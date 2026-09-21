"""Unit tests for result_database.py, the DuckDB writer layer."""

import os
import subprocess
import sys
import tempfile
import unittest

import duckdb_compat
import phenotype_scoring
import result_database


def leave_a_killed_writers_wal(db_path, table_name):
    """Leaves a real DuckDB write-ahead log next to ``db_path``, as a killed writer does.

    Writes ``table_name`` into a read-write database at ``db_path`` from a subprocess that
    then kills itself, so the log is never folded back into the database file.

    Args:
        db_path: The database path to write to and leave a log beside.
        table_name: The table the log will hold.
    """
    subprocess.run(
        [sys.executable, "-c",
         "import os, sys, duckdb\n"
         "connection = duckdb.connect(sys.argv[1])\n"
         "connection.execute('CREATE TABLE ' + sys.argv[2] + ' (a INTEGER)')\n"
         "connection.execute('INSERT INTO ' + sys.argv[2] + ' VALUES (1)')\n"
         "os.kill(os.getpid(), 9)\n",
         db_path, table_name],
        check=False)


def table_columns(connection, table_name):
    """Returns the ordered list of column names of ``table_name``."""
    return [column[0] for column in
            connection.execute(f"SELECT * FROM {table_name} LIMIT 0").description]


def table_names(connection):
    """Returns the set of table names in the database."""
    return duckdb_compat.list_tables(connection)


def primary_key_columns(connection, table_name):
    """Returns the columns declared PRIMARY KEY on ``table_name``."""
    row = connection.execute(
        "SELECT constraint_column_names FROM duckdb_constraints() "
        "WHERE table_name = ? AND constraint_type = 'PRIMARY KEY'",
        (table_name,)).fetchone()
    return list(row[0]) if row else []


class InferColumnTypesTests(unittest.TestCase):
    """DuckDB columns are statically typed, so every column needs a declared type."""

    def test_type_follows_the_widest_value_kind(self):
        rows = [
            {"count": 1, "score": 0.5, "name": "a", "mixed": 1},
            {"count": 2, "score": 2, "name": "b", "mixed": "not_available"},
        ]
        types = result_database.infer_column_types(
            ["count", "score", "name", "mixed"], rows)
        self.assertEqual(types["count"], "BIGINT")
        # One float among ints widens the column to DOUBLE.
        self.assertEqual(types["score"], "DOUBLE")
        self.assertEqual(types["name"], "VARCHAR")
        # A column carrying both numbers and a text marker has to be VARCHAR, or its
        # insert fails.
        self.assertEqual(types["mixed"], "VARCHAR")

    def test_all_null_column_falls_back_by_name(self):
        rows = [{"pLI": None, "GeneTableGeneSymbol": None}]
        types = result_database.infer_column_types(
            ["pLI", "GeneTableGeneSymbol"], rows)
        # pLI is compared to a number by the results server, and comparing a number
        # against a VARCHAR column is a binder error in DuckDB.
        self.assertEqual(types["pLI"], "DOUBLE")
        self.assertEqual(types["GeneTableGeneSymbol"], "VARCHAR")

    def test_a_real_value_overrides_the_name_fallback(self):
        types = result_database.infer_column_types(["pLI"], [{"pLI": "unknown"}])
        self.assertEqual(types["pLI"], "VARCHAR")


class WriteLociTableTests(unittest.TestCase):

    def setUp(self):
        self.connection = duckdb_compat.connect(":memory:")
        # A small subset of OUTPUT_COLUMNS, deliberately out of input-dict order.
        self.output_columns = [
            "LocusId", "Motif", "CanonicalMotif", "MotifSize",
            "FirstAffectedAlleleSize_AllAlleles", "Chrom",
        ]

    def tearDown(self):
        self.connection.close()

    def test_writes_columns_in_output_order_only_when_present(self):
        records = [
            {"LocusId": "1-1-10-A", "Motif": "A", "CanonicalMotif": "A",
             "MotifSize": 1, "Chrom": "chr1"},
            {"LocusId": "2-5-20-AG", "Motif": "AG", "CanonicalMotif": "AG",
             "MotifSize": 2, "Chrom": "chr2"},
        ]
        row_count = result_database.write_loci_table(
            self.connection, records, self.output_columns)

        self.assertEqual(row_count, 2)
        # FirstAffectedAlleleSize_AllAlleles is absent from every record, so it
        # must NOT appear; the rest keep OUTPUT_COLUMNS order (not dict order).
        self.assertEqual(table_columns(self.connection, "loci"),
                         ["LocusId", "Motif", "CanonicalMotif", "MotifSize", "Chrom"])

    def test_row_count_and_values_round_trip(self):
        records = [{"LocusId": "X-1-2-T", "Motif": "T", "MotifSize": 1}]
        result_database.write_loci_table(self.connection, records, self.output_columns)
        rows = self.connection.execute(
            "SELECT LocusId, Motif, MotifSize FROM loci").fetchall()
        self.assertEqual(rows, [("X-1-2-T", "T", 1)])

    def test_absent_key_in_a_present_column_is_null(self):
        # CanonicalMotif is present overall (second record has it) but the first
        # record omits it -> that cell must be NULL, not an error.
        records = [
            {"LocusId": "1-1-10-A", "Motif": "A"},
            {"LocusId": "2-1-10-A", "Motif": "A", "CanonicalMotif": "A"},
        ]
        result_database.write_loci_table(
            self.connection, records, self.output_columns)
        self.assertIn("CanonicalMotif", table_columns(self.connection, "loci"))
        values = self.connection.execute(
            "SELECT CanonicalMotif FROM loci ORDER BY LocusId").fetchall()
        self.assertEqual(values, [(None,), ("A",)])

    def test_mixed_value_column_round_trips_as_text(self):
        # An annotation column carrying both counts and a text marker must not fail
        # its insert on the way into a statically typed column.
        records = [
            {"LocusId": "1-1-10-A", "Motif": 4},
            {"LocusId": "2-1-10-A", "Motif": "not_available"},
        ]
        result_database.write_loci_table(self.connection, records, self.output_columns)
        values = self.connection.execute(
            "SELECT Motif FROM loci ORDER BY LocusId").fetchall()
        self.assertEqual(values, [("4",), ("not_available",)])

    def test_a_number_beside_a_null_in_a_text_column_keeps_its_spelling(self):
        # Same mixed column as above, with a record that omits it. The rows are loaded in
        # bulk, and a bulk loader that let pandas type the column would widen 4 to 4.0 and
        # write "4.0" into the VARCHAR column.
        records = [
            {"LocusId": "1-1-10-A", "Motif": 4},
            {"LocusId": "2-1-10-A"},
            {"LocusId": "3-1-10-A", "Motif": "not_available"},
        ]
        result_database.write_loci_table(self.connection, records, self.output_columns)
        values = self.connection.execute(
            "SELECT Motif FROM loci ORDER BY LocusId").fetchall()
        self.assertEqual(values, [("4",), (None,), ("not_available",)])

    def test_a_mixed_column_survives_a_marker_past_the_type_sample(self):
        # The same mixed column as above, at build scale. The loader hands DuckDB a chunk of
        # thousands of rows at a time and DuckDB used to type each column from a sample of
        # 1000 of them, so a marker at row 2500 was invisible: the column came out INT32 and
        # the insert failed with "Could not convert string 'not_available' to INT32".
        records = [{"LocusId": f"{i}-1-10-A", "Motif": i} for i in range(3000)]
        records[2500]["Motif"] = "not_available"
        result_database.write_loci_table(self.connection, records, self.output_columns)
        self.assertEqual(
            self.connection.execute("SELECT COUNT(*) FROM loci").fetchone()[0], 3000)
        self.assertEqual(
            self.connection.execute(
                "SELECT Motif FROM loci WHERE LocusId = '2500-1-10-A'").fetchone()[0],
            "not_available")
        self.assertEqual(
            self.connection.execute(
                "SELECT Motif FROM loci WHERE LocusId = '2499-1-10-A'").fetchone()[0], "2499")

    def test_all_null_numeric_column_is_comparable_to_a_number(self):
        # A build run without a gene table still writes pLI, all NULL, and the server
        # compares it to a number.
        records = [{"LocusId": "1-1-10-A", "pLI": None}]
        result_database.write_loci_table(
            self.connection, records, ["LocusId", "pLI"])
        self.assertEqual(
            self.connection.execute(
                "SELECT COUNT(*) FROM loci WHERE pLI > 0.9").fetchone()[0], 0)


class WriteSwimPlotTests(unittest.TestCase):

    def setUp(self):
        self.connection = duckdb_compat.connect(":memory:")

    def tearDown(self):
        self.connection.close()

    def _swim_row(self, sample_id, allele_size):
        return {
            "outlier_type": "AllAlleles", "outlier_rank": 1,
            "motif_category": "3bp", "SourceDb": "TRails",
            "allele_size": allele_size, "sample_id": sample_id,
            "family_id": "fam", "affected_status": "Affected",
            "analysis_status": "Unsolved", "sex": "male",
            "phenotype_description": None, "purity": None, "methylation": None,
            "FirstUnaffectedAlleleSize": None, "is_above_first_unaffected": 1,
            "LocusId": "1-1-10-AAG", "Motif": "AAG", "CanonicalMotif": "AAG",
            "MotifSize": 3, "gene_region": "intron",
            "GeneTableGeneSymbol": None, "IsInMendelianGene": 0, "IsKnownMotif": 0,
            "gene_id": None, "pLI": None, "NumRepeatsInReference": 3,
        }

    def test_row_count_and_column_order(self):
        rows = [self._swim_row("S1", 40), self._swim_row("S2", 39)]
        count = result_database.write_swim_plot(self.connection, rows)
        self.assertEqual(count, 2)
        self.assertEqual(
            self.connection.execute("SELECT COUNT(*) FROM swim_plot").fetchone()[0], 2)
        # Column order follows the producing dict's key order.
        self.assertEqual(table_columns(self.connection, "swim_plot")[0], "outlier_type")

    def test_empty_swim_rows_skips_table(self):
        count = result_database.write_swim_plot(self.connection, [])
        self.assertEqual(count, 0)
        self.assertNotIn("swim_plot", table_names(self.connection))

    def test_empty_swim_rows_with_columns_creates_the_table(self):
        # A build can legitimately find no outlier alleles; the server's swim-plot
        # endpoints still have to find a queryable table.
        count = result_database.write_swim_plot(
            self.connection, [], columns=list(self._swim_row("S1", 40).keys()))
        self.assertEqual(count, 0)
        self.assertIn("swim_plot", table_names(self.connection))
        self.assertEqual(
            self.connection.execute(
                "SELECT COUNT(*) FROM swim_plot WHERE allele_size > 10").fetchone()[0], 0)


class WritePhenotypeTablesTests(unittest.TestCase):

    def setUp(self):
        self.connection = duckdb_compat.connect(":memory:")

    def tearDown(self):
        self.connection.close()

    def test_empty_skips_both_tables(self):
        per_outlier, per_locus = result_database.write_phenotype_tables(
            self.connection, [], [])
        self.assertEqual((per_outlier, per_locus), (0, 0))
        self.assertNotIn("per_outlier_phenotype_scores", table_names(self.connection))
        self.assertNotIn("per_locus_phenotype_scores", table_names(self.connection))

    def _rows_from_the_producer(self):
        """Returns ``(per_outlier_rows, per_locus_rows)`` from a small real scoring run.

        The writer takes each table's column order from its first row dict, so a
        hand-written fixture would only pin the shape the fixture itself declares.
        phenotype_scoring.compute_phenotype_scores is the only producer of these
        rows in the build, so running it keeps this test tied to the schema the
        build actually writes.

        Returns:
            A ``(per_outlier_rows, per_locus_rows)`` tuple.
        """
        return phenotype_scoring.compute_phenotype_scores(
            [{"LocusId": "1-1-10-A", "gene_id": "ENSG001",
              "OutlierSampleIds_AllAlleles": "40x:S1,30x:S2"}],
            {"S1": {"HP:0001250"}, "S2": {"HP:0001250"}},
            {"ENSG001": {"gene_symbol": "FXN"}},
            {},
            {"S1": "affected", "S2": "affected"},
            {"S1": "unsolved", "S2": "unsolved"},
        )

    def test_writes_both_tables_when_non_empty(self):
        per_outlier_rows, per_locus_rows = self._rows_from_the_producer()
        per_outlier, per_locus = result_database.write_phenotype_tables(
            self.connection, per_outlier_rows, per_locus_rows)
        self.assertEqual((per_outlier, per_locus), (2, 1))
        # The expected names are spelled out here rather than derived from the rows,
        # so a change to either phenotype-score schema fails this test instead of
        # passing silently.
        self.assertEqual(
            table_columns(self.connection, "per_outlier_phenotype_scores"),
            ["locus_id", "sample_id", "outlier_type", "allele_size", "gene_symbol",
             "gene_phenotype_similarity", "gene_phenotype_overlap_count",
             "n_matching_diseases", "best_matching_disease",
             "best_disease_inheritance", "pairwise_similarity_to_next",
             "pairwise_shared_count_raw", "pairwise_shared_count_ic",
             "next_sample_id"])
        self.assertEqual(
            table_columns(self.connection, "per_locus_phenotype_scores"),
            ["locus_id", "outlier_type", "num_qualifying_samples",
             "sum_pairwise_similarity", "sum_pairwise_shared_raw",
             "sum_pairwise_shared_ic", "max_gene_phenotype_similarity",
             "qualifying_sample_ids"])
        self.assertEqual(
            self.connection.execute(
                "SELECT qualifying_sample_ids FROM per_locus_phenotype_scores").fetchone(),
            ("S1,S2",))


class WriteMendelianTablesTests(unittest.TestCase):

    def setUp(self):
        self.connection = duckdb_compat.connect(":memory:")

    def tearDown(self):
        self.connection.close()

    def test_empty_skips_tables(self):
        per_sample, per_motif = result_database.write_mendelian_tables(
            self.connection, [], [])
        self.assertEqual((per_sample, per_motif), (0, 0))
        self.assertNotIn("mendelian_violations", table_names(self.connection))

    def test_writes_tables_with_sample_id_primary_key(self):
        per_sample_rows = [{
            "sample_id": "CHILD1", "autosome_violations": 2, "autosome_total": 100,
            "total_violations": 2, "total_loci": 100,
        }]
        per_motif_rows = [{"sample_id": "CHILD1", "mv_A": 1, "total_A": 50}]
        per_sample, per_motif = result_database.write_mendelian_tables(
            self.connection, per_sample_rows, per_motif_rows)
        self.assertEqual((per_sample, per_motif), (1, 1))
        self.assertEqual(table_columns(self.connection, "mendelian_violations")[0], "sample_id")
        # sample_id declared as PRIMARY KEY.
        self.assertEqual(
            primary_key_columns(self.connection, "mendelian_violations"), ["sample_id"])
        self.assertEqual(
            self.connection.execute(
                "SELECT autosome_violations FROM mendelian_violations").fetchone(), (2,))


class OpenAndFinalizeTests(unittest.TestCase):

    def test_open_new_database_uses_tmp_path(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.duckdb")
            connection, tmp_path = result_database.open_new_database(final_path)
            self.assertEqual(tmp_path, final_path + ".tmp")
            self.assertTrue(os.path.exists(tmp_path))
            self.assertFalse(os.path.exists(final_path))
            connection.close()

    def test_open_new_database_removes_stale_tmp(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.duckdb")
            with open(final_path + ".tmp", "w") as stale:
                stale.write("stale garbage that is not a database file")
            # A stale write-ahead log belongs to the database that was just removed,
            # so it has to go with it.
            with open(final_path + ".tmp.wal", "w") as stale_wal:
                stale_wal.write("stale write-ahead log")
            connection, tmp_path = result_database.open_new_database(final_path)
            # A fresh, valid, empty database replaced the stale files.
            self.assertEqual(len(duckdb_compat.list_tables(connection)), 0)
            connection.close()

    def test_finalize_database_atomic_replace(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.duckdb")
            connection, tmp_path = result_database.open_new_database(final_path)
            result_database.write_loci_table(
                connection,
                [{"LocusId": "1-1-10-A", "Motif": "A"}],
                ["LocusId", "Motif"])
            result_database.finalize_database(connection, tmp_path, final_path)

            self.assertTrue(os.path.exists(final_path))
            self.assertFalse(os.path.exists(tmp_path))
            # The close folded the write-ahead log back in, so the finished database
            # is the single file the move just put in place.
            self.assertEqual(duckdb_compat.find_leftover_write_files(final_path), [])
            reopened = duckdb_compat.connect(final_path)
            self.assertEqual(
                reopened.execute("SELECT COUNT(*) FROM loci").fetchone()[0], 1)
            reopened.close()

    def test_finalize_overwrites_existing_final(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.duckdb")
            with open(final_path, "w") as existing:
                existing.write("an old database")
            connection, tmp_path = result_database.open_new_database(final_path)
            result_database.write_loci_table(
                connection, [{"LocusId": "X-1-2-T"}], ["LocusId"])
            result_database.finalize_database(connection, tmp_path, final_path)
            reopened = duckdb_compat.connect(final_path)
            self.assertEqual(
                reopened.execute("SELECT LocusId FROM loci").fetchone(), ("X-1-2-T",))
            reopened.close()

    def test_finalize_removes_the_destinations_stale_write_ahead_log(self):
        # DuckDB replays whatever "<db>.wal" it finds next to a database without checking
        # that the log belongs to it, so a log left by a killed writer would be folded into
        # the database the move just put there: the previous build's tables come back, and a
        # log that re-creates "loci" makes every open fail while replaying it.
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.duckdb")
            leave_a_killed_writers_wal(final_path, "old_table")
            self.assertTrue(os.path.exists(final_path + ".wal"),
                            "the killed writer left no write-ahead log to test against")

            connection, tmp_path = result_database.open_new_database(final_path)
            result_database.write_loci_table(
                connection, [{"LocusId": "X-1-2-T"}], ["LocusId"])
            result_database.finalize_database(connection, tmp_path, final_path)

            self.assertEqual(duckdb_compat.find_leftover_write_files(final_path), [])
            reopened = duckdb_compat.connect(final_path)
            try:
                self.assertEqual(duckdb_compat.list_tables(reopened), {"loci"})
                self.assertEqual(
                    reopened.execute("SELECT LocusId FROM loci").fetchone(), ("X-1-2-T",))
            finally:
                reopened.close()

    def test_finalize_keeps_the_destinations_wal_when_the_move_fails(self):
        # Until the move succeeds, that log is the missing tail of a database that is still
        # there, so the removal has to come after os.replace rather than before it.
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.duckdb")
            leave_a_killed_writers_wal(final_path, "old_table")
            connection, tmp_path = result_database.open_new_database(final_path)
            result_database.write_loci_table(
                connection, [{"LocusId": "X-1-2-T"}], ["LocusId"])
            os.remove(tmp_path)
            with self.assertRaises(OSError):
                result_database.finalize_database(connection, tmp_path, final_path)
            self.assertTrue(os.path.exists(final_path + ".wal"))


if __name__ == "__main__":
    unittest.main()
