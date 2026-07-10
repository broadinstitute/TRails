"""Unit tests for result_database.py — the pure SQLite writer layer."""

import os
import sqlite3
import tempfile
import unittest

import result_database


def table_columns(connection, table_name):
    """Returns the ordered list of column names of ``table_name``."""
    return [row[1] for row in connection.execute(
        f"SELECT * FROM pragma_table_info('{table_name}')")]


def table_names(connection):
    """Returns the set of table names in the database."""
    return {row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}


def index_names(connection, table_name):
    """Returns the list of index names defined on ``table_name``."""
    return [row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='index' AND tbl_name=?",
        (table_name,))]


class WriteLociTableTests(unittest.TestCase):

    def setUp(self):
        self.connection = sqlite3.connect(":memory:")
        # A small subset of OUTPUT_COLUMNS, deliberately out of input-dict order.
        self.output_columns = [
            "LocusId", "Motif", "CanonicalMotif", "MotifSize",
            "FirstAffectedAlleleSize_AllAlleles", "Chrom",
        ]

    def test_writes_columns_in_output_order_only_when_present(self):
        records = [
            {"LocusId": "1-1-10-A", "Motif": "A", "CanonicalMotif": "A",
             "MotifSize": 1, "Chrom": "chr1"},
            {"LocusId": "2-5-20-AG", "Motif": "AG", "CanonicalMotif": "AG",
             "MotifSize": 2, "Chrom": "chr2"},
        ]
        row_count, present = result_database.write_loci_table(
            self.connection, records, self.output_columns)

        self.assertEqual(row_count, 2)
        # FirstAffectedAlleleSize_AllAlleles is absent from every record, so it
        # must NOT appear; the rest keep OUTPUT_COLUMNS order (not dict order).
        self.assertEqual(present,
                         ["LocusId", "Motif", "CanonicalMotif", "MotifSize", "Chrom"])
        self.assertEqual(table_columns(self.connection, "loci"), present)

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
        _, present = result_database.write_loci_table(
            self.connection, records, self.output_columns)
        self.assertIn("CanonicalMotif", present)
        values = self.connection.execute(
            "SELECT CanonicalMotif FROM loci ORDER BY LocusId").fetchall()
        self.assertEqual(values, [(None,), ("A",)])

    def test_create_loci_indexes_guarded_by_present_columns(self):
        records = [{"LocusId": "1-1-10-A", "Motif": "A", "Chrom": "chr1"}]
        _, present = result_database.write_loci_table(
            self.connection, records, self.output_columns)
        created = result_database.create_loci_indexes(self.connection, present)
        # LocusId / Chrom present -> their indexes exist; gene_id absent -> none.
        self.assertIn("idx_loci_LocusId", created)
        self.assertIn("idx_loci_Chrom", created)
        self.assertNotIn("idx_loci_gene_id", created)
        self.assertIn("idx_loci_LocusId", index_names(self.connection, "loci"))


class WriteSwimPlotTests(unittest.TestCase):

    def setUp(self):
        self.connection = sqlite3.connect(":memory:")

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

    def test_row_count_and_index_present(self):
        rows = [self._swim_row("S1", 40), self._swim_row("S2", 39)]
        count = result_database.write_swim_plot(self.connection, rows)
        self.assertEqual(count, 2)
        self.assertEqual(
            self.connection.execute("SELECT COUNT(*) FROM swim_plot").fetchone()[0], 2)
        self.assertIn("idx_swim_allele_size", index_names(self.connection, "swim_plot"))
        # Column order follows the producing dict's key order.
        self.assertEqual(table_columns(self.connection, "swim_plot")[0], "outlier_type")

    def test_empty_swim_rows_skips_table(self):
        count = result_database.write_swim_plot(self.connection, [])
        self.assertEqual(count, 0)
        self.assertNotIn("swim_plot", table_names(self.connection))


class WriteSkinnyTablesTests(unittest.TestCase):

    def test_builds_three_skinny_tables(self):
        connection = sqlite3.connect(":memory:")
        connection.execute(
            "CREATE TABLE loci (LocusId, Motif, MotifSize, "
            "FirstAffectedAlleleSize_AllAlleles, HemizygousAlleleHistogram)")
        connection.execute(
            "INSERT INTO loci VALUES ('1-1-10-A', 'A', 1, 40, '3x:1')")
        connection.commit()
        names = result_database.write_skinny_tables(
            connection,
            ["LocusId", "Motif", "MotifSize",
             "FirstAffectedAlleleSize_AllAlleles", "HemizygousAlleleHistogram"])
        self.assertEqual(
            names, ["sk_AllAlleles", "sk_ShortAlleles", "sk_HemizygousAlleles"])
        self.assertEqual(
            connection.execute("SELECT COUNT(*) FROM sk_AllAlleles").fetchone()[0], 1)


class WritePhenotypeTablesTests(unittest.TestCase):

    def setUp(self):
        self.connection = sqlite3.connect(":memory:")

    def test_empty_skips_both_tables(self):
        per_outlier, per_locus = result_database.write_phenotype_tables(
            self.connection, [], [])
        self.assertEqual((per_outlier, per_locus), (0, 0))
        self.assertNotIn("per_outlier_phenotype_scores", table_names(self.connection))
        self.assertNotIn("per_locus_phenotype_scores", table_names(self.connection))

    def test_writes_both_tables_when_non_empty(self):
        per_outlier_rows = [{
            "locus_id": "1-1-10-A", "sample_id": "S1", "outlier_type": "AllAlleles",
            "allele_size": 40, "gene_symbol": "FXN",
            "gene_phenotype_similarity": 0.5, "gene_phenotype_overlap_count": 2.0,
            "n_matching_diseases": 1.0, "best_matching_disease": "FRDA",
            "best_disease_inheritance": "AR", "pairwise_similarity_to_next": None,
            "pairwise_shared_count_raw": None, "pairwise_shared_count_ic": None,
            "next_sample_id": None,
        }]
        per_locus_rows = [{
            "locus_id": "1-1-10-A", "outlier_type": "AllAlleles",
            "num_qualifying_samples": 1, "sum_pairwise_similarity": None,
            "sum_pairwise_shared_raw": None, "sum_pairwise_shared_ic": None,
            "max_gene_phenotype_similarity": 0.5, "qualifying_sample_ids": "S1",
            "IsKnownMotif": 0, "gene_region": "intron", "gene_region_rank": 5,
            "FirstAffectedAlleleSize": 40, "FirstUnaffectedAlleleSize": None,
            "NumRepeatsInReference": 1, "HPRC256_MaxAllele": None,
            "AoU1027_MaxAllele": None, "TenK10K_MaxAllele": None,
            "NumAffectedAboveUnaffected": 1, "NumAffectedFamiliesAboveUnaffected": 1,
            "MotifSize": 1,
        }]
        per_outlier, per_locus = result_database.write_phenotype_tables(
            self.connection, per_outlier_rows, per_locus_rows)
        self.assertEqual((per_outlier, per_locus), (1, 1))
        self.assertEqual(len(table_columns(self.connection, "per_outlier_phenotype_scores")), 14)
        self.assertEqual(len(table_columns(self.connection, "per_locus_phenotype_scores")), 20)
        self.assertIn("idx_locus_pheno_id",
                      index_names(self.connection, "per_locus_phenotype_scores"))


class WriteMendelianTablesTests(unittest.TestCase):

    def setUp(self):
        self.connection = sqlite3.connect(":memory:")

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
        pk_columns = [row[1] for row in self.connection.execute(
            "SELECT * FROM pragma_table_info('mendelian_violations')") if row[5]]
        self.assertEqual(pk_columns, ["sample_id"])
        self.assertEqual(
            self.connection.execute(
                "SELECT autosome_violations FROM mendelian_violations").fetchone(), (2,))


class OpenAndFinalizeTests(unittest.TestCase):

    def test_open_new_database_uses_tmp_path(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.db")
            connection, tmp_path = result_database.open_new_database(final_path)
            self.assertEqual(tmp_path, final_path + ".tmp")
            self.assertTrue(os.path.exists(tmp_path))
            self.assertFalse(os.path.exists(final_path))
            connection.close()

    def test_open_new_database_removes_stale_tmp(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.db")
            with open(final_path + ".tmp", "w") as stale:
                stale.write("stale garbage that is not a sqlite file")
            connection, tmp_path = result_database.open_new_database(final_path)
            # A fresh, valid, empty sqlite db replaced the stale file.
            self.assertEqual(
                connection.execute(
                    "SELECT COUNT(*) FROM sqlite_master").fetchone()[0], 0)
            connection.close()

    def test_finalize_database_atomic_replace(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.db")
            connection, tmp_path = result_database.open_new_database(final_path)
            result_database.write_loci_table(
                connection,
                [{"LocusId": "1-1-10-A", "Motif": "A"}],
                ["LocusId", "Motif"])
            result_database.finalize_database(connection, tmp_path, final_path)

            self.assertTrue(os.path.exists(final_path))
            self.assertFalse(os.path.exists(tmp_path))
            reopened = sqlite3.connect(final_path)
            self.assertEqual(
                reopened.execute("SELECT COUNT(*) FROM loci").fetchone()[0], 1)
            reopened.close()

    def test_finalize_overwrites_existing_final(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.db")
            with open(final_path, "w") as existing:
                existing.write("an old database")
            connection, tmp_path = result_database.open_new_database(final_path)
            result_database.write_loci_table(
                connection, [{"LocusId": "X-1-2-T"}], ["LocusId"])
            result_database.finalize_database(connection, tmp_path, final_path)
            reopened = sqlite3.connect(final_path)
            self.assertEqual(
                reopened.execute("SELECT LocusId FROM loci").fetchone(), ("X-1-2-T",))
            reopened.close()


if __name__ == "__main__":
    unittest.main()
