"""End-to-end test for build_database.build() — the pipeline orchestrator.

Synthesizes a tiny cohort (4 loci including one chrX, 5 samples including one
complete trio and 2 phenotyped samples) in a temp dir, runs the full build, and
asserts the resulting database has the expected tables, row counts, column order,
and populated outlier columns — and that a re-run is idempotent.
"""

import contextlib
import io
import os
import shutil
import tempfile
import unittest

import analysis_columns
import build_database
import duckdb_compat
import result_database
from tests.database_test_helpers import ordered_table_column_names


def _write_tsv(path, header, rows):
    """Write a tab-separated table with the given header and row lists."""
    with open(path, "w") as output_file:
        output_file.write("\t".join(header) + "\n")
        for row in rows:
            output_file.write("\t".join(str(cell) for cell in row) + "\n")


class BuildDatabaseEndToEndTests(unittest.TestCase):

    def setUp(self):
        self.directory = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.directory)

        # Five samples: a complete trio (CHILD/MOTHER/FATHER) plus two singletons,
        # two of which (CHILD, SINGLE1) carry phenotype terms.
        self.matrix_path = os.path.join(self.directory, "matrix.tsv")
        self.metadata_path = os.path.join(self.directory, "samples.tsv")
        self.phenotypes_path = os.path.join(self.directory, "phenotypes.tsv")
        self.db_path = os.path.join(self.directory, "result.duckdb")

        sample_ids = ["CHILD", "MOTHER", "FATHER", "SINGLE1", "SINGLE2"]

        # Four loci: 3 autosomal + 1 chrX. The chrX child genotype is hemizygous.
        # Genotypes are chosen so the CHILD carries the largest expanded allele
        # at the autosomal loci (so it lands in the affected outlier columns).
        matrix_rows = [
            # chr1 trinucleotide locus: CHILD expanded.
            ["chr1-1000-1009-CAG", "CAG",
             "5,60", "5,6", "5,7", "5,8", "5,9"],
            # chr2 dinucleotide locus.
            ["chr2-2000-2010-AT", "AT",
             "10,50", "10,11", "10,12", "10,13", "10,14"],
            # chr3 tetranucleotide locus.
            ["chr3-3000-3012-AAAT", "AAAT",
             "3,45", "3,4", "3,5", "3,6", "3,7"],
            # chrX locus: CHILD (male) hemizygous and expanded.
            ["chrX-4000-4006-CG", "CG",
             "55", "8,9", "10", "12,13", "14,15"],
        ]
        _write_tsv(self.matrix_path, ["trid", "motif"] + sample_ids, matrix_rows)

        _write_tsv(
            self.metadata_path,
            ["sample_id", "sex", "family_id", "maternal_id", "paternal_id",
             "affected_status", "analysis_status", "phenotype_description"],
            [
                ["CHILD", "male", "FAM1", "MOTHER", "FATHER", "affected", "unsolved", "Ataxia"],
                ["MOTHER", "female", "FAM1", "", "", "unaffected", "unaffected", ""],
                ["FATHER", "male", "FAM1", "", "", "unaffected", "unaffected", ""],
                ["SINGLE1", "female", "FAM2", "", "", "affected", "unsolved", "Seizures"],
                ["SINGLE2", "male", "FAM3", "", "", "unaffected", "unaffected", ""],
            ])

        _write_tsv(
            self.phenotypes_path,
            ["participant_id", "term_id"],
            [
                ["CHILD", "HP:0001251"],   # ataxia
                ["CHILD", "HP:0002066"],
                ["SINGLE1", "HP:0001250"],  # seizures
            ])

    def test_build_produces_expected_database(self):
        returned = build_database.build(
            self.matrix_path, self.metadata_path, self.db_path,
            phenotypes_table=self.phenotypes_path)
        self.assertEqual(returned, self.db_path)
        self.assertTrue(os.path.exists(self.db_path))
        self.assertFalse(os.path.exists(self.db_path + ".tmp"))

        connection = duckdb_compat.connect(self.db_path)
        try:
            tables = duckdb_compat.list_tables(connection)
            for required in ["loci", "swim_plot", "metadata"]:
                self.assertIn(required, tables)

            # The callset size recorded for the results server: the matrix's five
            # sample columns.
            self.assertEqual(result_database.read_sample_count(connection), 5)

            # Phenotype tables: two phenotyped affected/unsolved samples carry HPO
            # terms, so qualifying samples exist and the tables are written.
            self.assertIn("per_outlier_phenotype_scores", tables)
            self.assertIn("per_locus_phenotype_scores", tables)

            # Mendelian tables: one complete trio (CHILD/MOTHER/FATHER) exists.
            self.assertIn("mendelian_violations", tables)
            self.assertIn("mendelian_violations_per_motif", tables)
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM mendelian_violations").fetchone()[0], 1)
            self.assertEqual(
                connection.execute("SELECT sample_id FROM mendelian_violations").fetchone()[0],
                "CHILD")

            # loci row count == number of loci.
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM loci").fetchone()[0], 4)

            # loci columns are a subset of OUTPUT_COLUMNS, in OUTPUT_COLUMNS order.
            loci_columns = ordered_table_column_names(connection, "loci")
            self.assertTrue(set(loci_columns).issubset(set(analysis_columns.OUTPUT_COLUMNS)))
            self.assertEqual(
                loci_columns,
                [c for c in analysis_columns.OUTPUT_COLUMNS if c in set(loci_columns)])

            # Core locus columns must be present and populated.
            for required_column in ["LocusId", "Motif", "CanonicalMotif",
                                    "AllAlleleHistogram", "OutlierSampleIds_AllAlleles",
                                    "Chrom", "Start0Based", "End1Based"]:
                self.assertIn(required_column, loci_columns)

            # OutlierSampleIds populated for at least one locus.
            outlier_values = [row[0] for row in connection.execute(
                "SELECT OutlierSampleIds_AllAlleles FROM loci")]
            self.assertTrue(any(value for value in outlier_values),
                            "expected at least one populated OutlierSampleIds_AllAlleles")

            # The CHILD's expanded chr1 allele must be the first affected allele.
            first_affected = connection.execute(
                "SELECT FirstAffectedAlleleSize_AllAlleles FROM loci WHERE LocusId=?",
                ("chr1-1000-1009-CAG",)).fetchone()[0]
            self.assertEqual(first_affected, 60)

            # swim_plot has one row per outlier entry (non-empty here).
            self.assertGreater(
                connection.execute("SELECT COUNT(*) FROM swim_plot").fetchone()[0], 0)
        finally:
            connection.close()

    def test_build_is_idempotent(self):
        build_database.build(
            self.matrix_path, self.metadata_path, self.db_path,
            phenotypes_table=self.phenotypes_path)
        first = self._snapshot()

        build_database.build(
            self.matrix_path, self.metadata_path, self.db_path,
            phenotypes_table=self.phenotypes_path)
        second = self._snapshot()

        self.assertEqual(first, second)

    def _snapshot(self):
        """Return a deterministic snapshot of the built database's loci table."""
        connection = duckdb_compat.connect(self.db_path)
        try:
            columns = ordered_table_column_names(connection, "loci")
            rows = connection.execute(
                "SELECT * FROM loci ORDER BY LocusId").fetchall()
            return columns, rows
        finally:
            connection.close()

    def test_build_runs_without_optional_inputs(self):
        # No phenotypes, no gene table, no known-loci catalog: the build still
        # produces a valid database; phenotype + Mendelian tables behave per their
        # skip rules (Mendelian still written because the trio exists).
        build_database.build(self.matrix_path, self.metadata_path, self.db_path)
        connection = duckdb_compat.connect(self.db_path)
        try:
            tables = duckdb_compat.list_tables(connection)
            self.assertIn("loci", tables)
            self.assertIn("swim_plot", tables)
            self.assertIn("mendelian_violations", tables)
            # No phenotypes supplied -> no phenotype tables.
            self.assertNotIn("per_outlier_phenotype_scores", tables)
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM loci").fetchone()[0], 4)
        finally:
            connection.close()

    def test_a_build_without_trexplorer_annotations_warns_that_the_population_gate_is_off(self):
        # Without TRExplorer's 99th percentiles (and none in the matrix), no locus is gated, so
        # the #Aff / #Fam counts and phenotype scores come out larger rather than empty.
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            build_database.build(self.matrix_path, self.metadata_path, self.db_path)
        self.assertIn("WARNING: no TRExplorer annotations database was supplied", output.getvalue())
        self.assertIn("HPRC256_99thPercentile / AoU1027_99thPercentile", output.getvalue())


def read_counts_implied_by_affected_unsolved_outliers_table(connection):
    """Recounts NumAffectedUnsolved{Samples,Families}AboveUnaffected from the outlier table.

    Args:
        connection: An open connection to a built result database.

    Returns:
        A dict mapping (LocusId, outlier_type) to (sample count, distinct non-NULL family
        count), for every pair the table has rows for.
    """
    return {(locus_id, outlier_type): (sample_count, family_count)
            for locus_id, outlier_type, sample_count, family_count in connection.execute(
                f"SELECT LocusId, outlier_type, COUNT(*), COUNT(DISTINCT family_id) "
                f"FROM {result_database.AFFECTED_UNSOLVED_OUTLIERS_TABLE} "
                f"GROUP BY LocusId, outlier_type").fetchall()}


class AffectedUnsolvedOutliersTableTests(unittest.TestCase):
    """The outlier table must list exactly the samples and families the loci counts count."""

    def setUp(self):
        self.directory = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.directory)
        self.matrix_path = os.path.join(self.directory, "matrix.tsv")
        self.metadata_path = os.path.join(self.directory, "samples.tsv")
        self.db_path = os.path.join(self.directory, "result.duckdb")

        # A1 and A2 share a family, A4 has no family_id, GHOST has no metadata row (so it
        # counts as affected without a family), U1 is unaffected and U2 is solved.
        sample_ids = ["A1", "A2", "A3", "A4", "U1", "U2", "GHOST"]
        _write_tsv(
            self.matrix_path,
            ["trid", "motif", "HPRC256_99thPercentile", "AoU1027_99thPercentile",
             "TenK10K_99thPercentile"] + sample_ids,
            [
                # Every affected allele is above U1's 40, and above the long-read 99th
                # percentiles (max 30), but only A1, GHOST and A2 are above TenK10K's 52.
                # Every ShortAlleles allele is 5, a tie with U1, so that count is 0.
                ["chr1-1000-1030-CAG", "CAG", "20", "30", "52",
                 "5,60", "5,55", "5,50", "5,45", "5,40", "5,35", "5,58"],
                # U1 carries the largest allele, so no sample is above every unaffected.
                ["chr2-2000-2020-AT", "AT", "", "", "",
                 "10,30", "10,12", "10,12", "10,12", "10,31", "10,12", "10,12"],
                # Hemizygous A1 and A3 are above U1 and HPRC256's 25; only A1 is above TenK10K's 39.
                ["chrX-3000-3010-CG", "CG", "25", "", "39",
                 "40", "8,9", "38", "8,9", "20", "8,9", "8,9"],
            ])
        _write_tsv(
            self.metadata_path,
            ["sample_id", "sex", "family_id", "affected_status", "analysis_status"],
            [
                ["A1", "male", "FAM_A", "affected", "unsolved"],
                ["A2", "female", "FAM_A", "affected", "unsolved"],
                ["A3", "male", "FAM_B", "affected", "unsolved"],
                ["A4", "female", "", "affected", "partially solved"],
                ["U1", "male", "FAM_C", "unaffected", "unaffected"],
                ["U2", "female", "FAM_D", "affected", "solved"],
            ])

    def build_and_read(self, **build_arguments):
        """Builds the database and returns (stored counts, counts implied by the table, gate cohorts).

        The stored counts map (LocusId, outlier_type) to the loci table's
        (NumAffectedUnsolvedSamplesAboveUnaffected, NumAffectedUnsolvedFamiliesAboveUnaffected).
        """
        build_database.build(self.matrix_path, self.metadata_path, self.db_path, **build_arguments)
        connection = duckdb_compat.connect(self.db_path)
        try:
            stored_counts = {}
            for outlier_type in analysis_columns.OUTLIER_TYPES:
                for locus_id, sample_count, family_count in connection.execute(
                        f"SELECT LocusId, NumAffectedUnsolvedSamplesAboveUnaffected_{outlier_type}, "
                        f"NumAffectedUnsolvedFamiliesAboveUnaffected_{outlier_type} FROM loci").fetchall():
                    stored_counts[(locus_id, outlier_type)] = (sample_count, family_count)
            return (stored_counts, read_counts_implied_by_affected_unsolved_outliers_table(connection),
                    result_database.read_population_gate_cohorts(connection))
        finally:
            connection.close()

    def assert_table_agrees_with_counts(self, stored_counts, counts_implied_by_table):
        """Checks that every (locus, outlier type) count matches the rows of the outlier table."""
        self.assertEqual(len(stored_counts), 3 * len(analysis_columns.OUTLIER_TYPES))
        for key, (sample_count, family_count) in stored_counts.items():
            # A NULL count (no outlier list) and a count of 0 both mean no rows.
            self.assertEqual(counts_implied_by_table.get(key, (0, 0)),
                             (sample_count or 0, family_count or 0), key)
        self.assertLessEqual(set(counts_implied_by_table), set(stored_counts))

    def test_default_gate_table_agrees_with_counts(self):
        stored_counts, counts_implied_by_table, population_gate_cohorts = self.build_and_read()
        self.assert_table_agrees_with_counts(stored_counts, counts_implied_by_table)
        self.assertEqual(population_gate_cohorts, ("HPRC256", "AoU1027"))
        self.assertEqual(stored_counts[("chr1-1000-1030-CAG", "AllAlleles")], (5, 2))
        self.assertEqual(stored_counts[("chr1-1000-1030-CAG", "ShortAlleles")], (0, 0))
        self.assertEqual(stored_counts[("chr1-1000-1030-CAG", "HemizygousAlleles")], (None, None))
        self.assertEqual(stored_counts[("chr2-2000-2020-AT", "AllAlleles")], (0, 0))
        self.assertEqual(stored_counts[("chrX-3000-3010-CG", "HemizygousAlleles")], (2, 2))

    def test_tenk10k_gate_table_agrees_with_counts(self):
        stored_counts, counts_implied_by_table, population_gate_cohorts = self.build_and_read(
            population_gate_cohorts=("HPRC256", "AoU1027", "TenK10K"))
        self.assert_table_agrees_with_counts(stored_counts, counts_implied_by_table)
        self.assertEqual(population_gate_cohorts, ("HPRC256", "AoU1027", "TenK10K"))
        self.assertEqual(stored_counts[("chr1-1000-1030-CAG", "AllAlleles")], (3, 1))
        self.assertEqual(stored_counts[("chrX-3000-3010-CG", "HemizygousAlleles")], (1, 1))

    def test_table_rows_carry_family_and_largest_allele(self):
        self.build_and_read()
        connection = duckdb_compat.connect(self.db_path)
        try:
            rows = connection.execute(
                f"SELECT * FROM {result_database.AFFECTED_UNSOLVED_OUTLIERS_TABLE} "
                f"WHERE LocusId = 'chr1-1000-1030-CAG' ORDER BY sample_id").fetchall()
        finally:
            connection.close()
        self.assertEqual(rows, [
            ("chr1-1000-1030-CAG", "AllAlleles", "A1", "FAM_A", 60),
            ("chr1-1000-1030-CAG", "AllAlleles", "A2", "FAM_A", 55),
            ("chr1-1000-1030-CAG", "AllAlleles", "A3", "FAM_B", 50),
            ("chr1-1000-1030-CAG", "AllAlleles", "A4", None, 45),
            ("chr1-1000-1030-CAG", "AllAlleles", "GHOST", None, 58),
        ])


class StrchiveOnlyCatalogTests(unittest.TestCase):
    """A STRchive catalog on its own, with no known-loci catalog, still annotates known loci."""

    def setUp(self):
        self.directory = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.directory)
        self.matrix_path = os.path.join(self.directory, "matrix.tsv")
        self.metadata_path = os.path.join(self.directory, "samples.tsv")
        self.strchive_path = os.path.join(self.directory, "strchive.json")
        self.db_path = os.path.join(self.directory, "result.duckdb")

        _write_tsv(self.matrix_path, ["trid", "motif", "S1", "S2"],
                   [["chr2-499-530-CAG", "CAG", "10,11", "10,12"]])
        _write_tsv(self.metadata_path, ["sample_id", "sex", "affected_status"],
                   [["S1", "male", "affected"], ["S2", "female", "unaffected"]])
        with open(self.strchive_path, "w") as handle:
            handle.write(
                '[{"id": "STR_A", "disease": "DiseaseA", "chrom": "chr2", '
                '"start_hg38": 500, "stop_hg38": 530, '
                '"reference_motif_reference_orientation": ["CAG"]}]')

    def test_a_locus_in_the_strchive_catalog_gets_its_known_disease_locus(self):
        build_database.build(
            self.matrix_path, self.metadata_path, self.db_path,
            strchive_loci_json=self.strchive_path)
        connection = duckdb_compat.connect(self.db_path)
        try:
            value = connection.execute(
                "SELECT KnownDiseaseLocus FROM loci WHERE LocusId = 'chr2-499-530-CAG'"
            ).fetchone()[0]
        finally:
            connection.close()
        self.assertEqual(value, "STR_A")


if __name__ == "__main__":
    unittest.main()
