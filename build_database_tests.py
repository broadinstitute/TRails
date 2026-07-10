"""End-to-end test for build_database.build() — the pipeline orchestrator.

Synthesizes a tiny cohort (4 loci including one chrX, 5 samples including one
complete trio and 2 phenotyped samples) in a temp dir, runs the full build, and
asserts the resulting database has the expected tables, row counts, column order,
and populated outlier columns — and that a re-run is idempotent.
"""

import os
import sqlite3
import tempfile
import unittest

import analysis_columns
import build_database


def _write_tsv(path, header, rows):
    """Write a tab-separated table with the given header and row lists."""
    with open(path, "w") as output_file:
        output_file.write("\t".join(header) + "\n")
        for row in rows:
            output_file.write("\t".join(str(cell) for cell in row) + "\n")


def _table_names(connection):
    return {row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}


def _table_columns(connection, table_name):
    return [row[1] for row in connection.execute(
        f"SELECT * FROM pragma_table_info('{table_name}')")]


class BuildDatabaseEndToEndTests(unittest.TestCase):

    def setUp(self):
        self.directory = tempfile.mkdtemp()

        # Five samples: a complete trio (CHILD/MOTHER/FATHER) plus two singletons,
        # two of which (CHILD, SINGLE1) carry phenotype terms.
        self.matrix_path = os.path.join(self.directory, "matrix.tsv")
        self.metadata_path = os.path.join(self.directory, "samples.tsv")
        self.phenotypes_path = os.path.join(self.directory, "phenotypes.tsv")
        self.db_path = os.path.join(self.directory, "result.db")

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

    def tearDown(self):
        for name in os.listdir(self.directory):
            os.remove(os.path.join(self.directory, name))
        os.rmdir(self.directory)

    def test_build_produces_expected_database(self):
        returned = build_database.build(
            self.matrix_path, self.metadata_path, self.db_path,
            phenotypes_table=self.phenotypes_path)
        self.assertEqual(returned, self.db_path)
        self.assertTrue(os.path.exists(self.db_path))
        self.assertFalse(os.path.exists(self.db_path + ".tmp"))

        connection = sqlite3.connect(self.db_path)
        try:
            tables = _table_names(connection)
            for required in ["loci", "swim_plot", "sk_AllAlleles", "sk_ShortAlleles",
                             "sk_HemizygousAlleles"]:
                self.assertIn(required, tables)

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
            loci_columns = _table_columns(connection, "loci")
            self.assertTrue(set(loci_columns).issubset(set(analysis_columns.OUTPUT_COLUMNS)))
            self.assertEqual(
                loci_columns,
                [c for c in analysis_columns.OUTPUT_COLUMNS if c in set(loci_columns)])

            # Core Tier-A columns must be present and populated.
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
        connection = sqlite3.connect(self.db_path)
        try:
            columns = _table_columns(connection, "loci")
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
        connection = sqlite3.connect(self.db_path)
        try:
            tables = _table_names(connection)
            self.assertIn("loci", tables)
            self.assertIn("swim_plot", tables)
            self.assertIn("mendelian_violations", tables)
            # No phenotypes supplied -> no phenotype tables.
            self.assertNotIn("per_outlier_phenotype_scores", tables)
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM loci").fetchone()[0], 4)
        finally:
            connection.close()


if __name__ == "__main__":
    unittest.main()
