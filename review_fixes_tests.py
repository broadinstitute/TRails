"""Regression tests for the issues found in the 2026-06-21 code review.

Each test pins a specific bug so it cannot silently return. The X-ids refer to
the consolidated review findings.
"""

import os
import sqlite3
import tempfile
import unittest

import analysis_columns
import input_tables
import locus_annotations
import result_database
import results_server


class AnnotationColumnClassificationTests(unittest.TestCase):
    """X4 / X5 / X10: matrix annotation columns are not treated as samples."""

    def _write(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write(text)
        handle.close()
        self.addCleanup(os.remove, handle.name)
        return handle.name

    def test_annotation_columns_separated_from_samples(self):
        path = self._write(
            "trid\tmotif\tgene_id\tHPRC256_99thPercentile\tSAMPLE_A\tSAMPLE_B\n"
            "chr1-100-115-CAG\tCAG\tENSG001\t30\t12,40\t11,12\n")
        locus_rows, sample_id_list = input_tables.read_repeat_copy_numbers(path)
        # gene_id / HPRC256_* are NOT samples.
        self.assertEqual(sample_id_list, ["SAMPLE_A", "SAMPLE_B"])
        self.assertNotIn("gene_id", sample_id_list)
        row = locus_rows[0]
        self.assertEqual(set(row["genotypes"]), {"SAMPLE_A", "SAMPLE_B"})
        # Annotation columns are promoted (and numeric ones coerced).
        self.assertEqual(row["extra_columns"]["gene_id"], "ENSG001")
        self.assertEqual(row["extra_columns"]["HPRC256_99thPercentile"], 30)

    def test_gencode_alias_and_case_insensitive(self):
        path = self._write(
            "trid\tmotif\tGencodeGeneId\tGencode_Gene_Region\tS1\n"
            "chr1-100-115-CAG\tCAG\tENSG9\tCDS\t10,10\n")
        locus_rows, sample_id_list = input_tables.read_repeat_copy_numbers(path)
        self.assertEqual(sample_id_list, ["S1"])
        self.assertEqual(locus_rows[0]["extra_columns"]["GencodeGeneId"], "ENSG9")
        self.assertEqual(locus_rows[0]["extra_columns"]["GencodeGeneRegion"], "CDS")


class GeneTableColumnsTests(unittest.TestCase):
    """X3: add_gene_columns populates the 10 GeneTable* columns."""

    def test_gene_table_columns_populated(self):
        records = [{"gene_id": "ENSG1"}]
        gene_lookup = {"ENSG1": {
            "gene_symbol": "FMR1", "gene_aliases": "FRAXA", "pLI_v2": "0.9",
            "pLI_v4": "0.95", "lof_oe_ci_upper_v4": "0.3", "hgnc_gene_id": "HGNC:3",
            "inheritance": "XL", "disease_category": "ID", "LLM_phenotype_summary": "x",
            "sources": "OMIM",
        }}
        analysis_columns.add_gene_columns(records, gene_lookup)
        self.assertEqual(records[0]["GeneTableGeneSymbol"], "FMR1")
        self.assertEqual(records[0]["GeneTableInheritance"], "XL")
        self.assertEqual(records[0]["GeneTableLoeuf"], "0.3")

    def test_gene_table_columns_none_without_match(self):
        records = [{"gene_id": "missing"}]
        analysis_columns.add_gene_columns(records, {"ENSG1": {"gene_symbol": "X"}})
        # Every GeneTable* column is explicitly None so the loci schema is consistent.
        for column in analysis_columns.GENE_TABLE_COLUMN_BY_SOURCE_FIELD.values():
            self.assertIsNone(records[0][column])


class EmptyMotifTests(unittest.TestCase):
    """X8: an empty motif must not divide-by-zero."""

    def test_empty_motif_does_not_crash(self):
        records = [{"LocusId": "1-100-115-", "Motif": ""}]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["MotifSize"], 0)
        self.assertIsNone(records[0]["NumRepeatsInReference"])


class EmptyLociTableTests(unittest.TestCase):
    """X2: an empty record set still writes a valid loci table."""

    def test_empty_records_create_full_schema(self):
        connection = sqlite3.connect(":memory:")
        row_count, present_columns = result_database.write_loci_table(
            connection, [], analysis_columns.OUTPUT_COLUMNS)
        self.assertEqual(row_count, 0)
        self.assertEqual(present_columns, list(analysis_columns.OUTPUT_COLUMNS))
        # The table exists and is queryable (no 'CREATE TABLE loci ()' crash).
        self.assertEqual(connection.execute("SELECT COUNT(*) FROM loci").fetchone()[0], 0)
        connection.close()


class PermissiveSampleTableTests(unittest.TestCase):
    """X1 / X12: the server's load_sample_table accepts minimal metadata."""

    def _write(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write(text)
        handle.close()
        self.addCleanup(os.remove, handle.name)
        return handle.name

    def test_only_sample_id_does_not_crash(self):
        path = self._write("sample_id\nS1\nS2\n")
        sample_rows, affected_lookup, analysis_lookup = results_server.load_sample_table(path)
        self.assertEqual(set(sample_rows), {"S1", "S2"})
        self.assertEqual(analysis_lookup["S1"], "unknown")
        # Absent affected_status is falsy (treated as not-unaffected downstream).
        self.assertFalse(affected_lookup["S1"])

    def test_case_insensitive_id_and_extra_columns_preserved(self):
        path = self._write("Sample ID\tancestry\nS1\tEUR\n")
        sample_rows, _affected, _analysis = results_server.load_sample_table(path)
        self.assertIn("S1", sample_rows)
        self.assertEqual(sample_rows["S1"]["ancestry"], "EUR")  # extra column preserved

    def test_phenotype_strips_only_leading_na_prefix(self):
        path = self._write(
            "sample_id\taffected_status\tphenotype_description\n"
            "S1\tAffected\tNA; seizures; NA; ataxia\n")
        sample_rows, _affected, _analysis = results_server.load_sample_table(path)
        # Only the leading 'NA; ' is removed (matches the build), not the inner one.
        self.assertEqual(sample_rows["S1"]["phenotype_description"], "seizures; NA; ataxia")


if __name__ == "__main__":
    unittest.main()
