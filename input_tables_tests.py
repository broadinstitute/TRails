"""Unit tests for input_tables.py."""

import gzip
import os
import tempfile
import unittest

import pandas

import input_tables


class NormalizeColumnNameTests(unittest.TestCase):

    def test_collapses_case_underscore_and_space(self):
        self.assertEqual(input_tables.normalize_column_name("Sample_ID"), "sampleid")
        self.assertEqual(input_tables.normalize_column_name("sample id"), "sampleid")
        self.assertEqual(input_tables.normalize_column_name("  SAMPLEID "), "sampleid")
        self.assertEqual(input_tables.normalize_column_name("sampleid"), "sampleid")


class MatchColumnsTests(unittest.TestCase):

    def test_case_and_underscore_insensitive(self):
        df = pandas.DataFrame(columns=["TRID", "Motif", "Sample_A"])
        matches = input_tables.match_columns(df, ["trid", "motif", "missing"])
        self.assertEqual(matches["trid"], "TRID")
        self.assertEqual(matches["motif"], "Motif")
        self.assertNotIn("missing", matches)

    def test_first_match_wins(self):
        df = pandas.DataFrame(columns=["sample_id", "sampleid"])
        matches = input_tables.match_columns(df, ["sample_id"])
        self.assertEqual(matches["sample_id"], "sample_id")


def _write_tsv(rows, header):
    handle = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    handle.write("\t".join(header) + "\n")
    for row in rows:
        handle.write("\t".join(str(c) for c in row) + "\n")
    handle.close()
    return handle.name


class ReadRepeatCopyNumbersTests(unittest.TestCase):

    def test_basic_parse_and_sample_list(self):
        path = _write_tsv(
            [["chr1-100-110-AT", "AT", "12,40", "21"]],
            ["TRID", "Motif", "SampleA", "SampleB"],
        )
        try:
            locus_rows, sample_ids = input_tables.read_repeat_copy_numbers(path)
        finally:
            os.unlink(path)
        self.assertEqual(sample_ids, ["SampleA", "SampleB"])
        self.assertEqual(locus_rows[0]["trid"], "chr1-100-110-AT")
        self.assertEqual(locus_rows[0]["motif"], "AT")
        self.assertEqual(locus_rows[0]["genotypes"], {"SampleA": "12,40", "SampleB": "21"})

    def test_gzip_read(self):
        handle = tempfile.NamedTemporaryFile(suffix=".tsv.gz", delete=False)
        handle.close()
        with gzip.open(handle.name, "wt") as f:
            f.write("trid\tmotif\tS1\n")
            f.write("chr2-5-11-AAT\tAAT\t3\n")
        try:
            locus_rows, sample_ids = input_tables.read_repeat_copy_numbers(handle.name)
        finally:
            os.unlink(handle.name)
        self.assertEqual(sample_ids, ["S1"])
        self.assertEqual(locus_rows[0]["genotypes"], {"S1": "3"})

    def test_missing_trid_raises(self):
        path = _write_tsv([["AT", "12"]], ["Motif", "S1"])
        try:
            with self.assertRaises(ValueError):
                input_tables.read_repeat_copy_numbers(path)
        finally:
            os.unlink(path)


class ReadSampleMetadataTests(unittest.TestCase):

    def _read(self, rows, header):
        path = _write_tsv(rows, header)
        try:
            return input_tables.read_sample_metadata(path)
        finally:
            os.unlink(path)

    def test_case_insensitive_columns_and_minimal(self):
        sample_lookup, affected, analysis, df = self._read(
            [["S1"], ["S2"]], ["Sample_ID"]
        )
        self.assertIn("S1", sample_lookup)
        self.assertEqual(set(df.columns), {"sample_id"})

    def test_analysis_status_remap(self):
        _, _, analysis, _ = self._read(
            [["S1", "rncc"], ["S2", "rcpc"], ["S3", "s_kgfp"], ["S4", ""]],
            ["sample_id", "analysis_status"],
        )
        self.assertEqual(analysis["S1"], "unsolved")
        self.assertEqual(analysis["S2"], "unsolved")
        self.assertEqual(analysis["S3"], "solved")
        self.assertEqual(analysis["S4"], "unknown")

    def test_affected_status_possibly_affected_maps_to_affected(self):
        _, affected, _, _ = self._read(
            [["S1", "Possibly Affected"], ["S2", "Unaffected"]],
            ["sample_id", "affected_status"],
        )
        self.assertEqual(affected["S1"], "affected")
        self.assertEqual(affected["S2"], "unaffected")

    def test_sample_id_with_colon_rejected(self):
        with self.assertRaises(ValueError):
            self._read([["S:1"]], ["sample_id"])

    def test_dedup_keep_first(self):
        sample_lookup, _, _, df = self._read(
            [["S1", "male"], ["S1", "female"]],
            ["sample_id", "sex"],
        )
        self.assertEqual(len(df), 1)
        self.assertEqual(sample_lookup["S1"]["sex"], "male")

    def test_phenotype_na_prefix_stripped(self):
        sample_lookup, _, _, _ = self._read(
            [["S1", "NA; seizures"]],
            ["sample_id", "phenotype_description"],
        )
        self.assertEqual(sample_lookup["S1"]["phenotype_description"], "seizures")

    def test_invalid_sex_raises(self):
        with self.assertRaises(ValueError):
            self._read([["S1", "intersex"]], ["sample_id", "sex"])

    def test_full_metadata_lookup_fields(self):
        sample_lookup, _, _, _ = self._read(
            [["S1", "male", "F1", "M1", "P1", "ataxia", "affected", "unsolved"]],
            ["sample_id", "sex", "family_id", "maternal_id", "paternal_id",
             "phenotype_description", "affected_status", "analysis_status"],
        )
        row = sample_lookup["S1"]
        self.assertEqual(row["family_id"], "F1")
        self.assertEqual(row["maternal_id"], "M1")
        self.assertEqual(row["affected_status"], "affected")


class ReadPhenotypesTests(unittest.TestCase):

    def test_none_returns_empty(self):
        self.assertEqual(input_tables.read_phenotypes(None), {})

    def test_groups_hpo_terms(self):
        path = _write_tsv(
            [["S1", "HP:0001250"], ["S1", "HP:0000252"], ["S2", "not_hpo"]],
            ["participant_id", "term_id"],
        )
        try:
            result = input_tables.read_phenotypes(path)
        finally:
            os.unlink(path)
        self.assertEqual(result["S1"], {"HP:0001250", "HP:0000252"})
        self.assertNotIn("S2", result)

    def test_case_insensitive_header(self):
        path = _write_tsv(
            [["S1", "HP:0001250"]],
            ["Participant_ID", "Term_ID"],
        )
        try:
            result = input_tables.read_phenotypes(path)
        finally:
            os.unlink(path)
        self.assertEqual(result["S1"], {"HP:0001250"})


class ReadGeneTableTests(unittest.TestCase):

    def test_none_returns_empty(self):
        self.assertEqual(input_tables.read_gene_table(None), {})

    def test_keeps_known_columns(self):
        path = _write_tsv(
            [["ENSG1", "FMR1", "0.9", "0.95", "AD"]],
            ["gene_id", "gene_symbol", "pLI_v2", "pLI_v4", "inheritance"],
        )
        try:
            gene_lookup = input_tables.read_gene_table(path)
        finally:
            os.unlink(path)
        self.assertEqual(gene_lookup["ENSG1"]["gene_symbol"], "FMR1")
        self.assertEqual(gene_lookup["ENSG1"]["inheritance"], "AD")


class ReadGeneDiseasePhenotypesTests(unittest.TestCase):

    def test_none_returns_empty(self):
        self.assertEqual(input_tables.read_gene_disease_phenotypes(None), {})

    def test_splits_phenotype_and_inheritance(self):
        rows = [
            ["1", "FMR1", "HP:0001250", "x", "y", "OMIM:1"],
            ["1", "FMR1", "HP:0000006", "x", "y", "OMIM:1"],  # AD inheritance term
        ]
        path = _write_tsv(rows, ["ncbi", "gene", "hpo", "name", "freq", "disease"])
        try:
            data = input_tables.read_gene_disease_phenotypes(path)
        finally:
            os.unlink(path)
        self.assertEqual(data["FMR1"]["OMIM:1"]["phenotypes"], {"HP:0001250"})
        self.assertEqual(data["FMR1"]["OMIM:1"]["inheritance"], {"AD"})


if __name__ == "__main__":
    unittest.main()
