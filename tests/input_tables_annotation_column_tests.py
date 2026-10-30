"""Unit tests for input_tables.read_repeat_copy_numbers: matrix annotation columns are not samples."""

import os
import tempfile
import unittest

import input_tables


class AnnotationColumnClassificationTests(unittest.TestCase):
    """Recognized annotation columns annotate the locus instead of being read as samples."""

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


if __name__ == "__main__":
    unittest.main()
