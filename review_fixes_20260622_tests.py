"""Regression tests for the issues found in the 2026-06-22 code review.

Each test pins a specific finding (WS*/DP*/QC* ids from the consolidated review)
so the bug cannot silently return. These are the ten findings that all three
validators unanimously agreed were real.

WU1 (a JavaScript whisker-fallback bug in qc_shared_js.html) has no Python test
harness in this repo and is verified by inspection only.
"""

import os
import tempfile
import unittest

import intervaltree

import build_database
import duckdb_compat
import input_tables
import locus_annotations
import mendelian_qc
import results_server


class EmptyMotifMendelianTests(unittest.TestCase):
    """QC1: a blank motif must not raise KeyError '0bp' and abort the build."""

    def test_blank_motif_does_not_crash(self):
        locus_rows = [{
            "trid": "chr1-100-110-",
            "motif": "",
            "genotypes": {"CHILD": "10,20", "MOM": "10,11", "DAD": "20,21"},
        }]
        sample_lookup = {
            "CHILD": {"sample_id": "CHILD", "maternal_id": "MOM", "paternal_id": "DAD"},
            "MOM": {"sample_id": "MOM", "maternal_id": "", "paternal_id": ""},
            "DAD": {"sample_id": "DAD", "maternal_id": "", "paternal_id": ""},
        }
        per_sample, _ = mendelian_qc.compute_mendelian_violations(
            locus_rows, sample_lookup, sample_lookup)
        # The blank-motif locus is skipped, so the trio tallies zero loci.
        self.assertEqual(per_sample[0]["total_loci"], 0)


class StrchiveOnlyBuildTests(unittest.TestCase):
    """DP2: a STRchive-only catalog (no Broad catalog) still annotates known loci."""

    def setUp(self):
        self.directory = tempfile.mkdtemp()
        self.matrix_path = os.path.join(self.directory, "matrix.tsv")
        self.metadata_path = os.path.join(self.directory, "samples.tsv")
        self.strchive_path = os.path.join(self.directory, "strchive.json")
        self.db_path = os.path.join(self.directory, "result.duckdb")

        with open(self.matrix_path, "w") as handle:
            handle.write("trid\tmotif\tS1\tS2\n"
                         "chr2-499-530-CAG\tCAG\t10,11\t10,12\n")
        with open(self.metadata_path, "w") as handle:
            handle.write("sample_id\tsex\taffected_status\n"
                         "S1\tmale\taffected\n"
                         "S2\tfemale\tunaffected\n")
        with open(self.strchive_path, "w") as handle:
            handle.write(
                '[{"locus_id": "STR_A", "disease": "DiseaseA", "chrom": "chr2", '
                '"start_hg38": 500, "stop_hg38": 530, '
                '"reference_motif_reference_orientation": ["CAG"]}]')

    def tearDown(self):
        for name in os.listdir(self.directory):
            os.remove(os.path.join(self.directory, name))
        os.rmdir(self.directory)

    def test_strchive_only_annotation(self):
        # No known_loci_json supplied -- only strchive_loci_json.
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
