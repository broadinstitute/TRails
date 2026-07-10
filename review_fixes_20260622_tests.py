"""Regression tests for the issues found in the 2026-06-22 code review.

Each test pins a specific finding (WS*/DP*/QC* ids from the consolidated review)
so the bug cannot silently return. These are the ten findings that all three
validators unanimously agreed were real.

WU1 (a JavaScript whisker-fallback bug in qc_shared_js.html) has no Python test
harness in this repo and is verified by inspection only.
"""

import os
import sqlite3
import tempfile
import unittest

import intervaltree

import build_database
import input_tables
import locus_annotations
import mendelian_qc
import results_server


def _write(text, suffix=".tsv"):
    handle = tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False)
    handle.write(text)
    handle.close()
    return handle.name


class CoerceAnnotationValueTests(unittest.TestCase):
    """DP1: '.'/'NA'-style missing markers must not survive as strings."""

    def test_missing_markers_become_none(self):
        for marker in ("", ".", "NA", "N/A", "NaN", "nan", "none", "NULL", " . "):
            self.assertIsNone(
                input_tables._coerce_annotation_value(marker),
                msg=f"marker {marker!r} should coerce to None")

    def test_numeric_and_string_values_preserved(self):
        self.assertEqual(input_tables._coerce_annotation_value("30"), 30)
        self.assertEqual(input_tables._coerce_annotation_value("3.5"), 3.5)
        self.assertEqual(input_tables._coerce_annotation_value("CAG"), "CAG")

    def test_dot_in_population_column_does_not_survive(self):
        path = _write(
            "trid\tmotif\tHPRC256_99thPercentile\tS1\n"
            "chr1-100-115-CAG\tCAG\t.\t12,40\n")
        self.addCleanup(os.remove, path)
        locus_rows, _ = input_tables.read_repeat_copy_numbers(path)
        # A '.' population-stat cell must not reach the numeric downstream as a str.
        self.assertIsNone(locus_rows[0]["extra_columns"]["HPRC256_99thPercentile"])


class InheritanceHPOTermsTests(unittest.TestCase):
    """DP3: non-Mendelian inheritance HPO terms must not enter the inheritance set."""

    def test_unmapped_inheritance_term_dropped(self):
        path = _write(
            "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n"
            "1\tADD1\tHP:0001426\tNon-Mendelian\t-\tOMIM:145500\n"   # unmapped -> dropped
            "1\tADD1\tHP:0000006\tAD\t-\tOMIM:145500\n"             # mapped -> AD
            "1\tADD1\tHP:0004972\tphenotype\t-\tOMIM:145500\n",     # phenotype
            suffix=".txt")
        self.addCleanup(os.remove, path)
        data = input_tables.read_gene_disease_phenotypes(path)
        entry = data["ADD1"]["OMIM:145500"]
        self.assertEqual(entry["inheritance"], {"AD"})
        self.assertNotIn("HP:0001426", entry["inheritance"])
        self.assertEqual(entry["phenotypes"], {"HP:0004972"})


class AnnotationColumnPrefixTests(unittest.TestCase):
    """DP4: a sample column named like a cohort must not be swallowed as annotation."""

    def test_bare_cohort_sample_ids_are_samples(self):
        self.assertIsNone(input_tables._is_annotation_column("AoU_0001"))
        self.assertIsNone(input_tables._is_annotation_column("HG002"))

    def test_real_cohort_stat_columns_still_recognized(self):
        for column in ("HPRC256_99thPercentile", "AoU1027_MaxAllele",
                       "AoUPhase2HighCov_N", "TenK10K_MaxAllele", "TRExplorerMotif"):
            self.assertEqual(input_tables._is_annotation_column(column), column)

    def test_aou_sample_column_kept_as_sample(self):
        path = _write(
            "trid\tmotif\tAoU_0001\tHG002\n"
            "chr1-100-115-CAG\tCAG\t5,60\t5,6\n")
        self.addCleanup(os.remove, path)
        _, sample_id_list = input_tables.read_repeat_copy_numbers(path)
        self.assertEqual(sample_id_list, ["AoU_0001", "HG002"])


class PathogenicMotifsNullTests(unittest.TestCase):
    """WS2: PathogenicMotifs explicitly null must not raise TypeError."""

    def test_null_pathogenic_motifs_does_not_crash(self):
        tree = intervaltree.IntervalTree()
        tree.addi(100, 110, data={
            "RepeatUnit": "AT",
            "PathogenicMotifs": None,   # explicit null in the catalog JSON
            "LocusId": "DISEASE_AT",
            "Diseases": [{"Name": "DiseaseX"}],
        })
        lookups = {"locus": {}, "disease_trees": {"1": tree}, "strchive_trees": {}}
        row = {"LocusId": "q", "Chrom": "chr1", "Start0Based": 100,
               "End1Based": 110, "Motif": "AT"}
        result = results_server.compute_known_disease_info(row, lookups)
        self.assertIsNotNone(result)
        self.assertEqual(result["locus_id"], "DISEASE_AT")


class AffectedStatusDisplayTests(unittest.TestCase):
    """WS5: 'possibly affected' must be preserved for display, not collapsed."""

    def test_possibly_affected_preserved(self):
        sample_lookup = {"sA": {"affected_status": "possibly affected",
                                "family_id": "f1", "sex": "male"}}
        affected_lookup = {"sA": "affected"}  # collapsed for analysis logic only
        analysis_lookup = {"sA": "unsolved"}
        parsed = results_server.parse_outlier_samples(
            "60x:sA", sample_lookup, affected_lookup, analysis_lookup)
        self.assertEqual(parsed[0]["affected_status"], "Possibly Affected")

    def test_falls_back_to_lookup_when_row_missing(self):
        parsed = results_server.parse_outlier_samples(
            "60x:sB", {}, {"sB": "affected"}, {"sB": "unsolved"})
        self.assertEqual(parsed[0]["affected_status"], "Affected")


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
        self.db_path = os.path.join(self.directory, "result.db")

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
        connection = sqlite3.connect(self.db_path)
        try:
            value = connection.execute(
                "SELECT KnownDiseaseLocus FROM loci WHERE LocusId = 'chr2-499-530-CAG'"
            ).fetchone()[0]
        finally:
            connection.close()
        self.assertEqual(value, "STR_A")


def _build_db_with_unknown_category(path):
    """Minimal loci + swim_plot DB with a normal row and an 'Unknown' motif row."""
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE loci (LocusId TEXT, Motif TEXT, KnownDiseaseLocus TEXT)")
    conn.executemany("INSERT INTO loci VALUES (?, ?, ?)", [
        ("chr1-100-110-AT", "AT", ""),
        ("chrX-1-9-CTTTTT", "CTTTTT", ""),
    ])
    conn.execute("""CREATE TABLE swim_plot (
        outlier_type TEXT, motif_category TEXT, allele_size INTEGER, sample_id TEXT,
        family_id TEXT, affected_status TEXT, analysis_status TEXT, sex TEXT,
        phenotype_description TEXT, purity TEXT, methylation TEXT, LocusId TEXT,
        Motif TEXT, gene_region TEXT, GeneTableGeneSymbol TEXT)""")
    conn.executemany(
        "INSERT INTO swim_plot (outlier_type, motif_category, allele_size, sample_id, "
        "LocusId, Motif) VALUES (?, ?, ?, ?, ?, ?)",
        [
            ("AllAlleles", "2bp", 6, "sampleA", "chr1-100-110-AT", "AT"),
            ("AllAlleles", "Unknown", 99, "sampleZ", "chrX-1-9-CTTTTT", "CTTTTT"),
        ])
    conn.commit()
    conn.close()


class SchemaAndUnknownCategoryTests(unittest.TestCase):
    """WS3 + WS4: schema sample_ids fall back to swim_plot; Unknown rows served."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.db")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.db")
        _build_db_with_unknown_category(cls.db_path)
        # No sample_table -> the sample dropdown must fall back to swim_plot IDs.
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def test_schema_sample_ids_from_swim_plot(self):
        data = self.client.get("/api/v1/schema").get_json()
        self.assertEqual(data["sample_ids"], ["sampleA", "sampleZ"])

    def test_unknown_motif_category_returned(self):
        data = self.client.get("/api/v1/swim_plot_data?outlier_type=all").get_json()
        self.assertIn("Unknown", data["categories"])
        self.assertIn("sampleZ", {entry["sample_id"] for entry in data["data"]})


if __name__ == "__main__":
    unittest.main()
