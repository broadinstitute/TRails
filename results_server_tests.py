"""Unit tests for the single-database TRails Flask server (results_server.py).

These tests build a tiny temporary SQLite database with a minimal ``loci`` table
(plus a ``swim_plot`` table) and exercise the server through Flask's test client,
so no real TCP port is bound. They check that:

  * GET /api/v1/loci returns rows from the database,
  * a request carrying the legacy ``?source=`` query parameter does NOT 500
    (the parameter is ignored, not rejected),
  * there is no /readviz route (it returns 404), and
  * the kept pages/endpoints (schema, swim, export, qc) behave gracefully.
"""

import gzip
import json
import os
import sqlite3
import tempfile
import unittest

import results_server


def _build_minimal_db(path):
    """Create a tiny loci + swim_plot database for the AllAlleles outlier type."""
    conn = sqlite3.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT,
        Motif TEXT,
        CanonicalMotif TEXT,
        IsKnownMotif INTEGER,
        IsInMendelianGene INTEGER,
        AllAlleleHistogram TEXT,
        ShortAlleleHistogram TEXT,
        HemizygousAlleleHistogram TEXT,
        OutlierSampleIds_AllAlleles TEXT,
        OutlierSampleIds_ShortAlleles TEXT,
        OutlierSampleIds_HemizygousAlleles TEXT,
        Source TEXT,
        Chrom TEXT,
        Start0Based INTEGER,
        End1Based INTEGER,
        ReferenceRegion TEXT,
        KnownDiseaseLocus TEXT,
        MotifSize INTEGER,
        gene_id TEXT,
        gene_region TEXT,
        gene_region_rank INTEGER,
        NumRepeatsInReference INTEGER,
        pLI REAL,
        inheritance TEXT,
        FirstAffectedAlleleSize_AllAlleles INTEGER,
        SecondAffectedAlleleSize_AllAlleles INTEGER,
        ThirdAffectedAlleleSize_AllAlleles INTEGER,
        FirstUnaffectedAlleleSize_AllAlleles INTEGER,
        NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER,
        NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles INTEGER,
        FirstAffectedSampleId_AllAlleles TEXT,
        FirstAffectedPhenotype_AllAlleles TEXT,
        MaxGenePhenoSim_AllAlleles REAL,
        SumPairwiseSim_AllAlleles REAL,
        GeneTableGeneSymbol TEXT,
        GeneTableInheritance TEXT,
        GeneTableLLMPhenotypeSummary TEXT
    )""")
    rows = [
        ("chr1-100-110-AT", "AT", "AT", 1, 0, "5x:1,6x:2", "5x:1", None,
         "6x:sampleA,5x:sampleB", "5x:sampleB", None, "TRails", "chr1", 100, 110,
         "chr1:100-110", "", 2, "GENE1", "CDS", 1, 5, 0.9, "AD",
         6, 5, None, 4, 2, 1, "sampleA", "seizures", 0.5, 1.2,
         "GENE1", "AD", "summary"),
        ("chr2-200-209-AAG", "AAG", "AAG", 0, 1, "3x:4", "3x:4", None,
         "4x:sampleC", "4x:sampleC", None, "TRails", "chr2", 200, 209,
         "chr2:200-209", "", 3, "GENE2", "intron", 5, 3, 0.1, "AR",
         4, None, None, None, 1, 0, "sampleC", "ataxia", 0.0, 0.0,
         "GENE2", "AR", "summary2"),
    ]
    conn.executemany(
        "INSERT INTO loci VALUES (" + ",".join("?" * 37) + ")", rows
    )

    conn.execute("""CREATE TABLE swim_plot (
        outlier_type TEXT,
        outlier_rank INTEGER,
        motif_category TEXT,
        allele_size INTEGER,
        sample_id TEXT,
        family_id TEXT,
        affected_status TEXT,
        analysis_status TEXT,
        sex TEXT,
        phenotype_description TEXT,
        purity TEXT,
        methylation TEXT,
        FirstUnaffectedAlleleSize INTEGER,
        is_above_first_unaffected INTEGER,
        LocusId TEXT,
        Motif TEXT,
        CanonicalMotif TEXT,
        MotifSize INTEGER,
        gene_region TEXT,
        GeneTableGeneSymbol TEXT
    )""")
    conn.executemany(
        "INSERT INTO swim_plot VALUES (" + ",".join("?" * 20) + ")",
        [
            ("AllAlleles", 1, "2bp", 6, "sampleA", "famA", "Affected", "Unsolved",
             "male", "seizures", None, None, 4, 1, "chr1-100-110-AT", "AT", "AT", 2, "CDS", "GENE1"),
            ("AllAlleles", 2, "2bp", 5, "sampleB", "famB", "Unaffected", "Unaffected",
             "female", None, None, None, 4, 0, "chr1-100-110-AT", "AT", "AT", 2, "CDS", "GENE1"),
            ("AllAlleles", 1, "3bp", 4, "sampleC", "famC", "Affected", "Unsolved",
             "male", "ataxia", None, None, None, 1, "chr2-200-209-AAG", "AAG", "AAG", 3, "intron", "GENE2"),
        ],
    )
    conn.commit()
    conn.close()


class ResultsServerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.db")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.db")
        _build_minimal_db(cls.db_path)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    def test_loci_returns_rows(self):
        response = self.client.get("/api/v1/loci?outlier_type=all")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["total"], 2)
        self.assertEqual(len(data["results"]), 2)
        locus_ids = {r["LocusId"] for r in data["results"]}
        self.assertIn("chr1-100-110-AT", locus_ids)
        self.assertIn("chr2-200-209-AAG", locus_ids)

    def test_source_param_is_ignored_not_500(self):
        response = self.client.get("/api/v1/loci?outlier_type=all&source=foo")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["total"], 2)
        # The ignored source must not leak into filters_applied.
        self.assertNotIn("source", data["filters_applied"])

    def test_no_readviz_route(self):
        for url in ("/readviz/chr1-100-110-AT",
                    "/api/v1/readviz/chr1-100-110-AT",
                    "/api/v1/readviz/request",
                    "/api/v1/deepdive/request"):
            response = self.client.get(url)
            self.assertEqual(response.status_code, 404, f"{url} should not exist")

    def test_outlier_type_required(self):
        response = self.client.get("/api/v1/loci")
        self.assertEqual(response.status_code, 400)

    def test_locus_detail(self):
        response = self.client.get("/api/v1/loci/chr1-100-110-AT")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["locus"]["LocusId"], "chr1-100-110-AT")
        self.assertIn("AllAlleles", data["outlier_samples"])
        self.assertEqual(data["system_tags"], [])

    def test_schema_has_single_source(self):
        response = self.client.get("/api/v1/schema")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["total_loci"], 2)
        self.assertEqual(data["available_sources"], ["TRails"])
        # source must NOT be advertised as a filter in the single-DB server.
        self.assertNotIn("source", data["filters"])

    def test_swim_plot_data(self):
        response = self.client.get("/api/v1/swim_plot_data?outlier_type=all")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        sample_ids = {entry["sample_id"] for entry in data["data"]}
        self.assertIn("sampleA", sample_ids)

    def test_export_tsv(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=tsv")
        self.assertEqual(response.status_code, 200)
        self.assertIn(".tsv.gz", response.headers["Content-Disposition"])
        body = gzip.decompress(response.get_data()).decode("utf-8")
        self.assertIn("LocusId", body.splitlines()[0])
        self.assertIn("chr1-100-110-AT", body)

    def test_export_bed(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=bed")
        self.assertEqual(response.status_code, 200)
        body = response.get_data(as_text=True)
        self.assertIn("chr1\t100\t110\tchr1-100-110-AT", body)

    def test_export_json_nested_metadata_and_samples(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=json&min_pli=0.5")
        self.assertEqual(response.status_code, 200)
        self.assertIn(".json.gz", response.headers["Content-Disposition"])
        data = json.loads(gzip.decompress(response.get_data()))
        self.assertEqual(set(data.keys()), {"metadata", "loci"})
        self.assertEqual(data["metadata"]["outlier_type"], "all")
        self.assertEqual(data["metadata"]["total_loci"], 1)
        self.assertEqual(data["metadata"]["filters_applied"], {"min_pli": 0.5})
        self.assertNotIn("outlier_type", data["metadata"]["filters_applied"])

        self.assertEqual(len(data["loci"]), 1)
        locus = data["loci"][0]
        self.assertEqual(locus["LocusId"], "chr1-100-110-AT")
        # Every qualifying outlier sample must appear, not just the first one —
        # this is the whole point of the nested Samples list over the old flat
        # First/Second/ThirdAffected* columns.
        sample_ids = {s["Sample ID"] for s in locus["Samples"]["AllAlleleOutliers"]}
        self.assertEqual(sample_ids, {"sampleA", "sampleB"})

    def test_export_json_empty_result_is_well_formed(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=json&gene_id=NO_SUCH_GENE")
        self.assertEqual(response.status_code, 200)
        data = json.loads(gzip.decompress(response.get_data()))
        self.assertEqual(data["metadata"]["total_loci"], 0)
        self.assertEqual(data["loci"], [])

    def test_sample_qc_data(self):
        response = self.client.get("/api/v1/sample_qc_data")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertIn("rank1", data)
        self.assertIn("top10", data)

    def test_qc2_data_absent_is_graceful(self):
        # The minimal DB has no mendelian_violations table, so the QC2 endpoint
        # returns a controlled 500 (not an uncaught exception) and the page 404s.
        response = self.client.get("/api/v1/qc2_data")
        self.assertEqual(response.status_code, 500)
        self.assertIn("mendelian", response.get_json()["error"].lower())
        page = self.client.get("/qc2")
        self.assertEqual(page.status_code, 404)

    def test_sample_outlier_stats(self):
        response = self.client.get(
            "/api/v1/sample_outlier_stats?outlier_type=all&sample_id=sampleA&locus_id=chr1-100-110-AT"
        )
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["total_loci_above_unaffected"], 1)

    def test_annotation_tag_roundtrip(self):
        add = self.client.post(
            "/api/v1/annotations/chr1-100-110-AT/tags",
            json={"tag": "interesting"},
        )
        self.assertEqual(add.status_code, 200)
        self.assertIn("interesting", add.get_json()["tags"])
        # Filtering by the tag should now return exactly that locus.
        filtered = self.client.get("/api/v1/loci?outlier_type=all&tag=interesting")
        self.assertEqual(filtered.status_code, 200)
        self.assertEqual(filtered.get_json()["total"], 1)

    def test_index_page_renders(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)
        # The readviz request route is gone, so its handler script must not appear.
        self.assertNotIn("/readviz/request", response.get_data(as_text=True))


if __name__ == "__main__":
    unittest.main()
