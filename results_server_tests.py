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

    # The reference-region filter's backing index (see result_database.create_loci_indexes);
    # named through the server's constant so a rename cannot silently split the two.
    conn.execute(f"CREATE INDEX {results_server.REFERENCE_REGION_INDEX} "
                 f"ON loci(Chrom, Start0Based, End1Based)")

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

    # -- Reference region filter -------------------------------------------------

    def test_parse_reference_region_formats(self):
        for text, expected in [
            ("chr16:11579459-11579529", ("chr16", 11579459, 11579529)),
            ("chr16:11,579,459-11,579,529", ("chr16", 11579459, 11579529)),
            ("16:100-200", ("chr16", 100, 200)),
            ("  chrX:0-5  ", ("chrX", 0, 5)),
            # The prefix is rebuilt in lowercase whatever case was typed.
            ("Chr1:100-200", ("chr1", 100, 200)),
            ("CHR1:100-200", ("chr1", 100, 200)),
            # A bare position, and a zero-length interval, both mean that single base.
            ("chr1:100", ("chr1", 100, 101)),
            ("chr1:100-100", ("chr1", 100, 101)),
            # A chromosome on its own means the whole chromosome.
            ("chr1", ("chr1", 0, None)),
            ("chr1:", ("chr1", 0, None)),
        ]:
            region, error = results_server.parse_reference_region(text)
            self.assertIsNone(error, text)
            self.assertEqual(region, expected, text)

    def test_parse_reference_region_rejects_invalid(self):
        for text in ["chr1:200-100", "chr1:abc", "chr1:100-abc", "chr1:-200", "chr1:100-200-300", "ch r1:1-2",
                     # isdigit() accepts these but int() does not, and SQLite cannot bind
                     # a coordinate wider than a signed 64-bit integer.
                     "chr1:²", "chr1:9223372036854775808",
                     "chr1:99999999999999999999-99999999999999999999999",
                     # A dash means an end is coming; without one the region is truncated,
                     # not a bare position.
                     "chr1:100-", "chr1:-"]:
            region, error = results_server.parse_reference_region(text)
            self.assertIsNone(region, text)
            self.assertTrue(error, text)
        # Blank input is not an error, it just means "no region filter".
        self.assertEqual(results_server.parse_reference_region("  "), (None, None))

    def test_reference_region_filters_to_overlapping_loci(self):
        # loci: chr1:[100,110) and chr2:[200,209).
        for region, expected_ids in [
            ("chr1:105-106", {"chr1-100-110-AT"}),
            ("1:105-106", {"chr1-100-110-AT"}),
            # Case variants of the chromosome name still match the stored "chr1".
            ("Chr1:105-106", {"chr1-100-110-AT"}),
            ("CHR1:105-106", {"chr1-100-110-AT"}),
            ("chr1:105", {"chr1-100-110-AT"}),
            ("chr1:0-100000", {"chr1-100-110-AT"}),
            ("chr2", {"chr2-200-209-AAG"}),
            # Half-open: a region abutting the locus on either side does not overlap it.
            ("chr1:110-120", set()),
            ("chr1:90-100", set()),
            ("chr3:1-1000", set()),
        ]:
            response = self.client.get(f"/api/v1/loci?outlier_type=all&reference_region={region}")
            self.assertEqual(response.status_code, 200, region)
            data = response.get_json()
            self.assertEqual({r["LocusId"] for r in data["results"]}, expected_ids, region)
            self.assertEqual(data["total"], len(expected_ids), region)
            # The user-facing string is reported back; the parsed tuple stays internal.
            self.assertEqual(data["filters_applied"].get("reference_region"), region)
            self.assertNotIn("reference_region_parsed", data["filters_applied"])

    def test_reference_region_invalid_returns_400(self):
        # An out-of-range coordinate must be rejected up front, not surface as a 500 from
        # SQLite refusing to bind it.
        for region in ["chr1:200-100", "chr1:9223372036854775808"]:
            response = self.client.get(f"/api/v1/loci?outlier_type=all&reference_region={region}")
            self.assertEqual(response.status_code, 400, region)
            self.assertIn("reference_region", json.dumps(response.get_json()), region)

    def test_chromosome_name_variants(self):
        # Both naming conventions, both suffix cases, and the M/MT synonym.
        self.assertEqual(results_server.chromosome_name_variants("chr1"), ["1", "chr1"])
        self.assertEqual(results_server.chromosome_name_variants("chrx"),
                         ["X", "chrX", "chrx", "x"])
        self.assertEqual(results_server.chromosome_name_variants("chrMT"),
                         ["M", "MT", "chrM", "chrMT"])
        self.assertEqual(results_server.chromosome_name_variants("chrM"),
                         ["M", "MT", "chrM", "chrMT"])

    def test_reference_region_advertised_in_schema(self):
        filters = self.client.get("/api/v1/schema").get_json()["filters"]
        self.assertIn("reference_region", filters)

    def test_reference_region_export(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=tsv&reference_region=chr1:105-106")
        self.assertEqual(response.status_code, 200)
        body = gzip.decompress(response.get_data()).decode("utf-8")
        self.assertIn("chr1-100-110-AT", body)
        self.assertNotIn("chr2-200-209-AAG", body)

    def test_reference_region_query_bypasses_skinny_table(self):
        # Start0Based / End1Based are not projected into the skinny tables, so a region
        # query has to run against loci even when the skinny fast path is available.
        base = {"outlier_type": "all", "page": 1, "page_size": 50}
        original = results_server.app.config.get("HAS_SKINNY", False)
        results_server.app.config["HAS_SKINNY"] = True
        try:
            select_query = results_server.build_api_query(base)[0]
            self.assertIn("FROM sk_AllAlleles AS loci", select_query)

            region_params = dict(base, reference_region="chr1:105-106",
                                 reference_region_parsed=("chr1", 105, 106))
            select_query = results_server.build_api_query(region_params)[0]
            self.assertNotIn("sk_AllAlleles", select_query)
            self.assertIn("FROM loci AS loci", select_query)
        finally:
            results_server.app.config["HAS_SKINNY"] = original

    def test_reference_region_uses_the_composite_index(self):
        # The whole point of the filter's shape: an equality seek on Chrom followed by a
        # b-tree range scan over Start0Based, rather than a scan of every locus.
        _, count_query, _, _, sql_params, _ = results_server.build_api_query(
            {"outlier_type": "all", "page": 1, "page_size": 50,
             "reference_region_parsed": ("chr1", 105, 106)})
        conn = sqlite3.connect(self.db_path)
        try:
            plan = " ".join(str(row) for row in conn.execute(
                "EXPLAIN QUERY PLAN " + count_query, sql_params))
        finally:
            conn.close()
        self.assertIn(results_server.REFERENCE_REGION_INDEX, plan)


if __name__ == "__main__":
    unittest.main()
