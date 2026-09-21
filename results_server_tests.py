"""Unit tests for the single-database TRails Flask server (results_server.py).

These tests build a tiny temporary DuckDB database with a minimal ``loci`` table
(plus a ``swim_plot`` table) and exercise the server through Flask's test client,
so no real TCP port is bound. They check that:

  * GET /api/v1/loci returns rows from the database,
  * a request carrying the legacy ``?source=`` query parameter does NOT 500
    (the parameter is ignored, not rejected),
  * there is no /readviz route (it returns 404), and
  * the kept pages/endpoints (schema, swim, export, qc) behave gracefully.
"""

import gzip
import importlib.util
import json
import os
import sys
import tempfile
import threading
import unittest

import duckdb_compat
import locus_annotations
import result_database
import results_server
import swim_plot


def _build_minimal_db(path):
    """Create a tiny loci + swim_plot database for the AllAlleles outlier type."""
    conn = duckdb_compat.connect(path)
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

    # Build the fixture's swim_plot through the same writer the build uses, so it gets the real
    # column list rather than a hand-picked subset and the real per-column types. A filter that
    # reads a column the build actually writes, or compares one of its numeric columns to a
    # number, is then exercised here instead of failing only in a real deployment.
    swim_columns = swim_plot.SWIM_PLOT_COLUMNS

    def _swim_row(**values):
        """Return a full-width swim_plot row dict, defaulting every unnamed column to NULL."""
        unknown = set(values) - set(swim_columns)
        assert not unknown, "swim_plot has no column(s): %s" % sorted(unknown)
        return {column: values.get(column) for column in swim_columns}

    result_database.write_swim_plot(
        conn,
        [
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                      allele_size=6, sample_id="sampleA", family_id="famA",
                      affected_status="Affected", analysis_status="Unsolved", sex="male",
                      phenotype_description="seizures", FirstUnaffectedAlleleSize=4,
                      is_above_first_unaffected=1, LocusId="chr1-100-110-AT", Motif="AT",
                      CanonicalMotif="AT", MotifSize=2, gene_region="CDS",
                      GeneTableGeneSymbol="GENE1"),
            _swim_row(outlier_type="AllAlleles", outlier_rank=2, motif_category="2bp",
                      allele_size=5, sample_id="sampleB", family_id="famB",
                      affected_status="Unaffected", analysis_status="Unaffected", sex="female",
                      FirstUnaffectedAlleleSize=4, is_above_first_unaffected=0,
                      LocusId="chr1-100-110-AT", Motif="AT", CanonicalMotif="AT", MotifSize=2,
                      gene_region="CDS", GeneTableGeneSymbol="GENE1"),
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="3bp",
                      allele_size=4, sample_id="sampleC", family_id="famC",
                      affected_status="Affected", analysis_status="Unsolved", sex="male",
                      phenotype_description="ataxia", is_above_first_unaffected=1,
                      LocusId="chr2-200-209-AAG", Motif="AAG", CanonicalMotif="AAG", MotifSize=3,
                      gene_region="intron", GeneTableGeneSymbol="GENE2"),
        ],
        columns=swim_columns,
    )
    conn.commit()
    conn.close()


class ResultsServerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
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

    def test_schema_total_samples_is_null_without_the_metadata_table(self):
        # This fixture predates the metadata table, so the key is reported as null rather than
        # omitted: the index page checks for null before appending ", N samples" to the count.
        data = self.client.get("/api/v1/schema").get_json()
        self.assertIn("total_samples", data)
        self.assertIsNone(data["total_samples"])

    def test_variation_cluster_sort_is_accepted(self):
        # The sort key the "Variation Cluster" checkbox sends. This fixture has no
        # VariationClusterSizeDiff column, so the key drops out of the ORDER BY instead of
        # producing SQL that names a column the database does not have.
        self.assertIn("variation_cluster",
                      self.client.get("/api/v1/schema").get_json()["sort_options"])
        response = self.client.get("/api/v1/loci?outlier_type=all&sort_by=variation_cluster")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_json()["total"], 2)

    def test_variation_cluster_exclusions_without_the_column_exclude_nothing(self):
        # This fixture has no VariationClusterFilterReason column, so no locus was ever recorded
        # as filtered out and these exclusion boxes have nothing to exclude. Failing closed here
        # would empty the table the moment either box is ticked.
        for query in ("exclude_vc_depth_filtered=true",
                      "exclude_vc_size_filtered=true",
                      "exclude_vc_depth_filtered=true&exclude_vc_size_filtered=true"):
            response = self.client.get(f"/api/v1/loci?outlier_type=all&{query}")
            self.assertEqual(response.status_code, 200, query)
            self.assertEqual(response.get_json()["total"], 2, query)

    def test_swim_plot_data_with_variation_cluster_exclusions_without_the_column(self):
        # The swim-plot endpoint has no copy of these two flags (it ignores them), so it must
        # keep returning its data rather than going empty alongside the results table.
        response = self.client.get("/api/v1/swim_plot_data?outlier_type=all"
                                   "&exclude_vc_depth_filtered=true&exclude_vc_size_filtered=true")
        self.assertEqual(response.status_code, 200)
        sample_ids = {entry["sample_id"] for entry in response.get_json()["data"]}
        self.assertIn("sampleA", sample_ids)

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

    def test_export_expansion_hunter(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=expansion_hunter")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertIn("expansion_hunter.json", response.headers["Content-Disposition"])
        catalog = json.loads(response.get_data())
        self.assertEqual(len(catalog), 2)
        by_locus = {entry["LocusId"]: entry for entry in catalog}
        entry = by_locus["chr1-100-110-AT"]
        # The three fields ExpansionHunter reads, built from the loci row's coordinates and motif.
        self.assertEqual(entry["ReferenceRegion"], "chr1:100-110")
        self.assertEqual(entry["LocusStructure"], "(AT)*")
        self.assertEqual(entry["VariantType"], "Repeat")
        self.assertEqual(by_locus["chr2-200-209-AAG"]["ReferenceRegion"], "chr2:200-209")
        self.assertEqual(by_locus["chr2-200-209-AAG"]["LocusStructure"], "(AAG)*")
        # No sample table was supplied, so no outlier is known to be affected and no locus gets
        # a read-visualization threshold.
        self.assertNotIn("PlotReadVisualization", entry)

    def test_export_expansion_hunter_read_visualization_thresholds(self):
        # chr1's largest outlier (6x:sampleA) is affected and sits above the largest unaffected
        # outlier (5x:sampleB), so it earns a LongAllele threshold. The ShortAlleles outlier list
        # holds only the unaffected sampleB, so there is no ShortAllele threshold.
        lookups = results_server.app.config["LOOKUPS"]
        original = lookups["affected"]
        lookups["affected"] = {"sampleA": "affected", "sampleB": "unaffected"}
        try:
            response = self.client.get("/api/v1/export?outlier_type=all&format=expansion_hunter")
            self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
            catalog = {entry["LocusId"]: entry for entry in json.loads(response.get_data())}
        finally:
            lookups["affected"] = original
        self.assertEqual(catalog["chr1-100-110-AT"]["PlotReadVisualization"],
                         [{"If": "LongAllele", "Is": ">=", "Threshold": 6}])
        self.assertNotIn("PlotReadVisualization", catalog["chr2-200-209-AAG"])

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

    def test_concurrent_annotation_writes_all_succeed(self):
        # Flask serves with threaded=True and DuckDB raises instead of waiting when two
        # connections write the same row at once, so the annotation writes are serialized.
        # Without that, one of each overlapping pair came back as a 500.
        barrier = threading.Barrier(6)
        statuses = []
        lock = threading.Lock()

        def save_note(index):
            client = results_server.app.test_client()
            barrier.wait()
            response = client.put("/api/v1/annotations/chr1-100-110-AT/note",
                                  json={"note_text": f"note {index}"})
            with lock:
                statuses.append(response.status_code)

        threads = [threading.Thread(target=save_note, args=(i,)) for i in range(6)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(statuses, [200] * 6)
        # One of the six notes won, and it is the one both the cache and the database hold.
        stored = self.client.get("/api/v1/annotations/chr1-100-110-AT").get_json()["note"]
        reloaded = results_server.load_annotations(self.annotations_db)
        self.assertEqual(reloaded["notes"]["chr1-100-110-AT"]["note_text"], stored["note_text"])
        self.client.delete("/api/v1/annotations/chr1-100-110-AT/note")

    def test_a_failed_annotation_write_reports_the_error(self):
        # A write that cannot reach the database has to come back as a described error rather
        # than an uncaught exception, and must not leave the in-memory annotations claiming it
        # was saved.
        original = results_server.app.config["ANNOTATIONS_DB_PATH"]
        results_server.app.config["ANNOTATIONS_DB_PATH"] = "/nonexistent_dir/annotations.duckdb"
        try:
            response = self.client.put("/api/v1/annotations/chr2-200-209-AAG/note",
                                       json={"note_text": "will not save"})
        finally:
            results_server.app.config["ANNOTATIONS_DB_PATH"] = original
        self.assertEqual(response.status_code, 500)
        body = response.get_json()
        self.assertIn("note", body["error"].lower())
        self.assertTrue(body["detail"])
        self.assertNotIn("chr2-200-209-AAG", results_server.app.config["ANNOTATIONS"]["notes"])

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
                     # isdigit() accepts these but int() does not, and DuckDB cannot bind
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
        # DuckDB refusing to bind it.
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

    # -- Sample keyword filter ---------------------------------------------------

    def _sample_id_like_ids(self, keyword):
        response = self.client.get(f"/api/v1/loci?outlier_type=all&sample_id_like={keyword}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {r["LocusId"] for r in response.get_json()["results"]}

    def test_sample_id_like_does_not_match_allele_sizes(self):
        # OutlierSampleIds_AllAlleles is "6x:sampleA,5x:sampleB" at chr1 and "4x:sampleC" at chr2,
        # so a numeric keyword used to match the allele size in front of every sample id.
        self.assertEqual(self._sample_id_like_ids("6"), set())
        self.assertEqual(self._sample_id_like_ids("4"), set())
        # The "x:" delimiter is not part of a sample id either.
        self.assertEqual(self._sample_id_like_ids("x:"), set())

    def test_sample_id_like_still_matches_sample_ids(self):
        self.assertEqual(self._sample_id_like_ids("sampleA"), {"chr1-100-110-AT"})
        self.assertEqual(self._sample_id_like_ids("ampleC"), {"chr2-200-209-AAG"})
        self.assertEqual(self._sample_id_like_ids("sample"),
                         {"chr1-100-110-AT", "chr2-200-209-AAG"})

    def test_sample_id_like_treats_underscore_literally(self):
        # "_" is ILIKE's single-character wildcard, so without escaping this would match "sampleA".
        self.assertEqual(self._sample_id_like_ids("sample_A"), set())

    def test_phenotype_keyword_is_answered_from_swim_plot(self):
        # This fixture writes a swim_plot whose phenotype_description column is populated, so the
        # filter is answered from that table and names no AffectedPhenotype column at all. The
        # fallback for a database without swim_plot phenotypes is covered by
        # PhenotypeKeywordFallbackColumnSubsetTests below.
        self.assertTrue(results_server.app.config["HAS_SWIM_PLOT_PHENOTYPES"])
        response = self.client.get("/api/v1/loci?outlier_type=all&phenotype_keyword=seizures")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual({r["LocusId"] for r in response.get_json()["results"]},
                         {"chr1-100-110-AT"})

    # -- Merged search box: locus ids -------------------------------------------

    def test_pasted_locus_ids_become_one_set_lookup(self):
        # A predicate per pasted id costs terms x rows, which made a few hundred ids take
        # seconds. They belong in a single IN list instead, lower-cased on both sides so the
        # match stays case-insensitive. The behavior this replaced is asserted in
        # search_filter_tests.py (case-insensitive matching, a literal underscore); this test
        # pins the plan so the per-term chain cannot come back unnoticed.
        ids = ["chr1-100-110-AT", "CHR2-200-209-AAG", "chr1_alt-1-10-CAG"]
        select_query, _, _, _, sql_params, _ = results_server.build_api_query({
            "outlier_type": "all", "page": 1, "page_size": 50,
            "search_terms": [("locus_id", locus_id) for locus_id in ids],
        })
        self.assertEqual(select_query.count("lower(LocusId) IN ("), 1)
        self.assertIn("lower(LocusId) IN (?,?,?)", select_query)
        self.assertNotIn("LocusId ILIKE", select_query)
        self.assertEqual(sql_params[-3:], [locus_id.lower() for locus_id in ids])

    def test_pasted_locus_ids_keep_their_place_among_the_other_search_terms(self):
        # The ids are hoisted to the front of the OR, so their bind parameters have to move with
        # them: a gene symbol term binding out of order would filter on the wrong value.
        select_query, _, _, _, sql_params, _ = results_server.build_api_query({
            "outlier_type": "all", "page": 1, "page_size": 50,
            "search_terms": [("gene_symbol", "GENE1"), ("locus_id", "chr1-100-110-AT")],
        })
        self.assertIn("(lower(LocusId) IN (?) OR GeneTableGeneSymbol ILIKE ?)", select_query)
        self.assertEqual(sql_params[-2:], ["chr1-100-110-at", "%GENE1%"])

    def test_pasted_locus_ids_still_match_case_insensitively_end_to_end(self):
        response = self.client.get("/api/v1/loci",
                                   query_string={"outlier_type": "all",
                                                 "search": "CHR1-100-110-at,chr2-200-209-AAG"})
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual({r["LocusId"] for r in response.get_json()["results"]},
                         {"chr1-100-110-AT", "chr2-200-209-AAG"})

    def test_pasted_locus_ids_are_one_set_lookup_in_the_swim_plot_endpoint_too(self):
        response = self.client.get("/api/v1/swim_plot_data",
                                   query_string={"outlier_type": "all",
                                                 "search": "CHR1-100-110-at,chr2-200-209-AAG"})
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual({entry["LocusId"] for entry in response.get_json()["data"]},
                         {"chr1-100-110-AT", "chr2-200-209-AAG"})



def _build_db_for_sample_qc(path):
    """swim_plot DB where one sample holds two of a locus's first ten outlier entries."""
    conn = duckdb_compat.connect(path)
    conn.execute("CREATE TABLE loci (LocusId TEXT, Motif TEXT, KnownDiseaseLocus TEXT)")

    def _swim_row(**values):
        """Return a full-width swim_plot row dict, defaulting every unnamed column to NULL."""
        return {column: values.get(column) for column in swim_plot.SWIM_PLOT_COLUMNS}

    result_database.write_swim_plot(
        conn,
        [
            # Both alleles of one diploid sample at one locus, the way accumulate_locus records
            # them: two swim_plot rows for the same (sample, locus).
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                      allele_size=9, sample_id="sampleA", LocusId="chr1-100-110-AT",
                      Motif="AT", MotifSize=2),
            _swim_row(outlier_type="AllAlleles", outlier_rank=2, motif_category="2bp",
                      allele_size=8, sample_id="sampleA", LocusId="chr1-100-110-AT",
                      Motif="AT", MotifSize=2),
            # A second locus in the same motif-size bin, so the expected count is not 1 by default.
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                      allele_size=7, sample_id="sampleA", LocusId="chr2-200-210-AT",
                      Motif="AT", MotifSize=2),
            # Another sample's entry at the first locus, to keep the grouping honest.
            _swim_row(outlier_type="AllAlleles", outlier_rank=3, motif_category="2bp",
                      allele_size=6, sample_id="sampleB", LocusId="chr1-100-110-AT",
                      Motif="AT", MotifSize=2),
        ],
        columns=swim_plot.SWIM_PLOT_COLUMNS,
    )
    conn.commit()
    conn.close()


class SampleQcLocusCountTests(unittest.TestCase):
    """The sample-QC counts are counts of loci, not of swim_plot allele rows."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_for_sample_qc(cls.db_path)

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def test_both_alleles_at_one_locus_count_once(self):
        data = results_server.compute_sample_qc_data(self.db_path)
        top10 = {(entry["sample_id"], entry["bin"]): entry["count"] for entry in data["top10"]}
        # sampleA owns three of the four entries, but they sit at only two loci.
        self.assertEqual(top10[("sampleA", "2bp")], 2)
        self.assertEqual(top10[("sampleB", "2bp")], 1)

    def test_rank1_counts_distinct_loci(self):
        data = results_server.compute_sample_qc_data(self.db_path)
        rank1 = {(entry["sample_id"], entry["bin"]): entry["count"] for entry in data["rank1"]}
        self.assertEqual(rank1[("sampleA", "2bp")], 2)
        # sampleB is rank 3 everywhere, so it has no rank-1 row at all.
        self.assertNotIn(("sampleB", "2bp"), rank1)


STRCHIVE_TEST_LOCI = [
    {
        "id": "TEST_PATHOGENIC",
        "locus_id": "TEST_PATHOGENIC",
        "disease": "Test pathogenic-motif disease",
        # Every record in the shipped STRchive file stores inheritance as a list.
        "inheritance": ["AD"],
        "pathogenic_min": 60,
        "chrom": "chr9",
        # STRchive coordinates are 1-based, so this is the 0-based interval [100, 110).
        "start_hg38": 101,
        "stop_hg38": 110,
        "reference_motif_reference_orientation": ["AAAAG"],
        "pathogenic_motif_reference_orientation": ["AAGGG"],
    },
    {
        "id": "TEST_REFERENCE",
        "locus_id": "TEST_REFERENCE",
        "disease": "Test reference-motif disease",
        # A bare string, the shape an older cached STRchive file may still use.
        "inheritance": "AR",
        "pathogenic_min": 50,
        "chrom": "chr3",
        "start_hg38": 301,
        "stop_hg38": 310,
        "reference_motif_reference_orientation": ["CAG"],
        "pathogenic_motif_reference_orientation": [],
    },
]


def _build_db_for_strchive(path):
    """Two loci: one matching a STRchive pathogenic motif, one matching its reference motif.

    Their largest affected allele sizes straddle the STRchive pathogenic minimums (60 for
    TEST_PATHOGENIC, 50 for TEST_REFERENCE), so the pathogenic-threshold filter keeps the first
    and drops the second.
    """
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER,
        Motif TEXT, MotifSize INTEGER, KnownDiseaseLocus TEXT,
        OutlierSampleIds_AllAlleles TEXT,
        FirstAffectedAlleleSize_AllAlleles INTEGER)""")
    conn.executemany("INSERT INTO loci VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", [
        ("chr9-100-110-AAGGG", "chr9", 100, 110, "AAGGG", 5, "TEST_PATHOGENIC", "60x:sampleA", 60),
        ("chr3-300-310-CAG", "chr3", 300, 310, "CAG", 3, "TEST_REFERENCE", "40x:sampleB", 40),
    ])
    conn.commit()
    conn.close()


class StrchivePathogenicMotifTests(unittest.TestCase):
    """The locus-detail endpoint and the build agree on STRchive pathogenic-motif matches."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        cls.strchive_json = os.path.join(cls.tmpdir, "STRchive-loci.json")
        with open(cls.strchive_json, "w") as output_file:
            json.dump(STRCHIVE_TEST_LOCI, output_file)
        _build_db_for_strchive(cls.db_path)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db,
                                     strchive_loci_json=cls.strchive_json)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def _known_disease_locus(self, locus_id):
        response = self.client.get(f"/api/v1/loci/{locus_id}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return response.get_json()["annotations"]["known_disease_locus"]

    def test_detail_matches_the_build_for_a_pathogenic_motif(self):
        # The build flags this locus as known through the STRchive pathogenic motif; the detail
        # endpoint has to report the same locus rather than "no known disease locus".
        strchive_trees = results_server.app.config["LOOKUPS"]["strchive_trees"]
        self.assertEqual(
            locus_annotations.matches_disease_locus("chr9-100-110-AAGGG", {}, strchive_trees),
            "TEST_PATHOGENIC")
        known = self._known_disease_locus("chr9-100-110-AAGGG")
        self.assertIsNotNone(known)
        self.assertEqual(known["locus_id"], "TEST_PATHOGENIC")
        self.assertEqual(known["source"], "STRchive")
        self.assertEqual(known["pathogenic_min"], 60)

    def test_detail_matches_the_build_for_a_reference_motif(self):
        strchive_trees = results_server.app.config["LOOKUPS"]["strchive_trees"]
        self.assertEqual(
            locus_annotations.matches_disease_locus("chr3-300-310-CAG", {}, strchive_trees),
            "TEST_REFERENCE")
        self.assertEqual(self._known_disease_locus("chr3-300-310-CAG")["locus_id"],
                         "TEST_REFERENCE")

    def test_pathogenic_threshold_map_covers_strchive_only_loci(self):
        # KnownDiseaseLocus is set from the variant catalog OR the STRchive fallback, so the
        # threshold map has to be built from both. With only a STRchive catalog supplied, walking
        # the variant-catalog trees alone would leave it empty.
        self.assertEqual(
            results_server.app.config["KNOWN_DISEASE_LOCUS_THRESHOLDS"],
            {"chr9-100-110-AAGGG": 60, "chr3-300-310-CAG": 50})

    def test_pathogenic_threshold_filter_keeps_a_strchive_only_locus(self):
        # The filter and the locus-detail page must use the same pathogenic minimum: chr9 reaches
        # its 60 and survives, chr3 stops at 40 below its 50 and is dropped.
        self.assertEqual(self._known_disease_locus("chr9-100-110-AAGGG")["pathogenic_min"], 60)
        self.assertEqual(self._known_disease_locus("chr3-300-310-CAG")["pathogenic_min"], 50)
        response = self.client.get("/api/v1/loci?outlier_type=all&apply_pathogenic_threshold=true")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        data = response.get_json()
        self.assertEqual({r["LocusId"] for r in data["results"]}, {"chr9-100-110-AAGGG"})

    def test_strchive_inheritance_is_a_flat_list(self):
        # STRchive stores inheritance as a list, so reporting it must not nest it in another one.
        self.assertEqual(self._known_disease_locus("chr9-100-110-AAGGG")["inheritance"], ["AD"])
        # A bare string from an older cached file still comes back as a list of modes.
        self.assertEqual(self._known_disease_locus("chr3-300-310-CAG")["inheritance"], ["AR"])

    def test_strchive_inheritance_modes_shapes(self):
        self.assertEqual(results_server.strchive_inheritance_modes(["AD", "AR"]), ["AD", "AR"])
        self.assertEqual(results_server.strchive_inheritance_modes("XR"), ["XR"])
        self.assertIsNone(results_server.strchive_inheritance_modes([]))
        self.assertIsNone(results_server.strchive_inheritance_modes(None))


def _build_db_with_all_affected_phenotypes(path):
    """Loci carrying First/Second/Third AffectedPhenotype columns for the AllAlleles type."""
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER,
        Motif TEXT, MotifSize INTEGER, KnownDiseaseLocus TEXT,
        OutlierSampleIds_AllAlleles TEXT,
        FirstAffectedPhenotype_AllAlleles TEXT,
        SecondAffectedPhenotype_AllAlleles TEXT,
        ThirdAffectedPhenotype_AllAlleles TEXT)""")
    conn.executemany("INSERT INTO loci VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", [
        ("chr1-100-110-AT", "chr1", 100, 110, "AT", 2, "", "6x:sampleA,5x:sampleB",
         "seizures", "ataxia", None),
        ("chr2-200-209-AAG", "chr2", 200, 209, "AAG", 3, "", "4x:sampleC",
         "myopathy", None, "dystonia"),
    ])
    conn.commit()
    conn.close()


class PhenotypeKeywordAcrossAffectedOutliersTests(unittest.TestCase):
    """The fallback for a database with no swim_plot table: the summary phenotype columns.

    With swim_plot present the filter is answered from it (see the class below), because that is
    the table the swim plot itself matches. This fixture has no swim_plot at all, so the
    First/Second/ThirdAffectedPhenotype columns are the only phenotypes the database has and the
    filter searches all three of them rather than the first one only.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        _build_db_with_all_affected_phenotypes(cls.db_path)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def _locus_ids(self, keyword):
        response = self.client.get(f"/api/v1/loci?outlier_type=all&phenotype_keyword={keyword}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {r["LocusId"] for r in response.get_json()["results"]}

    def test_first_second_and_third_phenotypes_all_match(self):
        self.assertEqual(self._locus_ids("seizures"), {"chr1-100-110-AT"})
        # ataxia is the SECOND affected outlier's phenotype, and dystonia the third.
        self.assertEqual(self._locus_ids("ataxia"), {"chr1-100-110-AT"})
        self.assertEqual(self._locus_ids("dystonia"), {"chr2-200-209-AAG"})
        self.assertEqual(self._locus_ids("myopathy"), {"chr2-200-209-AAG"})

    def test_unmatched_keyword_returns_nothing(self):
        self.assertEqual(self._locus_ids("no_such_phenotype"), set())

    def test_keywords_are_or_combined_across_columns(self):
        self.assertEqual(self._locus_ids("ataxia,dystonia"),
                         {"chr1-100-110-AT", "chr2-200-209-AAG"})


def _build_db_with_swim_plot_phenotypes(path):
    """Loci plus the swim_plot rows the phenotypes really come from.

    chr1 has five outliers: three affected ones the loci table summarizes in its
    First/Second/ThirdAffectedPhenotype columns, a fourth affected one past the end of that
    summary, and an unaffected one. chr2 has a single affected outlier.
    """
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER,
        Motif TEXT, MotifSize INTEGER, KnownDiseaseLocus TEXT,
        OutlierSampleIds_AllAlleles TEXT,
        FirstAffectedPhenotype_AllAlleles TEXT,
        SecondAffectedPhenotype_AllAlleles TEXT,
        ThirdAffectedPhenotype_AllAlleles TEXT)""")
    conn.executemany("INSERT INTO loci VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", [
        ("chr1-100-110-AT", "chr1", 100, 110, "AT", 2, "",
         "9x:sampleA,8x:sampleB,7x:sampleC,6x:sampleD,5x:sampleE",
         "seizures", "ataxia", "dystonia"),
        ("chr2-200-209-AAG", "chr2", 200, 209, "AAG", 3, "", "4x:sampleF",
         "myopathy", None, None),
    ])

    def _swim_row(**values):
        """Return a full-width swim_plot row dict, defaulting every unnamed column to NULL."""
        return {column: values.get(column) for column in swim_plot.SWIM_PLOT_COLUMNS}

    result_database.write_swim_plot(
        conn,
        [
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                      allele_size=9, sample_id="sampleA", affected_status="Affected",
                      phenotype_description="seizures", LocusId="chr1-100-110-AT",
                      Motif="AT", CanonicalMotif="AT", MotifSize=2),
            _swim_row(outlier_type="AllAlleles", outlier_rank=2, motif_category="2bp",
                      allele_size=8, sample_id="sampleB", affected_status="Affected",
                      phenotype_description="ataxia", LocusId="chr1-100-110-AT",
                      Motif="AT", CanonicalMotif="AT", MotifSize=2),
            _swim_row(outlier_type="AllAlleles", outlier_rank=3, motif_category="2bp",
                      allele_size=7, sample_id="sampleC", affected_status="Affected",
                      phenotype_description="dystonia", LocusId="chr1-100-110-AT",
                      Motif="AT", CanonicalMotif="AT", MotifSize=2),
            # The fourth affected outlier: past the three ranks the loci table summarizes.
            _swim_row(outlier_type="AllAlleles", outlier_rank=4, motif_category="2bp",
                      allele_size=6, sample_id="sampleD", affected_status="Affected",
                      phenotype_description="neuropathy", LocusId="chr1-100-110-AT",
                      Motif="AT", CanonicalMotif="AT", MotifSize=2),
            # An unaffected outlier, which the summary columns never carry.
            _swim_row(outlier_type="AllAlleles", outlier_rank=5, motif_category="2bp",
                      allele_size=5, sample_id="sampleE", affected_status="Unaffected",
                      phenotype_description="headache", LocusId="chr1-100-110-AT",
                      Motif="AT", CanonicalMotif="AT", MotifSize=2),
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="3bp",
                      allele_size=4, sample_id="sampleF", affected_status="Affected",
                      phenotype_description="myopathy", LocusId="chr2-200-209-AAG",
                      Motif="AAG", CanonicalMotif="AAG", MotifSize=3),
        ],
        columns=swim_plot.SWIM_PLOT_COLUMNS,
    )
    conn.commit()
    conn.close()


class PhenotypeKeywordAgreesWithTheSwimPlotTests(unittest.TestCase):
    """A keyword selects a locus when ANY of its outliers carries that phenotype, in both views.

    The results table used to match the loci table's First/Second/ThirdAffectedPhenotype columns,
    which summarize at most the top three AFFECTED outliers, while the swim plot matched
    phenotype_description on every outlier row. A locus whose only matching phenotype belonged to
    its fourth affected outlier, or to an unaffected one, was therefore missing from the table
    while still being drawn in the swim plot.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        _build_db_with_swim_plot_phenotypes(cls.db_path)
        cls._saved_config = dict(results_server.app.config)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        results_server.app.config.clear()
        results_server.app.config.update(cls._saved_config)
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def _locus_ids(self, keyword):
        response = self.client.get(f"/api/v1/loci?outlier_type=all&phenotype_keyword={keyword}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {r["LocusId"] for r in response.get_json()["results"]}

    def _swim_locus_ids(self, keyword):
        response = self.client.get(
            f"/api/v1/swim_plot_data?outlier_type=all&phenotype_keyword={keyword}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {entry["LocusId"] for entry in response.get_json()["data"]}

    def test_the_top_three_affected_phenotypes_still_match(self):
        self.assertEqual(self._locus_ids("seizures"), {"chr1-100-110-AT"})
        self.assertEqual(self._locus_ids("ataxia"), {"chr1-100-110-AT"})
        self.assertEqual(self._locus_ids("dystonia"), {"chr1-100-110-AT"})
        self.assertEqual(self._locus_ids("myopathy"), {"chr2-200-209-AAG"})

    def test_a_fourth_affected_outlier_matches(self):
        # neuropathy belongs to the fourth affected outlier, which no summary column carries.
        self.assertEqual(self._locus_ids("neuropathy"), {"chr1-100-110-AT"})

    def test_an_unaffected_outlier_matches(self):
        self.assertEqual(self._locus_ids("headache"), {"chr1-100-110-AT"})

    def test_the_two_views_select_the_same_loci(self):
        for keyword in ("seizures", "ataxia", "dystonia", "neuropathy", "headache", "myopathy",
                        "no_such_phenotype"):
            self.assertEqual(self._locus_ids(keyword), self._swim_locus_ids(keyword), keyword)

    def test_keywords_are_or_combined(self):
        self.assertEqual(self._locus_ids("neuropathy,myopathy"),
                         {"chr1-100-110-AT", "chr2-200-209-AAG"})
        self.assertEqual(self._locus_ids("neuropathy,myopathy"),
                         self._swim_locus_ids("neuropathy,myopathy"))

    def test_matching_is_case_insensitive_and_partial(self):
        self.assertEqual(self._locus_ids("NEURO"), {"chr1-100-110-AT"})

    def test_an_unmatched_keyword_returns_nothing(self):
        self.assertEqual(self._locus_ids("no_such_phenotype"), set())


class ComputePlotThresholdsTests(unittest.TestCase):
    """The ExpansionHunter export's PlotReadVisualization thresholds."""

    AFFECTED = {"sampleA": "affected", "sampleB": "unaffected", "sampleC": "affected"}

    def test_largest_affected_above_the_largest_unaffected(self):
        self.assertEqual(
            results_server.compute_plot_thresholds({}, "6x:sampleA,5x:sampleB", "", self.AFFECTED),
            {"LongAllele": 6})

    def test_population_percentiles_raise_the_bar(self):
        # The comparison value is the largest of the unaffected outliers and the population
        # percentiles, so an affected allele below a percentile earns no threshold.
        self.assertIsNone(results_server.compute_plot_thresholds(
            {"HPRC256_99thPercentile": 10}, "6x:sampleA,5x:sampleB", "", self.AFFECTED))
        self.assertIsNone(results_server.compute_plot_thresholds(
            {"AoU1027_99thPercentile": 10}, "6x:sampleA,5x:sampleB", "", self.AFFECTED))
        # A percentile below the allele size leaves the threshold in place.
        self.assertEqual(
            results_server.compute_plot_thresholds(
                {"HPRC256_99thPercentile": 4, "AoU1027_99thPercentile": None},
                "6x:sampleA,5x:sampleB", "", self.AFFECTED),
            {"LongAllele": 6})

    def test_the_largest_outlier_has_to_be_affected(self):
        # sampleB is the largest here, and it is unaffected, so nothing qualifies even though a
        # smaller affected outlier follows it.
        self.assertIsNone(results_server.compute_plot_thresholds(
            {}, "7x:sampleB,6x:sampleA", "", self.AFFECTED))

    def test_short_allele_threshold(self):
        self.assertEqual(
            results_server.compute_plot_thresholds({}, "", "5x:sampleC", self.AFFECTED),
            {"ShortAllele": 5})
        # Both lists can qualify at once.
        self.assertEqual(
            results_server.compute_plot_thresholds(
                {}, "6x:sampleA,5x:sampleB", "5x:sampleC", self.AFFECTED),
            {"LongAllele": 6, "ShortAllele": 5})

    def test_no_outliers_and_unknown_samples(self):
        self.assertIsNone(results_server.compute_plot_thresholds({}, "", "", self.AFFECTED))
        # A sample with no metadata row is not known to be affected.
        self.assertIsNone(results_server.compute_plot_thresholds({}, "6x:sampleZ", "", {}))
        # Entries that do not parse are skipped rather than raising.
        self.assertIsNone(results_server.compute_plot_thresholds({}, "garbage", "", self.AFFECTED))


def _build_db_with_unknown_category(path):
    """Minimal loci + swim_plot DB with a normal row and an 'Unknown' motif row."""
    conn = duckdb_compat.connect(path)
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


class SchemaAndUnknownCategoryTests(unittest.TestCase):
    """WS3 + WS4: schema sample_ids fall back to swim_plot; Unknown rows served."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
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



def _build_db_with_sample_count(path, num_samples):
    """Minimal loci database that also carries the build's sample-count metadata row."""
    conn = duckdb_compat.connect(path)
    conn.execute("CREATE TABLE loci (LocusId TEXT, Motif TEXT, KnownDiseaseLocus TEXT)")
    conn.execute("INSERT INTO loci VALUES ('chr1-100-110-AT', 'AT', '')")
    result_database.write_sample_count(conn, num_samples)
    conn.commit()
    conn.close()


class SchemaSampleCountTests(unittest.TestCase):
    """/api/v1/schema reports the sample count build_database recorded in the database."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        _build_db_with_sample_count(cls.db_path, 2657)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def test_schema_reports_total_samples(self):
        data = self.client.get("/api/v1/schema").get_json()
        self.assertEqual(data["total_samples"], 2657)
        self.assertEqual(data["total_loci"], 1)

    def test_configure_app_reads_the_count_from_the_database(self):
        # The schema value comes from configure_app's one read at startup, not from a per-request
        # query, so pin the config entry the endpoint reads.
        self.assertEqual(results_server.app.config["TOTAL_SAMPLES"], 2657)


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


class _QueryRecordingConnection:
    """Wraps a DuckDB connection and records the SQL text of every execute() call."""

    def __init__(self, conn, recorded_queries):
        self._conn = conn
        self._recorded_queries = recorded_queries

    def execute(self, sql, *args, **kwargs):
        self._recorded_queries.append(sql)
        return self._conn.execute(sql, *args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._conn, name)


def _record_queries(test_case):
    """Route the server's get_db() through a recorder for the duration of one test.

    Args:
        test_case: The running TestCase, used to restore get_db() on cleanup.

    Returns:
        The list the executed SQL statements are appended to.
    """
    recorded_queries = []
    original_get_db = results_server.get_db

    def _recording_get_db():
        return _QueryRecordingConnection(original_get_db(), recorded_queries)

    results_server.get_db = _recording_get_db
    test_case.addCleanup(setattr, results_server, "get_db", original_get_db)
    return recorded_queries


def _build_db_with_phenotype_scores(path):
    """The minimal fixture plus the two phenotype-score tables, for both loci."""
    _build_minimal_db(path)
    conn = duckdb_compat.connect(path)
    try:
        result_database.write_phenotype_tables(
            conn,
            [
                {"locus_id": "chr1-100-110-AT", "sample_id": "sampleA",
                 "outlier_type": "AllAlleles", "allele_size": 6, "gene_symbol": "GENE1",
                 "gene_phenotype_similarity": 0.75, "gene_phenotype_overlap_count": 2,
                 "n_matching_diseases": 1, "best_matching_disease": "disease1",
                 "best_disease_inheritance": "AD", "pairwise_similarity_to_next": 0.25,
                 "pairwise_shared_count_raw": 3, "pairwise_shared_count_ic": 1.5,
                 "next_sample_id": "sampleB"},
                {"locus_id": "chr2-200-209-AAG", "sample_id": "sampleC",
                 "outlier_type": "AllAlleles", "allele_size": 4, "gene_symbol": "GENE2",
                 "gene_phenotype_similarity": 0.5, "gene_phenotype_overlap_count": 1,
                 "n_matching_diseases": 0, "best_matching_disease": None,
                 "best_disease_inheritance": None, "pairwise_similarity_to_next": None,
                 "pairwise_shared_count_raw": None, "pairwise_shared_count_ic": None,
                 "next_sample_id": None},
            ],
            [
                {"locus_id": "chr1-100-110-AT", "outlier_type": "AllAlleles",
                 "num_qualifying_samples": 1, "sum_pairwise_similarity": 0.25,
                 "sum_pairwise_shared_raw": 3, "sum_pairwise_shared_ic": 1.5,
                 "max_gene_phenotype_similarity": 0.75, "qualifying_sample_ids": "sampleA"},
            ],
        )
    finally:
        conn.close()


class ExportPhenotypeScoreQueryCountTests(unittest.TestCase):
    """Both enriching exports bulk-load the phenotype scores once, not once per locus."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        _build_db_with_phenotype_scores(cls.db_path)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def _score_queries(self, recorded_queries):
        """The recorded statements that read the per-outlier phenotype-score table."""
        return [sql for sql in recorded_queries if "per_outlier_phenotype_scores" in sql]

    def test_expansion_hunter_export_reads_the_score_table_once(self):
        # The enrichment branch used to call build_sample_details without phenotype_scores, so
        # every locus ran its own "WHERE locus_id = ?" query (plus a table listing).
        recorded_queries = _record_queries(self)
        response = self.client.get("/api/v1/export?outlier_type=all&format=expansion_hunter")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        catalog = {entry["LocusId"]: entry for entry in json.loads(response.get_data())}
        self.assertEqual(len(catalog), 2)
        score_queries = self._score_queries(recorded_queries)
        self.assertEqual(len(score_queries), 1, score_queries)
        self.assertNotIn("WHERE locus_id = ?", score_queries[0])

    def test_expansion_hunter_export_still_carries_the_scores(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=expansion_hunter")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        catalog = {entry["LocusId"]: entry for entry in json.loads(response.get_data())}
        samples = {s["Sample ID"]: s
                   for s in catalog["chr1-100-110-AT"]["Samples"]["AllAlleleOutliers"]}
        self.assertEqual(samples["sampleA"]["Gene-Pheno Sim"], 0.75)
        self.assertEqual(samples["sampleA"]["Pairwise Sim"], 0.25)
        self.assertEqual(
            {s["Sample ID"]: s["Gene-Pheno Sim"]
             for s in catalog["chr2-200-209-AAG"]["Samples"]["AllAlleleOutliers"]},
            {"sampleC": 0.5})

    def test_json_export_reads_the_score_table_once(self):
        # The path the expansion_hunter export was brought back into line with.
        recorded_queries = _record_queries(self)
        response = self.client.get("/api/v1/export?outlier_type=all&format=json")
        # The body is gzipped, so it is not decodable as text for an assertion message.
        self.assertEqual(response.status_code, 200)
        data = json.loads(gzip.decompress(response.get_data()))
        self.assertEqual(data["metadata"]["total_loci"], 2)
        self.assertEqual(len(self._score_queries(recorded_queries)), 1)


def _build_db_with_only_the_first_phenotype_column(path):
    """Loci carrying FirstAffectedPhenotype_AllAlleles only, and no swim_plot table."""
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER,
        Motif TEXT, MotifSize INTEGER, KnownDiseaseLocus TEXT,
        OutlierSampleIds_AllAlleles TEXT,
        FirstAffectedPhenotype_AllAlleles TEXT)""")
    conn.executemany("INSERT INTO loci VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", [
        ("chr1-100-110-AT", "chr1", 100, 110, "AT", 2, "", "6x:sampleA", "seizures"),
        ("chr2-200-209-AAG", "chr2", 200, 209, "AAG", 3, "", "4x:sampleC", "ataxia"),
    ])
    conn.commit()
    conn.close()


class PhenotypeKeywordFallbackColumnSubsetTests(unittest.TestCase):
    """The fallback narrowed to the phenotype columns a database actually declares.

    Without a swim_plot table the filter reads the summary columns, and it must name only the
    ones that exist: naming SecondAffectedPhenotype_AllAlleles here would be a SQL error on
    every request carrying a phenotype keyword.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        _build_db_with_only_the_first_phenotype_column(cls.db_path)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def test_the_swim_plot_branch_is_not_taken(self):
        self.assertFalse(results_server.app.config["HAS_SWIM_PLOT_PHENOTYPES"])

    def test_only_the_declared_column_is_named(self):
        select_query, _, _, _, sql_params, _ = results_server.build_api_query({
            "outlier_type": "all", "page": 1, "page_size": 50,
            "phenotype_keyword": ["seizures"],
        })
        self.assertIn("FirstAffectedPhenotype_AllAlleles ILIKE ?", select_query)
        self.assertNotIn("SecondAffectedPhenotype_AllAlleles", select_query)
        self.assertNotIn("ThirdAffectedPhenotype_AllAlleles", select_query)
        self.assertNotIn("swim_plot", select_query)
        self.assertEqual(sql_params[-1:], ["%seizures%"])

    def test_the_keyword_still_matches_through_the_endpoint(self):
        response = self.client.get("/api/v1/loci?outlier_type=all&phenotype_keyword=seizures")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual({r["LocusId"] for r in response.get_json()["results"]},
                         {"chr1-100-110-AT"})

    def test_an_unmatched_keyword_returns_nothing(self):
        response = self.client.get("/api/v1/loci?outlier_type=all&phenotype_keyword=myopathy")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual(response.get_json()["results"], [])


def _build_db_without_phenotype_columns(path):
    """Loci with no AffectedPhenotype column at all, and no swim_plot table."""
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER,
        Motif TEXT, MotifSize INTEGER, KnownDiseaseLocus TEXT,
        OutlierSampleIds_AllAlleles TEXT)""")
    conn.executemany("INSERT INTO loci VALUES (?, ?, ?, ?, ?, ?, ?, ?)", [
        ("chr1-100-110-AT", "chr1", 100, 110, "AT", 2, "", "6x:sampleA"),
    ])
    conn.commit()
    conn.close()


class PhenotypeKeywordWithoutAnyPhenotypeColumnTests(unittest.TestCase):
    """A phenotype keyword matches nothing when the database records no phenotypes at all."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        _build_db_without_phenotype_columns(cls.db_path)
        results_server.configure_app(cls.db_path, annotations_db=cls.annotations_db)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def test_the_filter_matches_no_rows_rather_than_being_ignored(self):
        select_query, _, _, _, _, _ = results_server.build_api_query({
            "outlier_type": "all", "page": 1, "page_size": 50,
            "phenotype_keyword": ["seizures"],
        })
        self.assertIn("1=0", select_query)
        self.assertNotIn("AffectedPhenotype", select_query)
        response = self.client.get("/api/v1/loci?outlier_type=all&phenotype_keyword=seizures")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual(response.get_json()["total"], 0)
        # Without the keyword the same fixture does return its one locus, so the empty result
        # above is the filter working and not an empty database.
        unfiltered = self.client.get("/api/v1/loci?outlier_type=all")
        self.assertEqual(unfiltered.get_json()["total"], 1)


def _load_results_server_without_locus_annotations():
    """Import a second copy of results_server.py with locus_annotations unimportable.

    This is the degraded mode the module's ImportError fallbacks exist for: locus_annotations
    imports intervaltree at module scope, so a missing intervaltree makes the import raise and
    the server falls back to its built-in stubs. Loading the file under a different module name
    leaves the results_server module the rest of these tests share untouched.

    Returns:
        The freshly executed module object, holding the fallback definitions.
    """
    spec = importlib.util.spec_from_file_location(
        "results_server_degraded_mode", results_server.__file__)
    module = importlib.util.module_from_spec(spec)
    missing = object()
    saved = sys.modules.get("locus_annotations", missing)
    sys.modules["locus_annotations"] = None  # makes "from locus_annotations import ..." raise
    try:
        spec.loader.exec_module(module)
    finally:
        if saved is missing:
            del sys.modules["locus_annotations"]
        else:
            sys.modules["locus_annotations"] = saved
    return module


class DegradedModeFallbackTests(unittest.TestCase):
    """The ImportError stubs must classify affected status the way locus_annotations does."""

    @classmethod
    def setUpClass(cls):
        cls.degraded = _load_results_server_without_locus_annotations()

    def test_the_fallback_really_is_a_stub(self):
        # If this fails the import did not raise and the assertions below would be testing the
        # real locus_annotations functions rather than the fallbacks.
        self.assertEqual(self.degraded.load_known_disease_loci("ignored"), ({}, {}, {}))

    def test_every_blank_spelling_normalizes_to_none(self):
        # Mapping only "nan" left an "NA" cell classified as the status "na", so the locus-detail
        # page showed "Na" instead of "Unknown" whenever the fallback was in use.
        for value in ("", "nan", "none", "na", "n/a", "null", "NA", " N/A ", "Null"):
            self.assertIsNone(self.degraded.normalize_affected_status_for_logic(value), value)

    def test_the_fallback_agrees_with_locus_annotations(self):
        for value in ("", "nan", "none", "na", "n/a", "null", "NA", " N/A ", "Null", None,
                      float("nan"), "Affected", "Possibly Affected", " UNAFFECTED "):
            self.assertEqual(
                self.degraded.normalize_affected_status_for_logic(value),
                locus_annotations.normalize_affected_status_for_logic(value),
                repr(value))
