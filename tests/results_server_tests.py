"""Unit tests for the TRails Flask server (results_server.py).

These tests build a tiny temporary DuckDB database with a minimal ``loci`` table
(plus a ``swim_plot`` table) and exercise the server through Flask's test client,
so no real TCP port is bound. They check that:

  * GET /api/v1/loci returns rows from the database, and
  * the pages and endpoints (schema, swim, export, qc) behave gracefully.
"""

import gzip
import importlib.util
import json
import os
import sys
import tempfile
import threading
import unittest
from unittest import mock

import intervaltree

import analysis_columns
import duckdb_compat
import input_tables
import locus_annotations
import result_database
import results_server
import swim_plot
from tests.server_test_fixtures import ServerTestCase, make_swim_plot_row


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
    result_database.write_swim_plot(
        conn,
        [
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                               allele_size=6, sample_id="sampleA", family_id="famA",
                               affected_status="Affected", analysis_status="Unsolved", sex="male",
                               phenotype_description="seizures", FirstUnaffectedAlleleSize=4,
                               is_above_first_unaffected=1, LocusId="chr1-100-110-AT", Motif="AT",
                               CanonicalMotif="AT", MotifSize=2, gene_region="CDS",
                               GeneTableGeneSymbol="GENE1"),
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=2, motif_category="2bp",
                               allele_size=5, sample_id="sampleB", family_id="famB",
                               affected_status="Unaffected", analysis_status="Unaffected", sex="female",
                               FirstUnaffectedAlleleSize=4, is_above_first_unaffected=0,
                               LocusId="chr1-100-110-AT", Motif="AT", CanonicalMotif="AT", MotifSize=2,
                               gene_region="CDS", GeneTableGeneSymbol="GENE1"),
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="3bp",
                               allele_size=4, sample_id="sampleC", family_id="famC",
                               affected_status="Affected", analysis_status="Unsolved", sex="male",
                               phenotype_description="ataxia", is_above_first_unaffected=1,
                               LocusId="chr2-200-209-AAG", Motif="AAG", CanonicalMotif="AAG", MotifSize=3,
                               gene_region="intron", GeneTableGeneSymbol="GENE2"),
        ],
        columns=swim_plot.SWIM_PLOT_COLUMNS,
    )
    conn.commit()
    conn.close()


class MinimalDatabaseTestCase(ServerTestCase):
    """Points the server at the _build_minimal_db fixture."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_minimal_db(cls.db_path)
        cls.configure_server(cls.db_path)


class ResultsServerTests(MinimalDatabaseTestCase):
    def test_loci_returns_rows(self):
        response = self.client.get("/api/v1/loci?outlier_type=all")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["total"], 2)
        self.assertEqual(len(data["results"]), 2)
        locus_ids = {r["LocusId"] for r in data["results"]}
        self.assertIn("chr1-100-110-AT", locus_ids)
        self.assertIn("chr2-200-209-AAG", locus_ids)

    def test_unknown_query_parameters_are_ignored(self):
        response = self.client.get("/api/v1/loci?outlier_type=all&unrecognized=foo")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["total"], 2)
        # An unrecognized parameter is ignored and not reported in filters_applied.
        self.assertNotIn("unrecognized", data["filters_applied"])

    def test_outlier_type_required(self):
        response = self.client.get("/api/v1/loci")
        self.assertEqual(response.status_code, 400)

    def test_locus_detail(self):
        response = self.client.get("/api/v1/loci/chr1-100-110-AT")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["locus"]["LocusId"], "chr1-100-110-AT")
        self.assertIn("AllAlleles", data["outlier_samples"])

    def test_schema_reports_the_database_total_loci(self):
        response = self.client.get("/api/v1/schema")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_json()["total_loci"], 2)

    def test_schema_total_samples_is_null_without_the_metadata_table(self):
        # This fixture has no metadata table, so the key is reported as null rather than
        # omitted: the index page checks for null before appending ", N samples" to the count.
        data = self.client.get("/api/v1/schema").get_json()
        self.assertIn("total_samples", data)
        self.assertIsNone(data["total_samples"])

    def test_schema_lists_the_population_metric_values_the_loci_endpoint_accepts(self):
        values = self.client.get("/api/v1/schema").get_json()["filters"]["population_metric"]["values"]
        self.assertEqual(values, ["99th", "max"])
        for value in values:
            self.assertEqual(self.client.get(f"/api/v1/loci?outlier_type=all&population_metric={value}").status_code,
                             200, value)
        self.assertEqual(self.client.get("/api/v1/loci?outlier_type=all&population_metric=median").status_code, 400)

    def test_the_page_escapes_the_gate_cohort_names(self):
        with mock.patch.dict(results_server.app.config, {"POPULATION_GATE_COHORTS": ('<img src=x onerror="f()">',)}):
            page = self.client.get("/").get_data(as_text=True)
        self.assertNotIn('<img src=x onerror="f()">', page)
        # The #Fam title attribute shows the name as text. The Group Results By help popup inserts
        # its data-html attribute's decoded value as HTML, so the name is escaped once more there.
        self.assertIn("built with (&lt;img src=x onerror=&#34;f()&#34;&gt;; #Aff counted by family)", page)
        self.assertIn("built with (&amp;lt;img src=x onerror=&amp;#34;f()&amp;#34;&amp;gt;). A row expands", page)

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
        # The swim plot applies these flags through the same locus filters as the results table,
        # so without the column it keeps returning its data too.
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
        records = [line.split("\t") for line in response.get_data(as_text=True).splitlines()]
        self.assertEqual(records, [["chr1", "100", "110", "chr1-100-110-AT", "0", "."],
                                   ["chr2", "200", "209", "chr2-200-209-AAG", "0", "."]])
        # Every record is valid BED6: an integer interval with start before end, a name, a score
        # from 0 to 1000 and a strand of +, - or .
        for chrom, start, end, name, score, strand in records:
            self.assertTrue(chrom and name)
            self.assertLess(int(start), int(end))
            self.assertTrue(0 <= int(score) <= 1000)
            self.assertIn(strand, ("+", "-", "."))

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
        # Every qualifying outlier sample appears in the nested Samples list, not just the first.
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
        # No sample table was supplied, so every outlier is of unknown status and, as #Aff counts
        # it, a candidate: the largest outlier of each list earns a read-visualization threshold.
        self.assertEqual(entry["PlotReadVisualization"], [{"If": "LongAllele", "Is": ">=", "Threshold": 6},
                                                          {"If": "ShortAllele", "Is": ">=", "Threshold": 5}])

    def test_export_expansion_hunter_read_visualization_thresholds(self):
        # chr1's largest outlier (6x:sampleA) is affected and sits above the largest unaffected
        # outlier (5x:sampleB), so it earns a LongAllele threshold. The ShortAlleles outlier list
        # holds only the unaffected sampleB, so there is no ShortAllele threshold. chr2's only
        # outlier, sampleC, is affected but solved, so it is part of the bar and earns nothing.
        with mock.patch.dict(results_server.app.config["LOOKUPS"], {
                "affected": {"sampleA": "affected", "sampleB": "unaffected", "sampleC": "affected"},
                "analysis": {"sampleC": "solved"}}):
            response = self.client.get("/api/v1/export?outlier_type=all&format=expansion_hunter")
            # The body streams, so it is read while the patched lookups are in place.
            self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
            catalog = {entry["LocusId"]: entry for entry in json.loads(response.get_data())}
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
        # The tag list the pages' menus are filled from now includes it
        tag_list = self.client.get("/api/v1/annotations/tags").get_json()
        self.assertEqual(["tags"], sorted(tag_list))
        self.assertIn("interesting", tag_list["tags"])

    def test_concurrent_annotation_writes_all_succeed(self):
        # Flask serves with threaded=True and DuckDB raises instead of waiting when two
        # connections write the same row at once, so the annotation writes are serialized and
        # overlapping writes all succeed.
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

    def test_the_in_memory_annotations_change_under_the_write_lock(self):
        # Otherwise two overlapping writes could commit in one order and update the in-memory copy
        # in the other, leaving it out of step with the database.
        lock_held_during_update = []

        def update_cached_annotations(annotations):
            lock_held_during_update.append(results_server.ANNOTATION_WRITE_LOCK.locked())
            return {"updated": True}

        self.assertEqual(
            results_server.execute_annotation_write(
                "DELETE FROM notes WHERE locus_id = ?", ("no-such-locus",), update_cached_annotations),
            (None, {"updated": True}))
        self.assertEqual(lock_held_during_update, [True])

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

    def test_initial_require_above_population_reaches_the_index_page(self):
        self.addCleanup(setattr, sys, "argv", sys.argv)
        sys.argv = ["results_server.py", "--db", "unused.duckdb", "--initial-require-above-population", "all"]
        self.assertEqual(results_server.parse_args().initial_require_above_population, "all")
        sys.argv = ["results_server.py", "--db", "unused.duckdb"]
        self.assertEqual(results_server.parse_args().initial_require_above_population, "")

        render_arguments = self.record_template_render_arguments()
        self.addCleanup(results_server.app.config.pop, "INITIAL_REQUIRE_ABOVE_POPULATION", None)
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(render_arguments["initial_require_above_population"], "")
        results_server.app.config["INITIAL_REQUIRE_ABOVE_POPULATION"] = "all"
        self.client.get("/")
        self.assertEqual(render_arguments["initial_require_above_population"], "all")

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

    def test_sample_id_like_does_not_match_allele_sizes(self):
        # OutlierSampleIds_AllAlleles is "6x:sampleA,5x:sampleB" at chr1 and "4x:sampleC" at chr2,
        # so a numeric keyword must not match the allele size in front of every sample id.
        self.assertEqual(self.get_locus_ids("sample_id_like=6"), set())
        self.assertEqual(self.get_locus_ids("sample_id_like=4"), set())
        # The "x:" delimiter is not part of a sample id either.
        self.assertEqual(self.get_locus_ids("sample_id_like=x:"), set())

    def test_sample_id_like_still_matches_sample_ids(self):
        self.assertEqual(self.get_locus_ids("sample_id_like=sampleA"), {"chr1-100-110-AT"})
        self.assertEqual(self.get_locus_ids("sample_id_like=ampleC"), {"chr2-200-209-AAG"})
        self.assertEqual(self.get_locus_ids("sample_id_like=sample"),
                         {"chr1-100-110-AT", "chr2-200-209-AAG"})

    def test_sample_id_like_treats_underscore_literally(self):
        # "_" is ILIKE's single-character wildcard, so without escaping this would match "sampleA".
        self.assertEqual(self.get_locus_ids("sample_id_like=sample_A"), set())

    def test_phenotype_keyword_is_answered_from_swim_plot(self):
        # This fixture writes a swim_plot whose phenotype_description column is populated, so the
        # filter is answered from that table and names no AffectedPhenotype column at all. The
        # fallback for a database without swim_plot phenotypes is covered by
        # PhenotypeKeywordFallbackColumnSubsetTests below.
        self.assertTrue(results_server.app.config["HAS_SWIM_PLOT_PHENOTYPES"])
        self.assertEqual(self.get_locus_ids("phenotype_keyword=seizures"), {"chr1-100-110-AT"})

    # -- Merged search box: locus ids -------------------------------------------

    def test_pasted_locus_ids_become_one_set_lookup(self):
        # A predicate per pasted id costs terms x rows, which for a few hundred ids takes
        # seconds, so the ids go in a single IN list, lower-cased on both sides so the match
        # stays case-insensitive. search_filter_tests.py checks the matching itself
        # (case-insensitive, a literal underscore); this test pins the query's shape.
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
    result_database.write_swim_plot(
        conn,
        [
            # Both alleles of one diploid sample at one locus, the way accumulate_locus records
            # them: two swim_plot rows for the same (sample, locus).
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                               allele_size=9, sample_id="sampleA", LocusId="chr1-100-110-AT",
                               Motif="AT", MotifSize=2),
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=2, motif_category="2bp",
                               allele_size=8, sample_id="sampleA", LocusId="chr1-100-110-AT",
                               Motif="AT", MotifSize=2),
            # A second locus in the same motif-size bin, so the expected count is not 1 by default.
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                               allele_size=7, sample_id="sampleA", LocusId="chr2-200-210-AT",
                               Motif="AT", MotifSize=2),
            # Another sample's entry at the first locus, to keep the grouping honest.
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=3, motif_category="2bp",
                               allele_size=6, sample_id="sampleB", LocusId="chr1-100-110-AT",
                               Motif="AT", MotifSize=2),
        ],
        columns=swim_plot.SWIM_PLOT_COLUMNS,
    )
    conn.commit()
    conn.close()


class ComputeOutlierWarningsTests(unittest.TestCase):
    """The reported bound is the smallest count that is flagged, so "≥ bound" says who is flagged."""

    def test_a_fractional_99th_percentile_is_reported_rounded_up(self):
        # 98 samples with 1 rank-1 locus, then 4 and 5: the 99th percentile is 4.01
        rank1_data = [{"sample_id": f"S{index}", "bin": "2bp", "count": 1} for index in range(98)]
        rank1_data += [{"sample_id": "FOUR", "bin": "2bp", "count": 4}, {"sample_id": "FIVE", "bin": "2bp", "count": 5}]
        warnings = results_server.compute_outlier_warnings(rank1_data)
        self.assertEqual(warnings, {"FIVE": {"outlier_bins": {"2bp": {"count": 5, "threshold": 5}}}})
        # Every flagged count is at least the reported bound, and every other count is below it
        bound = warnings["FIVE"]["outlier_bins"]["2bp"]["threshold"]
        for entry in rank1_data:
            self.assertEqual(entry["sample_id"] in warnings, entry["count"] >= bound, entry)


class SampleQcLocusCountTests(ServerTestCase):
    """The sample-QC counts are counts of loci, not of swim_plot allele rows."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_for_sample_qc(cls.db_path)

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
        "disease": "Test pathogenic-motif disease",
        # Every record in the shipped STRchive file stores inheritance as a list.
        "inheritance": ["AD"],
        "pathogenic_min": 60,
        "chrom": "chr9",
        # STRchive's start is 0-based and its stop exclusive, so this is the interval [100, 110).
        "start_hg38": 100,
        "stop_hg38": 110,
        "reference_motif_reference_orientation": ["AAAAG"],
        "pathogenic_motif_reference_orientation": ["AAGGG"],
    },
    {
        "id": "TEST_REFERENCE",
        "disease": "Test reference-motif disease",
        # A bare string, the shape an older cached STRchive file may still use.
        "inheritance": "AR",
        "pathogenic_min": 50,
        "chrom": "chr3",
        "start_hg38": 300,
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


class StrchivePathogenicMotifTests(ServerTestCase):
    """The locus-detail endpoint and the build agree on STRchive pathogenic-motif matches."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        cls.strchive_json = os.path.join(cls.tmpdir, "STRchive-loci.json")
        with open(cls.strchive_json, "w") as output_file:
            json.dump(STRCHIVE_TEST_LOCI, output_file)
        _build_db_for_strchive(cls.db_path)
        cls.configure_server(cls.db_path, strchive_loci_json=cls.strchive_json)

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


class PhenotypeKeywordAcrossAffectedOutliersTests(ServerTestCase):
    """The fallback for a database with no swim_plot table: the summary phenotype columns.

    With swim_plot present the filter is answered from it (see the class below), because that is
    the table the swim plot itself matches. This fixture has no swim_plot at all, so the
    First/Second/ThirdAffectedPhenotype columns are the only phenotypes the database has and the
    filter searches all three of them rather than the first one only.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_all_affected_phenotypes(cls.db_path)
        cls.configure_server(cls.db_path)

    def test_first_second_and_third_phenotypes_all_match(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=seizures"), {"chr1-100-110-AT"})
        # ataxia is the SECOND affected outlier's phenotype, and dystonia the third.
        self.assertEqual(self.get_locus_ids("phenotype_keyword=ataxia"), {"chr1-100-110-AT"})
        self.assertEqual(self.get_locus_ids("phenotype_keyword=dystonia"), {"chr2-200-209-AAG"})
        self.assertEqual(self.get_locus_ids("phenotype_keyword=myopathy"), {"chr2-200-209-AAG"})

    def test_unmatched_keyword_returns_nothing(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=no_such_phenotype"), set())

    def test_keywords_are_or_combined_across_columns(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=ataxia,dystonia"),
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
    result_database.write_swim_plot(
        conn,
        [
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                               allele_size=9, sample_id="sampleA", affected_status="Affected",
                               phenotype_description="seizures", LocusId="chr1-100-110-AT",
                               Motif="AT", CanonicalMotif="AT", MotifSize=2),
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=2, motif_category="2bp",
                               allele_size=8, sample_id="sampleB", affected_status="Affected",
                               phenotype_description="ataxia", LocusId="chr1-100-110-AT",
                               Motif="AT", CanonicalMotif="AT", MotifSize=2),
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=3, motif_category="2bp",
                               allele_size=7, sample_id="sampleC", affected_status="Affected",
                               phenotype_description="dystonia", LocusId="chr1-100-110-AT",
                               Motif="AT", CanonicalMotif="AT", MotifSize=2),
            # The fourth affected outlier: past the three ranks the loci table summarizes.
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=4, motif_category="2bp",
                               allele_size=6, sample_id="sampleD", affected_status="Affected",
                               phenotype_description="neuropathy", LocusId="chr1-100-110-AT",
                               Motif="AT", CanonicalMotif="AT", MotifSize=2),
            # An unaffected outlier, which the summary columns never carry.
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=5, motif_category="2bp",
                               allele_size=5, sample_id="sampleE", affected_status="Unaffected",
                               phenotype_description="headache", LocusId="chr1-100-110-AT",
                               Motif="AT", CanonicalMotif="AT", MotifSize=2),
            make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="3bp",
                               allele_size=4, sample_id="sampleF", affected_status="Affected",
                               phenotype_description="myopathy", LocusId="chr2-200-209-AAG",
                               Motif="AAG", CanonicalMotif="AAG", MotifSize=3),
        ],
        columns=swim_plot.SWIM_PLOT_COLUMNS,
    )
    conn.commit()
    conn.close()


class PhenotypeKeywordAgreesWithTheSwimPlotTests(ServerTestCase):
    """A keyword selects a locus when ANY of its outliers carries that phenotype, in both views.

    The results table matches every outlier's phenotype, as the swim plot does, rather than only
    the loci table's First/Second/ThirdAffectedPhenotype columns, which summarize at most the top
    three AFFECTED outliers. A locus whose only matching phenotype belongs to its fourth affected
    outlier, or to an unaffected one, is therefore in both views.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_swim_plot_phenotypes(cls.db_path)
        cls.configure_server(cls.db_path)

    def _swim_locus_ids(self, keyword):
        response = self.client.get(
            f"/api/v1/swim_plot_data?outlier_type=all&phenotype_keyword={keyword}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {entry["LocusId"] for entry in response.get_json()["data"]}

    def test_the_top_three_affected_phenotypes_still_match(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=seizures"), {"chr1-100-110-AT"})
        self.assertEqual(self.get_locus_ids("phenotype_keyword=ataxia"), {"chr1-100-110-AT"})
        self.assertEqual(self.get_locus_ids("phenotype_keyword=dystonia"), {"chr1-100-110-AT"})
        self.assertEqual(self.get_locus_ids("phenotype_keyword=myopathy"), {"chr2-200-209-AAG"})

    def test_a_fourth_affected_outlier_matches(self):
        # neuropathy belongs to the fourth affected outlier, which no summary column carries.
        self.assertEqual(self.get_locus_ids("phenotype_keyword=neuropathy"), {"chr1-100-110-AT"})

    def test_an_unaffected_outlier_matches(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=headache"), {"chr1-100-110-AT"})

    def test_the_two_views_select_the_same_loci(self):
        for keyword in ("seizures", "ataxia", "dystonia", "neuropathy", "headache", "myopathy",
                        "no_such_phenotype"):
            self.assertEqual(self.get_locus_ids(f"phenotype_keyword={keyword}"), self._swim_locus_ids(keyword),
                             keyword)

    def test_keywords_are_or_combined(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=neuropathy,myopathy"),
                         {"chr1-100-110-AT", "chr2-200-209-AAG"})
        self.assertEqual(self.get_locus_ids("phenotype_keyword=neuropathy,myopathy"),
                         self._swim_locus_ids("neuropathy,myopathy"))

    def test_matching_is_case_insensitive_and_partial(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=NEURO"), {"chr1-100-110-AT"})

    def test_an_unmatched_keyword_returns_nothing(self):
        self.assertEqual(self.get_locus_ids("phenotype_keyword=no_such_phenotype"), set())


class SplitOutlierEntriesTests(unittest.TestCase):
    """The one OutlierSampleIds parser behind the locus-detail samples and the read-visualization thresholds."""

    OUTLIER_STRING = "6x:sA,5:sB,garbage,x:sC,7x:sD:extra,4x:sE"

    def test_entries_that_do_not_parse_are_skipped_and_the_order_is_kept(self):
        self.assertEqual(results_server.split_outlier_entries(self.OUTLIER_STRING),
                         [(6, "sA"), (5, "sB"), (4, "sE")])
        self.assertEqual(results_server.split_outlier_entries(None), [])
        self.assertEqual(results_server.split_outlier_entries(""), [])

    def test_parse_outlier_samples_describes_each_sample(self):
        parsed = results_server.parse_outlier_samples(
            "6x:sA,5x:sB", {"sA": {"family_id": "F1", "sex": "male"}})
        self.assertEqual([(s["allele_size"], s["sample_id"], s["family_id"], s["sex"]) for s in parsed],
                         [(6, "sA", "F1", "male"), (5, "sB", None, None)])


class ComputePlotThresholdsTests(unittest.TestCase):
    """The ExpansionHunter export's PlotReadVisualization thresholds."""

    AFFECTED = {"sampleA": "affected", "sampleB": "unaffected", "sampleC": "affected"}
    LONG_READ_GATE = ("HPRC256", "AoU1027")

    def _thresholds(self, row_dict, long_outliers, short_outliers, affected_lookup, analysis_lookup=None,
                    population_gate_cohorts=LONG_READ_GATE):
        return results_server.compute_plot_thresholds(
            row_dict, long_outliers, short_outliers, affected_lookup, analysis_lookup or {},
            population_gate_cohorts)

    def test_largest_affected_above_the_largest_unaffected(self):
        self.assertEqual(self._thresholds({}, "6x:sampleA,5x:sampleB", "", self.AFFECTED), {"LongAllele": 6})

    def test_population_percentiles_raise_the_bar(self):
        # The comparison value is the largest of the unaffected outliers and the population
        # percentiles, so an affected allele below a percentile earns no threshold.
        self.assertIsNone(self._thresholds({"HPRC256_99thPercentile": 10}, "6x:sampleA,5x:sampleB", "", self.AFFECTED))
        self.assertIsNone(self._thresholds({"AoU1027_99thPercentile": 10}, "6x:sampleA,5x:sampleB", "", self.AFFECTED))
        # A percentile below the allele size leaves the threshold in place.
        self.assertEqual(
            self._thresholds({"HPRC256_99thPercentile": 4, "AoU1027_99thPercentile": None},
                             "6x:sampleA,5x:sampleB", "", self.AFFECTED),
            {"LongAllele": 6})

    def test_only_the_gate_cohorts_raise_the_bar(self):
        # TenK10K counts only when the build gated its counts with it.
        row_dict = {"HPRC256_99thPercentile": 4, "TenK10K_99thPercentile": 10}
        self.assertEqual(self._thresholds(row_dict, "6x:sampleA,5x:sampleB", "", self.AFFECTED), {"LongAllele": 6})
        self.assertIsNone(self._thresholds(row_dict, "6x:sampleA,5x:sampleB", "", self.AFFECTED,
                                           population_gate_cohorts=("HPRC256", "AoU1027", "TenK10K")))

    def test_the_largest_outlier_has_to_be_affected(self):
        # sampleB is the largest here, and it is unaffected, so nothing qualifies even though a
        # smaller affected outlier follows it.
        self.assertIsNone(self._thresholds({}, "7x:sampleB,6x:sampleA", "", self.AFFECTED))

    def test_short_allele_threshold(self):
        self.assertEqual(self._thresholds({}, "", "5x:sampleC", self.AFFECTED), {"ShortAllele": 5})
        # Both lists can qualify at once.
        self.assertEqual(self._thresholds({}, "6x:sampleA,5x:sampleB", "5x:sampleC", self.AFFECTED),
                         {"LongAllele": 6, "ShortAllele": 5})

    def test_no_outliers_and_entries_that_do_not_parse(self):
        self.assertIsNone(self._thresholds({}, "", "", self.AFFECTED))
        # Entries that do not parse are skipped rather than raising.
        self.assertIsNone(self._thresholds({}, "garbage", "", self.AFFECTED))

    def test_solved_samples_set_the_bar_and_unknown_samples_are_candidates(self):
        affected = {"SOLVED": "affected", "UNKNOWN": "unknown", "AFFECTED": "affected"}
        analysis = {"SOLVED": "solved", "AFFECTED": "unsolved"}
        # A solved affected sample is part of the bar, so the unsolved one below it earns nothing.
        self.assertIsNone(self._thresholds({}, "60x:SOLVED,50x:AFFECTED", "", affected, analysis))
        # A sample of unknown status is a candidate, and so is one with no metadata row.
        self.assertEqual(self._thresholds({}, "60x:UNKNOWN,50x:AFFECTED", "", affected, analysis), {"LongAllele": 60})
        self.assertEqual(self._thresholds({}, "6x:sampleZ", "", {}), {"LongAllele": 6})

    def test_a_threshold_is_set_exactly_when_the_build_counts_a_sample(self):
        # The build's #Aff count (find_affected_unsolved_samples_above_unaffected) and the
        # long-allele threshold have to agree on every locus: the threshold is the largest allele
        # the build counts, and there is none when the build counts no sample.
        affected = {"A1": "affected", "A2": "affected", "U1": "unaffected", "S1": "affected", "K1": "unknown"}
        analysis = {"A1": "unsolved", "A2": "partially solved", "S1": "probably solved"}
        row_dict = {"HPRC256_99thPercentile": 8}
        for long_outliers in ("10x:A1,9x:U1", "10x:U1,9x:A1", "10x:S1,9x:A1", "10x:K1,9x:A1", "10x:A1,10x:U1",
                              "9x:A2,7x:A1", "8x:A1", "12x:NOT_IN_THE_TABLE,11x:S1"):
            build_samples = analysis_columns.find_affected_unsolved_samples_above_unaffected(
                {"OutlierSampleIds_AllAlleles": long_outliers, **row_dict}, affected, analysis,
                {"AllAlleles": analysis_columns.parse_outlier_entries(long_outliers)},
                self.LONG_READ_GATE)[0]
            thresholds = self._thresholds(row_dict, long_outliers, "", affected, analysis)
            self.assertEqual(thresholds and thresholds["LongAllele"],
                             max(build_samples.values()) if build_samples else None, long_outliers)


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
        phenotype_description TEXT, LocusId TEXT,
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
    """Display keeps "possibly affected" rather than collapsing it to "affected"."""

    def test_possibly_affected_preserved(self):
        sample_lookup = {"sA": {"affected_status": "possibly affected", "analysis_status": "unsolved",
                                "family_id": "f1", "sex": "male"}}
        parsed = results_server.parse_outlier_samples("60x:sA", sample_lookup)
        self.assertEqual((parsed[0]["affected_status"], parsed[0]["analysis_status"]),
                         ("Possibly Affected", "Unsolved"))

    def test_a_sample_without_a_row_is_unknown(self):
        parsed = results_server.parse_outlier_samples("60x:sB", {})
        self.assertEqual((parsed[0]["affected_status"], parsed[0]["analysis_status"]), ("Unknown", "Unknown"))

    def test_every_blank_status_displays_as_unknown(self):
        # An "NA" cell is no status, so the locus-detail page shows "Unknown" rather than "Na".
        for value in ("", "nan", "none", "na", "n/a", "null", "NA", " N/A ", "Null"):
            described = results_server.describe_sample("S1", {"S1": {"affected_status": value}})
            self.assertEqual(described["affected_status"], "Unknown", repr(value))


class SchemaAndUnknownCategoryTests(ServerTestCase):
    """schema sample_ids fall back to swim_plot, and Unknown-category rows are served."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_unknown_category(cls.db_path)
        # No sample_table -> the sample dropdown must fall back to swim_plot IDs.
        cls.configure_server(cls.db_path)

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


class SchemaSampleCountTests(ServerTestCase):
    """/api/v1/schema reports the sample count build_database recorded in the database."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_sample_count(cls.db_path, 2657)
        cls.configure_server(cls.db_path)

    def test_schema_reports_total_samples(self):
        data = self.client.get("/api/v1/schema").get_json()
        self.assertEqual(data["total_samples"], 2657)
        self.assertEqual(data["total_loci"], 1)

    def test_configure_app_reads_the_count_from_the_database(self):
        # The schema value comes from configure_app's one read at startup, not from a per-request
        # query, so pin the config entry the endpoint reads.
        self.assertEqual(results_server.app.config["TOTAL_SAMPLES"], 2657)


class PermissiveSampleTableTests(unittest.TestCase):
    """load_sample_table accepts minimal metadata and reads it the way the build does."""

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
        # Absent statuses are falsy (treated as not-unaffected downstream) and display as Unknown.
        self.assertFalse(affected_lookup["S1"])
        self.assertFalse(analysis_lookup["S1"])
        described = results_server.describe_sample("S1", sample_rows)
        self.assertEqual((described["affected_status"], described["analysis_status"]), ("Unknown", "Unknown"))

    def test_statuses_are_shown_as_entered_with_the_analysis_status_remapped(self):
        # The sample row keeps "Possibly Affected" as entered and holds the analysis status after
        # read_sample_metadata's remapping, so the pages show what the analysis logic reads.
        path = self._write("sample_id\taffected_status\tanalysis_status\nS1\tPossibly Affected\tRNCC\n")
        sample_rows, _affected, _analysis = results_server.load_sample_table(path)
        described = results_server.describe_sample("S1", sample_rows)
        self.assertEqual((described["affected_status"], described["analysis_status"]),
                         ("Possibly Affected", "Unsolved"))

    def test_case_insensitive_id_and_extra_columns_preserved(self):
        path = self._write("Sample ID\tancestry\nS1\tEUR\n")
        sample_rows, _affected, _analysis = results_server.load_sample_table(path)
        self.assertIn("S1", sample_rows)
        self.assertEqual(sample_rows["S1"]["ancestry"], "EUR")  # extra column preserved

    def test_extra_columns_come_from_the_first_row_of_a_duplicated_sample(self):
        path = self._write("sample_id\tFamilyID\tancestry\nS1\tF1\tEUR\nS1\tF9\tAFR\n")
        sample_rows, _affected, _analysis = results_server.load_sample_table(path)
        self.assertEqual(sample_rows["S1"]["family_id"], "F1")
        self.assertEqual(sample_rows["S1"]["ancestry"], "EUR")
        self.assertNotIn("FamilyID", sample_rows["S1"])

    def test_phenotype_strips_only_leading_na_prefix(self):
        path = self._write(
            "sample_id\taffected_status\tphenotype_description\n"
            "S1\tAffected\tNA; seizures; NA; ataxia\n")
        sample_rows, _affected, _analysis = results_server.load_sample_table(path)
        # Only the leading 'NA; ' is removed (matches the build), not the inner one.
        self.assertEqual(sample_rows["S1"]["phenotype_description"], "seizures; NA; ataxia")

    def test_values_the_build_rejects_are_rejected(self):
        # The same file is forwarded to the build, which refuses an unknown status or sex.
        with self.assertRaises(ValueError):
            results_server.load_sample_table(self._write("sample_id\taffected_status\nS1\tmaybe\n"))
        with self.assertRaises(ValueError):
            results_server.load_sample_table(self._write("sample_id\tsex\nS1\tF\n"))

    def test_analysis_status_codes_and_blanks_are_normalized(self):
        path = self._write("sample_id\tanalysis_status\nS1\tRNCC\nS2\ts_kgfp\nS3\tNA\n")
        _rows, _affected, analysis_lookup = results_server.load_sample_table(path)
        self.assertEqual(analysis_lookup, {"S1": "unsolved", "S2": "solved", "S3": "unknown"})


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


class ExportPhenotypeScoreQueryCountTests(ServerTestCase):
    """Both enriching exports bulk-load the phenotype scores once, not once per locus."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_phenotype_scores(cls.db_path)
        cls.configure_server(cls.db_path)

    def _score_queries(self, recorded_queries):
        """The recorded statements that read the per-outlier phenotype-score table."""
        return [sql for sql in recorded_queries if "per_outlier_phenotype_scores" in sql]

    def test_expansion_hunter_export_reads_the_score_table_once(self):
        # The scores are loaded in one query for the whole export, not with a
        # "WHERE locus_id = ?" query per locus.
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
        recorded_queries = _record_queries(self)
        response = self.client.get("/api/v1/export?outlier_type=all&format=json")
        # The body is gzipped, so it is not decodable as text for an assertion message.
        self.assertEqual(response.status_code, 200)
        data = json.loads(gzip.decompress(response.get_data()))
        self.assertEqual(data["metadata"]["total_loci"], 2)
        self.assertEqual(len(self._score_queries(recorded_queries)), 1)

    def test_only_the_exported_loci_scores_are_loaded(self):
        # A score table can hold millions of rows, so only the exported loci's scores are loaded.
        conn = duckdb_compat.connect(self.db_path, read_only=True)
        try:
            scores = results_server.fetch_per_outlier_phenotype_scores_of_loci(
                conn, "SELECT * FROM loci WHERE LocusId = ? ORDER BY LocusId", ["chr2-200-209-AAG"])
        finally:
            conn.close()
        self.assertEqual(set(scores), {"chr2-200-209-AAG"})
        self.assertEqual(scores["chr2-200-209-AAG"]["AllAlleles"][0]["sample_id"], "sampleC")

    def test_the_exports_do_not_read_the_aggregate_score_table(self):
        # The exports list per-outlier scores only, so the per-locus aggregates are not loaded.
        for export_format in ("json", "expansion_hunter"):
            recorded_queries = _record_queries(self)
            response = self.client.get(f"/api/v1/export?outlier_type=all&format={export_format}")
            self.assertEqual(response.status_code, 200)
            response.get_data()
            self.assertFalse([sql for sql in recorded_queries if "per_locus_phenotype_scores" in sql],
                             export_format)

    def test_the_locus_detail_page_carries_both_kinds_of_score(self):
        scores = self.client.get("/api/v1/loci/chr1-100-110-AT").get_json()["phenotype_scores"]
        self.assertEqual([s["sample_id"] for s in scores["per_outlier"]["AllAlleles"]], ["sampleA"])
        self.assertEqual(scores["per_outlier"]["AllAlleles"][0]["gene_phenotype_similarity"], 0.75)
        self.assertEqual(scores["per_locus"]["AllAlleles"]["max_gene_phenotype_similarity"], 0.75)
        self.assertEqual(scores["per_locus"]["AllAlleles"]["qualifying_sample_ids"], "sampleA")
        # A locus with per-outlier scores and no aggregate row still gets both keys
        scores = self.client.get("/api/v1/loci/chr2-200-209-AAG").get_json()["phenotype_scores"]
        self.assertEqual(scores["per_locus"], {})
        self.assertIsNone(scores["per_outlier"]["AllAlleles"][0]["pairwise_similarity_to_next"])

    def test_a_filtered_json_export_carries_its_loci_scores(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=json&gene_symbol=GENE1")
        self.assertEqual(response.status_code, 200)
        loci = json.loads(gzip.decompress(response.get_data()))["loci"]
        self.assertEqual([locus["LocusId"] for locus in loci], ["chr1-100-110-AT"])
        samples = {s["Sample ID"]: s for s in loci[0]["Samples"]["AllAlleleOutliers"]}
        self.assertEqual(samples["sampleA"]["Gene-Pheno Sim"], 0.75)


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


class PhenotypeKeywordFallbackColumnSubsetTests(ServerTestCase):
    """The fallback narrowed to the phenotype columns a database actually declares.

    Without a swim_plot table the filter reads the summary columns, and it must name only the
    ones that exist: naming SecondAffectedPhenotype_AllAlleles here would be a SQL error on
    every request carrying a phenotype keyword.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_only_the_first_phenotype_column(cls.db_path)
        cls.configure_server(cls.db_path)

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
        self.assertEqual(self.get_locus_ids("phenotype_keyword=seizures"), {"chr1-100-110-AT"})

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


class PhenotypeKeywordWithoutAnyPhenotypeColumnTests(ServerTestCase):
    """A phenotype keyword matches nothing when the database records no phenotypes at all."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_without_phenotype_columns(cls.db_path)
        cls.configure_server(cls.db_path)

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


def _build_db_with_family_counts(path):
    """Two loci whose #Aff and #Fam counts put them in opposite orders."""
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER, Motif TEXT, MotifSize INTEGER,
        gene_region_rank INTEGER, FirstAffectedAlleleSize_AllAlleles INTEGER,
        NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER,
        NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles INTEGER,
        OutlierSampleIds_AllAlleles TEXT)""")
    conn.executemany("INSERT INTO loci VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", [
        # Three samples of one family
        ("chr1-100-110-AT", "chr1", 100, 110, "AT", 2, 1, 9, 3, 1, "9x:sampleA,8x:sampleB,7x:sampleC"),
        # Two samples of two families
        ("chr2-200-209-AAG", "chr2", 200, 209, "AAG", 3, 1, 9, 2, 2, "9x:sampleD,8x:sampleE"),
    ])
    conn.commit()
    conn.close()


class FamilyCountTests(ServerTestCase):
    """Each locus carries its #Fam count, and sort_by=family_count orders by it."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_family_counts(cls.db_path)
        cls.configure_server(cls.db_path)

    def _loci(self, query=""):
        response = self.client.get(f"/api/v1/loci?outlier_type=all{query}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return response.get_json()["results"]

    def test_each_locus_carries_the_family_count(self):
        self.assertEqual({r["LocusId"]: r["NumAffectedUnsolvedFamiliesAboveUnaffected"] for r in self._loci()},
                         {"chr1-100-110-AT": 1, "chr2-200-209-AAG": 2})

    def test_sort_by_family_count(self):
        self.assertEqual([r["LocusId"] for r in self._loci()], ["chr1-100-110-AT", "chr2-200-209-AAG"])
        self.assertEqual([r["LocusId"] for r in self._loci("&sort_by=family_count")],
                         ["chr2-200-209-AAG", "chr1-100-110-AT"])

    def test_the_tsv_export_has_the_family_count_column(self):
        response = self.client.get("/api/v1/export?outlier_type=all&format=tsv&sort_by=family_count")
        self.assertEqual(response.status_code, 200)
        lines = gzip.decompress(response.get_data()).decode().splitlines()
        rows = [dict(zip(lines[0].split("\t"), line.split("\t"))) for line in lines[1:]]
        self.assertEqual([(r["LocusId"], r["NumAffectedUnsolvedFamiliesAboveUnaffected"]) for r in rows],
                         [("chr2-200-209-AAG", "2"), ("chr1-100-110-AT", "1")])

    def test_a_filter_on_a_column_the_database_lacks_matches_nothing_instead_of_failing(self):
        # This database has only the AllAlleles outlier columns, and no NumRepeatsInReference,
        # gene_region or population columns
        for query, expected_total in (
                ("outlier_type=all&min_expansion=1", 0),
                ("outlier_type=all&gene_regions=cds", 0),
                ("outlier_type=all&require_above_population=all&include_loci_without_population_data=true", 2),
                # An exclusion has nothing to exclude without the column
                ("outlier_type=all&exclude_gene_regions=cds", 2),
                ("outlier_type=short&min_expansion=1", 0),
                ("outlier_type=short&sample_id_keyword=sampleA", 0),
                ("outlier_type=short&sample_id_like=sample", 0),
                ("outlier_type=short&apply_pathogenic_threshold=true", 0),
                ("outlier_type=short&require_above_population=all&include_loci_without_population_data=true", 0)):
            response = self.client.get(f"/api/v1/loci?{query}")
            self.assertEqual(response.status_code, 200, (query, response.get_data(as_text=True)))
            self.assertEqual(response.get_json()["total"], expected_total, query)


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
    """Without locus_annotations, no locus matches a known disease locus."""

    @classmethod
    def setUpClass(cls):
        cls.degraded = _load_results_server_without_locus_annotations()

    def test_the_fallback_really_is_a_stub(self):
        # If this fails the import did not raise and the assertions below would be testing the
        # real locus_annotations functions rather than the fallbacks.
        self.assertEqual(self.degraded.load_known_disease_loci("ignored"), ({}, {}))

    def test_a_locus_matches_no_known_disease_locus(self):
        lookups = {"disease_trees": {}, "strchive_trees": {}}
        row = {"LocusId": "chr1-100-110-AT", "Chrom": "chr1", "Start0Based": 100, "End1Based": 110, "Motif": "AT"}
        self.assertIsNone(self.degraded.compute_known_disease_info(row, lookups))
        self.assertIsNone(self.degraded.known_disease_pathogenic_min("chr1", 100, 110, "AT", {}, {}))


class RequestCheckTests(MinimalDatabaseTestCase):
    """No response carries CORS headers, and the Host check stops DNS rebinding whatever the bind
    address, while still answering other machines that use this one's address or names."""

    THIS_MACHINES_HOST_NAMES = {"labhost", "labhost.local", "labhost.example.org", "labhost.example.org.local"}

    def setUp(self):
        self.addCleanup(results_server.app.config.pop, "THIS_MACHINES_HOST_NAMES", None)
        results_server.app.config["THIS_MACHINES_HOST_NAMES"] = self.THIS_MACHINES_HOST_NAMES

    def _set_bind_host(self, bind_host):
        self.addCleanup(results_server.app.config.pop, "BIND_HOST", None)
        results_server.app.config["BIND_HOST"] = bind_host

    def _rebound_requests(self):
        """Read the results and write a note the way a page on rebound.example would.

        Returns:
            Tuple (read status, write status).
        """
        read = self.client.get("/api/v1/schema", base_url="http://rebound.example:5050")
        write = self.client.put("/api/v1/annotations/chr2-200-209-AAG/note", json={"note_text": "rebound"},
                                base_url="http://rebound.example:5050",
                                headers={"Origin": "http://rebound.example:5050"})
        self.assertNotIn("chr2-200-209-AAG", results_server.app.config["ANNOTATIONS"]["notes"])
        return read.status_code, write.status_code

    def test_no_response_carries_cors_headers(self):
        response = self.client.get("/api/v1/loci?outlier_type=all", headers={"Origin": "https://other.example"})
        self.assertEqual(response.status_code, 200)
        self.assertNotIn("Access-Control-Allow-Origin", response.headers)
        preflight = self.client.options("/api/v1/loci", headers={
            "Origin": "https://other.example", "Access-Control-Request-Method": "GET"})
        self.assertNotIn("Access-Control-Allow-Origin", preflight.headers)

    def test_loopback_host_names_are_answered(self):
        for bind_host in ("127.0.0.1", "localhost", "::1"):
            self._set_bind_host(bind_host)
            for host in ("localhost", "localhost:5050", "LOCALHOST:5050", "127.0.0.1:5050", "[::1]:5050"):
                for path in ("/", "/api/v1/schema"):
                    response = self.client.get(path, base_url=f"http://{host}")
                    self.assertEqual(response.status_code, 200, (bind_host, host, path))

    def test_another_host_name_is_refused_whatever_the_bind_address(self):
        # A page whose own domain was pointed at this machine connects from it, and sends that
        # domain as the Host and as the Origin of a write, so the two match and only the Host
        # check refuses it. Without a recorded bind address, restrict_writes treats the server as
        # bound to a non-loopback address.
        for bind_host in (None, "127.0.0.1", "localhost", "::1", "0.0.0.0", "192.168.1.5"):
            if bind_host is None:
                results_server.app.config.pop("BIND_HOST", None)
            else:
                self._set_bind_host(bind_host)
            self.assertEqual(self._rebound_requests(), (403, 403), bind_host)
            response = self.client.get("/", base_url="http://rebound.example")
            self.assertEqual(response.status_code, 403, bind_host)
            self.assertIn("labhost.local", response.get_json()["detail"])

    def test_names_that_only_look_like_this_machine_are_refused(self):
        for host in ("localhost.rebound.example", "labhost.rebound.example", "127.0.0.1.nip.io",
                     "rebound.example."):
            response = self.client.get("/api/v1/schema", base_url=f"http://{host}")
            self.assertEqual(response.status_code, 403, host)
        # A Host value that will not parse has no name, which is refused too
        self.assertEqual(results_server.parse_host_name("[::1"), "")
        self.assertEqual(results_server.parse_host_name("[::1]:5050"), "::1")

    def test_ip_addresses_and_this_machines_names_are_answered_when_bound_to_another_address(self):
        for bind_host in ("0.0.0.0", "192.168.1.5"):
            self._set_bind_host(bind_host)
            for host in ("192.168.1.5:5050", "10.0.0.7", "[fe80::1]:5050", "localhost:5050", "labhost:5050",
                         "LabHost.local:5050", "labhost.example.org", "labhost.example.org.:5050"):
                response = self.client.get("/api/v1/schema", base_url=f"http://{host}")
                self.assertEqual(response.status_code, 200, (bind_host, host))
            # Writes from another machine are still refused, by restrict_writes
            write = self.client.put("/api/v1/annotations/chr2-200-209-AAG/note", json={"note_text": "remote"},
                                    base_url="http://labhost.local:5050",
                                    environ_overrides={"REMOTE_ADDR": "192.168.1.9"})
            self.assertEqual(write.status_code, 403, bind_host)
            self.assertIn("non-loopback address", write.get_json()["detail"])

    def test_writes_from_other_machines_are_answered_with_allow_remote_writes(self):
        self._set_bind_host("0.0.0.0")
        self.addCleanup(results_server.app.config.pop, "ALLOW_REMOTE_WRITES", None)
        results_server.app.config["ALLOW_REMOTE_WRITES"] = True
        self.addCleanup(self.client.delete, "/api/v1/annotations/chr2-200-209-AAG/note")
        write = self.client.put("/api/v1/annotations/chr2-200-209-AAG/note", json={"note_text": "remote"},
                                base_url="http://labhost.local:5050",
                                environ_overrides={"REMOTE_ADDR": "192.168.1.9"})
        self.assertEqual(write.status_code, 200)
        # A rebinding page is refused all the same
        results_server.app.config["ANNOTATIONS"]["notes"].pop("chr2-200-209-AAG")
        self.assertEqual(self._rebound_requests(), (403, 403))

    def test_this_machines_host_names(self):
        with mock.patch.object(results_server.socket, "gethostname", return_value="LabHost"), \
                mock.patch.object(results_server.socket, "getfqdn", return_value="labhost.example.org."):
            self.assertEqual(results_server.find_this_machines_host_names(),
                             {"labhost", "labhost.local", "labhost.example.org", "labhost.example.org.local"})
        with mock.patch.object(results_server.socket, "gethostname", return_value="labhost.local"), \
                mock.patch.object(results_server.socket, "getfqdn", return_value="labhost.local"):
            self.assertEqual(results_server.find_this_machines_host_names(), {"labhost.local"})

    def test_a_named_input_file_that_does_not_exist_stops_the_server(self):
        existing_path = os.path.join(self.tmpdir, "results.duckdb")
        results_server.exit_if_a_supplied_input_file_is_missing(
            {"--sample-table": None, "--known-loci-json": existing_path, "--strchive-loci-json": None})
        with mock.patch("builtins.print") as printed, self.assertRaises(SystemExit) as raised:
            results_server.exit_if_a_supplied_input_file_is_missing(
                {"--sample-table": os.path.join(self.tmpdir, "no_such_samples.tsv"), "--known-loci-json": None})
        self.assertEqual(raised.exception.code, 1)
        self.assertIn("--sample-table file not found", printed.call_args.args[0])

    def test_filters_applied_is_the_same_in_the_list_and_the_json_export(self):
        query = "outlier_type=all&min_pli=0.5&motif_size=2&search=GENE1&reference_region=chr1:100-110&sort_by=size"
        expected = {"min_pli": 0.5, "motif_size": "2", "search": "GENE1", "reference_region": "chr1:100-110"}
        self.assertEqual(self.client.get(f"/api/v1/loci?{query}").get_json()["filters_applied"], expected)
        export = self.client.get(f"/api/v1/export?{query}&format=json")
        self.assertEqual(json.loads(gzip.decompress(export.get_data()))["metadata"]["filters_applied"], expected)


MENDELIAN_BINS = ["autosome", "chrX", "chrY"] + [f"motif_{category}" for category in (
    "1bp", "2bp", "3bp", "4bp", "5bp", "6bp", "7_24bp", "25plusbp")]


def _build_db_with_mendelian_violations(path):
    """Five trios from center ABC and one from XYZ, 1 violation in 100 calls per bin unless noted.

    ABC_1 has 20 autosome violations, above the warning threshold. ABC_2 has no chrX counts and
    no 1bp violation count, and ABC_3 no autosome total (NULLs).
    """
    conn = duckdb_compat.connect(path)
    conn.execute("CREATE TABLE loci (LocusId TEXT, Motif TEXT, KnownDiseaseLocus TEXT)")
    column_definitions = ", ".join(f"{name}_violations BIGINT, {name}_total BIGINT" for name in MENDELIAN_BINS)
    conn.execute(f"CREATE TABLE mendelian_violations (sample_id TEXT, {column_definitions})")
    overrides = {
        "ABC_1": {"autosome_violations": 20},
        "ABC_2": {"chrX_violations": None, "chrX_total": None, "motif_1bp_violations": None},
        "ABC_3": {"autosome_violations": 5, "autosome_total": None},
    }
    for sample_id in ("ABC_1", "ABC_2", "ABC_3", "ABC_4", "ABC_5", "XYZ_1"):
        values = {column: value for name in MENDELIAN_BINS
                  for column, value in ((f"{name}_violations", 1), (f"{name}_total", 100))}
        values.update(overrides.get(sample_id, {}))
        conn.execute(f"INSERT INTO mendelian_violations VALUES (?, {', '.join('?' * len(values))})",
                     [sample_id] + list(values.values()))
    conn.commit()
    conn.close()


class MendelianViolationTests(ServerTestCase):
    """The startup warnings and the QC2 page read the violation counts the same way, NULLs included."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_mendelian_violations(cls.db_path)
        cls.configure_server(cls.db_path)

    def test_the_startup_warnings_flag_only_the_sample_above_the_threshold(self):
        warnings = results_server.app.config["MENDELIAN_WARNINGS"]
        self.assertEqual(set(warnings), {"ABC_1"})
        self.assertEqual(warnings["ABC_1"]["autosome_rate"], 0.2)
        self.assertEqual(warnings["ABC_1"]["motif_rates"]["1bp"], {"rate": 0.01, "violations": 1, "total": 100})

    def test_the_qc2_page_reads_missing_counts_as_zero(self):
        data = self.client.get("/api/v1/qc2_data").get_json()
        by_chrom = {(r["sample_id"], r["bin"]): r for r in data["by_chrom"]}
        self.assertEqual(by_chrom[("ABC_2", "chrX")], {"sample_id": "ABC_2", "bin": "chrX", "fraction": 0,
                                                       "violations": 0, "total": 0})
        self.assertEqual(by_chrom[("ABC_3", "autosome")]["fraction"], 0)
        by_motif = {(r["sample_id"], r["bin"]): r for r in data["by_motif"]}
        self.assertEqual((by_motif[("ABC_2", "1bp")]["violations"], by_motif[("ABC_2", "1bp")]["fraction"]), (0, 0))
        self.assertEqual(by_motif[("ABC_4", "7-24bp")]["fraction"], 0.01)
        self.assertEqual(len(data["by_motif"]), 6 * 8)

    def test_only_centers_with_enough_trios_are_reported(self):
        data = self.client.get("/api/v1/qc2_data").get_json()
        by_center_chrom = {(r["center"], r["bin"]): r for r in data["by_center_chrom"]}
        self.assertEqual({center for center, _ in by_center_chrom}, {"ABC"})
        self.assertEqual(by_center_chrom[("ABC", "autosome")],
                         {"center": "ABC", "bin": "autosome", "fraction": 28 / 400, "violations": 28, "total": 400,
                          "trios": 5})
        self.assertEqual(by_center_chrom[("ABC", "chrX")]["total"], 400)
        self.assertEqual(len(data["by_center_motif"]), 8)

    def test_the_qc2_page_gets_the_trio_minimum_it_names(self):
        render_arguments = self.record_template_render_arguments()
        self.assertEqual(self.client.get("/qc2").status_code, 200)
        self.assertEqual(render_arguments["min_trios_per_center"], results_server.MIN_TRIOS_PER_CENTER)


def _build_db_for_min_expansion(path):
    """Loci whose second affected allele, per sample and per family, sits just above the unaffected one.

    A: reference 5, largest affected 6, second 5, unaffected 4, so it passes min_expansion=1 and
       require_above_unaffected=second on their own and together.
    B: the second affected allele ties the unaffected one.
    C: the largest affected allele equals the reference.
    """
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER, Motif TEXT, MotifSize INTEGER,
        OutlierSampleIds_AllAlleles TEXT, NumRepeatsInReference INTEGER,
        FirstAffectedAlleleSize_AllAlleles INTEGER, SecondAffectedAlleleSize_AllAlleles INTEGER,
        SecondAffectedAlleleSize_AllAlleles_ByFamily INTEGER, FirstUnaffectedAlleleSize_AllAlleles INTEGER)""")
    conn.executemany("INSERT INTO loci VALUES (?, 'chr1', ?, ?, 'AT', 2, '1x:sampleA', ?, ?, ?, ?, ?)", [
        ("A", 100, 110, 5, 6, 5, 5, 4),
        ("B", 200, 210, 5, 6, 4, 4, 4),
        ("C", 300, 310, 12, 12, 11, 11, 3),
    ])
    conn.commit()
    conn.close()


class MinExpansionTests(ServerTestCase):
    """min_expansion is a margin above the reference allele only, not above the unaffected one."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_for_min_expansion(cls.db_path)
        cls.configure_server(cls.db_path)

    def test_each_filter_on_its_own(self):
        self.assertEqual(self.get_locus_ids("min_expansion=1"), {"A", "B"})
        self.assertEqual(self.get_locus_ids("require_above_unaffected=second"), {"A", "C"})
        self.assertEqual(self.get_locus_ids("require_families_above_unaffected=second"), {"A", "C"})

    def test_combined_filters_keep_the_loci_that_pass_both(self):
        for unaffected_filter in ("require_above_unaffected", "require_families_above_unaffected"):
            for min_expansion in (0, 1):
                query = f"{unaffected_filter}=second&min_expansion={min_expansion}"
                self.assertEqual(self.get_locus_ids(query),
                                 self.get_locus_ids(f"{unaffected_filter}=second")
                                 & self.get_locus_ids(f"min_expansion={min_expansion}"),
                                 query)
            self.assertEqual(self.get_locus_ids(f"{unaffected_filter}=second&min_expansion=1"), {"A"})

    def test_a_missing_unaffected_column_matches_nothing(self):
        saved = results_server.app.config["DB_COLUMNS_SET"]
        self.addCleanup(results_server.app.config.__setitem__, "DB_COLUMNS_SET", saved)
        results_server.app.config["DB_COLUMNS_SET"] = saved - {"FirstUnaffectedAlleleSize_AllAlleles"}
        self.assertEqual(self.get_locus_ids("require_families_above_unaffected=second"), set())


class SwimPlotFilterTests(MinimalDatabaseTestCase):
    """The swim plot reads its parameters the way /api/v1/loci does.

    Its locus filters keep the rows of exactly the loci the results table returns, while the
    sample-level filters and SWIM_PLOT_ALLELE_LEVEL_FILTERS pick rows: each row is one allele.
    """

    def _swim_rows(self, query):
        response = self.client.get(f"/api/v1/swim_plot_data?outlier_type=all&{query}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {(entry["LocusId"], entry["sample_id"]) for entry in response.get_json()["data"]}

    def test_locus_filters_keep_the_loci_the_results_table_returns(self):
        self.client.post("/api/v1/annotations/chr2-200-209-AAG/tags", json={"tag": "swim"})
        self.addCleanup(self.client.delete, "/api/v1/annotations/chr2-200-209-AAG/tags/swim")
        for query in ("gene_regions=cds", "exclude_gene_regions=cds", "motif=TA", "motif=GAA", "motif_size=3",
                      "motif_size=-2", "mendelian_genes_only=true", "known_motifs_only=1", "min_pli=0.5",
                      "gene_symbol=GENE2", "gene_id=GENE1", "locus_id=chr1-100-110-AT", "chrom=chr2",
                      "reference_region=chr1:100-110", "search=GENE2", "search=chr1:105", "search=nothing",
                      "search=CHR2-200-209-aag", "tag=swim", "tag=no_such_tag", "known_loci_only=true",
                      "max_variation_cluster_size_diff=5", "min_repeats_threshold=5"):
            self.assertEqual({locus_id for locus_id, _ in self._swim_rows(query)}, self.get_locus_ids(query), query)

    def test_invalid_values_are_rejected_as_by_the_results_table(self):
        for query in ("gene_regions=foo", "exclude_gene_regions=foo", "motif=AXT", "motif_size=a-b",
                      "min_pli=high", "min_sigma_percentile=2", "max_variation_cluster_size_diff=abc",
                      "min_expansion=-1", "require_above_population=some", "population_metric=mean",
                      "search=chr1:200-100"):
            for path in ("/api/v1/swim_plot_data", "/api/v1/loci"):
                response = self.client.get(f"{path}?outlier_type=all&{query}")
                self.assertEqual(response.status_code, 400, (path, query))
        self.assertEqual(self.client.get("/api/v1/swim_plot_data").status_code, 400)

    def test_sample_level_filters_keep_the_rows_of_the_matching_samples(self):
        # The results table keeps the locus; the swim plot draws only that sample's allele there.
        self.assertEqual(self.get_locus_ids("sample_id_keyword=sampleB"), {"chr1-100-110-AT"})
        self.assertEqual(self._swim_rows("sample_id_keyword=SAMPLEB"), {("chr1-100-110-AT", "sampleB")})
        self.assertEqual(self._swim_rows("sample_id_like=ple_"), set())
        self.assertEqual(self._swim_rows("sample_id_like=pleC"), {("chr2-200-209-AAG", "sampleC")})
        self.assertEqual(self._swim_rows("phenotype_keyword=SEIZ"), {("chr1-100-110-AT", "sampleA")})

    def test_the_pathogenic_threshold_compares_each_allele(self):
        # chr1's largest affected allele (6) reaches the minimum, so the results table keeps the
        # locus, but sampleB's 5 does not, so the swim plot drops that allele.
        saved = results_server.app.config["KNOWN_DISEASE_LOCUS_THRESHOLDS"]
        self.addCleanup(results_server.app.config.__setitem__, "KNOWN_DISEASE_LOCUS_THRESHOLDS", saved)
        results_server.app.config["KNOWN_DISEASE_LOCUS_THRESHOLDS"] = {"chr1-100-110-AT": 6}
        self.assertEqual(self.get_locus_ids("apply_pathogenic_threshold=true"), {"chr1-100-110-AT"})
        self.assertEqual(self._swim_rows("apply_pathogenic_threshold=1"), {("chr1-100-110-AT", "sampleA")})

    def test_require_affected_above_unaffected(self):
        # sampleB is unaffected, and sampleC has no unaffected allele to be compared with.
        for value in ("1", "true"):
            self.assertEqual(self._swim_rows(f"require_affected_above_unaffected={value}"),
                             {("chr1-100-110-AT", "sampleA"), ("chr2-200-209-AAG", "sampleC")})

    def test_require_affected_above_unaffected_follows_the_builds_rule(self):
        # The build counts a solved affected sample with the unaffected ones, and a sample of
        # unknown status with the candidates. SOLVED is affected and solved, UNKNOWN has blank
        # statuses, AFFECTED is affected and unsolved, UNAFFECTED is unaffected.
        sample_table = os.path.join(self.tmpdir, "statuses.tsv")
        with open(sample_table, "w") as output_file:
            output_file.write("sample_id\taffected_status\tanalysis_status\n"
                              "SOLVED\taffected\tsolved\nUNKNOWN\t\t\n"
                              "AFFECTED\taffected\tunsolved\nUNAFFECTED\tunaffected\tunaffected\n")
        sample_lookup, affected_lookup, analysis_lookup = input_tables.read_sample_metadata(sample_table)
        records = [
            # UNKNOWN is above SOLVED; AFFECTED is not
            {"LocusId": "chr1-100-110-AT", "OutlierSampleIds_AllAlleles": "60x:UNKNOWN,55x:SOLVED,50x:AFFECTED"},
            # SOLVED is the largest, so no candidate is above it
            {"LocusId": "chr1-200-210-AT", "OutlierSampleIds_AllAlleles": "60x:SOLVED,50x:UNKNOWN"},
            {"LocusId": "chr1-300-310-AT", "OutlierSampleIds_AllAlleles": "40x:AFFECTED,30x:UNAFFECTED"},
        ]
        for record in records:
            start = int(record["LocusId"].split("-")[1])
            record.update({"Chrom": "chr1", "Start0Based": start, "End1Based": start + 10, "Motif": "AT",
                           "MotifSize": 2, "OutlierSampleIds_ShortAlleles": None,
                           "OutlierSampleIds_HemizygousAlleles": None})
        analysis_columns.add_all_outlier_columns(records, sample_lookup, affected_lookup, analysis_lookup)
        db_path = os.path.join(self.tmpdir, "statuses.duckdb")
        conn = duckdb_compat.connect(db_path)
        result_database.write_loci_table(conn, records, list(records[0]))
        result_database.write_swim_plot(conn, swim_plot.generate_swim_plot_table(records, sample_lookup),
                                        columns=swim_plot.SWIM_PLOT_COLUMNS)
        conn.commit()
        conn.close()
        self.configure_server(db_path, sample_table=sample_table)
        self.addCleanup(self.configure_server, self.db_path)

        self.assertEqual(self.get_locus_ids("require_above_unaffected=first"), {"chr1-100-110-AT", "chr1-300-310-AT"})
        swim_rows = self._swim_rows("require_affected_above_unaffected=1")
        self.assertEqual(swim_rows, {("chr1-100-110-AT", "UNKNOWN"), ("chr1-300-310-AT", "AFFECTED")})
        self.assertEqual({locus_id for locus_id, _ in swim_rows}, self.get_locus_ids("require_above_unaffected=first"))

    def test_the_matching_loci_are_looked_up_once_per_request(self):
        recorded_queries = _record_queries(self)
        self._swim_rows("gene_regions=cds&search=GENE1")
        self.assertEqual(len([sql for sql in recorded_queries if "FROM loci" in sql]), 1)
        recorded_queries.clear()
        self._swim_rows("")
        self.assertFalse([sql for sql in recorded_queries if "FROM loci" in sql])

    def test_min_delta_population_compares_the_gate_cohorts_as_the_results_table_does(self):
        # A: allele 20, long-read 99th percentiles 10, TenK10K 30. Its Δpop is 10 against the
        # long-read cohorts but -10 against every cohort, so min Δpop 5 keeps it only when the
        # cohorts are the long-read gate's. B: allele 20 with every percentile at 10.
        db_path = os.path.join(self.tmpdir, "population.duckdb")
        conn = duckdb_compat.connect(db_path)
        conn.execute("""CREATE TABLE loci (
            LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER, Motif TEXT, MotifSize INTEGER,
            OutlierSampleIds_AllAlleles TEXT, FirstAffectedAlleleSize_AllAlleles INTEGER,
            HPRC256_99thPercentile DOUBLE, AoU1027_99thPercentile DOUBLE, TenK10K_99thPercentile DOUBLE)""")
        population_percentiles = {"A": (10, 10, 30), "B": (10, 10, 10)}
        conn.executemany("INSERT INTO loci VALUES (?, 'chr1', ?, ?, 'AT', 2, ?, 20, ?, ?, ?)", [
            (locus_id, index * 100, index * 100 + 10, f"20x:sample{locus_id}", *percentiles)
            for index, (locus_id, percentiles) in enumerate(population_percentiles.items())])
        result_database.write_swim_plot(
            conn,
            [make_swim_plot_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp", allele_size=20,
                                sample_id=f"sample{locus_id}", LocusId=locus_id, Motif="AT", MotifSize=2,
                                HPRC256_99thPercentile=hprc, AoU1027_99thPercentile=aou,
                                TenK10K_99thPercentile=tenk10k)
             for locus_id, (hprc, aou, tenk10k) in population_percentiles.items()],
            columns=swim_plot.SWIM_PLOT_COLUMNS,
        )
        conn.commit()
        conn.close()
        self.configure_server(db_path)
        self.addCleanup(self.configure_server, self.db_path)

        for gate, expected_loci in (("long-read", {"A", "B"}), ("all", {"B"})):
            query = f"require_above_population={gate}&min_repeats_above_population=5"
            self.assertEqual(self.get_locus_ids(query), expected_loci, query)
            self.assertEqual({locus_id for locus_id, _ in self._swim_rows(query)}, expected_loci, query)


def _build_db_with_population_percentiles(path, population_gate_cohorts):
    """One locus whose only outlier, an affected 20, is above HPRC256 and AoU1027 but not TenK10K.

    Args:
        path: Where to write the database.
        population_gate_cohorts: The cohorts to record as the build's population gate, or None to
            record none, which reads back as the default long-read cohorts.
    """
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER, Motif TEXT, MotifSize INTEGER,
        OutlierSampleIds_AllAlleles TEXT, OutlierSampleIds_ShortAlleles TEXT,
        HPRC256_99thPercentile DOUBLE, AoU1027_99thPercentile DOUBLE, TenK10K_99thPercentile DOUBLE)""")
    conn.execute("INSERT INTO loci VALUES ('chr1-100-110-AT', 'chr1', 100, 110, 'AT', 2, '20x:sampleA', '', "
                 "10, 12, 30)")
    if population_gate_cohorts:
        result_database.write_population_gate_cohorts(conn, population_gate_cohorts)
    conn.commit()
    conn.close()


class ExportPopulationGateTests(ServerTestCase):
    """The exports' outlier lists and read-visualization thresholds use the build's population gate."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.sample_table = os.path.join(cls.tmpdir, "samples.tsv")
        with open(cls.sample_table, "w") as output_file:
            output_file.write("sample_id\taffected_status\tanalysis_status\nsampleA\taffected\tunsolved\n")

    def _exports(self, population_gate_cohorts):
        """Return the JSON export's locus and the ExpansionHunter export's entry under a gate."""
        db_path = os.path.join(self.tmpdir, f"results_{len(population_gate_cohorts or ())}.duckdb")
        _build_db_with_population_percentiles(db_path, population_gate_cohorts)
        self.configure_server(db_path, sample_table=self.sample_table)
        json_export = self.client.get("/api/v1/export?outlier_type=all&format=json")
        self.assertEqual(json_export.status_code, 200)
        expansion_hunter_export = self.client.get("/api/v1/export?outlier_type=all&format=expansion_hunter")
        self.assertEqual(expansion_hunter_export.status_code, 200)
        return (json.loads(gzip.decompress(json_export.get_data()))["loci"][0],
                json.loads(expansion_hunter_export.get_data())[0])

    def test_the_default_long_read_gate(self):
        locus, catalog_entry = self._exports(None)
        self.assertEqual([s["Sample ID"] for s in locus["Samples"]["AllAlleleOutliers"]], ["sampleA"])
        self.assertEqual(catalog_entry["PlotReadVisualization"], [{"If": "LongAllele", "Is": ">=", "Threshold": 20}])

    def test_a_gate_that_includes_tenk10k(self):
        # The build left the allele out of #Aff and the phenotype scores, so the exports leave it
        # out too.
        locus, catalog_entry = self._exports(("HPRC256", "AoU1027", "TenK10K"))
        self.assertNotIn("Samples", locus)
        self.assertNotIn("PlotReadVisualization", catalog_entry)


def _build_db_with_many_loci(path, num_loci):
    """num_loci loci, each with an affected outlier above an unaffected one and a phenotype score."""
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT, Chrom TEXT, Start0Based INTEGER, End1Based INTEGER, Motif TEXT, MotifSize INTEGER,
        OutlierSampleIds_AllAlleles TEXT, OutlierSampleIds_ShortAlleles TEXT,
        NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER)""")
    locus_ids = [f"chr1-{index * 100}-{index * 100 + 10}-AT" for index in range(num_loci)]
    conn.executemany("INSERT INTO loci VALUES (?, 'chr1', ?, ?, 'AT', 2, '20x:sampleA,10x:sampleB', '', 1)",
                     [(locus_id, index * 100, index * 100 + 10) for index, locus_id in enumerate(locus_ids)])
    result_database.write_phenotype_tables(
        conn,
        [{"locus_id": locus_id, "sample_id": "sampleA", "outlier_type": "AllAlleles", "allele_size": 20,
          "gene_symbol": "GENE1", "gene_phenotype_similarity": 0.5, "gene_phenotype_overlap_count": 1,
          "n_matching_diseases": 0, "best_matching_disease": None, "best_disease_inheritance": None,
          "pairwise_similarity_to_next": None, "pairwise_shared_count_raw": None,
          "pairwise_shared_count_ic": None, "next_sample_id": None} for locus_id in locus_ids],
        [],
    )
    conn.commit()
    conn.close()


class LargeExpansionHunterExportTests(ServerTestCase):
    """Every locus of an ExpansionHunter export is enriched, however many loci it has."""

    NUM_LOCI = 250

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path = os.path.join(cls.tmpdir, "results.duckdb")
        _build_db_with_many_loci(cls.db_path, cls.NUM_LOCI)
        sample_table = os.path.join(cls.tmpdir, "samples.tsv")
        with open(sample_table, "w") as output_file:
            output_file.write("sample_id\taffected_status\nsampleA\taffected\nsampleB\tunaffected\n")
        cls.configure_server(cls.db_path, sample_table=sample_table)

    def test_every_locus_carries_its_thresholds_and_samples(self):
        recorded_queries = _record_queries(self)
        response = self.client.get("/api/v1/export?outlier_type=all&format=expansion_hunter")
        self.assertEqual(response.status_code, 200)
        catalog = json.loads(response.get_data())
        self.assertEqual(len(catalog), self.NUM_LOCI)
        for entry in catalog:
            self.assertEqual(entry["PlotReadVisualization"], [{"If": "LongAllele", "Is": ">=", "Threshold": 20}],
                             entry["LocusId"])
            self.assertEqual([(s["Sample ID"], s["Gene-Pheno Sim"]) for s in entry["Samples"]["AllAlleleOutliers"]],
                             [("sampleA", 0.5), ("sampleB", None)], entry["LocusId"])
        # The scores are loaded in one query, not one per locus
        self.assertEqual(len([sql for sql in recorded_queries if "per_outlier_phenotype_scores" in sql]), 1)


class KnownDiseaseMatchTests(ServerTestCase):
    """The locus-detail page and the pathogenic-threshold filter take the same match and minimum."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        catalog_path = os.path.join(cls.tmpdir, "catalog.json")
        with open(catalog_path, "w") as output_file:
            json.dump([
                {"LocusId": "ZERO_MIN", "MainReferenceRegion": "chr1:100-110", "RepeatUnit": "AT",
                 "Diseases": [{"Name": "D1", "PathogenicMin": 0}, {"Name": "D2", "PathogenicMin": 40}]},
                {"LocusId": "NO_MIN", "MainReferenceRegion": "chr3:300-310", "RepeatUnit": "CAG",
                 "Diseases": [{"Name": "D3"}]},
            ], output_file)
        strchive_path = os.path.join(cls.tmpdir, "STRchive-loci.json")
        with open(strchive_path, "w") as output_file:
            json.dump(STRCHIVE_TEST_LOCI, output_file)
        cls.disease_trees, cls.strchive_trees = locus_annotations.load_known_disease_loci(
            catalog_path, strchive_filepath=strchive_path)

    def _detail(self, chrom, start, end, motif):
        lookups = {"disease_trees": self.disease_trees, "strchive_trees": self.strchive_trees}
        return results_server.compute_known_disease_info(
            {"LocusId": "query", "Chrom": chrom, "Start0Based": start, "End1Based": end, "Motif": motif}, lookups)

    def _pathogenic_min(self, chrom, start, end, motif):
        return results_server.known_disease_pathogenic_min(chrom, start, end, motif, self.disease_trees,
                                                           self.strchive_trees)

    def test_a_pathogenic_minimum_of_zero_is_kept_by_both(self):
        self.assertEqual(self._detail("chr1", 100, 110, "TA")["pathogenic_min"], 0)
        self.assertEqual(self._pathogenic_min("chr1", 100, 110, "TA"), 0)

    def test_a_catalog_match_without_a_minimum_hides_the_strchive_fallback_from_both(self):
        # STRchive's TEST_REFERENCE has a minimum of 50 at the same locus, but the catalog matched.
        detail = self._detail("chr3", 300, 310, "CAG")
        self.assertEqual((detail["locus_id"], detail["source"], detail["pathogenic_min"]),
                         ("NO_MIN", "variant_catalog", None))
        self.assertIsNone(self._pathogenic_min("chr3", 300, 310, "CAG"))

    def test_a_strchive_only_locus_matches_through_its_pathogenic_motif(self):
        detail = self._detail("chr9", 100, 110, "AAGGG")
        self.assertEqual((detail["locus_id"], detail["source"], detail["pathogenic_min"]),
                         ("TEST_PATHOGENIC", "STRchive", 60))
        self.assertEqual(self._pathogenic_min("chr9", 100, 110, "AAGGG"), 60)

    def test_too_little_overlap_matches_nothing(self):
        # [100, 110) against [100, 125) has a Jaccard index of 0.4
        self.assertIsNone(self._detail("chr1", 100, 125, "AT"))
        self.assertIsNone(self._pathogenic_min("chr1", 100, 125, "AT"))

    def test_a_null_pathogenic_motifs_list_still_matches_by_repeat_unit(self):
        tree = intervaltree.IntervalTree()
        tree.addi(100, 110, data={
            "RepeatUnit": "AT",
            "PathogenicMotifs": None,  # an explicit null in the catalog JSON
            "LocusId": "DISEASE_AT",
            "Diseases": [{"Name": "DiseaseX"}],
        })
        lookups = {"disease_trees": {"1": tree}, "strchive_trees": {}}
        result = results_server.compute_known_disease_info(
            {"LocusId": "q", "Chrom": "chr1", "Start0Based": 100, "End1Based": 110, "Motif": "AT"}, lookups)
        self.assertEqual(result["locus_id"], "DISEASE_AT")

    def test_the_server_matches_through_the_builds_function(self):
        # One definition, so the detail page and the filter cannot drift from the build's
        # KnownDiseaseLocus
        self.assertIs(results_server.find_matching_disease_locus, locus_annotations.find_matching_disease_locus)
        self.assertEqual(locus_annotations.matches_disease_locus("chr9-100-110-AAGGG", self.disease_trees,
                                                                 self.strchive_trees), "TEST_PATHOGENIC")
        self.assertEqual(locus_annotations.matches_disease_locus("chr3-300-310-CAG", self.disease_trees,
                                                                 self.strchive_trees), "NO_MIN")


if __name__ == "__main__":
    unittest.main()
