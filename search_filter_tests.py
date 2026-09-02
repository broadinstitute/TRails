"""Unit tests for the merged Search box and the Max Variation Cluster Size Diff filter.

The Search box replaced four separate filters (Gene Symbol(s), Gene ID, Locus Id(s) and
Reference Region), so the two things worth pinning down are:

  * classify_search_term(), which decides what each comma-separated term is from its shape, and
  * the WHERE clause build_api_query() produces, where terms OR together rather than AND.

A small throwaway loci table stands in for a real results database so the generated SQL
actually runs, including the numeric comparison against the TEXT VariationClusterSizeDiff
column that TRExplorer writes.
"""

import os
import sqlite3
import tempfile
import unittest

import results_server


# LocusId, Chrom, Start0Based, End1Based, CanonicalMotif, gene symbol, Ensembl id,
# VariationClusterSizeDiff (TEXT, as TRExplorer writes it), plus the columns every query touches.
TEST_LOCI = [
    ("1-100-110-AT", "chr1", 100, 110, "AT", "FMR1", "ENSG00000102081", "13",
     1, 1, "CDS", "6x:sampleA"),
    ("1-500-600-CAG", "chr1", 500, 600, "AGC", "HTT", "ENSG00000197386", "148",
     1, 1, "CDS", "9x:sampleB"),
    ("2-200-209-AAG", "chr2", 200, 209, "AAG", "ATXN1", "ENSG00000124788", None,
     1, 5, "intron", "4x:sampleC"),
    ("X-300-310-CCG", "chrX", 300, 310, "CCG", "AFF2", "ENSG00000155966", "4",
     1, 1, "CDS", "7x:sampleD"),
]


def _build_test_db(path):
    """Create a minimal loci table plus the reference-region index."""
    conn = sqlite3.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT,
        Chrom TEXT,
        Start0Based INTEGER,
        End1Based INTEGER,
        CanonicalMotif TEXT,
        GeneTableGeneSymbol TEXT,
        gene_id TEXT,
        VariationClusterSizeDiff TEXT,
        NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER,
        gene_region_rank INTEGER,
        gene_region TEXT,
        OutlierSampleIds_AllAlleles TEXT
    )""")
    conn.executemany("INSERT INTO loci VALUES (" + ",".join("?" * 12) + ")", TEST_LOCI)
    conn.execute(f"CREATE INDEX {results_server.REFERENCE_REGION_INDEX} "
                 f"ON loci(Chrom, Start0Based, End1Based)")
    conn.commit()
    conn.close()


class ClassifySearchTermTests(unittest.TestCase):

    def test_locus_ids(self):
        # Stored locus ids spell the motif in upper case, so a lower-case spelling is normalized
        # rather than matched with a case-insensitive comparison that would give up the index.
        for term, expected in [("2-89831737-89831752-CCATT", "2-89831737-89831752-CCATT"),
                               ("1-100-110-AT", "1-100-110-AT"),
                               ("X-300-310-ccg", "X-300-310-CCG")]:
            kind, value, error = results_server.classify_search_term(term)
            self.assertIsNone(error, term)
            self.assertEqual((kind, value), ("locus_id", expected), term)

    def test_gene_ids(self):
        # Stored gene ids are upper case with no version suffix; both variants normalize to it,
        # so they match instead of being classified as gene ids that can never match anything.
        for term in ["ENSG00000102081", "ENSG00000102081.5", "ensg00000102081"]:
            kind, value, error = results_server.classify_search_term(term)
            self.assertIsNone(error, term)
            self.assertEqual((kind, value), ("gene_id", "ENSG00000102081"), term)

    def test_regions(self):
        for term, expected in [
            ("chr16:11579459-11579529", ("chr16", 11579459, 11579529)),
            ("16:100-200", ("chr16", 100, 200)),
            ("chr16:11579459", ("chr16", 11579459, 11579460)),
            ("chr16", ("chr16", 0, None)),
            ("chrX", ("chrX", 0, None)),
            ("X", ("chrX", 0, None)),
        ]:
            kind, value, error = results_server.classify_search_term(term)
            self.assertIsNone(error, term)
            self.assertEqual((kind, value), ("region", expected), term)

    def test_gene_symbols(self):
        # Anything that is not one of the three recognized shapes is a gene symbol. This is the
        # case that would break if regions were recognized by handing every term to
        # parse_reference_region, which reads a bare alphanumeric string as a whole chromosome.
        for term in ["FMR1", "HTT", "C9orf72", "ATXN1", "RFC1"]:
            kind, value, error = results_server.classify_search_term(term)
            self.assertIsNone(error, term)
            self.assertEqual((kind, value), ("gene_symbol", term), term)

    def test_gene_id_normalization_helper(self):
        # The legacy ?gene_id= parameter runs through the same normalization as a search term, so
        # both paths match the stored spelling rather than one silently returning nothing.
        self.assertEqual(results_server.normalize_gene_id("ensg00000102081"), "ENSG00000102081")
        self.assertEqual(results_server.normalize_gene_id("ENSG00000102081.5"), "ENSG00000102081")
        self.assertEqual(results_server.normalize_gene_id("  ENSG00000102081  "), "ENSG00000102081")
        # Not an Ensembl gene id: callers fall back to the value as typed.
        self.assertIsNone(results_server.normalize_gene_id("FMR1"))
        self.assertIsNone(results_server.normalize_gene_id(""))
        self.assertIsNone(results_server.normalize_gene_id(None))

    def test_unparseable_region_reports_an_error(self):
        # Shaped like a region, so it is not silently demoted to a gene symbol.
        kind, value, error = results_server.classify_search_term("chr1:200-100")
        self.assertIsNone(kind)
        self.assertIsNone(value)
        self.assertTrue(error)


class SearchQueryTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "test_results.db")
        _build_test_db(cls.db_path)
        results_server.app.config["DB_PATH"] = cls.db_path
        results_server.app.config["DB_COLUMNS_SET"] = {
            row[1] for row in sqlite3.connect(cls.db_path).execute(
                "SELECT * FROM pragma_table_info('loci')")}
        results_server.app.config["HAS_SKINNY"] = False
        results_server.app.config.pop("MAX_LOCUS_SPAN", None)

    def _query(self, **extra_params):
        """Run the count and select queries for a set of params against the test database."""
        params = dict({"outlier_type": "all", "page": 1, "page_size": 50}, **extra_params)
        select_query, count_query, _, _, sql_params, sql_params_with_pagination = \
            results_server.build_api_query(params)
        conn = sqlite3.connect(self.db_path)
        try:
            locus_ids = {row[0] for row in conn.execute(select_query, sql_params_with_pagination)}
            total = conn.execute(count_query, sql_params).fetchone()[0]
        finally:
            conn.close()
        self.assertEqual(total, len(locus_ids))
        return locus_ids

    def _search(self, text):
        terms = []
        for term in text.split(","):
            term = term.strip()
            if not term:
                continue
            kind, value, error = results_server.classify_search_term(term)
            self.assertIsNone(error, term)
            terms.append((kind, value))
        return self._query(search=text, search_terms=terms)

    def test_single_term_of_each_kind(self):
        self.assertEqual(self._search("FMR1"), {"1-100-110-AT"})
        self.assertEqual(self._search("ENSG00000197386"), {"1-500-600-CAG"})
        self.assertEqual(self._search("2-200-209-AAG"), {"2-200-209-AAG"})
        self.assertEqual(self._search("chr1:105-106"), {"1-100-110-AT"})
        self.assertEqual(self._search("chrX"), {"X-300-310-CCG"})

    def test_terms_of_different_kinds_are_ored(self):
        # The whole point of the merged box: a locus matching any one term is kept, where the
        # four separate filters would have required a locus to satisfy all of them at once.
        self.assertEqual(
            self._search("FMR1, ENSG00000197386"),
            {"1-100-110-AT", "1-500-600-CAG"})
        self.assertEqual(
            self._search("FMR1, chr2:200-201, X-300-310-CCG"),
            {"1-100-110-AT", "2-200-209-AAG", "X-300-310-CCG"})
        # Terms that would have ANDed to nothing under the old filters.
        self.assertEqual(
            self._search("HTT, chr2"),
            {"1-500-600-CAG", "2-200-209-AAG"})

    def test_gene_symbol_matching_is_case_insensitive_and_partial(self):
        self.assertEqual(self._search("fmr1"), {"1-100-110-AT"})
        self.assertEqual(self._search("ATXN"), {"2-200-209-AAG"})

    def test_no_match(self):
        self.assertEqual(self._search("NOTAGENE"), set())

    def test_max_variation_cluster_size_diff(self):
        # Values are 13 (FMR1), 148 (HTT), 4 (AFF2) and NULL (ATXN1). A locus with no variation
        # cluster has no variation beyond the repeat itself, so it passes any maximum.
        self.assertEqual(
            self._query(max_variation_cluster_size_diff=20),
            {"1-100-110-AT", "X-300-310-CCG", "2-200-209-AAG"})
        self.assertEqual(
            self._query(max_variation_cluster_size_diff=4),
            {"X-300-310-CCG", "2-200-209-AAG"})
        self.assertEqual(
            self._query(max_variation_cluster_size_diff=1000),
            {locus[0] for locus in TEST_LOCI})
        # Compared numerically, not lexically: "148" must not sort below "20".
        self.assertNotIn("1-500-600-CAG", self._query(max_variation_cluster_size_diff=20))

    def test_variation_cluster_filter_combines_with_search(self):
        # Separate filters still AND with each other; only terms inside the search box OR.
        self.assertEqual(
            self._query(search="FMR1, HTT",
                        search_terms=[("gene_symbol", "FMR1"), ("gene_symbol", "HTT")],
                        max_variation_cluster_size_diff=20),
            {"1-100-110-AT"})

    def test_filters_degrade_when_the_database_lacks_the_column(self):
        # A database built without the gene table or without TRExplorer annotations must not
        # 500; the affected term drops out of the OR, and a filter left with nothing matches
        # nothing rather than everything.
        original = results_server.app.config["DB_COLUMNS_SET"]
        try:
            results_server.app.config["DB_COLUMNS_SET"] = original - {"GeneTableGeneSymbol"}
            self.assertEqual(self._search("FMR1"), set())
            self.assertEqual(self._search("FMR1, 2-200-209-AAG"), {"2-200-209-AAG"})
            results_server.app.config["DB_COLUMNS_SET"] = original - {"VariationClusterSizeDiff"}
            self.assertEqual(self._query(max_variation_cluster_size_diff=20), set())
        finally:
            results_server.app.config["DB_COLUMNS_SET"] = original


class SplitSearchTermsTests(unittest.TestCase):
    """Commas separate terms, except inside a region's coordinates."""

    def test_region_with_thousands_separators_stays_one_term(self):
        # parse_reference_region documents and accepts comma thousands separators, so splitting
        # the search string on every comma would shred a region copied out of a genome browser.
        self.assertEqual(results_server.split_search_terms("chr16:11,579,459-11,579,529"),
                         ["chr16:11,579,459-11,579,529"])
        self.assertEqual(results_server.split_search_terms("chr16:11,579,459"),
                         ["chr16:11,579,459"])
        self.assertEqual(
            results_server.split_search_terms("FMR1, chr16:11,579,459-11,579,529, ENSG00000102081"),
            ["FMR1", "chr16:11,579,459-11,579,529", "ENSG00000102081"])

    def test_ordinary_lists_still_split(self):
        self.assertEqual(results_server.split_search_terms("FMR1,HTT"), ["FMR1", "HTT"])
        self.assertEqual(results_server.split_search_terms("FMR1, HTT"), ["FMR1", "HTT"])
        # Bare chromosomes have no coordinate span, so they are not swallowed as one region.
        self.assertEqual(results_server.split_search_terms("1,2,3"), ["1", "2", "3"])
        self.assertEqual(results_server.split_search_terms("1, 2, 3"), ["1", "2", "3"])
        self.assertEqual(results_server.split_search_terms("1-100-110-AT, 2-200-209-AAG"),
                         ["1-100-110-AT", "2-200-209-AAG"])

    def test_whitespace_and_empty_terms(self):
        self.assertEqual(results_server.split_search_terms("  FMR1 ,  chrX  "), ["FMR1", "chrX"])
        self.assertEqual(results_server.split_search_terms(",,FMR1,,"), ["FMR1"])
        self.assertEqual(results_server.split_search_terms(""), [])
        self.assertEqual(results_server.split_search_terms("   "), [])
        self.assertEqual(results_server.split_search_terms(None), [])


class RequiresFullLociTableTests(unittest.TestCase):
    """Which filters can be answered from the narrow sk_{ot} tables."""

    def setUp(self):
        self._saved = results_server.app.config.get("SKINNY_COLUMNS")
        results_server.app.config["SKINNY_COLUMNS"] = {}

    def tearDown(self):
        if self._saved is None:
            results_server.app.config.pop("SKINNY_COLUMNS", None)
        else:
            results_server.app.config["SKINNY_COLUMNS"] = self._saved

    def test_reference_regions_always_fall_back(self):
        # Start0Based / End1Based are never projected into a skinny table.
        self.assertTrue(results_server.requires_full_loci_table(
            {"reference_region_parsed": ("chr1", 0, None)}, "AllAlleles"))
        self.assertTrue(results_server.requires_full_loci_table(
            {"search_terms": [("gene_symbol", "FMR1"), ("region", ("chr1", 0, None))]}, "AllAlleles"))

    def test_variation_cluster_filter_follows_the_skinny_projection(self):
        # Skinny tables built before the column was projected cannot answer the filter.
        self.assertTrue(results_server.requires_full_loci_table(
            {"max_variation_cluster_size_diff": 20}, "AllAlleles"))
        # Once the database is rebuilt with the column, the narrow table answers it instead.
        results_server.app.config["SKINNY_COLUMNS"] = {
            "AllAlleles": {"LocusId", "VariationClusterSizeDiff"}}
        self.assertFalse(results_server.requires_full_loci_table(
            {"max_variation_cluster_size_diff": 20}, "AllAlleles"))
        # ... and only for the outlier types whose skinny table actually has it.
        self.assertTrue(results_server.requires_full_loci_table(
            {"max_variation_cluster_size_diff": 20}, "ShortAlleles"))

    def test_filters_the_skinny_tables_can_answer(self):
        self.assertFalse(results_server.requires_full_loci_table({}, "AllAlleles"))
        self.assertFalse(results_server.requires_full_loci_table(
            {"search_terms": [("gene_symbol", "FMR1"), ("locus_id", "1-100-110-AT"),
                              ("gene_id", "ENSG00000102081")]}, "AllAlleles"))


if __name__ == "__main__":
    unittest.main()
