"""Unit tests for the merged Search box and the variation cluster filters.

The Search box replaced four separate filters (Gene Symbol(s), Gene ID, Locus Id(s) and
Reference Region), so the two things worth pinning down are:

  * classify_search_term(), which decides what each comma-separated term is from its shape, and
  * the WHERE clause build_api_query() produces, where terms OR together rather than AND.

A small throwaway loci table stands in for a real results database so the generated SQL
actually runs, including the numeric comparison against VariationClusterSizeDiff. That table
declares the column VARCHAR; the real build declares it BIGINT whenever every value it saw was
an integer, so VariationClusterTypeTests below covers that case through the real writer.
"""

import json
import os
import tempfile
import unittest

import duckdb_compat
import result_database
import results_server
import swim_plot


# LocusId, Chrom, Start0Based, End1Based, CanonicalMotif, gene symbol, Ensembl id,
# VariationClusterSizeDiff (VARCHAR here; see VariationClusterTypeTests), VariationClusterFilterReason, plus
# the columns every query touches. The two variation cluster fields are never both set on one
# locus: a size diff means a cluster was kept, a reason means one was computed and thrown away.
TEST_LOCI = [
    ("1-100-110-AT", "chr1", 100, 110, "AT", "FMR1", "ENSG00000102081", "13", None,
     1, 1, "CDS", "6x:sampleA"),
    ("1-500-600-CAG", "chr1", 500, 600, "AGC", "HTT", "ENSG00000197386", "148", None,
     1, 1, "CDS", "9x:sampleB"),
    ("2-200-209-AAG", "chr2", 200, 209, "AAG", "ATXN1", "ENSG00000124788", None, "DEPTH",
     1, 5, "intron", "4x:sampleC"),
    ("X-300-310-CCG", "chrX", 300, 310, "CCG", "AFF2", "ENSG00000155966", "4", None,
     1, 1, "CDS", "7x:sampleD"),
    ("3-700-712-AGGG", "chr3", 700, 712, "AGGG", "PABPN1", "ENSG00000100836", None, "EXTENSION",
     1, 1, "CDS", "5x:sampleE"),
]


def _build_test_db(path):
    """Create a minimal loci table."""
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT,
        Chrom TEXT,
        Start0Based INTEGER,
        End1Based INTEGER,
        CanonicalMotif TEXT,
        GeneTableGeneSymbol TEXT,
        gene_id TEXT,
        VariationClusterSizeDiff TEXT,
        VariationClusterFilterReason TEXT,
        NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER,
        gene_region_rank INTEGER,
        gene_region TEXT,
        OutlierSampleIds_AllAlleles TEXT
    )""")
    conn.executemany("INSERT INTO loci VALUES (" + ",".join("?" * 13) + ")", TEST_LOCI)
    conn.commit()
    conn.close()


class ClassifySearchTermTests(unittest.TestCase):

    def test_locus_ids(self):
        # The build stores the locus id exactly as the input matrix spelled it, so the term is
        # kept as typed (upper-casing the motif here would build an id the database may not
        # contain) and the endpoints compare it case-insensitively.
        for term in ["2-89831737-89831752-CCATT", "1-100-110-AT", "X-300-310-ccg",
                     "chr1_alt-1-10-CAG"]:
            kind, value, error = results_server.classify_search_term(term)
            self.assertIsNone(error, term)
            self.assertEqual((kind, value), ("locus_id", term), term)

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
        cls.db_path = os.path.join(cls.tmpdir, "test_results.duckdb")
        _build_test_db(cls.db_path)
        results_server.app.config["DB_PATH"] = cls.db_path
        conn = duckdb_compat.connect(cls.db_path, read_only=True)
        try:
            results_server.app.config["DB_COLUMNS_SET"] = duckdb_compat.table_columns(conn, "loci")
        finally:
            conn.close()

    def _query(self, **extra_params):
        """Run the count and select queries for a set of params against the test database."""
        params = dict({"outlier_type": "all", "page": 1, "page_size": 50}, **extra_params)
        select_query, count_query, _, _, sql_params, sql_params_with_pagination = \
            results_server.build_api_query(params)
        conn = duckdb_compat.connect(self.db_path, read_only=True)
        try:
            locus_ids = {row[0] for row in conn.execute(select_query, sql_params_with_pagination)}
            total = conn.execute(count_query, sql_params).fetchone()[0]
        finally:
            conn.close()
        self.assertEqual(total, len(locus_ids))
        return locus_ids

    def _search(self, text):
        # split_search_terms is what validate_params runs on the raw search string, so the term
        # list here is the one the server would really build. Splitting on every comma instead
        # would shred a comma-formatted region and test a case the server never produces.
        terms = []
        for term in results_server.split_search_terms(text):
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
        # Values are 13 (FMR1), 148 (HTT), 4 (AFF2) and NULL (ATXN1, PABPN1). A locus with no
        # variation cluster has no variation beyond the repeat itself, so it passes any maximum.
        self.assertEqual(
            self._query(max_variation_cluster_size_diff=20),
            {"1-100-110-AT", "X-300-310-CCG", "2-200-209-AAG", "3-700-712-AGGG"})
        self.assertEqual(
            self._query(max_variation_cluster_size_diff=4),
            {"X-300-310-CCG", "2-200-209-AAG", "3-700-712-AGGG"})
        self.assertEqual(
            self._query(max_variation_cluster_size_diff=1000),
            {locus[0] for locus in TEST_LOCI})
        # Compared numerically, not lexically: "148" must not sort below "20".
        self.assertNotIn("1-500-600-CAG", self._query(max_variation_cluster_size_diff=20))

    def test_exclude_variation_cluster_filter_reasons(self):
        # ATXN1 is DEPTH (the gold icon in the VC column), PABPN1 is EXTENSION (dark red). Each
        # box drops only its own reason, and a locus that was never filtered has no reason
        # recorded, so it survives both.
        never_filtered = {"1-100-110-AT", "1-500-600-CAG", "X-300-310-CCG"}
        self.assertEqual(
            self._query(exclude_vc_depth_filtered=True), never_filtered | {"3-700-712-AGGG"})
        self.assertEqual(
            self._query(exclude_vc_size_filtered=True), never_filtered | {"2-200-209-AAG"})
        self.assertEqual(
            self._query(exclude_vc_depth_filtered=True, exclude_vc_size_filtered=True),
            never_filtered)

    def test_exclude_variation_cluster_filter_reasons_with_the_size_filter(self):
        # The two kinds of variation cluster filter apply to disjoint sets of loci (a locus has
        # either a size diff or a reason, never both), so combining them keeps whatever passes both.
        self.assertEqual(
            self._query(exclude_vc_depth_filtered=True, max_variation_cluster_size_diff=20),
            {"1-100-110-AT", "X-300-310-CCG", "3-700-712-AGGG"})

    def test_exclude_variation_cluster_filter_reasons_exclude_nothing_without_the_column(self):
        # These are exclusion flags, so the same rule as the maximum below applies: a database
        # with no VariationClusterFilterReason column never recorded any locus as filtered out,
        # so there is nothing for them to exclude and every locus survives. Failing closed here
        # instead would empty the whole result the moment either box was ticked.
        all_loci = {locus[0] for locus in TEST_LOCI}
        original = results_server.app.config["DB_COLUMNS_SET"]
        try:
            results_server.app.config["DB_COLUMNS_SET"] = original - {"VariationClusterFilterReason"}
            self.assertEqual(self._query(exclude_vc_depth_filtered=True), all_loci)
            self.assertEqual(self._query(exclude_vc_size_filtered=True), all_loci)
            self.assertEqual(
                self._query(exclude_vc_depth_filtered=True, exclude_vc_size_filtered=True),
                all_loci)
        finally:
            results_server.app.config["DB_COLUMNS_SET"] = original

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
        finally:
            results_server.app.config["DB_COLUMNS_SET"] = original

    def test_variation_cluster_filter_passes_everything_without_the_column(self):
        # This filter is the exception to the rule above: its clause deliberately keeps loci whose
        # VariationClusterSizeDiff is NULL or "", so a database with no such column at all is the
        # same situation for every locus and every locus passes. Matching nothing instead would
        # contradict how the filter treats missing data when the column is present.
        original = results_server.app.config["DB_COLUMNS_SET"]
        try:
            results_server.app.config["DB_COLUMNS_SET"] = original - {"VariationClusterSizeDiff"}
            self.assertEqual(self._query(max_variation_cluster_size_diff=20),
                             self._query())
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


class VariationClusterSortTests(unittest.TestCase):
    """The ORDER BY the sort keys build: numeric size diffs, and a decided order for ties."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "sort_test_results.duckdb")
        conn = duckdb_compat.connect(cls.db_path)
        conn.execute("""CREATE TABLE loci (
            LocusId TEXT,
            VariationClusterSizeDiff TEXT,
            NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER,
            OutlierSampleIds_AllAlleles TEXT
        )""")
        # Stored as TEXT, the way TRExplorer writes them. Compared as strings "9" beats "20"
        # beats "148", which is the exact reverse of the intended widest-cluster-first order.
        conn.executemany("INSERT INTO loci VALUES (?, ?, ?, ?)", [
            ("1-100-110-AT", "9", 1, "6x:sampleA"),
            ("2-200-209-AAG", "20", 1, "6x:sampleB"),
            ("3-700-712-AGGG", "148", 1, "6x:sampleC"),
            ("X-300-310-CCG", None, 1, "6x:sampleD"),
        ])
        conn.commit()
        conn.close()

    @classmethod
    def tearDownClass(cls):
        for name in os.listdir(cls.tmpdir):
            os.remove(os.path.join(cls.tmpdir, name))
        os.rmdir(cls.tmpdir)

    def setUp(self):
        self._saved = {key: results_server.app.config.get(key)
                       for key in ("DB_PATH", "DB_COLUMNS_SET")}
        results_server.app.config["DB_PATH"] = self.db_path
        results_server.app.config["DB_COLUMNS_SET"] = {
            "LocusId", "VariationClusterSizeDiff",
            "NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles",
            "OutlierSampleIds_AllAlleles"}

    def tearDown(self):
        for key, value in self._saved.items():
            if value is None:
                results_server.app.config.pop(key, None)
            else:
                results_server.app.config[key] = value

    def _ordered_locus_ids(self, **extra_params):
        """Run the select query for a set of params and return the LocusIds in row order."""
        params = dict({"outlier_type": "all", "page": 1, "page_size": 50}, **extra_params)
        select_query, _, _, _, _, sql_params_with_pagination = \
            results_server.build_api_query(params)
        conn = duckdb_compat.connect(self.db_path, read_only=True)
        try:
            return [row[0] for row in conn.execute(select_query, sql_params_with_pagination)]
        finally:
            conn.close()

    def test_sorts_numerically_not_lexically(self):
        # 148 > 20 > 9, and the locus with no variation cluster at all sorts last.
        self.assertEqual(
            self._ordered_locus_ids(sort_by=["variation_cluster"]),
            ["3-700-712-AGGG", "2-200-209-AAG", "1-100-110-AT", "X-300-310-CCG"])

    def test_order_by_casts_the_text_column(self):
        # The cast is what makes the ordering above numeric, so pin the expression itself too.
        # TRY_CAST rather than CAST: DuckDB raises on a value that will not parse.
        order_by = results_server.build_api_order_by(
            {"outlier_type": "all", "sort_by": ["variation_cluster"]},
            results_server.app.config["DB_COLUMNS_SET"])
        self.assertIn("TRY_CAST(VariationClusterSizeDiff AS BIGINT) DESC", order_by)

    def test_ties_are_broken_by_locus_id(self):
        # DuckDB scans a table with several threads and returns rows in whatever order they
        # finish, so rows the sort keys leave tied come back in an arbitrary order: paging
        # through such a result can show one locus twice and another not at all. Every sort
        # ends with LocusId so no tie is left undecided.
        order_by = results_server.build_api_order_by(
            {"outlier_type": "all", "sort_by": ["count"]},
            results_server.app.config["DB_COLUMNS_SET"])
        self.assertTrue(order_by.endswith(", LocusId ASC"), order_by)
        # All four loci tie on the count, and this database has none of the other sort columns,
        # so LocusId alone decides the order.
        self.assertEqual(
            self._ordered_locus_ids(sort_by=["count"]),
            ["1-100-110-AT", "2-200-209-AAG", "3-700-712-AGGG", "X-300-310-CCG"])

    def test_sort_is_skipped_when_the_database_lacks_the_column(self):
        # An older database has no such column; the sort key drops out rather than producing
        # SQL that names a column the database does not have.
        results_server.app.config["DB_COLUMNS_SET"] = (
            results_server.app.config["DB_COLUMNS_SET"] - {"VariationClusterSizeDiff"})
        order_by = results_server.build_api_order_by(
            {"outlier_type": "all", "sort_by": ["variation_cluster"]},
            results_server.app.config["DB_COLUMNS_SET"])
        self.assertNotIn("VariationClusterSizeDiff", order_by)
        self.assertEqual(
            sorted(self._ordered_locus_ids(sort_by=["variation_cluster"])),
            ["1-100-110-AT", "2-200-209-AAG", "3-700-712-AGGG", "X-300-310-CCG"])


class VariationClusterTypeTests(unittest.TestCase):
    """The size-diff filter has to work on the column type the real build actually writes.

    result_database.infer_column_types declares VariationClusterSizeDiff BIGINT when every value
    the build saw was an integer, which is the ordinary case, and VARCHAR as soon as one is not.
    The fixture above hand-writes the column VARCHAR, so it never exercised the BIGINT case,
    where comparing the column to '' raised a ConversionException and 500'd every search using
    the filter. These tests go through result_database.write_loci_table so the declared type is
    the one the real build chooses, and run the query build_api_query actually emits.
    """

    COLUMNS = ["LocusId", "Chrom", "Start0Based", "End1Based", "CanonicalMotif",
               "GeneTableGeneSymbol", "gene_id", "VariationClusterSizeDiff",
               "VariationClusterFilterReason",
               "NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles", "gene_region_rank",
               "gene_region", "OutlierSampleIds_AllAlleles"]

    def _records(self, extra=None):
        """TEST_LOCI as record dicts, with the size diffs as ints so the column lands BIGINT."""
        records = []
        for locus in TEST_LOCI:
            record = dict(zip(self.COLUMNS, locus))
            size_diff = record["VariationClusterSizeDiff"]
            record["VariationClusterSizeDiff"] = None if size_diff is None else int(size_diff)
            records.append(record)
        if extra:
            records.append(dict(zip(self.COLUMNS, extra)))
        return records

    def _query(self, records, **extra_params):
        """Write records through the real writer, then run build_api_query's own select query."""
        db_path = os.path.join(tempfile.mkdtemp(), "typed.duckdb")
        conn = duckdb_compat.connect(db_path)
        try:
            result_database.write_loci_table(conn, records, self.COLUMNS)
            conn.commit()
            declared = conn.execute(
                "SELECT data_type FROM duckdb_columns() WHERE table_name = 'loci' "
                "AND column_name = 'VariationClusterSizeDiff'").fetchone()[0]
        finally:
            conn.close()

        saved = {key: results_server.app.config.get(key) for key in ("DB_PATH", "DB_COLUMNS_SET")}
        results_server.app.config["DB_PATH"] = db_path
        results_server.app.config["DB_COLUMNS_SET"] = set(self.COLUMNS)
        try:
            params = dict({"outlier_type": "all", "page": 1, "page_size": 50}, **extra_params)
            select_query, _c, _m, _g, _p, sql_params_with_pagination = \
                results_server.build_api_query(params)
            conn = duckdb_compat.connect(db_path, read_only=True)
            try:
                return declared, {row[0] for row in
                                  conn.execute(select_query, sql_params_with_pagination)}
            finally:
                conn.close()
        finally:
            for key, value in saved.items():
                if value is None:
                    results_server.app.config.pop(key, None)
                else:
                    results_server.app.config[key] = value

    def test_all_integer_column_is_written_bigint_and_still_filters(self):
        declared, kept = self._query(self._records(), max_variation_cluster_size_diff=20)
        self.assertEqual(declared, "BIGINT")
        # 148 (HTT) is over the maximum. The two loci with no cluster pass any maximum.
        self.assertEqual(
            kept, {"1-100-110-AT", "X-300-310-CCG", "2-200-209-AAG", "3-700-712-AGGG"})

    def test_mixed_column_is_written_varchar_and_still_filters(self):
        extra = ("4-1-2-AC", "chr4", 1, 2, "AC", "XYZ", "ENSG00000000001", "unknown", None,
                 1, 1, "CDS", "3x:sampleF")
        declared, kept = self._query(self._records(extra), max_variation_cluster_size_diff=20)
        self.assertEqual(declared, "VARCHAR")
        # A value that will not parse means the cluster size is unknown, so the locus is kept.
        self.assertEqual(
            kept, {"1-100-110-AT", "X-300-310-CCG", "2-200-209-AAG", "3-700-712-AGGG", "4-1-2-AC"})

    def test_sort_by_variation_cluster_on_a_bigint_column(self):
        # The ORDER BY casts the column; the cast must also be valid when it is already BIGINT.
        _declared, kept = self._query(self._records(), sort_by=["variation_cluster"])
        self.assertEqual(len(kept), len(TEST_LOCI))


# LocusId, Chrom, Start0Based, End1Based, Motif, CanonicalMotif, MotifSize, gene symbol,
# Ensembl id, gene region, its rank, the outlier count, the first affected phenotype and the
# outlier sample ids. Deliberately without a VariationClusterSizeDiff column: that is what makes
# this fixture the "no TRExplorer annotations" database the endpoint tests below need.
ENDPOINT_LOCI = [
    ("chr1-100-110-AT", "chr1", 100, 110, "AT", "AT", 2, "FMR1", "ENSG00000102081", "CDS", 1,
     2, "seizures", "6x:sampleA,5x:sampleB"),
    ("chr2-200-209-AAG", "chr2", 200, 209, "AAG", "AAG", 3, "ATXN1", "ENSG00000124788", "intron", 5,
     1, "ataxia", "4x:sampleC"),
]


def _build_endpoint_test_db(path):
    """Create a loci + swim_plot database with no variation-cluster columns."""
    conn = duckdb_compat.connect(path)
    conn.execute("""CREATE TABLE loci (
        LocusId TEXT,
        Chrom TEXT,
        Start0Based INTEGER,
        End1Based INTEGER,
        Motif TEXT,
        CanonicalMotif TEXT,
        MotifSize INTEGER,
        GeneTableGeneSymbol TEXT,
        gene_id TEXT,
        gene_region TEXT,
        gene_region_rank INTEGER,
        NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER,
        FirstAffectedPhenotype_AllAlleles TEXT,
        OutlierSampleIds_AllAlleles TEXT
    )""")
    conn.executemany("INSERT INTO loci VALUES (" + ",".join("?" * 14) + ")", ENDPOINT_LOCI)

    def _swim_row(**values):
        """Return a full-width swim_plot row dict, defaulting every unnamed column to NULL."""
        return {column: values.get(column) for column in swim_plot.SWIM_PLOT_COLUMNS}

    result_database.write_swim_plot(
        conn,
        [
            # The phenotypes match the loci table's FirstAffectedPhenotype_AllAlleles column: the
            # build fills both from the same sample metadata, and the phenotype_keyword filter is
            # answered from these rows.
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="2bp",
                      allele_size=6, sample_id="sampleA", affected_status="Affected",
                      phenotype_description="seizures",
                      LocusId="chr1-100-110-AT", Motif="AT", CanonicalMotif="AT", MotifSize=2),
            _swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="3bp",
                      allele_size=4, sample_id="sampleC", affected_status="Affected",
                      phenotype_description="ataxia",
                      LocusId="chr2-200-209-AAG", Motif="AAG", CanonicalMotif="AAG", MotifSize=3),
        ],
        columns=swim_plot.SWIM_PLOT_COLUMNS,
    )
    conn.commit()
    conn.close()


class EndpointConfigTestCase(unittest.TestCase):
    """Base class for the tests that go through configure_app and the Flask test client.

    configure_app rewrites the whole of app.config, which the query-building tests above share,
    so the config is snapshotted and put back afterwards rather than left pointing at this
    fixture.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "endpoint_results.duckdb")
        _build_endpoint_test_db(cls.db_path)
        cls._saved_config = dict(results_server.app.config)
        results_server.configure_app(
            cls.db_path, annotations_db=os.path.join(cls.tmpdir, "annotations.duckdb"))
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        results_server.app.config.clear()
        results_server.app.config.update(cls._saved_config)


class BlankFilterValueTests(EndpointConfigTestCase):
    """A filter whose value strips down to nothing must not become an empty SQL list.

    Each of these parameters is a comma-separated list. A value like "," is truthy as a request
    string but parses to no entries, and emitting "LocusId IN ()" or an empty "()" for it is a
    DuckDB syntax error, so the request 500'd. An empty list carries no information, so the
    filter is dropped instead and the search comes back unfiltered.
    """

    ALL_LOCUS_IDS = {locus[0] for locus in ENDPOINT_LOCI}

    def _loci(self, query):
        response = self.client.get(f"/api/v1/loci?outlier_type=all&{query}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return response.get_json()

    def test_blank_list_filters_are_dropped(self):
        for param in ("locus_id", "gene_symbol", "phenotype_keyword", "sample_id_keyword",
                      "sample_id_like"):
            for raw in (",", " , ", ",,"):
                data = self._loci(f"{param}={raw}")
                self.assertEqual({r["LocusId"] for r in data["results"]}, self.ALL_LOCUS_IDS,
                                 f"{param}={raw}")
                self.assertEqual(data["total"], len(ENDPOINT_LOCI), f"{param}={raw}")
                # Nothing was filtered, so nothing is reported as filtered either.
                self.assertNotIn(param, data["filters_applied"], f"{param}={raw}")

    def test_blank_list_filters_do_not_break_the_export(self):
        # /api/v1/export runs the same validate_params and build_api_query pair.
        response = self.client.get("/api/v1/export?outlier_type=all&format=bed&locus_id=,")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertIn("chr1-100-110-AT", response.get_data(as_text=True))

    def test_non_blank_list_filters_still_filter(self):
        # The guard drops only the empty case; a real value filters as before.
        self.assertEqual(
            {r["LocusId"] for r in self._loci("locus_id=chr1-100-110-AT")["results"]},
            {"chr1-100-110-AT"})
        self.assertEqual(
            {r["LocusId"] for r in self._loci("gene_symbol=ATXN1")["results"]},
            {"chr2-200-209-AAG"})
        self.assertEqual(
            {r["LocusId"] for r in self._loci("phenotype_keyword=ataxia")["results"]},
            {"chr2-200-209-AAG"})
        self.assertEqual(
            {r["LocusId"] for r in self._loci("sample_id_keyword=sampleA")["results"]},
            {"chr1-100-110-AT"})
        self.assertEqual(
            {r["LocusId"] for r in self._loci("sample_id_like=sampleC")["results"]},
            {"chr2-200-209-AAG"})


class VariationClusterFilterWithoutTheColumnTests(EndpointConfigTestCase):
    """The swim plot and the results table must answer this filter the same way.

    The clause keeps a locus whose variation cluster size is unknown, so a database with no
    VariationClusterSizeDiff column at all is that same situation for every locus and the filter
    adds no clause. The swim-plot endpoint used to append "1=0" instead, which emptied the chart
    for a request whose results table returned every locus.
    """

    def _swim_sample_ids(self, query):
        response = self.client.get(f"/api/v1/swim_plot_data?outlier_type=all&{query}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {entry["sample_id"] for entry in response.get_json()["data"]}

    def test_swim_plot_keeps_every_locus(self):
        unfiltered = self._swim_sample_ids("")
        self.assertEqual(unfiltered, {"sampleA", "sampleC"})
        self.assertEqual(self._swim_sample_ids("max_variation_cluster_size_diff=20"), unfiltered)

    def test_swim_plot_and_results_table_agree(self):
        response = self.client.get(
            "/api/v1/loci?outlier_type=all&max_variation_cluster_size_diff=20")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual({r["LocusId"] for r in response.get_json()["results"]},
                         {locus[0] for locus in ENDPOINT_LOCI})
        self.assertEqual(self._swim_sample_ids("max_variation_cluster_size_diff=20"),
                         {"sampleA", "sampleC"})

    def test_an_unparseable_maximum_is_still_rejected(self):
        # Dropping the "1=0" branch must not drop the 400 for a value that is not an integer.
        response = self.client.get(
            "/api/v1/swim_plot_data?outlier_type=all&max_variation_cluster_size_diff=abc")
        self.assertEqual(response.status_code, 400)


class DiseaseCatalogStartupTests(unittest.TestCase):
    """Either disease catalog is enough on its own.

    configure_app used to load both only when the Broad variant catalog was supplied, so a
    STRchive-only server silently had no disease details at all.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "catalog_results.duckdb")
        _build_endpoint_test_db(cls.db_path)
        cls.strchive_path = os.path.join(cls.tmpdir, "STRchive-loci.json")
        with open(cls.strchive_path, "wt") as output_file:
            json.dump([{
                "locus_id": "FMR1_test",
                "disease": "Test syndrome",
                "chrom": "chr1",
                "start_hg38": 101,
                "stop_hg38": 110,
                "reference_motif_reference_orientation": ["AT"],
                "pathogenic_motif_reference_orientation": ["AT"],
            }], output_file)

    def setUp(self):
        self._saved_config = dict(results_server.app.config)

    def tearDown(self):
        results_server.app.config.clear()
        results_server.app.config.update(self._saved_config)

    def _configure(self, **kwargs):
        results_server.configure_app(
            self.db_path, annotations_db=os.path.join(self.tmpdir, "annotations.duckdb"), **kwargs)
        return results_server.app.config["LOOKUPS"]

    def test_strchive_only_startup_loads_the_fallback(self):
        lookups = self._configure(strchive_loci_json=self.strchive_path)
        self.assertTrue(lookups["strchive_trees"])
        # The catalog trees stay empty: no Broad variant catalog was supplied.
        self.assertEqual(lookups["disease_trees"], {})
        # And the serving path actually resolves a locus through it, which is the point of
        # loading it at all.
        disease_info = results_server.compute_known_disease_info(
            {"LocusId": "chr1-100-110-AT", "Chrom": "chr1", "Start0Based": 100,
             "End1Based": 110, "Motif": "AT"},
            lookups)
        self.assertEqual(disease_info["source"], "STRchive")
        self.assertEqual(disease_info["locus_id"], "FMR1_test")

    def test_no_catalog_at_all_leaves_the_lookups_empty(self):
        lookups = self._configure()
        self.assertEqual(lookups["strchive_trees"], {})
        self.assertEqual(lookups["disease_trees"], {})
        self.assertEqual(lookups["locus"], {})


# LocusId, Chrom, Start0Based, End1Based, Motif, CanonicalMotif, MotifSize, gene symbol, Ensembl
# id, gene region, its rank, the outlier count and the outlier sample ids. The first locus spells
# its motif in lower case, the way a TRGT catalog may write it and the way build_database stores
# it; the second carries an underscore in its contig name, and the third is the same id with the
# underscore replaced, which is what an unescaped ILIKE wildcard would also match. The fourth sits
# at coordinates a genome browser writes with thousands separators.
SEARCH_LOCI_COLUMNS = ["LocusId", "Chrom", "Start0Based", "End1Based", "Motif", "CanonicalMotif",
                       "MotifSize", "GeneTableGeneSymbol", "gene_id", "gene_region",
                       "gene_region_rank",
                       "NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles",
                       "OutlierSampleIds_AllAlleles"]
SEARCH_LOCI = [
    ("chr1-1-10-cag", "chr1", 1, 10, "cag", "AGC", 3, "GENEA", "ENSG00000000001", "CDS", 1, 1,
     "6x:sampleA"),
    ("chr1_alt-1-10-CAG", "chr1_alt", 1, 10, "CAG", "AGC", 3, "GENEB", "ENSG00000000002", "CDS",
     1, 1, "6x:sampleB"),
    ("chr1Xalt-1-10-CAG", "chr1Xalt", 1, 10, "CAG", "AGC", 3, "GENEC", "ENSG00000000003", "CDS",
     1, 1, "6x:sampleC"),
    ("chr16-11579458-11579529-CCG", "chr16", 11579458, 11579529, "CCG", "CCG", 3, "GENED",
     "ENSG00000000004", "CDS", 1, 1, "7x:sampleD"),
]


def _build_search_endpoint_db(path):
    """Write SEARCH_LOCI and one swim_plot row per locus through the real writers."""
    conn = duckdb_compat.connect(path)
    try:
        result_database.write_loci_table(
            conn, [dict(zip(SEARCH_LOCI_COLUMNS, locus)) for locus in SEARCH_LOCI],
            SEARCH_LOCI_COLUMNS)

        def _swim_row(**values):
            """Return a full-width swim_plot row dict, defaulting every unnamed column to NULL."""
            return {column: values.get(column) for column in swim_plot.SWIM_PLOT_COLUMNS}

        result_database.write_swim_plot(
            conn,
            [_swim_row(outlier_type="AllAlleles", outlier_rank=1, motif_category="3bp",
                       allele_size=6 + index, sample_id=f"sample{letter}",
                       affected_status="Affected", LocusId=locus[0], Motif=locus[4],
                       CanonicalMotif=locus[5], MotifSize=locus[6])
             for index, (letter, locus) in enumerate(zip("ABCD", SEARCH_LOCI))],
            columns=swim_plot.SWIM_PLOT_COLUMNS,
        )
        conn.commit()
    finally:
        conn.close()


class SearchEndToEndTests(unittest.TestCase):
    """The ?search= parameter driven through the Flask client, on a database the writers wrote.

    The SQL-level tests above call build_api_query directly, so nothing exercised
    validate_params -> split_search_terms -> classify_search_term -> build_api_query as one path.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.db_path = os.path.join(cls.tmpdir, "search_results.duckdb")
        _build_search_endpoint_db(cls.db_path)
        cls._saved_config = dict(results_server.app.config)
        results_server.configure_app(
            cls.db_path, annotations_db=os.path.join(cls.tmpdir, "annotations.duckdb"))
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    @classmethod
    def tearDownClass(cls):
        results_server.app.config.clear()
        results_server.app.config.update(cls._saved_config)

    def _loci(self, search):
        response = self.client.get("/api/v1/loci", query_string={"outlier_type": "all",
                                                                "search": search})
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        data = response.get_json()
        # The raw string is what the UI shows as the applied filter; the classified terms stay
        # a SQL internal.
        self.assertEqual(data["filters_applied"].get("search"), search)
        self.assertNotIn("search_terms", data["filters_applied"])
        return {r["LocusId"] for r in data["results"]}

    def _swim_locus_ids(self, search):
        response = self.client.get("/api/v1/swim_plot_data",
                                   query_string={"outlier_type": "all", "search": search})
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {entry["LocusId"] for entry in response.get_json()["data"]}

    def test_a_lowercase_locus_id_is_found(self):
        # The stored id spells its motif in lower case, which is what the input matrix wrote.
        # Upper-casing the motif before comparing made this search return nothing at all.
        self.assertEqual(self._loci("chr1-1-10-cag"), {"chr1-1-10-cag"})
        self.assertEqual(self._swim_locus_ids("chr1-1-10-cag"), {"chr1-1-10-cag"})

    def test_the_locus_id_comparison_is_case_insensitive_both_ways(self):
        for spelling in ("CHR1-1-10-CAG", "chr1-1-10-CAG", "Chr1-1-10-Cag"):
            self.assertEqual(self._loci(spelling), {"chr1-1-10-cag"}, spelling)
            self.assertEqual(self._swim_locus_ids(spelling), {"chr1-1-10-cag"}, spelling)

    def test_an_underscore_in_a_locus_id_matches_literally(self):
        # "_" is ILIKE's single-character wildcard, so without escaping this would also return
        # the chr1Xalt locus.
        self.assertEqual(self._loci("chr1_alt-1-10-CAG"), {"chr1_alt-1-10-CAG"})
        self.assertEqual(self._swim_locus_ids("chr1_alt-1-10-cag"), {"chr1_alt-1-10-CAG"})

    def test_a_comma_formatted_region_stays_one_term(self):
        # Splitting the search string on every comma would turn this into five terms
        # ("chr16:11", "579", "459-11", ...) and match nothing.
        self.assertEqual(self._loci("chr16:11,579,459-11,579,529"),
                         {"chr16-11579458-11579529-CCG"})
        self.assertEqual(self._swim_locus_ids("chr16:11,579,459-11,579,529"),
                         {"chr16-11579458-11579529-CCG"})

    def test_a_comma_formatted_region_alongside_other_terms(self):
        self.assertEqual(
            self._loci("GENEA, chr16:11,579,459-11,579,529"),
            {"chr1-1-10-cag", "chr16-11579458-11579529-CCG"})

    def test_an_unparseable_search_term_is_a_400(self):
        response = self.client.get("/api/v1/loci",
                                   query_string={"outlier_type": "all", "search": "chr1:200-100"})
        self.assertEqual(response.status_code, 400)


if __name__ == "__main__":
    unittest.main()
