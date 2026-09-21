"""Unit tests for duckdb_compat, the sqlite3-shaped wrapper TRails runs on.

Everything here is a behavior the rest of TRails depends on and that plain duckdb does not
give you: rows addressable by name, SQLite's NULL ordering, a description that survives a
later query on the same connection, and iteration that streams rather than materializing a
whole export.
"""

import os
import shutil
import tempfile
import unittest
import unittest.mock

import duckdb_compat


class RecordingDuckdbConnection:
    """Stands in for a duckdb connection and records what the wrapper passes down."""

    def __init__(self):
        self.executemany_calls = []
        self.closed = False

    def executemany(self, sql, seq_of_params):
        self.executemany_calls.append((sql, seq_of_params))

    def close(self):
        self.closed = True


class RowTests(unittest.TestCase):
    """Row has to behave like sqlite3.Row: index by position AND by column name."""

    def setUp(self):
        self.conn = duckdb_compat.connect(":memory:")
        self.conn.row_factory = duckdb_compat.Row
        self.conn.execute("CREATE TABLE loci (LocusId TEXT, Motif TEXT, MotifSize INTEGER)")
        self.conn.execute("INSERT INTO loci VALUES ('1-100-110-AT', 'AT', 2), ('2-1-9-CAG', 'CAG', 3)")

    def tearDown(self):
        self.conn.close()

    def test_index_by_position_and_name(self):
        row = self.conn.execute("SELECT * FROM loci ORDER BY LocusId").fetchone()
        self.assertEqual(row[0], "1-100-110-AT")
        self.assertEqual(row["LocusId"], "1-100-110-AT")
        self.assertEqual(row["MotifSize"], 2)

    def test_dict_and_keys(self):
        row = self.conn.execute("SELECT * FROM loci ORDER BY LocusId").fetchone()
        self.assertEqual(dict(row), {"LocusId": "1-100-110-AT", "Motif": "AT", "MotifSize": 2})
        self.assertEqual(row.keys(), ["LocusId", "Motif", "MotifSize"])
        # The server writes `col in row.keys()` to test for an optional column.
        self.assertIn("Motif", row.keys())
        self.assertNotIn("Nonexistent", row.keys())

    def test_iterating_a_row_yields_values(self):
        # sqlite3.Row iterates over the values, not the column names.
        row = self.conn.execute("SELECT * FROM loci ORDER BY LocusId").fetchone()
        self.assertEqual(list(row), ["1-100-110-AT", "AT", 2])
        self.assertEqual(len(row), 3)

    def test_in_tests_the_values_the_way_sqlite3_does(self):
        # sqlite3.Row has no __contains__, so `x in row` falls through to iteration over the
        # values. A Row that answered for its column names instead would quietly invert the
        # test for whoever writes `column in row` rather than `column in row.keys()`.
        row = self.conn.execute("SELECT * FROM loci ORDER BY LocusId").fetchone()
        self.assertIn("AT", row)
        self.assertNotIn("Motif", row)

    def test_unknown_column_raises_index_error(self):
        row = self.conn.execute("SELECT * FROM loci ORDER BY LocusId LIMIT 1").fetchone()
        with self.assertRaises(IndexError):
            row["NoSuchColumn"]

    def test_no_row_factory_gives_tuples(self):
        self.conn.row_factory = None
        self.assertEqual(self.conn.execute("SELECT MotifSize FROM loci ORDER BY LocusId").fetchall(),
                         [(2,), (3,)])


class NullOrderingTests(unittest.TestCase):
    """Every ORDER BY in TRails was written against SQLite's NULL placement.

    SQLite sorts NULLs first ascending and last descending. DuckDB's own default is NULLS
    LAST in both directions, which would silently move every unannotated locus in a sorted
    result page.
    """

    def setUp(self):
        self.conn = duckdb_compat.connect(":memory:")
        self.conn.execute("CREATE TABLE t (id INTEGER, score DOUBLE)")
        self.conn.execute("INSERT INTO t VALUES (1, 5.0), (2, NULL), (3, 9.0)")

    def tearDown(self):
        self.conn.close()

    def test_nulls_first_ascending(self):
        self.assertEqual([r[0] for r in self.conn.execute("SELECT id FROM t ORDER BY score, id").fetchall()],
                         [2, 1, 3])

    def test_nulls_last_descending(self):
        self.assertEqual(
            [r[0] for r in self.conn.execute("SELECT id FROM t ORDER BY score DESC, id").fetchall()],
            [3, 1, 2])

    def test_explicit_nulls_last_still_wins(self):
        # The sigma-rank sorts in the server's ORDER BY builder write "ASC NULLS LAST".
        self.assertEqual(
            [r[0] for r in self.conn.execute(
                "SELECT id FROM t ORDER BY score ASC NULLS LAST, id").fetchall()],
            [1, 3, 2])


class CursorLifetimeTests(unittest.TestCase):
    """A duckdb connection holds one result; sqlite3 hands out independent cursors."""

    def setUp(self):
        self.conn = duckdb_compat.connect(":memory:")
        self.conn.execute("CREATE TABLE loci (LocusId TEXT, Motif TEXT)")
        self.conn.execute("INSERT INTO loci VALUES ('a', 'AT')")

    def tearDown(self):
        self.conn.close()

    def test_description_survives_a_later_query(self):
        # validate_database reads the loci column list off this cursor after listing tables.
        cursor = self.conn.execute("SELECT * FROM loci ORDER BY LocusId LIMIT 1")
        duckdb_compat.list_tables(self.conn)
        self.assertEqual([d[0] for d in cursor.description], ["LocusId", "Motif"])

    def test_fetching_a_replaced_result_raises(self):
        # The alternative is silently returning the newer query's rows.
        cursor = self.conn.execute("SELECT LocusId FROM loci")
        self.conn.execute("SELECT Motif FROM loci")
        with self.assertRaises(RuntimeError):
            cursor.fetchall()


class IterationTests(unittest.TestCase):

    def setUp(self):
        self.conn = duckdb_compat.connect(":memory:")
        self.conn.execute("CREATE TABLE t AS SELECT i AS a FROM range(5000) AS r(i)")

    def tearDown(self):
        self.conn.close()

    def test_iteration_yields_every_row(self):
        # More rows than one fetch batch, since the export endpoints stream whole tables.
        self.assertEqual(sum(row[0] for row in self.conn.execute("SELECT a FROM t")),
                         sum(range(5000)))

    def test_iteration_streams_rather_than_materializing(self):
        # Stopping early must leave the rest of the result unread. Reading one more row off
        # the cursor afterwards shows how far iteration actually got: one batch in, not all
        # 5,000 rows, which is what fetching the whole result up front would have left.
        cursor = self.conn.execute("SELECT a FROM t ORDER BY a")
        self.assertEqual(next(iter(cursor))[0], 0)
        self.assertEqual(cursor.fetchone()[0], duckdb_compat._ITERATION_BATCH_ROWS)


class ConnectionApiTests(unittest.TestCase):
    """The rest of the sqlite3 connection surface TRails is written against."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.tmpdir, "t.duckdb")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_executemany(self):
        conn = duckdb_compat.connect(self.db_path)
        try:
            conn.execute("CREATE TABLE t (a INTEGER, b TEXT)")
            conn.executemany("INSERT INTO t VALUES (?, ?)", [(1, "x"), (2, "y")])
            self.assertEqual(conn.execute("SELECT COUNT(*) FROM t").fetchone()[0], 2)
        finally:
            conn.close()

    def test_insert_or_replace_against_a_primary_key(self):
        # Three writers rewrite a row that is already there, each against a declared
        # primary key: write_sample_count INSERT OR REPLACEs metadata.key, and in the
        # annotations database upsert_note writes notes.locus_id (ON CONFLICT DO UPDATE)
        # and add_tag writes the (locus_id, tag) key of tags (INSERT OR IGNORE). The
        # metadata table's shape stands in for all three here.
        conn = duckdb_compat.connect(self.db_path)
        try:
            conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value VARCHAR)")
            conn.executemany("INSERT INTO metadata VALUES (?, ?)",
                             [("num_samples", "10"), ("other", "x")])
            conn.execute("INSERT OR REPLACE INTO metadata VALUES (?, ?)", ["num_samples", "20"])
            self.assertEqual(conn.execute("SELECT key, value FROM metadata ORDER BY key").fetchall(),
                             [("num_samples", "20"), ("other", "x")])
        finally:
            conn.close()

    def test_read_only_connection_sees_a_finished_database(self):
        writer = duckdb_compat.connect(self.db_path)
        writer.execute("CREATE TABLE t AS SELECT 1 AS a")
        writer.close()
        reader = duckdb_compat.connect(self.db_path, read_only=True)
        try:
            self.assertEqual(reader.execute("SELECT a FROM t").fetchone()[0], 1)
            with self.assertRaises(duckdb_compat.OperationalError):
                reader.execute("CREATE TABLE other (a INTEGER)")
        finally:
            reader.close()

    def test_read_only_connections_share_one_instance_and_get_their_own_temp_schema(self):
        # get_db() opens one of these per request, and the variation cluster filter builds a
        # temp table on it, so two concurrent requests must not collide over that name.
        writer = duckdb_compat.connect(self.db_path)
        writer.execute("CREATE TABLE t AS SELECT 1 AS a")
        writer.close()
        first = duckdb_compat.connect(self.db_path, read_only=True)
        second = duckdb_compat.connect(self.db_path, read_only=True)
        try:
            first.execute("CREATE OR REPLACE TEMP TABLE vc_filtered_loci AS SELECT * FROM t")
            second.execute("CREATE OR REPLACE TEMP TABLE vc_filtered_loci AS SELECT * FROM t")
            self.assertEqual(first.execute("SELECT COUNT(*) FROM vc_filtered_loci").fetchone()[0], 1)
            self.assertEqual(second.execute("SELECT COUNT(*) FROM vc_filtered_loci").fetchone()[0], 1)
        finally:
            first.close()
            second.close()

    def test_table_helpers(self):
        conn = duckdb_compat.connect(self.db_path)
        try:
            conn.execute("CREATE TABLE loci (LocusId TEXT, Motif TEXT)")
            conn.execute("CREATE TABLE swim_plot (sample_id TEXT)")
            self.assertEqual(duckdb_compat.list_tables(conn), {"loci", "swim_plot"})
            self.assertEqual(duckdb_compat.table_columns(conn, "loci"), {"LocusId", "Motif"})
            self.assertEqual(duckdb_compat.table_columns(conn, "nonexistent"), set())
        finally:
            conn.close()


    def test_atomic_replace_of_a_closed_database(self):
        # The build writes to <db>.tmp and os.replace()s it into place. DuckDB drops the
        # .wal sidecar on close, so the single remaining file is the whole database.
        tmp_path = self.db_path + ".tmp"
        conn = duckdb_compat.connect(tmp_path)
        conn.execute("CREATE TABLE loci (LocusId TEXT PRIMARY KEY, Motif TEXT)")
        conn.execute("INSERT INTO loci VALUES ('a', 'AT')")
        conn.commit()
        conn.close()
        self.assertFalse(os.path.exists(tmp_path + ".wal"))
        os.replace(tmp_path, self.db_path)
        reopened = duckdb_compat.connect(self.db_path, read_only=True)
        try:
            self.assertEqual(reopened.execute("SELECT Motif FROM loci").fetchone()[0], "AT")
        finally:
            reopened.close()

    def test_find_leftover_write_files(self):
        self.assertEqual(duckdb_compat.find_leftover_write_files(self.db_path), [])
        with open(self.db_path + ".wal", "w") as f:
            f.write("x" * 10)
        self.assertEqual(duckdb_compat.find_leftover_write_files(self.db_path),
                         [(self.db_path + ".wal", 10)])


class ConnectionWrapperTests(unittest.TestCase):
    """What Connection passes down to the duckdb connection it wraps."""

    def test_executemany_hands_duckdb_a_list_rather_than_a_generator(self):
        # duckdb 1.0.x refuses a generator ("executemany requires a list of parameter sets
        # to be provided"), and requirements.txt allows duckdb 1.0, so the rows have to be
        # materialized before they go down. The NaNs still become NULL on the way.
        def rows():
            yield ("a", float("nan"))
            yield ("b", 5.0)

        recorder = RecordingDuckdbConnection()
        duckdb_compat.Connection(recorder).executemany("INSERT INTO t VALUES (?, ?)", rows())
        sql, params = recorder.executemany_calls[0]
        self.assertEqual(sql, "INSERT INTO t VALUES (?, ?)")
        self.assertIsInstance(params, list)
        self.assertEqual(params, [["a", None], ["b", 5.0]])

    def test_close_closes_the_duckdb_connection(self):
        recorder = RecordingDuckdbConnection()
        duckdb_compat.Connection(recorder).close()
        self.assertTrue(recorder.closed)


class InsertRowsTests(unittest.TestCase):
    """insert_rows is how the build loads every table, a chunk of rows per statement.

    It replaced a per-row executemany, so what matters is that the values still land
    exactly as they did: pandas widens a column of ints holding a None to float64 on its
    own, and an int written into a VARCHAR column would then arrive as "1.0". A chunk also
    holds thousands of rows, more than the sample DuckDB would type an untyped column from,
    so the tests below cover a column whose exceptional value sits past that sample.
    """

    def setUp(self):
        self.conn = duckdb_compat.connect(":memory:")

    def tearDown(self):
        self.conn.close()

    def test_values_round_trip_through_each_declared_type(self):
        self.conn.execute("CREATE TABLE t (i BIGINT, d DOUBLE, s VARCHAR)")
        duckdb_compat.insert_rows(
            self.conn, "t", ["i", "d", "s"],
            [{"i": 1, "d": 1.5, "s": "a"}, {"i": 2, "d": 2.5, "s": "b"}])
        self.assertEqual(self.conn.execute("SELECT i, d, s FROM t ORDER BY i").fetchall(),
                         [(1, 1.5, "a"), (2, 2.5, "b")])

    def test_an_int_beside_a_null_stays_an_int(self):
        # The column is VARCHAR because another row holds text, and this row's 1 has to
        # arrive as "1". Left to itself pandas would make the chunk float64 and write "1.0".
        self.conn.execute("CREATE TABLE t (v VARCHAR)")
        duckdb_compat.insert_rows(self.conn, "t", ["v"],
                                  [{"v": 1}, {"v": None}, {"v": "not_available"}])
        self.assertEqual(self.conn.execute("SELECT v FROM t").fetchall(),
                         [("1",), (None,), ("not_available",)])

    def test_a_mixed_column_keeps_a_value_that_sits_past_the_type_sample(self):
        # The case infer_column_types declares VARCHAR for: an annotation column carrying
        # counts plus a text marker. DuckDB types an object-dtype column off a sample of the
        # frame (pandas_analyze_sample, 1000 rows), so with the marker at row 2500 of 3000
        # it used to see only ints, type the column INT32, and fail on the marker; a 5.5 in
        # the same position was silently stored as "6".
        self.conn.execute("CREATE TABLE t (i BIGINT, v VARCHAR)")
        rows = [{"i": i, "v": i} for i in range(3000)]
        rows[2500]["v"] = "not_available"
        rows[2600]["v"] = 5.5
        rows[2700]["v"] = True
        duckdb_compat.insert_rows(self.conn, "t", ["i", "v"], rows)
        self.assertEqual(
            self.conn.execute("SELECT v FROM t WHERE i IN (0, 2500, 2600, 2700) ORDER BY i")
            .fetchall(),
            [("0",), ("not_available",), ("5.5",), ("true",)])

    def test_a_big_int_past_the_type_sample_is_not_truncated(self):
        # The reverse hazard, on a column that really is all ints: typed off the first rows
        # the column came out INT32, and the one value wider than that raised "Value out of
        # range for type INT" instead of being stored.
        self.conn.execute("CREATE TABLE t (v BIGINT)")
        rows = [{"v": i} for i in range(3000)]
        rows[2500]["v"] = 9007199254740995
        duckdb_compat.insert_rows(self.conn, "t", ["v"], rows)
        self.assertEqual(self.conn.execute("SELECT MAX(v) FROM t").fetchone()[0],
                         9007199254740995)

    def test_an_all_null_column_keeps_its_declared_type(self):
        # Routine: a build with no gene table writes pLI and the GeneTable* columns as all
        # NULL, and the server compares pLI to a number, which is a binder error unless the
        # column really is numeric.
        self.conn.execute("CREATE TABLE t (i BIGINT, d DOUBLE, s VARCHAR)")
        duckdb_compat.insert_rows(self.conn, "t", ["i", "d", "s"],
                                  [{"i": None, "d": None, "s": None}, {}])
        self.assertEqual(self.conn.execute("SELECT COUNT(*), COUNT(i), COUNT(d), COUNT(s) "
                                           "FROM t").fetchone(), (2, 0, 0, 0))
        self.assertEqual(self.conn.execute("SELECT COUNT(*) FROM t WHERE d > 0.5").fetchone()[0], 0)

    def test_a_value_that_does_not_fit_the_declared_type_names_the_column(self):
        self.conn.execute("CREATE TABLE t (v BIGINT)")
        with self.assertRaises(ValueError) as caught:
            duckdb_compat.insert_rows(self.conn, "t", ["v"], [{"v": 1}, {"v": "abc"}])
        self.assertIn('Column "v" is declared BIGINT', str(caught.exception))

    def test_a_large_int_is_not_rounded(self):
        # Wider than float64 can hold exactly, which is what a pandas upcast would cost.
        self.conn.execute("CREATE TABLE t (v BIGINT)")
        duckdb_compat.insert_rows(self.conn, "t", ["v"], [{"v": 9007199254740995}, {"v": None}])
        self.assertEqual(self.conn.execute("SELECT v FROM t ORDER BY v").fetchall(),
                         [(None,), (9007199254740995,)])

    def test_an_absent_key_is_null(self):
        self.conn.execute("CREATE TABLE t (a VARCHAR, b VARCHAR)")
        duckdb_compat.insert_rows(self.conn, "t", ["a", "b"], [{"a": "x"}])
        self.assertEqual(self.conn.execute("SELECT a, b FROM t").fetchall(), [("x", None)])

    def test_nan_becomes_null(self):
        # Same reason execute() and executemany() map it: a stored NaN is not NULL, is
        # counted by COUNT(col), and compares as larger than every real number.
        self.conn.execute("CREATE TABLE t (v DOUBLE)")
        duckdb_compat.insert_rows(self.conn, "t", ["v"], [{"v": float("nan")}, {"v": 5.0}])
        self.assertEqual(self.conn.execute("SELECT COUNT(v) FROM t").fetchone()[0], 1)
        self.assertEqual(self.conn.execute("SELECT COUNT(*) FROM t WHERE v IS NULL").fetchone()[0], 1)

    def test_rows_spanning_several_chunks_all_arrive(self):
        self.conn.execute("CREATE TABLE t (a BIGINT, b VARCHAR)")
        rows = [{"a": i, "b": f"row{i}"} for i in range(500)]
        with unittest.mock.patch.object(duckdb_compat, "_INSERT_CHUNK_CELLS", 10):
            duckdb_compat.insert_rows(self.conn, "t", ["a", "b"], rows)
        self.assertEqual(self.conn.execute("SELECT COUNT(*) FROM t").fetchone()[0], 500)
        self.assertEqual(self.conn.execute("SELECT b FROM t ORDER BY a LIMIT 1").fetchone()[0], "row0")
        self.assertEqual(self.conn.execute("SELECT b FROM t ORDER BY a DESC LIMIT 1").fetchone()[0],
                         "row499")

    def test_a_generator_of_rows_is_consumed_a_chunk_at_a_time(self):
        # The whole point of chunking: the caller never has to hold a second copy of the table.
        self.conn.execute("CREATE TABLE t (a BIGINT)")
        with unittest.mock.patch.object(duckdb_compat, "_INSERT_CHUNK_CELLS", 10):
            duckdb_compat.insert_rows(self.conn, "t", ["a"], ({"a": i} for i in range(25)))
        self.assertEqual(self.conn.execute("SELECT SUM(a) FROM t").fetchone()[0], sum(range(25)))

    def test_no_rows_writes_nothing(self):
        self.conn.execute("CREATE TABLE t (a BIGINT)")
        duckdb_compat.insert_rows(self.conn, "t", ["a"], [])
        self.assertEqual(self.conn.execute("SELECT COUNT(*) FROM t").fetchone()[0], 0)


class NaNParameterTests(unittest.TestCase):
    """A NaN bound as a parameter has to land as NULL, the way sqlite3 bound it.

    Missing numbers reach the database writer as NaN. DuckDB stores a bound NaN as an actual
    NaN: it is not NULL, COUNT(col) counts it, and it compares as LARGER than every real
    number. A locus with no population statistic would then pass the server's "affected
    allele > population 99th percentile" gate as though the cohort had data for it.
    """

    def setUp(self):
        self.conn = duckdb_compat.connect(":memory:")
        self.conn.execute("CREATE TABLE t (LocusId VARCHAR, p99 DOUBLE)")

    def tearDown(self):
        self.conn.close()

    def assert_matches_sqlite(self):
        import sqlite3
        reference = sqlite3.connect(":memory:")
        reference.execute("CREATE TABLE t (LocusId TEXT, p99 REAL)")
        reference.executemany("INSERT INTO t VALUES (?, ?)",
                              [("a", float("nan")), ("b", 5.0)])
        for query in ("SELECT COUNT(*) FROM t WHERE p99 IS NULL",
                      "SELECT COUNT(p99) FROM t",
                      "SELECT COUNT(*) FROM t WHERE 3.0 > p99"):
            self.assertEqual(self.conn.execute(query).fetchone()[0],
                             reference.execute(query).fetchone()[0], query)

    def test_executemany(self):
        self.conn.executemany("INSERT INTO t VALUES (?, ?)", [("a", float("nan")), ("b", 5.0)])
        self.assert_matches_sqlite()

    def test_execute(self):
        self.conn.execute("INSERT INTO t VALUES (?, ?)", ["a", float("nan")])
        self.conn.execute("INSERT INTO t VALUES (?, ?)", ["b", 5.0])
        self.assert_matches_sqlite()

    def test_numpy_nan(self):
        import numpy as np
        self.conn.executemany("INSERT INTO t VALUES (?, ?)",
                              [("a", np.float64("nan")), ("b", np.float64(5.0))])
        self.assert_matches_sqlite()

    def test_other_values_are_untouched(self):
        self.conn.executemany("INSERT INTO t VALUES (?, ?)",
                              [("a", None), ("b", 0.0), ("", float("inf"))])
        self.assertEqual(self.conn.execute("SELECT LocusId, p99 FROM t ORDER BY LocusId").fetchall(),
                         [("", float("inf")), ("a", None), ("b", 0.0)])

    def test_executemany_accepts_a_generator(self):
        def rows():
            for i in range(3):
                yield (f"locus{i}", float("nan"))
        self.conn.executemany("INSERT INTO t VALUES (?, ?)", rows())
        self.assertEqual(self.conn.execute("SELECT COUNT(*) FROM t WHERE p99 IS NULL").fetchone()[0], 3)


class SqlDialectTests(unittest.TestCase):
    """The two dialect differences the port had to fix at every call site."""

    def setUp(self):
        self.conn = duckdb_compat.connect(":memory:")
        self.conn.execute("CREATE TABLE t (symbol TEXT, size_diff TEXT)")
        self.conn.execute("INSERT INTO t VALUES ('FMR1', '13'), ('HTT', 'unknown'), ('ATXN1', NULL)")

    def tearDown(self):
        self.conn.close()

    def test_ilike_is_the_case_insensitive_match(self):
        # SQLite's LIKE was case-insensitive; DuckDB's is not, so the filters use ILIKE.
        self.assertEqual(self.conn.execute("SELECT symbol FROM t WHERE symbol ILIKE '%fmr%'").fetchall(),
                         [("FMR1",)])
        self.assertEqual(self.conn.execute("SELECT symbol FROM t WHERE symbol LIKE '%fmr%'").fetchall(),
                         [])

    def test_try_cast_keeps_unparseable_text_from_raising(self):
        # SQLite turned 'unknown' into 0 here; a plain CAST raises in DuckDB.
        with self.assertRaises(duckdb_compat.Error):
            self.conn.execute("SELECT CAST(size_diff AS BIGINT) FROM t").fetchall()
        # ATXN1 has no value, FMR1 has 13, HTT has text that will not parse.
        self.assertEqual(
            self.conn.execute("SELECT COALESCE(TRY_CAST(size_diff AS BIGINT), 0) FROM t "
                              "ORDER BY symbol").fetchall(),
            [(0,), (13,), (0,)])


if __name__ == "__main__":
    unittest.main()
