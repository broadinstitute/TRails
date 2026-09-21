"""Unit tests for the sample-count metadata table written by result_database.py.

``build_database`` records the size of the callset once, and the results server reads it back
at startup to print "N loci from N genotyped samples" and to report ``total_samples`` through
/api/v1/schema. The three things worth pinning down are that the count round-trips, that
re-recording it replaces the row rather than adding a second one, and that a database built
before the table existed reads back as None instead of raising.
"""

import os
import tempfile
import unittest

import duckdb_compat
import result_database


class SampleCountTableTests(unittest.TestCase):

    def setUp(self):
        self.connection = duckdb_compat.connect(":memory:")

    def tearDown(self):
        self.connection.close()

    def test_database_without_metadata_table(self):
        # A database built before the metadata table existed has no such table at all, which is
        # not an error: the server shows the loci count on its own.
        self.assertIsNone(result_database.read_sample_count(self.connection))

    def test_write_then_read(self):
        result_database.write_sample_count(self.connection, 2657)
        self.assertEqual(result_database.read_sample_count(self.connection), 2657)

    def test_rewrite_replaces_value(self):
        result_database.write_sample_count(self.connection, 2136)
        result_database.write_sample_count(self.connection, 2657)
        self.assertEqual(result_database.read_sample_count(self.connection), 2657)
        self.assertEqual(
            self.connection.execute(
                f"SELECT COUNT(*) FROM {result_database.METADATA_TABLE}").fetchone()[0], 1)

    def test_metadata_table_without_a_sample_count_row(self):
        # The table exists but nothing recorded the count, e.g. a future build that writes some
        # other key first. The read must report "unknown", not blow up on the missing row.
        self.connection.execute(
            f"CREATE TABLE {result_database.METADATA_TABLE} (key TEXT PRIMARY KEY, value VARCHAR)")
        self.connection.execute(
            f"INSERT INTO {result_database.METADATA_TABLE} (key, value) VALUES (?, ?)",
            ("something_else", "1"))
        self.assertIsNone(result_database.read_sample_count(self.connection))

    def test_zero_samples_is_not_confused_with_absent(self):
        # 0 is falsy, so a naive read would report an empty callset the same way it reports a
        # database that never recorded the count.
        result_database.write_sample_count(self.connection, 0)
        self.assertEqual(result_database.read_sample_count(self.connection), 0)
        self.assertIsNotNone(result_database.read_sample_count(self.connection))


class SampleCountSurvivesFinalizeTests(unittest.TestCase):
    """write_sample_count leaves the commit to its caller, which is finalize_database."""

    def test_count_survives_the_atomic_move(self):
        with tempfile.TemporaryDirectory() as directory:
            final_path = os.path.join(directory, "result.duckdb")
            connection, tmp_path = result_database.open_new_database(final_path)
            result_database.write_loci_table(
                connection, [{"LocusId": "1-1-10-A", "Motif": "A"}], ["LocusId", "Motif"])
            result_database.write_sample_count(connection, 2657)
            result_database.finalize_database(connection, tmp_path, final_path)

            reopened = duckdb_compat.connect(final_path)
            try:
                self.assertEqual(result_database.read_sample_count(reopened), 2657)
            finally:
                reopened.close()


if __name__ == "__main__":
    unittest.main()
