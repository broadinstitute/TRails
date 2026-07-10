"""Unit tests for skinny_tables.py."""

import sqlite3
import unittest

import skinny_tables


def make_loci_connection(column_definitions, rows):
    """Creates an in-memory database with a minimal ``loci`` table.

    Args:
        column_definitions: An ordered list of ``(column_name, sql_type)`` pairs.
        rows: A list of value-tuples to insert (one per locus), aligned to
            ``column_definitions``.

    Returns:
        An open sqlite3 connection holding the populated ``loci`` table.
    """
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "CREATE TABLE loci ({})".format(
            ", ".join(f'"{name}" {sql_type}' for name, sql_type in column_definitions)))
    if rows:
        connection.executemany(
            "INSERT INTO loci VALUES ({})".format(
                ", ".join("?" for _ in column_definitions)),
            rows)
    connection.commit()
    return connection


def table_columns(connection, table_name):
    """Returns the ordered list of column names of ``table_name``."""
    return [row[1] for row in connection.execute(
        f"SELECT * FROM pragma_table_info('{table_name}')")]


def table_index_names(connection, table_name):
    """Returns the list of index names defined on ``table_name``."""
    return [row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='index' AND tbl_name=?",
        (table_name,))]


class TestColumnDefinitions(unittest.TestCase):

    def test_shared_columns_count_and_membership(self):
        shared = skinny_tables.shared_columns()
        # 21 entries, matching the reference add_skinny_tables.shared_columns()
        # exactly (the blueprint's "30" is an off-by-one; live source is source
        # of truth), with the AoU Phase 2 columns removed.
        self.assertEqual(len(shared), 21)
        self.assertEqual(len(set(shared)), 21)
        self.assertEqual(shared[0], "LocusId")
        for expected in ("Chrom", "Motif", "CanonicalMotif", "pLI",
                         "GeneTableGeneSymbol", "AoU1027_StdevPercentile"):
            self.assertIn(expected, shared)

    def test_per_outlier_columns_count_and_suffix(self):
        for outlier_type in skinny_tables.OUTLIER_TYPES:
            per_outlier = skinny_tables.per_outlier_columns(outlier_type)
            self.assertEqual(len(per_outlier), 13)
            self.assertEqual(len(set(per_outlier)), 13)
            for column in per_outlier:
                # Every per-outlier column carries the outlier type; the three
                # _ByFamily columns carry it before the _ByFamily suffix.
                self.assertIn(f"_{outlier_type}", column)
            self.assertIn(f"FirstAffectedAlleleSize_{outlier_type}", per_outlier)
            self.assertIn(f"OutlierSampleIds_{outlier_type}", per_outlier)
            self.assertIn(f"SecondAffectedAlleleSize_{outlier_type}_ByFamily", per_outlier)

    def test_all_alleles_total_is_34(self):
        self.assertEqual(
            len(skinny_tables.shared_columns()
                + skinny_tables.per_outlier_columns("AllAlleles")),
            34)


class TestBuildSkinnyTable(unittest.TestCase):

    def test_projects_only_columns_present_in_loci(self):
        # A subset of shared + AllAlleles per-outlier columns, plus an unrelated
        # column that must NOT appear in the skinny table.
        present = [
            "LocusId", "Chrom", "Motif", "pLI",
            "FirstAffectedAlleleSize_AllAlleles", "OutlierSampleIds_AllAlleles",
        ]
        column_definitions = [(name, "TEXT") for name in present] + [("UnrelatedColumn", "TEXT")]
        connection = make_loci_connection(
            column_definitions,
            [tuple(f"v{i}" for i in range(len(column_definitions)))])
        try:
            skinny_tables.build_skinny_table(
                connection, "AllAlleles", set(table_columns(connection, "loci")))
            self.assertEqual(table_columns(connection, "sk_AllAlleles"), present)
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM sk_AllAlleles").fetchone()[0], 1)
        finally:
            connection.close()

    def test_hemizygous_marker_is_one_or_null(self):
        column_definitions = [
            ("LocusId", "TEXT"),
            ("FirstAffectedAlleleSize_HemizygousAlleles", "INTEGER"),
            ("HemizygousAlleleHistogram", "BLOB"),
        ]
        rows = [
            ("locus_with_histogram", 10, "1:2|2:3"),
            ("locus_without_histogram", 20, None),
        ]
        connection = make_loci_connection(column_definitions, rows)
        try:
            skinny_tables.build_skinny_table(
                connection, "HemizygousAlleles", set(table_columns(connection, "loci")))
            self.assertIn("HemizygousAlleleHistogram",
                          table_columns(connection, "sk_HemizygousAlleles"))
            marker_by_locus = dict(connection.execute(
                "SELECT LocusId, HemizygousAlleleHistogram FROM sk_HemizygousAlleles"))
            self.assertEqual(marker_by_locus["locus_with_histogram"], 1)
            self.assertIsNone(marker_by_locus["locus_without_histogram"])
        finally:
            connection.close()

    def test_no_indexes_created(self):
        connection = make_loci_connection(
            [("LocusId", "TEXT"), ("Chrom", "TEXT")],
            [("locus1", "chr1")])
        try:
            for outlier_type in skinny_tables.OUTLIER_TYPES:
                skinny_tables.build_skinny_table(
                    connection, outlier_type, set(table_columns(connection, "loci")))
                self.assertEqual(
                    table_index_names(connection, f"sk_{outlier_type}"), [])
        finally:
            connection.close()

    def test_idempotent_rebuild(self):
        connection = make_loci_connection(
            [("LocusId", "TEXT")], [("locus1",), ("locus2",)])
        try:
            first = skinny_tables.build_skinny_table(
                connection, "ShortAlleles", set(table_columns(connection, "loci")))
            skinny_tables.build_skinny_table(
                connection, "ShortAlleles", set(table_columns(connection, "loci")))
            self.assertEqual(first, "sk_ShortAlleles")
            self.assertEqual(
                connection.execute("SELECT COUNT(*) FROM sk_ShortAlleles").fetchone()[0], 2)
        finally:
            connection.close()

    def test_hemizygous_marker_absent_when_histogram_column_absent(self):
        connection = make_loci_connection(
            [("LocusId", "TEXT")], [("locus1",)])
        try:
            skinny_tables.build_skinny_table(
                connection, "HemizygousAlleles", set(table_columns(connection, "loci")))
            self.assertNotIn("HemizygousAlleleHistogram",
                             table_columns(connection, "sk_HemizygousAlleles"))
        finally:
            connection.close()


if __name__ == "__main__":
    unittest.main()
