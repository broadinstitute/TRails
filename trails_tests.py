"""Unit tests for trails.py, the launcher.

The one thing here that a user cannot get back by re-running the build is their notes and
tags: they are typed in by hand in the results server, not derived from the input TSVs. The
DuckDB port renamed the file they live in, so these tests cover the migration that carries
them across.
"""

import os
import shutil
import sqlite3
import tempfile
import unittest

import duckdb_compat
import trails


def write_legacy_annotations(path, notes=(), tags=()):
    """Builds a SQLite annotations database shaped the way the pre-DuckDB TRails wrote it.

    Args:
        path: Where to create the file.
        notes: (locus_id, note_text, created_at, updated_at) tuples.
        tags: (locus_id, tag, created_at) tuples.
    """
    conn = sqlite3.connect(path)
    conn.execute("""CREATE TABLE notes (
        locus_id TEXT PRIMARY KEY,
        note_text TEXT NOT NULL,
        created_at TEXT NOT NULL,
        updated_at TEXT NOT NULL
    )""")
    conn.execute("""CREATE TABLE tags (
        locus_id TEXT NOT NULL,
        tag TEXT NOT NULL,
        created_at TEXT NOT NULL,
        PRIMARY KEY (locus_id, tag)
    )""")
    conn.executemany("INSERT INTO notes VALUES (?, ?, ?, ?)", list(notes))
    conn.executemany("INSERT INTO tags VALUES (?, ?, ?)", list(tags))
    conn.commit()
    conn.close()


def read_migrated(path):
    """Returns the (notes, tags) rows of a DuckDB annotations database."""
    conn = duckdb_compat.connect(path, read_only=True)
    try:
        return (conn.execute("SELECT locus_id, note_text, created_at, updated_at FROM notes "
                             "ORDER BY locus_id").fetchall(),
                conn.execute("SELECT locus_id, tag, created_at FROM tags "
                             "ORDER BY locus_id, tag").fetchall())
    finally:
        conn.close()


class BuildSourceFilesTests(unittest.TestCase):
    """The build fingerprint must cover the build's code and nothing else.

    A rebuild of the result database is expensive (the inputs are multi-gigabyte), so the
    fingerprint has to change when a module the build imports changes, and stay the same when
    a file the build never imports (the results server, a test file) changes.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.real_script_dir = trails.SCRIPT_DIR
        trails.SCRIPT_DIR = self.tmpdir
        # A miniature TRails directory: the two entry points, a module reached only through a
        # nested import two levels down, plus the two kinds of file the build never imports.
        self.write("trails.py", "import build_database\n")
        self.write("build_database.py", "def build():\n    import pipeline_step\n")
        self.write("pipeline_step.py", "from shared_helper import helper\n")
        self.write("shared_helper.py", "def helper():\n    return 1\n")
        self.write("results_server.py", "import flask\n")
        self.write("pipeline_step_tests.py", "import pipeline_step\n")

    def tearDown(self):
        trails.SCRIPT_DIR = self.real_script_dir
        shutil.rmtree(self.tmpdir)

    def write(self, name, text):
        """Writes text to name inside the fake script directory."""
        with open(os.path.join(self.tmpdir, name), "w") as f:
            f.write(text)

    def test_follows_imports_transitively_and_skips_everything_else(self):
        self.assertEqual(
            trails.build_source_files(),
            ["build_database.py", "pipeline_step.py", "shared_helper.py", "trails.py"])

    def test_editing_a_file_the_build_never_imports_does_not_change_the_fingerprint(self):
        before = trails.trails_code_version()
        self.write("results_server.py", "import flask\n\n# a server-only change\n")
        self.write("pipeline_step_tests.py", "import pipeline_step\n\n# another test\n")
        self.assertEqual(trails.trails_code_version(), before)

    def test_editing_a_build_module_changes_the_fingerprint(self):
        before = trails.trails_code_version()
        self.write("shared_helper.py", "def helper():\n    return 2\n")
        self.assertNotEqual(trails.trails_code_version(), before)

    def test_a_new_module_only_counts_once_a_build_module_imports_it(self):
        before = trails.trails_code_version()
        self.write("new_step.py", "def run():\n    return 3\n")
        self.assertEqual(trails.trails_code_version(), before)
        self.write("pipeline_step.py", "from shared_helper import helper\nimport new_step\n")
        after = trails.trails_code_version()
        self.assertIn("new_step.py", trails.build_source_files())
        self.assertNotEqual(after, before)

    def test_the_real_directory_covers_the_pipeline_and_excludes_server_and_tests(self):
        # Against the actual TRails checkout, not the fake directory above.
        trails.SCRIPT_DIR = self.real_script_dir
        names = trails.build_source_files()
        for expected in ["build_database.py", "input_tables.py", "mendelian_qc.py",
                         "motif_utilities.py", "result_database.py", "duckdb_compat.py"]:
            self.assertIn(expected, names)
        self.assertNotIn("results_server.py", names)
        self.assertEqual([n for n in names if n.endswith("_tests.py")], [])


class LegacyAnnotationsPathTests(unittest.TestCase):

    def test_derives_the_sqlite_path_the_previous_version_used(self):
        self.assertEqual(
            trails.legacy_annotations_db_path("/data/cohort.with_analysis_columns.duckdb"),
            "/data/cohort.with_analysis_columns.db.annotations.db")

    def test_a_db_path_that_is_already_sqlite_named_keeps_its_name(self):
        # --db can name anything; only the ".duckdb" the port introduced is swapped back.
        self.assertEqual(trails.legacy_annotations_db_path("/data/cohort.db"),
                         "/data/cohort.db.annotations.db")


class MigrateLegacyAnnotationsTests(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.annotations_db = os.path.join(self.tmpdir, "cohort.duckdb.annotations.duckdb")
        self.legacy_path = os.path.join(self.tmpdir, "cohort.db.annotations.db")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_notes_and_tags_arrive_intact(self):
        write_legacy_annotations(
            self.legacy_path,
            notes=[("1-100-110-AT", "candidate, check IGV", "2026-01-02 03:04:05",
                    "2026-01-09 10:11:12"),
                   ("2-1-9-CAG", "ruled out", "2026-02-02 02:02:02", "2026-02-02 02:02:02")],
            tags=[("1-100-110-AT", "followup", "2026-01-02 03:04:05"),
                  ("1-100-110-AT", "novel", "2026-01-03 03:04:05"),
                  ("2-1-9-CAG", "followup", "2026-02-02 02:02:02")])

        self.assertEqual(
            trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path), (2, 3))

        notes, tags = read_migrated(self.annotations_db)
        self.assertEqual(notes, [
            ("1-100-110-AT", "candidate, check IGV", "2026-01-02 03:04:05",
             "2026-01-09 10:11:12"),
            ("2-1-9-CAG", "ruled out", "2026-02-02 02:02:02", "2026-02-02 02:02:02")])
        self.assertEqual(tags, [
            ("1-100-110-AT", "followup", "2026-01-02 03:04:05"),
            ("1-100-110-AT", "novel", "2026-01-03 03:04:05"),
            ("2-1-9-CAG", "followup", "2026-02-02 02:02:02")])

    def test_the_legacy_file_is_left_untouched(self):
        write_legacy_annotations(
            self.legacy_path, notes=[("1-1-10-A", "keep me", "t1", "t2")])
        before = open(self.legacy_path, "rb").read()
        trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path)
        self.assertEqual(open(self.legacy_path, "rb").read(), before)

    def test_no_legacy_file_is_a_no_op(self):
        self.assertEqual(
            trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path), (0, 0))
        self.assertFalse(os.path.exists(self.annotations_db))

    def test_an_existing_annotations_database_is_never_overwritten(self):
        # The second run must not copy the old notes back over notes written since.
        write_legacy_annotations(self.legacy_path, notes=[("1-1-10-A", "old", "t1", "t2")])
        conn = duckdb_compat.connect(self.annotations_db)
        conn.execute("CREATE TABLE notes (locus_id TEXT PRIMARY KEY, note_text TEXT NOT NULL, "
                     "created_at TEXT NOT NULL, updated_at TEXT NOT NULL)")
        conn.execute("INSERT INTO notes VALUES ('1-1-10-A', 'new', 't3', 't4')")
        conn.commit()
        conn.close()

        self.assertEqual(
            trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path), (0, 0))
        conn = duckdb_compat.connect(self.annotations_db, read_only=True)
        try:
            self.assertEqual(conn.execute("SELECT note_text FROM notes").fetchone()[0], "new")
        finally:
            conn.close()

    def test_an_unreadable_legacy_file_does_not_stop_the_launcher(self):
        with open(self.legacy_path, "w") as f:
            f.write("this is not a SQLite database")
        self.assertEqual(
            trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path), (0, 0))
        self.assertFalse(os.path.exists(self.annotations_db))

    def test_an_empty_legacy_file_writes_no_database(self):
        write_legacy_annotations(self.legacy_path)
        self.assertEqual(
            trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path), (0, 0))
        self.assertFalse(os.path.exists(self.annotations_db))

    def test_tags_still_migrate_when_the_notes_table_is_missing(self):
        conn = sqlite3.connect(self.legacy_path)
        conn.execute("CREATE TABLE tags (locus_id TEXT NOT NULL, tag TEXT NOT NULL, "
                     "created_at TEXT NOT NULL, PRIMARY KEY (locus_id, tag))")
        conn.execute("INSERT INTO tags VALUES ('1-1-10-A', 'followup', 't1')")
        conn.commit()
        conn.close()

        self.assertEqual(
            trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path), (0, 1))
        notes, tags = read_migrated(self.annotations_db)
        self.assertEqual(notes, [])
        self.assertEqual(tags, [("1-1-10-A", "followup", "t1")])

    def test_a_legacy_path_with_uri_punctuation_is_read(self):
        # "#" and "?" end the path in a SQLite file: URI unless they are escaped, so an
        # unescaped path made SQLite open a different (empty) file and the user's notes and
        # tags were left behind. The destination is written to a plain directory here
        # because duckdb 1.5 cannot create a database under a directory whose name contains
        # "?": it splits the path at that character when it names the write-ahead log.
        awkward_dir = os.path.join(self.tmpdir, "run#3?draft")
        os.makedirs(awkward_dir)
        legacy_path = os.path.join(awkward_dir, "cohort.db.annotations.db")
        write_legacy_annotations(
            legacy_path,
            notes=[("1-100-110-AT", "candidate", "t1", "t2")],
            tags=[("1-100-110-AT", "followup", "t1")])

        self.assertEqual(
            trails.migrate_legacy_annotations(self.annotations_db, legacy_path), (1, 1))

        notes, tags = read_migrated(self.annotations_db)
        self.assertEqual(notes, [("1-100-110-AT", "candidate", "t1", "t2")])
        self.assertEqual(tags, [("1-100-110-AT", "followup", "t1")])

    def test_a_hash_in_the_data_directory_migrates_end_to_end(self):
        # The real layout: both databases sit in the user's data directory.
        awkward_dir = os.path.join(self.tmpdir, "run#3")
        os.makedirs(awkward_dir)
        legacy_path = os.path.join(awkward_dir, "cohort.db.annotations.db")
        annotations_db = os.path.join(awkward_dir, "cohort.duckdb.annotations.duckdb")
        write_legacy_annotations(
            legacy_path,
            notes=[("1-100-110-AT", "candidate", "t1", "t2")],
            tags=[("1-100-110-AT", "followup", "t1")])

        self.assertEqual(trails.migrate_legacy_annotations(annotations_db, legacy_path), (1, 1))

        notes, tags = read_migrated(annotations_db)
        self.assertEqual(notes, [("1-100-110-AT", "candidate", "t1", "t2")])
        self.assertEqual(tags, [("1-100-110-AT", "followup", "t1")])

    def test_no_temporary_file_is_left_behind(self):
        write_legacy_annotations(self.legacy_path, notes=[("1-1-10-A", "n", "t1", "t2")])
        trails.migrate_legacy_annotations(self.annotations_db, self.legacy_path)
        self.assertEqual(
            sorted(os.listdir(self.tmpdir)),
            sorted([os.path.basename(self.legacy_path), os.path.basename(self.annotations_db)]))


if __name__ == "__main__":
    unittest.main()
