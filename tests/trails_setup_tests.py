"""Unit tests for trails_setup, the reference-data download, with the network replaced by stubs."""

import contextlib
import gzip
import inspect
import io
import json
import os
import shutil
import tempfile
import unittest
from unittest import mock

import trails_launcher
import trails_setup
import trexplorer_annotations


class ResponseThatFailsAfterOneChunk(io.BytesIO):
    """A stub HTTP response that returns one chunk, then raises error.

    That is what a dropped connection, or Ctrl-C, does part way through a download.
    """

    def __init__(self, error):
        super().__init__(b"the first chunk of the table")
        self.error = error

    def read(self, size=-1):
        if self.tell():
            raise self.error
        return super().read(size)


class DownloadCleanupTests(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_a_download_that_fails_part_way_leaves_no_partial_file(self):
        with mock.patch.object(trails_setup.urllib.request, "urlopen",
                               return_value=ResponseThatFailsAfterOneChunk(ConnectionResetError("reset"))):
            with self.assertRaises(ConnectionResetError):
                trails_setup._download("https://example.org/table.tsv.gz", os.path.join(self.tmpdir, "table.tsv.gz"))
        self.assertEqual(os.listdir(self.tmpdir), [])

    def test_an_interrupted_trexplorer_download_leaves_no_partial_file(self):
        database_path = os.path.join(self.tmpdir, trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_NAME)
        with mock.patch.object(trails_setup.urllib.request, "urlopen",
                               return_value=ResponseThatFailsAfterOneChunk(KeyboardInterrupt())):
            with self.assertRaises(KeyboardInterrupt):
                trails_setup._download_trexplorer_annotations(trails_setup.TREXPLORER_ANNOTATIONS_TSV_URL,
                                                              database_path)
        self.assertEqual(os.listdir(self.tmpdir), [])


class GenelookupExportTests(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.dest = os.path.join(self.tmpdir, trails_setup.GENELOOKUP_GENE_TABLE_NAME)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def export(self, query_response, parts=()):
        """Runs the export against a stub service: query_response, then one gzipped part per URL."""
        responses = [io.BytesIO(json.dumps(query_response).encode())] + [
            io.BytesIO(gzip.compress(part.encode())) for part in parts]
        with mock.patch.object(trails_setup.urllib.request, "urlopen", side_effect=responses) as urlopen:
            trails_setup._export_genelookup_gene_table(trails_setup.GENELOOKUP_GENE_TABLE_SQL, self.dest)
        return urlopen

    def test_parts_are_joined_under_one_header(self):
        urlopen = self.export({"status": "success", "public_urls": ["https://a", "https://b"]},
                              ["gene_id\tgene_symbol\nENSG1\tA\n", "gene_id\tgene_symbol\nENSG2\tB\n"])
        with gzip.open(self.dest, "rt") as f:
            self.assertEqual(f.read(), "gene_id\tgene_symbol\nENSG1\tA\nENSG2\tB\n")
        query = urlopen.call_args_list[0].args[0]
        self.assertEqual(query.full_url, trails_setup.GENELOOKUP_QUERY_URL)
        self.assertEqual(json.loads(query.data), {"sql": trails_setup.GENELOOKUP_GENE_TABLE_SQL,
                                                  "export_to_file_format": "TSV"})

    def test_a_failed_export_raises_and_writes_nothing(self):
        with self.assertRaisesRegex(RuntimeError, "not permitted"):
            self.export({"error": "Export is not permitted for these columns"})
        self.assertFalse(os.path.exists(self.dest))

    def test_the_query_leaves_out_omim_columns(self):
        # genelookup refuses to export these; asking for one would fail the whole export
        self.assertNotIn("OMIM", trails_setup.GENELOOKUP_GENE_TABLE_SQL)


class FetchReferenceDataTests(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.fetched = []

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def fake_fetch(self, source, dest):
        self.fetched.append(source)
        with open(dest, "w") as f:
            f.write(source)

    def fetch(self, source, force=False):
        with mock.patch.object(trails_setup, "REFERENCE_FILES", [("table.tsv", source, self.fake_fetch, None)]):
            return trails_setup.fetch_reference_data(self.tmpdir, force=force)

    def test_each_source_is_fetched_once_until_it_changes_or_is_forced(self):
        self.assertEqual(self.fetch("v1"), {"table.tsv": os.path.join(self.tmpdir, "table.tsv")})
        self.fetch("v1")
        self.assertEqual(self.fetched, ["v1"])
        self.fetch("v2")
        self.fetch("v2", force=True)
        self.assertEqual(self.fetched, ["v1", "v2", "v2"])

    def test_a_failed_fetch_keeps_going_without_the_file(self):
        def failing_fetch(source, dest):
            raise OSError("offline")
        with mock.patch.object(trails_setup, "REFERENCE_FILES", [("table.tsv", "v1", failing_fetch, None)]):
            self.assertEqual(trails_setup.fetch_reference_data(self.tmpdir), {})

    def test_a_failed_fetch_names_a_trails_option_only_for_a_file_that_has_one(self):
        def failing_fetch(source, dest):
            raise OSError("offline")
        output = io.StringIO()
        with mock.patch.object(trails_setup, "REFERENCE_FILES",
                               [("table.tsv", "v1", failing_fetch, "--gene-table"),
                                ("database.duckdb", "v1", failing_fetch, None)]), \
                contextlib.redirect_stdout(output):
            trails_setup.fetch_reference_data(self.tmpdir)
        warnings = [line for line in output.getvalue().splitlines() if "[WARN]" in line]
        self.assertIn("trails.py --gene-table", warnings[0])
        self.assertNotIn("trails.py", warnings[1])

    def test_every_trails_option_named_is_one_trails_py_accepts(self):
        launcher_main_source = inspect.getsource(trails_launcher.main)
        for name, _source, _fetch, trails_option in trails_setup.REFERENCE_FILES:
            if trails_option:
                self.assertIn(f'"{trails_option}"', launcher_main_source, name)

    def test_the_downloads_include_trexplorer_and_genelookup(self):
        names = [name for name, _source, _fetch, _trails_option in trails_setup.REFERENCE_FILES]
        self.assertIn(trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_NAME, names)
        self.assertIn(trails_setup.GENELOOKUP_GENE_TABLE_NAME, names)
        self.assertIn(trails_setup.TREXPLORER_CATALOG_VERSION, trails_setup.TREXPLORER_ANNOTATIONS_TSV_URL)

    def test_every_downloaded_file_has_the_name_constant_trails_py_looks_it_up_by(self):
        # trails_launcher.main finds each file in fetch_reference_data's result by these names, and
        # a name that matched nothing there would silently turn that input off.
        self.assertEqual(
            sorted(name for name, _source, _fetch, _trails_option in trails_setup.REFERENCE_FILES),
            sorted([trails_setup.HPO_GENES_TO_PHENOTYPE_NAME, trails_setup.STR_ANALYSIS_KNOWN_LOCI_CATALOG_NAME,
                    trails_setup.STRCHIVE_LOCI_NAME, trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_NAME,
                    trails_setup.GENELOOKUP_GENE_TABLE_NAME]))


class TRExplorerAnnotationsDatabaseCacheTests(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_a_database_recorded_under_another_format_is_rebuilt_from_the_same_url(self):
        # A cached database recorded without the current database format may lack columns that
        # add_trexplorer_annotations reads, so it is rebuilt. The download URL stays the plain one.
        name, source, fetch, trails_option = next(entry for entry in trails_setup.REFERENCE_FILES
                                                  if entry[0] == trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_NAME)
        self.assertIn(f"database_format={trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_FORMAT_VERSION}", source)
        database_path = os.path.join(self.tmpdir, name)
        manifest_path = os.path.join(self.tmpdir, trails_setup.CACHE_MANIFEST_NAME)
        with open(database_path, "w") as f:
            f.write("a database built in an older format")
        with open(manifest_path, "w") as f:
            json.dump({name: trails_setup.TREXPLORER_ANNOTATIONS_TSV_URL}, f)

        with mock.patch.object(trails_setup, "REFERENCE_FILES", [(name, source, fetch, trails_option)]), \
                mock.patch.object(trails_setup, "_download") as download, \
                mock.patch.object(trexplorer_annotations, "build_trexplorer_annotations_database") as build:
            self.assertEqual(trails_setup.fetch_reference_data(self.tmpdir), {name: database_path})

        downloaded_table = database_path + ".download.tsv.gz"
        download.assert_called_once_with(trails_setup.TREXPLORER_ANNOTATIONS_TSV_URL, downloaded_table)
        build.assert_called_once_with(downloaded_table, database_path)
        with open(manifest_path) as f:
            self.assertEqual(json.load(f), {name: source})


if __name__ == "__main__":
    unittest.main()
