"""Unit tests for trails_setup, the reference-data download, with the network replaced by stubs."""

import gzip
import io
import json
import os
import shutil
import tempfile
import unittest
from unittest import mock

import trails_setup
import trexplorer_annotations


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
        with mock.patch.object(trails_setup, "REFERENCE_FILES", [("table.tsv", source, self.fake_fetch)]):
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
        with mock.patch.object(trails_setup, "REFERENCE_FILES", [("table.tsv", "v1", failing_fetch)]):
            self.assertEqual(trails_setup.fetch_reference_data(self.tmpdir), {})

    def test_the_downloads_include_trexplorer_and_genelookup(self):
        names = [name for name, _source, _fetch in trails_setup.REFERENCE_FILES]
        self.assertIn(trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_NAME, names)
        self.assertIn(trails_setup.GENELOOKUP_GENE_TABLE_NAME, names)
        self.assertIn(trails_setup.TREXPLORER_CATALOG_VERSION, trails_setup.TREXPLORER_ANNOTATIONS_TSV_URL)


class TRExplorerAnnotationsDatabaseCacheTests(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_a_database_recorded_under_another_format_is_rebuilt_from_the_same_url(self):
        # A cached database recorded without the current database format may lack columns that
        # add_trexplorer_annotations reads, so it is rebuilt. The download URL stays the plain one.
        name, source, fetch = next(entry for entry in trails_setup.REFERENCE_FILES
                                   if entry[0] == trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_NAME)
        self.assertIn(f"database_format={trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_FORMAT_VERSION}", source)
        database_path = os.path.join(self.tmpdir, name)
        manifest_path = os.path.join(self.tmpdir, trails_setup.CACHE_MANIFEST_NAME)
        with open(database_path, "w") as f:
            f.write("a database built in an older format")
        with open(manifest_path, "w") as f:
            json.dump({name: trails_setup.TREXPLORER_ANNOTATIONS_TSV_URL}, f)

        with mock.patch.object(trails_setup, "REFERENCE_FILES", [(name, source, fetch)]), \
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
