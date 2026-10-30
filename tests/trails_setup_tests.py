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


if __name__ == "__main__":
    unittest.main()
