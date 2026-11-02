"""Fixtures shared by the results-server test modules.

This module's name does not end in "tests.py", so test discovery imports it only through the test
modules that use it.
"""

import gzip
import os
import shutil
import tempfile
import unittest

import duckdb_compat
import results_server
import swim_plot


def read_gzipped_tsv_rows(response):
    """Return the rows of a gzipped TSV download, as dicts keyed by its header line.

    Args:
        response: A Flask test-client response whose body is a gzipped TSV with a header line.

    Returns:
        List of dicts mapping each header column to the row's value, as text.
    """
    lines = gzip.decompress(response.get_data()).decode("utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def make_swim_plot_row(**values):
    """Return a full-width swim_plot row dict, defaulting every unnamed column to NULL.

    Args:
        **values: Values of swim_plot columns, by column name.

    Returns:
        dict mapping every swim_plot.SWIM_PLOT_COLUMNS column to its value or None.

    Raises:
        AssertionError: If a name is not a swim_plot column, so a misspelled column fails the
            test instead of being written as NULL.
    """
    unknown_columns = set(values) - set(swim_plot.SWIM_PLOT_COLUMNS)
    assert not unknown_columns, f"swim_plot has no column(s): {sorted(unknown_columns)}"
    return {column: values.get(column) for column in swim_plot.SWIM_PLOT_COLUMNS}


class ServerTestCase(unittest.TestCase):
    """Base class for test classes that use a temporary directory or change the server's app.config.

    configure_app rewrites app.config, which every test module shares, so the configuration is
    snapshotted before the class's tests and put back after them. Each class also gets its own
    temporary directory, cls.tmpdir, which is removed after its tests.
    """

    @classmethod
    def setUpClass(cls):
        saved_config = dict(results_server.app.config)
        # Class cleanups run last-registered first: the directory goes, then the config is restored
        cls.addClassCleanup(results_server.app.config.update, saved_config)
        cls.addClassCleanup(results_server.app.config.clear)
        cls.tmpdir = tempfile.mkdtemp()
        cls.addClassCleanup(shutil.rmtree, cls.tmpdir)

    @classmethod
    def configure_server(cls, db_path, **configure_app_arguments):
        """Point the server at a test database through configure_app, the shared startup path.

        Sets cls.annotations_db (in cls.tmpdir) and cls.client, a Flask test client.

        Args:
            db_path: Path to the test database.
            **configure_app_arguments: Further configure_app arguments, e.g. sample_table.
        """
        cls.annotations_db = os.path.join(cls.tmpdir, "annotations.duckdb")
        results_server.configure_app(db_path, annotations_db=cls.annotations_db, **configure_app_arguments)
        results_server.app.config["TESTING"] = True
        cls.client = results_server.app.test_client()

    def get_locus_ids(self, query):
        """Return the LocusIds /api/v1/loci returns for outlier_type=all and further parameters.

        Args:
            query: The further query parameters, e.g. "phenotype_keyword=ataxia".

        Returns:
            Set of the LocusIds on the first page of results.
        """
        response = self.client.get(f"/api/v1/loci?outlier_type=all&{query}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {r["LocusId"] for r in response.get_json()["results"]}

    def get_swim_plot_rows(self, query):
        """Return the rows /api/v1/swim_plot_data returns for outlier_type=all and further parameters.

        Args:
            query: The further query parameters, e.g. "phenotype_keyword=ataxia", URL-encoded
                where a value needs it. "" for none.

        Returns:
            Set of the (LocusId, sample_id) pairs of the returned rows.
        """
        response = self.client.get(f"/api/v1/swim_plot_data?outlier_type=all&{query}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return {(entry["LocusId"], entry["sample_id"]) for entry in response.get_json()["data"]}

    def get_locus_ids_from_api_query(self, db_path, **params):
        """Return the LocusIds build_api_query's select query returns for outlier_type=all and params.

        Runs the query directly against db_path, with the server's current app.config, and checks
        that its count query agrees with the number of rows.

        Args:
            db_path: Path to the test database.
            **params: Further validated query parameters, e.g. sort_by=["variation_cluster"].

        Returns:
            List of the LocusIds, in the query's order.
        """
        select_query, count_query, _, _, sql_params = results_server.build_api_query(
            dict({"outlier_type": "all"}, **params))
        conn = duckdb_compat.connect(db_path, read_only=True)
        conn.row_factory = duckdb_compat.Row
        try:
            locus_ids = [row["LocusId"] for row in conn.execute(select_query, sql_params)]
            self.assertEqual(conn.execute(count_query, sql_params).fetchone()[0], len(locus_ids))
        finally:
            conn.close()
        return locus_ids

    def leave_columns_out_of_the_database_for_one_test(self, *columns):
        """Make the server treat its database as lacking some loci columns, for one test.

        Args:
            *columns: The loci columns to leave out of app.config["DB_COLUMNS_SET"].
        """
        saved_columns = results_server.app.config["DB_COLUMNS_SET"]
        self.addCleanup(results_server.app.config.__setitem__, "DB_COLUMNS_SET", saved_columns)
        results_server.app.config["DB_COLUMNS_SET"] = saved_columns - set(columns)

    def record_template_render_arguments(self):
        """Make the pages record their template's arguments instead of rendering, for one test.

        Returns:
            dict that each page render updates with its keyword arguments. The page body is "".
        """
        render_arguments = {}

        class RecordingTemplate:
            def render(self, **kwargs):
                render_arguments.update(kwargs)
                return ""

        original_get_template = results_server.jinja2_env.get_template
        results_server.jinja2_env.get_template = lambda name: RecordingTemplate()
        self.addCleanup(setattr, results_server.jinja2_env, "get_template", original_get_template)
        return render_arguments
