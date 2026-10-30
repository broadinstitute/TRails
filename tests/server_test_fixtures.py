"""Fixtures shared by the results-server test modules.

This module's name does not end in "tests.py", so test discovery imports it only through the test
modules that use it.
"""

import os
import shutil
import tempfile
import unittest

import results_server
import swim_plot


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
