"""Unit tests for trails_launcher.py, the launcher the top-level trails.py hands over to."""

import argparse
import json
import os
import shutil
import tempfile
import unittest
from unittest import mock

import trails_launcher


class BuildSourceFilesTests(unittest.TestCase):
    """The build fingerprint must cover the build's code and nothing else.

    A rebuild of the result database is expensive (the inputs are multi-gigabyte), so the
    fingerprint has to change when a module the build imports changes, and stay the same when
    a file the build never imports (the results server, a test file) changes.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.real_script_dir = trails_launcher.SCRIPT_DIR
        trails_launcher.SCRIPT_DIR = self.tmpdir
        # A miniature TRails directory: the two entry points, a module reached only through a
        # nested import two levels down, plus the two kinds of file the build never imports.
        self.write("trails_launcher.py", "import build_database\n")
        self.write("build_database.py", "def build():\n    import pipeline_step\n")
        self.write("pipeline_step.py", "from shared_helper import helper\n")
        self.write("shared_helper.py", "def helper():\n    return 1\n")
        self.write("results_server.py", "import flask\n")
        self.write("pipeline_step_tests.py", "import pipeline_step\n")

    def tearDown(self):
        trails_launcher.SCRIPT_DIR = self.real_script_dir
        shutil.rmtree(self.tmpdir)

    def write(self, name, text):
        """Writes text to name inside the fake script directory."""
        with open(os.path.join(self.tmpdir, name), "w") as f:
            f.write(text)

    def test_follows_imports_transitively_and_skips_everything_else(self):
        self.assertEqual(
            trails_launcher.build_source_files(),
            ["build_database.py", "pipeline_step.py", "shared_helper.py", "trails_launcher.py"])

    def test_editing_a_file_the_build_never_imports_does_not_change_the_fingerprint(self):
        before = trails_launcher.trails_code_version()
        self.write("results_server.py", "import flask\n\n# a server-only change\n")
        self.write("pipeline_step_tests.py", "import pipeline_step\n\n# another test\n")
        self.assertEqual(trails_launcher.trails_code_version(), before)

    def test_editing_a_build_module_changes_the_fingerprint(self):
        before = trails_launcher.trails_code_version()
        self.write("shared_helper.py", "def helper():\n    return 2\n")
        self.assertNotEqual(trails_launcher.trails_code_version(), before)

    def test_a_new_module_only_counts_once_a_build_module_imports_it(self):
        before = trails_launcher.trails_code_version()
        self.write("new_step.py", "def run():\n    return 3\n")
        self.assertEqual(trails_launcher.trails_code_version(), before)
        self.write("pipeline_step.py", "from shared_helper import helper\nimport new_step\n")
        after = trails_launcher.trails_code_version()
        self.assertIn("new_step.py", trails_launcher.build_source_files())
        self.assertNotEqual(after, before)

    def test_the_real_directory_covers_the_pipeline_and_excludes_server_and_tests(self):
        # Against the actual TRails checkout, not the fake directory above.
        trails_launcher.SCRIPT_DIR = self.real_script_dir
        names = trails_launcher.build_source_files()
        for expected in ["build_database.py", "input_tables.py", "mendelian_qc.py",
                         "motif_utilities.py", "result_database.py", "duckdb_compat.py"]:
            self.assertIn(expected, names)
        self.assertNotIn("results_server.py", names)
        self.assertEqual([n for n in names if n.endswith("_tests.py")], [])


class BuildManifestFlagsTests(unittest.TestCase):
    """The manifest decides whether a database is reused, so its flags must compare exactly."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.tmpdir, "result.duckdb")
        self.manifest_path = self.db_path + ".build_manifest.json"
        with open(self.db_path, "w") as f:
            f.write("a database")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def manifest_for(self, population_gate_cohorts):
        """Returns the manifest trails_launcher.main would build for this gate and no input files."""
        args = argparse.Namespace(n_loci=None, skip_phenotype_scores=False, source_label="",
                                  population_gate_cohorts=population_gate_cohorts)
        return trails_launcher.build_manifest({}, trails_launcher.build_flags(args))

    def write_manifest(self, manifest):
        with open(self.manifest_path, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)

    def test_an_unchanged_gate_reuses_the_database(self):
        self.write_manifest(self.manifest_for(("HPRC256", "AoU1027")))
        self.assertFalse(trails_launcher.needs_rebuild(
            self.db_path, self.manifest_path, self.manifest_for(("HPRC256", "AoU1027")), False))

    def test_a_different_gate_forces_a_rebuild(self):
        self.write_manifest(self.manifest_for(("HPRC256", "AoU1027")))
        self.assertTrue(trails_launcher.needs_rebuild(
            self.db_path, self.manifest_path, self.manifest_for(("HPRC256", "AoU1027", "TenK10K")), False))


class VirtualEnvironmentTests(unittest.TestCase):
    """The launcher always runs with the Python of TRails' own virtual environment."""

    def run_check(self, prefix, environment_works):
        """Returns the os.execv arguments the check would use, or None if it returned."""
        with mock.patch.object(trails_launcher.sys, "prefix", prefix), \
                mock.patch.object(trails_launcher, "virtual_environment_python_works", return_value=environment_works), \
                mock.patch.object(trails_launcher.subprocess, "run") as run, \
                mock.patch.object(trails_launcher.os, "execv") as execv:
            run.return_value.returncode = 0
            trails_launcher.ensure_running_in_trails_virtual_environment()
        self.created = run.called
        return execv.call_args.args if execv.called else None

    def test_already_inside_returns(self):
        self.assertIsNone(self.run_check(trails_launcher.VIRTUAL_ENVIRONMENT_DIRECTORY, True))
        self.assertFalse(self.created)

    def test_outside_reruns_with_the_environments_python(self):
        python, argv = self.run_check("/usr/local", True)
        self.assertEqual(python, trails_launcher.VIRTUAL_ENVIRONMENT_PYTHON)
        self.assertEqual(argv, [trails_launcher.VIRTUAL_ENVIRONMENT_PYTHON] + trails_launcher.sys.argv)
        self.assertFalse(self.created)

    def test_a_missing_or_broken_environment_is_created_first(self):
        self.assertIsNotNone(self.run_check("/usr/local", False))
        self.assertTrue(self.created)


class InitialRequireAbovePopulationTests(unittest.TestCase):
    """The launcher validates --initial-require-above-population itself and forwards it to the server."""

    def test_choices_match_the_server(self):
        # Imported here, not at the top: only this test needs flask.
        import results_server
        self.assertEqual(trails_launcher.INITIAL_REQUIRE_ABOVE_POPULATION_CHOICES,
                         list(results_server.POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION))

    def server_command_for(self, initial_require_above_population):
        """Returns the results_server.py command line serve() would run."""
        args = argparse.Namespace(sample_metadata_tsv="samples.tsv", port=None, host=None,
                                  initial_require_above_population=initial_require_above_population)
        real_run = trails_launcher.subprocess.run
        commands = []
        trails_launcher.subprocess.run = lambda cmd, check: commands.append(cmd)
        try:
            trails_launcher.serve("/data/cohort.duckdb", args)
        finally:
            trails_launcher.subprocess.run = real_run
        return commands[0]

    def test_a_chosen_value_is_forwarded(self):
        command = self.server_command_for("all")
        flag_position = command.index("--initial-require-above-population")
        self.assertEqual(command[flag_position + 1], "all")

    def test_no_value_forwards_no_flag(self):
        self.assertNotIn("--initial-require-above-population", self.server_command_for(None))


if __name__ == "__main__":
    unittest.main()
