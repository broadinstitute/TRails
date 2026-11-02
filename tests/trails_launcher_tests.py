"""Unit tests for trails_launcher.py, the launcher the top-level trails.py hands over to."""

import argparse
import contextlib
import importlib.util
import io
import json
import os
import shutil
import subprocess
import sys
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
        # A miniature TRails directory: the launcher, the build entry point, a module reached only
        # through a nested import two levels down, a download-only module that only the launcher
        # imports, plus the two kinds of file the build never imports.
        self.write("trails_launcher.py", "import reference_downloads\n\ndef run_build():\n    import build_database\n")
        self.write("reference_downloads.py", "def fetch():\n    return 'url'\n")
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

    def test_editing_download_code_only_the_launcher_imports_does_not_change_the_fingerprint(self):
        # The reference files it downloads are hashed one by one in the manifest's inputs.
        before = trails_launcher.trails_code_version()
        self.write("reference_downloads.py", "def fetch():\n    return 'another url'\n")
        self.assertEqual(trails_launcher.trails_code_version(), before)

    def test_editing_the_launcher_itself_changes_the_fingerprint(self):
        # It decides which inputs and flags reach the build.
        before = trails_launcher.trails_code_version()
        self.write("trails_launcher.py", "import reference_downloads\n\ndef run_build():\n"
                                         "    import build_database\n    return 1\n")
        self.assertNotEqual(trails_launcher.trails_code_version(), before)

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
                         "motif_utilities.py", "result_database.py", "duckdb_compat.py",
                         "trexplorer_annotations.py", "trails_launcher.py"]:
            self.assertIn(expected, names)
        self.assertNotIn("results_server.py", names)
        self.assertNotIn("trails_setup.py", names)
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
                mock.patch.object(trails_launcher, "virtual_environment_python_and_pip_work",
                                  return_value=environment_works), \
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

    def test_the_optional_install_marker_lives_inside_the_environment(self):
        # Recreating the environment (venv --clear) removes the optional packages, so it has to
        # remove the marker too, or no later run would try to install them again.
        self.assertEqual(os.path.dirname(trails_launcher.OPTIONAL_INSTALL_MARKER),
                         trails_launcher.VIRTUAL_ENVIRONMENT_DIRECTORY)


class VirtualEnvironmentHealthCheckTests(unittest.TestCase):
    """An environment counts as working only if its Python starts and can import pip."""

    def check(self, python_path):
        with mock.patch.object(trails_launcher, "VIRTUAL_ENVIRONMENT_PYTHON", python_path):
            return trails_launcher.virtual_environment_python_and_pip_work()

    def test_a_missing_python_does_not_work(self):
        self.assertFalse(self.check("/nonexistent/.venv/bin/python"))

    def test_an_environment_left_without_pip_does_not_work(self):
        # What `python3 -m venv` leaves on Debian or Ubuntu without python3-venv: the interpreter
        # is created, then installing pip fails. Such an environment has to be recreated, or
        # every later pip install in it fails with "No module named pip".
        directory = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, directory)
        environment = os.path.join(directory, ".venv")
        subprocess.run([sys.executable, "-m", "venv", "--without-pip", environment], check=True)
        python_path = os.path.join(environment, "bin", "python")
        self.assertEqual(subprocess.run([python_path, "-c", "import sys"]).returncode, 0)
        self.assertFalse(self.check(python_path))

    def test_a_python_that_can_import_pip_works(self):
        if importlib.util.find_spec("pip") is None:
            self.skipTest("the Python running the tests has no pip")
        self.assertTrue(self.check(sys.executable))


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
        with mock.patch.object(trails_launcher.subprocess, "run",
                               return_value=subprocess.CompletedProcess([], 0)) as run, \
                contextlib.redirect_stdout(io.StringIO()):
            trails_launcher.serve("/data/cohort.duckdb", args)
        return run.call_args.args[0]

    def test_a_chosen_value_is_forwarded(self):
        command = self.server_command_for("all")
        flag_position = command.index("--initial-require-above-population")
        self.assertEqual(command[flag_position + 1], "all")

    def test_no_value_forwards_no_flag(self):
        self.assertNotIn("--initial-require-above-population", self.server_command_for(None))


class DatabaseDirectoryTests(unittest.TestCase):
    """An explicit --db in a directory that does not exist yet must not fail only at the build's end."""

    def test_the_db_directory_is_created_before_the_reference_download(self):
        directory = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, directory)
        db_path = os.path.join(directory, "results", "cohort.duckdb")
        directory_existed_when = {}

        def fetch_reference_data(reference_data_dir):
            directory_existed_when["the reference download started"] = os.path.isdir(os.path.dirname(db_path))
            return {}

        def run_build(args, db_path, *other_arguments):
            with open(db_path, "w") as f:
                f.write("a built database")

        for input_name in ("cohort.tsv", "samples.tsv"):
            with open(os.path.join(directory, input_name), "w") as f:
                f.write("a header line\n")
        argv = ["trails.py", "--repeat-copy-numbers-tsv", os.path.join(directory, "cohort.tsv"),
                "--sample-metadata-tsv", os.path.join(directory, "samples.tsv"), "--db", db_path,
                "--work-dir", directory, "--no-serve"]
        with mock.patch.object(trails_launcher.sys, "argv", argv), \
                mock.patch.object(trails_launcher, "ensure_running_in_trails_virtual_environment"), \
                mock.patch.object(trails_launcher, "ensure_dependencies"), \
                mock.patch.object(trails_launcher.trails_setup, "fetch_reference_data", side_effect=fetch_reference_data), \
                mock.patch.object(trails_launcher, "run_build", side_effect=run_build), \
                contextlib.redirect_stdout(io.StringIO()):
            trails_launcher.main()
        self.assertEqual(directory_existed_when, {"the reference download started": True})
        self.assertTrue(os.path.exists(db_path + ".build_manifest.json"))


class InputFileTests(unittest.TestCase):
    """A mistyped input path stops trails.py before the virtual environment, pip and the download."""

    def run_main(self, *extra_arguments):
        """Runs main() with a real sample table, a missing matrix and extra_arguments.

        Returns:
            The SystemExit code, and whether any slow step ran.
        """
        directory = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, directory)
        samples_path = os.path.join(directory, "samples.tsv")
        with open(samples_path, "w") as f:
            f.write("sample_id\nS1\n")
        argv = ["trails.py", "--repeat-copy-numbers-tsv", os.path.join(directory, "cohort.tsv"),
                "--sample-metadata-tsv", samples_path, "--work-dir", os.path.join(directory, "work"),
                *extra_arguments]
        with mock.patch.object(trails_launcher.sys, "argv", argv), \
                mock.patch.object(trails_launcher, "ensure_running_in_trails_virtual_environment") as venv, \
                mock.patch.object(trails_launcher, "ensure_dependencies") as dependencies, \
                mock.patch.object(trails_launcher.trails_setup, "fetch_reference_data") as fetch:
            with self.assertRaises(SystemExit) as exit:
                trails_launcher.main()
        self.assertFalse(os.path.exists(os.path.join(directory, "work")))
        return exit.exception.code, venv.called or dependencies.called or fetch.called

    def test_a_missing_matrix_stops_before_any_slow_step_naming_the_option_and_the_path(self):
        message, a_slow_step_ran = self.run_main()
        self.assertFalse(a_slow_step_ran)
        self.assertRegex(message, r"^ERROR: --repeat-copy-numbers-tsv: .*cohort\.tsv is not an existing file$")

    def test_every_missing_file_is_named_including_the_optional_overrides(self):
        message, a_slow_step_ran = self.run_main("--gene-table", "/no/such/genes.tsv.gz")
        self.assertFalse(a_slow_step_ran)
        self.assertEqual([line.split(":")[1].strip() for line in message.splitlines()],
                         ["--repeat-copy-numbers-tsv", "--gene-table"])


class ServerExitTests(unittest.TestCase):
    """How trails.py ends when the results server stops: no Python traceback either way."""

    def serve(self, **run_behavior):
        """Runs serve() against a stubbed server; returns (SystemExit code or None, stdout, stderr)."""
        args = argparse.Namespace(sample_metadata_tsv="samples.tsv", port=None, host=None,
                                  initial_require_above_population=None)
        stdout, stderr = io.StringIO(), io.StringIO()
        exit_code = None
        with mock.patch.object(trails_launcher.subprocess, "run", **run_behavior), \
                contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                trails_launcher.serve("/data/cohort.duckdb", args)
            except SystemExit as exit:
                exit_code = exit.code
        return exit_code, stdout.getvalue(), stderr.getvalue()

    def test_ctrl_c_ends_quietly(self):
        # Ctrl-C stops the server and reaches the launcher too, as a KeyboardInterrupt.
        self.assertEqual(self.serve(side_effect=KeyboardInterrupt)[0::2], (None, ""))

    def test_a_server_that_stops_normally_ends_quietly(self):
        self.assertEqual(self.serve(return_value=subprocess.CompletedProcess([], 0))[0::2], (None, ""))

    def test_a_server_error_exits_with_its_status_and_one_line(self):
        # e.g. the port is already in use, or the database is missing; the server has said why.
        exit_code, _stdout, stderr = self.serve(return_value=subprocess.CompletedProcess([], 3))
        self.assertEqual(exit_code, 3)
        self.assertEqual(stderr.count("\n"), 1)
        self.assertIn("exit status 3", stderr)

    def test_a_server_killed_by_a_signal_exits_with_128_plus_the_signal_number(self):
        self.assertEqual(self.serve(return_value=subprocess.CompletedProcess([], -15))[0], 143)


if __name__ == "__main__":
    unittest.main()
