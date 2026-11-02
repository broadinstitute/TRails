"""TRails — one command from your two tables to a running server.

    python3 trails.py \\
        --repeat-copy-numbers-tsv your_cohort.tsv.gz \\
        --sample-metadata-tsv     your_samples.tsv.gz \\
        [--phenotypes-table your_phenotypes.tsv.gz]

In order, this:
  1. re-runs itself with the Python of TRails' own virtual environment, .venv/ (creating it
     if needed), and installs any missing Python dependencies there (requirements.txt),
  2. downloads any missing public (class-A) reference data into reference_data/: the HPO and
     known-disease-locus files, TRExplorer's population statistics and locus annotations, and
     genelookup's gene table,
  3. (re)builds the DuckDB result database if your inputs changed (sha256 build manifest),
  4. starts the local results server (unless --no-serve).

Only --repeat-copy-numbers-tsv and --sample-metadata-tsv are required; the phenotype table is
optional and skipped when absent. See docs/INPUT_FORMATS.md.

Note: this orchestrates build_database.build(), which reads the TSVs and populates the DB
directly, with no intermediate files, and results_server.py (serve).
"""

import argparse
import ast
import hashlib
import importlib.util
import json
import os
import subprocess
import sys

# analysis_columns imports only the standard library, so it can be imported here, before
# ensure_dependencies has installed anything; build_database cannot (see run_build).
import analysis_columns
import trails_setup

# The trails/ code directory, and the TRails directory above it that holds requirements.txt,
# reference_data/ and the top-level trails.py users run.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TRAILS_DIRECTORY = os.path.dirname(SCRIPT_DIR)

# Where the build's import graph starts. The code fingerprint in the build manifest covers
# these files and every TRails module they import, directly or indirectly, so a new pipeline
# module is picked up automatically as soon as an existing one imports it. Add an entry here
# only for a new way of starting a build that does not go through one of these files.
BUILD_ENTRY_POINTS = ["build_database.py"]

# Files the fingerprint covers without following their imports. This launcher decides which
# inputs and flags reach the build (build_flags, run_build), so it is covered itself. What it
# imports beyond the build is trails_setup, which only downloads the reference files, and
# build_manifest hashes each of those on its own, so editing download code rebuilds nothing.
BUILD_FILES_HASHED_WITHOUT_THEIR_IMPORTS = ["trails_launcher.py"]

# import names of the runtime deps in requirements.txt that the code actually imports.
REQUIRED_MODULES = ["flask", "pandas", "numpy", "intervaltree", "duckdb"]

# Not needed on any default code path, so a missing one is a warning, not an exit: pyhpo gives
# ontology-aware phenotype similarity, and phenotype_scoring falls back to Jaccard without it.
OPTIONAL_MODULES = ["pyhpo"]

# TRails' own Python virtual environment, which install.sh creates and installs the dependencies
# into. Recent Homebrew and Linux system Pythons refuse pip installs outside a virtual environment
# ("externally-managed-environment"), so TRails never installs into the Python it was started with.
VIRTUAL_ENVIRONMENT_DIRECTORY = os.path.join(TRAILS_DIRECTORY, ".venv")
VIRTUAL_ENVIRONMENT_PYTHON = os.path.join(VIRTUAL_ENVIRONMENT_DIRECTORY, "bin", "python")

# Written after the one optional-dependency install attempt, so a package that will not build here
# does not make every later run re-try (and re-wait for) the same failing pip invocation. It lives
# inside the virtual environment, so recreating the environment, which removes every installed
# package, also removes the marker, and the next run tries the optional install again.
OPTIONAL_INSTALL_MARKER = os.path.join(VIRTUAL_ENVIRONMENT_DIRECTORY, ".trails_optional_deps_attempted")

# The values results_server.py --initial-require-above-population accepts (the keys of its
# POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION). Copied rather than imported: results_server
# imports flask and numpy, which may not be installed yet when these arguments are parsed, before
# ensure_dependencies runs. trails_launcher_tests.py checks the two agree.
INITIAL_REQUIRE_ABOVE_POPULATION_CHOICES = ["long-read", "short-read", "all"]


def _has_module(name):
    """Return True if an import of name would resolve, without importing it."""
    try:
        return importlib.util.find_spec(name) is not None
    except ModuleNotFoundError:
        return False


def virtual_environment_python_and_pip_work():
    """Return True if TRails' virtual environment exists and its Python starts and can import pip.

    A virtual environment keeps pointing at the Python it was made from, so it stops working
    when that Python is removed (e.g. by a Homebrew upgrade). One whose creation stopped part way
    has a Python but no pip: on Debian or Ubuntu without the python3-venv package, `venv` creates
    the interpreter and then fails to install pip. Either kind is recreated.
    """
    if not os.path.exists(VIRTUAL_ENVIRONMENT_PYTHON):
        return False
    return subprocess.run([VIRTUAL_ENVIRONMENT_PYTHON, "-c", "import pip"],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def ensure_running_in_trails_virtual_environment():
    """Re-run this command with the Python of TRails' virtual environment, creating it if needed.

    Returns only when already running in that environment; otherwise replaces this process.
    """
    if os.path.realpath(sys.prefix) == os.path.realpath(VIRTUAL_ENVIRONMENT_DIRECTORY):
        return
    if not virtual_environment_python_and_pip_work():
        print(f"Creating a Python virtual environment in {VIRTUAL_ENVIRONMENT_DIRECTORY} ...")
        if subprocess.run([sys.executable, "-m", "venv", "--clear", VIRTUAL_ENVIRONMENT_DIRECTORY]).returncode != 0:
            sys.exit("ERROR: could not create a Python virtual environment. On Debian or Ubuntu, "
                     "install the python3-venv package (e.g. sudo apt install python3-venv) and re-run.")
    os.execv(VIRTUAL_ENVIRONMENT_PYTHON, [VIRTUAL_ENVIRONMENT_PYTHON] + sys.argv)


def ensure_dependencies():
    """pip-install the runtime dependencies if any are missing.

    The required set is installed from requirements.txt with check=True: without it nothing can
    run. The optional set is installed separately, by name, with check=False: pip fails the whole
    invocation if any listed package cannot be built, so installing it with the required set, or
    checking its result, would turn a soft dependency into a hard one.
    """
    if not all(_has_module(m) for m in REQUIRED_MODULES):
        print("Installing Python dependencies (requirements.txt) ...")
        subprocess.run([sys.executable, "-m", "pip", "install", "-r",
                        os.path.join(TRAILS_DIRECTORY, "requirements.txt")], check=True)
        still_missing = [m for m in REQUIRED_MODULES if not _has_module(m)]
        if still_missing:
            sys.exit(f"ERROR: dependencies still missing after install: {', '.join(still_missing)}")

    missing_optional = [m for m in OPTIONAL_MODULES if not _has_module(m)]
    if missing_optional and not os.path.exists(OPTIONAL_INSTALL_MARKER):
        print(f"Installing optional Python dependencies ({', '.join(missing_optional)}) ...")
        subprocess.run([sys.executable, "-m", "pip", "install"] + missing_optional, check=False)
        # Record the attempt either way. These packages are optional and one of them may simply
        # not build on this machine; retrying the install on every invocation would make every run
        # pay for a failure that is not going to resolve itself. Delete the marker to retry.
        try:
            with open(OPTIONAL_INSTALL_MARKER, "w") as f:
                f.write(" ".join(OPTIONAL_MODULES) + "\n")
        except OSError:
            pass
        missing_optional = [m for m in OPTIONAL_MODULES if not _has_module(m)]
    if missing_optional:
        print(f"NOTE: optional dependencies not installed: {', '.join(missing_optional)}. "
              f"Phenotype scoring falls back to Jaccard similarity without pyhpo. "
              f"Delete {OPTIONAL_INSTALL_MARKER} to retry.")


def hash_file(path):
    """Return the hex sha256 of a file, or None if path is empty or missing."""
    if not path or not os.path.exists(path):
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def module_names_imported_by(path):
    """Return the set of top-level module names the Python file at path imports.

    Both "import x" and "from x import y" count, at any nesting level, so a TRails module a
    build module imports inside a function rather than at the top of the file is still
    covered. Relative imports are ignored: TRails is a flat directory of modules with no
    package.

    Args:
        path: Path of the .py file to parse.

    Returns:
        A set of module-name strings, with any dotted suffix stripped.
    """
    with open(path) as f:
        tree = ast.parse(f.read(), filename=path)
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names.add(node.module.split(".")[0])
    return names


def build_source_files():
    """Return the sorted names of the TRails sources the built database depends on.

    Starts from BUILD_ENTRY_POINTS and follows imports that resolve to a .py file next to
    this script, transitively, so a module a pipeline module starts importing is picked up
    on its own and nobody has to maintain a list, then adds BUILD_FILES_HASHED_WITHOUT_THEIR_IMPORTS.
    Files the build never imports (trails_setup.py, results_server.py, which is launched as
    a subprocess, and the *_tests.py files) stay out, so editing them does not invalidate the
    manifest.

    Returns:
        A sorted list of file names (not paths), all of which exist in SCRIPT_DIR.
    """
    found = {name for name in BUILD_FILES_HASHED_WITHOUT_THEIR_IMPORTS
             if os.path.exists(os.path.join(SCRIPT_DIR, name))}
    seen = set()
    queue = list(BUILD_ENTRY_POINTS)
    while queue:
        name = queue.pop()
        if name in seen:
            continue
        seen.add(name)
        path = os.path.join(SCRIPT_DIR, name)
        if not os.path.exists(path):
            continue
        found.add(name)
        queue.extend(f"{module}.py" for module in module_names_imported_by(path))
    return sorted(found)


def trails_code_version():
    """Return a fingerprint of the TRails code for the build manifest.

    A sha256 over the TRails Python sources the build uses, so re-running after a code
    update correctly triggers a rebuild instead of reusing a stale database. The git
    commit is deliberately not used: it does not change when the working tree is
    edited, so in a checkout (the normal case while developing) an uncommitted
    change to the pipeline would silently reuse the database built by the old code.
    """
    digest = hashlib.sha256()
    for name in build_source_files():
        file_hash = hash_file(os.path.join(SCRIPT_DIR, name))
        if file_hash:
            digest.update(name.encode())
            digest.update(file_hash.encode())
    return digest.hexdigest()


def default_db_path(tsv_path, work_dir):
    """Return the default result-DB path derived from the genotype TSV's basename in work_dir."""
    base = os.path.basename(tsv_path)
    for ext in (".gz", ".tsv", ".txt"):
        if base.endswith(ext):
            base = base[:-len(ext)]
    return os.path.join(work_dir, base + ".with_analysis_columns.duckdb")


def build_flags(args):
    """Return {name: value} for every trails.py flag that affects the built DB, for the manifest.

    Every value is one that survives the manifest's JSON round trip unchanged, since
    needs_rebuild compares the manifest read back from disk with a freshly built one.

    Args:
        args: The parsed trails.py command line.
    """
    return {"n_loci": args.n_loci, "skip_phenotype_scores": args.skip_phenotype_scores,
            "source_label": args.source_label,
            # Comma-joined rather than a tuple: the manifest read back from JSON holds a list,
            # which never equals a tuple, so a tuple here would force a rebuild on every run.
            "population_gate_cohorts": ",".join(args.population_gate_cohorts)}


def build_manifest(inputs, flags):
    """Return the build manifest, which needs_rebuild compares to tell whether the DB is stale.

    It holds the input sha256s, the build flags, the phenotype-scoring backend and a fingerprint
    of the build's source files (trails_code_version).

    Args:
        inputs: {label: path} for every file whose contents affect the built DB.
        flags: {name: value} for every build flag that affects the DB.
    """
    # phenotype_scoring silently falls back to Jaccard similarity when pyhpo is missing, which
    # produces different scores from the same inputs. Record which backend was used, so installing
    # or removing pyhpo invalidates the cached database instead of leaving stale scores in place.
    return {"inputs": {label: hash_file(path) for label, path in inputs.items()},
            "flags": flags,
            "phenotype_backend": "pyhpo" if _has_module("pyhpo") else "jaccard",
            "trails_sha": trails_code_version()}


def needs_rebuild(db_path, manifest_path, manifest, force):
    """Return True if the DB must be (re)built: forced, missing, or inputs/flags changed."""
    if force or not os.path.exists(db_path) or not os.path.exists(manifest_path):
        return True
    try:
        with open(manifest_path) as f:
            return json.load(f) != manifest
    except (json.JSONDecodeError, OSError):
        return True


def run_build(args, db_path, gene_table, genes_to_phenotype, known_loci_json, strchive_loci_json,
              trexplorer_annotations_database):
    """Build the result DB in-process, directly from the two TSVs (no intermediate files)."""
    import build_database
    print(f"Building {db_path} from {args.repeat_copy_numbers_tsv} + {args.sample_metadata_tsv} ...")
    build_database.build(
        args.repeat_copy_numbers_tsv, args.sample_metadata_tsv, db_path,
        phenotypes_table=args.phenotypes_table, gene_table=gene_table,
        genes_to_phenotype=genes_to_phenotype, known_loci_json=known_loci_json,
        strchive_loci_json=strchive_loci_json, source_label=args.source_label,
        n_loci=args.n_loci, skip_phenotype_scores=args.skip_phenotype_scores,
        population_gate_cohorts=args.population_gate_cohorts,
        trexplorer_annotations_database=trexplorer_annotations_database)


def serve(db_path, args, known_loci_json=None, strchive_loci_json=None):
    """Start results_server.py against the single built DB (forwarding optional reference inputs)."""
    # No --annotations-db: the server keeps the notes and tags next to the database it is given,
    # one note/tag store per result database (results_server.default_annotations_db_path).
    cmd = [sys.executable, os.path.join(SCRIPT_DIR, "results_server.py"), "--db", db_path]
    for flag, value in [("--sample-table", args.sample_metadata_tsv),
                        ("--known-loci-json", known_loci_json),
                        ("--strchive-loci-json", strchive_loci_json),
                        ("--initial-require-above-population",
                         args.initial_require_above_population),
                        ("--port", args.port),
                        ("--host", args.host)]:
        if value:
            cmd += [flag, str(value)]
    print("+ " + " ".join(cmd))
    try:
        returncode = subprocess.run(cmd).returncode
    except KeyboardInterrupt:
        # Ctrl-C, the usual way to stop the server, reaches the server and this launcher alike.
        # The server shuts down on it, so there is nothing to report.
        return
    if returncode != 0:
        # The server has already printed why; a traceback here would only bury that. A server
        # killed by a signal exits with 128 + the signal number, as a shell would report it.
        print(f"ERROR: the results server stopped with exit status {returncode} (see its messages "
              f"above)", file=sys.stderr)
        sys.exit(returncode if returncode > 0 else 128 - returncode)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repeat-copy-numbers-tsv", required=True,
                        help="merged repeat-copy-numbers TSV (one row per locus, one column per sample)")
    parser.add_argument("--sample-metadata-tsv", required=True,
                        help="sample metadata TSV (one row per sample; only sample_id is required)")
    parser.add_argument("--phenotypes-table", help="optional phenotype TSV (HPO terms per sample)")
    parser.add_argument("--gene-table",
                        help="override the gene-disease table (default: genelookup's table, cached in reference_data/)")
    parser.add_argument("--genes-to-phenotype",
                        help="override the HPO genes_to_phenotype.txt (default: the cached reference copy)")
    parser.add_argument("--known-loci-json",
                        help="override the known-loci catalog (default: the cached reference copy)")
    parser.add_argument("--strchive-loci-json",
                        help="override the cached STRchive-loci.json (default: the cached reference copy)")
    parser.add_argument("--source-label", default="",
                        help="value written to the loci 'Source' column (e.g. the genotyping tool)")
    parser.add_argument("--db", help="result DB path (default: derived from the genotype TSV name)")
    parser.add_argument("--work-dir", default=".",
                        help="directory for the result database (default: current dir)")
    parser.add_argument("--reference-data-dir", default=trails_setup.DEFAULT_REFERENCE_DIR,
                        help="reference-data cache (default: reference_data/ in the TRails directory)")
    parser.add_argument("-n", "--n-loci", type=int, help="limit to the first N loci (for testing)")
    parser.add_argument("--skip-phenotype-scores", action="store_true",
                        help="skip phenotype scoring (faster build)")
    parser.add_argument("--population-gate-cohorts",
                        default=",".join(analysis_columns.POPULATION_P99_GATE_COHORTS),
                        help="comma-separated population cohorts whose 99th percentiles an allele must "
                             "exceed to count in the #Aff / #Fam outlier counts and the phenotype scores. "
                             f"Choices: {', '.join(analysis_columns.POPULATION_GATE_COHORT_CHOICES)} "
                             "(default: %(default)s). Changing it rebuilds the database")
    parser.add_argument("--rebuild", action="store_true", help="force a rebuild even if inputs are unchanged")
    parser.add_argument("--no-serve", action="store_true", help="build the DB but do not start the server")
    parser.add_argument("--initial-require-above-population", choices=INITIAL_REQUIRE_ABOVE_POPULATION_CHOICES,
                        help="population filter the page starts with when opened with no saved filters, "
                             "e.g. 'all' for a short-read callset (forwarded to results_server.py)")
    parser.add_argument("--port", help="server port (forwarded to results_server.py)")
    parser.add_argument("--host", help="server host (forwarded to results_server.py)")
    args = parser.parse_args()
    try:
        args.population_gate_cohorts = analysis_columns.parse_population_gate_cohorts(
            args.population_gate_cohorts)
    except ValueError as error:
        parser.error(str(error))
    # Checked before anything slow: otherwise a mistyped path surfaces only as a traceback after
    # the virtual environment, the pip install and the reference download (about 1 GB) have run.
    missing_input_files = [
        f"ERROR: {option}: {path} is not an existing file"
        for option, path in [("--repeat-copy-numbers-tsv", args.repeat_copy_numbers_tsv),
                             ("--sample-metadata-tsv", args.sample_metadata_tsv),
                             ("--phenotypes-table", args.phenotypes_table),
                             ("--gene-table", args.gene_table),
                             ("--genes-to-phenotype", args.genes_to_phenotype),
                             ("--known-loci-json", args.known_loci_json),
                             ("--strchive-loci-json", args.strchive_loci_json)]
        if path and not os.path.isfile(path)]
    if missing_input_files:
        sys.exit("\n".join(missing_input_files))

    os.makedirs(args.work_dir, exist_ok=True)
    db_path = args.db or default_db_path(args.repeat_copy_numbers_tsv, args.work_dir)
    # Created now: the build first writes into this directory at its very end, and the reference
    # download comes before that, so a missing directory would otherwise waste both.
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)

    # 1. dependencies, 2. public reference data.
    ensure_running_in_trails_virtual_environment()
    ensure_dependencies()
    print(f"TRails reference-data setup -> {args.reference_data_dir}")
    reference = trails_setup.fetch_reference_data(args.reference_data_dir)
    genes_to_phenotype = args.genes_to_phenotype or reference.get(trails_setup.HPO_GENES_TO_PHENOTYPE_NAME)
    known_loci_json = args.known_loci_json or reference.get(trails_setup.STR_ANALYSIS_KNOWN_LOCI_CATALOG_NAME)
    strchive_loci_json = args.strchive_loci_json or reference.get(trails_setup.STRCHIVE_LOCI_NAME)
    gene_table = args.gene_table or reference.get(trails_setup.GENELOOKUP_GENE_TABLE_NAME)
    trexplorer_annotations_database = reference.get(trails_setup.TREXPLORER_ANNOTATIONS_DATABASE_NAME)

    # 3. (re)build the DB if stale.
    manifest_path = db_path + ".build_manifest.json"
    manifest = build_manifest(
        inputs={"repeat_copy_numbers": args.repeat_copy_numbers_tsv,
                "sample_metadata": args.sample_metadata_tsv,
                "phenotypes": args.phenotypes_table,
                "gene_table": gene_table,
                "genes_to_phenotype": genes_to_phenotype,
                "known_loci_json": known_loci_json,
                "strchive_loci_json": strchive_loci_json,
                "trexplorer_annotations": trexplorer_annotations_database},
        flags=build_flags(args))

    if needs_rebuild(db_path, manifest_path, manifest, args.rebuild):
        run_build(args, db_path, gene_table, genes_to_phenotype, known_loci_json, strchive_loci_json,
                  trexplorer_annotations_database)
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
        print(f"Built {db_path}.")
    else:
        print(f"Database up to date ({db_path}); skipping rebuild. Use --rebuild to force.")

    # 4. serve.
    if args.no_serve:
        print("--no-serve set; not starting the server.")
        return
    serve(db_path, args, known_loci_json=known_loci_json, strchive_loci_json=strchive_loci_json)


if __name__ == "__main__":
    main()
