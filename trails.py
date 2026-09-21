"""TRails — one command from your two tables to a running server.

    python3 trails.py \
        --repeat-copy-numbers-tsv your_cohort.tsv.gz \
        --sample-metadata-tsv     your_samples.tsv.gz \
        [--phenotypes-table your_phenotypes.tsv.gz]

In order, this:
  1. installs any missing Python dependencies (requirements.txt),
  2. downloads any missing public (class-A) reference data into reference_data/,
  3. (re)builds the DuckDB result database if your inputs changed (sha256 build manifest),
  4. starts the local results server (unless --no-serve).

Only --repeat-copy-numbers-tsv and --sample-metadata-tsv are required. The phenotype table
and the licensed (class-B) gene-disease table (--gene-table) are optional; absent ones are
simply skipped. See docs/INPUT_FORMATS.md.

Note: this orchestrates build_database.build() — the ported analysis engine that reads the
TSVs and populates the DB directly, with no intermediate files — and results_server.py (serve).
"""

import argparse
import ast
import hashlib
import importlib.util
import json
import os
import subprocess
import sys

import trails_setup

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Where the build's import graph starts. The code fingerprint in the build manifest covers
# these files and every TRails module they import, directly or indirectly, so a new pipeline
# module is picked up automatically as soon as an existing one imports it. Add an entry here
# only for a new way of starting a build that does not go through one of these files.
BUILD_ENTRY_POINTS = ["trails.py", "build_database.py"]

# import names of the runtime deps in requirements.txt that the code actually imports.
REQUIRED_MODULES = ["flask", "pandas", "numpy", "intervaltree", "duckdb"]

# Neither of these is needed on any default code path, so a missing one is a warning, not an exit:
#   pyhpo    - ontology-aware phenotype similarity; phenotype_scoring falls back to Jaccard.
#   requests - imported only inside locus_annotations.fetch_strchive_loci, the opt-in network
#              fetch; trails_setup downloads reference data with urllib instead.
OPTIONAL_MODULES = ["pyhpo", "requests"]

# Written after the one optional-dependency install attempt, so a package that will not build here
# does not make every later run re-try (and re-wait for) the same failing pip invocation.
OPTIONAL_INSTALL_MARKER = os.path.join(SCRIPT_DIR, ".trails_optional_deps_attempted")


def _has_module(name):
    """Return True if an import of name would resolve, without importing it."""
    try:
        return importlib.util.find_spec(name) is not None
    except ModuleNotFoundError:
        return False


def ensure_dependencies():
    """pip-install the runtime dependencies if any are missing.

    The required set is installed with check=True: without it nothing can run. The optional set
    lives in its own requirements file and is installed with check=False, because pip fails the
    whole invocation if any listed package cannot be built, which would turn a soft dependency
    into a hard one.
    """
    if not all(_has_module(m) for m in REQUIRED_MODULES):
        print("Installing Python dependencies (requirements.txt) ...")
        subprocess.run([sys.executable, "-m", "pip", "install", "-r",
                        os.path.join(SCRIPT_DIR, "requirements.txt")], check=True)
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
              f"Phenotype scoring falls back to Jaccard similarity without pyhpo; the opt-in "
              f"STRchive fetch needs requests. Delete {OPTIONAL_INSTALL_MARKER} to retry.")


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

    Both "import x" and "from x import y" count, at any nesting level, because the build
    imports some modules inside a function rather than at the top of the file (this file
    imports build_database that way). Relative imports are ignored: TRails is a flat
    directory of modules with no package.

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
    on its own and nobody has to maintain a list. Files the build never imports
    (results_server.py, which is launched as a subprocess, and the *_tests.py files) stay
    out, so editing them does not invalidate the manifest.

    Returns:
        A sorted list of file names (not paths), all of which exist in SCRIPT_DIR.
    """
    found = set()
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


def build_manifest(inputs, flags):
    """Return the build manifest: input sha256s + build flags + the TRails commit.

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


def run_build(args, db_path, genes_to_phenotype, known_loci_json, strchive_loci_json):
    """Build the result DB in-process, directly from the two TSVs (no intermediate files)."""
    import build_database
    print(f"Building {db_path} from {args.repeat_copy_numbers_tsv} + {args.sample_metadata_tsv} ...")
    build_database.build(
        args.repeat_copy_numbers_tsv, args.sample_metadata_tsv, db_path,
        phenotypes_table=args.phenotypes_table, gene_table=args.gene_table,
        genes_to_phenotype=genes_to_phenotype, known_loci_json=known_loci_json,
        strchive_loci_json=strchive_loci_json, source_label=args.source_label,
        n_loci=args.n_loci, skip_phenotype_scores=args.skip_phenotype_scores)
    if not os.path.exists(db_path):
        sys.exit(f"ERROR: build_database did not produce {db_path}")


def legacy_annotations_db_path(db_path):
    """Return the SQLite annotations path a pre-DuckDB TRails would have written.

    Notes and tags lived in "<result database>.annotations.db" next to a result database
    whose own name ended in ".db". The result database now ends in ".duckdb", so swap the
    extension back to name the file the previous version left behind.
    """
    base = os.path.abspath(db_path)
    if base.endswith(".duckdb"):
        base = base[:-len(".duckdb")] + ".db"
    return base + ".annotations.db"


def migrate_legacy_annotations(annotations_db, legacy_path):
    """Copy notes and tags out of a pre-DuckDB SQLite annotations database, once.

    Notes and tags are typed in by hand in the results server, so they are the one thing in
    TRails a rebuild cannot reproduce. This runs only when the DuckDB annotations database
    does not exist yet, and it reads the SQLite file without writing to it, so the old copy
    survives the upgrade either way.

    Args:
        annotations_db: The DuckDB annotations database the server is about to open.
        legacy_path: The SQLite annotations database a previous version would have written.

    Returns:
        A ``(num_notes, num_tags)`` tuple; ``(0, 0)`` when there was nothing to migrate.
    """
    if os.path.exists(annotations_db) or not os.path.exists(legacy_path):
        return 0, 0

    # The only sqlite3 use left in TRails: the file being read here IS a SQLite database,
    # and reading it is this function's whole job. Imported here rather than at the top of
    # the module so it stays visibly local to the migration and nobody takes it for a sign
    # that TRails still runs on SQLite.
    import sqlite3
    import urllib.parse

    # The path goes into the URI percent-escaped: SQLite reads an unescaped "?" as the start
    # of the query parameters and an unescaped "#" as the start of the fragment, so a path
    # containing either would be truncated, and SQLite would quietly open (and create) a
    # different file than the one holding the user's notes.
    uri = "file:" + urllib.parse.quote(os.path.abspath(legacy_path)) + "?mode=ro"
    try:
        legacy = sqlite3.connect(uri, uri=True)
    except sqlite3.Error as error:
        print(f"NOTE: could not open the previous annotations database {legacy_path}: {error}")
        return 0, 0
    notes = []
    tags = []
    try:
        # Each table is read on its own: an unreadable or missing tags table must not cost
        # the user their notes, and a file from an older schema may have only one of them.
        for sql, destination in (
                ("SELECT locus_id, note_text, created_at, updated_at FROM notes", notes),
                ("SELECT locus_id, tag, created_at FROM tags", tags)):
            try:
                destination.extend(legacy.execute(sql).fetchall())
            except sqlite3.Error as error:
                print(f"NOTE: skipping part of {legacy_path}: {error}")
    finally:
        legacy.close()

    if not notes and not tags:
        return 0, 0

    import duckdb_compat

    # Built at a temporary path and moved into place, so an interrupted migration cannot
    # leave a half-copied database behind that the next run would take for a finished one.
    tmp_path = annotations_db + ".migration.tmp"
    for stale_path in (tmp_path, tmp_path + ".wal"):
        if os.path.exists(stale_path):
            os.remove(stale_path)
    try:
        connection = duckdb_compat.connect(tmp_path)
        try:
            # The same schema results_server.load_annotations creates.
            connection.execute(
                "CREATE TABLE notes (locus_id TEXT PRIMARY KEY, note_text TEXT NOT NULL, "
                "created_at TEXT NOT NULL, updated_at TEXT NOT NULL)")
            connection.execute(
                "CREATE TABLE tags (locus_id TEXT NOT NULL, tag TEXT NOT NULL, "
                "created_at TEXT NOT NULL, PRIMARY KEY (locus_id, tag))")
            if notes:
                connection.executemany("INSERT INTO notes VALUES (?, ?, ?, ?)", notes)
            if tags:
                connection.executemany("INSERT INTO tags VALUES (?, ?, ?)", tags)
            connection.commit()
        finally:
            connection.close()
        os.replace(tmp_path, annotations_db)
    except (duckdb_compat.Error, OSError) as error:
        for leftover in (tmp_path, tmp_path + ".wal"):
            if os.path.exists(leftover):
                os.remove(leftover)
        print(f"NOTE: could not migrate the notes and tags in {legacy_path}: {error}. They are "
              f"still in that file; delete {annotations_db} and re-run TRails to try again.")
        return 0, 0

    print(f"Migrated {len(notes):,} note(s) and {len(tags):,} tag(s) from {legacy_path} "
          f"into {annotations_db}. The original file is left untouched.")
    return len(notes), len(tags)


def serve(db_path, args, known_loci_json=None, strchive_loci_json=None):
    """Start results_server.py against the single built DB (forwarding optional reference inputs)."""
    cmd = [sys.executable, os.path.join(SCRIPT_DIR, "results_server.py"), "--db", db_path]
    # The server defaults --annotations-db to a bare relative name, which is resolved against the
    # server's cwd, so notes would follow the launch directory. Derive it from the full database
    # filename instead: anchoring on the directory alone would make every result database in that
    # directory share one note/tag store, and notes are keyed by LocusId only.
    annotations_db = os.path.abspath(db_path) + ".annotations.duckdb"
    # An upgrade from the SQLite version renamed the result database, and with it this
    # path. Bring the previous notes and tags across before the server creates an empty one.
    migrate_legacy_annotations(annotations_db, legacy_annotations_db_path(db_path))
    for flag, value in [("--sample-table", args.sample_metadata_tsv),
                        ("--annotations-db", annotations_db),
                        ("--known-loci-json", known_loci_json),
                        ("--strchive-loci-json", strchive_loci_json),
                        ("--port", args.port),
                        ("--host", args.host)]:
        if value:
            cmd += [flag, str(value)]
    print("+ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repeat-copy-numbers-tsv", required=True,
                        help="merged repeat-copy-numbers TSV (one row per locus, one column per sample)")
    parser.add_argument("--sample-metadata-tsv", required=True,
                        help="sample metadata TSV (one row per sample; only sample_id is required)")
    parser.add_argument("--phenotypes-table", help="optional phenotype TSV (HPO terms per sample)")
    parser.add_argument("--gene-table", help="optional gene-disease table (class B; license-gated)")
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
    parser.add_argument("--reference-data-dir", default=os.path.join(SCRIPT_DIR, "reference_data"),
                        help="reference-data cache (default: reference_data/ next to this script)")
    parser.add_argument("-n", "--n-loci", type=int, help="limit to the first N loci (for testing)")
    parser.add_argument("--skip-phenotype-scores", action="store_true",
                        help="skip phenotype scoring (faster build)")
    parser.add_argument("--rebuild", action="store_true", help="force a rebuild even if inputs are unchanged")
    parser.add_argument("--no-serve", action="store_true", help="build the DB but do not start the server")
    parser.add_argument("--port", help="server port (forwarded to results_server.py)")
    parser.add_argument("--host", help="server host (forwarded to results_server.py)")
    args = parser.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)

    # 1. dependencies, 2. public reference data.
    ensure_dependencies()
    print(f"TRails reference-data setup -> {args.reference_data_dir}")
    reference = trails_setup.fetch_reference_data(args.reference_data_dir)
    genes_to_phenotype = args.genes_to_phenotype or reference.get("genes_to_phenotype.txt")
    known_loci_json = args.known_loci_json or reference.get("variant_catalog_without_offtargets.GRCh38.json")
    strchive_loci_json = args.strchive_loci_json or reference.get("STRchive-loci.json")

    # 3. (re)build the DB if stale.
    db_path = args.db or default_db_path(args.repeat_copy_numbers_tsv, args.work_dir)
    manifest_path = db_path + ".build_manifest.json"
    manifest = build_manifest(
        inputs={"repeat_copy_numbers": args.repeat_copy_numbers_tsv,
                "sample_metadata": args.sample_metadata_tsv,
                "phenotypes": args.phenotypes_table,
                "gene_table": args.gene_table,
                "genes_to_phenotype": genes_to_phenotype,
                "known_loci_json": known_loci_json,
                "strchive_loci_json": strchive_loci_json},
        flags={"n_loci": args.n_loci, "skip_phenotype_scores": args.skip_phenotype_scores,
               "source_label": args.source_label})

    if needs_rebuild(db_path, manifest_path, manifest, args.rebuild):
        run_build(args, db_path, genes_to_phenotype, known_loci_json, strchive_loci_json)
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
