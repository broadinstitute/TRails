"""TRails — one command from your two tables to a running server.

    python3 trails.py \
        --repeat-copy-numbers-tsv your_cohort.tsv.gz \
        --sample-metadata-tsv     your_samples.tsv.gz \
        [--phenotypes-table your_phenotypes.tsv.gz]

In order, this:
  1. installs any missing Python dependencies (requirements.txt),
  2. downloads any missing public (class-A) reference data into reference_data/,
  3. (re)builds the SQLite result database if your inputs changed (sha256 build manifest),
  4. starts the local results server (unless --no-serve).

Only --repeat-copy-numbers-tsv and --sample-metadata-tsv are required. The phenotype table
and the licensed (class-B) gene-disease table (--gene-table) are optional; absent ones are
simply skipped. See docs/INPUT_FORMATS.md.

Note: this orchestrates build_database.build() — the ported analysis engine that reads the
TSVs and populates the DB directly, with no intermediate files — and results_server.py (serve).
Both land during bootstrap (see the refactor plan, step 0).
"""

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys

import trails_setup

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# import names of the runtime deps in requirements.txt.
REQUIRED_MODULES = ["flask", "pandas", "numpy", "tqdm", "msgpack", "intervaltree",
                    "pyhpo", "requests"]


def _has_module(name):
    """Return True if an import of name would resolve, without importing it."""
    try:
        return importlib.util.find_spec(name) is not None
    except ModuleNotFoundError:
        return False


def ensure_dependencies():
    """pip-install requirements.txt if any runtime dependency is missing."""
    if all(_has_module(m) for m in REQUIRED_MODULES):
        return
    print("Installing Python dependencies (requirements.txt) ...")
    subprocess.run([sys.executable, "-m", "pip", "install", "-r",
                    os.path.join(SCRIPT_DIR, "requirements.txt")], check=True)
    still_missing = [m for m in REQUIRED_MODULES if not _has_module(m)]
    if still_missing:
        sys.exit(f"ERROR: dependencies still missing after install: {', '.join(still_missing)}")


def hash_file(path):
    """Return the hex sha256 of a file, or None if path is empty or missing."""
    if not path or not os.path.exists(path):
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def trails_git_sha():
    """Return the TRails git commit, or 'unknown' outside a git checkout."""
    try:
        return subprocess.run(["git", "-C", SCRIPT_DIR, "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def default_db_path(tsv_path, work_dir):
    """Return the default result-DB path derived from the genotype TSV's basename in work_dir."""
    base = os.path.basename(tsv_path)
    for ext in (".gz", ".tsv", ".txt"):
        if base.endswith(ext):
            base = base[:-len(ext)]
    return os.path.join(work_dir, base + ".with_analysis_columns.db")


def build_manifest(inputs, flags):
    """Return the build manifest: input sha256s + build flags + the TRails commit.

    Args:
        inputs: {label: path} for every file whose contents affect the built DB.
        flags: {name: value} for every build flag that affects the DB.
    """
    return {"inputs": {label: hash_file(path) for label, path in inputs.items()},
            "flags": flags,
            "trails_sha": trails_git_sha()}


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


def serve(db_path, args, known_loci_json=None, strchive_loci_json=None):
    """Start results_server.py against the single built DB (forwarding optional reference inputs)."""
    cmd = [sys.executable, os.path.join(SCRIPT_DIR, "results_server.py"), "--db", db_path]
    for flag, value in [("--sample-table", args.sample_metadata_tsv),
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
