"""Download TRails' public (class-A) reference data — the single source of truth.

Both install.sh and trails.py call this, so the reference-data list lives in exactly one
place and the two can never drift. Besides the HPO and known-disease-locus files, this
fetches TRExplorer's population statistics and locus annotations, and genelookup's gene table
(without OMIM's licensed columns, which genelookup does not export). Your own (class-C) tables
are never fetched here.

CLI, from the TRails directory (the TRExplorer step needs the packages in .venv/):
    .venv/bin/python trails/trails_setup.py [fetch] [--dir reference_data] [--remote] [--force]
"""

import argparse
import contextlib
import gzip
import json
import os
import urllib.parse
import urllib.request

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# reference_data/ sits in the TRails directory, one level above this trails/ code directory
DEFAULT_REFERENCE_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "reference_data")

# Pinned public sources: bump these versions deliberately. The HPO, str-analysis and STRchive URLs
# below resolve to an immutable ref (a release tag or a commit sha), so two installs made at
# different times get the same files; a "main" URL would silently change the annotations.
# TRExplorer's table is pinned only by its release directory (v2.1): the file there can be replaced
# in place, and the bucket keeps no older copy to pin to, so an install made after a re-upload gets
# the new file while an existing install keeps its own until --force (TRAILS_FORCE=1). genelookup's
# gene table is not versioned at all (see below). Those two inputs can differ between two installs.
HPO_RELEASE = "v2025-05-06"   # Human Phenotype Ontology annotation release tag
STR_ANALYSIS_COMMIT = "15c642927de66703dcaf5ec92fb08ee8249dc745"   # broadinstitute/str-analysis
STRCHIVE_COMMIT = "d7961b2c3bbc13cbcb256c19ddbf8779f7e76ca6"       # dashnowlab/STRchive
TREXPLORER_CATALOG_VERSION = "v2.1"   # TRExplorer repeat catalog release (gs://tandem-repeat-catalog)

# The format of the database trexplorer_annotations.build_trexplorer_annotations_database builds
# from TRExplorer's table. Bump it whenever a change to trexplorer_annotations.py changes which
# columns that database keeps or how it builds them. It is part of the source the
# reference-versions manifest records for the database, so a bump makes existing installs rebuild
# their copy instead of reading columns it does not have.
TREXPLORER_ANNOTATIONS_DATABASE_FORMAT_VERSION = 1

TREXPLORER_ANNOTATIONS_TSV_URL = (
    f"https://storage.googleapis.com/tandem-repeat-catalog/{TREXPLORER_CATALOG_VERSION}/"
    f"TRExplorer.repeat_catalog_{TREXPLORER_CATALOG_VERSION}.hg38.1_to_1000bp_motifs.EH.with_annotations.tsv.gz")
TREXPLORER_ANNOTATIONS_DATABASE_NAME = f"TRExplorer_{TREXPLORER_CATALOG_VERSION}.hg38.annotations.duckdb"

# genelookup.broadinstitute.org's gene table, through the query service its pages use. The table is
# not versioned, so it is downloaded once and refreshed only by --force (TRAILS_FORCE=1).
GENELOOKUP_QUERY_URL = "https://us-central1-cmg-analysis.cloudfunctions.net/query_gene_lookup_db"
GENELOOKUP_GENE_TABLE_SQL = (
    "SELECT ensembl_gene_id AS gene_id, gene_symbol, gene_aliases, pLI_v2, pLI_v4, lof_oe_ci_upper_v4, "
    "hgnc_gene_id, inheritance, LLM_phenotype_summary, sources "
    "FROM `cmg-analysis.gene_lookup.combined_gene_disease_association_table`")
GENELOOKUP_GENE_TABLE_NAME = "genelookup_gene_table.tsv.gz"


@contextlib.contextmanager
def _removed_on_failure(path):
    """Remove path if the block raises anything, Ctrl-C included, then re-raise.

    A partial download would otherwise stay behind, and TRExplorer's table alone is about 900 MB,
    which matters most when the download failed because the disk filled up.
    """
    try:
        yield
    except BaseException:
        if os.path.exists(path):
            os.remove(path)
        raise


def _download(url, dest):
    """Download url to dest atomically (write a .tmp file, then rename)."""
    tmp = dest + ".tmp"
    with _removed_on_failure(tmp):
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "TRails-setup"})) as response, \
                open(tmp, "wb") as out:
            for chunk in iter(lambda: response.read(1 << 20), b""):
                out.write(chunk)
        os.replace(tmp, dest)


def _download_trexplorer_annotations(source, dest):
    """Download TRExplorer's annotated catalog table and keep the columns TRails uses in dest.

    The table is about 900 MB; only the much smaller database built from it is kept. source is
    the table's URL plus a "#database_format=N" fragment, which only the reference-versions
    manifest uses; it is dropped before downloading.
    """
    # Imported here so that the rest of this module, and --remote, work before the Python
    # dependencies (duckdb, pandas) are installed.
    import trexplorer_annotations

    url = urllib.parse.urldefrag(source).url
    downloaded_table = dest + ".download.tsv.gz"
    print(f"    downloading {url} (about 900 MB) ...")
    try:
        _download(url, downloaded_table)
        print(f"    keeping the columns TRails uses in {dest} ...")
        trexplorer_annotations.build_trexplorer_annotations_database(downloaded_table, dest)
    finally:
        if os.path.exists(downloaded_table):
            os.remove(downloaded_table)


def _export_genelookup_gene_table(sql, dest):
    """Export genelookup's gene table to dest (a gzipped TSV) through its query service.

    The service answers with links to one or more gzipped TSV files, each with its own header
    line, which are concatenated under a single header.
    """
    request = urllib.request.Request(
        GENELOOKUP_QUERY_URL, data=json.dumps({"sql": sql, "export_to_file_format": "TSV"}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "TRails-setup"})
    with urllib.request.urlopen(request, timeout=600) as response:
        result = json.load(response)
    if result.get("status") != "success":
        raise RuntimeError(f"genelookup export failed: {result.get('error', result)}")
    tmp = dest + ".tmp"
    with _removed_on_failure(tmp):
        with gzip.open(tmp, "wt") as out:
            for file_index, url in enumerate(result["public_urls"]):
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "TRails-setup"})) as part:
                    lines = gzip.decompress(part.read()).decode().splitlines(keepends=True)
                out.writelines(lines if file_index == 0 else lines[1:])
        os.replace(tmp, dest)


# The local file names of the other reference files, which trails_launcher looks them up by.
HPO_GENES_TO_PHENOTYPE_NAME = "genes_to_phenotype.txt"
STR_ANALYSIS_KNOWN_LOCI_CATALOG_NAME = "variant_catalog_without_offtargets.GRCh38.json"
STRCHIVE_LOCI_NAME = "STRchive-loci.json"

# (local filename, source, function that fetches source into the file, the trails.py option that
# takes a local copy of the file instead or None): one class-A file per entry. The source is what
# the reference-versions manifest below records. TRExplorer's database has no trails.py option:
# trails.py reads it only from the reference-data directory.
REFERENCE_FILES = [
    (HPO_GENES_TO_PHENOTYPE_NAME,
     f"https://github.com/obophenotype/human-phenotype-ontology/releases/download/{HPO_RELEASE}/genes_to_phenotype.txt",
     _download, "--genes-to-phenotype"),
    (STR_ANALYSIS_KNOWN_LOCI_CATALOG_NAME,
     f"https://raw.githubusercontent.com/broadinstitute/str-analysis/{STR_ANALYSIS_COMMIT}/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json",
     _download, "--known-loci-json"),
    (STRCHIVE_LOCI_NAME,
     f"https://raw.githubusercontent.com/dashnowlab/STRchive/{STRCHIVE_COMMIT}/data/STRchive-loci.json",
     _download, "--strchive-loci-json"),
    (TREXPLORER_ANNOTATIONS_DATABASE_NAME,
     f"{TREXPLORER_ANNOTATIONS_TSV_URL}#database_format={TREXPLORER_ANNOTATIONS_DATABASE_FORMAT_VERSION}",
     _download_trexplorer_annotations, None),
    (GENELOOKUP_GENE_TABLE_NAME, GENELOOKUP_GENE_TABLE_SQL, _export_genelookup_gene_table, "--gene-table"),
]


# Records which source each cached reference file came from (for TRExplorer's database, the table
# URL plus the database format). Without it, bumping a pin above would leave every existing install
# on its old copy (the filenames do not change), which is exactly the reproducibility the pins are
# meant to provide.
CACHE_MANIFEST_NAME = ".reference_versions.json"


def _read_cache_manifest(dest_dir):
    """Return {filename: source} for the cached reference files, or {} if unknown."""
    try:
        with open(os.path.join(dest_dir, CACHE_MANIFEST_NAME)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def fetch_reference_data(dest_dir=DEFAULT_REFERENCE_DIR, remote=False, force=False):
    """Ensure every class-A reference file is present in dest_dir; return {name: path}.

    Idempotent: a file already present, non-empty, and recorded as having come from the currently
    pinned URL is left as-is unless force=True. A file whose pin has since been bumped is
    re-downloaded, so an existing install picks up a deliberate version change.

    Args:
        dest_dir: directory the reference files are cached in (created if missing).
        remote: if True, don't download — just print the URLs (offline / remote mode).
        force: re-download even when a file is already present.
    """
    os.makedirs(dest_dir, exist_ok=True)
    cached_urls = _read_cache_manifest(dest_dir)
    paths = {}
    for name, source, fetch, trails_option in REFERENCE_FILES:
        dest = os.path.join(dest_dir, name)
        if remote:
            print(f"  [remote] {name} -> {source}")
            paths[name] = dest
        elif (not force and os.path.exists(dest) and os.path.getsize(dest) > 0
                and cached_urls.get(name, source) == source):
            # A cached file is kept when its recorded source matches the current pin, or when no
            # source is recorded for it, in which case the current pin is recorded for it now. Only
            # a recorded source that differs from the pin triggers a re-download. Keeping a file
            # with no record avoids a needless download and keeps the file usable offline.
            if name not in cached_urls:
                cached_urls[name] = source
            print(f"  [skip] {name} already present")
            paths[name] = dest
        else:
            print(f"  [fetch] {name}")
            # A reference-data download failure (e.g. offline / network error) must
            # NOT abort the build: the downstream pipeline tolerates any of these
            # class-A inputs being absent (KnownDiseaseLocus, phenotype scoring, the
            # population statistics and the gene table are simply left out). Warn clearly
            # and leave this reference unmapped.
            try:
                fetch(source, dest)
                paths[name] = dest
                cached_urls[name] = source
            except Exception as e:  # noqa: BLE001 - any network/IO failure is non-fatal here
                if os.path.exists(dest) and os.path.getsize(dest) > 0:
                    print(f"  [WARN] could not fetch {name} ({type(e).__name__}: {e}); "
                          f"using the copy already in {dest_dir}, which may predate the "
                          f"currently pinned version.")
                    paths[name] = dest
                else:
                    local_copy_advice = (f"Pass a local copy with trails.py {trails_option}, or re-run "
                                         if trails_option else "Re-run ")
                    print(f"  [WARN] could not fetch {name} ({type(e).__name__}: {e}); "
                          f"continuing without it. {local_copy_advice}trails_setup.py or "
                          f"install.sh later to retry.")
    if not remote:
        try:
            with open(os.path.join(dest_dir, CACHE_MANIFEST_NAME), "w") as f:
                json.dump(cached_urls, f, indent=2, sort_keys=True)
        except OSError as e:
            print(f"  [WARN] could not record reference versions ({type(e).__name__}: {e}); "
                  f"a future pin bump may not trigger a re-download.")
    return paths


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", nargs="?", default="fetch", choices=["fetch"],
                        help="only 'fetch' is supported (the default)")
    parser.add_argument("--dir", default=DEFAULT_REFERENCE_DIR,
                        help="reference-data directory (default: reference_data/ in the TRails directory)")
    parser.add_argument("--remote", action="store_true",
                        help="don't download — print the reference URLs instead (offline mode)")
    parser.add_argument("--force", action="store_true", help="re-download even if present")
    args = parser.parse_args()
    print(f"TRails reference-data setup -> {args.dir}")
    fetch_reference_data(args.dir, remote=args.remote, force=args.force)
    print("Done.")


if __name__ == "__main__":
    main()
