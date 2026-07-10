"""Download TRails' public (class-A) reference data — the single source of truth.

Both install.sh and trails.py call this, so the reference-data list lives in exactly one
place and the two can never drift. Class-B (licensed) and class-C (your own) inputs are never
fetched here — TRails runs without them; supply your own via CLI flags.

CLI:
    python3 trails_setup.py [fetch] [--dir reference_data] [--remote] [--force]
"""

import argparse
import os
import urllib.request

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_REFERENCE_DIR = os.path.join(SCRIPT_DIR, "reference_data")

# Pinned public sources — bump these versions deliberately.
HPO_RELEASE = "v2025-05-06"   # Human Phenotype Ontology annotation release tag

# (local filename, public URL) — one class-A file per entry.
REFERENCE_FILES = [
    ("genes_to_phenotype.txt",
     f"https://github.com/obophenotype/human-phenotype-ontology/releases/download/{HPO_RELEASE}/genes_to_phenotype.txt"),
    ("variant_catalog_without_offtargets.GRCh38.json",
     "https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
    ("STRchive-loci.json",
     "https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/STRchive-loci.json"),
]


def _download(url, dest):
    """Download url to dest atomically (write a .tmp file, then rename)."""
    tmp = dest + ".tmp"
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "TRails-setup"})) as response, \
            open(tmp, "wb") as out:
        for chunk in iter(lambda: response.read(1 << 20), b""):
            out.write(chunk)
    os.replace(tmp, dest)


def fetch_reference_data(dest_dir=DEFAULT_REFERENCE_DIR, remote=False, force=False):
    """Ensure every class-A reference file is present in dest_dir; return {name: path}.

    Idempotent: a file already present (and non-empty) is left as-is unless force=True.

    Args:
        dest_dir: directory the reference files are cached in (created if missing).
        remote: if True, don't download — just print the URLs (offline / remote mode).
        force: re-download even when a file is already present.
    """
    os.makedirs(dest_dir, exist_ok=True)
    paths = {}
    for name, url in REFERENCE_FILES:
        dest = os.path.join(dest_dir, name)
        if remote:
            print(f"  [remote] {name} -> {url}")
            paths[name] = dest
        elif not force and os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"  [skip] {name} already present")
            paths[name] = dest
        else:
            print(f"  [fetch] {name}")
            # A reference-data download failure (e.g. offline / network error) must
            # NOT abort the build: the downstream pipeline tolerates any of these
            # class-A inputs being absent (KnownDiseaseLocus / phenotype scoring are
            # simply skipped). Warn clearly and leave this reference unmapped.
            try:
                _download(url, dest)
                paths[name] = dest
            except Exception as e:  # noqa: BLE001 - any network/IO failure is non-fatal here
                print(f"  [WARN] could not fetch {name} ({type(e).__name__}: {e}); "
                      f"continuing without it. Pass a local copy on the trails.py "
                      f"command line, or re-run trails_setup.py later to retry.")
    return paths


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", nargs="?", default="fetch", choices=["fetch"],
                        help="only 'fetch' is supported (the default)")
    parser.add_argument("--dir", default=DEFAULT_REFERENCE_DIR,
                        help="reference-data directory (default: reference_data/ next to this script)")
    parser.add_argument("--remote", action="store_true",
                        help="don't download — print the reference URLs instead (offline mode)")
    parser.add_argument("--force", action="store_true", help="re-download even if present")
    args = parser.parse_args()
    print(f"TRails reference-data setup -> {args.dir}")
    fetch_reference_data(args.dir, remote=args.remote, force=args.force)
    print("Done.")


if __name__ == "__main__":
    main()
