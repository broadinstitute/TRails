#!/usr/bin/env bash
#
# install.sh — one-line installer / updater for TRails.
#
# Usage:
#   curl -fsSL https://raw.githubusercontent.com/broadinstitute/TRails/main/install.sh | bash
#   (or, from a checkout: ./install.sh)
#
# What it does, in order:
#   1. If not run from a checkout, downloads + extracts TRails into ./TRails (curl + tar,
#      no git). Re-running seamlessly UPDATES to the latest $TRAILS_VERSION; an up-to-date
#      copy is reused. Resumable (curl -C -) and idempotent.
#   2. Installs Python deps (requirements.txt).
#   3. Downloads PUBLIC (class-A) reference data into reference_data/ (skips files already
#      present; resumable). Licensed/controlled (class B) and your own (class C) data are
#      NEVER fetched — TRails runs without them; supply your own via CLI flags.
#      See docs/INPUT_FORMATS.md.
#
# Env knobs: TRAILS_REPO=owner/name  TRAILS_VERSION=<branch|tag|commit>  TRAILS_INSTALL_DIR=./TRails
#            TRAILS_FORCE=1 (force re-download)  TRAILS_REMOTE=1 (print ref-data URLs, don't download)

set -euo pipefail

# Overridable via env: TRAILS_REPO=owner/name TRAILS_VERSION=main TRAILS_INSTALL_DIR=./TRails
TRAILS_REPO="${TRAILS_REPO:-broadinstitute/TRails}"
TRAILS_VERSION="${TRAILS_VERSION:-main}"
INSTALL_DIR="${TRAILS_INSTALL_DIR:-./TRails}"

# Install directory may also be passed on the command line (takes precedence over the env
# var). When piping from curl, pass args after `--`, e.g.:
#   curl -fsSL .../install.sh | bash -s -- --dir /opt/TRails
#   curl -fsSL .../install.sh | bash -s -- /opt/TRails
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d|--dir)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a directory argument" >&2; exit 1; }
      INSTALL_DIR="$2"; INSTALL_DIR_EXPLICIT=1; shift 2 ;;
    --dir=*)   INSTALL_DIR="${1#*=}"; INSTALL_DIR_EXPLICIT=1; shift ;;
    -h|--help) echo "Usage: install.sh [--dir INSTALL_DIR]"; exit 0 ;;
    -*)        echo "ERROR: unknown option: $1" >&2; exit 1 ;;
    *)         INSTALL_DIR="$1"; INSTALL_DIR_EXPLICIT=1; shift ;;   # bare positional = install dir
  esac
done

# Detect mode: running from inside a checkout, or piped from `curl ... | bash`.
TRAILS_DIR=""
if [[ "${BASH_SOURCE[0]:-}" == *install.sh ]]; then
  _d="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  [[ -f "$_d/requirements.txt" ]] && TRAILS_DIR="$_d"
fi

# In checkout mode the install operates on the checkout itself, so an explicit
# install directory cannot be honored; reject it rather than silently ignoring it.
if [[ -n "$TRAILS_DIR" && "${INSTALL_DIR_EXPLICIT:-}" == "1" ]]; then
  echo "ERROR: --dir is not supported when running from a checkout (it installs into the checkout at $TRAILS_DIR)." >&2
  echo "       To install into a different directory, run the curl bootstrap instead:" >&2
  echo "         curl -fsSL https://raw.githubusercontent.com/$TRAILS_REPO/$TRAILS_VERSION/install.sh | bash -s -- --dir $INSTALL_DIR" >&2
  exit 1
fi

# Bootstrap mode: no local checkout -> download + extract the repo. Needs only curl + tar
# + python3 (no git, no jq). Idempotent, resumable, and self-updating:
#   - an existing up-to-date copy is reused (no work);
#   - if the remote $TRAILS_VERSION has advanced, it is re-downloaded seamlessly;
#   - a partial tarball download is resumed via curl -C -;
#   - TRAILS_FORCE=1 forces a fresh re-download.
# The installed commit is recorded in $INSTALL_DIR/.trails_version.
if [[ -z "$TRAILS_DIR" ]]; then
  command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is required" >&2; exit 1; }
  command -v curl    >/dev/null 2>&1 || { echo "ERROR: curl is required"    >&2; exit 1; }
  command -v tar     >/dev/null 2>&1 || { echo "ERROR: tar is required"     >&2; exit 1; }

  VERSION_FILE="$INSTALL_DIR/.trails_version"
  installed_sha=""
  [[ -f "$VERSION_FILE" ]] && installed_sha="$(cat "$VERSION_FILE" 2>/dev/null || true)"

  # Best-effort: resolve the remote commit SHA for $TRAILS_VERSION via the GitHub API.
  # Empty if offline / rate-limited / ref not found -> we fall back gracefully.
  remote_sha="$(curl -fsSL "https://api.github.com/repos/${TRAILS_REPO}/commits/${TRAILS_VERSION}" 2>/dev/null \
      | python3 -c "import sys,json
try: print(json.load(sys.stdin).get('sha',''))
except Exception: pass" 2>/dev/null || true)"

  need_download=0
  if [[ "${TRAILS_FORCE:-0}" == "1" ]]; then
    need_download=1
  elif [[ ! -f "$INSTALL_DIR/requirements.txt" ]]; then
    need_download=1   # nothing installed yet
  elif [[ -n "$remote_sha" && "$remote_sha" != "$installed_sha" ]]; then
    need_download=1
    echo "Update available for $TRAILS_VERSION (${installed_sha:-none} -> ${remote_sha}); updating ..."
  fi

  if [[ "$need_download" == "1" ]]; then
    mkdir -p "$INSTALL_DIR"
    ref_for_tarball="${remote_sha:-$TRAILS_VERSION}"
    if [[ -n "$remote_sha" ]]; then
      tarball_url="https://github.com/${TRAILS_REPO}/archive/${remote_sha}.tar.gz"
    else
      tarball_url="https://github.com/${TRAILS_REPO}/archive/refs/heads/${TRAILS_VERSION}.tar.gz"
    fi
    echo "Downloading TRails ($TRAILS_REPO@${ref_for_tarball}) -> $INSTALL_DIR ..."
    # Per-version tarball name so curl -C - only resumes the SAME version's partial.
    _tarball="$INSTALL_DIR/.trails_src.${ref_for_tarball}.tar.gz"
    curl -fL -C - --retry 3 -o "$_tarball" "$tarball_url"
    tar -xz -f "$_tarball" -C "$INSTALL_DIR" --strip-components=1
    rm -f "$INSTALL_DIR"/.trails_src.*.tar.gz   # clean up only after a successful extract
    printf '%s\n' "${remote_sha:-$TRAILS_VERSION}" > "$VERSION_FILE"
  else
    echo "TRails is up to date at $INSTALL_DIR (${installed_sha:-$TRAILS_VERSION})."
  fi
  TRAILS_DIR="$INSTALL_DIR"
  echo ""
fi

cd "$TRAILS_DIR"
TRAILS_DIR="$(pwd)"   # resolve to absolute so a relative INSTALL_DIR (e.g. ./TRails) doesn't double-nest below
REF_DIR="$TRAILS_DIR/reference_data"
mkdir -p "$REF_DIR"

# --- Python dependencies ------------------------------------------------------------
echo "Installing Python dependencies..."
python3 -m pip install -r "$TRAILS_DIR/requirements.txt"   # Flask, pandas, numpy, str-analysis, pyhpo, ...
echo ""

# --- Public (class-A) reference data ------------------------------------------------
# The class-A reference list + fetch logic live in trails_setup.py (one source of truth,
# shared with trails.py). TRAILS_REMOTE=1 prints URLs instead of downloading;
# TRAILS_FORCE=1 re-downloads.
REMOTE_FLAG=""; [[ "${TRAILS_REMOTE:-0}" == "1" ]] && REMOTE_FLAG="--remote"
FORCE_FLAG="";  [[ "${TRAILS_FORCE:-0}"  == "1" ]] && FORCE_FLAG="--force"
echo "Class A (public, fetched):"
python3 "$TRAILS_DIR/trails_setup.py" fetch --dir "$REF_DIR" $REMOTE_FLAG $FORCE_FLAG

echo ""
echo "Class B (licensed / controlled-access — NOT fetched; supply your own if licensed):"
echo "  - gene-disease table  -> pass via --gene-table PATH        (OMIM/HPO-derived)"
echo ""
echo "Class C (your own data — NOT fetched; see docs/INPUT_FORMATS.md):"
echo "  - repeat-copy-numbers TSV, sample-metadata TSV, optional phenotype TSV"
echo ""
echo "Done. Reference data in $REF_DIR"
echo ""
echo "Next:"
echo "  python3 $TRAILS_DIR/trails.py \\"
echo "      --repeat-copy-numbers-tsv <your_cohort.tsv.gz> \\"
echo "      --sample-metadata-tsv     <your_samples.tsv.gz>"
