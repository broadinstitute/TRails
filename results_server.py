#!/usr/bin/env python3
"""Flask-based read-only API server for TRails tandem-repeat outlier results.

Serves data from a single pre-computed SQLite database of locus-level outlier
statistics, enriched with sample metadata, gene-disease associations, and known
disease loci. This is the single-database TRails server: there is no source
selector, no read-visualization (readviz) integration, and no cloud dependency.

The server tolerates optional tables being absent. The ``loci`` table is the only
hard requirement (see ``REQUIRED_COLUMNS``); ``swim_plot``, the per-outlier-type
skinny tables (``sk_*``), the phenotype-score tables, and the Mendelian-violation
tables are all optional, and the corresponding pages/endpoints degrade gracefully
when they are missing.
"""

import argparse
import collections
from datetime import datetime
import json
import math
import os
import sqlite3
import sys
import traceback
import zlib

import msgpack
import flask
from flask import Flask, request, Response
from werkzeug.exceptions import HTTPException
import intervaltree
import jinja2
import numpy as np

# Standalone motif primitive (no str_analysis dependency, no sys.path insertion).
from motif_utilities import compute_canonical_motif

# Known-disease-locus matching and affected-status normalization live in the ported
# locus_annotations module. Import them defensively so the server stays importable
# even when that optional module is not present alongside it; the locus-detail page's
# known-disease annotations simply degrade to "no match" in that case.
try:
    from locus_annotations import (
        compute_jaccard,
        load_known_disease_loci,
        motifs_match,
        normalize_affected_status_for_logic,
    )
    HAVE_LOCUS_ANNOTATIONS = True
except ImportError:
    HAVE_LOCUS_ANNOTATIONS = False

    def normalize_affected_status_for_logic(value):
        """Lowercase + collapse 'possibly affected' to 'affected' (fallback)."""
        if value is None:
            return None
        if isinstance(value, float) and value != value:  # NaN
            return None
        normalized = str(value).strip().lower()
        if normalized == "nan":
            return None
        if normalized == "possibly affected":
            return "affected"
        return normalized

    def compute_jaccard(start1, end1, start2, end2):
        """Jaccard index for two half-open intervals (fallback)."""
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        if overlap_start >= overlap_end:
            return 0.0
        overlap_size = overlap_end - overlap_start
        union_size = (end1 - start1) + (end2 - start2) - overlap_size
        return overlap_size / union_size if union_size > 0 else 0.0

    def motifs_match(motif1, motif2):
        """Canonical match for <=6bp, length match for longer (fallback)."""
        if not motif1 or not motif2:
            return False
        if len(motif1) <= 6:
            return (compute_canonical_motif(motif1, include_reverse_complement=True)
                    == compute_canonical_motif(motif2, include_reverse_complement=True))
        return len(motif1) == len(motif2)

    def load_known_disease_loci(filepath, fetch_strchive=False, strchive_filepath=None,
                                build_locus_lookup=False):
        """Disabled fallback: no disease catalog available -> empty trees/lookup.

        Mirrors the real loader's (interval_trees, strchive_trees, locus_lookup)
        return shape so the startup call site unpacks identically.
        """
        return {}, {}, {}

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTLIER_TYPE_MAP = {"all": "AllAlleles", "short": "ShortAlleles", "hemi": "HemizygousAlleles"}

# Labels for export filenames (build_export_filename). Kept in sync with the "labels" advertised
# for outlier_type via /api/v1/schema ("Long Allele"/"Biallelic"/"Hemizygous").
FILENAME_OUTLIER_TYPE_LABELS = {"all": "long_allele", "short": "biallelic", "hemi": "hemizygous"}

# Valid sort_by values accepted by /api/v1/loci. Kept in sync with build_api_order_by's
# SORT_MAPPING and advertised verbatim via /api/v1/schema.
VALID_SORTS = (
    "region", "size", "count", "family_count",
    "pairwise_similarity", "gene_phenotype",
    "sigma_hprc_rank", "sigma_aou_rank",
)

REQUIRED_TABLE = "loci"
REQUIRED_COLUMNS = {"LocusId", "Chrom", "Start0Based", "End1Based", "Motif", "MotifSize"}

GENE_REGION_MAP = {
    "cds": ["CDS"],
    "promoter": ["promoter"],
    "utr": ["5' UTR", "3' UTR"],
    "intron": ["intron"],
    "exon": ["exon"],
    "intergenic": ["intergenic"],
}

# Columns to include in list response (static — not outlier-type-specific).
LIST_COLUMNS_STATIC = [
    "LocusId", "Chrom", "Start0Based", "End1Based",
    "Motif", "CanonicalMotif", "MotifSize", "NumRepeatsInReference",
    "Source", "ReferenceRegion",
    "gene_id", "gene_region", "gene_region_rank",
    "GeneTableGeneSymbol", "GeneTableInheritance", "GeneTableLLMPhenotypeSummary",
    "pLI", "inheritance",
    "IsKnownMotif", "IsInMendelianGene", "KnownDiseaseLocus",
    "AoU1027_MaxAllele", "AoU1027_99thPercentile",
    "AoU1027_Stdev", "AoU1027_StdevRankByMotif", "AoU1027_StdevRankTotalNumberByMotif",
    "HPRC256_MaxAllele", "HPRC256_99thPercentile",
    "HPRC256_Stdev", "HPRC256_StdevRankByMotif", "HPRC256_StdevRankTotalNumberByMotif",
    "HPRC256_StdevPercentile",
    "AoU1027_StdevPercentile",
    "TenK10K_MaxAllele", "TenK10K_99thPercentile",
    "TRExplorerSource", "TRExplorerReferenceRepeatPurity",
    "NonCodingAnnotations",
]

# Outlier-type-specific columns (the _{ot} suffix is stripped in the response).
LIST_COLUMNS_OT_SPECIFIC = [
    "FirstAffectedAlleleSize",
    "SecondAffectedAlleleSize",
    "ThirdAffectedAlleleSize",
    "FirstUnaffectedAlleleSize",
    "NumAffectedUnsolvedSamplesAboveUnaffected",
    "FirstAffectedSampleId",
    "FirstAffectedPhenotype",
    "MaxGenePhenoSim",
    "SumPairwiseSim",
]

# ---------------------------------------------------------------------------
# Flask app
# ---------------------------------------------------------------------------

app = Flask(__name__)


def sanitize_for_json(obj):
    """Recursively convert NaN/Infinity to None for JSON compliance."""
    if isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [sanitize_for_json(v) for v in obj]
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
    return obj


def resolve_tag_to_loci(tag):
    """Return the set of LocusIds matching a user-defined tag.

    Only user tags are supported in the single-DB server; the system tags
    (readviz / readviz pending) of the original multi-source server are gone.
    """
    return app.config["ANNOTATIONS"]["tag_to_loci"].get(tag, set())


def parse_motif_size_filter(motif_size_str):
    """Parse a motif_size filter string into SQL clauses and parameters.

    Supports comma-separated entries in these formats:
    - <int>: exact match (e.g., "3" means MotifSize = 3)
    - <int>-: open-ended range (e.g., "5-" means MotifSize >= 5)
    - -<int>: open-ended range (e.g., "-10" means MotifSize <= 10)
    - <int>-<int>: closed range (e.g., "3-6" means MotifSize BETWEEN 3 AND 6)

    Multiple entries are combined with OR.

    Args:
        motif_size_str: The filter string from the API parameter.

    Returns:
        Tuple of (sql_clause, params, error_message). sql_clause is a WHERE clause
        fragment, params is the list of parameter values, and error_message is a
        string when parsing failed (None otherwise).
    """
    if not motif_size_str or not motif_size_str.strip():
        return None, [], None

    parts = [p.strip() for p in motif_size_str.split(",") if p.strip()]
    if not parts:
        return None, [], None

    clauses = []
    params = []

    for part in parts:
        if "-" in part:
            idx = part.index("-")
            left = part[:idx].strip()
            right = part[idx + 1:].strip()
            if left and right:
                try:
                    clauses.append("MotifSize BETWEEN ? AND ?")
                    params.extend([int(left), int(right)])
                except ValueError:
                    return None, [], f"Invalid motif_size range: '{part}'"
            elif left and not right:
                try:
                    params.append(int(left))
                    clauses.append("MotifSize >= ?")
                except ValueError:
                    return None, [], f"Invalid motif_size value: '{part}'"
            elif not left and right:
                try:
                    params.append(int(right))
                    clauses.append("MotifSize <= ?")
                except ValueError:
                    return None, [], f"Invalid motif_size value: '{part}'"
            else:
                return None, [], f"Invalid motif_size format: '{part}'"
        else:
            try:
                params.append(int(part))
                clauses.append("MotifSize = ?")
            except ValueError:
                return None, [], f"Invalid motif_size value: '{part}'"

    if not clauses:
        return None, [], None

    return "(" + " OR ".join(clauses) + ")", params, None


# Override Flask's jsonify to handle NaN values.
_original_jsonify = flask.jsonify


def jsonify(*args, **kwargs):
    """Custom jsonify that converts NaN/Infinity to null."""
    if args:
        data = sanitize_for_json(args[0]) if len(args) == 1 else sanitize_for_json(args)
        return _original_jsonify(data, **kwargs)
    return _original_jsonify(**{k: sanitize_for_json(v) for k, v in kwargs.items()})


TEMPLATE_DIR = os.path.dirname(os.path.abspath(__file__))
jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader(TEMPLATE_DIR))


def _safe_tojson(value):
    """HTML-safe JSON for embedding values in <script> blocks.

    Escapes the characters that would otherwise let a value break out of the
    surrounding <script> tag or be reinterpreted as HTML. The escaped output is
    still valid JSON and parses identically in the browser.
    """
    return (
        json.dumps(value)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace(" ", "\\u2028")
        .replace(" ", "\\u2029")
    )


jinja2_env.filters["tojson"] = _safe_tojson

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description="Flask API server for TRails tandem-repeat outlier results.",
    )
    parser.add_argument(
        "--db",
        required=True,
        help="Path to the single TRails SQLite results database.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=5050,
        help="Server port (default: 5050).",
    )
    parser.add_argument(
        "--host",
        type=str,
        default="127.0.0.1",
        help="Bind address (default: 127.0.0.1).",
    )
    parser.add_argument(
        "--sample-table",
        default=None,
        help="Optional path to sample metadata table for outlier-sample enrichment.",
    )
    parser.add_argument(
        "--known-loci-json",
        default=None,
        help="Optional path to known disease loci JSON (variant catalog).",
    )
    parser.add_argument(
        "--strchive-loci-json",
        default=None,
        help="Optional path to a cached STRchive-loci.json (adds the detail-page STRchive fallback).",
    )
    parser.add_argument(
        "--annotations-db",
        default="annotations.db",
        help="Path to annotations SQLite database for notes/tags (default: annotations.db).",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run Flask in debug mode (auto-reloader + interactive debugger). Off by default; "
             "enabling it re-runs the full startup scan in the reloader child and exposes the debugger.",
    )
    return parser.parse_args()


def validate_database(path):
    """Validate that the database file exists and has the required table/columns.

    Args:
        path: Path to the SQLite database file.
    """
    if not os.path.exists(path):
        print(f"Error: database file not found: {path}")
        sys.exit(1)
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        tables = [row[0] for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()]
        if REQUIRED_TABLE not in tables:
            print(f"Error: database missing required table '{REQUIRED_TABLE}'. Found tables: {tables}")
            sys.exit(1)
        cursor = conn.execute(f"SELECT * FROM {REQUIRED_TABLE} LIMIT 1")
        db_columns = {desc[0] for desc in cursor.description}
        missing = REQUIRED_COLUMNS - db_columns
        if missing:
            print(f"Error: database table '{REQUIRED_TABLE}' missing required columns: {sorted(missing)}")
            sys.exit(1)
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Data loading functions
# ---------------------------------------------------------------------------


def _normalize_metadata_column(name):
    """Lowercase + drop underscores/spaces, for case/underscore-insensitive matching."""
    return str(name).strip().lower().replace("_", "").replace(" ", "")


def load_sample_table(filepath):
    """Load sample metadata and create lookup dictionaries.

    Permissive on purpose so it accepts exactly what the build accepts (the same
    file is forwarded by trails.py): only ``sample_id`` is required, column names
    are matched case/underscore-insensitively, and ``affected_status`` /
    ``analysis_status`` are optional (treated as unknown when absent). All other
    columns are preserved verbatim on each sample row (e.g. ``ancestry``).

    Returns:
        Tuple (sample_id_to_row, affected_status_lookup, analysis_status_lookup).
    """
    import pandas as pd

    df = pd.read_table(filepath, dtype=str, keep_default_na=False)

    # Canonicalize the columns we need, case/underscore-insensitively.
    column_by_normalized = {_normalize_metadata_column(c): c for c in df.columns}
    id_column = column_by_normalized.get("sampleid")
    if id_column is None:
        raise ValueError(
            f"Sample table {filepath} is missing a 'sample_id' column; "
            f"found columns: {list(df.columns)[:10]}")
    rename = {id_column: "sample_id"}
    for canonical in ("affected_status", "analysis_status", "phenotype_description"):
        actual = column_by_normalized.get(_normalize_metadata_column(canonical))
        if actual and actual != canonical:
            rename[actual] = canonical
    df = df.rename(columns=rename)

    analysis_status_remap = {"rncc": "unsolved", "rcpc": "unsolved", "s_kgfp": "solved"}
    if "analysis_status" in df.columns:
        df["analysis_status"] = (
            df["analysis_status"].astype(str).str.strip().str.lower().replace(analysis_status_remap)
        )
        df["analysis_status"] = df["analysis_status"].replace({"": "unknown", "nan": "unknown"})
    else:
        df["analysis_status"] = "unknown"
    if "affected_status" not in df.columns:
        df["affected_status"] = ""
    if "phenotype_description" in df.columns:
        # Strip only a leading "NA; " prefix (matches input_tables.read_sample_metadata).
        df["phenotype_description"] = df["phenotype_description"].apply(
            lambda p: p[len("NA; "):] if isinstance(p, str) and p.startswith("NA; ") else p
        )
    df = df.drop_duplicates(subset=["sample_id"], keep="first")
    sample_id_to_row = df.set_index("sample_id").to_dict(orient="index")
    affected_status_lookup = dict(
        zip(df.sample_id, df.affected_status.map(normalize_affected_status_for_logic))
    )
    analysis_status_lookup = dict(zip(df.sample_id, df.analysis_status))
    return sample_id_to_row, affected_status_lookup, analysis_status_lookup


def load_mendelian_warnings(db_path, threshold=0.10):
    """Load Mendelian violation data and flag samples with high violation rates.

    The Mendelian-violation tables now live in the SINGLE results database. If the
    table is absent (no trios were available at build time), returns {}.

    Args:
        db_path: Path to the results SQLite database.
        threshold: Violation-rate threshold (default 0.10 = 10%).

    Returns:
        dict mapping sample_id to warning info for samples exceeding the threshold
        on autosome OR chrX.
    """
    if not os.path.exists(db_path):
        return {}

    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        if not conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='mendelian_violations'"
        ).fetchone():
            return {}

        warnings = {}
        for row in conn.execute("SELECT * FROM mendelian_violations"):
            autosome_rate = row["autosome_violations"] / row["autosome_total"] if row["autosome_total"] > 0 else 0
            chrX_rate = row["chrX_violations"] / row["chrX_total"] if row["chrX_total"] > 0 else 0
            chrY_rate = row["chrY_violations"] / row["chrY_total"] if row["chrY_total"] > 0 else 0
            if autosome_rate > threshold or chrX_rate > threshold:
                motif_rates = {}
                for category in ["1bp", "2bp", "3bp", "4bp", "5bp", "6bp", "7_24bp", "25plusbp"]:
                    violations = row[f"motif_{category}_violations"]
                    total = row[f"motif_{category}_total"]
                    if total > 0:
                        motif_rates[category] = {
                            "rate": violations / total,
                            "violations": violations,
                            "total": total,
                        }
                warnings[row["sample_id"]] = {
                    "autosome_rate": autosome_rate,
                    "autosome_violations": row["autosome_violations"],
                    "autosome_total": row["autosome_total"],
                    "chrX_rate": chrX_rate,
                    "chrX_violations": row["chrX_violations"],
                    "chrX_total": row["chrX_total"],
                    "chrY_rate": chrY_rate,
                    "chrY_violations": row["chrY_violations"],
                    "chrY_total": row["chrY_total"],
                    "motif_rates": motif_rates,
                }
        return warnings
    finally:
        conn.close()


def compute_sample_qc_data(db_path):
    """Compute per-sample outlier counts grouped by motif size bin.

    Queries the swim_plot table filtered by outlier_type = 'AllAlleles', groups by
    sample_id and motif size bin, and returns counts for rank1 (largest outlier)
    and top10 (in the top 10 outliers).

    Returns:
        dict with {"rank1": [...], "top10": [...]}, or None if swim_plot / its
        outlier_rank column is missing.
    """
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        if not conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='swim_plot'"
        ).fetchone():
            return None
        columns = {desc[0] for desc in conn.execute("SELECT * FROM swim_plot LIMIT 1").description}
        if "outlier_rank" not in columns:
            return None

        bin_case = """
            CASE
                WHEN MotifSize = 1 THEN '1bp'
                WHEN MotifSize = 2 THEN '2bp'
                WHEN MotifSize >= 3 AND MotifSize <= 6 THEN '3-6bp'
                WHEN MotifSize >= 7 AND MotifSize <= 24 THEN '7-24bp'
                ELSE '25+bp'
            END
        """
        rank1_rows = conn.execute(
            f"SELECT sample_id, {bin_case} AS bin, COUNT(*) AS count FROM swim_plot "
            "WHERE outlier_type = 'AllAlleles' AND outlier_rank = 1 "
            "GROUP BY sample_id, bin ORDER BY sample_id, bin"
        ).fetchall()
        top10_rows = conn.execute(
            f"SELECT sample_id, {bin_case} AS bin, COUNT(*) AS count FROM swim_plot "
            "WHERE outlier_type = 'AllAlleles' AND outlier_rank <= 10 "
            "GROUP BY sample_id, bin ORDER BY sample_id, bin"
        ).fetchall()
        return {
            "rank1": [{"sample_id": r["sample_id"], "bin": r["bin"], "count": r["count"]} for r in rank1_rows],
            "top10": [{"sample_id": r["sample_id"], "bin": r["bin"], "count": r["count"]} for r in top10_rows],
        }
    finally:
        conn.close()


def compute_outlier_warnings(rank1_data):
    """Compute 99th-percentile warnings from rank1 outlier data.

    For each motif size bin, identifies samples at or above the 99th percentile of
    the per-sample largest-outlier count.

    Returns:
        dict mapping sample_id to {"outlier_bins": {bin: {"count", "threshold"}}}.
    """
    bin_counts = collections.defaultdict(list)
    for entry in rank1_data:
        bin_counts[entry["bin"]].append((entry["sample_id"], entry["count"]))

    bin_thresholds = {}
    for bin_name, samples in bin_counts.items():
        counts = [s[1] for s in samples]
        if counts:
            bin_thresholds[bin_name] = float(np.percentile(counts, 99))

    warnings = {}
    for bin_name, samples in bin_counts.items():
        threshold = bin_thresholds.get(bin_name, float("inf"))
        for sample_id, count in samples:
            if count >= threshold:
                warnings.setdefault(sample_id, {"outlier_bins": {}})
                warnings[sample_id]["outlier_bins"][bin_name] = {
                    "count": count,
                    "threshold": int(threshold),
                }
    return warnings


def load_annotations(db_path):
    """Load user annotations (notes/tags) from SQLite, creating the schema if needed.

    Returns:
        dict with keys: notes, tags, all_tags, tag_to_loci.
    """
    conn = sqlite3.connect(db_path)
    conn.execute("""CREATE TABLE IF NOT EXISTS notes (
        locus_id TEXT PRIMARY KEY,
        note_text TEXT NOT NULL,
        created_at TEXT NOT NULL DEFAULT (datetime('now')),
        updated_at TEXT NOT NULL DEFAULT (datetime('now'))
    )""")
    conn.execute("""CREATE TABLE IF NOT EXISTS tags (
        locus_id TEXT NOT NULL,
        tag TEXT NOT NULL,
        created_at TEXT NOT NULL DEFAULT (datetime('now')),
        PRIMARY KEY (locus_id, tag)
    )""")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_tags_tag ON tags(tag)")
    conn.commit()

    notes = {}
    for row in conn.execute("SELECT locus_id, note_text, updated_at FROM notes"):
        notes[row[0]] = {"note_text": row[1], "updated_at": row[2]}

    tags = {}
    tag_to_loci = collections.defaultdict(set)
    for row in conn.execute("SELECT locus_id, tag FROM tags"):
        tags.setdefault(row[0], []).append(row[1])
        tag_to_loci[row[1]].add(row[0])

    all_tags = sorted(tag_to_loci.keys())
    conn.close()
    return {
        "notes": notes,
        "tags": tags,
        "all_tags": all_tags,
        "tag_to_loci": dict(tag_to_loci),
    }


# ---------------------------------------------------------------------------
# Database connection helper
# ---------------------------------------------------------------------------


def get_db():
    """Get a read-only SQLite connection to the single results database."""
    db_path = app.config["DB_PATH"]
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


# ---------------------------------------------------------------------------
# CORS handler
# ---------------------------------------------------------------------------


@app.after_request
def add_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "GET, POST, PUT, DELETE, OPTIONS"
    return response


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def get_motif_variants(motif):
    """Return all rotations and reverse-complement rotations of a motif.

    Used as a fallback when the CanonicalMotif column doesn't exist, so motifs
    can be matched against the Motif column using all equivalent forms.
    """
    variants = set()
    motif = motif.upper()
    for i in range(len(motif)):
        variants.add(motif[i:] + motif[:i])
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rc = "".join(complement.get(b, b) for b in reversed(motif))
    for i in range(len(rc)):
        variants.add(rc[i:] + rc[:i])
    return variants


def validate_params():
    """Validate and parse query parameters from the request.

    The 'source' query parameter is ignored (single-database server): if present it
    is silently dropped rather than rejected, so old multi-source URLs keep working.

    Returns:
        Tuple (params_dict, error_response). If error_response is not None, the
        route handler should return it immediately.
    """
    params = {}
    errors = []

    # outlier_type (required)
    ot_raw = request.args.get("outlier_type")
    if not ot_raw:
        return None, (jsonify({"error": "Missing required parameter", "detail": "outlier_type is required (one of: all, short, hemi)"}), 400)
    if ot_raw not in OUTLIER_TYPE_MAP:
        return None, (jsonify({"error": "Invalid parameter", "detail": f"outlier_type must be one of: all, short, hemi — got '{ot_raw}'"}), 400)
    params["outlier_type"] = ot_raw

    # page (default: 1)
    page_raw = request.args.get("page", "1")
    try:
        params["page"] = int(page_raw)
        if params["page"] < 1:
            raise ValueError
    except ValueError:
        errors.append(f"page must be a positive integer, got '{page_raw}'")

    # page_size (default: 50, max: 500)
    page_size_raw = request.args.get("page_size", "50")
    try:
        params["page_size"] = int(page_size_raw)
        if params["page_size"] < 1 or params["page_size"] > 500:
            raise ValueError
    except ValueError:
        errors.append(f"page_size must be an integer between 1 and 500, got '{page_size_raw}'")

    if request.args.get("min_expansion"):
        try:
            params["min_expansion"] = int(request.args["min_expansion"])
            if params["min_expansion"] < 0:
                raise ValueError
        except ValueError:
            errors.append(f"min_expansion must be an integer >= 0, got '{request.args['min_expansion']}'")

    rau = request.args.get("require_above_unaffected")
    if rau:
        if rau not in ("first", "second", "third"):
            errors.append(f"require_above_unaffected must be one of: first, second, third — got '{rau}'")
        else:
            params["require_above_unaffected"] = rau

    rfau = request.args.get("require_families_above_unaffected")
    if rfau:
        if rfau not in ("first", "second", "third"):
            errors.append(f"require_families_above_unaffected must be one of: first, second, third — got '{rfau}'")
        else:
            params["require_families_above_unaffected"] = rfau

    rap = request.args.get("require_above_population")
    if rap:
        if rap not in ("long-read", "short-read"):
            errors.append(f"require_above_population must be one of: long-read, short-read — got '{rap}'")
        else:
            params["require_above_population"] = rap

    if request.args.get("include_loci_without_population_data") in ("1", "true"):
        params["include_loci_without_population_data"] = True

    pm = request.args.get("population_metric")
    if pm:
        if pm not in ("99th", "max"):
            errors.append(f"population_metric must be one of: 99th, max — got '{pm}'")
        else:
            params["population_metric"] = pm

    if request.args.get("motif"):
        motif_list = [m.strip() for m in request.args["motif"].split(",") if m.strip()]
        if motif_list:
            params["motif"] = motif_list

    if request.args.get("min_repeats_threshold"):
        try:
            params["min_repeats_threshold"] = int(request.args["min_repeats_threshold"])
        except ValueError:
            errors.append(f"min_repeats_threshold must be an integer, got '{request.args['min_repeats_threshold']}'")

    if request.args.get("gene_regions"):
        regions = [r.strip() for r in request.args["gene_regions"].split(",")]
        invalid = [r for r in regions if r not in GENE_REGION_MAP]
        if invalid:
            errors.append(f"gene_regions contains invalid values: {invalid}. Valid: {sorted(GENE_REGION_MAP)}")
        else:
            params["gene_regions"] = regions

    if request.args.get("exclude_gene_regions"):
        regions = [r.strip() for r in request.args["exclude_gene_regions"].split(",")]
        invalid = [r for r in regions if r not in GENE_REGION_MAP]
        if invalid:
            errors.append(f"exclude_gene_regions contains invalid values: {invalid}. Valid: {sorted(GENE_REGION_MAP)}")
        else:
            params["exclude_gene_regions"] = regions

    if request.args.get("min_pli"):
        try:
            params["min_pli"] = float(request.args["min_pli"])
        except ValueError:
            errors.append(f"min_pli must be a number, got '{request.args['min_pli']}'")

    if request.args.get("min_sigma_percentile"):
        try:
            val = float(request.args["min_sigma_percentile"])
            if not 0 <= val <= 1:
                raise ValueError
            params["min_sigma_percentile"] = val
        except ValueError:
            errors.append(f"min_sigma_percentile must be a number between 0 and 1, got '{request.args['min_sigma_percentile']}'")

    for flag_name in ("apply_pathogenic_threshold", "known_loci_only", "known_motifs_only", "mendelian_genes_only"):
        if request.args.get(flag_name, "").lower() in ("true", "1", "yes"):
            params[flag_name] = True

    if request.args.get("phenotype_keyword"):
        params["phenotype_keyword"] = [kw.strip() for kw in request.args["phenotype_keyword"].split(",") if kw.strip()]

    if request.args.get("sample_id_keyword"):
        params["sample_id_keyword"] = [kw.strip() for kw in request.args["sample_id_keyword"].split(",") if kw.strip()]

    if request.args.get("sample_id_like"):
        params["sample_id_like"] = [kw.strip() for kw in request.args["sample_id_like"].split(",") if kw.strip()]

    if request.args.get("locus_id"):
        params["locus_id"] = [lid.strip() for lid in request.args["locus_id"].split(",") if lid.strip()]

    if request.args.get("chrom"):
        params["chrom"] = request.args["chrom"]
    if request.args.get("gene_id"):
        params["gene_id"] = request.args["gene_id"]
    if request.args.get("gene_symbol"):
        params["gene_symbol"] = [gs.strip() for gs in request.args["gene_symbol"].split(",") if gs.strip()]

    if request.args.get("motif_size"):
        clause, params_list, error = parse_motif_size_filter(request.args["motif_size"])
        if error:
            errors.append(error)
        elif clause:
            params["motif_size"] = request.args["motif_size"].strip()
            params["motif_size_clause"] = clause
            params["motif_size_params"] = params_list

    if request.args.get("tag"):
        params["tag"] = request.args["tag"].strip()

    if request.args.get("sort_by"):
        sort_fields = []
        for sort_field in request.args["sort_by"].split(","):
            sort_field = sort_field.strip()
            if not sort_field:
                continue
            if sort_field not in VALID_SORTS:
                errors.append(f"sort_by value '{sort_field}' is invalid. Valid: {sorted(VALID_SORTS)}")
            elif sort_field not in sort_fields:
                sort_fields.append(sort_field)
        if sort_fields:
            params["sort_by"] = sort_fields

    if errors:
        return None, (jsonify({"error": "Invalid parameters", "detail": errors}), 400)

    return params, None


def build_api_order_by(params):
    """Build the ORDER BY clause from API parameters.

    Returns:
        Tuple (order_by_clause, needs_phenotype_join). The second value is always
        False (phenotype scores are denormalized onto loci); it is kept only for
        the caller's signature.
    """
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    source_columns = app.config.get("DB_COLUMNS_SET", set())
    sort_list = list(params.get("sort_by") or ["count"])

    SORT_MAPPING = {
        "region": ("gene_region_rank", "ASC"),
        "size": (f"FirstAffectedAlleleSize_{ot}", "DESC"),
        "count": (f"NumAffectedUnsolvedSamplesAboveUnaffected_{ot}", "DESC"),
        "family_count": (f"NumAffectedUnsolvedFamiliesAboveUnaffected_{ot}", "DESC"),
        "pairwise_similarity": (f"SumPairwiseSim_{ot}", "DESC"),
        "gene_phenotype": (f"MaxGenePhenoSim_{ot}", "DESC"),
        "sigma_hprc_rank": ("HPRC256_StdevPercentile", "ASC NULLS LAST"),
        "sigma_aou_rank": ("AoU1027_StdevPercentile", "ASC NULLS LAST"),
    }

    sequence = list(sort_list)
    for tiebreak in ("count", "region", "size"):
        if tiebreak not in sequence:
            sequence.append(tiebreak)

    order_parts = []
    for field in sequence:
        col, direction = SORT_MAPPING[field]
        if source_columns and col not in source_columns:
            continue
        order_parts.append(f"{col} {direction}")

    if not order_parts:
        return f" ORDER BY FirstAffectedAlleleSize_{ot} DESC", False
    return " ORDER BY " + ", ".join(order_parts), False


def build_api_query(params):
    """Build the SQL queries from validated API query parameters.

    Returns:
        Tuple (select_query, count_query, motif_count_query, gene_region_count_query,
        sql_params, sql_params_with_pagination).
    """
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    source_columns = app.config.get("DB_COLUMNS_SET", set())
    clauses = []
    sql_params = []
    min_expansion = params.get("min_expansion", 0)

    clauses.append(f"(loci.OutlierSampleIds_{ot} IS NOT NULL AND loci.OutlierSampleIds_{ot} != '')")

    if "min_expansion" in params:
        clauses.append(f"loci.FirstAffectedAlleleSize_{ot} >= loci.NumRepeatsInReference + ?")
        sql_params.append(params["min_expansion"])

    if "require_above_unaffected" in params:
        col_prefix = {"first": "First", "second": "Second", "third": "Third"}[params["require_above_unaffected"]]
        clauses.append(f"(loci.{col_prefix}AffectedAlleleSize_{ot} > loci.FirstUnaffectedAlleleSize_{ot} + ? OR (loci.FirstUnaffectedAlleleSize_{ot} IS NULL OR loci.FirstUnaffectedAlleleSize_{ot} = 0) AND loci.{col_prefix}AffectedAlleleSize_{ot} IS NOT NULL)")
        sql_params.append(min_expansion)

    if "require_families_above_unaffected" in params:
        col_prefix = {"first": "First", "second": "Second", "third": "Third"}[params["require_families_above_unaffected"]]
        by_family_col = f"loci.{col_prefix}AffectedAlleleSize_{ot}_ByFamily"
        if not source_columns or f"{col_prefix}AffectedAlleleSize_{ot}_ByFamily" in source_columns:
            clauses.append(f"({by_family_col} > loci.FirstUnaffectedAlleleSize_{ot} + ? OR (loci.FirstUnaffectedAlleleSize_{ot} IS NULL OR loci.FirstUnaffectedAlleleSize_{ot} = 0) AND {by_family_col} IS NOT NULL)")
            sql_params.append(min_expansion)

    if "require_above_population" in params:
        include_without_data = params.get("include_loci_without_population_data", False)
        m = "MaxAllele" if params.get("population_metric") == "max" else "99thPercentile"
        pop_datasets = (["HPRC256", "AoU1027"]
                        if params["require_above_population"] == "long-read" else ["TenK10K"])
        present = [d for d in pop_datasets if f"{d}_{m}" in source_columns]
        above_terms = [f"(loci.{d}_{m} IS NULL OR loci.FirstAffectedAlleleSize_{ot} > loci.{d}_{m})"
                       for d in present]
        if include_without_data:
            if above_terms:
                clauses.append("(" + " AND ".join(above_terms) + ")")
        else:
            if present:
                exists_term = " OR ".join(f"loci.{d}_{m} IS NOT NULL" for d in present)
                clauses.append("((" + exists_term + ") AND " + " AND ".join(above_terms) + ")")
            else:
                clauses.append("1=0")

    if "motif" in params:
        if "CanonicalMotif" in source_columns:
            canonical_motifs = {compute_canonical_motif(m, include_reverse_complement=True) for m in params["motif"]}
            clauses.append(f"CanonicalMotif IN ({','.join('?' * len(canonical_motifs))})")
            sql_params.extend(sorted(canonical_motifs))
        else:
            all_variants = set()
            for m in params["motif"]:
                all_variants.update(get_motif_variants(m))
            clauses.append(f"Motif IN ({','.join('?' * len(all_variants))})")
            sql_params.extend(sorted(all_variants))

    if "min_repeats_threshold" in params:
        clauses.append(f"loci.FirstAffectedAlleleSize_{ot} >= ?")
        sql_params.append(params["min_repeats_threshold"])

    if "gene_regions" in params:
        db_regions = []
        for region in params["gene_regions"]:
            db_regions.extend(GENE_REGION_MAP[region])
        clauses.append(f"loci.gene_region IN ({','.join('?' * len(db_regions))})")
        sql_params.extend(db_regions)

    if "exclude_gene_regions" in params:
        db_regions = []
        for region in params["exclude_gene_regions"]:
            db_regions.extend(GENE_REGION_MAP[region])
        clauses.append(f"loci.gene_region NOT IN ({','.join('?' * len(db_regions))})")
        sql_params.extend(db_regions)

    if "min_pli" in params:
        clauses.append("pLI >= ?")
        sql_params.append(params["min_pli"])

    if params.get("known_loci_only"):
        known_ids = app.config.get("KNOWN_DISEASE_LOCUS_IDS", set())
        if not known_ids:
            clauses.append("1=0")
        else:
            clauses.append(f"LocusId IN ({','.join('?' * len(known_ids))})")
            sql_params.extend(sorted(known_ids))

    if params.get("apply_pathogenic_threshold"):
        thresholds = app.config.get("KNOWN_DISEASE_LOCUS_THRESHOLDS", {})
        if not thresholds:
            clauses.append("1=0")
        else:
            threshold_clauses = []
            for locus_id, threshold in sorted(thresholds.items()):
                threshold_clauses.append(f"(loci.LocusId = ? AND loci.FirstAffectedAlleleSize_{ot} >= ?)")
                sql_params.extend([locus_id, threshold])
            clauses.append(f"({' OR '.join(threshold_clauses)})")

    if params.get("known_motifs_only"):
        clauses.append("loci.IsKnownMotif = 1" if "IsKnownMotif" in source_columns else "1=0")

    if params.get("mendelian_genes_only"):
        clauses.append("IsInMendelianGene = 1" if "IsInMendelianGene" in source_columns else "1=0")

    if "min_sigma_percentile" in params:
        threshold = 1 - params["min_sigma_percentile"]
        hprc_exists = "HPRC256_StdevPercentile" in source_columns
        aou_exists = "AoU1027_StdevPercentile" in source_columns
        if hprc_exists and aou_exists:
            clauses.append("((loci.HPRC256_StdevPercentile IS NOT NULL AND loci.HPRC256_StdevPercentile <= ?) OR (loci.AoU1027_StdevPercentile IS NOT NULL AND loci.AoU1027_StdevPercentile <= ?))")
            sql_params.extend([threshold, threshold])
        elif hprc_exists:
            clauses.append("(loci.HPRC256_StdevPercentile IS NOT NULL AND loci.HPRC256_StdevPercentile <= ?)")
            sql_params.append(threshold)
        elif aou_exists:
            clauses.append("(loci.AoU1027_StdevPercentile IS NOT NULL AND loci.AoU1027_StdevPercentile <= ?)")
            sql_params.append(threshold)

    if "phenotype_keyword" in params:
        kw_clauses = [f"FirstAffectedPhenotype_{ot} LIKE ? COLLATE NOCASE" for _ in params["phenotype_keyword"]]
        clauses.append(f"({' OR '.join(kw_clauses)})")
        sql_params.extend([f"%{kw}%" for kw in params["phenotype_keyword"]])

    if "sample_id_keyword" in params:
        kw_clauses = [
            f"((',' || OutlierSampleIds_{ot} || ',') LIKE ? COLLATE NOCASE"
            f" OR OutlierSampleIds_{ot} LIKE ? COLLATE NOCASE)"
            for _ in params["sample_id_keyword"]
        ]
        clauses.append(f"({' OR '.join(kw_clauses)})")
        for kw in params["sample_id_keyword"]:
            sql_params.extend([f"%:{kw},%", f"%:{kw}:%"])

    if "sample_id_like" in params:
        kw_clauses = [f"OutlierSampleIds_{ot} LIKE ? COLLATE NOCASE" for _ in params["sample_id_like"]]
        clauses.append(f"({' OR '.join(kw_clauses)})")
        sql_params.extend([f"%{kw}%" for kw in params["sample_id_like"]])

    if "locus_id" in params:
        clauses.append(f"LocusId IN ({','.join('?' * len(params['locus_id']))})")
        sql_params.extend(params["locus_id"])

    if "chrom" in params:
        clauses.append("Chrom = ?")
        sql_params.append(params["chrom"])

    if "gene_id" in params:
        clauses.append("gene_id = ?")
        sql_params.append(params["gene_id"])

    if "gene_symbol" in params:
        if "GeneTableGeneSymbol" in app.config["DB_COLUMNS_SET"]:
            gs_clauses = ["GeneTableGeneSymbol LIKE ? COLLATE NOCASE" for _ in params["gene_symbol"]]
            clauses.append(f"({' OR '.join(gs_clauses)})")
            sql_params.extend([f"%{gs}%" for gs in params["gene_symbol"]])
        else:
            # No gene table was built (column absent): a gene-symbol filter can
            # match nothing rather than crash the query.
            clauses.append("0")

    if "motif_size_clause" in params:
        clauses.append(params["motif_size_clause"].replace("MotifSize", "loci.MotifSize"))
        sql_params.extend(params["motif_size_params"])

    if "tag" in params:
        matching = resolve_tag_to_loci(params["tag"])
        if not matching:
            clauses.append("1=0")
        else:
            clauses.append(f"LocusId IN ({','.join('?' * len(matching))})")
            sql_params.extend(sorted(matching))

    if params["outlier_type"] == "hemi":
        clauses.append("Chrom IN ('chrX', 'chrY')")
        clauses.append("HemizygousAlleleHistogram IS NOT NULL")

    where = " AND ".join(clauses) if clauses else "1=1"
    order_by, _ = build_api_order_by(params)
    motif_count_col = "CanonicalMotif" if "CanonicalMotif" in source_columns else "Motif"

    if app.config.get("HAS_SKINNY", False):
        skinny = f"sk_{ot}"
        inner = f"SELECT LocusId FROM {skinny} AS loci WHERE {where}{order_by} LIMIT ? OFFSET ?"
        select_query = (f"SELECT loci.* FROM loci "
                        f"JOIN ({inner}) pick ON loci.LocusId = pick.LocusId{order_by}")
        count_query = f"SELECT COUNT(*) FROM {skinny} AS loci WHERE {where}"
        motif_count_query = f"SELECT {motif_count_col}, COUNT(*) as count FROM {skinny} AS loci WHERE {where} GROUP BY {motif_count_col} ORDER BY count DESC"
        gene_region_count_query = f"SELECT gene_region, COUNT(*) as count FROM {skinny} AS loci WHERE {where} GROUP BY gene_region"
    else:
        select_query = f"SELECT * FROM loci WHERE {where}{order_by} LIMIT ? OFFSET ?"
        count_query = f"SELECT COUNT(*) FROM loci WHERE {where}"
        motif_count_query = f"SELECT {motif_count_col}, COUNT(*) as count FROM loci WHERE {where} GROUP BY {motif_count_col} ORDER BY count DESC"
        gene_region_count_query = f"SELECT gene_region, COUNT(*) as count FROM loci WHERE {where} GROUP BY gene_region"

    sql_params_with_pagination = sql_params + [params["page_size"], (params["page"] - 1) * params["page_size"]]
    return select_query, count_query, motif_count_query, gene_region_count_query, sql_params, sql_params_with_pagination


def row_to_list_dict(row, ot):
    """Convert a sqlite3.Row to a dict for the list endpoint.

    Selects curated columns and renames outlier-type-specific columns by stripping
    the _{ot} suffix. NaN values become None.
    """
    d = {}
    row_dict = dict(row)
    for col in LIST_COLUMNS_STATIC:
        val = row_dict.get(col)
        if isinstance(val, float) and val != val:
            val = None
        d[col] = val
    for col_base in LIST_COLUMNS_OT_SPECIFIC:
        val = row_dict.get(f"{col_base}_{ot}")
        if isinstance(val, float) and val != val:
            val = None
        d[col_base] = val
    return d


def parse_outlier_samples(outlier_string, sample_lookup, affected_lookup, analysis_lookup):
    """Parse an OutlierSampleIds string and enrich with sample metadata.

    The format is comma-separated entries of: allele_size:sample_id[:purity:methylation].

    Returns:
        list of dicts with allele_size, sample_id, family_id, affected_status,
        analysis_status, sex, ancestry, phenotype_description, purity, methylation.
    """
    if not outlier_string:
        return []

    results = []
    for entry in outlier_string.split(","):
        parts = entry.strip().split(":")
        if len(parts) < 2:
            continue
        try:
            allele_size = int(parts[0].replace("x", ""))
        except ValueError:
            continue
        sample_id = parts[1]
        purity = None
        methylation = None
        if len(parts) > 2 and parts[2] != ".":
            try:
                purity = float(parts[2])
            except ValueError:
                pass
        if len(parts) > 3 and parts[3] != ".":
            try:
                methylation = float(parts[3])
            except ValueError:
                pass

        sample_row = sample_lookup.get(sample_id, {})

        # Prefer the raw sample-row affected_status for display/export so
        # "possibly affected" is preserved; affected_lookup collapses it to
        # "affected" for analysis logic only (see swim_plot._normalize_affected_status,
        # which mirrors this). Fall back to affected_lookup when the sample has no
        # usable raw status.
        affected_raw = str(sample_row.get("affected_status") or "").strip()
        if affected_raw.lower() in ("", "nan"):
            affected_raw = affected_lookup.get(sample_id, "unknown") or "unknown"
        if affected_raw.lower() == "affected":
            affected_display = "Affected"
        elif affected_raw.lower() == "unaffected":
            affected_display = "Unaffected"
        else:
            affected_display = affected_raw.title()

        analysis_raw = analysis_lookup.get(sample_id, "unknown")
        if analysis_raw == "solved":
            analysis_display = "Solved"
        elif analysis_raw == "unsolved":
            analysis_display = "Unsolved"
        elif analysis_raw == "unaffected":
            analysis_display = "Unaffected"
        else:
            analysis_display = analysis_raw.title() if analysis_raw else "Unknown"

        results.append({
            "allele_size": allele_size,
            "sample_id": sample_id,
            "family_id": sample_row.get("family_id"),
            "affected_status": affected_display,
            "analysis_status": analysis_display,
            "sex": sample_row.get("sex"),
            "ancestry": sample_row.get("ancestry"),
            "phenotype_description": sample_row.get("phenotype_description"),
            "purity": purity,
            "methylation": methylation,
        })
    return results


def build_sample_details(row_dict, raw_row, conn, lookups, phenotype_scores=None):
    """Build sample details for outliers exceeding the population 99th percentile.

    Args:
        phenotype_scores: Optional pre-fetched {'per_outlier': ..., 'per_locus': ...}
            dict for this locus (see fetch_all_phenotype_scores) — avoids one
            fetch_phenotype_scores() DB round trip per call when the caller already
            bulk-loaded scores for the whole result set. Falls back to fetching
            per-locus from conn when not supplied.

    Returns:
        dict with AllAlleleOutliers / ShortAlleleOutliers / HemizygousOutliers keys
        (only present if qualifying samples exist).
    """
    pop_max = max(
        row_dict.get("HPRC256_99thPercentile") or 0,
        row_dict.get("AoU1027_99thPercentile") or 0,
    )

    locus_id = row_dict.get("LocusId", "")
    if phenotype_scores is None:
        phenotype_scores = fetch_phenotype_scores(conn, locus_id) if conn else {}
    per_outlier_scores = phenotype_scores.get("per_outlier", {})

    outlier_configs = [
        ("OutlierSampleIds_AllAlleles", "AllAlleleOutliers", "AllAlleles"),
        ("OutlierSampleIds_ShortAlleles", "ShortAlleleOutliers", "ShortAlleles"),
        ("OutlierSampleIds_HemizygousAlleles", "HemizygousOutliers", "HemizygousAlleles"),
    ]

    result = {}
    motif = row_dict.get("Motif", "")
    for col_name, output_key, score_key in outlier_configs:
        outlier_str = raw_row[col_name] if col_name in raw_row.keys() else ""
        if not outlier_str:
            continue
        parsed = parse_outlier_samples(outlier_str, lookups["sample"], lookups["affected"], lookups["analysis"])
        score_lookup = {(s["sample_id"], s["allele_size"]): s for s in per_outlier_scores.get(score_key, [])}

        qualifying = []
        for sample in parsed:
            if sample["allele_size"] > pop_max:
                score_data = score_lookup.get((sample["sample_id"], sample["allele_size"]), {})
                qualifying.append({
                    "Allele": sample["allele_size"],
                    "Motif": motif,
                    "Sex": sample.get("sex") or "",
                    "Ancestry": sample.get("ancestry") or "",
                    "Sample ID": sample["sample_id"],
                    "Family": sample.get("family_id") or "",
                    "Status": sample.get("affected_status") or "",
                    "Analysis": sample.get("analysis_status") or "",
                    "Phenotype": sample.get("phenotype_description") or "",
                    "Purity": sample.get("purity"),
                    "Methylation": sample.get("methylation"),
                    "Gene-Pheno Sim": score_data.get("gene_phenotype_similarity"),
                    "Pairwise Sim": score_data.get("pairwise_similarity_to_next"),
                })
        if qualifying:
            result[output_key] = qualifying
    return result


def compute_known_disease_info(row, lookups):
    """Look up known-disease-locus info for a locus.

    Uses Jaccard > 0.66 overlap and length-dependent motif matching, with a STRchive
    fallback. Returns a dict with disease info or None when no match is found.
    """
    locus_id = row.get("LocusId", "")
    chrom = row.get("Chrom", "")
    chrom_key = chrom.replace("chr", "") if chrom else ""
    start_0based = row.get("Start0Based")
    end_1based = row.get("End1Based")
    motif = row.get("Motif", "")

    disease_info = None
    is_strchive = False

    if locus_id in lookups["locus"]:
        disease_info = lookups["locus"][locus_id]
    else:
        coord_key = f"{chrom_key}-{start_0based}-{end_1based}-{motif}"
        if coord_key in lookups["locus"]:
            disease_info = lookups["locus"][coord_key]

    if not disease_info and chrom_key in lookups["disease_trees"] and motif and start_0based is not None and end_1based is not None:
        for interval in lookups["disease_trees"][chrom_key].overlap(start_0based, end_1based):
            if not interval.data:
                continue
            if compute_jaccard(start_0based, end_1based, interval.begin, interval.end) <= 0.66:
                continue
            disease_motifs = [interval.data.get("RepeatUnit", "")] + (interval.data.get("PathogenicMotifs") or [])
            for dm in disease_motifs:
                if motifs_match(motif, dm):
                    disease_info = interval.data
                    break
            if disease_info:
                break

    if not disease_info and lookups.get("strchive_trees") and chrom_key in lookups["strchive_trees"] and motif and start_0based is not None and end_1based is not None:
        for interval in lookups["strchive_trees"][chrom_key].overlap(start_0based, end_1based):
            if not interval.data:
                continue
            if compute_jaccard(start_0based, end_1based, interval.begin, interval.end) <= 0.66:
                continue
            for strchive_motif in interval.data.get("reference_motif_reference_orientation", []):
                if motifs_match(motif, strchive_motif):
                    disease_info = interval.data
                    is_strchive = True
                    break
            if is_strchive:
                break

    if not disease_info:
        return None

    if is_strchive:
        return {
            "locus_id": disease_info.get("locus_id", disease_info.get("id")),
            "pathogenic_min": disease_info.get("pathogenic_min"),
            "inheritance": [disease_info.get("inheritance")] if disease_info.get("inheritance") else None,
            "diseases": [{
                "name": disease_info.get("disease"),
                "symbol": None,
                "inheritance": disease_info.get("inheritance"),
                "pathogenic_min": disease_info.get("pathogenic_min"),
            }],
            "source": "STRchive",
        }

    if not disease_info.get("Diseases"):
        return None

    pathogenic_min = None
    for disease in disease_info.get("Diseases", []):
        if disease.get("PathogenicMin") and (pathogenic_min is None or disease["PathogenicMin"] < pathogenic_min):
            pathogenic_min = disease["PathogenicMin"]

    inheritance_modes = {d["Inheritance"] for d in disease_info.get("Diseases", []) if d.get("Inheritance")}

    return {
        "locus_id": disease_info.get("LocusId"),
        "pathogenic_min": pathogenic_min,
        "inheritance": sorted(inheritance_modes) if inheritance_modes else None,
        "diseases": [{
            "name": d.get("Name"),
            "symbol": d.get("Symbol"),
            "inheritance": d.get("Inheritance"),
            "pathogenic_min": d.get("PathogenicMin"),
        } for d in disease_info.get("Diseases", [])],
        "source": "variant_catalog",
    }


def fetch_phenotype_scores(conn, locus_id):
    """Fetch phenotype scores for a locus, if the score tables exist.

    Returns:
        dict with 'per_outlier' and 'per_locus' keys, or {} when the tables are
        absent.
    """
    tables = {row[0] for row in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall()}
    if "per_outlier_phenotype_scores" not in tables:
        return {}

    result = {"per_outlier": {}, "per_locus": {}}
    rows = conn.execute(
        """SELECT sample_id, outlier_type, allele_size, gene_symbol,
                  gene_phenotype_similarity, gene_phenotype_overlap_count,
                  n_matching_diseases, best_matching_disease, best_disease_inheritance,
                  pairwise_similarity_to_next, pairwise_shared_count_raw,
                  pairwise_shared_count_ic, next_sample_id
           FROM per_outlier_phenotype_scores
           WHERE locus_id = ?
           ORDER BY outlier_type, allele_size DESC""",
        (locus_id,),
    ).fetchall()
    for row in rows:
        outlier_type = row[1]
        score_dict = {
            "sample_id": row[0],
            "allele_size": row[2],
            "gene_symbol": row[3],
            "gene_phenotype_similarity": row[4],
            "gene_phenotype_overlap_count": row[5],
            "n_matching_diseases": row[6],
            "best_matching_disease": row[7],
            "best_disease_inheritance": row[8],
            "pairwise_similarity_to_next": row[9],
            "pairwise_shared_count_raw": row[10],
            "pairwise_shared_count_ic": row[11],
            "next_sample_id": row[12],
        }
        for k, v in score_dict.items():
            if isinstance(v, float) and v != v:
                score_dict[k] = None
        result["per_outlier"].setdefault(outlier_type, []).append(score_dict)

    if "per_locus_phenotype_scores" in tables:
        rows = conn.execute(
            """SELECT outlier_type, num_qualifying_samples, sum_pairwise_similarity,
                      sum_pairwise_shared_raw, sum_pairwise_shared_ic,
                      max_gene_phenotype_similarity, qualifying_sample_ids
               FROM per_locus_phenotype_scores
               WHERE locus_id = ?""",
            (locus_id,),
        ).fetchall()
        for row in rows:
            score_dict = {
                "num_qualifying_samples": row[1],
                "sum_pairwise_similarity": row[2],
                "sum_pairwise_shared_raw": row[3],
                "sum_pairwise_shared_ic": row[4],
                "max_gene_phenotype_similarity": row[5],
                "qualifying_sample_ids": row[6],
            }
            for k, v in score_dict.items():
                if isinstance(v, float) and v != v:
                    score_dict[k] = None
            result["per_locus"][row[0]] = score_dict
    return result


def fetch_all_phenotype_scores(conn):
    """Bulk-load phenotype scores for every locus in one pass.

    Same per-locus output shape as fetch_phenotype_scores(), but issues two queries
    total instead of one query per locus. Used by exports that enrich every row with
    per-sample details, where calling fetch_phenotype_scores() per row would mean N
    SQL round trips for an N-locus export — the per_outlier/per_locus_phenotype_scores
    tables are small enough (low hundreds of thousands of rows) to hold in memory in
    full for the lifetime of one export request.

    Returns:
        dict mapping locus_id to the same {'per_outlier': ..., 'per_locus': ...}
        structure fetch_phenotype_scores() returns for a single locus.
    """
    tables = {row[0] for row in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall()}
    if "per_outlier_phenotype_scores" not in tables:
        return {}

    results = collections.defaultdict(lambda: {"per_outlier": {}, "per_locus": {}})

    rows = conn.execute(
        """SELECT locus_id, sample_id, outlier_type, allele_size, gene_symbol,
                  gene_phenotype_similarity, gene_phenotype_overlap_count,
                  n_matching_diseases, best_matching_disease, best_disease_inheritance,
                  pairwise_similarity_to_next, pairwise_shared_count_raw,
                  pairwise_shared_count_ic, next_sample_id
           FROM per_outlier_phenotype_scores
           ORDER BY locus_id, outlier_type, allele_size DESC"""
    ).fetchall()
    for row in rows:
        locus_id, outlier_type = row[0], row[2]
        score_dict = {
            "sample_id": row[1],
            "allele_size": row[3],
            "gene_symbol": row[4],
            "gene_phenotype_similarity": row[5],
            "gene_phenotype_overlap_count": row[6],
            "n_matching_diseases": row[7],
            "best_matching_disease": row[8],
            "best_disease_inheritance": row[9],
            "pairwise_similarity_to_next": row[10],
            "pairwise_shared_count_raw": row[11],
            "pairwise_shared_count_ic": row[12],
            "next_sample_id": row[13],
        }
        for k, v in score_dict.items():
            if isinstance(v, float) and v != v:
                score_dict[k] = None
        results[locus_id]["per_outlier"].setdefault(outlier_type, []).append(score_dict)

    if "per_locus_phenotype_scores" in tables:
        rows = conn.execute(
            """SELECT locus_id, outlier_type, num_qualifying_samples, sum_pairwise_similarity,
                      sum_pairwise_shared_raw, sum_pairwise_shared_ic,
                      max_gene_phenotype_similarity, qualifying_sample_ids
               FROM per_locus_phenotype_scores"""
        ).fetchall()
        for row in rows:
            locus_id, outlier_type = row[0], row[1]
            score_dict = {
                "num_qualifying_samples": row[2],
                "sum_pairwise_similarity": row[3],
                "sum_pairwise_shared_raw": row[4],
                "sum_pairwise_shared_ic": row[5],
                "max_gene_phenotype_similarity": row[6],
                "qualifying_sample_ids": row[7],
            }
            for k, v in score_dict.items():
                if isinstance(v, float) and v != v:
                    score_dict[k] = None
            results[locus_id]["per_locus"][outlier_type] = score_dict

    return dict(results)


def row_to_full_dict(row):
    """Convert a sqlite3.Row to a dict with all columns (NaN -> None)."""
    d = dict(row)
    for k, v in d.items():
        if isinstance(v, float) and v != v:
            d[k] = None
    return d


# ---------------------------------------------------------------------------
# Page routes
# ---------------------------------------------------------------------------


@app.route("/")
def index():
    """Serve the web UI via on-the-fly Jinja2 template rendering."""
    template = jinja2_env.get_template("index_page_template.html")
    return template.render(
        db_labels=app.config["DB_LABELS"],
        default_source=app.config["DB_LABELS"][0],
        mendelian_warnings=app.config.get("MENDELIAN_WARNINGS", {}),
        outlier_warnings_by_source={app.config["DB_LABELS"][0]: app.config.get("OUTLIER_WARNINGS", {})},
    )


@app.route("/swim")
def swim():
    """Serve the swim plot UI."""
    template = jinja2_env.get_template("swim_plot_template.html")
    return template.render(
        db_labels=app.config["DB_LABELS"],
        default_source=app.config["DB_LABELS"][0],
        mendelian_warnings=app.config.get("MENDELIAN_WARNINGS", {}),
        outlier_warnings_by_source={app.config["DB_LABELS"][0]: app.config.get("OUTLIER_WARNINGS", {})},
    )


@app.route("/sample_qc")
def sample_qc():
    """Serve the sample QC page."""
    template = jinja2_env.get_template("sample_qc_template.html")
    return template.render(
        db_labels=app.config["DB_LABELS"],
        default_source=app.config["DB_LABELS"][0],
        mendelian_warnings=app.config.get("MENDELIAN_WARNINGS", {}),
        outlier_warnings_by_source={app.config["DB_LABELS"][0]: app.config.get("OUTLIER_WARNINGS", {})},
    )


@app.route("/qc2")
def qc2():
    """Serve the QC2 (Mendelian violations) page.

    Hidden gracefully when the Mendelian tables are absent from the database.
    """
    if not app.config.get("HAS_MENDELIAN", False):
        return Response("Mendelian QC page is unavailable: the database has no Mendelian-violation tables.",
                        status=404, mimetype="text/plain")
    template = jinja2_env.get_template("qc2_template.html")
    return template.render(
        mendelian_warnings=app.config.get("MENDELIAN_WARNINGS", {}),
        outlier_warnings_by_source={app.config["DB_LABELS"][0]: app.config.get("OUTLIER_WARNINGS", {})},
    )


# ---------------------------------------------------------------------------
# Annotation CRUD endpoints
# ---------------------------------------------------------------------------


@app.route("/api/v1/annotations/tags", methods=["GET"])
def list_all_tags():
    """List all existing user tags (for autocomplete/dropdown)."""
    user_tags = app.config["ANNOTATIONS"]["all_tags"]
    return jsonify({"tags": user_tags, "system_tags": [], "all_tags": user_tags})


@app.route("/api/v1/annotations/<path:locus_id>", methods=["GET"])
def get_annotations(locus_id):
    """Get note + tags for a locus."""
    annotations = app.config["ANNOTATIONS"]
    result = {}
    if locus_id in annotations["notes"]:
        result["note"] = annotations["notes"][locus_id]
    if locus_id in annotations["tags"]:
        result["tags"] = annotations["tags"][locus_id]
    return jsonify(result)


@app.route("/api/v1/annotations/<path:locus_id>/note", methods=["PUT"])
def upsert_note(locus_id):
    """Upsert a note for a locus."""
    data = request.get_json(force=True)
    note_text = (data.get("note_text") or "").strip()
    if not note_text:
        return jsonify({"error": "note_text is required and cannot be empty"}), 400

    now = datetime.utcnow().isoformat(sep=" ", timespec="seconds")
    conn = sqlite3.connect(app.config["ANNOTATIONS_DB_PATH"])
    conn.execute(
        "INSERT INTO notes (locus_id, note_text, created_at, updated_at) VALUES (?, ?, ?, ?) "
        "ON CONFLICT(locus_id) DO UPDATE SET note_text=excluded.note_text, updated_at=excluded.updated_at",
        (locus_id, note_text, now, now),
    )
    conn.commit()
    conn.close()

    app.config["ANNOTATIONS"]["notes"][locus_id] = {"note_text": note_text, "updated_at": now}
    return jsonify({"locus_id": locus_id, "note": app.config["ANNOTATIONS"]["notes"][locus_id]})


@app.route("/api/v1/annotations/<path:locus_id>/note", methods=["DELETE"])
def delete_note(locus_id):
    """Delete a note for a locus."""
    conn = sqlite3.connect(app.config["ANNOTATIONS_DB_PATH"])
    conn.execute("DELETE FROM notes WHERE locus_id = ?", (locus_id,))
    conn.commit()
    conn.close()
    app.config["ANNOTATIONS"]["notes"].pop(locus_id, None)
    return jsonify({"locus_id": locus_id, "deleted": True})


@app.route("/api/v1/annotations/<path:locus_id>/tags", methods=["POST"])
def add_tag(locus_id):
    """Add a user tag to a locus."""
    data = request.get_json(force=True)
    tag = (data.get("tag") or "").strip()
    if not tag:
        return jsonify({"error": "tag is required and cannot be empty"}), 400

    now = datetime.utcnow().isoformat(sep=" ", timespec="seconds")
    conn = sqlite3.connect(app.config["ANNOTATIONS_DB_PATH"])
    conn.execute("INSERT OR IGNORE INTO tags (locus_id, tag, created_at) VALUES (?, ?, ?)", (locus_id, tag, now))
    conn.commit()
    conn.close()

    annotations = app.config["ANNOTATIONS"]
    annotations["tags"].setdefault(locus_id, [])
    if tag not in annotations["tags"][locus_id]:
        annotations["tags"][locus_id].append(tag)
    annotations["tag_to_loci"].setdefault(tag, set()).add(locus_id)
    annotations["all_tags"] = sorted(annotations["tag_to_loci"].keys())
    return jsonify({"locus_id": locus_id, "tags": annotations["tags"][locus_id]})


@app.route("/api/v1/annotations/<path:locus_id>/tags/<path:tag>", methods=["DELETE"])
def remove_tag(locus_id, tag):
    """Remove a user tag from a locus."""
    conn = sqlite3.connect(app.config["ANNOTATIONS_DB_PATH"])
    conn.execute("DELETE FROM tags WHERE locus_id = ? AND tag = ?", (locus_id, tag))
    conn.commit()
    conn.close()

    annotations = app.config["ANNOTATIONS"]
    if locus_id in annotations["tags"]:
        annotations["tags"][locus_id] = [t for t in annotations["tags"][locus_id] if t != tag]
        if not annotations["tags"][locus_id]:
            del annotations["tags"][locus_id]
    if tag in annotations["tag_to_loci"]:
        annotations["tag_to_loci"][tag].discard(locus_id)
        if not annotations["tag_to_loci"][tag]:
            del annotations["tag_to_loci"][tag]
    annotations["all_tags"] = sorted(annotations["tag_to_loci"].keys())
    return jsonify({"locus_id": locus_id, "tags": annotations["tags"].get(locus_id, [])})


# ---------------------------------------------------------------------------
# Loci endpoints
# ---------------------------------------------------------------------------


@app.route("/api/v1/loci")
def get_loci():
    """Query loci with filters, sorting, and pagination."""
    params, error = validate_params()
    if error:
        return error

    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    select_query, count_query, motif_count_query, gene_region_count_query, sql_params, sql_params_with_pagination = build_api_query(params)

    conn = get_db()
    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
        rows = conn.execute(select_query, sql_params_with_pagination).fetchall()
        motif_rows = conn.execute(motif_count_query, sql_params).fetchall()
        gene_region_rows = conn.execute(gene_region_count_query, sql_params).fetchall()
    finally:
        conn.close()

    motif_counts = {row[0]: row[1] for row in motif_rows if row[0]}
    gene_region_counts = {row[0]: row[1] for row in gene_region_rows if row[0]}

    total_pages = math.ceil(total / params["page_size"]) if total > 0 else 0
    results = [row_to_list_dict(row, ot) for row in rows]

    ann = app.config["ANNOTATIONS"]
    pathogenic_thresholds = app.config.get("KNOWN_DISEASE_LOCUS_THRESHOLDS", {})
    for r in results:
        lid = r.get("LocusId")
        if lid in ann["notes"]:
            r["has_note"] = True
        if lid in ann["tags"]:
            r["tags"] = ann["tags"][lid]
        if lid in pathogenic_thresholds:
            r["PathogenicMin"] = pathogenic_thresholds[lid]

    filters_applied = {k: v for k, v in params.items()
                       if k not in ("page", "page_size", "outlier_type", "sort_by", "motif_size_clause", "motif_size_params")}

    return jsonify({
        "total": total,
        "page": params["page"],
        "page_size": params["page_size"],
        "total_pages": total_pages,
        "outlier_type": ot,
        "filters_applied": filters_applied,
        "results": results,
        "motif_counts": motif_counts,
        "gene_region_counts": gene_region_counts,
    })


def gzip_stream(text_generator):
    """Compress a text-yielding generator into gzip byte chunks, one chunk per yield."""
    compressor = zlib.compressobj(9, zlib.DEFLATED, zlib.MAX_WBITS | 16)
    for chunk in text_generator:
        compressed = compressor.compress(chunk.encode("utf-8"))
        if compressed:
            yield compressed
    yield compressor.flush()


def build_export_filename(outlier_type_key, total, ext):
    """eg. tr_outliers.biallelic.33_loci.20260707_032555.json.gz"""
    search_type = FILENAME_OUTLIER_TYPE_LABELS[outlier_type_key]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"tr_outliers.{search_type}.{total}_loci.{timestamp}.{ext}"


def export_bed(conn, query, count_query, sql_params, params):
    """BED format: chrom, start, end, locus_id, motif, motif_size (no header). Streamed."""
    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
    except Exception:
        conn.close()
        raise
    def generate():
        try:
            for row in conn.execute(query, sql_params):
                r = dict(row)
                yield f"{r['Chrom']}\t{r['Start0Based']}\t{r['End1Based']}\t{r['LocusId']}\t{r['Motif']}\t{r['MotifSize']}\n"
        finally:
            conn.close()
    filename = build_export_filename(params["outlier_type"], total, "bed")
    return Response(generate(), mimetype="text/tab-separated-values",
                    headers={"Content-Disposition": f'attachment; filename="{filename}"'})


def export_tsv(conn, query, count_query, sql_params, ot, params):
    """TSV format: all list columns with header. Streamed, gzip-compressed."""
    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
    except Exception:
        conn.close()
        raise
    def generate():
        try:
            yield "\t".join(LIST_COLUMNS_STATIC + LIST_COLUMNS_OT_SPECIFIC) + "\n"
            for row in conn.execute(query, sql_params):
                row_dict = row_to_list_dict(row, ot)
                yield "\t".join("" if v is None else str(v).replace("\t", " ").replace("\n", " ") for v in row_dict.values()) + "\n"
        finally:
            conn.close()
    filename = build_export_filename(params["outlier_type"], total, "tsv.gz")
    return Response(gzip_stream(generate()), mimetype="application/gzip",
                    headers={"Content-Disposition": f'attachment; filename="{filename}"'})


def export_json(conn, query, count_query, sql_params, ot, params):
    """JSON format: {"metadata": {...}, "loci": [...]}.

    Every locus is always enriched with a "Samples" key listing every qualifying
    outlier sample (allele size, sample/family id, sex, ancestry, affected/analysis
    status, phenotype, purity, methylation, and phenotype-similarity scores) —
    regardless of result-set size. Phenotype scores are bulk-loaded once up front
    (fetch_all_phenotype_scores) rather than per row, so this stays cheap even for
    million-row exports. Rows still stream from the cursor so the full result set
    is never materialized in memory.
    """
    tags_dict = app.config["ANNOTATIONS"]["tags"]
    lookups = app.config["LOOKUPS"]

    non_filter_keys = {"source", "outlier_type", "sort_by", "page", "page_size",
                       "motif_size_clause", "motif_size_params"}
    filters_applied = {k: v for k, v in params.items() if k not in non_filter_keys}

    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
    except Exception:
        conn.close()
        raise

    metadata = {
        "source": app.config["DB_LABELS"][0],
        "outlier_type": params["outlier_type"],
        "sort_by": params.get("sort_by") or ["count"],
        "filters_applied": filters_applied,
        "total_loci": total,
        "generated_at": datetime.now().isoformat(),
    }

    def to_dict(row):
        row_dict = row_to_list_dict(row, ot)
        row_dict = {k: (None if isinstance(v, float) and (v != v) else v) for k, v in row_dict.items()}
        row_dict["Tags"] = ",".join(tags_dict.get(row_dict.get("LocusId", ""), []))
        return row_dict

    def generate():
        try:
            all_phenotype_scores = fetch_all_phenotype_scores(conn)
            yield '{\n  "metadata": ' + json.dumps(metadata, indent=2).replace("\n", "\n  ") + ',\n  "loci": [\n'
            first = True
            for raw_row in conn.execute(query, sql_params):
                row_dict = to_dict(raw_row)
                phenotype_scores = all_phenotype_scores.get(row_dict.get("LocusId", ""), {})
                sample_details = build_sample_details(row_dict, raw_row, conn, lookups, phenotype_scores=phenotype_scores)
                if sample_details:
                    row_dict["Samples"] = sample_details
                yield ("" if first else ",\n") + "    " + json.dumps(row_dict, indent=2).replace("\n", "\n    ")
                first = False
            yield "\n  ]\n}"
        finally:
            conn.close()
    filename = build_export_filename(params["outlier_type"], total, "json.gz")
    return Response(gzip_stream(generate()), mimetype="application/gzip",
                    headers={"Content-Disposition": f'attachment; filename="{filename}"'})


def compute_plot_thresholds(row_dict, all_allele_outliers_str, short_allele_outliers_str, affected_lookup):
    """Compute PlotReadVisualization thresholds for the ExpansionHunter export.

    Returns a dict with LongAllele/ShortAllele thresholds, or None if no criteria met.
    """
    def parse_outliers_simple(outlier_str):
        if not outlier_str:
            return []
        result = []
        for entry in outlier_str.split(","):
            parts = entry.strip().split(":")
            if len(parts) >= 2:
                try:
                    result.append((int(parts[0].replace("x", "")), parts[1]))
                except (ValueError, IndexError):
                    continue
        return result

    def find_threshold(outliers, lookup, comparison_value):
        if not outliers:
            return None
        allele_size, sample_id = outliers[0]
        if lookup.get(sample_id, "") == "affected" and allele_size > comparison_value:
            return allele_size
        return None

    long_outliers = parse_outliers_simple(all_allele_outliers_str)
    short_outliers = parse_outliers_simple(short_allele_outliers_str)
    if not long_outliers and not short_outliers:
        return None

    long_unaffected = [a for a, s in long_outliers if affected_lookup.get(s, "") != "affected"]
    short_unaffected = [a for a, s in short_outliers if affected_lookup.get(s, "") != "affected"]
    max_long_unaffected = max(long_unaffected) if long_unaffected else 0
    max_short_unaffected = max(short_unaffected) if short_unaffected else 0

    long_comparison = max(
        max_long_unaffected,
        row_dict.get("HPRC256_99thPercentile") or 0,
        row_dict.get("AoU1027_99thPercentile") or 0,
    )
    long_threshold = find_threshold(long_outliers, affected_lookup, long_comparison)
    short_threshold = find_threshold(short_outliers, affected_lookup, max_short_unaffected)

    if long_threshold is None and short_threshold is None:
        return None

    result = {}
    if long_threshold is not None:
        result["LongAllele"] = long_threshold
    if short_threshold is not None:
        result["ShortAllele"] = short_threshold
    return result


def export_expansion_hunter(conn, query, count_query, sql_params, ot, params):
    """ExpansionHunter variant catalog format. Enriches when <=200 loci, else streams."""
    tags_dict = app.config["ANNOTATIONS"]["tags"]
    affected_lookup = app.config["LOOKUPS"]["affected"]

    def to_dict(row):
        row_dict = row_to_list_dict(row, ot)
        row_dict = {k: (None if isinstance(v, float) and (v != v) else v) for k, v in row_dict.items()}
        row_dict["Tags"] = ",".join(tags_dict.get(row_dict.get("LocusId", ""), []))
        row_dict["ReferenceRegion"] = f"{row_dict['Chrom']}:{row_dict['Start0Based']}-{row_dict['End1Based']}"
        row_dict["LocusStructure"] = f"({row_dict['Motif']})*"
        row_dict["VariantType"] = "Repeat"
        return row_dict

    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
    except Exception:
        conn.close()
        raise

    if total <= 200:
        try:
            rows = conn.execute(query, sql_params).fetchall()
            results = [to_dict(row) for row in rows]
            lookups = app.config["LOOKUPS"]
            for row_dict, raw_row in zip(results, rows):
                thresholds = compute_plot_thresholds(
                    row_dict,
                    raw_row["OutlierSampleIds_AllAlleles"] or "",
                    raw_row["OutlierSampleIds_ShortAlleles"] or "",
                    affected_lookup,
                )
                if thresholds:
                    row_dict["PlotReadVisualization"] = [
                        {"If": allele, "Is": ">=", "Threshold": thresholds[allele]}
                        for allele in ("LongAllele", "ShortAllele") if allele in thresholds
                    ]
                sample_details = build_sample_details(row_dict, raw_row, conn, lookups)
                if sample_details:
                    row_dict["Samples"] = sample_details
        finally:
            conn.close()
        filename = build_export_filename(params["outlier_type"], total, "expansion_hunter.json")
        return Response(json.dumps(results, indent=2), mimetype="application/json",
                        headers={"Content-Disposition": f'attachment; filename="{filename}"'})

    def generate():
        try:
            yield "[\n"
            first = True
            for row in conn.execute(query, sql_params):
                yield ("" if first else ",\n") + json.dumps(to_dict(row), indent=2)
                first = False
            yield "\n]"
        finally:
            conn.close()
    filename = build_export_filename(params["outlier_type"], total, "expansion_hunter.json")
    return Response(generate(), mimetype="application/json",
                    headers={"Content-Disposition": f'attachment; filename="{filename}"'})


@app.route("/api/v1/export")
def export_loci():
    """Export loci in various formats: bed, tsv, json, expansion_hunter."""
    format_param = request.args.get("format", "").lower()
    valid_formats = {"bed", "tsv", "json", "expansion_hunter"}
    if format_param not in valid_formats:
        return jsonify({"error": f"format must be one of: {', '.join(sorted(valid_formats))}"}), 400

    params, error_response = validate_params()
    if error_response:
        return error_response

    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    select_query, count_query, _, _, sql_params, _ = build_api_query(params)
    base_query = select_query.replace(" LIMIT ? OFFSET ?", "")

    conn = get_db()
    if format_param == "bed":
        return export_bed(conn, base_query, count_query, sql_params, params)
    if format_param == "tsv":
        return export_tsv(conn, base_query, count_query, sql_params, ot, params)
    if format_param == "json":
        return export_json(conn, base_query, count_query, sql_params, ot, params)
    return export_expansion_hunter(conn, base_query, count_query, sql_params, ot, params)


@app.route("/api/v1/loci/<locus_id>")
def get_locus_detail(locus_id):
    """Get full detail for a single locus."""
    conn = get_db()
    try:
        row = conn.execute("SELECT * FROM loci WHERE LocusId = ?", (locus_id,)).fetchone()
        if row is None:
            return jsonify({"error": "Not found", "detail": f"No locus found with LocusId: {locus_id}"}), 404
        row_dict = row_to_full_dict(row)
        phenotype_scores = fetch_phenotype_scores(conn, locus_id)
    finally:
        conn.close()

    lookups = app.config["LOOKUPS"]
    outlier_samples = {}
    for ot_suffix in OUTLIER_TYPE_MAP.values():
        outlier_samples[ot_suffix] = parse_outlier_samples(
            row_dict.get(f"OutlierSampleIds_{ot_suffix}", ""),
            lookups["sample"], lookups["affected"], lookups["analysis"],
        )

    known_disease = compute_known_disease_info(row_dict, lookups)

    population_thresholds = {}
    for cohort in ("HPRC256", "AoU1027", "TenK10K"):
        for stat in ("MaxAllele", "99thPercentile", "90thPercentile", "Median", "Mode", "Stdev"):
            population_thresholds[f"{cohort}_{stat}"] = row_dict.get(f"{cohort}_{stat}")

    ann = app.config["ANNOTATIONS"]
    user_annotations = {}
    if locus_id in ann["notes"]:
        user_annotations["note"] = ann["notes"][locus_id]
    if locus_id in ann["tags"]:
        user_annotations["tags"] = ann["tags"][locus_id]

    return jsonify({
        "locus": row_dict,
        "system_tags": [],
        "outlier_samples": outlier_samples,
        "annotations": {
            "known_disease_locus": known_disease,
            "population_thresholds": population_thresholds,
        },
        "user_annotations": user_annotations,
        "phenotype_scores": phenotype_scores,
    })


@app.route("/api/v1/swim_plot_data")
def get_swim_plot_data():
    """Return swim-plot data: one entry per outlier sample per motif category."""
    outlier_type = request.args.get("outlier_type", "all")
    if outlier_type not in OUTLIER_TYPE_MAP:
        return jsonify({"error": "Invalid parameter", "detail": f"outlier_type must be one of: all, short, hemi — got '{outlier_type}'"}), 400
    ot = OUTLIER_TYPE_MAP[outlier_type]

    require_affected_above_unaffected = request.args.get("require_affected_above_unaffected", "")
    apply_pathogenic_threshold = request.args.get("apply_pathogenic_threshold", "")
    known_loci_only = request.args.get("known_loci_only", "")
    mendelian_only = request.args.get("mendelian_only", "")
    known_motifs_only = request.args.get("known_motifs_only", "")
    require_above_population = request.args.get("require_above_population", "")
    include_loci_without_population_data = request.args.get("include_loci_without_population_data", "") in ("1", "true")
    population_metric = "MaxAllele" if request.args.get("population_metric", "") == "max" else "99thPercentile"
    gene_regions_raw = request.args.get("gene_regions", "")
    exclude_gene_regions_raw = request.args.get("exclude_gene_regions", "")
    motif_size_raw = request.args.get("motif_size", "")
    motif_raw = request.args.get("motif", "")
    locus_id_raw = request.args.get("locus_id", "")
    gene_symbol_raw = request.args.get("gene_symbol", "")
    gene_id = request.args.get("gene_id", "")
    min_pli = request.args.get("min_pli", "")
    phenotype_keyword_raw = request.args.get("phenotype_keyword", "")
    sample_id_keyword_raw = request.args.get("sample_id_keyword", "")
    sample_id_like_raw = request.args.get("sample_id_like", "")
    min_expansion = request.args.get("min_expansion", "")
    min_sigma_percentile_raw = request.args.get("min_sigma_percentile", "")

    source_columns = app.config.get("DB_COLUMNS_SET", set())
    extra_clauses = []
    extra_params = []

    if require_affected_above_unaffected == "1":
        extra_clauses.append("""(
            affected_status IN ('Affected', 'Possibly Affected')
            AND allele_size > COALESCE(
                (SELECT MAX(inner_sp.allele_size)
                 FROM swim_plot inner_sp
                 WHERE inner_sp.LocusId = swim_plot.LocusId
                   AND inner_sp.outlier_type = swim_plot.outlier_type
                   AND inner_sp.affected_status = 'Unaffected'),
                0
            )
        )""")

    if known_loci_only == "1":
        known_ids = app.config.get("KNOWN_DISEASE_LOCUS_IDS", set())
        if not known_ids:
            extra_clauses.append("1=0")
        else:
            extra_clauses.append(f"LocusId IN ({','.join('?' * len(known_ids))})")
            extra_params.extend(sorted(known_ids))

    if apply_pathogenic_threshold == "1":
        thresholds = app.config.get("KNOWN_DISEASE_LOCUS_THRESHOLDS", {})
        if not thresholds:
            extra_clauses.append("1=0")
        else:
            threshold_clauses = []
            for locus_id, threshold in sorted(thresholds.items()):
                threshold_clauses.append("(LocusId = ? AND allele_size >= ?)")
                extra_params.extend([locus_id, threshold])
            extra_clauses.append(f"({' OR '.join(threshold_clauses)})")

    if mendelian_only == "1":
        extra_clauses.append("IsInMendelianGene = 1")

    if known_motifs_only == "1":
        extra_clauses.append("IsKnownMotif = 1")

    if require_above_population == "long-read":
        if include_loci_without_population_data:
            extra_clauses.append(f"""(
                (HPRC256_{population_metric} IS NULL OR allele_size > HPRC256_{population_metric}) AND
                (AoU1027_{population_metric} IS NULL OR allele_size > AoU1027_{population_metric})
            )""")
        else:
            extra_clauses.append(f"""(
                (HPRC256_{population_metric} IS NOT NULL OR AoU1027_{population_metric} IS NOT NULL) AND
                (HPRC256_{population_metric} IS NULL OR allele_size > HPRC256_{population_metric}) AND
                (AoU1027_{population_metric} IS NULL OR allele_size > AoU1027_{population_metric})
            )""")
    elif require_above_population == "short-read":
        if include_loci_without_population_data:
            extra_clauses.append(f"(TenK10K_{population_metric} IS NULL OR allele_size > TenK10K_{population_metric})")
        else:
            extra_clauses.append(f"(TenK10K_{population_metric} IS NOT NULL AND allele_size > TenK10K_{population_metric})")

    if gene_regions_raw:
        db_regions = []
        for r in (x.strip() for x in gene_regions_raw.split(",") if x.strip()):
            db_regions.extend(GENE_REGION_MAP.get(r, []))
        if db_regions:
            extra_clauses.append(f"gene_region IN ({','.join('?' * len(db_regions))})")
            extra_params.extend(db_regions)

    if exclude_gene_regions_raw:
        db_regions = []
        for r in (x.strip() for x in exclude_gene_regions_raw.split(",") if x.strip()):
            db_regions.extend(GENE_REGION_MAP.get(r, []))
        if db_regions:
            extra_clauses.append(f"gene_region NOT IN ({','.join('?' * len(db_regions))})")
            extra_params.extend(db_regions)

    if motif_size_raw:
        clause, params_list, error = parse_motif_size_filter(motif_size_raw)
        if error:
            return jsonify({"error": "Invalid parameter", "detail": error}), 400
        if clause:
            extra_clauses.append(clause)
            extra_params.extend(params_list)

    if motif_raw:
        motifs = [m.strip() for m in motif_raw.split(",") if m.strip()]
        if motifs:
            canonical_motifs = {compute_canonical_motif(m, include_reverse_complement=True) for m in motifs}
            extra_clauses.append(f"CanonicalMotif IN ({','.join('?' * len(canonical_motifs))})")
            extra_params.extend(sorted(canonical_motifs))

    if locus_id_raw:
        locus_ids = [lid.strip() for lid in locus_id_raw.split(",") if lid.strip()]
        if locus_ids:
            extra_clauses.append(f"LocusId IN ({','.join('?' * len(locus_ids))})")
            extra_params.extend(locus_ids)

    if gene_symbol_raw:
        symbols = [gs.strip() for gs in gene_symbol_raw.split(",") if gs.strip()]
        if symbols:
            gs_clauses = ["GeneTableGeneSymbol LIKE ? COLLATE NOCASE" for _ in symbols]
            extra_clauses.append(f"({' OR '.join(gs_clauses)})")
            extra_params.extend([f"%{gs}%" for gs in symbols])

    if gene_id:
        extra_clauses.append("gene_id = ?")
        extra_params.append(gene_id)

    if min_pli:
        try:
            extra_params.append(float(min_pli))
            extra_clauses.append("pLI >= ?")
        except ValueError:
            return jsonify({"error": "Invalid parameter", "detail": f"min_pli must be a number, got '{min_pli}'"}), 400

    if phenotype_keyword_raw:
        keywords = [kw.strip() for kw in phenotype_keyword_raw.split(",") if kw.strip()]
        if keywords:
            kw_clauses = ["phenotype_description LIKE ? COLLATE NOCASE" for _ in keywords]
            extra_clauses.append(f"({' OR '.join(kw_clauses)})")
            extra_params.extend([f"%{kw}%" for kw in keywords])

    if sample_id_keyword_raw:
        keywords = [kw.strip() for kw in sample_id_keyword_raw.split(",") if kw.strip()]
        if keywords:
            kw_clauses = ["sample_id = ? COLLATE NOCASE" for _ in keywords]
            extra_clauses.append(f"({' OR '.join(kw_clauses)})")
            extra_params.extend(keywords)

    if sample_id_like_raw:
        keywords = [kw.strip() for kw in sample_id_like_raw.split(",") if kw.strip()]
        if keywords:
            kw_clauses = ["sample_id LIKE ? COLLATE NOCASE" for _ in keywords]
            extra_clauses.append(f"({' OR '.join(kw_clauses)})")
            extra_params.extend([f"%{kw}%" for kw in keywords])

    if min_expansion:
        try:
            min_expansion_value = int(min_expansion)
            if min_expansion_value < 0:
                raise ValueError
        except ValueError:
            return jsonify({"error": "Invalid parameter", "detail": f"min_expansion must be an integer >= 0, got '{min_expansion}'"}), 400
        extra_clauses.append("allele_size >= NumRepeatsInReference + ?")
        extra_params.append(min_expansion_value)

    if min_sigma_percentile_raw:
        try:
            val = float(min_sigma_percentile_raw)
            if not 0 <= val <= 1:
                raise ValueError
            threshold = 1 - val
            hprc_exists = "HPRC256_StdevPercentile" in source_columns
            aou_exists = "AoU1027_StdevPercentile" in source_columns
            if hprc_exists and aou_exists:
                extra_clauses.append("((HPRC256_StdevPercentile IS NOT NULL AND HPRC256_StdevPercentile <= ?) OR (AoU1027_StdevPercentile IS NOT NULL AND AoU1027_StdevPercentile <= ?))")
                extra_params.extend([threshold, threshold])
            elif hprc_exists:
                extra_clauses.append("(HPRC256_StdevPercentile IS NOT NULL AND HPRC256_StdevPercentile <= ?)")
                extra_params.append(threshold)
            elif aou_exists:
                extra_clauses.append("(AoU1027_StdevPercentile IS NOT NULL AND AoU1027_StdevPercentile <= ?)")
                extra_params.append(threshold)
        except ValueError:
            return jsonify({"error": "Invalid parameter", "detail": f"min_sigma_percentile must be a number between 0 and 1, got '{min_sigma_percentile_raw}'"}), 400

    tag_raw = request.args.get("tag", "").strip()
    if tag_raw:
        matching = resolve_tag_to_loci(tag_raw)
        if not matching:
            extra_clauses.append("1=0")
        else:
            extra_clauses.append(f"LocusId IN ({','.join('?' * len(matching))})")
            extra_params.extend(sorted(matching))

    extra_where = (" AND " + " AND ".join(extra_clauses)) if extra_clauses else ""
    # "Unknown" is the motif_category swim_plot assigns rows whose MotifSize is
    # missing (swim_plot._compute_motif_category); include it so those rows are
    # not silently dropped from the endpoint.
    categories = [f"{i}bp" for i in range(1, 25)] + ["25+bp", "Unknown"]
    all_data = []

    conn = get_db()
    try:
        if not conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='swim_plot'"
        ).fetchone():
            return jsonify({"error": "swim_plot table not found", "detail": "The database has no swim_plot table."}), 500

        for motif_category in categories:
            query = ("""SELECT rowid, allele_size, motif_category, affected_status,
                       LocusId, sample_id, Motif, gene_region, GeneTableGeneSymbol,
                       purity, methylation, family_id, sex, analysis_status, phenotype_description
                FROM swim_plot
                WHERE outlier_type = ? AND motif_category = ?""" + extra_where + """
                ORDER BY allele_size DESC LIMIT 500""")
            for row in conn.execute(query, [ot, motif_category] + extra_params).fetchall():
                entry = {
                    "rowid": row["rowid"],
                    "allele_size": row["allele_size"],
                    "motif_category": row["motif_category"],
                    "affected_status": row["affected_status"],
                    "LocusId": row["LocusId"],
                    "sample_id": row["sample_id"],
                    "Motif": row["Motif"],
                    "gene_region": row["gene_region"],
                    "GeneTableGeneSymbol": row["GeneTableGeneSymbol"],
                    "purity": row["purity"],
                    "methylation": row["methylation"],
                    "family_id": row["family_id"],
                    "sex": row["sex"],
                    "analysis_status": row["analysis_status"],
                    "phenotype_description": row["phenotype_description"],
                }
                for key, val in entry.items():
                    if isinstance(val, float) and val != val:
                        entry[key] = None
                all_data.append(entry)
    finally:
        conn.close()

    ann = app.config["ANNOTATIONS"]
    for entry in all_data:
        lid = entry.get("LocusId")
        if lid in ann["notes"]:
            entry["has_note"] = True
        if lid in ann["tags"]:
            entry["tags"] = ann["tags"][lid]

    return jsonify({"categories": categories, "data": all_data})


@app.route("/api/v1/sample_qc_data")
def get_sample_qc_data():
    """Return per-sample outlier counts grouped by motif size bin (pre-computed at startup)."""
    # configure_app always sets SAMPLE_QC_CACHE to a dict (empty rank1/top10 when
    # there is no swim_plot table), so an empty payload is returned in that case.
    return jsonify(app.config.get("SAMPLE_QC_CACHE") or {"rank1": [], "top10": []})


@app.route("/api/v1/qc2_data")
def get_qc2_data():
    """Return Mendelian-violation fractions grouped by motif size, chromosome, and center.

    Reads the mendelian_violations table from the single results database. Returns a
    500 with a controlled message when that table is absent.
    """
    cached = app.config.get("QC2_CACHE")
    if cached is not None:
        return jsonify(cached)
    if not app.config.get("HAS_MENDELIAN", False):
        return jsonify({"error": "mendelian_violations table not found",
                        "detail": "The database does not contain the Mendelian-violation tables."}), 500

    conn = get_db()
    try:
        rows = conn.execute("SELECT * FROM mendelian_violations").fetchall()
    finally:
        conn.close()

    motif_columns = [
        ("motif_1bp_violations", "motif_1bp_total", "1bp"),
        ("motif_2bp_violations", "motif_2bp_total", "2bp"),
        ("motif_3bp_violations", "motif_3bp_total", "3bp"),
        ("motif_4bp_violations", "motif_4bp_total", "4bp"),
        ("motif_5bp_violations", "motif_5bp_total", "5bp"),
        ("motif_6bp_violations", "motif_6bp_total", "6bp"),
        ("motif_7_24bp_violations", "motif_7_24bp_total", "7-24bp"),
        ("motif_25plusbp_violations", "motif_25plusbp_total", "25+bp"),
    ]
    chrom_columns = [
        ("autosome_violations", "autosome_total", "autosome"),
        ("chrX_violations", "chrX_total", "chrX"),
        ("chrY_violations", "chrY_total", "chrY"),
    ]

    by_motif = []
    for row in rows:
        for viol_col, total_col, bin_name in motif_columns:
            violations = row[viol_col] or 0
            total = row[total_col] or 0
            by_motif.append({
                "sample_id": row["sample_id"], "bin": bin_name,
                "fraction": violations / total if total > 0 else 0,
                "violations": violations, "total": total,
            })

    by_chrom = []
    for row in rows:
        for viol_col, total_col, bin_name in chrom_columns:
            violations = row[viol_col] or 0
            total = row[total_col] or 0
            by_chrom.append({
                "sample_id": row["sample_id"], "bin": bin_name,
                "fraction": violations / total if total > 0 else 0,
                "violations": violations, "total": total,
            })

    def get_center(sample_id):
        return sample_id[:3].replace("_", "")

    center_trio_counts = {}
    for row in rows:
        center = get_center(row["sample_id"])
        center_trio_counts[center] = center_trio_counts.get(center, 0) + 1

    center_motif_agg = {}
    for row in rows:
        center = get_center(row["sample_id"])
        for viol_col, total_col, bin_name in motif_columns:
            agg = center_motif_agg.setdefault((center, bin_name), {"violations": 0, "total": 0})
            agg["violations"] += row[viol_col] or 0
            agg["total"] += row[total_col] or 0

    by_center_motif = []
    for (center, bin_name), agg in center_motif_agg.items():
        trios = center_trio_counts.get(center, 0)
        if trios < 5:
            continue
        by_center_motif.append({
            "center": center, "bin": bin_name,
            "fraction": agg["violations"] / agg["total"] if agg["total"] > 0 else 0,
            "violations": agg["violations"], "total": agg["total"], "trios": trios,
        })

    center_chrom_agg = {}
    for row in rows:
        center = get_center(row["sample_id"])
        for viol_col, total_col, bin_name in chrom_columns:
            agg = center_chrom_agg.setdefault((center, bin_name), {"violations": 0, "total": 0})
            agg["violations"] += row[viol_col] or 0
            agg["total"] += row[total_col] or 0

    by_center_chrom = []
    for (center, bin_name), agg in center_chrom_agg.items():
        trios = center_trio_counts.get(center, 0)
        if trios < 5:
            continue
        by_center_chrom.append({
            "center": center, "bin": bin_name,
            "fraction": agg["violations"] / agg["total"] if agg["total"] > 0 else 0,
            "violations": agg["violations"], "total": agg["total"], "trios": trios,
        })

    result = {
        "by_motif": by_motif,
        "by_chrom": by_chrom,
        "by_center_motif": by_center_motif,
        "by_center_chrom": by_center_chrom,
    }
    app.config["QC2_CACHE"] = result
    return jsonify(result)


@app.route("/api/v1/sample_outlier_stats")
def get_sample_outlier_stats():
    """Return counts of loci where a sample is an outlier above the first unaffected.

    Query params: sample_id, locus_id, outlier_type.
    """
    sample_id = request.args.get("sample_id")
    locus_id = request.args.get("locus_id")
    outlier_type_raw = request.args.get("outlier_type", "all")
    if outlier_type_raw not in OUTLIER_TYPE_MAP:
        return jsonify({"error": "Invalid parameter", "detail": "outlier_type must be one of: all, short, hemi"}), 400
    ot = OUTLIER_TYPE_MAP[outlier_type_raw]

    if not sample_id:
        return jsonify({"error": "Missing required parameter", "detail": "sample_id is required"}), 400

    conn = get_db()
    try:
        if not conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='swim_plot'"
        ).fetchone():
            return jsonify({"error": "swim_plot table not found", "detail": "The database has no swim_plot table."}), 500

        locus_row = None
        if locus_id:
            locus_row = conn.execute(
                "SELECT MotifSize, CanonicalMotif FROM swim_plot WHERE LocusId = ? AND outlier_type = ? LIMIT 1",
                (locus_id, ot),
            ).fetchone()

        motif_size = dict(locus_row)["MotifSize"] if locus_row else None
        canonical_motif = dict(locus_row)["CanonicalMotif"] if locus_row else None

        row = conn.execute(
            """SELECT
                COUNT(DISTINCT LocusId) as total,
                COUNT(DISTINCT CASE WHEN MotifSize = ? THEN LocusId END) as same_motif_size,
                COUNT(DISTINCT CASE WHEN CanonicalMotif = ? THEN LocusId END) as same_canonical_motif
            FROM swim_plot
            WHERE sample_id = ? AND outlier_type = ? AND is_above_first_unaffected = 1""",
            (motif_size, canonical_motif, sample_id, ot),
        ).fetchone()
        result = dict(row)
    finally:
        conn.close()

    return jsonify({
        "total_loci_above_unaffected": result["total"],
        "same_motif_size_loci": result["same_motif_size"],
        "same_canonical_motif_loci": result["same_canonical_motif"],
        "motif_size": motif_size,
        "canonical_motif": canonical_motif,
    })


@app.route("/api/v1/schema")
def get_schema():
    """Return available filter fields, sort options, and enum values."""
    label = app.config["DB_LABELS"][0]
    return jsonify({
        "source": label,
        "available_sources": app.config["DB_LABELS"],
        "filters": {
            "outlier_type": {
                "type": "enum",
                "required": True,
                "values": ["all", "short", "hemi"],
                "labels": {"all": "Long Allele", "short": "Biallelic", "hemi": "Hemizygous"},
            },
            "gene_regions": {
                "type": "multi_enum",
                "values": ["cds", "promoter", "utr", "intron", "exon", "intergenic"],
                "labels": {"cds": "CDS", "promoter": "Promoter", "utr": "UTR (5' & 3')", "intron": "Intron", "exon": "Exon", "intergenic": "Intergenic"},
            },
            "exclude_gene_regions": {
                "type": "multi_enum",
                "values": ["cds", "promoter", "utr", "intron", "exon", "intergenic"],
            },
            "require_above_unaffected": {"type": "enum", "values": ["first", "second", "third"]},
            "require_families_above_unaffected": {"type": "enum", "values": ["first", "second", "third"]},
            "require_above_population": {
                "type": "enum",
                "values": ["long-read", "short-read"],
                "labels": {"long-read": "Long-read (HPRC256, AoU1027)", "short-read": "Short-read (TenK10K)"},
            },
            "include_loci_without_population_data": {
                "type": "boolean",
                "description": "When require_above_population is set, also include loci without population data",
            },
            "min_expansion": {"type": "int", "min": 0},
            "min_repeats_threshold": {"type": "int", "min": 0},
            "min_pli": {"type": "float", "min": 0.0, "max": 1.0},
            "motif": {"type": "string"},
            "phenotype_keyword": {"type": "string"},
            "sample_id_keyword": {"type": "string", "description": "Comma-separated exact sample IDs from dropdown"},
            "sample_id_like": {"type": "string", "description": "Comma-separated keywords for partial sample ID matching, OR-combined"},
            "apply_pathogenic_threshold": {"type": "boolean"},
            "known_loci_only": {"type": "boolean"},
            "known_motifs_only": {"type": "boolean"},
            "mendelian_genes_only": {"type": "boolean"},
            "gene_id": {"type": "string"},
            "gene_symbol": {"type": "string"},
            "locus_id": {"type": "string"},
            "tag": {"type": "string", "description": "Filter to loci with this tag"},
            "motif_size": {"type": "string", "description": "Motif size filter. Comma-separated entries: exact (3), range (3-6), open-end (5-), open-start (-10)"},
        },
        "sort_options": list(VALID_SORTS),
        "total_loci": app.config.get("TOTAL_LOCI", 0),
        "database_columns": app.config.get("DB_COLUMNS", []),
        "all_tags": app.config["ANNOTATIONS"]["all_tags"],
        "sample_ids": app.config.get("SAMPLE_IDS", []),
    })


@app.errorhandler(400)
def bad_request(e):
    return jsonify({"error": "Bad request", "detail": str(e)}), 400


@app.errorhandler(404)
def not_found(e):
    return jsonify({"error": "Not found", "detail": str(e)}), 404


@app.errorhandler(500)
def internal_error(e):
    traceback.print_exc()
    return jsonify({"error": "Internal server error", "detail": str(e)}), 500


@app.errorhandler(Exception)
def handle_exception(e):
    # Let werkzeug HTTPExceptions (405, 403, 413, ...) keep their own status
    # instead of being flattened into a 500 with a stack trace. Only 400/404/500
    # have dedicated handlers above; everything else lands here.
    if isinstance(e, HTTPException):
        return e
    traceback.print_exc()
    return jsonify({"error": "Internal server error", "detail": str(e)}), 500


# ---------------------------------------------------------------------------
# Startup configuration
# ---------------------------------------------------------------------------


def configure_app(db_path, sample_table=None, known_loci_json=None, annotations_db="annotations.db",
                  strchive_loci_json=None):
    """Populate app.config for a single results database.

    This is the shared startup path used by both main() and the test client, so the
    config layout stays in one place. Reference-data inputs (sample table, known-loci
    JSON) are optional; when absent the corresponding enrichments are simply empty.
    """
    db_path = os.path.abspath(db_path)
    app.config["DB_PATH"] = db_path
    # The templates still expect a db_labels list (their source-selector block is
    # guarded by `db_labels|length > 1`, so a single entry hides the selector).
    app.config["DB_LABELS"] = ["TRails"]

    # Sample metadata (optional).
    if sample_table and os.path.exists(sample_table):
        sample_lookup, affected_lookup, analysis_lookup = load_sample_table(sample_table)
    else:
        sample_lookup, affected_lookup, analysis_lookup = {}, {}, {}

    # Known disease loci (optional, hermetic: no network fetch by default). A
    # cached STRchive-loci JSON, when supplied, adds the detail-page fallback.
    if known_loci_json and os.path.exists(known_loci_json):
        disease_trees, strchive_trees, locus_lookup = load_known_disease_loci(
            known_loci_json, fetch_strchive=False,
            strchive_filepath=strchive_loci_json if strchive_loci_json and os.path.exists(strchive_loci_json) else None,
            build_locus_lookup=True,
        )
    else:
        disease_trees, strchive_trees, locus_lookup = {}, {}, {}

    app.config["LOOKUPS"] = {
        "sample": sample_lookup,
        "affected": affected_lookup,
        "analysis": analysis_lookup,
        "disease_trees": disease_trees,
        "locus": locus_lookup,
        "strchive_trees": strchive_trees,
    }

    annotations = load_annotations(annotations_db)
    app.config["ANNOTATIONS"] = annotations
    app.config["ANNOTATIONS_DB_PATH"] = annotations_db

    # Inspect the database: loci columns, total count, present tables.
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        total_loci = conn.execute("SELECT COUNT(*) FROM loci").fetchone()[0]
        cursor = conn.execute("SELECT * FROM loci LIMIT 1")
        db_columns = [desc[0] for desc in cursor.description]
        tables = {row[0] for row in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()}
        # Sample IDs for the filter dropdown come from the sample table when one
        # is supplied; otherwise fall back to the distinct sample_ids in the
        # swim_plot table so the dropdown is still populated (those IDs are what
        # the sample_id_keyword filter matches against).
        if sample_lookup:
            sample_ids = sorted(sample_lookup.keys())
        elif "swim_plot" in tables:
            sample_ids = sorted({
                row[0] for row in conn.execute("SELECT DISTINCT sample_id FROM swim_plot") if row[0]
            })
        else:
            sample_ids = []
    finally:
        conn.close()

    db_columns_set = set(db_columns)
    app.config["TOTAL_LOCI"] = total_loci
    app.config["DB_COLUMNS"] = db_columns
    app.config["DB_COLUMNS_SET"] = db_columns_set
    app.config["DB_TABLES"] = tables
    app.config["HAS_SKINNY"] = {f"sk_{ot}" for ot in OUTLIER_TYPE_MAP.values()}.issubset(tables)
    app.config["HAS_MENDELIAN"] = "mendelian_violations" in tables
    app.config["SAMPLE_IDS"] = sample_ids

    # Known-disease-locus filter set + pathogenic-threshold map (from the loci table's
    # KnownDiseaseLocus column, and the disease catalog for thresholds when available).
    known_disease_locus_ids = set()
    known_disease_locus_thresholds = {}
    if "KnownDiseaseLocus" in db_columns_set:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        try:
            known_disease_locus_ids = {
                row[0] for row in conn.execute(
                    "SELECT LocusId FROM loci WHERE KnownDiseaseLocus IS NOT NULL AND KnownDiseaseLocus != ''"
                )
            }
            if disease_trees and known_disease_locus_ids:
                # Only known-disease loci can contribute a threshold (the overlap +
                # motif test below mirrors how KnownDiseaseLocus was set at build
                # time), so restrict the scan to that id set instead of every locus.
                placeholders = ",".join("?" * len(known_disease_locus_ids))
                for row in conn.execute(
                    f"SELECT LocusId, Chrom, Start0Based, End1Based, Motif FROM loci "
                    f"WHERE LocusId IN ({placeholders})",
                    tuple(known_disease_locus_ids)):
                    chrom_key = row[1].replace("chr", "") if row[1] else ""
                    if chrom_key not in disease_trees:
                        continue
                    overlaps = disease_trees[chrom_key].overlap(row[2], row[3])
                    if not (overlaps and row[4]):
                        continue
                    for interval in overlaps:
                        if compute_jaccard(row[2], row[3], interval.begin, interval.end) <= 0.66:
                            continue
                        disease_data = interval.data
                        disease_motifs = [disease_data.get("RepeatUnit", "")] + (disease_data.get("PathogenicMotifs") or [])
                        if not any(motifs_match(row[4], dm) for dm in disease_motifs):
                            continue
                        for disease in disease_data.get("Diseases", []):
                            pmin = disease.get("PathogenicMin")
                            if pmin is not None and (row[0] not in known_disease_locus_thresholds or pmin < known_disease_locus_thresholds[row[0]]):
                                known_disease_locus_thresholds[row[0]] = pmin
                        break
        finally:
            conn.close()
    app.config["KNOWN_DISEASE_LOCUS_IDS"] = known_disease_locus_ids
    app.config["KNOWN_DISEASE_LOCUS_THRESHOLDS"] = known_disease_locus_thresholds

    # Sample QC + outlier warnings (require swim_plot).
    qc_data = compute_sample_qc_data(db_path)
    if qc_data:
        app.config["SAMPLE_QC_CACHE"] = qc_data
        app.config["OUTLIER_WARNINGS"] = compute_outlier_warnings(qc_data["rank1"])
    else:
        app.config["SAMPLE_QC_CACHE"] = {"rank1": [], "top10": []}
        app.config["OUTLIER_WARNINGS"] = {}

    # Mendelian warnings (from the single DB; empty when the table is absent).
    app.config["MENDELIAN_WARNINGS"] = load_mendelian_warnings(db_path)
    app.config["QC2_CACHE"] = None

    return total_loci, db_columns_set, tables


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main():
    args = parse_args()
    validate_database(args.db)

    print(f"Loading database: {args.db}")
    if args.sample_table:
        print(f"Loading sample table: {args.sample_table}")
    if args.known_loci_json:
        print(f"Loading known disease loci: {args.known_loci_json}")
    print(f"Loading annotations: {args.annotations_db}")

    total_loci, _db_columns_set, tables = configure_app(
        args.db,
        sample_table=args.sample_table,
        known_loci_json=args.known_loci_json,
        annotations_db=args.annotations_db,
        strchive_loci_json=args.strchive_loci_json,
    )

    print(f"Loaded {total_loci:,d} loci")
    print(f"Loaded {len(app.config['LOOKUPS']['sample']):,d} samples")
    print(f"Loaded {len(app.config['ANNOTATIONS']['notes'])} notes, {len(app.config['ANNOTATIONS']['all_tags'])} unique tags")
    for col in ("CanonicalMotif", "IsKnownMotif", "IsInMendelianGene"):
        if col not in app.config["DB_COLUMNS_SET"]:
            print(f"Note: database missing optional column: {col} (using fallback)")
    print(f"Skinny fast-path: {'ON' if app.config['HAS_SKINNY'] else 'OFF (run the skinny-table build step)'}")
    print(f"Mendelian QC page: {'available' if app.config['HAS_MENDELIAN'] else 'hidden (no mendelian_violations table)'}")

    print(f"\nStarting server on {args.host}:{args.port}")
    print(f"Web UI: http://{args.host}:{args.port}/")
    app.run(host=args.host, port=args.port, debug=args.debug, use_reloader=args.debug, threaded=True)


if __name__ == "__main__":
    main()
