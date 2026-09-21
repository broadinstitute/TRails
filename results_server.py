#!/usr/bin/env python3
"""Flask-based read-only API server for TRails tandem-repeat outlier results.

Serves data from a single pre-computed DuckDB database of locus-level outlier
statistics, enriched with sample metadata, gene-disease associations, and known
disease loci. This is the single-database TRails server: there is no source
selector, no read-visualization (readviz) integration, and no cloud dependency.

The server tolerates optional tables being absent. The ``loci`` table is the only
hard requirement (see ``REQUIRED_COLUMNS``); ``swim_plot``, the phenotype-score
tables and the Mendelian-violation tables are all optional, and the corresponding
pages/endpoints degrade gracefully when they are missing.
"""

import argparse
import collections
from datetime import datetime
import json
import math
import os
import re
import sys
import threading
import traceback
import zlib

import flask
from flask import Flask, request, Response
from werkzeug.exceptions import HTTPException
import jinja2
import numpy as np

# The only module that imports duckdb; everything here goes through its sqlite3-shaped API.
import duckdb_compat

# Standalone motif primitive (no str_analysis dependency, no sys.path insertion).
from motif_utilities import COMPLEMENT, compute_canonical_motif

# The number of samples the database was built from, recorded in its metadata table.
from result_database import read_sample_count

# The bases compute_canonical_motif can canonicalize. Deriving this from motif_utilities.COMPLEMENT
# rather than spelling the alphabet out here keeps the two from drifting: anything outside the table
# makes reverse_complement raise KeyError, so the "motif" query parameter is validated against it.
VALID_MOTIF_PATTERN = re.compile(f"^[{''.join(sorted(COMPLEMENT))}]+$", re.IGNORECASE)

# Known-disease-locus matching and affected-status normalization live in the ported
# locus_annotations module. Import them defensively so the server stays importable
# even when that optional module is not present alongside it; the locus-detail page's
# known-disease annotations simply degrade to "no match" in that case.
try:
    from locus_annotations import (
        BLANK_STATUS_VALUES,
        compute_jaccard,
        load_known_disease_loci,
        motifs_match,
        normalize_affected_status_for_logic,
        strchive_locus_motifs,
    )
except ImportError as import_error:
    # A missing module is the supported degraded mode; a missing *name* means the two files have
    # drifted, which would otherwise disable the known-disease annotations without a word.
    print(f"WARNING: falling back to the built-in known-disease stubs: {import_error}")

    BLANK_STATUS_VALUES = {"", "nan", "none", "na", "n/a", "null"}

    def normalize_affected_status_for_logic(value):
        """Lowercase + collapse 'possibly affected' to 'affected' (fallback)."""
        if value is None:
            return None
        if isinstance(value, float) and value != value:  # NaN
            return None
        normalized = str(value).strip().lower()
        # Every blank spelling maps to None, exactly as locus_annotations does it. Checking only
        # "nan" here left a cell reading "NA" classified as the status "na", so the locus-detail
        # page showed "Na" instead of "Unknown" whenever this fallback was in use.
        if normalized in BLANK_STATUS_VALUES:
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

    def strchive_locus_motifs(locus_data):
        """Reference plus pathogenic STRchive motifs of one locus (fallback)."""
        return (list(locus_data.get("reference_motif_reference_orientation") or [])
                + list(locus_data.get("pathogenic_motif_reference_orientation") or []))

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
    "variation_cluster",
)

REQUIRED_TABLE = "loci"
REQUIRED_COLUMNS = {"LocusId", "Chrom", "Start0Based", "End1Based", "Motif", "MotifSize"}

# A coordinate is bound as a 64-bit integer, so a larger one cannot be passed as a query
# parameter. A reference-region start must stay strictly below the limit, since a bare
# position is widened to [start, start + 1).
MAX_REFERENCE_COORDINATE = 2 ** 63 - 1

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
    # The three variation-cluster columns are optional: they are carried over from the input
    # TSV only when it supplies them. A database predating them reports them as None
    # (row_to_list_dict tolerates a missing static column), which leaves the VC column blank.
    "VariationClusterSizeDiff", "VariationCluster", "VariationClusterFilterReason",
    "AoU1027_MaxAllele", "AoU1027_99thPercentile",
    "AoU1027_Stdev", "AoU1027_StdevRankByMotif", "AoU1027_StdevRankTotalNumberByMotif",
    "HPRC256_MaxAllele", "HPRC256_99thPercentile",
    "HPRC256_Stdev", "HPRC256_StdevRankByMotif", "HPRC256_StdevRankTotalNumberByMotif",
    "HPRC256_StdevPercentile",
    "AoU1027_StdevPercentile",
    "TenK10K_MaxAllele", "TenK10K_99thPercentile",
    "TRExplorerLocusId", "TRExplorerSource", "TRExplorerReferenceRepeatPurity",
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


# Integer query parameters are bound as BIGINT, which is signed 64-bit, so binding a larger Python
# int raises at query time and would surface as a 500. Every parsed integer that becomes a query
# parameter goes through parse_int64 so an oversized value fails validation with a 400 like any
# other bad input.
MIN_INT64 = -(2 ** 63)
MAX_INT64 = 2 ** 63 - 1


def escape_like_wildcards(text):
    """Escape SQL ILIKE metacharacters so text matches literally.

    Pair with ``ESCAPE '\\'`` in the ILIKE clause. Used where a value that is conceptually an exact
    match has to go through ILIKE because it is embedded in a larger delimited string.

    Args:
        text: The literal value to match.

    Returns:
        The value with ``\\``, ``%`` and ``_`` backslash-escaped.
    """
    return text.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def parse_int64(text):
    """Parse text as an integer that can be bound as a query parameter.

    Args:
        text: The raw parameter value.

    Returns:
        The parsed integer.

    Raises:
        ValueError: If text is not an integer, or is outside the signed 64-bit range.
    """
    value = int(text)
    if not MIN_INT64 <= value <= MAX_INT64:
        raise ValueError(f"integer out of the signed 64-bit range: '{text}'")
    return value


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
                    params.extend([parse_int64(left), parse_int64(right)])
                except ValueError:
                    return None, [], f"Invalid motif_size range: '{part}'"
            elif left and not right:
                try:
                    params.append(parse_int64(left))
                    clauses.append("MotifSize >= ?")
                except ValueError:
                    return None, [], f"Invalid motif_size value: '{part}'"
            elif not left and right:
                try:
                    params.append(parse_int64(right))
                    clauses.append("MotifSize <= ?")
                except ValueError:
                    return None, [], f"Invalid motif_size value: '{part}'"
            else:
                return None, [], f"Invalid motif_size format: '{part}'"
        else:
            try:
                params.append(parse_int64(part))
                clauses.append("MotifSize = ?")
            except ValueError:
                return None, [], f"Invalid motif_size value: '{part}'"

    if not clauses:
        return None, [], None

    return "(" + " OR ".join(clauses) + ")", params, None


# Shapes the merged search box routes on, checked in this order (see classify_search_term).
# A locus id is "{chrom}-{start}-{end}-{MOTIF}"; a gene id is an Ensembl accession; a region is
# a chromosome with an optional span. Everything else is treated as a gene symbol.
# The motif field accepts the same alphabet the motif filter does (motif_utilities.COMPLEMENT),
 # so a locus id whose motif carries an IUPAC ambiguity code is still recognized as a locus id
 # rather than falling through to a gene-symbol search.
SEARCH_LOCUS_ID_PATTERN = re.compile(
    r"^[0-9A-Za-z._]+-\d+-\d+-[%s%s]+$" % ("".join(sorted(COMPLEMENT)), "".join(sorted(COMPLEMENT)).lower()))
SEARCH_GENE_ID_PATTERN = re.compile(r"^(ENSG\d+)(?:\.\d+)?$", re.IGNORECASE)
# Commas do double duty in the search box: they separate terms, and they are the thousands
# separators parse_reference_region accepts inside a region ("chr16:11,579,459-11,579,529").
# Matching a region-with-span first consumes those internal commas as part of the term, leaving
# every other comma to separate terms. The span is required in this branch, so a bare chromosome
# is not swallowed and "1,2,3" still splits into three terms.
SEARCH_TERM_PATTERN = re.compile(
    r"\s*(?:(?P<region>(?:chr)?(?:\d{1,2}|X|Y|M|MT):[\d,]+(?:-[\d,]+)?)|(?P<other>[^,]+))\s*(?:,|$)",
    re.IGNORECASE)


def split_search_terms(text):
    """Split the search box's value into terms without breaking comma-formatted regions.

    Args:
        text: The raw search string.

    Returns:
        List of stripped, non-empty terms.
    """
    terms = []
    for match in SEARCH_TERM_PATTERN.finditer(text or ""):
        term = (match.group("region") or match.group("other") or "").strip()
        if term:
            terms.append(term)
    return terms
SEARCH_REGION_PATTERN = re.compile(r"^(?:chr)?(?:\d{1,2}|X|Y|M|MT)(?::[\d,]+(?:-[\d,]+)?)?$", re.IGNORECASE)


def normalize_gene_id(text):
    """Return an Ensembl gene id spelled the way the database stores it, else None.

    Stored gene ids are upper case and carry no ".<version>" suffix. Normalizing the queried value
    keeps the comparison a plain equality test rather than a case-insensitive one.

    Args:
        text: A candidate gene id.

    Returns:
        The normalized gene id, or None when text is not shaped like an Ensembl gene id.
    """
    match = SEARCH_GENE_ID_PATTERN.match((text or "").strip())
    return match.group(1).upper() if match else None


def classify_search_term(term):
    """Route one term from the merged search box to the filter it belongs to.

    The search box replaces the separate Gene Symbol(s), Gene ID, Locus Id(s) and Reference
    Region fields, so each term is classified by its shape:

      - "2-89831737-89831752-CCATT"           -> locus id
      - "ENSG00000102081"                     -> gene id
      - "chr16:11579459-11579529", "chr16"    -> reference region
      - anything else, e.g. "FMR1"            -> gene symbol

    Regions are recognized with a strict pattern rather than by handing every term to
    parse_reference_region, which reads any bare alphanumeric string as a whole chromosome and
    would therefore swallow every gene symbol. A term that carries a ":" is the one exception:
    gene symbols never contain one, so it is treated as a region attempt and its parse error is
    reported rather than silently searching for a gene of that name.

    Args:
        term: One stripped, non-empty term from the search box.

    Returns:
        Tuple (kind, value, error_message). kind is "locus_id", "gene_id", "region" or
        "gene_symbol". value is the parsed (chrom, start_0based, end_1based) tuple for a
        region and the term itself otherwise. On a region that will not parse, kind and value
        are None and error_message says why.
    """
    if SEARCH_LOCUS_ID_PATTERN.match(term):
        # The build stores the locus id exactly as the input matrix spelled it (build_database
        # copies the "trid" field verbatim), so a lower-case motif is a spelling the database can
        # really contain. The term is therefore kept as typed and both endpoints compare it
        # case-insensitively rather than normalizing it to a spelling that may not exist.
        return "locus_id", term, None
    normalized_gene_id = normalize_gene_id(term)
    if normalized_gene_id:
        return "gene_id", normalized_gene_id, None
    if SEARCH_REGION_PATTERN.match(term) or ":" in term:
        region, error = parse_reference_region(term)
        if region is None:
            return None, None, error or f"search term '{term}' is not a valid region"
        return "region", region, None
    return "gene_symbol", term, None


def parse_reference_region(region_string):
    """Parse a reference region string into a (chrom, start_0based, end_1based) tuple.

    The "chr" prefix is optional and comma thousands separators are stripped, so the
    lexical form of a region copied from IGV or the UCSC browser is accepted — but note
    that the start is read as 0-based, unlike those browsers' 1-based display. The
    returned chromosome always carries a lowercase "chr" prefix, whatever case was typed.

    - "chr1:100-200": the half-open interval [100, 200), i.e. the same 0-based
      start / 1-based end convention as the Start0Based and End1Based columns.
    - "chr1:100": the single base at 0-based position 100, i.e. [100, 101).
      "chr1:100-100" is read the same way, since a zero-length interval could
      never overlap any locus.
    - "chr1": the whole chromosome, returned with an end of None.

    Args:
        region_string: The filter string from the API parameter.

    Returns:
        Tuple (region, error_message). region is a (chrom, start_0based, end_1based)
        tuple, or None when region_string is blank or could not be parsed (in which
        case error_message describes the problem).
    """
    text = (region_string or "").strip().replace(",", "")
    if not text:
        return None, None

    expected = "expected 'chrom:start-end' with a 0-based start, e.g. chr16:11579459-11579529"
    chrom_text, _, span_text = text.partition(":")
    if not chrom_text or not all(c.isalnum() or c in "._" for c in chrom_text):
        return None, f"reference_region '{region_string}' has an invalid chromosome — {expected}"
    # Rebuild the prefix rather than keeping the typed one, so "Chr1" and "CHR1" both
    # come out as "chr1" and match the stored, case-sensitive chromosome names.
    suffix = chrom_text[len("chr"):] if chrom_text.lower().startswith("chr") else chrom_text
    chrom = f"chr{suffix}"

    span_text = span_text.strip()
    if not span_text:
        return (chrom, 0, None), None

    start_text, dash, end_text = span_text.partition("-")
    # isdecimal() rather than isdigit(): the latter also accepts characters like the
    # superscript "²", which int() then rejects with a ValueError. The end is required
    # whenever a dash was typed, so a truncated "chr1:100-" is an error rather than
    # quietly collapsing to the bare-position reading.
    if not start_text.isdecimal() or (dash and not end_text.isdecimal()):
        return None, f"reference_region '{region_string}' could not be parsed — {expected}"
    start = int(start_text)
    end = int(end_text) if dash else start
    if end < start:
        return None, f"reference_region '{region_string}' ends before it starts — {expected}"
    if start >= MAX_REFERENCE_COORDINATE or end > MAX_REFERENCE_COORDINATE:
        return None, (f"reference_region '{region_string}' has a coordinate above "
                      f"{MAX_REFERENCE_COORDINATE} — {expected}")

    return (chrom, start, max(end, start + 1)), None


def chromosome_name_variants(chrom):
    """Return the plausible spellings of a chromosome name, for matching a database.

    Databases differ in whether they store the "chr" prefix, users differ in how they
    capitalize the suffix, and the mitochondrial chromosome is written as either M or MT.
    Listing the alternatives as an IN set keeps the lookup an equality test, which DuckDB can
    evaluate against the Chrom column's per-row-group min/max; a case-insensitive comparison
    would match more spellings but would have to be evaluated row by row.

    Args:
        chrom: A chromosome name carrying the lowercase "chr" prefix.

    Returns:
        A sorted list of candidate names, each with and without the "chr" prefix.
    """
    suffix = chrom[len("chr"):]
    suffixes = {suffix, suffix.upper()}
    if suffix.upper() in ("M", "MT"):
        suffixes.update(("M", "MT"))
    return sorted(suffixes | {f"chr{s}" for s in suffixes})


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
        help="Path to the single TRails DuckDB results database.",
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
        "--allow-remote-writes",
        action="store_true",
        help="Allow notes and tags to be edited from other machines. Only meaningful with a "
             "non-loopback --host; the annotation routes have no authentication, so anyone who "
             "can reach the port could change them.",
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
        default="annotations.duckdb",
        help="Path to annotations DuckDB database for notes/tags (default: annotations.duckdb).",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run Flask in debug mode (auto-reloader + interactive debugger). Off by default; "
             "enabling it re-runs the full startup scan in the reloader child and exposes the debugger.",
    )
    return parser.parse_args()


def print_leftover_write_files_help(path, error, leftover):
    """Explain a read-only open that failed because of a leftover write-ahead log.

    Args:
        path: Path to the DuckDB database file.
        error: The DuckDB error that was raised.
        leftover: List of (filename, size in bytes) tuples from
            duckdb_compat.find_leftover_write_files.
    """
    abspath = os.path.abspath(path)
    print(f"Error: database could not be read: {error}")
    print(f"  This is almost certainly because of these leftover files next to {path}:")
    for filename, size in leftover:
        print(f"    {os.path.basename(filename)}  ({size:,d} bytes)")
    print("  They are left behind when a script writing to the database is killed part way through.")
    print("  DuckDB has to replay the write-ahead log into the database file before it will serve any")
    print("  read, and it can only do that with write access. This server opens databases read-only,")
    print("  so it cannot replay the log and even a SELECT fails.")
    print("  To recover, first check that nothing is still writing to the database:")
    print(f"    lsof '{abspath}'")
    print("  If that prints nothing, open it once read-write so DuckDB replays the log itself:")
    print(f"    python3 -c \"import duckdb; duckdb.connect('{abspath}').execute('SELECT 1')\"")
    print("  Whatever the interrupted write had already committed is kept and the rest is discarded,")
    print("  so whichever script was writing to the database has to be re-run.")


def validate_database(path):
    """Validate that the database file exists and has the required table/columns.

    Args:
        path: Path to the DuckDB database file.
    """
    if not os.path.exists(path):
        print(f"Error: database file not found: {path}")
        sys.exit(1)
    try:
        conn = duckdb_compat.connect(path, read_only=True)
    except duckdb_compat.OperationalError as e:
        leftover = duckdb_compat.find_leftover_write_files(path)
        if not leftover:
            raise
        print_leftover_write_files_help(path, e, leftover)
        sys.exit(1)
    try:
        tables = sorted(duckdb_compat.list_tables(conn))
        if REQUIRED_TABLE not in tables:
            print(f"Error: database missing required table '{REQUIRED_TABLE}'. Found tables: {tables}")
            sys.exit(1)
        missing = REQUIRED_COLUMNS - duckdb_compat.table_columns(conn, REQUIRED_TABLE)
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
    # Every column input_tables.read_sample_metadata canonicalizes for the build (input_tables.py's
    # match_columns list), so a "FamilyID"/"Sex" header reaches parse_outlier_samples under the
    # lower-case name it looks up instead of coming back as None.
    for canonical in ("sex", "family_id", "maternal_id", "paternal_id",
                      "affected_status", "analysis_status", "phenotype_description"):
        actual = column_by_normalized.get(_normalize_metadata_column(canonical))
        if actual and actual != canonical:
            rename[actual] = canonical
    df = df.rename(columns=rename)

    analysis_status_remap = {"rncc": "unsolved", "rcpc": "unsolved", "s_kgfp": "solved"}
    if "analysis_status" in df.columns:
        df["analysis_status"] = (
            df["analysis_status"].astype(str).str.strip().str.lower().replace(analysis_status_remap)
        )
        # The build maps every blank spelling to "unknown" (input_tables.read_sample_metadata);
        # covering only "" and "nan" here would leave a plain "NA" cell rendering as "Na".
        df["analysis_status"] = df["analysis_status"].apply(
            lambda v: "unknown" if v in BLANK_STATUS_VALUES else v)
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
        db_path: Path to the results DuckDB database.
        threshold: Violation-rate threshold (default 0.10 = 10%).

    Returns:
        dict mapping sample_id to warning info for samples exceeding the threshold
        on autosome OR chrX.
    """
    if not os.path.exists(db_path):
        return {}

    conn = duckdb_compat.connect(db_path, read_only=True)
    conn.row_factory = duckdb_compat.Row
    try:
        if "mendelian_violations" not in duckdb_compat.list_tables(conn):
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

    The rank1 rows are a subset of the top10 rows, so both counts come from one pass over
    outlier_rank <= 10 rather than a query each, which halves the work on a swim_plot table
    of tens of millions of rows.

    Returns:
        dict with {"rank1": [...], "top10": [...]}, or None if swim_plot / its
        outlier_rank column is missing.
    """
    conn = duckdb_compat.connect(db_path, read_only=True)
    conn.row_factory = duckdb_compat.Row
    try:
        if "swim_plot" not in duckdb_compat.list_tables(conn):
            return None
        if "outlier_rank" not in duckdb_compat.table_columns(conn, "swim_plot"):
            return None

        bin_case = """
            CASE
                WHEN MotifSize IS NULL OR MotifSize <= 0 THEN 'Unknown'
                WHEN MotifSize = 1 THEN '1bp'
                WHEN MotifSize = 2 THEN '2bp'
                WHEN MotifSize >= 3 AND MotifSize <= 6 THEN '3-6bp'
                WHEN MotifSize >= 7 AND MotifSize <= 24 THEN '7-24bp'
                ELSE '25+bp'
            END
        """
        # Both counts are counts of LOCI, which is what the QC page labels them. swim_plot holds
        # one row per outlier entry, and a diploid sample is recorded under both of its alleles, so
        # a sample can own two of a locus's first ten entries; COUNT(*) would report that locus
        # twice. Counting distinct LocusId is also what makes rank1_count independent of whether
        # two entries at one locus ever share rank 1.
        rows = conn.execute(
            f"SELECT sample_id, {bin_case} AS bin, "
            "COUNT(DISTINCT CASE WHEN outlier_rank = 1 THEN LocusId END) AS rank1_count, "
            "COUNT(DISTINCT LocusId) AS top10_count FROM swim_plot "
            "WHERE outlier_type = 'AllAlleles' AND outlier_rank <= 10 "
            "GROUP BY sample_id, bin ORDER BY sample_id, bin"
        ).fetchall()
        # A sample/bin with no rank-1 locus is dropped from rank1, which is what querying
        # outlier_rank = 1 on its own did; it stays in top10 with its full count.
        return {
            "rank1": [{"sample_id": r["sample_id"], "bin": r["bin"], "count": r["rank1_count"]}
                      for r in rows if r["rank1_count"]],
            "top10": [{"sample_id": r["sample_id"], "bin": r["bin"], "count": r["top10_count"]} for r in rows],
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
    """Load user annotations (notes/tags), creating the schema if needed.

    The timestamp columns carry no DEFAULT: every route that writes a row supplies
    created_at and updated_at itself.

    Returns:
        dict with keys: notes, tags, all_tags, tag_to_loci.
    """
    conn = duckdb_compat.connect(db_path)
    conn.execute("""CREATE TABLE IF NOT EXISTS notes (
        locus_id TEXT PRIMARY KEY,
        note_text TEXT NOT NULL,
        created_at TEXT NOT NULL,
        updated_at TEXT NOT NULL
    )""")
    conn.execute("""CREATE TABLE IF NOT EXISTS tags (
        locus_id TEXT NOT NULL,
        tag TEXT NOT NULL,
        created_at TEXT NOT NULL,
        PRIMARY KEY (locus_id, tag)
    )""")
    conn.commit()

    notes = {}
    for row in conn.execute("SELECT locus_id, note_text, updated_at FROM notes ORDER BY locus_id"):
        notes[row[0]] = {"note_text": row[1], "updated_at": row[2]}

    tags = {}
    tag_to_loci = collections.defaultdict(set)
    for row in conn.execute("SELECT locus_id, tag FROM tags ORDER BY locus_id, tag"):
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
    """Get a read-only DuckDB connection to the single results database."""
    conn = duckdb_compat.connect(app.config["DB_PATH"], read_only=True)
    conn.row_factory = duckdb_compat.Row
    return conn


# ---------------------------------------------------------------------------
# CORS handler
# ---------------------------------------------------------------------------


# Loopback addresses. A server bound to one of these is reachable only from this machine, which
# is what makes the unauthenticated annotation routes acceptable by default.
LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}


@app.before_request
def restrict_writes():
    """Reject state-changing requests that this server should not accept.

    Two separate exposures, both of which matter because the annotation routes (notes and tags)
    have no authentication:

      * Another origin. A POST with a plain content type is a "simple" request the browser sends
        before CORS gets a chance to hide anything, so the response headers below cannot stop
        another page the user has open from writing to this database.
      * Another machine. Bound to a non-loopback address the write routes are reachable by anyone
        who can route to the port, so they are refused unless --allow-remote-writes says otherwise.

    Read requests are never affected, and a local script with no Origin header still works.
    """
    if request.method in ("GET", "HEAD", "OPTIONS"):
        return None
    if not app.config.get("ALLOW_REMOTE_WRITES", False) \
            and app.config.get("BIND_HOST") not in LOOPBACK_HOSTS \
            and request.remote_addr not in LOOPBACK_HOSTS:
        return jsonify({
            "error": "annotation changes are limited to this machine",
            "detail": "This server is bound to a non-loopback address and has no authentication, "
                      "so notes and tags cannot be edited remotely. Restart it with "
                      "--allow-remote-writes if every client that can reach this port is trusted.",
        }), 403
    origin = request.headers.get("Origin")
    if origin and origin.rstrip("/") != request.host_url.rstrip("/"):
        return jsonify({"error": "cross-origin requests are not allowed"}), 403
    return None


@app.after_request
def add_cors_headers(response):
    """Allow cross-origin reads of the API, but never cross-origin writes.

    The bundled UI is same-origin (index_page_template.html builds API_BASE_URL from
    window.location.origin), so it needs no CORS headers at all. The wildcard exists only so a
    notebook or a separate viewer can GET the read-only endpoints; advertising PUT/POST/DELETE
    here would invite any page to mutate this database's notes and tags.
    """
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "GET, OPTIONS"
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
    # The full IUPAC table, not a local ACGT-only map: validate_params accepts every base in
    # COMPLEMENT, so an ambiguity code such as Y would otherwise be left uncomplemented here and
    # this fallback would disagree with the CanonicalMotif path it stands in for.
    rc = "".join(COMPLEMENT.get(base, base) for base in reversed(motif))
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
        params["page"] = parse_int64(page_raw)
        if params["page"] < 1:
            raise ValueError
    except ValueError:
        errors.append(f"page must be a positive integer, got '{page_raw}'")

    # page_size (default: 50, max: 500)
    page_size_raw = request.args.get("page_size", "50")
    try:
        params["page_size"] = parse_int64(page_size_raw)
        if params["page_size"] < 1 or params["page_size"] > 500:
            raise ValueError
    except ValueError:
        errors.append(f"page_size must be an integer between 1 and 500, got '{page_size_raw}'")

    # The derived OFFSET is bound, not `page`, so both have to be in range: page=2**62 with
    # page_size=50 passes the checks above and then overflows when bound.
    if "page" in params and "page_size" in params:
        if (params["page"] - 1) * params["page_size"] > MAX_INT64:
            errors.append(f"page is too large for page_size={params['page_size']}, "
                          f"got '{page_raw}'")

    if request.args.get("min_expansion"):
        try:
            params["min_expansion"] = parse_int64(request.args["min_expansion"])
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
        # compute_canonical_motif raises KeyError on any base outside the IUPAC table, which would
        # surface as a 500. Reject the value here so a typo gets the same 400 as every other filter.
        invalid_motifs = [m for m in motif_list if not VALID_MOTIF_PATTERN.match(m)]
        if invalid_motifs:
            errors.append(
                f"motif must contain only IUPAC bases ({''.join(sorted(COMPLEMENT))}); got "
                f"'{invalid_motifs[0]}'")
        elif motif_list:
            params["motif"] = motif_list

    if request.args.get("min_repeats_threshold"):
        try:
            params["min_repeats_threshold"] = parse_int64(request.args["min_repeats_threshold"])
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

    for flag_name in ("apply_pathogenic_threshold", "known_loci_only", "known_motifs_only", "mendelian_genes_only",
                      "exclude_vc_depth_filtered", "exclude_vc_size_filtered"):
        if request.args.get(flag_name, "").lower() in ("true", "1", "yes"):
            params[flag_name] = True

    # These list-valued filters are stored only when at least one entry survives the strip, the
    # way the motif filter above is. A value that is all separators and blanks (",", " , ") is
    # truthy as a request string but parses to nothing, and build_api_query would then emit
    # "IN ()" or an empty "()", which DuckDB rejects as a syntax error. An empty list is the
    # same as no filter at all, so drop it.
    phenotype_keywords = [kw.strip() for kw in request.args.get("phenotype_keyword", "").split(",") if kw.strip()]
    if phenotype_keywords:
        params["phenotype_keyword"] = phenotype_keywords

    sample_id_keywords = [kw.strip() for kw in request.args.get("sample_id_keyword", "").split(",") if kw.strip()]
    if sample_id_keywords:
        params["sample_id_keyword"] = sample_id_keywords

    sample_id_like_keywords = [kw.strip() for kw in request.args.get("sample_id_like", "").split(",") if kw.strip()]
    if sample_id_like_keywords:
        params["sample_id_like"] = sample_id_like_keywords

    # search (optional, comma-separated; each term routed by shape, all OR'ed together). The
    # separate locus_id / reference_region / gene_symbol / gene_id parameters below still work
    # for anyone calling the API directly; the UI sends this one instead.
    if request.args.get("search"):
        search_terms = []
        for term in split_search_terms(request.args["search"]):
            kind, value, error = classify_search_term(term)
            if error:
                errors.append(error)
            else:
                search_terms.append((kind, value))
        if search_terms:
            # Keep the user-facing string for filters_applied; the classified terms below are
            # a SQL internal consumed by build_api_query (excluded from filters_applied).
            params["search"] = request.args["search"].strip()
            params["search_terms"] = search_terms

    # max_variation_cluster_size_diff (optional, bp): keeps loci whose TRExplorer variation
    # cluster is at most this much larger than the repeat itself.
    if request.args.get("max_variation_cluster_size_diff"):
        try:
            params["max_variation_cluster_size_diff"] = parse_int64(request.args["max_variation_cluster_size_diff"])
        except ValueError:
            errors.append("max_variation_cluster_size_diff must be an integer number of base pairs, "
                          f"got '{request.args['max_variation_cluster_size_diff']}'")

    locus_ids = [lid.strip() for lid in request.args.get("locus_id", "").split(",") if lid.strip()]
    if locus_ids:
        params["locus_id"] = locus_ids

    if request.args.get("reference_region"):
        region, error = parse_reference_region(request.args["reference_region"])
        if error:
            errors.append(error)
        elif region:
            # Keep the user-facing string for filters_applied; the parsed tuple below is
            # a SQL internal consumed by build_api_query (excluded from filters_applied).
            params["reference_region"] = request.args["reference_region"].strip()
            params["reference_region_parsed"] = region

    if request.args.get("chrom"):
        params["chrom"] = request.args["chrom"]
    if request.args.get("gene_id"):
        # Normalized like a search-box term, so ?gene_id=ensg00000102081 and a versioned
        # ?gene_id=ENSG00000102081.5 match the stored spelling instead of silently returning
        # nothing. A value that is not an Ensembl gene id is passed through unchanged.
        params["gene_id"] = normalize_gene_id(request.args["gene_id"]) or request.args["gene_id"].strip()
    gene_symbols = [gs.strip() for gs in request.args.get("gene_symbol", "").split(",") if gs.strip()]
    if gene_symbols:
        params["gene_symbol"] = gene_symbols

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


def build_api_order_by(params, source_columns):
    """Build the ORDER BY clause from API parameters.

    Args:
        params: dict of validated query parameters.
        source_columns: this database's loci column set. Most sort keys read an optional
            annotation column (the sigma percentiles, the variation-cluster size, the
            phenotype scores), and a build whose input TSV did not supply one has no such
            column at all, so a sort key naming it is skipped rather than raising "no such
            column" for the whole search. An empty set means "not introspected", and every
            key is kept.

    Returns:
        The ORDER BY clause, as a string beginning with " ORDER BY".
    """
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
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
        # Widest variation cluster first. DESC also puts the loci annotated with no variation
        # cluster (NULL) last.
        "variation_cluster": ("VariationClusterSizeDiff", "DESC"),
    }

    # ORDER BY expression for the sort keys whose column cannot be ordered by as stored.
    # VariationClusterSizeDiff is VARCHAR whenever the build saw a value that is not an integer
    # (see result_database.infer_column_types), and ordering a VARCHAR by itself would compare
    # digit strings ("999" above "4847"). The cast is harmless on the BIGINT case. TRY_CAST
    # rather than CAST because DuckDB raises on a value that will not parse, where SQLite
    # silently yielded 0; a value that will not parse sorts as NULL, which DESC puts last
    # alongside the loci that have no variation cluster at all. The SORT_MAPPING column name
    # stays the real one so the missing-column check below still tests a column this database
    # either has or does not.
    SORT_EXPRESSIONS = {"variation_cluster": "TRY_CAST(VariationClusterSizeDiff AS BIGINT)"}

    sequence = list(sort_list)
    for tiebreak in ("count", "region", "size"):
        if tiebreak not in sequence:
            sequence.append(tiebreak)

    order_parts = []
    for field in sequence:
        col, direction = SORT_MAPPING[field]
        if source_columns and col not in source_columns:
            continue
        order_parts.append(f"{SORT_EXPRESSIONS.get(field, col)} {direction}")

    # LocusId last, so the order is a total one. Ties on every key above are common (count,
    # gene_region_rank and allele size are all small integers), and DuckDB scans with several
    # threads, so without a final unique key the tied rows come back in whatever order the
    # threads finish in. That is not cosmetic: paging through a result whose order is unstable
    # can show the same locus on two pages and never show another.
    if not order_parts:
        # Every mapped sort column is absent from this database. LocusId is in REQUIRED_COLUMNS, so
        # it is the one ordering that is always valid; naming a sort column here instead would
        # re-reference the very column the loop just skipped and raise "no such column".
        return " ORDER BY LocusId ASC"
    return " ORDER BY " + ", ".join(order_parts) + ", LocusId ASC"


# The VariationClusterFilterReason value behind each "hide these loci" flag, and the icon it
# draws in the results table's VC column. EXTENSION is the size filter, which is the name the UI
# uses, since the raw value says nothing about what was measured.
VC_FILTER_REASON_FLAGS = (
    ("exclude_vc_depth_filtered", "DEPTH"),       # gold crossed-out circle
    ("exclude_vc_size_filtered", "EXTENSION"),    # dark red crossed-out circle
)


def reference_region_clause(parsed_region):
    """Build the SQL that keeps the loci overlapping one reference region.

    Shared by the reference_region filter and by any region typed into the merged search box,
    so both get the same bounds.

    Args:
        parsed_region: (chrom, start_0based, end_1based) from parse_reference_region. The end
            is None for a whole chromosome.

    Returns:
        Tuple (clause, params) to append to a WHERE clause.
    """
    chrom, region_start, region_end = parsed_region
    # Chrom is matched against both naming conventions ("chr1" and "1") so the filter
    # works whichever one the database uses; either way it is an equality test DuckDB can
    # evaluate against the Chrom column's per-row-group min/max.
    chrom_names = chromosome_name_variants(chrom)
    region_clauses = [f"loci.Chrom IN ({','.join('?' * len(chrom_names))})"]
    region_params = list(chrom_names)
    if region_end is not None:
        # Half-open overlap: the locus starts before the region ends and ends after
        # the region starts. Both comparisons also drop a locus whose coordinate is NULL,
        # which is right for a coordinate filter.
        region_clauses.append("loci.Start0Based < ?")
        region_params.append(region_end)
        region_clauses.append("loci.End1Based > ?")
        region_params.append(region_start)
    else:
        # A whole chromosome has no coordinate comparison to drop the loci with no
        # Start0Based, so say so outright.
        region_clauses.append("loci.Start0Based IS NOT NULL")
    return " AND ".join(region_clauses), region_params


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

    # Guarded like every other filter below: a database built without this outlier type simply has
    # no rows for it, which is not the same thing as a SQL error. OUTLIER_TYPE_MAP advertises all
    # three types regardless of what the loaded database actually contains.
    if not source_columns or f"OutlierSampleIds_{ot}" in source_columns:
        clauses.append(f"(loci.OutlierSampleIds_{ot} IS NOT NULL AND loci.OutlierSampleIds_{ot} != '')")
    else:
        clauses.append("1=0")

    if "min_expansion" in params:
        clauses.append(f"loci.FirstAffectedAlleleSize_{ot} >= loci.NumRepeatsInReference + ?")
        sql_params.append(params["min_expansion"])

    if "require_above_unaffected" in params:
        col_prefix = {"first": "First", "second": "Second", "third": "Third"}[params["require_above_unaffected"]]
        if not source_columns or (f"{col_prefix}AffectedAlleleSize_{ot}" in source_columns and f"FirstUnaffectedAlleleSize_{ot}" in source_columns):
            # NULL means "this locus has no unaffected sample". 0 does not: it is a real allele
            # size (a full deletion of the repeat), and treating it as missing would pass every
            # affected allele at that locus without comparing it to anything.
            clauses.append(f"(loci.{col_prefix}AffectedAlleleSize_{ot} > loci.FirstUnaffectedAlleleSize_{ot} + ? OR loci.FirstUnaffectedAlleleSize_{ot} IS NULL AND loci.{col_prefix}AffectedAlleleSize_{ot} IS NOT NULL)")
            sql_params.append(min_expansion)
        else:
            clauses.append("1=0")

    if "require_families_above_unaffected" in params:
        col_prefix = {"first": "First", "second": "Second", "third": "Third"}[params["require_families_above_unaffected"]]
        by_family_col = f"loci.{col_prefix}AffectedAlleleSize_{ot}_ByFamily"
        if not source_columns or f"{col_prefix}AffectedAlleleSize_{ot}_ByFamily" in source_columns:
            # Same as above: only NULL means "no unaffected baseline"; 0 is a real allele size.
            clauses.append(f"({by_family_col} > loci.FirstUnaffectedAlleleSize_{ot} + ? OR loci.FirstUnaffectedAlleleSize_{ot} IS NULL AND {by_family_col} IS NOT NULL)")
            sql_params.append(min_expansion)
        else:
            # Backing column absent: the filter cannot be evaluated, so match no
            # rows rather than silently returning every locus (mirrors the
            # require_above_population "1=0" fallback below).
            clauses.append("1=0")

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
        if not source_columns or f"FirstAffectedAlleleSize_{ot}" in source_columns:
            clauses.append(f"loci.FirstAffectedAlleleSize_{ot} >= ?")
            sql_params.append(params["min_repeats_threshold"])
        else:
            clauses.append("1=0")

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
        # NULL NOT IN (...) is unknown in SQL, which would drop unannotated
        # (NULL gene_region) loci; keep them since they are not in the excluded set.
        clauses.append(f"(loci.gene_region IS NULL OR loci.gene_region NOT IN ({','.join('?' * len(db_regions))}))")
        sql_params.extend(db_regions)

    if "min_pli" in params:
        if not source_columns or "pLI" in source_columns:
            clauses.append("pLI >= ?")
            sql_params.append(params["min_pli"])
        else:
            clauses.append("1=0")

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
        else:
            # Neither sigma-percentile column is present: match no rows instead of
            # silently ignoring the filter.
            clauses.append("1=0")

    if "phenotype_keyword" in params:
        # A keyword selects a locus when ANY of its outlier samples carries that phenotype. The
        # per-outlier phenotypes live in swim_plot, one row per outlier entry, so the filter is
        # answered from there whenever the database has that table: the loci table's
        # First/Second/ThirdAffectedPhenotype columns are a three-deep display summary of the
        # affected outliers only (analysis_columns fills exactly three ranks), so matching them
        # instead would drop a locus whose only matching phenotype belongs to its fourth affected
        # outlier or to an unaffected one. The swim plot matches phenotype_description on every
        # outlier row, and this is the same predicate over the same rows, so the two views select
        # the same loci. Without a swim_plot table there is no swim plot to disagree with, and the
        # summary columns are the only phenotypes the database has, so they are used instead.
        if app.config.get("HAS_SWIM_PLOT_PHENOTYPES"):
            kw_clauses = ["phenotype_description ILIKE ?" for _ in params["phenotype_keyword"]]
            clauses.append("LocusId IN (SELECT LocusId FROM swim_plot WHERE outlier_type = ? "
                           f"AND ({' OR '.join(kw_clauses)}))")
            sql_params.append(ot)
            sql_params.extend(f"%{kw}%" for kw in params["phenotype_keyword"])
        else:
            phenotype_columns = [
                f"{col_prefix}AffectedPhenotype_{ot}" for col_prefix in ("First", "Second", "Third")
                if not source_columns or f"{col_prefix}AffectedPhenotype_{ot}" in source_columns
            ]
            if phenotype_columns:
                kw_clauses = ["(" + " OR ".join(f"{column} ILIKE ?" for column in phenotype_columns) + ")"
                              for _ in params["phenotype_keyword"]]
                clauses.append(f"({' OR '.join(kw_clauses)})")
                for kw in params["phenotype_keyword"]:
                    sql_params.extend([f"%{kw}%"] * len(phenotype_columns))
            else:
                # No phenotype column at all: match no rows rather than ignoring the filter.
                clauses.append("1=0")

    if "sample_id_keyword" in params:
        # These are exact sample ids picked from a dropdown, matched by ILIKE only because the ids
        # are packed into one "{n}x:{sample_id},..." string. Sample ids routinely contain "_",
        # which is ILIKE's single-character wildcard, so escape the metacharacters; otherwise
        # "sample_1" would also match the unrelated "sampleX1".
        kw_clauses = [
            f"((',' || OutlierSampleIds_{ot} || ',') ILIKE ? ESCAPE '\\'"
            f" OR OutlierSampleIds_{ot} ILIKE ? ESCAPE '\\')"
            for _ in params["sample_id_keyword"]
        ]
        clauses.append(f"({' OR '.join(kw_clauses)})")
        for kw in params["sample_id_keyword"]:
            escaped = escape_like_wildcards(kw)
            sql_params.extend([f"%:{escaped},%", f"%:{escaped}:%"])

    if "sample_id_like" in params:
        # ILIKE against the whole packed string would also match allele sizes and the "x:"
        # delimiters ("15" matching "15x:sampleA"), so split it into entries and match inside each
        # entry's sample-id field only, which is its second ":"-separated field. Wildcards are
        # escaped for the same reason as in sample_id_keyword above: a keyword such as "sample_1"
        # means a literal underscore, not ILIKE's single-character wildcard.
        kw_clauses = [
            f"len(list_filter(string_split(OutlierSampleIds_{ot}, ','),"
            f" entry -> split_part(entry, ':', 2) ILIKE ? ESCAPE '\\')) > 0"
            for _ in params["sample_id_like"]
        ]
        clauses.append(f"({' OR '.join(kw_clauses)})")
        sql_params.extend([f"%{escape_like_wildcards(kw)}%" for kw in params["sample_id_like"]])

    if "locus_id" in params:
        clauses.append(f"LocusId IN ({','.join('?' * len(params['locus_id']))})")
        sql_params.extend(params["locus_id"])

    if "reference_region_parsed" in params:
        clause, region_params = reference_region_clause(params["reference_region_parsed"])
        clauses.append(clause)
        sql_params.extend(region_params)

    # Merged search box: locus ids, gene ids, gene symbols and reference regions in one field,
    # matching a locus that satisfies any one of them. A term whose column the database lacks is
    # dropped from the OR, the same way the separate gene filters below degrade; if that leaves
    # no terms at all the filter matches nothing rather than everything.
    if "search_terms" in params:
        search_clauses = []
        search_params = []
        locus_id_values = []
        for kind, value in params["search_terms"]:
            if kind == "locus_id":
                # Collected into the one set lookup built below instead of a predicate per term.
                locus_id_values.append(value.lower())
            elif kind == "gene_id":
                if not source_columns or "gene_id" in source_columns:
                    search_clauses.append("gene_id = ?")
                    search_params.append(value)
            elif kind == "gene_symbol":
                if "GeneTableGeneSymbol" in app.config["DB_COLUMNS_SET"]:
                    search_clauses.append("GeneTableGeneSymbol ILIKE ?")
                    search_params.append(f"%{value}%")
            else:
                clause, region_params = reference_region_clause(value)
                search_clauses.append(f"({clause})")
                search_params.extend(region_params)
        if locus_id_values:
            # One "lower(LocusId) IN (...)" for every pasted id rather than one ILIKE each: a
            # chain of N ILIKEs costs N comparisons per row, while the IN list is a single hash
            # probe DuckDB can also answer from the column's statistics. Pasting a few hundred
            # ids is normal here, so the difference is seconds per request.
            # Both sides are lower-cased because the build stores a locus id as the input matrix
            # spelled it, which is what the ILIKE was there for. An IN list holds no wildcards,
            # so a contig name carrying "_" needs no escaping and matches literally.
            search_clauses.insert(0, f"lower(LocusId) IN ({','.join('?' * len(locus_id_values))})")
            search_params[:0] = locus_id_values
        clauses.append(f"({' OR '.join(search_clauses)})" if search_clauses else "0")
        sql_params.extend(search_params)

    # Max variation cluster size diff, compared numerically. The column's storage type depends on
    # the input: result_database.infer_column_types declares it BIGINT when every value the build
    # saw was an integer, and VARCHAR when any was not. TRY_CAST covers both, yielding NULL for a
    # blank or non-numeric value. Comparing the column to '' as well would raise on the BIGINT
    # case ("Could not convert string '' to INT64"), and is unnecessary: COALESCE sends a blank
    # to 0, which passes any maximum. A locus with no variation cluster has no variation beyond
    # the repeat itself, so it passes any maximum.
    if "max_variation_cluster_size_diff" in params:
        # The clause deliberately passes loci whose cluster size is unknown, so a database with no
        # VariationClusterSizeDiff column at all is the same situation for every locus: add no
        # clause rather than the "1=0" the filters below use, which would contradict this filter's
        # own treatment of missing data and return an empty result set.
        if not source_columns or "VariationClusterSizeDiff" in source_columns:
            clauses.append("(loci.VariationClusterSizeDiff IS NULL "
                           "OR COALESCE(TRY_CAST(loci.VariationClusterSizeDiff AS BIGINT), 0) <= ?)")
            sql_params.append(params["max_variation_cluster_size_diff"])

    # Filter: drop loci whose variation cluster was computed and then filtered out, which the
    # results table marks with a crossed-out circle in the VC column. A locus that was never
    # filtered has no reason recorded, so the NULL check keeps it. These are exclusions, so a
    # database with no VariationClusterFilterReason column recorded no locus as filtered out and
    # there is nothing to exclude: add no clause, the way the maximum above does. A "1=0" here
    # would empty the whole result set the moment either box is ticked.
    for flag_name, filter_reason in VC_FILTER_REASON_FLAGS:
        if not params.get(flag_name):
            continue
        if not source_columns or "VariationClusterFilterReason" in source_columns:
            clauses.append("(loci.VariationClusterFilterReason IS NULL "
                           "OR loci.VariationClusterFilterReason != ?)")
            sql_params.append(filter_reason)

    if "chrom" in params:
        clauses.append("Chrom = ?")
        sql_params.append(params["chrom"])

    if "gene_id" in params:
        if not source_columns or "gene_id" in source_columns:
            clauses.append("gene_id = ?")
            sql_params.append(params["gene_id"])
        else:
            clauses.append("1=0")

    if "gene_symbol" in params:
        if "GeneTableGeneSymbol" in app.config["DB_COLUMNS_SET"]:
            gs_clauses = ["GeneTableGeneSymbol ILIKE ?" for _ in params["gene_symbol"]]
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
        # The build writes "" (not NULL) for a locus with no hemizygous alleles
        # (allele_histograms.convert_counts_to_histogram_string returns "" for an empty count
        # dict), so IS NOT NULL alone matches every row and filters nothing. Guarded on the column
        # like the filters above, so an autosome-only database matches nothing instead of raising.
        if not source_columns or "HemizygousAlleleHistogram" in source_columns:
            clauses.append("HemizygousAlleleHistogram IS NOT NULL AND HemizygousAlleleHistogram != ''")
        else:
            clauses.append("1=0")

    where = " AND ".join(clauses) if clauses else "1=1"
    motif_count_col = "CanonicalMotif" if "CanonicalMotif" in source_columns else "Motif"
    # gene_region is an optional annotation column. Without it the region summary simply has
    # nothing to report; querying for it anyway would fail the whole /api/v1/loci request.
    has_gene_region = not source_columns or "gene_region" in source_columns

    order_by = build_api_order_by(params, source_columns)

    # One pass over loci, with no inner "pick the LocusIds first" subquery in front of it.
    # DuckDB is columnar: the filter and the sort read only the columns they name, and the
    # remaining ~180 are read only for the page's rows, so deferring the wide read by hand
    # would just make the table be scanned twice.
    select_query = f"SELECT * FROM loci WHERE {where}{order_by} LIMIT ? OFFSET ?"
    count_query = f"SELECT COUNT(*) FROM loci WHERE {where}"
    # The motif name breaks ties on count, so two motifs with the same number of loci keep a
    # fixed order in the response instead of swapping between identical requests.
    motif_count_query = (f"SELECT {motif_count_col}, COUNT(*) as count FROM loci WHERE {where} "
                         f"GROUP BY {motif_count_col} ORDER BY count DESC, {motif_count_col}")
    gene_region_count_query = (
        f"SELECT gene_region, COUNT(*) as count FROM loci WHERE {where} GROUP BY gene_region"
        if has_gene_region else None)

    sql_params_with_pagination = sql_params + [params["page_size"], (params["page"] - 1) * params["page_size"]]
    return select_query, count_query, motif_count_query, gene_region_count_query, sql_params, sql_params_with_pagination


def row_to_list_dict(row, ot):
    """Convert a duckdb_compat.Row to a dict for the list endpoint.

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
        if affected_raw.lower() in BLANK_STATUS_VALUES:
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


def strchive_inheritance_modes(value):
    """Return a STRchive record's inheritance modes as a flat list of strings, or None.

    The variant-catalog branch of compute_known_disease_info reports the top-level
    ``inheritance`` as a flat list of mode strings, so the STRchive branch has to as well.
    STRchive stores the value as a list already, so wrapping it unconditionally would return a
    list inside a list.

    Args:
        value: A STRchive record's ``inheritance`` value: a list of mode strings, a single mode
            string, or None.

    Returns:
        A list of mode strings, or None when nothing was recorded.
    """
    if not value:
        return None
    if isinstance(value, (list, tuple)):
        modes = [mode for mode in value if mode]
        return modes or None
    return [value]


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
            # Both motif lists, exactly as matches_disease_locus matches them during the build;
            # checking only the reference motifs here would show no disease details for a locus
            # the build flagged as known via a pathogenic motif.
            for strchive_motif in strchive_locus_motifs(interval.data):
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
            # Every record in the shipped STRchive file stores inheritance as a list, but an
            # older cached file may store a bare string, so accept both and always report the
            # flat list of modes the variant-catalog branch below returns.
            "inheritance": strchive_inheritance_modes(disease_info.get("inheritance")),
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


def known_disease_pathogenic_min(chrom, start_0based, end_1based, motif, disease_trees,
                                 strchive_trees):
    """Return the lowest pathogenic minimum recorded for one locus, or None.

    Matches the way locus_annotations.matches_disease_locus decides that a locus is a known
    disease locus: overlap with Jaccard above 0.66 plus a motif match, the variant catalog first
    and the STRchive trees as the fallback. Walking STRchive here too is what keeps the
    pathogenic-threshold filter from dropping a locus that is known only through STRchive, whose
    locus-detail page does report a pathogenic minimum.

    Args:
        chrom: The locus chromosome, with or without a "chr" prefix.
        start_0based: The locus start, 0-based.
        end_1based: The locus end, 1-based.
        motif: The locus motif.
        disease_trees: chrom -> IntervalTree of variant catalog disease loci.
        strchive_trees: chrom -> IntervalTree of STRchive disease loci (may be empty).

    Returns:
        The smallest pathogenic minimum of the matching disease locus, or None when nothing
        matches or the matched locus records no pathogenic minimum.
    """
    if not motif or start_0based is None or end_1based is None:
        return None
    chrom_key = chrom.replace("chr", "") if chrom else ""

    catalog_tree = disease_trees.get(chrom_key) if disease_trees else None
    if catalog_tree:
        for interval in catalog_tree.overlap(start_0based, end_1based):
            if not interval.data:
                continue
            if compute_jaccard(start_0based, end_1based, interval.begin, interval.end) <= 0.66:
                continue
            disease_motifs = [interval.data.get("RepeatUnit", "")] + (interval.data.get("PathogenicMotifs") or [])
            if not any(motifs_match(motif, dm) for dm in disease_motifs):
                continue
            pathogenic_mins = [d.get("PathogenicMin") for d in interval.data.get("Diseases", [])
                               if d.get("PathogenicMin") is not None]
            # The catalog matched, so stop here even when it records no threshold: the STRchive
            # trees are a fallback for loci the catalog does not know, exactly as in
            # matches_disease_locus and compute_known_disease_info.
            return min(pathogenic_mins) if pathogenic_mins else None

    strchive_tree = strchive_trees.get(chrom_key) if strchive_trees else None
    if strchive_tree:
        for interval in strchive_tree.overlap(start_0based, end_1based):
            if not interval.data:
                continue
            if compute_jaccard(start_0based, end_1based, interval.begin, interval.end) <= 0.66:
                continue
            # Both motif lists, the one definition of them, so the filter and the build agree on
            # which loci STRchive knows.
            if not any(motifs_match(motif, sm) for sm in strchive_locus_motifs(interval.data)):
                continue
            return interval.data.get("pathogenic_min")

    return None


def fetch_phenotype_scores(conn, locus_id):
    """Fetch phenotype scores for a locus, if the score tables exist.

    Returns:
        dict with 'per_outlier' and 'per_locus' keys, or {} when the tables are
        absent.
    """
    tables = duckdb_compat.list_tables(conn)
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
    tables = duckdb_compat.list_tables(conn)
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
    """Convert a duckdb_compat.Row to a dict with all columns (NaN -> None)."""
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
        has_mendelian=app.config.get("HAS_MENDELIAN", False),
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
        has_mendelian=app.config.get("HAS_MENDELIAN", False),
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
        has_mendelian=app.config.get("HAS_MENDELIAN", False),
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
        has_mendelian=True,
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


# Flask serves requests on several threads, so two annotation writes can overlap. DuckDB does not
# wait for the other writer the way sqlite3's file lock and busy timeout did: it raises on the
# second connection instead. Every write to the annotations database therefore runs under this one
# process-wide lock, which is enough because that database is only ever written from here.
ANNOTATION_WRITE_LOCK = threading.Lock()


def execute_annotation_write(statement, values):
    """Run one write against the annotations database, serialized against the other writers.

    Args:
        statement: The SQL statement to run.
        values: Its bind parameters.

    Returns:
        None when the write succeeded, or a string describing why it failed.
    """
    with ANNOTATION_WRITE_LOCK:
        try:
            conn = duckdb_compat.connect(app.config["ANNOTATIONS_DB_PATH"])
        except Exception as e:
            return str(e)
        try:
            conn.execute(statement, values)
            conn.commit()
        except Exception as e:
            return str(e)
        finally:
            conn.close()
    return None


@app.route("/api/v1/annotations/<path:locus_id>/note", methods=["PUT"])
def upsert_note(locus_id):
    """Upsert a note for a locus."""
    data = request.get_json(force=True)
    note_text = (data.get("note_text") or "").strip()
    if not note_text:
        return jsonify({"error": "note_text is required and cannot be empty"}), 400

    now = datetime.utcnow().isoformat(sep=" ", timespec="seconds")
    error = execute_annotation_write(
        "INSERT INTO notes (locus_id, note_text, created_at, updated_at) VALUES (?, ?, ?, ?) "
        "ON CONFLICT(locus_id) DO UPDATE SET note_text=excluded.note_text, updated_at=excluded.updated_at",
        (locus_id, note_text, now, now),
    )
    if error:
        return jsonify({"error": "Could not save the note", "detail": error}), 500

    app.config["ANNOTATIONS"]["notes"][locus_id] = {"note_text": note_text, "updated_at": now}
    return jsonify({"locus_id": locus_id, "note": app.config["ANNOTATIONS"]["notes"][locus_id]})


@app.route("/api/v1/annotations/<path:locus_id>/note", methods=["DELETE"])
def delete_note(locus_id):
    """Delete a note for a locus."""
    error = execute_annotation_write("DELETE FROM notes WHERE locus_id = ?", (locus_id,))
    if error:
        return jsonify({"error": "Could not delete the note", "detail": error}), 500
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
    error = execute_annotation_write(
        "INSERT OR IGNORE INTO tags (locus_id, tag, created_at) VALUES (?, ?, ?)",
        (locus_id, tag, now))
    if error:
        return jsonify({"error": "Could not add the tag", "detail": error}), 500

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
    error = execute_annotation_write(
        "DELETE FROM tags WHERE locus_id = ? AND tag = ?", (locus_id, tag))
    if error:
        return jsonify({"error": "Could not remove the tag", "detail": error}), 500

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
        gene_region_rows = (conn.execute(gene_region_count_query, sql_params).fetchall()
                            if gene_region_count_query else [])
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
                       if k not in ("page", "page_size", "outlier_type", "sort_by", "motif_size_clause",
                                    "motif_size_params", "reference_region_parsed",
                                    "search_terms")}

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
                       "motif_size_clause", "motif_size_params", "reference_region_parsed",
                       "search_terms"}
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
            # Bulk-loaded once for the whole export, the way export_json does it: passing the
            # per-locus scores in keeps build_sample_details from running its own
            # fetch_phenotype_scores() query (plus a list_tables() call) for every locus.
            all_phenotype_scores = fetch_all_phenotype_scores(conn)
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
                sample_details = build_sample_details(
                    row_dict, raw_row, conn, lookups,
                    phenotype_scores=all_phenotype_scores.get(row_dict.get("LocusId", ""), {}))
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
    search_raw = request.args.get("search", "")
    max_variation_cluster_size_diff = request.args.get("max_variation_cluster_size_diff", "")
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
        extra_clauses.append("IsInMendelianGene = 1" if (not source_columns or "IsInMendelianGene" in source_columns) else "1=0")

    if known_motifs_only == "1":
        extra_clauses.append("IsKnownMotif = 1" if (not source_columns or "IsKnownMotif" in source_columns) else "1=0")

    if require_above_population in ("long-read", "short-read"):
        # Only reference population-stat columns that actually exist in the DB
        # (a minimal-input build omits them); with none present the filter cannot
        # be evaluated, so match no rows. Mirrors build_api_query.
        datasets = ["HPRC256", "AoU1027"] if require_above_population == "long-read" else ["TenK10K"]
        present = [d for d in datasets if not source_columns or f"{d}_{population_metric}" in source_columns]
        above_terms = [f"({d}_{population_metric} IS NULL OR allele_size > {d}_{population_metric})" for d in present]
        if include_loci_without_population_data:
            if above_terms:
                extra_clauses.append("(" + " AND ".join(above_terms) + ")")
        else:
            if present:
                exists_term = " OR ".join(f"{d}_{population_metric} IS NOT NULL" for d in present)
                extra_clauses.append("((" + exists_term + ") AND " + " AND ".join(above_terms) + ")")
            else:
                extra_clauses.append("1=0")

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
            # Keep unannotated (NULL gene_region) rows: NULL NOT IN (...) is unknown
            # in SQL and would otherwise drop them.
            extra_clauses.append(f"(gene_region IS NULL OR gene_region NOT IN ({','.join('?' * len(db_regions))}))")
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
        # This endpoint parses its own parameters rather than going through validate_params, so it
        # repeats the IUPAC check here; without it compute_canonical_motif raises KeyError -> 500.
        invalid_motifs = [m for m in motifs if not VALID_MOTIF_PATTERN.match(m)]
        if invalid_motifs:
            return jsonify({
                "error": "Invalid parameter",
                "detail": f"motif must contain only IUPAC bases ({''.join(sorted(COMPLEMENT))}); "
                          f"got '{invalid_motifs[0]}'",
            }), 400
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
            gs_clauses = ["GeneTableGeneSymbol ILIKE ?" for _ in symbols]
            extra_clauses.append(f"({' OR '.join(gs_clauses)})")
            extra_params.extend([f"%{gs}%" for gs in symbols])

    if gene_id:
        if not source_columns or "gene_id" in source_columns:
            extra_clauses.append("gene_id = ?")
            extra_params.append(normalize_gene_id(gene_id) or gene_id.strip())
        else:
            extra_clauses.append("1=0")

    # Merged search box, the same OR over locus ids / gene ids / gene symbols / regions that
    # build_api_query applies. swim_plot carries no coordinates, so a region term is answered by
    # a subquery against loci in the same database rather than by materializing its locus ids,
    # which for a whole chromosome would be hundreds of thousands of bind parameters.
    if search_raw:
        search_clauses = []
        search_params = []
        locus_id_values = []
        for term in split_search_terms(search_raw):
            kind, value, error = classify_search_term(term)
            if error:
                return jsonify({"error": "Invalid parameter", "detail": error}), 400
            if kind == "locus_id":
                # Collected into the one set lookup built below instead of a predicate per term.
                locus_id_values.append(value.lower())
            elif kind == "gene_id":
                if not source_columns or "gene_id" in source_columns:
                    search_clauses.append("gene_id = ?")
                    search_params.append(value)
            elif kind == "gene_symbol":
                if not source_columns or "GeneTableGeneSymbol" in source_columns:
                    search_clauses.append("GeneTableGeneSymbol ILIKE ?")
                    search_params.append(f"%{value}%")
            else:
                clause, region_params = reference_region_clause(value)
                search_clauses.append(f"LocusId IN (SELECT LocusId FROM loci AS loci WHERE {clause})")
                search_params.extend(region_params)
        if locus_id_values:
            # Matched the same way the results table matches it (see build_api_query): one
            # case-insensitive set lookup over every pasted id, which needs no wildcard escaping.
            search_clauses.insert(0, f"lower(LocusId) IN ({','.join('?' * len(locus_id_values))})")
            search_params[:0] = locus_id_values
        extra_clauses.append(f"({' OR '.join(search_clauses)})" if search_clauses else "1=0")
        extra_params.extend(search_params)

    # Max variation cluster size diff. swim_plot does not copy VariationClusterSizeDiff, so this
    # is answered against loci too, the same way a region term is.
    # Set below when the variation-cluster filter is active, and consumed once the connection is
    # open to populate the vc_filtered_loci temp table the clause refers to.
    vc_filter_max = None
    if max_variation_cluster_size_diff:
        try:
            max_vc_diff = parse_int64(max_variation_cluster_size_diff)
        except ValueError:
            return jsonify({
                "error": "Invalid parameter",
                "detail": "max_variation_cluster_size_diff must be an integer number of base pairs, "
                          f"got '{max_variation_cluster_size_diff}'",
            }), 400
        # source_columns is the loci column set (configure_app fills DB_COLUMNS_SET from
        # "SELECT * FROM loci LIMIT 1"), and VariationClusterSizeDiff is only ever read from loci,
        # so the same set guards this filter. The clause below deliberately keeps a locus whose
        # cluster size is unknown, so a database with no such column at all is that same situation
        # for every locus: add no clause, exactly as build_api_query does. Matching nothing here
        # would empty the swim plot for a request whose results table still returns every locus.
        if not source_columns or "VariationClusterSizeDiff" in source_columns:
            # This endpoint runs one query per motif category, so an inline subquery over loci
            # would re-scan that column 25 times for a single chart. Materialize the matching
            # LocusIds into a temp table once, after the connection is opened, and join against
            # that instead.
            vc_filter_max = max_vc_diff
            extra_clauses.append("LocusId IN (SELECT LocusId FROM vc_filtered_loci)")

    if min_pli:
        try:
            pli_value = float(min_pli)
        except ValueError:
            return jsonify({"error": "Invalid parameter", "detail": f"min_pli must be a number, got '{min_pli}'"}), 400
        if not source_columns or "pLI" in source_columns:
            extra_params.append(pli_value)
            extra_clauses.append("pLI >= ?")
        else:
            extra_clauses.append("1=0")

    if phenotype_keyword_raw:
        keywords = [kw.strip() for kw in phenotype_keyword_raw.split(",") if kw.strip()]
        if keywords:
            kw_clauses = ["phenotype_description ILIKE ?" for _ in keywords]
            extra_clauses.append(f"({' OR '.join(kw_clauses)})")
            extra_params.extend([f"%{kw}%" for kw in keywords])

    if sample_id_keyword_raw:
        keywords = [kw.strip() for kw in sample_id_keyword_raw.split(",") if kw.strip()]
        if keywords:
            kw_clauses = ["lower(sample_id) = lower(?)" for _ in keywords]
            extra_clauses.append(f"({' OR '.join(kw_clauses)})")
            extra_params.extend(keywords)

    if sample_id_like_raw:
        keywords = [kw.strip() for kw in sample_id_like_raw.split(",") if kw.strip()]
        if keywords:
            # Escaped like the loci list's sample_id_like, so the same keyword selects the same
            # samples in both views.
            kw_clauses = ["sample_id ILIKE ? ESCAPE '\\'" for _ in keywords]
            extra_clauses.append(f"({' OR '.join(kw_clauses)})")
            extra_params.extend([f"%{escape_like_wildcards(kw)}%" for kw in keywords])

    if min_expansion:
        try:
            min_expansion_value = parse_int64(min_expansion)
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
            else:
                extra_clauses.append("1=0")
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
        if "swim_plot" not in duckdb_compat.list_tables(conn):
            return jsonify({"error": "swim_plot table not found", "detail": "The database has no swim_plot table."}), 500

        if vc_filter_max is not None:
            # A read-only main database still allows temp tables (they live in memory, in this
            # connection's own temp schema), so this costs one scan of loci's one relevant
            # column per request instead of one per category.
            conn.execute(
                "CREATE OR REPLACE TEMP TABLE vc_filtered_loci AS SELECT LocusId FROM loci "
                "WHERE VariationClusterSizeDiff IS NULL "
                "OR COALESCE(TRY_CAST(VariationClusterSizeDiff AS BIGINT), 0) <= ?",
                (vc_filter_max,))

        # LocusId and sample_id follow allele_size in the ORDER BY so the top 500 is a decided
        # set rather than an arbitrary one: hundreds of rows share an allele size, and DuckDB
        # reads the table with several threads, so ordering by allele_size alone would pick a
        # different 500 of the tied rows on each request.
        for motif_category in categories:
            query = ("""SELECT rowid, allele_size, motif_category, affected_status,
                       LocusId, sample_id, Motif, gene_region, GeneTableGeneSymbol,
                       purity, methylation, family_id, sex, analysis_status, phenotype_description
                FROM swim_plot
                WHERE outlier_type = ? AND motif_category = ?""" + extra_where + """
                ORDER BY allele_size DESC, LocusId, sample_id LIMIT 500""")
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
        if "swim_plot" not in duckdb_compat.list_tables(conn):
            return jsonify({"error": "swim_plot table not found", "detail": "The database has no swim_plot table."}), 500

        locus_row = None
        if locus_id:
            # One row per outlier sample matches, and they all carry the same MotifSize and
            # CanonicalMotif, but the LIMIT still needs an ORDER BY to name a single one of them.
            locus_row = conn.execute(
                "SELECT MotifSize, CanonicalMotif FROM swim_plot WHERE LocusId = ? AND outlier_type = ? "
                "ORDER BY sample_id, allele_size LIMIT 1",
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
            "min_sigma_percentile": {"type": "float", "min": 0.0, "max": 1.0, "description": "Keep loci whose HPRC256/AoU1027 stdev percentile is in the top (1 - value) fraction"},
            "chrom": {"type": "string", "description": "Filter to a single chromosome (e.g. chr1)"},
            "reference_region": {"type": "string", "description": "Keep loci overlapping a reference-genome region, as 'chrom:start-end' with a 0-based start and an exclusive end (e.g. chr16:11579459-11579529). The 'chr' prefix and comma separators are optional; a bare position or a bare chromosome is also accepted"},
            "motif": {"type": "string"},
            "phenotype_keyword": {"type": "string"},
            "sample_id_keyword": {"type": "string", "description": "Comma-separated exact sample IDs from dropdown"},
            "sample_id_like": {"type": "string", "description": "Comma-separated keywords for partial sample ID matching, OR-combined"},
            "apply_pathogenic_threshold": {"type": "boolean"},
            "known_loci_only": {"type": "boolean"},
            "known_motifs_only": {"type": "boolean"},
            "mendelian_genes_only": {"type": "boolean"},
            "search": {"type": "string", "description": "Comma-separated search terms, OR-combined. Each term is routed by its shape: a locus id (2-89831737-89831752-CCATT), an Ensembl gene id (ENSG00000102081), a reference region (chr16:11579459-11579529, chr16:11579459 or chr16), or otherwise a gene symbol matched as a substring. Commas inside a region's coordinates are treated as thousands separators, not term separators"},
            "max_variation_cluster_size_diff": {"type": "integer", "description": "Keep loci whose TRExplorer variation cluster is at most this many base pairs larger than the repeat itself. Loci with no variation cluster annotation are kept"},
            "exclude_vc_depth_filtered": {"type": "boolean", "description": "Drop loci where no variation cluster could be computed because the region had too little coverage in the population dataset the clusters were derived from (VariationClusterFilterReason = DEPTH), shown as a gold crossed-out circle in the VC column"},
            "exclude_vc_size_filtered": {"type": "boolean", "description": "Drop loci where a variation cluster was computed and then discarded because it was wider than the maximum allowed (VariationClusterFilterReason = EXTENSION), shown as a dark red crossed-out circle in the VC column"},
            # The legacy parameters below are superseded by 'search', which the UI sends instead; they still work for direct API callers.
            "gene_id": {"type": "string"},
            "gene_symbol": {"type": "string"},
            "locus_id": {"type": "string"},
            "tag": {"type": "string", "description": "Filter to loci with this tag"},
            "motif_size": {"type": "string", "description": "Motif size filter. Comma-separated entries: exact (3), range (3-6), open-end (5-), open-start (-10)"},
        },
        "sort_options": list(VALID_SORTS),
        "total_loci": app.config.get("TOTAL_LOCI", 0),
        # None for a database built before the pipeline started recording the sample count
        "total_samples": app.config.get("TOTAL_SAMPLES"),
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


def configure_app(db_path, sample_table=None, known_loci_json=None,
                  annotations_db="annotations.duckdb", strchive_loci_json=None):
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
    # cached STRchive-loci JSON, when supplied, adds the detail-page fallback. Either catalog is
    # enough on its own: load_known_disease_loci takes filepath=None to load only the STRchive
    # fallback, so a STRchive-only invocation still gets its disease details.
    known_loci_path = known_loci_json if known_loci_json and os.path.exists(known_loci_json) else None
    strchive_path = strchive_loci_json if strchive_loci_json and os.path.exists(strchive_loci_json) else None
    if known_loci_path or strchive_path:
        disease_trees, strchive_trees, locus_lookup = load_known_disease_loci(
            known_loci_path, fetch_strchive=False,
            strchive_filepath=strchive_path,
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
    conn = duckdb_compat.connect(db_path, read_only=True)
    try:
        total_loci = conn.execute("SELECT COUNT(*) FROM loci").fetchone()[0]
        # None for a database built before the pipeline started recording the sample count.
        total_samples = read_sample_count(conn)
        # LIMIT 0 rather than duckdb_compat.table_columns: the columns are advertised as an
        # ordered list via /api/v1/schema, and a query's description keeps the table's order
        # where the catalog lookup returns an unordered set. No row is read either way.
        db_columns = [desc[0] for desc in conn.execute("SELECT * FROM loci LIMIT 0").description]
        tables = duckdb_compat.list_tables(conn)
        # The loci-side phenotype_keyword filter is answered from swim_plot so it selects the same
        # loci the swim plot shows (see build_api_query), which needs both the table and that
        # column to be present.
        has_swim_plot_phenotypes = (
            "swim_plot" in tables
            and "phenotype_description" in duckdb_compat.table_columns(conn, "swim_plot"))
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
    app.config["TOTAL_SAMPLES"] = total_samples
    app.config["DB_COLUMNS"] = db_columns
    app.config["DB_COLUMNS_SET"] = db_columns_set
    app.config["HAS_MENDELIAN"] = "mendelian_violations" in tables
    app.config["HAS_SWIM_PLOT_PHENOTYPES"] = has_swim_plot_phenotypes
    app.config["SAMPLE_IDS"] = sample_ids

    # Known-disease-locus filter set + pathogenic-threshold map (from the loci table's
    # KnownDiseaseLocus column, and the disease catalog for thresholds when available).
    known_disease_locus_ids = set()
    known_disease_locus_thresholds = {}
    if "KnownDiseaseLocus" in db_columns_set:
        conn = duckdb_compat.connect(db_path, read_only=True)
        try:
            known_disease_locus_ids = {
                row[0] for row in conn.execute(
                    "SELECT LocusId FROM loci WHERE KnownDiseaseLocus IS NOT NULL AND KnownDiseaseLocus != ''"
                )
            }
            if (disease_trees or strchive_trees) and known_disease_locus_ids:
                # Only known-disease loci can contribute a threshold (known_disease_pathogenic_min
                # mirrors how KnownDiseaseLocus was set at build time, STRchive fallback
                # included), so restrict the scan to that id set instead of every locus.
                placeholders = ",".join("?" * len(known_disease_locus_ids))
                for row in conn.execute(
                    f"SELECT LocusId, Chrom, Start0Based, End1Based, Motif FROM loci "
                    f"WHERE LocusId IN ({placeholders})",
                    tuple(known_disease_locus_ids)):
                    pathogenic_min = known_disease_pathogenic_min(
                        row[1], row[2], row[3], row[4], disease_trees, strchive_trees)
                    if pathogenic_min is not None:
                        known_disease_locus_thresholds[row[0]] = pathogenic_min
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

    total_samples = app.config["TOTAL_SAMPLES"]
    if total_samples is None:
        print(f"Loaded {total_loci:,d} loci (sample count not recorded; rebuild the database to add it)")
    else:
        print(f"Loaded {total_loci:,d} loci from {total_samples:,d} genotyped samples")
    print(f"Loaded {len(app.config['LOOKUPS']['sample']):,d} samples")
    print(f"Loaded {len(app.config['ANNOTATIONS']['notes'])} notes, {len(app.config['ANNOTATIONS']['all_tags'])} unique tags")
    for col in ("CanonicalMotif", "IsKnownMotif", "IsInMendelianGene"):
        if col not in app.config["DB_COLUMNS_SET"]:
            print(f"Note: database missing optional column: {col} (using fallback)")
    print(f"Mendelian QC page: {'available' if app.config['HAS_MENDELIAN'] else 'hidden (no mendelian_violations table)'}")

    app.config["BIND_HOST"] = args.host
    app.config["ALLOW_REMOTE_WRITES"] = args.allow_remote_writes
    if args.host not in LOOPBACK_HOSTS:
        if args.allow_remote_writes:
            print(f"\nWARNING: bound to {args.host} with --allow-remote-writes: anyone who can "
                  f"reach this port can read the results and edit notes and tags.")
        else:
            print(f"\nNOTE: bound to {args.host}, so other machines can read the results. Notes "
                  f"and tags can only be edited from this machine; pass --allow-remote-writes to "
                  f"lift that (there is no authentication).")

    print(f"\nStarting server on {args.host}:{args.port}")
    print(f"Web UI: http://{args.host}:{args.port}/")
    app.run(host=args.host, port=args.port, debug=args.debug, use_reloader=args.debug, threaded=True)


if __name__ == "__main__":
    main()
