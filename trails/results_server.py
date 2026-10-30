#!/usr/bin/env python3
"""Flask-based read-only API server for TRails tandem-repeat outlier results.

Serves data from a single pre-computed DuckDB database of locus-level outlier
statistics, enriched with sample metadata, gene-disease associations, and known
disease loci.

The server tolerates optional tables being absent. The ``loci`` table is the only
hard requirement (see ``REQUIRED_COLUMNS``); ``swim_plot``, the phenotype-score
tables and the Mendelian-violation tables are all optional, and the corresponding
pages/endpoints degrade gracefully when they are missing.
"""

import argparse
import collections
from datetime import datetime
import ipaddress
import json
import math
import os
import re
import socket
import sys
import threading
import traceback
from urllib.parse import urlsplit
import zlib

import flask
from flask import Flask, request, Response
from werkzeug.exceptions import HTTPException
import jinja2
import numpy as np

# The only module that imports duckdb; everything here goes through its sqlite3-shaped API.
import duckdb_compat

# get_population_p99_threshold, the build's population bar, which the exports apply too, and
# BLANK_STATUS_VALUES, the spellings of a status cell that mean nothing was recorded
import analysis_columns

# Canonical motifs, and the bases they are defined for
from motif_utilities import COMPLEMENT, compute_canonical_motif

# The number of samples the database was built from and the cohorts that gated its #Aff / #Fam
# counts, both recorded in its metadata table, and the table of the (locus, sample) pairs #Aff
# counts, which the grouped view reads.
from result_database import AFFECTED_UNSOLVED_OUTLIERS_TABLE, read_population_gate_cohorts, read_sample_count

# The bases compute_canonical_motif can canonicalize. Deriving this from motif_utilities.COMPLEMENT
# rather than spelling the alphabet out here keeps the two from drifting: anything outside the table
# makes reverse_complement raise KeyError, so the "motif" query parameter is validated against it.
VALID_MOTIF_PATTERN = re.compile(f"^[{''.join(sorted(COMPLEMENT))}]+$", re.IGNORECASE)

# Known-disease-locus loading and matching live in the locus_annotations module, which needs
# intervaltree. Import them defensively so the server stays importable without it; the
# locus-detail page's known-disease annotations simply degrade to "no match" in that case.
try:
    from locus_annotations import find_matching_disease_locus, load_known_disease_loci
except ImportError as import_error:
    # A missing module is the supported degraded mode; a missing *name* means the two files have
    # drifted, which would otherwise disable the known-disease annotations without a word.
    print(f"WARNING: falling back to the built-in known-disease stubs: {import_error}")

    def load_known_disease_loci(filepath, strchive_filepath=None, build_locus_lookup=False):
        """Disabled fallback: no disease catalog available -> empty trees/lookup.

        Mirrors the real loader's (interval_trees, strchive_trees, locus_lookup)
        return shape so the startup call site unpacks identically.
        """
        return {}, {}, {}

    def find_matching_disease_locus(chrom, start, end, motif, interval_trees, strchive_trees=None):
        """Disabled fallback: with no disease catalog, no locus matches -> (None, None)."""
        return None, None

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTLIER_TYPE_MAP = {"all": "AllAlleles", "short": "ShortAlleles", "hemi": "HemizygousAlleles"}

# Population cohorts whose {cohort}_99thPercentile and {cohort}_MaxAllele columns the population
# filters and the page's Δpop column compare an affected allele against
LONG_READ_POPULATION_DATASETS = ("HPRC256", "AoU1027")
POPULATION_DATASETS = LONG_READ_POPULATION_DATASETS + ("TenK10K",)

# The cohorts each require_above_population value compares against ("all" is for short-read callsets, which are
# compared against both kinds of cohort)
POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION = {
    "long-read": LONG_READ_POPULATION_DATASETS,
    "short-read": ("TenK10K",),
    "all": POPULATION_DATASETS,
}

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

# The values of require_above_unaffected and require_families_above_unaffected, and the prefix of
# the affected-allele column each one compares
AFFECTED_ALLELE_COLUMN_PREFIXES = {"first": "First", "second": "Second", "third": "Third"}

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
    # TSV only when it supplies them. A database without them reports them as None
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
    "NumAffectedUnsolvedFamiliesAboveUnaffected",
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
    """Return the set of LocusIds matching a user-defined tag."""
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
    still valid JSON and parses identically in the browser. json.dumps already writes
    every non-ASCII character as a backslash-u escape, including the line separators
    U+2028 and U+2029, which a <script> block must not contain raw.
    """
    return (
        json.dumps(value)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
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
        "--initial-require-above-population",
        choices=[""] + list(POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION),
        default="",
        help="Population filter that the page starts with when it is opened with no saved filters, e.g. 'all' "
             "for a short-read callset (default: none).",
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


def load_sample_table(filepath):
    """Load sample metadata and create lookup dictionaries.

    Reads the table with input_tables.read_sample_metadata, the build's reader, so it accepts
    and normalizes exactly what the build does (trails.py forwards the same file to both): only
    ``sample_id`` is required, column names are matched case/underscore-insensitively, and an
    absent ``affected_status`` / ``analysis_status`` reads as unknown. The build keeps only the
    columns it uses, so every other column (e.g. ``ancestry``) is added to each sample's row
    here, verbatim, for the pages and exports to show.

    Args:
        filepath: Path to the sample-metadata TSV (optionally gzipped).

    Returns:
        Tuple (sample_id_to_row, affected_status_lookup, analysis_status_lookup).
    """
    # Imported here rather than at the top: input_tables imports locus_annotations, which needs
    # intervaltree, and the server runs without it as long as no sample table is given.
    import input_tables
    import pandas as pd

    sample_id_to_row, affected_status_lookup, analysis_status_lookup = input_tables.read_sample_metadata(filepath)

    raw_df = pd.read_table(filepath, dtype=str, keep_default_na=False)
    raw_sample_id_column = input_tables.match_columns(raw_df, ["sample_id"])["sample_id"]
    read_columns = {input_tables.normalize_column_name(column)
                    for sample_row in sample_id_to_row.values() for column in sample_row}
    other_columns = [column for column in raw_df.columns
                     if input_tables.normalize_column_name(column) not in read_columns]
    if other_columns:
        # The first row of a duplicated sample_id is the one read_sample_metadata keeps
        for row in raw_df.drop_duplicates(subset=[raw_sample_id_column], keep="first").to_dict(orient="records"):
            sample_id_to_row[row[raw_sample_id_column]].update({column: row[column] for column in other_columns})
    return sample_id_to_row, affected_status_lookup, analysis_status_lookup


def compute_violation_fraction(violations, total):
    """Return violations / total, with a missing (NULL) count read as 0, and 0 when there is no total.

    Shared by the startup Mendelian warnings and the QC2 page, so the two read the counts alike.
    """
    return (violations or 0) / total if total else 0


def load_mendelian_warnings(db_path, threshold=0.10):
    """Load Mendelian violation data and flag samples with high violation rates.

    The Mendelian-violation tables live in the results database. If the
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
            autosome_rate = compute_violation_fraction(row["autosome_violations"], row["autosome_total"])
            chrX_rate = compute_violation_fraction(row["chrX_violations"], row["chrX_total"])
            chrY_rate = compute_violation_fraction(row["chrY_violations"], row["chrY_total"])
            if autosome_rate > threshold or chrX_rate > threshold:
                motif_rates = {}
                for category in ["1bp", "2bp", "3bp", "4bp", "5bp", "6bp", "7_24bp", "25plusbp"]:
                    violations = row[f"motif_{category}_violations"] or 0
                    total = row[f"motif_{category}_total"] or 0
                    if total > 0:
                        motif_rates[category] = {
                            "rate": compute_violation_fraction(violations, total),
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
        # A sample/bin with no rank-1 locus is left out of rank1 and keeps its full count in top10.
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
# Request checks
# ---------------------------------------------------------------------------

# No response carries CORS headers, on purpose. The pages are served by this server, so their
# requests are same-origin, and scripts and notebooks are not subject to CORS at all. The only
# readers the headers would let in are pages on other sites the user has open, and those must not
# read sample ids, phenotypes and the other results.

# Loopback addresses. A server bound to one of these is reachable only from this machine, which
# is what makes the unauthenticated annotation routes acceptable by default.
LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}


def find_this_machines_host_names():
    """Return the names other machines may use to reach this one, lower-cased.

    These are socket.gethostname() and socket.getfqdn(), each also with the ".local" suffix
    under which macOS and many Linux machines announce themselves on the local network.

    Returns:
        A set of host names.
    """
    host_names = set()
    for name in (socket.gethostname(), socket.getfqdn()):
        name = name.lower().rstrip(".")
        if name:
            host_names.update((name, name if name.endswith(".local") else f"{name}.local"))
    return host_names


def parse_host_name(host):
    """Return the host name of a Host header value, lower-cased and without its port.

    Args:
        host: The Host header value, e.g. "localhost:5050" or "[::1]:5050".

    Returns:
        The host name, e.g. "localhost" or "::1", or "" when the value will not parse.
    """
    try:
        return (urlsplit(f"//{host}").hostname or "").rstrip(".")
    except ValueError:
        return ""


def is_ip_address(text):
    """Return whether text is an IPv4 or IPv6 address rather than a name."""
    try:
        ipaddress.ip_address(text)
    except ValueError:
        return False
    return True


@app.before_request
def reject_requests_for_another_host():
    """Refuse a request whose Host header names a host other than this machine.

    A web page on another site can point its own domain name at this machine's address (DNS
    rebinding) and then reach this server as if it were the page's own site, reading the results
    and editing notes and tags: the connection comes from this machine, and the Origin header
    matches the Host. The browser still names the page's domain in the Host header, so whatever
    address the server is bound to, it answers only requests that name it by an IP address (a
    rebinding page needs a domain name), by a loopback name, or by one of this machine's own
    names, which main() records in THIS_MACHINES_HOST_NAMES. Writes from other machines are
    refused separately, by restrict_writes.
    """
    host_name = parse_host_name(request.host)
    if (host_name in LOOPBACK_HOSTS or host_name in app.config.get("THIS_MACHINES_HOST_NAMES", ())
            or is_ip_address(host_name)):
        return None
    accepted_names = sorted(LOOPBACK_HOSTS | set(app.config.get("THIS_MACHINES_HOST_NAMES", ())))
    return jsonify({
        "error": "unexpected host name",
        "detail": f"This server only answers requests addressed to an IP address of this machine or to "
                  f"one of {', '.join(accepted_names)}, not '{request.host}'.",
    }), 403


@app.before_request
def restrict_writes():
    """Reject state-changing requests that this server should not accept.

    Two separate exposures, both of which matter because the annotation routes (notes and tags)
    have no authentication:

      * Another origin. A POST with a plain content type is a "simple" request the browser sends
        without asking this server first, so another page the user has open could write to this
        database if the Origin header were not checked.
      * Another machine. Bound to a non-loopback address the write routes are reachable by anyone
        who can route to the port, so they are refused unless --allow-remote-writes says otherwise.

    The Origin check relies on reject_requests_for_another_host: a page that rebinds its own
    domain to this machine sends an Origin matching its Host, and only the Host check refuses it.
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


def is_true_flag(name):
    """Return whether the request sets a boolean query parameter, as true, 1 or yes in any case.

    Args:
        name: The query parameter's name.

    Returns:
        True when the parameter is set to one of those values.
    """
    return request.args.get(name, "").lower() in ("true", "1", "yes")


def validate_params():
    """Validate and parse query parameters from the request.

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

    for name in ("require_above_unaffected", "require_families_above_unaffected"):
        value = request.args.get(name)
        if value:
            if value not in AFFECTED_ALLELE_COLUMN_PREFIXES:
                errors.append(f"{name} must be one of: {', '.join(AFFECTED_ALLELE_COLUMN_PREFIXES)}; got '{value}'")
            else:
                params[name] = value

    rap = request.args.get("require_above_population")
    if rap:
        if rap not in POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION:
            errors.append(f"require_above_population must be one of: "
                          f"{', '.join(POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION)}; got '{rap}'")
        else:
            params["require_above_population"] = rap

    if is_true_flag("include_loci_without_population_data"):
        params["include_loci_without_population_data"] = True

    # population_metric: which population column require_above_population and
    # min_repeats_above_population compare against, "99th" (default) or "max". The grouped view's
    # Δpop values and its exports use it too, whether or not a population filter is set.
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

    # min_repeats_above_population: the smallest Δpop the page's results table may show, i.e. the
    # affected allele minus the largest population value (see build_api_where_clause)
    if request.args.get("min_repeats_above_population"):
        try:
            params["min_repeats_above_population"] = parse_int64(request.args["min_repeats_above_population"])
        except ValueError:
            errors.append(f"min_repeats_above_population must be an integer, got '{request.args['min_repeats_above_population']}'")

    for name in ("gene_regions", "exclude_gene_regions"):
        if request.args.get(name):
            regions = [r.strip() for r in request.args[name].split(",")]
            invalid = [r for r in regions if r not in GENE_REGION_MAP]
            if invalid:
                errors.append(f"{name} contains invalid values: {invalid}. Valid: {sorted(GENE_REGION_MAP)}")
            else:
                params[name] = regions

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
        if is_true_flag(flag_name):
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
    # separate locus_id / reference_region / gene_symbol / gene_id parameters below are for
    # anyone calling the API directly; the UI sends this one instead.
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

    # Keep the loci where #Aff counts this sample, or one of this family's samples. A row of the
    # grouped view (/api/v1/outlier_samples) expands to these loci.
    for name in ("affected_unsolved_outlier_sample_id", "affected_unsolved_outlier_family_id"):
        if request.args.get(name, "").strip():
            params[name] = request.args[name].strip()

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


# validate_params entries that are not filters, or are the SQL form of a filter whose
# user-facing value is kept under its own name, so filters_applied leaves them out
NON_FILTER_PARAMS = ("page", "page_size", "outlier_type", "sort_by", "motif_size_clause", "motif_size_params",
                     "reference_region_parsed", "search_terms")


def build_filters_applied(params):
    """Return the filters of validated query parameters, as the request gave them.

    Args:
        params: dict of validated query parameters from validate_params().

    Returns:
        The params without NON_FILTER_PARAMS, for a response's filters_applied.
    """
    return {k: v for k, v in params.items() if k not in NON_FILTER_PARAMS}


def has_affected_unsolved_outliers_table():
    """Whether the database has the table the grouped view needs.

    A database without AFFECTED_UNSOLVED_OUTLIERS_TABLE gets it when it is rebuilt.

    Returns:
        True when the table exists.
    """
    return AFFECTED_UNSOLVED_OUTLIERS_TABLE in app.config.get("DB_TABLES", set())


def population_datasets_compared(require_above_population):
    """The population cohorts that Δpop and the population filters compare an allele against.

    Mirrors renderLocusRow in index_page_template.html: the cohorts of the
    require_above_population gate when it is set, and every cohort otherwise.

    Args:
        require_above_population: "long-read", "short-read", "all", or None when the gate is off.

    Returns:
        Tuple of cohort prefixes.
    """
    return POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION.get(require_above_population, POPULATION_DATASETS)


def compute_delta_population(allele_size, locus_row, params):
    """An allele minus the largest population value, as the page's Δpop column computes it.

    Uses the cohorts of population_datasets_compared and the population_metric parameter.

    Args:
        allele_size: The allele size to compare, in repeat units.
        locus_row: dict of the locus's loci-table columns.
        params: dict of validated query parameters.

    Returns:
        The difference, or None when the allele or every population value is missing.
    """
    metric = "MaxAllele" if params.get("population_metric") == "max" else "99thPercentile"
    values = []
    for dataset in population_datasets_compared(params.get("require_above_population")):
        value = locus_row.get(f"{dataset}_{metric}")
        if value is not None and not (isinstance(value, float) and math.isnan(value)):
            values.append(value)
    if allele_size is None or not values:
        return None
    return allele_size - max(values)


def build_api_order_by(params, source_columns):
    """Build the ORDER BY clause from API parameters.

    Args:
        params: dict of validated query parameters.
        source_columns: this database's loci column set. Most sort keys read an optional
            annotation column (the sigma percentiles, the variation-cluster size, the
            phenotype scores), and a build whose input TSV did not supply one has no such
            column at all, so a sort key naming it is skipped rather than raising "no such
            column" for the whole search.

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
    # rather than CAST because CAST raises on a value that will not parse; with TRY_CAST such a
    # value sorts as NULL, which DESC puts last alongside the loci that have no variation
    # cluster at all. The SORT_MAPPING column name stays the real one so the missing-column
    # check below still tests a column this database either has or does not.
    SORT_EXPRESSIONS = {"variation_cluster": "TRY_CAST(VariationClusterSizeDiff AS BIGINT)"}

    sequence = list(sort_list)
    for tiebreak in ("count", "region", "size"):
        if tiebreak not in sequence:
            sequence.append(tiebreak)

    order_parts = []
    for field in sequence:
        col, direction = SORT_MAPPING[field]
        if col not in source_columns:
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


def build_pathogenic_threshold_clause(table_name, allele_size_column):
    """Build the SQL that keeps the alleles at or above their known disease locus's pathogenic minimum.

    Shared by the apply_pathogenic_threshold filter of /api/v1/loci, which compares a locus's
    largest affected allele, and of /api/v1/swim_plot_data, which compares each allele it draws.

    Args:
        table_name: The table holding LocusId and the allele size column: "loci" or "swim_plot".
        allele_size_column: The allele size column compared.

    Returns:
        Tuple (clause, params) to append to a WHERE clause. The clause matches nothing when no
        known disease locus has a pathogenic minimum.
    """
    thresholds = app.config.get("KNOWN_DISEASE_LOCUS_THRESHOLDS", {})
    if not thresholds:
        return "1=0", []
    threshold_clauses = []
    threshold_params = []
    for locus_id, threshold in sorted(thresholds.items()):
        threshold_clauses.append(f"({table_name}.LocusId = ? AND {table_name}.{allele_size_column} >= ?)")
        threshold_params.extend([locus_id, threshold])
    return f"({' OR '.join(threshold_clauses)})", threshold_params


def build_require_above_population_clause(params, table_name, allele_size_column):
    """Build the SQL of the require_above_population filter for one allele size column.

    Keeps an allele above the population value (99th percentile, or the maximum when
    population_metric is "max") of every cohort the filter names that has data there. A locus
    with no data in any of those cohorts is dropped unless include_loci_without_population_data
    is set. A cohort whose column the loci table lacks is left out, and with none left nothing
    matches. Shared like build_pathogenic_threshold_clause.

    Args:
        params: dict of validated query parameters holding require_above_population.
        table_name: The table holding the allele size column and the {cohort}_{metric} columns:
            "loci", or "swim_plot", which copies them from loci.
        allele_size_column: The allele size column compared.

    Returns:
        The clause to append to a WHERE clause, or None when it would check nothing.
    """
    source_columns = app.config["DB_COLUMNS_SET"]
    m = "MaxAllele" if params.get("population_metric") == "max" else "99thPercentile"
    pop_datasets = POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION[params["require_above_population"]]
    present = [d for d in pop_datasets if f"{d}_{m}" in source_columns]
    above_terms = [f"({table_name}.{d}_{m} IS NULL OR {table_name}.{allele_size_column} > {table_name}.{d}_{m})"
                   for d in present]
    if params.get("include_loci_without_population_data"):
        return "(" + " AND ".join(above_terms) + ")" if above_terms else None
    if not present:
        return "1=0"
    exists_term = " OR ".join(f"{table_name}.{d}_{m} IS NOT NULL" for d in present)
    return "((" + exists_term + ") AND " + " AND ".join(above_terms) + ")"


def build_above_unaffected_clause(affected_allele_column, ot, source_columns):
    """Build the SQL that keeps the loci whose affected allele is above every unaffected one.

    Shared by require_above_unaffected and require_families_above_unaffected. The allele is
    compared with the largest unaffected allele directly: min_expansion is a margin above the
    reference allele only. A NULL unaffected allele means the locus has no unaffected sample, so
    any affected allele passes. 0 does not: it is a real allele size (a full deletion of the
    repeat), and treating it as missing would pass every affected allele at that locus without
    comparing it to anything.

    Args:
        affected_allele_column: The loci column holding the affected allele compared, e.g.
            "SecondAffectedAlleleSize_AllAlleles".
        ot: The outlier type's column suffix, e.g. "AllAlleles".
        source_columns: This database's loci column set.

    Returns:
        The clause to append to a WHERE clause. It matches no locus when the database lacks either
        column, since the filter cannot be evaluated, rather than returning every locus.
    """
    unaffected_allele_column = f"FirstUnaffectedAlleleSize_{ot}"
    if not {affected_allele_column, unaffected_allele_column} <= source_columns:
        return "1=0"
    return (f"(loci.{affected_allele_column} > loci.{unaffected_allele_column} "
            f"OR loci.{unaffected_allele_column} IS NULL AND loci.{affected_allele_column} IS NOT NULL)")


def build_api_where_clause(params):
    """Build the WHERE clause that picks the loci matching validated API query parameters.

    Shared by build_api_query (the locus list, its count, the charts and the export) and by
    /api/v1/outlier_samples, which groups the same loci by sample or family, so the two views
    always filter identically.

    Args:
        params: dict of validated query parameters from validate_params().

    Returns:
        Tuple (where, sql_params): a condition on the `loci` table and its parameters.
    """
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    source_columns = app.config["DB_COLUMNS_SET"]

    # The loci with an outlier of the requested type. Guarded like every filter: a database built
    # without this outlier type simply has no rows for it, which is not the same thing as a SQL
    # error. OUTLIER_TYPE_MAP advertises all three types regardless of what the loaded database
    # actually contains.
    if f"OutlierSampleIds_{ot}" in source_columns:
        clauses = [f"(loci.OutlierSampleIds_{ot} IS NOT NULL AND loci.OutlierSampleIds_{ot} != '')"]
    else:
        clauses = ["1=0"]

    if params["outlier_type"] == "hemi":
        clauses.append("Chrom IN ('chrX', 'chrY')")
        # The build writes "" (not NULL) for a locus with no hemizygous alleles
        # (allele_histograms.convert_counts_to_histogram_string returns "" for an empty count
        # dict), so IS NOT NULL alone matches every row and filters nothing. Guarded on the column
        # like the filters, so an autosome-only database matches nothing instead of raising.
        if "HemizygousAlleleHistogram" in source_columns:
            clauses.append("HemizygousAlleleHistogram IS NOT NULL AND HemizygousAlleleHistogram != ''")
        else:
            clauses.append("1=0")

    filter_clauses, sql_params = build_locus_filter_clauses(params)
    return " AND ".join(clauses + filter_clauses), sql_params


def build_locus_filter_clauses(params):
    """Build the SQL of every locus filter set in validated API query parameters.

    Only the filters: build_api_where_clause adds the restriction to the loci with an outlier of
    the requested type, which /api/v1/swim_plot_data does not need, since its rows each carry
    their own outlier type.

    Args:
        params: dict of validated query parameters from validate_params().

    Returns:
        Tuple (clauses, sql_params): conditions on the `loci` table, to be ANDed together, and
        their parameters in order.
    """
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    source_columns = app.config["DB_COLUMNS_SET"]
    clauses = []
    sql_params = []

    if "min_expansion" in params:
        clauses.append(f"loci.FirstAffectedAlleleSize_{ot} >= loci.NumRepeatsInReference + ?")
        sql_params.append(params["min_expansion"])

    if "require_above_unaffected" in params:
        clauses.append(build_above_unaffected_clause(
            f"{AFFECTED_ALLELE_COLUMN_PREFIXES[params['require_above_unaffected']]}AffectedAlleleSize_{ot}",
            ot, source_columns))

    if "require_families_above_unaffected" in params:
        clauses.append(build_above_unaffected_clause(
            f"{AFFECTED_ALLELE_COLUMN_PREFIXES[params['require_families_above_unaffected']]}AffectedAlleleSize_{ot}_ByFamily",
            ot, source_columns))

    if "require_above_population" in params:
        clause = build_require_above_population_clause(params, "loci", f"FirstAffectedAlleleSize_{ot}")
        if clause:
            clauses.append(clause)

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
        if f"FirstAffectedAlleleSize_{ot}" in source_columns:
            clauses.append(f"loci.FirstAffectedAlleleSize_{ot} >= ?")
            sql_params.append(params["min_repeats_threshold"])
        else:
            clauses.append("1=0")

    # Min repeats above population, i.e. a minimum on the page's Δpop column (see renderLocusRow in
    # index_page_template.html): the affected allele minus the largest value among the population
    # cohorts. Like Δpop, it uses the cohorts and population metric of the require_above_population
    # gate, or every cohort when that gate is off, and a locus with no population data has no Δpop,
    # so it never passes. Every cohort with data must be at least N - 0.5 below the allele, not N,
    # because the page rounds Δpop to a whole number and some TenK10K percentiles are fractional: a
    # difference rounds to N or more exactly when it is at least N - 0.5.
    if "min_repeats_above_population" in params:
        m = "MaxAllele" if params.get("population_metric") == "max" else "99thPercentile"
        pop_columns = [f"loci.{d}_{m}" for d in population_datasets_compared(params.get("require_above_population"))
                       if f"{d}_{m}" in source_columns]
        if pop_columns and f"FirstAffectedAlleleSize_{ot}" in source_columns:
            has_population_data = " OR ".join(f"{c} IS NOT NULL" for c in pop_columns)
            above_terms = " AND ".join(f"({c} IS NULL OR loci.FirstAffectedAlleleSize_{ot} - {c} >= ?)"
                                       for c in pop_columns)
            clauses.append(f"(({has_population_data}) AND {above_terms})")
            sql_params.extend([params["min_repeats_above_population"] - 0.5] * len(pop_columns))
        else:
            clauses.append("1=0")  # no population columns (or no affected allele) -> no locus has a Δpop

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
        if "pLI" in source_columns:
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
        clause, threshold_params = build_pathogenic_threshold_clause("loci", f"FirstAffectedAlleleSize_{ot}")
        clauses.append(clause)
        sql_params.extend(threshold_params)

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
                if f"{col_prefix}AffectedPhenotype_{ot}" in source_columns
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
        # are packed into one "{n}x:{sample_id},..." string. A comma is added at the end so every
        # sample id sits between ":" and ",". Sample ids routinely contain "_", which is ILIKE's
        # single-character wildcard, so escape the metacharacters; otherwise "sample_1" would also
        # match the unrelated "sampleX1".
        kw_clauses = [f"((OutlierSampleIds_{ot} || ',') ILIKE ? ESCAPE '\\')" for _ in params["sample_id_keyword"]]
        clauses.append(f"({' OR '.join(kw_clauses)})")
        sql_params.extend(f"%:{escape_like_wildcards(kw)},%" for kw in params["sample_id_keyword"])

    if "sample_id_like" in params:
        # ILIKE against the whole packed string would also match allele sizes and the "x:"
        # delimiters ("15" matching "15x:sampleA"), so split it into entries and match inside each
        # entry's sample-id field only, which is its second ":"-separated field. Wildcards are
        # escaped for the same reason as in sample_id_keyword above: a keyword such as "sample_1"
        # means a literal underscore, not ILIKE's single-character wildcard.
        kw_clauses = [
            f"len(list_filter(string_split(OutlierSampleIds_{ot}, ','),"
            f" lambda entry: split_part(entry, ':', 2) ILIKE ? ESCAPE '\\')) > 0"
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
                if "gene_id" in source_columns:
                    search_clauses.append("gene_id = ?")
                    search_params.append(value)
            elif kind == "gene_symbol":
                if "GeneTableGeneSymbol" in source_columns:
                    search_clauses.append("GeneTableGeneSymbol ILIKE ?")
                    search_params.append(f"%{value}%")
            else:
                clause, region_params = reference_region_clause(value)
                search_clauses.append(f"({clause})")
                search_params.extend(region_params)
        if locus_id_values:
            # One "lower(LocusId) IN (...)" for every pasted id: a single hash probe per row,
            # which DuckDB can also answer from the column's statistics, where a predicate per
            # id would cost one comparison per id per row. Pasting a few hundred ids is normal
            # here, so the difference is seconds per request.
            # Both sides are lower-cased because the build stores a locus id as the input matrix
            # spelled it and the match is case-insensitive. An IN list holds no wildcards, so a
            # contig name carrying "_" needs no escaping and matches literally.
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
        if "VariationClusterSizeDiff" in source_columns:
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
        if "VariationClusterFilterReason" in source_columns:
            clauses.append("(loci.VariationClusterFilterReason IS NULL "
                           "OR loci.VariationClusterFilterReason != ?)")
            sql_params.append(filter_reason)

    if "chrom" in params:
        clauses.append("Chrom = ?")
        sql_params.append(params["chrom"])

    if "gene_id" in params:
        if "gene_id" in source_columns:
            clauses.append("gene_id = ?")
            sql_params.append(params["gene_id"])
        else:
            clauses.append("1=0")

    if "gene_symbol" in params:
        if "GeneTableGeneSymbol" in source_columns:
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

    # The loci where #Aff counts one sample, or one family's samples (the loci a row of the grouped
    # view expands to). A database without the table cannot say, so nothing matches.
    for param_name, column in (("affected_unsolved_outlier_sample_id", "sample_id"),
                               ("affected_unsolved_outlier_family_id", "family_id")):
        if param_name not in params:
            continue
        if not has_affected_unsolved_outliers_table():
            clauses.append("1=0")
            continue
        clauses.append(f"loci.LocusId IN (SELECT LocusId FROM {AFFECTED_UNSOLVED_OUTLIERS_TABLE} "
                       f"WHERE outlier_type = ? AND {column} = ?)")
        sql_params.extend([ot, params[param_name]])

    return clauses, sql_params


def build_api_query(params):
    """Build the SQL queries from validated API query parameters.

    Returns:
        Tuple (select_query, count_query, motif_count_query, gene_region_count_query,
        sql_params, sql_params_with_pagination).
    """
    source_columns = app.config["DB_COLUMNS_SET"]
    where, sql_params = build_api_where_clause(params)
    motif_count_col = "CanonicalMotif" if "CanonicalMotif" in source_columns else "Motif"
    # gene_region is an optional annotation column. Without it the region summary simply has
    # nothing to report; querying for it anyway would fail the whole /api/v1/loci request.
    has_gene_region = "gene_region" in source_columns

    order_by = build_api_order_by(params, source_columns)

    # One pass over loci. DuckDB is columnar: the filter and the sort read only the columns they
    # name, and the remaining ~180 are read only for the page's rows, so picking the page's
    # LocusIds in a subquery first would only make the table be scanned twice.
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


def split_outlier_entries(outlier_string):
    """Split an OutlierSampleIds string into its entries, skipping entries that will not parse.

    The format is comma-separated "{allele}x:{sample_id}" entries, as
    allele_histograms.convert_sample_ids_to_string writes them.

    Args:
        outlier_string: An OutlierSampleIds_* value, or None.

    Returns:
        List of (allele_size, sample_id) tuples in the string's order, allele_size an int.
    """
    entries = []
    for entry in (outlier_string or "").split(","):
        parts = entry.strip().split(":")
        if len(parts) != 2:
            continue
        try:
            entries.append((int(parts[0].replace("x", "")), parts[1]))
        except ValueError:
            continue
    return entries


def parse_outlier_samples(outlier_string, sample_lookup, affected_lookup, analysis_lookup):
    """Parse an OutlierSampleIds string and enrich with sample metadata.

    The format is comma-separated "{allele}x:{sample_id}" entries (split_outlier_entries).

    Returns:
        list of dicts with allele_size, sample_id, family_id, affected_status,
        analysis_status, sex, ancestry and phenotype_description.
    """
    return [{"allele_size": allele_size, **describe_sample(sample_id, sample_lookup, affected_lookup, analysis_lookup)}
            for allele_size, sample_id in split_outlier_entries(outlier_string)]


def describe_sample(sample_id, sample_lookup, affected_lookup, analysis_lookup):
    """Return the sample-table fields the page shows for a sample, with statuses title-cased.

    Args:
        sample_id: Sample identifier.
        sample_lookup: dict mapping sample_id to row data.
        affected_lookup: dict mapping sample_id to affected status.
        analysis_lookup: dict mapping sample_id to analysis status.

    Returns:
        dict with sample_id, family_id, affected_status, analysis_status, sex, ancestry and
        phenotype_description.
    """
    sample_row = sample_lookup.get(sample_id, {})

    # Prefer the raw sample-row affected_status for display/export so
    # "possibly affected" is preserved; affected_lookup collapses it to
    # "affected" for analysis logic only (the affected_status handling in
    # swim_plot.generate_swim_plot_table mirrors this). Fall back to affected_lookup
    # when the sample has no usable raw status.
    affected_raw = str(sample_row.get("affected_status") or "").strip()
    if affected_raw.lower() in analysis_columns.BLANK_STATUS_VALUES:
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

    return {
        "sample_id": sample_id,
        "family_id": sample_row.get("family_id"),
        "affected_status": affected_display,
        "analysis_status": analysis_display,
        "sex": sample_row.get("sex"),
        "ancestry": sample_row.get("ancestry"),
        "phenotype_description": sample_row.get("phenotype_description"),
    }


def build_sample_details(row_dict, raw_row, lookups, per_outlier_scores, population_gate_cohorts):
    """Build sample details for outliers exceeding the population 99th percentile.

    Args:
        row_dict: dict of the locus's list columns (row_to_list_dict).
        raw_row: The locus's loci-table row, holding its OutlierSampleIds_* columns.
        lookups: app.config["LOOKUPS"].
        per_outlier_scores: dict mapping outlier type to the locus's per-outlier phenotype score
            dicts (see fetch_per_outlier_phenotype_scores_of_loci), {} when it has none.
        population_gate_cohorts: The cohorts whose 99th percentiles the build gated its counts
            and phenotype scores with; an outlier is listed only above the largest of them.

    Returns:
        dict with AllAlleleOutliers / ShortAlleleOutliers / HemizygousOutliers keys
        (only present if qualifying samples exist).
    """
    pop_max = analysis_columns.get_population_p99_threshold(row_dict, population_gate_cohorts) or 0

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


def lowest_catalog_pathogenic_min(disease_locus):
    """Return the smallest PathogenicMin among a variant catalog locus's diseases.

    Args:
        disease_locus: A variant catalog record.

    Returns:
        The smallest PathogenicMin, or None when none of its diseases records one.
    """
    pathogenic_mins = [disease["PathogenicMin"] for disease in disease_locus.get("Diseases", [])
                       if disease.get("PathogenicMin") is not None]
    return min(pathogenic_mins) if pathogenic_mins else None


def compute_known_disease_info(row, lookups):
    """Look up known-disease-locus info for a locus.

    Looks the locus up by LocusId and by coordinates in the variant catalog first, then matches
    it the way the build does (locus_annotations.find_matching_disease_locus). Returns a dict
    with disease info or None when no match is found.
    """
    locus_id = row.get("LocusId", "")
    chrom = row.get("Chrom", "")
    chrom_key = chrom.replace("chr", "") if chrom else ""
    start_0based = row.get("Start0Based")
    end_1based = row.get("End1Based")
    motif = row.get("Motif", "")

    disease_info = (lookups["locus"].get(locus_id)
                    or lookups["locus"].get(f"{chrom_key}-{start_0based}-{end_1based}-{motif}"))
    source = "variant_catalog"
    if not disease_info:
        disease_info, source = find_matching_disease_locus(
            chrom, start_0based, end_1based, motif, lookups["disease_trees"], lookups.get("strchive_trees"))

    if not disease_info:
        return None

    if source == "STRchive":
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

    inheritance_modes = {d["Inheritance"] for d in disease_info.get("Diseases", []) if d.get("Inheritance")}

    return {
        "locus_id": disease_info.get("LocusId"),
        "pathogenic_min": lowest_catalog_pathogenic_min(disease_info),
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

    The locus is matched with locus_annotations.find_matching_disease_locus, the same match the
    build and the locus-detail page use. Walking STRchive too is what keeps the pathogenic-threshold filter from dropping a
    locus that is known only through STRchive, whose locus-detail page does report a pathogenic
    minimum. A variant catalog match stops the search even when it records no minimum, since
    STRchive is only the fallback for loci the catalog does not know.

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
    disease_locus, source = find_matching_disease_locus(
        chrom, start_0based, end_1based, motif, disease_trees, strchive_trees)
    if disease_locus is None:
        return None
    if source == "STRchive":
        return disease_locus.get("pathogenic_min")
    return lowest_catalog_pathogenic_min(disease_locus)


# The per_outlier_phenotype_scores columns each per-outlier score dict holds
PER_OUTLIER_PHENOTYPE_SCORE_COLUMNS = (
    "sample_id", "allele_size", "gene_symbol",
    "gene_phenotype_similarity", "gene_phenotype_overlap_count",
    "n_matching_diseases", "best_matching_disease", "best_disease_inheritance",
    "pairwise_similarity_to_next", "pairwise_shared_count_raw",
    "pairwise_shared_count_ic", "next_sample_id",
)

# The per_locus_phenotype_scores columns each aggregate score dict holds
PER_LOCUS_PHENOTYPE_SCORE_COLUMNS = (
    "num_qualifying_samples", "sum_pairwise_similarity",
    "sum_pairwise_shared_raw", "sum_pairwise_shared_ic",
    "max_gene_phenotype_similarity", "qualifying_sample_ids",
)


def fetch_per_outlier_phenotype_scores_of_loci(conn, locus_ids_query, locus_ids_query_params):
    """Load the per-outlier phenotype scores of the loci a query returns, in one query.

    Used for one locus by the locus-detail page, and for every exported locus by the exports that
    list each locus's outlier samples. Only those loci are loaded: a per_outlier_phenotype_scores
    table can hold millions of rows, which take several GB as Python dicts.

    Args:
        conn: DuckDB database connection.
        locus_ids_query: A SELECT returning a LocusId column, e.g. an export's query over loci.
        locus_ids_query_params: The parameters of locus_ids_query.

    Returns:
        dict mapping locus_id to a dict mapping outlier_type to that type's score dicts
        (PER_OUTLIER_PHENOTYPE_SCORE_COLUMNS, NaN as None), largest allele first. {} when the
        database has no per_outlier_phenotype_scores table.
    """
    if "per_outlier_phenotype_scores" not in duckdb_compat.list_tables(conn):
        return {}

    scores = collections.defaultdict(dict)
    rows = conn.execute(
        f"SELECT locus_id, outlier_type, {', '.join(PER_OUTLIER_PHENOTYPE_SCORE_COLUMNS)} "
        f"FROM per_outlier_phenotype_scores "
        f"WHERE locus_id IN (SELECT LocusId FROM ({locus_ids_query}) AS selected_loci) "
        f"ORDER BY locus_id, outlier_type, allele_size DESC",
        list(locus_ids_query_params),
    ).fetchall()
    for locus_id, outlier_type, *values in rows:
        scores[locus_id].setdefault(outlier_type, []).append({
            column: None if isinstance(value, float) and value != value else value
            for column, value in zip(PER_OUTLIER_PHENOTYPE_SCORE_COLUMNS, values)})
    return dict(scores)


def fetch_phenotype_scores(conn, locus_id):
    """Fetch the phenotype scores of one locus, if the score tables exist.

    Args:
        conn: DuckDB database connection.
        locus_id: The locus's LocusId.

    Returns:
        dict with a 'per_outlier' key (outlier_type -> score dicts, see
        fetch_per_outlier_phenotype_scores_of_loci) and a 'per_locus' key (outlier_type -> that
        type's aggregate scores, PER_LOCUS_PHENOTYPE_SCORE_COLUMNS), each {} when the locus has
        none. {} when the database has no per_outlier_phenotype_scores table.
    """
    tables = duckdb_compat.list_tables(conn)
    if "per_outlier_phenotype_scores" not in tables:
        return {}

    result = {
        "per_outlier": fetch_per_outlier_phenotype_scores_of_loci(
            conn, "SELECT ? AS LocusId", [locus_id]).get(locus_id, {}),
        "per_locus": {},
    }
    if "per_locus_phenotype_scores" in tables:
        rows = conn.execute(
            f"SELECT outlier_type, {', '.join(PER_LOCUS_PHENOTYPE_SCORE_COLUMNS)} "
            f"FROM per_locus_phenotype_scores WHERE locus_id = ?",
            (locus_id,),
        ).fetchall()
        for outlier_type, *values in rows:
            result["per_locus"][outlier_type] = {
                column: None if isinstance(value, float) and value != value else value
                for column, value in zip(PER_LOCUS_PHENOTYPE_SCORE_COLUMNS, values)}
    return result


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
        has_mendelian=app.config.get("HAS_MENDELIAN", False),
        mendelian_warnings=app.config.get("MENDELIAN_WARNINGS", {}),
        outlier_warnings=app.config.get("OUTLIER_WARNINGS", {}),
        initial_require_above_population=app.config.get("INITIAL_REQUIRE_ABOVE_POPULATION", ""),
        population_gate_cohorts=app.config.get("POPULATION_GATE_COHORTS", ()),
    )


@app.route("/swim")
def swim():
    """Serve the swim plot UI."""
    template = jinja2_env.get_template("swim_plot_template.html")
    return template.render(
        has_mendelian=app.config.get("HAS_MENDELIAN", False),
        mendelian_warnings=app.config.get("MENDELIAN_WARNINGS", {}),
        outlier_warnings=app.config.get("OUTLIER_WARNINGS", {}),
    )


@app.route("/sample_qc")
def sample_qc():
    """Serve the sample QC page."""
    template = jinja2_env.get_template("sample_qc_template.html")
    return template.render(
        has_mendelian=app.config.get("HAS_MENDELIAN", False),
        mendelian_warnings=app.config.get("MENDELIAN_WARNINGS", {}),
        outlier_warnings=app.config.get("OUTLIER_WARNINGS", {}),
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
        min_trios_per_center=MIN_TRIOS_PER_CENTER,
    )


# ---------------------------------------------------------------------------
# Annotation CRUD endpoints
# ---------------------------------------------------------------------------


@app.route("/api/v1/annotations/tags", methods=["GET"])
def list_all_tags():
    """List all existing user tags (for autocomplete/dropdown)."""
    return jsonify({"tags": app.config["ANNOTATIONS"]["all_tags"]})


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


# Flask serves requests on several threads, so two annotation writes can overlap, and DuckDB raises
# on a second concurrent writer instead of waiting for it. Every write to the annotations database
# therefore runs under this one process-wide lock, which is enough because that database is only
# ever written from here.
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
        selected_alleles = fetch_selected_sample_or_family_alleles(conn, params, [dict(row)["LocusId"] for row in rows])
    finally:
        conn.close()

    motif_counts = {row[0]: row[1] for row in motif_rows if row[0]}
    gene_region_counts = {row[0]: row[1] for row in gene_region_rows if row[0]}

    total_pages = math.ceil(total / params["page_size"]) if total > 0 else 0
    results = [row_to_list_dict(row, ot) for row in rows]
    if selected_alleles is not None:
        for r in results:
            r["SelectedSampleOrFamilyAlleles"] = selected_alleles.get(r["LocusId"], [])

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

    return jsonify({
        "total": total,
        "page": params["page"],
        "page_size": params["page_size"],
        "total_pages": total_pages,
        "outlier_type": ot,
        "filters_applied": build_filters_applied(params),
        "results": results,
        "motif_counts": motif_counts,
        "gene_region_counts": gene_region_counts,
    })


def fetch_selected_sample_or_family_alleles(conn, params, locus_ids):
    """The selected sample's (or family's) own alleles at each locus of a restricted locus list.

    When /api/v1/loci is restricted to one sample or family (see build_api_where_clause), the
    page shows that sample's allele next to the locus's top affected allele.

    Args:
        conn: Open connection to the database.
        params: dict of validated query parameters.
        locus_ids: LocusIds of the page of results.

    Returns:
        None when the list is not restricted to a sample or family, else a dict mapping LocusId
        to a list of {"sample_id", "allele_size"} dicts, largest allele first.
    """
    for param_name, column in (("affected_unsolved_outlier_sample_id", "sample_id"),
                               ("affected_unsolved_outlier_family_id", "family_id")):
        if param_name in params:
            break
    else:
        return None
    if not locus_ids or not has_affected_unsolved_outliers_table():
        return {}
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    rows = conn.execute(
        f"SELECT LocusId, sample_id, allele_size FROM {AFFECTED_UNSOLVED_OUTLIERS_TABLE} "
        f"WHERE outlier_type = ? AND {column} = ? AND LocusId IN ({','.join('?' * len(locus_ids))}) "
        f"ORDER BY LocusId, allele_size DESC, sample_id",
        [ot, params[param_name]] + list(locus_ids)).fetchall()
    alleles = collections.defaultdict(list)
    for locus_id, sample_id, allele_size in rows:
        alleles[locus_id].append({"sample_id": sample_id, "allele_size": allele_size})
    return alleles


# group_sort values of /api/v1/outlier_samples: the groups CTE column each one sorts by, and the
# direction it sorts in unless group_sort_direction says otherwise
GROUPED_SORTS = {
    "num_loci": ("num_loci", "desc"),
    "num_disease_gene_loci": ("num_disease_gene_loci", "desc"),
    "best_gene_phenotype_similarity": ("best_gene_phenotype_similarity", "desc"),
    "largest_allele_size": ("largest_allele_size", "desc"),
    "group_id": ("group_id", "asc"),
}

# Filters that describe a sample rather than a locus. The grouped view applies them to each row's
# own sample (its id, or its own phenotype), where the locus view matches a locus's outlier list
# or the phenotypes of its outlier samples.
SAMPLE_LEVEL_FILTERS = ("sample_id_keyword", "sample_id_like", "phenotype_keyword")


def find_samples_passing_sample_level_filters(params, sample_lookup):
    """The samples that pass the grouped view's sample-level filters.

    Matching ignores case. sample_id_keyword matches whole sample ids, sample_id_like any part of
    one, and phenotype_keyword any part of the sample's phenotype description. Keywords of one
    filter are ORed, and the filters are ANDed.

    Args:
        params: dict of validated query parameters.
        sample_lookup: dict mapping sample_id to its sample-table row, or to a dict holding only
            its family_id for a sample in the outlier table that has no sample-table row.

    Returns:
        Sorted list of sample ids, or None when no sample-level filter is set.
    """
    if not any(name in params for name in SAMPLE_LEVEL_FILTERS):
        return None
    exact_ids = {kw.lower() for kw in params.get("sample_id_keyword", [])}
    passing = []
    for sample_id, sample_row in sample_lookup.items():
        if exact_ids and sample_id.lower() not in exact_ids:
            continue
        if "sample_id_like" in params and not any(kw.lower() in sample_id.lower()
                                                  for kw in params["sample_id_like"]):
            continue
        if "phenotype_keyword" in params:
            phenotype = sample_row.get("phenotype_description")
            phenotype = phenotype.lower() if isinstance(phenotype, str) else ""
            if not any(kw.lower() in phenotype for kw in params["phenotype_keyword"]):
                continue
        passing.append(sample_id)
    return sorted(passing)


def validate_grouped_params():
    """Validate the grouped view's parameters: every /api/v1/loci one, plus group_by, group_sort and
    group_sort_direction ("asc" or "desc", by default group_sort's own direction in GROUPED_SORTS).

    Returns:
        Tuple (params, group_by, group_sort, group_sort_direction, error_response).
        error_response is not None when the request is invalid or the database has no
        AFFECTED_UNSOLVED_OUTLIERS_TABLE.
    """
    def invalid(response):
        return None, None, None, None, response

    params, error = validate_params()
    if error:
        return invalid(error)
    group_by = request.args.get("group_by", "")
    if group_by not in ("sample", "family"):
        return invalid((jsonify({"error": "Invalid parameter", "detail": f"group_by must be one of: sample, family; got '{group_by}'"}), 400))
    group_sort = request.args.get("group_sort") or "num_loci"
    if group_sort not in GROUPED_SORTS:
        return invalid((jsonify({"error": "Invalid parameter", "detail": f"group_sort must be one of: {', '.join(GROUPED_SORTS)}; got '{group_sort}'"}), 400))
    group_sort_direction = request.args.get("group_sort_direction") or GROUPED_SORTS[group_sort][1]
    if group_sort_direction not in ("asc", "desc"):
        return invalid((jsonify({"error": "Invalid parameter", "detail": f"group_sort_direction must be one of: asc, desc; got '{group_sort_direction}'"}), 400))
    if not has_affected_unsolved_outliers_table():
        return invalid((jsonify({"error": "Not available", "detail": f"The database has no {AFFECTED_UNSOLVED_OUTLIERS_TABLE} "
                                 f"table. Rebuild it with trails.py and restart the server."}), 400))
    return params, group_by, group_sort, group_sort_direction, None


def build_grouped_ctes(params, group_by):
    """The SQL common to the grouped view and its export.

    A sample belongs under a locus when #Aff counts it there (AFFECTED_UNSOLVED_OUTLIERS_TABLE):
    affected status not "unaffected", analysis status not solved, and an allele above every
    unaffected or solved sample's and above the population 99th percentile. A family's rows are
    its samples' rows, so a sample without a family_id appears only when grouping by sample.
    Sample-level filters (SAMPLE_LEVEL_FILTERS) pick samples instead of loci, or in the family
    view the families with at least one such sample. A group's num_loci is therefore exactly the
    number of loci /api/v1/loci returns for it with the same filters, minus those, plus
    affected_unsolved_outlier_sample_id (or _family_id).

    Args:
        params: dict of validated query parameters.
        group_by: "sample" or "family".

    Returns:
        Tuple (ctes, ctes_params, group_column). ctes is a "WITH ..." prefix that defines
        matching_loci (the loci matching the locus filters, with their disease-gene flag, motif
        and gene region), pairs (the samples #Aff counts at them that pass the sample-level
        filters, with family_id and allele_size) and groups (one row per sample or family:
        group_id, num_loci, num_disease_gene_loci, largest_allele_size and
        best_gene_phenotype_similarity). group_column is the pairs column holding the group id.
    """
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    lookups = app.config["LOOKUPS"]
    where, where_params = build_api_where_clause({k: v for k, v in params.items() if k not in SAMPLE_LEVEL_FILTERS})
    source_columns = app.config["DB_COLUMNS_SET"]
    disease_gene_column = "loci.IsInMendelianGene" if "IsInMendelianGene" in source_columns else "NULL"
    motif_count_col = "CanonicalMotif" if "CanonicalMotif" in source_columns else "Motif"
    gene_region_column = "loci.gene_region" if "gene_region" in source_columns else "NULL"
    group_column = "pairs.sample_id" if group_by == "sample" else "pairs.family_id"

    # The (locus, sample) rows behind the view: the samples #Aff counts at the matching loci
    pair_conditions = ["p.outlier_type = ?"]
    pair_params = [ot]
    if group_by == "family":
        pair_conditions.append("p.family_id IS NOT NULL")
    # The sample-level filters see every sample of the outlier table, with the family the build
    # recorded for it, and every sample-table row (the only source of phenotypes)
    sample_rows = {sample_id: {"family_id": family_id}
                   for sample_id, family_id in app.config.get("OUTLIER_TABLE_FAMILY_ID_BY_SAMPLE", {}).items()}
    sample_rows.update(lookups["sample"])
    passing_samples = find_samples_passing_sample_level_filters(params, sample_rows)
    if passing_samples is not None:
        # In the family view the filters pick whole families: those with a sample that passes
        if group_by == "family":
            passing_column = "p.family_id"
            passing_ids = sorted({family_id for family_id in (sample_rows[s].get("family_id") for s in passing_samples)
                                  if isinstance(family_id, str) and family_id})
        else:
            passing_column, passing_ids = "p.sample_id", passing_samples
        if passing_ids:
            pair_conditions.append(f"{passing_column} IN ({','.join('?' * len(passing_ids))})")
            pair_params.extend(passing_ids)
        else:
            pair_conditions.append("1=0")
    if "per_outlier_phenotype_scores" in app.config.get("DB_TABLES", set()):
        best_similarity = "MAX(s.gene_phenotype_similarity)"
        phenotype_join = ("LEFT JOIN per_outlier_phenotype_scores s ON s.locus_id = pairs.LocusId "
                          "AND s.sample_id = pairs.sample_id AND s.outlier_type = ?")
        phenotype_join_params = [ot]
    else:
        best_similarity, phenotype_join, phenotype_join_params = "NULL", "", []
    ctes = (f"WITH matching_loci AS (SELECT loci.LocusId, {disease_gene_column} AS is_in_disease_gene, "
            f"loci.{motif_count_col} AS motif, {gene_region_column} AS gene_region FROM loci WHERE {where}), "
            f"pairs AS (SELECT p.LocusId, p.sample_id, p.family_id, p.allele_size, m.is_in_disease_gene "
            f"FROM {AFFECTED_UNSOLVED_OUTLIERS_TABLE} p JOIN matching_loci m ON m.LocusId = p.LocusId "
            f"WHERE {' AND '.join(pair_conditions)}), "
            f"groups AS (SELECT {group_column} AS group_id, "
            f"COUNT(DISTINCT pairs.LocusId) AS num_loci, "
            f"COUNT(DISTINCT CASE WHEN pairs.is_in_disease_gene = 1 THEN pairs.LocusId END) AS num_disease_gene_loci, "
            f"MAX(pairs.allele_size) AS largest_allele_size, "
            f"{best_similarity} AS best_gene_phenotype_similarity "
            f"FROM pairs {phenotype_join} GROUP BY {group_column}) ")
    return ctes, where_params + pair_params + phenotype_join_params, group_column


def grouped_order_by(group_sort, group_sort_direction):
    """ORDER BY terms (over the groups CTE's columns) that put the groups in group_sort's order.

    A group with no value (no gene-phenotype score) goes last in either direction. Ties are
    broken by #Loci, then by group id.

    Args:
        group_sort: A GROUPED_SORTS key.
        group_sort_direction: "asc" or "desc".

    Returns:
        The ORDER BY terms, without the "ORDER BY" keyword.
    """
    column = GROUPED_SORTS[group_sort][0]
    return f"{column} {group_sort_direction.upper()} NULLS LAST, num_loci DESC, group_id"


def describe_groups(conn, ctes, ctes_params, group_column, group_rows, group_by, params):
    """Complete rows of the groups CTE with what the grouped view shows for each.

    Adds each group's largest allele's locus, gene and Δpop (compute_delta_population), and
    the sample's sample-table fields, or a family's samples among the pairs.

    Args:
        conn: Open connection to the database.
        ctes, ctes_params, group_column: from build_grouped_ctes.
        group_rows: dicts of groups CTE rows.
        group_by: "sample" or "family".
        params: dict of validated query parameters.

    Returns:
        List of dicts, one per group row, in the same order.
    """
    lookups = app.config["LOOKUPS"]
    group_ids = [row["group_id"] for row in group_rows]
    largest_allele_loci = {}
    members = collections.defaultdict(set)
    if group_ids:
        in_groups = f"{group_column} IN ({','.join('?' * len(group_ids))})"
        # Each group's largest allele and the locus it is at (ties go to the smallest LocusId)
        for group_id, locus_id in conn.execute(
                f"{ctes}SELECT group_id, LocusId FROM (SELECT {group_column} AS group_id, pairs.LocusId, "
                f"ROW_NUMBER() OVER (PARTITION BY {group_column} ORDER BY pairs.allele_size DESC, pairs.LocusId) AS rank_in_group "
                f"FROM pairs WHERE {in_groups}) ranked WHERE rank_in_group = 1",
                ctes_params + group_ids).fetchall():
            largest_allele_loci[group_id] = locus_id
        if group_by == "family":
            for family_id, sample_id in conn.execute(
                    f"{ctes}SELECT DISTINCT pairs.family_id, pairs.sample_id FROM pairs WHERE {in_groups}",
                    ctes_params + group_ids).fetchall():
                members[family_id].add(sample_id)
    locus_rows = {}
    if largest_allele_loci:
        locus_ids = sorted(set(largest_allele_loci.values()))
        for row in conn.execute(f"SELECT * FROM loci WHERE LocusId IN ({','.join('?' * len(locus_ids))})", locus_ids):
            row = dict(row)
            locus_rows[row["LocusId"]] = row

    results = []
    for row in group_rows:
        locus_id = largest_allele_loci.get(row["group_id"])
        locus_row = locus_rows.get(locus_id, {})
        result = {
            "group_id": row["group_id"],
            "num_loci": row["num_loci"],
            "num_disease_gene_loci": row["num_disease_gene_loci"],
            "best_gene_phenotype_similarity": row["best_gene_phenotype_similarity"],
            "largest_allele_size": row["largest_allele_size"],
            "largest_allele_locus_id": locus_id,
            "largest_allele_gene_symbol": locus_row.get("GeneTableGeneSymbol"),
            "largest_allele_delta_population": compute_delta_population(row["largest_allele_size"], locus_row, params),
        }
        if group_by == "sample":
            result.update(describe_sample(row["group_id"], lookups["sample"], lookups["affected"], lookups["analysis"]))
            # The family the build recorded, which the family view groups by, rather than the sample
            # table's, which is empty for a sample with no row there
            result["family_id"] = app.config.get("OUTLIER_TABLE_FAMILY_ID_BY_SAMPLE", {}).get(row["group_id"])
        else:
            result["family_id"] = row["group_id"]
            result["members"] = [describe_sample(sample_id, lookups["sample"], lookups["affected"], lookups["analysis"])
                                 for sample_id in sorted(members[row["group_id"]])]
        results.append(result)
    return results


@app.route("/api/v1/outlier_samples")
def get_outlier_samples():
    """The loci matching the filters, grouped by affected unsolved sample or by family.

    Takes every /api/v1/loci parameter, plus group_by ("sample" or "family"), group_sort (a
    GROUPED_SORTS key, default "num_loci") and group_sort_direction ("asc" or "desc"). See
    build_grouped_ctes for which samples a group holds and how the filters apply.
    """
    params, group_by, group_sort, group_sort_direction, error = validate_grouped_params()
    if error:
        return error
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    ctes, ctes_params, group_column = build_grouped_ctes(params, group_by)
    totals_query = f"{ctes}SELECT COUNT(DISTINCT {group_column}), COUNT(DISTINCT pairs.LocusId) FROM pairs"
    page_query = f"{ctes}SELECT * FROM groups ORDER BY {grouped_order_by(group_sort, group_sort_direction)} LIMIT ? OFFSET ?"
    # The charts count the loci the rows span
    chart_loci = "FROM matching_loci WHERE LocusId IN (SELECT LocusId FROM pairs)"
    motif_count_query = f"{ctes}SELECT motif, COUNT(*) AS count {chart_loci} GROUP BY motif ORDER BY count DESC, motif"
    gene_region_count_query = f"{ctes}SELECT gene_region, COUNT(*) AS count {chart_loci} GROUP BY gene_region"

    conn = get_db()
    try:
        total_groups, total_loci = conn.execute(totals_query, ctes_params).fetchone()
        group_rows = [dict(row) for row in conn.execute(
            page_query, ctes_params + [params["page_size"], (params["page"] - 1) * params["page_size"]])]
        motif_rows = conn.execute(motif_count_query, ctes_params).fetchall()
        gene_region_rows = conn.execute(gene_region_count_query, ctes_params).fetchall()
        results = describe_groups(conn, ctes, ctes_params, group_column, group_rows, group_by, params)
    finally:
        conn.close()

    return jsonify({
        "group_by": group_by,
        "group_sort": group_sort,
        "group_sort_direction": group_sort_direction,
        "total": total_groups,
        "total_loci": total_loci,
        "page": params["page"],
        "page_size": params["page_size"],
        "total_pages": math.ceil(total_groups / params["page_size"]) if total_groups else 0,
        "outlier_type": ot,
        "results": results,
        "motif_counts": {row[0]: row[1] for row in motif_rows if row[0]},
        "gene_region_counts": {row[0]: row[1] for row in gene_region_rows if row[0]},
    })


# format values of /api/v1/outlier_samples/export: one row per sample (or family), or one row
# per sample (or family) and locus
GROUPED_EXPORT_FORMATS = ("groups_tsv", "group_loci_tsv")

# Sample-table fields the grouped exports write for each sample, after the group id
GROUPED_EXPORT_SAMPLE_FIELDS = ["sex", "affected_status", "analysis_status", "phenotype_description"]


def format_tsv_line(values):
    """One TSV line: a missing value (None or NaN) as an empty field, tabs and newlines as spaces."""
    return "\t".join("" if v is None or (isinstance(v, float) and math.isnan(v))
                     else str(v).replace("\t", " ").replace("\n", " ") for v in values) + "\n"


@app.route("/api/v1/outlier_samples/export")
def export_outlier_samples():
    """Download the grouped view as gzip-compressed TSV.

    Takes the /api/v1/outlier_samples parameters (page and page_size are ignored) and format:
      groups_tsv: one row per sample (or family), in the grouped view's order, with the columns
        the view shows.
      group_loci_tsv: one row per sample and locus (the pairs behind the view), groups in the
        view's order and, within a group, the largest allele first. Each row carries the
        sample's sample-table fields, its own allele and that allele's Δpop, then the columns of
        the locus TSV export.
    """
    params, group_by, group_sort, group_sort_direction, error = validate_grouped_params()
    if error:
        return error
    order_by = grouped_order_by(group_sort, group_sort_direction)
    export_format = request.args.get("format", "")
    if export_format not in GROUPED_EXPORT_FORMATS:
        return jsonify({"error": "Invalid parameter", "detail": f"format must be one of: {', '.join(GROUPED_EXPORT_FORMATS)}; got '{export_format}'"}), 400
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]
    ctes, ctes_params, group_column = build_grouped_ctes(params, group_by)
    group_noun = "samples" if group_by == "sample" else "families"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    lookups = app.config["LOOKUPS"]

    conn = get_db()
    if export_format == "groups_tsv":
        try:
            group_rows = [dict(row) for row in conn.execute(
                f"{ctes}SELECT * FROM groups ORDER BY {order_by}", ctes_params)]
            groups = describe_groups(conn, ctes, ctes_params, group_column, group_rows, group_by, params)
        finally:
            conn.close()
        group_columns = ["num_loci", "num_disease_gene_loci", "best_gene_phenotype_similarity", "largest_allele_size",
                         "largest_allele_delta_population", "largest_allele_locus_id", "largest_allele_gene_symbol"]
        if group_by == "sample":
            header = ["sample_id", "family_id"] + GROUPED_EXPORT_SAMPLE_FIELDS + group_columns
            lines = [format_tsv_line([g[c] for c in header]) for g in groups]
        else:
            header = ["family_id", "num_samples", "samples", "phenotypes"] + group_columns
            lines = [format_tsv_line(
                [g["family_id"], len(g["members"]),
                 "; ".join(f"{m['sample_id']} ({m['sex'] if isinstance(m['sex'], str) else ''}, "
                           f"{m['affected_status']}/{m['analysis_status']})" for m in g["members"]),
                 " | ".join(dict.fromkeys(m["phenotype_description"] for m in g["members"]
                                          if isinstance(m["phenotype_description"], str) and m["phenotype_description"]))]
                + [g[c] for c in group_columns]) for g in groups]
        filename = f"tr_outliers.{FILENAME_OUTLIER_TYPE_LABELS[params['outlier_type']]}.{len(groups)}_{group_noun}.{timestamp}.tsv.gz"
        return Response(gzip_stream(iter(["\t".join(header) + "\n"] + lines)), mimetype="application/gzip",
                        headers={"Content-Disposition": f'attachment; filename="{filename}"'})

    # group_loci_tsv, streamed: the pairs, ordered by their group's rank in the grouped view
    pairs_query = (f"{ctes}, ranked AS (SELECT group_id, ROW_NUMBER() OVER (ORDER BY {order_by}) AS group_rank "
                   f"FROM groups) "
                   f"SELECT pairs.sample_id AS pair_sample_id, pairs.family_id AS pair_family_id, "
                   f"pairs.allele_size AS pair_allele_size, loci.* "
                   f"FROM pairs JOIN ranked ON ranked.group_id = {group_column} JOIN loci ON loci.LocusId = pairs.LocusId "
                   f"ORDER BY ranked.group_rank, pairs.allele_size DESC, pairs.LocusId, pairs.sample_id")
    id_columns = ["sample_id", "family_id"] if group_by == "sample" else ["family_id", "sample_id"]
    try:
        num_pairs = conn.execute(f"{ctes}SELECT COUNT(*) FROM pairs", ctes_params).fetchone()[0]
    except Exception:
        conn.close()
        raise

    def generate():
        try:
            yield "\t".join(id_columns + GROUPED_EXPORT_SAMPLE_FIELDS + ["sample_allele_size", "sample_allele_delta_population"]
                            + LIST_COLUMNS_STATIC + LIST_COLUMNS_OT_SPECIFIC) + "\n"
            for row in conn.execute(pairs_query, ctes_params):
                row_dict = dict(row)
                sample = describe_sample(row_dict["pair_sample_id"], lookups["sample"], lookups["affected"], lookups["analysis"])
                # The family the rows were grouped by, as the build recorded it
                sample["family_id"] = row_dict["pair_family_id"]
                allele_size = row_dict["pair_allele_size"]
                yield format_tsv_line(
                    [sample[c] for c in id_columns + GROUPED_EXPORT_SAMPLE_FIELDS]
                    + [allele_size, compute_delta_population(allele_size, row_dict, params)]
                    + list(row_to_list_dict(row, ot).values()))
        finally:
            conn.close()
    filename = (f"tr_outliers.{FILENAME_OUTLIER_TYPE_LABELS[params['outlier_type']]}.{num_pairs}_"
                f"{'sample' if group_by == 'sample' else 'family'}_locus_pairs.{timestamp}.tsv.gz")
    return Response(gzip_stream(generate()), mimetype="application/gzip",
                    headers={"Content-Disposition": f'attachment; filename="{filename}"'})


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
    """BED6 format, streamed with no header: chrom, start, end, LocusId as the name, score 0, strand ".".

    Score and strand hold the values BED defines for "none": a tandem repeat locus has no score,
    and it has no strand. The LocusId carries the motif.
    """
    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
    except Exception:
        conn.close()
        raise
    def generate():
        try:
            for row in conn.execute(query, sql_params):
                r = dict(row)
                yield f"{r['Chrom']}\t{r['Start0Based']}\t{r['End1Based']}\t{r['LocusId']}\t0\t.\n"
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
                yield format_tsv_line(row_to_list_dict(row, ot).values())
        finally:
            conn.close()
    filename = build_export_filename(params["outlier_type"], total, "tsv.gz")
    return Response(gzip_stream(generate()), mimetype="application/gzip",
                    headers={"Content-Disposition": f'attachment; filename="{filename}"'})


def row_to_export_dict(row, ot, tags_by_locus):
    """Return a loci row's list columns plus its tags, as the JSON exports write them.

    Args:
        row: A loci-table row.
        ot: The outlier type's column suffix, e.g. "AllAlleles".
        tags_by_locus: dict mapping LocusId to the locus's tags.

    Returns:
        The row_to_list_dict dict, with a "Tags" key holding the locus's tags comma-separated.
    """
    row_dict = row_to_list_dict(row, ot)
    row_dict["Tags"] = ",".join(tags_by_locus.get(row_dict.get("LocusId", ""), []))
    return row_dict


def export_json(conn, query, count_query, sql_params, ot, params):
    """JSON format: {"metadata": {...}, "loci": [...]}.

    Every locus is always enriched with a "Samples" key listing every qualifying
    outlier sample (allele size, sample/family id, sex, ancestry, affected/analysis
    status, phenotype, and phenotype-similarity scores), whatever the size of the
    result set. The exported loci's phenotype scores are bulk-loaded
    once up front (fetch_per_outlier_phenotype_scores_of_loci) rather than per row. Rows
    still stream from the cursor.
    """
    tags_by_locus = app.config["ANNOTATIONS"]["tags"]
    lookups = app.config["LOOKUPS"]
    population_gate_cohorts = app.config["POPULATION_GATE_COHORTS"]

    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
    except Exception:
        conn.close()
        raise

    metadata = {
        "outlier_type": params["outlier_type"],
        "sort_by": params.get("sort_by") or ["count"],
        "filters_applied": build_filters_applied(params),
        "total_loci": total,
        "generated_at": datetime.now().isoformat(),
    }

    def generate():
        try:
            per_outlier_scores_by_locus = fetch_per_outlier_phenotype_scores_of_loci(conn, query, sql_params)
            yield '{\n  "metadata": ' + json.dumps(metadata, indent=2).replace("\n", "\n  ") + ',\n  "loci": [\n'
            first = True
            for raw_row in conn.execute(query, sql_params):
                row_dict = row_to_export_dict(raw_row, ot, tags_by_locus)
                sample_details = build_sample_details(
                    row_dict, raw_row, lookups, per_outlier_scores_by_locus.get(row_dict.get("LocusId", ""), {}),
                    population_gate_cohorts)
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


def compute_plot_thresholds(row_dict, all_allele_outliers_str, short_allele_outliers_str, affected_lookup,
                            population_gate_cohorts):
    """Compute PlotReadVisualization thresholds for the ExpansionHunter export.

    A locus's largest outlier earns a threshold when it is affected and above every unaffected
    outlier; for the long allele it also has to be above the largest 99th percentile of the
    population_gate_cohorts, the cohorts the build gated its counts with.

    Returns a dict with LongAllele/ShortAllele thresholds, or None if no criteria met.
    """
    def find_threshold(outliers, lookup, comparison_value):
        if not outliers:
            return None
        allele_size, sample_id = outliers[0]
        if lookup.get(sample_id, "") == "affected" and allele_size > comparison_value:
            return allele_size
        return None

    long_outliers = split_outlier_entries(all_allele_outliers_str)
    short_outliers = split_outlier_entries(short_allele_outliers_str)
    if not long_outliers and not short_outliers:
        return None

    long_unaffected = [a for a, s in long_outliers if affected_lookup.get(s, "") != "affected"]
    short_unaffected = [a for a, s in short_outliers if affected_lookup.get(s, "") != "affected"]
    max_long_unaffected = max(long_unaffected) if long_unaffected else 0
    max_short_unaffected = max(short_unaffected) if short_unaffected else 0

    long_comparison = max(
        max_long_unaffected,
        analysis_columns.get_population_p99_threshold(row_dict, population_gate_cohorts) or 0,
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
    """ExpansionHunter variant catalog format, streamed.

    Every locus carries its read-visualization thresholds (compute_plot_thresholds) and its
    qualifying outlier samples (build_sample_details), whatever the size of the export. As in
    export_json, the exported loci's phenotype scores are loaded once up front rather than per row.
    """
    tags_by_locus = app.config["ANNOTATIONS"]["tags"]
    lookups = app.config["LOOKUPS"]
    population_gate_cohorts = app.config["POPULATION_GATE_COHORTS"]

    try:
        total = conn.execute(count_query, sql_params).fetchone()[0]
    except Exception:
        conn.close()
        raise

    def generate():
        try:
            per_outlier_scores_by_locus = fetch_per_outlier_phenotype_scores_of_loci(conn, query, sql_params)
            yield "[\n"
            first = True
            for raw_row in conn.execute(query, sql_params):
                row_dict = row_to_export_dict(raw_row, ot, tags_by_locus)
                row_dict["ReferenceRegion"] = f"{row_dict['Chrom']}:{row_dict['Start0Based']}-{row_dict['End1Based']}"
                row_dict["LocusStructure"] = f"({row_dict['Motif']})*"
                row_dict["VariantType"] = "Repeat"
                # Every loci column, since a database may lack an outlier type's column
                raw_row_columns = dict(raw_row)
                thresholds = compute_plot_thresholds(
                    row_dict,
                    raw_row_columns.get("OutlierSampleIds_AllAlleles"),
                    raw_row_columns.get("OutlierSampleIds_ShortAlleles"),
                    lookups["affected"],
                    population_gate_cohorts,
                )
                if thresholds:
                    row_dict["PlotReadVisualization"] = [
                        {"If": allele, "Is": ">=", "Threshold": thresholds[allele]}
                        for allele in ("LongAllele", "ShortAllele") if allele in thresholds
                    ]
                sample_details = build_sample_details(
                    row_dict, raw_row, lookups, per_outlier_scores_by_locus.get(row_dict.get("LocusId", ""), {}),
                    population_gate_cohorts)
                if sample_details:
                    row_dict["Samples"] = sample_details
                yield ("" if first else ",\n") + json.dumps(row_dict, indent=2)
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
    for cohort in POPULATION_DATASETS:
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
        "outlier_samples": outlier_samples,
        "annotations": {
            "known_disease_locus": known_disease,
            "population_thresholds": population_thresholds,
        },
        "user_annotations": user_annotations,
        "phenotype_scores": phenotype_scores,
    })


# Filters the swim plot applies to each allele it draws, comparing that allele's own size, where
# /api/v1/loci compares a locus's largest affected allele
SWIM_PLOT_ALLELE_LEVEL_FILTERS = ("min_expansion", "apply_pathogenic_threshold", "require_above_population")


@app.route("/api/v1/swim_plot_data")
def get_swim_plot_data():
    """Return swim-plot data: one entry per outlier sample per motif category.

    Takes the /api/v1/loci parameters (see validate_params), plus require_affected_above_unaffected,
    which keeps only the affected alleles above every unaffected allele at their locus;
    mendelian_only is another name for mendelian_genes_only. Each row is one allele, so:

      * the sample-level filters (SAMPLE_LEVEL_FILTERS) keep the rows of the matching samples;
      * the SWIM_PLOT_ALLELE_LEVEL_FILTERS compare each row's own allele, where /api/v1/loci
        compares the locus's largest affected allele;
      * every other filter keeps the rows of the loci /api/v1/loci returns for the parameters
        without those two kinds (build_locus_filter_clauses). min_repeats_above_population
        therefore compares every population cohort, as /api/v1/loci does without
        require_above_population.
    """
    params, error = validate_params()
    if error:
        return error
    if is_true_flag("mendelian_only"):
        params["mendelian_genes_only"] = True
    ot = OUTLIER_TYPE_MAP[params["outlier_type"]]

    row_clauses = []
    row_params = []

    if is_true_flag("require_affected_above_unaffected"):
        row_clauses.append("""(
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

    if "min_expansion" in params:
        row_clauses.append("swim_plot.allele_size >= swim_plot.NumRepeatsInReference + ?")
        row_params.append(params["min_expansion"])

    if params.get("apply_pathogenic_threshold"):
        clause, threshold_params = build_pathogenic_threshold_clause("swim_plot", "allele_size")
        row_clauses.append(clause)
        row_params.extend(threshold_params)

    if "require_above_population" in params:
        clause = build_require_above_population_clause(params, "swim_plot", "allele_size")
        if clause:
            row_clauses.append(clause)

    if "phenotype_keyword" in params:
        row_clauses.append("(" + " OR ".join("swim_plot.phenotype_description ILIKE ?"
                                             for _ in params["phenotype_keyword"]) + ")")
        row_params.extend(f"%{kw}%" for kw in params["phenotype_keyword"])

    if "sample_id_keyword" in params:
        row_clauses.append("(" + " OR ".join("lower(swim_plot.sample_id) = lower(?)"
                                             for _ in params["sample_id_keyword"]) + ")")
        row_params.extend(params["sample_id_keyword"])

    if "sample_id_like" in params:
        # Escaped like the loci list's sample_id_like, so the same keyword selects the same samples
        # in both views.
        row_clauses.append("(" + " OR ".join("swim_plot.sample_id ILIKE ? ESCAPE '\\'"
                                             for _ in params["sample_id_like"]) + ")")
        row_params.extend(f"%{escape_like_wildcards(kw)}%" for kw in params["sample_id_like"])

    # The locus filters, answered from loci. swim_plot carries no coordinates or variation-cluster
    # columns, and this endpoint runs one query per motif category, so the matching loci are picked
    # once, into the swim_plot_matching_loci temp table, after the connection is opened.
    locus_filter_clauses, locus_filter_params = build_locus_filter_clauses(
        {k: v for k, v in params.items() if k not in SAMPLE_LEVEL_FILTERS + SWIM_PLOT_ALLELE_LEVEL_FILTERS})
    if locus_filter_clauses:
        row_clauses.append("swim_plot.LocusId IN (SELECT LocusId FROM swim_plot_matching_loci)")

    extra_where = (" AND " + " AND ".join(row_clauses)) if row_clauses else ""
    # "Unknown" is the motif_category swim_plot assigns rows whose MotifSize is
    # missing (swim_plot._compute_motif_category); include it so those rows are
    # not silently dropped from the endpoint.
    categories = [f"{i}bp" for i in range(1, 25)] + ["25+bp", "Unknown"]
    all_data = []

    conn = get_db()
    try:
        if "swim_plot" not in duckdb_compat.list_tables(conn):
            return jsonify({"error": "swim_plot table not found", "detail": "The database has no swim_plot table."}), 500

        if locus_filter_clauses:
            # A read-only main database still allows temp tables (they live in memory, in this
            # connection's own temp schema), so this costs one scan of loci per request instead
            # of one per category.
            conn.execute(
                "CREATE OR REPLACE TEMP TABLE swim_plot_matching_loci AS SELECT LocusId FROM loci "
                f"WHERE {' AND '.join(locus_filter_clauses)}",
                locus_filter_params)

        # LocusId and sample_id follow allele_size in the ORDER BY so the top 500 is a decided
        # set rather than an arbitrary one: hundreds of rows share an allele size, and DuckDB
        # reads the table with several threads, so ordering by allele_size alone would pick a
        # different 500 of the tied rows on each request.
        for motif_category in categories:
            query = ("""SELECT rowid, allele_size, motif_category, affected_status,
                       LocusId, sample_id, Motif, gene_region, GeneTableGeneSymbol,
                       family_id, sex, analysis_status, phenotype_description
                FROM swim_plot
                WHERE outlier_type = ? AND motif_category = ?""" + extra_where + """
                ORDER BY allele_size DESC, LocusId, sample_id LIMIT 500""")
            for row in conn.execute(query, [ot, motif_category] + row_params).fetchall():
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


# The fewest trios a sequencing center needs for the QC2 page to chart its violation fractions
MIN_TRIOS_PER_CENTER = 5


def compute_violation_fractions_by_sample(rows, column_triples):
    """Return each sample's Mendelian-violation fraction in each bin, for the QC2 page.

    Args:
        rows: mendelian_violations rows.
        column_triples: (violations column, total column, bin name) of each bin.

    Returns:
        List of {"sample_id", "bin", "fraction", "violations", "total"} dicts, by sample then bin,
        with a missing count read as 0.
    """
    return [{"sample_id": row["sample_id"], "bin": bin_name,
             "fraction": compute_violation_fraction(row[violations_column], row[total_column]),
             "violations": row[violations_column] or 0, "total": row[total_column] or 0}
            for row in rows for violations_column, total_column, bin_name in column_triples]


def compute_violation_fractions_by_center(rows, column_triples):
    """Return each sequencing center's Mendelian-violation fraction in each bin, for the QC2 page.

    A sample's center is the first three characters of its id, without underscores. Only centers
    with at least MIN_TRIOS_PER_CENTER trios (rows) are reported.

    Args:
        rows: mendelian_violations rows.
        column_triples: (violations column, total column, bin name) of each bin.

    Returns:
        List of {"center", "bin", "fraction", "violations", "total", "trios"} dicts, where
        violations and total are summed over the center's trios, with a missing count read as 0.
    """
    centers = [row["sample_id"][:3].replace("_", "") for row in rows]
    trios_by_center = collections.Counter(centers)
    sums = {}
    for row, center in zip(rows, centers):
        for violations_column, total_column, bin_name in column_triples:
            counts = sums.setdefault((center, bin_name), {"violations": 0, "total": 0})
            counts["violations"] += row[violations_column] or 0
            counts["total"] += row[total_column] or 0
    return [{"center": center, "bin": bin_name,
             "fraction": compute_violation_fraction(counts["violations"], counts["total"]),
             "violations": counts["violations"], "total": counts["total"], "trios": trios_by_center[center]}
            for (center, bin_name), counts in sums.items() if trios_by_center[center] >= MIN_TRIOS_PER_CENTER]


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

    # (violations column, total column, bin name) of each bin
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

    result = {
        "by_motif": compute_violation_fractions_by_sample(rows, motif_columns),
        "by_chrom": compute_violation_fractions_by_sample(rows, chrom_columns),
        "by_center_motif": compute_violation_fractions_by_center(rows, motif_columns),
        "by_center_chrom": compute_violation_fractions_by_center(rows, chrom_columns),
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
    return jsonify({
        "filters": {
            "outlier_type": {
                "type": "enum",
                "required": True,
                "values": list(OUTLIER_TYPE_MAP),
                "labels": {"all": "Long Allele", "short": "Biallelic", "hemi": "Hemizygous"},
            },
            "gene_regions": {
                "type": "multi_enum",
                "values": list(GENE_REGION_MAP),
                "labels": {"cds": "CDS", "promoter": "Promoter", "utr": "UTR (5' & 3')", "intron": "Intron", "exon": "Exon", "intergenic": "Intergenic"},
            },
            "exclude_gene_regions": {
                "type": "multi_enum",
                "values": list(GENE_REGION_MAP),
            },
            "require_above_unaffected": {"type": "enum", "values": list(AFFECTED_ALLELE_COLUMN_PREFIXES)},
            "require_families_above_unaffected": {"type": "enum", "values": list(AFFECTED_ALLELE_COLUMN_PREFIXES)},
            "require_above_population": {
                "type": "enum",
                "values": list(POPULATION_DATASETS_BY_REQUIRE_ABOVE_POPULATION),
                "labels": {"long-read": "Long-read (HPRC256, AoU1027)", "short-read": "Short-read (TenK10K)",
                           "all": "All (long-read and TenK10K)"},
            },
            "include_loci_without_population_data": {
                "type": "boolean",
                "description": "When require_above_population is set, also include loci without population data",
            },
            "min_expansion": {"type": "int", "min": 0},
            "min_repeats_threshold": {"type": "int", "min": 0},
            "min_repeats_above_population": {"type": "int", "description": "Keep loci whose affected allele is at least this many repeats above the largest population value (the page's Δpop column). Uses the cohorts and population_metric of require_above_population, or every cohort when it is not set. Loci without population data are dropped"},
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
            # The UI sends 'search' instead of the parameters below, which are for direct API callers.
            "gene_id": {"type": "string"},
            "gene_symbol": {"type": "string"},
            "locus_id": {"type": "string"},
            "tag": {"type": "string", "description": "Filter to loci with this tag"},
            "affected_unsolved_outlier_sample_id": {"type": "string", "description": "Keep the loci where #Aff counts this sample (affected and unsolved, allele above every unaffected or solved sample's and above the population 99th percentile). Each result then carries SelectedSampleOrFamilyAlleles"},
            "affected_unsolved_outlier_family_id": {"type": "string", "description": "Keep the loci where #Aff counts at least one sample of this family. Each result then carries SelectedSampleOrFamilyAlleles"},
            "motif_size": {"type": "string", "description": "Motif size filter. Comma-separated entries: exact (3), range (3-6), open-end (5-), open-start (-10)"},
        },
        # The extra parameters of /api/v1/outlier_samples (the "Group results by" view), which also
        # takes every filter above, and of its /api/v1/outlier_samples/export download
        "grouped_view_parameters": {
            "group_by": {"type": "enum", "required": True, "values": ["sample", "family"], "description": "One row per affected unsolved sample, or per family"},
            "group_sort": {"type": "enum", "values": list(GROUPED_SORTS), "description": "The column the rows are sorted by (default: num_loci)"},
            "group_sort_direction": {"type": "enum", "values": ["asc", "desc"], "description": "The sort direction (default: desc, or asc for group_id). Rows with no value sort last either way"},
            "format": {"type": "enum", "values": list(GROUPED_EXPORT_FORMATS), "description": "Export only, required: one row per sample or family (groups_tsv), or one per sample and locus (group_loci_tsv)"},
        },
        "sort_options": list(VALID_SORTS),
        # /api/v1/outlier_samples (the "Group results by" view) needs this table
        "has_affected_unsolved_outliers_table": has_affected_unsolved_outliers_table(),
        "group_sort_options": list(GROUPED_SORTS),
        "total_loci": app.config.get("TOTAL_LOCI", 0),
        # None for a database whose metadata table does not record the sample count
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


@app.errorhandler(Exception)
def handle_exception(e):
    # Every uncaught exception lands here. A werkzeug HTTPException (403, 405, 413, ...) keeps
    # its own status rather than being flattened into a 500; 400 and 404 are answered as JSON by
    # the handlers above. Anything else is a crash: print its traceback and answer 500 as JSON.
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

    # Sample metadata (optional).
    if sample_table and os.path.exists(sample_table):
        sample_lookup, affected_lookup, analysis_lookup = load_sample_table(sample_table)
    else:
        sample_lookup, affected_lookup, analysis_lookup = {}, {}, {}

    # Known disease loci (optional), read from local files only. A STRchive-loci JSON, when
    # supplied, adds the detail-page fallback. Either catalog is enough on its own:
    # load_known_disease_loci takes filepath=None to load only the STRchive fallback, so a
    # STRchive-only invocation still gets its disease details.
    known_loci_path = known_loci_json if known_loci_json and os.path.exists(known_loci_json) else None
    strchive_path = strchive_loci_json if strchive_loci_json and os.path.exists(strchive_loci_json) else None
    if known_loci_path or strchive_path:
        disease_trees, strchive_trees, locus_lookup = load_known_disease_loci(
            known_loci_path,
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
        # None for a database whose metadata table does not record the sample count.
        total_samples = read_sample_count(conn)
        # Named in the page's #Fam and Group Results By help text, and the population bar of the
        # exports' outlier sample lists (build_sample_details, compute_plot_thresholds)
        population_gate_cohorts = read_population_gate_cohorts(conn)
        # LIMIT 0 rather than duckdb_compat.table_columns: the columns are advertised as an
        # ordered list via /api/v1/schema, and a query's description keeps the table's order
        # where the catalog lookup returns an unordered set. No row is read either way.
        db_columns = [desc[0] for desc in conn.execute("SELECT * FROM loci LIMIT 0").description]
        tables = duckdb_compat.list_tables(conn)
        # The loci-side phenotype_keyword filter is answered from swim_plot so it selects the same
        # loci the swim plot shows (see build_locus_filter_clauses), which needs both the table and that
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
        # Each grouped-view sample's family as the build recorded it, so the grouped view's sample
        # filters and family ids also cover a sample with no sample-table row, or no --sample-table.
        if AFFECTED_UNSOLVED_OUTLIERS_TABLE in tables:
            outlier_table_family_id_by_sample = dict(conn.execute(
                f"SELECT DISTINCT sample_id, family_id FROM {AFFECTED_UNSOLVED_OUTLIERS_TABLE}").fetchall())
        else:
            outlier_table_family_id_by_sample = {}
    finally:
        conn.close()

    db_columns_set = set(db_columns)
    app.config["TOTAL_LOCI"] = total_loci
    app.config["TOTAL_SAMPLES"] = total_samples
    app.config["POPULATION_GATE_COHORTS"] = population_gate_cohorts
    app.config["DB_COLUMNS"] = db_columns
    app.config["DB_COLUMNS_SET"] = db_columns_set
    app.config["DB_TABLES"] = tables
    app.config["HAS_MENDELIAN"] = "mendelian_violations" in tables
    app.config["HAS_SWIM_PLOT_PHENOTYPES"] = has_swim_plot_phenotypes
    app.config["SAMPLE_IDS"] = sample_ids
    app.config["OUTLIER_TABLE_FAMILY_ID_BY_SAMPLE"] = outlier_table_family_id_by_sample

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

    configure_app(
        args.db,
        sample_table=args.sample_table,
        known_loci_json=args.known_loci_json,
        annotations_db=args.annotations_db,
        strchive_loci_json=args.strchive_loci_json,
    )

    total_loci = app.config["TOTAL_LOCI"]
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
    app.config["INITIAL_REQUIRE_ABOVE_POPULATION"] = args.initial_require_above_population
    app.config["THIS_MACHINES_HOST_NAMES"] = find_this_machines_host_names()
    if args.host not in LOOPBACK_HOSTS:
        if args.allow_remote_writes:
            print(f"\nWARNING: bound to {args.host} with --allow-remote-writes: anyone who can "
                  f"reach this port can read the results and edit notes and tags.")
        else:
            print(f"\nNOTE: bound to {args.host}, so other machines can read the results. Notes "
                  f"and tags can only be edited from this machine; pass --allow-remote-writes to "
                  f"lift that (there is no authentication).")
        print(f"Other machines have to address this one by IP address or as one of: "
              f"{', '.join(sorted(app.config['THIS_MACHINES_HOST_NAMES']))}.")

    print(f"\nStarting server on {args.host}:{args.port}")
    print(f"Web UI: http://{args.host}:{args.port}/")
    app.run(host=args.host, port=args.port, debug=args.debug, use_reloader=args.debug, threaded=True)


if __name__ == "__main__":
    main()
