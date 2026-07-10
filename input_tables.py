"""Readers for the three input tables that drive a TRails build.

TRails ingests up to three tab-separated input files plus two optional
annotation files:

  1. A repeat-copy-numbers matrix: one row per tandem-repeat locus, with a
     ``trid`` (locus id) column, a ``motif`` column, and one column per sample
     holding that sample's genotype as a comma-separated string of allele sizes
     (e.g. ``"12,40"`` for a diploid call, ``"21"`` for a hemizygous call, or an
     empty / ``.`` / ``./.`` cell for a no-call).
  2. A sample-metadata table: one row per sample, keyed by ``sample_id``, with
     optional ``sex``, ``family_id``, ``maternal_id``, ``paternal_id``,
     ``phenotype_description``, ``affected_status`` and ``analysis_status``.
  3. An optional phenotypes table mapping each ``participant_id`` to a set of
     HPO ``term_id`` values.

Plus two optional reference tables:

  4. A gene table (``gene_id`` keyed) supplying ``pLI`` inputs and inheritance.
  5. A gene-to-disease-phenotype table (HPO ``genes_to_phenotype.txt`` format)
     used by phenotype scoring.

Column-name matching across the input tables is case-insensitive AND
underscore/space-insensitive (see ``normalize_column_name``); ``sample_id``
*values* are still joined exactly. The readers here are deliberately permissive:
optional columns absent -> the field is simply omitted / NULL, and the build
continues (with a printed warning) rather than hard-failing.
"""

import gzip

import pandas


# Allowed normalized values for sample-metadata status columns. A value outside
# the corresponding set raises (mirrors validate_sample_metadata in the source).
ALLOWED_AFFECTED_STATUS = {"affected", "unaffected", "unknown", "possibly affected"}
ALLOWED_ANALYSIS_STATUS = {
    "solved",
    "unsolved",
    "unknown",
    "unaffected",
    "probably solved",
    "partially solved",
}
ALLOWED_SEX = {"male", "female", "unknown"}

# analysis_status codes that are remapped to canonical values before validation.
ANALYSIS_STATUS_REMAP = {
    "rncc": "unsolved",
    "rcpc": "unsolved",
    "s_kgfp": "solved",
}


def normalize_column_name(name):
    """Return a canonical key for a column name.

    Lowercases, strips surrounding whitespace, and removes underscores and
    spaces so that e.g. ``"Sample_ID"``, ``"sample id"`` and ``"sampleid"`` all
    collapse to the same key.

    Args:
        name: The raw column name.

    Returns:
        The normalized key string.
    """
    return str(name).strip().lower().replace("_", "").replace(" ", "")


def match_columns(df, logical_names):
    """Map each logical column name to the actual column present in ``df``.

    Matching is performed on normalized keys (see ``normalize_column_name``), so
    case and underscore/space differences are ignored. A logical name with no
    matching column is simply omitted from the returned mapping.

    Args:
        df: A pandas DataFrame whose columns are searched.
        logical_names: An iterable of desired logical column names.

    Returns:
        A dict mapping each matched logical name to the actual column name in
        ``df``. The first actual column matching a given normalized key wins.
    """
    normalized_to_actual = {}
    for actual_column in df.columns:
        normalized_to_actual.setdefault(normalize_column_name(actual_column), actual_column)

    matches = {}
    for logical_name in logical_names:
        actual = normalized_to_actual.get(normalize_column_name(logical_name))
        if actual is not None:
            matches[logical_name] = actual
    return matches


# Matrix columns whose header names a known per-locus annotation rather than a
# sample. When present they are promoted to the locus record (pass-through)
# instead of being treated as a sample genotype column — so they cannot pollute
# the allele histograms/outlier lists and they reach the downstream annotation
# stages. Matched case/underscore-insensitively; the value maps each normalized
# header to the canonical column name the record/output schema expects.
ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED = {normalize_column_name(name): name for name in [
    "gene_id", "GencodeGeneId", "gene_region", "GencodeGeneRegion", "gene_region_rank",
    "ReferenceRegion", "NumRepeatsInReference", "CanonicalMotif", "MotifSize",
    "Chrom", "Start0Based", "End1Based", "KnownDiseaseLocus", "IsKnownMotif",
    "IsInMendelianGene", "Source", "NonCodingAnnotations", "RepeatMaskerIntervals",
    "VariationClusterSizeDiff",
]}
# Per-cohort population-stat / TRExplorer annotation columns share these header
# prefixes (e.g. HPRC256_99thPercentile, AoU1027_Mode, TRExplorerMotif). They are
# passed through verbatim under their own header name. The prefixes are the
# specific cohort/version tokens actually emitted upstream (not the bare cohort
# abbreviation) so they can't swallow a genuine sample column whose id merely
# starts with a cohort name (e.g. a sample named "AoU_0001").
ANNOTATION_COLUMN_PREFIXES = ("hprc256", "aou1027", "tenk10k", "trexplorer")


def _coerce_annotation_value(value):
    """Coerce a raw string annotation cell to int/float when numeric, else str/None.

    Missing cells (empty, or a recognized TSV missing marker such as ``.`` /
    ``NA`` / ``N/A`` / ``NaN`` / ``None`` / ``null``, case-insensitively) become
    None; integer-looking values become int; other numeric values become float;
    everything else is left as the original string. Normalizing missing markers
    matters because the population-stat columns feed numeric comparisons
    downstream (e.g. the p99 thresholds in analysis_columns), where a leftover
    ``.``/``NA`` string would raise a TypeError against the numeric cohort values.
    """
    if value is None or str(value).strip().lower() in ("", ".", "na", "n/a", "nan", "none", "null"):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def _is_annotation_column(column_name):
    """Return the canonical name if column_name is a known annotation column, else None."""
    normalized = normalize_column_name(column_name)
    if normalized in ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED:
        return ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED[normalized]
    if normalized.startswith(ANNOTATION_COLUMN_PREFIXES):
        return column_name  # cohort-stat columns pass through under their own header
    return None


def read_repeat_copy_numbers(path):
    """Read the repeat-copy-numbers matrix into per-locus rows.

    The first two logical columns are ``trid`` (the locus id) and ``motif``,
    matched case/underscore-insensitively. Of the remaining columns, any whose
    header names a known per-locus annotation (see
    ``ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED`` /
    ``ANNOTATION_COLUMN_PREFIXES`` — e.g. ``gene_id``, ``GencodeGeneRegion``,
    ``HPRC256_99thPercentile``) is treated as an annotation pass-through column,
    NOT a sample; everything else is a sample whose value for a locus is that
    sample's genotype cell string. Separating annotation columns out keeps them
    from being mistaken for samples (which would otherwise pollute the allele
    histograms/outlier lists and never reach the annotation stages).

    Gzipped (``.gz``) files are supported.

    Args:
        path: Path to the matrix TSV (optionally gzipped).

    Returns:
        A tuple ``(locus_rows, sample_id_list)`` where ``locus_rows`` is a list
        of dicts, one per locus, each ``{"trid": ..., "motif": ...,
        "genotypes": {sample_id: cell_string}, "extra_columns": {canonical: value}}``
        and ``sample_id_list`` is the ordered list of sample-column names.
    """
    df = pandas.read_table(path, dtype=str, keep_default_na=False)

    matches = match_columns(df, ["trid", "motif"])
    if "trid" not in matches:
        raise ValueError(
            f"Repeat-copy-numbers table {path} is missing a 'trid' (locus id) column; "
            f"found columns: {list(df.columns)[:10]}"
        )
    if "motif" not in matches:
        raise ValueError(
            f"Repeat-copy-numbers table {path} is missing a 'motif' column; "
            f"found columns: {list(df.columns)[:10]}"
        )

    trid_column = matches["trid"]
    motif_column = matches["motif"]

    sample_id_list = []
    annotation_columns = {}  # actual header -> canonical record key
    for column in df.columns:
        if column in (trid_column, motif_column):
            continue
        canonical = _is_annotation_column(column)
        if canonical is None:
            sample_id_list.append(column)
        else:
            annotation_columns[column] = canonical

    locus_rows = []
    for record in df.to_dict(orient="records"):
        locus_rows.append({
            "trid": record[trid_column],
            "motif": record[motif_column],
            "genotypes": {sample_id: record[sample_id] for sample_id in sample_id_list},
            "extra_columns": {canonical: _coerce_annotation_value(record[actual])
                              for actual, canonical in annotation_columns.items()},
        })

    return locus_rows, sample_id_list


def _is_blank(value):
    """Return True if a metadata cell should be treated as missing."""
    if value is None:
        return True
    if isinstance(value, float):
        return pandas.isna(value)
    return str(value).strip().lower() in {"", "nan", "none", "na", "n/a", "null"}


def read_sample_metadata(path):
    """Read the sample-metadata table and build lookup dictionaries.

    Only ``sample_id`` is required (the original required eight columns; this is
    relaxed to an additive guard so a build can run with minimal metadata). All
    other recognized columns (``sex``, ``family_id``, ``maternal_id``,
    ``paternal_id``, ``phenotype_description``, ``affected_status``,
    ``analysis_status``) are matched case/underscore-insensitively and included
    when present.

    Normalizations applied (mirroring the source ``load_sample_table`` +
    ``validate_sample_metadata``):
      * ``analysis_status``: lowercased+stripped, remapped via
        ``ANALYSIS_STATUS_REMAP`` (``rncc``/``rcpc`` -> ``unsolved``,
        ``s_kgfp`` -> ``solved``), then blanks -> ``unknown``.
      * ``affected_status`` (for the affected lookup): lowercased+stripped, with
        ``possibly affected`` -> ``affected``.
      * ``phenotype_description``: a leading ``"NA; "`` prefix is stripped.
      * ``sample_id`` must not contain ``':'`` (it would break OutlierSampleIds
        parsing) -> raises ValueError.
      * Duplicate ``sample_id`` rows are dropped, keeping the first.
      * Values outside the allowed sets for affected/analysis/sex raise.

    Args:
        path: Path to the sample-metadata TSV (optionally gzipped).

    Returns:
        A tuple ``(sample_lookup, affected_lookup, analysis_lookup, sample_df)``
        where ``sample_lookup`` maps each sample_id to a dict of its metadata
        fields, ``affected_lookup`` maps sample_id -> normalized affected_status
        (or None), ``analysis_lookup`` maps sample_id -> normalized
        analysis_status, and ``sample_df`` is the cleaned DataFrame (with
        canonical-named columns: ``sample_id`` plus whichever optional columns
        were present).
    """
    raw_df = pandas.read_table(path, dtype=str, keep_default_na=False)

    matches = match_columns(raw_df, [
        "sample_id", "sex", "family_id", "maternal_id", "paternal_id",
        "phenotype_description", "affected_status", "analysis_status",
    ])
    if "sample_id" not in matches:
        raise ValueError(
            f"Sample-metadata table {path} is missing a required 'sample_id' column; "
            f"found columns: {list(raw_df.columns)[:10]}"
        )

    # Rebuild a DataFrame with canonical column names, only the columns present.
    df = pandas.DataFrame()
    for logical_name, actual_column in matches.items():
        df[logical_name] = raw_df[actual_column]

    # sample_id must not contain ':'.
    bad_sample_ids = [s for s in df["sample_id"] if ":" in str(s)]
    if bad_sample_ids:
        raise ValueError(
            f"sample_id values must not contain ':' (breaks OutlierSampleIds parsing): "
            f"{bad_sample_ids[:10]}"
        )

    # Normalize analysis_status.
    if "analysis_status" in df.columns:
        df["analysis_status"] = df["analysis_status"].apply(
            lambda v: ANALYSIS_STATUS_REMAP.get(str(v).strip().lower(), str(v).strip().lower())
        )
        df["analysis_status"] = df["analysis_status"].apply(
            lambda v: "unknown" if v in {"", "nan", "none", "na", "n/a", "null"} else v
        )

    # Strip leading "NA; " from phenotype_description.
    if "phenotype_description" in df.columns:
        df["phenotype_description"] = df["phenotype_description"].apply(
            lambda p: p[len("NA; "):] if isinstance(p, str) and p.startswith("NA; ") else p
        )

    # Drop duplicate sample_ids, keeping the first.
    df = df.drop_duplicates(subset=["sample_id"], keep="first").reset_index(drop=True)

    # Validate allowed value sets (raise if any value falls outside).
    _validate_status_column(df, "affected_status", ALLOWED_AFFECTED_STATUS)
    _validate_status_column(df, "analysis_status", ALLOWED_ANALYSIS_STATUS)
    _validate_status_column(df, "sex", ALLOWED_SEX)

    sample_lookup = {}
    affected_lookup = {}
    analysis_lookup = {}
    for record in df.to_dict(orient="records"):
        sample_id = record["sample_id"]
        sample_lookup[sample_id] = dict(record)
        affected_lookup[sample_id] = _normalize_affected_status_for_logic(record.get("affected_status"))
        analysis_lookup[sample_id] = record.get("analysis_status")

    return sample_lookup, affected_lookup, analysis_lookup, df


def _validate_status_column(df, column, allowed_set):
    """Raise ValueError if any non-blank value in ``column`` is outside ``allowed_set``."""
    if column not in df.columns:
        return
    errors = []
    for sample_id, value in zip(df["sample_id"], df[column]):
        if _is_blank(value):
            continue
        if str(value).strip().lower() not in allowed_set:
            errors.append(f"Unexpected {column} '{value}' for sample_id {sample_id}")
            if len(errors) >= 20:
                break
    if errors:
        raise ValueError(
            f"Found unexpected {column} values (showing up to 20):\n" + "\n".join(errors)
        )


def _normalize_affected_status_for_logic(value):
    """Lowercase + collapse 'possibly affected' to 'affected'; None for blanks."""
    if _is_blank(value):
        return None
    normalized = str(value).strip().lower()
    return "affected" if normalized == "possibly affected" else normalized


def read_phenotypes(path):
    """Read the phenotypes table into a participant -> HPO-term-set mapping.

    Expects ``participant_id`` and ``term_id`` columns (the format produced by
    generate_phenotype_table.py), matched case/underscore-insensitively. Only
    terms starting with ``"HP:"`` are retained.

    Args:
        path: Path to the phenotypes TSV (optionally gzipped), or None.

    Returns:
        A dict mapping each participant_id to a set of HPO term ids. Returns an
        empty dict if ``path`` is None.
    """
    if path is None:
        return {}

    participant_to_hpo = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as input_file:
        header = input_file.readline().rstrip("\n").split("\t")
        normalized_header = [normalize_column_name(c) for c in header]
        if "participantid" not in normalized_header or "termid" not in normalized_header:
            raise ValueError(
                f"Phenotypes table {path} must have 'participant_id' and 'term_id' columns; "
                f"found: {header}"
            )
        participant_index = normalized_header.index("participantid")
        term_index = normalized_header.index("termid")
        for line in input_file:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(participant_index, term_index):
                continue
            if fields[term_index].startswith("HP:"):
                participant_to_hpo.setdefault(fields[participant_index], set()).add(fields[term_index])

    return participant_to_hpo


def read_gene_table(path):
    """Read the gene table into a gene_id -> annotations mapping.

    Keeps the columns used downstream when present: ``gene_symbol``,
    ``gene_aliases``, ``pLI_v2``, ``pLI_v4``, ``lof_oe_ci_upper_v4``,
    ``hgnc_gene_id``, ``inheritance``, ``disease_category``,
    ``LLM_phenotype_summary``, ``sources``.

    Args:
        path: Path to the gene-table TSV (optionally gzipped), or None.

    Returns:
        A dict mapping each gene_id to a dict of its retained annotation fields.
        Returns an empty dict if ``path`` is None.
    """
    if path is None:
        return {}

    df = pandas.read_table(path)
    if "gene_id" not in df.columns:
        raise ValueError(
            f"Gene table {path} is missing a required 'gene_id' column; "
            f"found columns: {list(df.columns)[:10]}"
        )

    columns_to_keep = [c for c in [
        "gene_id", "gene_symbol", "gene_aliases", "pLI_v2", "pLI_v4",
        "lof_oe_ci_upper_v4", "hgnc_gene_id", "inheritance", "disease_category",
        "LLM_phenotype_summary", "sources",
    ] if c in df.columns]

    return df[columns_to_keep].set_index("gene_id").to_dict(orient="index")


# HPO term ids that encode an inheritance mode rather than a phenotype, and
# their mapping to short inheritance codes. The single source of truth for this
# classification (read_gene_disease_phenotypes uses it to route terms into a
# disease's inheritance vs. phenotypes set). Terms present here but absent from
# INHERITANCE_MAP (the non-Mendelian / polygenic group) carry no actionable
# inheritance code and are simply dropped from the inheritance set.
INHERITANCE_HPO_TERMS = {
    "HP:0000006",  # Autosomal dominant
    "HP:0000007",  # Autosomal recessive
    "HP:0001417",  # X-linked
    "HP:0001419",  # X-linked recessive
    "HP:0001423",  # X-linked dominant
    "HP:0001426",  # Non-Mendelian
    "HP:0001427",  # Mitochondrial
    "HP:0001450",  # Y-linked
    "HP:0010982",  # Polygenic
    "HP:0010983",  # Oligogenic
    "HP:0010984",  # Digenic
    "HP:0012275",  # Autosomal dominant with maternal imprinting
    "HP:0034341",  # Pseudoautosomal recessive
}
INHERITANCE_MAP = {
    "HP:0000006": "AD",
    "HP:0000007": "AR",
    "HP:0001417": "XL",
    "HP:0001419": "XR",
    "HP:0001423": "XD",
    "HP:0001427": "MT",
    "HP:0001450": "YL",
    "HP:0012275": "AD",  # AD with maternal imprinting -> autosomal dominant
    "HP:0034341": "AR",  # Pseudoautosomal recessive -> autosomal recessive
}


def read_gene_disease_phenotypes(path):
    """Read the HPO ``genes_to_phenotype`` table for phenotype scoring.

    The file is the HPO ``genes_to_phenotype.txt`` format: a header line
    followed by tab-separated rows whose 2nd field is the gene symbol, 3rd field
    is an HPO id, and 6th field is a disease id. Inheritance-mode HPO terms (see
    ``INHERITANCE_HPO_TERMS``) populate each disease's ``inheritance`` set;
    every other HPO id populates its ``phenotypes`` set.

    Args:
        path: Path to the genes_to_phenotype table (optionally gzipped), or None.

    Returns:
        A dict ``{gene_symbol: {disease_id: {"phenotypes": set,
        "inheritance": set}}}``. Returns an empty dict if ``path`` is None.
    """
    if path is None:
        return {}

    gene_disease_data = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as input_file:
        input_file.readline()  # skip header
        for line in input_file:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            disease_entry = gene_disease_data.setdefault(fields[1], {}).setdefault(
                fields[5], {"phenotypes": set(), "inheritance": set()}
            )
            if fields[2] in INHERITANCE_HPO_TERMS:
                # Only Mendelian inheritance modes carry an actionable short code;
                # non-Mendelian/polygenic terms (absent from INHERITANCE_MAP) are
                # dropped rather than added raw, so they can't make phenotype
                # scoring treat the disease as inheritance-restricted.
                if fields[2] in INHERITANCE_MAP:
                    disease_entry["inheritance"].add(INHERITANCE_MAP[fields[2]])
            else:
                disease_entry["phenotypes"].add(fields[2])

    return gene_disease_data
