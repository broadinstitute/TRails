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
continues rather than hard-failing.
"""

import gzip

import pandas

import analysis_columns


# Allowed normalized values for sample-metadata status columns. A value outside
# the corresponding set raises.
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
# sample. When present they are promoted to the locus record
# instead of being treated as a sample genotype column — so they cannot pollute
# the allele histograms/outlier lists and they reach the downstream annotation
# stages. Being recognized is not the same as being kept: Chrom, Start0Based,
# End1Based, MotifSize, CanonicalMotif and gene_region_rank are always recomputed
# from the LocusId, Motif and gene_region
# (locus_annotations.SOURCE_COLUMN_BY_RECOMPUTED_COLUMN), so for those the header
# only stays out of the sample list. KnownDiseaseLocus, IsKnownMotif and
# IsInMendelianGene are kept only when the build has no known-loci catalog or gene
# table to compute them from. The rest are kept as supplied. Matched
# case/underscore-insensitively; the value maps each normalized header to the
# canonical column name the record/output schema expects.
ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED = {normalize_column_name(name): name for name in [
    "gene_id", "GencodeGeneId", "gene_region", "GencodeGeneRegion", "gene_region_rank",
    "ReferenceRegion", "NumRepeatsInReference", "CanonicalMotif", "MotifSize",
    "Chrom", "Start0Based", "End1Based", "KnownDiseaseLocus", "IsKnownMotif",
    "IsInMendelianGene", "Source", "NonCodingAnnotations", "RepeatMaskerIntervals",
    "VariationClusterSizeDiff", "VariationCluster", "VariationClusterFilterReason",
]}
# Per-cohort population-stat / TRExplorer annotation columns share these header
# prefixes (e.g. HPRC256_99thPercentile, AoU1027_Mode, TenK10K_MaxAllele,
# TRExplorerMotif). A header matching a loci-table output column is stored under
# that column's canonical spelling (see _NORMALIZED_OUTPUT_COLUMNS below); any other
# prefix-matched header is passed through under its own name.
# Any header whose normalized form starts with one of these tokens is read as an
# annotation, never as a sample. "hprc256" and "aou1027" carry the cohort size, so
# a sample whose id merely starts with the cohort name (e.g. "AoU_0001") is still
# a sample. "tenk10k" and "trexplorer" are the bare cohort and tool names, so a
# sample column whose id starts with them (e.g. "TenK10K_S1") is also read as an
# annotation.
ANNOTATION_COLUMN_PREFIXES = ("hprc256", "aou1027", "tenk10k", "trexplorer")

# Normalized -> canonical spelling for every loci-table output column, so a
# prefix-matched cohort-stat header supplied in non-canonical casing/underscoring
# (e.g. "hprc256_99thpercentile") is stored under the exact OUTPUT_COLUMNS name the
# server reads ("HPRC256_99thPercentile") rather than as a separate raw column.
_NORMALIZED_OUTPUT_COLUMNS = {
    normalize_column_name(name): name for name in analysis_columns.OUTPUT_COLUMNS
}


def _coerce_annotation_value(value):
    """Coerce a raw string annotation cell to int/float when numeric, else str/None.

    Missing cells (empty, or a missing marker in analysis_columns.MISSING_VALUE_SPELLINGS
    such as ``.`` / ``NA`` / ``N/A`` / ``NaN`` / ``None`` / ``null``,
    case-insensitively) become None; integer-looking values become int; other
    numeric values become float; everything else is left as the original string.
    Normalizing missing markers matters because the population-stat columns feed
    numeric comparisons downstream (e.g. the p99 thresholds in analysis_columns),
    where a leftover ``.``/``NA`` string would raise a TypeError against the numeric
    cohort values.
    """
    if analysis_columns.blank_to_none(value) is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def _get_canonical_annotation_column_name_or_none(column_name):
    """Return the record key a matrix column's values are stored under, or None for a sample column.

    Args:
        column_name: A matrix column header.

    Returns:
        The canonical annotation column name for a recognized annotation header (the
        header itself for a prefix-matched one outside the loci-table schema), or None
        when the column is a sample.
    """
    normalized = normalize_column_name(column_name)
    if normalized in ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED:
        return ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED[normalized]
    if normalized.startswith(ANNOTATION_COLUMN_PREFIXES):
        # Canonicalize known cohort-stat columns to their OUTPUT_COLUMNS spelling;
        # genuinely unknown extras (not in the schema) pass through under their header.
        return _NORMALIZED_OUTPUT_COLUMNS.get(normalized, column_name)
    return None


def read_repeat_copy_numbers(path):
    """Read the repeat-copy-numbers matrix into per-locus rows.

    The first two logical columns are ``trid`` (the locus id) and ``motif``,
    matched case/underscore-insensitively. Of the remaining columns, any whose
    header names a known per-locus annotation (see
    ``ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED`` /
    ``ANNOTATION_COLUMN_PREFIXES`` — e.g. ``gene_id``, ``GencodeGeneRegion``,
    ``HPRC256_99thPercentile``) is treated as an annotation column,
    NOT a sample; everything else is a sample whose value for a locus is that
    sample's genotype cell string. Separating annotation columns out keeps them
    from being mistaken for samples (which would otherwise pollute the allele
    histograms/outlier lists and never reach the annotation stages). The returned
    values are the matrix's own (numeric ones coerced); the build later recomputes
    some of these columns (see the comment on
    ``ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED``).

    Gzipped (``.gz``) files are supported.

    Args:
        path: Path to the matrix TSV (optionally gzipped).

    Returns:
        A tuple ``(locus_rows, sample_id_list)`` where ``locus_rows`` is a list
        of dicts, one per locus, each ``{"trid": ..., "motif": ...,
        "genotypes": {sample_id: cell_string}, "extra_columns": {canonical: value}}``
        and ``sample_id_list`` is the ordered list of sample-column names.

    Raises:
        ValueError: If the ``trid`` or ``motif`` column is missing, or a sample
            column header contains ``':'`` or ``','``.
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
        canonical = _get_canonical_annotation_column_name_or_none(column)
        if canonical is None:
            sample_id_list.append(column)
        else:
            annotation_columns[column] = canonical
    _reject_sample_ids_containing_outlier_delimiters(sample_id_list, f"Sample column headers in {path}")

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


def _reject_sample_ids_containing_outlier_delimiters(sample_ids, description):
    """Raise a ValueError naming the sample ids that contain ':' or ','.

    The OutlierSampleIds_* strings pack entries as "{allele}x:{sample_id}" joined with ",",
    so either character in a sample id corrupts the parse and can silently merge two
    samples into one. The matrix's sample columns and the sample table's sample_id values
    both end up in those strings, so both readers apply this one rule.

    Args:
        sample_ids: The sample ids to check.
        description: What the ids are, named in the error, e.g. "Sample column headers in
            matrix.tsv".

    Raises:
        ValueError: If any id contains ':' or ','. The message lists up to 10 of them.
    """
    bad_sample_ids = [sample_id for sample_id in sample_ids if ":" in str(sample_id) or "," in str(sample_id)]
    if bad_sample_ids:
        raise ValueError(
            f"{description} must not contain ':' or ',' (both are OutlierSampleIds "
            f"delimiters): {bad_sample_ids[:10]}")


def read_sample_metadata(path):
    """Read the sample-metadata table and build lookup dictionaries.

    Only ``sample_id`` is required, so a build can run with minimal metadata. All
    other recognized columns (``sex``, ``family_id``, ``maternal_id``,
    ``paternal_id``, ``phenotype_description``, ``affected_status``,
    ``analysis_status``) are matched case/underscore-insensitively and included
    when present.

    Normalizations applied:
      * ``analysis_status``: lowercased+stripped, remapped via
        ``ANALYSIS_STATUS_REMAP`` (``rncc``/``rcpc`` -> ``unsolved``,
        ``s_kgfp`` -> ``solved``), then blanks -> ``unknown``.
      * ``affected_status`` (for the affected lookup): lowercased+stripped, with
        ``possibly affected`` -> ``affected``.
      * ``phenotype_description``: a leading ``"NA; "`` prefix is stripped.
      * ``sample_id`` must not contain ``':'`` or ``','`` (both are
        OutlierSampleIds delimiters) -> raises ValueError.
      * Duplicate ``sample_id`` rows are dropped, keeping the first.
      * Values outside the allowed sets for affected/analysis/sex raise.

    Args:
        path: Path to the sample-metadata TSV (optionally gzipped).

    Returns:
        A tuple ``(sample_lookup, affected_lookup, analysis_lookup)`` where
        ``sample_lookup`` maps each sample_id, in file order, to a dict of its
        metadata fields (``sample_id`` plus whichever optional columns were
        present, under their canonical names), ``affected_lookup`` maps
        sample_id -> normalized affected_status (or None), and
        ``analysis_lookup`` maps sample_id -> normalized analysis_status.
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

    _reject_sample_ids_containing_outlier_delimiters(df["sample_id"], f"sample_id values in {path}")

    # Normalize analysis_status.
    if "analysis_status" in df.columns:
        df["analysis_status"] = df["analysis_status"].apply(
            lambda v: ANALYSIS_STATUS_REMAP.get(str(v).strip().lower(), str(v).strip().lower())
        )
        df["analysis_status"] = df["analysis_status"].apply(
            lambda v: "unknown" if analysis_columns.blank_to_none(v) is None else v
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
        affected_lookup[sample_id] = analysis_columns.normalize_affected_status_for_logic(
            record.get("affected_status"))
        analysis_lookup[sample_id] = record.get("analysis_status")

    return sample_lookup, affected_lookup, analysis_lookup


def _validate_status_column(df, column, allowed_set):
    """Raise ValueError if any non-blank value in ``column`` is outside ``allowed_set``."""
    if column not in df.columns:
        return
    errors = []
    for sample_id, value in zip(df["sample_id"], df[column]):
        if analysis_columns.blank_to_none(value) is None:
            continue
        if str(value).strip().lower() not in allowed_set:
            errors.append(f"Unexpected {column} '{value}' for sample_id {sample_id}")
            if len(errors) >= 20:
                break
    if errors:
        raise ValueError(
            f"Found unexpected {column} values (showing up to 20):\n" + "\n".join(errors)
        )


def read_phenotypes(path):
    """Read the phenotypes table into a participant -> HPO-term-set mapping.

    Expects ``participant_id`` and ``term_id`` columns (the format described in
    docs/INPUT_FORMATS.md), matched case/underscore-insensitively. Only terms
    starting with ``"HP:"`` are retained.

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

    Keeps, when present, the fields analysis_columns.GENE_TABLE_COLUMN_BY_SOURCE_FIELD
    lists (gene_symbol, the pLI and LOEUF values, inheritance, sources, ...), the
    ones the downstream stages read.

    Args:
        path: Path to the gene-table TSV (optionally gzipped), or None.

    Returns:
        A dict mapping each gene_id to a dict of its retained annotation fields.
        Returns an empty dict if ``path`` is None.

    Raises:
        ValueError: If the table has no ``gene_id`` column, or lists a gene_id on more
            than one row.
    """
    if path is None:
        return {}

    # Missing markers: pandas' defaults ("", "NA", "N/A", "null", "nan", ...) plus TRails' own
    # ".", which the other inputs honor too (analysis_columns.MISSING_VALUE_SPELLINGS). Without it a "."
    # in pLI_v2/pLI_v4 arrives as the literal string and analysis_columns.add_gene_columns dies on
    # float(".").
    df = pandas.read_table(path, na_values=["."])

    # Matched with match_columns rather than by exact name, so this reader honors the
    # case/underscore-insensitive rule the module docstring states for every input table: a header
    # spelled "Gene_ID", "PLI_V4" or "gene symbol" would otherwise be silently dropped (or, for the
    # required gene_id column, rejected outright).
    matches = match_columns(df, ["gene_id"] + list(analysis_columns.GENE_TABLE_COLUMN_BY_SOURCE_FIELD))
    if "gene_id" not in matches:
        raise ValueError(
            f"Gene table {path} is missing a required 'gene_id' column; "
            f"found columns: {list(df.columns)[:10]}"
        )
    df = df.rename(columns={actual: logical for logical, actual in matches.items()})
    columns_to_keep = list(matches)

    if "hgnc_gene_id" in columns_to_keep:
        # A single blank cell upcasts this numeric-looking id column to float64, which would write
        # "1100.0" into GeneTableHgncGeneId. Render it back as an id string; blanks stay NaN so
        # analysis_columns.add_gene_columns (is_none_or_nan) still normalizes them to None.
        df["hgnc_gene_id"] = df["hgnc_gene_id"].apply(
            lambda value: value if analysis_columns.blank_to_none(value) is None
            else str(int(value)) if isinstance(value, float) else str(value))

    # One row per gene, as genelookup's table has. A table with one row per gene-disease pair can
    # disagree between a gene's rows (AD on one, AR on the next), so its rows are neither merged
    # nor dropped; without this check to_dict fails with a pandas error naming no file or gene.
    repeated_gene_ids = df.loc[df["gene_id"].duplicated(), "gene_id"].unique()
    if len(repeated_gene_ids):
        raise ValueError(
            f"Gene table {path} lists {len(repeated_gene_ids)} gene_id value(s) on more than one "
            f"row: {', '.join(str(gene_id) for gene_id in repeated_gene_ids[:10])}"
            f"{', ...' if len(repeated_gene_ids) > 10 else ''}. It needs one row per gene, so "
            f"combine each gene's rows into one")

    return df[columns_to_keep].set_index("gene_id").to_dict(orient="index")


# HPO term ids that encode an inheritance mode rather than a phenotype (descendants of
# HP:0000005 Mode of inheritance, including every one that genes_to_phenotype.txt v2025-05-06
# uses), and their mapping to short inheritance codes. The single source of truth for this
# classification (read_gene_disease_phenotypes uses it to route terms into a disease's
# inheritance vs. phenotypes set). Terms present here but absent from INHERITANCE_MAP (the
# non-Mendelian / polygenic group and qualifiers such as incomplete penetrance or genetic
# anticipation) carry no actionable inheritance code, so they are dropped: they enter neither
# the inheritance set nor the phenotype set that phenotype scoring compares patients against.
INHERITANCE_HPO_TERMS = {
    "HP:0000006",  # Autosomal dominant
    "HP:0000007",  # Autosomal recessive
    "HP:0001417",  # X-linked
    "HP:0001419",  # X-linked recessive
    "HP:0001423",  # X-linked dominant
    "HP:0001426",  # Non-Mendelian
    "HP:0001427",  # Mitochondrial
    "HP:0001442",  # Typified by somatic mosaicism
    "HP:0001450",  # Y-linked
    "HP:0001466",  # Contiguous gene syndrome
    "HP:0001470",  # Sex-limited expression
    "HP:0003743",  # Genetic anticipation
    "HP:0003744",  # Genetic anticipation with paternal anticipation bias
    "HP:0003745",  # Sporadic
    "HP:0003829",  # Typified by incomplete penetrance
    "HP:0003831",  # Typified by age-related disease onset
    "HP:0010982",  # Polygenic
    "HP:0010983",  # Oligogenic
    "HP:0010984",  # Digenic
    "HP:0012275",  # Autosomal dominant with maternal imprinting
    "HP:0025352",  # Typically de novo
    "HP:0032382",  # Uniparental disomy
    "HP:0034341",  # Pseudoautosomal recessive
    "HP:0034344",  # Female-limited expression
    "HP:4000158",  # Typified by high penetrance
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

    Gene symbols are upper-cased, because the gene table and HPO can spell the same
    symbol differently: HPO writes the HGNC symbol "C9orf72", while the default gene
    table (genelookup's) upper-cases it to "C9ORF72". Callers look a symbol up by its
    upper-cased form.

    Args:
        path: Path to the genes_to_phenotype table (optionally gzipped), or None.

    Returns:
        A dict ``{upper_cased_gene_symbol: {disease_id: {"phenotypes": set,
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
            disease_entry = gene_disease_data.setdefault(fields[1].upper(), {}).setdefault(
                fields[5], {"phenotypes": set(), "inheritance": set()}
            )
            if fields[2] in INHERITANCE_HPO_TERMS:
                # Only Mendelian inheritance modes carry an actionable short code;
                # the other inheritance terms (absent from INHERITANCE_MAP) are
                # dropped rather than added raw, so they can't make phenotype
                # scoring treat the disease as inheritance-restricted.
                if fields[2] in INHERITANCE_MAP:
                    disease_entry["inheritance"].add(INHERITANCE_MAP[fields[2]])
            else:
                disease_entry["phenotypes"].add(fields[2])

    return gene_disease_data
