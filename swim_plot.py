"""Build the swim-plot table: one row per outlier allele across all loci.

The swim plot visualizes every outlier sample's allele size against the locus's
unaffected baseline. This module flattens the per-locus loci records (each of
which carries the three ``OutlierSampleIds_*`` strings) into a long table with
exactly one row per outlier entry, per outlier type, enriched with the sample's
metadata (family / sex / affected status / analysis status / phenotype / purity /
methylation) and a fixed set of locus-level passthrough columns.

This is a faithful, standalone port of ``generate_swim_plot_table`` in the
reference ``analyze_results.py``. The reference iterated a pandas DataFrame; here
the input is the in-memory list of locus dicts produced by the earlier build
stages, so there is no pandas dependency.

Two behaviors carried over verbatim from the reference:

  * ``outlier_rank`` is 1-based in the order ``parse_outlier_entries`` yields,
    which is DESCENDING by allele size — so rank 1 is the largest allele.
  * ``motif_category`` bins motif size as ``"{size}bp"`` for sizes <= 24, else
    ``"25+bp"``; a missing motif size becomes ``"Unknown"``.

The reference re-split the raw outlier entry to recover purity/methylation; here
``parse_outlier_entries`` already returns them as the 3rd/4th tuple fields (with a
literal ``"."`` mapped to None), so this module consumes them directly.
"""

from analysis_columns import (
    OUTLIER_TYPES,
    is_above_unaffected,
    is_missing_outlier_value,
    parse_outlier_entries,
)


# Locus-level columns copied verbatim onto every swim-plot row (matching the
# reference's locus_columns list, in order).
SWIM_PLOT_LOCUS_COLUMNS = [
    "LocusId",
    "Motif",
    "CanonicalMotif",
    "MotifSize",
    "gene_region",
    "GeneTableGeneSymbol",
    "IsInMendelianGene",
    "IsKnownMotif",
    "gene_id",
    "pLI",
    "NumRepeatsInReference",
    "HPRC256_MaxAllele",
    "AoU1027_MaxAllele",
    "TenK10K_MaxAllele",
    "HPRC256_99thPercentile",
    "AoU1027_99thPercentile",
    "TenK10K_99thPercentile",
    "HPRC256_StdevPercentile",
    "AoU1027_StdevPercentile",
]


def _is_missing(value):
    """Return True if value is None or a float NaN (NaN is never equal to itself)."""
    return value is None or (isinstance(value, float) and value != value)


def _compute_motif_category(motif_size):
    """Bin a motif size into the swim-plot motif_category label.

    Args:
        motif_size: The locus motif size in bp (int-like), or a missing value.

    Returns:
        ``"Unknown"`` if the motif size is missing, ``"{size}bp"`` for sizes
        <= 24, else ``"25+bp"``.
    """
    if _is_missing(motif_size):
        return "Unknown"
    if int(motif_size) <= 24:
        return f"{int(motif_size)}bp"
    return "25+bp"


def _normalize_affected_status(raw_affected):
    """Title-case-normalize an affected_status for display, or "Unknown".

    Args:
        raw_affected: The sample's raw affected_status, or a missing value.

    Returns:
        "Affected"/"Unaffected" for those exact (case-insensitive) values,
        otherwise the Title-cased value, or "Unknown" when missing/blank.
    """
    if _is_missing(raw_affected):
        return "Unknown"
    normalized = str(raw_affected).strip().lower()
    if normalized == "affected":
        return "Affected"
    if normalized == "unaffected":
        return "Unaffected"
    return str(raw_affected).strip().title() if raw_affected else "Unknown"


def _normalize_analysis_status(raw_analysis):
    """Title-case-normalize an analysis_status for display, or "Unknown".

    Args:
        raw_analysis: The sample's raw analysis_status, or a missing value.

    Returns:
        "Solved"/"Unsolved"/"Unaffected" for those exact (case-insensitive)
        values, otherwise the Title-cased value, or "Unknown" when missing/blank.
    """
    if _is_missing(raw_analysis):
        return "Unknown"
    normalized = str(raw_analysis).strip().lower()
    if normalized == "solved":
        return "Solved"
    if normalized == "unsolved":
        return "Unsolved"
    if normalized == "unaffected":
        return "Unaffected"
    return str(raw_analysis).strip().title() if raw_analysis else "Unknown"


def _parse_stat_value(raw_value):
    """Convert a purity/methylation string to a float, or None.

    ``parse_outlier_entries`` already maps a literal ``"."`` (and absent fields)
    to None, so by the time a value reaches here it is either None or a numeric
    string. A non-numeric string degrades to None rather than raising.

    Args:
        raw_value: A purity or methylation string, or None.

    Returns:
        The float value, or None if missing / non-numeric.
    """
    if raw_value is None or raw_value == ".":
        return None
    try:
        return float(raw_value)
    except (ValueError, TypeError):
        return None


def generate_swim_plot_table(records, sample_lookup, affected_lookup, analysis_lookup):
    """Flatten locus records into one swim-plot row per outlier allele.

    Walks the three outlier types and, for each locus record carrying a non-empty
    ``OutlierSampleIds_{outlier_type}`` string, parses it (descending by allele
    size) and emits one row per entry. Each row carries the outlier rank, the
    motif_category bin, the sample's metadata (family / sex / affected /
    analysis / phenotype / purity / methylation), the locus's
    FirstUnaffectedAlleleSize for that outlier type, the is_above_first_unaffected
    flag, and the SWIM_PLOT_LOCUS_COLUMNS passthroughs.

    Args:
        records: List of locus record dicts (the build's in-memory representation),
            each carrying OutlierSampleIds_* plus the locus-level passthrough
            columns and FirstUnaffectedAlleleSize_*.
        sample_lookup: Mapping of sample_id -> sample row dict (family_id, sex,
            phenotype_description). A sample id absent here yields None metadata.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.

    Returns:
        list: One row-dict per outlier entry, in (outlier_type, locus, rank)
        iteration order.
    """
    rows = []

    for outlier_type in OUTLIER_TYPES:
        column = f"OutlierSampleIds_{outlier_type}"

        for record in records:
            outlier_value = record.get(column)
            if is_missing_outlier_value(outlier_value):
                continue

            motif_category = _compute_motif_category(record.get("MotifSize"))

            first_unaffected = record.get(f"FirstUnaffectedAlleleSize_{outlier_type}")
            first_unaffected = None if _is_missing(first_unaffected) else int(first_unaffected)

            for rank, (allele_size, sample_id, purity, methylation) in enumerate(
                parse_outlier_entries(outlier_value, locus_id=record.get("LocusId"), column=column),
                start=1,
            ):
                sample_row = sample_lookup.get(sample_id)
                phenotype_description = sample_row.get("phenotype_description") if sample_row else None
                if _is_missing(phenotype_description):
                    phenotype_description = None

                # Prefer the raw sample-row affected_status so "possibly affected" is
                # preserved for display; affected_lookup collapses it to "affected"
                # (that collapsed value is only for analysis logic). Fall back to
                # affected_lookup when the sample row or its field is missing.
                raw_affected = sample_row.get("affected_status") if sample_row else None
                if _is_missing(raw_affected):
                    raw_affected = affected_lookup.get(sample_id)

                row = {
                    "outlier_type": outlier_type,
                    "outlier_rank": rank,
                    "motif_category": motif_category,
                    "SourceDb": record.get("Source"),
                    "allele_size": allele_size,
                    "sample_id": sample_id,
                    "family_id": sample_row.get("family_id") if sample_row else None,
                    "affected_status": _normalize_affected_status(raw_affected),
                    "analysis_status": _normalize_analysis_status(analysis_lookup.get(sample_id)),
                    "sex": sample_row.get("sex") if sample_row else None,
                    "phenotype_description": phenotype_description,
                    "purity": _parse_stat_value(purity),
                    "methylation": _parse_stat_value(methylation),
                    "FirstUnaffectedAlleleSize": first_unaffected,
                    "is_above_first_unaffected": int(is_above_unaffected(allele_size, record, outlier_type)),
                }

                for locus_column in SWIM_PLOT_LOCUS_COLUMNS:
                    row[locus_column] = record.get(locus_column)

                rows.append(row)

    return rows
