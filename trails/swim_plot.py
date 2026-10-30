"""Build the swim-plot table: one row per outlier allele across all loci.

The swim plot visualizes every outlier sample's allele size against the locus's
unaffected baseline. This module flattens the per-locus loci records (each of
which carries the three ``OutlierSampleIds_*`` strings) into a long table with
exactly one row per outlier entry, per outlier type, enriched with the sample's
metadata (family / sex / affected status / analysis status / phenotype) and a
fixed set of locus-level passthrough columns.

The input is the in-memory list of locus dicts produced by the earlier build
stages, so there is no pandas dependency. Per row:

  * ``outlier_rank`` is 1-based in the order ``parse_outlier_entries`` yields,
    which is DESCENDING by allele size — so rank 1 is the largest allele.
  * ``motif_category`` bins motif size as ``"{size}bp"`` for sizes <= 24, else
    ``"25+bp"``; a missing motif size becomes ``"Unknown"``.
  * A blank family_id, sex or phenotype_description is stored as NULL, as in the
    other tables (``analysis_columns.blank_to_none``), and a blank affected or
    analysis status is shown as ``"Unknown"``.
"""

from analysis_columns import (
    OUTLIER_TYPES,
    blank_to_none,
    get_family_id,
    is_above_unaffected,
    is_missing_outlier_value,
    is_none_or_nan,
    parse_outlier_entries,
)


# Locus-level columns copied verbatim onto every swim-plot row, in order.
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

# The swim_plot table's columns, in order, matching the row dicts generate_swim_plot_table
# emits (the per-row keys below, then the copied locus columns).
# result_database.write_swim_plot builds the table from this list on every build, including
# one that yields zero outlier rows.
SWIM_PLOT_COLUMNS = [
    "outlier_type", "outlier_rank", "motif_category", "SourceDb", "allele_size",
    "sample_id", "family_id", "affected_status", "analysis_status", "sex",
    "phenotype_description", "FirstUnaffectedAlleleSize", "is_above_first_unaffected",
] + SWIM_PLOT_LOCUS_COLUMNS


def _compute_motif_category(motif_size):
    """Bin a motif size into the swim-plot motif_category label.

    Args:
        motif_size: The locus motif size in bp (int-like), or a missing value.

    Returns:
        ``"Unknown"`` if the motif size is missing or non-positive, ``"{size}bp"``
        for sizes <= 24, else ``"25+bp"``.

    A blank motif cell gives MotifSize 0 (``len("")``), which is a supported input
    class. Binning it as ``"0bp"`` would hide those rows: the swim-plot endpoint
    queries only ``1bp``..``24bp``, ``25+bp`` and ``Unknown``, so it would never
    ask for them.
    """
    if is_none_or_nan(motif_size) or int(motif_size) <= 0:
        return "Unknown"
    if int(motif_size) <= 24:
        return f"{int(motif_size)}bp"
    return "25+bp"


def _title_case_status(raw_status):
    """Title-case an affected_status or analysis_status for display, or "Unknown".

    A blank status ("", "NA", "N/A", "null", ...) means nothing was recorded, so it
    becomes "Unknown" rather than the display values "Na", "N/A" or "Null", which
    would show up as extra facet categories.

    Args:
        raw_status: The sample's raw status, or a missing value.

    Returns:
        "Unknown" when the status is missing or blank (``blank_to_none``), otherwise
        the stripped, Title-cased value, e.g. "Affected" or "Partially Solved".
    """
    status = blank_to_none(raw_status)
    return "Unknown" if status is None else str(status).strip().title()


def generate_swim_plot_table(records, sample_lookup, analysis_lookup):
    """Flatten locus records into one swim-plot row per outlier allele.

    Walks the three outlier types and, for each locus record carrying a non-empty
    ``OutlierSampleIds_{outlier_type}`` string, parses it (descending by allele
    size) and emits one row per entry. Each row carries the outlier rank, the
    motif_category bin, the sample's metadata (family / sex / affected /
    analysis / phenotype), the locus's
    FirstUnaffectedAlleleSize for that outlier type, the is_above_first_unaffected
    flag, and the SWIM_PLOT_LOCUS_COLUMNS passthroughs.

    Args:
        records: List of locus record dicts (the build's in-memory representation),
            each carrying OutlierSampleIds_* plus the locus-level passthrough
            columns and FirstUnaffectedAlleleSize_*.
        sample_lookup: Mapping of sample_id -> sample row dict (family_id, sex,
            affected_status, phenotype_description). A sample id absent here yields
            None metadata and an "Unknown" status.
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
            first_unaffected = None if is_none_or_nan(first_unaffected) else int(first_unaffected)

            for rank, (allele_size, sample_id) in enumerate(
                parse_outlier_entries(outlier_value, locus_id=record.get("LocusId"), column=column),
                start=1,
            ):
                sample_row = sample_lookup.get(sample_id)

                row = {
                    "outlier_type": outlier_type,
                    "outlier_rank": rank,
                    "motif_category": motif_category,
                    "SourceDb": record.get("Source"),
                    "allele_size": allele_size,
                    "sample_id": sample_id,
                    "family_id": get_family_id(sample_id, sample_lookup),
                    # The sample row's affected_status, so "possibly affected" is shown as it
                    # was entered; read_sample_metadata's affected_lookup collapses it to
                    # "affected" for the analysis logic.
                    "affected_status": _title_case_status(
                        sample_row.get("affected_status") if sample_row else None),
                    "analysis_status": _title_case_status(analysis_lookup.get(sample_id)),
                    "sex": blank_to_none(sample_row.get("sex")) if sample_row else None,
                    "phenotype_description": (
                        blank_to_none(sample_row.get("phenotype_description")) if sample_row else None),
                    "FirstUnaffectedAlleleSize": first_unaffected,
                    "is_above_first_unaffected": int(is_above_unaffected(allele_size, record, outlier_type)),
                }

                for locus_column in SWIM_PLOT_LOCUS_COLUMNS:
                    row[locus_column] = record.get(locus_column)

                rows.append(row)

    return rows
