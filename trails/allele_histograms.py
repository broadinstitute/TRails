"""Per-locus allele-histogram aggregation and outlier-sample selection.

This module computes, for each tandem-repeat locus, three allele-count
histograms (AllAlleles, ShortAlleles, HemizygousAlleles) plus the matching
lists of "outlier" sample ids (the samples carrying the largest alleles) from
an in-memory genotype matrix.

The genotype matrix is a mapping from locus to a per-sample cell value. Each
cell is one of:
  - "12,40" : a diploid genotype (two allele sizes, in repeat units)
  - "21"    : a hemizygous genotype (a single allele size)
  - ""/"."/"./." : a no-call (skipped)

Histogram routing (faithful port of
str_analysis/combine_single_sample_LPS_to_allele_histograms.py):
  - AllAlleleHistogram        : add the short allele; if diploid also add the long allele.
  - ShortAlleleHistogram      : add the short allele only.
  - HemizygousAlleleHistogram : add the short allele only when the genotype is hemizygous.

Histogram strings are sorted ASCENDING by allele size, formatted "{allele}x:{count}".
Outlier-sample-id strings are sorted DESCENDING by allele size, formatted
"{allele}x:{sample_id}", with the exact skip/stop rules transcribed from
convert_sample_ids_to_string in the reference source.
"""


NO_CALL_VALUES = {"", ".", "./.", "nan", "NaN", "NA"}


def parse_genotype(cell):
    """Parse a single genotype-matrix cell into a tuple of allele sizes.

    Args:
        cell: The raw cell value, e.g. "12,40" (diploid), "21" (hemizygous),
            or "" / "." / "./." (no-call). None is treated as a no-call.

    Returns:
        A tuple of ints (length 2 for diploid, length 1 for hemizygous), or
        None for a no-call / unparseable cell.
    """
    if cell is None:
        return None
    cell = str(cell).strip()
    if cell in NO_CALL_VALUES:
        return None
    parts = cell.split(",")
    if len(parts) not in (1, 2):
        return None
    try:
        return tuple(int(part) for part in parts)
    except ValueError:
        return None


def _record_sample_id(allele_to_sample_ids, allele, sample_id, occurrence_counts, n_outlier_sample_ids):
    """Track a sample id for an allele, capped at 2 * n_outlier_sample_ids occurrences.

    The occurrence count is incremented on every routed allele (so a homozygous
    sample, which routes the same allele twice, counts twice toward the cap),
    matching update_histograms in the reference source. The sample id is
    appended at most once per allele.
    """
    occurrence_counts[allele] += 1
    if occurrence_counts[allele] <= 2 * n_outlier_sample_ids:
        if sample_id not in allele_to_sample_ids[allele]:
            allele_to_sample_ids[allele].append(sample_id)


def accumulate_locus(genotypes_by_sample, n_outlier_sample_ids=10):
    """Accumulate the three histograms and outlier-sample maps for one locus.

    Args:
        genotypes_by_sample: Mapping of sample_id to raw genotype cell string.
        n_outlier_sample_ids: Sample-id occurrences per allele are tracked up to
            2 * n_outlier_sample_ids (homozygous samples count twice).

    Returns:
        A tuple of three (histogram_counts, allele_to_sample_ids) pairs in the
        order (AllAlleles, ShortAlleles, HemizygousAlleles). histogram_counts
        maps allele_size -> count; allele_to_sample_ids maps allele_size -> list
        of sample ids (insertion order preserved).
    """
    import collections

    all_counts = collections.defaultdict(int)
    all_sample_ids = collections.defaultdict(list)
    all_occurrences = collections.defaultdict(int)

    short_counts = collections.defaultdict(int)
    short_sample_ids = collections.defaultdict(list)
    short_occurrences = collections.defaultdict(int)

    hemizygous_counts = collections.defaultdict(int)
    hemizygous_sample_ids = collections.defaultdict(list)
    hemizygous_occurrences = collections.defaultdict(int)

    for sample_id, cell in genotypes_by_sample.items():
        allele_sizes = parse_genotype(cell)
        if allele_sizes is None:
            continue

        is_hemizygous = len(allele_sizes) == 1
        short_allele = min(allele_sizes)
        long_allele = max(allele_sizes)

        # AllAlleleHistogram: short allele always, long allele only if diploid.
        all_counts[short_allele] += 1
        _record_sample_id(all_sample_ids, short_allele, sample_id, all_occurrences, n_outlier_sample_ids)
        if not is_hemizygous:
            all_counts[long_allele] += 1
            _record_sample_id(all_sample_ids, long_allele, sample_id, all_occurrences, n_outlier_sample_ids)

        # ShortAlleleHistogram: short allele only.
        short_counts[short_allele] += 1
        _record_sample_id(short_sample_ids, short_allele, sample_id, short_occurrences, n_outlier_sample_ids)

        # HemizygousAlleleHistogram: short allele only if hemizygous.
        if is_hemizygous:
            hemizygous_counts[short_allele] += 1
            _record_sample_id(hemizygous_sample_ids, short_allele, sample_id, hemizygous_occurrences, n_outlier_sample_ids)

    return (
        (dict(all_counts), dict(all_sample_ids)),
        (dict(short_counts), dict(short_sample_ids)),
        (dict(hemizygous_counts), dict(hemizygous_sample_ids)),
    )


def convert_counts_to_histogram_string(allele_counts):
    """Format an allele-count dict as "{allele}x:{count}" entries, ASCENDING by allele.

    Args:
        allele_counts: Mapping of allele_size -> count.

    Returns:
        A comma-joined string, e.g. "9x:6,10x:4086,11x:178", or "" if empty.
    """
    return ",".join(f"{allele_size}x:{count}" for allele_size, count in sorted(allele_counts.items()))


def convert_sample_ids_to_string(allele_to_sample_ids, n_outlier_sample_ids=10):
    """Format outlier sample ids as "{allele}x:{sample_id}" entries, DESCENDING by allele.

    Faithful transcription of convert_sample_ids_to_string from the reference
    source: iterate alleles largest-first; break at the first allele carried by
    >= n_outlier_sample_ids samples (a "common" allele); emit sample ids in
    ascending order within each allele; stop once more than n_outlier_sample_ids
    ids have been emitted in total.

    Args:
        allele_to_sample_ids: Mapping of allele_size -> list of sample ids.
        n_outlier_sample_ids: The outlier cap controlling the skip/stop rules.

    Returns:
        A comma-joined string, or "" if no allele qualifies.
    """
    output_entries = []
    emitted_count = 0
    for allele, sample_list in sorted(allele_to_sample_ids.items(), key=lambda item: -item[0]):
        if len(sample_list) >= n_outlier_sample_ids:
            break
        for sample_id in sorted(sample_list):
            output_entries.append(f"{allele}x:{sample_id}")
            emitted_count += 1
        if emitted_count > n_outlier_sample_ids:
            break
    return ",".join(output_entries)


def build_histograms_and_outliers(locus_rows, sample_id_list, n_outlier_sample_ids=10):
    """Augment each locus row with histogram and outlier-sample-id columns.

    Args:
        locus_rows: List of dicts, one per locus. Each must carry a "genotypes"
            key mapping sample_id -> raw genotype cell string. Other keys (e.g.
            "trid", "motif") are preserved unchanged.
        sample_id_list: List of sample ids defining which matrix columns to
            consider for each locus. A sample id absent from a row's "genotypes"
            mapping is treated as a no-call for that locus.
        n_outlier_sample_ids: Outlier cap (default 10).

    Returns:
        The same list of locus-row dicts, each augmented in place with:
          AllAlleleHistogram, ShortAlleleHistogram, HemizygousAlleleHistogram,
          OutlierSampleIds_AllAlleles, OutlierSampleIds_ShortAlleles,
          OutlierSampleIds_HemizygousAlleles.
    """
    for locus_row in locus_rows:
        genotypes = locus_row.get("genotypes", {})
        all_pair, short_pair, hemizygous_pair = accumulate_locus(
            {sample_id: genotypes.get(sample_id) for sample_id in sample_id_list},
            n_outlier_sample_ids=n_outlier_sample_ids,
        )

        locus_row["AllAlleleHistogram"] = convert_counts_to_histogram_string(all_pair[0])
        locus_row["ShortAlleleHistogram"] = convert_counts_to_histogram_string(short_pair[0])
        locus_row["HemizygousAlleleHistogram"] = convert_counts_to_histogram_string(hemizygous_pair[0])

        locus_row["OutlierSampleIds_AllAlleles"] = convert_sample_ids_to_string(
            all_pair[1], n_outlier_sample_ids=n_outlier_sample_ids)
        locus_row["OutlierSampleIds_ShortAlleles"] = convert_sample_ids_to_string(
            short_pair[1], n_outlier_sample_ids=n_outlier_sample_ids)
        locus_row["OutlierSampleIds_HemizygousAlleles"] = convert_sample_ids_to_string(
            hemizygous_pair[1], n_outlier_sample_ids=n_outlier_sample_ids)

    return locus_rows
