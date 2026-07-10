"""The TRails outlier/affected analysis — the heart of the build pipeline.

This module turns the per-locus ``OutlierSampleIds_*`` strings (produced by
``allele_histograms``) plus the sample-metadata lookups into the full set of
analysis columns: the First/Second/Third affected and unaffected allele sizes
(both per-sample and per-family), the matching affected phenotypes and sample
ids, and the NumAffectedUnsolved{Samples,Families}AboveUnaffected counts. It
also adds the gene-derived ``pLI`` / ``inheritance`` columns, and exposes the
exact 129-name ordered ``OUTPUT_COLUMNS`` list that defines the loci table.

It is a faithful, standalone port of the corresponding logic in the reference
``analyze_results.py`` (parse_outlier_entries, is_unaffected_or_solved_status,
collect_samples_single_pass, the NumAffectedUnsolved* helpers, add_gene_columns,
add_all_outlier_columns, OUTPUT_COLUMNS). The reference operated on a pandas
DataFrame; here the in-memory representation is a list of dicts (one dict per
locus), so there is no pandas dependency and the algorithms are dict-native.

Two behaviors worth calling out, both load-bearing and matching the reference:

  * ``parse_outlier_entries`` sorts entries DESCENDING by allele size, so the
    "First" affected/unaffected sample is the one carrying the LARGEST allele.
  * ``is_unaffected_or_solved_status`` treats only affected_status="unaffected"
    or analysis_status in {solved, unaffected, probably solved} as
    unaffected/solved. Everything else — including affected="unknown",
    analysis="unsolved" and analysis="partially solved" — is grouped with the
    affecteds (the deliberate "include unknowns alongside affecteds" policy).
"""

import collections


# The three outlier types used throughout the analysis.
OUTLIER_TYPES = ["AllAlleles", "ShortAlleles", "HemizygousAlleles"]

# Outlier-field values that should be treated as "no outlier data".
MISSING_OUTLIER_STRINGS = {"", "NA", "N/A", "nan", "NaN", "None", "null"}


def _is_nan(value):
    """Return True if value is a float NaN (NaN is never equal to itself).

    Avoids a pandas/numpy dependency for the handful of NaN checks this module
    needs. Non-float values (including None and strings) return False.
    """
    return isinstance(value, float) and value != value


def _is_missing(value):
    """Return True if value is None or a float NaN."""
    return value is None or _is_nan(value)


def is_missing_outlier_value(value):
    """Return True if an outlier field should be treated as missing.

    Args:
        value: The raw outlier-column value (a string, None, or a float NaN).

    Returns:
        True for None, a float NaN, or a string that (stripped) is one of the
        recognized missing markers (``""``, ``"."``, ``"NA"``, ``"nan"``, ...).
    """
    if _is_missing(value):
        return True
    if isinstance(value, str) and value.strip() in MISSING_OUTLIER_STRINGS:
        return True
    if isinstance(value, str) and value.strip() == ".":
        return True
    return False


def parse_outlier_entries(value, locus_id=None, column=None):
    """Parse an outlier string into (allele_size, sample_id, purity, methylation) tuples.

    The outlier string is a comma-separated list of
    ``"{allele}x:{sample_id}[:{purity}[:{methylation}]]"`` entries. The trailing
    ``x`` on the allele is stripped before the int conversion. Purity and
    methylation are optional fourth/third fields; a literal ``"."`` (or an
    absent field) becomes None.

    Entries are sorted DESCENDING by allele size, so downstream consumers that
    rely on traversal order (First/Second/Third affected, the above-unaffected
    counts, consecutive-id dedup) always see the largest alleles first.

    Args:
        value: The outlier string, or a missing marker (see
            ``is_missing_outlier_value``).
        locus_id: Optional LocusId, used only to enrich error messages.
        column: Optional column name, used only to enrich error messages.

    Returns:
        A list of ``(allele_size, sample_id, purity, methylation)`` tuples sorted
        descending by allele_size. ``allele_size`` is an int; ``sample_id`` is a
        string; ``purity`` and ``methylation`` are strings or None.

    Raises:
        ValueError: If an entry has fewer than two ``:``-separated fields, or if
            the allele field is not an integer (after stripping ``x``).
    """
    if is_missing_outlier_value(value):
        return []

    if not isinstance(value, str):
        value = str(value)

    entries = []
    for raw_entry in [part for part in value.split(",") if part]:
        fields = raw_entry.split(":")
        if len(fields) < 2:
            raise ValueError(
                f"Malformed outlier entry (expected 'allele:sample_id:...')"
                f"{' in ' + column if column else ''}"
                f"{' for LocusId ' + str(locus_id) if locus_id else ''}: {raw_entry}"
            )
        sample_id = fields[1]
        try:
            allele_size = int(fields[0].replace("x", ""))
        except ValueError as exception:
            raise ValueError(
                f"Malformed allele size in outlier entry"
                f"{' in ' + column if column else ''}"
                f"{' for LocusId ' + str(locus_id) if locus_id else ''}: {raw_entry}"
            ) from exception

        purity = fields[2] if len(fields) > 2 and fields[2] != "." else None
        methylation = fields[3] if len(fields) > 3 and fields[3] != "." else None
        entries.append((allele_size, sample_id, purity, methylation))

    entries.sort(key=lambda entry: entry[0], reverse=True)
    return entries


def is_unaffected_or_solved_status(affected_status, analysis_status):
    """Return True if a sample should be treated as unaffected/solved.

    True iff affected_status is "unaffected" OR analysis_status is one of
    {solved, unaffected, probably solved} (case-insensitive, stripped).

    Intentional: samples with affected_status="unknown", analysis_status="unsolved",
    or analysis_status="partially solved" return False here and are therefore
    grouped with affected/unsolved samples by the outlier-count and ranking paths
    (First/Second/ThirdAffected*, NumAffectedUnsolved*). This is the deliberate
    "include unknowns alongside affecteds" policy — they are candidates that
    haven't been ruled out, so we want them visible in the affected columns.

    Args:
        affected_status: The sample's affected status (any case), or None.
        analysis_status: The sample's analysis status (any case), or None.

    Returns:
        bool: True if the sample is unaffected or solved.
    """
    affected_normalized = str(affected_status or "").strip().lower()
    analysis_normalized = str(analysis_status or "").strip().lower()
    return affected_normalized == "unaffected" or analysis_normalized in {
        "solved",
        "unaffected",
        "probably solved",
    }


def is_affected_unsolved(sample_id, affected_lookup, analysis_lookup):
    """Return True if a sample is affected AND not solved.

    Exact logical complement of ``is_unaffected_or_solved_status`` so the
    phenotype-scoring path and the outlier-count/ranking paths share the same
    affected/unsolved membership. Any sample that is not explicitly "unaffected"
    and not solved is treated as affected/unsolved — including
    affected_status="unknown" samples.

    Args:
        sample_id: The sample identifier to look up.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.

    Returns:
        bool: True if the sample is affected and unsolved.
    """
    affected_normalized = str(affected_lookup.get(sample_id, "")).strip().lower()
    analysis_normalized = str(analysis_lookup.get(sample_id, "")).strip().lower()
    if affected_normalized == "unaffected":
        return False
    if analysis_normalized in ("solved", "probably solved", "unaffected"):
        return False
    return True


def get_population_p99_threshold(row):
    """Return the max available population 99th-percentile threshold, or None.

    Reads whichever of HPRC256_99thPercentile, AoU1027_99thPercentile are
    present (not None / not NaN) and returns their maximum, mirroring the
    long-read cohorts in the server's above-population gate. Returns None if
    none present.

    Args:
        row: A locus record dict that may carry the percentile columns.

    Returns:
        float/int or None: The maximum p99 threshold, or None if unavailable.
    """
    candidate_values = [
        row.get("HPRC256_99thPercentile"),
        row.get("AoU1027_99thPercentile"),
    ]
    valid_values = [value for value in candidate_values if not _is_missing(value)]
    if not valid_values:
        return None
    return max(valid_values)


def is_above_unaffected(allele_size, row, outlier_type):
    """Return True if allele_size exceeds the largest unaffected/solved sample.

    With one row per LocusId, FirstUnaffectedAlleleSize is the max unaffected
    allele, so "above first unaffected" means "above all unaffected". When there
    is no unaffected sample (the column is NULL), any positive allele qualifies.

    Args:
        allele_size: The outlier's allele size.
        row: A locus record dict carrying FirstUnaffectedAlleleSize_* columns.
        outlier_type: One of AllAlleles / ShortAlleles / HemizygousAlleles.

    Returns:
        bool: True if the allele exceeds the first (largest) unaffected allele.
    """
    first_unaffected = row.get(f"FirstUnaffectedAlleleSize_{outlier_type}")
    if _is_missing(first_unaffected):
        return allele_size > 0
    return allele_size > first_unaffected


def is_above_population_p99(allele_size, row):
    """Return True if allele_size exceeds every available cohort's 99th percentile.

    If no cohort has data, the sample is included (returns True). Otherwise the
    allele must strictly exceed each available cohort's 99th percentile.

    Args:
        allele_size: The outlier's allele size.
        row: A locus record dict carrying the *_99thPercentile columns.

    Returns:
        bool: True if the allele exceeds the population p99 thresholds (or no
        data is available).
    """
    candidate_values = [
        row.get("HPRC256_99thPercentile"),
        row.get("AoU1027_99thPercentile"),
    ]
    valid_values = [value for value in candidate_values if not _is_missing(value)]
    if not valid_values:
        return True
    return all(allele_size > value for value in valid_values)


def compute_number_of_affected_unsolved_samples_above_unaffected(row, affected_lookup, analysis_lookup, parsed_by_type):
    """Count affected/unsolved samples above any unaffected/solved sample, per outlier type.

    For each outlier type, walking the entries largest-allele-first:
      1. Collect affected/unsolved sample ids per allele size (a set, so a
         diploid/homozygous sample that appears at two allele sizes is not
         double-counted), stopping at the FIRST unaffected/solved sample.
      2. Union affected sample ids over allele sizes from largest to smallest,
         stopping when an allele has an unaffected/solved sample OR when the
         allele size is <= the population p99 threshold.

    Args:
        row: A locus record dict carrying OutlierSampleIds_* and the percentile
            columns.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.
        parsed_by_type: Mapping of outlier_type -> pre-parsed outlier entries
            (as returned by ``parse_outlier_entries``).

    Returns:
        list: [AllAlleles_count, ShortAlleles_count, HemizygousAlleles_count],
        each an int count of distinct affected/unsolved samples above unaffected,
        or None when that outlier column is missing.
    """
    sample_count_for_each_outlier_type = [None, None, None]
    population_p99_threshold = get_population_p99_threshold(row)

    for outlier_type_index, outlier_type in enumerate(OUTLIER_TYPES):
        column = f"OutlierSampleIds_{outlier_type}"
        if is_missing_outlier_value(row.get(column)) or not row.get(column):
            sample_count_for_each_outlier_type[outlier_type_index] = None
            continue

        allele_size_to_affected_unsolved_samples = collections.defaultdict(set)
        allele_size_to_unaffected_or_solved_sample_count = collections.Counter()

        for allele_size, outlier_sample_id, _, _ in parsed_by_type[outlier_type]:
            if is_unaffected_or_solved_status(
                str(affected_lookup.get(outlier_sample_id, "") or ""),
                str(analysis_lookup.get(outlier_sample_id, "") or ""),
            ):
                allele_size_to_unaffected_or_solved_sample_count[allele_size] += 1
                break
            allele_size_to_affected_unsolved_samples[allele_size].add(outlier_sample_id)

        affected_unsolved_samples = set()
        for allele_size, sample_ids in sorted(
            allele_size_to_affected_unsolved_samples.items(), reverse=True
        ):
            if allele_size_to_unaffected_or_solved_sample_count[allele_size] > 0:
                break
            if population_p99_threshold is not None and allele_size <= population_p99_threshold:
                break
            affected_unsolved_samples.update(sample_ids)

        sample_count_for_each_outlier_type[outlier_type_index] = len(affected_unsolved_samples)

    return sample_count_for_each_outlier_type


def compute_number_of_affected_unsolved_families_above_unaffected(row, affected_lookup, analysis_lookup, sample_lookup, parsed_by_type):
    """Count distinct affected/unsolved families above any unaffected/solved sample.

    Same logic as ``compute_number_of_affected_unsolved_samples_above_unaffected``,
    but counts distinct family_ids instead of sample_ids. Entries whose sample has
    no metadata row, or whose family_id is missing/NaN, are skipped from the
    family set.

    Args:
        row: A locus record dict carrying OutlierSampleIds_* and the percentile
            columns.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.
        sample_lookup: Mapping of sample_id -> sample row dict (used for family_id).
        parsed_by_type: Mapping of outlier_type -> pre-parsed outlier entries.

    Returns:
        list: [AllAlleles_count, ShortAlleles_count, HemizygousAlleles_count],
        each an int count of distinct families, or None when that outlier column
        is missing.
    """
    family_count_for_each_outlier_type = [None, None, None]
    population_p99_threshold = get_population_p99_threshold(row)

    for outlier_type_index, outlier_type in enumerate(OUTLIER_TYPES):
        column = f"OutlierSampleIds_{outlier_type}"
        if is_missing_outlier_value(row.get(column)) or not row.get(column):
            continue

        allele_size_to_affected_families = collections.defaultdict(set)
        allele_size_to_unaffected_or_solved_sample_count = collections.Counter()

        for allele_size, outlier_sample_id, _, _ in parsed_by_type[outlier_type]:
            if is_unaffected_or_solved_status(
                str(affected_lookup.get(outlier_sample_id, "") or ""),
                str(analysis_lookup.get(outlier_sample_id, "") or ""),
            ):
                allele_size_to_unaffected_or_solved_sample_count[allele_size] += 1
                break

            sample_row = sample_lookup.get(outlier_sample_id)
            family_id = sample_row.get("family_id") if sample_row else None
            if family_id and not _is_missing(family_id):
                allele_size_to_affected_families[allele_size].add(family_id)

        affected_unsolved_families = set()
        for allele_size, family_ids in sorted(
            allele_size_to_affected_families.items(), reverse=True
        ):
            if allele_size_to_unaffected_or_solved_sample_count[allele_size] > 0:
                break
            if population_p99_threshold is not None and allele_size <= population_p99_threshold:
                break
            affected_unsolved_families.update(family_ids)

        family_count_for_each_outlier_type[outlier_type_index] = len(affected_unsolved_families)

    return family_count_for_each_outlier_type


def collect_samples_single_pass(outlier_entries, sample_lookup, affected_lookup, analysis_lookup):
    """Collect the first 3 affected and first 2 unaffected samples in one pass.

    The entries must already be sorted descending by allele size (as
    ``parse_outlier_entries`` returns them), so the first affected sample carries
    the largest affected allele. Two dedup modes are computed simultaneously:
    by sample_id (dedup on sample_id) and by family_id (dedup on family_id,
    skipping entries with a missing family_id). For each mode the first 3 affected
    and first 2 unaffected samples are kept.

    Args:
        outlier_entries: A list of ``(allele_size, sample_id, purity, methylation)``
            tuples sorted descending by allele size.
        sample_lookup: Mapping of sample_id -> sample row dict (carries family_id
            and phenotype_description). A sample id absent here is skipped.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.

    Returns:
        dict: ``{"sample_id": {"affected": [...], "unaffected": [...]},
        "family_id": {"affected": [...], "unaffected": [...]}}`` where each item
        is a dict with keys allele_size, sample_id, family_id, phenotype.
    """
    results = {
        "sample_id": {"affected": [], "unaffected": []},
        "family_id": {"affected": [], "unaffected": []},
    }
    seen_sample_ids = {"affected": set(), "unaffected": set()}
    seen_family_ids = {"affected": set(), "unaffected": set()}

    for allele_size, sample_id, _, _ in outlier_entries:
        # A matrix sample with no metadata row is kept (treated as Unknown), not
        # dropped: build_database admits such samples, the affected/unaffected
        # counters already count them, and the affected_lookup/analysis_lookup
        # gates below default to "" (-> affected). Only the family-level ranking
        # is skipped for them, since they have no family_id.
        sample_row = sample_lookup.get(sample_id) or {}
        family_id = sample_row.get("family_id")
        phenotype = sample_row.get("phenotype_description")
        if _is_missing(phenotype):
            phenotype = ""

        category = "unaffected" if is_unaffected_or_solved_status(
            str(affected_lookup.get(sample_id, "") or ""),
            str(analysis_lookup.get(sample_id, "") or ""),
        ) else "affected"
        max_needed = 2 if category == "unaffected" else 3

        entry = {
            "allele_size": allele_size,
            "sample_id": sample_id,
            "family_id": family_id,
            "phenotype": phenotype,
        }

        if sample_id not in seen_sample_ids[category]:
            seen_sample_ids[category].add(sample_id)
            if len(results["sample_id"][category]) < max_needed:
                results["sample_id"][category].append(entry)

        if family_id and not _is_missing(family_id):
            if family_id not in seen_family_ids[category]:
                seen_family_ids[category].add(family_id)
                if len(results["family_id"][category]) < max_needed:
                    results["family_id"][category].append(entry)

        if (len(results["sample_id"]["affected"]) >= 3 and
                len(results["sample_id"]["unaffected"]) >= 2 and
                len(results["family_id"]["affected"]) >= 3 and
                len(results["family_id"]["unaffected"]) >= 2):
            break

    return results


def add_all_outlier_columns(records, sample_lookup, affected_lookup, analysis_lookup):
    """Populate all outlier-analysis columns on each locus record, in place.

    For each outlier type (AllAlleles, ShortAlleles, HemizygousAlleles) and each
    record, parses the OutlierSampleIds_* string once and fills:
      * {First,Second,Third}AffectedAlleleSize_{ot} (+ _ByFamily)
      * {First,Second}UnaffectedAlleleSize_{ot} (+ _ByFamily), kept as nullable
        ints (None preserved, never coerced to 0)
      * {First,Second,Third}AffectedPhenotype_{ot} (+ _ByFamily)
      * {First,Second,Third}AffectedSampleId_{ot} (no _ByFamily)
      * NumAffectedUnsolved{Samples,Families}AboveUnaffected_{ot}

    Every target column is initialized to None on every record, so the loci table
    has a complete, consistent set of columns even where there are no outliers.

    Args:
        records: List of locus record dicts (mutated in place).
        sample_lookup: Mapping of sample_id -> sample row dict.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.

    Returns:
        The same list of records, each augmented with the outlier columns.
    """
    column_names = []
    for suffix in ["", "_ByFamily"]:
        for prefix in ["First", "Second", "Third"]:
            column_names.extend([f"{prefix}AffectedAlleleSize_{ot}{suffix}" for ot in OUTLIER_TYPES])
        for prefix in ["First", "Second"]:
            column_names.extend([f"{prefix}UnaffectedAlleleSize_{ot}{suffix}" for ot in OUTLIER_TYPES])
        for prefix in ["First", "Second", "Third"]:
            column_names.extend([f"{prefix}AffectedPhenotype_{ot}{suffix}" for ot in OUTLIER_TYPES])
    for prefix in ["First", "Second", "Third"]:
        column_names.extend([f"{prefix}AffectedSampleId_{ot}" for ot in OUTLIER_TYPES])
    column_names.extend([f"NumAffectedUnsolvedSamplesAboveUnaffected_{ot}" for ot in OUTLIER_TYPES])
    column_names.extend([f"NumAffectedUnsolvedFamiliesAboveUnaffected_{ot}" for ot in OUTLIER_TYPES])

    for record in records:
        for name in column_names:
            record[name] = None

        parsed_by_type = {
            ot: parse_outlier_entries(
                record.get(f"OutlierSampleIds_{ot}"),
                locus_id=record.get("LocusId"),
                column=f"OutlierSampleIds_{ot}",
            )
            for ot in OUTLIER_TYPES
        }

        num_affected_samples = compute_number_of_affected_unsolved_samples_above_unaffected(
            record, affected_lookup, analysis_lookup, parsed_by_type
        )
        num_affected_families = compute_number_of_affected_unsolved_families_above_unaffected(
            record, affected_lookup, analysis_lookup, sample_lookup, parsed_by_type
        )

        for outlier_type_index, outlier_type in enumerate(OUTLIER_TYPES):
            record[f"NumAffectedUnsolvedSamplesAboveUnaffected_{outlier_type}"] = num_affected_samples[outlier_type_index]
            record[f"NumAffectedUnsolvedFamiliesAboveUnaffected_{outlier_type}"] = num_affected_families[outlier_type_index]

            column = f"OutlierSampleIds_{outlier_type}"
            if is_missing_outlier_value(record.get(column)) or not record.get(column):
                continue

            collected = collect_samples_single_pass(
                parsed_by_type[outlier_type], sample_lookup, affected_lookup, analysis_lookup
            )
            affected_by_sample = collected["sample_id"]["affected"]
            affected_by_family = collected["family_id"]["affected"]
            unaffected_by_sample = collected["sample_id"]["unaffected"]
            unaffected_by_family = collected["family_id"]["unaffected"]

            for rank, prefix in enumerate(["First", "Second", "Third"]):
                record[f"{prefix}AffectedAlleleSize_{outlier_type}"] = (
                    affected_by_sample[rank]["allele_size"] if rank < len(affected_by_sample) else None
                )
                record[f"{prefix}AffectedAlleleSize_{outlier_type}_ByFamily"] = (
                    affected_by_family[rank]["allele_size"] if rank < len(affected_by_family) else None
                )
                record[f"{prefix}AffectedPhenotype_{outlier_type}"] = (
                    affected_by_sample[rank]["phenotype"] if rank < len(affected_by_sample) else None
                )
                record[f"{prefix}AffectedPhenotype_{outlier_type}_ByFamily"] = (
                    affected_by_family[rank]["phenotype"] if rank < len(affected_by_family) else None
                )
                record[f"{prefix}AffectedSampleId_{outlier_type}"] = (
                    affected_by_sample[rank]["sample_id"] if rank < len(affected_by_sample) else None
                )

            for rank, prefix in enumerate(["First", "Second"]):
                record[f"{prefix}UnaffectedAlleleSize_{outlier_type}"] = (
                    unaffected_by_sample[rank]["allele_size"] if rank < len(unaffected_by_sample) else None
                )
                record[f"{prefix}UnaffectedAlleleSize_{outlier_type}_ByFamily"] = (
                    unaffected_by_family[rank]["allele_size"] if rank < len(unaffected_by_family) else None
                )

    return records


# Maps each gene-table source field (as stored by input_tables.read_gene_table)
# to its GeneTable*-prefixed output column.
GENE_TABLE_COLUMN_BY_SOURCE_FIELD = {
    "gene_symbol": "GeneTableGeneSymbol",
    "gene_aliases": "GeneTableGeneAliases",
    "pLI_v2": "GeneTablepLI_v2",
    "pLI_v4": "GeneTablepLI_v4",
    "lof_oe_ci_upper_v4": "GeneTableLoeuf",
    "hgnc_gene_id": "GeneTableHgncGeneId",
    "inheritance": "GeneTableInheritance",
    "disease_category": "GeneTableDiseaseCategory",
    "LLM_phenotype_summary": "GeneTableLLMPhenotypeSummary",
    "sources": "GeneTableSources",
}


def add_gene_columns(records, gene_lookup):
    """Populate the gene-derived ``pLI`` and ``inheritance`` columns, in place.

    For each record, looks up its ``gene_id`` in ``gene_lookup``:
      * ``pLI`` = max of the gene's pLI_v2 and pLI_v4 (whichever are present),
        or None.
      * ``inheritance`` = the gene's inheritance mode, or None.

    Both columns are set to None when there is no gene_lookup, when the record's
    gene_id is missing, or when the gene_id is not in the lookup — so analysis
    runs without a gene table.

    Args:
        records: List of locus record dicts (mutated in place).
        gene_lookup: Mapping of gene_id -> gene info dict (with pLI_v2, pLI_v4,
            inheritance), or None / empty when no gene table was supplied.

    Returns:
        The same list of records, each with pLI and inheritance set.
    """
    for record in records:
        gene_id = record.get("gene_id")
        gene_row = gene_lookup.get(gene_id) if gene_lookup and not _is_missing(gene_id) else None

        if gene_row is None:
            record["pLI"] = None
            record["inheritance"] = None
            for column in GENE_TABLE_COLUMN_BY_SOURCE_FIELD.values():
                record[column] = None
            continue

        pli_v2 = gene_row.get("pLI_v2")
        pli_v4 = gene_row.get("pLI_v4")
        if not _is_missing(pli_v4) and not _is_missing(pli_v2):
            record["pLI"] = max(float(pli_v4), float(pli_v2))
        elif not _is_missing(pli_v4):
            record["pLI"] = float(pli_v4)
        elif not _is_missing(pli_v2):
            record["pLI"] = float(pli_v2)
        else:
            record["pLI"] = None

        record["inheritance"] = gene_row.get("inheritance")

        # Populate the 10 GeneTable* output columns from the gene row so they
        # exist in the loci table (the server filters/selects them directly).
        for source_field, column in GENE_TABLE_COLUMN_BY_SOURCE_FIELD.items():
            value = gene_row.get(source_field)
            record[column] = None if _is_missing(value) else value

    return records


# The exact 129-name ordered loci-table column list. The population-distribution-
# stat columns are added afterward by enrichment, not part of OUTPUT_COLUMNS.
OUTPUT_COLUMNS = [
    # Core locus info
    "LocusId",
    "Motif",
    "CanonicalMotif",
    "IsKnownMotif",
    "IsInMendelianGene",

    # Allele histograms
    "AllAlleleHistogram",
    "ShortAlleleHistogram",
    "HemizygousAlleleHistogram",

    # Outlier sample IDs
    "OutlierSampleIds_AllAlleles",
    "OutlierSampleIds_ShortAlleles",
    "OutlierSampleIds_HemizygousAlleles",

    # Source and coordinates
    "Source",
    "Chrom",
    "Start0Based",
    "End1Based",
    "ReferenceRegion",

    # Disease locus match
    "KnownDiseaseLocus",

    # Locus properties
    "MotifSize",
    "gene_id",
    "gene_region",
    "gene_region_rank",
    "NumRepeatsInReference",

    # Gene table columns (with GeneTable prefix)
    "GeneTableGeneSymbol",
    "GeneTableGeneAliases",
    "GeneTablepLI_v2",
    "GeneTablepLI_v4",
    "GeneTableLoeuf",
    "GeneTableHgncGeneId",
    "GeneTableInheritance",
    "GeneTableDiseaseCategory",
    "GeneTableLLMPhenotypeSummary",
    "GeneTableSources",

    # Combined pLI and inheritance
    "pLI",
    "inheritance",

    # Affected allele sizes (by sample_id)
    "FirstAffectedAlleleSize_AllAlleles",
    "FirstAffectedAlleleSize_ShortAlleles",
    "FirstAffectedAlleleSize_HemizygousAlleles",
    "SecondAffectedAlleleSize_AllAlleles",
    "SecondAffectedAlleleSize_ShortAlleles",
    "SecondAffectedAlleleSize_HemizygousAlleles",
    "ThirdAffectedAlleleSize_AllAlleles",
    "ThirdAffectedAlleleSize_ShortAlleles",
    "ThirdAffectedAlleleSize_HemizygousAlleles",

    # Affected allele sizes (by family_id)
    "FirstAffectedAlleleSize_AllAlleles_ByFamily",
    "FirstAffectedAlleleSize_ShortAlleles_ByFamily",
    "FirstAffectedAlleleSize_HemizygousAlleles_ByFamily",
    "SecondAffectedAlleleSize_AllAlleles_ByFamily",
    "SecondAffectedAlleleSize_ShortAlleles_ByFamily",
    "SecondAffectedAlleleSize_HemizygousAlleles_ByFamily",
    "ThirdAffectedAlleleSize_AllAlleles_ByFamily",
    "ThirdAffectedAlleleSize_ShortAlleles_ByFamily",
    "ThirdAffectedAlleleSize_HemizygousAlleles_ByFamily",

    # Unaffected allele sizes (by sample_id)
    "FirstUnaffectedAlleleSize_AllAlleles",
    "FirstUnaffectedAlleleSize_ShortAlleles",
    "FirstUnaffectedAlleleSize_HemizygousAlleles",
    "SecondUnaffectedAlleleSize_AllAlleles",
    "SecondUnaffectedAlleleSize_ShortAlleles",
    "SecondUnaffectedAlleleSize_HemizygousAlleles",

    # Unaffected allele sizes (by family_id)
    "FirstUnaffectedAlleleSize_AllAlleles_ByFamily",
    "FirstUnaffectedAlleleSize_ShortAlleles_ByFamily",
    "FirstUnaffectedAlleleSize_HemizygousAlleles_ByFamily",
    "SecondUnaffectedAlleleSize_AllAlleles_ByFamily",
    "SecondUnaffectedAlleleSize_ShortAlleles_ByFamily",
    "SecondUnaffectedAlleleSize_HemizygousAlleles_ByFamily",

    # Affected phenotypes (by sample_id)
    "FirstAffectedPhenotype_AllAlleles",
    "FirstAffectedPhenotype_ShortAlleles",
    "FirstAffectedPhenotype_HemizygousAlleles",
    "SecondAffectedPhenotype_AllAlleles",
    "SecondAffectedPhenotype_ShortAlleles",
    "SecondAffectedPhenotype_HemizygousAlleles",
    "ThirdAffectedPhenotype_AllAlleles",
    "ThirdAffectedPhenotype_ShortAlleles",
    "ThirdAffectedPhenotype_HemizygousAlleles",

    # Affected phenotypes (by family_id)
    "FirstAffectedPhenotype_AllAlleles_ByFamily",
    "FirstAffectedPhenotype_ShortAlleles_ByFamily",
    "FirstAffectedPhenotype_HemizygousAlleles_ByFamily",
    "SecondAffectedPhenotype_AllAlleles_ByFamily",
    "SecondAffectedPhenotype_ShortAlleles_ByFamily",
    "SecondAffectedPhenotype_HemizygousAlleles_ByFamily",
    "ThirdAffectedPhenotype_AllAlleles_ByFamily",
    "ThirdAffectedPhenotype_ShortAlleles_ByFamily",
    "ThirdAffectedPhenotype_HemizygousAlleles_ByFamily",

    # Affected sample IDs (no ByFamily variant)
    "FirstAffectedSampleId_AllAlleles",
    "FirstAffectedSampleId_ShortAlleles",
    "FirstAffectedSampleId_HemizygousAlleles",
    "SecondAffectedSampleId_AllAlleles",
    "SecondAffectedSampleId_ShortAlleles",
    "SecondAffectedSampleId_HemizygousAlleles",
    "ThirdAffectedSampleId_AllAlleles",
    "ThirdAffectedSampleId_ShortAlleles",
    "ThirdAffectedSampleId_HemizygousAlleles",

    # Num affected above unaffected
    "NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles",
    "NumAffectedUnsolvedSamplesAboveUnaffected_ShortAlleles",
    "NumAffectedUnsolvedSamplesAboveUnaffected_HemizygousAlleles",
    "NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles",
    "NumAffectedUnsolvedFamiliesAboveUnaffected_ShortAlleles",
    "NumAffectedUnsolvedFamiliesAboveUnaffected_HemizygousAlleles",

    # Phenotype scores (denormalized from per_locus_phenotype_scores)
    "MaxGenePhenoSim_AllAlleles",
    "MaxGenePhenoSim_ShortAlleles",
    "MaxGenePhenoSim_HemizygousAlleles",
    "SumPairwiseSim_AllAlleles",
    "SumPairwiseSim_ShortAlleles",
    "SumPairwiseSim_HemizygousAlleles",

    # Population reference data
    "AoU1027_99thPercentile",
    "AoU1027_MaxAllele",
    "HPRC256_99thPercentile",
    "HPRC256_MaxAllele",
    "TenK10K_99thPercentile",
    "TenK10K_MaxAllele",

    # Additional TRExplorer columns
    "AoU1027_OE_Length",
    "AoU1027_OE_LengthPercentile",
    "AoU1027_Stdev",
    "AoU1027_StdevRankByMotif",
    "AoU1027_StdevRankTotalNumberByMotif",
    "HPRC256_Stdev",
    "HPRC256_StdevRankByMotif",
    "HPRC256_StdevRankTotalNumberByMotif",
    "HPRC256_StdevPercentile",
    "AoU1027_StdevPercentile",
    "TenK10K_Stdev",

    # TRExplorer locus info
    "TRExplorerLocusId",
    "TRExplorerLocusJaccardSimilarity",
    "TRExplorerMotif",
    "TRExplorerReferenceRegion",
    "TRExplorerReferenceRepeatPurity",
    "TRExplorerSource",

    # Other annotations
    "NonCodingAnnotations",
    "RepeatMaskerIntervals",
    "VariationClusterSizeDiff",
]
