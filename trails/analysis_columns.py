"""The TRails outlier/affected analysis — the heart of the build pipeline.

This module turns the per-locus ``OutlierSampleIds_*`` strings (produced by
``allele_histograms``) plus the sample-metadata lookups into the full set of
analysis columns: the First/Second/Third affected and unaffected allele sizes
(both per-sample and per-family), the matching affected phenotypes and sample
ids, the NumAffectedUnsolved{Samples,Families}AboveUnaffected counts, and the
rows of the affected_unsolved_outliers_above_unaffected table that lists the
samples those counts count. It also adds the gene-derived ``pLI`` /
``inheritance`` columns, and exposes the exact 131-name ordered
``OUTPUT_COLUMNS`` list that defines the loci table.

The in-memory representation is a list of dicts (one dict per locus), so there is
no pandas dependency and the algorithms are dict-native.

Two behaviors worth calling out, both load-bearing:

  * ``parse_outlier_entries`` sorts entries DESCENDING by allele size, so the
    "First" affected/unaffected sample is the one carrying the LARGEST allele.
  * ``is_unaffected_or_solved_status`` treats only affected_status="unaffected"
    or analysis_status in {solved, unaffected, probably solved} as
    unaffected/solved. Everything else — including affected="unknown",
    analysis="unsolved" and analysis="partially solved" — is grouped with the
    affecteds (the deliberate "include unknowns alongside affecteds" policy).
"""


# The three outlier types used throughout the analysis.
OUTLIER_TYPES = ["AllAlleles", "ShortAlleles", "HemizygousAlleles"]

# Lower-cased, stripped cell values that mean "nothing was recorded", the one list for every
# input: a sample-metadata column (family_id, sex, affected_status, analysis_status,
# maternal_id, ...), an annotation column of the matrix (input_tables._coerce_annotation_value)
# and an OutlierSampleIds_* field. "." is TRails' own missing marker, as in the matrix's
# no-call cells and the gene table. Defined here, in a module with no third-party imports, so
# every module that reads these inputs can share it. Most code tests a value with blank_to_none
# instead.
MISSING_VALUE_SPELLINGS = {"", ".", "na", "n/a", "nan", "none", "null"}

# The cohorts whose 99th percentiles gate the NumAffectedUnsolved* counts and the phenotype
# scores by default (get_population_p99_threshold): the long-read
# ones, the same cohorts as the server's long-read above-population filter.
# build_database --population-gate-cohorts replaces it for one build, e.g. to add TenK10K
# (short-read ExpansionHunter calls) for a short-read callset. The build records the cohorts
# it used in the metadata table (result_database.write_population_gate_cohorts).
POPULATION_P99_GATE_COHORTS = ("HPRC256", "AoU1027")
POPULATION_GATE_COHORT_CHOICES = POPULATION_P99_GATE_COHORTS + ("TenK10K",)


def _is_nan(value):
    """Return True if value is a float NaN (NaN is never equal to itself).

    Avoids a pandas/numpy dependency for the handful of NaN checks this module
    needs. Non-float values (including None and strings) return False.
    """
    return isinstance(value, float) and value != value


def is_none_or_nan(value):
    """Return True if value is None or a float NaN."""
    return value is None or _is_nan(value)


def blank_to_none(value):
    """Return None for a missing or blank value, else the value unchanged.

    read_sample_metadata reads with keep_default_na=False, so an empty metadata cell
    arrives as '' and a cell written as NA, N/A, null, "." ... arrives as that string.
    This maps None, a float NaN and every MISSING_VALUE_SPELLINGS spelling to None, so every
    table stores a blank metadata field as NULL, and samples whose family_id cells all
    say "NA" are not counted as one family.

    Args:
        value: A metadata value (a string, None, or a float NaN).

    Returns:
        None for None, a float NaN, or a string whose stripped, lower-cased form is in
        MISSING_VALUE_SPELLINGS (which includes the empty string); otherwise value.
    """
    if is_none_or_nan(value):
        return None
    if isinstance(value, str) and value.strip().lower() in MISSING_VALUE_SPELLINGS:
        return None
    return value


def normalize_affected_status_for_logic(value):
    """Normalize an affected-status cell so comparisons don't depend on its spelling.

    The build and the server both read the sample table through
    input_tables.read_sample_metadata, which applies this, so they agree on what
    "Possibly Affected" means.

    Args:
        value: An affected-status cell, which may be None, NaN or any of the blank
            spellings in MISSING_VALUE_SPELLINGS.

    Returns:
        None for a blank value (``blank_to_none``), "affected" for "possibly affected",
        and the lower-cased stripped value otherwise.
    """
    if blank_to_none(value) is None:
        return None
    normalized = str(value).strip().lower()
    return "affected" if normalized == "possibly affected" else normalized


def is_missing_outlier_value(value):
    """Return True if an outlier field should be treated as missing.

    Args:
        value: The raw outlier-column value (a string, None, or a float NaN).

    Returns:
        True for None, a float NaN, or a string that is a MISSING_VALUE_SPELLINGS
        spelling (``""``, ``"."``, ``"NA"``, ``"nan"``, ... in any case).
    """
    return blank_to_none(value) is None


def parse_outlier_entries(value, locus_id=None, column=None):
    """Parse an outlier string into (allele_size, sample_id) pairs.

    The outlier string is a comma-separated list of ``"{allele}x:{sample_id}"``
    entries, the format ``allele_histograms.convert_sample_ids_to_string`` writes.
    The trailing ``x`` on the allele is stripped before the int conversion.

    Entries are sorted DESCENDING by allele size, so the consumers that rely on
    that order see the largest alleles first: the First/Second/Third affected and
    unaffected ranking (``collect_samples_single_pass``), the walk that stops at
    the first unaffected or sub-threshold allele
    (``find_affected_unsolved_samples_above_unaffected``), and keeping each sample
    at its largest allele.

    Args:
        value: The outlier string, or a missing marker (see
            ``is_missing_outlier_value``).
        locus_id: Optional LocusId, used only to enrich error messages.
        column: Optional column name, used only to enrich error messages.

    Returns:
        A list of ``(allele_size, sample_id)`` pairs sorted descending by
        allele_size. ``allele_size`` is an int; ``sample_id`` is a string.

    Raises:
        ValueError: If an entry does not have exactly two ``:``-separated fields
            (a sample_id cannot contain ``:``), or if the allele field is not an
            integer (after stripping ``x``).
    """
    if is_missing_outlier_value(value):
        return []

    if not isinstance(value, str):
        value = str(value)

    entries = []
    for raw_entry in [part for part in value.split(",") if part]:
        fields = raw_entry.split(":")
        if len(fields) != 2:
            raise ValueError(
                f"Malformed outlier entry (expected 'allele:sample_id')"
                f"{' in ' + column if column else ''}"
                f"{' for LocusId ' + str(locus_id) if locus_id else ''}: {raw_entry}"
            )
        allele_field, sample_id = fields
        try:
            allele_size = int(allele_field.replace("x", ""))
        except ValueError as exception:
            raise ValueError(
                f"Malformed allele size in outlier entry"
                f"{' in ' + column if column else ''}"
                f"{' for LocusId ' + str(locus_id) if locus_id else ''}: {raw_entry}"
            ) from exception
        entries.append((allele_size, sample_id))

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


def parse_population_gate_cohorts(comma_separated_cohorts):
    """Split a --population-gate-cohorts value into the cohort names it lists.

    Args:
        comma_separated_cohorts: Comma-separated cohort names, e.g. "HPRC256,AoU1027,TenK10K".
            Whitespace around a name is ignored.

    Returns:
        tuple: The cohort names, in the order given.

    Raises:
        ValueError: If a name is not one of POPULATION_GATE_COHORT_CHOICES.
    """
    population_gate_cohorts = tuple(cohort.strip() for cohort in comma_separated_cohorts.split(","))
    unknown_cohorts = set(population_gate_cohorts) - set(POPULATION_GATE_COHORT_CHOICES)
    if unknown_cohorts:
        raise ValueError(
            f"Unknown --population-gate-cohorts value(s): {', '.join(sorted(unknown_cohorts))}. "
            f"Choices: {', '.join(POPULATION_GATE_COHORT_CHOICES)}")
    return population_gate_cohorts


def get_population_p99_threshold(row, population_gate_cohorts=POPULATION_P99_GATE_COHORTS):
    """Return the max available population 99th-percentile threshold, or None.

    Reads the {cohort}_99thPercentile column of each of population_gate_cohorts (by
    default HPRC256 and AoU1027, the long-read cohorts of the server's above-population
    gate) and returns the maximum of those present (not None / not NaN). Returns None
    if none present. The same all-alleles percentiles gate every outlier type: the loci
    table has no ShortAlleles- or HemizygousAlleles-specific percentile columns, and the
    server's population filter compares against these same ones.

    Args:
        row: A locus record dict that may carry the percentile columns.
        population_gate_cohorts: The cohorts to take the 99th percentiles from.

    Returns:
        float/int or None: The maximum p99 threshold, or None if unavailable.
    """
    candidate_values = [row.get(f"{cohort}_99thPercentile") for cohort in population_gate_cohorts]
    valid_values = [value for value in candidate_values if not is_none_or_nan(value)]
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
    if is_none_or_nan(first_unaffected):
        return allele_size > 0
    return allele_size > first_unaffected


def parse_outlier_entries_by_type(record):
    """Parse each of a locus record's three OutlierSampleIds_* strings.

    Args:
        record: A locus record dict carrying OutlierSampleIds_* and LocusId.

    Returns:
        dict: outlier_type -> the ``parse_outlier_entries`` list for that type's column.
    """
    return {
        outlier_type: parse_outlier_entries(
            record.get(f"OutlierSampleIds_{outlier_type}"),
            locus_id=record.get("LocusId"),
            column=f"OutlierSampleIds_{outlier_type}",
        )
        for outlier_type in OUTLIER_TYPES
    }


def find_affected_unsolved_samples_above_unaffected(row, affected_lookup, analysis_lookup, parsed_by_type,
                                                    population_gate_cohorts=POPULATION_P99_GATE_COHORTS):
    """Find the affected/unsolved samples above any unaffected/solved sample, per outlier type.

    For each outlier type, the entries are walked largest-allele-first, keeping each
    affected/unsolved sample at the first (largest) allele it appears with, so a
    sample whose two alleles both pass counts once. The walk stops at the first
    allele size <= the population p99 threshold, or at the FIRST unaffected/solved
    sample, whose allele size is then dropped too: an affected sample tied with it
    is not "above" it.

    These are the samples NumAffectedUnsolvedSamplesAboveUnaffected counts, and their
    families are the ones NumAffectedUnsolvedFamiliesAboveUnaffected counts (see
    ``count_affected_unsolved_samples_and_families``). The
    affected_unsolved_outliers_above_unaffected table lists them (see
    ``list_affected_unsolved_outliers``), so the results server can group loci by
    sample or family, and phenotype_scoring scores the ones that have HPO terms.

    Args:
        row: A locus record dict carrying OutlierSampleIds_* and the percentile
            columns.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.
        parsed_by_type: Mapping of outlier_type -> pre-parsed outlier entries
            (as returned by ``parse_outlier_entries``).
        population_gate_cohorts: The cohorts whose 99th percentiles gate the alleles
            (see ``get_population_p99_threshold``).

    Returns:
        list: [AllAlleles, ShortAlleles, HemizygousAlleles], each a dict mapping
        sample_id to the sample's largest allele size above the thresholds, in
        descending allele order (samples sharing an allele size in the order the
        entries list them), or None when that outlier column is missing.
    """
    samples_for_each_outlier_type = [None, None, None]
    population_p99_threshold = get_population_p99_threshold(row, population_gate_cohorts)

    for outlier_type_index, outlier_type in enumerate(OUTLIER_TYPES):
        if is_missing_outlier_value(row.get(f"OutlierSampleIds_{outlier_type}")):
            continue

        sample_to_allele_size = {}
        for allele_size, outlier_sample_id in parsed_by_type[outlier_type]:
            if population_p99_threshold is not None and allele_size <= population_p99_threshold:
                break
            if is_unaffected_or_solved_status(
                affected_lookup.get(outlier_sample_id), analysis_lookup.get(outlier_sample_id),
            ):
                sample_to_allele_size = {sample_id: size for sample_id, size in sample_to_allele_size.items()
                                         if size != allele_size}
                break
            sample_to_allele_size.setdefault(outlier_sample_id, allele_size)

        samples_for_each_outlier_type[outlier_type_index] = sample_to_allele_size

    return samples_for_each_outlier_type


def get_family_id(sample_id, sample_lookup):
    """Return the sample's family_id from the sample metadata, or None when it has none.

    Args:
        sample_id: The sample identifier.
        sample_lookup: Mapping of sample_id -> sample row dict.

    Returns:
        The family_id, or None for a sample with no metadata row or a blank/NaN family_id
        (see ``blank_to_none``).
    """
    sample_row = sample_lookup.get(sample_id)
    return blank_to_none(sample_row.get("family_id")) if sample_row else None


def count_affected_unsolved_samples_and_families(samples_for_each_outlier_type, sample_lookup):
    """Count the samples and distinct families find_affected_unsolved_samples_above_unaffected found.

    These are the NumAffectedUnsolvedSamplesAboveUnaffected_* and
    NumAffectedUnsolvedFamiliesAboveUnaffected_* values. A sample without a family_id
    counts toward the samples but not the families.

    Args:
        samples_for_each_outlier_type: The list returned by
            ``find_affected_unsolved_samples_above_unaffected``.
        sample_lookup: Mapping of sample_id -> sample row dict (used for family_id).

    Returns:
        tuple: (sample counts, family counts), each a list
        [AllAlleles_count, ShortAlleles_count, HemizygousAlleles_count] with None
        where the outlier column is missing.
    """
    sample_counts = []
    family_counts = []
    for samples in samples_for_each_outlier_type:
        if samples is None:
            sample_counts.append(None)
            family_counts.append(None)
            continue
        sample_counts.append(len(samples))
        family_counts.append(len({get_family_id(sample_id, sample_lookup) for sample_id in samples} - {None}))
    return sample_counts, family_counts


def list_affected_unsolved_outliers(locus_id, samples_for_each_outlier_type, sample_lookup):
    """Return one locus's rows of the affected_unsolved_outliers_above_unaffected table.

    Args:
        locus_id: The locus's LocusId.
        samples_for_each_outlier_type: The list returned by
            ``find_affected_unsolved_samples_above_unaffected``.
        sample_lookup: Mapping of sample_id -> sample row dict (used for family_id).

    Returns:
        list: One dict per (outlier type, sample), keyed by the
        ``result_database.AFFECTED_UNSOLVED_OUTLIERS_COLUMN_TYPES`` columns (LocusId,
        outlier_type, sample_id, family_id, allele_size). family_id is None for a
        sample without one.
    """
    return [{"LocusId": locus_id, "outlier_type": outlier_type, "sample_id": sample_id,
             "family_id": get_family_id(sample_id, sample_lookup), "allele_size": allele_size}
            for outlier_type, samples in zip(OUTLIER_TYPES, samples_for_each_outlier_type)
            for sample_id, allele_size in sorted((samples or {}).items())]


def collect_samples_single_pass(outlier_entries, sample_lookup, affected_lookup, analysis_lookup):
    """Collect the first 3 affected and first 2 unaffected samples in one pass.

    The entries must already be sorted descending by allele size (as
    ``parse_outlier_entries`` returns them), so the first affected sample carries
    the largest affected allele. Two dedup modes are computed simultaneously:
    by sample_id (dedup on sample_id) and by family_id (dedup on family_id,
    skipping entries without one, as ``get_family_id`` decides). For each mode the
    first 3 affected and first 2 unaffected samples are kept.

    Args:
        outlier_entries: A list of ``(allele_size, sample_id)`` pairs sorted
            descending by allele size.
        sample_lookup: Mapping of sample_id -> sample row dict (carries family_id
            and phenotype_description). A sample id absent here is still ranked
            (treated as Unknown); it is only left out of the family-level ranking,
            because it has no family_id.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.

    Returns:
        dict: ``{"sample_id": {"affected": [...], "unaffected": [...]},
        "family_id": {"affected": [...], "unaffected": [...]}}`` where each item
        is a dict with keys allele_size, sample_id, family_id, phenotype. family_id
        and phenotype are None when blank (``get_family_id``, ``blank_to_none``).
    """
    results = {
        "sample_id": {"affected": [], "unaffected": []},
        "family_id": {"affected": [], "unaffected": []},
    }
    seen_sample_ids = {"affected": set(), "unaffected": set()}
    seen_family_ids = {"affected": set(), "unaffected": set()}

    for allele_size, sample_id in outlier_entries:
        # A matrix sample with no metadata row is kept (treated as Unknown), not
        # dropped: build_database admits such samples, the affected/unaffected
        # counters already count them, and the affected_lookup/analysis_lookup
        # gates below default to "" (-> affected). Only the family-level ranking
        # is skipped for them, since they have no family_id.
        sample_row = sample_lookup.get(sample_id) or {}
        family_id = get_family_id(sample_id, sample_lookup)
        phenotype = blank_to_none(sample_row.get("phenotype_description"))

        category = "unaffected" if is_unaffected_or_solved_status(
            affected_lookup.get(sample_id), analysis_lookup.get(sample_id),
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

        if family_id is not None:
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


def add_all_outlier_columns(records, sample_lookup, affected_lookup, analysis_lookup,
                            population_gate_cohorts=POPULATION_P99_GATE_COHORTS):
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
        population_gate_cohorts: The cohorts whose 99th percentiles gate the
            NumAffectedUnsolved* counts (see ``get_population_p99_threshold``).

    Returns:
        list: The affected_unsolved_outliers_above_unaffected rows of every record
        (see ``list_affected_unsolved_outliers``): one per (locus, outlier type,
        sample) that NumAffectedUnsolvedSamplesAboveUnaffected counts, so the table
        always agrees with the counts.
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

    affected_unsolved_outlier_rows = []
    for record in records:
        for name in column_names:
            record[name] = None

        parsed_by_type = parse_outlier_entries_by_type(record)

        # Find the samples NumAffectedUnsolved{Samples,Families}AboveUnaffected count, once for
        # all outlier types, and record them for the affected_unsolved_outliers_above_unaffected table
        affected_unsolved_samples = find_affected_unsolved_samples_above_unaffected(
            record, affected_lookup, analysis_lookup, parsed_by_type, population_gate_cohorts
        )
        num_affected_samples, num_affected_families = count_affected_unsolved_samples_and_families(
            affected_unsolved_samples, sample_lookup
        )
        affected_unsolved_outlier_rows.extend(list_affected_unsolved_outliers(
            record.get("LocusId"), affected_unsolved_samples, sample_lookup
        ))

        for outlier_type_index, outlier_type in enumerate(OUTLIER_TYPES):
            record[f"NumAffectedUnsolvedSamplesAboveUnaffected_{outlier_type}"] = num_affected_samples[outlier_type_index]
            record[f"NumAffectedUnsolvedFamiliesAboveUnaffected_{outlier_type}"] = num_affected_families[outlier_type_index]

            if is_missing_outlier_value(record.get(f"OutlierSampleIds_{outlier_type}")):
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

    return affected_unsolved_outlier_rows


# Maps each gene-table source field to its GeneTable*-prefixed output column. The fields are also
# the ones input_tables.read_gene_table keeps from the gene table.
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
    """Fill the GeneTable* columns, ``pLI`` and ``inheritance`` from the gene table, in place.

    For each record, looks up its ``gene_id`` in ``gene_lookup`` and sets:
      * the ten GENE_TABLE_COLUMN_BY_SOURCE_FIELD columns (``GeneTableGeneSymbol``,
        ``GeneTablepLI_v4``, ``GeneTableInheritance``, ...) to the gene row's
        fields, with a blank (NaN) field as None.
      * ``pLI`` = max of the gene's pLI_v2 and pLI_v4 (whichever are present),
        or None.
      * ``inheritance`` = the gene's inheritance mode, or None. It always equals
        ``GeneTableInheritance``.

    All of them are set to None when there is no gene_lookup, when the record's
    gene_id is missing, or when the gene_id is not in the lookup, so the loci
    table has the columns even for a build without a gene table.

    Args:
        records: List of locus record dicts (mutated in place).
        gene_lookup: Mapping of gene_id -> gene info dict (as read by
            input_tables.read_gene_table), or None / empty when no gene table was
            supplied.
    """
    for record in records:
        gene_id = record.get("gene_id")
        gene_row = gene_lookup.get(gene_id) if gene_lookup and not is_none_or_nan(gene_id) else None

        if gene_row is None:
            record["pLI"] = None
            for column in GENE_TABLE_COLUMN_BY_SOURCE_FIELD.values():
                record[column] = None
        else:
            pli_v2 = gene_row.get("pLI_v2")
            pli_v4 = gene_row.get("pLI_v4")
            if not is_none_or_nan(pli_v4) and not is_none_or_nan(pli_v2):
                record["pLI"] = max(float(pli_v4), float(pli_v2))
            elif not is_none_or_nan(pli_v4):
                record["pLI"] = float(pli_v4)
            elif not is_none_or_nan(pli_v2):
                record["pLI"] = float(pli_v2)
            else:
                record["pLI"] = None

            # Populate the 10 GeneTable* output columns from the gene row so they exist in the
            # loci table (the server filters/selects them directly). read_gene_table does not
            # pass keep_default_na=False, so a blank cell arrives as a float NaN; store it as None.
            for source_field, column in GENE_TABLE_COLUMN_BY_SOURCE_FIELD.items():
                value = gene_row.get(source_field)
                record[column] = None if is_none_or_nan(value) else value

        record["inheritance"] = record["GeneTableInheritance"]


# The ordered 131-name loci-table schema. Other per-locus annotation columns the matrix
# supplies (e.g. HPRC256_Mode) are appended after these by result_database.write_loci_table.
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
    "VariationCluster",
    "VariationClusterFilterReason",
]
