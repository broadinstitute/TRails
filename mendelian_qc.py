"""Trio Mendelian-violation QC computed directly from the in-memory genotype matrix.

This module ports the Mendelian-violation analysis of the original
``compute_mendelian_violations.py``. The only change versus the original is the
*data source*: rather than loading one locus-per-sample (LPS) genotype file per
trio member from disk/GCS, TRails reads the genotypes straight out of the
in-memory repeat-copy-numbers matrix (the ``locus_rows`` produced by
``input_tables.read_repeat_copy_numbers``) plus the trio relationships from the
sample-metadata table.

Everything else -- ``check_violation``, ``allele_matches_any``,
``get_chrom_category``, ``get_motif_size_category``, ``MOTIF_SIZE_CATEGORIES``,
the per-chromosome dispatch (chrY uses child+father, chrM uses child+mother,
autosome/chrX uses all three), the "require at least 2 distinct allele sizes"
filter, and the two output-table schemas -- is a faithful, verbatim port of the
original behavior.

Two TRails-specific additions handle the fact that a matrix cell can be a
no-call (the per-sample LPS files never were):
  * ``parse_genotype`` returns ``None`` for empty / ``.`` / ``./.`` cells (and
    for malformed cells), and such loci are skipped for the affected member.
  * The matrix is the single source of genotypes, so the "filter mother to loci
    present in the child" step of the original is automatic (we iterate the
    child's genotypes per locus).

All functions are pure -- they take their inputs as arguments and return values
without mutating shared module state.
"""

from collections import Counter

from motif_utilities import compute_canonical_motif, generate_all_canonical_motifs


# Motif size categories, in display/column order (verbatim from the source).
MOTIF_SIZE_CATEGORIES = ["1bp", "2bp", "3bp", "4bp", "5bp", "6bp", "7-24bp", "25+bp"]

# All canonical motifs of size 1-6bp, sorted for deterministic column ordering
# in the per-motif output table. (The source generated this as a sorted list;
# generate_all_canonical_motifs here returns a set, so we sort it once.)
ALL_CANONICAL_MOTIFS = sorted(generate_all_canonical_motifs(6))


def parse_genotype(cell):
    """Parse a matrix genotype cell into a tuple of integer allele sizes.

    Reproduces the LPS parse rule from the original ``load_lps_file``: the cell
    is split on ``,``; a 2-element split is a diploid call (a 2-tuple), a
    1-element split is a hemizygous call (a 1-tuple), and anything else is
    skipped. As a TRails addition, no-call cells (empty, ``.``, ``./.``) and
    cells that fail to parse as integers are skipped as well (the matrix can
    contain no-calls; the per-sample LPS files could not).

    Args:
        cell: The raw genotype string for one sample at one locus, e.g.
            ``"12,40"`` (diploid), ``"21"`` (hemizygous), or ``""`` / ``"."`` /
            ``"./."`` (no-call).

    Returns:
        A tuple of ints -- ``(short, long)`` for a diploid call or
        ``(allele,)`` for a hemizygous call -- or ``None`` for a no-call or
        malformed cell.
    """
    if cell is None:
        return None
    text = str(cell).strip()
    if text in ("", ".", "./."):
        return None

    allele_parts = text.split(",")
    if len(allele_parts) not in (1, 2):
        return None
    try:
        return tuple(int(part) for part in allele_parts)
    except ValueError:
        return None


def allele_matches_any(allele, genotype, threshold):
    """Return True if ``allele`` is within ``threshold`` of any allele in ``genotype``.

    A match means the two allele sizes differ by strictly fewer than
    ``threshold`` repeats (``abs(a - b) < threshold``).

    Args:
        allele: An integer allele size.
        genotype: A tuple of integer allele sizes to match against.
        threshold: Alleles match if they differ by strictly fewer than this many
            repeats.

    Returns:
        True if any allele in ``genotype`` matches ``allele``.
    """
    return any(abs(allele - other_allele) < threshold for other_allele in genotype)


def check_violation(child_genotype, mother_genotype, father_genotype, threshold):
    """Return True if the child's genotype violates Mendelian inheritance.

    A genotype is consistent with Mendelian inheritance if the child could have
    inherited one allele from the mother and one from the father. Alleles are
    considered matching if they differ by strictly fewer than ``threshold``
    repeats.

    Handles hemizygous calls (1-element tuples) for sex chromosomes:
      * Hemizygous child (e.g. a son on chrX): the single allele must come from
        the mother.
      * Hemizygous father on chrX: a daughter inherits his single allele.
    (chrY and chrM are handled by the per-chromosome dispatch in
    ``compute_mendelian_violations`` and never reach this function.)

    Args:
        child_genotype: Tuple of alleles for the child (1 or 2 elements).
        mother_genotype: Tuple of alleles for the mother (1 or 2 elements).
        father_genotype: Tuple of alleles for the father (1 or 2 elements).
        threshold: Alleles match if they differ by strictly fewer than this many
            repeats.

    Returns:
        True if the genotype is a Mendelian violation, False otherwise.
    """
    # Case 1: the child is hemizygous (a son on chrX; chrY is handled separately
    # and never reaches here), so the single X allele is inherited from the
    # mother regardless of the mother's ploidy.
    if len(child_genotype) == 1:
        return not allele_matches_any(child_genotype[0], mother_genotype, threshold)

    # Case 2: the child is diploid.
    first_child_allele, second_child_allele = child_genotype

    # Mother hemizygous (unusual, but handle it): the child must get the
    # mother's single allele, the other from the father.
    if len(mother_genotype) == 1:
        mother_allele = mother_genotype[0]
        return not (
            (abs(first_child_allele - mother_allele) < threshold
             and allele_matches_any(second_child_allele, father_genotype, threshold))
            or
            (abs(second_child_allele - mother_allele) < threshold
             and allele_matches_any(first_child_allele, father_genotype, threshold))
        )

    # Father hemizygous (e.g. chrX in a daughter): the child must get the
    # father's single allele, the other from the mother.
    if len(father_genotype) == 1:
        father_allele = father_genotype[0]
        return not (
            (abs(first_child_allele - father_allele) < threshold
             and allele_matches_any(second_child_allele, mother_genotype, threshold))
            or
            (abs(second_child_allele - father_allele) < threshold
             and allele_matches_any(first_child_allele, mother_genotype, threshold))
        )

    # Case 3: all diploid (autosomal / female chrX).
    first_mother_allele, second_mother_allele = mother_genotype
    first_father_allele, second_father_allele = father_genotype

    # Assignment 1: first child allele from mother, second from father.
    if (
        (abs(first_child_allele - first_mother_allele) < threshold
         or abs(first_child_allele - second_mother_allele) < threshold)
        and
        (abs(second_child_allele - first_father_allele) < threshold
         or abs(second_child_allele - second_father_allele) < threshold)
    ):
        return False

    # Assignment 2: first child allele from father, second from mother.
    if (
        (abs(first_child_allele - first_father_allele) < threshold
         or abs(first_child_allele - second_father_allele) < threshold)
        and
        (abs(second_child_allele - first_mother_allele) < threshold
         or abs(second_child_allele - second_mother_allele) < threshold)
    ):
        return False

    return True


def get_chrom_category(trid):
    """Categorize a locus id by chromosome type.

    The chromosome is taken from the first ``_``-delimited field of the trid, or
    the first ``-``-delimited field if there is no underscore, or the whole
    string otherwise (matching the original).

    Args:
        trid: Locus id, expected to start with a chromosome name (e.g.
            ``chr1_...`` or ``chrX-...``).

    Returns:
        One of ``"autosome"``, ``"chrX"``, ``"chrY"`` or ``"chrM"``.
    """
    if "_" in trid:
        chrom = trid.split("_")[0].lower()
    elif "-" in trid:
        chrom = trid.split("-")[0].lower()
    else:
        chrom = trid.lower()

    if chrom in ("chrx", "x"):
        return "chrX"
    if chrom in ("chry", "y"):
        return "chrY"
    if chrom in ("chrm", "m", "mt", "chrmt"):
        return "chrM"
    return "autosome"


def get_motif_size_category(motif):
    """Categorize a motif by its length.

    Args:
        motif: The motif string (e.g. ``"CAG"``, ``"AT"``).

    Returns:
        One of ``"1bp"`` .. ``"6bp"``, ``"7-24bp"`` or ``"25+bp"``.
    """
    size = len(motif)
    if size <= 6:
        return f"{size}bp"
    if size <= 24:
        return "7-24bp"
    return "25+bp"


def find_trios(sample_lookup_or_df):
    """Identify complete trios from the sample set.

    A complete trio is a child sample for which both the maternal_id and the
    paternal_id are non-blank and also present as sample_ids in the sample set.

    Accepts either the ``sample_lookup`` dict (sample_id -> metadata dict) or the
    cleaned ``sample_df`` DataFrame produced by
    ``input_tables.read_sample_metadata``; both expose ``maternal_id`` and
    ``paternal_id`` fields when present.

    Args:
        sample_lookup_or_df: Either a dict mapping sample_id to a metadata dict,
            or a pandas DataFrame with ``sample_id`` plus optional
            ``maternal_id`` / ``paternal_id`` columns.

    Returns:
        A list of ``(child_id, mother_id, father_id)`` tuples, one per complete
        trio, in the order the children appear in the input.

    Raises:
        ValueError: If any child_id would appear in more than one trio.
    """
    rows = _iter_sample_rows(sample_lookup_or_df)

    sample_ids = set(row["sample_id"] for row in rows)
    trios = []
    for row in rows:
        maternal_id = _blank_to_empty(row.get("maternal_id"))
        paternal_id = _blank_to_empty(row.get("paternal_id"))
        if not maternal_id or not paternal_id:
            continue
        if maternal_id in sample_ids and paternal_id in sample_ids:
            trios.append((row["sample_id"], maternal_id, paternal_id))

    duplicate_child_ids = sorted(
        child_id for child_id, count in Counter(trio[0] for trio in trios).items() if count > 1
    )
    if duplicate_child_ids:
        raise ValueError(f"Child ID(s) appear in multiple trios: {set(duplicate_child_ids)}")

    return trios


def _iter_sample_rows(sample_lookup_or_df):
    """Normalize a sample_lookup dict or sample_df into a list of row dicts."""
    if isinstance(sample_lookup_or_df, dict):
        return list(sample_lookup_or_df.values())
    # Assume a pandas DataFrame.
    return sample_lookup_or_df.to_dict(orient="records")


def _blank_to_empty(value):
    """Return a stripped string for ``value``, treating missing/NaN as ``""``."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in ("", "nan", "none", "na", "n/a", "null"):
        return ""
    return text


def _new_trio_stats():
    """Return a fresh per-trio stats accumulator (by-chromosome + by-motif-size)."""
    return {
        "by_chrom": {
            "autosome": {"violations": 0, "total": 0},
            "chrX": {"violations": 0, "total": 0},
            "chrY": {"violations": 0, "total": 0},
            "chrM": {"violations": 0, "total": 0},
        },
        "by_motif_size": {
            category: {"violations": 0, "total": 0} for category in MOTIF_SIZE_CATEGORIES
        },
    }


def _process_trio(child_id, mother_id, father_id, locus_rows, threshold):
    """Count Mendelian violations for a single trio over the in-memory matrix.

    Iterates the child's genotype at every locus, dispatching by chromosome
    category (chrY: child+father; chrM: child+mother; otherwise all three),
    requiring at least two distinct allele sizes across the present members, and
    tallying totals and violations by chromosome, by motif-size category, and by
    canonical 1-6bp motif.

    Args:
        child_id: Sample id of the child.
        mother_id: Sample id of the mother.
        father_id: Sample id of the father.
        locus_rows: The in-memory matrix rows (each with ``trid``, ``motif`` and
            a ``genotypes`` dict mapping sample_id -> cell string).
        threshold: Alleles match if they differ by strictly fewer than this many
            repeats.

    Returns:
        A tuple ``(stats, per_motif_stats)`` where ``stats`` has ``by_chrom`` and
        ``by_motif_size`` keys and ``per_motif_stats`` maps each canonical motif
        to ``{"violations": int, "total": int}``.
    """
    stats = _new_trio_stats()
    per_motif_stats = {motif: {"violations": 0, "total": 0} for motif in ALL_CANONICAL_MOTIFS}

    for locus_row in locus_rows:
        genotypes = locus_row["genotypes"]
        child_genotype = parse_genotype(genotypes.get(child_id))
        if child_genotype is None:
            continue

        motif = locus_row["motif"]
        # Skip loci with a blank motif (len 0 -> get_motif_size_category returns
        # "0bp", which is not a MOTIF_SIZE_CATEGORIES key and would KeyError in
        # _tally) or an N-containing motif.
        if not motif or "N" in motif.upper():
            continue

        chrom_category = get_chrom_category(locus_row["trid"])
        motif_size_category = get_motif_size_category(motif)
        canonical_motif = compute_canonical_motif(motif) if len(motif) <= 6 else None

        if chrom_category == "chrY":
            father_genotype = parse_genotype(genotypes.get(father_id))
            if father_genotype is None:
                continue
            if len(set(child_genotype) | set(father_genotype)) < 2:
                continue
            is_violation = any(
                not allele_matches_any(child_allele, father_genotype, threshold)
                for child_allele in child_genotype
            )
            _tally(stats, per_motif_stats, "chrY", motif_size_category, canonical_motif, is_violation)

        elif chrom_category == "chrM":
            mother_genotype = parse_genotype(genotypes.get(mother_id))
            if mother_genotype is None:
                continue
            if len(set(child_genotype) | set(mother_genotype)) < 2:
                continue
            is_violation = any(
                not allele_matches_any(child_allele, mother_genotype, threshold)
                for child_allele in child_genotype
            )
            _tally(stats, per_motif_stats, "chrM", motif_size_category, canonical_motif, is_violation)

        else:
            mother_genotype = parse_genotype(genotypes.get(mother_id))
            father_genotype = parse_genotype(genotypes.get(father_id))
            if mother_genotype is None or father_genotype is None:
                continue
            if len(set(child_genotype) | set(mother_genotype) | set(father_genotype)) < 2:
                continue
            _tally(
                stats, per_motif_stats, chrom_category, motif_size_category, canonical_motif,
                check_violation(child_genotype, mother_genotype, father_genotype, threshold),
            )

    return stats, per_motif_stats


def _tally(stats, per_motif_stats, chrom_category, motif_size_category, canonical_motif, is_violation):
    """Increment the total (and, if a violation, the violation) counters in place."""
    stats["by_chrom"][chrom_category]["total"] += 1
    stats["by_motif_size"][motif_size_category]["total"] += 1
    if canonical_motif is not None:
        per_motif_stats[canonical_motif]["total"] += 1
    if is_violation:
        stats["by_chrom"][chrom_category]["violations"] += 1
        stats["by_motif_size"][motif_size_category]["violations"] += 1
        if canonical_motif is not None:
            per_motif_stats[canonical_motif]["violations"] += 1


def per_sample_columns():
    """Return the ordered column names of the per-sample Mendelian-violation table.

    The layout matches the original ``write_output`` header: ``sample_id``, then
    ``{cat}_violations`` / ``{cat}_total`` for each of autosome, chrX, chrY,
    chrM, then ``motif_{safecat}_violations`` / ``motif_{safecat}_total`` for
    each motif-size category (with ``-`` -> ``_`` and ``+`` -> ``plus`` so
    ``7-24bp`` -> ``7_24bp`` and ``25+bp`` -> ``25plusbp``), then
    ``total_violations`` and ``total_loci``.

    Returns:
        A list of column-name strings.
    """
    columns = ["sample_id"]
    for chrom_category in ["autosome", "chrX", "chrY", "chrM"]:
        columns.append(f"{chrom_category}_violations")
        columns.append(f"{chrom_category}_total")
    for motif_size_category in MOTIF_SIZE_CATEGORIES:
        safe_category = motif_size_category.replace("-", "_").replace("+", "plus")
        columns.append(f"motif_{safe_category}_violations")
        columns.append(f"motif_{safe_category}_total")
    columns.append("total_violations")
    columns.append("total_loci")
    return columns


def per_motif_columns():
    """Return the ordered column names of the per-motif Mendelian-violation table.

    The layout matches the original ``write_per_motif_output`` header:
    ``sample_id``, then ``mv_{motif}`` / ``total_{motif}`` for each canonical
    1-6bp motif (in sorted order).

    Returns:
        A list of column-name strings.
    """
    columns = ["sample_id"]
    for motif in ALL_CANONICAL_MOTIFS:
        columns.append(f"mv_{motif}")
        columns.append(f"total_{motif}")
    return columns


def compute_mendelian_violations(locus_rows, sample_lookup, sample_df, threshold=2):
    """Compute per-sample and per-motif Mendelian-violation rows for all trios.

    Drives the computation entirely from the in-memory matrix (``locus_rows``)
    and the trio relationships in the sample metadata. Returns ``([], [])`` when
    no complete trio exists, so the caller can skip writing the Mendelian tables.

    Args:
        locus_rows: The in-memory matrix rows from
            ``input_tables.read_repeat_copy_numbers`` (each with ``trid``,
            ``motif`` and a ``genotypes`` dict).
        sample_lookup: The sample_id -> metadata-dict lookup (unused for the
            computation itself, accepted for a uniform call signature alongside
            the other build stages).
        sample_df: The cleaned sample-metadata DataFrame (or the sample_lookup
            dict) used by ``find_trios`` to discover trios.
        threshold: Alleles match if they differ by strictly fewer than this many
            repeats. Defaults to 2. The "match" test is strict
            (``abs(a - b) < threshold``), so a difference of exactly ``threshold``
            counts as a non-match (and can be a violation).

    Returns:
        A tuple ``(per_sample_rows, per_motif_rows)`` of two lists of dicts, one
        dict per trio child, keyed by the column names from
        ``per_sample_columns`` and ``per_motif_columns`` respectively. Both lists
        are empty when there are no trios.
    """
    trios = find_trios(sample_df if sample_df is not None else sample_lookup)
    if not trios:
        return [], []

    per_sample_rows = []
    per_motif_rows = []
    for child_id, mother_id, father_id in trios:
        stats, per_motif_stats = _process_trio(
            child_id, mother_id, father_id, locus_rows, threshold)

        per_sample_row = {"sample_id": child_id}
        total_violations = 0
        total_loci = 0
        for chrom_category in ["autosome", "chrX", "chrY", "chrM"]:
            per_sample_row[f"{chrom_category}_violations"] = stats["by_chrom"][chrom_category]["violations"]
            per_sample_row[f"{chrom_category}_total"] = stats["by_chrom"][chrom_category]["total"]
            total_violations += stats["by_chrom"][chrom_category]["violations"]
            total_loci += stats["by_chrom"][chrom_category]["total"]
        for motif_size_category in MOTIF_SIZE_CATEGORIES:
            safe_category = motif_size_category.replace("-", "_").replace("+", "plus")
            per_sample_row[f"motif_{safe_category}_violations"] = stats["by_motif_size"][motif_size_category]["violations"]
            per_sample_row[f"motif_{safe_category}_total"] = stats["by_motif_size"][motif_size_category]["total"]
        per_sample_row["total_violations"] = total_violations
        per_sample_row["total_loci"] = total_loci
        per_sample_rows.append(per_sample_row)

        per_motif_row = {"sample_id": child_id}
        for motif in ALL_CANONICAL_MOTIFS:
            per_motif_row[f"mv_{motif}"] = per_motif_stats[motif]["violations"]
            per_motif_row[f"total_{motif}"] = per_motif_stats[motif]["total"]
        per_motif_rows.append(per_motif_row)

    return per_sample_rows, per_motif_rows
