"""Trio Mendelian-violation QC computed directly from the in-memory genotype matrix.

The genotypes come straight out of the in-memory repeat-copy-numbers matrix (the
``locus_rows`` produced by ``input_tables.read_repeat_copy_numbers``) and the trio
relationships from the sample-metadata lookup. Each complete trio is compared at
every locus whose motif is made of A, C, G and T only:

  * Each matrix cell is parsed by ``allele_histograms.parse_genotype``. A no-call
    (empty, ``.``, ``./.``) or malformed cell of a member that the comparison needs
    skips the locus for that trio.
  * The per-chromosome dispatch (``_process_trio``): chrY compares the child with the
    father, chrM and a hemizygous (single-allele) chrX child compare the child with
    the mother, and the autosomes and a diploid chrX child compare the child with
    both parents (``check_violation``). A single-allele child call on an autosome may
    have come from either parent; only on chrX must it come from the mother.
  * A locus counts only when the compared members carry at least 2 distinct allele
    sizes. The filter is skipped for a haploid comparison (one allele on each side,
    as on chrY or chrM), where it would mean "the two disagree", would leave only
    discordant loci in the denominator, and would make ``chrY_violations /
    chrY_total`` something other than a violation rate.

Violations and totals are tallied by chromosome, by motif-size category and by
canonical 1-6 bp motif into the per-sample and per-motif rows (``_per_sample_row``
/ ``_per_motif_row``) that result_database writes as the ``mendelian_violations``
and ``mendelian_violations_per_motif`` tables.

All functions are pure -- they take their inputs as arguments and return values
without mutating shared module state.
"""

from allele_histograms import parse_genotype
from analysis_columns import blank_to_none
from motif_utilities import compute_canonical_motif, generate_all_canonical_motifs


# Motif size categories, in display/column order.
MOTIF_SIZE_CATEGORIES = ["1bp", "2bp", "3bp", "4bp", "5bp", "6bp", "7-24bp", "25+bp"]

# All canonical motifs of size 1-6bp, sorted for deterministic column ordering
# in the per-motif output table (generate_all_canonical_motifs returns a set).
ALL_CANONICAL_MOTIFS = sorted(generate_all_canonical_motifs(6))


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
    """Return True if the child's genotype violates Mendelian inheritance from both parents.

    A genotype is consistent with Mendelian inheritance if the child could have
    inherited one allele from the mother and one from the father. Alleles are
    considered matching if they differ by strictly fewer than ``threshold``
    repeats. A hemizygous father (on chrX, for a daughter) passes on his single
    allele.

    This is the comparison for the autosomes and a diploid chrX child. The loci
    inherited from one parent (chrY, chrM and a hemizygous son's chrX) are
    compared by ``_process_trio`` with that parent alone and never reach this
    function, so a single-allele child call here is just a one-allele call, which
    either parent could have supplied.

    Args:
        child_genotype: Tuple of alleles for the child (1 or 2 elements).
        mother_genotype: Tuple of alleles for the mother (1 or 2 elements).
        father_genotype: Tuple of alleles for the father (1 or 2 elements).
        threshold: Alleles match if they differ by strictly fewer than this many
            repeats.

    Returns:
        True if the genotype is a Mendelian violation, False otherwise.
    """
    # Case 1: a one-allele call, which either parent could have supplied.
    if len(child_genotype) == 1:
        return not (allele_matches_any(child_genotype[0], mother_genotype, threshold)
                    or allele_matches_any(child_genotype[0], father_genotype, threshold))

    # Case 2: the child is diploid. One allele must come from the mother and the
    # other from the father, in either assignment. allele_matches_any compares
    # against every allele a parent has, so the same rule covers a hemizygous
    # parent (e.g. the father on chrX for a daughter) and a diploid one.
    first_child_allele, second_child_allele = child_genotype
    return not (
        (allele_matches_any(first_child_allele, mother_genotype, threshold)
         and allele_matches_any(second_child_allele, father_genotype, threshold))
        or
        (allele_matches_any(first_child_allele, father_genotype, threshold)
         and allele_matches_any(second_child_allele, mother_genotype, threshold))
    )


def get_chrom_category(trid):
    """Categorize a locus id by chromosome type.

    The chromosome is the first ``-``-delimited field of the trid, as in the
    build's ``Chrom``: ``locus_annotations.parse_locus_id`` has already required
    the ``chrom-start-end-motif`` form.

    Args:
        trid: Locus id of the form ``chrom-start-end-motif`` (e.g. ``chrX-...``).

    Returns:
        One of ``"autosome"``, ``"chrX"``, ``"chrY"`` or ``"chrM"``.
    """
    chrom = trid.split("-")[0].lower()

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


def find_trios(sample_lookup):
    """Identify complete trios from the sample metadata.

    A complete trio is a child sample for which both the maternal_id and the
    paternal_id are non-blank and also present as sample_ids in the sample set.
    ``sample_lookup`` holds one entry per sample_id, so a child appears in at most
    one trio.

    Args:
        sample_lookup: Mapping of sample_id -> metadata dict (as returned by
            ``input_tables.read_sample_metadata``), with optional ``maternal_id`` /
            ``paternal_id`` fields.

    Returns:
        A list of ``(child_id, mother_id, father_id)`` tuples, one per complete
        trio, in the order the children appear in ``sample_lookup``.
    """
    trios = []
    for sample_id, row in sample_lookup.items():
        maternal_id = blank_to_none(row.get("maternal_id"))
        paternal_id = blank_to_none(row.get("paternal_id"))
        if maternal_id is None or paternal_id is None:
            continue
        maternal_id, paternal_id = str(maternal_id).strip(), str(paternal_id).strip()
        if maternal_id in sample_lookup and paternal_id in sample_lookup:
            trios.append((sample_id, maternal_id, paternal_id))
    return trios


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


def _keep_genotyped_trios(trios, locus_rows):
    """Drop trios whose members are not genotype columns of the matrix.

    The sample metadata can describe families whose samples were never
    genotyped. Such a trio is not comparable at any locus, so it would yield a
    per-sample row of zeros that reads as a clean trio. Only column presence is
    checked: a trio whose three members are genotype columns is kept even when
    none of its loci can be compared (e.g. every call is a no-call), and its row
    then has ``total_loci == 0``, which tells it apart from a clean trio (zero
    violations over a non-zero ``total_loci``).

    Args:
        trios: List of ``(child_id, mother_id, father_id)`` tuples.
        locus_rows: The in-memory matrix rows, each with a ``genotypes`` dict.
            ``input_tables.read_repeat_copy_numbers`` keys every row's dict by all
            of the matrix's sample columns, so the first row names them all.

    Returns:
        A list of the trios all three of whose members appear as genotype
        columns, in the input order.
    """
    genotyped = set(locus_rows[0]["genotypes"]) if locus_rows else set()

    kept = []
    for trio in trios:
        missing = [sample_id for sample_id in trio if sample_id not in genotyped]
        if missing:
            print(f"WARNING: skipping trio {trio} for Mendelian QC: "
                  f"not genotyped in the matrix: {missing}")
            continue
        kept.append(trio)
    return kept


def _prepare_loci(locus_rows):
    """Derive the per-locus values the trio loop needs, once for all trios.

    The chromosome category, the motif-size category and the canonical motif depend only on
    the locus, so they are computed here rather than inside ``_process_trio``, which runs
    once per trio over the same rows.

    Args:
        locus_rows: The in-memory matrix rows (each with ``trid``, ``motif`` and a
            ``genotypes`` dict mapping sample_id -> cell string).

    Returns:
        A list of ``(genotypes, chrom_category, motif_size_category, canonical_motif)``
        tuples, one per scorable locus, in the input order.
    """
    prepared = []
    for locus_row in locus_rows:
        motif = locus_row["motif"]
        # Skip loci with a blank motif (len 0 -> get_motif_size_category returns
        # "0bp", which is not a MOTIF_SIZE_CATEGORIES key and would KeyError in
        # _tally) or a motif carrying any base outside ACGT. compute_canonical_motif
        # accepts the whole IUPAC alphabet, so an ambiguity code survives
        # canonicalization, but per_motif_stats is keyed only by the ACGT canonical
        # motifs, so "N", "R", "Y" and friends would all KeyError in _tally.
        if not motif or not set(motif.upper()) <= set("ACGT"):
            continue
        prepared.append((
            locus_row["genotypes"],
            get_chrom_category(locus_row["trid"]),
            get_motif_size_category(motif),
            compute_canonical_motif(motif) if len(motif) <= 6 else None,
        ))
    return prepared


def _process_trio(child_id, mother_id, father_id, prepared_loci, threshold):
    """Count Mendelian violations for a single trio over the in-memory matrix.

    Iterates the child's genotype at every locus, dispatching by chromosome
    category (chrY: child+father; chrM and a hemizygous son's chrX: child+mother;
    otherwise all three), requiring at least two distinct allele sizes across the
    members that count unless both sides of a single-parent comparison are haploid
    (chrY, chrM), and tallying totals and violations by chromosome, by motif-size
    category, and by canonical 1-6bp motif.

    Args:
        child_id: Sample id of the child.
        mother_id: Sample id of the mother.
        father_id: Sample id of the father.
        prepared_loci: The ``_prepare_loci`` output for the matrix rows: one
            ``(genotypes, chrom_category, motif_size_category, canonical_motif)``
            tuple per scorable locus.
        threshold: Alleles match if they differ by strictly fewer than this many
            repeats.

    Returns:
        A tuple ``(stats, per_motif_stats)`` where ``stats`` has ``by_chrom`` and
        ``by_motif_size`` keys and ``per_motif_stats`` maps each canonical motif
        to ``{"violations": int, "total": int}``.
    """
    stats = _new_trio_stats()
    per_motif_stats = {motif: {"violations": 0, "total": 0} for motif in ALL_CANONICAL_MOTIFS}

    for genotypes, chrom_category, motif_size_category, canonical_motif in prepared_loci:
        child_genotype = parse_genotype(genotypes.get(child_id))
        if child_genotype is None:
            continue

        # chrY comes from the father, chrM from the mother, and a hemizygous son's single X allele
        # from the mother, so those loci are compared with that one parent; the other parent is
        # not consulted and need not have a call.
        if chrom_category == "chrY":
            single_parent_id = father_id
        elif chrom_category == "chrM" or (chrom_category == "chrX" and len(child_genotype) == 1):
            single_parent_id = mother_id
        else:
            single_parent_id = None

        if single_parent_id is not None:
            parent_genotype = parse_genotype(genotypes.get(single_parent_id))
            if parent_genotype is None:
                continue
            # The "at least 2 distinct alleles" variability filter degenerates on a haploid
            # comparison: with one allele on each side, "2 distinct" is exactly "parent and child
            # disagree", so applying it would leave only discordant loci in the denominator and
            # violations/total would stop being a violation rate. Keep it only when at least one
            # side is diploid, where it still means what it does on the autosomes.
            if (len(child_genotype) > 1 or len(parent_genotype) > 1) \
                    and len(set(child_genotype) | set(parent_genotype)) < 2:
                continue
            is_violation = any(
                not allele_matches_any(child_allele, parent_genotype, threshold)
                for child_allele in child_genotype
            )
        else:
            mother_genotype = parse_genotype(genotypes.get(mother_id))
            father_genotype = parse_genotype(genotypes.get(father_id))
            if mother_genotype is None or father_genotype is None:
                continue
            if len(set(child_genotype) | set(mother_genotype) | set(father_genotype)) < 2:
                continue
            is_violation = check_violation(child_genotype, mother_genotype, father_genotype, threshold)

        _tally(stats, per_motif_stats, chrom_category, motif_size_category, canonical_motif, is_violation)

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


def _per_sample_row(child_id, stats):
    """Return one trio child's row of the per-sample Mendelian-violation table.

    The keys, in column order, are ``sample_id``, then ``{cat}_violations`` /
    ``{cat}_total`` for each of autosome, chrX, chrY, chrM, then
    ``motif_{safecat}_violations`` / ``motif_{safecat}_total`` for each
    motif-size category (with ``-`` -> ``_`` and ``+`` -> ``plus`` so
    ``7-24bp`` -> ``7_24bp`` and ``25+bp`` -> ``25plusbp``), then
    ``total_violations`` and ``total_loci``. result_database.write_mendelian_tables
    takes the table's column order from these keys.

    Args:
        child_id: The trio child's sample id.
        stats: The by-chromosome / by-motif-size accumulator ``_process_trio`` returns.

    Returns:
        dict: The row, keyed in column order.
    """
    row = {"sample_id": child_id}
    for chrom_category, counts in stats["by_chrom"].items():
        row[f"{chrom_category}_violations"] = counts["violations"]
        row[f"{chrom_category}_total"] = counts["total"]
    for motif_size_category, counts in stats["by_motif_size"].items():
        safe_category = motif_size_category.replace("-", "_").replace("+", "plus")
        row[f"motif_{safe_category}_violations"] = counts["violations"]
        row[f"motif_{safe_category}_total"] = counts["total"]
    row["total_violations"] = sum(counts["violations"] for counts in stats["by_chrom"].values())
    row["total_loci"] = sum(counts["total"] for counts in stats["by_chrom"].values())
    return row


def _per_motif_row(child_id, per_motif_stats):
    """Return one trio child's row of the per-motif Mendelian-violation table.

    The keys, in column order, are ``sample_id``, then ``mv_{motif}`` /
    ``total_{motif}`` for each canonical 1-6bp motif in ALL_CANONICAL_MOTIFS
    (sorted) order.

    Args:
        child_id: The trio child's sample id.
        per_motif_stats: The canonical motif -> counts mapping ``_process_trio`` returns.

    Returns:
        dict: The row, keyed in column order.
    """
    row = {"sample_id": child_id}
    for motif in ALL_CANONICAL_MOTIFS:
        row[f"mv_{motif}"] = per_motif_stats[motif]["violations"]
        row[f"total_{motif}"] = per_motif_stats[motif]["total"]
    return row


def compute_mendelian_violations(locus_rows, sample_lookup, threshold=2):
    """Compute per-sample and per-motif Mendelian-violation rows for all trios.

    Drives the computation entirely from the in-memory matrix (``locus_rows``)
    and the trio relationships in the sample metadata. Trios whose child or
    parents are not genotype columns of the matrix are dropped (the metadata can
    describe samples the matrix does not carry). A trio whose members are all
    genotype columns is kept even when none of its loci could be compared; its
    row then has ``total_loci == 0`` (see ``_keep_genotyped_trios``). Returns
    ``([], [])`` when no genotyped complete trio exists, so the caller can skip
    writing the Mendelian tables.

    Args:
        locus_rows: The in-memory matrix rows from
            ``input_tables.read_repeat_copy_numbers`` (each with ``trid``,
            ``motif`` and a ``genotypes`` dict).
        sample_lookup: The sample_id -> metadata-dict lookup from
            ``input_tables.read_sample_metadata``, from which ``find_trios`` reads
            the trios.
        threshold: Alleles match if they differ by strictly fewer than this many
            repeats. Defaults to 2. The "match" test is strict
            (``abs(a - b) < threshold``), so a difference of exactly ``threshold``
            counts as a non-match (and can be a violation).

    Returns:
        A tuple ``(per_sample_rows, per_motif_rows)`` of two lists of dicts, one
        dict per trio child, laid out by ``_per_sample_row`` and ``_per_motif_row``
        respectively. Both lists are empty when there are no trios.
    """
    trios = _keep_genotyped_trios(find_trios(sample_lookup), locus_rows)
    if not trios:
        return [], []

    # Derived once here rather than per trio: none of these values depends on the trio.
    prepared_loci = _prepare_loci(locus_rows)

    per_sample_rows = []
    per_motif_rows = []
    for child_id, mother_id, father_id in trios:
        stats, per_motif_stats = _process_trio(
            child_id, mother_id, father_id, prepared_loci, threshold)
        per_sample_rows.append(_per_sample_row(child_id, stats))
        per_motif_rows.append(_per_motif_row(child_id, per_motif_stats))

    return per_sample_rows, per_motif_rows
