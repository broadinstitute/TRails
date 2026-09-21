"""Locus-level annotation helpers for the TRails build pipeline.

This module owns the two annotation stages that decorate each per-locus record
after the allele histograms have been built:

  * Stage 6 -- derived coordinate / motif columns parsed straight out of the
    ``LocusId`` (``add_derived_locus_columns``), plus the gene-region rank
    lookup.
  * Stage 7 -- known-disease-locus matching against a cached class-A variant
    catalog (and, optionally, a cached STRchive locus set), the canonical-motif
    column, and the ``IsKnownMotif`` / ``IsInMendelianGene`` flags.

It is a faithful, standalone port of the relevant pieces of the original
``tr_analysis_helpers.py`` (``compute_jaccard``, ``motifs_match``,
``load_known_disease_loci``, ``matches_disease_locus``) and the
``add_derived_locus_columns`` / known-motif logic from ``analyze_results.py``.

Design choices (per the TRails port blueprint):
  * No dependency on ``str_analysis``; motif primitives come from
    ``motif_utilities``.
  * No network access by default. ``load_known_disease_loci`` reads the cached
    class-A catalog JSON, and only builds STRchive interval trees when handed a
    cached STRchive-loci JSON path (``fetch_strchive`` triggers the optional
    network fetch, kept off by default).
  * The internal representation is a list of dicts (one dict per locus); these
    functions mutate each dict in place to add new columns.

All functions are pure with respect to global state (they only read their
arguments and, for the ``add_*`` helpers, mutate the passed-in records).
"""

import collections
import json

import intervaltree

# Result of load_known_disease_loci. interval_trees / strchive_trees are
# chrom -> IntervalTree (disease loci, payload in .data). locus_lookup maps a
# disease locus's LocusId and a "chrom-start0-end1-RepeatUnit" coordinate key to
# its catalog dict; it is built only when build_locus_lookup=True (the serving
# path needs it for the locus-detail page; the build path leaves it empty).
KnownDiseaseLoci = collections.namedtuple(
    "KnownDiseaseLoci", ["interval_trees", "strchive_trees", "locus_lookup"])

from motif_utilities import COMPLEMENT, compute_canonical_motif


# STRchive raw-data URL, used only when the optional network fetch is requested.
STRCHIVE_URL = (
    "https://raw.githubusercontent.com/dashnowlab/STRchive/"
    "refs/heads/main/data/STRchive-loci.json"
)


# Gene-region prioritization rank (lower number = more important / more likely
# functional). Matches GENE_REGION_RANK in the original analyze_results.py.
GENE_REGION_RANK = {
    "CDS": 1,
    "5' UTR": 2,
    "3' UTR": 3,
    "promoter": 4,
    "intron": 5,
    "exon": 6,
    "intergenic": 7,
}


def gene_region_rank(gene_region):
    """Return the numeric rank for a gene region, or None if unranked/missing.

    Args:
        gene_region: A gene-region string such as ``"CDS"`` or ``"intron"``,
            or None when the locus has no gene-region annotation.

    Returns:
        The integer rank from ``GENE_REGION_RANK``, or None if ``gene_region``
        is missing or not one of the known regions.
    """
    if gene_region is None:
        return None
    return GENE_REGION_RANK.get(gene_region)


def compute_jaccard(start1, end1, start2, end2):
    """Return the Jaccard index of two half-open intervals.

    Intervals are 0-based start, 1-based (exclusive) end. The Jaccard index is
    ``overlap / union``; it is 0.0 when the intervals do not overlap.

    Args:
        start1: 0-based start of the first interval.
        end1: 1-based (exclusive) end of the first interval.
        start2: 0-based start of the second interval.
        end2: 1-based (exclusive) end of the second interval.

    Returns:
        A float in ``[0.0, 1.0]``: 1.0 for identical intervals, 0.0 for
        disjoint ones, and the overlap-over-union fraction otherwise.
    """
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start >= overlap_end:
        return 0.0
    overlap_size = overlap_end - overlap_start
    union_size = (end1 - start1) + (end2 - start2) - overlap_size
    return overlap_size / union_size if union_size > 0 else 0.0


def motifs_match(motif1, motif2):
    """Return whether two motifs match under TRExplorer-style rules.

    For motifs of at most 6 bp, two motifs match when their canonical forms
    (including reverse complement) are equal. For longer motifs -- which are not
    meaningfully canonicalizable here -- two motifs match when they have the
    same length.

    Args:
        motif1: The first motif string (e.g. the query locus motif).
        motif2: The second motif string (e.g. a disease-locus repeat unit).

    Returns:
        True if the two motifs match under the size-dependent rule, else False.
        Returns False if either motif is empty/None.
    """
    if not motif1 or not motif2:
        return False
    if len(motif1) <= 6:
        return (compute_canonical_motif(motif1, include_reverse_complement=True)
                == compute_canonical_motif(motif2, include_reverse_complement=True))
    return len(motif1) == len(motif2)


def fetch_strchive_loci():
    """Fetch STRchive loci over the network and build chrom -> IntervalTree.

    This performs a live HTTP request and is therefore opt-in only (the build
    pipeline keeps it off by default for hermetic, reproducible runs). The
    ``requests`` import is local so the module has no hard network dependency.

    Returns:
        A dict mapping chromosome (without ``chr`` prefix) to an
        ``intervaltree.IntervalTree`` whose intervals carry the locus dict in
        their ``.data`` attribute.
    """
    import requests

    print(f"Fetching STRchive loci from {STRCHIVE_URL}")
    response = requests.get(STRCHIVE_URL)
    response.raise_for_status()
    return _build_strchive_trees(response.json())


def _build_strchive_trees(loci):
    """Build chrom -> IntervalTree from a list of STRchive locus dicts.

    Only loci with a disease, chromosome, and both hg38 coordinates are kept.
    The hg38 start is converted from 1-based to 0-based (``start_hg38 - 1``).

    Args:
        loci: A list of STRchive locus dicts (as parsed from STRchive-loci.json).

    Returns:
        A dict mapping chromosome (without ``chr`` prefix) to an IntervalTree of
        the qualifying loci.
    """
    interval_trees = collections.defaultdict(intervaltree.IntervalTree)
    count = 0
    for locus in loci:
        if not locus.get("disease") or not locus.get("chrom") or \
           not locus.get("start_hg38") or not locus.get("stop_hg38"):
            continue
        interval_trees[locus["chrom"].replace("chr", "")].addi(
            int(locus["start_hg38"]) - 1, int(locus["stop_hg38"]), data=locus)
        count += 1

    print(f"Loaded {count:,d} STRchive disease loci")
    return dict(interval_trees)


# Cell values that mean "nothing was recorded" in an affected_status / analysis_status column
# (matches input_tables._is_blank).
BLANK_STATUS_VALUES = {"", "nan", "none", "na", "n/a", "null"}


def normalize_affected_status_for_logic(value):
    """Normalize an affected-status cell so comparisons don't depend on its spelling.

    Both the build (input_tables.read_sample_metadata) and the server compare affected statuses,
    and they must agree on what "Possibly Affected" means, so this is the one definition.

    Args:
        value: An affected-status cell, which may be None, NaN or any of the blank spellings in
            BLANK_STATUS_VALUES.

    Returns:
        None for a blank value, "affected" for "possibly affected", and the lower-cased stripped
        value otherwise.
    """
    if value is None:
        return None
    if isinstance(value, float) and value != value:  # NaN
        return None
    normalized = str(value).strip().lower()
    if normalized in BLANK_STATUS_VALUES:
        return None
    return "affected" if normalized == "possibly affected" else normalized


def load_known_disease_loci(filepath=None, fetch_strchive=False, strchive_filepath=None,
                            build_locus_lookup=False):
    """Load known disease loci into chrom -> IntervalTree(s) for matching.

    Reads the cached class-A variant catalog JSON (e.g.
    ``variant_catalog_without_offtargets.GRCh38.json``) and keeps only loci that
    carry at least one named disease (``any(d["Name"] for d in Diseases)``) and a
    parseable ``MainReferenceRegion`` (``chrom:start-end``). The 1-based catalog
    start is converted to a 0-based start (``start - 1``); the end is already
    1-based. Each interval carries its full locus dict in ``.data``.

    The optional STRchive trees are built only when an explicit cached
    ``strchive_filepath`` is given, or when ``fetch_strchive`` is True (which
    triggers a live network fetch). Both default off, so the build does no
    network access by default.

    Args:
        filepath: Path to the variant catalog JSON file, or None to load only
            the STRchive fallback (``interval_trees`` is then empty).
        fetch_strchive: When True, fetch STRchive loci over the network. Kept
            False by default for hermetic builds.
        strchive_filepath: Optional path to a cached STRchive-loci.json file. If
            given, STRchive trees are built from it (no network access). Takes
            precedence over ``fetch_strchive``.
        build_locus_lookup: When True, also build ``locus_lookup`` (keyed by each
            locus's ``LocusId`` and by a ``"chrom-start0-end1-RepeatUnit"``
            coordinate key) for the serving path's locus-detail lookups. The
            build pipeline does not need it and leaves it ``{}`` (default).

    Returns:
        A ``KnownDiseaseLoci`` namedtuple ``(interval_trees, strchive_trees,
        locus_lookup)``. ``interval_trees`` and ``strchive_trees`` map a
        chromosome (without ``chr`` prefix) to an ``intervaltree.IntervalTree``;
        ``strchive_trees`` is ``{}`` when no STRchive source was provided and
        ``locus_lookup`` is ``{}`` unless ``build_locus_lookup`` was True.
    """
    loci = []
    if filepath is not None:
        with open(filepath, "rt") as input_file:
            loci = json.load(input_file)

    interval_trees = collections.defaultdict(intervaltree.IntervalTree)
    locus_lookup = {}
    count = 0
    for locus_data in loci:
        if build_locus_lookup and "LocusId" in locus_data:
            locus_lookup[locus_data["LocusId"]] = locus_data
        if not any(disease.get("Name") for disease in locus_data.get("Diseases", [])):
            continue
        reference_region = locus_data.get("MainReferenceRegion", "")
        if ":" not in reference_region or "-" not in reference_region:
            continue
        chrom_part, coords = reference_region.split(":")
        start_str, end_str = coords.split("-")
        chrom = chrom_part.replace("chr", "")
        start_0based = int(start_str) - 1
        end_1based = int(end_str)
        interval_trees[chrom].addi(start_0based, end_1based, data=locus_data)
        if build_locus_lookup:
            locus_lookup[f"{chrom}-{start_0based}-{end_1based}-{locus_data.get('RepeatUnit', '')}"] = locus_data
        count += 1

    print(f"Loaded {count:,d} variant catalog disease loci")

    strchive_trees = {}
    if strchive_filepath is not None:
        with open(strchive_filepath, "rt") as input_file:
            strchive_trees = _build_strchive_trees(json.load(input_file))
    elif fetch_strchive:
        strchive_trees = fetch_strchive_loci()

    return KnownDiseaseLoci(dict(interval_trees), strchive_trees, locus_lookup)


def strchive_locus_motifs(locus_data):
    """Return the motifs to match a STRchive locus against.

    A STRchive locus lists the reference motif and, separately, the motif(s) seen in
    the pathogenic allele, which for several loci differ (e.g. an interrupted reference
    unit expanding as a pure one). Matching only the reference motif would miss exactly
    the expansions this tool exists to surface, so both lists are returned. This is the
    one definition of that list: the build and the server's locus-detail page both call
    it, and they must agree on which loci are known disease loci.

    Args:
        locus_data: A STRchive locus dict (an interval's ``.data``).

    Returns:
        A list of motif strings: the reference motifs followed by the pathogenic ones.
    """
    return (list(locus_data.get("reference_motif_reference_orientation") or [])
            + list(locus_data.get("pathogenic_motif_reference_orientation") or []))


def matches_disease_locus(locus_id, interval_trees, strchive_trees=None):
    """Return the disease LocusId matched by a query locus, or None.

    Matching uses TRExplorer-style logic: a query locus matches a disease locus
    when they overlap, their interval Jaccard is greater than 0.66, and at least
    one of the disease locus's motifs (its ``RepeatUnit`` plus any
    ``PathogenicMotifs``) matches the query motif under ``motifs_match``. The
    variant catalog is consulted first; the optional STRchive trees are a
    fallback consulted only when the catalog yields no match.

    Args:
        locus_id: A ``"chrom-start-end-motif"`` string. ``chrom`` may carry a
            ``chr`` prefix; ``start`` is 0-based and ``end`` is 1-based.
        interval_trees: chrom -> IntervalTree from ``load_known_disease_loci``.
        strchive_trees: Optional chrom -> IntervalTree STRchive fallback (also
            from ``load_known_disease_loci``). When falsy, no fallback is used.

    Returns:
        The matching disease locus identifier (the catalog ``LocusId`` or, for
        STRchive, its ``locus_id``/``id``), or None if nothing matches.
    """
    chrom, start, end, motif = locus_id.split("-")
    start = int(start)
    end = int(end)
    chrom = chrom.replace("chr", "")

    for interval in interval_trees.get(chrom, intervaltree.IntervalTree()).overlap(start, end):
        if compute_jaccard(start, end, interval.begin, interval.end) <= 0.66:
            continue
        locus_data = interval.data
        for disease_motif in [locus_data["RepeatUnit"]] + (locus_data.get("PathogenicMotifs") or []):
            if motifs_match(motif, disease_motif):
                return locus_data["LocusId"]

    if strchive_trees:
        for interval in strchive_trees.get(chrom, intervaltree.IntervalTree()).overlap(start, end):
            if compute_jaccard(start, end, interval.begin, interval.end) <= 0.66:
                continue
            locus_data = interval.data
            for strchive_motif in strchive_locus_motifs(locus_data):
                if motifs_match(motif, strchive_motif):
                    return locus_data.get("locus_id", locus_data.get("id"))

    return None


def collect_known_disease_canonical_motifs(interval_trees, strchive_trees=None):
    """Return the set of canonical motifs associated with known disease loci.

    Walks every disease-locus interval and collects the canonical form
    (including reverse complement) of each locus's ``RepeatUnit`` and
    ``PathogenicMotifs``. When ``strchive_trees`` is supplied, each STRchive
    locus's motifs (``strchive_locus_motifs``) are collected too, so a
    locus that is recognized as a known disease locus only via the STRchive
    fallback (see ``matches_disease_locus``) still has its motif counted as
    known -- otherwise such a locus would get ``KnownDiseaseLocus`` set but
    ``IsKnownMotif=0``. Motifs containing ``N`` are intentionally skipped --
    they cannot be canonicalized into a single concrete sequence to compare
    against per-locus canonical motifs, so excluding them is by design (do not
    "fix" this).

    Args:
        interval_trees: chrom -> IntervalTree from ``load_known_disease_loci``.
        strchive_trees: chrom -> IntervalTree of STRchive loci from
            ``load_known_disease_loci`` (the fallback catalog), or None.

    Returns:
        A set of canonical motif strings.
    """
    known_canonical_motifs = set()
    for tree in interval_trees.values():
        for interval in tree:
            if not interval.data:
                continue
            for motif in [interval.data["RepeatUnit"]] + (interval.data.get("PathogenicMotifs") or []):
                if "N" in motif:
                    continue
                known_canonical_motifs.add(
                    compute_canonical_motif(motif, include_reverse_complement=True))
    for tree in (strchive_trees or {}).values():
        for interval in tree:
            if not interval.data:
                continue
            for motif in strchive_locus_motifs(interval.data):
                if "N" in motif:
                    continue
                known_canonical_motifs.add(
                    compute_canonical_motif(motif, include_reverse_complement=True))
    return known_canonical_motifs


def parse_locus_id(locus_id, repeat_copy_numbers_tsv=None):
    """Parse a ``chrom-start-end-motif`` locus id into its four fields.

    The locus id is the ``trid`` cell of the repeat-copy-numbers matrix, passed
    through verbatim by ``input_tables.read_repeat_copy_numbers``, so a matrix
    written by a tool that uses some other locus-id shape reaches this function
    unchanged. Rather than letting the tuple unpacking or ``int()`` raise a bare
    message, this names the column, the offending value and the file, in the
    same shape as the sample-id check in ``build_database.build``.

    Args:
        locus_id: The ``LocusId`` / ``trid`` value, e.g. ``"chr1-100-110-CAG"``.
        repeat_copy_numbers_tsv: Path of the matrix the value came from, used in
            the error message. None when the caller does not know it.

    Returns:
        A ``(chrom_field, start0_based, end1_based, motif_field)`` tuple, with
        the two coordinates as ints and the chromosome field exactly as written.

    Raises:
        ValueError: If the locus id does not have exactly four ``-``-separated
            fields, or its start/end fields are not integers.
    """
    source = repeat_copy_numbers_tsv or "the repeat-copy-numbers matrix"
    fields = locus_id.split("-")
    if len(fields) != 4:
        raise ValueError(
            f"The 'trid' (locus id) column in {source} must hold values of the form "
            f"'chrom-start-end-motif' (4 '-'-separated fields); got {locus_id!r} with "
            f"{len(fields)}")
    chrom_field, start_field, end_field, motif_field = fields
    try:
        start0_based = int(start_field)
        end1_based = int(end_field)
    except ValueError:
        raise ValueError(
            f"The 'trid' (locus id) column in {source} must hold values of the form "
            f"'chrom-start-end-motif' with integer start and end coordinates; got "
            f"{locus_id!r}") from None
    return chrom_field, start0_based, end1_based, motif_field


def validate_motif(motif, locus_id, repeat_copy_numbers_tsv=None):
    """Raise a descriptive ValueError if a motif is outside the IUPAC alphabet.

    ``compute_canonical_motif`` reverse-complements the motif through
    ``motif_utilities.COMPLEMENT``, which raises a bare ``KeyError`` naming only
    the offending character for anything outside the table (a decorated motif
    such as ``"(CAG)n"``, a whitespace-padded cell, an ``"N/A"`` placeholder).
    The build fails on such a motif rather than skipping the locus: the locus
    comes from a row the user put in the matrix, so dropping it would silently
    shrink the database, and ``CanonicalMotif`` is one of the columns the server
    groups and searches by, so the locus cannot be stored without it. (The
    Mendelian-QC stage skips rather than fails on the same class of motif
    because it computes a summary statistic over a subset of loci, where a
    skipped locus costs a little precision and nothing else.)

    Args:
        motif: The ``Motif`` value from the record, possibly empty.
        locus_id: The record's ``LocusId``, named in the error message.
        repeat_copy_numbers_tsv: Path of the matrix the value came from, used in
            the error message. None when the caller does not know it.

    Raises:
        ValueError: If the motif carries a character outside ``COMPLEMENT``.
    """
    # An empty motif cell is tolerated (see NumRepeatsInReference below), so only a
    # non-empty motif with a character outside the table is an error.
    if not motif or set(motif.upper()) <= set(COMPLEMENT):
        return
    raise ValueError(
        f"The 'motif' column in {repeat_copy_numbers_tsv or 'the repeat-copy-numbers matrix'} "
        f"must contain only IUPAC bases ({''.join(sorted(COMPLEMENT))}); locus {locus_id} has "
        f"motif {motif!r}")


def add_derived_locus_columns(locus_records, source_label="", repeat_copy_numbers_tsv=None):
    """Add coordinate, motif, source, gene-region, and canonical-motif columns.

    For each per-locus dict, parses the ``LocusId`` (``chrom-start-end-motif``)
    and adds:
      * ``Chrom`` -- chromosome with a ``chr`` prefix.
      * ``Start0Based`` / ``End1Based`` -- the parsed coordinates.
      * ``ReferenceRegion`` -- ``"chrom:start1based-end1based"`` (without the
        ``chr`` prefix), unless already present in the record.
      * ``MotifSize`` -- length of ``Motif``.
      * ``NumRepeatsInReference`` -- ``(End1Based - Start0Based) // len(Motif)``,
        unless already present in the record.
      * ``Source`` -- the supplied ``source_label`` (default ``""``).
      * ``CanonicalMotif`` -- canonical form of ``Motif`` (with reverse
        complement).
      * ``gene_id`` / ``gene_region`` -- pass-through from the record (also
        accepting the ``GencodeGeneId`` / ``GencodeGeneRegion`` aliases); None
        when absent.
      * ``gene_region_rank`` -- always the rank derived from ``gene_region``
        (None if the region is unranked or absent). A rank supplied by the input
        matrix is ignored, so the stored rank can never disagree with the
        region the locus actually has.

    Prints a one-line summary when the input matrix supplied any
    ``gene_region_rank`` values, so their replacement is not silent.

    Args:
        locus_records: A list of per-locus dicts, each with at least ``LocusId``
            and ``Motif`` keys. Mutated in place.
        source_label: The value to store in the ``Source`` column. Defaults to
            the empty string.
        repeat_copy_numbers_tsv: Path of the matrix the records came from, named
            in the error messages raised for a malformed ``LocusId`` or
            ``Motif``. None when the caller does not know it.

    Returns:
        The same ``locus_records`` list (mutated in place), for convenience.

    Raises:
        ValueError: If a ``LocusId`` is not ``chrom-start-end-motif`` with
            integer coordinates, or a ``Motif`` carries a character outside the
            IUPAC alphabet.
    """
    supplied_rank_count = 0

    for record in locus_records:
        chrom_field, start0_based, end1_based, _motif_field = parse_locus_id(
            record["LocusId"], repeat_copy_numbers_tsv)
        record["Chrom"] = f"chr{chrom_field.replace('chr', '')}"
        record["Start0Based"] = start0_based
        record["End1Based"] = end1_based

        if record.get("ReferenceRegion") is None:
            record["ReferenceRegion"] = (
                f"{record['Chrom'].replace('chr', '')}:"
                f"{record['Start0Based'] + 1}-{record['End1Based']}")

        record["MotifSize"] = len(record["Motif"])

        if record.get("NumRepeatsInReference") is None:
            # An empty motif (a malformed/blank cell in the required motif column)
            # has no repeat length; leave NumRepeatsInReference NULL rather than
            # dividing by zero and aborting the whole build.
            record["NumRepeatsInReference"] = (
                (record["End1Based"] - record["Start0Based"]) // record["MotifSize"]
                if record["MotifSize"] else None)

        # A Source column in the input matrix is a recognized pass-through annotation
        # (input_tables.ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED), so keep what it supplied and
        # fall back to the build-wide --source-label, matching the two columns above.
        if record.get("Source") is None:
            record["Source"] = source_label

        validate_motif(record["Motif"], record["LocusId"], repeat_copy_numbers_tsv)
        record["CanonicalMotif"] = compute_canonical_motif(
            record["Motif"], include_reverse_complement=True)

        if record.get("gene_id") is None:
            record["gene_id"] = record.get("GencodeGeneId")
        if record.get("gene_region") is None:
            record["gene_region"] = record.get("GencodeGeneRegion")

        # gene_region_rank is a recognized annotation column of the input matrix
        # (input_tables.ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED) so that it is not mistaken
        # for a sample column, but it is not a pass-through: unlike gene_id / gene_region /
        # Source, it carries no information of its own, being purely GENE_REGION_RANK applied
        # to gene_region. Deriving it unconditionally is what keeps it consistent with the
        # region stored next to it, which matters because the server's "region" ordering sorts
        # by this column while displaying gene_region (results_server SORT_MAPPING).
        if record.get("gene_region_rank") is not None:
            supplied_rank_count += 1
        record["gene_region_rank"] = gene_region_rank(record.get("gene_region"))

    if supplied_rank_count:
        print(f"  ignored the gene_region_rank supplied by the input matrix for "
              f"{supplied_rank_count:,d} of {len(locus_records):,d} loci; the rank is always "
              f"derived from gene_region")

    return locus_records


def compute_is_known_motif(canonical_motif, known_canonical_motifs):
    """Return 1 if a canonical motif is a known disease motif, else 0.

    Args:
        canonical_motif: A locus's canonical motif string.
        known_canonical_motifs: The set returned by
            ``collect_known_disease_canonical_motifs``.

    Returns:
        1 if ``canonical_motif`` is in the known set, else 0.
    """
    return 1 if canonical_motif in known_canonical_motifs else 0


def is_in_mendelian_gene(gene_id, gene_lookup):
    """Return 1 if a gene_id is in the Mendelian gene-disease table, else 0.

    Args:
        gene_id: The locus's gene id, or None when the locus has no gene
            annotation.
        gene_lookup: A dict keyed by gene id (the Mendelian gene-disease table).

    Returns:
        1 if ``gene_id`` is a non-None key of ``gene_lookup``, else 0.
    """
    return 1 if gene_id is not None and gene_id in gene_lookup else 0
