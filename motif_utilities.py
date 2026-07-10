"""Standalone tandem-repeat motif primitives for TRails.

This module reimplements the small set of motif-handling helpers that the
original pipeline borrowed from the ``str_analysis`` package, so that TRails has
no external dependency on ``str_analysis`` and inserts nothing onto ``sys.path``.

The behavior here is a faithful port of:
  * ``str_analysis/utils/canonical_repeat_unit.py`` (the doubling-rotation +
    reverse-complement-minimum canonicalization algorithm), and
  * ``str_analysis/utils/misc_utils.py`` (``reverse_complement``, the full IUPAC
    ``COMPLEMENT`` table, and ``parse_interval``).

Plus a small TRails-specific helper (``generate_all_canonical_motifs``) used
elsewhere in the build pipeline.

All functions are pure: they take their inputs as arguments and return values
without mutating shared state.
"""

import itertools


# Full IUPAC complement table (matches str_analysis.utils.misc_utils.COMPLEMENT
# exactly, including the ambiguity codes). Used by reverse_complement, which is
# in turn used by compute_canonical_motif.
COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "Y": "R",
    "R": "Y",
    "S": "S",
    "W": "W",
    "M": "K",
    "K": "M",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
}


def reverse_complement(dna):
    """Return the reverse complement of a DNA string.

    Uses the full IUPAC complement table, so ambiguity codes (e.g. ``Y``, ``R``)
    are complemented just as in the original ``str_analysis`` implementation. A
    base not present in the table raises ``KeyError`` (matching the original).

    Args:
        dna: A string of DNA bases such as ``"GAA"``.

    Returns:
        The reverse-complement string, e.g. ``reverse_complement("GAA")`` is
        ``"TTC"``.
    """
    return "".join([COMPLEMENT[base] for base in dna[::-1]])


def _alphabetically_first_motif_under_shift(motif):
    """Return the rotation of ``motif`` that sorts alphabetically first.

    Considers every cyclic rotation of the motif (implemented by scanning a
    doubled copy of the motif) and returns the smallest one.

    Args:
        motif: A repeat motif such as ``"CAG"``.

    Returns:
        The alphabetically first cyclic rotation of the motif.
    """
    minimal_motif = motif
    double_motif = motif + motif
    for offset in range(len(motif)):
        if double_motif[offset:offset + len(motif)] < minimal_motif:
            minimal_motif = double_motif[offset:offset + len(motif)]
    return minimal_motif


def compute_canonical_motif(motif, include_reverse_complement=True):
    """Return the canonical representation of a tandem-repeat motif.

    The canonical motif is the rearrangement of the motif's bases -- considering
    all cyclic rotations, and optionally the rotations of its reverse complement
    -- that is alphabetically first. For example, ``compute_canonical_motif("GAA")``
    returns ``"AAG"`` (the first among ``GAA``, ``AGA``, ``AAG``, and the
    reverse-complement rotations ``TTC``, ``TCT``, ``CTT``).

    Args:
        motif: A repeat motif such as ``"CAG"``. Case-insensitive (upcased
            internally).
        include_reverse_complement: When True (default), the reverse complement
            of the motif and its rotations are also considered.

    Returns:
        The canonical (alphabetically first) motif string.
    """
    motif = motif.upper()
    minimal_motif = _alphabetically_first_motif_under_shift(motif)

    if include_reverse_complement:
        minimal_motif_reverse_complement = _alphabetically_first_motif_under_shift(
            reverse_complement(motif))
        if minimal_motif_reverse_complement < minimal_motif:
            return minimal_motif_reverse_complement

    return minimal_motif


def parse_interval(interval_string):
    """Parse a ``"chrom:start-end"`` interval string into a 3-tuple.

    Supports super-contig names that themselves contain ``:`` by treating only
    the final ``:``-delimited field as the coordinate range (matching the
    original ``str_analysis`` implementation).

    Args:
        interval_string: A string like ``"chr1:100-200"``.

    Returns:
        A 3-tuple ``(chrom, start, end)`` where ``start`` and ``end`` are ints.
    """
    try:
        tokens = interval_string.split(":")
        chrom = ":".join(tokens[:-1])
        start, end = map(int, tokens[-1].split("-"))
    except Exception as exception:
        raise ValueError(f"Unable to parse interval: '{interval_string}': {exception}")

    return chrom, start, end


def generate_all_canonical_motifs(max_size=6):
    """Return the set of all canonical motifs of size 1 through ``max_size``.

    Enumerates every possible motif over the four canonical bases (A, C, G, T)
    of each length from 1 up to ``max_size``, canonicalizes each one (including
    reverse complement), and returns the deduplicated set. Used to build the
    per-motif Mendelian-violation table.

    Args:
        max_size: The largest motif length to enumerate. Defaults to 6.

    Returns:
        A set of canonical motif strings.
    """
    canonical_motifs = set()
    for motif_length in range(1, max_size + 1):
        for bases in itertools.product("ACGT", repeat=motif_length):
            canonical_motifs.add(compute_canonical_motif("".join(bases)))
    return canonical_motifs
