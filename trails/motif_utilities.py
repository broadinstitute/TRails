"""Tandem-repeat motif helpers for TRails.

  * ``reverse_complement`` over the full IUPAC ``COMPLEMENT`` table,
  * ``compute_canonical_motif``: the alphabetically first rotation of the motif or
    of its reverse complement,
  * ``compute_canonical_motifs_with_n_as_any_base``: the canonical motifs of every
    concrete motif an N-containing motif stands for, and
  * ``generate_all_canonical_motifs``: the set of all canonical motifs up to a length.

All functions are pure: they take their inputs as arguments and return values
without mutating shared state.
"""

import itertools


# Full IUPAC complement table, including the ambiguity codes. Used by
# reverse_complement, which is in turn used by compute_canonical_motif.
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
    are complemented too. A base not present in the table raises ``KeyError``.

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


def compute_canonical_motif(motif):
    """Return the canonical representation of a tandem-repeat motif.

    The canonical motif is the rearrangement of the motif's bases -- considering
    all cyclic rotations of the motif and of its reverse complement -- that is
    alphabetically first. For example, ``compute_canonical_motif("GAA")`` returns
    ``"AAG"`` (the first among ``GAA``, ``AGA``, ``AAG``, and the
    reverse-complement rotations ``TTC``, ``TCT``, ``CTT``).

    Args:
        motif: A repeat motif such as ``"CAG"``. Case-insensitive (upcased
            internally).

    Returns:
        The canonical (alphabetically first) motif string.
    """
    motif = motif.upper()
    return min(_alphabetically_first_motif_under_shift(motif),
               _alphabetically_first_motif_under_shift(reverse_complement(motif)))


def compute_canonical_motifs_with_n_as_any_base(motif):
    """Return the canonical motifs of every concrete motif that ``motif`` stands for.

    Each ``N`` stands for any of A, C, G and T, so the polyalanine catalog motif
    ``"GCN"`` stands for GCA, GCC, GCG and GCT. Two motifs match, N included, when
    their sets share a canonical motif. A motif without N gives just its own
    canonical motif. The number of concrete motifs is 4 to the power of the number
    of Ns, so callers expand only short motifs.

    Args:
        motif: A repeat motif such as ``"GCN"``. Case-insensitive.

    Returns:
        A set of canonical motif strings (``compute_canonical_motif``).
    """
    base_choices = [("A", "C", "G", "T") if base == "N" else (base,) for base in motif.upper()]
    return {compute_canonical_motif("".join(bases)) for bases in itertools.product(*base_choices)}


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
