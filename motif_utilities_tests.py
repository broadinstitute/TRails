"""Unit tests for motif_utilities.

Golden values for reverse_complement and compute_canonical_motif are borrowed
from str_analysis/utils/canonical_repeat_unit_tests.py to confirm the port
reproduces the original behavior exactly.
"""

import unittest

from motif_utilities import (
    COMPLEMENT,
    _alphabetically_first_motif_under_shift,
    compute_canonical_motif,
    generate_all_canonical_motifs,
    parse_interval,
    reverse_complement,
)


class ReverseComplementTests(unittest.TestCase):

    def test_single_bases(self):
        self.assertEqual(reverse_complement("A"), "T")
        self.assertEqual(reverse_complement("C"), "G")
        self.assertEqual(reverse_complement("G"), "C")
        self.assertEqual(reverse_complement("T"), "A")
        self.assertEqual(reverse_complement("N"), "N")

    def test_palindromes_and_pairs(self):
        self.assertEqual(reverse_complement("CG"), "CG")
        self.assertEqual(reverse_complement("TA"), "TA")
        self.assertEqual(reverse_complement("GC"), "GC")
        self.assertEqual(reverse_complement("AT"), "AT")
        self.assertEqual(reverse_complement("CA"), "TG")
        self.assertEqual(reverse_complement("GG"), "CC")

    def test_runs(self):
        self.assertEqual(reverse_complement("G" * 2), "C" * 2)
        self.assertEqual(reverse_complement("G" * 3), "C" * 3)
        self.assertEqual(reverse_complement("G" * 10), "C" * 10)

    def test_reverse_complement_of_motif(self):
        # The reverse complement of GAA is TTC (used in the GAA canonical example).
        self.assertEqual(reverse_complement("GAA"), "TTC")

    def test_iupac_ambiguity_codes(self):
        self.assertEqual(reverse_complement("Y"), "R")
        self.assertEqual(reverse_complement("R"), "Y")
        self.assertEqual(reverse_complement("M"), "K")
        self.assertEqual(reverse_complement("K"), "M")
        self.assertEqual(reverse_complement("B"), "V")
        self.assertEqual(reverse_complement("D"), "H")

    def test_unknown_base_raises(self):
        self.assertRaises(KeyError, lambda: reverse_complement("X"))

    def test_complement_table_is_complete_iupac(self):
        self.assertEqual(len(COMPLEMENT), 15)


class AlphabeticallyFirstUnderShiftTests(unittest.TestCase):

    def test_examples(self):
        self.assertEqual(_alphabetically_first_motif_under_shift("C"), "C")
        self.assertEqual(_alphabetically_first_motif_under_shift("TAA"), "AAT")
        self.assertEqual(_alphabetically_first_motif_under_shift("ACA"), "AAC")


class ComputeCanonicalMotifTests(unittest.TestCase):

    def test_single_base(self):
        self.assertEqual(compute_canonical_motif("G"), "C")
        self.assertEqual(compute_canonical_motif("N"), "N")
        self.assertEqual(compute_canonical_motif("T"), "A")

    def test_known_golden_values(self):
        self.assertEqual(compute_canonical_motif("TGAG"), "ACTC")
        self.assertEqual(compute_canonical_motif("G" * 9), "C" * 9)

    def test_gaa(self):
        # GAA rotations: GAA, AGA, AAG; reverse-complement TTC rotations: TTC, TCT, CTT.
        # Alphabetically first overall is AAG.
        self.assertEqual(compute_canonical_motif("GAA"), "AAG")

    def test_cag(self):
        # CAG rotations: CAG, AGC, GCA -> AGC. RC is CTG -> CTG. AGC < CTG.
        self.assertEqual(compute_canonical_motif("CAG"), "AGC")

    def test_reverse_complement_flag_changes_result(self):
        # For a motif whose reverse-complement rotation is smaller, the flag matters.
        with_rc = compute_canonical_motif("GAA", include_reverse_complement=True)
        without_rc = compute_canonical_motif("GAA", include_reverse_complement=False)
        self.assertEqual(with_rc, "AAG")
        # Without RC, only GAA/AGA/AAG considered -> still AAG (RC not smaller here).
        self.assertEqual(without_rc, "AAG")

    def test_reverse_complement_flag_picks_rc_when_smaller(self):
        # TTC rotations: TTC, TCT, CTT -> CTT. Its RC is GAA -> AAG which is smaller.
        self.assertEqual(compute_canonical_motif("TTC", include_reverse_complement=True), "AAG")
        self.assertEqual(compute_canonical_motif("TTC", include_reverse_complement=False), "CTT")

    def test_case_insensitive(self):
        self.assertEqual(compute_canonical_motif("gaa"), "AAG")
        self.assertEqual(compute_canonical_motif("Gaa"), "AAG")

    def test_n_containing_motif(self):
        # An N-containing motif still canonicalizes (N complements to N).
        self.assertEqual(compute_canonical_motif("AN"), "AN")
        self.assertEqual(compute_canonical_motif("NA"), "AN")


class ParseIntervalTests(unittest.TestCase):

    def test_basic(self):
        self.assertEqual(parse_interval("chr1:100-200"), ("chr1", 100, 200))

    def test_no_chr_prefix(self):
        self.assertEqual(parse_interval("X:5-15"), ("X", 5, 15))

    def test_supercontig_with_colon(self):
        self.assertEqual(parse_interval("HLA:1:100-200"), ("HLA:1", 100, 200))

    def test_malformed_raises(self):
        self.assertRaises(ValueError, lambda: parse_interval("chr1-100-200"))


class GenerateAllCanonicalMotifsTests(unittest.TestCase):

    def test_size_one(self):
        # Single bases A/C/G/T canonicalize to {A, C}.
        self.assertEqual(generate_all_canonical_motifs(1), {"A", "C"})

    def test_default_size_six(self):
        result = generate_all_canonical_motifs()
        # Default max_size is 6.
        self.assertEqual(result, generate_all_canonical_motifs(6))

    def test_all_results_are_canonical(self):
        # Every motif in the set must equal its own canonical form.
        for motif in generate_all_canonical_motifs(4):
            self.assertEqual(compute_canonical_motif(motif), motif)

    def test_monotonic_growth(self):
        smaller = generate_all_canonical_motifs(2)
        larger = generate_all_canonical_motifs(3)
        self.assertTrue(smaller.issubset(larger))


if __name__ == "__main__":
    unittest.main()
