"""Unit tests for mendelian_qc.py."""

import unittest

import pandas

import mendelian_qc
from mendelian_qc import (
    ALL_CANONICAL_MOTIFS,
    MOTIF_SIZE_CATEGORIES,
    allele_matches_any,
    check_violation,
    compute_mendelian_violations,
    find_trios,
    get_chrom_category,
    get_motif_size_category,
    parse_genotype,
    per_motif_columns,
    per_sample_columns,
)


class ParseGenotypeTests(unittest.TestCase):

    def test_diploid(self):
        self.assertEqual(parse_genotype("12,40"), (12, 40))

    def test_hemizygous(self):
        self.assertEqual(parse_genotype("21"), (21,))

    def test_no_call_variants(self):
        for cell in ("", ".", "./.", "  ", None):
            self.assertIsNone(parse_genotype(cell), cell)

    def test_malformed(self):
        self.assertIsNone(parse_genotype("1,2,3"))
        self.assertIsNone(parse_genotype("abc"))
        self.assertIsNone(parse_genotype("12,x"))

    def test_strips_whitespace(self):
        self.assertEqual(parse_genotype("  3 "), (3,))


class AlleleMatchesAnyTests(unittest.TestCase):

    def test_match_within_threshold(self):
        self.assertTrue(allele_matches_any(10, (11, 30), 2))

    def test_no_match(self):
        self.assertFalse(allele_matches_any(10, (30, 40), 2))

    def test_strict_threshold(self):
        # diff == threshold is NOT a match (strict <).
        self.assertFalse(allele_matches_any(10, (12,), 2))
        self.assertTrue(allele_matches_any(10, (11,), 2))


class CheckViolationTests(unittest.TestCase):

    def test_autosomal_consistent(self):
        # child (10, 20): 10 from mother (10,11), 20 from father (20,21).
        self.assertFalse(check_violation((10, 20), (10, 11), (20, 21), 2, "autosome"))

    def test_autosomal_consistent_swapped_assignment(self):
        # child (20, 10): 20 from father, 10 from mother (assignment 2).
        self.assertFalse(check_violation((20, 10), (10, 11), (20, 21), 2, "autosome"))

    def test_autosomal_violation(self):
        # child (50, 60): neither allele can come from a parent.
        self.assertTrue(check_violation((50, 60), (10, 11), (20, 21), 2, "autosome"))

    def test_autosomal_single_allele_child_can_come_from_either_parent(self):
        # A one-allele call on an autosome is not a hemizygous son: the allele may have come
        # from the father, so 30 matching the father (30, 31) is consistent, not a violation.
        self.assertFalse(check_violation((30,), (10, 11), (30, 31), 2, "autosome"))
        self.assertFalse(check_violation((10,), (10, 11), (30, 31), 2, "autosome"))
        # Matching neither parent is still a violation.
        self.assertTrue(check_violation((50,), (10, 11), (30, 31), 2, "autosome"))

    def test_chrx_hemizygous_child_consistent(self):
        # hemizygous son: single X allele must come from mother.
        self.assertFalse(check_violation((15,), (15, 30), (40,), 2, "chrX"))

    def test_chrx_hemizygous_child_violation(self):
        self.assertTrue(check_violation((15,), (30, 40), (15,), 2, "chrX"))

    def test_chrx_hemizygous_father_daughter_consistent(self):
        # daughter (12, 20): 12 from father's single allele, 20 from mother.
        self.assertFalse(check_violation((12, 20), (20, 21), (12,), 2, "chrX"))

    def test_chrx_hemizygous_father_daughter_violation(self):
        self.assertTrue(check_violation((50, 60), (20, 21), (12,), 2, "chrX"))

    def test_hemizygous_mother_consistent(self):
        # child (5, 30): 5 from mother's single allele, 30 from father.
        self.assertFalse(check_violation((5, 30), (5,), (30, 31), 2, "autosome"))

    def test_threshold_strictness_diff_equals_threshold_is_violation(self):
        # child allele 12 vs only possible source 10: diff == 2 == threshold,
        # which is NOT a match (strict <), so it's a violation.
        self.assertTrue(check_violation((12,), (10,), (99,), 2, "chrX"))
        # diff == 1 < threshold: a match, no violation.
        self.assertFalse(check_violation((11,), (10,), (99,), 2, "chrX"))


class ChromCategoryTests(unittest.TestCase):

    def test_autosome(self):
        self.assertEqual(get_chrom_category("chr1-100-110-AT"), "autosome")
        self.assertEqual(get_chrom_category("chr12_200_210_CAG"), "autosome")

    def test_chrx(self):
        self.assertEqual(get_chrom_category("chrX-1-2-A"), "chrX")
        self.assertEqual(get_chrom_category("X_1_2_A"), "chrX")

    def test_chry(self):
        self.assertEqual(get_chrom_category("chrY-1-2-A"), "chrY")

    def test_chrm(self):
        self.assertEqual(get_chrom_category("chrM-1-2-A"), "chrM")
        self.assertEqual(get_chrom_category("chrMT_1_2_A"), "chrM")


class MotifSizeCategoryTests(unittest.TestCase):

    def test_small_motifs(self):
        self.assertEqual(get_motif_size_category("A"), "1bp")
        self.assertEqual(get_motif_size_category("CAG"), "3bp")
        self.assertEqual(get_motif_size_category("ACGTAG"), "6bp")

    def test_medium_motif(self):
        self.assertEqual(get_motif_size_category("A" * 10), "7-24bp")

    def test_large_motif(self):
        self.assertEqual(get_motif_size_category("A" * 30), "25+bp")


class FindTriosTests(unittest.TestCase):

    def _df(self, rows):
        return pandas.DataFrame(rows)

    def test_complete_trio(self):
        df = self._df([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        self.assertEqual(find_trios(df), [("child", "mom", "dad")])

    def test_incomplete_missing_parent_in_set(self):
        df = self._df([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
        ])
        self.assertEqual(find_trios(df), [])

    def test_blank_parent_skipped(self):
        df = self._df([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": ""},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
        ])
        self.assertEqual(find_trios(df), [])

    def test_accepts_sample_lookup_dict(self):
        lookup = {
            "child": {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            "mom": {"sample_id": "mom"},
            "dad": {"sample_id": "dad"},
        }
        self.assertEqual(find_trios(lookup), [("child", "mom", "dad")])

    def test_duplicate_child_raises(self):
        # Two rows with the same child_id forming a trio -> duplicate.
        df = self._df([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        with self.assertRaises(ValueError):
            find_trios(df)


class ColumnLayoutTests(unittest.TestCase):

    def test_per_sample_columns_layout(self):
        columns = per_sample_columns()
        self.assertEqual(columns[0], "sample_id")
        self.assertIn("autosome_violations", columns)
        self.assertIn("autosome_total", columns)
        self.assertIn("chrM_total", columns)
        self.assertIn("motif_7_24bp_violations", columns)
        self.assertIn("motif_25plusbp_total", columns)
        self.assertEqual(columns[-2], "total_violations")
        self.assertEqual(columns[-1], "total_loci")
        # 1 + 4*2 + 8*2 + 2 = 27 columns.
        self.assertEqual(len(columns), 27)

    def test_per_motif_columns_layout(self):
        columns = per_motif_columns()
        self.assertEqual(columns[0], "sample_id")
        self.assertEqual(len(columns), 1 + 2 * len(ALL_CANONICAL_MOTIFS))
        self.assertIn("mv_AT", columns)
        self.assertIn("total_AT", columns)


class ComputeMendelianViolationsTests(unittest.TestCase):

    def test_no_trios_returns_empty(self):
        sample_df = pandas.DataFrame([{"sample_id": "s1"}, {"sample_id": "s2"}])
        per_sample, per_motif = compute_mendelian_violations(
            [], {}, sample_df, threshold=2)
        self.assertEqual(per_sample, [])
        self.assertEqual(per_motif, [])

    def test_end_to_end_three_loci_one_trio(self):
        # One complete trio, three loci:
        #   locus 1 (autosome, CAG): consistent.
        #   locus 2 (autosome, AT):  violation (child can't inherit either allele).
        #   locus 3 (autosome, single-distinct-allele): filtered out (<2 distinct).
        locus_rows = [
            {
                "trid": "chr1-100-110-CAG", "motif": "CAG",
                "genotypes": {
                    "child": "10,20", "mom": "10,11", "dad": "20,21", "sib": "5,5",
                },
            },
            {
                "trid": "chr2-200-210-AT", "motif": "AT",
                "genotypes": {
                    "child": "50,60", "mom": "10,11", "dad": "20,21", "sib": "",
                },
            },
            {
                "trid": "chr3-300-310-CAG", "motif": "CAG",
                "genotypes": {
                    # all members the same single allele -> <2 distinct -> skipped.
                    "child": "10,10", "mom": "10,10", "dad": "10,10", "sib": "10,10",
                },
            },
        ]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "sib", "maternal_id": "", "paternal_id": ""},
        ])
        sample_lookup = {row["sample_id"]: row for row in sample_df.to_dict(orient="records")}

        per_sample, per_motif = compute_mendelian_violations(
            locus_rows, sample_lookup, sample_df, threshold=2)

        self.assertEqual(len(per_sample), 1)
        self.assertEqual(len(per_motif), 1)
        row = per_sample[0]
        self.assertEqual(row["sample_id"], "child")
        # Two autosome loci counted (locus 3 filtered for <2 distinct alleles).
        self.assertEqual(row["autosome_total"], 2)
        self.assertEqual(row["autosome_violations"], 1)
        self.assertEqual(row["chrX_total"], 0)
        self.assertEqual(row["chrY_total"], 0)
        self.assertEqual(row["chrM_total"], 0)
        self.assertEqual(row["total_loci"], 2)
        self.assertEqual(row["total_violations"], 1)
        # by motif-size: CAG -> 3bp (locus 1 only, locus 3 filtered); AT -> 2bp (locus 2).
        self.assertEqual(row["motif_3bp_total"], 1)
        self.assertEqual(row["motif_3bp_violations"], 0)
        self.assertEqual(row["motif_2bp_total"], 1)
        self.assertEqual(row["motif_2bp_violations"], 1)

        # per-motif table: canonical(CAG) and canonical(AT).
        from motif_utilities import compute_canonical_motif
        cag = compute_canonical_motif("CAG")
        at = compute_canonical_motif("AT")
        motif_row = per_motif[0]
        self.assertEqual(motif_row["sample_id"], "child")
        self.assertEqual(motif_row[f"total_{cag}"], 1)
        self.assertEqual(motif_row[f"mv_{cag}"], 0)
        self.assertEqual(motif_row[f"total_{at}"], 1)
        self.assertEqual(motif_row[f"mv_{at}"], 1)

    def test_chry_uses_father_only(self):
        # chrY locus: child + father, mother irrelevant. Child allele not from
        # father -> violation.
        locus_rows = [{
            "trid": "chrY-1-10-A", "motif": "A",
            "genotypes": {"child": "30", "mom": "5,6", "dad": "10"},
        }]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["chrY_total"], 1)
        self.assertEqual(per_sample[0]["chrY_violations"], 1)
        self.assertEqual(per_sample[0]["autosome_total"], 0)

    def test_chry_concordant_locus_counts_toward_the_denominator(self):
        # A son whose chrY allele matches his father's is the normal case, and it has to be in
        # the denominator or chrY_violations / chrY_total is not a violation rate.
        locus_rows = [{
            "trid": "chrY-1-10-A", "motif": "A",
            "genotypes": {"child": "10", "mom": "5,6", "dad": "10"},
        }]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["chrY_total"], 1)
        self.assertEqual(per_sample[0]["chrY_violations"], 0)

    def test_chrm_uses_mother_only(self):
        # chrM locus: child + mother, father irrelevant. Child matches mother
        # -> no violation. Distinct allele provided by mother's second allele.
        locus_rows = [{
            "trid": "chrM-1-10-A", "motif": "A",
            "genotypes": {"child": "10", "mom": "10,40", "dad": "99"},
        }]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["chrM_total"], 1)
        self.assertEqual(per_sample[0]["chrM_violations"], 0)

    def test_chrx_hemizygous_child_evaluated_without_a_paternal_call(self):
        # A son's single X allele comes from the mother, so the comparison never needs the
        # father: the locus must be evaluated even though he has no call there.
        locus_rows = [
            {
                "trid": "chrX-1-10-A", "motif": "A",
                "genotypes": {"child": "10", "mom": "10,40", "dad": "."},
            },
            {
                "trid": "chrX-100-110-A", "motif": "A",
                "genotypes": {"child": "30", "mom": "10,11", "dad": "."},
            },
        ]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["chrX_total"], 2)
        self.assertEqual(per_sample[0]["chrX_violations"], 1)

    def test_chrx_diploid_child_still_requires_both_parents(self):
        # A daughter's two X alleles come one from each parent, so a missing paternal call
        # still means the locus cannot be evaluated.
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        without_father = [{
            "trid": "chrX-1-10-A", "motif": "A",
            "genotypes": {"child": "10,20", "mom": "10,11", "dad": "."},
        }]
        per_sample, _ = compute_mendelian_violations(without_father, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["chrX_total"], 0)

        with_father = [{
            "trid": "chrX-1-10-A", "motif": "A",
            "genotypes": {"child": "10,20", "mom": "10,11", "dad": "20,21"},
        }]
        per_sample, _ = compute_mendelian_violations(with_father, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["chrX_total"], 1)
        self.assertEqual(per_sample[0]["chrX_violations"], 0)

    def test_motif_with_n_skipped(self):
        locus_rows = [{
            "trid": "chr1-1-10-ANG", "motif": "ANG",
            "genotypes": {"child": "10,20", "mom": "5,6", "dad": "7,8"},
        }]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["total_loci"], 0)

    def test_motif_with_other_iupac_code_skipped(self):
        # compute_canonical_motif accepts the whole IUPAC alphabet, but the per-motif table is
        # keyed only by the ACGT canonical motifs, so an "R" motif used to KeyError here.
        locus_rows = [{
            "trid": "chr1-1-10-ARG", "motif": "ARG",
            "genotypes": {"child": "10,20", "mom": "5,6", "dad": "7,8"},
        }]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["total_loci"], 0)

    def test_motif_with_lowercase_bases_still_counted(self):
        # A lowercase motif is still an ACGT motif and must not be dropped by the
        # outside-ACGT skip.
        locus_rows = [{
            "trid": "chr1-1-10-cag", "motif": "cag",
            "genotypes": {"child": "10,20", "mom": "10,11", "dad": "20,21"},
        }]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["total_loci"], 1)

    def test_child_no_call_skipped(self):
        locus_rows = [{
            "trid": "chr1-1-10-A", "motif": "A",
            "genotypes": {"child": ".", "mom": "5,6", "dad": "7,8"},
        }]
        sample_df = pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])
        per_sample, _ = compute_mendelian_violations(locus_rows, {}, sample_df, threshold=2)
        self.assertEqual(per_sample[0]["total_loci"], 0)


class GenotypedTrioFilterTests(unittest.TestCase):
    """A trio described only in the metadata must not produce an all-zero row."""

    def _trio_metadata(self):
        return pandas.DataFrame([
            {"sample_id": "child", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])

    def test_trio_absent_from_the_matrix_is_dropped(self):
        locus_rows = [{
            "trid": "chr1-100-110-CAG", "motif": "CAG",
            "genotypes": {"other1": "10,20", "other2": "10,11"},
        }]
        per_sample, per_motif = compute_mendelian_violations(
            locus_rows, {}, self._trio_metadata(), threshold=2)
        self.assertEqual(per_sample, [])
        self.assertEqual(per_motif, [])

    def test_trio_with_one_ungenotyped_parent_is_dropped(self):
        locus_rows = [{
            "trid": "chr1-100-110-CAG", "motif": "CAG",
            "genotypes": {"child": "10,20", "mom": "10,11"},
        }]
        per_sample, per_motif = compute_mendelian_violations(
            locus_rows, {}, self._trio_metadata(), threshold=2)
        self.assertEqual(per_sample, [])
        self.assertEqual(per_motif, [])

    def test_genotyped_clean_trio_is_still_reported(self):
        # The clean trio must stay in the output and be distinguishable from a trio that was
        # never compared: zero violations over a non-zero denominator.
        locus_rows = [{
            "trid": "chr1-100-110-CAG", "motif": "CAG",
            "genotypes": {"child": "10,20", "mom": "10,11", "dad": "20,21"},
        }]
        per_sample, per_motif = compute_mendelian_violations(
            locus_rows, {}, self._trio_metadata(), threshold=2)
        self.assertEqual(len(per_sample), 1)
        self.assertEqual(len(per_motif), 1)
        self.assertEqual(per_sample[0]["sample_id"], "child")
        self.assertEqual(per_sample[0]["total_loci"], 1)
        self.assertEqual(per_sample[0]["total_violations"], 0)

    def test_genotyped_trio_with_only_no_calls_is_kept(self):
        # The members are genotype columns of the matrix, so the trio is real even though every
        # locus is skipped; the zero denominator is what marks it as uncompared.
        locus_rows = [{
            "trid": "chr1-100-110-CAG", "motif": "CAG",
            "genotypes": {"child": ".", "mom": ".", "dad": "."},
        }]
        per_sample, _ = compute_mendelian_violations(
            locus_rows, {}, self._trio_metadata(), threshold=2)
        self.assertEqual(len(per_sample), 1)
        self.assertEqual(per_sample[0]["total_loci"], 0)


class PerLocusDerivationTests(unittest.TestCase):
    """The per-locus values must be derived once per locus, not once per trio per locus."""

    def _two_trio_metadata(self):
        return pandas.DataFrame([
            {"sample_id": "child1", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "child2", "maternal_id": "mom", "paternal_id": "dad"},
            {"sample_id": "mom", "maternal_id": "", "paternal_id": ""},
            {"sample_id": "dad", "maternal_id": "", "paternal_id": ""},
        ])

    def _two_locus_rows(self):
        return [
            {"trid": "chr1-100-110-CAG", "motif": "CAG",
             "genotypes": {"child1": "10,20", "child2": "10,21", "mom": "10,11",
                           "dad": "20,21"}},
            {"trid": "chr2-100-110-AT", "motif": "AT",
             "genotypes": {"child1": "5,6", "child2": "5,7", "mom": "5,8", "dad": "6,7"}},
        ]

    def test_canonical_motif_is_computed_once_per_locus_for_two_trios(self):
        calls = []
        real_compute_canonical_motif = mendelian_qc.compute_canonical_motif

        def counting_compute_canonical_motif(motif):
            calls.append(motif)
            return real_compute_canonical_motif(motif)

        mendelian_qc.compute_canonical_motif = counting_compute_canonical_motif
        try:
            per_sample, _ = compute_mendelian_violations(
                self._two_locus_rows(), {}, self._two_trio_metadata(), threshold=2)
        finally:
            mendelian_qc.compute_canonical_motif = real_compute_canonical_motif

        self.assertEqual(len(per_sample), 2)
        self.assertEqual(sorted(calls), ["AT", "CAG"])

    def test_loci_with_an_unusable_motif_are_dropped_before_the_trio_loop(self):
        # A blank motif and one carrying a non-ACGT base are not scorable, and the filter now
        # lives in the shared per-locus pass, so they must not reach any trio's tally.
        locus_rows = self._two_locus_rows() + [
            {"trid": "chr3-1-10-", "motif": "",
             "genotypes": {"child1": "1,2", "child2": "1,3", "mom": "1,4", "dad": "2,3"}},
            {"trid": "chr3-20-30-CNG", "motif": "CNG",
             "genotypes": {"child1": "1,2", "child2": "1,3", "mom": "1,4", "dad": "2,3"}},
        ]
        per_sample, _ = compute_mendelian_violations(
            locus_rows, {}, self._two_trio_metadata(), threshold=2)
        self.assertEqual([row["total_loci"] for row in per_sample], [2, 2])

    def test_prepared_loci_carry_the_derived_values(self):
        prepared = mendelian_qc._prepare_loci(self._two_locus_rows())
        self.assertEqual([(chrom, size, canonical)
                          for _genotypes, chrom, size, canonical in prepared],
                         [("autosome", "3bp", "AGC"), ("autosome", "2bp", "AT")])


if __name__ == "__main__":
    unittest.main()
