"""Unit tests for mendelian_qc.py."""

import unittest

import pandas

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
        self.assertFalse(check_violation((10, 20), (10, 11), (20, 21), 2))

    def test_autosomal_consistent_swapped_assignment(self):
        # child (20, 10): 20 from father, 10 from mother (assignment 2).
        self.assertFalse(check_violation((20, 10), (10, 11), (20, 21), 2))

    def test_autosomal_violation(self):
        # child (50, 60): neither allele can come from a parent.
        self.assertTrue(check_violation((50, 60), (10, 11), (20, 21), 2))

    def test_chrx_hemizygous_child_consistent(self):
        # hemizygous son: single X allele must come from mother.
        self.assertFalse(check_violation((15,), (15, 30), (40,), 2))

    def test_chrx_hemizygous_child_violation(self):
        self.assertTrue(check_violation((15,), (30, 40), (15,), 2))

    def test_chrx_hemizygous_father_daughter_consistent(self):
        # daughter (12, 20): 12 from father's single allele, 20 from mother.
        self.assertFalse(check_violation((12, 20), (20, 21), (12,), 2))

    def test_chrx_hemizygous_father_daughter_violation(self):
        self.assertTrue(check_violation((50, 60), (20, 21), (12,), 2))

    def test_hemizygous_mother_consistent(self):
        # child (5, 30): 5 from mother's single allele, 30 from father.
        self.assertFalse(check_violation((5, 30), (5,), (30, 31), 2))

    def test_threshold_strictness_diff_equals_threshold_is_violation(self):
        # child allele 12 vs only possible source 10: diff == 2 == threshold,
        # which is NOT a match (strict <), so it's a violation.
        self.assertTrue(check_violation((12,), (10,), (99,), 2))
        # diff == 1 < threshold: a match, no violation.
        self.assertFalse(check_violation((11,), (10,), (99,), 2))


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


if __name__ == "__main__":
    unittest.main()
