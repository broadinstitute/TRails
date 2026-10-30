"""Unit tests for allele_histograms.py."""

import unittest

from allele_histograms import (
    parse_genotype,
    accumulate_locus,
    convert_counts_to_histogram_string,
    convert_sample_ids_to_string,
    build_histograms_and_outliers,
)


class ParseGenotypeTests(unittest.TestCase):

    def test_diploid(self):
        self.assertEqual(parse_genotype("12,40"), (12, 40))

    def test_hemizygous(self):
        self.assertEqual(parse_genotype("21"), (21,))

    def test_homozygous(self):
        self.assertEqual(parse_genotype("2,2"), (2, 2))

    def test_whitespace_is_stripped(self):
        self.assertEqual(parse_genotype("  10,11 "), (10, 11))

    def test_no_call_values_return_none(self):
        for cell in ("", ".", "./.", "nan", "NaN", "NA", None):
            self.assertIsNone(parse_genotype(cell), msg=cell)

    def test_unparseable_returns_none(self):
        self.assertIsNone(parse_genotype("abc"))
        self.assertIsNone(parse_genotype("1,2,3"))
        self.assertIsNone(parse_genotype("1.5"))


class HistogramStringTests(unittest.TestCase):

    def test_ascending_order(self):
        self.assertEqual(
            convert_counts_to_histogram_string({11: 178, 9: 6, 10: 4086}),
            "9x:6,10x:4086,11x:178",
        )

    def test_empty(self):
        self.assertEqual(convert_counts_to_histogram_string({}), "")


class AccumulateLocusTests(unittest.TestCase):

    def test_diploid_10_11(self):
        """Blueprint worked example: diploid 10,11."""
        all_pair, short_pair, hemizygous_pair = accumulate_locus({"s1": "10,11"})
        self.assertEqual(convert_counts_to_histogram_string(all_pair[0]), "10x:1,11x:1")
        self.assertEqual(convert_counts_to_histogram_string(short_pair[0]), "10x:1")
        self.assertEqual(convert_counts_to_histogram_string(hemizygous_pair[0]), "")

    def test_hemizygous_3(self):
        """Blueprint worked example: hemizygous 3 -> all three histograms 3x:1."""
        all_pair, short_pair, hemizygous_pair = accumulate_locus({"s1": "3"})
        self.assertEqual(convert_counts_to_histogram_string(all_pair[0]), "3x:1")
        self.assertEqual(convert_counts_to_histogram_string(short_pair[0]), "3x:1")
        self.assertEqual(convert_counts_to_histogram_string(hemizygous_pair[0]), "3x:1")

    def test_homozygous_double_counts_in_all_and_short_single(self):
        """Reference docstring example: 2,2 -> All '2x:2', Short '2x:1'."""
        all_pair, short_pair, hemizygous_pair = accumulate_locus({"s1": "2,2"})
        self.assertEqual(convert_counts_to_histogram_string(all_pair[0]), "2x:2")
        self.assertEqual(convert_counts_to_histogram_string(short_pair[0]), "2x:1")
        self.assertEqual(convert_counts_to_histogram_string(hemizygous_pair[0]), "")

    def test_no_call_cells_skipped(self):
        all_pair, short_pair, hemizygous_pair = accumulate_locus(
            {"s1": "10,11", "s2": "", "s3": ".", "s4": "./.", "s5": None})
        self.assertEqual(convert_counts_to_histogram_string(all_pair[0]), "10x:1,11x:1")

    def test_short_allele_is_min_regardless_of_order(self):
        all_pair, short_pair, _ = accumulate_locus({"s1": "40,12"})
        self.assertEqual(convert_counts_to_histogram_string(all_pair[0]), "12x:1,40x:1")
        self.assertEqual(convert_counts_to_histogram_string(short_pair[0]), "12x:1")

    def test_sample_id_tracking_per_allele(self):
        all_pair, short_pair, _ = accumulate_locus({"s1": "10,11", "s2": "10,12"})
        self.assertEqual(all_pair[1][10], ["s1", "s2"])
        self.assertEqual(all_pair[1][11], ["s1"])
        self.assertEqual(all_pair[1][12], ["s2"])

    def test_homozygous_sample_id_recorded_once_but_counts_twice_toward_cap(self):
        """A homozygous sample routes its allele twice (counts twice toward the
        2*n cap) but appears only once in the sample-id list."""
        all_pair, _, _ = accumulate_locus({"s1": "5,5"}, n_outlier_sample_ids=10)
        self.assertEqual(all_pair[0][5], 2)
        self.assertEqual(all_pair[1][5], ["s1"])

    def test_occurrence_cap_blocks_new_sample_ids(self):
        """With n_outlier_sample_ids=1, the cap is 2 occurrences for AllAlleles.
        Three homozygous samples on the same allele each add 2 occurrences, so
        only the first sample id is recorded (cap reached after sample 1)."""
        genotypes = {"s1": "5,5", "s2": "5,5", "s3": "5,5"}
        all_pair, _, _ = accumulate_locus(genotypes, n_outlier_sample_ids=1)
        self.assertEqual(all_pair[0][5], 6)
        self.assertEqual(all_pair[1][5], ["s1"])


class ConvertSampleIdsToStringTests(unittest.TestCase):

    def test_descending_by_allele(self):
        result = convert_sample_ids_to_string({10: ["b"], 12: ["a"], 11: ["c"]})
        self.assertEqual(result, "12x:a,11x:c,10x:b")

    def test_ascending_sample_ids_within_allele(self):
        result = convert_sample_ids_to_string({12: ["c", "a", "b"]})
        self.assertEqual(result, "12x:a,12x:b,12x:c")

    def test_skip_common_allele_at_or_above_n(self):
        """An allele with >= n samples causes a break (largest-first), so it and
        all smaller alleles are skipped."""
        allele_to_sample_ids = {
            20: ["x", "y", "z"],          # 3 >= n=3 -> break immediately
            10: ["a"],
        }
        self.assertEqual(convert_sample_ids_to_string(allele_to_sample_ids, n_outlier_sample_ids=3), "")

    def test_smaller_common_allele_stops_iteration(self):
        """The largest allele qualifies; a smaller common allele then breaks."""
        allele_to_sample_ids = {
            20: ["a"],                    # emitted
            15: ["p", "q", "r"],          # 3 >= n=3 -> break before emitting
            10: ["z"],
        }
        self.assertEqual(
            convert_sample_ids_to_string(allele_to_sample_ids, n_outlier_sample_ids=3),
            "20x:a",
        )

    def test_stop_once_more_than_n_emitted(self):
        """After emitting an allele's samples, if total emitted > n, stop."""
        allele_to_sample_ids = {
            20: ["a", "b"],
            19: ["c", "d", "e"],          # after this, 5 emitted > n=4 -> stop
            18: ["f"],
        }
        self.assertEqual(
            convert_sample_ids_to_string(allele_to_sample_ids, n_outlier_sample_ids=4),
            "20x:a,20x:b,19x:c,19x:d,19x:e",
        )

    def test_empty(self):
        self.assertEqual(convert_sample_ids_to_string({}), "")


class BuildHistogramsAndOutliersTests(unittest.TestCase):

    def test_augments_each_locus_row(self):
        locus_rows = [
            {"trid": "X-263540-263579-TTTA", "motif": "TTTA",
             "genotypes": {"sample1_xx": "10,11"}},
            {"trid": "X-95117-95209-GTC", "motif": "GTC",
             "genotypes": {"sample1_xy": "3"}},
        ]
        result = build_histograms_and_outliers(locus_rows, ["sample1_xx", "sample1_xy"])

        first = result[0]
        self.assertEqual(first["AllAlleleHistogram"], "10x:1,11x:1")
        self.assertEqual(first["ShortAlleleHistogram"], "10x:1")
        self.assertEqual(first["HemizygousAlleleHistogram"], "")
        # Outliers descending by allele: 11 then 10.
        self.assertEqual(first["OutlierSampleIds_AllAlleles"], "11x:sample1_xx,10x:sample1_xx")
        self.assertEqual(first["OutlierSampleIds_ShortAlleles"], "10x:sample1_xx")
        self.assertEqual(first["OutlierSampleIds_HemizygousAlleles"], "")

        second = result[1]
        self.assertEqual(second["AllAlleleHistogram"], "3x:1")
        self.assertEqual(second["ShortAlleleHistogram"], "3x:1")
        self.assertEqual(second["HemizygousAlleleHistogram"], "3x:1")
        self.assertEqual(second["OutlierSampleIds_AllAlleles"], "3x:sample1_xy")
        self.assertEqual(second["OutlierSampleIds_HemizygousAlleles"], "3x:sample1_xy")

    def test_sample_absent_from_row_treated_as_no_call(self):
        locus_rows = [{"genotypes": {"s1": "10,11"}}]
        result = build_histograms_and_outliers(locus_rows, ["s1", "s2"])
        self.assertEqual(result[0]["AllAlleleHistogram"], "10x:1,11x:1")

    def test_returns_same_list_mutated_in_place(self):
        locus_rows = [{"genotypes": {"s1": "5"}}]
        result = build_histograms_and_outliers(locus_rows, ["s1"])
        self.assertIs(result, locus_rows)
        self.assertIn("AllAlleleHistogram", locus_rows[0])


if __name__ == "__main__":
    unittest.main()
