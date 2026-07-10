"""Unit tests for analysis_columns.py — the TRails outlier/affected analysis."""

import unittest

import analysis_columns
from analysis_columns import (
    OUTLIER_TYPES,
    OUTPUT_COLUMNS,
    add_all_outlier_columns,
    add_gene_columns,
    collect_samples_single_pass,
    compute_number_of_affected_unsolved_families_above_unaffected,
    compute_number_of_affected_unsolved_samples_above_unaffected,
    get_population_p99_threshold,
    is_above_population_p99,
    is_above_unaffected,
    is_affected_unsolved,
    is_missing_outlier_value,
    is_unaffected_or_solved_status,
    parse_outlier_entries,
)


def make_outlier_string(pairs):
    """Build an "{allele}x:{sample_id}" comma string from (allele, sample_id) pairs."""
    return ",".join(f"{allele}x:{sample_id}" for allele, sample_id in pairs)


class ParseOutlierEntriesTests(unittest.TestCase):

    def test_descending_sort_and_x_strip(self):
        # Input is not pre-sorted; output must be descending by allele size.
        entries = parse_outlier_entries("10x:s1,40x:s2,25x:s3")
        self.assertEqual([allele for allele, _, _, _ in entries], [40, 25, 10])
        self.assertEqual([sid for _, sid, _, _ in entries], ["s2", "s3", "s1"])
        # The trailing 'x' was stripped and alleles are ints.
        self.assertTrue(all(isinstance(allele, int) for allele, _, _, _ in entries))

    def test_purity_and_methylation_parsing(self):
        entries = parse_outlier_entries("40x:s2:0.97:0.5,25x:s3:.:.,10x:s1")
        # Sorted descending: 40, 25, 10.
        self.assertEqual(entries[0], (40, "s2", "0.97", "0.5"))
        # A literal "." becomes None for both purity and methylation.
        self.assertEqual(entries[1], (25, "s3", None, None))
        # Missing fields become None.
        self.assertEqual(entries[2], (10, "s1", None, None))

    def test_missing_values(self):
        for missing in ["", ".", "NA", "N/A", "nan", "NaN", "None", "null", None, float("nan")]:
            self.assertTrue(is_missing_outlier_value(missing), missing)
            self.assertEqual(parse_outlier_entries(missing), [])

    def test_non_missing_values(self):
        self.assertFalse(is_missing_outlier_value("10x:s1"))
        self.assertFalse(is_missing_outlier_value("0x:s1"))

    def test_malformed_entry_raises(self):
        # Fewer than two colon-separated fields.
        with self.assertRaises(ValueError):
            parse_outlier_entries("10x")
        # Non-integer allele.
        with self.assertRaises(ValueError):
            parse_outlier_entries("abc:s1")

    def test_empty_segments_skipped(self):
        # Trailing comma / empty segment is dropped, not treated as malformed.
        entries = parse_outlier_entries("10x:s1,")
        self.assertEqual(len(entries), 1)


class StatusTruthTableTests(unittest.TestCase):

    def test_unaffected_or_solved_truth_table(self):
        # affected=="unaffected" -> True regardless of analysis.
        self.assertTrue(is_unaffected_or_solved_status("unaffected", "unsolved"))
        self.assertTrue(is_unaffected_or_solved_status("Unaffected", "anything"))
        # analysis in {solved, unaffected, probably solved} -> True.
        self.assertTrue(is_unaffected_or_solved_status("affected", "solved"))
        self.assertTrue(is_unaffected_or_solved_status("affected", "Probably Solved"))
        self.assertTrue(is_unaffected_or_solved_status("affected", "unaffected"))
        # Grouped with affecteds (False):
        self.assertFalse(is_unaffected_or_solved_status("affected", "unsolved"))
        self.assertFalse(is_unaffected_or_solved_status("unknown", "unknown"))
        self.assertFalse(is_unaffected_or_solved_status("affected", "partially solved"))
        self.assertFalse(is_unaffected_or_solved_status("unknown", "unsolved"))
        self.assertFalse(is_unaffected_or_solved_status(None, None))

    def test_is_affected_unsolved_is_exact_complement(self):
        affected_lookup = {
            "u": "unaffected", "solved": "affected", "ps": "affected",
            "unk": "unknown", "aff": "affected",
        }
        analysis_lookup = {
            "u": "unsolved", "solved": "solved", "ps": "partially solved",
            "unk": "unknown", "aff": "unsolved",
        }
        for sample_id in affected_lookup:
            self.assertEqual(
                is_affected_unsolved(sample_id, affected_lookup, analysis_lookup),
                not is_unaffected_or_solved_status(
                    affected_lookup[sample_id], analysis_lookup[sample_id]),
                sample_id,
            )
        # Spot check the partially-solved-affected case is treated as affected/unsolved.
        self.assertTrue(is_affected_unsolved("ps", affected_lookup, analysis_lookup))


class PopulationThresholdTests(unittest.TestCase):

    def test_get_population_p99_threshold_max_of_present(self):
        row = {"HPRC256_99thPercentile": 30, "AoU1027_99thPercentile": 45}
        self.assertEqual(get_population_p99_threshold(row), 45)

    def test_get_population_p99_threshold_none_when_absent(self):
        self.assertIsNone(get_population_p99_threshold({}))
        self.assertIsNone(get_population_p99_threshold(
            {"HPRC256_99thPercentile": None, "AoU1027_99thPercentile": float("nan")}))

    def test_is_above_population_p99(self):
        row = {"HPRC256_99thPercentile": 30, "AoU1027_99thPercentile": 45}
        self.assertFalse(is_above_population_p99(45, row))
        self.assertTrue(is_above_population_p99(46, row))
        # No cohort data -> include.
        self.assertTrue(is_above_population_p99(1, {}))

    def test_is_above_unaffected(self):
        self.assertTrue(is_above_unaffected(5, {}, "AllAlleles"))  # NULL first unaffected
        self.assertFalse(is_above_unaffected(0, {}, "AllAlleles"))
        row = {"FirstUnaffectedAlleleSize_AllAlleles": 20}
        self.assertTrue(is_above_unaffected(21, row, "AllAlleles"))
        self.assertFalse(is_above_unaffected(20, row, "AllAlleles"))


class CollectSamplesTests(unittest.TestCase):

    def setUp(self):
        # 4 affected samples and 2 unaffected, descending alleles already.
        self.sample_lookup = {
            "a1": {"family_id": "F1", "phenotype_description": "seizures"},
            "a2": {"family_id": "F1", "phenotype_description": "ataxia"},
            "a3": {"family_id": "F2", "phenotype_description": "tremor"},
            "a4": {"family_id": "F3", "phenotype_description": ""},
            "u1": {"family_id": "F4", "phenotype_description": ""},
            "u2": {"family_id": "F5", "phenotype_description": ""},
        }
        self.affected_lookup = {
            "a1": "affected", "a2": "affected", "a3": "affected", "a4": "affected",
            "u1": "unaffected", "u2": "unaffected",
        }
        self.analysis_lookup = {
            "a1": "unsolved", "a2": "unsolved", "a3": "unsolved", "a4": "unsolved",
            "u1": "unsolved", "u2": "unsolved",
        }

    def test_first_second_third_are_largest_affected(self):
        entries = parse_outlier_entries(make_outlier_string(
            [(50, "a1"), (40, "a2"), (30, "a3"), (20, "a4"), (15, "u1"), (10, "u2")]))
        collected = collect_samples_single_pass(
            entries, self.sample_lookup, self.affected_lookup, self.analysis_lookup)
        affected = collected["sample_id"]["affected"]
        self.assertEqual([entry["allele_size"] for entry in affected], [50, 40, 30])
        self.assertEqual([entry["sample_id"] for entry in affected], ["a1", "a2", "a3"])
        unaffected = collected["sample_id"]["unaffected"]
        self.assertEqual([entry["allele_size"] for entry in unaffected], [15, 10])

    def test_by_family_dedups_and_first_is_largest(self):
        # a1 and a2 share family F1 -> family-mode keeps only the first (largest, a1).
        entries = parse_outlier_entries(make_outlier_string(
            [(50, "a1"), (40, "a2"), (30, "a3"), (20, "a4")]))
        collected = collect_samples_single_pass(
            entries, self.sample_lookup, self.affected_lookup, self.analysis_lookup)
        by_family = collected["family_id"]["affected"]
        self.assertEqual([entry["family_id"] for entry in by_family], ["F1", "F2", "F3"])
        self.assertEqual([entry["allele_size"] for entry in by_family], [50, 30, 20])

    def test_by_family_skips_missing_family(self):
        sample_lookup = dict(self.sample_lookup)
        sample_lookup["a2"] = {"family_id": None, "phenotype_description": "x"}
        sample_lookup["a3"] = {"family_id": float("nan"), "phenotype_description": "x"}
        entries = parse_outlier_entries(make_outlier_string(
            [(50, "a1"), (40, "a2"), (30, "a3"), (20, "a4")]))
        collected = collect_samples_single_pass(
            entries, sample_lookup, self.affected_lookup, self.analysis_lookup)
        by_family = collected["family_id"]["affected"]
        # a2 (None) and a3 (NaN) skipped; only a1 (F1) and a4 (F3) remain.
        self.assertEqual([entry["sample_id"] for entry in by_family], ["a1", "a4"])

    def test_missing_metadata_row_kept_as_unknown_in_sample_ranking(self):
        # A matrix sample with no metadata row ("ghost") is kept (treated as
        # Unknown -> affected) in the sample-level ranking, but skipped in the
        # family-level ranking because it has no family_id. This matches the
        # affected/unsolved counters, which also count such samples.
        entries = parse_outlier_entries(make_outlier_string([(50, "ghost"), (40, "a1")]))
        collected = collect_samples_single_pass(
            entries, self.sample_lookup, self.affected_lookup, self.analysis_lookup)
        self.assertEqual([entry["sample_id"] for entry in collected["sample_id"]["affected"]],
                         ["ghost", "a1"])
        # ghost has no family_id -> only a1 (with a family) appears by family.
        self.assertEqual([entry["sample_id"] for entry in collected["family_id"]["affected"]],
                         ["a1"])


class NumAffectedUnsolvedSamplesTests(unittest.TestCase):

    def _parsed(self, row):
        return {ot: parse_outlier_entries(row.get(f"OutlierSampleIds_{ot}")) for ot in OUTLIER_TYPES}

    def test_stop_at_first_unaffected(self):
        # Descending: a1(50), a2(40), then an unaffected at 30 stops the walk.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(50, "a1"), (40, "a2"), (30, "u1"), (20, "a3")])}
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected", "u1": "unaffected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        counts = compute_number_of_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, analysis_lookup, self._parsed(row))
        # a1 and a2 counted; the break at u1 means a3 (below) is never reached.
        self.assertEqual(counts[0], 2)

    def test_stop_at_or_below_p99(self):
        # p99 = 40, so only the allele strictly above 40 (a1 at 50) is counted.
        row = {
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (40, "a2"), (30, "a3")]),
            "HPRC256_99thPercentile": 40,
        }
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        counts = compute_number_of_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, analysis_lookup, self._parsed(row))
        self.assertEqual(counts[0], 1)

    def test_missing_outlier_column_is_none(self):
        row = {"OutlierSampleIds_AllAlleles": ""}
        counts = compute_number_of_affected_unsolved_samples_above_unaffected(
            row, {}, {}, self._parsed(row))
        self.assertIsNone(counts[0])

    def test_distinct_sample_dedup_across_alleles(self):
        # The same sample at two allele sizes counts once.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (45, "a1")])}
        affected_lookup = {"a1": "affected"}
        analysis_lookup = {"a1": "unsolved"}
        counts = compute_number_of_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, analysis_lookup, self._parsed(row))
        self.assertEqual(counts[0], 1)

    def test_families_variant_dedups_by_family(self):
        # a1 and a2 share F1 -> one distinct family above unaffected.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (45, "a2")])}
        affected_lookup = {"a1": "affected", "a2": "affected"}
        analysis_lookup = {"a1": "unsolved", "a2": "unsolved"}
        sample_lookup = {"a1": {"family_id": "F1"}, "a2": {"family_id": "F1"}}
        counts = compute_number_of_affected_unsolved_families_above_unaffected(
            row, affected_lookup, analysis_lookup, sample_lookup, self._parsed(row))
        self.assertEqual(counts[0], 1)

    def test_families_variant_stop_at_unaffected(self):
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(50, "a1"), (40, "u1"), (30, "a2")])}
        affected_lookup = {"a1": "affected", "a2": "affected", "u1": "unaffected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        sample_lookup = {"a1": {"family_id": "F1"}, "a2": {"family_id": "F2"}, "u1": {"family_id": "F3"}}
        counts = compute_number_of_affected_unsolved_families_above_unaffected(
            row, affected_lookup, analysis_lookup, sample_lookup, self._parsed(row))
        self.assertEqual(counts[0], 1)


class AddAllOutlierColumnsTests(unittest.TestCase):

    def test_columns_populated_and_unaffected_null_preserved(self):
        records = [{
            "LocusId": "chr1-1-10-AAG",
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (40, "a2")]),
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": "",
        }]
        sample_lookup = {
            "a1": {"family_id": "F1", "phenotype_description": "seizures"},
            "a2": {"family_id": "F2", "phenotype_description": "ataxia"},
        }
        affected_lookup = {"a1": "affected", "a2": "affected"}
        analysis_lookup = {"a1": "unsolved", "a2": "unsolved"}

        add_all_outlier_columns(records, sample_lookup, affected_lookup, analysis_lookup)
        record = records[0]

        # First affected = largest (a1 at 50).
        self.assertEqual(record["FirstAffectedAlleleSize_AllAlleles"], 50)
        self.assertEqual(record["SecondAffectedAlleleSize_AllAlleles"], 40)
        self.assertIsNone(record["ThirdAffectedAlleleSize_AllAlleles"])
        self.assertEqual(record["FirstAffectedSampleId_AllAlleles"], "a1")
        self.assertEqual(record["FirstAffectedPhenotype_AllAlleles"], "seizures")
        # No unaffected -> NULL preserved (not 0).
        self.assertIsNone(record["FirstUnaffectedAlleleSize_AllAlleles"])
        # Counts present.
        self.assertEqual(record["NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles"], 2)
        # Missing outlier column -> count is None, allele sizes None.
        self.assertIsNone(record["NumAffectedUnsolvedSamplesAboveUnaffected_ShortAlleles"])
        self.assertIsNone(record["FirstAffectedAlleleSize_ShortAlleles"])

    def test_all_target_columns_present_on_record(self):
        records = [{
            "LocusId": "L",
            "OutlierSampleIds_AllAlleles": "",
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": "",
        }]
        add_all_outlier_columns(records, {}, {}, {})
        # Every analysis output column the loci table expects must be present.
        for outlier_type in OUTLIER_TYPES:
            for prefix in ["First", "Second", "Third"]:
                self.assertIn(f"{prefix}AffectedAlleleSize_{outlier_type}", records[0])
                self.assertIn(f"{prefix}AffectedAlleleSize_{outlier_type}_ByFamily", records[0])
                self.assertIn(f"{prefix}AffectedPhenotype_{outlier_type}", records[0])
                self.assertIn(f"{prefix}AffectedSampleId_{outlier_type}", records[0])
            for prefix in ["First", "Second"]:
                self.assertIn(f"{prefix}UnaffectedAlleleSize_{outlier_type}", records[0])
                self.assertIn(f"{prefix}UnaffectedAlleleSize_{outlier_type}_ByFamily", records[0])

    def test_by_family_dedup_in_full_pipeline(self):
        records = [{
            "LocusId": "L",
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (40, "a2"), (30, "a3")]),
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": "",
        }]
        sample_lookup = {
            "a1": {"family_id": "F1", "phenotype_description": ""},
            "a2": {"family_id": "F1", "phenotype_description": ""},
            "a3": {"family_id": "F2", "phenotype_description": ""},
        }
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        add_all_outlier_columns(records, sample_lookup, affected_lookup, analysis_lookup)
        # By sample: a1, a2, a3. By family: a1 (F1), a3 (F2).
        self.assertEqual(records[0]["FirstAffectedAlleleSize_AllAlleles_ByFamily"], 50)
        self.assertEqual(records[0]["SecondAffectedAlleleSize_AllAlleles_ByFamily"], 30)
        self.assertIsNone(records[0]["ThirdAffectedAlleleSize_AllAlleles_ByFamily"])


class AddGeneColumnsTests(unittest.TestCase):

    def test_pli_is_max_of_v2_v4(self):
        records = [{"gene_id": "G1"}, {"gene_id": "G2"}, {"gene_id": "G3"}, {"gene_id": "missing"}]
        gene_lookup = {
            "G1": {"pLI_v2": 0.3, "pLI_v4": 0.9, "inheritance": "AD"},
            "G2": {"pLI_v2": 0.5, "pLI_v4": float("nan"), "inheritance": "AR"},
            "G3": {"pLI_v2": float("nan"), "pLI_v4": float("nan"), "inheritance": "XR"},
        }
        add_gene_columns(records, gene_lookup)
        self.assertEqual(records[0]["pLI"], 0.9)
        self.assertEqual(records[0]["inheritance"], "AD")
        self.assertEqual(records[1]["pLI"], 0.5)
        self.assertIsNone(records[2]["pLI"])
        self.assertEqual(records[2]["inheritance"], "XR")
        # gene_id not in lookup -> both None.
        self.assertIsNone(records[3]["pLI"])
        self.assertIsNone(records[3]["inheritance"])

    def test_no_gene_lookup_yields_null(self):
        records = [{"gene_id": "G1"}, {"gene_id": None}]
        add_gene_columns(records, None)
        for record in records:
            self.assertIsNone(record["pLI"])
            self.assertIsNone(record["inheritance"])


class OutputColumnsTests(unittest.TestCase):

    def test_exact_count_and_uniqueness(self):
        self.assertEqual(len(OUTPUT_COLUMNS), 129)
        self.assertEqual(len(set(OUTPUT_COLUMNS)), 129)

    def test_anchor_positions(self):
        # Spot-check load-bearing positions from the blueprint's numbered schema.
        # Blueprint numbers are 1-based; subtract 1 for the list index.
        self.assertEqual(OUTPUT_COLUMNS[0], "LocusId")
        self.assertEqual(OUTPUT_COLUMNS[32], "pLI")  # col 33
        self.assertEqual(OUTPUT_COLUMNS[33], "inheritance")  # col 34
        self.assertEqual(OUTPUT_COLUMNS[34], "FirstAffectedAlleleSize_AllAlleles")  # col 35
        self.assertEqual(OUTPUT_COLUMNS[128], "VariationClusterSizeDiff")  # col 129


if __name__ == "__main__":
    unittest.main()
