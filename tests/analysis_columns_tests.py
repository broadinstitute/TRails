"""Unit tests for analysis_columns.py — the TRails outlier/affected analysis."""

import unittest

import analysis_columns
from analysis_columns import (
    OUTLIER_TYPES,
    OUTPUT_COLUMNS,
    add_all_outlier_columns,
    add_gene_columns,
    blank_to_none,
    collect_samples_single_pass,
    count_affected_unsolved_samples_and_families,
    find_affected_unsolved_samples_above_unaffected,
    get_family_id,
    get_population_p99_threshold,
    is_above_unaffected,
    is_missing_outlier_value,
    is_unaffected_or_solved_status,
    list_affected_unsolved_outliers,
    normalize_affected_status_for_logic,
    parse_outlier_entries,
    parse_population_gate_cohorts,
    title_case_status_for_display,
)


def make_outlier_string(pairs):
    """Build an "{allele}x:{sample_id}" comma string from (allele, sample_id) pairs."""
    return ",".join(f"{allele}x:{sample_id}" for allele, sample_id in pairs)


class ParseOutlierEntriesTests(unittest.TestCase):

    def test_descending_sort_and_x_strip(self):
        # Input is not pre-sorted; output must be descending by allele size.
        entries = parse_outlier_entries("10x:s1,40x:s2,25x:s3")
        # The trailing 'x' was stripped and alleles are ints.
        self.assertEqual(entries, [(40, "s2"), (25, "s3"), (10, "s1")])
        self.assertTrue(all(isinstance(allele, int) for allele, _ in entries))

    def test_missing_values(self):
        for missing in ["", ".", "NA", "N/A", "nan", "NaN", "None", "null", "NULL", None, float("nan")]:
            self.assertTrue(is_missing_outlier_value(missing), missing)
            self.assertEqual(parse_outlier_entries(missing), [])

    def test_non_missing_values(self):
        self.assertFalse(is_missing_outlier_value("10x:s1"))
        self.assertFalse(is_missing_outlier_value("0x:s1"))

    def test_malformed_entry_raises(self):
        # Fewer or more than two colon-separated fields (a sample_id cannot contain ':').
        for malformed in ("10x", "10x:s1:0.9"):
            with self.assertRaises(ValueError, msg=malformed):
                parse_outlier_entries(malformed)
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


class StatusNormalizationTests(unittest.TestCase):
    def test_affected_title_casing(self):
        self.assertEqual(title_case_status_for_display("affected"), "Affected")
        self.assertEqual(title_case_status_for_display("UNAFFECTED"), "Unaffected")
        self.assertEqual(title_case_status_for_display(" possibly affected "), "Possibly Affected")

    def test_analysis_title_casing(self):
        self.assertEqual(title_case_status_for_display("solved"), "Solved")
        self.assertEqual(title_case_status_for_display("UNSOLVED"), "Unsolved")
        self.assertEqual(title_case_status_for_display("partially solved"), "Partially Solved")

    def test_missing_and_blank_statuses_are_unknown(self):
        for blank in (None, float("nan"), "", "  ", "NA", "n/a", "null"):
            self.assertEqual(title_case_status_for_display(blank), "Unknown", repr(blank))


class BlankMetadataValueTests(unittest.TestCase):

    def test_blank_to_none(self):
        for blank in (None, float("nan"), "", "   ", ".", " . ", "NA", " na ", "N/A", "nan", "None", "NULL"):
            self.assertIsNone(blank_to_none(blank), repr(blank))
        self.assertEqual(blank_to_none("F1"), "F1")
        self.assertEqual(blank_to_none(" F1 "), " F1 ")
        self.assertEqual(blank_to_none("NAFAM"), "NAFAM")
        self.assertEqual(blank_to_none(".."), "..")

    def test_normalize_affected_status_for_logic(self):
        for blank in (None, float("nan"), "", ".", "NA", "null"):
            self.assertIsNone(normalize_affected_status_for_logic(blank), repr(blank))
        self.assertEqual(normalize_affected_status_for_logic(" Possibly Affected "), "affected")
        self.assertEqual(normalize_affected_status_for_logic("Affected"), "affected")
        self.assertEqual(normalize_affected_status_for_logic("UNAFFECTED"), "unaffected")
        self.assertEqual(normalize_affected_status_for_logic("Unknown"), "unknown")

    def test_one_list_of_missing_spellings(self):
        # Sample metadata, OutlierSampleIds_* fields and matrix annotation cells share one list.
        for spelling in analysis_columns.MISSING_VALUE_SPELLINGS:
            for value in (spelling, spelling.upper(), f" {spelling} "):
                self.assertIsNone(blank_to_none(value), repr(value))
                self.assertTrue(is_missing_outlier_value(value), repr(value))

    def test_get_family_id(self):
        sample_lookup = {"a": {"family_id": "F1"}, "b": {"family_id": ""},
                         "c": {"family_id": " "}, "d": {}, "e": {"family_id": "NA"}}
        self.assertEqual(get_family_id("a", sample_lookup), "F1")
        for sample_id in ("b", "c", "d", "e", "no_metadata_row"):
            self.assertIsNone(get_family_id(sample_id, sample_lookup), sample_id)

    def test_na_family_ids_do_not_merge_unrelated_samples(self):
        # Two unrelated samples whose family_id cells both say "NA" (read_sample_metadata
        # keeps that string) have no family: they are not one family in the counts, in the
        # _ByFamily ranking, or in the affected_unsolved_outliers_above_unaffected rows.
        records = [{
            "LocusId": "L",
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "S1"), (40, "S2"), (30, "S3")]),
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": "",
        }]
        sample_lookup = {"S1": {"family_id": "NA"}, "S2": {"family_id": "na"},
                         "S3": {"family_id": " "}}
        affected_lookup = {sample_id: "affected" for sample_id in sample_lookup}
        rows = add_all_outlier_columns(records, sample_lookup, affected_lookup, {})
        record = records[0]
        self.assertEqual(record["NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles"], 3)
        self.assertEqual(record["NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles"], 0)
        self.assertIsNone(record["FirstAffectedAlleleSize_AllAlleles_ByFamily"])
        self.assertEqual([row["family_id"] for row in rows], [None, None, None])

    def test_blank_phenotypes_are_stored_as_none(self):
        # The loci table's AffectedPhenotype columns follow the same blank rule as the swim plot's
        # phenotype_description; a sample without a metadata row has no phenotype either.
        records = [{
            "LocusId": "L",
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "S1"), (40, "S2"), (30, "S3")]),
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": "",
        }]
        sample_lookup = {"S1": {"phenotype_description": "NA"}, "S2": {"phenotype_description": ""}}
        add_all_outlier_columns(records, sample_lookup, {}, {})
        for prefix in ("First", "Second", "Third"):
            self.assertIsNone(records[0][f"{prefix}AffectedPhenotype_AllAlleles"], prefix)


class PopulationThresholdTests(unittest.TestCase):

    def test_get_population_p99_threshold_max_of_present(self):
        row = {"HPRC256_99thPercentile": 30, "AoU1027_99thPercentile": 45}
        self.assertEqual(get_population_p99_threshold(row), 45)

    def test_get_population_p99_threshold_none_when_absent(self):
        self.assertIsNone(get_population_p99_threshold({}))
        self.assertIsNone(get_population_p99_threshold(
            {"HPRC256_99thPercentile": None, "AoU1027_99thPercentile": float("nan")}))

    def test_population_gate_cohorts_choose_the_columns_read(self):
        row = {"HPRC256_99thPercentile": 30, "AoU1027_99thPercentile": 45,
               "TenK10K_99thPercentile": 60}
        # The default gate is the long-read cohorts only, so TenK10K is ignored.
        self.assertEqual(get_population_p99_threshold(row), 45)
        self.assertEqual(get_population_p99_threshold(row, ("HPRC256",)), 30)
        self.assertEqual(get_population_p99_threshold(row, ("HPRC256", "AoU1027", "TenK10K")), 60)

    def test_parse_population_gate_cohorts(self):
        self.assertEqual(parse_population_gate_cohorts("HPRC256,AoU1027"), ("HPRC256", "AoU1027"))
        self.assertEqual(parse_population_gate_cohorts("TenK10K, HPRC256"), ("TenK10K", "HPRC256"))
        with self.assertRaisesRegex(ValueError, r"value\(s\): HPRC\. Choices"):
            parse_population_gate_cohorts("HPRC,AoU1027")
        with self.assertRaises(ValueError):
            parse_population_gate_cohorts("")

    def test_is_above_unaffected(self):
        # With no unaffected sample (NULL first unaffected), every allele is above, even 0.
        self.assertTrue(is_above_unaffected(5, {}, "AllAlleles"))
        self.assertTrue(is_above_unaffected(0, {}, "AllAlleles"))
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

    def _counts(self, row, affected_lookup, analysis_lookup, sample_lookup=None,
                population_gate_cohorts=analysis_columns.POPULATION_P99_GATE_COHORTS):
        """Returns (sample counts, family counts) the way add_all_outlier_columns computes them."""
        samples = find_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, analysis_lookup, self._parsed(row), population_gate_cohorts)
        return count_affected_unsolved_samples_and_families(samples, sample_lookup or {})

    def test_stop_at_first_unaffected(self):
        # Descending: a1(50), a2(40), then an unaffected at 30 stops the walk.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(50, "a1"), (40, "a2"), (30, "u1"), (20, "a3")])}
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected", "u1": "unaffected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        sample_counts, _ = self._counts(row, affected_lookup, analysis_lookup)
        # a1 and a2 counted; the break at u1 means a3 (below) is never reached.
        self.assertEqual(sample_counts[0], 2)

    def test_tie_with_unaffected_is_not_above(self):
        # a2 shares allele 40 with the unaffected u1, so only a1 is above every unaffected.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(50, "a1"), (40, "a2"), (40, "u1")])}
        affected_lookup = {"a1": "affected", "a2": "affected", "u1": "unaffected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        sample_counts, _ = self._counts(row, affected_lookup, analysis_lookup)
        self.assertEqual(sample_counts[0], 1)

    def test_sample_also_at_the_tie_allele_keeps_its_larger_allele(self):
        # a1's second allele ties with u1, but its larger allele is above every unaffected.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(50, "a1"), (40, "a2"), (40, "a1"), (40, "u1")])}
        affected_lookup = {"a1": "affected", "a2": "affected", "u1": "unaffected"}
        samples = find_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, {}, self._parsed(row))
        self.assertEqual(samples[0], {"a1": 50})

    def test_stop_at_or_below_p99(self):
        # p99 = 40, so only the allele strictly above 40 (a1 at 50) is counted.
        row = {
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (40, "a2"), (30, "a3")]),
            "HPRC256_99thPercentile": 40,
        }
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        sample_counts, _ = self._counts(row, affected_lookup, analysis_lookup)
        self.assertEqual(sample_counts[0], 1)

    def test_population_gate_cohorts_choose_the_threshold(self):
        # TenK10K's 99th percentile of 45 only gates when it is one of the gate cohorts.
        row = {
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (40, "a2"), (30, "a3")]),
            "HPRC256_99thPercentile": 35,
            "TenK10K_99thPercentile": 45,
        }
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        self.assertEqual(self._counts(row, affected_lookup, analysis_lookup)[0][0], 2)
        self.assertEqual(self._counts(row, affected_lookup, analysis_lookup,
                                      population_gate_cohorts=("HPRC256", "TenK10K"))[0][0], 1)
        self.assertEqual(self._counts(row, affected_lookup, analysis_lookup,
                                      population_gate_cohorts=("AoU1027",))[0][0], 3)

    def test_missing_outlier_column_is_none(self):
        row = {"OutlierSampleIds_AllAlleles": ""}
        sample_counts, family_counts = self._counts(row, {}, {})
        self.assertIsNone(sample_counts[0])
        self.assertIsNone(family_counts[0])

    def test_distinct_sample_dedup_across_alleles(self):
        # The same sample at two allele sizes counts once, at its larger allele.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (45, "a1")])}
        affected_lookup = {"a1": "affected"}
        analysis_lookup = {"a1": "unsolved"}
        samples = find_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, analysis_lookup, self._parsed(row))
        self.assertEqual(samples[0], {"a1": 50})
        self.assertEqual(self._counts(row, affected_lookup, analysis_lookup)[0][0], 1)

    def test_samples_are_ordered_by_allele_then_entry_order(self):
        # phenotype_scoring compares consecutive samples in this order, so samples that share
        # an allele size must keep the order the outlier string lists them in.
        sample_ids = [f"s{index}" for index in range(8)]
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(40, "z")] + [(30, sample_id) for sample_id in sample_ids])}
        affected_lookup = {sample_id: "affected" for sample_id in sample_ids + ["z"]}
        samples = find_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, {}, self._parsed(row))
        self.assertEqual(list(samples[0]), ["z"] + sample_ids)

    def test_families_dedup_by_family(self):
        # a1 and a2 share F1 -> two samples but one distinct family above unaffected.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (45, "a2")])}
        affected_lookup = {"a1": "affected", "a2": "affected"}
        analysis_lookup = {"a1": "unsolved", "a2": "unsolved"}
        sample_lookup = {"a1": {"family_id": "F1"}, "a2": {"family_id": "F1"}}
        sample_counts, family_counts = self._counts(row, affected_lookup, analysis_lookup, sample_lookup)
        self.assertEqual(sample_counts[0], 2)
        self.assertEqual(family_counts[0], 1)

    def test_families_stop_at_unaffected(self):
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(50, "a1"), (40, "u1"), (30, "a2")])}
        affected_lookup = {"a1": "affected", "a2": "affected", "u1": "unaffected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        sample_lookup = {"a1": {"family_id": "F1"}, "a2": {"family_id": "F2"}, "u1": {"family_id": "F3"}}
        _, family_counts = self._counts(row, affected_lookup, analysis_lookup, sample_lookup)
        self.assertEqual(family_counts[0], 1)

    def test_sample_without_family_counts_toward_samples_only(self):
        # a2 has a blank family_id and "ghost" has no metadata row at all.
        row = {"OutlierSampleIds_AllAlleles": make_outlier_string(
            [(50, "a1"), (45, "a2"), (40, "ghost")])}
        affected_lookup = {"a1": "affected", "a2": "affected"}
        analysis_lookup = {"a1": "unsolved", "a2": "unsolved"}
        sample_lookup = {"a1": {"family_id": "F1"}, "a2": {"family_id": ""}}
        sample_counts, family_counts = self._counts(row, affected_lookup, analysis_lookup, sample_lookup)
        self.assertEqual(sample_counts[0], 3)
        self.assertEqual(family_counts[0], 1)


class ListAffectedUnsolvedOutliersTests(unittest.TestCase):

    def test_one_row_per_counted_sample_sorted_by_sample(self):
        row = {
            "OutlierSampleIds_AllAlleles": make_outlier_string(
                [(50, "a2"), (45, "a1"), (45, "a2"), (30, "u1"), (20, "a3")]),
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": make_outlier_string([(12, "a3")]),
        }
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected", "u1": "unaffected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}
        sample_lookup = {"a1": {"family_id": "F1"}, "a2": {"family_id": float("nan")}}
        parsed = {ot: parse_outlier_entries(row[f"OutlierSampleIds_{ot}"]) for ot in OUTLIER_TYPES}
        samples = find_affected_unsolved_samples_above_unaffected(
            row, affected_lookup, analysis_lookup, parsed)
        self.assertEqual(list_affected_unsolved_outliers("L1", samples, sample_lookup), [
            {"LocusId": "L1", "outlier_type": "AllAlleles", "sample_id": "a1",
             "family_id": "F1", "allele_size": 45},
            # a2's larger allele; its NaN family_id becomes None.
            {"LocusId": "L1", "outlier_type": "AllAlleles", "sample_id": "a2",
             "family_id": None, "allele_size": 50},
            {"LocusId": "L1", "outlier_type": "HemizygousAlleles", "sample_id": "a3",
             "family_id": None, "allele_size": 12},
        ])


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

    def test_returns_the_rows_the_counts_count_under_the_given_gate(self):
        records = [{
            "LocusId": "L",
            "OutlierSampleIds_AllAlleles": make_outlier_string([(50, "a1"), (40, "a2"), (30, "a3")]),
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": "",
            "TenK10K_99thPercentile": 35,
        }]
        sample_lookup = {"a1": {"family_id": "F1"}, "a2": {"family_id": "F1"}, "a3": {"family_id": "F2"}}
        affected_lookup = {"a1": "affected", "a2": "affected", "a3": "affected"}
        analysis_lookup = {k: "unsolved" for k in affected_lookup}

        rows = add_all_outlier_columns(records, sample_lookup, affected_lookup, analysis_lookup)
        self.assertEqual([row["sample_id"] for row in rows], ["a1", "a2", "a3"])
        self.assertEqual(records[0]["NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles"], 3)
        self.assertEqual(records[0]["NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles"], 2)

        # With TenK10K gating, a3's allele of 30 is at or below its 99th percentile.
        rows = add_all_outlier_columns(records, sample_lookup, affected_lookup, analysis_lookup,
                                       ("HPRC256", "AoU1027", "TenK10K"))
        self.assertEqual([row["sample_id"] for row in rows], ["a1", "a2"])
        self.assertEqual(records[0]["NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles"], 2)
        self.assertEqual(records[0]["NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles"], 1)


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
        self.assertEqual(records[0]["GeneTableInheritance"], "AD")
        self.assertEqual(records[1]["pLI"], 0.5)
        self.assertIsNone(records[2]["pLI"])
        self.assertEqual(records[2]["GeneTableInheritance"], "XR")
        # gene_id not in lookup -> both None.
        self.assertIsNone(records[3]["pLI"])
        self.assertIsNone(records[3]["GeneTableInheritance"])

    def test_pli_uses_whichever_values_are_present(self):
        gene_lookup = {
            "ONLY_V4": {"pLI_v2": None, "pLI_v4": 0.7},
            "ONLY_V4_NAN_V2": {"pLI_v2": float("nan"), "pLI_v4": 0.4},
            "NO_PLI_FIELDS": {"gene_symbol": "X"},
            "STRING_VALUES": {"pLI_v2": "0.2", "pLI_v4": "0.95"},
        }
        records = [{"gene_id": gene_id} for gene_id in gene_lookup]
        add_gene_columns(records, gene_lookup)
        self.assertEqual([record["pLI"] for record in records], [0.7, 0.4, None, 0.95])

    def test_no_gene_lookup_yields_null(self):
        records = [{"gene_id": "G1"}, {"gene_id": None}]
        add_gene_columns(records, None)
        for record in records:
            self.assertIsNone(record["pLI"])
            self.assertIsNone(record["GeneTableInheritance"])

    def test_blank_gene_table_inheritance_is_none(self):
        # read_gene_table gives a blank inheritance cell as a float NaN, stored as None. The
        # inheritance mode is stored only as GeneTableInheritance.
        records = [{"gene_id": "G1"}, {"gene_id": "G2"}, {"gene_id": "missing"}]
        add_gene_columns(records, {"G1": {"inheritance": "AD"}, "G2": {"inheritance": float("nan")}})
        self.assertEqual([record["GeneTableInheritance"] for record in records], ["AD", None, None])
        for record in records:
            self.assertNotIn("inheritance", record)


class OutputColumnsTests(unittest.TestCase):

    def test_exact_count_and_uniqueness(self):
        self.assertEqual(len(OUTPUT_COLUMNS), 136)
        self.assertEqual(len(set(OUTPUT_COLUMNS)), 136)

    def test_anchor_positions(self):
        # Spot-check fixed positions in OUTPUT_COLUMNS (the trailing comments give 1-based
        # column numbers).
        self.assertEqual(OUTPUT_COLUMNS[0], "LocusId")
        self.assertEqual(OUTPUT_COLUMNS[32], "pLI")  # col 33
        self.assertEqual(OUTPUT_COLUMNS[33], "FirstAffectedAlleleSize_AllAlleles")  # col 34
        self.assertEqual(OUTPUT_COLUMNS[133], "VariationClusterSizeDiff")  # col 134
        self.assertEqual(OUTPUT_COLUMNS[135], "VariationClusterFilterReason")  # col 136
        # The inheritance mode is GeneTableInheritance; no duplicate "inheritance" column.
        self.assertNotIn("inheritance", OUTPUT_COLUMNS)

    def test_each_cohort_has_median_and_mode_after_its_max_allele(self):
        for cohort in ("AoU1027", "HPRC256", "TenK10K"):
            position = OUTPUT_COLUMNS.index(f"{cohort}_MaxAllele")
            self.assertEqual(OUTPUT_COLUMNS[position + 1:position + 3],
                             [f"{cohort}_Median", f"{cohort}_Mode"], cohort)


class GeneTableColumnsTests(unittest.TestCase):
    """add_gene_columns copies a locus's gene-table row into the GeneTable* columns."""

    def test_gene_table_columns_populated(self):
        records = [{"gene_id": "ENSG1"}]
        gene_lookup = {"ENSG1": {
            "gene_symbol": "FMR1", "gene_aliases": "FRAXA", "pLI_v2": "0.9",
            "pLI_v4": "0.95", "lof_oe_ci_upper_v4": "0.3", "hgnc_gene_id": "HGNC:3",
            "inheritance": "XL", "disease_category": "ID", "LLM_phenotype_summary": "x",
            "sources": "OMIM",
        }}
        analysis_columns.add_gene_columns(records, gene_lookup)
        self.assertEqual(records[0]["GeneTableGeneSymbol"], "FMR1")
        self.assertEqual(records[0]["GeneTableInheritance"], "XL")
        self.assertEqual(records[0]["GeneTableLoeuf"], "0.3")

    def test_gene_table_columns_none_without_match(self):
        records = [{"gene_id": "missing"}]
        analysis_columns.add_gene_columns(records, {"ENSG1": {"gene_symbol": "X"}})
        # Every GeneTable* column is explicitly None so the loci schema is consistent.
        for column in analysis_columns.GENE_TABLE_COLUMN_BY_SOURCE_FIELD.values():
            self.assertIsNone(records[0][column])


if __name__ == "__main__":
    unittest.main()
