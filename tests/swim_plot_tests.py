"""Unit tests for swim_plot.generate_swim_plot_table and its helpers."""

import unittest

from swim_plot import (
    SWIM_PLOT_COLUMNS,
    SWIM_PLOT_LOCUS_COLUMNS,
    _compute_motif_category,
    _title_case_status,
    generate_swim_plot_table,
)


def _make_record(**overrides):
    """Build a minimal locus record with all swim-plot inputs, overridable."""
    record = {
        "LocusId": "1-100-110-AT",
        "Motif": "AT",
        "CanonicalMotif": "AT",
        "MotifSize": 2,
        "gene_region": "intron",
        "GeneTableGeneSymbol": "GENE1",
        "IsInMendelianGene": 0,
        "IsKnownMotif": 0,
        "gene_id": "ENSG1",
        "pLI": 0.5,
        "NumRepeatsInReference": 5,
        "HPRC256_MaxAllele": None,
        "AoU1027_MaxAllele": None,
        "TenK10K_MaxAllele": None,
        "HPRC256_99thPercentile": None,
        "AoU1027_99thPercentile": None,
        "TenK10K_99thPercentile": None,
        "HPRC256_StdevPercentile": None,
        "AoU1027_StdevPercentile": None,
        "Source": "TRails",
        "OutlierSampleIds_AllAlleles": None,
        "OutlierSampleIds_ShortAlleles": None,
        "OutlierSampleIds_HemizygousAlleles": None,
        "FirstUnaffectedAlleleSize_AllAlleles": None,
        "FirstUnaffectedAlleleSize_ShortAlleles": None,
        "FirstUnaffectedAlleleSize_HemizygousAlleles": None,
    }
    record.update(overrides)
    return record


class MotifCategoryTests(unittest.TestCase):
    def test_small_motif_size_bins(self):
        self.assertEqual(_compute_motif_category(2), "2bp")
        self.assertEqual(_compute_motif_category(24), "24bp")
        self.assertEqual(_compute_motif_category(1), "1bp")

    def test_large_motif_size_bins(self):
        self.assertEqual(_compute_motif_category(25), "25+bp")
        self.assertEqual(_compute_motif_category(100), "25+bp")

    def test_missing_motif_size(self):
        self.assertEqual(_compute_motif_category(None), "Unknown")
        self.assertEqual(_compute_motif_category(float("nan")), "Unknown")


class StatusNormalizationTests(unittest.TestCase):
    def test_affected_title_casing(self):
        self.assertEqual(_title_case_status("affected"), "Affected")
        self.assertEqual(_title_case_status("UNAFFECTED"), "Unaffected")
        self.assertEqual(_title_case_status(" possibly affected "), "Possibly Affected")

    def test_analysis_title_casing(self):
        self.assertEqual(_title_case_status("solved"), "Solved")
        self.assertEqual(_title_case_status("UNSOLVED"), "Unsolved")
        self.assertEqual(_title_case_status("partially solved"), "Partially Solved")

    def test_missing_and_blank_statuses_are_unknown(self):
        for blank in (None, float("nan"), "", "  ", "NA", "n/a", "null"):
            self.assertEqual(_title_case_status(blank), "Unknown", repr(blank))


class GenerateSwimPlotTableTests(unittest.TestCase):
    def _lookups(self):
        sample_lookup = {
            "S1": {"family_id": "F1", "sex": "male", "affected_status": "affected",
                   "phenotype_description": "seizures"},
            "S2": {"family_id": "F2", "sex": "female", "affected_status": "unaffected",
                   "phenotype_description": None},
            "S3": {"family_id": "F1", "sex": "male", "affected_status": "affected",
                   "phenotype_description": "ataxia"},
        }
        analysis_lookup = {"S1": "unsolved", "S2": "unaffected", "S3": "solved"}
        return sample_lookup, analysis_lookup

    def test_one_row_per_outlier_entry(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [
            _make_record(OutlierSampleIds_AllAlleles="40x:S1,30x:S2,20x:S3"),
        ]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(len(rows), 3)
        self.assertEqual([r["sample_id"] for r in rows], ["S1", "S2", "S3"])

    def test_rank_ordering_largest_is_rank_one(self):
        sample_lookup, analysis_lookup = self._lookups()
        # Provide entries out of order; parse sorts descending, so 40 -> rank 1.
        records = [_make_record(OutlierSampleIds_AllAlleles="20x:S3,40x:S1,30x:S2")]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(rows[0]["allele_size"], 40)
        self.assertEqual(rows[0]["outlier_rank"], 1)
        self.assertEqual(rows[1]["allele_size"], 30)
        self.assertEqual(rows[1]["outlier_rank"], 2)
        self.assertEqual(rows[2]["allele_size"], 20)
        self.assertEqual(rows[2]["outlier_rank"], 3)

    def test_row_keys_are_the_table_columns(self):
        # result_database.write_swim_plot declares the table from SWIM_PLOT_COLUMNS.
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(OutlierSampleIds_AllAlleles="40x:S1")]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(list(rows[0]), SWIM_PLOT_COLUMNS)

    def test_affected_and_analysis_title_cased_on_rows(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(OutlierSampleIds_AllAlleles="40x:S1,30x:S2")]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(rows[0]["affected_status"], "Affected")
        self.assertEqual(rows[0]["analysis_status"], "Unsolved")
        self.assertEqual(rows[1]["affected_status"], "Unaffected")
        self.assertEqual(rows[1]["analysis_status"], "Unaffected")

    def test_possibly_affected_preserved_from_sample_row(self):
        # sample_lookup rows carry the affected_status as entered, while read_sample_metadata's
        # affected_lookup collapses "possibly affected" to "affected" for the analysis logic. The
        # swim-plot row shows the entered "Possibly Affected".
        sample_lookup, analysis_lookup = self._lookups()
        sample_lookup["S1"]["affected_status"] = "possibly affected"
        records = [_make_record(OutlierSampleIds_AllAlleles="40x:S1")]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(rows[0]["affected_status"], "Possibly Affected")

    def test_is_above_first_unaffected_flag(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(
            OutlierSampleIds_AllAlleles="40x:S1,30x:S2,25x:S3",
            FirstUnaffectedAlleleSize_AllAlleles=30,
        )]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        by_sample = {r["sample_id"]: r for r in rows}
        self.assertEqual(by_sample["S1"]["is_above_first_unaffected"], 1)  # 40 > 30
        self.assertEqual(by_sample["S2"]["is_above_first_unaffected"], 0)  # 30 not > 30
        self.assertEqual(by_sample["S3"]["is_above_first_unaffected"], 0)  # 25 not > 30
        self.assertEqual(by_sample["S1"]["FirstUnaffectedAlleleSize"], 30)

    def test_is_above_first_unaffected_when_no_unaffected(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(
            OutlierSampleIds_AllAlleles="40x:S1",
            FirstUnaffectedAlleleSize_AllAlleles=None,
        )]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(rows[0]["is_above_first_unaffected"], 1)  # >0 with no baseline
        self.assertIsNone(rows[0]["FirstUnaffectedAlleleSize"])

    def test_motif_category_bins(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [
            _make_record(LocusId="L_small", MotifSize=3, OutlierSampleIds_AllAlleles="40x:S1"),
            _make_record(LocusId="L_big", MotifSize=30, OutlierSampleIds_AllAlleles="40x:S1"),
            _make_record(LocusId="L_unknown", MotifSize=None, OutlierSampleIds_AllAlleles="40x:S1"),
        ]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        by_locus = {r["LocusId"]: r["motif_category"] for r in rows}
        self.assertEqual(by_locus["L_small"], "3bp")
        self.assertEqual(by_locus["L_big"], "25+bp")
        self.assertEqual(by_locus["L_unknown"], "Unknown")

    def test_all_three_outlier_types_emitted(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(
            OutlierSampleIds_AllAlleles="40x:S1",
            OutlierSampleIds_ShortAlleles="20x:S2",
            OutlierSampleIds_HemizygousAlleles="15x:S3",
        )]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(len(rows), 3)
        self.assertEqual({r["outlier_type"] for r in rows},
                         {"AllAlleles", "ShortAlleles", "HemizygousAlleles"})

    def test_missing_outlier_value_skipped(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(OutlierSampleIds_AllAlleles=".")]
        self.assertEqual(generate_swim_plot_table(records, sample_lookup, analysis_lookup), [])

    def test_sample_missing_metadata_treated_as_none(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(OutlierSampleIds_AllAlleles="40x:UNKNOWN_SAMPLE")]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertEqual(len(rows), 1)
        self.assertIsNone(rows[0]["family_id"])
        self.assertIsNone(rows[0]["sex"])
        self.assertIsNone(rows[0]["phenotype_description"])
        self.assertEqual(rows[0]["affected_status"], "Unknown")
        self.assertEqual(rows[0]["analysis_status"], "Unknown")

    def test_locus_passthrough_columns_present(self):
        sample_lookup, analysis_lookup = self._lookups()
        records = [_make_record(OutlierSampleIds_AllAlleles="40x:S1")]
        row = generate_swim_plot_table(records, sample_lookup, analysis_lookup)[0]
        for column in SWIM_PLOT_LOCUS_COLUMNS:
            self.assertIn(column, row)
        self.assertEqual(row["SourceDb"], "TRails")
        self.assertEqual(row["GeneTableGeneSymbol"], "GENE1")
        self.assertEqual(row["LocusId"], "1-100-110-AT")

    def test_phenotype_nan_becomes_none(self):
        sample_lookup, analysis_lookup = self._lookups()
        sample_lookup["S9"] = {"family_id": "F9", "sex": "female", "phenotype_description": float("nan")}
        records = [_make_record(OutlierSampleIds_AllAlleles="40x:S9")]
        rows = generate_swim_plot_table(records, sample_lookup, analysis_lookup)
        self.assertIsNone(rows[0]["phenotype_description"])

    def test_blank_metadata_cells_become_none(self):
        # read_sample_metadata keeps an empty cell as '' and an "NA" cell as "NA", which must be
        # stored as NULL like the family_id of the affected_unsolved_outliers_above_unaffected
        # table.
        sample_lookup, analysis_lookup = self._lookups()
        sample_lookup["S9"] = {"family_id": "", "sex": " ", "phenotype_description": ""}
        sample_lookup["S8"] = {"family_id": "NA", "sex": "NA", "phenotype_description": "null"}
        records = [_make_record(OutlierSampleIds_AllAlleles="40x:S9,30x:S8")]
        for row in generate_swim_plot_table(records, sample_lookup, analysis_lookup):
            self.assertIsNone(row["family_id"], row["sample_id"])
            self.assertIsNone(row["sex"], row["sample_id"])
            self.assertIsNone(row["phenotype_description"], row["sample_id"])


if __name__ == "__main__":
    unittest.main()
