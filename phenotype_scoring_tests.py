"""Unit tests for phenotype_scoring.

All tests force the Jaccard fallback by monkeypatching
``phenotype_scoring.ensure_ontology`` to return False, so no pyhpo ontology is
built and the tests are deterministic regardless of whether pyhpo is installed.
"""

import unittest

import phenotype_scoring
from phenotype_scoring import (
    OUTLIER_TYPE_INHERITANCE,
    compute_combined_ic,
    compute_pairwise_shared_counts,
    compute_phenotype_scores,
    compute_similarity_pyhpo,
    get_qualifying_samples,
    score_patient_vs_gene_filtered,
)


class ForceJaccardTestCase(unittest.TestCase):
    """Base case that disables pyhpo so similarity uses the Jaccard fallback."""

    def setUp(self):
        self._saved_ensure_ontology = phenotype_scoring.ensure_ontology
        phenotype_scoring.ensure_ontology = lambda: False

    def tearDown(self):
        phenotype_scoring.ensure_ontology = self._saved_ensure_ontology


class SimilarityFallbackTests(ForceJaccardTestCase):
    def test_jaccard_similarity(self):
        # intersection {A,B}=2, union {A,B,C,D}=4 -> 0.5
        self.assertEqual(
            compute_similarity_pyhpo(["HP:A", "HP:B", "HP:C"], ["HP:A", "HP:B", "HP:D"]),
            0.5,
        )

    def test_jaccard_identical(self):
        self.assertEqual(compute_similarity_pyhpo(["HP:A"], ["HP:A"]), 1.0)

    def test_jaccard_disjoint(self):
        self.assertEqual(compute_similarity_pyhpo(["HP:A"], ["HP:B"]), 0.0)

    def test_jaccard_empty(self):
        self.assertEqual(compute_similarity_pyhpo([], ["HP:A"]), 0.0)
        self.assertEqual(compute_similarity_pyhpo([], []), 0.0)

    def test_shared_counts_flat_weight(self):
        # Without an ontology the IC weight is flat 1.0 per shared term, so the
        # IC-weighted count equals the raw count.
        raw_count, ic_count = compute_pairwise_shared_counts(
            ["HP:A", "HP:B", "HP:C"], ["HP:B", "HP:C", "HP:D"]
        )
        self.assertEqual(raw_count, 2)
        self.assertEqual(ic_count, 2.0)

    def test_shared_counts_none(self):
        self.assertEqual(compute_pairwise_shared_counts(["HP:A"], ["HP:B"]), (0, 0.0))


class CombinedICTests(unittest.TestCase):
    def test_average_when_both_positive(self):
        self.assertEqual(compute_combined_ic(_FakeTerm(2.0, 4.0)), 3.0)

    def test_single_positive(self):
        self.assertEqual(compute_combined_ic(_FakeTerm(5.0, 0.0)), 5.0)
        self.assertEqual(compute_combined_ic(_FakeTerm(0.0, 7.0)), 7.0)

    def test_default_weight_when_zero(self):
        self.assertEqual(compute_combined_ic(_FakeTerm(0.0, 0.0)), 1.0)

    def test_default_weight_on_error(self):
        self.assertEqual(compute_combined_ic(object()), 1.0)


class _FakeIC:
    def __init__(self, omim, gene):
        self.omim = omim
        self.gene = gene


class _FakeTerm:
    def __init__(self, omim, gene):
        self.information_content = _FakeIC(omim, gene)


# A tiny gene-disease catalog: gene GENEX has two diseases, one autosomal
# dominant (matches AllAlleles) and one autosomal recessive (matches
# ShortAlleles). Phenotype terms are deliberately disjoint between the two so
# the inheritance filter changes which disease is reachable.
GENE_DISEASE_DATA = {
    "GENEX": {
        "OMIM:1": {"phenotypes": {"HP:A", "HP:B"}, "inheritance": {"AD"}},
        "OMIM:2": {"phenotypes": {"HP:X", "HP:Y"}, "inheritance": {"AR"}},
    }
}


class ScorePatientVsGeneTests(ForceJaccardTestCase):
    def test_unknown_gene_has_no_annotation(self):
        result = score_patient_vs_gene_filtered({"HP:A"}, "NOPE", GENE_DISEASE_DATA, {"AD"})
        self.assertFalse(result["has_annotation"])
        self.assertEqual(result["n_matching_diseases"], 0)
        self.assertEqual(result["best_similarity"], 0.0)

    def test_inheritance_filter_selects_dominant_disease(self):
        # AD allowed -> only OMIM:1 (HP:A,HP:B) reachable; patient {HP:A,HP:B} -> 1.0
        result = score_patient_vs_gene_filtered({"HP:A", "HP:B"}, "GENEX", GENE_DISEASE_DATA, {"AD"})
        self.assertTrue(result["has_annotation"])
        self.assertEqual(result["n_matching_diseases"], 1)
        self.assertEqual(result["best_disease"], "OMIM:1")
        self.assertEqual(result["best_similarity"], 1.0)
        self.assertEqual(result["best_inheritance"], "AD")
        self.assertEqual(result["overlap_count"], 2)

    def test_inheritance_filter_excludes_dominant_for_recessive_outlier(self):
        # AR allowed -> only OMIM:2 (HP:X,HP:Y) reachable; patient {HP:A,HP:B} -> 0.0
        result = score_patient_vs_gene_filtered({"HP:A", "HP:B"}, "GENEX", GENE_DISEASE_DATA, {"AR"})
        self.assertEqual(result["best_disease"], "OMIM:2")
        self.assertEqual(result["best_similarity"], 0.0)
        self.assertEqual(result["n_matching_diseases"], 1)

    def test_annotated_gene_no_matching_inheritance(self):
        result = score_patient_vs_gene_filtered({"HP:A"}, "GENEX", GENE_DISEASE_DATA, {"XR"})
        self.assertTrue(result["has_annotation"])
        self.assertEqual(result["n_matching_diseases"], 0)
        self.assertIsNone(result["best_disease"])

    def test_outlier_type_inheritance_map(self):
        self.assertEqual(OUTLIER_TYPE_INHERITANCE["AllAlleles"], {"AD", "XD"})
        self.assertEqual(OUTLIER_TYPE_INHERITANCE["ShortAlleles"], {"AR"})
        self.assertEqual(OUTLIER_TYPE_INHERITANCE["HemizygousAlleles"], {"XR", "XL", "XD"})


class GetQualifyingSamplesTests(ForceJaccardTestCase):
    def _base_row(self, outlier_value):
        # No unaffected sample and no population percentiles -> gates 2 and 3
        # depend only on allele_size > 0 and "no cohort data -> include".
        return {
            "LocusId": "chr1-100-110-AC",
            "OutlierSampleIds_AllAlleles": outlier_value,
        }

    def test_gate_affected_unsolved(self):
        row = self._base_row("40x:s_affected,30x:s_unaffected")
        affected = {"s_affected": "affected", "s_unaffected": "unaffected"}
        analysis = {"s_affected": "unsolved", "s_unaffected": "unaffected"}
        hpo = {"s_affected": {"HP:A"}, "s_unaffected": {"HP:B"}}
        qualifying = get_qualifying_samples(row, "AllAlleles", affected, analysis, hpo)
        self.assertEqual([entry[0] for entry in qualifying], ["s_affected"])

    def test_gate_requires_hpo(self):
        row = self._base_row("40x:s1,30x:s2")
        affected = {"s1": "affected", "s2": "affected"}
        analysis = {"s1": "unsolved", "s2": "unsolved"}
        hpo = {"s1": {"HP:A"}}  # s2 has no HPO terms
        qualifying = get_qualifying_samples(row, "AllAlleles", affected, analysis, hpo)
        self.assertEqual([entry[0] for entry in qualifying], ["s1"])

    def test_gate_above_unaffected(self):
        row = self._base_row("40x:s_big,20x:s_small")
        row["FirstUnaffectedAlleleSize_AllAlleles"] = 30
        affected = {"s_big": "affected", "s_small": "affected"}
        analysis = {"s_big": "unsolved", "s_small": "unsolved"}
        hpo = {"s_big": {"HP:A"}, "s_small": {"HP:B"}}
        qualifying = get_qualifying_samples(row, "AllAlleles", affected, analysis, hpo)
        # only the 40 allele exceeds the unaffected threshold of 30
        self.assertEqual([entry[0] for entry in qualifying], ["s_big"])

    def test_gate_above_population_p99(self):
        row = self._base_row("40x:s_big,30x:s_small")
        row["HPRC256_99thPercentile"] = 35
        affected = {"s_big": "affected", "s_small": "affected"}
        analysis = {"s_big": "unsolved", "s_small": "unsolved"}
        hpo = {"s_big": {"HP:A"}, "s_small": {"HP:B"}}
        qualifying = get_qualifying_samples(row, "AllAlleles", affected, analysis, hpo)
        self.assertEqual([entry[0] for entry in qualifying], ["s_big"])

    def test_descending_order_and_dedup_keep_largest(self):
        # s_dup appears at two allele sizes (50 and 35); dedup keeps the 50 entry.
        row = self._base_row("50x:s_dup,45x:s_other,35x:s_dup")
        affected = {"s_dup": "affected", "s_other": "affected"}
        analysis = {"s_dup": "unsolved", "s_other": "unsolved"}
        hpo = {"s_dup": {"HP:A"}, "s_other": {"HP:B"}}
        qualifying = get_qualifying_samples(row, "AllAlleles", affected, analysis, hpo)
        self.assertEqual([(e[0], e[1]) for e in qualifying], [("s_dup", 50), ("s_other", 45)])

    def test_missing_outlier_value_returns_empty(self):
        self.assertEqual(get_qualifying_samples({"OutlierSampleIds_AllAlleles": ""},
                                                "AllAlleles", {}, {}, {"x": {"HP:A"}}), [])


class ComputePhenotypeScoresTests(ForceJaccardTestCase):
    def _records(self):
        # One locus in GENEX with three qualifying AllAlleles outliers.
        # HPO overlaps (Jaccard): s1{A,B} vs s2{A,B} -> 1.0 ; s2{A,B} vs s3{C} -> 0.0
        return [{
            "LocusId": "chr1-100-110-AC",
            "gene_id": "ENSG_X",
            "IsKnownMotif": 1,
            "gene_region": "CDS",
            "gene_region_rank": 1,
            "MotifSize": 2,
            "NumRepeatsInReference": 5,
            "FirstAffectedAlleleSize_AllAlleles": 50,
            "FirstUnaffectedAlleleSize_AllAlleles": None,
            "NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles": 3,
            "NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles": 2,
            "OutlierSampleIds_AllAlleles": "50x:s1,45x:s2,40x:s3",
            "OutlierSampleIds_ShortAlleles": "",
            "OutlierSampleIds_HemizygousAlleles": "",
        }]

    def _lookups(self):
        affected = {"s1": "affected", "s2": "affected", "s3": "affected"}
        analysis = {"s1": "unsolved", "s2": "unsolved", "s3": "unsolved"}
        return affected, analysis

    def test_empty_phenotypes_returns_empty(self):
        records = self._records()
        affected, analysis = self._lookups()
        per_outlier, per_locus = compute_phenotype_scores(
            records, {}, {}, {}, affected, analysis
        )
        self.assertEqual(per_outlier, [])
        self.assertEqual(per_locus, [])
        # No phenotype columns added to records.
        self.assertNotIn("MaxGenePhenoSim_AllAlleles", records[0])

    def test_per_locus_aggregation(self):
        records = self._records()
        affected, analysis = self._lookups()
        participant_to_hpo = {"s1": {"HP:A", "HP:B"}, "s2": {"HP:A", "HP:B"}, "s3": {"HP:C"}}
        gene_lookup = {"ENSG_X": {"gene_symbol": "GENEX"}}

        per_outlier, per_locus = compute_phenotype_scores(
            records, participant_to_hpo, gene_lookup, GENE_DISEASE_DATA, affected, analysis
        )

        all_alleles_locus_rows = [r for r in per_locus if r["outlier_type"] == "AllAlleles"]
        self.assertEqual(len(all_alleles_locus_rows), 1)
        locus_row = all_alleles_locus_rows[0]

        self.assertEqual(locus_row["num_qualifying_samples"], 3)
        self.assertEqual(locus_row["qualifying_sample_ids"], "s1,s2,s3")
        # pairwise: s1 vs s2 -> 1.0, s2 vs s3 -> 0.0 ; sum = 1.0
        self.assertEqual(locus_row["sum_pairwise_similarity"], 1.0)
        # gene-sim (AD diseases only): s1{A,B}->1.0, s2{A,B}->1.0, s3{C}->0.0 ; max=1.0
        self.assertEqual(locus_row["max_gene_phenotype_similarity"], 1.0)
        # shared raw: s1&s2={A,B}=2, s2&s3={}=0 ; sum=2
        self.assertEqual(locus_row["sum_pairwise_shared_raw"], 2)
        self.assertEqual(locus_row["sum_pairwise_shared_ic"], 2.0)
        # denormalized passthroughs
        self.assertEqual(locus_row["IsKnownMotif"], 1)
        self.assertEqual(locus_row["FirstAffectedAlleleSize"], 50)
        self.assertEqual(locus_row["NumAffectedAboveUnaffected"], 3)
        self.assertEqual(locus_row["NumAffectedFamiliesAboveUnaffected"], 2)

        # Records augmented with denormalized phenotype columns.
        self.assertEqual(records[0]["MaxGenePhenoSim_AllAlleles"], 1.0)
        self.assertEqual(records[0]["SumPairwiseSim_AllAlleles"], 1.0)

    def test_per_outlier_rows_structure(self):
        records = self._records()
        affected, analysis = self._lookups()
        participant_to_hpo = {"s1": {"HP:A", "HP:B"}, "s2": {"HP:A", "HP:B"}, "s3": {"HP:C"}}
        gene_lookup = {"ENSG_X": {"gene_symbol": "GENEX"}}

        per_outlier, _ = compute_phenotype_scores(
            records, participant_to_hpo, gene_lookup, GENE_DISEASE_DATA, affected, analysis
        )
        all_alleles_rows = [r for r in per_outlier if r["outlier_type"] == "AllAlleles"]
        self.assertEqual([r["sample_id"] for r in all_alleles_rows], ["s1", "s2", "s3"])
        # Last sample has no "next" sample -> pairwise fields are None.
        self.assertIsNone(all_alleles_rows[-1]["pairwise_similarity_to_next"])
        self.assertIsNone(all_alleles_rows[-1]["next_sample_id"])
        # First sample's next is s2.
        self.assertEqual(all_alleles_rows[0]["next_sample_id"], "s2")
        self.assertEqual(all_alleles_rows[0]["pairwise_similarity_to_next"], 1.0)
        # gene_symbol resolved via gene_id -> gene_lookup.
        self.assertEqual(all_alleles_rows[0]["gene_symbol"], "GENEX")

    def test_no_gene_symbol_leaves_gene_sim_none(self):
        records = self._records()
        records[0]["gene_id"] = None
        affected, analysis = self._lookups()
        participant_to_hpo = {"s1": {"HP:A", "HP:B"}, "s2": {"HP:A", "HP:B"}, "s3": {"HP:C"}}

        per_outlier, per_locus = compute_phenotype_scores(
            records, participant_to_hpo, {}, GENE_DISEASE_DATA, affected, analysis
        )
        all_alleles_rows = [r for r in per_outlier if r["outlier_type"] == "AllAlleles"]
        self.assertTrue(all(r["gene_phenotype_similarity"] is None for r in all_alleles_rows))
        locus_row = [r for r in per_locus if r["outlier_type"] == "AllAlleles"][0]
        self.assertIsNone(locus_row["max_gene_phenotype_similarity"])
        # pairwise still computed.
        self.assertEqual(locus_row["sum_pairwise_similarity"], 1.0)


if __name__ == "__main__":
    unittest.main()
