"""Unit tests for phenotype_scoring.

The ForceJaccardTestCase subclasses monkeypatch ``phenotype_scoring.ensure_ontology``
to return False, so they use the Jaccard fallback and give the same results whether
or not pyhpo is installed. CombinedICTests uses fake terms and needs no ontology.
PyhpoUnknownTermTests builds and uses the real pyhpo ontology, and is skipped when
pyhpo is not installed.
"""

import unittest

import analysis_columns
import input_tables
import phenotype_scoring
from phenotype_scoring import (
    OUTLIER_TYPE_INHERITANCE,
    allowed_inheritance_for_locus,
    compute_combined_ic,
    compute_pairwise_shared_counts,
    compute_phenotype_scores,
    compute_similarity_pyhpo,
    get_qualifying_samples,
    resolvable_hpo_terms,
    score_patient_vs_gene_filtered,
)
from tests.temp_file_helpers import write_temporary_file


def pyhpo_is_installed():
    """Return True if the optional pyhpo dependency can be imported."""
    try:
        import pyhpo  # noqa: F401
    except ImportError:
        return False
    return True


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
    def test_inheritance_filter_selects_dominant_disease(self):
        # AD allowed -> only OMIM:1 (HP:A,HP:B) reachable; patient {HP:A,HP:B} -> 1.0
        result = score_patient_vs_gene_filtered({"HP:A", "HP:B"}, "GENEX", GENE_DISEASE_DATA, {"AD"})
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
        self.assertEqual(result["n_matching_diseases"], 0)
        self.assertIsNone(result["best_disease"])

    def test_outlier_type_inheritance_map(self):
        self.assertEqual(OUTLIER_TYPE_INHERITANCE["AllAlleles"], {"AD", "XD", "MT", "YL"})
        self.assertEqual(OUTLIER_TYPE_INHERITANCE["ShortAlleles"], {"AR"})
        self.assertEqual(OUTLIER_TYPE_INHERITANCE["HemizygousAlleles"], {"XR", "XL", "XD", "MT", "YL"})

    def test_mitochondrial_and_y_linked_codes_are_allowed_only_on_their_chromosome(self):
        self.assertEqual(allowed_inheritance_for_locus("AllAlleles", "chr1"), {"AD", "XD"})
        self.assertEqual(allowed_inheritance_for_locus("AllAlleles", None), {"AD", "XD"})
        self.assertEqual(allowed_inheritance_for_locus("AllAlleles", "chrY"), {"AD", "XD", "YL"})
        self.assertEqual(allowed_inheritance_for_locus("HemizygousAlleles", "chrM"),
                         {"XR", "XL", "XD", "MT"})
        self.assertEqual(allowed_inheritance_for_locus("HemizygousAlleles", "chrMT"),
                         {"XR", "XL", "XD", "MT"})
        self.assertEqual(allowed_inheritance_for_locus("ShortAlleles", "chrY"), {"AR"})


class GetQualifyingSamplesTests(ForceJaccardTestCase):
    def _base_row(self, outlier_value):
        # No population percentiles, so only the affected/unaffected statuses and the HPO
        # terms decide which samples qualify.
        return {
            "LocusId": "chr1-100-110-AC",
            "OutlierSampleIds_AllAlleles": outlier_value,
        }

    def test_gate_affected_unsolved(self):
        row = self._base_row("40x:s_affected,30x:s_unaffected")
        affected = {"s_affected": "affected", "s_unaffected": "unaffected"}
        analysis = {"s_affected": "unsolved", "s_unaffected": "unaffected"}
        hpo = {"s_affected": {"HP:A"}, "s_unaffected": {"HP:B"}}
        qualifying = get_qualifying_samples(row, affected, analysis, hpo)["AllAlleles"]
        self.assertEqual([entry[0] for entry in qualifying], ["s_affected"])

    def test_gate_requires_hpo(self):
        row = self._base_row("40x:s1,30x:s2")
        affected = {"s1": "affected", "s2": "affected"}
        analysis = {"s1": "unsolved", "s2": "unsolved"}
        hpo = {"s1": {"HP:A"}}  # s2 has no HPO terms
        qualifying = get_qualifying_samples(row, affected, analysis, hpo)["AllAlleles"]
        self.assertEqual([entry[0] for entry in qualifying], ["s1"])

    def test_gate_above_unaffected(self):
        row = self._base_row("40x:s_big,30x:s_unaffected,20x:s_small")
        affected = {"s_big": "affected", "s_unaffected": "unaffected", "s_small": "affected"}
        analysis = {"s_big": "unsolved", "s_unaffected": "unaffected", "s_small": "unsolved"}
        hpo = {"s_big": {"HP:A"}, "s_small": {"HP:B"}}
        qualifying = get_qualifying_samples(row, affected, analysis, hpo)["AllAlleles"]
        # only the 40 allele exceeds the unaffected sample's 30
        self.assertEqual([entry[0] for entry in qualifying], ["s_big"])

    def test_gate_above_population_p99(self):
        row = self._base_row("40x:s_big,30x:s_small")
        row["HPRC256_99thPercentile"] = 35
        affected = {"s_big": "affected", "s_small": "affected"}
        analysis = {"s_big": "unsolved", "s_small": "unsolved"}
        hpo = {"s_big": {"HP:A"}, "s_small": {"HP:B"}}
        qualifying = get_qualifying_samples(row, affected, analysis, hpo)["AllAlleles"]
        self.assertEqual([entry[0] for entry in qualifying], ["s_big"])

    def test_gate_uses_the_population_gate_cohorts(self):
        # TenK10K's 99th percentile only gates when it is one of the gate cohorts.
        row = self._base_row("40x:s_big,30x:s_small")
        row["TenK10K_99thPercentile"] = 35
        affected = {"s_big": "affected", "s_small": "affected"}
        analysis = {"s_big": "unsolved", "s_small": "unsolved"}
        hpo = {"s_big": {"HP:A"}, "s_small": {"HP:B"}}
        qualifying = get_qualifying_samples(row, affected, analysis, hpo)["AllAlleles"]
        self.assertEqual([entry[0] for entry in qualifying], ["s_big", "s_small"])
        qualifying = get_qualifying_samples(row, affected, analysis, hpo,
                                            ("HPRC256", "AoU1027", "TenK10K"))["AllAlleles"]
        self.assertEqual([entry[0] for entry in qualifying], ["s_big"])

    def test_descending_order_and_dedup_keep_largest(self):
        # s_dup appears at two allele sizes (50 and 35); dedup keeps the 50 entry.
        row = self._base_row("50x:s_dup,45x:s_other,35x:s_dup")
        affected = {"s_dup": "affected", "s_other": "affected"}
        analysis = {"s_dup": "unsolved", "s_other": "unsolved"}
        hpo = {"s_dup": {"HP:A"}, "s_other": {"HP:B"}}
        qualifying = get_qualifying_samples(row, affected, analysis, hpo)["AllAlleles"]
        self.assertEqual([(e[0], e[1]) for e in qualifying], [("s_dup", 50), ("s_other", 45)])

    def test_missing_outlier_value_returns_empty(self):
        self.assertEqual(get_qualifying_samples({"OutlierSampleIds_AllAlleles": ""},
                                                {}, {}, {"x": {"HP:A"}}),
                         {"AllAlleles": [], "ShortAlleles": [], "HemizygousAlleles": []})

    def test_qualifying_samples_are_the_counted_samples_that_have_hpo_terms(self):
        # The phenotype scores and the NumAffectedUnsolved* counts must select the same
        # samples: a tie with the unaffected sample is not above it, the population gate
        # applies, and only the HPO filter is specific to phenotype scoring.
        record = {
            "LocusId": "chr1-100-110-AC",
            "OutlierSampleIds_AllAlleles": "60x:a1,50x:a2,45x:a3,40x:a4,40x:u1,30x:a5",
            "OutlierSampleIds_ShortAlleles": "45x:a3,40x:a4",
            "OutlierSampleIds_HemizygousAlleles": "",
            "AoU1027_99thPercentile": 44,
        }
        affected = {"a1": "affected", "a2": "unknown", "a3": "affected", "a4": "affected",
                    "a5": "affected", "u1": "unaffected"}
        analysis = {"a1": "unsolved", "a2": "unknown", "a3": "partially solved",
                    "a4": "unsolved", "a5": "unsolved", "u1": "unsolved"}
        hpo = {"a1": {"HP:A"}, "a2": {"HP:B"}, "a4": {"HP:C"}, "a5": {"HP:D"}}
        counted_rows = analysis_columns.add_all_outlier_columns(
            [dict(record)], {}, affected, analysis)
        qualifying = get_qualifying_samples(record, affected, analysis, hpo)
        for outlier_type in analysis_columns.OUTLIER_TYPES:
            counted_sample_ids_with_hpo_terms = {
                row["sample_id"] for row in counted_rows
                if row["outlier_type"] == outlier_type and row["sample_id"] in hpo}
            self.assertEqual({entry[0] for entry in qualifying[outlier_type]},
                             counted_sample_ids_with_hpo_terms, outlier_type)
        self.assertEqual([entry[0] for entry in qualifying["AllAlleles"]], ["a1", "a2"])
        self.assertEqual(qualifying["ShortAlleles"], [])


class ComputePhenotypeScoresTests(ForceJaccardTestCase):
    def _records(self):
        # One locus in GENEX with three qualifying AllAlleles outliers. The scorer reads only
        # these keys: it derives the qualifying samples from the outlier strings itself.
        # HPO overlaps (Jaccard): s1{A,B} vs s2{A,B} -> 1.0 ; s2{A,B} vs s3{C} -> 0.0
        return [{
            "LocusId": "chr1-100-110-AC",
            "gene_id": "ENSG_X",
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
        # The keys are absent rather than present-and-None, which is why the loci writer
        # (whose schema comes from the record keys) omits the columns entirely.
        self.assertNotIn("MaxGenePhenoSim_AllAlleles", records[0])
        self.assertNotIn("SumPairwiseSim_AllAlleles", records[0])

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
        # The row carries only the score columns the server selects; annotations such as
        # IsKnownMotif and MotifSize live in the loci table under the same LocusId.
        self.assertEqual(sorted(locus_row), [
            "locus_id", "max_gene_phenotype_similarity", "num_qualifying_samples",
            "outlier_type", "qualifying_sample_ids", "sum_pairwise_shared_ic",
            "sum_pairwise_shared_raw", "sum_pairwise_similarity",
        ])

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

    def test_mixed_case_gene_symbol_is_matched_case_insensitively(self):
        # HPO's genes_to_phenotype.txt writes the HGNC symbol "C9orf72", while genelookup's gene
        # table upper-cases it to "C9ORF72"; the ALS/FTD repeat locus must still be scored.
        gene_disease_data = input_tables.read_gene_disease_phenotypes(write_temporary_file(
            self,
            "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n"
            "203228\tC9orf72\tHP:0000006\tAutosomal dominant inheritance\t-\tOMIM:105550\n"
            "203228\tC9orf72\tHP:A\tphenotype A\t-\tOMIM:105550\n"
            "203228\tC9orf72\tHP:B\tphenotype B\t-\tOMIM:105550\n",
            suffix=".txt"))
        affected, analysis = self._lookups()
        participant_to_hpo = {"s1": {"HP:A", "HP:B"}, "s2": {"HP:A", "HP:B"}, "s3": {"HP:C"}}

        for gene_table_symbol in ("C9ORF72", "C9orf72"):
            records = self._records()
            per_outlier, _ = compute_phenotype_scores(
                records, participant_to_hpo, {"ENSG_X": {"gene_symbol": gene_table_symbol}},
                gene_disease_data, affected, analysis)
            all_alleles_rows = [r for r in per_outlier if r["outlier_type"] == "AllAlleles"]
            self.assertEqual([r["gene_phenotype_similarity"] for r in all_alleles_rows],
                             [1.0, 1.0, 0.0], gene_table_symbol)
            self.assertEqual(all_alleles_rows[0]["best_matching_disease"], "OMIM:105550")
            # The rows keep the gene table's spelling of the symbol.
            self.assertEqual(all_alleles_rows[0]["gene_symbol"], gene_table_symbol)
            self.assertEqual(records[0]["MaxGenePhenoSim_AllAlleles"], 1.0)

    def test_mitochondrial_and_y_linked_diseases_count_only_on_their_chromosome(self):
        # KDM5D's only HPO disease is Y-linked (OMIM:415000) and MT-ND5's is mitochondrial
        # (OMIM:540000). A single-allele chrY or chrM call is an outlier in all three types, and
        # the disease counts for AllAlleles and HemizygousAlleles, as an X-linked dominant one
        # does. HPO also tags SURF1's autosomal recessive Leigh syndrome (OMIM:220110) as
        # mitochondrial; on that chr9 gene only the AR code counts.
        gene_disease_data = {
            "KDM5D": {"OMIM:415000": {"phenotypes": {"HP:A", "HP:B"}, "inheritance": {"YL"}}},
            "MT-ND5": {"OMIM:540000": {"phenotypes": {"HP:A", "HP:B"}, "inheritance": {"MT"}}},
            "SURF1": {"OMIM:220110": {"phenotypes": {"HP:A", "HP:B"}, "inheritance": {"AR", "MT"}}},
        }
        gene_lookup = {gene_symbol: {"gene_symbol": gene_symbol} for gene_symbol in gene_disease_data}
        records = [{"LocusId": f"{chrom}-100-110-AC", "Chrom": chrom, "gene_id": gene_symbol,
                    **{f"OutlierSampleIds_{outlier_type}": "50x:s1"
                       for outlier_type in analysis_columns.OUTLIER_TYPES}}
                   for chrom, gene_symbol in [("chrY", "KDM5D"), ("chrM", "MT-ND5"), ("chr9", "SURF1")]]

        per_outlier, _ = compute_phenotype_scores(
            records, {"s1": {"HP:A", "HP:B"}}, gene_lookup, gene_disease_data,
            {"s1": "affected"}, {"s1": "unsolved"})

        self.assertEqual(
            {(row["gene_symbol"], row["outlier_type"]): row["best_disease_inheritance"]
             for row in per_outlier},
            {("KDM5D", "AllAlleles"): "YL", ("KDM5D", "ShortAlleles"): None,
             ("KDM5D", "HemizygousAlleles"): "YL",
             ("MT-ND5", "AllAlleles"): "MT", ("MT-ND5", "ShortAlleles"): None,
             ("MT-ND5", "HemizygousAlleles"): "MT",
             ("SURF1", "AllAlleles"): None, ("SURF1", "ShortAlleles"): "AR,MT",
             ("SURF1", "HemizygousAlleles"): None})


@unittest.skipUnless(pyhpo_is_installed(), "pyhpo is not installed")
class PyhpoUnknownTermTests(unittest.TestCase):
    """With pyhpo installed, an unknown HPO id must not zero a whole score.

    pyhpo's HPOSet.from_queries raises for the entire set when any single term is
    unknown, so a retired or misspelled id is dropped on its own and the remaining
    terms are still scored. These tests use the real ontology (HP:0001250 Seizure,
    HP:0001263 Global developmental delay) and skip when pyhpo is not installed.
    """

    @classmethod
    def setUpClass(cls):
        # resolvable_hpo_terms assumes ensure_ontology has already succeeded, so build the
        # ontology here rather than depending on an earlier test having built it.
        if not phenotype_scoring.ensure_ontology():
            raise unittest.SkipTest("the pyhpo ontology could not be initialized")

    def test_unknown_term_is_dropped(self):
        self.assertEqual(
            resolvable_hpo_terms(["HP:0001250", "HP:9999999", "not-an-hpo-id"]),
            ["HP:0001250"],
        )

    def test_unknown_term_does_not_zero_the_similarity(self):
        valid_only = compute_similarity_pyhpo(["HP:0001250"], ["HP:0001263"])
        self.assertGreater(valid_only, 0.0)
        with_unknown_term = compute_similarity_pyhpo(
            ["HP:0001250", "HP:9999999"], ["HP:0001263"])
        self.assertEqual(with_unknown_term, valid_only)

    def test_only_unknown_terms_scores_zero(self):
        self.assertEqual(compute_similarity_pyhpo(["HP:9999999"], ["HP:0001263"]), 0.0)


if __name__ == "__main__":
    unittest.main()
