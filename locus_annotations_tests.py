"""Unit tests for locus_annotations.py."""

import json
import os
import tempfile
import unittest

import intervaltree

import locus_annotations


class ComputeJaccardTests(unittest.TestCase):

    def test_disjoint_intervals_return_zero(self):
        self.assertEqual(locus_annotations.compute_jaccard(0, 10, 20, 30), 0.0)

    def test_touching_intervals_return_zero(self):
        # End-exclusive: [0,10) and [10,20) do not overlap.
        self.assertEqual(locus_annotations.compute_jaccard(0, 10, 10, 20), 0.0)

    def test_identical_intervals_return_one(self):
        self.assertEqual(locus_annotations.compute_jaccard(100, 200, 100, 200), 1.0)

    def test_partial_overlap(self):
        # [0,10) vs [5,15): overlap=5, union=15 -> 1/3.
        self.assertAlmostEqual(
            locus_annotations.compute_jaccard(0, 10, 5, 15), 5.0 / 15.0)

    def test_contained_interval(self):
        # [0,100) vs [25,75): overlap=50, union=100 -> 0.5.
        self.assertAlmostEqual(
            locus_annotations.compute_jaccard(0, 100, 25, 75), 0.5)


class MotifsMatchTests(unittest.TestCase):

    def test_short_motifs_match_by_canonical_form(self):
        # CAG and AGC are rotations -> same canonical -> match.
        self.assertTrue(locus_annotations.motifs_match("CAG", "AGC"))

    def test_short_motifs_match_by_reverse_complement(self):
        # GAA canonical AAG; CTT reverse-complements to AAG -> match.
        self.assertTrue(locus_annotations.motifs_match("GAA", "CTT"))

    def test_short_motifs_do_not_match(self):
        self.assertFalse(locus_annotations.motifs_match("CAG", "GCC"))

    def test_long_motifs_match_by_length(self):
        # >6 bp: matched purely on equal length, not sequence.
        self.assertTrue(locus_annotations.motifs_match("ACGTACG", "TTTTTTT"))

    def test_long_motifs_do_not_match_when_lengths_differ(self):
        self.assertFalse(locus_annotations.motifs_match("ACGTACG", "ACGTACGT"))

    def test_empty_motif_returns_false(self):
        self.assertFalse(locus_annotations.motifs_match("", "CAG"))
        self.assertFalse(locus_annotations.motifs_match("CAG", None))


class GeneRegionRankTests(unittest.TestCase):

    def test_known_regions(self):
        self.assertEqual(locus_annotations.gene_region_rank("CDS"), 1)
        self.assertEqual(locus_annotations.gene_region_rank("intergenic"), 7)
        self.assertEqual(locus_annotations.gene_region_rank("5' UTR"), 2)

    def test_unknown_region_returns_none(self):
        self.assertIsNone(locus_annotations.gene_region_rank("nonsense"))

    def test_missing_region_returns_none(self):
        self.assertIsNone(locus_annotations.gene_region_rank(None))


class AddDerivedLocusColumnsTests(unittest.TestCase):

    def test_parses_locus_id_and_computes_columns(self):
        records = [{"LocusId": "1-100000874-100000884-T", "Motif": "T"}]
        locus_annotations.add_derived_locus_columns(records, source_label="cohortA")
        record = records[0]
        self.assertEqual(record["Chrom"], "chr1")
        self.assertEqual(record["Start0Based"], 100000874)
        self.assertEqual(record["End1Based"], 100000884)
        self.assertEqual(record["ReferenceRegion"], "1:100000875-100000884")
        self.assertEqual(record["MotifSize"], 1)
        # (100000884 - 100000874) // 1 == 10
        self.assertEqual(record["NumRepeatsInReference"], 10)
        self.assertEqual(record["Source"], "cohortA")
        self.assertEqual(record["CanonicalMotif"], "A")
        self.assertIsNone(record["gene_id"])
        self.assertIsNone(record["gene_region"])
        self.assertIsNone(record["gene_region_rank"])

    def test_default_source_label_is_empty_string(self):
        records = [{"LocusId": "chr2-50-62-CAG", "Motif": "CAG"}]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["Source"], "")
        self.assertEqual(records[0]["Chrom"], "chr2")
        # (62 - 50) // 3 == 4
        self.assertEqual(records[0]["NumRepeatsInReference"], 4)
        self.assertEqual(records[0]["CanonicalMotif"], "AGC")

    def test_chr_prefix_in_locus_id_is_normalized(self):
        records = [{"LocusId": "chr3-10-22-AT", "Motif": "AT"}]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["Chrom"], "chr3")
        self.assertEqual(records[0]["ReferenceRegion"], "3:11-22")

    def test_preserves_existing_reference_region_and_num_repeats(self):
        records = [{
            "LocusId": "1-100-130-CAG", "Motif": "CAG",
            "ReferenceRegion": "1:custom", "NumRepeatsInReference": 999,
        }]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["ReferenceRegion"], "1:custom")
        self.assertEqual(records[0]["NumRepeatsInReference"], 999)

    def test_gene_columns_pass_through_and_aliases(self):
        records = [
            {"LocusId": "1-10-22-AT", "Motif": "AT", "gene_id": "ENSG1",
             "gene_region": "CDS"},
            {"LocusId": "1-30-42-AT", "Motif": "AT", "GencodeGeneId": "ENSG2",
             "GencodeGeneRegion": "intron"},
        ]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["gene_id"], "ENSG1")
        self.assertEqual(records[0]["gene_region"], "CDS")
        self.assertEqual(records[0]["gene_region_rank"], 1)
        self.assertEqual(records[1]["gene_id"], "ENSG2")
        self.assertEqual(records[1]["gene_region"], "intron")
        self.assertEqual(records[1]["gene_region_rank"], 5)


def _make_catalog_file(loci):
    """Write a catalog JSON to a temp file and return its path."""
    handle = tempfile.NamedTemporaryFile(
        mode="wt", suffix=".json", delete=False)
    json.dump(loci, handle)
    handle.close()
    return handle.name


class LoadKnownDiseaseLociTests(unittest.TestCase):

    def setUp(self):
        self.catalog = [
            {
                "LocusId": "DISEASE_A",
                "MainReferenceRegion": "chr1:1000-1030",
                "RepeatUnit": "CAG",
                "Diseases": [{"Name": "Disease A"}],
            },
            {
                # No named disease -> skipped.
                "LocusId": "NO_DISEASE",
                "MainReferenceRegion": "chr1:2000-2030",
                "RepeatUnit": "CAG",
                "Diseases": [{"Name": ""}],
            },
            {
                # Missing MainReferenceRegion -> skipped.
                "LocusId": "NO_REGION",
                "RepeatUnit": "CAG",
                "Diseases": [{"Name": "Disease C"}],
            },
        ]
        self.catalog_path = _make_catalog_file(self.catalog)

    def tearDown(self):
        os.remove(self.catalog_path)

    def test_only_named_disease_loci_with_regions_are_loaded(self):
        interval_trees, strchive_trees, locus_lookup = locus_annotations.load_known_disease_loci(
            self.catalog_path)
        self.assertEqual(strchive_trees, {})
        self.assertEqual(locus_lookup, {})  # not built unless build_locus_lookup=True
        self.assertIn("1", interval_trees)
        # Start is converted from 1-based 1000 to 0-based 999.
        intervals = sorted(interval_trees["1"])
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0].begin, 999)
        self.assertEqual(intervals[0].end, 1030)
        self.assertEqual(intervals[0].data["LocusId"], "DISEASE_A")

    def test_build_locus_lookup_keys_by_locus_id_and_coordinate(self):
        _trees, _strchive, locus_lookup = locus_annotations.load_known_disease_loci(
            self.catalog_path, build_locus_lookup=True)
        # Keyed by LocusId and by "chrom-start0-end1-RepeatUnit".
        self.assertIn("DISEASE_A", locus_lookup)
        self.assertIn("1-999-1030-CAG", locus_lookup)
        self.assertEqual(locus_lookup["DISEASE_A"]["LocusId"], "DISEASE_A")

    def test_no_network_fetch_by_default(self):
        # fetch_strchive defaults False; strchive_trees must be empty and no
        # network access should occur (a network call here would raise).
        _trees, strchive_trees, _lookup = locus_annotations.load_known_disease_loci(
            self.catalog_path)
        self.assertEqual(strchive_trees, {})

    def test_cached_strchive_file_builds_trees_without_network(self):
        strchive = [{
            "locus_id": "STR_A", "disease": "Some Disease", "chrom": "chr2",
            "start_hg38": 500, "stop_hg38": 530,
            "reference_motif_reference_orientation": ["CAG"],
        }]
        strchive_path = _make_catalog_file(strchive)
        try:
            _trees, strchive_trees, _lookup = locus_annotations.load_known_disease_loci(
                self.catalog_path, strchive_filepath=strchive_path)
            self.assertIn("2", strchive_trees)
            intervals = sorted(strchive_trees["2"])
            self.assertEqual(intervals[0].begin, 499)
            self.assertEqual(intervals[0].end, 530)
        finally:
            os.remove(strchive_path)


class MatchesDiseaseLocusTests(unittest.TestCase):

    def setUp(self):
        # Build a tiny catalog inline: one disease locus at chr1:1000-1030 (CAG).
        self.catalog = [{
            "LocusId": "ATXN_TEST",
            "MainReferenceRegion": "chr1:1000-1030",
            "RepeatUnit": "CAG",
            "PathogenicMotifs": ["CCG"],
            "Diseases": [{"Name": "Test ataxia"}],
        }]
        self.catalog_path = _make_catalog_file(self.catalog)
        self.interval_trees, _, _ = locus_annotations.load_known_disease_loci(
            self.catalog_path)
        # Stored interval is [999, 1030).

    def tearDown(self):
        os.remove(self.catalog_path)

    def test_exact_overlap_and_motif_match(self):
        # Query [999,1030) CAG: Jaccard 1.0, motif matches RepeatUnit.
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "1-999-1030-CAG", self.interval_trees),
            "ATXN_TEST")

    def test_matches_via_reverse_complement_motif(self):
        # CTG reverse-complements to CAG canonical -> match.
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "1-999-1030-CTG", self.interval_trees),
            "ATXN_TEST")

    def test_matches_pathogenic_motif(self):
        # CGG canonical matches the pathogenic motif CCG (rev-comp), same canonical.
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "1-999-1030-CGG", self.interval_trees),
            "ATXN_TEST")

    def test_motif_mismatch_returns_none(self):
        # Overlap is fine but AT does not match CAG or CCG.
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "1-999-1030-AT", self.interval_trees))

    def test_low_jaccard_overlap_returns_none(self):
        # Query [999, 1500): overlaps but Jaccard = 31/501 < 0.66 -> no match.
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "1-999-1500-CAG", self.interval_trees))

    def test_no_overlap_returns_none(self):
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "1-5000-5030-CAG", self.interval_trees))

    def test_wrong_chromosome_returns_none(self):
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "2-999-1030-CAG", self.interval_trees))

    def test_chr_prefix_in_query_is_normalized(self):
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "chr1-999-1030-CAG", self.interval_trees),
            "ATXN_TEST")

    def test_strchive_fallback(self):
        # No catalog match (wrong chrom), but STRchive has a match.
        strchive_trees = {
            "9": intervaltree.IntervalTree(),
        }
        strchive_trees["9"].addi(999, 1030, data={
            "locus_id": "STR_FALLBACK",
            "reference_motif_reference_orientation": ["CAG"],
        })
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "9-999-1030-CAG", self.interval_trees, strchive_trees),
            "STR_FALLBACK")


class KnownMotifAndMendelianGeneTests(unittest.TestCase):

    def setUp(self):
        self.catalog = [{
            "LocusId": "L1",
            "MainReferenceRegion": "chr1:1000-1030",
            "RepeatUnit": "CAG",
            "PathogenicMotifs": ["CCCCGG", "NGC"],
            "Diseases": [{"Name": "D"}],
        }]
        self.catalog_path = _make_catalog_file(self.catalog)
        self.interval_trees, _, _ = locus_annotations.load_known_disease_loci(
            self.catalog_path)

    def tearDown(self):
        os.remove(self.catalog_path)

    def test_collect_canonical_motifs_skips_N_containing(self):
        known = locus_annotations.collect_known_disease_canonical_motifs(
            self.interval_trees)
        from motif_utilities import compute_canonical_motif
        self.assertIn(compute_canonical_motif("CAG"), known)
        self.assertIn(compute_canonical_motif("CCCCGG"), known)
        # NGC contains N and is excluded.
        self.assertNotIn(compute_canonical_motif("NGC", include_reverse_complement=False), known)

    def test_compute_is_known_motif(self):
        known = locus_annotations.collect_known_disease_canonical_motifs(
            self.interval_trees)
        from motif_utilities import compute_canonical_motif
        self.assertEqual(
            locus_annotations.compute_is_known_motif(compute_canonical_motif("CAG"), known), 1)
        self.assertEqual(
            locus_annotations.compute_is_known_motif("ZZZ", known), 0)

    def test_is_in_mendelian_gene(self):
        gene_lookup = {"ENSG1": {}, "ENSG2": {}}
        self.assertEqual(locus_annotations.is_in_mendelian_gene("ENSG1", gene_lookup), 1)
        self.assertEqual(locus_annotations.is_in_mendelian_gene("ENSG_X", gene_lookup), 0)
        self.assertEqual(locus_annotations.is_in_mendelian_gene(None, gene_lookup), 0)


if __name__ == "__main__":
    unittest.main()
