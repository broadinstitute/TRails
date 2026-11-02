"""Unit tests for locus_annotations.py."""

import contextlib
import io
import json
import unittest

import intervaltree

import locus_annotations
import results_server
from tests.temp_file_helpers import write_temporary_file


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

    def test_n_in_a_short_motif_matches_any_base(self):
        # The catalogs write polyalanine loci as GCN or NGC; every alanine codon, in any
        # rotation and on either strand, matches, as does the other N spelling.
        for concrete_motif in ("GCG", "GCC", "GCA", "GCT", "CGC", "GGC", "ccg"):
            for catalog_motif in ("GCN", "NGC"):
                self.assertTrue(locus_annotations.motifs_match(concrete_motif, catalog_motif),
                                (concrete_motif, catalog_motif))
                self.assertTrue(locus_annotations.motifs_match(catalog_motif, concrete_motif),
                                (catalog_motif, concrete_motif))
        self.assertTrue(locus_annotations.motifs_match("NGC", "GCN"))

    def test_n_does_not_match_a_different_motif_or_length(self):
        for motif in ("AAT", "AAA", "TTA", "GCGG", "GC", "GGGGCC"):
            self.assertFalse(locus_annotations.motifs_match(motif, "GCN"), motif)


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
        self.assertEqual(record["ReferenceRegion"], "chr1:100000874-100000884")
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
        self.assertEqual(records[0]["ReferenceRegion"], "chr3:10-22")

    def test_reference_region_reads_back_as_the_locus_coordinates(self):
        # A derived ReferenceRegion pasted into the server's reference-region search must
        # select the locus's own interval: a 0-based start, like TRExplorer's
        # "chr1:10000-10108" for 1-10000-10108-TAACCC. A 1-based start would miss a 1 bp locus.
        records = [{"LocusId": "1-10000-10108-TAACCC", "Motif": "TAACCC"},
                   {"LocusId": "chrX-500-501-A", "Motif": "A"}]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["ReferenceRegion"], "chr1:10000-10108")
        for record in records:
            self.assertEqual(
                results_server.parse_reference_region(record["ReferenceRegion"]),
                ((record["Chrom"], record["Start0Based"], record["End1Based"]), None))

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

    def test_supplied_gene_region_rank_is_ignored_and_reported(self):
        # gene_region_rank is always derived from gene_region, so a rank the input matrix
        # supplied is dropped (here there is no gene_region, so the derived rank is None).
        records = [
            {"LocusId": "1-10-22-AT", "Motif": "AT", "gene_region_rank": 3},
            {"LocusId": "1-30-42-AT", "Motif": "AT"},
        ]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            locus_annotations.add_derived_locus_columns(records)
        self.assertIsNone(records[0]["gene_region_rank"])
        self.assertIsNone(records[1]["gene_region_rank"])
        self.assertIn("ignored the gene_region_rank supplied by the input matrix for 1 of 2 loci",
                      output.getvalue())

    def test_supplied_gene_region_rank_never_disagrees_with_gene_region(self):
        # The whole point of deriving the rank: the stored rank always matches the region the
        # locus actually has, even when the matrix supplied a rank for a different region.
        records = [{"LocusId": "1-10-22-AT", "Motif": "AT", "gene_region": "intron",
                    "gene_region_rank": 1}]
        with contextlib.redirect_stdout(io.StringIO()):
            locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["gene_region"], "intron")
        self.assertEqual(records[0]["gene_region_rank"], 5)

    def test_non_numeric_gene_region_rank_is_replaced_by_derived_rank(self):
        records = [{"LocusId": "1-10-22-AT", "Motif": "AT", "gene_region": "CDS",
                    "gene_region_rank": "high"}]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            locus_annotations.add_derived_locus_columns(records)
        # A string rank would make the loci-table column VARCHAR; the derived rank wins.
        self.assertEqual(records[0]["gene_region_rank"], 1)
        self.assertIn("ignored the gene_region_rank supplied by the input matrix for 1 of 1 loci",
                      output.getvalue())

    def test_supplied_coordinate_and_motif_columns_are_recomputed_and_reported(self):
        # Chrom, Start0Based, End1Based, MotifSize and CanonicalMotif are recognized matrix
        # columns but are always recomputed from the LocusId and Motif, so values that disagree
        # with them (here every one) are replaced, and each replacement is reported.
        records = [
            {"LocusId": "1-100-130-CAG", "Motif": "CAG", "Chrom": "chr2", "Start0Based": 99,
             "End1Based": 131, "MotifSize": 4, "CanonicalMotif": "CTG"},
            {"LocusId": "1-200-210-AT", "Motif": "AT", "Chrom": "chr1"},
        ]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(
            {column: records[0][column]
             for column in ("Chrom", "Start0Based", "End1Based", "MotifSize", "CanonicalMotif")},
            {"Chrom": "chr1", "Start0Based": 100, "End1Based": 130, "MotifSize": 3,
             "CanonicalMotif": "AGC"})
        self.assertEqual(output.getvalue().splitlines(), [
            "  ignored the Chrom supplied by the input matrix for 2 of 2 loci; it is always derived from LocusId",
            "  ignored the Start0Based supplied by the input matrix for 1 of 2 loci; it is always derived from LocusId",
            "  ignored the End1Based supplied by the input matrix for 1 of 2 loci; it is always derived from LocusId",
            "  ignored the MotifSize supplied by the input matrix for 1 of 2 loci; it is always derived from Motif",
            "  ignored the CanonicalMotif supplied by the input matrix for 1 of 2 loci; it is always derived from Motif",
        ])

    def test_no_supplied_gene_region_rank_prints_nothing(self):
        records = [{"LocusId": "1-10-22-AT", "Motif": "AT", "gene_region": "CDS"}]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["gene_region_rank"], 1)
        self.assertEqual(output.getvalue(), "")


class LocusIdValidationTests(unittest.TestCase):

    def test_locus_id_with_too_few_fields_names_column_value_and_file(self):
        records = [{"LocusId": "ATXN1", "Motif": "CAG"}]
        with self.assertRaises(ValueError) as raised:
            locus_annotations.add_derived_locus_columns(
                records, repeat_copy_numbers_tsv="/data/matrix.tsv")
        message = str(raised.exception)
        self.assertIn("'trid'", message)
        self.assertIn("'ATXN1'", message)
        self.assertIn("/data/matrix.tsv", message)

    def test_locus_id_with_too_many_fields_is_rejected(self):
        records = [{"LocusId": "chr1-100-110-CAG-extra", "Motif": "CAG"}]
        with self.assertRaises(ValueError) as raised:
            locus_annotations.add_derived_locus_columns(records)
        self.assertIn("'chr1-100-110-CAG-extra'", str(raised.exception))

    def test_locus_id_with_non_integer_coordinates_is_rejected(self):
        records = [{"LocusId": "chr1-start-end-CAG", "Motif": "CAG"}]
        with self.assertRaises(ValueError) as raised:
            locus_annotations.add_derived_locus_columns(
                records, repeat_copy_numbers_tsv="/data/matrix.tsv")
        message = str(raised.exception)
        self.assertIn("integer start and end", message)
        self.assertIn("'chr1-start-end-CAG'", message)
        self.assertIn("/data/matrix.tsv", message)

    def test_error_names_the_matrix_generically_when_the_path_is_unknown(self):
        with self.assertRaises(ValueError) as raised:
            locus_annotations.parse_locus_id("chr1_100_110_CAG")
        self.assertIn("the repeat-copy-numbers matrix", str(raised.exception))

    def test_well_formed_locus_id_parses(self):
        self.assertEqual(
            locus_annotations.parse_locus_id("chr1-100-110-CAG"), ("chr1", 100, 110, "CAG"))


class MotifValidationTests(unittest.TestCase):

    def test_decorated_motif_names_locus_column_value_and_file(self):
        records = [{"LocusId": "chr1-100-110-CAG", "Motif": "(CAG)n"}]
        with self.assertRaises(ValueError) as raised:
            locus_annotations.add_derived_locus_columns(
                records, repeat_copy_numbers_tsv="/data/matrix.tsv")
        message = str(raised.exception)
        self.assertIn("'motif'", message)
        self.assertIn("chr1-100-110-CAG", message)
        self.assertIn("'(CAG)n'", message)
        self.assertIn("/data/matrix.tsv", message)

    def test_iupac_ambiguity_codes_are_accepted(self):
        # compute_canonical_motif handles the whole IUPAC table, so a motif made of ambiguity
        # codes must not be rejected; only characters outside the table are.
        records = [{"LocusId": "chr1-100-110-CRG", "Motif": "crg"}]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["CanonicalMotif"], "CRG")

    def test_empty_motif_is_still_tolerated(self):
        # A blank motif cell leaves NumRepeatsInReference NULL rather than failing the build.
        records = [{"LocusId": "chr1-100-110-", "Motif": ""}]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["MotifSize"], 0)
        self.assertIsNone(records[0]["NumRepeatsInReference"])
        self.assertEqual(records[0]["CanonicalMotif"], "")


def _write_catalog_file(test_case, loci):
    """Write a catalog's loci as a temporary JSON file (see write_temporary_file) and return its path."""
    return write_temporary_file(test_case, json.dumps(loci), suffix=".json")


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
        self.catalog_path = _write_catalog_file(self, self.catalog)

    def test_only_named_disease_loci_with_regions_are_loaded(self):
        interval_trees, strchive_trees = locus_annotations.load_known_disease_loci(self.catalog_path)
        self.assertEqual(strchive_trees, {})
        self.assertIn("1", interval_trees)
        # MainReferenceRegion's start is already 0-based, so it is stored as it is.
        intervals = sorted(interval_trees["1"])
        self.assertEqual(len(intervals), 1)
        self.assertEqual(intervals[0].begin, 1000)
        self.assertEqual(intervals[0].end, 1030)
        self.assertEqual(intervals[0].data["LocusId"], "DISEASE_A")

    def test_null_pathogenic_motifs_still_match_by_repeat_unit(self):
        # Most loci of the real catalog (e.g. HTT, FMR1) have "PathogenicMotifs": null.
        catalog_path = _write_catalog_file(self, [{
            "LocusId": "DISEASE_AT", "MainReferenceRegion": "chr1:100-110", "RepeatUnit": "AT",
            "PathogenicMotifs": None, "Diseases": [{"Name": "Disease X"}],
        }])
        interval_trees, _ = locus_annotations.load_known_disease_loci(catalog_path)
        self.assertEqual(
            locus_annotations.matches_disease_locus("1-100-110-TA", interval_trees), "DISEASE_AT")
        self.assertEqual(
            locus_annotations.collect_known_disease_canonical_motifs(interval_trees), {"AT"})

    def test_cached_strchive_file_builds_trees_without_network(self):
        strchive = [{
            "id": "STR_A", "disease": "Some Disease", "chrom": "chr2",
            "start_hg38": 500, "stop_hg38": 530,
            "reference_motif_reference_orientation": ["CAG"],
        }]
        strchive_path = _write_catalog_file(self, strchive)
        _trees, strchive_trees = locus_annotations.load_known_disease_loci(
            self.catalog_path, strchive_filepath=strchive_path)
        self.assertIn("2", strchive_trees)
        intervals = sorted(strchive_trees["2"])
        self.assertEqual(intervals[0].begin, 500)
        self.assertEqual(intervals[0].end, 530)

    def test_strchive_pathogenic_motif_is_matched_and_collected(self):
        # Several STRchive loci expand as a motif that differs from the reference unit. Matching
        # only reference_motif_reference_orientation would miss exactly those expansions.
        strchive = [{
            "id": "STR_B", "disease": "Some Disease", "chrom": "chr2",
            "start_hg38": 500, "stop_hg38": 530,
            "reference_motif_reference_orientation": ["CAG"],
            "pathogenic_motif_reference_orientation": ["CCG"],
        }]
        strchive_path = _write_catalog_file(self, strchive)
        _trees, strchive_trees = locus_annotations.load_known_disease_loci(
            self.catalog_path, strchive_filepath=strchive_path)
        self.assertEqual(
            locus_annotations.matches_disease_locus("chr2-500-530-CCG", {}, strchive_trees),
            "STR_B")
        self.assertIn(
            locus_annotations.compute_canonical_motif("CCG"),
            locus_annotations.collect_known_disease_canonical_motifs({}, strchive_trees))

    def test_real_htt_locus_is_stored_at_its_locus_id_coordinates(self):
        # HTT's CAG repeat, as the GRCh38 variant catalog and STRchive list it. Both give the
        # 57 bp of 19 CAG repeats with a 0-based start, which TRExplorer and TRGT-style
        # catalogs call 4-3074876-3074933-CAG, so a locus with that id matches with Jaccard 1.
        catalog_path = _write_catalog_file(self, [{
            "LocusId": "HTT",
            "MainReferenceRegion": "chr4:3074876-3074933",
            "RepeatUnit": "CAG",
            "Diseases": [{"Symbol": "HD", "Name": "Huntington's disease"}],
        }])
        strchive_path = _write_catalog_file(self, [{
            "id": "HD_HTT", "disease": "Huntington disease", "chrom": "chr4",
            "start_hg38": 3074876, "stop_hg38": 3074933,
            "reference_motif_reference_orientation": ["CAG"],
            "pathogenic_motif_reference_orientation": ["CAG"],
        }])
        interval_trees, strchive_trees = locus_annotations.load_known_disease_loci(
            catalog_path, strchive_filepath=strchive_path)
        for trees in (interval_trees, strchive_trees):
            [interval] = trees["4"]
            self.assertEqual((interval.begin, interval.end), (3074876, 3074933))
            self.assertEqual((interval.end - interval.begin) % len("CAG"), 0)
        self.assertEqual(
            locus_annotations.matches_disease_locus("4-3074876-3074933-CAG", interval_trees), "HTT")
        self.assertEqual(
            locus_annotations.matches_disease_locus("4-3074876-3074933-CAG", {}, strchive_trees),
            "HD_HTT")


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
        self.interval_trees, _ = locus_annotations.load_known_disease_loci(
            _write_catalog_file(self, self.catalog))
        # Stored interval is [1000, 1030).

    def test_exact_overlap_and_motif_match(self):
        # Query [1000,1030) CAG: Jaccard 1.0, motif matches RepeatUnit.
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "1-1000-1030-CAG", self.interval_trees),
            "ATXN_TEST")

    def test_matches_via_reverse_complement_motif(self):
        # CTG reverse-complements to CAG canonical -> match.
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "1-1000-1030-CTG", self.interval_trees),
            "ATXN_TEST")

    def test_matches_pathogenic_motif(self):
        # CGG canonical matches the pathogenic motif CCG (rev-comp), same canonical.
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "1-1000-1030-CGG", self.interval_trees),
            "ATXN_TEST")

    def test_motif_mismatch_returns_none(self):
        # Overlap is fine but AT does not match CAG or CCG.
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "1-1000-1030-AT", self.interval_trees))

    def test_low_jaccard_overlap_returns_none(self):
        # Query [1000, 1500): overlaps but Jaccard = 30/500 < 0.66 -> no match.
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "1-1000-1500-CAG", self.interval_trees))

    def test_jaccard_of_exactly_the_threshold_matches(self):
        # [100, 133) against the disease locus [100, 150) has Jaccard 33/50 = 0.66 exactly, which
        # matches, as it does for the TRExplorer match; 32/50 = 0.64 does not.
        interval_trees = {"1": intervaltree.IntervalTree()}
        interval_trees["1"].addi(100, 150, data={"LocusId": "BOUNDARY_TEST", "RepeatUnit": "CAG"})
        self.assertEqual(locus_annotations.compute_jaccard(100, 133, 100, 150),
                         locus_annotations.MIN_JACCARD_SIMILARITY_OF_MATCHING_LOCI)
        self.assertEqual(
            locus_annotations.matches_disease_locus("1-100-133-CAG", interval_trees), "BOUNDARY_TEST")
        self.assertIsNone(locus_annotations.matches_disease_locus("1-100-132-CAG", interval_trees))

    def test_no_overlap_returns_none(self):
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "1-5000-5030-CAG", self.interval_trees))

    def test_wrong_chromosome_returns_none(self):
        self.assertIsNone(
            locus_annotations.matches_disease_locus(
                "2-1000-1030-CAG", self.interval_trees))

    def test_chr_prefix_in_query_is_normalized(self):
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "chr1-1000-1030-CAG", self.interval_trees),
            "ATXN_TEST")

    def test_strchive_fallback(self):
        # No catalog match (wrong chrom), but STRchive has a match.
        strchive_trees = {
            "9": intervaltree.IntervalTree(),
        }
        strchive_trees["9"].addi(999, 1030, data={
            "id": "STR_FALLBACK",
            "reference_motif_reference_orientation": ["CAG"],
        })
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "9-999-1030-CAG", self.interval_trees, strchive_trees),
            "STR_FALLBACK")

    def test_strchive_fallback_matches_a_pathogenic_motif(self):
        # The pathogenic motif differs from the reference one, which is the case this fallback
        # exists for: matching only the reference motif would miss the expansion.
        strchive_trees = {"9": intervaltree.IntervalTree()}
        strchive_trees["9"].addi(999, 1030, data={
            "id": "STR_PATHOGENIC",
            "reference_motif_reference_orientation": ["AAAAG"],
            "pathogenic_motif_reference_orientation": ["AAGGG"],
        })
        self.assertEqual(
            locus_annotations.matches_disease_locus(
                "9-999-1030-AAGGG", self.interval_trees, strchive_trees),
            "STR_PATHOGENIC")
        # The same locus's motif also counts as known, so IsKnownMotif agrees with
        # KnownDiseaseLocus.
        self.assertIn(
            locus_annotations.compute_canonical_motif("AAGGG"),
            locus_annotations.collect_known_disease_canonical_motifs({}, strchive_trees))


class StrchiveLocusMotifsTests(unittest.TestCase):
    """The one definition of a STRchive locus's motif list, shared by the build and the server."""

    def test_both_lists_are_returned_reference_first(self):
        self.assertEqual(
            locus_annotations.strchive_locus_motifs({
                "reference_motif_reference_orientation": ["AAAAG"],
                "pathogenic_motif_reference_orientation": ["AAGGG", "AGGGG"],
            }),
            ["AAAAG", "AAGGG", "AGGGG"])

    def test_missing_and_null_lists_are_treated_as_empty(self):
        self.assertEqual(locus_annotations.strchive_locus_motifs({}), [])
        self.assertEqual(
            locus_annotations.strchive_locus_motifs({
                "reference_motif_reference_orientation": None,
                "pathogenic_motif_reference_orientation": ["AAGGG"],
            }),
            ["AAGGG"])


class KnownMotifTests(unittest.TestCase):

    def setUp(self):
        self.catalog = [{
            "LocusId": "L1",
            "MainReferenceRegion": "chr1:1000-1030",
            "RepeatUnit": "CAG",
            "PathogenicMotifs": ["CCCCGG", "NGC"],
            "Diseases": [{"Name": "D"}],
        }]
        self.interval_trees, _ = locus_annotations.load_known_disease_loci(
            _write_catalog_file(self, self.catalog))

    def test_collect_canonical_motifs_expands_n_into_concrete_motifs(self):
        known = locus_annotations.collect_known_disease_canonical_motifs(
            self.interval_trees)
        from motif_utilities import compute_canonical_motif
        self.assertIn(compute_canonical_motif("CAG"), known)
        self.assertIn(compute_canonical_motif("CCCCGG"), known)
        # NGC stands for AGC, CGC, GGC and TGC; the N spelling itself is not a locus motif.
        for concrete_motif in ("AGC", "CGC", "GGC", "TGC"):
            self.assertIn(compute_canonical_motif(concrete_motif), known, concrete_motif)
        self.assertNotIn(compute_canonical_motif("NGC"), known)

    def test_long_motif_with_n_is_not_expanded(self):
        interval_trees, _ = locus_annotations.load_known_disease_loci(_write_catalog_file(self, [{
            "LocusId": "LONG", "MainReferenceRegion": "chr1:1000-1120", "RepeatUnit": "GGCTNNGGGNGC",
            "Diseases": [{"Name": "D"}]}]))
        self.assertEqual(locus_annotations.collect_known_disease_canonical_motifs(interval_trees), set())

    def test_compute_is_known_motif(self):
        known = locus_annotations.collect_known_disease_canonical_motifs(
            self.interval_trees)
        from motif_utilities import compute_canonical_motif
        self.assertEqual(
            locus_annotations.compute_is_known_motif(compute_canonical_motif("CAG"), known), 1)
        self.assertEqual(
            locus_annotations.compute_is_known_motif("ZZZ", known), 0)
        # A locus motif written with N is known when a concrete motif it stands for is.
        self.assertEqual(locus_annotations.compute_is_known_motif(compute_canonical_motif("GCN"), known), 1)
        self.assertEqual(locus_annotations.compute_is_known_motif(compute_canonical_motif("ANT"), known), 0)


class PolyalanineDiseaseLocusTests(unittest.TestCase):
    """Real ARX and PHOX2B entries, whose catalogs write the polyalanine repeat as NGC / GCN."""

    def setUp(self):
        # As in variant_catalog_without_offtargets.GRCh38.json and STRchive-loci.json.
        catalog_path = _write_catalog_file(self, [
            {"LocusId": "ARX_1", "MainReferenceRegion": "chrX:25013649-25013697", "RepeatUnit": "NGC",
             "PathogenicMotifs": None,
             "Diseases": [{"Name": "Developmental and epileptic encephalopathy-1"}]},
            {"LocusId": "PHOX2B", "MainReferenceRegion": "chr4:41745972-41746032", "RepeatUnit": "GCN",
             "PathogenicMotifs": None, "Diseases": [{"Name": "Congenital central hypoventilation"}]},
        ])
        strchive_path = _write_catalog_file(self, [{
            "id": "EIEE1_ARX", "disease": "Developmental and epileptic encephalopathy-1",
            "chrom": "chrX", "start_hg38": 25013649, "stop_hg38": 25013697,
            "reference_motif_reference_orientation": ["GCN"],
            "pathogenic_motif_reference_orientation": ["NGC"],
        }])
        self.interval_trees, self.strchive_trees = locus_annotations.load_known_disease_loci(
            catalog_path, strchive_filepath=strchive_path)

    def test_a_concrete_motif_matches_the_n_motif(self):
        # TRGT's pathogenic-repeat catalog calls ARX's repeat GCG; PHOX2B is often genotyped as GCC.
        self.assertEqual(locus_annotations.matches_disease_locus(
            "chrX-25013649-25013697-GCG", self.interval_trees, self.strchive_trees), "ARX_1")
        self.assertEqual(locus_annotations.matches_disease_locus(
            "chr4-41745972-41746032-GCC", self.interval_trees, self.strchive_trees), "PHOX2B")
        self.assertEqual(locus_annotations.matches_disease_locus(
            "chrX-25013649-25013697-GCG", {}, self.strchive_trees), "EIEE1_ARX")

    def test_a_different_motif_at_the_same_locus_does_not_match(self):
        for locus_id in ("chrX-25013649-25013697-AAT", "chr4-41745972-41746032-GGAA"):
            self.assertIsNone(locus_annotations.matches_disease_locus(
                locus_id, self.interval_trees, self.strchive_trees), locus_id)

    def test_is_known_motif_agrees_with_known_disease_locus(self):
        known = locus_annotations.collect_known_disease_canonical_motifs(
            self.interval_trees, self.strchive_trees)
        for motif, expected in (("GCG", 1), ("GCC", 1), ("AAT", 0)):
            self.assertEqual(locus_annotations.compute_is_known_motif(
                locus_annotations.compute_canonical_motif(motif), known), expected, motif)


class IsInMendelianGeneTests(unittest.TestCase):

    def test_is_in_mendelian_gene(self):
        # A table with neither an inheritance nor a sources column counts every gene it lists.
        gene_lookup = {"ENSG1": {}, "ENSG2": {}}
        self.assertEqual(locus_annotations.is_in_mendelian_gene("ENSG1", gene_lookup), 1)
        self.assertEqual(locus_annotations.is_in_mendelian_gene("ENSG_X", gene_lookup), 0)
        self.assertEqual(locus_annotations.is_in_mendelian_gene(None, gene_lookup), 0)

    def test_is_in_mendelian_gene_needs_mendelian_evidence_when_the_table_has_the_columns(self):
        # Rows shaped the way read_gene_table returns genelookup's table: a blank cell is NaN.
        nan = float("nan")
        gene_lookup = {
            "OMIM_GENE": {"inheritance": nan, "sources": "2: OMIM, pLI_v4:1.0"},
            "PANELAPP_GENE": {"inheritance": nan, "sources": "3: GWAS, PanelAppAU, dbNSFP"},
            "INHERITANCE_ONLY": {"inheritance": "AR", "sources": "2: Intolerome, dbNSFP"},
            "CONSTRAINED_ONLY": {"inheritance": nan, "sources": "1: pLI_v4:1.0"},
            "GWAS_AND_CLINVAR": {"inheritance": nan, "sources": "2: GWAS, ClinVar"},
            "BLANK_CELLS": {"inheritance": nan, "sources": nan},
            "NA_INHERITANCE": {"inheritance": "NA", "sources": "1: dbNSFP"},
        }
        self.assertEqual(
            {gene_id: locus_annotations.is_in_mendelian_gene(gene_id, gene_lookup) for gene_id in gene_lookup},
            {"OMIM_GENE": 1, "PANELAPP_GENE": 1, "INHERITANCE_ONLY": 1, "CONSTRAINED_ONLY": 0,
             "GWAS_AND_CLINVAR": 0, "BLANK_CELLS": 0, "NA_INHERITANCE": 0})


class StdevPercentileTests(unittest.TestCase):
    """add_population_stdev_percentile_columns: a locus's stdev rank over the number of ranked loci."""

    def test_rank_over_total(self):
        records = [{"HPRC256_StdevRankByMotif": 5, "HPRC256_StdevRankTotalNumberByMotif": 20},
                   {"AoU1027_StdevRankByMotif": 3, "AoU1027_StdevRankTotalNumberByMotif": 0,
                    "AoU1027_StdevPercentile": 0.4},
                   {"HPRC256_StdevRankByMotif": None, "HPRC256_StdevRankTotalNumberByMotif": 20}]
        locus_annotations.add_population_stdev_percentile_columns(records)
        self.assertEqual(records[0]["HPRC256_StdevPercentile"], 0.25)
        self.assertEqual(records[1]["AoU1027_StdevPercentile"], 0.4)  # no total: what it had stays
        self.assertNotIn("HPRC256_StdevPercentile", records[2])


if __name__ == "__main__":
    unittest.main()
