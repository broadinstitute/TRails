"""Unit tests for trexplorer_annotations: slimming TRExplorer's catalog table and matching loci to it."""

import gzip
import os
import shutil
import tempfile
import unittest

import duckdb_compat
import locus_annotations
import trexplorer_annotations

# Every column build_trexplorer_annotations_database reads, plus one it drops
CATALOG_COLUMNS = (
    ["LocusId", "ReferenceRegion", "LocusStructure", "Motif", "CanonicalMotif", "Source", "GencodeGeneId",
     "GencodeGeneRegion", "VariationCluster", "VariationClusterSizeDiff", "VariationClusterFilterReason",
     "ReferenceRepeatPurity"]
    + trexplorer_annotations.INTEGER_POPULATION_COLUMNS + trexplorer_annotations.FLOAT_POPULATION_COLUMNS)


def catalog_row(locus_id, **values):
    """Returns one catalog row for a locus id chrom-start-end-motif, with defaults for the rest."""
    chrom, start, end, motif = locus_id.split("-")
    row = {"LocusId": locus_id, "ReferenceRegion": f"chr{chrom}:{start}-{end}", "Motif": motif,
           "CanonicalMotif": locus_annotations.compute_canonical_motif(motif, include_reverse_complement=True),
           "Source": "TRExplorerV1:PerfectRepeatsInReference", "GencodeGeneRegion": "intergenic",
           "ReferenceRepeatPurity": "1.0"}
    row.update(values)
    return row


def user_locus(locus_id, **values):
    """Returns a locus record as add_derived_locus_columns leaves it."""
    chrom, start, end, motif = locus_id.split("-")
    record = {"LocusId": locus_id, "Chrom": f"chr{chrom.replace('chr', '')}", "Start0Based": int(start),
              "End1Based": int(end), "Motif": motif,
              "CanonicalMotif": locus_annotations.compute_canonical_motif(motif, include_reverse_complement=True)}
    record.update(values)
    return record


class TRExplorerAnnotationsTestCase(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def build_database(self, rows):
        """Writes rows as a gzipped catalog table and builds the annotations database from it."""
        table_path = os.path.join(self.tmpdir, "catalog.tsv.gz")
        with gzip.open(table_path, "wt") as f:
            f.write("\t".join(CATALOG_COLUMNS) + "\n")
            for row in rows:
                f.write("\t".join(str(row.get(column, "")) for column in CATALOG_COLUMNS) + "\n")
        database_path = os.path.join(self.tmpdir, "annotations.duckdb")
        self.assertEqual(trexplorer_annotations.build_trexplorer_annotations_database(table_path, database_path),
                         len(rows))
        return database_path

    def matched_locus_ids(self, rows, records):
        """Returns the TRExplorer locus id each record matched, or None."""
        matches = trexplorer_annotations._match_loci(records, self.build_database(rows))
        return [matches[i]["LocusId"] if i in matches else None for i in range(len(records))]


class BuildDatabaseTests(TRExplorerAnnotationsTestCase):

    def test_keeps_typed_columns_and_parses_the_region(self):
        database_path = self.build_database([catalog_row(
            "1-100-130-CAG", GencodeGeneId="ENSG00000197386", HPRC256_99thPercentile="12.5",
            HPRC256_MaxAllele="40", VariationClusterSizeDiff="4.0")])
        connection = duckdb_compat.connect(database_path, read_only=True)
        try:
            row = dict(zip(["Chrom", "Start0Based", "End1Based", "HPRC256_99thPercentile", "HPRC256_MaxAllele",
                            "VariationClusterSizeDiff", "GencodeGeneId", "AoU1027_Stdev", "VariationCluster"],
                           connection.execute(
                               "SELECT Chrom, Start0Based, End1Based, HPRC256_99thPercentile, HPRC256_MaxAllele, "
                               "VariationClusterSizeDiff, GencodeGeneId, AoU1027_Stdev, VariationCluster "
                               "FROM trexplorer_loci").fetchone()))
            columns = duckdb_compat.table_columns(connection, "trexplorer_loci")
        finally:
            connection.close()
        self.assertEqual(row, {"Chrom": "1", "Start0Based": 100, "End1Based": 130, "HPRC256_99thPercentile": 12.5,
                               "HPRC256_MaxAllele": 40, "VariationClusterSizeDiff": 4,
                               "GencodeGeneId": "ENSG00000197386", "AoU1027_Stdev": None, "VariationCluster": None})
        self.assertNotIn("LocusStructure", columns)
        self.assertFalse(os.path.exists(database_path + ".tmp"))


class MatchLociTests(TRExplorerAnnotationsTestCase):

    def test_the_locus_with_the_same_id_wins_over_a_more_variable_neighbor(self):
        rows = [catalog_row("1-100-130-CAG", HPRC256_Stdev="1.0"),
                catalog_row("1-101-130-CAG", HPRC256_Stdev="9.0")]
        self.assertEqual(self.matched_locus_ids(rows, [user_locus("chr1-100-130-CAG")]), ["1-100-130-CAG"])

    def test_otherwise_the_most_variable_overlapping_locus_wins(self):
        rows = [catalog_row("1-100-130-CAG", HPRC256_Stdev="1.0"),
                catalog_row("1-101-130-CAG", HPRC256_Stdev="9.0"),
                catalog_row("1-99-130-CAG", AoU1027_Stdev="50.0")]
        self.assertEqual(self.matched_locus_ids(rows, [user_locus("chr1-100-131-CAG")]), ["1-101-130-CAG"])

    def test_aou_stdev_then_purity_then_locus_id_break_the_rest_of_the_ties(self):
        rows = [catalog_row("1-100-130-CAG", ReferenceRepeatPurity="1.0"),
                catalog_row("1-101-130-CAG", AoU1027_Stdev="2.0")]
        self.assertEqual(self.matched_locus_ids(rows, [user_locus("1-100-131-CAG")]), ["1-101-130-CAG"])
        rows = [catalog_row("1-100-130-CAG", ReferenceRepeatPurity="0.9"),
                catalog_row("1-101-130-CAG", ReferenceRepeatPurity="0.8")]
        self.assertEqual(self.matched_locus_ids(rows, [user_locus("1-100-131-CAG")]), ["1-100-130-CAG"])
        rows = [catalog_row("1-100-130-CAG"), catalog_row("1-101-130-CAG")]
        self.assertEqual(self.matched_locus_ids(rows, [user_locus("1-100-131-CAG")]), ["1-101-130-CAG"])

    def test_too_little_overlap_matches_nothing(self):
        # Jaccard 20 / 40 = 0.5
        rows = [catalog_row("1-100-130-CAG")]
        self.assertEqual(self.matched_locus_ids(rows, [user_locus("1-110-140-CAG")]), [None])

    def test_the_motif_rule(self):
        rows = [catalog_row("1-100-130-CAG"), catalog_row("2-100-170-AAAAAAC"), catalog_row("3-100-160-AAAAAC")]
        records = [user_locus("1-100-130-AGC"),        # a rotation of CAG: same canonical motif
                   user_locus("1-100-130-CAA"),        # a different 3 bp motif
                   user_locus("2-100-170-GGGGGGT"),    # another 7 bp motif: same length is enough
                   user_locus("3-100-160-AAAAAAC")]    # 7 bp against 6 bp: never
        self.assertEqual(self.matched_locus_ids(rows, records), ["1-100-130-CAG", None, "2-100-170-AAAAAAC", None])

    def test_the_same_coordinates_on_another_chromosome_do_not_match(self):
        rows = [catalog_row("2-100-130-CAG"), catalog_row("X-100-130-CAG")]
        self.assertEqual(self.matched_locus_ids(rows, [user_locus("1-100-130-CAG"), user_locus("chrX-100-130-CAG")]),
                         [None, "X-100-130-CAG"])


class AddTRExplorerAnnotationsTests(TRExplorerAnnotationsTestCase):

    def test_fills_empty_columns_and_keeps_supplied_ones(self):
        database_path = self.build_database([catalog_row(
            "1-100-130-CAG", GencodeGeneId="ENSG00000197386", GencodeGeneRegion="CDS", HPRC256_99thPercentile="12.5",
            AoU1027_99thPercentile="14.0", VariationCluster="1-90-140", VariationClusterSizeDiff="20")])
        records = [user_locus("chr1-100-130-CAG", AoU1027_99thPercentile=99, gene_region=None, gene_region_rank=None),
                   user_locus("chr1-5000-5030-CAG")]
        self.assertEqual(trexplorer_annotations.add_trexplorer_annotations(records, database_path), (1, 0))
        annotated = records[0]
        self.assertEqual(annotated["HPRC256_99thPercentile"], 12.5)
        self.assertEqual(annotated["AoU1027_99thPercentile"], 99)  # supplied by the matrix, so kept
        self.assertEqual((annotated["gene_id"], annotated["gene_region"]), ("ENSG00000197386", "CDS"))
        self.assertEqual(annotated["gene_region_rank"], locus_annotations.gene_region_rank("CDS"))
        self.assertEqual((annotated["TRExplorerLocusId"], annotated["TRExplorerMotif"]), ("1-100-130-CAG", "CAG"))
        self.assertEqual(annotated["TRExplorerReferenceRegion"], "chr1:100-130")
        self.assertEqual(annotated["TRExplorerLocusJaccardSimilarity"], 1.0)
        self.assertEqual((annotated["VariationCluster"], annotated["VariationClusterSizeDiff"]), ("1-90-140", 20))
        self.assertNotIn("HPRC256_99thPercentile", records[1])

    def test_counts_overlap_matches_separately(self):
        database_path = self.build_database([catalog_row("1-100-130-CAG")])
        records = [user_locus("1-100-131-CAG")]
        self.assertEqual(trexplorer_annotations.add_trexplorer_annotations(records, database_path), (0, 1))
        self.assertEqual(records[0]["TRExplorerLocusJaccardSimilarity"], round(30 / 31, 4))


class StdevPercentileTests(unittest.TestCase):

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
