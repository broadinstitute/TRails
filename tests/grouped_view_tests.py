"""Unit tests for the "Group results by" view: the /api/v1/outlier_samples endpoint, its export,
and the /api/v1/loci restriction to one sample or family that a row of the view expands to.

A small throwaway database stands in for a real one. Its affected_unsolved_outliers_above_unaffected
table is written by hand here, with the rows the build derives from the loci table's outlier lists
(see LOCI), so these tests cover only how the server reads it.
"""

import gzip
import os
import unittest

import duckdb_compat
import results_server
from result_database import AFFECTED_UNSOLVED_OUTLIERS_TABLE, write_population_gate_cohorts
from tests.server_test_fixtures import ServerTestCase

# sample_id -> (family_id, affected status, analysis status, phenotype). S4 is solved and S5
# unaffected, so neither is ever counted; S6's analysis status is blank, which still counts as
# affected and unsolved, and S6 has no phenotype.
SAMPLES = {
    "S1": ("F1", "affected", "unsolved", "seizures, ataxia"),
    "S2": ("F1", "affected", "unsolved", "ataxia"),
    "S3": ("F2", "affected", "unsolved", "autism"),
    "S4": ("F3", "affected", "solved", "ataxia"),
    "S5": ("F4", "unaffected", "unaffected", ""),
    "S6": ("F5", "affected", "", ""),
}

# LocusId, OutlierSampleIds_AllAlleles, IsInMendelianGene, CanonicalMotif, gene_region,
# HPRC256_99thPercentile, and the stored #Aff and #Fam counts the outlier list implies:
#   L1: S1 and S2 are above S5 (unaffected) and the 99th percentile; S3 is below S5.
#   L2: S3 only; S4 is solved, which stops the walk before S1.
#   L3: S1 only; S6 ties the 99th percentile, which is not above it.
#   L4: S2 ties S5 (unaffected), so nothing counts.
#   L5: no outliers.
#   L6: no population data, so both samples count.
LOCI = [
    ("L1", "50x:S1,40x:S2,30x:S5,20x:S3", 1, "AGC", "CDS", 25.0, 2, 1),
    ("L2", "60x:S3,55x:S4,50x:S1", 0, "AT", "intron", 10.0, 1, 1),
    ("L3", "45x:S1,44x:S6", 1, "AGC", "CDS", 44.0, 1, 1),
    ("L4", "30x:S2,30x:S5", 0, "AAG", "intergenic", 10.0, 0, 0),
    ("L5", None, 0, "AAG", "intergenic", 10.0, None, None),
    ("L6", "70x:S6,20x:S1", 1, "CCG", "CDS", None, 2, 2),
]

# The rows the build writes for LOCI: LocusId, sample_id, family_id, allele_size
AFFECTED_UNSOLVED_OUTLIERS = [
    ("L1", "S1", "F1", 50), ("L1", "S2", "F1", 40), ("L2", "S3", "F2", 60),
    ("L3", "S1", "F1", 45), ("L6", "S1", "F1", 20), ("L6", "S6", "F5", 70),
]

# locus_id, sample_id, gene_phenotype_similarity
PHENOTYPE_SCORES = [("L1", "S1", 0.8), ("L3", "S1", 0.3), ("L2", "S3", 0.5)]


def build_test_database(directory, with_outliers_table=True, population_gate_cohorts=None):
    """Write the loci, phenotype score and sample tables, and optionally the grouped view's table.

    Args:
        directory: Directory to write the database and the sample table into.
        with_outliers_table: Whether to write AFFECTED_UNSOLVED_OUTLIERS_TABLE.
        population_gate_cohorts: Cohorts to record as the build's population gate, or None to
            record none, which reads back as the default long-read cohorts.

    Returns:
        Tuple (database path, sample table path).
    """
    db_path = os.path.join(directory, "grouped.duckdb")
    conn = duckdb_compat.connect(db_path)
    conn.execute("""CREATE TABLE loci (
        LocusId VARCHAR, Chrom VARCHAR, Motif VARCHAR, CanonicalMotif VARCHAR, gene_region VARCHAR,
        gene_region_rank INTEGER, GeneTableGeneSymbol VARCHAR, IsInMendelianGene INTEGER,
        HPRC256_99thPercentile DOUBLE, FirstAffectedAlleleSize_AllAlleles INTEGER,
        OutlierSampleIds_AllAlleles VARCHAR,
        NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles INTEGER,
        NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles INTEGER
    )""")
    for locus_id, outliers, mendelian, motif, region, p99, num_samples, num_families in LOCI:
        first_allele = int(outliers.split("x")[0]) if outliers else None
        conn.execute("INSERT INTO loci VALUES (?, 'chr1', ?, ?, ?, 1, ?, ?, ?, ?, ?, ?, ?)",
                     [locus_id, motif, motif, region, f"GENE_{locus_id}", mendelian, p99, first_allele,
                      outliers, num_samples, num_families])
    conn.execute("CREATE TABLE per_outlier_phenotype_scores (locus_id VARCHAR, sample_id VARCHAR, "
                 "outlier_type VARCHAR, gene_phenotype_similarity DOUBLE)")
    for locus_id, sample_id, similarity in PHENOTYPE_SCORES:
        conn.execute("INSERT INTO per_outlier_phenotype_scores VALUES (?, ?, 'AllAlleles', ?)",
                     [locus_id, sample_id, similarity])
    if with_outliers_table:
        conn.execute(f"CREATE TABLE {AFFECTED_UNSOLVED_OUTLIERS_TABLE} (LocusId VARCHAR, outlier_type VARCHAR, "
                     f"sample_id VARCHAR, family_id VARCHAR, allele_size INTEGER)")
        for locus_id, sample_id, family_id, allele_size in AFFECTED_UNSOLVED_OUTLIERS:
            conn.execute(f"INSERT INTO {AFFECTED_UNSOLVED_OUTLIERS_TABLE} VALUES (?, 'AllAlleles', ?, ?, ?)",
                         [locus_id, sample_id, family_id, allele_size])
    if population_gate_cohorts:
        write_population_gate_cohorts(conn, population_gate_cohorts)
    conn.commit()
    conn.close()

    sample_table_path = os.path.join(directory, "samples.tsv")
    with open(sample_table_path, "w") as f:
        f.write("sample_id\tfamily_id\tsex\taffected_status\tanalysis_status\tphenotype_description\n")
        for sample_id, (family_id, affected, analysis, phenotype) in SAMPLES.items():
            f.write(f"{sample_id}\t{family_id}\tfemale\t{affected}\t{analysis}\t{phenotype}\n")
    return db_path, sample_table_path


class GroupedViewTestCase(ServerTestCase):
    """Points the server at a test database through configure_app, the shared startup path."""

    WITH_OUTLIERS_TABLE = True
    WITH_SAMPLE_TABLE = True
    POPULATION_GATE_COHORTS = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.db_path, sample_table_path = build_test_database(cls.tmpdir, cls.WITH_OUTLIERS_TABLE,
                                                             cls.POPULATION_GATE_COHORTS)
        cls.configure_server(cls.db_path, sample_table=sample_table_path if cls.WITH_SAMPLE_TABLE else None)

    def _export(self, **params):
        response = self.client.get("/api/v1/outlier_samples/export", query_string=dict({"outlier_type": "all"}, **params))
        self.assertEqual(response.status_code, 200)
        lines = gzip.decompress(response.data).decode().splitlines()
        header = lines[0].split("\t")
        return response, [dict(zip(header, line.split("\t"))) for line in lines[1:]]


class GroupedViewTests(GroupedViewTestCase):

    def _groups(self, **params):
        response = self.client.get("/api/v1/outlier_samples", query_string=dict({"outlier_type": "all"}, **params))
        self.assertEqual(response.status_code, 200, response.get_json())
        return response.get_json()

    def _loci(self, **params):
        response = self.client.get("/api/v1/loci", query_string=dict({"outlier_type": "all"}, **params))
        self.assertEqual(response.status_code, 200, response.get_json())
        return response.get_json()

    def test_group_by_sample(self):
        data = self._groups(group_by="sample")
        self.assertEqual((data["total"], data["total_loci"]), (4, 4))
        self.assertEqual((data["group_sort"], data["group_sort_direction"]), ("num_loci", "desc"))
        rows = {r["group_id"]: r for r in data["results"]}
        self.assertEqual([r["group_id"] for r in data["results"]], ["S1", "S2", "S3", "S6"])
        self.assertEqual((rows["S1"]["num_loci"], rows["S1"]["num_disease_gene_loci"]), (3, 3))
        self.assertEqual((rows["S1"]["largest_allele_size"], rows["S1"]["largest_allele_locus_id"]), (50, "L1"))
        self.assertEqual(rows["S1"]["largest_allele_gene_symbol"], "GENE_L1")
        self.assertEqual(rows["S1"]["largest_allele_delta_population"], 25)
        self.assertIsNone(rows["S6"]["largest_allele_delta_population"])  # L6 has no population data
        self.assertEqual(rows["S1"]["best_gene_phenotype_similarity"], 0.8)
        self.assertIsNone(rows["S2"]["best_gene_phenotype_similarity"])
        self.assertEqual((rows["S6"]["affected_status"], rows["S6"]["analysis_status"]), ("Affected", "Unknown"))
        self.assertEqual((rows["S1"]["family_id"], rows["S1"]["sex"]), ("F1", "female"))
        self.assertEqual(data["motif_counts"], {"AGC": 2, "AT": 1, "CCG": 1})
        self.assertEqual(data["gene_region_counts"], {"CDS": 3, "intron": 1})

    def test_group_by_family(self):
        data = self._groups(group_by="family")
        rows = {r["group_id"]: r for r in data["results"]}
        self.assertEqual({f: r["num_loci"] for f, r in rows.items()}, {"F1": 3, "F2": 1, "F5": 1})
        self.assertEqual([m["sample_id"] for m in rows["F1"]["members"]], ["S1", "S2"])
        self.assertEqual(rows["F1"]["family_id"], "F1")

    def test_paging(self):
        data = self._groups(group_by="sample", page_size=3, page=2)
        self.assertEqual((data["total"], data["total_pages"]), (4, 2))
        self.assertEqual([r["group_id"] for r in data["results"]], ["S6"])

    def test_sorts(self):
        def order(group_sort):
            return [r["group_id"] for r in self._groups(group_by="sample", group_sort=group_sort)["results"]]
        self.assertEqual(order("best_gene_phenotype_similarity"), ["S1", "S3", "S2", "S6"])
        self.assertEqual(order("largest_allele_size"), ["S6", "S3", "S1", "S2"])
        self.assertEqual(order("num_disease_gene_loci"), ["S1", "S2", "S6", "S3"])
        self.assertEqual(order("group_id"), ["S1", "S2", "S3", "S6"])

    def test_sort_directions(self):
        def order(group_sort, direction):
            return [r["group_id"] for r in self._groups(
                group_by="sample", group_sort=group_sort, group_sort_direction=direction)["results"]]
        self.assertEqual(order("num_loci", "asc"), ["S2", "S3", "S6", "S1"])
        self.assertEqual(order("group_id", "desc"), ["S6", "S3", "S2", "S1"])
        self.assertEqual(order("largest_allele_size", "asc"), ["S2", "S1", "S3", "S6"])
        # Samples with no gene-phenotype score stay last either way
        self.assertEqual(order("best_gene_phenotype_similarity", "asc"), ["S3", "S1", "S2", "S6"])
        # Without a direction, each sort uses its own
        self.assertEqual(self._groups(group_by="sample", group_sort="group_id")["group_sort_direction"], "asc")
        response = self.client.get("/api/v1/outlier_samples", query_string={
            "outlier_type": "all", "group_by": "sample", "group_sort_direction": "up"})
        self.assertEqual(response.status_code, 400)

    def test_locus_filters_choose_the_loci(self):
        rows = {r["group_id"]: r["num_loci"] for r in self._groups(group_by="sample", mendelian_genes_only="true")["results"]}
        self.assertEqual(rows, {"S1": 3, "S2": 1, "S6": 1})
        rows = {r["group_id"]: r["num_loci"] for r in self._groups(group_by="sample", motif="CAG")["results"]}
        self.assertEqual(rows, {"S1": 2, "S2": 1})

    def test_sample_level_filters_choose_the_samples(self):
        def samples(**params):
            return [r["group_id"] for r in self._groups(group_by="sample", **params)["results"]]
        self.assertEqual(samples(phenotype_keyword="ATAXIA"), ["S1", "S2"])
        self.assertEqual(samples(sample_id_keyword="s3"), ["S3"])
        self.assertEqual(samples(sample_id_like="S"), ["S1", "S2", "S3", "S6"])
        self.assertEqual(samples(phenotype_keyword="nothing like this"), [])
        # In the family view they pick the families with a sample that passes, with all their loci
        rows = {r["group_id"]: r["num_loci"] for r in self._groups(group_by="family", sample_id_keyword="S2")["results"]}
        self.assertEqual(rows, {"F1": 3})

    def test_each_row_expands_to_exactly_its_loci(self):
        # The page drops the sample-level filters from a row's locus list, since they chose the row
        for group_by, param_name in (("sample", "affected_unsolved_outlier_sample_id"),
                                     ("family", "affected_unsolved_outlier_family_id")):
            for filters in ({}, {"mendelian_genes_only": "true"}, {"motif": "CAG"},
                            {"phenotype_keyword": "ataxia"}, {"sample_id_keyword": "S2", "motif": "CAG"},
                            {"min_repeats_above_population": "20"}):
                locus_filters = {k: v for k, v in filters.items() if k not in results_server.SAMPLE_LEVEL_FILTERS}
                rows = self._groups(group_by=group_by, **filters)["results"]
                self.assertTrue(rows, (group_by, filters))
                for row in rows:
                    loci = self._loci(**locus_filters, **{param_name: row["group_id"]})
                    self.assertEqual(loci["total"], row["num_loci"], (group_by, filters, row["group_id"]))

    def test_expanded_loci_carry_the_selected_alleles(self):
        loci = {r["LocusId"]: r["SelectedSampleOrFamilyAlleles"]
                for r in self._loci(affected_unsolved_outlier_family_id="F1")["results"]}
        self.assertEqual(loci["L1"], [{"sample_id": "S1", "allele_size": 50}, {"sample_id": "S2", "allele_size": 40}])
        self.assertEqual(set(loci), {"L1", "L3", "L6"})
        self.assertNotIn("SelectedSampleOrFamilyAlleles", self._loci()["results"][0])

    def test_loci_carry_the_unsolved_family_count(self):
        counts = {r["LocusId"]: r["NumAffectedUnsolvedFamiliesAboveUnaffected"] for r in self._loci()["results"]}
        self.assertEqual(counts["L1"], 1)
        self.assertEqual(counts["L6"], 2)

    def test_bad_parameters(self):
        for params in ({"group_by": "locus"}, {"group_by": "sample", "group_sort": "size"}, {}):
            response = self.client.get("/api/v1/outlier_samples", query_string=dict({"outlier_type": "all"}, **params))
            self.assertEqual(response.status_code, 400, params)
        response = self.client.get("/api/v1/outlier_samples", query_string={"group_by": "sample"})
        self.assertEqual(response.status_code, 400)  # outlier_type is required, as for /api/v1/loci
        response = self.client.get("/api/v1/outlier_samples/export",
                                   query_string={"outlier_type": "all", "group_by": "sample", "format": "bed"})
        self.assertEqual(response.status_code, 400)

    def test_the_page_names_the_gate_cohorts(self):
        # The test database records no gate, so it reads back as the build's default
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("the population cohorts the database was built with (HPRC256, AoU1027)", page)
        self.assertNotIn("population_gate_cohorts", page)

    def test_schema_advertises_the_view(self):
        data = self.client.get("/api/v1/schema").get_json()
        self.assertTrue(data["has_affected_unsolved_outliers_table"])
        self.assertEqual(data["group_sort_options"], list(results_server.GROUPED_SORTS))
        self.assertEqual(set(data["grouped_view_parameters"]), {"group_by", "group_sort", "group_sort_direction", "format"})
        for name in ("affected_unsolved_outlier_sample_id", "affected_unsolved_outlier_family_id",
                     "min_repeats_above_population"):
            self.assertIn(name, data["filters"])

    def test_export_one_row_per_sample(self):
        response, rows = self._export(group_by="sample", format="groups_tsv")
        self.assertIn("4_samples", response.headers["Content-Disposition"])
        self.assertEqual([(r["sample_id"], r["num_loci"], r["largest_allele_locus_id"]) for r in rows],
                         [("S1", "3", "L1"), ("S2", "1", "L1"), ("S3", "1", "L2"), ("S6", "1", "L6")])
        self.assertEqual(rows[0]["largest_allele_delta_population"], "25.0")
        self.assertEqual(rows[3]["phenotype_description"], "")  # no phenotype
        self.assertEqual(rows[3]["largest_allele_delta_population"], "")  # no population data
        # Every row, not one page, in the chosen order
        _, rows = self._export(group_by="sample", format="groups_tsv", group_sort="largest_allele_size", page_size=1)
        self.assertEqual([r["sample_id"] for r in rows], ["S6", "S3", "S1", "S2"])
        _, rows = self._export(group_by="sample", format="groups_tsv", group_sort="largest_allele_size",
                               group_sort_direction="asc")
        self.assertEqual([r["sample_id"] for r in rows], ["S2", "S1", "S3", "S6"])

    def test_export_one_row_per_family(self):
        response, rows = self._export(group_by="family", format="groups_tsv")
        self.assertIn("3_families", response.headers["Content-Disposition"])
        f1 = rows[0]
        self.assertEqual((f1["family_id"], f1["num_samples"], f1["num_loci"]), ("F1", "2", "3"))
        self.assertEqual(f1["samples"], "S1 (female, Affected/Unsolved); S2 (female, Affected/Unsolved)")
        self.assertEqual(f1["phenotypes"], "seizures, ataxia | ataxia")

    def test_export_one_row_per_sample_and_locus(self):
        response, rows = self._export(group_by="sample", format="group_loci_tsv")
        self.assertIn("6_sample_locus_pairs", response.headers["Content-Disposition"])
        # Groups in the view's order, then the largest allele first
        self.assertEqual([(r["sample_id"], r["LocusId"], r["sample_allele_size"]) for r in rows], [
            ("S1", "L1", "50"), ("S1", "L3", "45"), ("S1", "L6", "20"), ("S2", "L1", "40"),
            ("S3", "L2", "60"), ("S6", "L6", "70")])
        self.assertEqual(rows[0]["sample_allele_delta_population"], "25.0")
        self.assertEqual(rows[0]["NumAffectedUnsolvedFamiliesAboveUnaffected"], "1")
        # Family view: one row per sample of the family, family id first
        _, rows = self._export(group_by="family", format="group_loci_tsv", sample_id_keyword="S2")
        self.assertEqual([(r["family_id"], r["sample_id"], r["LocusId"]) for r in rows],
                         [("F1", "S1", "L1"), ("F1", "S1", "L3"), ("F1", "S2", "L1"), ("F1", "S1", "L6")])

    def test_export_applies_the_filters(self):
        _, rows = self._export(group_by="sample", format="group_loci_tsv", phenotype_keyword="ataxia", mendelian_genes_only="true")
        self.assertEqual({r["sample_id"] for r in rows}, {"S1", "S2"})
        self.assertNotIn("L2", {r["LocusId"] for r in rows})


class ComputeDeltaPopulationTests(unittest.TestCase):
    """The Δpop the grouped view and its export report uses the cohorts and metric the page does."""

    LOCUS_ROW = {"HPRC256_99thPercentile": 20.0, "AoU1027_99thPercentile": 25.0, "TenK10K_99thPercentile": 30.4,
                 "HPRC256_MaxAllele": 40.0, "AoU1027_MaxAllele": float("nan")}

    def _delta(self, **params):
        return results_server.compute_delta_population(50, self.LOCUS_ROW, dict({"outlier_type": "all"}, **params))

    def test_every_cohort_without_a_gate(self):
        self.assertAlmostEqual(self._delta(), 19.6)
        self.assertAlmostEqual(self._delta(require_above_population="all"), 19.6)

    def test_the_gate_cohorts(self):
        self.assertEqual(self._delta(require_above_population="long-read"), 25.0)
        self.assertAlmostEqual(self._delta(require_above_population="short-read"), 19.6)

    def test_the_max_metric_skips_missing_values(self):
        self.assertEqual(self._delta(population_metric="max"), 10.0)

    def test_no_value_or_no_allele(self):
        self.assertIsNone(results_server.compute_delta_population(50, {}, {"outlier_type": "all"}))
        self.assertIsNone(results_server.compute_delta_population(None, self.LOCUS_ROW, {"outlier_type": "all"}))


class NoSampleTableTests(GroupedViewTestCase):
    """Without a sample table every sample lacks a metadata row, as a callset sample missing from the
    sample table does. The sample-id filters and the family ids then come from the outlier table."""

    WITH_SAMPLE_TABLE = False

    def test_sample_id_filters_still_find_the_samples(self):
        def groups(**params):
            response = self.client.get("/api/v1/outlier_samples", query_string=dict({"outlier_type": "all"}, **params))
            self.assertEqual(response.status_code, 200, response.get_json())
            return {r["group_id"]: r for r in response.get_json()["results"]}
        self.assertEqual(list(groups(group_by="sample", sample_id_keyword="s2")), ["S2"])
        self.assertEqual(sorted(groups(group_by="sample", sample_id_like="S")), ["S1", "S2", "S3", "S6"])
        self.assertEqual({f: r["num_loci"] for f, r in groups(group_by="family", sample_id_keyword="S2").items()},
                         {"F1": 3})
        # Phenotypes come only from the sample table, so nothing matches one
        self.assertEqual(groups(group_by="sample", phenotype_keyword="ataxia"), {})
        self.assertEqual(groups(group_by="sample")["S1"]["family_id"], "F1")

    def test_exports_carry_the_recorded_family_ids(self):
        rows = self._export(group_by="family", format="group_loci_tsv")[1]
        self.assertEqual({(r["family_id"], r["sample_id"]) for r in rows},
                         {("F1", "S1"), ("F1", "S2"), ("F2", "S3"), ("F5", "S6")})
        rows = self._export(group_by="sample", format="groups_tsv")[1]
        self.assertEqual({r["sample_id"]: r["family_id"] for r in rows}, {"S1": "F1", "S2": "F1", "S3": "F2", "S6": "F5"})


class RecordedGateTests(GroupedViewTestCase):

    POPULATION_GATE_COHORTS = ("HPRC256", "AoU1027", "TenK10K")

    def test_the_page_names_the_recorded_cohorts(self):
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn("the population cohorts the database was built with (HPRC256, AoU1027, TenK10K)", page)


class MissingTableTests(GroupedViewTestCase):

    WITH_OUTLIERS_TABLE = False

    def test_grouped_view_is_unavailable_and_the_locus_restriction_matches_nothing(self):
        query = {"outlier_type": "all"}
        for path in ("/api/v1/outlier_samples", "/api/v1/outlier_samples/export"):
            response = self.client.get(path, query_string=dict(query, group_by="sample", format="groups_tsv"))
            self.assertEqual(response.status_code, 400)
            self.assertEqual(response.get_json()["error"], "Not available")
        data = self.client.get("/api/v1/loci", query_string=dict(query, affected_unsolved_outlier_sample_id="S1")).get_json()
        self.assertEqual(data["total"], 0)
        self.assertFalse(self.client.get("/api/v1/schema").get_json()["has_affected_unsolved_outliers_table"])


if __name__ == "__main__":
    unittest.main()
