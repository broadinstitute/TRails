"""Unit tests for analysis_columns.add_gene_columns, which fills the GeneTable* columns."""

import unittest

import analysis_columns


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
