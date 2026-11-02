"""Helpers shared by the database, swim plot and results-server test modules.

This module's name does not end in "tests.py", so test discovery imports it only through the test
modules that use it. It imports nothing, so a test module that uses it does not load the server.
"""


def make_locus_record_with_every_swim_plot_input(**overrides):
    """Returns a minimal locus record, as the build holds it in memory, with every swim-plot input.

    Args:
        **overrides: Values replacing the defaults, or further columns, by column name.

    Returns:
        dict mapping column names to values: one AT locus with no outliers and no population data,
        updated with overrides.
    """
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


def ordered_table_column_names(connection, table_name):
    """Returns the column names of ``table_name`` in table order.

    duckdb_compat.table_columns returns a set; the tests that check column order need the list.

    Args:
        connection: An open duckdb_compat connection.
        table_name: The table to describe.

    Returns:
        The list of column names, in the order the table declares them.
    """
    return [column[0] for column in
            connection.execute(f"SELECT * FROM {table_name} LIMIT 0").description]
