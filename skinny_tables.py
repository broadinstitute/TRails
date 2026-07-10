"""Build per-outlier-type "skinny" sort/filter tables in a TRails results database.

The results web server's locus list+export endpoint filters and sorts the very
wide ``loci`` table (158 columns, ~2KB per row). Any multi-column ``ORDER BY``
cannot be served by a single-column index, so it forces a full scan of the wide
table plus a temporary b-tree sort — slow on a multi-million-row database.

The remedy implemented here is a "deferred row lookup": for each outlier type we
materialize a narrow projection of ``loci`` that holds only the columns the
server actually filters and sorts on. The server runs its multi-key query over
this skinny table, takes the page's ``LocusId`` values, and then fetches the full
wide rows for just those few hundred ids. Because the skinny table keeps the same
column names as ``loci``, the server's existing WHERE/ORDER-BY SQL works verbatim
against it.

Three skinny tables are produced, one per outlier type:
``sk_AllAlleles``, ``sk_ShortAlleles``, and ``sk_HemizygousAlleles``.

The skinny tables are deliberately left without any indexes. A multi-column
``ORDER BY`` cannot be served by a single-column index, and with a selective
``WHERE`` an index on a sort/filter column pushes SQLite toward an index-walk
plan that random-accesses and re-evaluates the expensive predicate across a large
fraction of the table. A single sequential scan of the narrow table followed by a
temporary b-tree sort is faster — that narrow full scan is the whole point. The
outer page fetch joins back to ``loci`` through its own ``LocusId`` index, so the
skinny table needs no ``LocusId`` index either.

All functions are pure aside from the explicit database mutations performed by
``build_skinny_table`` against the connection it is handed.
"""

import time


# The three outlier types for which a skinny table is built.
OUTLIER_TYPES = ("AllAlleles", "ShortAlleles", "HemizygousAlleles")


def shared_columns():
    """Returns the shared (not outlier-type-specific) skinny-table column names.

    These are the 21 columns referenced by the server's filters and sorts that do
    not carry an outlier-type suffix. The order is preserved exactly so the
    projection is stable across runs.

    Returns:
        A list of 21 column-name strings.
    """
    return [
        "LocusId", "Chrom", "gene_id", "Motif", "CanonicalMotif",
        "gene_region", "gene_region_rank", "MotifSize", "NumRepeatsInReference",
        "pLI", "IsKnownMotif", "IsInMendelianGene", "GeneTableGeneSymbol",
        "HPRC256_99thPercentile", "HPRC256_MaxAllele",
        "AoU1027_99thPercentile", "AoU1027_MaxAllele",
        "TenK10K_99thPercentile", "TenK10K_MaxAllele",
        "HPRC256_StdevPercentile", "AoU1027_StdevPercentile",
    ]


def per_outlier_columns(outlier_type):
    """Returns the per-outlier-type skinny-table column names for one outlier type.

    These are the 13 columns referenced by the server's filters and sorts that
    carry an ``_{outlier_type}`` suffix.

    Args:
        outlier_type: One of "AllAlleles", "ShortAlleles", "HemizygousAlleles".

    Returns:
        A list of 13 column-name strings, each suffixed with ``outlier_type``.
    """
    return [
        f"FirstAffectedAlleleSize_{outlier_type}",
        f"SecondAffectedAlleleSize_{outlier_type}",
        f"ThirdAffectedAlleleSize_{outlier_type}",
        f"FirstAffectedAlleleSize_{outlier_type}_ByFamily",
        f"SecondAffectedAlleleSize_{outlier_type}_ByFamily",
        f"ThirdAffectedAlleleSize_{outlier_type}_ByFamily",
        f"FirstUnaffectedAlleleSize_{outlier_type}",
        f"NumAffectedUnsolvedSamplesAboveUnaffected_{outlier_type}",
        f"NumAffectedUnsolvedFamiliesAboveUnaffected_{outlier_type}",
        f"SumPairwiseSim_{outlier_type}",
        f"MaxGenePhenoSim_{outlier_type}",
        f"FirstAffectedPhenotype_{outlier_type}",
        f"OutlierSampleIds_{outlier_type}",
    ]


def build_skinny_table(connection, outlier_type, loci_columns):
    """Creates ``sk_{outlier_type}`` as a narrow projection of ``loci``.

    The projection is the intersection of the wanted columns (shared plus
    per-outlier-type) with the columns that actually exist in ``loci`` — different
    builds (e.g. VCF vs. LPS inputs) carry slightly different column sets, so the
    intersection keeps the build robust.

    For the ``HemizygousAlleles`` table, the hemizygous filter checks
    ``HemizygousAlleleHistogram IS NOT NULL``. Rather than copy that (potentially
    large) blob, a tiny ``1``/``NULL`` marker column of the same name is stored so
    the exact ``IS NOT NULL`` clause still works against the skinny table.

    The table is dropped and recreated, so this is idempotent. No indexes are
    created (see the module docstring for why).

    Args:
        connection: An open sqlite3 connection to the results database.
        outlier_type: One of "AllAlleles", "ShortAlleles", "HemizygousAlleles".
        loci_columns: A collection of the column names present in ``loci``.

    Returns:
        The name of the table that was created (``sk_{outlier_type}``).
    """
    table_name = f"sk_{outlier_type}"
    select_expressions = [
        column for column in shared_columns() + per_outlier_columns(outlier_type)
        if column in loci_columns
    ]
    if outlier_type == "HemizygousAlleles" and "HemizygousAlleleHistogram" in loci_columns:
        select_expressions.append(
            "CASE WHEN HemizygousAlleleHistogram IS NOT NULL THEN 1 END AS HemizygousAlleleHistogram")

    connection.execute(f"DROP TABLE IF EXISTS {table_name}")
    start_time = time.time()
    connection.execute(
        f"CREATE TABLE {table_name} AS SELECT {', '.join(select_expressions)} FROM loci")
    row_count = connection.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
    print(f"  {table_name}: {row_count:,} rows, {len(select_expressions)} columns "
          f"({time.time() - start_time:.1f}s)")
    return table_name
