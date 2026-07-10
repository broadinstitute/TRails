"""Pure SQLite writer for the TRails result database.

This module is the persistence layer of the TRails build pipeline. It owns every
table that ends up in the single ``*.db`` file the results server reads, and
nothing else: each function here takes already-computed Python data (lists of
record dicts produced by the in-memory analysis stages) plus an open
``sqlite3.Connection``, and turns it into tables and indexes. No analysis, no
TSV reading, no networking happens here — the orchestrator (``build_database``)
calls these writers in order.

The whole database is built into a temporary path (``<final>.tmp``) and then
atomically moved into place with ``os.replace`` so that a reader never observes
a half-written file and a failed build never clobbers a previous good database.

Tables written (mirroring the reference ``analyze_results.py`` SQLite output):

- ``loci``: the wide per-locus table, one row per locus, columns written in
  ``output_columns`` order (only the columns actually present in the records).
- ``swim_plot``: one row per outlier allele, with its 13 supporting indexes.
- ``sk_AllAlleles`` / ``sk_ShortAlleles`` / ``sk_HemizygousAlleles``: the narrow
  per-outlier-type sort/filter projections (delegated to ``skinny_tables``).
- ``per_outlier_phenotype_scores`` / ``per_locus_phenotype_scores``: optional,
  written only when phenotype scoring produced rows.
- ``mendelian_violations`` / ``mendelian_violations_per_motif``: optional,
  written into this same database only when at least one complete trio existed.

Column order for the secondary tables is taken from the insertion order of the
keys of the first row dict (Python dicts preserve insertion order), so the
producing modules define the schema and this writer reproduces it faithfully.

All functions are pure aside from the explicit mutations they perform on the
connection they are handed (and, for ``open_new_database`` / ``finalize_database``,
the filesystem move).
"""

import os
import sqlite3

import skinny_tables


# The three outlier types, in the canonical order used throughout TRails.
OUTLIER_TYPES = ("AllAlleles", "ShortAlleles", "HemizygousAlleles")


def open_new_database(db_path):
    """Opens a fresh SQLite database at a temporary sibling of ``db_path``.

    The database is created at ``db_path + '.tmp'`` (removing any stale temp
    file first) so the final path is only ever populated by an atomic move in
    ``finalize_database``.

    Args:
        db_path: The final database path the build is targeting.

    Returns:
        A ``(connection, tmp_path)`` tuple: an open sqlite3 connection to the
        temporary database and the temporary path it lives at.
    """
    tmp_path = db_path + ".tmp"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    return sqlite3.connect(tmp_path), tmp_path


def write_table_from_dicts(connection, table_name, rows, column_order=None,
                           primary_key_column=None, integer_columns=None):
    """Creates ``table_name`` and bulk-inserts a list of row dicts.

    The column set and order are taken from ``column_order`` when given,
    otherwise from the insertion order of the first row's keys (so the producing
    module's dict layout defines the schema). The table is dropped first, so the
    call is idempotent. ``NULL`` is written for any column absent from a given
    row dict.

    Args:
        connection: An open sqlite3 connection.
        table_name: The name of the table to (re)create and populate.
        rows: A list of row dicts. May be empty (an empty table is still created
            when ``column_order`` is supplied; otherwise nothing is written).
        column_order: Optional explicit list of column names defining the
            written column order. When omitted, the first row's keys are used.
        primary_key_column: Optional column name to declare ``PRIMARY KEY``.
        integer_columns: Optional collection of column names to declare with an
            ``INTEGER`` type affinity (the rest are created untyped, matching the
            dynamically typed values pandas' ``to_sql`` would store).

    Returns:
        The number of rows inserted.
    """
    if column_order is None:
        if not rows:
            return 0
        column_order = list(rows[0].keys())
    integer_columns = set(integer_columns or ())

    column_definitions = []
    for column in column_order:
        definition = column
        if column == primary_key_column:
            definition += " TEXT PRIMARY KEY"
        elif column in integer_columns:
            definition += " INTEGER"
        column_definitions.append(definition)

    connection.execute(f"DROP TABLE IF EXISTS {table_name}")
    connection.execute(f"CREATE TABLE {table_name} ({', '.join(column_definitions)})")

    placeholders = ", ".join("?" for _ in column_order)
    connection.executemany(
        f"INSERT INTO {table_name} VALUES ({placeholders})",
        [tuple(row.get(column) for column in column_order) for row in rows])
    return len(rows)


def write_loci_table(connection, records, output_columns):
    """Writes the wide ``loci`` table from the in-memory locus records.

    Columns are written in ``output_columns`` order, restricted to the columns
    that are actually present across the records (so an absent annotation column
    is simply not created rather than written as an all-NULL column). A locus
    record missing a present column gets ``NULL`` for it.

    Args:
        connection: An open sqlite3 connection.
        records: A list of per-locus record dicts (the build's in-memory rows).
        output_columns: The full ordered ``OUTPUT_COLUMNS`` schema; the written
            column order is this list filtered to columns present in ``records``.

    Returns:
        A ``(row_count, present_columns)`` tuple, where ``present_columns`` is the
        ordered list of columns actually written.
    """
    present = set()
    for record in records:
        present.update(record.keys())
    present_columns = [column for column in output_columns if column in present]
    # With no records there are no present columns; fall back to the full schema
    # so an empty cohort (e.g. -n 0 or a header-only matrix) still produces a
    # valid, queryable loci table instead of the invalid DDL `CREATE TABLE loci ()`.
    if not present_columns:
        present_columns = list(output_columns)

    row_count = write_table_from_dicts(
        connection, "loci", records, column_order=present_columns)
    print(f"  loci: {row_count:,} rows, {len(present_columns)} columns")
    return row_count, present_columns


def create_loci_indexes(connection, present_columns):
    """Creates the ``loci`` indexes the results server relies on.

    Mirrors the index set built by the reference pipeline, but every index is
    guarded by the presence of its column so a minimal database (without the
    optional gene / phenotype / population-stat columns) still indexes cleanly.

    Args:
        connection: An open sqlite3 connection holding a populated ``loci`` table.
        present_columns: The collection of column names present in ``loci``.

    Returns:
        The list of index names that were created.
    """
    present = set(present_columns)
    cursor = connection.cursor()
    created = []

    def add_index(index_name, column):
        if column in present:
            cursor.execute(
                f"CREATE INDEX IF NOT EXISTS {index_name} ON loci({column})")
            created.append(index_name)

    # Primary lookup + common single-column filter indexes.
    for column in [
        "LocusId", "gene_id", "Chrom", "KnownDiseaseLocus", "gene_region",
        "gene_region_rank", "pLI", "Motif", "CanonicalMotif", "IsKnownMotif",
        "IsInMendelianGene", "NumRepeatsInReference", "HPRC256_StdevPercentile",
        "AoU1027_StdevPercentile",
    ]:
        add_index(f"idx_loci_{column}", column)

    # Phenotype score indexes (present only when phenotype scoring ran).
    for outlier_type in OUTLIER_TYPES:
        add_index(f"idx_loci_MaxGenePhenoSim_{outlier_type}",
                  f"MaxGenePhenoSim_{outlier_type}")
        add_index(f"idx_loci_SumPairwiseSim_{outlier_type}",
                  f"SumPairwiseSim_{outlier_type}")

    # Affected / unaffected allele-size indexes for the above-unaffected filters.
    for prefix in ["First", "Second", "Third"]:
        for outlier_type in OUTLIER_TYPES:
            add_index(f"idx_loci_{prefix}AffectedAlleleSize_{outlier_type}",
                      f"{prefix}AffectedAlleleSize_{outlier_type}")
    for outlier_type in OUTLIER_TYPES:
        add_index(f"idx_loci_FirstUnaffectedAlleleSize_{outlier_type}",
                  f"FirstUnaffectedAlleleSize_{outlier_type}")

    # Family-level allele-size indexes.
    for prefix in ["First", "Second", "Third"]:
        for outlier_type in OUTLIER_TYPES:
            add_index(
                f"idx_loci_{prefix}AffectedAlleleSize_{outlier_type}_ByFamily",
                f"{prefix}AffectedAlleleSize_{outlier_type}_ByFamily")
    for outlier_type in OUTLIER_TYPES:
        add_index(
            f"idx_loci_FirstUnaffectedAlleleSize_{outlier_type}_ByFamily",
            f"FirstUnaffectedAlleleSize_{outlier_type}_ByFamily")

    # Sorting indexes used by the default ranked views.
    add_index("idx_loci_NumAffectedAboveUnaffected_AllAlleles",
              "NumAffectedUnsolvedSamplesAboveUnaffected_AllAlleles")
    add_index("idx_loci_NumAffectedFamiliesAboveUnaffected_AllAlleles",
              "NumAffectedUnsolvedFamiliesAboveUnaffected_AllAlleles")

    connection.commit()
    return created


def write_swim_plot(connection, swim_rows):
    """Writes the ``swim_plot`` table (one row per outlier allele) + its indexes.

    The 38-column schema and its order come from the keys of the first swim row
    (as produced by ``swim_plot.generate_swim_plot_table``). If ``swim_rows`` is
    empty the table is not created.

    Args:
        connection: An open sqlite3 connection.
        swim_rows: A list of swim-plot row dicts.

    Returns:
        The number of rows written.
    """
    if not swim_rows:
        print("  swim_plot: 0 rows (skipped)")
        return 0

    row_count = write_table_from_dicts(connection, "swim_plot", swim_rows)

    cursor = connection.cursor()
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_outlier_category ON swim_plot(outlier_type, motif_category)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_allele_size ON swim_plot(allele_size)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_sample_outlier ON swim_plot(sample_id, outlier_type, is_above_first_unaffected)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_locus_id ON swim_plot(LocusId, outlier_type)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_gene_symbol ON swim_plot(GeneTableGeneSymbol, outlier_type)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_gene_region ON swim_plot(gene_region)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_outlier_rank ON swim_plot(outlier_type, outlier_rank, MotifSize)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_known_motif ON swim_plot(IsKnownMotif)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_canonical_motif ON swim_plot(CanonicalMotif)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_pli ON swim_plot(pLI)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_motif_size ON swim_plot(MotifSize)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_gene_id ON swim_plot(gene_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_swim_mendelian ON swim_plot(IsInMendelianGene)")
    connection.commit()

    print(f"  swim_plot: {row_count:,} rows, 13 indexes")
    return row_count


def write_skinny_tables(connection, present_columns):
    """Builds the three per-outlier-type skinny sort/filter tables.

    Delegates to ``skinny_tables.build_skinny_table`` for each outlier type. The
    skinny tables intersect their wanted columns with ``present_columns`` so a
    minimal ``loci`` table still yields valid (narrower) projections.

    Args:
        connection: An open sqlite3 connection holding a populated ``loci`` table.
        present_columns: The collection of column names present in ``loci``.

    Returns:
        The list of skinny table names that were created.
    """
    present = set(present_columns)
    table_names = [
        skinny_tables.build_skinny_table(connection, outlier_type, present)
        for outlier_type in OUTLIER_TYPES
    ]
    connection.commit()
    return table_names


def write_phenotype_tables(connection, per_outlier_rows, per_locus_rows):
    """Writes the two phenotype-score tables + their indexes (when non-empty).

    Both tables are written only when ``per_outlier_rows`` is non-empty (matching
    the reference: the per-locus table is written alongside the per-outlier one).
    Column order for each table is taken from the first row dict.

    Args:
        connection: An open sqlite3 connection.
        per_outlier_rows: A list of per-outlier phenotype-score row dicts.
        per_locus_rows: A list of per-locus phenotype-score row dicts.

    Returns:
        A ``(per_outlier_count, per_locus_count)`` tuple. ``(0, 0)`` when skipped.
    """
    if not per_outlier_rows:
        print("  phenotype tables: skipped (no rows)")
        return 0, 0

    per_outlier_count = write_table_from_dicts(
        connection, "per_outlier_phenotype_scores", per_outlier_rows)
    cursor = connection.cursor()
    cursor.execute("CREATE INDEX idx_outlier_pheno_locus ON per_outlier_phenotype_scores(locus_id)")
    cursor.execute("CREATE INDEX idx_outlier_pheno_sample ON per_outlier_phenotype_scores(sample_id)")
    cursor.execute("CREATE INDEX idx_outlier_pheno_type ON per_outlier_phenotype_scores(outlier_type)")
    cursor.execute("CREATE INDEX idx_outlier_pheno_gene_sim ON per_outlier_phenotype_scores(gene_phenotype_similarity)")
    cursor.execute("CREATE INDEX idx_outlier_pheno_gene ON per_outlier_phenotype_scores(gene_symbol)")

    per_locus_count = write_table_from_dicts(
        connection, "per_locus_phenotype_scores", per_locus_rows)
    cursor.execute("CREATE INDEX idx_locus_pheno_id ON per_locus_phenotype_scores(locus_id)")
    cursor.execute("CREATE INDEX idx_locus_pheno_type ON per_locus_phenotype_scores(outlier_type)")
    cursor.execute("CREATE INDEX idx_locus_pheno_sum_sim ON per_locus_phenotype_scores(sum_pairwise_similarity)")
    cursor.execute("CREATE INDEX idx_locus_pheno_max_gene ON per_locus_phenotype_scores(max_gene_phenotype_similarity)")
    cursor.execute("CREATE INDEX idx_locus_pheno_num_samples ON per_locus_phenotype_scores(num_qualifying_samples)")
    cursor.execute("CREATE INDEX idx_locus_pheno_known_motif ON per_locus_phenotype_scores(IsKnownMotif)")
    cursor.execute("CREATE INDEX idx_locus_pheno_gene_region ON per_locus_phenotype_scores(gene_region)")
    cursor.execute("CREATE INDEX idx_locus_pheno_motif_size ON per_locus_phenotype_scores(MotifSize)")
    cursor.execute("CREATE INDEX idx_locus_pheno_pairwise_sort ON per_locus_phenotype_scores(outlier_type, sum_pairwise_similarity DESC, locus_id)")
    cursor.execute("CREATE INDEX idx_locus_pheno_gene_sort ON per_locus_phenotype_scores(outlier_type, max_gene_phenotype_similarity DESC, locus_id)")
    cursor.execute("CREATE INDEX idx_locus_pheno_gene_known ON per_locus_phenotype_scores(outlier_type, IsKnownMotif, max_gene_phenotype_similarity DESC, locus_id)")
    connection.commit()

    print(f"  per_outlier_phenotype_scores: {per_outlier_count:,} rows")
    print(f"  per_locus_phenotype_scores: {per_locus_count:,} rows")
    return per_outlier_count, per_locus_count


def write_mendelian_tables(connection, per_sample_rows, per_motif_rows):
    """Writes the two Mendelian-violation tables into this database (when present).

    These tables are written only when ``per_sample_rows`` is non-empty (i.e. at
    least one complete trio existed). ``sample_id`` is the primary key of each
    table and is forced to be the first column; the remaining columns and their
    order come from the first row dict (all declared ``INTEGER``, matching the
    reference per-trio count schema).

    Args:
        connection: An open sqlite3 connection.
        per_sample_rows: A list of per-trio-child row dicts for
            ``mendelian_violations``.
        per_motif_rows: A list of per-trio-child row dicts for
            ``mendelian_violations_per_motif``.

    Returns:
        A ``(per_sample_count, per_motif_count)`` tuple. ``(0, 0)`` when skipped.
    """
    if not per_sample_rows:
        print("  mendelian tables: skipped (no trios)")
        return 0, 0

    per_sample_count = _write_mendelian_table(
        connection, "mendelian_violations", per_sample_rows)
    per_motif_count = _write_mendelian_table(
        connection, "mendelian_violations_per_motif", per_motif_rows)
    connection.commit()

    print(f"  mendelian_violations: {per_sample_count:,} rows")
    print(f"  mendelian_violations_per_motif: {per_motif_count:,} rows")
    return per_sample_count, per_motif_count


def _write_mendelian_table(connection, table_name, rows):
    """Writes one Mendelian table with ``sample_id`` first as INTEGER-typed PK.

    Args:
        connection: An open sqlite3 connection.
        table_name: The table to (re)create.
        rows: A list of flat row dicts; every key other than ``sample_id`` holds
            an integer count.

    Returns:
        The number of rows written.
    """
    column_order = ["sample_id"] + [
        column for column in rows[0].keys() if column != "sample_id"]
    return write_table_from_dicts(
        connection, table_name, rows, column_order=column_order,
        primary_key_column="sample_id",
        integer_columns=[c for c in column_order if c != "sample_id"])


def finalize_database(connection, tmp_path, final_path):
    """Commits, closes, and atomically moves the temp database into place.

    Uses ``os.replace`` (not a shell ``mv``) so the move is atomic on the same
    filesystem and overwrites any existing final database in a single step.

    Args:
        connection: The open sqlite3 connection to the temporary database.
        tmp_path: The temporary database path (from ``open_new_database``).
        final_path: The destination path for the finished database.
    """
    connection.commit()
    connection.close()
    os.replace(tmp_path, final_path)
    print(f"  finalized database -> {final_path}")
