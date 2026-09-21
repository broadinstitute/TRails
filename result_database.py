"""Writer for the TRails DuckDB result database.

This module is the persistence layer of the TRails build pipeline. It owns every
table that ends up in the single ``*.duckdb`` file the results server reads, and
nothing else: each function here takes already-computed Python data (lists of
record dicts produced by the in-memory analysis stages) plus an open
``duckdb_compat`` connection, and turns it into tables. No analysis, no TSV
reading, no networking happens here: the orchestrator (``build_database``) calls
these writers in order.

The whole database is built into a temporary path (``<final>.tmp``) and then
atomically moved into place with ``os.replace`` so that a reader never observes
a half-written file and a failed build never clobbers a previous good database.

Tables written:

- ``loci``: the wide per-locus table, one row per locus, columns written in
  ``output_columns`` order (only the columns actually present in the records).
- ``swim_plot``: one row per outlier allele.
- ``per_outlier_phenotype_scores`` / ``per_locus_phenotype_scores``: optional,
  written only when phenotype scoring produced rows.
- ``mendelian_violations`` / ``mendelian_violations_per_motif``: optional,
  written into this same database only when at least one complete trio existed.
- ``metadata``: a small key/value table of build-time facts, currently just the
  number of samples the database was built from.

Column order for the secondary tables is taken from the insertion order of the
keys of the first row dict (Python dicts preserve insertion order), so the
producing modules define the schema and this writer reproduces it faithfully.

All functions are pure aside from the explicit mutations they perform on the
connection they are handed (and, for ``open_new_database`` / ``finalize_database``,
the filesystem move).
"""

import os

import duckdb_compat


# The three outlier types, in the canonical order used throughout TRails.
OUTLIER_TYPES = ("AllAlleles", "ShortAlleles", "HemizygousAlleles")

# The key/value metadata table and the key the sample count is stored under. The
# results server reads it through ``read_sample_count`` to show "N loci, N samples".
METADATA_TABLE = "metadata"
SAMPLE_COUNT_KEY = "num_samples"

# Columns that hold numbers rather than text. DuckDB columns are statically typed, so
# ``write_table_from_dicts`` reads each column's type off the values it is handed; this
# set decides the type only for a column whose rows are all NULL, where there is nothing
# to read. That case is routine rather than exotic: a build run without a gene table
# still writes ``pLI`` and the ``GeneTable*`` columns, all NULL, and the results server
# compares ``pLI`` to a number. Comparing a number against a VARCHAR column is a binder
# error in DuckDB, where SQLite's untyped column simply matched nothing. Any real value
# in the column overrides this set.
#
# The generated names are a superset of the schema (there is no
# ``ThirdUnaffectedAlleleSize_*``, and only AoU1027 carries the OE_* statistics), which
# keeps the patterns readable and costs nothing: a name no column ever uses is never
# looked up. The cohort-statistic patterns also cover the ``HPRC256_Mode`` /
# ``HPRC256_Median`` / ``HPRC256_90thPercentile`` style columns an input matrix can
# supply outside ``OUTPUT_COLUMNS``.
NUMERIC_COLUMNS = frozenset(
    [
        "IsKnownMotif", "IsInMendelianGene", "Start0Based", "End1Based", "MotifSize",
        "gene_region_rank", "NumRepeatsInReference", "pLI",
        "GeneTablepLI_v2", "GeneTablepLI_v4", "GeneTableLoeuf",
        "TRExplorerLocusJaccardSimilarity", "TRExplorerReferenceRepeatPurity",
        # swim_plot's own numeric columns. Its table is created empty (all columns, no
        # rows) when a build finds no outlier alleles, so it hits the same case.
        "outlier_rank", "allele_size", "purity", "methylation",
        "FirstUnaffectedAlleleSize", "is_above_first_unaffected",
    ]
    + [f"{prefix}AffectedAlleleSize_{outlier_type}{by_family}"
       for prefix in ("First", "Second", "Third")
       for outlier_type in OUTLIER_TYPES
       for by_family in ("", "_ByFamily")]
    + [f"{prefix}UnaffectedAlleleSize_{outlier_type}{by_family}"
       for prefix in ("First", "Second")
       for outlier_type in OUTLIER_TYPES
       for by_family in ("", "_ByFamily")]
    + [f"{name}_{outlier_type}"
       for name in ("NumAffectedUnsolvedSamplesAboveUnaffected",
                    "NumAffectedUnsolvedFamiliesAboveUnaffected",
                    "MaxGenePhenoSim", "SumPairwiseSim")
       for outlier_type in OUTLIER_TYPES]
    + [f"{cohort}_{statistic}"
       for cohort in ("AoU1027", "HPRC256", "TenK10K")
       for statistic in ("99thPercentile", "90thPercentile", "MaxAllele", "Mode",
                         "Median", "Stdev", "StdevPercentile", "StdevRankByMotif",
                         "StdevRankTotalNumberByMotif", "OE_Length",
                         "OE_LengthPercentile")])


def open_new_database(db_path):
    """Opens a fresh DuckDB database at a temporary sibling of ``db_path``.

    The database is created at ``db_path + '.tmp'`` (removing any stale temp
    file first) so the final path is only ever populated by an atomic move in
    ``finalize_database``. The temp file's ``.wal`` sidecar goes with it: a build
    killed part way through leaves one behind, and it belongs to a database that
    no longer exists.

    Args:
        db_path: The final database path the build is targeting.

    Returns:
        A ``(connection, tmp_path)`` tuple: an open connection to the temporary
        database and the temporary path it lives at.
    """
    tmp_path = db_path + ".tmp"
    for stale_path in (tmp_path, tmp_path + ".wal"):
        if os.path.exists(stale_path):
            os.remove(stale_path)
    return duckdb_compat.connect(tmp_path), tmp_path


def infer_column_types(column_order, rows):
    """Chooses the DuckDB type to declare each column with, from its values.

    A column is ``BIGINT`` when every value it holds is an int, ``DOUBLE`` when
    every value is numeric and at least one is a float, and ``VARCHAR`` as soon
    as one value is neither (which is what makes a legitimately mixed column,
    e.g. an annotation carrying both counts and a ``"not_available"`` marker,
    round-trip instead of failing its insert). A column with no values at all
    falls back to ``DOUBLE`` for a name in ``NUMERIC_COLUMNS`` and ``VARCHAR``
    otherwise.

    Args:
        column_order: The ordered list of column names being written.
        rows: The list of row dicts that will be inserted.

    Returns:
        A dict mapping each column name to its DuckDB type name.
    """
    # 0 = int, 1 = float, 2 = neither; the widest kind seen wins.
    widest_kind = {}
    for row in rows:
        for column, value in row.items():
            if value is None:
                continue
            kind = 0 if isinstance(value, int) else 1 if isinstance(value, float) else 2
            if kind > widest_kind.get(column, -1):
                widest_kind[column] = kind

    types = {}
    for column in column_order:
        kind = widest_kind.get(column)
        if kind == 0:
            types[column] = "BIGINT"
        elif kind == 1:
            types[column] = "DOUBLE"
        elif kind == 2:
            types[column] = "VARCHAR"
        else:
            types[column] = "DOUBLE" if column in NUMERIC_COLUMNS else "VARCHAR"
    return types


def write_table_from_dicts(connection, table_name, rows, column_order=None,
                           primary_key_column=None, integer_columns=None):
    """Creates ``table_name`` and bulk-inserts a list of row dicts.

    The rows go in through ``duckdb_compat.insert_rows``, a chunk per INSERT
    statement, because a per-row insert of a table as wide as ``loci`` spends
    most of the build's time in the database.

    The column set and order are taken from ``column_order`` when given,
    otherwise from the insertion order of the first row's keys (so the producing
    module's dict layout defines the schema). The table is dropped first, so the
    call is idempotent. ``NULL`` is written for any column absent from a given
    row dict. Column types come from ``infer_column_types`` unless the caller
    names the column explicitly.

    Args:
        connection: An open duckdb_compat connection.
        table_name: The name of the table to (re)create and populate.
        rows: A list of row dicts. May be empty (an empty table is still created
            when ``column_order`` is supplied; otherwise nothing is written).
        column_order: Optional explicit list of column names defining the
            written column order. When omitted, the first row's keys are used.
        primary_key_column: Optional column name to declare ``TEXT PRIMARY KEY``.
        integer_columns: Optional collection of column names to declare
            ``BIGINT`` regardless of the values present.

    Returns:
        The number of rows inserted.
    """
    if column_order is None:
        if not rows:
            return 0
        column_order = list(rows[0].keys())
    integer_columns = set(integer_columns or ())
    column_types = infer_column_types(column_order, rows)

    column_definitions = []
    for column in column_order:
        # Quote the identifier so annotation headers containing spaces or other
        # non-identifier characters (promoted verbatim from the input matrix) do
        # not corrupt the CREATE TABLE DDL.
        if column == primary_key_column:
            definition = f'"{column}" TEXT PRIMARY KEY'
        elif column in integer_columns:
            definition = f'"{column}" BIGINT'
        else:
            definition = f'"{column}" {column_types[column]}'
        column_definitions.append(definition)

    connection.execute(f"DROP TABLE IF EXISTS {table_name}")
    connection.execute(f"CREATE TABLE {table_name} ({', '.join(column_definitions)})")

    if rows:
        duckdb_compat.insert_rows(connection, table_name, column_order, rows)
    return len(rows)


def write_loci_table(connection, records, output_columns, extra_columns=None):
    """Writes the wide ``loci`` table from the in-memory locus records.

    Columns are written in ``output_columns`` order, restricted to the columns
    that are actually present across the records (so an absent annotation column
    is simply not created rather than written as an all-NULL column). A locus
    record missing a present column gets ``NULL`` for it.

    Args:
        connection: An open duckdb_compat connection.
        records: A list of per-locus record dicts (the build's in-memory rows).
        output_columns: The full ordered ``OUTPUT_COLUMNS`` schema; the written
            column order is this list filtered to columns present in ``records``.

    Returns:
        The number of rows written.
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
    # Recognized annotation columns the input supplied that are not part of the
    # fixed OUTPUT_COLUMNS schema (e.g. cohort HPRC256_Mode / _Median /
    # _90thPercentile that the locus-detail view reads) are appended after the
    # schema columns so they are persisted rather than silently dropped.
    if extra_columns:
        for column in sorted(extra_columns):
            if column in present and column not in output_columns:
                present_columns.append(column)

    row_count = write_table_from_dicts(
        connection, "loci", records, column_order=present_columns)
    print(f"  loci: {row_count:,} rows, {len(present_columns)} columns")
    return row_count


def write_sample_count(connection, num_samples):
    """Records how many samples the database was built from, in ``metadata``.

    The table is created if it does not exist yet and the row is written with
    ``INSERT OR REPLACE``, so rebuilding or re-recording the count is idempotent.
    The caller commits.

    Args:
        connection: An open duckdb_compat connection.
        num_samples: The number of samples in the callset, i.e. the number of
            sample genotype columns in the repeat-copy-numbers matrix the build
            read. This counts every sample column in the matrix, including ones
            with no metadata row and ones whose genotypes are all missing; it is
            not the number of rows in the sample-metadata TSV (which may list
            samples the matrix never genotyped).
    """
    connection.execute(
        f"CREATE TABLE IF NOT EXISTS {METADATA_TABLE} "
        f"(key TEXT PRIMARY KEY, value VARCHAR)")
    connection.execute(
        f"INSERT OR REPLACE INTO {METADATA_TABLE} (key, value) VALUES (?, ?)",
        (SAMPLE_COUNT_KEY, str(num_samples)))
    print(f"  {METADATA_TABLE}: {num_samples:,} samples")


def read_sample_count(connection):
    """Returns the sample count stored in ``metadata``, or None when absent.

    Args:
        connection: An open duckdb_compat connection to a results database.

    Returns:
        The recorded number of samples as an int, or ``None`` when the database
        predates the metadata table or has no sample-count row.
    """
    if METADATA_TABLE not in duckdb_compat.list_tables(connection):
        return None
    row = connection.execute(
        f"SELECT value FROM {METADATA_TABLE} WHERE key = ?", (SAMPLE_COUNT_KEY,)).fetchone()
    return int(row[0]) if row else None


def write_swim_plot(connection, swim_rows, columns=None):
    """Writes the ``swim_plot`` table (one row per outlier allele).

    The column schema and order come from the keys of the first swim row (as
    produced by ``swim_plot.generate_swim_plot_table``). When ``swim_rows`` is
    empty but ``columns`` (the canonical ``swim_plot.SWIM_PLOT_COLUMNS`` order) is
    supplied, an empty-but-queryable table is still created so the server's
    swim-plot endpoints return an empty result instead of 500-ing on a missing
    table (a valid build can legitimately yield zero outlier rows).

    Args:
        connection: An open duckdb_compat connection.
        swim_rows: A list of swim-plot row dicts.
        columns: Optional canonical column order used to create the table when
            ``swim_rows`` is empty.

    Returns:
        The number of rows written.
    """
    if not swim_rows and not columns:
        print("  swim_plot: 0 rows (skipped)")
        return 0

    column_order = None if swim_rows else list(columns)
    row_count = write_table_from_dicts(connection, "swim_plot", swim_rows, column_order=column_order)
    connection.commit()

    print(f"  swim_plot: {row_count:,} rows")
    return row_count


def write_phenotype_tables(connection, per_outlier_rows, per_locus_rows):
    """Writes the two phenotype-score tables (when non-empty).

    Both tables are written only when ``per_outlier_rows`` is non-empty (matching
    the reference: the per-locus table is written alongside the per-outlier one).
    Column order for each table is taken from the first row dict.

    Args:
        connection: An open duckdb_compat connection.
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
    per_locus_count = write_table_from_dicts(
        connection, "per_locus_phenotype_scores", per_locus_rows)
    connection.commit()

    print(f"  per_outlier_phenotype_scores: {per_outlier_count:,} rows")
    print(f"  per_locus_phenotype_scores: {per_locus_count:,} rows")
    return per_outlier_count, per_locus_count


def write_mendelian_tables(connection, per_sample_rows, per_motif_rows):
    """Writes the two Mendelian-violation tables into this database (when present).

    These tables are written only when ``per_sample_rows`` is non-empty (i.e. at
    least one complete trio existed). ``sample_id`` is the primary key of each
    table and is forced to be the first column; the remaining columns and their
    order come from the first row dict (all declared ``BIGINT``, matching the
    reference per-trio count schema).

    Args:
        connection: An open duckdb_compat connection.
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
    """Writes one Mendelian table: ``sample_id`` first as a TEXT primary key, every other
    column declared BIGINT.

    Args:
        connection: An open duckdb_compat connection.
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

    Closing before the move matters: DuckDB keeps a ``<path>.wal`` sidecar next to
    an open database and folds it back in on a clean close, so only a closed
    database is the single file that ``os.replace`` can move.

    Uses ``os.replace`` (not a shell ``mv``) so the move is atomic on the same
    filesystem and overwrites any existing final database in a single step.

    A ``<final_path>.wal`` left over from the database that was just overwritten is
    removed after the move. DuckDB does not check that a write-ahead log belongs to the
    database file next to it: it replays whatever it finds, so a log left behind by a
    killed writer (the recovery command the results server prints opens the database
    read-write, and so does the duckdb CLI) would be folded into the brand-new database,
    resurrecting the previous build's tables or, when it re-creates a table the new
    database already has, failing every open with ``Failure while replaying WAL file ...
    Table with name "loci" already exists!``. The removal comes after ``os.replace`` on
    purpose: until the move succeeds, that log is still the missing tail of a database
    that is still there.

    Args:
        connection: The open duckdb_compat connection to the temporary database.
        tmp_path: The temporary database path (from ``open_new_database``).
        final_path: The destination path for the finished database.
    """
    connection.commit()
    connection.close()
    os.replace(tmp_path, final_path)
    if os.path.exists(final_path + ".wal"):
        os.remove(final_path + ".wal")
    print(f"  finalized database -> {final_path}")
