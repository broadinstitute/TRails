#!/usr/bin/env python3
"""Thin DuckDB wrapper shaped like the sqlite3 module TRails used to run on.

TRails is a local app: a user installs it, builds a results database from their own TSVs,
and browses it in the results server. That database is DuckDB, and this module is the only
place in TRails that imports ``duckdb``. Everything else (``build_database``,
``result_database``, ``results_server`` and the tests) keeps calling the sqlite3 DB-API it
was written against, and this module presents that surface on top of DuckDB so three things
stay in one place:

  1. Row objects addressable by both column name and position, the way ``sqlite3.Row`` is.
     DuckDB returns plain tuples, and the query code reads columns by name.
  2. The NULL sort order. DuckDB defaults to NULLS LAST in both directions; SQLite puts
     NULLs first on ASC and last on DESC. Every ORDER BY in TRails was written against the
     SQLite order, so connect() sets default_null_order to match it.
  3. Read-only connections to the same file, which the results server opens once per
     request. DuckDB keeps one instance per file, so those share a cached connection and get
     an independent cursor each, which is DuckDB's supported way to query from threads.

The surface is exactly what TRails calls, and nothing more. Adding a name here that no
caller needs is how this module grows a second, subtly different DB-API, so the list is:

  - ``Error`` / ``OperationalError``: what a failed query raises.
  - ``connect(path, read_only=False)`` -> ``Connection``.
  - ``Connection``: ``execute``, ``executemany``, ``register``, ``unregister``, ``commit``,
    ``close``, and the ``row_factory`` attribute. There is deliberately no context manager:
    ``sqlite3.Connection.__exit__`` commits and leaves the connection OPEN, which is not
    what a reader of ``with duckdb_compat.connect(...)`` would get from a close, so every
    caller uses an explicit ``try: ... finally: conn.close()`` instead.
  - ``Cursor`` (returned by ``execute``, never constructed by a caller): ``description``,
    ``fetchone``, ``fetchall``, ``fetchmany`` and iteration.
  - ``Row``: the ``row_factory`` the results server sets. Indexing by name or position,
    ``keys()``, ``len()`` and iteration over the values, like ``sqlite3.Row``. ``key in
    row`` tests the VALUES there, so this defines no ``__contains__``; the server writes
    ``column in row.keys()`` to test for an optional column.
  - ``list_tables``, ``table_columns``, ``insert_rows``, ``find_leftover_write_files``.

Two SQL dialect differences are NOT hidden here, because hiding them would mean parsing
SQL. They are fixed at each call site instead:

  - LIKE is case-insensitive in SQLite and case-sensitive in DuckDB. Every
    "LIKE ? COLLATE NOCASE" became ILIKE.
  - CAST('abc' AS INTEGER) is 0 in SQLite and an error in DuckDB. Every cast of a text
    column that can hold non-numeric text became TRY_CAST, wrapped in COALESCE(..., 0)
    where the SQLite zero was load-bearing.
"""

import itertools
import os
import threading

import duckdb
import pandas as pd

# duckdb raises a different subclass per problem (BinderException for an unknown column,
# CatalogException for an unknown table, ConversionException for a bad cast), all under
# duckdb.Error. TRails only ever catches "the query did not run", so the callers that used
# to catch sqlite3.OperationalError catch this.
Error = duckdb.Error
OperationalError = duckdb.Error

# Applied to every connection. "nulls_first_on_asc_last_on_desc" is SQLite's ordering.
CONNECTION_CONFIG = {"default_null_order": "nulls_first_on_asc_last_on_desc"}

# Cells (rows x columns) per bulk INSERT in insert_rows, and the name the chunk is
# registered under while that statement runs. The budget is in cells rather than rows
# because what a chunk costs is cells: loci is ~180 columns wide and swim_plot is ~26.
_INSERT_CHUNK_CELLS = 2_000_000
_INSERT_VIEW_NAME = "_trails_insert_rows"

# Rows fetched per batch when iterating a cursor. Iterating has to stream, since the export
# endpoints iterate a query that can match every locus.
_ITERATION_BATCH_ROWS = 2048


def _bound_value(value):
    """Map a NaN parameter to NULL, the way sqlite3 binds it.

    Rows reaching the database writer carry missing numbers as NaN. sqlite3 stored that as
    NULL; DuckDB stores an actual NaN, which is not NULL, is counted by COUNT(col), and
    compares as LARGER than every real number. A locus whose population statistic is missing
    would then pass an "affected allele > population 99th percentile" gate as though it had
    data.

    `value != value` is true only for NaN, and works for numpy floats as well as Python's.
    """
    return None if value != value else value


def _bound_params(params):
    """Return the rows of parameters as a list, each row's NaNs mapped to NULL.

    ``params`` may be any iterable, including a generator, but the result is materialized:
    duckdb 1.0.x rejects a generator with "executemany requires a list of parameter sets to
    be provided", and requirements.txt allows duckdb 1.0. The only caller that matters is
    the notes/tags migration in trails.py, which passes one row per note or tag the user
    typed, so holding them all costs nothing. The bulk loads go through insert_rows, which
    streams a chunk at a time and never comes through here.
    """
    return [[_bound_value(value) for value in row] for row in params]


# One duckdb connection per read-only database file, shared process-wide. DuckDB refuses a
# second connection to a file that is already open with a different configuration, and
# re-opening a multi-gigabyte database per request would be wasteful anyway.
_read_only_connections = {}
_read_only_connections_lock = threading.Lock()


class Row:
    """A query result row addressable by column name or position, like sqlite3.Row.

    Supports row["Motif"], row[0], row.keys(), "Motif" in row.keys(), dict(row), len(row)
    and iteration over the values (all the shapes the TRails query code uses).
    """

    __slots__ = ("_column_indices", "_values")

    def __init__(self, column_indices, values):
        self._column_indices = column_indices
        self._values = values

    def __getitem__(self, key):
        if isinstance(key, str):
            try:
                return self._values[self._column_indices[key]]
            except KeyError:
                raise IndexError(f"No item with that key: {key!r}")
        return self._values[key]

    def keys(self):
        return list(self._column_indices)

    def __iter__(self):
        return iter(self._values)

    def __len__(self):
        return len(self._values)

    def __repr__(self):
        return f"Row({dict(zip(self._column_indices, self._values))!r})"


class Cursor:
    """Wraps a duckdb cursor so fetches honor the connection's row_factory.

    A duckdb connection holds one result at a time, where sqlite3 hands out independent
    cursors. Two consequences:

      - description is read and kept at construction, because the calling code reads it back
        off a cursor after running another query on the same connection.
      - the rows can only be fetched until the next execute() on the same connection.
        Fetching after that raises instead of silently returning the newer query's rows.
    """

    def __init__(self, connection, duckdb_cursor, generation, row_factory):
        self._connection = connection
        self._cursor = duckdb_cursor
        self._generation = generation
        self._row_factory = row_factory
        self.description = duckdb_cursor.description
        self._column_indices = {column[0]: i
                                for i, column in enumerate(self.description or ())}

    def _check_still_current(self):
        if self._connection.generation != self._generation:
            raise RuntimeError(
                "This result has been replaced by a later query on the same connection. "
                "DuckDB keeps one result per connection, so fetch a query's rows before "
                "running the next one, or open a second connection.")

    def _wrap_rows(self, rows):
        if self._row_factory is None:
            return rows
        return [self._row_factory(self._column_indices, row) for row in rows]

    def fetchone(self):
        self._check_still_current()
        values = self._cursor.fetchone()
        if self._row_factory is None or values is None:
            return values
        return self._row_factory(self._column_indices, values)

    def fetchall(self):
        self._check_still_current()
        return self._wrap_rows(self._cursor.fetchall())

    def fetchmany(self, size=1):
        self._check_still_current()
        return self._wrap_rows(self._cursor.fetchmany(size))

    def __iter__(self):
        """Yield rows a batch at a time, the way iterating a sqlite3 cursor does.

        The export endpoints iterate a query that can match every locus and stream each row
        straight into the response, so this must not materialize the whole result.
        """
        while True:
            rows = self.fetchmany(_ITERATION_BATCH_ROWS)
            if not rows:
                return
            for row in rows:
                yield row


class Connection:
    """A duckdb connection with sqlite3's execute/executemany/commit/close interface."""

    def __init__(self, duckdb_connection):
        self._connection = duckdb_connection
        # Bumped by every statement, so a Cursor can tell that its result is no longer the
        # one the connection holds. See Cursor._check_still_current.
        self.generation = 0
        self.row_factory = None

    def execute(self, sql, params=None):
        self.generation += 1
        # duckdb rejects an empty parameter list for a statement with no placeholders, and
        # callers pass () freely, so only forward parameters that exist.
        if params:
            cursor = self._connection.execute(sql, [_bound_value(v) for v in params])
        else:
            cursor = self._connection.execute(sql)
        return Cursor(self, cursor, self.generation, self.row_factory)

    def executemany(self, sql, seq_of_params):
        self.generation += 1
        self._connection.executemany(sql, _bound_params(seq_of_params))

    # Expose a DataFrame to SQL under a name, and take it away again. insert_rows is the
    # caller: it is how a chunk of rows is handed to DuckDB in one statement.
    def register(self, name, df):
        self._connection.register(name, df)

    def unregister(self, name):
        self._connection.unregister(name)

    def commit(self):
        # DuckDB autocommits each statement outside an explicit BEGIN, so there is nothing
        # to flush here. Kept so the build code reads the same as it did on sqlite3.
        pass

    def close(self):
        self._connection.close()


def connect(path, read_only=False):
    """Open a DuckDB database.

    Args:
        path: Path to the .duckdb file.
        read_only: Open without taking the write lock. Read-only connections to the same
            file share one underlying duckdb connection and get their own cursor from it.

    Returns:
        A Connection.
    """
    if not read_only:
        return Connection(duckdb.connect(path, config=CONNECTION_CONFIG))

    key = os.path.abspath(path)
    with _read_only_connections_lock:
        shared = _read_only_connections.get(key)
        if shared is None:
            shared = duckdb.connect(path, read_only=True, config=CONNECTION_CONFIG)
            _read_only_connections[key] = shared
    # A cursor of the shared connection is an independent connection to the same database:
    # safe to use from another thread, and with its own temp schema so two requests can each
    # create a temp table of the same name.
    return Connection(shared.cursor())


def list_tables(conn):
    """Return the set of table names in the database (SQLite's sqlite_master query)."""
    return {row[0] for row in conn.execute("SELECT table_name FROM duckdb_tables()").fetchall()}


def table_columns(conn, table_name):
    """Return the set of column names of a table, or an empty set if it does not exist."""
    return {row[0] for row in conn.execute(
        "SELECT column_name FROM duckdb_columns() WHERE table_name = ?", (table_name,)).fetchall()}


def _varchar_value(value):
    """Render one value the way DuckDB would render it when casting to VARCHAR.

    str() matches DuckDB's own cast for ints, floats (including "6.0", "1e+20" and "inf")
    and anything already a string. A bool is the one value the two disagree on: DuckDB
    writes "true", Python's str() writes "True".
    """
    if value is None or isinstance(value, str):
        return value
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _chunk_column(column_name, values, declared_type):
    """Build one chunk column with the pandas dtype its destination column is declared as.

    Args:
        column_name: The column's name, used only in the error message.
        values: The column's values for this chunk, NaNs already mapped to None.
        declared_type: The DuckDB type the destination column is declared as, or None when
            the table has no such column.

    Returns:
        A pandas array (or object-dtype Series for a type TRails does not declare itself).
    """
    try:
        if declared_type in ("BIGINT", "INTEGER"):
            return pd.array(values, dtype="Int64")
        if declared_type in ("DOUBLE", "FLOAT"):
            return pd.array(values, dtype="Float64")
        if declared_type == "VARCHAR":
            return pd.array([_varchar_value(value) for value in values], dtype="string")
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f'Column "{column_name}" is declared {declared_type} but holds a '
                         f"value that does not fit that type: {error}") from error
    # Any other declared type: hand the values over as they are and let DuckDB convert them.
    return pd.Series(values, dtype=object)


def _declared_column_types(conn, table_name):
    """Return {column name: declared DuckDB type} for an existing table."""
    return {row[0]: row[1] for row in conn.execute(f'DESCRIBE "{table_name}"').fetchall()}

def insert_rows(conn, table_name, column_order, rows):
    """Insert row dicts into an existing table, a chunk per statement.

    executemany runs one prepared statement per row, which is milliseconds per row on a
    table as wide as ``loci``. Registering a chunk as a DataFrame and inserting it with a
    single INSERT ... SELECT loads the same rows 50-90x faster (measured on 300,000
    swim_plot-shaped rows and 30,000 loci-shaped ones), and never holds more than one
    chunk in memory.

    Each column is built with the pandas dtype that matches the type its destination column
    is declared as, read off the table with DESCRIBE, so that nothing about a value is
    guessed. An object-dtype column would be: DuckDB types one from a sample of the rows
    (``pandas_analyze_sample``, 1000 by default) rather than from all of them, while a chunk
    here holds thousands. A column declared VARCHAR because it mixes counts with a
    "not_available" marker would then be typed INT32 off its first rows, and raise on the
    marker or silently round a later 5.5 to 6. Typing the column also keeps a float out of
    an int: left to itself pandas widens a column of ints holding a None to float64, and an
    int written into a VARCHAR column would land as "1.0" instead of the "1" executemany
    wrote.

    Args:
        conn: An open read-write Connection.
        table_name: An existing table, whose columns are ``column_order`` in that order.
        column_order: The ordered column names to read out of each row dict. A column
            absent from a row dict is inserted as NULL.
        rows: An iterable of row dicts. Consumed a chunk at a time, so a generator works.
    """
    declared_types = _declared_column_types(conn, table_name)
    chunk_rows = max(1, _INSERT_CHUNK_CELLS // len(column_order))
    row_iterator = iter(rows)
    while True:
        chunk = list(itertools.islice(row_iterator, chunk_rows))
        if not chunk:
            return
        frame = pd.DataFrame(
            {column: _chunk_column(column,
                                   [_bound_value(row.get(column)) for row in chunk],
                                   declared_types.get(column))
             for column in column_order},
            columns=column_order)
        conn.register(_INSERT_VIEW_NAME, frame)
        try:
            conn.execute(f'INSERT INTO "{table_name}" '
                         f'SELECT * FROM "{_INSERT_VIEW_NAME}"')
        finally:
            conn.unregister(_INSERT_VIEW_NAME)


def find_leftover_write_files(path):
    """Find a leftover write-ahead log sitting next to a DuckDB database.

    DuckDB keeps a "<db>.wal" beside the database while a read-write connection is open and
    removes it on a clean close, so one still sitting there means a process writing to the
    database was killed part way through. Its contents are not in the database file itself
    until a read-write open replays them, so the database is stale until then, and a copy of
    just the "<db>" file is missing whatever the .wal holds.

    Args:
        path: Path to the DuckDB database file.

    Returns:
        List of (filename, size in bytes) tuples for each leftover file that exists.
    """
    wal_path = path + ".wal"
    if not os.path.exists(wal_path):
        return []
    return [(wal_path, os.path.getsize(wal_path))]
