"""Helpers shared by the database test modules.

This module's name does not end in "tests.py", so test discovery imports it only through the test
modules that use it.
"""


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
