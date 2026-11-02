"""Temporary input files for the test modules.

This module's name does not end in "tests.py", so test discovery imports it only through the test
modules that use it.
"""

import gzip
import os
import tempfile


def write_temporary_file(test_case, text, suffix=".tsv"):
    """Write text to a temporary file that is removed when test_case finishes.

    Args:
        test_case: The running unittest.TestCase, which registers the cleanup.
        text: The file contents.
        suffix: The file name suffix; a name ending in ".gz" is written gzipped.

    Returns:
        The file's path.
    """
    handle = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    handle.close()
    test_case.addCleanup(os.remove, handle.name)
    with (gzip.open if suffix.endswith(".gz") else open)(handle.name, "wt") as output_file:
        output_file.write(text)
    return handle.name
