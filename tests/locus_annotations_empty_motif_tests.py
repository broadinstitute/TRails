"""Unit tests for locus_annotations.add_derived_locus_columns on a locus with an empty motif."""

import unittest

import locus_annotations


class EmptyMotifTests(unittest.TestCase):
    """An empty motif gives MotifSize 0 and no NumRepeatsInReference, not a division by zero."""

    def test_empty_motif_does_not_crash(self):
        records = [{"LocusId": "1-100-115-", "Motif": ""}]
        locus_annotations.add_derived_locus_columns(records)
        self.assertEqual(records[0]["MotifSize"], 0)
        self.assertIsNone(records[0]["NumRepeatsInReference"])


if __name__ == "__main__":
    unittest.main()
