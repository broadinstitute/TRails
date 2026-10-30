"""Unit tests for mendelian_qc.compute_mendelian_violations on a locus with a blank motif."""

import unittest

import mendelian_qc


class BlankMotifMendelianTests(unittest.TestCase):
    """A blank motif must not raise KeyError and abort the build; the locus is skipped."""

    def test_blank_motif_does_not_crash(self):
        locus_rows = [{
            "trid": "chr1-100-110-",
            "motif": "",
            "genotypes": {"CHILD": "10,20", "MOM": "10,11", "DAD": "20,21"},
        }]
        sample_lookup = {
            "CHILD": {"sample_id": "CHILD", "maternal_id": "MOM", "paternal_id": "DAD"},
            "MOM": {"sample_id": "MOM", "maternal_id": "", "paternal_id": ""},
            "DAD": {"sample_id": "DAD", "maternal_id": "", "paternal_id": ""},
        }
        per_sample, _ = mendelian_qc.compute_mendelian_violations(locus_rows, sample_lookup)
        # The blank-motif locus is skipped, so the trio tallies zero loci.
        self.assertEqual(per_sample[0]["total_loci"], 0)


if __name__ == "__main__":
    unittest.main()
