import unittest

import make_as_bed as m

LENGTHS = {"chr1": 10000, "chr2": 5000}


class TestMakeAsBed(unittest.TestCase):

    def test_pad_merge_and_names(self):
        rows = [("chr1", 1000, 1001, "DPYD"), ("chr1", 1500, 1600, "DPYD"), ("chr1", 3000, 3001, "MTHFR")]
        self.assertEqual(m.build(rows, 300, LENGTHS),
                         [("chr1", 700, 1900, "DPYD"), ("chr1", 2700, 3301, "MTHFR")])
        self.assertEqual(m.build(rows, 1000, LENGTHS), [("chr1", 0, 4001, "DPYD,MTHFR")])

    def test_clipped_to_chromosome_and_fai_order(self):
        rows = [("chr2", 4900, 4950, "B"), ("chr1", 50, 60, "A")]
        self.assertEqual(m.build(rows, 500, LENGTHS), [("chr1", 0, 560, "A"), ("chr2", 4400, 5000, "B")])

    def test_unknown_contig_stops(self):
        with self.assertRaises(SystemExit):
            m.build([("22", 1, 2, "CYP2D6")], 0, LENGTHS)


if __name__ == "__main__":
    unittest.main()
