import unittest
from tempfile import NamedTemporaryFile

import panel_coverage as pc

ONT = [100, 50, 30, 1]


class TestPanelCoverage(unittest.TestCase):

    def test_depth_colour_boundaries(self):
        cases = [(150, "green"), (100, "green"), (99, "palegreen"), (50, "palegreen"),
                 (49, "orange"), (30, "orange"), (29, "red"), (1, "red"), (0, "black")]
        for depth, expected in cases:
            self.assertEqual(pc.depth_colour(depth, ONT), expected, depth)

    def test_worst_colour(self):
        self.assertEqual(pc.worst_colour(["green", "orange", "palegreen"]), "orange")
        self.assertEqual(pc.worst_colour(["red", "black"]), "black")
        self.assertEqual(pc.worst_colour([]), "")

    def test_targets_and_depths(self):
        bed = NamedTemporaryFile("w", suffix=".bed", delete=False)
        bed.write("chr1\t10\t13\tTPMT|rs1\nchr1\t20\t21\tDPYD\nchr2\t5\t6\tNUDT15\n")
        bed.close()
        targets = pc.read_targets(bed.name)
        # BED is 0-based half-open: chr1:10-13 covers positions 11, 12, 13
        depth = ["chr1\t10\t500", "chr1\t11\t120", "chr1\t12\t40", "chr1\t13\t200", "chr1\t21\t0"]
        pc.add_depths(targets, depth)
        rows = {r["gene"]: r for r in pc.summarise(targets, ONT)}
        self.assertEqual(rows["TPMT"]["rsid"], "rs1")
        self.assertEqual(rows["TPMT"]["min_depth"], 40)
        self.assertEqual(rows["TPMT"]["mean_depth"], 120.0)
        self.assertEqual(rows["TPMT"]["colour"], "orange")
        self.assertEqual(rows["DPYD"]["colour"], "black")
        # no depth lines at all for a target counts as 0 reads
        self.assertEqual(rows["NUDT15"]["min_depth"], 0)

    def test_clinical_positions_decide_colour(self):
        bed = NamedTemporaryFile("w", suffix=".bed", delete=False)
        bed.write("chr1\t100\t110\tTPMT\nchr1\t200\t205\tHLA-B\n")
        bed.close()
        positions = NamedTemporaryFile("w", suffix=".bed", delete=False)
        positions.write("chr1\t104\t105\tTPMT|rs1\nchr1\t500\t501\tX|rs2\n")
        positions.close()
        targets = pc.read_targets(bed.name)
        pc.add_clinical_positions(targets, positions.name)
        # weak edge base (pos 101: 2 reads), clinical position 105 has 120 reads
        depth = [f"chr1\t{p}\t{2 if p == 101 else 120}" for p in range(101, 111)]
        depth += [f"chr1\t{p}\t{3 if p == 203 else 80}" for p in range(201, 206)]
        pc.add_depths(targets, depth)
        rows = {r["gene"]: r for r in pc.summarise(targets, ONT)}
        self.assertEqual(rows["TPMT"]["n_clinical_positions"], 1)
        self.assertEqual(rows["TPMT"]["basis"], "kliniska positioner")
        self.assertEqual(rows["TPMT"]["min_depth"], 2)
        self.assertEqual(rows["TPMT"]["colour"], "green")
        self.assertEqual(rows["TPMT"]["pct_ge_green"], 90.0)
        # no clinical positions inside -> worst base of the whole target
        self.assertEqual(rows["HLA-B"]["basis"], "hela målet")
        self.assertEqual(rows["HLA-B"]["colour"], "red")


if __name__ == "__main__":
    unittest.main()
