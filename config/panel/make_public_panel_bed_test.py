import json
import unittest
from tempfile import NamedTemporaryFile

import make_public_panel_bed as m

LENGTHS = {"chr1": 10000, "chr22": 50000}


def _file(text, suffix):
    fh = NamedTemporaryFile("w", suffix=suffix, delete=False)
    fh.write(text)
    fh.close()
    return fh.name


class TestPublicPanelBed(unittest.TestCase):

    def test_pharmcat_vcf_names_and_indel_length(self):
        vcf = _file("##x\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr1\t100\trs1\tG\tT\t.\tPASS\tPX=DPYD\n"
                    "chr1\t200\trs2\tTG\tT\t.\tPASS\tPX=DPYD\n"
                    "chr1\t900\trs3\tG\tA\t.\tPASS\tPOI\n", ".vcf")
        self.assertEqual(m.read_pharmcat_vcf(vcf),
                         [("chr1", 99, 100, "DPYD"), ("chr1", 199, 201, "DPYD"), ("chr1", 899, 900, "POI|rs3")])

    def test_merge_only_within_a_gene(self):
        rows = [("chr1", 99, 100, "DPYD"), ("chr1", 199, 201, "DPYD"), ("chr1", 250, 251, "TPMT")]
        self.assertEqual(m.build(rows, 100, LENGTHS),
                         [("chr1", 0, 301, "DPYD"), ("chr1", 150, 351, "TPMT")])

    def test_cyp2d6_locus_from_starphase(self):
        db = _file(json.dumps({"cyp2d6_config": {"cyp_coordinates": {
            "CYP2D6": {"chrom": "chr22", "start": 26259, "end": 32424},
            "REP6": {"chrom": "chr22", "start": 23191, "end": 25963},
            "CYP2D7": {"chrom": "chr22", "start": 39965, "end": 45903}}}}), ".json")
        self.assertEqual(m.read_cyp2d6_locus(db), [("chr22", 23191, 45903, "CYP2D6")])


if __name__ == "__main__":
    unittest.main()
