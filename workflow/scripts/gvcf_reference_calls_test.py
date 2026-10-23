import gzip
import unittest
from tempfile import NamedTemporaryFile

import gvcf_reference_calls as grc

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "##contig=<ID=chr2,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
)


def _gz(text):
    fh = NamedTemporaryFile(suffix=".vcf.gz", delete=False)
    fh.close()
    with gzip.open(fh.name, "wt") as out:
        out.write(text)
    return fh.name


class TestGvcfReferenceCalls(unittest.TestCase):

    def setUp(self):
        # variant VCF: one phased het SNV at chr1:10
        self.vcf = _gz(HEADER + "chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT:PS\t0|1:10\n")
        # expanded gVCF: per-base reference calls with MIN_DP, one low, one ./.
        gvcf = HEADER + "".join([
            "chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT:DP\t0/1:80\n",
            "chr1\t20\t.\tC\t.\t0\t.\t.\tGT:MIN_DP\t0/0:80\n",
            "chr1\t21\t.\tT\t.\t0\t.\t.\tGT:MIN_DP\t0/0:80\n",
            "chr1\t22\t.\tG\t.\t0\t.\t.\tGT:MIN_DP\t0/0:12\n",
            "chr1\t30\t.\tA\t.\t0\t.\t.\tGT:MIN_DP\t./.:60\n",
            "chr2\t5\t.\tC\t<*>\t0\t.\t.\tGT:MIN_DP\t0/0:45\n",
        ])
        self.bases = grc.read_gvcf_bases(_gz(gvcf), {"chr1", "chr2"})

    def _merge(self, positions):
        out = NamedTemporaryFile("w", suffix=".vcf", delete=False)
        out.close()
        n = grc.merge(self.vcf, self.bases, positions, 30, out.name)
        with open(out.name) as fh:
            text = fh.read()
        records = [l.split("\t") for l in text.splitlines() if not l.startswith("#")]
        return n, {(r[0], r[1]): r for r in records}, text

    def test_reference_depth_rule(self):
        n, rec, _ = self._merge({("chr1", 19, 20), ("chr1", 21, 22), ("chr1", 29, 30), ("chr2", 4, 5)})
        self.assertEqual(rec[("chr1", "20")][9].strip(), "0/0:80")
        self.assertEqual(rec[("chr1", "22")][9].strip(), "./.:12")   # depth below 30
        self.assertEqual(rec[("chr1", "30")][9].strip(), "./.:60")   # no confident call
        self.assertEqual(rec[("chr2", "5")][9].strip(), "0/0:45")
        self.assertEqual(n, 4)

    def test_variant_record_and_phasing_kept(self):
        _, rec, _ = self._merge({("chr1", 9, 10)})
        self.assertEqual(rec[("chr1", "10")][9].strip(), "0|1:10")
        self.assertEqual(len(rec), 1)

    def test_multibase_position(self):
        # chr1:20-22 (e.g. a deletion site): one weak base -> missing
        _, rec, _ = self._merge({("chr1", 19, 22)})
        self.assertEqual(rec[("chr1", "20")][3], "CTG")
        self.assertEqual(rec[("chr1", "20")][9].strip(), "./.:12")
        _, rec, _ = self._merge({("chr1", 19, 21)})
        self.assertEqual(rec[("chr1", "20")][3], "CT")
        self.assertEqual(rec[("chr1", "20")][9].strip(), "0/0:80")

    def test_position_missing_from_gvcf_is_skipped(self):
        n, rec, text = self._merge({("chr1", 499, 500)})
        self.assertEqual(n, 0)
        self.assertIn("##pgx_reference_calls=<min_depth=30", text)
        self.assertIn("##FORMAT=<ID=DP,", text)


if __name__ == "__main__":
    unittest.main()
