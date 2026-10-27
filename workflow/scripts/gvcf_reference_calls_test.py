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

    def test_pileup_parsing(self):
        # r1 normal with start marker, r2 mismatch, r3 normal + deletion after,
        # r4 deleted base, r5 reference skip, r6 normal with end marker
        reads = grc._pileup_reads("^].A,-2TC*>.$", ["r1", "r2", "r3", "r4", "r5", "r6"])
        self.assertEqual(reads, [("r1", True, False, False), ("r2", False, False, False),
                                 ("r3", True, False, True), ("r4", False, False, False),
                                 ("r6", True, False, False)])

    def _pileup(self, columns):
        """columns: {pos: [(read, base_token)]} on chr1 -> gzipped mpileup file."""
        lines = []
        for pos, reads in sorted(columns.items()):
            bases = "".join(tok for _, tok in reads)
            lines.append(f"chr1\t{pos}\tN\t{len(reads)}\t{bases}\t{'I' * len(reads)}\t"
                         + ",".join(r for r, _ in reads) + "\n")
        return _gz("".join(lines))

    def test_read_support_whole_span(self):
        # position chr1:30-31: 3 reads cover both bases; r2 has an insertion at 30;
        # r4 covers only 31 and is normal (not counted); r5 covers only 30 with a mismatch
        path = self._pileup({
            30: [("r1", "."), ("r2", ".+1A"), ("r3", ","), ("r5", "G")],
            31: [("r1", "."), ("r2", "."), ("r3", ","), ("r4", ".")],
        })
        support = grc.read_support(path, {("chr1", 29, 31)})
        self.assertEqual(support[("chr1", 29, 31)], (4, 2))

    def test_deletion_after_position_is_outside(self):
        # SNV position chr1:30: r2 deletes the next base (outside), r3 has an
        # insertion after it (counts), r4 deletes the next base of a 2-bp position
        path = self._pileup({
            30: [("r1", "."), ("r2", ".-1A"), ("r3", ".+1A")],
            31: [("r1", "."), ("r2", "*"), ("r3", "."), ("r4", ".-1C")],
        })
        support = grc.read_support(path, {("chr1", 29, 30), ("chr1", 30, 31)})
        self.assertEqual(support[("chr1", 29, 30)], (3, 2))
        self.assertEqual(support[("chr1", 30, 31)], (4, 3))

    def _support(self, reads, normal):
        return {("chr1", 29, 30): (reads, normal)}

    def test_read_rule(self):
        rule = {"min_reads": 30, "min_fraction": 0.95}
        out = NamedTemporaryFile("w", suffix=".vcf", delete=False)
        out.close()
        positions = {("chr1", 29, 30)}          # gVCF: ./. with 60 reads
        for (reads, normal), expected in [((40, 38), "0/0:40"),     # 95.0 %
                                          ((40, 37), "./.:60"),     # 92.5 %
                                          ((29, 29), "./.:60")]:    # too few reads
            grc.merge(self.vcf, self.bases, positions, 30, out.name, self._support(reads, normal), rule)
            with open(out.name) as fh:
                rec = [l.split("\t") for l in fh if l.startswith("chr1\t30\t")][0]
            self.assertEqual(rec[9].strip(), expected)
            if expected.startswith("0/0"):
                self.assertEqual(rec[7], "PGX_READRULE=0.950;PGX_READS=40")
        # without a rule the old behaviour stays
        grc.merge(self.vcf, self.bases, positions, 30, out.name)
        with open(out.name) as fh:
            self.assertNotIn("PGX_READRULE", fh.read())

    def test_refcall_without_genotype_is_replaced(self):
        # DeepVariant RefCall ./. at chr1:40 is decided by the read rule;
        # RefCall 0/0 (chr1:60, wanted) and RefCall away from wanted positions are kept
        vcf = _gz(HEADER + "chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT:PS\t0|1:10\n"
                           "chr1\t40\t.\tCT\tC\t0\tRefCall\t.\tGT:DP\t./.:70\n"
                           "chr1\t60\t.\tGT\tG\t0\tRefCall\t.\tGT:DP\t0/0:70\n"
                           "chr1\t90\t.\tGT\tG\t0\tRefCall\t.\tGT:DP\t./.:70\n")
        bases = grc.read_gvcf_bases(_gz(HEADER + "chr1\t40\t.\tCT\tC\t0\tRefCall\t.\tGT:DP\t./.:70\n"), {"chr1"})
        out = NamedTemporaryFile("w", suffix=".vcf", delete=False)
        out.close()
        rule = {"min_reads": 30, "min_fraction": 0.95}
        grc.merge(vcf, bases, {("chr1", 39, 41), ("chr1", 59, 60)}, 30, out.name,
                  {("chr1", 39, 41): (70, 68)}, rule)
        with open(out.name) as fh:
            rec = {(r[0], r[1]): r for r in (l.split("\t") for l in fh if not l.startswith("#"))}
        self.assertEqual(rec[("chr1", "40")][3:7], ["CT", ".", ".", "PASS"])
        self.assertEqual(rec[("chr1", "40")][9].strip(), "0/0:70")
        self.assertEqual(rec[("chr1", "60")][6], "RefCall")
        self.assertIn(("chr1", "90"), rec)   # RefCall away from wanted positions is kept

    def test_indel_position_gets_pharmcat_alleles(self):
        # PharmCAT rejects an indel reference record with ALT "."; chr1:20-21 is a
        # deletion site (CT -> C) in PharmCAT's positions file, chr2:5 a SNV
        positions_vcf = _gz(HEADER + "chr1\t20\trs1\tCT\tC\t.\t.\tPX=G\n"
                                     "chr2\t5\trs2\tC\tT\t.\t.\tPX=G\n")
        alleles = grc.read_indel_alleles(positions_vcf)
        self.assertEqual(alleles, {("chr1", 20, "CT"): "C"})
        out = NamedTemporaryFile("w", suffix=".vcf", delete=False)
        out.close()
        grc.merge(self.vcf, self.bases, {("chr1", 19, 21), ("chr2", 4, 5)}, 30, out.name, alleles=alleles)
        with open(out.name) as fh:
            rec = {(r[0], r[1]): r for r in (l.split("\t") for l in fh if not l.startswith("#"))}
        self.assertEqual(rec[("chr1", "20")][3:5], ["CT", "C"])
        self.assertEqual(rec[("chr1", "20")][9].strip(), "0/0:80")
        self.assertEqual(rec[("chr2", "5")][3:5], ["C", "."])

    def test_position_missing_from_gvcf_is_skipped(self):
        n, rec, text = self._merge({("chr1", 499, 500)})
        self.assertEqual(n, 0)
        self.assertIn("##pgx_reference_calls=<min_depth=30", text)
        self.assertIn("##FORMAT=<ID=DP,", text)


if __name__ == "__main__":
    unittest.main()
