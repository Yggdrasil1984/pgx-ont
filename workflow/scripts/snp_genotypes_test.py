import unittest
from tempfile import NamedTemporaryFile

import snp_genotypes as sg

ONT = [100, 50, 30, 1]

MARKERS_TSV = (
    "# comment\n"
    "gene\trsid\tchrom\tpos\tref\talts\tclinical_alt\ttranscript\thgvsc\thgvsp\thgvsg\tconsequence\tnote\tdefining_allele\tdefining_name\n"
    "F5\trs6025\tchr1\t100\tC\tA,T\tT\tNM_000130.5\tc.1601G>T,c.1601G>A\tp.Arg534Leu,p.Arg534Gln\t"
    "NC_000001.11:g.100C>T\tmissense_variant\t\tT\t\n"
    "MTHFR\trs1801133\tchr1\t200\tG\tA\tA\tNM_005957.5\tc.665C>T\tp.Ala222Val\t"
    "NC_000001.11:g.200G>A\tmissense_variant\t\tA\t\n"
    "F2\trs1799963\tchr11\t300\tG\tA\tA\tNM_000506.5\tc.*97G>A\t\tNC_000011.10:g.300G>A\t3_prime_UTR_variant\t\tA\t\n"
    "CYP3A5\trs41303343\tchr7\t500\t-\tA\tA\tNM_000777.5\tc.1035dup\tp.Thr346TyrfsTer3\t\tframeshift_variant\t\tA\t\n"
    "ALDH2\trs671\tchr12\t400\tG\tA\tA\tNM_000690.4\tc.1510G>A\tp.Glu504Lys\tNC_000012.12:g.400G>A\tmissense_variant\t\tA\t\n"
    "CYP3A4\trs2740574\tchr7\t600\tC\tT\tT\tNM_017460.5\tc.-392G>A\t\tNC_000007.14:g.600C>T\tupstream\t\tC\tCYP3A4*1B\n"
)


class TestSnpGenotypes(unittest.TestCase):

    def setUp(self):
        fh = NamedTemporaryFile("w", suffix=".tsv", delete=False)
        fh.write(MARKERS_TSV)
        fh.close()
        self.markers = sg.load_markers(fh.name)

    def _results(self, vcf_records, depths):
        rows = sg.genotype_markers(self.markers, vcf_records, depths, ONT)
        return {r["rsid"]: r for r in rows}

    def test_detected_allele_gets_its_own_hgvs(self):
        # F5: second alternative allele (T = Leiden, c.1601G>A) heterozygous
        rows = self._results({("chr1", 100): ("C", ["T"], "0|1")}, {("chr1", 100): 80})
        self.assertEqual(rows["rs6025"]["result"], sg.HETEROZYGOUS)
        self.assertEqual((rows["rs6025"]["genotype"], rows["rs6025"]["copies"]), ("C/T", 1))
        self.assertEqual(rows["rs6025"]["hgvsc"], "c.1601G>A")
        self.assertEqual(rows["rs6025"]["hgvsp"], "p.Arg534Gln")
        self.assertEqual(rows["rs6025"]["hgvsg"], "NC_000001.11:g.100C>T")

    def test_homozygous(self):
        rows = self._results({("chr1", 200): ("G", ["A"], "1/1")}, {("chr1", 200): 60})
        self.assertEqual(rows["rs1801133"]["result"], sg.HOMOZYGOUS)

    def test_no_variant_depth_decides(self):
        depths = {("chr11", 300): 30, ("chr12", 400): 29}
        rows = self._results({}, depths)
        self.assertEqual(rows["rs1799963"]["genotype"], "G/G")
        self.assertEqual(rows["rs1799963"]["copies"], 0)
        self.assertEqual(rows["rs1799963"]["colour"], "orange")
        self.assertEqual(rows["rs671"]["result"], sg.NOT_ASSESSABLE)
        # missing from the depth file counts as 0 reads
        self.assertEqual(rows["rs1801133"]["depth"], 0)
        self.assertEqual(rows["rs1801133"]["result"], sg.NOT_ASSESSABLE)

    def test_reference_call_counts_as_no_variant(self):
        rows = self._results({("chr1", 200): ("G", ["A"], "0/0")}, {("chr1", 200): 90})
        self.assertEqual(rows["rs1801133"]["genotype"], "G/G")
        self.assertEqual(rows["rs1801133"]["result"], sg.HOMOZYGOUS)
        # the clinical allele's HGVS is shown for a negative result
        self.assertEqual(rows["rs1801133"]["hgvsc"], "c.665C>T")

    def test_reference_base_defines_marker(self):
        # CYP3A4*1B: reference C is the defining allele; T/T = 0 copies of *1B
        rows = self._results({("chr7", 600): ("C", ["T"], "1/1")}, {("chr7", 600): 80})
        r = rows["rs2740574"]
        self.assertEqual((r["genotype"], r["result"], r["copies"]), ("T/T", sg.HOMOZYGOUS, 0))
        self.assertEqual(r["hgvsc"], "c.-392G>A")
        rows = self._results({}, {("chr7", 600): 80})  # reference on both: C/C = 2 copies
        self.assertEqual((rows["rs2740574"]["genotype"], rows["rs2740574"]["copies"]), ("C/C", 2))

    def test_explicit_nocall_not_assessable(self):
        rows = self._results({("chr1", 200): ("G", ["."], "./.")}, {("chr1", 200): 90})
        self.assertEqual(rows["rs1801133"]["result"], sg.NOT_ASSESSABLE)

    def test_indel_marker_not_assessable(self):
        rows = self._results({}, {("chr7", 500): 150})
        self.assertEqual(rows["rs41303343"]["result"], sg.NOT_ASSESSABLE_INDEL)
        self.assertIn("indel", rows["rs41303343"]["note"])

    def test_unexpected_allele_is_flagged(self):
        rows = self._results({("chr1", 200): ("G", ["C"], "0/1")}, {("chr1", 200): 90})
        self.assertEqual(rows["rs1801133"]["result"], sg.HETEROZYGOUS)
        self.assertEqual(rows["rs1801133"]["hgvsc"], "")
        self.assertEqual(rows["rs1801133"]["hgvsg"], "NC_000001.11:g.200G>C")
        self.assertIn("annan allel", rows["rs1801133"]["note"])


if __name__ == "__main__":
    unittest.main()
