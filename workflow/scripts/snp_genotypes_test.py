import unittest
from tempfile import NamedTemporaryFile

import snp_genotypes as sg

ONT = [100, 50, 30, 1]

MARKERS_TSV = (
    "# comment\n"
    "gene\trsid\tchrom\tpos\tref\talts\tclinical_alt\ttranscript\thgvsc\thgvsp\thgvsg\tconsequence\tnote\n"
    "F5\trs6025\tchr1\t100\tC\tA,T\tT\tNM_000130.5\tc.1601G>T,c.1601G>A\tp.Arg534Leu,p.Arg534Gln\t"
    "NC_000001.11:g.100C>T\tmissense_variant\t\n"
    "MTHFR\trs1801133\tchr1\t200\tG\tA\tA\tNM_005957.5\tc.665C>T\tp.Ala222Val\t"
    "NC_000001.11:g.200G>A\tmissense_variant\t\n"
    "F2\trs1799963\tchr11\t300\tG\tA\tA\tNM_000506.5\tc.*97G>A\t\tNC_000011.10:g.300G>A\t3_prime_UTR_variant\t\n"
    "CYP3A5\trs41303343\tchr7\t500\t-\tA\tA\tNM_000777.5\tc.1035dup\tp.Thr346TyrfsTer3\t\tframeshift_variant\t\n"
    "ALDH2\trs671\tchr12\t400\tG\tA\tA\tNM_000690.4\tc.1510G>A\tp.Glu504Lys\tNC_000012.12:g.400G>A\tmissense_variant\t\n"
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
        self.assertEqual(rows["rs6025"]["result"], sg.DETECTED_HET)
        self.assertEqual(rows["rs6025"]["hgvsc"], "c.1601G>A")
        self.assertEqual(rows["rs6025"]["hgvsp"], "p.Arg534Gln")
        self.assertEqual(rows["rs6025"]["hgvsg"], "NC_000001.11:g.100C>T")

    def test_homozygous(self):
        rows = self._results({("chr1", 200): ("G", ["A"], "1/1")}, {("chr1", 200): 60})
        self.assertEqual(rows["rs1801133"]["result"], sg.DETECTED_HOM)

    def test_no_variant_depth_decides(self):
        depths = {("chr11", 300): 30, ("chr12", 400): 29}
        rows = self._results({}, depths)
        self.assertEqual(rows["rs1799963"]["result"], sg.NOT_DETECTED)
        self.assertEqual(rows["rs1799963"]["colour"], "orange")
        self.assertEqual(rows["rs671"]["result"], sg.NOT_ASSESSABLE)
        # missing from the depth file counts as 0 reads
        self.assertEqual(rows["rs1801133"]["depth"], 0)
        self.assertEqual(rows["rs1801133"]["result"], sg.NOT_ASSESSABLE)

    def test_reference_call_counts_as_no_variant(self):
        rows = self._results({("chr1", 200): ("G", ["A"], "0/0")}, {("chr1", 200): 90})
        self.assertEqual(rows["rs1801133"]["result"], sg.NOT_DETECTED)
        # the clinical allele's HGVS is shown for a negative result
        self.assertEqual(rows["rs1801133"]["hgvsc"], "c.665C>T")

    def test_indel_marker_not_assessable(self):
        rows = self._results({}, {("chr7", 500): 150})
        self.assertEqual(rows["rs41303343"]["result"], sg.NOT_ASSESSABLE_INDEL)
        self.assertIn("indel", rows["rs41303343"]["note"])

    def test_unexpected_allele_is_flagged(self):
        rows = self._results({("chr1", 200): ("G", ["C"], "0/1")}, {("chr1", 200): 90})
        self.assertEqual(rows["rs1801133"]["result"], sg.DETECTED_HET)
        self.assertEqual(rows["rs1801133"]["hgvsc"], "")
        self.assertEqual(rows["rs1801133"]["hgvsg"], "NC_000001.11:g.200G>C")
        self.assertIn("annan allel", rows["rs1801133"]["note"])


if __name__ == "__main__":
    unittest.main()
