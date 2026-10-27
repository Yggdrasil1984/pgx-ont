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
    "IFNL4\trs11322783\tchr19\t50\tTT\tG,T\tG\tENST00000606380.2\tc.65_66delinsC,c.66del\t"
    "p.Glu22AlafsTer159,p.Gly23AlafsTer158\tNC_000019.10:g.50_51delinsG\tframeshift_variant\t\tG\tIFNL4 ΔG\n"
)

# Fake reference (1-based): CYP3A5 insertion site chr7:499 T, 500 A, 501 G;
# IFNL4 site chr19:46-52 GGCCTTC (TT at 50-51)
REF = {"chr7": "N" * 496 + "CATAGG", "chr19": "N" * 45 + "GGCCTTCGA"}


def fetch(chrom, start, end):
    return REF[chrom][start - 1:end]


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

    def test_read_rule_note(self):
        rows = self._results({("chr1", 200): ("G", ["."], "0/0", "PGX_READRULE=0.971;PGX_READS=68")},
                             {("chr1", 200): 68})
        r = rows["rs1801133"]
        self.assertEqual((r["genotype"], r["result"]), ("G/G", sg.HOMOZYGOUS))
        self.assertEqual(r["note"], "normal enligt läsandel (97 % av 68 läsningar)")
        rows = self._indel({("chr7", 499): ("TA", ["."], "0/0", "PGX_READRULE=0.960;PGX_READS=50")},
                           {("chr7", 500): 80})
        self.assertEqual(rows["rs41303343"]["genotype"], "-/-")
        self.assertIn("normal enligt läsandel", rows["rs41303343"]["note"])

    def test_explicit_nocall_not_assessable(self):
        rows = self._results({("chr1", 200): ("G", ["."], "./.")}, {("chr1", 200): 90})
        self.assertEqual(rows["rs1801133"]["result"], sg.NOT_ASSESSABLE_UNCERTAIN)
        rows = self._results({("chr1", 200): ("G", ["."], "./.")}, {("chr1", 200): 12})
        self.assertEqual(rows["rs1801133"]["result"], sg.NOT_ASSESSABLE)

    def test_indel_marker_without_reference_not_assessable(self):
        rows = self._results({}, {("chr7", 500): 150})
        self.assertEqual(rows["rs41303343"]["result"], sg.NOT_ASSESSABLE_INDEL)
        self.assertIn("indel", rows["rs41303343"]["note"])
        self.assertEqual(rows["rs41303343"]["defining_allele"], "A")
        self.assertEqual(rows["rs41303343"]["hgvsc"], "c.1035dup")

    def test_normalise(self):
        # VEP insertion before chr7:500 -> VCF anchor base at 499
        self.assertEqual(sg.normalise("chr7", 500, "-", "A", fetch), (499, "T", "TA"))
        # inserting A after 500 is the same event, left-aligned
        self.assertEqual(sg.normalise("chr7", 501, "-", "A", fetch), (499, "T", "TA"))
        self.assertEqual(sg.normalise("chr19", 50, "TT", "T", fetch), (49, "CT", "C"))
        self.assertEqual(sg.normalise("chr19", 49, "CTT", "CG", fetch), (50, "TT", "G"))
        self.assertIsNone(sg.normalise("chr19", 50, "TT", "TT", fetch))

    def _indel(self, vcf_records, depths):
        rows = sg.genotype_markers(self.markers, vcf_records, depths, ONT, fetch=fetch)
        return {r["rsid"]: r for r in rows}

    def test_indel_detected_in_vcf_style(self):
        # DeepVariant writes CYP3A5*7 with the anchor base
        rows = self._indel({("chr7", 499): ("T", ["TA"], "0|1")}, {("chr7", 500): 80})
        r = rows["rs41303343"]
        self.assertEqual((r["genotype"], r["result"], r["copies"]), ("-/A", sg.HETEROZYGOUS, 1))
        self.assertEqual((r["hgvsc"], r["hgvsp"]), ("c.1035dup", "p.Thr346TyrfsTer3"))
        # IFNL4 ΔG homozygous, written with an extra anchor base
        rows = self._indel({("chr19", 49): ("CTT", ["CG"], "1/1")}, {("chr19", 50): 70})
        r = rows["rs11322783"]
        self.assertEqual((r["genotype"], r["result"], r["copies"]), ("G/G", sg.HOMOZYGOUS, 2))
        self.assertEqual(r["hgvsc"], "c.65_66delinsC")
        self.assertEqual(r["hgvsg"], "NC_000019.10:g.50_51delinsG")

    def test_indel_multiallelic_record(self):
        # one chromosome ΔG, the other the TT>T deletion (non-clinical allele)
        rows = self._indel({("chr19", 49): ("CTT", ["CG", "CT"], "1/2")}, {("chr19", 50): 70})
        r = rows["rs11322783"]
        self.assertEqual((r["genotype"], r["result"], r["copies"]), ("G/T", sg.HETEROZYGOUS, 1))

    def test_indel_reference(self):
        # gVCF reference record over anchor + site -> ref/ref
        rows = self._indel({("chr7", 499): ("TA", ["."], "0/0")}, {("chr7", 500): 80})
        r = rows["rs41303343"]
        self.assertEqual((r["genotype"], r["result"], r["copies"]), ("-/-", sg.HOMOZYGOUS, 0))
        self.assertEqual(r["hgvsc"], "c.1035dup")   # clinical allele shown for a negative result
        # no record: depth decides, as for SNVs
        rows = self._indel({}, {("chr19", 50): 30, ("chr7", 500): 29})
        self.assertEqual((rows["rs11322783"]["genotype"], rows["rs11322783"]["copies"]), ("TT/TT", 0))
        self.assertEqual(rows["rs41303343"]["result"], sg.NOT_ASSESSABLE)

    def test_indel_uncertain_reference_record(self):
        rows = self._indel({("chr7", 499): ("TA", ["."], "./.")}, {("chr7", 500): 80})
        self.assertEqual(rows["rs41303343"]["result"], sg.NOT_ASSESSABLE_UNCERTAIN)

    def test_other_variant_at_indel_marker(self):
        # a different deletion overlapping the IFNL4 site is not called as ref
        rows = self._indel({("chr19", 49): ("CTTC", ["C"], "0/1")}, {("chr19", 50): 70})
        r = rows["rs11322783"]
        self.assertEqual((r["result"], r["genotype"]), (sg.NOT_ASSESSABLE_INDEL, ""))
        self.assertIn("annan variant", r["note"])

    def test_unexpected_allele_is_flagged(self):
        rows = self._results({("chr1", 200): ("G", ["C"], "0/1")}, {("chr1", 200): 90})
        self.assertEqual(rows["rs1801133"]["result"], sg.HETEROZYGOUS)
        self.assertEqual(rows["rs1801133"]["hgvsc"], "")
        self.assertEqual(rows["rs1801133"]["hgvsg"], "NC_000001.11:g.200G>C")
        self.assertIn("annan allel", rows["rs1801133"]["note"])


if __name__ == "__main__":
    unittest.main()
