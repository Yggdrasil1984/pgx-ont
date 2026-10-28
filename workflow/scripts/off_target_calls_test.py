import json
import unittest
import zipfile
from tempfile import NamedTemporaryFile

import off_target_calls as otc
import pharmcat_allele_positions as pap

THRESHOLDS = [100, 50, 30, 1]


def _rec(gt, depth, ref="G", alt=".", info="."):
    """A VCF record as read_vcf returns it."""
    alleles = [ref] + [a for a in alt.split(",") if a != "."]
    sep = "|" if "|" in gt else "/"
    return (alleles, [None if i == "." else int(i) for i in gt.replace("|", "/").split("/")], sep, depth, info)


class TestClassify(unittest.TestCase):

    def test_variant(self):
        self.assertEqual(otc.classify(10, _rec("0|1", 40, "G", "A"), None, (40, 20), 30),
                         ("G|A", 40, otc.VARIANT, "kontrollera manuellt"))

    def test_normal_and_read_rule(self):
        self.assertEqual(otc.classify(10, _rec("0/0", 80), None, (80, 79), 30)[:3], ("G/G", 80, otc.NORMAL))
        g, d, status, note = otc.classify(10, _rec("0/0", 80, info="PGX_READRULE=0.96;PGX_READS=80"),
                                          None, (80, 77), 30)
        self.assertEqual((g, status), ("G/G", otc.NORMAL))
        self.assertIn("normal enligt läsandel (77 av 80", note)

    def test_uncertain_with_enough_reads(self):
        self.assertEqual(otc.classify(10, _rec("./.", 70), None, (70, 60), 30)[2], otc.UNCERTAIN)

    def test_low_depth_shows_what_was_found(self):
        # refcalls step wrote ./. (12 reads < 30); DeepVariant's gVCF said 0/0
        g, d, status, note = otc.classify(10, _rec("./.", 12), _rec("0/0", 12), (12, 12), 30)
        self.assertEqual((g, d, status), ("G/G", 12, otc.LOW_DEPTH))
        self.assertEqual(note, "anropat G/G; 12 av 12 läsningar normala (100 %)")
        # a low-depth variant stays a variant
        self.assertEqual(otc.classify(10, _rec("0/1", 8, "G", "A"), None, (8, 4), 30)[2], otc.VARIANT)

    def test_reads_but_no_call_and_no_reads(self):
        self.assertEqual(otc.classify(10, None, None, (80, 79), 30)[2], otc.NOT_CALLED)
        self.assertEqual(otc.classify(10, None, None, (10, 10), 30)[2], otc.LOW_DEPTH)
        self.assertEqual(otc.classify(10, None, None, (0, 0), 30), ("", 0, otc.NOT_COVERED, ""))


class TestOffTargetCalls(unittest.TestCase):

    def test_only_outside_enriched_and_gene_rules(self):
        positions = [
            {"chrom": "chr1", "pos": 150, "end": 150, "gene": "TPMT", "rsid": "rs1", "ref": "G", "alleles": "*2"},
            {"chrom": "chr1", "pos": 900, "end": 900, "gene": "TPMT", "rsid": "rs2", "ref": "G", "alleles": "*3"},
            {"chrom": "chr2", "pos": 50, "end": 50, "gene": "RYR1", "rsid": "rs3", "ref": "G", "alleles": "x"},
            {"chrom": "chr3", "pos": 50, "end": 50, "gene": "CYP2D6", "rsid": "rs4", "ref": "G", "alleles": "*4"},
        ]
        enriched = {"chr1": [(100, 200)]}
        rows = otc.off_target_calls(positions, enriched, {"TPMT", "CYP2D6"}, {}, {}, {}, THRESHOLDS, {"CYP2D6"})
        # rs1 inside, RYR1 has no target and no reads, CYP2D6 is an outside call
        self.assertEqual([(r["rsid"], r["status"], r["colour"]) for r in rows], [("rs2", otc.NOT_COVERED, "black")])

    def test_tsv_round_trip(self):
        rows = [{"gene": "G6PD", "chrom": "chrX", "pos": 5, "rsid": "rs9", "alleles": "Aveiro", "genotype": "A/A",
                 "depth": 80, "colour": "palegreen", "status": otc.NORMAL, "note": ""}]
        out = NamedTemporaryFile("w", suffix=".tsv", delete=False)
        out.close()
        otc.write_tsv(rows, out.name)
        self.assertEqual(otc.read_tsv(out.name), rows)


class TestAllelePositions(unittest.TestCase):

    def test_alleles_per_position_from_jar(self):
        definition = {
            "gene": "TPMT",
            "variants": [{"chromosome": "chr6", "position": 100, "rsid": "rs1", "ref": "T"},
                         {"chromosome": "chr6", "position": 200, "rsid": None, "ref": "CA"}],
            "namedAlleles": [{"name": "*1", "reference": True, "alleles": ["T", "CA"]},
                             {"name": "*2", "reference": False, "alleles": ["C", None]},
                             {"name": "*3", "reference": False, "alleles": ["C", "C"]}],
        }
        jar = NamedTemporaryFile(suffix=".jar", delete=False)
        jar.close()
        with zipfile.ZipFile(jar.name, "w") as z:
            z.writestr(pap.PREFIX + "TPMT_translation.json", json.dumps(definition))
        out = NamedTemporaryFile("w", suffix=".tsv", delete=False)
        out.close()
        pap.main(jar.name, out.name)
        with open(out.name) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[1:], ["chr6\t100\t100\tTPMT\trs1\tT\t*2,*3", "chr6\t200\t201\tTPMT\t\tCA\t*3"])


if __name__ == "__main__":
    unittest.main()
