import json
import unittest
from tempfile import NamedTemporaryFile

import pharmcat_report_data as prd


def _gene(label, phenotypes, source="MATCHER", score=None):
    return {
        "callSource": source,
        "recommendationDiplotypes": [
            {"label": label, "phenotypes": phenotypes, "activityScore": score}
        ],
    }


PHENOTYPE_JSON = {
    "geneReports": {
        "CYP2D6": _gene("*4/*5", ["Poor Metabolizer"], source="OUTSIDE", score="0.0"),
        "TPMT": _gene("*1/*3A", ["Intermediate Metabolizer"]),
        "NUDT15": _gene("*1/*1", ["Normal Metabolizer"]),
        "DPYD": _gene("Reference/Reference", ["Normal Metabolizer"], score="2.0"),
        "VKORC1": _gene("rs9923231 reference (C)/rs9923231 variant (T)", []),
        "HLA-B": _gene("Unknown/Unknown", [], source="NONE"),
        "CYP2C19": _gene("*1/*1", ["Indeterminate"]),
    }
}

GUIDELINES_TSV = (
    "# comment line\n"
    "match\tphenotype_sv\tdrugs\ttext\tclassification\tsource\treplaces_english\n"
    "TPMT=Intermediate Metabolizer;NUDT15=Normal Metabolizer\tx\ttiopuriner\tSV-TPMT-IM\tStark\ttest\tja\n"
    "DPYD=AS:2.0\tx\tfluorouracil\tSV-DPYD-2\t\ttest\t\n"
    "DPYD=AS:1.5\tx\tfluorouracil\tSV-DPYD-1.5\t\ttest\tja\n"
    "CYP2D6=Poor Metabolizer\tx\tCYP2D6\tSV-CYP2D6-PM\t\ttest\tnej\n"
)


def _annotation(gene, text, population="general", classification="Strong"):
    return {
        "drugRecommendation": text,
        "classification": classification,
        "population": population,
        "implications": [],
        "genotypes": [{"diplotypes": [{"gene": gene}]}],
    }


REPORT_JSON = {
    "drugs": {
        "CPIC Guideline Annotation": {
            "codeine": {"guidelines": [{"annotations": [_annotation("CYP2D6", "Avoid codeine")]}]},
            "clopidogrel": {"guidelines": [{"annotations": [
                _annotation("CYP2C19", "Standard dose", population="cardiology"),
                _annotation("CYP2C19", "Standard dose", population="neurology"),
            ]}]},
            "azathioprine": {"guidelines": [{"annotations": [_annotation("TPMT", "Reduce dose")]}]},
        },
        "DPWG Guideline Annotation": {
            "codeine": {"guidelines": [{"annotations": [_annotation("CYP2D6", "DPWG codeine")]}]},
        },
        "FDA Label Annotation": {
            "codeine": {"guidelines": [{"annotations": [_annotation("CYP2D6", "FDA label")]}]},
        },
    }
}


def _write(content, suffix):
    fh = NamedTemporaryFile("w", suffix=suffix, delete=False, encoding="utf-8")
    fh.write(content if isinstance(content, str) else json.dumps(content))
    fh.close()
    return fh.name


class TestPharmcatReportData(unittest.TestCase):

    def setUp(self):
        self.calls = prd.load_gene_calls(
            _write(PHENOTYPE_JSON, ".json"), not_analysed={"HLA-B"})
        self.by_gene = {c["gene"]: c for c in self.calls}
        self.guidelines = prd.load_swedish_guidelines(_write(GUIDELINES_TSV, ".tsv"))

    def test_gene_calls(self):
        cyp2d6 = self.by_gene["CYP2D6"]
        self.assertEqual(cyp2d6["source"], "StarPhase")
        self.assertEqual(cyp2d6["status"], prd.STATUS_OK)
        self.assertEqual(cyp2d6["phenotype_sv"], "Långsam metaboliserare")
        self.assertEqual(cyp2d6["activity_score"], "0.0")
        # genotype without a CPIC phenotype still counts as called
        self.assertEqual(self.by_gene["VKORC1"]["status"], prd.STATUS_OK)
        self.assertEqual(self.by_gene["HLA-B"]["status"], prd.STATUS_NOT_ANALYSED)
        self.assertEqual(self.by_gene["CYP2C19"]["status"], prd.STATUS_NOT_DETERMINED)
        self.assertEqual(self.by_gene["CYP2C19"]["diplotype"], "")

    def test_gene_filter(self):
        calls = prd.load_gene_calls(_write(PHENOTYPE_JSON, ".json"), genes={"TPMT"})
        self.assertEqual([c["gene"] for c in calls], ["TPMT"])

    def test_swedish_selection(self):
        selected = {tuple(s["genes"]): s for s in prd.select_swedish(self.calls, self.guidelines)}
        self.assertEqual(selected[("NUDT15", "TPMT")]["guideline"]["text"], "SV-TPMT-IM")
        self.assertEqual(selected[("DPYD",)]["guideline"]["text"], "SV-DPYD-2")

    def test_replaces_english(self):
        selected = prd.select_swedish(self.calls, self.guidelines)
        # TPMT+NUDT15 "ja", DPYD empty = default "ja", CYP2D6 "nej"
        self.assertEqual(prd.genes_replacing_english(selected), {"TPMT", "NUDT15", "DPYD"})
        cyp2d6 = next(s for s in selected if s["genes"] == ["CYP2D6"])
        self.assertEqual(cyp2d6["guideline"]["text"], "SV-CYP2D6-PM")

    def test_swedish_no_match_is_manual(self):
        self.by_gene["TPMT"]["phenotype"] = "Possible Intermediate Metabolizer"
        selected = {tuple(s["genes"]): s for s in prd.select_swedish(self.calls, self.guidelines)}
        self.assertEqual(selected[("NUDT15", "TPMT")]["status"], prd.STATUS_MANUAL)
        self.assertIsNone(selected[("NUDT15", "TPMT")]["guideline"])

    def test_english_recommendations(self):
        recs = prd.english_recommendations(_write(REPORT_JSON, ".json"), skip_genes={"TPMT"})
        keys = sorted((r["drug"], r["source"]) for r in recs)
        # FDA left out, TPMT covered by Swedish text, clopidogrel merged
        self.assertEqual(keys, [("clopidogrel", "CPIC"), ("codeine", "CPIC"), ("codeine", "DPWG")])
        clopidogrel = next(r for r in recs if r["drug"] == "clopidogrel")
        self.assertEqual(clopidogrel["populations"], ["cardiology", "neurology"])


if __name__ == "__main__":
    unittest.main()
