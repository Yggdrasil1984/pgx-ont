import unittest

import generate_pharmcat_report as gpr
import pharmcat_report_data as prd


def _call(gene, status=prd.STATUS_OK):
    return {"gene": gene, "diplotype": "*1/*2", "phenotype": "Normal Metabolizer",
            "phenotype_sv": "Normal metaboliserare", "activity_score": "",
            "source": "PharmCAT", "status": status}


def _target(gene, colour, depth, n_clinical):
    return {"gene": gene, "colour": colour, "basis_depth": depth,
            "n_clinical_positions": n_clinical}


class TestAddCoverage(unittest.TestCase):

    def test_clinical_targets_decide_gene_colour(self):
        calls = gpr.add_coverage([_call("TPMT")], [
            _target("TPMT", "palegreen", 70, 4),
            _target("TPMT", "red", 1, 0),  # promoter region without clinical positions
        ])
        self.assertEqual(calls[0]["coverage"], "palegreen")
        self.assertEqual(calls[0]["coverage_label"], "70×")
        self.assertFalse(calls[0]["coverage_warning"])

    def test_gene_without_clinical_targets_uses_all(self):
        calls = gpr.add_coverage([_call("MTHFR")], [
            _target("MTHFR", "green", 120, 0),
            _target("MTHFR", "red", 12, 0),
        ])
        self.assertEqual(calls[0]["coverage"], "red")
        self.assertTrue(calls[0]["coverage_warning"])

    def test_gene_not_in_panel(self):
        calls = gpr.add_coverage([_call("CFTR")], [_target("TPMT", "green", 120, 4)])
        self.assertEqual(calls[0]["status"], gpr.STATUS_NOT_IN_PANEL)
        self.assertEqual(calls[0]["diplotype"], "")

    def test_not_analysed_gene_gets_no_coverage(self):
        calls = gpr.add_coverage([_call("HLA-B", status=prd.STATUS_NOT_ANALYSED)],
                                 [_target("HLA-B", "red", 4, 0)])
        self.assertNotIn("coverage", calls[0])
        self.assertEqual(calls[0]["status"], prd.STATUS_NOT_ANALYSED)


if __name__ == "__main__":
    unittest.main()
