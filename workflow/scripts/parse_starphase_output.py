#!/usr/bin/env python3
"""
parse_starphase_output.py
─────────────────────────
Reads a pb-StarPhase JSON output file and writes a small CSV summarising the
CYP2D6 diplotype result for consumption by generate_pgx_report.py.

Usage (called by Snakemake rule starphase_cyp2d6):
    python parse_starphase_output.py \\
        --json   results/starphase/{sample}_T_starphase.json \\
        --sample HD832_T \\
        --output results/starphase/{sample}_T_cyp2d6_result.csv

Output CSV (one row, no index):
    sample,gene,diplotype,hap1,hap2,activity_score,phenotype,cnv_status,cnv_note,raw_haplotypes
    HD832_T,CYP2D6,*5/*4,*5,*4,0.0,Poor Metaboliser,deletion,Gene deletion on one allele (*5)

Activity score rules (CPIC 2024):
    Each allele gets a function value; the two are summed.
    *1  → 1.0  (normal function)
    *2  → 1.0
    *3  → 0.0  (no function)
    *4  → 0.0
    *5  → 0.0  (gene deletion)
    *6  → 0.0
    *9  → 0.5  (decreased function)
    *10 → 0.25
    *17 → 0.5
    *29 → 0.5
    *41 → 0.5
    All others not listed → None (unknown)

Phenotype from activity score (CPIC):
    0.0         → Poor Metaboliser (PM)
    0.25–0.5    → Intermediate Metaboliser (IM)
    0.75–1.0    → Intermediate Metaboliser (IM)
    1.25–2.25   → Normal Metaboliser (NM)
    ≥ 2.5       → Ultrarapid Metaboliser (UM)
"""

import argparse
import json
import csv
import sys


# ── Activity values per CPIC 2024 ────────────────────────────────────────────
ACTIVITY_VALUE = {
    "*1":   1.0,
    "*2":   1.0,
    "*3":   0.0,
    "*4":   0.0,
    "*4.001": 0.0,  # sub-allele of *4
    "*5":   0.0,   # gene deletion
    "*6":   0.0,
    "*7":   0.0,
    "*8":   0.0,
    "*9":   0.5,
    "*10":  0.25,
    "*11":  0.0,
    "*12":  0.0,
    "*14":  0.0,
    "*17":  0.5,
    "*19":  0.0,
    "*29":  0.5,
    "*31":  0.0,
    "*36":  0.0,
    "*41":  0.5,
    "*56":  0.0,
    "*59":  0.0,
}


def allele_activity(allele):
    """Return activity value for an allele string, or None if unknown."""
    # Try exact match first, then strip sub-allele suffix (e.g. *4.001 → *4)
    if allele in ACTIVITY_VALUE:
        return ACTIVITY_VALUE[allele]
    base = allele.split(".")[0]
    return ACTIVITY_VALUE.get(base, None)


def activity_to_phenotype(score):
    """Convert numeric activity score to CPIC phenotype string."""
    if score is None:
        return "Indeterminate"
    if score == 0.0:
        return "Poor Metaboliser (PM)"
    if score <= 1.0:
        return "Intermediate Metaboliser (IM)"
    if score <= 2.25:
        return "Normal Metaboliser (NM)"
    return "Ultrarapid Metaboliser (UM)"


def cnv_status(hap1, hap2):
    """
    Detect copy number / structural variant status from the diplotype.
    Returns (cnv_status_code, human-readable note).
    """
    deletion_alleles = {"*5", "*13"}
    dup_alleles = {"*1xN", "*2xN", "*35xN", "*68xN"}   # simplified — StarPhase uses specific notation

    has_deletion = hap1 in deletion_alleles or hap2 in deletion_alleles
    has_duplication = any(
        "xN" in h or "+*" in h or "DUP" in h.upper()
        for h in (hap1, hap2)
    )

    if has_deletion and has_duplication:
        return "del+dup", "Gene deletion on one allele and duplication on the other"
    if has_deletion:
        deleted = hap1 if hap1 in deletion_alleles else hap2
        return "deletion", f"Gene deletion on one allele ({deleted})"
    if has_duplication:
        return "duplication", "Gene duplication detected (increased copy number)"
    return "normal", "Normal copy number (2 copies)"


def extract_cyp2d6(starphase_json: dict, sample_name: str) -> dict:
    """
    Pull CYP2D6 result from a StarPhase output JSON.

    StarPhase JSON structure (v1.5.0):
    {
      "sample_results": {
        "<sample>": {
          "CYP2D6": {
            "simple_diplotypes": [{"hap1": "*5", "hap2": "*4", "diplotype": "*5/*4"}],
            ...
          }
        }
      }
    }
    """
    result: dict = {
        "sample":         sample_name,
        "gene":           "CYP2D6",
        "diplotype":      "NO_RESULT",
        "hap1":           ".",
        "hap2":           ".",
        "activity_score": ".",
        "phenotype":      "Indeterminate",
        "cnv_status":     "unknown",
        "cnv_note":       "No result from StarPhase",
        "raw_haplotypes": "{}",
    }

    # v1.5.0: gene_details replaces the old sample_results wrapper
    gene_details = starphase_json.get("gene_details", {})
    if gene_details:
        gene_block = gene_details.get("CYP2D6")
    else:
        # Legacy fallback for StarPhase < 1.5.0
        sample_results = starphase_json.get("sample_results", {})
        if not sample_results:
            print("WARNING: No gene_details or sample_results found in StarPhase JSON", file=sys.stderr)
            return result
        sample_block = sample_results.get(sample_name) or next(iter(sample_results.values()), None)
        if sample_block is None:
            print(f"ERROR: Sample '{sample_name}' not found in StarPhase JSON", file=sys.stderr)
            return result
        gene_block = sample_block.get("CYP2D6")

    gene_block = gene_block  # v1.5.0 path sets gene_block directly
    if gene_block is None:
        print(f"WARNING: CYP2D6 not found in StarPhase output for {sample_name}", file=sys.stderr)
        result["cnv_note"] = "CYP2D6 not present in StarPhase result"
        return result

    simple = gene_block.get("simple_diplotypes", [])
    if not simple:
        print(f"WARNING: simple_diplotypes empty for CYP2D6 in {sample_name}", file=sys.stderr)
        result["cnv_note"] = "StarPhase produced no simple_diplotypes for CYP2D6"
        return result

    best = simple[0]   # StarPhase lists the highest-confidence call first
    hap1 = best.get("hap1", ".")
    hap2 = best.get("hap2", ".")
    diplotype = best.get("diplotype", f"{hap1}/{hap2}")

    av1 = allele_activity(hap1)
    av2 = allele_activity(hap2)
    if av1 is not None and av2 is not None:
        score = av1 + av2
        score_str = str(score)
    else:
        score = None
        score_str = "."

    phenotype = activity_to_phenotype(score)
    cnv_code, cnv_note = cnv_status(hap1, hap2)

    result.update(
        diplotype      = diplotype,
        hap1           = hap1,
        hap2           = hap2,
        activity_score = score_str,
        phenotype      = phenotype,
        cnv_status     = cnv_code,
        cnv_note       = cnv_note,
        raw_haplotypes = json.dumps(gene_block, separators=(",", ":")),
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract CYP2D6 result from StarPhase JSON")
    parser.add_argument("--json",   required=True, help="StarPhase output JSON file")
    parser.add_argument("--sample", required=True, help="Sample name (e.g. HD832_T)")
    parser.add_argument("--output", required=True, help="Output CSV path")
    args = parser.parse_args()

    with open(args.json) as fh:
        data = json.load(fh)

    row = extract_cyp2d6(data, args.sample)

    fieldnames = [
        "sample", "gene", "diplotype", "hap1", "hap2",
        "activity_score", "phenotype", "cnv_status", "cnv_note",
        "raw_haplotypes",
    ]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)

    print(f"CYP2D6 result for {args.sample}: {row['diplotype']} | {row['phenotype']}")


if __name__ == "__main__":
    main()
