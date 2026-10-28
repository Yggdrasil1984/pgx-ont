#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Swedish pgx report built from PharmCAT (all genes), StarPhase (CYP2D6
details), the Swedish guideline table and per-target panel coverage.

Rules for the gene summary:
- a PharmCAT gene with no target in the panel is "Ej analyserad"
  (PharmCAT reports something for every gene, but the panel never saw it)
- a gene's colour comes from its targets with clinical positions (all
  targets only if it has none); a called gene whose worst such target is
  red or black is flagged
  "Kontrollera läsdjup"; the genotype itself is never changed here
"""

import csv
import json
import os
import sys

import jinja2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pharmcat_report_data as prd  # noqa: E402
from panel_coverage import worst_colour  # noqa: E402
from generate_pgx_report import load_cyp2d6_result  # noqa: E402
from snp_genotypes import read_results as read_snp_results  # noqa: E402

STATUS_NOT_IN_PANEL = "Ej analyserad (ingår ej i panelen)"
WARN_COLOURS = ("red", "black")


def read_coverage(path):
    if not path:
        return []
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        for key in ("start", "end", "min_depth", "basis_depth", "n_clinical_positions"):
            row[key] = int(row[key])
    return rows


def add_coverage(gene_calls, coverage):
    """Add worst coverage colour per gene; genes without targets are not analysed."""
    if not coverage:
        return gene_calls
    by_gene = {}
    for row in coverage:
        by_gene.setdefault(row["gene"], []).append(row)
    for call in gene_calls:
        rows = by_gene.get(call["gene"])
        if call["status"] == prd.STATUS_NOT_ANALYSED:
            continue
        if not rows:
            call.update(status=STATUS_NOT_IN_PANEL, diplotype="", phenotype="",
                        phenotype_sv="", activity_score="", source="")
            continue
        # Targets with clinical positions decide the gene's colour; targets
        # without (e.g. promoter regions) only count if the gene has no others.
        clinical = [r for r in rows if r["n_clinical_positions"] > 0]
        rows = clinical or rows
        worst = worst_colour([r["colour"] for r in rows])
        worst_depth = min(r["basis_depth"] for r in rows if r["colour"] == worst)
        call["coverage"] = worst
        call["coverage_label"] = f"{worst_depth}×"
        call["coverage_warning"] = worst in WARN_COLOURS
    return gene_calls


def build_report(phenotype_json, report_json, match_json, swedish_tsv, starphase_csv, starphase_json,
                 coverage_tsv, template, sample, platform, thresholds,
                 not_analysed=("HLA-A", "HLA-B"), sources=("CPIC", "DPWG"), snp_tsv=None, vcf=None):
    gene_calls = prd.load_gene_calls(phenotype_json, not_analysed=set(not_analysed))
    coverage = read_coverage(coverage_tsv)
    gene_calls = add_coverage(gene_calls, coverage)
    uncallable = prd.load_uncallable(match_json)
    read_rule = prd.load_read_rule_calls(vcf, match_json) if vcf else {}
    for call in gene_calls:
        ok = call["status"] == prd.STATUS_OK
        call["uncallable"] = uncallable.get(call["gene"], []) if ok else []
        call["read_rule"] = read_rule.get(call["gene"], []) if ok else []

    guidelines = prd.load_swedish_guidelines(swedish_tsv)
    called = [c for c in gene_calls if c["status"] == prd.STATUS_OK]
    swedish = prd.select_swedish(called, guidelines)
    by_gene = {c["gene"]: c for c in gene_calls}
    for s in swedish:
        s["phenotypes"] = ", ".join(f"{g} {by_gene[g]['phenotype']}" for g in s["genes"])
    covered = prd.genes_replacing_english(swedish)
    reported = {c["gene"] for c in called}
    english = [r for r in prd.english_recommendations(report_json, sources=sources, skip_genes=covered)
               if set(r["genes"]) & reported]

    with open(report_json) as fh:
        meta = json.load(fh)
    with open(starphase_json) as fh:
        starphase_version = json.load(fh).get("pbstarphase_version", "")

    with open(template, encoding="utf-8") as fh:
        html = jinja2.Template(fh.read()).render(
            sample=sample,
            platform=platform,
            pharmcat_version=meta.get("pharmcatVersion", ""),
            pharmcat_data_version=meta.get("dataVersion", ""),
            starphase_version=starphase_version,
            thresholds=thresholds,
            genes=gene_calls,
            swedish=swedish,
            english=english,
            cyp2d6=load_cyp2d6_result(starphase_csv),
            coverage=coverage,
            snps=read_snp_results(snp_tsv) if snp_tsv else [],
        )
    return html


if __name__ == "__main__":
    html = build_report(
        phenotype_json=snakemake.input.phenotype,
        report_json=snakemake.input.report,
        match_json=snakemake.input.match,
        swedish_tsv=snakemake.input.swedish,
        starphase_csv=snakemake.input.starphase_csv,
        starphase_json=snakemake.input.starphase_json,
        coverage_tsv=getattr(snakemake.input, "coverage", None),
        template=snakemake.input.template,
        sample=f"{snakemake.wildcards.sample}_{snakemake.wildcards.type}",
        platform=snakemake.params.platform,
        thresholds=snakemake.params.thresholds,
        sources=snakemake.params.sources,
        snp_tsv=getattr(snakemake.input, "snps", None),
        vcf=getattr(snakemake.input, "vcf", None),
    )
    with open(snakemake.output.html, "w", encoding="utf-8") as fh:
        fh.write(html)
