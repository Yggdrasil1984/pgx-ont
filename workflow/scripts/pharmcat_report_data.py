#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read PharmCAT output and select report text for the pgx report.

- load_gene_calls():          one row per gene from PharmCAT's phenotype JSON
- load_swedish_guidelines():  the Swedish text table (kept outside the repo)
- select_swedish():           pick Swedish text per gene group, or flag
                              "Manuell bedömning" when no row matches
- english_recommendations():  CPIC/DPWG text from PharmCAT's report JSON

Swedish table format (tab-separated, '#' lines are comments):
    match           GENE=<PharmCAT phenotype>[;GENE2=<phenotype>]
                    or GENE=AS:<activity score>, e.g. DPYD=AS:1.5
    phenotype_sv    Swedish phenotype description
    drugs           drugs the text applies to
    text            Swedish recommendation text (HTML, <br> for new lines)
    classification  strength of recommendation per drug
    source          where the text comes from
    replaces_english  "ja" (default) hides CPIC/DPWG text for these genes;
                      "nej" shows both (e.g. a general CYP2D6 text)
"""

import csv
import json

# Swedish names used by the interpretation group for PharmCAT phenotypes.
PHENOTYPE_SV = {
    "Normal Metabolizer": "Normal metaboliserare",
    "Intermediate Metabolizer": "Intermediär metaboliserare",
    "Poor Metabolizer": "Långsam metaboliserare",
}

SOURCE_NAMES = {"OUTSIDE": "StarPhase", "MATCHER": "PharmCAT"}

STATUS_OK = "OK"
STATUS_NOT_ANALYSED = "Ej analyserad"
STATUS_NOT_DETERMINED = "Ej bestämd"
STATUS_MANUAL = "Manuell bedömning"


def _format_score(value):
    """Activity scores as '1.5'; None/'' stay ''."""
    if value in (None, ""):
        return ""
    try:
        return f"{float(value):.1f}"
    except ValueError:
        return str(value)


def load_gene_calls(phenotype_json, genes=None, not_analysed=()):
    """
    Return one dict per gene: gene, diplotype, phenotype, phenotype_sv,
    activity_score, source, status.

    genes:        only report these genes (None = all PharmCAT genes)
    not_analysed: genes deliberately not analysed (e.g. HLA-A, HLA-B)
    """
    with open(phenotype_json) as fh:
        reports = json.load(fh)["geneReports"]

    calls = []
    for gene in sorted(reports):
        if genes is not None and gene not in genes:
            continue
        report = reports[gene]
        diplotypes = report.get("recommendationDiplotypes") or []
        first = diplotypes[0] if diplotypes else {}
        phenotypes = first.get("phenotypes") or []
        phenotype = "; ".join(phenotypes)
        source = report.get("callSource")

        if gene in not_analysed:
            status = STATUS_NOT_ANALYSED
        elif (source not in SOURCE_NAMES or "Unknown" in first.get("label", "Unknown")
              or "Indeterminate" in phenotypes):
            # Some genes (e.g. CYP4F2, VKORC1) have a genotype but no CPIC
            # phenotype by design; those still count as called.
            status = STATUS_NOT_DETERMINED
        else:
            status = STATUS_OK

        calls.append({
            "gene": gene,
            "diplotype": first.get("label", "") if status == STATUS_OK else "",
            "phenotype": phenotype if status == STATUS_OK else "",
            "phenotype_sv": PHENOTYPE_SV.get(phenotype, phenotype) if status == STATUS_OK else "",
            "activity_score": _format_score(first.get("activityScore")) if status == STATUS_OK else "",
            "source": SOURCE_NAMES.get(source, ""),
            "status": status,
        })
    return calls


def load_uncallable(match_json):
    """{gene: [alleles PharmCAT could not test]} from PharmCAT's match JSON.

    An allele is uncallable when positions that define it are missing or
    not confidently called; a "*1/*1" result then does not exclude it.
    """
    with open(match_json) as fh:
        results = json.load(fh).get("results", [])
    return {r["gene"]: sorted(r.get("uncallableHaplotypes") or []) for r in results}


def _parse_match(match):
    """'TPMT=Normal Metabolizer;NUDT15=Poor Metabolizer' -> {gene: value}."""
    conditions = {}
    for part in match.split(";"):
        gene, value = part.split("=", 1)
        conditions[gene.strip()] = value.strip()
    return conditions


def load_swedish_guidelines(path):
    """Return a list of guideline rows, each with a parsed 'conditions' dict."""
    with open(path, encoding="utf-8", newline="") as fh:
        lines = [line for line in fh if not line.startswith("#")]
    rows = []
    for row in csv.DictReader(lines, delimiter="\t"):
        row["conditions"] = _parse_match(row["match"])
        rows.append(row)
    return rows


def _condition_met(call, value):
    if call is None or call["status"] != STATUS_OK:
        return False
    if value.startswith("AS:"):
        return call["activity_score"] != "" and float(call["activity_score"]) == float(value[3:])
    return call["phenotype"] == value


def select_swedish(gene_calls, guidelines):
    """
    For each gene group in the Swedish table (e.g. TPMT+NUDT15, DPYD), return
    {genes, status, guideline}: status OK with the matching row, or
    "Manuell bedömning" with guideline None when no row matches.
    Groups whose genes are not in gene_calls are skipped.
    """
    by_gene = {call["gene"]: call for call in gene_calls}
    groups = []
    for row in guidelines:
        group = tuple(sorted(row["conditions"]))
        if group not in groups:
            groups.append(group)

    selected = []
    for group in groups:
        if not all(gene in by_gene for gene in group):
            continue
        match = next(
            (row for row in guidelines
             if tuple(sorted(row["conditions"])) == group
             and all(_condition_met(by_gene[g], v) for g, v in row["conditions"].items())),
            None,
        )
        selected.append({
            "genes": list(group),
            "status": STATUS_OK if match else STATUS_MANUAL,
            "guideline": match,
        })
    return selected


def genes_replacing_english(selected):
    """Genes whose matched Swedish row replaces the CPIC/DPWG text."""
    return {
        gene
        for s in selected
        if s["status"] == STATUS_OK
        and (s["guideline"].get("replaces_english") or "ja").strip().lower() != "nej"
        for gene in s["genes"]
    }


def english_recommendations(report_json, sources=("CPIC", "DPWG"), skip_genes=()):
    """
    Return CPIC/DPWG recommendations from PharmCAT's report JSON as dicts:
    drug, genes, source, recommendation, classification, implications,
    populations. Annotations involving only genes in skip_genes (already
    covered by Swedish text) are left out. Identical texts for the same drug
    and source (PharmCAT lists one per patient population) are merged.
    """
    with open(report_json) as fh:
        drugs_by_source = json.load(fh)["drugs"]

    recs = {}
    for source_key, drugs in drugs_by_source.items():
        source = source_key.split()[0]
        if source not in sources:
            continue
        for drug in sorted(drugs):
            for guideline in drugs[drug].get("guidelines", []):
                for annotation in guideline.get("annotations", []):
                    text = annotation.get("drugRecommendation")
                    genes = sorted({
                        diplotype["gene"]
                        for genotype in annotation.get("genotypes", [])
                        for diplotype in genotype.get("diplotypes", [])
                    })
                    if not text or not genes or set(genes) <= set(skip_genes):
                        continue
                    key = (drug, source, text)
                    if key not in recs:
                        recs[key] = {
                            "drug": drug,
                            "genes": genes,
                            "source": source,
                            "recommendation": text,
                            "classification": annotation.get("classification", ""),
                            "implications": annotation.get("implications", []),
                            "populations": [],
                        }
                    population = annotation.get("population")
                    if population and population not in recs[key]["populations"]:
                        recs[key]["populations"].append(population)
    return list(recs.values())
