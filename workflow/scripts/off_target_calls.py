#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PharmCAT positions outside the enriched panel region ("Utanför panelens
mål"), reported separately from the PharmCAT result (decided by Lars-Ola
2026-10-28).

PharmCAT only gets positions inside the enriched region (panel targets
+/- off_target.enriched_pad bp), so its diplotypes rest on enriched
positions. Long reads still reach many positions outside it; those are
listed here per position, never fed into a diplotype:

  - a call other than 0/0                  -> "Variant", genotype in bases
  - 0/0 with >= min_depth reads (or the read rule) -> "Normal"
  - ./. with >= min_depth reads            -> "Ej bedömbar (osäkert anrop)"
  - fewer than min_depth reads: what was found is shown anyway (the caller's
    genotype if it gave one, and the share of normal reads), status
    "Lågt läsdjup" (red in the report)
  - reads but no call (outside nallo's calling region) -> "Ej anropad"
  - no reads                               -> "Ej täckt"

Genes in skip_genes (outside calls, e.g. CYP2D6 from StarPhase) are left
out. Genes without any panel target are listed only where reads exist.
"""

import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gvcf_reference_calls import read_support  # noqa: E402
from panel_coverage import depth_colour  # noqa: E402

VARIANT = "Variant"
NORMAL = "Normal"
UNCERTAIN = "Ej bedömbar (osäkert anrop)"
LOW_DEPTH = "Lågt läsdjup"
NOT_CALLED = "Ej anropad"
NOT_COVERED = "Ej täckt"

COLUMNS = ("gene", "chrom", "pos", "rsid", "alleles", "genotype", "depth", "colour", "status", "note")


def read_allele_positions(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            r["pos"], r["end"] = int(r["pos"]), int(r["end"])
            rows.append(r)
    return rows


def read_regions(bed_path, pad=0):
    """({chrom: [(start, end)]} padded, {genes with a target})."""
    regions, genes = {}, set()
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            regions.setdefault(f[0], []).append((max(0, int(f[1]) - pad), int(f[2]) + pad))
            if len(f) > 3 and f[3]:
                genes.add(f[3].split("|")[0])
    return regions, genes


def inside(regions, chrom, pos):
    return any(start < pos <= end for start, end in regions.get(chrom, ()))


def read_vcf(path):
    """{(chrom, pos): (alleles, gt, depth, info)}; gt as a list of indices (None = '.')."""
    records = {}
    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fields = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fields.get("GT", "./.")
            sep = "|" if "|" in gt else "/"
            indices = [None if i == "." else int(i) for i in gt.replace("|", "/").split("/")]
            depth = fields.get("DP") or fields.get("MIN_DP") or "0"
            alleles = [f[3]] + [a for a in f[4].split(",") if a not in (".", "<*>", "<NON_REF>")]
            records.setdefault((f[0], int(f[1])), (alleles, indices, sep, int(depth) if depth.isdigit() else 0, f[7]))
    return records


def _bases(alleles, indices, sep):
    return sep.join(alleles[i] if i is not None and i < len(alleles) else "." for i in indices)


def _share(normal, reads):
    return f"{normal} av {reads} läsningar normala ({100 * normal / reads:.0f} %)" if reads else ""


def classify(pos, record, raw, support, min_depth):
    """(genotype, depth, status, note) for one position."""
    reads, normal = support
    if record:
        alleles, indices, sep, depth, info = record
        if any(i not in (0, None) for i in indices):
            return _bases(alleles, indices, sep), depth, VARIANT, "kontrollera manuellt"
        if indices and all(i == 0 for i in indices):
            if "PGX_READRULE" in info:
                return _bases(alleles, indices, sep), reads or depth, NORMAL, "normal enligt läsandel (" + _share(normal, reads) + ")"
            if depth >= min_depth:
                return _bases(alleles, indices, sep), depth, NORMAL, ""
        if max(depth, reads) >= min_depth:
            return "", max(depth, reads), UNCERTAIN, _share(normal, reads)
    # fewer than min_depth reads, or no record at all
    if raw:
        alleles, indices, sep, depth, _ = raw
        found = _bases(alleles, indices, sep) if None not in indices else ""
        depth = max(depth, reads)
        if depth >= min_depth:
            return "", depth, UNCERTAIN, _share(normal, reads)
        if depth:
            note = ("anropat " + found + "; " if found else "inget säkert anrop; ") + _share(normal, reads)
            return found, depth, LOW_DEPTH, note.rstrip("; ")
    if reads:
        return "", reads, NOT_CALLED if reads >= min_depth else LOW_DEPTH, \
            "utanför anropsområdet; " + _share(normal, reads)
    return "", 0, NOT_COVERED, ""


def off_target_calls(positions, enriched, panel_genes, vcf, raw, support, thresholds, skip_genes=()):
    rows = []
    for p in positions:
        if p["gene"] in skip_genes or inside(enriched, p["chrom"], p["pos"]):
            continue
        key = (p["chrom"], p["pos"])
        span = (p["chrom"], p["pos"] - 1, p["end"])
        genotype, depth, status, note = classify(
            p["pos"], vcf.get(key), raw.get(key), support.get(span, (0, 0)), thresholds[2])
        if p["gene"] not in panel_genes and depth == 0:
            continue
        rows.append({
            "gene": p["gene"], "chrom": p["chrom"], "pos": p["pos"], "rsid": p["rsid"],
            "alleles": p["alleles"], "genotype": genotype, "depth": depth,
            "colour": depth_colour(depth, thresholds), "status": status, "note": note,
        })
    return sorted(rows, key=lambda r: (r["gene"], r["pos"]))


def write_tsv(rows, path):
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")


def read_tsv(path):
    """Rows as written by write_tsv (used by the report)."""
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
    for r in rows:
        r["pos"], r["depth"] = int(r["pos"]), int(r["depth"])
    return rows


if __name__ == "__main__":
    positions = read_allele_positions(snakemake.input.alleles)
    enriched, _ = read_regions(snakemake.input.enriched)
    _, panel_genes = read_regions(snakemake.input.targets)
    spans = [(p["chrom"], p["pos"] - 1, p["end"]) for p in positions]
    rows = off_target_calls(
        positions, enriched, panel_genes,
        read_vcf(snakemake.input.vcf), read_vcf(snakemake.input.gvcf),
        read_support(snakemake.input.pileup, spans),
        snakemake.params.thresholds, set(snakemake.params.skip_genes),
    )
    write_tsv(rows, snakemake.output.tsv)
