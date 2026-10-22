#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read depth per panel target, classified with the report's colour scale.

Input:  panel BED (BED4, name = GENE or GENE|rsID), clinical positions BED
        (e.g. PharmCAT's positions) and per-base depth from
        `samtools depth -a -b <bed>` (chrom, 1-based pos, depth).
Output: TSV with one row per target: chrom, start, end, gene, rsid,
        min_depth, mean_depth, pct_ge_green, n_clinical_positions, basis,
        basis_depth, colour.

The colour is based on the LOWEST depth at the clinical positions inside
the target (decided by Lars-Ola 2026-10-22), so weak bases elsewhere in a
large region don't decide the colour. Targets without clinical positions
fall back to the lowest depth in the whole target. Thresholds are given
per platform in the config as the lower limits for [green, pale green,
orange, red]; 0 reads is black.
"""

import csv
import gzip

COLOURS = ("green", "palegreen", "orange", "red")
NO_READS = "black"
COLOUR_ORDER = ("black", "red", "orange", "palegreen", "green")  # worst first


def depth_colour(depth, thresholds):
    """thresholds = [green, pale green, orange, red] lower limits, e.g. [100, 50, 30, 1]."""
    for colour, limit in zip(COLOURS, thresholds):
        if depth >= limit:
            return colour
    return NO_READS


def worst_colour(colours):
    """The worst colour in a list (black < red < orange < pale green < green)."""
    return min(colours, key=COLOUR_ORDER.index) if colours else ""


def read_targets(bed_path):
    targets = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            gene, _, rsid = name.partition("|")
            targets.append({"chrom": chrom, "start": int(start), "end": int(end),
                            "gene": gene, "rsid": rsid, "depths": [],
                            "clinical": set(), "clinical_depths": {}})
    return targets


def add_clinical_positions(targets, positions_bed):
    """Mark the 1-based positions from positions_bed that lie inside each target."""
    by_chrom = {}
    for target in targets:
        by_chrom.setdefault(target["chrom"], []).append(target)
    with open(positions_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            for pos in range(int(start) + 1, int(end) + 1):
                for target in by_chrom.get(chrom, []):
                    if target["start"] < pos <= target["end"]:
                        target["clinical"].add(pos)


def add_depths(targets, depth_lines):
    """Assign each (chrom, pos, depth) line to the targets containing pos."""
    by_chrom = {}
    for target in targets:
        by_chrom.setdefault(target["chrom"], []).append(target)
    for line in depth_lines:
        chrom, pos, depth = line.split("\t")[:3]
        pos, depth = int(pos), int(depth)
        for target in by_chrom.get(chrom, []):
            if target["start"] < pos <= target["end"]:
                target["depths"].append(depth)
                if pos in target["clinical"]:
                    target["clinical_depths"][pos] = depth


def summarise(targets, thresholds):
    rows = []
    for t in targets:
        depths = t["depths"] or [0]
        min_depth = min(depths)
        if t["clinical"]:
            basis = "kliniska positioner"
            basis_depth = min(t["clinical_depths"].get(pos, 0) for pos in t["clinical"])
        else:
            basis = "hela målet"
            basis_depth = min_depth
        rows.append({
            "chrom": t["chrom"], "start": t["start"], "end": t["end"],
            "gene": t["gene"], "rsid": t["rsid"],
            "min_depth": min_depth,
            "mean_depth": round(sum(depths) / len(depths), 1),
            "pct_ge_green": round(100 * sum(d >= thresholds[0] for d in depths) / len(depths), 1),
            "n_clinical_positions": len(t["clinical"]),
            "basis": basis,
            "basis_depth": basis_depth,
            "colour": depth_colour(basis_depth, thresholds),
        })
    return rows


def write_tsv(rows, path):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["chrom"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    targets = read_targets(snakemake.input.bed)
    add_clinical_positions(targets, snakemake.input.positions)
    with gzip.open(snakemake.input.depth, "rt") as fh:
        add_depths(targets, fh)
    write_tsv(summarise(targets, snakemake.params.thresholds), snakemake.output.tsv)
