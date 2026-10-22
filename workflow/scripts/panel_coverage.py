#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read depth per panel target, classified with the report's colour scale.

Input:  panel BED (BED4, name = GENE or GENE|rsID) and per-base depth from
        `samtools depth -a -b <bed>` (chrom, 1-based pos, depth).
Output: TSV with one row per target: chrom, start, end, gene, rsid,
        min_depth, mean_depth, colour.

The colour is based on the LOWEST depth within the target (the weakest
base decides). Thresholds are given per platform in the config as the
lower limits for [green, pale green, orange, red]; 0 reads is black.
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
                            "gene": gene, "rsid": rsid, "depths": []})
    return targets


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


def summarise(targets, thresholds):
    rows = []
    for t in targets:
        depths = t["depths"] or [0]
        min_depth = min(depths)
        rows.append({
            "chrom": t["chrom"], "start": t["start"], "end": t["end"],
            "gene": t["gene"], "rsid": t["rsid"],
            "min_depth": min_depth,
            "mean_depth": round(sum(depths) / len(depths), 1),
            "colour": depth_colour(min_depth, thresholds),
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
    with gzip.open(snakemake.input.depth, "rt") as fh:
        add_depths(targets, fh)
    write_tsv(summarise(targets, snakemake.params.thresholds), snakemake.output.tsv)
