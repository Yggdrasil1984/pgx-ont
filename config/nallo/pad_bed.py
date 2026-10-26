#!/usr/bin/env python3
"""Pad a BED file by N bp on each side and merge overlaps.

Usage: python3 pad_bed.py panel_targets.bed 500 > panel_targets_pad500.bed
"""
import sys


def main(bed_path, pad):
    rows = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            rows.append((chrom, max(0, int(start) - pad), int(end) + pad))
    merged = []
    for chrom, start, end in sorted(rows):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    for chrom, start, end in merged:
        print(f"{chrom}\t{start}\t{end}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
