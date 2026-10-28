#!/usr/bin/env python3
"""Pad BED files by N bp on each side and merge them into one region list.

One or more BED/padding pairs; overlapping regions are merged.

Usage: python3 pad_bed.py panel_targets.bed 500 > panel_targets_pad500.bed
       # nallo snv_call_regions: panel +/-500 bp, plus every PharmCAT position
       # and SNP marker (public_panel_targets.bed is already +/-100 bp)
       python3 pad_bed.py panel_targets.bed 500 public_panel_targets.bed 0 > call_regions.bed
"""
import sys


def read_padded(bed_path, pad):
    rows = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            rows.append((chrom, max(0, int(start) - pad), int(end) + pad))
    return rows


def merge(rows):
    merged = []
    for chrom, start, end in sorted(rows):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    return merged


def main(args):
    if not args or len(args) % 2:
        raise SystemExit(__doc__)
    rows = []
    for bed_path, pad in zip(args[::2], args[1::2]):
        rows += read_padded(bed_path, int(pad))
    for chrom, start, end in merge(rows):
        print(f"{chrom}\t{start}\t{end}")


if __name__ == "__main__":
    main(sys.argv[1:])
