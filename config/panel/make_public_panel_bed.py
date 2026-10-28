#!/usr/bin/env python3
"""Build a panel target BED from public sources only.

For labs without a capture-design BED (or where the design may not be shared).
Targets:
  - every PharmCAT position (pharmcat_positions.vcf from the PharmCAT container;
    gene from INFO/PX, "POI" for positions of interest without a gene)
  - every SNP marker in the marker table
  - the whole CYP2D6 locus as StarPhase defines it (REP6 .. CYP2D7, from the
    StarPhase database), named CYP2D6

Each position gets `pad` bp on both sides; overlapping windows of the same gene
are merged (different genes are kept apart so coverage is reported per gene).
Output: BED4 (name = GENE, or POI|rsID), in .fai contig order, the format
panel_coverage and nallo's target_regions expect. For nallo's
snv_call_regions, pad the result further with config/nallo/pad_bed.py (500).

Usage:
  apptainer exec pharmcat_3.4.0.sif cat /pharmcat/pharmcat_positions.vcf > pharmcat_positions.vcf
  python3 make_public_panel_bed.py --fai GRCh38_chr.fasta.fai \\
      --pharmcat-vcf pharmcat_positions.vcf --markers snp_markers.tsv \\
      --starphase-db pbstarphase_20251002.json > public_panel_targets.bed
"""
import argparse
import json


def read_fai(path):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            lengths[f[0]] = int(f[1])
    return lengths


def read_pharmcat_vcf(path):
    """[(chrom, start, end, name)] with end = start + len(REF)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.partition("=")[::2] for kv in f[7].split(";"))
            name = info["PX"] if info.get("PX") else f"POI|{f[2]}"
            start = int(f[1]) - 1
            rows.append((f[0], start, start + len(f[3]), name))
    return rows


def read_markers(path):
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            m = dict(zip(header, f))
            pos = int(m["pos"])
            rows.append((m["chrom"], pos - 1, pos, m["gene"]))
    return rows


def read_cyp2d6_locus(path):
    """The span of StarPhase's CYP2D6 coordinates (REP6 .. CYP2D7)."""
    with open(path) as fh:
        coords = json.load(fh)["cyp2d6_config"]["cyp_coordinates"].values()
    chroms = {c["chrom"] for c in coords}
    if len(chroms) != 1:
        raise SystemExit(f"StarPhase CYP2D6 coordinates on several contigs: {chroms}")
    return [(chroms.pop(), min(c["start"] for c in coords), max(c["end"] for c in coords), "CYP2D6")]


def build(rows, pad, lengths):
    """Pad, clip, merge within each name; sorted by .fai order."""
    order = {c: i for i, c in enumerate(lengths)}
    missing = sorted({r[0] for r in rows} - set(lengths))
    if missing:
        raise SystemExit(f"contigs not in the reference .fai: {', '.join(missing)}")
    padded = sorted(((c, max(0, s - pad), min(lengths[c], e + pad), n) for c, s, e, n in rows),
                    key=lambda r: (r[3], order[r[0]], r[1]))
    merged = []
    for chrom, start, end, name in padded:
        last = merged[-1] if merged else None
        if last and last[3] == name and last[0] == chrom and start <= last[2]:
            last[2] = max(last[2], end)
        else:
            merged.append([chrom, start, end, name])
    return sorted((tuple(r) for r in merged), key=lambda r: (order[r[0]], r[1], r[2]))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--fai", required=True)
    p.add_argument("--pharmcat-vcf", required=True)
    p.add_argument("--markers", help="SNP marker table")
    p.add_argument("--starphase-db", required=True)
    p.add_argument("--pad", type=int, default=100, help="bp on each side of a position (default 100)")
    a = p.parse_args()
    rows = read_pharmcat_vcf(a.pharmcat_vcf) + read_cyp2d6_locus(a.starphase_db)
    if a.markers:
        rows += read_markers(a.markers)
    # PharmCAT's own CYP2D6 positions lie inside the StarPhase locus
    for chrom, start, end, name in build(rows, a.pad, read_fai(a.fai)):
        print(f"{chrom}\t{start}\t{end}\t{name}")


if __name__ == "__main__":
    main()
