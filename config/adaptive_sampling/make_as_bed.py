#!/usr/bin/env python3
"""Build a MinKNOW adaptive-sampling BED from the regions the pgx workflow uses.

Regions: panel targets (BED, name "GENE" or "GENE|rsID"), PharmCAT positions
(BED, same naming) and SNP markers (snp_markers table: gene, rsid, chrom, pos).
Each region gets `flank` bp on both sides (MinKNOW decides from the start of a
read, so reads that start before a target must still be accepted); overlapping
regions are merged and column 4 lists the genes. Ends are clipped to the
chromosome length from the reference .fai, which must be the reference MinKNOW
uses (chr-prefixed GRCh38).

Output is a 4-column BED (both strands targeted), in .fai contig order.

Usage:
  python3 make_as_bed.py --flank 5000 --fai GRCh38_chr.fasta.fai \\
      --bed pgx_panel_v2_targets.bed --bed pharmcat_3.4.0_positions.bed \\
      --markers snp_markers.tsv > pgx_adaptive_sampling_flank5kb.bed
"""
import argparse


def read_fai(path):
    """{chrom: length}, in file order."""
    lengths = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            lengths[f[0]] = int(f[1])
    return lengths


def read_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            gene = f[3].split("|")[0] if len(f) > 3 and f[3] else ""
            rows.append((f[0], int(f[1]), int(f[2]), gene))
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


def build(rows, flank, lengths):
    """Pad, clip and merge; [(chrom, start, end, "GENE1,GENE2")]."""
    order = {c: i for i, c in enumerate(lengths)}
    missing = sorted({r[0] for r in rows} - set(lengths))
    if missing:
        raise SystemExit(f"contigs not in the reference .fai: {', '.join(missing)}")
    padded = sorted(((c, max(0, s - flank), min(lengths[c], e + flank), g) for c, s, e, g in rows),
                    key=lambda r: (order[r[0]], r[1], r[2]))
    merged = []
    for chrom, start, end, gene in padded:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end, []])
        if gene and gene not in merged[-1][3]:
            merged[-1][3].append(gene)
    return [(c, s, e, ",".join(g) or ".") for c, s, e, g in merged]


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--flank", type=int, required=True, help="bp added on each side")
    p.add_argument("--fai", required=True, help="reference .fai (same reference as in MinKNOW)")
    p.add_argument("--bed", action="append", default=[], help="target BED; may be repeated")
    p.add_argument("--markers", action="append", default=[], help="SNP marker table; may be repeated")
    a = p.parse_args()
    rows = [r for path in a.bed for r in read_bed(path)]
    rows += [r for path in a.markers for r in read_markers(path)]
    for chrom, start, end, genes in build(rows, a.flank, read_fai(a.fai)):
        print(f"{chrom}\t{start}\t{end}\t{genes}")


if __name__ == "__main__":
    main()
