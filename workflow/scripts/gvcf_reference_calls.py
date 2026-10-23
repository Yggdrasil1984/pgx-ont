#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add explicit reference genotypes from a gVCF to a (phased) variant VCF.

Without a gVCF, a PGx position missing from the VCF can mean "normal" or
"not sequenced", so PharmCAT reports such genes as "No Result". This script
takes the variant VCF unchanged (phasing kept) and adds one reference record
per wanted position that the variant VCF lacks:

  - every base of the position (several for indel positions) is a 0/0
    reference call in the gVCF with depth >= min_depth -> GT 0/0
  - otherwise (no call, ./., or depth below min_depth)  -> GT ./.
    (decided by Lars-Ola 2026-10-23: same depth rule as the SNP section)

Depth is MIN_DP for gVCF reference blocks, DP otherwise.

Inputs: variant VCF (gz), gVCF expanded per base (`bcftools convert
--gvcf2vcf`, gz), wanted positions as BED (0-based, end exclusive).
Output: plain VCF, sorted by the variant VCF's contig order.
"""

import gzip


def read_positions(bed_paths):
    positions = set()
    for path in bed_paths:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                chrom, start, end = line.split("\t")[:3]
                positions.add((chrom, int(start), int(end)))
    return positions


def read_gvcf_bases(gvcf_path, chroms):
    """{(chrom, pos): (ref_base, is_confident_ref, depth)} for single-base records."""
    bases = {}
    with gzip.open(gvcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in chroms:
                continue
            ref, alts = f[3], f[4].split(",")
            fields = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fields.get("GT", "./.").replace("|", "/")
            depth = fields.get("MIN_DP") or fields.get("DP") or "0"
            depth = int(depth) if depth.isdigit() else 0
            confident_ref = gt == "0/0"
            for offset, base in enumerate(ref):
                bases.setdefault((f[0], int(f[1]) + offset), (base, confident_ref, depth))
    return bases


def reference_record(chrom, start, end, bases, min_depth):
    """VCF fields for one wanted position, or None if the gVCF lacks a base."""
    span = [bases.get((chrom, pos)) for pos in range(start + 1, end + 1)]
    if any(b is None for b in span):
        return None
    ref = "".join(b[0] for b in span)
    depth = min(b[2] for b in span)
    confident = all(b[1] for b in span) and depth >= min_depth
    gt = "0/0" if confident else "./."
    return [chrom, str(start + 1), ".", ref, ".", ".", "PASS", ".", "GT:DP", f"{gt}:{depth}"]


def merge(variant_vcf, gvcf_bases, positions, min_depth, out_path):
    header, records, contigs = [], [], []
    with gzip.open(variant_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                header.append(line)
                if line.startswith("##contig=<ID="):
                    contigs.append(line[len("##contig=<ID="):].split(",")[0].rstrip(">\n"))
            elif line.startswith("#"):
                header.append(line)
            else:
                records.append(line.rstrip("\n").split("\t"))

    covered = set()
    for r in records:
        for offset in range(len(r[3])):
            covered.add((r[0], int(r[1]) + offset))
    added = []
    for chrom, start, end in sorted(positions):
        if any((chrom, pos) in covered for pos in range(start + 1, end + 1)):
            continue
        record = reference_record(chrom, start, end, gvcf_bases, min_depth)
        if record:
            added.append(record)

    order = {c: i for i, c in enumerate(contigs)}
    all_records = sorted(records + added, key=lambda r: (order.get(r[0], len(order)), int(r[1])))
    if not any(h.startswith("##FORMAT=<ID=DP,") for h in header):
        header.insert(-1, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    header.insert(-1, f"##pgx_reference_calls=<min_depth={min_depth},source=gVCF>\n")
    with open(out_path, "w") as out:
        out.writelines(header)
        for r in all_records:
            out.write("\t".join(r) + "\n")
    return len(added)


if __name__ == "__main__":
    positions = read_positions(snakemake.input.positions)
    chroms = {p[0] for p in positions}
    bases = read_gvcf_bases(snakemake.input.gvcf, chroms)
    n = merge(snakemake.input.vcf, bases, positions, snakemake.params.min_depth, snakemake.output.plain)
    with open(snakemake.log[0], "w") as log:
        log.write(f"added {n} reference records (min_depth={snakemake.params.min_depth})\n")
