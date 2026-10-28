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
  - otherwise, read rule (decided by Lars-Ola 2026-10-27): at least
    min_reads reads and at least min_fraction (0.95) of them normal over the
    whole position (reference at every base, no insertion or deletion)
    -> GT 0/0, marked INFO PGX_READRULE=<fraction>;PGX_READS=<reads>
  - otherwise (no call, ./., or depth below min_depth)  -> GT ./.
    (decided by Lars-Ola 2026-10-23: same depth rule as the SNP section)

Depth is MIN_DP for gVCF reference blocks, DP otherwise. DeepVariant's
RefCall records with GT 0/0 are kept as before; RefCall records without a
genotype (./.) that overlap a wanted position are replaced by the record
decided here, so the read rule can apply.

A 0/0 record that starts before a wanted position but spans it (e.g. a
rejected deletion AG -> A one base earlier) does not stand for the position:
PharmCAT only reads a record starting at its own position. Such a position
still gets its own record, and the bases after the first one of a multi-base
gVCF record are not a confident reference call (the deletion reads disagree
with the reference there), so the read rule decides.

At insertion/deletion positions the record carries PharmCAT's own ALT
alleles (from pharmcat_positions.vcf, e.g. TG -> T for DPYD*3): PharmCAT's
preprocessor discards an indel reference record whose ALT is ".".

Inputs: variant VCF (gz), gVCF expanded per base (`bcftools convert
--gvcf2vcf`, gz), wanted positions as BED (0-based, end exclusive),
optionally `samtools mpileup --output-QNAME` at the wanted positions and
PharmCAT's positions VCF.
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


def read_indel_alleles(vcf_path):
    """{(chrom, pos, ref): "alt1,alt2"} for insertion/deletion positions."""
    alleles = {}
    opener = gzip.open if vcf_path.endswith((".gz", ".bgz")) else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            ref, alts = f[3], f[4]
            if any(len(a) != len(ref) for a in alts.split(",")):
                alleles[(f[0], int(f[1]), ref)] = alts
    return alleles


def read_gvcf_bases(gvcf_path, chroms):
    """{(chrom, pos): (ref_base, is_confident_ref, depth, record_pos)} per base."""
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
                bases.setdefault((f[0], int(f[1]) + offset), (base, confident_ref, depth, int(f[1])))
    return bases


def _pileup_reads(bases, names):
    """[(read name, normal base?, insertion after it?, deletion after it?)]
    from one mpileup column.

    Reference skips (> <) are left out: the read does not cover the base."""
    reads, i, k = [], 0, 0
    while i < len(bases):
        c = bases[i]
        if c == "^":
            i += 2
        elif c == "$":
            i += 1
        elif c in "+-":
            j = i + 1
            while bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j])
            if reads and reads[-1][0] == names[k - 1]:
                name, normal, ins, dele = reads[-1]
                reads[-1] = (name, normal, ins or c == "+", dele or c == "-")
        else:
            if c not in "<>":
                reads.append((names[k], c in ".,", False, False))
            k += 1
            i += 1
    return reads


def read_support(pileup_path, positions):
    """{(chrom, start, end): (reads, normal_reads)} for the wanted positions.

    A read is normal if it covers every base of the position with the
    reference base and has no insertion or deletion there. A deletion that
    starts after the last base lies outside the position and is ignored; an
    insertion after the last base counts (PharmCAT insertion sites can be
    just the anchor base). Reads that only partly cover the position count
    only if they show a difference."""
    by_base = {}
    for chrom, start, end in positions:
        for pos in range(start + 1, end + 1):
            by_base.setdefault((chrom, pos), []).append((chrom, start, end))
    seen = {p: {} for p in positions}   # position -> read -> [bases covered, abnormal]
    opener = gzip.open if pileup_path.endswith(".gz") else open
    with opener(pileup_path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]))
            if key not in by_base or len(f) < 7 or f[3] == "0":
                continue
            for read, normal, ins, dele in _pileup_reads(f[4], f[6].split(",")):
                for p in by_base[key]:
                    state = seen[p].setdefault(read, [0, False])
                    state[0] += 1
                    state[1] = state[1] or not normal or ins or (dele and key[1] < p[2])
    support = {}
    for (chrom, start, end), reads in seen.items():
        length = end - start
        counted = [abnormal for covered, abnormal in reads.values() if abnormal or covered == length]
        support[(chrom, start, end)] = (len(counted), counted.count(False))
    return support


def reference_record(chrom, start, end, bases, min_depth, support=None, read_rule=None, alleles=None):
    """VCF fields for one wanted position, or None if the gVCF lacks a base."""
    span = [bases.get((chrom, pos)) for pos in range(start + 1, end + 1)]
    if any(b is None for b in span):
        return None
    ref = "".join(b[0] for b in span)
    depth = min(b[2] for b in span)
    # a base inside a multi-base record counts only if that record is this position
    first = [b[3] in (pos, start + 1) for pos, b in zip(range(start + 1, end + 1), span)]
    confident = all(b[1] for b in span) and all(first) and depth >= min_depth
    info = "."
    if not confident and support and read_rule:
        reads, normal = support.get((chrom, start, end), (0, 0))
        if reads >= read_rule["min_reads"] and normal >= read_rule["min_fraction"] * reads:
            confident, depth = True, reads
            info = f"PGX_READRULE={normal / reads:.3f};PGX_READS={reads}"
    gt = "0/0" if confident else "./."
    alt = (alleles or {}).get((chrom, start + 1, ref), ".")
    return [chrom, str(start + 1), ".", ref, alt, ".", "PASS", info, "GT:DP", f"{gt}:{depth}"]


def merge(variant_vcf, gvcf_bases, positions, min_depth, out_path, support=None, read_rule=None,
          alleles=None):
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

    def bases_of(r):
        return {(r[0], int(r[1]) + offset) for offset in range(len(r[3]))}

    wanted = {(chrom, pos) for chrom, start, end in positions for pos in range(start + 1, end + 1)}
    # DeepVariant RefCalls without a genotype at wanted positions are decided below instead
    records = [r for r in records
               if r[6] != "RefCall" or r[9].startswith("0/0") or not bases_of(r) & wanted]
    # A record at the position itself, or a call other than 0/0 overlapping it,
    # stands; a 0/0 record that only overlaps from an earlier base does not.
    starts = {(r[0], int(r[1])) for r in records}
    covered = set()
    for r in records:
        if not r[9].startswith("0/0") and not r[9].startswith("0|0"):
            covered |= bases_of(r)
    added = []
    for chrom, start, end in sorted(positions):
        if (chrom, start + 1) in starts or any((chrom, pos) in covered for pos in range(start + 1, end + 1)):
            continue
        record = reference_record(chrom, start, end, gvcf_bases, min_depth, support, read_rule, alleles)
        if record:
            added.append(record)

    order = {c: i for i, c in enumerate(contigs)}
    all_records = sorted(records + added, key=lambda r: (order.get(r[0], len(order)), int(r[1])))
    if not any(h.startswith("##FORMAT=<ID=DP,") for h in header):
        header.insert(-1, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    header.insert(-1, f"##pgx_reference_calls=<min_depth={min_depth},source=gVCF>\n")
    if read_rule:
        header.insert(-1, '##INFO=<ID=PGX_READRULE,Number=1,Type=Float,Description="Reference by read rule: '
                          f'fraction of normal reads (>= {read_rule["min_fraction"]} of >= {read_rule["min_reads"]} reads)">\n')
        header.insert(-1, '##INFO=<ID=PGX_READS,Number=1,Type=Integer,Description="Reads counted for PGX_READRULE">\n')
    with open(out_path, "w") as out:
        out.writelines(header)
        for r in all_records:
            out.write("\t".join(r) + "\n")
    return len(added)


if __name__ == "__main__":
    positions = read_positions(snakemake.input.positions)
    chroms = {p[0] for p in positions}
    bases = read_gvcf_bases(snakemake.input.gvcf, chroms)
    read_rule = snakemake.params.read_rule or None
    support = read_support(snakemake.input.pileup, positions) if read_rule else None
    n = merge(snakemake.input.vcf, bases, positions, snakemake.params.min_depth, snakemake.output.plain,
              support, read_rule, read_indel_alleles(snakemake.input.alleles))
    with open(snakemake.log[0], "w") as log:
        log.write(f"added {n} reference records (min_depth={snakemake.params.min_depth}, read_rule={read_rule})\n")
