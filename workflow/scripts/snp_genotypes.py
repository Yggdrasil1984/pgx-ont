#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plain SNP genotypes with HGVS names for markers outside PharmCAT.

Inputs:
  markers  TSV (kept outside the repo): gene, rsid, chrom, pos, ref, alts,
           clinical_alt, transcript, hgvsc, hgvsp, hgvsg, ... where alts,
           hgvsc and hgvsp list all alternative alleles in the same order
  vcf      sample VCF (bgzipped), e.g. nallo's DeepVariant VCF
  depth    per-base depth from `samtools depth -a` (gzipped)

Rules (decided by Lars-Ola 2026-10-22, while no gVCF is available):
  - a PASS variant at the marker -> heterozygot / homozygot, HGVS of the
    allele actually found; an allele not in the marker table is reported
    with its genomic position and flagged "annan allel"
  - no PASS variant and depth >= the orange threshold (30 for ONT)
    -> "Ingen variant påvisad"
  - no PASS variant and depth below it -> "Ej bedömbar (lågt läsdjup)"
  - insertion/deletion markers (ref or alt not a single base) are not
    matched against the VCF yet (VEP and VCF place them differently):
    always "Ej bedömbar", note "indel: granska manuellt"
"""

import csv
import gzip

from panel_coverage import depth_colour

DETECTED_HET = "Heterozygot"
DETECTED_HOM = "Homozygot"
NOT_DETECTED = "Ingen variant påvisad"
NOT_ASSESSABLE = "Ej bedömbar (lågt läsdjup)"
NOT_ASSESSABLE_INDEL = "Ej bedömbar"


def is_snv(marker):
    return len(marker["ref"]) == 1 and marker["ref"] in "ACGT" and all(
        len(a) == 1 and a in "ACGT" for a in marker["alts_list"])


def load_markers(path):
    with open(path, encoding="utf-8", newline="") as fh:
        lines = [line for line in fh if not line.startswith("#")]
    markers = []
    for row in csv.DictReader(lines, delimiter="\t"):
        row["pos"] = int(row["pos"])
        for key in ("alts", "hgvsc", "hgvsp"):
            row[key + "_list"] = row[key].split(",") if row[key] else []
        markers.append(row)
    return markers


def read_vcf(vcf_path, markers):
    """PASS records at the marker positions: {(chrom, pos): (ref, alts, gt)}."""
    wanted = {(m["chrom"], m["pos"]) for m in markers}
    found = {}
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, filt = f[0], int(f[1]), f[3], f[4].split(","), f[6]
            if (chrom, pos) not in wanted or filt not in ("PASS", "."):
                continue
            fmt, sample = f[8].split(":"), f[9].split(":")
            gt = dict(zip(fmt, sample)).get("GT", "./.")
            found[(chrom, pos)] = (ref, alts, gt)
    return found


def read_depths(depth_path, markers):
    wanted = {(m["chrom"], m["pos"]) for m in markers}
    depths = {}
    with gzip.open(depth_path, "rt") as fh:
        for line in fh:
            chrom, pos, depth = line.split("\t")[:3]
            if (chrom, int(pos)) in wanted:
                depths[(chrom, int(pos))] = int(depth)
    return depths


def _hgvs_for(marker, alt):
    """(c., p., g.) names for one alternative allele of a marker."""
    if alt in marker["alts_list"]:
        i = marker["alts_list"].index(alt)
        c = marker["hgvsc_list"][i] if i < len(marker["hgvsc_list"]) else ""
        p = marker["hgvsp_list"][i] if i < len(marker["hgvsp_list"]) else ""
    else:
        c = p = ""
    g = marker["hgvsg"].rsplit(">", 1)[0] + ">" + alt if marker["hgvsg"] else ""
    return c, p, g


def genotype_markers(markers, vcf_records, depths, thresholds):
    """One result dict per marker, following the rules in the module docstring."""
    min_assessable = thresholds[2]
    results = []
    for m in markers:
        key = (m["chrom"], m["pos"])
        depth = depths.get(key, 0)
        result = {
            "gene": m["gene"], "rsid": m["rsid"], "position": f"{m['chrom']}:{m['pos']}",
            "transcript": m["transcript"], "depth": depth,
            "colour": depth_colour(depth, thresholds), "note": "",
        }
        if not is_snv(m):
            result.update(result=NOT_ASSESSABLE_INDEL, note="indel: granska manuellt",
                          hgvsc=m["hgvsc"], hgvsp=m["hgvsp"], hgvsg=m["hgvsg"])
            results.append(result)
            continue
        record = vcf_records.get(key)
        called = []
        if record:
            _, alts, gt = record
            indices = [int(a) for a in gt.replace("|", "/").split("/") if a.isdigit()]
            called = [alts[i - 1] for i in indices if i > 0]
        if called:
            alt = called[0]
            result["result"] = DETECTED_HOM if len(called) == 2 and called[0] == called[1] else DETECTED_HET
            if len(set(called)) > 1:
                result["note"] = "två olika alternativa alleler"
            if alt not in m["alts_list"]:
                result["note"] = "annan allel än i markörtabellen"
        else:
            alt = m["clinical_alt"]
            result["result"] = NOT_DETECTED if depth >= min_assessable else NOT_ASSESSABLE
        result["hgvsc"], result["hgvsp"], result["hgvsg"] = _hgvs_for(m, alt)
        results.append(result)
    return results


def write_tsv(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["gene"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_results(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


if __name__ == "__main__":
    markers = load_markers(snakemake.input.markers)
    write_tsv(
        genotype_markers(markers,
                         read_vcf(snakemake.input.vcf, markers),
                         read_depths(snakemake.input.depth, markers),
                         snakemake.params.thresholds),
        snakemake.output.tsv,
    )
