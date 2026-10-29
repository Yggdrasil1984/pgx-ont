#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plain SNP genotypes with HGVS names for markers outside PharmCAT.

Inputs:
  markers  TSV (kept outside the repo): gene, rsid, chrom, pos, ref, alts,
           clinical_alt, transcript, hgvsc, hgvsp, hgvsg, defining_allele,
           defining_name, ... where alts, hgvsc and hgvsp list all alternative
           alleles in the same order. defining_allele is the allele that
           defines the marker; it may be the reference base (CYP3A4*1B).
  vcf      sample VCF (bgzipped), e.g. nallo's DeepVariant VCF
  depth    per-base depth from `samtools depth -a` (gzipped)

Rules (decided by Lars-Ola 2026-10-22 and 2026-10-26):
  - the report gives the actual genotype (e.g. T/T), its zygosity
    (homozygot = the same allele on both chromosomes, whichever allele) and
    the number of copies of the defining allele
  - genotype from the VCF record (variant, or 0/0 from gVCF reference calls);
    without a record, depth >= the orange threshold (30 for ONT) means
    reference on both chromosomes
  - a 0/0 record set by the read rule (INFO PGX_READRULE, see
    gvcf_reference_calls.py) gets the note "normal enligt läsandel"
  - depth below the threshold -> "Ej bedömbar (lågt läsdjup)"; ./. with
    enough depth (the caller is unsure) -> "Ej bedömbar (osäkert anrop)"
  - HGVS of the non-reference allele found, else of the clinical allele; an
    allele not in the marker table is flagged "annan allel"
  - insertion/deletion markers (ref or alt not a single base; the table
    uses VEP style, "-" for an empty allele, insertion before pos) are
    matched after normalising both the marker alleles and the VCF records
    the way `bcftools norm` does (trim shared bases, add the anchor base,
    left-align against the reference). A different variant overlapping the
    marker gives "Ej bedömbar", note "annan variant vid markören"
  - with off_target enabled (Lars-Ola 2026-10-29), markers outside the
    enriched region (panel targets +/- off_target.enriched_pad) are labelled
    enriched = "nej" and reported in their own table. Their depth comes from
    the read-rule pileup (samtools depth only covers the targets); below the
    depth limit what was found is shown ("Lågt läsdjup", red), and a variant
    gets "kontrollera manuellt", as for off-target PharmCAT positions
"""

import csv
import gzip

from panel_coverage import depth_colour

HETEROZYGOUS = "Heterozygot"
HOMOZYGOUS = "Homozygot"
NOT_ASSESSABLE = "Ej bedömbar (lågt läsdjup)"
NOT_ASSESSABLE_UNCERTAIN = "Ej bedömbar (osäkert anrop)"
NOT_ASSESSABLE_INDEL = "Ej bedömbar"


class Fasta:
    """Minimal reader for an uncompressed FASTA with a .fai index."""

    def __init__(self, path):
        self.path = path
        self.index = {}
        with open(path + ".fai") as fh:
            for line in fh:
                name, length, offset, line_bases, line_width = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(line_bases), int(line_width))

    def fetch(self, chrom, start, end):
        """Bases start..end, 1-based inclusive, upper case."""
        length, offset, line_bases, line_width = self.index[chrom]
        start, end = max(start, 1), min(end, length)
        first = offset + (start - 1) // line_bases * line_width + (start - 1) % line_bases
        last = offset + (end - 1) // line_bases * line_width + (end - 1) % line_bases
        with open(self.path) as fh:
            fh.seek(first)
            return fh.read(last - first + 1).replace("\n", "").replace("\r", "").upper()


def normalise(chrom, pos, ref, alt, fetch):
    """Left-aligned, parsimonious VCF form (pos, ref, alt), or None for ref == alt.

    fetch(chrom, start, end) gives reference bases (1-based, inclusive).
    "-" stands for an empty allele (VEP style)."""
    ref, alt = ref.replace("-", "").upper(), alt.replace("-", "").upper()
    if ref == alt:
        return None
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif not ref or not alt:
            pos -= 1
            base = fetch(chrom, pos, pos)
            ref, alt = base + ref, base + alt
        else:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


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


INDEL_WINDOW = 100


def read_vcf(vcf_path, markers):
    """PASS records at SNV marker positions and within INDEL_WINDOW bp of
    indel markers: {(chrom, pos): (ref, alts, gt)}."""
    wanted = {(m["chrom"], m["pos"]) for m in markers if is_snv(m)}
    windows = [(m["chrom"], m["pos"] - INDEL_WINDOW, m["pos"] + INDEL_WINDOW)
               for m in markers if not is_snv(m)]
    found = {}
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, filt = f[0], int(f[1]), f[3], f[4].split(","), f[6]
            near = (chrom, pos) in wanted or any(
                c == chrom and lo <= pos <= hi for c, lo, hi in windows)
            if not near or filt not in ("PASS", "."):
                continue
            fmt, sample = f[8].split(":"), f[9].split(":")
            gt = dict(zip(fmt, sample)).get("GT", "./.")
            found[(chrom, pos)] = (ref, alts, gt, f[7])
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


def _genotype_bases(marker, record, depth, min_assessable):
    """[allele1, allele2] or None when not assessable."""
    if record is None:
        return [marker["ref"]] * 2 if depth >= min_assessable else None
    ref, alts, gt = record[:3]
    indices = gt.replace("|", "/").split("/")
    if not all(i.isdigit() for i in indices) or len(indices) != 2:
        return None
    alleles = [ref] + alts
    return [alleles[int(i)] for i in indices]


READ_RULE_NOTE = "normal enligt läsandel"


def _read_rule_note(record):
    """Note for a 0/0 record set by the read rule, e.g. "... (97 % av 68 läsningar)"."""
    info = dict(kv.split("=", 1) for kv in (record[3] if len(record) > 3 else ".").split(";") if "=" in kv)
    if "PGX_READRULE" not in info:
        return ""
    return f"{READ_RULE_NOTE} ({float(info['PGX_READRULE']) * 100:.0f} % av {info.get('PGX_READS', '?')} läsningar)"


def _genotype_indel(marker, vcf_records, depth, min_assessable, fetch):
    """(allele labels or None, note) for an insertion/deletion marker.

    Labels are the marker table's alleles (ref, then alts), e.g. "-/A"."""
    chrom, labels = marker["chrom"], [marker["ref"]] + marker["alts_list"]
    read_rule_used = ""
    wanted = {}
    for label in marker["alts_list"]:
        norm = normalise(chrom, marker["pos"], marker["ref"], label, fetch)
        if norm:
            wanted[norm] = label
    span_lo = min(p for p, _, _ in wanted)
    span_hi = max(p + len(r) - 1 for p, r, _ in wanted)

    covering_ref = None
    for (c, pos), record in sorted(vcf_records.items()):
        ref, alts, gt = record[:3]
        end = pos + len(ref) - 1
        if c != chrom:
            continue
        real_alts = [a for a in alts if a not in (".", "<*>", "<NON_REF>", "*")]
        if not real_alts:
            if pos <= span_hi and end >= span_lo:
                covering_ref = gt
                read_rule_used = read_rule_used or _read_rule_note(record)
            continue
        names = {}
        for i, alt in enumerate(alts, start=1):
            norm = normalise(c, pos, ref, alt, fetch) if alt in real_alts else None
            if norm in wanted:
                names[i] = wanted[norm]
        if not names:
            if pos <= span_hi and end >= span_lo:
                return None, "annan variant vid markören: granska manuellt"
            continue
        indices = gt.replace("|", "/").split("/")
        if len(indices) != 2 or not all(i.isdigit() for i in indices):
            return None, ""
        if any(int(i) > 0 and int(i) not in names for i in indices):
            return None, "annan variant vid markören: granska manuellt"
        note = _read_rule_note(record) if all(i == "0" for i in indices) else ""
        return [labels[0] if i == "0" else names[int(i)] for i in indices], note

    if covering_ref is not None:
        confident = covering_ref.replace("|", "/") == "0/0"
        return ([labels[0]] * 2 if confident else None), (read_rule_used if confident else "")
    return ([labels[0]] * 2 if depth >= min_assessable else None), ""


def _not_assessable(depth, min_assessable):
    return NOT_ASSESSABLE if depth < min_assessable else NOT_ASSESSABLE_UNCERTAIN


def genotype_markers(markers, vcf_records, depths, thresholds, fetch=None):
    """One result dict per marker, following the rules in the module docstring.

    fetch(chrom, start, end) gives reference bases; it is needed for indel
    markers, which are reported as "Ej bedömbar" without it."""
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
        defining = m.get("defining_allele") or m["clinical_alt"]
        result.update(defining_allele=defining, defining_name=m.get("defining_name", ""))
        if not is_snv(m) and fetch is not None:
            bases, note = _genotype_indel(m, vcf_records, depth, min_assessable, fetch)
            if bases is None:
                result.update(result=NOT_ASSESSABLE_INDEL if note else _not_assessable(depth, min_assessable),
                              note=note, genotype="", copies="")
                alt = m["clinical_alt"]
            else:
                result.update(
                    result=HOMOZYGOUS if bases[0] == bases[1] else HETEROZYGOUS,
                    genotype="/".join(bases),
                    copies=bases.count(defining),
                    note=note,
                )
                found = [b for b in bases if b != m["ref"]]
                alt = found[0] if found else m["clinical_alt"]
            c, p, _ = _hgvs_for(m, alt)
            g = m["hgvsg"] if alt == m["clinical_alt"] else ""
            result.update(hgvsc=c, hgvsp=p, hgvsg=g)
            results.append(result)
            continue
        if not is_snv(m):
            c, p, g = _hgvs_for(m, m["clinical_alt"])
            result.update(result=NOT_ASSESSABLE_INDEL, note="indel: granska manuellt",
                          genotype="", copies="", hgvsc=c or m["hgvsc"], hgvsp=p or m["hgvsp"],
                          hgvsg=m["hgvsg"])
            results.append(result)
            continue
        bases = _genotype_bases(m, vcf_records.get(key), depth, min_assessable)
        if bases is None:
            result.update(result=_not_assessable(depth, min_assessable), genotype="", copies="")
            alt = m["clinical_alt"]
        else:
            result.update(
                result=HOMOZYGOUS if bases[0] == bases[1] else HETEROZYGOUS,
                genotype="/".join(bases),
                copies=bases.count(defining),
            )
            found = [b for b in bases if b != m["ref"]]
            alt = found[0] if found else m["clinical_alt"]
            if any(b not in m["alts_list"] for b in found):
                result["note"] = "annan allel än i markörtabellen"
            elif key in vcf_records:
                result["note"] = _read_rule_note(vcf_records[key])
        result["hgvsc"], result["hgvsp"], result["hgvsg"] = _hgvs_for(m, alt)
        results.append(result)
    return results


LOW_DEPTH = "Lågt läsdjup"
OFF_TARGET_VARIANT_NOTE = "utanför panelens mål: kontrollera manuellt"


def _span(marker):
    return (marker["chrom"], marker["pos"] - 1, marker["pos"])


def off_target_depths(markers, depths, enriched, support):
    """Depths with off-target markers taken from the read-rule pileup.

    samtools depth (panel_depth) only covers the panel targets, so a marker
    outside them would otherwise always show 0 reads."""
    from off_target_calls import inside
    merged = dict(depths)
    for m in markers:
        key = (m["chrom"], m["pos"])
        if not inside(enriched, *key):
            merged[key] = max(depths.get(key, 0), support.get(_span(m), (0, 0))[0])
    return merged


def mark_off_target(results, markers, enriched, vcf, raw, support, thresholds):
    """Label each result enriched ("ja"/"nej"); for off-target markers apply the
    off-target rule (Lars-Ola 2026-10-28): same rules, but below the depth limit
    show what was found ("Lågt läsdjup", red) instead of "Ej bedömbar", and
    flag a variant for manual review. vcf/raw as off_target_calls.read_vcf."""
    from off_target_calls import classify, inside, VARIANT, LOW_DEPTH as OT_LOW
    min_depth = thresholds[2]
    for m, r in zip(markers, results):
        key = (m["chrom"], m["pos"])
        if inside(enriched, *key):
            r["enriched"] = "ja"
            continue
        r["enriched"] = "nej"
        if r["result"] == NOT_ASSESSABLE:
            genotype, depth, status, note = classify(
                m["pos"], vcf.get(key), raw.get(key), support.get(_span(m), (0, 0)), min_depth)
            if status in (OT_LOW, VARIANT) and depth < min_depth:
                r.update(result=LOW_DEPTH, genotype=genotype, copies="", depth=depth,
                         colour=depth_colour(depth, thresholds),
                         note=(note + "; " if note else "") + ("kontrollera manuellt" if status == VARIANT else ""))
                r["note"] = r["note"].rstrip("; ")
        elif r["genotype"] and any(b != m["ref"] for b in r["genotype"].replace("|", "/").split("/")):
            r["note"] = (r["note"] + "; " if r["note"] else "") + OFF_TARGET_VARIANT_NOTE
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
    depths = read_depths(snakemake.input.depth, markers)
    enriched_bed = getattr(snakemake.input, "enriched", None)
    if enriched_bed:
        import off_target_calls as otc
        from gvcf_reference_calls import read_support
        enriched, _ = otc.read_regions(enriched_bed)
        support = read_support(snakemake.input.pileup, [_span(m) for m in markers])
        depths = off_target_depths(markers, depths, enriched, support)
    results = genotype_markers(markers,
                               read_vcf(snakemake.input.vcf, markers),
                               depths,
                               snakemake.params.thresholds,
                               fetch=Fasta(snakemake.input.ref).fetch)
    if enriched_bed:
        results = mark_off_target(results, markers, enriched,
                                  otc.read_vcf(snakemake.input.vcf), otc.read_vcf(snakemake.input.gvcf),
                                  support, snakemake.params.thresholds)
    write_tsv(results, snakemake.output.tsv)
