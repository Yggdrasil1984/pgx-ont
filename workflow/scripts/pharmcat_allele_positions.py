#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Which named alleles each PharmCAT position defines, from the allele
definitions bundled in pharmcat.jar (*_translation.json).

Output TSV: chrom, pos, end (0-based start + len(ref), as in a BED), gene,
rsid, ref, alleles (comma-separated; alleles whose definition differs from
the reference allele at this position).

Standard library only, so it runs with the container's own python3:
  apptainer exec pharmcat.sif python3 pharmcat_allele_positions.py /pharmcat/pharmcat.jar out.tsv
"""

import json
import sys
import zipfile

PREFIX = "org/pharmgkb/pharmcat/definition/alleles/"


def allele_positions(definition):
    """[(chrom, pos, end, gene, rsid, ref, [alleles])] for one gene's definition."""
    variants = definition["variants"]
    named = definition["namedAlleles"]
    reference = next((a for a in named if a.get("reference")), named[0])
    rows = []
    for i, v in enumerate(variants):
        ref_value = reference["alleles"][i]
        alleles = [a["name"] for a in named
                   if a is not reference and a["alleles"][i] is not None and a["alleles"][i] != ref_value]
        rows.append((v["chromosome"], v["position"], v["position"] - 1 + len(v["ref"]),
                     definition["gene"], v.get("rsid") or "", v["ref"], alleles))
    return rows


def main(jar_path, out_path):
    rows = []
    with zipfile.ZipFile(jar_path) as jar:
        for name in sorted(jar.namelist()):
            if name.startswith(PREFIX) and name.endswith("_translation.json"):
                rows += allele_positions(json.loads(jar.read(name)))
    with open(out_path, "w") as out:
        out.write("chrom\tpos\tend\tgene\trsid\tref\talleles\n")
        for chrom, pos, end, gene, rsid, ref, alleles in rows:
            out.write(f"{chrom}\t{pos}\t{end}\t{gene}\t{rsid}\t{ref}\t{','.join(alleles)}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
