# pharmcat_positions.smk
# ──────────────────────────────────────────────────────────────────────────────
# PharmCAT's own positions, taken from the container in use
# (/pharmcat/pharmcat_positions.vcf), so they always match the PharmCAT version.
#
# In the config, the value "pharmcat" for panel_coverage.clinical_positions or
# for an entry of reference_calls.positions means this generated BED:
#   chrom, pos-1, pos-1+len(REF), GENE|rsID   (GENE = "POI" for positions of
#   interest without a gene, e.g. rs12777823)
#
# OUTPUT:
#   results/gvcf/pharmcat_positions.vcf
#   results/gvcf/pharmcat_positions.bed
# ──────────────────────────────────────────────────────────────────────────────

PHARMCAT_POSITIONS_BED = "results/gvcf/pharmcat_positions.bed"


def positions_bed(path):
    """A configured positions BED, or the generated one for "pharmcat"."""
    return PHARMCAT_POSITIONS_BED if path == "pharmcat" else path


if config.get("pharmcat", {}).get("container"):

    rule pharmcat_positions_vcf:
        """PharmCAT's own positions/alleles, taken from the container in use."""
        output:
            vcf = "results/gvcf/pharmcat_positions.vcf",
        params:
            sif = config["pharmcat"]["container"],
        shell:
            "apptainer exec {params.sif} cat /pharmcat/pharmcat_positions.vcf > {output.vcf}"

    rule pharmcat_positions_bed:
        """PharmCAT positions as BED4 (name GENE|rsID), REF length for indels."""
        input:
            vcf = "results/gvcf/pharmcat_positions.vcf",
        output:
            bed = PHARMCAT_POSITIONS_BED,
        shell:
            """
            awk 'BEGIN {{ OFS = "\\t" }} !/^#/ {{
                gene = "POI"
                if (match($8, /PX=[^;]+/)) gene = substr($8, RSTART + 3, RLENGTH - 3)
                print $1, $2 - 1, $2 - 1 + length($4), gene "|" $3
            }}' {input.vcf} > {output.bed}
            """
