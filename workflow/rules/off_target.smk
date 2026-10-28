# off_target.smk
# ──────────────────────────────────────────────────────────────────────────────
# PharmCAT positions outside the enriched panel region (decided by Lars-Ola
# 2026-10-28). The panel is kept slim for MinION multiplexing, but long reads
# reach further than the baits; with nallo's calling region widened to every
# PharmCAT position (config/nallo/pad_bed.py), those positions get calls too.
#
#   - PharmCAT only gets the enriched region: panel_coverage.target_bed
#     +/- off_target.enriched_pad bp (500 = nallo's calling region decided
#     2026-10-26), so diplotypes rest on enriched positions only.
#   - Positions outside it are listed per position in their own report
#     section ("Utanför panelens mål"), never fed into a diplotype. Below the
#     depth limit the section still shows what was found, in red.
#
# PREREQUISITES: pharmcat, panel_coverage and reference_calls enabled.
#   off_target:
#     enabled: true
#     enriched_pad: 500
#
# OUTPUT:
#   results/pharmcat/enriched_regions.bed                  — targets +/- enriched_pad
#   results/pharmcat/{sample}_{type}.enriched.vcf.gz       — VCF given to PharmCAT
#   results/pharmcat/allele_positions.tsv                  — alleles per PharmCAT position
#   results/off_target/{sample}_{type}.off_target.tsv      — the report section
# ──────────────────────────────────────────────────────────────────────────────

_off_target_enabled = config.get("off_target", {}).get("enabled", False)


def pharmcat_input_vcf(wildcards):
    """VCF for PharmCAT: only the enriched region when off_target is enabled."""
    if _off_target_enabled:
        return f"results/pharmcat/{wildcards.sample}_{wildcards.type}.enriched.vcf.gz"
    return pharmcat_vcf(wildcards)


if _off_target_enabled:

    for _needed in ("pharmcat", "panel_coverage", "reference_calls"):
        if not config.get(_needed, {}).get("enabled", False):
            raise ValueError(f"off_target.enabled needs {_needed}.enabled")

    rule enriched_regions:
        """Panel targets +/- enriched_pad bp, merged."""
        input:
            bed = config["panel_coverage"]["target_bed"],
        output:
            bed = "results/pharmcat/enriched_regions.bed",
        params:
            pad    = config["off_target"].get("enriched_pad", 500),
            script = f"{workflow.basedir}/../config/nallo/pad_bed.py",
        shell:
            "python3 {params.script} {input.bed} {params.pad} > {output.bed}"

    rule pharmcat_enriched_vcf:
        """The PharmCAT input VCF restricted to the enriched region (by POS)."""
        input:
            vcf = pharmcat_vcf,
            bed = "results/pharmcat/enriched_regions.bed",
        output:
            vcf = "results/pharmcat/{sample}_{type}.enriched.vcf.gz",
        params:
            sif = config["pharmcat"]["container"],
        log:
            "logs/pharmcat/{sample}_{type}_enriched_vcf.log",
        shell:
            """
            apptainer exec {params.sif} bcftools view -T {input.bed} -Oz -o {output.vcf} {input.vcf} > {log} 2>&1
            """

    rule pharmcat_allele_positions:
        """Alleles defined at each PharmCAT position, from the container's jar."""
        output:
            tsv = "results/pharmcat/allele_positions.tsv",
        params:
            sif    = config["pharmcat"]["container"],
            script = f"{workflow.basedir}/scripts/pharmcat_allele_positions.py",
        log:
            "logs/pharmcat/allele_positions.log",
        shell:
            """
            apptainer exec {params.sif} python3 {params.script} /pharmcat/pharmcat.jar {output.tsv} > {log} 2>&1
            """

    rule off_target_calls:
        """Genotype, depth and status at PharmCAT positions outside the enriched region."""
        input:
            alleles  = "results/pharmcat/allele_positions.tsv",
            enriched = "results/pharmcat/enriched_regions.bed",
            targets  = config["panel_coverage"]["target_bed"],
            vcf      = pharmcat_vcf,
            gvcf     = "results/gvcf/{sample}_{type}.gvcf_expanded.vcf.gz",
            pileup   = "results/gvcf/{sample}_{type}.reference_pileup.tsv.gz",
        output:
            tsv = "results/off_target/{sample}_{type}.off_target.tsv",
        params:
            thresholds = panel_coverage_thresholds,
            skip_genes = config["pharmcat"].get("outside_call_genes", ["CYP2D6"]),
        log:
            "logs/off_target/{sample}_{type}.log",
        script:
            "../scripts/off_target_calls.py"
