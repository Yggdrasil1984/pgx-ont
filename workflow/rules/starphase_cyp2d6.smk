# starphase_cyp2d6.smk
# ──────────────────────────────────────────────────────────────────────────────
# Snakemake rules for CYP2D6 paralog-aware diplotyping using pb-StarPhase.
#
# These rules sit between the alignment/variant-calling layer (nallo or Clair3)
# and the clinical report generation step.
#
# PREREQUISITES (must exist in config.yaml):
#   starphase:
#     container:  "/home/your-username/singularity_images/starphase.sif"
#     database:   "/home/your-username/reference/starphase_db/pbstarphase_20251002.json"
#     reference:  "/home/your-username/reference/GRCh38_chr.fasta"  # full GRCh38, chr-prefixed
#
# INPUT FILES REQUIRED PER SAMPLE:
#   - A haplotagged BAM: produced by nallo (preferred) or WhatsHap haplotagging
#     of the Clair3/minimap2 BAM.  See config key starphase.bam_pattern below.
#   - A phased SNV VCF: produced by nallo DeepVariant+WhatsHap, or by the
#     standard pgx Clair3 pipeline (unphased is also accepted by StarPhase).
#
# OUTPUT:
#   results/starphase/{sample}_{type}_starphase.json    — full StarPhase JSON
#   results/starphase/{sample}_{type}_cyp2d6_result.csv — parsed CYP2D6 result
#     (columns: sample,gene,diplotype,hap1,hap2,activity_score,phenotype,
#                cnv_status,cnv_note)
# ──────────────────────────────────────────────────────────────────────────────


# ── Helper: resolve BAM path ──────────────────────────────────────────────────
# The BAM source depends on whether nallo was run (preferred haplotagged BAM)
# or the pipeline uses its own Clair3 BAM.  Configure in config.yaml:
#
#   starphase:
#     bam_source: "nallo"        # use nallo haplotagged BAM
#     nallo_outdir: "/path/to/nallo_results"
#   OR
#     bam_source: "pgx"          # use pgx alignment BAM (from ont_alignment.smk)
#   OR
#     bam_source: "custom"       # explicit files per sample, e.g. validation
#     custom_inputs:             # samples (HG00276, HG001) from elsewhere
#       HG00276_T:               # key = {sample}_{type}
#         bam: "/path/to/HG00276.haplotagged.bam"   # needs .bam.bai next to it
#         vcf: "/path/to/HG00276.filtered.vcf.gz"
#
# Default: pgx (no nallo required — works with current Chapter 2/3 setup).


# ── Enabled guard ────────────────────────────────────────────────────────────
_starphase_enabled = config.get("starphase", {}).get("enabled", False)


def starphase_custom_input(wildcards, key):
    """Return the bam/vcf path given under starphase.custom_inputs.{sample}_{type}."""
    sample_type = f"{wildcards.sample}_{wildcards.type}"
    custom = config["starphase"].get("custom_inputs", {})
    if sample_type not in custom:
        raise ValueError(
            f"starphase.bam_source is 'custom' but starphase.custom_inputs has no entry for {sample_type}"
        )
    return custom[sample_type][key]


def starphase_bam(wildcards):
    """Return the BAM path for StarPhase based on bam_source config."""
    source = config.get("starphase", {}).get("bam_source", "pgx")
    if source == "custom":
        return starphase_custom_input(wildcards, "bam")
    if source == "nallo":
        nallo_dir = config["starphase"]["nallo_outdir"]
        return (
            f"{nallo_dir}/aligned_reads/sample/"
            f"{wildcards.sample}_{wildcards.type}/"
            f"{wildcards.sample}_{wildcards.type}_haplotagged.bam"
        )
    # Default: pgx alignment BAM (from ont_alignment.smk)
    return f"alignment/samtools_merge_bam/{wildcards.sample}_{wildcards.type}.bam"


def starphase_bai(wildcards):
    return starphase_bam(wildcards) + ".bai"


def starphase_vcf(wildcards):
    """Return the VCF path for StarPhase."""
    source = config.get("starphase", {}).get("bam_source", "pgx")
    if source == "custom":
        return starphase_custom_input(wildcards, "vcf")
    if source == "nallo":
        nallo_dir = config["starphase"]["nallo_outdir"]
        return (
            f"{nallo_dir}/snvs/sample/"
            f"{wildcards.sample}_{wildcards.type}/"
            f"{wildcards.sample}_{wildcards.type}_deepvariant_snvs.vcf.gz"
        )
    # Default: pgx Clair3 merged VCF
    return (
        f"snv_indels/haplotypecaller/"
        f"{wildcards.sample}_{wildcards.type}.merged.vcf.gz"
    )


# ── Rule 1: run StarPhase diplotype ──────────────────────────────────────────
if _starphase_enabled:

    rule starphase_diplotype:
        """
        Run pb-StarPhase to call CYP2D6 (and other CPIC Level A genes) from a
        haplotagged BAM + phased VCF.

        StarPhase needs:
          --bam        haplotagged BAM (WhatsHap or nallo LongPhase)
          --vcf        phased SNV VCF (gzipped, tabix-indexed)
          --reference  GRCh38 chromosome FASTA (chr22 minimum for CYP2D6)
          --database   pre-built StarPhase JSON database
          --output-calls  JSON file to write the calls to

        StarPhase v1.5.0 writes the JSON directly to --output-calls (there is
        no --output-dir or --threads option in v1.5.0).
        """
        input:
            bam  = starphase_bam,
            bai  = starphase_bai,
            vcf  = starphase_vcf,
            ref  = config["starphase"]["reference"],
            db   = config["starphase"]["database"],
        output:
            json = "results/starphase/{sample}_{type}_starphase.json",
        params:
            sif      = config["starphase"]["container"],
        log:
            "logs/starphase/{sample}_{type}_diplotype.log",
        shell:
            """
            apptainer exec {params.sif} pbstarphase diplotype \\
                --bam          {input.bam}  \\
                --vcf          {input.vcf}  \\
                --reference    {input.ref}  \\
                --database     {input.db}   \\
                --output-calls {output.json} \\
                > {log} 2>&1
            """


    # ── Rule 2: parse StarPhase JSON → CYP2D6 CSV ────────────────────────────────
    rule starphase_parse_cyp2d6:
        """
        Extract CYP2D6 diplotype, phenotype, activity score, and CNV status
        from the StarPhase JSON into a one-row CSV for the report generator.
        """
        input:
            json = "results/starphase/{sample}_{type}_starphase.json",
        output:
            csv  = "results/starphase/{sample}_{type}_cyp2d6_result.csv",
        params:
            sample_name = "{sample}_{type}",
            script      = workflow.basedir + "/scripts/parse_starphase_output.py",
        log:
            "logs/starphase/{sample}_{type}_parse.log",
        shell:
            """
            python {params.script} \\
                --json   {input.json} \\
                --sample {params.sample_name} \\
                --output {output.csv} \\
                > {log} 2>&1
            """
