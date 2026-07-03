__author__ = "Lars-Ola Arvidsson"
__copyright__ = "Copyright 2026, Lars-Ola Arvidsson"
__email__ = "lars-ola.arvidsson@regionhalland.se"
__license__ = "GPL-3"



rule snv_indels_bed_split:
    input:
        bed = "pgx/reform_genomic_region/get_padded_bed/padded_bait_interval.bed"
    output:
        temp("snv_indels/bed_split/design_bedfile_{chr}.bed")
    log:
        "snv_indels/bed_split/bed_split_{chr}.bed.log"
    container:
        config["default_container"]
    message:
        "{rule}: generate bed file containing only {wildcards.chr}"
    shell:
        "(grep -P '^{wildcards.chr}\t' {input.bed} > {output}) 2> {log}"


rule snv_indels_haplotypecaller:
    input:
        bam   = "alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai   = "alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        fasta = config["reference"]["fasta"],
        bed   = "snv_indels/bed_split/design_bedfile_{chr}.bed"
    output:
        vcf = temp("snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf")
    params:
        model   = config["clair3"]["model_path"],
        sif     = config["clair3"]["container"],
        workdir = "snv_indels/clair3/{sample}_{type}_{chr}_workdir"
    log:
        "snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.log"
    threads: config.get("clair3", {}).get("threads", 8)
    message:
        "{rule}: call variants in {wildcards.chr} in {input.bam} using Clair3"
    shell:
        "(apptainer exec {params.sif} run_clair3.sh "
        "--bam_fn={input.bam} "
        "--ref_fn={input.fasta} "
        "--bed_fn={input.bed} "
        "--threads={threads} "
        "--platform=ont "
        "--model_path={params.model} "
        "--sample_name={wildcards.sample}_{wildcards.type} "
        "--output={params.workdir} && "
        "gunzip -c {params.workdir}/merge_output.vcf.gz > {output.vcf}) 2> {log}"


rule snv_indels_bgzip:
    input:
        "snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf"
    output:
        temp("snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.gz")
    log:
        "snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.gz.log"
    container:
        config["default_container"]
    message:
        "{rule}: bgzip {input}"
    shell:
        "(bgzip {input}) 2> {log}"


rule snv_indels_tabix:
    input:
        "snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.gz"
    output:
        temp("snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.gz.tbi")
    log:
        "snv_indels/haplotypecaller/{sample}_{type}_{chr}.vcf.gz.tbi.log"
    container:
        config["default_container"]
    message:
        "{rule}: index {input}"
    shell:
        "(tabix {input}) 2> {log}"


rule snv_indels_bcftools_concat:
    input:
        vcfs = lambda wildcards: expand(
            "snv_indels/haplotypecaller/{{sample}}_{{type}}_{chr}.vcf.gz",
            chr=extract_chr(
    "%s.fai" % (config["reference"]["fasta"]),
    filter_out=config.get("reference", {}).get("skip_chrs", [])
)
        ),
        tbis = lambda wildcards: expand(
            "snv_indels/haplotypecaller/{{sample}}_{{type}}_{chr}.vcf.gz.tbi",
            chr=extract_chr(
    "%s.fai" % (config["reference"]["fasta"]),
    filter_out=config.get("reference", {}).get("skip_chrs", [])
)
        )
    output:
        temp("snv_indels/haplotypecaller/{sample}_{type}.merged.vcf")
    log:
        "snv_indels/haplotypecaller/{sample}_{type}.merged.vcf.log"
    container:
        config["default_container"]
    message:
        "{rule}: concatenate {input.vcfs}"
    shell:
        "(bcftools concat {input.vcfs} -o {output}) 2> {log}"
