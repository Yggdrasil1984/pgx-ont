__author__ = "Lars-Ola Arvidsson"
__copyright__ = "Copyright 2026, Lars-Ola Arvidsson"
__email__ = "lars-ola.arvidsson@regionhalland.se"
__license__ = "GPL-3"


rule alignment_minimap2:
    input:
        fastq1 = lambda wildcards: units.xs(
            (wildcards.sample, wildcards.type),
            level=["sample", "type"]
        )["fastq1"].iloc[0],
        fastq2 = lambda wildcards: units.xs(
            (wildcards.sample, wildcards.type),
            level=["sample", "type"]
        )["fastq2"].iloc[0],
        ref = config["reference"]["fasta"],
        bwt = config["reference"].get("bwt", config["reference"]["fasta"] + ".bwt")
    output:
        bam = "alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai = "alignment/samtools_merge_bam/{sample}_{type}.bam.bai"
    params:
        platform = lambda wildcards: units.xs(
            (wildcards.sample, wildcards.type),
            level=["sample", "type"]
        )["platform"].iloc[0]
    log:
        "alignment/minimap2/{sample}_{type}.bam.log"
    threads: config.get("minimap2", {}).get("threads", 8)
    message:
        "{rule}: align {input.fastq1} to {input.ref} using {params.platform} mode"
    shell:
        """
        if [ "{params.platform}" = "Illumina" ]; then
            (bwa mem -t {threads} {input.ref} {input.fastq1} {input.fastq2} | \
            samtools sort -o {output.bam} && \
            samtools index {output.bam}) 2> {log}
        else
            (minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq1} | \
            samtools sort -o {output.bam} && \
            samtools index {output.bam}) 2> {log}
        fi
        """
