__author__ = "Lars-Ola Arvidsson"
__copyright__ = "Copyright 2026, Lars-Ola Arvidsson"
__email__ = "lars-ola.arvidsson@regionhalland.se"
__license__ = "GPL-3"


rule alignment_minimap2:
    input:
        fastq = lambda wildcards: units.xs(
            (wildcards.sample, wildcards.type),
            level=["sample", "type"]
        ).query("platform == 'ONT'")["fastq1"].iloc[0],
        ref   = config["reference"]["fasta"],
        bwt   = config["reference"].get("bwt", config["reference"]["fasta"] + ".bwt")
    output:
        bam = "alignment/samtools_merge_bam/{sample}_{type}_ONT.bam",
        bai = "alignment/samtools_merge_bam/{sample}_{type}_ONT.bam.bai"
    log:
        "alignment/minimap2/{sample}_{type}_ONT.bam.log"
    threads: config.get("minimap2", {}).get("threads", 8)
    message:
        "{rule}: align {input.fastq} to {input.ref} using ONT mode"
    shell:
        "(minimap2 -ax map-ont -t {threads} "
        "-R '@RG\\tID:{wildcards.sample}_{wildcards.type}\\tSM:{wildcards.sample}_{wildcards.type}\\tPL:ONT\\tLB:{wildcards.sample}_{wildcards.type}' "
        "{input.ref} {input.fastq} | "
        "/usr/bin/samtools sort -o {output.bam} && "
        "/usr/bin/samtools index {output.bam}) 2> {log}"


rule alignment_bwa:
    input:
        fastq1 = lambda wildcards: units.xs(
            (wildcards.sample, wildcards.type),
            level=["sample", "type"]
        ).query("platform == 'Illumina'")["fastq1"].iloc[0],
        fastq2 = lambda wildcards: units.xs(
            (wildcards.sample, wildcards.type),
            level=["sample", "type"]
        ).query("platform == 'Illumina'")["fastq2"].iloc[0],
        ref = config["reference"]["fasta"],
        bwt = config["reference"].get("bwt", config["reference"]["fasta"] + ".bwt")
    output:
        bam = "alignment/samtools_merge_bam/{sample}_{type}_Illumina.bam",
        bai = "alignment/samtools_merge_bam/{sample}_{type}_Illumina.bam.bai"
    log:
        "alignment/bwa/{sample}_{type}_Illumina.bam.log"
    threads: config.get("minimap2", {}).get("threads", 8)
    message:
        "{rule}: align {input.fastq1} to {input.ref} using Illumina mode"
    shell:
        "(bwa mem -t {threads} "
        "-R '@RG\\tID:{wildcards.sample}_{wildcards.type}\\tSM:{wildcards.sample}_{wildcards.type}\\tPL:ILLUMINA\\tLB:{wildcards.sample}_{wildcards.type}' "
        "{input.ref} {input.fastq1} {input.fastq2} | "
        "/usr/bin/samtools sort -o {output.bam} && "
        "/usr/bin/samtools index {output.bam}) 2> {log}"


rule alignment_merge_platform:
    input:
        bam = lambda wildcards: (
            "alignment/samtools_merge_bam/{sample}_{type}_ONT.bam"
            if units.xs(
                (wildcards.sample, wildcards.type),
                level=["sample", "type"]
            )["platform"].iloc[0] == "ONT"
            else
            "alignment/samtools_merge_bam/{sample}_{type}_Illumina.bam"
        ),
        bai = lambda wildcards: (
            "alignment/samtools_merge_bam/{sample}_{type}_ONT.bam.bai"
            if units.xs(
                (wildcards.sample, wildcards.type),
                level=["sample", "type"]
            )["platform"].iloc[0] == "ONT"
            else
            "alignment/samtools_merge_bam/{sample}_{type}_Illumina.bam.bai"
        )
    output:
        bam = "alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai = "alignment/samtools_merge_bam/{sample}_{type}.bam.bai"
    message:
        "{rule}: link platform BAM to unified path for {wildcards.sample}_{wildcards.type}"
    shell:
        "ln -sf $(basename {input.bam}) {output.bam} && "
        "ln -sf $(basename {input.bai}) {output.bai}"