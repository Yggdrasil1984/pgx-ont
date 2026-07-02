rule alignment_minimap2:
    input:
        fastq = lambda wildcards: units.xs(
    (wildcards.sample, wildcards.type),
    level=["sample", "type"]
)["fastq1"].iloc[0],
        ref   = config["reference"]["fasta"]
    output:
        bam = "alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai = "alignment/samtools_merge_bam/{sample}_{type}.bam.bai"
    log:
        "alignment/minimap2/{sample}_{type}.bam.log"
    threads: config.get("minimap2", {}).get("threads", 8)
    container:
        config.get("minimap2", {}).get("container", config["default_container"])
    message:
        "{rule}: align {input.fastq} to {input.ref} using minimap2"
    shell:
        "(minimap2 -ax map-ont -t {threads} {input.ref} {input.fastq} | "
        "samtools sort -o {output.bam} && "
        "samtools index {output.bam}) 2> {log}"
