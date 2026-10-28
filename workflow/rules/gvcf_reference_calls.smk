# gvcf_reference_calls.smk
# ──────────────────────────────────────────────────────────────────────────────
# Explicit reference genotypes for PharmCAT and the SNP section, from the
# DeepVariant gVCF that nallo writes per chunk.
#
# nallo does not publish the gVCF by default; publish it with a run config:
#   process { withName: '.*:CALL_SNVS:DEEPVARIANT_RUNDEEPVARIANT' {
#       publishDir = [ path: { "${params.outdir}/snvs/sample/${meta.id}/gvcf_chunks" },
#                      mode: 'copy', pattern: '*.g.vcf.gz*' ] } }
#
# Variant records (and phasing) come from nallo's phased VCF unchanged; the
# gVCF only adds reference records at PharmCAT positions and SNP-marker
# positions: 0/0 if every base is a reference call with depth >= the orange
# threshold (30 for ONT), otherwise ./. (decided by Lars-Ola 2026-10-23).
# StarPhase keeps using the phased VCF without these records.
#
# PREREQUISITES (config.yaml):
#   reference_calls:
#     enabled: true
#     positions: ["pharmcat"]   # "pharmcat" = generated from the container; SNP markers are added
#   starphase.bam_source: "nallo"
#   panel_coverage.enabled: true                        # thresholds
#
# Read rule (reference_calls.read_rule, Lars-Ola 2026-10-27): where the gVCF
# is not a confident reference call, the reads are counted from the BAM
# (samtools mpileup); >= 95% normal reads of >= 30 reads -> 0/0, marked
# INFO PGX_READRULE / PGX_READS.
#
# OUTPUT:
#   results/gvcf/{sample}_{type}.gvcf_expanded.vcf.gz   — gVCF, one record per base
#   results/gvcf/{sample}_{type}.with_refcalls.vcf.gz   — VCF for PharmCAT/SNPs
# ──────────────────────────────────────────────────────────────────────────────

import glob as _glob

_reference_calls_enabled = config.get("reference_calls", {}).get("enabled", False)
_read_rule = config.get("reference_calls", {}).get("read_rule", {})


def gvcf_chunks(wildcards):
    nallo_dir = config["starphase"]["nallo_outdir"]
    sid = nallo_sample_id(wildcards)
    chunks = sorted(_glob.glob(f"{nallo_dir}/snvs/sample/{sid}/gvcf_chunks/*.g.vcf.gz"))
    if not chunks:
        raise ValueError(
            f"No gVCF chunks for {sid} in {nallo_dir}/snvs/sample/{sid}/gvcf_chunks; "
            "run nallo with the gVCF publishDir config (see gvcf_reference_calls.smk)"
        )
    return chunks


def pharmcat_vcf(wildcards):
    """VCF for PharmCAT and the SNP section: with reference calls if enabled."""
    if _reference_calls_enabled:
        return f"results/gvcf/{wildcards.sample}_{wildcards.type}.with_refcalls.vcf.gz"
    return starphase_vcf(wildcards)


if _reference_calls_enabled:

    if config.get("starphase", {}).get("bam_source") != "nallo":
        raise ValueError("reference_calls.enabled needs starphase.bam_source: nallo (gVCF from nallo)")

    rule reference_positions:
        """Union of PharmCAT positions and SNP-marker positions as one BED."""
        input:
            beds    = [positions_bed(p) for p in config["reference_calls"]["positions"]],
            markers = config.get("snp_markers", {}).get("table", []) if config.get("snp_markers", {}).get("enabled", False) else [],
        output:
            bed = "results/gvcf/reference_positions.bed",
        shell:
            """
            ( cut -f1-3 {input.beds}
              for m in {input.markers}; do
                  grep -v '^#' "$m" | awk -F'\\t' 'NR > 1 {{
                      # SNV: the base itself; indel: also the anchor base before it
                      snv = ($5 ~ /^[ACGT]$/ && $6 ~ /^[ACGT](,[ACGT])*$/)
                      n = ($5 == "-") ? 1 : length($5)
                      if (snv) print $3 "\\t" $4 - 1 "\\t" $4
                      else     print $3 "\\t" $4 - 2 "\\t" $4 - 1 + n
                  }}'
              done ) | sort -k1,1 -k2,2n -u > {output.bed}
            """

    rule gvcf_expand:
        """Join nallo's gVCF chunks and expand reference blocks per base at the wanted positions."""
        input:
            chunks = gvcf_chunks,
            bed    = "results/gvcf/reference_positions.bed",
            ref    = config["starphase"]["reference"],
        output:
            gvcf     = temp("results/gvcf/{sample}_{type}.gvcf_sorted.vcf.gz"),
            expanded = "results/gvcf/{sample}_{type}.gvcf_expanded.vcf.gz",
        params:
            sif = config["pharmcat"]["container"],
        log:
            "logs/gvcf/{sample}_{type}_expand.log",
        shell:
            """
            (apptainer exec {params.sif} bcftools concat -a -Ou {input.chunks} \\
                | apptainer exec {params.sif} bcftools sort -Oz -o {output.gvcf} \\
             && apptainer exec {params.sif} bcftools index -t {output.gvcf} \\
             && apptainer exec {params.sif} bcftools convert --gvcf2vcf -f {input.ref} \\
                    -R {input.bed} -Oz -o {output.expanded} {output.gvcf}) > {log} 2>&1
            """

    rule reference_pileup:
        """Reads (with names) at the wanted positions, for the read rule."""
        input:
            bam = starphase_bam,
            bai = starphase_bai,
            bed = "results/gvcf/reference_positions.bed",
            ref = config["starphase"]["reference"],
        output:
            pileup = temp("results/gvcf/{sample}_{type}.reference_pileup.tsv.gz"),
        params:
            samtools = config.get("panel_coverage", {}).get("samtools", "samtools"),
            min_mapq = _read_rule.get("min_mapq", 5),
        log:
            "logs/gvcf/{sample}_{type}_reference_pileup.log",
        shell:
            """
            ({params.samtools} mpileup -B -Q 0 -q {params.min_mapq} -d 0 --output-QNAME \\
                -f {input.ref} -l {input.bed} {input.bam} | gzip > {output.pileup}) 2> {log}
            """

    rule add_reference_calls:
        """Phased variant VCF + reference records (depth rule, read rule), bgzipped and indexed."""
        input:
            vcf       = starphase_vcf,
            gvcf      = "results/gvcf/{sample}_{type}.gvcf_expanded.vcf.gz",
            positions = ["results/gvcf/reference_positions.bed"],
            pileup    = "results/gvcf/{sample}_{type}.reference_pileup.tsv.gz",
            alleles   = "results/gvcf/pharmcat_positions.vcf",   # ALT for indel reference records
        output:
            plain = temp("results/gvcf/{sample}_{type}.with_refcalls.vcf"),
        params:
            min_depth = lambda wildcards: panel_coverage_thresholds(wildcards)[2],
            read_rule = ({"min_reads": _read_rule.get("min_reads", 30),
                          "min_fraction": _read_rule.get("min_fraction", 0.95)}
                         if _read_rule.get("enabled", False) else {}),
        log:
            "logs/gvcf/{sample}_{type}_add_reference_calls.log",
        script:
            "../scripts/gvcf_reference_calls.py"

    rule bgzip_reference_calls:
        input:
            vcf = "results/gvcf/{sample}_{type}.with_refcalls.vcf",
        output:
            vcf = "results/gvcf/{sample}_{type}.with_refcalls.vcf.gz",
            tbi = "results/gvcf/{sample}_{type}.with_refcalls.vcf.gz.tbi",
        params:
            sif = config["pharmcat"]["container"],
        shell:
            """
            apptainer exec {params.sif} bgzip -c {input.vcf} > {output.vcf}
            apptainer exec {params.sif} tabix -p vcf {output.vcf}
            """
