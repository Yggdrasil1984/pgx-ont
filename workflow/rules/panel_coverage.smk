# panel_coverage.smk
# ──────────────────────────────────────────────────────────────────────────────
# Read depth per panel target with the report's colour scale.
#
# Uses the same BAM as StarPhase (starphase.bam_source decides which one).
# The colour is set by the lowest depth inside each target; thresholds are
# chosen by the sample's platform in units.tsv.
#
# PREREQUISITES (config.yaml):
#   panel_coverage:
#     enabled:    true
#     target_bed: "/path/to/panel_targets.bed"   # BED4, name = GENE or GENE|rsID
#     samtools:   "/usr/bin/samtools"
#     min_mapq:   0          # 0 = count all reads (as mosdepth's default)
#     thresholds:            # lower limits for [green, pale green, orange, red]
#       ONT: [100, 50, 30, 1]
#
# OUTPUT:
#   results/coverage/{sample}_{type}.panel_depth.tsv.gz — per-base depth
#   results/coverage/{sample}_{type}.panel_coverage.tsv — one row per target
# ──────────────────────────────────────────────────────────────────────────────


_panel_coverage_enabled = config.get("panel_coverage", {}).get("enabled", False)


def panel_coverage_thresholds(wildcards):
    """Thresholds for the sample's platform; fails if none are configured."""
    platform = units[
        (units["sample"] == wildcards.sample) & (units["type"] == wildcards.type)
    ]["platform"].iloc[0]
    thresholds = config["panel_coverage"]["thresholds"]
    if platform not in thresholds:
        raise ValueError(
            f"panel_coverage.thresholds has no entry for platform {platform!r} "
            f"(sample {wildcards.sample}_{wildcards.type})"
        )
    return thresholds[platform]


if _panel_coverage_enabled:

    rule panel_depth:
        """Per-base depth over the panel targets, including 0-depth positions."""
        input:
            bam = starphase_bam,
            bai = starphase_bai,
            bed = config["panel_coverage"]["target_bed"],
        output:
            depth = "results/coverage/{sample}_{type}.panel_depth.tsv.gz",
        params:
            samtools = config["panel_coverage"].get("samtools", "samtools"),
            min_mapq = config["panel_coverage"].get("min_mapq", 0),
        log:
            "logs/coverage/{sample}_{type}_panel_depth.log",
        shell:
            """
            ({params.samtools} depth -a -Q {params.min_mapq} -b {input.bed} {input.bam} \\
                | gzip > {output.depth}) 2> {log}
            """


    rule panel_coverage:
        """Summarise depth per target and assign the colour."""
        input:
            bed   = config["panel_coverage"]["target_bed"],
            depth = "results/coverage/{sample}_{type}.panel_depth.tsv.gz",
        output:
            tsv = "results/coverage/{sample}_{type}.panel_coverage.tsv",
        params:
            thresholds = panel_coverage_thresholds,
        log:
            "logs/coverage/{sample}_{type}_panel_coverage.log",
        script:
            "../scripts/panel_coverage.py"
