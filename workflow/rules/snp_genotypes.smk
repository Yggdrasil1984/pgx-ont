# snp_genotypes.smk
# ──────────────────────────────────────────────────────────────────────────────
# Plain SNP genotypes with HGVS names for panel markers outside PharmCAT
# (e.g. MTHFR, F2, F5, ALDH2, CYP1A2, CYP2C8).
#
# Uses the same VCF as StarPhase/PharmCAT (starphase.bam_source) and the
# per-base depth from panel_coverage. See scripts/snp_genotypes.py for how
# "Ingen variant påvisad" / "Ej bedömbar" is decided.
#
# PREREQUISITES (config.yaml):
#   snp_markers:
#     enabled: true
#     table:   "/path/to/snp_markers.tsv"   # kept outside the repo
#   panel_coverage:
#     enabled: true                          # provides depth and thresholds
#
# OUTPUT:
#   results/snp/{sample}_{type}.snp_genotypes.tsv
# ──────────────────────────────────────────────────────────────────────────────


_snp_markers_enabled = config.get("snp_markers", {}).get("enabled", False)

if _snp_markers_enabled:

    if not config.get("panel_coverage", {}).get("enabled", False):
        raise ValueError(
            "snp_markers.enabled is true but panel_coverage.enabled is false: "
            "the SNP genotypes need the per-base depth and thresholds"
        )

    rule snp_genotypes:
        """Genotype each marker from the VCF; depth decides negative results."""
        input:
            markers = config["snp_markers"]["table"],
            vcf     = pharmcat_vcf,
            depth   = "results/coverage/{sample}_{type}.panel_depth.tsv.gz",
        output:
            tsv = "results/snp/{sample}_{type}.snp_genotypes.tsv",
        params:
            thresholds = panel_coverage_thresholds,
        log:
            "logs/snp/{sample}_{type}_snp_genotypes.log",
        script:
            "../scripts/snp_genotypes.py"
