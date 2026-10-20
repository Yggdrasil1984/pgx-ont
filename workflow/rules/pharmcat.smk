# pharmcat.smk
# ──────────────────────────────────────────────────────────────────────────────
# Snakemake rules for PharmCAT (MPL-2.0): star alleles, phenotypes and
# CPIC/DPWG recommendations for all PharmCAT genes.
#
# CYP2D6 (and HLA, when present) is NOT called by PharmCAT itself. It comes
# from StarPhase as an "outside call" file (--pharmcat-tsv, see
# starphase_cyp2d6.smk). Other genes are called from the same VCF that
# StarPhase uses (starphase.bam_source decides which VCF that is).
#
# PREREQUISITES (config.yaml):
#   pharmcat:
#     enabled:   true
#     container: "/home/your-username/singularity_images/pharmcat_3.4.0.sif"
#     outside_call_genes: ["CYP2D6"]   # genes taken from StarPhase instead
#                                      # of PharmCAT's own VCF matcher
#   starphase:
#     enabled:   true        # PharmCAT needs the StarPhase outside-call file
#
# StarPhase's --pharmcat-tsv holds calls for ALL genes it types (~35, incl.
# HLA). An outside call REPLACES PharmCAT's own call for that gene, so only
# the genes in pharmcat.outside_call_genes are passed on. Adding a gene to
# that list makes StarPhase the caller of record for it (clinical decision).
#
# The container bundles its own GRCh38 reference and PGx position list, so no
# downloads happen at runtime.
#
# NOTE: the preprocessor is run WITHOUT -0 / --absent-to-ref /
# --unspecified-to-ref. Positions missing from the VCF stay "not called"
# instead of being assumed normal; turning them into "normal" is a clinical
# decision that belongs with the coverage QC.
#
# OUTPUT:
#   results/pharmcat/{sample}_{type}.preprocessed.vcf.bgz — PharmCAT-ready VCF
#   results/pharmcat/{sample}_{type}.missing_pgx_var.vcf  — PGx positions absent from the VCF
#   results/pharmcat/{sample}_{type}.outside_calls.tsv    — StarPhase calls passed to PharmCAT
#   results/pharmcat/{sample}_{type}.phenotype.json       — calls + phenotypes per gene
#   results/pharmcat/{sample}_{type}.report.json          — recommendations (read by the pgx report)
#   results/pharmcat/{sample}_{type}.report.html          — PharmCAT's own report (reference only)
# ──────────────────────────────────────────────────────────────────────────────


_pharmcat_enabled = config.get("pharmcat", {}).get("enabled", False)

if _pharmcat_enabled and not config.get("starphase", {}).get("enabled", False):
    raise ValueError(
        "pharmcat.enabled is true but starphase.enabled is false: PharmCAT "
        "needs the StarPhase outside-call file for CYP2D6"
    )


if _pharmcat_enabled:

    # ── Rule 1: prepare the VCF for PharmCAT ─────────────────────────────────
    rule pharmcat_preprocess:
        """
        Normalise the sample VCF to PharmCAT's expected positions and
        representation, and list PGx positions that are missing.
        """
        input:
            vcf = starphase_vcf,
        output:
            vcf     = "results/pharmcat/{sample}_{type}.preprocessed.vcf.bgz",
            missing = "results/pharmcat/{sample}_{type}.missing_pgx_var.vcf",
        params:
            outdir = "results/pharmcat",
            sif    = config["pharmcat"]["container"],
        log:
            "logs/pharmcat/{sample}_{type}_preprocess.log",
        shell:
            """
            apptainer exec {params.sif} pharmcat_vcf_preprocessor \\
                -vcf {input.vcf} \\
                -o   {params.outdir} \\
                -bf  {wildcards.sample}_{wildcards.type} \\
                > {log} 2>&1
            """


    # ── Rule 2: keep only the chosen genes from the StarPhase TSV ────────────
    rule pharmcat_outside_calls:
        """
        Filter the StarPhase PharmCAT TSV to pharmcat.outside_call_genes
        (default CYP2D6), keeping the header line. Fails if a listed gene
        is missing, so a gene is never silently left to the VCF matcher.
        """
        input:
            tsv = "results/starphase/{sample}_{type}_pharmcat.tsv",
        output:
            tsv = "results/pharmcat/{sample}_{type}.outside_calls.tsv",
        params:
            genes = ",".join(config.get("pharmcat", {}).get("outside_call_genes", ["CYP2D6"])),
        log:
            "logs/pharmcat/{sample}_{type}_outside_calls.log",
        shell:
            """
            awk -F'\\t' -v genes="{params.genes}" '
                BEGIN {{ n = split(genes, g, ","); for (i = 1; i <= n; i++) want[g[i]] = 1 }}
                /^#/ {{ print; next }}
                ($1 in want) {{ print; found[$1] = 1 }}
                END {{ for (k in want) if (!(k in found)) {{ print "missing outside call for " k > "/dev/stderr"; exit 1 }} }}
            ' {input.tsv} > {output.tsv} 2> {log}
            """


    # ── Rule 3: run PharmCAT with StarPhase as outside call ──────────────────
    rule pharmcat_run:
        """
        Run PharmCAT matcher, phenotyper and reporter. CYP2D6 is taken from
        the StarPhase outside-call file (-po).
        """
        input:
            vcf     = "results/pharmcat/{sample}_{type}.preprocessed.vcf.bgz",
            outside = "results/pharmcat/{sample}_{type}.outside_calls.tsv",
        output:
            phenotype = "results/pharmcat/{sample}_{type}.phenotype.json",
            report    = "results/pharmcat/{sample}_{type}.report.json",
            html      = "results/pharmcat/{sample}_{type}.report.html",
        params:
            outdir = "results/pharmcat",
            sif    = config["pharmcat"]["container"],
        log:
            "logs/pharmcat/{sample}_{type}_run.log",
        shell:
            """
            apptainer exec {params.sif} java -jar /pharmcat/pharmcat.jar \\
                -vcf {input.vcf} \\
                -po  {input.outside} \\
                -o   {params.outdir} \\
                -bf  {wildcards.sample}_{wildcards.type} \\
                -reporterJson -reporterHtml \\
                > {log} 2>&1
            """
