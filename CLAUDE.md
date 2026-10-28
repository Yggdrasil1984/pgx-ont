# CLAUDE.md — pgx ONT adaptation (Region Halland)

This file gives Claude Code the context for this project. Read it at the start of every session.

## Who I'm working with

- Lars-Ola, molecular biologist at Region Halland. No formal programming background.
- Explain technical concepts in plain language. Molecular biology analogies are welcome.
- Lars-Ola makes the decisions on direction, tool choice, and clinical interpretation. Claude does the hands-on terminal work.

## How to work

- **Go step by step and verify each step.** After each meaningful change, show what changed and confirm it worked before moving on.
- **Dry run before real run.** Always run `snakemake -n` (or the Nextflow equivalent, e.g. `-stub` / `-preview`) before a real execution.
- **Ask before anything hard to undo.** This includes deleting files, overwriting results, `git push`, `git reset --hard`, changing conda environments, or installing or upgrading system packages.
- **Stop and ask when a decision is clinical or architectural.** For example: which CYP2D6 caller to use, changing filtering thresholds, or anything that affects a reported genotype.
- **Keep commits small and well described**, one logical change per commit.
- When something fails, diagnose it and explain the cause in a sentence or two before fixing it.

## Data rules (important)

- **Never read, copy, or process patient data.** Work only with public reference samples (HG002, HG00276, GIAB, 1000 Genomes) and the test data listed below.
- If a path looks like it may contain patient samples, stop and ask.

## Licensing constraint (hard requirement)

The clinical pipeline must be redistributable through Genomic Medicine Sweden (GMS) repositories. That means:

- Only MIT, GPL-3, Apache, BSD-style, or MPL-2.0 licensed tools may be used in the clinical path. (MPL-2.0 approved by Lars-Ola on 2026-10-19 for PharmCAT.)
- **Do not use `epi2me-labs/wf-pgx` or Chinook in the clinical pipeline.** They are under the ONT Public License, which allows research use only and is not compatible with GMS redistribution. They may be used only for research comparison and validation.
- Check the license before suggesting any new tool.

## Project goal

Build an ONT long-read pharmacogenomics pipeline for clinical use, simple and cheap enough for smaller regional labs to run. Region Halland will use **ONT only**, so ONT compatibility is the priority. The pipeline should still **keep PacBio HiFi possible**: the user states the platform per sample, and that choice selects the upstream settings (nallo `--preset`: `ONT_R10`, `ONT_R10_AS` for adaptive sampling, `revio`/`pacbio` for HiFi). Everything after nallo should be platform-independent.

Planned clinical method (2027): a long-read-adapted Qiagen hybrid capture panel covering 49 biomarkers. Most are SNP regions; CYP2D6 has double tiling across the whole region.

## Target architecture

```
ONT (or HiFi) FASTQ/BAM
  → nallo     (minimap2 + DeepVariant + WhatsHap/LongPhase; haplotagged BAM, phased SNV VCF, Sniffles2 SV VCF)
  → StarPhase (CYP2D6 diplotyping incl. deletions, duplications, CYP2D6/CYP2D7 hybrids; also HLA) → --pharmcat-tsv
  → PharmCAT  (other genes from the VCF + StarPhase CYP2D6 as "outside call"; CPIC/DPWG phenotypes and recommendations)
  → pgx       (Swedish clinical HTML report, filled from PharmCAT JSON, plus pgx's own coverage QC)
```

A hybrid of nallo and pgx: nallo replaces pgx's own alignment/variant calling; pgx adds coverage QC over the panel targets (`panel_coverage`), gVCF reference calls, plain SNP genotypes and the Swedish report. pgx's original steps (depth_of_baits, missing positions, old report) are still there for the non-nallo path. In nallo, Paraphase and the rare-disease steps (VEP/CADD/genmod ranking, methylation, STRs, assembly) are not needed.

Licenses: nallo MIT, StarPhase BSD-3-Clause-Clear, PharmCAT MPL-2.0, pgx GPL-3.

## Repositories and locations

| What | Where | Notes |
|---|---|---|
| pgx pipeline | `~/projects/pgx` | Snakemake 7.18.0, conda env `pgx` (Python 3.8). Work on branch **`nanopore-adaptation`** |
| nallo | `~/projects/nallo` | Nextflow (26.04.6), Java 21, `dev` branch. Test run succeeded |
| DPYD amplicon pipeline | `~/projects/dpyd_amplicon/` | Standalone and separate from pgx. Validated on 6 blind samples |
| pgx test data | `~/projects/pgx/.tests/integration/` | |
| GRCh38 reference | `~/reference/GRCh38.fasta.gz` | |
| chr22 reference | `~/reference/chr22.fa` | Converted from Ensembl to `chr22` naming |
| GRCh38 with `chr` names | `~/reference/GRCh38_chr.fasta` | Required by StarPhase (full genome, chr-prefixed) |
| StarPhase database | `~/reference/starphase_db/pbstarphase_20251002.json` | Pre-built database |
| StarPhase container | `~/singularity_images/starphase.sif` | pbstarphase 1.5.0 |
| PharmCAT container | `~/singularity_images/pharmcat_3.4.0.sif` | From `docker://pgkb/pharmcat:3.4.0`; bundles its own reference, no downloads at runtime |
| HG00276 ONT test data | `~/projects/pgx/chapter4_validation/` | `HG00276.haplotagged.bam` + `HG00276.filtered.vcf.gz` (untracked, ~3 GB) |
| HG001 ONT data | `~/projects/pgx/chapter8_validation/` | CYP2D6 region extract from ont-open-data giab_2025.12 |
| CYP2D6 truth set | `~/Documents/CYP2D6 Coriell.xlsx` | GeT-RM consensus for 179 Coriell samples (Gaedigk et al. 2019). Always check truth here |
| nallo test run | `~/projects/nallo_runs/hg00276_panel/` | HG00276 panel reads; `params.yaml`, `local_resources.config`, `publish_gvcf.config` |
| pgx run on nallo output | `~/projects/pgx_runs/hg00276_nallo/` | Config `nallo_run.yaml` (+ `config/config.yaml`) |
| Clair3 container | `~/singularity_images/clair3.sif` | Model `/opt/models/r1041_e82_400bps_sup_v410` (inside the SIF) |


## pgx pipeline: what has been changed

- The hydra-genetics remote `module` imports for alignment and variant calling were replaced by local rule files:
  - `ont_alignment.smk`: bwa mem for Illumina, `minimap2 -ax map-ont` for ONT
  - `ont_varcall.smk`: Clair3 (via Apptainer) for both platforms, replacing GATK HaplotypeCaller
- `@RG` tags use `SM:{sample}_{type}` in both aligners. Clair3 gets `--sample_name {sample}_{type}` to match.
- `variant_filtration.py` reads `DP` from the FORMAT field (Clair3), not INFO (GATK).
- The platform is called "ONT" everywhere (not "nanopore"). The HTML report shows "Sekvenseringsplattform" in Provdetaljer.
- Known accepted limitation: one platform per sample per run.
- StarPhase (`starphase_cyp2d6.smk`) is switched on with `starphase.enabled`. `starphase.bam_source` picks the input: `"nallo"` (nallo output paths), `"custom"` (per-sample `bam`/`vcf` under `starphase.custom_inputs`, keyed `{sample}_{type}`), or `"pgx"` (does not work yet: pgx's Clair3 VCF is temporary, uncompressed and unphased, and its BAM isn't haplotagged).
- New rule files for the nallo path (all off by default in `config/config.yaml`; Snakefile include order: … `starphase_cyp2d6`, `panel_coverage`, `gvcf_reference_calls`, `snp_genotypes`, `pharmcat`):
  - `pharmcat.smk`: `pharmcat_preprocess`, `pharmcat_outside_calls` (StarPhase TSV filtered to `pharmcat.outside_call_genes`, default CYP2D6), `pharmcat_run`, `generate_pharmcat_report` → `results/report/{sample}_{type}.pgx_report.html` (script `generate_pharmcat_report.py`, template `config/templates/html_template_pharmcat.html`, texts via `pharmcat_report_data.py`).
  - `panel_coverage.smk`: `samtools depth -a` over the panel targets; colour per target from the lowest depth at clinical (PharmCAT) positions, thresholds per platform.
  - `gvcf_reference_calls.smk`: expands nallo's gVCF at PharmCAT and SNP-marker positions and adds 0/0 (≥ 30 reads, confident) or ./. records to the phased VCF. `pharmcat_vcf()` feeds PharmCAT and the SNP section; StarPhase keeps the phased VCF.
  - `snp_genotypes.smk`: markers outside PharmCAT from `snp_markers.table`; indel markers are normalised (anchor base, left-aligned) against the reference before matching.
- Validation samples need their own `samples.tsv`/`units.tsv` (passed with `--config samples=... units=...`), because `common.smk` only accepts sample names listed there.
- Documentation: `README_nanopore_adaptation.md`, `pgx_ubuntu_setup_guide.md`, `pgx_nanopore_chapter2_guide.md`, `pgx_nallo_chapter3_guide.md`.

## Known quirks and workarounds

- **samtools:** use `/usr/bin/samtools` (v1.22.1) explicitly. The `pgx` conda env contains an ancient samtools (v0.1.19) that shadows it.
- **Running pgx on a run folder:** from the run folder, `conda run -n pgx snakemake -s ~/projects/pgx/workflow/Snakefile --configfile config/config.yaml nallo_run.yaml --cores 2 --use-singularity results/report/<sample>_<type>.pgx_report.html` (dry run with `-n` first).
- **Snakemake 7 does not rerun after script edits**, and does not rebuild a missing intermediate file while the final report exists. Ask for the intermediate target directly (e.g. `results/gvcf/<s>_<t>.with_refcalls.vcf.gz`), then the report.
- **PharmCAT discards indel reference records with ALT "."**: write PharmCAT's own ALT from `/pharmcat/pharmcat_positions.vcf` (done in `gvcf_reference_calls.py`).
- **Unit tests:** in `workflow/scripts`, `conda run -n pgx python -m unittest discover -p '*_test.py'` (plain `python -m unittest` finds 0 tests).
- **`conda run` does not pass stdin:** put Python in a file (scratchpad) instead of a heredoc.
- **nallo:** the machine has 8 cores, so cap with `local_resources.config`; `skip_sex_check` needs a known sex in the samplesheet; needs `--par_regions` and `sv_callers: "sniffles"`; nallo does not publish the gVCF without `publish_gvcf.config`; the Nextflow JVM can hang after "Pipeline completed" (watch the log).
- **Clair3:** run as `apptainer exec ~/singularity_images/clair3.sif run_clair3.sh ...`
- **Running pgx:** use `--use-singularity` (as in the README). Without it, wrapper-based rules such as `depth_of_baits` fail with `No module named 'snakemake_wrapper_utils'`.
- **StarPhase v1.5.0 CLI:** output goes to `--output-calls <file.json>`. There is no `--output-dir` or `--threads`. `--pharmcat-tsv <file>` writes the PharmCAT outside-call file. A full run on a WGS-sized BAM takes ~15 min (most of it on CYP2D6).
- **PharmCAT:** `pharmcat_pipeline` has no outside-call flag, so run two stages: `pharmcat_vcf_preprocessor -vcf ...`, then `java -jar /pharmcat/pharmcat.jar -vcf <preprocessed> -po <starphase tsv> -reporterJson`. Run with `apptainer exec --bind /home/lars-ola-arvidsson,/tmp`. **Never** use `-0` / `--absent-to-ref` / `--unspecified-to-ref` (they turn "not seen" into "normal"; clinical decision).
- **The pgx conda env has no openpyxl.** Read `.xlsx` files with the Python standard library (zipfile + XML) instead of installing anything.
- **StarPhase `build` is broken** (SIGFPE in v2.2.0; v1.5.0 needs the PharmVar API, and pharmvar.org is blocked on the institution network). Use the pre-built databases from https://github.com/PacificBiosciences/pb-StarPhase/tree/main/data instead.
- **Paraphase is dropped from the clinical path** (decided 2026-10-19). It reports copy number and gene copies but no star alleles, StarPhase can't use its output, and it is HiFi-first. The code is parked on branch `paraphase-research`. StarPhase is the CYP2D6 caller.
- **Aldy4** has a documented hang on ONT CYP2D6 (GitHub issue #73). Avoid it.
- Past environment issues to be aware of: AppArmor restrictions, libdeflate conda conflicts, PEP 668 pip restrictions (don't `pip install` into the system Python), pulp version conflicts with Snakemake 7.

## Validation status

CYP2D6 regression pair (truth = GeT-RM consensus from the Coriell xlsx):

| Sample | GeT-RM truth | StarPhase | Status |
|---|---|---|---|
| HG00276 | \*4/\*5 | \*5/\*4, PM, AS 0.0 | ✅ also through Snakemake (`bam_source: custom`, 2026-10-19) and through PharmCAT as outside call |
| HG001 / NA12878 | \*3/(\*68)+\*4 | \*3/[\*68 + \*4] | ✅ correct. The Chapter 8 document wrongly used "\*2x2/\*4" as truth |

- Full pgx end-to-end run succeeded for both Illumina and ONT. Integration test passes with StarPhase disabled.
- **First full nallo → pgx run (2026-10-26)** on HG00276 panel reads (ONT, reads overlapping the panel; `~/projects/nallo_runs/hg00276_panel`, pgx side `~/projects/pgx_runs/hg00276_nallo`):
  - CYP2D6 \*5/\*4 PM.
  - With gVCF reference calls and SNV calling ±500 bp, **every PharmCAT panel gene matches the original WGS result**, incl. CYP3A4 \*1/\*2 and UGT1A1 \*1/\*80+\*28. Without padding, CYP3A4\*2 was missed (reported \*1/\*1) and UGT1A1\*28 was missing.
  - DPYD \*3 (c.1898delC) and \*7 (c.295_298delTCAT) were "not tested" until 2026-10-27. The cause was **not** DeepVariant: they had a confident 0/0 with 75/87 reads, but PharmCAT's preprocessor discards indel reference records whose ALT is "." (43 positions). Fixed by writing PharmCAT's own ALT (`pharmcat_positions.vcf` from the container). With that and the 95% read rule, 20 more alleles became testable (DPYD \*3/\*5/\*7, NUDT15 \*6/\*9/\*22, TPMT \*47, CYP2C9 \*6/\*25/\*85, 8 G6PD variants, …); no diplotype changed.
  - Still not testable in DPYD: c.1682G>T, c.3049G>A, \*9B, c.1615G>A, c.61C>T, c.62G>A, c.632A>G (82–90% normal reads, below 95%).
  - The DPYD amplicon pipeline only covers ~300 bp around rs75017182.
  - Note: the comparison is against the original WGS analysis of the same sample, not a truth set. Non-CYP2D6 genes still need GeT-RM/GIAB truth.
  - SNP section (2026-10-26): CYP3A4 rs2740574 T/T (0 copies of \*1B), MTHFR c.1286A>C heterozygous, CYP3A5\*7 -/- (71 reads). IFNL4 rs11322783 is "Ej bedömbar (osäkert anrop)": 64–66/68 reads show TT and none ΔG, but DeepVariant left chr19:39248515 as ./. (GQ 1).

## Coverage colour scale (decided by Lars-Ola, 2026-10-20)

Read depth per target/position in the report, for ONT:

| Depth | Colour |
|---|---|
| ≥ 100 | green |
| 50–99 | pale green |
| 30–49 | orange |
| 1–29 | red |
| 0 | black |

Reasoning: nallo users accept depth 30 for HiFi (Q > 30, 99.9% accurate reads). ONT reads average about Q20 (99%), range ~Q10–Q30, so roughly 100× is the ONT target. 24 samples per run will not reach that everywhere, hence the graded scale. HiFi can use lower thresholds, so the thresholds are configurable per platform (`panel_coverage.thresholds`, key = `platform` in units.tsv).

PacBio HiFi (decided by Lars-Ola, 2026-10-27), platform `PACBIO`: < 10 failed (red; 0 black), 10–19 warning (orange), 20–40 pass (pale green), > 40 over-sequenced (green). Config `[41, 20, 10, 1]`. Note: the orange limit is also the minimum depth for a 0/0 call and for "lågt läsdjup", so HiFi uses 10 reads there.

## Current focus (October 2026, updated 2026-10-27)

The chain nallo → StarPhase → PharmCAT → pgx report is connected and runs end to end (2026-10-26). It is enabled per run with `starphase`, `pharmcat`, `panel_coverage`, `snp_markers` and `reference_calls` in the config, and `bam_source: nallo`. nallo templates are in `config/nallo/` (params, gVCF publishDir, `pad_bed.py`).

Data kept outside the repo (drafts / local):
- `~/reference/guidelines/swedish_guidelines.tsv` — Swedish texts, draft. Format: `config/templates/swedish_guidelines.example.tsv`.
- `~/reference/panel/pgx_panel_v2_targets.bed` — named panel targets (Region Örebro län Twist design TE-98775480).
- `~/reference/panel/snp_markers_draft.tsv` — 37 non-PharmCAT markers with HGVS and `defining_allele`.
- `~/reference/panel/pharmcat_3.4.0_positions.bed`.
- `~/reference/panel/pgx_adaptive_sampling_flank5kb.bed` (capture libraries) and `..._flank20kb.bed` (native gDNA): MinKNOW adaptive-sampling targets, built with `config/adaptive_sampling/make_as_bed.py` from the panel BED + PharmCAT positions + SNP markers (2026-10-28). Use with the chr-prefixed GRCh38 in MinKNOW; nallo preset `ONT_R10_AS`. Comparison against hybrid capture planned.

Report rules decided by Lars-Ola:
- Depth colours (above). A gene's colour comes from its targets with clinical positions.
- 0/0 only with ≥ 30 reads (ONT; the platform's orange limit).
- **Read rule (2026-10-27):** where the caller is unsure, ≥ 95% normal reads with ≥ 30 reads counts as normal (`reference_calls.read_rule`). A read is normal if it has the reference at every base of the position and no insertion/deletion inside it. Marked INFO `PGX_READRULE`/`PGX_READS` in the VCF; the SNP section says "normal enligt läsandel (x % av n läsningar)".
- SNP section: genotype + zygosity + copies of the defining allele.
- "Ej testade alleler" shown per gene.
- Not assessable is labelled by cause: "Ej bedömbar (lågt läsdjup)" = < 30 reads; "Ej bedömbar (osäkert anrop)" = ≥ 30 reads but the caller gave ./. (technical labelling, not a new threshold).
- HLA not analysed (StarPhase HLA unvalidated).
- English CPIC/DPWG text where no Swedish text exists.
- Clinically relevant = ClinPGx level ≥ 3. All markers stay until the hybrid-capture runs show what can be called.

Open items:
- IFNL4 ΔG (rs11322783) on HG00276 is still not assessable: only 84% of reads are normal across all three bases (ONT noise next to a G run). Same for 7 DPYD positions (82–90%).
- PharmCAT report: show which genes/positions were called by the read rule (only the SNP section says so today)
- review the marker table: scope of markers. (Clinical alleles #6/#13/#20/#31 approved 2026-10-28; #13 rs55965422 labelled "CYP3A5*3.005, tidigare *5", suballele number still to check against PharmVar. Possible later: flag rs55965422 without rs776746 on the same haplotype.)
- HLA validation
- "journal version" of the report (genotypes without drug text)
- suggest to nallo (GMS) to publish the gVCF by default
- run on real hybrid-capture data (2027)

## Longer-term (separate project)

A new open-source ONT CYP2D6/CYP2D7 resolver (MIT or GPL-3, to be contributed to GMS), roughly 180–300 hours of work. This is only needed if StarPhase proves insufficient. So far StarPhase is 2/2 on the structural CYP2D6 cases (HG00276, HG001).
