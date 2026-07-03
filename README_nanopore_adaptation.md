# pgx — Nanopore Adaptation
### Branch: `nanopore-adaptation`
### Maintainer: Lars-Ola Arvidsson, Region Halland
### Based on: [genomic-medicine-sweden/pgx](https://github.com/genomic-medicine-sweden/pgx)

---

## What this branch does

This branch adapts the `genomic-medicine-sweden/pgx` pharmacogenomics pipeline — originally
designed for Illumina short-read WGS data — to support Oxford Nanopore Technology (ONT)
long-read FASTQ input. Both platforms can coexist in the same pipeline run, with each sample
routed to the correct aligner and variant caller based on the `platform` column in `units.tsv`.

---

## Key changes from the original pipeline

| Component | Original | This branch |
|---|---|---|
| Alignment | bwa mem (via hydra-genetics module) | minimap2 -ax map-ont (ONT) / bwa mem (Illumina) |
| Duplicate marking | Picard MarkDuplicates | Removed — not applicable to long reads |
| Variant calling | GATK HaplotypeCaller (via hydra-genetics module) | Clair3 v2.0.2 (via Apptainer container) |
| VCF depth filter | INFO/DP field | FORMAT/DP field (Clair3 format), with INFO/DP fallback |
| Sequencing platform | Illumina only | ONT + Illumina, platform-aware routing |

---

## Requirements

### Software
- Ubuntu 20.04+ (tested on Ubuntu 24)
- Conda / Miniconda3
- Apptainer (installed via `sudo apt install apptainer`)
- System samtools ≥ 1.17 (installed via `sudo apt install samtools`)
- AppArmor unprivileged userns enabled:
  ```bash
  echo 'kernel.apparmor_restrict_unprivileged_userns=0' | sudo tee /etc/sysctl.d/99-userns.conf
  sudo sysctl --system
  ```

### Conda environment
```bash
conda create -n pgx python=3.8
conda activate pgx
pip install -r requirements.txt --break-system-packages
conda install -c bioconda minimap2
conda install -c bioconda bwa
conda install -c conda-forge singularity
```

### Clair3 container
```bash
mkdir -p ~/singularity_images
singularity pull ~/singularity_images/clair3.sif docker://hkubal/clair3:latest
```

The container includes bundled models. For R10.4.1 + Dorado SUP basecalling:
```
model_path: "/opt/models/r1041_e82_400bps_sup_v410"
```

For R9.4.1 + Guppy SUP:
```
model_path: "/opt/models/r941_prom_sup_g5014"
```

---

## How to run

```bash
conda activate pgx
cd /path/to/pgx/.tests/integration   # or your own run folder
snakemake -s ../../workflow/Snakefile \
  --use-singularity \
  --cores 6 \
  --configfile config/config.yaml
```

---

## Input file structure

### `samples.tsv`
One row per patient/sample. The `sample` value becomes the identifier used in all output
file names.

```
sample          tumor_content
KA0822402033    1.0
KA0822402034    1.0
```

Allowed characters in `sample`: letters, numbers, and underscores. Hyphens may cause issues
with wildcard constraints — use underscores if your LID contains hyphens (e.g. `KA_0822402033`).

### `units.tsv`
One row per sequencing run. Multiple rows per sample are supported for different platforms.

```
sample          type  platform   barcode  machine  flowcell  lane  fastq1                          fastq2                          adapter
KA0822402033    T     Illumina   ACTGAC   HV123    FC001     L001  fastq/KA0822402033_R1.fastq.gz  fastq/KA0822402033_R2.fastq.gz  AGATCGG...
KA0822402034    T     ONT        NANA     ONT001   ONT001    L001  fastq/KA0822402034.fastq.gz     fastq/KA0822402034.fastq.gz     AGATCGG
```

Platform values: `Illumina` or `ONT` (case-sensitive).

For ONT samples, `fastq2` and `adapter` must still be present but can be set to the same
value as `fastq1` and a dummy adapter string — the pipeline ignores them for ONT samples.

---

## ⚠️ Known limitations

### 1. One platform per sample per run (IMPORTANT)

**If the same `sample`+`type` combination appears in both an Illumina row and an ONT row
in `units.tsv`, the two alignment jobs will overwrite each other's output BAM file.**

The current output path is:
```
alignment/samtools_merge_bam/{sample}_{type}.bam
```

This path does not encode the platform, so it is shared between Illumina and ONT runs
for the same sample.

**Workaround:** only include one platform per sample in `units.tsv` at a time. To compare
Illumina and ONT results for the same patient, run the pipeline twice — once with the
Illumina row, once with the ONT row — keeping the results in separate output directories.

**Planned fix:** encode platform in the BAM output path
(`{sample}_{type}_{platform}.bam`) and add a merge/handoff rule so downstream steps
remain platform-agnostic. This will be implemented in a future version.

### 2. Clair3 model must match basecalling chemistry

The Clair3 model path in `config.yaml` must match the chemistry and basecalling model
used to generate the ONT reads:

| Chemistry | Basecaller | Model |
|---|---|---|
| R10.4.1 | Dorado SUP | `r1041_e82_400bps_sup_v410` |
| R10.4.1 | Dorado HAC | `r1041_e82_400bps_hac_v410` |
| R9.4.1 | Guppy SUP | `r941_prom_sup_g5014` |

Using the wrong model will reduce variant calling accuracy significantly. Check your
basecalling log or run metadata to confirm which model was used.

### 3. System samtools used for alignment

The conda-installed samtools (0.1.19) in the `pgx` environment is too old to support
modern `samtools sort -o` syntax. The alignment rule explicitly calls `/usr/bin/samtools`
(system-installed v1.22.1). If you move this pipeline to a different machine, ensure
`/usr/bin/samtools` is version ≥ 1.17.

### 4. Clair3 is applied to Illumina samples too

The current variant calling rule uses Clair3 for all samples regardless of platform.
Clair3 is an ONT-optimised tool and will use the ONT model even for Illumina reads.
This is functionally acceptable for pipeline testing and cross-platform comparison, but
for production Illumina use, GATK HaplotypeCaller would give better results.

**Planned fix:** add platform-aware variant calling, routing Illumina samples back to
HaplotypeCaller and ONT samples to Clair3.

### 5. Pipeline currently covers TPMT only (chr6)

The reference genome, BED files, and clinical guidelines currently target only the TPMT
gene on chromosome 6. Expanding to other PGx genes requires:
- A full reference genome (not just chr6)
- Updated BED files for each gene region
- Updated haplotype definition and clinical guideline files

**Planned expansion:** CYP2D6 (chr22), CYP2C19 (chr10), DPYD (chr1)

### 6. No CYP2D6 CNV calling yet

CYP2D6 copy number variation (CNV) is not yet detected. This is planned as a major
feature addition using `Sniffles2` for structural variant calling on ONT BAMs, combined
with a new HTML report section showing CYP2D6 gene structure relative to CYP2D7.

---

## Output files

For each sample, the pipeline produces:

```
results/dna/pgx/
├── {sample}_T.pgx_report.html              # Main clinical HTML report
├── {sample}_T.pgx_clinical_recommendations.txt  # Text recommendations
├── {sample}_T.pgx_depth_table.xlsx         # Coverage table
└── additional_files/
    ├── {sample}_T.detected_variants.csv    # Detected PGx variants
    ├── {sample}_T.depth_of_baits.gdf       # Bait coverage
    └── {sample}_T.depth_at_missing_annotated.gdf  # Coverage at uncalled positions
```

---

## Planned features (roadmap)

- [ ] Platform-in-BAM-path fix for true parallel Illumina+ONT per sample
- [ ] Platform-aware variant calling (HaplotypeCaller for Illumina, Clair3 for ONT)
- [ ] Add "Sekvenseringsplattform" field to HTML report
- [ ] Expand gene panel to CYP2D6, CYP2C19, DPYD, CYP2C9
- [ ] CYP2D6 CNV calling with Sniffles2
- [ ] CYP2D6 HTML report section with gene structure visualisation vs CYP2D7
- [ ] Full reference genome support (currently chr6 only)
- [ ] Validation against known CYP2D6 genotypes from ONT open data

---

## Git branch history

| Commit | Description |
|---|---|
| Initial adaptation | Replace hydra-genetics modules with local ONT rules |
| Platform-aware routing | Single rule routing Illumina→bwa, ONT→minimap2 |
| Clair3 integration | Apptainer-based Clair3 replacing HaplotypeCaller |
| End-to-end validation | Full 45-step run successful on both platforms |
| ONT rename | Renamed "nanopore" to "ONT" throughout |