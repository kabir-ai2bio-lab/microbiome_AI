# microbiome_informed_disease_prediction

A framework for predicting disease from the gut microbiome:
**raw public sequencing data → genus-abundance matrix → statistics and ML → report.**

The first supported dataset is **PRJNA496408** (16S rRNA amplicon, 93 stool samples:
33 Alzheimer's disease, 32 MCI, 28 controls). Any paired-end Illumina 16S BioProject
on NCBI SRA can be added with a config file.

---

## Contents

```
Scripts/
  amplicon_16s/
    run_end_to_end.sh            ← one command: dry run → full run → analysis
    run_16s.sh                   ← run the preprocessing pipeline (activates conda envs)
    preprocess_16s_amplicon.py   ← SRA → QIIME 2 / DADA2 / VSEARCH-SILVA → genus tables
    analyze_genus.py             ← diversity, PERMANOVA, differential abundance, ML → HTML report
  data_selection.py              ← dataset selection
  preprocess_shotgun_data.py     ← shotgun metagenomics → MetaPhlAn / MVIB inputs
configs/
  PRJNA496408.yaml               ← all parameters + reference numbers for verification
envs/
  setup_local.sh                 ← one-time install (Linux / macOS / WSL)
  windows/setup_wsl.ps1          ← Windows: creates the Linux environment, then installs
  wslconfig.example              ← optional WSL2 CPU / RAM settings
docs/
  16S_PREPROCESSING.md           ← technical reference: steps, parameters, outputs, troubleshooting
```

## Quick start

```bash
git clone https://github.com/kabir-ai2bio-lab/microbiome_AI.git
cd microbiome_AI
bash envs/setup_local.sh                                       # 1. install (once, ~20-40 min)
bash Scripts/amplicon_16s/run_end_to_end.sh --dry-run-only     # 2. verify (~10-15 min)
bash Scripts/amplicon_16s/run_end_to_end.sh                    # 3. full run + analysis
```

Results are written to `results/PRJNA496408/`. Open
`results/PRJNA496408/full/analysis/analysis_report.html` in a browser.

---

## 1. What you need

| | Minimum | Recommended |
|---|---|---|
| OS | Linux, macOS (Intel or Apple Silicon), or Windows 10/11 (through WSL, set up by the helper script) | Linux |
| CPU | 4 cores | 8+ cores (every core is used automatically) |
| RAM | 16 GB | 32 GB |
| Disk | 30 GB free | |
| Internet | to download the data (NCBI SRA) and the SILVA reference | |

**Approximate run time for the full 93-sample dataset:**

| Machine | Time |
|---|---|
| 16 threads (Ryzen 7 5800X, Windows + WSL1) | about 1.5 h |
| 10 CPUs (Google Colab) | 3–4 h |
| 4 CPUs (laptop) | roughly 5–6 h (estimate) |

The slowest step is taxonomy assignment (VSEARCH against SILVA), which scales with the
number of cores. The dry run takes 10–15 minutes on any machine; about 5 of those
minutes are a one-time download of the SILVA reference.

## 2. Install

### Linux or macOS

```bash
git clone https://github.com/kabir-ai2bio-lab/microbiome_AI.git
cd microbiome_AI
bash envs/setup_local.sh
```

This installs into your home folder. No admin rights are needed, except `apt` packages
on Ubuntu and Rosetta 2 on Apple Silicon Macs (QIIME 2 only ships Intel macOS builds).
It creates:

| conda env | contents |
|---|---|
| `rachis-qiime2-2026.7` | QIIME 2 2026.7 with cutadapt, DADA2, feature-classifier, RESCRIPt, taxa |
| `biomefm-tools` | SRA Toolkit (`prefetch`, `fasterq-dump`), `pigz` |
| `figaro` | FIGARO (truncation-length advisor), pinned version, cloned to `~/opt/figaro` |

The script can be re-run safely. Anything already installed is skipped.

### Windows

QIIME 2 does not run natively on Windows, so it runs inside Linux (WSL). In PowerShell,
from the repository folder:

```powershell
powershell -ExecutionPolicy Bypass -File envs\windows\setup_wsl.ps1
```

The helper creates an Ubuntu environment called `BiomeFM` and installs everything inside it.
If CPU virtualization is turned off in the BIOS, it falls back to WSL1 automatically.
If Windows asks for a reboot, reboot and run the script again. Afterwards, open the
Linux shell and continue with step 3:

```powershell
wsl -d BiomeFM -u root
cd ~/microbiome_AI
```

## 3. Verify the installation (dry run)

```bash
bash Scripts/amplicon_16s/run_end_to_end.sh --dry-run-only
```

This runs **every pipeline step** on 3 samples (one per group, 5,000 read pairs each).
DADA2 is deterministic, so a correct installation reproduces these numbers **exactly**:

| Sample | Group | Reads in | Reads after DADA2 |
|---|---|---|---|
| SRR8061802 | MCI | 4,992 | 3,558 |
| SRR8061776 | AD | 4,990 | 3,203 |
| SRR8061756 | NC | 4,994 | 3,082 |

That gives **71 ASVs**. The log prints `REPRODUCIBILITY CHECK PASSED` when everything
matches. If the numbers differ, check that QIIME 2 2026.7 was installed.

## 4. Full run

```bash
bash Scripts/amplicon_16s/run_end_to_end.sh
```

This downloads all 93 runs and processes them in **one joint DADA2 run**, then assigns
taxonomy and builds the genus tables and the analysis report. It reproduces the original
Colab run exactly: 4,295,280 → 3,054,807 reads (median retention 71.96%), 2,002 ASVs.

* **Long runs:** start it inside `tmux` (or `screen`) so it survives closing the terminal.
* **Interrupted?** Run the same command again. Finished steps (downloads, QIIME
  artifacts, SILVA) are skipped automatically.
* **Options** (environment variables): `THREADS=8` limits CPU use, `SKIP_FIGARO=1`
  skips FIGARO, `EXPORT_DIR=/path` changes the results folder, and
  `COMPARE_TO=/path/to/matrix.csv` compares against a previous result.

Individual steps can also be run directly. See `bash Scripts/amplicon_16s/run_16s.sh --help`
and [docs/16S_PREPROCESSING.md](docs/16S_PREPROCESSING.md).

## 5. Results

`results/PRJNA496408/full/`:

| File | What it is |
|---|---|
| `genus/genus_relative_with_provisional_labels.csv` | **ML input**: 93 samples × genera (relative abundance) + `Diagnosis_provisional` |
| `genus/genus_counts.csv`, `genus_relative_abundance.csv` | the same without labels (counts / proportions) |
| `genus/*_filtered_single_project.csv` | genera present in ≥2 samples and ≥0.001% of reads |
| `analysis/analysis_report.html` | **report**: QC, composition, alpha/beta diversity, differential abundance, ML baseline |
| `analysis/*.csv` | every table in the report |
| `dada2_per_sample_stats.csv`, `pre_dada2_qc.csv` | read counts through each step, per sample |
| `asv_taxonomy_consensus.csv` | taxonomy for every ASV |
| `figaro/` | FIGARO truncation suggestions (advisory) |
| `*.qzv` | QIIME 2 quality plots, viewable at <https://view.qiime2.org> |
| `run_provenance.json` | every parameter and tool version used |

The analysis report contains:

1. **QC**: read depth and DADA2 retention by group
2. **Composition**: top genera per group
3. **Alpha diversity**: Shannon, Simpson, observed genera (rarefied), Kruskal-Wallis / Mann-Whitney
4. **Beta diversity**: Bray-Curtis PCoA, PERMANOVA (overall and pairwise)
5. **Differential abundance**: CLR transform, nonparametric tests, BH-FDR
6. **ML baseline**: random forest and logistic regression for AD vs NC, MCI vs NC,
   AD vs MCI and 3-class, with repeated stratified CV AUC and a label-permutation test
7. **Reproduction check** against a reference matrix (with `COMPARE_TO`)

To re-run only the analysis on an existing result:

```bash
conda activate rachis-qiime2-2026.7
python Scripts/amplicon_16s/analyze_genus.py --cohort-dir ~/biomefm_runs/PRJNA496408/cohorts/full
```

## 6. Adding another dataset

Copy `configs/PRJNA496408.yaml` to `configs/<BioProject>.yaml` and edit it: project ID,
expected sample and group counts, how labels are parsed from the SRA metadata, primers,
truncation lengths. Then:

```bash
CONFIG=configs/<BioProject>.yaml bash Scripts/amplicon_16s/run_end_to_end.sh
```

Run the dry run first and look at `pre_dada2_qc.csv` and the FIGARO report before choosing
truncation lengths. The pipeline stops if primers are not found in most reads.

## Caveats

* Diagnosis labels are **inferred** from NCBI `LibraryName` prefixes (AD/MCI/NC), not from
  verified clinical records.
* Primers and amplicon length were inferred from the data, not taken from a published protocol.
* Genus assignment uses a permissive 90% VSEARCH identity.
* This is a single BioProject without cross-study batch correction, so treat p-values and
  AUCs as exploratory.
