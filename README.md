# microbiome_informed_disease_prediction

A framework for predicting disease from the gut microbiome:
**raw public sequencing data → genus-abundance matrix → statistics and ML → report.**

It covers the **15 public 16S rRNA datasets** of Alzheimer's disease (AD), mild cognitive
impairment (MCI) and normal controls (NC) used by Wang et al. 2026 (*mSystems* 11:e01754-25),
1,700 participants in total. They are processed the way that paper describes:
SRA Toolkit → FIGARO → DADA2 (QIIME 2) → VSEARCH against SILVA 138 (99%) → genus tables →
merged table with the paper's genus filter.

---

## Contents

```
Scripts/amplicon_16s/
  run_end_to_end.sh            ← one command per dataset: dry run → full run → analysis report
  run_16s.sh                   ← the preprocessing pipeline (activates the conda envs)
  preprocess_16s_amplicon.py   ← SRA → QIIME 2 / DADA2 / VSEARCH-SILVA → genus tables
  analyze_genus.py             ← diversity, PERMANOVA, differential abundance, ML baseline → HTML report
  merge_genus_tables.py        ← combine finished datasets into one matrix (paper's genus filter)
configs/                       ← one YAML per dataset (16 files: PRJNA734525 has two sequencing batches)
labels/                        ← diagnosis labels you supply for the 3 datasets without labels in NCBI
envs/
  setup_local.sh               ← one-time install (Linux / macOS / WSL)
  windows/setup_wsl.ps1        ← Windows: creates the Linux environment, then installs
docs/16S_PREPROCESSING.md      ← technical reference: steps, options, troubleshooting
results/PRJNA496408/           ← finished reference results of the first dataset
```

## Quick start

```bash
git clone https://github.com/kabir-ai2bio-lab/microbiome_AI.git
cd microbiome_AI
bash envs/setup_local.sh                                                         # 1. install (once)
bash Scripts/amplicon_16s/run_end_to_end.sh --dry-run-only                       # 2. verify (~10-15 min)
CONFIG=configs/PRJNA554111.yaml bash Scripts/amplicon_16s/run_end_to_end.sh      # 3. one dataset
```

Each dataset's results go to `outputs/<dataset>/`. The ML input is
`outputs/<dataset>/full/genus/genus_relative_with_provisional_labels.csv` and the report is
`outputs/<dataset>/full/analysis/analysis_report.html`.

---

## 1. The 15 datasets

Every config reproduces the sample selection of the paper's Table 1: the pipeline checks the
run count and group counts against it and stops if NCBI gives anything else.

| # | Dataset (config) | Samples | AD / MCI / NC | Reads | Notes |
|---|---|---|---|---|---|
| 1 | PRJNA496408 | 93 | 33 / 32 / 28 | V3–V4, 2×300 | **done**, results in `results/` |
| 2 | PRJNA633959 | 171 | 100 / – / 71 | V3–V4, 2×300 | |
| 3 | PRJNA489760 | 90 | 30 / 30 / 30 | V3–V4, 2×300 | stool runs only (the other 90 runs are blood) |
| 4 | PRJNA611839 | 65 | 21 / – / 44 | V3–V4, 2×301 | ⚠ labels needed ([section 5](#5-the-3-datasets-without-labels-in-ncbi)) |
| 5 | PRJNA554111 | 86 | 43 / – / 43 | V3–V4, 2×301 | |
| 6 | PRJNA946900 | 438 | 225 / 38 / 175 | V4–V5, 2×238 | mixed read orientation, re-oriented automatically |
| 7 | PRJNA792014 | 26 | 13 / – / 13 | V3–V4, 2×250 | |
| 8 | PRJNA855571 | 165 | 125 / – / 40 | V3–V4, single (merged) | gut runs only (the other 171 are oral) |
| 9 | PRJNA690972 | 46 | – / 13 / 33 | V3–V4, single (merged) | ⚠ labels and run list needed |
| 10 | PRJNA811324 | 84 | 41 / – / 43 | V3–V4, single (merged) | |
| 11 | PRJNA734525_batch1 + _batch2 | 62 + 63 | 47 / 27 / 51 | V3–V4, 2×232 and 2×281 | two sequencing batches, each denoised on its own |
| 12 | PRJEB51982 | 50 | 25 / – / 25 | V4, 2×251 | |
| 13 | PRJNA801673 | 60 | – / 30 / 30 | V4, 2×153 | ⚠ labels needed |
| 14 | PRJNA770746 | 99 | 45 / – / 54 | V4, 2×251 | stool runs only (the other 99 are oral) |
| 15 | PRJEB86129 | 102 | 51 / – / 51 | V4, single 251 | stool runs only (saliva and a blank excluded) |

**12 datasets run as they are.** The other 3 need a label file first ([section 5](#5-the-3-datasets-without-labels-in-ncbi)).

How each dataset's samples and labels were found (all from NCBI metadata, checked against Table 1):
sample/library name prefixes (e.g. `AD12` / `CON7`), or BioSample attributes such as
`disease`, `host disease status`, `patient_group` or the BioSample title. Every config explains
its own choice in its header comment.

## 2. What you need

| | Minimum | Recommended |
|---|---|---|
| OS | Linux, macOS (Intel or Apple Silicon), or Windows 10/11 (through WSL, set up by the helper script) | Linux |
| CPU | 4 cores | 8+ cores (every core is used automatically) |
| RAM | 16 GB | 32 GB |
| Disk | 30 GB free for one dataset; ~150 GB for all 15 | |
| Internet | to download the data (NCBI SRA) and the SILVA reference | |

**Run time.** The slowest step is taxonomy assignment (VSEARCH against SILVA), which grows
with the number of distinct sequences (ASVs) and shrinks with more cores. Measured on 16 threads
(Ryzen 7 5800X): PRJNA496408 (93 samples) took 93 min, of which VSEARCH was 71 min. Expect
roughly 1–4 h per dataset on 16 threads, about two days for all 15, and 3–4× longer on a
4-core laptop. The SILVA reference is downloaded once (~5 min) and shared by all datasets.

## 3. Install

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
| `figaro` | FIGARO (chooses DADA2 truncation lengths), pinned version, cloned to `~/opt/figaro` |

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
Linux shell and continue there:

```powershell
wsl -d BiomeFM -u root
cd ~/microbiome_AI
```

### Verify the installation

```bash
bash Scripts/amplicon_16s/run_end_to_end.sh --dry-run-only
```

This runs **every pipeline step** on 3 samples of PRJNA496408 (5,000 read pairs each).
DADA2 is deterministic, so a correct installation reproduces these numbers **exactly**:

| Sample | Group | Reads in | Reads after DADA2 |
|---|---|---|---|
| SRR8061802 | MCI | 4,992 | 3,558 |
| SRR8061776 | AD | 4,990 | 3,203 |
| SRR8061756 | NC | 4,994 | 3,082 |

That gives **71 ASVs**, and the log prints `REPRODUCIBILITY CHECK PASSED`. If the numbers
differ, check that QIIME 2 2026.7 was installed.

## 4. Run the datasets

### One dataset

```bash
CONFIG=configs/PRJNA554111.yaml bash Scripts/amplicon_16s/run_end_to_end.sh --dry-run-only   # ~10 min check
CONFIG=configs/PRJNA554111.yaml bash Scripts/amplicon_16s/run_end_to_end.sh                  # full run + report
```

The dry run processes 2–3 samples (one per group) through every step and stops on any
problem (wrong sample counts, primers not found, too few reads kept), so start with it.
The full run then downloads every sample, denoises all of them in **one joint DADA2 run**,
assigns taxonomy and writes the genus tables and the analysis report.

### All datasets

```bash
for c in PRJNA792014 PRJEB51982 PRJNA811324 PRJNA554111 PRJNA734525_batch1 PRJNA734525_batch2 \
         PRJNA489760 PRJNA770746 PRJEB86129 PRJNA855571 PRJNA633959 PRJNA946900; do
  CONFIG=configs/$c.yaml bash Scripts/amplicon_16s/run_end_to_end.sh || echo "FAILED: $c"
done
```

This is the 12 datasets that have labels, smallest first. PRJNA496408 is already done; add it to
reproduce it. Add the other 3 after creating their label files ([section 5](#5-the-3-datasets-without-labels-in-ncbi)).

* **Long runs:** start them inside `tmux` (or `screen`) so they survive closing the terminal.
* **Interrupted?** Run the same command again. Finished steps (downloads, QIIME
  artifacts, SILVA) are skipped automatically.
* **Options** (environment variables): `THREADS=8` limits CPU use, `EXPORT_DIR=/path`
  changes the output folder (default `outputs/`), and `COMPARE_TO=/path/to/matrix.csv`
  compares the result with a previous one.

### Merge the datasets into one matrix

```bash
conda activate rachis-qiime2-2026.7
python Scripts/amplicon_16s/merge_genus_tables.py            # every finished dataset in outputs/
```

This writes `outputs/merged/`: one samples × genera table with `Dataset`, `Project` and
`Diagnosis_provisional` columns, before and after the paper's filter (keep genera present in at
least 2 samples and at least 0.001% of all reads). The two PRJNA734525 batches are joined here.
The paper then removes study effects with MMUPHin (R package); that step is not part of this
repository yet.

## 5. The 3 datasets without labels in NCBI

PRJNA611839, PRJNA801673 and PRJNA690972 have no diagnosis anywhere in their NCBI metadata.
Their configs are ready, but the pipeline stops with a clear message until you add a label file
(`labels/<dataset>.csv` with columns `Run,Diagnosis`) taken from the original publication or its
authors. See [labels/README.md](labels/README.md).

## 6. Results

`outputs/<dataset>/full/` (same layout as the reference results in `results/PRJNA496408/full/`):

| File | What it is |
|---|---|
| `genus/genus_relative_with_provisional_labels.csv` | **ML input**: samples × genera (relative abundance) + `Diagnosis_provisional` |
| `genus/genus_counts.csv`, `genus_relative_abundance.csv` | the same without labels (counts / proportions) |
| `genus/*_filtered_single_project.csv` | genera present in ≥2 samples and ≥0.001% of reads |
| `analysis/analysis_report.html` | **report**: QC, composition, alpha/beta diversity, differential abundance, ML baseline |
| `analysis/*.csv` | every table in the report |
| `dada2_per_sample_stats.csv`, `pre_dada2_qc.csv` | read counts through each step, per sample |
| `asv_taxonomy_consensus.csv` | taxonomy for every ASV |
| `figaro/` | FIGARO report; `dada2_params.json` = the truncation / maxEE DADA2 used |
| `*.qzv` | QIIME 2 quality plots, viewable at <https://view.qiime2.org> |
| `run_provenance.json` | every parameter and tool version used (`dada2_used` = truncation actually applied) |

The analysis report contains:

1. **QC**: read depth and DADA2 retention by group
2. **Composition**: top genera per group
3. **Alpha diversity**: Shannon, Simpson, observed genera (rarefied), Kruskal-Wallis / Mann-Whitney
4. **Beta diversity**: Bray-Curtis PCoA, PERMANOVA (overall and pairwise)
5. **Differential abundance**: CLR transform, nonparametric tests, BH-FDR
6. **ML baseline**: random forest and logistic regression for every pair of groups (and 3-class
   when there are three), with repeated stratified CV AUC and a label-permutation test
7. **Reproduction check** against a reference matrix (with `COMPARE_TO`)

To re-run only the analysis on an existing result:

```bash
conda activate rachis-qiime2-2026.7
python Scripts/amplicon_16s/analyze_genus.py --cohort-dir ~/biomefm_runs/<dataset>/cohorts/full
```

## Reference results (PRJNA496408)

Produced on 2026-10-22 by the full run (16 threads, 93 min) and committed in
[`results/PRJNA496408/full/`](results/PRJNA496408/full/). The genus matrix is
**identical** to the original Colab run: all 93 samples, every value and every label.

| | |
|---|---|
| Samples | 93 (AD 33, MCI 32, NC 28; provisional labels) |
| Reads | 4,295,280 raw pairs → 3,054,807 after DADA2 (median retention 71.96%, range 48.7–78.9%) |
| Features | 2,002 ASVs → 224 named genera (162 after the prevalence/abundance filter) |
| ML input | [`genus/genus_relative_with_provisional_labels.csv`](results/PRJNA496408/full/genus/genus_relative_with_provisional_labels.csv) (93 × 225 + label) |

| Analysis | Result |
|---|---|
| Alpha diversity | Shannon differs between groups (Kruskal-Wallis p = 0.015); observed genera p = 0.083 |
| Beta diversity | Bray-Curtis PERMANOVA R² = 0.063, p = 0.001; all pairwise comparisons p < 0.05 |
| Differential abundance | no genus significant after BH-FDR; top candidates *[Eubacterium] hallii group*, *Bacteroides*, *[Eubacterium] brachy group*, *Klebsiella*, *Alistipes* |
| ML baseline (5-fold CV × 10, logistic L2) | AD vs MCI AUC 0.80 ± 0.11 (perm p = 0.005); MCI vs NC 0.74 ± 0.12 (p = 0.02); AD vs NC 0.64 ± 0.16 (p = 0.08); 3-class 0.72 ± 0.08 (p = 0.005) |
| ML baseline (random forest) | AD vs MCI 0.74; MCI vs NC 0.69; AD vs NC 0.56; 3-class 0.64 |

Sequencing depth also differs between groups (Kruskal-Wallis p = 0.045), which can
contribute to the diversity and ML differences. With 93 samples and provisional labels,
these are exploratory baselines, not validated findings.

## 7. Adding a new dataset

Copy the config of a similar dataset (same read layout and primers) to
`configs/<BioProject>.yaml` and edit it. The keys:

| Key | Meaning |
|---|---|
| `project` / `dataset` | NCBI BioProject; optional `dataset` name for the output folder (e.g. one sequencing batch) |
| `layout` | `paired` (default) or `single` |
| `expected` | run count and group counts from the publication; the pipeline stops if the metadata gives anything else |
| `labels` | the metadata column with the diagnosis, and either `regex` (+ `map`, e.g. `CON` → `NC`) or `patterns` ({group: regex}); or `file` for labels you supply |
| `biosample_attributes` | `true` when the diagnosis or sample type is only in the NCBI BioSample record |
| `include` / `include_regex` | keep only matching runs, e.g. `isolation_source: feces` when a project also has oral or blood samples |
| `primers` | primers at the start of R1 / R2 (`r1`, `r2`), optional `r1_3prime` for reads spanning the whole amplicon, or `null` when the reads start after the primer |
| `reorient` | `true` for libraries sequenced in both orientations |
| `dada2.truncation` | `figaro` (paired-end: FIGARO's best truncation and maxEE, as in the paper) or `config` (`trunc_f` / `trunc_r` / `max_ee_*`, needed for single-end) |
| `figaro` | amplicon length without primers, and primer lengths including any spacer before them |
| `expected_results` | reference numbers; a matching run logs `REPRODUCIBILITY CHECK PASSED` |

Check the primers against a few raw reads (`fastq-dump -X 3000 --split-files <run>`), then
run the dry run and look at `pre_dada2_qc.csv` and `figaro/` before the full run.

## Caveats

* Diagnosis labels are **inferred** from NCBI metadata (sample names or BioSample attributes),
  not from verified clinical records. Their counts match the paper's Table 1.
* Primers and amplicon regions were inferred from the reads, not from the original protocols.
* The paper does not report DADA2 or VSEARCH settings. Paired-end datasets use FIGARO's choice;
  single-end datasets (FIGARO cannot model them) use full-length reads or a quality-based
  truncation stated in their config. PRJNA496408 keeps the fixed settings of the original
  Colab run (250/240), so its results stay identical to that run.
* Genus assignment uses a permissive 90% VSEARCH identity.
* Cross-study batch correction (MMUPHin in the paper) is not included yet, so treat
  single-dataset p-values and AUCs as exploratory.
