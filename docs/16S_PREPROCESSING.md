# 16S amplicon pipeline: technical reference

Installation and the normal workflow are covered in the [README](../README.md). This page
describes what each step does, every option, and how to troubleshoot.

```
SRA runs ─► prefetch / fasterq-dump ─► QIIME 2 import ─► Cutadapt (primers)
        ─► QC gate ─► [FIGARO] ─► joint DADA2 (all samples at once)
        ─► VSEARCH vs SILVA 138 NR99 ─► collapse to genus ─► CSV matrices ─► analysis report
```

## Steps

`preprocess_16s_amplicon.py` runs these steps in order. `--until STEP` stops after a step.

| Step | What it does | Gate / output |
|---|---|---|
| `metadata` | downloads NCBI RunInfo; parses provisional labels from `LibraryName` | stops unless 93 paired runs with groups AD 33 / MCI 32 / NC 28 |
| `select` | dry-run: 3 stratified runs incl. pilot SRR8061802; validation: 10; full: all | stops on mixed platform/strategy; `analysis_recipe.json` blocks mixing parameters |
| `download` | `prefetch` + `fasterq-dump` + `pigz`, 4 runs in parallel; dry-run keeps the first 5,000 pairs | gzip integrity check |
| `import` | QIIME 2 manifest import, `demux summarize` | `demux-summary.qzv` |
| `cutadapt` | removes primers (806R-like on R1, 338F-like on R2), discards pairs without primers | `trimmed-summary.qzv` |
| `qc` | per sample: primer retention ≥ 80%, ≥ 95% of pairs long enough for 250/240 truncation | stops if any sample fails |
| `figaro` | (`--figaro`) FIGARO on raw reads normalized to one length per direction | advisory only, DADA2 settings are not changed |
| `dada2` | **one** `denoise-paired` across all samples (trunc 250/240, maxEE 4/6) | |
| `dada2-stats` | per-sample retention table | stops if any sample keeps < 20% |
| `reference` | SILVA 138 SSURef NR99 via RESCRIPt, reverse-transcribed; cached | |
| `taxonomy` | `classify-consensus-vsearch` (identity 0.90, coverage 0.80, maxaccepts 10, both strands) | |
| `genus` | collapse to level 6, write count / relative / labelled / filtered matrices | |
| `finalize` | genus totals must equal DADA2 non-chimeric reads; compares with `expected_results`; provenance; export | |

All parameters are in `configs/PRJNA496408.yaml`.

## Options

```
run_16s.sh [--config FILE] [--mode dry-run|validation|full] [--threads N] [--until STEP]
           [--figaro] [--workdir DIR] [--reference-dir DIR] [--export-dir DIR] [--preflight]
```

| Option | Default | Meaning |
|---|---|---|
| `--mode` | `dry-run` | `dry-run`: 3 samples × 5,000 pairs; `validation`: 10 samples; `full`: all |
| `--threads` | all CPUs | used by fasterq-dump, pigz, Cutadapt, DADA2, VSEARCH, and the Python QC workers |
| `--workdir` | `~/biomefm_runs/<project>` | downloads, QIIME artifacts, logs |
| `--reference-dir` | `<workdir>/reference/silva138` | share one SILVA build between projects |
| `--export-dir` | none | copy final tables, reports and provenance here |
| `--preflight` | | check tools, plugins, CPU, RAM and disk; run nothing |

`run_end_to_end.sh` wraps these: preflight → dry run + analysis → full run → analysis. It
reads `CONFIG`, `WORKDIR`, `EXPORT_DIR` (default `outputs/`), `THREADS`, `SKIP_FIGARO` and
`COMPARE_TO` from the environment.

## Resuming and reproducibility

* Rerun the same command after an interruption. Finished outputs are detected and skipped.
  QIIME writes to `*.partial.qza` and renames only on success, so a killed step never
  looks finished.
* `analysis_recipe.json` stops the pipeline from reusing a cohort folder made with
  different samples or parameters. Move `cohorts/<mode>/` aside to change them.
* `expected_results` in the config holds reference numbers (dry run and full run).
  `finalize` writes `expected_results_check.json` and logs `REPRODUCIBILITY CHECK PASSED`
  or a warning listing the differences.
* `step_timings.json` records seconds per step.

## Measured run times (16 threads, WSL1, 32 GB RAM)

| Step | Dry run (3 × 5k pairs) | Full (93 samples, 4.3 M pairs) |
|---|---|---|
| download | < 1 min | 2.5 min |
| import + summaries | < 1 min | 3 min |
| Cutadapt + summary | < 1 min | 7 min |
| QC gate + FIGARO | < 1 min | < 1 min |
| DADA2 | 37 s | 8.4 min |
| SILVA reference (once) | 5 min | cached |
| VSEARCH taxonomy | 3 min | 71 min |
| genus + finalize | 1 min | 1 min |
| **total** | **about 12 min** | **93 min** (+ 3 min analysis) |

VSEARCH is an exhaustive search of about 2,000 ASVs against about 500k SILVA sequences,
so it dominates the run and scales almost linearly with cores.

## Outputs

See the README for the main files. Everything lives under
`<workdir>/cohorts/<mode>/`: QIIME artifacts (`*.qza`), exported tables, `logs/` (one
log per command plus a pipeline log), `figaro/`, `genus/`, and `analysis/` when
`analyze_genus.py` runs.

## Platform notes

* **Windows:** runs inside WSL. `envs/windows/setup_wsl.ps1` picks WSL2 when CPU
  virtualization is enabled in the BIOS/UEFI (*SVM* on AMD, *VT-x* on Intel), and
  otherwise WSL1. Keep data inside the Linux filesystem (`~/...`); `/mnt/c` is much
  slower. For WSL2, `envs/wslconfig.example` raises the default limits (half the RAM).
  WSL1 always sees all CPUs and RAM.
* **macOS:** QIIME 2 2026.7 provides `osx-64` builds only. On Apple Silicon,
  `setup_local.sh` installs Rosetta 2 and creates the env with `CONDA_SUBDIR=osx-64`.
* **Linux:** tested on Ubuntu 24.04.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `wsl --import` / `--install` error 0x80370102 | virtualization is off in the BIOS; the helper uses WSL1 instead |
| `WSL1 is not supported` | the Windows feature is not enabled yet; run `setup_wsl.ps1` and reboot |
| preflight: `qiime MISSING` | use `run_16s.sh` (it activates the env) or `conda activate rachis-qiime2-2026.7` |
| pip / PyPI `NameResolutionError` inside WSL | temporary DNS issue; rerun. FIGARO's dependencies come from conda for this reason |
| `Partial ... outputs` | a previous run left some outputs of a multi-output step; delete the listed files |
| QC gate fails | inspect `pre_dada2_qc.csv`; check primers/orientation, don't just lower thresholds |
| VSEARCH killed (out of memory) | raise WSL2 `memory`, or use `--threads` to lower memory use |
| `REPRODUCIBILITY CHECK` warning | compare tool versions in `run_provenance.json` with QIIME 2 2026.7 |
