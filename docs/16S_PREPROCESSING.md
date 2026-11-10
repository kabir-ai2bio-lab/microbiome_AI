# 16S amplicon pipeline: technical reference

Installation and the normal workflow are covered in the [README](../README.md). This page
describes what each step does, every option, and how to troubleshoot.

```
SRA runs ─► prefetch / fasterq-dump ─► QIIME 2 import ─► Cutadapt (primers)
        ─► FIGARO (truncation) ─► QC gate ─► joint DADA2 (all samples at once)
        ─► VSEARCH vs SILVA 138 NR99 ─► collapse to genus ─► CSV matrices ─► analysis report
```

## Steps

`preprocess_16s_amplicon.py` runs these steps in order. `--until STEP` stops after a step.

| Step | What it does | Gate / output |
|---|---|---|
| `metadata` | downloads NCBI RunInfo (+ BioSample attributes if `biosample_attributes: true`); keeps paired runs and the `include` filter; parses labels with `labels.regex` | stops unless run and group counts equal `expected` (the paper's Table 1) |
| `select` | dry-run: `n_samples` stratified runs (2–3, one per group); validation: 10; full: all | stops on mixed platform/strategy; `analysis_recipe.json` blocks mixing parameters |
| `download` | `prefetch` + `fasterq-dump --split-3` (unpaired reads go to a separate file, so R1/R2 stay in step) + `pigz`, 4 runs in parallel; dry-run keeps the first 5,000 pairs; `reorient: true` swaps mates of pairs sequenced in the reverse orientation | gzip integrity check; share of swapped pairs logged |
| `import` | QIIME 2 manifest import, `demux summarize` | `demux-summary.qzv` |
| `cutadapt` | removes the config's primers, discards pairs without primers; skipped when `primers: null` (reads start after the primer) | `trimmed-summary.qzv` |
| `figaro` | (`--figaro`) FIGARO on raw reads normalized to one length per direction (the 1st-percentile read length; the 5th percentile when a few short-tailed runs would leave too little overlap, e.g. PRJNA734525_batch1) | with `dada2.truncation: figaro` its best choice sets DADA2 truncation (minus primer length) and maxEE; otherwise advisory |
| `qc` | per sample: primer retention ≥ 80%, ≥ 95% of pairs long enough for the truncation lengths | stops if any sample fails |
| `dada2` | **one** `denoise-paired` across all samples (truncation/maxEE from FIGARO or the config) | values used are in `run_provenance.json` → `dada2_used` |
| `dada2-stats` | per-sample retention table | stops if any sample keeps < 20% |
| `reference` | SILVA 138 SSURef NR99 via RESCRIPt, reverse-transcribed; cached | |
| `taxonomy` | `classify-consensus-vsearch` (identity 0.90, coverage 0.80, maxaccepts 10, both strands). With `asv_filter: {min_samples: n}` it first keeps only ASVs found in ≥ n samples (`*_min{n}samples` table/rep-seqs) | with `asv_filter`: `asv_filter_report.csv` (ASVs and reads kept per sample), `asv_filter.json` |
| `genus` | collapse to level 6, write count / relative / labelled / filtered matrices | |
| `finalize` | genus totals must equal DADA2 non-chimeric reads (the reads kept by `asv_filter`, when used); compares with `expected_results`; provenance; export | |

All parameters are in `configs/<dataset>.yaml` (keys explained in the README, section 7).

**Single-end datasets** (`layout: single`) use the same steps with single-end QIIME 2 types:
`cutadapt trim-single` (optionally a 3' primer via `r1_3prime`), no FIGARO (it models read-pair overlap),
`dada2 denoise-single` with the config's `trunc_f` (0 = keep full length, used for reads the
submitters already merged) and optional `trim_left`.

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
| `--reference-dir` | `<workdir>/../reference/silva138` | SILVA build, shared by all projects |
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
so it dominates the run and scales almost linearly with cores. It costs about 3.3–3.6 s per ASV
for V3–V4 reads on 16 threads, so the dry run's seconds per ASV predict the full run well. The
full runs of the 14 labelled datasets took about 63 h of compute on one 16-thread PC (0.4–12.5 h
each; taxonomy was most of every run). Per-dataset times are in `results/<dataset>/full/step_timings.json`.

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
| `PermissionError: [Errno 13] ... /tmp/qiime2/.../processes/` (WSL1) | a Windows process (usually the virus scanner) briefly locked a file in QIIME's cache; the pipeline retries this automatically (3 attempts). If it still fails, rerun the same command, or exclude the WSL folder from real-time scanning |
| `Partial ... outputs` | a previous run left some outputs of a multi-output step; delete the listed files |
| QC gate fails | inspect `pre_dada2_qc.csv`; check primers/orientation, don't just lower thresholds |
| VSEARCH killed (out of memory) | raise WSL2 `memory`, or use `--threads` to lower memory use |
| `REPRODUCIBILITY CHECK` warning | compare tool versions in `run_provenance.json` with QIIME 2 2026.7 |
