# User-supplied diagnosis labels

Three of the 15 datasets have **no diagnosis in their NCBI metadata**, so the pipeline cannot
label them on its own:

| Dataset | Runs needed | Groups (Wang et al. 2026, Table 1) | File to create |
|---|---|---|---|
| PRJNA611839 | 65 | AD 21 / NC 44 | `PRJNA611839.csv` |
| PRJNA801673 | 60 | MCI 30 / NC 30 | `PRJNA801673.csv` |
| PRJNA690972 | 46 of the 80 runs | MCI 13 / NC 33 | `PRJNA690972.csv` |

Get the labels from the original publication (supplementary tables) or from its authors, and
save them as a CSV in this folder:

```csv
Run,Diagnosis
SRR11287963,AD
SRR11287964,NC
```

* `Run` is the SRA run accession. To join on another metadata column instead (e.g. `SampleName`),
  change `labels.join_column` in the dataset's config.
* `Diagnosis` must be one of the groups in the config (`AD`, `MCI`, `NC`).
* Only runs listed in the file are used. For PRJNA690972, list exactly the 46 runs the paper used.

The pipeline checks the run count and group counts against the config's `expected` block and
stops if they differ, so a wrong or incomplete file cannot slip through.
