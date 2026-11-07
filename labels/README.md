# Diagnosis labels from the original publications

Three of the 15 datasets have **no diagnosis in their NCBI metadata**. Their labels come from the
original papers (the studies Wang et al. 2026 cite for them):

| Dataset | File | Source | Result |
|---|---|---|---|
| PRJNA801673 | `PRJNA801673.csv` | McLeod et al. 2023, *PLoS One* 18:e0280211, S1 Table (`RRC_ID` = SRA SampleName, `MCI` Y/N) | 30 MCI / 30 NC, exactly as in Wang et al. Table 1 |
| PRJNA690972 | `PRJNA690972.csv` | Duan et al. 2021, *Front Cell Infect Microbiol* 11:592842, Supplementary Table 1 (sample → NC / SCD / aMCI), matched to SRA names (G03A → G3) | 13 aMCI → MCI, 33 NC, exactly as in Wang et al.; SCD and 8 runs absent from the table are not used |
| PRJNA611839 | `PRJNA611839.csv` (**missing**) | Xi et al. 2021, *BMC Microbiol* 21:226: no per-sample table is published | needs labels from the authors |

## File format

```csv
Run,Diagnosis
SRR11287963,AD
SRR11287964,NC
```

* The join column is set by `labels.join_column` in the dataset's config (`Run` or `SampleName`).
* `Diagnosis` must be one of the groups in the config (`AD`, `MCI`, `NC`).
* Only runs listed in the file are used.

The pipeline checks the run count and group counts against the config's `expected` block (the
paper's Table 1) and stops if they differ, so a wrong or incomplete file cannot slip through.

## PRJNA611839

NCBI stores only numeric sample codes (e.g. `238439`), and neither NCBI, ENA nor the paper's
supplementary files link them to AD or control. The labels (21 AD / 44 NC) have to be requested
from the corresponding author of Xi et al. 2021. Save them as `PRJNA611839.csv` with columns
`Run,Diagnosis` (or the numeric codes with `labels.join_column: SampleName`).
