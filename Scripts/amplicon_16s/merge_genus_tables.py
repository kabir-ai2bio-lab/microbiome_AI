#!/usr/bin/env python3
"""
Merge the genus tables of several finished datasets into one matrix, as in Wang et al. 2026:
union of genera (missing = 0), then drop genera present in only one sample and genera below
0.001% of all reads. Batch correction (the paper uses MMUPHin) is a separate, later step.

    python Scripts/amplicon_16s/merge_genus_tables.py                     # every dataset in outputs/
    python Scripts/amplicon_16s/merge_genus_tables.py --datasets PRJNA554111 PRJNA633959 \
        --results-dir outputs --out-dir outputs/merged

Datasets processed in batches (e.g. PRJNA734525_batch1 / _batch2) are listed separately and keep
their BioProject in the `Project` column. When a dataset exists both as <name> and <name>_figaro,
only the FIGARO run is used (unless --datasets names them explicitly).

Outputs (in --out-dir):
    merged_genus_counts.csv            samples x genera (counts) + Dataset, Project, Diagnosis_provisional
    merged_genus_relative.csv          same, relative abundance (per sample, before filtering)
    merged_genus_counts_filtered.csv   after the paper's prevalence / abundance filter
    merged_genus_relative_filtered.csv
    merged_summary.csv                 samples and groups per dataset
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

LABEL = "Diagnosis_provisional"
TABLE = Path("full") / "genus" / "genus_counts_with_provisional_labels.csv"


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, default=repo / "outputs",
                    help="Folder holding <dataset>/full/genus/ (default: outputs/)")
    ap.add_argument("--datasets", nargs="*", help="Dataset names (default: every finished one)")
    ap.add_argument("--out-dir", type=Path, default=None, help="Default: <results-dir>/merged")
    ap.add_argument("--min-samples", type=int, default=2, help="Keep genera present in >= N samples (paper: 2)")
    ap.add_argument("--min-fraction", type=float, default=0.00001,
                    help="Keep genera with >= this fraction of all reads (paper: 0.001%% = 0.00001)")
    a = ap.parse_args()

    names = a.datasets or sorted(p.parent.parent.parent.name for p in a.results_dir.glob(f"*/{TABLE}"))
    if not a.datasets:
        # The same samples run twice (e.g. PRJNA496408 and PRJNA496408_figaro): keep the FIGARO run,
        # which follows the paper's method.
        for n in [n for n in names if n + "_figaro" in names]:
            print(f"skipping {n}: {n}_figaro (the paper's FIGARO method) has the same samples")
            names.remove(n)
    if not names:
        print(f"No finished datasets under {a.results_dir} (looked for */{TABLE})", file=sys.stderr)
        return 1
    frames = []
    for name in names:
        f = a.results_dir / name / TABLE
        if not f.exists():
            print(f"missing: {f}", file=sys.stderr)
            return 1
        t = pd.read_csv(f, index_col=0)
        t.index = name + ":" + t.index.astype(str)          # run IDs stay unique across datasets
        t.insert(0, "Project", name.split("_")[0])           # PRJNA734525_batch1 -> PRJNA734525
        t.insert(0, "Dataset", name)
        frames.append(t)
        print(f"{name:<22} {len(t):>4} samples  {t[LABEL].value_counts().to_dict()}")

    meta_cols = ["Dataset", "Project", LABEL]
    merged = pd.concat(frames, axis=0, sort=True)
    genera = sorted(c for c in merged.columns if c not in meta_cols)
    counts = merged[genera].fillna(0).astype("int64")
    meta = merged[meta_cols]

    present = (counts > 0).sum(axis=0)
    fraction = counts.sum(axis=0) / counts.to_numpy().sum()
    keep = (present >= a.min_samples) & (fraction >= a.min_fraction)
    if "Unassigned" in keep:
        keep["Unassigned"] = False                          # not a genus

    out = a.out_dir or a.results_dir / "merged"
    out.mkdir(parents=True, exist_ok=True)
    rel = counts.div(counts.sum(axis=1), axis=0)
    meta.join(counts).to_csv(out / "merged_genus_counts.csv")
    meta.join(rel).to_csv(out / "merged_genus_relative.csv")
    kept = counts.loc[:, keep]
    meta.join(kept).to_csv(out / "merged_genus_counts_filtered.csv")
    meta.join(kept.div(kept.sum(axis=1), axis=0)).to_csv(out / "merged_genus_relative_filtered.csv")
    summary = meta.groupby(["Dataset", LABEL]).size().unstack(fill_value=0)
    summary["total"] = summary.sum(axis=1)
    summary.to_csv(out / "merged_summary.csv")
    print(f"\n{len(counts)} samples from {len(names)} datasets | {len(genera)} genera -> {int(keep.sum())} "
          f"after filter (>= {a.min_samples} samples, >= {a.min_fraction:.5%} of reads)")
    print(f"Groups: {meta[LABEL].value_counts().to_dict()}")
    print(f"Written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
