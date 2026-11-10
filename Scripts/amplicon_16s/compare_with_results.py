#!/usr/bin/env python3
"""
Check that your runs reproduce the committed results: compares every dataset you finished in
outputs/ with the same dataset in results/ (and outputs/merged with results/merged).

    python Scripts/amplicon_16s/compare_with_results.py                 # every dataset in outputs/
    python Scripts/amplicon_16s/compare_with_results.py PRJNA554111     # one dataset

For each dataset it compares the samples, labels, ASV count, DADA2 read counts per sample and the
genus count table value by value. The pipeline is deterministic (DADA2, VSEARCH and the analysis use
fixed seeds; the thread count does not matter), so a correct installation gives IDENTICAL.
Exit code 0 = everything compared is identical.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

LABEL = "Diagnosis_provisional"


def table(path: Path) -> pd.DataFrame:
    t = pd.read_csv(path, index_col=0)
    return t.sort_index().sort_index(axis=1)


def compare_dataset(mine: Path, ref: Path) -> list:
    """Return a list of differences (empty = identical)."""
    diffs = []
    a = table(mine / "genus" / "genus_counts_with_provisional_labels.csv")
    b = table(ref / "genus" / "genus_counts_with_provisional_labels.csv")
    if list(a.index) != list(b.index):
        diffs.append(f"samples differ: {len(a)} here vs {len(b)} in results/")
    if not a[LABEL].reindex(b.index).equals(b[LABEL]):
        diffs.append("diagnosis labels differ")
    if set(a.columns) != set(b.columns):
        diffs.append(f"genera differ: {len(set(a.columns) - set(b.columns))} only here, "
                     f"{len(set(b.columns) - set(a.columns))} only in results/")
    common = a.columns.intersection(b.columns).drop(LABEL)
    idx = a.index.intersection(b.index)
    n_diff = int((a.loc[idx, common] != b.loc[idx, common]).to_numpy().sum())
    if n_diff:
        diffs.append(f"{n_diff} genus count values differ")
    na = len(pd.read_csv(mine / "asv_taxonomy_consensus.csv"))
    nb = len(pd.read_csv(ref / "asv_taxonomy_consensus.csv"))
    if na != nb:
        diffs.append(f"ASVs: {na} here vs {nb} in results/")
    sa = pd.read_csv(mine / "dada2_per_sample_stats.csv", index_col=0)
    sb = pd.read_csv(ref / "dada2_per_sample_stats.csv", index_col=0)
    col = "non-chimeric"
    if col in sa and col in sb and not sa[col].sort_index().equals(sb[col].sort_index()):
        diffs.append("DADA2 non-chimeric reads per sample differ")
    return diffs


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("datasets", nargs="*", help="Dataset names (default: every finished one in --outputs)")
    ap.add_argument("--outputs", type=Path, default=repo / "outputs", help="Your results (default: outputs/)")
    ap.add_argument("--results", type=Path, default=repo / "results", help="Committed results (default: results/)")
    a = ap.parse_args()

    names = a.datasets or sorted(p.parent.name for p in a.outputs.glob("*/full") if (p / "genus").is_dir())
    if not names and not (a.outputs / "merged").is_dir():
        print(f"Nothing to compare: no finished dataset in {a.outputs}", file=sys.stderr)
        return 1
    failed = 0
    for name in names:
        mine, ref = a.outputs / name / "full", a.results / name / "full"
        if not ref.is_dir():
            print(f"{name:<22} no committed results to compare with ({ref})")
            continue
        if not (mine / "genus").is_dir():
            print(f"{name:<22} not finished in {mine}")
            failed += 1
            continue
        diffs = compare_dataset(mine, ref)
        print(f"{name:<22} {'IDENTICAL' if not diffs else 'DIFFERENT: ' + '; '.join(diffs)}")
        failed += bool(diffs)

    if not a.datasets and (a.outputs / "merged").is_dir() and (a.results / "merged").is_dir():
        f = "merged_genus_counts_filtered.csv"
        same = table(a.outputs / "merged" / f).equals(table(a.results / "merged" / f))
        print(f"{'merged':<22} {'IDENTICAL' if same else 'DIFFERENT (expected if not all 14 datasets were run)'}")
        failed += not same
    print("\nAll compared results are identical." if not failed else
          "\nSome results differ: check that QIIME 2 2026.7 is installed and the dry-run check passed.")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
