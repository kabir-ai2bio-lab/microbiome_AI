#!/usr/bin/env bash
# End-to-end 16S workflow: preflight -> dry run (verified against expected numbers)
# -> full run -> analysis report. Every stage resumes if rerun.
#
#   bash Scripts/amplicon_16s/run_end_to_end.sh                 # everything
#   bash Scripts/amplicon_16s/run_end_to_end.sh --dry-run-only  # just verify the setup
#   EXPORT_DIR=/mnt/e/results THREADS=8 bash Scripts/amplicon_16s/run_end_to_end.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
CONFIG="${CONFIG:-$REPO/configs/PRJNA496408.yaml}"
PROJECT="$(awk '/^project:/{print $2}' "$CONFIG")"
WORKDIR="${WORKDIR:-$HOME/biomefm_runs/$PROJECT}"
EXPORT_DIR="${EXPORT_DIR:-$REPO/outputs}"
FIGARO_FLAG="--figaro"; [[ -n "${SKIP_FIGARO:-}" ]] && FIGARO_FLAG=""
EXTRA=(); [[ -n "${THREADS:-}" ]] && EXTRA+=(--threads "$THREADS")
CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"

# ${EXTRA[@]+...} keeps macOS's bash 3.2 happy with an empty array under set -u
run() { bash "$HERE/run_16s.sh" --config "$CONFIG" --workdir "$WORKDIR" --export-dir "$EXPORT_DIR" ${EXTRA[@]+"${EXTRA[@]}"} "$@"; }
analyze() {
  source "$CONDA_BASE/etc/profile.d/conda.sh"; set +u; conda activate "${QIIME_ENV:-rachis-qiime2-2026.7}"; set -u
  python "$HERE/analyze_genus.py" "$@"
}
banner() { printf '\n######## %s ########\n' "$*"; }

banner "1/5 preflight"
run --mode dry-run $FIGARO_FLAG --preflight

banner "2/5 dry run (3 samples x 5,000 read pairs)"
run --mode dry-run $FIGARO_FLAG
analyze --cohort-dir "$WORKDIR/cohorts/dry-run" --out-dir "$EXPORT_DIR/$PROJECT/dry-run/analysis" --n-permutations 99
[[ "${1:-}" == "--dry-run-only" ]] && { banner "dry run finished: $EXPORT_DIR/$PROJECT/dry-run"; exit 0; }

banner "3/5 full run (all samples, one joint DADA2)"
run --mode full $FIGARO_FLAG

banner "4/5 analysis"
analyze --cohort-dir "$WORKDIR/cohorts/full" --out-dir "$EXPORT_DIR/$PROJECT/full/analysis" ${COMPARE_TO:+--compare-to "$COMPARE_TO"}

banner "5/5 done"
echo "Results:  $EXPORT_DIR/$PROJECT/full"
echo "Report:   $EXPORT_DIR/$PROJECT/full/analysis/analysis_report.html"
echo "ML input: $EXPORT_DIR/$PROJECT/full/genus/genus_relative_with_provisional_labels.csv"
