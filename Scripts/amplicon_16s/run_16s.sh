#!/usr/bin/env bash
# Activate the QIIME 2 env, add SRA tools / FIGARO paths, then run the pipeline.
# All arguments are passed through, e.g.:
#   bash Scripts/amplicon_16s/run_16s.sh --mode dry-run --preflight
#   bash Scripts/amplicon_16s/run_16s.sh --mode full --figaro --export-dir /mnt/e/BiomeFM/outputs
set -euo pipefail

CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
QIIME_ENV="${QIIME_ENV:-rachis-qiime2-2026.7}"
TOOLS_ENV="${TOOLS_ENV:-biomefm-tools}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

source "$CONDA_BASE/etc/profile.d/conda.sh"
set +u  # some conda activation scripts reference unset variables
conda activate "$QIIME_ENV"
set -u
# Append (not prepend) so the QIIME env's python stays first on PATH.
export PATH="$PATH:$CONDA_BASE/envs/$TOOLS_ENV/bin"
export FIGARO_HOME="${FIGARO_HOME:-$HOME/opt/figaro}"
export FIGARO_PYTHON="${FIGARO_PYTHON:-$CONDA_BASE/envs/figaro/bin/python}"

if [[ " $* " != *" --config "* ]]; then
  set -- --config "$REPO/configs/PRJNA496408.yaml" "$@"
fi
exec python "$HERE/preprocess_16s_amplicon.py" "$@"
