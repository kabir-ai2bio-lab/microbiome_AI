#!/usr/bin/env bash
# One-time setup for the 16S pipeline on Linux, macOS (Intel or Apple Silicon) or
# Windows WSL. Installs Miniforge, QIIME 2 2026.7 (amplicon plugins), SRA Toolkit +
# pigz, and FIGARO. Safe to re-run: anything already installed is skipped.
#
#   bash envs/setup_local.sh
#   SKIP_FIGARO=1 bash envs/setup_local.sh       # without FIGARO
set -euo pipefail

CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
QIIME_ENV="${QIIME_ENV:-rachis-qiime2-2026.7}"
TOOLS_ENV="${TOOLS_ENV:-biomefm-tools}"
FIGARO_HOME="${FIGARO_HOME:-$HOME/opt/figaro}"
FIGARO_ENV="${FIGARO_ENV:-figaro}"
FIGARO_REV="0916f0d593eb6ee8e7fd210d7e48171882f7c219"
QIIME_URL="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released"

say() { printf '\n==> %s\n' "$*"; }

OS="$(uname -s)"; ARCH="$(uname -m)"
case "$OS" in
  Linux)  QIIME_YML="$QIIME_URL/rachis-qiime2-linux-64-conda.yml"
          MINIFORGE="Miniforge3-Linux-x86_64.sh" ;;
  Darwin) QIIME_YML="$QIIME_URL/rachis-qiime2-osx-64-conda.yml"
          MINIFORGE="Miniforge3-MacOSX-${ARCH}.sh"
          # QIIME 2 only ships Intel macOS builds; Apple Silicon runs them under Rosetta 2.
          if [[ "$ARCH" == "arm64" ]]; then
            export CONDA_SUBDIR=osx-64
            /usr/bin/pgrep -q oahd || softwareupdate --install-rosetta --agree-to-license
          fi ;;
  *) echo "Unsupported OS: $OS (on Windows, run this inside WSL)"; exit 1 ;;
esac
sed_inplace() { if [[ "$OS" == Darwin ]]; then sed -i '' "$@"; else sed -i "$@"; fi; }

if [[ "$OS" == Linux ]] && command -v apt-get >/dev/null; then
  say "System packages (git, wget, build tools)"
  SUDO=""; [[ $EUID -ne 0 ]] && SUDO="sudo"
  $SUDO apt-get update -qq
  DEBIAN_FRONTEND=noninteractive $SUDO apt-get install -y -qq git wget curl ca-certificates build-essential tmux
fi

if [[ ! -x "$CONDA_BASE/bin/conda" ]]; then
  say "Installing Miniforge to $CONDA_BASE"
  curl -fsSL -o /tmp/miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/$MINIFORGE"
  bash /tmp/miniforge.sh -b -p "$CONDA_BASE"
  rm -f /tmp/miniforge.sh
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

if ! env_exists "$QIIME_ENV"; then
  say "Creating QIIME 2 env $QIIME_ENV (10-30 min)"
  conda env create --name "$QIIME_ENV" --file "$QIIME_YML" --yes
  [[ -n "${CONDA_SUBDIR:-}" ]] && conda run -n "$QIIME_ENV" conda config --env --set subdir osx-64
fi

if ! env_exists "$TOOLS_ENV"; then
  say "Creating $TOOLS_ENV env (sra-tools, pigz)"
  conda create -y -n "$TOOLS_ENV" -c conda-forge -c bioconda sra-tools pigz
fi

if [[ -z "${SKIP_FIGARO:-}" ]]; then
  say "FIGARO at $FIGARO_HOME (rev ${FIGARO_REV:0:7})"
  if [[ ! -d "$FIGARO_HOME/.git" ]]; then
    mkdir -p "$(dirname "$FIGARO_HOME")"
    git clone -q https://github.com/Zymo-Research/figaro.git "$FIGARO_HOME"
  fi
  git -C "$FIGARO_HOME" checkout -q "$FIGARO_REV"
  # Python 3.10 removed collections.Iterable; FIGARO still uses it.
  compat="$FIGARO_HOME/figaro/easyMultiprocessing.py"
  sed_inplace 's/collections\.Iterable/collections.abc.Iterable/g' "$compat"
  grep -q '^import collections.abc' "$compat" || { printf 'import collections.abc\n' | cat - "$compat" > "$compat.tmp" && mv "$compat.tmp" "$compat"; }
  # FIGARO's requirements.txt is only numpy/scipy/matplotlib; take them from conda-forge
  # rather than pip (pip's PyPI lookups fail on some WSL1 DNS setups).
  if ! env_exists "$FIGARO_ENV"; then
    conda create -y -n "$FIGARO_ENV" -c conda-forge python=3.10 "numpy>=1.13.1" "scipy>=1.2.1" "matplotlib>=3.0.2"
  else
    conda install -y -n "$FIGARO_ENV" -c conda-forge "numpy>=1.13.1" "scipy>=1.2.1" "matplotlib>=3.0.2"
  fi
fi

say "Verifying"
set +u
conda activate "$QIIME_ENV"
set -u
export PATH="$PATH:$CONDA_BASE/envs/$TOOLS_ENV/bin"
qiime --version | head -n 1
prefetch --version | grep -m1 -o 'prefetch : [0-9.]*' || prefetch --version | head -n 2
pigz --version
python -c "import os; print('CPUs visible:', os.cpu_count())"
say "Setup complete. Next: bash Scripts/amplicon_16s/run_16s.sh --mode dry-run --preflight"
