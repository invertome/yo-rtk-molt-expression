#!/bin/bash
# Run the whole R side of the pipeline for one configuration file, on the current machine.
#
# Usage: bash run_pipeline.sh <config.yaml>
#
# The three R scripts run in order: dge_panel.R (the fits and the panel tables), supp_contrasts.R
# (the pairwise tables and the PCA), validate_outputs.R (consistency checks). The script stops at
# the first failure and prints "PIPELINE DONE" only when all three succeeded.
#
# Machine-specific settings come from a site file, ./site.env in the current directory or the file
# named by DGE_SITE_ENV (template: site.env.example). It may name the conda environment that holds
# R and the BLAS setting. Without a site file the Rscript found on PATH is used and the BLAS kernel
# is left unpinned, which the script says out loud.
#
# Environment variables that change where things are, all optional:
#   DGE_PIPELINE_DIR   directory holding the R scripts (default: the directory of this script)
#   BASE               overrides the config's base: line, so the same config can run against data
#                      kept somewhere else (paths in the config are relative to BASE)
#   PIPELINE_OUT       overrides the config's output_dir: line (through a temporary copy of the
#                      config written under BASE/tmp), to write a second set of tables without
#                      overwriting the first
# base: and output_dir: must be plain top-level lines of the YAML, as in configs/example.yaml,
# because this script reads them with sed rather than a YAML parser.
set -eo pipefail
CFG="${1:?usage: bash run_pipeline.sh <config.yaml>}"; [ -f "$CFG" ] || { echo "config not found: $CFG" >&2; exit 2; }
CFG="$(cd "$(dirname "$CFG")" && pwd)/$(basename "$CFG")"
export DGE_PIPELINE_DIR="${DGE_PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SITE="${DGE_SITE_ENV:-./site.env}"
if [ -f "$SITE" ]; then set -a; . "$SITE"; set +a; fi
# Activate the R environment if the site file names one. This happens before `set -u` because
# conda's own activation scripts read variables that may be unset and would abort under it.
if [ -n "${CONDA_ENV:-}" ]; then
  sh="${CONDA_SH:?CONDA_SH must name conda.sh when CONDA_ENV is set}"; sh="${sh/#\~/$HOME}"; sh="${sh/#\$HOME/$HOME}"
  . "$sh" && conda activate "$CONDA_ENV"
fi
set -u
# BLAS setting (see assert_blas_pin in lib_dge.R): either OPENBLAS_CORETYPE pins one kernel, or
# DGE_BLAS_PIN=none opts out on purpose. With a site file that sets neither, stop. With no site
# file at all, run unpinned and say so.
if [ -n "${OPENBLAS_CORETYPE:-}" ]; then export OPENBLAS_CORETYPE
elif [ "${DGE_BLAS_PIN:-}" = none ]; then export DGE_BLAS_PIN
elif [ ! -f "$SITE" ]; then export DGE_BLAS_PIN=none; echo "BLAS kernel not pinned (no site file); set OPENBLAS_CORETYPE in site.env to reproduce a pinned run"
else echo "set OPENBLAS_CORETYPE or DGE_BLAS_PIN=none in $SITE" >&2; exit 2; fi
# Read base: and name: from the config (quotes stripped), apply the BASE override, and if
# PIPELINE_OUT is set write a temporary config whose output_dir: points there.
base="$(sed -n 's/^base:[[:space:]]*//p' "$CFG" | head -1 | tr -d "'\"")"; name="$(sed -n 's/^name:[[:space:]]*//p' "$CFG" | head -1 | tr -d "'\"")"
[ -n "$base" ] || { echo "config has no top-level base: line: $CFG" >&2; exit 2; }
export BASE="${BASE:-$base}"
if [ -n "${PIPELINE_OUT:-}" ]; then
  mkdir -p "$BASE/tmp"; tmpcfg="$BASE/tmp/${name:-config}_$$.yaml"
  sed "s#^output_dir:.*#output_dir: $PIPELINE_OUT#" "$CFG" > "$tmpcfg"; CFG="$tmpcfg"
fi
# The R scripts run from BASE so that every relative path in the config resolves there.
cd "$BASE"
for s in dge_panel.R supp_contrasts.R validate_outputs.R; do echo "=== [$(date)] $s ==="; Rscript "$DGE_PIPELINE_DIR/$s" "$CFG"; done
echo "[$(date)] PIPELINE DONE"
