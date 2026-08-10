#!/bin/bash
# ============================================================
# refresh_env.sh — serial, guarded refresh of the HPC py-cuda13 pixi env
# ============================================================
# Run on the HPC LOGIN NODE only. Install on the login node is allowed, but
# never while jobs are running:
#   - conda (pixi install) and remotes (pixi run setup) both write into
#     .pixi/envs/py-cuda13/lib/R/library; a mutation racing with running array
#     tasks corrupts the shared library (observed: digest/Meta/package.rds
#     missing -> 'there is no package called digest'; mime lazyload DB missing).
#   - `pixi install` does NOT verify installed package files, so a corrupt
#     package dir survives re-installs; the setup task now runs an integrity
#     check and fails loudly with the wipe-and-reinstall repair.
#   - an SSH drop can kill a running install mid-write -> run inside tmux.
#
# Usage (HPC login node):
#   tmux new -s env-refresh
#   src/utils/bash/refresh_env.sh
#
# This script refuses to start while ANY of your Slurm jobs are active.
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
cd "${PROJECT_ROOT}"

PIXI_BIN="${HOME}/.pixi/bin/pixi"

# --- Guard: no active jobs while mutating the env ---------------------------
ACTIVE_JOBS="$(squeue -u "${USER}" -h -o "%j" 2>/dev/null || true)"
if [[ -n "${ACTIVE_JOBS}" ]]; then
  echo "ERROR: active Slurm jobs detected — env refresh must run while no jobs are running" >&2
  echo "       (concurrent installs corrupt the shared R library). Active jobs:" >&2
  squeue -u "${USER}" >&2
  exit 1
fi

echo "=== [1/3] pixi install -e py-cuda13 (conda + pypi deps) ==="
"${PIXI_BIN}" install --environment py-cuda13

echo "=== [2/3] pixi run -e py-cuda13 setup (R source packages + integrity check) ==="
"${PIXI_BIN}" run -e py-cuda13 setup

echo "=== [3/3] smoke check: critical packages load ==="
"${PIXI_BIN}" run -e py-cuda13 Rscript --vanilla -e '
  invisible(lapply(c("digest", "SeuratObject", "Seurat", "anndataR"), library, character.only = TRUE))
  cat("All critical packages load OK\n")
'

echo "Environment refresh complete. You can now submit jobs."
