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

# --- Lock guard: serialize env mutations (shared with setup_env_sbatch.sh) ---
# Both entry points write into .pixi/envs/py-cuda13 (conda + R library). The
# lock is an atomic directory containing PID + timestamp metadata.  It is
# fail-closed after interruption; remove it only after verifying no writer is
# active, rather than risking concurrent prefix mutation.
ENV_LOCK_FILE="${LOGS_DIR}/env_refresh.lock"
source "${SCRIPT_DIR}/env_mutation_lock.sh"

# --- Guard: no active jobs while mutating the env ---------------------------
acquire_env_lock
ACTIVE_JOBS="$(squeue -u "${USER}" -h -o "%j" 2>/dev/null || true)"
if [[ -n "${ACTIVE_JOBS}" ]]; then
  echo "ERROR: active Slurm jobs detected — env refresh must run while no jobs are running" >&2
  echo "       (concurrent installs corrupt the shared R library). Active jobs:" >&2
  squeue -u "${USER}" >&2
  exit 1
fi

# --- Failure diagnostics: printed only when pixi run setup fails ---------------
# Identical output contract as setup_env_sbatch.sh's run_env_diagnostics, so
# login-node and worker-node failures are consistent. Single-quoted shell block
# -> R uses DOUBLE quotes here (opposite constraint from the pixi.toml
# [tasks.setup] TOML block).
run_env_diagnostics() {
  "${PIXI_BIN}" run -e py-cuda13 Rscript --vanilla -e '
    cat("R version:", R.version.string, "\n")
    cat("Library paths:\n")
    for (lib in .libPaths()) cat("  ", lib, "\n")
    for (lib in .libPaths()) {
      dirs <- list.dirs(lib, recursive = FALSE)
      has_desc <- vapply(dirs, function(d) file.exists(file.path(d, "DESCRIPTION")), logical(1))
      dirs <- dirs[has_desc]
      cat("Library ", lib, ": ", length(dirs), " dirs with DESCRIPTION\n", sep = "")
      for (p in dirs) {
        if (!file.exists(file.path(p, "Meta", "package.rds"))) {
          pkg_line <- grep("^Package:", readLines(file.path(p, "DESCRIPTION"), n = 20), value = TRUE)
          cat("Missing Meta/package.rds:", p, paste(pkg_line, collapse = "; "), "\n")
        }
      }
    }
    cat("Hint: the only legitimate case is .../translations (R message catalogs, false positive -> update pixi.toml); anything else means the env is genuinely corrupt -> wipe-and-reinstall (rm -rf .pixi/envs/py-cuda13 && pixi install -e py-cuda13 && pixi run -e py-cuda13 setup).\n")
  ' || true
}

echo "=== [1/3] pixi install -e py-cuda13 (conda + pypi deps) ==="
"${PIXI_BIN}" install --environment py-cuda13

echo "=== [2/3] pixi run -e py-cuda13 setup (R source packages + integrity check) ==="
if ! "${PIXI_BIN}" run -e py-cuda13 setup; then
  echo "ERROR: pixi run setup failed — running env diagnostics..." >&2
  run_env_diagnostics
  exit 1
fi

echo "=== [3/3] smoke check: notebook + worker import lists + env_check.R ==="
# Sources src/utils/imports.R (the notebook loader attach list), the slim
# worker subsets (imports_worker_core.R / imports_worker_transzeroimp.R,
# attached by the R benchmark / transzeroimp workers) and the repo-wide
# src/utils/env_check.R (requireNamespace check: attach ∪ namespaced-only ∪
# worker/annotation packages — covers MOFA2/scITD/GloScope/annotation extras
# no longer attached by imports.R); a missing or broken package stops here
# with the same error a worker would produce (set -euo pipefail turns it into
# a FAILED refresh). Run from the project root (done above), so the relative
# source paths resolve.
"${PIXI_BIN}" run -e py-cuda13 Rscript --vanilla -e '
  source("src/utils/imports.R")
  source("src/utils/imports_worker_core.R")
  source("src/utils/imports_worker_transzeroimp.R")
  source("src/utils/env_check.R")
  cat("All packages in src/utils/imports.R + worker subsets + env_check.R load OK\n")
'

echo "Environment refresh complete. You can now submit jobs."
