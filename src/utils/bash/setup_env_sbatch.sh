#!/bin/bash
#SBATCH --job-name=setup_env
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=shared-cpu
#SBATCH --mail-type=END,FAIL
# Logs: hardcoded to the current HPC repo location (single-user repo). SLURM
# does NOT expand ${HOME} inside #SBATCH lines; the script mkdir -p's LOGS_DIR.
#SBATCH --output=[REDACTED_PATH]/ECODA_paper/logs/setup_env_%j.out
#SBATCH --error=[REDACTED_PATH]/ECODA_paper/logs/setup_env_%j.err
# ============================================================
# setup_env_sbatch.sh — full py-cuda13 env build on a worker node
# ============================================================
# Submit from the HPC login node (defaults: 16 cpus / 64G / 2h, overridable):
#   sbatch src/utils/bash/setup_env_sbatch.sh
#   sbatch --cpus-per-task=32 --mem=128G src/utils/bash/setup_env_sbatch.sh
#
# Worker nodes share $HOME/scratch and have internet; R source packages are
# compiled with parallel make (MAKEFLAGS=-jN via R_SETUP_JOBS) so the build is
# much faster than on the 2-core login node and survives SSH drops.
#
# Serialization contract with src/utils/bash/refresh_env.sh (login-node path):
# both acquire the shared atomic directory lock
# ${LOGS_DIR}/env_refresh.lock (PID + timestamp).  The lock fails closed after
# interruption; remove it only after verifying no writer is active.  Both
# entry points also refuse to start while other jobs are active — a mutation
# racing with running array tasks corrupts the shared R library (observed:
# digest/Meta/package.rds missing, mime lazyload DB missing).
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch the script is copied to /var/spool/slurmd/job<id>/slurm_script,
# so BASH_SOURCE no longer points at the repo. Recover the real path from the
# job record (Command= field); BASH_SOURCE fallback for login-node execution.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
cd "${PROJECT_ROOT}"
mkdir -p "${LOGS_DIR}"

PIXI_BIN="${HOME}/.pixi/bin/pixi"

# [tasks.setup] reads this for MAKEFLAGS=-jN (parallel R source compilation).
export R_SETUP_JOBS="${SLURM_CPUS_PER_TASK:-16}"

# --- Lock guard: serialize env mutations (shared with refresh_env.sh) --------
# Both entry points write into .pixi/envs/py-cuda13 (conda + R library). The
# lock is an atomic directory containing PID + timestamp metadata.  It is
# fail-closed after interruption; remove it only after verifying no writer is
# active, rather than risking concurrent prefix mutation.
ENV_LOCK_FILE="${LOGS_DIR}/env_refresh.lock"
source "${SCRIPT_DIR}/env_mutation_lock.sh"

# --- Guard: no other active jobs while mutating the env ----------------------
check_no_active_jobs() {
  ecoda_require_no_active_jobs "${SLURM_JOB_ID:-}"
}

# --- Toolchain preflight (conda compilers, GCCcore module fallback) ----------
check_toolchain() {
  if "${PIXI_BIN}" run -e py-cuda13 bash -c 'command -v gcc >/dev/null 2>&1 && command -v make >/dev/null 2>&1'; then
    echo "OK: gcc and make available in the py-cuda13 env (conda compilers)."
    return 0
  fi
  echo "WARNING: gcc/make not found in the py-cuda13 env — falling back to the GCCcore module..." >&2
  module load GCCcore/12.2.0 >/dev/null 2>&1 || true
  if command -v gcc >/dev/null 2>&1 && command -v make >/dev/null 2>&1; then
    echo "OK: gcc and make found via the GCCcore/12.2.0 module."
    return 0
  fi
  echo "ERROR: no C/C++ toolchain available for R source compilation." >&2
  echo "       Run src/utils/bash/refresh_env.sh once on the login node, or verify" >&2
  echo "       the compilers/make pixi deps were installed into the py-cuda13 env." >&2
  return 1
}

# --- Failure diagnostics: printed only when pixi run setup fails ---------------
# Prints R version, .libPaths(), per-library DESCRIPTION-dir counts and every
# dir missing Meta/package.rds (full path + Package: field), plus the
# translations hint. NOTE: single-quoted shell block -> R uses DOUBLE quotes
# here (opposite constraint from the pixi.toml [tasks.setup] TOML block).
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

echo "=== [0/4] guards (env lock + no other active jobs) ==="
acquire_env_lock
check_no_active_jobs
echo "OK: env lock acquired, no other jobs active."

start="$(date +%s)"

echo "=== [1/4] pixi install -e py-cuda13 (conda + pypi deps) ==="
"${PIXI_BIN}" install --environment py-cuda13
echo "pixi install done: $(( $(date +%s) - start ))s elapsed."

echo "=== [2/4] toolchain preflight (gcc/make for R source compilation) ==="
check_toolchain

echo "=== [3/4] pixi run -e py-cuda13 setup (R source packages, R_SETUP_JOBS=${R_SETUP_JOBS} parallel make) ==="
if ! "${PIXI_BIN}" run -e py-cuda13 setup; then
  echo "ERROR: pixi run setup failed — running env diagnostics..." >&2
  run_env_diagnostics
  exit 1
fi
echo "setup done: $(( $(date +%s) - start ))s elapsed."

echo "=== [4/4] smoke check: notebook + worker import lists + env_check.R ==="
# Sources src/utils/imports.R (the notebook loader attach list), the slim
# worker subsets (imports_worker_core.R / imports_worker_transzeroimp.R,
# attached by the R benchmark / transzeroimp workers) and the repo-wide
# src/utils/env_check.R (requireNamespace check: attach ∪ namespaced-only ∪
# worker/annotation packages — covers MOFA2/scITD/GloScope/annotation extras
# no longer attached by imports.R); a missing or broken package stops here
# with the same error a worker would produce (set -euo pipefail turns it into
# a FAILED build). Run from the project root (done above), so the relative
# source paths resolve.
"${PIXI_BIN}" run -e py-cuda13 Rscript --vanilla -e '
  source("src/utils/imports.R")
  source("src/utils/imports_worker_core.R")
  source("src/utils/imports_worker_transzeroimp.R")
  source("src/utils/env_check.R")
  cat("All packages in src/utils/imports.R + worker subsets + env_check.R load OK\n")
'

echo "Environment build complete ($(( $(date +%s) - start ))s total). You can now submit jobs."
