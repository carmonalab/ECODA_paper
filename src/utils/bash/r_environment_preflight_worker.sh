#!/bin/bash
# Read-only R package smoke check on one allocated compute node.
set -euo pipefail

SCRIPT_DIR=""
if [[ -n "${PROJECT_ROOT:-}" &&
      -f "${PROJECT_ROOT}/src/slurm_config.sh" ]]; then
  SCRIPT_DIR="${PROJECT_ROOT}/src/utils/bash"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" &&
        -f "${SLURM_SUBMIT_DIR}/src/slurm_config.sh" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}/src/utils/bash"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  if [[ -n "${SLURM_JOB_ID:-}" ]] &&
     command -v scontrol >/dev/null 2>&1; then
    submitted_command="$(scontrol show job "${SLURM_JOB_ID}" -o 2>/dev/null |
      sed -n 's/.* Command=\([^ ]*\).*/\1/p' | head -1 || true)"
    submitted_dir="$(dirname "${submitted_command}")"
    if [[ -n "${submitted_command}" &&
          -f "${submitted_dir}/../../slurm_config.sh" ]]; then
      SCRIPT_DIR="$(cd "${submitted_dir}" && pwd)"
    fi
  fi
fi
if [[ -z "${SCRIPT_DIR}" || ! -f "${SCRIPT_DIR}/../../slurm_config.sh" ]]; then
  echo "ERROR: could not recover the repository source directory." >&2
  exit 1
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
ecoda_runtime_reexec_worker "${ECODA_RUNTIME_PROFILE:-stage5}" \
  "${SCRIPT_DIR}/r_environment_preflight_worker.sh" || exit 1
cd "${PROJECT_ROOT}"
R_COMMAND="${R_ENV_PREFLIGHT_RSCRIPT:-${PIXI_RSCRIPT:-}}"
[[ -n "${R_COMMAND}" ]] || {
  echo "ERROR: configured PIXI_RSCRIPT is unavailable." >&2
  exit 1
}

${R_COMMAND} -e '
  cat("node=", Sys.info()[["nodename"]], "\n", sep="")
  cat("R=", R.version.string, "\n", sep="")
  cat("HOME=", Sys.getenv("HOME"), "\n", sep="")
  cat("libPaths=", paste(.libPaths(), collapse=";"), "\n", sep="")
  packages <- c("arrow", "DelayedArray", "DESeq2", "EPIC", "GloScope", "MOFA2", "scITD")
  for (pkg in packages) {
    location <- find.package(pkg, quiet=TRUE)
    rdb <- if (length(location)) file.path(location, "R", paste0(pkg, ".rdb")) else ""
    cat("package=", pkg, " path=", if (length(location)) location else "<missing>",
        " rdb=", if (nzchar(rdb)) rdb else "<missing>",
        " exists=", if (nzchar(rdb)) file.exists(rdb) else FALSE, "\n", sep="")
  }
  suppressPackageStartupMessages({
    library(arrow)
    library(DESeq2)
    library(EPIC)
    library(GloScope)
    library(MOFA2)
    library(scITD)
  })
  source("src/utils/imports_worker_core.R")
  cat("benchmark R environment preflight OK\n")
'
