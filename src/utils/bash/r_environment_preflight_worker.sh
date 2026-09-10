#!/bin/bash
# Read-only R package smoke check on one allocated compute node.
set -euo pipefail

source_root="${ECODA_SOURCE_ROOT:-}"
source_manifest="${ECODA_SOURCE_MANIFEST:-}"
source_required="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
run_root="${R_ENV_PREFLIGHT_RUN_ROOT:-${ECODA_RUN_ROOT:-}}"
run_id="${R_ENV_PREFLIGHT_RUN_ID:-${ECODA_RUN_ID:-${run_root##*/}}}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${run_root}" ]] || {
  echo "ERROR: R environment preflight run root binding mismatch" >&2
  exit 1
}
[[ -z "${ECODA_RUN_ID:-}" || "${ECODA_RUN_ID}" == "${run_id}" ]] || {
  echo "ERROR: R environment preflight run ID binding mismatch" >&2
  exit 1
}
[[ -z "${R_ENV_PREFLIGHT_RUN_ID:-}" || "${R_ENV_PREFLIGHT_RUN_ID}" == "${run_id}" ]] || {
  echo "ERROR: R environment preflight scheduler run ID binding mismatch" >&2
  exit 1
}

[[ "${source_required}" == "1" ]] || {
  echo "ERROR: legacy_source_unpinned: R environment preflight requires an immutable source snapshot." >&2
  exit 1
}
[[ "${source_root}" = /* && "${source_manifest}" = /* ]] || {
  echo "ERROR: R environment preflight source root and manifest are required." >&2
  exit 1
}
[[ "${run_root}" = /* ]] || {
  echo "ERROR: R environment preflight run root is required." >&2
  exit 1
}
[[ "${run_id}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ &&
   "${run_root##*/}" == "${run_id}" ]] || {
  echo "ERROR: invalid R environment preflight run identity" >&2
  exit 1
}
[[ -f "${source_root}/src/slurm_config.sh" &&
   -r "${source_root}/src/slurm_config.sh" ]] || {
  echo "ERROR: immutable R environment preflight source root is incomplete: ${source_root}" >&2
  exit 1
}

SCRIPT_DIR="${source_root%/}/src/utils/bash"
source "${source_root}/src/slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
[[ "${ECODA_RUNS_ROOT:-}" = /* &&
   "${run_root}" == "${ECODA_RUNS_ROOT%/}/${run_id}" ]] || {
  echo "ERROR: R environment preflight run root is not the global run root for its ID" >&2
  exit 1
}
worker_script="${source_root%/}/src/utils/bash/r_environment_preflight_worker.sh"
worker_script="$(ecoda_require_source_script_path "${worker_script}" "${source_root}")" ||
  exit 1

export ECODA_RUN_ROOT="${run_root}" ECODA_RUN_ID="${run_id}"
export ECODA_SOURCE_ROOT="${source_root}"
export ECODA_SOURCE_MANIFEST="${source_manifest}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED="${source_required}"
export ECODA_RUNTIME_IMAGE="${ECODA_RUNTIME_IMAGE:-}"
export ECODA_RUNTIME_MANIFEST="${ECODA_RUNTIME_MANIFEST:-}"
export ECODA_HOST_ENV_PREFIX="${ECODA_HOST_ENV_PREFIX:-}"
export ECODA_SCRATCH_ROOT="${ECODA_SCRATCH_ROOT:-${HPC_SCRATCH_DIR:-}}"
export ECODA_LOGS_DIR="${ECODA_LOGS_DIR:-${LOGS_DIR:-}}"
export ECODA_AUX_ROOT="${ECODA_AUX_ROOT:-${source_root%/}/aux}"
export R_ENV_PREFLIGHT_RUN_ROOT="${run_root}" R_ENV_PREFLIGHT_RUN_ID="${run_id}"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" &&
      "${ECODA_RUNTIME_MODE:-host}" == "host" ]]; then
  ecoda_runtime_validate_bound_run || exit 1
fi
ecoda_runtime_reexec_worker "${ECODA_RUNTIME_PROFILE:-stage5}" \
  "${worker_script}" || exit 1
cd "${PROJECT_ROOT}"
if [[ -n "${R_ENV_PREFLIGHT_RSCRIPT:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" &&
      "${R_ENV_PREFLIGHT_RSCRIPT}" != "${PIXI_RSCRIPT}" ]]; then
  echo "ERROR: container R environment preflight cannot override its in-image Rscript." >&2
  exit 1
fi
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
