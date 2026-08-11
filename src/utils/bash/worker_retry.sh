#!/bin/bash
#
# Worker-node self-healing helpers: transient-failure requeue + R library
# staging + thread pinning.
#
# Source AFTER slurm_config.sh, from an array worker running on a compute
# node (uses PROJECT_ROOT, HPC_SCRATCH_DIR, LOGS_DIR):
#   source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"
#   (adjust relative path based on script depth: annotation/preprocess workers
#   sit one level below src/, benchmark workers two levels).
#
# Provided:
#   TRANSIENT_REQEX
#       grep -Ei pattern of transient failures (stale BeeGFS client-cache
#       views, env-layout races, missing files/imports). Anything matching
#       makes a task eligible for self-requeue.
#   worker_retry_count_file
#       Prints the per-task retry counter path
#       ${HPC_SCRATCH_DIR}/_worker_retries/<jobid>_<taskid>.count (jobid =
#       SLURM_ARRAY_JOB_ID if set, else SLURM_JOB_ID; taskid = array task id
#       or 0). Filenames embed job+task ids, so counters never collide across
#       submissions and stale files (task killed after success) are harmless.
#   worker_bump_retry_count <max>
#       Reads the counter (0 if missing); returns 1 when count >= max, else
#       writes count+1 atomically (tmp + mv) and echoes the new count.
#   worker_clear_retry_count
#       rm -f the counter file (called on the success path).
#   worker_requeue_if_transient <err_file> [max]
#       Greps <err_file> (the Slurm --error file for this task; must exist and
#       be non-empty) with TRANSIENT_REQEX. On match, bumps the retry counter
#       (default max ${WORKER_MAX_RETRIES:-3}); if the cap is not reached:
#       echoes "Transient failure detected; requeueing (attempt N)", runs
#       `scontrol requeue <array_job>_<task>` (or <job> for plain jobs),
#       sleeps 2 s, returns 0 (the caller exits 0; the task restarts, likely
#       on another node). Otherwise returns 1 (the caller propagates the real
#       exit code, failing the task fail-closed).
#   stage_env_rlib [prefix]
#       Copies the pixi R library (${PROJECT_ROOT}/.pixi/envs/py-cuda13/lib/R/
#       library) to a per-task dir under node-local storage so R package
#       loading is immune to stale BeeGFS client-cache views. Guard: skipped
#       (warn, return 0) when WORKER_STAGE_R_LIB != 1 (default 1) or the
#       library exceeds WORKER_R_LIB_MAX_MB (default 10240; measured with
#       du -sm — the size guard keeps metadata load off BeeGFS). Target root:
#       ${WORKER_STAGE_ROOT:-/scratch} (falls back to /tmp when /scratch is
#       absent), then ${STAGE_ROOT}/${prefix}_R_library under the per-task dir
#       <jobid>_<taskid>. On success exports R_LIBS with the staged copy
#       FIRST (default env library remains the fallback). On copy failure
#       returns 1 — the caller's transient grep covers stale-view "No such
#       file or directory" errors. Node-shared /srv/share/users/... staging is
#       a documented follow-up; the staging root is pluggable via
#       WORKER_STAGE_ROOT.
#   export_worker_thread_env
#       Pins BLAS/OMP/MKL/NUMEXPR/veclib to 1 thread so CPU time ~= wall time
#       (annotation worker only — benchmark workers stay hardware-pinned for
#       runtime comparability, preprocess is intentionally multi-threaded).
#
# All functions are safe under `set -euo pipefail` when invoked in the
# standard conditional patterns (function calls in `if` conditions disable
# errexit inside the body); every fallible command is guarded or explicit.

TRANSIENT_REQEX='No such file or directory|cannot open file|cannot open shared object|package or namespace load failed|there is no package called|cannot open connection|failed to load|No module named'

worker_retry_count_file() {
  local JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-unknown}}"
  local TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
  printf '%s' "${HPC_SCRATCH_DIR:-/tmp}/_worker_retries/${JOB_ID}_${TASK_ID}.count"
}

worker_bump_retry_count() {
  local MAX="${1:-${WORKER_MAX_RETRIES:-3}}"
  local FILE
  FILE="$(worker_retry_count_file)"
  local COUNT=0
  if [[ -f "${FILE}" ]]; then
    COUNT="$(cat "${FILE}" 2>/dev/null || true)"
    COUNT="${COUNT:-0}"
    [[ "${COUNT}" =~ ^[0-9]+$ ]] || COUNT=0
  fi
  if (( COUNT >= MAX )); then
    return 1
  fi
  local NEW=$(( COUNT + 1 ))
  mkdir -p "$(dirname "${FILE}")"
  local TMP="${FILE}.tmp.$$"
  printf '%s\n' "${NEW}" > "${TMP}"
  mv -f "${TMP}" "${FILE}"
  printf '%s' "${NEW}"
  return 0
}

worker_clear_retry_count() {
  rm -f "$(worker_retry_count_file)"
}

worker_requeue_if_transient() {
  local ERR_FILE="${1:-}"
  local MAX="${2:-${WORKER_MAX_RETRIES:-3}}"
  if [[ -z "${ERR_FILE}" || ! -s "${ERR_FILE}" ]]; then
    echo "Task ${SLURM_ARRAY_TASK_ID:-?}: no non-empty err file (${ERR_FILE:-unset}); not requeueing."
    return 1
  fi
  if ! grep -Eiq "${TRANSIENT_REQEX}" "${ERR_FILE}"; then
    echo "Task ${SLURM_ARRAY_TASK_ID:-?}: no transient signature in ${ERR_FILE}; not requeueing."
    return 1
  fi
  local N=""
  if ! N="$(worker_bump_retry_count "${MAX}")"; then
    echo "Task ${SLURM_ARRAY_TASK_ID:-?}: transient failure detected but retry limit (${MAX}) reached; not requeueing."
    return 1
  fi
  echo "Transient failure detected; requeueing (attempt ${N})"
  if [[ -n "${SLURM_ARRAY_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    scontrol requeue "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
  else
    scontrol requeue "${SLURM_JOB_ID}"
  fi
  sleep 2
  return 0
}

stage_env_rlib() {
  local PREFIX="${1:-env}"
  if [[ "${WORKER_STAGE_R_LIB:-1}" != "1" ]]; then
    echo "WORKER_STAGE_R_LIB=${WORKER_STAGE_R_LIB:-1}; skipping R library staging (BeeGFS fallback; retry covers flakes)."
    return 0
  fi
  local SRC="${PROJECT_ROOT}/.pixi/envs/py-cuda13/lib/R/library"
  if [[ ! -d "${SRC}" ]]; then
    echo "WARNING: R library source ${SRC} not found; skipping staging (BeeGFS fallback)."
    return 0
  fi
  local MAX_MB="${WORKER_R_LIB_MAX_MB:-10240}"
  local SIZE_MB=0
  if command -v du >/dev/null 2>&1; then
    SIZE_MB="$(du -sm "${SRC}" 2>/dev/null | awk '{print $1}' || true)"
    SIZE_MB="${SIZE_MB:-0}"
    [[ "${SIZE_MB}" =~ ^[0-9]+$ ]] || SIZE_MB=0
    echo "R library size: ${SIZE_MB} MB (staging guard ${MAX_MB} MB)"
    if (( SIZE_MB > MAX_MB )); then
      echo "WARNING: R library size ${SIZE_MB} MB exceeds WORKER_R_LIB_MAX_MB=${MAX_MB}; skipping staging (BeeGFS fallback)."
      return 0
    fi
  else
    echo "WARNING: du not available on this node; skipping the size guard."
  fi

  local STAGE_ROOT="${WORKER_STAGE_ROOT:-}"
  if [[ -z "${STAGE_ROOT}" ]]; then
    if [[ -d "/scratch" && -w "/scratch" ]]; then
      STAGE_ROOT="/scratch"
    else
      STAGE_ROOT="/tmp"
    fi
  fi
  local JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-unknown}}"
  local TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
  local STAGE_DIR="${STAGE_ROOT}/${JOB_ID}_${TASK_ID}"
  local TARGET="${STAGE_DIR}/${PREFIX}_R_library"
  local START=0 END=0
  if ! mkdir -p "${STAGE_DIR}"; then
    echo "ERROR: cannot create staging dir ${STAGE_DIR}; falling back to BeeGFS (task may fail on stale views)."
    return 1
  fi
  if [[ -d "${TARGET}" ]]; then
    echo "Staged R library already exists; reusing ${TARGET}"
  else
    START="$(date +%s)"
    echo "Staging R library ${SRC} -> ${TARGET} ..."
    if ! cp -a "${SRC}" "${TARGET}"; then
      echo "ERROR: cp -a staging failed (${SRC} -> ${TARGET}); falling back to BeeGFS (transient grep will requeue if stale-view)."
      return 1
    fi
    END="$(date +%s)"
    echo "Staging complete: ${TARGET} (${SIZE_MB} MB, $(( END - START )) s)"
  fi
  export R_LIBS="${TARGET}${R_LIBS:+:${R_LIBS}}"
  echo "R_LIBS=${R_LIBS}"
  return 0
}

export_worker_thread_env() {
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
         NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
}
