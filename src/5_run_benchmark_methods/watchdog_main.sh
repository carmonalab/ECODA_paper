#!/bin/bash
#SBATCH --job-name=benchmark_watchdog
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-type=END,FAIL

# Compute-node watchdog for one benchmark method array: owns the terminal wait
# + per-task gate + OOM auto-escalation (re-submitting only OOM'd datasets
# with doubled --mem, clamped to BENCHMARK_MEM_MAX) that the login submitter
# tail used to run inline. The login tail (benchmark_wait_watchdog) only
# waits for THIS job and reads its status file — so an SSH drop of the tail
# (SIGHUP) can no longer interrupt an escalation chain (observed 2026-08-14:
# 4313942_3 scitd OOM'd at 128G, no retry ever submitted, no sync email).
#
# Invoked by benchmark_submit_watchdog (benchmark_submit_common.sh) as
#   watchdog_main.sh <array_id> <label> <manifest> <mode> -- \
#     <partition> <throttle> <log_prefix> <worker_script> <runtime_export> \
#     [flags...]
# where mode is `strict` (method arrays) or `soft-gate` (prepare_pseudobulk:
# artifact-completeness pass first, strict OOM-aware gate fallback — the
# former benchmark_wait_prep_array logic). METHOD/BENCHMARK_MANIFEST/
# FORCE_BENCHMARK and the explicit runtime export propagate via sbatch
# inheritance (submitter -> watchdog -> retry workers).
#
# Protocol: writes ${WATCHDOG_STATUS_DIR}/<SLURM_JOB_ID>.status (self-named;
# the id is unknowable at submit time) before exiting — STATE=OK|FAIL, LABEL=,
# one JOB_REPORT= line per gated array, FAIL_REASON=/REPORT= on failure.
# No emailing, no NAS access, no exec-log merging (login tail only).
# Pure bash + slurm CLI (squeue/sacct/sbatch/scontrol) — no pixi/R/Python.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch the script is copied to /var/spool/slurmd/job<id>/slurm_script,
# so BASH_SOURCE no longer points at the repo. Recover the real path from the
# job record (Command= field); BASH_SOURCE fallback for login-node execution.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
source "${SCRIPT_DIR}/benchmark_submit_common.sh"

# Retry-array spec forwarded by benchmark_submit_watchdog (after the --):
# partition/throttle/log-prefix/worker-script/runtime-export + per-method
# sbatch flags (--gpus/--constraint/--cpus-per-task pins preserved).
WD_PARTITION=""
WD_THROTTLE=""
WD_LOG_PREFIX=""
WD_WORKER_SCRIPT=""
WD_RUNTIME_EXPORT=""
WD_FLAGS=()

# Retry-array submission closure (mirrors the submitters'
# submit_*_method_array_retry): writes the reduced manifest (one dataset per
# line), exports METHOD/BENCHMARK_MANIFEST, sbatch's the reduced array with
# the watchdog's forwarded spec + --mem + log patterns, echoes ONLY the new
# array job id on stdout.
watchdog_resubmit() {
  local LABEL="$1"
  local DS_CSV="$2"
  local MEM="$3"
  local MANIFEST="$4"
  local DS_LIST
  IFS=',' read -r -a DS_LIST <<< "${DS_CSV}"
  local MANIFEST_TMP="${MANIFEST}.build.$$"
  : > "${MANIFEST_TMP}"
  local name
  for name in "${DS_LIST[@]}"; do
    echo "${name}" >> "${MANIFEST_TMP}"
  done
  [[ -s "${MANIFEST_TMP}" ]] || return 1
  mv -f "${MANIFEST_TMP}" "${MANIFEST}"
  if [[ -n "${ANALYSIS_PASS:-}" ]]; then
    export ANALYSIS_MANIFEST="${MANIFEST}"
    unset BENCHMARK_MANIFEST
  else
    export BENCHMARK_MANIFEST="${MANIFEST}"
  fi
  # METHOD must be exported too: sbatch propagates only the exported
  # environment, and the workers (1.1_run_worker.sh) hard-require it.
  export METHOD="${LABEL}"

  local runtime_export="${WD_RUNTIME_EXPORT:-}"
  local validated_runtime_export
  [[ -n "${runtime_export}" ]] || {
    echo "ERROR: watchdog retry is missing the runtime export." >&2
    return 1
  }
  export ECODA_RUNTIME_PROFILE="${ECODA_RUNTIME_PROFILE:-stage5}"
  ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE:-host}" || {
    echo "ERROR: watchdog retry runtime validation failed; refusing retry submission." >&2
    return 1
  }
  validated_runtime_export="$(ecoda_runtime_export_csv \
    "${ECODA_RUNTIME_PROFILE}" "${ECODA_APPTAINER_NV:-0}")" || return 1
  [[ "${validated_runtime_export}" == "${runtime_export}" ]] || {
    echo "ERROR: watchdog runtime export differs from inherited worker export." >&2
    return 1
  }

  echo "Watchdog resubmitting ${LABEL} array (${#DS_LIST[@]} datasets, partition=${WD_PARTITION}, " >&2
  echo "  flags: ${WD_FLAGS[*]}, mem=${MEM}, throttle=${WD_THROTTLE})" >&2

  local SUBMIT_MSG
  SUBMIT_MSG=$(sbatch \
      --array="1-${#DS_LIST[@]}%${WD_THROTTLE}" \
      --partition="${WD_PARTITION}" \
      "${WD_FLAGS[@]}" \
      --time="${BENCHMARK_WORKER_TIME_LIMIT:-12:00:00}" \
      --mem="${MEM}" \
      --output="${LOGS_DIR}/${WD_LOG_PREFIX}_%A_%a.log" \
      --error="${LOGS_DIR}/${WD_LOG_PREFIX}_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      --export="ALL,${runtime_export}" \
      "${WD_WORKER_SCRIPT}")

  local ARRAY_JOB_ID
  ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "  ${LABEL} retry array job ID: ${ARRAY_JOB_ID}" >&2
  echo "${ARRAY_JOB_ID}"
}

watchdog_main() {
  local ARRAY_ID="$1"
  local LABEL="$2"
  local MANIFEST="$3"
  local MODE="$4"
  shift 4
  if [[ "${1:-}" == "--" ]]; then
    shift
  fi
  WD_PARTITION="$1"
  WD_THROTTLE="$2"
  WD_LOG_PREFIX="$3"
  WD_WORKER_SCRIPT="$4"
  WD_RUNTIME_EXPORT="$5"
  [[ -n "${WD_RUNTIME_EXPORT}" ]] || {
    echo "ERROR: watchdog runtime export is required." >&2
    return 1
  }
  shift 5
  WD_FLAGS=("$@")

  local STATUS_FILE="${WATCHDOG_STATUS_DIR}/${SLURM_JOB_ID}.status"
  mkdir -p "${WATCHDOG_STATUS_DIR}"
  echo "Watchdog ${SLURM_JOB_ID}: ${LABEL} array ${ARRAY_ID} (mode=${MODE}, manifest=${MANIFEST})"
  WATCHDOG_SCHEDULER_IDS=("${ARRAY_ID}")
  [[ -n "${SLURM_JOB_ID:-}" ]] && WATCHDOG_SCHEDULER_IDS+=("${SLURM_JOB_ID}")

  case "${MODE}" in
    strict)
      benchmark_wait_oom_retry "${ARRAY_ID}" "${LABEL}" watchdog_resubmit "${MANIFEST}" "${STATUS_FILE}"
      exit 0
      ;;
    soft-gate)
      # prepare_pseudobulk: artifact-based gate (all PB_VARIANT_NAMES
      # variants present in benchmark/pseudobulks/ per dataset — the
      # mofa/pseudobulk workers only consume those files, with an on-the-fly
      # fallback for missing ones) with strict OOM-aware task-state fallback.
      # Under --force the prep tasks MUST recompute, so the strict gate
      # always applies. Dataset list comes from the manifest (one per line).
      benchmark_wait_array_terminal "${ARRAY_ID}" "${LABEL}"
      if [[ "${FORCE_BENCHMARK:-0}" == "1" ]]; then
        echo "prepare_pseudobulk watchdog: --force requested; applying the strict OOM-aware task-state gate."
        benchmark_wait_oom_retry "${ARRAY_ID}" "${LABEL}" watchdog_resubmit "${MANIFEST}" "${STATUS_FILE}"
        exit 0
      fi
      local PB_DIR="${HPC_SCRATCH_DIR}/benchmark/pseudobulks"
      local PB_VARIANTS
      PB_VARIANTS="$(benchmark_pb_variant_names "${SCRIPT_DIR}/benchmark_hpc_utils.R")"
      if [[ -z "${PB_VARIANTS}" ]]; then
        echo "WARNING: could not parse PB_VARIANT_NAMES from benchmark_hpc_utils.R; applying the OOM-aware task-state gate." >&2
        benchmark_wait_oom_retry "${ARRAY_ID}" "${LABEL}" watchdog_resubmit "${MANIFEST}" "${STATUS_FILE}"
        exit 0
      fi
      local MISSING=() PREP_DS=() DS V
      while IFS= read -r DS; do
        PREP_DS+=("${DS}")
      done < "${MANIFEST}"
      for DS in "${PREP_DS[@]}"; do
        for V in ${PB_VARIANTS}; do
          if [[ ! -f "${PB_DIR}/${DS}_pseudobulk_${V}.rds" ]]; then
            MISSING+=("${DS}/${V}")
          fi
        done
      done
      if [[ ${#MISSING[@]} -eq 0 ]]; then
        echo "prepare_pseudobulk: all $(wc -w <<< "${PB_VARIANTS}" | tr -d ' ') variants x ${#PREP_DS[@]} datasets present in ${PB_DIR}; skipping the task-state gate."
        # Visibility only: show non-COMPLETED tasks without failing
        # (best-effort; sacct may be empty or the master line may itself be
        # COMPLETED).
        sacct -j "${ARRAY_ID}" --format=JobID,State,ExitCode -n 2>/dev/null \
          | grep -v 'COMPLETED' | sed 's/^/  note: /' || true
        WATCHDOG_GATED_REPORTS=("${LABEL}|${ARRAY_ID}|$(array_wall_time "${ARRAY_ID}")")
        benchmark_write_status_file "${STATUS_FILE}" OK "${LABEL}"
        exit 0
      fi
      echo "prepare_pseudobulk: missing variant file(s): ${MISSING[*]}; applying the OOM-aware task-state gate."
      benchmark_wait_oom_retry "${ARRAY_ID}" "${LABEL}" watchdog_resubmit "${MANIFEST}" "${STATUS_FILE}"
      exit 0
      ;;
    *)
      echo "ERROR: unknown watchdog mode '${MODE}' (expected strict or soft-gate)." >&2
      exit 1
      ;;
  esac
}

watchdog_main "$@"
