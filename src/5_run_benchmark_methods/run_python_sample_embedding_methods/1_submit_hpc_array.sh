#!/bin/bash
set -euo pipefail

# Submit Python benchmark method arrays (MrVI/scPoli/PILOT-GM-VAE on GPU,
# PILOT/QOT on CPU)
#
# Usage:
#   ./1_submit_hpc_array.sh                                    # all methods, all benchmark datasets
#   ./1_submit_hpc_array.sh --ds_name _debug --methods mrvi    # single dataset, single method
#   ./1_submit_hpc_array.sh --methods mrvi,scpoli --force      # recompute existing feathers
#   ./1_submit_hpc_array.sh --partition debug-cpu              # override per-method partitions (drops the constraint pin)
#   ./1_submit_hpc_array.sh --sync-only 12345,12346            # resume: skip submission, re-check + sync
#                                                              # (repeat the original --ds_name/--methods flags;
#                                                              #  watchdog ids resume via their status files)
#
# One SLURM array per method; array task IDs map 1:1 to lines of the
# per-method manifest ${HPC_SCRATCH_DIR}/benchmark_manifest_<method>_<pid>.txt
# (written by this script, one dataset per line; PID suffix so overlapping
# submissions do not clobber queued arrays' manifests). Hardware is PINNED for
# runtime comparability (see slurm_config.sh benchmark vars): GPU methods
# (mrvi, scpoli, pilotgm) on ${SLURM_PARTITION_BENCHMARK_GPU} with
# ${BENCHMARK_GPU_CONSTRAINT}, PILOT/QOT on ${SLURM_PARTITION_BENCHMARK_CPU} with
# ${BENCHMARK_CPU_CONSTRAINT}. Flags are passed on the sbatch command line
# (SLURM directives do not expand env vars). An explicit --partition <P>
# override (e.g. --partition ${SLURM_PARTITION_PRIVATE} for _debug runs on the
# private node, or an ad-hoc --partition shared-gpu for debug) DROPS the
# method's --constraint flag — an explicit partition
# choice means the user accepts non-pinned hardware; keeping the constraint
# would otherwise hang jobs PENDING forever on nodes whose CPU/GPU differ.
# After all arrays complete:
# NAS reachability check -> merge per-task exec logs (job-id scoped,
# existing-log continuity) -> sync to NAS -> only then delete this run's
# per-task logs.
#
# Each method array gets its own compute-node WATCHDOG (watchdog_main.sh via
# benchmark_submit_watchdog) that owns the terminal wait + per-task gate +
# OOM auto-escalation, so an SSH drop of this login tail can never interrupt
# an escalation chain; this tail only waits for the watchdog jobs
# (benchmark_wait_watchdog) and then syncs. OOM escalation: an
# OUT_OF_MEMORY task cannot self-requeue (the task is dead), so the watchdog
# re-submits only the OOM'd tasks' datasets with doubled --mem (128G -> 256G
# -> 500G, clamped to the BENCHMARK_MEM_MAX ceiling) via its own resubmit
# closure (same per-method partition/--gpus/--constraint/throttle), before
# failing closed with an OOM report in its status file; the tail emails on
# FAIL. Non-OOM failures fail closed exactly as before.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/../benchmark_submit_common.sh"
cd "${PROJECT_ROOT}"

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi

DS_NAME_ARG=""
METHODS_ARG=""
PARTITION_ARG=""
FORCE_ARG=0
SYNC_ONLY_IDS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ds_name)
      DS_NAME_ARG="${2:-}"
      shift 2
      ;;
    --ds_name=*)
      DS_NAME_ARG="${1#*=}"
      shift
      ;;
    --methods)
      METHODS_ARG="${2:-}"
      shift 2
      ;;
    --methods=*)
      METHODS_ARG="${1#*=}"
      shift
      ;;
    --partition)
      PARTITION_ARG="${2:-}"
      shift 2
      ;;
    --partition=*)
      PARTITION_ARG="${1#*=}"
      shift
      ;;
    --force)
      FORCE_ARG=1
      shift
      ;;
    --sync-only)
      SYNC_ONLY_IDS="${2:-}"
      if [[ -z "${SYNC_ONLY_IDS}" ]]; then
        echo "ERROR: --sync-only requires at least one job id (comma-separated)."
        exit 1
      fi
      shift 2
      ;;
    --sync-only=*)
      SYNC_ONLY_IDS="${1#*=}"
      if [[ -z "${SYNC_ONLY_IDS}" ]]; then
        echo "ERROR: --sync-only requires at least one job id (comma-separated)."
        exit 1
      fi
      shift
      ;;
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

if [[ -n "${SYNC_ONLY_IDS}" && ${FORCE_ARG} -eq 1 ]]; then
  echo "ERROR: --sync-only cannot be combined with --force."
  exit 1
fi

# Passed to workers via the environment (sbatch propagates the submit
# script's environment); 1.1_run_worker.sh forwards it to the py script.
export FORCE_BENCHMARK="${FORCE_ARG}"

# ---------------------------------------------------------------------------
# Resolve datasets (shared helper; skips `_*` keys unless --ds_name).
# ---------------------------------------------------------------------------
benchmark_resolve_datasets "${DS_NAME_ARG}"

METHODS=(mrvi scpoli pilot qot pilotgm)
if [[ -n "${METHODS_ARG}" ]]; then
  IFS=',' read -r -a METHODS <<< "${METHODS_ARG}"
fi

# ---------------------------------------------------------------------------
# Submit one array per method
# ---------------------------------------------------------------------------
if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Sync-only resume mode: jobs ${SYNC_ONLY_IDS} (no submission) ==="
  IFS=',' read -r -a ARRAY_JOB_IDS <<< "${SYNC_ONLY_IDS}"
  # Watchdog ids resume via their status files (wait + read STATE); array ids
  # keep the strict all-rows gate.
  for i in "${!ARRAY_JOB_IDS[@]}"; do
    JOB_ID="${ARRAY_JOB_IDS[$i]}"
    STATUS_FILE="${WATCHDOG_STATUS_DIR}/${JOB_ID}.status"
    if [[ -f "${STATUS_FILE}" ]]; then
      SYNC_LABEL="$(grep -E '^LABEL=' "${STATUS_FILE}" | head -1 | cut -d= -f2- | tr -d '[:space:]')"
      echo "  Job ${JOB_ID}: watchdog status file found; gating via the watchdog."
      benchmark_wait_watchdog "${JOB_ID}" "${SYNC_LABEL:-watchdog}"
    else
      benchmark_wait_for_array "${JOB_ID}" ""
    fi
  done
else
echo "=== Submitting Python benchmark method arrays ==="

mkdir -p "${LOGS_DIR}"

# Prints ONLY the array job id on stdout (progress goes to stderr) so the
# caller can capture the id with $(...) — a multi-line capture would break
# the sacct gates downstream.
#
# Retry-capable: the OOM auto-escalation loop inside the per-method watchdog
# (watchdog_main.sh -> benchmark_wait_oom_retry) re-invokes ITS OWN closure
# with a comma-separated REDUCED dataset list, a doubled --mem and a fresh
# manifest path (only the OOM'd datasets; array task ids map 1:1 to its lines
# via the existing sed -n mechanism — no worker changes needed). This
# submit-side function is used only for the normal full-dataset submission;
# per-method partition/--gpus/--constraint/throttle come from py_method_spec
# below (same as the watchdog's forwarded spec, PARTITION_ARG override
# included).

# Per-method submission spec (partition/throttle/flags) into the globals
# PY_METHOD_PARTITION/PY_METHOD_THROTTLE/PY_METHOD_FLAGS. An explicit
# --partition override drops the --constraint pin (kept:
# --gpus/--cpus-per-task/--mem).
py_method_spec() {
  local METHOD="$1"
  PY_METHOD_PARTITION=""
  PY_METHOD_THROTTLE=""
  PY_METHOD_FLAGS=()
  case "${METHOD}" in
    mrvi|scpoli|pilotgm)
      PY_METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_GPU}"
      PY_METHOD_FLAGS=(--gpus="${BENCHMARK_GPU_COUNT}"
                       --constraint="${BENCHMARK_GPU_CONSTRAINT}"
                       --cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}")
      PY_METHOD_THROTTLE="${BENCHMARK_GPU_ARRAY_THROTTLE}"
      ;;
    pilot|qot)
      PY_METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
      PY_METHOD_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}"
                       --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
      PY_METHOD_THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
      ;;
    *)
      echo "ERROR: Unknown method '${METHOD}' (expected mrvi, scpoli, pilot, qot or pilotgm)." >&2
      exit 1
      ;;
  esac
  if [[ -n "${PARTITION_ARG}" ]]; then
    PY_METHOD_PARTITION="${PARTITION_ARG}"
    # Explicit partition choice = user accepts non-pinned hardware: drop the
    # --constraint pin (kept: --gpus/--cpus-per-task/--mem). Without this,
    # e.g. --partition private-carmona-gpu would sit PENDING forever — its
    # CPU/GPU never match the pinned BENCHMARK_*_CONSTRAINT.
    local FILTERED_FLAGS=() FLAG
    for FLAG in "${PY_METHOD_FLAGS[@]}"; do
      case "${FLAG}" in
        --constraint|--constraint=*)
          ;;
        *)
          FILTERED_FLAGS+=("${FLAG}")
          ;;
      esac
    done
    PY_METHOD_FLAGS=("${FILTERED_FLAGS[@]}")
    echo "  NOTE: --partition override drops the --constraint hardware pin" >&2
  fi
}

submit_python_method_array_retry() {
  local METHOD="$1"
  local DS_CSV="$2"
  local MEM="$3"
  local MANIFEST="$4"
  py_method_spec "${METHOD}"

  # Per-method manifest: one dataset per line; rebuilt every run. The name
  # carries the submit PID so an overlapping second submission cannot
  # clobber the manifest of already-submitted (queued) arrays.
  : > "${MANIFEST}"
  IFS=',' read -r -a DS_LIST <<< "${DS_CSV}"
  for name in "${DS_LIST[@]}"; do
    echo "${name}" >> "${MANIFEST}"
  done
  export BENCHMARK_MANIFEST="${MANIFEST}"
  # METHOD must be exported too: sbatch propagates only the exported
  # environment, and 1.1_run_worker.sh hard-requires it.
  export METHOD="${METHOD}"

  echo "Submitting ${METHOD} array (${#DS_LIST[@]} datasets, partition=${PY_METHOD_PARTITION}, " >&2
  echo "  flags: ${PY_METHOD_FLAGS[*]}, mem=${MEM}, throttle=${PY_METHOD_THROTTLE})" >&2

  local SUBMIT_MSG
  SUBMIT_MSG=$(sbatch \
      --array="1-${#DS_LIST[@]}%${PY_METHOD_THROTTLE}" \
      --partition="${PY_METHOD_PARTITION}" \
      "${PY_METHOD_FLAGS[@]}" \
      --mem="${MEM}" \
      --output="${LOGS_DIR}/5_benchmark_${METHOD}_%A_%a.log" \
      --error="${LOGS_DIR}/5_benchmark_${METHOD}_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      "${SCRIPT_DIR}/1.1_run_worker.sh")

  local ARRAY_JOB_ID
  ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "  ${METHOD} array job ID: ${ARRAY_JOB_ID}" >&2
  echo "${ARRAY_JOB_ID}"
}

ARRAY_JOB_IDS=()
ARRAY_JOB_METHODS=()
WATCHDOG_JOB_IDS=()
for METHOD in "${METHODS[@]}"; do
  DS_CSV=""
  for name in "${DATASET_NAMES[@]}"; do
    [[ -n "${DS_CSV}" ]] && DS_CSV+=","
    DS_CSV+="${name}"
  done
  NEW_ARRAY_ID="$(submit_python_method_array_retry "${METHOD}" "${DS_CSV}" "${BENCHMARK_MEM}" \
    "${HPC_SCRATCH_DIR}/benchmark_manifest_${METHOD}_$$.txt")"
  ARRAY_JOB_IDS+=("${NEW_ARRAY_ID}")
  ARRAY_JOB_METHODS+=("${METHOD}")
  # Per-method watchdog: owns the terminal wait + gate + OOM escalation on a
  # compute node (survives SSH drops of this tail); per-method partition /
  # throttle / GPU+CPU flags preserved for its retry-array submissions.
  py_method_spec "${METHOD}"
  WATCHDOG_JOB_IDS+=("$(benchmark_submit_watchdog \
    "${NEW_ARRAY_ID}" "${METHOD}" \
    "${HPC_SCRATCH_DIR}/benchmark_manifest_${METHOD}_$$.txt" \
    strict "${PY_METHOD_PARTITION}" "${PY_METHOD_THROTTLE}" \
    "5_benchmark_${METHOD}" "${SCRIPT_DIR}/1.1_run_worker.sh" \
    "${PY_METHOD_FLAGS[@]}")")
done
fi

# ---------------------------------------------------------------------------
# Monitor & verify & sync results back to NAS (shared tail).
# Each array is gated via its watchdog: an OUT_OF_MEMORY task's dataset is
# re-submitted by the watchdog with doubled --mem (ceiling BENCHMARK_MEM_MAX);
# non-OOM failures fail closed as before. The watchdog runs on a compute
# node, so an SSH drop of this tail cannot interrupt the escalation.
# ---------------------------------------------------------------------------
if [[ -z "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Monitoring job completion ==="
  for i in "${!ARRAY_JOB_IDS[@]}"; do
    benchmark_wait_watchdog "${WATCHDOG_JOB_IDS[$i]}" "${ARRAY_JOB_METHODS[$i]}"
  done
fi

# Labels for the exec-log merge = the submitted methods.
benchmark_merge_sync_cleanup "${METHODS[@]}"
