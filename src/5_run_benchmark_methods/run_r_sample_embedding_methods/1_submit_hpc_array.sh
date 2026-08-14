#!/bin/bash
set -euo pipefail

# Submit R benchmark method arrays (GloScope, MOFA, Pseudobulk, scITD) on CPU
#
# Usage:
#   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh   # all methods, all benchmark datasets
#   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --ds_name _debug --methods prepare_pseudobulk,pseudobulk
#   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --methods mofa --force
#   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --partition debug-cpu
#   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --sync-only 12345,12346   # resume: skip submission, re-check + sync
#                                                                                                              # (repeat the original --ds_name/--methods flags;
#                                                                                                              #  watchdog ids resume via their status files)
#
# One SLURM array per method; array task IDs map 1:1 to lines of the
# per-method manifest ${HPC_SCRATCH_DIR}/benchmark_manifest_<method>_<pid>.txt
# (written by this script, one dataset per line; PID suffix so overlapping
# submissions do not clobber queued arrays' manifests). All methods are PINNED
# to the CPU benchmark class (${SLURM_PARTITION_BENCHMARK_CPU} with
# ${BENCHMARK_CPU_CONSTRAINT} EPYC-7742, ${BENCHMARK_CPU_CPUS_PER_TASK} cpus,
# ${BENCHMARK_MEM}) — the same class as PILOT, so cross-method runtime
# comparisons stay valid. Flags are passed on the sbatch command line (SLURM
# directives do not expand env vars). An explicit --partition <P> override
# DROPS the --constraint pin (explicit partition choice = user accepts
# non-pinned hardware; keeping the constraint would hang jobs PENDING forever
# on nodes whose CPU differs).
#
# Submit order: prepare_pseudobulk array FIRST with its own compute-node
# WATCHDOG (soft-gate mode: gated on ARTIFACT COMPLETENESS — all
# PB_VARIANT_NAMES variants present in benchmark/pseudobulks/ per dataset —
# instead of strict task states, so a prep task that failed on a transient
# node issue must not block the method arrays when its variants already exist
# on disk; the strict OOM-aware sacct gate applies only under --force or when
# variant files are actually missing), waited via its watchdog id. THEN the
# mofa/pseudobulk arrays (and gloscope/scitd) are submitted, each with its own
# strict watchdog, and waited the same way. If mofa or pseudobulk is requested
# without prepare_pseudobulk it is auto-prepended. Every gate + OOM escalation
# runs in the watchdog jobs (compute nodes) so an SSH drop of this login tail
# can never interrupt an escalation; this tail only waits for the watchdog
# jobs and then runs the shared merge/sync/cleanup tail (NAS reachability
# check -> RDS integrity sidecar -> merge per-task exec logs -> sync to NAS ->
# only then delete per-task logs) — see benchmark_submit_common.sh.
#
# OOM auto-escalation: an OUT_OF_MEMORY task (e.g. scitd on a large dataset
# via process_scitd_fig's full-matrix scITD container) can not self-requeue
# (the task is dead), so each watchdog (watchdog_main.sh, submitted via
# benchmark_submit_watchdog) re-submits only the OOM'd tasks' datasets with
# doubled --mem (128G -> 256G -> 500G, clamped to the BENCHMARK_MEM_MAX
# ceiling), before failing closed with an OOM report in its status file; the
# tail's benchmark_wait_watchdog reads the file and emails on FAIL.
# Non-OOM failures fail closed exactly as before.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/../benchmark_submit_common.sh"
cd "${PROJECT_ROOT}"

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
# script's environment); 1.1_run_worker.sh forwards it to the R scripts.
export FORCE_BENCHMARK="${FORCE_ARG}"

# Resolve datasets (shared helper; skips `_*` keys unless --ds_name).
benchmark_resolve_datasets "${DS_NAME_ARG}"

METHODS=(gloscope mofa pseudobulk scitd)
if [[ -n "${METHODS_ARG}" ]]; then
  IFS=',' read -r -a METHODS <<< "${METHODS_ARG}"
fi

# Validate against the known method set; dedupe preserving order.
VALID_METHODS="prepare_pseudobulk|gloscope|mofa|pseudobulk|scitd"
UNIQUE_METHODS=()
for M in "${METHODS[@]}"; do
  if ! [[ "${M}" =~ ^(${VALID_METHODS})$ ]]; then
    echo "ERROR: Unknown method '${M}' (expected prepare_pseudobulk, gloscope, mofa, pseudobulk or scitd)."
    exit 1
  fi
  if [[ ! " ${UNIQUE_METHODS[*]} " =~ " ${M} " ]]; then
    UNIQUE_METHODS+=("${M}")
  fi
done
METHODS=("${UNIQUE_METHODS[@]}")

# mofa/pseudobulk consume the shared DESeq2 pseudobulks: auto-prepend the
# prepare_pseudobulk prep step if not explicitly listed.
NEEDS_PREP=0
for M in "${METHODS[@]}"; do
  if [[ "${M}" == "mofa" || "${M}" == "pseudobulk" ]]; then
    NEEDS_PREP=1
  fi
done
if [[ ${NEEDS_PREP} -eq 1 ]]; then
  HAS_PREP=0
  for M in "${METHODS[@]}"; do
    if [[ "${M}" == "prepare_pseudobulk" ]]; then
      HAS_PREP=1
    fi
  done
  if [[ ${HAS_PREP} -eq 0 ]]; then
    echo "NOTE: mofa/pseudobulk requested without prepare_pseudobulk; auto-prepending it."
    METHODS=("prepare_pseudobulk" "${METHODS[@]}")
  fi
fi

# ---------------------------------------------------------------------------
# Hardware: all R benchmark methods pinned to the CPU benchmark class
# (EPYC-7742, BENCHMARK_CPU_CPUS_PER_TASK cpus, BENCHMARK_MEM) for runtime
# comparability with PILOT; --partition override drops the constraint pin.
# ---------------------------------------------------------------------------
PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
EXTRA_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}"
             --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
if [[ -n "${PARTITION_ARG}" ]]; then
  PARTITION="${PARTITION_ARG}"
  # Explicit partition choice = user accepts non-pinned hardware: drop the
  # --constraint pin (kept: --cpus-per-task/--mem). Without this, e.g.
  # --partition private-carmona-gpu would sit PENDING forever — its CPU never
  # matches the pinned BENCHMARK_CPU_CONSTRAINT.
  FILTERED_FLAGS=()
  for FLAG in "${EXTRA_FLAGS[@]}"; do
    case "${FLAG}" in
      --constraint|--constraint=*)
        ;;
      *)
        FILTERED_FLAGS+=("${FLAG}")
        ;;
    esac
  done
  EXTRA_FLAGS=("${FILTERED_FLAGS[@]}")
  echo "  NOTE: --partition override drops the --constraint hardware pin"
fi

mkdir -p "${LOGS_DIR}"

# ---------------------------------------------------------------------------
# Submit one array per method (prepare_pseudobulk first, then the rest), each
# with its own compute-node watchdog (benchmark_submit_watchdog): the
# watchdog owns the terminal wait + gate + OOM escalation (survives SSH drops
# of this login tail), and this tail gates on the watchdog id via
# benchmark_wait_watchdog.
# ---------------------------------------------------------------------------
if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Sync-only resume mode: jobs ${SYNC_ONLY_IDS} (no submission, no prepare_pseudobulk wait) ==="
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
echo "=== Submitting R benchmark method arrays ==="

# Prints ONLY the array job id on stdout (progress goes to stderr) so the
# caller can capture the id with $(...) — a multi-line capture would break
# the sacct gates downstream.
#
# Retry-capable: the OOM auto-escalation loop inside the per-method watchdog
# (watchdog_main.sh -> benchmark_wait_oom_retry) re-invokes ITS OWN closure
# with a comma-separated REDUCED dataset list, a doubled --mem and a fresh
# manifest path (only the OOM'd datasets; array task ids map 1:1 to its lines
# via the existing sed -n mechanism — no worker changes needed). This
# submit-side function is used only for the normal full-dataset submission
# (submit_method_array below); flags/partition/throttle mirror the watchdog's
# forwarded spec.
submit_method_array_retry() {
  local METHOD="$1"
  local DS_CSV="$2"
  local MEM="$3"
  local MANIFEST="$4"
  IFS=',' read -r -a DS_LIST <<< "${DS_CSV}"
  : > "${MANIFEST}"
  for name in "${DS_LIST[@]}"; do
    echo "${name}" >> "${MANIFEST}"
  done
  export BENCHMARK_MANIFEST="${MANIFEST}"
  # METHOD must be exported too: sbatch propagates only the exported
  # environment, and 1.1_run_worker.sh hard-requires it.
  export METHOD="${METHOD}"

  echo "Submitting ${METHOD} array (${#DS_LIST[@]} datasets, partition=${PARTITION}, " >&2
  echo "  flags: ${EXTRA_FLAGS[*]}, mem=${MEM}, throttle=${MAX_NUM_CHUNKS_PARALLEL})" >&2

  local SUBMIT_MSG
  SUBMIT_MSG=$(sbatch \
      --array="1-${#DS_LIST[@]}%${MAX_NUM_CHUNKS_PARALLEL}" \
      --partition="${PARTITION}" \
      "${EXTRA_FLAGS[@]}" \
      --mem="${MEM}" \
      --output="${LOGS_DIR}/5_benchmark_r_${METHOD}_%A_%a.log" \
      --error="${LOGS_DIR}/5_benchmark_r_${METHOD}_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      "${SCRIPT_DIR}/1.1_run_worker.sh")

  local ARRAY_JOB_ID
  ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "  ${METHOD} array job ID: ${ARRAY_JOB_ID}" >&2
  echo "${ARRAY_JOB_ID}"
}

# Normal full-dataset submission: delegates to submit_method_array_retry with
# the resolved DATASET_NAMES and the default BENCHMARK_MEM.
submit_method_array() {
  local METHOD="$1"
  local DS_CSV=""
  for name in "${DATASET_NAMES[@]}"; do
    [[ -n "${DS_CSV}" ]] && DS_CSV+=","
    DS_CSV+="${name}"
  done
  submit_method_array_retry "${METHOD}" "${DS_CSV}" "${BENCHMARK_MEM}" \
    "${HPC_SCRATCH_DIR}/benchmark_manifest_${METHOD}_$$.txt"
}

PREP_JOB_ID=""
ARRAY_JOB_IDS=()
ARRAY_JOB_METHODS=()
WATCHDOG_JOB_IDS=()

for METHOD in "${METHODS[@]}"; do
  if [[ "${METHOD}" == "prepare_pseudobulk" ]]; then
    PREP_JOB_ID="$(submit_method_array "${METHOD}")"
    break
  fi
done

# prepare_pseudobulk must complete before the mofa/pseudobulk arrays start.
# Its soft-gate watchdog (artifact completeness pass first, strict OOM-aware
# fallback) owns the gate; this tail blocks on the watchdog id. FORCE_ARG is
# exported as FORCE_BENCHMARK above, so the watchdog applies the strict gate
# under --force.
if [[ -n "${PREP_JOB_ID}" ]]; then
  PREP_WATCHDOG_ID="$(benchmark_submit_watchdog \
    "${PREP_JOB_ID}" "prepare_pseudobulk" \
    "${HPC_SCRATCH_DIR}/benchmark_manifest_prepare_pseudobulk_$$.txt" \
    soft-gate "${PARTITION}" "${MAX_NUM_CHUNKS_PARALLEL}" \
    "5_benchmark_r_prepare_pseudobulk" "${SCRIPT_DIR}/1.1_run_worker.sh" \
    "${EXTRA_FLAGS[@]}")"
  benchmark_wait_watchdog "${PREP_WATCHDOG_ID}" "prepare_pseudobulk"
fi

for METHOD in "${METHODS[@]}"; do
  if [[ "${METHOD}" == "prepare_pseudobulk" ]]; then
    continue
  fi
  NEW_ARRAY_ID="$(submit_method_array "${METHOD}")"
  ARRAY_JOB_IDS+=("${NEW_ARRAY_ID}")
  ARRAY_JOB_METHODS+=("${METHOD}")
  WATCHDOG_JOB_IDS+=("$(benchmark_submit_watchdog \
    "${NEW_ARRAY_ID}" "${METHOD}" \
    "${HPC_SCRATCH_DIR}/benchmark_manifest_${METHOD}_$$.txt" \
    strict "${PARTITION}" "${MAX_NUM_CHUNKS_PARALLEL}" \
    "5_benchmark_r_${METHOD}" "${SCRIPT_DIR}/1.1_run_worker.sh" \
    "${EXTRA_FLAGS[@]}")")
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
  for i in "${!ARRAY_JOB_IDS[@]}"; do
    benchmark_wait_watchdog "${WATCHDOG_JOB_IDS[$i]}" "${ARRAY_JOB_METHODS[$i]}"
  done
fi

# Labels for the exec-log merge = the submitted methods (includes the
# auto-prepended prepare_pseudobulk, whose worker also writes a per-task log).
benchmark_merge_sync_cleanup "${METHODS[@]}"
