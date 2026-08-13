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
#                                                                                                              # (repeat the original --ds_name/--methods flags)
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
# Submit order: prepare_pseudobulk array FIRST, waited until it leaves the
# scheduler, then gated on ARTIFACT COMPLETENESS (all PB_VARIANT_NAMES
# variants present in benchmark/pseudobulks/ per dataset — the mofa/pseudobulk
# workers only consume those files, with an on-the-fly fallback for missing
# ones) instead of strict task states: a prep task that failed on a transient
# node issue (stale BeeGFS view -> "missing from the pixi environment") must
# not block the method arrays when its variants already exist on disk. The
# strict fail-closed sacct gate applies only under --force (prep must
# recompute) or when variant files are actually missing — both strict paths
# are OOM-AWARE (see below). THEN the mofa/pseudobulk arrays (and
# gloscope/scitd) are submitted, waited + gated the same way. If mofa or
# pseudobulk is requested without prepare_pseudobulk it is auto-prepended.
# After all arrays complete the shared merge/sync/cleanup tail runs (NAS
# reachability check -> RDS integrity sidecar -> fail-closed sacct gate ->
# merge per-task exec logs -> sync to NAS -> only then delete per-task
# logs) — see benchmark_submit_common.sh.
#
# OOM auto-escalation: an OUT_OF_MEMORY task (e.g. scitd on a large dataset
# via process_scitd_fig's full-matrix scITD container) can not self-requeue
# (the task is dead), so the gates here run through
# benchmark_wait_oom_retry: only the OOM'd tasks' datasets are re-submitted
# with doubled --mem (128G -> 256G -> 512G) via submit_method_array_retry,
# up to the BENCHMARK_MEM_MAX ceiling, before failing closed with an OOM
# report. Non-OOM failures fail closed exactly as before.

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
# prepare_pseudobulk gate: artifact-based (soft) with strict fallback.
#
# The mofa/pseudobulk workers consume ONLY the shared variant files in
# ${HPC_SCRATCH_DIR}/benchmark/pseudobulks/ (with an on-the-fly fallback for
# missing variants), so a prep task that failed on a transient node issue
# (stale BeeGFS view -> "missing from the pixi environment") while all its
# variants already exist on disk must NOT block the method arrays. Wait for
# the array to leave the scheduler (same squeue poll as
# benchmark_wait_for_array), then check artifact completeness; if every
# variant exists per dataset, proceed without the strict task-state gate.
# Under --force the prep tasks MUST recompute, so the strict gate always
# applies; when variant files are missing (or PB_VARIANT_NAMES cannot be
# parsed) the OOM-aware gate applies instead of the plain strict one:
# OUT_OF_MEMORY prep tasks are re-submitted with doubled --mem (only those
# datasets), all-COMPLETED passes, non-OOM failures fail closed. The soft
# gate is NOT used in --sync-only mode (all provided job ids are gated
# strictly there).
# PB_VARIANT_NAMES (benchmark_hpc_utils.R) is the single source of truth;
# parsed here so no second list has to be maintained.
# ---------------------------------------------------------------------------
benchmark_wait_prep_array() {
  local JOB_ID="$1"
  local PREP_MANIFEST="${HPC_SCRATCH_DIR}/benchmark_manifest_prepare_pseudobulk_$$.txt"
  echo "=== Monitoring prepare_pseudobulk array ${JOB_ID} ==="
  # Block until the job leaves the scheduler (same poll as the shared gate).
  while squeue -u "$USER" -h -o "%A" 2>/dev/null | grep -qx "${JOB_ID}"; do
    sleep 60
  done
  echo "prepare_pseudobulk array ${JOB_ID} left the scheduler."

  if [[ ${FORCE_ARG} -eq 1 ]]; then
    echo "prepare_pseudobulk: --force requested; applying the strict task-state gate."
    benchmark_wait_for_array "${JOB_ID}" "prepare_pseudobulk"
    return 0
  fi

  local PB_DIR="${HPC_SCRATCH_DIR}/benchmark/pseudobulks"
  local PB_VARIANTS
  PB_VARIANTS="$(sed -n '/^PB_VARIANT_NAMES <- c(/,/^)/p' \
    "${SCRIPT_DIR}/../benchmark_hpc_utils.R" | grep -oE '"[a-zA-Z0-9_.]+"' | tr -d '"')"
  if [[ -z "${PB_VARIANTS}" ]]; then
    echo "WARNING: could not parse PB_VARIANT_NAMES from benchmark_hpc_utils.R; applying the OOM-aware task-state gate."
    benchmark_wait_oom_retry "${JOB_ID}" "prepare_pseudobulk" submit_method_array_retry "${PREP_MANIFEST}"
    return 0
  fi

  local MISSING=()
  local DS V
  for DS in "${DATASET_NAMES[@]}"; do
    for V in ${PB_VARIANTS}; do
      if [[ ! -f "${PB_DIR}/${DS}_pseudobulk_${V}.rds" ]]; then
        MISSING+=("${DS}/${V}")
      fi
    done
  done
  if [[ ${#MISSING[@]} -eq 0 ]]; then
    echo "prepare_pseudobulk: all $(wc -w <<< "${PB_VARIANTS}" | tr -d ' ') variants x ${NUM_DATASETS} datasets present in ${PB_DIR}; skipping the task-state gate."
    # Visibility only: show non-COMPLETED tasks without failing (best-effort;
    # sacct may be empty or the master line may itself be COMPLETED).
    sacct -j "${JOB_ID}" --format=JobID,State,ExitCode -n 2>/dev/null \
      | grep -v 'COMPLETED' | sed 's/^/  note: /' || true
    return 0
  fi
  echo "prepare_pseudobulk: missing variant file(s): ${MISSING[*]}; applying the OOM-aware task-state gate."
  benchmark_wait_oom_retry "${JOB_ID}" "prepare_pseudobulk" submit_method_array_retry "${PREP_MANIFEST}"
}

# ---------------------------------------------------------------------------
# Submit one array per method (prepare_pseudobulk first, then the rest)
# ---------------------------------------------------------------------------
if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Sync-only resume mode: jobs ${SYNC_ONLY_IDS} (no submission, no prepare_pseudobulk wait) ==="
  IFS=',' read -r -a ARRAY_JOB_IDS <<< "${SYNC_ONLY_IDS}"
else
echo "=== Submitting R benchmark method arrays ==="

# Prints ONLY the array job id on stdout (progress goes to stderr) so the
# caller can capture the id with $(...) — a multi-line capture would break
# the sacct gate in benchmark_wait_oom_retry.
#
# Retry-capable: the OOM auto-escalation loop (benchmark_wait_oom_retry)
# re-invokes this with a comma-separated REDUCED dataset list, a doubled
# --mem and a fresh manifest path (only the OOM'd datasets; array task ids
# map 1:1 to its lines via the existing sed -n mechanism — no worker changes
# needed). Flags/partition/throttle are identical to the normal submission.
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

for METHOD in "${METHODS[@]}"; do
  if [[ "${METHOD}" == "prepare_pseudobulk" ]]; then
    PREP_JOB_ID="$(submit_method_array "${METHOD}")"
    break
  fi
done

# prepare_pseudobulk must complete before the mofa/pseudobulk arrays start;
# the gate is artifact-based (soft) with an OOM-aware fail-closed fallback.
if [[ -n "${PREP_JOB_ID}" ]]; then
  benchmark_wait_prep_array "${PREP_JOB_ID}"
fi

for METHOD in "${METHODS[@]}"; do
  if [[ "${METHOD}" == "prepare_pseudobulk" ]]; then
    continue
  fi
  ARRAY_JOB_IDS+=("$(submit_method_array "${METHOD}")")
  ARRAY_JOB_METHODS+=("${METHOD}")
done
fi

# ---------------------------------------------------------------------------
# Monitor & verify & sync results back to NAS (shared tail).
# Each array is gated OOM-aware: an OUT_OF_MEMORY task's dataset is
# re-submitted with doubled --mem (submit_method_array_retry, ceiling
# BENCHMARK_MEM_MAX); non-OOM failures fail closed as before.
# ---------------------------------------------------------------------------
for i in "${!ARRAY_JOB_IDS[@]}"; do
  benchmark_wait_oom_retry "${ARRAY_JOB_IDS[$i]}" "${ARRAY_JOB_METHODS[$i]}" \
    submit_method_array_retry "${HPC_SCRATCH_DIR}/benchmark_manifest_${ARRAY_JOB_METHODS[$i]}_$$.txt"
done

# Labels for the exec-log merge = the submitted methods (includes the
# auto-prepended prepare_pseudobulk, whose worker also writes a per-task log).
benchmark_merge_sync_cleanup "${METHODS[@]}"
