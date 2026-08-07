#!/bin/bash
set -euo pipefail

# Submit the ECODA transformation analysis and zero-imputation analysis
# arrays (Pipeline B) on CPU
#
# Usage:
#   ./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh   # both analyses, all benchmark datasets
#   ./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --ds_name _debug --analysis trans,zeroimp
#   ./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --analysis trans --force
#   ./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --partition debug-cpu
#
# Two SLURM arrays (one per analysis: trans, zeroimp); array task IDs map 1:1
# to lines of the per-analysis manifest
# ${HPC_SCRATCH_DIR}/benchmark_manifest_<analysis>_<pid>.txt (one dataset per
# line; PID suffix so overlapping submissions do not clobber queued arrays'
# manifests). Partition ${SLURM_PARTITION} defaults to shared-cpu with 4 cpus
# and 32G (--cpus-per-task/--mem/--time passed on the sbatch command line —
# SLURM directives do not expand env vars); these two analyses are NOT part
# of the pinned runtime comparisons, so no --constraint pin is applied.
# After all arrays complete the shared merge/sync/cleanup tail runs (NAS
# reachability check -> RDS integrity sidecar -> fail-closed sacct gate ->
# merge per-task exec logs -> sync to NAS -> only then delete this run's
# per-task logs) — see benchmark_submit_common.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/../benchmark_submit_common.sh"
cd "${PROJECT_ROOT}"

DS_NAME_ARG=""
ANALYSIS_ARG=""
PARTITION_ARG=""
FORCE_ARG=0
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
    --analysis)
      ANALYSIS_ARG="${2:-}"
      shift 2
      ;;
    --analysis=*)
      ANALYSIS_ARG="${1#*=}"
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
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Passed to workers via the environment (sbatch propagates the submit
# script's environment); 1.1_run_worker.sh forwards it to the R scripts.
export FORCE_BENCHMARK="${FORCE_ARG}"

# Resolve datasets (shared helper; skips `_*` keys unless --ds_name).
benchmark_resolve_datasets "${DS_NAME_ARG}"

ANALYSES=(trans zeroimp)
if [[ -n "${ANALYSIS_ARG}" ]]; then
  IFS=',' read -r -a ANALYSES <<< "${ANALYSIS_ARG}"
fi
UNIQUE_ANALYSES=()
for A in "${ANALYSES[@]}"; do
  if [[ "${A}" != "trans" && "${A}" != "zeroimp" ]]; then
    echo "ERROR: Unknown analysis '${A}' (expected trans or zeroimp)."
    exit 1
  fi
  if [[ ! " ${UNIQUE_ANALYSES[*]} " =~ " ${A} " ]]; then
    UNIQUE_ANALYSES+=("${A}")
  fi
done
ANALYSES=("${UNIQUE_ANALYSES[@]}")

PARTITION="${SLURM_PARTITION}"
EXTRA_FLAGS=()
if [[ -n "${PARTITION_ARG}" ]]; then
  PARTITION="${PARTITION_ARG}"
fi

mkdir -p "${LOGS_DIR}"

# ---------------------------------------------------------------------------
# Submit one array per analysis
# ---------------------------------------------------------------------------
echo "=== Submitting transformation/zero-imputation analysis arrays ==="

ARRAY_JOB_IDS=()
for ANALYSIS in "${ANALYSES[@]}"; do
  MANIFEST="${HPC_SCRATCH_DIR}/benchmark_manifest_${ANALYSIS}_$$.txt"
  : > "${MANIFEST}"
  for name in "${DATASET_NAMES[@]}"; do
    echo "${name}" >> "${MANIFEST}"
  done
  export BENCHMARK_MANIFEST="${MANIFEST}"
  # ANALYSIS must be exported too: sbatch propagates only the exported
  # environment, and 1.1_run_worker.sh hard-requires it.
  export ANALYSIS="${ANALYSIS}"

  echo "Submitting ${ANALYSIS} array (${NUM_DATASETS} datasets, partition=${PARTITION}, "
  echo "  flags: ${EXTRA_FLAGS[*]}, mem=32G, time=02:00:00, throttle=${MAX_NUM_CHUNKS_PARALLEL})"

  SUBMIT_MSG=$(sbatch \
      --array="1-${NUM_DATASETS}%${MAX_NUM_CHUNKS_PARALLEL}" \
      --partition="${PARTITION}" \
      "${EXTRA_FLAGS[@]}" \
      --cpus-per-task=4 \
      --mem=32G \
      --time=02:00:00 \
      --output="${LOGS_DIR}/5_transzeroimp_${ANALYSIS}_%A_%a.log" \
      --error="${LOGS_DIR}/5_transzeroimp_${ANALYSIS}_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      "${SCRIPT_DIR}/1.1_run_worker.sh")

  ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "  ${ANALYSIS} array job ID: ${ARRAY_JOB_ID}"
  ARRAY_JOB_IDS+=("${ARRAY_JOB_ID}")
done

# ---------------------------------------------------------------------------
# Monitor & verify & sync results back to NAS (shared tail)
# ---------------------------------------------------------------------------
for JOB_ID in "${ARRAY_JOB_IDS[@]}"; do
  benchmark_wait_for_array "${JOB_ID}" "analysis"
done

benchmark_merge_sync_cleanup "${ARRAY_JOB_IDS[@]}"
