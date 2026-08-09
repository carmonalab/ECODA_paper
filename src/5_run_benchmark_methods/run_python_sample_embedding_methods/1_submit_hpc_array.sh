#!/bin/bash
set -euo pipefail

# Submit Python benchmark method arrays (MrVI/scPoli on GPU, PILOT on CPU)
#
# Usage:
#   ./1_submit_hpc_array.sh                                    # all methods, all benchmark datasets
#   ./1_submit_hpc_array.sh --ds_name _debug --methods mrvi    # single dataset, single method
#   ./1_submit_hpc_array.sh --methods mrvi,scpoli --force      # recompute existing feathers
#   ./1_submit_hpc_array.sh --partition debug-cpu              # override per-method partitions (drops the constraint pin)
#   ./1_submit_hpc_array.sh --sync-only 12345,12346            # resume: skip submission, re-check + sync
#                                                              # (repeat the original --ds_name/--methods flags)
#
# One SLURM array per method; array task IDs map 1:1 to lines of the
# per-method manifest ${HPC_SCRATCH_DIR}/benchmark_manifest_<method>_<pid>.txt
# (written by this script, one dataset per line; PID suffix so overlapping
# submissions do not clobber queued arrays' manifests). Hardware is PINNED for
# runtime comparability (see slurm_config.sh benchmark vars): GPU methods
# (mrvi, scpoli) on ${SLURM_PARTITION_BENCHMARK_GPU} with
# ${BENCHMARK_GPU_CONSTRAINT}, PILOT on ${SLURM_PARTITION_BENCHMARK_CPU} with
# ${BENCHMARK_CPU_CONSTRAINT}. Flags are passed on the sbatch command line
# (SLURM directives do not expand env vars). An explicit --partition <P>
# override (e.g. --partition ${SLURM_PARTITION_PRIVATE} for _debug runs on the
# private node, or an ad-hoc --partition shared-gpu for debug) DROPS the
# method's --constraint flag — an explicit partition
# choice means the user accepts non-pinned hardware; keeping the constraint
# would otherwise hang jobs PENDING forever on nodes whose CPU/GPU differ.
# After all arrays complete:
# NAS reachability check -> fail-closed sacct gate -> merge per-task exec
# logs (job-id scoped, existing-log continuity) -> sync to NAS -> only then
# delete this run's per-task logs.

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

METHODS=(mrvi scpoli pilot)
if [[ -n "${METHODS_ARG}" ]]; then
  IFS=',' read -r -a METHODS <<< "${METHODS_ARG}"
fi

# ---------------------------------------------------------------------------
# Submit one array per method
# ---------------------------------------------------------------------------
if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Sync-only resume mode: jobs ${SYNC_ONLY_IDS} (no submission) ==="
  IFS=',' read -r -a ARRAY_JOB_IDS <<< "${SYNC_ONLY_IDS}"
else
echo "=== Submitting Python benchmark method arrays ==="

mkdir -p "${LOGS_DIR}"
ARRAY_JOB_IDS=()

for METHOD in "${METHODS[@]}"; do
  case "${METHOD}" in
    mrvi|scpoli)
      PARTITION="${SLURM_PARTITION_BENCHMARK_GPU}"
      EXTRA_FLAGS=(--gpus="${BENCHMARK_GPU_COUNT}"
                   --constraint="${BENCHMARK_GPU_CONSTRAINT}"
                   --cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}")
      THROTTLE="${BENCHMARK_GPU_ARRAY_THROTTLE}"
      ;;
    pilot)
      PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
      EXTRA_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}"
                   --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
      THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
      ;;
    *)
      echo "ERROR: Unknown method '${METHOD}' (expected mrvi, scpoli or pilot)."
      exit 1
      ;;
  esac
  if [[ -n "${PARTITION_ARG}" ]]; then
    PARTITION="${PARTITION_ARG}"
    # Explicit partition choice = user accepts non-pinned hardware: drop the
    # --constraint pin (kept: --gpus/--cpus-per-task/--mem). Without this,
    # e.g. --partition private-carmona-gpu would sit PENDING forever — its
    # CPU/GPU never match the pinned BENCHMARK_*_CONSTRAINT.
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

  # Per-method manifest: one dataset per line; rebuilt every run. The name
  # carries the submit PID so an overlapping second submission cannot
  # clobber the manifest of already-submitted (queued) arrays.
  MANIFEST="${HPC_SCRATCH_DIR}/benchmark_manifest_${METHOD}_$$.txt"
  : > "${MANIFEST}"
  for name in "${DATASET_NAMES[@]}"; do
    echo "${name}" >> "${MANIFEST}"
  done
  export BENCHMARK_MANIFEST="${MANIFEST}"
  # METHOD must be exported too: sbatch propagates only the exported
  # environment, and 1.1_run_worker.sh hard-requires it.
  export METHOD="${METHOD}"

  echo "Submitting ${METHOD} array (${NUM_DATASETS} datasets, partition=${PARTITION}, "
  echo "  flags: ${EXTRA_FLAGS[*]}, mem=${BENCHMARK_MEM}, throttle=${THROTTLE})"

  SUBMIT_MSG=$(sbatch \
      --array="1-${NUM_DATASETS}%${THROTTLE}" \
      --partition="${PARTITION}" \
      "${EXTRA_FLAGS[@]}" \
      --mem="${BENCHMARK_MEM}" \
      --output="${LOGS_DIR}/5_benchmark_${METHOD}_%A_%a.log" \
      --error="${LOGS_DIR}/5_benchmark_${METHOD}_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      "${SCRIPT_DIR}/1.1_run_worker.sh")

  ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "  ${METHOD} array job ID: ${ARRAY_JOB_ID}"
  ARRAY_JOB_IDS+=("${ARRAY_JOB_ID}")
done
fi

# ---------------------------------------------------------------------------
# Monitor & verify & sync results back to NAS (shared tail)
# ---------------------------------------------------------------------------
echo "=== Monitoring job completion ==="
for JOB_ID in "${ARRAY_JOB_IDS[@]}"; do
  benchmark_wait_for_array "${JOB_ID}" "benchmark"
done

# Labels for the exec-log merge = the submitted methods.
benchmark_merge_sync_cleanup "${METHODS[@]}"
