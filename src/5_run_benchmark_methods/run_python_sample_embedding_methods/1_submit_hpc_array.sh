#!/bin/bash
set -euo pipefail

# Submit Python benchmark method arrays (MrVI/scPoli on GPU, PILOT on CPU)
#
# Usage:
#   ./1_submit_hpc_array.sh                                    # all methods, all benchmark datasets
#   ./1_submit_hpc_array.sh --ds_name _debug --methods mrvi    # single dataset, single method
#   ./1_submit_hpc_array.sh --methods mrvi,scpoli --force      # recompute existing feathers
#   ./1_submit_hpc_array.sh --partition debug-cpu              # override per-method partitions
#
# One SLURM array per method; array task IDs map 1:1 to lines of the
# per-method manifest ${HPC_SCRATCH_DIR}/benchmark_manifest_<method>_<pid>.txt
# (written by this script, one dataset per line; PID suffix so overlapping
# submissions do not clobber queued arrays' manifests). Hardware is PINNED for
# runtime comparability (see slurm_config.sh benchmark vars): GPU methods
# (mrvi, scpoli) on ${SLURM_PARTITION_BENCHMARK_GPU} with
# ${BENCHMARK_GPU_CONSTRAINT}, PILOT on ${SLURM_PARTITION_BENCHMARK_CPU} with
# ${BENCHMARK_CPU_CONSTRAINT}. Flags are passed on the sbatch command line
# (SLURM directives do not expand env vars). After all arrays complete:
# NAS reachability check -> fail-closed sacct gate -> merge per-task exec
# logs (job-id scoped, existing-log continuity) -> sync to NAS -> only then
# delete this run's per-task logs.

source "$(dirname "${BASH_SOURCE[0]}")/../../slurm_config.sh"
cd "${PROJECT_ROOT}"

module load jq/1.6

DS_NAME_ARG=""
METHODS_ARG=""
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
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Passed to workers via the environment (sbatch propagates the submit
# script's environment); 1.1_run_worker.sh forwards it to the py script.
export FORCE_BENCHMARK="${FORCE_ARG}"

# ---------------------------------------------------------------------------
# Resolve datasets: use_for_benchmark == true AND a benchmark_analysis view.
# Convention: default-all skips keys starting with "_" (e.g. _debug) unless
# explicitly requested via --ds_name.
# ---------------------------------------------------------------------------
DATASET_NAMES=()
if [[ -n "${DS_NAME_ARG}" ]]; then
  if ! jq -e --arg ds "${DS_NAME_ARG}" 'has($ds)' "${DATASETS_JSON_FILE}" > /dev/null 2>&1; then
    echo "ERROR: '${DS_NAME_ARG}' is not a dataset in ${DATASETS_JSON_FILE}."
    exit 1
  fi
  DATASET_NAMES+=("${DS_NAME_ARG}")
else
  while IFS= read -r name; do
    DATASET_NAMES+=("$name")
  done < <(jq -r 'to_entries[] |
    select(.value.use_for_benchmark == true) |
    select(.value.views.benchmark_analysis != null) |
    .key | select(startswith("_") | not)' "${DATASETS_JSON_FILE}")
fi

NUM_DATASETS=${#DATASET_NAMES[@]}
if [[ ${NUM_DATASETS} -eq 0 ]]; then
  echo "ERROR: No benchmark datasets found in ${DATASETS_JSON_FILE}."
  exit 1
fi

echo "Found ${NUM_DATASETS} benchmark datasets."

METHODS=(mrvi scpoli pilot)
if [[ -n "${METHODS_ARG}" ]]; then
  IFS=',' read -r -a METHODS <<< "${METHODS_ARG}"
fi

# ---------------------------------------------------------------------------
# Submit one array per method
# ---------------------------------------------------------------------------
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
      "$(dirname "${BASH_SOURCE[0]}")/1.1_run_worker.sh")

  ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "  ${METHOD} array job ID: ${ARRAY_JOB_ID}"
  ARRAY_JOB_IDS+=("${ARRAY_JOB_ID}")
done

# ---------------------------------------------------------------------------
# Monitor & verify & sync results back to NAS
# ---------------------------------------------------------------------------
echo "=== Monitoring job completion ==="
for JOB_ID in "${ARRAY_JOB_IDS[@]}"; do
  while squeue -u "$USER" 2>/dev/null | grep -q "${JOB_ID}"; do
    sleep 60
  done
done

# Give sacct a moment to record final states (mirrors 1_submit_hpc.sh)
sleep 30

echo "Arrays finished. Checking task states..."
for JOB_ID in "${ARRAY_JOB_IDS[@]}"; do
  STATES="$(sacct -j "${JOB_ID}" --format=State -n 2>/dev/null || true)"
  if [[ -z "${STATES//[[:space:]]/}" ]]; then
    echo "ERROR: sacct returned no states for Array Job ${JOB_ID}; NOT syncing to NAS."
    exit 1
  fi
  # Fail-closed: every row (array master + tasks + batch steps) must be COMPLETED.
  if grep -qvE '^ *COMPLETED *$' <<< "${STATES}"; then
    echo "ERROR: Array Job ${JOB_ID} had non-COMPLETED tasks; NOT syncing to NAS."
    sacct -j "${JOB_ID}" --format=JobID,JobName,State,ExitCode
    exit 1
  fi
  echo "Array Job ${JOB_ID}: all tasks COMPLETED."
done

# NAS must be reachable BEFORE the merge: the merge with --no-cleanup keeps
# the per-task logs until after the rsync, but a merge-then-fail would
# otherwise leave the pipeline unable to sync anything without a --force
# recompute.
echo "Checking NAS reachability..."
if ! ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable."
    exit 1
fi
mkdir -p "${NAS_TARGET_DIR}/benchmark"

echo "All tasks completed successfully. Merging execution-time logs..."
# --no-cleanup: per-task logs are deleted only AFTER the rsync below
# succeeds. --job_ids scopes the merge to THIS run's arrays (task ids are
# per-array, job ids are unique) so stale logs from previous failed runs
# never leak in. --existing-log preserves the NAS log across partial
# (e.g. --ds_name _debug) runs instead of overwriting it with subset rows.
MERGE_ARGS=(--output_dir "${HPC_SCRATCH_DIR}/benchmark/embeddings"
            --no-cleanup
            --job_ids "${ARRAY_JOB_IDS[@]}")
EXISTING_LOG="${NAS_TARGET_DIR}/benchmark/execution_times.feather"
if [[ -f "${EXISTING_LOG}" ]]; then
    MERGE_ARGS+=(--existing-log "${EXISTING_LOG}")
fi
"${PYTHON_BIN}" \
  "$(dirname "${BASH_SOURCE[0]}")/1.1.2_merge_execution_times.py" \
  "${MERGE_ARGS[@]}"

echo "Merged logs. Syncing results to NAS..."
rsync -rlptDv "${HPC_SCRATCH_DIR}/benchmark/" "${NAS_TARGET_DIR}/benchmark/"
echo "Results synchronized to ${NAS_TARGET_DIR}/benchmark/"

# Per-task logs may be deleted only now that the sync has succeeded.
for JOB_ID in "${ARRAY_JOB_IDS[@]}"; do
    rm -f "${HPC_SCRATCH_DIR}/benchmark/embeddings"/execution_times_task_${JOB_ID}_*.feather
done
echo "Deleted per-task execution-time logs for job ids: ${ARRAY_JOB_IDS[*]}"
