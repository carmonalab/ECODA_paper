#!/bin/bash
#
# Shared helper functions for the benchmark SLURM submitters
# (src/5_run_benchmark_methods/run_python_sample_embedding_methods/,
#  run_r_sample_embedding_methods/, run_transformation_zeroimp_analysis/).
#
# Source AFTER slurm_config.sh and `cd "${PROJECT_ROOT}"` (functions use the
# exported slurm_config vars at call time):
#   source "$(dirname "${BASH_SOURCE[0]}")/benchmark_submit_common.sh"
#
# Provided:
#   benchmark_resolve_datasets <ds_name_arg>
#       Fills the global DATASET_NAMES/NUM_DATASETS (jq over datasets.json:
#       use_for_benchmark == true AND a benchmark_analysis view; `_*` keys
#       skipped unless explicitly requested). Errors out on unknown --ds_name
#       or zero datasets.
#   benchmark_wait_for_array <job_id> <label>
#       Blocks on `scontrol wait <job_id>` (event-driven; exit code ignored),
#       then polls sacct (bounded, 5s interval) until every state row is
#       terminal, then runs a fail-closed sacct gate (every state row must be
#       COMPLETED; aborts without syncing on any non-COMPLETED state or empty
#       sacct output).
#   benchmark_merge_sync_cleanup <labels...>
#       NAS reachability check FIRST (fail before any destructive merge
#       work), then writes the RDS integrity sidecar (benchmark/checksums.md5
#       over results/, pseudobulks/, gloscope_dists/ — verified by the
#       notebook's load_hpc_benchmark_results() before readRDS), merges the
#       per-task exec logs via 1.1.2_merge_execution_times.py (--no-cleanup,
#       --labels/--datasets-scoped over DATASET_NAMES, --existing-log NAS
#       continuity), rsyncs ${HPC_SCRATCH_DIR}/benchmark/ ->
#       ${NAS_TARGET_DIR}/benchmark/, and only then deletes this run's
#       per-task logs (scoped to the run's label x dataset cross product,
#       plus a separate legacy execution_times_task_* sweep).
#
# The three submitters also support `--sync-only <id1,id2,...>` resume mode
# (skip the sbatch submission loops, re-check the provided job ids via
# benchmark_wait_for_array, then run benchmark_merge_sync_cleanup): these
# helpers are reused as-is; the submitters only branch on the flag.
# Sync-status emails are sent by the helpers via notify_sync_status
# (src/utils/bash/sync_status_email.sh, sourced below) — one per gate
# failure ("NOT synced — reason") and one after a successful rsync.
# ============================================================================

# Path to the shared exec-log merge script, resolved from THIS file's location
# (BASH_SOURCE[0] inside a sourced file is the sourced file's path).
BENCHMARK_MERGE_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py"

# Sync-status email helper (best-effort; requires USER_EMAIL from slurm_config.sh).
source "$(dirname "${BASH_SOURCE[0]}")/../utils/bash/sync_status_email.sh"

# ---------------------------------------------------------------------------
# Dataset resolution (see header)
# ---------------------------------------------------------------------------
benchmark_resolve_datasets() {
  local DS_NAME_ARG="$1"
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
}

# ---------------------------------------------------------------------------
# Monitor an array to completion, then fail-closed sacct gate
# ---------------------------------------------------------------------------
benchmark_wait_for_array() {
  local JOB_ID="$1"
  local LABEL="$2"
  echo "=== Monitoring ${LABEL} array ${JOB_ID} ==="
  # Event-driven block until the job leaves the scheduler (no polling).
  # Exit code deliberately ignored: the fail-closed sacct gate below is the
  # authoritative check (covers cancellation, failure, purged controller
  # records).
  scontrol wait "${JOB_ID}" > /dev/null 2>&1 || true
  # sacct may lag a few seconds behind the job leaving the scheduler; poll
  # (bounded) until every state row is terminal instead of a blind fixed sleep.
  # The 180-iteration cap (15 min) plus a 60-iteration grace window (5 min)
  # covers pathological SlurmDBD accounting lag (scheduler said done, sacct
  # still reports RUNNING); the fail-closed gate below is unchanged.
  local TAIL_ITER=0
  local STATES=""
  while (( TAIL_ITER < 180 )); do  # max 15 min at 5s
    STATES="$(sacct -j "${JOB_ID}" --format=State -n 2>/dev/null || true)"
    if [[ -n "${STATES//[[:space:]]/}" ]] \
       && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
      break
    fi
    sleep 5
    TAIL_ITER=$((TAIL_ITER + 1))
  done
  if (( TAIL_ITER >= 180 )); then
    echo "WARNING: sacct still reporting non-terminal states after 15 min; extending wait by a 5 min grace window..." >&2
    TAIL_ITER=0
    while (( TAIL_ITER < 60 )); do
      STATES="$(sacct -j "${JOB_ID}" --format=State -n 2>/dev/null || true)"
      if [[ -n "${STATES//[[:space:]]/}" ]] \
         && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
        break
      fi
      sleep 5
      TAIL_ITER=$((TAIL_ITER + 1))
    done
  fi
  echo "${LABEL} array ${JOB_ID} finished. Checking task states..."
  STATES="$(sacct -j "${JOB_ID}" --format=State -n 2>/dev/null || true)"
  if [[ -z "${STATES//[[:space:]]/}" ]]; then
    echo "ERROR: sacct returned no states for Array Job ${JOB_ID}; NOT syncing to NAS."
    notify_sync_status \
      "ECODA: benchmark NOT synced (job ${JOB_ID})" \
      "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID} (datasets: ${DATASET_NAMES[*]}): sacct returned no states (job purged or unknown id)."
    exit 1
  fi
  # Fail-closed: every row (array master + tasks + batch steps) must be COMPLETED.
  if grep -qvE '^ *COMPLETED *$' <<< "${STATES}"; then
    echo "ERROR: Array Job ${JOB_ID} had non-COMPLETED tasks; NOT syncing to NAS."
    sacct -j "${JOB_ID}" --format=JobID,JobName,State,ExitCode
    notify_sync_status \
      "ECODA: benchmark NOT synced (job ${JOB_ID})" \
      "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID} (datasets: ${DATASET_NAMES[*]}): non-COMPLETED tasks.
$(sacct -j "${JOB_ID}" --format=JobID,JobName,State,ExitCode -n 2>/dev/null || true)"
    exit 1
  fi
  echo "Array Job ${JOB_ID} (${LABEL}): all tasks COMPLETED."
}

# ---------------------------------------------------------------------------
# NAS check -> RDS integrity sidecar -> merge exec logs -> rsync -> cleanup
# ---------------------------------------------------------------------------
benchmark_merge_sync_cleanup() {
  local LABELS=("$@")

  # NAS must be reachable BEFORE the merge: the merge with --no-cleanup keeps
  # the per-task logs until after the rsync, but a merge-then-fail would
  # otherwise leave the pipeline unable to sync anything without a --force
  # recompute.
  echo "Checking NAS reachability..."
  if ! ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
      echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable."
      notify_sync_status \
        "ECODA: benchmark NOT synced (no NAS access)" \
        "Benchmark sync to NAS skipped (datasets: ${DATASET_NAMES[*]}, labels: ${LABELS[*]}): NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
      exit 1
  fi
  mkdir -p "${NAS_TARGET_DIR}/benchmark"

  # Write an md5 checksum sidecar over the RDS result bundles so the notebook
  # (load_hpc_benchmark_results) can verify NAS-loaded bundles before
  # deserializing them. Best-effort: only if at least one RDS exists.
  if ls "${HPC_SCRATCH_DIR}/benchmark"/results/*.rds \
         "${HPC_SCRATCH_DIR}/benchmark"/pseudobulks/*.rds \
         "${HPC_SCRATCH_DIR}/benchmark"/gloscope_dists/*.rds > /dev/null 2>&1; then
    FIND_DIRS=()
    for d in results pseudobulks gloscope_dists; do
      [[ -d "${HPC_SCRATCH_DIR}/benchmark/${d}" ]] && FIND_DIRS+=("${d}")
    done
    (cd "${HPC_SCRATCH_DIR}/benchmark" && \
       find "${FIND_DIRS[@]}" -type f -name '*.rds' -exec md5sum {} + > checksums.md5)
    echo "Wrote benchmark/checksums.md5 (RDS bundle integrity sidecar)."
  fi

  echo "All tasks completed successfully. Merging execution-time logs..."
  # --no-cleanup: per-task logs are deleted only AFTER the rsync below
  # succeeds. --labels x --datasets scopes the merge to THIS run's
  # (method/analysis x dataset) cross product so stale logs from previous
  # failed runs never leak in. --existing-log preserves the NAS log across
  # partial (e.g. --ds_name _debug) runs instead of overwriting it with
  # subset rows. DATASET_NAMES is filled by benchmark_resolve_datasets
  # (called before this function by every submitter).
  local MERGE_ARGS=(--output_dir "${HPC_SCRATCH_DIR}/benchmark/embeddings"
                    --no-cleanup
                    --labels "${LABELS[@]}"
                    --datasets "${DATASET_NAMES[@]}")
  # The rsync below copies ${HPC_SCRATCH_DIR}/benchmark/ wholesale, so the
  # merged log lives at benchmark/embeddings/execution_times.feather.
  local EXISTING_LOG="${NAS_TARGET_DIR}/benchmark/embeddings/execution_times.feather"
  if [[ -f "${EXISTING_LOG}" ]]; then
      MERGE_ARGS+=(--existing-log "${EXISTING_LOG}")
  fi
  "${PYTHON_BIN}" "${BENCHMARK_MERGE_SCRIPT}" "${MERGE_ARGS[@]}"

  echo "Merged logs. Syncing results to NAS..."
  rsync -rlptDv "${HPC_SCRATCH_DIR}/benchmark/" "${NAS_TARGET_DIR}/benchmark/"
  echo "Results synchronized to ${NAS_TARGET_DIR}/benchmark/"
  notify_sync_status \
    "ECODA: benchmark synced to NAS" \
    "Benchmark results synced to ${NAS_TARGET_DIR}/benchmark/ (datasets: ${DATASET_NAMES[*]}, labels: ${LABELS[*]})."

  # Per-task logs may be deleted only now that the sync has succeeded.
  # Scoped to THIS run's (label x dataset) cross product so an overlapping
  # submission's not-yet-merged logs are never deleted; the legacy
  # execution_times_task_* sweep is a separate glob (no current worker
  # produces that naming, so it can only hit stale files). Neither glob
  # matches the merged execution_times.feather (no "_" suffix).
  for LABEL in "${LABELS[@]}"; do
    for DS in "${DATASET_NAMES[@]}"; do
      rm -f "${HPC_SCRATCH_DIR}/benchmark/embeddings"/"execution_times_${LABEL}_${DS}.feather"
    done
  done
  rm -f "${HPC_SCRATCH_DIR}/benchmark/embeddings"/execution_times_task_*.feather
  echo "Deleted per-task execution-time logs."
}
