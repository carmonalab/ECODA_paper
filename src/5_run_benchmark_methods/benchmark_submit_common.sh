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
#       Blocks on an exact-id `squeue -j` poll (60s interval; `scontrol` has
#       no plain wait command), then waits for sacct rows to become terminal.
#       The watchdog time limit is the outer bound; two consecutive empty
#       sacct responses fail closed. Every state row must be COMPLETED before
#       syncing.
#   benchmark_wait_oom_retry <job_id> <label> <resubmit_fn> <manifest> [status_file]
#       OOM-auto-escalating variant of the wait+gate for the benchmark
#       submitters' own arrays: waits for the array to leave the scheduler
#       (shared benchmark_wait_array_terminal poll), then gates on per-TASK
#       states only (`<job_id>_<n>` rows; .batch/.extern/master rows excluded).
#       All-COMPLETED -> JOB_REPORTS record + return 0. Any non-COMPLETED,
#       non-OUT_OF_MEMORY state -> fail closed (task report + sync-status
#       email, exit 1, no sync). OUT_OF_MEMORY tasks -> re-submit ONLY those
#       tasks' datasets (mapped from <manifest> via sed -n <task>p) with
#       doubled --mem (benchmark_bump_mem, CLAMPED to the BENCHMARK_MEM_MAX
#       ceiling so a retry never exceeds the nodes' RAM), via
#       ${resubmit_fn} <label> <comma-separated datasets> <new
#       mem> <new manifest path> — which must write the manifest, sbatch the
#       reduced array (same flags/throttle as the normal submission) and echo
#       ONLY the new array job id on stdout. Loop back with the new id/manifest;
#       belt-and-braces attempt cap (4). At the ceiling or on non-memory
#       failures: fail closed with an OOM report incl. per-task MaxRSS.
#       Optional 5th arg <status_file>: watchdog mode — every terminal path
#       writes the status file (STATE=OK|FAIL, LABEL=, one JOB_REPORT= line per
#       gated array incl. the final retry id, FAIL_REASON=, REPORT=) instead of
#       emailing (compute nodes have no mail CLI; the login tail emails from
#       the file via benchmark_wait_watchdog) and instead of appending to the
#       JOB_REPORTS global (the tail merges the JOB_REPORT= lines from the
#       file). Behavior without the arg is unchanged.
#   benchmark_submit_watchdog <array_id> <label> <manifest> <mode> <partition>
#       <throttle> <log_prefix> <worker_script> [flags...]
#       Submits one compute-node watchdog job (watchdog_main.sh, 1 cpu/2G/
#       WATCHDOG_TIME_LIMIT, default partition of the method — no constraint
#       pin) that owns the terminal wait + OOM escalation for <array_id> and
#       writes its status file `${WATCHDOG_STATUS_DIR}/<watchdog_job_id>.status`
#       (self-named from SLURM_JOB_ID; unknowable at submit time). <mode> is
#       strict (method arrays) or soft-gate (prepare_pseudobulk artifact gate).
#       <throttle>/<log_prefix>/<worker_script>/<flags> are forwarded for the
#       watchdog's retry-array submissions. Echoes ONLY the watchdog job id.
#   benchmark_wait_watchdog <watchdog_id> <label>
#       Login-tail counterpart: waits for the watchdog job to leave the
#       scheduler, polls for its status file (<=2 min grace), then parses
#       STATE: OK -> merge its JOB_REPORT= lines into JOB_REPORTS, return 0;
#       FAIL -> print the report, notify_sync_status, exit 1; watchdog job
#       non-COMPLETED without a status file -> fail closed with its sacct
#       State,ExitCode + a pointer to its logs; COMPLETED without a status
#       file -> fail closed ("exited without a status file").
#   benchmark_bump_mem <mem> / benchmark_mem_ge <a> <b>
#       Mem-string helpers for the OOM escalation: <N>G/<N>T -> 2N (same
#       suffix; non-zero exit on unparseable input) / truthy when a >= b.
#   benchmark_merge_sync_cleanup <labels...>
#       NAS reachability check FIRST (fail before any destructive merge
#       work), merges execution logs atomically, writes a scope-specific
#       checksum sidecar over final RDS/Feather outputs, rsyncs the selected
#       canonical root, verifies the remote checksum, and only then deletes
#       this run's per-task logs.
#
# The canonical 1_submit_hpc_array.sh owns the matrix submission, watchdog
# aggregation, and sync-only RUN_ID resume path. The family entrypoints are
# compatibility delegators; these helpers remain available to the watchdog
# regression tests and any narrowly scoped legacy worker.
# Sync-status emails are sent by the helpers via notify_sync_status
# (src/utils/bash/sync_status_email.sh, sourced below) — one per gate
# failure ("NOT synced — reason") and one after a successful rsync. Gate
# failures carry a per-task report (task -> DATASET_NAMES[i-1], state,
# elapsed, exit code) + array wall time; the final success/NAS-unreachable
# emails carry a "Job durations" block (label, job id, array wall time) for
# every array gated during this run, accumulated in the global JOB_REPORTS
# by benchmark_wait_for_array.
# ============================================================================

# Path to the shared exec-log merge script, resolved from THIS file's location
# (BASH_SOURCE[0] inside a sourced file is the sourced file's path).
ANALYSIS_MERGE_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py"

# Compute-node watchdog entry script (same directory as this file).
WATCHDOG_MAIN_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/watchdog_main.sh"

# Watchdog status files: one `<watchdog_job_id>.status` per watchdog, written
# by watchdog_main.sh (self-named from SLURM_JOB_ID) and read by the login
# tail's benchmark_wait_watchdog. Outside benchmark/ so it is never rsync'd.
if [[ -n "${ANALYSIS_PASS:-}" ]]; then
  WATCHDOG_STATUS_DIR="${HPC_SCRATCH_DIR}/_batch_effect_watchdog/${ANALYSIS_PASS}"
else
  WATCHDOG_STATUS_DIR="${HPC_SCRATCH_DIR}/_benchmark_watchdog"
fi

# Per-array job duration records for the final email: one "<label>|<job id>|<array wall time>"
# entry per gated array, appended by benchmark_wait_for_array /
# benchmark_wait_oom_retry in submission order (only the FINAL, successful
# retry of an OOM-escalated array is recorded — intermediate OOM'd attempts
# are never appended).
JOB_REPORTS=()

# Same records, collected for the status file in watchdog mode: every array
# gated by the watchdog run (original + each retry) so the login tail's final
# email shows the whole escalation chain. Only used when benchmark_wait_oom_retry
# is called with a status_file arg; the tail merges these into JOB_REPORTS.
WATCHDOG_GATED_REPORTS=()
WATCHDOG_SCHEDULER_IDS=()

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

# Validate authoritative processed h5ads before submitting any benchmark
# workers. This is intentionally independent of result-file existence: a
# metadata/PCA-only local mirror must never be accepted as a worker input.
benchmark_validate_h5ad_inputs() {
  local VIEW="$1"
  local METHOD="$2"
  shift 2
  local DS_NAME OUTPUT_FILE H5AD_PATH
  for DS_NAME in "$@"; do
    OUTPUT_FILE="$(jq -r --arg ds "${DS_NAME}" --arg view "${VIEW}" \
      '.[$ds].views[$view].output_file_name // empty' \
      "${DATASETS_JSON_FILE}")"
    if [[ -z "${OUTPUT_FILE}" ]]; then
      echo "ERROR: no output_file_name for ${DS_NAME}/${VIEW} in ${DATASETS_JSON_FILE}."
      return 1
    fi
    H5AD_PATH="${HPC_SCRATCH_DIR}/${DS_NAME}/output/${OUTPUT_FILE}"
    echo "Validating h5ad contract: ${H5AD_PATH}" >&2
    "${PYTHON_BIN}" \
      "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
      --path "${H5AD_PATH}" \
      --view "${VIEW}" \
      --method "${METHOD}" || return 1
  done
}

# ---------------------------------------------------------------------------
# Per-task email report for a gated array (task -> DATASET_NAMES[i-1], state,
# elapsed, exit code) + array wall time. Prints an n/a line when sacct is
# unavailable or purged. DATASET_NAMES must be filled (benchmark_resolve_datasets
# runs before every wait loop, in submission AND --sync-only modes, so the
# mapping is valid in both). An optional second argument <manifest> maps task
# ids from that file instead (used by the OOM retry loop, whose retried
# arrays' manifests hold only the re-submitted datasets).
# ---------------------------------------------------------------------------
benchmark_task_report() {
  local JOB_ID="$1"
  local MANIFEST="${2:-}"
  local TASKS
  TASKS="$(array_task_report "${JOB_ID}")"
  if [[ -z "${TASKS}" ]]; then
    printf '%s' "Per-task report: n/a (sacct unavailable or job purged)."
    return 0
  fi
  printf '%s\n' "Per-task report (task -> dataset, state, elapsed, exit code):"
  local task state elapsed exitcode extra ds_name
  while IFS=$'\t' read -r task state elapsed exitcode; do
    ds_name="?"
    if [[ -n "${MANIFEST}" && -f "${MANIFEST}" ]]; then
      ds_name="$(sed -n "${task}p" "${MANIFEST}")"
    elif [[ ${task} -ge 1 && ${task} -le ${NUM_DATASETS} ]]; then
      ds_name="${DATASET_NAMES[$((task - 1))]}"
    fi
    extra=""
    [[ "${state}" == "COMPLETED" ]] || extra=" (${exitcode})"
    printf '  %s -> %-30s %-14s%s  %s\n' "${task}" "${ds_name}" "${state}" "${extra}" "${elapsed}"
  done <<< "${TASKS}"
  printf 'Array wall time: %s\n' "$(array_wall_time "${JOB_ID}")"
}

# ---------------------------------------------------------------------------
# "Job durations" block for the final emails: label, job id, array wall time
# for every array gated during this run (JOB_REPORTS, in submission order).
# ---------------------------------------------------------------------------
benchmark_job_durations_block() {
  if (( ${#JOB_REPORTS[@]} == 0 )); then
    printf '%s' "Job durations: n/a (no gated arrays)."
    return 0
  fi
  printf '%s\n' "Job durations (label, job id, array wall time):"
  local line label jid wall
  for line in "${JOB_REPORTS[@]}"; do
    IFS='|' read -r label jid wall <<< "${line}"
    printf '  %-30s %s  %s\n' "${label}" "${jid}" "${wall}"
  done
}

# ---------------------------------------------------------------------------
# Memory escalation helpers (OOM retry)
# ---------------------------------------------------------------------------

# benchmark_bump_mem <mem>: parse <N>G/<N>T, echo 2N with the same suffix.
# Non-zero exit on unparseable input (callers fail closed).
benchmark_bump_mem() {
  local MEM="$1"
  if [[ "${MEM}" =~ ^([0-9]+)([GT])$ ]]; then
    printf '%d%s' $(( ${BASH_REMATCH[1]} * 2 )) "${BASH_REMATCH[2]}"
    return 0
  fi
  return 1
}

# benchmark_mem_ge <a> <b>: exit 0 when a >= b (G/T suffix aware), 1 otherwise
# or on unparseable input.
benchmark_mem_ge() {
  local A="$1" B="$2"
  local a_num b_num a_suf b_suf
  [[ "${A}" =~ ^([0-9]+)([GT])$ ]] || return 1
  a_num=${BASH_REMATCH[1]}
  a_suf=${BASH_REMATCH[2]}
  [[ "${B}" =~ ^([0-9]+)([GT])$ ]] || return 1
  b_num=${BASH_REMATCH[1]}
  b_suf=${BASH_REMATCH[2]}
  if [[ "${a_suf}" == "${b_suf}" ]]; then
    (( a_num >= b_num ))
  elif [[ "${a_suf}" == "T" ]]; then
    (( a_num * 1024 >= b_num ))
  else
    (( a_num >= b_num * 1024 ))
  fi
}

# ---------------------------------------------------------------------------
# PB_VARIANT_NAMES parse (benchmark_hpc_utils.R is the single source of
# truth): extracts the c("...", ...) list via sed range + grep, one variant
# per line on stdout; prints nothing (non-zero-ish empty output, callers fail
# closed) when the list cannot be parsed. Used by the watchdog's soft-gate
# mode (prepare_pseudobulk artifact gate).
# ---------------------------------------------------------------------------
benchmark_pb_variant_names() {
  local HPC_UTILS="$1"
  sed -n '/^PB_VARIANT_NAMES <- c(/,/^)/p' \
    "${HPC_UTILS}" | grep -oE '"[a-zA-Z0-9_.]+"' | tr -d '"'
}

# ---------------------------------------------------------------------------
 # Shared array terminal wait: exact-id squeue poll (60s) + sacct
 # poll-until-terminal. Used by benchmark_wait_for_array (all states) and
 # benchmark_wait_oom_retry (per-task states) — behavior of the wait itself is
 # identical; the gate semantics differ per caller.
# ---------------------------------------------------------------------------
benchmark_wait_array_terminal() {
  local JOB_ID="$1"
  local LABEL="$2"
  echo "=== Monitoring ${LABEL} array ${JOB_ID} ==="
  # Block until the job leaves the scheduler. scontrol has NO plain `wait`
  # command (only `wait_job`, which waits for node-ready — not completion —
  # and is documented as unusable with SLURM_ARRAY_JOB_ID), so query squeue
  # for this array ID directly. A user-wide query can omit a still-running
  # array task under scheduler/permission lag and falsely release the gate.
  # The fail-closed gates downstream are the authoritative check (covers
  # cancellation, failure, purged controller records).
  # A transient squeue failure (non-zero exit / non-empty stderr) or a clean
  # response still containing the id resets the miss counter. Only 2
  # CONSECUTIVE clean responses without the id count as "left the scheduler".
  # This avoids a transient empty/error squeue response (stale BeeGFS view,
  # scheduler hiccup) being misread as an early exit; the bounded sacct poll
  # below is the final safety net, not the primary signal. After 3 failures,
  # a terminal sacct response is an explicit fallback; failure from both
  # commands still stops the wait closed.
  local MISS=0 SQUEUE_FAILURES=0 OUT=""
  while :; do
    if ! OUT="$(timeout 30 squeue -j "${JOB_ID}" -h -o "%A" 2>/dev/null)"; then
      SQUEUE_FAILURES=$((SQUEUE_FAILURES + 1))
      if (( SQUEUE_FAILURES >= 3 )); then
        local FALLBACK_STATES
        FALLBACK_STATES="$(timeout 30 sacct -j "${JOB_ID}" --format=State -n 2>/dev/null || true)"
        if [[ -z "${FALLBACK_STATES//[[:space:]]/}" ]]; then
          echo "ERROR: squeue and sacct could not report terminal state for ${JOB_ID}; stopping wait." >&2
          return 1
        elif ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${FALLBACK_STATES}"; then
          echo "WARNING: squeue unavailable for ${JOB_ID}; sacct reports terminal state, continuing." >&2
          break
        fi
        echo "WARNING: squeue unavailable for ${JOB_ID}; sacct still reports a non-terminal state." >&2
      else
        echo "WARNING: squeue query failed while waiting for ${JOB_ID}; keeping polling." >&2
      fi
      MISS=0
    elif grep -qx "${JOB_ID}" <<< "${OUT}"; then
      SQUEUE_FAILURES=0
      MISS=0
    else
      SQUEUE_FAILURES=0
      MISS=$((MISS + 1))
    fi
    (( MISS >= 2 )) && break
    sleep 60
  done
  echo "${LABEL} array ${JOB_ID} left the scheduler."
  # sacct may lag a few seconds after the scheduler query stops returning the
  # array. Wait until every accounting row is terminal; do not impose a short
  # fixed cap because a legitimate large-dataset task can run longer than the
  # old 20-minute grace window. The watchdog's Slurm time limit remains the
  # outer fail-closed ceiling. Two consecutive empty accounting responses are
  # still treated as an accounting failure and released to the caller's
  # fail-closed check.
  echo "Waiting for sacct to record terminal states for job ${JOB_ID}..."
  local EMPTY=0 STATES="" EMPTY_GRACE="${ECODA_ACCOUNTING_EMPTY_GRACE:-3}"
  while :; do
    STATES="$(timeout 30 sacct -j "${JOB_ID}" --format=State -n 2>/dev/null || true)"
    if [[ -z "${STATES//[[:space:]]/}" ]]; then
      EMPTY=$((EMPTY + 1))
      if (( EMPTY >= EMPTY_GRACE )); then
        echo "WARNING: sacct returned no states ${EMPTY_GRACE} times for ${JOB_ID}; stopping wait." >&2
        break
      fi
    elif ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
      break
    else
      EMPTY=0
    fi
    sleep 5
  done
  echo "${LABEL} array ${JOB_ID} finished. Checking task states..."
}

# ---------------------------------------------------------------------------
# OOM per-task email report: task -> dataset (mapped from the run's manifest
# when given, else DATASET_NAMES), state, MaxRSS, elapsed + array wall time.
# Prints an n/a line when sacct is unavailable or purged.
# ---------------------------------------------------------------------------
benchmark_oom_task_report() {
  local JOB_ID="$1"
  local MANIFEST="${2:-}"
  local ROWS
  ROWS="$(sacct -j "${JOB_ID}" -n --parsable2 --format=JobID,State,MaxRSS,Elapsed 2>/dev/null || true)"
  if [[ -z "${ROWS//[[:space:]]/}" ]]; then
    printf '%s' "OOM per-task report: n/a (sacct unavailable or job purged)."
    return 0
  fi
  printf '%s\n' "OOM per-task report (task -> dataset, state, MaxRSS, elapsed):"
  local jid state maxrss elapsed task ds_name
  while IFS='|' read -r jid state maxrss elapsed; do
    jid="$(tr -d '[:space:]' <<< "${jid}")"
    [[ "${jid}" =~ ^${JOB_ID}_[0-9]+$ ]] || continue
    task="${jid#${JOB_ID}_}"
    ds_name="?"
    if [[ -n "${MANIFEST}" && -f "${MANIFEST}" ]]; then
      ds_name="$(sed -n "${task}p" "${MANIFEST}")"
    elif [[ ${task} -ge 1 && ${task} -le ${NUM_DATASETS} ]]; then
      ds_name="${DATASET_NAMES[$((task - 1))]}"
    fi
    printf '  %s -> %-30s %-14s %-10s  %s\n' "${task}" "${ds_name}" "$(tr -d '[:space:]' <<< "${state}")" "${maxrss}" "${elapsed}"
  done <<< "${ROWS}"
  printf 'Array wall time: %s\n' "$(array_wall_time "${JOB_ID}")"
}

# ---------------------------------------------------------------------------
# Watchdog status file writer (watchdog mode of benchmark_wait_oom_retry and
# the watchdog's soft-gate OK path): atomic write (tmp + mv) of
#   STATE=OK|FAIL
#   LABEL=<label>
#   JOB_REPORT=<label>|<id>|<wall>        (one per gated array, from
#                                          WATCHDOG_GATED_REPORTS)
#   FAIL_REASON=<reason>                  (FAIL only)
#   REPORT=<multiline>                    (FAIL only; per-task report)
# ---------------------------------------------------------------------------
benchmark_write_status_file() {
  local STATUS_FILE="$1"
  local STATE="$2"
  local LABEL="$3"
  local FAIL_REASON="${4:-}"
  local REPORT="${5:-}"
  local TMP_FILE="${STATUS_FILE}.tmp"
  mkdir -p "$(dirname "${STATUS_FILE}")"
  {
    printf 'STATE=%s\n' "${STATE}"
    printf 'LABEL=%s\n' "${LABEL}"
    local scheduler_id
    for scheduler_id in "${WATCHDOG_SCHEDULER_IDS[@]}"; do
      printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"
    done
    local line
    for line in "${WATCHDOG_GATED_REPORTS[@]}"; do
      printf 'JOB_REPORT=%s\n' "${line}"
    done
    if [[ -n "${FAIL_REASON}" ]]; then
      printf 'FAIL_REASON=%s\n' "${FAIL_REASON}"
      printf 'REPORT=\n%s\n' "${REPORT}"
    fi
  } > "${TMP_FILE}"
  mv "${TMP_FILE}" "${STATUS_FILE}"
  echo "Watchdog status written: ${STATUS_FILE} (STATE=${STATE})." >&2
}

# ---------------------------------------------------------------------------
# Monitor an array to completion with OOM auto-escalation (see header):
# per-task sacct gate; OUT_OF_MEMORY tasks are re-submitted with doubled
# --mem (ceiling BENCHMARK_MEM_MAX) via ${resubmit_fn}; all-COMPLETED passes;
# any non-COMPLETED, non-OOM state fails closed exactly like the strict gate.
# Optional 5th arg <status_file>: watchdog mode — every terminal path writes
# the status file (and skips notify_sync_status + the JOB_REPORTS append; the
# login tail emails/merges from the file).
# ---------------------------------------------------------------------------
benchmark_wait_oom_retry() {
  local JOB_ID="$1"
  local LABEL="$2"
  local RESUBMIT_FN="$3"
  local MANIFEST="$4"
  local STATUS_FILE="${5:-}"
  local MEM="${BENCHMARK_MEM}"
  local MAX_ATTEMPTS=4
  local ATTEMPT=0
  local TASK_STATES OOM_TASKS=() BAD_TASK="" DS_CSV="" NEW_MEM="" NEW_MANIFEST="" NEW_ID=""
  local JID STATE MASTER_STATE="" TASK_ROWS_FOUND=0 t ds_name CLAMPED=0
  if [[ -n "${STATUS_FILE}" ]]; then
    WATCHDOG_GATED_REPORTS=()
    WATCHDOG_SCHEDULER_IDS=("${JOB_ID}")
    [[ -n "${SLURM_JOB_ID:-}" ]] && WATCHDOG_SCHEDULER_IDS+=("${SLURM_JOB_ID}")
  fi

  while (( ATTEMPT < MAX_ATTEMPTS )); do
    benchmark_wait_array_terminal "${JOB_ID}" "${LABEL}"
    if [[ -n "${STATUS_FILE}" ]]; then
      WATCHDOG_GATED_REPORTS+=("${LABEL}|${JOB_ID}|$(array_wall_time "${JOB_ID}")")
    fi
    # Task-row-only states: `--parsable2` (full, un-truncated values) keeps
    # OUT_OF_MEMORY readable; rows are <job_id>_<n> (batch/extern/master rows
    # filtered out below).
    TASK_STATES="$(sacct -j "${JOB_ID}" -n --parsable2 --format=JobID,State 2>/dev/null || true)"
    if [[ -z "${TASK_STATES//[[:space:]]/}" ]]; then
      echo "ERROR: sacct returned no states for Array Job ${JOB_ID}; NOT syncing to NAS."
      if [[ -n "${STATUS_FILE}" ]]; then
        benchmark_write_status_file "${STATUS_FILE}" FAIL "${LABEL}" \
          "sacct returned no states for array ${JOB_ID} (job purged or unknown id)" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID}: sacct returned no states."
      else
        notify_sync_status \
          "ECODA: benchmark NOT synced (job ${JOB_ID})" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID} (datasets: ${DATASET_NAMES[*]}): sacct returned no states (job purged or unknown id)."
      fi
      exit 1
    fi
    OOM_TASKS=()
    BAD_TASK=""
    MASTER_STATE=""
    TASK_ROWS_FOUND=0
    while IFS='|' read -r JID STATE; do
      JID="$(tr -d '[:space:]' <<< "${JID}")"
      if [[ "${JID}" == "${JOB_ID}" ]]; then
        MASTER_STATE="$(tr -d '[:space:]' <<< "${STATE}")"
        continue
      fi
      [[ "${JID}" =~ ^${JOB_ID}_[0-9]+$ ]] || continue
      TASK_ROWS_FOUND=1
      STATE="$(tr -d '[:space:]' <<< "${STATE}")"
      case "${STATE}" in
        COMPLETED)
          ;;
        OUT_OF_MEMORY)
          OOM_TASKS+=("${JID#${JOB_ID}_}")
          ;;
        *)
          BAD_TASK="${JID} (${STATE})"
          break
          ;;
      esac
    done <<< "${TASK_STATES}"
    # Degenerate case: no task rows at all — rely on the master row so a
    # non-COMPLETED (e.g. CANCELLED) empty array still fails closed.
    if [[ ${TASK_ROWS_FOUND} -eq 0 && "${MASTER_STATE}" != "COMPLETED" ]]; then
      BAD_TASK="${JOB_ID} (master row: ${MASTER_STATE:-unknown})"
    fi
    if [[ -n "${BAD_TASK}" ]]; then
      echo "ERROR: Array Job ${JOB_ID} had non-COMPLETED, non-OOM tasks; NOT syncing to NAS."
      sacct -j "${JOB_ID}" --format=JobID,JobName,State,ExitCode
      if [[ -n "${STATUS_FILE}" ]]; then
        benchmark_write_status_file "${STATUS_FILE}" FAIL "${LABEL}" \
          "non-COMPLETED, non-OOM tasks in array ${JOB_ID}" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID}: non-COMPLETED, non-OOM tasks.
$(benchmark_task_report "${JOB_ID}" "${MANIFEST}")"
      else
        notify_sync_status \
          "ECODA: benchmark NOT synced (job ${JOB_ID})" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID} (datasets: ${DATASET_NAMES[*]}): non-COMPLETED, non-OOM tasks.
$(benchmark_task_report "${JOB_ID}" "${MANIFEST}")"
      fi
      exit 1
    fi
    if [[ ${#OOM_TASKS[@]} -eq 0 ]]; then
      echo "Array Job ${JOB_ID} (${LABEL}): all tasks COMPLETED."
      if [[ -n "${STATUS_FILE}" ]]; then
        benchmark_write_status_file "${STATUS_FILE}" OK "${LABEL}"
      else
        JOB_REPORTS+=("${LABEL}|${JOB_ID}|$(array_wall_time "${JOB_ID}")")
      fi
      return 0
    fi
    # OOM tasks: escalate memory, or fail closed at the ceiling. The
    # MaxRSS report (sacct --format=JobID,State,MaxRSS,Elapsed) documents
    # how far the OOM'd attempt got. The doubled value is CLAMPED to the
    # ceiling so a retry can never request more than the nodes fit (e.g.
    # 512G = 524288 MB would never schedule on the 512000 MB shared-cpu
    # nodes and the squeue poll has no timeout for a PENDING-forever retry);
    # when the current mem is already at/above the ceiling, fail closed.
    if benchmark_mem_ge "${MEM}" "${BENCHMARK_MEM_MAX}"; then
      echo "ERROR: Array Job ${JOB_ID} (${LABEL}) OOM'd at the ${BENCHMARK_MEM_MAX} memory ceiling; NOT syncing to NAS."
      sacct -j "${JOB_ID}" --format=JobID,JobName,State,MaxRSS,Elapsed
      if [[ -n "${STATUS_FILE}" ]]; then
        benchmark_write_status_file "${STATUS_FILE}" FAIL "${LABEL}" \
          "OUT_OF_MEMORY tasks at the ${BENCHMARK_MEM_MAX} ceiling (BENCHMARK_MEM_MAX)" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID}: OUT_OF_MEMORY tasks at the ${BENCHMARK_MEM_MAX} ceiling.
$(benchmark_oom_task_report "${JOB_ID}" "${MANIFEST}")"
      else
        notify_sync_status \
          "ECODA: benchmark NOT synced (OOM at mem ceiling)" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID}: OUT_OF_MEMORY tasks at the ${BENCHMARK_MEM_MAX} ceiling (BENCHMARK_MEM_MAX).
$(benchmark_oom_task_report "${JOB_ID}" "${MANIFEST}")"
      fi
      exit 1
    fi
    if ! NEW_MEM="$(benchmark_bump_mem "${MEM}")"; then
      echo "ERROR: Array Job ${JOB_ID} (${LABEL}) has unparseable mem '${MEM}'; NOT syncing to NAS."
      if [[ -n "${STATUS_FILE}" ]]; then
        benchmark_write_status_file "${STATUS_FILE}" FAIL "${LABEL}" \
          "unparseable BENCHMARK_MEM '${MEM}' (not <N>G/<N>T)" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID}: BENCHMARK_MEM '${MEM}' is not <N>G/<N>T."
      else
        notify_sync_status \
          "ECODA: benchmark NOT synced (unparseable mem)" \
          "Benchmark sync to NAS skipped for ${LABEL} job ${JOB_ID}: BENCHMARK_MEM '${MEM}' is not <N>G/<N>T."
      fi
      exit 1
    fi
    CLAMPED=0
    if benchmark_mem_ge "${NEW_MEM}" "${BENCHMARK_MEM_MAX}"; then
      NEW_MEM="${BENCHMARK_MEM_MAX}"
      CLAMPED=1
    fi
    DS_CSV=""
    for t in "${OOM_TASKS[@]}"; do
      ds_name="$(sed -n "${t}p" "${MANIFEST}")"
      if [[ -z "${ds_name}" ]]; then
        echo "ERROR: No manifest entry for OOM task ${t} in ${MANIFEST}; NOT syncing to NAS." >&2
        if [[ -n "${STATUS_FILE}" ]]; then
          benchmark_write_status_file "${STATUS_FILE}" FAIL "${LABEL}" \
            "no manifest entry for OOM task ${t} in ${MANIFEST}"
        fi
        exit 1
      fi
      [[ -n "${DS_CSV}" ]] && DS_CSV+=","
      DS_CSV+="${ds_name}"
    done
    if [[ ${CLAMPED} -eq 1 ]]; then
      echo "OOM escalation (${LABEL}): task(s) ${OOM_TASKS[*]} -> dataset(s) ${DS_CSV}; retrying with mem ${MEM} -> ${NEW_MEM} (clamped to BENCHMARK_MEM_MAX=${BENCHMARK_MEM_MAX}) (attempt $((ATTEMPT + 1)) of ${MAX_ATTEMPTS})."
    else
      echo "OOM escalation (${LABEL}): task(s) ${OOM_TASKS[*]} -> dataset(s) ${DS_CSV}; retrying with mem ${MEM} -> ${NEW_MEM} (attempt $((ATTEMPT + 1)) of ${MAX_ATTEMPTS})."
    fi
    if [[ -n "${ANALYSIS_PASS:-}" ]]; then
      NEW_MANIFEST="${ANALYSIS_ROOT}/manifests/batch_effect_${ANALYSIS_PASS}_manifest_${LABEL}_retry_$$.txt"
    else
      NEW_MANIFEST="${HPC_SCRATCH_DIR}/benchmark_manifest_${LABEL}_retry_$$.txt"
    fi
    NEW_ID="$("${RESUBMIT_FN}" "${LABEL}" "${DS_CSV}" "${NEW_MEM}" "${NEW_MANIFEST}")"
    if [[ -z "${NEW_ID}" ]]; then
      echo "ERROR: ${RESUBMIT_FN} returned no array job id for the ${LABEL} retry; NOT syncing to NAS." >&2
      if [[ -n "${STATUS_FILE}" ]]; then
        benchmark_write_status_file "${STATUS_FILE}" FAIL "${LABEL}" \
          "resubmit function returned no array job id for the ${LABEL} retry"
      fi
      exit 1
    fi
    if [[ -n "${STATUS_FILE}" ]]; then
      WATCHDOG_SCHEDULER_IDS+=("${NEW_ID}")
    fi
    JOB_ID="${NEW_ID}"
    MANIFEST="${NEW_MANIFEST}"
    MEM="${NEW_MEM}"
    ATTEMPT=$((ATTEMPT + 1))
  done
  echo "ERROR: ${LABEL} exceeded ${MAX_ATTEMPTS} OOM retry attempts; NOT syncing to NAS." >&2
  if [[ -n "${STATUS_FILE}" ]]; then
    benchmark_write_status_file "${STATUS_FILE}" FAIL "${LABEL}" \
      "exceeded ${MAX_ATTEMPTS} OOM retry attempts"
  fi
  exit 1
}

# ---------------------------------------------------------------------------
# Monitor an array to completion, then fail-closed sacct gate (all rows)
# ---------------------------------------------------------------------------
benchmark_wait_for_array() {
  local JOB_ID="$1"
  local LABEL="$2"
  benchmark_wait_array_terminal "${JOB_ID}" "${LABEL}"
  local STATES
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
$(benchmark_task_report "${JOB_ID}")"
    exit 1
  fi
  echo "Array Job ${JOB_ID} (${LABEL}): all tasks COMPLETED."
  JOB_REPORTS+=("${LABEL}|${JOB_ID}|$(array_wall_time "${JOB_ID}")")
}

# ---------------------------------------------------------------------------
# Submit one compute-node watchdog job (watchdog_main.sh) that owns the
# terminal wait + OOM escalation for a method array (see header). The
# watchdog runs on the method's partition WITHOUT the pinned constraint class
# (modest 1 cpu/2G job; it must never compete with the pinned workers or
# distort the benchmark resource picture) and logs to
# ${LOGS_DIR}/5_benchmark_watchdog_<label>_<id>.log/.err. The retry arrays it
# submits reuse the given throttle/log_prefix/worker_script/flags (the flags
# carry the per-method --gpus/--constraint/--cpus-per-task pins). Echoes ONLY
# the watchdog job id on stdout (progress to stderr) so the caller can capture
# it with $(...) — a multi-line capture would break the gates downstream.
# ---------------------------------------------------------------------------
benchmark_submit_watchdog() {
  local ARRAY_ID="$1"
  local LABEL="$2"
  local MANIFEST="$3"
  local MODE="$4"
  local PARTITION="$5"
  local THROTTLE="$6"
  local LOG_PREFIX="$7"
  local WORKER_SCRIPT="$8"
  shift 8
  local WATCHDOG_FLAGS=("$@")

  mkdir -p "${WATCHDOG_STATUS_DIR}"
  echo "Submitting ${LABEL} watchdog for array ${ARRAY_ID} (mode=${MODE}, partition=${PARTITION}, " >&2
  echo "  time=${WATCHDOG_TIME_LIMIT}, flags for retries: ${WATCHDOG_FLAGS[*]})" >&2

  local WD_JOB_NAME WD_OUTPUT_PREFIX
  if [[ -n "${ANALYSIS_PASS:-}" ]]; then
    WD_JOB_NAME="batch_effect_watchdog_${ANALYSIS_PASS}_${LABEL}"
    WD_OUTPUT_PREFIX="5_batch_effect_watchdog_${ANALYSIS_PASS}_${LABEL}"
  else
    WD_JOB_NAME="benchmark_watchdog_${LABEL}"
    WD_OUTPUT_PREFIX="5_benchmark_watchdog_${LABEL}"
  fi
  local SUBMIT_MSG
  SUBMIT_MSG=$(sbatch \
      --job-name="${WD_JOB_NAME}" \
      --ntasks=1 --cpus-per-task=1 --mem=2G \
      --time="${WATCHDOG_TIME_LIMIT}" \
      --partition="${PARTITION}" \
      --output="${LOGS_DIR}/${WD_OUTPUT_PREFIX}_%A.log" \
      --error="${LOGS_DIR}/${WD_OUTPUT_PREFIX}_%A.err" \
      --mail-user="${USER_EMAIL}" \
      "${WATCHDOG_MAIN_SCRIPT}" \
      "${ARRAY_ID}" "${LABEL}" "${MANIFEST}" "${MODE}" -- \
      "${PARTITION}" "${THROTTLE}" "${LOG_PREFIX}" "${WORKER_SCRIPT}" "${WATCHDOG_FLAGS[@]}")

  local WATCHDOG_ID
  WATCHDOG_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "  ${LABEL} watchdog job ID: ${WATCHDOG_ID} (status file: ${WATCHDOG_STATUS_DIR}/${WATCHDOG_ID}.status)" >&2
  echo "${WATCHDOG_ID}"
}

# ---------------------------------------------------------------------------
# Login-tail counterpart of the watchdog (see header): wait for the watchdog
# job to leave the scheduler, poll for its status file (<=2 min grace), then
# branch on its STATE:
#   OK    -> merge the JOB_REPORT= lines into JOB_REPORTS, return 0.
#   FAIL  -> print the report, notify_sync_status ("NOT synced — watchdog
#            failed"), exit 1.
#   watchdog non-COMPLETED without a status file (FAILED/TIMEOUT/PREEMPTED/
#   CANCELLED, e.g. node panic before the file was written) -> fail closed
#   with its sacct State,ExitCode + a pointer to its logs.
#   COMPLETED but no status file after the grace -> fail closed ("watchdog
#   exited without a status file").
# ---------------------------------------------------------------------------
benchmark_wait_watchdog() {
  local WATCHDOG_ID="$1"
  local LABEL="$2"
  local STATUS_FILE="${WATCHDOG_STATUS_DIR}/${WATCHDOG_ID}.status"
  local WSTATE REASON REPORT WD_STATE WD_EXIT line GRACE

  benchmark_wait_array_terminal "${WATCHDOG_ID}" "${LABEL} watchdog"

  # Grace for the status file: the watchdog writes it right before exiting;
  # allow up to 2 min for scheduler/sacct lag and file visibility
  # (WATCHDOG_STATUS_GRACE_ITERS overridable for tests).
  GRACE=0
  while [[ ! -s "${STATUS_FILE}" && ${GRACE} -lt ${WATCHDOG_STATUS_GRACE_ITERS:-24} ]]; do
    sleep 5
    GRACE=$((GRACE + 1))
  done

  if [[ -s "${STATUS_FILE}" ]]; then
    WSTATE="$(grep -E '^STATE=' "${STATUS_FILE}" | head -1 | cut -d= -f2- | tr -d '[:space:]')"
    case "${WSTATE}" in
      OK)
        echo "${LABEL} watchdog ${WATCHDOG_ID}: STATE=OK."
        while IFS= read -r line; do
          JOB_REPORTS+=("${line#JOB_REPORT=}")
        done < <(grep '^JOB_REPORT=' "${STATUS_FILE}" || true)
        return 0
        ;;
      FAIL)
        REASON="$(grep -E '^FAIL_REASON=' "${STATUS_FILE}" | head -1 | cut -d= -f2-)"
        REPORT="$(sed -n '/^REPORT=$/,$p' "${STATUS_FILE}" | tail -n +2)"
        echo "ERROR: ${LABEL} watchdog ${WATCHDOG_ID}: STATE=FAIL (${REASON:-unknown}); NOT syncing to NAS." >&2
        notify_sync_status \
          "ECODA: benchmark NOT synced (watchdog failed)" \
          "Benchmark sync to NAS skipped for ${LABEL}: watchdog job ${WATCHDOG_ID} failed (${REASON:-unknown}).
${REPORT}"
        exit 1
        ;;
      *)
        echo "ERROR: ${LABEL} watchdog ${WATCHDOG_ID}: unparseable status file (STATE='${WSTATE:-}'); NOT syncing to NAS." >&2
        notify_sync_status \
          "ECODA: benchmark NOT synced (unparseable watchdog status)" \
          "Benchmark sync to NAS skipped for ${LABEL}: watchdog job ${WATCHDOG_ID} wrote an unparseable status file (${STATUS_FILE})."
        exit 1
        ;;
    esac
  fi

  # No status file after the grace: fail closed on the watchdog job's own
  # state (a non-COMPLETED watchdog — FAILED/TIMEOUT/PREEMPTED/CANCELLED, e.g.
  # node panic before the file was written — points at its logs; a COMPLETED
  # watchdog that wrote nothing is a bug in the status protocol).
  WD_STATE="$(sacct -j "${WATCHDOG_ID}" -X -n --format=State 2>/dev/null | head -1 | tr -d '[:space:]' || true)"
  WD_EXIT="$(sacct -j "${WATCHDOG_ID}" -X -n --format=ExitCode 2>/dev/null | head -1 | tr -d '[:space:]' || true)"
  echo "ERROR: ${LABEL} watchdog ${WATCHDOG_ID} exited without a status file (sacct State=${WD_STATE:-n/a}, ExitCode=${WD_EXIT:-n/a}); NOT syncing to NAS." >&2
  echo "  Check ${LOGS_DIR}/5_benchmark_watchdog_${LABEL}_${WATCHDOG_ID}.log/.err" >&2
  notify_sync_status \
    "ECODA: benchmark NOT synced (watchdog lost)" \
    "Benchmark sync to NAS skipped for ${LABEL}: watchdog job ${WATCHDOG_ID} exited without a status file (sacct State=${WD_STATE:-n/a}, ExitCode=${WD_EXIT:-n/a}).
Check ${LOGS_DIR}/5_benchmark_watchdog_${LABEL}_${WATCHDOG_ID}.log/.err; recover with --sync-only ${WATCHDOG_ID} or a re-run (idempotent)."
  exit 1
}

# ---------------------------------------------------------------------------
# NAS check -> RDS integrity sidecar -> merge exec logs -> rsync -> cleanup
# ---------------------------------------------------------------------------

benchmark_sync_artifacts_for() {
  local ds="$1" label="$2" n suffix stem
  local root="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
  SYNC_ARTIFACTS=()
  if [[ "${label}" == prepare_pseudobulk ]]; then
    if [[ -n "${ANALYSIS_PASS:-}" ]]; then
      SYNC_ARTIFACTS+=("${root}/pseudobulks/${ds}_batch_effect_${ANALYSIS_PASS}_pseudobulk_hvg2000.rds")
    else
      for suffix in schvg2000 hvg2000 hvg500 hvg2000_bl hvg1000 hvg3000; do
        SYNC_ARTIFACTS+=("${root}/pseudobulks/${ds}_pseudobulk_${suffix}.rds")
      done
    fi
    return 0
  fi
  case "${label}" in
    mrvi)
      if [[ -n "${ANALYSIS_PASS:-}" ]]; then
        SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_batch_effect_${ANALYSIS_PASS}_hvg2000_highres_mrvi_dists.feather")
      else
        for n in 1000 2000 3000; do
          SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_mrvi_dists.feather")
        done
      fi
      ;;
    scpoli)
      SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_lowres_scpoli_dims15_embs.feather")
      for n in 1000 3000; do
        SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_highres_scpoli_dims15_embs.feather")
      done
      for n in 2 3 5 10 15; do
        SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_highres_scpoli_dims${n}_embs.feather")
      done
      ;;
    pilot|qot|pilotgm)
      suffix="${label}"
      if [[ -n "${ANALYSIS_PASS:-}" ]]; then
        SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_batch_effect_${ANALYSIS_PASS}_hvg2000_highres_${suffix}_dists.feather")
      else
        SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_lowres_${suffix}_dists.feather")
        for n in 1000 2000 3000; do
          SYNC_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_highres_${suffix}_dists.feather")
        done
      fi
      ;;
    trans|zeroimp)
      SYNC_ARTIFACTS+=("${root}/results/${ds}_${label}.rds")
      ;;
    gloscope|mofa|pseudobulk|scitd)
      stem="${ds}"
      [[ -n "${ANALYSIS_PASS:-}" ]] && stem="${ds}_batch_effect_${ANALYSIS_PASS}"
      SYNC_ARTIFACTS+=("${root}/results/${stem}_${label}.rds")
      ;;
    composition)
      stem="${ds}"
      [[ -n "${ANALYSIS_PASS:-}" ]] && stem="${ds}_batch_effect_${ANALYSIS_PASS}"
      SYNC_ARTIFACTS+=("${root}/results/${stem}_composition.rds")
      SYNC_ARTIFACTS+=("${root}/results/${stem}_metadata.rds")
      ;;
    *)
      echo "ERROR: unsupported benchmark sync label '${label}'." >&2
      return 1
      ;;
  esac
}
analysis_merge_sync_cleanup() (
  local LABELS=("$@")
  local LOCAL_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
  local REMOTE_ROOT="${ANALYSIS_NAS_ROOT:-${NAS_TARGET_DIR}/benchmark}"
  local LOG_PREFIX="${ANALYSIS_LOG_PREFIX:-execution_times_}"
  local RUN_ROOT="${ECODA_RUN_ROOT:-}"
  local RUN_ID="${ECODA_RUN_ID:-}"
  local RUN_LOG_DIR="${EXECUTION_LOG_DIR:-${RUN_ROOT}/logs}"
  local KIND="benchmark" KIND_CAP="Benchmark"
  local SYNC_OWNER_DIR="" SYNC_LOCK_DIR="" SYNC_FINAL_STATE="FAIL"
  local SYNC_FILES SYNC_FILES_TMP NO_CHECKSUM_FILES_TMP CLEANUP_MANIFEST
  local CHECKSUM_TMP REMOTE_CHECKSUM_TMP EXISTING_LOG
  local ds view row_label label path rel line selected_seen=""
  local -a SELECTED_RELS=()
  sync_fail() {
    local reason="$1"
    SYNC_FINAL_STATE="FAIL"
    echo "ERROR: ${reason}" >&2
    notify_sync_status \
      "ECODA: ${KIND} NOT synced" \
      "${KIND_CAP} sync to NAS skipped (datasets: ${DATASET_NAMES[*]}, labels: ${LABELS[*]}): ${reason}" || true
    exit 1
  }
  sync_cleanup() {
    local rc="$?"
    if [[ -n "${SYNC_OWNER_DIR}" ]]; then
      if ! ecoda_owner_set_state "${SYNC_OWNER_DIR}" "${SYNC_FINAL_STATE}" \
          "Stage 5 synchronization terminal state"; then
        rc=1
      fi
    fi
    if [[ -n "${SYNC_LOCK_DIR}" ]]; then
      rmdir "${SYNC_LOCK_DIR}" 2>/dev/null || rc=1
    fi
    exit "${rc}"
  }
  trap sync_cleanup EXIT

  [[ -n "${RUN_ROOT}" && -n "${RUN_ID}" ]] || sync_fail "Stage 5 run root/id is not set"
  if [[ "${LOCAL_ROOT}" != "${HPC_SCRATCH_DIR}/benchmark" ]]; then
    KIND="analysis"
    KIND_CAP="Analysis"
  fi
  echo "Checking NAS reachability before any destructive merge work..."
  [[ -d "${NAS_TARGET_DIR}" ]] || sync_fail "NAS path ${REMOTE_ROOT} is unreachable"
  if ! mkdir -p "${REMOTE_ROOT}" "${LOCAL_ROOT}/embeddings" "${RUN_ROOT}/manifests"; then
    sync_fail "failed to create benchmark sync directories"
  fi
  LOCAL_ROOT="$(cd "${LOCAL_ROOT}" && pwd)" ||
    sync_fail "failed to canonicalize local benchmark root"
  REMOTE_ROOT="$(cd "${REMOTE_ROOT}" && pwd)" ||
    sync_fail "failed to canonicalize remote benchmark root"
  RUN_ROOT="$(cd "${RUN_ROOT}" && pwd)" ||
    sync_fail "failed to canonicalize Stage 5 run root"
  RUN_LOG_DIR="$(cd "${RUN_LOG_DIR}" && pwd)" ||
    sync_fail "failed to canonicalize Stage 5 run log directory"
  ECODA_RUN_ROOT="${RUN_ROOT}"
  SYNC_OWNER_DIR="$(ecoda_owner_acquire stage5 "sync/${LOCAL_ROOT}" "${RUN_ID}" 0 0)" || {
    sync_fail "shared Stage 5 sync owner is unavailable"
  }
  SYNC_LOCK_DIR="${RUN_ROOT}/sync.lock"
  mkdir "${SYNC_LOCK_DIR}" 2>/dev/null || sync_fail "run sync lock already exists: ${SYNC_LOCK_DIR}"

  SYNC_FILES="${RUN_ROOT}/manifests/sync_files.tsv"
  SYNC_FILES_TMP="${SYNC_FILES}.build.$$"
  NO_CHECKSUM_FILES_TMP="${RUN_ROOT}/manifests/sync_files.no_checksum.build.$$"
  CLEANUP_MANIFEST="${RUN_ROOT}/manifests/cleanup.tsv"
  if ! : > "${SYNC_FILES_TMP}"; then
    sync_fail "failed to create selected sync manifest"
  fi
  SELECTED_RELS=()

  add_sync_artifact() {
    local artifact="$1" artifact_rel
    [[ "${artifact}" == "${LOCAL_ROOT}/"* ]] || sync_fail "selected artifact is outside analysis root: ${artifact}"
    artifact_rel="${artifact#${LOCAL_ROOT}/}"
    case " ${selected_seen} " in
      *" ${artifact_rel} "*) return 0 ;;
    esac
    [[ -s "${artifact}" ]] || sync_fail "selected artifact is missing or empty: ${artifact}"
    ecoda_validate_checksum "${artifact}" || sync_fail "selected artifact checksum failed: ${artifact}"
    selected_seen="${selected_seen} ${artifact_rel}"
    SELECTED_RELS+=("${artifact_rel}")
    printf '%s\n' "${artifact_rel}" >> "${SYNC_FILES_TMP}"
    printf '%s\n' "${artifact_rel}.md5" >> "${SYNC_FILES_TMP}"
  }

  if [[ -r "${ECODA_SELECTION_MANIFEST:-}" ]]; then
    while IFS=$'\t' read -r ds view row_label; do
      [[ -n "${ds}" && -n "${view}" ]] || continue
      SYNC_LABELS=("${LABELS[@]}")
      if [[ "${ECODA_EXACT_SELECTION:-0}" == "1" && -z "${ANALYSIS_PASS:-}" ]]; then
        SYNC_LABELS=("${row_label}")
      fi
      for label in "${SYNC_LABELS[@]}"; do
        benchmark_sync_artifacts_for "${ds}" "${label}" || sync_fail "cannot resolve ${ds}/${view}/${label}"
        for path in "${SYNC_ARTIFACTS[@]}"; do add_sync_artifact "${path}"; done
      done
    done < "${ECODA_SELECTION_MANIFEST}"
  else
    for ds in "${DATASET_NAMES[@]}"; do
      for label in "${LABELS[@]}"; do
        benchmark_sync_artifacts_for "${ds}" "${label}" || sync_fail "cannot resolve ${ds}/${label}"
        for path in "${SYNC_ARTIFACTS[@]}"; do add_sync_artifact "${path}"; done
      done
    done
  fi
  [[ ${#SELECTED_RELS[@]} -gt 0 ]] || sync_fail "no selected benchmark artifacts"
  if [[ -s "${REMOTE_ROOT}/checksums.md5" ]]; then
    (cd "${REMOTE_ROOT}" && md5sum -c checksums.md5 >/dev/null 2>&1) || {
      sync_fail "existing remote checksum manifest failed validation"
    }
  fi

  EXISTING_LOG="${REMOTE_ROOT}/embeddings/execution_times.feather"
  if [[ -f "${EXISTING_LOG}" && ! -f "${EXISTING_LOG}.md5" ]]; then
    "${PYTHON_BIN}" "${ANALYSIS_MERGE_SCRIPT}" \
      --migrate-existing-log "${EXISTING_LOG}" || {
      sync_fail "existing remote execution log failed guarded sidecar migration"
    }
  fi
  if [[ -e "${EXISTING_LOG}" || -e "${EXISTING_LOG}.md5" ]]; then
    ecoda_validate_checksum_remote "${EXISTING_LOG}" "${EXISTING_LOG}.md5" || {
      sync_fail "existing remote execution log checksum failed"
    }
  fi
  MERGE_ARGS=(--output_dir "${LOCAL_ROOT}/embeddings"
              --log-dir "${RUN_LOG_DIR}"
              --no-cleanup
              --filename_prefix "${LOG_PREFIX}"
              --labels "${LABELS[@]}"
              --datasets "${DATASET_NAMES[@]}"
              --cleanup-manifest "${CLEANUP_MANIFEST}")
  [[ -s "${EXISTING_LOG}" ]] && MERGE_ARGS+=(--existing-log "${EXISTING_LOG}")
  "${PYTHON_BIN}" "${ANALYSIS_MERGE_SCRIPT}" "${MERGE_ARGS[@]}" || {
    sync_fail "execution-log merge failed"
  }
  add_sync_artifact "${LOCAL_ROOT}/embeddings/execution_times.feather"

  CHECKSUM_TMP="${LOCAL_ROOT}/checksums.md5.build.$$"
  if ! : > "${CHECKSUM_TMP}"; then
    sync_fail "failed to create benchmark checksum manifest"
  fi
  if [[ -s "${REMOTE_ROOT}/checksums.md5" ]]; then
    while IFS= read -r line || [[ -n "${line}" ]]; do
      rel="$(printf '%s\n' "${line}" | awk '{p=$2; sub(/^\*/, "", p); print p}')"
      case " ${selected_seen} " in
        *" ${rel} "*) ;;
        *) printf '%s\n' "${line}" >> "${CHECKSUM_TMP}" ||
             sync_fail "failed to preserve remote checksum entry" ;;
      esac
    done < "${REMOTE_ROOT}/checksums.md5"
  fi
  for rel in "${SELECTED_RELS[@]}"; do
    (cd "${LOCAL_ROOT}" && md5sum "${rel}") >> "${CHECKSUM_TMP}" || {
      rm -f "${CHECKSUM_TMP}"
      sync_fail "cannot checksum selected relative path: ${rel}"
    }
  done
  [[ -s "${CHECKSUM_TMP}" ]] || {
    rm -f "${CHECKSUM_TMP}"
    sync_fail "final benchmark checksum manifest is empty"
  }
  if ! mv -f "${CHECKSUM_TMP}" "${LOCAL_ROOT}/checksums.md5"; then
    rm -f "${CHECKSUM_TMP}"
    sync_fail "cannot install local benchmark checksum manifest"
  fi
  printf 'checksums.md5\n' >> "${SYNC_FILES_TMP}" ||
    sync_fail "cannot add checksum manifest to selected sync"
  [[ -s "${SYNC_FILES_TMP}" ]] || sync_fail "selected sync file manifest is empty"
  if ! mv -f "${SYNC_FILES_TMP}" "${SYNC_FILES}"; then
    rm -f "${SYNC_FILES_TMP}"
    sync_fail "cannot install selected sync manifest"
  fi

  if ! : > "${NO_CHECKSUM_FILES_TMP}"; then
    sync_fail "failed to create no-checksum sync manifest"
  fi
  while IFS= read -r rel || [[ -n "${rel}" ]]; do
    if [[ "${rel}" != "checksums.md5" ]]; then
      printf '%s\n' "${rel}" >> "${NO_CHECKSUM_FILES_TMP}" ||
        sync_fail "failed to build no-checksum sync manifest"
    fi
  done < "${SYNC_FILES}"
  rsync -rlptDv --files-from="${NO_CHECKSUM_FILES_TMP}" \
    "${LOCAL_ROOT}/" "${REMOTE_ROOT}/" || sync_fail "selected benchmark rsync failed"
  REMOTE_CHECKSUM_TMP="${REMOTE_ROOT}/.checksums.md5.tmp.${RUN_ID}"
  cp "${LOCAL_ROOT}/checksums.md5" "${REMOTE_CHECKSUM_TMP}" || {
    sync_fail "cannot stage remote checksum manifest"
  }
  mv -f "${REMOTE_CHECKSUM_TMP}" "${REMOTE_ROOT}/checksums.md5" || {
    sync_fail "cannot atomically install remote checksum manifest"
  }
  (cd "${REMOTE_ROOT}" && md5sum -c checksums.md5) || sync_fail "remote checksum verification failed"
  if [[ -s "${CLEANUP_MANIFEST}" ]]; then
    while IFS= read -r path || [[ -n "${path}" ]]; do
      [[ "${path}" == "${RUN_LOG_DIR}/"* ]] ||
        sync_fail "cleanup manifest escapes run log directory: ${path}"
      rm -f "${path}" "${path}.md5" ||
        sync_fail "failed to clean up run log artifact: ${path}"
    done < "${CLEANUP_MANIFEST}"
  fi
  notify_sync_status \
    "ECODA: ${KIND} synced to NAS" \
    "${KIND_CAP} selected results synced to ${REMOTE_ROOT}/ (datasets: ${DATASET_NAMES[*]}, labels: ${LABELS[*]}).
$(benchmark_job_durations_block)" || true
  SYNC_FINAL_STATE="OK"
)

# Existing submitters retain their public helper name; batch-effect submitters
# call the generic implementation directly.
benchmark_merge_sync_cleanup() {
  analysis_merge_sync_cleanup "$@"
}
