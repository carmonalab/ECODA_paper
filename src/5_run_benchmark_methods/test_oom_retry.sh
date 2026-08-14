#!/usr/bin/env bash
#
# Mocked unit tests for the benchmark OOM auto-escalation + compute-node
# watchdog machinery:
#   - src/5_run_benchmark_methods/benchmark_submit_common.sh (sourced):
#     benchmark_bump_mem / benchmark_mem_ge / benchmark_pb_variant_names /
#     benchmark_wait_oom_retry (status-file mode) /
#     benchmark_submit_watchdog / benchmark_wait_watchdog
#   - src/5_run_benchmark_methods/watchdog_main.sh (functions extracted via
#     awk, then eval'd): watchdog_main (strict + soft-gate) / watchdog_resubmit
#
# The slurm CLI (squeue/sacct/sbatch) and notify_sync_status are stubbed:
# squeue says "not in scheduler" (terminal immediately), sacct serves canned
# rows from the per-test globals SACCT_ROWS/SACCT_MEM_ROWS/SACCT_XEXIT/
# SACCT_WALL, sbatch echoes "Submitted batch job <id>" and captures its args
# to ${CAPTURE_DIR}/sbatch_calls.txt, notify_sync_status appends to
# ${CAPTURE_DIR}/notify_calls.txt. Captures live on DISK because the
# functions under test run in subshells (their exit paths call `exit`), so
# plain shell globals would be lost. The status-file grace poll is
# short-circuited via WATCHDOG_STATUS_GRACE_ITERS=0 and sleep is stubbed.
#
# Run: bash src/5_run_benchmark_methods/test_oom_retry.sh
# (bash >= 3.2; no bashisms newer than 3.2 — associative arrays / mapfile /
# negative array indices are NOT used.)

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export USER_EMAIL="test@example.com"

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda_watchdog_test.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT

# shellcheck source=/dev/null
source "${REPO_ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
# Override the HPC paths BEFORE sourcing the common file (WATCHDOG_STATUS_DIR
# is derived from HPC_SCRATCH_DIR at source time; slurm_config.sh must come
# first as it assigns HPC_SCRATCH_DIR/LOGS_DIR unconditionally).
export HPC_SCRATCH_DIR="${TMP_DIR}/scratch"
export LOGS_DIR="${TMP_DIR}/logs"
export DATASETS_JSON_FILE="${REPO_ROOT}/datasets.json"
export BENCHMARK_MEM="${BENCHMARK_MEM:-128G}"
export BENCHMARK_MEM_MAX="${BENCHMARK_MEM_MAX:-500G}"
export FORCE_BENCHMARK=0
export WATCHDOG_STATUS_GRACE_ITERS=0

# shellcheck source=/dev/null
source "${REPO_ROOT}/src/5_run_benchmark_methods/benchmark_submit_common.sh"

# watchdog_main.sh executes `watchdog_main "$@"` at the bottom, so its
# functions are text-extracted (awk stops at the first bare-} line) and
# eval'd instead of sourcing the whole script.
extract_fn() {
  local FILE="$1" FN="$2"
  awk -v fn="${FN}" '
    index($0, fn "() {") == 1 { capture = 1 }
    capture { print }
    capture && /^}$/ { exit }
  ' "${FILE}"
}
WATCHDOG_MAIN_SH="${REPO_ROOT}/src/5_run_benchmark_methods/watchdog_main.sh"
SCRIPT_DIR="$(dirname "${WATCHDOG_MAIN_SH}")"
eval "$(extract_fn "${WATCHDOG_MAIN_SH}" watchdog_resubmit)"
eval "$(extract_fn "${WATCHDOG_MAIN_SH}" watchdog_main)"

# ---------------------------------------------------------------------------
# Stubs (captures on disk; see header)
# ---------------------------------------------------------------------------
CAPTURE_DIR="${TMP_DIR}/captures"
reset_captures() {
  rm -rf "${CAPTURE_DIR}"
  mkdir -p "${CAPTURE_DIR}"
  echo "22222" > "${CAPTURE_DIR}/sbatch_next.txt"
  : > "${CAPTURE_DIR}/sbatch_calls.txt"
  : > "${CAPTURE_DIR}/sbatch_env.txt"
  : > "${CAPTURE_DIR}/notify_calls.txt"
  echo "22222" > "${CAPTURE_DIR}/resubmit_retry_id.txt"
  : > "${CAPTURE_DIR}/resubmit_capture.txt"
}

squeue() { return 1; }   # nothing in the scheduler -> terminal immediately
sleep() { :; }           # no real waits (poll loops otherwise take minutes)

SACCT_ROWS=""            # "<jid>:<STATE>" entries (master + task rows)
SACCT_MEM_ROWS=""        # "<jid>:<MaxRSS>" entries (optional)
SACCT_XEXIT="0:0"
SACCT_WALL="00:05:00"
sacct() {
  local jid="" fmt="" XS=0 entries=() e m mem
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -j) jid="$2"; shift 2 ;;
      -X) XS=1; shift ;;
      --format=*) fmt="${1#--format=}"; shift ;;
      *) shift ;;
    esac
  done
  # sacct -j <id> returns the master row AND the <id>_<n> task rows.
  for e in ${SACCT_ROWS}; do
    case "${e%%:*}" in
      ${jid}|${jid}_*) entries+=("${e}") ;;
    esac
  done
  [[ ${#entries[@]} -eq 0 ]] && return 0
  case "${fmt}" in
    JobID,State)
      for e in "${entries[@]}"; do printf '%s|%s\n' "${e%%:*}" "${e#*:}"; done
      ;;
    State)
      if [[ ${XS} -eq 1 ]]; then
        printf '%s\n' "${entries[0]#*:}"
      else
        for e in "${entries[@]}"; do printf '%s\n' "${e#*:}"; done
      fi
      ;;
    JobID,State,MaxRSS,Elapsed)
      for e in "${entries[@]}"; do
        mem="n/a"
        for m in ${SACCT_MEM_ROWS}; do
          [[ "${m%%:*}" == "${e%%:*}" ]] && mem="${m#*:}"
        done
        printf '%s|%s|%s|00:04:00\n' "${e%%:*}" "${e#*:}" "${mem}"
      done
      ;;
    JobID,State,Elapsed,ExitCode)
      for e in "${entries[@]}"; do printf '%s|%s|00:04:00|0:0\n' "${e%%:*}" "${e#*:}"; done
      ;;
    ExitCode)
      printf '%s\n' "${SACCT_XEXIT}"
      ;;
    Elapsed)
      printf '%s\n' "${SACCT_WALL}"
      ;;
    *)
      for e in "${entries[@]}"; do printf '%s|%s|0:0\n' "${e%%:*}" "${e#*:}"; done
      ;;
  esac
}

sbatch() {
  echo "$*" >> "${CAPTURE_DIR}/sbatch_calls.txt"
  # The sbatch stub runs in the caller's process, so it can observe the env
  # exports the closures make for the workers (METHOD/BENCHMARK_MANIFEST).
  printf 'METHOD=%s\nMANIFEST=%s\n' "${METHOD:-}" "${BENCHMARK_MANIFEST:-}" \
    >> "${CAPTURE_DIR}/sbatch_env.txt"
  local next
  next="$(cat "${CAPTURE_DIR}/sbatch_next.txt")"
  echo "Submitted batch job ${next}"
  echo $((next + 1)) > "${CAPTURE_DIR}/sbatch_next.txt"
}

notify_sync_status() {
  # Flatten newlines so the capture stays one BODY= line per call.
  printf 'SUBJECT=%s\nBODY=%s\n---\n' "$1" "$(printf '%s' "$2" | tr '\n' ' ')" \
    >> "${CAPTURE_DIR}/notify_calls.txt"
}
notify_count() {
  local c
  c="$(grep -c '^SUBJECT=' "${CAPTURE_DIR}/notify_calls.txt" 2>/dev/null)"
  echo "${c:-0}"
}
notify_last_subject() {
  grep '^SUBJECT=' "${CAPTURE_DIR}/notify_calls.txt" | tail -1 | cut -d= -f2-
}
notify_last_body() {
  grep '^BODY=' "${CAPTURE_DIR}/notify_calls.txt" | tail -1 | cut -d= -f2-
}

# Resubmit closure stub: records its args to disk, writes the retry manifest,
# echoes the configured retry array id ("" -> the "no id" fail-closed path).
# resubmit_retry_id.txt holds one id per line; the stub POPS the first so a
# test can script a chain of retries (e.g. for the attempt-cap path).
resubmit_stub() {
  printf 'LABEL=%s\nDS=%s\nMEM=%s\nMANIFEST=%s\n' "$1" "$2" "$3" "$4" \
    > "${CAPTURE_DIR}/resubmit_capture.txt"
  : > "$4"
  local IFS=',' ds
  for ds in $2; do
    echo "${ds}" >> "$4"
  done
  local ids rest
  ids="$(head -1 "${CAPTURE_DIR}/resubmit_retry_id.txt")"
  rest="$(tail -n +2 "${CAPTURE_DIR}/resubmit_retry_id.txt")"
  echo "${ids}" > "${CAPTURE_DIR}/resubmit_retry_id.txt"
  printf '%s\n' "${rest}" >> "${CAPTURE_DIR}/resubmit_retry_id.txt"
  echo "${ids}"
}

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------
TESTS=0
FAILURES=0
assert_eq() {
  local desc="$1" expected="$2" got="$3"
  TESTS=$((TESTS + 1))
  if [[ "${expected}" == "${got}" ]]; then
    echo "  ok: ${desc}"
  else
    echo "  FAIL: ${desc} (expected '${expected}', got '${got}')"
    FAILURES=$((FAILURES + 1))
  fi
}
assert_contains() {
  local desc="$1" needle="$2" haystack="$3"
  TESTS=$((TESTS + 1))
  if [[ "${haystack}" == *"${needle}"* ]]; then
    echo "  ok: ${desc}"
  else
    echo "  FAIL: ${desc} ('${haystack}' does not contain '${needle}')"
    FAILURES=$((FAILURES + 1))
  fi
}
assert_exit() {
  local desc="$1" expected="$2"
  shift 2
  TESTS=$((TESTS + 1))
  ( "$@" ) >/dev/null 2>&1
  local rc=$?
  if [[ ${rc} -eq ${expected} ]]; then
    echo "  ok: ${desc}"
  else
    echo "  FAIL: ${desc} (expected exit ${expected}, got ${rc})"
    FAILURES=$((FAILURES + 1))
  fi
}

setup_oom_scenario() {
  mkdir -p "${HPC_SCRATCH_DIR}"
  reset_captures
  SACCT_ROWS=""
  SACCT_MEM_ROWS=""
  SACCT_XEXIT="0:0"
  SACCT_WALL="00:05:00"
  JOB_REPORTS=()
  WATCHDOG_GATED_REPORTS=()
  export BENCHMARK_MEM="128G"
  export BENCHMARK_MEM_MAX="500G"
  export FORCE_BENCHMARK=0
  export SLURM_JOB_ID=""
}

# Runs watchdog_main in a subshell (it calls `exit`, which would kill the
# harness), recording its exit code via an EXIT trap (traps defined inside a
# subshell fire when the subshell exits).
run_watchdog() {
  (
    trap 'echo "RC=$?" > "${CAPTURE_DIR}/wd_rc.txt"' EXIT
    watchdog_main "$@"
  ) >/dev/null 2>&1
}
wd_rc() {
  cut -d= -f2 "${CAPTURE_DIR}/wd_rc.txt" 2>/dev/null
}

# Shared 3-dataset manifest for the OOM mapping assertions.
MANIFEST_FILE=""
make_manifest() {
  MANIFEST_FILE="${HPC_SCRATCH_DIR}/manifest_scitd.txt"
  printf 'Gongsharma\nStephenson\nWu\n' > "${MANIFEST_FILE}"
}

echo "=== benchmark_bump_mem / benchmark_mem_ge ==="
assert_eq "bump 128G -> 256G" "256G" "$(benchmark_bump_mem 128G)"
assert_eq "bump 256G -> 512G" "512G" "$(benchmark_bump_mem 256G)"
assert_eq "bump 1T -> 2T" "2T" "$(benchmark_bump_mem 1T)"
assert_exit "bump unparseable fails" 1 benchmark_bump_mem foo
assert_exit "mem_ge 256G >= 128G" 0 benchmark_mem_ge 256G 128G
assert_exit "mem_ge 500G >= 500G" 0 benchmark_mem_ge 500G 500G
assert_exit "mem_ge 128G >= 256G" 1 benchmark_mem_ge 128G 256G
assert_exit "mem_ge 1T >= 500G" 0 benchmark_mem_ge 1T 500G
assert_exit "mem_ge 500G >= 1T" 1 benchmark_mem_ge 500G 1T
assert_exit "mem_ge unparseable fails" 1 benchmark_mem_ge foo 128G

echo "=== benchmark_pb_variant_names ==="
PB_VARIANTS="$(benchmark_pb_variant_names "${SCRIPT_DIR}/benchmark_hpc_utils.R")"
assert_eq "PB variant count" "6" "$(wc -w <<< "${PB_VARIANTS}" | tr -d ' ')"
assert_contains "first variant schvg2000" "schvg2000" "${PB_VARIANTS}"
assert_contains "last variant hvg3000" "hvg3000" "${PB_VARIANTS}"

echo "=== benchmark_wait_oom_retry: OK path with status file (OOM -> retry -> COMPLETED) ==="
setup_oom_scenario
make_manifest
SACCT_ROWS="11111:COMPLETED 11111_1:COMPLETED 11111_2:COMPLETED 11111_3:OUT_OF_MEMORY 22222:COMPLETED 22222_1:COMPLETED"
STATUS_FILE="${WATCHDOG_STATUS_DIR}/55555.status"
mkdir -p "${WATCHDOG_STATUS_DIR}"
benchmark_wait_oom_retry 11111 scitd resubmit_stub "${MANIFEST_FILE}" "${STATUS_FILE}"
assert_eq "exit 0 on OK" "0" "$?"
RESUBMIT_CAP="$(cat "${CAPTURE_DIR}/resubmit_capture.txt")"
assert_contains "resubmit got only the OOM'd dataset" "DS=Wu" "${RESUBMIT_CAP}"
assert_contains "resubmit mem doubled" "MEM=256G" "${RESUBMIT_CAP}"
STATUS_CONTENT="$(cat "${STATUS_FILE}")"
assert_contains "status STATE=OK" "STATE=OK" "${STATUS_CONTENT}"
assert_contains "status has JOB_REPORT for original id" "JOB_REPORT=scitd|11111|" "${STATUS_CONTENT}"
assert_contains "status has JOB_REPORT incl. FINAL retry id" "JOB_REPORT=scitd|22222|" "${STATUS_CONTENT}"
assert_eq "no notify_sync_status in status-file mode" "0" "$(notify_count)"
assert_eq "retry manifest reduced" "Wu" "$(cat "$(grep '^MANIFEST=' "${CAPTURE_DIR}/resubmit_capture.txt" | cut -d= -f2-)")"

echo "=== benchmark_wait_oom_retry: clamp 256G -> 500G (not 512G) ==="
setup_oom_scenario
export BENCHMARK_MEM="256G"
SACCT_ROWS="33333:COMPLETED 33333_1:OUT_OF_MEMORY 44444:COMPLETED 44444_1:COMPLETED"
echo "44444" > "${CAPTURE_DIR}/resubmit_retry_id.txt"
STATUS_FILE="${WATCHDOG_STATUS_DIR}/clamp.status"
benchmark_wait_oom_retry 33333 gloscope resubmit_stub "${MANIFEST_FILE}" "${STATUS_FILE}"
assert_eq "exit 0 after clamped retry" "0" "$?"
assert_contains "clamped mem is 500G" "MEM=500G" "$(cat "${CAPTURE_DIR}/resubmit_capture.txt")"
assert_contains "status STATE=OK after clamped retry" "STATE=OK" "$(cat "${STATUS_FILE}")"

echo "=== benchmark_wait_oom_retry: ceiling OOM -> FAIL with MaxRSS report ==="
setup_oom_scenario
export BENCHMARK_MEM="500G"
SACCT_ROWS="55555:COMPLETED 55555_1:OUT_OF_MEMORY"
SACCT_MEM_ROWS="55555_1:127505612K"
STATUS_FILE="${WATCHDOG_STATUS_DIR}/ceiling.status"
assert_exit "exit 1 at ceiling" 1 benchmark_wait_oom_retry 55555 scitd resubmit_stub "${MANIFEST_FILE}" "${STATUS_FILE}"
STATUS_CONTENT="$(cat "${STATUS_FILE}")"
assert_contains "STATE=FAIL" "STATE=FAIL" "${STATUS_CONTENT}"
assert_contains "FAIL_REASON mentions ceiling" "FAIL_REASON=OUT_OF_MEMORY tasks at the 500G ceiling" "${STATUS_CONTENT}"
assert_contains "report carries MaxRSS" "127505612K" "${STATUS_CONTENT}"
assert_eq "no email at ceiling (status-file mode)" "0" "$(notify_count)"

echo "=== benchmark_wait_oom_retry: non-OOM task -> FAIL ==="
setup_oom_scenario
SACCT_ROWS="66666:COMPLETED 66666_1:COMPLETED 66666_2:FAILED"
STATUS_FILE="${WATCHDOG_STATUS_DIR}/badtask.status"
assert_exit "exit 1 on non-OOM fail" 1 benchmark_wait_oom_retry 66666 mofa resubmit_stub "${MANIFEST_FILE}" "${STATUS_FILE}"
STATUS_CONTENT="$(cat "${STATUS_FILE}")"
assert_contains "STATE=FAIL" "STATE=FAIL" "${STATUS_CONTENT}"
assert_contains "FAIL_REASON mentions non-OOM" "FAIL_REASON=non-COMPLETED, non-OOM tasks" "${STATUS_CONTENT}"
assert_contains "report maps task 2 -> dataset" "2 ->" "${STATUS_CONTENT}"

echo "=== benchmark_wait_oom_retry: empty sacct -> FAIL ==="
setup_oom_scenario
SACCT_ROWS=""
STATUS_FILE="${WATCHDOG_STATUS_DIR}/emptysacct.status"
assert_exit "exit 1 on empty sacct" 1 benchmark_wait_oom_retry 77777 pseudobulk resubmit_stub "${MANIFEST_FILE}" "${STATUS_FILE}"
STATUS_CONTENT="$(cat "${STATUS_FILE}")"
assert_contains "STATE=FAIL" "STATE=FAIL" "${STATUS_CONTENT}"
assert_contains "FAIL_REASON mentions no states" "FAIL_REASON=sacct returned no states" "${STATUS_CONTENT}"

echo "=== benchmark_wait_oom_retry: unparseable mem -> FAIL ==="
setup_oom_scenario
export BENCHMARK_MEM="foo"
SACCT_ROWS="88888:COMPLETED 88888_1:OUT_OF_MEMORY"
STATUS_FILE="${WATCHDOG_STATUS_DIR}/badmem.status"
assert_exit "exit 1 on unparseable mem" 1 benchmark_wait_oom_retry 88888 pilot resubmit_stub "${MANIFEST_FILE}" "${STATUS_FILE}"
STATUS_CONTENT="$(cat "${STATUS_FILE}")"
assert_contains "STATE=FAIL" "STATE=FAIL" "${STATUS_CONTENT}"
assert_contains "FAIL_REASON mentions unparseable" "FAIL_REASON=unparseable BENCHMARK_MEM 'foo'" "${STATUS_CONTENT}"

echo "=== benchmark_wait_oom_retry: attempt cap -> FAIL ==="
setup_oom_scenario
export BENCHMARK_MEM="32G"
SACCT_ROWS="99999:COMPLETED 99999_1:OUT_OF_MEMORY 10000:COMPLETED 10000_1:OUT_OF_MEMORY 10001:COMPLETED 10001_1:OUT_OF_MEMORY 10002:COMPLETED 10002_1:OUT_OF_MEMORY 10003:COMPLETED 10003_1:OUT_OF_MEMORY"
printf '10000\n10001\n10002\n10003\n' > "${CAPTURE_DIR}/resubmit_retry_id.txt"
STATUS_FILE="${WATCHDOG_STATUS_DIR}/cap.status"
assert_exit "exit 1 at attempt cap" 1 benchmark_wait_oom_retry 99999 scpoli resubmit_stub "${MANIFEST_FILE}" "${STATUS_FILE}"
assert_contains "STATE=FAIL" "STATE=FAIL" "$(cat "${STATUS_FILE}")"
assert_contains "FAIL_REASON mentions attempt cap" "FAIL_REASON=exceeded 4 OOM retry attempts" "$(cat "${STATUS_FILE}")"

echo "=== benchmark_wait_oom_retry: regression without status file (email + JOB_REPORTS) ==="
setup_oom_scenario
SACCT_ROWS="11111:COMPLETED 11111_1:COMPLETED 11111_2:COMPLETED 11111_3:OUT_OF_MEMORY 22222:COMPLETED 22222_1:COMPLETED"
benchmark_wait_oom_retry 11111 scitd resubmit_stub "${MANIFEST_FILE}"
assert_eq "exit 0 on OK (no status file)" "0" "$?"
assert_eq "no email on OK" "0" "$(notify_count)"
assert_eq "JOB_REPORTS records only the final retry" "scitd|22222|00:05:00" "${JOB_REPORTS[0]}"
SACCT_ROWS="66666:COMPLETED 66666_1:FAILED"
assert_exit "exit 1 on non-OOM fail (no status file)" 1 benchmark_wait_oom_retry 66666 mofa resubmit_stub "${MANIFEST_FILE}"
assert_eq "email sent on fail without status file" "1" "$(notify_count)"
assert_contains "email subject says NOT synced" "NOT synced" "$(notify_last_subject)"

echo "=== benchmark_submit_watchdog: sbatch construction + id-only stdout ==="
setup_oom_scenario
mkdir -p "${WATCHDOG_STATUS_DIR}"
WATCHDOG_ID="$(benchmark_submit_watchdog \
  11111 scitd "${MANIFEST_FILE}" strict shared-cpu 1000 \
  "5_benchmark_r_scitd" "${REPO_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh" \
  --constraint=EPYC-7742 --cpus-per-task=16)"
assert_eq "watchdog id is the sbatch-returned id" "22222" "${WATCHDOG_ID}"
SBATCH_ARGS="$(head -1 "${CAPTURE_DIR}/sbatch_calls.txt")"
assert_contains "job name" "--job-name=benchmark_watchdog_scitd" "${SBATCH_ARGS}"
assert_contains "1 task / 1 cpu" "--ntasks=1 --cpus-per-task=1" "${SBATCH_ARGS}"
assert_contains "modest mem 2G" "--mem=2G" "${SBATCH_ARGS}"
assert_contains "WATCHDOG_TIME_LIMIT default 12h" "--time=12:00:00" "${SBATCH_ARGS}"
assert_contains "method partition (no pinned constraint)" "--partition=shared-cpu" "${SBATCH_ARGS}"
assert_contains "watchdog log naming" "5_benchmark_watchdog_scitd_%A.log" "${SBATCH_ARGS}"
assert_contains "array id forwarded" " 11111 scitd" "${SBATCH_ARGS}"
assert_contains "mode forwarded" " strict -- " "${SBATCH_ARGS}"
assert_contains "throttle forwarded" " shared-cpu 1000 5_benchmark_r_scitd" "${SBATCH_ARGS}"
assert_contains "worker script forwarded" "1.1_run_worker.sh" "${SBATCH_ARGS}"
assert_contains "flags preserved for retries" "--constraint=EPYC-7742 --cpus-per-task=16" "${SBATCH_ARGS}"

echo "=== benchmark_wait_watchdog: OK status -> pass + JOB_REPORTS merged ==="
setup_oom_scenario
mkdir -p "${WATCHDOG_STATUS_DIR}"
printf 'STATE=OK\nLABEL=scitd\nJOB_REPORT=scitd|11111|00:05:00\nJOB_REPORT=scitd|22222|00:05:00\n' \
  > "${WATCHDOG_STATUS_DIR}/55555.status"
SACCT_ROWS="55555:COMPLETED"
benchmark_wait_watchdog 55555 scitd
assert_eq "exit 0 on watchdog OK" "0" "$?"
assert_eq "JOB_REPORTS merged from file" "2" "${#JOB_REPORTS[@]}"
assert_eq "merged entry 1" "scitd|11111|00:05:00" "${JOB_REPORTS[0]}"
assert_eq "merged entry 2 (final retry)" "scitd|22222|00:05:00" "${JOB_REPORTS[1]}"
assert_eq "no email on OK" "0" "$(notify_count)"

echo "=== benchmark_wait_watchdog: FAIL status -> email + exit 1 ==="
setup_oom_scenario
mkdir -p "${WATCHDOG_STATUS_DIR}"
printf 'STATE=FAIL\nLABEL=scitd\nFAIL_REASON=OUT_OF_MEMORY tasks at the 500G ceiling\nREPORT=\nOOM per-task report\n' \
  > "${WATCHDOG_STATUS_DIR}/55556.status"
SACCT_ROWS="55556:COMPLETED"
assert_exit "exit 1 on watchdog FAIL" 1 benchmark_wait_watchdog 55556 scitd
assert_eq "email sent on watchdog FAIL" "1" "$(notify_count)"
assert_contains "email subject" "watchdog failed" "$(notify_last_subject)"
assert_contains "email body carries reason" "500G ceiling" "$(notify_last_body)"
assert_contains "email body carries report" "OOM per-task report" "$(notify_last_body)"

echo "=== benchmark_wait_watchdog: watchdog non-COMPLETED, no status file -> fail closed ==="
setup_oom_scenario
rm -f "${WATCHDOG_STATUS_DIR}/55557.status"
SACCT_ROWS="55557:FAILED"
SACCT_XEXIT="1:0"
assert_exit "exit 1 on lost watchdog" 1 benchmark_wait_watchdog 55557 scitd
assert_eq "email sent" "1" "$(notify_count)"
assert_contains "email subject" "watchdog lost" "$(notify_last_subject)"
assert_contains "body has sacct state" "State=FAILED" "$(notify_last_body)"
assert_contains "body has log pointer" "5_benchmark_watchdog_scitd_55557.log" "$(notify_last_body)"

echo "=== benchmark_wait_watchdog: COMPLETED but no status file -> fail closed ==="
setup_oom_scenario
rm -f "${WATCHDOG_STATUS_DIR}/55558.status"
SACCT_ROWS="55558:COMPLETED"
assert_exit "exit 1 on missing status file" 1 benchmark_wait_watchdog 55558 scitd
assert_eq "email sent" "1" "$(notify_count)"
assert_contains "body mentions missing status file" "without a status file" "$(notify_last_body)"

echo "=== watchdog_main (extracted): strict mode OOM escalation ==="
setup_oom_scenario
export SLURM_JOB_ID="77777"
WD_MANIFEST="${HPC_SCRATCH_DIR}/manifest_wd.txt"
printf 'Gongsharma\nStephenson\nWu\n' > "${WD_MANIFEST}"
SACCT_ROWS="11111:COMPLETED 11111_1:COMPLETED 11111_2:COMPLETED 11111_3:OUT_OF_MEMORY 22222:COMPLETED 22222_1:COMPLETED"
run_watchdog 11111 scitd "${WD_MANIFEST}" strict -- shared-cpu 1000 5_benchmark_r_scitd \
  "${REPO_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh" \
  --constraint=EPYC-7742 --cpus-per-task=16
assert_eq "watchdog strict exit 0" "0" "$(wd_rc)"
WD_STATUS="${WATCHDOG_STATUS_DIR}/77777.status"
assert_contains "status STATE=OK" "STATE=OK" "$(cat "${WD_STATUS}")"
assert_contains "status incl. retry id" "JOB_REPORT=scitd|22222|" "$(cat "${WD_STATUS}")"
RETRY_ARGS="$(head -1 "${CAPTURE_DIR}/sbatch_calls.txt")"
assert_contains "retry array size 1" "--array=1-1%1000" "${RETRY_ARGS}"
assert_contains "retry --mem doubled to 256G" "--mem=256G" "${RETRY_ARGS}"
assert_contains "retry partition" "--partition=shared-cpu" "${RETRY_ARGS}"
assert_contains "retry flags preserved" "--constraint=EPYC-7742 --cpus-per-task=16" "${RETRY_ARGS}"
assert_contains "retry log naming" "5_benchmark_r_scitd_%A_%a.log" "${RETRY_ARGS}"
assert_contains "retry worker script" "1.1_run_worker.sh" "${RETRY_ARGS}"
RETRY_MANIFEST="$(ls "${HPC_SCRATCH_DIR}"/benchmark_manifest_scitd_retry_*.txt 2>/dev/null | head -1)"
assert_eq "retry manifest reduced to OOM'd dataset" "Wu" "$(cat "${RETRY_MANIFEST}")"
assert_eq "METHOD exported for the retry workers" "scitd" "$(grep '^METHOD=' "${CAPTURE_DIR}/sbatch_env.txt" | head -1 | cut -d= -f2-)"
assert_eq "BENCHMARK_MANIFEST exported for the retry workers" "${RETRY_MANIFEST}" "$(grep '^MANIFEST=' "${CAPTURE_DIR}/sbatch_env.txt" | head -1 | cut -d= -f2-)"

echo "=== watchdog_main (extracted): strict clamp 256G -> 500G ==="
setup_oom_scenario
export SLURM_JOB_ID="77778"
export BENCHMARK_MEM="256G"
SACCT_ROWS="33333:COMPLETED 33333_1:OUT_OF_MEMORY 44444:COMPLETED 44444_1:COMPLETED"
echo "44444" > "${CAPTURE_DIR}/sbatch_next.txt"
run_watchdog 33333 gloscope "${WD_MANIFEST}" strict -- shared-cpu 1000 5_benchmark_r_gloscope \
  "${REPO_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
assert_eq "watchdog strict (clamp) exit 0" "0" "$(wd_rc)"
assert_contains "retry --mem clamped to 500G" "--mem=500G" "$(head -1 "${CAPTURE_DIR}/sbatch_calls.txt")"
assert_contains "status STATE=OK" "STATE=OK" "$(cat "${WATCHDOG_STATUS_DIR}/77778.status")"

echo "=== watchdog_main (extracted): soft-gate, all variants present -> OK without gate ==="
setup_oom_scenario
export SLURM_JOB_ID="88888"
mkdir -p "${HPC_SCRATCH_DIR}/benchmark/pseudobulks"
for ds in Gongsharma Stephenson Wu; do
  for v in ${PB_VARIANTS}; do
    touch "${HPC_SCRATCH_DIR}/benchmark/pseudobulks/${ds}_pseudobulk_${v}.rds"
  done
done
SACCT_ROWS="11111:COMPLETED 11111_1:COMPLETED 11111_2:COMPLETED 11111_3:FAILED"
run_watchdog 11111 prepare_pseudobulk "${WD_MANIFEST}" soft-gate -- shared-cpu 1000 5_benchmark_r_prepare_pseudobulk \
  "${REPO_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
assert_eq "soft-gate exit 0" "0" "$(wd_rc)"
assert_eq "no resubmit in soft-gate pass" "0" "$(wc -l < "${CAPTURE_DIR}/sbatch_calls.txt" | tr -d ' ')"
WD_STATUS="${WATCHDOG_STATUS_DIR}/88888.status"
assert_contains "status STATE=OK" "STATE=OK" "$(cat "${WD_STATUS}")"
assert_contains "status JOB_REPORT for prep array" "JOB_REPORT=prepare_pseudobulk|11111|" "$(cat "${WD_STATUS}")"

echo "=== watchdog_main (extracted): soft-gate, variants missing -> strict gate ==="
setup_oom_scenario
export SLURM_JOB_ID="88889"
rm -rf "${HPC_SCRATCH_DIR}/benchmark/pseudobulks"
mkdir -p "${HPC_SCRATCH_DIR}/benchmark/pseudobulks"
touch "${HPC_SCRATCH_DIR}/benchmark/pseudobulks/Gongsharma_pseudobulk_schvg2000.rds"
SACCT_ROWS="11111:COMPLETED 11111_1:COMPLETED 11111_2:COMPLETED 11111_3:COMPLETED"
run_watchdog 11111 prepare_pseudobulk "${WD_MANIFEST}" soft-gate -- shared-cpu 1000 5_benchmark_r_prepare_pseudobulk \
  "${REPO_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
assert_eq "soft-gate fallback exit 0 (all COMPLETED)" "0" "$(wd_rc)"
assert_contains "status STATE=OK via strict gate" "STATE=OK" "$(cat "${WATCHDOG_STATUS_DIR}/88889.status")"

echo "=== watchdog_main (extracted): soft-gate under --force -> strict gate ==="
setup_oom_scenario
export SLURM_JOB_ID="88890"
export FORCE_BENCHMARK=1
SACCT_ROWS="11111:COMPLETED 11111_1:COMPLETED 11111_2:COMPLETED 11111_3:COMPLETED"
run_watchdog 11111 prepare_pseudobulk "${WD_MANIFEST}" soft-gate -- shared-cpu 1000 5_benchmark_r_prepare_pseudobulk \
  "${REPO_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
assert_eq "soft-gate --force exit 0 (all COMPLETED)" "0" "$(wd_rc)"
assert_contains "status STATE=OK" "STATE=OK" "$(cat "${WATCHDOG_STATUS_DIR}/88890.status")"

echo "=== watchdog_main (extracted): unknown mode -> exit 1, no status file ==="
setup_oom_scenario
export SLURM_JOB_ID="77779"
assert_exit "unknown mode exit 1" 1 watchdog_main 11111 scitd "${WD_MANIFEST}" bogus -- shared-cpu 1000 5_benchmark_r_scitd \
  "${REPO_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
assert_eq "no status file on unknown mode" "" "$(cat "${WATCHDOG_STATUS_DIR}/77779.status" 2>/dev/null)"

# ---------------------------------------------------------------------------
echo
if [[ ${FAILURES} -eq 0 ]]; then
  echo "ALL ${TESTS} ASSERTIONS PASSED"
  exit 0
else
  echo "${FAILURES} OF ${TESTS} ASSERTIONS FAILED"
  exit 1
fi
