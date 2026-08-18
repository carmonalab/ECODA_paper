#!/bin/bash
# ==============================================================================
# download_datasets_hpc.sh -- Phase-5 downloads: HPC scratch -> NAS sync tail
# ==============================================================================
# Login-node submitter for the 9 PILOT-GM-VAE study dataset downloads
# (onboarding plan .kilo/plans/1786899069337-onboard-new-datasets-phase5.md).
#
# ROUTE (user decision 2026-08-17 -- replaces the Mac->NAS direct route): the
# Mac NAS mount is unstable and the Mac is short on disk, so downloads land in
# BeeGFS scratch FIRST, then a login-node tail rsyncs them to the NAS folder
# ${NAS_SC_DIR}/JooM_2025_41097818/output/ (where the onboarding check
# notebooks read them from the Mac). The original Mac->NAS
# `download_datasets.sh` stays as a NAS-mounted fallback ONLY.
#
# FLOW
#   * Smoke-tests compute-node egress (srun curl to Zenodo + CellxGene on
#     debug-cpu); on failure it falls back to quiet login-node downloads
#     (`nice -n 19` + `curl --limit-rate`, the documented high-bandwidth
#     transfer path) instead of submitting an array.
#   * Submits ONE SLURM array job -- one task per requested key, concurrency
#     throttled (DOWNLOAD_ARRAY_THROTTLE, default 3) -- onto shared-cpu. Each
#     task runs run_download_worker.sh: `curl -L -C -` resumable + verified per
#     key (Zenodo md5s / CellxGene size-verified via HEAD content-length -- no
#     `.h5ad.md5` sidecar exists; computed md5 recorded as informational) +
#     selective tar extraction + tar deletion, into ${HPC_SCRATCH_DIR}/_downloads/
#     (BeeGFS scratch: 1.1 PB, no per-user size quota, NOT $HOME itself).
#   * Waits + sacct-gates every task (fail-closed: any non-COMPLETED task, or
#     any key without a recorded STATUS=OK, blocks the NAS sync), then
#     rsyncs to NAS -- ONLY the files recorded in a STATUS=OK status file
#     (leftover tars / _tar_tmp.* leaks / partial h5ads of failed keys are
#     cleaned from scratch first and never synced; scratch copy KEPT for later
#     pipeline use) -- verifies the NAS md5s against the worker records
#     (Zenodo files: expected md5; CellxGene files: worker-recorded
#     MD5_RECORDED; tar-extracted files: worker-recorded MD5_<file>=), and
#     appends the per-key report (+ tar
#     listings + the `--sync-only <job-id>` resume command) to
#     notebooks/dataset_onboarding/download_log.md (commit from the Mac).
#   * Concurrency: a flock on ${DOWNLOAD_DIR}/.submit.lock serializes
#     submitter runs (multiple login terminals), and the array path refuses to
#     submit while an `onboard_download` job is already queued -- two
#     concurrent runs would corrupt the shared partial files via interleaved
#     `curl -C -` and double-sync to the NAS. Never run two copies.
#   * Cancellation is safe at any point (Ctrl-C / scancel); re-running the
#     same command resumes partial files via `curl -C -`.
#
# Usage (run from the HPC login node; `git pull` the HPC clone first):
#   ./notebooks/dataset_onboarding/download_datasets_hpc.sh                      # all 8 keys
#   ./notebooks/dataset_onboarding/download_datasets_hpc.sh --only breast        # one key
#   ./notebooks/dataset_onboarding/download_datasets_hpc.sh --login-node         # login-node fallback
#   ./notebooks/dataset_onboarding/download_datasets_hpc.sh --sync-only <job-id> # resume: gate + NAS sync tail only
#
# Env overrides: DOWNLOAD_ARRAY_THROTTLE (default 3), DOWNLOAD_LIMIT_RATE
# (login-node mode, default 2m), DOWNLOAD_SKIP_NAS_MD5=1 (skip the NAS-side
# md5 re-verification -- doubles read traffic, several minutes on the tail).
# ==============================================================================
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../../src/slurm_config.sh"
cd "${PROJECT_ROOT}"

KEYS=(alzheimer breast covid_lupus_tar diabetes kidney lung_tar myocardial parkinson)
DOWNLOAD_DIR="${HPC_SCRATCH_DIR}/_downloads"
STATUS_DIR="${DOWNLOAD_DIR}/_status"
NAS_DEST="${NAS_SC_DIR}/JooM_2025_41097818/output"
ONBOARD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER="${ONBOARD_DIR}/run_download_worker.sh"
LOG_FILE="${ONBOARD_DIR}/download_log.md"
THROTTLE="${DOWNLOAD_ARRAY_THROTTLE:-3}"

ONLY=""
SYNC_ONLY=""
LOGIN_NODE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)
      ONLY="${2:-}"
      shift 2
      ;;
    --sync-only)
      SYNC_ONLY="${2:-}"
      if [[ -z "${SYNC_ONLY}" ]]; then
        echo "ERROR: --sync-only requires a job id."
        exit 1
      fi
      shift 2
      ;;
    --login-node)
      LOGIN_NODE=1
      shift
      ;;
    *)
      echo "ERROR: unknown argument $1 (accepted: --only <key> | --sync-only <id> | --login-node)." >&2
      exit 1
      ;;
  esac
done

requested_keys() {
  if [[ -n "${ONLY}" ]]; then
    local ok=0 k
    for k in "${KEYS[@]}"; do
      [[ "${k}" == "${ONLY}" ]] && ok=1
    done
    if [[ ${ok} -ne 1 ]]; then
      echo "ERROR: --only '${ONLY}' invalid. Valid keys: ${KEYS[*]}" >&2
      exit 1
    fi
    printf '%s\n' "${ONLY}"
  else
    printf '%s\n' "${KEYS[@]}"
  fi
}
REQ_KEYS=()
while IFS= read -r k; do REQ_KEYS+=("${k}"); done < <(requested_keys)

# ---------------------------------------------------------------------------
# Concurrency guard: two concurrent submitters (or a submitter racing an array
# from an interrupted earlier run) would make two tasks `curl -C -` into the
# same partial file (interleaved appends corrupt it) and double-sync to the
# NAS. The flock covers concurrent login-node runs; the squeue check below
# covers an orphaned array from a dropped submitter session.
# ---------------------------------------------------------------------------
mkdir -p "${DOWNLOAD_DIR}"
LOCK_FILE="${DOWNLOAD_DIR}/.submit.lock"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  echo "ERROR: another download_datasets_hpc.sh run holds ${LOCK_FILE}; refusing to run concurrently (concurrent runs corrupt the shared partial files)." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Log: per-key status + md5s + tar listings + resume command -> download_log.md
# (defined early: called from every exit path below)
# ---------------------------------------------------------------------------
ARRAY_JOB_ID="${ARRAY_JOB_ID:-}"
RUN_MODE="${RUN_MODE:-}"
append_log() {
  local header="$1"
  {
    printf '\n### %s\n\n- **HPC download run:** %s\n' "$(date '+%Y-%m-%d %H:%M')" "${header}"
    printf -- '- **Job:** %s\n- **Mode:** %s\n- **Keys:** %s\n' \
      "${ARRAY_JOB_ID}" "${RUN_MODE}" "${REQ_KEYS[*]}"
    for k in "${KEYS[@]}"; do
      local sf="${STATUS_DIR}/${k}.status"
      if [[ -n "${ONLY}" && "${k}" != "${ONLY}" ]]; then
        continue
      fi
      if [[ -f "${sf}" ]]; then
        local st f sz ma me note
        st="$(grep '^STATUS=' "${sf}" | cut -d= -f2- | head -1 || true)"
        f="$(grep '^FILE=' "${sf}" | cut -d= -f2- | head -1 || true)"
        sz="$(grep '^SIZE=' "${sf}" | cut -d= -f2- | head -1 || true)"
        ma="$(grep '^MD5_ACTUAL=' "${sf}" | cut -d= -f2- | head -1 || true)"
        me="$(grep '^MD5_EXPECTED=' "${sf}" | cut -d= -f2- | head -1 || true)"
        note="$(grep -E '^(SKIPPED|FILES|ERROR|VERIFY|MD5_RECORDED)=' "${sf}" | cut -d= -f2- | paste -sd ', ' - || true)"
        printf -- '- **Key:** %s **Status:** %s\n  - File: %s (%s)\n  - md5: %s (expected %s)\n' \
          "${k}" "${st:-?}" "${f:-?}" "${sz:-?}" "${ma:-n/a}" "${me:-n/a}"
        if [[ -n "${note}" ]]; then
          printf '  - Note: %s\n' "${note}"
        fi
      else
        printf -- '- **Key:** %s **Status:** n/a (no status file)\n' "${k}"
      fi
      local tl="${STATUS_DIR}/${k}.tar_listing"
      if [[ -f "${tl}" ]]; then
        printf '  - Tar members:\n'
        sed 's/^/      /' "${tl}"
      fi
    done
    printf -- '- **Resume:** re-run the same command, or `--sync-only %s` for the gate+sync tail (repeat the original `--only` flags)\n' "${ARRAY_JOB_ID}"
  } >>"${LOG_FILE}"
}

echo "=== Phase-5 dataset downloads (keys: ${REQ_KEYS[*]}) ==="
echo "Download dir (BeeGFS scratch): ${DOWNLOAD_DIR}"
echo "NAS destination: ${NAS_DEST}"
echo ""
echo "Scratch free space (${HOME}/scratch):"
df -h "${HOME}/scratch" || true

# ---------------------------------------------------------------------------
# Download stage: compute-node array (default) or login-node fallback
# ---------------------------------------------------------------------------
if [[ -n "${SYNC_ONLY}" ]]; then
  RUN_MODE="sync-only"
  ARRAY_JOB_ID="${SYNC_ONLY}"
  echo "=== Sync-only resume mode: job ${ARRAY_JOB_ID} (no submission, no smoke test) ==="
elif [[ ${LOGIN_NODE} -eq 1 ]]; then
  RUN_MODE="login-node"
  ARRAY_JOB_ID="login-node-run"
  echo "=== Login-node fallback mode (documented high-bandwidth transfer path) ==="
else
  echo "=== Smoke-testing compute-node egress (debug-cpu) ==="
  # NOTE: probe REAL object URLs -- a HEAD on the CellxGene bucket root
  # returns 403 (S3 denies bucket-root HEAD), which falsely triggered the
  # login-node fallback on 2026-08-17. The SEA-AD .h5ad object is the
  # canonical CellxGene probe. srun stderr is NOT swallowed: a genuine
  # srun/partition failure must be visible in the output.
  smoke_host() {
    local url="$1" out code
    out="$(srun -p debug-cpu --time=00:05:00 curl -sIL -o /dev/null -w '%{http_code}' "${url}" 2>&1 || true)"
    code="$(printf '%s\n' "${out}" | tail -1)"
    [[ "${code}" =~ ^(2|3)[0-9][0-9]$ ]]
  }
  SMOKE_URLS=(https://zenodo.org \
    https://datasets.cellxgene.cziscience.com/c2b49431-9288-4d94-8ca5-f6723b72217e.h5ad)
  SMOKE_FAIL=0
  for u in "${SMOKE_URLS[@]}"; do
    if smoke_host "${u}"; then
      echo "OK: compute nodes reach ${u}"
    else
      echo "FAIL: compute nodes could not reach ${u} (see srun/curl output above)" >&2
      SMOKE_FAIL=1
    fi
  done
  if [[ ${SMOKE_FAIL} -eq 0 ]]; then
    echo "OK: compute-node egress confirmed -- submitting the download array."
    RUN_MODE="array"
    ARRAY_JOB_ID=""
  else
    echo "WARNING: compute-node egress check failed -- falling back to login-node downloads (nice + --limit-rate)." >&2
    RUN_MODE="login-node"
    ARRAY_JOB_ID="login-node-run"
  fi
fi

if [[ "${RUN_MODE}" == "array" ]]; then
  SPEC_ARGS=()
  for k in "${REQ_KEYS[@]}"; do
    for i in "${!KEYS[@]}"; do
      if [[ "${KEYS[$i]}" == "${k}" ]]; then SPEC_ARGS+=("$((i + 1))"); fi
    done
  done
  ARRAY_SPEC="$(IFS=,; echo "${SPEC_ARGS[*]}")"
  if squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -qx 'onboard_download'; then
    echo "ERROR: an onboard_download array job is already in the queue; refusing to submit a second one (interleaved curl -C - corrupts the shared partial files)." >&2
    echo "Wait for it to finish, or resume its sync tail with --sync-only <job-id>." >&2
    exit 1
  fi
  mkdir -p "${LOGS_DIR}"
  export DOWNLOAD_KEYS="${KEYS[*]}"
  export DOWNLOAD_LIMIT_RATE=""
  echo "=== Submitting download array job (tasks: ${ARRAY_SPEC}, throttle ${THROTTLE}) ==="
  SUBMIT_MSG="$(sbatch \
    --array="${ARRAY_SPEC}%${THROTTLE}" \
    --partition="${SLURM_PARTITION_BENCHMARK_CPU}" \
    --cpus-per-task=1 --mem=4G --time=12:00:00 \
    --output="${LOGS_DIR}/onboard_download_%A_%a.log" \
    --error="${LOGS_DIR}/onboard_download_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    --mail-type=END,FAIL \
    "${WORKER}")"
  ARRAY_JOB_ID="$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')"
  echo "Array Job ID allocated: ${ARRAY_JOB_ID}"
fi

if [[ "${RUN_MODE}" == "login-node" ]]; then
  export DOWNLOAD_KEYS="${KEYS[*]}"
  export DOWNLOAD_LIMIT_RATE="${DOWNLOAD_LIMIT_RATE:-2m}"
  FAILED_KEYS=""
  echo "=== Downloading on the login node (nice -n 19, --limit-rate ${DOWNLOAD_LIMIT_RATE}) ==="
  for k in "${REQ_KEYS[@]}"; do
    echo ""
    echo ">>> ${k}"
    if ! nice -n 19 bash "${WORKER}" --key "${k}"; then
      FAILED_KEYS="${FAILED_KEYS} ${k}"
    fi
  done
  if [[ -n "${FAILED_KEYS}" ]]; then
    echo "ERROR: login-node downloads failed for:${FAILED_KEYS}" >&2
    append_log "FAILED -- NOT synced to NAS (login-node download stage incomplete)"
    exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Monitor + gate (array path: squeue wait; array AND sync-only: sacct gate with
# fail-closed unknown/purged handling per AGENTS.md)
# ---------------------------------------------------------------------------
if [[ "${RUN_MODE}" == "array" ]]; then
  echo ""
  echo "=== Monitoring job completion ==="
  while squeue -u "$USER" -h -o "%A" 2>/dev/null | grep -qx "${ARRAY_JOB_ID}"; do
    sleep 60
  done
  echo "Array Job ${ARRAY_JOB_ID} left the scheduler."
fi

if [[ "${RUN_MODE}" == "array" || "${RUN_MODE}" == "sync-only" ]]; then
  echo ""
  echo "Checking sacct terminal states for job ${ARRAY_JOB_ID} (bounded, max 15 min)..."
  TAIL_ITER=0
  while ((TAIL_ITER < 180)); do
    STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
    if [[ -n "${STATES//[[:space:]]/}" ]] \
      && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<<"${STATES}"; then
      break
    fi
    sleep 5
    TAIL_ITER=$((TAIL_ITER + 1))
  done
  STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
  if [[ -z "${STATES//[[:space:]]/}" ]]; then
    echo "ERROR: sacct returned no states for Array Job ${ARRAY_JOB_ID}; NOT syncing (unknown/purged id)." >&2
    append_log "FAILED -- NOT synced (sacct records purged/unknown for ${ARRAY_JOB_ID})"
    echo "Resume: --sync-only ${ARRAY_JOB_ID}" >&2
    exit 1
  fi
  if grep -qvE '^ *COMPLETED *$' <<<"${STATES}"; then
    echo "ERROR: Array Job ${ARRAY_JOB_ID} had non-COMPLETED tasks; NOT syncing." >&2
    sacct -j "${ARRAY_JOB_ID}" --format=JobID,State,ExitCode
    append_log "FAILED -- NOT synced (non-COMPLETED tasks; partial files in ${DOWNLOAD_DIR})"
    echo "Partial files persist in ${DOWNLOAD_DIR}; re-run the same command (resumes via curl -C -) or fix + --sync-only ${ARRAY_JOB_ID}." >&2
    exit 1
  fi
fi

# Per-key status gate (covers both array and login-node runs).
echo ""
echo "=== Download status per key (${STATUS_DIR}/) ==="
STATUS_FAILED=()
for k in "${REQ_KEYS[@]}"; do
  sf="${STATUS_DIR}/${k}.status"
  if [[ -f "${sf}" ]] && grep -q '^STATUS=OK$' "${sf}"; then
    f="$(grep '^FILE=' "${sf}" | cut -d= -f2- | head -1 || true)"
    sz="$(grep '^SIZE=' "${sf}" | cut -d= -f2- | head -1 || true)"
    echo "  OK   ${k}  (${f:-?}, ${sz:-?})"
  else
    echo "  FAIL ${k}"
    STATUS_FAILED+=("${k}")
  fi
done
if [[ ${#STATUS_FAILED[@]} -gt 0 ]]; then
  echo "ERROR: ${#STATUS_FAILED[@]} key(s) not OK (${STATUS_FAILED[*]}); NOT syncing to NAS." >&2
  append_log "FAILED -- NOT synced to NAS (download stage incomplete; partial files in ${DOWNLOAD_DIR})"
  echo "Partial files persist in ${DOWNLOAD_DIR}; re-run the same command to resume." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# NAS sync tail: rsync scratch -> NAS, keep scratch copy, delete tar archives
# ---------------------------------------------------------------------------
echo ""
echo "=== Syncing _downloads -> NAS (${NAS_DEST}) ==="
if [[ ! -d "$(dirname "${NAS_DEST}")" ]]; then
  echo "ERROR: NAS mount missing ($(dirname "${NAS_DEST}")) -- is the NAS reachable from the login node?" >&2
  append_log "FAILED -- NAS unreachable at sync (files remain in scratch: ${DOWNLOAD_DIR})"
  exit 1
fi
mkdir -p "${NAS_DEST}"

# Cleanup BEFORE the rsync: leftover tar.gz archives + extraction temp dirs are
# removed from scratch and never reach the NAS (they are also excluded below as
# a belt-and-braces guard).
find "${DOWNLOAD_DIR}" -maxdepth 1 -name '*.tar.gz' -delete
find "${DOWNLOAD_DIR}" -maxdepth 1 -type d -name '_tar_tmp.*' -exec rm -rf {} + 2>/dev/null || true

# Sync ONLY the files recorded in a STATUS=OK status file -- never leftover
# tars, _tar_tmp.* leaks, or partial h5ads of failed keys.
SYNC_LIST="${STATUS_DIR}/.sync_list"
: > "${SYNC_LIST}"
for sf in "${STATUS_DIR}"/*.status; do
  [[ -f "${sf}" ]] || continue
  grep -q '^STATUS=OK$' "${sf}" || continue
  f="$(grep '^FILE=' "${sf}" | cut -d= -f2- | head -1 || true)"
  [[ -n "${f}" ]] && printf '%s\n' "${f}" >> "${SYNC_LIST}"
  for nf in $(grep '^FILES=' "${sf}" 2>/dev/null | cut -d= -f2- | tr ',' ' '); do
    [[ -n "${nf}" ]] && printf '%s\n' "${nf}" >> "${SYNC_LIST}"
  done
done
sort -u "${SYNC_LIST}" -o "${SYNC_LIST}"
rsync -rlptDv \
  --exclude '_status/' \
  --exclude '_tar_tmp.*' \
  --include '*/' \
  --include-from "${SYNC_LIST}" \
  --exclude '*' \
  "${DOWNLOAD_DIR}/" "${NAS_DEST}/"

# Verify NAS md5s against the worker records (belt-and-braces; rsync already
# guarantees byte-identical copies). Skippable: doubles read traffic.
if [[ "${DOWNLOAD_SKIP_NAS_MD5:-0}" == "1" ]]; then
  echo "Skipping NAS md5 verification (DOWNLOAD_SKIP_NAS_MD5=1)."
fi
NAS_VERIFY_FAIL=0
for k in "${REQ_KEYS[@]}"; do
  sf="${STATUS_DIR}/${k}.status"
  # Tar keys: the canonical outputs carry worker-recorded md5s (MD5_<file>=) --
  # verify the NAS copy against them; fall back to scratch==NAS identity for
  # statuses written by older workers without recorded extracted md5s.
  files_list="$(grep '^FILES=' "${sf}" 2>/dev/null | cut -d= -f2- | tr ',' ' ' || true)"
  if [[ -n "${files_list}" ]]; then
    for nf in ${files_list}; do
      if [[ "${DOWNLOAD_SKIP_NAS_MD5:-0}" == "1" ]]; then
        echo "  ${nf}: NAS size $(stat -c %s "${NAS_DEST}/${nf}" 2>/dev/null || echo '?') (md5 verify skipped)"
        continue
      fi
      exp_md5="$(grep "^MD5_${nf}=" "${sf}" 2>/dev/null | cut -d= -f2- | head -1 || true)"
      if [[ -n "${exp_md5}" ]]; then
        na="$(md5sum "${NAS_DEST}/${nf}" | awk '{print $1}')"
        if [[ "${na}" == "${exp_md5}" ]]; then
          echo "  OK   ${nf} (NAS md5 ${na})"
        else
          echo "  FAIL ${nf}: NAS md5 ${na} != worker-verified ${exp_md5}" >&2
          NAS_VERIFY_FAIL=1
        fi
      else
        sc="$(md5sum "${DOWNLOAD_DIR}/${nf}" | awk '{print $1}')"
        na="$(md5sum "${NAS_DEST}/${nf}" | awk '{print $1}')"
        if [[ "${sc}" == "${na}" ]]; then
          echo "  OK   ${nf} (scratch == NAS md5 ${na})"
        else
          echo "  FAIL ${nf}: scratch md5 ${sc} != NAS md5 ${na}" >&2
          NAS_VERIFY_FAIL=1
        fi
      fi
    done
    continue
  fi
  f="$(grep '^FILE=' "${sf}" 2>/dev/null | cut -d= -f2- | head -1 || true)"
  exp="$(grep '^MD5_ACTUAL=' "${sf}" 2>/dev/null | cut -d= -f2- | head -1 || true)"
  [[ -z "${f}" || -z "${exp}" ]] && continue
  if [[ "${DOWNLOAD_SKIP_NAS_MD5:-0}" == "1" ]]; then
    echo "  ${f}: NAS size $(stat -c %s "${NAS_DEST}/${f}" 2>/dev/null || echo '?') (md5 verify skipped)"
    continue
  fi
  got="$(md5sum "${NAS_DEST}/${f}" | awk '{print $1}')"
  if [[ "${got}" == "${exp}" ]]; then
    echo "  OK   ${f} (NAS md5 ${got})"
  else
    echo "  FAIL ${f}: NAS md5 ${got} != expected ${exp}" >&2
    NAS_VERIFY_FAIL=1
  fi
done
if [[ ${NAS_VERIFY_FAIL} -ne 0 ]]; then
  echo "ERROR: NAS md5 verification failed; re-run the same command (idempotent rsync) and check the NAS." >&2
  append_log "FAILED -- NAS md5 verification failed (downloads OK in scratch: ${DOWNLOAD_DIR})"
  exit 1
fi

append_log "completed -- synced to NAS (scratch copy kept in ${DOWNLOAD_DIR})"

echo ""
echo "=== Downloads complete + synced to NAS ==="
echo "  NAS:    ${NAS_DEST}"
echo "  Scratch: ${DOWNLOAD_DIR}  (kept for later pipeline use)"
echo "  Log:    ${LOG_FILE}  (commit from the Mac)"
echo ""
echo "Next: check the files from the Mac:"
echo "  ls -lh /Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets/JooM_2025_41097818/output/"