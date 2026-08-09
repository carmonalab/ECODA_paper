#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Dispatcher for dataset-specific preprocessing steps (src/2_dataset_specific_preprocessing/).
#
# Submits every per-step sbatch script in this folder (1.*_submit_*.sh) IN
# PARALLEL, waits for all of them, then reports per-job state via sacct and
# exits non-zero if any step failed.
#
# IMPORTANT: Steps are submitted in parallel and therefore MUST be mutually
# independent. If a future step depends on another step, it must be submitted
# after the wait loop (or via `--dependency`) — do NOT rely on submission
# order alone.
#
# Run AFTER src/1_stage_data/1_stage_data.sh and BEFORE
# src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh.
#
# Usage:
#   ./1_submit_hpc.sh                    # submit + wait + report
#   ./1_submit_hpc.sh --sync-only <id1,id2,...>   # resume: skip submission, re-check the given jobs
#
# NOTE: This dispatcher never syncs to NAS (dataset-specific preprocessing
# writes to scratch only; the preprocess array syncs afterwards), so
# --sync-only here only re-runs the wait + sacct gate + summary — no status
# email is sent.
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

STEP_SCRIPTS=( "${SCRIPT_DIR}"/1.*_submit_*.sh )

SYNC_ONLY_IDS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
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

if [[ ${#STEP_SCRIPTS[@]} -eq 0 && -z "${SYNC_ONLY_IDS}" ]]; then
  echo "ERROR: No step scripts found in ${SCRIPT_DIR}."
  exit 1
fi

mkdir -p "${LOGS_DIR}"

JOB_IDS=()
if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Sync-only resume mode: jobs ${SYNC_ONLY_IDS} (no submission) ==="
  IFS=',' read -r -a JOB_IDS <<< "${SYNC_ONLY_IDS}"
else
for step_script in "${STEP_SCRIPTS[@]}"; do
  echo "=== Submitting $(basename "${step_script}") ==="
  # --output/--error/--partition passed on the sbatch command line (not
  # #SBATCH lines): SLURM directives do not expand environment variables.
  STEP_LOG_STEM="${LOGS_DIR}/$(basename "${step_script}" .sh)_%j"
  JOB_IDS+=("$(sbatch --parsable \
      --partition="${SLURM_PARTITION}" \
      --output="${STEP_LOG_STEM}.log" \
      --error="${STEP_LOG_STEM}.err" \
      "${step_script}")")
done
fi

echo "Submitted jobs: ${JOB_IDS[*]}"

# Wait for all jobs to leave the queue (squeue poll, 60s interval; exit
# codes ignored — the per-job sacct COMPLETED check below is the
# authoritative gate). scontrol has no plain `wait` command (only
# `wait_job`, which waits for node-ready — not completion), so poll squeue
# for the exact job id (`-o %A` prints the array master id for every task,
# or the job id for plain jobs).
for job_id in "${JOB_IDS[@]}"; do
  while squeue -u "$USER" -h -o "%A" 2>/dev/null | grep -qx "${job_id}"; do
    sleep 60
  done
  echo "Job ${job_id} left the scheduler."
done

# sacct may lag a few seconds behind jobs leaving the scheduler; poll
# (bounded) until every job's master state is terminal instead of a blind
# fixed sleep (same -X master-only granularity as the COMPLETED check below).
# The 180-iteration cap (15 min) plus a 60-iteration grace window (5 min)
# covers pathological SlurmDBD accounting lag (scheduler said done, sacct
# still reports RUNNING); the per-job COMPLETED gate below is unchanged.
echo "Waiting for sacct to record terminal states (bounded, max 20 min)..."
TAIL_ITER=0
while (( TAIL_ITER < 180 )); do  # max 15 min at 5s
  settled=1
  for job_id in "${JOB_IDS[@]}"; do
    state="$(sacct -j "${job_id}" -X -n -o "State" 2>/dev/null | tr -d ' \n' || true)"
    if [[ -z "${state}" ]] || [[ "${state}" =~ PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING ]]; then
      settled=0
    fi
  done
  if [[ ${settled} -eq 1 ]]; then
    break
  fi
  sleep 5
  TAIL_ITER=$((TAIL_ITER + 1))
done
if (( TAIL_ITER >= 180 )); then
  echo "WARNING: sacct still reporting non-terminal states after 15 min; extending wait by a 5 min grace window..." >&2
  TAIL_ITER=0
  while (( TAIL_ITER < 60 )); do
    settled=1
    for job_id in "${JOB_IDS[@]}"; do
      state="$(sacct -j "${job_id}" -X -n -o "State" 2>/dev/null | tr -d ' \n' || true)"
      if [[ -z "${state}" ]] || [[ "${state}" =~ PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING ]]; then
        settled=0
      fi
    done
    if [[ ${settled} -eq 1 ]]; then
      break
    fi
    sleep 5
    TAIL_ITER=$((TAIL_ITER + 1))
  done
fi

echo "=== Job summary ==="
sacct -j "$(IFS=,; echo "${JOB_IDS[*]}")" --format=JobID,JobName,State,ExitCode

# Exit non-zero if any step failed
failed=0
for job_id in "${JOB_IDS[@]}"; do
  state=$(sacct -j "${job_id}" -X -n -o "State" 2>/dev/null | tr -d ' \n')
  if [[ "${state}" != "COMPLETED" ]]; then
    echo "ERROR: Job ${job_id} ended in state '${state}'."
    failed=1
  fi
done

if [[ ${failed} -ne 0 ]]; then
  exit 1
fi

echo "All dataset-specific preprocessing steps completed."
