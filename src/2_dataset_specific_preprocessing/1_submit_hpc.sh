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
# independent. The one exception is wired explicitly via --dependency: the
# GongSharma cap step (1.1_submit_gongsharma.sh, when present) is submitted
# first and the CombinedPBMC step (1.2) is gated behind it with
# --dependency=afterok (1.1 overwrites the staged SoundLife h5ads IN PLACE,
# which 1.2 reads in backed mode — a race would nondeterminize the CombinedPBMC
# dataset). For any future step dependency, use the same pattern (submit first
# + --dependency) — do NOT rely on submission order alone.
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
# --sync-only here only re-runs the wait + sacct gate + summary. A best-effort
# completion email (per-step state + elapsed, via
# src/utils/bash/sync_status_email.sh) is sent at the end — success or
# failure — and --mail-user="${USER_EMAIL}" is passed on every sbatch call so
# Slurm's own END/FAIL job emails reach the user instead of the cluster
# default (the step scripts only carry #SBATCH --mail-type=END,FAIL).
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/sync_status_email.sh"
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
  # GongSharma cap step (1.1) must finish before the CombinedPBMC step (1.2):
  # 1.1 overwrites the staged SoundLife h5ads IN PLACE and 1.2 reads the same
  # staged files in backed mode — running them in parallel would nondeterminize
  # the CombinedPBMC dataset. Submit 1.1 first, then gate 1.2 behind it with
  # --dependency=afterok (fail-closed: a failed cap means 1.2 is never
  # submitted and the sacct gate below reports non-COMPLETED).
  CAP_STEP_SCRIPT="${SCRIPT_DIR}/1.1_submit_gongsharma.sh"
  if [[ -f "${CAP_STEP_SCRIPT}" ]]; then
    echo "=== Submitting $(basename "${CAP_STEP_SCRIPT}") (prerequisite for CombinedPBMC) ==="
    STEP_LOG_STEM="${LOGS_DIR}/$(basename "${CAP_STEP_SCRIPT}" .sh)_%j"
    CAP_JOB_ID="$(sbatch --parsable \
        --partition="${SLURM_PARTITION}" \
        --output="${STEP_LOG_STEM}.log" \
        --error="${STEP_LOG_STEM}.err" \
        --mail-user="${USER_EMAIL}" \
        "${CAP_STEP_SCRIPT}")"
    JOB_IDS+=("${CAP_JOB_ID}")
  fi
  for step_script in "${STEP_SCRIPTS[@]}"; do
    step_name="$(basename "${step_script}")"
    # Skip the cap step in the loop: the glob above (1.*_submit_*.sh) matches
    # 1.1_submit_gongsharma.sh, which was ALREADY submitted as the cap
    # prerequisite (two concurrent cap jobs would race on the deterministic
    # *.capped_tmp.h5ad temp path + os.replace, and CombinedPBMC would be
    # gated only on the first).
    if [[ "${step_name}" == "1.1_submit_gongsharma.sh" ]]; then
      echo "=== Skipping ${step_name} (already submitted as cap prerequisite above) ==="
      continue
    fi
    echo "=== Submitting ${step_name} ==="
    # --output/--error/--partition passed on the sbatch command line (not
    # #SBATCH lines): SLURM directives do not expand environment variables.
    STEP_LOG_STEM="${LOGS_DIR}/$(basename "${step_script}" .sh)_%j"
    SBATCH_ARGS=(
      --partition="${SLURM_PARTITION}"
      --output="${STEP_LOG_STEM}.log"
      --error="${STEP_LOG_STEM}.err"
      --mail-user="${USER_EMAIL}"
    )
    if [[ "${step_name}" == "1.2_submit_combinedpbmc.sh" && -n "${CAP_JOB_ID}" ]]; then
      SBATCH_ARGS+=(--dependency="afterok:${CAP_JOB_ID}")
    fi
    JOB_IDS+=("$(sbatch --parsable "${SBATCH_ARGS[@]}" "${step_script}")")
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

# Per-step report for the completion email: one line per job (bare job id row
# only — excludes .batch/.extern step rows), JobName + State + Elapsed +
# ExitCode. Plain jobs, no task->dataset mapping needed.
stage2_report() {
  local job_id jid jname state elapsed exitcode
  printf '%s\n' "Per-step report (step, state, elapsed, exit code):"
  for job_id in "${JOB_IDS[@]}"; do
    while IFS='|' read -r jid jname state elapsed exitcode; do
      jid="${jid//[[:space:]]/}"
      exitcode="${exitcode%|}"
      [[ "${jid}" == "${job_id}" ]] || continue
      printf '  %-40s %-14s %s (%s)\n' "${jname}" "${state}" "${elapsed}" "${exitcode}"
    done < <(sacct -j "${job_id}" -n -P --format=JobID,JobName,State,Elapsed,ExitCode 2>/dev/null || true)
  done
}

# Exit non-zero if any step failed
failed=0
FAILED_STEP_NAMES=()
for job_id in "${JOB_IDS[@]}"; do
  state=$(sacct -j "${job_id}" -X -n -o "State" 2>/dev/null | tr -d ' \n')
  if [[ "${state}" != "COMPLETED" ]]; then
    echo "ERROR: Job ${job_id} ended in state '${state}'."
    failed=1
    FAILED_STEP_NAME="$(sacct -j "${job_id}" -X -n -P --format=JobName 2>/dev/null | head -1 | tr -d '[:space:]|' || true)"
    [[ -n "${FAILED_STEP_NAME}" ]] || FAILED_STEP_NAME="${job_id}"
    FAILED_STEP_NAMES+=("${FAILED_STEP_NAME}")
  fi
done

if [[ ${failed} -ne 0 ]]; then
  notify_sync_status \
    "ECODA: stage-2 steps FAILED (${#FAILED_STEP_NAMES[@]}/${#JOB_IDS[@]})" \
    "Stage-2 dataset-specific preprocessing steps FAILED (${#FAILED_STEP_NAMES[@]}/${#JOB_IDS[@]}). Failed steps: ${FAILED_STEP_NAMES[*]}.
$(stage2_report)"
  exit 1
fi

notify_sync_status \
  "ECODA: stage-2 steps COMPLETED (${#JOB_IDS[@]} steps)" \
  "Stage-2 dataset-specific preprocessing completed (${#JOB_IDS[@]} steps; no NAS sync by design).
$(stage2_report)"

echo "All dataset-specific preprocessing steps completed."
