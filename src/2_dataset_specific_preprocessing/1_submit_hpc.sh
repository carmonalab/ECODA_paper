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
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

STEP_SCRIPTS=( "${SCRIPT_DIR}"/1.*_submit_*.sh )

if [[ ${#STEP_SCRIPTS[@]} -eq 0 ]]; then
  echo "ERROR: No step scripts found in ${SCRIPT_DIR}."
  exit 1
fi

JOB_IDS=()
for step_script in "${STEP_SCRIPTS[@]}"; do
  echo "=== Submitting $(basename "${step_script}") ==="
  JOB_IDS+=("$(sbatch --parsable "${step_script}")")
done

echo "Submitted jobs: ${JOB_IDS[*]}"

# Wait for all jobs to leave the queue
while :; do
  pending=0
  for job_id in "${JOB_IDS[@]}"; do
    if squeue -u "$USER" -j "${job_id}" -h -o "%i" 2>/dev/null | grep -q "${job_id}"; then
      pending=1
    fi
  done
  if [[ ${pending} -eq 0 ]]; then
    break
  fi
  sleep 60
done

# Give sacct a moment to record final states
sleep 30

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
