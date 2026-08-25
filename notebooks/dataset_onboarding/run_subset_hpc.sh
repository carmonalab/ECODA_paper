#!/bin/bash
# ==============================================================================
# run_subset_hpc.sh -- Submit parallel full-file registry audits + diagnostic subsets
# ==============================================================================
# Submits a SLURM array job to worker nodes so all 9 datasets run in parallel
# (~2-3 minutes total) on BeeGFS scratch instead of sequentially on login nodes.
#
# Usage (from HPC login node, repo root):
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh                  # all 9 datasets in parallel
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh --only breast    # one dataset
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh --debug-cpu     # submit to debug-cpu (max 15 min)
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh --direct        # direct execution (debug only)
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh --sync-nas      # copy subsets to NAS too
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/src/slurm_config.sh"
cd "${PROJECT_ROOT}"

KEYS=(alzheimer breast covid19 diabetes kidney lung lupus myocardial parkinson)
IN_DIR="${HPC_SCRATCH_DIR}/_downloads"
OUT_DIR="${IN_DIR}/subsets"
NAS_SUBSETS="${NAS_SC_DIR}/JooM_2025_41097818/subsets"
SUBSET_LOGS_DIR="${LOGS_DIR}/onboard_subsets"

ONLY=""
PARTITION="shared-cpu"
TIME_LIMIT="00:45:00"
DIRECT=0
SYNC_NAS=0
SKIP_EXISTING=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)
      ONLY="${2:-}"
      shift 2
      ;;
    --partition)
      PARTITION="${2:-}"
      shift 2
      ;;
    --time)
      TIME_LIMIT="${2:-}"
      shift 2
      ;;
    --debug-cpu)
      PARTITION="debug-cpu"
      TIME_LIMIT="00:15:00"
      shift
      ;;
    --skip-existing)
      SKIP_EXISTING=1
      shift
      ;;
    --direct)
      DIRECT=1
      shift
      ;;
    --sync-nas)
      SYNC_NAS=1
      shift
      ;;
    *)
      echo "ERROR: Unknown option $1 (accepted: --only <key> | --partition <p> | --time <t> | --skip-existing | --debug-cpu | --direct | --sync-nas)" >&2
      exit 1
      ;;
  esac
done

mkdir -p "${OUT_DIR}" "${SUBSET_LOGS_DIR}"

EXTRA_WORKER_ARG=""
if [[ "${SKIP_EXISTING}" -eq 1 ]]; then
  EXTRA_WORKER_ARG="--skip-existing"
fi

if [[ "${DIRECT}" -eq 1 ]]; then
  echo "=== Running direct subsetting in current shell ==="
  WORKER_SCRIPT="${SCRIPT_DIR}/run_subset_worker.sh"
  if [[ -n "${ONLY}" ]]; then
    "${WORKER_SCRIPT}" "${IN_DIR}" "${OUT_DIR}" "${ONLY}" "${EXTRA_WORKER_ARG}"
  else
    for i in "${!KEYS[@]}"; do
      SLURM_ARRAY_TASK_ID="${i}" "${WORKER_SCRIPT}" "${IN_DIR}" "${OUT_DIR}" "" "${EXTRA_WORKER_ARG}"
    done
  fi
  exit 0
fi

# Build array range
if [[ -n "${ONLY}" ]]; then
  TARGET_INDEX=-1
  for i in "${!KEYS[@]}"; do
    if [[ "${KEYS[$i]}" == "${ONLY}" ]]; then
      TARGET_INDEX="${i}"
      break
    fi
  done
  if [[ "${TARGET_INDEX}" -eq -1 ]]; then
    echo "ERROR: Unknown dataset key '${ONLY}'. Valid keys: ${KEYS[*]}" >&2
    exit 1
  fi
  ARRAY_SPEC="${TARGET_INDEX}"
  N_TASKS=1
else
  ARRAY_SPEC="0-$(( ${#KEYS[@]} - 1 ))"
  N_TASKS="${#KEYS[@]}"
fi

echo "=============================================================================="
echo "Submitting HPC Subsetting SLURM Array"
echo "=============================================================================="
echo "Datasets:      ${ONLY:-all 9 datasets in parallel}"
echo "Array spec:    ${ARRAY_SPEC} (${N_TASKS} tasks)"
echo "Partition:     ${PARTITION}"
echo "Time limit:    ${TIME_LIMIT}"
echo "CPUs per task: 4"
echo "Memory:        32G per task"
echo "Scratch dir:   ${IN_DIR}"
echo "Subsets dir:   ${OUT_DIR}"
echo "Logs dir:      ${SUBSET_LOGS_DIR}"
echo "Skip existing: ${SKIP_EXISTING}"
echo "=============================================================================="

# Submit sbatch array
SBATCH_CMD=(
  sbatch
  --parsable
  --job-name="onboard_subsets"
  --partition="${PARTITION}"
  --time="${TIME_LIMIT}"
  --cpus-per-task=4
  --mem=32G
  --array="${ARRAY_SPEC}"
  --output="${SUBSET_LOGS_DIR}/subset_%A_%a.out"
  --error="${SUBSET_LOGS_DIR}/subset_%A_%a.err"
  "${SCRIPT_DIR}/run_subset_worker.sh" "${IN_DIR}" "${OUT_DIR}" "" "${EXTRA_WORKER_ARG}"
)

JOB_ID="$("${SBATCH_CMD[@]}")"
echo "Submitted SLURM Array Job ID: ${JOB_ID}"
echo ""

# Non-blocking progress monitoring
echo "Monitoring task progress (polling every 5s)..."
START_TIME=$(date +%s)

while true; do
  ACTIVE=$(squeue -h -j "${JOB_ID}" 2>/dev/null | wc -l || true)
  ELAPSED=$(( $(date +%s) - START_TIME ))
  
  if [[ "${ACTIVE}" -eq 0 ]]; then
    break
  fi
  
  RUNNING=$(squeue -h -j "${JOB_ID}" -t R 2>/dev/null | wc -l || true)
  PENDING=$(squeue -h -j "${JOB_ID}" -t PD 2>/dev/null | wc -l || true)
  
  printf "\r[%3ds] Array %s: %d total active (%d running, %d pending)..." "${ELAPSED}" "${JOB_ID}" "${ACTIVE}" "${RUNNING}" "${PENDING}"
  sleep 5
done

echo ""
echo "All tasks finished after ${ELAPSED}s. Checking task statuses via sacct..."
echo ""

# Report per-task exit status
FAILED=0
while read -r line; do
  [[ -z "${line}" ]] && continue
  echo "  ${line}"
  if echo "${line}" | grep -Eq 'FAILED|CANCELLED|OUT_OF_MEMORY|TIMEOUT'; then
    FAILED=1
  fi
done < <(sacct -j "${JOB_ID}" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS --noheader 2>/dev/null | grep -v 'batch\|extern' || true)

# Check per-task logs for errors
for log in "${SUBSET_LOGS_DIR}"/subset_"${JOB_ID}"_*.err; do
  if [[ -f "${log}" && -s "${log}" ]]; then
    TASK_ID="$(basename "${log}" | sed "s/subset_${JOB_ID}_\([0-9]*\)\.err/\1/")"
    KEY_NAME="${KEYS[${TASK_ID}]:-unknown}"
    if grep -iq "error\|traceback" "${log}" 2>/dev/null; then
      echo ""
      echo "--- ERROR in task ${TASK_ID} (${KEY_NAME}) log: ${log} ---"
      tail -n 25 "${log}"
      FAILED=1
    fi
  fi
done

if [[ "${SYNC_NAS}" -eq 1 ]]; then
  echo ""
  echo "=== Syncing subsets to NAS (${NAS_SUBSETS}) ==="
  mkdir -p "${NAS_SUBSETS}"
  rsync -av "${OUT_DIR}/" "${NAS_SUBSETS}/"
fi

echo ""
echo "=============================================================================="
if [[ "${FAILED}" -eq 0 ]]; then
  echo "SUCCESS: All full-file registry audits and diagnostic subsets generated in:"
  echo "  ${OUT_DIR}"
  echo ""
  echo "To sync to your Mac, run from your local Mac repo root:"
  echo "  mkdir -p data/new_dataset_checks/subsets"
  echo "  rsync -avP bamboo:${OUT_DIR}/ data/new_dataset_checks/subsets/"
else
  echo "WARNING: Some tasks encountered errors. Review logs in ${SUBSET_LOGS_DIR}/"
fi
echo "=============================================================================="

exit "${FAILED}"
