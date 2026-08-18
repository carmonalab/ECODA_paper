#!/bin/bash
# ==============================================================================
# run_subset_worker.sh -- Worker task for HPC dataset subsetting
# ==============================================================================
# Invoked by run_subset_hpc.sh as a SLURM array task or direct srun.
# Processes one dataset key from the catalog using ${PYTHON_BIN}.
# ==============================================================================
set -euo pipefail

# Spool-safe SCRIPT_DIR recovery (AGENTS.md HPC invariant)
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/src/slurm_config.sh"
cd "${PROJECT_ROOT}"

KEYS=(alzheimer breast covid19 diabetes kidney lung lupus myocardial parkinson)

IN_DIR="${1:-${HPC_SCRATCH_DIR}/_downloads}"
OUT_DIR="${2:-${IN_DIR}/subsets}"
TASK_INDEX="${SLURM_ARRAY_TASK_ID:-0}"

if [[ -n "${3:-}" ]]; then
    KEY="$3"
else
    KEY="${KEYS[${TASK_INDEX}]}"
fi

echo "=============================================================================="
echo "Worker starting for dataset key: ${KEY} (Task ID: ${TASK_INDEX})"
echo "Host: $(hostname) | Date: $(date) | CPUs: $(nproc)"
echo "Input dir:  ${IN_DIR}"
echo "Output dir: ${OUT_DIR}"
echo "=============================================================================="

mkdir -p "${OUT_DIR}"

PYTHON_EXEC="${PYTHON_BIN}"
if [[ ! -x "${PYTHON_EXEC}" ]]; then
    PYTHON_EXEC="$(which python3)"
fi

START_SEC=$(date +%s)
"${PYTHON_EXEC}" "${PROJECT_ROOT}/notebooks/dataset_onboarding/create_subsets_hpc.py" \
    --in-dir "${IN_DIR}" \
    --out-dir "${OUT_DIR}" \
    --only "${KEY}"
END_SEC=$(date +%s)
ELAPSED=$((END_SEC - START_SEC))

echo "=============================================================================="
echo "Worker finished for ${KEY} in ${ELAPSED}s with exit code 0"
echo "=============================================================================="
