#!/bin/bash
# ==============================================================================
# run_subset_hpc.sh -- Run dataset subsetting directly on HPC scratch
# ==============================================================================
# Executes create_subsets_hpc.py using BeeGFS scratch files to generate
# lightweight diagnostic subsets (~15-40 MB each, ~200 MB total) and metadata
# summaries in seconds.
#
# Usage (from HPC login node or compute session; repo root):
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh                  # all 9 datasets
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh --only breast    # one dataset
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh --debug-cpu     # via srun on debug-cpu
#   ./notebooks/dataset_onboarding/run_subset_hpc.sh --sync-nas      # copy subsets to NAS too
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_ROOT}/src/slurm_config.sh"
cd "${PROJECT_ROOT}"

IN_DIR="${HPC_SCRATCH_DIR}/_downloads"
OUT_DIR="${IN_DIR}/subsets"
NAS_SUBSETS="${NAS_SC_DIR}/JooM_2025_41097818/subsets"

ONLY=""
DEBUG_CPU=0
SYNC_NAS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)
      ONLY="${2:-}"
      shift 2
      ;;
    --debug-cpu)
      DEBUG_CPU=1
      shift
      ;;
    --sync-nas)
      SYNC_NAS=1
      shift
      ;;
    *)
      echo "ERROR: Unknown option $1 (accepted: --only <key> | --debug-cpu | --sync-nas)" >&2
      exit 1
      ;;
  esac
done

mkdir -p "${OUT_DIR}"

CMD_ARGS=("notebooks/dataset_onboarding/create_subsets_hpc.py" "--in-dir" "${IN_DIR}" "--out-dir" "${OUT_DIR}")
if [[ -n "${ONLY}" ]]; then
  CMD_ARGS+=("--only" "${ONLY}")
fi

echo "=== Running HPC dataset subsetting ==="
echo "Input directory:  ${IN_DIR}"
echo "Output directory: ${OUT_DIR}"

if [[ "${DEBUG_CPU}" -eq 1 ]]; then
  echo "Executing via srun on partition '${SLURM_PARTITION_CPU_DEBUG}'..."
  srun -p "${SLURM_PARTITION_CPU_DEBUG}" --time=00:30:00 --mem=32G --cpus-per-task=4 \
    pixi run --as-is -e default python "${CMD_ARGS[@]}"
else
  pixi run --as-is -e default python "${CMD_ARGS[@]}"
fi

if [[ "${SYNC_NAS}" -eq 1 ]]; then
  echo ""
  echo "=== Syncing subsets to NAS (${NAS_SUBSETS}) ==="
  mkdir -p "${NAS_SUBSETS}"
  rsync -av "${OUT_DIR}/" "${NAS_SUBSETS}/"
fi

echo ""
echo "=============================================================================="
echo "Subsets generated successfully in: ${OUT_DIR}"
echo ""
echo "To pull the subsets directly to your local Mac, run this from your Mac repo root:"
echo "  mkdir -p data/new_dataset_checks/subsets"
echo "  rsync -avP bamboo:${OUT_DIR}/ data/new_dataset_checks/subsets/"
echo "=============================================================================="
