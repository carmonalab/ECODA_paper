#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${ROOT}/src/slurm_config.sh"
EXPECTED_RSCRIPT="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/Rscript --vanilla"
[[ "${PIXI_RSCRIPT}" == "${EXPECTED_RSCRIPT}" ]] || {
  echo "ERROR: workers must use the direct py-cuda13 Rscript binary." >&2
  exit 1
}

PROJECT_ROOT="" SLURM_JOB_ID=999999 SLURM_SUBMIT_DIR="${ROOT}" \
  R_ENV_PREFLIGHT_RSCRIPT=true \
  bash "${ROOT}/src/utils/bash/r_environment_preflight_worker.sh"
echo "R environment preflight bootstrap: OK"
