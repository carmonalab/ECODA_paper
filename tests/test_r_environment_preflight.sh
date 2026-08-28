#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="" SLURM_JOB_ID=999999 SLURM_SUBMIT_DIR="${ROOT}" \
  R_ENV_PREFLIGHT_RSCRIPT=true \
  bash "${ROOT}/src/utils/bash/r_environment_preflight_worker.sh"
echo "R environment preflight bootstrap: OK"
