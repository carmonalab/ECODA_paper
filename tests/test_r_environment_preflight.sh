#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export ECODA_RUNTIME_MODE=host
export ECODA_RUNTIME_IN_CONTAINER=0
source "${ROOT}/src/slurm_config.sh"
EXPECTED_PYTHON="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python"
EXPECTED_RSCRIPT="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/Rscript --vanilla"
[[ "${PYTHON_BIN}" == "${EXPECTED_PYTHON}" ]] || {
  echo "ERROR: workers must use the direct py-cuda13 Python binary." >&2
  exit 1
}
[[ "${PIXI_RSCRIPT}" == "${EXPECTED_RSCRIPT}" ]] || {
  echo "ERROR: workers must use the direct py-cuda13 Rscript binary." >&2
  exit 1
}
[[ "${PIXI_RSCRIPT}" != *"pixi run"* ]] || {
  echo "ERROR: worker R must not use pixi run." >&2
  exit 1
}

# ECODA_RUNTIME_MODE=apptainer is a submission-time choice; host-side
# validators still use the direct host binaries until a worker re-execs.
(
  export ECODA_RUNTIME_MODE=apptainer
  export ECODA_RUNTIME_IN_CONTAINER=0
  source "${ROOT}/src/slurm_config.sh"
  [[ "${PYTHON_BIN}" == "${EXPECTED_PYTHON}" ]] || {
    echo "ERROR: apptainer mode changed the host Python path before re-exec." >&2
    exit 1
  }
  [[ "${PIXI_RSCRIPT}" == "${EXPECTED_RSCRIPT}" ]] || {
    echo "ERROR: apptainer mode changed the host R path before re-exec." >&2
    exit 1
  }
)

# The container branch is configuration-only here: no SIF is needed to prove
# that direct in-image paths and host import overrides are sanitized.
RUNTIME_TEST_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-r-preflight.XXXXXX")"
cleanup_runtime_test() {
  rm -rf "${RUNTIME_TEST_ROOT}"
}
trap cleanup_runtime_test EXIT
mkdir -p "${RUNTIME_TEST_ROOT}/bin" "${RUNTIME_TEST_ROOT}/lib/R"
touch "${RUNTIME_TEST_ROOT}/bin/python" "${RUNTIME_TEST_ROOT}/bin/Rscript"
chmod +x "${RUNTIME_TEST_ROOT}/bin/python" "${RUNTIME_TEST_ROOT}/bin/Rscript"
(
  export ECODA_RUNTIME_MODE=apptainer
  export ECODA_RUNTIME_IN_CONTAINER=1
  export ECODA_RUNTIME_PREFIX="${RUNTIME_TEST_ROOT}"
  export PYTHONHOME=/tmp/host-python
  export PYTHONPATH=/tmp/host-pythonpath
  export R_LIBS_USER=/tmp/host-r-user
  export R_LIBS_SITE=/tmp/host-r-site
  source "${ROOT}/src/slurm_config.sh"
  [[ "${PYTHON_BIN}" == "${RUNTIME_TEST_ROOT}/bin/python" ]] || exit 1
  [[ "${PIXI_RSCRIPT}" == "${RUNTIME_TEST_ROOT}/bin/Rscript --vanilla" ]] || exit 1
  [[ "${R_HOME}" == "${RUNTIME_TEST_ROOT}/lib/R" ]] || exit 1
  [[ "${RETICULATE_PYTHON}" == "${RUNTIME_TEST_ROOT}/bin/python" ]] || exit 1
  [[ -z "${PYTHONHOME:-}" && -z "${PYTHONPATH:-}" ]] || exit 1
  [[ -z "${R_LIBS_USER:-}" && -z "${R_LIBS_SITE:-}" ]] || exit 1
)

PROJECT_ROOT="" SLURM_JOB_ID=999999 SLURM_SUBMIT_DIR="${ROOT}" \
  R_ENV_PREFLIGHT_RSCRIPT=true \
  bash "${ROOT}/src/utils/bash/r_environment_preflight_worker.sh"
echo "R environment preflight bootstrap: OK"
