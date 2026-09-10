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
RUNTIME_TEST_ROOT="$(cd "${RUNTIME_TEST_ROOT}" && pwd -P)"
cleanup_runtime_test() {
  chmod -R u+w "${RUNTIME_TEST_ROOT}" >/dev/null 2>&1 || true
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

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

# Keep the final worker smoke independent of the mutable checkout.  The
# fixture is intentionally small, but contains every file and identity field
# required by the FORMAT=1 source snapshot validator.
SNAPSHOT_COMMIT="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SNAPSHOT_ROOT="${RUNTIME_TEST_ROOT}/snapshots/${SNAPSHOT_COMMIT}"
SOURCE_ROOT="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY}/source.tar"
mkdir -p "${SOURCE_ROOT}/src/utils/bash" "${SOURCE_ROOT}/aux" \
  "${SOURCE_IDENTITY}"
cp "${ROOT}/src/slurm_config.sh" "${SOURCE_ROOT}/src/slurm_config.sh"
for source_file in ecoda_run_common.sh ecoda_runtime.sh \
  r_environment_preflight_worker.sh; do
  cp "${ROOT}/src/utils/bash/${source_file}" \
    "${SOURCE_ROOT}/src/utils/bash/${source_file}"
done
printf '{}\n' > "${SOURCE_ROOT}/datasets.json"
printf 'fixture config\n' > "${SOURCE_ROOT}/config_helper.R"
printf 'fixture pixi toml\n' > "${SOURCE_ROOT}/pixi.toml"
printf 'fixture pixi lock\n' > "${SOURCE_ROOT}/pixi.lock"
printf 'fixture scGate database\n' > "${SOURCE_ROOT}/aux/scGateDB.rds"
printf 'fixture genes blocklist\n' > "${SOURCE_ROOT}/aux/genes.blocklist.rds"
printf 'fixture Ensembl map\n' \
  > "${SOURCE_ROOT}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_ROOT}" .
SOURCE_ARCHIVE_SHA256="$(sha256_file "${SOURCE_ARCHIVE}")"
SOURCE_CONFIG_SHA256="$(sha256_file "${SOURCE_ROOT}/config_helper.R")"
SOURCE_DATASETS_SHA256="$(sha256_file "${SOURCE_ROOT}/datasets.json")"
SOURCE_TOML_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.toml")"
SOURCE_LOCK_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.lock")"
cat > "${SOURCE_MANIFEST}" <<EOF
FORMAT=1
SOURCE_ROOT=${SOURCE_ROOT}
SOURCE_COMMIT=${SNAPSHOT_COMMIT}
SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}
SOURCE_ARCHIVE_SHA256=${SOURCE_ARCHIVE_SHA256}
CONFIG_HELPER_SHA256=${SOURCE_CONFIG_SHA256}
DATASETS_SHA256=${SOURCE_DATASETS_SHA256}
PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}
PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}
AUX_ROOT=${SOURCE_ROOT}/aux
SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4
EOF
printf 'COMPLETE\n' > "${SNAPSHOT_ROOT}/COMPLETE"
chmod -R a-w "${SNAPSHOT_ROOT}"

# The host environment has the same direct interpreter layout as the real
# py-cuda13 prefix.  Its Rscript stub asserts the bindings that the worker
# must retain after validating the immutable source/runtime pair.
HOST_ENV_PREFIX="${RUNTIME_TEST_ROOT}/host/.pixi/envs/py-cuda13"
SCRATCH_ROOT="${RUNTIME_TEST_ROOT}/scratch/ECODA_paper"
LOGS_ROOT="${RUNTIME_TEST_ROOT}/logs"
TMP_ROOT="${RUNTIME_TEST_ROOT}/tmp"
mkdir -p "${HOST_ENV_PREFIX}/bin" "${HOST_ENV_PREFIX}/lib" \
  "${SCRATCH_ROOT}" "${LOGS_ROOT}" "${TMP_ROOT}"
cat > "${HOST_ENV_PREFIX}/bin/python" <<'STUB'
#!/bin/bash
exit 0
STUB
cat > "${HOST_ENV_PREFIX}/bin/Rscript" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ "${1:-}" == "--vanilla" && "${2:-}" == "-e" ]] || exit 10
[[ "${PROJECT_ROOT:-}" == "${R_PREFLIGHT_EXPECTED_SOURCE_ROOT}" ]] || exit 11
[[ "${ECODA_SOURCE_ROOT:-}" == "${R_PREFLIGHT_EXPECTED_SOURCE_ROOT}" ]] || exit 12
[[ "${ECODA_SOURCE_MANIFEST:-}" == "${R_PREFLIGHT_EXPECTED_SOURCE_MANIFEST}" ]] || exit 13
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-}" == "1" ]] || exit 14
[[ "${ECODA_AUX_ROOT:-}" == "${R_PREFLIGHT_EXPECTED_SOURCE_ROOT}/aux" ]] || exit 15
[[ -r "${ECODA_AUX_ROOT}/scGateDB.rds" ]] || exit 16
[[ "${SCGATE_DB_PATH:-}" == "${R_PREFLIGHT_EXPECTED_SOURCE_ROOT}/aux/scGateDB.rds" ]] || exit 17
[[ "${ECODA_HOST_ENV_PREFIX:-}" == "${R_PREFLIGHT_EXPECTED_HOST_PREFIX}" ]] || exit 18
[[ "${ECODA_HOST_PYTHON_BIN:-}" == "${R_PREFLIGHT_EXPECTED_HOST_PREFIX}/bin/python" ]] || exit 19
[[ "${ECODA_HOST_PIXI_RSCRIPT:-}" == "${R_PREFLIGHT_EXPECTED_HOST_PREFIX}/bin/Rscript --vanilla" ]] || exit 20
[[ "${PYTHON_BIN:-}" == "${R_PREFLIGHT_EXPECTED_HOST_PREFIX}/bin/python" ]] || exit 21
[[ "${PIXI_RSCRIPT:-}" == "${R_PREFLIGHT_EXPECTED_HOST_PREFIX}/bin/Rscript --vanilla" ]] || exit 22
[[ "${ECODA_RUN_ROOT:-}" == "${R_PREFLIGHT_EXPECTED_RUN_ROOT}" ]] || exit 23
[[ "${ECODA_RUN_ID:-}" == "${R_PREFLIGHT_EXPECTED_RUN_ID}" ]] || exit 24
[[ "${R_ENV_PREFLIGHT_RUN_ROOT:-}" == "${R_PREFLIGHT_EXPECTED_RUN_ROOT}" ]] || exit 25
[[ "${R_ENV_PREFLIGHT_RUN_ID:-}" == "${R_PREFLIGHT_EXPECTED_RUN_ID}" ]] || exit 26
[[ "${ECODA_RUNTIME_IMAGE:-}" == "${R_PREFLIGHT_EXPECTED_RUNTIME_IMAGE}" ]] || exit 27
[[ "${ECODA_RUNTIME_MANIFEST:-}" == "${R_PREFLIGHT_EXPECTED_RUNTIME_MANIFEST}" ]] || exit 28
printf 'host Rscript smoke\n' > "${R_PREFLIGHT_MARKER:?}"
STUB
chmod +x "${HOST_ENV_PREFIX}/bin/python" "${HOST_ENV_PREFIX}/bin/Rscript"
HOST_PYTHON_SHA256="$(sha256_file "${HOST_ENV_PREFIX}/bin/python")"
HOST_RSCRIPT_SHA256="$(sha256_file "${HOST_ENV_PREFIX}/bin/Rscript")"

RUNTIME_ID="r-preflight-runtime"
RUNTIME_DIR="${SCRATCH_ROOT}/_ecoda_runtime/${RUNTIME_ID}"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
mkdir -p "${RUNTIME_DIR}"
printf 'deterministic FORMAT=2 runtime fixture\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA256="$(sha256_file "${RUNTIME_IMAGE}")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_PATH=${RUNTIME_IMAGE}
IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.49.0
APPTAINER_VERSION=1.3.2
IMAGE_BUILD_GIT_REVISION=runtime-build-fixture
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}
EOF
RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")"
RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
chmod -R a-w "${RUNTIME_DIR}"

RUN_ID="r-environment-preflight"
RUN_ROOT="${SCRATCH_ROOT}/_ecoda_runs/${RUN_ID}"
mkdir -p "${RUN_ROOT}/manifests" "${RUN_ROOT}/logs"
cp "${SOURCE_MANIFEST}" "${RUN_ROOT}/manifests/source.manifest"
cat > "${RUN_ROOT}/manifests/runtime.identity" <<EOF
RUNTIME_IMAGE=${RUNTIME_IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
RUNTIME_IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}
RUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA256}
RUNTIME_IMAGE_SIZE=${RUNTIME_IMAGE_SIZE}
RUNTIME_MANIFEST_SIZE=${RUNTIME_MANIFEST_SIZE}
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}
EOF
chmod 444 "${RUN_ROOT}/manifests/source.manifest" \
  "${RUN_ROOT}/manifests/runtime.identity"
cat > "${RUN_ROOT}/metadata" <<EOF
STAGE=stage5
RUN_ID=${RUN_ID}
STATE=ACTIVE
SOURCE_MANIFEST=${SOURCE_MANIFEST}
SOURCE_MANIFEST_COPY=${RUN_ROOT}/manifests/source.manifest
SOURCE_ROOT=${SOURCE_ROOT}
SOURCE_COMMIT=${SNAPSHOT_COMMIT}
SOURCE_AUX_ROOT=${SOURCE_ROOT}/aux
RUNTIME_IDENTITY=${RUN_ROOT}/manifests/runtime.identity
RUNTIME_IMAGE=${RUNTIME_IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
EOF

[[ "$(wc -l < "${SOURCE_MANIFEST}" | tr -d '[:space:]')" == "11" ]] || exit 1
[[ "$(sed -n '1p' "${SOURCE_MANIFEST}")" == "FORMAT=1" ]] || exit 1
[[ "$(sed -n '1p' "${RUNTIME_MANIFEST}")" == "FORMAT=2" ]] || exit 1
[[ "$(wc -l < "${RUN_ROOT}/manifests/runtime.identity" | tr -d '[:space:]')" == "8" ]] || exit 1
case "$(cat "${RUN_ROOT}/metadata")" in
  *"SOURCE_ROOT=${SOURCE_ROOT}"*"RUNTIME_IMAGE=${RUNTIME_IMAGE}"*) ;;
  *) echo "ERROR: preflight run metadata omitted source/runtime identity." >&2; exit 1 ;;
esac

export HOME="${RUNTIME_TEST_ROOT}/home"
export TMPDIR="${TMP_ROOT}"
export HPC_SCRATCH_DIR="${SCRATCH_ROOT}"
export ECODA_SCRATCH_ROOT="${SCRATCH_ROOT}"
export ECODA_LOGS_DIR="${LOGS_ROOT}"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_AUX_ROOT="${SOURCE_ROOT}/aux"
export ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}"
export ECODA_HOST_PYTHON_BIN="${HOST_ENV_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${HOST_ENV_PREFIX}/bin/Rscript --vanilla"
export ECODA_HOST_PYTHON_SHA256="${HOST_PYTHON_SHA256}"
export ECODA_HOST_RSCRIPT_SHA256="${HOST_RSCRIPT_SHA256}"
export ECODA_RUNTIME_MODE=host
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_PROFILE=stage5
export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export ECODA_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
export ECODA_RUN_ROOT="${RUN_ROOT}"
export ECODA_RUN_ID="${RUN_ID}"
export R_ENV_PREFLIGHT_RUN_ROOT="${RUN_ROOT}"
export R_ENV_PREFLIGHT_RUN_ID="${RUN_ID}"
export R_ENV_PREFLIGHT_RSCRIPT="${HOST_ENV_PREFIX}/bin/Rscript --vanilla"
export R_PREFLIGHT_EXPECTED_SOURCE_ROOT="${SOURCE_ROOT}"
export R_PREFLIGHT_EXPECTED_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export R_PREFLIGHT_EXPECTED_HOST_PREFIX="${HOST_ENV_PREFIX}"
export R_PREFLIGHT_EXPECTED_RUN_ROOT="${RUN_ROOT}"
export R_PREFLIGHT_EXPECTED_RUN_ID="${RUN_ID}"
export R_PREFLIGHT_EXPECTED_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export R_PREFLIGHT_EXPECTED_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export R_PREFLIGHT_MARKER="${RUNTIME_TEST_ROOT}/rscript.marker"
export SLURM_JOB_ID=999999
export PROJECT_ROOT=""
export SLURM_SUBMIT_DIR="${ROOT}"

# Snapshot-backed workers fail closed when the immutable-source requirement is
# absent; this negative check keeps that boundary covered before the success
# smoke below.
if ECODA_SOURCE_SNAPSHOT_REQUIRED=0 \
  bash "${SOURCE_ROOT}/src/utils/bash/r_environment_preflight_worker.sh" \
  >/dev/null 2>&1; then
  echo "ERROR: preflight accepted an unpinned source." >&2
  exit 1
fi
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
bash "${SOURCE_ROOT}/src/utils/bash/r_environment_preflight_worker.sh"
[[ -s "${R_PREFLIGHT_MARKER}" ]] || {
  echo "ERROR: host Rscript smoke did not execute." >&2
  exit 1
}
echo "R environment preflight bootstrap: OK"
