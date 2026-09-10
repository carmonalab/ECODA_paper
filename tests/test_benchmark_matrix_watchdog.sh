#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(cd "$(mktemp -d "${TMPDIR:-/tmp}/ecoda-matrix-watchdog.XXXXXX")" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
HPC_ROOT="${TMP_DIR}/home/scratch/ECODA_paper"
SOURCE_COMMIT="bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
SNAPSHOT_ROOT="${HPC_ROOT}/_ecoda_source_snapshots/${SOURCE_COMMIT}"
SOURCE_TREE="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY_DIR="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY_DIR}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY_DIR}/source.tar"
HOST_PREFIX="${TMP_DIR}/host/.pixi/envs/py-cuda13"
RUNTIME_DIR="${HPC_ROOT}/_ecoda_runtime/stage5-test"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
RUN_ROOT="${HPC_ROOT}/_ecoda_runs/watchdog-run"
TEST_LOGS="${TMP_DIR}/ecoda-logs"
TEST_TMP="${TMP_DIR}/node-tmp"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${HOST_PREFIX}/bin" \
  "${HOST_PREFIX}/lib/R" "${SOURCE_TREE}" "${SOURCE_IDENTITY_DIR}" \
  "${RUNTIME_DIR}" "${RUN_ROOT}/manifests" "${RUN_ROOT}/status/watchdogs" \
  "${RUN_ROOT}/logs" "${TEST_LOGS}" "${TEST_TMP}" \
  "${HPC_ROOT}/batch_effect/uncorrected/embeddings" \
  "${TMP_DIR}/nas/batch_effect/uncorrected/embeddings"
cp -R "${ROOT}/src" "${SOURCE_TREE}/src"
cp "${ROOT}/datasets.json" "${SOURCE_TREE}/datasets.json"
cp "${ROOT}/config_helper.R" "${SOURCE_TREE}/config_helper.R"
cp "${ROOT}/pixi.toml" "${SOURCE_TREE}/pixi.toml"
cp "${ROOT}/pixi.lock" "${SOURCE_TREE}/pixi.lock"
mkdir -p "${SOURCE_TREE}/aux"
cp "${ROOT}/aux/scGateDB.rds" "${SOURCE_TREE}/aux/scGateDB.rds"
cp "${ROOT}/aux/genes.blocklist.rds" "${SOURCE_TREE}/aux/genes.blocklist.rds"
cp "${ROOT}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz" \
  "${SOURCE_TREE}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
printf '#!/bin/bash\nexit 0\n' > "${HOST_PREFIX}/bin/python"
printf '#!/bin/bash\nexit 0\n' > "${HOST_PREFIX}/bin/Rscript"
chmod +x "${HOST_PREFIX}/bin/python" "${HOST_PREFIX}/bin/Rscript"
sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_TREE}" .
SOURCE_ARCHIVE_SHA="$(sha256_file "${SOURCE_ARCHIVE}")"
SOURCE_CONFIG_SHA="$(sha256_file "${SOURCE_TREE}/config_helper.R")"
SOURCE_DATASETS_SHA="$(sha256_file "${SOURCE_TREE}/datasets.json")"
SOURCE_TOML_SHA="$(sha256_file "${SOURCE_TREE}/pixi.toml")"
SOURCE_LOCK_SHA="$(sha256_file "${SOURCE_TREE}/pixi.lock")"
printf 'FORMAT=1\nSOURCE_ROOT=%s\nSOURCE_COMMIT=%s\nSOURCE_ARCHIVE_PATH=%s\nSOURCE_ARCHIVE_SHA256=%s\nCONFIG_HELPER_SHA256=%s\nDATASETS_SHA256=%s\nPIXI_TOML_SHA256=%s\nPIXI_LOCK_SHA256=%s\nAUX_ROOT=%s\nSCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4\n' \
  "${SOURCE_TREE}" "${SOURCE_COMMIT}" "${SOURCE_ARCHIVE}" \
  "${SOURCE_ARCHIVE_SHA}" "${SOURCE_CONFIG_SHA}" "${SOURCE_DATASETS_SHA}" \
  "${SOURCE_TOML_SHA}" "${SOURCE_LOCK_SHA}" "${SOURCE_TREE}/aux" \
  > "${SOURCE_MANIFEST}"
chmod -R a-w "${SOURCE_TREE}"
touch "${SNAPSHOT_ROOT}/COMPLETE"
printf 'format-2-test-image\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA="$(sha256_file "${RUNTIME_IMAGE}")"
printf 'FORMAT=2\nIMAGE_BUILD_GIT_REVISION=build-commit-b\nIMAGE_PATH=%s\nIMAGE_SHA256=%s\nRUNTIME_ENV=py-cuda13\nRUNTIME_LAYOUT=relocated\nCONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13\nBASE_IMAGE=rockylinux:9\nPIXITAINER_VERSION=0.8.3\nPIXI_VERSION=0.49.0\nAPPTAINER_VERSION=1.3.2\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
  "${RUNTIME_IMAGE}" "${RUNTIME_IMAGE_SHA}" "${SOURCE_TOML_SHA}" \
  "${SOURCE_LOCK_SHA}" > "${RUNTIME_MANIFEST}"
chmod 444 "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}"
chmod 555 "${RUNTIME_DIR}"
cp "${SOURCE_MANIFEST}" "${RUN_ROOT}/manifests/source.manifest"
RUNTIME_MANIFEST_SHA="$(sha256_file "${RUNTIME_MANIFEST}")"
RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
printf 'RUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA}" \
  "${RUNTIME_MANIFEST_SHA}" "${RUNTIME_IMAGE_SIZE}" \
  "${RUNTIME_MANIFEST_SIZE}" "${SOURCE_TOML_SHA}" "${SOURCE_LOCK_SHA}" \
  > "${RUN_ROOT}/manifests/runtime.identity"
chmod 600 "${RUN_ROOT}/manifests/source.manifest" "${RUN_ROOT}/manifests/runtime.identity"
cat > "${TMP_DIR}/bin/sacct" <<'STUB'
#!/bin/bash
case "$*" in
  *"-j 1001"*) printf '1001|COMPLETED|0:0\n1001_1|COMPLETED|0:0\n1001_2|OUT_OF_MEMORY|0:0\n' ;;
  *"-j 1002"*) printf '1002|COMPLETED|0:0\n1002_1|COMPLETED|0:0\n' ;;
  *) printf 'COMPLETED\n' ;;
esac
STUB
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf 'ENV ECODA_SOURCE_ROOT=%s ECODA_SOURCE_MANIFEST=%s ECODA_RUNTIME_IMAGE=%s ECODA_RUNTIME_MANIFEST=%s ECODA_RUN_ROOT=%s ECODA_RUN_ID=%s ARGS %s\n' \
  "${ECODA_SOURCE_ROOT:-}" "${ECODA_SOURCE_MANIFEST:-}" \
  "${ECODA_RUNTIME_IMAGE:-}" "${ECODA_RUNTIME_MANIFEST:-}" \
  "${ECODA_RUN_ROOT:-}" "${ECODA_RUN_ID:-}" "$*" >> "${CAPTURE}"
printf '1002\n'
STUB
chmod +x "${TMP_DIR}/bin/sacct" "${TMP_DIR}/bin/sbatch"
export CAPTURE="${TMP_DIR}/sbatch.calls"
export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_SCRATCH_ROOT="${HPC_ROOT}"
export NAS_TARGET_DIR="${TMP_DIR}/nas"
export ECODA_LOGS_DIR="${TEST_LOGS}" TMPDIR="${TEST_TMP}"
export ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_AUX_ROOT="${SOURCE_TREE}/aux"
export ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
export ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" ECODA_RUNTIME_PROFILE=stage5
export ECODA_APPTAINER_NV=1 APPTAINER_BIN="${TMP_DIR}/bin/apptainer"
cat > "${TMP_DIR}/bin/apptainer" <<'STUB'
#!/bin/bash
exit 0
STUB
chmod +x "${TMP_DIR}/bin/apptainer"
export ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID=watchdog-run
export ECODA_LOGS_DIR="${RUN_ROOT}/logs"
RUNTIME_EXPORT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_SCRATCH_ROOT="${HPC_ROOT}" \
  ECODA_LOGS_DIR="${RUN_ROOT}/logs" TMPDIR="${TEST_TMP}" \
  ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
  ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" \
  ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" \
  ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" ECODA_RUNTIME_PROFILE=stage5 \
  ECODA_APPTAINER_NV=1 APPTAINER_BIN="${TMP_DIR}/bin/apptainer" \
  ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID=watchdog-run \
  bash -c '
    set -e
    source "$1/src/slurm_config.sh"
    source "$1/src/utils/bash/ecoda_runtime.sh"
    ecoda_runtime_validate_bound_run >/dev/null
    ecoda_runtime_export_csv stage5 1
  ' _ "${ROOT}"
)"
MANIFEST="${RUN_ROOT}/manifests/mrvi.tsv"
printf 'Joanito\tbatch_effect_uncorrected\tmrvi\thvg2000\nStephenson\tbatch_effect_uncorrected\tmrvi\thvg1000\n' > "${MANIFEST}"
WATCHDOG_OUTPUT="$(HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" ANALYSIS_PASS=uncorrected \
  PROJECT_ROOT="${SOURCE_TREE}" \
  ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" \
  ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_PROFILE=stage5 ECODA_APPTAINER_NV=1 \
  METHOD_TIME_LIMIT=03:00:00 MATRIX_WATCHDOG_MAX_POLLS=1 \
  SLURM_JOB_ID=999999 SLURM_SUBMIT_DIR="${ROOT}" \
  bash "${ROOT}/src/5_run_benchmark_methods/matrix_watchdog.sh" "${RUN_ROOT}" mrvi "${MANIFEST}" 1001 128G 256G shared-gpu 4 \
  "${SOURCE_TREE}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh" \
  "${RUNTIME_EXPORT}" --gpus=1)"
STATUS="${RUN_ROOT}/status/watchdogs/mrvi.status"
[[ "$(grep '^STATE=' "${STATUS}")" == "STATE=OK" ]]
RETRY="${RUN_ROOT}/manifests/mrvi.retry_1.tsv"
[[ "$(cat "${RETRY}")" == $'Stephenson\tbatch_effect_uncorrected\tmrvi\thvg1000' ]]
[[ "$(grep -c '^SCHEDULER_ID=' "${STATUS}")" == 2 ]]
[[ "$(grep -c '^SCHEDULER_ID=1001$' "${STATUS}")" == 1 ]]
[[ "$(grep -c '^SCHEDULER_ID=1002$' "${STATUS}")" == 1 ]]
case "${WATCHDOG_OUTPUT}" in *"BATCH_EFFECT_RETRY_ARRAY_JOB_ID=1002"*) ;; *) echo "batch retry marker missing" >&2; exit 1 ;; esac
if grep -q 'BENCHMARK_MANIFEST=' "${CAPTURE}"; then
  echo "batch retry exported BENCHMARK_MANIFEST" >&2
  exit 1
fi
case "$(cat "${CAPTURE}")" in *"MATRIX_RETRY_MANIFEST=${RETRY}"*) ;; *) echo "matrix retry manifest export missing" >&2; exit 1 ;; esac
RETRY_CALLS="$(cat "${CAPTURE}")"
case "${RETRY_CALLS}" in
  *"${SOURCE_TREE}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"*) ;;
  *) echo "matrix retry worker escaped the immutable source root" >&2; exit 1 ;;
esac
for required_identity in \
  "ECODA_SOURCE_ROOT=${SOURCE_TREE}" \
  "ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}" \
  "ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}" \
  "ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}" \
  "ECODA_RUN_ROOT=${RUN_ROOT}" \
  "ECODA_RUN_ID=watchdog-run"; do
  case "${RETRY_CALLS}" in
    *"${required_identity}"*) ;;
    *) echo "matrix retry omitted immutable run identity: ${required_identity}" >&2; exit 1 ;;
  esac
done
OUTSIDE_SCRIPT="${TMP_DIR}/outside-matrix-worker.sh"
printf '#!/bin/bash\nexit 0\n' > "${OUTSIDE_SCRIPT}"
chmod +x "${OUTSIDE_SCRIPT}"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" \
  bash -c '
    set -e
    source "$1/src/slurm_config.sh"
    source "$1/src/utils/bash/ecoda_run_common.sh"
    ecoda_require_source_script_path "$2" "$3"
  ' _ "${ROOT}" "${OUTSIDE_SCRIPT}" "${SOURCE_TREE}"; then
  echo "mutable retry worker script escaped source containment" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]

OWNERSHIP_SELECTION="${RUN_ROOT}/manifests/retry-ownership.tsv"
printf 'Joanito\tbatch_effect_uncorrected\tmrvi\n' > "${OWNERSHIP_SELECTION}"
RETRY_ARTIFACT="${HPC_ROOT}/batch_effect/uncorrected/embeddings/Joanito_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"
mkdir -p "${TMP_DIR}/nas/batch_effect/uncorrected/embeddings"
: > "${CAPTURE}"
(
  set -e
  export HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}"
  export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
  export ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID=watchdog-run
  export NAS_TARGET_DIR="${TMP_DIR}/nas"
  source "${ROOT}/src/slurm_config.sh"
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  ecoda_artifact_owner_acquire "${RETRY_ARTIFACT}" stage5 watchdog-run 1 0 0 >/dev/null
  ecoda_validate_output_ownership stage5 "${OWNERSHIP_SELECTION}" watchdog-run
)
[[ ! -s "${CAPTURE}" ]]
rm -rf "${HPC_ROOT}/_ecoda_owners/artifact"
: > "${CAPTURE}"
(
  set -e
  export HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}"
  export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
  export ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID=watchdog-run
  export NAS_TARGET_DIR="${TMP_DIR}/nas"
  source "${ROOT}/src/slurm_config.sh"
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  ecoda_artifact_owner_acquire "${RETRY_ARTIFACT}" stage5 other-run 1 0 0 >/dev/null
  if ecoda_validate_output_ownership stage5 "${OWNERSHIP_SELECTION}" watchdog-run; then
    exit 1
  fi
)
[[ ! -s "${CAPTURE}" ]]
rm -rf "${HPC_ROOT}/_ecoda_owners/artifact"
case "${RETRY_CALLS}" in *"--time=03:00:00"*) ;; *) echo "matrix retry worker time limit missing" >&2; exit 1 ;; esac
case "${RETRY_CALLS}" in *"ECODA_RUNTIME_PROFILE=stage5"*"ECODA_APPTAINER_NV=1"*) ;; *) echo "matrix retry runtime export missing" >&2; exit 1 ;; esac
SCHEDULER_MANIFEST="${RUN_ROOT}/manifests/scheduler_ids.tsv"
printf 'ARRAY\t1001\nWATCHDOG\t1003\n' > "${SCHEDULER_MANIFEST}"
PROJECT_ROOT="${SOURCE_TREE}" SLURM_JOB_ID=1004 HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID=watchdog-run \
  bash "${SOURCE_TREE}/src/5_run_benchmark_methods/matrix_gate.sh" \
  "${RUN_ROOT}" mrvi "${SCHEDULER_MANIFEST}"
[[ "$(grep '^STATE=' "${RUN_ROOT}/status/aggregate")" == "STATE=OK" ]]
[[ "$(grep -c '^SCHEDULER_ID=' "${RUN_ROOT}/status/aggregate")" == 4 ]]
[[ "$(grep -c '^SCHEDULER_ID=1002$' "${RUN_ROOT}/status/aggregate")" == 1 ]]
BAD_GATE_ROOT="${HPC_ROOT}/_ecoda_runs/unpinned-gate"
mkdir -p "${BAD_GATE_ROOT}/manifests" "${BAD_GATE_ROOT}/status/watchdogs"
cp "${RUN_ROOT}/manifests/source.manifest" "${BAD_GATE_ROOT}/manifests/source.manifest"
cp "${RUN_ROOT}/manifests/runtime.identity" "${BAD_GATE_ROOT}/manifests/runtime.identity"
cp "${RUN_ROOT}/status/watchdogs/mrvi.status" "${BAD_GATE_ROOT}/status/watchdogs/mrvi.status"
BAD_SCHEDULER_MANIFEST="${BAD_GATE_ROOT}/manifests/scheduler_ids.tsv"
printf 'ARRAY\t1001\nWATCHDOG\t1003\n' > "${BAD_SCHEDULER_MANIFEST}"
if PROJECT_ROOT="${ROOT}" SLURM_JOB_ID=1005 HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  ECODA_SOURCE_SNAPSHOT_REQUIRED=0 \
  ECODA_RUN_ROOT="${BAD_GATE_ROOT}" ECODA_RUN_ID=unpinned-gate \
  bash "${ROOT}/src/5_run_benchmark_methods/matrix_gate.sh" \
  "${BAD_GATE_ROOT}" mrvi "${BAD_SCHEDULER_MANIFEST}"; then
  echo "mutable/unpinned matrix gate was accepted" >&2
  exit 1
fi
[[ ! -e "${BAD_GATE_ROOT}/status/aggregate" ]]
