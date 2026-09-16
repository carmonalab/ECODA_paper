#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage5-corrected-retry.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

HPC_ROOT="${TMP_DIR}/home/scratch/ECODA_paper"
NAS_ROOT="${TMP_DIR}/nas/project"
SOURCE_COMMIT="cccccccccccccccccccccccccccccccccccccccc"
SNAPSHOT_ROOT="${HPC_ROOT}/_ecoda_source_snapshots/${SOURCE_COMMIT}"
SOURCE_TREE="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY_DIR="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY_DIR}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY_DIR}/source.tar"
HOST_PREFIX="${TMP_DIR}/host/.pixi/envs/py-cuda13"
RUNTIME_DIR="${HPC_ROOT}/_ecoda_runtime/stage5-corrected-retry"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
CORRECTED_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/corrected-retry"
LEGACY_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/legacy-retry"
TEST_TMP="${TMP_DIR}/node-tmp"
CAPTURE="${TMP_DIR}/sbatch.calls"

mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${HOST_PREFIX}/bin" \
  "${HOST_PREFIX}/lib/R" "${SOURCE_TREE}" "${SOURCE_IDENTITY_DIR}" \
  "${RUNTIME_DIR}" "${CORRECTED_RUN_ROOT}/manifests" \
  "${CORRECTED_RUN_ROOT}/logs" "${LEGACY_RUN_ROOT}/manifests" \
  "${LEGACY_RUN_ROOT}/logs" "${TEST_TMP}" \
  "${HPC_ROOT}/batch_effect/corrected_final" \
  "${HPC_ROOT}/batch_effect/benchmark" \
  "${NAS_ROOT}/batch_effect/corrected_final" \
  "${NAS_ROOT}/batch_effect/benchmark"

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

printf 'corrected-retry-runtime-image\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA="$(sha256_file "${RUNTIME_IMAGE}")"
printf 'FORMAT=2\nIMAGE_BUILD_GIT_REVISION=build-c\nIMAGE_PATH=%s\nIMAGE_SHA256=%s\nRUNTIME_ENV=py-cuda13\nRUNTIME_LAYOUT=relocated\nCONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13\nBASE_IMAGE=rockylinux:9\nPIXITAINER_VERSION=0.8.3\nPIXI_VERSION=0.49.0\nAPPTAINER_VERSION=1.3.2\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
  "${RUNTIME_IMAGE}" "${RUNTIME_IMAGE_SHA}" "${SOURCE_TOML_SHA}" \
  "${SOURCE_LOCK_SHA}" > "${RUNTIME_MANIFEST}"
RUNTIME_MANIFEST_SHA="$(sha256_file "${RUNTIME_MANIFEST}")"
RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
chmod 444 "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}"
chmod 555 "${RUNTIME_DIR}"

write_run_identity() {
  local run_root="$1"
  mkdir -p "${run_root}/manifests" "${run_root}/logs"
  cp "${SOURCE_MANIFEST}" "${run_root}/manifests/source.manifest"
  printf 'RUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
    "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA}" \
    "${RUNTIME_MANIFEST_SHA}" "${RUNTIME_IMAGE_SIZE}" \
    "${RUNTIME_MANIFEST_SIZE}" "${SOURCE_TOML_SHA}" "${SOURCE_LOCK_SHA}" \
    > "${run_root}/manifests/runtime.identity"
  chmod 600 "${run_root}/manifests/source.manifest" \
    "${run_root}/manifests/runtime.identity"
}
write_run_identity "${CORRECTED_RUN_ROOT}"
write_run_identity "${LEGACY_RUN_ROOT}"

cat > "${TMP_DIR}/bin/sacct" <<'STUB'
#!/bin/bash
set -euo pipefail
job=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -j) job="${2:-}"; shift 2 ;;
    *) shift ;;
  esac
done
case "${job}" in
  7001) printf '7001_1|OUT_OF_MEMORY|0:0\n' ;;
  7002) printf '7002_1|COMPLETED|0:0\n' ;;
  7101) printf '7101_1|OUT_OF_MEMORY|0:0\n' ;;
  7102) printf '7102_1|COMPLETED|0:0\n' ;;
  *) printf '%s_1|COMPLETED|0:0\n' "${job}" ;;
esac
STUB
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
printf '%s\n' "${SBATCH_RESULT_ID}"
STUB
chmod +x "${TMP_DIR}/bin/sacct" "${TMP_DIR}/bin/sbatch"

runtime_export_for() {
  local run_root="$1"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_SCRATCH_ROOT="${HPC_ROOT}" \
    NAS_TARGET_DIR="${NAS_ROOT}" ECODA_LOGS_DIR="${run_root}/logs" \
    TMPDIR="${TEST_TMP}" ECODA_SOURCE_ROOT="${SOURCE_TREE}" \
    ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
    ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_TREE}/aux" \
    ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" ECODA_RUNTIME_MODE=apptainer \
    ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" \
    ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
    ECODA_RUNTIME_PROFILE=stage5 ECODA_APPTAINER_NV=0 \
    ECODA_RUN_ROOT="${run_root}" ECODA_RUN_ID="${run_root##*/}" \
    bash -c '
      set -e
      source "$1/src/slurm_config.sh"
      source "$1/src/utils/bash/ecoda_runtime.sh"
      ecoda_runtime_validate_bound_run >/dev/null
      ecoda_runtime_export_csv stage5 0
    ' _ "${SOURCE_TREE}"
}

export CAPTURE HOME="${TMP_DIR}/home" HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_SCRATCH_ROOT="${HPC_ROOT}"
export NAS_TARGET_DIR="${NAS_ROOT}" TMPDIR="${TEST_TMP}"
export PATH="${TMP_DIR}/bin:${PATH}"
export PROJECT_ROOT="${SOURCE_TREE}" SLURM_SUBMIT_DIR="${SOURCE_TREE}"
export ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_TREE}/aux"
export ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export ECODA_RUNTIME_PROFILE=stage5 ECODA_APPTAINER_NV=0
export ECODA_RUNTIME_IN_CONTAINER=0 USER_EMAIL="test@example.invalid"
export MATRIX_WATCHDOG_POLL_SECONDS=0 ECODA_ACCOUNTING_EMPTY_GRACE=1
export METHOD_TIME_LIMIT=03:00:00 FORCE_BENCHMARK=0
export WORKER_SCRIPT="${SOURCE_TREE}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"

CORRECTED_MANIFEST="${CORRECTED_RUN_ROOT}/manifests/composition.tsv"
printf 'Joanito\tbatch_effect_corrected\tcomposition\n' > "${CORRECTED_MANIFEST}"
CORRECTED_RUNTIME_EXPORT="$(runtime_export_for "${CORRECTED_RUN_ROOT}")"
export ECODA_RUN_ROOT="${CORRECTED_RUN_ROOT}" ECODA_RUN_ID=corrected-retry
export ECODA_LOGS_DIR="${CORRECTED_RUN_ROOT}/logs"
export ANALYSIS_VARIANT=corrected_final ANALYSIS_PASS=corrected
export ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION=recovery_35row
export ANALYSIS_ROOT="${HPC_ROOT}/batch_effect/corrected_final/recovery_35row"
export ANALYSIS_NAS_ROOT="${NAS_ROOT}/batch_effect/corrected_final/recovery_35row"
export ANALYSIS_LOG_PREFIX="execution_times_batch_effect_corrected_final_"
export SBATCH_RESULT_ID=7002
: > "${CAPTURE}"
CORRECTED_OUTPUT="$(bash "${SOURCE_TREE}/src/5_run_benchmark_methods/matrix_watchdog.sh" \
  "${CORRECTED_RUN_ROOT}" composition "${CORRECTED_MANIFEST}" 7001 128G 256G \
  shared-cpu 1 "${WORKER_SCRIPT}" "${CORRECTED_RUNTIME_EXPORT}" --cpus-per-task=1)"
CORRECTED_RETRY_MANIFEST="${CORRECTED_RUN_ROOT}/manifests/composition_corrected_final.retry_1.tsv"
CORRECTED_RETRY_LOG_PREFIX="${CORRECTED_RUN_ROOT}/logs/5_matrix_composition_corrected_final_retry1"
[[ "$(sed -n 's/^STATE=//p' "${CORRECTED_RUN_ROOT}/status/watchdogs/composition.status")" == "OK" ]]
case "${CORRECTED_OUTPUT}" in
  *"BATCH_EFFECT_RETRY_ARRAY_JOB_ID=7002"*) ;;
  *) echo "corrected-final retry marker missing" >&2; exit 1 ;;
esac
[[ -s "${CORRECTED_RETRY_MANIFEST}" ]]
CORRECTED_CALLS="$(<"${CAPTURE}")"
for required_identity in \
  "ANALYSIS_VARIANT=corrected_final" \
  "ANALYSIS_PASS=corrected" \
  "ANALYSIS_ROOT=${HPC_ROOT}/batch_effect/corrected_final/recovery_35row" \
  "ANALYSIS_NAS_ROOT=${NAS_ROOT}/batch_effect/corrected_final/recovery_35row" \
  "ANALYSIS_LOG_PREFIX=execution_times_batch_effect_corrected_final_" \
  "ECODA_SOURCE_ROOT=${SOURCE_TREE}" \
  "ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}" \
  "ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}" \
  "ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}" \
  "ECODA_RUNTIME_MODE=apptainer" \
  "ECODA_RUNTIME_PROFILE=stage5" \
  "ECODA_APPTAINER_NV=0" \
  "ECODA_SOURCE_SNAPSHOT_REQUIRED=1" \
  "ECODA_RUN_ID=corrected-retry"; do
  case "${CORRECTED_CALLS}" in
    *"${required_identity}"*) ;;
    *) echo "corrected-final retry omitted ${required_identity}" >&2; exit 1 ;;
  esac
done
for required_path in \
  "MATRIX_RETRY_MANIFEST=${CORRECTED_RETRY_MANIFEST}" \
  "ANALYSIS_MANIFEST=${CORRECTED_RETRY_MANIFEST}" \
  "JOB_LOG_PREFIX=${CORRECTED_RETRY_LOG_PREFIX}" \
  "--output=${CORRECTED_RETRY_LOG_PREFIX}_%A_%a.log" \
  "--error=${CORRECTED_RETRY_LOG_PREFIX}_%A_%a.err"; do
  case "${CORRECTED_CALLS}" in
    *"${required_path}"*) ;;
    *) echo "corrected-final retry omitted ${required_path}" >&2; exit 1 ;;
  esac
done

LEGACY_MANIFEST="${LEGACY_RUN_ROOT}/manifests/composition.tsv"
printf 'Joanito\tbatch_effect_uncorrected\tcomposition\n' > "${LEGACY_MANIFEST}"
LEGACY_RUNTIME_EXPORT="$(runtime_export_for "${LEGACY_RUN_ROOT}")"
unset ANALYSIS_VARIANT ANALYSIS_PASS
export ECODA_RUN_ROOT="${LEGACY_RUN_ROOT}" ECODA_RUN_ID=legacy-retry
export ECODA_LOGS_DIR="${LEGACY_RUN_ROOT}/logs"
export ANALYSIS_ROOT="${HPC_ROOT}/batch_effect/benchmark"
export ANALYSIS_NAS_ROOT="${NAS_ROOT}/batch_effect/benchmark"
export ANALYSIS_LOG_PREFIX="execution_times_"
export SBATCH_RESULT_ID=7102
: > "${CAPTURE}"
LEGACY_OUTPUT="$(bash "${SOURCE_TREE}/src/5_run_benchmark_methods/matrix_watchdog.sh" \
  "${LEGACY_RUN_ROOT}" composition "${LEGACY_MANIFEST}" 7101 128G 256G \
  shared-cpu 1 "${WORKER_SCRIPT}" "${LEGACY_RUNTIME_EXPORT}" --cpus-per-task=1)"
LEGACY_RETRY_MANIFEST="${LEGACY_RUN_ROOT}/manifests/composition.retry_1.tsv"
LEGACY_RETRY_LOG_PREFIX="${LEGACY_RUN_ROOT}/logs/5_matrix_composition_retry1"
[[ "$(sed -n 's/^STATE=//p' "${LEGACY_RUN_ROOT}/status/watchdogs/composition.status")" == "OK" ]]
case "${LEGACY_OUTPUT}" in
  *"MATRIX_RETRY_ARRAY_JOB_ID=7102"*) ;;
  *) echo "legacy retry marker missing" >&2; exit 1 ;;
esac
[[ -s "${LEGACY_RETRY_MANIFEST}" ]]
LEGACY_CALLS="$(<"${CAPTURE}")"
case "${LEGACY_CALLS}" in
  *corrected_final*) echo "legacy retry gained corrected-final identity" >&2; exit 1 ;;
esac
for required_path in \
  "MATRIX_RETRY_MANIFEST=${LEGACY_RETRY_MANIFEST}" \
  "ANALYSIS_MANIFEST=${LEGACY_RETRY_MANIFEST}" \
  "JOB_LOG_PREFIX=${LEGACY_RETRY_LOG_PREFIX}" \
  "--output=${LEGACY_RETRY_LOG_PREFIX}_%A_%a.log" \
  "--error=${LEGACY_RETRY_LOG_PREFIX}_%A_%a.err"; do
  case "${LEGACY_CALLS}" in
    *"${required_path}"*) ;;
    *) echo "legacy retry omitted ${required_path}" >&2; exit 1 ;;
  esac
done
