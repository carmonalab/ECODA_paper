#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
unset HPC_SCRATCH_DIR ECODA_SOURCE_ROOT ECODA_SOURCE_MANIFEST \
  ECODA_SOURCE_SNAPSHOT_REQUIRED ECODA_RUNTIME_IMAGE ECODA_RUNTIME_MANIFEST \
  ECODA_RUNTIME_IDENTITY ECODA_RUN_ROOT ECODA_RUN_ID \
  ANALYSIS_VARIANT ANALYSIS_PASS ANALYSIS_ROOT ANALYSIS_NAS_ROOT \
  ANALYSIS_LOG_PREFIX
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-selection.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
HPC_ROOT="${TMP_DIR}/home/scratch/ECODA_paper"
SOURCE_COMMIT="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SNAPSHOT_ROOT="${HPC_ROOT}/_ecoda_source_snapshots/${SOURCE_COMMIT}"
SOURCE_TREE="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY_DIR="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY_DIR}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY_DIR}/source.tar"
HOST_PREFIX="${TMP_DIR}/host/.pixi/envs/py-cuda13"
RUNTIME_DIR="${HPC_ROOT}/_ecoda_runtime/stage5-test"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
TEST_TMP="${TMP_DIR}/node-tmp"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${HOST_PREFIX}/bin" \
  "${HOST_PREFIX}/lib/R" "${SOURCE_TREE}" "${SOURCE_IDENTITY_DIR}" \
  "${RUNTIME_DIR}" "${TEST_TMP}"
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
chmod 444 "${SOURCE_MANIFEST}" "${SOURCE_ARCHIVE}" "${SNAPSHOT_ROOT}/COMPLETE"
chmod 555 "${SOURCE_IDENTITY_DIR}" "${SNAPSHOT_ROOT}"
printf 'format-2-test-image\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA="$(sha256_file "${RUNTIME_IMAGE}")"
printf 'FORMAT=2\nIMAGE_BUILD_GIT_REVISION=build-commit-b\nIMAGE_PATH=%s\nIMAGE_SHA256=%s\nRUNTIME_ENV=py-cuda13\nRUNTIME_LAYOUT=relocated\nCONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13\nBASE_IMAGE=rockylinux:9\nPIXITAINER_VERSION=0.8.3\nPIXI_VERSION=0.49.0\nAPPTAINER_VERSION=1.3.2\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
  "${RUNTIME_IMAGE}" "${RUNTIME_IMAGE_SHA}" "${SOURCE_TOML_SHA}" \
  "${SOURCE_LOCK_SHA}" > "${RUNTIME_MANIFEST}"
chmod 444 "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}"
chmod 555 "${RUNTIME_DIR}"
cat > "${TMP_DIR}/bin/apptainer" <<'STUB'
#!/bin/bash
case "${1:-}" in
  inspect|--version) exit 0 ;;
  *) exit 0 ;;
esac
STUB
chmod +x "${TMP_DIR}/bin/apptainer"
TEST_LOGS="${TMP_DIR}/logs"
mkdir -p "${TEST_LOGS}"
export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_SCRATCH_ROOT="${HPC_ROOT}"
export ECODA_LOGS_DIR="${TEST_LOGS}" TMPDIR="${TEST_TMP}"
export ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_TREE}/aux"
export ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
export ECODA_HOST_PYTHON_BIN="${HOST_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${HOST_PREFIX}/bin/Rscript --vanilla"
export ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" ECODA_RUNTIME_PROFILE=stage5
export ECODA_APPTAINER_NV=0 APPTAINER_BIN="${TMP_DIR}/bin/apptainer"

CAPTURE="${TMP_DIR}/calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '79000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
printf 'Adams\tbenchmark_analysis\tmrvi\nBassez\tbenchmark_analysis\tgloscope\nKfoury\ttrans\ttrans\nKim\tzeroimp\tzeroimp\n' > "${TMP_DIR}/selection.tsv"
OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
  bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" --selection-file "${TMP_DIR}/selection.tsv"
)"
RUN_ID="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^BENCHMARK_RUN_ID=//p')"
test -n "${RUN_ID}"
RUN_ROOT="${HPC_ROOT}/_ecoda_runs/${RUN_ID}"
SELECTION="${RUN_ROOT}/manifests/selection.tsv"
EXPECTED_SELECTION=$'Adams\tbenchmark_analysis\tmrvi\nBassez\tbenchmark_analysis\tgloscope\nKfoury\tbenchmark_analysis\ttrans\nKim\tbenchmark_analysis\tzeroimp'
test "$(wc -l < "${SELECTION}" | tr -d '[:space:]')" = 4
test "$(cat "${SELECTION}")" = "${EXPECTED_SELECTION}"
cmp -s "${RUN_ROOT}/manifests/source.manifest" "${SOURCE_MANIFEST}"
RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
test -s "${RUNTIME_IDENTITY}"
test "$(sed -n 's/^FORMAT=//p' "${SOURCE_MANIFEST}")" = 1
test "$(sed -n 's/^FORMAT=//p' "${RUNTIME_MANIFEST}")" = 2
test "$(wc -l < "${RUNTIME_IDENTITY}" | tr -d '[:space:]')" = 8
test "$(sed -n 's/^RUNTIME_IMAGE=//p' "${RUNTIME_IDENTITY}")" = "${RUNTIME_IMAGE}"
test "$(sed -n 's/^RUNTIME_MANIFEST=//p' "${RUNTIME_IDENTITY}")" = "${RUNTIME_MANIFEST}"

MRVI_GPU="${RUN_ROOT}/manifests/matrix_benchmark_analysis_mrvi_default_gpu.tsv"
MRVI_CPU="${RUN_ROOT}/manifests/matrix_benchmark_analysis_mrvi_cpu.tsv"
GLOSCOPE_CPU="${RUN_ROOT}/manifests/matrix_benchmark_analysis_gloscope_cpu.tsv"
TRANS="${RUN_ROOT}/manifests/matrix_benchmark_analysis_trans.tsv"
ZEROIMP="${RUN_ROOT}/manifests/matrix_benchmark_analysis_zeroimp.tsv"
test "$(wc -l < "${MRVI_GPU}" | tr -d '[:space:]')" = 1
test "$(cat "${MRVI_GPU}")" = $'Adams\tbenchmark_analysis\tmrvi\thvg2000'
test "$(wc -l < "${MRVI_CPU}" | tr -d '[:space:]')" = 2
test "$(cat "${MRVI_CPU}")" = $'Adams\tbenchmark_analysis\tmrvi\thvg1000\nAdams\tbenchmark_analysis\tmrvi\thvg3000'
test "$(wc -l < "${GLOSCOPE_CPU}" | tr -d '[:space:]')" = 5
test "$(cat "${GLOSCOPE_CPU}")" = $'Bassez\tbenchmark_analysis\tgloscope\thvg2000_pcadims10\nBassez\tbenchmark_analysis\tgloscope\thvg2000_pcadims30\nBassez\tbenchmark_analysis\tgloscope\thvg2000_pcadims50\nBassez\tbenchmark_analysis\tgloscope\thvg1000_pcadims30\nBassez\tbenchmark_analysis\tgloscope\thvg3000_pcadims30'
test "$(wc -l < "${TRANS}" | tr -d '[:space:]')" = 1
test "$(cat "${TRANS}")" = $'Kfoury\tbenchmark_analysis\ttrans'
test "$(wc -l < "${ZEROIMP}" | tr -d '[:space:]')" = 1
test "$(cat "${ZEROIMP}")" = $'Kim\tbenchmark_analysis\tzeroimp'

test "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" = 11
test "$(grep -c -- '--array=1-1' "${CAPTURE}")" = 3
test "$(grep -c -- '--array=1-2' "${CAPTURE}")" = 1
test "$(grep -c -- '--array=1-5' "${CAPTURE}")" = 1
test "$(grep -c 'matrix_gate.sh' "${CAPTURE}")" = 1
echo "exact benchmark selection and aliases: OK"
expect_submit_failure() {
  local label="$1"
  shift
  local before output
  before="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
  if output="$(
    HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
      BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
      bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
      "$@" 2>&1
  )"; then
    echo "expected Stage 5 selection failure: ${label}" >&2
    exit 1
  fi
  if [[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" != "${before}" ]]; then
    echo "invalid Stage 5 selection reached scheduler: ${label}" >&2
    printf '%s\n' "${output}" >&2
    exit 1
  fi
}

NAS_ROOT="${TMP_DIR}/nas/project"
mkdir -p "${NAS_ROOT}"
export NAS_TARGET_DIR="${NAS_ROOT}"
FINAL_METHODS="prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot"
FINAL_SELECTION="${TMP_DIR}/final-selection.tsv"
printf 'Covid19_PBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nKidney_KPMP_full\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n' \
  > "${FINAL_SELECTION}"

FINAL_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${FINAL_SELECTION}" \
    --pass uncorrected \
    --analysis-variant final \
    --methods "${FINAL_METHODS}"
)"
FINAL_RUN_ID="$(printf '%s\n' "${FINAL_OUTPUT}" | sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
test -n "${FINAL_RUN_ID}"
FINAL_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/${FINAL_RUN_ID}"
FINAL_RUN_SELECTION="${FINAL_RUN_ROOT}/manifests/selection.tsv"
test "$(cat "${FINAL_RUN_SELECTION}")" = "$(cat "${FINAL_SELECTION}")"
test "$(wc -l < "${FINAL_RUN_SELECTION}" | tr -d '[:space:]')" = 5
FINAL_METADATA="${FINAL_RUN_ROOT}/metadata"
grep -q '^ANALYSIS_VARIANT=final$' "${FINAL_METADATA}"
grep -q '^ANALYSIS_ROOT=.*/batch_effect/uncorrected_final$' "${FINAL_METADATA}"
grep -q "^ANALYSIS_NAS_ROOT=${NAS_ROOT}/batch_effect/uncorrected_final$" "${FINAL_METADATA}"
grep -q '^ANALYSIS_PASS=uncorrected$' "${FINAL_METADATA}"
grep -q '^PASS=uncorrected$' "${FINAL_METADATA}"
grep -q '^ROOT=.*/batch_effect/uncorrected_final$' "${FINAL_METADATA}"
grep -q '^ANALYSIS_LOG_PREFIX=execution_times_batch_effect_uncorrected_final_$' \
  "${FINAL_METADATA}"
FINAL_PENDING="$(sed -n 's/^PENDING_SELECTION=//p' "${FINAL_METADATA}")"
test -s "${FINAL_PENDING}"
test "$(wc -l < "${FINAL_PENDING}" | tr -d '[:space:]')" = 35
grep -q 'ANALYSIS_ROOT=.*/batch_effect/uncorrected_final' "${CAPTURE}"
if grep -Eq 'ANALYSIS_ROOT=.*/batch_effect/uncorrected(,|$)' "${CAPTURE}"; then
  echo "final worker selection leaked the legacy analysis root" >&2
  exit 1
fi
expect_submit_failure "final variant with broad force" \
  --selection-file "${FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant final --methods "${FINAL_METHODS}" --force
expect_submit_failure "final variant without explicit selection" \
  --pass uncorrected --analysis-variant final --methods "${FINAL_METHODS}"
expect_submit_failure "final variant without explicit method suite" \
  --selection-file "${FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant final
expect_submit_failure "final variant with broad dataset selection" \
  --datasets Covid19_PBMC --pass uncorrected --analysis-variant final \
  --methods "${FINAL_METHODS}"
expect_submit_failure "final variant with corrected pass" \
  --selection-file "${FINAL_SELECTION}" --pass corrected \
  --analysis-variant final --methods "${FINAL_METHODS}"
expect_submit_failure "final variant with ordinary selection" \
  --selection-file "${TMP_DIR}/selection.tsv" --analysis-variant final \
  --methods "${FINAL_METHODS}"
expect_submit_failure "final variant with forbidden method" \
  --selection-file "${FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant final \
  --methods "prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,mofa"
expect_submit_failure "final variant with ordinary analyses" \
  --selection-file "${FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant final --analyses trans
expect_submit_failure "final variant with exact historical mode" \
  --selection-file "${FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant final --exact-batch-selection
CORRECTED_FINAL_SELECTION="${TMP_DIR}/corrected-final-selection.tsv"
printf 'Breast_cancer\tbatch_effect_corrected\tbatch_effect_corrected\nJoanito\tbatch_effect_corrected\tbatch_effect_corrected\nStephenson\tbatch_effect_corrected\tbatch_effect_corrected\nCovid19_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nKidney_KPMP_full\tbatch_effect_corrected\tbatch_effect_corrected\nDiabetes\tbatch_effect_corrected\tbatch_effect_corrected\nLupus_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nLung\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${CORRECTED_FINAL_SELECTION}"
EXPECTED_CORRECTED_FINAL_SELECTION=$'Breast_cancer\tbatch_effect_corrected\tbatch_effect_corrected\nJoanito\tbatch_effect_corrected\tbatch_effect_corrected\nStephenson\tbatch_effect_corrected\tbatch_effect_corrected\nCovid19_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nKidney_KPMP_full\tbatch_effect_corrected\tbatch_effect_corrected\nDiabetes\tbatch_effect_corrected\tbatch_effect_corrected\nLupus_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nLung\tbatch_effect_corrected\tbatch_effect_corrected'
test "$(cat "${CORRECTED_FINAL_SELECTION}")" = \
  "${EXPECTED_CORRECTED_FINAL_SELECTION}"
test "$(wc -l < "${CORRECTED_FINAL_SELECTION}" | tr -d '[:space:]')" = 8

CORRECTED_FINAL_METHOD_MATRIX="${TMP_DIR}/corrected-final-method-matrix.tsv"
printf 'Breast_cancer\tbatch_effect_corrected\tprepare_pseudobulk\nBreast_cancer\tbatch_effect_corrected\tpseudobulk\nBreast_cancer\tbatch_effect_corrected\tgloscope\nBreast_cancer\tbatch_effect_corrected\tcomposition\nBreast_cancer\tbatch_effect_corrected\tmrvi\nBreast_cancer\tbatch_effect_corrected\tpilot\nBreast_cancer\tbatch_effect_corrected\tqot\nJoanito\tbatch_effect_corrected\tprepare_pseudobulk\nJoanito\tbatch_effect_corrected\tpseudobulk\nJoanito\tbatch_effect_corrected\tgloscope\nJoanito\tbatch_effect_corrected\tcomposition\nStephenson\tbatch_effect_corrected\tprepare_pseudobulk\nStephenson\tbatch_effect_corrected\tpseudobulk\nStephenson\tbatch_effect_corrected\tgloscope\nStephenson\tbatch_effect_corrected\tcomposition\nCovid19_PBMC\tbatch_effect_corrected\tprepare_pseudobulk\nCovid19_PBMC\tbatch_effect_corrected\tpseudobulk\nCovid19_PBMC\tbatch_effect_corrected\tgloscope\nCovid19_PBMC\tbatch_effect_corrected\tcomposition\nKidney_KPMP_full\tbatch_effect_corrected\tprepare_pseudobulk\nKidney_KPMP_full\tbatch_effect_corrected\tpseudobulk\nKidney_KPMP_full\tbatch_effect_corrected\tgloscope\nKidney_KPMP_full\tbatch_effect_corrected\tcomposition\nDiabetes\tbatch_effect_corrected\tprepare_pseudobulk\nDiabetes\tbatch_effect_corrected\tpseudobulk\nDiabetes\tbatch_effect_corrected\tgloscope\nDiabetes\tbatch_effect_corrected\tcomposition\nLupus_PBMC\tbatch_effect_corrected\tprepare_pseudobulk\nLupus_PBMC\tbatch_effect_corrected\tpseudobulk\nLupus_PBMC\tbatch_effect_corrected\tgloscope\nLupus_PBMC\tbatch_effect_corrected\tcomposition\nLung\tbatch_effect_corrected\tprepare_pseudobulk\nLung\tbatch_effect_corrected\tpseudobulk\nLung\tbatch_effect_corrected\tgloscope\nLung\tbatch_effect_corrected\tcomposition\n' \
  > "${CORRECTED_FINAL_METHOD_MATRIX}"
EXPECTED_CORRECTED_FINAL_METHOD_MATRIX=$'Breast_cancer\tbatch_effect_corrected\tprepare_pseudobulk\nBreast_cancer\tbatch_effect_corrected\tpseudobulk\nBreast_cancer\tbatch_effect_corrected\tgloscope\nBreast_cancer\tbatch_effect_corrected\tcomposition\nBreast_cancer\tbatch_effect_corrected\tmrvi\nBreast_cancer\tbatch_effect_corrected\tpilot\nBreast_cancer\tbatch_effect_corrected\tqot\nJoanito\tbatch_effect_corrected\tprepare_pseudobulk\nJoanito\tbatch_effect_corrected\tpseudobulk\nJoanito\tbatch_effect_corrected\tgloscope\nJoanito\tbatch_effect_corrected\tcomposition\nStephenson\tbatch_effect_corrected\tprepare_pseudobulk\nStephenson\tbatch_effect_corrected\tpseudobulk\nStephenson\tbatch_effect_corrected\tgloscope\nStephenson\tbatch_effect_corrected\tcomposition\nCovid19_PBMC\tbatch_effect_corrected\tprepare_pseudobulk\nCovid19_PBMC\tbatch_effect_corrected\tpseudobulk\nCovid19_PBMC\tbatch_effect_corrected\tgloscope\nCovid19_PBMC\tbatch_effect_corrected\tcomposition\nKidney_KPMP_full\tbatch_effect_corrected\tprepare_pseudobulk\nKidney_KPMP_full\tbatch_effect_corrected\tpseudobulk\nKidney_KPMP_full\tbatch_effect_corrected\tgloscope\nKidney_KPMP_full\tbatch_effect_corrected\tcomposition\nDiabetes\tbatch_effect_corrected\tprepare_pseudobulk\nDiabetes\tbatch_effect_corrected\tpseudobulk\nDiabetes\tbatch_effect_corrected\tgloscope\nDiabetes\tbatch_effect_corrected\tcomposition\nLupus_PBMC\tbatch_effect_corrected\tprepare_pseudobulk\nLupus_PBMC\tbatch_effect_corrected\tpseudobulk\nLupus_PBMC\tbatch_effect_corrected\tgloscope\nLupus_PBMC\tbatch_effect_corrected\tcomposition\nLung\tbatch_effect_corrected\tprepare_pseudobulk\nLung\tbatch_effect_corrected\tpseudobulk\nLung\tbatch_effect_corrected\tgloscope\nLung\tbatch_effect_corrected\tcomposition'
test "$(cat "${CORRECTED_FINAL_METHOD_MATRIX}")" = \
  "${EXPECTED_CORRECTED_FINAL_METHOD_MATRIX}"
test "$(wc -l < "${CORRECTED_FINAL_METHOD_MATRIX}" | tr -d '[:space:]')" = 35
test "$(sort -u "${CORRECTED_FINAL_METHOD_MATRIX}" | wc -l | tr -d '[:space:]')" = 35
for corrected_matrix_dataset in Breast_cancer Joanito Stephenson Covid19_PBMC \
  Kidney_KPMP_full Diabetes Lupus_PBMC Lung; do
  test "$(grep -F -c "${corrected_matrix_dataset}" \
    "${CORRECTED_FINAL_METHOD_MATRIX}")" -ge 4
done
test "$(grep -F -c $'Breast_cancer\tbatch_effect_corrected\t' \
  "${CORRECTED_FINAL_METHOD_MATRIX}")" = 7
test "$(grep -F -c $'Joanito\tbatch_effect_corrected\t' \
  "${CORRECTED_FINAL_METHOD_MATRIX}")" = 4

# A valid non-Breast cache is independently owned and recorded so the matrix
# run must skip exactly that one declared row while retaining every Breast row.
rm -rf "${HPC_ROOT}/_ecoda_owners"
MATRIX_VALID_PRODUCER_RUN_ID="corrected-final-valid-prep"
MATRIX_VALID_PRODUCER_ROOT="${HPC_ROOT}/_ecoda_runs/${MATRIX_VALID_PRODUCER_RUN_ID}"
MATRIX_VALID_CACHE="${HPC_ROOT}/batch_effect/corrected_final/recovery_35row/pseudobulks/Joanito_batch_effect_corrected_final_pseudobulk_hvg2000.rds"
mkdir -p "${MATRIX_VALID_PRODUCER_ROOT}/manifests" "$(dirname "${MATRIX_VALID_CACHE}")"
printf 'valid corrected-final prepare cache fixture\n' > "${MATRIX_VALID_CACHE}"
(
  set -euo pipefail
  export HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}"
  export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
  export ECODA_RUN_ROOT="${MATRIX_VALID_PRODUCER_ROOT}" \
    ECODA_RUN_ID="${MATRIX_VALID_PRODUCER_RUN_ID}"
  export ANALYSIS_VARIANT=corrected_final PASS_ARG=corrected \
    ANALYSIS_PASS=corrected \
    ANALYSIS_ROOT="${HPC_ROOT}/batch_effect/corrected_final/recovery_35row" \
    ANALYSIS_NAS_ROOT="${NAS_ROOT}/batch_effect/corrected_final/recovery_35row"
  source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  source "${ROOT}/src/utils/bash/ecoda_stage5_policy.sh"
  ecoda_write_checksum "${MATRIX_VALID_CACHE}" >/dev/null
  ecoda_write_artifact_record "${MATRIX_VALID_CACHE}" \
    stage5_prepare_pseudobulk_hvg2000 "${MATRIX_VALID_PRODUCER_RUN_ID}" >/dev/null
  ecoda_artifact_owner_acquire "${MATRIX_VALID_CACHE}" stage5 \
    "${MATRIX_VALID_PRODUCER_RUN_ID}" 0 0 0 >/dev/null
  ecoda_artifact_owner_set_state "${MATRIX_VALID_CACHE}" OK \
    "valid corrected-final matrix fixture"
)
CORRECTED_FINAL_CAPTURE_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
CORRECTED_FINAL_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${CORRECTED_FINAL_SELECTION}" \
    --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" \
    --pass corrected \
    --analysis-variant corrected_final \
    --methods "${FINAL_METHODS}"
)"
CORRECTED_FINAL_RUN_ID="$(printf '%s\n' "${CORRECTED_FINAL_OUTPUT}" |
  sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
test -n "${CORRECTED_FINAL_RUN_ID}"
CORRECTED_FINAL_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/${CORRECTED_FINAL_RUN_ID}"
CORRECTED_FINAL_RUN_SELECTION="${CORRECTED_FINAL_RUN_ROOT}/manifests/selection.tsv"
CORRECTED_FINAL_RUN_MATRIX="${CORRECTED_FINAL_RUN_ROOT}/manifests/method_matrix.tsv"
test "$(cat "${CORRECTED_FINAL_RUN_SELECTION}")" = \
  "${EXPECTED_CORRECTED_FINAL_SELECTION}"
test "$(cat "${CORRECTED_FINAL_RUN_MATRIX}")" = \
  "${EXPECTED_CORRECTED_FINAL_METHOD_MATRIX}"
test -s "${CORRECTED_FINAL_RUN_MATRIX}.md5"
test "$(wc -l < "${CORRECTED_FINAL_RUN_MATRIX}" | tr -d '[:space:]')" = 35
test "$(sort -u "${CORRECTED_FINAL_RUN_MATRIX}" | wc -l | tr -d '[:space:]')" = 35
CORRECTED_FINAL_METADATA="${CORRECTED_FINAL_RUN_ROOT}/metadata"
grep -q '^ANALYSIS_VARIANT=corrected_final$' "${CORRECTED_FINAL_METADATA}"
grep -q '^ANALYSIS_ROOT=.*/batch_effect/corrected_final/recovery_35row$' \
  "${CORRECTED_FINAL_METADATA}"
grep -q "^ANALYSIS_NAS_ROOT=${NAS_ROOT}/batch_effect/corrected_final/recovery_35row$" \
  "${CORRECTED_FINAL_METADATA}"
grep -q '^ANALYSIS_PASS=corrected$' "${CORRECTED_FINAL_METADATA}"
grep -q '^PASS=corrected$' "${CORRECTED_FINAL_METADATA}"
grep -q '^METHODS='"${FINAL_METHODS}"'$' "${CORRECTED_FINAL_METADATA}"
grep -q '^METHOD_MATRIX=.*/manifests/method_matrix.tsv$' \
  "${CORRECTED_FINAL_METADATA}"
grep -q "^METHOD_MATRIX_SOURCE=${CORRECTED_FINAL_METHOD_MATRIX}$" \
  "${CORRECTED_FINAL_METADATA}"
test "${CORRECTED_FINAL_RUN_MATRIX}" != "${CORRECTED_FINAL_METHOD_MATRIX}"
grep -q '^ANALYSIS_ROOT_VERSION=recovery_35row$' \
  "${CORRECTED_FINAL_METADATA}"
grep -q '^ANALYSIS_ROOT_IDENTITY=corrected_final/recovery_35row$' \
  "${CORRECTED_FINAL_METADATA}"
grep -q '^METHOD_MATRIX_MD5=[0-9a-f]\{32\}$' "${CORRECTED_FINAL_METADATA}"
grep -q '^METHOD_MATRIX_SIZE=[1-9][0-9]*$' "${CORRECTED_FINAL_METADATA}"
CORRECTED_FINAL_MATRIX_SHA="$(sha256_file "${CORRECTED_FINAL_RUN_MATRIX}")"
grep -q "^METHOD_MATRIX_SHA256=${CORRECTED_FINAL_MATRIX_SHA}$" \
  "${CORRECTED_FINAL_METADATA}"
grep -q "^METHOD_MATRIX_IDENTITY=${CORRECTED_FINAL_MATRIX_SHA}$" \
  "${CORRECTED_FINAL_METADATA}"
CORRECTED_FINAL_MATRIX_MD5="$(sed -n 's/^METHOD_MATRIX_MD5=//p' \
  "${CORRECTED_FINAL_METADATA}")"
CORRECTED_FINAL_MATRIX_SIZE="$(sed -n 's/^METHOD_MATRIX_SIZE=//p' \
  "${CORRECTED_FINAL_METADATA}")"
test "$(sed -n 's/^MD5=//p' "${CORRECTED_FINAL_RUN_MATRIX}.md5")" = \
  "${CORRECTED_FINAL_MATRIX_MD5}"
test "${CORRECTED_FINAL_MATRIX_SIZE}" = \
  "$(wc -c < "${CORRECTED_FINAL_RUN_MATRIX}" | tr -d '[:space:]')"
grep -q '^DECLARED_METHOD_ROWS=35$' "${CORRECTED_FINAL_METADATA}"
grep -q '^PENDING_METHOD_ROWS=34$' "${CORRECTED_FINAL_METADATA}"
CORRECTED_FINAL_PENDING="$(sed -n 's/^PENDING_SELECTION=//p' \
  "${CORRECTED_FINAL_METADATA}")"
test -s "${CORRECTED_FINAL_PENDING}"
test -s "${CORRECTED_FINAL_PENDING}.md5"
test "$(wc -l < "${CORRECTED_FINAL_PENDING}" | tr -d '[:space:]')" = 34
EXPECTED_CORRECTED_FINAL_PENDING="$(grep -v -F \
  $'Joanito\tbatch_effect_corrected\tprepare_pseudobulk' \
  "${CORRECTED_FINAL_METHOD_MATRIX}")"
test "$(cat "${CORRECTED_FINAL_PENDING}")" = \
  "${EXPECTED_CORRECTED_FINAL_PENDING}"
if grep -F -q \
    $'Joanito\tbatch_effect_corrected\tprepare_pseudobulk' \
    "${CORRECTED_FINAL_PENDING}"; then
  echo "valid non-Breast corrected-final row was not skipped" >&2
  exit 1
fi
grep -F -q $'Breast_cancer\tbatch_effect_corrected\tprepare_pseudobulk' \
  "${CORRECTED_FINAL_PENDING}"
for corrected_pending_method in prepare_pseudobulk pseudobulk gloscope \
  composition mrvi pilot qot; do
  grep -F -q $'\t'"${corrected_pending_method}" \
    "${CORRECTED_FINAL_RUN_MATRIX}"
done
for corrected_pending_method in prepare_pseudobulk pseudobulk gloscope \
  composition mrvi pilot qot; do
  corrected_worker_matrix="${CORRECTED_FINAL_RUN_ROOT}/manifests/matrix_batch_effect_corrected_${corrected_pending_method}.tsv"
  test -s "${corrected_worker_matrix}"
  test -s "${corrected_worker_matrix}.md5"
  case "${corrected_pending_method}" in
    prepare_pseudobulk) corrected_expected_rows=7 ;;
    pseudobulk|gloscope|composition) corrected_expected_rows=8 ;;
    mrvi|pilot|qot) corrected_expected_rows=1 ;;
  esac
  test "$(wc -l < "${corrected_worker_matrix}" | tr -d '[:space:]')" = \
    "${corrected_expected_rows}"
  while IFS= read -r corrected_worker_row; do
    grep -F -q "${corrected_worker_row}" "${CORRECTED_FINAL_PENDING}"
  done < "${corrected_worker_matrix}"
done
CORRECTED_FINAL_CAPTURE="${TMP_DIR}/corrected-final-calls"
sed -n "$((CORRECTED_FINAL_CAPTURE_BEFORE + 1)),\$p" "${CAPTURE}" \
  > "${CORRECTED_FINAL_CAPTURE}"
grep -q 'ANALYSIS_ROOT=.*/batch_effect/corrected_final' \
  "${CORRECTED_FINAL_CAPTURE}"
grep -q 'ANALYSIS_VARIANT=corrected_final' "${CORRECTED_FINAL_CAPTURE}"
if grep -Eq 'ANALYSIS_ROOT=.*/batch_effect/(uncorrected|uncorrected_final)' \
    "${CORRECTED_FINAL_CAPTURE}"; then
  echo "corrected-final worker selection leaked an uncorrected analysis root" >&2
  exit 1
fi
if grep -Eq 'ANALYSIS_ROOT=.*/batch_effect/corrected_final(,|$)' \
    "${CORRECTED_FINAL_CAPTURE}"; then
  echo "corrected-final matrix worker leaked the historical root" >&2
  exit 1
fi
grep -q "ECODA_STAGE5_METHOD_MATRIX=${CORRECTED_FINAL_RUN_MATRIX}" \
  "${CORRECTED_FINAL_CAPTURE}"
echo "corrected-final explicit 35-row method matrix contract: OK"

CORRECTED_FINAL_MATRIX_DUPLICATE="${TMP_DIR}/corrected-final-matrix-duplicate.tsv"
cp "${CORRECTED_FINAL_METHOD_MATRIX}" "${CORRECTED_FINAL_MATRIX_DUPLICATE}"
printf 'Breast_cancer\tbatch_effect_corrected\tprepare_pseudobulk\n' \
  >> "${CORRECTED_FINAL_MATRIX_DUPLICATE}"
expect_submit_failure "corrected-final matrix duplicate row" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_MATRIX_DUPLICATE}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_MATRIX_UNKNOWN="${TMP_DIR}/corrected-final-matrix-unknown.tsv"
sed '1s/Breast_cancer/Unknown/' "${CORRECTED_FINAL_METHOD_MATRIX}" \
  > "${CORRECTED_FINAL_MATRIX_UNKNOWN}"
expect_submit_failure "corrected-final matrix unknown dataset" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_MATRIX_UNKNOWN}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_MATRIX_FORBIDDEN="${TMP_DIR}/corrected-final-matrix-alzheimer.tsv"
sed '1s/Breast_cancer/Alzheimer/' "${CORRECTED_FINAL_METHOD_MATRIX}" \
  > "${CORRECTED_FINAL_MATRIX_FORBIDDEN}"
expect_submit_failure "corrected-final matrix forbidden Alzheimer dataset" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_MATRIX_FORBIDDEN}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_MATRIX_WRONG_VIEW="${TMP_DIR}/corrected-final-matrix-view.tsv"
sed '1s/batch_effect_corrected/batch_effect_uncorrected/' \
  "${CORRECTED_FINAL_METHOD_MATRIX}" > "${CORRECTED_FINAL_MATRIX_WRONG_VIEW}"
expect_submit_failure "corrected-final matrix wrong view" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_MATRIX_WRONG_VIEW}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_MATRIX_BAD_METHOD="${TMP_DIR}/corrected-final-matrix-method.tsv"
sed '7s/qot/mofa/' "${CORRECTED_FINAL_METHOD_MATRIX}" \
  > "${CORRECTED_FINAL_MATRIX_BAD_METHOD}"
expect_submit_failure "corrected-final matrix unsupported method" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_MATRIX_BAD_METHOD}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_SELECTION_DUPLICATE="${TMP_DIR}/corrected-final-selection-duplicate.tsv"
sed '8s/Lung/Breast_cancer/' "${CORRECTED_FINAL_SELECTION}" \
  > "${CORRECTED_FINAL_SELECTION_DUPLICATE}"
expect_submit_failure "corrected-final selection duplicate dataset" \
  --selection-file "${CORRECTED_FINAL_SELECTION_DUPLICATE}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_SELECTION_FORBIDDEN="${TMP_DIR}/corrected-final-selection-forbidden.tsv"
sed '2s/Joanito/Alzheimer/' "${CORRECTED_FINAL_SELECTION}" \
  > "${CORRECTED_FINAL_SELECTION_FORBIDDEN}"
expect_submit_failure "corrected-final selection forbidden dataset" \
  --selection-file "${CORRECTED_FINAL_SELECTION_FORBIDDEN}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_SELECTION_DISABLED="${TMP_DIR}/corrected-final-selection-disabled.tsv"
sed '3s/Stephenson/Myocardial_infarction/' \
  "${CORRECTED_FINAL_SELECTION}" > "${CORRECTED_FINAL_SELECTION_DISABLED}"
expect_submit_failure "corrected-final selection disabled dataset" \
  --selection-file "${CORRECTED_FINAL_SELECTION_DISABLED}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_SELECTION_DEBUG="${TMP_DIR}/corrected-final-selection-debug.tsv"
sed '4s/Covid19_PBMC/_debug/' "${CORRECTED_FINAL_SELECTION}" \
  > "${CORRECTED_FINAL_SELECTION_DEBUG}"
expect_submit_failure "corrected-final selection debug dataset" \
  --selection-file "${CORRECTED_FINAL_SELECTION_DEBUG}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_SELECTION_UNKNOWN="${TMP_DIR}/corrected-final-selection-unknown.tsv"
sed '5s/Kidney_KPMP_full/Unknown/' "${CORRECTED_FINAL_SELECTION}" \
  > "${CORRECTED_FINAL_SELECTION_UNKNOWN}"
expect_submit_failure "corrected-final selection unknown dataset" \
  --selection-file "${CORRECTED_FINAL_SELECTION_UNKNOWN}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final matrix without explicit selection" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "method matrix with legacy ordinary mode" \
  --selection-file "${TMP_DIR}/selection.tsv" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}"
expect_submit_failure "corrected-final matrix with uncorrected pass" \
  --selection-file "${FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" --pass uncorrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final matrix with uncorrected final variant" \
  --selection-file "${FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" --pass uncorrected \
  --analysis-variant final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final matrix without explicit method suite" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" --pass corrected \
  --analysis-variant corrected_final
expect_submit_failure "corrected-final matrix with broad dataset selection" \
  --datasets Breast_cancer --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final matrix with ordinary analyses" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}" --analyses trans
expect_submit_failure "corrected-final matrix with historical exact mode" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}" \
  --exact-batch-selection
expect_submit_failure "corrected-final matrix with force" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" \
  --method-matrix "${CORRECTED_FINAL_METHOD_MATRIX}" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}" --force

# Alzheimer remains a separate ordinary corrected-final one-row follow-up.
ALZHEIMER_CORRECTED_FINAL_SELECTION="${TMP_DIR}/corrected-final-alzheimer-selection.tsv"
printf 'Alzheimer\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${ALZHEIMER_CORRECTED_FINAL_SELECTION}"
ALZHEIMER_CAPTURE_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
ALZHEIMER_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${ALZHEIMER_CORRECTED_FINAL_SELECTION}" \
    --pass corrected \
    --analysis-variant corrected_final \
    --methods "${FINAL_METHODS}"
)"
ALZHEIMER_RUN_ID="$(printf '%s\n' "${ALZHEIMER_OUTPUT}" |
  sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
test -n "${ALZHEIMER_RUN_ID}"
ALZHEIMER_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/${ALZHEIMER_RUN_ID}"
ALZHEIMER_METADATA="${ALZHEIMER_RUN_ROOT}/metadata"
grep -q '^ANALYSIS_VARIANT=corrected_final$' "${ALZHEIMER_METADATA}"
grep -q '^PASS=corrected$' "${ALZHEIMER_METADATA}"
ALZHEIMER_SELECTION="${ALZHEIMER_RUN_ROOT}/manifests/selection.tsv"
test "$(cat "${ALZHEIMER_SELECTION}")" = \
  "$(cat "${ALZHEIMER_CORRECTED_FINAL_SELECTION}")"
ALZHEIMER_PENDING="$(sed -n 's/^PENDING_SELECTION=//p' \
  "${ALZHEIMER_METADATA}")"
test -s "${ALZHEIMER_PENDING}"
test "$(wc -l < "${ALZHEIMER_PENDING}" | tr -d '[:space:]')" = 7
for alzheimer_method in prepare_pseudobulk pseudobulk gloscope composition \
  mrvi pilot qot; do
  test "$(grep -c $'\t'"${alzheimer_method}"$ "${ALZHEIMER_PENDING}")" = 1
  alzheimer_matrix="${ALZHEIMER_RUN_ROOT}/manifests/matrix_batch_effect_corrected_${alzheimer_method}.tsv"
  test -s "${alzheimer_matrix}"
  test "$(wc -l < "${alzheimer_matrix}" | tr -d '[:space:]')" = 1
  test "$(awk -F $'\t' '$1 != "Alzheimer" { bad=1 } END { print bad + 0 }' \
    "${alzheimer_matrix}")" = 0
done
ALZHEIMER_CAPTURE="${TMP_DIR}/alzheimer-calls"
sed -n "$((ALZHEIMER_CAPTURE_BEFORE + 1)),\$p" "${CAPTURE}" \
  > "${ALZHEIMER_CAPTURE}"
test "$(grep -c -- '--array=1-1' "${ALZHEIMER_CAPTURE}")" = 7
if grep -Eq -- '--array=1-(8|9)' "${ALZHEIMER_CAPTURE}"; then
  echo "Alzheimer corrected-final follow-up expanded beyond its one-row scope" >&2
  exit 1
fi
echo "corrected-final Alzheimer follow-up selection contract: OK"

echo "corrected-final benchmark variant selection contract: OK"

echo "final benchmark variant selection contract: OK"
TARGETED_PREP_PRODUCER_RUN_ID="targeted-final-prep-cache"
TARGETED_PREP_ROOT="${HPC_ROOT}/_ecoda_runs/${TARGETED_PREP_PRODUCER_RUN_ID}"
mkdir -p "${TARGETED_PREP_ROOT}/manifests"
for targeted_ds in Covid19_PBMC Diabetes Joanito Lung; do
  targeted_prep="${HPC_ROOT}/batch_effect/uncorrected_final/pseudobulks/${targeted_ds}_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds"
  mkdir -p "$(dirname "${targeted_prep}")"
  printf 'valid targeted final pseudobulk cache\n' > "${targeted_prep}"
  (
    export HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}"
    export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
    export ECODA_RUN_ROOT="${TARGETED_PREP_ROOT}" \
      ECODA_RUN_ID="${TARGETED_PREP_PRODUCER_RUN_ID}"
    source "${ROOT}/src/slurm_config.sh"
    source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
    source "${ROOT}/src/utils/bash/ecoda_stage5_policy.sh"
    ecoda_write_checksum "${targeted_prep}" >/dev/null
    ecoda_write_artifact_record "${targeted_prep}" \
      stage5_prepare_pseudobulk_hvg2000 "${TARGETED_PREP_PRODUCER_RUN_ID}" >/dev/null
    ecoda_artifact_owner_acquire "${targeted_prep}" stage5 \
      "${TARGETED_PREP_PRODUCER_RUN_ID}" 0 0 0 >/dev/null
    ecoda_artifact_owner_set_state "${targeted_prep}" OK \
      "targeted final dependency cache fixture"
  )
done
TARGETED_FINAL_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${FINAL_SELECTION}" --pass uncorrected \
    --analysis-variant final --target-methods gloscope,composition,mrvi
)"
TARGETED_FINAL_RUN_ID="$(printf '%s\n' "${TARGETED_FINAL_OUTPUT}" |
  sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
test -n "${TARGETED_FINAL_RUN_ID}"
TARGETED_FINAL_ROOT="${HPC_ROOT}/_ecoda_runs/${TARGETED_FINAL_RUN_ID}"
TARGETED_FINAL_PENDING="$(sed -n 's/^PENDING_SELECTION=//p' \
  "${TARGETED_FINAL_ROOT}/metadata")"
EXPECTED_TARGETED_FINAL_PENDING=$'Kidney_KPMP_full\tbatch_effect_uncorrected\tprepare_pseudobulk\nCovid19_PBMC\tbatch_effect_uncorrected\tgloscope\nDiabetes\tbatch_effect_uncorrected\tgloscope\nJoanito\tbatch_effect_uncorrected\tgloscope\nLung\tbatch_effect_uncorrected\tgloscope\nKidney_KPMP_full\tbatch_effect_uncorrected\tgloscope\nCovid19_PBMC\tbatch_effect_uncorrected\tcomposition\nDiabetes\tbatch_effect_uncorrected\tcomposition\nJoanito\tbatch_effect_uncorrected\tcomposition\nLung\tbatch_effect_uncorrected\tcomposition\nKidney_KPMP_full\tbatch_effect_uncorrected\tcomposition\nCovid19_PBMC\tbatch_effect_uncorrected\tmrvi\nDiabetes\tbatch_effect_uncorrected\tmrvi\nJoanito\tbatch_effect_uncorrected\tmrvi\nLung\tbatch_effect_uncorrected\tmrvi\nKidney_KPMP_full\tbatch_effect_uncorrected\tmrvi'
test "$(cat "${TARGETED_FINAL_PENDING}")" = "${EXPECTED_TARGETED_FINAL_PENDING}"
test "$(wc -l < "${TARGETED_FINAL_PENDING}" | tr -d '[:space:]')" = 16
test "$(grep -c $'\tgloscope$' "${TARGETED_FINAL_PENDING}")" = 5
test "$(grep -c $'\tcomposition$' "${TARGETED_FINAL_PENDING}")" = 5
test "$(grep -c $'\tmrvi$' "${TARGETED_FINAL_PENDING}")" = 5
echo "targeted final five-row repair selection contract: OK"
rm -rf "${HPC_ROOT}/_ecoda_owners/stage5"
COMPOSITION_ONLY_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${FINAL_SELECTION}" --pass uncorrected \
    --analysis-variant final --target-methods composition
)"
COMPOSITION_ONLY_RUN_ID="$(printf '%s\n' "${COMPOSITION_ONLY_OUTPUT}" |
  sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
test -n "${COMPOSITION_ONLY_RUN_ID}"
COMPOSITION_ONLY_ROOT="${HPC_ROOT}/_ecoda_runs/${COMPOSITION_ONLY_RUN_ID}"
COMPOSITION_ONLY_PENDING="$(sed -n 's/^PENDING_SELECTION=//p' \
  "${COMPOSITION_ONLY_ROOT}/metadata")"
EXPECTED_COMPOSITION_ONLY_PENDING=$'Kidney_KPMP_full\tbatch_effect_uncorrected\tprepare_pseudobulk\nCovid19_PBMC\tbatch_effect_uncorrected\tcomposition\nDiabetes\tbatch_effect_uncorrected\tcomposition\nJoanito\tbatch_effect_uncorrected\tcomposition\nLung\tbatch_effect_uncorrected\tcomposition\nKidney_KPMP_full\tbatch_effect_uncorrected\tcomposition'
test "$(cat "${COMPOSITION_ONLY_PENDING}")" = \
  "${EXPECTED_COMPOSITION_ONLY_PENDING}"
test "$(wc -l < "${COMPOSITION_ONLY_PENDING}" | tr -d '[:space:]')" = 6
test "$(grep -c $'\tcomposition$' "${COMPOSITION_ONLY_PENDING}")" = 5
COMPOSITION_ONLY_CAPTURE="${TMP_DIR}/composition-only-calls"
printf '%s\n' "${COMPOSITION_ONLY_OUTPUT}" |
  grep -E 'BATCH_EFFECT_(ARRAY|WATCHDOG|AGGREGATE_GATE)_JOB_ID=' \
  > "${COMPOSITION_ONLY_CAPTURE}"
echo "targeted composition-only recovery selection contract: OK"
