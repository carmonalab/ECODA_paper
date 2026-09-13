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
FOUR_ROW_FINAL_SELECTION="${TMP_DIR}/four-row-final-selection.tsv"
printf 'Covid19_PBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n' \
  > "${FOUR_ROW_FINAL_SELECTION}"
expect_submit_failure "final variant with retired four-row selection" \
  --selection-file "${FOUR_ROW_FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant final --methods "${FINAL_METHODS}"
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
CAPTURE_BEFORE_CORRECTED="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
expect_submit_failure "final variant with exact historical mode" \
  --selection-file "${FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant final --exact-batch-selection
CORRECTED_FINAL_SELECTION="${TMP_DIR}/corrected-final-selection.tsv"
printf 'Joanito\tbatch_effect_corrected\tbatch_effect_corrected\nStephenson\tbatch_effect_corrected\tbatch_effect_corrected\nAlzheimer\tbatch_effect_corrected\tbatch_effect_corrected\nBreast_cancer\tbatch_effect_corrected\tbatch_effect_corrected\nCovid19_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nKidney_KPMP_full\tbatch_effect_corrected\tbatch_effect_corrected\nDiabetes\tbatch_effect_corrected\tbatch_effect_corrected\nLupus_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nLung\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${CORRECTED_FINAL_SELECTION}"
CORRECTED_FINAL_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${CORRECTED_FINAL_SELECTION}" \
    --pass corrected \
    --analysis-variant corrected_final \
    --methods "${FINAL_METHODS}"
)"
CORRECTED_FINAL_RUN_ID="$(printf '%s\n' "${CORRECTED_FINAL_OUTPUT}" | sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
test -n "${CORRECTED_FINAL_RUN_ID}"
CORRECTED_FINAL_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/${CORRECTED_FINAL_RUN_ID}"
CORRECTED_FINAL_RUN_SELECTION="${CORRECTED_FINAL_RUN_ROOT}/manifests/selection.tsv"
test "$(cat "${CORRECTED_FINAL_RUN_SELECTION}")" = "$(cat "${CORRECTED_FINAL_SELECTION}")"
test "$(wc -l < "${CORRECTED_FINAL_RUN_SELECTION}" | tr -d '[:space:]')" = 9
CORRECTED_FINAL_METADATA="${CORRECTED_FINAL_RUN_ROOT}/metadata"
grep -q '^ANALYSIS_VARIANT=corrected_final$' "${CORRECTED_FINAL_METADATA}"
grep -q '^ANALYSIS_ROOT=.*/batch_effect/corrected_final$' "${CORRECTED_FINAL_METADATA}"
grep -q "^ANALYSIS_NAS_ROOT=${NAS_ROOT}/batch_effect/corrected_final$" \
  "${CORRECTED_FINAL_METADATA}"
grep -q '^ANALYSIS_PASS=corrected$' "${CORRECTED_FINAL_METADATA}"
grep -q '^PASS=corrected$' "${CORRECTED_FINAL_METADATA}"
CORRECTED_CAPTURE="${TMP_DIR}/corrected-calls"
sed -n "$((CAPTURE_BEFORE_CORRECTED + 1)),\$p" "${CAPTURE}" > "${CORRECTED_CAPTURE}"
grep -q 'ANALYSIS_ROOT=.*/batch_effect/corrected_final' "${CORRECTED_CAPTURE}"
grep -q 'ANALYSIS_VARIANT=corrected_final' "${CORRECTED_CAPTURE}"
grep -q 'ANALYSIS_PASS=corrected' "${CORRECTED_CAPTURE}"
grep -q 'ANALYSIS_NAS_ROOT=.*/batch_effect/corrected_final' "${CORRECTED_CAPTURE}"
grep -q 'ANALYSIS_LOG_PREFIX=execution_times_batch_effect_corrected_final_' \
  "${CORRECTED_CAPTURE}"
if grep -Eq 'ANALYSIS_ROOT=.*/batch_effect/(uncorrected|uncorrected_final)(,|$)' \
    "${CORRECTED_CAPTURE}"; then
  echo "corrected-final worker selection leaked an uncorrected analysis root" >&2
  exit 1
fi
CORRECTED_FINAL_SHORT_SELECTION="${TMP_DIR}/corrected-final-short-selection.tsv"
printf 'Joanito\tbatch_effect_corrected\tbatch_effect_corrected\nStephenson\tbatch_effect_corrected\tbatch_effect_corrected\nAlzheimer\tbatch_effect_corrected\tbatch_effect_corrected\nBreast_cancer\tbatch_effect_corrected\tbatch_effect_corrected\nCovid19_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nKidney_KPMP_full\tbatch_effect_corrected\tbatch_effect_corrected\nDiabetes\tbatch_effect_corrected\tbatch_effect_corrected\nLupus_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${CORRECTED_FINAL_SHORT_SELECTION}"
expect_submit_failure "corrected-final variant with incomplete nine-row selection" \
  --selection-file "${CORRECTED_FINAL_SHORT_SELECTION}" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}"
CORRECTED_FINAL_REORDERED_SELECTION="${TMP_DIR}/corrected-final-reordered-selection.tsv"
printf 'Stephenson\tbatch_effect_corrected\tbatch_effect_corrected\nJoanito\tbatch_effect_corrected\tbatch_effect_corrected\nAlzheimer\tbatch_effect_corrected\tbatch_effect_corrected\nBreast_cancer\tbatch_effect_corrected\tbatch_effect_corrected\nCovid19_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nKidney_KPMP_full\tbatch_effect_corrected\tbatch_effect_corrected\nDiabetes\tbatch_effect_corrected\tbatch_effect_corrected\nLupus_PBMC\tbatch_effect_corrected\tbatch_effect_corrected\nLung\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${CORRECTED_FINAL_REORDERED_SELECTION}"
expect_submit_failure "corrected-final variant with wrong config order" \
  --selection-file "${CORRECTED_FINAL_REORDERED_SELECTION}" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final variant with force" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}" --force
expect_submit_failure "corrected-final variant without explicit selection" \
  --pass corrected --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final variant with uncorrected pass" \
  --selection-file "${FINAL_SELECTION}" --pass uncorrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final variant with broad dataset selection" \
  --datasets Joanito --pass corrected --analysis-variant corrected_final \
  --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final variant with ordinary selection" \
  --selection-file "${TMP_DIR}/selection.tsv" --pass corrected \
  --analysis-variant corrected_final --methods "${FINAL_METHODS}"
expect_submit_failure "corrected-final variant with forbidden method" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" --pass corrected \
  --analysis-variant corrected_final \
  --methods "prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,mofa"
expect_submit_failure "corrected-final variant with exact historical mode" \
  --selection-file "${CORRECTED_FINAL_SELECTION}" --pass corrected \
  --analysis-variant corrected_final --exact-batch-selection

echo "corrected-final benchmark variant selection contract: OK"

echo "final benchmark variant selection contract: OK"
md5_file_for_fixture() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

# Exercise the validator-only Kidney legacy inventory independently of the
# all-missing five-row final fixture above.  The arbitrary payload is enough
# for the existing Rscript validator stub; no H5AD/RDS computation is run.
rm -rf "${HPC_ROOT}/_ecoda_owners"
LEGACY_KIDNEY_ROOT="${HPC_ROOT}/batch_effect/uncorrected"
LEGACY_KIDNEY_PREP="${LEGACY_KIDNEY_ROOT}/pseudobulks/Kidney_KPMP_full_batch_effect_uncorrected_pseudobulk_hvg2000.rds"
mkdir -p "$(dirname "${LEGACY_KIDNEY_PREP}")"
printf 'legacy Kidney prepare pseudobulk fixture\n' > "${LEGACY_KIDNEY_PREP}"
(
  set -euo pipefail
  source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  ecoda_write_checksum "${LEGACY_KIDNEY_PREP}" >/dev/null
)
cp "${LEGACY_KIDNEY_PREP}" "${TMP_DIR}/kidney-legacy-prep.before"
cp "${LEGACY_KIDNEY_PREP}.md5" "${TMP_DIR}/kidney-legacy-prep.before.md5"
FIXTURE_CAPTURE_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
KIDNEY_FIXTURE_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${FINAL_SELECTION}" \
    --pass uncorrected \
    --analysis-variant final \
    --methods "${FINAL_METHODS}"
)"
KIDNEY_FIXTURE_RUN_ID="$(printf '%s\n' "${KIDNEY_FIXTURE_OUTPUT}" |
  sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
test -n "${KIDNEY_FIXTURE_RUN_ID}"
KIDNEY_FIXTURE_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/${KIDNEY_FIXTURE_RUN_ID}"
KIDNEY_FIXTURE_METADATA="${KIDNEY_FIXTURE_RUN_ROOT}/metadata"
KIDNEY_FIXTURE_PENDING="$(sed -n 's/^PENDING_SELECTION=//p' \
  "${KIDNEY_FIXTURE_METADATA}")"
test -s "${KIDNEY_FIXTURE_PENDING}"
test "$(wc -l < "${KIDNEY_FIXTURE_PENDING}" | tr -d '[:space:]')" = 34
if grep -F -q $'Kidney_KPMP_full\tbatch_effect_uncorrected\tprepare_pseudobulk' \
    "${KIDNEY_FIXTURE_PENDING}"; then
  echo "validated legacy Kidney preparation appeared in pending selection" >&2
  exit 1
fi
grep -F -q $'Kidney_KPMP_full\tbatch_effect_uncorrected\tqot' \
  "${KIDNEY_FIXTURE_PENDING}"

KIDNEY_FIXTURE_OWNERS="${KIDNEY_FIXTURE_RUN_ROOT}/manifests/owners.tsv"
test -s "${KIDNEY_FIXTURE_OWNERS}"
if grep -F -q \
    $'uncorrected/Kidney_KPMP_full/batch_effect_uncorrected/prepare_pseudobulk' \
    "${KIDNEY_FIXTURE_OWNERS}"; then
  echo "validated legacy Kidney preparation appeared in owner expansion" >&2
  exit 1
fi
grep -F -q \
  $'uncorrected/Kidney_KPMP_full/batch_effect_uncorrected/qot' \
  "${KIDNEY_FIXTURE_OWNERS}"
KIDNEY_FIXTURE_PREP_MATRIX="${KIDNEY_FIXTURE_RUN_ROOT}/manifests/matrix_batch_effect_uncorrected_prepare_pseudobulk.tsv"
KIDNEY_FIXTURE_QOT_MATRIX="${KIDNEY_FIXTURE_RUN_ROOT}/manifests/matrix_batch_effect_uncorrected_qot.tsv"
test "$(wc -l < "${KIDNEY_FIXTURE_PREP_MATRIX}" | tr -d '[:space:]')" = 4
if grep -F -q $'Kidney_KPMP_full\tbatch_effect_uncorrected\tprepare_pseudobulk' \
    "${KIDNEY_FIXTURE_PREP_MATRIX}"; then
  echo "validated legacy Kidney preparation appeared in matrix owner expansion" >&2
  exit 1
fi
test "$(wc -l < "${KIDNEY_FIXTURE_QOT_MATRIX}" | tr -d '[:space:]')" = 5
grep -F -q $'Kidney_KPMP_full\tbatch_effect_uncorrected\tqot' \
  "${KIDNEY_FIXTURE_QOT_MATRIX}"

KIDNEY_FIXTURE_INVENTORY="${KIDNEY_FIXTURE_RUN_ROOT}/manifests/kidney_legacy_inventory.tsv"
test -s "${KIDNEY_FIXTURE_INVENTORY}"
test -s "${KIDNEY_FIXTURE_INVENTORY}.md5"
[[ ! -L "${KIDNEY_FIXTURE_INVENTORY}" &&
   ! -L "${KIDNEY_FIXTURE_INVENTORY}.md5" ]]
test "$(awk -F $'\t' '$1 != "Kidney_KPMP_full" { bad=1 } END { print bad + 0 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = 0
test "$(awk -F $'\t' \
  '{ printf "%s%s", $2, (NR == 7 ? "\n" : ",") }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = "${FINAL_METHODS}"
test "$(wc -l < "${KIDNEY_FIXTURE_INVENTORY}" | tr -d '[:space:]')" = 7
test "$(awk -F $'\t' 'NF != 4 { bad=1 } END { print bad + 0 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = 0
test "$(awk -F $'\t' '$2 == "prepare_pseudobulk" { print $3 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = valid
test "$(awk -F $'\t' '$3 == "valid" { count++ } END { print count + 0 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = 1
test "$(awk -F $'\t' '$3 == "missing" { count++ } END { print count + 0 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = 6
test "$(awk -F $'\t' '$2 == "qot" { print $3 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = missing
test "$(awk -F $'\t' '$2 == "prepare_pseudobulk" { print $4 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = "${LEGACY_KIDNEY_PREP}"
test "$(awk -F $'\t' '$2 == "qot" { print $4 }' \
  "${KIDNEY_FIXTURE_INVENTORY}")" = \
  "${LEGACY_KIDNEY_ROOT}/embeddings/Kidney_KPMP_full_batch_effect_uncorrected_hvg2000_highres_qot_dists.feather"
KIDNEY_FIXTURE_INVENTORY_MD5="$(md5_file_for_fixture "${KIDNEY_FIXTURE_INVENTORY}")"
KIDNEY_FIXTURE_INVENTORY_SIZE="$(wc -c < "${KIDNEY_FIXTURE_INVENTORY}" |
  tr -d '[:space:]')"
test "$(sed -n 's/^MD5=//p' "${KIDNEY_FIXTURE_INVENTORY}.md5")" = \
  "${KIDNEY_FIXTURE_INVENTORY_MD5}"
test "$(sed -n 's/^SIZE=//p' "${KIDNEY_FIXTURE_INVENTORY}.md5")" = \
  "${KIDNEY_FIXTURE_INVENTORY_SIZE}"
test "$(sed -n 's/^PATH=//p' "${KIDNEY_FIXTURE_INVENTORY}.md5")" = \
  "${KIDNEY_FIXTURE_INVENTORY}"
test "$(sed -n 's/^KIDNEY_LEGACY_INVENTORY=//p' \
  "${KIDNEY_FIXTURE_METADATA}")" = "${KIDNEY_FIXTURE_INVENTORY}"
test "$(sed -n 's/^KIDNEY_LEGACY_INVENTORY_MD5=//p' \
  "${KIDNEY_FIXTURE_METADATA}")" = "${KIDNEY_FIXTURE_INVENTORY_MD5}"
test "$(sed -n 's/^KIDNEY_LEGACY_INVENTORY_SIZE=//p' \
  "${KIDNEY_FIXTURE_METADATA}")" = "${KIDNEY_FIXTURE_INVENTORY_SIZE}"
test "$(sed -n 's/^KIDNEY_LEGACY_INVENTORY_STATUS=//p' \
  "${KIDNEY_FIXTURE_METADATA}")" = \
  "prepare_pseudobulk=valid;pseudobulk=missing;gloscope=missing;composition=missing;mrvi=missing;pilot=missing;qot=missing"
test "$(sed -n 's/^KIDNEY_LEGACY_VALID_METHODS=//p' \
  "${KIDNEY_FIXTURE_METADATA}")" = prepare_pseudobulk
test "$(sed -n 's/^KIDNEY_LEGACY_MISSING_METHODS=//p' \
  "${KIDNEY_FIXTURE_METADATA}")" = "pseudobulk gloscope composition mrvi pilot qot"
test "$(sed -n 's/^KIDNEY_LEGACY_INVALID_METHODS=//p' \
  "${KIDNEY_FIXTURE_METADATA}")" = ""

# The missing qot row must resolve to the final stem, while the valid legacy
# preparation must remain outside both final roots and retain its source bytes.
QOT_FINAL_PATH="$(
  set -euo pipefail
  source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  PASS_ARG=uncorrected
  ANALYSIS_PASS=uncorrected
  ANALYSIS_VARIANT=final
  ANALYSIS_ROOT="${HPC_ROOT}/batch_effect/uncorrected_final"
  ANALYSIS_NAS_ROOT="${NAS_ROOT}/batch_effect/uncorrected_final"
  _ecoda_stage5_artifacts_for \
    Kidney_KPMP_full batch_effect_uncorrected qot
  test "${#ECODA_BENCHMARK_ARTIFACTS[@]}" = 1
  printf '%s\n' "${ECODA_BENCHMARK_ARTIFACTS[0]}"
)"
test "${QOT_FINAL_PATH}" = \
  "${HPC_ROOT}/batch_effect/uncorrected_final/embeddings/Kidney_KPMP_full_batch_effect_uncorrected_final_hvg2000_highres_qot_dists.feather"
test "${QOT_FINAL_PATH}" != \
  "${LEGACY_KIDNEY_ROOT}/embeddings/Kidney_KPMP_full_batch_effect_uncorrected_hvg2000_highres_qot_dists.feather"
[[ ! -e "${HPC_ROOT}/batch_effect/uncorrected_final/pseudobulks/Kidney_KPMP_full_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds" &&
   ! -L "${HPC_ROOT}/batch_effect/uncorrected_final/pseudobulks/Kidney_KPMP_full_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds" ]]
[[ ! -e "${NAS_ROOT}/batch_effect/uncorrected_final/pseudobulks/Kidney_KPMP_full_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds" &&
   ! -L "${NAS_ROOT}/batch_effect/uncorrected_final/pseudobulks/Kidney_KPMP_full_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds" ]]
cmp -s "${TMP_DIR}/kidney-legacy-prep.before" "${LEGACY_KIDNEY_PREP}"
cmp -s "${TMP_DIR}/kidney-legacy-prep.before.md5" "${LEGACY_KIDNEY_PREP}.md5"
KIDNEY_FIXTURE_CAPTURE="${TMP_DIR}/kidney-fixture-calls"
sed -n "$((FIXTURE_CAPTURE_BEFORE + 1)),\$p" "${CAPTURE}" \
  > "${KIDNEY_FIXTURE_CAPTURE}"
if grep -F -q "${LEGACY_KIDNEY_PREP}" "${KIDNEY_FIXTURE_CAPTURE}"; then
  echo "legacy Kidney artifact leaked into final scheduler payload" >&2
  exit 1
fi
echo "validator-only Kidney legacy inventory reuse: OK"
# Final targeted recovery may name only the failed method classes across the
# exact five-row final selection.  The legacy Kidney composition and MRVI
# artifacts are valid, so only the five GloScope plus four composition and four
# MRVI rows may enter pending selection.
LEGACY_KIDNEY_COMPOSITION="${LEGACY_KIDNEY_ROOT}/results/Kidney_KPMP_full_batch_effect_uncorrected_composition.rds"
LEGACY_KIDNEY_METADATA="${LEGACY_KIDNEY_ROOT}/results/Kidney_KPMP_full_batch_effect_uncorrected_metadata.rds"
LEGACY_KIDNEY_MRVI="${LEGACY_KIDNEY_ROOT}/embeddings/Kidney_KPMP_full_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"
LEGACY_KIDNEY_MRVI_RUNTIME="${LEGACY_KIDNEY_MRVI}.runtime.json"
for legacy_path in "${LEGACY_KIDNEY_COMPOSITION}" "${LEGACY_KIDNEY_METADATA}" \
  "${LEGACY_KIDNEY_MRVI}" "${LEGACY_KIDNEY_MRVI_RUNTIME}"; do
  mkdir -p "$(dirname "${legacy_path}")"
  printf 'valid legacy targeted fixture\n' > "${legacy_path}"
  digest="$(md5_file_for_fixture "${legacy_path}")"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" \
    "$(wc -c < "${legacy_path}" | tr -d '[:space:]')" "${legacy_path}" \
    > "${legacy_path}.md5"
done
printf '{"dataset":"Kidney_KPMP_full","method":"MrVI_hvg2000"}\n' \
  > "${LEGACY_KIDNEY_MRVI_RUNTIME}"
digest="$(md5_file_for_fixture "${LEGACY_KIDNEY_MRVI_RUNTIME}")"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" \
  "$(wc -c < "${LEGACY_KIDNEY_MRVI_RUNTIME}" | tr -d '[:space:]')" \
  "${LEGACY_KIDNEY_MRVI_RUNTIME}" > "${LEGACY_KIDNEY_MRVI_RUNTIME}.md5"
rm -rf "${HPC_ROOT}/_ecoda_owners"
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
EXPECTED_TARGETED_FINAL_PENDING=$'Covid19_PBMC\tbatch_effect_uncorrected\tgloscope\nDiabetes\tbatch_effect_uncorrected\tgloscope\nJoanito\tbatch_effect_uncorrected\tgloscope\nLung\tbatch_effect_uncorrected\tgloscope\nKidney_KPMP_full\tbatch_effect_uncorrected\tgloscope\nCovid19_PBMC\tbatch_effect_uncorrected\tcomposition\nDiabetes\tbatch_effect_uncorrected\tcomposition\nJoanito\tbatch_effect_uncorrected\tcomposition\nLung\tbatch_effect_uncorrected\tcomposition\nCovid19_PBMC\tbatch_effect_uncorrected\tmrvi\nDiabetes\tbatch_effect_uncorrected\tmrvi\nJoanito\tbatch_effect_uncorrected\tmrvi\nLung\tbatch_effect_uncorrected\tmrvi'
test "$(cat "${TARGETED_FINAL_PENDING}")" = "${EXPECTED_TARGETED_FINAL_PENDING}"
test "$(wc -l < "${TARGETED_FINAL_PENDING}" | tr -d '[:space:]')" = 13
test "$(grep -c $'\tgloscope$' "${TARGETED_FINAL_PENDING}")" = 5
test "$(grep -c $'\tcomposition$' "${TARGETED_FINAL_PENDING}")" = 4
test "$(grep -c $'\tmrvi$' "${TARGETED_FINAL_PENDING}")" = 4
! grep -Eq $'\t(prepare_pseudobulk|pseudobulk|pilot|qot)$' \
  "${TARGETED_FINAL_PENDING}"
echo "targeted final five-row repair selection contract: OK"
