#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
unset HPC_SCRATCH_DIR ECODA_SOURCE_ROOT ECODA_SOURCE_MANIFEST \
  ECODA_SOURCE_SNAPSHOT_REQUIRED ECODA_RUNTIME_IMAGE ECODA_RUNTIME_MANIFEST \
  ECODA_RUNTIME_IDENTITY ECODA_RUN_ROOT ECODA_RUN_ID
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
