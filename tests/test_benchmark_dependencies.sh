#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
unset HPC_SCRATCH_DIR ECODA_SOURCE_ROOT ECODA_SOURCE_MANIFEST \
  ECODA_SOURCE_SNAPSHOT_REQUIRED ECODA_RUNTIME_IMAGE ECODA_RUNTIME_MANIFEST \
  ECODA_RUNTIME_IDENTITY ECODA_RUN_ROOT ECODA_RUN_ID
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-benchmark-deps.XXXXXX")"
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
TEST_LOGS="${TMP_DIR}/ecoda-logs"
TEST_TMP="${TMP_DIR}/node-tmp"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${TMP_DIR}/nas" "${HOST_PREFIX}/bin" \
  "${HOST_PREFIX}/lib/R" "${SOURCE_TREE}" "${SOURCE_IDENTITY_DIR}" \
  "${RUNTIME_DIR}" "${TEST_LOGS}" "${TEST_TMP}"
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
CAPTURE="${TMP_DIR}/calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '78000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_SCRATCH_ROOT="${HPC_ROOT}"
export NAS_PREFIX="${TMP_DIR}/nas" \
  NAS_BASE_DIR="${TMP_DIR}/nas/DataCollections" \
  NAS_SC_DIR="${TMP_DIR}/nas/DataCollections/Standardized_SingleCell_Datasets" \
  NAS_TARGET_DIR="${TMP_DIR}/nas/Projects/ECODA_paper"
export ECODA_LOGS_DIR="${TEST_LOGS}" TMPDIR="${TEST_TMP}"
export ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_TREE}/aux"
export ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
export ECODA_HOST_PYTHON_BIN="${HOST_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${HOST_PREFIX}/bin/Rscript --vanilla"
export ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" ECODA_RUNTIME_PROFILE=stage5
export ECODA_APPTAINER_NV=0 APPTAINER_BIN="${TMP_DIR}/bin/apptainer"
printf '_debug\tbenchmark_analysis\tselected\n_debug\tbatch_effect_uncorrected\tselected\nAdams\tbenchmark_analysis\tselected\nBassez\tbenchmark_analysis\tselected\n' > "${TMP_DIR}/selection.tsv"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid BENCHMARK_MATRIX_TEST=1 \
  bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
  --selection-file "${TMP_DIR}/selection.tsv" --methods mofa,gloscope,mrvi >/dev/null
CALLS="${CAPTURE}"
[[ "$(wc -l < "${CALLS}" | tr -d '[:space:]')" == 19 ]]
RUN_ROOTS=("${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/"*)
[[ ${#RUN_ROOTS[@]} -eq 1 ]]
RUN_ROOT="${RUN_ROOTS[0]}"
MANIFEST_NAMES=(
  matrix_benchmark_analysis_prepare_pseudobulk.tsv
  matrix_batch_effect_uncorrected_prepare_pseudobulk.tsv
  matrix_benchmark_analysis_mofa.tsv
  matrix_batch_effect_uncorrected_mofa.tsv
  matrix_benchmark_analysis_gloscope_cpu.tsv
  matrix_batch_effect_uncorrected_gloscope.tsv
  matrix_benchmark_analysis_mrvi_default_gpu.tsv
  matrix_benchmark_analysis_mrvi_cpu.tsv
  matrix_batch_effect_uncorrected_mrvi.tsv
)
MANIFEST_COUNTS=(2 1 2 1 15 1 3 6 1)
MANIFEST_COLUMNS=(3 3 3 3 4 3 4 4 3)
for manifest_idx in "${!MANIFEST_NAMES[@]}"; do
  manifest="${RUN_ROOT}/manifests/${MANIFEST_NAMES[${manifest_idx}]}"
  [[ -s "${manifest}" ]]
  [[ "$(wc -l < "${manifest}" | tr -d '[:space:]')" == "${MANIFEST_COUNTS[${manifest_idx}]}" ]]
  awk -F '\t' -v expected="${MANIFEST_COLUMNS[${manifest_idx}]}" \
    'NF != expected {exit 1}' "${manifest}"
done
MOFA_BENCH="$(sed -n '5p' "${CALLS}")"
MOFA_BATCH="$(sed -n '7p' "${CALLS}")"
case "${MOFA_BENCH}" in *"--dependency=afterok:780002"*) ;; *) echo "benchmark MOFA dependency missing" >&2; exit 1 ;; esac
case "${MOFA_BATCH}" in *"--dependency=afterok:780004"*) ;; *) echo "batch-view MOFA dependency missing" >&2; exit 1 ;; esac
for line in 9 11 13 15 17; do
  if sed -n "${line}p" "${CALLS}" | grep -q -- '--dependency=afterok:'; then
    echo "independent Stage 5 method carried pseudobulk dependency" >&2
    exit 1
  fi
done
if sed -n '1,18p' "${CALLS}" | grep -q -- '--wait'; then
  echo "an independent array waited before aggregate submission" >&2
  exit 1
fi
case "$(sed -n '19p' "${CALLS}")" in *"matrix_gate.sh"*) ;; *) echo "aggregate gate was not last" >&2; exit 1 ;; esac
echo "benchmark dependency edges: OK"
