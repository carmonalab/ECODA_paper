#!/bin/bash
# Deterministic Stage 5 matrix contract: grouped view/label arrays, one
# aggregate gate, run-owned manifests/logs, and no hidden dependency waits.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
unset HPC_SCRATCH_DIR
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-matrix-stage.XXXXXX")"
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
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${HOST_PREFIX}/bin" \
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
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf 'ENV ECODA_SOURCE_ROOT=%s ECODA_SOURCE_MANIFEST=%s ECODA_RUNTIME_IMAGE=%s ECODA_RUNTIME_MANIFEST=%s ECODA_RUN_ROOT=%s ECODA_RUN_ID=%s ARGS %s\n' \
  "${ECODA_SOURCE_ROOT:-}" "${ECODA_SOURCE_MANIFEST:-}" \
  "${ECODA_RUNTIME_IMAGE:-}" "${ECODA_RUNTIME_MANIFEST:-}" \
  "${ECODA_RUN_ROOT:-}" "${ECODA_RUN_ID:-}" "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '70000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
export HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_SCRATCH_ROOT="${HPC_ROOT}"
export ECODA_LOGS_DIR="${TEST_LOGS}" TMPDIR="${TEST_TMP}"
export ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}"
export ECODA_HOST_PYTHON_BIN="${HOST_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${HOST_PREFIX}/bin/Rscript --vanilla"
export ECODA_AUX_ROOT="${SOURCE_TREE}/aux"
export ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" ECODA_RUNTIME_PROFILE=stage5
export ECODA_APPTAINER_NV=0 APPTAINER_BIN="${TMP_DIR}/bin/apptainer"
OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams,Bassez,Kfoury,Kim --methods mrvi,gloscope
)"
RUN_ID="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^BENCHMARK_RUN_ID=//p')"
[[ -n "${RUN_ID}" ]]
RUN_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}"
[[ -s "${RUN_ROOT}/manifests/selection.tsv" ]]
MANIFEST_NAMES=(
  matrix_benchmark_analysis_mrvi_default_gpu.tsv
  matrix_benchmark_analysis_mrvi_cpu.tsv
  matrix_benchmark_analysis_gloscope_cpu.tsv
)
MANIFEST_COUNTS=(4 8 20)
for manifest_idx in 0 1 2; do
  manifest_name="${MANIFEST_NAMES[${manifest_idx}]}"
  manifest="${RUN_ROOT}/manifests/${manifest_name}"
  [[ "$(wc -l < "${manifest}" | tr -d '[:space:]')" == "${MANIFEST_COUNTS[${manifest_idx}]}" ]]
  while IFS=$'\t' read -r ds view label combo; do
    [[ -n "${ds}" && "${view}" == benchmark_analysis && -n "${combo}" ]]
    case "${label}" in mrvi|gloscope) ;; *) exit 1 ;; esac
  done < "${manifest}"
done
[[ "$(wc -l < "${RUN_ROOT}/manifests/scheduler_ids.tsv" | tr -d '[:space:]')" == 7 ]]
CALLS="$(cat "${CAPTURE}")"
[[ "$(printf '%s\n' "${CALLS}" | wc -l | tr -d '[:space:]')" == 7 ]]
FIRST_ARRAY_CALL="$(printf '%s\n' "${CALLS}" | grep -- '--array=' | sed -n '1p')"
case "${FIRST_ARRAY_CALL}" in
  *"--export=ALL,"*"METHOD=mrvi"*"ANALYSIS_MANIFEST="*) ;;
  *) echo "Stage 5 array omitted worker environment export" >&2; exit 1 ;;
esac
case "${FIRST_ARRAY_CALL}" in
  *"${SOURCE_TREE}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"*) ;;
  *) echo "Stage 5 array omitted immutable worker script" >&2; exit 1 ;;
esac
for required_path in \
  "${SOURCE_TREE}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh" \
  "${SOURCE_TREE}/src/5_run_benchmark_methods/matrix_watchdog.sh" \
  "${SOURCE_TREE}/src/5_run_benchmark_methods/matrix_gate.sh"; do
  case "${CALLS}" in
    *"${required_path}"*) ;;
    *) echo "scheduler command escaped immutable source root: ${required_path}" >&2; exit 1 ;;
  esac
done
for required_identity in \
  "ECODA_SOURCE_ROOT=${SOURCE_TREE}" \
  "ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}" \
  "ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}" \
  "ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}" \
  "ECODA_RUN_ROOT=${RUN_ROOT}"; do
  case "${CALLS}" in
    *"${required_identity}"*) ;;
    *) echo "scheduler command omitted immutable run identity: ${required_identity}" >&2; exit 1 ;;
  esac
done
GATE_CALL="$(printf '%s\n' "${CALLS}" | grep 'matrix_gate.sh' | tail -1)"
case "${GATE_CALL}" in
  *"${SOURCE_TREE}/src/5_run_benchmark_methods/matrix_gate.sh"*) ;;
  *) echo "aggregate gate was not rooted in the immutable source snapshot" >&2; exit 1 ;;
esac
for gate_identity in \
  "ECODA_SOURCE_ROOT=${SOURCE_TREE}" \
  "ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}" \
  "ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}" \
  "ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}" \
  "ECODA_RUN_ROOT=${RUN_ROOT}"; do
  case "${GATE_CALL}" in
    *"${gate_identity}"*) ;;
    *) echo "aggregate gate omitted immutable run identity: ${gate_identity}" >&2; exit 1 ;;
  esac
done
[[ "$(grep '^SOURCE_ROOT=' "${RUN_ROOT}/manifests/source.manifest")" == "SOURCE_ROOT=${SOURCE_TREE}" ]]
[[ "$(grep '^RUNTIME_IMAGE=' "${RUN_ROOT}/manifests/runtime.identity")" == "RUNTIME_IMAGE=${RUNTIME_IMAGE}" ]]
[[ "$(grep '^RUNTIME_MANIFEST=' "${RUN_ROOT}/manifests/runtime.identity")" == "RUNTIME_MANIFEST=${RUNTIME_MANIFEST}" ]]
grep -q "^SOURCE_SOURCE_ROOT=${SOURCE_TREE}$" "${RUN_ROOT}/metadata"
grep -q "^RUNTIME_IMAGE=${RUNTIME_IMAGE}$" "${RUN_ROOT}/metadata"
grep -q "^RUNTIME_MANIFEST=${RUNTIME_MANIFEST}$" "${RUN_ROOT}/metadata"
PREFLIGHT_LINE="$(awk '/^stage5_compute_h5ad_preflight \|\|/{print NR; exit}' "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh")"
IDENTITY_LINE="$(awk '/^stage5_prepare_source_identity \|\|/{print NR; exit}' "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh")"
[[ -n "${PREFLIGHT_LINE}" && -n "${IDENTITY_LINE}" && ${PREFLIGHT_LINE} -lt ${IDENTITY_LINE} ]] ||
  { echo "Stage 5 source identity was trusted before strict H5AD preflight" >&2; exit 1; }
case "${CALLS}" in *"matrix_gate.sh"*) ;; *) echo "aggregate gate was not submitted" >&2; exit 1 ;; esac
case "${CALLS}" in *"--dependency=afterany:"*) ;; *) echo "aggregate gate dependency missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--array=1-4"*) ;; *) echo "default MrVI shard array missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--array=1-8"*) ;; *) echo "CPU MrVI shard array missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--array=1-20"*) ;; *) echo "GloScope shard array missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--partition=${SLURM_PARTITION_BENCHMARK_GPU}"*"--constraint=${BENCHMARK_GPU_CONSTRAINT}"*) ;; *) echo "default GPU method resource class missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--partition=${SLURM_PARTITION_BENCHMARK_CPU}"*"METHOD=mrvi"*) ;; *) echo "CPU MrVI shard resource class missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--time=${BENCHMARK_GPU_DEFAULT_TIME_LIMIT}"*) ;; *) echo "default GPU method time limit missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--partition=${SLURM_PARTITION_BENCHMARK_CPU}"*) ;; *) echo "CPU method resource class missing" >&2; exit 1 ;; esac
while IFS= read -r call; do
  case "${call}" in
    *matrix_watchdog*)
      case "${call}" in
        *"--partition=${SLURM_PARTITION_BENCHMARK_CPU}"*) ;;
        *) echo "matrix watchdog was scheduled outside CPU partition" >&2; exit 1 ;;
      esac
      case "${call}" in
        *"--dependency=afterany:"*) ;;
        *) echo "matrix watchdog dependency on its array is missing" >&2; exit 1 ;;
      esac
      ;;
  esac
done <<< "${CALLS}"
case "${CALLS}" in *"_ecoda_runs/${RUN_ID}/logs/"*) ;; *) echo "scheduler logs escaped run root" >&2; exit 1 ;; esac

: > "${CAPTURE}"
SYNC_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --sync-only "${RUN_ID}"
)"
case "${SYNC_OUTPUT}" in *"BENCHMARK_RUN_ID=${RUN_ID}"*) ;; *) echo "sync-only recovery failed" >&2; exit 1 ;; esac
[[ ! -s "${CAPTURE}" ]]
: > "${CAPTURE}"
ANY_OUTPUT="$(
  BENCHMARK_GPU_ANY_VRAM_PER_GPU=80G \
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams --methods scpoli --gpu-policy any
)"
ANY_CALLS="$(cat "${CAPTURE}")"
case "${ANY_CALLS}" in *"--partition=${BENCHMARK_GPU_ANY_PARTITION}"*) ;; *) echo "any-GPU policy partition missing" >&2; exit 1 ;; esac
case "${ANY_CALLS}" in *"--constraint=${BENCHMARK_GPU_CONSTRAINT}"*) echo "any-GPU policy retained H200 constraint" >&2; exit 1 ;; esac
case "${ANY_CALLS}" in *"--time=${BENCHMARK_GPU_ANY_TIME_LIMIT}"*) ;; *) echo "any-GPU policy time limit missing" >&2; exit 1 ;; esac
case "${ANY_CALLS}" in *"ECODA_APPTAINER_NV=1"*) ;; *) echo "any-GPU policy runtime NV flag missing" >&2; exit 1 ;; esac
case "${ANY_CALLS}" in *"--gres=gpu:${BENCHMARK_GPU_COUNT},VramPerGpu:80G"*) ;; *) echo "any-GPU VRAM request missing" >&2; exit 1 ;; esac
case "${ANY_CALLS}" in *"--gpus="*) echo "any-GPU VRAM request retained --gpus" >&2; exit 1 ;; esac
: > "${CAPTURE}"
PILOTGM_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams --methods pilotgm
)"
PILOTGM_CALLS="$(cat "${CAPTURE}")"
case "${PILOTGM_CALLS}" in *"--partition=${SLURM_PARTITION_BENCHMARK_CPU}"*) ;; *) echo "pilotgm was not moved to CPU resources" >&2; exit 1 ;; esac
case "${PILOTGM_CALLS}" in *"--gpus="*) echo "pilotgm retained a GPU flag" >&2; exit 1 ;; esac
case "${PILOTGM_CALLS}" in *"ECODA_APPTAINER_NV=0"*) ;; *) echo "pilotgm runtime NV flag missing" >&2; exit 1 ;; esac
case "${PILOTGM_CALLS}" in *"--array=1-1"*) ;; *) echo "pilotgm parameter screening was not reduced to default-only" >&2; exit 1 ;; esac

PY_CALL_LOG="${TMP_DIR}/python.call"
FAKE_PREFIX="${TMP_DIR}/worker-prefix"
mkdir -p "${FAKE_PREFIX}/bin" "${FAKE_PREFIX}/lib/R"
cat > "${FAKE_PREFIX}/bin/python" <<'STUB'
#!/bin/bash
printf '%s\n' "$*" > "${PY_CALL_LOG:?}"
exit 0
STUB
cat > "${FAKE_PREFIX}/bin/Rscript" <<'STUB'
#!/bin/bash
exit 0
STUB
chmod +x "${FAKE_PREFIX}/bin/python" "${FAKE_PREFIX}/bin/Rscript"
WORKER_SCRATCH="${TMP_DIR}/worker-scratch"
WORKER_RUN_ID="stage5-direct-worker"
WORKER_RUN_ROOT="${WORKER_SCRATCH}/_ecoda_runs/${WORKER_RUN_ID}"
WORKER_RUNTIME_MANIFEST_SHA="$(sha256_file "${RUNTIME_MANIFEST}")"
WORKER_RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
WORKER_RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
WORKER_ROOT="${WORKER_RUN_ROOT}/analysis"
mkdir -p "${WORKER_RUN_ROOT}/manifests" "${WORKER_ROOT}"
cp "${SOURCE_MANIFEST}" "${WORKER_RUN_ROOT}/manifests/source.manifest"
printf 'RUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA}" \
  "${WORKER_RUNTIME_MANIFEST_SHA}" "${WORKER_RUNTIME_IMAGE_SIZE}" \
  "${WORKER_RUNTIME_MANIFEST_SIZE}" "${SOURCE_TOML_SHA}" "${SOURCE_LOCK_SHA}" \
  > "${WORKER_RUN_ROOT}/manifests/runtime.identity"
chmod 444 "${WORKER_RUN_ROOT}/manifests/source.manifest" \
  "${WORKER_RUN_ROOT}/manifests/runtime.identity"
WORKER_MANIFEST="${WORKER_RUN_ROOT}/manifests/worker.tsv"
printf 'Adams\tbenchmark_analysis\tmrvi\n' > "${WORKER_MANIFEST}"
PY_CALL_LOG="${PY_CALL_LOG}" \
HOME="${TMP_DIR}/home" PATH="/usr/bin:/bin" TMPDIR="${TMP_DIR}" \
HPC_SCRATCH_DIR="${WORKER_SCRATCH}" LOGS_DIR="${TMP_DIR}/worker-logs" \
ECODA_RUNTIME_MODE=host ECODA_RUNTIME_IN_CONTAINER=1 \
ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_TREE}/aux" \
ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" \
ECODA_RUN_ROOT="${WORKER_RUN_ROOT}" ECODA_RUN_ID="${WORKER_RUN_ID}" \
ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" \
ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
ECODA_RUNTIME_PREFIX="${FAKE_PREFIX}" ECODA_APPTAINER_NV=1 \
ECODA_RUNTIME_PROFILE=stage5 METHOD_GPU_POLICY=default METHOD=mrvi \
ANALYSIS_MANIFEST="${WORKER_MANIFEST}" ANALYSIS_ROOT="${WORKER_ROOT}" \
EXECUTION_LOG_DIR="${WORKER_ROOT}/embeddings" \
SLURM_ARRAY_TASK_ID=1 SLURM_ARRAY_JOB_ID=91001 \
  bash "${SOURCE_TREE}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
case "$(cat "${PY_CALL_LOG}")" in
  *"--device cuda"*) ;;
  *) echo "GPU worker did not force --device cuda" >&2; exit 1 ;;
esac


: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  ECODA_RUNTIME_IMAGE="${TMP_DIR}/missing-runtime.sif" \
  ECODA_RUNTIME_MANIFEST="${TMP_DIR}/missing-runtime.sif.manifest" \
  BENCHMARK_MATRIX_TEST=1 \
  bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
  --datasets Adams --methods mrvi; then
  echo "apptainer mode accepted a missing image" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
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
  echo "mutable matrix worker script escaped source containment" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  HPC_SCRATCH_DIR="${HPC_ROOT}" ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" \
  bash -c '
    set -e
    source "$1/src/slurm_config.sh"
    source "$1/src/utils/bash/ecoda_run_common.sh"
    ecoda_require_source_script_path "$2" "$3"
  ' _ "${ROOT}" \
  "${SOURCE_TREE}/src/5_run_benchmark_methods/matrix_gate.sh" "${SOURCE_TREE}"
BATCH_SELECTION="${TMP_DIR}/batch-matrix.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nKidney_KPMP_full\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n' > "${BATCH_SELECTION}"
: > "${CAPTURE}"
PASS_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${BATCH_SELECTION}" --exact-batch-selection --pass uncorrected \
    --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
)"
case "${PASS_OUTPUT}" in *"BATCH_EFFECT_RUN_ID="*) ;; *) echo "batch run marker missing" >&2; exit 1 ;; esac
if printf '%s\n' "${PASS_OUTPUT}" | grep -q '^BENCHMARK_'; then
  echo "ordinary BENCHMARK marker leaked into batch mode" >&2
  exit 1
fi
PASS_CALLS="$(cat "${CAPTURE}")"
case "${PASS_CALLS}" in *pilotgm*) echo "batch PILOT-GM-VAE leaked into scheduler calls" >&2; exit 1 ;; esac
case "${PASS_CALLS}" in *"ANALYSIS_MANIFEST="*) ;; *) echo "batch ANALYSIS_MANIFEST export missing" >&2; exit 1 ;; esac
case "${PASS_CALLS}" in *"BENCHMARK_MANIFEST="*) echo "batch BENCHMARK_MANIFEST export leaked" >&2; exit 1 ;; esac
case "${PASS_CALLS}" in *"/batch_effect/uncorrected"*) ;; *) echo "batch analysis root was not pass-scoped" >&2; exit 1 ;; esac
case "${PASS_CALLS}" in *"--array=1-12"*) ;; *) echo "batch matrix arrays did not use twelve rows" >&2; exit 1 ;; esac
case "${PASS_CALLS}" in *"--partition=${BENCHMARK_GPU_ANY_PARTITION}"*) ;; *) echo "batch GPU method did not use any-GPU partition" >&2; exit 1 ;; esac
case "${PASS_CALLS}" in *"--constraint=${BENCHMARK_GPU_CONSTRAINT}"*) echo "batch GPU method retained H200 constraint" >&2; exit 1 ;; esac
case "${PASS_CALLS}" in *"--time=${BENCHMARK_GPU_ANY_TIME_LIMIT}"*) ;; *) echo "batch GPU method time limit missing" >&2; exit 1 ;; esac
rm -rf "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners"
TARGET_SELECTION="${TMP_DIR}/batch-target.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n' > "${TARGET_SELECTION}"
: > "${CAPTURE}"
TARGET_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${TARGET_SELECTION}" --pass uncorrected \
    --target-methods gloscope --partition shared-bigmem --mem 500G --max-mem 500G
)"
TARGET_RUN_ID="$(printf '%s\n' "${TARGET_OUTPUT}" | sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
[[ -n "${TARGET_RUN_ID}" ]]
TARGET_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${TARGET_RUN_ID}"
[[ "$(cat "${TARGET_ROOT}/manifests/pending_selection.tsv")" == $'Alzheimer\tbatch_effect_uncorrected\tgloscope' ]]
TARGET_MATRIX="${TARGET_ROOT}/manifests/matrix_batch_effect_uncorrected_gloscope.tsv"
[[ -s "${TARGET_MATRIX}" ]]
[[ "$(cat "${TARGET_MATRIX}")" == $'Alzheimer\tbatch_effect_uncorrected\tgloscope' ]]
[[ ! -e "${TARGET_ROOT}/manifests/matrix_batch_effect_uncorrected_prepare_pseudobulk.tsv" ]]
TARGET_CALLS="$(cat "${CAPTURE}")"
case "${TARGET_CALLS}" in *"--partition=shared-bigmem"*"METHOD=gloscope"*) ;; *) echo "targeted GloScope did not use explicit bigmem partition" >&2; exit 1 ;; esac
case "${TARGET_CALLS}" in *"prepare_pseudobulk"*) echo "targeted GloScope scheduled pseudobulk preparation" >&2; exit 1 ;; esac
case "${TARGET_CALLS}" in *"--array=1-1"*) ;; *) echo "targeted GloScope array did not contain one selected row" >&2; exit 1 ;; esac
rm -rf "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners"
# Only the changed high-resolution consumers and the new current Kidney key
# remain absent.  Valid unaffected batch rows must stay outside this targeted
# recomputation selection.
for ds in Alzheimer Breast_cancer Covid19_PBMC Myocardial_infarction Diabetes Lung Joanito Stephenson CombinedPBMC; do
  path="${HPC_ROOT}/batch_effect/uncorrected/embeddings/${ds}_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"
  printf 'valid-%s\n' "${ds}" > "${path}"
  digest="$(md5sum "${path}" | cut -d' ' -f1)"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
done
: > "${CAPTURE}"
TARGETED_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${BATCH_SELECTION}" --pass uncorrected --target-methods mrvi
)"
TARGETED_RUN_ID="$(printf '%s\n' "${TARGETED_OUTPUT}" | sed -n 's/^BATCH_EFFECT_RUN_ID=//p')"
[[ -n "${TARGETED_RUN_ID}" ]]
TARGETED_ROOT="${HPC_ROOT}/_ecoda_runs/${TARGETED_RUN_ID}"
[[ "$(cat "${TARGETED_ROOT}/manifests/pending_selection.tsv")" == $'Kidney_KPMP_full\tbatch_effect_uncorrected\tmrvi\nLupus_PBMC\tbatch_effect_uncorrected\tmrvi\nParkinson\tbatch_effect_uncorrected\tmrvi' ]]
for unaffected in Alzheimer Breast_cancer Covid19_PBMC Myocardial_infarction Diabetes Lung Joanito Stephenson CombinedPBMC; do
  ! grep -q "^${unaffected}" "${TARGETED_ROOT}/manifests/pending_selection.tsv"
done
TARGETED_MATRIX="${TARGETED_ROOT}/manifests/matrix_batch_effect_uncorrected_mrvi.tsv"
[[ "$(cat "${TARGETED_MATRIX}")" == $'Kidney_KPMP_full\tbatch_effect_uncorrected\tmrvi\nLupus_PBMC\tbatch_effect_uncorrected\tmrvi\nParkinson\tbatch_effect_uncorrected\tmrvi' ]]
TARGETED_CALLS="$(cat "${CAPTURE}")"
[[ "$(printf '%s\n' "${TARGETED_CALLS}" | wc -l | tr -d '[:space:]')" == 3 ]]
case "${TARGETED_CALLS}" in *"--array=1-3"*) ;; *) echo "targeted changed-input rows did not form a three-row array" >&2; exit 1 ;; esac
case "${TARGETED_CALLS}" in *"METHOD=mrvi"*) ;; *) echo "targeted changed-input method was not submitted" >&2; exit 1 ;; esac
case "${TARGETED_CALLS}" in *"METHOD=prepare_pseudobulk"*) echo "unselected batch method leaked into targeted rerun" >&2; exit 1 ;; esac
rm -rf "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners"

BAD_SCOPE_SELECTION="${TMP_DIR}/batch-wrong-scope.tsv"
sed '1s/batch_effect_uncorrected$/wrong_scope/' "${BATCH_SELECTION}" > "${BAD_SCOPE_SELECTION}"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${BAD_SCOPE_SELECTION}" --pass uncorrected \
    --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot; then
  echo "wrong batch scope was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
BAD_EXACT_SELECTION="${TMP_DIR}/batch-malformed-exact.tsv"
sed '1s/^Alzheimer/Breast_cancer/' "${BATCH_SELECTION}" > "${BAD_EXACT_SELECTION}"
RUNS_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${BAD_EXACT_SELECTION}" --exact-batch-selection \
    --pass uncorrected \
    --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot; then
  echo "malformed exact batch selection was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]
BAD_LEGACY_SELECTION="${TMP_DIR}/batch-legacy-key.tsv"
sed 's/^Kidney_KPMP_full/Kidney_KPMP/' "${BATCH_SELECTION}" > "${BAD_LEGACY_SELECTION}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${BAD_LEGACY_SELECTION}" --exact-batch-selection \
    --pass uncorrected \
    --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot; then
  echo "legacy Kidney_KPMP exact batch key was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]


rm -rf "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners"
HPC_ROOT="${TMP_DIR}/home/scratch/ECODA_paper"
mkdir -p "${HPC_ROOT}/benchmark/embeddings"
for n in 1000 2000 3000; do
  selected="${HPC_ROOT}/benchmark/embeddings/Adams_hvg${n}_mrvi_dists.feather"
  printf 'selected\n' > "${selected}"
  digest="$(md5sum "${selected}" | cut -d' ' -f1)"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${selected}" | tr -d '[:space:]')" "${selected}" > "${selected}.md5"
done
unrelated="${HPC_ROOT}/benchmark/embeddings/Unrelated_hvg1000_mrvi_dists.feather"
printf 'unrelated\n' > "${unrelated}"
digest="$(md5sum "${unrelated}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${unrelated}" | tr -d '[:space:]')" "${unrelated}" > "${unrelated}.md5"
: > "${CAPTURE}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams --methods mrvi --force >/dev/null
for n in 1000 2000 3000; do
  [[ ! -e "${HPC_ROOT}/benchmark/embeddings/Adams_hvg${n}_mrvi_dists.feather.md5" ]]
done
[[ -s "${unrelated}.md5" ]]
rm -rf "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners"
for n in 1000 2000 3000; do
  selected="${HPC_ROOT}/benchmark/embeddings/Adams_hvg${n}_mrvi_dists.feather"
  pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"s1":[1.0],"s2":[0.0]},index=["s1"]).to_feather(sys.argv[1])' "${selected}"
  digest="$(md5sum "${selected}" | cut -d' ' -f1)"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${selected}" | tr -d '[:space:]')" "${selected}" > "${selected}.md5"
done
ARTIFACT_PATH="${HPC_ROOT}/benchmark/embeddings/Adams_hvg1000_mrvi_dists.feather"
ARTIFACT_RECORD="$(
  HOME="${TMP_DIR}/home" HPC_SCRATCH_DIR="${HPC_ROOT}" \
  ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" ECODA_RUN_ROOT="${RUN_ROOT}" \
  ECODA_RUN_ID="${RUN_ID}" PATH="${TMP_DIR}/bin:${PATH}" \
  bash -c '
    set -e
    source "$1/src/slurm_config.sh"
    source "$1/src/utils/bash/ecoda_run_common.sh"
    ecoda_write_artifact_record "$2" stage5_mrvi "$3"
  ' _ "${ROOT}" "${ARTIFACT_PATH}" "${RUN_ID}"
)"
[[ -s "${ARTIFACT_RECORD}" ]]
grep -q "^PATH=${ARTIFACT_PATH}$" "${ARTIFACT_RECORD}"
grep -q '^PRODUCER=stage5_mrvi$' "${ARTIFACT_RECORD}"
grep -q "^RUN_ID=${RUN_ID}$" "${ARTIFACT_RECORD}"
HOME="${TMP_DIR}/home" HPC_SCRATCH_DIR="${HPC_ROOT}" \
  ECODA_HOST_ENV_PREFIX="${HOST_PREFIX}" ECODA_RUN_ROOT="${RUN_ROOT}" \
  ECODA_RUN_ID="${RUN_ID}" PATH="${TMP_DIR}/bin:${PATH}" \
  bash -c '
    set -e
    source "$1/src/slurm_config.sh"
    source "$1/src/utils/bash/ecoda_run_common.sh"
    ecoda_validate_artifact_record "$2" stage5_mrvi "$3" >/dev/null
  ' _ "${ROOT}" "${ARTIFACT_PATH}" "${RUN_ID}"

# A terminally failed/stale gate is evidence only.  It must not add rows to a
# later explicit selection.  prepare_pseudobulk is a registered non-default
# method and is submitted only because it is explicitly named here, while
# the valid baseline MrVI rows remain skipped.
STALE_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/stale_failed_gate"
mkdir -p "${STALE_RUN_ROOT}/manifests" "${STALE_RUN_ROOT}/status"
printf 'STATE=FAIL\nREASON=stale historical gate\n' > "${STALE_RUN_ROOT}/status/terminal"
printf 'Bassez\tbenchmark_analysis\tprepare_pseudobulk\n' > "${STALE_RUN_ROOT}/manifests/selection.tsv"
POST_BASELINE_SELECTION="${TMP_DIR}/ordinary-target.tsv"
printf 'Adams\tbenchmark_analysis\tmrvi\nAdams\tbenchmark_analysis\tprepare_pseudobulk\n' > \
  "${POST_BASELINE_SELECTION}"
: > "${CAPTURE}"
POST_BASELINE_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --selection-file "${POST_BASELINE_SELECTION}"
)"
POST_BASELINE_RUN_ID="$(printf '%s\n' "${POST_BASELINE_OUTPUT}" | sed -n 's/^BENCHMARK_RUN_ID=//p')"
[[ -n "${POST_BASELINE_RUN_ID}" ]]
POST_BASELINE_ROOT="${HPC_ROOT}/_ecoda_runs/${POST_BASELINE_RUN_ID}"
[[ "$(cat "${POST_BASELINE_ROOT}/manifests/pending_selection.tsv")" == $'Adams\tbenchmark_analysis\tprepare_pseudobulk' ]]
POST_BASELINE_CALLS="$(cat "${CAPTURE}")"
[[ "$(printf '%s\n' "${POST_BASELINE_CALLS}" | wc -l | tr -d '[:space:]')" == 3 ]]
case "${POST_BASELINE_CALLS}" in *"METHOD=mrvi"*) echo "valid baseline MrVI row was resubmitted" >&2; exit 1 ;; esac
case "${POST_BASELINE_CALLS}" in *"METHOD=prepare_pseudobulk"*) ;; *) echo "explicit non-default row was not submitted" >&2; exit 1 ;; esac
case "${POST_BASELINE_CALLS}" in *"FORCE_BENCHMARK=1"*) echo "implicit force leaked into explicit non-default selection" >&2; exit 1 ;; esac
case "$(cat "${POST_BASELINE_ROOT}/manifests/pending_selection.tsv")" in *Bassez*) echo "stale gate added a selection row" >&2; exit 1 ;; esac
: > "${CAPTURE}"
NOOP_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams --methods mrvi
)"
case "${NOOP_OUTPUT}" in *"NOOP_VALIDATED=1"*) ;; *) echo "valid baseline row did not take the no-op path" >&2; exit 1 ;; esac
[[ "$(grep '^STATE=' "${HPC_ROOT}/_ecoda_runs/$(printf '%s\n' "${NOOP_OUTPUT}" | sed -n 's/^BENCHMARK_RUN_ID=//p')/status/report")" == "STATE=NOOP_VALIDATED" ]]
[[ ! -s "${CAPTURE}" ]]
printf 'MD5=00000000000000000000000000000000\n' > "${HPC_ROOT}/benchmark/embeddings/Adams_hvg2000_mrvi_dists.feather.md5"
: > "${CAPTURE}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams --methods mrvi >/dev/null
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == 3 ]]
[[ ! -e "${HPC_ROOT}/benchmark/embeddings/Adams_hvg2000_mrvi_dists.feather.md5" ]]
HOSTILE_OUTPUT="${TMP_DIR}/outside-stage5-output"
mkdir -p "${HOSTILE_OUTPUT}" "${HPC_ROOT}/Adams"
ln -s "${HOSTILE_OUTPUT}" "${HPC_ROOT}/Adams/output"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
  --datasets Adams --methods mrvi >/dev/null 2>&1; then
  echo "Stage 5 accepted an output-root symlink escape" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
rm -f "${HPC_ROOT}/Adams/output"
HOSTILE_SELECTION="${TMP_DIR}/hostile-dataset.tsv"
printf '../Adams\tbenchmark_analysis\tmrvi\n' > "${HOSTILE_SELECTION}"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
  --selection-file "${HOSTILE_SELECTION}" >/dev/null 2>&1; then
  echo "Stage 5 accepted a traversal dataset selection" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
MALFORMED_SELECTION="${TMP_DIR}/malformed-selection.tsv"
printf 'Adams\tbenchmark_analysis\tmrvi\textra\n' > "${MALFORMED_SELECTION}"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
  --selection-file "${MALFORMED_SELECTION}" >/dev/null 2>&1; then
  echo "Stage 5 accepted a four-column selection row" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
echo "benchmark matrix submitter: OK"
