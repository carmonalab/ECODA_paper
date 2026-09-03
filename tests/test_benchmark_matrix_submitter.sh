#!/bin/bash
# Deterministic Stage 5 matrix contract: grouped view/label arrays, one
# aggregate gate, run-owned manifests/logs, and no hidden dependency waits.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
unset HPC_SCRATCH_DIR
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-matrix-stage.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '70000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
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
WORKER_MANIFEST="${TMP_DIR}/worker.tsv"
printf 'Adams\tbenchmark_analysis\tmrvi\n' > "${WORKER_MANIFEST}"
WORKER_SCRATCH="${TMP_DIR}/worker-scratch"
WORKER_ROOT="${TMP_DIR}/worker-root"
mkdir -p "${WORKER_SCRATCH}" "${WORKER_ROOT}"
PY_CALL_LOG="${PY_CALL_LOG}" \
HOME="${TMP_DIR}/home" PATH="/usr/bin:/bin" TMPDIR="${TMP_DIR}" \
HPC_SCRATCH_DIR="${WORKER_SCRATCH}" LOGS_DIR="${TMP_DIR}/worker-logs" \
ECODA_RUNTIME_MODE=host ECODA_RUNTIME_IN_CONTAINER=1 \
ECODA_RUNTIME_PREFIX="${FAKE_PREFIX}" ECODA_APPTAINER_NV=1 \
ECODA_RUNTIME_PROFILE=stage5 METHOD_GPU_POLICY=default METHOD=mrvi \
ANALYSIS_MANIFEST="${WORKER_MANIFEST}" ANALYSIS_ROOT="${WORKER_ROOT}" \
EXECUTION_LOG_DIR="${WORKER_ROOT}/embeddings" \
SLURM_ARRAY_TASK_ID=1 SLURM_ARRAY_JOB_ID=91001 \
  bash "${ROOT}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
case "$(cat "${PY_CALL_LOG}")" in
  *"--device cuda"*) ;;
  *) echo "GPU worker did not force --device cuda" >&2; exit 1 ;;
esac


: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  ECODA_RUNTIME_MODE=apptainer BENCHMARK_MATRIX_TEST=1 \
  bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
  --datasets Adams --methods mrvi; then
  echo "apptainer mode accepted a missing image" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]

BATCH_SELECTION="${TMP_DIR}/batch-matrix.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nKidney_KPMP\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n' > "${BATCH_SELECTION}"
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
: > "${CAPTURE}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams --methods mrvi >/dev/null
[[ ! -s "${CAPTURE}" ]]
printf 'MD5=00000000000000000000000000000000\n' > "${HPC_ROOT}/benchmark/embeddings/Adams_hvg2000_mrvi_dists.feather.md5"
: > "${CAPTURE}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  BENCHMARK_MATRIX_TEST=1 bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
    --datasets Adams --methods mrvi >/dev/null
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == 3 ]]
[[ ! -e "${HPC_ROOT}/benchmark/embeddings/Adams_hvg2000_mrvi_dists.feather.md5" ]]
echo "benchmark matrix submitter: OK"
