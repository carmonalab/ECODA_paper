#!/bin/bash
set -euo pipefail

SUBMIT_MODE="${1:-all}"
ARRAY_OVERRIDE="${2:-}"
case "${SUBMIT_MODE}" in
  all|cpu|gpu) ;;
  *)
    echo "Usage: $0 [all|cpu|gpu] [array_spec]" >&2
    exit 1
    ;;
esac
if [[ "${SUBMIT_MODE}" == "all" && -n "${ARRAY_OVERRIDE}" ]]; then
  echo "ERROR: array_spec override requires mode 'cpu' or 'gpu'" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
source "${PROJECT_ROOT}/src/slurm_config.sh"
cd "${PROJECT_ROOT}"

RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/sf15_classic_profile_20261126"
MANIFEST_ALL="${SCRIPT_DIR}/manifest_16rows.tsv"
MANIFEST_CPU="${RUN_ROOT}/manifests/manifest_cpu_11rows.tsv"
MANIFEST_GPU="${RUN_ROOT}/manifests/manifest_gpu_5rows.tsv"
WORKER_SH="${SCRIPT_DIR}/task_worker.sh"

mkdir -p "${RUN_ROOT}/manifests" "${RUN_ROOT}/rows" "${RUN_ROOT}/outputs" "${RUN_ROOT}/logs"

awk -F'\t' 'NR > 1 && $4 == "cpu" { print $0 }' "${MANIFEST_ALL}" > "${MANIFEST_CPU}"
awk -F'\t' 'NR > 1 && $4 == "gpu" { print $0 }' "${MANIFEST_ALL}" > "${MANIFEST_GPU}"

CPU_COUNT="$(wc -l < "${MANIFEST_CPU}" | tr -d ' ')"
GPU_COUNT="$(wc -l < "${MANIFEST_GPU}" | tr -d ' ')"
if [[ "${CPU_COUNT}" != "11" || "${GPU_COUNT}" != "5" ]]; then
  echo "ERROR: expected 11 CPU and 5 GPU rows, got cpu=${CPU_COUNT} gpu=${GPU_COUNT}" >&2
  exit 1
fi

echo "RUN_ROOT=${RUN_ROOT}"
echo "SUBMIT_MODE=${SUBMIT_MODE}"

if [[ "${SUBMIT_MODE}" == "all" || "${SUBMIT_MODE}" == "cpu" ]]; then
  CPU_ARRAY_SPEC="${ARRAY_OVERRIDE:-1-${CPU_COUNT}}"
  CPU_JOB_ID="$(
    sbatch --parsable \
      --job-name=sf15_prof_cpu \
      --partition="${SLURM_PARTITION_BENCHMARK_CPU}" \
      --constraint="${BENCHMARK_CPU_CONSTRAINT}" \
      --nodes=1 \
      --ntasks=1 \
      --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}" \
      --mem="${BENCHMARK_MEM}" \
      --time=04:00:00 \
      --array="${CPU_ARRAY_SPEC}" \
      --output="${RUN_ROOT}/logs/cpu_%A_%a.out" \
      --error="${RUN_ROOT}/logs/cpu_%A_%a.err" \
      --export=ALL,ECODA_PROJECT_ROOT="${PROJECT_ROOT}",SF15_RUN_ROOT="${RUN_ROOT}",SF15_MANIFEST_FILE="${MANIFEST_CPU}" \
      "${WORKER_SH}"
  )"
  echo "CPU_JOB_ID=${CPU_JOB_ID} (array=${CPU_ARRAY_SPEC})"
fi

if [[ "${SUBMIT_MODE}" == "all" || "${SUBMIT_MODE}" == "gpu" ]]; then
  GPU_ARRAY_SPEC="${ARRAY_OVERRIDE:-1-${GPU_COUNT}}"
  GPU_JOB_ID="$(
    sbatch --parsable \
      --job-name=sf15_prof_gpu \
      --partition="${BENCHMARK_GPU_DEFAULT_PARTITION}" \
      --constraint="${BENCHMARK_GPU_CONSTRAINT}" \
      --gres="gpu:${BENCHMARK_GPU_COUNT}" \
      --nodes=1 \
      --ntasks=1 \
      --cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}" \
      --mem="${BENCHMARK_MEM}" \
      --time=04:00:00 \
      --array="${GPU_ARRAY_SPEC}" \
      --output="${RUN_ROOT}/logs/gpu_%A_%a.out" \
      --error="${RUN_ROOT}/logs/gpu_%A_%a.err" \
      --export=ALL,ECODA_PROJECT_ROOT="${PROJECT_ROOT}",SF15_RUN_ROOT="${RUN_ROOT}",SF15_MANIFEST_FILE="${MANIFEST_GPU}" \
      "${WORKER_SH}"
  )"
  echo "GPU_JOB_ID=${GPU_JOB_ID} (array=${GPU_ARRAY_SPEC})"
fi
