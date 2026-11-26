#!/bin/bash
set -euo pipefail

# Recover PROJECT_ROOT cleanly even when Slurm runs from /var/spool/slurmd/.
if [[ -n "${ECODA_PROJECT_ROOT:-}" ]]; then
  PROJECT_ROOT="${ECODA_PROJECT_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/src/slurm_config.sh" ]]; then
  PROJECT_ROOT="${SLURM_SUBMIT_DIR}"
else
  PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
fi

source "${PROJECT_ROOT}/src/slurm_config.sh"
cd "${PROJECT_ROOT}"

MANIFEST_FILE="${SF15_MANIFEST_FILE:?SF15_MANIFEST_FILE is required}"
RUN_ROOT="${SF15_RUN_ROOT:?SF15_RUN_ROOT is required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

LINE="$(sed -n "${TASK_ID}p" "${MANIFEST_FILE}")"
[[ -n "${LINE}" ]] || {
  echo "ERROR: missing manifest line ${TASK_ID} in ${MANIFEST_FILE}" >&2
  exit 1
}

IFS=$'\t' read -r DS_NAME METHOD_NAME LANG_MODE RESOURCE_MODE <<< "${LINE}"
[[ -n "${DS_NAME}" && -n "${METHOD_NAME}" && -n "${LANG_MODE}" && -n "${RESOURCE_MODE}" ]] || {
  echo "ERROR: malformed manifest line ${TASK_ID}: ${LINE}" >&2
  exit 1
}

export ECODA_PROF_NODE="$(hostname -s)"
export ECODA_PROF_CPU_MODEL="$(lscpu 2>/dev/null | awk -F: '/Model name/ {gsub(/^[ \t]+/, "", $2); print $2; exit}' || true)"
if [[ "${RESOURCE_MODE}" == "gpu" ]]; then
  export ECODA_PROF_GPU_MODEL="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n 1 || echo "unavailable")"
else
  export ECODA_PROF_GPU_MODEL="none"
fi

SAFE_METHOD="${METHOD_NAME//[^A-Za-z0-9_.-]/_}"
OUT_CSV="${RUN_ROOT}/rows/${RESOURCE_MODE}_${TASK_ID}_${DS_NAME}__${SAFE_METHOD}.csv"
WORK_DIR="${RUN_ROOT}/outputs/${DS_NAME}_${SAFE_METHOD}"
mkdir -p "${RUN_ROOT}/rows" "${WORK_DIR}"

echo "SF15_TASK_START task=${TASK_ID} ds=${DS_NAME} method=${METHOD_NAME} lang=${LANG_MODE} mode=${RESOURCE_MODE} node=${ECODA_PROF_NODE} cpu=${ECODA_PROF_CPU_MODEL} gpu=${ECODA_PROF_GPU_MODEL}"

HARNESS_DIR="${PROJECT_ROOT}/src/5_run_benchmark_methods/sf15_classic_profile_16rows"
TIME_CMD=""
if command -v /usr/bin/time >/dev/null 2>&1; then
  TIME_CMD="/usr/bin/time -v"
fi

if [[ "${LANG_MODE}" == "R" ]]; then
  ${TIME_CMD} ${PIXI_RSCRIPT} "${HARNESS_DIR}/profile_r_classic.R" \
    --ds_name "${DS_NAME}" \
    --method "${METHOD_NAME}" \
    --out_csv "${OUT_CSV}" \
    --work_dir "${WORK_DIR}"
elif [[ "${LANG_MODE}" == "py" ]]; then
  DEVICE_ARG="cpu"
  if [[ "${RESOURCE_MODE}" == "gpu" ]]; then
    DEVICE_ARG="cuda"
  fi
  ${TIME_CMD} "${PYTHON_BIN}" "${HARNESS_DIR}/profile_py_classic.py" \
    --ds_name "${DS_NAME}" \
    --method "${METHOD_NAME}" \
    --out_csv "${OUT_CSV}" \
    --work_dir "${WORK_DIR}" \
    --device "${DEVICE_ARG}"
else
  echo "ERROR: unsupported lang=${LANG_MODE}" >&2
  exit 1
fi
