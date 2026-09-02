#!/bin/bash
# Validate one existing H5AD on an allocated compute node.
set -euo pipefail

SCRIPT_DIR=""
if [[ -n "${PROJECT_ROOT:-}" &&
      -f "${PROJECT_ROOT}/src/slurm_config.sh" ]]; then
  SCRIPT_DIR="${PROJECT_ROOT}/src/utils/bash"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" &&
        -f "${SLURM_SUBMIT_DIR}/src/slurm_config.sh" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}/src/utils/bash"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  if [[ -n "${SLURM_JOB_ID:-}" ]] &&
     command -v scontrol >/dev/null 2>&1; then
    submitted_command="$(scontrol show job "${SLURM_JOB_ID}" -o 2>/dev/null |
      sed -n 's/.* Command=\([^ ]*\).*/\1/p' | head -1 || true)"
    submitted_dir="$(dirname "${submitted_command}")"
    if [[ -n "${submitted_command}" &&
          -f "${submitted_dir}/../../slurm_config.sh" ]]; then
      SCRIPT_DIR="$(cd "${submitted_dir}" && pwd)"
    fi
  fi
fi
if [[ -z "${SCRIPT_DIR}" || ! -f "${SCRIPT_DIR}/../../slurm_config.sh" ]]; then
  echo "ERROR: could not recover the repository source directory." >&2
  exit 1
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
ecoda_runtime_reexec_worker "${ECODA_RUNTIME_PROFILE:-stage3}" \
  "${SCRIPT_DIR}/h5ad_preflight_worker.sh" || exit 1
cd "${PROJECT_ROOT}"
if [[ -n "${H5AD_PREFLIGHT_PYTHON_BIN:-}" ]]; then
  if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" &&
        "${H5AD_PREFLIGHT_PYTHON_BIN}" != "${PYTHON_BIN}" ]]; then
    echo "ERROR: container H5AD preflight cannot override its in-image Python." >&2
    exit 1
  fi
  [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]] ||
    PYTHON_BIN="${H5AD_PREFLIGHT_PYTHON_BIN}"
fi

manifest="${H5AD_PREFLIGHT_MANIFEST:?H5AD_PREFLIGHT_MANIFEST is required}"
status_dir="${H5AD_PREFLIGHT_STATUS_DIR:?H5AD_PREFLIGHT_STATUS_DIR is required}"
run_root="${H5AD_PREFLIGHT_RUN_ROOT:-}"
run_id="${H5AD_PREFLIGHT_RUN_ID:-${run_root##*/}}"
mode="${H5AD_PREFLIGHT_MODE:-require}"
task_id="${SLURM_ARRAY_TASK_ID:-${H5AD_PREFLIGHT_TASK_ID:-}}"

[[ "${run_id}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || {
  echo "ERROR: invalid H5AD preflight run ID" >&2
  exit 1
}
if [[ -n "${run_root}" ]]; then
  [[ "${run_root##*/}" == "${run_id}" ]] || {
    echo "ERROR: H5AD preflight run ID does not match run root" >&2
    exit 1
  }
fi
[[ "${task_id}" =~ ^[0-9]+$ && ${task_id} -gt 0 ]] || {
  echo "ERROR: invalid H5AD preflight task ID" >&2
  exit 1
}
[[ -r "${manifest}" ]] || {
  echo "ERROR: H5AD preflight manifest is unreadable: ${manifest}" >&2
  exit 1
}
if [[ -n "${run_root}" ]]; then
  ecoda_validate_run_owned_path "${manifest}" "${run_root}" || {
    echo "ERROR: H5AD preflight manifest is not run-owned: ${manifest}" >&2
    exit 1
  }
fi
case "${mode}" in
  require|classify) ;;
  *) echo "ERROR: invalid H5AD preflight mode: ${mode}" >&2; exit 1 ;;
esac

line="$(sed -n "${task_id}p" "${manifest}")"
IFS=$'\t' read -r dataset view path extra <<< "${line}"
[[ -n "${dataset}" && -n "${view}" && -n "${path}" && -z "${extra}" ]] || {
  echo "ERROR: malformed H5AD preflight row ${task_id}" >&2
  exit 1
}

safe="$(printf '%s__%s' "${dataset}" "${view}" | tr '/:,\t |' '______')"
status_file="${status_dir}/${safe}.status"
mkdir -p "${status_dir}"

state="OK"
contract_rc=0
checksum_rc=0
if [[ ! -s "${path}" ]]; then
  contract_rc=1
  checksum_rc=1
else
  set +e
  "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
    --path "${path}" --view "${view}" --method "H5AD compute preflight" >/dev/null 2>&1
  contract_rc=$?
  ecoda_validate_checksum "${path}"
  checksum_rc=$?
  set -e
fi
if [[ ${contract_rc} -ne 0 || ${checksum_rc} -ne 0 ]]; then
  if [[ "${mode}" == "classify" ]]; then
    state="REBUILD"
  else
    state="FAIL"
  fi
fi

if ! ecoda_atomic_write "${status_file}" \
  "STATE=${state}\nRUN_ID=${run_id}\nDATASET=${dataset}\nVIEW=${view}\nTASK_ID=${task_id}\n"; then
  echo "ERROR: could not persist H5AD preflight status: ${status_file}" >&2
  exit 1
fi
if [[ "${state}" == "FAIL" ]]; then
  exit 1
fi
