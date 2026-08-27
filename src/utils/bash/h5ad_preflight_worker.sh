#!/bin/bash
# Validate one existing H5AD on an allocated compute node.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
  submitted_dir="$(scontrol show job "${SLURM_JOB_ID}" 2>/dev/null |
    awk -F= '/Command=/ {print $2}' | xargs dirname 2>/dev/null || true)"
  [[ -n "${submitted_dir}" && -d "${submitted_dir}" ]] &&
    SCRIPT_DIR="$(cd "${submitted_dir}" && pwd)"
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"
if [[ -n "${H5AD_PREFLIGHT_PYTHON_BIN:-}" ]]; then
  PYTHON_BIN="${H5AD_PREFLIGHT_PYTHON_BIN}"
fi

manifest="${H5AD_PREFLIGHT_MANIFEST:?H5AD_PREFLIGHT_MANIFEST is required}"
status_dir="${H5AD_PREFLIGHT_STATUS_DIR:?H5AD_PREFLIGHT_STATUS_DIR is required}"
run_root="${H5AD_PREFLIGHT_RUN_ROOT:-}"
mode="${H5AD_PREFLIGHT_MODE:-require}"
task_id="${SLURM_ARRAY_TASK_ID:-${H5AD_PREFLIGHT_TASK_ID:-}}"

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
  "STATE=${state}\nDATASET=${dataset}\nVIEW=${view}\nTASK_ID=${task_id}\n"; then
  echo "ERROR: could not persist H5AD preflight status: ${status_file}" >&2
  exit 1
fi
if [[ "${state}" == "FAIL" ]]; then
  exit 1
fi
