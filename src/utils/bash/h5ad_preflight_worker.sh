#!/bin/bash
# Validate one existing H5AD on an allocated compute node.
set -euo pipefail

source_root="${ECODA_SOURCE_ROOT:-}"
source_manifest="${ECODA_SOURCE_MANIFEST:-}"
source_required="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
run_root="${H5AD_PREFLIGHT_RUN_ROOT:-${ECODA_RUN_ROOT:-}}"
run_id="${H5AD_PREFLIGHT_RUN_ID:-${ECODA_RUN_ID:-${run_root##*/}}}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${run_root}" ]] || {
  echo "ERROR: H5AD preflight run root binding mismatch" >&2
  exit 1
}
[[ -z "${ECODA_RUN_ID:-}" || "${ECODA_RUN_ID}" == "${run_id}" ]] || {
  echo "ERROR: H5AD preflight run ID binding mismatch" >&2
  exit 1
}
[[ -z "${H5AD_PREFLIGHT_RUN_ID:-}" || "${H5AD_PREFLIGHT_RUN_ID}" == "${run_id}" ]] || {
  echo "ERROR: H5AD preflight scheduler run ID binding mismatch" >&2
  exit 1
}
worker_script=""

[[ "${source_required}" == "1" ]] || {
  echo "ERROR: legacy_source_unpinned: H5AD preflight requires an immutable source snapshot." >&2
  exit 1
}
[[ "${source_root}" = /* && "${source_manifest}" = /* ]] || {
  echo "ERROR: H5AD preflight source root and manifest are required." >&2
  exit 1
}
[[ "${run_root}" = /* ]] || {
  echo "ERROR: H5AD preflight run root is required." >&2
  exit 1
}
[[ "${run_id}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || {
  echo "ERROR: invalid H5AD preflight run ID" >&2
  exit 1
}
[[ "${run_root##*/}" == "${run_id}" ]] || {
  echo "ERROR: H5AD preflight run ID does not match run root" >&2
  exit 1
}
[[ -f "${source_root}/src/slurm_config.sh" &&
   -r "${source_root}/src/slurm_config.sh" ]] || {
  echo "ERROR: immutable H5AD preflight source root is incomplete: ${source_root}" >&2
  exit 1
}

SCRIPT_DIR="${source_root%/}/src/utils/bash"
source "${source_root}/src/slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
[[ "${ECODA_RUNS_ROOT:-}" = /* &&
   "${run_root}" == "${ECODA_RUNS_ROOT%/}/${run_id}" ]] || {
  echo "ERROR: H5AD preflight run root is not the global run root for its ID" >&2
  exit 1
}
worker_script="${source_root%/}/src/utils/bash/h5ad_preflight_worker.sh"
worker_script="$(ecoda_require_source_script_path "${worker_script}" "${source_root}")" ||
  exit 1

export ECODA_RUN_ROOT="${run_root}" ECODA_RUN_ID="${run_id}"
export ECODA_SOURCE_ROOT="${source_root}"
export ECODA_SOURCE_MANIFEST="${source_manifest}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED="${source_required}"
export ECODA_RUNTIME_IMAGE="${ECODA_RUNTIME_IMAGE:-}"
export ECODA_RUNTIME_MANIFEST="${ECODA_RUNTIME_MANIFEST:-}"
export ECODA_HOST_ENV_PREFIX="${ECODA_HOST_ENV_PREFIX:-}"
export ECODA_SCRATCH_ROOT="${ECODA_SCRATCH_ROOT:-${HPC_SCRATCH_DIR:-}}"
export ECODA_LOGS_DIR="${ECODA_LOGS_DIR:-${LOGS_DIR:-}}"
export ECODA_AUX_ROOT="${ECODA_AUX_ROOT:-${source_root%/}/aux}"
export H5AD_PREFLIGHT_RUN_ROOT="${run_root}" H5AD_PREFLIGHT_RUN_ID="${run_id}"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" &&
      "${ECODA_RUNTIME_MODE:-host}" == "host" ]]; then
  ecoda_runtime_validate_bound_run || exit 1
fi
ecoda_runtime_reexec_worker "${ECODA_RUNTIME_PROFILE:-stage3}" \
  "${worker_script}" || exit 1
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
case "${manifest}" in
  "${run_root}"/*) ;;
  *) echo "ERROR: H5AD preflight manifest escapes its run root: ${manifest}" >&2; exit 1 ;;
esac
ecoda_validate_run_owned_path "${manifest}" "${run_root}" || {
  echo "ERROR: H5AD preflight manifest is not run-owned: ${manifest}" >&2
  exit 1
}
case "${status_dir}" in
  "${run_root}"/*) ;;
  *) echo "ERROR: H5AD preflight status directory escapes its run root: ${status_dir}" >&2; exit 1 ;;
esac
mkdir -p "${status_dir}"
ecoda_validate_run_owned_path "${status_dir}" "${run_root}" || {
  echo "ERROR: H5AD preflight status directory is not run-owned: ${status_dir}" >&2
  exit 1
}
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

if [[ "${state}" == "OK" ]]; then
  case "${ECODA_RUNTIME_PROFILE:-stage3}" in
    stage3) producer="stage3_preflight" ;;
    stage5) producer="stage5_preflight" ;;
    *) producer="h5ad_preflight" ;;
  esac
  if ! ecoda_validate_artifact_record "${path}" "${producer}" "${run_id}" >/dev/null 2>&1; then
    if ! ecoda_write_artifact_record "${path}" "${producer}" "${run_id}" >/dev/null; then
      if [[ "${mode}" == "classify" ]]; then
        state="REBUILD"
      else
        state="FAIL"
      fi
    fi
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
