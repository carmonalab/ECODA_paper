#!/bin/bash
# Validator-only evidence: no scheduler, repair, fit, or NAS publication.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
CONTROL=""; CONTROL_MD5=""; REPORT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --control) CONTROL="${2:-}"; shift 2 ;;
    --control-md5) CONTROL_MD5="${2:-}"; shift 2 ;;
    --report) REPORT="${2:-}"; shift 2 ;;
    *) echo 'Usage: validate_bl_science.sh --control ABS --control-md5 MD5 --report ABS' >&2; exit 1 ;;
  esac
done
[[ -n "${CONTROL}" && -n "${REPORT}" && "${CONTROL_MD5}" =~ ^[a-f0-9]{32}$ &&
   "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 ]] || exit 1
ecoda_require_source_script_path "${SCRIPT_DIR}/validate_bl_science.sh" "${ECODA_SOURCE_ROOT}" >/dev/null
reader="${PROJECT_ROOT}/src/utils/r/validate_bl_science.R"
ecoda_require_source_script_path "${reader}" "${ECODA_SOURCE_ROOT}" >/dev/null
[[ ! -e "${REPORT}" && ! -L "${REPORT}" && ! -e "${REPORT}.md5" && ! -L "${REPORT}.md5" ]] || {
  echo 'ERROR: existing science evidence is immutable' >&2; exit 1;
}
export RETICULATE_PYTHON="${PYTHON_BIN}"
audit_root="${REPORT%/artifacts/evidence/bl_science.json}"
scratch_root="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
[[ "${audit_root}" != "${REPORT}" && "${audit_root%/*}" == "${scratch_root}/_ecoda_runs" ]] || exit 1
export ECODA_RUN_ROOT="${audit_root}" ECODA_RUN_ID="${audit_root##*/}"
record="$(ecoda_artifact_record_path "${REPORT}" "${ECODA_RUN_ID}")"
[[ ! -e "${record}" && ! -L "${record}" ]] || exit 1
${PIXI_RSCRIPT} --vanilla "${reader}" "${CONTROL}" "${CONTROL_MD5}" "${REPORT}" "${PROJECT_ROOT}"
ecoda_write_checksum "${REPORT}"
ecoda_write_artifact_record "${REPORT}" bl_science_validation "${ECODA_RUN_ID}" >/dev/null
ecoda_validate_artifact_record "${REPORT}" bl_science_validation "${ECODA_RUN_ID}" >/dev/null
printf 'NOOP_VALIDATED=%s\nPUBLICATION_COMPLETED=0\n' "${REPORT}"
