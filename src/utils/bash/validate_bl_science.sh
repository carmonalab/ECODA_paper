#!/bin/bash
# Validator-only evidence: no scheduler, repair, fit, or NAS publication.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
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
   "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && -n "${ECODA_RUN_ID:-}" ]] || exit 1
ecoda_validate_run_id "${ECODA_RUN_ID}" || exit 1
ecoda_require_source_script_path "${SCRIPT_DIR}/validate_bl_science.sh" "${ECODA_SOURCE_ROOT}" >/dev/null
reader="${PROJECT_ROOT}/src/utils/r/validate_bl_science.R"
ecoda_require_source_script_path "${reader}" "${ECODA_SOURCE_ROOT}" >/dev/null
[[ ! -e "${REPORT}" && ! -L "${REPORT}" && ! -e "${REPORT}.md5" && ! -L "${REPORT}.md5" ]] || {
  echo 'ERROR: existing science evidence is immutable' >&2; exit 1;
}
scratch_root="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
audit_root="${scratch_root}/_ecoda_runs/${ECODA_RUN_ID}"
[[ "${REPORT}" == "${audit_root}/artifacts/evidence/bl_science.json" &&
   ! -e "${audit_root}" && ! -L "${audit_root}" ]] || exit 1
[[ -z "${ECODA_RUN_ROOT:-}" || "$(_ecoda_canonical_path "${ECODA_RUN_ROOT}")" == "${audit_root}" ]] || exit 1
[[ "${CONTROL}" == "${scratch_root}/_ecoda_inputs/"* && -f "${CONTROL}" && ! -L "${CONTROL}" &&
   "$(ecoda_realpath_existing "${CONTROL}")" == "${CONTROL}" &&
   "$(ecoda_md5_file "${CONTROL}")" == "${CONTROL_MD5}" ]] || exit 1
ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
export ECODA_RUN_ROOT="${audit_root}" ECODA_RUN_ID RETICULATE_PYTHON="${PYTHON_BIN}"
report_owned=0
validation_failed() {
  ecoda_set_run_state FAIL bl_science_validation_failed || true
  if [[ ${report_owned} -eq 1 ]]; then
    ecoda_artifact_owner_set_state "${REPORT}" FAIL bl_science_validation_failed >/dev/null || true
  fi
}
trap validation_failed EXIT
ecoda_atomic_copy "${CONTROL}" "${ECODA_RUN_ROOT}/manifests/control.json"
[[ "$(ecoda_md5_file "${ECODA_RUN_ROOT}/manifests/control.json")" == "${CONTROL_MD5}" &&
   "$(ecoda_md5_file "${CONTROL}")" == "${CONTROL_MD5}" ]] || exit 1
ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
export ECODA_RUNTIME_PROFILE=stage5
ecoda_runtime_validate_submission apptainer
for name in control.json source.manifest runtime.identity; do
  ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/${name}"
  ecoda_write_artifact_record "${ECODA_RUN_ROOT}/manifests/${name}" bl_science_scope "${ECODA_RUN_ID}" >/dev/null
done
ecoda_runtime_validate_bound_run
mkdir -p "${REPORT%/*}"
ecoda_artifact_owner_acquire "${REPORT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
report_owned=1
${PIXI_RSCRIPT} --vanilla "${reader}" "${ECODA_RUN_ROOT}/manifests/control.json" "${CONTROL_MD5}" "${REPORT}" "${PROJECT_ROOT}"
ecoda_runtime_validate_bound_run
for name in control.json source.manifest runtime.identity; do
  ecoda_validate_artifact_record "${ECODA_RUN_ROOT}/manifests/${name}" bl_science_scope "${ECODA_RUN_ID}" >/dev/null
done
ecoda_write_checksum "${REPORT}"
ecoda_write_artifact_record "${REPORT}" bl_science_validation "${ECODA_RUN_ID}" >/dev/null
ecoda_validate_artifact_record "${REPORT}" bl_science_validation "${ECODA_RUN_ID}" >/dev/null
ecoda_artifact_owner_set_state "${REPORT}" OK bl_science_validated_publication_pending >/dev/null
ecoda_set_run_state NOOP_VALIDATED original_science_validated_original_gate_failed_publication_pending
trap - EXIT
printf 'NOOP_VALIDATED=%s\nPUBLICATION_COMPLETED=0\n' "${REPORT}"
