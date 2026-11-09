#!/bin/bash
# Fixed reviewed science publication; no scheduler, fit, or old-run write.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
CONTROL=""; CONTROL_MD5=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --control) CONTROL="${2:-}"; shift 2 ;;
    --control-md5) CONTROL_MD5="${2:-}"; shift 2 ;;
    *) echo 'Usage: publish_bl_science.sh --control ABS --control-md5 MD5' >&2; exit 1 ;;
  esac
done
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && -n "${ECODA_RUN_ID:-}" &&
   "${CONTROL_MD5}" =~ ^[a-f0-9]{32}$ ]] || exit 1
ecoda_validate_run_id "${ECODA_RUN_ID}"
helper="${PROJECT_ROOT}/src/utils/py/bl_science_publication.py"
nas_helper="${SCRIPT_DIR}/ecoda_nas_backup.sh"
for script in "${SCRIPT_DIR}/publish_bl_science.sh" "${helper}" "${nas_helper}"; do
  ecoda_require_source_script_path "${script}" "${ECODA_SOURCE_ROOT}" >/dev/null
done
scratch="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
expected_run="${scratch}/_ecoda_runs/${ECODA_RUN_ID}"
[[ "${CONTROL}" == "${scratch}/_ecoda_inputs/"* && -f "${CONTROL}" && ! -L "${CONTROL}" &&
   "$(ecoda_realpath_existing "${CONTROL}")" == "${CONTROL}" &&
   "$(ecoda_md5_file "${CONTROL}")" == "${CONTROL_MD5}" &&
   ! -e "${expected_run}" && ! -L "${expected_run}" &&
   ( -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${expected_run}" ) &&
   -d "${NAS_TARGET_DIR}" && ! -L "${NAS_TARGET_DIR}" &&
   "$(ecoda_realpath_existing "${NAS_TARGET_DIR}")" == "${NAS_TARGET_DIR}" ]] || exit 1
ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
export ECODA_RUN_ROOT="${expected_run}" ECODA_RUN_ID
OWNED=(); SYNC_OWNER=""; SYNC_LOCK=""
failed() {
  ecoda_set_run_state FAIL bl_science_publication_failed || true
  for path in "${OWNED[@]:-}"; do
    [[ -n "${path}" ]] && ecoda_artifact_owner_set_state "${path}" FAIL bl_release_incomplete >/dev/null || true
  done
  [[ -z "${SYNC_OWNER}" ]] || ecoda_owner_set_state "${SYNC_OWNER}" FAIL bl_release_incomplete || true
  [[ -z "${SYNC_LOCK}" ]] || rmdir "${SYNC_LOCK}" 2>/dev/null || true
}
trap failed EXIT
ecoda_atomic_copy "${CONTROL}" "${ECODA_RUN_ROOT}/manifests/control.json"
[[ "$(ecoda_md5_file "${CONTROL}")" == "${CONTROL_MD5}" &&
   "$(ecoda_md5_file "${ECODA_RUN_ROOT}/manifests/control.json")" == "${CONTROL_MD5}" ]]
CONTROL="${ECODA_RUN_ROOT}/manifests/control.json"
campaign_schema="$("${PYTHON_BIN}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["schema"])' "${CONTROL}")"
publication_context=benchmark
case "${campaign_schema}" in
  ecoda_bl_publication_control_v1)
    COPY_PRODUCER=stage5_bl_science_release_copy; SCOPE_PRODUCER=bl_publication_scope
    REPORT_PRODUCER=bl_science_publication; REPORT_NAME=bl_publication.json
    VALIDATION_MARKER=BL_PUBLICATION_VALIDATED ;;
  ecoda_bassez_hr22_publication_control_v1)
    COPY_PRODUCER=stage5_bassez_hr22_release_copy; SCOPE_PRODUCER=bassez_hr22_publication_scope
    REPORT_PRODUCER=bassez_hr22_science_publication; REPORT_NAME=bassez_hr22_publication.json
    VALIDATION_MARKER=BASSEZ_HR22_PUBLICATION_VALIDATED ;;
  ecoda_main_v2_login_reboot_publication_control_v1)
    COPY_PRODUCER=stage5_main_v2_reboot_science_release_copy; SCOPE_PRODUCER=main_v2_reboot_publication_scope
    REPORT_PRODUCER=main_v2_reboot_science_publication; REPORT_NAME=main_v2_reboot_science_publication.json
    VALIDATION_MARKER=MAIN_V2_REBOOT_SCIENCE_PUBLICATION_VALIDATED
    publication_context="$("${PYTHON_BIN}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["context"])' "${CONTROL}")"
    case "${publication_context}" in benchmark|batch_effect/uncorrected_final|batch_effect/corrected_final) ;; *) exit 1 ;; esac ;;
  ecoda_main_v2_science_publication_control_v1)
    COPY_PRODUCER=stage5_main_v2_science_release_copy; SCOPE_PRODUCER=main_v2_publication_scope
    REPORT_PRODUCER=main_v2_science_publication; REPORT_NAME=main_v2_science_publication.json
    VALIDATION_MARKER=MAIN_V2_SCIENCE_PUBLICATION_VALIDATED
    publication_context="$("${PYTHON_BIN}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["context"])' "${CONTROL}")"
    case "${publication_context}" in benchmark|batch_effect/uncorrected_final|batch_effect/corrected_final) ;; *) exit 1 ;; esac ;;
  *) echo 'ERROR: unsupported fixed publication campaign' >&2; exit 1 ;;
esac
ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
export ECODA_RUNTIME_PROFILE=stage5
ecoda_runtime_validate_submission apptainer
ecoda_runtime_validate_bound_run
INVENTORY="${ECODA_RUN_ROOT}/manifests/publication_inventory.tsv"
ANALYSIS_ROOT="${ECODA_RUN_ROOT}/artifacts/${publication_context}"
DESTINATION="$("${PYTHON_BIN}" "${helper}" plan --control "${CONTROL}" --control-md5 "${CONTROL_MD5}" \
  --run-root "${ECODA_RUN_ROOT}" --nas-root "${NAS_TARGET_DIR}" --inventory "${INVENTORY}")"
bash "${nas_helper}" preflight --nas-root "${NAS_TARGET_DIR}" --report "${ECODA_RUN_ROOT}/manifests/nas_preflight.json"
for name in control.json source.manifest runtime.identity publication_inventory.tsv nas_preflight.json; do
  ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/${name}"
  ecoda_write_artifact_record "${ECODA_RUN_ROOT}/manifests/${name}" "${SCOPE_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
done
SYNC_OWNER="$(ecoda_owner_acquire stage5 "sync/${ANALYSIS_ROOT}" "${ECODA_RUN_ID}" 0 1)"
mkdir "${SYNC_OWNER}.publication.lock"
SYNC_LOCK="${SYNC_OWNER}.publication.lock"
_ecoda_validate_path_ancestors "${DESTINATION}" "${NAS_TARGET_DIR}"
mkdir -p "${DESTINATION%/*}" "${ANALYSIS_ROOT}"
mkdir "${DESTINATION}"
REPORT="${ECODA_RUN_ROOT}/artifacts/evidence/${REPORT_NAME}"
mkdir -p "${REPORT%/*}"
ecoda_artifact_owner_acquire "${REPORT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
OWNED+=("${REPORT}")
CHECKSUMS="${ECODA_RUN_ROOT}/manifests/nas_checksums.md5"
: > "${CHECKSUMS}"
normalize_nas() {
  local file="$1" directory
  _ecoda_validate_path_ancestors "${file}" "${NAS_TARGET_DIR}" || return 1
  [[ -f "${file}" && ! -L "${file}" ]] || return 1
  chmod 666 "${file}" || return 1
  [[ "$(_ecoda_artifact_mode_bits "${file}")" == 666 ]] || return 1
  directory="${file%/*}"
  while [[ "${directory}" != "${NAS_TARGET_DIR}" ]]; do
    [[ "${directory}" == "${NAS_TARGET_DIR}/"* && -d "${directory}" && ! -L "${directory}" ]] || return 1
    _ecoda_validate_path_ancestors "${directory}/entry" "${NAS_TARGET_DIR}" || return 1
    chmod 777 "${directory}" || return 1
    [[ "$(_ecoda_artifact_mode_bits "${directory}")" == 777 ]] || return 1
    directory="${directory%/*}"
  done
}
while IFS=$'\t' read -r dataset kind source relative md5 size producer extra; do
  [[ -z "${extra}" ]]
  copied="${ANALYSIS_ROOT}/${relative}"; published="${DESTINATION}/${relative}"
  for target in "${copied}" "${published}"; do
    boundary="${ECODA_RUN_ROOT}"
    [[ "${target}" != "${published}" ]] || boundary="${NAS_TARGET_DIR}"
    _ecoda_validate_path_ancestors "${target}" "${boundary}"
    mkdir -p "${target%/*}"
    [[ ! -e "${target}" && ! -L "${target}" && ! -e "${target}.md5" && ! -L "${target}.md5" ]]
    ecoda_artifact_owner_acquire "${target}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
    OWNED+=("${target}")
    copy_source="${source}"
    [[ "${target}" != "${published}" ]] || copy_source="${copied}"
    "${PYTHON_BIN}" "${helper}" copy --source "${copy_source}" --destination "${target}" \
      --boundary "${boundary}" --md5 "${md5}" --size "${size}"
    ecoda_write_checksum "${target}"
    ecoda_write_artifact_record "${target}" "${COPY_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
    ecoda_validate_artifact_record "${target}" "${COPY_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
    if [[ "${target}" == "${published}" ]]; then
      normalize_nas "${target}" || exit 1
      normalize_nas "${target}.md5" || exit 1
    fi
  done
  printf '%s  %s\n' "${md5}" "${relative}" >> "${CHECKSUMS}"
  printf '%s  %s.md5\n' "$(ecoda_md5_file "${published}.md5")" "${relative}" >> "${CHECKSUMS}"
done < "${INVENTORY}"
ecoda_write_checksum "${CHECKSUMS}"
ecoda_write_artifact_record "${CHECKSUMS}" "${SCOPE_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
"${PYTHON_BIN}" "${helper}" copy --source "${CHECKSUMS}" --destination "${DESTINATION}/checksums.md5" \
  --boundary "${NAS_TARGET_DIR}" --md5 "$(ecoda_md5_file "${CHECKSUMS}")" --size "$(wc -c < "${CHECKSUMS}" | tr -d '[:space:]')"
normalize_nas "${DESTINATION}/checksums.md5" || exit 1
(cd "${DESTINATION}" && md5sum -c checksums.md5)
bash "${nas_helper}" preflight --nas-root "${NAS_TARGET_DIR}" --report "${ECODA_RUN_ROOT}/manifests/nas_postflight.json"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/nas_postflight.json"
ecoda_write_artifact_record "${ECODA_RUN_ROOT}/manifests/nas_postflight.json" "${SCOPE_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
ecoda_runtime_validate_bound_run
for name in control.json source.manifest runtime.identity publication_inventory.tsv nas_preflight.json nas_postflight.json nas_checksums.md5; do
  ecoda_validate_artifact_record "${ECODA_RUN_ROOT}/manifests/${name}" "${SCOPE_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
done
"${PYTHON_BIN}" "${helper}" verify --control "${CONTROL}" --control-md5 "${CONTROL_MD5}" \
  --run-root "${ECODA_RUN_ROOT}" --nas-root "${NAS_TARGET_DIR}" --inventory "${INVENTORY}"
for path in "${OWNED[@]}"; do
  [[ "${path}" == "${REPORT}" ]] || ecoda_artifact_owner_set_state "${path}" OK reviewed_bytes_released >/dev/null
done
"${PYTHON_BIN}" "${helper}" finish --control "${CONTROL}" --control-md5 "${CONTROL_MD5}" \
  --run-root "${ECODA_RUN_ROOT}" --nas-root "${NAS_TARGET_DIR}" --inventory "${INVENTORY}" --report "${REPORT}"
ecoda_write_checksum "${REPORT}"
ecoda_write_artifact_record "${REPORT}" "${REPORT_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
ecoda_validate_artifact_record "${REPORT}" "${REPORT_PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
ecoda_artifact_owner_set_state "${REPORT}" OK reviewed_publication_complete >/dev/null
ecoda_owner_set_state "${SYNC_OWNER}" OK reviewed_publication_complete
rmdir "${SYNC_LOCK}"; SYNC_LOCK=""
if [[ "${campaign_schema}" == ecoda_main_v2_login_reboot_publication_control_v1 ]]; then
  ecoda_set_run_state NOOP_VALIDATED reviewed_science_published_without_refit_original_reboot_stop_preserved
else
  ecoda_set_run_state NOOP_VALIDATED reviewed_science_published_without_refit_original_gate_failed
fi
trap - EXIT
printf '%s=%s\nPUBLICATION_COMPLETED=1\nNEW_MODEL_FIT=0\n' "${VALIDATION_MARKER}" "${REPORT}"
