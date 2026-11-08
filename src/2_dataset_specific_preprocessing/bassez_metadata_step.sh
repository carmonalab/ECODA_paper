#!/bin/bash
# Sourced by the immutable Stage 2 submitter, worker, watchdog, and auditor.
# One explicitly selected metadata repair; no canonical scientific destinations.

_stage2_bassez_sha() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

_stage2_bassez_reference() {
  local reference="$1" input md5 mode proof parent producer proof_sha extra expected_name
  [[ "${DATASETS_JSON_FILE}" == "${ECODA_SOURCE_ROOT:-${PROJECT_ROOT}}/datasets.json" ]] || return 1
  jq -e '.Bassez.use_for_benchmark == true and (.Bassez.views | all(.[]; has("columns") | not))' "${DATASETS_JSON_FILE}" >/dev/null || return 1
  [[ -f "${reference}" && ! -L "${reference}" ]] || return 1
  awk -F '\t' 'NF != 7 || NR > 1 {exit 1} END {if (NR != 1) exit 1}' "${reference}" || return 1
  IFS=$'\t' read -r input md5 mode proof parent producer proof_sha extra < "${reference}"
  [[ -n "${input}" && -n "${proof}" && -n "${parent}" && -n "${producer}" &&
     -z "${extra}" && "${md5}" =~ ^[0-9a-f]{32}$ && "${proof_sha}" =~ ^[0-9a-f]{64}$ ]] || return 1
  [[ "${input}" = /* && -f "${input}" && ! -L "${input}" &&
     "${proof}" = /* && -s "${proof}" && ! -L "${proof}" ]] || return 1
  [[ "$(ecoda_realpath_existing "${input}")" == "${input}" &&
     "$(ecoda_realpath_existing "${proof}")" == "${proof}" ]] || return 1
  case "${input}" in "${HPC_SCRATCH_DIR%/}/"*) ;; *) return 1 ;; esac
  expected_name="$(ecoda_view_output_name Bassez benchmark_analysis)" || return 1
  [[ "${input##*/}" == "${expected_name}" && "$(_stage2_bassez_sha "${proof}")" == "${proof_sha}" ]] || return 1
  case "${mode}" in
    legacy) [[ "${parent}" == - && "${producer}" == - ]] || return 1 ;;
    producer) ecoda_validate_run_id "${parent}" || return 1; [[ "${producer}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || return 1 ;;
    *) return 1 ;;
  esac
  BASSEZ_METADATA_INPUT="${input}" BASSEZ_METADATA_MD5="${md5}"
  BASSEZ_METADATA_MODE="${mode}" BASSEZ_METADATA_PROOF="${proof}"
  BASSEZ_METADATA_PARENT_RUN="${parent}" BASSEZ_METADATA_PRODUCER="${producer}"
}

stage2_bassez_metadata_bind() {
  local reference="$1" release="$2" target binding temporary output source_root nas_base
  [[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || return 1
  source_root="${ECODA_SOURCE_ROOT:?}"
  [[ "${ECODA_RUN_ROOT}" == "${HPC_SCRATCH_DIR%/}/_ecoda_runs/${ECODA_RUN_ID}" ]] || return 1
  [[ "$(ecoda_realpath_existing "${ECODA_RUN_ROOT}")" == "${ECODA_RUN_ROOT}" &&
     "$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")" == "${HPC_SCRATCH_DIR}" ]] || return 1
  [[ "${NAS_TARGET_DIR:-}" = /* && -d "${NAS_TARGET_DIR}" && ! -L "${NAS_TARGET_DIR}" ]] || return 1
  nas_base="$(ecoda_realpath_existing "${NAS_TARGET_DIR}")" || return 1
  [[ "${nas_base}" == "${NAS_TARGET_DIR}" ]] || return 1
  _stage2_bassez_reference "${reference}" || return 1
  target="${ECODA_RUN_ROOT}/manifests/bassez_metadata.input.tsv"
  binding="${ECODA_RUN_ROOT}/manifests/bassez_metadata.json"
  [[ ! -e "${target}" && ! -L "${target}" && ! -e "${binding}" && ! -L "${binding}" ]] || return 1
  ecoda_atomic_copy "${reference}" "${target}" || return 1
  ecoda_write_checksum "${target}" || return 1
  output="${ECODA_RUN_ROOT}/artifacts/datasets/Bassez/output/$(ecoda_view_output_name Bassez benchmark_analysis)" || return 1
  temporary="$(mktemp "${binding}.build.XXXXXX")" || return 1
  if ! jq -n --arg schema ecoda_bassez_metadata_run_v1 --arg run "${ECODA_RUN_ID}" \
    --arg root "${ECODA_RUN_ROOT}" --arg scratch "${HPC_SCRATCH_DIR}" --arg nas "${nas_base}" \
    --arg release "${release}" --arg output "${output}" --arg reference "${target}" \
    --arg reference_sha "$(_stage2_bassez_sha "${target}")" \
    --arg source_sha "$(_stage2_bassez_sha "${ECODA_RUN_ROOT}/manifests/source.manifest")" \
    --arg runtime_sha "$(_stage2_bassez_sha "${ECODA_RUN_ROOT}/manifests/runtime.identity")" \
    --arg config_sha "$(_stage2_bassez_sha "${source_root}/datasets.json")" \
    '{schema:$schema,run_id:$run,run_root:$root,scratch_root:$scratch,nas_base:$nas,
      release_id:$release,output:$output,reference:$reference,reference_sha256:$reference_sha,
      source_sha256:$source_sha,runtime_sha256:$runtime_sha,config_sha256:$config_sha}' \
    > "${temporary}"; then rm -f "${temporary}"; return 1; fi
  if ! ln "${temporary}" "${binding}"; then rm -f "${temporary}"; return 1; fi
  rm -f "${temporary}"
  ecoda_write_checksum "${binding}" || return 1
  stage2_bassez_metadata_input_check
}

stage2_bassez_metadata_load() {
  local binding="${ECODA_RUN_ROOT:?}/manifests/bassez_metadata.json" reference output release source_root
  source_root="${ECODA_SOURCE_ROOT:-${PROJECT_ROOT}}"
  [[ "${ECODA_RUN_ROOT}" == "${HPC_SCRATCH_DIR%/}/_ecoda_runs/${ECODA_RUN_ID}" ]] || return 1
  [[ -f "${binding}" && ! -L "${binding}" ]] || return 1
  ecoda_validate_checksum "${binding}" >/dev/null || return 1
  reference="${ECODA_RUN_ROOT}/manifests/bassez_metadata.input.tsv"
  ecoda_validate_checksum "${reference}" >/dev/null || return 1
  output="${ECODA_RUN_ROOT}/artifacts/datasets/Bassez/output/$(ecoda_view_output_name Bassez benchmark_analysis)" || return 1
  jq -e --arg run "${ECODA_RUN_ID}" --arg root "${ECODA_RUN_ROOT}" \
    --arg scratch "${HPC_SCRATCH_DIR}" --arg nas "${NAS_TARGET_DIR}" --arg output "${output}" \
    --arg reference "${reference}" --arg reference_sha "$(_stage2_bassez_sha "${reference}")" \
    --arg source_sha "$(_stage2_bassez_sha "${ECODA_RUN_ROOT}/manifests/source.manifest")" \
    --arg runtime_sha "$(_stage2_bassez_sha "${ECODA_RUN_ROOT}/manifests/runtime.identity")" \
    --arg config_sha "$(_stage2_bassez_sha "${source_root}/datasets.json")" \
    '.schema=="ecoda_bassez_metadata_run_v1" and .run_id==$run and .run_root==$root and
     .scratch_root==$scratch and .nas_base==$nas and .output==$output and .reference==$reference and
     .reference_sha256==$reference_sha and .source_sha256==$source_sha and
     .runtime_sha256==$runtime_sha and .config_sha256==$config_sha' "${binding}" >/dev/null || return 1
  release="$(jq -er '.release_id' "${binding}")" || return 1
  [[ "${release}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || return 1
  _stage2_bassez_reference "${reference}" || return 1
  BASSEZ_METADATA_OUTPUT="${output}"
  BASSEZ_METADATA_NAS="${NAS_TARGET_DIR%/}/releases/${release}/${ECODA_RUN_ID}/datasets/Bassez/output/${output##*/}"
  BASSEZ_METADATA_RELEASE="${release}"
}

stage2_bassez_metadata_output() { stage2_bassez_metadata_load || return 1; printf '%s\n' "${BASSEZ_METADATA_OUTPUT}"; }
stage2_bassez_metadata_outputs() { stage2_bassez_metadata_load || return 1; printf '%s;%s.repair.json\n' "${BASSEZ_METADATA_OUTPUT}" "${BASSEZ_METADATA_OUTPUT}"; }
stage2_bassez_metadata_script() { printf '%s/src/2_dataset_specific_preprocessing/1.submit.sh\n' "${ECODA_SOURCE_ROOT:-${PROJECT_ROOT}}"; }
stage2_bassez_metadata_nas_output() {
  stage2_bassez_metadata_load || return 1
  case "${1:-${BASSEZ_METADATA_OUTPUT}}" in
    "${BASSEZ_METADATA_OUTPUT}") printf '%s\n' "${BASSEZ_METADATA_NAS}" ;;
    "${BASSEZ_METADATA_OUTPUT}.repair.json") printf '%s.repair.json\n' "${BASSEZ_METADATA_NAS}" ;;
    *) return 1 ;;
  esac
}
stage2_bassez_metadata_nas_path() { stage2_bassez_metadata_nas_output "$@"; }

stage2_bassez_metadata_input_check() {
  stage2_bassez_metadata_load || return 1
  ecoda_require_input_ownership "${BASSEZ_METADATA_INPUT}" "${ECODA_RUN_ID}" || return 1
  if [[ "${BASSEZ_METADATA_MODE}" == producer ]]; then
    ecoda_validate_input_artifact "${BASSEZ_METADATA_INPUT}" "${BASSEZ_METADATA_PRODUCER}" "${BASSEZ_METADATA_PARENT_RUN}" || return 1
  fi
}

_stage2_bassez_metadata_args() {
  stage2_bassez_metadata_load || return 1
  BASSEZ_METADATA_ARGS=(--input-file "${BASSEZ_METADATA_INPUT}" --expected-md5 "${BASSEZ_METADATA_MD5}"
    --scratch-root "${HPC_SCRATCH_DIR}" --run-root "${ECODA_RUN_ROOT}" --output-file "${BASSEZ_METADATA_OUTPUT}")
  if [[ "${BASSEZ_METADATA_MODE}" == legacy ]]; then
    BASSEZ_METADATA_ARGS+=(--legacy-input-proof "${BASSEZ_METADATA_PROOF}")
  else
    BASSEZ_METADATA_ARGS+=(--parent-record "${BASSEZ_METADATA_PROOF}" --parent-run-id "${BASSEZ_METADATA_PARENT_RUN}" --parent-producer "${BASSEZ_METADATA_PRODUCER}")
  fi
}

stage2_bassez_metadata_run() {
  local command="$1" worker rc
  [[ "${command}" == create || "${command}" == validate ]] || return 1
  stage2_bassez_metadata_input_check || return 1
  _stage2_bassez_metadata_args || return 1
  worker="${ECODA_SOURCE_ROOT:-${PROJECT_ROOT}}/src/2_dataset_specific_preprocessing/bassez_h5ad_metadata.py"
  [[ -f "${worker}" && ! -L "${worker}" && "${PYTHON_BIN}" = /* && -x "${PYTHON_BIN}" ]] || return 1
  if "${PYTHON_BIN}" "${worker}" "${command}" "${BASSEZ_METADATA_ARGS[@]}"; then rc=0; else rc=$?; fi
  stage2_bassez_metadata_input_check || return 1
  return "${rc}"
}

stage2_bassez_metadata_validate() {
  local path="$1" source nas
  stage2_bassez_metadata_load || return 1
  case "${path}" in
    "${BASSEZ_METADATA_OUTPUT}"|"${BASSEZ_METADATA_OUTPUT}.repair.json") stage2_bassez_metadata_run validate ;;
    "${BASSEZ_METADATA_NAS}"|"${BASSEZ_METADATA_NAS}.repair.json")
      stage2_bassez_metadata_run validate || return 1
      source="${BASSEZ_METADATA_OUTPUT}"; nas="${BASSEZ_METADATA_NAS}"
      [[ "${path}" != *.repair.json ]] || { source="${source}.repair.json"; nas="${nas}.repair.json"; }
      ecoda_validate_checksum "${source}" || return 1
      ecoda_validate_checksum "${nas}" || return 1
      [[ "$(ecoda_md5_file "${source}")" == "${ECODA_CHECKSUM_MD5}" ]]
      ;;
    *) return 1 ;;
  esac
}

stage2_bassez_metadata_sync() {
  local source destination md5 size directory mode sidecar temporary
  local -a directories
  [[ "${NAS_TARGET_DIR:-}" = /* && -d "${NAS_TARGET_DIR}" && ! -L "${NAS_TARGET_DIR}" &&
     "$(ecoda_realpath_existing "${NAS_TARGET_DIR}")" == "${NAS_TARGET_DIR}" ]] || return 1
  stage2_bassez_metadata_run validate || return 1
  directory="${BASSEZ_METADATA_NAS%/*}"
  ecoda_validate_run_owned_path "${BASSEZ_METADATA_OUTPUT}" "${ECODA_RUN_ROOT}" || return 1
  _ecoda_validate_path_ancestors "${directory}" "${NAS_TARGET_DIR}" || return 1
  directories=("${NAS_TARGET_DIR}/releases" "${NAS_TARGET_DIR}/releases/${BASSEZ_METADATA_RELEASE}"
    "${NAS_TARGET_DIR}/releases/${BASSEZ_METADATA_RELEASE}/${ECODA_RUN_ID}"
    "${NAS_TARGET_DIR}/releases/${BASSEZ_METADATA_RELEASE}/${ECODA_RUN_ID}/datasets"
    "${directory%/*}" "${directory}")
  mkdir -p "${directory}" || return 1
  for source in "${BASSEZ_METADATA_OUTPUT}" "${BASSEZ_METADATA_OUTPUT}.repair.json"; do
    destination="$(stage2_bassez_metadata_nas_output "${source}")" || return 1
    ecoda_validate_checksum "${source}" || return 1
    md5="${ECODA_CHECKSUM_MD5}"; size="${ECODA_CHECKSUM_SIZE}"
    [[ ! -L "${destination}" && ! -L "${destination}.md5" ]] || return 1
    if [[ -e "${destination}" ]]; then
      [[ "$(ecoda_md5_file "${destination}")" == "${md5}" ]] || return 1
    else
      rsync -rlpt --ignore-existing --chmod=ugo+rwX "${source}" "${destination}" || return 1
      [[ "$(ecoda_md5_file "${destination}")" == "${md5}" ]] || return 1
    fi
    sidecar="${destination}.md5"
    if [[ ! -e "${sidecar}" ]]; then
      temporary="$(mktemp "${sidecar}.build.XXXXXX")" || return 1
      printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${md5}" "${size}" "${destination}" > "${temporary}" || { rm -f "${temporary}"; return 1; }
      ln "${temporary}" "${sidecar}" || { rm -f "${temporary}"; return 1; }
      rm -f "${temporary}"
    fi
    chmod 666 "${destination}" "${sidecar}" || return 1
    ecoda_validate_checksum "${destination}" || return 1
    [[ "$(_ecoda_artifact_mode_bits "${destination}")" == 666 && "$(_ecoda_artifact_mode_bits "${sidecar}")" == 666 ]] || return 1
  done
  for directory in "${directories[@]}"; do
    chmod 777 "${directory}" || return 1
    [[ "$(_ecoda_artifact_mode_bits "${directory}")" == 777 ]] || return 1
  done
}

# Login-node publication follows the compute watchdog's scratch-only proof.
stage2_bassez_metadata_publish() {
  local path owner step_owner
  stage2_bassez_metadata_load || return 1
  _ecoda_expand_output_selection stage2 "${ECODA_RUN_ROOT}/manifests/steps.tsv" || return 1
  [[ ${#ECODA_OUTPUT_PATHS[@]} -eq 4 ]] || return 1
  for path in "${ECODA_OUTPUT_PATHS[@]}"; do
    owner="$(ecoda_artifact_owner_dir "${path}")" || return 1
    _ecoda_artifact_owner_validate_dir "${owner}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" &&
       "${ECODA_ARTIFACT_OWNER_STAGE}" == stage2 ]] || return 1
    case "${path}" in
      "${BASSEZ_METADATA_OUTPUT}"|"${BASSEZ_METADATA_OUTPUT}.repair.json")
        [[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
        ecoda_validate_artifact_record "${path}" bassez_h5ad_metadata "${ECODA_RUN_ID}" >/dev/null || return 1 ;;
      "${BASSEZ_METADATA_NAS}"|"${BASSEZ_METADATA_NAS}.repair.json")
        [[ "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE || "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1 ;;
      *) return 1 ;;
    esac
  done
  step_owner="$(ecoda_owner_dir stage2 bassez_h5ad_metadata)" || return 1
  [[ "$(ecoda_owner_run "${step_owner}")" == "${ECODA_RUN_ID}" &&
     "$(ecoda_owner_field "${step_owner}" STAGE)" == stage2 ]] || return 1
  stage2_bassez_metadata_sync || return 1
  for path in "${BASSEZ_METADATA_NAS}" "${BASSEZ_METADATA_NAS}.repair.json"; do
    stage2_bassez_metadata_validate "${path}" || return 1
    ecoda_write_artifact_record "${path}" bassez_h5ad_metadata "${ECODA_RUN_ID}" >/dev/null || return 1
    # NAS publication permissions differ from read-only scratch artifacts.
    chmod 666 "${path}" "${path}.md5" || return 1
    [[ "$(_ecoda_artifact_mode_bits "${path}")" == 666 &&
       "$(_ecoda_artifact_mode_bits "${path}.md5")" == 666 ]] || return 1
    ecoda_validate_artifact_record "${path}" bassez_h5ad_metadata "${ECODA_RUN_ID}" >/dev/null || return 1
    ecoda_artifact_owner_set_state "${path}" OK "validated and published by Stage 2 login controller" || return 1
  done
  ecoda_owner_set_state "${step_owner}" OK "scratch watchdog and NAS publication validated"
}

stage2_bassez_metadata_fail_owners() {
  local reason="$1" owner state path rc=0
  for owner in "${ECODA_ACQUIRED_OWNERS[@]:-}"; do
    [[ -n "${owner}" ]] || continue
    path="$(ecoda_owner_field "${owner}" PATH 2>/dev/null || true)"
    state="$(ecoda_owner_state "${owner}" 2>/dev/null || true)"
    # Successful scratch evidence remains reusable after publication failure.
    case "${path}" in "${ECODA_RUN_ROOT}/artifacts/datasets/Bassez/output/"*)
      [[ "${state}" != OK ]] || continue ;;
    esac
    ecoda_owner_set_state "${owner}" FAIL "${reason}" || rc=1
  done
  return "${rc}"
}
