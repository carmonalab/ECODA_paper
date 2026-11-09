#!/bin/bash
# Explicit derivative mode for run-owned outputs and read-only input references.
# Source after slurm_config.sh and ecoda_run_common.sh; Bash 3.2 compatible.

_ecoda_derivative_cli() {
  local helper="${PROJECT_ROOT}/src/utils/py/derivative_artifacts.py"
  [[ -n "${PYTHON_BIN:-}" && -x "${PYTHON_BIN}" && -r "${helper}" ]] || {
    echo "ERROR: configured derivative-artifact interpreter/helper is unavailable" >&2
    return 1
  }
  "${PYTHON_BIN}" "${helper}" "$@" --config "${DATASETS_JSON_FILE}"
}

ecoda_derivative_bind() {
  local stage="$1" release="$2" selection="$3" inputs="$4" binding
  local -a extra=()
  [[ -n "${ECODA_RUN_ROOT:-}" && -n "${NAS_TARGET_DIR:-}" ]] || return 1
  if [[ "${stage}" == stage4 ]]; then
    [[ -n "${STAGE4_INPUT_PRODUCER_RUN_ID:-}" ]] || return 1
    extra=(--producer-run "${STAGE4_INPUT_PRODUCER_RUN_ID}")
  fi
  binding="$(_ecoda_derivative_cli bind --stage "${stage}" \
    --run-root "${ECODA_RUN_ROOT}" --release-id "${release}" \
    --nas-root "${NAS_TARGET_DIR}" --selection "${selection}" \
    --input-artifacts "${inputs}" ${extra[@]+"${extra[@]}"})" || return 1
  export ECODA_ARTIFACT_BINDING="${binding}"
  ecoda_derivative_load "${stage}"
}

ecoda_stage4_owner_key() {
  if [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
    ecoda_derivative_load stage4 || return 1
    printf 'derivative/%s/%s' "${ECODA_RUN_ID}" "$1"
  else
    printf '%s' "$1"
  fi
}

ecoda_stage4_annotation_config() {
  ecoda_derivative_load stage4 || return 1
  [[ "$1" == Stephenson ]] || return 1
  TISSUE_TYPE="$(jq -er --arg ds "$1" '.[$ds].tissue' "${DATASETS_JSON_FILE}")" || return 1
  NORMAL_TISSUE="$(jq -er --arg ds "$1" '.[$ds].normal_tissue' "${DATASETS_JSON_FILE}")" || return 1
  [[ "${TISSUE_TYPE}" == Blood && "${NORMAL_TISSUE}" == true ]] || return 1
  export TISSUE_TYPE NORMAL_TISSUE
}

ecoda_stage4_path() {
  local dataset="$1" view="$2" kind="${3:-output}" name
  if [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
    ecoda_derivative_load stage4 || return 1
    case "${kind}" in
      input) ecoda_derivative_input "${dataset}" "${view}" ;;
      output) ecoda_derivative_dataset_output "${dataset}" "${view}" ;;
      nas) ecoda_derivative_dataset_output "${dataset}" "${view}" nas ;;
      *) return 1 ;;
    esac
  else
    name="$(ecoda_view_output_name "${dataset}" "${view}")" || return 1
    printf '%s/%s/output/%s' "${HPC_SCRATCH_DIR}" "${dataset}" "${name}"
  fi
}

ecoda_stage4_parent_check() {
  local path parent
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]] || return 0
  ecoda_derivative_load stage4 || return 1
  parent="$(jq -er '.input_producer_run' "${ECODA_ARTIFACT_BINDING}")" || return 1
  path="$(ecoda_derivative_input "$1" "$2" file content)" || return 1
  ecoda_require_input_ownership "${path}" "${ECODA_RUN_ID}" &&
    ecoda_validate_input_artifact "${path}" stage3 "${parent}"
}

ecoda_stage4_prepare_copy() {
  local dataset="$1" view="$2" rc receipt
  ecoda_stage4_parent_check "${dataset}" "${view}" || return 1
  if _ecoda_derivative_cli copy-stage4 --binding "${ECODA_ARTIFACT_BINDING}" \
      --dataset "${dataset}" --view "${view}" >/dev/null; then rc=0; else rc=$?; fi
  ecoda_stage4_parent_check "${dataset}" "${view}" || return 1
  [[ ${rc} -eq 0 ]] || return "${rc}"
  receipt="${ECODA_RUN_ROOT}/manifests/stage4_input_${dataset}_${view}.json"
  if [[ -e "${receipt}.md5" || -L "${receipt}.md5" ]]; then
    ecoda_validate_checksum "${receipt}" || return 1
  else
    ecoda_write_checksum "${receipt}" || return 1
  fi
  ecoda_write_artifact_record "${receipt}" stage4_input_copy "${ECODA_RUN_ID}" >/dev/null
}

ecoda_stage4_validate_copy_receipt() {
  local dataset="$1" view="$2" receipt source output parent
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]] || return 0
  ecoda_derivative_load stage4 || return 1
  receipt="${ECODA_RUN_ROOT}/manifests/stage4_input_${dataset}_${view}.json"
  ecoda_validate_checksum "${receipt}" &&
    ecoda_validate_artifact_record "${receipt}" stage4_input_copy "${ECODA_RUN_ID}" >/dev/null || return 1
  source="$(ecoda_stage4_path "${dataset}" "${view}" input)" || return 1
  output="$(ecoda_stage4_path "${dataset}" "${view}")" || return 1
  parent="$(jq -er '.input_producer_run' "${ECODA_ARTIFACT_BINDING}")" || return 1
  ecoda_validate_checksum "${source}" || return 1
  jq -e --arg source "${source}" --arg output "${output}" --arg parent "${parent}" \
    --arg binding "$(ecoda_sha256_file "${ECODA_ARTIFACT_BINDING}")" \
    --arg md5 "${ECODA_CHECKSUM_MD5}" --argjson size "${ECODA_CHECKSUM_SIZE}" \
    '.schema == "ecoda_stage4_input_copy_v1" and .source == $source and .output == $output and
     .input_producer_run == $parent and .binding_sha256 == $binding and .source_md5 == $md5 and .size == $size' \
    "${receipt}" >/dev/null
}

ecoda_stage4_validate_nas() {
  local dataset="$1" view="$2" source destination directory
  source="$(ecoda_stage4_path "${dataset}" "${view}")" || return 1
  destination="$(ecoda_stage4_path "${dataset}" "${view}" nas)" || return 1
  ecoda_validate_artifact_record "${source}" stage4_merge "${ECODA_RUN_ID}" >/dev/null &&
    ecoda_validate_checksum "${source}" && ecoda_validate_checksum "${destination}" &&
    cmp -s "${source}" "${destination}" &&
    ecoda_validate_artifact_record "${destination}" stage4_merge "${ECODA_RUN_ID}" >/dev/null || return 1
  [[ "$(_ecoda_artifact_mode_bits "${destination}")" == 666 &&
     "$(_ecoda_artifact_mode_bits "${destination}.md5")" == 666 ]] || return 1
  directory="${destination%/*}"
  while [[ "${directory}" != "${NAS_TARGET_DIR}" ]]; do
    [[ "${directory}" == "${NAS_TARGET_DIR}/"* &&
       "$(_ecoda_artifact_mode_bits "${directory}")" == 777 ]] || return 1
    directory="${directory%/*}"
  done
}

ecoda_stage4_publish() {
  local dataset="$1" view="$2" source destination directory temporary size md5
  [[ -d "${NAS_TARGET_DIR:-}" && ! -L "${NAS_TARGET_DIR}" ]] || return 1
  ecoda_stage4_parent_check "${dataset}" "${view}" || return 1
  ecoda_stage4_validate_copy_receipt "${dataset}" "${view}" || return 1
  source="$(ecoda_stage4_path "${dataset}" "${view}")" || return 1
  destination="$(ecoda_stage4_path "${dataset}" "${view}" nas)" || return 1
  ecoda_validate_artifact_record "${source}" stage4_merge "${ECODA_RUN_ID}" >/dev/null || return 1
  ecoda_validate_checksum "${source}" || return 1
  md5="${ECODA_CHECKSUM_MD5}"; size="${ECODA_CHECKSUM_SIZE}"
  _ecoda_validate_path_ancestors "${destination}" "${NAS_TARGET_DIR}" || return 1
  [[ ! -L "${destination}" && ! -L "${destination}.md5" ]] || return 1
  mkdir -p "${destination%/*}" || return 1
  if [[ -e "${destination}" ]]; then
    cmp -s "${source}" "${destination}" || return 1
  else
    rsync -rlpt --ignore-existing --chmod=ugo+rwX "${source}" "${destination}" || return 1
  fi
  if [[ ! -e "${destination}.md5" ]]; then
    temporary="$(mktemp "${destination}.md5.build.XXXXXX")" || return 1
    printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${md5}" "${size}" "${destination}" > "${temporary}" || { rm -f "${temporary}"; return 1; }
    ln "${temporary}" "${destination}.md5" || { rm -f "${temporary}"; return 1; }
    rm -f "${temporary}"
  fi
  ecoda_validate_checksum "${destination}" || return 1
  ecoda_write_artifact_record "${destination}" stage4_merge "${ECODA_RUN_ID}" >/dev/null || return 1
  chmod 666 "${destination}" "${destination}.md5" || return 1
  directory="${destination%/*}"
  while [[ "${directory}" != "${NAS_TARGET_DIR}" ]]; do
    [[ "${directory}" == "${NAS_TARGET_DIR}/"* ]] || return 1
    chmod 777 "${directory}" || return 1
    directory="${directory%/*}"
  done
  ecoda_stage4_validate_nas "${dataset}" "${view}"
}

ecoda_stage4_finalize_failure() {
  local reason="$1" owner path state rc=0
  for owner in "${ECODA_ACQUIRED_OWNERS[@]:-}"; do
    [[ -n "${owner}" ]] || continue
    path="$(ecoda_owner_field "${owner}" PATH 2>/dev/null || true)"
    state="$(ecoda_owner_state "${owner}" 2>/dev/null || true)"
    case "${path}" in "${ECODA_RUN_ROOT}/artifacts/datasets/"*)
      [[ "${state}" != OK ]] || continue ;;
    esac
    ecoda_owner_set_state "${owner}" FAIL "${reason}" || rc=1
  done
  return "${rc}"
}

ecoda_derivative_load() {
  local stage="${1:-${ECODA_GATE_STAGE:-}}" expected
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]] || return 0
  expected="${ECODA_RUN_ROOT:-}/manifests/artifact_bindings.json"
  [[ "${ECODA_RUN_ROOT:-}" = /* && "${ECODA_ARTIFACT_BINDING}" == "${expected}" ]] || {
    echo "ERROR: derivative binding is not owned by the active run" >&2
    return 1
  }
  _ecoda_derivative_cli validate --binding "${expected}" \
    --run-root "${ECODA_RUN_ROOT}" --nas-root "${NAS_TARGET_DIR}" \
    --stage "${stage}" >/dev/null
}

ecoda_derivative_root() {
  local context="$1" destination="${2:-scratch}"
  local -a extra=()
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]] || return 1
  [[ "${destination}" == scratch || "${destination}" == nas ]] || return 1
  [[ "${destination}" == nas ]] && extra=(--nas)
  _ecoda_derivative_cli root --binding "${ECODA_ARTIFACT_BINDING}" \
    --run-root "${ECODA_RUN_ROOT}" --nas-root "${NAS_TARGET_DIR}" \
    --context "${context}" ${extra[@]+"${extra[@]}"}
}

ecoda_derivative_input() {
  local dataset="$1" view="$2" kind="${3:-file}" verify="${4:-metadata}"
  local -a extra=()
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]] || return 1
  case "${kind}" in directory) extra+=(--directory);; file) ;; *) return 1;; esac
  case "${verify}" in content) extra+=(--verify-content);; metadata) ;; *) return 1;; esac
  _ecoda_derivative_cli input --binding "${ECODA_ARTIFACT_BINDING}" \
    --run-root "${ECODA_RUN_ROOT}" --nas-root "${NAS_TARGET_DIR}" \
    --dataset "${dataset}" --view "${view}" ${extra[@]+"${extra[@]}"}
}

ecoda_derivative_dataset_output() {
  local dataset="$1" view="$2" destination="${3:-scratch}" root name path
  root="$(ecoda_derivative_root datasets "${destination}")" || return 1
  name="$(ecoda_view_output_name "${dataset}" "${view}")" || return 1
  [[ "${dataset}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
     "${name}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || return 1
  path="${root}/${dataset}/output/${name}"
  _ecoda_validate_path_ancestors "${path}" "${root}" || return 1
  printf '%s' "${path}"
}

ecoda_derivative_stage5_context() {
  case "${1-${ANALYSIS_PASS:-}}:${2-${ANALYSIS_VARIANT:-}}" in
    :) printf 'benchmark' ;;
    uncorrected:final) printf 'batch_effect/uncorrected_final' ;;
    corrected:corrected_final) printf 'batch_effect/corrected_final' ;;
    *) echo 'ERROR: derivative Stage 5 requires benchmark or an explicit final batch view' >&2; return 1 ;;
  esac
}

ecoda_derivative_stage5_validate() {
  local context expected_root expected_nas
  ecoda_derivative_load stage5 || return 1
  context="$(ecoda_derivative_stage5_context "$@")" || return 1
  expected_root="$(ecoda_derivative_root "${context}")" || return 1
  expected_nas="$(ecoda_derivative_root "${context}" nas)" || return 1
  [[ "${ANALYSIS_ROOT:-}" == "${expected_root}" &&
     "${ANALYSIS_NAS_ROOT:-}" == "${expected_nas}" ]] || {
    echo 'ERROR: Stage 5 roots disagree with the immutable derivative binding' >&2
    return 1
  }
}
