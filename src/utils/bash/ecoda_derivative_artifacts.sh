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
  [[ -n "${ECODA_RUN_ROOT:-}" && -n "${NAS_TARGET_DIR:-}" ]] || return 1
  binding="$(_ecoda_derivative_cli bind --stage "${stage}" \
    --run-root "${ECODA_RUN_ROOT}" --release-id "${release}" \
    --nas-root "${NAS_TARGET_DIR}" --selection "${selection}" \
    --input-artifacts "${inputs}")" || return 1
  export ECODA_ARTIFACT_BINDING="${binding}"
  ecoda_derivative_load "${stage}"
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
