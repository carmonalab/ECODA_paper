#!/bin/bash
# Canonical manifest-driven Stage 3 preprocessing gate.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ECODA_RUNTIME_MODE_REQUESTED="${ECODA_RUNTIME_MODE:-}"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
export ECODA_GATE_STAGE=stage3
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/h5ad_preflight_submit.sh"
cd "${PROJECT_ROOT}"

DATASETS_ARG=""
DATASETS_SET=0
VIEWS_ARG=""
VIEWS_SET=0
SELECTION_FILE_ARG=""
SELECTION_FILE_SET=0
EXACT_BATCH_SELECTION=0
COMBINED_BATCH_SELECTION=0
FORCE_ARG=0
SYNC_ONLY_RUN=""
SYNC_ONLY_SET=0
MEMORY="128G"
MAX_MEMORY="500G"
PARTITION="${SLURM_PARTITION}"
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
RUNTIME_EXPORT=""
STAGE3_SELECTION_CLASSIFICATION=""
STAGE3_UNCORRECTED_BATCH_SELECTION=0
STAGE3_COMBINED_BATCH_SELECTION=0
STAGE3_CORRECTED_SELECTION=0
STAGE3_COVID_PREFLIGHT_REQUIRED=0
STAGE3_COVID_PREFLIGHT_ROOT_REQUESTED="${STAGE3_COVID_PREFLIGHT_ROOT:-${ECODA_COVID_PREFLIGHT_ROOT:-}}"
STAGE3_COVID_PREFLIGHT_ROOT=""

usage() {
  cat <<'EOF'
Usage: 1_submit_hpc_array.sh [--datasets LIST] [--views LIST]
       [--selection-file TSV] [--combined-batch-selection]
       [--exact-batch-selection] [--force]
       [--sync-only RUN_ID] [--mem VALUE] [--max-mem VALUE]
       [--partition NAME] [--throttle N]

Each manifest row is DATASET<TAB>VIEW. --ds_name and --view remain accepted
as compatibility aliases for one dataset/view selection. Exact batch mode
requires the immutable twelve-row uncorrected selection file. A corrected-only
selection must contain every current non-underscore use_for_batch_effect
dataset's corrected row in config order. Combined batch mode requires the four
target uncorrected rows followed by every configured non-underscore
use_for_batch_effect dataset in corrected mode, in config order. An
all-corrected manifest remains a separate corrected-only release path.
EOF
}


while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets) DATASETS_ARG="${2:-}"; DATASETS_SET=1; shift 2 ;;
    --datasets=*) DATASETS_ARG="${1#*=}"; DATASETS_SET=1; shift ;;
    --ds_name) DATASETS_ARG="${2:-}"; DATASETS_SET=1; shift 2 ;;
    --ds_name=*) DATASETS_ARG="${1#*=}"; DATASETS_SET=1; shift ;;
    --views) VIEWS_ARG="${2:-}"; VIEWS_SET=1; shift 2 ;;
    --views=*) VIEWS_ARG="${1#*=}"; VIEWS_SET=1; shift ;;
    --view) VIEWS_ARG="${2:-}"; VIEWS_SET=1; shift 2 ;;
    --view=*) VIEWS_ARG="${1#*=}"; VIEWS_SET=1; shift ;;
    --selection-file) SELECTION_FILE_ARG="${2:-}"; SELECTION_FILE_SET=1; shift 2 ;;
    --selection-file=*) SELECTION_FILE_ARG="${1#*=}"; SELECTION_FILE_SET=1; shift ;;
    --combined-batch-selection) COMBINED_BATCH_SELECTION=1; shift ;;
    --exact-batch-selection) EXACT_BATCH_SELECTION=1; shift ;;
    --force) FORCE_ARG=1; shift ;;
    --sync-only) SYNC_ONLY_RUN="${2:-}"; SYNC_ONLY_SET=1; shift 2 ;;
    --sync-only=*) SYNC_ONLY_RUN="${1#*=}"; SYNC_ONLY_SET=1; shift ;;
    --mem) MEMORY="${2:-}"; shift 2 ;;
    --mem=*) MEMORY="${1#*=}"; shift ;;
    --max-mem) MAX_MEMORY="${2:-}"; shift 2 ;;
    --max-mem=*) MAX_MEMORY="${1#*=}"; shift ;;
    --partition) PARTITION="${2:-}"; shift 2 ;;
    --partition=*) PARTITION="${1#*=}"; shift ;;
    --throttle) THROTTLE="${2:-}"; shift 2 ;;
    --throttle=*) THROTTLE="${1#*=}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -n "${SYNC_ONLY_RUN}" && ${FORCE_ARG} -eq 1 ]]; then
  echo "ERROR: --sync-only cannot be combined with --force." >&2
  exit 1
fi
if [[ ${DATASETS_SET} -eq 1 && -z "${DATASETS_ARG}" ]]; then
  echo "ERROR: --datasets must not be empty." >&2
  exit 1
fi
if [[ ${VIEWS_SET} -eq 1 && -z "${VIEWS_ARG}" ]]; then
  echo "ERROR: --view/--views must not be empty." >&2
  exit 1
fi
if [[ ${SELECTION_FILE_SET} -eq 1 && -z "${SELECTION_FILE_ARG}" ]]; then
  echo "ERROR: --selection-file must not be empty." >&2
  exit 1
fi
if [[ ${SYNC_ONLY_SET} -eq 1 && -z "${SYNC_ONLY_RUN}" ]]; then
  echo "ERROR: --sync-only requires a scheduler ID or run ID." >&2
  exit 1
fi
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required for Stage 3 selection." >&2
  exit 1
fi
if [[ "${ECODA_SOURCE_ROOT:-}" = */tree &&
      -f "${ECODA_SOURCE_ROOT}/datasets.json" &&
      ! -L "${ECODA_SOURCE_ROOT}/datasets.json" &&
      -r "${ECODA_SOURCE_ROOT}/datasets.json" ]]; then
  # Selection validation must use the immutable snapshot config whenever the
  # durable wrapper has already bound one; the run-bound loader repeats this
  # check after the manifest is copied into the run root.
  DATASETS_JSON_FILE="${ECODA_SOURCE_ROOT}/datasets.json"
  export DATASETS_JSON_FILE
fi
if [[ ${COMBINED_BATCH_SELECTION} -eq 1 && ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
  echo "ERROR: --combined-batch-selection cannot be combined with --exact-batch-selection." >&2
  exit 1
fi
if [[ ${COMBINED_BATCH_SELECTION} -eq 1 &&
      ${SYNC_ONLY_SET} -eq 0 &&
      ${SELECTION_FILE_SET} -eq 0 ]]; then
  echo "ERROR: --combined-batch-selection requires --selection-file." >&2
  exit 1
fi
if [[ ${COMBINED_BATCH_SELECTION} -eq 1 &&
      (${DATASETS_SET} -eq 1 || ${VIEWS_SET} -eq 1) ]]; then
  echo "ERROR: --combined-batch-selection cannot combine --datasets/--views." >&2
  exit 1
fi

stage3_manifest_value() {
  local manifest="$1" key="$2" value
  [[ -f "${manifest}" && ! -L "${manifest}" && -r "${manifest}" ]] || return 1
  value="$(awk -v wanted="${key}" '
    index($0, wanted "=") == 1 {
      count++
      result=substr($0, length(wanted) + 2)
    }
    END {
      if (count != 1 || result == "") exit 1
      print result
    }
  ' "${manifest}")" || return 1
  printf '%s\n' "${value}"
}

stage3_copy_atomic() {
  local source="$1" destination="$2" temporary
  [[ -f "${source}" && ! -L "${source}" && -r "${source}" ]] || return 1
  mkdir -p "$(dirname "${destination}")" || return 1
  temporary="${destination}.build.$$"
  rm -f "${temporary}"
  cp "${source}" "${temporary}" || {
    rm -f "${temporary}"
    return 1
  }
  chmod 600 "${temporary}" || {
    rm -f "${temporary}"
    return 1
  }
  mv -f "${temporary}" "${destination}" || {
    rm -f "${temporary}"
    return 1
  }
}

stage3_load_bound_run() {
  local source_copy="${ECODA_RUN_ROOT:-}/manifests/source.manifest"
  local runtime_identity="${ECODA_RUN_ROOT:-}/manifests/runtime.identity"
  local source_root source_manifest_original runtime_image runtime_manifest
  local snapshot_root source_format runtime_format identity_count
  [[ -n "${ECODA_RUN_ROOT:-}" && -d "${ECODA_RUN_ROOT}" ]] || return 2
  [[ -s "${source_copy}" && ! -L "${source_copy}" && -r "${source_copy}" ]] || return 2
  [[ -s "${runtime_identity}" && ! -L "${runtime_identity}" && -r "${runtime_identity}" ]] || return 2
  ecoda_validate_run_owned_path "${source_copy}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_run_owned_path "${runtime_identity}" "${ECODA_RUN_ROOT}" || return 1
  source_format="$(stage3_manifest_value "${source_copy}" FORMAT)" || return 1
  [[ "${source_format}" == "1" ]] || return 1
  source_root="$(stage3_manifest_value "${source_copy}" SOURCE_ROOT)" || return 1
  [[ "${source_root}" = /* && "${source_root}" == */tree ]] || return 1
  snapshot_root="${source_root%/tree}"
  source_manifest_original="${snapshot_root}/identity/source.manifest"
  [[ -f "${source_manifest_original}" && ! -L "${source_manifest_original}" &&
     -r "${source_manifest_original}" ]] || return 1
  cmp -s "${source_copy}" "${source_manifest_original}" || {
    echo "ERROR: run-owned Stage 3 source manifest differs from immutable snapshot." >&2
    return 1
  }
  runtime_image="$(stage3_manifest_value "${runtime_identity}" RUNTIME_IMAGE)" || return 1
  runtime_manifest="$(stage3_manifest_value "${runtime_identity}" RUNTIME_MANIFEST)" || return 1
  [[ "${runtime_image}" = /* && "${runtime_manifest}" = /* ]] || return 1
  identity_count="$(wc -l < "${runtime_identity}" | tr -d '[:space:]')" || return 1
  runtime_format="$(_ecoda_runtime_manifest_value "${runtime_manifest}" FORMAT 2>/dev/null || true)"
  case "${runtime_format}" in
    1) [[ "${identity_count}" == "6" ]] || return 1 ;;
    2) [[ "${identity_count}" == "8" ]] || return 1 ;;
    *) return 1 ;;
  esac
  SOURCE_ROOT="${source_root}"
  SOURCE_MANIFEST_ORIGINAL="${source_manifest_original}"
  SOURCE_MANIFEST_RUN="${source_copy}"
  RUNTIME_IDENTITY="${runtime_identity}"
  export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
  export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST_ORIGINAL}"
  export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
  export ECODA_AUX_ROOT="${SOURCE_ROOT%/}/aux"
  export ECODA_RUNTIME_IMAGE="${runtime_image}"
  export ECODA_RUNTIME_MANIFEST="${runtime_manifest}"
  export ECODA_RUNTIME_IDENTITY="${RUNTIME_IDENTITY}"
  export ECODA_RUN_ID="${RUN_ID}"
  export ECODA_RUNTIME_MODE="${ECODA_RUNTIME_MODE:-apptainer}"
  PROJECT_ROOT="${SOURCE_ROOT}"
  DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
  SCRIPT_DIR="${PROJECT_ROOT}/src/3_scrnaseq_preprocessing"
  export PROJECT_ROOT DATASETS_JSON_FILE
  LOGS_DIR="${ECODA_LOGS_DIR:-${LOGS_DIR:-${ECODA_RUN_ROOT}/logs}}"
  export LOGS_DIR ECODA_LOGS_DIR="${LOGS_DIR}"
  return 0
}

stage3_record_identity_metadata() {
  local source_commit source_archive source_archive_sha source_config source_datasets
  local source_toml source_lock source_aux source_branch key value tmp
  local runtime_image runtime_manifest image_sha manifest_sha image_size manifest_size
  local image_toml image_lock
  source_commit="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SOURCE_COMMIT)" || return 1
  source_archive="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SOURCE_ARCHIVE_PATH)" || return 1
  source_archive_sha="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SOURCE_ARCHIVE_SHA256)" || return 1
  source_config="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" CONFIG_HELPER_SHA256)" || return 1
  source_datasets="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" DATASETS_SHA256)" || return 1
  source_toml="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" PIXI_TOML_SHA256)" || return 1
  source_lock="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" PIXI_LOCK_SHA256)" || return 1
  source_aux="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" AUX_ROOT)" || return 1
  source_branch="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SCGATE_DB_BRANCH)" || return 1
  runtime_image="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE)" || return 1
  runtime_manifest="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST)" || return 1
  image_sha="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE_SHA256)" || return 1
  manifest_sha="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SHA256)" || return 1
  image_size="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE_SIZE)" || return 1
  manifest_size="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SIZE)" || return 1
  image_toml="$(stage3_manifest_value "${RUNTIME_IDENTITY}" IMAGE_PIXI_TOML_SHA256 2>/dev/null || true)"
  image_lock="$(stage3_manifest_value "${RUNTIME_IDENTITY}" IMAGE_PIXI_LOCK_SHA256 2>/dev/null || true)"
  tmp="${ECODA_RUN_ROOT}/metadata.build.$$"
  {
    cat "${ECODA_RUN_ROOT}/metadata"
    printf 'SOURCE_MANIFEST=%s\nSOURCE_MANIFEST_COPY=%s\nSOURCE_ROOT=%s\nSOURCE_COMMIT=%s\nSOURCE_ARCHIVE_PATH=%s\nSOURCE_ARCHIVE_SHA256=%s\nSOURCE_CONFIG_HELPER_SHA256=%s\nSOURCE_DATASETS_SHA256=%s\nSOURCE_PIXI_TOML_SHA256=%s\nSOURCE_PIXI_LOCK_SHA256=%s\nSOURCE_AUX_ROOT=%s\nSOURCE_SCGATE_DB_BRANCH=%s\n' \
      "${SOURCE_MANIFEST_ORIGINAL}" "${SOURCE_MANIFEST_RUN}" "${SOURCE_ROOT}" \
      "${source_commit}" "${source_archive}" "${source_archive_sha}" \
      "${source_config}" "${source_datasets}" "${source_toml}" "${source_lock}" \
      "${source_aux}" "${source_branch}"
    printf 'RUNTIME_IDENTITY=%s\nRUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\n' \
      "${RUNTIME_IDENTITY}" "${runtime_image}" "${runtime_manifest}" \
      "${image_sha}" "${manifest_sha}" "${image_size}" "${manifest_size}"
    [[ -n "${image_toml}" ]] && printf 'IMAGE_PIXI_TOML_SHA256=%s\n' "${image_toml}"
    [[ -n "${image_lock}" ]] && printf 'IMAGE_PIXI_LOCK_SHA256=%s\n' "${image_lock}"
  } > "${tmp}" || {
    rm -f "${tmp}"
    return 1
  }
  mv -f "${tmp}" "${ECODA_RUN_ROOT}/metadata" || {
    rm -f "${tmp}"
    return 1
  }
}

stage3_artifact_record_valid() {
  local path="$1" producer="$2" record
  record="$(ecoda_artifact_record_path "${path}" "${RUN_ID}" 2>/dev/null || true)"
  [[ -n "${record}" && -f "${record}" && ! -L "${record}" ]] || return 1
  ecoda_validate_artifact_record "${path}" "${producer}" "${RUN_ID}"
}

stage3_artifact_record_any() {
  local path="$1"
  stage3_artifact_record_valid "${path}" stage3 ||
    stage3_artifact_record_valid "${path}" stage3_preflight
}

stage3_require_new_snapshot() {
  [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" &&
     "${ECODA_SOURCE_ROOT:-}" = /* &&
     "${ECODA_SOURCE_MANIFEST:-}" = /* ]] || {
    echo "legacy_source_unpinned" >&2
    return 1
  }
  [[ -f "${ECODA_SOURCE_MANIFEST}" && ! -L "${ECODA_SOURCE_MANIFEST}" &&
     -r "${ECODA_SOURCE_MANIFEST}" ]] || return 1
}

stage3_install_source_manifest() {
  local incoming="${ECODA_SOURCE_MANIFEST:-}"
  local destination="${ECODA_RUN_ROOT}/manifests/source.manifest"
  stage3_require_new_snapshot || return 1
  stage3_copy_atomic "${incoming}" "${destination}" || return 1
  ecoda_validate_run_owned_path "${destination}" "${ECODA_RUN_ROOT}" || return 1
  [[ "$(stage3_manifest_value "${destination}" FORMAT)" == "1" ]] || return 1
  [[ "$(stage3_manifest_value "${destination}" SOURCE_ROOT)" == "${ECODA_SOURCE_ROOT}" ]] ||
    return 1
  cmp -s "${destination}" "${incoming}" || return 1
}

stage3_source_script() {
  local relative="$1"
  printf '%s/%s' "${SOURCE_ROOT:-${ECODA_SOURCE_ROOT:-${SCRIPT_DIR}}}" "${relative}"
}

stage3_require_source_script() {
  local candidate="$1"
  [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]] || return 1
  [[ -n "${SOURCE_ROOT:-${ECODA_SOURCE_ROOT:-}}" ]] || return 1
  ecoda_require_source_script_path "${candidate}" "${SOURCE_ROOT:-${ECODA_SOURCE_ROOT}}"
}


if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
  [[ -n "${SELECTION_FILE_ARG}" ]] || {
    echo "ERROR: --exact-batch-selection requires --selection-file." >&2
    exit 1
  }
  [[ -r "${SELECTION_FILE_ARG}" ]] || {
    echo "ERROR: exact batch selection file is unreadable: ${SELECTION_FILE_ARG}" >&2
    exit 1
  }
  ecoda_validate_exact_batch_selection "${SELECTION_FILE_ARG}" 2 || exit 1
fi

output_path_for() {
  local ds="$1" view="$2" output
  output="$(ecoda_view_output_name "${ds}" "${view}")"
  [[ -n "${output}" ]] || return 1
  printf '%s/%s/output/%s' "${HPC_SCRATCH_DIR}" "${ds}" "${output}"
}
validate_external_selection() {
  local ds view row seen=""
  local selection="${1:-}"
  [[ -r "${selection}" ]] || return 1
  ecoda_validate_manifest "${selection}" 2 || return 1
  while IFS=$'\t' read -r ds view; do
    ecoda_dataset_exists "${ds}" || return 1
    ecoda_view_exists "${ds}" "${view}" || return 1
    [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" &&
      -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] || return 1
    row="${ds}/${view}"
    case " ${seen} " in *" ${row} "*) return 1 ;; esac
    seen="${seen} ${row}"
  done < "${selection}"
}

stage3_corrected_dataset_list() {
  local config="${1:-${DATASETS_JSON_FILE:-}}"
  [[ -r "${config}" && ! -L "${config}" ]] || return 1
  jq -r -e '
    if type != "object" then
      error("datasets configuration must be a JSON object")
    elif any(to_entries[]; (.value | type) != "object") then
      error("dataset entries must be JSON objects")
    else
      to_entries[]
      | select((.key | startswith("_") | not) and
               (.value.use_for_batch_effect == true))
      | .key
    end
  ' "${config}"
}

stage3_validate_selection_row() {
  local ds="$1" view="$2"
  [[ -n "${ds}" && -n "${view}" &&
     "${ds}" != *$'\t'* && "${ds}" != *$'\n'* &&
     "${view}" != *$'\t'* && "${view}" != *$'\n'* ]] || return 1
  ecoda_dataset_exists "${ds}" || return 1
  ecoda_view_exists "${ds}" "${view}" || return 1
  [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" &&
     -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] || return 1
}

stage3_validate_configured_corrected_row() {
  local config="${1:-${DATASETS_JSON_FILE:-}}" ds="${2:-}" view="${3:-}"
  local sample_col label_col input_name output_name
  [[ "${view}" == "batch_effect_corrected" && "${ds}" != _* ]] || {
    echo "ERROR: corrected Stage 3 row is not a production corrected view: ${ds}/${view}" >&2
    return 1
  }
  stage3_validate_selection_row "${ds}" "${view}" || return 1
  jq -e --arg ds "${ds}" --arg view "${view}" '
    .[$ds] as $entry
    | (($entry.columns // {}) * ($entry.views[$view].columns // {})) as $columns
    | ($columns.sample | type) == "string"
    and (($columns.sample | test("[^[:space:]]")))
    and (($columns.label | type) == "string")
    and (($columns.label | test("[^[:space:]]")))
    and ($columns.sample != $columns.label)
    and (($entry.use_for_batch_effect // false) == true)
    and (($entry.views[$view].input_file_name // $entry.views[$view].input_file)
         | type) == "string"
    and (($entry.views[$view].output_file_name // $entry.views[$view].output_file)
         | type) == "string"
  ' "${config}" >/dev/null || {
    echo "ERROR: corrected Stage 3 row lacks a complete configured sample/label/input/output contract: ${ds}/${view}" >&2
    return 1
  }
  ecoda_validate_corrected_batch_columns "${config}" "${ds}" "${view}" || return 1
  sample_col="$(jq -r --arg ds "${ds}" --arg view "${view}" \
    '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).sample' \
    "${config}")" || return 1
  label_col="$(jq -r --arg ds "${ds}" --arg view "${view}" \
    '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).label' \
    "${config}")" || return 1
  input_name="$(jq -r --arg ds "${ds}" --arg view "${view}" \
    '.[$ds].views[$view].input_file_name // .[$ds].views[$view].input_file // empty' \
    "${config}")" || return 1
  output_name="$(jq -r --arg ds "${ds}" --arg view "${view}" \
    '.[$ds].views[$view].output_file_name // .[$ds].views[$view].output_file // empty' \
    "${config}")" || return 1
  [[ -n "${sample_col}" && -n "${label_col}" &&
     -n "${input_name}" && -n "${output_name}" ]]
}

stage3_selection_is_uncorrected_four_candidate() {
  local selection="${1:-}" ds view extra
  local count=0
  local expected_datasets=(
    Covid19_PBMC
    Diabetes
    Joanito
    Lung
  )
  [[ -r "${selection}" ]] || return 1
  ecoda_validate_manifest "${selection}" 2 || return 1
  while IFS=$'\t' read -r ds view extra; do
    count=$((count + 1))
    [[ ${count} -le 4 &&
       "${ds}" == "${expected_datasets[$((count - 1))]}" &&
       "${view}" == "batch_effect_uncorrected" &&
       -z "${extra}" ]] || return 1
  done < "${selection}"
  [[ ${count} -eq 4 ]]
}

stage3_validate_uncorrected_four_selection() {
  local selection="${1:-}" config="${2:-${DATASETS_JSON_FILE:-}}"
  local ds view count=0
  [[ -r "${selection}" && -r "${config}" && ! -L "${config}" ]] || return 1
  stage3_selection_is_uncorrected_four_candidate "${selection}" || return 1
  while IFS=$'\t' read -r ds view; do
    stage3_validate_selection_row "${ds}" "${view}" || return 1
    jq -e --arg ds "${ds}" \
      '.[$ds].use_for_batch_effect == true' "${config}" >/dev/null || {
      echo "ERROR: uncorrected Stage 3 target row is not batch-enabled: ${ds}" >&2
      return 1
    }
    count=$((count + 1))
  done < "${selection}"
  [[ ${count} -eq 4 ]] || return 1
  STAGE3_SELECTION_CLASSIFICATION="uncorrected_four"
  STAGE3_UNCORRECTED_BATCH_SELECTION=1
  STAGE3_COMBINED_BATCH_SELECTION=0
  STAGE3_CORRECTED_SELECTION=0
}

stage3_validate_combined_batch_selection() {
  local selection="${1:-}" config="${2:-${DATASETS_JSON_FILE:-}}"
  local corrected_rows corrected_count expected_count count=0
  local ds view expected_ds expected_view line corrected_index
  local expected_datasets=(
    Covid19_PBMC
    Diabetes
    Joanito
    Lung
  )
  local expected_views=(
    batch_effect_uncorrected
    batch_effect_uncorrected
    batch_effect_uncorrected
    batch_effect_uncorrected
  )
  [[ -r "${selection}" ]] || return 1
  ecoda_validate_manifest "${selection}" 2 || return 1
  corrected_rows="$(stage3_corrected_dataset_list "${config}")" || return 1
  corrected_count="$(printf '%s\n' "${corrected_rows}" | awk 'NF {count++} END {print count + 0}')"
  [[ "${corrected_count}" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: configured corrected Stage 3 selection is empty." >&2
    return 1
  }
  expected_count=$((4 + corrected_count))
  while IFS= read -r line || [[ -n "${line}" ]]; do
    IFS=$'\t' read -r ds view <<< "${line}"
    count=$((count + 1))
    if [[ ${count} -le 4 ]]; then
      expected_ds="${expected_datasets[$((count - 1))]}"
      expected_view="${expected_views[$((count - 1))]}"
      [[ "${ds}" == "${expected_ds}" && "${view}" == "${expected_view}" ]] || {
        echo "ERROR: combined Stage 3 row ${count} must be ${expected_ds}/${expected_view}." >&2
        return 1
      }
      stage3_validate_selection_row "${ds}" "${view}" || {
        echo "ERROR: combined Stage 3 uncorrected row is not configured: ${ds}/${view}" >&2
        return 1
      }
      jq -e --arg ds "${ds}" '.[$ds].use_for_batch_effect == true' \
        "${config}" >/dev/null || {
        echo "ERROR: combined Stage 3 uncorrected row is not batch-enabled: ${ds}" >&2
        return 1
      }
    else
      corrected_index=$((count - 4))
      expected_ds="$(printf '%s\n' "${corrected_rows}" |
        sed -n "${corrected_index}p")"
      [[ -n "${expected_ds}" &&
         "${ds}" == "${expected_ds}" &&
         "${view}" == "batch_effect_corrected" ]] || {
        echo "ERROR: combined Stage 3 corrected row ${corrected_index} is stale or out of order." >&2
        return 1
      }
      stage3_validate_configured_corrected_row \
        "${config}" "${ds}" "${view}" || return 1
    fi
  done < "${selection}"
  [[ ${count} -eq ${expected_count} ]] || {
    echo "ERROR: combined Stage 3 selection requires exactly four uncorrected rows followed by ${corrected_count} current-config corrected rows." >&2
    return 1
  }
  STAGE3_SELECTION_CLASSIFICATION="combined"
  STAGE3_UNCORRECTED_BATCH_SELECTION=0
  STAGE3_COMBINED_BATCH_SELECTION=1
  STAGE3_CORRECTED_SELECTION=1
}

stage3_selection_is_combined_candidate() {
  local selection="${1:-}" ds view extra line expected_ds expected_view
  local total=0 uncorrected=0 corrected=0 other=0 prefix_ok=1
  local expected_datasets=(
    Covid19_PBMC
    Diabetes
    Joanito
    Lung
  )
  [[ -r "${selection}" ]] || return 1
  while IFS= read -r line || [[ -n "${line}" ]]; do
    IFS=$'\t' read -r ds view extra <<< "${line}"
    [[ -n "${ds}" && -n "${view}" && -z "${extra}" ]] || return 1
    total=$((total + 1))
    if [[ ${total} -le 4 ]]; then
      expected_ds="${expected_datasets[$((total - 1))]}"
      expected_view="batch_effect_uncorrected"
      [[ "${ds}" == "${expected_ds}" &&
         "${view}" == "${expected_view}" ]] || prefix_ok=0
    fi
    case "${view}" in
      batch_effect_uncorrected) uncorrected=$((uncorrected + 1)) ;;
      batch_effect_corrected) corrected=$((corrected + 1)) ;;
      *) other=1 ;;
    esac
  done < "${selection}"
  if [[ ${total} -gt 4 && ${prefix_ok} -eq 1 && ${corrected} -gt 0 ]]; then
    return 0
  fi
  [[ ${total} -gt 0 && ${uncorrected} -gt 0 &&
     ${corrected} -gt 0 && ${other} -eq 0 ]]
}

stage3_selection_contains_corrected() {
  local selection="${1:-}" ds view extra
  [[ -r "${selection}" ]] || return 1
  while IFS=$'\t' read -r ds view extra; do
    [[ -n "${ds}" && -n "${view}" && -z "${extra}" ]] || return 1
    [[ "${view}" == "batch_effect_corrected" ]] && return 0
  done < "${selection}"
  return 1
}

stage3_validate_corrected_selection() {
  local selection="${1:-}" config="${2:-${DATASETS_JSON_FILE:-}}"
  local corrected_rows corrected_count expected_ds ds view extra count=0
  [[ -r "${selection}" && -r "${config}" && ! -L "${config}" ]] || return 1
  ecoda_validate_manifest "${selection}" 2 || return 1
  corrected_rows="$(stage3_corrected_dataset_list "${config}")" || return 1
  corrected_count="$(printf '%s\n' "${corrected_rows}" |
    awk 'NF {count++} END {print count + 0}')"
  [[ "${corrected_count}" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: configured corrected Stage 3 selection is empty." >&2
    return 1
  }
  while IFS=$'\t' read -r ds view extra; do
    count=$((count + 1))
    expected_ds="$(printf '%s\n' "${corrected_rows}" | sed -n "${count}p")"
    [[ -n "${expected_ds}" && "${ds}" == "${expected_ds}" &&
       "${view}" == "batch_effect_corrected" && -z "${extra}" ]] || {
      echo "ERROR: corrected-only Stage 3 selection must contain every current configured corrected dataset in config order." >&2
      return 1
    }
    stage3_validate_configured_corrected_row \
      "${config}" "${ds}" "${view}" || return 1
  done < "${selection}"
  [[ ${count} -eq ${corrected_count} ]] || {
    echo "ERROR: corrected-only Stage 3 selection is partial; expected exactly ${corrected_count} configured corrected rows." >&2
    return 1
  }
  STAGE3_SELECTION_CLASSIFICATION="corrected_only"
  STAGE3_UNCORRECTED_BATCH_SELECTION=0
  STAGE3_COMBINED_BATCH_SELECTION=0
  STAGE3_CORRECTED_SELECTION=1
}

stage3_selection_contains_target_uncorrected() {
  local selection="${1:-}" ds view extra
  [[ -r "${selection}" ]] || return 1
  while IFS=$'\t' read -r ds view extra; do
    [[ -n "${ds}" && -n "${view}" && -z "${extra}" ]] || return 1
    if [[ "${view}" == "batch_effect_uncorrected" &&
          ( "${ds}" == "Covid19_PBMC" || "${ds}" == "Diabetes" ||
            "${ds}" == "Joanito" || "${ds}" == "Lung" ) ]]; then
      return 0
    fi
  done < "${selection}"
  return 1
}

stage3_classify_selection() {
  local selection="${1:-}" config="${2:-${DATASETS_JSON_FILE:-}}"
  stage3_reset_selection_classification
  # --exact-batch-selection has already validated the immutable historical
  # twelve-row matrix before this function is reached.  Keep it outside the
  # corrected release modes and therefore outside all source preflights.
  if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
    STAGE3_SELECTION_CLASSIFICATION="historical_exact"
    return 0
  fi
  if [[ ${COMBINED_BATCH_SELECTION} -eq 1 ]]; then
    stage3_validate_combined_batch_selection "${selection}" "${config}" || return 1
  elif stage3_selection_is_uncorrected_four_candidate "${selection}"; then
    stage3_validate_uncorrected_four_selection "${selection}" "${config}" || return 1
  elif stage3_selection_contains_target_uncorrected "${selection}"; then
    # Target uncorrected rows are released only as the exact four-row mode or
    # the exact four-row prefix plus the complete dynamic corrected tail.
    stage3_selection_is_combined_candidate "${selection}" &&
      stage3_validate_combined_batch_selection "${selection}" "${config}" || {
        echo "ERROR: partial target uncorrected Stage 3 selection is not an approved four-row mode." >&2
        return 1
      }
  elif stage3_selection_contains_corrected "${selection}"; then
    stage3_validate_corrected_selection "${selection}" "${config}" || return 1
  fi
}

stage3_reset_selection_classification() {
  STAGE3_SELECTION_CLASSIFICATION=""
  STAGE3_UNCORRECTED_BATCH_SELECTION=0
  STAGE3_COMBINED_BATCH_SELECTION=0
  STAGE3_CORRECTED_SELECTION=0
}

stage3_validate_covid_obs_report() {
  local view="$1" report expected_rule actual_rule input_path configured_input_path actual_md5 actual_size
  local preflight_root="${STAGE3_COVID_PREFLIGHT_ROOT:-}"
  local config="${ECODA_SOURCE_ROOT:-}/datasets.json"
  [[ "${ECODA_RUN_ROOT:-}" = /* && -d "${ECODA_RUN_ROOT}" &&
     ! -L "${ECODA_RUN_ROOT}" ]] || return 1
  [[ "${preflight_root}" = "${ECODA_RUN_ROOT}/preflight" &&
     -d "${preflight_root}" && ! -L "${preflight_root}" ]] || {
    echo "ERROR: final Stage 3 Covid obs preflight root is not run-bound." >&2
    return 1
  }
  ecoda_validate_run_owned_path "${preflight_root}" "${ECODA_RUN_ROOT}" || {
    echo "ERROR: Covid obs preflight root is not run-owned." >&2
    return 1
  }
  report="${preflight_root}/Covid19_PBMC_${view}.json"
  [[ -f "${report}" && ! -L "${report}" && -r "${report}" && -s "${report}" &&
     -f "${report}.md5" && ! -L "${report}.md5" && -r "${report}.md5" &&
     -s "${report}.md5" ]] || {
    echo "ERROR: missing Covid ${view} obs-only preflight report/checksum: ${report}" >&2
    return 1
  }
  ecoda_validate_run_owned_path "${report}" "${ECODA_RUN_ROOT}" || {
    echo "ERROR: Covid ${view} obs-only report is not run-owned." >&2
    return 1
  }
  ecoda_validate_run_owned_path "${report}.md5" "${ECODA_RUN_ROOT}" || {
    echo "ERROR: Covid ${view} obs-only checksum is not run-owned." >&2
    return 1
  }
  ecoda_validate_checksum "${report}" || {
    echo "ERROR: invalid Covid ${view} obs-only preflight checksum: ${report}" >&2
    return 1
  }
  [[ -r "${config}" && ! -L "${config}" ]] || return 1
  expected_rule="$(jq -cS --arg view "${view}" \
    '.Covid19_PBMC.views[$view].subset_vars' "${config}")" || return 1
  actual_rule="$(jq -cS --arg view "${view}" \
    '.datasets[0].subset_vars // null' "${report}")" || return 1
  [[ "${expected_rule}" != "null" && "${actual_rule}" == "${expected_rule}" ]] || {
    echo "ERROR: Covid ${view} obs-only report does not record the configured subset rule." >&2
    return 1
  }
  input_path="${HPC_SCRATCH_DIR}/Covid19_PBMC/data/Covid19_Ren2021.h5ad"
  configured_input_path="${input_path}"
  [[ -f "${configured_input_path}" && ! -L "${configured_input_path}" &&
     -r "${configured_input_path}" ]] || {
    echo "ERROR: Covid direct H5AD input is missing or unsafe: ${configured_input_path}" >&2
    return 1
  }
  input_path="$(ecoda_realpath_existing "${configured_input_path}")" || {
    echo "ERROR: Covid direct H5AD input could not be canonicalized: ${configured_input_path}" >&2
    return 1
  }
  [[ -f "${input_path}" && ! -L "${input_path}" && -r "${input_path}" ]] || {
    echo "ERROR: Covid direct H5AD input is missing or unsafe: ${input_path}" >&2
    return 1
  }
  actual_md5="$(ecoda_md5_file "${input_path}")" || return 1
  actual_size="$(wc -c < "${input_path}" | tr -d '[:space:]')" || return 1
  jq -e --arg view "${view}" \
    --arg source "${ECODA_SOURCE_ROOT}" \
    --arg manifest "${ECODA_SOURCE_MANIFEST}" \
    --arg input "${input_path}" \
    --arg md5 "${actual_md5}" --argjson size "${actual_size}" '
      .obs_only == true and
      .view == $view and
      (.datasets | length) == 1 and
      .datasets[0].dataset == "Covid19_PBMC" and
      .datasets[0].view == $view and
      .datasets[0].sample_column == "sampleID" and
      .datasets[0].input_identity.path == $input and
      .datasets[0].input_identity.md5 == $md5 and
      .datasets[0].input_identity.size == $size and
      .datasets[0].split_sample_count == 0 and
      .datasets[0].configured_cardinalities.sampleID != null and
      .datasets[0].configured_cardinalities.PatientID != null and
      .datasets[0].sampling_day_audit.column == "Sampling day (Days after symptom onset)" and
      (.datasets[0].sampling_day_audit.raw_unique_values | type) == "array" and
      (.provenance.source_root == $source) and
      (.provenance.source_manifest.path == $manifest) and
      (.provenance.runtime_identity.path | type) == "string" and
      (.provenance.runtime_manifest.path | type) == "string" and
      (.provenance.runtime_image.path | type) == "string"
    ' "${report}" >/dev/null || {
      echo "ERROR: Covid ${view} obs-only preflight identity/predicate audit failed." >&2
      return 1
    }
}

stage3_validate_covid_obs_reports() {
  stage3_validate_covid_obs_report batch_effect_uncorrected || return 1
  stage3_validate_covid_obs_report batch_effect_corrected || return 1
}

stage3_install_covid_obs_evidence() {
  local source_root="${STAGE3_COVID_PREFLIGHT_ROOT:-}"
  local destination="${ECODA_RUN_ROOT}/manifests/covid_obs_preflight"
  local manifest="${ECODA_RUN_ROOT}/manifests/covid_obs_preflight.tsv"
  local tmp="${manifest}.build.$$"
  local view report target digest size
  [[ ${STAGE3_COVID_PREFLIGHT_REQUIRED} -eq 1 ]] || return 0
  [[ "${source_root}" = /* && -d "${source_root}" && ! -L "${source_root}" ]] || return 1
  mkdir -p "${destination}" || return 1
  : > "${tmp}" || return 1
  for view in batch_effect_uncorrected batch_effect_corrected; do
    report="${source_root}/Covid19_PBMC_${view}.json"
    target="${destination}/$(basename "${report}")"
    stage3_copy_atomic "${report}" "${target}" || {
      rm -f "${tmp}"
      return 1
    }
    ecoda_write_checksum "${target}" || {
      rm -f "${tmp}"
      return 1
    }
    ecoda_validate_run_owned_path "${target}" "${ECODA_RUN_ROOT}" || {
      rm -f "${tmp}"
      return 1
    }
    digest="${ECODA_CHECKSUM_MD5}"
    size="${ECODA_CHECKSUM_SIZE}"
    printf '%s\t%s\t%s\t%s\n' "${view}" "${target}" "${digest}" "${size}" >> "${tmp}"
  done
  ecoda_atomic_install_manifest "${tmp}" "${manifest}" 4 || {
    rm -f "${tmp}"
    return 1
  }
  rm -f "${tmp}"
  ecoda_write_checksum "${manifest}"
}

stage3_selection_contains_covid() {
  local selection="${1:-}" ds view extra
  [[ -r "${selection}" ]] || return 1
  while IFS=$'\t' read -r ds view extra; do
    [[ -n "${ds}" && -n "${view}" && -z "${extra}" ]] || return 1
    [[ "${ds}" == "Covid19_PBMC" ]] && return 0
  done < "${selection}"
  return 1
}

stage3_run_covid_obs_preflight() {
  local selection="${1:-}" audit_manifest="${ECODA_RUN_ROOT}/manifests/h5ad_obs_audit.tsv"
  local audit_tmp="${audit_manifest}.build.$$"
  local status_dir="${ECODA_RUN_ROOT}/status/h5ad_obs_audit"
  local input_path="${HPC_SCRATCH_DIR}/Covid19_PBMC/data/Covid19_Ren2021.h5ad"
  local expected_root="${ECODA_RUN_ROOT}/preflight"
  local view path status safe state status_run status_dataset status_view
  local status_task status_path preflight_script preflight_id preflight_rc
  local status_count=0
  case "${STAGE3_SELECTION_CLASSIFICATION}" in
    uncorrected_four|combined) ;;
    *) return 0 ;;
  esac
  stage3_selection_contains_covid "${selection}" || return 0
  STAGE3_COVID_PREFLIGHT_REQUIRED=1
  if [[ -n "${STAGE3_COVID_PREFLIGHT_ROOT_REQUESTED}" &&
        "${STAGE3_COVID_PREFLIGHT_ROOT_REQUESTED}" != "${expected_root}" ]]; then
    echo "ERROR: Stage 3 Covid preflight root override is not run-bound." >&2
    return 1
  fi
  STAGE3_COVID_PREFLIGHT_ROOT="${expected_root}"
  mkdir -p "${STAGE3_COVID_PREFLIGHT_ROOT}" "${status_dir}" || return 1
  [[ ! -L "${STAGE3_COVID_PREFLIGHT_ROOT}" && ! -L "${status_dir}" ]] || return 1
  ecoda_validate_run_owned_path "${STAGE3_COVID_PREFLIGHT_ROOT}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_run_owned_path "${status_dir}" "${ECODA_RUN_ROOT}" || return 1
  [[ -f "${input_path}" && ! -L "${input_path}" && -r "${input_path}" ]] || {
    echo "ERROR: Covid direct H5AD input is missing or unsafe: ${input_path}" >&2
    return 1
  }
  : > "${audit_tmp}" || return 1
  printf 'Covid19_PBMC\tbatch_effect_uncorrected\t%s\n' "${input_path}" >> "${audit_tmp}"
  printf 'Covid19_PBMC\tbatch_effect_corrected\t%s\n' "${input_path}" >> "${audit_tmp}"
  ecoda_atomic_install_manifest "${audit_tmp}" "${audit_manifest}" 3 || {
    rm -f "${audit_tmp}"
    return 1
  }
  rm -f "${audit_tmp}"
  ecoda_write_checksum "${audit_manifest}" || return 1
  rm -f "${status_dir}"/*.status
  preflight_script="$(stage3_require_source_script \
    "${SCRIPT_DIR}/../utils/bash/h5ad_obs_audit_worker.sh")" || return 1
  set +e
  preflight_id="$(
    ecoda_submit_h5ad_preflight "${audit_manifest}" "${status_dir}" \
      "${ECODA_RUN_ROOT}" require "${PARTITION}" "${MEMORY}" "${THROTTLE}" \
      "${LOGS_DIR}" stage3 "${preflight_script}" "${RUNTIME_EXPORT}"
  )"
  preflight_rc=$?
  set -e
  if [[ "${preflight_id}" =~ ^[0-9]+$ ]]; then
    stage3_install_scheduler_record PREFLIGHT "${preflight_id}" || return 1
  fi
  [[ "${preflight_id}" =~ ^[0-9]+$ && ${preflight_rc} -eq 0 ]] || {
    echo "ERROR: Covid obs-only preflight worker failed: job=${preflight_id:-unknown} rc=${preflight_rc}" >&2
    return 1
  }
  ecoda_wait_h5ad_preflight_status_files "${audit_manifest}" "${status_dir}" || {
    echo "ERROR: Covid obs-only preflight statuses did not settle." >&2
    return 1
  }
  while IFS=$'\t' read -r ds view path; do
    status_count=$((status_count + 1))
    safe="$(_ecoda_safe_component "${ds}__${view}")"
    status="${status_dir}/${safe}.status"
    [[ -s "${status}" && ! -L "${status}" ]] || return 1
    state="$(sed -n 's/^STATE=//p' "${status}" | head -1)"
    status_run="$(sed -n 's/^RUN_ID=//p' "${status}" | head -1)"
    status_dataset="$(sed -n 's/^DATASET=//p' "${status}" | head -1)"
    status_view="$(sed -n 's/^VIEW=//p' "${status}" | head -1)"
    status_task="$(sed -n 's/^TASK_ID=//p' "${status}" | head -1)"
    status_path="$(sed -n 's/^INPUT_FILE=//p' "${status}" | head -1)"
    [[ "${state}" == "OK" && "${status_run}" == "${ECODA_RUN_ID}" &&
       "${status_dataset}" == "${ds}" && "${status_view}" == "${view}" &&
       "${status_task}" == "${status_count}" && "${status_path}" == "${path}" ]] || {
      echo "ERROR: malformed Covid obs-only preflight status: ${status}" >&2
      return 1
    }
  done < "${audit_manifest}"
  [[ ${status_count} -eq 2 ]] || return 1
  stage3_validate_covid_obs_reports || return 1
  stage3_install_covid_obs_evidence || return 1
}

stage3_validate_corrected_source_release() {
  local selection="${1:-}" config="${2:-${DATASETS_JSON_FILE:-}}"
  local manifest="${ECODA_RUN_ROOT}/manifests/corrected_batch_preflight.tsv"
  local manifest_tmp="${manifest}.build.$$"
  local status_dir="${ECODA_RUN_ROOT}/status/corrected_batch_preflight"
  local ds view input_name source_path source_real report report_tmp raw_tmp safe
  local sample_col label_col batch_json subset_json status source_type
  local source_md5 source_sha256 source_size
  local source_commit source_archive source_archive_sha source_config source_datasets
  local source_toml source_lock source_aux source_branch source_format
  local runtime_image runtime_manifest runtime_image_sha runtime_manifest_sha
  local runtime_image_size runtime_manifest_size rscript_bin source_script
  [[ ${STAGE3_CORRECTED_SELECTION} -eq 1 ]] || return 0
  [[ -r "${selection}" && "${ECODA_RUN_ROOT:-}" = /* &&
     -d "${ECODA_RUN_ROOT}" && ! -L "${ECODA_RUN_ROOT}" ]] || return 1
  mkdir -p "${status_dir}" || return 1
  [[ ! -L "${status_dir}" ]] || return 1
  ecoda_validate_run_owned_path "${status_dir}" "${ECODA_RUN_ROOT}" || return 1
  : > "${manifest_tmp}" || return 1
  source_format="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" FORMAT)" || return 1
  source_commit="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SOURCE_COMMIT)" || return 1
  source_archive="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SOURCE_ARCHIVE_PATH)" || return 1
  source_archive_sha="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SOURCE_ARCHIVE_SHA256)" || return 1
  source_config="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" CONFIG_HELPER_SHA256)" || return 1
  source_datasets="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" DATASETS_SHA256)" || return 1
  source_toml="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" PIXI_TOML_SHA256)" || return 1
  source_lock="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" PIXI_LOCK_SHA256)" || return 1
  source_aux="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" AUX_ROOT)" || return 1
  source_branch="$(stage3_manifest_value "${SOURCE_MANIFEST_RUN}" SCGATE_DB_BRANCH)" || return 1
  runtime_image="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE)" || return 1
  runtime_manifest="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST)" || return 1
  runtime_image_sha="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE_SHA256)" || return 1
  runtime_manifest_sha="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SHA256)" || return 1
  runtime_image_size="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE_SIZE)" || return 1
  runtime_manifest_size="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SIZE)" || return 1
  while IFS=$'\t' read -r ds view; do
    [[ "${view}" == "batch_effect_corrected" ]] || continue
    stage3_validate_configured_corrected_row \
      "${config}" "${ds}" "${view}" || return 1
    sample_col="$(jq -r --arg ds "${ds}" --arg view "${view}" \
      '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).sample' \
      "${config}")" || return 1
    label_col="$(jq -r --arg ds "${ds}" --arg view "${view}" \
      '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).label' \
      "${config}")" || return 1
    batch_json="$(jq -c --arg ds "${ds}" --arg view "${view}" \
      '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).batch' \
      "${config}")" || return 1
    subset_json="$(jq -c --arg ds "${ds}" --arg view "${view}" \
      '.[$ds].views[$view].subset_vars // {}' "${config}")" || return 1
    input_name="$(jq -r --arg ds "${ds}" --arg view "${view}" \
      '.[$ds].views[$view].input_file_name // .[$ds].views[$view].input_file // empty' \
      "${config}")" || return 1
    case "${input_name}" in
      /*) source_path="${input_name}" ;;
      *"/"*|*$'\n'*|*$'\t'*)
        echo "ERROR: corrected source input name is unsafe: ${ds}/${view}" >&2
        return 1
        ;;
      *) source_path="${HPC_SCRATCH_DIR}/${ds}/data/${input_name}" ;;
    esac
    source_real="${source_path}"
    if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" != "1" ]]; then
      [[ -f "${source_path}" && ! -L "${source_path}" &&
         -r "${source_path}" ]] || {
        echo "ERROR: corrected Stage 3 source input is missing or unsafe: ${source_path}" >&2
        return 1
      }
      source_real="$(ecoda_realpath_existing "${source_path}")" || return 1
      [[ -f "${source_real}" && ! -L "${source_real}" && -r "${source_real}" ]] || return 1
      source_md5="$(ecoda_md5_file "${source_real}")" || return 1
      source_sha256="$(ecoda_sha256_file "${source_real}")" || return 1
      source_size="$(wc -c < "${source_real}" | tr -d '[:space:]')" || return 1
    fi
    report="${ECODA_RUN_ROOT}/manifests/corrected_batch_preflight/$(
      _ecoda_safe_component "${ds}__${view}"
    ).json"
    mkdir -p "$(dirname "${report}")" || return 1
    ecoda_validate_run_owned_path "$(dirname "${report}")" "${ECODA_RUN_ROOT}" || return 1
    report_tmp="${report}.build.$$"
    raw_tmp="${report}.raw.$$"
    rm -f "${report_tmp}" "${raw_tmp}"
    status=""
    case "${source_real}" in
      *.rds|*.RDS)
        source_type="rds"
        status="SOURCE_METADATA_VALIDATED_RDS"
        if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" == "1" ]]; then
          status="CONFIG_ONLY_TEST"
          jq -n --arg dataset "${ds}" --arg view "${view}" \
            --arg source "${source_real}" --arg sample "${sample_col}" \
            --arg label "${label_col}" --arg status "${status}" \
            --arg source_root "${SOURCE_ROOT}" \
            --arg source_manifest "${SOURCE_MANIFEST_ORIGINAL}" \
            --arg runtime_identity "${RUNTIME_IDENTITY}" \
            --arg runtime_manifest "${runtime_manifest}" \
            --arg runtime_image "${runtime_image}" \
            --argjson batch "${batch_json}" --argjson subset "${subset_json}" \
            '{
              schema_version: 1, dataset: $dataset, view: $view, status: $status,
              source_type: "rds", source_path: $source,
              config: {dataset: $dataset, view: $view, sample_column: $sample,
                label_column: $label, batch_keys: $batch, subset_vars: $subset},
              provenance: {source_root: $source_root,
                source_manifest: {path: $source_manifest},
                runtime_identity: {path: $runtime_identity},
                runtime_manifest: {path: $runtime_manifest},
                runtime_image: {path: $runtime_image}}
            }' > "${report_tmp}" || return 1
          mv -f "${report_tmp}" "${report}" || return 1
        else
          source_script="$(stage3_require_source_script \
            "${SCRIPT_DIR}/../utils/r/audit_corrected_source_metadata.R")" || return 1
          rscript_bin="${PIXI_RSCRIPT%% *}"
          [[ -n "${rscript_bin}" ]] || rscript_bin="${ECODA_HOST_ENV_PREFIX:-}/bin/Rscript"
          [[ -x "${rscript_bin}" ]] || {
            echo "ERROR: immutable runtime Rscript is unavailable for corrected RDS preflight." >&2
            return 1
          }
          (
            cd "${SOURCE_ROOT}" || exit 1
            "${rscript_bin}" --vanilla "${source_script}" \
              --config "${config}" --input-file "${source_real}" \
              --output "${report}" --dataset "${ds}" \
              --view "${view}" --source-root "${SOURCE_ROOT}" \
              --source-manifest "${SOURCE_MANIFEST_ORIGINAL}" \
              --runtime-identity "${RUNTIME_IDENTITY}" \
              --run-root "${ECODA_RUN_ROOT}"
          ) || {
            rm -f "${report}" "${report}.md5"
            echo "ERROR: corrected Stage 3 RDS source metadata preflight failed: ${ds}/${view}" >&2
            return 1
          }
        fi
        ;;
      *.h5ad|*.H5AD)
        source_type="h5ad"
        status="SOURCE_METADATA_VALIDATED_H5AD"
        if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" == "1" ]]; then
          status="CONFIG_ONLY_TEST"
          jq -n --arg dataset "${ds}" --arg view "${view}" \
            --arg source "${source_real}" --arg sample "${sample_col}" \
            --arg label "${label_col}" --arg status "${status}" \
            --arg source_root "${SOURCE_ROOT}" \
            --arg source_manifest "${SOURCE_MANIFEST_ORIGINAL}" \
            --arg runtime_identity "${RUNTIME_IDENTITY}" \
            --arg runtime_manifest "${runtime_manifest}" \
            --arg runtime_image "${runtime_image}" \
            --argjson batch "${batch_json}" --argjson subset "${subset_json}" \
            '{
              schema_version: 1, dataset: $dataset, view: $view, status: $status,
              source_type: "h5ad", source_path: $source,
              config: {dataset: $dataset, view: $view, sample_column: $sample,
                label_column: $label, batch_keys: $batch, subset_vars: $subset},
              provenance: {source_root: $source_root,
                source_manifest: {path: $source_manifest},
                runtime_identity: {path: $runtime_identity},
                runtime_manifest: {path: $runtime_manifest},
                runtime_image: {path: $runtime_image}}
            }' > "${report_tmp}" || return 1
          mv -f "${report_tmp}" "${report}" || return 1
        else
          source_script="$(stage3_require_source_script \
            "${SCRIPT_DIR}/../utils/py/audit_corrected_h5ad_source.py")" || return 1
          (
            cd "${SOURCE_ROOT}" || exit 1
            export PYTHONPATH="${SOURCE_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
            "${PYTHON_BIN}" "${source_script}" \
              --config "${config}" \
              --input-file "${source_real}" \
              --output "${report}" \
              --dataset "${ds}" \
              --view "${view}" \
              --source-root "${SOURCE_ROOT}" \
              --source-manifest "${SOURCE_MANIFEST_ORIGINAL}" \
              --runtime-identity "${RUNTIME_IDENTITY}" \
              --run-root "${ECODA_RUN_ROOT}" \
              --sample-column "${sample_col}" \
              --label-column "${label_col}" \
              --batch-keys-json "${batch_json}" \
              --subset-vars-json "${subset_json}"
          ) || {
            rm -f "${report}" "${report}.md5"
            echo "ERROR: corrected Stage 3 H5AD source metadata preflight failed: ${ds}/${view}" >&2
            return 1
          }
        fi
        ;;
      *)
        echo "ERROR: corrected source metadata validation requires an RDS or H5AD input: ${ds}/${view}" >&2
        return 1
        ;;
    esac
    [[ -s "${report}" && ! -L "${report}" && -r "${report}" ]] || return 1
    ecoda_validate_run_owned_path "${report}" "${ECODA_RUN_ROOT}" || return 1
    ecoda_write_checksum "${report}" || return 1
    printf '%s\t%s\t%s\t%s\t%s\n' "${ds}" "${view}" "${status}" \
      "${source_real}" "${report}" >> "${manifest_tmp}" || return 1
    safe="$(_ecoda_safe_component "${ds}__${view}")"
    ecoda_atomic_write "${status_dir}/${safe}.status" \
      "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\nDATASET=${ds}\nVIEW=${view}\nSTATUS=${status}\nSOURCE_TYPE=${source_type}\nSOURCE_PATH=${source_real}\nREPORT=${report}\n" ||
      return 1
  done < "${selection}"
  ecoda_atomic_install_manifest "${manifest_tmp}" "${manifest}" 5 || {
    rm -f "${manifest_tmp}"
    return 1
  }
  rm -f "${manifest_tmp}"
  ecoda_write_checksum "${manifest}" || return 1
  ecoda_validate_manifest "${manifest}" 5
}

stage3_validate_corrected_source_evidence() {
  local selection="${1:-}"
  local manifest="${ECODA_RUN_ROOT}/manifests/corrected_batch_preflight.tsv"
  local status_dir="${ECODA_RUN_ROOT}/status/corrected_batch_preflight"
  local ds view extra row_ds row_view row_status row_source row_report
  local safe status expected_count=0 evidence_count=0 expected_report
  local source_real source_md5 source_sha256 source_size source_type
  local sample_col label_col batch_json batch_expected_json subset_json
  local runtime_manifest runtime_image runtime_image_sha runtime_manifest_sha
  local runtime_image_size runtime_manifest_size
  [[ ${STAGE3_CORRECTED_SELECTION} -eq 1 ]] || return 0
  [[ -r "${selection}" && -f "${manifest}" && ! -L "${manifest}" &&
     -r "${manifest}" && -s "${manifest}" ]] || return 1
  ecoda_validate_run_owned_path "${manifest}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_checksum "${manifest}" || return 1
  ecoda_validate_manifest "${manifest}" 5 || return 1
  [[ -d "${status_dir}" && ! -L "${status_dir}" ]] || return 1
  ecoda_validate_run_owned_path "${status_dir}" "${ECODA_RUN_ROOT}" || return 1
  runtime_manifest="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST)" || return 1
  runtime_image="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE)" || return 1
  runtime_image_sha="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE_SHA256)" || return 1
  runtime_manifest_sha="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SHA256)" || return 1
  runtime_image_size="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_IMAGE_SIZE)" || return 1
  runtime_manifest_size="$(stage3_manifest_value "${RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SIZE)" || return 1
  while IFS=$'\t' read -r ds view extra; do
    [[ -n "${ds}" && -n "${view}" && -z "${extra}" ]] || return 1
    [[ "${view}" == "batch_effect_corrected" ]] || continue
    expected_count=$((expected_count + 1))
    row="$(sed -n "$((expected_count))p" "${manifest}")"
    IFS=$'\t' read -r row_ds row_view row_status row_source row_report extra <<< "${row}"
    [[ "${row_ds}" == "${ds}" && "${row_view}" == "${view}" &&
       -n "${row_status}" && "${row_source}" = /* &&
       -n "${row_report}" && -z "${extra}" ]] || return 1
    expected_report="${ECODA_RUN_ROOT}/manifests/corrected_batch_preflight/$(
      _ecoda_safe_component "${ds}__${view}"
    ).json"
    [[ "${row_report}" == "${expected_report}" ]] || return 1
    case "${row_status}" in
      CONFIG_ONLY_TEST) [[ "${PREPROCESS_SUBMITTER_TEST:-0}" == "1" ]] || return 1 ;;
      SOURCE_METADATA_VALIDATED_H5AD)
        source_type="h5ad"
        ;;
      SOURCE_METADATA_VALIDATED_RDS)
        source_type="rds"
        ;;
      SOURCE_METADATA_VALIDATED_H5AD_CACHE)
        source_type="h5ad"
        ;;
      *) return 1 ;;
    esac
    [[ -f "${row_report}" && ! -L "${row_report}" &&
       -r "${row_report}" && -s "${row_report}" ]] || return 1
    ecoda_validate_run_owned_path "${row_report}" "${ECODA_RUN_ROOT}" || return 1
    ecoda_validate_checksum "${row_report}" || return 1
    safe="$(_ecoda_safe_component "${ds}__${view}")"
    status="${status_dir}/${safe}.status"
    [[ -f "${status}" && ! -L "${status}" && -r "${status}" ]] || return 1
    ecoda_validate_run_owned_path "${status}" "${ECODA_RUN_ROOT}" || return 1
    [[ "$(sed -n 's/^STATE=//p' "${status}" | head -1)" == "OK" &&
       "$(sed -n 's/^RUN_ID=//p' "${status}" | head -1)" == "${ECODA_RUN_ID}" &&
       "$(sed -n 's/^DATASET=//p' "${status}" | head -1)" == "${ds}" &&
       "$(sed -n 's/^VIEW=//p' "${status}" | head -1)" == "${view}" &&
       "$(sed -n 's/^STATUS=//p' "${status}" | head -1)" == "${row_status}" &&
       "$(sed -n 's/^REPORT=//p' "${status}" | head -1)" == "${row_report}" ]] ||
      return 1
    sample_col="$(jq -r --arg ds "${ds}" --arg view "${view}" \
      '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).sample' \
      "${DATASETS_JSON_FILE}")" || return 1
    label_col="$(jq -r --arg ds "${ds}" --arg view "${view}" \
      '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).label' \
      "${DATASETS_JSON_FILE}")" || return 1
    batch_json="$(jq -c --arg ds "${ds}" --arg view "${view}" \
      '((.[$ds].columns // {}) * (.[$ds].views[$view].columns // {})).batch' \
      "${DATASETS_JSON_FILE}")" || return 1
    subset_json="$(jq -c --arg ds "${ds}" --arg view "${view}" \
      '.[$ds].views[$view].subset_vars // {}' "${DATASETS_JSON_FILE}")" || return 1
    batch_expected_json="$(printf '%s\n' "${batch_json}" | jq -c '
      if type == "array" then .
      elif type == "string" then [.]
      else error("configured batch keys must be a string or array")
      end
    ')" || return 1
    if [[ "${row_status}" == "CONFIG_ONLY_TEST" ]]; then
      jq -e --arg dataset "${ds}" --arg view "${view}" --arg status "${row_status}" \
        --arg source "${row_source}" --arg sample "${sample_col}" \
        --arg label "${label_col}" --argjson batch "${batch_json}" \
        '.schema_version == 1 and .dataset == $dataset and .view == $view and
         .status == $status and .source_path == $source and
         .config.dataset == $dataset and .config.view == $view and
         .config.sample_column == $sample and .config.label_column == $label and
         .config.batch_keys == $batch' "${row_report}" >/dev/null || return 1
    else
      source_real="$(ecoda_realpath_existing "${row_source}")" || return 1
      [[ "${source_real}" == "${row_source}" ]] || return 1
      source_md5="$(ecoda_md5_file "${source_real}")" || return 1
      source_sha256="$(ecoda_sha256_file "${source_real}")" || return 1
      source_size="$(wc -c < "${source_real}" | tr -d '[:space:]')" || return 1
      jq -e --arg dataset "${ds}" --arg view "${view}" --arg status "${row_status}" \
        --arg source "${row_source}" --arg source_type "${source_type}" \
        --arg source_md5 "${source_md5}" --arg source_sha256 "${source_sha256}" \
        --argjson source_size "${source_size}" --arg sample "${sample_col}" \
        --arg label "${label_col}" --argjson batch "${batch_json}" \
        --argjson batch_expected "${batch_expected_json}" \
        --argjson subset "${subset_json}" --arg source_root "${SOURCE_ROOT}" \
        --arg source_manifest "${SOURCE_MANIFEST_ORIGINAL}" \
        --arg runtime_identity "${RUNTIME_IDENTITY}" \
        --arg runtime_manifest "${runtime_manifest}" --arg runtime_image "${runtime_image}" \
        --arg runtime_image_sha "${runtime_image_sha}" \
        --arg runtime_manifest_sha "${runtime_manifest_sha}" \
        --arg runtime_image_size "${runtime_image_size}" \
        --arg runtime_manifest_size "${runtime_manifest_size}" '
        .schema_version == 1 and .dataset == $dataset and .view == $view and
        .status == $status and .source_type == $source_type and
        .source_path == $source and
        .source_identity.path == $source and
        .source_identity.md5 == $source_md5 and
        .source_identity.sha256 == $source_sha256 and
        .source_identity.size == $source_size and
        .config.dataset == $dataset and .config.view == $view and
        .config.sample_column == $sample and .config.label_column == $label and
        (if $source_type == "h5ad" then
           .config.batch_keys == $batch_expected
         else
           (.config.batch_keys == $batch or
            .config.batch_keys == $batch_expected)
         end) and
        (.subset_audit as $subset_audit |
          ($subset_audit | type == "object") and
          $subset_audit.sample_column == $sample and
          ($subset_audit.total_cells | type == "number") and
          ($subset_audit.retained_cells | type == "number") and
          ($subset_audit.dropped_cells | type == "number") and
          ($subset_audit.total_samples | type == "number") and
          ($subset_audit.retained_samples | type == "number") and
          ($subset_audit.dropped_samples | type == "number") and
          ($subset_audit.total_cells >= 1) and
          ($subset_audit.retained_cells >= 1) and
          ($subset_audit.dropped_cells >= 0) and
          ($subset_audit.total_samples >= 1) and
          ($subset_audit.retained_samples >= 1) and
          ($subset_audit.dropped_samples >= 0) and
          ($subset_audit.retained_cells + $subset_audit.dropped_cells ==
            $subset_audit.total_cells) and
          ($subset_audit.retained_samples + $subset_audit.dropped_samples ==
            $subset_audit.total_samples) and
          ($subset_audit.retained_sample_ids | type == "array") and
          ($subset_audit.dropped_sample_ids | type == "array") and
          ($subset_audit.split_sample_ids | type == "array") and
          ($subset_audit.retained_sample_ids | length) == $subset_audit.retained_samples and
          ($subset_audit.dropped_sample_ids | length) == $subset_audit.dropped_samples and
          ($subset_audit.split_sample_ids | length) == 0 and
          $subset_audit.split_sample_count == 0 and
          (.source_type != "h5ad" or
            ($subset_audit.configured == $subset and
             $subset_audit.validation_scope == "full_source_metadata"))
        ) and
        (.validation_summary as $validation |
          ($validation | type == "object") and
          (($validation | keys | sort) ==
            ["batch_keys", "biological_column", "composite_design_columns",
             "composite_design_rank", "design_columns", "design_rank",
             "estimable", "fingerprint", "key_level_counts",
             "key_near_unique_fraction", "method_id", "model_id", "n_cells",
             "n_samples", "near_unique_fraction", "sample_column", "valid"]) and
          $validation.valid == true and
          $validation.estimable == true and
          $validation.sample_column == $sample and
          $validation.biological_column == $label and
          ($validation.batch_keys == $batch_expected or
           ($source_type != "h5ad" and $validation.batch_keys == $batch)) and
          ($validation.n_cells | type == "number") and
          ($validation.n_samples | type == "number") and
          ($validation.n_cells == .subset_audit.retained_cells) and
          ($validation.n_samples == .subset_audit.retained_samples) and
          ($validation.key_level_counts | type == "object") and
          ($validation.key_near_unique_fraction | type == "object") and
          (($validation.key_level_counts | keys | sort) ==
            ($batch_expected | sort)) and
          (($validation.key_near_unique_fraction | keys | sort) ==
            ($batch_expected | sort)) and
          ($validation.near_unique_fraction | type == "number") and
          ($validation.near_unique_fraction > 0 and
            $validation.near_unique_fraction < 1) and
          ($validation.design_rank | type == "number") and
          ($validation.design_columns | type == "number") and
          ($validation.composite_design_rank | type == "number") and
          ($validation.composite_design_columns | type == "number") and
          ($validation.fingerprint |
            if type == "string" then test("^[0-9a-f]{64}$") else false end) and
          $validation.method_id == "preprocess" and
          $validation.model_id == "hvg_composite_v1"
        ) and
        .provenance.source_root == $source_root and
        .provenance.source_manifest.path == $source_manifest and
        .provenance.runtime_identity.path == $runtime_identity and
        .provenance.runtime_manifest.path == $runtime_manifest and
        .provenance.runtime_manifest.sha256 == $runtime_manifest_sha and
        (.provenance.runtime_manifest.size | tostring) == $runtime_manifest_size and
        .provenance.runtime_image.path == $runtime_image and
        .provenance.runtime_image.sha256 == $runtime_image_sha and
        (.provenance.runtime_image.size | tostring) == $runtime_image_size' \
        "${row_report}" >/dev/null || {
          echo "ERROR: corrected Stage 3 source preflight report identity/content failed: ${ds}/${view}" >&2
          return 1
        }
    fi
    evidence_count=$((evidence_count + 1))
  done < "${selection}"
  [[ ${evidence_count} -gt 0 && "${evidence_count}" == "$(
    wc -l < "${manifest}" | tr -d '[:space:]'
  )" ]]
}


validate_h5ad() {
  local ds="$1" view="$2" path="$3"
  [[ -s "${path}" ]] || return 1
  "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
    --path "${path}" --view "${view}" --method "Stage 3 preprocessing" >/dev/null 2>&1
}
STAGE3_PREFLIGHT_ACTIVE=0
STAGE3_PREFLIGHT_MANIFEST=""
STAGE3_PREFLIGHT_STATUS_DIR=""

stage3_compute_validate_existing() {
  local preflight_manifest="${ECODA_RUN_ROOT}/manifests/h5ad_preflight.tsv"
  local preflight_tmp="${preflight_manifest}.build.$$"
  local status_dir="${ECODA_RUN_ROOT}/status/h5ad_preflight"
  local ds view path safe status state status_run status_dataset status_view status_task
  local preflight_id preflight_rc preflight_script count=0 status_count=0
  STAGE3_PREFLIGHT_MANIFEST="${preflight_manifest}"
  STAGE3_PREFLIGHT_STATUS_DIR="${status_dir}"
  [[ "${PREPROCESS_SUBMITTER_TEST:-0}" == 1 || ${FORCE_ARG} -eq 1 ]] && return 0
  : > "${preflight_tmp}" || return 1
  while IFS=$'\t' read -r ds view; do
    path="$(output_path_for "${ds}" "${view}")" || return 1
    if [[ -s "${path}" ]]; then
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${path}" >> "${preflight_tmp}" || return 1
      count=$((count + 1))
    fi
  done < "${MANIFEST}"
  if [[ ${count} -eq 0 ]]; then
    rm -f "${preflight_tmp}"
    return 0
  fi
  ecoda_atomic_install_manifest "${preflight_tmp}" "${preflight_manifest}" 3 || {
    rm -f "${preflight_tmp}"
    return 1
  }
  rm -f "${preflight_tmp}"
  ecoda_write_checksum "${preflight_manifest}" || return 1
  mkdir -p "${status_dir}" "${LOGS_DIR}" || return 1
  rm -f "${status_dir}"/*.status
  preflight_script="$(stage3_require_source_script \
    "${SCRIPT_DIR}/../utils/bash/h5ad_preflight_worker.sh")" || return 1
  set +e
  preflight_id="$(
    ecoda_submit_h5ad_preflight "${preflight_manifest}" "${status_dir}" \
      "${ECODA_RUN_ROOT}" classify "${PARTITION}" "${MEMORY}" "${THROTTLE}" \
      "${LOGS_DIR}" stage3 "${preflight_script}" \
      "${RUNTIME_EXPORT}"
  )"
  preflight_rc=$?
  set -e
  if [[ "${preflight_id}" =~ ^[0-9]+$ ]]; then
    stage3_install_scheduler_record PREFLIGHT "${preflight_id}" || return 1
  fi
  if [[ ${preflight_rc} -ne 0 || ! "${preflight_id}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Stage 3 H5AD preflight scheduler wait failed: job=${preflight_id:-unknown} rc=${preflight_rc}" >&2
    return 1
  fi
  ecoda_wait_h5ad_preflight_status_files "${preflight_manifest}" "${status_dir}" || {
    echo "ERROR: Stage 3 H5AD preflight statuses did not settle within ${H5AD_PREFLIGHT_STATUS_GRACE_SECONDS:-60}s" >&2
    return 1
  }
  while IFS=$'\t' read -r ds view path; do
    safe="$(_ecoda_safe_component "${ds}__${view}")"
    status="${status_dir}/${safe}.status"
    [[ -s "${status}" ]] || return 1
    status_run="$(sed -n 's/^RUN_ID=//p' "${status}" | head -1)"
    status_dataset="$(sed -n 's/^DATASET=//p' "${status}" | head -1)"
    status_view="$(sed -n 's/^VIEW=//p' "${status}" | head -1)"
    status_task="$(sed -n 's/^TASK_ID=//p' "${status}" | head -1)"
    [[ "${status_run}" == "${ECODA_RUN_ID}" &&
       "${status_dataset}" == "${ds}" &&
       "${status_view}" == "${view}" &&
       "${status_task}" == "$((status_count + 1))" ]] || return 1
    state="$(sed -n 's/^STATE=//p' "${status}" | head -1)"
    case "${state}" in
      OK)
        validate_h5ad "${ds}" "${view}" "${path}" || return 1
        stage3_artifact_record_valid "${path}" stage3_preflight || return 1
        ;;
      REBUILD) ;;
      *) return 1 ;;
    esac
    status_count=$((status_count + 1))
  done < "${preflight_manifest}"
  [[ "${status_count}" == "${count}" ]] || return 1
  STAGE3_PREFLIGHT_ACTIVE=1
}

stage3_existing_output_valid() {
  local ds="$1" view="$2" path="$3" safe status state
  if [[ ${FORCE_ARG} -ne 0 || ! -s "${path}" ]]; then
    return 1
  fi
  [[ ${STAGE3_PREFLIGHT_ACTIVE} -eq 1 ]] || return 2
  safe="$(_ecoda_safe_component "${ds}__${view}")"
  status="${STAGE3_PREFLIGHT_STATUS_DIR}/${safe}.status"
  [[ -s "${status}" ]] || return 2
  state="$(sed -n 's/^STATE=//p' "${status}" | head -1)"
  case "${state}" in
    OK)
      validate_h5ad "${ds}" "${view}" "${path}" || return 1
      stage3_artifact_record_valid "${path}" stage3_preflight || return 2
      return 0
      ;;
    REBUILD) return 1 ;;
    *) return 2 ;;
  esac
}


stage3_finalize_owner_manifest() {
  local state="$1" reason="$2" owners_file="${ECODA_RUN_ROOT:-}/manifests/owners.tsv"
  local row owner rc=0
  [[ -r "${owners_file}" ]] || return 1
  [[ -s "${owners_file}" ]] || return 0
  while IFS=$'\t' read -r row owner; do
    [[ -n "${row}" && -n "${owner}" ]] || { rc=1; continue; }
    if ! ecoda_owner_set_state "${owner}" "${state}" "${reason}"; then
      rc=1
    fi
  done < "${owners_file}"
  return "${rc}"
}

stage3_abort() {
  local reason="$1"
  local rc=0
  ecoda_owner_finalize_tracked FAIL "${reason}" || rc=1
  if [[ -n "${ECODA_RUN_ROOT:-}" && -r "${ECODA_RUN_ROOT}/manifests/owners.tsv" ]]; then
    stage3_finalize_owner_manifest FAIL "${reason}" || rc=1
  fi
  if [[ -n "${ECODA_RUN_ROOT:-}" ]]; then
    ecoda_set_run_state FAIL "${reason}" || rc=1
  fi
  echo "ERROR: ${reason}" >&2
  exit 1
}
stage3_validate_scheduler_manifest() {
  local manifest="$1" require_jobs="${2:-1}"
  local kind scheduler_id seen="" count=0 array_count=0 watchdog_count=0
  [[ -r "${manifest}" ]] || return 1
  ecoda_validate_run_owned_path "${manifest}" "${ECODA_RUN_ROOT}" || return 1
  if [[ ! -s "${manifest}" ]]; then
    [[ "${require_jobs}" == "0" ]] || return 1
    return 0
  fi
  ecoda_validate_manifest "${manifest}" 2 || return 1
  while IFS=$'\t' read -r kind scheduler_id; do
    [[ -n "${kind}" && "${scheduler_id}" =~ ^[0-9]+$ ]] || return 1
    case "${kind}" in
      ARRAY) array_count=$((array_count + 1)) ;;
      WATCHDOG) watchdog_count=$((watchdog_count + 1)) ;;
      STATUS|PREFLIGHT) ;;
      *) return 1 ;;
    esac
    case " ${seen} " in
      *" ${scheduler_id} "*) return 1 ;;
    esac
    seen="${seen} ${scheduler_id}"
    count=$((count + 1))
  done < "${manifest}"
  [[ "${require_jobs}" == "0" || ${array_count} -gt 0 ]] || return 1
  [[ "${require_jobs}" == "0" || ${watchdog_count} -gt 0 ]]
}

stage3_install_scheduler_record() {
  local kind="$1" scheduler_id="$2"
  local manifest="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
  local tmp="${manifest}.record.$$" existing_kind existing_id
  [[ "${kind}" == "ARRAY" || "${kind}" == "WATCHDOG" ||
     "${kind}" == "STATUS" || "${kind}" == "PREFLIGHT" ]] ||
    return 1
  [[ "${scheduler_id}" =~ ^[0-9]+$ ]] || return 1
  if [[ -s "${manifest}" ]]; then
    while IFS=$'\t' read -r existing_kind existing_id; do
      [[ -n "${existing_kind}" && -n "${existing_id}" ]] || return 1
      if [[ "${existing_id}" == "${scheduler_id}" ]]; then
        [[ "${existing_kind}" == "${kind}" || "${kind}" == "STATUS" ]] || return 1
        return 0
      fi
    done < "${manifest}"
  fi
  if [[ -s "${manifest}" ]]; then
    cp "${manifest}" "${tmp}" || return 1
  else
    : > "${tmp}" || return 1
  fi
  printf '%s\t%s\n' "${kind}" "${scheduler_id}" >> "${tmp}" || {
    rm -f "${tmp}"
    return 1
  }
  if ! ecoda_atomic_install_manifest "${tmp}" "${manifest}" 2; then
    rm -f "${tmp}"
    return 1
  fi
  rm -f "${tmp}"
  stage3_validate_scheduler_manifest "${manifest}" 0
}

stage3_record_watchdog_status_ids() {
  local status_file="${ECODA_RUN_ROOT}/status/watchdog" status_line scheduler_id
  [[ -r "${status_file}" ]] || return 0
  while IFS= read -r status_line; do
    case "${status_line}" in
      SCHEDULER_ID=*)
        scheduler_id="${status_line#*=}"
        stage3_install_scheduler_record STATUS "${scheduler_id}" || return 1
        ;;
    esac
  done < "${status_file}"
}

stage3_watchdog_terminal_ok() {
  local scheduler_manifest="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
  local status_file="${ECODA_RUN_ROOT}/status/watchdog"
  local state status_run watchdog_id rows exitcode
  if [[ -s "${status_file}" ]]; then
    state="$(sed -n 's/^STATE=//p' "${status_file}" | head -1)"
    status_run="$(sed -n 's/^RUN_ID=//p' "${status_file}" | head -1)"
    [[ "${status_run}" == "${ECODA_RUN_ID}" ]] || return 1
    case "${state}" in
      OK) stage3_record_watchdog_status_ids || return 1; return 0 ;;
      FAIL) return 1 ;;
      ACTIVE|"") ;;
      *) return 1 ;;
    esac
  fi
  watchdog_id="$(awk -F '\t' '$1 == "WATCHDOG" {print $2; exit}' "${scheduler_manifest}")"
  [[ "${watchdog_id}" =~ ^[0-9]+$ ]] || return 1
  ecoda_wait_scalar_accounting "${watchdog_id}" \
    "${STAGE3_WATCHDOG_POLL_SECONDS:-30}" || return 1
  rows="${ECODA_ACCOUNTING_ROWS:-}"
  exitcode="$(printf '%s\n' "${rows}" | awk -F '|' 'NR == 1 {print $3}')"
  [[ "${ECODA_ACCOUNTING_STATE:-}" == "COMPLETED" &&
     "${exitcode}" == 0:0* ]]
}



sync_selected() {
  local manifest="$1" ds view path output remote_dir sync_lock lock_root
  local transfer_manifest verify_manifest expected_output local_digest local_size
  local rc=0
  [[ -d "${NAS_TARGET_DIR}" ]] || {
    echo "ERROR: NAS path is unreachable: ${NAS_TARGET_DIR}" >&2
    return 1
  }
  lock_root="${SYNC_LOCK_ROOT:-${ECODA_RUN_ROOT:-${TMPDIR:-/tmp}/ecoda-stage3-sync}}"
  mkdir -p "${lock_root}" || return 1
  sync_lock="${lock_root}/sync.lock"
  mkdir "${sync_lock}" 2>/dev/null || {
    echo "ERROR: another Stage 3 sync owns ${sync_lock}" >&2
    return 1
  }
  transfer_manifest="${sync_lock}/files.$$"
  verify_manifest="${sync_lock}/verify.$$"
  : > "${transfer_manifest}" || rc=1
  : > "${verify_manifest}" || rc=1
  while IFS=$'\t' read -r ds view; do
    [[ -n "${ds}" && -n "${view}" ]] || { rc=1; continue; }
    if ! path="$(output_path_for "${ds}" "${view}")"; then
      rc=1
      continue
    fi
    output="$(basename "${path}")"
    remote_dir="${NAS_TARGET_DIR}/${ds}/output"
    expected_output="${ds}/output/${output}"
    if ! mkdir -p "${remote_dir}" || ! validate_h5ad "${ds}" "${view}" "${path}"; then
      rc=1
      continue
    fi
    if [[ -n "${ECODA_RUN_ROOT:-}" &&
          -s "${ECODA_RUN_ROOT}/manifests/source.manifest" ]]; then
      stage3_artifact_record_any "${path}" || {
        rc=1
        continue
      }
    else
      ecoda_validate_checksum "${path}" || {
        rc=1
        continue
      }
    fi
    local_digest="${ECODA_CHECKSUM_MD5}"
    local_size="${ECODA_CHECKSUM_SIZE}"
    printf '%s\n' "${expected_output}" >> "${transfer_manifest}" || rc=1
    printf '%s\n' "${expected_output}.md5" >> "${transfer_manifest}" || rc=1
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${ds}" "${view}" "${path}" "${output}" \
      "${local_digest}" "${local_size}" >> "${verify_manifest}" || rc=1
  done < "${manifest}"
  if [[ ${rc} -eq 0 && -s "${transfer_manifest}" ]]; then
    rsync -rlptDv --files-from="${transfer_manifest}" \
      "${HPC_SCRATCH_DIR}/" "${NAS_TARGET_DIR}/" || rc=1
  elif [[ ${rc} -eq 0 ]]; then
    rc=1
  fi
  if [[ ${rc} -eq 0 ]]; then
    while IFS=$'\t' read -r ds view path output local_digest local_size; do
      remote_dir="${NAS_TARGET_DIR}/${ds}/output"
      ecoda_compare_checksum_remote "${path}" \
        "${remote_dir}/${output}" "${remote_dir}/${output}.md5" \
        "${local_digest}" "${local_size}" || rc=1
    done < "${verify_manifest}"
  fi
  rm -f "${transfer_manifest}" "${verify_manifest}"
  rmdir "${sync_lock}" 2>/dev/null || rc=1
  return "${rc}"
}

build_recovery_selection() {
  local target="$1" ds view
  : > "${target}" || return 1
  if [[ -n "${SELECTION_FILE_ARG}" ]]; then
    cp "${SELECTION_FILE_ARG}" "${target}" || return 1
  else
    [[ -n "${DATASETS_ARG}" && -n "${VIEWS_ARG}" ]] || {
      echo "ERROR: numeric --sync-only requires --datasets and --view/--views." >&2
      return 1
    }
    ecoda_split_csv "${DATASETS_ARG}" || return 1
    DATASET_NAMES=("${ECODA_ARRAY[@]}")
    ecoda_assert_unique_items "${DATASET_NAMES[@]}" || return 1
    ecoda_split_csv "${VIEWS_ARG}" || return 1
    RECOVERY_VIEWS=("${ECODA_ARRAY[@]}")
    ecoda_assert_unique_items "${RECOVERY_VIEWS[@]}" || return 1
    for ds in "${DATASET_NAMES[@]}"; do
      ecoda_dataset_exists "${ds}" || return 1
      for view in "${RECOVERY_VIEWS[@]}"; do
        ecoda_view_exists "${ds}" "${view}" || return 1
        printf '%s\t%s\n' "${ds}" "${view}" >> "${target}" || return 1
      done
    done
  fi
  ecoda_validate_manifest "${target}" 2 || return 1
  if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
    ecoda_validate_exact_batch_selection "${target}" 2 || return 1
  fi
}

gate_recovery_scheduler_id() {
  local scheduler_id="$1" expected="$2" rows
  [[ "${scheduler_id}" =~ ^[0-9]+$ ]] || return 1
  ECODA_ACCOUNTING_ROWS=""
  if ecoda_wait_array_accounting "${scheduler_id}" "${expected}" \
      "${STAGE3_WATCHDOG_POLL_SECONDS:-30}"; then
    return 0
  fi
  rows="${ECODA_ACCOUNTING_ROWS:-}"
  if [[ "${rows}" == *"${scheduler_id}_"* ]]; then
    return 1
  fi
  ecoda_wait_scalar_accounting "${scheduler_id}" "${STAGE3_WATCHDOG_POLL_SECONDS:-30}"
}

numeric_sync_only() {
  local ids="$1" recovery_manifest expected scheduler_id path ds view failed=0
  recovery_manifest="${TMPDIR:-/tmp}/ecoda_stage3_sync_${$}.tsv"
  build_recovery_selection "${recovery_manifest}" || {
    rm -f "${recovery_manifest}"
    return 1
  }
  expected="$(wc -l < "${recovery_manifest}" | tr -d '[:space:]')"
  ecoda_split_csv "${ids}" || { rm -f "${recovery_manifest}"; return 1; }
  for scheduler_id in "${ECODA_ARRAY[@]}"; do
    gate_recovery_scheduler_id "${scheduler_id}" "${expected}" || {
      echo "ERROR: scheduler recovery gate failed for ${scheduler_id}." >&2
      rm -f "${recovery_manifest}"
      return 1
    }
  done
  while IFS=$'\t' read -r ds view; do
    path="$(output_path_for "${ds}" "${view}")" || failed=1
    validate_h5ad "${ds}" "${view}" "${path}" || failed=1
    ecoda_validate_checksum "${path}" || failed=1
  done < "${recovery_manifest}"
  if [[ ${failed} -ne 0 ]]; then
    rm -f "${recovery_manifest}"
    return 1
  fi
  if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" != "1" ]]; then
    SYNC_LOCK_ROOT="${TMPDIR:-/tmp}/ecoda_stage3_sync_${$}"
    export SYNC_LOCK_ROOT
    sync_selected "${recovery_manifest}" || {
      rm -f "${recovery_manifest}"
      return 1
    }
  fi
  rm -f "${recovery_manifest}"
  printf 'PREPROCESS_SYNC_ONLY_IDS=%s\n' "${ids}"
}

if [[ -n "${SYNC_ONLY_RUN}" &&
      "${SYNC_ONLY_RUN}" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
  numeric_sync_only "${SYNC_ONLY_RUN}" || exit 1
  exit 0
fi

if [[ -n "${SYNC_ONLY_RUN}" ]]; then
  ecoda_open_run "${SYNC_ONLY_RUN}" || exit 1
  RUN_ID="${SYNC_ONLY_RUN}"
  export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT
  set +e
  stage3_load_bound_run
  bound_rc=$?
  set -e
  if [[ ${bound_rc} -eq 2 ]]; then
    stage3_abort "legacy_source_unpinned"
  fi
  [[ ${bound_rc} -eq 0 ]] ||
    stage3_abort "Stage 3 run-bound source/runtime identity is invalid"
  ecoda_runtime_validate_bound_run ||
    stage3_abort "Stage 3 run-bound runtime validation failed"
  export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST_ORIGINAL}"
  RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage3 0)" ||
    stage3_abort "Stage 3 run-bound runtime export construction failed"
  RUNTIME_EXPORT="${RUNTIME_EXPORT},ECODA_RUNTIME_IDENTITY=${RUNTIME_IDENTITY},ECODA_SOURCE_MANIFEST_RUN=${SOURCE_MANIFEST_RUN},ECODA_RUNTIME_RUN_ID=${RUN_ID}"
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  PENDING_MANIFEST="${ECODA_RUN_ROOT}/manifests/pending.tsv"
  SCHEDULER_IDS_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
  STATUS_FILE="${ECODA_RUN_ROOT}/status/watchdog"
  ecoda_validate_run_owned_path "${MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage3_abort "Stage 3 selection manifest is not run-owned"
  ecoda_validate_manifest "${MANIFEST}" 2 ||
    stage3_abort "Stage 3 selection manifest is invalid"
  ecoda_validate_checksum "${MANIFEST}" ||
    stage3_abort "Stage 3 selection checksum is invalid"
  stage3_classify_selection "${MANIFEST}" "${DATASETS_JSON_FILE}" ||
    stage3_abort "Stage 3 selection classification is invalid"
  case "${STAGE3_SELECTION_CLASSIFICATION}" in
    uncorrected_four|combined)
      if [[ -n "${STAGE3_COVID_PREFLIGHT_ROOT_REQUESTED}" &&
            "${STAGE3_COVID_PREFLIGHT_ROOT_REQUESTED}" != "${ECODA_RUN_ROOT}/preflight" ]]; then
        stage3_abort "Stage 3 sync-only Covid preflight root override is not run-bound"
      fi
      STAGE3_COVID_PREFLIGHT_ROOT="${ECODA_RUN_ROOT}/preflight"
      stage3_validate_covid_obs_reports ||
        stage3_abort "Stage 3 sync-only Covid obs preflight validation failed"
      ;;
  esac
  stage3_validate_corrected_source_evidence "${MANIFEST}" ||
    stage3_abort "Stage 3 corrected source preflight evidence is missing or invalid"
  [[ -r "${PENDING_MANIFEST}" ]] ||
    stage3_abort "Stage 3 pending manifest is missing"
  ecoda_validate_run_owned_path "${PENDING_MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage3_abort "Stage 3 pending manifest is not run-owned"
  if [[ -s "${PENDING_MANIFEST}" ]]; then
    ecoda_validate_manifest "${PENDING_MANIFEST}" 2 ||
      stage3_abort "Stage 3 pending manifest is invalid"
    require_scheduler_jobs=1
  else
    require_scheduler_jobs=0
  fi
  stage3_validate_scheduler_manifest "${SCHEDULER_IDS_FILE}" \
    "${require_scheduler_jobs}" ||
    stage3_abort "Stage 3 scheduler ID manifest is missing or invalid"
  if [[ "${require_scheduler_jobs}" -eq 1 ]]; then
    stage3_watchdog_terminal_ok ||
      stage3_abort "Stage 3 watchdog has not reached terminal success"
    stage3_record_watchdog_status_ids ||
      stage3_abort "Stage 3 watchdog scheduler records are invalid"
    stage3_validate_scheduler_manifest "${SCHEDULER_IDS_FILE}" 1 ||
      stage3_abort "Stage 3 scheduler ID manifest is incomplete"
    if [[ -s "${STATUS_FILE}" ]]; then
      grep -q '^STATE=OK$' "${STATUS_FILE}"
    else
      [[ -s "${ECODA_RUN_ROOT}/status/terminal" ]] &&
        grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/terminal"
    fi ||
      stage3_abort "all-skipped Stage 3 run lacks terminal success"
  fi
  failed=0
  while IFS=$'\t' read -r ds view; do
    [[ -n "${ds}" && -n "${view}" ]] || { failed=1; continue; }
    path="$(output_path_for "${ds}" "${view}")" || { failed=1; continue; }
    validate_h5ad "${ds}" "${view}" "${path}" || failed=1
    stage3_artifact_record_any "${path}" || failed=1
  done < "${MANIFEST}"
  [[ ${failed} -eq 0 ]] ||
    stage3_abort "Stage 3 sync-only h5ad contract/checksum failed"
  if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" != "1" ]]; then
    sync_selected "${MANIFEST}" ||
      stage3_abort "selected Stage 3 sync failed"
  fi
  stage3_finalize_owner_manifest OK "Stage 3 sync-only completed" ||
    stage3_abort "failed to finalize Stage 3 owners after sync"
  ecoda_set_run_state OK "sync-only Stage 3 validation and selected sync passed" ||
    stage3_abort "failed to write Stage 3 terminal OK state"
  echo "PREPROCESS_RUN_ID=${SYNC_ONLY_RUN}"
  exit 0
fi

DATASET_NAMES=()
if [[ -n "${SELECTION_FILE_ARG}" ]]; then
  [[ -r "${SELECTION_FILE_ARG}" ]] || { echo "ERROR: selection file is unreadable: ${SELECTION_FILE_ARG}" >&2; exit 1; }
elif [[ -n "${DATASETS_ARG}" ]]; then
  ecoda_split_csv "${DATASETS_ARG}"
  DATASET_NAMES=("${ECODA_ARRAY[@]}")
  ecoda_assert_unique_items "${DATASET_NAMES[@]}"
else
  while IFS= read -r ds; do DATASET_NAMES+=("${ds}"); done < <(jq -r 'keys[] | select(startswith("_") | not)' "${DATASETS_JSON_FILE}")
fi

if [[ -z "${SELECTION_FILE_ARG}" && ${#DATASET_NAMES[@]} -eq 0 ]]; then
  echo "ERROR: no Stage 3 datasets selected." >&2
  exit 1
fi
if [[ -z "${SELECTION_FILE_ARG}" ]]; then
  for ds in "${DATASET_NAMES[@]}"; do
    ecoda_dataset_exists "${ds}" || { echo "ERROR: unknown dataset '${ds}'." >&2; exit 1; }
  done
fi
stage3_build_generated_selection() {
  local target="$1" ds view
  : > "${target}" || return 1
  for ds in "${DATASET_NAMES[@]}"; do
    ecoda_dataset_exists "${ds}" || return 1
    if [[ -n "${VIEWS_ARG}" ]]; then
      ecoda_split_csv "${VIEWS_ARG}" || return 1
      for view in "${ECODA_ARRAY[@]}"; do
        ecoda_view_exists "${ds}" "${view}" || return 1
        [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" &&
           -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] || return 1
        printf '%s\t%s\n' "${ds}" "${view}" >> "${target}" || return 1
      done
    else
      while IFS= read -r view; do
        [[ -n "${view}" ]] || continue
        [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" &&
           -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] || return 1
        printf '%s\t%s\n' "${ds}" "${view}" >> "${target}" || return 1
      done < <(jq -r --arg ds "${ds}" '.[$ds].views // {} | keys[]' "${DATASETS_JSON_FILE}")
    fi
  done
  ecoda_validate_manifest "${target}" 2
}

if [[ -z "${SELECTION_FILE_ARG}" ]]; then
  PREVALIDATION_SELECTION="${TMPDIR:-/tmp}/ecoda_stage3_selection_${$}.tsv"
  stage3_build_generated_selection "${PREVALIDATION_SELECTION}" || {
    rm -f "${PREVALIDATION_SELECTION}"
    echo "ERROR: generated Stage 3 selection is malformed." >&2
    exit 1
  }
  stage3_classify_selection "${PREVALIDATION_SELECTION}" "${DATASETS_JSON_FILE}" || {
    rm -f "${PREVALIDATION_SELECTION}"
    echo "ERROR: generated Stage 3 selection classification is invalid." >&2
    exit 1
  }
  rm -f "${PREVALIDATION_SELECTION}"
fi

if [[ -n "${SELECTION_FILE_ARG}" ]]; then
  validate_external_selection "${SELECTION_FILE_ARG}" || {
    echo "ERROR: Stage 3 selection file is malformed or semantically invalid." >&2
    exit 1
  }
  stage3_classify_selection "${SELECTION_FILE_ARG}" "${DATASETS_JSON_FILE}" || {
    echo "ERROR: Stage 3 selection classification is invalid." >&2
    exit 1
  }
fi

stage3_require_new_snapshot || {
  echo "ERROR: legacy_source_unpinned" >&2
  exit 1
}
if [[ "${ECODA_RUNTIME_MODE_REQUESTED}" == "host" ]]; then
  echo "ERROR: host_snapshot_requires_apptainer" >&2
  exit 1
fi
RUN_ID="${ECODA_RUN_ID:-$(ecoda_new_run_id stage3)}"
ecoda_init_run stage3 "${RUN_ID}" >/dev/null || {
  echo "ERROR: Stage 3 run root already exists or could not be initialized." >&2
  exit 1
}
export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT
SOURCE_ROOT="${ECODA_SOURCE_ROOT}"
SOURCE_MANIFEST_ORIGINAL="${ECODA_SOURCE_MANIFEST}"
SOURCE_MANIFEST_RUN="${ECODA_RUN_ROOT}/manifests/source.manifest"
RUNTIME_IDENTITY="${ECODA_RUN_ROOT}/manifests/runtime.identity"
stage3_install_source_manifest ||
  stage3_abort "failed to copy immutable Stage 3 source manifest"
if [[ -z "${ECODA_RUNTIME_MODE_REQUESTED}" ]]; then
  export ECODA_RUNTIME_MODE=apptainer
else
  export ECODA_RUNTIME_MODE="${ECODA_RUNTIME_MODE_REQUESTED}"
fi
MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
MANIFEST_TMP="${MANIFEST}.build.$$"
: > "${MANIFEST_TMP}"
echo "PREPROCESS_RUN_ID=${RUN_ID}"
echo "PREPROCESS_DATASET_MANIFEST=${MANIFEST}"

append_selection() {
  local ds="$1" view
  ecoda_dataset_exists "${ds}" || { echo "ERROR: unknown dataset '${ds}'." >&2; return 1; }
  if [[ -n "${VIEWS_ARG}" ]]; then
    ecoda_split_csv "${VIEWS_ARG}"
    for view in "${ECODA_ARRAY[@]}"; do
      ecoda_view_exists "${ds}" "${view}" || { echo "ERROR: ${ds}/${view} is not declared in datasets.json." >&2; return 1; }
      [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" && -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] || { echo "ERROR: ${ds}/${view} has no input/output contract." >&2; return 1; }
      printf '%s\t%s\n' "${ds}" "${view}" >> "${MANIFEST_TMP}"
    done
  else
    while IFS= read -r view; do
      [[ -n "${view}" ]] || continue
      [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" && -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] || { echo "ERROR: ${ds}/${view} has no input/output contract." >&2; return 1; }
      printf '%s\t%s\n' "${ds}" "${view}" >> "${MANIFEST_TMP}"
    done < <(jq -r --arg ds "${ds}" '.[$ds].views // {} | keys[]' "${DATASETS_JSON_FILE}")
  fi
}

if [[ -n "${SELECTION_FILE_ARG}" ]]; then
  cp "${SELECTION_FILE_ARG}" "${MANIFEST_TMP}" ||
    stage3_abort "failed to copy Stage 3 selection file"
  ecoda_validate_manifest "${MANIFEST_TMP}" 2 ||
    stage3_abort "malformed Stage 3 selection file"
else
  for ds in "${DATASET_NAMES[@]}"; do
    append_selection "${ds}" || stage3_abort "invalid Stage 3 selection"
  done
  ecoda_validate_manifest "${MANIFEST_TMP}" 2 ||
    stage3_abort "empty Stage 3 selection"
fi
if ! ecoda_atomic_install_manifest "${MANIFEST_TMP}" "${MANIFEST}" 2; then
  stage3_abort "failed to install Stage 3 selection atomically"
fi
rm -f "${MANIFEST_TMP}"
ecoda_write_checksum "${MANIFEST}" || stage3_abort "failed to checksum Stage 3 selection"

# Reject duplicate dataset/view owners without relying on JSON key order.
SEEN_ROWS=""
while IFS=$'\t' read -r ds view; do
  row="${ds}/${view}"
  case " ${SEEN_ROWS} " in
    *" ${row} "*) stage3_abort "duplicate selection row ${row}" ;;
  esac
  SEEN_ROWS="${SEEN_ROWS} ${row}"
done < "${MANIFEST}"
export ECODA_RUNTIME_PROFILE=stage3
ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE:-apptainer}" ||
  stage3_abort "Stage 3 immutable runtime validation failed"
[[ -s "${RUNTIME_IDENTITY}" && ! -L "${RUNTIME_IDENTITY}" ]] ||
  stage3_abort "Stage 3 runtime submission did not write runtime.identity"
stage3_load_bound_run ||
  stage3_abort "Stage 3 run-bound source/runtime identity is invalid"
stage3_record_identity_metadata ||
  stage3_abort "failed to record Stage 3 source/runtime identity"
stage3_classify_selection "${MANIFEST}" "${DATASETS_JSON_FILE}" ||
  stage3_abort "Stage 3 selection classification is invalid"
stage3_validate_corrected_source_release "${MANIFEST}" "${DATASETS_JSON_FILE}" ||
  stage3_abort "Stage 3 corrected source metadata/configuration release failed"
stage3_validate_corrected_source_evidence "${MANIFEST}" ||
  stage3_abort "Stage 3 corrected source preflight evidence is missing or invalid"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage3 0)" ||
  stage3_abort "Stage 3 runtime export construction failed"
RUNTIME_EXPORT="${RUNTIME_EXPORT},ECODA_RUNTIME_IDENTITY=${RUNTIME_IDENTITY},ECODA_SOURCE_MANIFEST_RUN=${SOURCE_MANIFEST_RUN}"
SCHEDULER_IDS_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
ecoda_atomic_write "${SCHEDULER_IDS_FILE}" "" ||
  stage3_abort "failed to initialize Stage 3 scheduler manifest"
stage3_run_covid_obs_preflight "${MANIFEST}" ||
  stage3_abort "Covid obs-only preflight failed; Stage 3 release is blocked"
stage3_compute_validate_existing ||
  stage3_abort "Stage 3 compute-node existing-output preflight failed"

PENDING_MANIFEST="${ECODA_RUN_ROOT}/manifests/pending.tsv"
PENDING_TMP="${PENDING_MANIFEST}.build.$$"
OWNERS_MANIFEST="${ECODA_RUN_ROOT}/manifests/owners.tsv"
OWNERS_TMP="${OWNERS_MANIFEST}.build.$$"
PENDING_COUNT=0
ecoda_owner_clear_tracked
while IFS=$'\t' read -r ds view; do
  path="$(output_path_for "${ds}" "${view}")" ||
    stage3_abort "missing output contract for ${ds}/${view}"
  valid=0
  if [[ ${FORCE_ARG} -eq 0 && -s "${path}" ]]; then
    set +e
    stage3_existing_output_valid "${ds}" "${view}" "${path}"
    preflight_rc=$?
    set -e
    if [[ ${preflight_rc} -eq 0 ]]; then
      valid=1
    elif [[ ${preflight_rc} -ne 1 ]]; then
      stage3_abort "missing or malformed Stage 3 preflight status for ${ds}/${view}"
    fi
  fi
  if [[ ${valid} -eq 1 ]]; then
    echo "Skipping validated Stage 3 artifact ${ds}/${view}."
    continue
  fi
  set +e
  owner_dir="$(ecoda_owner_acquire stage3 "${ds}/${view}" "${RUN_ID}" "${FORCE_ARG}" 0)"
  owner_rc=$?
  set -e
  if [[ ${owner_rc} -ne 0 ]]; then
    stage3_abort "ownership conflict for ${ds}/${view}"
  fi
  ecoda_owner_track "${owner_dir}" ||
    stage3_abort "failed to track Stage 3 owner for ${ds}/${view}"
  if [[ ${FORCE_ARG} -eq 1 ]]; then
    ecoda_invalidate_artifact "${path}" ||
      stage3_abort "failed to invalidate Stage 3 artifact ${path}"
  fi
  printf '%s\t%s\n' "${ds}" "${view}" >> "${PENDING_TMP}"
  printf '%s\t%s\n' "${ds}/${view}" "${owner_dir}" >> "${OWNERS_TMP}"
  PENDING_COUNT=$((PENDING_COUNT + 1))
done < "${MANIFEST}"

if [[ ${PENDING_COUNT} -gt 0 ]]; then
  ecoda_atomic_install_manifest "${PENDING_TMP}" "${PENDING_MANIFEST}" 2 ||
    stage3_abort "failed to install Stage 3 pending manifest atomically"
  ecoda_atomic_install_manifest "${OWNERS_TMP}" "${OWNERS_MANIFEST}" 2 ||
    stage3_abort "failed to install Stage 3 owner manifest atomically"
else
  ecoda_atomic_write "${PENDING_MANIFEST}" "" ||
    stage3_abort "failed to create empty Stage 3 pending manifest"
  ecoda_atomic_write "${OWNERS_MANIFEST}" "" ||
    stage3_abort "failed to create empty Stage 3 owner manifest"
fi
rm -f "${PENDING_TMP}" "${OWNERS_TMP}"

if [[ ${PENDING_COUNT} -eq 0 ]]; then
  ecoda_validate_output_ownership stage3 "${MANIFEST}" "${RUN_ID}" 1 ||
    stage3_abort "Stage 3 no-op output ownership validation failed"
  if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" != "1" ]] &&
     ! sync_selected "${MANIFEST}"; then
    stage3_abort "selected Stage 3 sync failed"
  fi
  ecoda_owner_finalize_tracked OK "Stage 3 no-op artifacts validated" ||
    stage3_abort "failed to finalize Stage 3 no-op output owners"

  NOOP_REPORT="${ECODA_RUN_REPORT:-${ECODA_RUN_ROOT}/status/noop}"
  case "${NOOP_REPORT}" in
    "${ECODA_RUN_ROOT}"/*) ;;
    *) NOOP_REPORT="${ECODA_RUN_ROOT}/status/noop" ;;
  esac
  ecoda_atomic_write "${NOOP_REPORT}" \
    "STATE=NOOP_VALIDATED\nRUN_ID=${RUN_ID}\nREASON=all selected Stage 3 artifacts already validated\n" ||
    stage3_abort "failed to write Stage 3 validator-only report"
  ecoda_set_run_state NOOP_VALIDATED "all selected Stage 3 artifacts already validated" ||
    stage3_abort "failed to write Stage 3 NOOP_VALIDATED state"
  echo "NOOP_VALIDATED=${RUN_ID}"
  echo "PREPROCESS_RUN_ID=${RUN_ID}"
  exit 0
fi

ecoda_validate_output_ownership stage3 "${PENDING_MANIFEST}" "${RUN_ID}" 1 ||
  stage3_abort "Stage 3 output ownership validation failed before array submission"
mkdir -p "${LOGS_DIR}" || stage3_abort "failed to create Stage 3 log directory"
export FORCE_PREPROCESS="${FORCE_ARG}"
export PREPROCESS_SELECTION_FILE="${PENDING_MANIFEST}"
export PREPROCESS_RUN_ROOT="${ECODA_RUN_ROOT}"
export PREPROCESS_ERROR_PREFIX="${LOGS_DIR}/3_scrnaseq_preprocessing"
worker_script="$(stage3_require_source_script "${SCRIPT_DIR}/1.1_run_worker.sh")" ||
  stage3_abort "Stage 3 worker script escaped immutable source root"
set +e
ARRAY_MSG="$(sbatch --parsable --array="1-${PENDING_COUNT}%${THROTTLE}" \
  --mem="${MEMORY}" --partition="${PARTITION}" \
  --output="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.log" \
  --error="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.err" \
  --mail-user="${USER_EMAIL}" \
  --export="ALL,PREPROCESS_SELECTION_FILE=${PENDING_MANIFEST},PREPROCESS_RUN_ROOT=${ECODA_RUN_ROOT},FORCE_PREPROCESS=${FORCE_ARG},PREPROCESS_ERROR_PREFIX=${LOGS_DIR}/3_scrnaseq_preprocessing,${RUNTIME_EXPORT}" \
  "${worker_script}")"
array_rc=$?
set -e
[[ ${array_rc} -eq 0 ]] || stage3_abort "sbatch rejected Stage 3 preprocessing array"
ARRAY_ID="${ARRAY_MSG%%;*}"
[[ "${ARRAY_ID}" =~ ^[0-9]+$ ]] || stage3_abort "invalid Stage 3 array id"
echo "PREPROCESS_ARRAY_JOB_ID=${ARRAY_ID}"
stage3_install_scheduler_record ARRAY "${ARRAY_ID}" ||
  stage3_abort "failed to persist Stage 3 array scheduler ID"
stage3_validate_scheduler_manifest "${SCHEDULER_IDS_FILE}" 0 ||
  stage3_abort "Stage 3 scheduler ID manifest is invalid after array submission"
watchdog_script="$(stage3_require_source_script "${SCRIPT_DIR}/1.2_preprocess_watchdog.sh")" ||
  stage3_abort "Stage 3 watchdog script escaped immutable source root"
set +e
WATCHDOG_MSG="$(sbatch --parsable --wait --dependency="afterany:${ARRAY_ID}" \
  --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=2G \
  --time="${STAGE3_WATCHDOG_TIME_LIMIT:-12:00:00}" \
  --output="${LOGS_DIR}/3_scrnaseq_preprocessing_watchdog_%j.log" \
  --error="${LOGS_DIR}/3_scrnaseq_preprocessing_watchdog_%j.err" \
  --mail-user="${USER_EMAIL}" \
  --export="ALL,PREPROCESS_RUN_ROOT=${ECODA_RUN_ROOT},PREPROCESS_PENDING_MANIFEST=${PENDING_MANIFEST},${RUNTIME_EXPORT}" \
  "${watchdog_script}" "${RUN_ID}" "${MANIFEST}" \
  "${ARRAY_ID}" "${MEMORY}" "${MAX_MEMORY}" "${PARTITION}" "${THROTTLE}")"
watchdog_rc=$?
set -e
WATCHDOG_ID="${WATCHDOG_MSG%%;*}"
[[ "${WATCHDOG_ID}" =~ ^[0-9]+$ ]] ||
  stage3_abort "invalid Stage 3 watchdog id"
echo "PREPROCESS_WATCHDOG_JOB_ID=${WATCHDOG_ID}"
stage3_install_scheduler_record WATCHDOG "${WATCHDOG_ID}" ||
  stage3_abort "failed to persist Stage 3 watchdog scheduler ID"
if [[ ${watchdog_rc} -ne 0 ]]; then
  stage3_record_watchdog_status_ids ||
    stage3_abort "failed to preserve Stage 3 watchdog scheduler IDs"
  stage3_abort "Stage 3 watchdog job failed"
fi

if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" == "1" ]]; then
  ecoda_set_run_state OK "submitter test mode; scheduler calls validated" ||
    stage3_abort "failed to write Stage 3 terminal OK state"
  exit 0
fi
if [[ ! -s "${ECODA_RUN_ROOT}/status/watchdog" ]] ||
   ! grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/watchdog"; then
  stage3_abort "Stage 3 watchdog did not report OK"
fi
while IFS= read -r status_line; do
  case "${status_line}" in
    SCHEDULER_ID=*|ARRAY_JOB_ID=*)
      printf 'PREPROCESS_SCHEDULER_ID=%s:%s\n' "${RUN_ID}" "${status_line#*=}"
      ;;
  esac
done < "${ECODA_RUN_ROOT}/status/watchdog"
printf 'PREPROCESS_SCHEDULER_ID=%s:%s\n' "${RUN_ID}" "${WATCHDOG_ID}"
stage3_record_watchdog_status_ids ||
  stage3_abort "failed to install complete Stage 3 scheduler manifest"
stage3_validate_scheduler_manifest "${SCHEDULER_IDS_FILE}" 1 ||
  stage3_abort "Stage 3 scheduler manifest is incomplete"
if ! sync_selected "${MANIFEST}"; then
  stage3_abort "selected Stage 3 sync failed"
fi
stage3_finalize_owner_manifest OK "Stage 3 sync completed" ||
  stage3_abort "failed to finalize Stage 3 owners"
ecoda_set_run_state OK "Stage 3 preprocessing, validation, and selected sync completed" ||
  stage3_abort "failed to write Stage 3 terminal OK state"
