#!/bin/bash
# Read-only H5AD obs audit/metadata export.  This worker never writes to the
# H5AD source or creates an artifact ownership record for it.
set -euo pipefail

BOOTSTRAP_SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
BOOTSTRAP_SOURCE_VALID=0
if [[ -n "${BOOTSTRAP_SOURCE_ROOT}" &&
      "${BOOTSTRAP_SOURCE_ROOT}" = /* &&
      "${BOOTSTRAP_SOURCE_ROOT}" != *$'\n'* &&
      "${BOOTSTRAP_SOURCE_ROOT}" != *$'\t'* &&
      "${BOOTSTRAP_SOURCE_ROOT##*/}" == "tree" ]]; then
  BOOTSTRAP_SNAPSHOT_ROOT="${BOOTSTRAP_SOURCE_ROOT%/tree}"
  if [[ "${BOOTSTRAP_SNAPSHOT_ROOT##*/}" =~ ^[[:xdigit:]]{40}$ &&
        -d "${BOOTSTRAP_SOURCE_ROOT}" &&
        ! -L "${BOOTSTRAP_SOURCE_ROOT}" &&
        -f "${BOOTSTRAP_SOURCE_ROOT}/src/slurm_config.sh" &&
        ! -L "${BOOTSTRAP_SOURCE_ROOT}/src/slurm_config.sh" &&
        -r "${BOOTSTRAP_SOURCE_ROOT}/src/slurm_config.sh" ]]; then
    BOOTSTRAP_SOURCE_VALID=1
  fi
fi

if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
  [[ "${BOOTSTRAP_SOURCE_VALID}" == "1" ]] || {
    echo "ERROR: container audit worker requires a valid immutable ECODA_SOURCE_ROOT." >&2
    exit 1
  }
  SCRIPT_DIR="${BOOTSTRAP_SOURCE_ROOT%/}/src/utils/bash"
elif [[ "${BOOTSTRAP_SOURCE_VALID}" == "1" ]]; then
  SCRIPT_DIR="${BOOTSTRAP_SOURCE_ROOT%/}/src/utils/bash"
elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
  command -v scontrol >/dev/null 2>&1 || {
    echo "ERROR: scontrol is required to recover the immutable worker path." >&2
    exit 1
  }
  SCRIPT_DIR="$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/ {print $2}' | xargs dirname)"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
[[ -n "${SCRIPT_DIR}" ]] || {
  echo "ERROR: could not recover the immutable worker directory." >&2
  exit 1
}
SLURM_CONFIG_FILE="${SCRIPT_DIR}/../../slurm_config.sh"
[[ -f "${SLURM_CONFIG_FILE}" &&
   ! -L "${SLURM_CONFIG_FILE}" &&
   -r "${SLURM_CONFIG_FILE}" ]] || {
  echo "ERROR: immutable Slurm config is missing or unsafe: ${SLURM_CONFIG_FILE}" >&2
  exit 1
}
source "${SLURM_CONFIG_FILE}"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"

AUDIT_MODE="${H5AD_OBS_AUDIT_MODE:-audit}"
case "${AUDIT_MODE}" in
  audit) AUDIT_RUNTIME_STAGE=stage3 ;;
  metadata) AUDIT_RUNTIME_STAGE=stage5 ;;
  *) echo "ERROR: unsupported obs worker mode: ${AUDIT_MODE}" >&2; exit 1 ;;
esac
SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
SOURCE_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
RUN_ROOT="${H5AD_OBS_AUDIT_RUN_ROOT:-${ECODA_RUN_ROOT:-}}"
RUN_ID="${H5AD_OBS_AUDIT_RUN_ID:-${ECODA_RUN_ID:-${RUN_ROOT##*/}}}"
fail() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ "${SOURCE_REQUIRED}" == "1" ]] ||
  fail "read-only H5AD obs audit requires an immutable source snapshot"
[[ "${SOURCE_ROOT}" = /* && "${SOURCE_MANIFEST}" = /* ]] ||
  fail "immutable source root and manifest are required"
[[ "${SOURCE_ROOT##*/}" == "tree" ]] ||
  fail "immutable source root is not a snapshot tree"
SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
[[ "${SNAPSHOT_ROOT##*/}" =~ ^[[:xdigit:]]{40}$ ]] ||
  fail "immutable source root is not commit keyed"
EXPECTED_SOURCE_MANIFEST="${SNAPSHOT_ROOT}/identity/source.manifest"
[[ "${SOURCE_MANIFEST}" == "${EXPECTED_SOURCE_MANIFEST}" ]] ||
  fail "source manifest is not bound to immutable source tree"
[[ -d "${SOURCE_ROOT}" && ! -L "${SOURCE_ROOT}" &&
   -f "${SOURCE_MANIFEST}" && ! -L "${SOURCE_MANIFEST}" &&
   -r "${SOURCE_MANIFEST}" ]] ||
  fail "immutable source identity is missing or unsafe"

[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${RUN_ROOT}" ]] ||
  fail "inherited run root does not match audit run root"
[[ -z "${ECODA_RUN_ID:-}" || "${ECODA_RUN_ID}" == "${RUN_ID}" ]] ||
  fail "inherited run ID does not match audit run ID"
[[ -z "${H5AD_OBS_AUDIT_RUN_ID:-}" ||
   "${H5AD_OBS_AUDIT_RUN_ID}" == "${RUN_ID}" ]] ||
  fail "obs audit scheduler run ID does not match audit run ID"
[[ "${RUN_ROOT}" = /* && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] ||
  fail "read-only H5AD obs audit requires an existing run root"
[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ &&
   "${RUN_ROOT##*/}" == "${RUN_ID}" ]] ||
  fail "read-only H5AD obs audit run identity is invalid"
[[ "${ECODA_RUNS_ROOT:-}" = /* &&
   "${ECODA_RUNS_ROOT}" != *$'\n'* &&
   "${ECODA_RUNS_ROOT}" != *$'\t'* &&
   -d "${ECODA_RUNS_ROOT}" && ! -L "${ECODA_RUNS_ROOT}" ]] ||
  fail "global ECODA run root is missing or unsafe"
EXPECTED_RUN_ROOT="${ECODA_RUNS_ROOT%/}/${RUN_ID}"
[[ "${RUN_ROOT}" == "${EXPECTED_RUN_ROOT}" ]] ||
  fail "audit run root is not the global run root for its ID"
_ecoda_validate_path_ancestors "${RUN_ROOT}" "${ECODA_RUNS_ROOT}" ||
  fail "audit run root has a symlinked ancestor"
RUNS_ROOT_REAL="$(realpath -e "${ECODA_RUNS_ROOT}" 2>/dev/null ||
  realpath "${ECODA_RUNS_ROOT}" 2>/dev/null)" ||
  fail "could not canonicalize global ECODA run root"
[[ "${RUNS_ROOT_REAL}" == "${ECODA_RUNS_ROOT}" ]] ||
  fail "global ECODA run root is not canonical"
RUN_ROOT_REAL="$(realpath -e "${RUN_ROOT}" 2>/dev/null || realpath "${RUN_ROOT}" 2>/dev/null)" ||
  fail "could not canonicalize audit run root"
[[ "${RUN_ROOT_REAL}" == "${RUN_ROOT}" ]] ||
  fail "audit run root is not canonical"
_ecoda_validate_run_owned_directory_for_create() {
  local candidate="${1:-}"
  local boundary="${2:-${RUN_ROOT}}"
  local boundary_real current parent candidate_real current_real
  [[ "${candidate}" = /* && "${candidate}" != *$'\n'* &&
     "${candidate}" != *$'\t'* ]] || return 1
  [[ "${boundary}" = /* && "${boundary}" != *$'\n'* &&
     "${boundary}" != *$'\t'* ]] || return 1
  case "${candidate}" in
    *"/.."|*"/../"*) return 1 ;;
  esac
  case "${boundary}" in
    *"/.."|*"/../"*) return 1 ;;
  esac
  case "${candidate}" in
    "${boundary}"|"${boundary}"/*) ;;
    *) return 1 ;;
  esac
  if [[ "${boundary}" == "${RUN_ROOT}" ]]; then
    boundary_real="${RUN_ROOT_REAL}"
  else
    boundary_real="$(realpath -e "${boundary}" 2>/dev/null ||
      realpath "${boundary}" 2>/dev/null)" || return 1
  fi
  [[ -n "${boundary_real}" ]] || return 1
  [[ ! -L "${candidate}" ]] || return 1
  _ecoda_validate_path_ancestors "${candidate}" "${boundary}" ||
    return 1

  current="${candidate}"
  while [[ ! -e "${current}" && ! -L "${current}" ]]; do
    parent="$(dirname "${current}")"
    [[ "${parent}" != "${current}" ]] || return 1
    current="${parent}"
  done
  [[ -d "${current}" && ! -L "${current}" ]] || return 1
  current_real="$(realpath -e "${current}" 2>/dev/null ||
    realpath "${current}" 2>/dev/null)" || return 1
  case "${current_real}" in
    "${boundary_real}"|"${boundary_real}"/*) ;;
    *) return 1 ;;
  esac

  if [[ -e "${candidate}" ]]; then
    [[ -d "${candidate}" ]] || return 1
    candidate_real="$(realpath -e "${candidate}" 2>/dev/null ||
      realpath "${candidate}" 2>/dev/null)" || return 1
    case "${candidate_real}" in
      "${boundary_real}"|"${boundary_real}"/*) ;;
      *) return 1 ;;
    esac
  fi
}


RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"

RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"
RUN_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
[[ -f "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" &&
   -r "${RUN_SOURCE_MANIFEST}" &&
   -f "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" &&
   -r "${RUN_RUNTIME_IDENTITY}" ]] ||
  fail "run-bound source/runtime identity manifests are missing or unsafe"
ecoda_validate_run_owned_path "${RUN_SOURCE_MANIFEST}" "${RUN_ROOT}" ||
  fail "run-bound source manifest is not run-owned"
ecoda_validate_run_owned_path "${RUN_RUNTIME_IDENTITY}" "${RUN_ROOT}" ||
  fail "run-bound runtime identity is not run-owned"
cmp -s "${RUN_SOURCE_MANIFEST}" "${SOURCE_MANIFEST}" ||
  fail "run-bound source manifest differs from immutable source manifest"

IDENTITY_IMAGE="$(_ecoda_runtime_require_identity_value \
  "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE)" ||
  fail "run runtime identity lacks RUNTIME_IMAGE"
IDENTITY_MANIFEST="$(_ecoda_runtime_require_identity_value \
  "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST)" ||
  fail "run runtime identity lacks RUNTIME_MANIFEST"
[[ "${IDENTITY_IMAGE}" = /* && "${IDENTITY_MANIFEST}" = /* &&
   "${ECODA_RUNTIME_IMAGE:-}" == "${IDENTITY_IMAGE}" &&
   "${ECODA_RUNTIME_MANIFEST:-}" == "${IDENTITY_MANIFEST}" ]] ||
  fail "runtime paths do not match run-bound runtime identity"
[[ -f "${IDENTITY_IMAGE}" && ! -L "${IDENTITY_IMAGE}" &&
   -f "${IDENTITY_MANIFEST}" && ! -L "${IDENTITY_MANIFEST}" ]] ||
  fail "run-bound runtime image or manifest is missing or unsafe"

export ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID="${RUN_ID}"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_IMAGE="${IDENTITY_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${IDENTITY_MANIFEST}"
export ECODA_RUNTIME_IDENTITY="${RUN_RUNTIME_IDENTITY}"
export ECODA_AUX_ROOT="${SOURCE_ROOT%/}/aux"
export ECODA_RUNTIME_PROFILE="${AUDIT_RUNTIME_STAGE}"
PROJECT_ROOT="${SOURCE_ROOT}"
DATASETS_JSON_FILE="${SOURCE_ROOT}/datasets.json"
export PROJECT_ROOT DATASETS_JSON_FILE

if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" &&
      "${ECODA_RUNTIME_MODE:-host}" == "host" ]]; then
  ecoda_runtime_validate_bound_run ||
    fail "run-bound runtime validation failed"
fi

WORKER_SCRIPT="${SOURCE_ROOT}/src/utils/bash/h5ad_obs_audit_worker.sh"
WORKER_SCRIPT="$(ecoda_require_source_script_path "${WORKER_SCRIPT}" "${SOURCE_ROOT}")" ||
  fail "obs audit worker escaped immutable source root"
ecoda_runtime_reexec_worker "${AUDIT_RUNTIME_STAGE}" "${WORKER_SCRIPT}" ||
  fail "could not enter immutable audit runtime"
if [[ "${AUDIT_MODE}" == metadata ]]; then
  METADATA_MANIFEST="${H5AD_METADATA_EXPORT_MANIFEST:-}"
  STATUS_DIR="${H5AD_METADATA_EXPORT_STATUS_DIR:-}"
  TASK_ID="${SLURM_ARRAY_TASK_ID:-${H5AD_METADATA_EXPORT_TASK_ID:-}}"
  if [[ -n "${H5AD_METADATA_EXPORT_VARIANT:-}" &&
        -n "${ANALYSIS_VARIANT:-}" &&
        "${H5AD_METADATA_EXPORT_VARIANT}" != "${ANALYSIS_VARIANT}" ]]; then
    fail "metadata export variant disagrees with Stage 5 analysis variant"
  fi
  if [[ -n "${H5AD_METADATA_EXPORT_PASS:-}" &&
        -n "${ANALYSIS_PASS:-}" &&
        "${H5AD_METADATA_EXPORT_PASS}" != "${ANALYSIS_PASS}" ]]; then
    fail "metadata export pass disagrees with Stage 5 analysis pass"
  fi
  METADATA_VARIANT="${H5AD_METADATA_EXPORT_VARIANT:-${ANALYSIS_VARIANT:-}}"
  METADATA_PASS="${H5AD_METADATA_EXPORT_PASS:-${ANALYSIS_PASS:-}}"
  METADATA_ROOT="${H5AD_METADATA_EXPORT_ROOT:-${ANALYSIS_ROOT:-}}"
  METADATA_NAS_ROOT="${H5AD_METADATA_EXPORT_NAS_ROOT:-${ANALYSIS_NAS_ROOT:-}}"
  case "${METADATA_VARIANT}" in
    final)
      EXPECTED_METADATA_PASS="uncorrected"
      EXPECTED_METADATA_ROOT="${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final"
      EXPECTED_METADATA_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/uncorrected_final"
      ;;
    corrected_final)
      EXPECTED_METADATA_PASS="corrected"
      if [[ "${METADATA_ROOT}" == */batch_effect/corrected_final/recovery_35row ||
            "${METADATA_NAS_ROOT}" == */batch_effect/corrected_final/recovery_35row ]]; then
        EXPECTED_METADATA_ROOT="${HPC_SCRATCH_DIR}/batch_effect/corrected_final/recovery_35row"
        EXPECTED_METADATA_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/corrected_final/recovery_35row"
      else
        EXPECTED_METADATA_ROOT="${HPC_SCRATCH_DIR}/batch_effect/corrected_final"
        EXPECTED_METADATA_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/corrected_final"
      fi
      ;;
    *)
      fail "metadata export requires analysis variant final or corrected_final"
      ;;
  esac
  [[ -z "${METADATA_PASS}" ||
     "${METADATA_PASS}" == "${EXPECTED_METADATA_PASS}" ]] ||
    fail "metadata export variant/pass identity mismatch"
  METADATA_PASS="${EXPECTED_METADATA_PASS}"
  [[ "${METADATA_ROOT}" == "${EXPECTED_METADATA_ROOT}" ]] ||
    fail "metadata export variant/root identity mismatch"
  if [[ -n "${METADATA_NAS_ROOT}" &&
        "${METADATA_NAS_ROOT}" != "${EXPECTED_METADATA_NAS_ROOT}" ]]; then
    fail "metadata export NAS root does not match analysis variant"
  fi
  export ANALYSIS_VARIANT="${METADATA_VARIANT}"
  export ANALYSIS_PASS="${METADATA_PASS}"
  export ANALYSIS_ROOT="${METADATA_ROOT}"
  [[ -n "${METADATA_MANIFEST}" && -r "${METADATA_MANIFEST}" &&
     ! -L "${METADATA_MANIFEST}" ]] ||
    fail "metadata export manifest is missing or unsafe"
  [[ -n "${STATUS_DIR}" && "${STATUS_DIR}" = /* &&
     ! -L "${STATUS_DIR}" &&
     "${TASK_ID}" =~ ^[0-9]+$ && ${TASK_ID} -gt 0 ]] ||
    fail "metadata export requires a run-owned status directory and task ID"
  ecoda_validate_run_owned_path "${METADATA_MANIFEST}" "${RUN_ROOT}" ||
    fail "metadata export manifest is not run-owned"
  ecoda_validate_checksum "${METADATA_MANIFEST}" ||
    fail "metadata export manifest checksum is invalid"
  ecoda_validate_manifest "${METADATA_MANIFEST}" 4 ||
    fail "metadata export manifest is malformed"
  row=""
  manifest_rows=0
  seen_metadata_rows=""
  manifest_line=""
  while IFS= read -r manifest_line || [[ -n "${manifest_line}" ]]; do
    IFS=$'\t' read -r row_dataset row_view row_input row_output row_extra <<< "${manifest_line}"
    manifest_rows=$((manifest_rows + 1))
    [[ -n "${row_dataset}" &&
       "${row_view}" == "batch_effect_${METADATA_PASS}" &&
       "${row_input}" = /* && -n "${row_output}" &&
       -z "${row_extra}" ]] ||
      fail "metadata export manifest row ${manifest_rows} is malformed"
    [[ "${row_dataset}" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ ]] ||
      fail "metadata export dataset is not a safe path component"
    ecoda_dataset_exists "${row_dataset}" ||
      fail "metadata export dataset is not configured: ${row_dataset}"
    ecoda_view_exists "${row_dataset}" "${row_view}" ||
      fail "metadata export view is not configured: ${row_dataset}/${row_view}"
    expected_input_name="$(ecoda_view_output_name "${row_dataset}" "${row_view}")" ||
      fail "metadata export output name is unavailable: ${row_dataset}/${row_view}"
    [[ -n "${expected_input_name}" &&
       "${expected_input_name}" != */* &&
       "${expected_input_name}" != *$'\n'* &&
       "${expected_input_name}" != *$'\t'* ]] ||
      fail "metadata export configured H5AD output name is unsafe"
    expected_input="${HPC_SCRATCH_DIR%/}/${row_dataset}/output/${expected_input_name}"
    [[ "${row_input}" == "${expected_input}" ]] ||
      fail "metadata export input is not the configured ${METADATA_PASS} H5AD path"
    case " ${seen_metadata_rows} " in
      *" ${row_dataset}|${row_view} "*)
        fail "metadata export manifest contains duplicate dataset/view" ;;
      *) seen_metadata_rows="${seen_metadata_rows} ${row_dataset}|${row_view}" ;;
    esac
    expected_output="${METADATA_ROOT%/}/metadata/${row_dataset}_sample_metadata.feather"
    [[ "${row_output}" == "${expected_output}" ]] ||
      fail "metadata export row ${manifest_rows} is not bound to the exact variant output"
    [[ -f "${row_input}" && ! -L "${row_input}" && -r "${row_input}" ]] ||
      fail "metadata export H5AD input is missing or unsafe: ${row_input}"
    if [[ ${manifest_rows} -eq ${TASK_ID} ]]; then
      row="${row_dataset}"$'\t'"${row_view}"$'\t'"${row_input}"$'\t'"${row_output}"
    fi
  done < "${METADATA_MANIFEST}"
  [[ ${manifest_rows} -gt 0 && ${TASK_ID} -le ${manifest_rows} ]] ||
    fail "metadata export task ID is outside manifest rows"
  [[ -n "${row}" ]] || fail "metadata export selected row is missing"
  IFS=$'\t' read -r row_dataset row_view row_input row_output row_extra <<< "${row}"
  metadata_output_parent="$(dirname "${row_output}")"
  _ecoda_validate_run_owned_directory_for_create \
    "${metadata_output_parent}" "${HPC_SCRATCH_DIR}" ||
    fail "metadata export output parent is not a safe final path"
  _ecoda_validate_run_owned_directory_for_create "${STATUS_DIR}" "${RUN_ROOT}" ||
    fail "metadata export status directory is not a safe run-owned path"
  [[ ! -L "${row_output}" && ! -L "${row_output}.md5" ]] ||
    fail "metadata export output is a symlink"
  mkdir -p "${STATUS_DIR}"
  ecoda_validate_run_owned_path "${STATUS_DIR}" "${RUN_ROOT}" ||
    fail "metadata export status directory is not run-owned"
  EXPORTER="${SOURCE_ROOT}/src/utils/py/export_h5ad_sample_metadata.py"
  EXPORTER="$(ecoda_require_source_script_path "${EXPORTER}" "${SOURCE_ROOT}")" ||
    fail "metadata exporter escaped immutable source root"
  mkdir -p "${metadata_output_parent}"
  status_kind="EXPORTED"
  if "${PYTHON_BIN}" "${EXPORTER}" \
      --config "${DATASETS_JSON_FILE}" --dataset "${row_dataset}" \
      --view "${row_view}" --analysis-variant "${METADATA_VARIANT}" \
      --input-file "${row_input}" --output "${row_output}" --check >/dev/null 2>&1; then
    status_kind="NOOP_VALIDATED"
  else
    "${PYTHON_BIN}" "${EXPORTER}" \
      --config "${DATASETS_JSON_FILE}" --dataset "${row_dataset}" \
      --view "${row_view}" --analysis-variant "${METADATA_VARIANT}" \
      --input-file "${row_input}" --output "${row_output}" ||
      fail "metadata export failed for ${row_dataset}"
  fi
  ecoda_validate_checksum "${row_output}" ||
    fail "metadata export checksum is invalid: ${row_output}"
  safe="$(_ecoda_safe_component "${row_dataset}__${row_view}")"
  ecoda_atomic_write "${STATUS_DIR}/${safe}.status" \
    "STATE=OK\nSTATUS=${status_kind}\nRUN_ID=${RUN_ID}\nANALYSIS_VARIANT=${METADATA_VARIANT}\nANALYSIS_PASS=${METADATA_PASS}\nANALYSIS_ROOT=${METADATA_ROOT}\nDATASET=${row_dataset}\nVIEW=${row_view}\nTASK_ID=${TASK_ID}\nINPUT_FILE=${row_input}\nOUTPUT_FILE=${row_output}\n"
  printf 'H5AD_METADATA_EXPORT=%s\n' "${row_output}"
  exit 0
fi


# A scheduler array supplies one manifest row per view.  A direct invocation
# without a manifest remains useful for diagnostics and runs both views, but
# never does so when a row-isolated task binding is present.
INPUT_FILE="${H5AD_OBS_AUDIT_INPUT_FILE:-${HPC_SCRATCH_DIR}/Covid19_PBMC/data/Covid19_Ren2021.h5ad}"
EXPECTED_INPUT="${HPC_SCRATCH_DIR}/Covid19_PBMC/data/Covid19_Ren2021.h5ad"
[[ "${INPUT_FILE}" == "${EXPECTED_INPUT}" ]] ||
  fail "Covid obs audit input path is not the configured direct source"
OUTPUT_ROOT="${H5AD_OBS_AUDIT_OUTPUT_ROOT:-${RUN_ROOT}/preflight}"
[[ "${OUTPUT_ROOT}" = /* ]] || fail "obs audit output root must be absolute"
[[ "${OUTPUT_ROOT}" == "${RUN_ROOT}/preflight" ]] ||
  fail "obs audit output root must be the run-owned preflight directory"
_ecoda_validate_run_owned_directory_for_create "${OUTPUT_ROOT}" "${RUN_ROOT}" ||
  fail "obs audit output root is not a safe run-owned path"
for audit_view in batch_effect_uncorrected batch_effect_corrected; do
  audit_report="${OUTPUT_ROOT}/Covid19_PBMC_${audit_view}.json"
  [[ ! -L "${audit_report}" && ! -L "${audit_report}.md5" ]] ||
    fail "obs audit report or checksum is a symlink: ${audit_report}"
  _ecoda_validate_run_owned_directory_for_create \
    "$(dirname "${audit_report}")" "${RUN_ROOT}" ||
    fail "obs audit report parent is not run-owned: ${audit_report}"
  _ecoda_validate_run_owned_directory_for_create \
    "$(dirname "${audit_report}.md5")" "${RUN_ROOT}" ||
    fail "obs audit checksum parent is not run-owned: ${audit_report}.md5"
done

AUDIT_MANIFEST="${H5AD_OBS_AUDIT_MANIFEST:-${H5AD_PREFLIGHT_MANIFEST:-}}"
STATUS_DIR="${H5AD_OBS_AUDIT_STATUS_DIR:-${H5AD_PREFLIGHT_STATUS_DIR:-}}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-${H5AD_OBS_AUDIT_TASK_ID:-${H5AD_PREFLIGHT_TASK_ID:-}}}"
VIEWS=()
STATUS_FILE=""
if [[ -n "${AUDIT_MANIFEST}" ]]; then
  [[ -n "${STATUS_DIR}" && "${TASK_ID}" =~ ^[0-9]+$ && ${TASK_ID} -gt 0 ]] ||
    fail "row-isolated obs audit requires manifest, status directory, and task ID"
  [[ -r "${AUDIT_MANIFEST}" && ! -L "${AUDIT_MANIFEST}" ]] ||
    fail "obs audit manifest is missing or unsafe"
  ecoda_validate_run_owned_path "${AUDIT_MANIFEST}" "${RUN_ROOT}" ||
    fail "obs audit manifest is not run-owned"
  ecoda_validate_manifest "${AUDIT_MANIFEST}" 3 ||
    fail "obs audit manifest is malformed"
  [[ "${STATUS_DIR}" = /* && ! -L "${STATUS_DIR}" ]] ||
    fail "obs audit status directory is missing or unsafe"
  _ecoda_validate_run_owned_directory_for_create "${STATUS_DIR}" "${RUN_ROOT}" ||
    fail "obs audit status directory is not a safe run-owned path"
  row=""
  manifest_rows=0
  seen_audit_views=""
  manifest_line=""
  while IFS= read -r manifest_line || [[ -n "${manifest_line}" ]]; do
    IFS=$'\t' read -r row_dataset row_view row_input row_extra <<< "${manifest_line}"
    manifest_rows=$((manifest_rows + 1))
    [[ "${row_dataset}" == "Covid19_PBMC" &&
       ( "${row_view}" == "batch_effect_uncorrected" ||
         "${row_view}" == "batch_effect_corrected" ) &&
       "${row_input}" == "${INPUT_FILE}" && -z "${row_extra}" ]] ||
      fail "obs audit manifest row ${manifest_rows} is not the approved direct Covid binding"
    case " ${seen_audit_views} " in
      *" ${row_view} "*)
        fail "obs audit manifest contains duplicate view: ${row_view}" ;;
      *) seen_audit_views="${seen_audit_views} ${row_view}" ;;
    esac
    if [[ ${manifest_rows} -eq ${TASK_ID} ]]; then
      row="${row_dataset}"$'\t'"${row_view}"$'\t'"${row_input}"
    fi
  done < "${AUDIT_MANIFEST}"
  [[ ${manifest_rows} -gt 0 && ${TASK_ID} -le ${manifest_rows} ]] ||
    fail "obs audit task ID is outside manifest rows"
  [[ -n "${row}" ]] ||
    fail "obs audit selected row is missing"
  IFS=$'\t' read -r row_dataset row_view row_input row_extra <<< "${row}"
  VIEWS=("${row_view}")
  safe="$(_ecoda_safe_component "${row_dataset}__${row_view}")"
  STATUS_FILE="${STATUS_DIR}/${safe}.status"
else
  VIEWS=(batch_effect_uncorrected batch_effect_corrected)
fi

mkdir -p "${OUTPUT_ROOT}"
[[ ! -L "${OUTPUT_ROOT}" ]] || fail "obs audit output root must not be a symlink"
ecoda_validate_run_owned_path "${OUTPUT_ROOT}" "${RUN_ROOT}" ||
  fail "obs audit output root is not run-owned"
if [[ -n "${STATUS_FILE}" ]]; then
  mkdir -p "${STATUS_DIR}"
  ecoda_validate_run_owned_path "${STATUS_DIR}" "${RUN_ROOT}" ||
    fail "obs audit status directory is not run-owned"
fi
[[ -x "${PYTHON_BIN}" ]] || fail "immutable runtime Python is unavailable: ${PYTHON_BIN}"
for view in "${VIEWS[@]}"; do
  report="${OUTPUT_ROOT}/Covid19_PBMC_${view}.json"
  "${PYTHON_BIN}" "${SOURCE_ROOT}/src/3_scrnaseq_preprocessing/1.0_audit_input_views.py" \
    --config "${SOURCE_ROOT}/datasets.json" \
    --input-file "${INPUT_FILE}" \
    --output-root "${OUTPUT_ROOT}" \
    --output "${report}" \
    --view "${view}" --ds-name Covid19_PBMC --obs-only
  [[ -s "${report}" && ! -L "${report}" && -s "${report}.md5" &&
     ! -L "${report}.md5" ]] || fail "obs audit report or checksum is missing: ${report}"
  ecoda_validate_run_owned_path "${report}" "${RUN_ROOT}" ||
    fail "obs audit report is not run-owned: ${report}"
  ecoda_validate_run_owned_path "${report}.md5" "${RUN_ROOT}" ||
    fail "obs audit checksum is not run-owned: ${report}.md5"
  ecoda_validate_checksum "${report}" ||
    fail "obs audit report checksum is invalid: ${report}"
  if [[ -n "${STATUS_FILE}" ]]; then
    ecoda_atomic_write "${STATUS_FILE}" \
      "STATE=OK\nRUN_ID=${RUN_ID}\nDATASET=Covid19_PBMC\nVIEW=${view}\nTASK_ID=${TASK_ID}\nINPUT_FILE=${INPUT_FILE}\nREPORT=${report}\n"
  fi
done
for view in "${VIEWS[@]}"; do
  printf 'H5AD_OBS_AUDIT_REPORT=%s\n' "${OUTPUT_ROOT}/Covid19_PBMC_${view}.json"
done
