#!/bin/bash
# Canonical manifest-driven Stage 3 preprocessing gate.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
FORCE_ARG=0
SYNC_ONLY_RUN=""
SYNC_ONLY_SET=0
MEMORY="128G"
MAX_MEMORY="500G"
PARTITION="${SLURM_PARTITION}"
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
RUNTIME_EXPORT=""

usage() {
  cat <<'EOF'
Usage: 1_submit_hpc_array.sh [--datasets LIST] [--views LIST]
       [--selection-file TSV] [--exact-batch-selection] [--force]
       [--sync-only RUN_ID] [--mem VALUE] [--max-mem VALUE]
       [--partition NAME] [--throttle N]

Each manifest row is DATASET<TAB>VIEW. --ds_name and --view remain accepted
as compatibility aliases for one dataset/view selection. Exact batch mode
requires the immutable twelve-row uncorrected selection file.
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
mkdir -p "${LOGS_DIR}" || {
  echo "ERROR: could not create Stage 3 log directory." >&2
  exit 1
}
export ECODA_RUNTIME_PROFILE=stage3
ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE}" || {
  echo "ERROR: Stage 3 immutable runtime validation failed." >&2
  exit 1
}
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage3 0)" || {
  echo "ERROR: Stage 3 runtime export construction failed." >&2
  exit 1
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
  local preflight_id preflight_rc count=0 status_count=0
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
  set +e
  preflight_id="$(
    ecoda_submit_h5ad_preflight "${preflight_manifest}" "${status_dir}" \
      "${ECODA_RUN_ROOT}" classify "${PARTITION}" "${MEMORY}" "${THROTTLE}" \
      "${LOGS_DIR}" stage3 "${SCRIPT_DIR}/../utils/bash/h5ad_preflight_worker.sh" \
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
      OK|REBUILD) ;;
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
    OK) return 0 ;;
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
    if ! mkdir -p "${remote_dir}" ||
       ! ecoda_validate_checksum "${path}"; then
      rc=1
      continue
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
    ecoda_validate_checksum "${path}" || failed=1
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
if [[ -n "${SELECTION_FILE_ARG}" ]]; then
  validate_external_selection "${SELECTION_FILE_ARG}" || {
    echo "ERROR: Stage 3 selection file is malformed or semantically invalid." >&2
    exit 1
  }
fi

RUN_ID="${ECODA_RUN_ID:-$(ecoda_new_run_id stage3)}"
ecoda_init_run stage3 "${RUN_ID}" >/dev/null
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
SCHEDULER_IDS_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
ecoda_atomic_write "${SCHEDULER_IDS_FILE}" "" ||
  stage3_abort "failed to initialize Stage 3 scheduler manifest"
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
  if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" != "1" ]] &&
     ! sync_selected "${MANIFEST}"; then
    stage3_abort "selected Stage 3 sync failed"
  fi
  ecoda_set_run_state OK "all selected Stage 3 artifacts already validated and synced" ||
    stage3_abort "failed to write Stage 3 terminal OK state"
  echo "PREPROCESS_RUN_ID=${RUN_ID}"
  exit 0
fi

mkdir -p "${LOGS_DIR}" || stage3_abort "failed to create Stage 3 log directory"
export FORCE_PREPROCESS="${FORCE_ARG}"
export PREPROCESS_SELECTION_FILE="${PENDING_MANIFEST}"
export PREPROCESS_RUN_ROOT="${ECODA_RUN_ROOT}"
export PREPROCESS_ERROR_PREFIX="${LOGS_DIR}/3_scrnaseq_preprocessing"
set +e
ARRAY_MSG="$(sbatch --parsable --array="1-${PENDING_COUNT}%${THROTTLE}" \
  --mem="${MEMORY}" --partition="${PARTITION}" \
  --output="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.log" \
  --error="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.err" \
  --mail-user="${USER_EMAIL}" \
  --export="ALL,PREPROCESS_SELECTION_FILE=${PENDING_MANIFEST},PREPROCESS_RUN_ROOT=${ECODA_RUN_ROOT},FORCE_PREPROCESS=${FORCE_ARG},PREPROCESS_ERROR_PREFIX=${LOGS_DIR}/3_scrnaseq_preprocessing,${RUNTIME_EXPORT}" \
  "${SCRIPT_DIR}/1.1_run_worker.sh")"
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
set +e
WATCHDOG_MSG="$(sbatch --parsable --wait --dependency="afterany:${ARRAY_ID}" \
  --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=2G \
  --time="${STAGE3_WATCHDOG_TIME_LIMIT:-12:00:00}" \
  --output="${LOGS_DIR}/3_scrnaseq_preprocessing_watchdog_%j.log" \
  --error="${LOGS_DIR}/3_scrnaseq_preprocessing_watchdog_%j.err" \
  --mail-user="${USER_EMAIL}" \
  --export="ALL,PREPROCESS_RUN_ROOT=${ECODA_RUN_ROOT},PREPROCESS_PENDING_MANIFEST=${PENDING_MANIFEST},${RUNTIME_EXPORT}" \
  "${SCRIPT_DIR}/1.2_preprocess_watchdog.sh" "${RUN_ID}" "${MANIFEST}" \
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
