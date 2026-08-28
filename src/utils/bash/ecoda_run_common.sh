#!/bin/bash
# Shared run/selection/ownership primitives for ECODA Pipelines 2-5.
# Source after slurm_config.sh. Bash 3.2-compatible: no namerefs, mapfile, or
# associative arrays. Functions mutate documented global variables.

ECODA_RUNS_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs"
ECODA_OWNERS_ROOT="${HPC_SCRATCH_DIR}/_ecoda_owners"
ECODA_RUN_ROOT=""
ECODA_RUN_ID=""
ECODA_ARRAY=()

_ecoda_die() {
  echo "ERROR: $*" >&2
  return 1
}

ecoda_validate_run_id() {
  local run_id="${1:-}"
  [[ "${run_id}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || {
    _ecoda_die "invalid run ID; expected one safe path component: ${run_id}"
    return 1
  }
}

ecoda_realpath_existing() {
  local path="$1"
  [[ -e "${path}" || -L "${path}" ]] || return 1
  command -v realpath >/dev/null 2>&1 || {
    _ecoda_die "realpath is required for run-owned path validation"
    return 1
  }
  realpath "${path}" 2>/dev/null
}

ecoda_validate_run_owned_path() {
  local candidate="${1:-}"
  local run_root="${2:-}"
  local candidate_real root_real
  [[ -n "${candidate}" && -n "${run_root}" ]] || {
    _ecoda_die "run-owned path validation requires candidate and run root"
    return 1
  }
  root_real="$(ecoda_realpath_existing "${run_root}")" || {
    _ecoda_die "run root is missing or cannot be canonicalized: ${run_root}"
    return 1
  }
  candidate_real="$(ecoda_realpath_existing "${candidate}")" || {
    _ecoda_die "run-owned path is missing or cannot be canonicalized: ${candidate}"
    return 1
  }
  case "${candidate_real}" in
    "${root_real}"|"${root_real}"/*) return 0 ;;
    *) _ecoda_die "run-owned path escapes ${root_real}: ${candidate}"; return 1 ;;
  esac
}

ECODA_ACQUIRED_OWNERS=()

ecoda_owner_clear_tracked() {
  ECODA_ACQUIRED_OWNERS=()
}

ecoda_owner_track() {
  local owner_dir="$1" owner
  [[ -n "${owner_dir}" ]] || {
    _ecoda_die "cannot track an empty owner"
    return 1
  }
  if [[ ${#ECODA_ACQUIRED_OWNERS[@]} -gt 0 ]]; then
    for owner in "${ECODA_ACQUIRED_OWNERS[@]}"; do
      [[ "${owner}" == "${owner_dir}" ]] && return 0
    done
  fi
  ECODA_ACQUIRED_OWNERS+=("${owner_dir}")
}

ecoda_owner_finalize_tracked() {
  local state="$1" reason="${2:-}" owner rc=0
  if [[ ${#ECODA_ACQUIRED_OWNERS[@]} -gt 0 ]]; then
    for owner in "${ECODA_ACQUIRED_OWNERS[@]}"; do
      if ! ecoda_owner_set_state "${owner}" "${state}" "${reason}"; then
        rc=1
      fi
    done
  fi
  return "${rc}"
}

_ecoda_safe_component() {
  printf '%s' "$1" | tr '/:,\t |' '______'
}

ecoda_atomic_write() {
  local destination="$1"
  local content="$2"
  local parent tmp
  parent="$(dirname "${destination}")"
  mkdir -p "${parent}"
  tmp="${destination}.tmp.$$"
  umask 077
  printf '%b' "${content}" > "${tmp}"
  mv -f "${tmp}" "${destination}"
}

ecoda_atomic_install_manifest() {
  local source="$1"
  local destination="$2"
  local columns="$3"
  local parent tmp
  [[ -r "${source}" ]] || { _ecoda_die "manifest source is unreadable: ${source}"; return 1; }
  parent="$(dirname "${destination}")"
  mkdir -p "${parent}"
  tmp="${destination}.tmp.$$"
  cp "${source}" "${tmp}"
  if ! ecoda_validate_manifest "${tmp}" "${columns}"; then
    rm -f "${tmp}"
    return 1
  fi
  mv -f "${tmp}" "${destination}"
}

ecoda_owner_field() {
  local owner_dir="$1"
  local field="$2"
  [[ -r "${owner_dir}/owner" ]] || return 1
  sed -n "s/^${field}=//p" "${owner_dir}/owner" | head -1
}

ecoda_owner_state() {
  ecoda_owner_field "$1" STATE
}

ecoda_owner_run() {
  ecoda_owner_field "$1" RUN_ID
}

ecoda_owner_reclaim_terminal() {
  local owner_dir="$1"
  local force="${2:-0}"
  local artifact_valid="${3:-1}"
  local state
  state="$(ecoda_owner_state "${owner_dir}" 2>/dev/null || true)"
  case "${state}" in
    OK|FAIL)
      ;;
    ACTIVE)
      _ecoda_die "active owner cannot be reclaimed: ${owner_dir}"
      return 1
      ;;
    *)
      _ecoda_die "owner state is missing or invalid: ${owner_dir}"
      return 1
      ;;
  esac
  if [[ "${force}" != "1" && "${artifact_valid}" != "0" ]]; then
    _ecoda_die "valid terminal owner must be skipped or forced: ${owner_dir}"
    return 1
  fi
  rm -f "${owner_dir}/owner"
  rmdir "${owner_dir}" 2>/dev/null || {
    _ecoda_die "cannot reclaim terminal owner: ${owner_dir}"
    return 1
  }
}

ecoda_new_run_id() {
  local stage="$1"
  printf '%s_%s_%s' "${stage}" "$(date +%Y%m%d%H%M%S)" "$$"
}

ecoda_init_run() {
  local stage="$1"
  local requested_id="${2:-}"
  local run_id="${requested_id:-$(ecoda_new_run_id "${stage}")}"
  ecoda_validate_run_id "${run_id}" || return 1
  local root="${ECODA_RUNS_ROOT}/${run_id}"
  if [[ -e "${root}" ]]; then
    _ecoda_die "run root already exists: ${root}"
    return 1
  fi
  mkdir -p "${root}/manifests" "${root}/status" "${root}/logs"
  ECODA_RUN_ID="${run_id}"
  ECODA_RUN_ROOT="${root}"
  ecoda_owner_clear_tracked
  ecoda_atomic_write "${root}/metadata" \
    "STAGE=${stage}\nRUN_ID=${run_id}\nSTATE=ACTIVE\nPID=$$\nCREATED=$(date -u +%Y-%m-%dT%H:%M:%SZ)\n"
  printf '%s' "${run_id}"
}

ecoda_open_run() {
  local run_id="$1"
  ecoda_validate_run_id "${run_id}" || return 1
  local root="${ECODA_RUNS_ROOT}/${run_id}"
  [[ -d "${root}" ]] || { _ecoda_die "run root does not exist: ${root}"; return 1; }
  [[ -r "${root}/metadata" ]] || { _ecoda_die "run metadata missing: ${root}/metadata"; return 1; }
  ECODA_RUN_ID="${run_id}"
  ECODA_RUN_ROOT="${root}"
  ecoda_owner_clear_tracked
}

ecoda_set_run_state() {
  local state="$1"
  local reason="${2:-}"
  [[ -n "${ECODA_RUN_ROOT}" ]] || { _ecoda_die "run root is not open"; return 1; }
  ecoda_atomic_write "${ECODA_RUN_ROOT}/status/terminal" \
    "STATE=${state}\nRUN_ID=${ECODA_RUN_ID}\nREASON=${reason}\nTIME=$(date -u +%Y-%m-%dT%H:%M:%SZ)\n"
}

ecoda_split_csv() {
  local csv="${1:-}"
  ECODA_ARRAY=()
  [[ -n "${csv}" ]] || { _ecoda_die "selection list must not be empty"; return 1; }
  local old_ifs="${IFS}" item
  IFS=','
  read -r -a ECODA_ARRAY <<< "${csv}"
  IFS="${old_ifs}"
  [[ ${#ECODA_ARRAY[@]} -gt 0 ]] || { _ecoda_die "selection list must not be empty"; return 1; }
  for item in "${ECODA_ARRAY[@]}"; do
    [[ -n "${item}" ]] || { _ecoda_die "selection list contains an empty item: ${csv}"; return 1; }
    [[ "${item}" != *$'\t'* && "${item}" != *$'\n'* ]] || {
      _ecoda_die "selection item contains a tab/newline: ${item}"; return 1;
    }
  done
}

ecoda_assert_unique_items() {
  local seen="" item
  for item in "$@"; do
    case " ${seen} " in
      *" ${item} "*) _ecoda_die "duplicate selection item: ${item}"; return 1 ;;
    esac
    seen="${seen} ${item}"
  done
}

ecoda_dataset_exists() {
  local ds="$1"
  jq -e --arg ds "${ds}" 'has($ds)' "${DATASETS_JSON_FILE}" >/dev/null 2>&1
}

ecoda_view_field() {
  local ds="$1" view="$2" field="$3"
  jq -r --arg ds "${ds}" --arg view "${view}" --arg field "${field}" \
    '.[$ds].views[$view][$field] // empty' "${DATASETS_JSON_FILE}"
}

ecoda_view_exists() {
  local ds="$1" view="$2"
  jq -e --arg ds "${ds}" --arg view "${view}" \
    '.[$ds].views[$view] != null' "${DATASETS_JSON_FILE}" >/dev/null 2>&1
}

ecoda_view_input_name() {
  local value
  value="$(ecoda_view_field "$1" "$2" input_file_name)"
  if [[ -z "${value}" ]]; then
    value="$(ecoda_view_field "$1" "$2" input_file)"
  fi
  printf '%s' "${value}"
}

ecoda_view_output_name() {
  local value
  value="$(ecoda_view_field "$1" "$2" output_file_name)"
  if [[ -z "${value}" ]]; then
    value="$(ecoda_view_field "$1" "$2" output_file)"
  fi
  printf '%s' "${value}"
}

ecoda_md5_file() {
  local path="$1"
  [[ -s "${path}" ]] || return 1
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "${path}" | cut -d' ' -f1
  elif command -v md5 >/dev/null 2>&1; then
    md5 -q "${path}"
  else
    _ecoda_die "neither md5sum nor md5 is available"
    return 1
  fi
}

ecoda_write_checksum() {
  local path="$1" sidecar="${2:-${1}.md5}" digest size
  [[ -s "${path}" ]] || { _ecoda_die "cannot checksum missing/empty artifact: ${path}"; return 1; }
  digest="$(ecoda_md5_file "${path}")" || return 1
  size="$(wc -c < "${path}" | tr -d '[:space:]')"
  ecoda_atomic_write "${sidecar}" "MD5=${digest}\nSIZE=${size}\nPATH=${path}\n"
}

ecoda_invalidate_artifact() {
  local path="$1"
  shift
  rm -f "${path}.md5"
  local marker
  for marker in "$@"; do
    rm -f "${marker}"
  done
}

ecoda_validate_checksum() {
  local path="$1" sidecar="${2:-${1}.md5}" expected actual expected_size actual_size recorded_path
  [[ -s "${path}" ]] || return 1
  [[ -s "${sidecar}" ]] || return 1
  expected="$(sed -n 's/^MD5=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  expected_size="$(sed -n 's/^SIZE=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  recorded_path="$(sed -n 's/^PATH=//p' "${sidecar}" | head -1)"
  [[ "${recorded_path}" == "${path}" ]] || return 1
  [[ "${expected}" =~ ^[[:xdigit:]]{32}$ ]] || return 1
  actual="$(ecoda_md5_file "${path}")" || return 1
  [[ "${actual}" == "${expected}" ]] || return 1
  if [[ -n "${expected_size}" ]]; then
    actual_size="$(wc -c < "${path}" | tr -d '[:space:]')"
    [[ "${actual_size}" == "${expected_size}" ]] || return 1
  fi
}

ecoda_validate_checksum_remote() {
  local path="$1" sidecar="${2:-${1}.md5}" expected actual expected_size actual_size
  [[ -s "${path}" && -s "${sidecar}" ]] || return 1
  expected="$(sed -n 's/^MD5=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  expected_size="$(sed -n 's/^SIZE=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  [[ "${expected}" =~ ^[[:xdigit:]]{32}$ ]] || return 1
  actual="$(ecoda_md5_file "${path}")" || return 1
  [[ "${actual}" == "${expected}" ]] || return 1
  if [[ -n "${expected_size}" ]]; then
    actual_size="$(wc -c < "${path}" | tr -d '[:space:]')"
    [[ "${actual_size}" == "${expected_size}" ]] || return 1
  fi
}

ecoda_compare_checksum_remote() {
  local local_path="$1"
  local remote_path="$2"
  local remote_sidecar="${3:-${remote_path}.md5}"
  local local_digest remote_digest remote_digest_actual local_size remote_size expected_size recorded_path
  ecoda_validate_checksum "${local_path}" || return 1
  [[ -s "${remote_path}" && -s "${remote_sidecar}" ]] || return 1
  remote_digest="$(sed -n 's/^MD5=//p' "${remote_sidecar}" | head -1 | tr -d '[:space:]')"
  expected_size="$(sed -n 's/^SIZE=//p' "${remote_sidecar}" | head -1 | tr -d '[:space:]')"
  recorded_path="$(sed -n 's/^PATH=//p' "${remote_sidecar}" | head -1)"
  [[ "${remote_digest}" =~ ^[[:xdigit:]]{32}$ ]] || return 1
  [[ "${expected_size}" =~ ^[0-9]+$ ]] || return 1
  local_digest="$(ecoda_md5_file "${local_path}")" || return 1
  remote_digest_actual="$(ecoda_md5_file "${remote_path}")" || return 1
  local_size="$(wc -c < "${local_path}" | tr -d '[:space:]')"
  remote_size="$(wc -c < "${remote_path}" | tr -d '[:space:]')"
  [[ "${local_digest}" == "${remote_digest_actual}" ]] || return 1
  [[ "${remote_digest_actual}" == "${remote_digest}" ]] || return 1
  [[ "${local_size}" == "${remote_size}" && "${remote_size}" == "${expected_size}" ]] || return 1
  [[ -n "${recorded_path}" ]] || return 1
}

ecoda_owner_dir() {
  local stage="$1" key="$2"
  printf '%s/%s/%s' "${ECODA_OWNERS_ROOT}" "${stage}" "$(_ecoda_safe_component "${key}")"
}

ecoda_owner_acquire() {
  local stage="$1" key="$2" run_id="$3" force="${4:-0}"
  local artifact_valid="${5:-1}"
  local owner_dir state owner_run
  owner_dir="$(ecoda_owner_dir "${stage}" "${key}")"
  mkdir -p "$(dirname "${owner_dir}")"
  if mkdir "${owner_dir}" 2>/dev/null; then
    ecoda_atomic_write "${owner_dir}/owner" \
      "RUN_ID=${run_id}\nSTATE=ACTIVE\nSTAGE=${stage}\nKEY=${key}\nPID=$$\n"
    printf '%s' "${owner_dir}"
    return 0
  fi

  owner_run="$(ecoda_owner_run "${owner_dir}" 2>/dev/null || true)"
  state="$(ecoda_owner_state "${owner_dir}" 2>/dev/null || true)"
  if [[ "${state}" == "ACTIVE" ]]; then
    _ecoda_die "active owner ${owner_run:-unknown} already owns ${stage}/${key}"
    return 1
  fi
  case "${state}" in
    OK|FAIL)
      ;;
    *)
      _ecoda_die "owner state is missing or invalid for ${stage}/${key}"
      return 1
      ;;
  esac
  if [[ "${force}" != "1" && "${artifact_valid}" != "0" ]]; then
    return 2
  fi
  ecoda_owner_reclaim_terminal "${owner_dir}" "${force}" "${artifact_valid}" || return 1
  mkdir "${owner_dir}" 2>/dev/null || {
    _ecoda_die "owner was claimed concurrently for ${stage}/${key}"
    return 1
  }
  ecoda_atomic_write "${owner_dir}/owner" \
    "RUN_ID=${run_id}\nSTATE=ACTIVE\nSTAGE=${stage}\nKEY=${key}\nPID=$$\n"
  printf '%s' "${owner_dir}"
}

ecoda_owner_set_state() {
  local owner_dir="$1" state="$2" reason="${3:-}"
  [[ -d "${owner_dir}" && -r "${owner_dir}/owner" ]] || {
    _ecoda_die "owner directory or state file missing: ${owner_dir}"
    return 1
  }
  case "${state}" in
    ACTIVE|OK|FAIL) ;;
    *) _ecoda_die "invalid owner state: ${state}"; return 1 ;;
  esac
  local run_id stage key
  run_id="$(ecoda_owner_field "${owner_dir}" RUN_ID)"
  stage="$(ecoda_owner_field "${owner_dir}" STAGE)"
  key="$(ecoda_owner_field "${owner_dir}" KEY)"
  [[ -n "${run_id}" && -n "${stage}" && -n "${key}" ]] || {
    _ecoda_die "owner state metadata is incomplete: ${owner_dir}"
    return 1
  }
  ecoda_atomic_write "${owner_dir}/owner" \
    "RUN_ID=${run_id}\nSTATE=${state}\nSTAGE=${stage}\nKEY=${key}\nREASON=${reason}\nTIME=$(date -u +%Y-%m-%dT%H:%M:%SZ)\n"
}

ecoda_validate_manifest() {
  local manifest="$1" columns="$2" line expected_count=0
  [[ -s "${manifest}" ]] || { _ecoda_die "manifest is missing or empty: ${manifest}"; return 1; }
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -n "${line}" ]] || { _ecoda_die "manifest contains a blank row: ${manifest}"; return 1; }
    case "${line}" in
      \#*) _ecoda_die "manifest must be headerless: ${manifest}"; return 1 ;;
    esac
    expected_count=$((expected_count + 1))
  done < "${manifest}"
  [[ ${expected_count} -gt 0 ]] || { _ecoda_die "manifest has no rows: ${manifest}"; return 1; }
  if [[ "${columns}" =~ ^[2-9][0-9]*$ ]]; then
    awk -F '\t' -v expected="${columns}" 'NF != expected { exit 1 }' "${manifest}" || {
      _ecoda_die "manifest row has unexpected column count: ${manifest}"
      return 1
    }
  else
    _ecoda_die "unsupported manifest column count: ${columns}"
    return 1
  fi
}

# Validate the immutable twelve-cohort uncorrected batch selection. The helper
# is deliberately stricter than ecoda_validate_manifest: it checks row order,
# dataset identity, view identity, and (for Stage 5) the third label field.
ecoda_validate_exact_batch_selection() {
  local manifest="$1"
  local columns="$2"
  local expected_count=12
  local count=0
  local ds view label
  local expected
  local expected_datasets=(
    Alzheimer
    Breast_cancer
    Covid19_PBMC
    Kidney_KPMP
    Myocardial_infarction
    Diabetes
    Lupus_PBMC
    Lung
    Parkinson
    Joanito
    Stephenson
    CombinedPBMC
  )

  [[ "${columns}" == "2" || "${columns}" == "3" ]] || {
    _ecoda_die "exact batch selection requires two or three columns"
    return 1
  }
  ecoda_validate_manifest "${manifest}" "${columns}" || return 1

  while IFS=$'\t' read -r ds view label; do
    count=$((count + 1))
    expected="${expected_datasets[$((count - 1))]}"
    [[ "${ds}" == "${expected}" ]] || {
      _ecoda_die "exact batch selection row ${count} must be ${expected}, got ${ds}"
      return 1
    }
    [[ "${view}" == "batch_effect_uncorrected" ]] || {
      _ecoda_die "exact batch selection row ${count} must use batch_effect_uncorrected"
      return 1
    }
    if [[ "${columns}" == "3" && "${label}" != "batch_effect_uncorrected" ]]; then
      _ecoda_die "exact batch selection row ${count} label must be batch_effect_uncorrected"
      return 1
    fi
  done < "${manifest}"

  [[ ${count} -eq ${expected_count} ]] || {
    _ecoda_die "exact batch selection requires exactly ${expected_count} rows"
    return 1
  }
}

_ecoda_accounting_active() {
  case "$1" in
    PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING) return 0 ;;
    *) return 1 ;;
  esac
}

ecoda_wait_scalar_accounting() {
  local job="$1" poll_seconds="${2:-30}" rows state empty=0 unresolved=0
  while :; do
    rows="$(sacct -j "${job}" -X -n -P --format=JobIDRaw,State,ExitCode 2>/dev/null || true)"
    if [[ -z "${rows//[[:space:]]/}" ]]; then
      empty=$((empty + 1))
      [[ ${empty} -lt ${ECODA_ACCOUNTING_EMPTY_GRACE:-3} ]] || return 1
    else
      empty=0
      state="$(printf '%s\n' "${rows}" | awk -F '|' 'NR == 1 {print $2}')"
      [[ -n "${state}" ]] || state="$(printf '%s\n' "${rows}" | awk 'NR == 1 {print $1}')"
      state="${state%%+*}"
      if [[ -n "${state}" ]] && ! _ecoda_accounting_active "${state}"; then
        ECODA_ACCOUNTING_STATE="${state}"
        ECODA_ACCOUNTING_ROWS="${rows}"

        return 0
      fi
      if [[ -z "${state}" ]]; then
        unresolved=$((unresolved + 1))
        [[ ${unresolved} -lt ${ECODA_ACCOUNTING_EMPTY_GRACE:-3} ]] || return 1
      else
        unresolved=0
      fi
    fi
    sleep "${poll_seconds}"
  done
}

ecoda_wait_array_accounting() {
  local job="$1" expected="$2" poll_seconds="${3:-30}"
  local rows jid state found pending empty=0 missing=0 scheduler_active active_jobs
  while :; do
    rows="$(sacct -j "${job}" -n -P --format=JobIDRaw,State,ExitCode 2>/dev/null || true)"
    scheduler_active=0
    if command -v squeue >/dev/null 2>&1; then
      active_jobs="$(squeue -j "${job}" -h -o "%A" 2>/dev/null || true)"
      case " ${active_jobs} " in
        *" ${job} "*) scheduler_active=1 ;;
      esac
    fi
    if [[ -z "${rows//[[:space:]]/}" ]]; then
      if [[ ${scheduler_active} -eq 1 ]]; then
        empty=0
        missing=0
        sleep "${poll_seconds}"
        continue
      fi
      empty=$((empty + 1))
      [[ ${empty} -lt ${ECODA_ACCOUNTING_EMPTY_GRACE:-3} ]] || return 1
    else
      empty=0
    fi
    found=0
    pending=0
    while IFS='|' read -r jid state exitcode; do
      [[ "${jid}" =~ ^${job}_[0-9]+$ ]] || continue
      state="${state%%+*}"
      found=$((found + 1))
      _ecoda_accounting_active "${state}" && pending=1
    done <<< "${rows}"
    if [[ ${found} -lt ${expected} ]]; then
      if [[ ${scheduler_active} -eq 1 ]]; then
        missing=0
      else
        missing=$((missing + 1))
        [[ ${missing} -lt ${ECODA_ACCOUNTING_EMPTY_GRACE:-3} ]] || return 1
      fi
    else
      missing=0
    fi
    if [[ ${found} -ge ${expected} && ${pending} -eq 0 ]]; then
      ECODA_ACCOUNTING_ROWS="${rows}"
      return 0
    fi
    sleep "${poll_seconds}"
  done
}

# Return 0 when every path has a valid nonempty checksum; this is deliberately
# independent from schema validators so callers can fail closed on either.
ecoda_validate_artifacts() {
  local path
  for path in "$@"; do
    ecoda_validate_checksum "${path}" || return 1
  done
}
# Validate Stage 2 derived outputs semantically after checksum validation.
# These checks are intentionally delegated to the pinned Python/R runtimes so
# a stale or structurally plausible artifact cannot satisfy a prerequisite.
ecoda_validate_stage2_output() {
  local step="$1"
  local path="$2"
  case "${step}" in
    myocardial_counts)
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/derived_prerequisite_contract.py" \
        --path "${path}" --kind myocardial >/dev/null 2>&1
      ;;
    combinedpbmc)
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/derived_prerequisite_contract.py" \
        --path "${path}" --kind combinedpbmc >/dev/null 2>&1
      ;;
    joanito)
      case "${path}" in
        *.rds)
          joanito_check="${TMPDIR:-/tmp}/ecoda_joanito_check.$$.R"
          if ! cat > "${joanito_check}" <<'RSCRIPT'
p <- commandArgs(trailingOnly = TRUE)[1]
x <- readRDS(p)
md <- if ("meta.data" %in% slotNames(x)) x@meta.data else x
needed <- c("dataset", "cell.type", "iCMS", "seqtec", "cell.type_new")
if (!is.data.frame(md) || !all(needed %in% colnames(md))) {
  stop("Joanito RDS lacks raw/derived metadata columns")
}
expected_seqtec <- ifelse(
  as.character(md$dataset) %in% c("CRC-SG1", "KUL5"),
  "5' seq", "3' seq"
)
base <- as.character(md$cell.type)
expected_new <- base
has_icms <- !is.na(md$iCMS)
expected_new[has_icms] <- paste0(
  base[has_icms], "_",
  ifelse(as.character(md$iCMS[has_icms]) == "Normal",
         "Normal", "Cancer")
)
same <- function(actual, expected) {
  actual <- as.character(actual)
  (is.na(actual) & is.na(expected)) |
    (!is.na(actual) & !is.na(expected) & actual == expected)
}
if (!all(same(md$seqtec, expected_seqtec)) ||
    !all(same(md$cell.type_new, expected_new))) {
  stop("Joanito derived metadata is stale")
}
for (nm in c("seqtec", "cell.type_new")) {
  v <- as.character(md[[nm]])
  if (!any(!is.na(v) & nzchar(trimws(v)))) stop("Joanito derived column is empty")
}
RSCRIPT
          then
            return 1
          fi
          if ${PIXI_RSCRIPT} "${joanito_check}" "${path}" >/dev/null 2>&1; then
            joanito_rc=0
          else
            joanito_rc=$?
          fi
          rm -f "${joanito_check}"
          return "${joanito_rc}"
          ;;
        *.h5ad)
          "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/derived_prerequisite_contract.py" \
            --path "${path}" --kind joanito-debug >/dev/null 2>&1
          ;;
        *) return 1 ;;
      esac
      ;;
    *)
      return 0
      ;;
  esac
}
