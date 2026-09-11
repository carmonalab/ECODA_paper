#!/bin/bash
# Shared run/selection/ownership primitives for ECODA Pipelines 2-5.
# Source after slurm_config.sh. Bash 3.2-compatible: no namerefs, mapfile, or
# associative arrays. Functions mutate documented global variables.

ECODA_RUNS_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs"
ECODA_OWNERS_ROOT="${HPC_SCRATCH_DIR}/_ecoda_owners"
ECODA_RUN_ROOT="${ECODA_RUN_ROOT-}"
ECODA_RUN_ID="${ECODA_RUN_ID-}"
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
_ecoda_validate_path_ancestors() {
  local candidate="${1:-}" boundary="${2:-}" current parent
  local boundary_root=""
  [[ "${candidate}" = /* && "${candidate}" != *$'\n'* &&
     "${candidate}" != *$'\t'* ]] || {
    _ecoda_die "path must be absolute and free of record delimiters: ${candidate}"
    return 1
  }
  if [[ -n "${boundary}" ]]; then
    [[ "${boundary}" = /* && "${boundary}" != *$'\n'* &&
       "${boundary}" != *$'\t'* ]] || return 1
    boundary_root="${boundary%/}"
    [[ -n "${boundary_root}" ]] || boundary_root="/"
    case "${candidate}" in
      "${boundary_root}"|"${boundary_root}"/*)
        [[ "${candidate}" == "${boundary_root}" ]] && return 0
        ;;
      *) boundary_root="" ;;
    esac
  fi
  # The final component may be a symlink for read-only callers that need its
  # canonical target.  When a trusted configured root is supplied, reject
  # redirects below that root while allowing system prefixes such as
  # macOS's /var -> /private/var.
  current="$(dirname "${candidate}")"
  while :; do
    if [[ -n "${boundary_root}" && "${current}" == "${boundary_root}" ]]; then
      break
    fi
    if [[ -n "${boundary_root}" ]]; then
      [[ ! -L "${current}" ]] || {
        _ecoda_die "path has a symlinked parent component: ${current}"
        return 1
      }
    fi
    [[ ! -e "${current}" || -d "${current}" ]] || {
      _ecoda_die "path parent component is not a directory: ${current}"
      return 1
    }
    [[ "${current}" == "/" ]] && break
    parent="$(dirname "${current}")"
    [[ "${parent}" != "${current}" ]] || break
    current="${parent}"
  done
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
  declare -p ECODA_ACQUIRED_OWNERS >/dev/null 2>&1 ||
    ECODA_ACQUIRED_OWNERS=()
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
  _ecoda_init_output_arrays
  declare -p ECODA_ACQUIRED_OWNERS >/dev/null 2>&1 ||
    ECODA_ACQUIRED_OWNERS=()
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

ECODA_RECLAIM_TOMBSTONE=""

_ecoda_owner_reclaim_to_tombstone() {
  local owner_dir="${1:-}" tombstone suffix=0
  ECODA_RECLAIM_TOMBSTONE=""
  [[ -d "${owner_dir}" && ! -L "${owner_dir}" ]] || {
    _ecoda_die "owner directory is not a regular directory: ${owner_dir}"
    return 1
  }
  tombstone="${owner_dir}.reclaim.$$"
  while [[ -e "${tombstone}" || -L "${tombstone}" ]]; do
    suffix=$((suffix + 1))
    tombstone="${owner_dir}.reclaim.$$.${suffix}"
  done
  mv "${owner_dir}" "${tombstone}" 2>/dev/null || {
    _ecoda_die "cannot move owner to reclaim tombstone: ${owner_dir}"
    return 1
  }
  [[ ! -e "${owner_dir}" && ! -L "${owner_dir}" &&
     -d "${tombstone}" && ! -L "${tombstone}" ]] || {
    _ecoda_die "owner reclaim rename was unsafe: ${owner_dir}"
    return 1
  }
  ECODA_RECLAIM_TOMBSTONE="${tombstone}"
}

_ecoda_owner_reclaim_tombstone_clean() {
  local tombstone="${1:-}"
  [[ -d "${tombstone}" && ! -L "${tombstone}" ]] || return 1
  rm -f "${tombstone}/owner" || return 1
  rmdir "${tombstone}" 2>/dev/null || {
    _ecoda_die "cannot clean owner reclaim tombstone: ${tombstone}"
    return 1
  }
}

_ecoda_owner_reclaim_tombstone_restore() {
  local owner_dir="${1:-}" tombstone="${2:-}"
  [[ -d "${tombstone}" && ! -L "${tombstone}" ]] || return 1
  [[ ! -e "${owner_dir}" && ! -L "${owner_dir}" ]] || return 1
  mkdir "${owner_dir}" 2>/dev/null || return 1
  mv "${tombstone}/owner" "${owner_dir}/owner" 2>/dev/null &&
    rmdir "${tombstone}" 2>/dev/null
}

ecoda_owner_reclaim_terminal() {
  local owner_dir="$1"
  local force="${2:-0}"
  local artifact_valid="${3:-1}"
  local state tombstone
  [[ -d "${owner_dir}" && ! -L "${owner_dir}" ]] || {
    _ecoda_die "owner directory is not a regular directory: ${owner_dir}"
    return 1
  }
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
  _ecoda_owner_reclaim_to_tombstone "${owner_dir}" || return 1
  tombstone="${ECODA_RECLAIM_TOMBSTONE}"
  _ecoda_owner_reclaim_tombstone_clean "${tombstone}"
}


ecoda_new_run_id() {
  local stage="$1"
  printf '%s_%s_%s' "${stage}" "$(date +%Y%m%d%H%M%S)" "$$"
}
_ecoda_runs_root_canonical() {
  local create="${1:-0}" root="${ECODA_RUNS_ROOT:-}" canonical boundary
  [[ "${root}" = /* && "${root}" != *$'\n'* &&
     "${root}" != *$'\t'* ]] || {
    _ecoda_die "ECODA_RUNS_ROOT must be an absolute path: ${root}"
    return 1
  }
  boundary="${HPC_SCRATCH_DIR:-}"
  if [[ -z "${boundary}" || "${boundary}" != /* ]]; then
    boundary="$(dirname "${root}")"
  else
    boundary="${boundary%/}"
    [[ -n "${boundary}" ]] || boundary="/"
    case "${root}" in
      "${boundary}"|"${boundary}"/*) ;;
      *) boundary="$(dirname "${root}")" ;;
    esac
  fi
  _ecoda_validate_path_ancestors "${root}" "${boundary}" || return 1
  if [[ ! -e "${root}" && ! -L "${root}" ]]; then
    [[ "${create}" == "1" ]] || {
      _ecoda_die "ECODA_RUNS_ROOT does not exist: ${root}"
      return 1
    }
    mkdir -p "${root}" || {
      _ecoda_die "could not create ECODA_RUNS_ROOT: ${root}"
      return 1
    }
  fi
  [[ -d "${root}" && ! -L "${root}" ]] || {
    _ecoda_die "ECODA_RUNS_ROOT is not a regular directory: ${root}"
    return 1
  }
  _ecoda_validate_path_ancestors "${root}" "${boundary}" || return 1
  canonical="$(ecoda_realpath_existing "${root}")" || {
    _ecoda_die "ECODA_RUNS_ROOT cannot be canonicalized: ${root}"
    return 1
  }
  printf '%s' "${canonical}"
}

_ecoda_validate_run_metadata() {
  local metadata="$1" expected_run="$2" expected_stage="${3:-}"
  [[ -f "${metadata}" && ! -L "${metadata}" && -r "${metadata}" &&
     -s "${metadata}" ]] || {
    _ecoda_die "run metadata is missing or not a regular readable file: ${metadata}"
    return 1
  }
  stage_count="$(sed -n 's/^STAGE=//p' "${metadata}" | wc -l | tr -d '[:space:]')"
  run_count="$(sed -n 's/^RUN_ID=//p' "${metadata}" | wc -l | tr -d '[:space:]')"
  [[ "${stage_count}" == "1" && "${run_count}" == "1" ]] || {
    _ecoda_die "run metadata must contain one STAGE and RUN_ID: ${metadata}"
    return 1
  }
  metadata_stage="$(sed -n 's/^STAGE=//p' "${metadata}" | head -1)"
  metadata_run="$(sed -n 's/^RUN_ID=//p' "${metadata}" | head -1)"
  case "${metadata_stage}" in
    stage2|stage3|stage4|stage5) ;;
    *)
      _ecoda_die "run metadata has an invalid stage: ${metadata}"
      return 1
      ;;
  esac
  [[ "${metadata_run}" == "${expected_run}" ]] || {
    _ecoda_die "run metadata RUN_ID does not match requested run: ${metadata}"
    return 1
  }
  if [[ -n "${expected_stage}" ]]; then
    [[ "${expected_stage}" == "${metadata_stage}" ]] || {
      _ecoda_die "run metadata stage does not match requested stage: ${metadata}"
      return 1
    }
  fi
}

ecoda_init_run() {
  local stage="$1"
  local requested_id="${2:-}"
  local run_id="${requested_id:-$(ecoda_new_run_id "${stage}")}"
  local runs_root_real root root_real
  [[ -n "${stage}" && "${stage}" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]] || {
    _ecoda_die "invalid run stage: ${stage}"
    return 1
  }
  ecoda_validate_run_id "${run_id}" || return 1
  runs_root_real="$(_ecoda_runs_root_canonical 1)" || return 1
  root="${ECODA_RUNS_ROOT}/${run_id}"
  if [[ -e "${root}" || -L "${root}" ]]; then
    _ecoda_die "run root already exists: ${root}"
    return 1
  fi
  # Claim the exact run directory atomically.  The existence check above is
  # diagnostic only; mkdir is the synchronization boundary.
  mkdir "${root}" 2>/dev/null || {
    _ecoda_die "run root already exists or could not be initialized: ${root}"
    return 1
  }
  root_real="$(ecoda_realpath_existing "${root}")" || {
    rmdir "${root}" 2>/dev/null || true
    _ecoda_die "new run root cannot be canonicalized: ${root}"
    return 1
  }
  [[ "${root_real}" == "${runs_root_real}/${run_id}" &&
     -d "${root}" && ! -L "${root}" ]] || {
    rmdir "${root}" 2>/dev/null || true
    _ecoda_die "new run root has an invalid canonical layout: ${root}"
    return 1
  }
  mkdir "${root}/manifests" "${root}/status" "${root}/logs" || {
    _ecoda_die "could not create run subdirectories: ${root}"
    return 1
  }
  ecoda_atomic_write "${root}/metadata" \
    "STAGE=${stage}\nRUN_ID=${run_id}\nSTATE=ACTIVE\nPID=$$\nCREATED=$(date -u +%Y-%m-%dT%H:%M:%SZ)\n" || return 1
  ECODA_RUN_ID="${run_id}"
  ECODA_RUN_ROOT="${root}"
  ecoda_owner_clear_tracked
  printf '%s' "${run_id}"
}

ecoda_open_run() {
  local run_id="$1"
  local expected_stage="${2:-${ECODA_GATE_STAGE:-}}"
  local runs_root_real root root_real
  ecoda_validate_run_id "${run_id}" || return 1
  runs_root_real="$(_ecoda_runs_root_canonical 0)" || return 1
  root="${ECODA_RUNS_ROOT}/${run_id}"
  [[ -d "${root}" && ! -L "${root}" ]] || {
    _ecoda_die "run root does not exist or is not a regular directory: ${root}"
    return 1
  }
  _ecoda_validate_path_ancestors "${root}" "${ECODA_RUNS_ROOT}" || return 1
  root_real="$(ecoda_realpath_existing "${root}")" || {
    _ecoda_die "run root cannot be canonicalized: ${root}"
    return 1
  }
  [[ "${root_real}" == "${runs_root_real}/${run_id}" ]] || {
    _ecoda_die "run root has an invalid canonical layout: ${root}"
    return 1
  }
  _ecoda_validate_run_metadata "${root_real}/metadata" "${run_id}" \
    "${expected_stage}" || return 1
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

ECODA_CHECKSUM_PATH=""
ECODA_CHECKSUM_MD5=""
ECODA_CHECKSUM_SIZE=""
ECODA_CHECKSUM_STRICT=0

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
  ECODA_CHECKSUM_PATH=""
  ECODA_CHECKSUM_MD5=""
  ECODA_CHECKSUM_SIZE=""
  ECODA_CHECKSUM_STRICT=0
  [[ -s "${path}" ]] || { _ecoda_die "cannot checksum missing/empty artifact: ${path}"; return 1; }
  digest="$(ecoda_md5_file "${path}")" || return 1
  size="$(wc -c < "${path}" | tr -d '[:space:]')" || return 1
  ecoda_atomic_write "${sidecar}" "MD5=${digest}\nSIZE=${size}\nPATH=${path}\n" || return 1
  ECODA_CHECKSUM_PATH="${path}"
  ECODA_CHECKSUM_MD5="${digest}"
  ECODA_CHECKSUM_STRICT=1
  ECODA_CHECKSUM_SIZE="${size}"
}
_ecoda_artifact_mode_bits() {
  local path="${1:-}" mode=""
  mode="$(stat -c '%a' "${path}" 2>/dev/null || true)"
  if [[ ! "${mode}" =~ ^[0-7]+$ ]]; then
    mode="$(stat -f '%Lp' "${path}" 2>/dev/null || true)"
  fi
  [[ "${mode}" =~ ^[0-7]{3,4}$ ]] || {
    _ecoda_die "could not inspect artifact permissions: ${path}"
    return 1
  }
  printf '%s' "${mode: -3}"
}

_ecoda_artifact_is_nonwritable() {
  local path="${1:-}" mode
  [[ -f "${path}" && ! -L "${path}" ]] || return 1
  mode="$(_ecoda_artifact_mode_bits "${path}")" || return 2
  [[ "${mode}" != *[2367]* ]]
}

_ecoda_artifact_make_nonwritable() {
  local path="${1:-}"
  [[ -f "${path}" && ! -L "${path}" ]] || {
    _ecoda_die "published artifact is not a regular file: ${path}"
    return 1
  }
  chmod a-w "${path}" || {
    _ecoda_die "could not make published artifact read-only: ${path}"
    return 1
  }
  _ecoda_artifact_is_nonwritable "${path}" || {
    _ecoda_die "published artifact remains writable: ${path}"
    return 1
  }
}



ecoda_invalidate_artifact() {
  local path="$1" canonical marker
  shift
  if [[ -e "${path}" || -L "${path}" ]]; then
    canonical="$(_ecoda_canonical_path "${path}")" || return 1
    if [[ -f "${canonical}" && ! -L "${canonical}" ]]; then
      chmod u+w "${canonical}" || {
        _ecoda_die "could not reopen artifact for recomputation: ${canonical}"
        return 1
      }
    fi
  fi
  rm -f "${path}.md5" || return 1
  for marker in "$@"; do
    rm -f "${marker}" || return 1
  done
}

_ecoda_checksum_path_matches() {
  local recorded_path="$1" requested_path="$2" recorded_real requested_real
  [[ "${recorded_path}" = /* && "${requested_path}" = /* &&
     "${recorded_path}" != *$'\n'* && "${requested_path}" != *$'\n'* &&
     "${recorded_path}" != *$'\t'* && "${requested_path}" != *$'\t'* ]] ||
    return 1
  [[ "${recorded_path}" == "${requested_path}" ]] && return 0
  recorded_real="$(ecoda_realpath_existing "${recorded_path}" 2>/dev/null)" ||
    return 1
  requested_real="$(ecoda_realpath_existing "${requested_path}" 2>/dev/null)" ||
    return 1
  [[ "${recorded_real}" == "${requested_real}" ]]
}

ecoda_validate_checksum() {
  local path="$1" sidecar="${2:-${1}.md5}" expected actual expected_size actual_size recorded_path
  ECODA_CHECKSUM_PATH=""
  ECODA_CHECKSUM_MD5=""
  ECODA_CHECKSUM_SIZE=""
  ECODA_CHECKSUM_STRICT=0
  [[ -s "${path}" && -s "${sidecar}" ]] || return 1
  expected="$(sed -n 's/^MD5=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  expected_size="$(sed -n 's/^SIZE=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  recorded_path="$(sed -n 's/^PATH=//p' "${sidecar}" | head -1)"
  _ecoda_checksum_path_matches "${recorded_path}" "${path}" || return 1
  [[ "${expected}" =~ ^[[:xdigit:]]{32}$ ]] || return 1
  actual="$(ecoda_md5_file "${path}")" || return 1
  actual_size="$(wc -c < "${path}" | tr -d '[:space:]')" || return 1
  [[ "${actual}" == "${expected}" ]] || return 1
  [[ "${expected_size}" =~ ^[1-9][0-9]*$ &&
     "${actual_size}" == "${expected_size}" ]] || return 1
  ECODA_CHECKSUM_PATH="${path}"
  ECODA_CHECKSUM_MD5="${actual}"
  ECODA_CHECKSUM_STRICT=1
  ECODA_CHECKSUM_SIZE="${actual_size}"
}

# Validate strict sidecar fields against a digest/size computed immediately
# before this no-write confirmation; never use this as a standalone checksum.
ecoda_validate_checksum_record() {
  local path="$1" expected_digest="${2:-}" expected_size="${3:-}"
  local sidecar="${4:-${1}.md5}" recorded_digest recorded_size recorded_path actual_size
  ECODA_CHECKSUM_PATH=""
  ECODA_CHECKSUM_STRICT=0
  ECODA_CHECKSUM_MD5=""
  ECODA_CHECKSUM_SIZE=""
  [[ -s "${path}" && -s "${sidecar}" ]] || return 1
  [[ "${expected_digest}" =~ ^[[:xdigit:]]{32}$ ]] || return 1
  [[ "${expected_size}" =~ ^[1-9][0-9]*$ ]] || return 1
  recorded_digest="$(sed -n 's/^MD5=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  recorded_size="$(sed -n 's/^SIZE=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  recorded_path="$(sed -n 's/^PATH=//p' "${sidecar}" | head -1)"
  _ecoda_checksum_path_matches "${recorded_path}" "${path}" || return 1
  [[ "${recorded_digest}" =~ ^[[:xdigit:]]{32}$ &&
     "${recorded_digest}" == "${expected_digest}" ]] || return 1
  [[ "${recorded_size}" =~ ^[1-9][0-9]*$ &&
     "${recorded_size}" == "${expected_size}" ]] || return 1
  actual_size="$(wc -c < "${path}" | tr -d '[:space:]')" || return 1
  [[ "${actual_size}" == "${expected_size}" ]] || return 1
  ECODA_CHECKSUM_PATH="${path}"
  ECODA_CHECKSUM_MD5="${expected_digest}"
  ECODA_CHECKSUM_SIZE="${expected_size}"
}

ecoda_validate_checksum_remote() {
  local path="$1" sidecar="${2:-${1}.md5}" expected actual expected_size actual_size
  [[ -s "${path}" && -s "${sidecar}" ]] || return 1
  expected="$(sed -n 's/^MD5=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  expected_size="$(sed -n 's/^SIZE=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
  [[ "${expected}" =~ ^[[:xdigit:]]{32}$ &&
     "${expected_size}" =~ ^[1-9][0-9]*$ ]] || return 1
  actual="$(ecoda_md5_file "${path}")" || return 1
  [[ "${actual}" == "${expected}" ]] || return 1
  actual_size="$(wc -c < "${path}" | tr -d '[:space:]')" || return 1
  [[ "${actual_size}" == "${expected_size}" ]] || return 1
}

ecoda_compare_checksum_remote() {
  local local_path="$1"
  local remote_path="$2"
  local remote_sidecar="${3:-${remote_path}.md5}"
  local known_digest="${4:-}" known_size="${5:-}"
  local local_digest remote_digest remote_digest_actual local_size remote_size expected_size recorded_path
  if [[ -n "${known_digest}" || -n "${known_size}" ]]; then
    [[ -n "${known_digest}" && -n "${known_size}" ]] || return 1
    ecoda_validate_checksum_record "${local_path}" "${known_digest}" "${known_size}" || return 1
    local_digest="${known_digest}"
    local_size="${known_size}"
  else
    ecoda_validate_checksum "${local_path}" || return 1
    local_digest="${ECODA_CHECKSUM_MD5}"
    local_size="${ECODA_CHECKSUM_SIZE}"
  fi
  [[ -s "${remote_path}" && -s "${remote_sidecar}" ]] || return 1
  remote_digest="$(sed -n 's/^MD5=//p' "${remote_sidecar}" | head -1 | tr -d '[:space:]')"
  expected_size="$(sed -n 's/^SIZE=//p' "${remote_sidecar}" | head -1 | tr -d '[:space:]')"
  recorded_path="$(sed -n 's/^PATH=//p' "${remote_sidecar}" | head -1)"
  [[ "${remote_digest}" =~ ^[[:xdigit:]]{32}$ ]] || return 1
  [[ "${expected_size}" =~ ^[0-9]+$ ]] || return 1
  remote_digest_actual="$(ecoda_md5_file "${remote_path}")" || return 1
  remote_size="$(wc -c < "${remote_path}" | tr -d '[:space:]')"
  [[ "${local_digest}" == "${remote_digest_actual}" ]] || return 1
  [[ "${remote_digest_actual}" == "${remote_digest}" ]] || return 1
  [[ "${local_size}" == "${remote_size}" && "${remote_size}" == "${expected_size}" ]] || return 1
  [[ -n "${recorded_path}" ]] || return 1
  if [[ "${recorded_path}" != "${remote_path}" ]]; then
    ecoda_atomic_write "${remote_sidecar}" \
      "MD5=${remote_digest_actual}\nSIZE=${remote_size}\nPATH=${remote_path}\n" ||
      return 1
    ecoda_validate_checksum_record "${remote_path}" "${remote_digest_actual}" \
      "${remote_size}" "${remote_sidecar}" || return 1
  fi
}
# SHA-256 helpers are intentionally local to this shared shell library.  The
# runtime helper has a similarly portable implementation; keeping this one
# here avoids requiring callers to source a second file merely to derive the
# bounded artifact/owner key.
_ecoda_sha256_text() {
  local value="${1:-}" digest
  if command -v sha256sum >/dev/null 2>&1; then
    digest="$(printf '%s' "${value}" | sha256sum | awk '{print $1}')"
  elif command -v shasum >/dev/null 2>&1; then
    digest="$(printf '%s' "${value}" | shasum -a 256 | awk '{print $1}')"
  else
    _ecoda_die "sha256sum or shasum is required for artifact identity"
    return 1
  fi
  [[ "${digest}" =~ ^[[:xdigit:]]{64}$ ]] || {
    _ecoda_die "could not derive a SHA-256 artifact identity"
    return 1
  }
  printf '%s' "${digest}" | tr '[:upper:]' '[:lower:]'
}

# Canonicalize an absolute path while permitting the final artifact itself to
# be not-yet-created.  The containing directory must already exist; this is
# deliberate so ownership cannot be reserved for an unresolved path.
_ecoda_canonical_path() {
  local candidate="${1:-}" parent base parent_real boundary root
  local suffix="" current="${1:-}" part
  [[ -n "${candidate}" && "${candidate}" = /* ]] || return 1
  [[ "${candidate}" != *$'\n'* && "${candidate}" != *$'\t'* ]] || return 1
  # Apply symlink-component checks below the most-specific configured root,
  # while allowing trusted system prefixes and configured-root symlinks.
  for root in "${ECODA_RUNS_ROOT:-}" "${ECODA_OWNERS_ROOT:-}" \
              "${NAS_TARGET_DIR:-}" "${HPC_SCRATCH_DIR:-}"; do
    [[ -n "${root}" && "${root}" = /* ]] || continue
    boundary="${root%/}"
    [[ -n "${boundary}" ]] || boundary="/"
    case "${candidate}" in
      "${boundary}"|"${boundary}"/*)
        _ecoda_validate_path_ancestors "${candidate}" "${boundary}" || return 1
        break
        ;;
    esac
  done
  if [[ -e "${candidate}" || -L "${candidate}" ]]; then
    ecoda_realpath_existing "${candidate}"
    return $?
  fi
  # Resolve the nearest existing ancestor, then append the missing suffix.
  # This lets an owner be reserved before a producer creates its output
  # directory while rejecting symlinked parent components.
  while [[ ! -e "${current}" && ! -L "${current}" ]]; do
    part="$(basename "${current}")"
    [[ -n "${part}" && "${part}" != "." && "${part}" != ".." ]] || return 1
    if [[ -n "${suffix}" ]]; then
      suffix="${part}/${suffix}"
    else
      suffix="${part}"
    fi
    parent="$(dirname "${current}")"
    [[ "${parent}" != "${current}" ]] || return 1
    current="${parent}"
  done
  [[ -d "${current}" && ! -L "${current}" ]] || {
    _ecoda_die "artifact path has no existing directory ancestor: ${candidate}"
    return 1
  }
  parent_real="$(ecoda_realpath_existing "${current}")" || return 1
  if [[ -n "${suffix}" ]]; then
    printf '%s/%s' "${parent_real%/}" "${suffix}"
  else
    printf '%s' "${parent_real}"
  fi
}

ecoda_canonical_path() {
  _ecoda_canonical_path "$1"
}

_ecoda_artifact_record_load() {
  local record="$1" expected_path="$2" expected_producer="$3" expected_run="$4"
  local line key value index=0 expected_key
  local keys=(PATH SIZE MD5 RUN_ID PRODUCER STATE)
  [[ -f "${record}" && ! -L "${record}" && -r "${record}" && -s "${record}" ]] || {
    _ecoda_die "artifact record is missing or not a regular readable file: ${record}"
    return 1
  }
  [[ "$(tail -c 1 "${record}" 2>/dev/null; printf '\001')" == $'\n\001' ]] || {
    _ecoda_die "artifact record must end with a newline: ${record}"
    return 1
  }
  ECODA_ARTIFACT_RECORD_PATH="${record}"
  ECODA_ARTIFACT_RECORD_CANONICAL_PATH=""
  ECODA_ARTIFACT_RECORD_SIZE=""
  ECODA_ARTIFACT_RECORD_MD5=""
  ECODA_ARTIFACT_RECORD_RUN_ID=""
  ECODA_ARTIFACT_RECORD_PRODUCER=""
  ECODA_ARTIFACT_RECORD_STATE=""
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    [[ ${index} -le ${#keys[@]} ]] || {
      _ecoda_die "artifact record has extra fields: ${record}"
      return 1
    }
    expected_key="${keys[$((index - 1))]}"
    [[ "${line}" == "${expected_key}="* ]] || {
      _ecoda_die "artifact record field ${index} must be ${expected_key}: ${record}"
      return 1
    }
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" && -n "${value}" ]] || {
      _ecoda_die "artifact record has an empty or malformed ${expected_key}: ${record}"
      return 1
    }
    case "${key}" in
      PATH) ECODA_ARTIFACT_RECORD_CANONICAL_PATH="${value}" ;;
      SIZE) ECODA_ARTIFACT_RECORD_SIZE="${value}" ;;
      MD5) ECODA_ARTIFACT_RECORD_MD5="${value}" ;;
      RUN_ID) ECODA_ARTIFACT_RECORD_RUN_ID="${value}" ;;
      PRODUCER) ECODA_ARTIFACT_RECORD_PRODUCER="${value}" ;;
      STATE) ECODA_ARTIFACT_RECORD_STATE="${value}" ;;
    esac
  done < "${record}"
  [[ ${index} -eq ${#keys[@]} ]] || {
    _ecoda_die "artifact record has the wrong number of fields: ${record}"
    return 1
  }
  [[ "${ECODA_ARTIFACT_RECORD_CANONICAL_PATH}" == "${expected_path}" &&
     "${ECODA_ARTIFACT_RECORD_PRODUCER}" == "${expected_producer}" &&
     "${ECODA_ARTIFACT_RECORD_RUN_ID}" == "${expected_run}" &&
     "${ECODA_ARTIFACT_RECORD_STATE}" == "PUBLISHED" ]] || {
    _ecoda_die "artifact record identity/state mismatch: ${record}"
    return 1
  }
  [[ "${ECODA_ARTIFACT_RECORD_CANONICAL_PATH}" = /* &&
     "${ECODA_ARTIFACT_RECORD_CANONICAL_PATH}" != *$'\n'* &&
     "${ECODA_ARTIFACT_RECORD_CANONICAL_PATH}" != *$'\t'* ]] || return 1
  [[ "${ECODA_ARTIFACT_RECORD_SIZE}" =~ ^[1-9][0-9]*$ &&
     "${ECODA_ARTIFACT_RECORD_MD5}" =~ ^[[:xdigit:]]{32}$ ]] || {
    _ecoda_die "artifact record digest or size is malformed: ${record}"
    return 1
  }
  [[ "${ECODA_ARTIFACT_RECORD_RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || return 1
}

ecoda_artifact_record_path() {
  local path="${1:-}" run_id="${2:-}" canonical digest root root_real runs_root_real
  ecoda_validate_run_id "${run_id}" || return 1
  runs_root_real="$(_ecoda_runs_root_canonical 0)" || return 1
  root="${ECODA_RUNS_ROOT}/${run_id}"
  [[ -d "${root}" && ! -L "${root}" ]] || {
    _ecoda_die "producer run root is missing or is not a regular directory: ${root}"
    return 1
  }
  _ecoda_validate_path_ancestors "${root}" "${ECODA_RUNS_ROOT}" || return 1
  root_real="$(ecoda_realpath_existing "${root}")" || return 1
  [[ "${root_real}" == "${runs_root_real}/${run_id}" ]] || {
    _ecoda_die "producer run root has an invalid canonical layout: ${root}"
    return 1
  }
  canonical="$(_ecoda_canonical_path "${path}")" || {
    _ecoda_die "artifact path is missing or cannot be canonicalized: ${path}"
    return 1
  }
  digest="$(_ecoda_sha256_text "${canonical}")" || return 1
  ECODA_ARTIFACT_CANONICAL_PATH="${canonical}"
  ECODA_ARTIFACT_RECORD_PATH="${root}/manifests/artifacts/${digest:0:32}.record"
  printf '%s' "${ECODA_ARTIFACT_RECORD_PATH}"
}

ecoda_write_artifact_record() {
  local path="${1:-}" producer="${2:-}" run_id="${3:-}"
  local canonical record digest size content reused=0
  ecoda_validate_run_id "${run_id}" || return 1
  [[ -n "${producer}" && "${producer}" != *$'\n'* &&
     "${producer}" != *$'\t'* && "${producer}" != *'='* ]] || {
    _ecoda_die "artifact producer is empty or contains a record delimiter"
    return 1
  }
  canonical="$(_ecoda_canonical_path "${path}")" || {
    _ecoda_die "artifact path is missing or cannot be canonicalized: ${path}"
    return 1
  }
  # A worker may have just performed the strict full-file preflight.  Reuse
  # that immediately preceding digest only after checking its sidecar fields
  # and current size; otherwise this function is itself the first strict
  # publication boundary.
  if [[ "${ECODA_CHECKSUM_STRICT:-0}" == "1" &&
        ( "${ECODA_CHECKSUM_PATH:-}" == "${path}" ||
          "${ECODA_CHECKSUM_PATH:-}" == "${canonical}" ) &&
        "${ECODA_CHECKSUM_MD5:-}" =~ ^[[:xdigit:]]{32}$ &&
        "${ECODA_CHECKSUM_SIZE:-}" =~ ^[1-9][0-9]*$ ]]; then
    if ecoda_validate_checksum_record "${path}" \
        "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}"; then
      reused=1
    elif [[ "${path}" != "${canonical}" ]] &&
         ecoda_validate_checksum_record "${canonical}" \
           "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}" \
           "${path}.md5"; then
      reused=1
    fi
  fi
  if [[ ${reused} -eq 0 ]]; then
    ecoda_validate_checksum "${path}" || {
      _ecoda_die "artifact cannot be published without a strict checksum: ${path}"
      return 1
    }
  fi
  digest="${ECODA_CHECKSUM_MD5}"
  size="${ECODA_CHECKSUM_SIZE}"
  record="$(ecoda_artifact_record_path "${canonical}" "${run_id}")" || return 1
  content="PATH=${canonical}\nSIZE=${size}\nMD5=${digest}\nRUN_ID=${run_id}\nPRODUCER=${producer}\nSTATE=PUBLISHED\n"
  if [[ -e "${record}" || -L "${record}" ]]; then
    _ecoda_artifact_record_load "${record}" "${canonical}" "${producer}" "${run_id}" || return 1
    [[ "${ECODA_ARTIFACT_RECORD_MD5}" == "${digest}" &&
       "${ECODA_ARTIFACT_RECORD_SIZE}" == "${size}" ]] || {
      _ecoda_die "artifact record already exists with different content: ${record}"
      return 1
    }
    _ecoda_artifact_make_nonwritable "${canonical}" || return 1
    printf '%s' "${record}"
    return 0
  fi
  ecoda_atomic_write "${record}" "${content}" || return 1
  if ! _ecoda_artifact_make_nonwritable "${canonical}"; then
    if ! rm -f "${record}" >/dev/null 2>&1; then
      _ecoda_die "could not roll back writable artifact record: ${record}"
    fi
    return 1
  fi
  printf '%s' "${record}"
}

ecoda_validate_artifact_record() {
  local path="${1:-}" producer="${2:-}" run_id="${3:-}"
  local canonical record sidecar_ok=1 immutable=0
  ecoda_validate_run_id "${run_id}" || return 1
  [[ -n "${producer}" && "${producer}" != *$'\n'* &&
     "${producer}" != *$'\t'* && "${producer}" != *'='* ]] || return 1
  canonical="$(_ecoda_canonical_path "${path}")" || {
    _ecoda_die "artifact path is missing or cannot be canonicalized: ${path}"
    return 1
  }
  [[ -f "${canonical}" && ! -L "${canonical}" && -s "${canonical}" ]] || {
    _ecoda_die "artifact is missing or empty: ${canonical}"
    return 1
  }
  record="$(ecoda_artifact_record_path "${canonical}" "${run_id}")" || return 1
  _ecoda_artifact_record_load "${record}" "${canonical}" "${producer}" "${run_id}" || return 1
  if _ecoda_artifact_is_nonwritable "${canonical}" >/dev/null 2>&1; then
    immutable=1
  fi
  # Existing checksum sidecars retain their strict PATH semantics.  For a
  # symlinked caller path, accept either the caller spelling or its canonical
  # spelling, but never skip digest/size validation.
  if [[ ${immutable} -eq 1 ]]; then
    if ecoda_validate_checksum_record "${path}" \
        "${ECODA_ARTIFACT_RECORD_MD5}" "${ECODA_ARTIFACT_RECORD_SIZE}"; then
      sidecar_ok=0
    elif [[ "${path}" != "${canonical}" ]] &&
         ecoda_validate_checksum_record "${canonical}" \
           "${ECODA_ARTIFACT_RECORD_MD5}" "${ECODA_ARTIFACT_RECORD_SIZE}"; then
      sidecar_ok=0
    elif [[ "${path}" != "${canonical}" ]] &&
         ecoda_validate_checksum_record "${canonical}" \
           "${ECODA_ARTIFACT_RECORD_MD5}" "${ECODA_ARTIFACT_RECORD_SIZE}" \
           "${path}.md5"; then
      sidecar_ok=0
    fi
  else
    # Writable/legacy artifacts are not eligible for record-only reuse.  A
    # fresh full checksum catches same-size mutations before trusting the
    # published digest, while immutable artifacts retain the fast path above.
    if ecoda_validate_checksum "${path}" &&
       [[ "${ECODA_CHECKSUM_MD5}" == "${ECODA_ARTIFACT_RECORD_MD5}" &&
          "${ECODA_CHECKSUM_SIZE}" == "${ECODA_ARTIFACT_RECORD_SIZE}" ]]; then
      sidecar_ok=0
    elif [[ "${path}" != "${canonical}" ]] &&
         ecoda_validate_checksum "${canonical}" &&
         [[ "${ECODA_CHECKSUM_MD5}" == "${ECODA_ARTIFACT_RECORD_MD5}" &&
            "${ECODA_CHECKSUM_SIZE}" == "${ECODA_ARTIFACT_RECORD_SIZE}" ]]; then
      sidecar_ok=0
    elif [[ "${path}" != "${canonical}" ]] &&
         ecoda_validate_checksum "${canonical}" "${path}.md5" &&
         [[ "${ECODA_CHECKSUM_MD5}" == "${ECODA_ARTIFACT_RECORD_MD5}" &&
            "${ECODA_CHECKSUM_SIZE}" == "${ECODA_ARTIFACT_RECORD_SIZE}" ]]; then
      sidecar_ok=0
    fi
  fi
  [[ ${sidecar_ok} -eq 0 ]] || {
    _ecoda_die "artifact checksum sidecar does not match published record: ${canonical}"
    return 1
  }
  ECODA_ARTIFACT_CANONICAL_PATH="${canonical}"
  printf '%s' "${record}"
}

ecoda_require_source_script_path() {
  local candidate="${1:-}" source_root="${2:-}" candidate_real root_real
  [[ "${candidate}" = /* && "${source_root}" = /* ]] || {
    _ecoda_die "source script and source root must be absolute"
    return 1
  }
  [[ -f "${candidate}" && ! -L "${candidate}" && -r "${candidate}" ]] || {
    _ecoda_die "source script is missing, unreadable, or symlinked: ${candidate}"
    return 1
  }
  [[ -d "${source_root}" ]] || {
    _ecoda_die "source root is missing: ${source_root}"
    return 1
  }
  candidate_real="$(ecoda_realpath_existing "${candidate}")" || return 1
  root_real="$(ecoda_realpath_existing "${source_root}")" || return 1
  case "${candidate_real}" in
    "${root_real}"/*) ;;
    *) _ecoda_die "source script escapes immutable source root: ${candidate}"; return 1 ;;
  esac
  printf '%s' "${candidate_real}"
}

ecoda_artifact_owner_key() {
  local path="${1:-}" canonical digest key key_size basename safe_basename
  canonical="$(_ecoda_canonical_path "${path}")" || {
    _ecoda_die "artifact owner path is missing or cannot be canonicalized: ${path}"
    return 1
  }
  [[ "${canonical}" = /* ]] || return 1
  basename="$(basename "${canonical}")"
  [[ -n "${basename}" && "${basename}" != "." && "${basename}" != ".." ]] || return 1
  safe_basename="$(_ecoda_safe_component "${basename}")"
  [[ -n "${safe_basename}" ]] || return 1
  digest="$(_ecoda_sha256_text "${canonical}")" || return 1
  key="${safe_basename}_${digest:0:32}"
  key_size="$(printf '%s' "${key}" | wc -c | tr -d '[:space:]')" || return 1
  if [[ "${key_size}" =~ ^[0-9]+$ && ${key_size} -le 100 ]]; then
    printf '%s' "${key}"
    return 0
  fi
  printf 'artifact_%s' "${digest:0:32}"
}

ecoda_artifact_owner_dir() {
  local path="${1:-}" canonical key
  canonical="$(_ecoda_canonical_path "${path}")" || return 1
  key="$(ecoda_artifact_owner_key "${canonical}")" || return 1
  printf '%s/artifact/%s' "${ECODA_OWNERS_ROOT}" "${key}"
}

ECODA_ARTIFACT_OWNER_DIR=""
ECODA_ARTIFACT_OWNER_CANONICAL_PATH=""
ECODA_ARTIFACT_OWNER_ACQUIRED=0
ECODA_ARTIFACT_OWNER_STATE=""
ECODA_OUTPUT_PATHS=()
ECODA_OUTPUT_WRITE_FLAGS=()
ECODA_OUTPUT_OWNER_DIRS=()
_ecoda_init_output_arrays() {
  declare -p ECODA_OUTPUT_PATHS >/dev/null 2>&1 || ECODA_OUTPUT_PATHS=()
  declare -p ECODA_OUTPUT_WRITE_FLAGS >/dev/null 2>&1 || ECODA_OUTPUT_WRITE_FLAGS=()
  declare -p ECODA_OUTPUT_OWNER_DIRS >/dev/null 2>&1 || ECODA_OUTPUT_OWNER_DIRS=()
}

_ecoda_artifact_owner_validate_path() {
  local owner_dir="${1:-}" allow_missing="${2:-0}"
  local owners_root="${ECODA_OWNERS_ROOT:-}" artifact_root artifact_root_real owner_real
  artifact_root="${owners_root%/}/artifact"
  [[ "${owners_root}" = /* && "${owners_root}" != *$'\n'* &&
     "${owners_root}" != *$'\t'* ]] || {
    _ecoda_die "ECODA_OWNERS_ROOT must be an absolute path: ${owners_root}"
    return 1
  }
  [[ "${owner_dir}" = /* && "${owner_dir}" != *$'\n'* &&
     "${owner_dir}" != *$'\t'* ]] || {
    _ecoda_die "global artifact owner path is malformed: ${owner_dir}"
    return 1
  }
  case "${owner_dir}" in
    "${artifact_root}"/*) ;;
    *)
      _ecoda_die "global artifact owner is outside ${artifact_root}: ${owner_dir}"
      return 1
      ;;
  esac
  _ecoda_validate_path_ancestors "${owner_dir}" "${owners_root}" || return 1
  if [[ ! -e "${owner_dir}" && ! -L "${owner_dir}" ]]; then
    [[ "${allow_missing}" == "1" ]] && return 0
    _ecoda_die "global artifact owner directory is missing: ${owner_dir}"
    return 1
  fi
  [[ -d "${owner_dir}" && ! -L "${owner_dir}" ]] || {
    _ecoda_die "global artifact owner is not a regular directory: ${owner_dir}"
    return 1
  }
  artifact_root_real="$(ecoda_realpath_existing "${artifact_root}")" || {
    _ecoda_die "global artifact owner root cannot be canonicalized: ${artifact_root}"
    return 1
  }
  owner_real="$(ecoda_realpath_existing "${owner_dir}")" || {
    _ecoda_die "global artifact owner cannot be canonicalized: ${owner_dir}"
    return 1
  }
  case "${owner_real}" in
    "${artifact_root_real}"/*) ;;
    *)
      _ecoda_die "global artifact owner escapes ${artifact_root_real}: ${owner_dir}"
      return 1
      ;;
  esac
}

_ecoda_artifact_owner_validate_dir() {
  local owner_dir="${1:-}" expected_path="${2:-}"
  local line key value expected_key index=0 reason_seen=0 time_seen=0
  local owner_run owner_stage owner_path owner_state owner_pid
  local keys=(RUN_ID STAGE PATH STATE PID)
  _ecoda_artifact_owner_validate_path "${owner_dir}" || return 1
  [[ -f "${owner_dir}/owner" &&
     ! -L "${owner_dir}/owner" && -r "${owner_dir}/owner" ]] || {
    _ecoda_die "global artifact owner is missing or malformed: ${owner_dir}"
    return 1
  }
  [[ "$(tail -c 1 "${owner_dir}/owner" 2>/dev/null; printf '\001')" == $'\n\001' ]] || return 1
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    if [[ ${index} -le ${#keys[@]} ]]; then
      expected_key="${keys[$((index - 1))]}"
      [[ "${line}" == "${expected_key}="* ]] || return 1
    else
      case "${line}" in
        REASON=*) [[ ${reason_seen} -eq 0 ]] || return 1; reason_seen=1 ;;
        TIME=*) [[ ${time_seen} -eq 0 ]] || return 1; time_seen=1 ;;
        *) return 1 ;;
      esac
      expected_key="${line%%=*}"
    fi
    key="${line%%=*}"
    value="${line#*=}"
    if [[ "${key}" == REASON ]]; then
      [[ "${value}" != *$'\n'* && "${value}" != *$'\t'* &&
         "${value}" != *'='* ]] || return 1
      continue
    fi
    [[ -n "${value}" ]] || return 1
    case "${key}" in
      RUN_ID) owner_run="${value}" ;;
      STAGE) owner_stage="${value}" ;;
      PATH) owner_path="${value}" ;;
      STATE) owner_state="${value}" ;;
      PID) owner_pid="${value}" ;;
      TIME) [[ "${value}" =~ ^[0-9T:Z-]+$ ]] || return 1 ;;
      *) return 1 ;;
    esac
  done < "${owner_dir}/owner"
  [[ ${index} -ge ${#keys[@]} && ${index} -le 7 ]] || return 1
  [[ "${owner_run}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ &&
     -n "${owner_stage}" && "${owner_stage}" != *$'\n'* &&
     "${owner_stage}" != *$'\t'* && "${owner_stage}" != *'='* ]] || return 1
  [[ "${owner_path}" = /* && "${owner_path}" != *$'\n'* &&
     "${owner_path}" != *$'\t'* ]] || return 1
  [[ "${owner_state}" == ACTIVE || "${owner_state}" == OK ||
     "${owner_state}" == FAIL ]] || return 1
  [[ "${owner_pid}" =~ ^[0-9]+$ ]] || return 1
  [[ -n "${expected_path}" ]] || expected_path="${owner_path}"
  [[ "${owner_path}" == "${expected_path}" ]] || return 1
  ECODA_ARTIFACT_OWNER_DIR="${owner_dir}"
  ECODA_ARTIFACT_OWNER_CANONICAL_PATH="${owner_path}"
  ECODA_ARTIFACT_OWNER_STATE="${owner_state}"
  ECODA_ARTIFACT_OWNER_RUN="${owner_run}"
  ECODA_ARTIFACT_OWNER_STAGE="${owner_stage}"
  ECODA_ARTIFACT_OWNER_PID="${owner_pid}"
}

_ecoda_artifact_owner_set_state_dir() {
  local owner_dir="${1:-}" state="${2:-}" reason="${3:-}"
  local owner_run owner_stage owner_path owner_pid
  [[ "${state}" == "ACTIVE" || "${state}" == "OK" ||
     "${state}" == "FAIL" ]] || {
    _ecoda_die "invalid global artifact owner state: ${state}"
    return 1
  }
  _ecoda_artifact_owner_validate_dir "${owner_dir}" || return 1
  owner_run="${ECODA_ARTIFACT_OWNER_RUN}"
  owner_stage="${ECODA_ARTIFACT_OWNER_STAGE}"
  owner_path="${ECODA_ARTIFACT_OWNER_CANONICAL_PATH}"
  owner_pid="${ECODA_ARTIFACT_OWNER_PID}"
  [[ "${reason}" != *$'\n'* && "${reason}" != *$'\t'* &&
     "${reason}" != *'='* ]] || {
    _ecoda_die "global artifact owner reason contains a record delimiter"
    return 1
  }
  ecoda_atomic_write "${owner_dir}/owner" \
    "RUN_ID=${owner_run}\nSTAGE=${owner_stage}\nPATH=${owner_path}\nSTATE=${state}\nPID=${owner_pid}\nREASON=${reason}\nTIME=$(date -u +%Y-%m-%dT%H:%M:%SZ)\n"
  ECODA_ARTIFACT_OWNER_STATE="${state}"
}

ecoda_artifact_owner_acquire() {
  local path="${1:-}" stage="${2:-}" run_id="${3:-}"
  local allow_same_active="${4:-0}" reclaim_fail="${5:-1}" reclaim_terminal="${6:-0}"
  local canonical owner_dir state owner_run owner_stage owner_pid tombstone
  ECODA_ARTIFACT_OWNER_DIR=""
  ECODA_ARTIFACT_OWNER_CANONICAL_PATH=""
  ECODA_ARTIFACT_OWNER_ACQUIRED=0
  ECODA_ARTIFACT_OWNER_STATE=""
  ecoda_validate_run_id "${run_id}" || return 1
  [[ -n "${stage}" && "${stage}" != *$'\n'* &&
     "${stage}" != *$'\t'* && "${stage}" != *'='* ]] || {
    _ecoda_die "global artifact owner stage is empty or malformed"
    return 1
  }
  canonical="$(_ecoda_canonical_path "${path}")" || {
    _ecoda_die "global artifact path is missing or cannot be canonicalized: ${path}"
    return 1
  }
  owner_dir="$(ecoda_artifact_owner_dir "${canonical}")" || return 1
  _ecoda_artifact_owner_validate_path "${owner_dir}" 1 || return 1
  mkdir -p "$(dirname "${owner_dir}")" || return 1
  if mkdir "${owner_dir}" 2>/dev/null; then
    _ecoda_artifact_owner_validate_path "${owner_dir}" || {
      rmdir "${owner_dir}" 2>/dev/null || true
      return 1
    }
    if ! ecoda_atomic_write "${owner_dir}/owner" \
      "RUN_ID=${run_id}\nSTAGE=${stage}\nPATH=${canonical}\nSTATE=ACTIVE\nPID=$$\n"; then
      rm -f "${owner_dir}/owner"
      rmdir "${owner_dir}" 2>/dev/null || true
      return 1
    fi
    ECODA_ARTIFACT_OWNER_DIR="${owner_dir}"
    ECODA_ARTIFACT_OWNER_CANONICAL_PATH="${canonical}"
    ECODA_ARTIFACT_OWNER_ACQUIRED=1
    ECODA_ARTIFACT_OWNER_STATE="ACTIVE"
    ECODA_ARTIFACT_OWNER_RUN="${run_id}"
    ECODA_ARTIFACT_OWNER_STAGE="${stage}"
    ECODA_ARTIFACT_OWNER_PID="$$"
    printf '%s' "${owner_dir}"
    return 0
  fi
  _ecoda_artifact_owner_validate_dir "${owner_dir}" "${canonical}" || return 1
  state="${ECODA_ARTIFACT_OWNER_STATE}"
  owner_run="${ECODA_ARTIFACT_OWNER_RUN}"
  owner_stage="${ECODA_ARTIFACT_OWNER_STAGE}"
  owner_pid="${ECODA_ARTIFACT_OWNER_PID}"
  if [[ "${state}" == "ACTIVE" ]]; then
    if [[ "${allow_same_active}" == "1" && "${owner_run}" == "${run_id}" ]]; then
      ECODA_ARTIFACT_OWNER_DIR="${owner_dir}"
      ECODA_ARTIFACT_OWNER_CANONICAL_PATH="${canonical}"
      ECODA_ARTIFACT_OWNER_ACQUIRED=0
      ECODA_ARTIFACT_OWNER_STATE="${state}"
      ECODA_ARTIFACT_OWNER_RUN="${owner_run}"
      ECODA_ARTIFACT_OWNER_STAGE="${owner_stage}"
      ECODA_ARTIFACT_OWNER_PID="${owner_pid}"
      printf '%s' "${owner_dir}"
      return 0
    fi
    _ecoda_die "active global artifact owner ${owner_run} already owns ${canonical}"
    return 1
  fi
  if [[ "${state}" == "OK" && "${reclaim_terminal}" != "1" ]]; then
    _ecoda_die "terminal OK global artifact owner already exists: ${canonical}"
    return 2
  fi
  if [[ "${state}" == "FAIL" && "${reclaim_fail}" != "1" ]]; then
    _ecoda_die "terminal failed global artifact owner requires explicit reclaim: ${canonical}"
    return 2
  fi
  if [[ "${state}" == "OK" || "${state}" == "FAIL" ]]; then
    _ecoda_owner_reclaim_to_tombstone "${owner_dir}" || return 1
    tombstone="${ECODA_RECLAIM_TOMBSTONE}"
    if mkdir "${owner_dir}" 2>/dev/null; then
      if ! _ecoda_artifact_owner_validate_path "${owner_dir}"; then
        rmdir "${owner_dir}" 2>/dev/null || true
        if [[ ! -e "${owner_dir}" && ! -L "${owner_dir}" ]]; then
          _ecoda_owner_reclaim_tombstone_restore "${owner_dir}" "${tombstone}" ||
            _ecoda_owner_reclaim_tombstone_clean "${tombstone}" || true
        else
          _ecoda_owner_reclaim_tombstone_clean "${tombstone}" || true
        fi
        return 1
      fi
      if ! ecoda_atomic_write "${owner_dir}/owner" \
        "RUN_ID=${run_id}\nSTAGE=${stage}\nPATH=${canonical}\nSTATE=ACTIVE\nPID=$$\n"; then
        rm -f "${owner_dir}/owner.tmp.$$"
        if [[ ! -e "${owner_dir}/owner" && ! -L "${owner_dir}/owner" ]]; then
          rmdir "${owner_dir}" 2>/dev/null || true
        fi
        if [[ ! -e "${owner_dir}" && ! -L "${owner_dir}" ]]; then
          _ecoda_owner_reclaim_tombstone_restore "${owner_dir}" "${tombstone}" ||
            _ecoda_owner_reclaim_tombstone_clean "${tombstone}" || true
        else
          _ecoda_owner_reclaim_tombstone_clean "${tombstone}" || true
        fi
        return 1
      fi
      _ecoda_owner_reclaim_tombstone_clean "${tombstone}" || return 1
      ECODA_ARTIFACT_OWNER_DIR="${owner_dir}"
      ECODA_ARTIFACT_OWNER_CANONICAL_PATH="${canonical}"
      ECODA_ARTIFACT_OWNER_ACQUIRED=1
      ECODA_ARTIFACT_OWNER_STATE="ACTIVE"
      ECODA_ARTIFACT_OWNER_RUN="${run_id}"
      ECODA_ARTIFACT_OWNER_STAGE="${stage}"
      ECODA_ARTIFACT_OWNER_PID="$$"
    else
      if [[ -e "${owner_dir}" || -L "${owner_dir}" ]]; then
        _ecoda_owner_reclaim_tombstone_clean "${tombstone}" || true
      elif ! _ecoda_owner_reclaim_tombstone_restore \
          "${owner_dir}" "${tombstone}"; then
        _ecoda_owner_reclaim_tombstone_clean "${tombstone}" || true
      fi
      _ecoda_die "global artifact owner was claimed concurrently: ${owner_dir}"
      return 1
    fi
  fi
  printf '%s' "${owner_dir}"
}

ecoda_artifact_owner_set_state() {
  local path="${1:-}" state="${2:-}" reason="${3:-}" canonical owner_dir
  canonical="$(_ecoda_canonical_path "${path}")" || return 1
  owner_dir="$(ecoda_artifact_owner_dir "${canonical}")" || return 1
  _ecoda_artifact_owner_set_state_dir "${owner_dir}" "${state}" "${reason}"
}

ecoda_artifact_owner_validate() {
  local path="${1:-}" expected_run="${2:-}" canonical owner_dir
  canonical="$(_ecoda_canonical_path "${path}")" || return 1
  owner_dir="$(ecoda_artifact_owner_dir "${canonical}")" || return 1
  _ecoda_artifact_owner_validate_dir "${owner_dir}" "${canonical}" || return 1
  if [[ -n "${expected_run}" && "${ECODA_ARTIFACT_OWNER_RUN}" != "${expected_run}" ]]; then
    _ecoda_die "global artifact owner run mismatch for ${canonical}"
    return 1
  fi
  printf '%s' "${owner_dir}"
}

ecoda_artifact_owner_release() {
  ecoda_artifact_owner_set_state "$1" "${2:-OK}" "${3:-}"
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
  case "${owner_dir}" in
    "${ECODA_OWNERS_ROOT}/artifact/"*)
      _ecoda_artifact_owner_set_state_dir "${owner_dir}" "${state}" "${reason}"
      return $?
      ;;
  esac
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
_ecoda_output_add_path() {
  _ecoda_init_output_arrays
  local path="${1:-}" write_flag="${2:-1}"
  local canonical existing root root_real boundary lexically_rooted=0 rooted=0
  [[ -n "${path}" && "${path}" = /* &&
     "${path}" != *$'\n'* && "${path}" != *$'\t'* ]] || return 1
  [[ "${write_flag}" == 0 || "${write_flag}" == 1 ]] || return 1
  # Validate lexical components below each configured root before resolving
  # them.  The root itself may be a trusted symlink (for example a macOS
  # temporary directory), but an output subdirectory may not redirect away.
  for root in "${HPC_SCRATCH_DIR:-}" "${NAS_TARGET_DIR:-}"; do
    [[ -n "${root}" && "${root}" = /* && -d "${root}" ]] || continue
    boundary="${root%/}"
    [[ -n "${boundary}" ]] || boundary="/"
    case "${path}" in
      "${boundary}"|"${boundary}"/*)
        lexically_rooted=1
        _ecoda_validate_path_ancestors "${path}" "${boundary}" || return 1
        ;;
    esac
  done
  [[ ${lexically_rooted} -eq 1 ]] || {
    _ecoda_die "selected output is not lexically under a configured scratch/NAS root: ${path}"
    return 1
  }
  canonical="$(_ecoda_canonical_path "${path}")" || return 1
  [[ ! -d "${canonical}" ]] || return 1
  for root in "${HPC_SCRATCH_DIR:-}" "${NAS_TARGET_DIR:-}"; do
    [[ -n "${root}" && "${root}" = /* && -d "${root}" ]] || continue
    root_real="$(ecoda_realpath_existing "${root}")" || continue
    case "${canonical}" in
      "${root_real}"/*)
        _ecoda_validate_path_ancestors "${canonical}" "${root_real}" || return 1
        rooted=1
        break
        ;;
    esac
  done
  [[ ${rooted} -eq 1 ]] || {
    _ecoda_die "selected output escapes configured scratch/NAS roots: ${canonical}"
    return 1
  }
  if [[ ${#ECODA_OUTPUT_PATHS[@]} -gt 0 ]]; then
    for existing in "${ECODA_OUTPUT_PATHS[@]}"; do
      if [[ "${existing}" == "${canonical}" ]]; then
        _ecoda_die "selection expands to duplicate artifact path: ${canonical}"
        return 1
      fi
    done
  fi
  ECODA_OUTPUT_PATHS+=("${canonical}")
  ECODA_OUTPUT_WRITE_FLAGS+=("${write_flag}")
}

_ecoda_output_add_scratch_nas_pair() {
  local scratch_path="${1:-}" nas_path
  _ecoda_output_add_path "${scratch_path}" || return 1
  if [[ -n "${NAS_TARGET_DIR:-}" ]]; then
    [[ "${NAS_TARGET_DIR}" = /* ]] || return 1
    nas_path="${2:-}"
    [[ -n "${nas_path}" ]] || return 1
    _ecoda_output_add_path "${nas_path}" || return 1
  fi
}

_ecoda_stage5_artifacts_for() {
  local ds="$1" view="$2" label="$3" pass="${PASS_ARG:-${ANALYSIS_PASS:-}}"
  local root nas_root stem suffix n
  if [[ -z "${pass}" ]]; then
    case "${view}" in
      batch_effect_uncorrected) pass="uncorrected" ;;
      batch_effect_corrected) pass="corrected" ;;
    esac
  fi
  ECODA_BENCHMARK_ARTIFACTS=()
  if [[ -n "${ANALYSIS_ROOT:-}" ]]; then
    root="${ANALYSIS_ROOT}"
  elif [[ -n "${pass}" ]]; then
    root="${HPC_SCRATCH_DIR}/batch_effect/${pass}"
  else
    root="${HPC_SCRATCH_DIR}/benchmark"
  fi
  if [[ -n "${ANALYSIS_NAS_ROOT:-}" ]]; then
    nas_root="${ANALYSIS_NAS_ROOT}"
  elif [[ -n "${pass}" && -n "${NAS_TARGET_DIR:-}" ]]; then
    nas_root="${NAS_TARGET_DIR}/batch_effect/${pass}"
  elif [[ -n "${NAS_TARGET_DIR:-}" ]]; then
    nas_root="${NAS_TARGET_DIR}/benchmark"
  else
    nas_root=""
  fi
  if [[ "${label}" == "prepare_pseudobulk" ]]; then
    if [[ -n "${pass}" ]]; then
      ECODA_BENCHMARK_ARTIFACTS+=("${root}/pseudobulks/${ds}_batch_effect_${pass}_pseudobulk_hvg2000.rds")
    else
      for stem in schvg2000 hvg2000 hvg500 hvg2000_bl hvg1000 hvg3000; do
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/pseudobulks/${ds}_pseudobulk_${stem}.rds")
      done
    fi
  else
    case "${label}" in
      mrvi)
        if [[ -n "${pass}" ]]; then
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_batch_effect_${pass}_hvg2000_highres_mrvi_dists.feather")
        else
          for n in 1000 2000 3000; do
            ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_mrvi_dists.feather")
          done
        fi
        ;;
      scpoli)
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_lowres_scpoli_dims15_embs.feather")
        for n in 1000 3000; do
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_highres_scpoli_dims15_embs.feather")
        done
        for n in 2 3 5 10 15; do
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_highres_scpoli_dims${n}_embs.feather")
        done
        ;;
      pilot|qot)
        suffix="${label}"
        if [[ -n "${pass}" ]]; then
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_batch_effect_${pass}_hvg2000_highres_${suffix}_dists.feather")
        else
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_lowres_${suffix}_dists.feather")
          for n in 1000 2000 3000; do
            ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_highres_${suffix}_dists.feather")
          done
        fi
        ;;
      pilotgm)
        [[ -z "${pass}" ]] || return 1
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_highres_pilotgm_dists.feather")
        ;;
      trans|zeroimp)
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${ds}_${label}.rds")
        ;;
      gloscope|mofa|pseudobulk|scitd)
        stem="${ds}"
        [[ -z "${pass}" ]] || stem="${ds}_batch_effect_${pass}"
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_${label}.rds")
        ;;
      composition)
        stem="${ds}"
        [[ -z "${pass}" ]] || stem="${ds}_batch_effect_${pass}"
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_composition.rds")
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_metadata.rds")
        ;;
      *)
        # Explicit post-baseline methods use one method-specific result key.
        # The submitter's method registry remains authoritative for whether
        # the row is runnable; this fallback gives the ownership layer a
        # deterministic path without broad filesystem discovery.
        stem="${ds}"
        [[ -z "${pass}" ]] || stem="${ds}_batch_effect_${pass}"
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_${label}.rds")
        ;;
    esac
  fi
  ECODA_BENCHMARK_ARTIFACT_NAS=()
  if [[ -n "${nas_root}" ]]; then
    for stem in "${ECODA_BENCHMARK_ARTIFACTS[@]}"; do
      case "${stem}" in
        "${root}"/*) ECODA_BENCHMARK_ARTIFACT_NAS+=("${nas_root}/${stem#${root}/}") ;;
        *) return 1 ;;
      esac
    done
  fi
}

_ecoda_expand_output_selection() {
  local stage="${1:-}" selection="${2:-}" step script ds view label outputs extra
  local path name pass nas_path artifact_index write_flag
  local old_ifs
  _ecoda_init_output_arrays

  [[ -r "${selection}" && ! -L "${selection}" ]] || {
    _ecoda_die "selection manifest is missing or unreadable: ${selection}"
    return 1
  }
  ECODA_OUTPUT_PATHS=()
  ECODA_OUTPUT_WRITE_FLAGS=()
  ECODA_OUTPUT_OWNER_DIRS=()
  case "${stage}" in
    stage2)
      ecoda_validate_manifest "${selection}" 5 || return 1
      while IFS=$'\t' read -r step script outputs dependency owner extra; do
        [[ "${step}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ &&
           -n "${script}" && -n "${outputs}" && -n "${dependency}" &&
           -n "${owner}" && -z "${extra}" ]] || return 1
        if type step_script >/dev/null 2>&1; then
          expected_script="$(step_script "${step}" 2>/dev/null || true)"
          [[ -n "${expected_script}" && "${script}" == "${expected_script}" ]] || {
            _ecoda_die "Stage 2 script contract mismatch for ${step}"
            return 1
          }
        fi
        if type step_outputs >/dev/null 2>&1; then
          expected_outputs="$(step_outputs "${step}" 2>/dev/null || true)"
          [[ -n "${expected_outputs}" && "${outputs}" == "${expected_outputs}" ]] || {
            _ecoda_die "Stage 2 output contract mismatch for ${step}"
            return 1
          }
        fi
        write_flag=1
        [[ "${owner}" == "-" ]] && write_flag=0
        old_ifs="${IFS}"
        IFS=';'
        read -r -a ECODA_SELECTION_OUTPUTS <<< "${outputs}"
        IFS="${old_ifs}"
        [[ ${#ECODA_SELECTION_OUTPUTS[@]} -gt 0 ]] || return 1
        for path in "${ECODA_SELECTION_OUTPUTS[@]}"; do
          [[ -n "${path}" ]] || return 1
          _ecoda_output_add_path "${path}" "${write_flag}" || return 1
        done
      done < "${selection}"
      ;;
    stage3|stage4)
      ecoda_validate_manifest "${selection}" 2 || return 1
      [[ -r "${DATASETS_JSON_FILE:-}" ]] || {
        _ecoda_die "configured datasets.json is missing for ownership expansion"
        return 1
      }
      while IFS=$'\t' read -r ds view extra; do
        [[ -n "${ds}" && -n "${view}" && -z "${extra}" ]] || return 1
        ecoda_dataset_exists "${ds}" && ecoda_view_exists "${ds}" "${view}" || return 1
        name="$(ecoda_view_output_name "${ds}" "${view}")"
        [[ -n "${name}" && "${name}" != */* &&
           "${name}" != *$'\n'* && "${name}" != *$'\t'* ]] || return 1
        _ecoda_output_add_scratch_nas_pair \
          "${HPC_SCRATCH_DIR}/${ds}/output/${name}" \
          "${NAS_TARGET_DIR:-}/${ds}/output/${name}" || return 1
      done < "${selection}"
      ;;
    stage5)
      ecoda_validate_manifest "${selection}" 3 || return 1
      [[ -r "${DATASETS_JSON_FILE:-}" ]] || {
        _ecoda_die "configured datasets.json is missing for ownership expansion"
        return 1
      }
      while IFS=$'\t' read -r ds view label extra; do
        [[ -n "${ds}" && -n "${view}" && -n "${label}" && -z "${extra}" &&
           "${label}" =~ ^[A-Za-z0-9_.-]+$ ]] || return 1
        ecoda_dataset_exists "${ds}" && ecoda_view_exists "${ds}" "${view}" || return 1
        _ecoda_stage5_artifacts_for "${ds}" "${view}" "${label}" || return 1
        artifact_index=0
        for path in "${ECODA_BENCHMARK_ARTIFACTS[@]}"; do
          nas_path=""
          if [[ ${#ECODA_BENCHMARK_ARTIFACT_NAS[@]} -gt 0 ]]; then
            nas_path="${ECODA_BENCHMARK_ARTIFACT_NAS[${artifact_index}]}"
          fi
          artifact_index=$((artifact_index + 1))
          if [[ -n "${nas_path}" ]]; then
            _ecoda_output_add_scratch_nas_pair "${path}" "${nas_path}" || return 1
          else
            _ecoda_output_add_path "${path}" || return 1
          fi
        done
      done < "${selection}"
      ;;
    *) _ecoda_die "unsupported stage for ownership expansion: ${stage}"; return 1 ;;
  esac
  [[ ${#ECODA_OUTPUT_PATHS[@]} -gt 0 ]] || {
    _ecoda_die "selection expands to no output artifacts: ${selection}"
    return 1
  }
}

_ecoda_validate_output_owner_availability() {
  local run_id="$1" path owner_dir state owner_run
  _ecoda_init_output_arrays

  if [[ ${#ECODA_OUTPUT_PATHS[@]} -gt 0 ]]; then
    for path in "${ECODA_OUTPUT_PATHS[@]}"; do
      owner_dir="$(ecoda_artifact_owner_dir "${path}")" || return 1
      if [[ -e "${owner_dir}" || -L "${owner_dir}" ]]; then
        _ecoda_artifact_owner_validate_dir "${owner_dir}" "${path}" || return 1
        state="${ECODA_ARTIFACT_OWNER_STATE}"
        owner_run="${ECODA_ARTIFACT_OWNER_RUN}"
        if [[ "${state}" == "ACTIVE" && "${owner_run}" != "${run_id}" ]]; then
          _ecoda_die "active global artifact owner ${owner_run} overlaps ${path}"
          return 1
        fi
      fi
    done
  fi
}

ecoda_validate_output_ownership() {
  local stage="${1:-}" selection="${2:-}" run_id="${3:-}"
  local reclaim_terminal="${4:-0}"
  local path owner_dir acquired_dir existing_owner write_flag rc=0 owner_listed=0
  local output_index=0
  local newly_acquired=()
  _ecoda_init_output_arrays

  [[ "${reclaim_terminal}" == "0" || "${reclaim_terminal}" == "1" ]] || {
    _ecoda_die "artifact owner terminal reclaim flag must be 0 or 1"
    return 1
  }
  ecoda_validate_run_id "${run_id}" || return 1
  [[ "${stage}" == stage2 || "${stage}" == stage3 ||
     "${stage}" == stage4 || "${stage}" == stage5 ]] || {
    _ecoda_die "unsupported stage for output ownership: ${stage}"
    return 1
  }
  [[ "${selection}" = /* ]] || {
    _ecoda_die "ownership selection must be an absolute path"
    return 1
  }
  _ecoda_expand_output_selection "${stage}" "${selection}" || return 1
  _ecoda_validate_output_owner_availability "${run_id}" || return 1
  # Reserve only after the complete selection has been expanded and checked.
  if [[ ${#ECODA_OUTPUT_PATHS[@]} -gt 0 ]]; then
    for path in "${ECODA_OUTPUT_PATHS[@]}"; do
      write_flag="${ECODA_OUTPUT_WRITE_FLAGS[${output_index}]:-1}"
      output_index=$((output_index + 1))
      [[ "${write_flag}" == "1" ]] || continue
      if ! ecoda_artifact_owner_acquire "${path}" "${stage}" "${run_id}" 1 \
          "${reclaim_terminal}" "${reclaim_terminal}" >/dev/null; then
        rc=1
        break
      fi
      owner_dir="${ECODA_ARTIFACT_OWNER_DIR}"
      [[ -n "${owner_dir}" ]] || { rc=1; break; }
      if [[ "${ECODA_ARTIFACT_OWNER_ACQUIRED}" == "1" ]]; then
        newly_acquired+=("${owner_dir}")
      fi
      ecoda_owner_track "${owner_dir}" || { rc=1; break; }
      owner_listed=0
      for existing_owner in "${ECODA_OUTPUT_OWNER_DIRS[@]:-}"; do
        [[ -n "${existing_owner}" ]] || continue
        if [[ "${existing_owner}" == "${owner_dir}" ]]; then
          owner_listed=1
          break
        fi
      done
      [[ ${owner_listed} -eq 1 ]] ||
        ECODA_OUTPUT_OWNER_DIRS+=("${owner_dir}")
    done
  else
    _ecoda_die "selection expands to no output artifacts: ${selection}"
    return 1
  fi
  if [[ ${rc} -ne 0 ]]; then
    if [[ ${#newly_acquired[@]} -gt 0 ]]; then
      for acquired_dir in "${newly_acquired[@]}"; do
        ecoda_owner_set_state "${acquired_dir}" FAIL \
          "output ownership acquisition failed" >/dev/null 2>&1 || true
      done
    fi
    return 1
  fi
}

ecoda_require_input_ownership() {
  local path="${1:-}" run_id="${2:-}" canonical owner_dir
  ecoda_validate_run_id "${run_id}" || return 1
  canonical="$(_ecoda_canonical_path "${path}")" || {
    _ecoda_die "input path is missing or cannot be canonicalized: ${path}"
    return 1
  }
  owner_dir="$(ecoda_artifact_owner_dir "${canonical}")" || return 1
  [[ -e "${owner_dir}" || -L "${owner_dir}" ]] || return 0
  _ecoda_artifact_owner_validate_dir "${owner_dir}" "${canonical}" || return 1
  if [[ "${ECODA_ARTIFACT_OWNER_STATE}" == "ACTIVE" ]]; then
    _ecoda_die "input path has an active writer ${ECODA_ARTIFACT_OWNER_RUN}: ${canonical}"
    return 1
  fi
}

_ecoda_validate_input_schema_contract() {
  local path="$1" validator python_bin
  case "${path}" in
    *.h5ad)
      validator="${ECODA_ARTIFACT_CONTRACT_PY:-${PROJECT_ROOT:-}/src/utils/py/artifact_contract.py}"
      python_bin="${PYTHON_BIN:-}"
      [[ -r "${validator}" && -n "${python_bin}" ]] || {
        _ecoda_die "H5AD schema validator is unavailable: ${path}"
        return 1
      }
      "${python_bin}" "${validator}" --path "${path}" --kind h5ad >/dev/null 2>&1 || return 1
      ;;
    *) ;;
  esac
}

ecoda_validate_input_artifact() {
  local path="${1:-}" producer="${2:-}" run_id="${3:-}" canonical owner_dir
  ecoda_validate_run_id "${run_id}" || return 1
  canonical="$(_ecoda_canonical_path "${path}")" || return 1
  owner_dir="$(ecoda_artifact_owner_dir "${canonical}")" || return 1
  [[ -d "${owner_dir}" ]] || {
    _ecoda_die "input artifact has no global owner: ${canonical}"
    return 1
  }
  _ecoda_artifact_owner_validate_dir "${owner_dir}" "${canonical}" || return 1
  [[ "${ECODA_ARTIFACT_OWNER_STATE}" == "OK" ]] || {
    _ecoda_die "input artifact writer is not terminal OK: ${canonical}"
    return 1
  }
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${run_id}" ]] || {
    _ecoda_die "input artifact owner run does not match producer run: ${canonical}"
    return 1
  }
  ecoda_validate_artifact_record "${path}" "${producer}" "${run_id}" || return 1
  _ecoda_validate_input_schema_contract "${canonical}" || return 1
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
    Kidney_KPMP_full
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
    rows="$(sacct -j "${job}" -n -P --format=JobID,State,ExitCode 2>/dev/null || true)"
    scheduler_active=0
    if command -v squeue >/dev/null 2>&1; then
      active_jobs="$(squeue -j "${job}" -h -o "%A" 2>/dev/null || true)"
      while IFS= read -r active_id; do
        case "${active_id}" in
          "${job}"|"${job}"_*) scheduler_active=1; break ;;
        esac
      done <<< "${active_jobs}"
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
