#!/usr/bin/env bash
# Self-verifying, commit-keyed source snapshots for ECODA workers.
#
# This file deliberately has no dependency on the mutable checkout.  The copy
# in a published snapshot is the bootstrap used by subsequent workers.

set -u -o pipefail

PROGRAM="${0##*/}"
SOURCE_MANIFEST_FORMAT="1"
DEFAULT_SCGATE_DB_BRANCH="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4"

usage() {
  cat >&2 <<'EOF'
Usage:
  ecoda_source_snapshot.sh create \
    --source-root ABSOLUTE_REPOSITORY \
    --snapshot-parent ABSOLUTE_PARENT \
    --commit FULL_COMMIT

  ecoda_source_snapshot.sh exec \
    --source-root ABSOLUTE_SNAPSHOT_TREE \
    --source-manifest ABSOLUTE_MANIFEST \
    --host-env-prefix ABSOLUTE_ENV_PREFIX \
    --runtime-image ABSOLUTE_IMAGE \
    --runtime-manifest ABSOLUTE_IMAGE_MANIFEST \
    --run-id RUN_ID \
    --scratch-root ABSOLUTE_SCRATCH_ROOT \
    --logs-root ABSOLUTE_LOG_ROOT \
    --script RELATIVE_SCRIPT -- [SCRIPT_ARGS...]
EOF
  return 2
}

error() {
  printf 'ERROR: %s\n' "$*" >&2
  return 1
}

is_safe_value() {
  local value="${1:-}"
  [[ -n "${value}" ]] || return 1
  [[ "${value}" != *$'\n'* && "${value}" != *$'\r'* && \
    "${value}" != *$'\t'* && "${value}" != *[[:space:]]* && \
    "${value}" != *"="* ]] || return 1
  return 0
}

require_absolute_argument() {
  local label="${1:-path}"
  local value="${2:-}"
  [[ "${value}" = /* ]] || {
    error "${label} must be an absolute path: ${value}"
    return 1
  }
  is_safe_value "${value}" || {
    error "${label} contains unsafe whitespace, control characters, or '=': ${value}"
    return 1
  }
  case "${value}" in
    */../*|*/..|../*|..)
      error "${label} contains a parent-directory escape: ${value}"
      return 1
      ;;
  esac
  return 0
}

no_symlink_components() {
  local path="${1:-}"
  local component prefix="/"
  local old_ifs="${IFS}"
  local parts=()
  [[ "${path}" = /* ]] || return 1
  IFS='/' read -r -a parts <<< "${path#/}"
  IFS="${old_ifs}"
  for component in "${parts[@]}"; do
    [[ -n "${component}" && "${component}" != "." ]] || continue
    [[ "${component}" != ".." ]] || return 1
    prefix="${prefix%/}/${component}"
    [[ ! -L "${prefix}" ]] || return 1
  done
  return 0
}

canonical_existing() {
  local path="${1:-}"
  local resolved
  no_symlink_components "${path}" || {
    error "path contains a symlink or escape: ${path}"
    return 1
  }
  if command -v realpath >/dev/null 2>&1; then
    if resolved="$(realpath -e "${path}" 2>/dev/null)" && [[ -n "${resolved}" ]]; then
      printf '%s\n' "${resolved}"
      return 0
    fi
    if resolved="$(realpath "${path}" 2>/dev/null)" && [[ -e "${resolved}" ]]; then
      printf '%s\n' "${resolved}"
      return 0
    fi
  fi
  if [[ -d "${path}" ]] && resolved="$(cd -- "${path}" 2>/dev/null && pwd -P)" && [[ -n "${resolved}" ]]; then
    printf '%s\n' "${resolved}"
    return 0
  fi
  error "path does not exist or cannot be canonicalized: ${path}"
  return 1
}


sha256_file() {
  local path="${1:-}"
  [[ -f "${path}" && ! -L "${path}" ]] || {
    error "cannot hash missing or non-regular file: ${path}"
    return 1
  }
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${path}" | awk '{print $1}'
    return "${PIPESTATUS[0]}"
  fi
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${path}" | awk '{print $1}'
    return "${PIPESTATUS[0]}"
  fi
  error "sha256sum or shasum is required"
  return 1
}

file_mode() {
  local path="${1:-}"
  local mode
  if mode="$(stat -c '%a' "${path}" 2>/dev/null)" && [[ -n "${mode}" ]]; then
    printf '%s\n' "${mode}"
    return 0
  fi
  if mode="$(stat -f '%Lp' "${path}" 2>/dev/null)" && [[ -n "${mode}" ]]; then
    printf '%s\n' "${mode}"
    return 0
  fi
  error "cannot inspect permissions: ${path}"
  return 1
}

mode_is_readonly() {
  local mode="${1:-}"
  local last_three="${mode: -3}"
  [[ "${last_three}" != *[2367]* ]]
}

check_readonly_tree() {
  local root="${1:-}"
  local path mode
  [[ -d "${root}" && ! -L "${root}" ]] || {
    error "read-only check requires a real directory: ${root}"
    return 1
  }
  while IFS= read -r path; do
    [[ ! -L "${path}" ]] || {
      error "snapshot contains a symlink: ${path}"
      return 1
    }
    [[ -r "${path}" ]] || {
      error "published snapshot is not readable: ${path}"
      return 1
    }
    if [[ -d "${path}" ]]; then
      [[ -x "${path}" ]] || {
        error "published snapshot directory is not searchable: ${path}"
        return 1
      }
    fi
    mode="$(file_mode "${path}")" || return 1
    mode_is_readonly "${mode}" || {
      error "published snapshot is writable: ${path}"
      return 1
    }
  done < <(find "${root}" -print 2>/dev/null)
  return 0
}

check_no_symlinks_or_git() {
  local root="${1:-}"
  local found
  found="$(find "${root}" -type l -print -quit 2>/dev/null)"
  [[ -z "${found}" ]] || {
    error "snapshot contains a symlink: ${found}"
    return 1
  }
  found="$(find "${root}" -name .git -print -quit 2>/dev/null)"
  [[ -z "${found}" ]] || {
    error "snapshot contains a .git entry: ${found}"
    return 1
  }
  return 0
}

path_is_within() {
  local candidate="${1:-}"
  local parent="${2:-}"
  [[ "${candidate}" == "${parent}" || "${candidate}" == "${parent}/"* ]]
}
path_identity() {
  local path="${1:-}"
  local identity
  [[ -e "${path}" && ! -L "${path}" ]] || {
    error "cannot inspect missing or symlinked path identity: ${path}"
    return 1
  }
  if identity="$(stat -c '%d:%i' "${path}" 2>/dev/null)" &&
     [[ "${identity}" =~ ^[0-9]+:[0-9]+$ ]]; then
    printf '%s\n' "${identity}"
    return 0
  fi
  if identity="$(stat -f '%d:%i' "${path}" 2>/dev/null)" &&
     [[ "${identity}" =~ ^[0-9]+:[0-9]+$ ]]; then
    printf '%s\n' "${identity}"
    return 0
  fi
  error "cannot inspect path identity: ${path}"
  return 1
}

verify_snapshot_identity_now() {
  local source_root="${1:-}"
  local snapshot_root="${2:-}"
  local snapshot_parent="${3:-}"
  local source_manifest="${4:-}"
  local expected_parent_id="${5:-}"
  local expected_snapshot_id="${6:-}"
  local expected_tree_id="${7:-}"
  local expected_identity_dir_id="${8:-}"
  local expected_manifest_id="${9:-}"
  local expected_archive_id="${10:-}"
  local expected_marker_id="${11:-}"
  local snapshot_id="${snapshot_root##*/}"
  local canonical_parent canonical_snapshot canonical_tree canonical_manifest
  local current_parent_id current_snapshot_id current_tree_id current_identity_dir_id
  local current_manifest_id current_archive_id current_marker_id
  local identity_dir="${snapshot_root}/identity"
  local archive="${identity_dir}/source.tar"
  local marker="${snapshot_root}/COMPLETE"

  [[ "${snapshot_parent}" != "/" && -n "${snapshot_parent}" &&
    "${snapshot_root}" == "${snapshot_parent}/${snapshot_id}" &&
    "${source_root}" == "${snapshot_root}/tree" &&
    "${snapshot_id}" =~ ^[0-9a-fA-F]{40}$ &&
    "${source_manifest}" == "${identity_dir}/source.manifest" ]] || {
    error "snapshot identity is not an exact parent/commit/tree layout"
    return 1
  }
  no_symlink_components "${snapshot_parent}" || {
    error "snapshot parent ancestry contains a symlink: ${snapshot_parent}"
    return 1
  }
  no_symlink_components "${snapshot_root}" || {
    error "snapshot ancestry contains a symlink: ${snapshot_root}"
    return 1
  }
  no_symlink_components "${source_root}" || {
    error "snapshot tree ancestry contains a symlink: ${source_root}"
    return 1
  }
  no_symlink_components "${source_manifest}" || {
    error "snapshot manifest ancestry contains a symlink: ${source_manifest}"
    return 1
  }
  [[ -d "${snapshot_parent}" && ! -L "${snapshot_parent}" &&
    -d "${snapshot_root}" && ! -L "${snapshot_root}" &&
    -d "${source_root}" && ! -L "${source_root}" &&
    -d "${identity_dir}" && ! -L "${identity_dir}" &&
    -f "${source_manifest}" && ! -L "${source_manifest}" &&
    -f "${archive}" && ! -L "${archive}" &&
    -f "${marker}" && ! -L "${marker}" ]] || {
    error "snapshot identity disappeared or changed type"
    return 1
  }
  canonical_parent="$(canonical_existing "${snapshot_parent}")" || return 1
  canonical_snapshot="$(canonical_existing "${snapshot_root}")" || return 1
  canonical_tree="$(canonical_existing "${source_root}")" || return 1
  canonical_manifest="$(canonical_existing "${source_manifest}")" || return 1
  [[ "${canonical_parent}" == "${snapshot_parent}" &&
    "${canonical_snapshot}" == "${snapshot_root}" &&
    "${canonical_tree}" == "${source_root}" &&
    "${canonical_manifest}" == "${source_manifest}" ]] || {
    error "snapshot identity canonical path changed"
    return 1
  }
  current_parent_id="$(path_identity "${snapshot_parent}")" || return 1
  current_snapshot_id="$(path_identity "${snapshot_root}")" || return 1
  current_tree_id="$(path_identity "${source_root}")" || return 1
  current_identity_dir_id="$(path_identity "${identity_dir}")" || return 1
  current_manifest_id="$(path_identity "${source_manifest}")" || return 1
  current_archive_id="$(path_identity "${archive}")" || return 1
  current_marker_id="$(path_identity "${marker}")" || return 1
  [[ "${current_parent_id}" == "${expected_parent_id}" &&
    "${current_snapshot_id}" == "${expected_snapshot_id}" &&
    "${current_tree_id}" == "${expected_tree_id}" &&
    "${current_identity_dir_id}" == "${expected_identity_dir_id}" &&
    "${current_manifest_id}" == "${expected_manifest_id}" &&
    "${current_archive_id}" == "${expected_archive_id}" &&
    "${current_marker_id}" == "${expected_marker_id}" ]] || {
    error "snapshot identity was replaced during validation"
    return 1
  }
  return 0
}
run_snapshot_with_parent_lock() {
  local snapshot_parent="${1:-}"
  local snapshot_parent_mode="${2:-}"
  local source_root="${3:-}"
  local snapshot_root="${4:-}"
  local source_manifest="${5:-}"
  local script_path="${6:-}"
  local expected_parent_id="${7:-}"
  local expected_snapshot_id="${8:-}"
  local expected_tree_id="${9:-}"
  local expected_identity_dir_id="${10:-}"
  local expected_manifest_id="${11:-}"
  local expected_archive_id="${12:-}"
  local expected_marker_id="${13:-}"
  local expected_script_id="${14:-}"
  shift 14

  # The parent stays writable between snapshot creations, but not while a
  # worker is being launched.  Locking only for this invocation closes the
  # ordinary rename/unlink race without permanently disabling future create
  # operations.  Components above this parent are never followed through
  # symlinks and remain the caller's immutable-ancestry responsibility.
  # Entering the verified directory and chmod'ing "." avoids following a
  # last-moment replacement symlink; exact path and inode checks are repeated
  # after the lock and immediately before launch.
  (
    snapshot_parent_lock_path="${snapshot_parent}"
    snapshot_parent_lock_mode="${snapshot_parent_mode}"
    snapshot_parent_lock_dir="${snapshot_parent_lock_path}/.ecoda-exec-lock"
    snapshot_parent_lock_claimed=0
    snapshot_parent_lock_acquired=0
    trap 'if [[ "${snapshot_parent_lock_claimed}" -eq 1 ]]; then
      if [[ "${snapshot_parent_lock_acquired}" -eq 1 ]] &&
         ! chmod "${snapshot_parent_lock_mode}" . >/dev/null 2>&1; then
        printf "ERROR: could not restore snapshot parent permissions: %s\n" "${snapshot_parent_lock_path}" >&2
      fi
      if ! rmdir .ecoda-exec-lock >/dev/null 2>&1; then
        printf "ERROR: could not release snapshot parent lock: %s\n" "${snapshot_parent_lock_dir}" >&2
      fi
    fi' EXIT

    local current_parent_mode locked_parent_mode current_parent_id current_script_id
    cd -- "${snapshot_parent_lock_path}" || {
      error "cannot enter snapshot parent before locking: ${snapshot_parent_lock_path}"
      exit 1
    }
    current_parent_mode="$(file_mode .)" || exit 1
    current_parent_id="$(path_identity .)" || exit 1
    [[ "${current_parent_mode}" == "${snapshot_parent_lock_mode}" &&
      "${current_parent_id}" == "${expected_parent_id}" ]] || {
      error "snapshot parent changed during validation: ${snapshot_parent_lock_path}"
      exit 1
    }
    verify_snapshot_identity_now "${source_root}" "${snapshot_root}" \
      "${snapshot_parent}" "${source_manifest}" \
      "${expected_parent_id}" "${expected_snapshot_id}" "${expected_tree_id}" \
      "${expected_identity_dir_id}" "${expected_manifest_id}" \
      "${expected_archive_id}" "${expected_marker_id}" || exit 1
    no_symlink_components "${snapshot_parent_lock_path}" || {
      error "snapshot parent ancestry changed before execution: ${snapshot_parent_lock_path}"
      exit 1
    }
    [[ ! -e ".ecoda-exec-lock" && ! -L ".ecoda-exec-lock" ]] || {
      error "snapshot parent execution lock is already held: ${snapshot_parent_lock_dir}"
      exit 1
    }
    if mode_is_readonly "${current_parent_mode}"; then
      :
    else
      mkdir .ecoda-exec-lock || {
        error "cannot claim snapshot parent execution lock: ${snapshot_parent_lock_dir}"
        exit 1
      }
      snapshot_parent_lock_claimed=1
      chmod a-w . || {
        error "cannot temporarily lock snapshot parent for execution: ${snapshot_parent_lock_path}"
        exit 1
      }
      snapshot_parent_lock_acquired=1
    fi
    locked_parent_mode="$(file_mode .)" || exit 1
    mode_is_readonly "${locked_parent_mode}" || {
      error "snapshot parent remained writable during execution setup: ${snapshot_parent_lock_path}"
      exit 1
    }
    current_parent_id="$(path_identity .)" || exit 1
    [[ "${current_parent_id}" == "${expected_parent_id}" ]] || {
      error "snapshot parent changed while being locked: ${snapshot_parent_lock_path}"
      exit 1
    }
    verify_snapshot_identity_now "${source_root}" "${snapshot_root}" \
      "${snapshot_parent}" "${source_manifest}" \
      "${expected_parent_id}" "${expected_snapshot_id}" "${expected_tree_id}" \
      "${expected_identity_dir_id}" "${expected_manifest_id}" \
      "${expected_archive_id}" "${expected_marker_id}" || exit 1
    validate_snapshot_paths_and_manifest "${snapshot_root}" \
      "${snapshot_root##*/}" || exit 1
    no_symlink_components "${script_path}" || {
      error "snapshot worker path changed before execution: ${script_path}"
      exit 1
    }
    current_script_id="$(path_identity "${script_path}")" || exit 1
    [[ "${current_script_id}" == "${expected_script_id}" ]] || {
      error "snapshot worker path was replaced before execution: ${script_path}"
      exit 1
    }
    verify_snapshot_identity_now "${source_root}" "${snapshot_root}" \
      "${snapshot_parent}" "${source_manifest}" \
      "${expected_parent_id}" "${expected_snapshot_id}" "${expected_tree_id}" \
      "${expected_identity_dir_id}" "${expected_manifest_id}" \
      "${expected_archive_id}" "${expected_marker_id}" || exit 1
    current_parent_id="$(path_identity .)" || exit 1
    [[ "${current_parent_id}" == "${expected_parent_id}" ]] || {
      error "snapshot parent changed immediately before execution: ${snapshot_parent_lock_path}"
      exit 1
    }
    /bin/bash "${script_path}" "$@"
  )
}



compare_extracted_trees() {
  local left="${1:-}"
  local right="${2:-}"
  local left_list right_list
  left_list="$(cd -- "${left}" && find . -print | LC_ALL=C sort)" || return 1
  right_list="$(cd -- "${right}" && find . -print | LC_ALL=C sort)" || return 1
  [[ "${left_list}" == "${right_list}" ]] || {
    error "archive extraction file lists differ"
    return 1
  }
  diff -r -q "${left}" "${right}" >/dev/null 2>&1 || {
    error "archive extraction contents differ"
    return 1
  }
  return 0
}

archive_entries_safe() {
  local archive="${1:-}"
  local entries line first
  entries="$(tar -tf "${archive}" 2>/dev/null)" || {
    error "cannot list source archive: ${archive}"
    return 1
  }
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -n "${line}" ]] || continue
    [[ "${line}" != *$'\n'* && "${line}" != *$'\r'* ]] || {
      error "source archive contains an unsafe filename"
      return 1
    }
    line="${line%/}"
    [[ "${line}" != /* ]] || {
      error "source archive contains an absolute path: ${line}"
      return 1
    }
    case "/${line}/" in
      */../*)
        error "source archive contains a parent-directory escape: ${line}"
        return 1
        ;;
    esac
  done <<< "${entries}"

  # A symlink in an archive could make a later archive member escape while it
  # is being extracted.  Git trees used by create reject symlinks too, but
  # exec must defend against a retained archive that was modified later.
  entries="$(tar -tvf "${archive}" 2>/dev/null)" || {
    error "cannot inspect source archive entries: ${archive}"
    return 1
  }
  while IFS= read -r line || [[ -n "${line}" ]]; do
    first="${line#"${line%%[![:space:]]*}"}"
    case "${first}" in
      l*)
        error "source archive contains a symlink"
        return 1
        ;;
    esac
  done <<< "${entries}"
  return 0
}

extract_and_compare_archive() {
  local archive="${1:-}"
  local tree="${2:-}"
  local fresh
  archive_entries_safe "${archive}" || return 1
  fresh="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-source-verify.XXXXXXXX")" || {
    error "cannot create temporary archive extraction directory"
    return 1
  }
  if ! tar -xf "${archive}" -C "${fresh}" 2>/dev/null; then
    chmod -R u+w "${fresh}" >/dev/null 2>&1 || true
    rm -rf "${fresh}"
    error "cannot extract source archive"
    return 1
  fi
  if ! check_no_symlinks_or_git "${fresh}" || ! compare_extracted_trees "${tree}" "${fresh}"; then
    chmod -R u+w "${fresh}" >/dev/null 2>&1 || true
    rm -rf "${fresh}"
    return 1
  fi
  chmod -R u+w "${fresh}" >/dev/null 2>&1 || true
  rm -rf "${fresh}"
  return 0
}

source_manifest_load() {
  local manifest="${1:-}"
  local expected_key line key value index=0
  local keys=(
    FORMAT SOURCE_ROOT SOURCE_COMMIT SOURCE_ARCHIVE_PATH
    SOURCE_ARCHIVE_SHA256 CONFIG_HELPER_SHA256 DATASETS_SHA256
    PIXI_TOML_SHA256 PIXI_LOCK_SHA256 AUX_ROOT SCGATE_DB_BRANCH
  )
  [[ -f "${manifest}" && ! -L "${manifest}" && -r "${manifest}" ]] || {
    error "source manifest is missing, unreadable, or not regular: ${manifest}"
    return 1
  }
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    [[ "${index}" -le "${#keys[@]}" ]] || {
      error "source manifest has extra fields: ${manifest}"
      return 1
    }
    expected_key="${keys[$((index - 1))]}"
    [[ "${line}" == "${expected_key}="* ]] || {
      error "source manifest field ${index} must be ${expected_key}: ${manifest}"
      return 1
    }
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" && -n "${value}" ]] || {
      error "source manifest has an empty or malformed ${expected_key}: ${manifest}"
      return 1
    }
    is_safe_value "${value}" || {
      error "source manifest has an unsafe value for ${expected_key}: ${manifest}"
      return 1
    }
    case "${expected_key}" in
      FORMAT) SOURCE_FORMAT="${value}" ;;
      SOURCE_ROOT) SOURCE_ROOT_VALUE="${value}" ;;
      SOURCE_COMMIT) SOURCE_COMMIT_VALUE="${value}" ;;
      SOURCE_ARCHIVE_PATH) SOURCE_ARCHIVE_PATH_VALUE="${value}" ;;
      SOURCE_ARCHIVE_SHA256) SOURCE_ARCHIVE_SHA256_VALUE="${value}" ;;
      CONFIG_HELPER_SHA256) CONFIG_HELPER_SHA256_VALUE="${value}" ;;
      DATASETS_SHA256) DATASETS_SHA256_VALUE="${value}" ;;
      PIXI_TOML_SHA256) PIXI_TOML_SHA256_VALUE="${value}" ;;
      PIXI_LOCK_SHA256) PIXI_LOCK_SHA256_VALUE="${value}" ;;
      AUX_ROOT) AUX_ROOT_VALUE="${value}" ;;
      SCGATE_DB_BRANCH) SCGATE_DB_BRANCH_VALUE="${value}" ;;
    esac
  done < "${manifest}"
  [[ "${index}" -eq "${#keys[@]}" ]] || {
    error "source manifest has the wrong number of fields: ${manifest}"
    return 1
  }
  [[ "${SOURCE_FORMAT}" == "${SOURCE_MANIFEST_FORMAT}" ]] || {
    error "unsupported source manifest FORMAT: ${SOURCE_FORMAT}"
    return 1
  }
  [[ "${SOURCE_COMMIT_VALUE}" =~ ^[0-9a-fA-F]{40}$ ]] || {
    error "source manifest SOURCE_COMMIT is not a full commit: ${manifest}"
    return 1
  }
  for value in "${SOURCE_ARCHIVE_SHA256_VALUE}" "${CONFIG_HELPER_SHA256_VALUE}" \
    "${DATASETS_SHA256_VALUE}" "${PIXI_TOML_SHA256_VALUE}" "${PIXI_LOCK_SHA256_VALUE}"; do
    [[ "${value}" =~ ^[0-9a-fA-F]{64}$ ]] || {
      error "source manifest contains an invalid SHA-256 digest: ${manifest}"
      return 1
    }
  done
  return 0
}

runtime_manifest_load() {
  local manifest="${1:-}"
  local line key value
  local seen_format=0 seen_image_path=0 seen_image_sha=0
  local seen_build_revision=0 seen_env=0 seen_layout=0 seen_prefix=0
  local seen_base=0 seen_pixitainer=0 seen_pixi=0 seen_apptainer=0
  local seen_image_pixi=0 seen_image_lock=0
  [[ -f "${manifest}" && ! -L "${manifest}" && -r "${manifest}" && -s "${manifest}" ]] || {
    error "runtime manifest is missing, unreadable, or empty: ${manifest}"
    return 1
  }
  while IFS= read -r line || [[ -n "${line}" ]]; do
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${line}" == *=* && "${key}" =~ ^[A-Z][A-Z0-9_]*$ && -n "${value}" ]] || {
      error "runtime manifest has invalid syntax: ${manifest}"
      return 1
    }
    is_safe_value "${value}" || {
      error "runtime manifest contains an unsafe value: ${manifest}"
      return 1
    }
    case "${key}" in
      FORMAT)
        ((seen_format++ == 0)) || { error "runtime manifest duplicates FORMAT"; return 1; }
        RUNTIME_FORMAT="${value}" ;;
      IMAGE_PATH)
        ((seen_image_path++ == 0)) || { error "runtime manifest duplicates IMAGE_PATH"; return 1; }
        RUNTIME_IMAGE_PATH="${value}" ;;
      IMAGE_SHA256)
        ((seen_image_sha++ == 0)) || { error "runtime manifest duplicates IMAGE_SHA256"; return 1; }
        RUNTIME_IMAGE_SHA256_VALUE="${value}" ;;
      IMAGE_BUILD_GIT_REVISION)
        ((seen_build_revision++ == 0)) || { error "runtime manifest duplicates IMAGE_BUILD_GIT_REVISION"; return 1; }
        RUNTIME_BUILD_REVISION="${value}" ;;
      RUNTIME_ENV)
        ((seen_env++ == 0)) || { error "runtime manifest duplicates RUNTIME_ENV"; return 1; }
        RUNTIME_ENV_VALUE="${value}" ;;
      RUNTIME_LAYOUT)
        ((seen_layout++ == 0)) || { error "runtime manifest duplicates RUNTIME_LAYOUT"; return 1; }
        RUNTIME_LAYOUT_VALUE="${value}" ;;
      CONTAINER_ENV_PREFIX)
        ((seen_prefix++ == 0)) || { error "runtime manifest duplicates CONTAINER_ENV_PREFIX"; return 1; }
        RUNTIME_PREFIX_VALUE="${value}" ;;
      BASE_IMAGE)
        ((seen_base++ == 0)) || { error "runtime manifest duplicates BASE_IMAGE"; return 1; }
        RUNTIME_BASE_VALUE="${value}" ;;
      PIXITAINER_VERSION)
        ((seen_pixitainer++ == 0)) || { error "runtime manifest duplicates PIXITAINER_VERSION"; return 1; }
        RUNTIME_PIXITAINER_VALUE="${value}" ;;
      PIXI_VERSION)
        ((seen_pixi++ == 0)) || { error "runtime manifest duplicates PIXI_VERSION"; return 1; }
        RUNTIME_PIXI_VALUE="${value}" ;;
      APPTAINER_VERSION)
        ((seen_apptainer++ == 0)) || { error "runtime manifest duplicates APPTAINER_VERSION"; return 1; }
        RUNTIME_APPTAINER_VALUE="${value}" ;;
      IMAGE_PIXI_TOML_SHA256)
        ((seen_image_pixi++ == 0)) || { error "runtime manifest duplicates IMAGE_PIXI_TOML_SHA256"; return 1; }
        RUNTIME_IMAGE_PIXI_VALUE="${value}" ;;
      IMAGE_PIXI_LOCK_SHA256)
        ((seen_image_lock++ == 0)) || { error "runtime manifest duplicates IMAGE_PIXI_LOCK_SHA256"; return 1; }
        RUNTIME_IMAGE_LOCK_VALUE="${value}" ;;
    esac
  done < "${manifest}"
  [[ "${seen_format}" -eq 1 && "${seen_image_path}" -eq 1 && "${seen_image_sha}" -eq 1 && \
    "${seen_build_revision}" -eq 1 && "${seen_env}" -eq 1 && "${seen_layout}" -eq 1 && \
    "${seen_prefix}" -eq 1 && "${seen_base}" -eq 1 && "${seen_pixitainer}" -eq 1 && \
    "${seen_pixi}" -eq 1 && "${seen_apptainer}" -eq 1 && "${seen_image_pixi}" -eq 1 && \
    "${seen_image_lock}" -eq 1 ]] || {
    error "runtime manifest is missing required FORMAT=2 identity fields: ${manifest}"
    return 1
  }
  [[ "${RUNTIME_FORMAT}" == "2" ]] || {
    error "snapshot execution requires runtime manifest FORMAT=2: ${manifest}"
    return 1
  }
  [[ "${RUNTIME_IMAGE_SHA256_VALUE}" =~ ^[0-9a-fA-F]{64}$ && \
    "${RUNTIME_IMAGE_PIXI_VALUE}" =~ ^[0-9a-fA-F]{64}$ && \
    "${RUNTIME_IMAGE_LOCK_VALUE}" =~ ^[0-9a-fA-F]{64}$ ]] || {
    error "runtime manifest has an invalid SHA-256 identity: ${manifest}"
    return 1
  }
  [[ "${RUNTIME_PREFIX_VALUE}" = /* ]] || {
    error "runtime manifest CONTAINER_ENV_PREFIX must be absolute"
    return 1
  }
  case "${RUNTIME_LAYOUT_VALUE}" in
    relocated|path-preserving) ;;
    *)
      error "runtime manifest RUNTIME_LAYOUT is unrecognized: ${RUNTIME_LAYOUT_VALUE}"
      return 1
      ;;
  esac
  return 0
}

check_commit_tree_membership() {
  local source_root="${1:-}"
  local commit="${2:-}"
  local rel type record meta path mode
  local required=(
    src config_helper.R datasets.json pixi.toml pixi.lock
    aux/scGateDB.rds aux/genes.blocklist.rds
    aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz
  )
  for rel in "${required[@]}"; do
    git -C "${source_root}" cat-file -e "${commit}:${rel}" 2>/dev/null || {
      error "requested commit is missing required path: ${rel}"
      return 1
    }
    type="$(git -C "${source_root}" cat-file -t "${commit}:${rel}" 2>/dev/null)" || {
      error "cannot inspect requested commit path: ${rel}"
      return 1
    }
    if [[ "${rel}" == "src" ]]; then
      [[ "${type}" == tree ]] || {
        error "requested commit src is not a tree"
        return 1
      }
    else
      [[ "${type}" == blob ]] || {
        error "requested commit required path is not a regular file: ${rel}"
        return 1
      }
    fi
  done

  while IFS= read -r -d '' record; do
    meta="${record%%$'\t'*}"
    path="${record#*$'\t'}"
    [[ "${path}" != *$'\n'* && "${path}" != *$'\r'* && "${path}" != *$'\t'* ]] || {
      error "requested commit contains an unsafe filename: ${path}"
      return 1
    }
    case "/${path}/" in
      */.git|*/.git/*)
        error "requested commit contains a .git path: ${path}"
        return 1
        ;;
    esac
    mode="${meta%% *}"
    case "${mode}" in
      120000|160000)
        error "requested commit contains an unsupported symlink or submodule: ${path}"
        return 1
        ;;
    esac
  done < <(git -C "${source_root}" ls-tree -r -z --full-tree "${commit}" 2>/dev/null)
  return 0
}

ensure_commit_files_extracted() {
  local source_root="${1:-}"
  local commit="${2:-}"
  local tree="${3:-}"
  local record meta path mode
  while IFS= read -r -d '' record; do
    meta="${record%%$'\t'*}"
    path="${record#*$'\t'}"
    mode="${meta%% *}"
    case "${mode}" in
      120000|160000) continue ;;
    esac
    [[ -f "${tree}/${path}" && ! -L "${tree}/${path}" ]] || {
      error "git archive omitted committed path: ${path}"
      return 1
    }
  done < <(git -C "${source_root}" ls-tree -r -z --full-tree "${commit}" 2>/dev/null)
  return 0
}

extract_scgate_branch() {
  local tree="${1:-}"
  local config="${tree}/src/slurm_config.sh"
  local branch count=0 line
  if [[ -f "${config}" && ! -L "${config}" ]]; then
    while IFS= read -r line || [[ -n "${line}" ]]; do
      case "${line}" in
        *SCGATE_DB_BRANCH=*)
          branch="${line#*SCGATE_DB_BRANCH=}"
          branch="${branch%%[[:space:]]*}"
          branch="${branch#\"}"
          branch="${branch%\"}"
          branch="${branch#\'}"
          branch="${branch%\'}"
          if [[ -n "${branch}" ]]; then
            count=$((count + 1))
            SCGATE_DB_BRANCH_VALUE="${branch}"
          fi
          ;;
      esac
    done < "${config}"
  fi
  if [[ "${count}" -eq 0 ]]; then
    SCGATE_DB_BRANCH_VALUE="${DEFAULT_SCGATE_DB_BRANCH}"
  elif [[ "${count}" -ne 1 ]]; then
    error "snapshot slurm_config.sh has duplicate SCGATE_DB_BRANCH values"
    return 1
  fi
  is_safe_value "${SCGATE_DB_BRANCH_VALUE}" || {
    error "snapshot SCGATE_DB_BRANCH is unsafe"
    return 1
  }
  return 0
}

write_source_manifest() {
  local manifest="${1:-}"
  local source_root="${2:-}"
  local commit="${3:-}"
  local archive="${4:-}"
  local archive_sha="${5:-}"
  local tree="${6:-}"
  local config_sha datasets_sha pixi_toml_sha pixi_lock_sha
  config_sha="$(sha256_file "${tree}/config_helper.R")" || return 1
  datasets_sha="$(sha256_file "${tree}/datasets.json")" || return 1
  pixi_toml_sha="$(sha256_file "${tree}/pixi.toml")" || return 1
  pixi_lock_sha="$(sha256_file "${tree}/pixi.lock")" || return 1
  extract_scgate_branch "${tree}" || return 1
  {
    printf '%s\n' \
      "FORMAT=${SOURCE_MANIFEST_FORMAT}" \
      "SOURCE_ROOT=${source_root}" \
      "SOURCE_COMMIT=${commit}" \
      "SOURCE_ARCHIVE_PATH=${archive}" \
      "SOURCE_ARCHIVE_SHA256=${archive_sha}" \
      "CONFIG_HELPER_SHA256=${config_sha}" \
      "DATASETS_SHA256=${datasets_sha}" \
      "PIXI_TOML_SHA256=${pixi_toml_sha}" \
      "PIXI_LOCK_SHA256=${pixi_lock_sha}" \
      "AUX_ROOT=${source_root}/aux" \
      "SCGATE_DB_BRANCH=${SCGATE_DB_BRANCH_VALUE}"
  } > "${manifest}" || {
    error "cannot write source manifest: ${manifest}"
    return 1
  }
  return 0
}

write_complete_marker() {
  local marker="${1:-}"
  printf 'COMPLETE\n' > "${marker}" || {
    error "cannot write snapshot completion marker: ${marker}"
    return 1
  }
}

validate_snapshot_paths_and_manifest() {
  local snapshot_root="${1:-}"
  local expected_commit="${2:-}"
  local source_root="${snapshot_root}/tree"
  local identity="${snapshot_root}/identity"
  local manifest="${identity}/source.manifest"
  local archive="${identity}/source.tar"
  local marker="${snapshot_root}/COMPLETE"
  local actual_archive_sha config_sha datasets_sha pixi_toml_sha pixi_lock_sha manifest_branch
  source_manifest_load "${manifest}" || return 1
  [[ "${SOURCE_ROOT_VALUE}" == "${source_root}" ]] || {
    error "source manifest SOURCE_ROOT does not match snapshot tree"
    return 1
  }
  [[ "${SOURCE_COMMIT_VALUE}" == "${expected_commit}" ]] || {
    error "source manifest SOURCE_COMMIT does not match snapshot directory"
    return 1
  }
  [[ "${SOURCE_ARCHIVE_PATH_VALUE}" == "${archive}" ]] || {
    error "source manifest SOURCE_ARCHIVE_PATH does not match snapshot identity"
    return 1
  }
  [[ "${AUX_ROOT_VALUE}" == "${source_root}/aux" ]] || {
    error "source manifest AUX_ROOT is not the sole snapshot auxiliary root"
    return 1
  }
  [[ -d "${source_root}" && -d "${identity}" && -f "${archive}" && \
    -f "${marker}" && ! -L "${source_root}" && ! -L "${identity}" && \
    ! -L "${archive}" && ! -L "${marker}" ]] || {
    error "snapshot is incomplete: ${snapshot_root}"
    return 1
  }
  [[ "$(cat "${marker}" 2>/dev/null)" == "COMPLETE" ]] || {
    error "snapshot completion marker is invalid: ${marker}"
    return 1
  }
  actual_archive_sha="$(sha256_file "${archive}")" || return 1
  [[ "${actual_archive_sha}" == "${SOURCE_ARCHIVE_SHA256_VALUE}" ]] || {
    error "source archive SHA-256 does not match source manifest"
    return 1
  }
  check_no_symlinks_or_git "${source_root}" || return 1
  extract_and_compare_archive "${archive}" "${source_root}" || return 1
  for rel in config_helper.R datasets.json pixi.toml pixi.lock; do
    [[ -f "${source_root}/${rel}" && ! -L "${source_root}/${rel}" ]] || {
      error "snapshot is missing required source file: ${rel}"
      return 1
    }
  done
  config_sha="$(sha256_file "${source_root}/config_helper.R")" || return 1
  datasets_sha="$(sha256_file "${source_root}/datasets.json")" || return 1
  pixi_toml_sha="$(sha256_file "${source_root}/pixi.toml")" || return 1
  pixi_lock_sha="$(sha256_file "${source_root}/pixi.lock")" || return 1
  [[ "${config_sha}" == "${CONFIG_HELPER_SHA256_VALUE}" && \
    "${datasets_sha}" == "${DATASETS_SHA256_VALUE}" && \
    "${pixi_toml_sha}" == "${PIXI_TOML_SHA256_VALUE}" && \
    "${pixi_lock_sha}" == "${PIXI_LOCK_SHA256_VALUE}" ]] || {
    error "source manifest dependency/file digests do not match snapshot tree"
    return 1
  }
  manifest_branch="${SCGATE_DB_BRANCH_VALUE}"
  extract_scgate_branch "${source_root}" || return 1
  [[ "${manifest_branch}" == "${SCGATE_DB_BRANCH_VALUE}" ]] || {
    error "source manifest SCGATE_DB_BRANCH does not match snapshot tree"
    return 1
  }
  check_readonly_tree "${snapshot_root}" || return 1
  return 0
}

create_snapshot() {
  local source_input="${1:-}"
  local parent_input="${2:-}"
  local commit_input="${3:-}"
  local source_root parent resolved_commit git_root final tmp archive tree identity
  local archive_sha
  require_absolute_argument "--source-root" "${source_input}" || return 1
  require_absolute_argument "--snapshot-parent" "${parent_input}" || return 1
  [[ "${commit_input}" =~ ^[0-9a-fA-F]{40}$ ]] || {
    error "--commit must be a full 40-hex commit"
    return 1
  }
  source_root="$(canonical_existing "${source_input}")" || return 1
  [[ -d "${source_root}" ]] || {
    error "source root is not a directory: ${source_root}"
    return 1
  }
  parent="$(canonical_existing "${parent_input}")" || return 1
  [[ "${parent}" != "/" && -d "${parent}" && -w "${parent}" ]] || {
    error "snapshot parent must be a writable non-root directory: ${parent}"
    return 1
  }
  no_symlink_components "${source_root}" || {
    error "source root contains a symlink: ${source_root}"
    return 1
  }
  no_symlink_components "${parent}" || {
    error "snapshot parent contains a symlink: ${parent}"
    return 1
  }
  git_root="$(git -C "${source_root}" rev-parse --show-toplevel 2>/dev/null)" || {
    error "source root is not an existing Git work tree: ${source_root}"
    return 1
  }
  git_root="$(canonical_existing "${git_root}")" || return 1
  [[ "${git_root}" == "${source_root}" ]] || {
    error "--source-root must be the Git repository root: ${source_root}"
    return 1
  }
  [[ -z "$(git -C "${source_root}" status --porcelain --untracked-files=all 2>/dev/null)" ]] || {
    error "source checkout is dirty or has untracked files"
    return 1
  }
  resolved_commit="$(git -C "${source_root}" rev-parse --verify "${commit_input}^{commit}" 2>/dev/null)" || {
    error "requested commit cannot be resolved: ${commit_input}"
    return 1
  }
  [[ "${resolved_commit}" =~ ^[0-9a-fA-F]{40}$ ]] || {
    error "Git returned an invalid full commit for: ${commit_input}"
    return 1
  }
  resolved_commit="$(printf '%s' "${resolved_commit}" | tr '[:upper:]' '[:lower:]')" || return 1
  final="${parent}/${resolved_commit}"
  if [[ -e "${final}" || -L "${final}" ]]; then
    [[ -d "${final}" && ! -L "${final}" ]] || {
      error "commit-keyed snapshot identity conflicts with an existing path: ${final}"
      return 1
    }
    validate_snapshot_paths_and_manifest "${final}" "${resolved_commit}" || {
      error "existing commit-keyed snapshot is invalid or conflicts: ${final}"
      return 1
    }
    printf 'Source snapshot already exists and is verified: %s\n' "${final}"
    return 0
  fi

  check_commit_tree_membership "${source_root}" "${resolved_commit}" || return 1
  tmp="$(mktemp -d "${parent}/.ecoda-source-${resolved_commit}.XXXXXXXX")" || {
    error "cannot create same-filesystem temporary snapshot parent under ${parent}"
    return 1
  }
  local cleanup_tmp=1
  cleanup_snapshot_tmp() {
    if [[ "${cleanup_tmp}" -eq 1 && -n "${tmp}" && -d "${tmp}" ]]; then
      chmod -R u+w "${tmp}" >/dev/null 2>&1 || true
      rm -rf "${tmp}"
    fi
  }
  trap cleanup_snapshot_tmp EXIT

  tree="${tmp}/tree"
  identity="${tmp}/identity"
  archive="${identity}/source.tar"
  mkdir -p "${tree}" "${identity}" || {
    error "cannot create temporary snapshot layout"
    return 1
  }
  git -C "${source_root}" archive --format=tar --output="${archive}" "${resolved_commit}" 2>/dev/null || {
    error "cannot archive requested commit: ${resolved_commit}"
    return 1
  }
  archive_sha="$(sha256_file "${archive}")" || return 1
  archive_entries_safe "${archive}" || return 1
  tar -xf "${archive}" -C "${tree}" 2>/dev/null || {
    error "cannot extract requested source archive"
    return 1
  }
  check_no_symlinks_or_git "${tree}" || return 1
  ensure_commit_files_extracted "${source_root}" "${resolved_commit}" "${tree}" || return 1
  extract_and_compare_archive "${archive}" "${tree}" || return 1
  write_source_manifest "${tmp}/identity/source.manifest" "${final}/tree" "${resolved_commit}" \
    "${final}/identity/source.tar" "${archive_sha}" "${tree}" || return 1
  write_complete_marker "${tmp}/COMPLETE" || return 1

  # The content/list comparison above intentionally happened while extraction
  # still had normal archive modes.  Only after that proof do we make the
  # complete temporary publication immutable.
  chmod -R a-w "${tmp}" || {
    error "cannot make temporary source snapshot read-only"
    return 1
  }
  if mv -T "${tmp}" "${final}" 2>/dev/null; then
    :
  else
    if [[ -e "${final}" || -L "${final}" ]]; then
      error "commit-keyed snapshot appeared during atomic publication: ${final}"
      return 1
    fi
    mv "${tmp}" "${final}" || {
      error "cannot atomically publish source snapshot: ${final}"
      return 1
    }
  fi
  cleanup_tmp=0
  tmp=""
  trap - EXIT
  validate_snapshot_paths_and_manifest "${final}" "${resolved_commit}" || {
    error "published source snapshot failed final verification: ${final}"
    return 1
  }
  printf 'Created verified source snapshot: %s\n' "${final}"
  return 0
}

prepare_logs_root() {
  local logs_input="${1:-}"
  local logs_parent logs_base logs_parent_real logs_real
  require_absolute_argument "--logs-root" "${logs_input}" || return 1
  no_symlink_components "${logs_input}" 2>/dev/null || {
    # A missing final component is fine; canonicalize its existing parent.
    :
  }
  logs_base="${logs_input##*/}"
  logs_parent="${logs_input%/*}"
  [[ -n "${logs_base}" && -n "${logs_parent}" ]] || {
    error "--logs-root must name a directory below an existing parent"
    return 1
  }
  logs_parent_real="$(canonical_existing "${logs_parent}")" || return 1
  [[ -d "${logs_parent_real}" ]] || {
    error "logs parent is not a directory: ${logs_parent_real}"
    return 1
  }
  logs_real="${logs_parent_real%/}/${logs_base}"
  if [[ -e "${logs_real}" || -L "${logs_real}" ]]; then
    [[ -d "${logs_real}" && ! -L "${logs_real}" ]] || {
      error "logs root is not a real directory: ${logs_input}"
      return 1
    }
  else
    [[ -w "${logs_parent_real}" ]] || {
      error "logs parent is not writable: ${logs_parent_real}"
      return 1
    }
  fi
  PREPARED_LOGS_ROOT="${logs_real}"
  return 0
}

validate_runtime_identity_paths() {
  local image_input="${1:-}"
  local manifest_input="${2:-}"
  local image manifest version_dir runtime_collection runtime_id image_base
  require_absolute_argument "--runtime-image" "${image_input}" || return 1
  require_absolute_argument "--runtime-manifest" "${manifest_input}" || return 1
  image="$(canonical_existing "${image_input}")" || return 1
  manifest="$(canonical_existing "${manifest_input}")" || return 1
  [[ -f "${image}" && -r "${image}" && -s "${image}" && ! -L "${image}" ]] || {
    error "runtime image is missing, empty, unreadable, or a symlink: ${image_input}"
    return 1
  }
  [[ -f "${manifest}" && -r "${manifest}" && -s "${manifest}" && ! -L "${manifest}" ]] || {
    error "runtime manifest is missing, empty, unreadable, or a symlink: ${manifest_input}"
    return 1
  }
  [[ "${manifest}" == "${image}.manifest" ]] || {
    error "runtime manifest must be adjacent to the exact runtime image: ${manifest_input}"
    return 1
  }
  image_base="${image##*/}"
  [[ "${image_base}" == *.sif ]] || {
    error "runtime image must be a versioned .sif file: ${image}"
    return 1
  }
  version_dir="${image%/*}"
  runtime_id="${version_dir##*/}"
  runtime_collection="${version_dir%/*}"
  [[ "${runtime_collection##*/}" == "_ecoda_runtime" && \
    "${runtime_id}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || {
    error "runtime image/manifest must be under _ecoda_runtime/RUNTIME_ID: ${image}"
    return 1
  }
  no_symlink_components "${image}" || {
    error "runtime image path contains a symlink: ${image}"
    return 1
  }
  no_symlink_components "${manifest}" || {
    error "runtime manifest path contains a symlink: ${manifest}"
    return 1
  }
  RUNTIME_IMAGE_CANONICAL="${image}"
  RUNTIME_MANIFEST_CANONICAL="${manifest}"
  return 0
}

exec_snapshot() {
  local source_input="${1:-}" manifest_input="${2:-}" host_input="${3:-}"
  local image_input="${4:-}" runtime_manifest_input="${5:-}" run_id="${6:-}"
  local scratch_input="${7:-}" logs_input="${8:-}" script="${9:-}"
  local source_root source_manifest snapshot_root snapshot_parent snapshot_id archive marker aux
  local snapshot_parent_mode snapshot_parent_id snapshot_identity_id source_tree_id
  local identity_dir_id source_manifest_id archive_id marker_id script_id
  local scratch_root host_prefix script_path actual_archive_sha manifest_branch
  require_absolute_argument "--source-root" "${source_input}" || return 1
  require_absolute_argument "--source-manifest" "${manifest_input}" || return 1
  require_absolute_argument "--host-env-prefix" "${host_input}" || return 1
  require_absolute_argument "--runtime-image" "${image_input}" || return 1
  require_absolute_argument "--runtime-manifest" "${runtime_manifest_input}" || return 1
  require_absolute_argument "--scratch-root" "${scratch_input}" || return 1
  require_absolute_argument "--logs-root" "${logs_input}" || return 1
  [[ "${run_id}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || {
    error "--run-id must be a nonempty safe identifier"
    return 1
  }
  [[ -n "${script}" && "${script}" != /* && "${script}" != *$'\n'* && \
    "${script}" != *$'\r'* && "${script}" != *$'\t'* ]] || {
    error "--script must be a nonempty relative path"
    return 1
  }
  case "${script}" in
    ../*|*/../*|*/..|..|*//*)
      error "--script contains an unsafe relative path: ${script}"
      return 1
      ;;
  esac

  source_root="$(canonical_existing "${source_input}")" || return 1
  [[ -d "${source_root}" && ! -L "${source_root}" ]] || {
    error "snapshot source root is not a directory: ${source_input}"
    return 1
  }
  snapshot_root="${source_root%/}"
  [[ "${source_root##*/}" == tree ]] || {
    error "--source-root must be the tree directory of a commit snapshot"
    return 1
  }
  snapshot_root="${source_root%/tree}"
  snapshot_id="${snapshot_root##*/}"
  [[ "${snapshot_id}" =~ ^[0-9a-fA-F]{40}$ ]] || {
    error "snapshot tree parent is not a full commit identity: ${snapshot_root}"
    return 1
  }
  snapshot_root="$(canonical_existing "${snapshot_root}")" || return 1
  source_root="${snapshot_root}/tree"
  source_manifest="$(canonical_existing "${manifest_input}")" || return 1
  [[ "${source_manifest}" == "${snapshot_root}/identity/source.manifest" ]] || {
    error "source manifest is not the identity manifest for this snapshot"
    return 1
  }
  source_manifest_load "${source_manifest}" || return 1
  [[ "${SOURCE_ROOT_VALUE}" == "${source_root}" && \
    "${SOURCE_COMMIT_VALUE}" == "${snapshot_id}" && \
    "${SOURCE_ARCHIVE_PATH_VALUE}" == "${snapshot_root}/identity/source.tar" && \
    "${AUX_ROOT_VALUE}" == "${source_root}/aux" ]] || {
    error "source manifest does not identify this exact commit snapshot"
    return 1
  }
  [[ "${SOURCE_COMMIT_VALUE}" =~ ^[0-9a-fA-F]{40}$ ]] || {
    error "source manifest has an invalid SOURCE_COMMIT"
    return 1
  }
  archive="${snapshot_root}/identity/source.tar"
  marker="${snapshot_root}/COMPLETE"
  aux="${source_root}/aux"
  snapshot_parent="${snapshot_root%/*}"
  [[ "${snapshot_parent}" != "/" && -n "${snapshot_parent}" ]] || {
    error "snapshot commit identity has no usable parent: ${snapshot_root}"
    return 1
  }
  snapshot_parent_id="$(path_identity "${snapshot_parent}")" || return 1
  snapshot_identity_id="$(path_identity "${snapshot_root}")" || return 1
  source_tree_id="$(path_identity "${source_root}")" || return 1
  identity_dir_id="$(path_identity "${snapshot_root}/identity")" || return 1
  source_manifest_id="$(path_identity "${source_manifest}")" || return 1
  archive_id="$(path_identity "${archive}")" || return 1
  marker_id="$(path_identity "${marker}")" || return 1
  [[ -f "${marker}" && ! -L "${marker}" ]] || {
    error "snapshot is missing COMPLETE: ${marker}"
    return 1
  }
  [[ "$(cat "${marker}" 2>/dev/null)" == "COMPLETE" ]] || {
    error "snapshot COMPLETE marker is invalid"
    return 1
  }
  [[ -f "${archive}" && ! -L "${archive}" && -r "${archive}" ]] || {
    error "snapshot source archive is missing or unsafe: ${archive}"
    return 1
  }
  actual_archive_sha="$(sha256_file "${archive}")" || return 1
  [[ "${actual_archive_sha}" == "${SOURCE_ARCHIVE_SHA256_VALUE}" ]] || {
    error "snapshot source archive SHA-256 mismatch"
    return 1
  }
  check_no_symlinks_or_git "${snapshot_root}" || return 1
  extract_and_compare_archive "${archive}" "${source_root}" || return 1
  for rel in config_helper.R datasets.json pixi.toml pixi.lock; do
    [[ -f "${source_root}/${rel}" && ! -L "${source_root}/${rel}" ]] || {
      error "snapshot is missing required source file: ${rel}"
      return 1
    }
  done
  [[ -d "${aux}" && ! -L "${aux}" ]] || {
    error "snapshot auxiliary root is missing: ${aux}"
    return 1
  }
  for rel in aux/scGateDB.rds aux/genes.blocklist.rds aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
    [[ -f "${source_root}/${rel}" && ! -L "${source_root}/${rel}" ]] || {
      error "snapshot is missing required auxiliary file: ${rel}"
      return 1
    }
  done
  check_readonly_tree "${snapshot_root}" || return 1
  manifest_branch="${SCGATE_DB_BRANCH_VALUE}"
  extract_scgate_branch "${source_root}" || return 1
  [[ "${manifest_branch}" == "${SCGATE_DB_BRANCH_VALUE}" ]] || {
    error "source manifest SCGATE_DB_BRANCH does not match snapshot tree"
    return 1
  }

  host_prefix="$(canonical_existing "${host_input}")" || return 1
  [[ -d "${host_prefix}" && ! -L "${host_prefix}" ]] || {
    error "host environment prefix is not an existing directory: ${host_input}"
    return 1
  }
  scratch_root="$(canonical_existing "${scratch_input}")" || return 1
  [[ -d "${scratch_root}" && ! -L "${scratch_root}" && "${scratch_root}" != "/" ]] || {
    error "scratch root is not a safe existing directory: ${scratch_input}"
    return 1
  }
  path_is_within "${scratch_root}" "${source_root}" && {
    error "scratch root cannot be inside the immutable source tree"
    return 1
  }
  prepare_logs_root "${logs_input}" || return 1
  [[ "${PREPARED_LOGS_ROOT}" != "${source_root}" && \
    ! "${PREPARED_LOGS_ROOT}" == "${source_root}/"* && \
    "${PREPARED_LOGS_ROOT}" != "${snapshot_root}" && \
    ! "${PREPARED_LOGS_ROOT}" == "${snapshot_root}/"* ]] || {
    error "logs root cannot be inside the immutable source snapshot"
    return 1
  }
  validate_runtime_identity_paths "${image_input}" "${runtime_manifest_input}" || return 1
  runtime_manifest_load "${RUNTIME_MANIFEST_CANONICAL}" || return 1
  [[ "${RUNTIME_IMAGE_PATH}" == "${RUNTIME_IMAGE_CANONICAL}" ]] || {
    error "runtime manifest IMAGE_PATH does not match the supplied image"
    return 1
  }
  local runtime_parent image_mode manifest_mode parent_mode
  runtime_parent="${RUNTIME_IMAGE_CANONICAL%/*}"
  image_mode="$(file_mode "${RUNTIME_IMAGE_CANONICAL}")" || return 1
  manifest_mode="$(file_mode "${RUNTIME_MANIFEST_CANONICAL}")" || return 1
  parent_mode="$(file_mode "${runtime_parent}")" || return 1
  mode_is_readonly "${image_mode}" || {
    error "runtime image is writable: ${RUNTIME_IMAGE_CANONICAL}"
    return 1
  }
  mode_is_readonly "${manifest_mode}" || {
    error "runtime manifest is writable: ${RUNTIME_MANIFEST_CANONICAL}"
    return 1
  }
  mode_is_readonly "${parent_mode}" || {
    error "runtime image parent is writable: ${runtime_parent}"
    return 1
  }
  [[ "${RUNTIME_IMAGE_PIXI_VALUE}" == "${PIXI_TOML_SHA256_VALUE}" && \
    "${RUNTIME_IMAGE_LOCK_VALUE}" == "${PIXI_LOCK_SHA256_VALUE}" ]] || {
    error "runtime image dependency identity does not match the source snapshot"
    return 1
  }
  script_path="${source_root}/${script}"
  [[ -f "${script_path}" && -r "${script_path}" && ! -L "${script_path}" ]] || {
    error "script is not a readable regular file under the immutable source tree: ${script}"
    return 1
  }
  no_symlink_components "${script_path}" || {
    error "script path contains a symlink: ${script}"
    return 1
  }
  script_id="$(path_identity "${script_path}")" || return 1
  snapshot_parent_mode="$(file_mode "${snapshot_parent}")" || return 1

  # This is the only directory this executor may create.  Its parent must
  # already exist; no scratch run root, scheduler job, or worker directory is
  # created here.
  if [[ ! -e "${PREPARED_LOGS_ROOT}" ]]; then
    mkdir "${PREPARED_LOGS_ROOT}" || {
      error "cannot create explicit logs root: ${PREPARED_LOGS_ROOT}"
      return 1
    }
  fi
  [[ -d "${PREPARED_LOGS_ROOT}" && -w "${PREPARED_LOGS_ROOT}" ]] || {
    error "explicit logs root is not writable: ${PREPARED_LOGS_ROOT}"
    return 1
  }

  export ECODA_SOURCE_ROOT="${source_root}"
  export ECODA_SOURCE_MANIFEST="${source_manifest}"
  export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
  export ECODA_HOST_ENV_PREFIX="${host_prefix}"
  export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE_CANONICAL}"
  export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST_CANONICAL}"
  export ECODA_RUNTIME_MODE=apptainer
  export ECODA_RUN_ID="${run_id}"
  export HPC_SCRATCH_DIR="${scratch_root}"
  export ECODA_SCRATCH_ROOT="${scratch_root}"
  export ECODA_LOGS_DIR="${PREPARED_LOGS_ROOT}"
  export ECODA_AUX_ROOT="${aux}"
  export SCGATE_DB_BRANCH="${SCGATE_DB_BRANCH_VALUE}"
  export SCGATE_DB_PATH="${aux}/scGateDB.rds"
  export PYTHONDONTWRITEBYTECODE=1

  unset BASH_ENV ENV
  run_snapshot_with_parent_lock \
    "${snapshot_parent}" "${snapshot_parent_mode}" \
    "${source_root}" "${snapshot_root}" "${source_manifest}" "${script_path}" \
    "${snapshot_parent_id}" "${snapshot_identity_id}" "${source_tree_id}" \
    "${identity_dir_id}" "${source_manifest_id}" "${archive_id}" \
    "${marker_id}" "${script_id}" "${@:10}"
}

parse_create() {
  local source_root="" snapshot_parent="" commit="" arg
  shift
  while [[ "$#" -gt 0 ]]; do
    arg="$1"
    case "${arg}" in
      --source-root)
        [[ "$#" -ge 2 ]] || { error "--source-root requires a value"; return 2; }
        source_root="$2"; shift 2 ;;
      --snapshot-parent)
        [[ "$#" -ge 2 ]] || { error "--snapshot-parent requires a value"; return 2; }
        snapshot_parent="$2"; shift 2 ;;
      --commit)
        [[ "$#" -ge 2 ]] || { error "--commit requires a value"; return 2; }
        commit="$2"; shift 2 ;;
      -h|--help) usage ;;
      *) error "unknown create argument: ${arg}"; usage; return 2 ;;
    esac
  done
  [[ -n "${source_root}" && -n "${snapshot_parent}" && -n "${commit}" ]] || {
    error "create requires --source-root, --snapshot-parent, and --commit"
    usage
    return 2
  }
  create_snapshot "${source_root}" "${snapshot_parent}" "${commit}"
}
parse_exec() {
  local source_root="" source_manifest="" host_prefix="" runtime_image=""
  local runtime_manifest="" run_id="" scratch_root="" logs_root="" script=""
  local arg
  local separator=0
  local -a script_args=()
  shift
  while [[ "$#" -gt 0 ]]; do
    arg="$1"
    if [[ "${arg}" == "--" ]]; then
      separator=1
      shift
      script_args=("$@")
      break
    fi
    case "${arg}" in
      --source-root)
        [[ "$#" -ge 2 ]] || { error "--source-root requires a value"; return 2; }
        source_root="$2"; shift 2 ;;
      --source-manifest)
        [[ "$#" -ge 2 ]] || { error "--source-manifest requires a value"; return 2; }
        source_manifest="$2"; shift 2 ;;
      --host-env-prefix)
        [[ "$#" -ge 2 ]] || { error "--host-env-prefix requires a value"; return 2; }
        host_prefix="$2"; shift 2 ;;
      --runtime-image)
        [[ "$#" -ge 2 ]] || { error "--runtime-image requires a value"; return 2; }
        runtime_image="$2"; shift 2 ;;
      --runtime-manifest)
        [[ "$#" -ge 2 ]] || { error "--runtime-manifest requires a value"; return 2; }
        runtime_manifest="$2"; shift 2 ;;
      --run-id)
        [[ "$#" -ge 2 ]] || { error "--run-id requires a value"; return 2; }
        run_id="$2"; shift 2 ;;
      --scratch-root)
        [[ "$#" -ge 2 ]] || { error "--scratch-root requires a value"; return 2; }
        scratch_root="$2"; shift 2 ;;
      --logs-root)
        [[ "$#" -ge 2 ]] || { error "--logs-root requires a value"; return 2; }
        logs_root="$2"; shift 2 ;;
      --script)
        [[ "$#" -ge 2 ]] || { error "--script requires a value"; return 2; }
        script="$2"; shift 2 ;;
      -h|--help) usage ;;
      *) error "unknown exec argument (expected -- before script arguments): ${arg}"; usage; return 2 ;;
    esac
  done
  [[ "${separator}" -eq 1 ]] || {
    error "exec requires '--' before script arguments"
    usage
    return 2
  }
  [[ -n "${source_root}" && -n "${source_manifest}" && -n "${host_prefix}" && \
    -n "${runtime_image}" && -n "${runtime_manifest}" && -n "${run_id}" && \
    -n "${scratch_root}" && -n "${logs_root}" && -n "${script}" ]] || {
    error "exec requires every documented option and a relative --script"
    usage
    return 2
  }
  # Preserve arbitrary script arguments exactly.  exec_snapshot expects them
  # in positions ten onward, after its nine parsed values.
  exec_snapshot "${source_root}" "${source_manifest}" "${host_prefix}" \
    "${runtime_image}" "${runtime_manifest}" "${run_id}" "${scratch_root}" \
    "${logs_root}" "${script}" "${script_args[@]}"
}

main() {
  [[ "$#" -gt 0 ]] || { usage; return 2; }
  case "$1" in
    create) parse_create "$@" ;;
    exec) parse_exec "$@" ;;
    -h|--help) usage ;;
    *) error "unknown subcommand: $1"; usage; return 2 ;;
  esac
}

main "$@"
