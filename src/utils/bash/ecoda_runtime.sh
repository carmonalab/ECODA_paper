#!/bin/bash
# Shared immutable worker-runtime boundary for ECODA pipelines 2-5.
# Source after slurm_config.sh.  Submitters and host-side watchdogs stay on the
# host Pixi environment; only scientific workers cross into an Apptainer image.
# Bash 3.2-compatible: indexed arrays only, no namerefs or associative arrays.

ECODA_RUNTIME_BIND_ARGS=()
ECODA_RUNTIME_BIND_DESTS=""
ECODA_RUNTIME_LAYOUT=""
ECODA_RUNTIME_CONTAINER_PREFIX=""

_ecoda_runtime_die() {
  echo "ERROR: $*" >&2
  return 1
}

_ecoda_runtime_mode() {
  local mode="${1:-${ECODA_RUNTIME_MODE:-host}}"
  case "${mode}" in
    host|apptainer) printf '%s\n' "${mode}" ;;
    *) _ecoda_runtime_die "invalid ECODA_RUNTIME_MODE (expected host or apptainer): ${mode}"; return 1 ;;
  esac
}

_ecoda_runtime_profile() {
  case "${1:-}" in
    default|stage2|stage3|stage4|stage5) return 0 ;;
    *) _ecoda_runtime_die "invalid immutable runtime bind profile: ${1:-}"; return 1 ;;
  esac
}

_ecoda_runtime_sha256() {
  local path="${1:-}"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${path}" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${path}" | awk '{print $1}'
  else
    _ecoda_runtime_die "sha256sum or shasum is required for runtime validation"
    return 1
  fi
}

_ecoda_runtime_realpath_existing() {
  local path="${1:-}"
  [[ -n "${path}" ]] || {
    _ecoda_runtime_die "cannot resolve an empty runtime bind source"
    return 1
  }
  command -v realpath >/dev/null 2>&1 || {
    _ecoda_runtime_die "realpath is required for immutable runtime bind validation"
    return 1
  }
  local resolved
  if resolved="$(realpath -e "${path}" 2>/dev/null)"; then
    printf '%s\n' "${resolved}"
    return 0
  fi
  if [[ -e "${path}" || -L "${path}" ]] &&
     resolved="$(realpath "${path}" 2>/dev/null)"; then
    printf '%s\n' "${resolved}"
    return 0
  fi
  _ecoda_runtime_die "runtime bind source is missing or cannot be canonicalized: ${path}"
  return 1
}

_ecoda_runtime_manifest_value() {
  local manifest="${1:-}"
  local key="${2:-}"
  [[ "${key}" =~ ^[A-Z][A-Z0-9_]*$ ]] || return 1
  awk -v wanted="${key}" '
    index($0, wanted "=") == 1 {
      count++
      value = substr($0, length(wanted) + 2)
    }
    END {
      if (count != 1) exit 1
      print value
    }
  ' "${manifest}"
}

_ecoda_runtime_require_manifest_value() {
  local manifest="${1:-}"
  local key="${2:-}"
  local value
  if ! value="$(_ecoda_runtime_manifest_value "${manifest}" "${key}")"; then
    _ecoda_runtime_die "runtime manifest is missing or duplicates ${key}: ${manifest}"
    return 1
  fi
  [[ -n "${value}" ]] || {
    _ecoda_runtime_die "runtime manifest value is empty: ${key}"
    return 1
  }
  printf '%s\n' "${value}"
}

_ecoda_runtime_validate_manifest_shape() {
  local manifest="${1:-}"
  awk '
    BEGIN { ok = 1 }
    $0 !~ /^[A-Z][A-Z0-9_]*=[^[:space:]]+$/ { ok = 0 }
    {
      key = $0
      sub(/=.*/, "", key)
      if (++seen[key] > 1) ok = 0
    }
    END { exit(ok ? 0 : 1) }
  ' "${manifest}" || {
    _ecoda_runtime_die "runtime manifest has invalid syntax or duplicate keys: ${manifest}"
    return 1
  }
}
_ecoda_runtime_file_size() {
  local path="${1:-}"
  [[ -f "${path}" && -r "${path}" ]] || {
    _ecoda_runtime_die "runtime file is missing or unreadable: ${path}"
    return 1
  }
  wc -c < "${path}" | tr -d '[:space:]'
}

_ecoda_runtime_mode_bits() {
  local path="${1:-}"
  local mode=""
  mode="$(stat -c '%a' "${path}" 2>/dev/null || true)"
  if [[ ! "${mode}" =~ ^[0-7]+$ ]]; then
    mode="$(stat -f '%Lp' "${path}" 2>/dev/null || true)"
  fi
  [[ "${mode}" =~ ^[0-7]{3,4}$ ]] || {
    _ecoda_runtime_die "could not inspect runtime permissions: ${path}"
    return 1
  }
  printf '%s\n' "${mode: -3}"
}

_ecoda_runtime_require_nonwritable() {
  local path="${1:-}"
  local mode
  [[ -e "${path}" && ! -L "${path}" ]] || {
    _ecoda_runtime_die "runtime path is missing or is a symlink: ${path}"
    return 1
  }
  mode="$(_ecoda_runtime_mode_bits "${path}")" || return 1
  [[ "${mode}" != *[2367]* ]] || {
    _ecoda_runtime_die "runtime path is writable: ${path}"
    return 1
  }
}

_ecoda_runtime_require_tree_nonwritable() {
  local root="${1:-}"
  local write_bit
  _ecoda_runtime_require_nonwritable "${root}" || return 1
  for write_bit in 200 020 002; do
    [[ -z "$(find "${root}" -perm -"${write_bit}" -print -quit 2>/dev/null)" ]] || {
      _ecoda_runtime_die "immutable source tree contains writable paths: ${root}"
      return 1
    }
  done
}

_ecoda_runtime_validate_source_archive() {
  local source_root="${1:-}"
  local archive="${2:-}"
  local expected_sha="${3:-}"
  local actual_sha tmp entry
  [[ -d "${source_root}" && -r "${source_root}" ]] || {
    _ecoda_runtime_die "immutable source root is missing or unreadable: ${source_root}"
    return 1
  }
  [[ -f "${archive}" && ! -L "${archive}" && -r "${archive}" && -s "${archive}" ]] || {
    _ecoda_runtime_die "immutable source archive is missing or unsafe: ${archive}"
    return 1
  }
  [[ "${expected_sha}" =~ ^[[:xdigit:]]{64}$ ]] || {
    _ecoda_runtime_die "source archive SHA-256 is invalid: ${archive}"
    return 1
  }
  actual_sha="$(_ecoda_runtime_sha256 "${archive}")" || return 1
  [[ "${actual_sha}" == "${expected_sha}" ]] || {
    _ecoda_runtime_die "immutable source archive SHA-256 mismatch: ${archive}"
    return 1
  }
  tar -tf "${archive}" >/dev/null 2>&1 || {
    _ecoda_runtime_die "could not inspect immutable source archive: ${archive}"
    return 1
  }
  while IFS= read -r entry; do
    entry="${entry%/}"
    [[ -z "${entry}" || "${entry}" == "." ]] && continue
    case "${entry}" in
      /*|..|../*|*/../*) _ecoda_runtime_die "source archive contains an escaping path: ${entry}"; return 1 ;;
    esac
  done < <(tar -tf "${archive}" 2>/dev/null)
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-source-verify.XXXXXX")" || {
    _ecoda_runtime_die "could not create source archive verification directory"
    return 1
  }
  if ! tar -xf "${archive}" -C "${tmp}" 2>/dev/null; then
    rm -rf "${tmp}"
    _ecoda_runtime_die "could not extract immutable source archive: ${archive}"
    return 1
  fi
  if ! diff -qr "${tmp}" "${source_root}" >/dev/null 2>&1; then
    rm -rf "${tmp}"
    _ecoda_runtime_die "immutable source tree differs from its retained archive"
    return 1
  fi
  rm -rf "${tmp}"
}

_ecoda_runtime_validate_host_environment() {
  local prefix="${ECODA_HOST_ENV_PREFIX:-}"
  [[ "${prefix}" = /* && "${prefix}" == */.pixi/envs/py-cuda13 ]] || {
    _ecoda_runtime_die "ECODA_HOST_ENV_PREFIX must be an absolute .pixi/envs/py-cuda13 prefix"
    return 1
  }
  [[ -d "${prefix}" && -x "${prefix}/bin/python" && -x "${prefix}/bin/Rscript" ]] || {
    _ecoda_runtime_die "recorded host environment prefix is missing Python/Rscript: ${prefix}"
    return 1
  }
  if [[ -n "${ECODA_HOST_PYTHON_BIN:-}" ]]; then
    [[ "${ECODA_HOST_PYTHON_BIN}" == "${prefix}/bin/python" ]] || {
      _ecoda_runtime_die "ECODA_HOST_PYTHON_BIN does not use ECODA_HOST_ENV_PREFIX"
      return 1
    }
  fi
  if [[ -n "${ECODA_HOST_PIXI_RSCRIPT:-}" ]]; then
    [[ "${ECODA_HOST_PIXI_RSCRIPT}" == "${prefix}/bin/Rscript --vanilla" ]] || {
      _ecoda_runtime_die "ECODA_HOST_PIXI_RSCRIPT does not use ECODA_HOST_ENV_PREFIX"
      return 1
    }
  fi
}

_ecoda_runtime_host_binary_digests() {
  local require_recorded="${1:-0}"
  local expected_python="${ECODA_HOST_PYTHON_SHA256:-}"
  local expected_rscript="${ECODA_HOST_RSCRIPT_SHA256:-}"
  local actual_python actual_rscript
  case "${require_recorded}" in
    0|1) ;;
    *) _ecoda_runtime_die "host binary identity mode must be 0 or 1"; return 1 ;;
  esac
  _ecoda_runtime_validate_host_environment || return 1
  if [[ "${require_recorded}" == "1" ||
        -n "${expected_python}" || -n "${expected_rscript}" ]]; then
    [[ "${expected_python}" =~ ^[[:xdigit:]]{64}$ &&
       "${expected_rscript}" =~ ^[[:xdigit:]]{64}$ ]] || {
      _ecoda_runtime_die "run-bound host binary identity is missing or invalid"
      return 1
    }
  fi
  actual_python="$(_ecoda_runtime_sha256 "${ECODA_HOST_ENV_PREFIX}/bin/python")" || return 1
  actual_rscript="$(_ecoda_runtime_sha256 "${ECODA_HOST_ENV_PREFIX}/bin/Rscript")" || return 1
  if [[ -n "${expected_python}" ]]; then
    [[ "${actual_python}" == "${expected_python}" ]] || {
      _ecoda_runtime_die "host Python binary SHA-256 changed"
      return 1
    }
    [[ "${actual_rscript}" == "${expected_rscript}" ]] || {
      _ecoda_runtime_die "host Rscript binary SHA-256 changed"
      return 1
    }
  fi
  export ECODA_HOST_PYTHON_SHA256="${actual_python}"
  export ECODA_HOST_RSCRIPT_SHA256="${actual_rscript}"
  printf '%s\n%s\n' "${actual_python}" "${actual_rscript}"
}

_ecoda_runtime_validate_source_snapshot() {
  local runtime_manifest="${1:-}"
  local source_manifest source_root source_archive aux_root
  local manifest_root manifest_aux snapshot_root identity_dir complete
  local expected_archive expected_config expected_datasets expected_toml expected_lock
  local expected_commit expected_branch actual source_field_count
  local required_path required_sha
  [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]] || {
    _ecoda_runtime_die "format-2 runtime requires ECODA_SOURCE_SNAPSHOT_REQUIRED=1"
    return 1
  }
  source_manifest="${ECODA_SOURCE_MANIFEST:-}"
  source_root="${ECODA_SOURCE_ROOT:-}"
  [[ "${source_manifest}" = /* && -f "${source_manifest}" && ! -L "${source_manifest}" && -r "${source_manifest}" ]] || {
    _ecoda_runtime_die "ECODA_SOURCE_MANIFEST is missing or unsafe"
    return 1
  }
  [[ "${source_root}" = /* && -d "${source_root}" ]] || {
    _ecoda_runtime_die "ECODA_SOURCE_ROOT is missing or not a directory"
    return 1
  }
  source_field_count="$(wc -l < "${source_manifest}" | tr -d '[:space:]')" || return 1
  [[ "${source_field_count}" == "11" ]] || {
    _ecoda_runtime_die "source manifest must contain exactly its FORMAT=1 identity fields"
    return 1
  }
  [[ "$(_ecoda_runtime_require_manifest_value "${source_manifest}" FORMAT)" == "1" ]] || {
    _ecoda_runtime_die "source manifest FORMAT must be 1"
    return 1
  }
  expected_commit="$(_ecoda_runtime_require_manifest_value "${source_manifest}" SOURCE_COMMIT)" || return 1
  [[ "${expected_commit}" =~ ^[[:xdigit:]]{40}$ ]] || {
    _ecoda_runtime_die "source manifest SOURCE_COMMIT is not a full Git commit"
    return 1
  }
  manifest_root="$(_ecoda_runtime_require_manifest_value "${source_manifest}" SOURCE_ROOT)" || return 1
  manifest_aux="$(_ecoda_runtime_require_manifest_value "${source_manifest}" AUX_ROOT)" || return 1
  source_archive="$(_ecoda_runtime_require_manifest_value "${source_manifest}" SOURCE_ARCHIVE_PATH)" || return 1
  [[ "${source_archive}" = /* && -f "${source_archive}" && ! -L "${source_archive}" && -r "${source_archive}" ]] || {
    _ecoda_runtime_die "source archive is missing or unsafe"
    return 1
  }
  expected_archive="$(_ecoda_runtime_require_manifest_value "${source_manifest}" SOURCE_ARCHIVE_SHA256)" || return 1
  expected_config="$(_ecoda_runtime_require_manifest_value "${source_manifest}" CONFIG_HELPER_SHA256)" || return 1
  expected_datasets="$(_ecoda_runtime_require_manifest_value "${source_manifest}" DATASETS_SHA256)" || return 1
  expected_toml="$(_ecoda_runtime_require_manifest_value "${source_manifest}" PIXI_TOML_SHA256)" || return 1
  expected_lock="$(_ecoda_runtime_require_manifest_value "${source_manifest}" PIXI_LOCK_SHA256)" || return 1
  expected_branch="$(_ecoda_runtime_require_manifest_value "${source_manifest}" SCGATE_DB_BRANCH)" || return 1
  [[ "$(_ecoda_runtime_realpath_existing "${source_root}")" == "$(_ecoda_runtime_realpath_existing "${manifest_root}")" ]] || {
    _ecoda_runtime_die "source manifest SOURCE_ROOT does not match ECODA_SOURCE_ROOT"
    return 1
  }
  [[ "${manifest_aux}" == "${source_root%/}/aux" ]] || {
    _ecoda_runtime_die "source manifest AUX_ROOT is not the frozen source aux directory"
    return 1
  }
  [[ "${ECODA_AUX_ROOT:-${source_root%/}/aux}" == "${source_root%/}/aux" ]] || {
    _ecoda_runtime_die "ECODA_AUX_ROOT does not match the frozen source aux directory"
    return 1
  }
  source_root="$(_ecoda_runtime_realpath_existing "${source_root}")" || return 1
  source_manifest="$(_ecoda_runtime_realpath_existing "${source_manifest}")" || return 1
  source_archive="$(_ecoda_runtime_realpath_existing "${source_archive}")" || return 1
  snapshot_root="$(dirname "${source_root}")"
  identity_dir="${snapshot_root}/identity"
  complete="${snapshot_root}/COMPLETE"
  [[ "${source_root}" == "${snapshot_root}/tree" ]] || {
    _ecoda_runtime_die "ECODA_SOURCE_ROOT is not a commit-keyed snapshot tree"
    return 1
  }
  [[ "${source_manifest}" == "${identity_dir}/source.manifest" &&
     "${source_archive}" == "${identity_dir}/source.tar" &&
     -f "${complete}" ]] || {
    _ecoda_runtime_die "source manifest/archive are outside the verified snapshot identity"
    return 1
  }
  [[ ! -e "${source_root}/.git" && -d "${source_root}/aux" ]] || {
    _ecoda_runtime_die "snapshot source tree has an invalid Git or aux layout"
    return 1
  }
  for required_path in \
    "${source_root}/config_helper.R" "${source_root}/datasets.json" \
    "${source_root}/pixi.toml" "${source_root}/pixi.lock" \
    "${source_root}/aux/scGateDB.rds" "${source_root}/aux/genes.blocklist.rds" \
    "${source_root}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"; do
    [[ -f "${required_path}" && -r "${required_path}" ]] || {
      _ecoda_runtime_die "snapshot source file is missing or unreadable: ${required_path}"
      return 1
    }
  done
  for required_sha in \
    "${expected_archive}" "${expected_config}" "${expected_datasets}" \
    "${expected_toml}" "${expected_lock}"; do
    [[ "${required_sha}" =~ ^[[:xdigit:]]{64}$ ]] || {
      _ecoda_runtime_die "source manifest contains an invalid SHA-256 value"
      return 1
    }
  done
  actual="$(_ecoda_runtime_sha256 "${source_root}/config_helper.R")" || return 1
  [[ "${actual}" == "${expected_config}" ]] || { _ecoda_runtime_die "snapshot config_helper.R digest mismatch"; return 1; }
  actual="$(_ecoda_runtime_sha256 "${source_root}/datasets.json")" || return 1
  [[ "${actual}" == "${expected_datasets}" ]] || { _ecoda_runtime_die "snapshot datasets.json digest mismatch"; return 1; }
  actual="$(_ecoda_runtime_sha256 "${source_root}/pixi.toml")" || return 1
  [[ "${actual}" == "${expected_toml}" ]] || { _ecoda_runtime_die "snapshot pixi.toml digest mismatch"; return 1; }
  actual="$(_ecoda_runtime_sha256 "${source_root}/pixi.lock")" || return 1
  [[ "${actual}" == "${expected_lock}" ]] || { _ecoda_runtime_die "snapshot pixi.lock digest mismatch"; return 1; }
  if [[ -n "${SCGATE_DB_BRANCH:-}" ]]; then
    [[ "${SCGATE_DB_BRANCH}" == "${expected_branch}" ]] || {
      _ecoda_runtime_die "snapshot SCGATE_DB_BRANCH mismatch"
      return 1
    }
  fi
  _ecoda_runtime_validate_source_archive "${source_root}" "${source_archive}" "${expected_archive}" || return 1
  _ecoda_runtime_require_tree_nonwritable "${source_root}" || return 1
  if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
    _ecoda_runtime_validate_host_environment || return 1
  fi
  printf '%s\n' "${expected_toml}" "${expected_lock}"
}
_ecoda_runtime_require_run_identity() {
  local run_root="${ECODA_RUN_ROOT:-}"
  local run_root_base manifests_path identity
  local run_root_real manifests_real identity_real expected_identity
  [[ "${run_root}" = /* && -d "${run_root}" ]] || {
    _ecoda_runtime_die "bound runtime validation requires an absolute existing ECODA_RUN_ROOT"
    return 1
  }
  run_root_base="${run_root%/}"
  [[ -n "${run_root_base}" ]] || run_root_base="/"
  manifests_path="${run_root_base}/manifests"
  identity="${manifests_path}/runtime.identity"
  [[ -d "${manifests_path}" && ! -L "${manifests_path}" ]] || {
    _ecoda_runtime_die "run manifests directory is missing or unsafe"
    return 1
  }
  [[ -f "${identity}" && ! -L "${identity}" && -r "${identity}" && -s "${identity}" ]] || {
    _ecoda_runtime_die "run-bound runtime.identity is missing or unsafe: ${identity}"
    return 1
  }
  run_root_real="$(_ecoda_runtime_realpath_existing "${run_root}")" || return 1
  manifests_real="$(_ecoda_runtime_realpath_existing "${manifests_path}")" || return 1
  identity_real="$(_ecoda_runtime_realpath_existing "${identity}")" || return 1
  expected_identity="${run_root_real%/}/manifests/runtime.identity"
  [[ "${manifests_real}" == "${run_root_real%/}/manifests" &&
     "${identity_real}" == "${expected_identity}" ]] || {
    _ecoda_runtime_die "run-bound runtime.identity is not contained by ECODA_RUN_ROOT"
    return 1
  }
  printf '%s\n' "${identity_real}"
}

_ecoda_runtime_require_run_source_manifest() {
  local run_root="${ECODA_RUN_ROOT:-}"
  local run_root_base manifests_path run_manifest env_manifest
  local run_root_real manifests_real run_real env_real source_root source_root_real
  local expected_run_manifest expected_source_manifest
  [[ "${run_root}" = /* && -d "${run_root}" ]] || {
    _ecoda_runtime_die "snapshot-backed run requires an absolute existing ECODA_RUN_ROOT"
    return 1
  }
  run_root_base="${run_root%/}"
  [[ -n "${run_root_base}" ]] || run_root_base="/"
  manifests_path="${run_root_base}/manifests"
  run_manifest="${manifests_path}/source.manifest"
  env_manifest="${ECODA_SOURCE_MANIFEST:-}"
  [[ -d "${manifests_path}" && ! -L "${manifests_path}" ]] || {
    _ecoda_runtime_die "run manifests directory is missing or unsafe"
    return 1
  }
  [[ -f "${run_manifest}" && ! -L "${run_manifest}" && -r "${run_manifest}" && -s "${run_manifest}" ]] || {
    _ecoda_runtime_die "run-bound source.manifest is missing or unsafe: ${run_manifest}"
    return 1
  }
  [[ "${env_manifest}" = /* && -f "${env_manifest}" && ! -L "${env_manifest}" && -r "${env_manifest}" ]] || {
    _ecoda_runtime_die "ECODA_SOURCE_MANIFEST is missing or unsafe"
    return 1
  }
  [[ "${ECODA_SOURCE_ROOT:-}" = /* && -d "${ECODA_SOURCE_ROOT}" ]] || {
    _ecoda_runtime_die "snapshot-backed run requires an absolute existing ECODA_SOURCE_ROOT"
    return 1
  }
  run_root_real="$(_ecoda_runtime_realpath_existing "${run_root}")" || return 1
  manifests_real="$(_ecoda_runtime_realpath_existing "${manifests_path}")" || return 1
  run_real="$(_ecoda_runtime_realpath_existing "${run_manifest}")" || return 1
  env_real="$(_ecoda_runtime_realpath_existing "${env_manifest}")" || return 1
  source_root="${ECODA_SOURCE_ROOT}"
  source_root_real="$(_ecoda_runtime_realpath_existing "${source_root}")" || return 1
  expected_run_manifest="${run_root_real%/}/manifests/source.manifest"
  [[ "${manifests_real}" == "${run_root_real%/}/manifests" &&
     "${run_real}" == "${expected_run_manifest}" ]] || {
    _ecoda_runtime_die "run-bound source.manifest is not contained by ECODA_RUN_ROOT"
    return 1
  }
  case "${source_root_real}" in
    */tree) ;;
    *) _ecoda_runtime_die "ECODA_SOURCE_ROOT is not a commit-keyed snapshot tree"; return 1 ;;
  esac
  expected_source_manifest="${source_root_real%/tree}/identity/source.manifest"
  [[ "${env_real}" == "${expected_source_manifest}" ]] || {
    _ecoda_runtime_die "ECODA_SOURCE_MANIFEST is not the exact snapshot identity manifest"
    return 1
  }
  cmp -s "${run_real}" "${env_real}" || {
    _ecoda_runtime_die "run-bound source.manifest differs from ECODA_SOURCE_MANIFEST"
    return 1
  }
}


_ecoda_runtime_validate_source_identity() {
  local manifest="${1:-}"
  local expected_lock="${2:-}"
  local expected_revision="${3:-}"
  local expected_toml="${4:-}"
  local format current_lock current_toml current_revision source_values
  format="$(_ecoda_runtime_require_manifest_value "${manifest}" FORMAT)" || return 1
  case "${format}" in
    1|2) ;;
    *) _ecoda_runtime_die "unsupported immutable runtime manifest FORMAT: ${format}"; return 1 ;;
  esac
  if [[ "${format}" == "1" && "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" != "1" ]]; then
    [[ -f "${PROJECT_ROOT}/pixi.lock" ]] || {
      _ecoda_runtime_die "pixi.lock is missing for immutable runtime identity validation"
      return 1
    }
    current_lock="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.lock")" || return 1
    [[ "${current_lock}" == "${expected_lock}" ]] || {
      _ecoda_runtime_die "PIXI_LOCK_SHA256 mismatch for immutable runtime image"
      return 1
    }
    command -v git >/dev/null 2>&1 || {
      _ecoda_runtime_die "git is required for legacy immutable runtime validation"
      return 1
    }
    current_revision="$(git -C "${PROJECT_ROOT}" rev-parse HEAD 2>/dev/null)" || {
      _ecoda_runtime_die "could not resolve source revision for legacy immutable runtime validation"
      return 1
    }
    [[ -n "${current_revision}" && "${current_revision}" == "${expected_revision}" ]] || {
      _ecoda_runtime_die "GIT_REVISION mismatch for immutable runtime image"
      return 1
    }
    return 0
  fi
  if [[ "${ECODA_RUNTIME_BUILD_VALIDATION:-0}" == "1" ]]; then
    [[ -f "${PROJECT_ROOT}/pixi.toml" && -f "${PROJECT_ROOT}/pixi.lock" ]] || {
      _ecoda_runtime_die "build checkout is missing Pixi dependency files"
      return 1
    }
    current_toml="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.toml")" || return 1
    current_lock="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.lock")" || return 1
    [[ -n "${expected_toml}" && "${current_toml}" == "${expected_toml}" ]] || {
      _ecoda_runtime_die "IMAGE_PIXI_TOML_SHA256 does not match the build checkout"
      return 1
    }
    [[ "${current_lock}" == "${expected_lock}" ]] || {
      _ecoda_runtime_die "IMAGE_PIXI_LOCK_SHA256 does not match the build checkout"
      return 1
    }
    return 0
  fi
  _ecoda_runtime_require_run_source_manifest || return 1
  source_values="$(_ecoda_runtime_validate_source_snapshot "${manifest}")" || return 1
  current_toml="$(printf '%s\n' "${source_values}" | sed -n '1p')"
  current_lock="$(printf '%s\n' "${source_values}" | sed -n '2p')"
  if [[ "${format}" == "2" ]]; then
    [[ -n "${expected_toml}" && "${current_toml}" == "${expected_toml}" ]] || {
      _ecoda_runtime_die "source PIXI_TOML_SHA256 does not match immutable runtime image"
      return 1
    }
  fi
  [[ "${current_lock}" == "${expected_lock}" ]] || {
    _ecoda_runtime_die "source PIXI_LOCK_SHA256 does not match immutable runtime image"
    return 1
  }
}
_ecoda_runtime_write_identity() {
  local image="${1:-}"
  local manifest="${2:-}"
  local image_sha="${3:-}"
  local format="${4:-}"
  local identity identity_tmp manifest_sha image_size manifest_size image_toml image_lock
  [[ -n "${ECODA_RUN_ROOT:-}" ]] || {
    [[ "${format}" == "2" ]] || return 0
    _ecoda_runtime_die "format-2 submission requires ECODA_RUN_ROOT for runtime.identity"
    return 1
  }
  identity="${ECODA_RUN_ROOT}/manifests/runtime.identity"
  [[ -d "${ECODA_RUN_ROOT}/manifests" && -w "${ECODA_RUN_ROOT}/manifests" ]] || {
    _ecoda_runtime_die "run manifests directory is not writable: ${ECODA_RUN_ROOT}/manifests"
    return 1
  }
  manifest_sha="$(_ecoda_runtime_sha256 "${manifest}")" || return 1
  image_size="$(_ecoda_runtime_file_size "${image}")" || return 1
  manifest_size="$(_ecoda_runtime_file_size "${manifest}")" || return 1
  identity_tmp="${identity}.tmp.$$"
  umask 077
  {
    printf '%s\n' \
      "RUNTIME_IMAGE=${image}" \
      "RUNTIME_MANIFEST=${manifest}" \
      "RUNTIME_IMAGE_SHA256=${image_sha}" \
      "RUNTIME_MANIFEST_SHA256=${manifest_sha}" \
      "RUNTIME_IMAGE_SIZE=${image_size}" \
      "RUNTIME_MANIFEST_SIZE=${manifest_size}"
    if [[ "${format}" == "2" ]]; then
      image_toml="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PIXI_TOML_SHA256)" || {
        rm -f "${identity_tmp}"
        return 1
      }
      image_lock="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PIXI_LOCK_SHA256)" || {

        rm -f "${identity_tmp}"
        return 1
      }
      printf '%s\n' \
        "IMAGE_PIXI_TOML_SHA256=${image_toml}" \
        "IMAGE_PIXI_LOCK_SHA256=${image_lock}"
    fi
  } > "${identity_tmp}" || {
    rm -f "${identity_tmp}"
    _ecoda_runtime_die "could not write runtime identity"
    return 1
  }
  chmod 600 "${identity_tmp}" || {
    rm -f "${identity_tmp}"
    return 1
  }
  mv -f "${identity_tmp}" "${identity}" || {
    rm -f "${identity_tmp}"
    _ecoda_runtime_die "could not install runtime identity: ${identity}"
    return 1
  }
}

_ecoda_runtime_require_identity_value() {
  local identity="${1:-}"
  local key="${2:-}"
  _ecoda_runtime_require_manifest_value "${identity}" "${key}"
}

_ecoda_runtime_add_bind() {
  local source="${1:-}"
  local destination="${2:-}"
  local mode="${3:-}"
  [[ -n "${source}" && -n "${destination}" ]] || {
    _ecoda_runtime_die "runtime bind requires source and destination"
    return 1
  }
  case "${mode}" in
    ro|rw) ;;
    *) _ecoda_runtime_die "runtime bind mode must be ro or rw: ${mode}"; return 1 ;;
  esac
  [[ "${source}" = /* && "${destination}" = /* ]] || {
    _ecoda_runtime_die "runtime bind source and destination must be absolute: ${source}:${destination}"
    return 1
  }
  case "${destination}" in
    "${ECODA_RUNTIME_CONTAINER_PREFIX}"|"${ECODA_RUNTIME_CONTAINER_PREFIX}"/*)
      _ecoda_runtime_die "runtime bind destination collides with embedded environment: ${destination}"
      return 1
      ;;
  esac
  case $'\n'"${ECODA_RUNTIME_BIND_DESTS}"$'\n' in
    *$'\n'"${destination}"$'\n'*)
      _ecoda_runtime_die "duplicate runtime bind destination: ${destination}"
      return 1
      ;;
  esac
  ECODA_RUNTIME_BIND_DESTS="${ECODA_RUNTIME_BIND_DESTS}${ECODA_RUNTIME_BIND_DESTS:+$'\n'}${destination}"
  ECODA_RUNTIME_BIND_ARGS+=("${source}:${destination}:${mode}")
}

ecoda_runtime_build_bind_args() {
  local profile="${1:-}"
  local project_source source_root source_aux
  local scratch_source logs_source tmp_source ref_source
  local scratch_dest logs_dest tmp_dest source_datasets source_config
  local source_manifest_path snapshot_identity snapshot_root run_source_manifest run_runtime_identity
  local source_runtime_format="${ECODA_RUNTIME_FORMAT:-1}"
  _ecoda_runtime_profile "${profile}" || return 1
  [[ -n "${ECODA_RUNTIME_LAYOUT:-}" ]] || {
    _ecoda_runtime_die "runtime layout is unavailable; validate the image manifest first"
    return 1
  }
  [[ -n "${ECODA_RUNTIME_CONTAINER_PREFIX:-}" ]] || {
    _ecoda_runtime_die "container environment prefix is unavailable; validate the image manifest first"
    return 1
  }

  ECODA_RUNTIME_BIND_ARGS=()
  ECODA_RUNTIME_BIND_DESTS=""
  scratch_dest="${HPC_SCRATCH_DIR:-}"
  logs_dest="${LOGS_DIR:-}"
  tmp_dest="${TMPDIR:-/tmp}"
  if [[ "${source_runtime_format}" == "2" ]]; then
    source_root="${ECODA_SOURCE_ROOT:-}"
    source_aux="${ECODA_AUX_ROOT:-${source_root%/}/aux}"
    scratch_dest="${ECODA_SCRATCH_ROOT:-${HPC_SCRATCH_DIR:-}}"
    logs_dest="${ECODA_LOGS_DIR:-${LOGS_DIR:-}}"
    [[ "${source_root}" = /* && "${source_aux}" == "${source_root%/}/aux" ]] || {
      _ecoda_runtime_die "format-2 binds require the immutable source and aux roots"
      return 1
    }
    source_manifest_path="${ECODA_SOURCE_MANIFEST:-}"
    [[ "${source_manifest_path}" = /* && -f "${source_manifest_path}" && ! -L "${source_manifest_path}" ]] || {
      _ecoda_runtime_die "format-2 binds require the canonical source manifest"
      return 1
    }
    source_manifest_path="$(_ecoda_runtime_realpath_existing "${source_manifest_path}")" || return 1
    snapshot_identity="$(dirname "${source_manifest_path}")"
    snapshot_root="$(dirname "${snapshot_identity}")"
    [[ "${snapshot_identity##*/}" == "identity" &&
       -f "${snapshot_identity}/source.tar" &&
       -f "${snapshot_root}/COMPLETE" ]] || {
      _ecoda_runtime_die "format-2 source manifest is not a verified snapshot identity"
      return 1
    }
    [[ -n "${ECODA_RUN_ROOT:-}" && "${ECODA_RUN_ROOT}" = /* &&
       -d "${ECODA_RUN_ROOT}/manifests" ]] || {
      _ecoda_runtime_die "format-2 binds require the run manifests directory"
      return 1
    }
    run_source_manifest="${ECODA_RUN_ROOT}/manifests/source.manifest"
    run_runtime_identity="${ECODA_RUN_ROOT}/manifests/runtime.identity"
    source_root="$(_ecoda_runtime_realpath_existing "${source_root}")" || return 1
    source_aux="$(_ecoda_runtime_realpath_existing "${source_aux}")" || return 1
  else
    project_source="$(_ecoda_runtime_realpath_existing "${PROJECT_ROOT}")" || return 1
  fi
  scratch_source="$(_ecoda_runtime_realpath_existing "${scratch_dest}")" || return 1
  logs_source="$(_ecoda_runtime_realpath_existing "${logs_dest}")" || return 1
  tmp_source="$(_ecoda_runtime_realpath_existing "${tmp_dest}")" || return 1

  if [[ "${source_runtime_format}" == "2" ]]; then
    _ecoda_runtime_add_bind "${scratch_source}" "${scratch_dest}" rw || return 1
    _ecoda_runtime_add_bind "${logs_source}" "${logs_dest}" rw || return 1
    _ecoda_runtime_add_bind "${tmp_source}" "${tmp_dest}" rw || return 1
    if [[ "${profile}" == "stage4" ]]; then
      ref_source="$(_ecoda_runtime_realpath_existing "${HOME_REF_DIR}")" || return 1
      _ecoda_runtime_add_bind "${ref_source}" "${HOME_REF_DIR}" ro || return 1
    fi
    _ecoda_runtime_add_bind "${snapshot_root}" "${snapshot_root}" ro || return 1
    _ecoda_runtime_add_bind "${snapshot_identity}" "${snapshot_identity}" ro || return 1
    _ecoda_runtime_add_bind "${run_source_manifest}" "${run_source_manifest}" ro || return 1
    _ecoda_runtime_add_bind "${run_runtime_identity}" "${run_runtime_identity}" ro || return 1
    _ecoda_runtime_add_bind "${source_root}" "${source_root}" ro || return 1
    _ecoda_runtime_add_bind "${source_aux}" "${source_root}/aux" ro || return 1
    if [[ -x /usr/bin/scontrol && -r /usr/lib64/slurm/libslurmfull.so &&
          -r /etc/slurm/slurm.conf && -d /etc/slurm/slurm.d ]]; then
      _ecoda_runtime_add_bind /usr/bin/scontrol /usr/bin/scontrol ro || return 1
      _ecoda_runtime_add_bind /usr/lib64/slurm/libslurmfull.so /usr/lib64/slurm/libslurmfull.so ro || return 1
      _ecoda_runtime_add_bind /etc/slurm/slurm.conf /etc/slurm/slurm.conf ro || return 1
      _ecoda_runtime_add_bind /etc/slurm/slurm.d /etc/slurm/slurm.d ro || return 1
    fi
    return 0
  fi

  case "${ECODA_RUNTIME_LAYOUT}" in
    relocated)
      _ecoda_runtime_add_bind "${project_source}" "${PROJECT_ROOT}" ro || return 1
      ;;
    path-preserving)
      source_datasets="$(_ecoda_runtime_realpath_existing "${PROJECT_ROOT}/datasets.json")" || return 1
      source_config="$(_ecoda_runtime_realpath_existing "${PROJECT_ROOT}/config_helper.R")" || return 1
      source_aux="$(_ecoda_runtime_realpath_existing "${PROJECT_ROOT}/aux")" || return 1
      _ecoda_runtime_add_bind "${project_source}/src" "${PROJECT_ROOT}/src" ro || return 1
      _ecoda_runtime_add_bind "${source_datasets}" "${PROJECT_ROOT}/datasets.json" ro || return 1
      _ecoda_runtime_add_bind "${source_config}" "${PROJECT_ROOT}/config_helper.R" ro || return 1
      _ecoda_runtime_add_bind "${source_aux}" "${PROJECT_ROOT}/aux" ro || return 1
      ;;
    *)
      _ecoda_runtime_die "unrecognized runtime layout: ${ECODA_RUNTIME_LAYOUT}"
      return 1
      ;;
  esac

  _ecoda_runtime_add_bind "${scratch_source}" "${scratch_dest}" rw || return 1
  _ecoda_runtime_add_bind "${logs_source}" "${logs_dest}" rw || return 1
  _ecoda_runtime_add_bind "${tmp_source}" "${tmp_dest}" rw || return 1

  if [[ "${profile}" == "stage4" ]]; then
    ref_source="$(_ecoda_runtime_realpath_existing "${HOME_REF_DIR}")" || return 1
    _ecoda_runtime_add_bind "${ref_source}" "${HOME_REF_DIR}" ro || return 1
  fi
}

ecoda_runtime_validate_submission() {
  local mode="${1:-${ECODA_RUNTIME_MODE:-host}}"
  local image manifest image_path image_sha actual_sha
  local format runtime_env layout prefix base pixitainer pixi_version apptainer_version
  local git_revision lock_sha image_toml image_lock image_build_revision
  local manifest_project source_root
  _ecoda_runtime_mode "${mode}" >/dev/null || return 1
  case "${mode}" in
    host)
      if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
        _ecoda_runtime_require_run_source_manifest || return 1
        _ecoda_runtime_validate_source_snapshot "" >/dev/null || return 1
      fi
      return 0
      ;;
    apptainer)
      image="${ECODA_RUNTIME_IMAGE:-}"
      manifest="${ECODA_RUNTIME_MANIFEST:-${image}.manifest}"
      [[ "${image}" = /* ]] || {
        _ecoda_runtime_die "ECODA_RUNTIME_IMAGE must be an absolute path"
        return 1
      }
      [[ -f "${image}" && -r "${image}" && -s "${image}" ]] || {
        _ecoda_runtime_die "immutable runtime image is missing, unreadable, or empty: ${image}"
        return 1
      }
      [[ -f "${manifest}" && -r "${manifest}" && -s "${manifest}" ]] || {
        _ecoda_runtime_die "immutable runtime manifest is missing, unreadable, or empty: ${manifest}"
        return 1
      }
      _ecoda_runtime_validate_manifest_shape "${manifest}" || return 1
      image_path="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PATH)" || return 1
      image_sha="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_SHA256)" || return 1
      format="$(_ecoda_runtime_require_manifest_value "${manifest}" FORMAT)" || return 1
      runtime_env="$(_ecoda_runtime_require_manifest_value "${manifest}" RUNTIME_ENV)" || return 1
      layout="$(_ecoda_runtime_require_manifest_value "${manifest}" RUNTIME_LAYOUT)" || return 1
      prefix="$(_ecoda_runtime_require_manifest_value "${manifest}" CONTAINER_ENV_PREFIX)" || return 1
      base="$(_ecoda_runtime_require_manifest_value "${manifest}" BASE_IMAGE)" || return 1
      pixitainer="$(_ecoda_runtime_require_manifest_value "${manifest}" PIXITAINER_VERSION)" || return 1
      pixi_version="$(_ecoda_runtime_require_manifest_value "${manifest}" PIXI_VERSION)" || return 1
      apptainer_version="$(_ecoda_runtime_require_manifest_value "${manifest}" APPTAINER_VERSION)" || return 1
      case "${format}" in
        1)
          git_revision="$(_ecoda_runtime_require_manifest_value "${manifest}" GIT_REVISION)" || return 1
          lock_sha="$(_ecoda_runtime_require_manifest_value "${manifest}" PIXI_LOCK_SHA256)" || return 1
          image_toml=""
          image_build_revision="${git_revision}"
          ;;
        2)
          image_build_revision="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_BUILD_GIT_REVISION)" || return 1
          image_toml="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PIXI_TOML_SHA256)" || return 1
          lock_sha="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PIXI_LOCK_SHA256)" || return 1
          [[ "${image_build_revision}" =~ ^[^,[:space:]]+$ ]] || {
            _ecoda_runtime_die "IMAGE_BUILD_GIT_REVISION is invalid"
            return 1
          }
          [[ "${image_toml}" =~ ^[[:xdigit:]]{64}$ &&
             "${lock_sha}" =~ ^[[:xdigit:]]{64}$ ]] || {
            _ecoda_runtime_die "format-2 dependency identity is not SHA-256"
            return 1
          }
          ;;
        *)
          _ecoda_runtime_die "unsupported immutable runtime manifest FORMAT: ${format}"
          return 1
          ;;
      esac
      [[ "${image_path}" == "${image}" ]] || {
        _ecoda_runtime_die "IMAGE_PATH does not match ECODA_RUNTIME_IMAGE"
        return 1
      }
      [[ "${image_sha}" =~ ^[[:xdigit:]]{64}$ ]] || {
        _ecoda_runtime_die "IMAGE_SHA256 is not a SHA-256 digest"
        return 1
      }
      [[ "${runtime_env}" == "py-cuda13" ]] || {
        _ecoda_runtime_die "runtime manifest RUNTIME_ENV is not py-cuda13"
        return 1
      }
      case "${layout}" in
        relocated)
          [[ "${prefix}" == "/opt/ecoda/py-cuda13" ]] || {
            _ecoda_runtime_die "relocated runtime has an unexpected container prefix: ${prefix}"
            return 1
          }
          ;;
        path-preserving)
          [[ "${prefix}" == "${PROJECT_ROOT}/.pixi/envs/py-cuda13" ]] || {
            _ecoda_runtime_die "path-preserving runtime prefix does not match PROJECT_ROOT"
            return 1
          }
          manifest_project="$(_ecoda_runtime_require_manifest_value "${manifest}" CONTAINER_PROJECT_ROOT)" || return 1
          if [[ "${format}" == "2" ]]; then
            source_root="${ECODA_SOURCE_ROOT:-}"
            [[ "${ECODA_RUNTIME_BUILD_VALIDATION:-0}" == "1" ]] && source_root="${PROJECT_ROOT}"
            [[ -n "${source_root}" && "${manifest_project}" == "${source_root}" ]] || {
              _ecoda_runtime_die "format-2 path-preserving runtime source root differs from its build root"
              return 1
            }
          else
            [[ "${manifest_project}" == "${PROJECT_ROOT}" ]] || {
              _ecoda_runtime_die "path-preserving runtime was built for a different project root"
              return 1
            }
          fi
          ;;
        *)
          _ecoda_runtime_die "runtime manifest RUNTIME_LAYOUT is unrecognized: ${layout}"
          return 1
          ;;
      esac
      [[ "${base}" == "rockylinux:9" ]] || {
        _ecoda_runtime_die "runtime manifest BASE_IMAGE is not rockylinux:9"
        return 1
      }
      [[ "${pixitainer}" == "0.8.3" ]] || {
        _ecoda_runtime_die "runtime manifest PIXITAINER_VERSION is not 0.8.3"
        return 1
      }
      [[ -n "${pixi_version}" && -n "${apptainer_version}" ]] || {
        _ecoda_runtime_die "runtime manifest is missing toolchain identity"
        return 1
      }
      _ecoda_runtime_validate_source_identity \
        "${manifest}" "${lock_sha}" "${image_build_revision}" "${image_toml}" || return 1
      actual_sha="$(_ecoda_runtime_sha256 "${image}")" || return 1
      [[ "${actual_sha}" == "${image_sha}" ]] || {
        _ecoda_runtime_die "immutable runtime image SHA-256 mismatch"
        return 1
      }
      if [[ "${format}" == "2" && "${ECODA_RUNTIME_BUILD_VALIDATION:-0}" != "1" ]]; then
        _ecoda_runtime_require_nonwritable "${image}" || return 1
        _ecoda_runtime_require_nonwritable "${manifest}" || return 1
        _ecoda_runtime_require_nonwritable "$(dirname "${image}")" || return 1
      fi
      command -v "${APPTAINER_BIN:-apptainer}" >/dev/null 2>&1 || {
        _ecoda_runtime_die "apptainer is required for immutable runtime validation"
        return 1
      }
      "${APPTAINER_BIN:-apptainer}" inspect "${image}" >/dev/null 2>&1 || {
        _ecoda_runtime_die "apptainer inspect failed for immutable runtime image"
        return 1
      }
      ECODA_RUNTIME_FORMAT="${format}"
      ECODA_RUNTIME_LAYOUT="${layout}"
      ECODA_RUNTIME_CONTAINER_PREFIX="${prefix}"
      _ecoda_runtime_profile "${ECODA_RUNTIME_PROFILE:-default}" || return 1
      if [[ "${ECODA_RUNTIME_BUILD_VALIDATION:-0}" != "1" ]]; then
        ecoda_runtime_build_bind_args "${ECODA_RUNTIME_PROFILE:-default}" || return 1
        _ecoda_runtime_write_identity "${image}" "${manifest}" "${image_sha}" "${format}" || return 1
      fi
      ;;
  esac
}

ecoda_runtime_validate_bound_run() {
  local identity image manifest recorded_image_sha recorded_manifest_sha
  local format image_path image_sha runtime_env layout prefix base pixitainer
  local pixi_version apptainer_version image_toml image_lock image_build_revision
  local identity_toml identity_lock identity_field_count
  local manifest_project source_root
  [[ $# -eq 0 ]] || {
    _ecoda_runtime_die "ecoda_runtime_validate_bound_run takes no arguments"
    return 1
  }
  [[ -n "${ECODA_RUN_ROOT:-}" && "${ECODA_RUN_ROOT}" = /* ]] || {
    _ecoda_runtime_die "bound runtime validation requires ECODA_RUN_ROOT"
    return 1
  }
  identity="$(_ecoda_runtime_require_run_identity)" || return 1
  _ecoda_runtime_validate_manifest_shape "${identity}" || return 1
  recorded_image_sha="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_IMAGE_SHA256)" || return 1
  recorded_manifest_sha="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_MANIFEST_SHA256)" || return 1
  recorded_image_size="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_IMAGE_SIZE)" || return 1
  recorded_manifest_size="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_MANIFEST_SIZE)" || return 1
  image="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_IMAGE)" || return 1
  manifest="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_MANIFEST)" || return 1
  [[ "${image}" = /* && "${manifest}" = /* ]] || {
    _ecoda_runtime_die "run-bound runtime paths must be absolute"
    return 1
  }
  [[ "${recorded_image_sha}" =~ ^[[:xdigit:]]{64}$ &&
     "${recorded_manifest_sha}" =~ ^[[:xdigit:]]{64}$ &&
     "${recorded_image_size}" =~ ^[0-9]+$ &&
     "${recorded_manifest_size}" =~ ^[0-9]+$ ]] || {
    _ecoda_runtime_die "run-bound runtime identity has invalid digest or size"
    return 1
  }
  [[ -f "${image}" && -r "${image}" && -s "${image}" ]] || {
    _ecoda_runtime_die "run-bound runtime image is missing or unreadable: ${image}"
    return 1
  }
  [[ -f "${manifest}" && -r "${manifest}" && -s "${manifest}" ]] || {
    _ecoda_runtime_die "run-bound runtime manifest is missing or unreadable: ${manifest}"
    return 1
  }
  image_size="$(_ecoda_runtime_file_size "${image}")" || return 1
  manifest_size="$(_ecoda_runtime_file_size "${manifest}")" || return 1
  [[ "${image_size}" == "${recorded_image_size}" ]] || {
    _ecoda_runtime_die "run-bound runtime image size changed"
    return 1
  }
  [[ "${manifest_size}" == "${recorded_manifest_size}" ]] || {
    _ecoda_runtime_die "run-bound runtime manifest size changed"
    return 1
  }
  [[ "$(_ecoda_runtime_sha256 "${manifest}")" == "${recorded_manifest_sha}" ]] || {
    _ecoda_runtime_die "run-bound runtime manifest SHA-256 changed"
    return 1
  }
  _ecoda_runtime_require_nonwritable "${image}" || return 1
  _ecoda_runtime_require_nonwritable "${manifest}" || return 1
  _ecoda_runtime_require_nonwritable "$(dirname "${image}")" || return 1
  _ecoda_runtime_validate_manifest_shape "${manifest}" || return 1
  image_path="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PATH)" || return 1
  image_sha="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_SHA256)" || return 1
  format="$(_ecoda_runtime_require_manifest_value "${manifest}" FORMAT)" || return 1
  identity_field_count="$(wc -l < "${identity}" | tr -d '[:space:]')" || return 1
  case "${format}" in
    1)
      [[ "${identity_field_count}" == "6" ]] || {
        _ecoda_runtime_die "FORMAT=1 runtime.identity has unexpected fields"
        return 1
      }
      ;;
    2)
      [[ "${identity_field_count}" == "8" ]] || {
        _ecoda_runtime_die "FORMAT=2 runtime.identity has unexpected fields"
        return 1
      }
      identity_toml="$(_ecoda_runtime_require_identity_value "${identity}" IMAGE_PIXI_TOML_SHA256)" || return 1
      identity_lock="$(_ecoda_runtime_require_identity_value "${identity}" IMAGE_PIXI_LOCK_SHA256)" || return 1
      ;;
    *) _ecoda_runtime_die "unsupported run-bound runtime FORMAT: ${format}"; return 1 ;;
  esac
  runtime_env="$(_ecoda_runtime_require_manifest_value "${manifest}" RUNTIME_ENV)" || return 1
  layout="$(_ecoda_runtime_require_manifest_value "${manifest}" RUNTIME_LAYOUT)" || return 1
  prefix="$(_ecoda_runtime_require_manifest_value "${manifest}" CONTAINER_ENV_PREFIX)" || return 1
  base="$(_ecoda_runtime_require_manifest_value "${manifest}" BASE_IMAGE)" || return 1
  pixitainer="$(_ecoda_runtime_require_manifest_value "${manifest}" PIXITAINER_VERSION)" || return 1
  pixi_version="$(_ecoda_runtime_require_manifest_value "${manifest}" PIXI_VERSION)" || return 1
  apptainer_version="$(_ecoda_runtime_require_manifest_value "${manifest}" APPTAINER_VERSION)" || return 1
  [[ "${format}" == "2" ]] || {
    export ECODA_RUNTIME_IMAGE="${image}" ECODA_RUNTIME_MANIFEST="${manifest}"
    ecoda_runtime_validate_submission apptainer
    return $?
  }
  _ecoda_runtime_require_run_source_manifest || return 1
  [[ "${image}" == "${image_path}" && "${image_sha}" == "${recorded_image_sha}" ]] || {
    _ecoda_runtime_die "run-bound runtime identity does not match its image manifest"
    return 1
  }
  [[ "${manifest}" == "${image}.manifest" ]] || {
    _ecoda_runtime_die "format-2 runtime manifest is not beside its versioned image"
    return 1
  }
  case "${image}" in
    */_ecoda_runtime/*/*.sif) ;;
    *) _ecoda_runtime_die "format-2 runtime image is not versioned under _ecoda_runtime"; return 1 ;;
  esac
  image_build_revision="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_BUILD_GIT_REVISION)" || return 1
  image_toml="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PIXI_TOML_SHA256)" || return 1
  image_lock="$(_ecoda_runtime_require_manifest_value "${manifest}" IMAGE_PIXI_LOCK_SHA256)" || return 1
  [[ "${identity_toml}" == "${image_toml}" && "${identity_lock}" == "${image_lock}" ]] || {
    _ecoda_runtime_die "run-bound dependency identity does not match its image manifest"
    return 1
  }
  [[ "${identity_toml}" =~ ^[[:xdigit:]]{64}$ && "${identity_lock}" =~ ^[[:xdigit:]]{64}$ &&
     "${image_toml}" =~ ^[[:xdigit:]]{64}$ && "${image_lock}" =~ ^[[:xdigit:]]{64}$ ]] || {
    _ecoda_runtime_die "format-2 dependency identity is not SHA-256"
    return 1
  }
  [[ "${runtime_env}" == "py-cuda13" && "${base}" == "rockylinux:9" &&
     "${pixitainer}" == "0.8.3" && -n "${pixi_version}" && -n "${apptainer_version}" ]] || {
    _ecoda_runtime_die "format-2 runtime manifest has invalid toolchain identity"
    return 1
  }
  case "${layout}" in
    relocated)
      [[ "${prefix}" == "/opt/ecoda/py-cuda13" ]] || {
        _ecoda_runtime_die "format-2 relocated runtime prefix is invalid"
        return 1
      }
      ;;
    path-preserving)
      manifest_project="$(_ecoda_runtime_require_manifest_value "${manifest}" CONTAINER_PROJECT_ROOT)" || return 1
      source_root="${ECODA_SOURCE_ROOT:-}"
      [[ -n "${source_root}" && "${manifest_project}" == "${source_root}" ]] || {
        _ecoda_runtime_die "format-2 path-preserving runtime source root differs from its build root"
        return 1
      }
      ;;
    *) _ecoda_runtime_die "format-2 runtime layout is invalid: ${layout}"; return 1 ;;
  esac
  export ECODA_RUNTIME_IMAGE="${image}" ECODA_RUNTIME_MANIFEST="${manifest}"
  export ECODA_RUNTIME_FORMAT=2 ECODA_RUNTIME_LAYOUT="${layout}" ECODA_RUNTIME_CONTAINER_PREFIX="${prefix}"
  _ecoda_runtime_validate_source_identity \
    "${manifest}" "${image_lock}" "${image_build_revision}" "${image_toml}" || return 1
  if [[ "${ECODA_RUNTIME_MODE:-host}" == "host" &&
        "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
    _ecoda_runtime_host_binary_digests 1 >/dev/null || return 1
  fi
  _ecoda_runtime_profile "${ECODA_RUNTIME_PROFILE:-default}" || return 1
  ecoda_runtime_build_bind_args "${ECODA_RUNTIME_PROFILE:-default}" || return 1
}

ecoda_runtime_export_csv() {
  local profile="${1:-}"
  local nv="${2:-}"
  local mode="${ECODA_RUNTIME_MODE:-host}"
  local image="${ECODA_RUNTIME_IMAGE:-}"
  local manifest="${ECODA_RUNTIME_MANIFEST:-}"
  local format=""
  local identity=""
  local image_sha manifest_sha image_size manifest_size image_toml image_lock
  local host_python_sha host_rscript_sha host_binary_values
  local csv_value
  local -a csv_values
  local image_ro=0 manifest_ro=0 parent_ro=0
  local source_root="${ECODA_SOURCE_ROOT:-}"
  local source_manifest="${ECODA_SOURCE_MANIFEST:-}"
  local source_required="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
  local run_id="${ECODA_RUN_ID:-}"
  _ecoda_runtime_profile "${profile}" || return 1
  case "${nv}" in
    0|1) ;;
    *) _ecoda_runtime_die "ECODA_APPTAINER_NV must be 0 or 1: ${nv}"; return 1 ;;
  esac
  _ecoda_runtime_mode "${mode}" >/dev/null || return 1
  if [[ -n "${manifest}" && -f "${manifest}" ]]; then
    format="$(_ecoda_runtime_manifest_value "${manifest}" FORMAT 2>/dev/null || true)"
  fi
  if [[ "${format}" == "2" ]]; then
    identity="$(_ecoda_runtime_require_run_identity)" || return 1
    _ecoda_runtime_require_run_source_manifest || return 1
    _ecoda_runtime_validate_manifest_shape "${identity}" || return 1
    image="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_IMAGE)" || return 1
    manifest="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_MANIFEST)" || return 1
    image_sha="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_IMAGE_SHA256)" || return 1
    manifest_sha="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_MANIFEST_SHA256)" || return 1
    image_size="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_IMAGE_SIZE)" || return 1
    manifest_size="$(_ecoda_runtime_require_identity_value "${identity}" RUNTIME_MANIFEST_SIZE)" || return 1
    image_toml="$(_ecoda_runtime_require_identity_value "${identity}" IMAGE_PIXI_TOML_SHA256)" || return 1
    image_lock="$(_ecoda_runtime_require_identity_value "${identity}" IMAGE_PIXI_LOCK_SHA256)" || return 1
    [[ -n "${source_root}" && -n "${source_manifest}" && -n "${run_id}" ]] || {
      _ecoda_runtime_die "format-2 runtime export requires source and run identity"
      return 1
    }
    if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
      host_python_sha="${ECODA_HOST_PYTHON_SHA256:-}"
      host_rscript_sha="${ECODA_HOST_RSCRIPT_SHA256:-}"
      [[ "${host_python_sha}" =~ ^[[:xdigit:]]{64}$ &&
         "${host_rscript_sha}" =~ ^[[:xdigit:]]{64}$ ]] || {
        _ecoda_runtime_die "container runtime export requires recorded host binary identity"
        return 1
      }
    else
      host_binary_values="$(_ecoda_runtime_host_binary_digests)" || return 1
      host_python_sha="$(printf '%s\n' "${host_binary_values}" | sed -n '1p')"
      host_rscript_sha="$(printf '%s\n' "${host_binary_values}" | sed -n '2p')"
    fi
    export ECODA_HOST_PYTHON_SHA256="${host_python_sha}"
    export ECODA_HOST_RSCRIPT_SHA256="${host_rscript_sha}"
    _ecoda_runtime_require_nonwritable "${image}" >/dev/null 2>&1 && image_ro=1 || true
    _ecoda_runtime_require_nonwritable "${manifest}" >/dev/null 2>&1 && manifest_ro=1 || true
    _ecoda_runtime_require_nonwritable "$(dirname "${image}")" >/dev/null 2>&1 && parent_ro=1 || true
    csv_values=(
      "${mode}" "${image}" "${manifest}" "${profile}" "${nv}"
      "${source_root}" "${source_manifest}" "${source_required}" "${run_id}"
      "${image_sha}" "${manifest_sha}" "${image_size}" "${manifest_size}"
      "${image_toml}" "${image_lock}" "${image_ro}" "${manifest_ro}" "${parent_ro}"
      "${host_python_sha}" "${host_rscript_sha}"
    )
    for csv_value in "${csv_values[@]}"; do
      case "${csv_value}" in
        *[,]*|*$'\n'*) _ecoda_runtime_die "runtime export values cannot contain commas or newlines"; return 1 ;;
      esac
    done
    printf 'ECODA_RUNTIME_MODE=%s,ECODA_RUNTIME_IMAGE=%s,ECODA_RUNTIME_MANIFEST=%s,ECODA_RUNTIME_PROFILE=%s,ECODA_APPTAINER_NV=%s,ECODA_SOURCE_ROOT=%s,ECODA_SOURCE_MANIFEST=%s,ECODA_SOURCE_SNAPSHOT_REQUIRED=%s,ECODA_RUN_ID=%s,ECODA_RUNTIME_IMAGE_SHA256=%s,ECODA_RUNTIME_MANIFEST_SHA256=%s,ECODA_RUNTIME_IMAGE_SIZE=%s,ECODA_RUNTIME_MANIFEST_SIZE=%s,ECODA_IMAGE_PIXI_TOML_SHA256=%s,ECODA_IMAGE_PIXI_LOCK_SHA256=%s,ECODA_RUNTIME_IMAGE_READONLY=%s,ECODA_RUNTIME_MANIFEST_READONLY=%s,ECODA_RUNTIME_PARENT_READONLY=%s,ECODA_HOST_PYTHON_SHA256=%s,ECODA_HOST_RSCRIPT_SHA256=%s\n' \
      "${mode}" "${image}" "${manifest}" "${profile}" "${nv}" \
      "${source_root}" "${source_manifest}" "${source_required}" "${run_id}" \
      "${image_sha}" "${manifest_sha}" "${image_size}" "${manifest_size}" \
      "${image_toml}" "${image_lock}" "${image_ro}" "${manifest_ro}" "${parent_ro}" \
      "${host_python_sha}" "${host_rscript_sha}"
    return 0
  fi
  case "${mode}${image}${manifest}${profile}" in
    *[,]*|*$'\n'*) _ecoda_runtime_die "runtime export values cannot contain commas or newlines"; return 1 ;;
  esac
  [[ -n "${image}" && -n "${manifest}" ]] || {
    _ecoda_runtime_die "runtime export requires image and manifest paths"
    return 1
  }
  printf 'ECODA_RUNTIME_MODE=%s,ECODA_RUNTIME_IMAGE=%s,ECODA_RUNTIME_MANIFEST=%s,ECODA_RUNTIME_PROFILE=%s,ECODA_APPTAINER_NV=%s\n' \
    "${mode}" "${image}" "${manifest}" "${profile}" "${nv}"
}

_ecoda_runtime_require_source_script() {
  local script="${1:-}"
  local source_root="${2:-}"
  local script_real root_real
  [[ "${script}" = /* && -f "${script}" && -r "${script}" ]] || {
    _ecoda_runtime_die "worker source script must be an absolute readable file: ${script}"
    return 1
  }
  root_real="$(_ecoda_runtime_realpath_existing "${source_root}")" || return 1
  script_real="$(_ecoda_runtime_realpath_existing "${script}")" || return 1
  case "${script_real}" in
    "${root_real}"/*) printf '%s\n' "${script_real}" ;;
    *) _ecoda_runtime_die "worker source script is outside the immutable source root: ${script}"; return 1 ;;
  esac
}

ecoda_runtime_reexec_worker() {
  local profile="${1:-}"
  local script="${2:-}"
  local mode="${ECODA_RUNTIME_MODE:-host}"
  local image manifest prefix nv apptainer_bin
  local source_root source_manifest source_aux scratch_dest logs_dest
  local runtime_manifest_probe runtime_format_probe identity_probe
  local container_project data_file scgate_path
  local env_declaration env_name env_value bind_arg
  local -a apptainer_args
  _ecoda_runtime_profile "${profile}" || return 1
  _ecoda_runtime_mode "${mode}" >/dev/null || return 1
  runtime_manifest_probe="${ECODA_RUNTIME_MANIFEST:-}"
  identity_probe="${ECODA_RUN_ROOT:-}/manifests/runtime.identity"
  if [[ -n "${ECODA_RUN_ROOT:-}" &&
        ( -e "${identity_probe}" || -L "${identity_probe}" ) ]]; then
    identity_probe="$(_ecoda_runtime_require_run_identity)" || return 1
    runtime_manifest_probe="$(_ecoda_runtime_manifest_value "${identity_probe}" RUNTIME_MANIFEST 2>/dev/null || true)"
  fi
  runtime_format_probe=""
  if [[ -n "${runtime_manifest_probe}" && -f "${runtime_manifest_probe}" ]]; then
    runtime_format_probe="$(_ecoda_runtime_manifest_value "${runtime_manifest_probe}" FORMAT 2>/dev/null || true)"
  fi
  [[ "${runtime_format_probe}" != "2" || "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]] || {
    _ecoda_runtime_die "FORMAT=2 worker reexec requires ECODA_SOURCE_SNAPSHOT_REQUIRED=1"
    return 1
  }

  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    source_root="${ECODA_SOURCE_ROOT:-}"
    source_manifest="${ECODA_SOURCE_MANIFEST:-}"
    [[ "${source_root}" = /* && "${source_manifest}" = /* ]] || {
      _ecoda_runtime_die "snapshot-backed worker reexec requires source root and manifest"
      return 1
    }
    script="$(_ecoda_runtime_require_source_script "${script}" "${source_root}")" || return 1
  else
    [[ "${script}" = /* && -f "${script}" && -r "${script}" ]] || {
      _ecoda_runtime_die "worker source script must be an absolute readable file: ${script}"
      return 1
    }
  fi

  if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
    if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
      export ECODA_RUNTIME_PROFILE="${profile}"
      ecoda_runtime_validate_bound_run || return 1
    fi
    return 0
  fi
  if [[ "${mode}" == "host" ]]; then
    if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
      export ECODA_RUNTIME_PROFILE="${profile}"
      ecoda_runtime_validate_bound_run || return 1
    fi
    return 0
  fi

  export ECODA_RUNTIME_PROFILE="${profile}"
  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    ecoda_runtime_validate_bound_run || return 1
  else
    ecoda_runtime_validate_submission apptainer || return 1
  fi
  image="${ECODA_RUNTIME_IMAGE}"
  manifest="${ECODA_RUNTIME_MANIFEST}"
  prefix="${ECODA_RUNTIME_CONTAINER_PREFIX}"
  nv="${ECODA_APPTAINER_NV:-0}"
  apptainer_bin="${APPTAINER_BIN:-apptainer}"
  ecoda_runtime_build_bind_args "${profile}" || return 1
  container_project="${PROJECT_ROOT}"
  data_file="${DATASETS_JSON_FILE}"
  scratch_dest="${HPC_SCRATCH_DIR}"
  logs_dest="${LOGS_DIR}"
  scgate_path="${SCGATE_DB_PATH:-}"
  if [[ "${ECODA_RUNTIME_FORMAT:-1}" == "2" ]]; then
    source_root="${ECODA_SOURCE_ROOT}"
    source_aux="${ECODA_AUX_ROOT:-${source_root%/}/aux}"
    container_project="${source_root}"
    data_file="${source_root}/datasets.json"
    scratch_dest="${ECODA_SCRATCH_ROOT:-${HPC_SCRATCH_DIR}}"
    logs_dest="${ECODA_LOGS_DIR:-${LOGS_DIR}}"
    scgate_path="${source_aux}/scGateDB.rds"
  fi

  apptainer_args=(
    exec
    --containall
    --no-home
    --no-mount home,cwd,hostfs,bind-paths
  )
  while IFS= read -r env_declaration; do
    env_name="$(printf '%s\n' "${env_declaration}" |
      sed -n 's/^declare -x \([A-Za-z_][A-Za-z0-9_]*\).*/\1/p')"
    [[ -n "${env_name}" ]] || continue
    case "${env_name}" in
      PROJECT_ROOT|DATASETS_JSON_FILE|HPC_SCRATCH_DIR|LOGS_DIR|HOME_REF_DIR|\
      NAS_PREFIX|NAS_SC_DIR|NAS_TARGET_DIR|NAS_REF_DIR|SCGATE_DB_PATH|\
      SCGATE_DB_BRANCH|SCGATE_MODEL_CACHE_DIR|SCGATE_ONTOLOGY_BRANCH|\
      SAMPLE_COLNAME|PATH|LD_LIBRARY_PATH|TMPDIR|HOME|\
      PWD|OLDPWD|SHLVL|_|PYTHON_BIN|PIXI_RSCRIPT|R_HOME|RETICULATE_PYTHON|\
      PYTHONHOME|PYTHONPATH|PYTHONNOUSERSITE|R_LIBS_*|R_ENVIRON_USER|\
      R_PROFILE_USER|ECODA_HOST_*|ECODA_SOURCE_*|ECODA_AUX_ROOT|\
      ECODA_LOGS_DIR|ECODA_SCRATCH_ROOT|ECODA_RUNTIME_*|APPTAINER_*|\
      APPTAINERENV_*|SINGULARITY_*|SINGULARITYENV_*|BASH*|SHELLOPTS|\
      BASHOPTS|EUID|UID|PPID)
        continue
        ;;
    esac
    eval "env_value=\${${env_name}:-}"
    case "${env_value}" in
      *$'\n'*) _ecoda_runtime_die "exported runtime variable contains a newline: ${env_name}"; return 1 ;;
    esac
    apptainer_args+=(--env "${env_name}=${env_value}")
  done < <(export -p)
  if [[ "${nv}" == "1" ]]; then
    apptainer_args+=(--nv)
  fi
  for bind_arg in "${ECODA_RUNTIME_BIND_ARGS[@]}"; do
    apptainer_args+=(--bind "${bind_arg}")
  done
  apptainer_args+=(
    --env "PROJECT_ROOT=${container_project}"
    --env "DATASETS_JSON_FILE=${data_file}"
    --env "HPC_SCRATCH_DIR=${scratch_dest}"
    --env "LOGS_DIR=${logs_dest}"
    --env "HOME_REF_DIR=${HOME_REF_DIR:-}"
    --env "NAS_PREFIX=${NAS_PREFIX:-}"
    --env "NAS_SC_DIR=${NAS_SC_DIR:-}"
    --env "NAS_TARGET_DIR=${NAS_TARGET_DIR:-}"
    --env "NAS_REF_DIR=${NAS_REF_DIR:-}"
    --env "SCGATE_DB_PATH=${scgate_path}"
    --env "SCGATE_DB_BRANCH=${SCGATE_DB_BRANCH:-}"
    --env "SCGATE_MODEL_CACHE_DIR=${SCGATE_MODEL_CACHE_DIR:-}"
    --env "SCGATE_ONTOLOGY_BRANCH=${SCGATE_ONTOLOGY_BRANCH:-}"
    --env "SAMPLE_COLNAME=${SAMPLE_COLNAME:-}"
    --env "USER_EMAIL=${USER_EMAIL:-}"
    --env "TMPDIR=${TMPDIR:-/tmp}"
    --env "ECODA_RUNTIME_MODE=apptainer"
    --env "ECODA_RUNTIME_IMAGE=${image}"
    --env "ECODA_RUNTIME_MANIFEST=${manifest}"
    --env "ECODA_RUNTIME_PROFILE=${profile}"
    --env "ECODA_APPTAINER_NV=${nv}"
    --env "ECODA_RUNTIME_IN_CONTAINER=1"
    --env "ECODA_RUNTIME_PREFIX=${prefix}"
    --env "PYTHONDONTWRITEBYTECODE=1"
  )
  if [[ "${ECODA_RUNTIME_FORMAT:-1}" == "2" ]]; then
    apptainer_args+=(
      --env "ECODA_SOURCE_ROOT=${ECODA_SOURCE_ROOT}"
      --env "ECODA_SOURCE_MANIFEST=${ECODA_SOURCE_MANIFEST}"
      --env "ECODA_SOURCE_SNAPSHOT_REQUIRED=1"
      --env "ECODA_HOST_ENV_PREFIX=${ECODA_HOST_ENV_PREFIX}"
      --env "ECODA_HOST_PYTHON_SHA256=${ECODA_HOST_PYTHON_SHA256:-}"
      --env "ECODA_HOST_RSCRIPT_SHA256=${ECODA_HOST_RSCRIPT_SHA256:-}"
      --env "ECODA_AUX_ROOT=${ECODA_AUX_ROOT:-${ECODA_SOURCE_ROOT%/}/aux}"
      --env "ECODA_LOGS_DIR=${logs_dest}"
      --env "ECODA_SCRATCH_ROOT=${scratch_dest}"
      --env "ECODA_RUN_ID=${ECODA_RUN_ID:-}"
    )
  fi
  apptainer_args+=(
    "${image}"
    /bin/bash
    "${script}"
  )
  shift 2
  apptainer_args+=("$@")
  exec "${apptainer_bin}" "${apptainer_args[@]}"
}
