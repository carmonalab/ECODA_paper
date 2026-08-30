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

_ecoda_runtime_validate_source_identity() {
  local manifest="${1:-}"
  local expected_lock="${2:-}"
  local expected_revision="${3:-}"
  local current_lock current_revision
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
    _ecoda_runtime_die "git is required for immutable runtime source identity validation"
    return 1
  }
  current_revision="$(git -C "${PROJECT_ROOT}" rev-parse HEAD 2>/dev/null)" || {
    _ecoda_runtime_die "could not resolve source revision for immutable runtime validation"
    return 1
  }
  [[ -n "${current_revision}" && "${current_revision}" == "${expected_revision}" ]] || {
    _ecoda_runtime_die "GIT_REVISION mismatch for immutable runtime image"
    return 1
  }
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
  local project_source scratch_source logs_source tmp_source ref_source
  local source_datasets source_config source_aux
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
  project_source="$(_ecoda_runtime_realpath_existing "${PROJECT_ROOT}")" || return 1
  scratch_source="$(_ecoda_runtime_realpath_existing "${HPC_SCRATCH_DIR}")" || return 1
  logs_source="$(_ecoda_runtime_realpath_existing "${LOGS_DIR}")" || return 1
  tmp_source="$(_ecoda_runtime_realpath_existing "${TMPDIR:-/tmp}")" || return 1

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

  _ecoda_runtime_add_bind "${scratch_source}" "${HPC_SCRATCH_DIR}" rw || return 1
  _ecoda_runtime_add_bind "${logs_source}" "${LOGS_DIR}" rw || return 1
  _ecoda_runtime_add_bind "${tmp_source}" "${TMPDIR:-/tmp}" rw || return 1

  if [[ "${profile}" == "stage4" ]]; then
    ref_source="$(_ecoda_runtime_realpath_existing "${HOME_REF_DIR}")" || return 1
    _ecoda_runtime_add_bind "${ref_source}" "${HOME_REF_DIR}" ro || return 1
  fi
}

ecoda_runtime_validate_submission() {
  local mode="${1:-${ECODA_RUNTIME_MODE:-host}}"
  local image manifest image_path image_sha actual_sha
  local format runtime_env layout prefix base pixitainer pixi_version apptainer_version
  local git_revision lock_sha manifest_project
  _ecoda_runtime_mode "${mode}" >/dev/null || return 1
  case "${mode}" in
    host)
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
      git_revision="$(_ecoda_runtime_require_manifest_value "${manifest}" GIT_REVISION)" || return 1
      lock_sha="$(_ecoda_runtime_require_manifest_value "${manifest}" PIXI_LOCK_SHA256)" || return 1
      [[ "${image_path}" == "${image}" ]] || {
        _ecoda_runtime_die "IMAGE_PATH does not match ECODA_RUNTIME_IMAGE"
        return 1
      }
      [[ "${image_sha}" =~ ^[[:xdigit:]]{64}$ ]] || {
        _ecoda_runtime_die "IMAGE_SHA256 is not a SHA-256 digest"
        return 1
      }
      actual_sha="$(_ecoda_runtime_sha256 "${image}")" || return 1
      [[ "${actual_sha}" == "${image_sha}" ]] || {
        _ecoda_runtime_die "immutable runtime image SHA-256 mismatch"
        return 1
      }
      [[ "${format}" == "1" ]] || {
        _ecoda_runtime_die "unsupported immutable runtime manifest FORMAT: ${format}"
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
          [[ "${manifest_project}" == "${PROJECT_ROOT}" ]] || {
            _ecoda_runtime_die "path-preserving runtime was built for a different project root"
            return 1
          }
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
      _ecoda_runtime_validate_source_identity "${manifest}" "${lock_sha}" "${git_revision}" || return 1
      command -v "${APPTAINER_BIN:-apptainer}" >/dev/null 2>&1 || {
        _ecoda_runtime_die "apptainer is required for immutable runtime validation"
        return 1
      }
      "${APPTAINER_BIN:-apptainer}" inspect "${image}" >/dev/null 2>&1 || {
        _ecoda_runtime_die "apptainer inspect failed for immutable runtime image"
        return 1
      }
      ECODA_RUNTIME_LAYOUT="${layout}"
      ECODA_RUNTIME_CONTAINER_PREFIX="${prefix}"
      _ecoda_runtime_profile "${ECODA_RUNTIME_PROFILE:-default}" || return 1
      ecoda_runtime_build_bind_args "${ECODA_RUNTIME_PROFILE:-default}" || return 1
      ;;
  esac
}

ecoda_runtime_export_csv() {
  local profile="${1:-}"
  local nv="${2:-}"
  local mode="${ECODA_RUNTIME_MODE:-host}"
  local image="${ECODA_RUNTIME_IMAGE:-}"
  local manifest="${ECODA_RUNTIME_MANIFEST:-}"
  _ecoda_runtime_profile "${profile}" || return 1
  case "${nv}" in
    0|1) ;;
    *) _ecoda_runtime_die "ECODA_APPTAINER_NV must be 0 or 1: ${nv}"; return 1 ;;
  esac
  _ecoda_runtime_mode "${mode}" >/dev/null || return 1
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

ecoda_runtime_reexec_worker() {
  local profile="${1:-}"
  local script="${2:-}"
  local image manifest prefix nv
  local apptainer_bin
  local -a apptainer_args
  _ecoda_runtime_profile "${profile}" || return 1
  [[ "${script}" = /* && -f "${script}" && -r "${script}" ]] || {
    _ecoda_runtime_die "worker source script must be an absolute readable file: ${script}"
    return 1
  }
  _ecoda_runtime_mode "${ECODA_RUNTIME_MODE:-host}" >/dev/null || return 1
  if [[ "${ECODA_RUNTIME_MODE:-host}" == "host" || "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
    return 0
  fi

  export ECODA_RUNTIME_PROFILE="${profile}"
  ecoda_runtime_validate_submission apptainer || return 1
  image="${ECODA_RUNTIME_IMAGE}"
  manifest="${ECODA_RUNTIME_MANIFEST}"
  prefix="${ECODA_RUNTIME_CONTAINER_PREFIX}"
  nv="${ECODA_APPTAINER_NV:-0}"
  apptainer_bin="${APPTAINER_BIN:-apptainer}"
  ecoda_runtime_build_bind_args "${profile}" || return 1


  apptainer_args=(
    exec
    --containall
    --no-home
    --no-mount home,cwd,hostfs,bind-paths
  )
  local env_declaration env_name env_value
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
      R_PROFILE_USER|ECODA_HOST_*|ECODA_RUNTIME_*|APPTAINER_*|\
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
  local bind_arg
  for bind_arg in "${ECODA_RUNTIME_BIND_ARGS[@]}"; do
    apptainer_args+=(--bind "${bind_arg}")
  done
  apptainer_args+=(
    --env "PROJECT_ROOT=${PROJECT_ROOT}"
    --env "DATASETS_JSON_FILE=${DATASETS_JSON_FILE}"
    --env "HPC_SCRATCH_DIR=${HPC_SCRATCH_DIR}"
    --env "LOGS_DIR=${LOGS_DIR}"
    --env "HOME_REF_DIR=${HOME_REF_DIR:-}"
    --env "NAS_PREFIX=${NAS_PREFIX:-}"
    --env "NAS_SC_DIR=${NAS_SC_DIR:-}"
    --env "NAS_TARGET_DIR=${NAS_TARGET_DIR:-}"
    --env "NAS_REF_DIR=${NAS_REF_DIR:-}"
    --env "SCGATE_DB_PATH=${SCGATE_DB_PATH:-}"
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
    "${image}"
    /bin/bash
    "${script}"
  )
  shift 2
  apptainer_args+=("$@")
  exec "${apptainer_bin}" "${apptainer_args[@]}"
}
