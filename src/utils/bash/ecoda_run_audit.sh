#!/bin/bash
# Validator-only audit for one immutable ECODA run.
# This script deliberately performs no scheduler submission, owner repair, or
# filesystem discovery outside the paths named by the caller.
set -euo pipefail

PROGRAM="${0##*/}"
RUN_ROOT_ARG=""
STAGE_ARG=""
SELECTION_ARG=""
SOURCE_MANIFEST_ARG=""
RUNTIME_IDENTITY_ARG=""

AUDIT_METADATA_VARIANT=""
AUDIT_METADATA_ANALYSIS_ROOT=""
AUDIT_METADATA_ANALYSIS_NAS_ROOT=""
AUDIT_METADATA_ANALYSIS_PASS=""
AUDIT_METADATA_ANALYSIS_LOG_PREFIX=""
AUDIT_METADATA_ANALYSIS_ROOT_VERSION=""
AUDIT_METADATA_ANALYSIS_ROOT_IDENTITY=""
AUDIT_METADATA_METHOD_MATRIX=""
AUDIT_METADATA_METHOD_MATRIX_MD5=""
AUDIT_METADATA_METHOD_MATRIX_SIZE=""
AUDIT_METADATA_METHOD_MATRIX_SHA256=""
AUDIT_METADATA_METHOD_MATRIX_IDENTITY=""
AUDIT_METADATA_DECLARED_METHOD_ROWS=""
AUDIT_METADATA_PENDING_METHOD_ROWS=""
AUDIT_METADATA_METHOD_MATRIX_MODE=0
AUDIT_METADATA_EXPORT_MANIFEST=""
AUDIT_METADATA_EXPORT_STATUS=""
AUDIT_METADATA_METHODS=""
AUDIT_METADATA_PENDING_SELECTION=""
AUDIT_METADATA_PENDING_MD5=""
AUDIT_METADATA_PENDING_SIZE=""
AUDIT_METADATA_DISPATCH_SELECTION=""
AUDIT_METADATA_DISPATCH_MD5=""
AUDIT_METADATA_DISPATCH_SIZE=""
AUDIT_METADATA_DISPATCH_ROWS=""
AUDIT_METADATA_BATCH_CONTRACT_MANIFEST=""
AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_MD5=""
AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SIZE=""
AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SHA256=""

_audit_die() {
  echo "ERROR: $*" >&2
  return 1
}

_audit_usage() {
  cat <<'EOF'
Usage: ecoda_run_audit.sh --run-root ABSOLUTE_RUN_ROOT --stage STAGE \
       --selection ABSOLUTE_SELECTION \
       --source-manifest ABSOLUTE_SOURCE_MANIFEST \
       --runtime-identity ABSOLUTE_RUNTIME_IDENTITY
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-root)
      RUN_ROOT_ARG="${2:-}"
      shift 2
      ;;
    --run-root=*)
      RUN_ROOT_ARG="${1#*=}"
      shift
      ;;
    --stage)
      STAGE_ARG="${2:-}"
      shift 2
      ;;
    --stage=*)
      STAGE_ARG="${1#*=}"
      shift
      ;;
    --selection)
      SELECTION_ARG="${2:-}"
      shift 2
      ;;
    --selection=*)
      SELECTION_ARG="${1#*=}"
      shift
      ;;
    --source-manifest)
      SOURCE_MANIFEST_ARG="${2:-}"
      shift 2
      ;;
    --source-manifest=*)
      SOURCE_MANIFEST_ARG="${1#*=}"
      shift
      ;;
    --runtime-identity)
      RUNTIME_IDENTITY_ARG="${2:-}"
      shift 2
      ;;
    --runtime-identity=*)
      RUNTIME_IDENTITY_ARG="${1#*=}"
      shift
      ;;
    -h|--help)
      _audit_usage
      exit 0
      ;;
    *)
      _audit_usage >&2
      _audit_die "unknown argument: $1"
      exit 1
      ;;
  esac
done

[[ -n "${RUN_ROOT_ARG}" && -n "${STAGE_ARG}" &&
   -n "${SELECTION_ARG}" && -n "${SOURCE_MANIFEST_ARG}" &&
   -n "${RUNTIME_IDENTITY_ARG}" ]] || {
  _audit_usage >&2
  _audit_die "all five options are required"
  exit 1
}
[[ "${RUN_ROOT_ARG}" = /* && "${SELECTION_ARG}" = /* &&
   "${SOURCE_MANIFEST_ARG}" = /* && "${RUNTIME_IDENTITY_ARG}" = /* ]] || {
  _audit_die "run root, selection, source manifest, and runtime identity must be absolute"
  exit 1
}
[[ "${STAGE_ARG}" == stage2 || "${STAGE_ARG}" == stage3 ||
   "${STAGE_ARG}" == stage4 || "${STAGE_ARG}" == stage5 ]] || {
  _audit_die "unsupported stage: ${STAGE_ARG}"
  exit 1
}

command -v realpath >/dev/null 2>&1 || {
  _audit_die "realpath is required for run audit"
  exit 1
}
RUN_ROOT_REAL="$(realpath "${RUN_ROOT_ARG}" 2>/dev/null)" || {
  _audit_die "run root is missing or cannot be canonicalized: ${RUN_ROOT_ARG}"
  exit 1
}
[[ -d "${RUN_ROOT_REAL}" && ! -L "${RUN_ROOT_REAL}" ]] || {
  _audit_die "run root is not a regular directory: ${RUN_ROOT_ARG}"
  exit 1
}
case "${RUN_ROOT_REAL}" in
  */_ecoda_runs/*)
    SCRATCH_ROOT="${RUN_ROOT_REAL%/_ecoda_runs/*}"
    ;;
  *)
    _audit_die "run root is not under an exact _ecoda_runs root: ${RUN_ROOT_REAL}"
    exit 1
    ;;
esac
[[ -n "${SCRATCH_ROOT}" && "${SCRATCH_ROOT}" = /* &&
   -d "${SCRATCH_ROOT}" ]] || {
  _audit_die "run scratch root is missing: ${SCRATCH_ROOT}"
  exit 1
}
RUN_ID="${RUN_ROOT_REAL##*/}"
[[ "${RUN_ROOT_REAL}" == "${SCRATCH_ROOT}/_ecoda_runs/${RUN_ID}" ]] || {
  _audit_die "run root has an invalid canonical layout: ${RUN_ROOT_REAL}"
  exit 1
}

# Set these before sourcing the shared library so its roots are derived from
# this exact run, never from a caller's mutable checkout or current run root.
export HPC_SCRATCH_DIR="${SCRATCH_ROOT}"
export ECODA_RUNS_ROOT="${SCRATCH_ROOT}/_ecoda_runs"
export ECODA_OWNERS_ROOT="${SCRATCH_ROOT}/_ecoda_owners"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
ecoda_validate_run_id "${RUN_ID}" || exit 1

_audit_regular_file() {
  local path="$1"
  [[ -f "${path}" && ! -L "${path}" && -r "${path}" && -s "${path}" ]] || {
    _audit_die "required regular readable file is missing or empty: ${path}"
    return 1
  }
}

_audit_sha256_file() {
  local path="$1" digest
  if command -v sha256sum >/dev/null 2>&1; then
    digest="$(sha256sum "${path}" | awk '{print $1}')"
  elif command -v shasum >/dev/null 2>&1; then
    digest="$(shasum -a 256 "${path}" | awk '{print $1}')"
  else
    _audit_die "sha256sum or shasum is required for run audit"
    return 1
  fi
  [[ "${digest}" =~ ^[[:xdigit:]]{64}$ ]] || return 1
  printf '%s' "${digest}" | tr '[:upper:]' '[:lower:]'
}

_audit_same_bytes() {
  cmp -s "$1" "$2" || {
    _audit_die "run-owned manifest copy differs from supplied identity: $2"
    return 1
  }
}

_audit_source_manifest() {
  local manifest="$1"
  local line key value index=0 expected_key
  local keys=(FORMAT SOURCE_ROOT SOURCE_COMMIT SOURCE_ARCHIVE_PATH
    SOURCE_ARCHIVE_SHA256 CONFIG_HELPER_SHA256 DATASETS_SHA256
    PIXI_TOML_SHA256 PIXI_LOCK_SHA256 AUX_ROOT SCGATE_DB_BRANCH)
  local source_root source_archive aux_root snapshot_root identity_dir
  local source_file source_path expected_sha aux_file fresh entries entry details first
  _audit_regular_file "${manifest}" || return 1
  [[ "$(tail -c 1 "${manifest}" 2>/dev/null; printf '\001')" == $'\n\001' ]] || return 1
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    [[ ${index} -le ${#keys[@]} ]] || {
      _audit_die "source manifest has extra fields: ${manifest}"
      return 1
    }
    expected_key="${keys[$((index - 1))]}"
    [[ "${line}" == "${expected_key}="* ]] || {
      _audit_die "source manifest field ${index} must be ${expected_key}: ${manifest}"
      return 1
    }
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" && -n "${value}" &&
       "${value}" != *$'\n'* && "${value}" != *$'\t'* ]] || return 1
    case "${key}" in
      FORMAT) source_format="${value}" ;;
      SOURCE_ROOT) source_root="${value}" ;;
      SOURCE_COMMIT) source_commit="${value}" ;;
      SOURCE_ARCHIVE_PATH) source_archive="${value}" ;;
      SOURCE_ARCHIVE_SHA256) source_archive_sha="${value}" ;;
      CONFIG_HELPER_SHA256) source_config_sha="${value}" ;;
      DATASETS_SHA256) source_datasets_sha="${value}" ;;
      PIXI_TOML_SHA256) source_toml_sha="${value}" ;;
      PIXI_LOCK_SHA256) source_lock_sha="${value}" ;;
      AUX_ROOT) aux_root="${value}" ;;
      SCGATE_DB_BRANCH) source_branch="${value}" ;;
    esac
  done < "${manifest}"
  [[ ${index} -eq ${#keys[@]} && "${source_format}" == 1 &&
     "${source_commit}" =~ ^[[:xdigit:]]{40}$ ]] || {
    _audit_die "source manifest format or commit is invalid: ${manifest}"
    return 1
  }
  for value in "${source_archive_sha}" "${source_config_sha}" \
    "${source_datasets_sha}" "${source_toml_sha}" "${source_lock_sha}"; do
    [[ "${value}" =~ ^[[:xdigit:]]{64}$ ]] || {
      _audit_die "source manifest digest is invalid: ${manifest}"
      return 1
    }
  done
  [[ "${source_root}" = /* && "${aux_root}" = /* &&
     "${source_archive}" = /* && "${aux_root}" == "${source_root%/}/aux" ]] || {
    _audit_die "source manifest roots are invalid: ${manifest}"
    return 1
  }
  snapshot_root="${source_root%/tree}"
  [[ "${source_root##*/}" == tree &&
     "${snapshot_root##*/}" == "${source_commit}" ]] || {
    _audit_die "source manifest SOURCE_ROOT is not the commit snapshot tree: ${manifest}"
    return 1
  }
  [[ -d "${source_root}" && ! -L "${source_root}" ]] || return 1
  [[ -d "${aux_root}" && ! -L "${aux_root}" ]] || return 1
  [[ -f "${source_archive}" && ! -L "${source_archive}" &&
     -s "${source_archive}" ]] || return 1
  identity_dir="${snapshot_root}/identity"
  [[ "$(realpath "${source_root}")" == "${source_root}" &&
     "$(realpath "${aux_root}")" == "${aux_root}" &&
     "$(realpath "${source_archive}")" == "${source_archive}" &&
     "${manifest}" == "${identity_dir}/source.manifest" &&
     "${source_archive}" == "${identity_dir}/source.tar" &&
     -f "${snapshot_root}/COMPLETE" && ! -L "${snapshot_root}/COMPLETE" &&
     "$(cat "${snapshot_root}/COMPLETE" 2>/dev/null)" == COMPLETE &&
     ! -e "${source_root}/.git" ]] || {
    _audit_die "source manifest is not bound to one verified snapshot: ${manifest}"
    return 1
  }
  [[ "$(_audit_sha256_file "${source_archive}")" == "${source_archive_sha}" ]] || {
    _audit_die "source archive digest mismatch: ${source_archive}"
    return 1
  }
  for source_file in config_helper.R datasets.json pixi.toml pixi.lock; do
    source_path="${source_root}/${source_file}"
    _audit_regular_file "${source_path}" || return 1
    case "${source_file}" in
      config_helper.R) expected_sha="${source_config_sha}" ;;
      datasets.json) expected_sha="${source_datasets_sha}" ;;
      pixi.toml) expected_sha="${source_toml_sha}" ;;
      pixi.lock) expected_sha="${source_lock_sha}" ;;
    esac
    [[ "$(_audit_sha256_file "${source_path}")" == "${expected_sha}" ]] || {
      _audit_die "source snapshot digest mismatch: ${source_path}"
      return 1
    }
  done
  for aux_file in scGateDB.rds genes.blocklist.rds EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
    _audit_regular_file "${aux_root}/${aux_file}" || return 1
  done
  entries="$(tar -tf "${source_archive}" 2>/dev/null)" || {
    _audit_die "cannot list source archive: ${source_archive}"
    return 1
  }
  [[ -n "${entries}" ]] || return 1
  while IFS= read -r entry || [[ -n "${entry}" ]]; do
    entry="${entry%/}"
    [[ -n "${entry}" && "${entry}" != /* &&
       "${entry}" != ../* && "${entry}" != */../* &&
       "${entry}" != .. ]] || return 1
  done <<< "${entries}"
  details="$(tar -tvf "${source_archive}" 2>/dev/null)" || return 1
  while IFS= read -r entry || [[ -n "${entry}" ]]; do
    first="${entry#"${entry%%[![:space:]]*}"}"
    case "${first}" in
      l*) return 1 ;;
    esac
  done <<< "${details}"
  fresh="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-audit-source.XXXXXXXX")" || return 1
  if ! tar -xf "${source_archive}" -C "${fresh}" 2>/dev/null; then
    rm -rf "${fresh}"
    return 1
  fi
  if [[ -n "$(find "${fresh}" -type l -print -quit 2>/dev/null)" ||
        -n "$(find "${fresh}" -name .git -print -quit 2>/dev/null)" ]]; then
    rm -rf "${fresh}"
    return 1
  fi
  left_list="$(cd "${source_root}" && find . -print | LC_ALL=C sort)" || {
    rm -rf "${fresh}"
    return 1
  }
  right_list="$(cd "${fresh}" && find . -print | LC_ALL=C sort)" || {
    rm -rf "${fresh}"
    return 1
  }
  if [[ "${left_list}" != "${right_list}" ]] ||
     ! diff -r -q "${source_root}" "${fresh}" >/dev/null 2>&1; then
    rm -rf "${fresh}"
    _audit_die "source archive/tree contents differ: ${source_root}"
    return 1
  fi
  rm -rf "${fresh}"
  while IFS= read -r source_path; do
    [[ ! -w "${source_path}" ]] || {
      _audit_die "source snapshot is writable: ${source_path}"
      return 1
    }
  done < <(find "${source_root}" -print 2>/dev/null)
  AUDIT_SOURCE_ROOT="${source_root}"
  AUDIT_SOURCE_ARCHIVE="${source_archive}"
  AUDIT_SOURCE_SNAPSHOT="${snapshot_root}"
  AUDIT_SOURCE_TOML_SHA="${source_toml_sha}"
  AUDIT_SOURCE_LOCK_SHA="${source_lock_sha}"
}

_audit_runtime_identity() {
  local identity="$1"
  local line key value index=0 expected_key
  local keys=(RUNTIME_IMAGE RUNTIME_MANIFEST RUNTIME_IMAGE_SHA256
    RUNTIME_MANIFEST_SHA256 RUNTIME_IMAGE_SIZE RUNTIME_MANIFEST_SIZE)
  local image manifest image_sha manifest_sha image_size manifest_size
  local image_toml image_lock manifest_image_path manifest_image_sha
  _audit_regular_file "${identity}" || return 1
  [[ "$(tail -c 1 "${identity}" 2>/dev/null; printf '\001')" == $'\n\001' ]] || return 1
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    if [[ ${index} -le ${#keys[@]} ]]; then
      expected_key="${keys[$((index - 1))]}"
      [[ "${line}" == "${expected_key}="* ]] || {
        _audit_die "runtime identity field ${index} must be ${expected_key}: ${identity}"
        return 1
      }
    elif [[ ${index} -eq 7 ]]; then
      [[ "${line}" == IMAGE_PIXI_TOML_SHA256=* ]] || return 1
      expected_key=IMAGE_PIXI_TOML_SHA256
    elif [[ ${index} -eq 8 ]]; then
      [[ "${line}" == IMAGE_PIXI_LOCK_SHA256=* ]] || return 1
      expected_key=IMAGE_PIXI_LOCK_SHA256
    else
      _audit_die "runtime identity has extra fields: ${identity}"
      return 1
    fi
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" && -n "${value}" &&
       "${value}" != *$'\n'* && "${value}" != *$'\t'* ]] || return 1
    case "${key}" in
      RUNTIME_IMAGE) image="${value}" ;;
      RUNTIME_MANIFEST) manifest="${value}" ;;
      RUNTIME_IMAGE_SHA256) image_sha="${value}" ;;
      RUNTIME_MANIFEST_SHA256) manifest_sha="${value}" ;;
      RUNTIME_IMAGE_SIZE) image_size="${value}" ;;
      RUNTIME_MANIFEST_SIZE) manifest_size="${value}" ;;
      IMAGE_PIXI_TOML_SHA256) image_toml="${value}" ;;
      IMAGE_PIXI_LOCK_SHA256) image_lock="${value}" ;;
    esac
  done < "${identity}"
  [[ ${index} -eq 6 || ${index} -eq 8 ]] || return 1
  [[ "${image}" = /* && "${manifest}" = /* &&
     "${manifest}" == "${image}.manifest" &&
     "${image_sha}" =~ ^[[:xdigit:]]{64}$ &&
     "${manifest_sha}" =~ ^[[:xdigit:]]{64}$ &&
     "${image_size}" =~ ^[1-9][0-9]*$ &&
     "${manifest_size}" =~ ^[1-9][0-9]*$ ]] || {
    _audit_die "runtime identity values are malformed: ${identity}"
    return 1
  }
  if [[ ${index} -eq 8 ]]; then
    [[ "${image_toml}" =~ ^[[:xdigit:]]{64}$ &&
       "${image_lock}" =~ ^[[:xdigit:]]{64}$ &&
       "${image_toml}" == "${AUDIT_SOURCE_TOML_SHA}" &&
       "${image_lock}" == "${AUDIT_SOURCE_LOCK_SHA}" &&
       "${image}" == */_ecoda_runtime/*/ecoda-py-cuda13.sif ]] || return 1
  fi
  _audit_regular_file "${image}" || return 1
  _audit_regular_file "${manifest}" || return 1
  if [[ ${index} -eq 8 ]]; then
    [[ "$(realpath "${image}")" == "${image}" &&
       "$(realpath "${manifest}")" == "${manifest}" &&
       ! -w "${image}" && ! -w "${manifest}" &&
       ! -w "$(dirname "${image}")" &&
       ! -w "$(dirname "${manifest}")" ]] || {
      _audit_die "format-2 runtime image/manifest is writable: ${identity}"
      return 1
    }
  fi
  [[ "$(wc -c < "${image}" | tr -d '[:space:]')" == "${image_size}" &&
     "$(wc -c < "${manifest}" | tr -d '[:space:]')" == "${manifest_size}" ]] || {
    _audit_die "runtime identity file size mismatch: ${identity}"
    return 1
  }
  [[ "$(_audit_sha256_file "${manifest}")" == "${manifest_sha}" ]] || {
    _audit_die "runtime manifest digest mismatch: ${manifest}"
    return 1
  }
  manifest_image_path="$(sed -n 's/^IMAGE_PATH=//p' "${manifest}" | sed -n '1p')"
  manifest_image_sha="$(sed -n 's/^IMAGE_SHA256=//p' "${manifest}" | sed -n '1p')"
  [[ "${manifest_image_path}" == "${image}" &&
     "${manifest_image_sha}" == "${image_sha}" ]] || {
    _audit_die "runtime identity does not match runtime manifest: ${identity}"
    return 1
  }
  AUDIT_RUNTIME_IMAGE="${image}"
  AUDIT_RUNTIME_MANIFEST="${manifest}"
}

_audit_metadata_count() {
  local metadata="$1" field="$2"
  sed -n "s/^${field}=//p" "${metadata}" | wc -l | tr -d '[:space:]'
}

_audit_metadata_value() {
  local metadata="$1" field="$2"
  sed -n "s/^${field}=//p" "${metadata}" | sed -n '1p'
}

_audit_stage5_identity() {
  local variant="${AUDIT_METADATA_VARIANT:-}"
  local pass="${AUDIT_METADATA_PASS:-}"
  local analysis_root="${AUDIT_METADATA_ANALYSIS_ROOT:-}"
  local analysis_nas_root="${AUDIT_METADATA_ANALYSIS_NAS_ROOT:-}"
  local analysis_pass="${AUDIT_METADATA_ANALYSIS_PASS:-}"
  local expected_suffix expected_pass expected_root expected_nas
  local expected_log_prefix
  local analysis_log_prefix="${AUDIT_METADATA_ANALYSIS_LOG_PREFIX:-}"
  local scratch_root nas_target nas_base

  if [[ "${STAGE_ARG}" != stage5 && -n "${variant}" ]]; then
    _audit_die "ANALYSIS_VARIANT is only valid for Stage 5 run metadata"
    return 1
  fi
  [[ "${STAGE_ARG}" == stage5 ]] || return 0
  case "${pass}" in
    ""|uncorrected|corrected) ;;
    *)
      _audit_die "invalid Stage 5 PASS metadata: ${pass}"
      return 1
      ;;
  esac

  # Do not let caller-provided analysis or matrix state influence an audit.
  # Run metadata and the run-owned matrix copy are the only identity sources.
  unset ANALYSIS_VARIANT ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_PASS \
    ANALYSIS_LOG_PREFIX PASS_ARG ECODA_STAGE5_METHOD_MATRIX METHOD_MATRIX \
    ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION

  if [[ -n "${variant}" ]]; then
    case "${variant}" in
      final)
        expected_suffix="uncorrected_final"
        expected_pass="uncorrected"
        ;;
      corrected_final)
        expected_suffix="corrected_final/recovery_35row"
        expected_pass="corrected"
        expected_log_prefix="execution_times_batch_effect_corrected_final_"
        ;;
      *)
        _audit_die "unsupported Stage 5 analysis variant: ${variant}"
        return 1
        ;;
    esac
    scratch_root="${HPC_SCRATCH_DIR%/}"
    [[ -n "${scratch_root}" ]] || scratch_root="/"
    if [[ "${scratch_root}" == "/" ]]; then
      expected_root="/batch_effect/${expected_suffix}"
    else
      expected_root="${scratch_root}/batch_effect/${expected_suffix}"
    fi
    [[ "${pass}" == "${expected_pass}" &&
       "${analysis_pass}" == "${expected_pass}" &&
       "${AUDIT_METADATA_ROOT}" == "${analysis_root}" &&
       "${analysis_root}" == "${expected_root}" &&
       "${analysis_root}" = /* &&
       "${analysis_nas_root}" = /* &&
       "${analysis_log_prefix}" == "${expected_log_prefix}" ]] || {
      _audit_die "Stage 5 variant metadata has mismatched pass/root/log identity"
      return 1
    }

    case "${analysis_nas_root}" in
      */batch_effect/${expected_suffix})
        nas_base="${analysis_nas_root%/batch_effect/${expected_suffix}}"
        [[ -n "${nas_base}" ]] || nas_base="/"
        ;;
      *)
        _audit_die "Stage 5 variant NAS root is not variant-qualified"
        return 1
        ;;
    esac
    nas_target="${NAS_TARGET_DIR:-}"
    if [[ -n "${nas_target}" ]]; then
      nas_target="${nas_target%/}"
      [[ -n "${nas_target}" ]] || nas_target="/"
      [[ "${nas_target}" == "${nas_base}" ]] || {
        _audit_die "Stage 5 variant NAS root disagrees with configured NAS root"
        return 1
      }
    else
      nas_target="${nas_base}"
    fi
    [[ "${nas_target}" = /* && -d "${nas_target}" ]] || {
      _audit_die "Stage 5 variant NAS root parent is missing: ${nas_target}"
      return 1
    }
    expected_nas="${nas_target%/}"
    [[ -n "${expected_nas}" ]] || expected_nas="/"
    if [[ "${expected_nas}" == "/" ]]; then
      expected_nas="/batch_effect/${expected_suffix}"
    else
      expected_nas="${expected_nas}/batch_effect/${expected_suffix}"
    fi
    [[ "${analysis_nas_root}" == "${expected_nas}" ]] || {
      _audit_die "Stage 5 variant NAS root has the wrong suffix"
      return 1
    }
    export NAS_TARGET_DIR="${nas_target}"
    export ANALYSIS_VARIANT="${variant}"
    export ANALYSIS_ROOT="${analysis_root}"
    export ANALYSIS_NAS_ROOT="${analysis_nas_root}"
    export ANALYSIS_PASS="${analysis_pass}"
    export ANALYSIS_LOG_PREFIX="${analysis_log_prefix}"
    export PASS_ARG="${pass}"
    if [[ "${variant}" == corrected_final ]]; then
      export ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION=recovery_35row
    fi
    return 0
  fi

  # Legacy runs do not record ANALYSIS_* fields.  Recreate the historical
  # defaults from ROOT/PASS without allowing stale variant state to leak in.
  if [[ -n "${AUDIT_METADATA_ROOT}" ]]; then
    [[ "${AUDIT_METADATA_ROOT}" = /* ]] || {
      _audit_die "legacy Stage 5 ROOT metadata is not absolute"
      return 1
    }
    export ANALYSIS_ROOT="${AUDIT_METADATA_ROOT}"
  elif [[ -n "${pass}" ]]; then
    export ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/${pass}"
  else
    export ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/benchmark"
  fi
  if [[ -n "${pass}" ]]; then
    export ANALYSIS_PASS="${pass}"
    export PASS_ARG="${pass}"
    export ANALYSIS_LOG_PREFIX="execution_times_batch_effect_${pass}_"
    if [[ -n "${NAS_TARGET_DIR:-}" ]]; then
      if [[ -n "${AUDIT_METADATA_ROOT}" &&
            "${AUDIT_METADATA_ROOT}" == "${HPC_SCRATCH_DIR}"/* ]]; then
        export ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR%/}/${AUDIT_METADATA_ROOT#${HPC_SCRATCH_DIR}/}"
      else
        export ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR%/}/batch_effect/${pass}"
      fi
    fi
  else
    unset ANALYSIS_PASS PASS_ARG
    export ANALYSIS_LOG_PREFIX="execution_times_"
    if [[ -n "${NAS_TARGET_DIR:-}" ]]; then
      if [[ -n "${AUDIT_METADATA_ROOT}" &&
            "${AUDIT_METADATA_ROOT}" == "${HPC_SCRATCH_DIR}"/* ]]; then
        export ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR%/}/${AUDIT_METADATA_ROOT#${HPC_SCRATCH_DIR}/}"
      fi
    fi
  fi
}

_audit_run_metadata() {
  local metadata="${RUN_ROOT_REAL}/metadata" metadata_stage metadata_run
  local field field_count variant_count identity_count
  local root_version_count root_identity_count matrix_field_count
  local -a identity_fields=(
    ANALYSIS_VARIANT ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_PASS
    ANALYSIS_LOG_PREFIX
  )
  _audit_regular_file "${metadata}" || return 1
  metadata_stage="$(_audit_metadata_value "${metadata}" STAGE)"
  metadata_run="$(_audit_metadata_value "${metadata}" RUN_ID)"
  [[ "${metadata_stage}" == "${STAGE_ARG}" &&
     "${metadata_run}" == "${RUN_ID}" ]] || {
    _audit_die "run metadata does not match requested stage/run: ${metadata}"
    return 1
  }
  [[ "$(_audit_metadata_count "${metadata}" STAGE)" == 1 &&
     "$(_audit_metadata_count "${metadata}" RUN_ID)" == 1 ]] || {
    _audit_die "run metadata has duplicate stage/run identity fields: ${metadata}"
    return 1
  }
  for field in ROOT PASS; do
    field_count="$(_audit_metadata_count "${metadata}" "${field}")"
    [[ "${field_count}" =~ ^[01]$ ]] || {
      _audit_die "run metadata has duplicate ${field} fields: ${metadata}"
      return 1
    }
  done
  AUDIT_METADATA_ROOT="$(_audit_metadata_value "${metadata}" ROOT)"
  AUDIT_METADATA_PASS="$(_audit_metadata_value "${metadata}" PASS)"
  AUDIT_METADATA_METHODS="$(_audit_metadata_value "${metadata}" METHODS)"
  AUDIT_METADATA_PENDING_SELECTION="$(
    _audit_metadata_value "${metadata}" PENDING_SELECTION
  )"
  AUDIT_METADATA_PENDING_MD5="$(
    _audit_metadata_value "${metadata}" PENDING_SELECTION_MD5
  )"
  AUDIT_METADATA_PENDING_SIZE="$(
    _audit_metadata_value "${metadata}" PENDING_SELECTION_SIZE
  )"
  AUDIT_METADATA_VARIANT="$(_audit_metadata_value "${metadata}" ANALYSIS_VARIANT)"
  AUDIT_METADATA_ANALYSIS_ROOT="$(
    _audit_metadata_value "${metadata}" ANALYSIS_ROOT
  )"
  AUDIT_METADATA_ANALYSIS_NAS_ROOT="$(
    _audit_metadata_value "${metadata}" ANALYSIS_NAS_ROOT
  )"
  AUDIT_METADATA_ANALYSIS_PASS="$(
    _audit_metadata_value "${metadata}" ANALYSIS_PASS
  )"
  AUDIT_METADATA_ANALYSIS_LOG_PREFIX="$(
    _audit_metadata_value "${metadata}" ANALYSIS_LOG_PREFIX
  )"
  AUDIT_METADATA_ANALYSIS_ROOT_VERSION="$(
    _audit_metadata_value "${metadata}" ANALYSIS_ROOT_VERSION
  )"
  AUDIT_METADATA_ANALYSIS_ROOT_IDENTITY="$(
    _audit_metadata_value "${metadata}" ANALYSIS_ROOT_IDENTITY
  )"
  AUDIT_METADATA_METHOD_MATRIX="$(
    _audit_metadata_value "${metadata}" METHOD_MATRIX
  )"
  AUDIT_METADATA_METHOD_MATRIX_MD5="$(
    _audit_metadata_value "${metadata}" METHOD_MATRIX_MD5
  )"
  AUDIT_METADATA_METHOD_MATRIX_SIZE="$(
    _audit_metadata_value "${metadata}" METHOD_MATRIX_SIZE
  )"
  AUDIT_METADATA_METHOD_MATRIX_SHA256="$(
    _audit_metadata_value "${metadata}" METHOD_MATRIX_SHA256
  )"
  AUDIT_METADATA_METHOD_MATRIX_IDENTITY="$(
    _audit_metadata_value "${metadata}" METHOD_MATRIX_IDENTITY
  )"
  AUDIT_METADATA_DECLARED_METHOD_ROWS="$(
    _audit_metadata_value "${metadata}" DECLARED_METHOD_ROWS
  )"
  AUDIT_METADATA_PENDING_METHOD_ROWS="$(
    _audit_metadata_value "${metadata}" PENDING_METHOD_ROWS
  )"
  AUDIT_METADATA_DISPATCH_SELECTION="$(
    _audit_metadata_value "${metadata}" DISPATCH_SELECTION
  )"
  AUDIT_METADATA_DISPATCH_MD5="$(
    _audit_metadata_value "${metadata}" DISPATCH_SELECTION_MD5
  )"
  AUDIT_METADATA_DISPATCH_SIZE="$(
    _audit_metadata_value "${metadata}" DISPATCH_SELECTION_SIZE
  )"
  AUDIT_METADATA_DISPATCH_ROWS="$(
    _audit_metadata_value "${metadata}" DISPATCH_SELECTION_ROWS
  )"

  variant_count="$(_audit_metadata_count "${metadata}" ANALYSIS_VARIANT)"
  for field in "${identity_fields[@]}"; do
    field_count="$(_audit_metadata_count "${metadata}" "${field}")"
    [[ "${field_count}" =~ ^[01]$ ]] || {
      _audit_die "run metadata has duplicate ${field} fields: ${metadata}"
      return 1
    }
  done
  if [[ ${variant_count} -eq 1 ]]; then
    [[ -n "${AUDIT_METADATA_VARIANT}" ]] || {
      _audit_die "run metadata has an empty ANALYSIS_VARIANT: ${metadata}"
      return 1
    }
    for field in ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_PASS \
      ANALYSIS_LOG_PREFIX; do
      identity_count="$(_audit_metadata_count "${metadata}" "${field}")"
      [[ "${identity_count}" == 1 &&
         -n "$(_audit_metadata_value "${metadata}" "${field}")" ]] || {
        _audit_die "variant run metadata is missing ${field}: ${metadata}"
        return 1
      }
    done
  else
    for field in ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_PASS \
      ANALYSIS_LOG_PREFIX; do
      identity_count="$(_audit_metadata_count "${metadata}" "${field}")"
      [[ "${identity_count}" == 0 ]] || {
        _audit_die "legacy run metadata contains partial variant identity: ${metadata}"
        return 1
      }
    done
  fi
  # Matrix metadata is a separate contract from the dataset/view selection.
  # Corrected-final metadata always binds the recovery_35row root; matrix
  # metadata adds the independent method-scope contract.
  matrix_count="$(_audit_metadata_count "${metadata}" METHOD_MATRIX)"
  [[ "${matrix_count}" == 0 || "${matrix_count}" == 1 ]] || {
    _audit_die "run metadata has duplicate METHOD_MATRIX fields: ${metadata}"
    return 1
  }
  AUDIT_METADATA_METHOD_MATRIX_MODE=0
  for field in METHOD_MATRIX METHOD_MATRIX_MD5 METHOD_MATRIX_SIZE \
    METHOD_MATRIX_SHA256 METHOD_MATRIX_IDENTITY DECLARED_METHOD_ROWS \
    PENDING_METHOD_ROWS ANALYSIS_ROOT_VERSION ANALYSIS_ROOT_IDENTITY; do
    matrix_field_count="$(_audit_metadata_count "${metadata}" "${field}")"
    [[ "${matrix_field_count}" =~ ^[01]$ ]] || {
      _audit_die "run metadata has duplicate ${field} fields: ${metadata}"
      return 1
    }
  done
  if [[ "${matrix_count}" == 1 ]]; then
    AUDIT_METADATA_METHOD_MATRIX_MODE=1
    [[ "${STAGE_ARG}" == stage5 &&
       "${AUDIT_METADATA_VARIANT}" == corrected_final &&
       "${AUDIT_METADATA_PASS}" == corrected ]] || {
      _audit_die "METHOD_MATRIX is only valid for corrected-final Stage 5 runs"
      return 1
    }
    for field in METHOD_MATRIX METHOD_MATRIX_MD5 METHOD_MATRIX_SIZE \
      METHOD_MATRIX_SHA256 METHOD_MATRIX_IDENTITY DECLARED_METHOD_ROWS \
      PENDING_METHOD_ROWS ANALYSIS_ROOT_VERSION ANALYSIS_ROOT_IDENTITY; do
      matrix_field_count="$(_audit_metadata_count "${metadata}" "${field}")"
      [[ "${matrix_field_count}" == 1 &&
         -n "$(_audit_metadata_value "${metadata}" "${field}")" ]] || {
        _audit_die "matrix run metadata is missing ${field}: ${metadata}"
        return 1
      }
    done
    [[ "${AUDIT_METADATA_METHOD_MATRIX_MD5}" =~ ^[[:xdigit:]]{32}$ &&
       "${AUDIT_METADATA_METHOD_MATRIX_SIZE}" =~ ^[1-9][0-9]*$ &&
       "${AUDIT_METADATA_METHOD_MATRIX_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
       "${AUDIT_METADATA_METHOD_MATRIX_IDENTITY}" =~ ^[[:xdigit:]]{64}$ &&
       "${AUDIT_METADATA_METHOD_MATRIX_IDENTITY}" == "${AUDIT_METADATA_METHOD_MATRIX_SHA256}" &&
       "${AUDIT_METADATA_DECLARED_METHOD_ROWS}" =~ ^[1-9][0-9]*$ &&
       "${AUDIT_METADATA_PENDING_METHOD_ROWS}" =~ ^(0|[1-9][0-9]*)$ &&
       "${AUDIT_METADATA_ANALYSIS_ROOT_VERSION}" == recovery_35row &&
       "${AUDIT_METADATA_ANALYSIS_ROOT_IDENTITY}" == corrected_final/recovery_35row ]] || {
      _audit_die "corrected-final matrix metadata values are invalid: ${metadata}"
      return 1
    }
    matrix_field_count="$(_audit_metadata_count "${metadata}" PENDING_SELECTION)"
    [[ "${matrix_field_count}" == 1 &&
       -n "${AUDIT_METADATA_PENDING_SELECTION}" ]] || {
      _audit_die "matrix run metadata is missing PENDING_SELECTION: ${metadata}"
      return 1
    }
    for field in PENDING_SELECTION_MD5 PENDING_SELECTION_SIZE; do
      matrix_field_count="$(_audit_metadata_count "${metadata}" "${field}")"
      [[ "${matrix_field_count}" =~ ^[01]$ ]] || {
        _audit_die "matrix run metadata duplicates ${field}: ${metadata}"
        return 1
      }
    done
    if [[ "${AUDIT_METADATA_PENDING_METHOD_ROWS}" -gt 0 ]]; then
      [[ "$(_audit_metadata_count "${metadata}" PENDING_SELECTION_MD5)" == 1 &&
         "$(_audit_metadata_count "${metadata}" PENDING_SELECTION_SIZE)" == 1 &&
         -n "${AUDIT_METADATA_PENDING_MD5}" &&
         -n "${AUDIT_METADATA_PENDING_SIZE}" &&
         "${AUDIT_METADATA_PENDING_MD5}" =~ ^[[:xdigit:]]{32}$ &&
         "${AUDIT_METADATA_PENDING_SIZE}" =~ ^[1-9][0-9]*$ ]] || {
        _audit_die "matrix pending-selection metadata values are invalid: ${metadata}"
        return 1
      }
    else
      [[ -z "${AUDIT_METADATA_PENDING_MD5}" &&
         -z "${AUDIT_METADATA_PENDING_SIZE}" ]] || {
        _audit_die "NOOP matrix metadata must not claim pending bytes: ${metadata}"
        return 1
      }
    fi
  else
    for field in METHOD_MATRIX_MD5 METHOD_MATRIX_SIZE METHOD_MATRIX_SHA256 \
      METHOD_MATRIX_IDENTITY DECLARED_METHOD_ROWS PENDING_METHOD_ROWS; do
      matrix_field_count="$(_audit_metadata_count "${metadata}" "${field}")"
      [[ "${matrix_field_count}" == 0 ]] || {
        _audit_die "run metadata contains partial METHOD_MATRIX identity: ${metadata}"
        return 1
      }
    done
    root_version_count="$(_audit_metadata_count "${metadata}" ANALYSIS_ROOT_VERSION)"
    root_identity_count="$(_audit_metadata_count "${metadata}" ANALYSIS_ROOT_IDENTITY)"
    if [[ "${AUDIT_METADATA_VARIANT}" == corrected_final ]]; then
      case "${AUDIT_METADATA_ANALYSIS_ROOT}" in
        */batch_effect/corrected_final)
          [[ "${root_version_count}" == 0 &&
             "${root_identity_count}" == 0 ]] || {
            _audit_die "direct corrected-final root has unexpected replacement identity metadata"
            return 1
          }
          ;;
        */batch_effect/corrected_final/recovery_35row)
          [[ "${root_version_count}" == 1 &&
             "${root_identity_count}" == 1 &&
             "${AUDIT_METADATA_ANALYSIS_ROOT_VERSION}" == recovery_35row &&
             "${AUDIT_METADATA_ANALYSIS_ROOT_IDENTITY}" == corrected_final/recovery_35row ]] || {
            _audit_die "corrected-final root identity metadata is missing or invalid"
            return 1
          }
          ;;
        *)
          _audit_die "corrected-final metadata names an unsupported analysis root"
          return 1
          ;;
      esac
    else
      [[ "${root_version_count}" == 0 &&
         "${root_identity_count}" == 0 ]] || {
        _audit_die "run metadata contains an unexpected root identity"
        return 1
      }
    fi
  fi


  AUDIT_METADATA_EXPORT_MANIFEST="$(
    _audit_metadata_value "${metadata}" METADATA_EXPORT_MANIFEST
  )"
  AUDIT_METADATA_EXPORT_STATUS="$(
    _audit_metadata_value "${metadata}" METADATA_EXPORT_STATUS
  )"
  export_count="$(_audit_metadata_count "${metadata}" METADATA_EXPORT_MANIFEST)"
  [[ "${export_count}" == "$(_audit_metadata_count \
    "${metadata}" METADATA_EXPORT_STATUS)" ]] || {
    _audit_die "run metadata has an incomplete metadata-export identity: ${metadata}"
    return 1
  }
  if [[ -n "${AUDIT_METADATA_VARIANT}" ]]; then
    [[ "${export_count}" == 1 &&
       -n "${AUDIT_METADATA_EXPORT_MANIFEST}" &&
       -n "${AUDIT_METADATA_EXPORT_STATUS}" ]] || {
      _audit_die "variant run metadata is missing metadata-export identity: ${metadata}"
      return 1
    }
  else
    [[ "${export_count}" == 0 ]] || {
      _audit_die "legacy run metadata contains metadata-export identity: ${metadata}"
      return 1
    }
  fi


  AUDIT_METADATA_BATCH_CONTRACT_MANIFEST=""
  AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_MD5=""
  AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SIZE=""
  AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SHA256=""
  if [[ "${AUDIT_METADATA_PASS}" == corrected ]]; then
    for field in BATCH_CONTRACT_MANIFEST BATCH_CONTRACT_MANIFEST_MD5 \
      BATCH_CONTRACT_MANIFEST_SIZE BATCH_CONTRACT_MANIFEST_SHA256; do
      field_count="$(_audit_metadata_count "${metadata}" "${field}")"
      [[ "${field_count}" == 1 ]] || {
        _audit_die "corrected run metadata is missing or duplicates ${field}: ${metadata}"
        return 1
      }
    done
    AUDIT_METADATA_BATCH_CONTRACT_MANIFEST="$(
      _audit_metadata_value "${metadata}" BATCH_CONTRACT_MANIFEST
    )"
    AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_MD5="$(
      _audit_metadata_value "${metadata}" BATCH_CONTRACT_MANIFEST_MD5
    )"
    AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SIZE="$(
      _audit_metadata_value "${metadata}" BATCH_CONTRACT_MANIFEST_SIZE
    )"
    AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SHA256="$(
      _audit_metadata_value "${metadata}" BATCH_CONTRACT_MANIFEST_SHA256
    )"
  fi
  _audit_stage5_identity || return 1
}

_audit_terminal_status() {
  local terminal="${RUN_ROOT_REAL}/status/terminal" line state status_run
  _audit_regular_file "${terminal}" || return 1
  state="$(sed -n 's/^STATE=//p' "${terminal}" | sed -n '1p')"
  status_run="$(sed -n 's/^RUN_ID=//p' "${terminal}" | sed -n '1p')"
  [[ "${status_run}" == "${RUN_ID}" &&
     ( "${state}" == OK || "${state}" == NOOP_VALIDATED ) ]] || {
    _audit_die "run terminal status is not successful: ${terminal}"
    return 1
  }
  AUDIT_TERMINAL_STATE="${state}"
}

_audit_selection() {
  local columns rows
  case "${STAGE_ARG}" in
    stage2) columns=5 ;;
    stage3|stage4) columns=2 ;;
    stage5) columns=3 ;;
  esac
  ecoda_validate_run_owned_path "${SELECTION_ARG}" "${RUN_ROOT_REAL}" || {
    _audit_die "selection is not run-owned: ${SELECTION_ARG}"
    return 1
  }
  ecoda_validate_manifest "${SELECTION_ARG}" "${columns}" || return 1
  ecoda_validate_checksum "${SELECTION_ARG}" || {
    _audit_die "selection checksum is missing or invalid: ${SELECTION_ARG}"
    return 1
  }
  rows="$(wc -l < "${SELECTION_ARG}" | tr -d '[:space:]')"
  [[ "${rows}" =~ ^[1-9][0-9]*$ ]] || return 1
  AUDIT_SELECTION_ROWS="${rows}"
}

_audit_stage5_method_matrix() {
  local matrix="${AUDIT_METADATA_METHOD_MATRIX:-}"
  local matrix_real matrix_md5 matrix_size matrix_sha
  local pending="${AUDIT_METADATA_PENDING_SELECTION:-}"
  local pending_md5 pending_size pending_actual_md5 pending_actual_size
  local row_dataset row_view row_method extra
  local pending_dataset pending_view pending_method pending_extra
  local selection_dataset selection_view selection_label selection_extra
  local key matrix_keys="" selection_keys="" pending_keys="" index=0
  local pending_count=0 selection_count=0

  [[ "${STAGE_ARG}" == stage5 ]] || return 0
  if [[ "${AUDIT_METADATA_METHOD_MATRIX_MODE:-0}" != 1 ]]; then
    unset ECODA_STAGE5_METHOD_MATRIX METHOD_MATRIX
    return 0
  fi
  [[ "${AUDIT_METADATA_VARIANT:-}" == corrected_final &&
     "${AUDIT_METADATA_PASS:-}" == corrected ]] || {
    _audit_die "METHOD_MATRIX requires the corrected-final Stage 5 pass"
    return 1
  }
  [[ "${AUDIT_METADATA_METHODS}" == "prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot" ]] || {
    _audit_die "corrected-final matrix metadata has the wrong method suite"
    return 1
  }
  [[ "${matrix}" == "${RUN_ROOT_REAL}/manifests/method_matrix.tsv" ]] || {
    _audit_die "METHOD_MATRIX is not the canonical run-owned matrix path"
    return 1
  }
  _audit_regular_file "${matrix}" || return 1
  _audit_regular_file "${matrix}.md5" || return 1
  matrix_real="$(ecoda_realpath_existing "${matrix}")" || return 1
  [[ "${matrix_real}" == "${matrix}" ]] || {
    _audit_die "METHOD_MATRIX path is not canonical: ${matrix}"
    return 1
  }
  ecoda_validate_run_owned_path "${matrix}" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_run_owned_path "${matrix}.md5" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_manifest "${matrix}" 3 || return 1
  ecoda_validate_checksum "${matrix}" || {
    _audit_die "METHOD_MATRIX checksum is invalid: ${matrix}"
    return 1
  }
  matrix_md5="${ECODA_CHECKSUM_MD5}"
  matrix_size="${ECODA_CHECKSUM_SIZE}"
  matrix_sha="$(_audit_sha256_file "${matrix}")" || return 1
  [[ "${AUDIT_METADATA_METHOD_MATRIX_MD5}" == "${matrix_md5}" &&
     "${AUDIT_METADATA_METHOD_MATRIX_SIZE}" == "${matrix_size}" &&
     "${AUDIT_METADATA_METHOD_MATRIX_SHA256}" == "${matrix_sha}" &&
     "${AUDIT_METADATA_METHOD_MATRIX_IDENTITY}" == "${matrix_sha}" ]] || {
    _audit_die "METHOD_MATRIX checksum/size metadata mismatches: ${matrix}"
    return 1
  }

  # The matrix is an ordered declaration of authorized method rows.  Its
  # dataset/method cardinality is intentionally supplied by the manifest.
  while IFS=$'\t' read -r row_dataset row_view row_method extra; do
    [[ -n "${row_dataset}" && -n "${row_view}" && -n "${row_method}" &&
       "${row_view}" == batch_effect_corrected &&
       "${row_dataset}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${row_method}" =~ ^[A-Za-z0-9_.-]+$ &&
       -z "${extra}" ]] || {
      _audit_die "METHOD_MATRIX row ${index} is malformed"
      return 1
    }
    case ",prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot," in
      *,"${row_method}",*) ;;
      *) _audit_die "METHOD_MATRIX contains an unsupported method: ${row_method}"; return 1 ;;
    esac
    key="${row_dataset}|${row_view}|${row_method}"
    case " ${matrix_keys} " in
      *" ${key} "*)
        _audit_die "METHOD_MATRIX contains a duplicate triple: ${key}"
        return 1
        ;;
    esac
    matrix_keys="${matrix_keys} ${key}"
    index=$((index + 1))
  done < "${matrix}"
  [[ ${index} -gt 0 ]] || {
    _audit_die "METHOD_MATRIX is empty"
    return 1
  }
  [[ "${AUDIT_METADATA_DECLARED_METHOD_ROWS}" == "${index}" ]] || {
    _audit_die "DECLARED_METHOD_ROWS disagrees with METHOD_MATRIX"
    return 1
  }

  [[ "${SELECTION_ARG}" != "${matrix}" ]] || {
    _audit_die "dataset selection and METHOD_MATRIX cannot be the same file"
    return 1
  }
  while IFS=$'\t' read -r selection_dataset selection_view \
    selection_label selection_extra; do
    [[ -n "${selection_dataset}" &&
       "${selection_view}" == batch_effect_corrected &&
       "${selection_label}" == batch_effect_corrected &&
       -z "${selection_extra}" ]] || {
      _audit_die "corrected-final dataset selection is malformed"
      return 1
    }
    key="${selection_dataset}|${selection_view}"
    case " ${selection_keys} " in
      *" ${key} "*)
        _audit_die "corrected-final dataset selection contains a duplicate row"
        return 1
        ;;
    esac
    selection_keys="${selection_keys} ${key}"
    case " ${matrix_keys} " in
      *" ${selection_dataset}|${selection_view}|"*) ;;
      *)
        _audit_die "METHOD_MATRIX does not cover selection row: ${key}"
        return 1
        ;;
    esac
    selection_count=$((selection_count + 1))
  done < "${SELECTION_ARG}"
  [[ ${selection_count} -gt 0 ]] || {
    _audit_die "corrected-final matrix selection is empty"
    return 1
  }
  while IFS=$'\t' read -r row_dataset row_view row_method extra; do
    key="${row_dataset}|${row_view}"
    case " ${selection_keys} " in
      *" ${key} "*) ;;
      *) _audit_die "METHOD_MATRIX escapes the dataset selection: ${key}"; return 1 ;;
    esac
  done < "${matrix}"

  # PENDING_SELECTION is the derived method subset.  It must be run-owned,
  # checksummed independently, and contain only declared matrix triples.
  [[ "${pending}" == "${RUN_ROOT_REAL}/manifests/pending_selection.tsv" ]] || {
    _audit_die "matrix PENDING_SELECTION is not the canonical run-owned path"
    return 1
  }
  [[ -f "${pending}" && ! -L "${pending}" && -r "${pending}" ]] || {
    _audit_die "matrix PENDING_SELECTION is missing or unsafe: ${pending}"
    return 1
  }
  if [[ -s "${pending}" ]]; then
    _audit_regular_file "${pending}.md5" || return 1
    ecoda_validate_manifest "${pending}" 3 || return 1
    ecoda_validate_checksum "${pending}" || {
      _audit_die "matrix PENDING_SELECTION checksum is invalid: ${pending}"
      return 1
    }
    pending_actual_md5="${ECODA_CHECKSUM_MD5}"
    pending_actual_size="${ECODA_CHECKSUM_SIZE}"
    pending_md5="${AUDIT_METADATA_PENDING_MD5}"
    pending_size="${AUDIT_METADATA_PENDING_SIZE}"
    [[ "${pending_md5}" == "${pending_actual_md5}" &&
       "${pending_size}" == "${pending_actual_size}" ]] || {
      _audit_die "PENDING_SELECTION checksum/size metadata mismatches: ${pending}"
      return 1
    }
  else
    [[ "${AUDIT_TERMINAL_STATE}" == NOOP_VALIDATED &&
       "${AUDIT_METADATA_PENDING_METHOD_ROWS}" == 0 &&
       -z "${AUDIT_METADATA_PENDING_MD5}" &&
       -z "${AUDIT_METADATA_PENDING_SIZE}" ]] || {
      _audit_die "empty matrix PENDING_SELECTION is only valid for a NOOP run"
      return 1
    }
  fi
  if [[ -s "${pending}" ]]; then
    while IFS=$'\t' read -r pending_dataset pending_view \
      pending_method pending_extra; do
      [[ -n "${pending_dataset}" && -n "${pending_view}" &&
         -n "${pending_method}" && -z "${pending_extra}" ]] || {
        _audit_die "matrix PENDING_SELECTION row is malformed"
        return 1
      }
      key="${pending_dataset}|${pending_view}|${pending_method}"
      case " ${matrix_keys} " in
        *" ${key} "*) ;;
        *)
          _audit_die "PENDING_SELECTION escapes METHOD_MATRIX: ${key}"
          return 1
          ;;
      esac
      case " ${pending_keys} " in
        *" ${key} "*)
          _audit_die "PENDING_SELECTION contains a duplicate triple: ${key}"
          return 1
          ;;
      esac
      pending_keys="${pending_keys} ${key}"
      pending_count=$((pending_count + 1))
    done < "${pending}"
  fi
  [[ "${AUDIT_METADATA_PENDING_METHOD_ROWS}" == "${pending_count}" ]] || {
    _audit_die "PENDING_METHOD_ROWS disagrees with PENDING_SELECTION"
    return 1
  }
  AUDIT_MATRIX_PENDING_ROWS="${pending_count}"
}

_audit_dispatch_selection() {
  local dispatch="${AUDIT_METADATA_DISPATCH_SELECTION:-}"
  local dispatch_md5 dispatch_size dispatch_rows actual_md5 actual_size actual_rows
  local row_dataset row_view row_method extra key seen="" selection_keys=""
  [[ "${STAGE_ARG}" == stage5 ]] || return 0
  [[ -n "${dispatch}" ]] || return 0
  [[ "${dispatch}" == "${RUN_ROOT_REAL}/manifests/dispatch_selection.tsv" ]] || {
    _audit_die "dispatch selection is not the canonical run-owned path"
    return 1
  }
  _audit_regular_file "${dispatch}" || return 1
  _audit_regular_file "${dispatch}.md5" || return 1
  ecoda_validate_run_owned_path "${dispatch}" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_run_owned_path "${dispatch}.md5" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_manifest "${dispatch}" 3 || return 1
  ecoda_validate_checksum "${dispatch}" || return 1
  actual_md5="${ECODA_CHECKSUM_MD5}"
  actual_size="${ECODA_CHECKSUM_SIZE}"
  actual_rows="$(awk 'END { print NR }' "${dispatch}")" || return 1
  dispatch_md5="${AUDIT_METADATA_DISPATCH_MD5}"
  dispatch_size="${AUDIT_METADATA_DISPATCH_SIZE}"
  dispatch_rows="${AUDIT_METADATA_DISPATCH_ROWS}"
  [[ "${dispatch_md5}" =~ ^[[:xdigit:]]{32}$ &&
     "${dispatch_size}" =~ ^[1-9][0-9]*$ &&
     "${dispatch_rows}" =~ ^[1-9][0-9]*$ &&
     "${dispatch_md5}" == "${actual_md5}" &&
     "${dispatch_size}" == "${actual_size}" &&
     "${dispatch_rows}" == "${actual_rows}" ]] || {
    _audit_die "dispatch selection metadata does not match its manifest"
    return 1
  }
  while IFS=$'\t' read -r selection_dataset selection_view \
    _selection_label selection_extra; do
    [[ -n "${selection_dataset}" && -n "${selection_view}" &&
       -z "${selection_extra}" ]] || return 1
    selection_keys="${selection_keys} ${selection_dataset}|${selection_view}"
  done < "${SELECTION_ARG}"
  while IFS=$'\t' read -r row_dataset row_view row_method extra; do
    [[ -n "${row_dataset}" && -n "${row_view}" && -n "${row_method}" &&
       -z "${extra}" &&
       "${row_dataset}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${row_view}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${row_method}" =~ ^[A-Za-z0-9_.-]+$ ]] || return 1
    key="${row_dataset}|${row_view}|${row_method}"
    case " ${seen} " in
      *" ${key} "*) _audit_die "dispatch selection contains a duplicate row"; return 1 ;;
    esac
    seen="${seen} ${key}"
    case " ${selection_keys} " in
      *" ${row_dataset}|${row_view} "*) ;;
      *) _audit_die "dispatch selection escapes dataset selection"; return 1 ;;
    esac
    if [[ "${AUDIT_METADATA_METHOD_MATRIX_MODE:-0}" == 1 ]]; then
      awk -F '\t' -v ds="${row_dataset}" -v view="${row_view}" \
        -v method="${row_method}" \
        '$1 == ds && $2 == view && $3 == method && NF == 3 { found=1 }
         END { exit(found ? 0 : 1) }' \
        "${AUDIT_METADATA_METHOD_MATRIX}" || {
        _audit_die "dispatch selection escapes METHOD_MATRIX"
        return 1
      }
    else
      case ",prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot,mofa,scitd,scpoli,pilotgm,trans,zeroimp," in
        *,"${row_method}",*) ;;
        *) _audit_die "dispatch selection contains an unsupported method"; return 1 ;;
      esac
    fi
  done < "${dispatch}"
}

_audit_batch_contract_identity_path() {
  local dataset="${1:-}" view="${2:-}" method="${3:-}"
  local index row_dataset row_view row_method row_path row_md5 row_size extra
  local found=0
  index="${#AUDIT_BATCH_CONTRACT_PATHS[@]}"
  while [[ ${index} -gt 0 ]]; do
    index=$((index - 1))
    if [[ "${AUDIT_BATCH_CONTRACT_DATASETS[${index}]}" == "${dataset}" &&
          "${AUDIT_BATCH_CONTRACT_VIEWS[${index}]}" == "${view}" &&
          "${AUDIT_BATCH_CONTRACT_METHODS[${index}]}" == "${method}" ]]; then
      [[ ${found} -eq 0 ]] || return 1
      found=1
      printf '%s' "${AUDIT_BATCH_CONTRACT_PATHS[${index}]}"
    fi
  done
  [[ ${found} -eq 1 ]]
}
_audit_batch_contract_identity_for_label() {
  local dataset="${1:-}" view="${2:-}" label="${3:-}"
  local source_path="${4:-}" method_id model_id
  [[ -n "${dataset}" && -n "${view}" && -n "${label}" ]] || {
    _audit_die "corrected batch-contract identity requires dataset, view, and method"
    return 1
  }
  # Manifest and matrix rows use the stable shell labels.  The Python/R
  # contract builders accept only the canonical semantic method IDs, so map
  # each label together with its method-specific model before invoking them.
  ecoda_corrected_batch_method_policy "${label}" || {
    _audit_die "unsupported corrected batch method label: ${label}"
    return 1
  }
  method_id="${ECODA_CORRECTED_BATCH_METHOD_ID}"
  model_id="${ECODA_CORRECTED_BATCH_MODEL_ID}"
  [[ -n "${method_id}" && -n "${model_id}" ]] || {
    _audit_die "corrected batch method policy produced an empty identity: ${label}"
    return 1
  }
  if [[ "${AUDIT_METADATA_VARIANT:-}" == corrected_final ]]; then
    ecoda_batch_contract_identity "${DATASETS_JSON_FILE}" "${dataset}" \
      "${view}" "${method_id}" "${model_id}"
  else
    [[ -n "${source_path}" ]] || {
      _audit_die "direct-root corrected batch contract requires a source H5AD: ${label}"
      return 1
    }
    ecoda_batch_contract_identity "${DATASETS_JSON_FILE}" "${dataset}" \
      "${view}" "${method_id}" "${model_id}" "${source_path}"
  fi
}


_audit_batch_contract_manifest() {
  local manifest="${AUDIT_METADATA_BATCH_CONTRACT_MANIFEST:-}"
  local row_dataset row_view row_method row_path row_md5 row_size extra
  local expected_path expected_identity safe key actual_count=0 expected_count=0
  local manifest_real row_path_real
  local source_path output_name
  local duplicate_keys="" dataset view label method contract_method found index
  local matrix_mode="${AUDIT_METADATA_METHOD_MATRIX_MODE:-0}"
  local -a configured_methods=()
  local -a expected_methods=()
  local -a scope_datasets=() scope_views=()
  local -a actual_datasets=() actual_views=() actual_methods=() actual_paths=()
  [[ "${STAGE_ARG}" == stage5 && "${AUDIT_METADATA_PASS:-}" == corrected ]] ||
    return 0
  export ECODA_SOURCE_ROOT="${AUDIT_SOURCE_ROOT}"
  PYTHON_BIN="$(_audit_python_binary)" || return 1
  export PYTHON_BIN
  export PROJECT_ROOT="${AUDIT_SOURCE_ROOT}"
  export DATASETS_JSON_FILE="${AUDIT_SOURCE_ROOT}/datasets.json"
  [[ -n "${manifest}" &&
     -f "${manifest}" && ! -L "${manifest}" && -r "${manifest}" ]] || {
    _audit_die "corrected run batch-contract manifest is missing or unsafe"
    return 1
  }
  manifest_real="$(ecoda_realpath_existing "${manifest}")" || return 1
  [[ "${manifest_real}" == "${RUN_ROOT_REAL}/manifests/batch_contract.tsv" ]] || {
    _audit_die "corrected run batch-contract manifest is not run-owned"
    return 1
  }
  ecoda_validate_run_owned_path "${manifest}" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_manifest "${manifest}" 6 || return 1
  ecoda_validate_checksum "${manifest}" || {
    _audit_die "corrected run batch-contract manifest checksum is invalid"
    return 1
  }
  [[ "${AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_MD5}" == "${ECODA_CHECKSUM_MD5}" &&
     "${AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SIZE}" == "${ECODA_CHECKSUM_SIZE}" ]] || {
    _audit_die "corrected run batch-contract manifest checksum metadata mismatches"
    return 1
  }
  [[ "${AUDIT_METADATA_BATCH_CONTRACT_MANIFEST_SHA256}" == "$(ecoda_sha256_file "${manifest}")" ]] || {
    _audit_die "corrected run batch-contract manifest SHA-256 mismatches"
    return 1
  }
  [[ -n "${AUDIT_METADATA_METHODS}" ]] || {
    _audit_die "corrected run metadata has no selected methods"
    return 1
  }
  ecoda_split_csv "${AUDIT_METADATA_METHODS}" || return 1
  for method in "${ECODA_ARRAY[@]}"; do
    [[ "${method}" != _ecoda_none_ ]] || return 1
    ecoda_corrected_batch_method_policy "${method}" || return 1
    configured_methods+=("${method}")
  done
  while IFS=$'\t' read -r dataset view label extra; do
    [[ -n "${dataset}" && -n "${view}" && -n "${label}" && -z "${extra}" &&
       "${dataset}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${view}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${label}" == batch_effect_corrected ]] || {
      _audit_die "corrected selection has an invalid dataset/view row"
      return 1
    }
    key="${dataset}|${view}"
    case " ${duplicate_keys} " in
      *" ${key} "*)
        _audit_die "corrected selection contains a duplicate dataset/view"
        return 1
        ;;
    esac
    duplicate_keys="${duplicate_keys} ${key}"
    scope_datasets+=("${dataset}")
    scope_views+=("${view}")
  done < "${SELECTION_ARG}"
  [[ ${#scope_datasets[@]} -gt 0 ]] || return 1
  for index in "${!scope_datasets[@]}"; do
    ecoda_validate_corrected_batch_columns \
      "${DATASETS_JSON_FILE}" "${scope_datasets[${index}]}" \
      "${scope_views[${index}]}" || return 1
  done
  AUDIT_BATCH_CONTRACT_DATASETS=()
  AUDIT_BATCH_CONTRACT_VIEWS=()
  AUDIT_BATCH_CONTRACT_METHODS=()
  AUDIT_BATCH_CONTRACT_PATHS=()
  while IFS=$'\t' read -r row_dataset row_view row_method row_path row_md5 row_size extra; do
    [[ -n "${row_dataset}" && -n "${row_view}" && -n "${row_method}" &&
       -n "${row_path}" && -n "${row_md5}" && -n "${row_size}" && -z "${extra}" &&
       "${row_dataset}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${row_view}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${row_method}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${row_path}" = /* && "${row_path}" != *$'\n'* &&
       "${row_path}" != *$'\t'* &&
       "${row_md5}" =~ ^[[:xdigit:]]{32}$ &&
       "${row_size}" =~ ^[1-9][0-9]*$ ]] || {
      _audit_die "corrected batch-contract manifest row is malformed"
      return 1
    }
    key="${row_dataset}|${row_view}|${row_method}"
    case " ${duplicate_keys} " in
      *" ${key} "*)
        _audit_die "corrected batch-contract manifest contains duplicate ${key}"
        return 1
        ;;
    esac
    duplicate_keys="${duplicate_keys} ${key}"
    safe="$(_ecoda_safe_component "${row_dataset}__${row_view}__${row_method}")" ||
      return 1
    expected_path="${RUN_ROOT_REAL}/manifests/batch_contracts/${safe}.json"
    row_path_real="$(ecoda_realpath_existing "${row_path}")" || return 1
    [[ "${row_path_real}" == "${expected_path}" ]] || {
      _audit_die "corrected batch-contract identity path is not run-owned: ${key}"
      return 1
    }
    ecoda_validate_checksum "${row_path}" || {
      _audit_die "corrected batch-contract identity checksum is invalid: ${key}"
      return 1
    }
    [[ "${row_md5}" == "${ECODA_CHECKSUM_MD5}" &&
       "${row_size}" == "${ECODA_CHECKSUM_SIZE}" ]] || return 1
    output_name="$(ecoda_view_output_name "${row_dataset}" "${row_view}")" ||
      return 1
    source_path="${HPC_SCRATCH_DIR}/${row_dataset}/output/${output_name}"
    [[ -s "${source_path}" ]] || {
      _audit_die "corrected source H5AD is missing or empty: ${source_path}"
      return 1
    }
    expected_identity="$(
      _audit_batch_contract_identity_for_label \
        "${row_dataset}" "${row_view}" "${row_method}" "${source_path}"
    )" || return 1
    cmp -s "${row_path}" <(printf '%s\n' "${expected_identity}") || {
      _audit_die "corrected batch-contract identity or validated summary mismatches ${key}"
      return 1
    }
    actual_datasets+=("${row_dataset}")
    actual_views+=("${row_view}")
    actual_methods+=("${row_method}")
    actual_paths+=("${row_path}")
    actual_count=$((actual_count + 1))
  done < "${manifest}"
  for index in "${!scope_datasets[@]}"; do
    dataset="${scope_datasets[${index}]}"
    view="${scope_views[${index}]}"
    if [[ "${matrix_mode}" == 1 ]]; then
      expected_methods=(preprocess)
      while IFS=$'\t' read -r row_dataset row_view row_method extra; do
        [[ "${row_dataset}" == "${dataset}" &&
           "${row_view}" == "${view}" ]] &&
          expected_methods+=("${row_method}")
      done < "${AUDIT_METADATA_METHOD_MATRIX}"
      [[ ${#expected_methods[@]} -gt 1 ]] || {
        _audit_die "METHOD_MATRIX has no methods for ${dataset}/${view}"
        return 1
      }
    else
      expected_methods=(preprocess "${configured_methods[@]}")
    fi
    for contract_method in "${expected_methods[@]}"; do
      found=0
      for method_index in "${!actual_datasets[@]}"; do
        if [[ "${actual_datasets[${method_index}]}" == "${dataset}" &&
              "${actual_views[${method_index}]}" == "${view}" &&
              "${actual_methods[${method_index}]}" == "${contract_method}" ]]; then
          found=1
          break
        fi
      done
      [[ ${found} -eq 1 ]] || {
        _audit_die "corrected batch-contract identity is missing ${dataset}/${view}/${contract_method}"
        return 1
      }
      expected_count=$((expected_count + 1))
    done
  done
  [[ ${actual_count} -eq ${expected_count} ]] || {
    _audit_die "corrected batch-contract manifest has unexpected rows"
    return 1
  }
  AUDIT_BATCH_CONTRACT_DATASETS=("${actual_datasets[@]}")
  AUDIT_BATCH_CONTRACT_VIEWS=("${actual_views[@]}")
  AUDIT_BATCH_CONTRACT_METHODS=("${actual_methods[@]}")
  AUDIT_BATCH_CONTRACT_PATHS=("${actual_paths[@]}")
}

_audit_stage5_variant_view() {
  local view="$1"
  case "${AUDIT_METADATA_VARIANT:-}" in
    "") ;;
    final)
      [[ "${view}" == batch_effect_uncorrected ]] || {
        _audit_die "final Stage 5 selection must use batch_effect_uncorrected"
        return 1
      }
      ;;
    corrected_final)
      [[ "${view}" == batch_effect_corrected ]] || {
        _audit_die "corrected_final Stage 5 selection must use batch_effect_corrected"
        return 1
      }
      ;;
    *)
      _audit_die "unsupported Stage 5 analysis variant"
      return 1
      ;;
  esac
}

_audit_stage5_scope_selection() {
  local dataset view label extra scope
  local pending_dataset pending_view pending_label pending_extra pending_key
  local view_rows=0 method_rows=0 scope_keys="" key
  while IFS=$'\t' read -r dataset view label extra; do
    [[ -n "${dataset}" && -n "${view}" && -n "${label}" && -z "${extra}" &&
       "${dataset}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${view}" =~ ^[A-Za-z0-9_.-]+$ ]] || {
      _audit_die "Stage 5 selection row is malformed: ${SELECTION_ARG}"
      return 1
    }
    ecoda_dataset_exists "${dataset}" &&
      ecoda_view_exists "${dataset}" "${view}" || {
      _audit_die "Stage 5 selection dataset/view is not declared: ${dataset}/${view}"
      return 1
    }
    _audit_stage5_variant_view "${view}" || return 1
    case "${label}" in
      benchmark_analysis|batch_effect_uncorrected|batch_effect_corrected)
        [[ "${label}" == "${view}" ]] || {
          _audit_die "Stage 5 view scope does not match its view column: ${SELECTION_ARG}"
          return 1
        }
        key="${dataset}|${view}"
        case " ${scope_keys} " in
          *" ${key} "*)
            _audit_die "Stage 5 selection contains a duplicate view scope: ${dataset}/${view}"
            return 1
            ;;
        esac
        scope_keys="${scope_keys} ${key}"
        view_rows=1
        ;;
      *) method_rows=1 ;;
    esac
  done < "${SELECTION_ARG}"
  [[ ${view_rows} -eq 0 || ${method_rows} -eq 0 ]] || {
    _audit_die "Stage 5 selection mixes view and method scope rows: ${SELECTION_ARG}"
    return 1
  }
  if [[ ${view_rows} -eq 1 ]]; then
    if [[ -z "${AUDIT_METADATA_PENDING_SELECTION}" &&
          "${AUDIT_TERMINAL_STATE}" != NOOP_VALIDATED ]]; then
      _audit_die "Stage 5 pending selection metadata is missing: ${RUN_ROOT_REAL}/metadata"
      return 1
    fi
    scope="${AUDIT_METADATA_PENDING_SELECTION:-${RUN_ROOT_REAL}/manifests/pending_selection.tsv}"
    [[ -f "${scope}" && ! -L "${scope}" && -r "${scope}" ]] || {
      _audit_die "Stage 5 pending selection is missing or unsafe: ${scope}"
      return 1
    }
    ecoda_validate_run_owned_path "${scope}" "${RUN_ROOT_REAL}" || return 1
    if [[ -s "${scope}" ]]; then
      [[ -f "${scope}.md5" && ! -L "${scope}.md5" &&
         -r "${scope}.md5" ]] || {
        _audit_die "Stage 5 pending selection checksum is missing or unsafe: ${scope}.md5"
        return 1
      }
      ecoda_validate_checksum "${scope}" || {
        _audit_die "Stage 5 pending selection checksum is invalid: ${scope}"
        return 1
      }
      if [[ -n "${AUDIT_METADATA_PENDING_SELECTION}" ]]; then
        [[ "${AUDIT_METADATA_PENDING_MD5}" == "${ECODA_CHECKSUM_MD5}" &&
           "${AUDIT_METADATA_PENDING_SIZE}" == "${ECODA_CHECKSUM_SIZE}" ]] || {
          _audit_die "Stage 5 pending selection metadata checksum does not match: ${scope}"
          return 1
        }
      fi
      ecoda_validate_manifest "${scope}" 3 || return 1
      while IFS=$'\t' read -r pending_dataset pending_view pending_label pending_extra; do
        [[ -n "${pending_dataset}" && -n "${pending_view}" &&
           -n "${pending_label}" && -z "${pending_extra}" &&
           "${pending_dataset}" =~ ^[A-Za-z0-9_.-]+$ &&
           "${pending_view}" =~ ^[A-Za-z0-9_.-]+$ &&
           "${pending_label}" =~ ^[A-Za-z0-9_.-]+$ ]] || {
          _audit_die "Stage 5 pending selection row is malformed: ${scope}"
          return 1
        }
        ecoda_dataset_exists "${pending_dataset}" &&
          ecoda_view_exists "${pending_dataset}" "${pending_view}" || {
          _audit_die "Stage 5 pending dataset/view is not declared: ${pending_dataset}/${pending_view}"
          return 1
        }
        _audit_stage5_variant_view "${pending_view}" || return 1
        case "${pending_label}" in
          benchmark_analysis|batch_effect_uncorrected|batch_effect_corrected)
            _audit_die "Stage 5 pending selection uses a view as method scope: ${scope}"
            return 1
            ;;
        esac
        pending_key="${pending_dataset}|${pending_view}"
        case " ${scope_keys} " in
          *" ${pending_key} "*) ;;
          *)
            _audit_die "Stage 5 pending selection escapes supplied view scope: ${scope}"
            return 1
            ;;
        esac
      done < "${scope}"
    else
      [[ "${AUDIT_TERMINAL_STATE}" == NOOP_VALIDATED ]] || {
        _audit_die "Stage 5 pending selection is empty for a non-NOOP run: ${scope}"
        return 1
      }
    fi
  else
    scope="${SELECTION_ARG}"
  fi
  AUDIT_STAGE5_SCOPE_SELECTION="${scope}"
}

_audit_scheduler_ids() {
  local scheduler_file="${RUN_ROOT_REAL}/manifests/scheduler_ids.tsv"
  local kind scheduler_id extra id_list="" rows row
  local field1 field2 field3 field4 field5 display_id raw_id state exit_code
  local normalized_state match_id requested_id found_id found root_row
  local manifest_seen="" accounting_seen=""
  local requested_ids=() found_ids=() root_success_ids="" root_failed_ids=""
  [[ -f "${scheduler_file}" && ! -L "${scheduler_file}" &&
     -r "${scheduler_file}" ]] || {
    _audit_die "scheduler ID manifest is missing: ${scheduler_file}"
    return 1
  }
  if [[ ! -s "${scheduler_file}" ]]; then
    [[ "${AUDIT_TERMINAL_STATE}" == NOOP_VALIDATED ||
       "${AUDIT_TERMINAL_STATE}" == OK ]] || return 1
    return 0
  fi
  while IFS=$'\t' read -r kind scheduler_id extra; do
    [[ -n "${kind}" && "${kind}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${scheduler_id}" =~ ^[0-9]+$ && -z "${extra}" ]] || {
      _audit_die "scheduler ID manifest row is malformed: ${scheduler_file}"
      return 1
    }
    case " ${manifest_seen} " in
      *" ${scheduler_id} "*)
        _audit_die "scheduler ID manifest contains a duplicate: ${scheduler_id}"
        return 1
        ;;
    esac
    manifest_seen="${manifest_seen} ${scheduler_id}"
    requested_ids+=("${scheduler_id}")
    if [[ -z "${id_list}" ]]; then
      id_list="${scheduler_id}"
    else
      id_list+=",${scheduler_id}"
    fi
  done < "${scheduler_file}"
  [[ ${#requested_ids[@]} -gt 0 ]] || return 1
  if ! rows="$(sacct -n -P -X -j "${id_list}" \
    --format=JobID,JobIDRaw,State,ExitCode 2>/dev/null)"; then
    _audit_die "scheduler accounting query failed for requested IDs: ${id_list}"
    return 1
  fi
  [[ -n "${rows//[[:space:]]/}" ]] || {
    _audit_die "scheduler accounting is missing for requested IDs: ${id_list}"
    return 1
  }
  while IFS= read -r row; do
    [[ -n "${row}" ]] || {
      _audit_die "scheduler accounting returned a blank row"
      return 1
    }
    IFS='|' read -r field1 field2 field3 field4 field5 <<< "${row}"
    [[ -z "${field5}" ]] || {
      _audit_die "scheduler accounting row has extra fields: ${row}"
      return 1
    }
    if [[ -n "${field4}" ]]; then
      display_id="${field1}"
      raw_id="${field2}"
      state="${field3}"
      exit_code="${field4}"
    else
      # Keep compatibility with three-column local accounting stubs. Production
      # Slurm output uses JobID plus JobIDRaw so array children can be tied to
      # their requested array root without mistaking their numeric raw IDs for
      # independent jobs.
      display_id="${field1}"
      raw_id="${field1%%_*}"
      state="${field2}"
      exit_code="${field3}"
    fi
    [[ "${display_id}" =~ ^[0-9]+(_[0-9]+)?$ &&
       "${raw_id}" =~ ^[0-9]+$ &&
       -n "${state}" && -n "${exit_code}" ]] || {
      _audit_die "scheduler accounting row is malformed: ${row}"
      return 1
    }
    case " ${accounting_seen} " in
      *" ${display_id} "*)
        _audit_die "scheduler accounting returned a duplicate job row: ${display_id}"
        return 1
        ;;
    esac
    accounting_seen="${accounting_seen} ${display_id}"
    match_id=""
    root_row=0
    for requested_id in "${requested_ids[@]}"; do
      if [[ "${raw_id}" == "${requested_id}" ]]; then
        match_id="${requested_id}"
        root_row=1
        break
      fi
    done
    if [[ -z "${match_id}" ]]; then
      for requested_id in "${requested_ids[@]}"; do
        if [[ "${display_id}" == "${requested_id}" ||
              "${display_id}" =~ ^${requested_id}_[0-9]+$ ]]; then
          match_id="${requested_id}"
          break
        fi
      done
    fi
    [[ -n "${match_id}" ]] || {
      _audit_die "scheduler accounting returned an unexpected job row: ${display_id}"
      return 1
    }
    if [[ ${root_row} -eq 0 && "${display_id}" == "${match_id}" ]]; then
      root_row=1
    fi
    if [[ ${root_row} -eq 0 ]]; then
      # Array task rows are descendants of a requested root. Their state is
      # retained by Slurm accounting and intentionally is not promoted to a
      # separate scheduler requirement: an OOM task may have been superseded
      # by the watchdog's recorded retry array.
      continue
    fi
    found_ids+=("${match_id}")
    normalized_state="${state%%+*}"
    case "${normalized_state}" in
      COMPLETED)
        if [[ "${exit_code}" == 0:0* ]]; then
          case " ${root_success_ids} " in
            *" ${match_id} "*) ;;
            *) root_success_ids="${root_success_ids} ${match_id}" ;;
          esac
        else
          case " ${root_failed_ids} " in
            *" ${match_id} "*) ;;
            *) root_failed_ids="${root_failed_ids} ${match_id}" ;;
          esac
        fi
        ;;
      *)
        case " ${root_failed_ids} " in
          *" ${match_id} "*) ;;
          *) root_failed_ids="${root_failed_ids} ${match_id}" ;;
        esac
        ;;
    esac
  done <<< "${rows}"
  for requested_id in "${requested_ids[@]}"; do
    found=0
    for found_id in "${found_ids[@]}"; do
      if [[ "${found_id}" == "${requested_id}" ]]; then
        found=1
        break
      fi
    done
    [[ ${found} -eq 1 ]] || {
      _audit_die "scheduler accounting is missing a successful root row for ${requested_id}"
      return 1
    }
    case " ${root_failed_ids} " in
      *" ${requested_id} "*)
        _audit_die "scheduler ID is not a successful terminal job: ${requested_id}"
        return 1
        ;;
    esac
    case " ${root_success_ids} " in
      *" ${requested_id} "*) ;;
      *)
        _audit_die "scheduler ID has no successful terminal root row: ${requested_id}"
        return 1
        ;;
    esac
  done
}

_audit_record_for_path() {
  local path="$1" owner_run="${2:-${RUN_ID}}" record producer
  record="$(ecoda_artifact_record_path "${path}" "${owner_run}")" || return 1
  _audit_regular_file "${record}" || return 1
  producer="$(sed -n 's/^PRODUCER=//p' "${record}" | sed -n '1p')"
  [[ -n "${producer}" ]] || return 1
  ecoda_validate_artifact_record "${path}" "${producer}" "${owner_run}" >/dev/null || {
    _audit_die "selected artifact record is invalid: ${record}"
    return 1
  }
}

_audit_python_binary() {
  local binary="${PYTHON_BIN:-}"
  if [[ -z "${binary}" ]]; then
    if command -v python3 >/dev/null 2>&1; then
      binary=python3
    elif command -v python >/dev/null 2>&1; then
      binary=python
    else
      _audit_die "python is required for selected artifact semantic validation"
      return 1
    fi
  fi
  command -v "${binary}" >/dev/null 2>&1 || {
    _audit_die "configured Python interpreter is unavailable: ${binary}"
    return 1
  }
  printf '%s' "${binary}"
}

_audit_benchmark_h5ad() {
  local path="$1" view="$2" method="$3" identity_path="${4:-}"
  local python_bin validator variant="${AUDIT_METADATA_VARIANT:-}"
  local -a validator_args
  python_bin="$(_audit_python_binary)" || return 1
  validator="${AUDIT_SOURCE_ROOT}/src/utils/py/benchmark_h5ad_contract.py"
  [[ -r "${validator}" ]] || {
    _audit_die "benchmark H5AD validator is missing: ${validator}"
    return 1
  }
  validator_args=(--path "${path}" --view "${view}" --method "${method}")
  [[ -n "${identity_path}" ]] &&
    validator_args+=(--expected-batch-contract "${identity_path}")
  if [[ "${variant}" == corrected_final &&
        "${view}" == batch_effect_corrected &&
        "${method}" == "Stage 3 preprocessing" ]]; then
    validator_args+=(--allow-missing-corrected-summary)
  fi
  "${python_bin}" "${validator}" "${validator_args[@]}" >/dev/null 2>&1 || {
    _audit_die "selected H5AD contract is invalid: ${path}"
    return 1
  }
}

_audit_corrected_source_contracts() {
  local dataset view source_path identity_path output_name key seen=""
  local index
  [[ "${STAGE_ARG}" == stage5 &&
     "${AUDIT_METADATA_PASS:-}" == corrected ]] || return 0
  for index in "${!AUDIT_BATCH_CONTRACT_DATASETS[@]}"; do
    dataset="${AUDIT_BATCH_CONTRACT_DATASETS[${index}]}"
    view="${AUDIT_BATCH_CONTRACT_VIEWS[${index}]}"
    key="${dataset}|${view}"
    case " ${seen} " in *" ${key} "*) continue ;; esac
    seen="${seen} ${key}"
    output_name="$(ecoda_view_output_name "${dataset}" "${view}")" || return 1
    source_path="${HPC_SCRATCH_DIR}/${dataset}/output/${output_name}"
    identity_path="$(
      _audit_batch_contract_identity_path "${dataset}" "${view}" preprocess
    )" || return 1
    ecoda_validate_checksum "${source_path}" || {
      _audit_die "corrected source H5AD checksum is invalid: ${source_path}"
      return 1
    }
    if [[ "${AUDIT_METADATA_VARIANT:-}" == corrected_final ]]; then
      _audit_benchmark_h5ad "${source_path}" "${view}" \
        "Stage 3 preprocessing" "${identity_path}" || return 1
    else
      _audit_benchmark_h5ad "${source_path}" "${view}" \
        "Stage 5 corrected source" "${identity_path}" || return 1
    fi
  done
}

_audit_metadata_export_manifest() {
  local manifest="${AUDIT_METADATA_EXPORT_MANIFEST:-}"
  local status_report="${AUDIT_METADATA_EXPORT_STATUS:-}"
  local expected_view="batch_effect_${AUDIT_METADATA_ANALYSIS_PASS:-}"
  local manifest_real status_real
  local dataset view input output extra key name expected_input expected_output
  local row_count=0 actual_count=0 found
  local status_field
  local expected_keys="" actual_keys=""
  [[ "${STAGE_ARG}" == stage5 &&
     -n "${AUDIT_METADATA_VARIANT:-}" ]] || return 0
  export ECODA_SOURCE_ROOT="${AUDIT_SOURCE_ROOT}"
  export PROJECT_ROOT="${AUDIT_SOURCE_ROOT}"
  export DATASETS_JSON_FILE="${AUDIT_SOURCE_ROOT}/datasets.json"
  [[ "${manifest}" == "${RUN_ROOT_REAL}/manifests/metadata_export.tsv" &&
     "${status_report}" == "${RUN_ROOT_REAL}/status/metadata_export.report" ]] || {
    _audit_die "metadata-export paths are not run-owned"
    return 1
  }
  [[ -f "${manifest}" && ! -L "${manifest}" && -r "${manifest}" ]] || {
    _audit_die "metadata-export manifest is missing or unsafe: ${manifest}"
    return 1
  }
  [[ -f "${status_report}" && ! -L "${status_report}" &&
     -r "${status_report}" ]] || {
    _audit_die "metadata-export status is missing or unsafe: ${status_report}"
    return 1
  }
  _audit_regular_file "${manifest}.md5" || return 1
  _audit_regular_file "${status_report}.md5" || return 1
  manifest_real="$(ecoda_realpath_existing "${manifest}")" || return 1
  status_real="$(ecoda_realpath_existing "${status_report}")" || return 1
  [[ "${manifest_real}" == "${manifest}" &&
     "${status_real}" == "${status_report}" ]] || {
    _audit_die "metadata-export identity paths are not canonical"
    return 1
  }
  ecoda_validate_run_owned_path "${manifest}" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_run_owned_path "${manifest}.md5" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_run_owned_path "${status_report}" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_run_owned_path "${status_report}.md5" "${RUN_ROOT_REAL}" || return 1
  ecoda_validate_manifest "${manifest}" 4 || return 1
  ecoda_validate_checksum "${manifest}" || {
    _audit_die "metadata-export manifest checksum is invalid: ${manifest}"
    return 1
  }

  # The exporter is declared for each selected dataset/view, independently of
  # pending method rows.  Build that exact set from the caller's selection;
  # never discover additional rows from the run root.
  while IFS=$'\t' read -r dataset view _label extra; do
    [[ -n "${dataset}" && -n "${view}" && -z "${extra}" ]] || return 1
    _audit_stage5_variant_view "${view}" || return 1
    key="${dataset}|${view}"
    case " ${expected_keys} " in
      *" ${key} "*)
        continue
        ;;
    esac
    expected_keys="${expected_keys} ${key}"
    row_count=$((row_count + 1))
  done < "${SELECTION_ARG}"
  [[ ${row_count} -gt 0 ]] || return 1

  while IFS=$'\t' read -r dataset view input output extra; do
    [[ -n "${dataset}" && -n "${view}" && -n "${input}" &&
       -n "${output}" && -z "${extra}" &&
       "${dataset}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${view}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${input}" = /* && "${output}" = /* ]] || {
      _audit_die "metadata-export manifest row is malformed: ${manifest}"
      return 1
    }
    _audit_stage5_variant_view "${view}" || return 1
    key="${dataset}|${view}"
    case " ${expected_keys} " in
      *" ${key} "*) ;;
      *)
        _audit_die "metadata-export row escapes selected scope: ${key}"
        return 1
        ;;
    esac
    case " ${actual_keys} " in
      *" ${key} "*)
        _audit_die "metadata-export manifest contains a duplicate: ${key}"
        return 1
        ;;
    esac
    actual_keys="${actual_keys} ${key}"
    name="$(ecoda_view_output_name "${dataset}" "${view}")" || return 1
    expected_input="${HPC_SCRATCH_DIR}/${dataset}/output/${name}"
    expected_output="${ANALYSIS_ROOT}/metadata/${dataset}_sample_metadata.feather"
    [[ "${input}" == "${expected_input}" &&
       "${output}" == "${expected_output}" ]] || {
      _audit_die "metadata-export row is not bound to the selected variant path: ${key}"
      return 1
    }
    _audit_regular_file "${input}" || return 1
    _audit_regular_file "${output}" || return 1
    ecoda_validate_checksum "${input}" || {
      _audit_die "metadata-export input H5AD checksum is invalid: ${input}"
      return 1
    }
    ecoda_validate_checksum "${output}" || {
      _audit_die "metadata-export Feather checksum is invalid: ${output}"
      return 1
    }
    actual_count=$((actual_count + 1))
  done < "${manifest}"
  [[ ${actual_count} -eq ${row_count} ]] || {
    _audit_die "metadata-export manifest does not match selected datasets"
    return 1
  }

  for status_field in STATE RUN_ID ANALYSIS_VARIANT ANALYSIS_ROOT \
    ANALYSIS_NAS_ROOT ANALYSIS_PASS ANALYSIS_LOG_PREFIX MANIFEST COUNT PENDING; do
    [[ "$(_audit_metadata_count "${status_report}" "${status_field}")" == 1 ]] || {
      _audit_die "metadata-export status is missing or duplicates ${status_field}"
      return 1
    }
  done
  [[ "$(_audit_metadata_value "${status_report}" STATE)" == OK &&
     "$(_audit_metadata_value "${status_report}" RUN_ID)" == "${RUN_ID}" &&
     "$(_audit_metadata_value "${status_report}" ANALYSIS_VARIANT)" == \
       "${AUDIT_METADATA_VARIANT}" &&
     "$(_audit_metadata_value "${status_report}" ANALYSIS_ROOT)" == \
       "${AUDIT_METADATA_ANALYSIS_ROOT}" &&
     "$(_audit_metadata_value "${status_report}" ANALYSIS_NAS_ROOT)" == \
       "${AUDIT_METADATA_ANALYSIS_NAS_ROOT}" &&
     "$(_audit_metadata_value "${status_report}" ANALYSIS_PASS)" == \
       "${AUDIT_METADATA_ANALYSIS_PASS}" &&
     "$(_audit_metadata_value "${status_report}" ANALYSIS_LOG_PREFIX)" == \
       "${AUDIT_METADATA_ANALYSIS_LOG_PREFIX}" &&
     "$(_audit_metadata_value "${status_report}" MANIFEST)" == "${manifest}" &&
     "$(_audit_metadata_value "${status_report}" COUNT)" == "${row_count}" &&
     "$(_audit_metadata_value "${status_report}" PENDING)" =~ ^[0-9]+$ ]] || {
    _audit_die "metadata-export status identity mismatches run metadata"
    return 1
  }
  ecoda_validate_checksum "${status_report}" || {
    _audit_die "metadata-export status checksum is invalid: ${status_report}"
    return 1
  }
}


_audit_annotation_artifact() {
  local path="$1" option="$2" python_bin validator
  python_bin="$(_audit_python_binary)" || return 1
  validator="${AUDIT_SOURCE_ROOT}/src/utils/py/annotation_contract.py"
  [[ -r "${validator}" ]] || {
    _audit_die "annotation artifact validator is missing: ${validator}"
    return 1
  }
  "${python_bin}" "${validator}" "--${option}" "${path}" \
    --sidecar-validated >/dev/null 2>&1 || {
    _audit_die "selected annotation contract is invalid: ${path}"
    return 1
  }
}

_audit_matrix_feather() {
  local path="$1" owner_run="$2" identity_path="${3:-}"
  local python_bin validator variant="${AUDIT_METADATA_VARIANT:-}"
  local -a validator_args
  python_bin="$(_audit_python_binary)" || return 1
  validator="${AUDIT_SOURCE_ROOT}/src/5_run_benchmark_methods/matrix_artifact_validator.py"
  [[ -r "${validator}" ]] || {
    _audit_die "matrix artifact validator is missing: ${validator}"
    return 1
  }
  validator_args=(--artifact "${path}" --producer-run-id "${owner_run}")
  [[ -n "${variant}" ]] &&
    validator_args+=(--analysis-variant "${variant}")
  if [[ "${AUDIT_METADATA_PASS:-}" == corrected ]]; then
    validator_args+=(--batch --batch-pass corrected)
  fi
  [[ -n "${identity_path}" ]] &&
    validator_args+=(--expected-batch-contract "${identity_path}")
  "${python_bin}" "${validator}" "${validator_args[@]}" >/dev/null 2>&1 || {
    _audit_die "selected Feather contract is invalid: ${path}"
    return 1
  }
}
_audit_benchmark_rds() {
  local path="$1" dataset="$2" view="$3" label="$4" identity_path="${5:-}"
  local rscript="${PIXI_RSCRIPT:-Rscript}"
  local validator="${AUDIT_SOURCE_ROOT}/src/5_run_benchmark_methods/validate_benchmark_rds_contract.R"
  local variant="${AUDIT_METADATA_VARIANT:-}"
  local -a rscript_cmd rds_args
  rscript_cmd=()
  rds_args=()
  read -r -a rscript_cmd <<< "${rscript}"
  [[ ${#rscript_cmd[@]} -gt 0 ]] || {
    _audit_die "Rscript is required for selected RDS semantic validation"
    return 1
  }
  command -v "${rscript_cmd[0]}" >/dev/null 2>&1 || {
    _audit_die "configured Rscript interpreter is unavailable: ${rscript_cmd[0]}"
    return 1
  }
  [[ -r "${validator}" ]] || {
    _audit_die "benchmark RDS validator is missing: ${validator}"
    return 1
  }
  ecoda_validate_checksum "${path}" || {
    _audit_die "selected RDS checksum is invalid: ${path}"
    return 1
  }
  rds_args=(--artifact "${path}" --method "${label}" --dataset "${dataset}"
    --view "${view}" --input-root "${HPC_SCRATCH_DIR}"
    --config "${DATASETS_JSON_FILE}")
  [[ -n "${AUDIT_METADATA_PASS:-}" ]] &&
    rds_args+=(--batch-pass "${AUDIT_METADATA_PASS}")
  [[ -n "${variant}" ]] &&
    rds_args+=(--analysis-variant "${variant}")
  [[ -n "${identity_path}" ]] &&
    rds_args+=(--expected-batch-contract "${identity_path}")
  [[ "${path}" == *_metadata.rds ]] && rds_args+=(--metadata)
  "${rscript_cmd[@]}" "${validator}" "${rds_args[@]}" >/dev/null 2>&1 || {
    _audit_die "selected RDS contract is invalid: ${path}"
    return 1
  }
}

_audit_semantic_artifact() {
  local path="$1" dataset="$2" view="$3" label="$4" owner_run="$5"
  local identity_path=""
  if [[ "${STAGE_ARG}" == stage5 &&
        "${AUDIT_METADATA_PASS:-}" == corrected &&
        "${label}" != batch_effect_corrected &&
        "${label}" != batch_effect_uncorrected ]]; then
    identity_path="$(
      _audit_batch_contract_identity_path "${dataset}" "${view}" "${label}"
    )" || return 1
  fi
  case "${STAGE_ARG}" in
    stage3)
      case "${path}" in
        *.h5ad)
          ecoda_validate_checksum "${path}" || return 1
          _audit_benchmark_h5ad "${path}" "${view}" "Stage 3 audit" || return 1
          ;;
      esac
      ;;
    stage4)
      case "${path}" in
        *.h5ad)
          ecoda_validate_checksum "${path}" || return 1
          _audit_annotation_artifact "${path}" h5ad || return 1
          ;;
        *.feather)
          ecoda_validate_checksum "${path}" || return 1
          _audit_annotation_artifact "${path}" path || return 1
          ;;
      esac
      ;;
    stage5)
      case "${path}" in
        *.rds)
          _audit_benchmark_rds "${path}" "${dataset}" "${view}" "${label}" \
            "${identity_path}" || return 1
          ;;
        *.feather)
          _audit_matrix_feather "${path}" "${owner_run}" "${identity_path}" || return 1
          ;;
        *.h5ad)
          ecoda_validate_checksum "${path}" || return 1
          _audit_benchmark_h5ad "${path}" "${view}" "${label}" \
            "${identity_path}" || return 1
          ;;
      esac
      ;;
  esac
}
_audit_output_metadata_add() {
  local raw_path="$1" dataset="$2" view="$3" label="$4" canonical
  canonical="$(_ecoda_canonical_path "${raw_path}")" || return 1
  AUDIT_OUTPUT_DATASETS+=("${dataset}")
  AUDIT_OUTPUT_VIEWS+=("${view}")
  AUDIT_OUTPUT_LABELS+=("${label}")
  AUDIT_OUTPUT_CANONICAL_PATHS+=("${canonical}")
}

_audit_prepare_output_semantics() {
  local scope="$1" dataset view label extra name raw_path nas_path
  local artifact_index
  AUDIT_OUTPUT_DATASETS=()
  AUDIT_OUTPUT_VIEWS=()
  AUDIT_OUTPUT_LABELS=()
  AUDIT_OUTPUT_CANONICAL_PATHS=()
  case "${STAGE_ARG}" in
    stage3|stage4)
      while IFS=$'\t' read -r dataset view extra; do
        name="$(ecoda_view_output_name "${dataset}" "${view}")" || return 1
        raw_path="${HPC_SCRATCH_DIR}/${dataset}/output/${name}"
        _audit_output_metadata_add "${raw_path}" "${dataset}" "${view}" "" ||
          return 1
        if [[ -n "${NAS_TARGET_DIR:-}" ]]; then
          nas_path="${NAS_TARGET_DIR}/${dataset}/output/${name}"
          _audit_output_metadata_add "${nas_path}" "${dataset}" "${view}" "" ||
            return 1
        fi
      done < "${scope}"
      ;;
    stage5)
      while IFS=$'\t' read -r dataset view label extra; do
        _ecoda_stage5_artifacts_for "${dataset}" "${view}" "${label}" || return 1
        artifact_index=0
        for raw_path in "${ECODA_BENCHMARK_ARTIFACTS[@]}"; do
          _audit_output_metadata_add "${raw_path}" "${dataset}" "${view}" "${label}" ||
            return 1
          if [[ ${#ECODA_BENCHMARK_ARTIFACT_NAS[@]} -gt 0 ]]; then
            nas_path="${ECODA_BENCHMARK_ARTIFACT_NAS[${artifact_index}]}"
            _audit_output_metadata_add "${nas_path}" "${dataset}" "${view}" "${label}" ||
              return 1
          fi
          artifact_index=$((artifact_index + 1))
        done
      done < "${scope}"
      ;;
  esac
  [[ ${#AUDIT_OUTPUT_CANONICAL_PATHS[@]} -eq ${#ECODA_OUTPUT_PATHS[@]} ]] || {
    _audit_die "selected artifact semantic scope does not match expanded output paths"
    return 1
  }
}

_audit_selected_artifacts() {
  local path owner_dir owner_run metadata_root metadata_pass scope_selection
  local scope_rows
  local semantic_index=0 semantic_dataset="" semantic_view="" semantic_label=""
  local artifact_marker=""
  scope_selection="${SELECTION_ARG}"
  export DATASETS_JSON_FILE="${AUDIT_SOURCE_ROOT}/datasets.json"
  export PROJECT_ROOT="${AUDIT_SOURCE_ROOT}"
  if [[ "${STAGE_ARG}" == stage5 ]]; then
    if [[ "${AUDIT_METADATA_METHOD_MATRIX_MODE:-0}" == 1 ]]; then
      # The dataset selection is the scope source. The run-owned matrix is an
      # independent authorization input consumed for each pending method row.
      export ECODA_STAGE5_METHOD_MATRIX="${AUDIT_METADATA_METHOD_MATRIX}"
      unset METHOD_MATRIX
    else
      unset ECODA_STAGE5_METHOD_MATRIX METHOD_MATRIX
    fi
  fi
  metadata_root="${AUDIT_METADATA_ROOT:-}"
  metadata_pass="${AUDIT_METADATA_PASS:-}"
  if [[ -n "${AUDIT_METADATA_VARIANT:-}" ]]; then
    ecoda_stage5_validate_identity "${metadata_pass}" \
      "${AUDIT_METADATA_VARIANT}" || return 1
  fi
  if [[ "${STAGE_ARG}" == stage5 ]]; then
    case "${metadata_pass}" in
      ""|uncorrected|corrected) ;;
      *) _audit_die "invalid Stage 5 PASS metadata: ${metadata_pass}"; return 1 ;;
    esac
    if [[ -n "${metadata_root}" ]]; then
      [[ "${metadata_root}" = /* ]] || return 1
      export ANALYSIS_ROOT="${metadata_root}"
      if [[ -n "${NAS_TARGET_DIR:-}" &&
            "${metadata_root}" == "${HPC_SCRATCH_DIR}"/* ]]; then
        export ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/${metadata_root#${HPC_SCRATCH_DIR}/}"
      fi
    elif [[ -n "${metadata_pass}" ]]; then
      export ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/${metadata_pass}"
      if [[ -n "${NAS_TARGET_DIR:-}" ]]; then
        export ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/${metadata_pass}"
      fi
    fi
  fi
  if [[ "${STAGE_ARG}" == stage5 &&
        -n "${AUDIT_METADATA_VARIANT:-}" ]]; then
    ecoda_stage5_validate_identity "${metadata_pass}" \
      "${AUDIT_METADATA_VARIANT}" || return 1
  fi
  if [[ "${STAGE_ARG}" == stage5 ]]; then
    _audit_stage5_scope_selection || return 1
    scope_selection="${AUDIT_STAGE5_SCOPE_SELECTION}"
    scope_rows="$(wc -l < "${scope_selection}" | tr -d '[:space:]')"
    [[ "${scope_rows}" == "${AUDIT_MATRIX_PENDING_ROWS:-${scope_rows}}" ]] || {
      _audit_die "expanded Stage 5 scope disagrees with PENDING_METHOD_ROWS"
      return 1
    }
    [[ -s "${scope_selection}" ]] || return 0
  fi
  if [[ -n "${NAS_TARGET_DIR:-}" ]]; then
    export NAS_TARGET_DIR
  fi
  _ecoda_expand_output_selection "${STAGE_ARG}" "${scope_selection}" || {
    _audit_die "selected output expansion failed"
    return 1
  }
  case "${STAGE_ARG}" in
    stage3|stage4|stage5)
      _audit_prepare_output_semantics "${scope_selection}" || return 1
      ;;
  esac
  for path in "${ECODA_OUTPUT_PATHS[@]}"; do
    if [[ "${STAGE_ARG}" == stage3 ||
          "${STAGE_ARG}" == stage4 ||
          "${STAGE_ARG}" == stage5 ]]; then
      [[ "${path}" == "${AUDIT_OUTPUT_CANONICAL_PATHS[${semantic_index}]}" ]] || {
        _audit_die "selected artifact semantic path is not run-bound: ${path}"
        return 1
      }
      semantic_dataset="${AUDIT_OUTPUT_DATASETS[${semantic_index}]}"
      semantic_view="${AUDIT_OUTPUT_VIEWS[${semantic_index}]}"
      semantic_label="${AUDIT_OUTPUT_LABELS[${semantic_index}]}"
      semantic_index=$((semantic_index + 1))
      if [[ "${STAGE_ARG}" == stage5 &&
            -n "${AUDIT_METADATA_VARIANT:-}" ]]; then
        case "${AUDIT_METADATA_VARIANT}" in
          final) artifact_marker="_batch_effect_uncorrected_final_" ;;
          corrected_final) artifact_marker="_batch_effect_corrected_final_" ;;
          *) return 1 ;;
        esac
        [[ "${path##*/}" == *"${artifact_marker}"* ]] || {
          _audit_die "selected Stage 5 artifact stem is not variant-qualified: ${path}"
          return 1
        }
      fi
    else
      semantic_dataset=""
      semantic_view=""
      semantic_label=""
    fi
    owner_dir="$(ecoda_artifact_owner_dir "${path}")" || return 1
    [[ -d "${owner_dir}" ]] || {
      _audit_die "selected artifact owner is missing: ${path}"
      return 1
    }
    _ecoda_artifact_owner_validate_dir "${owner_dir}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || {
      _audit_die "selected artifact owner is not terminal OK: ${path}"
      return 1
    }
    owner_run="${ECODA_ARTIFACT_OWNER_RUN}"
    if [[ -n "${NAS_TARGET_DIR:-}" &&
          "${path}" == "${NAS_TARGET_DIR}"/* ]]; then
      continue
    fi
    _audit_record_for_path "${path}" "${owner_run}" || return 1
    _audit_semantic_artifact "${path}" "${semantic_dataset}" \
      "${semantic_view}" "${semantic_label}" "${owner_run}" || return 1
  done
  if [[ "${STAGE_ARG}" == stage3 ||
        "${STAGE_ARG}" == stage4 ||
        "${STAGE_ARG}" == stage5 ]]; then
    [[ ${semantic_index} -eq ${#AUDIT_OUTPUT_CANONICAL_PATHS[@]} ]] || {
      _audit_die "selected artifact semantic paths were not compared in full"
      return 1
    }
  fi

}

# Validate the supplied source/runtime identities before looking at any
# selected artifact; every path below is derived from these exact arguments or
# this one run root.  No glob, scheduler submission, or repair is performed.
_audit_run_metadata || exit 1
_audit_regular_file "${SOURCE_MANIFEST_ARG}" || exit 1
_audit_regular_file "${RUNTIME_IDENTITY_ARG}" || exit 1
SOURCE_COPY="${RUN_ROOT_REAL}/manifests/source.manifest"
RUNTIME_COPY="${RUN_ROOT_REAL}/manifests/runtime.identity"
_audit_regular_file "${SOURCE_COPY}" || exit 1
_audit_regular_file "${RUNTIME_COPY}" || exit 1
ecoda_validate_run_owned_path "${SOURCE_COPY}" "${RUN_ROOT_REAL}" || exit 1
ecoda_validate_run_owned_path "${RUNTIME_COPY}" "${RUN_ROOT_REAL}" || exit 1
_audit_same_bytes "${SOURCE_MANIFEST_ARG}" "${SOURCE_COPY}" || exit 1
_audit_same_bytes "${RUNTIME_IDENTITY_ARG}" "${RUNTIME_COPY}" || exit 1
_audit_source_manifest "${SOURCE_MANIFEST_ARG}" || exit 1
_audit_batch_contract_manifest || exit 1
_audit_corrected_source_contracts || exit 1
_audit_runtime_identity "${RUNTIME_IDENTITY_ARG}" || exit 1
_audit_selection || exit 1
_audit_terminal_status || exit 1
_audit_stage5_method_matrix || exit 1
_audit_dispatch_selection || exit 1
_audit_metadata_export_manifest || exit 1
_audit_scheduler_ids || exit 1
_audit_selected_artifacts || exit 1
printf 'ECODA_RUN_AUDIT_OK=%s\n' "${RUN_ID}"
