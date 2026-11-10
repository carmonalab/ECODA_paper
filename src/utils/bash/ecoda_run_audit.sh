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
SCHEDULER_ACCOUNTING_ARG=""
STAGE3_REUSE_REPORT_ARG=""
REUSE_SELECTION_ARG=""
TERMINAL_INSPECT_ARG=""
FAILED_GATE_MANIFEST_ARG=""


_audit_die() {
  echo "ERROR: $*" >&2
  return 1
}

_audit_usage() {
  cat <<'EOF'
Usage: ecoda_run_audit.sh --run-root ABSOLUTE_RUN_ROOT --stage STAGE \
       --selection ABSOLUTE_SELECTION \
       --source-manifest ABSOLUTE_SOURCE_MANIFEST \
       --runtime-identity ABSOLUTE_RUNTIME_IDENTITY \
       [--scheduler-accounting ABSOLUTE_RUN_OWNED_ACCOUNTING_TSV]
Stage 3 validator-only subset reuse (preserves the original failed gate):
       --stage3-reuse-report NEW_ABSOLUTE_JSON --reuse-selection ABSOLUTE_SUBSET_TSV
       --terminal-inspect ORIGINAL_INSPECT_JSON --failed-gate-manifest ORIGINAL_GATE_JSON
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
    --scheduler-accounting)
      [[ -n "${2:-}" ]] || { _audit_die "--scheduler-accounting requires a path"; exit 1; }
      SCHEDULER_ACCOUNTING_ARG="$2"
      shift 2
      ;;
    --scheduler-accounting=*)
      SCHEDULER_ACCOUNTING_ARG="${1#*=}"
      [[ -n "${SCHEDULER_ACCOUNTING_ARG}" ]] || { _audit_die "--scheduler-accounting requires a path"; exit 1; }
      shift
      ;;
    --stage3-reuse-report) STAGE3_REUSE_REPORT_ARG="${2:-}"; shift 2 ;;
    --reuse-selection) REUSE_SELECTION_ARG="${2:-}"; shift 2 ;;
    --terminal-inspect) TERMINAL_INSPECT_ARG="${2:-}"; shift 2 ;;
    --failed-gate-manifest) FAILED_GATE_MANIFEST_ARG="${2:-}"; shift 2 ;;
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

if [[ -n "${STAGE3_REUSE_REPORT_ARG}${REUSE_SELECTION_ARG}${TERMINAL_INSPECT_ARG}${FAILED_GATE_MANIFEST_ARG}" ]]; then
  [[ "${STAGE_ARG}" == stage3 && -z "${SCHEDULER_ACCOUNTING_ARG}" &&
     "${STAGE3_REUSE_REPORT_ARG}" = /* && "${REUSE_SELECTION_ARG}" = /* &&
     "${TERMINAL_INSPECT_ARG}" = /* && "${FAILED_GATE_MANIFEST_ARG}" = /* ]] || {
    _audit_die 'Stage 3 reuse requires all four explicit absolute proof/report options'; exit 1;
  }
fi

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

_audit_run_metadata() {
  local metadata="${RUN_ROOT_REAL}/metadata"
  local metadata_stage metadata_run field field_count
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
}

_audit_selection() {
  local columns rows
  case "${STAGE_ARG}" in
    stage2) columns=5 ;;
    stage3|stage4) columns=2 ;;
    *) _audit_die "shared selection audit does not implement stage policy: ${STAGE_ARG}"; return 1 ;;
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

# These derived-data producers have finite operation-specific inventories,
# rather than benchmark family outputs. Their source-pinned read-only validator
# owns that schema; the generic identity and accounting checks still apply.
_audit_analysis_release() {
  local purpose operation expected_operation validator
  local expected_scheduler_kind=ARRAY expected_scalar=1
  local scope_producer=analysis_projection_scope submission_producer=analysis_projection_submission
  local validator_name=project_analysis_release.sh
  local metadata="${RUN_ROOT_REAL}/metadata"
  purpose="$(_audit_metadata_value "${metadata}" PURPOSE)"
  [[ "$(_audit_metadata_count "${metadata}" PURPOSE)" == 1 ]] || return 1
  case "${purpose}" in
    analysis_metadata_projection) expected_operation=projection ;;
    analysis_pca_render) expected_operation=render-pca ;;
    singlecell_breast_unc_pilot)
      expected_operation=singlecell-breast-unc-pilot
      scope_producer=singlecell_pilot_scope
      submission_producer=singlecell_pilot_submission
      validator_name=singlecell_batch_geometry.sh
      ;;
    singlecell_batch_main_bulk)
      expected_operation=singlecell-batch-main-bulk
      scope_producer=singlecell_bulk_scope
      submission_producer=singlecell_bulk_submission
      validator_name=singlecell_batch_geometry_bulk.sh
      expected_scalar=0
      ;;
    saved_cluster_counts_screen)
      expected_operation=saved-clusters
      scope_producer=saved_cluster_scope
      submission_producer=saved_cluster_submission
      validator_name=saved_cluster_screen.sh
      ;;
    counts_only_cell_depth_screen)
      expected_operation=cell-depth-screen
      scope_producer=cell_depth_scope
      submission_producer=cell_depth_submission
      validator_name=cell_depth_screen.sh
      expected_scheduler_kind=worker
      ;;
    counts_only_zeroimp_screen)
      expected_operation=zeroimp-count-screen
      scope_producer=zeroimp_scope
      submission_producer=zeroimp_submission
      validator_name=zeroimp_count_screen.sh
      ;;
    legacy_analysis_compatibility_audit)
      expected_operation=legacy-audit
      scope_producer=legacy_analysis_audit_scope
      submission_producer=legacy_analysis_audit_submission
      validator_name=legacy_analysis_audit.sh
      ;;
    orphaned_main_dependency_observer)
      expected_operation=orphaned-main-observer
      scope_producer=orphaned_main_observer_scope
      submission_producer=orphaned_main_observer_submission
      validator_name=orphaned_main_observer.sh
      ;;
    *) _audit_die "unsupported Stage 5 PURPOSE: ${purpose}"; return 1 ;;
  esac
  [[ "${SELECTION_ARG}" == "${RUN_ROOT_REAL}/manifests/selection.tsv" &&
     "${AUDIT_TERMINAL_STATE}" == OK ]] || return 1
  _audit_regular_file "${RUN_ROOT_REAL}/manifests/operation" || return 1
  operation="$(cat "${RUN_ROOT_REAL}/manifests/operation")"
  [[ "${operation}" == "${expected_operation}" ]] || {
    _audit_die 'analysis PURPOSE and operation disagree'; return 1;
  }
  ecoda_validate_artifact_record "${RUN_ROOT_REAL}/manifests/operation" \
    "${scope_producer}" "${RUN_ID}" >/dev/null || return 1
  ecoda_validate_artifact_record "${SELECTION_ARG}" \
    "${scope_producer}" "${RUN_ID}" >/dev/null || return 1
  ecoda_validate_artifact_record "${RUN_ROOT_REAL}/manifests/scheduler_ids.tsv" \
    "${submission_producer}" "${RUN_ID}" >/dev/null || return 1
  if [[ ${expected_scalar} -eq 1 ]]; then
    [[ "$(wc -l < "${RUN_ROOT_REAL}/manifests/scheduler_ids.tsv" | tr -d '[:space:]')" == 1 ]] || return 1
    awk -F '\t' -v expected="${expected_scheduler_kind}" 'NF != 2 || $1 != expected || $2 !~ /^[0-9]+$/ { exit 1 }' \
      "${RUN_ROOT_REAL}/manifests/scheduler_ids.tsv" || return 1
  else
    # Only the fixed11 purpose owns two initial arrays, one watchdog and retries.
    awk -F '\t' 'NF!=2 || $2!~/^[1-9][0-9]*$/ || seen[$2]++ {exit 1}
      $1=="ARRAY"{a++;next} $1=="WATCHDOG"{w++;next} {exit 1}
      END{if(a<2 || a>6 || w!=1)exit 1}' "${RUN_ROOT_REAL}/manifests/scheduler_ids.tsv" || return 1
  fi
  _audit_scheduler_ids || return 1
  validator="$(ecoda_require_source_script_path \
    "${AUDIT_SOURCE_ROOT}/src/utils/bash/${validator_name}" \
    "${AUDIT_SOURCE_ROOT}")" || return 1
  PROJECT_ROOT="${AUDIT_SOURCE_ROOT}" ECODA_SOURCE_ROOT="${AUDIT_SOURCE_ROOT}" \
    ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST_ARG}" ECODA_SOURCE_SNAPSHOT_REQUIRED=1 \
    bash "${validator}" --validate
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
  if [[ -n "${SCHEDULER_ACCOUNTING_ARG}" ]]; then
    # Durable-gate inspect captures its one real query before artifact audits.
    # Consume those exact bytes; never issue another query as a fallback.
    [[ "${SCHEDULER_ACCOUNTING_ARG}" == "${RUN_ROOT_REAL}/manifests/scheduler_accounting.tsv" ]] || {
      _audit_die "scheduler accounting capture must belong to this exact run"
      return 1
    }
    _audit_regular_file "${SCHEDULER_ACCOUNTING_ARG}" || return 1
    [[ "$(realpath "${SCHEDULER_ACCOUNTING_ARG}")" == "${SCHEDULER_ACCOUNTING_ARG}" ]] || return 1
    _audit_regular_file "${SCHEDULER_ACCOUNTING_ARG}.md5" || return 1
    [[ "$(realpath "${SCHEDULER_ACCOUNTING_ARG}.md5")" == "${SCHEDULER_ACCOUNTING_ARG}.md5" ]] || return 1
    ecoda_validate_checksum "${SCHEDULER_ACCOUNTING_ARG}" || return 1
    rows="$(cat "${SCHEDULER_ACCOUNTING_ARG}")" || return 1
  elif ! rows="$(sacct -n -P -X -j "${id_list}" \
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
      # Gate captures use array-qualified JobID in three columns. Standalone
      # four-column queries also include numeric JobIDRaw. Raw child IDs alone
      # cannot be routed to their array root.
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
  local python_bin validator
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
  "${python_bin}" "${validator}" "${validator_args[@]}" >/dev/null 2>&1 || {
    _audit_die "selected H5AD contract is invalid: ${path}"
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


_audit_output_metadata_add() {
  local raw_path="$1" dataset="$2" view="$3" label="$4" canonical
  canonical="$(_ecoda_canonical_path "${raw_path}")" || return 1
  AUDIT_OUTPUT_DATASETS+=("${dataset}")
  AUDIT_OUTPUT_VIEWS+=("${view}")
  AUDIT_OUTPUT_LABELS+=("${label}")
  AUDIT_OUTPUT_CANONICAL_PATHS+=("${canonical}")
}

_audit_semantic_artifact() {
  local path="$1" dataset="$2" view="$3" input_dir python_bin
  case "${STAGE_ARG}" in
    stage2)
      if [[ -f "${RUN_ROOT_REAL}/manifests/bassez_metadata.json" ]]; then
        ecoda_validate_checksum "${path}" || return 1
        case "${path}" in
          *.h5ad) ecoda_validate_stage2_output bassez_h5ad_metadata "${path}" || return 1 ;;
        esac
      fi
      ;;
    stage3)
      case "${path}" in
        *.h5ad)
          ecoda_validate_checksum "${path}" || return 1
          _audit_benchmark_h5ad "${path}" "${view}" "Stage 3 audit" || return 1
          if [[ -n "${STAGE3_REUSE_REPORT_ARG:-}" ]]; then
            # Reuse validates this same cohort/obs policy in the report helper
            # against the immutable production source references. It does not
            # re-read historical raw expression files to recompute their MD5.
            return 0
          fi
          if [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
            input_dir="$(ecoda_derivative_input "${dataset}" "${view}" directory)" || return 1
          else
            input_dir="${HPC_SCRATCH_DIR}/${dataset}/data"
          fi
          python_bin="$(_audit_python_binary)" || return 1
          "${python_bin}" "${AUDIT_SOURCE_ROOT}/src/utils/py/cohort_selection.py" \
            --config "${AUDIT_SOURCE_ROOT}/datasets.json" --dataset "${dataset}" \
            --view "${view}" --input-dir "${input_dir}" --path "${path}" \
            >/dev/null 2>&1 || return 1
          ;;
      esac
      ;;
    stage4)
      if [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
        ecoda_stage4_parent_check "${dataset}" "${view}" || return 1
        ecoda_stage4_validate_copy_receipt "${dataset}" "${view}" || return 1
      fi
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
  esac
}

_audit_prepare_output_semantics() {
  local scope="$1" dataset view extra name raw_path nas_path
  AUDIT_OUTPUT_DATASETS=()
  AUDIT_OUTPUT_VIEWS=()
  AUDIT_OUTPUT_LABELS=()
  AUDIT_OUTPUT_CANONICAL_PATHS=()
  case "${STAGE_ARG}" in
    stage3|stage4)
      while IFS=$'\t' read -r dataset view extra; do
        name="$(ecoda_view_output_name "${dataset}" "${view}")" || return 1
        if [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
          raw_path="$(ecoda_derivative_dataset_output "${dataset}" "${view}")" || return 1
        else
          raw_path="${HPC_SCRATCH_DIR}/${dataset}/output/${name}"
        fi
        _audit_output_metadata_add "${raw_path}" "${dataset}" "${view}" "" ||
          return 1
        if [[ -n "${NAS_TARGET_DIR:-}" ]]; then
          if [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
            nas_path="$(ecoda_derivative_dataset_output "${dataset}" "${view}" nas)" || return 1
          else
            nas_path="${NAS_TARGET_DIR}/${dataset}/output/${name}"
          fi
          _audit_output_metadata_add "${nas_path}" "${dataset}" "${view}" "" ||
            return 1
        fi
      done < "${scope}"
      ;;
    *) return 1 ;;
  esac
  [[ ${#AUDIT_OUTPUT_CANONICAL_PATHS[@]} -eq ${#ECODA_OUTPUT_PATHS[@]} ]] || {
    _audit_die "selected artifact semantic scope does not match expanded output paths"
    return 1
  }
}

_audit_bassez_nas_artifact() {
  local path="$1" directory
  ecoda_validate_artifact_record "${path}" bassez_h5ad_metadata "${RUN_ID}" || {
    _audit_die "Bassez NAS producer record is missing or invalid: ${path}"
    return 1
  }
  ecoda_validate_stage2_output bassez_h5ad_metadata "${path}" || {
    _audit_die "Bassez NAS publication validation failed: ${path}"
    return 1
  }
  [[ "$(_ecoda_artifact_mode_bits "${path}")" == 666 &&
     "$(_ecoda_artifact_mode_bits "${path}.md5")" == 666 ]] || {
    _audit_die "Bassez NAS files must have mode 666: ${path}"
    return 1
  }
  directory="$(dirname "${path}")"
  while [[ "${directory}" != "${NAS_TARGET_DIR%/}" ]]; do
    [[ "${directory}" == "${NAS_TARGET_DIR%/}/"* &&
       "$(_ecoda_artifact_mode_bits "${directory}")" == 777 ]] || {
      _audit_die "Bassez selected NAS directories must have mode 777: ${directory}"
      return 1
    }
    directory="$(dirname "${directory}")"
  done
}

_audit_selected_artifacts() {
  local path owner_dir owner_run scope_selection
  local semantic_index=0 semantic_dataset="" semantic_view="" semantic_label=""
  scope_selection="${SELECTION_ARG}"
  export DATASETS_JSON_FILE="${AUDIT_SOURCE_ROOT}/datasets.json"
  export PROJECT_ROOT="${AUDIT_SOURCE_ROOT}"
  _ecoda_expand_output_selection "${STAGE_ARG}" "${scope_selection}" || {
    _audit_die "selected output expansion failed"
    return 1
  }
  case "${STAGE_ARG}" in
    stage3|stage4)
      _audit_prepare_output_semantics "${scope_selection}" || return 1
      ;;
  esac
  for path in "${ECODA_OUTPUT_PATHS[@]}"; do
    if [[ "${STAGE_ARG}" == stage3 ||
          "${STAGE_ARG}" == stage4 ]]; then
      [[ "${path}" == "${AUDIT_OUTPUT_CANONICAL_PATHS[${semantic_index}]}" ]] || {
        _audit_die "selected artifact semantic path is not run-bound: ${path}"
        return 1
      }
      semantic_dataset="${AUDIT_OUTPUT_DATASETS[${semantic_index}]}"
      semantic_view="${AUDIT_OUTPUT_VIEWS[${semantic_index}]}"
      semantic_label="${AUDIT_OUTPUT_LABELS[${semantic_index}]}"
      semantic_index=$((semantic_index + 1))
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
      if [[ "${STAGE_ARG}" == stage4 && -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
        ecoda_stage4_validate_nas "${semantic_dataset}" "${semantic_view}" || return 1
      fi
      if [[ "${STAGE_ARG}" == stage2 &&
            -f "${RUN_ROOT_REAL}/manifests/bassez_metadata.json" ]]; then
        _audit_bassez_nas_artifact "${path}" || return 1
      fi
      continue
    fi
    _audit_record_for_path "${path}" "${owner_run}" || return 1
    _audit_semantic_artifact "${path}" "${semantic_dataset}" \
      "${semantic_view}" "${semantic_label}" || return 1
  done
  if [[ "${STAGE_ARG}" == stage3 ||
        "${STAGE_ARG}" == stage4 ]]; then
    [[ ${semantic_index} -eq ${#AUDIT_OUTPUT_CANONICAL_PATHS[@]} ]] || {
      _audit_die "selected artifact semantic paths were not compared in full"
      return 1
    }
  fi
}



# Validate the supplied source/runtime identities before looking at any
# selected artifact; every path below is derived from these exact arguments or
# this one run root.  No glob, scheduler submission, or repair is performed.
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
export ECODA_RUN_ROOT="${RUN_ROOT_REAL}" ECODA_RUN_ID="${RUN_ID}"
if [[ "${STAGE_ARG}" == stage2 &&
      ( -e "${RUN_ROOT_REAL}/manifests/bassez_metadata.json" ||
        -e "${RUN_ROOT_REAL}/manifests/bassez_metadata.input.tsv" ) ]]; then
  export PROJECT_ROOT="${AUDIT_SOURCE_ROOT}" ECODA_SOURCE_ROOT="${AUDIT_SOURCE_ROOT}"
  export DATASETS_JSON_FILE="${AUDIT_SOURCE_ROOT}/datasets.json"
  export PYTHON_BIN="$(_audit_python_binary)"
  ecoda_load_stage2_bassez_metadata_helper || exit 1
fi
if [[ -f "${RUN_ROOT_REAL}/manifests/artifact_bindings.json" ||
      -n "${ECODA_ARTIFACT_BINDING:-}" ]]; then
  expected_binding="${RUN_ROOT_REAL}/manifests/artifact_bindings.json"
  [[ -z "${ECODA_ARTIFACT_BINDING:-}" ||
     "${ECODA_ARTIFACT_BINDING}" == "${expected_binding}" ]] || exit 1
  export ECODA_ARTIFACT_BINDING="${expected_binding}"
  export PROJECT_ROOT="${AUDIT_SOURCE_ROOT}"
  export DATASETS_JSON_FILE="${AUDIT_SOURCE_ROOT}/datasets.json"
  export PYTHON_BIN="$(_audit_python_binary)"
  source "${AUDIT_SOURCE_ROOT}/src/utils/bash/ecoda_derivative_artifacts.sh"
  ecoda_derivative_load "${STAGE_ARG}" || exit 1
fi
if [[ "${STAGE_ARG}" == stage5 ]]; then
  STAGE5_POLICY_SCRIPT="$(
    ecoda_require_source_script_path \
      "${AUDIT_SOURCE_ROOT}/src/utils/bash/ecoda_stage5_policy.sh" \
      "${AUDIT_SOURCE_ROOT}"
  )" || exit 1
  source "${STAGE5_POLICY_SCRIPT}"
  STAGE5_AUDIT_POLICY_SCRIPT="$(
    ecoda_require_source_script_path \
      "${AUDIT_SOURCE_ROOT}/src/utils/bash/ecoda_run_audit_stage5.sh" \
      "${AUDIT_SOURCE_ROOT}"
  )" || exit 1
  source "${STAGE5_AUDIT_POLICY_SCRIPT}"
fi
_audit_run_metadata || exit 1
_audit_runtime_identity "${RUNTIME_IDENTITY_ARG}" || exit 1
if [[ "${STAGE_ARG}" == stage5 &&
      "$(_audit_metadata_count "${RUN_ROOT_REAL}/metadata" PURPOSE)" == 1 &&
      "$(_audit_metadata_value "${RUN_ROOT_REAL}/metadata" PURPOSE)" == singlecell_batch_main_bulk ]]; then
  # The producer's frozen Stage5 policy predates this operation's four-column
  # row/dataset/representation/memory schema. Keep that original source pinned;
  # only this exact purpose uses its own fixed11 selection before delegation.
  [[ "${SELECTION_ARG}" == "${RUN_ROOT_REAL}/manifests/selection.tsv" ]] || exit 1
  ecoda_validate_run_owned_path "${SELECTION_ARG}" "${RUN_ROOT_REAL}" || exit 1
  ecoda_validate_manifest "${SELECTION_ARG}" 4 || exit 1
  ecoda_validate_checksum "${SELECTION_ARG}" || exit 1
  AUDIT_SELECTION_ROWS="$(wc -l < "${SELECTION_ARG}" | tr -d '[:space:]')"
  [[ "${AUDIT_SELECTION_ROWS}" == 11 ]] || { _audit_die 'singlecell bulk selection requires exactly11 rows'; exit 1; }
else
  _audit_selection || exit 1
fi
_audit_terminal_status || exit 1
if [[ "${STAGE_ARG}" == stage5 &&
      ( "$(_audit_metadata_count "${RUN_ROOT_REAL}/metadata" PURPOSE)" != 0 ||
        -e "${RUN_ROOT_REAL}/manifests/operation" ||
        -L "${RUN_ROOT_REAL}/manifests/operation" ) ]]; then
  _audit_analysis_release || exit 1
  printf 'ECODA_RUN_AUDIT_OK=%s\n' "${RUN_ID}"
  exit 0
fi
if [[ "${STAGE_ARG}" == stage5 ]]; then
  _audit_stage5_method_matrix || exit 1
  _audit_dispatch_selection || exit 1
  _audit_metadata_export_manifest || exit 1
  _audit_batch_contract_manifest || exit 1
  _audit_corrected_source_contracts || exit 1
fi
if [[ -n "${STAGE3_REUSE_REPORT_ARG}" ]]; then
  # Audit the original full-run identities above. Expensive source/H5AD/NAS
  # checks use only this exact subset; the complete failed accounting is kept.
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" ]] || {
    _audit_die 'Stage 3 subset reuse requires an immutable derivative binding'; exit 1;
  }
  reuse_helper="$(realpath "${SCRIPT_DIR}/../py/stage3_reuse_report.py")" || exit 1
  _audit_regular_file "${reuse_helper}" || exit 1
  reuse_python="$(_audit_python_binary)" || exit 1
  reuse_args=(--run-root "${RUN_ROOT_REAL}" --selection "${SELECTION_ARG}"
    --reuse-selection "${REUSE_SELECTION_ARG}" --source-manifest "${SOURCE_MANIFEST_ARG}"
    --runtime-identity "${RUNTIME_IDENTITY_ARG}" --terminal-inspect "${TERMINAL_INSPECT_ARG}"
    --failed-gate-manifest "${FAILED_GATE_MANIFEST_ARG}")
  "${reuse_python}" "${reuse_helper}" "${reuse_args[@]}" --check-history || exit 1
  SELECTION_ARG="${REUSE_SELECTION_ARG}"
  _audit_selected_artifacts || exit 1
  "${reuse_python}" "${reuse_helper}" "${reuse_args[@]}" \
    --stage3-reuse-report "${STAGE3_REUSE_REPORT_ARG}" --audit-entrypoint "${SCRIPT_DIR}/ecoda_run_audit.sh" || exit 1
  printf 'ECODA_STAGE3_NOOP_VALIDATED=%s\n' "${STAGE3_REUSE_REPORT_ARG}"
  exit 0
fi
_audit_scheduler_ids || exit 1
_audit_selected_artifacts || exit 1
printf 'ECODA_RUN_AUDIT_OK=%s\n' "${RUN_ID}"
