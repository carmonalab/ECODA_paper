#!/bin/bash
#
# Shared helper functions for the benchmark SLURM submitters
# (src/5_run_benchmark_methods/run_python_sample_embedding_methods/,
#  run_r_sample_embedding_methods/, run_transformation_zeroimp_analysis/).
#
# Source AFTER slurm_config.sh and `cd "${PROJECT_ROOT}"` (functions use the
# exported slurm_config vars at call time):
#   source "$(dirname "${BASH_SOURCE[0]}")/benchmark_submit_common.sh"
#
# Provided:
#   benchmark_resolve_datasets <ds_name_arg>
#       Fills the global DATASET_NAMES/NUM_DATASETS from datasets.json.
#   benchmark_validate_h5ad_inputs <view> <method> <datasets...>
#       Validates authoritative processed H5AD inputs before worker submission.
#   benchmark_sync_artifacts_for <dataset> <label>
#       Resolves the guarded artifact set for one dataset and method.
#   benchmark_validate_runtime_metadata <metadata> <output> <md5>
#       Validates the runtime metadata bound to one output artifact.
#   benchmark_merge_sync_cleanup <labels...>
#       Merges execution logs, verifies checksums, synchronizes selected
#       artifacts, and cleans up the run-owned execution logs.
#
# The canonical Stage 5 launcher (`stage5_dispatcher.sh`) owns matrix
# submission, `matrix_watchdog.sh` owns per-method retries and task gates,
# and `matrix_gate.sh` owns the aggregate gate. Compatibility entrypoints
# delegate to that launcher.
# ============================================================================

# Path to the shared exec-log merge script, resolved from THIS file's location
# (BASH_SOURCE[0] inside a sourced file is the sourced file's path).
ANALYSIS_MERGE_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py"

if ! command -v ecoda_stage5_validate_identity >/dev/null 2>&1; then
  source "$(dirname "${BASH_SOURCE[0]}")/../utils/bash/ecoda_stage5_policy.sh"
fi
benchmark_stage5_identity_guard() {
  local requested_pass="${1:-${PASS_ARG:-${ANALYSIS_PASS:-}}}"
  local scratch_root="${HPC_SCRATCH_DIR:-}" nas_root="${NAS_TARGET_DIR:-}"
  local expected_root expected_nas expected_pass expected_suffix
  local variant="${ANALYSIS_VARIANT:-}"
  if command -v ecoda_stage5_validate_identity >/dev/null 2>&1; then
    ecoda_stage5_validate_identity "${requested_pass}"
    return $?
  fi
  case "${variant}" in
    "")
      return 0
      ;;
    final)
      expected_pass="uncorrected"
      expected_suffix="uncorrected_final"
      ;;
    corrected_final)
      expected_pass="corrected"
      expected_suffix="corrected_final"
      if [[ "${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION:-}" == "recovery_35row" ||
            "${ANALYSIS_ROOT:-}" == */batch_effect/corrected_final/recovery_35row ||
            "${ANALYSIS_NAS_ROOT:-}" == */batch_effect/corrected_final/recovery_35row ]]; then
        expected_suffix="corrected_final/recovery_35row"
      fi
      ;;
    *)
      echo "ERROR: unsupported Stage 5 analysis variant: ${variant}" >&2
      return 1
      ;;
  esac
  [[ "${requested_pass}" == "${expected_pass}" ]] || {
    echo "ERROR: ${variant} Stage 5 pass identity disagrees with ${expected_pass}." >&2
    return 1
  }
  [[ -z "${PASS_ARG:-}" || "${PASS_ARG}" == "${requested_pass}" ]] || {
    echo "ERROR: final Stage 5 PASS_ARG disagrees with pass identity." >&2
    return 1
  }
  [[ -z "${ANALYSIS_PASS:-}" || "${ANALYSIS_PASS}" == "${requested_pass}" ]] || {
    echo "ERROR: final Stage 5 ANALYSIS_PASS disagrees with pass identity." >&2
    return 1
  }
  [[ "${scratch_root}" = /* && "${nas_root}" = /* ]] || {
    echo "ERROR: final Stage 5 identity is incomplete or not absolute." >&2
    return 1
  }
  scratch_root="${scratch_root%/}"
  nas_root="${nas_root%/}"
  [[ -n "${scratch_root}" ]] || scratch_root="/"
  [[ -n "${nas_root}" ]] || nas_root="/"
  if [[ "${scratch_root}" == "/" ]]; then
    expected_root="/batch_effect/${expected_suffix}"
  else
    expected_root="${scratch_root}/batch_effect/${expected_suffix}"
  fi
  if [[ "${nas_root}" == "/" ]]; then
    expected_nas="/batch_effect/${expected_suffix}"
  else
    expected_nas="${nas_root}/batch_effect/${expected_suffix}"
  fi
  [[ "${ANALYSIS_ROOT:-}" == "${expected_root}" &&
     "${ANALYSIS_NAS_ROOT:-}" == "${expected_nas}" ]] || {
    echo "ERROR: ${variant} Stage 5 roots are not bound to ${expected_suffix}." >&2
    return 1
  }
}

benchmark_stage5_method_guard() {
  local method="${1:-}"
  if command -v ecoda_stage5_validate_final_method >/dev/null 2>&1; then
    ecoda_stage5_validate_final_method "${method}"
    return $?
  fi
  [[ -n "${ANALYSIS_VARIANT:-}" ]] || return 0
  case "${method}" in
    prepare_pseudobulk|pseudobulk|gloscope|composition|mrvi|pilot|qot) ;;
    *) echo "ERROR: method is not permitted in the final Stage 5 suite: ${method}" >&2; return 1 ;;
  esac
}

# Per-array job duration records for optional synchronization reports.
JOB_REPORTS=()

# Sync-status email helper (best-effort; requires USER_EMAIL from slurm_config.sh).
source "$(dirname "${BASH_SOURCE[0]}")/../utils/bash/sync_status_email.sh"

benchmark_source_script_path() {
  local relative="$1"
  local root="${ECODA_SOURCE_ROOT:-${PROJECT_ROOT}}"
  printf '%s/%s' "${root%/}" "${relative}"
}

benchmark_require_source_script_path() {
  local candidate="$1"
  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    [[ -n "${ECODA_SOURCE_ROOT:-}" ]] || return 1
    ecoda_require_source_script_path "${candidate}" "${ECODA_SOURCE_ROOT}"
  else
    [[ -f "${candidate}" && -r "${candidate}" && ! -L "${candidate}" ]]
  fi
}

benchmark_validate_bound_runtime() {
  benchmark_stage5_identity_guard || return 1
  if [[ -n "${ECODA_RUN_ROOT:-}" ||
        "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    ecoda_runtime_validate_bound_run
  else
    ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE:-host}"
  fi
}

benchmark_validate_output_scope() {
  local selection="${ECODA_SELECTION_MANIFEST:-${ECODA_RUN_ROOT:-}/manifests/selection.tsv}"
  local ownership_selection="${selection}" ownership_tmp="" rc
  benchmark_stage5_identity_guard || return 1
  [[ -n "${ECODA_RUN_ID:-}" ]] || return 0
  [[ -r "${selection}" ]] || return 1
  command -v ecoda_stage5_validate_output_ownership >/dev/null 2>&1 || return 1
  if awk -F '\t' '
      NF != 3 { saw_extended=1 }
      END { exit(saw_extended ? 0 : 1) }
    ' "${selection}"; then
    ownership_tmp="${selection}.ownership.$$"
    ownership_selection="${ownership_tmp}"
  fi
  ecoda_stage5_validate_output_ownership "${ownership_selection}" "${ECODA_RUN_ID}"
  rc=$?
  [[ -n "${ownership_tmp}" ]] && rm -f "${ownership_tmp}"
  return "${rc}"
}


# ---------------------------------------------------------------------------
# Dataset resolution (see header)
# ---------------------------------------------------------------------------
benchmark_resolve_datasets() {
  local DS_NAME_ARG="$1"
  DATASET_NAMES=()
  if [[ -n "${DS_NAME_ARG}" ]]; then
    if ! jq -e --arg ds "${DS_NAME_ARG}" 'has($ds)' "${DATASETS_JSON_FILE}" > /dev/null 2>&1; then
      echo "ERROR: '${DS_NAME_ARG}' is not a dataset in ${DATASETS_JSON_FILE}."
      exit 1
    fi
    DATASET_NAMES+=("${DS_NAME_ARG}")
  else
    while IFS= read -r name; do
      DATASET_NAMES+=("$name")
    done < <(jq -r 'to_entries[] |
      select(.value.use_for_benchmark == true) |
      select(.value.views.benchmark_analysis != null) |
      .key | select(startswith("_") | not)' "${DATASETS_JSON_FILE}")
  fi

  NUM_DATASETS=${#DATASET_NAMES[@]}
  if [[ ${NUM_DATASETS} -eq 0 ]]; then
    echo "ERROR: No benchmark datasets found in ${DATASETS_JSON_FILE}."
    exit 1
  fi

  echo "Found ${NUM_DATASETS} benchmark datasets."
}

# Validate authoritative processed h5ads before submitting any benchmark
# workers. This is intentionally independent of result-file existence: a
# metadata/PCA-only local mirror must never be accepted as a worker input.
benchmark_validate_h5ad_inputs() {
  local VIEW="$1"
  local METHOD="$2"
  shift 2
  local DS_NAME OUTPUT_FILE H5AD_PATH
  for DS_NAME in "$@"; do
    OUTPUT_FILE="$(jq -r --arg ds "${DS_NAME}" --arg view "${VIEW}" \
      '.[$ds].views[$view].output_file_name // empty' \
      "${DATASETS_JSON_FILE}")"
    if [[ -z "${OUTPUT_FILE}" ]]; then
      echo "ERROR: no output_file_name for ${DS_NAME}/${VIEW} in ${DATASETS_JSON_FILE}."
      return 1
    fi
    H5AD_PATH="${HPC_SCRATCH_DIR}/${DS_NAME}/output/${OUTPUT_FILE}"
    echo "Validating h5ad contract: ${H5AD_PATH}" >&2
    "${PYTHON_BIN}" \
      "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
      --path "${H5AD_PATH}" \
      --view "${VIEW}" \
      --method "${METHOD}" || return 1
  done
}


# ---------------------------------------------------------------------------
# "Job durations" block for the final emails: label, job id, array wall time
# for every array gated during this run (JOB_REPORTS, in submission order).
# ---------------------------------------------------------------------------
benchmark_job_durations_block() {
  if (( ${#JOB_REPORTS[@]} == 0 )); then
    printf '%s' "Job durations: n/a (no gated arrays)."
    return 0
  fi
  printf '%s\n' "Job durations (label, job id, array wall time):"
  local line label jid wall
  for line in "${JOB_REPORTS[@]}"; do
    IFS='|' read -r label jid wall <<< "${line}"
    printf '  %-30s %s  %s\n' "${label}" "${jid}" "${wall}"
  done
}

# ---------------------------------------------------------------------------
# NAS check -> RDS integrity sidecar -> merge exec logs -> rsync -> cleanup
# ---------------------------------------------------------------------------

benchmark_sync_artifacts_for() {
  local ds="$1" label="$2" pass="${ANALYSIS_PASS:-${PASS_ARG:-}}"
  local root view saved_root had_root=0 runtime_count runtime_index
  SYNC_ARTIFACTS=()
  benchmark_stage5_identity_guard "${pass}" || return 1
  benchmark_stage5_method_guard "${label}" || return 1
  if [[ -n "${ECODA_STAGE5_METHOD_MATRIX:-${METHOD_MATRIX:-}}" ]] &&
     ! ecoda_stage5_method_matrix_allows "${ds}" \
       "batch_effect_${pass}" "${label}"; then
    # A declared dataset row can have a narrower method scope than the
    # global fixed suite.  Treat unauthorized methods as absent from sync,
    # never reconstruct or claim their paths.
    SYNC_ARTIFACTS=()
    return 0
  fi
  root="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
  if [[ -n "${pass}" && -z "${ANALYSIS_ROOT:-}" ]]; then
    if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
      root="${HPC_SCRATCH_DIR}/batch_effect/$(ecoda_stage5_analysis_root_suffix)"
    else
      root="${HPC_SCRATCH_DIR}/batch_effect/${pass}"
    fi
  fi
  case "${pass}" in
    uncorrected|corrected) view="batch_effect_${pass}" ;;
    *) view="benchmark_analysis" ;;
  esac
  if [[ -n "${ANALYSIS_ROOT+x}" ]]; then
    had_root=1
    saved_root="${ANALYSIS_ROOT}"
  fi
  ANALYSIS_ROOT="${root}"
  if ! _ecoda_stage5_artifacts_for "${ds}" "${view}" "${label}"; then
    if [[ ${had_root} -eq 1 ]]; then
      ANALYSIS_ROOT="${saved_root}"
    else
      unset ANALYSIS_ROOT
    fi
    return 1
  fi
  SYNC_ARTIFACTS=("${ECODA_BENCHMARK_ARTIFACTS[@]}")
  if [[ ${had_root} -eq 1 ]]; then
    ANALYSIS_ROOT="${saved_root}"
  else
    unset ANALYSIS_ROOT
  fi
  case "${label}" in
    mrvi|scpoli|pilot|qot|pilotgm|trans|zeroimp)
      # Keep the payload selection semantics above unchanged, then require its
      # runtime record as an artifact in the same guarded selection. The
      # generic add_sync_artifact caller validates each JSON and automatically
      # carries the JSON's .md5 sidecar into the sync manifests.
      runtime_count=${#SYNC_ARTIFACTS[@]}
      runtime_index=0
      while [[ ${runtime_index} -lt ${runtime_count} ]]; do
        SYNC_ARTIFACTS+=("${SYNC_ARTIFACTS[${runtime_index}]}.runtime.json")
        runtime_index=$((runtime_index + 1))
      done
      ;;
  esac
}
# Validate the runtime metadata bound to a selected output. The JSON parser is
# deliberately standard-library-only so this guard remains usable with the
# immutable Python interpreter configured by the submitter.
benchmark_validate_runtime_metadata() {
  local metadata_path="$1" output_path="$2" output_md5="$3"
  benchmark_stage5_identity_guard || return 1
  if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
    case "${output_path}" in
      "${ANALYSIS_ROOT}"/*|"${ANALYSIS_NAS_ROOT}"/*) ;;
      *) echo "ERROR: ${ANALYSIS_VARIANT} runtime output is outside its analysis roots." >&2; return 1 ;;
    esac
  fi
  "${PYTHON_BIN}" - "${metadata_path}" "${output_path}" "${output_md5}" <<'PY'
import json
import math
import re
import sys


def reject_duplicate_keys(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON field")
        value[key] = item
    return value


def reject_nonfinite_json_constant(value):
    raise ValueError("non-finite JSON constant: {}".format(value))


def valid_nonnegative_number(value, allow_none=False):
    if value is None and allow_none:
        return True
    if type(value) not in (int, float):
        return False
    try:
        number = float(value)
    except (OverflowError, TypeError, ValueError):
        return False
    return math.isfinite(number) and number >= 0


metadata_path, output_path, output_md5 = sys.argv[1:]
required = {
    "schema_version",
    "artifact_path",
    "artifact_md5",
    "dataset",
    "method",
    "time_secs",
    "mem_GB",
}
try:
    with open(metadata_path, "r", encoding="utf-8") as handle:
        payload = json.load(
            handle,
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonfinite_json_constant,
        )
    if not isinstance(payload, dict) or set(payload) != required:
        raise ValueError("runtime metadata must have exactly seven fields")
    if type(payload["schema_version"]) is not int or payload["schema_version"] != 1:
        raise ValueError("invalid schema_version")
    if (
        type(payload["artifact_path"]) is not str
        or payload["artifact_path"] != output_path
    ):
        raise ValueError("runtime metadata artifact_path does not match output")
    artifact_md5 = payload["artifact_md5"]
    if (
        type(artifact_md5) is not str
        or re.fullmatch(r"[0-9a-f]{32}", artifact_md5) is None
        or artifact_md5 != output_md5
    ):
        raise ValueError("runtime metadata artifact_md5 does not match output")
    for field in ("dataset", "method"):
        value = payload[field]
        if type(value) is not str or not value.strip():
            raise ValueError("runtime metadata {} is blank or not a string".format(field))
    if not valid_nonnegative_number(payload["time_secs"]):
        raise ValueError("runtime metadata time_secs is invalid")
    if not valid_nonnegative_number(payload["mem_GB"], allow_none=True):
        raise ValueError("runtime metadata mem_GB is invalid")
except (OSError, UnicodeError, ValueError, TypeError, json.JSONDecodeError):
    sys.exit(1)
PY
}

analysis_merge_sync_cleanup() (
  local LABELS=("$@")
  local STAGE5_PASS="${ANALYSIS_PASS:-${PASS_ARG:-}}"
  local LOCAL_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
  local REMOTE_ROOT="${ANALYSIS_NAS_ROOT:-${NAS_TARGET_DIR}/benchmark}"
  local ds view row_label label path rel line selected_seen=""
  local metadata_manifest metadata_ds metadata_view metadata_input metadata_output metadata_extra expected_metadata_view
  benchmark_stage5_identity_guard "${STAGE5_PASS}" || exit 1
  if [[ -z "${ANALYSIS_ROOT:-}" && -n "${STAGE5_PASS}" ]]; then
    if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
      LOCAL_ROOT="${HPC_SCRATCH_DIR}/batch_effect/$(ecoda_stage5_analysis_root_suffix)"
    else
      LOCAL_ROOT="${HPC_SCRATCH_DIR}/batch_effect/${STAGE5_PASS}"
    fi
  fi
  if [[ -z "${ANALYSIS_NAS_ROOT:-}" && -n "${STAGE5_PASS}" ]]; then
    if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
      REMOTE_ROOT="${NAS_TARGET_DIR}/batch_effect/$(ecoda_stage5_analysis_root_suffix)"
    else
      REMOTE_ROOT="${NAS_TARGET_DIR}/batch_effect/${STAGE5_PASS}"
    fi
  fi
  local LOG_PREFIX="${ANALYSIS_LOG_PREFIX:-execution_times_}"
  local RUN_ROOT="${ECODA_RUN_ROOT:-}"
  local RUN_ID="${ECODA_RUN_ID:-}"
  local RUN_LOG_DIR="${EXECUTION_LOG_DIR:-${RUN_ROOT}/logs}"
  local KIND="benchmark" KIND_CAP="Benchmark"
  local SYNC_OWNER_DIR="" SYNC_LOCK_DIR="" SYNC_FINAL_STATE="FAIL"
  local sync_owner_artifact_valid=0
  local recovery_sync_owner_dir="" recovery_sync_owner_state=""
  local SYNC_FILES SYNC_FILES_TMP NO_CHECKSUM_FILES_TMP CLEANUP_MANIFEST
  local CHECKSUM_TMP REMOTE_CHECKSUM_TMP EXISTING_LOG
  local merge_script="${ANALYSIS_MERGE_SCRIPT}"
  local artifact selected_index
  # Worker artifacts are terminal and immutable during this sync; retain each
  # strict validation digest for the local checksums manifest instead of
  # hashing every selected payload a second time.
  local -a SELECTED_RELS=() SELECTED_DIGESTS=()
  sync_fail() {
    local reason="$1"
    SYNC_FINAL_STATE="FAIL"
    echo "ERROR: ${reason}" >&2
    notify_sync_status \
      "ECODA: ${KIND} NOT synced" \
      "${KIND_CAP} sync to NAS skipped (datasets: ${DATASET_NAMES[*]}, labels: ${LABELS[*]}): ${reason}" || true
    exit 1
  }
  sync_cleanup() {
    local rc="$?"
    if [[ -n "${SYNC_OWNER_DIR}" ]]; then
      if ! ecoda_owner_set_state "${SYNC_OWNER_DIR}" "${SYNC_FINAL_STATE}" \
          "Stage 5 synchronization terminal state"; then
        rc=1
      fi
    fi
    if [[ -n "${SYNC_LOCK_DIR}" ]]; then
      rmdir "${SYNC_LOCK_DIR}" 2>/dev/null || rc=1
    fi
    exit "${rc}"
  }
  trap sync_cleanup EXIT

  [[ -n "${RUN_ROOT}" && -n "${RUN_ID}" ]] || sync_fail "Stage 5 run root/id is not set"
  if [[ "${LOCAL_ROOT}" != "${HPC_SCRATCH_DIR}/benchmark" ]]; then
    KIND="analysis"
    KIND_CAP="Analysis"
  fi
  echo "Checking NAS reachability before any destructive merge work..."
  [[ -d "${NAS_TARGET_DIR}" ]] || sync_fail "NAS path ${REMOTE_ROOT} is unreachable"
  if ! mkdir -p "${REMOTE_ROOT}" "${LOCAL_ROOT}/embeddings" "${RUN_ROOT}/manifests"; then
    sync_fail "failed to create benchmark sync directories"
  fi
  LOCAL_ROOT="$(cd "${LOCAL_ROOT}" && pwd)" ||
    sync_fail "failed to canonicalize local benchmark root"
  REMOTE_ROOT="$(cd "${REMOTE_ROOT}" && pwd)" ||
    sync_fail "failed to canonicalize remote benchmark root"
  RUN_ROOT="$(cd "${RUN_ROOT}" && pwd)" ||
    sync_fail "failed to canonicalize Stage 5 run root"
  RUN_LOG_DIR="$(cd "${RUN_LOG_DIR}" && pwd)" ||
    sync_fail "failed to canonicalize Stage 5 run log directory"
  [[ -n "${ECODA_SOURCE_ROOT:-}" ]] &&
    merge_script="$(benchmark_source_script_path src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py)"
  ECODA_RUN_ROOT="${RUN_ROOT}"
  # The replacement corrected-final lane is a terminal publication root.
  # Reject a successful owner without touching it, but preserve the
  # established terminal-failure recovery path. Direct legacy roots
  # intentionally retain their historical reentry behavior.
  if [[ "${LOCAL_ROOT}" == */batch_effect/corrected_final/recovery_35row &&
        "${REMOTE_ROOT}" == */batch_effect/corrected_final/recovery_35row ]]; then
    sync_owner_artifact_valid=1
    recovery_sync_owner_dir="$(ecoda_owner_dir stage5 "sync/${LOCAL_ROOT}")" ||
      sync_fail "cannot resolve shared Stage 5 sync owner"
    if [[ -e "${recovery_sync_owner_dir}" || -L "${recovery_sync_owner_dir}" ]]; then
      recovery_sync_owner_state="$(
        ecoda_owner_state "${recovery_sync_owner_dir}" 2>/dev/null || true
      )"
      case "${recovery_sync_owner_state}" in
        FAIL)
          # A failed synchronization may be reclaimed by a later run.
          sync_owner_artifact_valid=0
          ;;
        OK)
          sync_fail "shared Stage 5 sync owner is already finalized successfully"
          ;;
        ACTIVE)
          sync_fail "shared Stage 5 sync owner is active"
          ;;
        *)
          sync_fail "shared Stage 5 sync owner has an invalid state"
          ;;
      esac
    fi
  fi
  SYNC_OWNER_DIR="$(ecoda_owner_acquire stage5 "sync/${LOCAL_ROOT}" "${RUN_ID}" 0 \
    "${sync_owner_artifact_valid}")" || {
    sync_fail "shared Stage 5 sync owner is unavailable"
  }
  SYNC_LOCK_DIR="${RUN_ROOT}/sync.lock"
  mkdir "${SYNC_LOCK_DIR}" 2>/dev/null || sync_fail "run sync lock already exists: ${SYNC_LOCK_DIR}"

  SYNC_FILES="${RUN_ROOT}/manifests/sync_files.tsv"
  SYNC_FILES_TMP="${SYNC_FILES}.build.$$"
  NO_CHECKSUM_FILES_TMP="${RUN_ROOT}/manifests/sync_files.no_checksum.build.$$"
  CLEANUP_MANIFEST="${RUN_ROOT}/manifests/cleanup.tsv"
  if ! : > "${SYNC_FILES_TMP}"; then
    sync_fail "failed to create selected sync manifest"
  fi
  SELECTED_RELS=()
  SELECTED_DIGESTS=()

  add_sync_artifact() {
    local artifact="$1" artifact_rel runtime_output
    [[ "${artifact}" == "${LOCAL_ROOT}/"* ]] || sync_fail "selected artifact is outside analysis root: ${artifact}"
    artifact_rel="${artifact#${LOCAL_ROOT}/}"
    case " ${selected_seen} " in
      *" ${artifact_rel} "*) return 0 ;;
    esac
    [[ -s "${artifact}" ]] || sync_fail "selected artifact is missing or empty: ${artifact}"
    if [[ "${artifact}" == *.runtime.json ]]; then
      runtime_output="${artifact%.runtime.json}"
      ecoda_validate_checksum "${runtime_output}" ||
        sync_fail "selected runtime output checksum failed: ${runtime_output}"
      benchmark_validate_runtime_metadata \
        "${artifact}" "${runtime_output}" "${ECODA_CHECKSUM_MD5}" ||
        sync_fail "selected runtime metadata validation failed: ${artifact}"
    fi
    ecoda_validate_checksum "${artifact}" || sync_fail "selected artifact checksum failed: ${artifact}"
    selected_seen="${selected_seen} ${artifact_rel}"
    SELECTED_RELS+=("${artifact_rel}")
    SELECTED_DIGESTS+=("${ECODA_CHECKSUM_MD5}")
    printf '%s\n' "${artifact_rel}" >> "${SYNC_FILES_TMP}"
    printf '%s\n' "${artifact_rel}.md5" >> "${SYNC_FILES_TMP}"
  }

  if [[ -r "${ECODA_SELECTION_MANIFEST:-}" ]]; then
    while IFS=$'\t' read -r ds view row_label; do
      [[ -n "${ds}" && -n "${view}" ]] || continue
      SYNC_LABELS=("${LABELS[@]}")
      if [[ "${ECODA_EXACT_SELECTION:-0}" == "1" && -z "${STAGE5_PASS}" ]]; then
        SYNC_LABELS=("${row_label}")
      fi
      for label in "${SYNC_LABELS[@]}"; do
        benchmark_sync_artifacts_for "${ds}" "${label}" || sync_fail "cannot resolve ${ds}/${view}/${label}"
        for path in "${SYNC_ARTIFACTS[@]}"; do add_sync_artifact "${path}"; done
      done
    done < "${ECODA_SELECTION_MANIFEST}"
  else
    for ds in "${DATASET_NAMES[@]}"; do
      for label in "${LABELS[@]}"; do
        benchmark_sync_artifacts_for "${ds}" "${label}" || sync_fail "cannot resolve ${ds}/${label}"
        for path in "${SYNC_ARTIFACTS[@]}"; do add_sync_artifact "${path}"; done
      done
    done
  fi
  if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
    metadata_manifest="${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
    expected_metadata_view="batch_effect_${STAGE5_PASS}"
    [[ -s "${metadata_manifest}" ]] ||
      sync_fail "${ANALYSIS_VARIANT} metadata export manifest is missing"
    while IFS=$'\t' read -r metadata_ds metadata_view metadata_input metadata_output metadata_extra; do
      [[ -n "${metadata_ds}" && "${metadata_view}" == "${expected_metadata_view}" &&
         -n "${metadata_input}" && -n "${metadata_output}" &&
         "${metadata_output}" == "${LOCAL_ROOT}/metadata/"* &&
         -z "${metadata_extra}" ]] ||
        sync_fail "malformed ${ANALYSIS_VARIANT} metadata export manifest row"
      add_sync_artifact "${metadata_output}"
    done < "${metadata_manifest}"
  fi
  [[ ${#SELECTED_RELS[@]} -gt 0 ]] || sync_fail "no selected benchmark artifacts"
  if [[ -s "${REMOTE_ROOT}/checksums.md5" ]]; then
    (cd "${REMOTE_ROOT}" && md5sum -c checksums.md5 >/dev/null 2>&1) || {
      sync_fail "existing remote checksum manifest failed validation"
    }
  fi

  normalize_external_checksum_path() {
    local path="$1" sidecar="${2:-${1}.md5}"
    local md5_line size_line path_line line_count tmp
    [[ -s "${path}" && -s "${sidecar}" ]] || return 1
    line_count="$(awk 'END { print NR }' "${sidecar}")" || return 1
    [[ "${line_count}" == "3" ]] || return 1
    md5_line="$(sed -n '1p' "${sidecar}")"
    size_line="$(sed -n '2p' "${sidecar}")"
    path_line="$(sed -n '3p' "${sidecar}")"
    [[ "${md5_line}" == MD5=* && "${size_line}" == SIZE=* &&
       "${path_line}" == PATH=* && -n "${path_line#PATH=}" ]] || return 1
    [[ "${path_line}" == "PATH=${path}" ]] && return 0
    tmp="${sidecar}.tmp.$$"
    if ! printf '%s\n%s\nPATH=%s\n' \
        "${md5_line}" "${size_line}" "${path}" > "${tmp}"; then
      rm -f "${tmp}"
      return 1
    fi
    if ! mv -f "${tmp}" "${sidecar}"; then
      rm -f "${tmp}"
      return 1
    fi
  }

  EXISTING_LOG="${REMOTE_ROOT}/embeddings/execution_times.feather"
  if [[ -f "${EXISTING_LOG}" && ! -f "${EXISTING_LOG}.md5" ]]; then
    "${PYTHON_BIN}" "${merge_script}" \
      --migrate-existing-log "${EXISTING_LOG}" || {
      sync_fail "existing remote execution log failed guarded sidecar migration"
    }
  fi
  if [[ -e "${EXISTING_LOG}" || -e "${EXISTING_LOG}.md5" ]]; then
    ecoda_validate_checksum_remote "${EXISTING_LOG}" "${EXISTING_LOG}.md5" || {
      sync_fail "existing remote execution log checksum failed"
    }
    case "${EXISTING_LOG}" in
      "${RUN_ROOT}"/*) ;;
      *)
        normalize_external_checksum_path "${EXISTING_LOG}" ||
          sync_fail "existing remote execution log sidecar schema failed"
        ;;
    esac
  fi
  MERGE_ARGS=(--output_dir "${LOCAL_ROOT}/embeddings"
              --log-dir "${RUN_LOG_DIR}"
              --no-cleanup
              --filename_prefix "${LOG_PREFIX}"
              --labels "${LABELS[@]}"
              --datasets "${DATASET_NAMES[@]}"
              --cleanup-manifest "${CLEANUP_MANIFEST}")
  [[ -s "${EXISTING_LOG}" ]] && MERGE_ARGS+=(--existing-log "${EXISTING_LOG}")
  "${PYTHON_BIN}" "${merge_script}" "${MERGE_ARGS[@]}" || {
    sync_fail "execution-log merge failed"
  }
  add_sync_artifact "${LOCAL_ROOT}/embeddings/execution_times.feather"

  CHECKSUM_TMP="${LOCAL_ROOT}/checksums.md5.build.$$"
  if ! : > "${CHECKSUM_TMP}"; then
    sync_fail "failed to create benchmark checksum manifest"
  fi
  if [[ -s "${REMOTE_ROOT}/checksums.md5" ]]; then
    while IFS= read -r line || [[ -n "${line}" ]]; do
      rel="$(printf '%s\n' "${line}" | awk '{p=$2; sub(/^\*/, "", p); print p}')"
      case " ${selected_seen} " in
        *" ${rel} "*) ;;
        *) printf '%s\n' "${line}" >> "${CHECKSUM_TMP}" ||
             sync_fail "failed to preserve remote checksum entry" ;;
      esac
    done < "${REMOTE_ROOT}/checksums.md5"
  fi
  [[ ${#SELECTED_RELS[@]} -eq ${#SELECTED_DIGESTS[@]} ]] ||
    sync_fail "selected checksum records are incomplete"
  selected_index=0
  for rel in "${SELECTED_RELS[@]}"; do
    printf '%s  %s\n' "${SELECTED_DIGESTS[${selected_index}]}" "${rel}" >> "${CHECKSUM_TMP}" || {
      rm -f "${CHECKSUM_TMP}"
      sync_fail "cannot write selected checksum record: ${rel}"
    }
    selected_index=$((selected_index + 1))
  done
  [[ -s "${CHECKSUM_TMP}" ]] || {
    rm -f "${CHECKSUM_TMP}"
    sync_fail "final benchmark checksum manifest is empty"
  }
  if ! mv -f "${CHECKSUM_TMP}" "${LOCAL_ROOT}/checksums.md5"; then
    rm -f "${CHECKSUM_TMP}"
    sync_fail "cannot install local benchmark checksum manifest"
  fi
  printf 'checksums.md5\n' >> "${SYNC_FILES_TMP}" ||
    sync_fail "cannot add checksum manifest to selected sync"
  [[ -s "${SYNC_FILES_TMP}" ]] || sync_fail "selected sync file manifest is empty"
  if ! mv -f "${SYNC_FILES_TMP}" "${SYNC_FILES}"; then
    rm -f "${SYNC_FILES_TMP}"
    sync_fail "cannot install selected sync manifest"
  fi

  if ! : > "${NO_CHECKSUM_FILES_TMP}"; then
    sync_fail "failed to create no-checksum sync manifest"
  fi
  while IFS= read -r rel || [[ -n "${rel}" ]]; do
    if [[ "${rel}" != "checksums.md5" ]]; then
      printf '%s\n' "${rel}" >> "${NO_CHECKSUM_FILES_TMP}" ||
        sync_fail "failed to build no-checksum sync manifest"
    fi
  done < "${SYNC_FILES}"
  rsync -rlptDv --files-from="${NO_CHECKSUM_FILES_TMP}" \
    "${LOCAL_ROOT}/" "${REMOTE_ROOT}/" || sync_fail "selected benchmark rsync failed"
  REMOTE_CHECKSUM_TMP="${REMOTE_ROOT}/.checksums.md5.tmp.${RUN_ID}"
  cp "${LOCAL_ROOT}/checksums.md5" "${REMOTE_CHECKSUM_TMP}" || {
    sync_fail "cannot stage remote checksum manifest"
  }
  mv -f "${REMOTE_CHECKSUM_TMP}" "${REMOTE_ROOT}/checksums.md5" || {
    sync_fail "cannot atomically install remote checksum manifest"
  }
  (cd "${REMOTE_ROOT}" && md5sum -c checksums.md5) || sync_fail "remote checksum verification failed"
  if [[ -s "${CLEANUP_MANIFEST}" ]]; then
    while IFS= read -r path || [[ -n "${path}" ]]; do
      [[ "${path}" == "${RUN_LOG_DIR}/"* ]] ||
        sync_fail "cleanup manifest escapes run log directory: ${path}"
      rm -f "${path}" "${path}.md5" ||
        sync_fail "failed to clean up run log artifact: ${path}"
    done < "${CLEANUP_MANIFEST}"
  fi
  notify_sync_status \
    "ECODA: ${KIND} synced to NAS" \
    "${KIND_CAP} selected results synced to ${REMOTE_ROOT}/ (datasets: ${DATASET_NAMES[*]}, labels: ${LABELS[*]}).
$(benchmark_job_durations_block)" || true
  SYNC_FINAL_STATE="OK"
)

# Existing submitters retain their public helper name; batch-effect submitters
# call the generic implementation directly.
benchmark_merge_sync_cleanup() {
  analysis_merge_sync_cleanup "$@"
}
