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
# and `matrix_gate.sh` owns the aggregate gate. Invoke the root launcher
# directly; family-specific compatibility entrypoints were removed.
# ============================================================================

# Path to the shared exec-log merge script, resolved from THIS file's location
# (BASH_SOURCE[0] inside a sourced file is the sourced file's path).
ANALYSIS_MERGE_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py"

if ! command -v ecoda_stage5_validate_identity >/dev/null 2>&1; then
  source "$(dirname "${BASH_SOURCE[0]}")/../utils/bash/ecoda_stage5_policy.sh"
fi

# Bind local evidence to the immutable composite, not its replaceable sidecar.
stage5_inventory_binds_file() {
  local receipt="$1" path="$2" root inventory relative digest size
  root="$(jq -er '.source_root' "${receipt}")" || return 1
  inventory="$(jq -er '.inventory_path' "${receipt}")" || return 1
  case "${path}" in "${root}/"*) relative="${path#"${root}/"}" ;; *) return 1 ;; esac
  [[ -f "${path}" && ! -L "${path}" ]] || return 1
  if [[ -s "${path}" ]]; then
    digest="$(ecoda_md5_file "${path}")" || return 1
  else
    digest=d41d8cd98f00b204e9800998ecf8427e
  fi
  size="$(wc -c < "${path}" | tr -d '[:space:]')" || return 1
  awk -F '\t' -v p="${relative}" -v d="${digest}" -v s="${size}" \
    '$1==p {n++; if($2!=d || $3!=s) bad=1} END {exit(n==1 && !bad ? 0 : 1)}' "${inventory}"
}

stage5_bind_precompute_evidence() {
  local receipt="$1" expected="$2" promotion child inputs nas
  promotion="$(jq -er '.source_root' "${receipt}")" || return 1
  for child in inputs canonical archive; do
    stage5_inventory_binds_file "${receipt}" "${promotion}/manifests/${child}.receipt.json" || return 1
  done
  inputs="${promotion}/manifests/inputs.receipt.json"
  nas="$(jq -er '.nas_root' "${receipt}")" || return 1
  # New releases require the normal dataset/output paths. Utility verification
  # enforces the exact configured H5AD/sidecar set for this canonical root;
  # historical completed-promotion readers remain read-only and separate.
  [[ "$(jq -er '.destination_root' "${inputs}")" == "${nas}" &&
     "$(jq -er '.nas_root' "${inputs}")" == "${nas}" &&
     "$(jq -er '.selection_sha256' "${inputs}")" == "$(jq -er '.selection_sha256' "${receipt}")" ]] || return 1
  stage5_inventory_binds_file "${receipt}" "${promotion}/manifests/expected_pending.tsv" || return 1
  [[ -f "${expected}" && ! -L "${expected}" ]] || return 1
  cmp -s "${expected}" "${promotion}/manifests/expected_pending.tsv"
}

benchmark_require_nas_mode() {
  "${PYTHON_BIN}" - "$1" "$2" <<'PY'
import os, stat, sys
p, mode = sys.argv[1:]
sys.exit(0 if not os.path.islink(p) and stat.S_IMODE(os.stat(p).st_mode) == int(mode, 8) else 1)
PY
}
benchmark_stage5_identity_guard() {
  local requested_pass="${1:-${PASS_ARG:-${ANALYSIS_PASS:-}}}"
  local scratch_root="${HPC_SCRATCH_DIR:-}" nas_root="${NAS_TARGET_DIR:-}"
  local expected_root expected_nas expected_pass expected_suffix
  local variant="${ANALYSIS_VARIANT:-}"
  if [[ -n "${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION:-}" ||
        -n "${ANALYSIS_ROOT_VERSION:-}" ||
        -n "${ECODA_STAGE5_CORRECTED_FINAL_ROOT_IDENTITY:-}" ||
        -n "${ANALYSIS_ROOT_IDENTITY:-}" ]]; then
    echo "ERROR: historical Stage 5 root identity is not permitted for a new run." >&2
    return 1
  fi
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
  local selection="${ECODA_SELECTION_MANIFEST:-}"
  local matrix="${ECODA_STAGE5_METHOD_MATRIX:-${METHOD_MATRIX:-}}"
  local nas_utility receipt="${ECODA_PRECOMPUTE_BACKUP_RECEIPT:-}" promotion_root=""
  local archive_receipt="" publication_identical=1
  local receipt_hash expected_hash child remote_dir next_dir
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
  local GLOBAL_SYNC_LOCK=""
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
    if [[ -n "${GLOBAL_SYNC_LOCK}" ]]; then
      rmdir "${GLOBAL_SYNC_LOCK}" 2>/dev/null || rc=1
    fi
    exit "${rc}"
  }
  trap sync_cleanup EXIT

  [[ -n "${RUN_ROOT}" && -n "${RUN_ID}" ]] || sync_fail "Stage 5 run root/id is not set"
  ecoda_validate_run_owned_path "${selection}" "${RUN_ROOT}" &&
    ecoda_validate_manifest "${selection}" 3 &&
    ecoda_validate_checksum "${selection}" ||
    sync_fail "explicit run-owned checksum-valid selection is required"
  if [[ "${ECODA_STAGE5_METHOD_MATRIX_MODE:-0}" == "1" || -n "${matrix}" ]]; then
    ecoda_validate_run_owned_path "${matrix}" "${RUN_ROOT}" &&
      ecoda_validate_manifest "${matrix}" 3 &&
      ecoda_validate_checksum "${matrix}" ||
      sync_fail "explicit run-owned checksum-valid method matrix is required"
  fi
  if [[ "${ANALYSIS_VARIANT:-}" == corrected_final ]]; then
    [[ -n "${receipt}" && -n "${ECODA_EXPECTED_PENDING:-}" &&
       -n "${ECODA_PRECOMPUTE_BACKUP_RECEIPT_SHA256:-}" &&
       -n "${ECODA_EXPECTED_PENDING_SHA256:-}" ]] ||
      sync_fail "corrected publication requires bound precompute receipt and pending scope"
    receipt_hash="$(shasum -a 256 "${receipt}")" || sync_fail "cannot hash precompute receipt"
    expected_hash="$(shasum -a 256 "${ECODA_EXPECTED_PENDING}")" || sync_fail "cannot hash expected pending"
    [[ "${receipt_hash%% *}" == "${ECODA_PRECOMPUTE_BACKUP_RECEIPT_SHA256}" &&
       "${expected_hash%% *}" == "${ECODA_EXPECTED_PENDING_SHA256}" ]] ||
      sync_fail "precompute receipt or expected pending binding changed"
  fi
  if [[ "${LOCAL_ROOT}" != "${HPC_SCRATCH_DIR}/benchmark" ]]; then
    KIND="analysis"
    KIND_CAP="Analysis"
  fi
  nas_utility="$(benchmark_source_script_path src/utils/bash/ecoda_nas_backup.sh)"
  benchmark_require_source_script_path "${nas_utility}" ||
    sync_fail "verified NAS utility is unavailable"
  if [[ -n "${receipt}" ]]; then
    bash "${nas_utility}" verify --receipt "${receipt}" ||
      sync_fail "precompute backup receipt verification failed"
    if [[ "${ANALYSIS_VARIANT:-}" == corrected_final ]]; then
      stage5_bind_precompute_evidence "${receipt}" "${ECODA_EXPECTED_PENDING}" ||
        sync_fail "precompute composite evidence binding changed"
    fi
    [[ "$(jq -r '.nas_root' "${receipt}")" == "${NAS_TARGET_DIR}" ]] ||
      sync_fail "precompute backup NAS root differs"
    promotion_root="$(jq -er '.source_root' "${receipt}")" ||
      sync_fail "precompute promotion source is missing"
    for child in inputs canonical archive; do
      bash "${nas_utility}" verify --receipt "${promotion_root}/manifests/${child}.receipt.json" ||
        sync_fail "precompute ${child} receipt verification failed"
    done
    archive_receipt="${promotion_root}/manifests/archive.receipt.json"
  fi
  bash "${nas_utility}" preflight --nas-root "${NAS_TARGET_DIR}" \
    --report "${RUN_ROOT}/manifests/sync_nas_preflight.json" ||
    sync_fail "NAS mount/write/rename/effective-mode preflight failed"
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
  # Inspect the one global key before any merge. Reclaim is deferred until the
  # complete intended publication exists and its archived predecessor is proved.
  recovery_sync_owner_dir="$(ecoda_owner_dir stage5 "sync/${LOCAL_ROOT}")" ||
    sync_fail "cannot resolve shared Stage 5 sync owner"
  mkdir -p "$(dirname "${recovery_sync_owner_dir}")" ||
    sync_fail "cannot create shared sync lock parent"
  if ! mkdir "${recovery_sync_owner_dir}.publication.lock" 2>/dev/null; then
    sync_fail "shared Stage 5 publication lock is already held"
  fi
  GLOBAL_SYNC_LOCK="${recovery_sync_owner_dir}.publication.lock"
  if [[ -e "${recovery_sync_owner_dir}" || -L "${recovery_sync_owner_dir}" ]]; then
    [[ -d "${recovery_sync_owner_dir}" && ! -L "${recovery_sync_owner_dir}" &&
       -f "${recovery_sync_owner_dir}/owner" && ! -L "${recovery_sync_owner_dir}/owner" ]] ||
      sync_fail "shared Stage 5 sync owner is malformed"
    "${PYTHON_BIN}" - "${recovery_sync_owner_dir}/owner" <<'PY' ||
import pathlib, re, sys
try:
    fields = {}
    for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
        key, value = line.split("=", 1)
        if key in fields:
            raise ValueError("duplicate owner field")
        fields[key] = value
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", fields["RUN_ID"]):
        raise ValueError("invalid owner run")
    if fields["STATE"] not in ("OK", "FAIL"):
        raise ValueError("nonterminal owner")
except (OSError, ValueError, KeyError):
    sys.exit(1)
PY
      sync_fail "shared Stage 5 sync owner is active or malformed"
    recovery_sync_owner_state="$(ecoda_owner_state "${recovery_sync_owner_dir}" 2>/dev/null)" ||
      sync_fail "shared Stage 5 sync owner is malformed"
    [[ "$(ecoda_owner_field "${recovery_sync_owner_dir}" STAGE)" == stage5 &&
       "$(ecoda_owner_field "${recovery_sync_owner_dir}" KEY)" == "sync/${LOCAL_ROOT}" &&
       -n "$(ecoda_owner_run "${recovery_sync_owner_dir}")" ]] ||
      sync_fail "shared Stage 5 sync owner identity is malformed"
    case "${recovery_sync_owner_state}" in
      OK|FAIL) ;;
      *) sync_fail "shared Stage 5 sync owner is active or malformed" ;;
    esac
  fi
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

  while IFS=$'\t' read -r ds view row_label; do
    if [[ -n "${STAGE5_PASS}" ]]; then
      [[ "${view}" == "batch_effect_${STAGE5_PASS}" ]] ||
        sync_fail "selection view disagrees with synchronization pass"
    else
      [[ "${view}" == benchmark_analysis ]] ||
        sync_fail "selection view disagrees with benchmark synchronization"
    fi
    SYNC_LABELS=("${LABELS[@]}")
    if [[ "${ECODA_EXACT_SELECTION:-0}" == "1" && -z "${STAGE5_PASS}" ]]; then
      SYNC_LABELS=("${row_label}")
    fi
    for label in "${SYNC_LABELS[@]}"; do
      if [[ "${ECODA_STAGE5_METHOD_MATRIX_MODE:-0}" == "1" || -n "${matrix}" ]]; then
        ecoda_stage5_method_matrix_allows "${ds}" "${view}" "${label}" || continue
      fi
      benchmark_sync_artifacts_for "${ds}" "${label}" || sync_fail "cannot resolve ${ds}/${view}/${label}"
      for path in "${SYNC_ARTIFACTS[@]}"; do add_sync_artifact "${path}"; done
    done
  done < "${selection}"
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
      awk -F '\t' -v ds="${metadata_ds}" -v view="${metadata_view}" \
        '$1 == ds && $2 == view { found=1 } END { exit !found }' "${selection}" ||
        sync_fail "metadata artifact is outside the explicit selection"
      add_sync_artifact "${metadata_output}"
    done < "${metadata_manifest}"
  fi
  [[ ${#SELECTED_RELS[@]} -gt 0 ]] || sync_fail "no selected benchmark artifacts"
  local retry_publication=0 retry_continuity=""
  local transition="${RUN_ROOT}/manifests/sync_transition.json"
  if [[ "${ANALYSIS_VARIANT:-}" == corrected_final && "${recovery_sync_owner_state}" == FAIL &&
        "$(ecoda_owner_run "${recovery_sync_owner_dir}")" == "${RUN_ID}" &&
        ( -e "${transition}" || -L "${transition}" ) ]]; then
    ecoda_validate_run_owned_path "${transition}" "${RUN_ROOT}" &&
      ecoda_validate_checksum "${transition}" &&
      ecoda_validate_artifact_record "${transition}" stage5_sync_transition "${RUN_ID}" ||
      sync_fail "publication transition evidence is invalid"
    [[ "${recovery_sync_owner_state}" == FAIL &&
       "$(ecoda_owner_run "${recovery_sync_owner_dir}")" == "${RUN_ID}" ]] ||
      sync_fail "publication retry requires this run's terminal FAIL owner"
    retry_continuity="$("${PYTHON_BIN}" - "${transition}" "${receipt}" "${RUN_ID}" \
      "${LOCAL_ROOT}" "${REMOTE_ROOT}" "${RUN_ROOT}" "${recovery_sync_owner_dir}/owner" <<'PY'
import hashlib, json, os, pathlib, sys
transition, receipt, run_id, local, remote, run_root, owner = sys.argv[1:]
def identity(path):
    p = pathlib.Path(path)
    if p.is_symlink():
        raise ValueError("symlink publication evidence")
    if not p.exists():
        return None
    if not p.is_file():
        raise ValueError("nonregular publication evidence")
    h = hashlib.md5()
    with p.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return [h.hexdigest(), str(p.stat().st_size)]
def selected_path(rel):
    p = pathlib.PurePosixPath(rel)
    if p.is_absolute() or ".." in p.parts or str(p) != rel:
        raise ValueError("unsafe publication path")
    target = pathlib.Path(remote) / p
    if target.resolve() != target:
        raise ValueError("symlink publication ancestor")
    return target
try:
    t = json.loads(pathlib.Path(transition).read_text())
    if (t["run_id"], t["local_root"], t["remote_root"], t["owner"]) != (run_id, local, remote, owner):
        raise ValueError("transition identity differs")
    if t["composite_sha256"] != hashlib.sha256(pathlib.Path(receipt).read_bytes()).hexdigest():
        raise ValueError("transition composite differs")
    if set(t["original"]) != set(t["publication"]):
        raise ValueError("transition selection differs")
    for rel, intended in t["publication"].items():
        if identity(selected_path(rel)) not in (t["original"][rel], intended):
            raise ValueError("selected destination is neither original nor intended")
    for rel, original in t["unselected"].items():
        if rel in t["publication"] or identity(selected_path(rel)) != original:
            raise ValueError("unselected canonical destination changed")
    continuity = t["continuity"]
    if continuity:
        p = pathlib.Path(continuity["path"])
        if p.parent != pathlib.Path(run_root) / "manifests" or p.resolve() != p:
            raise ValueError("continuity outside run")
        if identity(p) != continuity["identity"]:
            raise ValueError("continuity changed")
        print(p)
except (OSError, ValueError, KeyError, TypeError):
    sys.exit(1)
PY
    )" || sync_fail "partial publication differs from sealed transition"
    retry_publication=1
  fi
  if [[ ${retry_publication} -eq 0 && -s "${REMOTE_ROOT}/checksums.md5" ]]; then
    (cd "${REMOTE_ROOT}" && md5sum -c checksums.md5 >/dev/null 2>&1) || {
      sync_fail "existing remote checksum manifest failed validation"
    }
  fi

  EXISTING_LOG="${REMOTE_ROOT}/embeddings/execution_times.feather"
  if [[ ${retry_publication} -eq 1 ]]; then
    EXISTING_LOG="${retry_continuity}"
  elif [[ -e "${EXISTING_LOG}" || -e "${EXISTING_LOG}.md5" ]]; then
    ecoda_validate_checksum_remote "${EXISTING_LOG}" "${EXISTING_LOG}.md5" ||
      sync_fail "existing remote execution log checksum failed"
    # NAS sidecars intentionally retain their source PATH. Give the merger a
    # run-owned continuity copy rather than rewriting historical NAS evidence.
    local continuity_md5 continuity_log
    continuity_md5="$(ecoda_md5_file "${EXISTING_LOG}")" ||
      sync_fail "cannot identify existing execution log"
    continuity_log="${RUN_ROOT}/manifests/execution_times_existing_${continuity_md5}.feather"
    if [[ -e "${continuity_log}" || -L "${continuity_log}" ]]; then
      ecoda_validate_run_owned_path "${continuity_log}" "${RUN_ROOT}" &&
        ecoda_validate_checksum "${continuity_log}" &&
        cmp -s "${EXISTING_LOG}" "${continuity_log}" ||
        sync_fail "existing run-owned continuity copy differs"
    else
      ecoda_atomic_copy "${EXISTING_LOG}" "${continuity_log}" &&
        cmp -s "${EXISTING_LOG}" "${continuity_log}" &&
        ecoda_write_checksum "${continuity_log}" ||
        sync_fail "cannot stage run-owned execution-log continuity copy"
    fi
    ecoda_write_artifact_record "${continuity_log}" stage5_execution_log "${RUN_ID}" >/dev/null ||
      sync_fail "cannot publish continuity log record"
    EXISTING_LOG="${continuity_log}"
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
  # Byte-identical publication is validation-only: do not reclaim its owner,
  # rewrite payloads/sidecars/checksums, or normalize already-correct modes.
  while IFS= read -r rel; do
    [[ ! -L "${LOCAL_ROOT}/${rel}" && ! -L "${REMOTE_ROOT}/${rel}" ]] ||
      sync_fail "symlink in selected publication: ${rel}"
    if ! cmp -s "${LOCAL_ROOT}/${rel}" "${REMOTE_ROOT}/${rel}"; then
      publication_identical=0
    fi
  done < "${SYNC_FILES}"
  if [[ "${publication_identical}" == 1 && ${retry_publication} -eq 0 ]]; then
    benchmark_require_nas_mode "${REMOTE_ROOT}" 777 || sync_fail "NAS root effective mode is not 777"
    while IFS= read -r rel; do
      benchmark_require_nas_mode "${REMOTE_ROOT}/${rel}" 666 ||
        sync_fail "NAS effective file mode is not 666: ${rel}"
      remote_dir="$(dirname "${REMOTE_ROOT}/${rel}")"
      while [[ "${remote_dir}" != "${REMOTE_ROOT}" ]]; do
        benchmark_require_nas_mode "${remote_dir}" 777 ||
          sync_fail "NAS effective directory mode is not 777: ${remote_dir}"
        remote_dir="$(dirname "${remote_dir}")"
      done
    done < "${SYNC_FILES}"
    echo "Validated identical intended NAS publication; no publication or owner rewrite."
    exit 0
  fi
  if [[ "${ANALYSIS_VARIANT:-}" == corrected_final ]]; then
    if [[ -e "${transition}" || -L "${transition}" ]]; then
      ecoda_validate_run_owned_path "${transition}" "${RUN_ROOT}" &&
        ecoda_validate_checksum "${transition}" &&
        ecoda_validate_artifact_record "${transition}" stage5_sync_transition "${RUN_ID}" ||
        sync_fail "publication transition evidence is invalid"
    fi
    [[ -n "${archive_receipt}" ]] || sync_fail "terminal shared sync transition requires archive receipt"
    "${PYTHON_BIN}" - "${receipt}" "${archive_receipt}" "${promotion_root}" \
      "${recovery_sync_owner_dir}/owner" "${SYNC_FILES}" "${LOCAL_ROOT}" "${REMOTE_ROOT}" \
      "${transition}" "${RUN_ID}" "${EXISTING_LOG}" <<'PY' ||
import csv, hashlib, json, os, pathlib, sys
receipt_path, archive_path, promotion, owner, files, local, remote, transition_path, run_id, continuity = sys.argv[1:]
def identity(path):
    p = pathlib.Path(path)
    if p.is_symlink() or not p.is_file():
        raise ValueError("nonregular evidence")
    digest = hashlib.md5()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest(), str(p.stat().st_size)
def inventory(path):
    rows = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) != 3 or row[0] in rows:
                raise ValueError("invalid inventory")
            rows[row[0]] = (row[1], row[2])
    return rows
try:
    main = json.load(open(receipt_path))
    composite = inventory(main["inventory_path"])
    required = ["manifests/global_sync_owner.tsv", "manifests/archive_mapping.tsv"]
    for kind in ("inputs", "canonical", "archive"):
        child_rel = "manifests/" + kind + ".receipt.json"
        required.append(child_rel)
        child = json.load(open(os.path.join(promotion, child_rel)))
        required.append(os.path.relpath(child["inventory_path"], promotion))
    for rel in required:
        if rel not in composite or identity(os.path.join(promotion, rel)) != composite[rel]:
            raise ValueError("prerequisite missing from composite")
    archive = json.load(open(archive_path))
    archived = inventory(archive["inventory_path"])
    mapping = {}
    with open(os.path.join(promotion, "manifests/archive_mapping.tsv")) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != ["original_path", "archive_path", "md5", "size"]:
            raise ValueError("invalid archive mapping header")
        for row in reader:
            if row["original_path"] in mapping:
                raise ValueError("duplicate archive mapping")
            mapping[row["original_path"]] = row
    def require_archived(path):
        row = mapping[path]
        target = row["archive_path"]
        rel = os.path.relpath(target, archive["destination_root"])
        if rel == ".." or rel.startswith("../") or rel not in archived:
            raise ValueError("publication evidence not archived")
        if (row["md5"], row["size"]) != archived[rel] or identity(path) != archived[rel]:
            raise ValueError("publication evidence changed since archive")
        if identity(target) != archived[rel]:
            raise ValueError("archive evidence differs")
        return target
    with open(os.path.join(promotion, "manifests/global_sync_owner.tsv")) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    matches = [r for r in rows if r["path"] == owner]
    if len(matches) > 1:
        raise ValueError("duplicate owner authorization")
    row = matches[0] if matches else None
    transition = pathlib.Path(transition_path)
    intent = {
        rel: identity(os.path.join(local, rel))
        for rel in pathlib.Path(files).read_text().splitlines()
    }
    authorization = dict(
        run_id=run_id, owner=owner, local_root=local, remote_root=remote,
        composite_sha256=hashlib.sha256(pathlib.Path(receipt_path).read_bytes()).hexdigest(),
        predecessor=row, publication={k: list(v) for k, v in intent.items()},
    )
    retry = transition.exists()
    if retry:
        stored = json.loads(transition.read_text())
        if transition.is_symlink() or any(stored.get(k) != v for k, v in authorization.items()):
            raise ValueError("publication transition intent changed")
    elif row is None:
        if os.path.lexists(owner):
            raise ValueError("unrecorded owner appeared")
    elif identity(owner) != (row["md5"], row["size"]):
        raise ValueError("owner changed")
    fields = {}
    if os.path.lexists(owner):
        for line in pathlib.Path(owner).read_text().splitlines():
            key, value = line.split("=", 1)
            if key in fields:
                raise ValueError("duplicate owner field")
            fields[key] = value
    if retry:
        if fields.get("STATE") != "FAIL" or fields.get("RUN_ID") != run_id:
            raise ValueError("retry requires this run's terminal FAIL owner")
        if row is not None:
            predecessor = row["archive_path"]
            if identity(predecessor) != (row["md5"], row["size"]):
                raise ValueError("archived predecessor changed")
            archived_fields = dict(line.split("=", 1) for line in pathlib.Path(predecessor).read_text().splitlines())
            if archived_fields["STATE"] != row["state"] or archived_fields["RUN_ID"] != row["run_id"]:
                raise ValueError("archived predecessor identity differs")
            mapped = mapping[owner]
            rel = os.path.relpath(predecessor, archive["destination_root"])
            if (mapped["archive_path"] != predecessor or rel not in archived
                    or archived[rel] != identity(predecessor)
                    or (mapped["md5"], mapped["size"]) != (row["md5"], row["size"])):
                raise ValueError("predecessor archive mapping differs")
    elif row is not None:
        if fields["STATE"] not in ("OK", "FAIL") or fields["STATE"] != row["state"] or fields["RUN_ID"] != row["run_id"]:
            raise ValueError("owner is not the archived terminal owner")
        if require_archived(owner) != row["archive_path"]:
            raise ValueError("owner archive mapping differs")
    for rel in pathlib.Path(files).read_text().splitlines():
        target = os.path.join(remote, rel)
        if os.path.lexists(target) and identity(target) != identity(os.path.join(local, rel)):
            require_archived(target)
    if not retry:
        authorization["original"] = {
            rel: list(identity(os.path.join(remote, rel))) if os.path.lexists(os.path.join(remote, rel)) else None
            for rel in intent
        }
        authorization["unselected"] = {}
        checksums = pathlib.Path(remote) / "checksums.md5"
        if checksums.exists():
            for line in checksums.read_text().splitlines():
                digest, rel = line.split(None, 1)
                rel = rel.lstrip("*")
                p = pathlib.PurePosixPath(rel)
                if p.is_absolute() or ".." in p.parts or str(p) != rel:
                    raise ValueError("unsafe canonical checksum path")
                if rel not in intent:
                    actual = identity(os.path.join(remote, rel))
                    if actual[0] != digest or rel in authorization["unselected"]:
                        raise ValueError("unselected canonical checksum differs")
                    authorization["unselected"][rel] = list(actual)
        authorization["continuity"] = (
            dict(path=continuity, identity=list(identity(continuity)))
            if continuity and os.path.isfile(continuity) else None
        )
        with transition.open("x") as stream:
            json.dump(authorization, stream, sort_keys=True)
            stream.write("\n")
except (OSError, ValueError, KeyError, TypeError):
    sys.exit(1)
PY
      sync_fail "terminal shared sync owner/publication archive authorization failed"
    if [[ ! -e "${transition}.md5" ]]; then
      ecoda_write_checksum "${transition}" &&
        ecoda_write_artifact_record "${transition}" stage5_sync_transition "${RUN_ID}" >/dev/null ||
        sync_fail "cannot seal publication transition"
    fi
    sync_owner_artifact_valid=0
  elif [[ -n "${recovery_sync_owner_state}" ]]; then
    [[ "${recovery_sync_owner_state}" == FAIL ]] ||
      sync_fail "shared Stage 5 sync owner is already finalized successfully"
    sync_owner_artifact_valid=0
  fi
  SYNC_OWNER_DIR="$(ecoda_owner_acquire stage5 "sync/${LOCAL_ROOT}" "${RUN_ID}" 0 \
    "${sync_owner_artifact_valid}")" || sync_fail "shared Stage 5 sync owner is unavailable"

  if ! : > "${NO_CHECKSUM_FILES_TMP}"; then
    sync_fail "failed to create no-checksum sync manifest"
  fi
  while IFS= read -r rel || [[ -n "${rel}" ]]; do
    if [[ "${rel}" != "checksums.md5" ]]; then
      printf '%s\n' "${rel}" >> "${NO_CHECKSUM_FILES_TMP}" ||
        sync_fail "failed to build no-checksum sync manifest"
    fi
  done < "${SYNC_FILES}"
  rsync -rlptDv --checksum --chmod=ugo+rwX --files-from="${NO_CHECKSUM_FILES_TMP}" \
    "${LOCAL_ROOT}/" "${REMOTE_ROOT}/" || sync_fail "selected benchmark rsync failed"
  REMOTE_CHECKSUM_TMP="${REMOTE_ROOT}/.checksums.md5.tmp.${RUN_ID}"
  cp "${LOCAL_ROOT}/checksums.md5" "${REMOTE_CHECKSUM_TMP}" || {
    sync_fail "cannot stage remote checksum manifest"
  }
  chmod 666 "${REMOTE_CHECKSUM_TMP}" || {
    sync_fail "cannot set NAS checksum manifest permissions"
  }
  benchmark_require_nas_mode "${REMOTE_CHECKSUM_TMP}" 666 ||
    sync_fail "NAS checksum manifest effective mode is not 666"
  mv -f "${REMOTE_CHECKSUM_TMP}" "${REMOTE_ROOT}/checksums.md5" || {
    sync_fail "cannot atomically install remote checksum manifest"
  }
  (cd "${REMOTE_ROOT}" && md5sum -c checksums.md5) || sync_fail "remote checksum verification failed"
  chmod 777 "${REMOTE_ROOT}" ||
    sync_fail "cannot set NAS benchmark root permissions"
  benchmark_require_nas_mode "${REMOTE_ROOT}" 777 ||
    sync_fail "NAS root effective mode is not 777"
  while IFS= read -r rel || [[ -n "${rel}" ]]; do
    [[ -n "${rel}" ]] || continue
    chmod 666 "${REMOTE_ROOT}/${rel}" ||
      sync_fail "cannot set NAS file permissions: ${rel}"
    benchmark_require_nas_mode "${REMOTE_ROOT}/${rel}" 666 ||
      sync_fail "NAS effective file mode is not 666: ${rel}"
    remote_dir="${REMOTE_ROOT}/$(dirname "${rel}")"
    while [[ "${remote_dir}" != "${REMOTE_ROOT}" ]]; do
      chmod 777 "${remote_dir}" ||
        sync_fail "cannot set NAS directory permissions: ${remote_dir}"
      benchmark_require_nas_mode "${remote_dir}" 777 ||
        sync_fail "NAS effective directory mode is not 777: ${remote_dir}"
      next_dir="$(dirname "${remote_dir}")"
      [[ "${next_dir}" != "${remote_dir}" ]] ||
        sync_fail "NAS permission path did not converge: ${remote_dir}"
      remote_dir="${next_dir}"
    done
  done < "${SYNC_FILES}"
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
