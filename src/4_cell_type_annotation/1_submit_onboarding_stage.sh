#!/bin/bash
# Canonical Stage 4 gate: preparation and merges are dataset-parallel; view
# updates remain serial within one dataset because they share one annotation set.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
export ECODA_GATE_STAGE=stage4
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"

DATASETS_ARG=""
DATASETS_SET=0
VIEWS_ARG=""
VIEWS_SET=0
SELECTION_FILE_ARG=""
SELECTION_FILE_SET=0
EXACT_BATCH_SELECTION=0
FORCE_ARG=0
SKIP_PREPARE=0
REUSE_RUN_ARG=""
SYNC_ONLY_RUN=""
SYNC_ONLY_SET=0
MEMORY="32G"
MAX_MEMORY="128G"
PARTITION="${SLURM_PARTITION_BENCHMARK_CPU:-shared-cpu}"
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"

usage() {
  cat <<'EOF'
Usage: 1_submit_onboarding_stage.sh [--datasets LIST] [--views LIST]
       [--selection-file TSV] [--exact-batch-selection] [--skip-prepare
       --reuse-run RUN_ID] [--force] [--sync-only RUN_ID]
       [--mem VALUE] [--max-mem VALUE] [--partition NAME] [--throttle N]

Selection rows are DATASET<TAB>VIEW. --view is accepted as a compatibility
alias for --views. Exact batch mode requires the immutable twelve-row
uncorrected selection file.
EOF
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets) DATASETS_ARG="${2:-}"; DATASETS_SET=1; shift 2 ;;
    --datasets=*) DATASETS_ARG="${1#*=}"; DATASETS_SET=1; shift ;;
    --view|--views) VIEWS_ARG="${2:-}"; VIEWS_SET=1; shift 2 ;;
    --view=*|--views=*) VIEWS_ARG="${1#*=}"; VIEWS_SET=1; shift ;;
    --selection-file) SELECTION_FILE_ARG="${2:-}"; SELECTION_FILE_SET=1; shift 2 ;;
    --selection-file=*) SELECTION_FILE_ARG="${1#*=}"; SELECTION_FILE_SET=1; shift ;;
    --exact-batch-selection) EXACT_BATCH_SELECTION=1; shift ;;
    --reuse-run) REUSE_RUN_ARG="${2:-}"; shift 2 ;;
    --reuse-run=*) REUSE_RUN_ARG="${1#*=}"; shift ;;
    --skip-prepare) SKIP_PREPARE=1; shift ;;
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
if [[ ${SKIP_PREPARE} -eq 1 && -z "${REUSE_RUN_ARG}" ]]; then
  echo "ERROR: --skip-prepare requires --reuse-run RUN_ID; no scheduler submission made." >&2
  exit 1
fi
if [[ -n "${REUSE_RUN_ARG}" && ${SKIP_PREPARE} -eq 0 ]]; then
  echo "ERROR: --reuse-run is valid only with --skip-prepare." >&2
  exit 1
fi
if [[ -n "${REUSE_RUN_ARG}" && -n "${SYNC_ONLY_RUN}" ]]; then
  echo "ERROR: --reuse-run cannot be combined with --sync-only." >&2
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
  echo "ERROR: --sync-only requires a run ID." >&2
  exit 1
fi
command -v jq >/dev/null 2>&1 || { echo "ERROR: jq is required for Stage 4 selection." >&2; exit 1; }

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
  local ds="$1" view="$2" name
  name="$(ecoda_view_output_name "${ds}" "${view}")"
  [[ -n "${name}" ]] || return 1
  printf '%s/%s/output/%s' "${HPC_SCRATCH_DIR}" "${ds}" "${name}"
}

stage4_finalize_owner_manifest() {
  local state="$1" reason="$2" owners_file="${ECODA_RUN_ROOT}/manifests/owners.tsv"
  local ds owner rc=0
  [[ -r "${owners_file}" ]] || return 1
  [[ -s "${owners_file}" ]] || return 0
  while IFS=$'\t' read -r ds owner; do
    [[ -n "${ds}" && -n "${owner}" ]] || { rc=1; continue; }
    if ! ecoda_owner_set_state "${owner}" "${state}" "${reason}"; then
      rc=1
    fi
  done < "${owners_file}"
  return "${rc}"
}

stage4_abort() {
  local reason="$1"
  local rc=0
  ecoda_owner_finalize_tracked FAIL "${reason}" || rc=1
  if [[ -n "${ECODA_RUN_ROOT:-}" && -r "${ECODA_RUN_ROOT}/manifests/owners.tsv" ]]; then
    stage4_finalize_owner_manifest FAIL "${reason}" || rc=1
  fi
  if [[ -n "${ECODA_RUN_ROOT:-}" ]]; then
    ecoda_set_run_state FAIL "${reason}" || rc=1
  fi
  echo "ERROR: ${reason}" >&2
  exit 1
}
stage4_record_scheduler() {
  local kind="$1" scheduler_id="$2"
  local manifest="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
  local tmp="${manifest}.record.$$" existing_kind existing_id
  [[ "${kind}" == "SCGATE" || "${kind}" == "ARRAY" ||
     "${kind}" == "WATCHDOG" || "${kind}" == "STATUS" ]] || return 1
  [[ "${scheduler_id}" =~ ^[0-9]+$ ]] || return 1
  if [[ -s "${manifest}" ]]; then
    while IFS=$'\t' read -r existing_kind existing_id; do
      [[ -n "${existing_kind}" && "${existing_id}" =~ ^[0-9]+$ ]] || return 1
      [[ "${existing_id}" == "${scheduler_id}" ]] && return 0
    done < "${manifest}"
    cp "${manifest}" "${tmp}" || return 1
  else
    : > "${tmp}" || return 1
  fi
  printf '%s\t%s\n' "${kind}" "${scheduler_id}" >> "${tmp}" || {
    rm -f "${tmp}"
    return 1
  }
  mv -f "${tmp}" "${manifest}" || {
    rm -f "${tmp}"
    return 1
  }
}


validate_reuse_preparation() {
  local prep_manifest="${ECODA_RUN_ROOT}/manifests/preparation.tsv"
  local chunk_manifest="${ECODA_RUN_ROOT}/manifests/chunks.tsv"
  local ds views run_root chunk feather_dir chunk_num feather union_path sample expected_args
  [[ -s "${prep_manifest}" && -s "${chunk_manifest}" ]] || {
    echo "ERROR: --skip-prepare reuse run lacks immutable preparation/chunk manifests." >&2
    return 1
  }
  ecoda_validate_run_owned_path "${prep_manifest}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_run_owned_path "${chunk_manifest}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_manifest "${prep_manifest}" 3 || return 1
  ecoda_validate_manifest "${chunk_manifest}" 3 || return 1
  ecoda_validate_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv" || return 1
  while IFS=$'\t' read -r ds views run_root; do
    selected_dataset="$(awk -F '\t' -v target="${ds}" '$1 == target {found=1} END {print found+0}' "${MANIFEST}")"
    [[ "${selected_dataset}" == "1" ]] || return 1
    unsuitable="$(jq -r --arg ds "${ds}" '.[$ds].not_suitable_for_auto_annotation // [] | (index("hitme") != null and index("scatomic") != null)' "${DATASETS_JSON_FILE}")"
    [[ "${unsuitable}" != "true" ]] || return 1
    [[ "${run_root}" == "${ECODA_RUN_ROOT}" ]] || {
      echo "ERROR: reuse preparation row is not owned by ${ECODA_RUN_ROOT}: ${ds}" >&2
      return 1
    }
    export DS_NAME="${ds}" ANNOTATION_VIEWS="${views}" ANNOTATION_RUN_ROOT="${run_root}" ANNOTATION_RUN_ID="${RUN_ID}"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/1.1_prepare_chunks.py" \
      --views "${views}" --run-root "${run_root}" --validate-only || return 1
  done < "${prep_manifest}"
  while IFS=$'\t' read -r ds chunk feather_dir; do
    [[ -n "${ds}" && -n "${chunk}" && -n "${feather_dir}" ]] || return 1
    expected_chunk_dir="${ECODA_RUN_ROOT}/datasets/${ds}/chunks"
    expected_feather_dir="${ECODA_RUN_ROOT}/datasets/${ds}/annotations"
    expected_union="${ECODA_RUN_ROOT}/datasets/${ds}/union/union.h5ad"
    [[ "${chunk}" == "${expected_chunk_dir}/chunk_"*.txt ]] || return 1
    [[ "${feather_dir}" == "${expected_feather_dir}" ]] || return 1
    ecoda_validate_run_owned_path "${chunk}" "${ECODA_RUN_ROOT}" || return 1
    chunk_num="${chunk##*/chunk_}"
    chunk_num="${chunk_num%.txt}"
    [[ "${chunk_num}" =~ ^[1-9][0-9]*$ ]] || return 1
    feather="${feather_dir}/annotations_chunk_${chunk_num}.feather"
    union_path="$(sed -n '1p' "${chunk}")"
    [[ "${union_path}" == "${expected_union}" ]] || return 1
    ecoda_validate_run_owned_path "${union_path}" "${ECODA_RUN_ROOT}" || return 1
    ecoda_validate_checksum "${union_path}" || return 1
    expected_args=(--expected-union "${union_path}")
    while IFS= read -r sample; do
      [[ -n "${sample}" ]] || return 1
      expected_args+=(--expected-sample "${sample}")
    done < <(sed -n '2,$p' "${chunk}")
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
      --path "${feather}" --require-sidecar "${expected_args[@]}" || {
      echo "ERROR: reuse annotation feather failed validation: ${feather}" >&2
      return 1
    }
  done < "${chunk_manifest}"
}

sync_one_dataset() {
  local ds="$1" manifest="$2" sync_dir="${ECODA_RUN_ROOT}/status/sync"
  local path view output remote_dir dataset_root lock
  local transfer_manifest verify_manifest rc=0 selected=0
  mkdir -p "${sync_dir}" || return 1
  lock="${sync_dir}/${ds}.lock"
  mkdir "${lock}" 2>/dev/null || return 1
  remote_dir="${NAS_TARGET_DIR}/${ds}/output"
  transfer_manifest="${lock}/files"
  verify_manifest="${lock}/verify"
  if ! mkdir -p "${remote_dir}" ||
     ! : > "${transfer_manifest}" ||
     ! : > "${verify_manifest}"; then
    rc=1
  fi
  if [[ ${rc} -eq 0 ]]; then
    while IFS=$'\t' read -r row_ds view; do
      [[ "${row_ds}" == "${ds}" ]] || continue
      selected=1
      path="$(output_path_for "${row_ds}" "${view}")" || { rc=1; continue; }
      output="$(basename "${path}")"
      if ! ecoda_validate_checksum "${path}" ||
         ! "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
           --h5ad "${path}" --require-sidecar >/dev/null 2>&1; then
        rc=1
        continue
      fi
      printf '%s\n' "${output}" >> "${transfer_manifest}" || rc=1
      printf '%s\t%s\t%s\n' "${view}" "${path}" "${output}" >> "${verify_manifest}" || rc=1
      printf '%s\n' "${output}.md5" >> "${transfer_manifest}" || rc=1
    done < "${manifest}"
  fi
  [[ ${selected} -eq 1 ]] || rc=1
  if [[ ${rc} -eq 0 && -s "${transfer_manifest}" ]]; then
    rsync -rlptDv --files-from="${transfer_manifest}" \
      "${HPC_SCRATCH_DIR}/${ds}/output/" "${remote_dir}/" || rc=1
  elif [[ ${rc} -eq 0 ]]; then
    rc=1
  fi
  if [[ ${rc} -eq 0 ]]; then
    while IFS=$'\t' read -r view path output; do
      ecoda_compare_checksum_remote "${path}" \
        "${remote_dir}/${output}" "${remote_dir}/${output}.md5" || rc=1
    done < "${verify_manifest}"
  fi
  if [[ ${rc} -eq 0 ]]; then
    dataset_root="${ECODA_RUN_ROOT}/datasets/${ds}"
    rm -rf "${dataset_root}/chunks" "${dataset_root}/union" "${dataset_root}/annotations/annotation_tmp" ||
      rc=1
  fi
  if [[ ${rc} -eq 0 ]]; then
    printf 'OK\n' > "${sync_dir}/${ds}.status" || rc=1
  fi
  rm -f "${transfer_manifest}" "${verify_manifest}"
  rmdir "${lock}" 2>/dev/null || rc=1
  return "${rc}"
}

if [[ -n "${SYNC_ONLY_RUN}" ]]; then
  ecoda_open_run "${SYNC_ONLY_RUN}" || exit 1
  RUN_ID="${SYNC_ONLY_RUN}"
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  ecoda_validate_run_owned_path "${MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage4_abort "Stage 4 selection manifest is not run-owned"
  ecoda_validate_manifest "${MANIFEST}" 2 ||
    stage4_abort "Stage 4 selection manifest is invalid"
  ecoda_validate_checksum "${MANIFEST}" ||
    stage4_abort "Stage 4 selection checksum is invalid"
  if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
    ecoda_validate_exact_batch_selection "${MANIFEST}" 2 ||
      stage4_abort "Stage 4 exact selection is invalid"
  fi
  [[ -s "${ECODA_RUN_ROOT}/status/merge_watchdog" ]] &&
    grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/merge_watchdog" ||
    stage4_abort "Stage 4 merge watchdog has not reported terminal OK"
  OWNERS_MANIFEST="${ECODA_RUN_ROOT}/manifests/owners.tsv"
  ecoda_validate_run_owned_path "${OWNERS_MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage4_abort "Stage 4 owner manifest is missing or not run-owned"
  ecoda_validate_manifest "${OWNERS_MANIFEST}" 2 ||
    stage4_abort "Stage 4 owner manifest is invalid"
  SYNC_MANIFEST="${TMPDIR:-/tmp}/ecoda_stage4_sync_${$}.tsv"
  : > "${SYNC_MANIFEST}"
  skip_seen=""
  runnable_count=0
  while IFS=$'\t' read -r ds view; do
    unsuitable="$(jq -r --arg ds "${ds}" \
      '.[$ds].not_suitable_for_auto_annotation // [] |
       (index("hitme") != null and index("scatomic") != null)' \
      "${DATASETS_JSON_FILE}")"
    if [[ "${unsuitable}" == "true" ]]; then
      skip_seen="${skip_seen} ${ds}"
    else
      printf '%s\t%s\n' "${ds}" "${view}" >> "${SYNC_MANIFEST}"
      runnable_count=$((runnable_count + 1))
    fi
  done < "${MANIFEST}"
  if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
    [[ ${runnable_count} -eq 9 && "${skip_seen}" == " Alzheimer Diabetes Parkinson" ]] ||
      { rm -f "${SYNC_MANIFEST}"; stage4_abort "Stage 4 exact sync set is not nine runnable and three excluded cohorts"; }
  fi
  owner_failed=0
  sync_datasets=""
  while IFS=$'\t' read -r ds owner; do
    expected_owner="$(ecoda_owner_dir stage4 "${ds}")"
    [[ "${owner}" == "${expected_owner}" ]] || { owner_failed=1; continue; }
    owner_state="$(ecoda_owner_state "${owner}" 2>/dev/null || true)"
    owner_run="$(ecoda_owner_run "${owner}" 2>/dev/null || true)"
    [[ "${owner_state}" != "ACTIVE" || "${owner_run}" == "${RUN_ID}" ]] || { owner_failed=1; continue; }
    [[ "${owner_state}" == "OK" ]] || { owner_failed=1; continue; }
    case " ${sync_datasets} " in *" ${ds} "*) ;; *) sync_datasets="${sync_datasets} ${ds}" ;; esac
  done < "${OWNERS_MANIFEST}"
  [[ ${owner_failed} -eq 0 ]] || { rm -f "${SYNC_MANIFEST}"; stage4_abort "Stage 4 sync ownership validation failed"; }
  failed=0
  for ds in ${sync_datasets}; do
    sync_one_dataset "${ds}" "${SYNC_MANIFEST}" || failed=1
  done
  rm -f "${SYNC_MANIFEST}"
  [[ ${failed} -eq 0 ]] || stage4_abort "Stage 4 sync-only selected sync failed"
  stage4_finalize_owner_manifest OK "Stage 4 sync-only selected sync completed" ||
    stage4_abort "failed to finalize Stage 4 owners after sync"
  ecoda_set_run_state OK "Stage 4 sync-only selected sync passed" ||
    stage4_abort "failed to write Stage 4 terminal OK state"
  echo "ANNOTATION_RUN_ID=${SYNC_ONLY_RUN}"
  exit 0
fi

DATASET_NAMES=()
if [[ -n "${REUSE_RUN_ARG}" ]]; then
  ecoda_open_run "${REUSE_RUN_ARG}" || exit 1
  RUN_ID="${REUSE_RUN_ARG}"
  stage="$(sed -n 's/^STAGE=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  [[ "${stage}" == "stage4" ]] || stage4_abort "reuse run is not a Stage 4 run: ${RUN_ID}"
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  ecoda_validate_run_owned_path "${MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage4_abort "reuse selection manifest is not run-owned"
  ecoda_validate_manifest "${MANIFEST}" 2 ||
    stage4_abort "reuse run selection manifest is invalid"
  ecoda_validate_checksum "${MANIFEST}" ||
    stage4_abort "reuse run selection checksum is invalid"
  if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
    ecoda_validate_exact_batch_selection "${MANIFEST}" 2 ||
      stage4_abort "reuse run exact selection is invalid"
  fi
  if [[ -n "${SELECTION_FILE_ARG}" ]]; then
    selection_tmp="${ECODA_RUN_ROOT}/manifests/reuse-selection.check.$$"
    cp "${SELECTION_FILE_ARG}" "${selection_tmp}" ||
      stage4_abort "failed to copy supplied reuse selection"
    if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
      ecoda_validate_exact_batch_selection "${selection_tmp}" 2 ||
        stage4_abort "supplied reuse exact selection is invalid"
    else
      ecoda_validate_manifest "${selection_tmp}" 2 ||
        stage4_abort "supplied reuse selection is invalid"
    fi
    cmp -s "${selection_tmp}" "${MANIFEST}" ||
      stage4_abort "supplied selection differs from immutable reuse selection"
    rm -f "${selection_tmp}"
  fi
  validate_reuse_preparation ||
    stage4_abort "reuse preparation/chunk artifacts failed validation"
else
  if [[ -n "${SELECTION_FILE_ARG}" ]]; then
    [[ -r "${SELECTION_FILE_ARG}" ]] || { echo "ERROR: unreadable selection file." >&2; exit 1; }
  elif [[ -n "${DATASETS_ARG}" ]]; then
    ecoda_split_csv "${DATASETS_ARG}"
    DATASET_NAMES=("${ECODA_ARRAY[@]}")
    ecoda_assert_unique_items "${DATASET_NAMES[@]}"
  else
    while IFS= read -r ds; do DATASET_NAMES+=("${ds}"); done < <(
      jq -r 'keys[] | select(startswith("_") | not)' "${DATASETS_JSON_FILE}"
    )
  fi
  if [[ -z "${SELECTION_FILE_ARG}" ]]; then
    for ds in "${DATASET_NAMES[@]}"; do
      ecoda_dataset_exists "${ds}" || { echo "ERROR: unknown dataset '${ds}'." >&2; exit 1; }
    done
  fi

  RUN_ID="${ECODA_RUN_ID:-$(ecoda_new_run_id stage4)}"
  ecoda_init_run stage4 "${RUN_ID}" >/dev/null
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  MANIFEST_TMP="${MANIFEST}.build.$$"
  : > "${MANIFEST_TMP}"
  if [[ -n "${SELECTION_FILE_ARG}" ]]; then
    cp "${SELECTION_FILE_ARG}" "${MANIFEST_TMP}"
    ecoda_validate_manifest "${MANIFEST_TMP}" 2 || {
      ecoda_set_run_state FAIL "malformed Stage 4 selection"
      exit 1
    }
  else
    for ds in "${DATASET_NAMES[@]}"; do
      resolve_dataset_views "${ds}" >> "${MANIFEST_TMP}" || {
        ecoda_set_run_state FAIL "invalid Stage 4 view selection"
        exit 1
      }
    done
    ecoda_validate_manifest "${MANIFEST_TMP}" 2 || {
      ecoda_set_run_state FAIL "no Stage 4 selected views"
      exit 1
    }
  fi
  ecoda_atomic_install_manifest "${MANIFEST_TMP}" "${MANIFEST}" 2 || {
    ecoda_set_run_state FAIL "failed to install Stage 4 selection atomically"
    exit 1
  }
  rm -f "${MANIFEST_TMP}"
ecoda_write_checksum "${MANIFEST}"
fi
echo "ANNOTATION_RUN_ID=${RUN_ID}"
echo "ANNOTATION_SELECTION_MANIFEST=${MANIFEST}"
ecoda_validate_manifest "${MANIFEST}" 2 || stage4_abort "Stage 4 selection manifest is invalid"
ecoda_validate_checksum "${MANIFEST}" || stage4_abort "Stage 4 selection checksum is invalid"
if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
  ecoda_validate_exact_batch_selection "${MANIFEST}" 2 ||
    stage4_abort "Stage 4 exact selection is invalid"
fi
FULL_SELECTION_SEEN=""
while IFS=$'\t' read -r row_ds row_view; do
  ecoda_dataset_exists "${row_ds}" || stage4_abort "unknown Stage 4 dataset: ${row_ds}"
  ecoda_view_exists "${row_ds}" "${row_view}" ||
    stage4_abort "undeclared Stage 4 view: ${row_ds}/${row_view}"
  [[ -n "$(ecoda_view_input_name "${row_ds}" "${row_view}")" &&
     -n "$(ecoda_view_output_name "${row_ds}" "${row_view}")" ]] ||
    stage4_abort "missing Stage 4 input/output contract: ${row_ds}/${row_view}"
  full_row="${row_ds}/${row_view}"
  case " ${FULL_SELECTION_SEEN} " in
    *" ${full_row} "*) stage4_abort "duplicate Stage 4 selection row ${full_row}" ;;
  esac
  FULL_SELECTION_SEEN="${FULL_SELECTION_SEEN} ${full_row}"
done < "${MANIFEST}"

# Remove explicitly skipped cohorts from the runnable manifest and reject
# duplicate dataset/view owners.
RUNNABLE_MANIFEST="${ECODA_RUN_ROOT}/manifests/runnable_selection.tsv"
RUNNABLE_TMP="${RUNNABLE_MANIFEST}.build.$$"
OWNERS_MANIFEST="${ECODA_RUN_ROOT}/manifests/owners.tsv"
: > "${RUNNABLE_TMP}"
SEEN=""

# Derive exclusions once from the immutable selection, before the runnable
# manifest is constructed. Duplicate dataset/view rows are rejected while
# legitimate multi-view rows still produce one skip record per dataset.
SKIP_FILE="${ECODA_RUN_ROOT}/status/skipped"
SKIP_TMP="${SKIP_FILE}.build.$$"
SKIP_SEEN=""
SKIP_ROWS=""
SKIP_COUNT=0
while IFS=$'\t' read -r row_ds row_view; do
  unsuitable="$(jq -r --arg ds "${row_ds}" \
    '.[$ds].not_suitable_for_auto_annotation // [] |
     (index("hitme") != null and index("scatomic") != null)' \
    "${DATASETS_JSON_FILE}")"
  [[ "${unsuitable}" == "true" ]] || continue
  case " ${SKIP_SEEN} " in
    *" ${row_ds} "*) continue ;;
  esac
  SKIP_SEEN="${SKIP_SEEN} ${row_ds}"
  SKIP_COUNT=$((SKIP_COUNT + 1))
  printf '%s\tSKIP_NOT_SUITABLE\n' "${row_ds}" >> "${SKIP_TMP}"
done < "${MANIFEST}"
if [[ ${SKIP_COUNT} -gt 0 ]]; then
  ecoda_atomic_install_manifest "${SKIP_TMP}" "${SKIP_FILE}" 2 ||
    stage4_abort "failed to install Stage 4 skip records"
else
  rm -f "${SKIP_TMP}" "${SKIP_FILE}"
fi
if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
  [[ ${SKIP_COUNT} -eq 3 ]] ||
    stage4_abort "exact batch selection did not yield three exclusions"
  [[ "${SKIP_SEEN}" == " Alzheimer Diabetes Parkinson" ]] ||
    stage4_abort "unexpected Stage 4 exclusion set"
fi

while IFS=$'\t' read -r ds view; do
  unsuitable="$(jq -r --arg ds "${ds}" \
    '.[$ds].not_suitable_for_auto_annotation // [] |
     (index("hitme") != null and index("scatomic") != null)' \
    "${DATASETS_JSON_FILE}")"
  if [[ "${unsuitable}" == "true" ]]; then
    continue
  fi
  row="${ds}/${view}"
  case " ${SEEN} " in
    *" ${row} "*) stage4_abort "duplicate Stage 4 row ${row}" ;;
  esac
  SEEN="${SEEN} ${row}"
  input_path="$(output_path_for "${ds}" "${view}")" ||
    stage4_abort "missing Stage 4 output contract ${row}"
  if [[ "${ANNOTATION_SUBMITTER_TEST:-0}" != "1" ]]; then
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
      --path "${input_path}" --view "${view}" --method "Stage 4 annotation" >/dev/null 2>&1 ||
      stage4_abort "invalid Stage 4 h5ad ${row}"
    ecoda_validate_checksum "${input_path}" ||
      stage4_abort "invalid Stage 4 checksum ${row}"
  fi
  if [[ ${FORCE_ARG} -eq 1 ]]; then
    ecoda_invalidate_artifact "${ECODA_RUN_ROOT}/datasets/${ds}/merge.ok" ||
      stage4_abort "failed to invalidate run-owned merge marker ${row}"
  fi
  printf '%s\t%s\n' "${ds}" "${view}" >> "${RUNNABLE_TMP}"
done < "${MANIFEST}"
if [[ ! -s "${RUNNABLE_TMP}" ]]; then
  ecoda_atomic_write "${RUNNABLE_MANIFEST}" "" ||
    stage4_abort "failed to create empty Stage 4 runnable manifest"
  ecoda_atomic_write "${OWNERS_MANIFEST}" "" ||
    stage4_abort "failed to create empty Stage 4 owner manifest"
  ecoda_atomic_write "${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv" "" ||
    stage4_abort "failed to create empty Stage 4 scheduler manifest"
  rm -f "${RUNNABLE_TMP}"
  ecoda_set_run_state OK "all selected datasets explicitly exempt from auto-annotation" ||
    stage4_abort "failed to write Stage 4 terminal OK state"
  exit 0
fi
if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
  runnable_rows="$(wc -l < "${RUNNABLE_TMP}" | tr -d '[:space:]')"
  [[ "${runnable_rows}" == "9" ]] ||
    stage4_abort "exact batch selection did not yield nine runnable rows"
fi
ecoda_atomic_install_manifest "${RUNNABLE_TMP}" "${RUNNABLE_MANIFEST}" 2 ||
  stage4_abort "failed to install Stage 4 runnable manifest atomically"
rm -f "${RUNNABLE_TMP}"
WORK_MANIFEST="${RUNNABLE_MANIFEST}"
if [[ ${FORCE_ARG} -eq 0 && "${ANNOTATION_SUBMITTER_TEST:-0}" != "1" ]]; then
  VALIDATION_MANIFEST="${ECODA_RUN_ROOT}/manifests/runnable_validated.tsv"
  VALIDATION_TMP="${VALIDATION_MANIFEST}.build.$$"
  : > "${VALIDATION_TMP}"
  SCAN_DATASETS=""
  while IFS=$'\t' read -r ds view; do
    case " ${SCAN_DATASETS} " in *" ${ds} "*) ;; *) SCAN_DATASETS="${SCAN_DATASETS} ${ds}" ;; esac
  done < "${WORK_MANIFEST}"
  for ds in ${SCAN_DATASETS}; do
    dataset_valid=1
    while IFS=$'\t' read -r row_ds view; do
      [[ "${row_ds}" == "${ds}" ]] || continue
      input_path="$(output_path_for "${ds}" "${view}")"
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
        --h5ad "${input_path}" --require-sidecar >/dev/null 2>&1 || dataset_valid=0
      ecoda_validate_checksum "${input_path}" || dataset_valid=0
    done < "${WORK_MANIFEST}"
    if [[ ${dataset_valid} -eq 1 ]]; then
      printf '%s\tALREADY_ANNOTATED\n' "${ds}" >> "${ECODA_RUN_ROOT}/status/skipped"
    else
      while IFS=$'\t' read -r row_ds view; do
        [[ "${row_ds}" == "${ds}" ]] && printf '%s\t%s\n' "${row_ds}" "${view}" >> "${VALIDATION_TMP}"
      done < "${WORK_MANIFEST}"
    fi
  done
  if [[ -s "${VALIDATION_TMP}" ]]; then
    ecoda_atomic_install_manifest "${VALIDATION_TMP}" "${VALIDATION_MANIFEST}" 2 || {
      ecoda_set_run_state FAIL "failed to install Stage 4 validated manifest atomically"
      exit 1
    }
    WORK_MANIFEST="${VALIDATION_MANIFEST}"
  else
    WORK_MANIFEST=""
  fi
  rm -f "${VALIDATION_TMP}"
fi
if [[ ! -s "${WORK_MANIFEST}" ]]; then
  ecoda_set_run_state OK "all selected annotation artifacts already validated" ||
    stage4_abort "failed to write Stage 4 terminal OK state"
  exit 0
fi

PREP_MANIFEST="${ECODA_RUN_ROOT}/manifests/preparation.tsv"
OWNERS_MANIFEST="${ECODA_RUN_ROOT}/manifests/owners.tsv"
ecoda_owner_clear_tracked
if [[ ${SKIP_PREPARE} -eq 1 ]]; then
  ecoda_validate_run_owned_path "${PREP_MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage4_abort "reuse preparation manifest is not run-owned"
  ecoda_validate_manifest "${PREP_MANIFEST}" 3 ||
    stage4_abort "reuse preparation manifest is invalid"
  OWNERS_TMP="${OWNERS_MANIFEST}.build.$$"
  : > "${OWNERS_TMP}"
  DATASET_GROUPS=""
  while IFS=$'\t' read -r ds view; do
    case " ${DATASET_GROUPS} " in *" ${ds} "*) ;; *) DATASET_GROUPS="${DATASET_GROUPS} ${ds}" ;; esac
  done < "${WORK_MANIFEST}"
  for ds in ${DATASET_GROUPS}; do
    set +e
    owner_dir="$(ecoda_owner_acquire stage4 "${ds}" "${RUN_ID}" "${FORCE_ARG}" 0)"
    owner_rc=$?
    set -e
    [[ ${owner_rc} -eq 0 ]] || stage4_abort "ownership conflict for ${ds}"
    ecoda_owner_track "${owner_dir}" || stage4_abort "failed to track owner for ${ds}"
    printf '%s\t%s\n' "${ds}" "${owner_dir}" >> "${OWNERS_TMP}"
  done
  ecoda_atomic_install_manifest "${OWNERS_TMP}" "${OWNERS_MANIFEST}" 2 ||
    stage4_abort "failed to install Stage 4 reuse owner manifest atomically"
  rm -f "${OWNERS_TMP}"
else
  PREP_TMP="${PREP_MANIFEST}.build.$$"
  OWNERS_TMP="${OWNERS_MANIFEST}.build.$$"
  : > "${PREP_TMP}"
  : > "${OWNERS_TMP}"
  DATASET_GROUPS=""
  while IFS=$'\t' read -r ds view; do
    case " ${DATASET_GROUPS} " in *" ${ds} "*) ;; *) DATASET_GROUPS="${DATASET_GROUPS} ${ds}" ;; esac
  done < "${WORK_MANIFEST}"
  for ds in ${DATASET_GROUPS}; do
    views=""
    while IFS=$'\t' read -r row_ds view; do
      [[ "${row_ds}" == "${ds}" ]] || continue
      [[ -n "${views}" ]] && views="${views},"
      views="${views}${view}"
    done < "${WORK_MANIFEST}"
    set +e
    owner_dir="$(ecoda_owner_acquire stage4 "${ds}" "${RUN_ID}" "${FORCE_ARG}" 0)"
    owner_rc=$?
    set -e
    [[ ${owner_rc} -eq 0 ]] || stage4_abort "ownership conflict for ${ds}"
    ecoda_owner_track "${owner_dir}" || stage4_abort "failed to track owner for ${ds}"
    printf '%s\t%s\t%s\n' "${ds}" "${views}" "${ECODA_RUN_ROOT}" >> "${PREP_TMP}"
    printf '%s\t%s\n' "${ds}" "${owner_dir}" >> "${OWNERS_TMP}"
  done
  ecoda_atomic_install_manifest "${PREP_TMP}" "${PREP_MANIFEST}" 3 ||
    stage4_abort "failed to install Stage 4 preparation manifest atomically"
  ecoda_atomic_install_manifest "${OWNERS_TMP}" "${OWNERS_MANIFEST}" 2 ||
    stage4_abort "failed to install Stage 4 owner manifest atomically"
  rm -f "${PREP_TMP}" "${OWNERS_TMP}"
fi
if [[ -z "${DATASET_GROUPS:-}" ]]; then
  DATASET_GROUPS=""
  while IFS=$'\t' read -r ds view; do
    case " ${DATASET_GROUPS} " in *" ${ds} "*) ;; *) DATASET_GROUPS="${DATASET_GROUPS} ${ds}" ;; esac
  done < "${WORK_MANIFEST}"
fi
SCHEDULER_IDS_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
if [[ ! -e "${SCHEDULER_IDS_FILE}" ]]; then
  ecoda_atomic_write "${SCHEDULER_IDS_FILE}" "" ||
    stage4_abort "failed to initialize Stage 4 scheduler manifest"
fi

mkdir -p "${LOGS_DIR}"
emit_watchdog_status_ids() {
  local kind="$1" status_file="$2" status_line
  [[ -s "${status_file}" ]] || {
    echo "ERROR: ${kind} watchdog status is missing: ${status_file}" >&2
    return 1
  }
  while IFS= read -r status_line; do
    case "${status_line}" in
      SCHEDULER_ID=*|ARRAY_JOB_ID=*)
        printf 'ANNOTATION_SCHEDULER_ID=%s:%s\n' "${kind}" "${status_line#*=}"
        ;;
    esac
  done < "${status_file}"
}
if [[ "${ANNOTATION_SUBMITTER_TEST:-0}" != "1" ]]; then
  if ! "${SCRIPT_DIR}/1.0_stage_reference_maps.sh"; then
    stage4_abort "Stage 4 reference-map staging failed"
  fi
  cache_valid=0
  if ${PIXI_RSCRIPT} "${SCRIPT_DIR}/2.0_create_scgate_db.R" --validate-only >/dev/null 2>&1; then
    cache_valid=1
  fi
  if [[ -e "${SCGATE_DB_PATH}" && ${cache_valid} -eq 0 && ${FORCE_ARG} -eq 0 ]]; then
    stage4_abort "existing scGate cache failed validation; use --force to rebuild"
  fi
  if [[ ! -e "${SCGATE_DB_PATH}" || ${FORCE_ARG} -eq 1 ]]; then
    SCGATE_FORCE_ARG=()
    [[ ${FORCE_ARG} -eq 1 ]] && SCGATE_FORCE_ARG=(--force)
    set +e
    scgate_msg="$(sbatch --parsable --wait --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=4G \
      --time="${ANNOTATION_WATCHDOG_TIME_LIMIT:-12:00:00}" --output="${LOGS_DIR}/4_scgate_%j.log" \
      --error="${LOGS_DIR}/4_scgate_%j.err" --mail-user="${USER_EMAIL}" \
      --export="ALL" --wrap="${PIXI_RSCRIPT} ${SCRIPT_DIR}/2.0_create_scgate_db.R ${SCGATE_FORCE_ARG[*]}")"
    scgate_rc=$?
    set -e
    scgate_job="${scgate_msg%%;*}"
    if [[ "${scgate_job}" =~ ^[0-9]+$ ]]; then
      stage4_record_scheduler SCGATE "${scgate_job}" ||
        stage4_abort "failed to persist scGate scheduler ID"
    fi
    [[ ${scgate_rc} -eq 0 && "${scgate_job}" =~ ^[0-9]+$ && -s "${SCGATE_DB_PATH}" ]] ||
      stage4_abort "scGate cache preflight failed"
    ${PIXI_RSCRIPT} "${SCRIPT_DIR}/2.0_create_scgate_db.R" --validate-only >/dev/null 2>&1 ||
      stage4_abort "scGate cache failed post-build validation"
    echo "ANNOTATION_SCGATE_JOB_ID=${scgate_job}"
  fi
fi
if [[ ${SKIP_PREPARE} -eq 0 ]]; then
  prep_count="$(wc -l < "${PREP_MANIFEST}" | tr -d '[:space:]')"
  [[ "${prep_count}" =~ ^[1-9][0-9]*$ ]] || stage4_abort "Stage 4 preparation manifest is empty"
  export FORCE_ANNOTATION="${FORCE_ARG}" ANNOTATION_RUN_ID="${RUN_ID}"
  set +e
  prep_msg="$(sbatch --parsable --array="1-${prep_count}%${THROTTLE}" --mem="${MEMORY}" --partition="${PARTITION}" \
    --output="${LOGS_DIR}/4_annotation_prepare_%A_%a.log" --error="${LOGS_DIR}/4_annotation_prepare_%A_%a.err" \
    --mail-user="${USER_EMAIL}" --export="ALL,ANNOTATION_PREP_MANIFEST=${PREP_MANIFEST},ANNOTATION_RUN_ID=${RUN_ID},FORCE_ANNOTATION=${FORCE_ARG}" \
    "${SCRIPT_DIR}/1.2_prepare_chunks_worker.sh")"
  prep_rc=$?
  set -e
  PREP_ARRAY="${prep_msg%%;*}"
  if [[ "${PREP_ARRAY}" =~ ^[0-9]+$ ]]; then
    stage4_record_scheduler ARRAY "${PREP_ARRAY}" ||
      stage4_abort "failed to persist Stage 4 preparation array ID"
  fi
  [[ ${prep_rc} -eq 0 && "${PREP_ARRAY}" =~ ^[0-9]+$ ]] ||
    stage4_abort "Stage 4 preparation array submission failed"
  echo "ANNOTATION_PREP_ARRAY_JOB_ID=${PREP_ARRAY}"
  set +e
  prep_watchdog_msg="$(sbatch --parsable --wait --dependency="afterany:${PREP_ARRAY}" --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=2G \
    --time="${ANNOTATION_WATCHDOG_TIME_LIMIT:-12:00:00}" --output="${LOGS_DIR}/4_annotation_prepare_watchdog_%j.log" \
    --error="${LOGS_DIR}/4_annotation_prepare_watchdog_%j.err" --mail-user="${USER_EMAIL}" \
    --export="ALL,ANNOTATION_RUN_ID=${RUN_ID}" "${SCRIPT_DIR}/1.3_prepare_chunks_watchdog.sh" "${RUN_ID}" "${PREP_MANIFEST}" \
    "${PREP_ARRAY}" "${MEMORY}" "${MAX_MEMORY}" "${PARTITION}" "${THROTTLE}")"
  prep_watchdog_rc=$?
  set -e
  PREP_WATCHDOG="${prep_watchdog_msg%%;*}"
  if [[ "${PREP_WATCHDOG}" =~ ^[0-9]+$ ]]; then
    stage4_record_scheduler WATCHDOG "${PREP_WATCHDOG}" ||
      stage4_abort "failed to persist Stage 4 preparation watchdog ID"
  fi
  [[ ${prep_watchdog_rc} -eq 0 && "${PREP_WATCHDOG}" =~ ^[0-9]+$ ]] ||
    stage4_abort "Stage 4 preparation watchdog submission or execution failed"
  echo "ANNOTATION_PREP_WATCHDOG_JOB_ID=${PREP_WATCHDOG}"
  if [[ "${ANNOTATION_SUBMITTER_TEST:-0}" != "1" ]]; then
    emit_watchdog_status_ids preparation "${ECODA_RUN_ROOT}/status/preparation_watchdog" ||
      stage4_abort "Stage 4 preparation watchdog status missing"
  fi
fi

CHUNK_MANIFEST="${ECODA_RUN_ROOT}/manifests/chunks.tsv"
if [[ ${SKIP_PREPARE} -eq 1 ]]; then
  ecoda_validate_run_owned_path "${CHUNK_MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage4_abort "reuse chunk manifest is not run-owned"
  ecoda_validate_manifest "${CHUNK_MANIFEST}" 3 ||
    stage4_abort "reuse chunk manifest is invalid"
  TOTAL_CHUNKS="$(wc -l < "${CHUNK_MANIFEST}" | tr -d '[:space:]')"
else
  CHUNK_TMP="${CHUNK_MANIFEST}.build.$$"
  : > "${CHUNK_TMP}"
  TOTAL_CHUNKS=0
  for ds in ${DATASET_GROUPS}; do
    chunk_dir="${ECODA_RUN_ROOT}/datasets/${ds}/chunks"
    feather_dir="${ECODA_RUN_ROOT}/datasets/${ds}/annotations"
    shopt -s nullglob
    chunks=("${chunk_dir}"/chunk_*.txt)
    shopt -u nullglob
    for chunk in "${chunks[@]-}"; do
      printf '%s\t%s\t%s\n' "${ds}" "${chunk}" "${feather_dir}" >> "${CHUNK_TMP}"
      TOTAL_CHUNKS=$((TOTAL_CHUNKS + 1))
    done
  done
  if [[ ${TOTAL_CHUNKS} -eq 0 && "${ANNOTATION_SUBMITTER_TEST:-0}" == "1" ]]; then
    for ds in ${DATASET_GROUPS}; do
      mkdir -p "${ECODA_RUN_ROOT}/datasets/${ds}/chunks" "${ECODA_RUN_ROOT}/datasets/${ds}/annotations"
      fake_chunk="${ECODA_RUN_ROOT}/datasets/${ds}/chunks/chunk_1.txt"
      printf '%s\nfake_sample\n' "${ECODA_RUN_ROOT}/datasets/${ds}/union/union.h5ad" > "${fake_chunk}"
      printf '%s\t%s\t%s\n' "${ds}" "${fake_chunk}" "${ECODA_RUN_ROOT}/datasets/${ds}/annotations" >> "${CHUNK_TMP}"
      TOTAL_CHUNKS=$((TOTAL_CHUNKS + 1))
    done
  fi
  [[ ${TOTAL_CHUNKS} -gt 0 ]] || stage4_abort "no annotation chunks prepared"
  ecoda_atomic_install_manifest "${CHUNK_TMP}" "${CHUNK_MANIFEST}" 3 ||
    stage4_abort "failed to install Stage 4 chunk manifest atomically"
  rm -f "${CHUNK_TMP}"
fi
[[ ${TOTAL_CHUNKS} -gt 0 ]] || stage4_abort "reuse chunk manifest is empty"
set +e
annot_msg="$(sbatch --parsable --array="1-${TOTAL_CHUNKS}%${THROTTLE}" --mem="${MEMORY}" --partition="${PARTITION}" \
  --output="${LOGS_DIR}/4_cell_type_annotation_%A_%a.log" --error="${LOGS_DIR}/4_cell_type_annotation_%A_%a.err" \
  --mail-user="${USER_EMAIL}" --export="ALL,CHUNKS_MANIFEST=${CHUNK_MANIFEST},ANNOTATION_RUN_ID=${RUN_ID},ANNOTATION_ERROR_PREFIX=${LOGS_DIR}/4_cell_type_annotation" \
  "${SCRIPT_DIR}/2.1_run_worker.sh")"
annot_rc=$?
set -e
ANNOT_ARRAY="${annot_msg%%;*}"
if [[ "${ANNOT_ARRAY}" =~ ^[0-9]+$ ]]; then
  stage4_record_scheduler ARRAY "${ANNOT_ARRAY}" ||
    stage4_abort "failed to persist Stage 4 annotation array ID"
fi
[[ ${annot_rc} -eq 0 && "${ANNOT_ARRAY}" =~ ^[0-9]+$ ]] ||
  stage4_abort "Stage 4 annotation array submission failed"
echo "ANNOTATION_ARRAY_JOB_ID=${ANNOT_ARRAY}"
set +e
annot_watchdog_msg="$(sbatch --parsable --wait --dependency="afterany:${ANNOT_ARRAY}" --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=2G \
  --time="${ANNOTATION_WATCHDOG_TIME_LIMIT:-12:00:00}" --output="${LOGS_DIR}/4_cell_type_annotation_watchdog_%j.log" \
  --error="${LOGS_DIR}/4_cell_type_annotation_watchdog_%j.err" --mail-user="${USER_EMAIL}" --export="ALL,ANNOTATION_RUN_ID=${RUN_ID}" \
  "${SCRIPT_DIR}/1.2_annotation_watchdog.sh" "${RUN_ID}" "${CHUNK_MANIFEST}" "${ANNOT_ARRAY}" \
  "${MEMORY}" "${MAX_MEMORY}" "${PARTITION}" "${THROTTLE}")"
annot_watchdog_rc=$?
set -e
ANNOT_WATCHDOG="${annot_watchdog_msg%%;*}"
if [[ "${ANNOT_WATCHDOG}" =~ ^[0-9]+$ ]]; then
  stage4_record_scheduler WATCHDOG "${ANNOT_WATCHDOG}" ||
    stage4_abort "failed to persist Stage 4 annotation watchdog ID"
fi
[[ ${annot_watchdog_rc} -eq 0 && "${ANNOT_WATCHDOG}" =~ ^[0-9]+$ ]] ||
  stage4_abort "Stage 4 annotation watchdog submission or execution failed"
echo "ANNOTATION_WATCHDOG_JOB_ID=${ANNOT_WATCHDOG}"
if [[ "${ANNOTATION_SUBMITTER_TEST:-0}" != "1" ]]; then
  emit_watchdog_status_ids annotation "${ECODA_RUN_ROOT}/status/annotation_watchdog" ||
    stage4_abort "Stage 4 annotation watchdog status missing"
fi
MERGE_MANIFEST="${ECODA_RUN_ROOT}/manifests/merge.tsv"
MERGE_TMP="${MERGE_MANIFEST}.build.$$"
: > "${MERGE_TMP}"
for ds in ${DATASET_GROUPS}; do
  views=""
  while IFS=$'\t' read -r row_ds view; do
    [[ "${row_ds}" == "${ds}" ]] || continue
    [[ -n "${views}" ]] && views="${views},"
    views="${views}${view}"
  done < "${WORK_MANIFEST}"
  printf '%s\t%s\t%s\n' "${ds}" "${views}" "${ECODA_RUN_ROOT}" >> "${MERGE_TMP}"
done
MERGE_COUNT="$(wc -l < "${MERGE_TMP}" | tr -d '[:space:]')"
[[ "${MERGE_COUNT}" =~ ^[1-9][0-9]*$ ]] || stage4_abort "Stage 4 merge manifest is empty"
ecoda_atomic_install_manifest "${MERGE_TMP}" "${MERGE_MANIFEST}" 3 ||
  stage4_abort "failed to install Stage 4 merge manifest atomically"
rm -f "${MERGE_TMP}"
set +e
merge_msg="$(sbatch --parsable --array="1-${MERGE_COUNT}%${THROTTLE}" --mem="${MEMORY}" --partition="${PARTITION}" \
  --output="${LOGS_DIR}/4_annotation_merge_%A_%a.log" --error="${LOGS_DIR}/4_annotation_merge_%A_%a.err" \
  --mail-user="${USER_EMAIL}" --export="ALL,ANNOTATION_MERGE_MANIFEST=${MERGE_MANIFEST},ANNOTATION_RUN_ID=${RUN_ID},FORCE_ANNOTATION=${FORCE_ARG}" \
  "${SCRIPT_DIR}/3.2_merge_worker.sh")"
merge_rc=$?
set -e
MERGE_ARRAY="${merge_msg%%;*}"
if [[ "${MERGE_ARRAY}" =~ ^[0-9]+$ ]]; then
  stage4_record_scheduler ARRAY "${MERGE_ARRAY}" ||
    stage4_abort "failed to persist Stage 4 merge array ID"
fi
[[ ${merge_rc} -eq 0 && "${MERGE_ARRAY}" =~ ^[0-9]+$ ]] ||
  stage4_abort "Stage 4 merge array submission failed"
echo "ANNOTATION_MERGE_ARRAY_JOB_ID=${MERGE_ARRAY}"
set +e
merge_watchdog_msg="$(sbatch --parsable --wait --dependency="afterany:${MERGE_ARRAY}" --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=2G \
  --time="${ANNOTATION_WATCHDOG_TIME_LIMIT:-12:00:00}" --output="${LOGS_DIR}/4_annotation_merge_watchdog_%j.log" \
  --error="${LOGS_DIR}/4_annotation_merge_watchdog_%j.err" --mail-user="${USER_EMAIL}" --export="ALL,ANNOTATION_RUN_ID=${RUN_ID}" \
  "${SCRIPT_DIR}/3.3_merge_watchdog.sh" "${RUN_ID}" "${MERGE_MANIFEST}" "${MERGE_ARRAY}" \
  "${MEMORY}" "${MAX_MEMORY}" "${PARTITION}" "${THROTTLE}")"
merge_watchdog_rc=$?
set -e
MERGE_WATCHDOG="${merge_watchdog_msg%%;*}"
if [[ "${MERGE_WATCHDOG}" =~ ^[0-9]+$ ]]; then
  stage4_record_scheduler WATCHDOG "${MERGE_WATCHDOG}" ||
    stage4_abort "failed to persist Stage 4 merge watchdog ID"
fi
[[ ${merge_watchdog_rc} -eq 0 && "${MERGE_WATCHDOG}" =~ ^[0-9]+$ ]] ||
  stage4_abort "Stage 4 merge watchdog submission or execution failed"
echo "ANNOTATION_MERGE_WATCHDOG_JOB_ID=${MERGE_WATCHDOG}"
if [[ "${ANNOTATION_SUBMITTER_TEST:-0}" != "1" ]]; then
  emit_watchdog_status_ids merge "${ECODA_RUN_ROOT}/status/merge_watchdog" ||
    stage4_abort "Stage 4 merge watchdog status missing"
fi
SCHEDULER_IDS_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
[[ -r "${SCHEDULER_IDS_FILE}" ]] ||
  stage4_abort "Stage 4 scheduler manifest is missing"
if [[ "${ANNOTATION_SUBMITTER_TEST:-0}" != "1" ]]; then
  COMPLETE_SCHEDULER_TMP="${SCHEDULER_IDS_FILE}.complete.build.$$"
  cp "${SCHEDULER_IDS_FILE}" "${COMPLETE_SCHEDULER_TMP}" ||
    stage4_abort "failed to copy complete Stage 4 scheduler manifest"
  SCHEDULER_ID_SEEN=""
  while IFS=$'\t' read -r scheduler_kind scheduler_id; do
    [[ "${scheduler_kind}" =~ ^(SCGATE|ARRAY|WATCHDOG|STATUS)$ &&
       "${scheduler_id}" =~ ^[0-9]+$ ]] ||
      stage4_abort "invalid Stage 4 scheduler manifest row"
    case " ${SCHEDULER_ID_SEEN} " in
      *" ${scheduler_id} "*) stage4_abort "duplicate Stage 4 scheduler ID" ;;
    esac
    SCHEDULER_ID_SEEN="${SCHEDULER_ID_SEEN} ${scheduler_id}"
  done < "${SCHEDULER_IDS_FILE}"
  STATUS_FILES=(
    "${ECODA_RUN_ROOT}/status/annotation_watchdog"
    "${ECODA_RUN_ROOT}/status/merge_watchdog"
  )
  [[ -n "${PREP_ARRAY:-}" ]] &&
    STATUS_FILES=("${ECODA_RUN_ROOT}/status/preparation_watchdog" "${STATUS_FILES[@]}")
  for status_file in "${STATUS_FILES[@]}"; do
    while IFS= read -r status_line; do
      case "${status_line}" in
        SCHEDULER_ID=*)
          scheduler_id="${status_line#*=}"
          [[ "${scheduler_id}" =~ ^[0-9]+$ ]] ||
            stage4_abort "invalid Stage 4 watchdog scheduler ID"
          case " ${SCHEDULER_ID_SEEN} " in
            *" ${scheduler_id} "*) ;;
            *)
              printf 'STATUS\t%s\n' "${scheduler_id}" >> "${COMPLETE_SCHEDULER_TMP}"
              SCHEDULER_ID_SEEN="${SCHEDULER_ID_SEEN} ${scheduler_id}"
              ;;
          esac
          ;;
      esac
    done < "${status_file}"
  done
  if ! ecoda_atomic_install_manifest "${COMPLETE_SCHEDULER_TMP}" "${SCHEDULER_IDS_FILE}" 2; then
    stage4_abort "failed to install complete Stage 4 scheduler ID manifest"
  fi
  rm -f "${COMPLETE_SCHEDULER_TMP}"
else
  while IFS=$'\t' read -r scheduler_kind scheduler_id; do
    [[ "${scheduler_kind}" =~ ^(SCGATE|ARRAY|WATCHDOG|STATUS)$ &&
       "${scheduler_id}" =~ ^[0-9]+$ ]] ||
      stage4_abort "invalid Stage 4 scheduler manifest row"
  done < "${SCHEDULER_IDS_FILE}"
fi
if [[ "${ANNOTATION_SUBMITTER_TEST:-0}" == "1" ]]; then
  ecoda_set_run_state OK "submitter test mode; dataset arrays submitted" ||
    stage4_abort "failed to write Stage 4 terminal OK state"
  exit 0
fi
[[ -d "${NAS_TARGET_DIR}" ]] || stage4_abort "NAS unreachable"
mkdir -p "${ECODA_RUN_ROOT}/status/sync" ||
  stage4_abort "failed to create Stage 4 sync status directory"
SYNC_PIDS=()
for ds in ${DATASET_GROUPS}; do
  sync_one_dataset "${ds}" "${WORK_MANIFEST}" &
  SYNC_PIDS+=("$!")
done
sync_failed=0
for pid in "${SYNC_PIDS[@]}"; do
  wait "${pid}" || sync_failed=1
done
[[ ${sync_failed} -eq 0 ]] || stage4_abort "Stage 4 selected NAS sync failed"
stage4_finalize_owner_manifest OK "Stage 4 merge and selected sync completed" ||
  stage4_abort "failed to finalize Stage 4 owners"
ecoda_set_run_state OK "Stage 4 preparation, annotation, merge, and selected sync completed" ||
  stage4_abort "failed to write Stage 4 terminal OK state"
