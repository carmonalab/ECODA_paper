#!/bin/bash
# Canonical Pipeline 5 matrix wrapper. It submits every independent selected
# dataset/method row before waiting, then owns one aggregate gate and one final
# merge/checksum/NAS synchronization.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
export ECODA_GATE_STAGE=stage5
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/h5ad_preflight_submit.sh"
cd "${PROJECT_ROOT}"

DATASETS_ARG=""
DATASETS_SET=0
METHODS_ARG=""
METHODS_SET=0
ANALYSES_ARG=""
ANALYSES_SET=0
SELECTION_FILE_ARG=""
SELECTION_FILE_SET=0
PASS_ARG=""
PASS_SET=0
EXACT_BATCH_SELECTION=0
FORCE_ARG=0
SYNC_ONLY_RUN=""
SYNC_ONLY_SET=0
PARTITION_ARG=""
MEMORY="${BENCHMARK_MEM}"
MAX_MEMORY="${BENCHMARK_MEM_MAX}"
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"

usage() {
  cat <<'EOF'
Usage: 1_submit_hpc_array.sh [--datasets LIST] [--methods LIST]
       [--analyses trans,zeroimp] [--selection-file TSV]
       [--exact-batch-selection] [--pass uncorrected|corrected]
       [--force] [--sync-only RUN_ID]
       [--partition NAME] [--mem VALUE] [--max-mem VALUE] [--throttle N]

Selection-file rows are DATASET<TAB>VIEW<TAB>LABEL. Ordinary methods use the
benchmark_analysis view; batch mode uses the selected explicit pass view.
Exact batch mode requires the immutable twelve-row uncorrected matrix.
EOF
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets) DATASETS_ARG="${2:-}"; DATASETS_SET=1; shift 2 ;;
    --datasets=*) DATASETS_ARG="${1#*=}"; DATASETS_SET=1; shift ;;
    --ds_name) DATASETS_ARG="${2:-}"; DATASETS_SET=1; shift 2 ;;
    --ds_name=*) DATASETS_ARG="${1#*=}"; DATASETS_SET=1; shift ;;
    --methods) METHODS_ARG="${2:-}"; METHODS_SET=1; shift 2 ;;
    --methods=*) METHODS_ARG="${1#*=}"; METHODS_SET=1; shift ;;
    --analyses|--analysis) ANALYSES_ARG="${2:-}"; ANALYSES_SET=1; shift 2 ;;
    --analyses=*|--analysis=*) ANALYSES_ARG="${1#*=}"; ANALYSES_SET=1; shift ;;
    --selection-file) SELECTION_FILE_ARG="${2:-}"; SELECTION_FILE_SET=1; shift 2 ;;
    --selection-file=*) SELECTION_FILE_ARG="${1#*=}"; SELECTION_FILE_SET=1; shift ;;
    --exact-batch-selection) EXACT_BATCH_SELECTION=1; shift ;;
    --pass|--analysis-pass) PASS_ARG="${2:-}"; PASS_SET=1; shift 2 ;;
    --pass=*|--analysis-pass=*) PASS_ARG="${1#*=}"; PASS_SET=1; shift ;;
    --force) FORCE_ARG=1; shift ;;
    --sync-only) SYNC_ONLY_RUN="${2:-}"; SYNC_ONLY_SET=1; shift 2 ;;
    --sync-only=*) SYNC_ONLY_RUN="${1#*=}"; SYNC_ONLY_SET=1; shift ;;
    --partition) PARTITION_ARG="${2:-}"; shift 2 ;;
    --partition=*) PARTITION_ARG="${1#*=}"; shift ;;
    --mem) MEMORY="${2:-}"; shift 2 ;;
    --mem=*) MEMORY="${1#*=}"; shift ;;
    --max-mem) MAX_MEMORY="${2:-}"; shift 2 ;;
    --max-mem=*) MAX_MEMORY="${1#*=}"; shift ;;
    --throttle) THROTTLE="${2:-}"; shift 2 ;;
    --throttle=*) THROTTLE="${1#*=}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done
if [[ -n "${SYNC_ONLY_RUN}" && ${FORCE_ARG} -eq 1 ]]; then echo "ERROR: --sync-only cannot use --force." >&2; exit 1; fi
if [[ -n "${PASS_ARG}" && "${PASS_ARG}" != uncorrected && "${PASS_ARG}" != corrected ]]; then echo "ERROR: --pass must be uncorrected or corrected." >&2; exit 1; fi
command -v jq >/dev/null 2>&1 || { echo "ERROR: jq is required for benchmark selection." >&2; exit 1; }
if [[ ${DATASETS_SET} -eq 1 && -z "${DATASETS_ARG}" ]]; then
  echo "ERROR: --datasets must not be empty." >&2
  exit 1
fi
if [[ ${METHODS_SET} -eq 1 && -z "${METHODS_ARG}" ]]; then
  echo "ERROR: --methods must not be empty." >&2
  exit 1
fi
if [[ ${ANALYSES_SET} -eq 1 && -z "${ANALYSES_ARG}" ]]; then
  echo "ERROR: --analyses must not be empty." >&2
  exit 1
fi
if [[ ${SELECTION_FILE_SET} -eq 1 && -z "${SELECTION_FILE_ARG}" ]]; then
  echo "ERROR: --selection-file must not be empty." >&2
  exit 1
fi
if [[ ${PASS_SET} -eq 1 && -z "${PASS_ARG}" ]]; then
  echo "ERROR: --pass must not be empty." >&2
  exit 1
fi
if [[ ${SYNC_ONLY_SET} -eq 1 && -z "${SYNC_ONLY_RUN}" ]]; then
  echo "ERROR: --sync-only requires a run ID." >&2
  exit 1
fi
EXPECTED_BATCH_METHODS="prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,pilotgm,qot"
if [[ -n "${PASS_ARG}" && ${METHODS_SET} -eq 1 &&
      "${METHODS_ARG}" != "${EXPECTED_BATCH_METHODS}" ]]; then
  echo "ERROR: batch-effect pass requires the fixed ordered method suite: ${EXPECTED_BATCH_METHODS}" >&2
  exit 1
fi
if [[ -n "${PASS_ARG}" && -n "${SELECTION_FILE_ARG}" ]]; then
  [[ -r "${SELECTION_FILE_ARG}" ]] || { echo "ERROR: selection file is unreadable." >&2; exit 1; }
  selection_invalid=0
  while IFS=$'\t' read -r selection_ds selection_view selection_label; do
    [[ -n "${selection_ds}" && -n "${selection_view}" && -n "${selection_label}" ]] || selection_invalid=1
    [[ "${selection_view}" == "batch_effect_${PASS_ARG}" ]] || selection_invalid=1
    [[ "${selection_label}" == "batch_effect_${PASS_ARG}" ]] || selection_invalid=1
  done < "${SELECTION_FILE_ARG}"
  [[ ${selection_invalid} -eq 0 ]] || {
    echo "ERROR: pass-mode selection rows must use batch_effect_${PASS_ARG}." >&2
    exit 1
  }
fi

# Exact validation is a preflight: reject malformed input before any run-root,
# pending-manifest, owner, or scheduler state can be created.
if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
  [[ -n "${SELECTION_FILE_ARG}" ]] || {
    echo "ERROR: --exact-batch-selection requires --selection-file." >&2
    exit 1
  }
  [[ "${PASS_ARG}" == "uncorrected" ]] || {
    echo "ERROR: --exact-batch-selection requires --pass uncorrected." >&2
    exit 1
  }
  [[ -r "${SELECTION_FILE_ARG}" ]] || {
    echo "ERROR: exact batch selection file is unreadable: ${SELECTION_FILE_ARG}" >&2
    exit 1
  }
  ecoda_validate_exact_batch_selection "${SELECTION_FILE_ARG}" 3 || exit 1
fi

if [[ -z "${PASS_ARG}" && -n "${SYNC_ONLY_RUN}" ]]; then
  sync_metadata="${HPC_SCRATCH_DIR}/_ecoda_runs/${SYNC_ONLY_RUN}/metadata"
  if [[ -r "${sync_metadata}" ]]; then
    sync_pass="$(sed -n 's/^PASS=//p' "${sync_metadata}" | head -1 || true)"
    case "${sync_pass}" in
      uncorrected|corrected) PASS_ARG="${sync_pass}" ;;
    esac
  fi
fi

# Pass-sensitive helper paths (notably watchdog status roots) are resolved
# only after --pass has been parsed and validated.
if [[ -n "${PASS_ARG}" ]]; then
  export ANALYSIS_PASS="${PASS_ARG}"
else
  unset ANALYSIS_PASS
fi
if [[ -n "${PASS_ARG}" ]]; then
  unset BENCHMARK_MANIFEST
fi
source "${SCRIPT_DIR}/benchmark_submit_common.sh"

stage5_finalize_owner_manifest() {
  local state="$1" reason="$2" owner_file="${ECODA_RUN_ROOT:-}/manifests/owners.tsv"
  local owner_key owner rc=0
  [[ -r "${owner_file}" ]] || return 1
  [[ -s "${owner_file}" ]] || return 0
  while IFS=$'\t' read -r owner_key owner; do
    [[ -n "${owner_key}" && -n "${owner}" ]] || { rc=1; continue; }
    if ! ecoda_owner_set_state "${owner}" "${state}" "${reason}"; then
      rc=1
    fi
  done < "${owner_file}"
  return "${rc}"
}

stage5_abort() {
  local reason="$1"
  local rc=0
  ecoda_owner_finalize_tracked FAIL "${reason}" || rc=1
  if [[ -n "${ECODA_RUN_ROOT:-}" && -r "${ECODA_RUN_ROOT}/manifests/owners.tsv" ]]; then
    stage5_finalize_owner_manifest FAIL "${reason}" || rc=1
  fi
  if [[ -n "${ECODA_RUN_ROOT:-}" ]]; then
    ecoda_set_run_state FAIL "${reason}" || rc=1
  fi
  echo "ERROR: ${reason}" >&2
  exit 1
}
stage5_record_scheduler() {
  local kind="$1" scheduler_id="$2"
  local tmp="${SCHEDULER_FILE}.record.$$" existing_kind existing_id
  [[ "${kind}" == "ARRAY" || "${kind}" == "WATCHDOG" ||
     "${kind}" == "STATUS" || "${kind}" == "AGGREGATE_GATE" ||
     "${kind}" == "PREFLIGHT" ]] || return 1
  [[ "${scheduler_id}" =~ ^[0-9]+$ ]] || return 1
  if [[ -s "${SCHEDULER_FILE}" ]]; then
    while IFS=$'\t' read -r existing_kind existing_id; do
      [[ -n "${existing_kind}" && "${existing_id}" =~ ^[0-9]+$ ]] || return 1
      [[ "${existing_id}" == "${scheduler_id}" ]] && return 0
    done < "${SCHEDULER_FILE}"
    cp "${SCHEDULER_FILE}" "${tmp}" || return 1
  else
    : > "${tmp}" || return 1
  fi
  printf '%s\t%s\n' "${kind}" "${scheduler_id}" >> "${tmp}" || {
    rm -f "${tmp}"
    return 1
  }
  mv -f "${tmp}" "${SCHEDULER_FILE}" || {
    rm -f "${tmp}"
    return 1
  }
}


method_spec() {
  local method="$1"
  METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
  METHOD_THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
  METHOD_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}" --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
  METHOD_WORKER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
  case "${method}" in
    mrvi|pilotgm|scpoli)
      METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_GPU}"
      METHOD_THROTTLE="${BENCHMARK_GPU_ARRAY_THROTTLE}"
      METHOD_FLAGS=(--gpus="${BENCHMARK_GPU_COUNT}" --constraint="${BENCHMARK_GPU_CONSTRAINT}" --cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}")
      METHOD_WORKER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
      ;;
    pilot|qot)
      METHOD_WORKER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
      ;;
    gloscope|mofa|pseudobulk|composition|scitd|prepare_pseudobulk)
      ;;
    trans|zeroimp)
      METHOD_WORKER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh"
      METHOD_FLAGS=(--cpus-per-task=4)
      ;;
    *) echo "ERROR: unsupported benchmark method/analysis '${method}'." >&2; return 1 ;;
  esac
  if [[ -n "${PARTITION_ARG}" ]]; then
    METHOD_PARTITION="${PARTITION_ARG}"
    filtered=(); flag=""
    for flag in "${METHOD_FLAGS[@]}"; do
      case "${flag}" in --constraint|--constraint=*) ;; *) filtered+=("${flag}");; esac
    done
    METHOD_FLAGS=("${filtered[@]}")
  fi
}

validate_input_row() {
  local ds="$1" view="$2" name path
  name="$(ecoda_view_output_name "${ds}" "${view}")"; [[ -n "${name}" ]] || return 1
  path="${HPC_SCRATCH_DIR}/${ds}/output/${name}"
  # Full persisted-content validation runs in the compute-node preflight
  # array.  The login submitter only checks presence here; source identity
  # creation below validates the sidecar, size, MD5, and Sample column.
  [[ -s "${path}" ]]
}

stage5_repair_missing_h5ad_sidecars() {
  local repair_manifest="${ECODA_RUN_ROOT}/manifests/h5ad_checksum_repairs.tsv"
  local repair_tmp="${repair_manifest}.build.$$"
  local ds view path sidecar digest size root repairs=0
  local scratch_path nas_path scratch_digest nas_digest source_key seen_sources=""
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0
  : > "${repair_tmp}" || return 1
  while IFS=$'\t' read -r ds view _scope; do
    source_key="${ds}/${view}"
    case " ${seen_sources} " in
      *" ${source_key} "*) continue ;;
    esac
    seen_sources="${seen_sources} ${source_key}"
    scratch_path="${HPC_SCRATCH_DIR}/${ds}/output/$(ecoda_view_output_name "${ds}" "${view}")"
    nas_path="${NAS_TARGET_DIR}/${ds}/output/$(ecoda_view_output_name "${ds}" "${view}")"
    [[ -s "${scratch_path}" && -s "${nas_path}" ]] || {
      rm -f "${repair_tmp}"
      return 1
    }
    scratch_digest="$(ecoda_md5_file "${scratch_path}")" || {
      rm -f "${repair_tmp}"
      return 1
    }
    nas_digest="$(ecoda_md5_file "${nas_path}")" || {
      rm -f "${repair_tmp}"
      return 1
    }
    [[ "${scratch_digest}" == "${nas_digest}" ]] || {
      rm -f "${repair_tmp}"
      return 1
    }
    for root in "${HPC_SCRATCH_DIR}" "${NAS_TARGET_DIR}"; do
      path="${root}/${ds}/output/$(ecoda_view_output_name "${ds}" "${view}")"
      sidecar="${path}.md5"
      if [[ -e "${sidecar}" || -L "${sidecar}" ]]; then
        ecoda_validate_checksum "${path}" || {
          rm -f "${repair_tmp}"
          return 1
        }
        continue
      fi
      # A missing sidecar is repaired only after the persisted H5AD content
      # contract passes. Existing invalid sidecars are never overwritten.
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
        --path "${path}" --view "${view}" \
        --method "Stage 5 source checksum repair" >/dev/null 2>&1 || {
        rm -f "${repair_tmp}"
        return 1
      }
      ecoda_write_checksum "${path}" || {
        rm -f "${repair_tmp}"
        return 1
      }
      ecoda_validate_checksum "${path}" || {
        rm -f "${repair_tmp}"
        return 1
      }
      digest="$(sed -n 's/^MD5=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
      size="$(sed -n 's/^SIZE=//p' "${sidecar}" | head -1 | tr -d '[:space:]')"
      printf '%s\t%s\t%s\t%s\t%s\n' "${root}" "${ds}" "${view}" \
        "${path}" "${digest}:${size}" >> "${repair_tmp}" || {
        rm -f "${repair_tmp}"
        return 1
      }
      repairs=$((repairs + 1))
    done
  done < "${MANIFEST}"
  if [[ ${repairs} -gt 0 ]]; then
    ecoda_atomic_install_manifest "${repair_tmp}" "${repair_manifest}" 5 || {
      rm -f "${repair_tmp}"
      return 1
    }
    rm -f "${repair_tmp}"
    ecoda_write_checksum "${repair_manifest}" || return 1
  else
    rm -f "${repair_tmp}"
  fi
  return 0
}

stage5_prepare_source_identity() {
  stage5_repair_missing_h5ad_sidecars || return 1
  local identity="${SOURCE_IDENTITY}" identity_sidecar="${SOURCE_IDENTITY}.md5"
  if [[ -e "${identity}" || -L "${identity}" ]]; then
    ecoda_validate_run_owned_path "${identity}" "${ECODA_RUN_ROOT}" || return 1
  fi
  if [[ -e "${identity_sidecar}" || -L "${identity_sidecar}" ]]; then
    ecoda_validate_run_owned_path "${identity_sidecar}" "${ECODA_RUN_ROOT}" || return 1
  fi
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0
  if [[ -s "${identity}" && -s "${identity_sidecar}" ]]; then
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/h5ad_source_identity.py" \
      --identity "${identity}" --selection "${MANIFEST}" \
      --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" ||
      return 1
  else
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/h5ad_source_identity.py" \
      --output "${identity}" --selection "${MANIFEST}" \
      --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" ||
      return 1
    ecoda_write_checksum "${identity}" || return 1
  fi
  ecoda_validate_checksum "${identity}" || return 1
}

stage5_compute_h5ad_preflight() {
  local preflight_manifest="${ECODA_RUN_ROOT}/manifests/h5ad_preflight.tsv"
  local preflight_tmp="${preflight_manifest}.build.$$"
  local expected_tmp="${preflight_manifest}.expected.$$"
  local status_dir="${ECODA_RUN_ROOT}/status/h5ad_preflight"
  local preflight_logs="${ECODA_RUN_ROOT}/logs"
  local ds view source_path safe status state preflight_row seen_rows
  local preflight_id preflight_rc count=0
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0

  if [[ -n "${SYNC_ONLY_RUN}" && ! -s "${preflight_manifest}" ]]; then
    # Runs created before the compute boundary have no preflight record.  Do
    # not resubmit during recovery; validate the immutable selected sources
    # locally and fail closed if any source is invalid.
    while IFS=$'\t' read -r ds view _scope; do
      source_path="${HPC_SCRATCH_DIR}/${ds}/output/$(ecoda_view_output_name "${ds}" "${view}")"
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
        --path "${source_path}" --view "${view}" --method "Stage 5 sync recovery" >/dev/null 2>&1 ||
        return 1
      ecoda_validate_checksum "${source_path}" || return 1
    done < "${MANIFEST}"
    return 0
  fi

  if [[ -n "${SYNC_ONLY_RUN}" ]]; then
    ecoda_validate_run_owned_path "${preflight_manifest}" "${ECODA_RUN_ROOT}" || return 1
    ecoda_validate_manifest "${preflight_manifest}" 3 || return 1
    ecoda_validate_checksum "${preflight_manifest}" || return 1
    : > "${expected_tmp}" || return 1
    seen_rows=""
    while IFS=$'\t' read -r ds view _scope; do
      preflight_row="${ds}/${view}"
      case " ${seen_rows} " in *" ${preflight_row} "*) continue ;; esac
      seen_rows="${seen_rows} ${preflight_row}"
      source_path="${HPC_SCRATCH_DIR}/${ds}/output/$(ecoda_view_output_name "${ds}" "${view}")"
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${source_path}" >> "${expected_tmp}" || return 1
    done < "${MANIFEST}"
    if ! cmp -s "${expected_tmp}" "${preflight_manifest}"; then
      rm -f "${expected_tmp}"
      return 1
    fi
    rm -f "${expected_tmp}"
  else
    : > "${preflight_tmp}" || return 1
    seen_rows=""
    while IFS=$'\t' read -r ds view _scope; do
      preflight_row="${ds}/${view}"
      case " ${seen_rows} " in *" ${preflight_row} "*) continue ;; esac
      seen_rows="${seen_rows} ${preflight_row}"
      source_path="${HPC_SCRATCH_DIR}/${ds}/output/$(ecoda_view_output_name "${ds}" "${view}")"
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${source_path}" >> "${preflight_tmp}" || return 1
    done < "${MANIFEST}"
    ecoda_atomic_install_manifest "${preflight_tmp}" "${preflight_manifest}" 3 || {
      rm -f "${preflight_tmp}"
      return 1
    }
    rm -f "${preflight_tmp}"
    ecoda_write_checksum "${preflight_manifest}" || return 1
    mkdir -p "${status_dir}" "${preflight_logs}" || return 1
    rm -f "${status_dir}"/*.status
    set +e
    preflight_id="$(
      ecoda_submit_h5ad_preflight "${preflight_manifest}" "${status_dir}" \
        "${ECODA_RUN_ROOT}" require "${PARTITION_ARG:-${SLURM_PARTITION_BENCHMARK_CPU}}" \
        "${MEMORY}" "${THROTTLE}" "${preflight_logs}" stage5 \
        "${SCRIPT_DIR}/../utils/bash/h5ad_preflight_worker.sh"
    )"
    preflight_rc=$?
    set -e
    if [[ "${preflight_id}" =~ ^[0-9]+$ ]]; then
      stage5_record_scheduler PREFLIGHT "${preflight_id}" || return 1
    fi
    [[ ${preflight_rc} -eq 0 && "${preflight_id}" =~ ^[0-9]+$ ]] || return 1
  fi

  while IFS=$'\t' read -r ds view _scope; do
    safe="$(_ecoda_safe_component "${ds}__${view}")"
    status="${status_dir}/${safe}.status"
    [[ -s "${status}" ]] || return 1
    state="$(sed -n 's/^STATE=//p' "${status}" | head -1)"
    [[ "${state}" == OK ]] || return 1
    count=$((count + 1))
  done < "${MANIFEST}"
  [[ ${count} -gt 0 ]] || return 1
}

# Resolve the immutable dataset/view selection before any scheduler submission.
RUN_ID="${ECODA_RUN_ID:-$(ecoda_new_run_id stage5)}"
if [[ -n "${SYNC_ONLY_RUN}" ]]; then
  ecoda_open_run "${SYNC_ONLY_RUN}" || exit 1
  RUN_ID="${SYNC_ONLY_RUN}"
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  ecoda_validate_run_owned_path "${MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage5_abort "Stage 5 selection manifest is not run-owned"
  ecoda_validate_manifest "${MANIFEST}" 3 ||
    stage5_abort "Stage 5 selection manifest is invalid"
  ecoda_validate_checksum "${MANIFEST}" ||
    stage5_abort "Stage 5 selection checksum is invalid"
  methods_meta="$(sed -n 's/^METHODS=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  [[ -n "${METHODS_ARG}" ]] || METHODS_ARG="${methods_meta}"
  pass_meta="$(sed -n 's/^PASS=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  [[ -n "${PASS_ARG}" ]] || PASS_ARG="${pass_meta}"
  [[ -n "${PASS_ARG}" ]] && unset BENCHMARK_MANIFEST
else
  ecoda_init_run stage5 "${RUN_ID}" >/dev/null
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  MANIFEST_TMP="${MANIFEST}.build.$$"
  : > "${MANIFEST_TMP}"
  if [[ -n "${SELECTION_FILE_ARG}" ]]; then
    [[ -r "${SELECTION_FILE_ARG}" ]] || stage5_abort "selection file is unreadable"
    while IFS=$'\t' read -r ds input_view row_label; do
      [[ -n "${ds}" && -n "${input_view}" && -n "${row_label}" ]] ||
        stage5_abort "selection rows require DATASET<TAB>VIEW<TAB>LABEL"
      if [[ -n "${PASS_ARG}" ]]; then
        [[ "${input_view}" == "batch_effect_${PASS_ARG}" ]] ||
          stage5_abort "pass-mode selection has the wrong view: ${input_view}"
        view="${input_view}"
      else
        view="${input_view}"
        case "${view}" in trans|zeroimp) view="benchmark_analysis" ;; esac
      fi
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${row_label}" >> "${MANIFEST_TMP}"
    done < "${SELECTION_FILE_ARG}"
  else
    DATASET_NAMES_TMP=()
    if [[ -n "${DATASETS_ARG}" ]]; then
      ecoda_split_csv "${DATASETS_ARG}" || stage5_abort "invalid benchmark dataset selection"
      DATASET_NAMES_TMP=("${ECODA_ARRAY[@]}")
      ecoda_assert_unique_items "${DATASET_NAMES_TMP[@]}" ||
        stage5_abort "duplicate benchmark dataset selection"
    elif [[ -n "${PASS_ARG}" ]]; then
      while IFS= read -r ds; do DATASET_NAMES_TMP+=("${ds}"); done < <(
        jq -r 'to_entries[] | select(.value.use_for_batch_effect == true) |
          .key | select(startswith("_") | not)' "${DATASETS_JSON_FILE}"
      )
    else
      while IFS= read -r ds; do DATASET_NAMES_TMP+=("${ds}"); done < <(
        jq -r 'to_entries[] | select(.value.use_for_benchmark == true) |
          select(.value.views.benchmark_analysis != null) |
          .key | select(startswith("_") | not)' "${DATASETS_JSON_FILE}"
      )
    fi
    [[ ${#DATASET_NAMES_TMP[@]} -gt 0 ]] || stage5_abort "no benchmark datasets selected"
    for ds in "${DATASET_NAMES_TMP[@]}"; do
      ecoda_dataset_exists "${ds}" || stage5_abort "unknown dataset ${ds}"
      if [[ -n "${PASS_ARG}" ]]; then
        view="batch_effect_${PASS_ARG}"
        ecoda_view_exists "${ds}" "${view}" ||
          stage5_abort "${ds}/${view} is not declared"
      else
        view="benchmark_analysis"
      fi
      [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" &&
        -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] ||
        stage5_abort "${ds}/${view} has no input/output"
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${view}" >> "${MANIFEST_TMP}"
    done
  fi
  ecoda_atomic_install_manifest "${MANIFEST_TMP}" "${MANIFEST}" 3 ||
    stage5_abort "failed to install Stage 5 selection atomically"
  rm -f "${MANIFEST_TMP}"
  ecoda_write_checksum "${MANIFEST}" || stage5_abort "failed to checksum Stage 5 selection"
fi
SOURCE_IDENTITY="${ECODA_RUN_ROOT}/manifests/source_identity.json"
SCHEDULER_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
if [[ -z "${SYNC_ONLY_RUN}" ]]; then
  ecoda_atomic_write "${SCHEDULER_FILE}" "" ||
    stage5_abort "failed to initialize Stage 5 scheduler manifest"
fi
if [[ -n "${PASS_ARG}" && -n "${METHODS_ARG}" ]]; then
  [[ "${METHODS_ARG}" == "${EXPECTED_BATCH_METHODS}" ]] ||
    stage5_abort "batch-effect pass requires the fixed ordered method suite"
fi

# Validate row syntax, owners, and source h5ads. No biological label is used as
# a processing covariate; LABEL is only a scheduler/output grouping token.
ecoda_validate_manifest "${MANIFEST}" 3 ||
  stage5_abort "invalid Stage 5 selection"
if [[ ${EXACT_BATCH_SELECTION} -eq 1 ]]; then
  ecoda_validate_exact_batch_selection "${MANIFEST}" 3 ||
    stage5_abort "invalid exact Stage 5 selection"
fi
DATASET_NAMES=(); SEEN_DS=""; SEEN_ROW=""
while IFS=$'\t' read -r ds view row_label; do
  ecoda_dataset_exists "${ds}" || stage5_abort "unknown dataset ${ds}"
  ecoda_view_exists "${ds}" "${view}" ||
    stage5_abort "${ds}/${view} is not declared"
  if [[ -n "${PASS_ARG}" ]]; then
    [[ "${view}" == "batch_effect_${PASS_ARG}" ]] ||
      stage5_abort "pass-mode selection has the wrong view for ${ds}"
  fi
  if [[ "${PASS_ARG}" == corrected ]] && [[ -z "$(jq -r --arg ds "${ds}" '.[$ds].columns.batch // empty' "${DATASETS_JSON_FILE}")" ]]; then
    stage5_abort "corrected batch-effect view requires a confirmed columns.batch"
  fi
  row="${ds}/${view}/${row_label}"
  case " ${SEEN_ROW} " in *" ${row} "*) stage5_abort "duplicate Stage 5 selection ${row}" ;; esac
  SEEN_ROW="${SEEN_ROW} ${row}"
  case " ${SEEN_DS} " in *" ${ds} "*) ;; *) DATASET_NAMES+=("${ds}"); SEEN_DS="${SEEN_DS} ${ds}" ;; esac
  if [[ "${BENCHMARK_MATRIX_TEST:-0}" != 1 ]]; then
    validate_input_row "${ds}" "${view}" || stage5_abort "invalid Stage 5 h5ad ${row}"
  fi
done < "${MANIFEST}"
stage5_prepare_source_identity ||
  stage5_abort "failed to build or verify Stage 5 source identity"
stage5_compute_h5ad_preflight ||
  stage5_abort "Stage 5 compute-node H5AD preflight failed"
METHODS=(gloscope mofa pseudobulk composition scitd mrvi scpoli pilot qot pilotgm)
ANALYSES=(_ecoda_none_)
ANALYSES_SELECTED=0
EXACT_SELECTION=0
[[ ${EXACT_BATCH_SELECTION} -eq 1 ]] && EXACT_SELECTION=1
BATCH_EFFECT_METHODS=(prepare_pseudobulk pseudobulk gloscope composition mrvi pilot pilotgm qot)
if [[ -n "${PASS_ARG}" ]]; then METHODS=("${BATCH_EFFECT_METHODS[@]}"); fi
if [[ -n "${METHODS_ARG}" ]]; then
  ecoda_split_csv "${METHODS_ARG}" || stage5_abort "invalid benchmark method selection"
  METHODS=("${ECODA_ARRAY[@]}")
fi
if [[ -n "${PASS_ARG}" ]]; then
  methods_csv="$(IFS=,; echo "${METHODS[*]}")"
  [[ "${methods_csv}" == "${EXPECTED_BATCH_METHODS}" ]] ||
    stage5_abort "batch-effect pass requires the fixed ordered method suite"
fi
if [[ -n "${ANALYSES_ARG}" ]]; then
  ecoda_split_csv "${ANALYSES_ARG}" || stage5_abort "invalid benchmark analysis selection"
  ANALYSES=("${ECODA_ARRAY[@]}")
  ANALYSES_SELECTED=1
  [[ -n "${METHODS_ARG}" || -n "${PASS_ARG}" ]] || METHODS=(_ecoda_none_)
fi
if [[ -n "${SELECTION_FILE_ARG}" && -z "${METHODS_ARG}" && -z "${ANALYSES_ARG}" && -z "${PASS_ARG}" ]]; then
  EXACT_SELECTION=1
  METHODS=(_ecoda_none_); ANALYSES=(_ecoda_none_); ANALYSES_SELECTED=0
  while IFS=$'\t' read -r ds view row_label; do
    case "${row_label}" in
      trans|zeroimp)
        [[ "${ANALYSES[*]}" == *"${row_label}"* ]] || ANALYSES+=("${row_label}")
        ANALYSES_SELECTED=1
        ;;
      *)
        [[ "${METHODS[*]}" == *"${row_label}"* ]] || METHODS+=("${row_label}")
        ;;
    esac
  done < "${MANIFEST}"
  if [[ ${#METHODS[@]} -gt 1 ]]; then METHODS=("${METHODS[@]:1}"); fi
  if [[ ${#ANALYSES[@]} -gt 1 ]]; then ANALYSES=("${ANALYSES[@]:1}"); fi
fi

if [[ -n "${PASS_ARG}" ]]; then
  [[ ${ANALYSES_SELECTED} -eq 0 ]] ||
    stage5_abort "batch-effect pass does not accept ordinary analyses"
  for method in "${METHODS[@]}"; do
    [[ "${method}" == _ecoda_none_ ]] && continue
    case " ${BATCH_EFFECT_METHODS[*]} " in
      *" ${method} "*) ;;
      *) stage5_abort "unsupported batch-effect method ${method}" ;;
    esac
  done
fi
if [[ "${METHODS[*]}" == *"_ecoda_none_"* && ${ANALYSES_SELECTED} -eq 0 ]]; then
  stage5_abort "no benchmark methods or analyses selected"
fi
if [[ ${ANALYSES_SELECTED} -eq 1 ]]; then
  ecoda_assert_unique_items "${METHODS[@]}" "${ANALYSES[@]}" ||
    stage5_abort "duplicate benchmark methods or analyses"
else
  ecoda_assert_unique_items "${METHODS[@]}" || stage5_abort "duplicate benchmark methods"
fi
for method in "${METHODS[@]}"; do
  [[ "${method}" == _ecoda_none_ ]] && continue
  method_spec "${method}" || stage5_abort "unsupported method"
done
for method in "${ANALYSES[@]}"; do
  [[ "${method}" == _ecoda_none_ ]] && continue
  method_spec "${method}" || stage5_abort "unsupported analysis"
done
NEEDS_PREP=0
for method in "${METHODS[@]}"; do case "${method}" in mofa|pseudobulk|composition) NEEDS_PREP=1;; esac; done
if [[ ${NEEDS_PREP} -eq 1 ]]; then case " ${METHODS[*]} " in *" prepare_pseudobulk "*) ;; *) METHODS=(prepare_pseudobulk "${METHODS[@]}");; esac; fi
export ECODA_SELECTION_MANIFEST="${MANIFEST}"
export ECODA_EXACT_SELECTION="${EXACT_SELECTION}"

if [[ -n "${PASS_ARG}" ]]; then
  ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/${PASS_ARG}"
else
  ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/benchmark"
fi

benchmark_artifacts_for() {
  local ds="$1" view="$2" label="$3"
  local stem suffix n
  ARTIFACT_PATHS=()
  if [[ "${label}" == prepare_pseudobulk ]]; then
    if [[ -n "${PASS_ARG}" ]]; then
      ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/pseudobulks/${ds}_batch_effect_${PASS_ARG}_pseudobulk_hvg2000.rds")
    else
      for stem in schvg2000 hvg2000 hvg500 hvg2000_bl hvg1000 hvg3000; do
        ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/pseudobulks/${ds}_pseudobulk_${stem}.rds")
      done
    fi
    return 0
  fi
  case "${label}" in
    mrvi)
      if [[ -n "${PASS_ARG}" ]]; then
        ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_batch_effect_${PASS_ARG}_hvg2000_highres_mrvi_dists.feather")
      else
        for n in 1000 2000 3000; do
          ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_hvg${n}_mrvi_dists.feather")
        done
      fi
      ;;
    scpoli)
      for n in 2000; do ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_hvg${n}_lowres_scpoli_dims15_embs.feather"); done
      for n in 1000 3000; do ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_hvg${n}_highres_scpoli_dims15_embs.feather"); done
      for n in 2 3 5 10 15; do ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_hvg2000_highres_scpoli_dims${n}_embs.feather"); done
      ;;
    pilot|qot|pilotgm)
      suffix="${label}"
      if [[ -n "${PASS_ARG}" ]]; then
        ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_batch_effect_${PASS_ARG}_hvg2000_highres_${suffix}_dists.feather")
      else
        for n in 2000; do ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_hvg${n}_lowres_${suffix}_dists.feather"); done
        for n in 1000 2000 3000; do ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/embeddings/${ds}_hvg${n}_highres_${suffix}_dists.feather"); done
      fi
      ;;
    trans|zeroimp)
      ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/results/${ds}_${label}.rds")
      ;;
    gloscope|mofa|pseudobulk|scitd)
      stem="${ds}"
      [[ -n "${PASS_ARG}" ]] && stem="${ds}_batch_effect_${PASS_ARG}"
      ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/results/${stem}_${label}.rds")
      ;;
    composition)
      stem="${ds}"
      [[ -n "${PASS_ARG}" ]] && stem="${ds}_batch_effect_${PASS_ARG}"
      ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/results/${stem}_composition.rds")
      ARTIFACT_PATHS+=("${ANALYSIS_ROOT}/results/${stem}_metadata.rds")
      ;;
    *)
      echo "ERROR: unsupported benchmark artifact label '${label}'." >&2
      return 1
      ;;
  esac
}

RDS_PREFLIGHT_DONE=""
RDS_PREFLIGHT_FAILED=""

benchmark_rds_group_valid() {
  local ds="$1" view="$2" key="${PASS_ARG:-ordinary}/${ds}/${view}"
  local safe list list_tmp label path metadata group_rc
  local selected_labels="" seen_label=""
  local has_rds=0
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 2
  case " ${RDS_PREFLIGHT_DONE} " in *" ${key} "*) return 0 ;; esac
  case " ${RDS_PREFLIGHT_FAILED} " in *" ${key} "*) return 1 ;; esac
  safe="$(_ecoda_safe_component "${key}")"
  list="${ECODA_RUN_ROOT}/manifests/rds_preflight_${safe}.tsv"
  list_tmp="${list}.build.$$"
  : > "${list_tmp}" || return 1
  for label in "${METHODS[@]}" "${ANALYSES[@]}"; do
    [[ "${label}" == _ecoda_none_ ]] && continue
    case " ${seen_label} " in *" ${label} "*) continue ;; esac
    seen_label="${seen_label} ${label}"
    case "${label}" in
      gloscope|mofa|pseudobulk|composition|scitd|prepare_pseudobulk|trans|zeroimp) ;;
      *) continue ;;
    esac
    benchmark_artifacts_for "${ds}" "${view}" "${label}" || {
      rm -f "${list_tmp}"
      return 1
    }
    for path in "${ARTIFACT_PATHS[@]}"; do
      case "${path}" in
        *.rds)
          [[ -s "${path}" && -s "${path}.md5" ]] || {
            rm -f "${list_tmp}"
            return 2
          }
          metadata=0
          [[ "${path}" == *_metadata.rds ]] && metadata=1
          printf '%s\t%s\t%s\t%s\t%s\n' \
            "${path}" "${label}" "${ds}" "${view}" "${metadata}" >> "${list_tmp}" || {
            rm -f "${list_tmp}"
            return 1
          }
          has_rds=1
          ;;
      esac
    done
  done
  [[ ${has_rds} -eq 1 ]] || { rm -f "${list_tmp}"; return 2; }
  ecoda_atomic_install_manifest "${list_tmp}" "${list}" 5 || {
    rm -f "${list_tmp}"
    return 1
  }
  rm -f "${list_tmp}"
  ecoda_write_checksum "${list}" || return 1
  rds_args=(--artifact-list "${list}" --config "${DATASETS_JSON_FILE}" \
    --input-root "${HPC_SCRATCH_DIR}" --source-identity "${SOURCE_IDENTITY}" \
    --source-identity-verified)
  [[ -n "${PASS_ARG}" ]] && rds_args+=(--batch-pass "${PASS_ARG}")
  set +e
  ${PIXI_RSCRIPT} "${SCRIPT_DIR}/validate_benchmark_rds_contract.R" \
    "${rds_args[@]}" >/dev/null 2>&1
  group_rc=$?
  set -e
  if [[ ${group_rc} -ne 0 ]]; then
    RDS_PREFLIGHT_FAILED="${RDS_PREFLIGHT_FAILED} ${key}"
    return 1
  fi
  RDS_PREFLIGHT_DONE="${RDS_PREFLIGHT_DONE} ${key}"
  rm -f "${list}" "${list}.md5"
  return 0
}

benchmark_selected_artifacts_valid() {
  local ds="$1" view="$2" label="$3" path artifact_check
  local has_feather=0 rds_grouped=0 group_rc
  case "${label}" in
    gloscope|mofa|pseudobulk|composition|scitd|prepare_pseudobulk|trans|zeroimp)
      if benchmark_rds_group_valid "${ds}" "${view}"; then
        rds_grouped=1
      else
        group_rc=$?
        [[ ${group_rc} -eq 2 ]] || return 1
      fi
      ;;
  esac
  benchmark_artifacts_for "${ds}" "${view}" "${label}" || return 1
  [[ ${#ARTIFACT_PATHS[@]} -gt 0 ]] || return 1
  for path in "${ARTIFACT_PATHS[@]}"; do
    ecoda_validate_checksum "${path}" || return 1
    case "${path}" in
      *.feather)
        has_feather=1
        ;;
      *.rds)
        [[ ${rds_grouped} -eq 1 ]] && continue
        rds_args=(--artifact "${path}" --method "${label}" --dataset "${ds}" --view "${view}" \
          --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}")
        [[ -n "${PASS_ARG}" ]] && rds_args+=(--batch-pass "${PASS_ARG}")
        [[ "${path}" == *_metadata.rds ]] && rds_args+=(--metadata)
        [[ -s "${SOURCE_IDENTITY}" ]] && rds_args+=(--source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
        ${PIXI_RSCRIPT} "${SCRIPT_DIR}/validate_benchmark_rds_contract.R" \
          "${rds_args[@]}" >/dev/null 2>&1 || return 1
        ;;
      *) return 1 ;;
    esac
  done
  if [[ ${has_feather} -eq 1 ]]; then
    [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0
    if [[ -n "${PASS_ARG}" ]]; then
      artifact_check="${ECODA_RUN_ROOT}/manifests/artifact_check_${ds}_${label}.tsv"
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${label}" > "${artifact_check}.build.$$"
      mv -f "${artifact_check}.build.$$" "${artifact_check}"
      ecoda_write_checksum "${artifact_check}" || {
        rm -f "${artifact_check}"
        return 1
      }
      "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" \
        --root "${ANALYSIS_ROOT}" --selection "${artifact_check}" \
        --labels "${label}" --batch --batch-pass "${PASS_ARG}" \
        --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" \
        --source-identity "${SOURCE_IDENTITY}" --source-identity-verified >/dev/null 2>&1 || {
        rm -f "${artifact_check}" "${artifact_check}.md5"
        return 1
      }
      rm -f "${artifact_check}" "${artifact_check}.md5"
    else
      artifact_check="${ECODA_RUN_ROOT}/manifests/artifact_check_${ds}_${label}.tsv"
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${label}" > "${artifact_check}.build.$$"
      mv -f "${artifact_check}.build.$$" "${artifact_check}"
      ecoda_write_checksum "${artifact_check}" || {
        rm -f "${artifact_check}"
        return 1
      }
      "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" \
        --root "${ANALYSIS_ROOT}" --selection "${artifact_check}" \
        --labels "${label}" --input-root "${HPC_SCRATCH_DIR}" \
        --config "${DATASETS_JSON_FILE}" --source-identity "${SOURCE_IDENTITY}" \
        --source-identity-verified >/dev/null 2>&1 || {
        rm -f "${artifact_check}" "${artifact_check}.md5"
        return 1
      }
      rm -f "${artifact_check}" "${artifact_check}.md5"
    fi
  fi
}

OWNERS_FILE="${ECODA_RUN_ROOT}/manifests/owners.tsv"
if [[ -z "${SYNC_ONLY_RUN}" ]]; then
  OWNERS_TMP="${OWNERS_FILE}.build.$$"
  PENDING_SELECTION="${ECODA_RUN_ROOT}/manifests/pending_selection.tsv"
  PENDING_SELECTION_TMP="${PENDING_SELECTION}.build.$$"
  : > "${OWNERS_TMP}"
  : > "${PENDING_SELECTION_TMP}"
  OWNER_SEEN=""
  ecoda_owner_clear_tracked
  for method in "${METHODS[@]}" "${ANALYSES[@]}"; do
    [[ "${method}" == _ecoda_none_ ]] && continue
    while IFS=$'\t' read -r ds view row_label; do
      if [[ ${EXACT_SELECTION} -eq 1 && -z "${PASS_ARG}" ]]; then
        case "${method}:${row_label}" in
          prepare_pseudobulk:mofa|prepare_pseudobulk:pseudobulk|prepare_pseudobulk:composition) ;;
          prepare_pseudobulk:*) continue ;;
          *:"${method}") ;;
          *) continue ;;
        esac
      fi
      owner_key="${PASS_ARG:-ordinary}/${ds}/${view}/${method}"
      case " ${OWNER_SEEN} " in
        *" ${owner_key} "*) continue ;;
      esac
      OWNER_SEEN="${OWNER_SEEN} ${owner_key}"
      if benchmark_selected_artifacts_valid "${ds}" "${view}" "${method}" &&
         [[ ${FORCE_ARG} -eq 0 ]]; then
        echo "Skipping validated Stage 5 artifact ${ds}/${view}/${method}."
        continue
      fi
      benchmark_artifacts_for "${ds}" "${view}" "${method}" ||
        stage5_abort "cannot resolve output contract for ${ds}/${view}/${method}"
      set +e
      owner_dir="$(ecoda_owner_acquire stage5 "${owner_key}" "${RUN_ID}" "${FORCE_ARG}" 0)"
      owner_rc=$?
      set -e
      [[ ${owner_rc} -eq 0 ]] || stage5_abort "ownership conflict for ${owner_key}"
      ecoda_owner_track "${owner_dir}" ||
        stage5_abort "failed to track owner for ${owner_key}"
      for path in "${ARTIFACT_PATHS[@]}"; do
        ecoda_invalidate_artifact "${path}" ||
          stage5_abort "failed to invalidate Stage 5 artifact ${path}"
      done
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${method}" >> "${PENDING_SELECTION_TMP}"
      printf '%s\t%s\n' "${owner_key}" "${owner_dir}" >> "${OWNERS_TMP}"
    done < "${MANIFEST}"
  done
  if [[ -s "${PENDING_SELECTION_TMP}" ]]; then
    ecoda_atomic_install_manifest "${PENDING_SELECTION_TMP}" "${PENDING_SELECTION}" 3 ||
      stage5_abort "failed to install Stage 5 pending manifest atomically"
    ecoda_atomic_install_manifest "${OWNERS_TMP}" "${OWNERS_FILE}" 2 ||
      stage5_abort "failed to install Stage 5 owner manifest atomically"
  else
    ecoda_atomic_write "${PENDING_SELECTION}" "" ||
      stage5_abort "failed to create empty Stage 5 pending manifest"
    ecoda_atomic_write "${OWNERS_FILE}" "" ||
      stage5_abort "failed to create empty Stage 5 owner manifest"
  fi
  rm -f "${PENDING_SELECTION_TMP}" "${OWNERS_TMP}"
else
  ecoda_validate_run_owned_path "${OWNERS_FILE}" "${ECODA_RUN_ROOT}" ||
    stage5_abort "Stage 5 owner manifest is missing or not run-owned"
  [[ ! -s "${OWNERS_FILE}" ]] || ecoda_validate_manifest "${OWNERS_FILE}" 2 ||
    stage5_abort "Stage 5 owner manifest is invalid"
fi

if [[ -n "${SYNC_ONLY_RUN}" ]]; then
  [[ -s "${ECODA_RUN_ROOT}/status/aggregate" ]] &&
    grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/aggregate" ||
    stage5_abort "Stage 5 aggregate status is not OK"
else
  methods_csv=""
  for method in "${METHODS[@]}"; do
    [[ "${method}" == _ecoda_none_ ]] || {
      [[ -n "${methods_csv}" ]] && methods_csv="${methods_csv},"
      methods_csv="${methods_csv}${method}"
    }
  done
  analyses_csv=""
  if [[ ${ANALYSES_SELECTED} -eq 1 ]]; then analyses_csv="$(IFS=,; echo "${ANALYSES[*]}")"; fi
  RUN_METADATA="STAGE=stage5\nRUN_ID=${RUN_ID}\nSTATE=ACTIVE\nMETHODS=${methods_csv}\nANALYSES=${analyses_csv}\nPASS=${PASS_ARG}\nEXACT_SELECTION=${EXACT_SELECTION}\nSOURCE_IDENTITY=${SOURCE_IDENTITY}\nH5AD_PREFLIGHT=${ECODA_RUN_ROOT}/manifests/h5ad_preflight.tsv\nROOT=${HPC_SCRATCH_DIR}/$([[ -n "${PASS_ARG}" ]] && printf 'batch_effect/%s' "${PASS_ARG}" || printf 'benchmark')\n"
  ecoda_atomic_write "${ECODA_RUN_ROOT}/metadata" "${RUN_METADATA}" ||
    stage5_abort "failed to write Stage 5 run metadata"
fi

if [[ -n "${PASS_ARG}" ]]; then
  ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/${PASS_ARG}"
  ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/${PASS_ARG}"
  ANALYSIS_LOG_PREFIX="execution_times_batch_effect_${PASS_ARG}_"
  export ANALYSIS_PASS="${PASS_ARG}" ANALYSIS_HIGH_RES_ONLY=1
else
  ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/benchmark"
  ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/benchmark"
  ANALYSIS_LOG_PREFIX="execution_times_"
  unset ANALYSIS_PASS ANALYSIS_HIGH_RES_ONLY
fi
RUN_LOG_DIR="${ECODA_RUN_ROOT}/logs"
export ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_LOG_PREFIX FORCE_BENCHMARK="${FORCE_ARG}"
export ECODA_RUN_ROOT ECODA_RUN_ID EXECUTION_LOG_DIR="${RUN_LOG_DIR}" LOGS_DIR="${RUN_LOG_DIR}"
unset BENCHMARK_MANIFEST
mkdir -p "${ANALYSIS_ROOT}/embeddings" "${ANALYSIS_ROOT}/results" \
  "${ANALYSIS_ROOT}/pseudobulks" "${ANALYSIS_ROOT}/gloscope_dists" "${RUN_LOG_DIR}"

if [[ -z "${SYNC_ONLY_RUN}" ]]; then
  WATCHDOG_LABELS=()
  WATCHDOG_IDS=()
  ARRAY_IDS=()
  LABELS=()
  GROUP_KEYS=()
  GROUP_VIEWS=()
  GROUP_LABELS=()
  GROUP_MANIFESTS=()
  GROUP_TMP_FILES=()
  PREP_VIEWS=()
  PREP_WATCHDOGS=()

  group_add_row() {
    local ds="$1" view="$2" label="$3" key="${2}|${3}"
    local idx=0 safe
    while [[ ${idx} -lt ${#GROUP_KEYS[@]} && "${GROUP_KEYS[${idx}]}" != "${key}" ]]; do
      idx=$((idx + 1))
    done
    if [[ ${idx} -eq ${#GROUP_KEYS[@]} ]]; then
      safe="$(printf '%s' "${key}" | tr '/:,\t |' '______')"
      GROUP_KEYS+=("${key}")
      GROUP_VIEWS+=("${view}")
      GROUP_LABELS+=("${label}")
      GROUP_MANIFESTS+=("${ECODA_RUN_ROOT}/manifests/matrix_${safe}.tsv")
      GROUP_TMP_FILES+=("${ECODA_RUN_ROOT}/manifests/matrix_${safe}.build.$$")
      : > "${GROUP_TMP_FILES[${idx}]}"
    fi
    printf '%s\t%s\t%s\n' "${ds}" "${view}" "${label}" >> "${GROUP_TMP_FILES[${idx}]}"
  }

  submit_matrix() {
    local group_label="$1" manifest="$2" dependency="$3" method="$4" view="$5"
    method_spec "${method}" || return 1
    local worker="${METHOD_WORKER}" worker_partition="${METHOD_PARTITION}" \
      watchdog_partition="${SLURM_PARTITION_BENCHMARK_CPU}" throttle="${METHOD_THROTTLE}"
    local safe="$(printf '%s' "${group_label}" | tr '/:,\t ' '_____')"
    local array_msg array_id array_rc wd_msg wd_id wd_rc
    local worker_env="METHOD=${method},ANALYSIS=${method},ANALYSIS_MANIFEST=${manifest},ANALYSIS_VIEW=${view},ANALYSIS_ROOT=${ANALYSIS_ROOT},EXECUTION_LOG_DIR=${RUN_LOG_DIR},ECODA_RUN_ROOT=${ECODA_RUN_ROOT},ECODA_RUN_ID=${RUN_ID},FORCE_BENCHMARK=${FORCE_ARG},JOB_LOG_PREFIX=${RUN_LOG_DIR}/5_matrix_${safe}"
    if [[ -n "${PASS_ARG}" ]]; then
      worker_env="${worker_env},ANALYSIS_PASS=${PASS_ARG}"
    else
      worker_env="${worker_env},BENCHMARK_MANIFEST=${manifest}"
    fi
    array_args=(--parsable --array="1-$(wc -l < "${manifest}" | tr -d '[:space:]')%${throttle}" --partition="${worker_partition}" "${METHOD_FLAGS[@]}" --mem="${MEMORY}" \
      --output="${RUN_LOG_DIR}/5_matrix_${safe}_%A_%a.log" --error="${RUN_LOG_DIR}/5_matrix_${safe}_%A_%a.err" --mail-user="${USER_EMAIL}")
    [[ -n "${dependency}" ]] && array_args+=(--dependency="afterok:${dependency}")
    array_args+=(--export="ALL,${worker_env}" "${worker}")
    set +e
    array_msg="$(sbatch "${array_args[@]}")"
    array_rc=$?
    set -e
    array_id="${array_msg%%;*}"
    if [[ "${array_id}" =~ ^[0-9]+$ ]]; then
      stage5_record_scheduler ARRAY "${array_id}" || return 1
    else
      return 1
    fi
    [[ ${array_rc} -eq 0 ]] || return 1
    set +e
    wd_msg="$(sbatch --parsable --partition="${watchdog_partition}" --ntasks=1 --cpus-per-task=1 --mem=2G --time="${WATCHDOG_TIME_LIMIT}" \
      --output="${RUN_LOG_DIR}/5_matrix_watchdog_${safe}_%A.log" --error="${RUN_LOG_DIR}/5_matrix_watchdog_${safe}_%A.err" --mail-user="${USER_EMAIL}" \
      --export="ALL,${worker_env},MATRIX_WATCHDOG_ROOT=${ECODA_RUN_ROOT}" \
      "${SCRIPT_DIR}/matrix_watchdog.sh" "${ECODA_RUN_ROOT}" "${group_label}" "${manifest}" "${array_id}" "${MEMORY}" "${MAX_MEMORY}" "${worker_partition}" "${throttle}" "${worker}" "${METHOD_FLAGS[@]}")"
    wd_rc=$?
    set -e
    wd_id="${wd_msg%%;*}"
    if [[ "${wd_id}" =~ ^[0-9]+$ ]]; then
      stage5_record_scheduler WATCHDOG "${wd_id}" || return 1
    else
      return 1
    fi
    [[ ${wd_rc} -eq 0 ]] || return 1
    LAST_ARRAY_ID="${array_id}"
    LAST_WATCHDOG_ID="${wd_id}"
    ARRAY_IDS+=("${array_id}")
    WATCHDOG_IDS+=("${wd_id}")
    WATCHDOG_LABELS+=("${group_label}")
    if [[ -n "${PASS_ARG}" ]]; then
      echo "BATCH_EFFECT_ARRAY_JOB_ID=${group_label}:${array_id}"
      echo "BATCH_EFFECT_WATCHDOG_JOB_ID=${group_label}:${wd_id}"
    else
      echo "BENCHMARK_ARRAY_JOB_ID=${group_label}:${array_id}"
      echo "BENCHMARK_WATCHDOG_JOB_ID=${group_label}:${wd_id}"
    fi
  }

  if [[ -s "${PENDING_SELECTION}" ]]; then
    while IFS=$'\t' read -r ds view method; do
      [[ -n "${ds}" && -n "${view}" && -n "${method}" ]] || continue
      group_add_row "${ds}" "${view}" "${method}"
    done < "${PENDING_SELECTION}"
  fi
  for idx in "${!GROUP_KEYS[@]}"; do
    ecoda_atomic_install_manifest "${GROUP_TMP_FILES[${idx}]}" "${GROUP_MANIFESTS[${idx}]}" 3 ||
      stage5_abort "failed to install Stage 5 matrix manifest atomically"
    rm -f "${GROUP_TMP_FILES[${idx}]}"
  done

  for idx in "${!GROUP_KEYS[@]}"; do
    view="${GROUP_VIEWS[${idx}]}"
    method="${GROUP_LABELS[${idx}]}"
    if [[ "${method}" == prepare_pseudobulk ]]; then
      submit_matrix "${view}__${method}" "${GROUP_MANIFESTS[${idx}]}" "" "${method}" "${view}" ||
        stage5_abort "submission failed for ${method}/${view}"
      PREP_VIEWS+=("${view}")
      PREP_WATCHDOGS+=("${LAST_WATCHDOG_ID}")
    fi
  done
  for idx in "${!GROUP_KEYS[@]}"; do
    view="${GROUP_VIEWS[${idx}]}"
    method="${GROUP_LABELS[${idx}]}"
    [[ "${method}" == prepare_pseudobulk ]] && continue
    dependency=""
    case "${method}" in
      mofa|pseudobulk|composition)
        for prep_idx in "${!PREP_VIEWS[@]}"; do
          [[ "${PREP_VIEWS[${prep_idx}]}" == "${view}" ]] && dependency="${PREP_WATCHDOGS[${prep_idx}]}"
        done
        ;;
    esac
    submit_matrix "${view}__${method}" "${GROUP_MANIFESTS[${idx}]}" "${dependency}" "${method}" "${view}" ||
      stage5_abort "submission failed for ${method}/${view}"
  done
  if [[ ${#WATCHDOG_IDS[@]} -gt 0 ]]; then
    watchdog_ids_colon="$(IFS=:; echo "${WATCHDOG_IDS[*]}")"
    watchdog_labels_csv="$(IFS=,; echo "${WATCHDOG_LABELS[*]}")"
    set +e
    gate_msg="$(sbatch --parsable --wait --dependency="afterany:${watchdog_ids_colon}" --partition="${SLURM_PARTITION_BENCHMARK_CPU}" \
      --ntasks=1 --cpus-per-task=1 --mem=2G --time="${WATCHDOG_TIME_LIMIT}" --output="${RUN_LOG_DIR}/5_matrix_gate_%j.log" \
      --error="${RUN_LOG_DIR}/5_matrix_gate_%j.err" --mail-user="${USER_EMAIL}" \
      "${SCRIPT_DIR}/matrix_gate.sh" "${ECODA_RUN_ROOT}" "${watchdog_labels_csv}" "${SCHEDULER_FILE}")"
    gate_rc=$?
    set -e
    GATE_ID="${gate_msg%%;*}"
    if [[ "${GATE_ID}" =~ ^[0-9]+$ ]]; then
      stage5_record_scheduler AGGREGATE_GATE "${GATE_ID}" ||
        stage5_abort "failed to persist Stage 5 aggregate gate ID"
    fi
    [[ ${gate_rc} -eq 0 ]] ||
      stage5_abort "benchmark aggregate gate job failed"
    [[ "${GATE_ID}" =~ ^[0-9]+$ ]] ||
      stage5_abort "invalid benchmark aggregate gate id"
    if [[ -n "${PASS_ARG}" ]]; then
      echo "BATCH_EFFECT_AGGREGATE_GATE_JOB_ID=${GATE_ID}"
    else
      echo "BENCHMARK_AGGREGATE_GATE_JOB_ID=${GATE_ID}"
    fi
    if [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]]; then
      mkdir -p "${ECODA_RUN_ROOT}/status/watchdogs"
      for idx in "${!WATCHDOG_LABELS[@]}"; do
        safe="$(printf '%s' "${WATCHDOG_LABELS[${idx}]}" | tr '/:,\t |' '______')"
        ecoda_atomic_write "${ECODA_RUN_ROOT}/status/watchdogs/${safe}.status" \
          "STATE=OK\nLABEL=${WATCHDOG_LABELS[${idx}]}\nARRAY_JOB_ID=${ARRAY_IDS[${idx}]}\nWATCHDOG_JOB_ID=${WATCHDOG_IDS[${idx}]}\nSCHEDULER_ID=${ARRAY_IDS[${idx}]}\nSCHEDULER_ID=${WATCHDOG_IDS[${idx}]}\n"
      done
      ecoda_atomic_write "${ECODA_RUN_ROOT}/status/aggregate" \
        "STATE=OK\nWATCHDOG_LABELS=${watchdog_labels_csv}\nSCHEDULER_ID=${GATE_ID}\n"
    fi
  else
    ecoda_atomic_write "${ECODA_RUN_ROOT}/status/aggregate" "STATE=OK\nWATCHDOG_LABELS=\n" ||
      stage5_abort "failed to write empty Stage 5 aggregate status"
  fi
  if [[ ! -s "${ECODA_RUN_ROOT}/status/aggregate" ]] ||
     ! grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/aggregate"; then
    stage5_abort "benchmark aggregate gate did not report OK"
  fi
  for idx in "${!WATCHDOG_LABELS[@]}"; do
    safe="$(printf '%s' "${WATCHDOG_LABELS[${idx}]}" | tr '/:,\t |' '______')"
    status="${ECODA_RUN_ROOT}/status/watchdogs/${safe}.status"
    [[ -s "${status}" ]] ||
      stage5_abort "missing Stage 5 watchdog status ${WATCHDOG_LABELS[${idx}]}"
    while IFS= read -r status_line; do
      case "${status_line}" in
        SCHEDULER_ID=*|ARRAY_JOB_ID=*|WATCHDOG_JOB_ID=*)
          if [[ "${status_line}" == ARRAY_JOB_ID=* && -n "${PASS_ARG}" &&
                "${status_line#*=}" != "${ARRAY_IDS[${idx}]}" ]]; then
            printf 'BATCH_EFFECT_RETRY_ARRAY_JOB_ID=%s:%s\n' \
              "${WATCHDOG_LABELS[${idx}]}" "${status_line#*=}"
          fi
          if [[ -n "${PASS_ARG}" ]]; then
            printf 'BATCH_EFFECT_SCHEDULER_ID=%s:%s\n' "${WATCHDOG_LABELS[${idx}]}" "${status_line#*=}"
          else
            printf 'BENCHMARK_SCHEDULER_ID=%s:%s\n' "${WATCHDOG_LABELS[${idx}]}" "${status_line#*=}"
          fi
          ;;
      esac
    done < "${status}"
  done
  if [[ -n "${SCHEDULER_FILE:-}" && -s "${SCHEDULER_FILE}" ]]; then
    COMPLETE_SCHEDULER_TMP="${SCHEDULER_FILE}.complete.build.$$"
    : > "${COMPLETE_SCHEDULER_TMP}"
    SCHEDULER_SEEN=""
    append_scheduler_record() {
      local record_kind="$1" record_id="$2" record_key="${record_id}"
      [[ "${record_kind}" =~ ^(ARRAY|WATCHDOG|STATUS|AGGREGATE_GATE|PREFLIGHT)$ &&
         "${record_id}" =~ ^[0-9]+$ ]] ||
        stage5_abort "invalid Stage 5 scheduler record"
      case " ${SCHEDULER_SEEN} " in *" ${record_key} "*) return 0 ;; esac
      SCHEDULER_SEEN="${SCHEDULER_SEEN} ${record_key}"
      printf '%s\t%s\n' "${record_kind}" "${record_id}" >> "${COMPLETE_SCHEDULER_TMP}" ||
        stage5_abort "failed to write Stage 5 scheduler record"
    }
    while IFS=$'\t' read -r record_kind record_id; do
      append_scheduler_record "${record_kind}" "${record_id}"
    done < "${SCHEDULER_FILE}"
    for idx in "${!WATCHDOG_LABELS[@]}"; do
      safe="$(printf '%s' "${WATCHDOG_LABELS[${idx}]}" | tr '/:,\t |' '______')"
      status="${ECODA_RUN_ROOT}/status/watchdogs/${safe}.status"
      while IFS= read -r status_line; do
        case "${status_line}" in
          SCHEDULER_ID=*) append_scheduler_record STATUS "${status_line#*=}" ;;
        esac
      done < "${status}"
    done
    if ! ecoda_atomic_install_manifest "${COMPLETE_SCHEDULER_TMP}" "${SCHEDULER_FILE}" 2; then
      stage5_abort "failed to install complete Stage 5 scheduler ID manifest"
    fi
    rm -f "${COMPLETE_SCHEDULER_TMP}"
  fi
else
  EXACT_SELECTION=0
  exact_line="$(sed -n 's/^EXACT_SELECTION=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  [[ "${exact_line}" == 1 ]] && EXACT_SELECTION=1
  METHODS=(_ecoda_none_)
  ANALYSES=(_ecoda_none_)
  ANALYSES_SELECTED=0
  methods_line="$(sed -n 's/^METHODS=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  if [[ -n "${methods_line}" ]]; then ecoda_split_csv "${methods_line}"; METHODS=("${ECODA_ARRAY[@]}"); fi
  labels_line="$(sed -n 's/^ANALYSES=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  if [[ -n "${labels_line}" ]]; then ecoda_split_csv "${labels_line}"; ANALYSES=("${ECODA_ARRAY[@]}"); ANALYSES_SELECTED=1; fi
  [[ -s "${ECODA_RUN_ROOT}/status/aggregate" ]] &&
    grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/aggregate" ||
    stage5_abort "Stage 5 aggregate status is not OK"
fi

# Exactly one shared merge/checksum/sync owner for this run. Tests stop before
# NAS access; real durable-gate invocations run this tail only after the
# aggregate watchdog status is complete.
if [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]]; then
  ecoda_set_run_state OK "benchmark matrix submission/gate test completed"
  if [[ -n "${PASS_ARG}" ]]; then
    echo "BATCH_EFFECT_RUN_ID=${RUN_ID}"
  else
    echo "BENCHMARK_RUN_ID=${RUN_ID}"
  fi
  exit 0
fi
stage5_prepare_source_identity ||
  stage5_abort "Stage 5 source identity changed before final validation"
LABELS=()
for method in "${METHODS[@]}"; do
  [[ "${method}" == _ecoda_none_ ]] || LABELS+=("${method}")
done
if [[ ${ANALYSES_SELECTED} -eq 1 ]]; then LABELS+=("${ANALYSES[@]}"); fi
validation_args=(--root "${ANALYSIS_ROOT}" --selection "${MANIFEST}" --labels "${LABELS[@]}" \
  --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" \
  --source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
[[ -n "${PASS_ARG}" ]] && validation_args+=(--batch --batch-pass "${PASS_ARG}")
[[ ${EXACT_SELECTION} -eq 1 ]] && validation_args+=(--exact)
if ! "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" "${validation_args[@]}"; then
  stage5_abort "Stage 5 matrix artifact validation failed"
fi
RDS_LABELS=()
for label in "${LABELS[@]}"; do
  case "${label}" in
    gloscope|mofa|pseudobulk|composition|scitd|prepare_pseudobulk|trans|zeroimp)
      RDS_LABELS+=("${label}")
      ;;
  esac
done
if [[ ${#RDS_LABELS[@]} -gt 0 ]]; then
  rds_args=(--root "${ANALYSIS_ROOT}" --selection "${MANIFEST}" \
    --labels "$(IFS=,; echo "${RDS_LABELS[*]}")" \
    --config "${DATASETS_JSON_FILE}" --input-root "${HPC_SCRATCH_DIR}" \
    --source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
  [[ -n "${PASS_ARG}" ]] && rds_args+=(--batch-pass "${PASS_ARG}")
  [[ ${EXACT_SELECTION} -eq 1 ]] && rds_args+=(--exact)
  if ! ${PIXI_RSCRIPT} "${SCRIPT_DIR}/validate_benchmark_rds_contract.R" "${rds_args[@]}"; then
    stage5_abort "Stage 5 RDS artifact validation failed"
  fi
fi
if ! benchmark_merge_sync_cleanup "${LABELS[@]}"; then
  stage5_abort "Stage 5 benchmark synchronization failed"
fi
stage5_finalize_owner_manifest OK "benchmark matrix sync completed" ||
  stage5_abort "failed to finalize Stage 5 owners"
ecoda_set_run_state OK "benchmark matrix aggregate, artifact validation, and selected sync completed" ||
  stage5_abort "failed to write Stage 5 terminal OK state"
if [[ -n "${PASS_ARG}" ]]; then
  echo "BATCH_EFFECT_RUN_ID=${RUN_ID}"
else
  echo "BENCHMARK_RUN_ID=${RUN_ID}"
fi
