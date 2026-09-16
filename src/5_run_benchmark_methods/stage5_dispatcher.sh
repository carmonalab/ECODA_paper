#!/bin/bash
# Canonical Pipeline 5 manifest dispatcher. The public launcher delegates here;
# this implementation owns one run-bound dispatch manifest, scheduler wave,
# aggregate gate, and final merge/checksum/NAS synchronization.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
export ECODA_GATE_STAGE=stage5
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_stage5_policy.sh"
source "${SCRIPT_DIR}/../utils/bash/h5ad_preflight_submit.sh"
cd "${PROJECT_ROOT}"

DATASETS_ARG=""
DATASETS_SET=0
METHODS_ARG=""
METHODS_SET=0
ANALYSES_ARG=""
ANALYSES_SET=0
TARGET_METHODS_ARG=""
TARGET_METHODS_SET=0
TARGET_METHODS=()
SELECTION_FILE_ARG=""
METHOD_MATRIX_ARG=""
METHOD_MATRIX_SET=0
METHOD_MATRIX_MODE=0
METHOD_MATRIX=""
METHOD_MATRIX_SOURCE_PATH=""
METHOD_MATRIX_MD5=""
METHOD_MATRIX_SIZE=""
METHOD_MATRIX_SHA256=""
METHOD_MATRIX_IDENTITY=""
METHOD_MATRIX_DECLARED_COUNT=0
METHOD_MATRIX_PENDING_COUNT=0
METHOD_MATRIX_ROOT_VERSION=""
METHOD_MATRIX_METADATA_MD5=""
METHOD_MATRIX_METADATA_SIZE=""
METHOD_MATRIX_METADATA_SHA256=""
METHOD_MATRIX_METADATA_DECLARED_COUNT=""
METHOD_MATRIX_METADATA_IDENTITY=""
METHOD_MATRIX_METADATA_PENDING_COUNT=""
SELECTION_FILE_SET=0
PASS_ARG=""
ANALYSIS_VARIANT_ARG=""
ANALYSIS_VARIANT_SET=0
PASS_SET=0
EXACT_BATCH_SELECTION=0
FORCE_ARG=0
FORCE_TARGETED_ARG=0
FORCE_REASON_ARG=""
FORCE_REASON_SET=0
SYNC_ONLY_RUN=""
SYNC_ONLY_SET=0
PARTITION_ARG=""
MEMORY="${BENCHMARK_MEM}"
MAX_MEMORY="${BENCHMARK_MEM_MAX}"
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
RUNTIME_EXPORT=""
unset H5AD_ALLOW_MISSING_SUMMARY

GPU_POLICY="auto"
BASELINE_METHODS=(gloscope mofa pseudobulk composition scitd mrvi scpoli pilot qot pilotgm)
STAGE5_INPUT_PRODUCER_RUN_ID="${STAGE5_INPUT_PRODUCER_RUN_ID:-${STAGE4_RUN_ID:-${ANNOTATION_RUN_ID:-${PREPROCESS_RUN_ID:-${STAGE3_RUN_ID:-${INPUT_PRODUCER_RUN_ID:-${ECODA_ARTIFACT_PRODUCER_RUN_ID:-${ECODA_PRODUCER_RUN_ID:-}}}}}}}}"

usage() {
  cat <<'EOF'
Usage: 1_submit_hpc_array.sh [--datasets LIST] [--methods LIST]
       [--target-methods LIST] [--method-matrix TSV]
       [--analyses trans,zeroimp] [--selection-file TSV]
       [--exact-batch-selection] [--pass uncorrected|corrected]
       [--analysis-variant final|corrected_final]
       [--gpu-policy auto|default|any] [--force]
       [--force-targeted --force-reason REASON] [--sync-only RUN_ID]
       [--partition NAME] [--mem VALUE] [--max-mem VALUE] [--throttle N]

Selection-file rows are DATASET<TAB>VIEW<TAB>LABEL. The run-owned
dispatch_selection.tsv records the exact DATASET<TAB>VIEW<TAB>METHOD rows
that the dispatcher may submit; no release scope is inferred from dataset
names or row counts. Ordinary methods use benchmark_analysis; batch mode uses
the selected explicit pass view. --method-matrix supplies an explicit
corrected-final dispatch scope and may contain any non-empty manifest-defined
row set. --target-methods is an explicit selection-file-scoped partial batch
recovery. --force-targeted force-reclaims only those target methods and
requires --target-methods, --pass, --selection-file, and --force-reason.
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
    --method-matrix) METHOD_MATRIX_ARG="${2:-}"; METHOD_MATRIX_SET=1; shift 2 ;;
    --method-matrix=*) METHOD_MATRIX_ARG="${1#*=}"; METHOD_MATRIX_SET=1; shift ;;
    --target-methods) TARGET_METHODS_ARG="${2:-}"; TARGET_METHODS_SET=1; shift 2 ;;
    --target-methods=*) TARGET_METHODS_ARG="${1#*=}"; TARGET_METHODS_SET=1; shift ;;
    --exact-batch-selection) EXACT_BATCH_SELECTION=1; shift ;;
    --pass|--analysis-pass) PASS_ARG="${2:-}"; PASS_SET=1; shift 2 ;;
    --pass=*|--analysis-pass=*) PASS_ARG="${1#*=}"; PASS_SET=1; shift ;;
    --analysis-variant) ANALYSIS_VARIANT_ARG="${2:-}"; ANALYSIS_VARIANT_SET=1; shift 2 ;;
    --analysis-variant=*) ANALYSIS_VARIANT_ARG="${1#*=}"; ANALYSIS_VARIANT_SET=1; shift ;;
    --force) FORCE_ARG=1; shift ;;
    --force-targeted) FORCE_TARGETED_ARG=1; shift ;;
    --force-reason) FORCE_REASON_ARG="${2:-}"; FORCE_REASON_SET=1; shift 2 ;;
    --force-reason=*) FORCE_REASON_ARG="${1#*=}"; FORCE_REASON_SET=1; shift ;;
    --gpu-policy) GPU_POLICY="${2:-}"; shift 2 ;;
    --gpu-policy=*) GPU_POLICY="${1#*=}"; shift ;;
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
if [[ ${FORCE_TARGETED_ARG} -eq 1 ]]; then
  [[ ${FORCE_ARG} -eq 0 ]] || {
    echo "ERROR: --force-targeted cannot be combined with --force." >&2
    exit 1
  }
  [[ ${TARGET_METHODS_SET} -eq 1 && -n "${TARGET_METHODS_ARG}" ]] || {
    echo "ERROR: --force-targeted requires --target-methods." >&2
    exit 1
  }
  [[ ${PASS_SET} -eq 1 && "${PASS_ARG}" == "uncorrected" ]] || {
    echo "ERROR: --force-targeted requires --pass uncorrected." >&2
    exit 1
  }
  [[ ${SELECTION_FILE_SET} -eq 1 && -n "${SELECTION_FILE_ARG}" ]] || {
    echo "ERROR: --force-targeted requires --selection-file." >&2
    exit 1
  }
  [[ ${METHODS_SET} -eq 0 && ${ANALYSES_SET} -eq 0 ]] || {
    echo "ERROR: --force-targeted cannot be combined with --methods or --analyses." >&2
    exit 1
  }
  [[ ${EXACT_BATCH_SELECTION} -eq 0 ]] || {
    echo "ERROR: --force-targeted cannot be combined with --exact-batch-selection." >&2
    exit 1
  }
  [[ ${SYNC_ONLY_SET} -eq 0 ]] || {
    echo "ERROR: --force-targeted cannot be combined with --sync-only." >&2
    exit 1
  }
  [[ ${FORCE_REASON_SET} -eq 1 && -n "${FORCE_REASON_ARG}" ]] || {
    echo "ERROR: --force-targeted requires --force-reason." >&2
    exit 1
  }
  case "${FORCE_REASON_ARG}" in
    *$'\n'*|*$'\r'*|*$'\t'*|*=*|*\\*)
      echo "ERROR: --force-reason contains a record delimiter." >&2
      exit 1
      ;;
  esac
elif [[ ${FORCE_REASON_SET} -eq 1 ]]; then
  echo "ERROR: --force-reason requires --force-targeted." >&2
  exit 1
fi
case "${GPU_POLICY}" in
  auto|default|any) ;;
  *) echo "ERROR: --gpu-policy must be auto, default, or any." >&2; exit 1 ;;
esac
if [[ -n "${BENCHMARK_GPU_ANY_VRAM_PER_GPU}" &&
      ! "${BENCHMARK_GPU_ANY_VRAM_PER_GPU}" =~ ^[0-9]+(G|GB)$ ]]; then
  echo "ERROR: BENCHMARK_GPU_ANY_VRAM_PER_GPU must be an integer G/GB value." >&2
  exit 1
fi
for gpu_value in "${BENCHMARK_GPU_ANY_PARTITION}" \
                 "${BENCHMARK_GPU_DEFAULT_TIME_LIMIT}" \
                 "${BENCHMARK_GPU_ANY_TIME_LIMIT}"; do
  [[ -n "${gpu_value}" && "${gpu_value}" != *$'\n'* && "${gpu_value}" != *' '* ]] || {
    echo "ERROR: GPU resource configuration contains a blank/newline value." >&2
    exit 1
  }
done
if [[ -n "${SYNC_ONLY_RUN}" && ${FORCE_ARG} -eq 1 ]]; then echo "ERROR: --sync-only cannot use --force." >&2; exit 1; fi
if [[ -n "${PASS_ARG}" && "${PASS_ARG}" != uncorrected && "${PASS_ARG}" != corrected ]]; then echo "ERROR: --pass must be uncorrected or corrected." >&2; exit 1; fi
command -v jq >/dev/null 2>&1 || { echo "ERROR: jq is required for benchmark selection." >&2; exit 1; }
mkdir -p "${LOGS_DIR}" || {
  echo "ERROR: could not create Stage 5 log directory." >&2
  exit 1
}
export ECODA_RUNTIME_PROFILE=stage5
RUNTIME_EXPORT=""
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
if [[ ${METHOD_MATRIX_SET} -eq 1 && -z "${METHOD_MATRIX_ARG}" ]]; then
  echo "ERROR: --method-matrix must not be empty." >&2
  exit 1
fi
if [[ ${PASS_SET} -eq 1 && -z "${PASS_ARG}" ]]; then
  echo "ERROR: --pass must not be empty." >&2
  exit 1
fi
if [[ ${ANALYSIS_VARIANT_SET} -eq 1 && -z "${ANALYSIS_VARIANT_ARG}" ]]; then
  echo "ERROR: --analysis-variant must not be empty." >&2
  exit 1
fi
case "${ANALYSIS_VARIANT_ARG}" in
  "") ;;
  final|corrected_final) ;;
  *) echo "ERROR: --analysis-variant must be final or corrected_final." >&2; exit 1 ;;
esac
if [[ ${SYNC_ONLY_SET} -eq 1 && -z "${SYNC_ONLY_RUN}" ]]; then
  echo "ERROR: --sync-only requires a run ID." >&2
  exit 1
fi
EXPECTED_BATCH_METHODS="prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot"
METHOD_MATRIX_ROOT_VERSION="recovery_35row"

# Scope is declared by the caller's manifests.  The dispatcher validates
# syntax, configured dataset/view identity, uniqueness, and method membership;
# it does not contain a release-specific dataset inventory.
stage5_validate_method_matrix_source() {
  local matrix="${METHOD_MATRIX_ARG:-}" row_ds row_view row_method extra
  local selection_ds selection_view selection_label selection_extra
  local actual_md5 actual_size actual_sha256 index=0 seen="" selected=0
  [[ -n "${matrix}" ]] || return 0
  [[ "${ANALYSIS_VARIANT_ARG:-}" == corrected_final &&
     "${PASS_ARG:-}" == corrected && ${PASS_SET} -eq 1 &&
     ${SELECTION_FILE_SET} -eq 1 && -n "${SELECTION_FILE_ARG}" &&
     ${METHODS_SET} -eq 1 &&
     "${METHODS_ARG}" == "${EXPECTED_BATCH_METHODS}" ]] || {
    echo "ERROR: --method-matrix requires corrected_final, corrected pass, an explicit selection file, and the fixed batch method suite." >&2
    return 1
  }
  [[ ${DATASETS_SET} -eq 0 && ${ANALYSES_SET} -eq 0 &&
     ${TARGET_METHODS_SET} -eq 0 && ${EXACT_BATCH_SELECTION} -eq 0 &&
     ${FORCE_ARG} -eq 0 && ${FORCE_TARGETED_ARG} -eq 0 ]] || {
    echo "ERROR: --method-matrix rejects broad, ordinary, exact, targeted, and force modes." >&2
    return 1
  }
  [[ -f "${matrix}" && ! -L "${matrix}" && -r "${matrix}" ]] || {
    echo "ERROR: method matrix is missing, unreadable, or symlinked: ${matrix}" >&2
    return 1
  }
  [[ "${matrix}" != *$'\n'* && "${matrix}" != *$'\r'* &&
     "${matrix}" != *$'\t'* ]] || {
    echo "ERROR: method matrix path contains a record delimiter." >&2
    return 1
  }
  ecoda_validate_manifest "${matrix}" 3 || {
    echo "ERROR: method matrix must be a non-empty headerless three-column TSV." >&2
    return 1
  }
  ecoda_validate_manifest "${SELECTION_FILE_ARG}" 3 || {
    echo "ERROR: method-matrix selection must be a non-empty three-column TSV." >&2
    return 1
  }
  while IFS=$'\t' read -r row_ds row_view row_method extra; do
    [[ -n "${row_ds}" && -n "${row_view}" && -n "${row_method}" &&
       -z "${extra}" ]] || {
      echo "ERROR: method matrix contains an empty or extended row." >&2
      return 1
    }
    [[ "${row_view}" == "batch_effect_corrected" &&
       "${row_ds}" =~ ^[A-Za-z0-9_.-]+$ &&
       "${row_method}" =~ ^[A-Za-z0-9_.-]+$ ]] || {
      echo "ERROR: method matrix contains an unsafe or non-corrected row." >&2
      return 1
    }
    case ",${EXPECTED_BATCH_METHODS}," in
      *,"${row_method}",*) ;;
      *) echo "ERROR: method matrix contains unsupported method ${row_method}." >&2; return 1 ;;
    esac
    ecoda_dataset_exists "${row_ds}" &&
      ecoda_view_exists "${row_ds}" "${row_view}" || {
      echo "ERROR: method matrix references an unknown dataset/view: ${row_ds}/${row_view}" >&2
      return 1
    }
    case " ${seen} " in
      *" ${row_ds}/${row_view}/${row_method} "*)
        echo "ERROR: method matrix contains a duplicate row." >&2
        return 1
        ;;
    esac
    seen="${seen} ${row_ds}/${row_view}/${row_method}"
    selected=0
    while IFS=$'\t' read -r selection_ds selection_view selection_label selection_extra; do
      [[ "${selection_ds}" == "${row_ds}" &&
         "${selection_view}" == "${row_view}" &&
         -z "${selection_extra}" ]] && selected=1
    done < "${SELECTION_FILE_ARG}"
    [[ ${selected} -eq 1 ]] || {
      echo "ERROR: method matrix row is outside the dataset selection: ${row_ds}/${row_view}" >&2
      return 1
    }
    index=$((index + 1))
  done < "${matrix}"
  [[ ${index} -gt 0 ]] || {
    echo "ERROR: method matrix has no rows." >&2
    return 1
  }
  actual_md5="$(ecoda_md5_file "${matrix}")" || return 1
  actual_size="$(wc -c < "${matrix}" | tr -d '[:space:]')" || return 1
  actual_sha256="$(ecoda_sha256_file "${matrix}")" || return 1
  [[ "${actual_md5}" =~ ^[[:xdigit:]]{32}$ &&
     "${actual_size}" =~ ^[1-9][0-9]*$ &&
     "${actual_sha256}" =~ ^[[:xdigit:]]{64}$ ]] || return 1
  METHOD_MATRIX_MODE=1
  METHOD_MATRIX_SOURCE_PATH="${matrix}"
  METHOD_MATRIX_MD5="${actual_md5}"
  METHOD_MATRIX_SIZE="${actual_size}"
  METHOD_MATRIX_SHA256="${actual_sha256}"
  METHOD_MATRIX_IDENTITY="${actual_sha256}"
  METHOD_MATRIX_DECLARED_COUNT=${index}
}

if [[ -n "${PASS_ARG}" && ${METHODS_SET} -eq 1 &&
      ${TARGET_METHODS_SET} -eq 0 &&
      "${METHODS_ARG}" != "${EXPECTED_BATCH_METHODS}" ]]; then
  echo "ERROR: batch-effect pass requires the fixed ordered method suite: ${EXPECTED_BATCH_METHODS}" >&2
  exit 1
fi
STAGE5_MATRIX_METHODS=()
stage5_method_matrix_methods_for() {
  local ds="${1:-}" view="${2:-}" row_ds row_view row_method extra
  local matrix="${METHOD_MATRIX:-${METHOD_MATRIX_SOURCE_PATH:-}}"
  STAGE5_MATRIX_METHODS=()
  [[ ${METHOD_MATRIX_MODE} -eq 1 ]] || return 0
  [[ -s "${matrix}" ]] || return 1
  while IFS=$'\t' read -r row_ds row_view row_method extra; do
    [[ -z "${extra}" ]] || return 1
    if [[ "${row_ds}" == "${ds}" && "${row_view}" == "${view}" ]]; then
      STAGE5_MATRIX_METHODS+=("${row_method}")
    fi
  done < "${matrix}"
  [[ ${#STAGE5_MATRIX_METHODS[@]} -gt 0 ]]
}
if [[ ${TARGET_METHODS_SET} -eq 1 ]]; then
  [[ -n "${TARGET_METHODS_ARG}" ]] || {
    echo "ERROR: --target-methods must not be empty." >&2
    exit 1
  }
  [[ -n "${PASS_ARG}" && ${SELECTION_FILE_SET} -eq 1 ]] || {
    echo "ERROR: --target-methods requires --pass and --selection-file." >&2
    exit 1
  }
  [[ ${EXACT_BATCH_SELECTION} -eq 0 ]] || {
    echo "ERROR: --target-methods cannot be combined with --exact-batch-selection." >&2
    exit 1
  }
  [[ ${METHODS_SET} -eq 0 && ${ANALYSES_SET} -eq 0 ]] || {
    echo "ERROR: --target-methods cannot be combined with --methods or --analyses." >&2
    exit 1
  }
  ecoda_split_csv "${TARGET_METHODS_ARG}" || exit 1
  ecoda_assert_unique_items "${ECODA_ARRAY[@]}" || exit 1
  TARGET_METHODS=("${ECODA_ARRAY[@]}")
  for target_method in "${ECODA_ARRAY[@]}"; do
    case ",${EXPECTED_BATCH_METHODS}," in
      *,"${target_method}",*) ;;
      *) echo "ERROR: unsupported targeted batch-effect method: ${target_method}" >&2; exit 1 ;;
    esac
  done
fi
if [[ -n "${PASS_ARG}" && -n "${SELECTION_FILE_ARG}" ]]; then
  [[ -r "${SELECTION_FILE_ARG}" ]] || { echo "ERROR: selection file is unreadable." >&2; exit 1; }
  ecoda_validate_manifest "${SELECTION_FILE_ARG}" 3 || {
    echo "ERROR: selection file must contain exactly three columns per row." >&2
    exit 1
  }
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
stage5_validate_final_selection() {
  local variant="${ANALYSIS_VARIANT_ARG:-}" expected_pass expected_view
  local row_count=0 ds view label extra final_line row_key matrix_key selected=""
  case "${variant}" in
    "") return 0 ;;
    final)
      expected_pass="uncorrected"
      expected_view="batch_effect_uncorrected"
      ;;
    corrected_final)
      expected_pass="corrected"
      expected_view="batch_effect_corrected"
      ;;
    *)
      echo "ERROR: unsupported Stage 5 analysis variant." >&2
      return 1
      ;;
  esac
  [[ ${PASS_SET} -eq 1 && "${PASS_ARG}" == "${expected_pass}" ]] || {
    echo "ERROR: ${variant} analysis variant requires --pass ${expected_pass}." >&2
    return 1
  }
  [[ ${SELECTION_FILE_SET} -eq 1 && -n "${SELECTION_FILE_ARG}" ]] || {
    echo "ERROR: ${variant} analysis variant requires an explicit --selection-file." >&2
    return 1
  }
  [[ ${DATASETS_SET} -eq 0 && ${ANALYSES_SET} -eq 0 &&
     ${EXACT_BATCH_SELECTION} -eq 0 ]] || {
    echo "ERROR: ${variant} analysis variant rejects broad/default or ordinary selection modes." >&2
    return 1
  }
  if [[ ${FORCE_ARG} -eq 1 ]]; then
    echo "ERROR: ${variant} analysis variant does not support broad --force." >&2
    return 1
  fi
  if [[ ${TARGET_METHODS_SET} -eq 1 ]]; then
    [[ ${METHODS_SET} -eq 0 && ${FORCE_ARG} -eq 0 ]] || {
      echo "ERROR: targeted final recovery cannot include a full method selection or --force." >&2
      return 1
    }
  else
    [[ ${METHODS_SET} -eq 1 &&
       "${METHODS_ARG}" == "${EXPECTED_BATCH_METHODS}" ]] || {
      echo "ERROR: ${variant} analysis variant requires explicit --methods ${EXPECTED_BATCH_METHODS}." >&2
      return 1
    }
  fi
  [[ -r "${SELECTION_FILE_ARG}" ]] || {
    echo "ERROR: ${variant} selection file is unreadable." >&2
    return 1
  }
  ecoda_validate_manifest "${SELECTION_FILE_ARG}" 3 || {
    echo "ERROR: ${variant} selection file must contain exactly three columns per row." >&2
    return 1
  }
  while IFS= read -r final_line || [[ -n "${final_line}" ]]; do
    IFS=$'\t' read -r ds view label extra <<< "${final_line}"
    [[ -n "${ds}" && "${view}" == "${expected_view}" &&
       "${label}" == "${expected_view}" && -z "${extra}" ]] || {
      echo "ERROR: ${variant} selection rows must use ${expected_view}." >&2
      return 1
    }
    row_key="${ds}/${view}"
    case " ${selected} " in
      *" ${row_key} "*) echo "ERROR: ${variant} selection contains duplicate ${row_key}." >&2; return 1 ;;
    esac
    selected="${selected} ${row_key}"
    row_count=$((row_count + 1))
  done < "${SELECTION_FILE_ARG}"
  [[ ${row_count} -gt 0 ]] || {
    echo "ERROR: ${variant} selection is empty." >&2
    return 1
  }
  if [[ ${METHOD_MATRIX_MODE} -eq 1 ]]; then
    while IFS=$'\t' read -r ds view label extra; do
      matrix_key="${ds}/${view}"
      case " ${selected} " in
        *" ${matrix_key} "*) ;;
        *) echo "ERROR: method matrix scope is outside the selection manifest." >&2; return 1 ;;
      esac
    done < "${METHOD_MATRIX_SOURCE_PATH}"
    while IFS=$'\t' read -r ds view label; do
      stage5_method_matrix_methods_for "${ds}" "${view}" || {
        echo "ERROR: selection row has no method scope in the method matrix: ${ds}/${view}" >&2
        return 1
      }
    done < "${SELECTION_FILE_ARG}"
  fi
}
if [[ -n "${SYNC_ONLY_RUN}" ]]; then
  ecoda_validate_run_id "${SYNC_ONLY_RUN}" || exit 1
  sync_metadata="${HPC_SCRATCH_DIR}/_ecoda_runs/${SYNC_ONLY_RUN}/metadata"
  [[ -r "${sync_metadata}" ]] || {
    echo "ERROR: sync-only run metadata is missing or unreadable." >&2
    exit 1
  }
  sync_pass="$(sed -n 's/^PASS=//p' "${sync_metadata}" | head -1 || true)"
  if [[ -n "${sync_pass}" ]]; then
    [[ -z "${PASS_ARG}" || "${PASS_ARG}" == "${sync_pass}" ]] || {
      echo "ERROR: sync-only pass does not match run metadata." >&2
      exit 1
    }
    PASS_ARG="${sync_pass}"
    PASS_SET=1
  fi
  sync_variant="$(sed -n 's/^ANALYSIS_VARIANT=//p' "${sync_metadata}" | head -1 || true)"
  [[ -z "${sync_variant}" || "${sync_variant}" == final ||
     "${sync_variant}" == corrected_final ]] || {
    echo "ERROR: sync-only run has an unknown analysis variant." >&2
    exit 1
  }
  if [[ -n "${sync_variant}" ]]; then
    [[ -z "${ANALYSIS_VARIANT_ARG}" ||
       "${ANALYSIS_VARIANT_ARG}" == "${sync_variant}" ]] || {
      echo "ERROR: sync-only analysis variant does not match run metadata." >&2
      exit 1
    }
    ANALYSIS_VARIANT_ARG="${sync_variant}"
  fi
  stored_root_version="$(sed -n 's/^ANALYSIS_ROOT_VERSION=//p' "${sync_metadata}" | head -1 || true)"
  stored_root_identity="$(sed -n 's/^ANALYSIS_ROOT_IDENTITY=//p' "${sync_metadata}" | head -1 || true)"
  if [[ "${ANALYSIS_VARIANT_ARG:-}" == corrected_final ]]; then
    [[ "${stored_root_version}" == "${METHOD_MATRIX_ROOT_VERSION}" &&
       "${stored_root_identity}" == "corrected_final/${METHOD_MATRIX_ROOT_VERSION}" ]] || {
      echo "ERROR: sync-only corrected-final root-version identity is invalid." >&2
      exit 1
    }
  fi
  if [[ "${ANALYSIS_VARIANT_ARG:-}" == final ||
        "${ANALYSIS_VARIANT_ARG:-}" == corrected_final ]]; then
    expected_selection="${HPC_SCRATCH_DIR}/_ecoda_runs/${SYNC_ONLY_RUN}/manifests/selection.tsv"
    if [[ -z "${SELECTION_FILE_ARG}" ]]; then
      SELECTION_FILE_ARG="${expected_selection}"
      SELECTION_FILE_SET=1
    else
      [[ "${SELECTION_FILE_ARG}" == "${expected_selection}" ]] || {
        echo "ERROR: final sync-only selection is not the run-bound selection." >&2
        exit 1
      }
    fi
    stored_method_matrix="$(sed -n 's/^METHOD_MATRIX=//p' "${sync_metadata}" | head -1 || true)"
    if [[ -n "${stored_method_matrix}" ]]; then
      expected_method_matrix="${ECODA_RUNS_ROOT}/${SYNC_ONLY_RUN}/manifests/method_matrix.tsv"
      [[ "${stored_method_matrix}" == "${expected_method_matrix}" ]] || {
        echo "ERROR: sync-only method matrix is not the run-bound manifest." >&2
        exit 1
      }
      if [[ ${METHOD_MATRIX_SET} -eq 1 ]]; then
        [[ "${METHOD_MATRIX_ARG}" == "${stored_method_matrix}" ]] || {
          echo "ERROR: sync-only method matrix does not match run metadata." >&2
          exit 1
        }
      else
        METHOD_MATRIX_ARG="${stored_method_matrix}"
        METHOD_MATRIX_SET=1
      fi
      METHOD_MATRIX_METADATA_MD5="$(sed -n 's/^METHOD_MATRIX_MD5=//p' "${sync_metadata}" | head -1 || true)"
      METHOD_MATRIX_METADATA_SIZE="$(sed -n 's/^METHOD_MATRIX_SIZE=//p' "${sync_metadata}" | head -1 || true)"
      METHOD_MATRIX_METADATA_SHA256="$(sed -n 's/^METHOD_MATRIX_SHA256=//p' "${sync_metadata}" | head -1 || true)"
      METHOD_MATRIX_METADATA_DECLARED_COUNT="$(sed -n 's/^DECLARED_METHOD_ROWS=//p' "${sync_metadata}" | head -1 || true)"
      METHOD_MATRIX_METADATA_IDENTITY="$(sed -n 's/^METHOD_MATRIX_IDENTITY=//p' "${sync_metadata}" | head -1 || true)"
      METHOD_MATRIX_METADATA_PENDING_COUNT="$(sed -n 's/^PENDING_METHOD_ROWS=//p' "${sync_metadata}" | head -1 || true)"
    fi
    stored_methods="$(sed -n 's/^METHODS=//p' "${sync_metadata}" | head -1 || true)"
    stored_target_methods="$(sed -n 's/^TARGET_METHODS=//p' "${sync_metadata}" | head -1 || true)"
    if [[ -n "${stored_target_methods}" ]]; then
      [[ ${METHODS_SET} -eq 0 &&
         ( ${TARGET_METHODS_SET} -eq 0 ||
           "${TARGET_METHODS_ARG}" == "${stored_target_methods}" ) ]] || {
        echo "ERROR: final sync-only target methods do not match run metadata." >&2
        exit 1
      }
      TARGET_METHODS_ARG="${stored_target_methods}"
      TARGET_METHODS_SET=1
      METHODS_ARG=""
      METHODS_SET=0
      ecoda_split_csv "${TARGET_METHODS_ARG}" || exit 1
      TARGET_METHODS=("${ECODA_ARRAY[@]}")
      ecoda_assert_unique_items "${TARGET_METHODS[@]}" || exit 1
      for stored_target_method in "${TARGET_METHODS[@]}"; do
        case ",${EXPECTED_BATCH_METHODS}," in
          *,"${stored_target_method}",*) ;;
          *) echo "ERROR: final sync-only metadata contains an unsupported target method." >&2; exit 1 ;;
        esac
      done
    elif [[ -z "${METHODS_ARG}" ]]; then
      METHODS_ARG="${stored_methods}"
      METHODS_SET=1
    fi
  fi
fi
stage5_validate_method_matrix_source || exit 1
stage5_validate_final_selection || exit 1


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


# Corrected configuration is a validator-only boundary.  Resolve the same
# dataset scope that the later selection builder will use, but do not create
# run, selection, scheduler, worker, or preflight state while checking it.
stage5_validate_corrected_batch_selection() {
  local dataset selection_manifest=""
  local old_ifs contract_method
  local -a corrected_datasets=()
  local -a corrected_contract_methods=(
    preprocess prepare_pseudobulk pseudobulk gloscope composition mrvi pilot qot
  )

  [[ -r "${DATASETS_JSON_FILE}" ]] || {
    _ecoda_die "corrected batch configuration is unreadable: ${DATASETS_JSON_FILE}"
    return 1
  }
  jq -e 'type == "object"' "${DATASETS_JSON_FILE}" >/dev/null 2>&1 || {
    _ecoda_die "corrected batch configuration is not a JSON object"
    return 1
  }

  if [[ -n "${SELECTION_FILE_ARG}" ]]; then
    while IFS=$'\t' read -r dataset _selection_view _selection_label; do
      [[ -n "${dataset}" ]] || {
        _ecoda_die "corrected selection contains an empty dataset"
        return 1
      }
      corrected_datasets+=("${dataset}")
    done < "${SELECTION_FILE_ARG}"
  elif [[ -n "${DATASETS_ARG}" ]]; then
    [[ "${DATASETS_ARG}" != ,* && "${DATASETS_ARG}" != *, &&
       "${DATASETS_ARG}" != *,,* ]] || {
      _ecoda_die "corrected dataset selection contains an empty item"
      return 1
    }
    old_ifs="${IFS}"
    IFS=','
    read -r -a corrected_datasets <<< "${DATASETS_ARG}"
    IFS="${old_ifs}"
  elif [[ -n "${SYNC_ONLY_RUN}" ]]; then
    selection_manifest="${ECODA_RUNS_ROOT}/${SYNC_ONLY_RUN}/manifests/selection.tsv"
    if [[ -r "${selection_manifest}" ]]; then
      while IFS=$'\t' read -r dataset _selection_view _selection_label; do
        [[ -n "${dataset}" ]] || {
          _ecoda_die "corrected run selection contains an empty dataset"
          return 1
        }
        corrected_datasets+=("${dataset}")
      done < "${selection_manifest}"
    fi
  fi

  if [[ ${#corrected_datasets[@]} -eq 0 ]]; then
    while IFS= read -r dataset; do
      [[ -n "${dataset}" ]] && corrected_datasets+=("${dataset}")
    done < <(
      jq -r '
        to_entries[]
        | select((.key | startswith("_") | not)
                 and (.value | type == "object")
                 and (.value.use_for_batch_effect == true))
        | .key
      ' "${DATASETS_JSON_FILE}"
    )
  fi
  [[ ${#corrected_datasets[@]} -gt 0 ]] || {
    _ecoda_die "corrected batch selection contains no datasets"
    return 1
  }

  for dataset in "${corrected_datasets[@]}"; do
    ecoda_validate_corrected_batch_columns \
      "${DATASETS_JSON_FILE}" "${dataset}" batch_effect_corrected || return 1
    for contract_method in "${corrected_contract_methods[@]}"; do
      ecoda_corrected_batch_method_policy "${contract_method}" || return 1
      ecoda_batch_contract_identity \
        "${DATASETS_JSON_FILE}" "${dataset}" batch_effect_corrected \
        "${ECODA_CORRECTED_BATCH_METHOD_ID}" \
        "${ECODA_CORRECTED_BATCH_MODEL_ID}" >/dev/null || {
        _ecoda_die "invalid corrected batch contract identity for ${dataset}/${contract_method}"
        return 1
      }
    done
  done
}

if [[ "${PASS_ARG}" == corrected ]]; then
  stage5_validate_corrected_batch_selection || exit 1
fi

# Pass-sensitive helper paths are resolved only after --pass and
# --analysis-variant have been parsed and validated.
stage5_configure_analysis_context() {
  case "${ANALYSIS_VARIANT_ARG:-}" in
    final)
      [[ "${PASS_ARG}" == uncorrected ]] || {
        echo "ERROR: final analysis variant requires uncorrected pass." >&2
        return 1
      }
      unset ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION
      export ANALYSIS_VARIANT=final
      ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final"
      ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/uncorrected_final"
      ANALYSIS_LOG_PREFIX="execution_times_batch_effect_uncorrected_final_"
      ;;
    corrected_final)
      [[ "${PASS_ARG}" == corrected ]] || {
        echo "ERROR: corrected_final analysis variant requires corrected pass." >&2
        return 1
      }
      export ANALYSIS_VARIANT=corrected_final
      export ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION="${METHOD_MATRIX_ROOT_VERSION}"
      ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/corrected_final/${METHOD_MATRIX_ROOT_VERSION}"
      ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/corrected_final/${METHOD_MATRIX_ROOT_VERSION}"
      ANALYSIS_LOG_PREFIX="execution_times_batch_effect_corrected_final_"
      ;;
    "")
      unset ANALYSIS_VARIANT ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION
      if [[ -n "${PASS_ARG}" ]]; then
        ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/${PASS_ARG}"
        ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/${PASS_ARG}"
        ANALYSIS_LOG_PREFIX="execution_times_batch_effect_${PASS_ARG}_"
      else
        ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/benchmark"
        ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/benchmark"
        ANALYSIS_LOG_PREFIX="execution_times_"
      fi
      ;;
    *)
      echo "ERROR: unsupported analysis variant." >&2
      return 1
      ;;
  esac
  export ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_LOG_PREFIX ROOT="${ANALYSIS_ROOT}"
  if [[ -n "${PASS_ARG}" ]]; then
    export ANALYSIS_PASS="${PASS_ARG}"
  else
    unset ANALYSIS_PASS
  fi
}
stage5_configure_analysis_context || exit 1
if [[ -n "${PASS_ARG}" ]]; then
  unset BENCHMARK_MANIFEST
fi
source "${SCRIPT_DIR}/benchmark_submit_common.sh"

stage5_bind_method_matrix() {
  local destination="${ECODA_RUN_ROOT}/manifests/method_matrix.tsv"
  local source="${METHOD_MATRIX_SOURCE_PATH:-${METHOD_MATRIX_ARG:-}}"
  local source_sha256 destination_sha256 pending_rows
  [[ ${METHOD_MATRIX_MODE} -eq 1 ]] || {
    unset ECODA_STAGE5_METHOD_MATRIX METHOD_MATRIX
    return 0
  }
  if [[ -n "${SYNC_ONLY_RUN}" ]]; then
    [[ "${METHOD_MATRIX_ARG}" == "${destination}" ]] || return 1
    [[ -f "${destination}" && ! -L "${destination}" && -r "${destination}" ]] || return 1
  else
    [[ -f "${source}" && ! -L "${source}" && -r "${source}" ]] || return 1
    source_sha256="$(ecoda_sha256_file "${source}")" || return 1
    ecoda_atomic_install_manifest "${source}" "${destination}" 3 || return 1
    destination_sha256="$(ecoda_sha256_file "${destination}")" || return 1
    [[ "${source_sha256}" == "${destination_sha256}" ]] || {
      echo "ERROR: method matrix changed while being copied." >&2
      return 1
    }
  fi
  if [[ -z "${SYNC_ONLY_RUN}" ]]; then
    ecoda_write_checksum "${destination}" || return 1
  fi
  ecoda_validate_run_owned_path "${destination}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_manifest "${destination}" 3 || return 1
  ecoda_validate_checksum "${destination}" || return 1
  METHOD_MATRIX="${destination}"
  METHOD_MATRIX_MD5="${ECODA_CHECKSUM_MD5}"
  METHOD_MATRIX_SIZE="${ECODA_CHECKSUM_SIZE}"
  METHOD_MATRIX_SHA256="$(ecoda_sha256_file "${destination}")" || return 1
  METHOD_MATRIX_IDENTITY="${METHOD_MATRIX_SHA256}"
  METHOD_MATRIX_DECLARED_COUNT="$(awk 'END { print NR }' "${destination}")" || return 1
  export METHOD_MATRIX ECODA_STAGE5_METHOD_MATRIX="${METHOD_MATRIX}"
  export ECODA_STAGE5_METHOD_MATRIX_MODE=1
  if [[ -n "${SYNC_ONLY_RUN}" ]]; then
    [[ "${METHOD_MATRIX_METADATA_MD5}" =~ ^[[:xdigit:]]{32}$ &&
       "${METHOD_MATRIX_METADATA_SIZE}" =~ ^[1-9][0-9]*$ &&
       "${METHOD_MATRIX_METADATA_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
       "${METHOD_MATRIX_METADATA_IDENTITY}" =~ ^[[:xdigit:]]{64}$ &&
       "${METHOD_MATRIX_METADATA_DECLARED_COUNT}" =~ ^[1-9][0-9]*$ &&
       "${METHOD_MATRIX_METADATA_PENDING_COUNT}" =~ ^[0-9]+$ ]] || return 1
  fi
  [[ -z "${METHOD_MATRIX_METADATA_MD5}" ||
     "${METHOD_MATRIX_METADATA_MD5}" == "${METHOD_MATRIX_MD5}" ]] || return 1
  [[ -z "${METHOD_MATRIX_METADATA_SIZE}" ||
     "${METHOD_MATRIX_METADATA_SIZE}" == "${METHOD_MATRIX_SIZE}" ]] || return 1
  [[ -z "${METHOD_MATRIX_METADATA_SHA256}" ||
     "${METHOD_MATRIX_METADATA_SHA256}" == "${METHOD_MATRIX_SHA256}" ]] || return 1
  [[ -z "${METHOD_MATRIX_METADATA_DECLARED_COUNT}" ||
     "${METHOD_MATRIX_METADATA_DECLARED_COUNT}" == "${METHOD_MATRIX_DECLARED_COUNT}" ]] || return 1
  [[ -z "${METHOD_MATRIX_METADATA_IDENTITY}" ||
     "${METHOD_MATRIX_METADATA_IDENTITY}" == "${METHOD_MATRIX_IDENTITY}" ]] || return 1
  if [[ -n "${METHOD_MATRIX_METADATA_PENDING_COUNT}" ]]; then
    pending_rows=0
    if [[ -r "${ECODA_RUN_ROOT}/manifests/pending_selection.tsv" ]]; then
      pending_rows="$(awk 'END { print NR }' \
        "${ECODA_RUN_ROOT}/manifests/pending_selection.tsv")" || return 1
    fi
    [[ "${METHOD_MATRIX_METADATA_PENDING_COUNT}" =~ ^[0-9]+$ &&
       "${pending_rows}" == "${METHOD_MATRIX_METADATA_PENDING_COUNT}" ]] || return 1
  fi
}

stage5_source_script() {
  local relative="$1"
  local source_root="${ECODA_SOURCE_ROOT:-${PROJECT_ROOT}}"
  printf '%s/%s' "${source_root%/}" "${relative}"
}

stage5_require_source_script() {
  local candidate="$1" source_root="${ECODA_SOURCE_ROOT:-}"
  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    [[ -n "${source_root}" ]] || return 1
    ecoda_require_source_script_path "${candidate}" "${source_root}"
  else
    [[ -f "${candidate}" && -r "${candidate}" && ! -L "${candidate}" ]]
  fi
}

stage5_install_source_manifest() {
  local incoming="${ECODA_SOURCE_MANIFEST:-}"
  local destination="${ECODA_RUN_ROOT}/manifests/source.manifest"
  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" != "1" ]]; then
    [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0


    echo "ERROR: new Stage 5 runs require an immutable source snapshot." >&2
    return 1
  fi
  ecoda_atomic_copy "${incoming}" "${destination}" || return 1
  ecoda_validate_run_owned_path "${destination}" "${ECODA_RUN_ROOT}" || return 1
  [[ "$(_ecoda_runtime_require_manifest_value "${destination}" FORMAT)" == "1" ]] || return 1
  [[ "$(_ecoda_runtime_require_manifest_value "${destination}" SOURCE_ROOT)" == "${ECODA_SOURCE_ROOT}" ]] ||
    return 1
}
stage5_validate_output_ownership() {
  local manifest="$1" reclaim_terminal="${2:-0}"
  local ownership_manifest="${1}" ownership_tmp="" rc
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  [[ "${reclaim_terminal}" == "0" || "${reclaim_terminal}" == "1" ]] || return 1
  if awk -F '\t' '
      NF != 3 { saw_extended=1 }
      END { exit(saw_extended ? 0 : 1) }
    ' "${manifest}"; then
    ownership_tmp="${manifest}.ownership.$$"
    awk -F '\t' 'NF >= 3 { print $1 "\t" $2 "\t" $3 }' \
      "${manifest}" > "${ownership_tmp}" || return 1
    ownership_manifest="${ownership_tmp}"
  fi
  ecoda_stage5_validate_output_ownership "${ownership_manifest}" "${RUN_ID}" \
    "${reclaim_terminal}"
  rc=$?
  [[ -n "${ownership_tmp}" ]] && rm -f "${ownership_tmp}"
  return "${rc}"
}

stage5_bind_run_identity() {
  local source_copy="${ECODA_RUN_ROOT}/manifests/source.manifest"
  local runtime_identity="${ECODA_RUN_ROOT}/manifests/runtime.identity"
  local source_root source_manifest image manifest
  [[ -s "${source_copy}" && -s "${runtime_identity}" ]] || {
    echo "legacy_source_unpinned" >&2
    return 1
  }
  ecoda_validate_run_owned_path "${source_copy}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_run_owned_path "${runtime_identity}" "${ECODA_RUN_ROOT}" || return 1
  source_root="$(_ecoda_runtime_require_manifest_value "${source_copy}" SOURCE_ROOT)" || return 1
  [[ "${source_root}" = /* ]] || return 1
  source_manifest="${source_root%/tree}/identity/source.manifest"
  [[ -f "${source_manifest}" && -r "${source_manifest}" ]] || return 1
  cmp -s "${source_copy}" "${source_manifest}" || {
    echo "ERROR: run-owned source manifest differs from the immutable snapshot manifest." >&2
    return 1
  }
  image="$(_ecoda_runtime_require_manifest_value "${runtime_identity}" RUNTIME_IMAGE)" || return 1
  manifest="$(_ecoda_runtime_require_manifest_value "${runtime_identity}" RUNTIME_MANIFEST)" || return 1
  [[ "${image}" = /* && "${manifest}" = /* ]] || return 1
  export ECODA_SOURCE_ROOT="${source_root}"
  export ECODA_SOURCE_MANIFEST="${source_manifest}"
  export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
  export ECODA_AUX_ROOT="${source_root%/}/aux"
  export ECODA_RUNTIME_IMAGE="${image}"
  export ECODA_RUNTIME_MANIFEST="${manifest}"
  export ECODA_RUN_ID="${RUN_ID}"
  PROJECT_ROOT="${source_root}"
  DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
  export PROJECT_ROOT DATASETS_JSON_FILE
  SCRIPT_DIR="${PROJECT_ROOT}/src/5_run_benchmark_methods"
  ANALYSIS_MERGE_SCRIPT="${SCRIPT_DIR}/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py"
}

stage5_validate_bound_runtime() {
  if [[ -n "${ECODA_RUN_ROOT:-}" ||
        "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    ecoda_runtime_validate_bound_run
  else
    ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE:-host}"
  fi
}

stage5_record_run_identity_metadata() {
  local source_copy="${ECODA_RUN_ROOT}/manifests/source.manifest"
  local runtime_identity="${ECODA_RUN_ROOT}/manifests/runtime.identity"
  local key value
  [[ -s "${source_copy}" && -s "${runtime_identity}" ]] || return 0
  printf 'SOURCE_MANIFEST=%s\nRUNTIME_IDENTITY=%s\n' "${source_copy}" "${runtime_identity}"
  for key in FORMAT SOURCE_ROOT SOURCE_COMMIT SOURCE_ARCHIVE_PATH SOURCE_ARCHIVE_SHA256 \
    CONFIG_HELPER_SHA256 DATASETS_SHA256 PIXI_TOML_SHA256 PIXI_LOCK_SHA256 AUX_ROOT \
    SCGATE_DB_BRANCH; do
    value="$(_ecoda_runtime_require_manifest_value "${source_copy}" "${key}")" || return 1
    printf '%s=%s\n' "SOURCE_${key}" "${value}"
  done
  while IFS='=' read -r key value; do
    [[ -n "${key}" && -n "${value}" ]] || continue
    printf '%s=%s\n' "${key}" "${value}"
  done < "${runtime_identity}"
}

ECODA_BATCH_CONTRACT_MANIFEST=""
ECODA_BATCH_CONTRACT_MANIFEST_MD5=""
ECODA_BATCH_CONTRACT_MANIFEST_SIZE=""
ECODA_BATCH_CONTRACT_MANIFEST_SHA256=""

stage5_batch_contract_identity_path() {
  local ds="${1:-}" view="${2:-}" method="${3:-}"
  local manifest="${ECODA_BATCH_CONTRACT_MANIFEST:-}"
  local row_ds row_view row_method row_path row_md5 row_size extra found=0
  [[ -n "${manifest}" && -r "${manifest}" &&
     -n "${ds}" && -n "${view}" && -n "${method}" ]] || return 1
  while IFS=$'\t' read -r row_ds row_view row_method row_path row_md5 row_size extra; do
    [[ -z "${extra}" ]] || return 1
    if [[ "${row_ds}" == "${ds}" && "${row_view}" == "${view}" &&
          "${row_method}" == "${method}" ]]; then
      [[ ${found} -eq 0 ]] || return 1
      found=1
      printf '%s' "${row_path}"
    fi
  done < "${manifest}"
  [[ ${found} -eq 1 ]]
}
stage5_build_batch_contract_identity() {
  local ds="${1:-}" view="${2:-}" method_id="${3:-}" model_id="${4:-}"
  local source_path="${5:-}"
  [[ -n "${ds}" && -n "${view}" && -n "${method_id}" &&
     -n "${model_id}" ]] || return 1
  if [[ "${ANALYSIS_VARIANT_ARG:-${ANALYSIS_VARIANT:-}}" == corrected_final ]]; then
    ecoda_batch_contract_identity \
      "${DATASETS_JSON_FILE}" "${ds}" "${view}" "${method_id}" \
      "${model_id}"
  else
    [[ -n "${source_path}" ]] || return 1
    ecoda_batch_contract_identity \
      "${DATASETS_JSON_FILE}" "${ds}" "${view}" "${method_id}" \
      "${model_id}" "${source_path}"
  fi
}


stage5_record_batch_contract_metadata() {
  [[ "${PASS_ARG:-}" == corrected ]] || return 0
  [[ -n "${ECODA_BATCH_CONTRACT_MANIFEST}" &&
     "${ECODA_BATCH_CONTRACT_MANIFEST_MD5}" =~ ^[[:xdigit:]]{32}$ &&
     "${ECODA_BATCH_CONTRACT_MANIFEST_SIZE}" =~ ^[1-9][0-9]*$ &&
     "${ECODA_BATCH_CONTRACT_MANIFEST_SHA256}" =~ ^[[:xdigit:]]{64}$ ]] || {
    echo "ERROR: corrected Stage 5 batch-contract manifest metadata is incomplete." >&2
    return 1
  }
  printf 'BATCH_CONTRACT_MANIFEST=%s\nBATCH_CONTRACT_MANIFEST_MD5=%s\nBATCH_CONTRACT_MANIFEST_SIZE=%s\nBATCH_CONTRACT_MANIFEST_SHA256=%s\n' \
    "${ECODA_BATCH_CONTRACT_MANIFEST}" \
    "${ECODA_BATCH_CONTRACT_MANIFEST_MD5}" \
    "${ECODA_BATCH_CONTRACT_MANIFEST_SIZE}" \
    "${ECODA_BATCH_CONTRACT_MANIFEST_SHA256}"
}


stage5_validate_batch_contract_manifest() {
  local manifest="${ECODA_BATCH_CONTRACT_MANIFEST:-}"
  local row_ds row_view row_method row_path row_md5 row_size extra
  local source_path
  local expected_path expected_identity actual_count=0 duplicate_key=""
  local safe key found expected_count=0
  local -a actual_ds=() actual_views=() actual_methods=() actual_paths=()
  local -a contract_methods=(preprocess)
  local contract_method ds view method
  local -a expected_contract_methods=()
  [[ "${PASS_ARG:-}" == corrected ]] || return 0
  [[ -n "${manifest}" && "${manifest}" == "${ECODA_RUN_ROOT}/manifests/batch_contract.tsv" &&
     -f "${manifest}" && ! -L "${manifest}" && -r "${manifest}" ]] || {
    echo "ERROR: corrected Stage 5 batch-contract manifest is missing or unsafe." >&2
    return 1
  }
  ecoda_validate_run_owned_path "${manifest}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_manifest "${manifest}" 6 || return 1
  ecoda_validate_checksum "${manifest}" || return 1
  [[ "${ECODA_BATCH_CONTRACT_MANIFEST_MD5}" == "${ECODA_CHECKSUM_MD5}" &&
     "${ECODA_BATCH_CONTRACT_MANIFEST_SIZE}" == "${ECODA_CHECKSUM_SIZE}" ]] || {
    echo "ERROR: corrected Stage 5 batch-contract manifest checksum metadata mismatches." >&2
    return 1
  }
  [[ "${ECODA_BATCH_CONTRACT_MANIFEST_SHA256}" == "$(ecoda_sha256_file "${manifest}")" ]] || {
    echo "ERROR: corrected Stage 5 batch-contract manifest SHA-256 mismatches." >&2
    return 1
  }
  for method in "${METHODS[@]:-}"; do
    [[ "${method}" == _ecoda_none_ ]] || contract_methods+=("${method}")
  done
  while IFS=$'\t' read -r row_ds row_view row_method row_path row_md5 row_size extra; do
    [[ -n "${row_ds}" && -n "${row_view}" && -n "${row_method}" &&
       -n "${row_path}" && -n "${row_md5}" && -n "${row_size}" &&
       -z "${extra}" &&
       "${row_ds}" != *$'\n'* && "${row_ds}" != *$'\t'* &&
       "${row_view}" != *$'\n'* && "${row_view}" != *$'\t'* &&
       "${row_method}" != *$'\n'* && "${row_method}" != *$'\t'* &&
       "${row_path}" = /* && "${row_path}" != *$'\n'* &&
       "${row_path}" != *$'\t'* &&
       "${row_md5}" =~ ^[[:xdigit:]]{32}$ &&
       "${row_size}" =~ ^[1-9][0-9]*$ ]] || return 1
    key="${row_ds}|${row_view}|${row_method}"
    case " ${duplicate_key} " in
      *" ${key} "*)
        echo "ERROR: corrected Stage 5 batch-contract manifest has duplicate key ${key}." >&2
        return 1
        ;;
    esac
    duplicate_key="${duplicate_key} ${key}"
    safe="$(_ecoda_safe_component "${row_ds}__${row_view}__${row_method}")" || return 1
    expected_path="${ECODA_RUN_ROOT}/manifests/batch_contracts/${safe}.json"
    [[ "${row_path}" == "${expected_path}" ]] || return 1
    ecoda_validate_run_owned_path "${row_path}" "${ECODA_RUN_ROOT}" || return 1
    ecoda_validate_checksum "${row_path}" || return 1
    [[ "${row_md5}" == "${ECODA_CHECKSUM_MD5}" &&
       "${row_size}" == "${ECODA_CHECKSUM_SIZE}" ]] || return 1
    ecoda_corrected_batch_method_policy "${row_method}" || return 1
    source_path="$(stage5_input_path "${row_ds}" "${row_view}")" || return 1
    if [[ "${BENCHMARK_MATRIX_TEST:-0}" != "1" ]]; then
      [[ -s "${source_path}" ]] || return 1
    fi
    expected_identity="$(
      stage5_build_batch_contract_identity "${row_ds}" "${row_view}" \
        "${ECODA_CORRECTED_BATCH_METHOD_ID}" \
        "${ECODA_CORRECTED_BATCH_MODEL_ID}" "${source_path}"
    )" || return 1

    cmp -s "${row_path}" <(printf '%s\n' "${expected_identity}") || {
      echo "ERROR: corrected Stage 5 batch-contract identity mismatches ${key}." >&2
      return 1
    }
    actual_ds+=("${row_ds}")
    actual_views+=("${row_view}")
    actual_methods+=("${row_method}")
    actual_paths+=("${row_path}")
    actual_count=$((actual_count + 1))
  done < "${manifest}"
  while IFS=$'\t' read -r ds view _row_label; do
    [[ -n "${ds}" && -n "${view}" ]] || return 1
    expected_contract_methods=(preprocess)
    if [[ ${METHOD_MATRIX_MODE} -eq 1 ]]; then
      stage5_method_matrix_methods_for "${ds}" "${view}" || return 1
      expected_contract_methods+=("${STAGE5_MATRIX_METHODS[@]}")
    else
      expected_contract_methods=("${contract_methods[@]}")
    fi
    for contract_method in "${expected_contract_methods[@]}"; do
      found=0
      for index in "${!actual_ds[@]}"; do
        if [[ "${actual_ds[${index}]}" == "${ds}" &&
              "${actual_views[${index}]}" == "${view}" &&
              "${actual_methods[${index}]}" == "${contract_method}" ]]; then
          found=1
          break
        fi
      done
      [[ ${found} -eq 1 ]] || {
        echo "ERROR: corrected Stage 5 batch-contract identity is missing ${ds}/${view}/${contract_method}." >&2
        return 1
      }
      expected_count=$((expected_count + 1))
    done
  done < "${MANIFEST}"
  [[ ${actual_count} -eq ${expected_count} ]] || {
    echo "ERROR: corrected Stage 5 batch-contract manifest has unexpected rows." >&2
    return 1
  }
}

stage5_create_batch_contract_manifest() {
  local manifest="${ECODA_RUN_ROOT}/manifests/batch_contract.tsv"
  local contracts_dir="${ECODA_RUN_ROOT}/manifests/batch_contracts"
  local manifest_tmp="${manifest}.build.$$"
  local identity_path identity_tmp identity_json identity_md5 identity_size
  local source_path
  local ds view row_label contract_method safe fingerprint key
  local seen_views="" seen_contracts=""
  local -a contract_methods=(preprocess)
  [[ "${PASS_ARG:-}" == corrected ]] || return 0
  local -a contract_methods_for=()
  for method in "${METHODS[@]:-}"; do
    [[ "${method}" == _ecoda_none_ ]] || contract_methods+=("${method}")
  done
  if [[ -n "${SYNC_ONLY_RUN:-}" ]]; then
    ECODA_BATCH_CONTRACT_MANIFEST="$(
      sed -n 's/^BATCH_CONTRACT_MANIFEST=//p' \
        "${ECODA_RUN_ROOT}/metadata" | head -1 || true
    )"
    ECODA_BATCH_CONTRACT_MANIFEST_MD5="$(
      sed -n 's/^BATCH_CONTRACT_MANIFEST_MD5=//p' \
        "${ECODA_RUN_ROOT}/metadata" | head -1 || true
    )"
    ECODA_BATCH_CONTRACT_MANIFEST_SIZE="$(
      sed -n 's/^BATCH_CONTRACT_MANIFEST_SIZE=//p' \
        "${ECODA_RUN_ROOT}/metadata" | head -1 || true
    )"
    ECODA_BATCH_CONTRACT_MANIFEST_SHA256="$(
      sed -n 's/^BATCH_CONTRACT_MANIFEST_SHA256=//p' \
        "${ECODA_RUN_ROOT}/metadata" | head -1 || true
    )"
    stage5_validate_batch_contract_manifest
    return $?
  fi
  mkdir -p "${contracts_dir}" || return 1
  : > "${manifest_tmp}" || return 1
  while IFS=$'\t' read -r ds view row_label; do
    key="${ds}|${view}"
    case " ${seen_views} " in
      *" ${key} "*) continue ;;
    esac
    source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
    if [[ "${BENCHMARK_MATRIX_TEST:-0}" != "1" ]]; then
      [[ -s "${source_path}" ]] || return 1
    fi
    contract_methods_for=(preprocess)
    if [[ ${METHOD_MATRIX_MODE} -eq 1 ]]; then
      stage5_method_matrix_methods_for "${ds}" "${view}" || return 1
      contract_methods_for+=("${STAGE5_MATRIX_METHODS[@]}")
    else
      contract_methods_for=("${contract_methods[@]}")
    fi
    for contract_method in "${contract_methods_for[@]}"; do
      key="${ds}|${view}|${contract_method}"
      case " ${seen_contracts} " in
        *" ${key} "*) return 1 ;;
      esac
      seen_contracts="${seen_contracts} ${key}"
      ecoda_corrected_batch_method_policy "${contract_method}" || return 1
      identity_json="$(
        stage5_build_batch_contract_identity "${ds}" "${view}" \
          "${ECODA_CORRECTED_BATCH_METHOD_ID}" \
          "${ECODA_CORRECTED_BATCH_MODEL_ID}" "${source_path}"
      )" || return 1
      [[ -n "${identity_json}" ]] || return 1
      safe="$(_ecoda_safe_component "${ds}__${view}__${contract_method}")" || return 1
      identity_path="${contracts_dir}/${safe}.json"
      identity_tmp="${identity_path}.build.$$"
      printf '%s\n' "${identity_json}" > "${identity_tmp}" || {
        rm -f "${identity_tmp}"
        return 1
      }
      chmod 600 "${identity_tmp}" || {
        rm -f "${identity_tmp}"
        return 1
      }
      mv -f "${identity_tmp}" "${identity_path}" || {
        rm -f "${identity_tmp}"
        return 1
      }
      ecoda_write_checksum "${identity_path}" || return 1
      identity_md5="${ECODA_CHECKSUM_MD5}"
      identity_size="${ECODA_CHECKSUM_SIZE}"
      fingerprint="$(jq -er '.fingerprint' "${identity_path}")" || return 1
      [[ "${fingerprint}" =~ ^[[:xdigit:]]{64}$ ]] || return 1
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${ds}" "${view}" "${contract_method}" "${identity_path}" \
        "${identity_md5}" "${identity_size}" >> "${manifest_tmp}" || return 1
    done
  done < "${MANIFEST}"
  ecoda_atomic_install_manifest "${manifest_tmp}" "${manifest}" 6 || {
    rm -f "${manifest_tmp}"
    return 1
  }
  rm -f "${manifest_tmp}"
  ecoda_write_checksum "${manifest}" || return 1
  ECODA_BATCH_CONTRACT_MANIFEST="${manifest}"
  ECODA_BATCH_CONTRACT_MANIFEST_MD5="${ECODA_CHECKSUM_MD5}"
  ECODA_BATCH_CONTRACT_MANIFEST_SIZE="${ECODA_CHECKSUM_SIZE}"
  ECODA_BATCH_CONTRACT_MANIFEST_SHA256="$(ecoda_sha256_file "${manifest}")" || return 1
  stage5_validate_batch_contract_manifest
}

stage5_watchdog_status_value() {
  local status="$1" field="$2"
  awk -F= -v wanted="${field}" '
    $1 == wanted {
      count++
      value=substr($0, length(wanted) + 2)
    }
    END {
      if (count != 1 || value == "" ||
          value ~ /[\t\r\n]/) {
        exit 1
      }
      print value
    }
  ' "${status}"
}

stage5_watchdog_safe_label() {
  printf '%s' "$1" | tr '/:,\t |' '______'
}

# Resolve one method owner to the exact run-owned watchdog status that covers
# its dataset/view row.  A method is successful only when every matching
# watchdog record is a terminal STATE=OK record; an absent, malformed, or
# failed record remains fail-closed.
stage5_watchdog_state_for_method() {
  local ds="$1" view="$2" method="$3"
  local status_dir="${ECODA_RUN_ROOT:-}/status/watchdogs"
  local matrix status status_label expected_label safe
  local matrix_base matrix_prefix matrix_suffix matrix_columns matrix_match state
  local row_ds row_view row_method row_extra
  local expected_count=0 saw_ok=0 saw_fail=0
  STAGE5_WATCHDOG_STATE=FAIL
  [[ -d "${status_dir}" && ! -L "${status_dir}" ]] || return 0
  matrix_prefix="${view}_${method}"
  for matrix in "${ECODA_RUN_ROOT:-}"/manifests/matrix_*.tsv; do
    [[ -e "${matrix}" ]] || continue
    [[ -f "${matrix}" && ! -L "${matrix}" && -r "${matrix}" ]] || continue
    ecoda_validate_run_owned_path "${matrix}" "${ECODA_RUN_ROOT}" \
      >/dev/null 2>&1 || continue
    matrix_columns=0
    if ecoda_validate_manifest "${matrix}" 3 >/dev/null 2>&1; then
      matrix_columns=3
    elif ecoda_validate_manifest "${matrix}" 4 >/dev/null 2>&1; then
      matrix_columns=4
    else
      continue
    fi
    matrix_match=0
    while IFS=$'\t' read -r row_ds row_view row_method row_extra; do
      if [[ "${row_ds}" == "${ds}" && "${row_view}" == "${view}" &&
            "${row_method}" == "${method}" ]]; then
        matrix_match=1
      fi
    done < "${matrix}"
    [[ ${matrix_match} -eq 1 ]] || continue
    matrix_base="${matrix##*/}"
    matrix_base="${matrix_base#matrix_}"
    matrix_base="${matrix_base%.tsv}"
    if [[ "${matrix_base}" == "${matrix_prefix}" ]]; then
      expected_label="${view}__${method}"
    elif [[ "${matrix_base}" == "${matrix_prefix}_"* ]]; then
      matrix_suffix="${matrix_base#${matrix_prefix}_}"
      case "${matrix_suffix}" in
        cpu|default_gpu|any_gpu) ;;
        *) saw_fail=1; continue ;;
      esac
      expected_label="${view}__${method}__${matrix_suffix}"
    else
      saw_fail=1
      continue
    fi
    expected_count=$((expected_count + 1))
    safe="$(stage5_watchdog_safe_label "${expected_label}")" || {
      saw_fail=1
      continue
    }
    status="${status_dir}/${safe}.status"
    [[ -s "${status}" && -f "${status}" && ! -L "${status}" &&
       -r "${status}" ]] || {
      saw_fail=1
      continue
    }
    ecoda_validate_run_owned_path "${status}" "${ECODA_RUN_ROOT}" \
      >/dev/null 2>&1 || {
      saw_fail=1
      continue
    }
    if ! status_label="$(stage5_watchdog_status_value "${status}" LABEL \
        2>/dev/null)"; then
      saw_fail=1
      continue
    fi
    [[ "${status_label}" == "${expected_label}" ]] || {
      saw_fail=1
      continue
    }
    if ! state="$(stage5_watchdog_status_value "${status}" STATE \
        2>/dev/null)"; then
      saw_fail=1
      continue
    fi
    case "${state}" in
      OK) saw_ok=$((saw_ok + 1)) ;;
      FAIL) saw_fail=1 ;;
      *) saw_fail=1 ;;
    esac
  done
  if [[ ${expected_count} -gt 0 && ${saw_ok} -eq ${expected_count} &&
        ${saw_fail} -eq 0 ]]; then
    STAGE5_WATCHDOG_STATE=OK
  fi
}

stage5_watchdog_state_for_owner_key() {
  local owner_key="$1"
  local scope ds view method extra
  STAGE5_WATCHDOG_STATE=FAIL
  IFS='/' read -r scope ds view method extra <<< "${owner_key}"
  case "${scope}" in
    ordinary|uncorrected|corrected) ;;
    *) return 0 ;;
  esac
  [[ -n "${ds}" && -n "${view}" && -n "${method}" &&
     -z "${extra}" ]] || return 0
  stage5_watchdog_state_for_method "${ds}" "${view}" "${method}"
}

# Global artifact owners are keyed only by their canonical path.  Resolve that
# path through the exact stage5 owner manifest and centralized artifact
# contract before consulting its method watchdog state; never infer success
# from a method name or artifact presence.
stage5_watchdog_state_for_artifact_path() {
  local target_path="$1"
  local owner_file="${ECODA_RUN_ROOT:-}/manifests/owners.tsv"
  local owner_key owner extra owner_metadata_key expected_owner
  local scope ds view method owner_path expected_path
  local found=0 saw_ok=0 saw_fail=0
  STAGE5_WATCHDOG_STATE=FAIL
  [[ -r "${owner_file}" && ! -L "${owner_file}" ]] || return 0
  ecoda_validate_run_owned_path "${owner_file}" "${ECODA_RUN_ROOT}" \
    >/dev/null 2>&1 || return 0
  while IFS=$'\t' read -r owner_key owner extra; do
    [[ -n "${owner_key}" && -n "${owner}" && -z "${extra}" ]] || continue
    expected_owner="$(ecoda_owner_dir stage5 "${owner_key}" 2>/dev/null ||
      true)"
    owner_metadata_key="$(ecoda_owner_field "${owner}" KEY 2>/dev/null ||
      true)"
    [[ "${owner}" == "${expected_owner}" &&
       "${owner_metadata_key}" == "${owner_key}" ]] || continue
    IFS='/' read -r scope ds view method extra <<< "${owner_key}"
    case "${scope}" in
      ordinary|uncorrected|corrected) ;;
      *) continue ;;
    esac
    [[ -n "${ds}" && -n "${view}" && -n "${method}" &&
       -z "${extra}" ]] || continue
    _ecoda_stage5_artifacts_for "${ds}" "${view}" "${method}" \
      >/dev/null 2>&1 || continue
    for owner_path in "${ECODA_BENCHMARK_ARTIFACTS[@]}" \
                      "${ECODA_BENCHMARK_ARTIFACT_NAS[@]}"; do
      [[ -n "${owner_path}" ]] || continue
      if ! expected_path="$(ecoda_canonical_path "${owner_path}" \
          2>/dev/null)"; then
        continue
      fi
      [[ "${expected_path}" == "${target_path}" ]] || continue
      found=1
      stage5_watchdog_state_for_method "${ds}" "${view}" "${method}"
      if [[ "${STAGE5_WATCHDOG_STATE}" == OK ]]; then
        saw_ok=1
      else
        saw_fail=1
      fi
    done
  done < "${owner_file}"
  if [[ ${found} -eq 1 && ${saw_ok} -eq 1 && ${saw_fail} -eq 0 ]]; then
    STAGE5_WATCHDOG_STATE=OK
  fi
}

stage5_finalize_tracked_owner_states() {
  local reason="$1"
  local owner owner_state owner_key owner_path owner_run owner_stage expected_owner
  local current_run="${RUN_ID:-${ECODA_RUN_ID:-${ECODA_RUN_ROOT##*/}}}"
  local owner_reason owner_mutable rc=0
  declare -p ECODA_ACQUIRED_OWNERS >/dev/null 2>&1 ||
    ECODA_ACQUIRED_OWNERS=()
  if [[ ${#ECODA_ACQUIRED_OWNERS[@]} -eq 0 ]]; then
    return 0
  fi
  for owner in "${ECODA_ACQUIRED_OWNERS[@]}"; do
    [[ -n "${owner}" ]] || { rc=1; continue; }
    owner_state=FAIL
    owner_mutable=0
    case "${owner}" in
      "${ECODA_OWNERS_ROOT:-}/artifact/"*)
        if _ecoda_artifact_owner_validate_dir "${owner}" \
            >/dev/null 2>&1; then
          owner_run="${ECODA_ARTIFACT_OWNER_RUN:-}"
          owner_stage="${ECODA_ARTIFACT_OWNER_STAGE:-}"
          owner_path="${ECODA_ARTIFACT_OWNER_CANONICAL_PATH:-}"
          if [[ -n "${owner_path}" && "${owner_run}" == "${current_run}" &&
                "${owner_stage}" == stage5 ]]; then
            owner_mutable=1
            stage5_watchdog_state_for_artifact_path "${owner_path}"
            owner_state="${STAGE5_WATCHDOG_STATE}"
          fi
        fi
        ;;
      *)
        owner_key="$(ecoda_owner_field "${owner}" KEY 2>/dev/null || true)"
        expected_owner="$(ecoda_owner_dir stage5 "${owner_key}" \
          2>/dev/null || true)"
        owner_run="$(ecoda_owner_field "${owner}" RUN_ID 2>/dev/null ||
          true)"
        owner_stage="$(ecoda_owner_field "${owner}" STAGE 2>/dev/null ||
          true)"
        if [[ -n "${owner_key}" && "${owner}" == "${expected_owner}" &&
              "${owner_run}" == "${current_run}" &&
              "${owner_stage}" == stage5 ]]; then
          owner_mutable=1
          stage5_watchdog_state_for_owner_key "${owner_key}"
          owner_state="${STAGE5_WATCHDOG_STATE}"
        fi
        ;;
    esac
    if [[ ${owner_mutable} -ne 1 ]]; then
      rc=1
      continue
    fi
    if [[ "${owner_state}" == OK ]]; then
      owner_reason="method watchdog terminal OK; ${reason}"
    else
      owner_reason="${reason}"
    fi
    ecoda_owner_set_state "${owner}" "${owner_state}" "${owner_reason}" ||
      rc=1
  done
  return "${rc}"
}
stage5_finalize_owner_manifest() {
  local state="$1" reason="$2" owner_file="${ECODA_RUN_ROOT:-}/manifests/owners.tsv"
  local owner_key owner extra owner_metadata_key expected_owner owner_run owner_stage owner_state
  local current_run="${RUN_ID:-${ECODA_RUN_ID:-${ECODA_RUN_ROOT##*/}}}"
  local owner_reason owner_mutable rc=0
  [[ "${state}" == "OK" || "${state}" == "FAIL" ]] || return 1
  [[ -r "${owner_file}" && ! -L "${owner_file}" ]] || return 1
  ecoda_validate_run_owned_path "${owner_file}" "${ECODA_RUN_ROOT}" \
    >/dev/null 2>&1 || return 1
  while IFS=$'\t' read -r owner_key owner extra; do
    if [[ -z "${owner_key}" || -z "${owner}" || -n "${extra}" ]]; then
      rc=1
      continue
    fi
    expected_owner="$(ecoda_owner_dir stage5 "${owner_key}" 2>/dev/null ||
      true)"
    owner_metadata_key="$(ecoda_owner_field "${owner}" KEY 2>/dev/null ||
      true)"
    owner_run="$(ecoda_owner_field "${owner}" RUN_ID 2>/dev/null || true)"
    owner_stage="$(ecoda_owner_field "${owner}" STAGE 2>/dev/null || true)"
    owner_state="${state}"
    owner_mutable=0
    if [[ "${state}" == "FAIL" ]]; then
      if [[ "${owner}" == "${expected_owner}" &&
            "${owner_metadata_key}" == "${owner_key}" &&
            "${owner_run}" == "${current_run}" &&
            "${owner_stage}" == stage5 ]]; then
        owner_mutable=1
        stage5_watchdog_state_for_owner_key "${owner_key}"
        owner_state="${STAGE5_WATCHDOG_STATE}"
      else
        rc=1
      fi
    elif [[ "${owner}" == "${expected_owner}" &&
            "${owner_metadata_key}" == "${owner_key}" &&
            "${owner_run}" == "${current_run}" &&
            "${owner_stage}" == stage5 ]]; then
      owner_mutable=1
    else
      rc=1
    fi
    if [[ "${state}" == "FAIL" && "${owner_state}" == OK ]]; then
      owner_reason="method watchdog terminal OK; ${reason}"
    else
      owner_reason="${reason}"
    fi
    if [[ ${owner_mutable} -eq 1 ]]; then
      ecoda_owner_set_state "${owner}" "${owner_state}" "${owner_reason}" ||
        rc=1
    fi
  done < "${owner_file}"
  return "${rc}"
}

stage5_abort() {
  local reason="$1"
  local rc=0
  stage5_finalize_tracked_owner_states "${reason}" || rc=1
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
     "${kind}" == "PREFLIGHT" || "${kind}" == "METADATA_EXPORT" ]] || return 1
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


gpu_method_is_default() {
  case ",${BENCHMARK_GPU_DEFAULT_METHODS}," in
    *,"$1",*) return 0 ;;
    *) return 1 ;;
  esac
}

method_spec() {
  local method="$1" effective_gpu_policy="${GPU_POLICY}"
  local gpu_request=()
  METHOD_RUNTIME_NV=0
  METHOD_TIME_LIMIT="${BENCHMARK_CPU_TIME_LIMIT}"
  METHOD_GPU_POLICY=cpu
  METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
  METHOD_THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
  METHOD_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}" --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
  METHOD_WORKER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
  case "${method}" in
    mrvi|scpoli)
      METHOD_RUNTIME_NV=1
      if [[ "${method}" == mrvi && -n "${PASS_ARG}" ]]; then
        effective_gpu_policy=cpu
      elif [[ "${effective_gpu_policy}" == auto ]]; then
        if ! gpu_method_is_default "${method}"; then
          effective_gpu_policy=any
        else
          effective_gpu_policy=default
        fi
      fi
      case "${effective_gpu_policy}" in
        cpu)
          METHOD_RUNTIME_NV=0
          METHOD_GPU_POLICY=cpu
          METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
          METHOD_THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
          METHOD_TIME_LIMIT="${BENCHMARK_CPU_TIME_LIMIT}"
          METHOD_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}" --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
          ;;
        default)
          METHOD_GPU_POLICY=default
          METHOD_PARTITION="${BENCHMARK_GPU_DEFAULT_PARTITION}"
          METHOD_THROTTLE="${BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE}"
          METHOD_TIME_LIMIT="${BENCHMARK_GPU_DEFAULT_TIME_LIMIT}"
          METHOD_FLAGS=(--gpus="${BENCHMARK_GPU_COUNT}" --constraint="${BENCHMARK_GPU_CONSTRAINT}" --cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}")
          ;;
        any)
          METHOD_GPU_POLICY=any
          METHOD_PARTITION="${BENCHMARK_GPU_ANY_PARTITION}"
          METHOD_THROTTLE="${BENCHMARK_GPU_ANY_ARRAY_THROTTLE}"
          METHOD_TIME_LIMIT="${BENCHMARK_GPU_ANY_TIME_LIMIT}"
          gpu_request=(--gpus="${BENCHMARK_GPU_COUNT}")
          if [[ -n "${BENCHMARK_GPU_ANY_VRAM_PER_GPU}" ]]; then
            gpu_request=(--gres="gpu:${BENCHMARK_GPU_COUNT},VramPerGpu:${BENCHMARK_GPU_ANY_VRAM_PER_GPU}")
          fi
          METHOD_FLAGS=("${gpu_request[@]}" --cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}")
          ;;
        *) echo "ERROR: unsupported effective GPU policy '${effective_gpu_policy}'." >&2; return 1 ;;
      esac
      METHOD_WORKER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
      ;;
    pilot|pilotgm|qot)
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
  if [[ "${method}" == mrvi && -n "${PASS_ARG}" && -n "${PARTITION_ARG}" ]]; then
    echo "ERROR: batch-view MRVI uses the deterministic CPU resource policy and rejects --partition overrides." >&2
    return 1
  fi
  if [[ -n "${PARTITION_ARG}" ]]; then
    METHOD_PARTITION="${PARTITION_ARG}"
    filtered=(); flag=""
    for flag in "${METHOD_FLAGS[@]}"; do
      case "${flag}" in --constraint|--constraint=*) ;; *) filtered+=("${flag}");; esac
    done
    METHOD_FLAGS=("${filtered[@]}")
  fi
}

# Heavy methods are split only where each parameter run is expensive enough to
# amortize loading one dataset into a fresh worker.  Small R methods remain
# one task per dataset; batch-effect runs are intentionally unsharded.
benchmark_method_shard_rows() {
  local method="$1" view="$2"
  [[ -z "${PASS_ARG}" && "${view}" == "benchmark_analysis" ]] || return 0
  case "${method}" in
    gloscope)
      printf '%s\n' \
        $'hvg2000_pcadims10\tcpu' \
        $'hvg2000_pcadims30\tcpu' \
        $'hvg2000_pcadims50\tcpu' \
        $'hvg1000_pcadims30\tcpu' \
        $'hvg3000_pcadims30\tcpu'
      ;;
    mrvi)
      printf '%s\n' \
        $'hvg2000\tdefault_gpu' \
        $'hvg1000\tcpu' \
        $'hvg3000\tcpu'
      ;;
    scpoli)
      printf '%s\n' \
        $'hvg2000_highres_dims15\tdefault_gpu' \
        $'hvg2000_lowres_dims15\tany_gpu' \
        $'hvg1000_highres_dims15\tany_gpu' \
        $'hvg2000_highres_dims2\tany_gpu' \
        $'hvg2000_highres_dims3\tany_gpu' \
        $'hvg2000_highres_dims5\tany_gpu' \
        $'hvg2000_highres_dims10\tany_gpu' \
        $'hvg3000_highres_dims15\tany_gpu'
      ;;
    pilotgm)
      # PILOT-GM-VAE is benchmark-default-only: no parameter screening.
      printf '%s\n' $'hvg2000_highres\tcpu'
      ;;
  esac
}

method_spec_for_group() {
  local method="$1" resource_class="$2" saved_gpu_policy="${GPU_POLICY}" rc
  case "${resource_class}" in
    base)
      method_spec "${method}"
      rc=$?
      ;;
    default_gpu)
      if [[ "${GPU_POLICY}" == "auto" ]]; then GPU_POLICY=default; fi
      method_spec "${method}"
      rc=$?
      ;;
    any_gpu)
      if [[ "${GPU_POLICY}" == "auto" ]]; then GPU_POLICY=any; fi
      method_spec "${method}"
      rc=$?
      ;;
    cpu)
      if [[ "${method}" == "mrvi" && "${GPU_POLICY}" == "auto" ]]; then
        METHOD_RUNTIME_NV=0
        METHOD_GPU_POLICY=cpu
        METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
        METHOD_THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
        METHOD_TIME_LIMIT="${BENCHMARK_CPU_TIME_LIMIT}"
        METHOD_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}" --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
        METHOD_WORKER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
        if [[ -n "${PARTITION_ARG}" ]]; then
          METHOD_PARTITION="${PARTITION_ARG}"
          METHOD_FLAGS=(--cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
        fi
        rc=0
      else
        method_spec "${method}"
        rc=$?
      fi
      ;;
    *)
      echo "ERROR: unsupported benchmark resource class '${resource_class}'." >&2
      rc=1
      ;;
  esac
  GPU_POLICY="${saved_gpu_policy}"
  return "${rc}"
}

stage5_validate_path_component() {
  local component="${1:-}" kind="${2:-path}"
  case "${component}" in
    ""|"."|".."|/*|*/*|*$'\n'*|*$'\r'*|*$'\t'*)
      echo "ERROR: invalid ${kind} path component." >&2
      return 1
      ;;
  esac
  [[ "${component}" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ ]] || {
    echo "ERROR: invalid ${kind} path component." >&2
    return 1
  }
}
stage5_method_is_forced() {
  local method="$1" target_method
  [[ ${FORCE_ARG} -eq 1 ]] && return 0
  [[ ${FORCE_TARGETED_ARG} -eq 1 ]] || return 1
  for target_method in "${TARGET_METHODS[@]}"; do
    [[ "${target_method}" == "${method}" ]] && return 0
  done
  return 1
}
# Corrected-final targeted manifests rebuild stale prepare caches and their
# dependent method rows instead of reusing terminal final artifacts.
stage5_corrected_final_recovery_mode() {
  [[ "${ANALYSIS_VARIANT:-${ANALYSIS_VARIANT_ARG:-}}" == corrected_final &&
     "${PASS_ARG:-}" == corrected &&
     ${TARGET_METHODS_SET} -eq 1 ]]
}




stage5_configured_output_name() {
  local ds="$1" view="$2"
  [[ -r "${DATASETS_JSON_FILE:-}" ]] || return 1
  jq -er --arg ds "${ds}" --arg view "${view}" '
    def nonempty_string:
      (type == "string" and length > 0);
    .[$ds].views[$view] as $spec
    | if ($spec | type) != "object" then empty
      elif ($spec.output_file_name | nonempty_string) then $spec.output_file_name
      elif (($spec.output_file_name == null
             or (($spec.output_file_name | type) == "string"
                 and ($spec.output_file_name | length) == 0))
            and ($spec.output_file | nonempty_string)) then $spec.output_file
      else empty
      end
    | select(. != "." and . != "..")
    | select(test("^[A-Za-z0-9_][A-Za-z0-9_.-]*$"))
  ' "${DATASETS_JSON_FILE}" || return 1
}

stage5_canonical_output_path() {
  local root="${1:-}" ds="${2:-}" view="${3:-}" name
  local output_root candidate root_real output_root_real candidate_real
  local expected_output_root
  STAGE5_OUTPUT_PATH=""
  STAGE5_OUTPUT_CANONICAL_PATH=""
  [[ "${root}" = /* &&
      "${root}" != *$'\n'* && "${root}" != *$'\r'* && "${root}" != *$'\t'* ]] ||
    return 1
  stage5_validate_path_component "${ds}" dataset || return 1
  stage5_validate_path_component "${view}" view || return 1
  name="$(stage5_configured_output_name "${ds}" "${view}")" || return 1
  stage5_validate_path_component "${name}" output || return 1
  output_root="${root%/}/${ds}/output"
  candidate="${output_root}/${name}"
  [[ ! -e "${root}" && ! -L "${root}" || -d "${root}" ]] || return 1
  [[ ! -e "${output_root}" && ! -L "${output_root}" || -d "${output_root}" ]] || return 1
  root_real="$(ecoda_canonical_path "${root}")" || return 1
  output_root_real="$(ecoda_canonical_path "${output_root}")" || return 1
  candidate_real="$(ecoda_canonical_path "${candidate}")" || return 1
  expected_output_root="${root_real%/}/${ds}/output"
  [[ "${output_root_real}" == "${expected_output_root}" ]] || return 1
  case "${candidate_real}" in
    "${output_root_real}"/*) ;;
    *) return 1 ;;
  esac
  STAGE5_OUTPUT_PATH="${candidate}"
  STAGE5_OUTPUT_CANONICAL_PATH="${candidate_real}"
  printf '%s' "${candidate}"
}

stage5_input_path() {
  stage5_canonical_output_path "${HPC_SCRATCH_DIR}" "$1" "$2"
}
stage5_validate_corrected_h5ad_contract() {
  local validator="${1:-}" path="${2:-}" view="${3:-}" identity_path="${4:-}"
  local -a validator_args=(
    --path "${path}"
    --view "${view}"
    --method "Stage 3 preprocessing"
    --expected-batch-contract "${identity_path}"
  )
  [[ -r "${validator}" && -s "${path}" && -s "${identity_path}" ]] || return 1
  if [[ "${ANALYSIS_VARIANT_ARG:-${ANALYSIS_VARIANT:-}}" == corrected_final ]]; then
    validator_args+=(--allow-missing-corrected-summary)
  fi
  "${PYTHON_BIN}" "${validator}" "${validator_args[@]}"
}


validate_input_row() {
  local ds="$1" view="$2" path
  path="$(stage5_input_path "${ds}" "${view}")" || return 1
  # Full persisted-content validation runs in the compute-node preflight
  # array.  The login submitter only checks presence here; source identity
  # creation below validates the sidecar, size, MD5, and Sample column.
  [[ -s "${path}" ]]
}

stage5_validate_corrected_source_contracts() {
  local ds view source_path identity_path row seen_sources=""
  local validator="${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py"
  [[ "${PASS_ARG:-}" == corrected ]] || return 0
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  [[ -r "${validator}" ]] || return 1
  while IFS=$'\t' read -r ds view _scope; do
    row="${ds}/${view}"
    case " ${seen_sources} " in
      *" ${row} "*) continue ;;
    esac
    seen_sources="${seen_sources} ${row}"
    source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
    identity_path="$(
      stage5_batch_contract_identity_path "${ds}" "${view}" preprocess
    )" || return 1
    stage5_validate_corrected_h5ad_contract \
      "${validator}" "${source_path}" "${view}" "${identity_path}" \
      >/dev/null 2>&1 || return 1

  done < "${MANIFEST}"
}


stage5_output_path_for_root() {
  stage5_canonical_output_path "$1" "$2" "$3"
}

stage5_require_input_ownerships() {
  local ds view source_path row seen_sources=""
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  while IFS=$'\t' read -r ds view _scope; do
    row="${ds}/${view}"
    case " ${seen_sources} " in *" ${row} "*) continue ;; esac
    seen_sources="${seen_sources} ${row}"
    source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
    [[ -s "${source_path}" ]] || return 1
    ecoda_require_input_ownership "${source_path}" "${RUN_ID}" || return 1
  done < "${MANIFEST}"
}

stage5_validate_input_provenance() {
  local ds view source_path row seen_sources="" owner_stage producer_run
  local producer owner_dir producer_ok
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]] || return 0
  command -v ecoda_stage5_artifact_owner_validate >/dev/null 2>&1 || return 1
  command -v ecoda_validate_input_artifact >/dev/null 2>&1 || return 1
  while IFS=$'\t' read -r ds view _scope; do
    row="${ds}/${view}"
    case " ${seen_sources} " in
      *" ${row} "*) continue ;;
    esac
    seen_sources="${seen_sources} ${row}"
    source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
    ECODA_ARTIFACT_OWNER_CANONICAL_PATH=""
    ecoda_artifact_owner_validate "${source_path}" >/dev/null 2>&1 || return 1
    owner_dir="$(ecoda_artifact_owner_dir \
      "${ECODA_ARTIFACT_OWNER_CANONICAL_PATH}")" || return 1
    owner_stage="${ECODA_ARTIFACT_OWNER_STAGE:-}"
    [[ -n "${owner_dir}" && "${ECODA_ARTIFACT_OWNER_STATE:-}" == "OK" ]] ||
      return 1
    case "${owner_stage}" in
      stage3)
        producer_run="${STAGE5_INPUT_PRODUCER_RUN_ID:-${ECODA_ARTIFACT_OWNER_RUN:-}}"
        producer_ok=0
        for producer in stage3 stage3_preflight; do
          if ecoda_validate_input_artifact \
              "${source_path}" "${producer}" "${producer_run}" >/dev/null 2>&1; then
            producer_ok=1
            break
          fi
        done
        ;;
      stage4)
        producer_run="${STAGE5_INPUT_PRODUCER_RUN_ID:-${ECODA_ARTIFACT_OWNER_RUN:-}}"
        producer_ok=0
        if ecoda_validate_input_artifact \
            "${source_path}" stage4_merge "${producer_run}" >/dev/null 2>&1; then
          producer_ok=1
        fi
        ;;
      *)
        return 1
        ;;
    esac
    ecoda_validate_run_id "${producer_run}" || return 1
    [[ ${producer_ok} -eq 1 ]] || return 1
  done < "${MANIFEST}"
}

stage5_validate_source_artifact_records() {
  local ds view source_path identity_path row seen_sources="" producer="stage5_preflight"
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  while IFS=$'\t' read -r ds view _scope; do
    row="${ds}/${view}"
    case " ${seen_sources} " in *" ${row} "*) continue ;; esac
    seen_sources="${seen_sources} ${row}"
    source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
    ecoda_validate_artifact_record "${source_path}" "${producer}" "${RUN_ID}" ||
      return 1
    if [[ "${PASS_ARG:-}" == corrected ]]; then
      identity_path="$(
        stage5_batch_contract_identity_path "${ds}" "${view}" preprocess
      )" || return 1
      stage5_validate_corrected_h5ad_contract \
        "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
        "${source_path}" "${view}" "${identity_path}" || return 1

    else
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
        --path "${source_path}" --view "${view}" \
        --method "Stage 5 source preflight" || return 1
    fi
  done < "${MANIFEST}"
}


stage5_repair_missing_h5ad_sidecars() {
  local repair_manifest="${ECODA_RUN_ROOT}/manifests/h5ad_checksum_repairs.tsv"
  local repair_tmp="${repair_manifest}.build.$$"
  local ds view path sidecar digest size root repairs=0
  local scratch_path nas_path scratch_digest scratch_size source_key seen_sources=""
  local identity_path
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0
  : > "${repair_tmp}" || return 1
  while IFS=$'\t' read -r ds view _scope; do
    source_key="${ds}/${view}"
    case " ${seen_sources} " in
      *" ${source_key} "*) continue ;;
    esac
    seen_sources="${seen_sources} ${source_key}"
    scratch_path="$(stage5_input_path "${ds}" "${view}")" || {
      rm -f "${repair_tmp}"
      return 1
    }
    nas_path="$(stage5_output_path_for_root "${NAS_TARGET_DIR}" "${ds}" "${view}")" || {
      rm -f "${repair_tmp}"
      return 1
    }
    [[ -s "${scratch_path}" && -s "${nas_path}" ]] || {
      rm -f "${repair_tmp}"
      return 1
    }
    ecoda_require_input_ownership "${scratch_path}" "${RUN_ID}" || return 1
    ecoda_require_input_ownership "${nas_path}" "${RUN_ID}" || return 1
    echo "Stage 5 source checksum validation: ${ds}/${view}" >&2
    for root in "${HPC_SCRATCH_DIR}" "${NAS_TARGET_DIR}"; do
      if [[ "${root}" == "${HPC_SCRATCH_DIR}" ]]; then
        path="${scratch_path}"
      else
        path="${nas_path}"
      fi
      sidecar="${path}.md5"
      if [[ -e "${sidecar}" || -L "${sidecar}" ]]; then
        # Validate the existing sidecar once and reuse its strict digest/size
        # record for the scratch/NAS identity comparison below.
        ecoda_validate_checksum "${path}" || {
          rm -f "${repair_tmp}"
          return 1
        }
        digest="${ECODA_CHECKSUM_MD5}"
        size="${ECODA_CHECKSUM_SIZE}"
      else
        # A missing sidecar is repaired only after the persisted H5AD content
        # contract passes. Existing invalid sidecars are never overwritten.
        if [[ "${PASS_ARG:-}" == corrected ]]; then
          identity_path="$(
            stage5_batch_contract_identity_path "${ds}" "${view}" preprocess
          )" || {
            rm -f "${repair_tmp}"
            return 1
          }
          stage5_validate_corrected_h5ad_contract \
            "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
            "${path}" "${view}" "${identity_path}" >/dev/null 2>&1 || {
            rm -f "${repair_tmp}"
            return 1
          }

        else
          "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
            --path "${path}" --view "${view}" \
            --method "Stage 5 source checksum repair" >/dev/null 2>&1 || {
            rm -f "${repair_tmp}"
            return 1
          }
        fi
        ecoda_write_checksum "${path}" || {
          rm -f "${repair_tmp}"
          return 1
        }
        digest="${ECODA_CHECKSUM_MD5}"
        size="${ECODA_CHECKSUM_SIZE}"
        ecoda_validate_checksum_record "${path}" "${digest}" "${size}" || {
          rm -f "${repair_tmp}"
          return 1
        }
        printf '%s\t%s\t%s\t%s\t%s\n' "${root}" "${ds}" "${view}" \
          "${path}" "${digest}:${size}" >> "${repair_tmp}" || {
          rm -f "${repair_tmp}"
          return 1
        }
        repairs=$((repairs + 1))
      fi
      if [[ "${root}" == "${HPC_SCRATCH_DIR}" ]]; then
        scratch_digest="${digest}"
        scratch_size="${size}"
      else
        [[ "${scratch_digest}" == "${digest}" &&
           "${scratch_size}" == "${size}" ]] || {
          rm -f "${repair_tmp}"
          return 1
        }
      fi
    done
    echo "Stage 5 source checksum validation passed: ${ds}/${view}" >&2
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
  local identity="${SOURCE_IDENTITY}" identity_sidecar="${SOURCE_IDENTITY}.md5"
  local identity_script
  if [[ -e "${identity}" || -L "${identity}" ]]; then
    ecoda_validate_run_owned_path "${identity}" "${ECODA_RUN_ROOT}" || return 1
  fi
  if [[ -e "${identity_sidecar}" || -L "${identity_sidecar}" ]]; then
    ecoda_validate_run_owned_path "${identity_sidecar}" "${ECODA_RUN_ROOT}" || return 1
  fi
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0
  identity_script="$(stage5_source_script src/utils/py/h5ad_source_identity.py)"
  if [[ -s "${identity}" && -s "${identity_sidecar}" ]]; then
    "${PYTHON_BIN}" "${identity_script}" \
      --identity "${identity}" --selection "${MANIFEST}" \
      --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" \
      --validated-source-sidecars ||
      return 1
  else
    "${PYTHON_BIN}" "${identity_script}" \
      --output "${identity}" --selection "${MANIFEST}" \
      --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" \
      --validated-source-sidecars ||
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
  local ds view source_path safe status state status_run status_dataset status_view status_task
  local preflight_id preflight_rc count=0 preflight_expected
  local identity_path
  local preflight_runtime_export
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0

  if [[ -n "${SYNC_ONLY_RUN}" && ! -s "${preflight_manifest}" ]]; then
    # Runs created before the compute boundary have no preflight record.  Do
    # not resubmit during recovery; validate the immutable selected sources
    # locally and fail closed if any source is invalid.
    while IFS=$'\t' read -r ds view _scope; do
      source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
      if [[ "${PASS_ARG:-}" == corrected ]]; then
        identity_path="$(
          stage5_batch_contract_identity_path "${ds}" "${view}" preprocess
        )" || return 1
        stage5_validate_corrected_h5ad_contract \
          "$(stage5_source_script src/utils/py/benchmark_h5ad_contract.py)" \
          "${source_path}" "${view}" "${identity_path}" \
          >/dev/null 2>&1 || return 1

      else
        "${PYTHON_BIN}" "$(stage5_source_script src/utils/py/benchmark_h5ad_contract.py)" \
          --path "${source_path}" --view "${view}" --method "Stage 5 sync recovery" \
          >/dev/null 2>&1 || return 1
      fi
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
      source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
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
      source_path="$(stage5_input_path "${ds}" "${view}")" || return 1
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
    preflight_worker="$(stage5_source_script src/utils/bash/h5ad_preflight_worker.sh)"
    stage5_validate_bound_runtime || return 1
    stage5_require_source_script "${preflight_worker}" || return 1
    preflight_runtime_export="${RUNTIME_EXPORT}"
    if [[ "${PASS_ARG:-}" == corrected ]]; then
      preflight_runtime_export="${preflight_runtime_export},ECODA_BATCH_CONTRACT_MANIFEST=${ECODA_BATCH_CONTRACT_MANIFEST}"
      if [[ "${ANALYSIS_VARIANT_ARG:-${ANALYSIS_VARIANT:-}}" == corrected_final ]]; then
        preflight_runtime_export="${preflight_runtime_export},H5AD_ALLOW_MISSING_SUMMARY=1"
      else
        preflight_runtime_export="${preflight_runtime_export},H5AD_ALLOW_MISSING_SUMMARY=0"
      fi
    fi
    set +e
    preflight_id="$(
      ecoda_submit_h5ad_preflight "${preflight_manifest}" "${status_dir}" \
        "${ECODA_RUN_ROOT}" require "${PARTITION_ARG:-${SLURM_PARTITION_BENCHMARK_CPU}}" \
        "${MEMORY}" "${THROTTLE}" "${preflight_logs}" stage5 \
        "${preflight_worker}" \
        "${preflight_runtime_export}"
    )"
    preflight_rc=$?
    set -e
    if [[ "${preflight_id}" =~ ^[0-9]+$ ]]; then
      stage5_record_scheduler PREFLIGHT "${preflight_id}" || return 1
    fi
    if [[ ${preflight_rc} -ne 0 || ! "${preflight_id}" =~ ^[0-9]+$ ]]; then
      echo "ERROR: Stage 5 H5AD preflight scheduler wait failed: job=${preflight_id:-unknown} rc=${preflight_rc}" >&2
      return 1
    fi
    preflight_expected="$(awk 'END { print NR }' "${preflight_manifest}")" || return 1
    [[ "${preflight_expected}" =~ ^[1-9][0-9]*$ ]] || {
      echo "ERROR: Stage 5 H5AD preflight manifest has no rows" >&2
      return 1
    }
    if ! ecoda_wait_array_accounting "${preflight_id}" "${preflight_expected}" \
        "${H5AD_PREFLIGHT_ACCOUNTING_POLL_SECONDS:-30}"; then
      echo "ERROR: Stage 5 H5AD preflight array did not settle: job=${preflight_id}" >&2
      return 1
    fi
  ecoda_wait_h5ad_preflight_status_files "${preflight_manifest}" "${status_dir}" || {
    echo "ERROR: Stage 5 H5AD preflight statuses did not settle within ${H5AD_PREFLIGHT_STATUS_GRACE_SECONDS:-60}s" >&2
    return 1
  }
  while IFS=$'\t' read -r ds view _scope; do
    safe="$(_ecoda_safe_component "${ds}__${view}")"
    status="${status_dir}/${safe}.status"
    [[ -s "${status}" ]] || {
      echo "ERROR: Stage 5 H5AD preflight status is missing: ${status}" >&2
      return 1
    }
    status_run="$(sed -n 's/^RUN_ID=//p' "${status}" | head -1)"
    status_dataset="$(sed -n 's/^DATASET=//p' "${status}" | head -1)"
    status_view="$(sed -n 's/^VIEW=//p' "${status}" | head -1)"
    status_task="$(sed -n 's/^TASK_ID=//p' "${status}" | head -1)"
    state="$(sed -n 's/^STATE=//p' "${status}" | head -1)"
    [[ "${state}" == OK &&
       "${status_run}" == "${ECODA_RUN_ID}" &&
       "${status_dataset}" == "${ds}" &&
       "${status_view}" == "${view}" &&
       "${status_task}" == "$((count + 1))" ]] || {
      echo "ERROR: Stage 5 H5AD preflight status mismatch: ${status}" >&2
      return 1
    }
    count=$((count + 1))
  done < "${preflight_manifest}"
  [[ ${count} -gt 0 ]] || return 1
  fi
}

# Resolve the immutable dataset/view selection before any scheduler submission.
RUN_ID="${ECODA_RUN_ID:-$(ecoda_new_run_id stage5)}"
if [[ -n "${SYNC_ONLY_RUN}" ]]; then
  ecoda_open_run "${SYNC_ONLY_RUN}" || exit 1
  RUN_ID="${SYNC_ONLY_RUN}"
  export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT
  stage5_bind_run_identity ||
    { echo "legacy_source_unpinned" >&2; exit 1; }
  stage5_validate_bound_runtime || exit 1
  RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)" || {
    echo "ERROR: Stage 5 bound runtime export construction failed." >&2
    exit 1
  }
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  ecoda_validate_run_owned_path "${MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage5_abort "Stage 5 selection manifest is not run-owned"
  ecoda_validate_manifest "${MANIFEST}" 3 ||
    stage5_abort "Stage 5 selection manifest is invalid"
  ecoda_validate_checksum "${MANIFEST}" ||
    stage5_abort "Stage 5 selection checksum is invalid"
  methods_meta="$(sed -n 's/^METHODS=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  target_methods_meta="$(sed -n 's/^TARGET_METHODS=//p' \
    "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  if [[ -z "${target_methods_meta}" && -z "${METHODS_ARG}" ]]; then
    METHODS_ARG="${methods_meta}"
    METHODS_SET=1
  fi
  pass_meta="$(sed -n 's/^PASS=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  [[ -n "${PASS_ARG}" ]] || PASS_ARG="${pass_meta}"
  [[ -n "${PASS_ARG}" ]] && unset BENCHMARK_MANIFEST
  variant_meta="$(sed -n 's/^ANALYSIS_VARIANT=//p' \
    "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  if [[ -z "${ANALYSIS_VARIANT_ARG}" ]]; then
    ANALYSIS_VARIANT_ARG="${variant_meta}"
  fi
  if [[ -n "${variant_meta}" && "${ANALYSIS_VARIANT_ARG}" != "${variant_meta}" ]]; then
    stage5_abort "sync-only analysis variant does not match run metadata"
  fi
  if [[ -n "${ANALYSIS_VARIANT_ARG}" ]]; then
    ANALYSIS_VARIANT_SET=1
    PASS_SET=1
    if [[ -n "${target_methods_meta}" ]]; then
      [[ ${METHODS_SET} -eq 0 ]] || {
        stage5_abort "final sync-only metadata contains target methods and a full METHODS selection"
      }
      if [[ ${TARGET_METHODS_SET} -eq 1 ]]; then
        [[ "${TARGET_METHODS_ARG}" == "${target_methods_meta}" ]] || {
          stage5_abort "final sync-only target methods do not match run metadata"
        }
      else
        TARGET_METHODS_ARG="${target_methods_meta}"
        TARGET_METHODS_SET=1
      fi
      METHODS_ARG=""
      METHODS_SET=0
      ecoda_split_csv "${TARGET_METHODS_ARG}" ||
        stage5_abort "invalid stored target methods"
      ecoda_assert_unique_items "${ECODA_ARRAY[@]}" ||
        stage5_abort "duplicate stored target methods"
      TARGET_METHODS=("${ECODA_ARRAY[@]}")
      for stored_target_method in "${TARGET_METHODS[@]}"; do
        case ",${EXPECTED_BATCH_METHODS}," in
          *,"${stored_target_method}",*) ;;
          *) stage5_abort "final sync-only metadata contains an unsupported target method" ;;
        esac
      done
    elif [[ -z "${METHODS_ARG}" ]]; then
      METHODS_ARG="${methods_meta}"
      METHODS_SET=1
    fi
    if [[ -z "${target_methods_meta}" && ${TARGET_METHODS_SET} -eq 1 ]]; then
      stage5_abort "final sync-only target methods are missing from run metadata"
    fi
    stage5_validate_final_selection || stage5_abort "stored final selection is invalid"
    stage5_configure_analysis_context ||
      stage5_abort "stored Stage 5 analysis context is invalid"
  fi
else
  ecoda_init_run stage5 "${RUN_ID}" >/dev/null
  export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT
  stage5_install_source_manifest ||
    stage5_abort "new Stage 5 run has no immutable source manifest"
  runtime_submission_mode="${ECODA_RUNTIME_MODE:-host}"
  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    runtime_submission_mode=apptainer
    export ECODA_RUNTIME_MODE=apptainer
  fi
  ecoda_runtime_validate_submission "${runtime_submission_mode}" || {
    stage5_abort "Stage 5 immutable runtime validation failed"
  }
  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]]; then
    [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]] ||
      stage5_abort "Stage 5 runtime submission did not write runtime.identity"
    stage5_bind_run_identity ||
      stage5_abort "Stage 5 run-bound source/runtime identity is invalid"
  fi
  RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)" || {
    stage5_abort "Stage 5 runtime export construction failed"
  }
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  MANIFEST_TMP="${MANIFEST}.build.$$"
  : > "${MANIFEST_TMP}"
  if [[ -n "${SELECTION_FILE_ARG}" ]]; then
    [[ -r "${SELECTION_FILE_ARG}" ]] || stage5_abort "selection file is unreadable"
    ecoda_validate_manifest "${SELECTION_FILE_ARG}" 3 ||
      stage5_abort "selection file must contain exactly three columns per row"
    while IFS=$'\t' read -r ds input_view row_label; do
      [[ -n "${ds}" && -n "${input_view}" && -n "${row_label}" ]] ||
        stage5_abort "selection rows require DATASET<TAB>VIEW<TAB>LABEL"
      stage5_validate_path_component "${ds}" dataset ||
        stage5_abort "invalid Stage 5 dataset selection"
      stage5_validate_path_component "${input_view}" view ||
        stage5_abort "invalid Stage 5 view selection"
      stage5_validate_path_component "${row_label}" selection ||
        stage5_abort "invalid Stage 5 selection label"
      ecoda_dataset_exists "${ds}" ||
        stage5_abort "unknown dataset ${ds}"
      if [[ -n "${PASS_ARG}" ]]; then
        [[ "${input_view}" == "batch_effect_${PASS_ARG}" ]] ||
          stage5_abort "pass-mode selection has the wrong view: ${input_view}"
        view="${input_view}"
      else
        view="${input_view}"
        case "${view}" in trans|zeroimp) view="benchmark_analysis" ;; esac
      fi
      ecoda_view_exists "${ds}" "${view}" ||
        stage5_abort "${ds}/${view} is not declared"
      stage5_input_path "${ds}" "${view}" >/dev/null ||
        stage5_abort "invalid Stage 5 H5AD path contract for ${ds}/${view}"
      [[ -n "$(ecoda_view_input_name "${ds}" "${view}")" &&
        -n "$(ecoda_view_output_name "${ds}" "${view}")" ]] ||
        stage5_abort "${ds}/${view} has no input/output"
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
      stage5_validate_path_component "${ds}" dataset ||
        stage5_abort "invalid Stage 5 dataset selection"
      ecoda_dataset_exists "${ds}" || stage5_abort "unknown dataset ${ds}"
      if [[ -n "${PASS_ARG}" ]]; then
        view="batch_effect_${PASS_ARG}"
        ecoda_view_exists "${ds}" "${view}" ||
          stage5_abort "${ds}/${view} is not declared"
      else
        view="benchmark_analysis"
      fi
      stage5_validate_path_component "${view}" view ||
        stage5_abort "invalid Stage 5 view selection"
      stage5_input_path "${ds}" "${view}" >/dev/null ||
        stage5_abort "invalid Stage 5 H5AD path contract for ${ds}/${view}"
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
stage5_bind_method_matrix ||
  stage5_abort "failed to bind the run-owned method matrix"
SOURCE_IDENTITY="${ECODA_RUN_ROOT}/manifests/source_identity.json"
SCHEDULER_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
if [[ -z "${SYNC_ONLY_RUN}" ]]; then
  ecoda_atomic_write "${SCHEDULER_FILE}" "" ||
    stage5_abort "failed to initialize Stage 5 scheduler manifest"
fi
if [[ -n "${PASS_ARG}" && -n "${METHODS_ARG}" &&
      ${TARGET_METHODS_SET} -eq 0 ]]; then
  [[ "${METHODS_ARG}" == "${EXPECTED_BATCH_METHODS}" ]] ||
    stage5_abort "batch-effect pass requires the fixed ordered method suite"
fi

# Validate row syntax, owners, and source h5ads. No biological label is used as
# a processing covariate; LABEL is only a scheduler/output grouping token.
ecoda_validate_manifest "${MANIFEST}" 3 ||
  stage5_abort "invalid Stage 5 selection"
stage5_validate_input_provenance ||
  stage5_abort "Stage 5 source H5AD lacks a terminal upstream owner/record"
DATASET_NAMES=(); SEEN_DS=""; SEEN_ROW=""
while IFS=$'\t' read -r ds view row_label; do
  stage5_validate_path_component "${ds}" dataset ||
    stage5_abort "invalid Stage 5 dataset selection"
  stage5_validate_path_component "${view}" view ||
    stage5_abort "invalid Stage 5 view selection"
  stage5_validate_path_component "${row_label}" selection ||
    stage5_abort "invalid Stage 5 selection label"
  ecoda_dataset_exists "${ds}" || stage5_abort "unknown dataset ${ds}"
  ecoda_view_exists "${ds}" "${view}" ||
    stage5_abort "${ds}/${view} is not declared"
  stage5_input_path "${ds}" "${view}" >/dev/null ||
    stage5_abort "invalid Stage 5 H5AD path contract for ${ds}/${view}"
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
stage5_require_input_ownerships ||
  stage5_abort "Stage 5 source input has an ACTIVE writer"
METHODS=("${BASELINE_METHODS[@]}")
ANALYSES=(_ecoda_none_)
ANALYSES_SELECTED=0
EXACT_SELECTION=0
[[ ${EXACT_BATCH_SELECTION} -eq 1 ]] && EXACT_SELECTION=1
BATCH_EFFECT_METHODS=(prepare_pseudobulk pseudobulk gloscope composition mrvi pilot qot)
if [[ ${TARGET_METHODS_SET} -eq 1 ]]; then
  METHODS=("${TARGET_METHODS[@]}")
elif [[ -n "${PASS_ARG}" ]]; then
  METHODS=("${BATCH_EFFECT_METHODS[@]}")
fi
if [[ -n "${METHODS_ARG}" ]]; then
  ecoda_split_csv "${METHODS_ARG}" ||
    stage5_abort "invalid benchmark method selection"
  METHODS=("${ECODA_ARRAY[@]}")
fi
if [[ -n "${PASS_ARG}" && ${TARGET_METHODS_SET} -eq 0 ]]; then
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
R_ENV_PREFLIGHT_REQUIRED=0
for method in "${METHODS[@]}" "${ANALYSES[@]}"; do
  case "${method}" in
    gloscope|mofa|pseudobulk|composition|scitd|prepare_pseudobulk|trans|zeroimp)
      R_ENV_PREFLIGHT_REQUIRED=1
      ;;
  esac
done
export ECODA_SELECTION_MANIFEST="${MANIFEST}"
export ECODA_EXACT_SELECTION="${EXACT_SELECTION}"

if [[ "${PASS_ARG}" == corrected ]]; then
  stage5_create_batch_contract_manifest ||
    stage5_abort "failed to create or validate corrected Stage 5 batch-contract manifest"
fi
stage5_validate_corrected_source_contracts ||
  stage5_abort "Stage 5 corrected source batch-contract validation failed"
DISPATCH_MANIFEST="${ECODA_RUN_ROOT}/manifests/dispatch_selection.tsv"
DISPATCH_MANIFEST_MD5=""
DISPATCH_MANIFEST_SIZE=""
stage5_build_dispatch_manifest() {
  local tmp="${DISPATCH_MANIFEST}.build.$$"
  local ds view row_label method extra key seen=""
  [[ -n "${DISPATCH_MANIFEST}" ]] || return 1
  if [[ -n "${SYNC_ONLY_RUN}" ]]; then
    [[ -f "${DISPATCH_MANIFEST}" && ! -L "${DISPATCH_MANIFEST}" &&
       -r "${DISPATCH_MANIFEST}" ]] || return 1
  else
    : > "${tmp}" || return 1
    if [[ ${METHOD_MATRIX_MODE} -eq 1 ]]; then
      while IFS=$'\t' read -r ds view method extra; do
        [[ -n "${ds}" && -n "${view}" && -n "${method}" && -z "${extra}" ]] ||
          return 1
        key="${ds}/${view}/${method}"
        case " ${seen} " in *" ${key} "*) continue ;; esac
        seen="${seen} ${key}"
        printf '%s\t%s\t%s\n' "${ds}" "${view}" "${method}" >> "${tmp}" || return 1
      done < "${METHOD_MATRIX}"
    else
      for method in "${METHODS[@]}" "${ANALYSES[@]}"; do
        [[ "${method}" == _ecoda_none_ ]] && continue
        while IFS=$'\t' read -r ds view row_label; do
          if [[ ${EXACT_SELECTION} -eq 1 && -z "${PASS_ARG}" ]]; then
            case "${method}:${row_label}" in
              prepare_pseudobulk:mofa|prepare_pseudobulk:pseudobulk|prepare_pseudobulk:composition|prepare_pseudobulk:prepare_pseudobulk) ;;
              prepare_pseudobulk:*) continue ;;
              *:"${method}") ;;
              *) continue ;;
            esac
          fi
          key="${ds}/${view}/${method}"
          case " ${seen} " in *" ${key} "*) continue ;; esac
          seen="${seen} ${key}"
          printf '%s\t%s\t%s\n' "${ds}" "${view}" "${method}" >> "${tmp}" || return 1
        done < "${MANIFEST}"
      done
    fi
    [[ -s "${tmp}" ]] || {
      rm -f "${tmp}"
      return 1
    }
    ecoda_atomic_install_manifest "${tmp}" "${DISPATCH_MANIFEST}" 3 || {
      rm -f "${tmp}"
      return 1
    }
    rm -f "${tmp}"
    ecoda_write_checksum "${DISPATCH_MANIFEST}" || return 1
  fi
  ecoda_validate_run_owned_path "${DISPATCH_MANIFEST}" "${ECODA_RUN_ROOT}" || return 1
  ecoda_validate_manifest "${DISPATCH_MANIFEST}" 3 || return 1
  ecoda_validate_checksum "${DISPATCH_MANIFEST}" || return 1
  DISPATCH_MANIFEST_MD5="${ECODA_CHECKSUM_MD5}"
  DISPATCH_MANIFEST_SIZE="${ECODA_CHECKSUM_SIZE}"
  export ECODA_DISPATCH_SELECTION="${DISPATCH_MANIFEST}"
}
stage5_build_dispatch_manifest ||
  stage5_abort "failed to build or validate the Stage 5 dispatch manifest"





benchmark_artifacts_for() {
  local ds="$1" view="$2" label="$3"
  ARTIFACT_PATHS=()
  _ecoda_stage5_artifacts_for "${ds}" "${view}" "${label}" || return 1
  ARTIFACT_PATHS=("${ECODA_BENCHMARK_ARTIFACTS[@]}")
  [[ ${#ARTIFACT_PATHS[@]} -gt 0 ]]
}


RDS_PREFLIGHT_DONE=""
RDS_PREFLIGHT_FAILED=""

benchmark_rds_group_valid() {
  local ds="$1" view="$2" pass="${PASS_ARG:-${ANALYSIS_PASS:-}}"
  local key="${pass:-ordinary}/${ds}/${view}"
  local safe list list_tmp label path metadata group_rc
  local selected_labels="" seen_label=""
  local has_rds=0
  ecoda_stage5_validate_identity "${pass}" || return 1
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 2
  [[ "${pass:-}" == corrected ]] && return 2
  case " ${RDS_PREFLIGHT_DONE} " in *" ${key} "*) return 0 ;; esac
  case " ${RDS_PREFLIGHT_FAILED} " in *" ${key} "*) return 1 ;; esac
  safe="$(_ecoda_safe_component "${key}")"
  list="${ECODA_RUN_ROOT}/manifests/rds_preflight_${safe}.tsv"
  list_tmp="${list}.build.$$"
  : > "${list_tmp}" || return 1
  for label in "${METHODS[@]}" "${ANALYSES[@]}"; do
    [[ "${label}" == _ecoda_none_ ]] && continue
    if [[ ${METHOD_MATRIX_MODE} -eq 1 ]] &&
       ! ecoda_stage5_method_matrix_allows "${ds}" "${view}" "${label}"; then
      continue
    fi
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
    --input-root "${HPC_SCRATCH_DIR}")
  [[ -s "${SOURCE_IDENTITY}" ]] &&
    rds_args+=(--source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
  [[ -n "${pass}" ]] && rds_args+=(--batch-pass "${pass}")
  if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
    rds_args+=(--analysis-variant "${ANALYSIS_VARIANT}")
  fi
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

stage5_validate_reusable_artifact() {
  local path="$1" producer="$2" record="" recorded_producer
  local owner_run owner_stage owner_state
  if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
    # Final artifacts are reusable across runs only after the immutable path,
    # terminal global owner, and producer-run record all agree.  A record in
    # the current recovery run is not required; the owner identifies the
    # prior terminal producer run.
    ecoda_stage5_validate_identity "${PASS_ARG:-${ANALYSIS_PASS:-}}" || return 1
    ecoda_stage5_validate_artifact_path "${path}" || return 1
    ecoda_require_input_ownership "${path}" "${RUN_ID}" || return 1
    ecoda_stage5_artifact_owner_validate "${path}" >/dev/null 2>&1 || return 1
    owner_state="${ECODA_ARTIFACT_OWNER_STATE:-}"
    owner_stage="${ECODA_ARTIFACT_OWNER_STAGE:-}"
    [[ "${owner_state}" == "OK" && "${owner_stage}" == "stage5" ]] || return 1
    owner_run="${ECODA_ARTIFACT_OWNER_RUN:-}"
    ecoda_validate_run_id "${owner_run}" || return 1
    record="$(ecoda_artifact_record_path "${path}" "${owner_run}" 2>/dev/null || true)"
    [[ -n "${record}" && -f "${record}" && ! -L "${record}" ]] || return 1
    recorded_producer="$(stage5_recorded_producer "${record}")" || return 1
    stage5_producer_allowed "${path}" "${producer}" "${recorded_producer}" ||
      return 1
    ecoda_validate_artifact_record \
      "${path}" "${recorded_producer}" "${owner_run}" || return 1
    return 0
  fi
  # A checksum/record is not enough to declare a reusable artifact: a live
  # writer from another run may be mutating the same path.  Missing owners
  # remain valid for legacy baseline artifacts.
  ecoda_require_input_ownership "${path}" "${RUN_ID}" || return 1
  if command -v ecoda_artifact_record_path >/dev/null 2>&1; then
    record="$(ecoda_artifact_record_path "${path}" "${RUN_ID}" 2>/dev/null || true)"
    if [[ -n "${record}" && -e "${record}" ]]; then
      recorded_producer="$(stage5_recorded_producer "${record}")" || return 1
      stage5_producer_allowed "${path}" "${producer}" "${recorded_producer}" ||
        return 1
      ecoda_validate_artifact_record \
        "${path}" "${recorded_producer}" "${RUN_ID}" || return 1
      return 0
    fi
  fi
  # Run-bound artifacts may be reused only with a run-owned record.  The
  # sidecar-only fallback is reserved for validator-only legacy mode.
  if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" &&
        "${BENCHMARK_MATRIX_TEST:-0}" != "1" ]]; then
    return 1
  fi
  ecoda_validate_checksum "${path}"
}
stage5_prepare_pseudobulk_valid() {
  local ds="$1" view="$2" path owner_dir producer_run record recorded_producer expected_cache
  local prepare_rds_args=()
  [[ -n "${PASS_ARG}" ]] || return 1
  stage5_corrected_final_recovery_mode && return 1
  benchmark_artifacts_for "${ds}" "${view}" prepare_pseudobulk || return 1
  [[ ${#ARTIFACT_PATHS[@]} -eq 1 ]] || return 1
  path="${ARTIFACT_PATHS[0]}"
  expected_cache="$(ecoda_stage5_batch_stem "${ds}" "${PASS_ARG}" "${ANALYSIS_VARIANT:-}")" ||
    return 1
  [[ "${path}" == "${ANALYSIS_ROOT}/pseudobulks/${expected_cache}_pseudobulk_hvg2000.rds" ]] ||
    return 1

  # Targeted recovery runs have a new RUN_ID.  Validate the immutable cache
  # against its terminal global owner first, then resolve the record under
  # that producer run rather than the recovery run.
  ecoda_validate_checksum "${path}" || return 1
  ecoda_stage5_artifact_owner_validate "${path}" >/dev/null 2>&1 || return 1
  owner_dir="${ECODA_ARTIFACT_OWNER_DIR:-}"
  [[ -n "${owner_dir}" &&
     "${ECODA_ARTIFACT_OWNER_STATE:-}" == "OK" &&
     "${ECODA_ARTIFACT_OWNER_STAGE:-}" == "stage5" ]] || return 1
  producer_run="${ECODA_ARTIFACT_OWNER_RUN:-}"
  ecoda_validate_run_id "${producer_run}" || return 1
  record="$(ecoda_artifact_record_path "${path}" "${producer_run}" 2>/dev/null || true)"
  [[ -n "${record}" && -f "${record}" && ! -L "${record}" ]] || return 1
  recorded_producer="$(stage5_recorded_producer "${record}")" || return 1
  stage5_producer_allowed "${path}" prepare_pseudobulk "${recorded_producer}" ||
    return 1
  ecoda_validate_artifact_record \
    "${path}" "${recorded_producer}" "${producer_run}" >/dev/null || return 1

  # The matrix fixture does not have a real R runtime, but it must still pass
  # every checksum, owner, and producer-record check above.
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  prepare_rds_args=(--artifact "${path}" --method prepare_pseudobulk
    --dataset "${ds}" --view "${view}" --input-root "${HPC_SCRATCH_DIR}"
    --config "${DATASETS_JSON_FILE}")
  prepare_rds_args+=(--batch-pass "${PASS_ARG}")
  if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
    prepare_rds_args+=(--analysis-variant "${ANALYSIS_VARIANT}")
  fi
  [[ -s "${SOURCE_IDENTITY}" ]] &&
    prepare_rds_args+=(--source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
  ${PIXI_RSCRIPT} "${SCRIPT_DIR}/validate_benchmark_rds_contract.R" \
    "${prepare_rds_args[@]}" >/dev/null 2>&1
}

benchmark_selected_artifacts_valid() {
  local ds="$1" view="$2" label="$3" path artifact_check batch_identity
  local has_feather=0 rds_grouped=0 group_rc
  local artifact_validator_args=()
  batch_identity=""
  if [[ "${PASS_ARG:-}" == corrected ]]; then
    batch_identity="$(
      stage5_batch_contract_identity_path "${ds}" "${view}" "${label}"
    )" || return 1
  fi
  case "${label}" in
    gloscope|mofa|pseudobulk|composition|scitd|prepare_pseudobulk|trans|zeroimp)
      if benchmark_rds_group_valid "${ds}" "${view}"; then
        rds_grouped=1
      else
        group_rc=$?
        [[ ${group_rc} -eq 2 ]] || return 1
      fi
  esac
  benchmark_artifacts_for "${ds}" "${view}" "${label}" || return 1
  [[ ${#ARTIFACT_PATHS[@]} -gt 0 ]] || return 1
  for path in "${ARTIFACT_PATHS[@]}"; do
    stage5_validate_reusable_artifact "${path}" "${label}" || return 1
    case "${path}" in
      *.feather)
        has_feather=1
        ;;
      *.rds)
        [[ ${rds_grouped} -eq 1 ]] && continue
        rds_args=(--artifact "${path}" --method "${label}" --dataset "${ds}" --view "${view}" \
          --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}")
        if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
          rds_args+=(--analysis-variant "${ANALYSIS_VARIANT}")
        fi
        [[ "${path}" == *_metadata.rds ]] && rds_args+=(--metadata)
        [[ -s "${SOURCE_IDENTITY}" ]] && rds_args+=(--source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
        [[ -n "${batch_identity}" ]] &&
          rds_args+=(--expected-batch-contract "${batch_identity}")
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
      artifact_validator_args=(--root "${ANALYSIS_ROOT}" --selection "${artifact_check}" \
        --labels "${label}" --batch --batch-pass "${PASS_ARG}" \
        --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}")
      if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
        artifact_validator_args+=(--analysis-variant "${ANALYSIS_VARIANT}")
      fi
      [[ -s "${SOURCE_IDENTITY}" ]] &&
        artifact_validator_args+=(--source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
      [[ -n "${batch_identity}" ]] &&
        artifact_validator_args+=(--expected-batch-contract "${batch_identity}")
      "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" \
        "${artifact_validator_args[@]}" >/dev/null 2>&1 || {
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
      artifact_validator_args=(--root "${ANALYSIS_ROOT}" --selection "${artifact_check}" \
        --labels "${label}" --input-root "${HPC_SCRATCH_DIR}" \
        --config "${DATASETS_JSON_FILE}")
      [[ -s "${SOURCE_IDENTITY}" ]] &&
        artifact_validator_args+=(--source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
      "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" \
        "${artifact_validator_args[@]}" >/dev/null 2>&1 || {
        rm -f "${artifact_check}" "${artifact_check}.md5"
        return 1
      }
      rm -f "${artifact_check}" "${artifact_check}.md5"
    fi
  fi
}

stage5_publish_output_records() {
  local ds view label path runtime_path
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  command -v ecoda_write_artifact_record >/dev/null 2>&1 || return 1
  while IFS=$'\t' read -r ds view label; do
    [[ -n "${ds}" && -n "${view}" && -n "${label}" ]] || return 1
    benchmark_artifacts_for "${ds}" "${view}" "${label}" || return 1
    for path in "${ARTIFACT_PATHS[@]}"; do
      stage5_publish_or_validate_artifact_record "${path}" "${label}" ||
        return 1
      case "${label}" in
        mrvi|scpoli|pilot|qot|pilotgm|trans|zeroimp)
          runtime_path="${path}.runtime.json"
          [[ -s "${runtime_path}" ]] || return 1
          stage5_publish_or_validate_artifact_record "${runtime_path}" "${label}" ||
            return 1
          ;;
      esac
    done
  done < "${PENDING_SELECTION}"
}

stage5_recorded_producer() {
  local record="$1" producer="" producer_count=0 key value
  [[ -f "${record}" && ! -L "${record}" && -r "${record}" ]] || return 1
  while IFS='=' read -r key value; do
    if [[ "${key}" == PRODUCER ]]; then
      producer_count=$((producer_count + 1))
      producer="${value}"
    fi
  done < "${record}"
  [[ ${producer_count} -eq 1 && -n "${producer}" &&
     "${producer}" =~ ^[A-Za-z0-9_.-]+$ ]] || return 1
  printf '%s' "${producer}"
}

stage5_allowed_producers_for() {
  local path="$1" label="$2" basename variant
  STAGE5_ALLOWED_PRODUCERS=()
  case "${label}" in
    prepare_pseudobulk)
      basename="${path##*/}"
      variant="${basename##*_pseudobulk_}"
      variant="${variant%.rds}"
      stage5_validate_path_component "${variant}" pseudobulk ||
        return 1
      STAGE5_ALLOWED_PRODUCERS=("stage5_prepare_pseudobulk_${variant}")
      ;;
    gloscope)
      if [[ -n "${PASS_ARG}" ]]; then
        STAGE5_ALLOWED_PRODUCERS=(stage5_gloscope)
      else
        STAGE5_ALLOWED_PRODUCERS=(stage5_gloscope_consolidate)
      fi
      ;;
    *)
      stage5_validate_path_component "${label}" producer ||
        return 1
      STAGE5_ALLOWED_PRODUCERS=("stage5_${label}")
      ;;
  esac
}

stage5_producer_allowed() {
  local path="$1" label="$2" producer="$3" allowed
  stage5_allowed_producers_for "${path}" "${label}" || return 1
  for allowed in "${STAGE5_ALLOWED_PRODUCERS[@]}"; do
    [[ "${producer}" == "${allowed}" ]] && return 0
  done
  return 1
}

stage5_publish_or_validate_artifact_record() {
  local path="$1" fallback_producer="$2" record producer canonical_producer
  stage5_allowed_producers_for "${path}" "${fallback_producer}" || return 1
  canonical_producer="${STAGE5_ALLOWED_PRODUCERS[0]}"
  record="$(ecoda_artifact_record_path "${path}" "${RUN_ID}" 2>/dev/null || true)"
  [[ -n "${record}" ]] || return 1
  if [[ -e "${record}" || -L "${record}" ]]; then
    producer="$(stage5_recorded_producer "${record}")" || return 1
    stage5_producer_allowed "${path}" "${fallback_producer}" "${producer}" ||
      return 1
    ecoda_validate_artifact_record "${path}" "${producer}" "${RUN_ID}" >/dev/null ||
      return 1
    return 0
  fi
  ecoda_write_artifact_record "${path}" "${canonical_producer}" "${RUN_ID}" >/dev/null
}

stage5_validate_pending_rds_records() {
  local ds view label path record producer
  while IFS=$'\t' read -r ds view label; do
    [[ -n "${ds}" && -n "${view}" && -n "${label}" ]] || return 1
    benchmark_artifacts_for "${ds}" "${view}" "${label}" || return 1
    for path in "${ARTIFACT_PATHS[@]}"; do
      case "${path}" in
        *.rds)
          record="$(ecoda_artifact_record_path \
            "${path}" "${RUN_ID}" 2>/dev/null || true)"
          [[ -n "${record}" && -e "${record}" ]] || return 1
          producer="$(stage5_recorded_producer "${record}")" || return 1
          stage5_producer_allowed "${path}" "${label}" "${producer}" || return 1
          ecoda_validate_artifact_record \
            "${path}" "${producer}" "${RUN_ID}" >/dev/null || return 1
          ;;
      esac
    done
  done < "${PENDING_SELECTION}"
}

stage5_matrix_validation_selection() {
  local label="$1" selection tmp
  STAGE5_MATRIX_VALIDATION_SELECTION="${MANIFEST}"
  STAGE5_MATRIX_VALIDATION_SELECTION_TEMP=""
  if [[ ${EXACT_SELECTION} -eq 1 && -z "${PASS_ARG}" ]]; then
    selection="${ECODA_RUN_ROOT}/manifests/matrix_validation_${label}.tsv"
    tmp="${selection}.build.$$"
    awk -F '\t' -v wanted="${label}" '$3 == wanted { print }' \
      "${MANIFEST}" > "${tmp}" || return 1
    [[ -s "${tmp}" ]] || {
      rm -f "${tmp}"
      return 2
    }
    ecoda_atomic_install_manifest "${tmp}" "${selection}" 3 || {
      rm -f "${tmp}"
      return 1
    }
    rm -f "${tmp}"
    ecoda_write_checksum "${selection}" || return 1
    STAGE5_MATRIX_VALIDATION_SELECTION="${selection}"
    STAGE5_MATRIX_VALIDATION_SELECTION_TEMP="${selection}"
  fi
}

stage5_validate_corrected_matrix_rows() {
  local label="$1" ds view row_label identity_path safe selection tmp
  local -a validation_args
  [[ "${PASS_ARG:-}" == corrected ]] || return 1
  while IFS=$'\t' read -r ds view row_label; do
    if [[ ${METHOD_MATRIX_MODE} -eq 1 ]] &&
       ! ecoda_stage5_method_matrix_allows "${ds}" "${view}" "${label}"; then
      continue
    fi
    identity_path="$(
      stage5_batch_contract_identity_path "${ds}" "${view}" "${label}"
    )" || return 1
    safe="$(_ecoda_safe_component "${ds}__${view}__${label}")" || return 1
    selection="${ECODA_RUN_ROOT}/manifests/matrix_validation_${safe}.tsv"
    tmp="${selection}.build.$$"
    printf '%s\t%s\t%s\n' "${ds}" "${view}" "${row_label}" > "${tmp}" || return 1
    ecoda_atomic_install_manifest "${tmp}" "${selection}" 3 || {
      rm -f "${tmp}"
      return 1
    }
    rm -f "${tmp}"
    ecoda_write_checksum "${selection}" || {
      rm -f "${selection}"
      return 1
    }
    validation_args=(--root "${ANALYSIS_ROOT}" --selection "${selection}" \
      --labels "${label}" --batch --batch-pass corrected \
      --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" \
      --source-identity "${SOURCE_IDENTITY}" --source-identity-verified \
      --expected-batch-contract "${identity_path}" \
      --producer "stage5_${label}" --producer-run-id "${RUN_ID}")
    [[ -n "${ANALYSIS_VARIANT:-}" ]] &&
      validation_args+=(--analysis-variant "${ANALYSIS_VARIANT}")
    if ! "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" \
        "${validation_args[@]}" >/dev/null 2>&1; then
      rm -f "${selection}" "${selection}.md5"
      return 1
    fi
    rm -f "${selection}" "${selection}.md5"
  done < "${MANIFEST}"
}

stage5_track_pending_artifact_owners() {
  local owner_dir reclaim_terminal=0
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  [[ -s "${PENDING_SELECTION}" ]] || return 0
  stage5_corrected_final_recovery_mode && reclaim_terminal=1
  ecoda_stage5_validate_output_ownership "${PENDING_SELECTION}" "${RUN_ID}" \
    "${reclaim_terminal}" ||
    return 1
  for owner_dir in "${ECODA_OUTPUT_OWNER_DIRS[@]:-}"; do
    [[ -n "${owner_dir}" ]] || return 1
    _ecoda_artifact_owner_validate_dir "${owner_dir}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_STATE:-}" == "ACTIVE" &&
       "${ECODA_ARTIFACT_OWNER_RUN:-}" == "${RUN_ID}" ]] || return 1
    ecoda_owner_track "${owner_dir}" || return 1
  done
}

stage5_selection_has_pending_rows() {
  local ds view method extra
  while IFS=$'\t' read -r ds view method extra; do
    [[ -n "${ds}" && -n "${view}" && -n "${method}" && -z "${extra}" ]] ||
      return 1
    if [[ ${METHOD_MATRIX_MODE} -eq 1 && "${ds}" == "Breast_cancer" ]]; then
      return 0
    fi
    if stage5_corrected_final_recovery_mode; then
      return 0
    fi
    if stage5_method_is_forced "${method}"; then
      return 0
    fi
    if [[ ${FORCE_TARGETED_ARG} -eq 1 &&
          "${method}" == "prepare_pseudobulk" ]]; then
      stage5_prepare_pseudobulk_valid "${ds}" "${view}" || return 0
    elif ! benchmark_selected_artifacts_valid "${ds}" "${view}" "${method}"; then
      return 0
    fi
  done < "${DISPATCH_MANIFEST}"
  return 1
}

stage5_variant_metadata() {
  if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
    ecoda_stage5_validate_identity "${PASS_ARG:-${ANALYSIS_PASS:-}}" || return 1
    printf 'ANALYSIS_VARIANT=%s\nANALYSIS_ROOT=%s\nANALYSIS_NAS_ROOT=%s\nANALYSIS_PASS=%s\nANALYSIS_LOG_PREFIX=%s\nMETADATA_EXPORT_MANIFEST=%s\nMETADATA_EXPORT_STATUS=%s\n' \
      "${ANALYSIS_VARIANT}" "${ANALYSIS_ROOT}" "${ANALYSIS_NAS_ROOT}" \
      "${ANALYSIS_PASS:-${PASS_ARG}}" "${ANALYSIS_LOG_PREFIX}" \
      "${ECODA_RUN_ROOT}/manifests/metadata_export.tsv" \
      "${ECODA_RUN_ROOT}/status/metadata_export.report"
    if [[ "${ANALYSIS_VARIANT}" == corrected_final &&
          "${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION:-}" == "recovery_35row" ]]; then
      printf 'ANALYSIS_ROOT_VERSION=%s\nANALYSIS_ROOT_IDENTITY=%s\n' \
        "${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION}" \
        "corrected_final/${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION}"
    fi
    if [[ ${METHOD_MATRIX_MODE} -eq 1 ]]; then
      [[ -s "${METHOD_MATRIX}" && -s "${METHOD_MATRIX}.md5" ]] || return 1
      METHOD_MATRIX_PENDING_COUNT=0
      if [[ -r "${PENDING_SELECTION:-}" ]]; then
        METHOD_MATRIX_PENDING_COUNT="$(awk 'END { print NR }' \
          "${PENDING_SELECTION}")" || return 1
      fi
      [[ "${METHOD_MATRIX_PENDING_COUNT}" =~ ^[0-9]+$ ]] || return 1
      printf 'METHOD_MATRIX=%s\nMETHOD_MATRIX_SOURCE=%s\nMETHOD_MATRIX_MD5=%s\nMETHOD_MATRIX_SIZE=%s\nMETHOD_MATRIX_SHA256=%s\nMETHOD_MATRIX_IDENTITY=%s\nDECLARED_METHOD_ROWS=%s\nPENDING_METHOD_ROWS=%s\n' \
        "${METHOD_MATRIX}" "${METHOD_MATRIX_SOURCE_PATH}" \
        "${METHOD_MATRIX_MD5}" "${METHOD_MATRIX_SIZE}" \
        "${METHOD_MATRIX_SHA256}" "${METHOD_MATRIX_IDENTITY}" \
        "${METHOD_MATRIX_DECLARED_COUNT}" "${METHOD_MATRIX_PENDING_COUNT}"
    fi
    if [[ "${ANALYSIS_VARIANT}" == corrected_final &&
          -s "${ECODA_RUN_ROOT}/manifests/corrected_final_consumer_contract.json" ]]; then
      printf 'CORRECTED_FINAL_CONSUMER_CONTRACT=%s\n' \
        "${ECODA_RUN_ROOT}/manifests/corrected_final_consumer_contract.json"
    fi
  fi
}

stage5_export_final_sample_metadata() {
  local manifest="${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
  local manifest_tmp="${manifest}.build.$$"
  local pending_manifest="${ECODA_RUN_ROOT}/manifests/metadata_export_pending.tsv"
  local pending_tmp="${pending_manifest}.build.$$"
  local status_dir="${ECODA_RUN_ROOT}/status/metadata_export"
  local status_report="${ECODA_RUN_ROOT}/status/metadata_export.report"
  local exporter worker export_runtime export_msg export_id export_rc
  local ds view input output safe state status_task pending_count=0 pending_index=0
  local row extra seen="" required_view="batch_effect_${PASS_ARG:-}"
  [[ -n "${ANALYSIS_VARIANT:-}" ]] || return 0
  ecoda_stage5_validate_identity "${PASS_ARG:-${ANALYSIS_PASS:-}}" || return 1
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == "1" ]] && return 0
  exporter="$(stage5_source_script src/utils/py/export_h5ad_sample_metadata.py)"
  stage5_require_source_script "${exporter}" || return 1
  if [[ -n "${SYNC_ONLY_RUN}" ]]; then
    ecoda_validate_run_owned_path "${manifest}" "${ECODA_RUN_ROOT}" || return 1
    ecoda_validate_manifest "${manifest}" 4 || return 1
    ecoda_validate_checksum "${manifest}" || return 1
    while IFS=$'\t' read -r ds view input output extra; do
      [[ -n "${ds}" && "${view}" == "${required_view}" &&
         -n "${input}" && -n "${output}" && -z "${extra}" ]] || return 1
      ecoda_validate_checksum "${output}" || return 1
    done < "${manifest}"
    return 0
  fi
  : > "${manifest_tmp}" || return 1
  : > "${pending_tmp}" || {
    rm -f "${manifest_tmp}"
    return 1
  }
  while IFS=$'\t' read -r ds view _row_label extra; do
    [[ -n "${ds}" && "${view}" == "${required_view}" &&
       -z "${extra}" ]] || {
      rm -f "${manifest_tmp}" "${pending_tmp}"
      return 1
    }
    case " ${seen} " in *" ${ds} "*) continue ;; esac
    seen="${seen} ${ds}"
    input="$(stage5_input_path "${ds}" "${view}")" || {
      rm -f "${manifest_tmp}" "${pending_tmp}"
      return 1
    }
    output="${ANALYSIS_ROOT}/metadata/${ds}_sample_metadata.feather"
    printf '%s\t%s\t%s\t%s\n' "${ds}" "${view}" "${input}" "${output}" >> "${manifest_tmp}" ||
      return 1
    if ! "${PYTHON_BIN}" "${exporter}" \
        --config "${DATASETS_JSON_FILE}" --dataset "${ds}" --view "${view}" \
        --input-file "${input}" --output "${output}" --check >/dev/null 2>&1; then
      printf '%s\t%s\t%s\t%s\n' "${ds}" "${view}" "${input}" "${output}" >> "${pending_tmp}" ||
        return 1
      pending_count=$((pending_count + 1))
    fi
  done < "${MANIFEST}"
  [[ -s "${manifest_tmp}" ]] || {
    rm -f "${manifest_tmp}" "${pending_tmp}"
    return 1
  }
  ecoda_atomic_install_manifest "${manifest_tmp}" "${manifest}" 4 || {
    rm -f "${manifest_tmp}" "${pending_tmp}"
    return 1
  }
  rm -f "${manifest_tmp}"
  ecoda_write_checksum "${manifest}" || {
    rm -f "${pending_tmp}"
    return 1
  }
  mkdir -p "${status_dir}" || {
    rm -f "${pending_tmp}"
    return 1
  }
  ecoda_validate_run_owned_path "${status_dir}" "${ECODA_RUN_ROOT}" || return 1
  rm -f "${status_dir}"/*.status
  if [[ ${pending_count} -gt 0 ]]; then
    ecoda_atomic_install_manifest "${pending_tmp}" "${pending_manifest}" 4 || {
      rm -f "${pending_tmp}"
      return 1
    }
    rm -f "${pending_tmp}"
    ecoda_write_checksum "${pending_manifest}" || return 1
    worker="$(stage5_source_script src/utils/bash/h5ad_obs_audit_worker.sh)"
    stage5_validate_bound_runtime || return 1
    stage5_require_source_script "${worker}" || return 1
    export_runtime="${RUNTIME_EXPORT}"
    [[ -n "${ECODA_SOURCE_ROOT:-}" ]] &&
      export_runtime="${export_runtime},ECODA_SOURCE_ROOT=${ECODA_SOURCE_ROOT},ECODA_SOURCE_MANIFEST=${ECODA_SOURCE_MANIFEST},ECODA_SOURCE_SNAPSHOT_REQUIRED=1"
    export_runtime="${export_runtime},ANALYSIS_VARIANT=${ANALYSIS_VARIANT},ANALYSIS_PASS=${PASS_ARG},ANALYSIS_ROOT=${ANALYSIS_ROOT},ANALYSIS_NAS_ROOT=${ANALYSIS_NAS_ROOT},ANALYSIS_LOG_PREFIX=${ANALYSIS_LOG_PREFIX}"
    set +e
    export_msg="$(sbatch --parsable --wait \
      --array="1-${pending_count}%${THROTTLE}" \
      --partition="${PARTITION_ARG:-${SLURM_PARTITION_BENCHMARK_CPU}}" \
      --ntasks=1 --cpus-per-task=1 --mem="${MEMORY}" \
      --time="${WATCHDOG_TIME_LIMIT}" \
      --output="${ECODA_RUN_ROOT}/logs/metadata_export_%A_%a.log" \
      --error="${ECODA_RUN_ROOT}/logs/metadata_export_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      --export="ALL,H5AD_OBS_AUDIT_MODE=metadata,H5AD_METADATA_EXPORT_MANIFEST=${pending_manifest},H5AD_METADATA_EXPORT_STATUS_DIR=${status_dir},${export_runtime}" \
      "${worker}")"
    export_rc=$?
    set -e
    export_id="${export_msg%%;*}"
    [[ "${export_id}" =~ ^[0-9]+$ ]] || return 1
    stage5_record_scheduler METADATA_EXPORT "${export_id}" || return 1
    [[ ${export_rc} -eq 0 ]] || return 1
  else
    rm -f "${pending_tmp}" "${pending_manifest}" "${pending_manifest}.md5"
  fi
  pending_index=0
  while IFS=$'\t' read -r ds view input output extra; do
    [[ -n "${ds}" && "${view}" == "${required_view}" &&
       -n "${input}" && -n "${output}" && -z "${extra}" ]] || return 1
    safe="$(_ecoda_safe_component "${ds}__${view}")"
    status="${status_dir}/${safe}.status"
    if "${PYTHON_BIN}" "${exporter}" \
        --config "${DATASETS_JSON_FILE}" --dataset "${ds}" --view "${view}" \
        --input-file "${input}" --output "${output}" --check >/dev/null 2>&1; then
      ecoda_validate_checksum "${output}" || return 1
      if [[ ${pending_count} -gt 0 && -s "${pending_manifest}" ]]; then
        if grep -F -q "${ds}" "${pending_manifest}"; then
          pending_index=$((pending_index + 1))
          [[ -s "${status}" ]] || return 1
          state="$(sed -n 's/^STATE=//p' "${status}" | head -1 || true)"
          status_task="$(sed -n 's/^TASK_ID=//p' "${status}" | head -1 || true)"
          [[ "${state}" == OK && "${status_task}" == "${pending_index}" ]] || return 1
        fi
      else
        ecoda_atomic_write "${status}" \
          "STATE=OK\nSTATUS=NOOP_VALIDATED\nRUN_ID=${RUN_ID}\nDATASET=${ds}\nVIEW=${view}\nTASK_ID=0\nINPUT_FILE=${input}\nOUTPUT_FILE=${output}\n" || return 1
      fi
    else
      return 1
    fi
  done < "${manifest}"
  [[ ${pending_count} -eq 0 || ${pending_index} -eq ${pending_count} ]] || return 1
  ecoda_atomic_write "${status_report}" \
    "STATE=OK\nANALYSIS_VARIANT=${ANALYSIS_VARIANT}\nANALYSIS_ROOT=${ANALYSIS_ROOT}\nANALYSIS_NAS_ROOT=${ANALYSIS_NAS_ROOT}\nANALYSIS_PASS=${PASS_ARG}\nANALYSIS_LOG_PREFIX=${ANALYSIS_LOG_PREFIX}\nRUN_ID=${RUN_ID}\nMANIFEST=${manifest}\nCOUNT=$(awk 'END { print NR }' "${manifest}")\nPENDING=${pending_count}\n" || return 1
  ecoda_write_checksum "${status_report}" || return 1
  return 0
}

stage5_export_final_sample_metadata ||
  stage5_abort "final sample metadata export failed"


PENDING_SELECTION="${ECODA_RUN_ROOT}/manifests/pending_selection.tsv"
DISPATCH_ROWS="$(awk 'END { print NR }' "${DISPATCH_MANIFEST}")" ||
  stage5_abort "failed to count Stage 5 dispatch manifest rows"
dispatch_metadata="DISPATCH_SELECTION=${DISPATCH_MANIFEST}\nDISPATCH_SELECTION_MD5=${DISPATCH_MANIFEST_MD5}\nDISPATCH_SELECTION_SIZE=${DISPATCH_MANIFEST_SIZE}\nDISPATCH_SELECTION_ROWS=${DISPATCH_ROWS}\n"
if [[ -z "${SYNC_ONLY_RUN}" ]] &&
   ! stage5_selection_has_pending_rows; then
  ecoda_atomic_write "${PENDING_SELECTION}" "" ||
    stage5_abort "failed to create empty Stage 5 pending manifest"
  ecoda_atomic_write "${ECODA_RUN_ROOT}/manifests/owners.tsv" "" ||
    stage5_abort "failed to create empty Stage 5 owner manifest"
  identity_metadata="$(stage5_record_run_identity_metadata)" ||
    stage5_abort "failed to read Stage 5 source/runtime identity"
  [[ -z "${identity_metadata}" ]] ||
    identity_metadata="${identity_metadata}"$'\n'
  batch_contract_metadata="$(stage5_record_batch_contract_metadata)" ||
    stage5_abort "failed to read corrected Stage 5 batch-contract metadata"
  batch_contract_metadata_suffix=""
  [[ "${PASS_ARG:-}" == corrected ]] &&
    batch_contract_metadata_suffix="${batch_contract_metadata}\n"
  methods_csv="$(IFS=,; echo "${METHODS[*]}")"
  analyses_csv=""
  [[ ${ANALYSES_SELECTED} -eq 1 ]] &&
    analyses_csv="$(IFS=,; echo "${ANALYSES[*]}")"
  target_methods_csv=""
  if [[ ${TARGET_METHODS_SET} -eq 1 ]]; then
    target_methods_csv="$(IFS=,; echo "${TARGET_METHODS[*]}")"
  fi
  target_methods_metadata=""
  if [[ -n "${ANALYSIS_VARIANT:-}" && ${TARGET_METHODS_SET} -eq 1 ]]; then
    target_methods_metadata="TARGET_METHODS=${target_methods_csv}\n"
  fi
  analysis_variant_metadata="$(stage5_variant_metadata)"
  [[ -z "${analysis_variant_metadata}" ]] ||
    analysis_variant_metadata="${analysis_variant_metadata}"$'\n'
  root_metadata="ROOT=${ANALYSIS_ROOT}\n"
  ecoda_atomic_write "${ECODA_RUN_ROOT}/metadata" \
    "STAGE=stage5\nRUN_ID=${RUN_ID}\nSTATE=ACTIVE\nMETHODS=${methods_csv}\nANALYSES=${analyses_csv}\nPASS=${PASS_ARG}\nEXACT_SELECTION=${EXACT_SELECTION}\nFORCE_TARGETED=${FORCE_TARGETED_ARG}\nFORCE_REASON=${FORCE_REASON_ARG}\nSOURCE_IDENTITY=${SOURCE_IDENTITY}\n${dispatch_metadata}${analysis_variant_metadata}${target_methods_metadata}${root_metadata}${identity_metadata}\n${batch_contract_metadata_suffix}" ||
    stage5_abort "failed to write Stage 5 NOOP metadata"
  ecoda_atomic_write "${ECODA_RUN_ROOT}/status/report" \
    "STATE=NOOP_VALIDATED\nRUN_ID=${RUN_ID}\nFORCE_TARGETED=${FORCE_TARGETED_ARG}\nFORCE_REASON=${FORCE_REASON_ARG}\nREASON=all selected benchmark artifacts are valid; no rerun selected\n" ||
    stage5_abort "failed to write Stage 5 NOOP_VALIDATED report"
  ecoda_set_run_state OK "NOOP_VALIDATED: all selected benchmark artifacts are valid" ||
    stage5_abort "failed to write Stage 5 NOOP_VALIDATED terminal state"
  echo "NOOP_VALIDATED=1"
  if [[ -n "${PASS_ARG}" ]]; then
    echo "BATCH_EFFECT_RUN_ID=${RUN_ID}"
  else
    echo "BENCHMARK_RUN_ID=${RUN_ID}"
  fi
  exit 0
fi


if [[ -z "${SYNC_ONLY_RUN}" ]]; then
  stage5_repair_missing_h5ad_sidecars ||
    stage5_abort "Stage 5 source H5AD sidecar validation failed"
fi
stage5_compute_h5ad_preflight ||
  stage5_abort "Stage 5 compute-node H5AD preflight failed"
if [[ -z "${SYNC_ONLY_RUN}" ||
      -s "${ECODA_RUN_ROOT}/manifests/h5ad_preflight.tsv" ]]; then
  stage5_validate_source_artifact_records ||
    stage5_abort "Stage 5 H5AD preflight artifact record validation failed"
fi
stage5_prepare_source_identity ||
  stage5_abort "failed to build or verify Stage 5 source identity"

stage5_run_r_environment_preflight() {
  local log_dir="${ECODA_RUN_ROOT}/logs/r_environment_preflight"
  local preflight_msg preflight_id preflight_rc preflight_worker preflight_export
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0
  mkdir -p "${log_dir}" || return 1
  preflight_worker="$(stage5_source_script src/utils/bash/r_environment_preflight_worker.sh)"
  preflight_export="${RUNTIME_EXPORT},ECODA_RUN_ROOT=${ECODA_RUN_ROOT},ECODA_RUN_ID=${RUN_ID},R_ENV_PREFLIGHT_RUN_ROOT=${ECODA_RUN_ROOT}"
  [[ -n "${ECODA_SOURCE_ROOT:-}" ]] &&
    preflight_export="${preflight_export},ECODA_SOURCE_ROOT=${ECODA_SOURCE_ROOT},ECODA_SOURCE_MANIFEST=${ECODA_SOURCE_MANIFEST},ECODA_SOURCE_SNAPSHOT_REQUIRED=1"
  stage5_validate_bound_runtime || return 1
  stage5_require_source_script "${preflight_worker}" || return 1
  set +e
  preflight_msg="$(sbatch --parsable --wait \
    --partition="${SLURM_PARTITION_BENCHMARK_CPU}" \
    --ntasks=1 --cpus-per-task=1 --mem=2G \
    --time="${WATCHDOG_TIME_LIMIT}" \
    --output="${log_dir}/r_environment_preflight_%j.log" \
    --error="${log_dir}/r_environment_preflight_%j.err" \
    --mail-user="${USER_EMAIL}" --export="ALL,${preflight_export}" \
    "${preflight_worker}")"
  preflight_rc=$?
  set -e
  preflight_id="${preflight_msg%%;*}"
  [[ "${preflight_id}" =~ ^[0-9]+$ ]] || return 1
  stage5_record_scheduler PREFLIGHT "${preflight_id}" || return 1
  [[ ${preflight_rc} -eq 0 ]] || return 1
  echo "BENCHMARK_R_ENV_PREFLIGHT_JOB_ID=${preflight_id}"
}

if [[ -z "${SYNC_ONLY_RUN}" && ${R_ENV_PREFLIGHT_REQUIRED} -eq 1 ]]; then
  stage5_run_r_environment_preflight ||
    stage5_abort "Stage 5 compute-node R environment preflight failed"
fi

stage5_run_corrected_final_consumer_barrier() {
  local validator report
  local -a consumer_args=()
  [[ "${BENCHMARK_MATRIX_TEST:-0}" == 1 ]] && return 0
  [[ "${ANALYSIS_VARIANT:-}" == corrected_final &&
     -z "${SYNC_ONLY_RUN}" ]] || return 0
  validator="$(
    stage5_source_script \
      src/5_run_benchmark_methods/validate_corrected_final_consumer_contracts.R
  )" || return 1
  stage5_require_source_script "${validator}" || return 1
  report="${ECODA_RUN_ROOT}/manifests/corrected_final_consumer_contract.json"
  stage5_validate_bound_runtime || return 1
  consumer_args=(
    --config "${DATASETS_JSON_FILE}"
    --selection "${MANIFEST}"
    --analysis-root "${ANALYSIS_ROOT}"
    --input-root "${HPC_SCRATCH_DIR}"
    --output "${report}"
  )
  if [[ ${METHOD_MATRIX_MODE} -eq 1 ]]; then
    consumer_args+=(
      --method-matrix "${METHOD_MATRIX}"
      --method-matrix-md5 "${METHOD_MATRIX_MD5}"
      --method-matrix-size "${METHOD_MATRIX_SIZE}"
      --method-matrix-sha256 "${METHOD_MATRIX_SHA256}"
      --method-matrix-identity "${METHOD_MATRIX_IDENTITY}"
    )
  fi
  ${PIXI_RSCRIPT} "${validator}" "${consumer_args[@]}" >/dev/null 2>&1 || return 1
  [[ -s "${report}" && -s "${report}.md5" ]] || return 1
  ecoda_validate_checksum "${report}" || return 1
}

stage5_run_corrected_final_consumer_barrier ||
  stage5_abort "corrected-final consumer contract barrier failed"

OWNERS_FILE="${ECODA_RUN_ROOT}/manifests/owners.tsv"
PENDING_SELECTION_MD5=""
PENDING_SELECTION_SIZE=""
if [[ -z "${SYNC_ONLY_RUN}" ]]; then
  OWNERS_TMP="${OWNERS_FILE}.build.$$"
  PENDING_SELECTION="${ECODA_RUN_ROOT}/manifests/pending_selection.tsv"
  PENDING_SELECTION_TMP="${PENDING_SELECTION}.build.$$"
  : > "${OWNERS_TMP}"
  : > "${PENDING_SELECTION_TMP}"
  OWNER_SEEN=""
  ecoda_owner_clear_tracked
  stage5_add_pending_row() {
    local ds="$1" view="$2" method="$3"
    local owner_key owner_dir owner_rc method_force=0
    local matrix_breast_pending=0 path
    [[ "${method}" == _ecoda_none_ ]] && return 0
    stage5_method_is_forced "${method}" && method_force=1
    [[ ${METHOD_MATRIX_MODE} -eq 1 && "${ds}" == "Breast_cancer" ]] &&
      matrix_breast_pending=1
    owner_key="${PASS_ARG:-ordinary}/${ds}/${view}/${method}"
    case " ${OWNER_SEEN} " in
      *" ${owner_key} "*) return 0 ;;
    esac
    OWNER_SEEN="${OWNER_SEEN} ${owner_key}"
    if [[ "${method}" == "prepare_pseudobulk" &&
          ${method_force} -eq 0 && ${matrix_breast_pending} -eq 0 ]] &&
       [[ ${FORCE_TARGETED_ARG} -eq 1 || ${TARGET_METHODS_SET} -eq 1 ]]; then
      if stage5_corrected_final_recovery_mode; then
        :
      elif stage5_prepare_pseudobulk_valid "${ds}" "${view}"; then
        echo "Skipping validated Stage 5 pseudobulk cache ${ds}/${view}/${method}."
        return 0
      fi
    elif [[ ${method_force} -eq 0 &&
            ${matrix_breast_pending} -eq 0 ]] &&
         ! stage5_corrected_final_recovery_mode &&
         benchmark_selected_artifacts_valid "${ds}" "${view}" "${method}"; then
      echo "Skipping validated Stage 5 artifact ${ds}/${view}/${method}."
      return 0
    fi
    benchmark_artifacts_for "${ds}" "${view}" "${method}" ||
      stage5_abort "cannot resolve output contract for ${ds}/${view}/${method}"
    if [[ ${METHOD_MATRIX_MODE} -eq 1 &&
          "${ANALYSIS_VARIANT:-}" == corrected_final ]]; then
      for path in "${ARTIFACT_PATHS[@]}"; do
        for existing_target in "${path}" "${path}.md5" \
          "${path}.runtime.json" "${path}.runtime.json.md5"; do
          [[ ! -e "${existing_target}" && ! -L "${existing_target}" ]] || {
            stage5_abort "corrected-final method matrix target already exists; refusing to invalidate ${existing_target}"
          }
        done
      done
    fi
    set +e
    owner_dir="$(ecoda_owner_acquire stage5 "${owner_key}" "${RUN_ID}" "${method_force}" 0)"
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
  }
  while IFS=$'\t' read -r ds view method extra; do
    [[ -n "${ds}" && -n "${view}" && -n "${method}" && -z "${extra}" ]] ||
      stage5_abort "invalid Stage 5 dispatch manifest row"
    stage5_add_pending_row "${ds}" "${view}" "${method}"
  done < "${DISPATCH_MANIFEST}"
  if [[ -s "${PENDING_SELECTION_TMP}" ]]; then
    ecoda_atomic_install_manifest "${PENDING_SELECTION_TMP}" "${PENDING_SELECTION}" 3 ||
      stage5_abort "failed to install Stage 5 pending manifest atomically"
    ecoda_write_checksum "${PENDING_SELECTION}" ||
      stage5_abort "failed to checksum Stage 5 pending manifest"
    PENDING_SELECTION_MD5="${ECODA_CHECKSUM_MD5}"
    PENDING_SELECTION_SIZE="${ECODA_CHECKSUM_SIZE}"
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
  target_methods_csv=""
  if [[ ${TARGET_METHODS_SET} -eq 1 ]]; then
    target_methods_csv="$(IFS=,; echo "${TARGET_METHODS[*]}")"
  fi
  analysis_variant_metadata="$(stage5_variant_metadata)"
  [[ -z "${analysis_variant_metadata}" ]] ||
    analysis_variant_metadata="${analysis_variant_metadata}"$'\n'
  target_methods_metadata=""
  [[ -n "${ANALYSIS_VARIANT:-}" ]] &&
    target_methods_metadata="TARGET_METHODS=${target_methods_csv}\n"
  identity_metadata="$(stage5_record_run_identity_metadata)" ||
    stage5_abort "failed to read Stage 5 source/runtime identity"
  [[ -z "${identity_metadata}" ]] ||
    identity_metadata="${identity_metadata}"$'\n'
  batch_contract_metadata="$(stage5_record_batch_contract_metadata)" ||
    stage5_abort "failed to read corrected Stage 5 batch-contract metadata"
  batch_contract_metadata_suffix=""
  [[ "${PASS_ARG:-}" == corrected ]] &&
    batch_contract_metadata_suffix="${batch_contract_metadata}\n"
  dispatch_rows="$(awk 'END { print NR }' "${DISPATCH_MANIFEST}")" ||
    stage5_abort "failed to count Stage 5 dispatch manifest rows"
  dispatch_metadata="DISPATCH_SELECTION=${DISPATCH_MANIFEST}\nDISPATCH_SELECTION_MD5=${DISPATCH_MANIFEST_MD5}\nDISPATCH_SELECTION_SIZE=${DISPATCH_MANIFEST_SIZE}\nDISPATCH_SELECTION_ROWS=${dispatch_rows}\n"
  root_metadata="ROOT=${ANALYSIS_ROOT}\n"
  RUN_METADATA="STAGE=stage5\nRUN_ID=${RUN_ID}\nSTATE=ACTIVE\nMETHODS=${methods_csv}\nANALYSES=${analyses_csv}\nPASS=${PASS_ARG}\nEXACT_SELECTION=${EXACT_SELECTION}\nFORCE_TARGETED=${FORCE_TARGETED_ARG}\nFORCE_REASON=${FORCE_REASON_ARG}\nSOURCE_IDENTITY=${SOURCE_IDENTITY}\nH5AD_PREFLIGHT=${ECODA_RUN_ROOT}/manifests/h5ad_preflight.tsv\nPENDING_SELECTION=${PENDING_SELECTION}\nPENDING_SELECTION_MD5=${PENDING_SELECTION_MD5}\nPENDING_SELECTION_SIZE=${PENDING_SELECTION_SIZE}\n${dispatch_metadata}${analysis_variant_metadata}${target_methods_metadata}${root_metadata}${identity_metadata}\n${batch_contract_metadata_suffix}"
  ecoda_atomic_write "${ECODA_RUN_ROOT}/metadata" "${RUN_METADATA}" ||
    stage5_abort "failed to write Stage 5 run metadata"
fi


stage5_configure_analysis_context ||
  stage5_abort "Stage 5 analysis context is not bound to the selected variant"
if [[ -n "${PASS_ARG}" ]]; then
  export ANALYSIS_HIGH_RES_ONLY=1
else
  unset ANALYSIS_HIGH_RES_ONLY
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
  GROUP_STATUS_LABELS=()
  GROUP_RESOURCE_CLASSES=()
  GROUP_MANIFEST_COLUMNS=()
  GROUP_MANIFESTS=()
  GROUP_TMP_FILES=()
  PREP_VIEWS=()
  PREP_WATCHDOGS=()

  group_add_row() {
    local ds="$1" view="$2" label="$3" combo="${4:-}" resource_class="${5:-base}"
    local key idx=0 safe status_label
    if [[ "${resource_class}" == base ]]; then
      key="${view}|${label}"
    else
      key="${view}|${label}|${resource_class}"
    fi
    while [[ ${idx} -lt ${#GROUP_KEYS[@]} && "${GROUP_KEYS[${idx}]}" != "${key}" ]]; do
      idx=$((idx + 1))
    done
    if [[ ${idx} -eq ${#GROUP_KEYS[@]} ]]; then
      safe="$(printf '%s' "${key}" | tr '/:,\t |' '______')"
      status_label="${view}__${label}"
      [[ "${resource_class}" == base ]] ||
        status_label="${status_label}__${resource_class}"
      GROUP_KEYS+=("${key}")
      GROUP_VIEWS+=("${view}")
      GROUP_LABELS+=("${label}")
      GROUP_STATUS_LABELS+=("${status_label}")
      GROUP_RESOURCE_CLASSES+=("${resource_class}")
      GROUP_MANIFEST_COLUMNS+=("$([[ -n "${combo}" ]] && printf 4 || printf 3)")
      GROUP_MANIFESTS+=("${ECODA_RUN_ROOT}/manifests/matrix_${safe}.tsv")
      GROUP_TMP_FILES+=("${ECODA_RUN_ROOT}/manifests/matrix_${safe}.build.$$")
      : > "${GROUP_TMP_FILES[${idx}]}"
    fi
    if [[ -n "${combo}" ]]; then
      [[ "${GROUP_MANIFEST_COLUMNS[${idx}]}" == 4 ]] || return 1
      printf '%s\t%s\t%s\t%s\n' "${ds}" "${view}" "${label}" "${combo}" \
        >> "${GROUP_TMP_FILES[${idx}]}"
    else
      [[ "${GROUP_MANIFEST_COLUMNS[${idx}]}" == 3 ]] || return 1
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${label}" \
        >> "${GROUP_TMP_FILES[${idx}]}"
    fi
  }

  submit_matrix() {
    local group_label="$1" manifest="$2" dependency="$3" method="$4"
    local view="$5" resource_class="$6"
    method_spec_for_group "${method}" "${resource_class}" || return 1
    local method_force=0
    stage5_method_is_forced "${method}" && method_force=1
    local worker="${METHOD_WORKER}" worker_partition="${METHOD_PARTITION}" \
      watchdog_partition="${SLURM_PARTITION_BENCHMARK_CPU}" throttle="${METHOD_THROTTLE}" \
      method_time_limit="${METHOD_TIME_LIMIT}" method_gpu_policy="${METHOD_GPU_POLICY}"
    local method_runtime_export watchdog_script ownership_manifest ownership_tmp
    method_runtime_export="$(ecoda_runtime_export_csv stage5 "${METHOD_RUNTIME_NV}")" || return 1
    watchdog_script="$(stage5_source_script src/5_run_benchmark_methods/matrix_watchdog.sh)"
    local safe="$(printf '%s' "${group_label}" | tr '/:,\t ' '_____')"
    local array_msg array_id array_rc wd_msg wd_id wd_rc
    local worker_env="METHOD=${method},ANALYSIS=${method},ANALYSIS_MANIFEST=${manifest},ANALYSIS_VIEW=${view},ANALYSIS_ROOT=${ANALYSIS_ROOT},EXECUTION_LOG_DIR=${RUN_LOG_DIR},ECODA_RUN_ROOT=${ECODA_RUN_ROOT},ECODA_RUN_ID=${RUN_ID},ECODA_SELECTION_MANIFEST=${manifest},FORCE_BENCHMARK=${method_force},METHOD_TIME_LIMIT=${method_time_limit},METHOD_GPU_POLICY=${method_gpu_policy},ECODA_ARTIFACT_PRODUCER=stage5_${method},JOB_LOG_PREFIX=${RUN_LOG_DIR}/5_matrix_${safe}"
    if [[ -n "${PASS_ARG}" ]]; then
      worker_env="${worker_env},ANALYSIS_PASS=${PASS_ARG}"
      if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
        worker_env="${worker_env},ANALYSIS_VARIANT=${ANALYSIS_VARIANT},ANALYSIS_NAS_ROOT=${ANALYSIS_NAS_ROOT},ANALYSIS_LOG_PREFIX=${ANALYSIS_LOG_PREFIX}"
    if [[ ${METHOD_MATRIX_MODE} -eq 0 &&
          -n "${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION:-}" ]]; then
      worker_env="${worker_env},ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION=${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION}"
    fi
      fi
      [[ "${PASS_ARG}" == corrected ]] &&
        worker_env="${worker_env},ECODA_BATCH_CONTRACT_MANIFEST=${ECODA_BATCH_CONTRACT_MANIFEST}"
    else
      worker_env="${worker_env},BENCHMARK_MANIFEST=${manifest}"
    fi
    if [[ ${METHOD_MATRIX_MODE} -eq 1 ]]; then
      worker_env="${worker_env},METHOD_MATRIX=${METHOD_MATRIX},ECODA_STAGE5_METHOD_MATRIX=${METHOD_MATRIX},ECODA_STAGE5_METHOD_MATRIX_MODE=1,METHOD_MATRIX_IDENTITY=${METHOD_MATRIX_IDENTITY},METHOD_MATRIX_MD5=${METHOD_MATRIX_MD5},METHOD_MATRIX_SIZE=${METHOD_MATRIX_SIZE},ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION=${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION:-}"
    fi
    worker_env="${worker_env},${method_runtime_export}"
    array_args=(--parsable --array="1-$(wc -l < "${manifest}" | tr -d '[:space:]')%${throttle}" --partition="${worker_partition}" "${METHOD_FLAGS[@]}" --time="${method_time_limit}" --mem="${MEMORY}" \
      --output="${RUN_LOG_DIR}/5_matrix_${safe}_%A_%a.log" --error="${RUN_LOG_DIR}/5_matrix_${safe}_%A_%a.err" --mail-user="${USER_EMAIL}" \
      --export="ALL,${worker_env}" "${worker}")
    [[ -n "${dependency}" ]] && array_args+=(--dependency="afterok:${dependency}")
    stage5_validate_bound_runtime || return 1
    stage5_validate_output_ownership "${manifest}" 1 || return 1
    stage5_require_source_script "${worker}" || return 1
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
    stage5_validate_bound_runtime || return 1
    stage5_validate_output_ownership "${manifest}" || return 1
    stage5_require_source_script "${watchdog_script}" || return 1
    set +e
    wd_msg="$(sbatch --parsable --dependency="afterany:${array_id}" --partition="${watchdog_partition}" --ntasks=1 --cpus-per-task=1 --mem=2G --time="${WATCHDOG_TIME_LIMIT}" \
      --output="${RUN_LOG_DIR}/5_matrix_watchdog_${safe}_%A.log" --error="${RUN_LOG_DIR}/5_matrix_watchdog_${safe}_%A.err" --mail-user="${USER_EMAIL}" \
      --export="ALL,${worker_env},MATRIX_WATCHDOG_ROOT=${ECODA_RUN_ROOT}" \
      "${watchdog_script}" "${ECODA_RUN_ROOT}" "${group_label}" "${manifest}" "${array_id}" "${MEMORY}" "${MAX_MEMORY}" "${worker_partition}" "${throttle}" "${worker}" "${method_runtime_export}" "${METHOD_FLAGS[@]}")"
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
      shard_rows="$(benchmark_method_shard_rows "${method}" "${view}")"
      if [[ -n "${shard_rows}" ]]; then
        while IFS=$'\t' read -r shard_token resource_class; do
          [[ -n "${shard_token}" && -n "${resource_class}" ]] ||
            stage5_abort "invalid parameter-shard row for ${method}"
          group_add_row "${ds}" "${view}" "${method}" \
            "${shard_token}" "${resource_class}" ||
            stage5_abort "inconsistent parameter-shard manifest for ${method}"
        done <<< "${shard_rows}"
      else
        group_add_row "${ds}" "${view}" "${method}" "" base ||
          stage5_abort "failed to build Stage 5 matrix manifest"
      fi
    done < "${PENDING_SELECTION}"
  fi
  for idx in "${!GROUP_KEYS[@]}"; do
    ecoda_atomic_install_manifest "${GROUP_TMP_FILES[${idx}]}" "${GROUP_MANIFESTS[${idx}]}" \
      "${GROUP_MANIFEST_COLUMNS[${idx}]}" ||
      stage5_abort "failed to install Stage 5 matrix manifest"
    ecoda_write_checksum "${GROUP_MANIFESTS[${idx}]}" ||
      stage5_abort "failed to checksum Stage 5 matrix manifest"
    rm -f "${GROUP_TMP_FILES[${idx}]}"
  done

  for idx in "${!GROUP_KEYS[@]}"; do
    view="${GROUP_VIEWS[${idx}]}"
    method="${GROUP_LABELS[${idx}]}"
    group_label="${GROUP_STATUS_LABELS[${idx}]}"
    resource_class="${GROUP_RESOURCE_CLASSES[${idx}]}"
    if [[ "${method}" == prepare_pseudobulk ]]; then
      submit_matrix "${group_label}" "${GROUP_MANIFESTS[${idx}]}" "" \
        "${method}" "${view}" "${resource_class}" ||
        stage5_abort "submission failed for ${method}/${view}"
      PREP_VIEWS+=("${view}")
      PREP_WATCHDOGS+=("${LAST_WATCHDOG_ID}")
    fi
  done
  for idx in "${!GROUP_KEYS[@]}"; do
    view="${GROUP_VIEWS[${idx}]}"
    method="${GROUP_LABELS[${idx}]}"
    [[ "${method}" == prepare_pseudobulk ]] && continue
    group_label="${GROUP_STATUS_LABELS[${idx}]}"
    resource_class="${GROUP_RESOURCE_CLASSES[${idx}]}"
    dependency=""
    case "${method}" in
      mofa|pseudobulk|composition)
        for prep_idx in "${!PREP_VIEWS[@]}"; do
          [[ "${PREP_VIEWS[${prep_idx}]}" == "${view}" ]] &&
            dependency="${PREP_WATCHDOGS[${prep_idx}]}"
        done
        ;;
    esac
    submit_matrix "${group_label}" "${GROUP_MANIFESTS[${idx}]}" "${dependency}" \
      "${method}" "${view}" "${resource_class}" ||
      stage5_abort "submission failed for ${method}/${view}"
  done
  if [[ ${#WATCHDOG_IDS[@]} -gt 0 ]]; then
    watchdog_ids_colon="$(IFS=:; echo "${WATCHDOG_IDS[*]}")"
    watchdog_labels_csv="$(IFS=,; echo "${WATCHDOG_LABELS[*]}")"
    gate_script="$(stage5_source_script src/5_run_benchmark_methods/matrix_gate.sh)"
    stage5_validate_bound_runtime ||
      stage5_abort "Stage 5 bound runtime validation failed before aggregate gate"
    stage5_validate_output_ownership "${PENDING_SELECTION}" ||
      stage5_abort "Stage 5 output ownership changed before aggregate gate"
    stage5_require_source_script "${gate_script}" ||
      stage5_abort "Stage 5 aggregate gate script is outside the immutable source root"
    set +e
    gate_msg="$(sbatch --parsable --wait --dependency="afterany:${watchdog_ids_colon}" --partition="${SLURM_PARTITION_BENCHMARK_CPU}" \
      --ntasks=1 --cpus-per-task=1 --mem=2G --time="${WATCHDOG_TIME_LIMIT}" --output="${RUN_LOG_DIR}/5_matrix_gate_%j.log" \
      --error="${RUN_LOG_DIR}/5_matrix_gate_%j.err" --mail-user="${USER_EMAIL}" \
      "${gate_script}" "${ECODA_RUN_ROOT}" "${watchdog_labels_csv}" "${SCHEDULER_FILE}")"
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
      [[ "${record_kind}" =~ ^(ARRAY|WATCHDOG|STATUS|AGGREGATE_GATE|PREFLIGHT|METADATA_EXPORT)$ &&
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
  target_methods_line="$(sed -n 's/^TARGET_METHODS=//p' \
    "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  methods_line="$(sed -n 's/^METHODS=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  if [[ -n "${target_methods_line}" ]]; then
    ecoda_split_csv "${target_methods_line}"
    METHODS=("${ECODA_ARRAY[@]}")
  elif [[ -n "${methods_line}" ]]; then
    ecoda_split_csv "${methods_line}"
    METHODS=("${ECODA_ARRAY[@]}")
  fi
  labels_line="$(sed -n 's/^ANALYSES=//p' "${ECODA_RUN_ROOT}/metadata" | head -1 || true)"
  if [[ -n "${labels_line}" ]]; then ecoda_split_csv "${labels_line}"; ANALYSES=("${ECODA_ARRAY[@]}"); ANALYSES_SELECTED=1; fi
  [[ -s "${ECODA_RUN_ROOT}/status/aggregate" ]] &&
    grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/aggregate" ||
    stage5_abort "Stage 5 aggregate status is not OK"
fi

stage5_track_pending_artifact_owners ||
  stage5_abort "Stage 5 artifact ownership state is missing or foreign"

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

if [[ -z "${PASS_ARG}" ]]; then
  GLOSCOPE_SHARD_MANIFEST="${ECODA_RUN_ROOT}/manifests/matrix_benchmark_analysis_gloscope_cpu.tsv"
  if [[ -s "${GLOSCOPE_SHARD_MANIFEST}" ]]; then
    ${PIXI_RSCRIPT} "${SCRIPT_DIR}/consolidate_gloscope_results.R" \
      --manifest "${GLOSCOPE_SHARD_MANIFEST}" \
      --results_dir "${ANALYSIS_ROOT}/results" ||
      stage5_abort "Stage 5 GloScope shard consolidation failed"
  fi
fi

stage5_prepare_source_identity ||
  stage5_abort "Stage 5 source identity changed before final validation"
LABELS=()
for method in "${METHODS[@]}"; do
  [[ "${method}" == _ecoda_none_ ]] || LABELS+=("${method}")
done
if [[ ${ANALYSES_SELECTED} -eq 1 ]]; then LABELS+=("${ANALYSES[@]}"); fi
FEATHER_LABELS=()
for label in "${LABELS[@]}"; do
  case "${label}" in
    mrvi|scpoli|pilot|qot|pilotgm)
      FEATHER_LABELS+=("${label}")
      ;;
  esac
done
stage5_validate_final_matrix_rows() {
  local label ds view row_label selection tmp
  local -a validation_args
  [[ "${ANALYSIS_VARIANT:-}" == final ]] || return 1
  [[ ${#FEATHER_LABELS[@]} -gt 0 ]] || return 0
  for label in "${FEATHER_LABELS[@]}"; do
    while IFS=$'\t' read -r ds view row_label; do
      selection="${ECODA_RUN_ROOT}/manifests/matrix_validation_final_${label}.tsv"
      tmp="${selection}.build.$$"
      printf '%s\t%s\t%s\n' "${ds}" "${view}" "${row_label}" > "${tmp}" ||
        return 1
      ecoda_atomic_install_manifest "${tmp}" "${selection}" 3 || {
        rm -f "${tmp}"
        return 1
      }
      rm -f "${tmp}"
      ecoda_write_checksum "${selection}" || return 1
      validation_args=(
        --root "${ANALYSIS_ROOT}"
        --selection "${selection}"
        --labels "${label}"
        --input-root "${HPC_SCRATCH_DIR}"
        --config "${DATASETS_JSON_FILE}"
        --source-identity "${SOURCE_IDENTITY}"
        --source-identity-verified
        --producer "stage5_${label}"
        --producer-run-id "${RUN_ID}"
        --batch
        --batch-pass uncorrected
        --analysis-variant final
      )
      if ! "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" \
          "${validation_args[@]}" >/dev/null 2>&1; then
        rm -f "${selection}" "${selection}.md5"
        return 1
      fi
      rm -f "${selection}" "${selection}.md5"
    done < "${MANIFEST}"
  done
}

if [[ "${ANALYSIS_VARIANT:-}" == final ]]; then
  stage5_validate_final_matrix_rows ||
    stage5_abort "Stage 5 final matrix artifact validation failed"
else
  if [[ ${#FEATHER_LABELS[@]} -gt 0 ]]; then
    for label in "${FEATHER_LABELS[@]}"; do
    if [[ "${PASS_ARG:-}" == corrected ]]; then
      stage5_validate_corrected_matrix_rows "${label}" ||
        stage5_abort "Stage 5 corrected ${label} matrix artifact validation failed"
      continue
    fi
    if stage5_matrix_validation_selection "${label}"; then
      :
    else
      selection_rc=$?
      [[ ${selection_rc} -eq 2 ]] && continue
      stage5_abort "failed to build Stage 5 ${label} validation selection"
    fi
    validation_args=(--root "${ANALYSIS_ROOT}" \
      --selection "${STAGE5_MATRIX_VALIDATION_SELECTION}" --labels "${label}" \
      --input-root "${HPC_SCRATCH_DIR}" --config "${DATASETS_JSON_FILE}" \
      --source-identity "${SOURCE_IDENTITY}" --source-identity-verified \
      --producer "stage5_${label}" --producer-run-id "${RUN_ID}")
    [[ -n "${PASS_ARG}" ]] && validation_args+=(--batch --batch-pass "${PASS_ARG}")
    [[ -n "${ANALYSIS_VARIANT:-}" ]] &&
      validation_args+=(--analysis-variant "${ANALYSIS_VARIANT}")
    [[ ${EXACT_SELECTION} -eq 1 ]] && validation_args+=(--exact)
    if ! "${PYTHON_BIN}" "${SCRIPT_DIR}/matrix_artifact_validator.py" \
        "${validation_args[@]}"; then
      [[ -n "${STAGE5_MATRIX_VALIDATION_SELECTION_TEMP}" ]] &&
        rm -f "${STAGE5_MATRIX_VALIDATION_SELECTION_TEMP}" \
          "${STAGE5_MATRIX_VALIDATION_SELECTION_TEMP}.md5"
      stage5_abort "Stage 5 ${label} matrix artifact validation failed"
    fi
    [[ -n "${STAGE5_MATRIX_VALIDATION_SELECTION_TEMP}" ]] &&
      rm -f "${STAGE5_MATRIX_VALIDATION_SELECTION_TEMP}" \
        "${STAGE5_MATRIX_VALIDATION_SELECTION_TEMP}.md5"
    done
  fi
fi
RDS_LABELS=()
for label in "${LABELS[@]}"; do
  case "${label}" in
    gloscope|mofa|pseudobulk|composition|scitd|prepare_pseudobulk|trans|zeroimp)
      RDS_LABELS+=("${label}")
      ;;
  esac
done

stage5_validate_corrected_rds_rows() {
  local ds view row_label label path identity_path
  local -a corrected_rds_args
  [[ "${PASS_ARG:-}" == corrected ]] || return 1
  [[ ${#RDS_LABELS[@]} -gt 0 ]] || return 0
  while IFS=$'\t' read -r ds view row_label; do
    for label in "${RDS_LABELS[@]}"; do
      if [[ ${METHOD_MATRIX_MODE} -eq 1 ]] &&
         ! ecoda_stage5_method_matrix_allows "${ds}" "${view}" "${label}"; then
        continue
      fi
      identity_path="$(
        stage5_batch_contract_identity_path "${ds}" "${view}" "${label}"
      )" || return 1
      benchmark_artifacts_for "${ds}" "${view}" "${label}" || return 1
      for path in "${ARTIFACT_PATHS[@]}"; do
        corrected_rds_args=(--artifact "${path}" --method "${label}" \
          --dataset "${ds}" --view "${view}" --config "${DATASETS_JSON_FILE}" \
          --input-root "${HPC_SCRATCH_DIR}" --batch-pass corrected \
          --source-identity "${SOURCE_IDENTITY}" --source-identity-verified \
          --expected-batch-contract "${identity_path}")
        [[ -n "${ANALYSIS_VARIANT:-}" ]] &&
          corrected_rds_args+=(--analysis-variant "${ANALYSIS_VARIANT}")
        [[ "${path}" == *_metadata.rds ]] &&
          corrected_rds_args+=(--metadata)
        ${PIXI_RSCRIPT} "${SCRIPT_DIR}/validate_benchmark_rds_contract.R" \
          "${corrected_rds_args[@]}" >/dev/null 2>&1 || return 1
      done
    done
  done < "${MANIFEST}"
}

stage5_validate_final_rds_rows() {
  local label ds view row_label path metadata
  local list tmp
  local -a rds_args
  [[ "${ANALYSIS_VARIANT:-}" == final ]] || return 1
  [[ ${#RDS_LABELS[@]} -gt 0 ]] || return 0
  for label in "${RDS_LABELS[@]}"; do
    list="${ECODA_RUN_ROOT}/manifests/rds_validation_final_${label}.tsv"
    tmp="${list}.build.$$"
    : > "${tmp}" || return 1
    while IFS=$'\t' read -r ds view row_label; do
      benchmark_artifacts_for "${ds}" "${view}" "${label}" || {
        rm -f "${tmp}"
        return 1
      }
      for path in "${ARTIFACT_PATHS[@]}"; do
        metadata=0
        [[ "${path}" == *_metadata.rds ]] && metadata=1
        printf '%s\t%s\t%s\t%s\t%s\n' \
          "${path}" "${label}" "${ds}" "${view}" "${metadata}" >> "${tmp}" ||
          return 1
      done
    done < "${MANIFEST}"
    if [[ ! -s "${tmp}" ]]; then
      rm -f "${tmp}"
      continue
    fi
    ecoda_atomic_install_manifest "${tmp}" "${list}" 5 || {
      rm -f "${tmp}"
      return 1
    }
    rm -f "${tmp}"
    ecoda_write_checksum "${list}" || return 1
    rds_args=(
      --artifact-list "${list}"
      --config "${DATASETS_JSON_FILE}"
      --input-root "${HPC_SCRATCH_DIR}"
      --batch-pass uncorrected
      --source-identity "${SOURCE_IDENTITY}"
      --source-identity-verified
    )
    if ! ${PIXI_RSCRIPT} "${SCRIPT_DIR}/validate_benchmark_rds_contract.R" \
        "${rds_args[@]}" >/dev/null 2>&1; then
      rm -f "${list}" "${list}.md5"
      return 1
    fi
    rm -f "${list}" "${list}.md5"
  done
}

if [[ "${PASS_ARG:-}" == corrected ]]; then
  stage5_validate_corrected_rds_rows ||
    stage5_abort "Stage 5 corrected RDS artifact validation failed"
elif [[ ${#RDS_LABELS[@]} -gt 0 ]]; then
  stage5_validate_pending_rds_records ||
    stage5_abort "Stage 5 RDS artifact record validation failed"
  if [[ "${ANALYSIS_VARIANT:-}" == final ]]; then
    stage5_validate_final_rds_rows ||
      stage5_abort "Stage 5 final RDS artifact validation failed"
  else
    rds_args=(--root "${ANALYSIS_ROOT}" --selection "${MANIFEST}" \
      --labels "$(IFS=,; echo "${RDS_LABELS[*]}")" \
      --config "${DATASETS_JSON_FILE}" --input-root "${HPC_SCRATCH_DIR}" \
      --source-identity "${SOURCE_IDENTITY}" --source-identity-verified)
    [[ -n "${PASS_ARG}" ]] && rds_args+=(--batch-pass "${PASS_ARG}")
    [[ ${EXACT_SELECTION} -eq 1 ]] && rds_args+=(--exact)
    if ! ${PIXI_RSCRIPT} "${SCRIPT_DIR}/validate_benchmark_rds_contract.R" \
        "${rds_args[@]}"; then
      stage5_abort "Stage 5 RDS artifact validation failed"
    fi
  fi
fi
stage5_publish_output_records ||
  stage5_abort "Stage 5 artifact record publication failed"
if ! ( benchmark_merge_sync_cleanup "${LABELS[@]}" ); then
  stage5_abort "Stage 5 benchmark synchronization failed"
fi
stage5_track_pending_artifact_owners ||
  stage5_abort "Stage 5 artifact ownership changed before finalization"
ecoda_owner_finalize_tracked OK "benchmark matrix sync completed" ||
  stage5_abort "failed to finalize global Stage 5 artifact owners"
stage5_finalize_owner_manifest OK "benchmark matrix sync completed" ||
  stage5_abort "failed to finalize Stage 5 owners"
if [[ -n "${PASS_ARG}" ]]; then
  echo "BATCH_EFFECT_RUN_ID=${RUN_ID}"
else
  echo "BENCHMARK_RUN_ID=${RUN_ID}"
fi
