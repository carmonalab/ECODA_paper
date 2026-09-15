#!/bin/bash
#SBATCH --job-name=benchmark_r_worker
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SOURCE_RELATIVE="src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
SOURCE_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
RUNTIME_IMAGE="${ECODA_RUNTIME_IMAGE:-}"
RUNTIME_MANIFEST="${ECODA_RUNTIME_MANIFEST:-}"
RUN_ID="${ECODA_RUN_ID:-}"
RUN_ROOT="${ECODA_RUN_ROOT:-}"
[[ "${SOURCE_REQUIRED}" == "1" && "${SOURCE_ROOT}" = /* &&
   "${SOURCE_MANIFEST}" = /* ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 worker requires an immutable source snapshot." >&2
  exit 1
}
[[ -n "${RUNTIME_IMAGE}" && "${RUNTIME_IMAGE}" = /* &&
   -n "${RUNTIME_MANIFEST}" && "${RUNTIME_MANIFEST}" = /* ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 worker requires a bound runtime image and manifest." >&2
  exit 1
}
[[ -n "${RUN_ID}" && "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ &&
   -n "${RUN_ROOT}" && "${RUN_ROOT}" = /* ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 worker requires a valid bound run root and ID." >&2
  exit 1
}
[[ -d "${SOURCE_ROOT}" && -f "${SOURCE_MANIFEST}" &&
   -r "${SOURCE_MANIFEST}" && ! -L "${SOURCE_MANIFEST}" ]] || {
  echo "ERROR: immutable source manifest is missing or unreadable: ${SOURCE_MANIFEST}" >&2
  exit 1
}
SNAPSHOT_WORKER="${SOURCE_ROOT%/}/${WORKER_SOURCE_RELATIVE}"
[[ -f "${SNAPSHOT_WORKER}" && -r "${SNAPSHOT_WORKER}" && ! -L "${SNAPSHOT_WORKER}" ]] || {
  echo "ERROR: immutable Stage 5 worker is missing or unreadable: ${SNAPSHOT_WORKER}" >&2
  exit 1
}
command -v realpath >/dev/null 2>&1 || {
  echo "ERROR: realpath is required for immutable Stage 5 workers." >&2
  exit 1
}
SNAPSHOT_WORKER_REAL="$(realpath "${SNAPSHOT_WORKER}" 2>/dev/null || true)"
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  command -v scontrol >/dev/null 2>&1 || {
    echo "ERROR: scontrol is required to verify the immutable worker command." >&2
    exit 1
  }
  SCHEDULED_COMMAND="$(scontrol show job "${SLURM_JOB_ID}" -o 2>/dev/null || true)"
  SCHEDULED_SCRIPT="${SCHEDULED_COMMAND#*Command=}"
  SCHEDULED_SCRIPT="${SCHEDULED_SCRIPT%% *}"
  SCHEDULED_SCRIPT_REAL="$(realpath "${SCHEDULED_SCRIPT}" 2>/dev/null || true)"
  [[ -n "${SCHEDULED_SCRIPT_REAL}" && "${SCHEDULED_SCRIPT_REAL}" == "${SNAPSHOT_WORKER_REAL}" ]] || {
    echo "ERROR: mutable or non-snapshot Stage 5 worker command rejected." >&2
    exit 1
  }
else
  CURRENT_SCRIPT_REAL="$(realpath "${BASH_SOURCE[0]}" 2>/dev/null || true)"
  [[ -n "${CURRENT_SCRIPT_REAL}" && "${CURRENT_SCRIPT_REAL}" == "${SNAPSHOT_WORKER_REAL}" ]] || {
    echo "ERROR: mutable or non-snapshot Stage 5 worker path rejected." >&2
    exit 1
  }
fi
source "${SOURCE_ROOT%/}/src/slurm_config.sh"
source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_runtime.sh"
PREEXEC_ROOT_VERSION="${ANALYSIS_ROOT_VERSION:-}"
PREEXEC_ROOT_IDENTITY="${ANALYSIS_ROOT_IDENTITY:-}"
PREEXEC_ROOT_ALIAS="${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION:-}"
if [[ -n "${PREEXEC_ROOT_ALIAS}" &&
      -n "${PREEXEC_ROOT_VERSION}" &&
      "${PREEXEC_ROOT_ALIAS}" != "${PREEXEC_ROOT_VERSION}" ]]; then
  echo "ERROR: Stage 5 root-version aliases disagree before runtime reexec." >&2
  exit 1
fi
if [[ -n "${PREEXEC_ROOT_ALIAS}" ]]; then
  [[ "${PREEXEC_ROOT_ALIAS}" == "recovery_35row" ]] || {
    echo "ERROR: unsupported corrected-final root version before runtime reexec: ${PREEXEC_ROOT_ALIAS}" >&2
    exit 1
  }
  PREEXEC_ROOT_VERSION="${PREEXEC_ROOT_ALIAS}"
fi
[[ -z "${PREEXEC_ROOT_VERSION}" ||
   "${PREEXEC_ROOT_VERSION}" == "recovery_35row" ]] || {
  echo "ERROR: unsupported ANALYSIS_ROOT_VERSION before runtime reexec: ${PREEXEC_ROOT_VERSION}" >&2
  exit 1
}
if [[ "${PREEXEC_ROOT_VERSION}" == "recovery_35row" ]]; then
  [[ -z "${PREEXEC_ROOT_IDENTITY}" ||
     "${PREEXEC_ROOT_IDENTITY}" == "corrected_final/recovery_35row" ]] || {
    echo "ERROR: replacement root identity disagrees before runtime reexec." >&2
    exit 1
  }
  PREEXEC_ROOT_IDENTITY="corrected_final/recovery_35row"
fi
ANALYSIS_ROOT_VERSION="${PREEXEC_ROOT_VERSION}"
ANALYSIS_ROOT_IDENTITY="${PREEXEC_ROOT_IDENTITY}"
export ANALYSIS_VARIANT="${ANALYSIS_VARIANT:-}" \
  ANALYSIS_PASS="${ANALYSIS_PASS:-}" \
  ANALYSIS_ROOT="${ANALYSIS_ROOT:-}" \
  ANALYSIS_NAS_ROOT="${ANALYSIS_NAS_ROOT:-}" \
  ANALYSIS_LOG_PREFIX="${ANALYSIS_LOG_PREFIX:-}"
if [[ -n "${ANALYSIS_ROOT_VERSION}" ]]; then
  export ANALYSIS_ROOT_VERSION
else
  unset ANALYSIS_ROOT_VERSION
fi
if [[ -n "${ANALYSIS_ROOT_IDENTITY}" ]]; then
  export ANALYSIS_ROOT_IDENTITY
else
  unset ANALYSIS_ROOT_IDENTITY
fi
PREEXEC_MATRIX="${ECODA_STAGE5_METHOD_MATRIX:-${METHOD_MATRIX:-}}"
if [[ -n "${ECODA_STAGE5_METHOD_MATRIX:-}" &&
      -n "${METHOD_MATRIX:-}" &&
      "${ECODA_STAGE5_METHOD_MATRIX}" != "${METHOD_MATRIX}" ]]; then
  echo "ERROR: Stage 5 method-matrix aliases disagree before runtime reexec." >&2
  exit 1
fi
if [[ -n "${PREEXEC_MATRIX}" ]]; then
  export METHOD_MATRIX="${PREEXEC_MATRIX}"
fi
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export ECODA_RUN_ID="${RUN_ID}"
export ECODA_RUN_ROOT="${RUN_ROOT}"
EXPECTED_RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
[[ "${RUN_ROOT}" == "${EXPECTED_RUN_ROOT}" && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 5 worker run root is not the exact bound run root: ${RUN_ROOT}" >&2
  exit 1
}
RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"
RUN_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
[[ -f "${RUN_SOURCE_MANIFEST}" && -r "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" &&
   -f "${RUN_RUNTIME_IDENTITY}" && -r "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 run-bound source/runtime manifests are missing." >&2
  exit 1
}
cmp -s "${SOURCE_MANIFEST}" "${RUN_SOURCE_MANIFEST}" || {
  echo "ERROR: Stage 5 run source manifest differs from the immutable source manifest." >&2
  exit 1
}
ecoda_runtime_reexec_worker stage5 \
  "${SNAPSHOT_WORKER}" || exit 1
SCRIPT_DIR="${SOURCE_ROOT%/}/src/5_run_benchmark_methods/run_r_sample_embedding_methods"
source "${SOURCE_ROOT%/}/src/slurm_config.sh"
source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_runtime.sh"
source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_run_common.sh"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUN_ID="${RUN_ID}"
export ECODA_RUN_ROOT="${RUN_ROOT}"
if [[ "${ECODA_RUNTIME_MODE:-host}" == "host" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  ecoda_runtime_validate_bound_run || {
    echo "ERROR: Stage 5 run-bound runtime validation failed." >&2
    exit 1
  }
fi
ecoda_validate_run_id "${ECODA_RUN_ID}" || exit 1
[[ "${ECODA_RUN_ROOT}" == "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 5 run root changed during immutable worker reexec." >&2
  exit 1
}
cd "${PROJECT_ROOT}"
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"
[[ -n "${METHOD:-}" ]] || { echo "ERROR: METHOD is not set." >&2; exit 1; }
case "${METHOD}" in
  gloscope|mofa|pseudobulk|scitd|composition|trans|zeroimp)
    export ECODA_ARTIFACT_PRODUCER="stage5_${METHOD}"
    ;;
esac
MANIFEST_PATH="${MATRIX_RETRY_MANIFEST:-${ANALYSIS_MANIFEST:-${BENCHMARK_MANIFEST:-}}}"
[[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: benchmark manifest is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
ecoda_validate_run_owned_path "${MANIFEST_PATH}" "${ECODA_RUN_ROOT}" || {
  echo "ERROR: benchmark manifest is outside the bound Stage 5 run root: ${MANIFEST_PATH}" >&2
  exit 1
}
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
[[ -n "${line}" ]] || { echo "ERROR: no benchmark row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME ROW_VIEW ROW_LABEL ROW_COMBO ROW_EXTRA <<< "${line}"
[[ -n "${DS_NAME}" ]] || { echo "ERROR: malformed benchmark row" >&2; exit 1; }
[[ -z "${ROW_EXTRA:-}" ]] || {
  echo "ERROR: benchmark row has more than four tab-separated fields" >&2
  exit 1
}
if [[ -n "${ROW_COMBO:-}" ]]; then
  [[ "${METHOD}" == gloscope ]] || {
    echo "ERROR: combo token is only supported for method gloscope" >&2
    exit 1
  }
  [[ -z "${ANALYSIS_PASS:-}" ]] || {
    echo "ERROR: GloScope combo shards are only supported for ordinary runs" >&2
    exit 1
  }
  case "${ROW_COMBO}" in
    hvg2000_pcadims10|hvg2000_pcadims30|hvg2000_pcadims50|\
      hvg1000_pcadims30|hvg3000_pcadims30) ;;
    *)
      echo "ERROR: invalid GloScope combo token: ${ROW_COMBO}" >&2
      exit 1
      ;;
  esac
fi
export DS_NAME
ANALYSIS_VIEW="${ROW_VIEW:-${ANALYSIS_VIEW:-benchmark_analysis}}"
export ANALYSIS_VIEW
ANALYSIS_VARIANT="${ANALYSIS_VARIANT:-}"
ANALYSIS_PASS="${ANALYSIS_PASS:-}"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
ANALYSIS_NAS_ROOT="${ANALYSIS_NAS_ROOT:-}"
ANALYSIS_LOG_PREFIX="${ANALYSIS_LOG_PREFIX:-}"
ANALYSIS_ROOT_VERSION="${ANALYSIS_ROOT_VERSION:-}"
ANALYSIS_ROOT_IDENTITY="${ANALYSIS_ROOT_IDENTITY:-}"
ROOT_VERSION_ALIAS="${ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION:-}"
if [[ -n "${ROOT_VERSION_ALIAS}" ]]; then
  [[ "${ROOT_VERSION_ALIAS}" == "recovery_35row" ]] || {
    echo "ERROR: unsupported corrected-final root version: ${ROOT_VERSION_ALIAS}" >&2
    exit 1
  }
  if [[ -n "${ANALYSIS_ROOT_VERSION}" &&
        "${ANALYSIS_ROOT_VERSION}" != "${ROOT_VERSION_ALIAS}" ]]; then
    echo "ERROR: Stage 5 root-version aliases disagree." >&2
    exit 1
  fi
  ANALYSIS_ROOT_VERSION="${ROOT_VERSION_ALIAS}"
fi
[[ -z "${ANALYSIS_ROOT_VERSION}" ||
   "${ANALYSIS_ROOT_VERSION}" == "recovery_35row" ]] || {
  echo "ERROR: unsupported Stage 5 ANALYSIS_ROOT_VERSION: ${ANALYSIS_ROOT_VERSION}" >&2
  exit 1
}
EXPECTED_VARIANT_PASS=""
EXPECTED_ANALYSIS_ROOT=""
EXPECTED_ANALYSIS_NAS_ROOT=""
EXPECTED_LOG_PREFIX=""
case "${ANALYSIS_VARIANT}" in
  "")
    [[ -z "${ANALYSIS_ROOT_VERSION}" &&
       -z "${ANALYSIS_ROOT_IDENTITY}" &&
       -z "${ROOT_VERSION_ALIAS}" ]] || {
      echo "ERROR: root-version identity requires an explicit Stage 5 variant." >&2
      exit 1
    }
    ;;
  final)
    EXPECTED_VARIANT_PASS="uncorrected"
    EXPECTED_ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final"
    EXPECTED_ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/uncorrected_final"
    EXPECTED_LOG_PREFIX="execution_times_batch_effect_uncorrected_final_"
    [[ -z "${ANALYSIS_ROOT_VERSION}" &&
       ( -z "${ANALYSIS_ROOT_IDENTITY}" ||
         "${ANALYSIS_ROOT_IDENTITY}" == "uncorrected_final" ) ]] || {
      echo "ERROR: final Stage 5 workers cannot use a corrected-final root identity." >&2
      exit 1
    }
    ;;
  corrected_final)
    EXPECTED_VARIANT_PASS="corrected"
    EXPECTED_LOG_PREFIX="execution_times_batch_effect_corrected_final_"
    if [[ -n "${ANALYSIS_ROOT_VERSION}" ||
          "${ANALYSIS_ROOT_IDENTITY}" == "corrected_final/recovery_35row" ]]; then
      ANALYSIS_ROOT_VERSION="${ANALYSIS_ROOT_VERSION:-recovery_35row}"
      [[ "${ANALYSIS_ROOT_VERSION}" == "recovery_35row" ]] || {
        echo "ERROR: replacement corrected-final root requires recovery_35row." >&2
        exit 1
      }
      [[ -z "${ANALYSIS_ROOT_IDENTITY}" ||
         "${ANALYSIS_ROOT_IDENTITY}" == "corrected_final/recovery_35row" ]] || {
        echo "ERROR: corrected-final root identity disagrees with its version." >&2
        exit 1
      }
      ANALYSIS_ROOT_IDENTITY="corrected_final/recovery_35row"
      EXPECTED_ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/corrected_final/recovery_35row"
      EXPECTED_ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/corrected_final/recovery_35row"
    else
      [[ -z "${ANALYSIS_ROOT_IDENTITY}" ||
         "${ANALYSIS_ROOT_IDENTITY}" == "corrected_final" ]] || {
        echo "ERROR: invalid direct corrected-final root identity." >&2
        exit 1
      }
      EXPECTED_ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/corrected_final"
      EXPECTED_ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/corrected_final"
      case "${ANALYSIS_ROOT}:${ANALYSIS_NAS_ROOT}" in
        */batch_effect/corrected_final/recovery_35row:*)
          echo "ERROR: replacement corrected-final root is missing its bound identity." >&2
          exit 1
          ;;
        *:*/batch_effect/corrected_final/recovery_35row)
          echo "ERROR: replacement corrected-final NAS root is missing its bound identity." >&2
          exit 1
          ;;
      esac
      ANALYSIS_ROOT_IDENTITY="corrected_final"
    fi
    ;;
  *)
    echo "ERROR: unsupported Stage 5 R analysis variant: ${ANALYSIS_VARIANT}" >&2
    exit 1
    ;;
esac
if [[ -n "${ANALYSIS_VARIANT}" ]]; then
  [[ "${ANALYSIS_PASS}" == "${EXPECTED_VARIANT_PASS}" &&
     "${ANALYSIS_ROOT}" == "${EXPECTED_ANALYSIS_ROOT}" ]] || {
    echo "ERROR: Stage 5 R worker variant/pass/root identity mismatch." >&2
    exit 1
  }
  if [[ -n "${ANALYSIS_NAS_ROOT}" &&
        "${ANALYSIS_NAS_ROOT}" != "${EXPECTED_ANALYSIS_NAS_ROOT}" ]]; then
    echo "ERROR: Stage 5 R worker NAS root identity mismatch." >&2
    exit 1
  fi
  if [[ -n "${ANALYSIS_LOG_PREFIX}" &&
        "${ANALYSIS_LOG_PREFIX}" != "${EXPECTED_LOG_PREFIX}" ]]; then
    echo "ERROR: Stage 5 R worker log identity mismatch." >&2
    exit 1
  fi
  ANALYSIS_LOG_PREFIX="${EXPECTED_LOG_PREFIX}"
fi
if [[ -n "${ANALYSIS_ROOT_VERSION}" ]]; then
  ROOT_VERSION_ALIAS="${ANALYSIS_ROOT_VERSION}"
  export ANALYSIS_ROOT_VERSION
  export ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION="${ROOT_VERSION_ALIAS}"
else
  unset ANALYSIS_ROOT_VERSION ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION
fi
if [[ -n "${ANALYSIS_ROOT_IDENTITY}" ]]; then
  export ANALYSIS_ROOT_IDENTITY
else
  unset ANALYSIS_ROOT_IDENTITY
fi
validate_method_matrix_environment() {
  local matrix="${ECODA_STAGE5_METHOD_MATRIX:-${METHOD_MATRIX:-}}"
  local legacy_matrix="${METHOD_MATRIX:-}"
  local expected_md5="${METHOD_MATRIX_MD5:-}"
  local expected_size="${METHOD_MATRIX_SIZE:-}"
  local expected_sha256="${METHOD_MATRIX_SHA256:-}"
  local expected_identity="${METHOD_MATRIX_IDENTITY:-}"
  local expected_count="${METHOD_MATRIX_DECLARED_COUNT:-}"
  local actual_sha256 actual_count
  if [[ -n "${ECODA_STAGE5_METHOD_MATRIX:-}" &&
        -n "${legacy_matrix}" &&
        "${ECODA_STAGE5_METHOD_MATRIX}" != "${legacy_matrix}" ]]; then
    echo "ERROR: Stage 5 method-matrix aliases disagree." >&2
    return 1
  fi
  [[ -z "${matrix}" ]] && return 0
  [[ "${matrix}" = /* && "${matrix}" != *$'\n'* &&
     "${matrix}" != *$'\r'* && "${matrix}" != *$'\t'* &&
     -f "${matrix}" && ! -L "${matrix}" && -r "${matrix}" ]] || {
    echo "ERROR: Stage 5 method matrix is missing or unsafe: ${matrix}" >&2
    return 1
  }
  ecoda_validate_run_owned_path "${matrix}" "${ECODA_RUN_ROOT}" || {
    echo "ERROR: Stage 5 method matrix is outside the bound run root: ${matrix}" >&2
    return 1
  }
  ecoda_validate_manifest "${matrix}" 3 || {
    echo "ERROR: Stage 5 method matrix is not a valid three-column TSV." >&2
    return 1
  }
  ecoda_validate_checksum "${matrix}" || {
    echo "ERROR: Stage 5 method matrix checksum is invalid: ${matrix}" >&2
    return 1
  }
  [[ -z "${expected_md5}" || "${expected_md5}" == "${ECODA_CHECKSUM_MD5}" ]] || return 1
  [[ -z "${expected_size}" || "${expected_size}" == "${ECODA_CHECKSUM_SIZE}" ]] || return 1
  actual_sha256="$(ecoda_sha256_file "${matrix}")" || return 1
  [[ -z "${expected_sha256}" || "${expected_sha256}" == "${actual_sha256}" ]] || return 1
  [[ -z "${expected_identity}" || "${expected_identity}" == "${actual_sha256}" ]] || return 1
  actual_count="$(wc -l < "${matrix}" | tr -d '[:space:]')" || return 1
  [[ -z "${expected_count}" || "${expected_count}" == "${actual_count}" ]] || return 1
  export METHOD_MATRIX="${matrix}" ECODA_STAGE5_METHOD_MATRIX="${matrix}"
  export ECODA_STAGE5_METHOD_MATRIX_MODE=1
}
validate_method_matrix_environment || exit 1
export ANALYSIS_VARIANT ANALYSIS_PASS ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_LOG_PREFIX
PASS="${ANALYSIS_PASS}"
ROOT="${ANALYSIS_ROOT}"
export PASS ROOT
RESULTS_DIR="${ANALYSIS_ROOT}/results"
PSEUDOBULK_DIR="${ANALYSIS_ROOT}/pseudobulks"
GLOSCOPE_DIR="${ANALYSIS_ROOT}/gloscope_dists"
EMBED_DIR="${ANALYSIS_ROOT}/embeddings"
EXECUTION_LOG_DIR="${EXECUTION_LOG_DIR:-${EMBED_DIR}}"
if [[ -n "${ROW_COMBO:-}" ]]; then
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_${METHOD}_${DS_NAME}_${ROW_COMBO}.feather"
elif [[ -n "${ANALYSIS_VARIANT}" ]]; then
  LOG_FILE="${EXECUTION_LOG_DIR}/${ANALYSIS_LOG_PREFIX}${METHOD}_${DS_NAME}.feather"
elif [[ -n "${ANALYSIS_PASS}" ]]; then
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_batch_effect_${ANALYSIS_PASS}_${METHOD}_${DS_NAME}.feather"
else
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_${METHOD}_${DS_NAME}.feather"
fi
FORCE_FLAG=()
[[ "${FORCE_BENCHMARK:-0}" == 1 ]] && FORCE_FLAG=(--force)
ANALYSIS_PASS_FLAG=()
[[ -n "${ANALYSIS_PASS:-}" ]] && ANALYSIS_PASS_FLAG=(--analysis_pass "${ANALYSIS_PASS}")
COMBO_FLAG=()
[[ -n "${ROW_COMBO:-}" ]] && COMBO_FLAG=(--combo "${ROW_COMBO}")
if [[ "${METHOD}" == prepare_pseudobulk ]]; then
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_prepare_pseudobulk.R"
else
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_run_benchmark_methods_r.R"
fi
printf 'ECODA_WORKER_DISPATCH METHOD=%s R_SCRIPT=%s ANALYSIS_PASS=%s FORCE_BENCHMARK=%s\n' \
  "${METHOD}" "${R_SCRIPT}" "${ANALYSIS_PASS:-}" "${FORCE_BENCHMARK:-0}"

set +e
R_ARGS=(
  "${R_SCRIPT}"
  --config_path "${DATASETS_JSON_FILE}" --ds_name "${DS_NAME}" --view "${ANALYSIS_VIEW}"
  --method "${METHOD}" --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output"
  --results_dir "${RESULTS_DIR}" --pseudobulk_dir "${PSEUDOBULK_DIR}"
  --gloscope_cache_dir "${GLOSCOPE_DIR}" --log_file "${LOG_FILE}"
)
if [[ ${#FORCE_FLAG[@]} -gt 0 ]]; then R_ARGS+=("${FORCE_FLAG[@]}"); fi
if [[ ${#ANALYSIS_PASS_FLAG[@]} -gt 0 ]]; then R_ARGS+=("${ANALYSIS_PASS_FLAG[@]}"); fi
if [[ ${#COMBO_FLAG[@]} -gt 0 ]]; then R_ARGS+=("${COMBO_FLAG[@]}"); fi
${PIXI_RSCRIPT} "${R_ARGS[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then worker_clear_retry_count; exit 0; fi
ERR_PREFIX="${JOB_LOG_PREFIX:-5_benchmark_r_${METHOD}}"
ERR_FILE="${LOGS_DIR}/${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then exit 0; fi
exit ${RC}
