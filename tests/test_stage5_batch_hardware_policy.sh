#!/bin/bash
# Actual dispatcher resource/submission functions and watchdog retry command;
# no Slurm, runtime, scientific artifacts or HPC access.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-batch-hardware.XXXXXX")"
trap 'rm -rf "${TMP}"' EXIT
DISPATCHER="${ROOT}/src/5_run_benchmark_methods/stage5_dispatcher.sh"
WATCHDOG="${ROOT}/src/5_run_benchmark_methods/matrix_watchdog.sh"
for name in gpu_method_is_default method_spec method_spec_for_group; do
  eval "$(awk -v name="${name}" '$0 == name "() {" {copy=1} copy {print} copy && /^}$/ {exit}' "${DISPATCHER}")"
done
eval "$(awk '/^  submit_matrix\(\) \{/ {copy=1} copy {print} copy && /^  }$/ {exit}' "${DISPATCHER}")"
RETRY_COMMAND="$(awk '/^  retry_msg="\$\(sbatch / {copy=1} copy {print} copy && /^  retry_rc=\$\?$/ {exit}' "${WATCHDOG}")"
[[ -n "${RETRY_COMMAND}" ]]
PROJECT_ROOT="${ROOT}"
GPU_POLICY=auto; PASS_ARG=""; PARTITION_ARG=""
EXCLUDED_NODES="${TEST_EXCLUDED_NODES:-}"; EXCLUDE_FLAG=""
[[ -z "${EXCLUDED_NODES}" ]] || EXCLUDE_FLAG="--exclude=${EXCLUDED_NODES}"
BENCHMARK_CPU_CONSTRAINT=EPYC-7742; BENCHMARK_GPU_CONSTRAINT=nvidia_h200_nvl
BENCHMARK_CPU_CPUS_PER_TASK=16; BENCHMARK_GPU_CPUS_PER_TASK=8; BENCHMARK_GPU_COUNT=1
BENCHMARK_GPU_DEFAULT_METHODS=mrvi,scpoli
SLURM_PARTITION_BENCHMARK_CPU=shared-cpu; MAX_NUM_CHUNKS_PARALLEL=12
BENCHMARK_GPU_DEFAULT_PARTITION=shared-gpu; BENCHMARK_GPU_ANY_PARTITION=shared-gpu,private-gpu
BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE=4; BENCHMARK_GPU_ANY_ARRAY_THROTTLE=3
BENCHMARK_CPU_TIME_LIMIT=72:00:00; BENCHMARK_GPU_DEFAULT_TIME_LIMIT=12:00:00; BENCHMARK_GPU_ANY_TIME_LIMIT=18:00:00
BENCHMARK_GPU_ANY_VRAM_PER_GPU=16G
MEMORY=4G; MAX_MEMORY=8G; WATCHDOG_TIME_LIMIT=01:00:00; USER_EMAIL=test@example.invalid
ECODA_RECOVERY_ACTIVE=0; METHOD_MATRIX_MODE=0; ANALYSIS_ROOT="${TMP}/analysis"
ANALYSIS_VARIANT=""; ECODA_BATCH_CONTRACT_MANIFEST="${TMP}/batch-contract.json"
ARRAY_IDS=(); WATCHDOG_IDS=(); WATCHDOG_LABELS=()
stage5_method_is_forced() { return 1; }
ecoda_stage5_runtime_export_csv() { printf 'RUNTIME_NV=%s' "$2"; }
stage5_source_script() { printf '%s/%s' "${ROOT}" "$1"; }
stage5_validate_bound_runtime() { :; }
stage5_validate_output_ownership() { :; }
stage5_require_source_script() { [[ -f "$1" ]]; }
stage5_record_scheduler() { printf '%s\t%s\n' "$1" "$2" >> "${CASE_ROOT}/ids"; }
sbatch() {
  local arg seen_worker=0 position=0
  printf '%s\n' "$*" >> "${CASE_ROOT}/calls"
  for arg in "$@"; do
    if [[ ${seen_worker} -eq 1 ]]; then
      position=$((position + 1))
      [[ ${position} -le 10 ]] || printf '%s\n' "${arg}" >> "${CASE_ROOT}/retry-flags"
    elif [[ "${arg}" == */matrix_watchdog.sh ]]; then seen_worker=1; fi
  done
  printf '900%s\n' "$(wc -l < "${CASE_ROOT}/calls" | tr -d '[:space:]')"
}
exercise() {
  local view="$1" method="$2" class="$3" flag expected_flags=() actual_flags=()
  CASE_ROOT="${TMP}/${view}_${method}_${class}"
  mkdir -p "${CASE_ROOT}/manifests" "${CASE_ROOT}/logs"
  ECODA_RUN_ROOT="${CASE_ROOT}"; RUN_ID="${view}_${method}_${class}"; RUN_LOG_DIR="${CASE_ROOT}/logs"
  case "${view}" in
    benchmark_analysis) PASS_ARG="" ;;
    batch_effect_uncorrected) PASS_ARG=uncorrected ;;
    batch_effect_corrected) PASS_ARG=corrected ;;
  esac
  # Same pre-existing policy/context; only the actual view controls model pins.
  method_spec_for_group "${method}" "${class}" benchmark_analysis
  local expected_partition="${METHOD_PARTITION}" expected_time="${METHOD_TIME_LIMIT}" expected_nv="${METHOD_RUNTIME_NV}" expected_policy="${METHOD_GPU_POLICY}"
  for flag in "${METHOD_FLAGS[@]}"; do
    case "${view}:${flag}" in batch_effect_*:--constraint=*) ;; *) expected_flags+=("${flag}");; esac
  done
  method_spec_for_group "${method}" "${class}" "${view}"
  [[ "${METHOD_FLAGS[*]}" == "${expected_flags[*]}" && "${METHOD_PARTITION}" == "${expected_partition}" &&
     "${METHOD_TIME_LIMIT}" == "${expected_time}" && "${METHOD_RUNTIME_NV}" == "${expected_nv}" && "${METHOD_GPU_POLICY}" == "${expected_policy}" ]]
  local manifest="${CASE_ROOT}/manifests/selection.tsv"
  printf 'D\t%s\t%s\n' "${view}" "${method}" > "${manifest}"
  submit_matrix "${view}__${method}" "${manifest}" "" "${method}" "${view}" "${class}" >/dev/null
  [[ "$(wc -l < "${CASE_ROOT}/calls" | tr -d '[:space:]')" == 2 ]]
  while IFS= read -r flag; do actual_flags+=("${flag}"); done < "${CASE_ROOT}/retry-flags"
  [[ "${actual_flags[*]}" == "${expected_flags[*]}" ]]
  local initial
  initial="$(sed -n '1p' "${CASE_ROOT}/calls")"
  case "${initial}" in *"--mem=4G"*"RUNTIME_NV=${expected_nv}"*) ;; *) exit 1;; esac
  # Execute the watchdog's exact OOM retry submission with its captured flags.
  WORKER_FLAGS=("${actual_flags[@]}"); PARTITION="${expected_partition}"; THROTTLE="${METHOD_THROTTLE}"
  WORKER_SCRIPT="${METHOD_WORKER}"; WORKER_TIME_LIMIT="${expected_time}"; NEXT_MEMORY=8G; retry_count=1
  RETRY_LOG_PREFIX="${CASE_ROOT}/logs/retry"; retry_export="ALL,RUNTIME_NV=${expected_nv}"
  eval "${RETRY_COMMAND}"
  [[ ${retry_rc} -eq 0 && "$(wc -l < "${CASE_ROOT}/calls" | tr -d '[:space:]')" == 3 ]]
  local retry
  retry="$(sed -n '3p' "${CASE_ROOT}/calls")"
  case "${retry}" in *"--mem=8G"*"RUNTIME_NV=${expected_nv}"*) ;; *) exit 1;; esac
  for flag in "${expected_flags[@]}"; do
    case "${initial}" in *"${flag}"*) ;; *) exit 1;; esac
    case "${retry}" in *"${flag}"*) ;; *) exit 1;; esac
  done
  if [[ "${view}" == batch_effect_* ]]; then
    if rg -q -- '--constraint(=| )' "${CASE_ROOT}/calls"; then echo 'batch model constraint survived initial/watchdog/retry' >&2; exit 1; fi
  fi
  if [[ -n "${EXCLUDED_NODES}" ]]; then
    [[ $(grep -c -- "--exclude=${EXCLUDED_NODES}" "${CASE_ROOT}/calls") == 3 ]]
  elif grep -q -- '--exclude=' "${CASE_ROOT}/calls"; then
    echo 'default submissions unexpectedly excluded nodes' >&2; exit 1
  fi
}
for view in benchmark_analysis batch_effect_uncorrected batch_effect_corrected; do
  for method in prepare_pseudobulk pseudobulk composition gloscope mofa scitd pilot pilotgm qot trans zeroimp; do exercise "${view}" "${method}" base; done
  exercise "${view}" mrvi base
  exercise "${view}" mrvi cpu
  exercise "${view}" scpoli default_gpu
  exercise "${view}" scpoli any_gpu
done
rg -q -- '--constraint=EPYC-7742' "${TMP}/benchmark_analysis_composition_base/calls"
rg -q -- '--constraint=nvidia_h200_nvl' "${TMP}/benchmark_analysis_mrvi_base/calls"
rg -q -- '--gpus=1' "${TMP}/batch_effect_corrected_scpoli_default_gpu/calls"
rg -q -- '--gres=gpu:1,VramPerGpu:16G' "${TMP}/batch_effect_uncorrected_scpoli_any_gpu/calls"
if method_spec_for_group composition base foreign_view >/dev/null 2>&1; then exit 1; fi
echo 'PASS: benchmark model pins preserved; both batch views model-free with identical allocations and actual initial/watchdog/retry flags'
