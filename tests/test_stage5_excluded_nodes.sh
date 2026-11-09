#!/bin/bash
# Native Slurm arguments only; no scheduler, H5AD, fits or environment mutation.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(cd "$(mktemp -d "${TMPDIR:-/tmp}/ecoda-excluded-nodes.XXXXXX")" && pwd -P)"
trap 'rm -rf "${TMP}"' EXIT
HPC_SCRATCH_DIR="${TMP}/scratch"; mkdir -p "${HPC_SCRATCH_DIR}"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/utils/bash/h5ad_preflight_submit.sh"
NODES=cpu046,cpu051,cpu052
for value in cpu046 "${NODES}"; do ecoda_validate_excluded_nodes "${value}"; done
for value in '' cpu046,cpu046 'cpu[046-052]' 'cpu046, cpu051' '--help' 'cpu046,' 'cpu046;false'; do
  if ecoda_validate_excluded_nodes "${value}"; then echo 'unsafe excluded-node list accepted' >&2; exit 1; fi
done
DISPATCHER="${ROOT}/src/5_run_benchmark_methods/stage5_dispatcher.sh"
for name in stage5_run_r_environment_preflight stage5_validate_bound_runtime; do
  eval "$(awk -v name="${name}" '$0 == name "() {" {copy=1} copy {print} copy && /^}$/ {exit}' "${DISPATCHER}")"
done
METADATA_COMMAND="$(awk '/^    export_msg="\$\(sbatch / {copy=1} copy {print} copy && /^    export_rc=\$\?$/ {exit}' "${DISPATCHER}")"
AGGREGATE_COMMAND="$(awk '/^    gate_msg="\$\(sbatch / {copy=1} copy {print} copy && /^    gate_rc=\$\?$/ {exit}' "${DISPATCHER}")"
[[ -n "${METADATA_COMMAND}" && -n "${AGGREGATE_COMMAND}" ]]
# Use real run-owned path/checksum validation; isolate runtime and scheduler.
ecoda_runtime_validate_bound_run() { :; }
ecoda_runtime_validate_submission() { :; }
ecoda_require_source_script_path() { printf '%s' "$1"; }
stage5_source_script() { printf '%s/%s' "${ROOT}" "$1"; }
stage5_require_source_script() { [[ -f "$1" ]]; }
stage5_record_scheduler() { :; }
sbatch() { printf '%s\n' "$@" >> "${CAPTURE}"; printf '90001\n'; }
USER_EMAIL=test@example.invalid; MEMORY=4G; THROTTLE=3; WATCHDOG_TIME_LIMIT=00:10:00
SLURM_PARTITION_BENCHMARK_CPU=shared-cpu; RUNTIME_EXPORT=RUNTIME_NV=0; ECODA_RECOVERY_ACTIVE=0
ECODA_SOURCE_ROOT="${ROOT}"; ECODA_SOURCE_MANIFEST="${TMP}/source.manifest"; ECODA_SOURCE_SNAPSHOT_REQUIRED=1
ECODA_HOST_ENV_PREFIX="${TMP}/host"; ECODA_RUNTIME_IMAGE="${TMP}/runtime.sif"; ECODA_RUNTIME_MANIFEST="${TMP}/runtime.manifest"
ECODA_AUX_ROOT="${TMP}/aux"; ECODA_LOGS_DIR="${TMP}/logs"
for excluded in '' "${NODES}"; do
  EXCLUDED_NODES="${excluded}"; EXCLUDE_FLAG=""
  [[ -z "${excluded}" ]] || EXCLUDE_FLAG="--exclude=${excluded}"
  RUN_ID="case_${excluded:+excluded}"; ECODA_RUN_ID="${RUN_ID}"; ECODA_RUN_ROOT="${ECODA_RUNS_ROOT}/${RUN_ID}"
  mkdir -p "${ECODA_RUN_ROOT}/manifests" "${ECODA_RUN_ROOT}/logs"
  if [[ -n "${excluded}" ]]; then
    printf '%s\n' "${excluded}" > "${ECODA_RUN_ROOT}/manifests/excluded_nodes.tsv"
    ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/excluded_nodes.tsv"
  fi
  manifest="${ECODA_RUN_ROOT}/manifests/input.tsv"; printf 'D\tbenchmark_analysis\t/input.h5ad\n' > "${manifest}"
  CAPTURE="${ECODA_RUN_ROOT}/calls"; : > "${CAPTURE}"
  ecoda_submit_h5ad_preflight "${manifest}" "${ECODA_RUN_ROOT}/status/preflight" "${ECODA_RUN_ROOT}" require shared-cpu 4G 3 \
    "${ECODA_RUN_ROOT}/logs" test "${ROOT}/src/utils/bash/h5ad_preflight_worker.sh" "${RUNTIME_EXPORT}" "${excluded}" >/dev/null
  stage5_run_r_environment_preflight >/dev/null
  pending_count=1; PARTITION_ARG=""; pending_manifest="${manifest}"; status_dir="${ECODA_RUN_ROOT}/status/metadata"
  export_runtime=RUNTIME_NV=0; worker="${ROOT}/src/utils/bash/h5ad_obs_audit_worker.sh"
  eval "${METADATA_COMMAND}"; [[ ${export_rc} == 0 ]]
  watchdog_ids_colon=90001; watchdog_labels_csv=test; RUN_LOG_DIR="${ECODA_RUN_ROOT}/logs"
  gate_script="${ROOT}/src/5_run_benchmark_methods/matrix_gate.sh"; SCHEDULER_FILE="${manifest}"
  eval "${AGGREGATE_COMMAND}"; [[ ${gate_rc} == 0 ]]
  [[ $(grep -c '^--parsable$' "${CAPTURE}") == 4 ]]
  if [[ -n "${excluded}" ]]; then
    [[ $(grep -cx -- "--exclude=${excluded}" "${CAPTURE}") == 4 ]]
    stage5_validate_bound_runtime
    printf 'cpu999\n' > "${ECODA_RUN_ROOT}/manifests/excluded_nodes.tsv"
    if stage5_validate_bound_runtime >/dev/null 2>&1; then echo 'changed exclusion manifest accepted' >&2; exit 1; fi
    ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/excluded_nodes.tsv"
    if stage5_validate_bound_runtime >/dev/null 2>&1; then echo 'resealed different exclusion list accepted' >&2; exit 1; fi
  elif grep -q '^--exclude' "${CAPTURE}"; then echo 'default submissions unexpectedly excluded nodes' >&2; exit 1; fi
done
# Actual method/initial-watchdog/retry arguments, preserving resource policies.
TEST_EXCLUDED_NODES="${NODES}" /bin/bash "${ROOT}/tests/test_stage5_batch_hardware_policy.sh"
echo 'PASS: finite native exclusion list on preflight/Renv/metadata/method/watchdog/aggregate/retry; default and hardware policies unchanged'
