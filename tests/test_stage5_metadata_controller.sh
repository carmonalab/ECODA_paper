#!/bin/bash
# Exercise the actual coordinator; stub only Slurm and the semantic reader.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-metadata-controller.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT

run_case() (
  export FAULT="$1" PASS_ARG="$2" HPC_SCRATCH_DIR="${TMP}/$1-$2"
  mkdir -p "${HPC_SCRATCH_DIR}"
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  eval "$(awk '/^stage5_export_final_sample_metadata\(\)/ {copy=1} copy {print} copy && /^}$/ {exit}' \
    "${ROOT}/src/5_run_benchmark_methods/stage5_dispatcher.sh")"
  eval "$(awk '/^stage5_record_scheduler\(\)/ {copy=1} copy {print} copy && /^}$/ {exit}' \
    "${ROOT}/src/5_run_benchmark_methods/stage5_dispatcher.sh")"
  export ECODA_RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/fixture" RUN_ID=fixture
  export ANALYSIS_PASS="${PASS_ARG}"
  if [[ "${PASS_ARG}" == corrected ]]; then export ANALYSIS_VARIANT=corrected_final;
  else export ANALYSIS_VARIANT=final; fi
  export ANALYSIS_ROOT="${ECODA_RUN_ROOT}/artifacts/batch_effect/${PASS_ARG}_final"
  export ANALYSIS_NAS_ROOT="${TMP}/nas/batch_effect/${PASS_ARG}_final"
  export DATASETS_JSON_FILE="${ROOT}/datasets.json"
  ANALYSIS_LOG_PREFIX=fixture; ECODA_RECOVERY_ACTIVE=0; SYNC_ONLY_RUN=""
  BENCHMARK_MATRIX_TEST=0; PARTITION_ARG=fixture; THROTTLE=8; MEMORY=1G
  WATCHDOG_TIME_LIMIT=00:01:00; USER_EMAIL=fixture; RUNTIME_EXPORT=fixture
  MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  SCHEDULER_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
  mkdir -p "${ECODA_RUN_ROOT}/manifests"
  : > "${SCHEDULER_FILE}"
  for DS in Joanito Stephenson Alzheimer Breast_cancer Covid19_PBMC Kidney_KPMP_full Lupus_PBMC Lung; do
    printf '%s\tbatch_effect_%s\tall\n' "${DS}" "${PASS_ARG}" >> "${MANIFEST}"
  done
  stage5_source_script() { printf '%s/%s\n' "${ROOT}" "$1"; }
  stage5_require_source_script() { [[ -f "$1" ]]; }
  ecoda_stage5_validate_identity() { [[ "$1" == "${PASS_ARG}" ]]; }
  stage5_validate_bound_runtime() { return 0; }
  stage5_input_path() { printf '%s/never-opened-%s.h5ad\n' "${HPC_SCRATCH_DIR}" "$1"; }
  PYTHON_BIN="${HPC_SCRATCH_DIR}/semantic-check"
  cat > "${PYTHON_BIN}" <<'SH'
#!/bin/bash
set -euo pipefail
shift
OUTPUT=""; DATASET=""
while [[ $# -gt 0 ]]; do
  case "$1" in --output) OUTPUT="$2"; shift 2 ;; --dataset) DATASET="$2"; shift 2 ;; --check) shift ;; *) shift 2 ;; esac
done
[[ "${ANALYSIS_VARIANT}" == final || "${ANALYSIS_VARIANT}" == corrected_final ]]
[[ -s "${OUTPUT}" ]] || exit 1
if [[ "${FAULT}" == semantic && "${DATASET}" == Joanito ]]; then
  echo 'semantic reader detail preserved' >&2; exit 9
fi
SH
  chmod +x "${PYTHON_BIN}"
  sbatch() {
    local DS VIEW INPUT OUTPUT EXTRA INDEX=0 STATUS
    while IFS=$'\t' read -r DS VIEW INPUT OUTPUT EXTRA; do
      INDEX=$((INDEX + 1)); mkdir -p "$(dirname "${OUTPUT}")"
      printf 'preserved tiny metadata %s\n' "${DS}" > "${OUTPUT}"
      ecoda_write_checksum "${OUTPUT}"
      STATUS="${ECODA_RUN_ROOT}/status/metadata_export/${DS}__${VIEW}.status"
      printf 'STATE=OK\nTASK_ID=%s\n' "${INDEX}" > "${STATUS}"
      if [[ "${DS}" == Joanito ]]; then
        case "${FAULT}" in
          missing-status) rm "${STATUS}" ;;
          state) printf 'STATE=FAIL\nTASK_ID=1\n' > "${STATUS}" ;;
          task) printf 'STATE=OK\nTASK_ID=2\n' > "${STATUS}" ;;
        esac
      fi
    done < "${ECODA_RUN_ROOT}/manifests/metadata_export_pending.tsv"
    if [[ "${FAULT}" == parse ]]; then printf 'invalid-id\n'; else printf '123456;fixture\n'; fi
    [[ "${FAULT}" != wait ]] || return 7
  }
  if [[ "${FAULT}" == record ]]; then stage5_record_scheduler() { return 1; }; fi
  if [[ "${FAULT}" == checksum ]]; then
    eval "$(declare -f ecoda_validate_checksum | sed '1s/ecoda_validate_checksum/real_checksum/')"
    ecoda_validate_checksum() { [[ "$1" != *Joanito_sample_metadata.feather ]] && real_checksum "$@"; }
  fi
  if [[ "${FAULT}" == report || "${FAULT}" == report-checksum ]]; then
    eval "$(declare -f ecoda_atomic_write | sed '1s/ecoda_atomic_write/real_write/')"
    ecoda_atomic_write() { [[ "${FAULT}" != report || "$1" != */metadata_export.report ]] && real_write "$@"; }
    eval "$(declare -f ecoda_write_checksum | sed '1s/ecoda_write_checksum/real_write_checksum/')"
    ecoda_write_checksum() { [[ "${FAULT}" != report-checksum || "$1" != */metadata_export.report ]] && real_write_checksum "$@"; }
  fi
  if [[ "${FAULT}" == count ]]; then
    grep() { [[ "$*" != *Lung* ]] && command grep "$@"; }
  fi
  RC=0; stage5_export_final_sample_metadata 0 || RC=$?
  for DS in Joanito Stephenson Alzheimer Breast_cancer Covid19_PBMC Kidney_KPMP_full Lupus_PBMC Lung; do
    OUTPUT="${ANALYSIS_ROOT}/metadata/${DS}_sample_metadata.feather"
    [[ "$(cat "${OUTPUT}")" == "preserved tiny metadata ${DS}" ]]
  done
  if [[ "${FAULT}" == valid ]]; then
    [[ ${RC} -eq 0 ]]
    command grep -q '^COUNT=8$' "${ECODA_RUN_ROOT}/status/metadata_export.report"
    command grep -q '^PENDING=8$' "${ECODA_RUN_ROOT}/status/metadata_export.report"
    ecoda_validate_checksum "${ECODA_RUN_ROOT}/status/metadata_export.report"
  else [[ ${RC} -ne 0 ]]; fi
)

for LANE in uncorrected corrected; do run_case valid "${LANE}"; done
expect_failure() {
  run_case "$1" corrected > "${TMP}/$1.log" 2>&1
  grep -F -q "$2" "${TMP}/$1.log" || { cat "${TMP}/$1.log" >&2; return 1; }
}
expect_failure parse 'scheduler ID is invalid'
expect_failure record 'could not record metadata export scheduler ID'
expect_failure wait 'scheduler client failed for 123456 (exit 7)'
expect_failure semantic 'semantic validation failed for Joanito/batch_effect_corrected'
grep -q 'semantic reader detail preserved' "${TMP}/semantic.log"
expect_failure checksum 'checksum validation failed for Joanito/batch_effect_corrected'
expect_failure missing-status 'status is missing for Joanito/batch_effect_corrected'
expect_failure state 'STATE=FAIL TASK_ID=1, expected OK/1'
expect_failure task 'STATE=OK TASK_ID=2, expected OK/1'
expect_failure count 'pending count mismatch: validated 7, expected 8'
expect_failure report 'could not write metadata export report'
expect_failure report-checksum 'could not checksum metadata export report'
echo 'Stage 5 metadata coordinator diagnostics: PASS'
