#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-common.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
export HPC_SCRATCH_DIR="${TMP_DIR}/scratch"
export LOGS_DIR="${TMP_DIR}/logs"
export DATASETS_JSON_FILE="${ROOT}/datasets.json"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
RUN_ID="test_run_$$"
ecoda_init_run test "${RUN_ID}" >/dev/null
MANIFEST_SOURCE="${TMP_DIR}/manifest.source"
MANIFEST_DEST="${TMP_DIR}/manifest.dest"
printf 'Adams\tbenchmark_analysis\n' > "${MANIFEST_SOURCE}"
ecoda_atomic_install_manifest "${MANIFEST_SOURCE}" "${MANIFEST_DEST}" 2
[[ "$(cat "${MANIFEST_DEST}")" == $'Adams\tbenchmark_analysis' ]]
printf 'Bassez\tbenchmark_analysis\npartial\n' > "${MANIFEST_SOURCE}"
set +e
ecoda_atomic_install_manifest "${MANIFEST_SOURCE}" "${MANIFEST_DEST}" 2 >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
[[ "$(cat "${MANIFEST_DEST}")" == $'Adams\tbenchmark_analysis' ]]
BATCH_MANIFEST="${TMP_DIR}/batch-selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\nKidney_KPMP\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\n' > "${BATCH_MANIFEST}"
ecoda_validate_exact_batch_selection "${BATCH_MANIFEST}" 2
BATCH_MATRIX="${TMP_DIR}/batch-matrix.tsv"
sed 's/$/\tbatch_effect_uncorrected/' "${BATCH_MANIFEST}" > "${BATCH_MATRIX}"
ecoda_validate_exact_batch_selection "${BATCH_MATRIX}" 3
for replacement in batch_effect_analysis batch_effect_corrected; do
  invalid="${TMP_DIR}/${replacement}.tsv"
  sed "1s/batch_effect_uncorrected/${replacement}/" "${BATCH_MANIFEST}" > "${invalid}"
  set +e
  ecoda_validate_exact_batch_selection "${invalid}" 2 >/dev/null 2>&1
  RC=$?
  set -e
  [[ ${RC} -ne 0 ]]
done
invalid="${TMP_DIR}/duplicate.tsv"
sed '12d' "${BATCH_MANIFEST}" | sed '2s/.*/Alzheimer\tbatch_effect_uncorrected/' > "${invalid}"
set +e
ecoda_validate_exact_batch_selection "${invalid}" 2 >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
OWNER="$(ecoda_owner_acquire test dataset/view "${RUN_ID}" 0)"
[[ -d "${OWNER}" ]]
set +e
ecoda_owner_acquire test dataset/view other_run 0 >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -eq 1 ]]
ecoda_owner_set_state "${OWNER}" OK done
set +e
ecoda_owner_acquire test dataset/view other_run 0 >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -eq 2 ]]
OWNER2="$(ecoda_owner_acquire test dataset/view other_run 0 0)"
[[ "${OWNER2}" == "${OWNER}" ]]
ecoda_owner_set_state "${OWNER2}" FAIL retry
OWNER3="$(ecoda_owner_acquire test dataset/view final_run 1)"
[[ "${OWNER3}" == "${OWNER}" ]]
[[ "$(ecoda_owner_state "${OWNER3}")" == "ACTIVE" ]]
printf 'complete artifact\n' > "${TMP_DIR}/artifact"
ecoda_write_checksum "${TMP_DIR}/artifact"
ecoda_validate_checksum "${TMP_DIR}/artifact"
printf 'corrupt\n' >> "${TMP_DIR}/artifact"
set +e
ecoda_validate_checksum "${TMP_DIR}/artifact"
RC=$?
set -e
[[ ${RC} -ne 0 ]]
sacct() { return 0; }
ECODA_ACCOUNTING_EMPTY_GRACE=2
set +e
ecoda_wait_array_accounting missing_job 1 0 >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
sacct() {
  ACCOUNTING_CALLS=$((ACCOUNTING_CALLS + 1))
  if [[ ${ACCOUNTING_CALLS} -eq 1 ]]; then
    printf 'active_job_1|PENDING|0:0\n'
  else
    printf 'active_job_1|COMPLETED|0:0\nactive_job_2|COMPLETED|0:0\n'
  fi
}
squeue() { printf 'active_job_[1-2%%2]\n'; }
ACCOUNTING_CALLS=0
ecoda_wait_array_accounting active_job 2 0 >/dev/null
[[ "${ECODA_ACCOUNTING_ROWS}" == *"active_job_2|COMPLETED|0:0"* ]]
echo "common ownership/checksum: OK"
