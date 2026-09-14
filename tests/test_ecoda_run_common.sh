#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-common.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
NAS_TEST_ROOT="$(cd "${TMP_DIR}" && pwd -P)"
NAS_TEST_HOST_ENV="${NAS_TEST_ROOT}/fake-host-env"
mkdir -p "${NAS_TEST_HOST_ENV}/bin"
printf '#!/bin/sh\n' > "${NAS_TEST_HOST_ENV}/bin/python"
printf '#!/bin/sh\n' > "${NAS_TEST_HOST_ENV}/bin/Rscript"
chmod +x "${NAS_TEST_HOST_ENV}/bin/python" "${NAS_TEST_HOST_ENV}/bin/Rscript"

run_nas_config() {
  local cluster="$1"
  local nas_prefix="${2:-}"
  local nas_sc_dir="${3:-}"
  local nas_target_dir="${4:-}"
  env -i \
    PATH="/usr/bin:/bin" \
    USER="ecoda-test" \
    ECODA_HPC_CLUSTER="${cluster}" \
    ECODA_HOST_ENV_PREFIX="${NAS_TEST_HOST_ENV}" \
    HPC_SCRATCH_DIR="${NAS_TEST_ROOT}/nas-scratch" \
    HOME="${NAS_TEST_ROOT}/nas-home" \
    ECODA_NAS_PREFIX="${nas_prefix}" \
    ECODA_NAS_SC_DIR="${nas_sc_dir}" \
    ECODA_NAS_TARGET_DIR="${nas_target_dir}" \
    bash -c '
      set -euo pipefail
      unset NAS_PREFIX NAS_BASE_DIR NAS_SC_DIR NAS_TARGET_DIR SLURM_CLUSTER_NAME
      [[ -n "${ECODA_NAS_PREFIX:-}" ]] || unset ECODA_NAS_PREFIX
      [[ -n "${ECODA_NAS_SC_DIR:-}" ]] || unset ECODA_NAS_SC_DIR
      [[ -n "${ECODA_NAS_TARGET_DIR:-}" ]] || unset ECODA_NAS_TARGET_DIR
      source "$1"
      printf "%s\t%s\t%s\n" \
        "${NAS_PREFIX}" "${NAS_SC_DIR}" "${NAS_TARGET_DIR}"
    ' _ "${ROOT}/src/slurm_config.sh"
}

BAMBOO_NAS="$(run_nas_config bamboo)"
[[ "${BAMBOO_NAS}" == \
  $'/srv/smednas515.unige.ch/carmona_smb\t/srv/smednas515.unige.ch/carmona_smb/DataCollections/Standardized_SingleCell_Datasets\t/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper' ]]

if YGGDRASIL_NAS="$(run_nas_config yggdrasil 2>&1)"; then
  echo "Yggdrasil configuration unexpectedly inherited a NAS default." >&2
  exit 1
else
  YGGDRASIL_RC=$?
fi
[[ ${YGGDRASIL_RC} -ne 0 ]]
[[ "${YGGDRASIL_NAS}" == *"explicit"* ]]
[[ "${YGGDRASIL_NAS}" == *"configuration"* ]]

YGG_PREFIX="${NAS_TEST_ROOT}/nas/ygg-prefix"
YGG_SC_DIR="${NAS_TEST_ROOT}/nas/ygg-sc"
YGG_TARGET_DIR="${NAS_TEST_ROOT}/nas/ygg-target"
YGGDRASIL_OVERRIDE_NAS="$(
  run_nas_config yggdrasil "${YGG_PREFIX}" "${YGG_SC_DIR}" "${YGG_TARGET_DIR}"
)"
[[ "${YGGDRASIL_OVERRIDE_NAS}" == \
  "${YGG_PREFIX}"$'\t'"${YGG_SC_DIR}"$'\t'"${YGG_TARGET_DIR}" ]]

for relative_nas_var in ECODA_NAS_PREFIX ECODA_NAS_SC_DIR ECODA_NAS_TARGET_DIR; do
  relative_prefix="${YGG_PREFIX}"
  relative_sc_dir="${YGG_SC_DIR}"
  relative_target_dir="${YGG_TARGET_DIR}"
  case "${relative_nas_var}" in
    ECODA_NAS_PREFIX) relative_prefix="relative/nas-prefix" ;;
    ECODA_NAS_SC_DIR) relative_sc_dir="relative/nas-sc" ;;
    ECODA_NAS_TARGET_DIR) relative_target_dir="relative/nas-target" ;;
  esac
  if relative_nas_output="$(
    run_nas_config yggdrasil \
      "${relative_prefix}" "${relative_sc_dir}" "${relative_target_dir}" \
      2>&1
  )"; then
    echo "${relative_nas_var} unexpectedly accepted a relative path." >&2
    exit 1
  else
    relative_nas_rc=$?
  fi
  [[ ${relative_nas_rc} -ne 0 ]]
  [[ "${relative_nas_output}" == *"absolute"* ]]
done

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
printf 'Alzheimer\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\nKidney_KPMP_full\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\n' > "${BATCH_MANIFEST}"
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
CORRECTED_BATCH_BASE="${TMP_DIR}/corrected-batch-base.json"
cat > "${CORRECTED_BATCH_BASE}" <<'JSON'
{"Fixture":{"columns":{"sample":"sample_id","label":"label","batch":"batch_a"},"views":{"batch_effect_corrected":{"columns":{}}}}}
JSON
CORRECTED_BATCH_STATE="${TMP_DIR}/corrected-batch-state"
CORRECTED_BATCH_SBATCH="${TMP_DIR}/corrected-batch.sbatch.calls"
: > "${CORRECTED_BATCH_SBATCH}"
for corrected_case in null object empty duplicate blank label_overlap sample_overlap reserved_sample; do
  case "${corrected_case}" in
    null) raw_batch='null' ;;
    object) raw_batch='{"name":"batch_a"}' ;;
    empty) raw_batch='[]' ;;
    duplicate) raw_batch='["batch_a","batch_a"]' ;;
    blank) raw_batch='["   "]' ;;
    label_overlap) raw_batch='["label"]' ;;
    sample_overlap) raw_batch='["sample_id"]' ;;
    reserved_sample) raw_batch='["Sample"]' ;;
  esac
  corrected_fixture="${TMP_DIR}/corrected-${corrected_case}.json"
  jq --argjson batch "${raw_batch}" \
    '.Fixture.columns.batch = $batch' "${CORRECTED_BATCH_BASE}" \
    > "${corrected_fixture}"
  corrected_run_root="${CORRECTED_BATCH_STATE}/_ecoda_runs/${corrected_case}"
  corrected_runs_before="$(printf '%s\n' "${ECODA_RUNS_ROOT}"/*)"
  set +e
  ecoda_validate_corrected_batch_columns \
    "${corrected_fixture}" Fixture batch_effect_corrected >/dev/null 2>&1
  RC=$?
  set -e
  [[ ${RC} -ne 0 ]]
  [[ ! -s "${CORRECTED_BATCH_SBATCH}" ]]
  [[ "${corrected_runs_before}" == "$(printf '%s\n' "${ECODA_RUNS_ROOT}"/*)" ]]
  [[ ! -e "${CORRECTED_BATCH_STATE}" ]]
  [[ ! -e "${corrected_run_root}" ]]
  [[ ! -e "${corrected_run_root}/manifests/scheduler_ids.tsv" ]]
  [[ ! -e "${corrected_run_root}/manifests/selection.tsv" ]]
  [[ ! -e "${corrected_run_root}/manifests/pending.tsv" ]]
  [[ ! -e "${corrected_run_root}/status/compute" ]]
done
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
GLOBAL_ARTIFACT="${TMP_DIR}/global-artifact"
printf 'global artifact\n' > "${GLOBAL_ARTIFACT}"
GLOBAL_OWNER="$(ecoda_artifact_owner_dir "${GLOBAL_ARTIFACT}")"
mkdir -p "$(dirname "${GLOBAL_OWNER}")"
SYMLINK_OWNER_TARGET="${TMP_DIR}/global-owner-target"
mkdir -p "${SYMLINK_OWNER_TARGET}"
GLOBAL_CANONICAL="$(ecoda_canonical_path "${GLOBAL_ARTIFACT}")"
printf 'RUN_ID=symlink-owner\nSTAGE=stage5\nPATH=%s\nSTATE=ACTIVE\nPID=%s\n' \
  "${GLOBAL_CANONICAL}" "$$" > "${SYMLINK_OWNER_TARGET}/owner"
ln -s "${SYMLINK_OWNER_TARGET}" "${GLOBAL_OWNER}"
SYMLINK_OWNER_BEFORE="$(cat "${SYMLINK_OWNER_TARGET}/owner")"
set +e
ecoda_artifact_owner_acquire "${GLOBAL_ARTIFACT}" stage5 symlink-owner 1 1 0 \
  >/dev/null 2>&1
SYMLINK_ACQUIRE_RC=$?
ecoda_artifact_owner_validate "${GLOBAL_ARTIFACT}" symlink-owner \
  >/dev/null 2>&1
SYMLINK_VALIDATE_RC=$?
ecoda_artifact_owner_set_state "${GLOBAL_ARTIFACT}" OK should-reject \
  >/dev/null 2>&1
SYMLINK_SET_STATE_RC=$?
set -e
[[ ${SYMLINK_ACQUIRE_RC} -ne 0 ]]
[[ ${SYMLINK_VALIDATE_RC} -ne 0 ]]
[[ ${SYMLINK_SET_STATE_RC} -ne 0 ]]
[[ "$(cat "${SYMLINK_OWNER_TARGET}/owner")" == "${SYMLINK_OWNER_BEFORE}" ]]
printf 'complete artifact\n' > "${TMP_DIR}/artifact"
ecoda_write_checksum "${TMP_DIR}/artifact"
ecoda_validate_checksum "${TMP_DIR}/artifact"
printf 'corrupt\n' >> "${TMP_DIR}/artifact"
set +e
ecoda_validate_checksum "${TMP_DIR}/artifact"
RC=$?
set -e
[[ ${RC} -ne 0 ]]
printf 'complete artifact\n' > "${TMP_DIR}/artifact"
ecoda_write_checksum "${TMP_DIR}/artifact"
remote_artifact="${TMP_DIR}/remote-artifact"
cp "${TMP_DIR}/artifact" "${remote_artifact}"
cp "${TMP_DIR}/artifact.md5" "${remote_artifact}.md5"
ecoda_compare_checksum_remote "${TMP_DIR}/artifact" "${remote_artifact}" "${remote_artifact}.md5"
[[ "$(sed -n 's/^PATH=//p' "${remote_artifact}.md5")" == "${remote_artifact}" ]]
ecoda_validate_checksum "${remote_artifact}" "${remote_artifact}.md5"
SYMLINK_INPUT_REAL="${TMP_DIR}/input-real"
SYMLINK_INPUT_LINK="${TMP_DIR}/input-link"
mkdir -p "${SYMLINK_INPUT_REAL}"
ln -s "${SYMLINK_INPUT_REAL}" "${SYMLINK_INPUT_LINK}"
SYMLINK_INPUT="${SYMLINK_INPUT_LINK}/artifact.bin"
printf 'symlinked input artifact\n' > "${SYMLINK_INPUT}"
ecoda_write_checksum "${SYMLINK_INPUT}"
SYMLINK_PRODUCER_RUN="symlink_input_$$"
ecoda_init_run stage3 "${SYMLINK_PRODUCER_RUN}" >/dev/null
ecoda_write_artifact_record "${SYMLINK_INPUT}" stage3 "${SYMLINK_PRODUCER_RUN}" >/dev/null
SYMLINK_INPUT_OWNER="$(ecoda_artifact_owner_acquire \
  "${SYMLINK_INPUT}" stage3 "${SYMLINK_PRODUCER_RUN}" 0)"
ecoda_artifact_owner_set_state "${SYMLINK_INPUT}" OK \
  "symlinked input regression" >/dev/null
ecoda_validate_input_artifact \
  "${SYMLINK_INPUT}" stage3 "${SYMLINK_PRODUCER_RUN}" >/dev/null
SYMLINK_INPUT_CANONICAL="$(ecoda_canonical_path "${SYMLINK_INPUT}")"
ecoda_validate_checksum "${SYMLINK_INPUT_CANONICAL}" \
  "${SYMLINK_INPUT}.md5"
ecoda_validate_checksum_record "${SYMLINK_INPUT_CANONICAL}" \
  "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}" \
  "${SYMLINK_INPUT}.md5"
ecoda_validate_input_artifact \
  "${SYMLINK_INPUT_CANONICAL}" stage3 "${SYMLINK_PRODUCER_RUN}" >/dev/null
sacct() { return 0; }
ECODA_ACCOUNTING_EMPTY_GRACE=2
set +e
ecoda_wait_array_accounting missing_job 1 0 >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
ACCOUNTING_STATE="${TMP_DIR}/accounting.calls"
printf '0\n' > "${ACCOUNTING_STATE}"
sacct() {
  local calls
  calls="$(cat "${ACCOUNTING_STATE}")"
  calls=$((calls + 1))
  printf '%s\n' "${calls}" > "${ACCOUNTING_STATE}"
  if [[ ${calls} -eq 1 ]]; then
    printf 'active_job_1|PENDING|0:0\n'
  else
    printf 'active_job_1|COMPLETED|0:0\nactive_job_2|COMPLETED|0:0\n'
  fi
}
squeue() { printf 'active_job_[1-2%%2]\n'; }
ecoda_wait_array_accounting active_job 2 0 >/dev/null
[[ "${ECODA_ACCOUNTING_ROWS}" == *"active_job_2|COMPLETED|0:0"* ]]
