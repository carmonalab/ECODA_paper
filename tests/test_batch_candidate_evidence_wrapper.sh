#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-evidence-wrapper.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
emit_csv() {
  local path="$1" dataset="$2"
  local header='dataset,method,method_available,method_applicable,method_reason,artifact,candidate,present,completeness,levels,samples_per_level,nmi_biology,constant_candidate,sample_unique_candidate,perfect_confounded,marginal_r2,marginal_p,joint_r2,joint_p,warnings,marginal_p_holm,joint_p_holm'
  printf '%s\n' "${header}" > "${path}"
  printf '%s\n' "${dataset},test,TRUE,TRUE,test,artifact,candidate,TRUE,1,2,A:1;B:1,0,FALSE,FALSE,FALSE,0,0,0,0,,0,0" >> "${path}"
}
emit_outputs() {
  local output_dir="${EVIDENCE_OUTPUT_DIR:?}"
  local ds path digest size
  mkdir -p "${output_dir}"
  for ds in Alzheimer Breast_cancer Covid19_PBMC Kidney_KPMP Myocardial_infarction Diabetes Lupus_PBMC Lung Parkinson Joanito Stephenson CombinedPBMC; do
    path="${output_dir}/${ds}_batch_candidate_evidence.csv"
    emit_csv "${path}" "${ds}"
    digest="$(md5sum "${path}" | cut -d' ' -f1)"
    size="$(wc -c < "${path}" | tr -d '[:space:]')"
    printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "${size}" "${path}" > "${path}.md5"
  done
  path="${output_dir}/batch_candidate_review.csv"
  emit_csv "${path}" Alzheimer
  digest="$(md5sum "${path}" | cut -d' ' -f1)"
  size="$(wc -c < "${path}" | tr -d '[:space:]')"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "${size}" "${path}" > "${path}.md5"
}
if [[ "${SBATCH_FAIL:-0}" == "1" ]]; then
  printf '99002\n'
  exit 7
fi
if [[ "${SBATCH_EMIT_EVIDENCE:-0}" == "1" ]]; then
  emit_outputs
fi
printf '%s\n' "${SBATCH_ID:-99001}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
cat > "${TMP_DIR}/bin/rsync" <<'STUB'
#!/bin/bash
set -euo pipefail
if [[ "${RSYNC_FAIL:-0}" == "1" ]]; then
  exit 9
fi
n="$#"
eval "source_file=\${$((n - 2))}"
eval "source_sidecar=\${$((n - 1))}"
eval "destination=\${$n}"
mkdir -p "${destination}"
cp "${source_file}" "${source_sidecar}" "${destination}"
STUB
chmod +x "${TMP_DIR}/bin/rsync"
SELECTION="${TMP_DIR}/selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\nKidney_KPMP\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\n' > "${SELECTION}"
digest="$(md5sum "${SELECTION}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${SELECTION}" | tr -d '[:space:]')" "${SELECTION}" > "${SELECTION}.md5"
OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  EVIDENCE_SUBMITTER_TEST=1 bash "${ROOT}/notebooks/dataset_onboarding/submit_batch_candidate_evidence.sh" \
    --selection-file "${SELECTION}" \
    --analysis-root "${TMP_DIR}/analysis" --input-root "${TMP_DIR}/input" \
    --output-dir "${TMP_DIR}/evidence" --run-id evidence-test
)"
case "${OUTPUT}" in *"BATCH_EFFECT_EVIDENCE_SCHEDULER_ID=99001"*) ;; *) echo "evidence scheduler marker missing" >&2; exit 1 ;; esac
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == 1 ]]
case "$(cat "${CAPTURE}")" in *"--wait"*) ;; *) echo "evidence wrapper did not wait for job" >&2; exit 1 ;; esac
case "$(cat "${CAPTURE}")" in *"build_batch_candidate_evidence.R"*) ;; *) echo "evidence worker command missing" >&2; exit 1 ;; esac
: > "${CAPTURE}"
BAD="${TMP_DIR}/bad.tsv"
sed '1s/batch_effect_uncorrected/batch_effect_corrected/' "${SELECTION}" > "${BAD}"
set +e
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  EVIDENCE_SUBMITTER_TEST=1 bash "${ROOT}/notebooks/dataset_onboarding/submit_batch_candidate_evidence.sh" \
    --selection-file "${BAD}" --analysis-root "${TMP_DIR}/analysis" --input-root "${TMP_DIR}/input" \
    --output-dir "${TMP_DIR}/evidence-bad" --run-id evidence-bad >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
[[ ! -s "${CAPTURE}" ]]
set +e
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  SBATCH_FAIL=1 bash "${ROOT}/notebooks/dataset_onboarding/submit_batch_candidate_evidence.sh" \
    --selection-file "${SELECTION}" --analysis-root "${TMP_DIR}/analysis" --input-root "${TMP_DIR}/input" \
    --output-dir "${TMP_DIR}/evidence-failed" --run-id evidence-failed >/dev/null 2>&1
FAIL_RC=$?
set -e
[[ ${FAIL_RC} -ne 0 ]]
FAIL_RUN="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/evidence-failed"
[[ "$(grep '^STATE=' "${FAIL_RUN}/status/terminal")" == "STATE=FAIL" ]]
[[ "$(cat "${FAIL_RUN}/manifests/scheduler_ids.tsv")" == $'EVIDENCE\t99002' ]]

NAS_ROOT="${TMP_DIR}/nas"
mkdir -p "${NAS_ROOT}"
SUCCESS_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  EVIDENCE_SUBMITTER_TEST_SYNC=1 EVIDENCE_TEST_NAS_TARGET_DIR="${NAS_ROOT}" \
  EVIDENCE_OUTPUT_DIR="${TMP_DIR}/evidence-success" SBATCH_EMIT_EVIDENCE=1 SBATCH_ID=99003 \
  bash "${ROOT}/notebooks/dataset_onboarding/submit_batch_candidate_evidence.sh" \
    --selection-file "${SELECTION}" --analysis-root "${TMP_DIR}/analysis" --input-root "${TMP_DIR}/input" \
    --output-dir "${TMP_DIR}/evidence-success" --run-id evidence-success
)"
case "${SUCCESS_OUTPUT}" in *"EVIDENCE_RUN_ID=evidence-success"*) ;; *) echo "evidence sync success marker missing" >&2; exit 1 ;; esac
SUCCESS_RUN="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/evidence-success"
[[ "$(grep '^STATE=' "${SUCCESS_RUN}/status/terminal")" == "STATE=OK" ]]
SUCCESS_REMOTE="${NAS_ROOT}/batch_effect/uncorrected/evidence/evidence-success"
for ds in Alzheimer Breast_cancer Covid19_PBMC Kidney_KPMP Myocardial_infarction Diabetes Lupus_PBMC Lung Parkinson Joanito Stephenson CombinedPBMC; do
  [[ -s "${SUCCESS_REMOTE}/${ds}_batch_candidate_evidence.csv" ]]
  [[ -s "${SUCCESS_REMOTE}/${ds}_batch_candidate_evidence.csv.md5" ]]
done
[[ -s "${SUCCESS_REMOTE}/batch_candidate_review.csv" ]]
[[ -s "${SUCCESS_REMOTE}/batch_candidate_review.csv.md5" ]]
[[ -s "${SUCCESS_REMOTE}/checksums.md5" ]]
[[ -s "${SUCCESS_REMOTE}/checksums.md5.md5" ]]
[[ "$(printf '%s\n' "${SUCCESS_REMOTE}"/* | wc -l | tr -d '[:space:]')" == 28 ]]

set +e
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  EVIDENCE_SUBMITTER_TEST_SYNC=1 EVIDENCE_TEST_NAS_TARGET_DIR="${NAS_ROOT}" \
  EVIDENCE_OUTPUT_DIR="${TMP_DIR}/evidence-transfer-failed" \
  SBATCH_EMIT_EVIDENCE=1 SBATCH_ID=99004 RSYNC_FAIL=1 \
  bash "${ROOT}/notebooks/dataset_onboarding/submit_batch_candidate_evidence.sh" \
    --selection-file "${SELECTION}" --analysis-root "${TMP_DIR}/analysis" --input-root "${TMP_DIR}/input" \
    --output-dir "${TMP_DIR}/evidence-transfer-failed" --run-id evidence-transfer-failed >/dev/null 2>&1
TRANSFER_FAIL_RC=$?
set -e
[[ ${TRANSFER_FAIL_RC} -ne 0 ]]
TRANSFER_FAIL_RUN="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/evidence-transfer-failed"
[[ "$(grep '^STATE=' "${TRANSFER_FAIL_RUN}/status/terminal")" == "STATE=FAIL" ]]
[[ "$(cat "${TRANSFER_FAIL_RUN}/manifests/scheduler_ids.tsv")" == $'EVIDENCE\t99004' ]]
[[ -s "${TMP_DIR}/evidence-transfer-failed/Alzheimer_batch_candidate_evidence.csv" ]]
[[ ! -e "${NAS_ROOT}/batch_effect/uncorrected/evidence/evidence-transfer-failed/Alzheimer_batch_candidate_evidence.csv" ]]

echo "evidence wrapper: OK"
