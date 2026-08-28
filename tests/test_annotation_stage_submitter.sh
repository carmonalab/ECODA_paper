#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-annotation-stage.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '72000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
SELECTION="${TMP_DIR}/selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\nKidney_KPMP\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\n' > "${SELECTION}"
OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  ANNOTATION_SUBMITTER_TEST=1 \
  bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
    --selection-file "${SELECTION}" --exact-batch-selection
)"
RUN_ID="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^ANNOTATION_RUN_ID=//p')"
[[ -n "${RUN_ID}" ]]
RUN_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}"
SKIPS="${RUN_ROOT}/status/skipped"
[[ "$(cat "${SKIPS}")" == $'Alzheimer\tSKIP_NOT_SUITABLE\nDiabetes\tSKIP_NOT_SUITABLE\nParkinson\tSKIP_NOT_SUITABLE' ]]
RUNNABLE="${RUN_ROOT}/manifests/runnable_selection.tsv"
[[ "$(wc -l < "${RUNNABLE}" | tr -d '[:space:]')" == 9 ]]
! grep -qE '^(Alzheimer|Diabetes|Parkinson)\t' "${RUNNABLE}"
[[ "$(wc -l < "${RUN_ROOT}/manifests/scheduler_ids.tsv" | tr -d '[:space:]')" == 6 ]]
CALLS="$(cat "${CAPTURE}")"
[[ "$(printf '%s\n' "${CALLS}" | wc -l | tr -d '[:space:]')" == 6 ]]
case "${CALLS}" in *"--array=1-9%1000"*) ;; *) echo "runnable arrays did not use nine rows" >&2; exit 1 ;; esac
case "${CALLS}" in *"1.2_prepare_chunks_worker.sh"*) ;; *) echo "preparation worker array missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"2.1_run_worker.sh"*) ;; *) echo "annotation worker array missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"3.2_merge_worker.sh"*) ;; *) echo "merge worker array missing" >&2; exit 1 ;; esac
: > "${CAPTURE}"
BAD="${TMP_DIR}/bad.tsv"
sed '1s/batch_effect_uncorrected/batch_effect_corrected/' "${SELECTION}" > "${BAD}"
[[ ! -e "${RUN_ROOT}/datasets/Alzheimer" ]]
[[ ! -e "${RUN_ROOT}/datasets/Diabetes" ]]
[[ ! -e "${RUN_ROOT}/datasets/Parkinson" ]]
set +e
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  ANNOTATION_SUBMITTER_TEST=1 bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  --selection-file "${BAD}" --exact-batch-selection >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
[[ ! -s "${CAPTURE}" ]]
for bad_kind in legacy missing; do
  : > "${CAPTURE}"
  if [[ "${bad_kind}" == legacy ]]; then
    sed '1s/batch_effect_uncorrected/batch_effect_analysis/' "${SELECTION}" > "${BAD}"
  else
    sed '12d' "${SELECTION}" > "${BAD}"
  fi
  set +e
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
    ANNOTATION_SUBMITTER_TEST=1 bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
    --selection-file "${BAD}" --exact-batch-selection >/dev/null 2>&1
  RC=$?
  set -e
  [[ ${RC} -ne 0 ]]
  [[ ! -s "${CAPTURE}" ]]
done
echo "annotation stage submitter: OK"
