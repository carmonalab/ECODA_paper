#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
reference_stage_line="$(awk '/1\.0_stage_reference_maps\.sh/{print NR; exit}' "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh")"
runtime_validation_line="$(awk '/ecoda_runtime_validate_submission/{print NR; exit}' "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh")"
[[ -n "${reference_stage_line}" && -n "${runtime_validation_line}" &&
   ${reference_stage_line} -lt ${runtime_validation_line} ]] || {
  echo "Stage 4 runtime validation must follow reference-map staging." >&2
  exit 1
}
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
[[ "$(printf '%s\n' "${CALLS}" | wc -l | tr -d '[:space:]')" == "6" ]]
case "${CALLS}" in *"--array=1-9%1000"*) ;; *) echo "runnable arrays did not use nine rows" >&2; exit 1 ;; esac
case "${CALLS}" in *"1.2_prepare_chunks_worker.sh"*) ;; *) echo "preparation worker array missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--time=02:00:00"*"2.1_run_worker.sh"*) ;; *) echo "annotation worker time limit missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--time=02:00:00"*"3.2_merge_worker.sh"*) ;; *) echo "merge worker time limit missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_MODE=host"*"ECODA_RUNTIME_PROFILE=stage4"*) ;; *) echo "Stage 4 runtime export missing" >&2; exit 1 ;; esac
: > "${CAPTURE}"
DATASET_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  ANNOTATION_SUBMITTER_TEST=1 \
  bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
    --datasets _debug --views benchmark_analysis,batch_effect_uncorrected
)"
DATASET_MANIFEST="$(printf '%s\n' "${DATASET_OUTPUT}" | sed -n 's/^ANNOTATION_SELECTION_MANIFEST=//p')"
[[ "$(wc -l < "${DATASET_MANIFEST}" | tr -d '[:space:]')" == 2 ]]
[[ "$(sed -n '1p' "${DATASET_MANIFEST}")" == $'_debug\tbenchmark_analysis' ]]
[[ "$(sed -n '2p' "${DATASET_MANIFEST}")" == $'_debug\tbatch_effect_uncorrected' ]]
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  ECODA_RUNTIME_MODE=apptainer ANNOTATION_SUBMITTER_TEST=1 \
  bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  --selection-file "${SELECTION}" --exact-batch-selection >/dev/null 2>&1; then
  echo "Stage 4 accepted missing immutable runtime image" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
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
