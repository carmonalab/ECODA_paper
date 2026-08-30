#!/bin/bash
# Focused contract test for the canonical manifest-driven preprocessing gate.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-preprocess-stage.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
CALL_COUNT="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
case "${CALL_COUNT}" in
  1) printf '600001\n' ;;
  2) printf '600002\n' ;;
  *) printf 'unexpected sbatch call count\n' >&2; exit 1 ;;
esac
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
SELECTION="${TMP_DIR}/selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\nKidney_KPMP\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\n' > "${SELECTION}"
OUTPUT="$(
  HOME="${TMP_DIR}/home" \
  PATH="${TMP_DIR}/bin:${PATH}" \
  USER_EMAIL="test@example.invalid" \
  PREPROCESS_SUBMITTER_TEST=1 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
    --selection-file "${SELECTION}" --exact-batch-selection
)"
case "${OUTPUT}" in *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;; *) echo "missing array marker" >&2; exit 1 ;; esac
case "${OUTPUT}" in *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) ;; *) echo "missing watchdog marker" >&2; exit 1 ;; esac
MANIFEST="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
[[ -s "${MANIFEST}" ]]
[[ "$(wc -l < "${MANIFEST}" | tr -d '[:space:]')" == 12 ]]
[[ "$(sed -n '1p' "${MANIFEST}")" == $'Alzheimer\tbatch_effect_uncorrected' ]]
[[ "$(sed -n '12p' "${MANIFEST}")" == $'CombinedPBMC\tbatch_effect_uncorrected' ]]
SCHEDULER_MANIFEST="$(dirname "${MANIFEST}")/scheduler_ids.tsv"
[[ "$(wc -l < "${SCHEDULER_MANIFEST}" | tr -d '[:space:]')" == 2 ]]
CALLS="$(cat "${CAPTURE}")"
case "${CALLS}" in *"--array=1-12%1000"*) ;; *) echo "array was not submitted with all selected rows" >&2; exit 1 ;; esac
case "${CALLS}" in *"--dependency=afterany:600001"*) ;; *) echo "watchdog dependency missing" >&2; exit 1 ;; esac
case "${CALLS}" in *PREPROCESS_SELECTION_FILE=*pending.tsv*) ;; *) echo "pending manifest was not exported" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_MODE=host"*"ECODA_RUNTIME_PROFILE=stage3"*) ;; *) echo "Stage 3 runtime export missing" >&2; exit 1 ;; esac
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  ECODA_RUNTIME_MODE=apptainer PREPROCESS_SUBMITTER_TEST=1 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
  --selection-file "${SELECTION}" --exact-batch-selection >/dev/null 2>&1; then
  echo "Stage 3 accepted missing immutable runtime image" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
: > "${CAPTURE}"
BAD="${TMP_DIR}/bad.tsv"
for bad_kind in legacy corrected missing; do
  if [[ "${bad_kind}" == missing ]]; then
    sed '12d' "${SELECTION}" > "${BAD}"
  else
    bad_view="batch_effect_corrected"
    [[ "${bad_kind}" == legacy ]] && bad_view="batch_effect_analysis"
    sed "1s/batch_effect_uncorrected/${bad_view}/" "${SELECTION}" > "${BAD}"
  fi
  RUNS_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs"
  BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
  set +e
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
    PREPROCESS_SUBMITTER_TEST=1 bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
    --selection-file "${BAD}" --exact-batch-selection >/dev/null 2>&1
  RC=$?
  set -e
  [[ ${RC} -ne 0 ]]
  [[ ! -s "${CAPTURE}" ]]
  [[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]
done
BAD_NONEXACT="${TMP_DIR}/bad-nonexact.tsv"
printf 'Adams\tmissing_view\n' > "${BAD_NONEXACT}"
: > "${CAPTURE}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
set +e
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
  --selection-file "${BAD_NONEXACT}" >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]

MISSING_HOME="${TMP_DIR}/missing-home"
mkdir -p "${MISSING_HOME}"
: > "${CAPTURE}"
set +e
HOME="${MISSING_HOME}" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
  --datasets Adams --views benchmark_analysis >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
MISSING_CALLS="$(cat "${CAPTURE}")"
case "${MISSING_CALLS}" in *"--array=1-1%1000"*) ;; *) echo "missing-output path did not submit preprocessing work" >&2; exit 1 ;; esac
case "${MISSING_CALLS}" in *"h5ad_preflight"*) echo "missing-output path submitted an unnecessary preflight" >&2; exit 1 ;; esac

echo "preprocessing stage submitter: OK"
