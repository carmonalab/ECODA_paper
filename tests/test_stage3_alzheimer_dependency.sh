#!/bin/bash
# Focused producer boundary, without Slurm/runtime submission or cohort reads.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-alzheimer-dependency.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'rm -rf "${TMP}"' EXIT
export HPC_SCRATCH_DIR="${TMP}/scratch"
mkdir -p "${HPC_SCRATCH_DIR}/Alzheimer/data" "${TMP}/source"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
eval "$(sed -n '/^stage3_validate_alzheimer_input_dependency() {/,/^}/p' \
  "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh")"
stage3_selection_contains_alzheimer() { awk -F '\t' '$1=="Alzheimer" {found=1} END {exit(!found)}' "$1"; }
_ecoda_runtime_manifest_value() { sed -n "s/^${2}=//p" "$1"; }
ecoda_view_input_name() { printf 'SEAAD_Alzheimer_donor_assay.h5ad'; }
ecoda_validate_input_artifact() { [[ "$1" == "${INPUT}" && "$2" == alzheimer_donor_assay && "$3" == producer ]]; }
stage3_record_alzheimer_input_binding() { printf '%s\t%s\n' "${STAGE3_ALZHEIMER_INPUT_VIEW}" "${STAGE3_ALZHEIMER_INPUT_PATH}" > "${TMP}/record"; }
ecoda_derivative_input() {
  if [[ "$2" == batch_effect_uncorrected && "${DIFFERENT_SECOND_INPUT:-0}" == 1 ]]; then
    printf '%s' "${TMP}/other.h5ad"
  else
    printf '%s' "${INPUT}"
  fi
}
INPUT="${HPC_SCRATCH_DIR}/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad"
printf 'immutable fixture\n' > "${INPUT}"
printf 'different raw input\n' > "${TMP}/other.h5ad"
ln -s "${HPC_SCRATCH_DIR}" "${TMP}/scratch-alias"
PRODUCER_ROOT="${ECODA_RUNS_ROOT}/producer"
mkdir -p "${PRODUCER_ROOT}/manifests" "${PRODUCER_ROOT}/status"
printf 'STAGE=stage2\nSTATE=ACTIVE\n' > "${PRODUCER_ROOT}/metadata"
printf 'STATE=OK\n' > "${PRODUCER_ROOT}/status/terminal"
printf 'alzheimer_donor_assay\tfixture\t%s\t-\t-\n' \
  "${TMP}/scratch-alias/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad" > "${PRODUCER_ROOT}/manifests/steps.tsv"
MANIFEST="${TMP}/selection.tsv"
printf 'Adams\tbenchmark_analysis\nAlzheimer\tbatch_effect_corrected\nAlzheimer\tbatch_effect_uncorrected\n' > "${MANIFEST}"
DATASETS_JSON_FILE="${TMP}/datasets.json"
printf '{"Alzheimer":{"columns":{"sample":"donor_id_assay"}}}\n' > "${DATASETS_JSON_FILE}"
ECODA_SOURCE_SNAPSHOT_REQUIRED=1
ECODA_SOURCE_ROOT="${TMP}/source"
ECODA_SOURCE_MANIFEST="${TMP}/source.manifest"
STAGE3_ALZHEIMER_INPUT_PRODUCER=alzheimer_donor_assay
STAGE3_ALZHEIMER_INPUT_NAME=SEAAD_Alzheimer_donor_assay.h5ad
STAGE3_ALZHEIMER_SAMPLE_COLUMN=donor_id_assay
if stage3_validate_alzheimer_input_dependency > "${TMP}/missing.log" 2>&1; then
  echo 'missing selected Alzheimer producer was accepted' >&2; exit 1
fi
STAGE3_INPUT_PRODUCER_RUN_ID=producer
ECODA_ARTIFACT_BINDING=fixture
stage3_validate_alzheimer_input_dependency
[[ "$(cat "${TMP}/record")" == "$(printf 'batch_effect_corrected\t%s' "${INPUT}")" ]]
DIFFERENT_SECOND_INPUT=1
if stage3_validate_alzheimer_input_dependency > "${TMP}/different.log" 2>&1; then
  echo 'second Alzheimer view escaped its producer proof' >&2; exit 1
fi
unset STAGE3_INPUT_PRODUCER_RUN_ID
printf 'Adams\tbenchmark_analysis\n' > "${MANIFEST}"
stage3_validate_alzheimer_input_dependency
printf 'PASS: selected Alzheimer producer required, mixed selection and canonical alias checked\n'
