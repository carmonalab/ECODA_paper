#!/bin/bash
# Submit exactly one uncorrected batch-candidate evidence job and own its
# terminal validation/synchronization. The durable ECODA gate wraps this
# command; this wrapper itself waits for the one evidence scheduler job.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../src/slurm_config.sh"
export ECODA_GATE_STAGE=evidence
source "${SCRIPT_DIR}/../../src/utils/bash/ecoda_run_common.sh"
if [[ "${EVIDENCE_SUBMITTER_TEST_SYNC:-0}" == "1" ]]; then
  [[ -n "${EVIDENCE_TEST_NAS_TARGET_DIR:-}" ]] || {
    echo "ERROR: evidence sync test requires EVIDENCE_TEST_NAS_TARGET_DIR." >&2
    exit 1
  }
  NAS_TARGET_DIR="${EVIDENCE_TEST_NAS_TARGET_DIR}"
  export NAS_TARGET_DIR
fi

cd "${PROJECT_ROOT}"

SELECTION_FILE_ARG=""
ANALYSIS_ROOT_ARG=""
INPUT_ROOT_ARG=""
OUTPUT_DIR_ARG=""
RUN_ID_ARG=""
PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"

usage() {
  cat <<'EOF'
Usage: submit_batch_candidate_evidence.sh --selection-file TSV \
  --analysis-root DIR --input-root DIR --output-dir DIR --run-id ID
EOF
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --selection-file) SELECTION_FILE_ARG="${2:-}"; shift 2 ;;
    --selection-file=*) SELECTION_FILE_ARG="${1#*=}"; shift ;;
    --analysis-root) ANALYSIS_ROOT_ARG="${2:-}"; shift 2 ;;
    --analysis-root=*) ANALYSIS_ROOT_ARG="${1#*=}"; shift ;;
    --input-root) INPUT_ROOT_ARG="${2:-}"; shift 2 ;;
    --input-root=*) INPUT_ROOT_ARG="${1#*=}"; shift ;;
    --output-dir) OUTPUT_DIR_ARG="${2:-}"; shift 2 ;;
    --output-dir=*) OUTPUT_DIR_ARG="${1#*=}"; shift ;;
    --run-id) RUN_ID_ARG="${2:-}"; shift 2 ;;
    --run-id=*) RUN_ID_ARG="${1#*=}"; shift ;;
    --partition) PARTITION="${2:-}"; shift 2 ;;
    --partition=*) PARTITION="${1#*=}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "${SELECTION_FILE_ARG}" && -r "${SELECTION_FILE_ARG}" ]] || {
  echo "ERROR: --selection-file is required and unreadable." >&2
  exit 1
}
[[ -n "${ANALYSIS_ROOT_ARG}" && -n "${INPUT_ROOT_ARG}" &&
   -n "${OUTPUT_DIR_ARG}" && -n "${RUN_ID_ARG}" ]] || {
  echo "ERROR: --analysis-root, --input-root, --output-dir, and --run-id are required." >&2
  exit 1
}
RUN_ID="${RUN_ID_ARG}"
ecoda_validate_run_id "${RUN_ID}" || exit 1
ecoda_validate_exact_batch_selection "${SELECTION_FILE_ARG}" 2 || exit 1
ecoda_validate_checksum "${SELECTION_FILE_ARG}" || {
  echo "ERROR: evidence selection checksum is missing or invalid." >&2
  exit 1
}

ecoda_init_run evidence "${RUN_ID}" >/dev/null
EVIDENCE_TERMINAL=0
EVIDENCE_FAILURE_REASON=""
evidence_fail() {
  local reason="$1"
  EVIDENCE_FAILURE_REASON="${reason}"
  EVIDENCE_TERMINAL=1
  ecoda_set_run_state FAIL "${reason}" || true
  echo "ERROR: ${reason}" >&2
  exit 1
}
evidence_exit_trap() {
  local rc="$?"
  if [[ ${rc} -ne 0 && ${EVIDENCE_TERMINAL} -eq 0 ]]; then
    ecoda_set_run_state FAIL "${EVIDENCE_FAILURE_REASON:-evidence wrapper aborted}" || true
  fi
  exit "${rc}"
}
trap evidence_exit_trap EXIT

RUN_SELECTION="${ECODA_RUN_ROOT}/manifests/selection.tsv"
if ! cp "${SELECTION_FILE_ARG}" "${RUN_SELECTION}.build.$$"; then
  evidence_fail "failed to copy evidence selection"
fi
if ! ecoda_validate_exact_batch_selection "${RUN_SELECTION}.build.$$" 2; then
  rm -f "${RUN_SELECTION}.build.$$"
  evidence_fail "immutable evidence selection validation failed"
fi
if ! mv -f "${RUN_SELECTION}.build.$$" "${RUN_SELECTION}"; then
  evidence_fail "failed to install immutable evidence selection"
fi
if ! ecoda_write_checksum "${RUN_SELECTION}"; then
  evidence_fail "failed to checksum immutable evidence selection"
fi

EVIDENCE_DATASETS=(
  Alzheimer Breast_cancer Covid19_PBMC Kidney_KPMP Myocardial_infarction
  Diabetes Lupus_PBMC Lung Parkinson Joanito Stephenson CombinedPBMC
)
EVIDENCE_CSVS=()
for ds in "${EVIDENCE_DATASETS[@]}"; do
  EVIDENCE_CSVS+=("${OUTPUT_DIR_ARG}/${ds}_batch_candidate_evidence.csv")
done
REVIEW_CSV="${OUTPUT_DIR_ARG}/batch_candidate_review.csv"
CHECKSUMS="${OUTPUT_DIR_ARG}/checksums.md5"

if ! mkdir -p "${OUTPUT_DIR_ARG}" "${ECODA_RUN_ROOT}/logs"; then
  evidence_fail "failed to create evidence output/log directories"
fi
JOB_WRAP="${PIXI_RSCRIPT} ${SCRIPT_DIR}/build_batch_candidate_evidence.R --selection-file ${RUN_SELECTION} --analysis-root ${ANALYSIS_ROOT_ARG} --input-root ${INPUT_ROOT_ARG} --output-dir ${OUTPUT_DIR_ARG} --config ${DATASETS_JSON_FILE}"
set +e
JOB_ID="$(sbatch --parsable --wait --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=8G \
  --time="${EVIDENCE_TIME_LIMIT:-12:00:00}" \
  --output="${ECODA_RUN_ROOT}/logs/evidence_%j.log" \
  --error="${ECODA_RUN_ROOT}/logs/evidence_%j.err" \
  --mail-user="${USER_EMAIL}" --export="ALL,EVIDENCE_RUN_ROOT=${ECODA_RUN_ROOT}" \
  --wrap="${JOB_WRAP}")"
JOB_RC=$?
set -e
JOB_ID="${JOB_ID%%;*}"
if [[ "${JOB_ID}" =~ ^[0-9]+$ ]]; then
  if ! ecoda_atomic_write "${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv" \
      "EVIDENCE	${JOB_ID}\n"; then
    evidence_fail "failed to persist evidence scheduler ID"
  fi
  echo "BATCH_EFFECT_EVIDENCE_SCHEDULER_ID=${JOB_ID}"
  echo "EVIDENCE_SCHEDULER_ID=${JOB_ID}"
else
  evidence_fail "evidence scheduler returned an invalid job ID: ${JOB_ID}"
fi
[[ ${JOB_RC} -eq 0 ]] || evidence_fail "evidence scheduler submission/wait failed: ${JOB_ID}"

if [[ "${EVIDENCE_SUBMITTER_TEST:-0}" == "1" ]]; then
  if ! ecoda_set_run_state OK "evidence submitter scheduler stub validated"; then
    evidence_fail "failed to write evidence terminal OK state"
  fi
  EVIDENCE_TERMINAL=1
  echo "EVIDENCE_RUN_ID=${RUN_ID}"
  exit 0
fi

CSV_VALIDATOR="${ECODA_RUN_ROOT}/logs/validate_evidence_csv.py"
if ! cat > "${CSV_VALIDATOR}" <<'PY'
import csv
import math
import sys

EXPECTED = [
    "dataset", "method", "method_available", "method_applicable",
    "method_reason", "artifact", "candidate", "present", "completeness",
    "levels", "samples_per_level", "nmi_biology", "constant_candidate",
    "sample_unique_candidate", "perfect_confounded", "marginal_r2",
    "marginal_p", "joint_r2", "joint_p", "warnings", "marginal_p_holm",
    "joint_p_holm",
]
NUMERIC = {
    "completeness", "levels", "nmi_biology", "marginal_r2", "marginal_p",
    "joint_r2", "joint_p", "marginal_p_holm", "joint_p_holm",
}
BOOLEAN = {
    "method_available", "method_applicable", "present", "constant_candidate",
    "sample_unique_candidate", "perfect_confounded",
}

def fail(message):
    raise SystemExit(message)

path, expected_dataset = sys.argv[1], sys.argv[2]
with open(path, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames != EXPECTED:
        fail("unexpected evidence CSV schema: " + path)
    rows = list(reader)
if not rows:
    fail("evidence CSV has no data rows: " + path)
for row in rows:
    if row["dataset"] != expected_dataset and expected_dataset != "__review__":
        fail("evidence CSV dataset does not match filename: " + path)
    if expected_dataset == "__review__" and row["dataset"] not in {
        "Alzheimer", "Breast_cancer", "Covid19_PBMC", "Kidney_KPMP",
        "Myocardial_infarction", "Diabetes", "Lupus_PBMC", "Lung",
        "Parkinson", "Joanito", "Stephenson", "CombinedPBMC",
    }:
        fail("review CSV contains an unknown dataset: " + path)
    if not row["method"].strip() or not row["candidate"].strip():
        fail("evidence CSV has blank method/candidate: " + path)
    for column in BOOLEAN:
        if row[column] not in {"TRUE", "FALSE", ""}:
            fail("evidence CSV has invalid boolean " + column + ": " + path)
    for column in NUMERIC:
        if row[column] == "":
            continue
        try:
            value = float(row[column])
        except ValueError:
            fail("evidence CSV has nonnumeric " + column + ": " + path)
        if not math.isfinite(value):
            fail("evidence CSV has nonfinite " + column + ": " + path)
    if row["completeness"] and not 0.0 <= float(row["completeness"]) <= 1.0:
        fail("evidence CSV completeness is outside [0,1]: " + path)
    if row["levels"] and float(row["levels"]) < 0 or (
        row["levels"] and float(row["levels"]) != int(float(row["levels"]))
    ):
        fail("evidence CSV levels is invalid: " + path)
PY
then
  evidence_fail "failed to install evidence CSV validator"
fi

evidence_allowed_name() {
  local name="$1" ds
  [[ "${name}" == "batch_candidate_review.csv" ||
     "${name}" == "batch_candidate_review.csv.md5" ||
     "${name}" == "checksums.md5" ||
     "${name}" == "checksums.md5.md5" ]] && return 0
  for ds in "${EVIDENCE_DATASETS[@]}"; do
    [[ "${name}" == "${ds}_batch_candidate_evidence.csv" ||
       "${name}" == "${ds}_batch_candidate_evidence.csv.md5" ]] && return 0
  done
  return 1
}

validate_csv() {
  local path="$1" expected_dataset="$2"
  [[ -s "${path}" ]] || return 1
  "${PYTHON_BIN}" "${CSV_VALIDATOR}" "${path}" "${expected_dataset}"
}

validate_outputs() {
  local ds path
  for ds_index in "${!EVIDENCE_DATASETS[@]}"; do
    ds="${EVIDENCE_DATASETS[${ds_index}]}"
    path="${EVIDENCE_CSVS[${ds_index}]}"
    validate_csv "${path}" "${ds}" || return 1
    ecoda_validate_checksum "${path}" || return 1
  done
  validate_csv "${REVIEW_CSV}" "__review__" || return 1
  ecoda_validate_checksum "${REVIEW_CSV}" || return 1
}

for path in "${OUTPUT_DIR_ARG}"/*; do
  [[ -e "${path}" ]] || continue
  evidence_allowed_name "$(basename "${path}")" ||
    evidence_fail "unexpected evidence output: ${path}"
done
validate_outputs || evidence_fail "evidence outputs failed twelve-cohort validation"

CHECKSUMS_TMP="${CHECKSUMS}.tmp.$$"
: > "${CHECKSUMS_TMP}" || evidence_fail "failed to create evidence checksum manifest"
for path in "${EVIDENCE_CSVS[@]}" "${REVIEW_CSV}"; do
  digest="$(ecoda_md5_file "${path}")" ||
    evidence_fail "failed to checksum evidence output: ${path}"
  printf '%s  %s\n' "${digest}" "${path}" >> "${CHECKSUMS_TMP}" ||
    evidence_fail "failed to write evidence checksum manifest"
done
if ! mv -f "${CHECKSUMS_TMP}" "${CHECKSUMS}" ||
   ! ecoda_write_checksum "${CHECKSUMS}" ||
   ! ecoda_validate_checksum "${CHECKSUMS}"; then
  rm -f "${CHECKSUMS_TMP}"
  evidence_fail "evidence checksum manifest validation failed"
fi

REMOTE_DIR="${NAS_TARGET_DIR}/batch_effect/uncorrected/evidence/${RUN_ID}"
[[ -d "${NAS_TARGET_DIR}" ]] ||
  evidence_fail "NAS evidence target is unreachable"
mkdir -p "${REMOTE_DIR}" ||
  evidence_fail "failed to create remote evidence destination"
for path in "${EVIDENCE_CSVS[@]}" "${REVIEW_CSV}" "${CHECKSUMS}"; do
  ecoda_validate_checksum "${path}" ||
    evidence_fail "local evidence checksum validation failed: ${path}"
  rsync -rlptDv "${path}" "${path}.md5" "${REMOTE_DIR}/" ||
    evidence_fail "evidence rsync failed: ${path}"
done
for path in "${EVIDENCE_CSVS[@]}" "${REVIEW_CSV}" "${CHECKSUMS}"; do
  remote_path="${REMOTE_DIR}/$(basename "${path}")"
  ecoda_validate_checksum_remote "${remote_path}" "${remote_path}.md5" ||
    evidence_fail "remote evidence checksum validation failed: ${remote_path}"
  cmp -s "${path}.md5" "${remote_path}.md5" ||
    evidence_fail "remote evidence sidecar differs: ${remote_path}"
done
for path in "${REMOTE_DIR}"/*; do
  [[ -e "${path}" ]] || continue
  evidence_allowed_name "$(basename "${path}")" ||
    evidence_fail "unexpected remote evidence output: ${path}"
done
if ! ecoda_set_run_state OK "twelve-cohort evidence validated and synchronized"; then
  evidence_fail "failed to write evidence terminal OK state"
fi
EVIDENCE_TERMINAL=1
echo "EVIDENCE_RUN_ID=${RUN_ID}"
