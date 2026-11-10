#!/bin/bash
# Fixed scalar Breast UNC pilot; immutable executor and durable ecoda-benchmark only.
set -euo pipefail
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && "${ECODA_SOURCE_ROOT:-}" == /* ]] || exit 1
SCRIPT_DIR="${ECODA_SOURCE_ROOT%/}/src/utils/bash"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
SELF="${SCRIPT_DIR}/singlecell_batch_geometry.sh"
RUNNER="${PROJECT_ROOT}/src/notebooks/run_singlecell_batch_geometry.py"
PURPOSE=singlecell_breast_unc_pilot
PRODUCER=singlecell_batch_geometry
MODE=controller
if [[ "${1:-}" == --worker || "${1:-}" == --validate ]]; then MODE="${1#--}"; shift; fi
CONTROL=""; CONTROL_SHA=""
while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || exit 1
  case "$1" in
    --control) [[ -z "${CONTROL}" ]] || exit 1; CONTROL="$2" ;;
    --control-sha256) [[ -z "${CONTROL_SHA}" ]] || exit 1; CONTROL_SHA="$2" ;;
    *) echo 'ERROR: unknown fixed-pilot argument' >&2; exit 1 ;;
  esac
  shift 2
done
ecoda_validate_run_id "${ECODA_RUN_ID:-}"
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
ecoda_require_source_script_path "${RUNNER}" "${ECODA_SOURCE_ROOT}" >/dev/null
SCRATCH="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
RUN="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${RUN}" ]] || exit 1
if [[ "${MODE}" == controller ]]; then
  [[ ! -e "${RUN}" && ! -L "${RUN}" && "${CONTROL}" == "${SCRATCH}/_ecoda_inputs/"* &&
     "${CONTROL_SHA}" =~ ^[a-f0-9]{64}$ && -f "${CONTROL}" && ! -L "${CONTROL}" &&
     "$(ecoda_realpath_existing "${CONTROL}")" == "${CONTROL}" &&
     "$(ecoda_sha256_file "${CONTROL}")" == "${CONTROL_SHA}" ]] || exit 1
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  ecoda_atomic_write "${RUN}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
else
  [[ -z "${CONTROL}${CONTROL_SHA}" ]] || exit 1
  ecoda_open_run "${ECODA_RUN_ID}" stage5
fi
export ECODA_RUN_ROOT="${RUN}" ECODA_RUNTIME_PROFILE=stage5 PYTHONDONTWRITEBYTECODE=1
M="${RUN}/manifests"; OUTPUT="${RUN}/artifacts/singlecell_batch_geometry"; STATUS="${RUN}/status/singlecell_pilot.status"
SCOPE=(source.manifest runtime.identity operation selection.tsv control.json control.sha256 inputs.tsv)
OWNED=()
failed() {
  local path owner
  ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\n" || true
  for path in "${OWNED[@]:-}" "${OUTPUT}"; do
    [[ -n "${path}" ]] || continue
    owner="$(ecoda_artifact_owner_dir "${path}")" || continue
    if _ecoda_artifact_owner_validate_dir "${owner}" "${path}" &&
       [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]; then
      ecoda_artifact_owner_set_state "${path}" FAIL singlecell_pilot_failed >/dev/null || true
    fi
  done
  if [[ "${MODE}" == controller ]]; then ecoda_set_run_state FAIL singlecell_pilot_failed || true; fi
}
if [[ "${MODE}" != validate ]]; then trap failed EXIT; fi
if [[ "${MODE}" == controller ]]; then
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest"
  ecoda_runtime_validate_submission apptainer
  ecoda_atomic_copy "${CONTROL}" "${M}/control.json"
  ecoda_atomic_write "${M}/control.sha256" "${CONTROL_SHA}\n"
  ecoda_atomic_write "${M}/operation" "singlecell-breast-unc-pilot\n"
  ecoda_atomic_write "${M}/selection.tsv" "Breast_cancer\tbatch_effect_uncorrected\tunc_raw_pilot\n"
  # Input inventory contains metadata/provenance references; parent H5 is NEVER rehashed here.
  "${PYTHON_BIN}" -B - "${M}" <<'PY'
import hashlib,json,pathlib,sys
m=pathlib.Path(sys.argv[1]);c=json.loads((m/'control.json').read_text())
paths=[m/'control.json',pathlib.Path(c['stage3_receipt']['path']),pathlib.Path(c['stage3_review']['path'])]
with (m/'inputs.tsv').open('x') as f:
 for p in paths:f.write(str(p)+'\t'+hashlib.md5(p.read_bytes()).hexdigest()+'\n')
PY
  for NAME in "${SCOPE[@]}"; do
    ecoda_write_checksum "${M}/${NAME}"
    ecoda_write_artifact_record "${M}/${NAME}" singlecell_pilot_scope "${ECODA_RUN_ID}" >/dev/null
  done
fi
pilot_python() {
  "${PYTHON_BIN}" -B "${RUNNER}" --control "${M}/control.json" --control-sha256 "$(cat "${M}/control.sha256")" "$@"
}
validate_scope() {
  local name
  ecoda_runtime_validate_bound_run || return 1
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest" || return 1
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" &&
     "$(cat "${M}/operation")" == singlecell-breast-unc-pilot &&
     "$(cat "${M}/selection.tsv")" == "$(printf 'Breast_cancer\tbatch_effect_uncorrected\tunc_raw_pilot')" ]] || return 1
  for name in "${SCOPE[@]}"; do ecoda_validate_artifact_record "${M}/${name}" singlecell_pilot_scope "${ECODA_RUN_ID}" >/dev/null || return 1; done
  pilot_python --check-control
}
validate_report() {
  local path md5 extra owner
  validate_scope || return 1
  [[ "$(sed -n 's/^STATE=//p' "${STATUS}")" == OK && "$(sed -n 's/^RUN_ID=//p' "${STATUS}")" == "${ECODA_RUN_ID}" ]] || return 1
  ecoda_validate_artifact_record "${OUTPUT}/files.tsv" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
  while IFS=$'\t' read -r path md5 extra; do
    [[ "${path}" == "${OUTPUT}/"* && "${md5}" =~ ^[a-f0-9]{32}$ && -z "${extra}" ]] || return 1
    ecoda_validate_artifact_record "${path}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
    [[ "${ECODA_ARTIFACT_RECORD_MD5}" == "${md5}" ]] || return 1
    owner="$(ecoda_artifact_owner_dir "${path}")" || return 1
    _ecoda_artifact_owner_validate_dir "${owner}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  done < "${OUTPUT}/files.tsv"
  for path in "${OUTPUT}" "${OUTPUT}/files.tsv"; do
    owner="$(ecoda_artifact_owner_dir "${path}")" || return 1; _ecoda_artifact_owner_validate_dir "${owner}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  done
  pilot_python --output "${OUTPUT}" --run-id "${ECODA_RUN_ID}" --validate-output
}
if [[ "${MODE}" == validate ]]; then
  validate_report
  [[ "$(sed -n 's/^STATE=//p' "${RUN}/status/terminal")" == OK &&
     "$(sed -n 's/^RUN_ID=//p' "${RUN}/status/terminal")" == "${ECODA_RUN_ID}" ]] || exit 1
  for NAME in scheduler_ids.tsv sbatch_submission.txt; do ecoda_validate_artifact_record "${M}/${NAME}" singlecell_pilot_submission "${ECODA_RUN_ID}" >/dev/null; done
  [[ "$(wc -l < "${M}/scheduler_ids.tsv" | tr -d '[:space:]')" == 1 ]] || exit 1
  awk -F '\t' 'NF != 2 || $1 != "ARRAY" || $2 !~ /^[1-9][0-9]*$/ {exit 1}' "${M}/scheduler_ids.tsv"
  echo "SINGLECELL_BREAST_UNC_PILOT_VALIDATED=${ECODA_RUN_ID}"; exit 0
fi
if [[ "${MODE}" == worker ]]; then
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker
  cd "${PROJECT_ROOT}"; validate_scope
  [[ ! -e "${OUTPUT}" && ! -L "${OUTPUT}" ]] || exit 1
  OWNER="$(ecoda_artifact_owner_dir "${OUTPUT}")"; _ecoda_artifact_owner_validate_dir "${OWNER}" "${OUTPUT}"
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]] || exit 1
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" NUMBA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
  export NUMBA_CACHE_DIR="${RUN}/cache/numba" MPLCONFIGDIR="${RUN}/cache/matplotlib"; mkdir -p "${NUMBA_CACHE_DIR}" "${MPLCONFIGDIR}"
  pilot_python --output "${OUTPUT}" --run-id "${ECODA_RUN_ID}"
  validate_scope
  while IFS=$'\t' read -r PATH_IN MD5_IN EXTRA; do
    [[ "${PATH_IN}" == "${OUTPUT}/"* && "${MD5_IN}" =~ ^[a-f0-9]{32}$ && -z "${EXTRA}" && "$(ecoda_md5_file "${PATH_IN}")" == "${MD5_IN}" ]] || exit 1
    ecoda_artifact_owner_acquire "${PATH_IN}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null; OWNED+=("${PATH_IN}")
  done < "${OUTPUT}/files.tsv"
  ecoda_artifact_owner_acquire "${OUTPUT}/files.tsv" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null; OWNED+=("${OUTPUT}/files.tsv")
  for PATH_IN in "${OWNED[@]}"; do ecoda_write_checksum "${PATH_IN}"; ecoda_write_artifact_record "${PATH_IN}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null; done
  pilot_python --output "${OUTPUT}" --run-id "${ECODA_RUN_ID}" --validate-output
  for PATH_IN in "${OWNED[@]}" "${OUTPUT}"; do ecoda_artifact_owner_set_state "${PATH_IN}" OK singlecell_pilot_complete >/dev/null; done
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\n"; validate_report; trap - EXIT; exit 0
fi
validate_scope
mkdir -p "${OUTPUT%/*}"; ecoda_artifact_owner_acquire "${OUTPUT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)"
set +e
JOB="$(sbatch --parsable --wait --partition="${SLURM_PARTITION_BENCHMARK_CPU}" --ntasks=1 --cpus-per-task=8 --mem=64G --time=02:00:00 \
  --output="${RUN}/logs/singlecell_pilot_%j.log" --error="${RUN}/logs/singlecell_pilot_%j.err" \
  --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" "${SELF}" --worker)"
RC=$?; set -e
ecoda_atomic_write "${M}/sbatch_submission.txt" "${JOB}\n"; JOB="${JOB%%;*}"
[[ "${JOB}" =~ ^[0-9]+$ ]] || { echo 'ERROR: ambiguous pilot submission; preserve and inspect' >&2; exit 1; }
echo "SINGLECELL_PILOT_JOB_ID=${JOB}"
ecoda_atomic_write "${M}/scheduler_ids.tsv" "ARRAY\t${JOB}\n"
for NAME in scheduler_ids.tsv sbatch_submission.txt; do ecoda_write_checksum "${M}/${NAME}"; ecoda_write_artifact_record "${M}/${NAME}" singlecell_pilot_submission "${ECODA_RUN_ID}" >/dev/null; done
[[ ${RC} -eq 0 ]] || exit "${RC}"
validate_report; ecoda_set_run_state OK singlecell_pilot_complete; trap - EXIT
# NAS publication is a separate existing finite ecoda_nas_backup.sh operation after Root review.
