#!/bin/bash
# One counts-only diagnostic job. Invoke through the durable ecoda-benchmark gate.
set -euo pipefail
if [[ -n "${ECODA_SOURCE_ROOT:-}" ]]; then
  SCRIPT_DIR="${ECODA_SOURCE_ROOT}/src/utils/bash"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
SELF="${PROJECT_ROOT}/src/utils/bash/zeroimp_count_screen.sh"
RUNNER="${PROJECT_ROOT}/src/notebooks/run_zeroimp_count_screen.R"
PURPOSE=counts_only_zeroimp_screen
PRODUCER=zeroimp_count_screen
MODE=controller
if [[ "${1:-}" == --worker || "${1:-}" == --validate ]]; then MODE="${1#--}"; shift; fi
CONTROL=""; CONTROL_MD5=""; DATASETS=""
while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || { echo 'ERROR: missing zeroimp count-screen argument value' >&2; exit 1; }
  case "$1" in
    --analysis_release) CONTROL="$2" ;;
    --analysis_release_md5) CONTROL_MD5="$2" ;;
    --datasets) DATASETS="$2" ;;
    *) echo "ERROR: unknown zeroimp count-screen argument: $1" >&2; exit 1 ;;
  esac
  shift 2
done
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && -n "${ECODA_RUN_ID:-}" ]] || {
  echo 'ERROR: zeroimp count-screen requires an immutable snapshot and run ID' >&2; exit 1;
}
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
ecoda_require_source_script_path "${RUNNER}" "${ECODA_SOURCE_ROOT}" >/dev/null
if [[ "${MODE}" == controller ]]; then
  [[ -n "${CONTROL}" && "${CONTROL_MD5}" =~ ^[a-f0-9]{32}$ && -n "${DATASETS}" ]] || {
    echo 'ERROR: explicit checked --analysis_release and --datasets are required' >&2; exit 1;
  }
  [[ ! -e "${ECODA_RUNS_ROOT}/${ECODA_RUN_ID}" && ! -L "${ECODA_RUNS_ROOT}/${ECODA_RUN_ID}" ]] || { echo 'ERROR: zeroimp run root is create-only' >&2; exit 1; }
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  ecoda_atomic_write "${ECODA_RUN_ROOT}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
  trap 'ecoda_set_run_state FAIL zeroimp_count_screen_binding_failed || true' EXIT
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
  export ECODA_RUNTIME_PROFILE=stage5
  ecoda_runtime_validate_submission apptainer
else
  [[ -z "${CONTROL}${CONTROL_MD5}${DATASETS}" ]] || {
    echo 'ERROR: worker/validate scope comes only from the frozen run' >&2; exit 1;
  }
  ecoda_open_run "${ECODA_RUN_ID}" stage5
fi
RUN="${ECODA_RUN_ROOT}"
export ECODA_RUN_ROOT ECODA_RUN_ID
MANIFESTS="${RUN}/manifests"
OUTPUT="${RUN}/artifacts/zeroimp_count_screen/analysis"
STATUS="${RUN}/status/zeroimp_count_screen.status"
SCOPE_FILES=(source.manifest runtime.identity operation selection.tsv analysis_release.rds inputs.tsv checkpoint_reuse.json screen_identity.json)
OUTPUT_FILES=()

if [[ "${MODE}" == controller ]]; then
  # Read the reviewed control with configured R; inspect compact HR counts only.
  # Freeze every selected artifact, receipt and review evidence before submission.
  ${PIXI_RSCRIPT} --vanilla "${RUNNER}" --project_root "${PROJECT_ROOT}" \
    --analysis_release "${CONTROL}" --analysis_release_md5 "${CONTROL_MD5}" \
    --datasets "${DATASETS}" --config_path "${DATASETS_JSON_FILE}" \
    --run_id "${ECODA_RUN_ID}" --bind "${MANIFESTS}"
  ecoda_atomic_write "${MANIFESTS}/operation" "zeroimp-count-screen\n"
  ecoda_atomic_copy "${CONTROL}" "${MANIFESTS}/analysis_release.rds"
  "${PYTHON_BIN}" - "${RUN}" "${PURPOSE}" "${CONTROL_MD5}" "${DATASETS_JSON_FILE}" <<'PY'
import hashlib,json,pathlib,sys
run=pathlib.Path(sys.argv[1]);m=run/'manifests'
names=['source.manifest','runtime.identity','operation','selection.tsv','analysis_release.rds','inputs.tsv','checkpoint_reuse.json']
identity={'schema':1,'purpose':sys.argv[2],'run_id':run.name,'control_md5':sys.argv[3],
          'config_path':sys.argv[4],'config_sha256':hashlib.sha256(pathlib.Path(sys.argv[4]).read_bytes()).hexdigest(),
          'identities':{n:hashlib.sha256((m/n).read_bytes()).hexdigest() for n in names}}
if hashlib.md5((m/'analysis_release.rds').read_bytes()).hexdigest()!=identity['control_md5']: raise SystemExit('ERROR: control changed during binding')
(m/'screen_identity.json').write_text(json.dumps(identity,sort_keys=True)+'\n')
PY
  for NAME in "${SCOPE_FILES[@]}"; do
    ecoda_write_checksum "${MANIFESTS}/${NAME}"
    ecoda_write_artifact_record "${MANIFESTS}/${NAME}" zeroimp_scope "${ECODA_RUN_ID}" >/dev/null
  done
  mkdir -p "${OUTPUT%/*}"
fi

ecoda_validate_artifact_record "${MANIFESTS}/selection.tsv" zeroimp_scope "${ECODA_RUN_ID}" >/dev/null
OUTPUT_FILES=(screen_manifest.json status.json)
while IFS=$'\t' read -r DS VIEW METHOD; do
  [[ "${DS}" =~ ^[A-Za-z0-9_]+$ && "${VIEW}" == benchmark_analysis && "${METHOD}" == zeroimp_count_screen ]] || exit 1
  OUTPUT_FILES+=("${DS}_zeroimp_counts.rds" "${DS}_zeroimp_counts.csv")
done < "${MANIFESTS}/selection.tsv"

validate_scope() {
  local NAME PATH_IN MD5_IN OWNER PARENT RECORD BIRTH_PRODUCER
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" ]] || return 1
  [[ "$(cat "${MANIFESTS}/operation")" == zeroimp-count-screen ]] || return 1
  ecoda_runtime_validate_bound_run || return 1
  for NAME in "${SCOPE_FILES[@]}"; do
    ecoda_validate_artifact_record "${MANIFESTS}/${NAME}" zeroimp_scope "${ECODA_RUN_ID}" >/dev/null || return 1
  done
  "${PYTHON_BIN}" - "${RUN}" "${PURPOSE}" "${HPC_SCRATCH_DIR}" "${DATASETS_JSON_FILE}" <<'PY' || return 1
import hashlib,json,pathlib,re,sys
run=pathlib.Path(sys.argv[1]);m=run/'manifests';x=json.loads((m/'screen_identity.json').read_text())
expected={'source.manifest','runtime.identity','operation','selection.tsv','analysis_release.rds','inputs.tsv','checkpoint_reuse.json'}
if x['schema']!=1 or x['purpose']!=sys.argv[2] or x['run_id']!=run.name or set(x['identities'])!=expected or x['config_path']!=sys.argv[4]: raise SystemExit('ERROR: zeroimp count-screen identity differs')
if hashlib.sha256(pathlib.Path(sys.argv[4]).read_bytes()).hexdigest()!=x['config_sha256']: raise SystemExit('ERROR: config changed')
for n,sha in x['identities'].items():
    p=m/n
    if p.resolve()!=p or hashlib.sha256(p.read_bytes()).hexdigest()!=sha: raise SystemExit('ERROR: frozen scope changed: '+n)
rows=[r.split('\t') for r in (m/'selection.tsv').read_text().splitlines()]
allowed={'Bassez','Stephenson','Adams','Gongsharma_cmv_young_males','Kfoury','Kim','Lee','Pelka','Smillie','Wu','Zhang'}
if not rows or any(len(r)!=3 or r[1:]!=['benchmark_analysis','zeroimp_count_screen'] or r[0] not in allowed for r in rows) or len({r[0] for r in rows})!=len(rows): raise SystemExit('ERROR: invalid frozen selection')
if hashlib.md5((m/'analysis_release.rds').read_bytes()).hexdigest()!=x['control_md5']: raise SystemExit('ERROR: control checksum differs')
seen=set();scratch=pathlib.Path(sys.argv[3]).resolve()
for row in (m/'inputs.tsv').read_text().splitlines():
    fields=row.split('\t')
    if len(fields)!=2: raise SystemExit('ERROR: invalid input inventory')
    name,md5=fields;p=pathlib.Path(name)
    if name in seen or not p.is_absolute() or p.resolve()!=p or not p.is_file() or scratch not in p.parents or not re.fullmatch('[a-f0-9]{32}',md5): raise SystemExit('ERROR: unsafe selected input: '+name)
    seen.add(name);digest=hashlib.md5()
    with p.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    if digest.hexdigest()!=md5:raise SystemExit('ERROR: selected input changed: '+name)
PY
  # Never fabricate producer provenance for reviewed legacy inputs. A run-owned
  # artifact must additionally have its actual published producer and OK owner.
  while IFS=$'\t' read -r PATH_IN MD5_IN; do
    ecoda_require_input_ownership "${PATH_IN}" "${ECODA_RUN_ID}" || return 1
    case "${PATH_IN}" in "${ECODA_RUNS_ROOT}/"*)
      OWNER="$(ecoda_artifact_owner_dir "${PATH_IN}")" || return 1
      _ecoda_artifact_owner_validate_dir "${OWNER}" "${PATH_IN}" || return 1
      PARENT="${PATH_IN#${ECODA_RUNS_ROOT}/}"; PARENT="${PARENT%%/*}"
      [[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK && "${ECODA_ARTIFACT_OWNER_RUN}" == "${PARENT}" ]] || return 1
      RECORD="$(ecoda_artifact_record_path "${PATH_IN}" "${PARENT}")" || return 1
      BIRTH_PRODUCER="$(sed -n 's/^PRODUCER=//p' "${RECORD}")" || return 1
      ecoda_validate_artifact_record "${PATH_IN}" "${BIRTH_PRODUCER}" "${PARENT}" >/dev/null || return 1
      ;; esac
  done < "${MANIFESTS}/inputs.tsv"
}
screen_r() {
  local ACTION="$1"; shift
  local BOUND_MD5 BOUND_CSV
  BOUND_MD5="$("${PYTHON_BIN}" -c 'import json,sys;print(json.load(open(sys.argv[1]))["control_md5"])' "${MANIFESTS}/screen_identity.json")" || return 1
  BOUND_CSV="$(cut -f1 "${MANIFESTS}/selection.tsv" | paste -sd, -)" || return 1
  export RETICULATE_PYTHON="${PYTHON_BIN}" ECODA_DERIVED_ROOT="${OUTPUT%/*}"
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
  local FLAG=--validate
  [[ "${ACTION}" != produce ]] || FLAG=--produce
  ${PIXI_RSCRIPT} --vanilla "${RUNNER}" --project_root "${PROJECT_ROOT}" \
    --analysis_release "${MANIFESTS}/analysis_release.rds" --analysis_release_md5 "${BOUND_MD5}" \
    --datasets "${BOUND_CSV}" --config_path "${DATASETS_JSON_FILE}" \
    --run_id "${ECODA_RUN_ID}" "${FLAG}" "${OUTPUT}"

}
validate_report() {
  local NAME OWNER
  validate_scope || return 1
  [[ -s "${STATUS}" && "$(sed -n 's/^STATE=//p' "${STATUS}")" == OK &&
     "$(sed -n 's/^RUN_ID=//p' "${STATUS}")" == "${ECODA_RUN_ID}" ]] || return 1
  for NAME in "${OUTPUT_FILES[@]}"; do
    ecoda_validate_artifact_record "${OUTPUT}/${NAME}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
    OWNER="$(ecoda_artifact_owner_dir "${OUTPUT}/${NAME}")" || return 1
    _ecoda_artifact_owner_validate_dir "${OWNER}" "${OUTPUT}/${NAME}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  done
  screen_r validate || return 1
}
fail_active_outputs() {
  local NAME OWNER
  for NAME in "${OUTPUT_FILES[@]}"; do
    OWNER="$(ecoda_artifact_owner_dir "${OUTPUT}/${NAME}" 2>/dev/null)" || continue
    if _ecoda_artifact_owner_validate_dir "${OWNER}" "${OUTPUT}/${NAME}" &&
       [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]; then
      ecoda_artifact_owner_set_state "${OUTPUT}/${NAME}" FAIL zeroimp_count_screen_failed >/dev/null || true
    fi
  done
}
if [[ "${MODE}" == validate ]]; then
  validate_report || exit 1
  [[ "$(sed -n 's/^STATE=//p' "${RUN}/status/terminal")" == OK && "$(sed -n 's/^RUN_ID=//p' "${RUN}/status/terminal")" == "${ECODA_RUN_ID}" ]] || exit 1
  for NAME in scheduler_id scheduler_ids.tsv sbatch_submission.txt; do
    ecoda_validate_artifact_record "${MANIFESTS}/${NAME}" zeroimp_submission "${ECODA_RUN_ID}" >/dev/null || exit 1
  done
  [[ "$(cat "${MANIFESTS}/scheduler_ids.tsv")" == "$(printf 'ARRAY\t%s' "$(cat "${MANIFESTS}/scheduler_id")")" &&
     "$(cat "${MANIFESTS}/scheduler_id")" =~ ^[0-9]+$ ]] || exit 1
  echo "ZEROIMP_COUNT_SCREEN_VALIDATED=${ECODA_RUN_ID}"
  exit 0
fi
if [[ "${MODE}" == worker ]]; then
  failed() {
    ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\n" || true
    fail_active_outputs
  }
  trap failed EXIT
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker
  cd "${PROJECT_ROOT}"
  validate_scope || exit 1
  [[ ! -e "${OUTPUT}" && ! -L "${OUTPUT}" ]] || { echo 'ERROR: zeroimp count-screen output is create-only' >&2; exit 1; }
  screen_r produce
  validate_scope || exit 1
  screen_r validate
  for NAME in "${OUTPUT_FILES[@]}"; do
    # The runner writes the scientific pair's truthful birth-path sidecars.
    case "${NAME}" in *.rds|*.csv) ecoda_validate_checksum "${OUTPUT}/${NAME}" ;;
      *) ecoda_write_checksum "${OUTPUT}/${NAME}" ;; esac
    ecoda_write_artifact_record "${OUTPUT}/${NAME}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
    ecoda_artifact_owner_set_state "${OUTPUT}/${NAME}" OK zeroimp_count_screen_complete >/dev/null
  done
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\n"
  validate_report || exit 1
  trap - EXIT
  exit 0
fi
controller_failed() {
  ecoda_set_run_state FAIL zeroimp_count_screen_failed || true
  fail_active_outputs
}
trap controller_failed EXIT
for NAME in "${OUTPUT_FILES[@]}"; do
  ecoda_artifact_owner_acquire "${OUTPUT}/${NAME}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
done
validate_scope || exit 1
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)"
set +e
JOB="$(sbatch --parsable --wait --partition="${SLURM_PARTITION}" --ntasks=1 --cpus-per-task=1 \
  --mem=64G --time=02:00:00 --output="${RUN}/logs/zeroimp_count_screen_%j.log" --error="${RUN}/logs/zeroimp_count_screen_%j.err" \
  --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" "${SELF}" --worker)"
RC=$?
set -e
ecoda_atomic_write "${MANIFESTS}/sbatch_submission.txt" "${JOB}\n"
JOB="${JOB%%;*}"
[[ "${JOB}" =~ ^[0-9]+$ ]] || { echo 'ERROR: ambiguous zeroimp count-screen submission; preserve run and inspect it' >&2; exit 1; }
printf 'ZEROIMP_COUNT_SCREEN_JOB_ID=%s\n' "${JOB}"
ecoda_atomic_write "${MANIFESTS}/scheduler_id" "${JOB}\n"
ecoda_atomic_write "${MANIFESTS}/scheduler_ids.tsv" "ARRAY\t${JOB}\n"
for NAME in scheduler_id scheduler_ids.tsv sbatch_submission.txt; do
  ecoda_write_checksum "${MANIFESTS}/${NAME}"
  ecoda_write_artifact_record "${MANIFESTS}/${NAME}" zeroimp_submission "${ECODA_RUN_ID}" >/dev/null
done
[[ ${RC} -eq 0 ]] || exit "${RC}"
validate_report || exit 1
ecoda_set_run_state OK zeroimp_count_screen_complete
trap - EXIT
