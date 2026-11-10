#!/bin/bash
# Fixed remaining11 main representations; pilot stays a separate immutable run.
set -euo pipefail
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && "${ECODA_SOURCE_ROOT:-}" == /* ]] || exit 1
SCRIPT_DIR="${ECODA_SOURCE_ROOT%/}/src/utils/bash"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
SELF="${SCRIPT_DIR}/singlecell_batch_geometry_bulk.sh"
RUNNER="${PROJECT_ROOT}/src/notebooks/run_singlecell_batch_geometry.py"
MODE=controller; CONTROL=""; CONTROL_SHA=""; MANIFEST=""; MEMORY=""
case "${1:-}" in --worker|--watchdog|--validate) MODE="${1#--}";shift ;; esac
while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || exit 1
  case "$1" in
    --control) [[ -z "${CONTROL}" ]] || exit 1;CONTROL="$2" ;;
    --control-sha256) [[ -z "${CONTROL_SHA}" ]] || exit 1;CONTROL_SHA="$2" ;;
    --manifest) [[ -z "${MANIFEST}" ]] || exit 1;MANIFEST="$2" ;;
    --memory) [[ -z "${MEMORY}" ]] || exit 1;MEMORY="$2" ;;
    *) exit 1 ;;
  esac;shift 2
done
ecoda_validate_run_id "${ECODA_RUN_ID:-}"
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
ecoda_require_source_script_path "${RUNNER}" "${ECODA_SOURCE_ROOT}" >/dev/null
SCRATCH="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")";RUN="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${RUN}" ]] || exit 1
export ECODA_RUN_ROOT="${RUN}" ECODA_RUNTIME_PROFILE=stage5 PYTHONDONTWRITEBYTECODE=1
M="${RUN}/manifests";PRODUCER=singlecell_batch_geometry;PURPOSE=singlecell_batch_main_bulk
SCOPE=(source.manifest runtime.identity operation selection.tsv control.json control.sha256 inputs.tsv initial_64.tsv initial_128.tsv)
expected_selection() {
  printf '1\tAlzheimer\tunc_raw\t128\n2\tAlzheimer\tcor_harmony\t128\n3\tBreast_cancer\tcor_harmony\t64\n4\tCovid19_PBMC\tunc_raw\t64\n5\tCovid19_PBMC\tcor_harmony\t64\n6\tJoanito\tunc_raw\t64\n7\tJoanito\tcor_harmony\t64\n8\tLupus_PBMC\tunc_raw\t64\n9\tLupus_PBMC\tcor_harmony\t64\n10\tStephenson\tunc_raw\t64\n11\tStephenson\tcor_harmony\t64\n'
}
publish() { ecoda_write_checksum "$1";ecoda_write_artifact_record "$1" "$2" "${ECODA_RUN_ID}" >/dev/null; }
python_run() { "${PYTHON_BIN}" -B "${RUNNER}" --control "${M}/control.json" --control-sha256 "$(cat "${M}/control.sha256")" "$@"; }
scope_validate() {
  local n
  ecoda_runtime_validate_bound_run || return 1
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest" || return 1
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" && "$(cat "${M}/operation")" == singlecell-batch-main-bulk && "$(cat "${M}/selection.tsv")" == "$(expected_selection)" ]] || return 1
  for n in "${SCOPE[@]}";do ecoda_validate_artifact_record "${M}/${n}" singlecell_bulk_scope "${ECODA_RUN_ID}" >/dev/null || return 1;done
  python_run --check-control
}
manifest_validate() {
  ecoda_validate_run_owned_path "$1" "${RUN}" && ecoda_validate_artifact_record "$1" singlecell_bulk_scope "${ECODA_RUN_ID}" >/dev/null &&
    ecoda_validate_manifest "$1" 4 && awk -F '\t' 'NR==FNR {a[$0]=1;next} !($0 in a)||seen[$0]++ {exit 1}' "${M}/selection.tsv" "$1"
}
owner_ok() {
  local owner;owner="$(ecoda_artifact_owner_dir "$1")" || return 1
  _ecoda_artifact_owner_validate_dir "${owner}" "$1" && [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]]
}
row_paths() {
  ROW="$1";DS="$2";REP="$3";INITIAL="$4";MEM="$5"
  [[ "${INITIAL}:${MEM}" == 64:64 || "${INITIAL}:${MEM}" == 64:128 || "${INITIAL}:${MEM}" == 64:256 || "${INITIAL}:${MEM}" == 128:128 || "${INITIAL}:${MEM}" == 128:256 || "${INITIAL}:${MEM}" == 128:512 ]] || return 1
  OUTPUT="${RUN}/artifacts/singlecell_batch_geometry/${DS}/${REP}/attempt_${MEM}"
  STATUS="${RUN}/status/row_${ROW}.status"
}
row_validate() {
  local p md5 extra
  ecoda_validate_artifact_record "${STATUS}" singlecell_bulk_row "${ECODA_RUN_ID}" >/dev/null && owner_ok "${STATUS}" || return 1
  [[ "$(cat "${STATUS}")" == "$(printf 'STATE=OK\nRUN_ID=%s\nROW=%s\nMEMORY=%s\nOUTPUT=%s\n' "${ECODA_RUN_ID}" "${ROW}" "${MEM}" "${OUTPUT}")" ]] || return 1
  ecoda_validate_artifact_record "${OUTPUT}/files.tsv" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null && owner_ok "${OUTPUT}/files.tsv" && owner_ok "${OUTPUT}" || return 1
  [[ "$(cut -f1 "${OUTPUT}/files.tsv")" == "$(printf '%s\n' "${OUTPUT}/embedding.h5ad" "${OUTPUT}/metrics.json" "${OUTPUT}/display.feather" "${OUTPUT}/phases.json")" ]] || return 1
  while IFS=$'\t' read -r p md5 extra;do
    [[ "${md5}" =~ ^[a-f0-9]{32}$ && -z "${extra}" ]] && ecoda_validate_artifact_record "${p}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null && [[ "${ECODA_ARTIFACT_RECORD_MD5}" == "${md5}" ]] && owner_ok "${p}" || return 1
  done < "${OUTPUT}/files.tsv"
}
all_rows_validate() {
  local row ds rep initial mem
  while IFS=$'\t' read -r row ds rep initial;do
    mem="$(sed -n 's/^MEMORY=//p' "${RUN}/status/row_${row}.status")" || return 1
    row_paths "${row}" "${ds}" "${rep}" "${initial}" "${mem}" && row_validate || return 1
  done < "${M}/selection.tsv"
}
submission_record() {
  local kind="$1" job="$2" manifest="$3" mem="$4" tmp="${M}/attempts.tmp.$$"
  [[ "${job}" =~ ^[1-9][0-9]*$ ]] || return 1
  [[ ! -e "${M}/submissions/${job}.tsv" && ! -L "${M}/submissions/${job}.tsv" ]] || return 1
  ecoda_atomic_write "${M}/submissions/${job}.tsv" "${kind}\t${job}\t${manifest}\t${mem}\n" && publish "${M}/submissions/${job}.tsv" singlecell_bulk_submission || return 1
  : > "${tmp}" || return 1
  if [[ -f "${M}/attempts.tsv" ]];then ecoda_validate_checksum "${M}/attempts.tsv" && cat "${M}/attempts.tsv" > "${tmp}" || return 1;fi
  printf '%s\t%s\t%s\t%s\n' "${kind}" "${job}" "${manifest}" "${mem}" >> "${tmp}" || return 1
  awk -F '\t' 'NF!=4||seen[$2]++{exit 1}' "${tmp}" || return 1
  mv "${tmp}" "${M}/attempts.tsv" && ecoda_write_checksum "${M}/attempts.tsv" || return 1
  ecoda_atomic_write "${M}/scheduler_ids.tsv" "$(cut -f1,2 "${M}/attempts.tsv")\n" && ecoda_write_checksum "${M}/scheduler_ids.tsv"
}
submit_array() {
  local manifest="$1" mem="$2" count raw rc
  JOB="";manifest_validate "${manifest}" || return 1;count="$(wc -l < "${manifest}" | tr -d '[:space:]')"
  set +e
  raw="$(sbatch --parsable --partition="${SLURM_PARTITION_BENCHMARK_CPU}" --array="1-${count}" --ntasks=1 --cpus-per-task=8 --mem="${mem}G" --time=02:00:00 --output="${RUN}/logs/geometry_%A_%a.log" --error="${RUN}/logs/geometry_%A_%a.err" --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" "${SELF}" --worker --manifest "${manifest}" --memory "${mem}")";rc=$?
  set -e
  ecoda_atomic_write "${manifest}.submission.txt" "${raw}\n" && publish "${manifest}.submission.txt" singlecell_bulk_submission || return 1
  JOB="${raw%%;*}";[[ "${JOB}" =~ ^[1-9][0-9]*$ ]] || return 1
  submission_record ARRAY "${JOB}" "${manifest}" "${mem}" || return 1;echo "SINGLECELL_BULK_ARRAY_ID=${JOB}"
  [[ ${rc} -eq 0 ]]
}
if [[ "${MODE}" == controller ]];then
  [[ -z "${MANIFEST}${MEMORY}" && ! -e "${RUN}" && ! -L "${RUN}" && "${CONTROL}" == "${SCRATCH}/_ecoda_inputs/"* && "${CONTROL_SHA}" =~ ^[a-f0-9]{64}$ && -f "${CONTROL}" && ! -L "${CONTROL}" && "$(ecoda_realpath_existing "${CONTROL}")" == "${CONTROL}" && "$(ecoda_sha256_file "${CONTROL}")" == "${CONTROL_SHA}" ]] || exit 1
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  trap 'ecoda_set_run_state FAIL singlecell_bulk_controller_failed || true' EXIT
  ecoda_atomic_write "${RUN}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest";ecoda_runtime_validate_submission apptainer
  ecoda_atomic_copy "${CONTROL}" "${M}/control.json";ecoda_atomic_write "${M}/control.sha256" "${CONTROL_SHA}\n"
  ecoda_atomic_write "${M}/operation" "singlecell-batch-main-bulk\n";expected_selection > "${M}/selection.tsv"
  awk -F '\t' '$4==64' "${M}/selection.tsv" > "${M}/initial_64.tsv";awk -F '\t' '$4==128' "${M}/selection.tsv" > "${M}/initial_128.tsv"
  "${PYTHON_BIN}" -B - "${M}" <<'PY'
import hashlib,json,pathlib,sys
m=pathlib.Path(sys.argv[1]);c=json.loads((m/'control.json').read_text())
assert c['schema']=='ecoda_singlecell_batch_main_v2'
with (m/'inputs.tsv').open('x')as f:
 for p in [m/'control.json',pathlib.Path(c['stage3_receipt']['path']),pathlib.Path(c['stage3_review']['path'])]:f.write(str(p)+'\t'+hashlib.md5(p.read_bytes()).hexdigest()+'\n')
PY
  for NAME in "${SCOPE[@]}";do publish "${M}/${NAME}" singlecell_bulk_scope;done
else
  [[ -z "${CONTROL}${CONTROL_SHA}" ]] || exit 1;ecoda_open_run "${ECODA_RUN_ID}" stage5
fi
scope_validate
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)"
if [[ "${MODE}" == worker ]];then
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker --manifest "${MANIFEST}" --memory "${MEMORY}"
  manifest_validate "${MANIFEST}";[[ "${SLURM_ARRAY_TASK_ID:-}" =~ ^[1-9][0-9]*$ ]] || exit 1
  LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")";[[ -n "${LINE}" ]] || exit 1
  IFS=$'\t' read -r ROW DS REP INITIAL <<< "${LINE}";row_paths "${ROW}" "${DS}" "${REP}" "${INITIAL}" "${MEMORY}"
  if [[ -e "${STATUS}" || -L "${STATUS}" ]];then row_validate;exit 0;fi
  [[ ! -e "${OUTPUT}" && ! -L "${OUTPUT}" ]] || exit 1
  mkdir -p "${OUTPUT%/*}";OWNED=("${OUTPUT}");ecoda_artifact_owner_acquire "${OUTPUT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
  failed_worker() { local p owner;for p in "${OWNED[@]}";do owner="$(ecoda_artifact_owner_dir "${p}")";if _ecoda_artifact_owner_validate_dir "${owner}" "${p}" && [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]];then ecoda_artifact_owner_set_state "${p}" FAIL singlecell_bulk_attempt_failed >/dev/null || true;fi;done; }
  trap failed_worker EXIT
  export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8 NUMBA_CACHE_DIR="${RUN}/cache/${ROW}_${MEM}/numba" MPLCONFIGDIR="${RUN}/cache/${ROW}_${MEM}/matplotlib"
  mkdir -p "${NUMBA_CACHE_DIR}" "${MPLCONFIGDIR}";cd "${PROJECT_ROOT}"
  python_run --dataset "${DS}" --representation "${REP}" --output "${OUTPUT}" --run-id "${ECODA_RUN_ID}"
  python_run --dataset "${DS}" --representation "${REP}" --output "${OUTPUT}" --run-id "${ECODA_RUN_ID}" --validate-output
  [[ "$(cut -f1 "${OUTPUT}/files.tsv")" == "$(printf '%s\n' "${OUTPUT}/embedding.h5ad" "${OUTPUT}/metrics.json" "${OUTPUT}/display.feather" "${OUTPUT}/phases.json")" ]]
  while IFS=$'\t' read -r P MD5 EXTRA;do
    [[ "${MD5}" =~ ^[a-f0-9]{32}$ && -z "${EXTRA}" && "$(ecoda_md5_file "${P}")" == "${MD5}" ]]
    ecoda_artifact_owner_acquire "${P}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null;OWNED+=("${P}")
  done < "${OUTPUT}/files.tsv"
  ecoda_artifact_owner_acquire "${OUTPUT}/files.tsv" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null;OWNED+=("${OUTPUT}/files.tsv")
  for P in "${OWNED[@]:1}";do publish "${P}" "${PRODUCER}";done
  scope_validate
  for P in "${OWNED[@]}";do ecoda_artifact_owner_set_state "${P}" OK singlecell_bulk_attempt_complete >/dev/null;done
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\nROW=${ROW}\nMEMORY=${MEM}\nOUTPUT=${OUTPUT}\n"
  ecoda_artifact_owner_acquire "${STATUS}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null;publish "${STATUS}" singlecell_bulk_row;ecoda_artifact_owner_set_state "${STATUS}" OK singlecell_bulk_row_complete >/dev/null
  row_validate;trap - EXIT;exit 0
fi
if [[ "${MODE}" == watchdog ]];then
  [[ -z "${MANIFEST}${MEMORY}" && -n "${SLURM_JOB_ID:-}" ]] || exit 1
  # Controller is the sole initial ledger writer; watchdog takes over after this marker.
  for N in {1..60};do [[ ! -f "${M}/submission.ready" ]] || break;sleep 1;done
  ecoda_validate_artifact_record "${M}/submission.ready" singlecell_bulk_submission "${ECODA_RUN_ID}" >/dev/null
  ecoda_validate_checksum "${M}/attempts.tsv";ecoda_validate_checksum "${M}/scheduler_ids.tsv"
  HAD_FAILURE=0
  watchdog_fail() { ecoda_atomic_write "${RUN}/status/singlecell_bulk_watchdog.status" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\n" || true; }
  trap watchdog_fail EXIT
  STARTS=(64 128);MEMORIES=(64 128);CEILINGS=(256 512)
  CURRENTS=("${M}/initial_64.tsv" "${M}/initial_128.tsv");JOBS=()
  for CURRENT in "${CURRENTS[@]}";do JOBS+=("$(awk -F '\t' -v p="${CURRENT}" '$1=="ARRAY"&&$3==p{print $2}' "${M}/attempts.tsv")");done
  # Each round submits both independent OOM groups before waiting on retries.
  for ROUND in 0 1 2;do
    for INDEX in 0 1;do
      START="${STARTS[${INDEX}]}";CURRENT="${CURRENTS[${INDEX}]}";MEM="${MEMORIES[${INDEX}]}";JOB="${JOBS[${INDEX}]}"
      if [[ -z "${JOB}" ]];then continue;fi
      if [[ ! "${JOB}" =~ ^[1-9][0-9]*$ ]];then HAD_FAILURE=1;JOBS[${INDEX}]="";continue;fi
      manifest_validate "${CURRENT}";COUNT="$(wc -l < "${CURRENT}" | tr -d '[:space:]')"
      if ! ecoda_wait_array_accounting "${JOB}" "${COUNT}" "${ECODA_SINGLECELL_POLL_SECONDS:-30}";then HAD_FAILURE=1;JOBS[${INDEX}]="";continue;fi
      manifest_validate "${CURRENT}"
      RAW="${M}/accounting_${JOB}.txt";ecoda_atomic_write "${RAW}" "${ECODA_ACCOUNTING_ROWS}\n";publish "${RAW}" singlecell_bulk_accounting
      OOM=();SEEN=" "
      while IFS='|' read -r J STATE EXITCODE;do
        [[ "${J}" =~ ^${JOB}_([0-9]+)$ ]] || continue;TASK="${BASH_REMATCH[1]}"
        [[ ${TASK} -ge 1 && ${TASK} -le ${COUNT} && "${SEEN}" != *" ${TASK} "* ]] || exit 1;SEEN="${SEEN}${TASK} "
        case "${STATE%%+*}" in
          COMPLETED) [[ "${EXITCODE}" == 0:0 ]] || HAD_FAILURE=1 ;;
          OUT_OF_MEMORY) OOM+=("${TASK}") ;;
          *) HAD_FAILURE=1 ;;
        esac
      done <<< "${ECODA_ACCOUNTING_ROWS}"
      [[ "$(wc -w <<< "${SEEN}" | tr -d '[:space:]')" == "${COUNT}" ]] || exit 1
      JOBS[${INDEX}]=""
      [[ ${#OOM[@]} -gt 0 ]] || continue
      if [[ ${MEM} -ge ${CEILINGS[${INDEX}]} ]];then HAD_FAILURE=1;continue;fi
      NEXT=$((MEM * 2));ATTEMPT=$((ROUND + 1));RETRY="${M}/retry_${START}_${ATTEMPT}.tsv";TMP="${RETRY}.build.$$";: > "${TMP}"
      for TASK in "${OOM[@]}";do sed -n "${TASK}p" "${CURRENT}" >> "${TMP}";done
      RETRY_MD5="$(ecoda_md5_file "${TMP}")"
      [[ ! -e "${RETRY}" && ! -L "${RETRY}" ]];ecoda_atomic_install_manifest "${TMP}" "${RETRY}" 4;rm -f "${TMP}";publish "${RETRY}" singlecell_bulk_scope
      manifest_validate "${RETRY}";[[ "$(ecoda_md5_file "${RETRY}")" == "${RETRY_MD5}" ]]
      if ! submit_array "${RETRY}" "${NEXT}";then HAD_FAILURE=1;continue;fi
      CURRENTS[${INDEX}]="${RETRY}";MEMORIES[${INDEX}]="${NEXT}";JOBS[${INDEX}]="${JOB}"
    done
  done
  [[ ! -e "${M}/initial_submission.failure" ]] || HAD_FAILURE=1
  [[ ${HAD_FAILURE} -eq 0 ]];all_rows_validate
  ecoda_atomic_write "${RUN}/status/singlecell_bulk_watchdog.status" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\n";trap - EXIT;exit 0
fi
[[ -z "${MANIFEST}${MEMORY}" ]] || exit 1
if [[ "${MODE}" == validate ]];then
  all_rows_validate
  [[ "$(cat "${RUN}/status/singlecell_bulk_watchdog.status")" == "$(printf 'STATE=OK\nRUN_ID=%s\n' "${ECODA_RUN_ID}")" && "$(sed -n 's/^STATE=//p' "${RUN}/status/terminal")" == OK ]]
  for NAME in scheduler_ids.tsv attempts.tsv;do ecoda_validate_artifact_record "${M}/${NAME}" singlecell_bulk_submission "${ECODA_RUN_ID}" >/dev/null;done
  [[ "$(cut -f1,2 "${M}/attempts.tsv")" == "$(cat "${M}/scheduler_ids.tsv")" ]]
  while IFS=$'\t' read -r KIND ID PATH_IN MEM_IN;do
    ecoda_validate_artifact_record "${M}/submissions/${ID}.tsv" singlecell_bulk_submission "${ECODA_RUN_ID}" >/dev/null
    [[ "$(cat "${M}/submissions/${ID}.tsv")" == "$(printf '%s\t%s\t%s\t%s' "${KIND}" "${ID}" "${PATH_IN}" "${MEM_IN}")" ]]
    if [[ "${KIND}" == ARRAY ]];then manifest_validate "${PATH_IN}";fi
  done < "${M}/attempts.tsv"
  awk -F '\t' -v m="${M}" 'NF!=4||$2!~/^[1-9][0-9]*$/||seen[$2]++{exit 1}
    $1=="WATCHDOG"{if($3!="-"||$4!=0)exit 1;w++;next}
    $1!="ARRAY"||paths[$3]++{exit 1}
    $3==m"/initial_64.tsv"&&$4==64{a++;next}
    $3==m"/initial_128.tsv"&&$4==128{b++;next}
    $3==m"/retry_64_1.tsv"&&$4==128{next}
    $3==m"/retry_64_2.tsv"&&$4==256{next}
    $3==m"/retry_128_1.tsv"&&$4==256{next}
    $3==m"/retry_128_2.tsv"&&$4==512{next}
    {exit 1} END{if(a!=1||b!=1||w!=1)exit 1}' "${M}/attempts.tsv"
  exit 0
fi
DEPENDENCIES=""
for MEM in 64 128;do
  if submit_array "${M}/initial_${MEM}.tsv" "${MEM}";then :;else ecoda_atomic_write "${M}/initial_submission.failure" "STATE=FAIL\n";fi
  if [[ "${JOB:-}" =~ ^[1-9][0-9]*$ ]];then DEPENDENCIES="${DEPENDENCIES:+${DEPENDENCIES}:}${JOB}";fi
done
[[ -n "${DEPENDENCIES}" ]] || exit 1
RAW="$(sbatch --parsable --partition="${SLURM_PARTITION_BENCHMARK_CPU}" --dependency="afterany:${DEPENDENCIES}" --ntasks=1 --cpus-per-task=1 --mem=2G --time="${WATCHDOG_TIME_LIMIT}" --output="${RUN}/logs/watchdog_%j.log" --error="${RUN}/logs/watchdog_%j.err" --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" "${SELF}" --watchdog)"
ecoda_atomic_write "${M}/watchdog_submission.txt" "${RAW}\n";WD="${RAW%%;*}";[[ "${WD}" =~ ^[1-9][0-9]*$ ]]
submission_record WATCHDOG "${WD}" - 0;echo "SINGLECELL_BULK_WATCHDOG_ID=${WD}"
ecoda_atomic_write "${M}/submission.ready" "STATE=READY\n";publish "${M}/submission.ready" singlecell_bulk_submission
WAIT_RC=0;ecoda_wait_scalar_accounting "${WD}" "${ECODA_SINGLECELL_POLL_SECONDS:-30}" || WAIT_RC=$?
[[ ${WAIT_RC} -eq 0 ]]
case "${ECODA_ACCOUNTING_STATE}" in COMPLETED|FAILED|CANCELLED*|OUT_OF_MEMORY|TIMEOUT|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE) : ;; *) exit 1 ;; esac
# Mutable checkpoints become PUBLISHED only after the sole watchdog writer ends.
for NAME in attempts.tsv scheduler_ids.tsv;do publish "${M}/${NAME}" singlecell_bulk_submission;done
[[ "${ECODA_ACCOUNTING_STATE}" == COMPLETED ]]
awk -F '|' -v job="${WD}" '$1==job{n++;if($2!="COMPLETED"||$3!="0:0")bad=1}END{exit(n==1&&!bad?0:1)}' <<< "${ECODA_ACCOUNTING_ROWS}"
scope_validate;all_rows_validate
ecoda_set_run_state OK singlecell_bulk_complete;trap - EXIT
# Finite NAS publication remains a separate reviewed existing helper operation.
