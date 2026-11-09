#!/bin/bash
# One dependency barrier after the recorded login reboot. No old-run outcome
# inference, accounting query, scientific work, synchronization or old writes.
set -euo pipefail
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && "${ECODA_SOURCE_ROOT:-}" == /* ]] || exit 1
# Slurm executes its spool copy; recover the pinned source before configuration.
SCRIPT_DIR="${ECODA_SOURCE_ROOT%/}/src/utils/bash"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
SELF="${SCRIPT_DIR}/orphaned_main_observer.sh"
PURPOSE=orphaned_main_dependency_observer
PRODUCER=orphaned_main_dependency_observation
MODE=controller
if [[ "${1:-}" == --worker || "${1:-}" == --validate ]]; then MODE="${1#--}"; shift; fi
CONTROL=""; CONTROL_SHA=""
while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || exit 1
  case "$1" in
    --control) [[ -z "${CONTROL}" ]] || exit 1; CONTROL="$2" ;;
    --control-sha256) [[ -z "${CONTROL_SHA}" ]] || exit 1; CONTROL_SHA="$2" ;;
    *) echo 'Usage: orphaned_main_observer.sh --control ABS_JSON --control-sha256 SHA256' >&2; exit 1 ;;
  esac
  shift 2
done
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 ]] || exit 1
ecoda_validate_run_id "${ECODA_RUN_ID:?}"
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
SCRATCH="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
RUN="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${RUN}" ]] || exit 1
if [[ "${MODE}" == controller ]]; then
  [[ ! -e "${RUN}" && ! -L "${RUN}" && "${CONTROL}" == "${SCRATCH}/_ecoda_inputs/"* && "${CONTROL_SHA}" =~ ^[a-f0-9]{64}$ ]] || exit 1
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  ecoda_atomic_write "${RUN}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
else
  [[ -z "${CONTROL}${CONTROL_SHA}" ]] || exit 1
  ecoda_open_run "${ECODA_RUN_ID}" stage5
fi
export ECODA_RUN_ROOT="${RUN}" ECODA_RUNTIME_PROFILE=stage5
M="${RUN}/manifests"
REPORT="${RUN}/artifacts/evidence/dependency_observation.json"
STATUS="${RUN}/status/orphaned_main_observer.status"
SCOPE_FILES=(control.json original_gate.json incident_observation.json scheduler_catalog.tsv original_source.manifest original_runtime.identity source.manifest runtime.identity operation selection.tsv)
failed() {
  local owner
  ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\n" || true
  [[ "${MODE}" != controller ]] || ecoda_set_run_state FAIL observer_failed_no_old_outcome_claim || true
  owner="$(ecoda_artifact_owner_dir "${REPORT}")" || return 0
  if _ecoda_artifact_owner_validate_dir "${owner}" "${REPORT}" &&
     [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]; then
    ecoda_artifact_owner_set_state "${REPORT}" FAIL observer_failed >/dev/null || true
  fi
}
[[ "${MODE}" == validate ]] || trap failed EXIT
metadata() {
  "${PYTHON_BIN}" - "$1" "${RUN}" "${CONTROL}" "${CONTROL_SHA}" <<'PY'
import csv,datetime,hashlib,json,os,pathlib,re,shlex,sys,tempfile
mode,run,control,expected_sha=sys.argv[1:];run=pathlib.Path(run);m=run/'manifests'
def need(ok,message):
    if not ok: raise SystemExit('ERROR: '+message)
def sha(data):return hashlib.sha256(data).hexdigest()
def read(name):
    p=pathlib.Path(name);need(not any(c in str(p) for c in '\t\r\n') and p.is_absolute() and p.is_file() and not p.is_symlink() and p.resolve()==p and p.stat().st_size<=2*1024*1024,'unsafe/oversized observer metadata');return p.read_bytes()
def put(path,data):
    fd,tmp=tempfile.mkstemp(prefix='.observer-',dir=path.parent)
    try:
        with os.fdopen(fd,'wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
        os.link(tmp,path)
    finally:os.unlink(tmp)
def fields(data):
    lines=[line.split('=',1) for line in data.decode().splitlines()];need(all(len(p)==2 for p in lines) and len({p[0] for p in lines})==len(lines),'invalid identity fields');return dict(lines)
refs={'original_gate':'original_gate.json','incident_observation':'incident_observation.json','scheduler_catalog':'scheduler_catalog.tsv','original_source':'original_source.manifest','original_runtime':'original_runtime.identity'}
if mode=='bind':
    raw=read(control);need(sha(raw)==expected_sha,'control SHA differs');put(m/'control.json',raw)
q=json.loads(read(m/'control.json'))
need(set(q)=={'schema','active_dependencies',*refs} and q['schema']=='ecoda_orphaned_main_observer_v1','unknown observer control')
for key,name in refs.items():
    proof=q[key];need(set(proof)=={'path','sha256'} and re.fullmatch('[a-f0-9]{64}',proof['sha256']),'invalid observer reference')
    if mode=='bind':
        data=read(proof['path']);need(sha(data)==proof['sha256'],'reference SHA differs');put(m/name,data)
    need(sha(read(m/name))==proof['sha256'],'frozen reference changed')
gate=json.loads(read(m/'original_gate.json'));incident=json.loads(read(m/'incident_observation.json'))
source=fields(read(m/'original_source.manifest'));runtime=fields(read(m/'original_runtime.identity'))
old_id='ecoda_stage5_main_confounds_20261109_v2';parent_id='stage5_main_confounds_20261109_v2';scratch=run.parent.parent
need(gate['gate_id']==old_id and gate['state']=='PRELAUNCH_STOP' and gate.get('terminal_status') is None and gate['profile']=='durable-hpc-gate-ecoda' and gate['serialization_group']=='ecoda-benchmark' and re.fullmatch('[a-f0-9]{32}',gate['event_generation']) and gate['command_digest']==sha(gate['exact_command'].encode()),'original gate is not the exact stopped orchestration')
argv=shlex.split(gate['exact_command'])
def option(flag):
    found=[argv[i+1] for i,v in enumerate(argv[:-1]) if v==flag];need(len(found)==1,'ambiguous original wrapper option');return found[0]
need(option('--run-id')==parent_id and option('--script')=='src/5_run_benchmark_methods/stage5_campaign.sh' and option('--scratch-root')==str(scratch) and option('--source-root')==source['SOURCE_ROOT'] and source['SOURCE_COMMIT']=='9ca2453eea7320a721ef6dfd458f9db108abecac','wrong original source/wrapper')
need(q['original_source']['path']==str(scratch/'_ecoda_runs'/parent_id/'manifests/source.manifest') and q['original_runtime']['path']==str(scratch/'_ecoda_runs'/parent_id/'manifests/runtime.identity') and option('--runtime-image')==runtime['RUNTIME_IMAGE'] and option('--runtime-manifest')==runtime['RUNTIME_MANIFEST'],'wrong original identity namespace')
current_runtime=fields(read(m/'runtime.identity'))
need(all(runtime[k]==current_runtime[k] for k in ('RUNTIME_IMAGE_SHA256','RUNTIME_MANIFEST_SHA256')),'observer runtime differs from original pinned runtime')
need(incident['schema']=='ecoda_login_reboot_operational_observation_v1' and incident['hostname']=='login1.bamboo' and incident['boot_utc']=='2026-11-09T14:47:27+00:00' and incident['no_accounting_query'] is True and incident['no_scientific_payload_reads'] is True and incident['no_remote_writes'] is True,'wrong reboot observation')
need(datetime.datetime.fromisoformat(gate['launch_timestamp'])<datetime.datetime.fromisoformat(incident['boot_utc'])<=datetime.datetime.fromisoformat(incident['observed_at']),'reboot observation chronology differs')
need(incident['tmux']['exit']==1 and not incident['tmux']['stdout'] and incident['processes']['exit']==0 and gate['remote_runner'] not in incident['processes']['stdout'] and parent_id not in incident['processes']['stdout'],'runner absence was not established')
files={x['path']:x for x in incident['files']}
need(all(files.get(p,{}).get('exists') is False for p in (gate['remote_status'],gate['remote_failure_status'],gate['completion_channel'][5:])),'original completion evidence was not absent')
# These are snapshots, not assertions that live worker logs remain unchanged.
children={'benchmark':'stage5_main_benchmark_20261109_v2','uncorrected':'stage5_main_uncorrected_20261109_v2','corrected':'stage5_main_corrected_20261109_v2'}
rows=list(csv.reader(read(m/'scheduler_catalog.tsv').decode().splitlines(),delimiter='\t'))
need(len(rows)==73 and all(len(r)==4 and children.get(r[0])==r[1] and r[2] in ('PREFLIGHT','ARRAY','WATCHDOG','AGGREGATE_GATE') and re.fullmatch('[1-9][0-9]*',r[3]) for r in rows) and len({r[3] for r in rows})==73,'initial 73-ID catalog differs')
observed_catalog=set()
for context,child in children.items():
    snapshot=files.get(str(scratch/'_ecoda_runs'/child/'manifests/scheduler_ids.tsv'),{})
    need(snapshot.get('exists') is True and sha(snapshot['text'].encode())==snapshot['sha256'],'original scheduler snapshot differs')
    entries=list(csv.reader(snapshot['text'].splitlines(),delimiter='\t'))
    need(entries and all(len(r)==2 for r in entries),'original scheduler snapshot rows differ')
    observed_catalog.update((context,child,*r) for r in entries)
observed_catalog.update((ctx,children[ctx],'AGGREGATE_GATE',j) for ctx,j in (('benchmark','4665612'),('uncorrected','4665594'),('corrected','4665633')))
need({tuple(r) for r in rows}==observed_catalog,'catalog differs from original scheduler snapshots and three identified aggregate roots')
by_id={r[3]:r for r in rows}
need(all(by_id.get(j)==[ctx,children[ctx],'AGGREGATE_GATE',j] for ctx,j in (('benchmark','4665612'),('uncorrected','4665594'),('corrected','4665633'))),'aggregate context identity differs')
active=('4665568','4665570','4665612');descriptors=incident['active_job_descriptors']
need(set(descriptors)==set(active) and [r['job_id'] for r in q['active_dependencies']]==list(active) and all(set(r)=={'job_id','descriptor_sha256'} for r in q['active_dependencies']),'observer dependency scope differs')
for row,kind in zip(q['active_dependencies'],('ARRAY','WATCHDOG','AGGREGATE_GATE')):
    j=row['job_id'];d=descriptors[j]
    need(by_id.get(j)==['benchmark',children['benchmark'],kind,j] and d['argv']==['scontrol','show','job','-o',j] and d['exit']==0 and sha(d['stdout'].encode())==row['descriptor_sha256'] and d['stdout'].startswith('JobId='+j+' ') and 'JobState=PENDING' in d['stdout'],'active dependency descriptor differs')
need('Dependency=afterany:4665568' in descriptors['4665570']['stdout'] and 'Dependency=afterany:4665570' in descriptors['4665612']['stdout'],'remaining watchdog/aggregate graph differs')
selection=''.join(j+'\tdependency\tobserve_afterany\n' for j in active).encode()
if mode=='bind':put(m/'selection.tsv',selection);put(m/'operation',b'orphaned-main-observer\n')
else:need(read(m/'selection.tsv')==selection and read(m/'operation')==b'orphaned-main-observer\n','observer scope changed')
if mode in ('write','report'):
    expected=dict(schema='ecoda_orphaned_main_dependency_observation_v1',run_id=run.name,original_gate_id=old_id,original_gate_state='PRELAUNCH_STOP',original_command_digest=gate['command_digest'],original_event_generation=gate['event_generation'],control_sha256=sha(read(m/'control.json')),dependency_roots=list(active),dependency_barrier_observed=True,original_outcomes='UNKNOWN_UNTIL_ORIGINAL_INSPECT',original_wrapper_exit_code=None,initial_scheduler_catalog_size=73,final_emitted_catalog_established=False,scientifically_approved_rows=0,model_fits=0,accounting_queries=0,old_run_writes=0)
    path=run/'artifacts/evidence/dependency_observation.json'
    if mode=='write':
        job=os.environ.get('SLURM_JOB_ID','');need(re.fullmatch('[1-9][0-9]*',job) is not None,'observer must run in its own Slurm allocation');expected['observer_job_id']=job;put(path,(json.dumps(expected,sort_keys=True)+'\n').encode())
    else:
        report=json.loads(read(path));jobs=read(m/'scheduler_ids.tsv').decode().splitlines();need(len(jobs)==1,'observer must have one scheduler ID');job=jobs[0].split('\t');need(len(job)==2 and job[0]=='ARRAY' and re.fullmatch('[1-9][0-9]*',job[1]),'malformed observer scheduler ID');expected['observer_job_id']=job[1];need(report==expected,'observer receipt differs')
PY
}
if [[ "${MODE}" == controller ]]; then
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest"
  ecoda_runtime_validate_submission apptainer
  metadata bind
  for NAME in "${SCOPE_FILES[@]}"; do ecoda_write_checksum "${M}/${NAME}"; ecoda_write_artifact_record "${M}/${NAME}" orphaned_main_observer_scope "${ECODA_RUN_ID}" >/dev/null; done
fi
validate_scope() {
  local name
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" ]] || return 1
  ecoda_runtime_validate_bound_run || return 1
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest" || return 1
  for name in "${SCOPE_FILES[@]}"; do ecoda_validate_artifact_record "${M}/${name}" orphaned_main_observer_scope "${ECODA_RUN_ID}" >/dev/null || return 1; done
  metadata check
}
validate_report() {
  local owner
  validate_scope || return 1
  [[ "$(sed -n 's/^STATE=//p' "${STATUS}")" == OK && "$(sed -n 's/^RUN_ID=//p' "${STATUS}")" == "${ECODA_RUN_ID}" ]] || return 1
  ecoda_validate_artifact_record "${REPORT}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
  owner="$(ecoda_artifact_owner_dir "${REPORT}")" || return 1
  _ecoda_artifact_owner_validate_dir "${owner}" "${REPORT}" || return 1
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == OK && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 ]] || return 1
  metadata report
}
if [[ "${MODE}" == validate ]]; then
  validate_report || exit 1
  [[ "$(sed -n 's/^STATE=//p' "${RUN}/status/terminal")" == OK && "$(sed -n 's/^RUN_ID=//p' "${RUN}/status/terminal")" == "${ECODA_RUN_ID}" ]] || exit 1
  for NAME in scheduler_ids.tsv sbatch_submission.txt; do ecoda_validate_artifact_record "${M}/${NAME}" orphaned_main_observer_submission "${ECODA_RUN_ID}" >/dev/null; done
  echo "ORPHANED_MAIN_OBSERVER_VALIDATED=${ECODA_RUN_ID}"; exit 0
fi
if [[ "${MODE}" == worker ]]; then
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker
  validate_scope || exit 1
  OWNER="$(ecoda_artifact_owner_dir "${REPORT}")"
  _ecoda_artifact_owner_validate_dir "${OWNER}" "${REPORT}"
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE && ! -e "${REPORT}" && ! -L "${REPORT}" ]] || exit 1
  metadata write
  ecoda_write_checksum "${REPORT}"
  ecoda_write_artifact_record "${REPORT}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
  ecoda_artifact_owner_set_state "${REPORT}" OK dependency_barrier_only_old_outcomes_unknown >/dev/null
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\n"
  trap - EXIT; exit 0
fi
validate_scope || exit 1
mkdir -p "${REPORT%/*}"
ecoda_artifact_owner_acquire "${REPORT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)"
set +e
JOB="$(sbatch --parsable --wait --dependency=afterany:4665568:4665570:4665612 --kill-on-invalid-dep=yes --partition="${SLURM_PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:05:00 --output="${RUN}/logs/dependency_observer_%j.log" --error="${RUN}/logs/dependency_observer_%j.err" --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" "${SELF}" --worker)"
RC=$?
set -e
ecoda_atomic_write "${M}/sbatch_submission.txt" "${JOB}\n"
JOB="${JOB%%;*}"; [[ "${JOB}" =~ ^[1-9][0-9]*$ ]] || exit 1
printf 'ORPHANED_MAIN_OBSERVER_JOB_ID=%s\n' "${JOB}"
ecoda_atomic_write "${M}/scheduler_ids.tsv" "ARRAY\t${JOB}\n"
for NAME in scheduler_ids.tsv sbatch_submission.txt; do ecoda_write_checksum "${M}/${NAME}"; ecoda_write_artifact_record "${M}/${NAME}" orphaned_main_observer_submission "${ECODA_RUN_ID}" >/dev/null; done
[[ ${RC} -eq 0 ]] || exit "${RC}"
validate_report || exit 1
ecoda_set_run_state OK observer_completed_old_outcomes_unknown
trap - EXIT
printf 'DEPENDENCY_OBSERVATION_JSON=%s\nORIGINAL_OUTCOMES=UNKNOWN_UNTIL_ORIGINAL_INSPECT\n' "${REPORT}"
