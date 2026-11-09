#!/bin/bash
# Finite compatibility evidence only: no fits, scores, approval or NAS sync.
set -euo pipefail
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && "${ECODA_SOURCE_ROOT:-}" == /* ]] || exit 1
# Slurm executes its spool copy; recover the pinned source before configuration.
SCRIPT_DIR="${ECODA_SOURCE_ROOT%/}/src/utils/bash"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
SELF="${SCRIPT_DIR}/legacy_analysis_audit.sh"
PURPOSE=legacy_analysis_compatibility_audit
PRODUCER=legacy_analysis_compatibility
MODE=controller
if [[ "${1:-}" == --worker || "${1:-}" == --validate ]]; then MODE="${1#--}"; shift; fi
CONTROL=""; CONTROL_SHA=""
while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || exit 1
  case "$1" in
    --control) [[ -z "${CONTROL}" ]] || exit 1; CONTROL="$2" ;;
    --control-sha256) [[ -z "${CONTROL_SHA}" ]] || exit 1; CONTROL_SHA="$2" ;;
    *) echo 'Usage: legacy_analysis_audit.sh --control ABS_JSON --control-sha256 SHA256' >&2; exit 1 ;;
  esac
  shift 2
done
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && -n "${ECODA_RUN_ID:-}" ]] || exit 1
ecoda_validate_run_id "${ECODA_RUN_ID}"
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
ecoda_require_source_script_path "${PROJECT_ROOT}/src/notebooks/audit_legacy_analysis_inputs.R" "${ECODA_SOURCE_ROOT}" >/dev/null
SCRATCH="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
RUN="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${RUN}" ]] || exit 1
if [[ "${MODE}" == controller ]]; then
  [[ ! -e "${RUN}" && ! -L "${RUN}" && "${CONTROL}" == "${SCRATCH}/_ecoda_inputs/"* &&
     "${CONTROL_SHA}" =~ ^[a-f0-9]{64}$ ]] || exit 1
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  ecoda_atomic_write "${RUN}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
else
  [[ -z "${CONTROL}${CONTROL_SHA}" ]] || exit 1
  ecoda_open_run "${ECODA_RUN_ID}" stage5
fi
export ECODA_RUN_ROOT="${RUN}" ECODA_RUNTIME_PROFILE=stage5 RETICULATE_PYTHON="${PYTHON_BIN}"
M="${RUN}/manifests"
REPORT="${RUN}/artifacts/evidence/legacy_analysis_compatibility.json"
STATUS="${RUN}/status/legacy_analysis_audit.status"
SCOPE_FILES=(request.json control.json rows.209.tsv files.106.tsv author_hr_metadata_pairs.9.tsv
  avg_legacy.10.tsv zeroimp_nas_counterparts.9.tsv observation.json zeroimp_schema.json zeroimp_policy.json
  feather_headers.json input_review.json input_witness.json parent_gate.json parent_inspect.json
  parent_accounting.tsv source.manifest runtime.identity operation selection.tsv inputs.tsv scope.json)
failed() {
  local owner
  ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\n" || true
  if [[ "${MODE}" == controller ]]; then ecoda_set_run_state FAIL legacy_compatibility_audit_failed || true; fi
  owner="$(ecoda_artifact_owner_dir "${REPORT}")" || return 0
  if _ecoda_artifact_owner_validate_dir "${owner}" "${REPORT}" &&
     [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]; then
    ecoda_artifact_owner_set_state "${REPORT}" FAIL legacy_compatibility_audit_failed >/dev/null || true
  fi
}
if [[ "${MODE}" != validate ]]; then trap failed EXIT; fi

# Metadata and integer stat checks only; the R worker alone opens result RDSs
# and configured H5AD obs. No current H5AD content hash is calculated here.
metadata() {
  "${PYTHON_BIN}" - "$1" "${RUN}" "${PROJECT_ROOT}" "${CONTROL}" "${CONTROL_SHA}" <<'PY'
import csv,hashlib,io,json,os,pathlib,re,sys,tempfile
mode,root,project,request,request_sha=sys.argv[1:];root=pathlib.Path(root);m=root/'manifests';project=pathlib.Path(project)
def need(ok,message):
    if not ok: raise SystemExit('ERROR: '+message)
def sha(data): return hashlib.sha256(data).hexdigest()
def regular(name):
    p=pathlib.Path(name);need(not any(c in str(p) for c in '\t\r\n') and p.is_absolute() and p.is_file() and not p.is_symlink() and p.resolve()==p,'noncanonical file: '+str(p));return p
def small(name):
    p=regular(name);need(p.stat().st_size<=16*1024*1024,'oversized metadata: '+str(p));return p.read_bytes()
def ref(value):
    need(isinstance(value,dict) and set(value)=={'path','sha256'} and re.fullmatch('[a-f0-9]{64}',value['sha256']),'invalid explicit SHA reference')
    data=small(value['path']);need(sha(data)==value['sha256'],'metadata SHA mismatch: '+value['path']);return data
def put(name,data):
    fd,tmp=tempfile.mkstemp(prefix='.legacy-scope-',dir=m)
    try:
        with os.fdopen(fd,'wb') as f: f.write(data);f.flush();os.fsync(f.fileno())
        os.link(tmp,m/name)
    finally: os.unlink(tmp)
def js(data): return json.loads(data)
def table(name): return list(csv.DictReader(io.StringIO(small(m/name).decode()),delimiter='\t'))
def stat(name):
    s=regular(name).stat();return dict(dev=s.st_dev,ino=s.st_ino,size=s.st_size,mode=s.st_mode&0o7777,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)
def fields(data):
    lines=data.decode().splitlines();pairs=[x.split('=',1) for x in lines]
    need(all(len(x)==2 for x in pairs) and len({x[0] for x in pairs})==len(pairs),'duplicate/malformed evidence fields');return dict(pairs)
entries=('science_control','observation','zeroimp_schema','zeroimp_policy','feather_headers','input_review')
if mode=='bind':
    raw=small(request);need(re.fullmatch('[a-f0-9]{64}',request_sha) and sha(raw)==request_sha,'request SHA mismatch');q=js(raw)
    need(set(q)=={'schema',*entries} and q['schema']=='ecoda_legacy_analysis_audit_gate_v1','unknown gate request fields')
    put('request.json',raw)
    for key in entries: put('control.json' if key=='science_control' else key+'.json',ref(q[key]))
    control=js(small(m/'control.json'))
    expected={'rows.209.tsv','files.106.tsv','author_hr_metadata_pairs.9.tsv','avg_legacy.10.tsv','zeroimp_nas_counterparts.9.tsv'}
    need(set(control['controls'])==expected,'unexpected control file union')
    for name in expected: put(name,ref(dict(path=str(pathlib.Path(q['science_control']['path']).parent/name),sha256=control['controls'][name])))
    review=js(small(m/'input_review.json'))
    for key,name in [('observation','input_witness.json'),('gate_manifest','parent_gate.json'),('terminal_inspect','parent_inspect.json'),('accounting','parent_accounting.tsv')]: put(name,ref(review[key]))
else: q=js(small(m/'request.json'))
control=js(small(m/'control.json'));rows=table('rows.209.tsv');files=table('files.106.tsv');pairs=table('author_hr_metadata_pairs.9.tsv')
need(control['schema']=='ecoda_p6_legacy_science_audit_proposal_v4' and control['total_logical_rows']==209 and control['physical_files']==106 and control['author_hr_metadata_pairs']==9 and control['scientifically_approved_rows']==0,'wrong finite science control')
need((len(rows),len(files),len(pairs))==(209,106,9),'finite manifest counts differ')
keys=[(r['dataset'],r['view'],r['method_key']) for r in rows]
need(len(set(keys))==209 and all(r['view']=='benchmark_analysis' and r['origin']=='legacy_reuse' and r['status']=='BLOCKED_UNTIL_SCIENTIFIC_REVIEW' for r in rows),'logical scope/origin differs')
need(all(all(v and not any(c in v for c in '\t\r\n') for v in k) for k in keys),'invalid selection token')
config=js(small(project/'datasets.json'))
need(all(config.get(r['dataset'],{}).get('use_for_benchmark') is True for r in rows),'disabled/unconfigured cohort')
inputs={}
for r in rows:
    v={k:r[k] for k in ('current_input_path','current_input_md5','current_input_size','current_input_proof')}
    need(r['dataset'] not in inputs or inputs[r['dataset']]==v,'conflicting input references');inputs[r['dataset']]=v
need(len(inputs)==10 and len({r['artifact_path'] for r in files})==106,'input/physical union differs')
review=js(small(m/'input_review.json'));w=js(small(m/'input_witness.json'));gate=js(small(m/'parent_gate.json'));inspect=js(small(m/'parent_inspect.json'))
need(review['schema']=='ecoda_p6_preflight_witness_review_v1' and review['decision']=='APPROVED_INPUT_CONTENT_WITNESS_ONLY' and review['scientifically_approved_rows']==0 and review['formal_pipeline_dependency'] is False and review['review_complete'] is True,'pending/unapproved input witness')
for key,name in [('observation','input_witness.json'),('gate_manifest','parent_gate.json'),('terminal_inspect','parent_inspect.json'),('accounting','parent_accounting.tsv')]: need(sha(small(m/name))==review[key]['sha256'],'review evidence differs')
need(w['schema']=='ecoda_p6_preflight_witness_observation_v1' and w['scientifically_approved'] is False and w['h5ad_payload_bytes_read']==0 and len(w['rows'])==10,'invalid original input observation')
need(gate['gate_id']==inspect['gate_id']==w['parent_gate_id'] and gate['state']==inspect['state']==review['parent_terminal_state'] and gate['state'] in ('FAILED','COMPLETED'),'parent not terminal or review differs')
a=inspect['audit']['accounting'];raw=small(m/'parent_accounting.tsv')
need(gate['schema_version']=='1.0' and inspect['operation']=='inspect' and gate['audit']['accounting']==a,'gate and original inspection schemas/evidence differ')
need(a['attempted'] is True and a['query_count']==1 and a['actual_exit']==0 and raw==a['stdout'].encode() and b'[output truncated]' not in raw,'original complete terminal accounting required')
requested=gate['scheduler_array_ids']+gate['scheduler_watchdog_ids']
need(requested and len(set(requested))==len(requested) and all(re.fullmatch('[0-9]+',v) for v in requested) and a['scheduler_ids']==requested and a['command']=='sacct -n -P -X -j '+','.join(requested)+' --format=JobIDRaw,State,ExitCode','original accounting query scope differs')
account={};parsed=[]
need(raw.endswith(b'\n'),'unterminated original accounting')
for number,line in enumerate(raw.decode().splitlines(),1):
    p=line.split('|');need(len(p)==3 and re.fullmatch('[0-9]+',p[0]) and p[0] not in account,'malformed/duplicate original accounting');account[p[0]]=p[1:]
    parsed.append(dict(job_id=p[0],state=p[1],exit_code=p[2],line=number))
need(a['format']=='JobIDRaw|State|ExitCode' and a['rows']==parsed and a['malformed_rows']==[],'inspect parsed accounting differs from preserved raw rows')
selected={(r['dataset'],r['fullhash_array_id'],r['task_id'],r['current_input_path'],r['current_input_md5'],str(r['current_input_size'])) for r in w['rows']}
approved={(r['dataset'],r['fullhash_array_id'],r['task_id'],r['current_input_path'],r['current_input_md5'],str(r['current_input_size'])) for r in review['selected_preflight_tasks'] if r['status']=='PASS'}
need(len(selected)==10 and selected==approved and len(review['selected_preflight_tasks'])==10,'selected preflight review scope differs')
metadata_paths=[]
for key in ('source','runtime','input_manifest'):
    v=w[key];data=small(v['path']);need(sha(data)==v['sha256'] and data.decode()==v['text'],'prior source/runtime/input manifest changed');metadata_paths.append(v['path'])
parent=pathlib.Path(w['source']['path']).parent.parent
need(parent.name==w['parent_run_id'] and str(parent.parent)==str(root.parent),'foreign preflight run root')
need(w['source']['path']==str(parent/'manifests/source.manifest') and w['runtime']['path']==str(parent/'manifests/runtime.identity') and w['input_manifest']['path']==str(parent/'manifests/h5ad_preflight.tsv'),'noncanonical prior identity paths')
preflight=[line.split('\t') for line in w['input_manifest']['text'].splitlines()]
need(len(preflight)==11 and all(len(row)==3 and row[1]=='benchmark_analysis' and row[2].startswith('/') for row in preflight) and len({row[0] for row in preflight})==11 and {row[0] for row in preflight}==set(inputs)|{'Stephenson'},'prior full-hash physical manifest differs')
for r in w['rows']:
    need(r['dataset'] in inputs and all(str(r[k])==v for k,v in inputs[r['dataset']].items()),'witness/current-input union differs')
    need(type(r['task_id']) is int and 1<=r['task_id']<=len(preflight) and preflight[r['task_id']-1]==[r['dataset'],'benchmark_analysis',r['current_input_path']],'preflight task does not match its exact physical manifest row')
    need(r['view']=='benchmark_analysis' and r['no_change_timestamps_since_record'] is True and account.get(r['fullhash_array_id'])==['COMPLETED','0:0'],'selected full-hash array did not succeed')
    need(r['fullhash_array_id'] in gate['scheduler_array_ids'],'unrecorded full-hash array')
    current=stat(r['current_input_path']);need(current==r['current_stat'] and current['size']==int(r['current_input_size']),'current input changed since observation')
    for key in ('proof','record','status','log'):
        v=r[key];data=small(v['path']);need(sha(data)==v['sha256'] and data.decode()==v['text'] and stat(v['path'])==v['stat'],'preflight evidence changed: '+key);metadata_paths.append(v['path'])
    expected_record=parent/'manifests/artifacts'/(sha(r['current_input_path'].encode())[:32]+'.record')
    need(r['record']['path']==str(expected_record),'noncanonical preflight record')
    need(fields(small(expected_record))==dict(PATH=r['current_input_path'],SIZE=str(r['current_input_size']),MD5=r['current_input_md5'],RUN_ID=w['parent_run_id'],PRODUCER='stage5_preflight',STATE='PUBLISHED'),'wrong full-hash record')
    status=fields(small(r['status']['path']))
    need(r['status']['path']==str(parent/'status/h5ad_preflight'/(r['dataset']+'__benchmark_analysis.status')) and status==dict(STATE='OK',RUN_ID=w['parent_run_id'],DATASET=r['dataset'],VIEW='benchmark_analysis',TASK_ID=str(r['task_id'])),'wrong preflight task status')
    need(r['log']['path']==str(parent/'logs'/('h5ad_preflight_stage5_'+r['fullhash_array_id']+'_'+str(r['task_id'])+'.log')),'wrong preflight log identity')
    need(fields(small(r['current_input_proof']))==dict(MD5=r['current_input_md5'],SIZE=str(r['current_input_size']),PATH=r['current_input_path']),'wrong expected input proof')
    need(current['mtime_ns']<=r['record']['stat']['mtime_ns'] and current['ctime_ns']<=r['record']['stat']['mtime_ns'],'input timestamps postdate full-hash record')
paths=list(dict.fromkeys([r['artifact_path'] for r in files]+[v['current_input_path'] for v in inputs.values()]+metadata_paths))
selection=(''.join('\t'.join(k)+'\n' for k in keys)).encode()
if mode=='bind':
    put('operation',b'legacy-audit\n');put('selection.tsv',selection);put('inputs.tsv',(''.join(p+'\n' for p in paths)).encode())
    names=[p.name for p in m.iterdir() if p.is_file()]
    identity=dict(schema=1,run_id=root.name,purpose='legacy_analysis_compatibility_audit',config_sha256=sha(small(project/'datasets.json')),metadata_sha256={n:sha(small(m/n)) for n in names},input_stats={p:stat(p) for p in paths})
    put('scope.json',(json.dumps(identity,sort_keys=True)+'\n').encode())
else:
    identity=js(small(m/'scope.json'));need(identity['run_id']==root.name and identity['purpose']=='legacy_analysis_compatibility_audit' and identity['config_sha256']==sha(small(project/'datasets.json')),'scope/source config changed')
    need(small(m/'selection.tsv')==selection and small(m/'inputs.tsv')==(''.join(p+'\n' for p in paths)).encode(),'selection/input inventory changed')
    for n,h in identity['metadata_sha256'].items(): need(sha(small(m/n))==h,'frozen metadata changed: '+n)
    need(identity['input_stats']=={p:stat(p) for p in paths},'selected input/evidence stat changed during audit')
if mode=='report':
    x=js(small(root/'artifacts/evidence/legacy_analysis_compatibility.json'))
    need(x['schema']=='ecoda_legacy_analysis_compatibility_v1' and x['physical_files']==106 and x['logical_rows']==209 and len(x['rows'])==209 and len(x['author_hr_pairs'])==9 and len(x['cohorts'])==10,'report scope differs')
    need(x['scientifically_approved_rows']==0 and x['model_fits']==0 and x['score_calculations']==0 and x['h5ad_expression_values_read'] is False and x['h5ad_content_rehashed_here'] is False and x['prior_validated_content_witness_required'] is True,'report exceeds evidence scope')
    need(set(x['rows'])=={r['dataset']+'/'+r['method_key'] for r in rows},'report omitted selected logical rows')
    v=x['validation_identities'];pairs=(('control','control.json'),('source-identity','source.manifest'),('runtime-identity','runtime.identity'),('config',str(project/'datasets.json')),('observation','observation.json'),('zeroimp-schema','zeroimp_schema.json'),('zeroimp-policy','zeroimp_policy.json'),('feather-headers','feather_headers.json'))
    for key,name in pairs:
        p=pathlib.Path(name) if name.startswith('/') else m/name;digest='-md5' if key in ('source-identity','runtime-identity','config') else '-sha256'
        need(v['--'+key]==str(p) and v['--'+key+digest]==(hashlib.md5(small(p)).hexdigest() if digest=='-md5' else sha(small(p))),'worker identity differs: '+key)
    need(x['compatible'] is True and x['status']=='COMPATIBILITY_EVIDENCE_COMPLETE_PENDING_REVIEW','compatibility mismatches require separate review')
PY
}
if [[ "${MODE}" == controller ]]; then
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest"
  ecoda_runtime_validate_submission apptainer
  metadata bind
  for NAME in "${SCOPE_FILES[@]}"; do
    ecoda_write_checksum "${M}/${NAME}"
    ecoda_write_artifact_record "${M}/${NAME}" legacy_analysis_audit_scope "${ECODA_RUN_ID}" >/dev/null
  done
fi
validate_scope() {
  local name path
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" && "$(cat "${M}/operation")" == legacy-audit ]] || return 1
  ecoda_runtime_validate_bound_run || return 1
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest" || return 1
  for name in "${SCOPE_FILES[@]}"; do ecoda_validate_artifact_record "${M}/${name}" legacy_analysis_audit_scope "${ECODA_RUN_ID}" >/dev/null || return 1; done
  while IFS= read -r path; do ecoda_require_input_ownership "${path}" "${ECODA_RUN_ID}" || return 1; done < "${M}/inputs.tsv"
  metadata check
}
validate_report() {
  local owner
  validate_scope || return 1
  [[ "$(sed -n 's/^STATE=//p' "${STATUS}")" == OK && "$(sed -n 's/^RUN_ID=//p' "${STATUS}")" == "${ECODA_RUN_ID}" ]] || return 1
  ecoda_validate_artifact_record "${REPORT}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
  owner="$(ecoda_artifact_owner_dir "${REPORT}")" || return 1
  _ecoda_artifact_owner_validate_dir "${owner}" "${REPORT}" || return 1
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  metadata report
}
if [[ "${MODE}" == validate ]]; then
  validate_report || exit 1
  [[ "$(sed -n 's/^STATE=//p' "${RUN}/status/terminal")" == OK && "$(sed -n 's/^RUN_ID=//p' "${RUN}/status/terminal")" == "${ECODA_RUN_ID}" ]] || exit 1
  for NAME in scheduler_ids.tsv sbatch_submission.txt; do ecoda_validate_artifact_record "${M}/${NAME}" legacy_analysis_audit_submission "${ECODA_RUN_ID}" >/dev/null; done
  echo "LEGACY_ANALYSIS_AUDIT_VALIDATED=${ECODA_RUN_ID}"
  exit 0
fi
if [[ "${MODE}" == worker ]]; then
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker
  cd "${PROJECT_ROOT}"
  validate_scope || exit 1
  [[ ! -e "${REPORT}" && ! -L "${REPORT}" ]] || exit 1
  args=(--control "${M}/control.json" --control-sha256 "$(ecoda_sha256_file "${M}/control.json")"
    --config "${PROJECT_ROOT}/datasets.json" --config-md5 "$(ecoda_md5_file "${PROJECT_ROOT}/datasets.json")"
    --source-identity "${M}/source.manifest" --source-identity-md5 "$(ecoda_md5_file "${M}/source.manifest")"
    --runtime-identity "${M}/runtime.identity" --runtime-identity-md5 "$(ecoda_md5_file "${M}/runtime.identity")")
  for KEY in observation zeroimp-schema zeroimp-policy feather-headers; do
    NAME="${KEY//-/_}.json"; args+=("--${KEY}" "${M}/${NAME}" "--${KEY}-sha256" "$(ecoda_sha256_file "${M}/${NAME}")")
  done
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
  set +e
  ${PIXI_RSCRIPT} --vanilla "${PROJECT_ROOT}/src/notebooks/audit_legacy_analysis_inputs.R" "${args[@]}" --report "${REPORT}"
  RC=$?
  set -e
  validate_scope || exit 1
  if [[ -f "${REPORT}" && ! -L "${REPORT}" ]]; then
    ecoda_write_checksum "${REPORT}"
    ecoda_write_artifact_record "${REPORT}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
  fi
  [[ ${RC} -eq 0 ]] || exit "${RC}"
  metadata report || exit 1
  ecoda_artifact_owner_set_state "${REPORT}" OK compatibility_evidence_pending_scientific_review >/dev/null
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\n"
  validate_report || exit 1
  trap - EXIT
  exit 0
fi
validate_scope || exit 1
mkdir -p "${REPORT%/*}"
ecoda_artifact_owner_acquire "${REPORT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)"
set +e
JOB="$(sbatch --parsable --wait --partition="${SLURM_PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=64G --time=01:00:00 \
  --output="${RUN}/logs/legacy_analysis_%j.log" --error="${RUN}/logs/legacy_analysis_%j.err" \
  --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" "${SELF}" --worker)"
RC=$?
set -e
ecoda_atomic_write "${M}/sbatch_submission.txt" "${JOB}\n"
JOB="${JOB%%;*}"; [[ "${JOB}" =~ ^[0-9]+$ ]] || { echo 'ERROR: ambiguous legacy audit submission; preserve run' >&2; exit 1; }
printf 'LEGACY_ANALYSIS_AUDIT_JOB_ID=%s\n' "${JOB}"
ecoda_atomic_write "${M}/scheduler_ids.tsv" "ARRAY\t${JOB}\n"
for NAME in scheduler_ids.tsv sbatch_submission.txt; do
  ecoda_write_checksum "${M}/${NAME}"
  ecoda_write_artifact_record "${M}/${NAME}" legacy_analysis_audit_submission "${ECODA_RUN_ID}" >/dev/null
done
[[ ${RC} -eq 0 ]] || exit "${RC}"
validate_report || exit 1
ecoda_set_run_state OK compatibility_evidence_pending_scientific_review
trap - EXIT
printf 'LEGACY_ANALYSIS_COMPATIBILITY_JSON=%s\nSCIENTIFICALLY_APPROVED_ROWS=0\n' "${REPORT}"
