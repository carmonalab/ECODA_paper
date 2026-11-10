#!/bin/bash
# Real controller/worker/common records; runtime, scheduler and scientific R
# execution are isolated here. The accepted R fixture tests the science reader.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# Execute the actual production stat helper/check without filesystem mutation.
"${ROOT}/.pixi/envs/default/bin/python" - "${ROOT}/src/utils/bash/legacy_analysis_audit.sh" <<'PY_STAT'
import ast,json,pathlib,sys
source=pathlib.Path(sys.argv[1]).read_text().split("<<'PY'\n",1)[1].split('\nPY\n}',1)[0];tree=ast.parse(source)
ns=dict(json=json,selected_beegfs=lambda path:None,stopped=True,review={'reviewed_device_transition':{'old_dev':66,'new_dev':68,'reason':'OBSERVED_LOGIN_REBOOT_REMOUNT','paths':['/selected/input']}})
for node in tree.body:
 if isinstance(node,ast.FunctionDef) and node.name in ('need','current_witness_stat','portable_stat'):
  exec(compile(ast.Module(body=[node],type_ignores=[]),'<actual-stat-helper>','exec'),ns)
# Actual e8cf historical STOP snapshot has no terminal_status field.
condition=next(n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='need' and len(n.args)>1 and isinstance(n.args[1],ast.Constant) and n.args[1].value=='historical stopped parent identity differs')
check=ast.Module(body=[ast.Expr(value=condition)],type_ignores=[]);ast.fix_missing_locations(check)
ns['w']={'parent_gate_id':'g'};base=dict(schema_version='1.0',gate_id='g',state='PRELAUNCH_STOP',remote_exit_code=None,audit_completed=False,audit_state='NOT_STARTED')
for gate in [base,dict(base,terminal_status=None)]:
 ns['gate']=gate;exec(compile(check,'<actual-stopped-parent-condition>','exec'),ns)
for gate in [dict(base,terminal_status={'state':'FAILED'}),dict(base,terminal_status=''),dict(base,remote_exit_code=1),dict(base,state='COMPLETED'),dict(base,audit_completed=True)]:
 ns['gate']=gate
 try:exec(compile(check,'<actual-stopped-parent-condition>','exec'),ns)
 except SystemExit:pass
 else:raise AssertionError('non-STOP/null parent accepted')
prior=dict(dev=66,ino=100,size=25,mode=0o600,mtime_ns=10,ctime_ns=11)
expected=dict(prior,dev=68);assert ns['current_witness_stat'](prior,'/selected/input')==expected
for bad,path in [(dict(prior,dev=67),'/selected/input'),(prior,'/foreign'),(dict(prior,ino=True),'/selected/input')]:
 try:ns['current_witness_stat'](bad,path)
 except SystemExit:pass
 else:raise AssertionError('invalid stat accepted')
# Exercise the exact production input-stat comparison for each mismatch.
loop=next(n for n in tree.body if isinstance(n,ast.For) and ast.unparse(n.target)=='r' and any(isinstance(v,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='current' for t in v.targets) for v in n.body))
check=next(v for v in loop.body if isinstance(v,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='current' for t in v.targets))
check_index=loop.body.index(check);comparison=loop.body[check_index+1]
ns['r']=dict(current_input_path='/selected/input',current_stat=prior,current_input_size='25')
def run(current):
 ns['stat']=lambda path:current
 exec(compile(ast.Module(body=[check,comparison],type_ignores=[]),'<actual-input-stat-check>','exec'),ns)
run(expected)
for key in expected:
 changed=dict(expected);changed[key]+=1
 if key=='dev':run(changed);continue
 try:run(changed)
 except SystemExit:pass
 else:raise AssertionError('stat drift accepted: '+key)
ns['stopped']=False;assert ns['current_witness_stat'](prior,'/foreign')==prior
print('PASS actual stat boundary: reviewed historical66→68 facts retained; portabledev/all other fields/terminal strict')
PY_STAT
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-legacy-audit.XXXXXX")"; TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
TREE="${TMP}/snapshot/tree"; SCRATCH="${TMP}/scratch"; INPUTS="${SCRATCH}/_ecoda_inputs/audit"
mkdir -p "${TREE}/src/utils/bash" "${TREE}/src/notebooks" "${INPUTS}" "${TMP}/bin"
cp "${ROOT}/src/utils/bash/legacy_analysis_audit.sh" "${ROOT}/src/utils/bash/ecoda_run_common.sh" "${TREE}/src/utils/bash/"
# Only the test-owned copy simulates the selected BeeGFS mount on macOS.
"${ROOT}/.pixi/envs/default/bin/python" - "${TREE}/src/utils/bash/legacy_analysis_audit.sh" "${SCRATCH}" <<'PY_MOUNT'
import pathlib,sys
p=pathlib.Path(sys.argv[1]);text=p.read_text().replace("pathlib.Path('/proc/self/mountinfo').read_text()",repr('10 1 0:68 / '+sys.argv[2]+' rw - beegfs beegfs_nodev rw\n'))
# The fixture worker simulates a different kernel's device number; every
# inode, byte, mode and timestamp remains the same real temporary file.
text=text.replace('s.st_dev',"s.st_dev+(1000 if os.environ.get('LEGACY_TEST_PINNED_WORKER')=='1' else 0)")
p.write_text(text)
PY_MOUNT
cp "${ROOT}/src/notebooks/audit_legacy_analysis_inputs.R" "${TREE}/src/notebooks/"
cp "${ROOT}/datasets.json" "${TREE}/datasets.json"
printf 'immutable fixture source\n' > "${TREE}/source.manifest"
cat > "${TREE}/src/slurm_config.sh" <<'SH'
export PROJECT_ROOT="${TEST_TREE}" HPC_SCRATCH_DIR="${TEST_SCRATCH}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners"
export PYTHON_BIN="${TEST_PYTHON}" PIXI_RSCRIPT="${TEST_RSCRIPT}"
SLURM_PARTITION=shared-cpu
SH
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'SH'
ecoda_runtime_validate_submission() { [[ "$1" == apptainer ]]; printf 'pinned runtime\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run() { [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]]; cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"; }
ecoda_runtime_export_csv() { printf ECODA_RUNTIME_PROFILE=stage5; }
ecoda_runtime_reexec_worker() {
  [[ "$1" == stage5 && "$2" == "${TEST_TREE}/src/utils/bash/legacy_analysis_audit.sh" && "$3" == --worker ]]
  export LEGACY_TEST_PINNED_WORKER=1
  echo PINNED_WORKER
}
SH
cat > "${TMP}/bin/sbatch" <<'SH'
#!/bin/bash
set -euo pipefail
[[ "$*" == *'--parsable --wait'* && "$*" == *'--cpus-per-task=1'* && "$*" == *'--mem=64G'* && "$*" == *'--time=01:00:00'* && "$*" != *'--gpus'* ]]
printf '%s\n' "$*" >> "${TEST_SCRATCH}/scheduler.calls"
args=("$@")
set +e
# Slurm runs a copied spool script with no configuration siblings.
spool="${TEST_SCRATCH}/spool/${ECODA_RUN_ID}"
mkdir -p "${spool}"
cp "${args[${#args[@]}-2]}" "${spool}/slurm_script"
SLURM_JOB_ID=765432 bash "${spool}/slurm_script" "${args[${#args[@]}-1]}" > "${ECODA_RUN_ROOT}/logs/legacy_analysis_765432.log" 2> "${ECODA_RUN_ROOT}/logs/legacy_analysis_765432.err"
rc=$?
set -e
[[ ${rc} -eq 0 ]] || cat "${ECODA_RUN_ROOT}/logs/legacy_analysis_765432.err" >&2
printf '765432;fixture\n'
exit "${rc}"
SH
cat > "${TMP}/bin/Rscript" <<'SH'
#!/bin/bash
exec "${TEST_PYTHON}" "${TEST_TMP}/emit_report.py" "$@"
SH
for name in sbatch Rscript; do chmod +x "${TMP}/bin/${name}"; done
for name in sacct squeue rsync; do printf '#!/bin/bash\necho forbidden >> "%s"\nexit 99\n' "${TMP}/forbidden" > "${TMP}/bin/${name}"; chmod +x "${TMP}/bin/${name}"; done
# Read exact historical producer blobs; never execute them in the fixture.
OLD_SOURCE="${SCRATCH}/_ecoda_source_snapshots/native-9ca/tree"
mkdir -p "${OLD_SOURCE}/src/utils/bash"
for name in h5ad_preflight_worker.sh ecoda_run_common.sh; do
  git -C "${ROOT}" show "9ca2453eea7320a721ef6dfd458f9db108abecac:src/utils/bash/${name}" > "${OLD_SOURCE}/src/utils/bash/${name}"
done
export LEGACY_TEST_OLD_SOURCE="${OLD_SOURCE}"
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_TMP="${TMP}"
export TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_RSCRIPT="${TMP}/bin/Rscript"
export PATH="${TMP}/bin:${PATH}" ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
"${TEST_PYTHON}" - "${INPUTS}" "${SCRATCH}" <<'PY'
import copy,csv,hashlib,json,os,pathlib,sys
p,scratch=map(pathlib.Path,sys.argv[1:]);parent=scratch/'_ecoda_runs/old_preflight';m=parent/'manifests';m.mkdir(parents=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def stat(p):
 s=p.stat();return dict(dev=s.st_dev,ino=s.st_ino,size=s.st_size,mode=s.st_mode&0o7777,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)
def evidence(path,text):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text);return dict(path=str(path),sha256=sha(path),stat=stat(path),text=text)
def save(path,x):path.write_text(json.dumps(x));return dict(path=str(path),sha256=sha(path))
cohorts=['Adams','Bassez','Gongsharma_cmv_young_males','Kfoury','Kim','Lee','Pelka','Smillie','Wu','Zhang']
source=evidence(m/'source.manifest','SOURCE_COMMIT=9ca2453eea7320a721ef6dfd458f9db108abecac\nSOURCE_ROOT='+os.environ['LEGACY_TEST_OLD_SOURCE']+'\nDATASETS_SHA256='+sha(pathlib.Path(os.environ['TEST_TREE'])/'datasets.json')+'\n');runtime=evidence(m/'runtime.identity','prior pinned runtime\n')
w=[];files=[];rows=[];pairs=[]
for n,ds in enumerate(cohorts):
 task=2 if ds=='Bassez' else 3 if ds=='Adams' else n+2
 path=scratch/ds/'input.h5ad';path.parent.mkdir();path.write_bytes(b'fixture H5AD content: no reader permitted in wrapper');os.chmod(path,0o600)
 md5=hashlib.md5(path.read_bytes()).hexdigest();size=str(path.stat().st_size);proof=evidence(p/(ds+'.md5'),f'MD5={md5}\nSIZE={size}\nPATH={path}\n')
 record=evidence(m/'artifacts'/(hashlib.sha256(str(path).encode()).hexdigest()[:32]+'.record'),f'PATH={path}\nSIZE={size}\nMD5={md5}\nRUN_ID=old_preflight\nPRODUCER=stage5_preflight\nSTATE=PUBLISHED\n')
 status=evidence(parent/'status/h5ad_preflight'/(ds+'__benchmark_analysis.status'),f'STATE=OK\nRUN_ID=old_preflight\nDATASET={ds}\nVIEW=benchmark_analysis\nTASK_ID={task}\n')
 log=evidence(parent/'logs'/f'h5ad_preflight_stage5_4665493_{task}.log','')
 w.append(dict(dataset=ds,view='benchmark_analysis',current_input_path=str(path),current_input_md5=md5,current_input_size=size,current_input_proof=proof['path'],current_stat=stat(path),fullhash_array_id='4665493',task_id=task,no_change_timestamps_since_record=True,proof=proof,record=record,status=status,log=log))
 artifacts=[]
 for j in range(11 if n<6 else 10):
  f=p/f'{ds}_{j}.rds';f.write_text('old result fixture');artifacts.append(f);files.append(dict(dataset=ds,artifact_path=str(f),md5=hashlib.md5(f.read_bytes()).hexdigest(),size=str(f.stat().st_size),needed_bundle_keys=''))
 for j in range(21 if n<9 else 20):
  f=artifacts[j%len(artifacts)];method=f'M{j}'
  rows.append(dict(dataset=ds,view='benchmark_analysis',method_key=method,kind='standalone_result',artifact_path=str(f),bundle_key='',sample_id_col='',md5=hashlib.md5(f.read_bytes()).hexdigest(),size=str(f.stat().st_size),current_input_path=str(path),current_input_md5=md5,current_input_size=size,current_input_proof=proof['path'],origin='legacy_reuse',status='BLOCKED_UNTIL_SCIENTIFIC_REVIEW'))
 if n<9:pairs.append(dict(dataset=ds,author_path=str(artifacts[0]),author_md5=rows[-1]['md5'],author_bundle_key='ECODA_authors_HR',metadata_path=str(artifacts[1]),metadata_md5=rows[-1]['md5'],current_input_path=str(path),current_input_md5=md5))
def tsv(name,values):
 with (p/name).open('w') as f:
  out=csv.DictWriter(f,fieldnames=list(values[0]),delimiter='\t',lineterminator='\n');out.writeheader();out.writerows(values)
for name,values in [('rows.209.tsv',rows),('files.106.tsv',files),('author_hr_metadata_pairs.9.tsv',pairs),('avg_legacy.10.tsv',[dict(dataset=x) for x in cohorts]),('zeroimp_nas_counterparts.9.tsv',[dict(dataset=x) for x in cohorts[:9]])]:tsv(name,values)
science=dict(schema='ecoda_p6_legacy_science_audit_proposal_v4',total_logical_rows=209,physical_files=106,author_hr_metadata_pairs=9,scientifically_approved_rows=0,controls={n:sha(p/n) for n in ('rows.209.tsv','files.106.tsv','author_hr_metadata_pairs.9.tsv','avg_legacy.10.tsv','zeroimp_nas_counterparts.9.tsv')})
preflight='Stephenson\tbenchmark_analysis\t'+str(scratch/'unselected_Stephenson.h5ad')+'\n'+''.join(r['dataset']+'\tbenchmark_analysis\t'+r['current_input_path']+'\n' for r in sorted(w,key=lambda r:r['task_id']))
input_manifest=evidence(m/'h5ad_preflight.tsv',preflight)
witness=save(p/'witness.json',dict(schema='ecoda_p6_preflight_witness_observation_v1',parent_run_id='old_preflight',parent_gate_id='old_failed_gate',scientifically_approved=False,terminal_accounting_review_pending=True,h5ad_payload_bytes_read=0,source=source,runtime=runtime,input_manifest=input_manifest,rows=w))
account=evidence(p/'accounting.tsv','4665493|COMPLETED|0:0\n999999|FAILED|1:0\n')
# Consumed gate fields and the complete accounting shape are anchored to
# .gate/ecoda_stage5_main_confounds_20261109_v1{,.inspect}.json (one inspect).
# Only IDs, fixture evidence and states are synthetic; no .gate runtime input.
def audit_capture(text):
 parsed=[dict(job_id=v[0],state=v[1],exit_code=v[2],line=n) for n,line in enumerate(text.splitlines(),1) for v in [line.split('|')]]
 failed=[dict(r,scheduler_id=r['job_id']) for r in parsed if r['state']!='COMPLETED' or r['exit_code']!='0:0']
 accounting=dict(actual_exit=0,attempted=True,command='sacct -n -P -X -j 4665493,999999 --format=JobIDRaw,State,ExitCode',failed_rows=failed,format='JobIDRaw|State|ExitCode',malformed_rows=[],matches={key:[r for r in parsed if r['job_id']==key] for key in ('4665493','999999')},missing_scheduler_ids=[],ok=not failed,query_count=1,rows=parsed,scheduler_ids=['4665493','999999'],stdout=text)
 return dict(accounting=accounting,artifact_contracts=[],commands=[],immutable_fingerprints=[],passed=False,performed_at='2026-11-09T00:00:00Z')
audit=audit_capture(account['text']);reviewer=dict(approved=False,approved_at=None,required=True,reviewer=None)
gate=save(p/'parent_gate.json',dict(schema_version='1.0',gate_id='old_failed_gate',state='FAILED',scheduler_array_ids=['4665493','999999'],scheduler_watchdog_ids=[],audit=audit,audit_completed=True,audit_state='COMPLETED',reviewer=reviewer))
inspect=save(p/'inspect.json',dict(operation='inspect',gate_id='old_failed_gate',state='FAILED',ok=True,release_eligible=False,reviewer=reviewer,audit=audit))
review=save(p/'review.json',dict(schema='ecoda_p6_preflight_witness_review_v1',decision='APPROVED_INPUT_CONTENT_WITNESS_ONLY',review_complete=True,scientifically_approved_rows=0,formal_pipeline_dependency=False,parent_terminal_state='FAILED',observation=witness,gate_manifest=gate,terminal_inspect=inspect,accounting={k:account[k] for k in ('path','sha256')},selected_preflight_tasks=[dict(**{k:r[k] for k in ('dataset','fullhash_array_id','task_id','current_input_path','current_input_md5','current_input_size')},status='PASS') for r in w]))
request=dict(schema='ecoda_legacy_analysis_audit_gate_v1',science_control=save(p/'science_control.json',science),observation=save(p/'observation.json',dict(status='OBSERVED')),zeroimp_schema=save(p/'zero_schema.json',{}),zeroimp_policy=save(p/'zero_policy.json',{}),feather_headers=save(p/'headers.json',{}),input_review=review)
save(p/'request.json',request)
for state in ('COMPLETED','SUCCEEDED'):
 g=json.loads((p/'parent_gate.json').read_text());g['state']=state
 i=json.loads((p/'inspect.json').read_text());i['state']=state;i['audit']=audit_capture(i['audit']['accounting']['stdout'].replace('999999|FAILED|1:0','999999|COMPLETED|0:0'));g['audit']=i['audit']
 a=p/(state+'_accounting.tsv');a.write_text(i['audit']['accounting']['stdout'])
 r=json.loads((p/'review.json').read_text());r['parent_terminal_state']=state;r['gate_manifest']=save(p/(state+'_gate.json'),g);r['terminal_inspect']=save(p/(state+'_inspect.json'),i);r['accounting']=dict(path=str(a),sha256=sha(a))
 q=dict(request);q['input_review']=save(p/(state+'_review.json'),r);save(p/(state+'_request.json'),q)
pending=json.loads((p/'review.json').read_text());pending['review_complete']=False
request['input_review']=save(p/'pending_review.json',pending);save(p/'pending_request.json',request)
wrong=copy.deepcopy(json.loads((p/'witness.json').read_text()));wrong['rows'][0]['task_id']=2
r=json.loads((p/'review.json').read_text());r['observation']=save(p/'wrong_task_witness.json',wrong);r['selected_preflight_tasks'][0]['task_id']=2
q=json.loads((p/'request.json').read_text());q['input_review']=save(p/'wrong_task_review.json',r);save(p/'wrong_task_request.json',q)

# Stopped parent is a frozen historical snapshot, never a fabricated inspection.
g=json.loads((p/'parent_gate.json').read_text());g.update(state='PRELAUNCH_STOP',terminal_status=None,remote_exit_code=None,audit_completed=False,audit_state='NOT_STARTED');g.pop('audit');g.pop('terminal_status')
label_rows=[]
config=json.loads((pathlib.Path(os.environ['TEST_TREE'])/'datasets.json').read_text())
def full_stat(value):
 return {'st_'+k:v for k,v in value.items() if k!='mode'}|{'st_mode':0o100000|value['mode']}
for row in w:
 ds=row['dataset'];columns=config[ds]['columns']
 label_rows.append(dict(dataset=ds,view='benchmark_analysis',path=row['current_input_path'],expected_source_md5=row['current_input_md5'],size=int(row['current_input_size']),canonical_sample_column='Sample',configured_sample_column=columns['sample'],label_column=columns['label'],sample_ids=['s2','s1'],named_labels={'s1':'B','s2':'A'},ordered_primary_labels=['A','B'],observed_first_cell_sample_ids=['s2','s1'],source_identity_order_matches_observed=True,n_samples=2,stat_before=full_stat(row['current_stat']),stat_after=full_stat(row['current_stat']),preflight_record=dict(path=row['record']['path'],sha256=row['record']['sha256'],stat=full_stat(row['record']['stat'])),preflight_status={k:row['status'][k] for k in ('path','sha256')},proof={k:row['proof'][k] for k in ('path','sha256')}))
label_rows.append(dict(dataset='Stephenson')) # Unselected labels are not consumed.
labels=dict(schema='ecoda_old56_current_benchmark_labels_observation_v1',scientific_reuse_approved=False,expression_values_read=False,full_h5ad_hashes=0,fits=0,scores_computed=0,normalization_policy='exact_text_storage_v2_no_numeric_or_logical_conversion',original_source_commit='9ca2453eea7320a721ef6dfd458f9db108abecac',collector_sha256='a'*64,rows=label_rows)
label_ref=save(p/'labels.json',labels)
label_review=dict(schema='ecoda_current_primary_label_witness_review_v1',status='APPROVED_CURRENT_PRIMARY_LABEL_WITNESS',scientific_result_reuse_approved=False,observation_sha256=label_ref['sha256'],collector_sha256=labels['collector_sha256'],datasets=11)
r=dict(reviewed_device_transition=None,reviewer='fixture input-only reviewer',reviewed_at='2026-11-09T22:00:00Z',schema='ecoda_p6_stopped_parent_input_review_v1',decision='APPROVED_PUBLISHED_INPUT_CONTENT_WITNESS_ONLY',review_complete=True,scientifically_approved_rows=0,formal_pipeline_dependency=False,parent_gate_approved=False,scheduler_exit_established=False,method_science_approved=False,parent_gate_state='PRELAUNCH_STOP',completion_basis='PUBLISHED_FULL_HASH_RECORD_AND_NATIVE_TASK_OK',observation=witness,gate_manifest=save(p/'stopped_gate.json',g),current_labels=label_ref,current_labels_review=save(p/'label_review.json',label_review),selected_preflight_tasks=json.loads((p/'review.json').read_text())['selected_preflight_tasks'])
for key,name in [('native_preflight_worker','h5ad_preflight_worker.sh'),('native_checksum_helper','ecoda_run_common.sh')]:
 path=pathlib.Path(os.environ['LEGACY_TEST_OLD_SOURCE'])/'src/utils/bash'/name;r[key]=dict(path=str(path),sha256=sha(path))
def stopped_request(name,review):
 q=json.loads((p/'request.json').read_text());q['input_review']=save(p/(name+'_review.json'),review);save(p/(name+'_request.json'),q)
stopped_request('stopped',r)
for name,key,value in [('stopped_pending','review_complete',False),('stopped_science','method_science_approved',True),('stopped_exit_claim','scheduler_exit_established',True),('stopped_unknown','schema','unknown'),('stopped_accounting','accounting',{}),('stopped_dependency','formal_pipeline_dependency',True)]:
 bad=copy.deepcopy(r);bad[key]=value;stopped_request(name,bad)
for name,key,value in [('stopped_exit','remote_exit_code',1),('stopped_success','state','COMPLETED'),('stopped_audited','audit_completed',True),('stopped_terminal','terminal_status',{'state':'FAILED'})]:
 bad=copy.deepcopy(r);gate=copy.deepcopy(g);gate[key]=value;bad['gate_manifest']=save(p/(name+'_gate.json'),gate);stopped_request(name,bad)
for name,key in [('stopped_source','source'),('stopped_runtime','runtime'),('stopped_record','record'),('stopped_proof','proof'),('stopped_status','status')]:
 bad=copy.deepcopy(r);wx=json.loads((p/'witness.json').read_text());target=wx[key] if key in ('source','runtime') else wx['rows'][0][key];target['text']+='TAMPER\n';target['sha256']=hashlib.sha256(target['text'].encode()).hexdigest();bad['observation']=save(p/(name+'_witness.json'),wx);stopped_request(name,bad)
bad=copy.deepcopy(r);bad['native_preflight_worker']['sha256']='b'*64;stopped_request('stopped_producer',bad)
bad=copy.deepcopy(r);wx=json.loads((p/'witness.json').read_text());wx['rows'][0]['task_id']=2;bad['observation']=save(p/'stopped_task_witness.json',wx);bad['selected_preflight_tasks'][0]['task_id']=2;stopped_request('stopped_task',bad)
for name,change in [('stopped_labels',lambda row:row['ordered_primary_labels'].reverse()),('stopped_stat',lambda row:row['stat_after'].__setitem__('st_ino',999)),('stopped_label_proof',lambda row:row['proof'].__setitem__('sha256','c'*64))]:
 bad=copy.deepcopy(r);lx=copy.deepcopy(labels);change(lx['rows'][0]);bad['current_labels']=save(p/(name+'_labels.json'),lx);lr=copy.deepcopy(label_review);lr['observation_sha256']=bad['current_labels']['sha256'];bad['current_labels_review']=save(p/(name+'_label_review.json'),lr);stopped_request(name,bad)
PY
cat > "${TMP}/emit_report.py" <<'PY'
import csv,hashlib,json,os,pathlib,sys
assert os.environ.get('LEGACY_TEST_PINNED_WORKER')=='1'
assert sys.argv[1]=='--vanilla' and sys.argv[2].endswith('/src/notebooks/audit_legacy_analysis_inputs.R')
a=sys.argv[3:];v=dict(zip(a[::2],a[1::2]));assert len(v)==17 and len(a)==34
m=pathlib.Path(v['--control']).parent
rows=list(csv.DictReader((m/'rows.209.tsv').open(),delimiter='\t'))
failed=os.environ.get('LEGACY_TEST_FAIL')=='1'
x=dict(schema='ecoda_legacy_analysis_compatibility_v1',physical_files=106,logical_rows=209,rows={r['dataset']+'/'+r['method_key']:{} for r in rows},author_hr_pairs={str(i):{} for i in range(9)},cohorts={str(i):{} for i in range(10)},scientifically_approved_rows=0,model_fits=0,score_calculations=0,h5ad_expression_values_read=False,h5ad_content_rehashed_here=False,prior_validated_content_witness_required=True,validation_identities={k:w for k,w in v.items() if k!='--report'},compatible=not failed,status='COMPATIBILITY_MISMATCHES' if failed else 'COMPATIBILITY_EVIDENCE_COMPLETE_PENDING_REVIEW')
with open(v['--report'],'x') as f:json.dump(x,f)
print('startup noise is a log, never report content')
sys.exit(1 if failed else 0)
PY
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
invoke() { ECODA_RUN_ID="$1" ECODA_RUN_ROOT="${SCRATCH}/_ecoda_runs/$1" bash "${TREE}/src/utils/bash/legacy_analysis_audit.sh" --control "$2" --control-sha256 "$(ecoda_sha256_file "$2")"; }
reject() { if "$@" > "${TMP}/negative.log" 2>&1; then echo 'expected legacy audit rejection' >&2; exit 1; fi; }
protected_before="$("${TEST_PYTHON}" - "${SCRATCH}" <<'PY'
import hashlib,pathlib,sys
p=pathlib.Path(sys.argv[1]);files=sorted([*p.glob('*/input.h5ad'),*p.glob('_ecoda_runs/old_preflight/**/*')]);print([(str(f),hashlib.sha256(f.read_bytes()).hexdigest(),f.stat().st_mode) for f in files if f.is_file()])
PY
)"
invoke good "${INPUTS}/request.json" > "${TMP}/good.log"
run="${SCRATCH}/_ecoda_runs/good"; report="${run}/artifacts/evidence/legacy_analysis_compatibility.json"
grep -q '^STATE=OK$' "${run}/status/terminal"
grep -q '^LEGACY_ANALYSIS_AUDIT_JOB_ID=765432$' "${TMP}/good.log"
[[ "$(wc -l < "${run}/manifests/selection.tsv" | tr -d ' ')" == 209 ]]
ecoda_validate_artifact_record "${report}" legacy_analysis_compatibility good >/dev/null
ECODA_RUN_ID=good ECODA_RUN_ROOT="${run}" bash "${TREE}/src/utils/bash/legacy_analysis_audit.sh" --validate > "${TMP}/validated.log"
cp "${run}/status/terminal" "${TMP}/terminal.saved"
sed 's/STATE=OK/STATE=FAIL/' "${TMP}/terminal.saved" > "${run}/status/terminal"
before="$(ecoda_sha256_file "${run}/status/terminal")"
ECODA_RUN_ID=good ECODA_RUN_ROOT="${run}" reject bash "${TREE}/src/utils/bash/legacy_analysis_audit.sh" --validate
[[ "$(ecoda_sha256_file "${run}/status/terminal")" == "${before}" ]]
cp "${TMP}/terminal.saved" "${run}/status/terminal"
invoke completed "${INPUTS}/COMPLETED_request.json" > "${TMP}/completed.log"
grep -q '^STATE=OK$' "${SCRATCH}/_ecoda_runs/completed/status/terminal"
invoke stopped "${INPUTS}/stopped_request.json" > "${TMP}/stopped.log"
stopped_run="${SCRATCH}/_ecoda_runs/stopped"
[[ ! -e "${stopped_run}/manifests/parent_inspect.json" && ! -e "${stopped_run}/manifests/parent_accounting.tsv" ]]
"${TEST_PYTHON}" - "${stopped_run}/manifests/input_basis.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]));assert x['parent_snapshot_state']=='PRELAUNCH_STOP' and x['parent_snapshot_terminal_status'] is None and x['parent_snapshot_remote_exit_code'] is None
assert x['scheduler_exit_established'] is False and x['method_science_approved'] is False and x['parent_gate_approved'] is False and x['formal_pipeline_dependency'] is False and x['scientifically_approved_rows']==0
PY
ECODA_RUN_ID=stopped ECODA_RUN_ROOT="${stopped_run}" bash "${TREE}/src/utils/bash/legacy_analysis_audit.sh" --validate > "${TMP}/stopped-validated.log"
cp "${stopped_run}/status/terminal" "${TMP}/stopped-terminal.saved"
sed 's/STATE=OK/STATE=FAIL/' "${TMP}/stopped-terminal.saved" > "${stopped_run}/status/terminal"
before="$(ecoda_sha256_file "${stopped_run}/status/terminal")"
ECODA_RUN_ID=stopped ECODA_RUN_ROOT="${stopped_run}" reject bash "${TREE}/src/utils/bash/legacy_analysis_audit.sh" --validate
[[ "$(ecoda_sha256_file "${stopped_run}/status/terminal")" == "${before}" ]]
cp "${TMP}/stopped-terminal.saved" "${stopped_run}/status/terminal"
for name in pending science exit_claim unknown accounting dependency exit success audited terminal source runtime record proof status producer task labels stat label_proof; do
  reject invoke "stopped_${name}" "${INPUTS}/stopped_${name}_request.json"
done
reject invoke invented_state "${INPUTS}/SUCCEEDED_request.json"
reject invoke wrong_task "${INPUTS}/wrong_task_request.json"
grep -q 'preflight task does not match its exact physical manifest row' "${TMP}/negative.log"
reject invoke pending "${INPUTS}/pending_request.json"
grep -q '^STATE=FAIL$' "${SCRATCH}/_ecoda_runs/pending/status/terminal"
LEGACY_TEST_FAIL=1 reject invoke mismatch "${INPUTS}/request.json"
failed_run="${SCRATCH}/_ecoda_runs/mismatch"
grep -q '^LEGACY_ANALYSIS_AUDIT_JOB_ID=765432$' "${TMP}/negative.log"
grep -q '^STATE=FAIL$' "${failed_run}/status/terminal"
[[ -f "${failed_run}/artifacts/evidence/legacy_analysis_compatibility.json" && "$(cat "${failed_run}/manifests/scheduler_ids.tsv")" == $'ARRAY\t765432' ]]
protected_after="$("${TEST_PYTHON}" - "${SCRATCH}" <<'PY'
import hashlib,pathlib,sys
p=pathlib.Path(sys.argv[1]);files=sorted([*p.glob('*/input.h5ad'),*p.glob('_ecoda_runs/old_preflight/**/*')]);print([(str(f),hashlib.sha256(f.read_bytes()).hexdigest(),f.stat().st_mode) for f in files if f.is_file()])
PY
)"
[[ "${protected_before}" == "${protected_after}" ]]
printf changed >> "${SCRATCH}/Adams/input.h5ad"
reject invoke changed "${INPUTS}/request.json"
reject invoke stopped_changed "${INPUTS}/stopped_request.json"
[[ "$(wc -l < "${SCRATCH}/scheduler.calls" | tr -d ' ')" == 4 && ! -e "${TMP}/forbidden" ]]
echo 'PASS: exact finite legacy scope, terminal/stopped typed input-content review/stat/label/source witness, pinned CPU job, failed ID/mismatch report preservation, read-only terminal validation and immutable original sources'
