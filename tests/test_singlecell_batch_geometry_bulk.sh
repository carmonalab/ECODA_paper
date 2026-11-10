#!/bin/bash
# Native manifests/owners/spool/accounting; all scheduler and numerical effects fake.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
"${ROOT}/.pixi/envs/default/bin/python" -B - "${ROOT}" <<'PY'
import hashlib,json,os,pathlib,shutil,subprocess,sys,tempfile
repo=pathlib.Path(sys.argv[1])
with tempfile.TemporaryDirectory(prefix='ecoda-bulk-shell-')as tmp:
 p=pathlib.Path(tmp).resolve();tree=p/'snapshot/tree';u=tree/'src/utils/bash';u.mkdir(parents=True);(tree/'src/notebooks').mkdir()
 for n in ['singlecell_batch_geometry_bulk.sh','ecoda_run_common.sh']:shutil.copy2(repo/'src/utils/bash'/n,u/n)
 (tree/'source.manifest').write_text('SOURCE_ROOT='+str(tree)+'\n')
 (tree/'src/slurm_config.sh').write_text('''export PROJECT_ROOT="${TEST_TREE}" HPC_SCRATCH_DIR="${TEST_SCRATCH}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners" PYTHON_BIN="${TEST_PYTHON}"
SLURM_PARTITION_BENCHMARK_CPU=shared-cpu
''')
 (u/'ecoda_runtime.sh').write_text('''ecoda_runtime_validate_submission(){ printf 'fixture runtime\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run(){ test -s "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_export_csv(){ printf ECODA_RUNTIME_PROFILE=stage5; }
ecoda_runtime_reexec_worker(){ :; }
''')
 (tree/'src/notebooks/run_singlecell_batch_geometry.py').write_text('''import argparse,hashlib,json,os,pathlib
p=argparse.ArgumentParser();p.add_argument('--control');p.add_argument('--control-sha256');p.add_argument('--dataset');p.add_argument('--representation');p.add_argument('--output');p.add_argument('--run-id');p.add_argument('--check-control',action='store_true');p.add_argument('--validate-output',action='store_true');a=p.parse_args()
assert hashlib.sha256(pathlib.Path(a.control).read_bytes()).hexdigest()==a.control_sha256
if a.check_control:raise SystemExit(0)
out=pathlib.Path(a.output)
if a.validate_output:
 assert len((out/'files.tsv').read_text().splitlines())==4;raise SystemExit(0)
assert (a.dataset,a.representation)!=('Breast_cancer','unc_raw')
out.mkdir();files=['embedding.h5ad','metrics.json','display.feather','phases.json']
for name in files:(out/name).write_text(a.dataset+' '+a.representation+' '+name+' immutable')
with (out/'files.tsv').open('x')as f:
 for name in files:f.write(str(out/name)+'\\t'+hashlib.md5((out/name).read_bytes()).hexdigest()+'\\n')
''')
 binpath=p/'bin';binpath.mkdir()
 fake=p/'slurm.py';fake.write_text(r'''import json,os,pathlib,subprocess,sys,time
root=pathlib.Path(os.environ['TEST_STATE']);root.mkdir(exist_ok=True);args=sys.argv[2:];mode=sys.argv[1]
def save(job,rows):
 q=root/(job+'.tmp');q.write_text(rows);q.replace(root/(job+'.state'))
if mode=='wd':
 job=args[0];env=dict(os.environ,SLURM_JOB_ID=job);r=subprocess.run(args[1:],env=env,stdout=open(root/(job+'.log'),'w'),stderr=open(root/(job+'.err'),'w'));save(job,job+'|'+('COMPLETED|0:0'if r.returncode==0 else'FAILED|1:0')+'\n');raise SystemExit(0)
if mode=='sbatch':
 count=root/'counter';n=int(count.read_text())+1 if count.exists()else 900001;count.write_text(str(n));job=str(n)
 flags=[x for x in args if x.startswith('--')];assert '--partition=shared-cpu'in flags and not any(x.startswith(('--constraint','--gpus'))for x in flags)
 script=next(i for i,x in enumerate(args)if x.endswith('singlecell_batch_geometry_bulk.sh'));command=args[script:]
 spool=root/('spool_'+job);shutil=__import__('shutil');shutil.copy2(command[0],spool);command[0]=str(spool)
 if '--watchdog'in command:
  assert any(x.startswith('--dependency=afterany:')for x in flags)
  assert sum(1 for x in root.glob('*.state'))==2,'Both independent arrays must be submitted first'
  subprocess.Popen([sys.executable,__file__,'wd',job,'/bin/bash',*command],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
 else:
  assert '--cpus-per-task=8'in flags and '--time=02:00:00'in flags
  manifest=pathlib.Path(command[command.index('--manifest')+1]);mem=command[command.index('--memory')+1]
  assert '--mem='+mem+'G'in flags and manifest.is_file()and pathlib.Path(str(manifest)+'.md5').is_file()
  (root/(job+'.manifest')).write_text(manifest.name)
  rows=manifest.read_text().splitlines();assert '--array=1-'+str(len(rows))in flags
  outcomes=[]
  for i,line in enumerate(rows,1):
   row=line.split('\t')[0];scenario=os.environ['TEST_SCENARIO'];state=None
   if row=='4'and (mem=='64'or scenario=='ceiling'):state='OUT_OF_MEMORY'
   if row=='2'and (mem=='128'or scenario=='ceiling'):state='OUT_OF_MEMORY'
   if row=='5'and scenario=='nonoom':state='FAILED'
   if state is None:
    e=dict(os.environ,SLURM_ARRAY_TASK_ID=str(i),SLURM_ARRAY_JOB_ID=job,SLURM_CPUS_PER_TASK='8')
    result=subprocess.run(['/bin/bash',*command],env=e,stdout=open(root/(job+'_'+str(i)+'.log'),'w'),stderr=open(root/(job+'_'+str(i)+'.err'),'w'))
    state='COMPLETED'if result.returncode==0 else'FAILED'
   outcomes.append(job+'_'+str(i)+'|'+state+'|'+('0:0'if state=='COMPLETED'else'1:0'))
  save(job,'\n'.join(outcomes)+'\n')
 print(job+';fixture');raise SystemExit(0)
if mode=='sacct':
 job=args[args.index('-j')+1];q=root/(job+'.state');meta=root/(job+'.manifest')
 if meta.exists()and meta.read_text()=='retry_64_1.tsv':assert 'retry_128_1.tsv'in [f.read_text()for f in root.glob('*.manifest')],'Independent retry groups serialized'
 print(q.read_text()if q.exists()else job+'|RUNNING|0:0');raise SystemExit(0)
if mode=='squeue':raise SystemExit(0)
raise AssertionError(mode)
''')
 for n in ['sbatch','sacct','squeue']:
  q=binpath/n;q.write_text('#!/bin/bash\nexec '+sys.executable+' -B '+str(fake)+' '+n+' "$@"\n');q.chmod(0o755)
 scratch=p/'scratch';inputs=scratch/'_ecoda_inputs';inputs.mkdir(parents=True)
 for n in ['receipt.json','review.json']:(inputs/n).write_text('{}')
 control=inputs/'control.json';control.write_text(json.dumps(dict(schema='ecoda_singlecell_batch_main_v2',stage3_receipt={'path':str(inputs/'receipt.json')},stage3_review={'path':str(inputs/'review.json')})))
 common=dict(os.environ,TEST_TREE=str(tree),TEST_SCRATCH=str(scratch),TEST_PYTHON=sys.executable,ECODA_SOURCE_ROOT=str(tree),ECODA_SOURCE_MANIFEST=str(tree/'source.manifest'),ECODA_SOURCE_SNAPSHOT_REQUIRED='1',PATH=str(binpath)+os.pathsep+os.environ['PATH'],ECODA_SINGLECELL_POLL_SECONDS='.02')
 for scenario in ['mixed','nonoom','ceiling']:
  runid='bulk_'+scenario;run=scratch/'_ecoda_runs'/runid;env=dict(common,ECODA_RUN_ID=runid,ECODA_RUN_ROOT=str(run),TEST_SCENARIO=scenario,TEST_STATE=str(p/('slurm_'+scenario)))
  result=subprocess.run(['/bin/bash',str(u/'singlecell_batch_geometry_bulk.sh'),'--control',str(control),'--control-sha256',hashlib.sha256(control.read_bytes()).hexdigest()],env=env,text=True,capture_output=True,timeout=180)
  if (scenario=='mixed')!=(result.returncode==0):
   print(result.stdout,result.stderr)
   for q in pathlib.Path(env['TEST_STATE']).glob('*.err'):print(q.name,q.read_text())
   raise AssertionError((scenario,result.returncode))
  m=run/'manifests';attempts=[s.split('\t')for s in (m/'attempts.tsv').read_text().splitlines()];assert sum(s[0]=='WATCHDOG'for s in attempts)==1
  assert attempts[0][3]=='64'and attempts[1][3]=='128'and attempts[2][0]=='WATCHDOG'
  assert (m/'retry_64_1.tsv').read_text().startswith('4\tCovid19_PBMC\tunc_raw\t64\n')
  assert (m/'retry_128_1.tsv').read_text().startswith('2\tAlzheimer\tcor_harmony\t128\n')
  if scenario!='ceiling':assert not (m/'retry_128_2.tsv').exists(),'No broad/group rerun'
  assert not any(x[3]not in ['0','64','128','256','512']for x in attempts)
  statuses=list((run/'status').glob('row_*.status'));assert len(statuses)==({'mixed':11,'nonoom':10,'ceiling':9}[scenario])
  # Initial successful row never appears in an OOM retry; original owner stays OK.
  assert all(not line.startswith('3\t')for q in m.glob('retry_*.tsv')for line in q.read_text().splitlines())
  original=run/'artifacts/singlecell_batch_geometry/Breast_cancer/cor_harmony/attempt_64/metrics.json';before=original.read_bytes()
  if scenario=='mixed':
   validated=subprocess.run(['/bin/bash',str(u/'singlecell_batch_geometry_bulk.sh'),'--validate'],env=env,capture_output=True,text=True);assert validated.returncode==0,validated.stderr
   # Exercise actual generic purpose branch with native checks/validator. Accounting itself is stubbed only here.
   source=(repo/'src/utils/bash/ecoda_run_audit.sh').read_text();body=source[source.index('_audit_analysis_release() {'):source.index('\n\n\n_audit_scheduler_ids()')]
   audit=p/'audit_function.sh';audit.write_text('source "$1/src/slurm_config.sh"\nsource "$1/src/utils/bash/ecoda_run_common.sh"\n'+body+'''
_audit_metadata_value(){ sed -n "s/^$2=//p" "$1"; }
_audit_metadata_count(){ grep -c "^$2=" "$1"; }
_audit_regular_file(){ test -f "$1" && test ! -L "$1"; }
_audit_scheduler_ids(){ return 0; }
RUN_ROOT_REAL="$ECODA_RUN_ROOT";RUN_ID="$ECODA_RUN_ID";SELECTION_ARG="$RUN_ROOT_REAL/manifests/selection.tsv";AUDIT_TERMINAL_STATE=OK;AUDIT_SOURCE_ROOT="$ECODA_SOURCE_ROOT";SOURCE_MANIFEST_ARG="$ECODA_SOURCE_MANIFEST"
_audit_analysis_release
''')
   r=subprocess.run(['/bin/bash','-eu',str(audit),str(tree)],env=env,capture_output=True,text=True);assert r.returncode==0,r.stderr
   # Wrong-stage owner is rejected read-only.
   owner=next(q for q in (scratch/'_ecoda_owners/artifact').glob('*/owner')if 'PATH='+str(original)+'\n'in q.read_text());owner.write_text(owner.read_text().replace('STAGE=stage5','STAGE=stage3'));bad=owner.read_bytes()
   r=subprocess.run(['/bin/bash',str(u/'singlecell_batch_geometry_bulk.sh'),'--validate'],env=env,capture_output=True,text=True);assert r.returncode!=0 and owner.read_bytes()==bad
  if scenario=='ceiling':
   assert (m/'retry_64_2.tsv').exists()and not (m/'retry_64_3.tsv').exists()
   assert (m/'retry_128_2.tsv').exists()and not (m/'retry_128_3.tsv').exists()
   assert [x[3]for x in attempts if x[2].endswith('retry_128_2.tsv')]==['512']
  if scenario=='nonoom':assert not any(line.startswith('5\t')for q in m.glob('retry_*.tsv')for line in q.read_text().splitlines())
  assert original.read_bytes()==before
  print('PASS',scenario,'all emitted IDs retained; exact OOM rows/memory; valid rows immutable',flush=True)
 print('PASS native bulk shell/generic audit scope/spool/owner/non-OOM/ceiling contracts; no live scheduler')
PY
