#!/bin/bash
# Real observer/common records; only Slurm/runtime are stubbed. Metadata only.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-orphaned-observer.XXXXXX")"; TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
TREE="${TMP}/snapshot/tree"; SCRATCH="${TMP}/scratch"; INPUTS="${SCRATCH}/_ecoda_inputs/observer"
mkdir -p "${TREE}/src/utils/bash" "${INPUTS}" "${TMP}/bin"
cp "${ROOT}/src/utils/bash/orphaned_main_observer.sh" "${ROOT}/src/utils/bash/ecoda_run_common.sh" "${TREE}/src/utils/bash/"
printf 'pinned fixture source\n' > "${TREE}/source.manifest"
cat > "${TREE}/src/slurm_config.sh" <<'BASH'
export PROJECT_ROOT="${TEST_TREE}" HPC_SCRATCH_DIR="${TEST_SCRATCH}" PYTHON_BIN="${TEST_PYTHON}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners"
SLURM_PARTITION=shared-cpu
BASH
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'BASH'
ecoda_runtime_validate_submission() { [[ "$1" == apptainer ]] || return 1; printf 'RUNTIME_IMAGE_SHA256=%064d\nRUNTIME_MANIFEST_SHA256=%064d\n' 1 2 > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run() { cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"; }
ecoda_runtime_export_csv() { [[ "$1" == stage5 && "$2" == 0 ]] || return 1; printf ECODA_RUNTIME_PROFILE=stage5; }
ecoda_runtime_reexec_worker() { [[ "$1" == stage5 && "$2" == "${TEST_TREE}/src/utils/bash/orphaned_main_observer.sh" && "$3" == --worker ]] || return 1; echo PINNED_OBSERVER_WORKER; }
BASH
cat > "${TMP}/bin/sbatch" <<'BASH'
#!/bin/bash
set -euo pipefail
[[ "$*" == *'--parsable --wait --dependency=afterany:4665568:4665570:4665612 --kill-on-invalid-dep=yes'* && "$*" == *'--cpus-per-task=1 --mem=2G --time=00:05:00'* && "$*" != *'--array'* && "$*" != *'--gpus'* ]] || exit 99
printf '%s\n' "$*" >> "${TEST_TMP}/sbatch.calls"
args=("$@")
rc=1
if [[ "${OBSERVER_TEST_FAIL:-0}" != 1 ]]; then
  mkdir -p "${TEST_TMP}/spool"
  cp "${args[${#args[@]}-2]}" "${TEST_TMP}/spool/slurm_script"
  set +e
  SLURM_JOB_ID=765432 bash "${TEST_TMP}/spool/slurm_script" "${args[${#args[@]}-1]}" > "${ECODA_RUN_ROOT}/logs/observer.log" 2> "${ECODA_RUN_ROOT}/logs/observer.err"
  rc=$?
  set -e
fi
printf '765432;fixture\n'
exit "${rc}"
BASH
chmod +x "${TMP}/bin/sbatch"
for command in sacct squeue scontrol rsync; do printf '#!/bin/bash\necho forbidden >> "%s"\nexit 99\n' "${TMP}/forbidden" > "${TMP}/bin/${command}"; chmod +x "${TMP}/bin/${command}"; done
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_TMP="${TMP}" TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python"
export ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1 PATH="${TMP}/bin:${PATH}"
"${TEST_PYTHON}" - "${INPUTS}" "${SCRATCH}" <<'PY'
import hashlib,json,pathlib,shlex,sys
p,scratch=map(pathlib.Path,sys.argv[1:]);parent=scratch/'_ecoda_runs/stage5_main_confounds_20261109_v2';(parent/'manifests').mkdir(parents=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(path,data):
    path.write_bytes(data);return dict(path=str(path),sha256=sha(data))
def js(path,value):return save(path,(json.dumps(value)+'\n').encode())
source=dict(SOURCE_ROOT='/fixture/old/tree',SOURCE_COMMIT='9ca2453eea7320a721ef6dfd458f9db108abecac')
runtime=dict(RUNTIME_IMAGE='/fixture/runtime.sif',RUNTIME_MANIFEST='/fixture/runtime.sif.manifest',RUNTIME_IMAGE_SHA256=f'{1:064d}',RUNTIME_MANIFEST_SHA256=f'{2:064d}')
source_ref=save(parent/'manifests/source.manifest',''.join(f'{k}={v}\n' for k,v in source.items()).encode())
runtime_ref=save(parent/'manifests/runtime.identity',''.join(f'{k}={v}\n' for k,v in runtime.items()).encode())
argv=['/bin/bash','/fixture/snapshot.sh','exec','--run-id',parent.name,'--scratch-root',str(scratch),'--source-root',source['SOURCE_ROOT'],'--runtime-image',runtime['RUNTIME_IMAGE'],'--runtime-manifest',runtime['RUNTIME_MANIFEST'],'--script','src/5_run_benchmark_methods/stage5_campaign.sh']
gate=dict(gate_id='ecoda_stage5_main_confounds_20261109_v2',state='PRELAUNCH_STOP',profile='durable-hpc-gate-ecoda',serialization_group='ecoda-benchmark',exact_command=shlex.join(argv),launch_timestamp='2026-11-09T08:20:00+00:00',event_generation='a'*32,remote_runner=str(scratch/'gates/main.runner.sh'),remote_status=str(scratch/'gates/main.status.json'),remote_failure_status=str(scratch/'gates/main.status.failure'),completion_channel='file:'+str(scratch/'gates/main.done'))
gate['command_digest']=sha(gate['exact_command'].encode())
children={'benchmark':'stage5_main_benchmark_20261109_v2','uncorrected':'stage5_main_uncorrected_20261109_v2','corrected':'stage5_main_corrected_20261109_v2'}
files=[dict(path=gate[k],exists=False) for k in ('remote_status','remote_failure_status')]+[dict(path=gate['completion_channel'][5:],exists=False)]
rows=[];counter=900000
for ctx,n in (('benchmark',38),('uncorrected',17),('corrected',16)):
    entries=[['ARRAY','4665568'],['WATCHDOG','4665570']] if ctx=='benchmark' else ([['AGGREGATE_GATE','4665594']] if ctx=='uncorrected' else [])
    while len(entries)<n:counter+=1;entries.append(['ARRAY',str(counter)])
    text=''.join('\t'.join(r)+'\n' for r in entries)
    files.append(dict(path=str(scratch/'_ecoda_runs'/children[ctx]/'manifests/scheduler_ids.tsv'),exists=True,text=text,sha256=sha(text.encode())))
    rows.extend((ctx,children[ctx],*r) for r in entries)
for ctx,j in (('benchmark','4665612'),('corrected','4665633')):rows.append((ctx,children[ctx],'AGGREGATE_GATE',j))
active=['4665568','4665570','4665612'];descriptors={}
for j,dep in zip(active,['',' Dependency=afterany:4665568_*(unfulfilled)',' Dependency=afterany:4665570(unfulfilled)']):
    descriptors[j]=dict(argv=['scontrol','show','job','-o',j],exit=0,stdout=f'JobId={j} JobState=PENDING{dep} ')
# An embedded live log snapshot is evidence at observation time, not a promise
# that the old running worker can never append to its own log.
live=parent/'live-worker.log';live.write_text('observed log\n');files.append(dict(path=str(live),exists=True,text=live.read_text(),sha256=sha(live.read_bytes())))
incident=dict(schema='ecoda_login_reboot_operational_observation_v1',hostname='login1.bamboo',boot_utc='2026-11-09T14:47:27+00:00',observed_at='2026-11-09T16:05:11+00:00',no_accounting_query=True,no_scientific_payload_reads=True,no_remote_writes=True,tmux=dict(exit=1,stdout=''),processes=dict(exit=0,stdout='collector only'),files=files,active_job_descriptors=descriptors)
q=dict(schema='ecoda_orphaned_main_observer_v1',original_gate=js(p/'old-gate.json',gate),incident_observation=js(p/'incident.json',incident),scheduler_catalog=save(p/'catalog.tsv',''.join('\t'.join(r)+'\n' for r in rows).encode()),original_source=source_ref,original_runtime=runtime_ref,active_dependencies=[dict(job_id=j,descriptor_sha256=sha(descriptors[j]['stdout'].encode())) for j in active])
js(p/'control.json',q)
for fault in ('gate','dependencies','descriptor','catalog'):
    bad=json.loads(json.dumps(q))
    if fault=='gate':value=dict(gate,state='RUNNING');bad['original_gate']=js(p/'bad-gate.json',value)
    elif fault=='dependencies':bad['active_dependencies'].append(dict(job_id='999999',descriptor_sha256='0'*64))
    elif fault=='descriptor':bad['active_dependencies'][0]['descriptor_sha256']='0'*64
    else:
        badrows=list(rows);badrows[0]=('corrected',children['corrected'],badrows[0][2],badrows[0][3]);bad['scheduler_catalog']=save(p/'bad-catalog.tsv',''.join('\t'.join(r)+'\n' for r in badrows).encode())
    js(p/f'{fault}.control.json',bad)
live.write_text('legitimate later worker progress\n')
PY
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
invoke() { ECODA_RUN_ID="$1" ECODA_RUN_ROOT="${SCRATCH}/_ecoda_runs/$1" bash "${TREE}/src/utils/bash/orphaned_main_observer.sh" --control "$2" --control-sha256 "$(ecoda_sha256_file "$2")"; }
reject() { if "$@" > "${TMP}/negative.log" 2>&1; then echo 'expected observer rejection' >&2; exit 1; fi; }
protected_before="$("${TEST_PYTHON}" - "${SCRATCH}/_ecoda_runs/stage5_main_confounds_20261109_v2" <<'PY'
import hashlib,pathlib,sys
print([(str(p),hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mode) for p in sorted(pathlib.Path(sys.argv[1]).rglob('*')) if p.is_file()])
PY
)"
if ! invoke good "${INPUTS}/control.json" > "${TMP}/good.log" 2>&1; then cat "${TMP}/good.log" >&2; cat "${SCRATCH}/_ecoda_runs/good/logs/observer.err" >&2; exit 1; fi
run="${SCRATCH}/_ecoda_runs/good"; report="${run}/artifacts/evidence/dependency_observation.json"
grep -q '^STATE=OK$' "${run}/status/terminal"
grep -q '^ORPHANED_MAIN_OBSERVER_JOB_ID=765432$' "${TMP}/good.log"
[[ "$(cat "${run}/manifests/scheduler_ids.tsv")" == $'ARRAY\t765432' ]]
[[ $(wc -l < "${run}/manifests/selection.tsv" | tr -d ' ') == 3 ]]
ecoda_validate_artifact_record "${report}" orphaned_main_dependency_observation good >/dev/null
"${TEST_PYTHON}" - "${report}" <<'PY'
import json,pathlib,sys
x=json.loads(pathlib.Path(sys.argv[1]).read_text());assert x['dependency_barrier_observed'] is True and x['original_outcomes']=='UNKNOWN_UNTIL_ORIGINAL_INSPECT' and x['original_wrapper_exit_code'] is None and x['final_emitted_catalog_established'] is False and x['initial_scheduler_catalog_size']==73 and x['scientifically_approved_rows']==x['accounting_queries']==x['model_fits']==x['old_run_writes']==0
PY
ECODA_RUN_ID=good ECODA_RUN_ROOT="${run}" bash "${TREE}/src/utils/bash/orphaned_main_observer.sh" --validate > "${TMP}/validated.log"
cp "${run}/status/terminal" "${TMP}/terminal.saved"
sed 's/STATE=OK/STATE=FAIL/' "${TMP}/terminal.saved" > "${run}/status/terminal"
before="$(ecoda_sha256_file "${run}/status/terminal")"
ECODA_RUN_ID=good ECODA_RUN_ROOT="${run}" reject bash "${TREE}/src/utils/bash/orphaned_main_observer.sh" --validate
[[ "$(ecoda_sha256_file "${run}/status/terminal")" == "${before}" ]]
cp "${TMP}/terminal.saved" "${run}/status/terminal"
for fault in gate dependencies descriptor catalog; do reject invoke "bad_${fault}" "${INPUTS}/${fault}.control.json"; done
OBSERVER_TEST_FAIL=1 reject invoke failed "${INPUTS}/control.json"
grep -q '^ORPHANED_MAIN_OBSERVER_JOB_ID=765432$' "${TMP}/negative.log"
grep -q '^STATE=FAIL$' "${SCRATCH}/_ecoda_runs/failed/status/terminal"
[[ "$(cat "${SCRATCH}/_ecoda_runs/failed/manifests/scheduler_ids.tsv")" == $'ARRAY\t765432' ]]
protected_after="$("${TEST_PYTHON}" - "${SCRATCH}/_ecoda_runs/stage5_main_confounds_20261109_v2" <<'PY'
import hashlib,pathlib,sys
print([(str(p),hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mode) for p in sorted(pathlib.Path(sys.argv[1]).rglob('*')) if p.is_file()])
PY
)"
[[ "${protected_before}" == "${protected_after}" && ! -e "${TMP}/forbidden" && $(wc -l < "${TMP}/sbatch.calls" | tr -d ' ') == 2 ]]
echo 'PASS: one source-pinned 1CPU/2G/5min afterany observer, exact active3/initial73 snapshot, no old writes/queries/fits/outcome claims, readonly validation and failure IDs'
