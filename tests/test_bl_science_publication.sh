#!/bin/bash
# Actual release adapter/common records; finite fake bytes, runtime/NAS stubs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-bl-publication.XXXXXX")"; TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
TREE="${TMP}/snapshot/tree"; SCRATCH="${TMP}/scratch"; NAS="${TMP}/nas"
mkdir -p "${TREE}/src/utils/bash" "${TREE}/src/utils/py" "${SCRATCH}/_ecoda_inputs" "${SCRATCH}/_ecoda_runs" "${NAS}" "${TMP}/bin"
cp "${ROOT}/src/utils/bash/publish_bl_science.sh" "${ROOT}/src/utils/bash/ecoda_run_common.sh" "${TREE}/src/utils/bash/"
cp "${ROOT}/src/utils/py/bl_science_publication.py" "${TREE}/src/utils/py/"
printf 'SOURCE_ROOT=%s\nSOURCE_COMMIT=%040d\n' "${TREE}" 1 > "${TREE}/source.manifest"
cat > "${TREE}/src/slurm_config.sh" <<'SH'
export PROJECT_ROOT="${TEST_TREE}" HPC_SCRATCH_DIR="${TEST_SCRATCH}" NAS_TARGET_DIR="${TEST_NAS}"
export PYTHON_BIN="${TEST_PYTHON}"
SH
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'SH'
ecoda_runtime_validate_submission() {
  [[ "$1" == apptainer && "${ECODA_RUNTIME_PROFILE}" == stage5 ]]
  printf 'RUNTIME_IMAGE_SHA256=new-bound-runtime\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"
}
ecoda_runtime_validate_bound_run() {
  grep -q '^RUNTIME_IMAGE_SHA256=new-bound-runtime$' "${ECODA_RUN_ROOT}/manifests/runtime.identity"
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
}
SH
cat > "${TREE}/src/utils/bash/ecoda_nas_backup.sh" <<'SH'
#!/bin/bash
set -euo pipefail
[[ "$1" == preflight && "$2" == --nas-root && "$3" == "${TEST_NAS}" && "$4" == --report ]]
[[ "${TEST_MISSING_MOUNT:-0}" != 1 ]] || exit 7
printf '{"mount_source":"fixture-mounted-NAS","nas_root":"%s"}\n' "${TEST_NAS}" > "$5"
SH
cat > "${TMP}/bin/python" <<'SH'
#!/bin/bash
set -euo pipefail
if [[ "${TEST_COPY_FAILURE:-0}" == 1 && " ${*} " == *' copy '* && " ${*} " == *" --destination ${TEST_NAS}/"* ]]; then exit 7; fi
exec "${REAL_PYTHON}" "$@"
SH
cat > "${TMP}/bin/chmod" <<'SH'
#!/bin/bash
set -euo pipefail
if [[ "${TEST_MODE_FAILURE:-0}" == 1 && "$1" == 666 && "$2" == "${TEST_NAS}/"* ]]; then exit 0; fi
exec /bin/chmod "$@"
SH
for command in sbatch sacct squeue Rscript; do
  printf '#!/bin/bash\nprintf forbidden >> "%s"\nexit 99\n' "${TMP}/forbidden.called" > "${TMP}/bin/${command}"
done
chmod +x "${TMP}/bin/"*
export REAL_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_PYTHON="${TMP}/bin/python"
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_NAS="${NAS}" TEST_TMP="${TMP}"
export ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export PATH="${TMP}/bin:${PATH}"
"${REAL_PYTHON}" - <<'PY'
import hashlib,json,os,pathlib
t=pathlib.Path(os.environ['TEST_TMP']);scratch=pathlib.Path(os.environ['TEST_SCRATCH'])
old=scratch/'_ecoda_runs/stage5_benchmark_bl10_20261109_v2';validation=scratch/'_ecoda_runs/science_validation'
for r in (old,validation):(r/'manifests/artifacts').mkdir(parents=True)
(old/'logs').mkdir();(validation/'artifacts/evidence').mkdir(parents=True)
def write(p,text='fixture\n'):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text);return p
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def proof(p):return dict(path=str(p),sha256=sha(p))
def record(p,run,producer):
 md5=hashlib.md5(p.read_bytes()).hexdigest();size=str(p.stat().st_size);value=dict(PATH=str(p),SIZE=size,MD5=md5,RUN_ID=run.name,PRODUCER=producer,STATE='PUBLISHED')
 side=write(pathlib.Path(str(p)+'.md5'),f'MD5={md5}\nSIZE={size}\nPATH={p}\n')
 key=hashlib.sha256(str(p).encode()).hexdigest()[:32];rp=write(run/'manifests/artifacts'/f'{key}.record',''.join(f'{k}={v}\n' for k,v in value.items()))
 return dict(path=str(p),size=size,md5=md5,record=value,record_path=str(rp),record_sha256=sha(rp),sidecar_sha256=sha(side))
def owner(p,run,state):
 key=p.name+'_'+hashlib.sha256(str(p).encode()).hexdigest()[:32]
 if len(key.encode())>100:key='artifact_'+hashlib.sha256(str(p).encode()).hexdigest()[:32]
 values=dict(RUN_ID=run.name,STAGE='stage5',PATH=str(p),STATE=state,PID='123')
 p=write(scratch/'_ecoda_owners/artifact'/key/'owner',''.join(f'{k}={v}\n' for k,v in values.items()))
 return dict(**proof(p),fields=values,policy='Original state preserved')
datasets=[f'Fixture{i}' for i in range(10)];artifacts={};timings={};workers=[]
for i,ds in enumerate(datasets,1):
 cache=record(write(old/f'artifacts/benchmark/pseudobulks/{ds}_pseudobulk_hvg2000_bl.rds'),old,'stage5_prepare_pseudobulk_hvg2000_bl')
 result=record(write(old/f'artifacts/benchmark/results/{ds}_Pseudobulk_hvg2000_bl.rds'),old,'stage5_pseudobulk')
 artifacts[ds]=dict(cache=cache,result=result,owners=[owner(pathlib.Path(cache['path']),old,'FAIL'),owner(pathlib.Path(result['path']),old,'ACTIVE')])
 timings[ds]=dict(proof=record(write(old/f'logs/execution_times_pseudobulk_{ds}_hvg2000_bl.feather'),old,'stage5_execution_log'))
 workers.append(dict(stdout=proof(write(old/f'logs/worker_{i}.log')),stderr=proof(write(old/f'logs/worker_{i}.err',''))))
success=['4663549','4663561','4663563']+[str(5000000+i) for i in range(18)]
parts=[[i,'COMPLETED','0:0'] for i in success]+[['4663573','FAILED','1:0'],['4663574','FAILED','1:0']]
raw=''.join('|'.join(r)+'\n' for r in parts);rows=[dict(line=i,job_id=r[0],state=r[1],exit_code=r[2]) for i,r in enumerate(parts,1)]
accounting=dict(stdout=raw,rows=rows);audit=dict(accounting=accounting)
ids=success[:3]+['4663573','4663574']
evidence={k:proof(write(old/'manifests'/name,json.dumps(dict(state='FAILED',audit=audit,scheduler_array_ids=ids[:3],scheduler_watchdog_ids=ids[3:])))) for k,name in [('gate_manifest','gate.json'),('inspect','inspect.json')]}
matrix=write(old/'manifests/matrix.tsv',''.join(f'{d}\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n' for d in datasets));evidence['worker_manifest']=proof(matrix)
rawpath=write(old/'manifests/accounting.raw.tsv',raw)
failures=[dict(stdout=proof(write(old/f'logs/{i}.log')),stderr=proof(write(old/f'logs/{i}.err',''))) for i in ['4663573','4663574']]
h5ad=write(t/'never_read.h5ad','MUST_NOT_BE_CONSUMED');h5ad.chmod(0)
inputproof=write(old/'manifests/input.expected.md5','reference only\n');inputrecord=write(old/'manifests/input.record');status=write(old/'status/preflight.status')
science=dict(schema='ecoda_bl_science_receipt_v1',status='NOOP_VALIDATED',science_validated=True,publication_completed=False,original_gate_status='FAILED',original_run_root=str(old),audit_root=str(validation),
 original_evidence=evidence,validation_source=proof(write(validation/'manifests/source.manifest')),validation_runtime=proof(write(validation/'manifests/runtime.identity')),
 production_validator=proof(write(t/'original_validator.R')),original_accounting=dict(query_count=1,worker_array_id='4663563',scheduler_ids=ids,rows=rows,raw_path=str(rawpath),sha256=sha(rawpath),raw_stdout=raw,orchestration_failures=failures),
 input_proofs={ds:dict(record_path=str(inputrecord),record_sha256=sha(inputrecord),reference=[ds,'benchmark_analysis',str(h5ad),'0'*32,'20',str(inputproof)],external_proof_sha256=sha(inputproof),status=proof(status)) for ds in datasets},worker_evidence=workers,artifacts=artifacts,timings=timings)
receipt=write(validation/'artifacts/evidence/bl_science.json',json.dumps(science));record(receipt,validation,'bl_science_validation');owner(receipt,validation,'OK')
review=write(scratch/'_ecoda_inputs/review.json',json.dumps(dict(status='APPROVED_SCIENTIFIC_ARTIFACT_REUSE',receipt_sha256=sha(receipt),datasets=10,rds_artifacts=20,timing_feathers=10,original_gate_status='FAILED',publication_completed=False)))
write(t/'base.json',json.dumps(dict(schema='ecoda_bl_publication_control_v1',release_id='release',nas_root=os.environ['TEST_NAS'],science_receipt=proof(receipt),review=proof(review))))
# Snapshot all original finite evidence and owner bytes/modes; no H5AD reads.
paths=[p for root in (old,validation,scratch/'_ecoda_owners') for p in root.rglob('*') if p.is_file()]
write(t/'original_snapshot.json',json.dumps({str(p):[sha(p),p.stat().st_mode] for p in paths}))
PY
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
invoke() {
  local id="$1" control="${SCRATCH}/_ecoda_inputs/$1.json" md5
  "${REAL_PYTHON}" - "${TMP}/base.json" "${control}" "${id}" <<'PY'
import hashlib,json,pathlib,sys
x=json.loads(pathlib.Path(sys.argv[1]).read_text());x['run_id']=sys.argv[3]
# Bind the currently staged review so explicit review corruption tests reach validation.
p=pathlib.Path(x['review']['path']);x['review']['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
pathlib.Path(sys.argv[2]).write_text(json.dumps(x))
PY
  md5="$(ecoda_md5_file "${control}")"
  ECODA_RUN_ID="${id}" ECODA_RUN_ROOT="${SCRATCH}/_ecoda_runs/${id}" \
    bash "${TREE}/src/utils/bash/publish_bl_science.sh" --control "${control}" --control-md5 "${md5}"
}
expect_fail() { if "$@" > "${TMP}/failure.log" 2>&1; then echo 'Expected publication failure' >&2; exit 1; fi; }
invoke good > "${TMP}/success.log" 2>&1 || { cat "${TMP}/success.log"; exit 1; }
report="${SCRATCH}/_ecoda_runs/good/artifacts/evidence/bl_publication.json"
ecoda_validate_artifact_record "${report}" bl_science_publication good >/dev/null
grep -q '^STATE=NOOP_VALIDATED$' "${SCRATCH}/_ecoda_runs/good/status/terminal"
"${REAL_PYTHON}" - "${report}" <<'PY'
import json,pathlib,sys
x=json.loads(pathlib.Path(sys.argv[1]).read_text());assert x['publication_completed'] and not x['new_model_fit'] and x['original_gate_status']=='FAILED'
assert len(x['mappings'])==30 and len(x['original_accounting']['rows'])==23
assert len(list(pathlib.Path(x['destination_root']).rglob('*.*')))==61
PY
before="$(ecoda_md5_file "${report}")"; expect_fail invoke good
[[ "$(ecoda_md5_file "${report}")" == "${before}" ]]
control="${SCRATCH}/_ecoda_runs/good/manifests/control.json"
validate_release() {
  "${REAL_PYTHON}" "${TREE}/src/utils/py/bl_science_publication.py" validate --control "${control}" \
    --control-md5 "$(ecoda_md5_file "${control}")" --run-root "${SCRATCH}/_ecoda_runs/good" --nas-root "${NAS}" \
    --inventory "${SCRATCH}/_ecoda_runs/good/manifests/publication_inventory.tsv" --report "${report}"
}
validate_release > "${TMP}/validation.log"
grep -q '^BL_PUBLICATION_CONTRACT_OK$' "${TMP}/validation.log"
nas_payload="${NAS}/releases/release/good/benchmark/results/Fixture0_Pseudobulk_hvg2000_bl.rds"
# This sandbox strips setuid/setgid; sticky bits remain observable on files.
chmod 1666 "${nas_payload}"; expect_fail validate_release; chmod 666 "${nas_payload}"
chmod 1777 "${nas_payload%/*}"; expect_fail validate_release; chmod 777 "${nas_payload%/*}"
nas_checksums="${NAS}/releases/release/good/benchmark/checksums.md5"
cp "${nas_checksums}" "${TMP}/checksums.saved"; printf tamper > "${nas_checksums}"
expect_fail validate_release
cp "${TMP}/checksums.saved" "${nas_checksums}"
expect_fail "${REAL_PYTHON}" "${TREE}/src/utils/py/bl_science_publication.py" copy \
  --source "${TMP}/never_read.h5ad" --destination "${SCRATCH}/_ecoda_runs/good/forbidden.h5ad" \
  --boundary "${SCRATCH}/_ecoda_runs/good" --md5 00000000000000000000000000000000 --size 20
[[ ! -e "${SCRATCH}/_ecoda_runs/good/forbidden.h5ad" ]]
TEST_MISSING_MOUNT=1 expect_fail invoke missing_mount
[[ ! -e "${NAS}/releases/release/missing_mount" ]]
TEST_MODE_FAILURE=1 expect_fail invoke bad_mode
grep -q '^STATE=FAIL$' "${SCRATCH}/_ecoda_runs/bad_mode/status/terminal"
TEST_COPY_FAILURE=1 expect_fail invoke failed_copy
grep -q '^STATE=FAIL$' "${SCRATCH}/_ecoda_runs/failed_copy/status/terminal"
mkdir -p "${NAS}/releases/release/conflict/benchmark"; printf foreign > "${NAS}/releases/release/conflict/benchmark/foreign"
expect_fail invoke conflict
[[ "$(cat "${NAS}/releases/release/conflict/benchmark/foreign")" == foreign ]]
cp "${SCRATCH}/_ecoda_inputs/review.json" "${TMP}/review.saved"
printf '{"status":"REJECTED"}\n' > "${SCRATCH}/_ecoda_inputs/review.json"
expect_fail invoke bad_review
cp "${TMP}/review.saved" "${SCRATCH}/_ecoda_inputs/review.json"
old="${SCRATCH}/_ecoda_runs/stage5_benchmark_bl10_20261109_v2/artifacts/benchmark/results/Fixture0_Pseudobulk_hvg2000_bl.rds"
cp "${old}" "${TMP}/payload.saved"; printf changed > "${old}"
expect_fail invoke bad_payload
cp "${TMP}/payload.saved" "${old}"
owner_file="$(ecoda_artifact_owner_dir "${old}")/owner"
cp "${owner_file}" "${TMP}/owner.saved"; printf 'STATE=OK\n' > "${owner_file}"
expect_fail invoke bad_owner
cp "${TMP}/owner.saved" "${owner_file}"
"${REAL_PYTHON}" - "${TMP}/original_snapshot.json" <<'PY'
import hashlib,json,pathlib,sys
for name,(sha,mode) in json.loads(pathlib.Path(sys.argv[1]).read_text()).items():
 p=pathlib.Path(name);assert hashlib.sha256(p.read_bytes()).hexdigest()==sha and p.stat().st_mode==mode,name
PY
"${REAL_PYTHON}" - "${TREE}/src/utils/py/bl_science_publication.py" "${TMP}" <<'PY'
import importlib.util,json,pathlib,sys
from unittest.mock import patch
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('bl_publication',sys.argv[1]);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
root=pathlib.Path(sys.argv[2]);target=root/'atomic.json';module.write_receipt(target,{'complete':True})
assert json.loads(target.read_text())=={'complete':True}
try:module.write_receipt(target,{'replace':True})
except FileExistsError:pass
else:raise AssertionError('existing receipt replaced')
assert json.loads(target.read_text())=={'complete':True}
failed=root/'failed-atomic.json'
with patch.object(module.os,'link',side_effect=OSError('injected install failure')):
 try:module.write_receipt(failed,{'partial':True})
 except OSError:pass
 else:raise AssertionError('injected failure ignored')
assert not failed.exists() and not list(root.glob('.bl-receipt-*'))
PY
[[ ! -e "${TMP}/forbidden.called" ]]
echo 'PASS: finite reviewed BL publication, exact30 payloads/60checksums, new owners/runtime/source, original FAILED evidence unchanged, create-only and mount/mode/copy/review/payload/owner failures'
