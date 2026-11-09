#!/usr/bin/env python3
"""Real audit entrypoint, tiny sparse H5AD, and preserved bounded OOM evidence."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py import stage3_reuse_report as reuse


def setup(base, source, scratch, runtime, source_manifest):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from src.utils.py.cohort_selection import cohort_configuration, resolve_cohort_selection, cohort_record, input_file_identities
    from src.utils.py.datasets_io import read_datasets_json
    run = scratch / '_ecoda_runs/reuse_fixture'; manifests = run / 'manifests'; manifests.mkdir(parents=True, exist_ok=True)
    pairs = [(d, v) for d in ('Alzheimer','Breast_cancer','Covid19_PBMC','Joanito','Kidney_KPMP_full','Lupus_PBMC','Lung','Stephenson') for v in ('batch_effect_uncorrected','batch_effect_corrected')]
    pairs += [('Stephenson','benchmark_analysis')]
    (manifests/'selection.tsv').write_text(''.join('\t'.join(r)+'\n' for r in pairs))
    (base/'subset.tsv').write_text('Stephenson\tbenchmark_analysis\n')
    for i in (1, 2): (manifests/f'selection.retry_{i}.tsv').write_text(''.join('\t'.join(r)+'\n' for r in pairs[:2]))
    config = json.loads((source/'datasets.json').read_text()); refs=[]
    for d,v in pairs:
        name=config[d]['views'][v]['input_file_name']; assert isinstance(name,str)
        p=base/'raw'/d/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'tiny source identity fixture')
        md5=hashlib.md5(p.read_bytes()).hexdigest();proof=Path(str(p)+'.md5');proof.write_text(f'MD5={md5}\nSIZE={p.stat().st_size}\nPATH={p}\n')
        refs.append(f'{d}\t{v}\t{p}\t{md5}\t{p.stat().st_size}\t{proof}\n')
    (base/'refs.tsv').write_text(''.join(refs))
    entry=read_datasets_json(source/'datasets.json')['Stephenson']; view='benchmark_analysis'
    obs=pd.DataFrame({'Sample':np.repeat(['good_a','good_b'],500),'Status':np.repeat(['Healthy','Covid'],500),'Site':'Ncl'},index=[f'cell_{i}' for i in range(1000)])
    cfg=cohort_configuration(dataset='Stephenson',view=view,sample_column=entry['sample_col'],label_column=entry['label_col'],batch_keys=entry['batch_col'],subset_vars=entry['views'][view]['subset_vars'],batch_metadata_policy=entry['batch_metadata_policy'],evaluation_labels=entry['evaluation_labels'])
    mask,_,audit=resolve_cohort_selection(obs,cfg); assert mask.all()
    counts=sp.csr_matrix((np.ones(1000,dtype=np.float32),(np.arange(1000),np.arange(1000))),shape=(1000,3000))
    data=ad.AnnData(X=counts.copy(),obs=obs,var=pd.DataFrame({'hvg_rank':np.arange(3000,dtype=float)},index=[f'gene_{i}' for i in range(3000)]));data.layers['counts']=counts
    for key in ('X_pca_benchmark_analysis_hvg1000','X_pca_benchmark_analysis_hvg2000','X_pca_benchmark_analysis_hvg3000','X_pca_harmony_benchmark_analysis_hvg2000'):data.obsm[key]=np.zeros((1000,2),dtype=np.float32)
    data.uns['cohort_selection']=cohort_record(audit,input_file_identities(entry['views'][view]['input_file'],base/'raw/Stephenson'),obs)
    output=run/'artifacts/datasets/Stephenson/output'/config['Stephenson']['views'][view]['output_file_name'];output.parent.mkdir(parents=True);data.write_h5ad(output)
    (base/'output.txt').write_text(str(output))
    ids=['101','102','103','104','105']
    (manifests/'scheduler_ids.tsv').write_text(''.join(f'{kind}\t{job}\n' for kind,job in zip(['PREFLIGHT','ARRAY','WATCHDOG','STATUS','STATUS'],ids)))
    (manifests/'h5ad_obs_audit.tsv').write_text('Covid19_PBMC\tbatch_effect_uncorrected\tfixture\nCovid19_PBMC\tbatch_effect_corrected\tfixture\n')
    logs=scratch/'_ecoda_run_logs'/run.name;logs.mkdir(parents=True)
    witnesses=[('3_scrnaseq_preprocessing_102_1.err','201'),('3_scrnaseq_preprocessing_102_2.err','202'),('3_scrnaseq_preprocessing_retry1_104_1.err','204'),('3_scrnaseq_preprocessing_retry1_104_2.err','104')]
    for name,raw in witnesses:(logs/name).write_text(f'slurmstepd: error: Detected 1 oom_kill event in StepId={raw}.batch\n')
    (run/'status').mkdir(exist_ok=True)
    (run/'status/watchdog').write_text('STATE=OK\nRUN_ID=reuse_fixture\nRETRY_INDEX=2\nARRAY_JOB_ID=105\n'+''.join(f'SCHEDULER_ID={i}\n' for i in ids[1:]))
    (run/'status/terminal').write_text('STATE=OK\nRUN_ID=reuse_fixture\nREASON=preprocessing validation and selected sync completed\n')
    cmd=shlex.join(['immutable_executor','--run-id',run.name,'--scratch-root',str(scratch),'--script','src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh','--source-manifest',str(source_manifest),'--selection-file',str(manifests/'selection.tsv'),'--logs-root',str(logs),'--mem','128G','--max-mem','500G'])
    raw_ids=['101','301','102','103','201','202']+[str(i) for i in range(302,316)]+['104','204','105','316'];assert len(raw_ids)==24
    raw=[[i,'OUT_OF_MEMORY' if i in ('201','202','204','104') else 'COMPLETED','0:125' if i in ('201','202','204','104') else '0:0'] for i in raw_ids]
    accounting={'query_count':1,'actual_exit':0,'attempted':True,'malformed_rows':[],'missing_scheduler_ids':[],'scheduler_ids':ids,'command':'sacct -n -P -X -j '+','.join(ids)+' --format=JobIDRaw,State,ExitCode','stdout':''.join('|'.join(r)+'\n' for r in raw),'rows':[dict(job_id=r[0],state=r[1],exit_code=r[2]) for r in raw]}
    audit={'accounting':accounting,'commands':[],'artifact_contracts':[],'immutable_fingerprints':[]}
    digest=hashlib.sha256(cmd.encode()).hexdigest();terminal={'state':'COMPLETED','exit_code':0,'gate_id':'failed_fixture','exact_command':cmd,'command_digest':digest,'event_generation':'generation'}
    gate={'state':'FAILED','gate_id':'failed_fixture','audit_completed':True,'audit_state':'COMPLETED','audit':audit,'terminal_status':terminal,'exact_command':cmd,'command_digest':digest,'event_generation':'generation','scheduler_array_ids':[ids[i] for i in (0,1,3,4)],'scheduler_watchdog_ids':['103'],'discrepancies':[{'kind':'terminal_audit_failed'}]}
    (base/'gate.json').write_text(json.dumps(gate));(base/'inspect.json').write_text(json.dumps({'state':'FAILED','gate_id':'failed_fixture','audit':audit}))


def verify(base, source, scratch, source_manifest):
    run=scratch/'_ecoda_runs/reuse_fixture'
    args=['--run-root',str(run),'--stage','stage3','--selection',str(run/'manifests/selection.tsv'),'--source-manifest',str(source_manifest),'--runtime-identity',str(run/'manifests/runtime.identity'),'--stage3-reuse-report',str(base/'reuse.json'),'--reuse-selection',str(base/'subset.tsv'),'--terminal-inspect',str(base/'inspect.json'),'--failed-gate-manifest',str(base/'gate.json')]
    audit=ROOT/'src/utils/bash/ecoda_run_audit.sh'
    def invoke(options=args):return subprocess.run(['bash',str(audit),*options],capture_output=True,text=True)
    protected=list(run.rglob('*'))+list((base/'raw').rglob('*'))+list((base/'nas').rglob('*'))+[base/'gate.json',base/'inspect.json']
    before={str(p):(hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mode) for p in protected if p.is_file()}
    selected_raw=base/'raw/Stephenson'/json.loads((source/'datasets.json').read_text())['Stephenson']['views']['benchmark_analysis']['input_file_name']
    raw_mode=selected_raw.stat().st_mode
    selected_raw.chmod(0)  # Actual child process must never read raw contents.
    before.pop(str(selected_raw))
    success=invoke();assert success.returncode==0,(success.stdout,success.stderr)
    assert selected_raw.stat().st_mode & 0o777 == 0
    selected_raw.chmod(raw_mode)
    assert before=={p:(hashlib.sha256(Path(p).read_bytes()).hexdigest(),Path(p).stat().st_mode) for p in before}
    ordinary=invoke(args[:args.index('--stage3-reuse-report')])
    assert ordinary.returncode!=0 and 'scheduler accounting query failed' in ordinary.stderr,ordinary.stderr
    report=json.loads((base/'reuse.json').read_text());assert report['status']=='NOOP_VALIDATED' and report['failed_gate_released'] is False and len(report['history']['accounting_rows'])==24 and len(report['artifacts'])==1 and report['artifacts'][0]['n_samples']==2
    assert len(report['history']['oom_witnesses'])==4
    receipt=reuse.checksum(base/'reuse.json');assert receipt['md5']==hashlib.md5((base/'reuse.json').read_bytes()).hexdigest()
    assert (base/'reuse.json').stat().st_mode & 0o777 == (base/'reuse.json.md5').stat().st_mode & 0o777 == 0o444
    assert 'no fresh raw content hash' in report['artifacts'][0]['input']['verification']
    # No unselected output exists, so success proves subset-only expensive reads.
    assert len(list((run/'artifacts/datasets').glob('*/output/*.h5ad')))==1
    assert invoke().returncode!=0  # Create-only report.
    original_gate=(base/'gate.json').read_bytes();original_inspect=(base/'inspect.json').read_bytes()
    namespace=argparse.Namespace(run_root=str(run),selection=str(run/'manifests/selection.tsv'),reuse_selection=str(base/'subset.tsv'),source_manifest=str(source_manifest),runtime_identity=str(run/'manifests/runtime.identity'),terminal_inspect=str(base/'inspect.json'),failed_gate_manifest=str(base/'gate.json'))
    def reject_change(change):
        g=json.loads(original_gate);change(g);(base/'gate.json').write_text(json.dumps(g));(base/'inspect.json').write_text(json.dumps({'state':'FAILED','gate_id':'failed_fixture','audit':g['audit']}))
        try:reuse.history(namespace)
        except ValueError:pass
        else:raise AssertionError('invalid lineage accepted')
        (base/'gate.json').write_bytes(original_gate);(base/'inspect.json').write_bytes(original_inspect)
    reject_change(lambda g:g['audit']['accounting'].update(stdout=g['audit']['accounting']['stdout']+'999|COMPLETED|0:0\n'))
    def non_oom(g):
        a=g['audit']['accounting'];a['stdout']=a['stdout'].replace('201|OUT_OF_MEMORY','201|FAILED')
        for row in a['rows']:
            if row['job_id']=='201':row['state']='FAILED'
    reject_change(non_oom)
    witness=scratch/'_ecoda_run_logs/reuse_fixture/3_scrnaseq_preprocessing_retry1_104_1.err';saved=witness.read_bytes();witness.write_text('StepId=999.batch oom_kill')
    try:reuse.history(namespace)
    except ValueError:pass
    else:raise AssertionError('foreign OOM witness accepted')
    witness.write_bytes(saved)
    before=(base/'subset.tsv').read_bytes();(base/'subset.tsv').write_text('Unknown\tbenchmark_analysis\n')
    try:reuse.history(namespace)
    except ValueError:pass
    else:raise AssertionError('foreign subset accepted')
    (base/'subset.tsv').write_bytes(before)
    nas=Path(report['artifacts'][0]['nas']['path']);old=nas.read_bytes();nas.write_bytes(old+b'X')
    fresh=args.copy();fresh[fresh.index('--stage3-reuse-report')+1]=str(base/'fault.json')
    failure=invoke(fresh);assert failure.returncode!=0 and 'checksum mismatch' in failure.stderr,failure.stderr
    nas.write_bytes(old);os.chmod(nas,0o644)
    failure=invoke(fresh);assert failure.returncode!=0 and 'mode differs' in failure.stderr,failure.stderr
    os.chmod(nas,0o666);assert not (base/'fault.json').exists()
    assert (base/'gate.json').read_bytes()==original_gate and (base/'inspect.json').read_bytes()==original_inspect
    print('PASS: actual Stage3 subset audit, full witnessed OOM lineage, selected-only scientific/NAS reads, create-only receipt and immutable failed gate')


if __name__=='__main__':
    import argparse
    if len(sys.argv)>1:
        mode,*paths=sys.argv[1:];paths=list(map(Path,paths))
        if mode=='--fixture':setup(*paths)
        else:verify(*paths)
    else:
        prefix=(ROOT/'tests/test_benchmark_matrix_submitter.sh').read_text().split('if [[ "${ECODA_MATRIX_REQUEST_ONLY:-0}" == 1 ]]; then',1)[0]
        prefix=prefix.replace('touch "${SNAPSHOT_ROOT}/COMPLETE"', "printf 'COMPLETE\\n' > \"${SNAPSHOT_ROOT}/COMPLETE\"")
        prefix=prefix.replace('ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"','ROOT='+shlex.quote(str(ROOT)))
        body=r'''
export PROJECT_ROOT="${SOURCE_TREE}" DATASETS_JSON_FILE="${SOURCE_TREE}/datasets.json" NAS_TARGET_DIR="${TMP_DIR}/nas"
mkdir -p "${NAS_TARGET_DIR}"
export PYTHON_BIN="${TEST_REAL_PYTHON}" ECODA_GATE_STAGE=stage3
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/utils/bash/ecoda_runtime.sh"
source "${ROOT}/src/utils/bash/ecoda_derivative_artifacts.sh"
ecoda_init_run stage3 reuse_fixture >/dev/null
cp "${SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
ecoda_runtime_validate_submission apptainer
"${TEST_REAL_PYTHON}" "${ROOT}/tests/test_stage3_reuse_report.py" --fixture "${TMP_DIR}" "${SOURCE_TREE}" "${HPC_ROOT}" "${RUNTIME_MANIFEST}" "${SOURCE_MANIFEST}"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_derivative_bind stage3 test-release "${ECODA_RUN_ROOT}/manifests/selection.tsv" "${TMP_DIR}/refs.tsv"
output="$(cat "${TMP_DIR}/output.txt")"
nas="${NAS_TARGET_DIR}/releases/test-release/reuse_fixture/datasets/Stephenson/output/${output##*/}"
mkdir -p "${nas%/*}"
cp "${output}" "${nas}"
for path in "${output}" "${nas}"; do
 ecoda_write_checksum "${path}"
 ecoda_artifact_owner_acquire "${path}" stage3 reuse_fixture 0 0 0 >/dev/null
 ecoda_artifact_owner_set_state "${path}" OK fixture >/dev/null
done
ecoda_write_artifact_record "${output}" stage3 reuse_fixture >/dev/null
chmod 666 "${nas}" "${nas}.md5"
"${TEST_REAL_PYTHON}" -c 'import pathlib,sys; p=pathlib.Path(sys.argv[1]); [d.chmod(0o777) for d in [p,*p.rglob("*")] if d.is_dir()]' "${NAS_TARGET_DIR}/releases"
# Any accidental scheduler query/submission fails deterministically.
printf '#!/bin/bash\necho FORBIDDEN-SCHEDULER >&2\nexit 88\n' > "${TMP_DIR}/bin/sacct"
chmod +x "${TMP_DIR}/bin/sacct"
export PATH="${TMP_DIR}/bin:${PATH}"
"${TEST_REAL_PYTHON}" "${ROOT}/tests/test_stage3_reuse_report.py" --verify "${TMP_DIR}" "${SOURCE_TREE}" "${HPC_ROOT}" "${SOURCE_MANIFEST}"
'''
        with tempfile.TemporaryDirectory(prefix='ecoda-reuse-test-') as temporary:
            p=Path(temporary)/'run.sh';p.write_text(prefix+body)
            result=subprocess.run(['bash',str(p)],cwd=ROOT,text=True,capture_output=True)
            print(result.stdout,end='');print(result.stderr,end='',file=sys.stderr)
            raise SystemExit(result.returncode)
