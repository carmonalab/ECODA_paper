"""Fixed four-fit release: real transfer/owners, immutable fit evidence, no jobs."""
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import numpy as np
import pandas as pd

from test_main_v2_science_publication import build, write, proof, sha
from bl_science_publication import UNC4_SCHEMA, UNC4_SOURCE, UNC4_RUN, UNC4_DATASETS, UNC4_SAMPLES, owner_file, verified_scope
ROOT=Path(__file__).resolve().parents[1]
CONTEXT='batch_effect/uncorrected_final'


def fixture(base):
    x=build(base,CONTEXT); scratch=x['scratch']; old=scratch/'_ecoda_runs'/UNC4_RUN; parent=old.parent/'stage5_failed16_20261109_v1'
    manifest=old/'manifests'; manifest.mkdir(parents=True)
    source=x['source']; config=x['config']
    sm=write(manifest/'source.manifest',f'SOURCE_ROOT={source}\nSOURCE_COMMIT={UNC4_SOURCE}\nAUX_ROOT={source}/aux\nDATASETS_SHA256={sha(config)}\n')
    original_manifest=write(source/'unc4.source.manifest',sm.read_text())
    runtime=write(manifest/'runtime.identity',x['runtime'].read_text())
    selection=write(manifest/'selection.tsv',''.join(f'{ds}\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n' for ds in UNC4_DATASETS))
    matrix=write(manifest/'matrix_batch_effect_uncorrected_mrvi.tsv',''.join(f'{ds}\tbatch_effect_uncorrected\tmrvi\n' for ds in UNC4_DATASETS))
    dispatch=write(manifest/'dispatch_selection.tsv',matrix.read_text())
    def sealed(path,producer):
        md5=hashlib.md5(path.read_bytes()).hexdigest();size=str(path.stat().st_size)
        side=write(Path(str(path)+'.md5'),f'MD5={md5}\nSIZE={size}\nPATH={path}\n')
        rec=write(manifest/'artifacts'/(hashlib.sha256(str(path).encode()).hexdigest()[:32]+'.record'),f'PATH={path}\nSIZE={size}\nMD5={md5}\nRUN_ID={old.name}\nPRODUCER={producer}\nSTATE=PUBLISHED\n')
        return md5,size,side,rec
    sealed(matrix,'stage5_matrix');sealed(selection,'stage5_selection')
    child=write(manifest/'scheduler_ids.tsv','PREFLIGHT\t4666178\nARRAY\t4666187\nWATCHDOG\t4666191\nAGGREGATE_GATE\t4666192\n');sealed(child,'stage5_schedulers')
    scheduler_rows=[['benchmark','stage5_failed16_benchmark_20261109_v1','PREFLIGHT','4666176'],['benchmark','stage5_failed16_benchmark_20261109_v1','PREFLIGHT','4666180'],['benchmark','stage5_failed16_benchmark_20261109_v1','ARRAY','4666183'],['benchmark','stage5_failed16_benchmark_20261109_v1','WATCHDOG','4666184'],['benchmark','stage5_failed16_benchmark_20261109_v1','AGGREGATE_GATE','4666185']]
    scheduler_rows += [['uncorrected',old.name,*row.split('\t')] for row in child.read_text().splitlines()]
    scheduler_rows += [['corrected','stage5_failed16_corrected_20261109_v1','PREFLIGHT','4666196']]
    parent_ids=write(parent/'manifests/campaign_scheduler_ids.tsv',''.join('\t'.join(row)+'\n' for row in scheduler_rows))
    ids=[row[3] for row in scheduler_rows]
    raw_ids=['4666176','4666179','4666181','4666182','4666178','4666180','4666183','4666184','4666185','4666188','4666189','4666190','4666187','4666191','4666192','4666197','4666199','4666200','4666202','4666203','4666204','4666205','4666196']
    raw=''.join(job+('|FAILED|0:53\n' if job in {'4666183','4666196','4666204','4666205'} else '|FAILED|1:0\n' if job in {'4666184','4666185'} else '|COMPLETED|0:0\n') for job in raw_ids)
    accounting=dict(query_count=1,attempted=True,actual_exit=0,scheduler_ids=ids,command='sacct -n -P -X -j '+','.join(ids)+' --format=JobIDRaw,State,ExitCode',stdout=raw,rows=[dict(line=i,job_id=r[0],state=r[1],exit_code=r[2]) for i,r in enumerate((line.split('|') for line in raw.splitlines()),1)])
    command=shlex.join(['executor','--run-id',parent.name,'--source-root',str(source),'--source-manifest',str(original_manifest),'--scratch-root',str(scratch),'--runtime-image',str(x['image']),'--runtime-manifest',str(x['image_manifest']),'--script','src/5_run_benchmark_methods/stage5_campaign.sh'])
    gate=dict(gate_id='ecoda_stage5_failed16_20261109_v1',state='FAILED',profile='fixture',profile_digest=sha(x['profile']),exact_command=command,command_digest=hashlib.sha256(command.encode()).hexdigest(),event_generation='real-generation',scheduler_array_ids=ids,scheduler_watchdog_ids=[],audit_completed=True,audit_state='COMPLETED',audit=dict(passed=False,accounting=accounting))
    gate['terminal_status']=dict(schema_version='1.0',state='FAILED',exit_code=1,started_at='2026-11-09T16:41:50Z',finished_at='2026-11-09T17:38:05Z',**{key:gate[key] for key in ('gate_id','exact_command','command_digest','event_generation')})
    gate_path=write(base/'evidence/unc4.gate.json',json.dumps(gate));inspect=write(base/'evidence/unc4.inspect.json',json.dumps(dict(operation='inspect',state='FAILED',gate_id=gate['gate_id'],audit=gate['audit'])))
    metadata=write(old/'metadata','STATE=ACTIVE\nTARGET_METHODS=mrvi\n')
    terminal=write(old/'status/terminal','STATE=FAIL\nREASON=Stage 5 final matrix artifact validation failed\n')
    aggregate=write(old/'status/aggregate','STATE=OK\nWATCHDOG_LABELS=batch_effect_uncorrected__mrvi\nSCHEDULER_ID=4666178\nSCHEDULER_ID=4666187\nSCHEDULER_ID=4666191\nSCHEDULER_ID=4666192\n')
    wd=write(old/'status/watchdogs/batch_effect_uncorrected__mrvi.status','STATE=OK\nLABEL=batch_effect_uncorrected__mrvi\nARRAY_JOB_ID=4666187\nRETRY_INDEX=0\nSCHEDULER_ID=4666187\n')
    refs=[];identities=[];distances=[];companions=[];distance_owners=[];row_owners=[];owner_rows=[];workers=[]
    registry=json.loads(config.read_text())
    for i,(ds,n) in enumerate(zip(UNC4_DATASETS,UNC4_SAMPLES),1):
        source_h5ad=write(base/'inputs'/ds/registry[ds]['views']['batch_effect_uncorrected']['output_file_name'],'DO NOT READ EXPRESSION\n')
        imd5=hashlib.md5(source_h5ad.read_bytes()).hexdigest();isize=source_h5ad.stat().st_size;ip=write(Path(str(source_h5ad)+'.md5'),f'MD5={imd5}\nSIZE={isize}\nPATH={source_h5ad}\n');source_h5ad.chmod(0)
        refs.append(f'{ds}\tbatch_effect_uncorrected\t{source_h5ad}\t{imd5}\t{isize}\t{ip}\n')
        sample_ids=[f's{k}' for k in range(n)];identities.append(dict(dataset=ds,view='batch_effect_uncorrected',path=str(source_h5ad),md5=imd5,size=isize,sample_ids=sample_ids))
        path=old/'artifacts'/CONTEXT/'embeddings'/f'{ds}_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather';path.parent.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(np.ones((n,n))-np.eye(n),index=sample_ids,columns=sample_ids).to_feather(path)
        md5,size,_,rec=sealed(path,'stage5_mrvi')
        values=dict(RUN_ID=old.name,STAGE='stage5',PATH=str(path),STATE='OK');of=write(owner_file(path,old),''.join(f'{k}={v}\n' for k,v in values.items()));distance_owners.append(dict(**proof(of),fields=values))
        key=f'derivative/{old.name}/uncorrected/{ds}/batch_effect_uncorrected/mrvi';row_values=dict(RUN_ID=old.name,STAGE='stage5',KEY=key,STATE='OK');row_owner=write(scratch/'_ecoda_owners/stage5'/re.sub(r'[/,:\t |]','_',key)/'owner',''.join(f'{k}={v}\n' for k,v in row_values.items()));row_owners.append(dict(**proof(row_owner),fields=row_values));owner_rows.append(f'{key}\t{row_owner.parent}\n')
        preflight=write(old/f'status/h5ad_preflight/{ds}__batch_effect_uncorrected.status',f'STATE=OK\nRUN_ID={old.name}\nDATASET={ds}\nVIEW=batch_effect_uncorrected\nTASK_ID={i}\n')
        distances.append(dict(dataset=ds,path=str(path),md5=md5,size=int(size),n_samples=n,frame_shape=[n,n],input_md5=imd5,input_preflight=proof(preflight),producer_record=dict(PATH=str(path),SIZE=size,MD5=md5,RUN_ID=old.name,PRODUCER='stage5_mrvi',STATE='PUBLISHED')))
        payload=dict(schema_version=1,artifact_path=str(path),artifact_md5=md5,dataset=ds,method='MrVI_hvg2000',time_secs=float(i),mem_GB=1.0)
        rt=write(Path(str(path)+'.runtime.json'),json.dumps(payload));rm,rs,rside,rrec=sealed(rt,'stage5_mrvi')
        timing=old/'logs'/f'execution_times_batch_effect_uncorrected_final_mrvi_{ds}.feather';timing.parent.mkdir(exist_ok=True)
        frame=pd.DataFrame([dict(dataset=ds,method='MrVI_hvg2000',time_secs=float(i),mem_GB=1.0)]);frame.to_feather(timing);tm,ts,tside,trec=sealed(timing,'stage5_execution_log')
        def observed(p,m,s,side,record):return dict(path=str(p),exists=True,md5=m,size=int(s),sidecar=proof(side),producer_record=dict(**proof(record),fields=dict(line.split('=',1) for line in record.read_text().splitlines())),owner=dict(path=str(owner_file(p,old)),exists=False))
        runtime_row=observed(rt,rm,rs,rside,rrec);runtime_row['payload']=payload
        timing_row=observed(timing,tm,ts,tside,trec);timing_row.update(columns=list(frame.columns),rows=json.loads(frame.to_json(orient='records')))
        companions.append(dict(dataset=ds,runtime=runtime_row,timing=timing_row))
        workers.append(dict(dataset=ds,task_index=i,stdout=proof(write(old/f'logs/5_matrix_batch_effect_uncorrected__mrvi_4666187_{i}.log','fitted\n')),stderr=proof(write(old/f'logs/5_matrix_batch_effect_uncorrected__mrvi_4666187_{i}.err',''))))
    owners=write(manifest/'owners.tsv',''.join(owner_rows));inputs=write(manifest/'input_artifacts.tsv',''.join(refs));identity=write(manifest/'source_identity.json',json.dumps(dict(entries=identities)))
    binding=write(manifest/'artifact_bindings.json',json.dumps(dict(schema='ecoda_derivative_artifacts_v1',stage='stage5',run_id=old.name,run_root=str(old),selection=str(selection),input_artifacts=str(inputs),identities={str(p):sha(p) for p in (config,sm,runtime,selection,inputs)})))
    patched=base/'patched/tree';patched_validator=write(patched/'src/5_run_benchmark_methods/matrix_artifact_validator.py',(ROOT/'src/5_run_benchmark_methods/matrix_artifact_validator.py').read_text())
    patched_manifest=write(base/'patched/identity/source.manifest',f'SOURCE_ROOT={patched}\nSOURCE_COMMIT=7a6c513f7f3d3bc4b1f0710029e12aef2d26aedd\n')
    control_root=scratch/'_ecoda_inputs/unc4'
    validation=write(control_root/'validation.json',json.dumps(dict(schema='ecoda_failed16_unc4_readonly_validation_v1',status='VALIDATED_CURRENT_FOUR_MRVI_ARTIFACTS',run_id=old.name,original_source_commit=UNC4_SOURCE,no_source_writes=True,no_accounting_query=True,no_model_fit=True,no_h5ad_payload_read=True,validation_source=dict(**proof(patched_manifest),commit='7a6c513f7f3d3bc4b1f0710029e12aef2d26aedd'),validator=proof(patched_validator),config=proof(config),binding=proof(binding),unchanged=[proof(p) for p in (binding,identity,selection,metadata,terminal,aggregate)],payloads=distances)))
    companion_report=write(control_root/'companions.json',json.dumps(dict(schema='ecoda_failed16_unc4_companion_audit_v1',run_id=old.name,no_source_writes=True,no_accounting_query=True,no_scheduler_query=True,no_model_fit=True,no_h5ad_payload_read=True,rows=companions)))
    evidence=dict(gate_manifest=gate_path,inspect=inspect,profile=x['profile'],source_manifest=sm,runtime_identity=runtime,binding=binding,source_identity=identity,input_artifacts=inputs,selection=selection,dispatch=dispatch,child_scheduler_ids=child,metadata=metadata,owners_manifest=owners,terminal=terminal,aggregate=aggregate,watchdog=wd,matrix=matrix,parent_scheduler_ids=parent_ids)
    science=write(control_root/'science.json',json.dumps(dict(schema='ecoda_failed16_unc4_science_receipt_v1',status='NOOP_VALIDATED',science_validated=True,publication_completed=False,original_gate_status='FAILED',source_commit=UNC4_SOURCE,original_run_root=str(old),worker_exit_not_established=True,validation=proof(validation),companions=proof(companion_report),original_evidence={key:proof(value) for key,value in evidence.items()},original_accounting=accounting,worker_evidence=workers,distance_owners=distance_owners,row_owners=row_owners)))
    review=write(control_root/'review.json',json.dumps(dict(status='APPROVED_SCIENTIFIC_ARTIFACT_REUSE',receipt_sha256=sha(science),validation_sha256=sha(validation),companions_sha256=sha(companion_report),datasets=4,distance_feathers=4,runtime_jsons=4,timing_feathers=4,payloads=12,original_gate_status='FAILED',publication_completed=False,worker_exit_not_established=True)))
    def invoke(name):
        control=write(control_root/(name+'.json'),json.dumps(dict(schema=UNC4_SCHEMA,run_id=name,release_id='unc4-release',nas_root=str(x['nas']),science_receipt=proof(science),review=proof(review))))
        return subprocess.run(['bash',str(x['tree']/'src/utils/bash/publish_bl_science.sh'),'--control',str(control),'--control-md5',hashlib.md5(control.read_bytes()).hexdigest()],env=dict(x['env'],ECODA_RUN_ID=name,ECODA_RUN_ROOT=str(scratch/'_ecoda_runs'/name)),capture_output=True,text=True)
    return locals()


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-unc4-publication-') as temporary:
        x=fixture(Path(temporary).resolve());old=x['old'];before={p:(sha(p),p.stat().st_mode) for root in (old,x['scratch']/'_ecoda_owners') for p in root.rglob('*') if p.is_file()}
        result=x['invoke']('published');assert result.returncode==0,result.stdout+result.stderr
        run_root=x['scratch']/'_ecoda_runs/published';report_path=run_root/'artifacts/evidence/failed16_unc4_science_publication.json';report=json.loads(report_path.read_text())
        assert len(report['mappings'])==12 and report['original_gate_status']=='FAILED' and not report['new_model_fit']
        assert report['sample_metadata_published'] is False and report['runtime_companion_policy']=='BYTE_IDENTICAL_ORIGINAL_FIT_EVIDENCE_NOT_REBASED'
        assert len((Path(report['destination_root'])/'checksums.md5').read_text().splitlines())==24
        for row in report['mappings']:
            if row['kind']=='runtime':assert json.loads(Path(row['scratch']).read_text())['artifact_path'].startswith(str(old))
        for path,value in before.items():assert (sha(path),path.stat().st_mode)==value,str(path)
        control=run_root/'manifests/control.json'
        command=[sys.executable,str(x['x']['tree']/'src/utils/py/bl_science_publication.py'),'validate','--control',str(control),'--control-md5',hashlib.md5(control.read_bytes()).hexdigest(),'--run-root',str(run_root),'--nas-root',str(x['x']['nas']),'--inventory',str(run_root/'manifests/publication_inventory.tsv'),'--report',str(report_path)]
        result=subprocess.run(command,env=x['x']['env'],capture_output=True,text=True);assert result.returncode==0,result.stdout+result.stderr
        def check():return verified_scope(SimpleNamespace(control=str(control),control_md5=hashlib.md5(control.read_bytes()).hexdigest(),run_root=str(run_root),nas_root=str(x['x']['nas'])))
        check()
        # Mutations retain real original identity/proof guards; never recompute or repair.
        for path,transform in ((x['matrix'],lambda t:t.replace('Lung','Foreign')),(x['wd'],lambda t:t.replace('STATE=OK','STATE=FAIL')),(Path(x['distance_owners'][0]['path']),lambda t:t.replace('STATE=OK','STATE=ACTIVE').replace(UNC4_RUN,'foreign'))):
            previous=path.read_bytes();path.write_text(transform(previous.decode()))
            try:check()
            except ValueError:pass
            else:raise AssertionError('bad original scope/WD/owner accepted')
            path.write_bytes(previous)
        def resealed_status_rejected(path, transform):
            previous=path.read_bytes();science_bytes=x['science'].read_bytes();review_bytes=x['review'].read_bytes();control_bytes=control.read_bytes();mode=control.stat().st_mode
            path.write_text(transform(previous.decode()))
            science=json.loads(science_bytes)
            for key,item in science['original_evidence'].items():
                if item['path']==str(path):science['original_evidence'][key]=proof(path)
            write(x['science'],json.dumps(science));review=json.loads(review_bytes);review['receipt_sha256']=sha(x['science']);write(x['review'],json.dumps(review))
            ctrl=json.loads(control_bytes);ctrl.update(science_receipt=proof(x['science']),review=proof(x['review']));control.chmod(0o600);write(control,json.dumps(ctrl))
            try:
                try:check()
                except ValueError:pass
                else:raise AssertionError('resealed foreign watchdog/aggregate lineage accepted')
            finally:
                path.write_bytes(previous);x['science'].write_bytes(science_bytes);x['review'].write_bytes(review_bytes);control.write_bytes(control_bytes);control.chmod(mode)
        resealed_status_rejected(x['wd'],lambda t:t.replace('ARRAY_JOB_ID=4666187','ARRAY_JOB_ID=999'))
        resealed_status_rejected(x['wd'],lambda t:t.replace('LABEL=batch_effect_uncorrected__mrvi','LABEL=foreign'))
        resealed_status_rejected(x['aggregate'],lambda t:t.replace('WATCHDOG_LABELS=batch_effect_uncorrected__mrvi','WATCHDOG_LABELS=foreign'))
        def blank_generation(t):
            gate=json.loads(t);gate['event_generation']='';gate['terminal_status']['event_generation']='';return json.dumps(gate)
        def nonterminal_exit(t):
            gate=json.loads(t);gate['terminal_status']['state']='RUNNING';gate['terminal_status']['exit_code']=None;return json.dumps(gate)
        resealed_status_rejected(x['gate_path'],blank_generation)
        resealed_status_rejected(x['gate_path'],nonterminal_exit)
        # Direct JSON copying cannot bypass the fixed typed control.
        runtime=Path(x['companions'][0]['runtime']['path']);destination=run_root/'artifacts/forbidden.json'
        direct=[sys.executable,str(x['x']['tree']/'src/utils/py/bl_science_publication.py'),'copy','--source',str(runtime),'--destination',str(destination),'--boundary',str(run_root),'--md5',hashlib.md5(runtime.read_bytes()).hexdigest(),'--size',str(runtime.stat().st_size)]
        assert subprocess.run(direct,capture_output=True).returncode!=0 and not destination.exists()
        assert not (x['base']/'forbidden').exists()
    print('PASS fixed unc4 dist/runtime/timing12 publication; original failed23 outcomes/owners preserved; JSON birthPATH unchanged')


if __name__=='__main__':run()
