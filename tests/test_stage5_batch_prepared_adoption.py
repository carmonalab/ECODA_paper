"""Strict32 mixed-source adoption: actual tiny schema2 caches and metadata."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src/utils/py'))
import stage5_prepared_adoption as adoption
spec=importlib.util.spec_from_file_location('metadata_fixture',ROOT/'tests/test_stage5_prepared_metadata.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
from src.utils.py.derivative_artifacts import bind
UNC='batch_effect/uncorrected_final';CORR='batch_effect/corrected_final'


def witness(path):return dict(path=str(path),sha256=m.digest(path))
def put(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value)+'\n');return path


def approve(f):
    value=f['value'];artifacts=value['artifacts'];prepared=[a for a in artifacts if a['kind'] in ('pseudobulk','timing')]
    obs=put(f['base']/'cache-bytes.json',dict(rows=[dict(dataset=ds,cache=dict(path=next(a['source'] for a in prepared if a['dataset']==ds and a['kind']=='pseudobulk'),md5=next(a['md5'] for a in prepared if a['dataset']==ds and a['kind']=='pseudobulk'),size=next(a['size'] for a in prepared if a['dataset']==ds and a['kind']=='pseudobulk')),timing=dict(path=next(a['source'] for a in prepared if a['dataset']==ds and a['kind']=='timing'),md5=next(a['md5'] for a in prepared if a['dataset']==ds and a['kind']=='timing'),size=next(a['size'] for a in prepared if a['dataset']==ds and a['kind']=='timing'))) for ds in m.BATCH_DATASETS]))
    cache=put(f['base']/'cache-semantics.json',dict(schema='ecoda_unc8_cache_semantics_audit_v1',cache_semantics_validated=True,timing_semantics_validated=True,actual_cli=dict(status=0),reuse_approved=False,original_gate_state='PRELAUNCH_STOP',original_wrapper_exit_code=None,original_run_root=str(f['prep_run']),metadata=witness(obs),rows=[dict(dataset=row['dataset'],identity={key:row['cache'][key] for key in ('md5','size')},timing_file=dict(identity={key:row['timing'][key] for key in ('md5','size')})) for row in json.loads(obs.read_text())['rows']]))
    by_dataset={row['dataset']:row for row in json.loads(obs.read_text())['rows']}
    canonical=sorted([dict(dataset=row['dataset'],kind=kind,path=by_dataset[row['dataset']]['cache' if kind=='pseudobulk' else 'timing']['path'],md5=identity['md5'],size=identity['size'],sha256=m.digest(Path(by_dataset[row['dataset']]['cache' if kind=='pseudobulk' else 'timing']['path']))) for row in json.loads(cache.read_text())['rows'] for kind,identity in [('pseudobulk',row['identity']),('timing',row['timing_file']['identity'])]],key=lambda row:(row['dataset'],row['kind'],row['path']))
    cache_value=json.loads(cache.read_text())
    for row in cache_value['rows']:
        row['identity']['sha256']=m.digest(Path(by_dataset[row['dataset']]['cache']['path']))
        row['timing_file']['identity']['sha256']=m.digest(Path(by_dataset[row['dataset']]['timing']['path']))
    put(cache,cache_value)
    observation=put(f['base']/'scoped-observation.json',dict(schema='ecoda_unc8_scoped_writer_observation_v1',status='READONLY_CANDIDATE',artifact_identity_stable=True,blockers=[],scoped_no_stale_writers_approved=False,scientific_reuse_approved=False,source_run_root=str(f['prep_run']),original_gate_state='PRELAUNCH_STOP',wrapper_exit_code=None,canonical_artifacts=canonical,canonical_artifacts_sha256=adoption.scope_digest(canonical),control_proofs=[witness(obs),witness(cache)]))
    writers=put(f['base']/'writers.json',dict(schema='ecoda_stage5_batch_prepared_writer_review_v1',status='APPROVED_SCOPED_NO_STALE_WRITERS',source_run_root=str(f['prep_run']),original_gate_state='PRELAUNCH_STOP',wrapper_exit_code=None,scoped_no_stale_writers=True,artifacts_sha256=adoption.scope_digest(prepared),evidence=[witness(observation)]))
    value['evidence']=[witness(cache),f['metadata_evidence'],witness(writers)]
    review=put(f['base']/'review.json',dict(schema='ecoda_stage5_batch_prepared_review_v1',status='APPROVED_COPY_ONLY',config_sha256=value['config_sha256'],sources_sha256=adoption.scope_digest(value['sources']),artifacts_sha256=adoption.scope_digest(artifacts),evidence=value['evidence'],source_gate_states={k:v[2] for k,v in adoption.BATCH_ROLES.items()},cache_artifacts=8,timing_artifacts=8,metadata_artifacts=16,cache_semantics_validated=True,timing_semantics_validated=True,metadata_semantics_validated=True,scoped_no_stale_writers=True))
    value['reuse_review']=witness(review);m.save_receipt(f)


def fixture(base):
    run_ids={UNC:('stage5_main_uncorrected_20261109_v1',adoption.BATCH_TARGETS[UNC]),CORR:('stage5_main_corrected_20261109_v1',adoption.BATCH_TARGETS[CORR])}
    f=m.build_fixture(base,run_ids=run_ids,old_source_commit=adoption.BATCH_ROLES['uncorrected_metadata'][3]);v=f['value']
    artifacts=[a for a in v['artifacts'] if a['kind']=='metadata'];sources={}
    for context,role in ((UNC,'uncorrected_metadata'),(CORR,'corrected_metadata')):
        sources[role]=dict(v['sources'][context],context=context,gate_state='FAILED',gate_manifest=witness(put(f['base']/(role+'.gate.json'),dict(gate_id='ecoda_stage5_main_confounds_20261109_v1',state='FAILED'))))
        for a in artifacts:
            if a['context']==context:a['source_role']=role
    prep=f['scratch']/'_ecoda_runs'/adoption.BATCH_ROLES['uncorrected_preparation'][1];manifest=prep/'manifests';manifest.mkdir(parents=True)
    old_unc=f['contexts'][UNC]
    m.fields(manifest/'source.manifest',dict(FORMAT=1,SOURCE_ROOT=f['config'].parent,SOURCE_COMMIT=adoption.BATCH_ROLES['uncorrected_preparation'][3],DATASETS_SHA256=m.digest(f['config'])))
    (manifest/'runtime.identity').write_bytes((old_unc['old_run']/'manifests/runtime.identity').read_bytes())
    selection=manifest/'selection.tsv';selection.write_bytes((old_unc['old_run']/'manifests/selection.tsv').read_bytes())
    binding=bind(run_root=prep,release_id='fixture-release',nas_root=f['nas'],stage='stage5',config_path=f['config'],selection=selection,input_artifacts=old_unc['inputs'])
    sources['uncorrected_preparation']=dict(run_root=str(prep),context=UNC,gate_state='PRELAUNCH_STOP',binding=witness(binding),config=witness(f['config']),source_manifest=witness(manifest/'source.manifest'),runtime_identity=witness(manifest/'runtime.identity'),gate_manifest=witness(put(f['base']/'preparation.gate.json',dict(gate_id='ecoda_stage5_main_confounds_20261109_v2',state='PRELAUNCH_STOP',remote_exit_code=None,terminal_status=None))))
    script=f['base']/'produce_tiny_caches.R'
    script.write_text(r'''
a<-commandArgs(TRUE);root<-a[1];old<-a[2]
e<-new.env(parent=globalenv())
for(p in c('src/utils/pseudobulk.R','src/utils/artifact_records.R','src/5_run_benchmark_methods/benchmark_pipeline.R','src/5_run_benchmark_methods/benchmark_hpc_utils.R'))sys.source(file.path(root,p),e)
genes<-c('retained_A','retained_B','retained_C')
e$fit_pseudobulk_deseq2<-function(counts,...)list(norm_matrix=counts,variance_order=rownames(counts))
e$exec_time<-function(expr){force(expr);as.difftime(1,units='secs')};e$peak_rss_gb<-function()0.25
refs<-read.delim(file.path(old,'manifests/input_artifacts.tsv'),header=FALSE)
for(i in seq_len(nrow(refs))){ds<-refs[i,1];ids<-paste0(ds,c('_s1','_s2'));counts<-matrix(seq_len(6),nrow=3,dimnames=list(genes,ids));metadata<-data.frame(Sample=ids,row.names=ids)
 e$load_h5ad_sample_aggregate<-function(...,metadata_columns)list(counts=counts,metadata=metadata,sample_ids=ids)
 stem<-paste0(ds,'_batch_effect_uncorrected_final')
 values<-e$prepare_pseudobulks_hpc(refs[i,3],hvg_rank_genes=genes,variants='hvg2000',cache_stem=stem,run_id=basename(old))
 e$save_rds_atomic(values$hvg2000,file.path(old,'artifacts/batch_effect/uncorrected_final/pseudobulks',paste0(stem,'_pseudobulk_hvg2000.rds')),producer='stage5_prepare_pseudobulk_hvg2000')
 e$emit_pseudobulk_timing_rows(values,ds,log_file=file.path(old,'logs',paste0('execution_times_batch_effect_uncorrected_final_prepare_pseudobulk_',ds,'.feather')))
}
''')
    env=dict(old_unc['env'],ECODA_RUN_ROOT=str(prep),ECODA_RUN_ID=prep.name,ECODA_AUX_ROOT=str(ROOT/'aux'))
    result=subprocess.run([str(ROOT/'.pixi/envs/default/bin/Rscript'),str(script),str(ROOT),str(prep)],env=env,text=True,capture_output=True)
    assert result.returncode==0,result.stdout+result.stderr
    for ds in m.BATCH_DATASETS:
        for kind,p in [('pseudobulk',prep/'artifacts'/UNC/'pseudobulks'/f'{ds}_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds'),('timing',prep/'logs'/f'execution_times_batch_effect_uncorrected_final_prepare_pseudobulk_{ds}.feather')]:
            record=prep/'manifests/artifacts'/(hashlib.sha256(str(p).encode()).hexdigest()[:32]+'.record');owner_sha=None
            if kind=='pseudobulk':
                m.common(env,'ecoda_artifact_owner_acquire "$1" stage5 "$2" 0 0 0 >/dev/null\necoda_artifact_owner_set_state "$1" FAIL original-failed-source >/dev/null\n',p,prep.name)
                owner_sha=m.digest(adoption.owner_path(p,f['scratch']))
            artifacts.append(dict(context=UNC,kind=kind,dataset=ds,view='batch_effect_uncorrected',source_role='uncorrected_preparation',source=str(p),md5=hashlib.md5(p.read_bytes()).hexdigest(),size=p.stat().st_size,sidecar_sha256=m.digest(Path(str(p)+'.md5')),record_sha256=m.digest(record),owner_sha256=owner_sha,status=None))
    for context in (UNC,CORR):
        row=f['contexts'][context];refs=[line.split('\t') for line in row['inputs'].read_text().splitlines()]
        identity=row['new_run']/'manifests/source_identity.json';identity.write_text(json.dumps(dict(schema=1,entries=[dict(dataset=a[0],view=a[1],path=a[2],md5=a[3],size=a[4],sample_ids=[a[0]+'_s1',a[0]+'_s2']) for a in refs])));m.proof(identity)
        row['env']['SOURCE_IDENTITY']=str(identity)
    f.update(prep_run=prep,metadata_evidence=v['evidence'][1],observation=v['evidence'][2],artifacts=artifacts)
    f['value']=dict(schema=adoption.BATCH_SCHEMA,state='APPROVED_COPY_ONLY',config_sha256=m.digest(f['config']),target_run_ids=dict(adoption.BATCH_TARGETS),sources=sources,artifacts=artifacts)
    approve(f);return f


def plan(f,context,valid):
    m.save_receipt(f)
    r=subprocess.run([sys.executable,str(ROOT/'src/utils/py/stage5_prepared_adoption.py'),'plan','--receipt',str(f['receipt']),'--sha256',m.digest(f['receipt']),'--config',str(f['config']),'--context',context],env=f['contexts'][context]['env'],text=True,capture_output=True)
    assert (r.returncode==0)==valid,(r.stdout,r.stderr);return r


def campaign_scope(f):
    source=(ROOT/'src/5_run_benchmark_methods/stage5_campaign.sh').read_text()
    start=source.index('import hashlib, json, os, pathlib, sys\n')
    body=source[start:source.index('\nPY\n',start)]
    parent=f['scratch']/'_ecoda_runs/campaign_preflight';manifest=parent/'manifests';manifest.mkdir(parents=True)
    rows=[]
    for name,context in (('benchmark','benchmark'),('uncorrected',UNC),('corrected',CORR)):
        current=f['contexts'][context];view=current['view'];path=f['base']/(name+'.public.tsv')
        datasets=('Stephenson',) if context=='benchmark' else m.BATCH_DATASETS
        path.write_text(''.join(f"{ds}\t{view}\t{('composition' if context=='benchmark' else view)}\n" for ds in datasets))
        child='scope_benchmark' if context=='benchmark' else adoption.BATCH_TARGETS[context]
        rows.append(f"{name}\t{child}\t{path}\t{current['inputs']}\n")
    (manifest/'campaign_requested.tsv').write_text(''.join(rows))
    (manifest/'campaign_prepared_adoption.json').write_bytes(f['receipt'].read_bytes())
    held=[]
    try:
        for context in (UNC,CORR):
            path=f['contexts'][context]['new_run'];other=path.with_name(path.name+'.held');path.rename(other);held.append((path,other))
        result=subprocess.run([sys.executable,'-',str(parent),str(ROOT),str(f['scratch']/'_ecoda_runs'),'auto'],input=body,env=dict(os.environ,NAS_TARGET_DIR=str(f['nas'])),text=True,capture_output=True)
        assert result.returncode==0,result.stdout+result.stderr
    finally:
        for path,other in held:other.rename(path)
    assert len((manifest/'campaign_contexts.tsv').read_text().splitlines())==3
    start=source.index('    adoption_context=benchmark');end=source.index('    [[ "${context}" == benchmark',start);forward=source[start:end]
    for context in ('benchmark','uncorrected','corrected'):
        code='args=(base)\n'+forward+'\nprintf "%s\n" "${args[@]}"\n'
        result=subprocess.run(['bash','-ec',code],env=dict(os.environ,context=context,PREPARED_ADOPTION=str(f['receipt']),PREPARED_ADOPTION_SHA=m.digest(f['receipt'])),text=True,capture_output=True)
        assert result.returncode==0,result.stderr
        assert ('--prepared-adoption' in result.stdout)==(context!='benchmark'),result.stdout


def main():
    with tempfile.TemporaryDirectory(prefix='ecoda-batch-prepared-') as tmp:
        f=fixture(Path(tmp));original=copy.deepcopy(f['value'])
        assert len(plan(f,UNC,True).stdout.splitlines())==24
        assert len(plan(f,CORR,True).stdout.splitlines())==8
        campaign_scope(f)
        for fault in ('role','missing','extra','global_failed','target','unapproved_review','cache_unvalidated','writers_unapproved','active_cache','ok_cache','foreign_cache','missing_cache_owner','timing_owner','status'):
            f['value']=copy.deepcopy(original);a=next(a for a in f['value']['artifacts'] if a['kind']=='pseudobulk');restore=None
            if fault=='role':a['source_role']='uncorrected_metadata'
            elif fault=='missing':f['value']['artifacts'].pop()
            elif fault=='extra':f['value']['artifacts'].append(copy.deepcopy(a))
            elif fault=='global_failed':f['value']['original_gate_state']='FAILED'
            elif fault=='target':f['value']['target_run_ids'][UNC]='foreign'
            elif fault=='unapproved_review':
                p=Path(f['value']['reuse_review']['path']);restore=(p,p.read_bytes());p.write_text(p.read_text().replace('APPROVED_COPY_ONLY','DRAFT'));f['value']['reuse_review']=witness(p)
            elif fault in ('cache_unvalidated','writers_unapproved'):
                index=0 if fault=='cache_unvalidated' else 2;p=Path(f['value']['evidence'][index]['path']);restore=(p,p.read_bytes())
                value=json.loads(p.read_text())
                if fault=='cache_unvalidated':value['cache_semantics_validated']=False
                else:value['status']='DRAFT_UNREVIEWED'
                put(p,value);f['value']['evidence'][index]=witness(p)
                review=Path(f['value']['reuse_review']['path']);value=json.loads(review.read_text());value['evidence']=f['value']['evidence'];put(review,value);f['value']['reuse_review']=witness(review)
            elif fault in ('active_cache','ok_cache','foreign_cache'):
                p=adoption.owner_path(Path(a['source']),f['scratch']);restore=(p,p.read_bytes())
                replacement='STATE=ACTIVE' if fault=='active_cache' else 'STATE=OK' if fault=='ok_cache' else 'STATE=FAIL'
                text=p.read_text().replace('STATE=FAIL',replacement)
                if fault=='foreign_cache':text=text.replace('RUN_ID='+f['prep_run'].name,'RUN_ID=foreign')
                p.write_text(text);a['owner_sha256']=m.digest(p)
            elif fault=='missing_cache_owner':a['owner_sha256']=None
            elif fault=='timing_owner':
                t=next(a for a in f['value']['artifacts'] if a['kind']=='timing');p=adoption.owner_path(Path(t['source']),f['scratch']);p.parent.mkdir(parents=True);m.fields(p,dict(PATH=t['source'],RUN_ID=f['prep_run'].name,STAGE='stage5',STATE='OK'));restore=(p,None);t['owner_sha256']=m.digest(p)
            else:
                meta=next(a for a in f['value']['artifacts'] if a['kind']=='metadata');p=Path(meta['status']['path']);restore=(p,p.read_bytes());p.write_text(p.read_text().replace('STATUS=EXPORTED','STATUS=NOOP_VALIDATED'));meta['status']=witness(p)
            if fault not in ('unapproved_review','cache_unvalidated','writers_unapproved','global_failed','target','role','missing','extra'):approve(f)
            plan(f,UNC,False)
            if restore:
                p,data=restore;p.unlink() if data is None else p.write_bytes(data)
        f['value']=copy.deepcopy(original);approve(f)
        protected=m.protected_sources(f);protected.update({p:(p.read_bytes(),p.stat().st_mode) for p in f['prep_run'].rglob('*') if p.is_file()})
        for context,count in ((UNC,24),(CORR,8)):
            result=m.adopt(f,context,validate_coordinator=True);assert result.returncode==0,result.stdout+result.stderr
            report=json.loads((f['contexts'][context]['new_run']/'manifests/prepared_adoption.json').read_text())
            assert len(report['artifacts'])==count and 'original_gate_state' not in report
            assert report['source_gate_states']['uncorrected_preparation']=='PRELAUNCH_STOP'
            for row in report['artifacts']:assert Path(row['source']).read_bytes()==Path(row['target']).read_bytes() and row['source_role'] in adoption.BATCH_ROLES
        # The current native classifier must skip all8 adopted prep rows, not fit them again.
        lines=(ROOT/'src/5_run_benchmark_methods/stage5_dispatcher.sh').read_text().splitlines();functions=[]
        for name in ('stage5_validate_path_component','benchmark_artifacts_for','benchmark_rds_group_valid','stage5_validate_reusable_artifact','benchmark_selected_artifacts_valid','stage5_recorded_producer','stage5_allowed_producers_for','stage5_producer_allowed'):
            start=lines.index(name+'() {');end=next(i for i in range(start+1,len(lines)) if lines[i]=='}');functions.extend(lines[start:end+1])
        script='source "$PROJECT_ROOT/src/utils/bash/ecoda_derivative_artifacts.sh"\necoda_derivative_load stage5\nsource "$PROJECT_ROOT/src/utils/bash/ecoda_stage5_policy.sh"\nsource "$PROJECT_ROOT/src/5_run_benchmark_methods/stage5_pending.sh"\n'
        script+='SCRIPT_DIR="$PROJECT_ROOT/src/5_run_benchmark_methods"; PASS_ARG=uncorrected; ECODA_RECOVERY_ACTIVE=0; METHOD_MATRIX_MODE=0; BENCHMARK_MATRIX_TEST=0; RDS_PREFLIGHT_DONE=""; RDS_PREFLIGHT_FAILED=""\nMETHODS=(prepare_pseudobulk pseudobulk composition pilot qot); ANALYSES=(_ecoda_none_); stage5_method_is_forced() { return 1; }\n'
        script+='\n'.join(functions)+'\n'
        script+=': > "$ECODA_RUN_ROOT/manifests/test_dispatch.tsv"\nfor method in prepare_pseudobulk pseudobulk composition pilot qot; do\nfor ds in '+ ' '.join(m.BATCH_DATASETS)+'; do printf "%s\tbatch_effect_uncorrected\t%s\n" "$ds" "$method" >> "$ECODA_RUN_ROOT/manifests/test_dispatch.tsv"; done; done\n'
        script+='stage5_classify_pending_rows "$ECODA_RUN_ROOT/manifests/test_dispatch.tsv" "$ECODA_RUN_ROOT/manifests/test_pending.tsv"\ncount=$(wc -l < "$ECODA_RUN_ROOT/manifests/test_pending.tsv" | tr -d " "); echo "NATIVE_PENDING_ROWS=$count" >&2\n[[ "$count" == 32 ]]\n! grep -q prepare_pseudobulk "$ECODA_RUN_ROOT/manifests/test_pending.tsv"\n'
        m.common(f['contexts'][UNC]['env'],script)
        for p,before in protected.items():assert (p.read_bytes(),p.stat().st_mode)==before,p
    print('PASS strict32 mixed-source adoption: real8 caches/timings+16metadata, oldFAIL/ABSENT owners/states preserved, native prep8 reusable/32 science pending; no benchmark fits/exports')


if __name__=='__main__':main()
