"""Current UNC32 checked-copy boundary: no jobs, H5AD opens or scientific fits."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace

from test_main_v2_science_publication import build, write, proof, sha
import bl_science_publication as pub

IDS = '4666239 4666241 4666242 4666244 4666247 4666249 4666253 4666256 4666259 4666261 4666265 4666269 4666286 4666288 4666297 4666306 4666315 4666277 4666287 4666325 4666334 4666336 4666345 4666347 4666351 4666357 4666243 4666246 4666248 4666252 4666255 4666258 4666260 4666264 4666267 4666268 4666296 4666305 4666314 4666323 4666324 4666333 4666335 4666344 4666346 4666349 4666356 4666365 4666366'.split()


def fixture(base):
    x = build(base, pub.UNC32_CONTEXT); old = x['scratch'] / '_ecoda_runs' / pub.UNC32_RUN
    m = old / 'manifests'; expected = pub.unc32_paths(old)
    control_root = x['scratch'] / '_ecoda_inputs/unc32'
    sm = write(m / 'source.manifest', f"SOURCE_ROOT={x['source']}\nSOURCE_COMMIT={pub.UNC32_SOURCE}\nAUX_ROOT={x['source']}/aux\n")
    runtime = write(m / 'runtime.identity', x['runtime'].read_text())
    terminal = write(old / 'status/terminal', 'STATE=FAIL\nRUN_ID='+old.name+'\nREASON=Stage 5 final RDS artifact validation failed\n')
    prep_review = write(control_root/'prep-review.json', '{}')
    prep = write(control_root/'prepared_adoption.json', json.dumps(dict(state='APPROVED_COPY_ONLY',reuse_review=proof(prep_review),evidence=[proof(prep_review)])))
    refs = []; exports = []; input_stats = []
    def stat(path): return {k:getattr(path.stat(),k) for k in ('st_dev','st_ino','st_size','st_mode','st_mtime_ns','st_ctime_ns')}
    def seal(path, producer, owned=False):
        md5 = hashlib.md5(path.read_bytes()).hexdigest(); size = str(path.stat().st_size)
        write(Path(str(path)+'.md5'), f'PATH={path}\nMD5={md5}\nSIZE={size}\n')
        rec = write(m/'artifacts'/(hashlib.sha256(str(path).encode()).hexdigest()[:32]+'.record'), f'PATH={path}\nMD5={md5}\nSIZE={size}\nRUN_ID={old.name}\nPRODUCER={producer}\nSTATE=PUBLISHED\n')
        own = pub.owner_file(path, old)
        if owned: write(own,f'PATH={path}\nRUN_ID={old.name}\nSTAGE=stage5\nSTATE=FAIL\n')
        def snap(p): return dict(path=str(p),exists=p.exists(),**(dict(sha256=sha(p),stat=stat(p)) if p.exists() else {}))
        return dict(payload=dict(path=str(path),exists=True,stat=stat(path)), sidecar=snap(Path(str(path)+'.md5')), record=snap(rec), owner=snap(own))
    observed_rows = []; qualification_rows = []
    for item in expected:
        path = Path(item['path'])
        if item['kind'] == 'runtime':
            core = str(path).removesuffix('.runtime.json')
            value = dict(schema_version=1,artifact_path=core,artifact_md5=hashlib.md5(Path(core).read_bytes()).hexdigest(),dataset=item['dataset'],method=item['method'].upper()+'_hvg2000_highres',time_secs=2.0,mem_GB=None)
            write(path,json.dumps(value))
        else: write(path,'original exported bytes\n' if item['kind'] == 'sample_metadata' else 'small derived fixture; never deserialized\n')
        observed_rows.append(seal(path,item['producer'],item['kind'] in {'core','sample_metadata'}))
        qualification_rows.append(dict(**{k:item[k] for k in ('path','kind','method','dataset')},md5=hashlib.md5(path.read_bytes()).hexdigest(),schema_valid=True,current_input_valid=True,sample_order_valid=True,timing_valid=True,named_labels_valid=True,geometry_valid=True,geometry_disposition='NOT_APPLICABLE_SAMPLE_METADATA' if str(path).endswith('_metadata.rds') else 'VALIDATED_NAMED_GEOMETRY',scoring_disposition='UNKNOWN_STORED_POLICY_REQUIRES_DISTANCE_RESCORE',original_export_lineage_verified=True))
    for i,ds in enumerate(pub.UNC32_DATASETS,1):
        h5 = write(base/'inputs'/f'{ds}.h5ad','MUST NEVER BE OPENED'); h5.chmod(0)
        input_stats.append(dict(path=str(h5),stat=stat(h5)))
        refs.append([ds,'batch_effect_uncorrected',str(h5),'a'*32,str(h5.stat().st_size),str(h5)+'.md5'])
        oldexport=write(x['scratch']/'_ecoda_runs/stage5_main_uncorrected_20261109_v1/artifacts'/pub.UNC32_CONTEXT/'metadata'/f'{ds}_sample_metadata.feather','original exported bytes\n')
        md5=hashlib.md5(oldexport.read_bytes()).hexdigest();size=oldexport.stat().st_size
        side=write(Path(str(oldexport)+'.md5'),f'PATH={oldexport}\nMD5={md5}\nSIZE={size}\n')
        status=write(x['scratch']/'_ecoda_runs/stage5_main_uncorrected_20261109_v1/status/metadata_export'/f'{ds}__batch_effect_uncorrected.status',f'STATE=OK\nSTATUS=EXPORTED\nTASK_ID={i}\nRUN_ID=stage5_main_uncorrected_20261109_v1\nOUTPUT_FILE={oldexport}\nINPUT_FILE={h5}\n')
        exports.append(dict(dataset=ds,source=str(oldexport),md5=md5,size=size,sidecar_sha256=sha(side),status=proof(status)))
    native = write(control_root/'native.json',json.dumps(dict(source_commit=pub.UNC32_SOURCE,run_id=old.name,status='EXACT32_RDS_NATIVE_DIAGNOSTIC_PASSED_NO_REUSE_APPROVAL',original_owners_and_statuses_unchanged=True,metadata={str(p):dict(stat=stat(p),sha256=sha(p)) for p in (sm,runtime,terminal)},inputs=input_stats,artifacts=[{}]*32,commands=[{}]*6)))
    scope = write(control_root/'scope.json',json.dumps(dict(source_commit=pub.UNC32_SOURCE,source_manifest_sha256=sha(sm),source_root=str(x['source']),runtime_identity_sha256=sha(runtime),config_sha256=sha(x['config']),strict32_receipt=proof(prep),contexts=[dict(context=pub.UNC32_CONTEXT,run_root=str(old),run_id=old.name,input_references=refs,adopted_source_rows=exports)],artifacts=[dict(r,context=pub.UNC32_CONTEXT) for r in expected])))
    observed = write(control_root/'metadata.json',json.dumps(dict(source_commit=pub.UNC32_SOURCE,issues=[],changes=[],files=[],rows=observed_rows)))
    pins = dict(scope=sha(scope),metadata=sha(observed),native32=sha(native),source_manifest=sha(sm))
    qualified=write(control_root/'qualification.json',json.dumps(dict(schema='ecoda_remaining97_unc32_science_validation_v1',run_id=old.name,source_commit=pub.UNC32_SOURCE,status='VALIDATED_CURRENT_UNC32_NO_PUBLICATION',native32_sha256=pins['native32'],new_model_fit=False,evidence=[proof(native)],preparation_lineage=[proof(prep),proof(prep_review)],rows=qualification_rows)))
    template='sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode'
    raw=[[v,'COMPLETED','0:0'] for v in IDS+[str(9000000+i) for i in range(92)]]
    accounting=dict(stdout=''.join('|'.join(r)+'\n' for r in raw),scheduler_ids=IDS,command=template.replace('{scheduler_ids}',','.join(IDS)),query_count=1,attempted=True,actual_exit=0,rows=[dict(line=i,job_id=r[0],state=r[1],exit_code=r[2]) for i,r in enumerate(raw,1)])
    profile=write(control_root/'profile.json',json.dumps(dict(name='ecoda',policy=dict(accounting_command=template))))
    command=f'executor --source-root {x["source"]} --run-id stage5_remaining95_20261109_v1'
    gate=dict(gate_id='ecoda_stage5_remaining97_20261109_v1',state='FAILED',audit_completed=True,audit_state='COMPLETED',exact_command=command,command_digest=hashlib.sha256(command.encode()).hexdigest(),event_generation='generation',profile='ecoda',profile_digest=sha(profile),scheduler_array_ids=IDS[:26],scheduler_watchdog_ids=IDS[26:],audit=dict(passed=True,accounting=accounting))
    gate['terminal_status']=dict(schema_version='1.0',state='FAILED',exit_code=1,started_at='2026-11-09T21:40:10Z',finished_at='2026-11-10T00:30:20Z',**{k:gate[k] for k in ('gate_id','exact_command','command_digest','event_generation')})
    gatefile=write(control_root/'gate.json',json.dumps(gate));inspect=write(control_root/'inspect.json',json.dumps(dict(gate_id=gate['gate_id'],state='FAILED',operation='inspect',audit=gate['audit'])))
    science=write(control_root/'science.json',json.dumps(dict(schema='ecoda_remaining97_unc32_science_receipt_v1',original_run_root=str(old),status='NOOP_VALIDATED',science_validated=True,publication_completed=False,source_commit=pub.UNC32_SOURCE,original_gate_status='FAILED',worker_exit_not_established=True,preparation_lineage=[proof(prep),proof(prep_review)],scope=proof(scope),metadata=proof(observed),native32=proof(native),qualification=proof(qualified),gate_manifest=proof(gatefile),inspect=proof(inspect),profile=proof(profile),original_accounting=accounting)))
    review=write(control_root/'review.json',json.dumps(dict(status='APPROVED_SCIENTIFIC_ARTIFACT_REUSE',receipt_sha256=sha(science),qualification_sha256=sha(qualified),datasets=8,physical_rows=32,payloads=104,original_gate_status='FAILED',publication_completed=False,selected_writers_quiescent=True,metadata_sha256=pins['metadata'],native32_sha256=pins['native32'])))
    # Replace fixed real evidence identities ONLY in the isolated synthetic executor.
    helper=x['tree']/'src/utils/py/bl_science_publication.py';code=helper.read_text()
    for key,value in pins.items():code=code.replace(pub.UNC32_PINS[key],value)
    helper.write_text(code)
    def invoke(name):
        control=write(control_root/(name+'.json'),json.dumps(dict(schema=pub.UNC32_SCHEMA,run_id=name,release_id='unc32-release',nas_root=str(x['nas']),science_receipt=proof(science),review=proof(review))))
        return subprocess.run(['bash',str(x['tree']/'src/utils/bash/publish_bl_science.sh'),'--control',str(control),'--control-md5',hashlib.md5(control.read_bytes()).hexdigest()],env=dict(x['env'],ECODA_RUN_ID=name,ECODA_RUN_ROOT=str(x['scratch']/'_ecoda_runs'/name)),capture_output=True,text=True)
    return locals()


class CurrentUNC32(unittest.TestCase):
    def test_scope(self):
        scope=pub.unc32_paths(Path('/scratch/_ecoda_runs')/pub.UNC32_RUN)
        from collections import Counter
        self.assertEqual(Counter(r['kind'] for r in scope),dict(core=40,runtime=16,timing=40,sample_metadata=8))
        self.assertFalse(any('/pseudobulks/' in r['path'] for r in scope))
        self.assertEqual(pub.analysis_context(dict(schema=pub.UNC32_SCHEMA)),pub.UNC32_CONTEXT)
        self.assertNotEqual(pub.campaign(dict(schema=pub.UNC32_SCHEMA)),pub.campaign(dict(schema=pub.UNC4_SCHEMA)))

    def test_real_copy_and_readonly_rejections(self):
        with tempfile.TemporaryDirectory(prefix='ecoda-unc32-publication-') as temp:
            x=fixture(Path(temp).resolve());old=x['old'];before={p:(sha(p),p.stat().st_mode) for root in (old,x['x']['scratch']/'_ecoda_owners') for p in root.rglob('*') if p.is_file()}
            r=x['invoke']('published_unc32');self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            run=x['x']['scratch']/'_ecoda_runs/published_unc32';report=run/'artifacts/evidence/remaining97_unc32_publication.json';value=json.loads(report.read_text())
            self.assertEqual(len(value['mappings']),104);self.assertEqual(len(value['runtime_companions']),16)
            self.assertFalse(value['preparation_caches_copied']);self.assertTrue(value['sample_metadata_published']);self.assertEqual(value['original_gate_status'],'FAILED')
            self.assertEqual(len((Path(value['destination_root'])/'checksums.md5').read_text().splitlines()),208)
            for item in value['runtime_companions']:
                self.assertTrue(json.loads(Path(item['copied_runtime']['scratch']).read_text())['artifact_path'].startswith(str(old)))
            for p,identity in before.items():self.assertEqual((sha(p),p.stat().st_mode),identity)
            def validate():
                c=run/'manifests/control.json'
                return subprocess.run([sys.executable,str(x['helper']),'validate','--control',str(c),'--control-md5',hashlib.md5(c.read_bytes()).hexdigest(),'--run-root',str(run),'--nas-root',str(x['x']['nas']),'--inventory',str(run/'manifests/publication_inventory.tsv'),'--report',str(report)],capture_output=True,text=True)
            self.assertEqual(validate().returncode,0)
            for path,change in ((old/'status/terminal',lambda b:b.replace(b'STATE=FAIL',b'STATE=OK')),(pub.owner_file(Path(x['expected'][0]['path']),old),lambda b:b.replace(old.name.encode(),b'foreign_run')),(Path(x['expected'][0]['path']),lambda b:b+b'corruption')):
                prior=path.read_bytes();path.write_bytes(change(prior));self.assertNotEqual(validate().returncode,0);path.write_bytes(prior)
            self.assertFalse((x['x']['base']/'forbidden').exists() if 'base' in x['x'] else (Path(temp)/'forbidden').exists())

    def test_exact_scientific_qualification(self):
        with tempfile.TemporaryDirectory() as temp:
            x=fixture(Path(temp).resolve());q=json.loads(x['qualified'].read_text());original=pub.UNC32_PINS.copy()
            try:
                pub.UNC32_PINS.update(x['pins']);science=dict(preparation_lineage=q['preparation_lineage'])
                self.assertEqual(len(pub.unc32_qualification(q,x['expected'],science)),104)
                for mutate in (lambda z:z['rows'].pop(),lambda z:z['rows'].append(z['rows'][0]),lambda z:z['rows'][0].update(named_labels_valid=False),lambda z:z['rows'][-1].update(original_export_lineage_verified=False),lambda z:z.update(source_commit=pub.MAIN_SOURCE)):
                    bad=copy.deepcopy(q);mutate(bad)
                    with self.assertRaises(ValueError):pub.unc32_qualification(bad,x['expected'],science)
            finally:pub.UNC32_PINS.clear();pub.UNC32_PINS.update(original)

    def test_current_copy_consumed_without_new_fit_path_claim(self):
        import importlib.util
        import os
        import pandas as pd
        from test_oldmain54_publication import R_P6
        with tempfile.TemporaryDirectory(prefix='ecoda-unc32-p6-') as temp:
            base=Path(temp).resolve();old=base/'_ecoda_runs'/pub.UNC32_RUN;run=base/'_ecoda_runs/new_publication'
            run.mkdir(parents=True);nas=base/'nas';nas.mkdir()
            item=next(r for r in pub.unc32_paths(old) if r['dataset']=='Alzheimer' and r['kind']=='runtime' and r['method']=='pilot')
            runtime=Path(item['path']);core=Path(str(runtime).removesuffix('.runtime.json'));core.parent.mkdir(parents=True)
            pd.DataFrame({'Sample':['s1','s2'],'s1':[0.,1.],'s2':[1.,0.]}).to_feather(core)
            md5=lambda p:hashlib.md5(p.read_bytes()).hexdigest()
            write(runtime,json.dumps(dict(schema_version=1,artifact_path=str(core),artifact_md5=md5(core),dataset='Alzheimer',method='PILOT_hvg2000_highres',time_secs=2.,mem_GB=None)))
            control=write(run/'manifests/control.json',json.dumps(dict(schema=pub.UNC32_SCHEMA,run_id=run.name)))
            inventory=[];mappings=[]
            for source in (core,runtime):
                rel=source.relative_to(old/'artifacts'/pub.UNC32_CONTEXT);targets=[]
                for boundary,target in ((run,run/'artifacts'/pub.UNC32_CONTEXT/rel),(nas,nas/rel)):
                    target.parent.mkdir(parents=True,exist_ok=True)
                    pub.copy_payload(SimpleNamespace(source=str(source),destination=str(target),boundary=str(boundary),md5=md5(source),size=source.stat().st_size,control=str(control),control_md5=md5(control)))
                    write(Path(str(target)+'.md5'),f'MD5={md5(target)}\nSIZE={target.stat().st_size}\nPATH={target}\n');targets.append(str(target))
                row=dict(source=str(source),kind='runtime' if source==runtime else 'core',md5=md5(source));inventory.append(row);mappings.append(dict(row,scratch=targets[0],nas=targets[1]))
            lineage=pub.audit54_runtime_lineage(inventory,mappings);self.assertEqual(len(lineage),1)
            copied=Path(mappings[0]['scratch'])
            # UNC has no corrected batch-contract JSON; keep the original path
            # explicitly rather than inventing a rebased new-fit runtime record.
            self.assertNotEqual(json.loads(Path(mappings[1]['scratch']).read_text())['artifact_path'],str(copied))
            evidence=write(base/'lineage.json',json.dumps(lineage));root=Path(__file__).resolve().parents[1]
            script=R_P6.replace('batch_effect_corrected','batch_effect_uncorrected').replace('MrVI_hvg2000','PILOT_hvg2000_highres').replace("origin='legacy_reuse'","origin='new'").replace("receipt$origin,'legacy_reuse'","receipt$origin,'new'")
            script=script.replace('scope$require_new<-FALSE','scope$require_new<-TRUE').replace('bad$require_new<-TRUE','bad$current_input_md5<-strrep("b",32)')
            rp=write(base/'p6.R',script)
            r=subprocess.run([str(root/'.pixi/envs/default/bin/Rscript'),'--vanilla',str(rp),str(root),str(base),str(copied),str(evidence)],env=dict(os.environ,R_HOME=str(root/'.pixi/envs/default/lib/R')),capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stdout+r.stderr);self.assertIn('CURRENT_P6_READER_PROJECTION_OK',r.stdout)
            arbitrary=write(base/'arbitrary.runtime.json','{}')
            with self.assertRaisesRegex(ValueError,'outside its fixed'):
                pub.copy_payload(SimpleNamespace(source=str(arbitrary),destination=str(nas/'bad.runtime.json'),boundary=str(nas),md5=md5(arbitrary),size=2,control=str(control),control_md5=md5(control)))

    def test_accounting_preserves_failed_wrapper(self):
        with tempfile.TemporaryDirectory() as temp:
            x=fixture(Path(temp).resolve());g=json.loads(x['gatefile'].read_text());i=json.loads(x['inspect'].read_text());p=json.loads(x['profile'].read_text())
            self.assertEqual(len(pub.unc32_accounting(g,i,p)['rows']),141)
            for change in (lambda z:z.update(state='COMPLETED'),lambda z:z['audit'].update(passed=False),lambda z:z['audit']['accounting'].update(query_count=2),lambda z:z['audit']['accounting']['rows'].pop(),lambda z:z['scheduler_array_ids'].__setitem__(0,'999999'),lambda z:(z.update(event_generation=''),z['terminal_status'].update(event_generation='')),lambda z:z['terminal_status'].pop('started_at'),lambda z:z['terminal_status'].update(finished_at=123)):
                bad=copy.deepcopy(g);change(bad)
                with self.assertRaises(ValueError):pub.unc32_accounting(bad,i,p)


if __name__=='__main__':unittest.main()
