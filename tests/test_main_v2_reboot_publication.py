"""One fixed login-reboot publication; no fabricated wrapper/job outcomes."""
import copy
import hashlib
import json
import os
import shlex
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_main_v2_science_publication import build, proof, sha, write
from bl_science_publication import MAIN_RUNS, REBOOT_SCHEMA, verified_scope


def fixture(base, context):
    x = build(base, context)
    old, parent, scratch = (x[key] for key in ('old', 'parent', 'scratch'))
    names = dict(benchmark=MAIN_RUNS['benchmark'], uncorrected=MAIN_RUNS['batch_effect/uncorrected_final'], corrected=MAIN_RUNS['batch_effect/corrected_final'])
    snapshots = []
    ledgers = []
    initial = []
    serial = 10000
    for label, child in names.items():
        path = scratch / '_ecoda_runs' / child / 'manifests/scheduler_ids.tsv'
        selected = child == old.name
        rows = [['ARRAY', '101'], ['WATCHDOG', '102'], ['ARRAY', '103'], ['WATCHDOG', '104'], ['AGGREGATE_GATE', '105']] if selected else []
        target = dict(benchmark=38, uncorrected=17, corrected=16)[label]
        while len(rows) < target:
            rows.append(['STATUS', str(serial)]); serial += 1
        write(path, ''.join('\t'.join(row) + '\n' for row in rows))
        write(Path(str(path) + '.md5'), f'PATH={path}\nMD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\n')
        ledgers.append(dict(context=label, **proof(path), sidecar_sha256=sha(Path(str(path) + '.md5'))))
        initial.extend([[label, child, *row] for row in rows])
        if selected:
            for wd in ('102', '104'):
                log = write(old / f'logs/5_matrix_watchdog_method{wd}_{wd}.log', 'MATRIX_RETRY_ARRAY_JOB_ID=990001\n' if wd == '102' else '')
                err = write(old / f'logs/5_matrix_watchdog_method{wd}_{wd}.err', '')
                status = write(old / f'status/watchdogs/method{wd}.status', f'STATE=FAIL\nLABEL=method{wd}\nSCHEDULER_ID={101 if wd == "102" else 103}\n' + ('SCHEDULER_ID=990001\n' if wd == '102' else ''))
                snapshots.append(dict(context=label, watchdog_id=wd, stdout=dict(**proof(log), exists=True), stderr=dict(**proof(err), exists=True), status=dict(**proof(status), exists=True)))
    initial += [['benchmark', names['benchmark'], 'AGGREGATE_GATE', '4665612'], ['corrected', names['corrected'], 'AGGREGATE_GATE', '4665633']]
    root = scratch / '_ecoda_inputs/reboot'
    initial_path = write(root / 'initial.tsv', ''.join('\t'.join(row) + '\n' for row in initial))
    initial_inventory = write(root / 'initial.json', json.dumps(dict(schema='ecoda_orphaned_main_initial_scheduler_inventory_v1', initial_id_count=73, final_emitted_catalog_complete=False, original_files_written=False, context_counts=dict(benchmark=39, uncorrected=17, corrected=17), catalog_sha256=sha(initial_path), original_child_ledgers=ledgers)))
    label = next(key for key, value in names.items() if value == old.name)
    emitted = initial + [[label, old.name, 'RETRY_ARRAY', '990001']]
    emitted_path = write(root / 'emitted.tsv', ''.join('\t'.join(row) + '\n' for row in emitted))
    ids = [row[3] for row in emitted]
    final_inventory = write(root / 'emitted.json', json.dumps(dict(schema='ecoda_orphaned_main_emitted_scheduler_inventory_v1', final_emitted_catalog_complete=True, original_files_written=False, catalog_sha256=sha(emitted_path), initial_catalog=proof(initial_path), initial_inventory=proof(initial_inventory), single_inspect_scheduler_ids=ids, original_child_ledgers=ledgers, watchdog_evidence=snapshots)))
    gate = json.loads(x['gate_path'].read_text())
    # Real old dispatcher argv can name only external controls, with RUN_ID in env.
    external = []
    for line in x['context_file'].read_text().splitlines():
        label, child, selected_path, refs_path = line.split('\t')
        selection_path = write(root / 'original-controls' / (label + '.selection.tsv'), Path(selected_path).read_text())
        references_path = write(root / 'original-controls' / (label + '.inputs.tsv'), Path(refs_path).read_text())
        external.append('\t'.join((label, child, str(selection_path), str(references_path))) + '\n')
    original_contexts = write(root / 'original-controls/contexts.tsv', ''.join(external))
    gate['exact_command'] = shlex.join([*shlex.split(gate['exact_command']), '--contexts', str(original_contexts)])
    gate['command_digest'] = hashlib.sha256(gate['exact_command'].encode()).hexdigest()
    gate.update(state='PRELAUNCH_STOP', terminal_status=None, remote_exit_code=None, remote_host='bamboo', launch_timestamp='2026-11-09T08:18:46Z', tmux_session='ecoda-old-main-v2')
    for key, suffix in (('remote_manifest', '.manifest.json'), ('remote_runner', '.runner.sh'), ('remote_status', '.status.json'), ('remote_failure_status', '.status.json.failure')):
        gate[key] = str(scratch / 'gates' / (gate['gate_id'] + suffix))
    launch = {key: gate[key] for key in ('gate_id', 'exact_command', 'command_digest', 'event_generation', 'remote_host')}
    write(Path(gate['remote_manifest']), json.dumps(launch)); write(Path(gate['remote_runner']), gate['event_generation'] + '\n' + gate['exact_command'] + '\n')
    files = [dict(**proof(Path(gate[key])), text=Path(gate[key]).read_text(), exists=True, truncated=False) for key in ('remote_manifest', 'remote_runner')]
    files += [dict(path=gate[key], exists=False) for key in ('remote_status', 'remote_failure_status')]
    files.append(dict(path=str(scratch / 'gates' / (gate['gate_id'] + '.done')), exists=False))
    obs = dict(schema='ecoda_login_reboot_operational_observation_v1', observed_at='2026-11-09T16:05:11.164862+00:00', hostname='login1.bamboo', boot_epoch=1794235647, boot_utc='2026-11-09T14:47:27+00:00', no_accounting_query=True, no_scientific_payload_reads=True, no_remote_writes=True, tmux=dict(exit=1, stdout=''), processes=dict(exit=0, stdout='123 1 python -\n'), queue=dict(exit=0, stdout=''), files=files, active_job_descriptors={'4665612': dict(exit=0, stdout=f'JobId=4665612 matrix_gate.sh {names["benchmark"]}')})
    observation = write(root / 'observation.json', json.dumps(obs))
    quiet = copy.deepcopy(obs); quiet['observed_at'] = '2026-11-10T04:00:00+00:00'; quiet['tmux'] = dict(exit=0, stdout='ecoda-fresh-recovery: 1 windows\n')
    quiet['processes']['stdout'] = f'123 1 bash {base}/fresh0314/tree/src/5_run_benchmark_methods/stage5_dispatcher.sh --selection-file {base}/fresh-controls/benchmark.selection.tsv\n'
    quiescence = write(root / 'quiescence.json', json.dumps(quiet))
    agg = write(scratch / '_ecoda_runs' / names['corrected'] / 'logs/5_matrix_gate_4665633.log', '')
    aggregate = write(root / 'aggregate.json', json.dumps([dict(**proof(agg), exists=True)]))
    incident = write(root / 'incident.json', json.dumps(dict(schema='ecoda_main_v2_login_reboot_incident_v1', orchestration_outcome='LOGIN_REBOOT_LOSS', wrapper_exit_code=None, observation=proof(observation), aggregate_evidence=proof(aggregate), initial_catalog=proof(initial_path), initial_inventory=proof(initial_inventory), emitted_catalog=proof(emitted_path), emitted_inventory=proof(final_inventory), quiescence_observation=proof(quiescence))))
    raw = ''.join(f'{job}|{"COMPLETED|0:0" if job in {"101", "103"} else "FAILED|0:53"}\n' for job in ids)
    accounting = dict(attempted=True, query_count=1, actual_exit=0, scheduler_ids=ids, command=f'sacct -n -P -X -j {",".join(ids)} --format=JobIDRaw,State,ExitCode', stdout=raw, rows=[dict(line=i, job_id=row[0], state=row[1], exit_code=row[2]) for i, row in enumerate((line.split('|') for line in raw.splitlines()), 1)])
    gate.update(scheduler_array_ids=ids, scheduler_watchdog_ids=[], audit=dict(passed=False, accounting=accounting))
    write(x['gate_path'], json.dumps(gate)); write(x['inspection'], json.dumps(dict(operation='inspect', state='PRELAUNCH_STOP', gate_id=gate['gate_id'], audit=gate['audit'])))
    write(x['raw_path'], raw)
    # No invented parent scheduler ledger or original wrapper terminal status.
    x['campaign_sched'].unlink()
    for name in ('terminal', 'aggregate'):
        (old / 'status' / name).unlink()
    science = json.loads(x['receipt'].read_text())
    science.update(schema='ecoda_main_v2_reboot_science_receipt_v1', original_gate_status='PRELAUNCH_STOP', original_accounting=accounting, incident_sha256=sha(incident), orchestration_outcome='LOGIN_REBOOT_LOSS', wrapper_exit_code=None)
    science['original_evidence'].update(gate_manifest=proof(x['gate_path']), inspect=proof(x['inspection']), accounting_raw=proof(x['raw_path']), campaign_scheduler_ids=proof(emitted_path), child_scheduler_ids=proof(old / 'manifests/scheduler_ids.tsv'))
    science['status_snapshots'] = [dict(path=str(old / 'status' / name), exists=False) for name in ('terminal', 'aggregate')]
    write(x['receipt'], json.dumps(science))
    review = json.loads(x['review'].read_text()); review.update(receipt_sha256=sha(x['receipt']), original_gate_status='PRELAUNCH_STOP', incident_sha256=sha(incident), orchestration_outcome='LOGIN_REBOOT_LOSS', wrapper_exit_code=None)
    write(x['review'], json.dumps(review))

    def invoke(name):
        control = write(root / (name + '.json'), json.dumps(dict(schema=REBOOT_SCHEMA, context=context, run_id=name, release_id='reboot-release', nas_root=str(x['nas']), science_receipt=proof(x['receipt']), review=proof(x['review']), incident=proof(incident))))
        return subprocess.run(['bash', str(x['tree'] / 'src/utils/bash/publish_bl_science.sh'), '--control', str(control), '--control-md5', hashlib.md5(control.read_bytes()).hexdigest()], env=dict(x['env'], ECODA_RUN_ID=name, ECODA_RUN_ROOT=str(scratch / '_ecoda_runs' / name)), capture_output=True, text=True)
    return dict(x, **{key: value for key, value in locals().items() if key not in x}) | {'invoke': invoke}


def source_precondition_tests(base):
    """Real file/proof/stat boundaries, with synthetic source bytes explicitly isolated."""
    import bl_science_publication as mod
    x = fixture(base, 'benchmark')
    for snap in x['snapshots']:
        for kind in ('stdout','stderr','status'):
            Path(snap[kind]['path']).unlink()
    scratch=x['scratch']; roots={short:scratch/'_ecoda_runs'/child for short,child in x['names'].items()}
    original_sources=mod.REBOOT_RETRY_SOURCES
    # The fixed production allowlist is20 files; fixture bytes never stand in for real9ca approval.
    assert len(original_sources)==20
    synthetic_sources={}
    for relative in original_sources:
        path=x['source']/relative
        if not path.exists(): write(path, '# bounded synthetic source fixture\n'+relative+'\n')
        synthetic_sources[relative]=sha(path)
    mod.REBOOT_RETRY_SOURCES=synthetic_sources
    try:
        pairs=[]; initial=[]; ledgers=[]; serial=20000
        for short,root in roots.items():
            context=next(k for k,v in MAIN_RUNS.items() if v==root.name)
            rows=[['PREFLIGHT',str(serial)],['PREFLIGHT',str(serial+1)]];serial+=2
            for index in range(18 if short=='benchmark' else 7):
                if short=='benchmark' and index<2: array,wd=str(101+2*index),str(102+2*index)
                else: array,wd=str(serial),str(serial+1);serial+=2
                label=(('benchmark_analysis' if short=='benchmark' else 'batch_effect_'+short)+'__fixture'+str(index))
                matrix=write(root/'manifests'/('matrix_'+label.replace('__','_')+'.tsv'),'Stephenson\tbenchmark_analysis\tgloscope\n')
                write(Path(str(matrix)+'.md5'),f'MD5={hashlib.md5(matrix.read_bytes()).hexdigest()}\nSIZE={matrix.stat().st_size}\nPATH={matrix}\n')
                pair=dict(context=context,run_root=str(root),label=label,array_id=array,watchdog_id=wd,matrix=str(matrix));pairs.append(pair)
                rows.extend([['ARRAY',array],['WATCHDOG',wd]])
                #14 exact native terminal0 statuses,18 absent records as in the incident.
                if len(pairs)<=14:
                    write(root/'status/watchdogs'/(label+'.status'),f'STATE=FAIL\nLABEL={label}\nREASON=fixture\nARRAY_JOB_ID={array}\nRETRY_INDEX=0\nSCHEDULER_ID={array}\n')
                    write(root/'logs'/('5_matrix_watchdog_'+label+'_'+wd+'.log'),'')
                    write(root/'logs'/('5_matrix_watchdog_'+label+'_'+wd+'.err'),'')
            if short=='uncorrected':rows.append(['AGGREGATE_GATE','105'])
            path=write(root/'manifests/scheduler_ids.tsv',''.join('\t'.join(row)+'\n' for row in rows))
            write(Path(str(path)+'.md5'),f'PATH={path}\nMD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\n')
            ledgers.append(dict(context=short,**proof(path),sidecar_sha256=sha(Path(str(path)+'.md5'))))
            initial.extend([[short,root.name,*r] for r in rows])
            for filename in ('source.manifest','runtime.identity'):
                if root!=x['old']:write(root/'manifests'/filename,(x['old']/'manifests'/filename).read_text())
            log=''.join(f'BENCHMARK_ARRAY_JOB_ID={p["label"]}:{p["array_id"]}\nBENCHMARK_WATCHDOG_JOB_ID={p["label"]}:{p["watchdog_id"]}\n' for p in pairs if p['run_root']==str(root))
            write(x['parent']/'logs'/(short+'.log'),log)
        initial += [['benchmark',roots['benchmark'].name,'AGGREGATE_GATE','4665612'],['corrected',roots['corrected'].name,'AGGREGATE_GATE','4665633']]
        assert len(initial)==73
        write(x['initial_path'],''.join('\t'.join(r)+'\n' for r in initial))
        initial_proof=json.loads(x['initial_inventory'].read_text());initial_proof.update(catalog_sha256=sha(x['initial_path']),original_child_ledgers=ledgers)
        write(x['initial_inventory'],json.dumps(initial_proof));write(x['emitted_path'],x['initial_path'].read_text())
        def snapshot(path):
            if not path.exists():return dict(path=str(path),exists=False)
            st=path.stat();return dict(**proof(path),exists=True,size=st.st_size,text=path.read_text(),stat={k:getattr(st,k) for k in ('st_dev','st_ino','st_mode','st_size','st_mtime_ns','st_ctime_ns')})
        snapshots={}; groups=[];watches=[]
        def save(path):
            item=snapshot(path);snapshots[str(path)]=item;return item
        for root in roots.values():
            for filename in ('source.manifest','runtime.identity','scheduler_ids.tsv','scheduler_ids.tsv.md5'):save(root/'manifests'/filename)
        for short in roots:save(x['parent']/'logs'/(short+'.log'))
        for pair in pairs:
            root=Path(pair['run_root']);label=pair['label'];wd=pair['watchdog_id'];save(Path(pair['matrix']));save(Path(pair['matrix']+'.md5'))
            files={k:save(root/('status/watchdogs/'+label+'.status' if k=='status' else 'logs/5_matrix_watchdog_'+label+'_'+wd+suffix)) for k,suffix in [('status',''),('stdout','.log'),('stderr','.err')]}
            groups.append(dict(pair=pair,files={k:v['path'] for k,v in files.items()}))
            watches.append(dict(context=next(k for k,v in roots.items() if str(v)==pair['run_root']),watchdog_id=wd,**files))
            for i in range(1,5):
                for suffix in ('','_corrected_final'):
                    path=root/'manifests'/(label+suffix+'.retry_'+str(i)+'.tsv');save(path);save(Path(str(path)+'.md5'))
        directories=[]
        for root in roots.values():
            for name in ('logs','manifests','status/watchdogs'):
                path=root/name;path.mkdir(parents=True,exist_ok=True)
                entries=[dict(name=p.name,mode=p.lstat().st_mode,size=p.lstat().st_size,symlink=p.is_symlink()) for p in path.iterdir()]
                directories.append(dict(path=str(path),exists=True,entries=sorted(entries,key=lambda a:a['name'])))
        observation=write(x['root']/'coverage.observed.json',json.dumps(dict(schema='ecoda_old_main_catalog32_observation_v1',status='OBSERVED_NOT_APPROVED',initial_ids=73,collector_sha256='9867940459b07962504c35d66a6fcd82aad2de8836f4248f6c79509fab2387f8',control_sha256='ef5d687d045890b31d769bee25bd04de09d4919726bf4c38b8409712dbec9a6a',stable=True,final_emitted_catalog_complete=False,scientific_reuse_approved=False,no_accounting_query=True,no_queue_query=True,no_binary_payload_reads=True,no_remote_writes=True,observed_at_epoch=1,files=list(snapshots.values()),groups=groups,directories=directories)))
        common=dict(original_source_commit=mod.MAIN_SOURCE,original_run_roots=[str(r) for r in roots.values()],observation_sha256=sha(observation),initial_catalog_sha256=sha(x['initial_path']),original_gate_state='PRELAUNCH_STOP',wrapper_exit_code=None,reviewed_at_epoch=2)
        preserve=write(x['root']/'preservation.json',json.dumps(dict(common,schema='ecoda_main_v2_retry_manifest_preservation_review_v1',status='APPROVED_SOURCE_PRECONDITION_NO_RETRY_COVERAGE',source_pins_sha256=hashlib.sha256(json.dumps(synthetic_sources,sort_keys=True,separators=(',',':')).encode()).hexdigest(),source_preservation_confirmed=True,no_retry_manifest_deletion_or_controller_reexecution=True)))
        writer=write(x['root']/'writers.json',json.dumps(dict(common,schema='ecoda_main_v2_scoped_quiescence_review_v1',status='APPROVED_NO_OLD_ROOT_WRITERS',quiescence_observation_sha256=sha(x['quiescence']),no_old_root_writers=True,no_unknown_writers=True,allowed_disjoint_source_commit='32769f1cde54346e8b003ad096ca9e5b7490a517',allowed_disjoint_run_roots=[str(scratch/'_ecoda_runs'/('stage5_remaining95_'+part+'_20261109_v1')) for part in ('benchmark','uncorrected','corrected')]+[str(scratch/'_ecoda_runs/stage5_remaining95_20261109_v1')],evidence=[proof(x['quiescence'])])))
        coverage_path=write(x['root']/'coverage.json',json.dumps(dict(schema='ecoda_main_v2_no_retry_coverage_v1',status='APPROVED_NO_NEW_RETRY_IDS',source_commit=mod.MAIN_SOURCE,initial_catalog_sha256=sha(x['initial_path']),observation=proof(observation),preservation_review=proof(preserve),writer_review=proof(writer))))
        final=json.loads(x['final_inventory'].read_text());final.update(catalog_sha256=sha(x['emitted_path']),initial_catalog=proof(x['initial_path']),initial_inventory=proof(x['initial_inventory']),single_inspect_scheduler_ids=[r[3] for r in initial],original_child_ledgers=ledgers,watchdog_evidence=watches);write(x['final_inventory'],json.dumps(final))
        incident=json.loads(x['incident'].read_text());incident.update(initial_catalog=proof(x['initial_path']),initial_inventory=proof(x['initial_inventory']),emitted_catalog=proof(x['emitted_path']),emitted_inventory=proof(x['final_inventory']),no_retry_coverage=proof(coverage_path));write(x['incident'],json.dumps(incident))
        def check():
            ctrl=dict(incident=proof(x['incident']));sc=dict(incident_sha256=sha(x['incident']),orchestration_outcome='LOGIN_REBOOT_LOSS',wrapper_exit_code=None)
            return mod.verify_login_reboot(ctrl,sc,sc,json.loads(x['gate_path'].read_text()),{},dict(campaign_scheduler_ids=x['emitted_path']),initial,[r[3] for r in initial],x['old'])
        # No file/status/owner mutation on the diagnostic acceptance path.
        before={p:(sha(p),p.stat().st_mode) for root in roots.values() for p in root.rglob('*') if p.is_file()}
        covered=check();assert len(covered['_watchdogs'])==32
        assert all((sha(p),p.stat().st_mode)==v for p,v in before.items())
        def negative(path,change,expected):
            data=path.read_bytes();change()
            try:check()
            except (ValueError,KeyError,OSError) as e:assert expected in str(e),(expected,str(e))
            else:raise AssertionError('accepted '+expected)
            finally:path.write_bytes(data)
        # Reseal review references so failures exercise semantic approval checks.
        old_coverage=coverage_path.read_bytes();old_incident=x['incident'].read_bytes()
        def reseal_review(path,changes):
            value=json.loads(path.read_text());value.update(changes);write(path,json.dumps(value));c=json.loads(coverage_path.read_text());c['writer_review' if path==writer else 'preservation_review']=proof(path);write(coverage_path,json.dumps(c));i=json.loads(x['incident'].read_text());i['no_retry_coverage']=proof(coverage_path);write(x['incident'],json.dumps(i))
        for path,changes,expected in [(preserve,dict(source_preservation_confirmed=False),'preservation approval'),(writer,dict(no_old_root_writers=False),'no-old-writer approval')]:
            data=path.read_bytes();negative(path,lambda:reseal_review(path,changes),expected);coverage_path.write_bytes(old_coverage);x['incident'].write_bytes(old_incident)
        changed_source=x['source']/next(iter(synthetic_sources));negative(changed_source,lambda:changed_source.write_text('altered'), 'retry source changed')
        pair=pairs[-1];new_retry=Path(pair['run_root'])/'manifests'/(pair['label']+'.retry_1.tsv');write(new_retry,'new\n')
        try:
            try:check()
            except ValueError as e:assert 'absent no-retry path appeared' in str(e)
            else:raise AssertionError('new retry accepted')
        finally:new_retry.unlink()
        # Even a resealed approval cannot contradict the native retry counter.
        status_path=Path(pairs[0]['run_root'])/'status/watchdogs'/(pairs[0]['label']+'.status')
        old_status=status_path.read_bytes();old_observation=observation.read_bytes();old_preserve=preserve.read_bytes();old_writer=writer.read_bytes()
        status_path.write_text(status_path.read_text().replace('RETRY_INDEX=0','RETRY_INDEX=1'))
        value=json.loads(observation.read_text())
        for idx,item in enumerate(value['files']):
            if item['path']==str(status_path):value['files'][idx]=snapshot(status_path)
        write(observation,json.dumps(value))
        for path in (preserve,writer):
            value=json.loads(path.read_text());value['observation_sha256']=sha(observation);write(path,json.dumps(value))
        value=json.loads(coverage_path.read_text());value.update(observation=proof(observation),preservation_review=proof(preserve),writer_review=proof(writer));write(coverage_path,json.dumps(value))
        value=json.loads(x['incident'].read_text());value['no_retry_coverage']=proof(coverage_path);write(x['incident'],json.dumps(value))
        try:
            try:check()
            except ValueError as e:assert 'Native watchdog contradicts' in str(e),str(e)
            else:raise AssertionError('contradictory native retry accepted')
        finally:
            status_path.write_bytes(old_status)
            # Restore the exact baseline stat witness after the intentionally modified file.
            observation.write_bytes(old_observation);preserve.write_bytes(old_preserve);writer.write_bytes(old_writer)
            coverage_path.write_bytes(old_coverage);x['incident'].write_bytes(old_incident)
        # Explicit producer-design arm: real PUBLISHED record/proof and family owner, no fabricated global artifact owner.
        original=x['old'];family_key='derivative/'+original.name+'/ordinary/Stephenson/benchmark_analysis/gloscope'
        import re
        family=write(scratch/'_ecoda_owners/stage5'/re.sub(r'[/,:\t |]','_',family_key)/'owner',f'RUN_ID={original.name}\nSTAGE=stage5\nKEY={family_key}\nSTATE=ACTIVE\n')
        names={'Stephenson_GloScope_'+t+'.rds' for t in ('hvg2000_pcadims10','hvg2000_pcadims50','hvg1000_pcadims30','hvg3000_pcadims30')}
        review=write(x['root']/'ownerless.json',json.dumps(dict(schema='ecoda_main_v2_ownerless_gloscope_review_v1',status='APPROVED_EXACT_GLOSCOPE_OWNERLESS_PRODUCER_DESIGN',original_source_commit=mod.MAIN_SOURCE,original_run_root=str(original),paths=[str(original/'artifacts/benchmark/results'/n) for n in names],producer='stage5_gloscope',no_original_owner_created=True,writer_review_sha256=sha(writer))))
        covered['ownerless_gloscope_review']=proof(review)
        for name in names:
            path=write(original/'artifacts/benchmark/results'/name,'opaque already-science-validated synthetic bundle\n');md5=hashlib.md5(path.read_bytes()).hexdigest();size=str(path.stat().st_size)
            write(Path(str(path)+'.md5'),f'PATH={path}\nMD5={md5}\nSIZE={size}\n');record=write(original/'manifests/artifacts'/(hashlib.sha256(str(path).encode()).hexdigest()[:32]+'.record'),f'PATH={path}\nSIZE={size}\nMD5={md5}\nRUN_ID={original.name}\nPRODUCER=stage5_gloscope\nSTATE=PUBLISHED\n')
            mod.record(path,original,'stage5_gloscope')
            row=dict(source=str(path),kind='result',original_producer='stage5_gloscope');item=dict(owner=None);checked=dict(row=['Stephenson','benchmark_analysis','gloscope',path.stem.removeprefix('Stephenson_GloScope_')],sample_ids=['s'+str(i) for i in range(53)],dropped_sample_ids=[],artifact_contract_valid=True,current_labels_valid=True,sample_order_valid=True,timing_valid=True)
            mod.verify_ownerless_gloscope(row,item,'benchmark',original,checked,covered)
            for key in ('artifact_contract_valid','current_labels_valid','sample_order_valid','timing_valid'):
                bad=dict(checked,**{key:False})
                try:mod.verify_ownerless_gloscope(row,item,'benchmark',original,bad,covered)
                except ValueError:pass
                else:raise AssertionError('bad science accepted '+key)
            owner_path=mod.owner_file(path,original);write(owner_path,'RUN_ID=foreign\nSTATE=ACTIVE\n')
            try:
                try:mod.verify_ownerless_gloscope(row,item,'benchmark',original,checked,covered)
                except ValueError:pass
                else:raise AssertionError('foreign owner accepted')
            finally:owner_path.unlink()
            bad=dict(row,source=str(path).replace('pcadims10','pcadims30') if 'pcadims10' in name else str(path)+'other')
            try:mod.verify_ownerless_gloscope(bad,item,'benchmark',original,checked,covered)
            except ValueError:pass
            else:raise AssertionError('unlisted ownerless path accepted')
        assert family.read_text().endswith('STATE=ACTIVE\n')
        print('source-precondition32 and exact ownerlessGlo4 tests passed')
    finally:
        mod.REBOOT_RETRY_SOURCES=original_sources


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-reboot-publication-') as temporary:
        base = Path(temporary).resolve()
        source_precondition_tests(base / 'source-precondition')
        for number, context in enumerate(MAIN_RUNS):
            x = fixture(base / str(number), context)
            preserved = {p: (sha(p), p.stat().st_mode) for root in (x['old'].parent, x['scratch'] / '_ecoda_owners', x['scratch'] / 'gates') for p in root.rglob('*') if p.is_file()}
            result = x['invoke']('reboot_good')
            assert result.returncode == 0, result.stdout + result.stderr
            run_root = x['scratch'] / '_ecoda_runs/reboot_good'
            report = json.loads((run_root / 'artifacts/evidence/main_v2_reboot_science_publication.json').read_text())
            assert report['original_gate_status'] == 'PRELAUNCH_STOP' and report['wrapper_exit_code'] is None
            assert report['orchestration_outcome'] == 'LOGIN_REBOOT_LOSS' and report['worker_exit_not_established'] is True
            assert not report['sample_metadata_published'] and not x['campaign_sched'].exists()
            for path, previous in preserved.items():
                assert (sha(path), path.stat().st_mode) == previous, str(path)
            assert not (x['base'] / 'forbidden').exists()
            control = run_root / 'manifests/control.json'
            def check():
                return verified_scope(SimpleNamespace(control=str(control), control_md5=hashlib.md5(control.read_bytes()).hexdigest(), run_root=str(run_root), nas_root=str(x['nas'])))
            check()
            # Reseal only the proposed review/control: structural guards must still reject.
            def rejected(path, mutate, message):
                original = path.read_bytes(); cmode = control.stat().st_mode; control.chmod(0o600); cbytes = control.read_bytes(); ibytes = x['incident'].read_bytes(); sbytes=x['receipt'].read_bytes(); rbytes=x['review'].read_bytes()
                value = json.loads(original); mutate(value); write(path, json.dumps(value))
                incident = json.loads(x['incident'].read_text())
                for key, item in incident.items():
                    if isinstance(item, dict) and item.get('path') == str(path): incident[key] = proof(path)
                write(x['incident'], json.dumps(incident))
                science = json.loads(x['receipt'].read_text()); science['incident_sha256'] = sha(x['incident'])
                for key, item in science['original_evidence'].items():
                    if item['path'] == str(path): science['original_evidence'][key] = proof(path)
                if path == x['gate_path']:
                    inspection=json.loads(x['inspection'].read_text()); inspection['state']=value['state']; inspection['audit']=value['audit']; write(x['inspection'],json.dumps(inspection)); science['original_evidence']['inspect']=proof(x['inspection'])
                write(x['receipt'],json.dumps(science)); review=json.loads(x['review'].read_text()); review.update(receipt_sha256=sha(x['receipt']),incident_sha256=sha(x['incident'])); write(x['review'],json.dumps(review))
                ctrl=json.loads(control.read_text());ctrl.update(science_receipt=proof(x['receipt']),review=proof(x['review']),incident=proof(x['incident']));write(control,json.dumps(ctrl))
                try:
                    check()
                except (ValueError, KeyError, OSError): pass
                else: raise AssertionError(message)
                path.write_bytes(original);x['incident'].write_bytes(ibytes);x['receipt'].write_bytes(sbytes);x['review'].write_bytes(rbytes);control.write_bytes(cbytes); control.chmod(cmode)
            rejected(x['quiescence'], lambda v:v['queue'].update(stdout='990001 RUNNING\n'), 'queued retry accepted')
            rejected(x['quiescence'], lambda v:v['processes'].update(stdout=MAIN_RUNS[context]), 'live controller accepted')
            old_argv = f"123 1 bash {x['source']}/src/5_run_benchmark_methods/stage5_dispatcher.sh --selection-file {x['root']}/original-controls/benchmark.selection.tsv\n"
            assert not any(run_id in old_argv for run_id in (*MAIN_RUNS.values(), 'stage5_main_confounds_20261109_v2'))
            rejected(x['quiescence'], lambda v:v['processes'].update(stdout=old_argv), 'env-only old dispatcher accepted')
            control_argv = f"123 1 bash /some-other-entry.sh --selection-file {x['root']}/original-controls/benchmark.selection.tsv\n"
            rejected(x['quiescence'], lambda v:v['processes'].update(stdout=control_argv), 'original control writer accepted')
            rejected(x['quiescence'], lambda v:v['tmux'].update(stdout='ecoda-old-main-v2: 1 windows\n'), 'old tmux accepted')
            rejected(x['observation'], lambda v:v.update(boot_epoch=1), 'wrong boot accepted')
            rejected(x['observation'], lambda v:v['files'][-1].update(exists=True), 'fabricated done accepted')
            rejected(x['final_inventory'], lambda v:v.update(final_emitted_catalog_complete=False), 'initial catalog promoted')
            rejected(x['final_inventory'], lambda v:v.update(watchdog_evidence=[]), 'retry evidence omitted')
            missing_log = Path(x['snapshots'][0]['stdout']['path']); held_log = missing_log.with_suffix('.held')
            missing_log.rename(held_log)
            try:
                rejected(x['final_inventory'], lambda v:v['watchdog_evidence'][0]['stdout'].update(exists=False), 'absent watchdog stdout claimed complete retry coverage')
            finally:
                held_log.rename(missing_log)
            rejected(x['incident'], lambda v:v.update(wrapper_exit_code=1), 'invented wrapper exit accepted')
            rejected(x['gate_path'], lambda v:v.update(state='FAILED'), 'normal failure substituted')
            # Restore inspect after the intentionally resealed gate-state negative.
            write(x['inspection'], json.dumps(dict(operation='inspect',state='PRELAUNCH_STOP',gate_id=x['gate']['gate_id'],audit=x['gate']['audit'])))
            assert not (x['base'] / 'forbidden').exists()
    print('main-v2 login-reboot publication tests passed')


if __name__ == '__main__':
    run()
