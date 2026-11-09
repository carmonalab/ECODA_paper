"""One fixed login-reboot publication; no fabricated wrapper/job outcomes."""
import copy
import hashlib
import json
import os
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


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-reboot-publication-') as temporary:
        base = Path(temporary).resolve()
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
