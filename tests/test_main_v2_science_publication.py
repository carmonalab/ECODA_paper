"""Fixed main-v2 copy-only release; real publisher/owners, no Slurm/RNA/fits."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/utils/py'))
from bl_science_publication import MAIN_RUNS, MAIN_SCHEMA, MAIN_SOURCE


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value='fixture\n'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value)
    return path


def proof(path):
    return dict(path=str(path), sha256=sha(path))


def build(base, context):
    scratch = base / 'scratch'; nas = base / 'nas'; tree = base / 'publication/tree'; source = base / 'original/tree'
    old = scratch / '_ecoda_runs' / MAIN_RUNS[context]; parent = old.parent / 'stage5_main_confounds_20261109_v2'
    for path in (scratch / '_ecoda_inputs', nas, tree / 'src/utils/bash', tree / 'src/utils/py', base / 'bin'):
        path.mkdir(parents=True)
    for name in ('publish_bl_science.sh', 'ecoda_run_common.sh'):
        shutil.copyfile(ROOT / 'src/utils/bash' / name, tree / 'src/utils/bash' / name)
    for name in ('bl_science_publication.py', 'derivative_artifacts.py'):
        shutil.copyfile(ROOT / 'src/utils/py' / name, tree / 'src/utils/py' / name)
    write(tree / 'source.manifest', f'SOURCE_ROOT={tree}\n')
    write(tree / 'src/slurm_config.sh', 'export PROJECT_ROOT="$TEST_TREE" HPC_SCRATCH_DIR="$TEST_SCRATCH" NAS_TARGET_DIR="$TEST_NAS"\nexport PYTHON_BIN="$TEST_PYTHON" PYTHONDONTWRITEBYTECODE=1\n')
    write(tree / 'src/utils/bash/ecoda_runtime.sh', '''ecoda_runtime_validate_submission() {
 [[ "$1" == apptainer && "$ECODA_RUNTIME_PROFILE" == stage5 ]]
 printf 'RUNTIME_IMAGE_SHA256=new-runtime\n' > "$ECODA_RUN_ROOT/manifests/runtime.identity"
}
ecoda_runtime_validate_bound_run() { grep -qx RUNTIME_IMAGE_SHA256=new-runtime "$ECODA_RUN_ROOT/manifests/runtime.identity"; }
''')
    write(tree / 'src/utils/bash/ecoda_nas_backup.sh', '''#!/bin/bash
set -euo pipefail
[[ "$1" == preflight && "$3" == "$TEST_NAS" ]] || exit 7
[[ "${TEST_MISSING_MOUNT:-0}" != 1 ]] || exit 7
printf '{"mount_source":"fixture","nas_root":"%s"}\n' "$TEST_NAS" > "$5"
''')
    for command in ('sbatch', 'sacct', 'squeue', 'Rscript'):
        write(base / 'bin' / command, f'#!/bin/bash\ntouch "{base}/forbidden"\nexit 99\n').chmod(0o755)
    write(base / 'bin/chmod', '''#!/bin/bash
if [[ "${TEST_MODE_FAILURE:-0}" == 1 && "$1" == 666 && "$2" == "$TEST_NAS/"* ]]; then exit 0; fi
exec /bin/chmod "$@"
''').chmod(0o755)
    config = write(source / 'datasets.json', (ROOT / 'datasets.json').read_text())
    contracts = ('src/utils/bash/ecoda_run_common.sh', 'src/utils/bash/ecoda_stage5_policy.sh',
                 'src/utils/bash/ecoda_derivative_artifacts.sh', 'src/utils/py/derivative_artifacts.py',
                 'src/utils/py/benchmark_request_scope.py')
    for relative in contracts:
        target = source / relative; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    source_manifest = write(old / 'manifests/source.manifest', f'SOURCE_ROOT={source}\nSOURCE_COMMIT={MAIN_SOURCE}\nAUX_ROOT={source}/aux\nDATASETS_SHA256={sha(config)}\n')
    image = write(base / 'runtime/image.sif'); image_manifest = write(base / 'runtime/image.sif.manifest', f'IMAGE_SHA256={sha(image)}\n')
    runtime = write(old / 'manifests/runtime.identity', f'RUNTIME_IMAGE={image}\nRUNTIME_MANIFEST={image_manifest}\nRUNTIME_MANIFEST_SHA256={sha(image_manifest)}\nRUNTIME_IMAGE_SHA256={sha(image)}\nRUNTIME_IMAGE_SIZE={image.stat().st_size}\n')
    registry = json.loads(config.read_text()); ds = 'Stephenson'
    view = 'benchmark_analysis' if context == 'benchmark' else 'batch_effect_' + ('corrected' if context.endswith('/corrected_final') else 'uncorrected')
    h5ad = write(base / 'inputs' / registry[ds]['views'][view]['output_file_name'], 'NEVER READ H5AD PAYLOAD\n')
    input_md5 = hashlib.md5(h5ad.read_bytes()).hexdigest(); size = h5ad.stat().st_size
    checksum = write(base / 'inputs/expected.md5', f'PATH={h5ad}\nMD5={input_md5}\nSIZE={size}\n')
    h5ad.chmod(0)
    inputs = write(old / 'manifests/input_artifacts.tsv', f'{ds}\t{view}\t{h5ad}\t{input_md5}\t{size}\t{checksum}\n')
    identity = write(old / 'manifests/source_identity.json', json.dumps(dict(schema='ecoda_h5ad_source_identity_v1', entries=[dict(dataset=ds, view=view, path=str(h5ad), md5=input_md5, size=size, sample_ids=['s1', 's2', 's3'])])))
    methods = ('scitd', 'qot') if context == 'benchmark' else ('composition', 'qot')
    selection = write(old / 'manifests/selection.tsv', ''.join(f'{ds}\t{view}\t{method}\n' for method in methods) if context == 'benchmark' else f'{ds}\t{view}\t{view}\n')
    dispatch = write(old / 'manifests/dispatch_selection.tsv', ''.join(f'{ds}\t{view}\t{m}\n' for m in methods))
    binding = dict(schema='ecoda_derivative_artifacts_v1', stage='stage5', run_id=old.name, run_root=str(old),
                   artifact_root=str(old / 'artifacts'), nas_base=str(nas), release_id='original-release', nas_root=str(nas / 'releases/original-release' / old.name),
                   selection=str(selection), input_artifacts=str(inputs), identities={str(p): sha(p) for p in (config, source_manifest, runtime, selection, inputs)})
    if context == 'benchmark':
        requested = write(old / 'manifests/requested_selection.tsv', selection.read_text()); binding['requested_selection'] = str(requested); binding['identities'][str(requested)] = sha(requested)
    binding_path = write(old / 'manifests/artifact_bindings.json', json.dumps(binding))
    write(Path(str(binding_path) + '.md5'), f'PATH={binding_path}\nMD5={hashlib.md5(binding_path.read_bytes()).hexdigest()}\nSIZE={binding_path.stat().st_size}\n')
    preflight = write(old / 'manifests/h5ad_preflight.tsv', f'{ds}\t{view}\t{h5ad}\n')
    input_status = write(old / f'status/h5ad_preflight/{ds}__{view}.status', f'STATE=OK\nRUN_ID={old.name}\nDATASET={ds}\nVIEW={view}\nTASK_ID=1\n')
    child_sched = write(old / 'manifests/scheduler_ids.tsv', 'ARRAY\t101\nWATCHDOG\t102\nARRAY\t103\nWATCHDOG\t104\nAGGREGATE_GATE\t105\n')
    names = dict(benchmark=MAIN_RUNS['benchmark'], uncorrected=MAIN_RUNS['batch_effect/uncorrected_final'], corrected=MAIN_RUNS['batch_effect/corrected_final'])
    selected_name = next(key for key, value in names.items() if value == old.name)
    contexts = []; schedulers = []
    for label, child in names.items():
        parent_selection = write(parent / f'manifests/campaign_{label}.selection.tsv', (old / 'manifests/requested_selection.tsv').read_text() if label == selected_name and context == 'benchmark' else selection.read_text())
        parent_inputs = write(parent / f'manifests/campaign_{label}.inputs.tsv', inputs.read_text())
        contexts.append(f'{label}\t{child}\t{parent_selection}\t{parent_inputs}\n')
        if child == old.name:
            schedulers.extend(f'{label}\t{child}\t{row}\n' for row in child_sched.read_text().splitlines())
        else:
            schedulers.append(f'{label}\t{child}\tAGGREGATE_GATE\t{201 if label == "benchmark" else 202 if label == "uncorrected" else 203}\n')
    context_file = write(parent / 'manifests/campaign_contexts.tsv', ''.join(contexts))
    campaign_sched = write(parent / 'manifests/campaign_scheduler_ids.tsv', ''.join(schedulers))
    ids = [r.split('\t')[-1] for r in schedulers]
    ids = [i.strip() for i in ids]
    profile = write(base / 'evidence/profile.json', json.dumps(dict(name='fixture', policy=dict(accounting_command='sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode'))))
    raw = ''.join(f'{i}|{ "COMPLETED|0:0" if i in {"101", "103"} else "FAILED|1:0"}\n' for i in ids) + '901|FAILED|1:0\n902|COMPLETED|0:0\n'
    accounting = dict(attempted=True, query_count=1, actual_exit=0, scheduler_ids=ids, command=f'sacct -n -P -X -j {",".join(ids)} --format=JobIDRaw,State,ExitCode', stdout=raw, rows=[dict(line=i, job_id=r[0], state=r[1], exit_code=r[2]) for i, r in enumerate((line.split('|') for line in raw.splitlines()), 1)])
    immutable_source = write(source / 'source.manifest', source_manifest.read_text())
    argv = ['executor', '--run-id', parent.name, '--scratch-root', str(scratch), '--source-root', str(source), '--source-manifest', str(immutable_source), '--runtime-image', str(image), '--runtime-manifest', str(image_manifest), '--script', 'src/5_run_benchmark_methods/stage5_campaign.sh']
    gate = dict(gate_id='ecoda_stage5_main_confounds_20261109_v2', state='FAILED', audit_completed=True, audit_state='COMPLETED', audit=dict(passed=False, accounting=accounting), profile='fixture', profile_digest=sha(profile), exact_command=shlex.join(argv), command_digest=hashlib.sha256(shlex.join(argv).encode()).hexdigest(), event_generation='fixture-generation', scheduler_array_ids=ids, scheduler_watchdog_ids=[])
    gate['terminal_status'] = dict(schema_version='1.0', state='FAILED', exit_code=1, started_at='2026-11-09T08:00:00Z', finished_at='2026-11-09T09:00:00Z', **{key: gate[key] for key in ('gate_id', 'exact_command', 'command_digest', 'event_generation')})
    gate_path = write(base / 'evidence/gate.json', json.dumps(gate))
    inspection = write(base / 'evidence/inspect.json', json.dumps(dict(operation='inspect', state='FAILED', gate_id=gate['gate_id'], audit=gate['audit'])))
    raw_path = write(base / 'evidence/accounting.raw.tsv', raw)
    terminal = write(old / 'status/terminal', 'STATE=FAIL\n'); aggregate = write(old / 'status/aggregate', 'STATE=FAIL\nREASON=watchdog startup failed\n')
    missing_status = old / 'status/watchdogs/absent.status'
    row_owners = []; owner_lines = []; selected = []; checks = []; inventory = []; payloads = []

    def sealed(path, producer):
        md5 = hashlib.md5(path.read_bytes()).hexdigest(); size = str(path.stat().st_size)
        side = write(Path(str(path) + '.md5'), f'MD5={md5}\nSIZE={size}\nPATH={path}\n')
        record = write(old / 'manifests/artifacts' / (hashlib.sha256(str(path).encode()).hexdigest()[:32] + '.record'), f'PATH={path}\nMD5={md5}\nSIZE={size}\nRUN_ID={old.name}\nPRODUCER={producer}\nSTATE=PUBLISHED\n')
        return md5, size, proof(side), proof(record)

    for i, method in enumerate(methods):
        matrix = write(old / f'manifests/matrix_{view}_{method}.tsv', f'{ds}\t{view}\t{method}\n'); sealed(matrix, 'matrix_fixture')
        selected.append(dict(dataset=ds, view=view, method=method, parameter='', array_id=str(101 + i * 2), task_index=1, manifest=proof(matrix), stdout=proof(write(old / f'logs/{method}_{101 + i * 2}_1.log', 'completed worker\n')), stderr=proof(write(old / f'logs/{method}_{101 + i * 2}_1.err', ''))))
        scope = 'ordinary' if context == 'benchmark' else 'corrected' if context.endswith('/corrected_final') else 'uncorrected'
        owner_key = f'derivative/{old.name}/{scope}/{ds}/{view}/{method}'
        row_owner = write(scratch / '_ecoda_owners/stage5' / re.sub(r'[/,:\t |]', '_', owner_key) / 'owner', f'RUN_ID={old.name}\nSTAGE=stage5\nKEY={owner_key}\nSTATE=FAIL\n')
        row_owners.append(dict(**proof(row_owner), fields=dict(RUN_ID=old.name, STAGE='stage5', KEY=owner_key, STATE='FAIL'))); owner_lines.append(f'{owner_key}\t{row_owner.parent}\n')
        row_paths = []; row_timings = []
        stem = ds if context == 'benchmark' else ds + '_batch_effect_' + ('corrected' if context.endswith('/corrected_final') else 'uncorrected') + '_final'
        if method == 'scitd':
            relatives = [('result', f'results/{ds}_scitd.rds')]
        elif method == 'composition':
            relatives = [('result', f'results/{stem}_composition.rds'), ('metadata', f'results/{stem}_metadata.rds')]
            if context.endswith('/corrected_final'):
                relatives += [('result', f'results/{stem}_{name}.rds') for name in ('ECODA_authors_HR', 'ECODA_authors_HR_NULL', 'ecoda_unsupervised_harmony', 'ecoda_unsupervised_harmony_limma')]
        else:
            relatives = [('result', f'embeddings/{ds}_{token}_qot_dists.feather') for token in ('hvg2000_lowres', 'hvg1000_highres', 'hvg2000_highres', 'hvg3000_highres')] if context == 'benchmark' else [('result', f'embeddings/{stem}_hvg2000_highres_qot_dists.feather')]
        relatives.append(('timing', f'embeddings/timing_sources/{method}.feather'))
        for kind, relative in relatives:
            path = write(old / ('logs/' + Path(relative).name if kind == 'timing' else 'artifacts/' + context + '/' + relative), 'opaque reviewed synthetic science\n')
            producer = 'stage5_execution_log' if kind == 'timing' else 'stage5_' + method
            md5, size, side, record = sealed(path, producer)
            row = dict(dataset=ds, kind=kind, source=str(path), relative=relative, md5=md5, size=size, original_producer=producer)
            owner = None
            if kind != 'timing':
                name = path.name + '_' + hashlib.sha256(str(path).encode()).hexdigest()[:32]
                if len(name.encode()) > 100: name = 'artifact_' + hashlib.sha256(str(path).encode()).hexdigest()[:32]
                state = 'ACTIVE' if i else 'FAIL'
                values = dict(RUN_ID=old.name, STAGE='stage5', PATH=str(path), STATE=state, PID='123')
                owned = write(scratch / '_ecoda_owners/artifact' / name / 'owner', ''.join(f'{k}={v}\n' for k, v in values.items()))
                owner = dict(**proof(owned), fields=values)
            inventory.append(row); payloads.append(dict(**row, producer_record=record, sidecar=side, owner=owner, runtime_evidence=[proof(runtime)]))
            (row_timings if kind == 'timing' else row_paths).append(str(path))
        dropped = ['s2'] if method == 'scitd' else []
        checks.append(dict(row=[ds, view, method, ''], input_md5=input_md5, completion_criterion='VALIDATED_PUBLISHED_ARTIFACTS', worker_exit_not_established=True, artifact_contract_valid=True, current_labels_valid=True, sample_order_valid=True, timing_valid=True, sample_ids=[s for s in ['s1', 's2', 's3'] if s not in dropped], dropped_sample_ids=dropped, labels=['a', 'b'] if dropped else ['a', 'a', 'b'], expected_artifacts=row_paths, timing_artifacts=row_timings))
    owners = write(old / 'manifests/owners.tsv', ''.join(owner_lines))
    control_root = scratch / '_ecoda_inputs/main-v2-science'
    catalog = write(control_root / 'catalog.json', json.dumps(dict(schema='ecoda_main_v2_publication_catalog_v1', context=context, rows=selected)))
    inv_keys = ('dataset', 'kind', 'source', 'relative', 'md5', 'size', 'original_producer')
    inv_path = write(control_root / 'inventory.tsv', ''.join('\t'.join(row[k] for k in inv_keys) + '\n' for row in inventory))
    def rds_method(row):
        if context.endswith('/corrected_final'):
            for arm in ('ECODA_authors_HR', 'ECODA_authors_HR_NULL', 'ecoda_unsupervised_harmony', 'ecoda_unsupervised_harmony_limma'):
                if row['source'].endswith('_' + arm + '.rds'): return arm
        return methods[0]
    rds_inputs = write(base / 'validation/rds.inputs.tsv', ''.join(f'{row["source"]}\t{rds_method(row)}\t{ds}\t{view}\t{int(row["kind"] == "metadata")}\n' for row in inventory if row['source'].endswith('.rds')))
    runs = []
    for kind, name, marker in (('rds', 'validate_benchmark_rds_contract.R', 'benchmark RDS artifact-list contract OK'), ('feather', 'matrix_artifact_validator.py', 'matrix artifact contract OK')):
        script = write(source / 'src/5_run_benchmark_methods' / name)
        paths = [row['source'] for row in inventory if row['kind'] != 'timing' and row['source'].endswith('.feather')] if kind == 'feather' else [str(rds_inputs)]
        for number, path in enumerate(paths):
            item = dict(**proof(script), kind=kind, actual_exit=0, argv=[str(script), '--config', str(config), '--source-identity', str(identity), '--artifact-list' if kind == 'rds' else '--artifact', path], stdout=proof(write(base / f'validation/{kind}{number}.out', marker + '\n')), stderr=proof(write(base / f'validation/{kind}{number}.err', '')))
            if context != 'benchmark': item['argv'] += ['--batch-pass', 'corrected' if context.endswith('/corrected_final') else 'uncorrected', '--analysis-variant', 'corrected_final' if context.endswith('/corrected_final') else 'final']
            if context.endswith('/corrected_final') and kind == 'feather':
                contract = write(old / 'manifests/batch_contract_qot.json', '{"fixture":"expected QOT contract"}\n')
                item['expected_batch_contract'] = proof(contract); item['argv'] += ['--expected-batch-contract', str(contract)]
            if kind == 'rds': item['input_list'] = proof(rds_inputs)
            runs.append(item)
    validation = write(control_root / 'validation.json', json.dumps(dict(schema='ecoda_main_v2_artifact_validation_v1', catalog_sha256=sha(catalog), source_identity_sha256=sha(identity), config_sha256=sha(config), new_model_fit=False, contract_sources=[proof(source / relative) for relative in contracts], execution_evidence=[proof(write(base / 'validation/audit.R')), proof(runtime)], validator_runs=runs, rows=checks)))
    evidence = dict(gate_manifest=gate_path, inspect=inspection, profile=profile, accounting_raw=raw_path, campaign_contexts=context_file, campaign_scheduler_ids=campaign_sched, child_scheduler_ids=child_sched, source_manifest=source_manifest, runtime_identity=runtime, artifact_binding=binding_path, input_artifacts=inputs, source_identity=identity, selection=selection, dispatch_selection=dispatch, owners_manifest=owners, preflight_manifest=preflight)
    science = dict(schema='ecoda_main_v2_science_receipt_v1', context=context, status='NOOP_VALIDATED', science_validated=True, publication_completed=False, original_gate_status='FAILED', source_commit=MAIN_SOURCE, original_run_root=str(old), original_evidence={k: proof(v) for k, v in evidence.items()}, original_accounting=accounting, catalog=proof(catalog), scientific_validation=proof(validation), inventory=proof(inv_path), status_snapshots=[dict(**proof(p), exists=True) for p in (terminal, aggregate)] + [dict(path=str(missing_status), exists=False)], row_owners=row_owners, payloads=payloads, input_statuses=[dict(dataset=ds, view=view, status=proof(input_status))])
    receipt = write(control_root / 'science.json', json.dumps(science))
    review = write(control_root / 'review.json', json.dumps(dict(status='APPROVED_SCIENTIFIC_ARTIFACT_REUSE', receipt_sha256=sha(receipt), catalog_sha256=sha(catalog), inventory_sha256=sha(inv_path), context=context, selected_rows=2, payloads=len(inventory), worker_exit_not_established=True, original_gate_status='FAILED', publication_completed=False)))
    env = dict(os.environ, TEST_TREE=str(tree), TEST_SCRATCH=str(scratch), TEST_NAS=str(nas), TEST_PYTHON=sys.executable, ECODA_SOURCE_ROOT=str(tree), ECODA_SOURCE_MANIFEST=str(tree / 'source.manifest'), ECODA_SOURCE_SNAPSHOT_REQUIRED='1', PYTHONDONTWRITEBYTECODE='1', PATH=str(base / 'bin') + ':' + os.environ['PATH'])

    def invoke(name, **overrides):
        control = write(control_root / (name + '.json'), json.dumps(dict(schema=MAIN_SCHEMA, context=context, run_id=name, release_id='release', nas_root=str(nas), science_receipt=proof(receipt), review=proof(review))))
        return subprocess.run(['bash', str(tree / 'src/utils/bash/publish_bl_science.sh'), '--control', str(control), '--control-md5', hashlib.md5(control.read_bytes()).hexdigest()], env=dict(env, ECODA_RUN_ID=name, ECODA_RUN_ROOT=str(scratch / '_ecoda_runs' / name), **overrides), capture_output=True, text=True)

    return locals()


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-main-v2-publication-') as temporary:
        base = Path(temporary).resolve()
        for index, context in enumerate(MAIN_RUNS):
            x = build(base / str(index), context); invoke = x['invoke']
            immutable = {p: (sha(p), p.stat().st_mode) for root in (x['old'], x['parent'], x['source'], x['scratch'] / '_ecoda_owners') for p in root.rglob('*') if p.is_file()}
            result = invoke('good')
            assert result.returncode == 0, result.stdout + result.stderr
            run_root = x['scratch'] / '_ecoda_runs/good'; report = run_root / 'artifacts/evidence/main_v2_science_publication.json'
            published = json.loads(report.read_text())
            assert published['analysis_root'] == str(run_root / 'artifacts' / context)
            assert len(published['mappings']) == len(x['inventory']) and not published['new_model_fit'] and published['original_gate_status'] == 'FAILED'
            assert published['worker_exit_not_established'] and published['completion_criterion'] == 'VALIDATED_PUBLISHED_ARTIFACTS'
            assert published['sample_metadata_published'] is False
            assert len((Path(published['destination_root']) / 'checksums.md5').read_text().splitlines()) == 2 * len(x['inventory'])
            assert all('PRODUCER=stage5_main_v2_science_release_copy\n' in Path(row['scratch_record']['record_path']).read_text() for row in published['mappings'])
            control = run_root / 'manifests/control.json'
            command = [sys.executable, str(x['tree'] / 'src/utils/py/bl_science_publication.py'), 'validate', '--control', str(control), '--control-md5', hashlib.md5(control.read_bytes()).hexdigest(), '--run-root', str(run_root), '--nas-root', str(x['nas']), '--inventory', str(run_root / 'manifests/publication_inventory.tsv'), '--report', str(report)]
            assert subprocess.run(command, env=x['env'], capture_output=True).returncode == 0
            assert invoke('good').returncode != 0
            for p, before in immutable.items(): assert (sha(p), p.stat().st_mode) == before, p
            assert not (x['base'] / 'forbidden').exists()
            if context != 'batch_effect/corrected_final': continue
            assert invoke('mount_failed', TEST_MISSING_MOUNT='1').returncode != 0
            assert invoke('mode_failed', TEST_MODE_FAILURE='1').returncode != 0
            # Changed old owners/statuses are rejected, never repaired.
            for path in (Path(x['payloads'][0]['owner']['path']), x['terminal'], x['checksum']):
                saved = path.read_bytes(); path.write_bytes(saved + b'changed\n')
                assert invoke('changed_' + path.name.replace('.', '_')).returncode != 0
                path.write_bytes(saved)
            saved = x['raw_path'].read_bytes(); x['raw_path'].write_bytes(saved.replace(b'902|COMPLETED|0:0', b'902|RUNNING|0:0'))
            assert invoke('accounting_changed').returncode != 0; x['raw_path'].write_bytes(saved)
            # A resealed/root-approved receipt cannot bypass current scientific
            # validation or claim an unproved zero worker exit code.
            receipt_saved = x['receipt'].read_bytes(); review_saved = x['review'].read_bytes(); validation_saved = x['validation'].read_bytes()
            gate_saved = x['gate_path'].read_bytes()
            for fault in ('missing', 'nonterminal', 'stale', 'exit_mismatch', 'partial_audit'):
                gate = json.loads(gate_saved)
                if fault == 'missing': gate.pop('terminal_status')
                elif fault == 'nonterminal': gate['terminal_status']['state'] = 'RUNNING'
                elif fault == 'stale': gate['terminal_status']['event_generation'] = 'different-generation'
                elif fault == 'partial_audit': gate['audit_completed'] = False
                else: gate['terminal_status']['exit_code'] = 0
                x['gate_path'].write_text(json.dumps(gate)); science = copy.deepcopy(x['science']); science['original_evidence']['gate_manifest'] = proof(x['gate_path'])
                x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
                assert invoke('wrapper_' + fault).returncode != 0
            # A zero wrapper exit does not turn its separately failed audit
            # or original scheduler outcomes into a successful gate.
            gate = json.loads(gate_saved); gate['terminal_status'].update(state='COMPLETED', exit_code=0)
            x['gate_path'].write_text(json.dumps(gate)); science = copy.deepcopy(x['science']); science['original_evidence']['gate_manifest'] = proof(x['gate_path'])
            x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
            result = invoke('wrapper_completed_audit_failed'); assert result.returncode == 0, result.stderr
            x['gate_path'].write_bytes(gate_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            for field in ('current_labels_valid', 'sample_order_valid', 'timing_valid', 'artifact_contract_valid', 'worker_exit_not_established'):
                changed = json.loads(validation_saved); changed['rows'][0][field] = False
                x['validation'].write_text(json.dumps(changed)); science = copy.deepcopy(x['science']); science['scientific_validation'] = proof(x['validation']); x['receipt'].write_text(json.dumps(science))
                review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
                assert invoke('invalid_' + field).returncode != 0
            x['validation'].write_bytes(validation_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            # A resealed approval of a foreign ACTIVE owner never changes
            # the original-run ownership requirement.
            foreign_owner = Path(x['payloads'][0]['owner']['path']); owner_saved = foreign_owner.read_bytes()
            foreign_owner.write_bytes(owner_saved.replace(x['old'].name.encode(), b'foreign_run').replace(b'STATE=FAIL', b'STATE=ACTIVE'))
            science = copy.deepcopy(x['science']); science['payloads'][0]['owner'] = dict(**proof(foreign_owner), fields=dict(line.split('=', 1) for line in foreign_owner.read_text().splitlines()))
            x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
            rejected = invoke('foreign_active_owner'); assert rejected.returncode != 0 and 'Wrong artifact owner' in rejected.stderr, rejected.stderr
            foreign_owner.write_bytes(owner_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            # A forged but fully resealed audit/review cannot truncate one
            # corrected composition arm: the original contract expands all6.
            inv_saved = x['inv_path'].read_bytes(); rds_saved = x['rds_inputs'].read_bytes()
            removed = next(row['source'] for row in x['inventory'] if row['source'].endswith('_ecoda_unsupervised_harmony_limma.rds'))
            x['inv_path'].write_text(''.join(line + '\n' for line in inv_saved.decode().splitlines() if removed not in line))
            x['rds_inputs'].write_text(''.join(line + '\n' for line in rds_saved.decode().splitlines() if removed not in line))
            changed = json.loads(validation_saved); changed['rows'][0]['expected_artifacts'].remove(removed)
            changed['validator_runs'][0]['input_list'] = proof(x['rds_inputs'])
            x['validation'].write_text(json.dumps(changed)); science = copy.deepcopy(x['science'])
            science['inventory'] = proof(x['inv_path']); science['scientific_validation'] = proof(x['validation'])
            science['payloads'] = [p for p in science['payloads'] if p['source'] != removed]
            x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved)
            review.update(receipt_sha256=sha(x['receipt']), inventory_sha256=sha(x['inv_path']), payloads=len(science['payloads']))
            x['review'].write_text(json.dumps(review)); rejected = invoke('truncated_family')
            assert rejected.returncode != 0 and 'Incomplete per-row production artifact contract' in rejected.stderr, rejected.stderr
            x['inv_path'].write_bytes(inv_saved); x['rds_inputs'].write_bytes(rds_saved)
            x['validation'].write_bytes(validation_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            # Missing actual validator coverage fails even when all claimed
            # per-row booleans and expected artifact names remain complete.
            x['rds_inputs'].write_text(''.join(line + '\n' for line in rds_saved.decode().splitlines() if removed not in line))
            changed = json.loads(validation_saved); changed['validator_runs'][0]['input_list'] = proof(x['rds_inputs'])
            x['validation'].write_text(json.dumps(changed)); science = copy.deepcopy(x['science']); science['scientific_validation'] = proof(x['validation'])
            x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
            rejected = invoke('truncated_validator_inputs')
            assert rejected.returncode != 0 and 'Actual production validator input union differs' in rejected.stderr, rejected.stderr
            x['rds_inputs'].write_bytes(rds_saved); x['validation'].write_bytes(validation_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            for kind in ('rds', 'feather'):
                changed = json.loads(validation_saved); invocation = next(r for r in changed['validator_runs'] if r['kind'] == kind)
                argv = invocation['argv']; i = argv.index('--batch-pass'); del argv[i:i+2]
                x['validation'].write_text(json.dumps(changed)); science = copy.deepcopy(x['science']); science['scientific_validation'] = proof(x['validation'])
                x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
                rejected = invoke('missing_mode_' + kind)
                assert rejected.returncode != 0 and 'batch/variant mode differs' in rejected.stderr, rejected.stderr
            x['validation'].write_bytes(validation_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            for field in (1, 3, 4):
                lines = [line.split('\t') for line in rds_saved.decode().splitlines()]
                if field == 1:
                    row = next(r for r in lines if r[0].endswith('_ecoda_unsupervised_harmony.rds')); row[field] = 'ECODA_authors_HR'
                else:
                    lines[0][field] = 'benchmark_analysis' if field == 3 else '1'
                x['rds_inputs'].write_text(''.join('\t'.join(r) + '\n' for r in lines))
                changed = json.loads(validation_saved); changed['validator_runs'][0]['input_list'] = proof(x['rds_inputs'])
                x['validation'].write_text(json.dumps(changed)); science = copy.deepcopy(x['science']); science['scientific_validation'] = proof(x['validation'])
                x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
                rejected = invoke('wrong_rds_column_' + str(field))
                assert rejected.returncode != 0 and 'dataset/view/method/metadata columns differ' in rejected.stderr, rejected.stderr
            x['rds_inputs'].write_bytes(rds_saved); x['validation'].write_bytes(validation_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            changed = json.loads(validation_saved); changed['rows'][0].update(sample_ids=['s1', 's3'], labels=['a', 'b'], dropped_sample_ids=['s2'])
            x['validation'].write_text(json.dumps(changed)); science = copy.deepcopy(x['science']); science['scientific_validation'] = proof(x['validation'])
            x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
            rejected = invoke('non_scitd_drop'); assert rejected.returncode != 0 and 'only scITD may drop' in rejected.stderr, rejected.stderr
            x['validation'].write_bytes(validation_saved); x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            # Prep adoption is explicit immutable lineage, not a new fit.
            adopted = x['payloads'][0]['source']
            adoption = write(x['old'] / 'manifests/prepared_adoption.tsv', f'pseudobulk\tStephenson\t{x["view"]}\t/original/cache.rds\t{adopted}\tstage5_composition\t0\t1\n')
            receipt_lineage = write(x['old'] / 'manifests/prepared_adoption.receipt.json', '{"schema":"fixture-original-adoption"}\n')
            report_lineage = write(x['old'] / 'manifests/prepared_adoption.json', json.dumps(dict(operation='COPY_WITHOUT_REFIT', original_gate_state='FAILED', receipt_path=str(receipt_lineage), receipt_sha256=sha(receipt_lineage), binding_sha256=sha(x['binding_path']))))
            x['sealed'](receipt_lineage, 'stage5_prepared_adoption_receipt'); x['sealed'](report_lineage, 'stage5_prepared_adoption')
            rejected = invoke('missing_adoption_lineage'); assert rejected.returncode != 0
            science = copy.deepcopy(x['science']); science['adoption_lineage'] = [proof(p) for p in (adoption, receipt_lineage, report_lineage)]
            x['receipt'].write_text(json.dumps(science)); review = json.loads(review_saved); review['receipt_sha256'] = sha(x['receipt']); x['review'].write_text(json.dumps(review))
            assert invoke('adoption_lineage').returncode == 0
            x['receipt'].write_bytes(receipt_saved); x['review'].write_bytes(review_saved)
            for p in (adoption, receipt_lineage, report_lineage):
                p.unlink()
                Path(str(p) + '.md5').unlink(missing_ok=True)
                (x['old'] / 'manifests/artifacts' / (hashlib.sha256(str(p).encode()).hexdigest()[:32] + '.record')).unlink(missing_ok=True)
            saved = Path(x['payloads'][0]['source']).read_bytes(); Path(x['payloads'][0]['source']).write_bytes(saved + b'corrupt')
            assert invoke('payload_corrupt').returncode != 0; Path(x['payloads'][0]['source']).write_bytes(saved)
            x['missing_status'].parent.mkdir(parents=True, exist_ok=True); x['missing_status'].write_text('STATE=OK\n')
            assert invoke('fabricated_missing_status').returncode != 0; x['missing_status'].unlink()
            for p, before in immutable.items(): assert (sha(p), p.stat().st_mode) == before, p
            assert not (x['base'] / 'forbidden').exists()
        # An actual fourth-column CSV is one physical task with exactly the
        # original requested HR parameters, without the unrequested LR output.
        x = build(base / 'parameters', 'benchmark')
        tokens = ('hvg1000_highres', 'hvg2000_highres', 'hvg3000_highres')
        catalog = json.loads(x['catalog'].read_text()); catalog['rows'][1]['parameter'] = ','.join(tokens)
        matrix = Path(catalog['rows'][1]['manifest']['path'])
        matrix.write_text(f'Stephenson\tbenchmark_analysis\tqot\t{",".join(tokens)}\n'); x['sealed'](matrix, 'matrix_fixture')
        catalog['rows'][1]['manifest'] = proof(matrix); x['catalog'].write_text(json.dumps(catalog))
        requested = x['old'] / 'manifests/requested_selection.tsv'
        requested.write_text('Stephenson\tbenchmark_analysis\tscitd\n' + ''.join(f'Stephenson\tbenchmark_analysis\tqot\t{token}\n' for token in tokens))
        Path(x['contexts'][0].split('\t')[2]).write_bytes(requested.read_bytes())
        binding = json.loads(x['binding_path'].read_text()); binding['identities'][str(requested)] = sha(requested)
        x['binding_path'].write_text(json.dumps(binding))
        write(Path(str(x['binding_path']) + '.md5'), f'PATH={x["binding_path"]}\nMD5={hashlib.md5(x["binding_path"].read_bytes()).hexdigest()}\nSIZE={x["binding_path"].stat().st_size}\n')
        removed = next(r['source'] for r in x['inventory'] if 'hvg2000_lowres' in r['source'])
        inventory = [r for r in x['inventory'] if r['source'] != removed]
        keys = ('dataset', 'kind', 'source', 'relative', 'md5', 'size', 'original_producer')
        x['inv_path'].write_text(''.join('\t'.join(r[k] for k in keys) + '\n' for r in inventory))
        validation = json.loads(x['validation'].read_text()); validation['catalog_sha256'] = sha(x['catalog'])
        validation['rows'][1]['row'][3] = ','.join(tokens); validation['rows'][1]['expected_artifacts'].remove(removed)
        validation['validator_runs'] = [r for r in validation['validator_runs'] if removed not in r['argv']]
        x['validation'].write_text(json.dumps(validation)); science = copy.deepcopy(x['science'])
        science['catalog'] = proof(x['catalog']); science['inventory'] = proof(x['inv_path']); science['scientific_validation'] = proof(x['validation'])
        science['payloads'] = [p for p in science['payloads'] if p['source'] != removed]
        science['original_evidence'] = {k: proof(Path(p['path'])) for k, p in science['original_evidence'].items()}
        x['receipt'].write_text(json.dumps(science)); review = json.loads(x['review'].read_text())
        review.update(receipt_sha256=sha(x['receipt']), catalog_sha256=sha(x['catalog']), inventory_sha256=sha(x['inv_path']), payloads=len(inventory))
        x['review'].write_text(json.dumps(review)); result = x['invoke']('parameter_csv')
        assert result.returncode == 0, result.stdout + result.stderr
        report = json.loads((x['scratch'] / '_ecoda_runs/parameter_csv/artifacts/evidence/main_v2_science_publication.json').read_text())
        assert not any('lowres' in row['source'] for row in report['mappings']) and len(report['mappings']) == len(inventory)
        print('PASS: fixed main-v2 three-context copy-only publications; old FAIL/ACTIVE/status proofs immutable, current science/scope/accounting/owner/NAS guards')


if __name__ == '__main__':
    run()
