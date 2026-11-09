"""Fixed BL10 and Bassez HR22 releases; never read H5ADs or fit models.

Control schema ecoda_bl_publication_control_v1 contains run_id, release_id,
nas_root, science_receipt and review (the latter two are {path, sha256}).
The review is APPROVED_SCIENTIFIC_ARTIFACT_REUSE, binds receipt_sha256 and
declares datasets=10, rds_artifacts=20, timing_feathers=10,
original_gate_status=FAILED, publication_completed=false.
The second fixed schema is ecoda_bassez_hr22_publication_control_v1, with
the same keys and review fields but datasets=1, rds_artifacts=23, timings=1.
The fixed main-v2 arm adds context to those control keys. A reviewed finite
catalog, actual validator input coverage and independently expanded production
contracts qualify PUBLISHED artifacts. Complete failed terminal accounting is
preserved; physical worker exit is explicitly not established. Global batch
sample-metadata exports remain outside this selected-method publication.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import secrets
import shlex
import shutil
import subprocess
import sys

STAT_KEYS = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
ORIGINAL = 'stage5_benchmark_bl10_20261109_v2'
PRODUCER = 'stage5_bl_science_release_copy'
HR_SCHEMA = 'ecoda_bassez_hr22_publication_control_v1'
MAIN_SCHEMA = 'ecoda_main_v2_science_publication_control_v1'
REBOOT_SCHEMA = 'ecoda_main_v2_login_reboot_publication_control_v1'
MAIN_SCHEMAS = {MAIN_SCHEMA, REBOOT_SCHEMA}
MAIN_SOURCE = '9ca2453eea7320a721ef6dfd458f9db108abecac'
UNC4_SCHEMA = 'ecoda_failed16_unc4_publication_control_v1'
UNC4_SOURCE = '0314f207dfe17f70c2ffdfcc31c4ecd00379b41c'
UNC4_RUN = 'stage5_failed16_uncorrected_20261109_v1'
UNC4_DATASETS = ('Covid19_PBMC', 'Kidney_KPMP_full', 'Lupus_PBMC', 'Lung')
UNC4_SAMPLES = (44, 80, 94, 80)
MAIN_RUNS = dict(benchmark='stage5_main_benchmark_20261109_v2',
                 **{'batch_effect/uncorrected_final': 'stage5_main_uncorrected_20261109_v2',
                    'batch_effect/corrected_final': 'stage5_main_corrected_20261109_v2'})


def analysis_context(control):
    if control['schema'] == UNC4_SCHEMA:
        return 'batch_effect/uncorrected_final'
    return control['context'] if control['schema'] in MAIN_SCHEMAS else 'benchmark'


def campaign(control):
    if control['schema'] == UNC4_SCHEMA:
        return ('stage5_failed16_unc4_science_release_copy', 'failed16_unc4_science_publication.json',
                'ecoda_failed16_unc4_science_publication_receipt_v1', 'failed16_unc4_science_publication')
    if control['schema'] == REBOOT_SCHEMA:
        return ('stage5_main_v2_reboot_science_release_copy', 'main_v2_reboot_science_publication.json',
                'ecoda_main_v2_reboot_science_publication_receipt_v1', 'main_v2_reboot_science_publication')
    if control['schema'] in MAIN_SCHEMAS:
        return ('stage5_main_v2_science_release_copy', 'main_v2_science_publication.json',
                'ecoda_main_v2_science_publication_receipt_v1', 'main_v2_science_publication')
    if control['schema'] == HR_SCHEMA:
        return ('stage5_bassez_hr22_release_copy', 'bassez_hr22_publication.json',
                'ecoda_bassez_hr22_publication_receipt_v1', 'bassez_hr22_science_publication')
    require(control['schema'] == 'ecoda_bl_publication_control_v1', 'Unknown publication campaign')
    return (PRODUCER, 'bl_publication.json', 'ecoda_bl_publication_receipt_v1', 'bl_science_publication')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def safe(value, exists=True):
    require(isinstance(value, str) and value.startswith('/') and not any(c in value for c in '\n\r\t'), 'Unsafe absolute path')
    p = Path(value)
    require(str(p) == value and '..' not in p.parts and p.resolve() == p, 'Aliased/noncanonical path: ' + value)
    require(not any(x.is_symlink() for x in (p, *p.parents)), 'Symlink path: ' + value)
    if exists:
        require(p.exists(), 'Missing evidence: ' + value)
    return p


def digest(p, algorithm='sha256'):
    p = safe(str(p)); before = p.stat(); h = hashlib.new(algorithm)
    require(p.is_file(), 'Not a regular evidence file')
    require(p.suffix.lower() not in {'.h5ad', '.h5'}, 'Publication must not read expression inputs')
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    after = p.stat()
    require(all(getattr(before, k) == getattr(after, k) for k in STAT_KEYS), 'Evidence changed during read')
    return h.hexdigest()


def unique(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON/proof key')
        result[key] = value
    return result


def read_json(p):
    return json.loads(safe(str(p)).read_text(), object_pairs_hook=unique)


def fields(p):
    lines = safe(str(p)).read_text().splitlines()
    require(lines and all('=' in line for line in lines), 'Malformed proof fields')
    return unique(line.split('=', 1) for line in lines)


def proof(item):
    p = safe(item['path'])
    require(re.fullmatch('[a-f0-9]{64}', item['sha256']) and digest(p) == item['sha256'], 'Proof digest differs: ' + str(p))
    return p


def record(p, run, producer):
    p = safe(str(p)); md5 = digest(p, 'md5'); size = str(p.stat().st_size)
    expected = dict(PATH=str(p), MD5=md5, SIZE=size)
    require(fields(str(p) + '.md5') == expected and int(size) > 0, 'Strict sidecar differs')
    key = hashlib.sha256(str(p).encode()).hexdigest()[:32]
    path = run / 'manifests/artifacts' / (key + '.record')
    require(fields(path) == dict(**expected, RUN_ID=run.name, PRODUCER=producer, STATE='PUBLISHED'), 'Producer record differs')
    return dict(path=str(p), md5=md5, size=size, record_path=str(path), record_sha256=digest(path))


def owner_file(p, run):
    key = re.sub(r'[/,:\t |]', '_', p.name) + '_' + hashlib.sha256(str(p).encode()).hexdigest()[:32]
    if len(key.encode()) > 100:
        key = 'artifact_' + hashlib.sha256(str(p).encode()).hexdigest()[:32]
    return run.parent.parent / '_ecoda_owners/artifact' / key / 'owner'


def owner(p, run, states):
    path = owner_file(p, run)
    value = fields(path)
    require(value.get('RUN_ID') == run.name and value.get('STAGE') == 'stage5' and value.get('PATH') == str(p)
            and value.get('STATE') in states, 'Wrong artifact owner')
    return dict(path=str(path), sha256=digest(path), fields=value)


def main_v2_expected_artifacts(original, source, config, binding, context, row):
    """Expand existing production contracts; never discover files or fit work."""
    ds, view, method, parameter = (row[k] for k in ('dataset', 'view', 'method', 'parameter'))
    root = original / 'artifacts' / context
    if parameter:
        require(context == 'benchmark', 'Batch parameter publication is unsupported')
        path = Path(source['SOURCE_ROOT']) / 'src/utils/py/benchmark_request_scope.py'
        namespace = {'__name__': '_main_v2_original_requests'}
        exec(compile(path.read_bytes(), str(path), 'exec'), namespace)
        tokens = parameter.split(',')
        require(len(tokens) == len(set(tokens)) and all(token in namespace['TOKENS'].get(method, ()) for token in tokens), 'Unknown original parameter scope')
        paths = {str(root / relative) for relative in namespace['selected_files'](ds, method, tuple(tokens))}
        return paths, {(path, method, ds, view, str(int(path.endswith('_metadata.rds')))) for path in paths if path.endswith('.rds')}
    env = dict(os.environ, PROJECT_ROOT=source['SOURCE_ROOT'], DATASETS_JSON_FILE=str(config), PYTHONDONTWRITEBYTECODE='1',
               PYTHON_BIN=sys.executable, HPC_SCRATCH_DIR=str(original.parent.parent),
               NAS_TARGET_DIR=binding['nas_base'], ECODA_RUN_ROOT=str(original), ECODA_RUN_ID=original.name,
               ECODA_ARTIFACT_BINDING=str(original / 'manifests/artifact_bindings.json'),
               ANALYSIS_ROOT=str(root), ANALYSIS_NAS_ROOT=str(Path(binding['nas_root']) / context),
               PASS_ARG='' if context == 'benchmark' else 'corrected' if context.endswith('/corrected_final') else 'uncorrected',
               ANALYSIS_VARIANT='' if context == 'benchmark' else 'corrected_final' if context.endswith('/corrected_final') else 'final')
    env['ANALYSIS_PASS'] = env['PASS_ARG']
    for name in ('METHOD_MATRIX', 'ECODA_STAGE5_METHOD_MATRIX', 'ECODA_STAGE5_METHOD_MATRIX_MODE'):
        env.pop(name, None)
    command = ['bash', '-c', 'set -euo pipefail; source "$1/src/utils/bash/ecoda_run_common.sh"; '
               'source "$1/src/utils/bash/ecoda_stage5_policy.sh"; '
               '_ecoda_stage5_artifacts_for "$2" "$3" "$4" || exit 1; '
               'for path in "${ECODA_BENCHMARK_ARTIFACTS[@]}"; do '
               'printf "%s\\t%s\\n" "$path" "$(ecoda_stage5_rds_artifact_method "$4" "$path")"; done',
               'publication-contract', source['SOURCE_ROOT'], ds, view, method]
    result = subprocess.run(command, env=env, capture_output=True, text=True)
    require(result.returncode == 0 and result.stdout.strip(), 'Original artifact contract expansion failed: ' + result.stderr)
    paths = [line.split('\t') for line in result.stdout.splitlines()]
    require(all(len(item) == 2 for item in paths), 'Malformed original artifact contract expansion')
    return {item[0] for item in paths}, {(path, label, ds, view, str(int(path.endswith('_metadata.rds'))))
                                      for path, label in paths if path.endswith('.rds')}


def cli_value(argv, name):
    require(isinstance(argv, list) and all(isinstance(word, str) for word in argv)
            and not any(word.startswith(name + '=') for word in argv), 'Ambiguous validator option ' + name)
    indices = [i for i, word in enumerate(argv) if word == name]
    require(len(indices) <= 1 and (not indices or indices[0] + 1 < len(argv)), 'Ambiguous validator option ' + name)
    return argv[indices[0] + 1] if indices else None


def verified_scope(args):
    control_path = safe(args.control)
    require(digest(control_path, 'md5') == args.control_md5, 'Control MD5 differs')
    control = read_json(control_path)
    expected_keys = {'schema', 'run_id', 'release_id', 'nas_root', 'science_receipt', 'review'}
    if control.get('schema') in MAIN_SCHEMAS:
        expected_keys.add('context')
    if control.get('schema') == REBOOT_SCHEMA:
        expected_keys.add('incident')
    require(set(control) == expected_keys
            and control['schema'] in {'ecoda_bl_publication_control_v1', HR_SCHEMA, UNC4_SCHEMA, *MAIN_SCHEMAS}, 'Unknown publication control')
    run = safe(args.run_root); nas = safe(args.nas_root)
    require(run.is_dir() and run.parent.name == '_ecoda_runs' and control['run_id'] == run.name
            and control['nas_root'] == str(nas) and nas.is_dir()
            and re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*', control['release_id']), 'Publication roots differ')
    require(control_path == run / 'manifests/control.json', 'Publication control must be run-owned')
    receipt_path = proof(control['science_receipt']); review = read_json(proof(control['review'])); science = read_json(receipt_path)
    if control['schema'] == UNC4_SCHEMA:
        return verified_unc4_scope(control, science, review, run, nas)
    if control['schema'] in MAIN_SCHEMAS:
        return verified_main_v2_scope(control, science, review, run, nas)
    if control['schema'] == HR_SCHEMA:
        return verified_hr22_scope(control, science, review, receipt_path, run, nas)
    require(review['status'] == 'APPROVED_SCIENTIFIC_ARTIFACT_REUSE'
            and review['receipt_sha256'] == control['science_receipt']['sha256']
            and review['datasets'] == 10 and review['rds_artifacts'] == 20 and review['timing_feathers'] == 10
            and review['original_gate_status'] == 'FAILED' and review['publication_completed'] is False, 'Scientific review does not approve this exact receipt')
    require(science['schema'] == 'ecoda_bl_science_receipt_v1' and science['status'] == 'NOOP_VALIDATED'
            and science['science_validated'] is True and science['publication_completed'] is False
            and science['original_gate_status'] == 'FAILED', 'Science receipt is not the reviewed failed-gate validation')
    original = safe(science['original_run_root']); validation = safe(science['audit_root'])
    require(original.name == ORIGINAL and original.parent == run.parent and validation.parent == run.parent
            and run not in (original, validation) and receipt_path == validation / 'artifacts/evidence/bl_science.json', 'Wrong original/validation run')
    record(receipt_path, validation, 'bl_science_validation'); owner(receipt_path, validation, {'OK'})
    for item in science['original_evidence'].values():
        proof(item)
    proof(science['validation_source']); proof(science['validation_runtime']); proof(science['production_validator'])
    raw = science['original_accounting']; rows = raw['rows']
    require(raw['query_count'] == 1 and raw['worker_array_id'] == '4663563' and len(rows) == 23
            and len({row['job_id'] for row in rows}) == 23, 'Original accounting scope differs')
    require(proof(dict(path=raw['raw_path'], sha256=raw['sha256'])).read_bytes() == raw['raw_stdout'].encode(), 'Original accounting bytes differ')
    parsed = [line.split('|') for line in raw['raw_stdout'].splitlines()]
    require(raw['raw_stdout'].endswith('\n') and len(parsed) == 23
            and rows == [dict(line=i, job_id=r[0], state=r[1], exit_code=r[2]) for i, r in enumerate(parsed, 1)]
            and all(re.fullmatch('[1-9][0-9]*', r[0]) and r[1:] == (['FAILED', '1:0'] if r[0] in {'4663573', '4663574'} else ['COMPLETED', '0:0']) for r in parsed)
            and {'4663573', '4663574'} <= {r[0] for r in parsed}, 'All23 original outcomes must remain visible')
    gate = read_json(proof(science['original_evidence']['gate_manifest']))
    inspected = read_json(proof(science['original_evidence']['inspect']))
    scheduler_ids = gate['scheduler_array_ids'] + gate['scheduler_watchdog_ids']
    require(gate['state'] == inspected['state'] == 'FAILED' and gate['audit'] == inspected['audit']
            and inspected['audit']['accounting']['stdout'] == raw['raw_stdout']
            and inspected['audit']['accounting']['rows'] == rows
            and scheduler_ids == raw['scheduler_ids'] and len(set(scheduler_ids)) == 5
            and set(scheduler_ids) <= {r[0] for r in parsed}, 'Original gate evidence changed')
    for failure in raw['orchestration_failures']:
        proof(failure['stdout']); proof(failure['stderr'])
    for item in science['input_proofs'].values():
        proof(dict(path=item['record_path'], sha256=item['record_sha256']))
        proof(dict(path=item['reference'][5], sha256=item['external_proof_sha256'])); proof(item['status'])
    for worker in science['worker_evidence']:
        proof(worker['stdout']); proof(worker['stderr'])
    dataset_order = [r.split('\t')[0] for r in proof(science['original_evidence']['worker_manifest']).read_text().splitlines()]
    require(len(dataset_order) == len(set(dataset_order)) == 10
            and set(dataset_order) == set(science['artifacts']) == set(science['timings']) == set(science['input_proofs']), 'Exact physical BL10 coverage differs')
    inventory = []
    for ds in dataset_order:
        require(re.fullmatch('[A-Za-z0-9_]+', ds), 'Unsafe dataset ID')
        artifacts = science['artifacts'][ds]
        specifications = [('cache', 'stage5_prepare_pseudobulk_hvg2000_bl', f'pseudobulks/{ds}_pseudobulk_hvg2000_bl.rds'),
                          ('result', 'stage5_pseudobulk', f'results/{ds}_Pseudobulk_hvg2000_bl.rds'),
                          ('timing', 'stage5_execution_log', f'embeddings/timing_sources/execution_times_pseudobulk_{ds}_hvg2000_bl.feather')]
        for kind, producer, relative in specifications:
            item = science['timings'][ds]['proof'] if kind == 'timing' else artifacts[kind]
            path = safe(item['path'])
            expected = original / ('logs/' + Path(relative).name if kind == 'timing' else 'artifacts/benchmark/' + relative)
            require(path == expected, 'Original payload path differs')
            current = record(path, original, producer)
            require(all(current[k] == item[k] for k in current) and fields(item['record_path']) == item['record']
                    and digest(str(path) + '.md5') == item['sidecar_sha256'], 'Approved payload/proof bytes changed')
            if kind != 'timing':
                previous_owner = artifacts['owners'][0 if kind == 'cache' else 1]
                require(owner(path, original, {'OK', 'FAIL', 'ACTIVE'}) == {k: previous_owner[k] for k in ('path', 'sha256', 'fields')}, 'Original owner changed')
            inventory.append(dict(dataset=ds, kind=kind, source=str(path), relative=relative, md5=item['md5'], size=item['size'], original_producer=producer))
    destination = safe(str(nas / 'releases' / control['release_id'] / run.name / 'benchmark'), False)
    return control, science, inventory, run, nas, destination


def main_v2_writer_tokens(gate):
    """Exact old entrypoints/controls also identify env-only run-ID controllers."""
    argv = shlex.split(gate['exact_command'])
    source = safe(cli_value(argv, '--source-root'))
    contexts = safe(cli_value(argv, '--contexts'))
    require(contexts.suffix == '.tsv' and contexts.stat().st_size <= 65536, 'Wrong original contexts control')
    before = digest(contexts)
    rows = [line.split('\t') for line in contexts.read_text().splitlines()]
    children = {'benchmark': MAIN_RUNS['benchmark'], 'uncorrected': MAIN_RUNS['batch_effect/uncorrected_final'],
                'corrected': MAIN_RUNS['batch_effect/corrected_final']}
    require(len(rows) == 3 and {row[0] for row in rows} == set(children)
            and all(len(row) == 4 and children.get(row[0]) == row[1] for row in rows)
            and digest(contexts) == before, 'Wrong/changing original context identities')
    return [str(source / 'src/5_run_benchmark_methods' / script) for script in ('stage5_campaign.sh', 'stage5_dispatcher.sh')] + [
        str(contexts), *(str(safe(path)) for row in rows for path in row[2:])]


REBOOT_RETRY_SOURCES = {
    'src/5_run_benchmark_methods/matrix_watchdog.sh': 'a58132d2e63aaddd3622a2b64b6f0ecceb061ae1da2866dba038815712e68ccf',
    'src/5_run_benchmark_methods/stage5_dispatcher.sh': 'f1c243f80a0c578eba67a5b6b88976d04f113b332feb34010ede27f8354b05a3',
    'src/5_run_benchmark_methods/benchmark_submit_common.sh': '1729ce35a5002fe85d9ed670fd12addf7e9e41ec7ea62a5270133b0841280844',
    'src/5_run_benchmark_methods/stage5_campaign.sh': 'b85f11a1f9d976b5263a4a3beb9947bddf9ba4a1751ae7151b2c6197be740cca',
    'src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh': 'a2591ef9fa889df0685f2103e1133737d15a4da7cd95b9d5b5f1e30279fb23e2',
    'src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh': '08b13e30d2fc5928860007b4fed4492cc58939354ca983e2d4bfde5b142da17f',
    'src/5_run_benchmark_methods/matrix_gate.sh': '5130fba6a99e6b22cd6cebd0db125372418a2741149448208a01f7f7bbf48cfd',
    'src/utils/bash/ecoda_run_common.sh': '0ad37f4feaeca5d561d47da60bfd85db38c3bcb9d33c6f02d9284b61ec844cb1',
    'src/utils/bash/worker_retry.sh': '404ba975dafa7ecf67302715fc072b398b79457746aea223563f53594782fb0e',
    'src/utils/bash/ecoda_runtime.sh': '540fb2fd2f9462053ba81382c5890a68da256da3975f392557192f6fb81b2cdd',
    'src/utils/bash/ecoda_stage5_policy.sh': '8fe00e36ef8a1649a0e37720f89ecc8046f906039615627b574f5e48d7a39e6f',
    'src/utils/bash/h5ad_preflight_worker.sh': '352a601a435136101b07f4161111543885cb89c10851a5e4d0f25249796c147f',
    'src/utils/bash/h5ad_obs_audit_worker.sh': 'a109ab63b3b4eea5d77aa15823d5ccf557e661c5ce5db0064021fb00b8eaf4e3',
    'src/slurm_config.sh': 'c96e7fb8422801f2a76cd2b2df007b81f69bbc6fb0592750a9ceb7ceff7b1a9a',
    'src/utils/bash/r_environment_preflight_worker.sh': '3ade9d24036b6b74c04d3144b958dc0143b1e699f8072dcd5b42208f3604a1ec',
    'src/utils/bash/h5ad_preflight_submit.sh': 'e3dd16ebf2e7ddf10d3400d652808cee50c988c1c6c2f3b9fcfab5d67be0afe3',
    'src/5_run_benchmark_methods/stage5_pending.sh': '9fe6b1403f874ff2d8cf6d396aa8cc5264e0dd3cf40074853d9f8f777060595c',
    'src/5_run_benchmark_methods/stage5_validated_recovery.sh': '5984b7c6c1dfbb2295aaac93cd5ac88fcf5b45c40dbb90dafc064e682f44ce6b',
    'src/5_run_benchmark_methods/stage5_prepared_adoption.sh': 'b580a272165e1d1535bb2c7689fef37d6b076f9a4a6d2592ac2d221c76cef57b',
    'src/utils/bash/ecoda_derivative_artifacts.sh': '2a621812cabbe8fee9fa084fc1cc73725bce767685552ae84fc11bfac7c659ba',
}


def verify_no_retry_coverage(incident, paths, initial, original):
    """Fixed9ca alternative: a retry cannot be submitted before its preserved manifest."""
    coverage = read_json(proof(incident['no_retry_coverage']))
    observed_path = proof(coverage['observation']); observed = read_json(observed_path)
    roots = {str(original.parent / child) for child in MAIN_RUNS.values()}
    initial_sha = digest(paths['initial_catalog'])
    require(coverage['schema'] == 'ecoda_main_v2_no_retry_coverage_v1'
            and coverage['status'] == 'APPROVED_NO_NEW_RETRY_IDS'
            and coverage['source_commit'] == MAIN_SOURCE and coverage['initial_catalog_sha256'] == initial_sha
            and observed['schema'] == 'ecoda_old_main_catalog32_observation_v1'
            and observed['status'] == 'OBSERVED_NOT_APPROVED' and observed['initial_ids'] == 73
            and observed['collector_sha256'] == '9867940459b07962504c35d66a6fcd82aad2de8836f4248f6c79509fab2387f8'
            and observed['control_sha256'] == 'ef5d687d045890b31d769bee25bd04de09d4919726bf4c38b8409712dbec9a6a'
            and observed['stable'] is True and observed['final_emitted_catalog_complete'] is False
            and observed['scientific_reuse_approved'] is False
            and all(observed[k] is True for k in ('no_accounting_query', 'no_queue_query', 'no_binary_payload_reads', 'no_remote_writes')),
            'No-retry coverage is not its reviewed bounded observation')
    source_root = Path(fields(original / 'manifests/source.manifest')['SOURCE_ROOT'])
    for relative, expected in REBOOT_RETRY_SOURCES.items():
        require(digest(source_root / relative) == expected, 'Exact9ca retry source changed: ' + relative)
    source_digest = hashlib.sha256(json.dumps(REBOOT_RETRY_SOURCES, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    preservation = read_json(proof(coverage['preservation_review']))
    writers = read_json(proof(coverage['writer_review']))
    require(preservation['schema'] == 'ecoda_main_v2_retry_manifest_preservation_review_v1'
            and preservation['status'] == 'APPROVED_SOURCE_PRECONDITION_NO_RETRY_COVERAGE'
            and preservation['source_pins_sha256'] == source_digest
            and preservation['source_preservation_confirmed'] is True
            and preservation['no_retry_manifest_deletion_or_controller_reexecution'] is True,
            'Missing exact source/persistent-manifest preservation approval')
    require(writers['schema'] == 'ecoda_main_v2_scoped_quiescence_review_v1'
            and writers['status'] == 'APPROVED_NO_OLD_ROOT_WRITERS'
            and writers['quiescence_observation_sha256'] == digest(paths['quiescence_observation'])
            and datetime.datetime.fromisoformat(read_json(paths['quiescence_observation'])['observed_at']).timestamp() >= observed['observed_at_epoch']
            and writers['no_old_root_writers'] is True and writers['no_unknown_writers'] is True
            and writers['allowed_disjoint_source_commit'] == '32769f1cde54346e8b003ad096ca9e5b7490a517'
            and set(writers['allowed_disjoint_run_roots']) == {str(original.parent / ('stage5_remaining95_' + part + '_20261109_v1')) for part in ('benchmark', 'uncorrected', 'corrected')} | {str(original.parent / 'stage5_remaining95_20261109_v1')}
            and writers['evidence'] and all(proof(item) for item in writers['evidence']),
            'Missing fresh scoped no-old-writer approval')
    for review in (preservation, writers):
        require(review['original_source_commit'] == MAIN_SOURCE and set(review['original_run_roots']) == roots
                and review['observation_sha256'] == digest(observed_path)
                and review['initial_catalog_sha256'] == initial_sha
                and review['original_gate_state'] == 'PRELAUNCH_STOP' and review['wrapper_exit_code'] is None
                and type(review['reviewed_at_epoch']) in (float, int) and review['reviewed_at_epoch'] >= observed['observed_at_epoch'],
                'No-retry review scope/freshness/unknown wrapper differs')
    files = {row['path']: row for row in observed['files']}
    require(len(files) == len(observed['files']), 'Duplicate no-retry snapshot')
    expected_paths = set()
    def unchanged(path):
        path = str(path); expected_paths.add(path)
        require(path in files and type(files[path]['exists']) is bool, 'Missing explicit no-retry path snapshot')
        item = files[path]; current = safe(path, item['exists'])
        if item['exists']:
            proof(item)
            require(set(item['stat']) == set(STAT_KEYS) | {'st_mode'}
                    and current.is_file() and current.stat().st_size == item['size']
                    and hashlib.sha256(item['text'].encode()).hexdigest() == item['sha256']
                    and all(getattr(current.stat(), key) == value for key, value in item['stat'].items()),
                    'No-retry observed text/stat changed')
        else:
            require(not current.exists(), 'An absent no-retry path appeared')
        return item
    pairs = observed['groups']; by_watchdog = {}
    require(len(pairs) == 32 and len({g['pair']['watchdog_id'] for g in pairs}) == 32, 'Incomplete32 native/precondition coverage')
    for root in roots:
        for name in ('source.manifest', 'runtime.identity', 'scheduler_ids.tsv', 'scheduler_ids.tsv.md5'):
            unchanged(Path(root) / 'manifests' / name)
    parent = original.parent / 'stage5_main_confounds_20261109_v2'
    controller_pairs = {}
    for short in ('benchmark', 'uncorrected', 'corrected'):
        log = unchanged(parent / 'logs' / (short + '.log'))
        require(log['exists'], 'Missing original controller pairing evidence')
        arrays = re.findall(r'(?:BENCHMARK|BATCH_EFFECT)_ARRAY_JOB_ID=([^\s:]+):([1-9][0-9]*)', log['text'])
        watches = re.findall(r'(?:BENCHMARK|BATCH_EFFECT)_WATCHDOG_JOB_ID=([^\s:]+):([1-9][0-9]*)', log['text'])
        require(len(dict(arrays)) == len(arrays) and len(dict(watches)) == len(watches)
                and dict(arrays).keys() == dict(watches).keys(), 'Ambiguous original controller pairings')
        for label, wd in watches:
            require(wd not in controller_pairs, 'Duplicate original watchdog pairing')
            controller_pairs[wd] = (short, label, dict(arrays)[label])
    for group in pairs:
        pair = group['pair']; child = Path(pair['run_root']); label = pair['label']; wd = pair['watchdog_id']
        require(str(child) in roots and re.fullmatch('[A-Za-z0-9_]+', label)
                and MAIN_RUNS.get(pair['context']) == child.name, 'Foreign no-retry root')
        short = next(k for k,v in {'benchmark': MAIN_RUNS['benchmark'], 'uncorrected': MAIN_RUNS['batch_effect/uncorrected_final'], 'corrected': MAIN_RUNS['batch_effect/corrected_final']}.items() if v == child.name)
        require(controller_pairs.get(wd) == (short, label, pair['array_id'])
                and [short, child.name, 'WATCHDOG', wd] in initial and [short, child.name, 'ARRAY', pair['array_id']] in initial
                and pair['matrix'] == str(child / 'manifests' / ('matrix_' + label.replace('__', '_') + '.tsv')),
                'No-retry pair differs from original ledger/matrix')
        matrix = unchanged(pair['matrix']); unchanged(pair['matrix'] + '.md5')
        require(matrix['exists'], 'Original matrix disappeared')
        snapshots = {}
        for kind, suffix in (('status', ''), ('stdout', '.log'), ('stderr', '.err')):
            path = child / ('status/watchdogs/' + label + '.status' if kind == 'status' else 'logs/5_matrix_watchdog_' + label + '_' + wd + suffix)
            require(group['files'][kind] == str(path), 'Foreign watchdog evidence path')
            snapshots[kind] = unchanged(path)
        status = snapshots['status']
        if status['exists']:
            lines = status['text'].splitlines(); ids = [line[13:] for line in lines if line.startswith('SCHEDULER_ID=')]
            ordinary = [line for line in lines if not line.startswith('SCHEDULER_ID=')]
            values = unique(line.split('=', 1) for line in ordinary)
            require(set(values) == {'STATE','LABEL','REASON','ARRAY_JOB_ID','RETRY_INDEX'}
                    and values['STATE'] in {'OK','FAIL'} and values['LABEL'] == label
                    and values['ARRAY_JOB_ID'] == pair['array_id'] and values['RETRY_INDEX'] == '0' and ids == [pair['array_id']],
                    'Native watchdog contradicts no-new-retry coverage')
        for item in snapshots.values():
            require(not item['exists'] or not re.search(r'^(?:MATRIX|BATCH_EFFECT)_RETRY_ARRAY_JOB_ID=', item['text'], re.M), 'Recorded retry contradicts coverage')
        for index in range(1,5):
            for suffix in ('', '_corrected_final'):
                path = child / 'manifests' / (label + suffix + '.retry_' + str(index) + '.tsv')
                require(not unchanged(path)['exists'] and not unchanged(str(path) + '.md5')['exists'], 'A retry pre-submission manifest/proof exists')
        by_watchdog[wd] = snapshots
    require(set(controller_pairs) == set(by_watchdog) and len(controller_pairs) == 32, 'Incomplete original controller pairing coverage')
    require(set(files) == expected_paths, 'Unexpected/missing finite no-retry metadata paths')
    directories = observed['directories']
    require(len(directories) == 9 and {row['path'] for row in directories} == {str(Path(root)/name) for root in roots for name in ('logs','manifests','status/watchdogs')}, 'Wrong bounded directory scope')
    for row in directories:
        directory = safe(row['path']); entries = []
        for entry in directory.iterdir():
            st = entry.lstat()
            require('retry' not in entry.name.lower(), 'Unexpected retry filename contradicts coverage')
            if directory.name == 'watchdogs' or entry.name.startswith('5_matrix_watchdog_'):
                require(str(entry) in expected_paths, 'Unclassified native watchdog evidence')
            entries.append(dict(name=entry.name, mode=st.st_mode, size=st.st_size, symlink=entry.is_symlink()))
            require(len(entries) <= 10000, 'Unbounded original metadata directory')
        require(row['exists'] is True and sorted(entries, key=lambda item:item['name']) == row['entries'], 'No-retry directory listing changed')
    coverage['_watchdogs'] = by_watchdog
    return coverage


def verify_ownerless_gloscope(row, item, context, original, checked, coverage):
    require(coverage is not None and 'ownerless_gloscope_review' in coverage, 'Missing explicit ownerless GloScope review')
    review = read_json(proof(coverage['ownerless_gloscope_review']))
    names = {'Stephenson_GloScope_' + token + '.rds' for token in ('hvg2000_pcadims10','hvg2000_pcadims50','hvg1000_pcadims30','hvg3000_pcadims30')}
    paths = {str(original / 'artifacts/benchmark/results' / name) for name in names}
    family_key = 'derivative/' + original.name + '/ordinary/Stephenson/benchmark_analysis/gloscope'
    family_path = original.parent.parent / '_ecoda_owners/stage5' / re.sub(r'[/,:\t |]', '_', family_key) / 'owner'
    family = fields(family_path)
    require(family['RUN_ID'] == original.name and family['STAGE'] == 'stage5'
            and family['KEY'] == family_key and family['STATE'] in {'OK', 'FAIL', 'ACTIVE'}, 'Wrong original GloScope family owner')
    safe(str(owner_file(Path(row['source']), original)), False)
    require(context == 'benchmark' and original.name == MAIN_RUNS['benchmark']
            and review['schema'] == 'ecoda_main_v2_ownerless_gloscope_review_v1'
            and review['status'] == 'APPROVED_EXACT_GLOSCOPE_OWNERLESS_PRODUCER_DESIGN'
            and review['original_source_commit'] == MAIN_SOURCE and review['original_run_root'] == str(original)
            and set(review['paths']) == paths and review['producer'] == 'stage5_gloscope'
            and review['no_original_owner_created'] is True and review['writer_review_sha256'] == coverage['writer_review']['sha256']
            and row['source'] in paths and row['kind'] == 'result' and row['original_producer'] == 'stage5_gloscope'
            and item['owner'] is None and not owner_file(Path(row['source']), original).exists()
            and not owner_file(Path(row['source']), original).is_symlink()
            and checked['row'] == ['Stephenson','benchmark_analysis','gloscope',Path(row['source']).stem.removeprefix('Stephenson_GloScope_')]
            and len(checked['sample_ids']) == 53 and checked['dropped_sample_ids'] == []
            and all(checked[k] is True for k in ('artifact_contract_valid','current_labels_valid','sample_order_valid','timing_valid')),
            'Ownerless permission is limited to four reviewed53-sample GloScope shards')


def verify_login_reboot(control, science, review, gate, inspected, evidence, schedulers, ids, original):
    """One fixed reboot incident; never turn missing wrapper status into success."""
    incident = read_json(proof(control['incident']))
    require(set(incident) - {'no_retry_coverage'} == {'schema', 'orchestration_outcome', 'wrapper_exit_code', 'observation',
                             'aggregate_evidence', 'initial_catalog', 'initial_inventory', 'emitted_catalog',
                             'emitted_inventory', 'quiescence_observation'}
            and incident['schema'] == 'ecoda_main_v2_login_reboot_incident_v1'
            and incident['orchestration_outcome'] == 'LOGIN_REBOOT_LOSS'
            and incident['wrapper_exit_code'] is None, 'Wrong fixed reboot incident')
    for value in (science, review):
        require(value['incident_sha256'] == control['incident']['sha256']
                and value['orchestration_outcome'] == 'LOGIN_REBOOT_LOSS'
                and value['wrapper_exit_code'] is None, 'Reboot scientific review/receipt differs')
    require(gate.get('terminal_status') is None and gate['remote_exit_code'] is None
            and gate['remote_host'] == 'bamboo' and gate['launch_timestamp'] == '2026-11-09T08:18:46Z',
            'Reboot must preserve the absent actual wrapper outcome')
    paths = {key: proof(incident[key]) for key in ('observation', 'aggregate_evidence', 'initial_catalog',
             'initial_inventory', 'emitted_catalog', 'emitted_inventory', 'quiescence_observation')}
    require(paths['emitted_catalog'] == evidence['campaign_scheduler_ids']
            and original.parent / 'stage5_main_confounds_20261109_v2' not in paths['emitted_catalog'].parents,
            'Reboot inventory must be additive, outside original parent manifests')
    for key, suffix in (('remote_manifest', '.manifest.json'), ('remote_runner', '.runner.sh'),
                        ('remote_status', '.status.json'), ('remote_failure_status', '.status.json.failure')):
        require(gate[key] == str(original.parent.parent / 'gates' / (gate['gate_id'] + suffix)), 'Foreign original gate path')
    observation = read_json(paths['observation']); quiescence = read_json(paths['quiescence_observation'])
    boot = datetime.datetime.fromisoformat('2026-11-09T14:47:27+00:00')
    writer_tokens = [*MAIN_RUNS.values(), 'stage5_main_confounds_20261109_v2', *main_v2_writer_tokens(gate)]
    for item in (observation, quiescence):
        require(item['schema'] == 'ecoda_login_reboot_operational_observation_v1'
                and item['hostname'] == 'login1.bamboo' and item['boot_epoch'] == 1794235647
                and datetime.datetime.fromisoformat(item['boot_utc']) == boot
                and datetime.datetime.fromisoformat(item['observed_at']) > boot
                and all(item[key] is True for key in ('no_accounting_query', 'no_scientific_payload_reads', 'no_remote_writes'))
                and item['tmux']['exit'] in {0, 1}
                and (item['tmux']['exit'] == 0 or not item['tmux']['stdout'])
                and gate['tmux_session'] not in item['tmux']['stdout']
                and item['processes']['exit'] == 0 and item['queue']['exit'] == 0,
                'Missing reboot/host/absence observation')
        processes = item['processes']['stdout']
        require(not any(name in processes for name in writer_tokens),
                'An original controller/runner is still observed')
        files = {row['path']: row for row in item['files']}
        require(len(files) == len(item['files']), 'Duplicate reboot snapshot')
        for path in (gate['remote_status'], gate['remote_failure_status'],
                     str(Path(gate['remote_status']).with_name(gate['gate_id'] + '.done'))):
            require(path in files and files[path]['exists'] is False and not safe(path, False).exists(),
                    'Original terminal/event evidence was created or was not observed absent')
    require(observation['tmux']['exit'] == 1 and not observation['tmux']['stdout'], 'Original incident did not observe tmux loss')
    require(datetime.datetime.fromisoformat(quiescence['observed_at']) >= datetime.datetime.fromisoformat(observation['observed_at']),
            'Stale quiescence observation')
    # Preserve actual pre-loss command/generation evidence, not an invented exit.
    observed = {row['path']: row for row in observation['files']}
    for key in ('remote_manifest', 'remote_runner'):
        row = observed[gate[key]]
        proof(row)
        require(row['exists'] is True and row['truncated'] is False
                and hashlib.sha256(row['text'].encode()).hexdigest() == row['sha256'], 'Truncated/invalid original runner snapshot')
    launch = json.loads(observed[gate['remote_manifest']]['text'], object_pairs_hook=unique)
    require(all(launch[key] == gate[key] for key in ('gate_id', 'exact_command', 'command_digest', 'event_generation', 'remote_host'))
            and gate['exact_command'] in observed[gate['remote_runner']]['text']
            and gate['event_generation'] in observed[gate['remote_runner']]['text'], 'Reboot launch identity differs')
    lines = lambda path: [line.split('\t') for line in path.read_text().splitlines()]
    initial = lines(paths['initial_catalog']); initial_proof = read_json(paths['initial_inventory'])
    final = read_json(paths['emitted_inventory'])
    require(initial_proof['schema'] == 'ecoda_orphaned_main_initial_scheduler_inventory_v1'
            and initial_proof['catalog_sha256'] == digest(paths['initial_catalog'])
            and initial_proof['initial_id_count'] == len(initial) == 73
            and initial_proof['final_emitted_catalog_complete'] is False
            and initial_proof['original_files_written'] is False
            and initial_proof['context_counts'] == {'benchmark': 39, 'uncorrected': 17, 'corrected': 17}
            and {key: sum(row[0] == key for row in initial) for key in initial_proof['context_counts']} == initial_proof['context_counts']
            and len({row[3] for row in initial}) == len(initial) and set(map(tuple, initial)) <= set(map(tuple, schedulers)),
            'Wrong initial73 scheduler evidence (not a final catalog)')
    require(final['schema'] == 'ecoda_orphaned_main_emitted_scheduler_inventory_v1'
            and final['final_emitted_catalog_complete'] is True and final['original_files_written'] is False
            and final['catalog_sha256'] == digest(paths['emitted_catalog'])
            and final['initial_catalog'] == incident['initial_catalog'] and final['initial_inventory'] == incident['initial_inventory']
            and final['single_inspect_scheduler_ids'] == ids, 'Final emitted scheduler inventory is not sealed to the sole inspect')
    children = {'benchmark': MAIN_RUNS['benchmark'], 'uncorrected': MAIN_RUNS['batch_effect/uncorrected_final'],
                'corrected': MAIN_RUNS['batch_effect/corrected_final']}
    ledgers = initial_proof['original_child_ledgers']; expected = set()
    require(len(ledgers) == 3 and {row['context'] for row in ledgers} == set(children)
            and final['original_child_ledgers'] == ledgers, 'All three immutable original ledgers are required')
    for row in ledgers:
        path = proof(row); child = original.parent / children[row['context']]
        require(path == child / 'manifests/scheduler_ids.tsv'
                and digest(str(path) + '.md5') == row['sidecar_sha256']
                and fields(str(path) + '.md5') == dict(PATH=str(path), MD5=digest(path, 'md5'), SIZE=str(path.stat().st_size)),
                'Original child ledger/proof changed')
        values = lines(path)
        require(values and all(len(value) == 2 and value[0] in {'ARRAY', 'WATCHDOG', 'STATUS', 'AGGREGATE_GATE', 'PREFLIGHT', 'METADATA_EXPORT'}
                and re.fullmatch('[1-9][0-9]*', value[1]) for value in values), 'Invalid original child ledger')
        expected.update((row['context'], child.name, *value) for value in values)
    additional = {('benchmark', children['benchmark'], 'AGGREGATE_GATE', '4665612'),
                  ('corrected', children['corrected'], 'AGGREGATE_GATE', '4665633')}
    require(expected | additional == set(map(tuple, initial)), 'Initial catalog omits/adds original emitted IDs')
    descriptors = observation['active_job_descriptors']
    require(any(item.get('exit') == 0 and 'JobId=4665612' in item.get('stdout', '') and children['benchmark'] in item['stdout']
                and 'matrix_gate.sh' in item['stdout'] for item in descriptors.values()),
            'Benchmark aggregate descriptor witness missing')
    aggregates = read_json(paths['aggregate_evidence'])
    require(any(item.get('exists') is True and item.get('path') == str(original.parent / children['corrected'] / 'logs/5_matrix_gate_4665633.log') for item in aggregates), 'Corrected aggregate log witness missing')
    for item in aggregates:
        if item.get('path', '').endswith('/5_matrix_gate_4665633.log') and item.get('exists'):
            proof(item)
    expected |= additional
    watchdogs = {row[3]: row[:2] for row in initial if row[2] == 'WATCHDOG'}
    coverage = verify_no_retry_coverage(incident, paths, initial, original) if 'no_retry_coverage' in incident else None
    watches = final['watchdog_evidence']
    require(len(watches) == len(watchdogs) and {item['watchdog_id'] for item in watches} == set(watchdogs),
            'Final catalog lacks complete watchdog/retry evidence')
    for item in watches:
        context, child_name = watchdogs[item['watchdog_id']]
        require(item['context'] == context, 'Foreign watchdog evidence')
        child = original.parent / child_name
        require(item['stdout']['exists'] is True or coverage is not None, 'Complete retry coverage requires every original watchdog stdout')
        if coverage is not None:
            require(item['watchdog_id'] in coverage['_watchdogs']
                    and all(item[kind] == coverage['_watchdogs'][item['watchdog_id']][kind] for kind in ('stdout','stderr','status')),
                    'Final watchdog evidence differs from approved precondition coverage')
        retry_ids = set()
        log_match = re.fullmatch(r'5_matrix_watchdog_(.+)_' + item['watchdog_id'] + r'\.log', Path(item['stdout']['path']).name)
        require(log_match is not None, 'Watchdog stdout is not its production log')
        label = log_match[1]
        for kind in ('stdout', 'stderr', 'status'):
            snapshot = item[kind]; path = safe(snapshot['path'], snapshot['exists'])
            expected_path = child / ('status/watchdogs/' + label + '.status' if kind == 'status'
                                    else 'logs/5_matrix_watchdog_' + label + '_' + item['watchdog_id'] + ('.log' if kind == 'stdout' else '.err'))
            require(path == expected_path, 'Foreign watchdog snapshot')
            if snapshot['exists']:
                proof(snapshot)
                text = path.read_text()
                if kind == 'status':
                    require('LABEL=' + label in text.splitlines(), 'Wrong watchdog status label')
                    retry_ids.update(re.findall(r'^SCHEDULER_ID=([1-9][0-9]*)$', text, re.M))
                else:
                    retry_ids.update(re.findall(r'^(?:BATCH_EFFECT|MATRIX)_RETRY_ARRAY_JOB_ID=([1-9][0-9]*)$', text, re.M))
            else:
                require(not path.exists(), 'An absent watchdog snapshot changed')
        for job in retry_ids:
            if job not in {row[3] for row in initial}:
                expected.add((context, child_name, 'RETRY_ARRAY', job))
    require(expected == set(map(tuple, schedulers)), 'Final emitted catalog differs from complete original/retry evidence')
    if coverage is not None:
        require(expected == set(map(tuple, initial)), 'No-new-retry coverage cannot add scheduler IDs')
    queued = quiescence['queue']['stdout']
    require(not any(re.search(r'(?<![0-9])' + job + r'(?:_|\b)', queued) for job in ids), 'An original scheduler job remains queued')
    return coverage


def verified_main_v2_scope(control, science, review, run, nas):
    """Fixed failed main-v2 input; the reviewed finite audit defines the rows.

    This checks the audit's bytes and complete original terminal accounting,
    not method internals. Existing RDS/Feather validators must have established
    the expected payload union, current labels/order and scITD drops first.
    """
    context = control['context']
    reboot = control['schema'] == REBOOT_SCHEMA
    gate_state = 'PRELAUNCH_STOP' if reboot else 'FAILED'
    receipt_schema = 'ecoda_main_v2_reboot_science_receipt_v1' if reboot else 'ecoda_main_v2_science_receipt_v1'
    require(context in MAIN_RUNS and science['schema'] == receipt_schema
            and science['context'] == context and science['status'] == 'NOOP_VALIDATED'
            and science['science_validated'] is True and science['publication_completed'] is False
            and science['original_gate_status'] == gate_state and science['source_commit'] == MAIN_SOURCE,
            'Wrong fixed main-v2 science receipt')
    original = safe(science['original_run_root'])
    require(original.parent == run.parent and original.name == MAIN_RUNS[context] and original != run,
            'Wrong main-v2 original child')
    evidence = {key: proof(value) for key, value in science['original_evidence'].items()}
    required = {'gate_manifest', 'inspect', 'profile', 'accounting_raw', 'campaign_contexts', 'campaign_scheduler_ids',
                'child_scheduler_ids', 'source_manifest', 'runtime_identity', 'artifact_binding', 'input_artifacts',
                'source_identity', 'selection', 'dispatch_selection', 'owners_manifest', 'preflight_manifest'}
    require(set(evidence) == required, 'Incomplete main-v2 original evidence')
    for key, basename in (('child_scheduler_ids', 'scheduler_ids.tsv'), ('source_manifest', 'source.manifest'),
                          ('runtime_identity', 'runtime.identity'), ('artifact_binding', 'artifact_bindings.json'),
                          ('input_artifacts', 'input_artifacts.tsv'), ('selection', 'selection.tsv'),
                          ('dispatch_selection', 'dispatch_selection.tsv'), ('owners_manifest', 'owners.tsv'),
                          ('preflight_manifest', 'h5ad_preflight.tsv')):
        require(evidence[key] == original / 'manifests' / basename, 'Wrong original child proof namespace')
    require(evidence['source_identity'] == original / 'manifests/source_identity.json', 'Wrong original sample identity namespace')
    source = fields(evidence['source_manifest']); runtime = fields(evidence['runtime_identity'])
    config = safe(source['SOURCE_ROOT'] + '/datasets.json')
    configuration = read_json(config)
    require(source['SOURCE_COMMIT'] == MAIN_SOURCE and source['AUX_ROOT'] == source['SOURCE_ROOT'] + '/aux'
            and digest(config) == source['DATASETS_SHA256'], 'Original main-v2 source/config differs')
    require(digest(runtime['RUNTIME_MANIFEST']) == runtime['RUNTIME_MANIFEST_SHA256']
            and fields(runtime['RUNTIME_MANIFEST'])['IMAGE_SHA256'] == runtime['RUNTIME_IMAGE_SHA256']
            and re.fullmatch('[a-f0-9]{64}', runtime['RUNTIME_IMAGE_SHA256'])
            and safe(runtime['RUNTIME_IMAGE']).stat().st_size == int(runtime['RUNTIME_IMAGE_SIZE']), 'Original main-v2 runtime differs')
    binding = read_json(evidence['artifact_binding'])
    require(binding['schema'] == 'ecoda_derivative_artifacts_v1' and binding['stage'] == 'stage5'
            and binding['run_id'] == original.name and binding['run_root'] == str(original)
            and binding['selection'] == str(evidence['selection']) and binding['input_artifacts'] == str(evidence['input_artifacts']),
            'Original derivative binding differs')
    for path in (config, evidence['source_manifest'], evidence['runtime_identity'], evidence['selection'], evidence['input_artifacts']):
        require(binding['identities'].get(str(path)) == digest(path), 'Original binding proof changed')
    gate = read_json(evidence['gate_manifest']); inspected = read_json(evidence['inspect']); profile = read_json(evidence['profile'])
    require(gate['gate_id'] == 'ecoda_stage5_main_confounds_20261109_v2'
            and gate['state'] == inspected['state'] == gate_state and inspected['operation'] == 'inspect'
            and inspected['gate_id'] == gate['gate_id'] and gate['audit'] == inspected['audit']
            and gate['audit_completed'] is True and gate['audit_state'] == 'COMPLETED'
            and gate['audit']['passed'] is False and gate['profile'] == profile['name']
            and gate['profile_digest'] == digest(evidence['profile'])
            and gate['command_digest'] == hashlib.sha256(gate['exact_command'].encode()).hexdigest(), 'Original failed gate/inspect/profile differs')
    if not reboot:
        wrapper_status = gate.get('terminal_status')
        require(isinstance(wrapper_status, dict) and wrapper_status.get('schema_version') == '1.0'
                and wrapper_status.get('state') in {'FAILED', 'COMPLETED'}
                and type(wrapper_status.get('exit_code')) is int
                and (wrapper_status['exit_code'] == 0) == (wrapper_status['state'] == 'COMPLETED')
                and all(wrapper_status.get(key) == gate[key] for key in ('gate_id', 'exact_command', 'command_digest', 'event_generation'))
                and isinstance(gate['event_generation'], str) and bool(gate['event_generation'])
                and all(isinstance(wrapper_status.get(key), str) and wrapper_status[key] for key in ('started_at', 'finished_at')),
                'Original wrapper terminal status is missing, nonterminal or stale')
    argv = shlex.split(gate['exact_command'])
    def option(name):
        values = [argv[i+1] for i, word in enumerate(argv[:-1]) if word == name]
        require(len(values) == 1, 'Ambiguous original wrapper option ' + name)
        return values[0]
    require(option('--run-id') == 'stage5_main_confounds_20261109_v2'
            and option('--scratch-root') == str(original.parent.parent)
            and option('--source-root') == source['SOURCE_ROOT']
            and safe(option('--source-manifest')).read_bytes() == evidence['source_manifest'].read_bytes()
            and option('--runtime-image') == runtime['RUNTIME_IMAGE']
            and option('--runtime-manifest') == runtime['RUNTIME_MANIFEST']
            and option('--script') == 'src/5_run_benchmark_methods/stage5_campaign.sh', 'Wrong original main-v2 wrapper')
    ids = gate['scheduler_array_ids'] + gate['scheduler_watchdog_ids']
    require(ids and len(set(ids)) == len(ids) and all(re.fullmatch('[1-9][0-9]*', i) for i in ids), 'Incomplete original scheduler union')
    raw = inspected['audit']['accounting']; template = 'sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode'
    require(profile['policy']['accounting_command'] == template and raw['attempted'] is True
            and type(raw['query_count']) is int and raw['query_count'] == 1
            and type(raw['actual_exit']) is int and raw['actual_exit'] == 0 and raw['scheduler_ids'] == ids
            and raw['command'] == template.replace('{scheduler_ids}', ','.join(ids))
            and evidence['accounting_raw'].read_bytes() == raw['stdout'].encode()
            and raw['stdout'].endswith('\n'), 'Original single complete accounting query differs')
    parsed = [line.split('|') for line in raw['stdout'].splitlines()]
    terminal = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED', 'SPECIAL_EXIT'}
    require(parsed and all(len(r) == 3 and re.fullmatch('[1-9][0-9]*(?:_[1-9][0-9]*)?', r[0])
            and r[1].rstrip('+') in terminal and re.fullmatch('[0-9]+:[0-9]+', r[2]) for r in parsed)
            and len({r[0] for r in parsed}) == len(parsed) and set(ids) <= {r[0] for r in parsed}
            and raw['rows'] == [dict(line=i, job_id=r[0], state=r[1], exit_code=r[2]) for i, r in enumerate(parsed, 1)]
            and science['original_accounting'] == raw, 'Original accounting is incomplete/nonterminal')
    rows = lambda path: [line.split('\t') for line in path.read_text().splitlines()]
    schedulers = rows(evidence['campaign_scheduler_ids']); child_ids = rows(evidence['child_scheduler_ids'])
    contexts = rows(evidence['campaign_contexts'])
    expected_children = {'benchmark': MAIN_RUNS['benchmark'], 'uncorrected': MAIN_RUNS['batch_effect/uncorrected_final'], 'corrected': MAIN_RUNS['batch_effect/corrected_final']}
    require(len(contexts) == 3 and all(len(r) == 4 and expected_children.get(r[0]) == r[1] for r in contexts)
            and {r[0] for r in contexts} == set(expected_children)
            and all(len(r) == 4 and expected_children.get(r[0]) == r[1] for r in schedulers)
            and len({r[3] for r in schedulers}) == len(schedulers) and {r[3] for r in schedulers} == set(ids)
            and child_ids and all(len(r) == 2 for r in child_ids)
            and (set(map(tuple, child_ids)) <= {tuple(r[2:]) for r in schedulers if r[1] == original.name} if reboot
                 else child_ids == [r[2:] for r in schedulers if r[1] == original.name]), 'Original campaign/child scheduler scope differs')
    parent = original.parent / 'stage5_main_confounds_20261109_v2'
    require(evidence['campaign_contexts'] == parent / 'manifests/campaign_contexts.tsv', 'Wrong original parent context namespace')
    reboot_coverage = None
    if reboot:
        reboot_coverage = verify_login_reboot(control, science, review, gate, inspected, evidence, schedulers, ids, original)
    else:
        require(evidence['campaign_scheduler_ids'] == parent / 'manifests/campaign_scheduler_ids.tsv', 'Wrong original parent scheduler namespace')
    selected_context = next(r for r in contexts if r[1] == original.name)
    selected_scope = evidence['selection']
    if context == 'benchmark':
        selected_scope = safe(binding['requested_selection'])
        require(selected_scope == original / 'manifests/requested_selection.tsv'
                and binding['identities'].get(str(selected_scope)) == digest(selected_scope), 'Original requested-scope binding differs')
    require(safe(selected_context[2]).read_bytes() == selected_scope.read_bytes()
            and safe(selected_context[3]).read_bytes() == evidence['input_artifacts'].read_bytes(), 'Original parent/child scientific scope differs')
    catalog_path = proof(science['catalog']); catalog = read_json(catalog_path)
    validation_path = proof(science['scientific_validation']); validation = read_json(validation_path)
    require(catalog['schema'] == 'ecoda_main_v2_publication_catalog_v1' and catalog['context'] == context
            and validation['schema'] == 'ecoda_main_v2_artifact_validation_v1'
            and validation['catalog_sha256'] == science['catalog']['sha256']
            and validation['source_identity_sha256'] == digest(evidence['source_identity'])
            and validation['config_sha256'] == digest(config) and validation['new_model_fit'] is False,
            'Scientific validation is not bound to this current catalog/input/config')
    for item in validation['execution_evidence']:
        proof(item)
    require(validation['execution_evidence'], 'Scientific validator source/runtime/log evidence missing')
    contracts = {'src/utils/bash/ecoda_run_common.sh', 'src/utils/bash/ecoda_stage5_policy.sh',
                 'src/utils/bash/ecoda_derivative_artifacts.sh', 'src/utils/py/derivative_artifacts.py',
                 'src/utils/py/benchmark_request_scope.py'}
    require({str(proof(item).relative_to(source['SOURCE_ROOT'])) for item in validation['contract_sources']} == contracts,
            'Original artifact contract sources missing or foreign')
    inventory_path = proof(science['inventory']); inventory_rows = rows(inventory_path)
    keys = ('dataset', 'kind', 'source', 'relative', 'md5', 'size', 'original_producer')
    require(inventory_rows and all(len(r) == 7 and all(r) for r in inventory_rows), 'Malformed finite publication inventory')
    inventory = [dict(zip(keys, r)) for r in inventory_rows]
    require(len({r['source'] for r in inventory}) == len(inventory) and len({r['relative'] for r in inventory}) == len(inventory), 'Duplicate publication payload')
    needed_validators = {'rds' if Path(r['source']).suffix == '.rds' else 'feather'
                         for r in inventory if r['kind'] != 'timing'}
    observed_validators = set(); validator_inputs = set(); rds_inputs = set()
    for invocation in validation['validator_runs']:
        kind = invocation['kind']; script = proof(invocation)
        expected = {'rds': ('validate_benchmark_rds_contract.R', 'benchmark RDS artifact-list contract OK'),
                    'feather': ('matrix_artifact_validator.py', 'matrix artifact contract OK')}
        require(kind in expected and script == Path(source['SOURCE_ROOT']) / 'src/5_run_benchmark_methods' / expected[kind][0]
                and type(invocation['actual_exit']) is int and invocation['actual_exit'] == 0
                and str(script) in invocation['argv']
                and cli_value(invocation['argv'], '--config') == str(config)
                and cli_value(invocation['argv'], '--source-identity') == str(evidence['source_identity'])
                and expected[kind][1] in proof(invocation['stdout']).read_text().splitlines(), 'Missing/failed current production artifact validator')
        proof(invocation['stderr']); observed_validators.add(kind)
        argv = invocation['argv']
        require(not any(word.startswith('--') and 'historical' in word for word in argv), 'Historical validator relaxation is prohibited')
        expected_pass = None if context == 'benchmark' else 'corrected' if context.endswith('/corrected_final') else 'uncorrected'
        expected_variant = None if context == 'benchmark' else 'corrected_final' if context.endswith('/corrected_final') else 'final'
        require(cli_value(argv, '--batch-pass') == expected_pass and cli_value(argv, '--analysis-variant') == expected_variant,
                'Production validator batch/variant mode differs from publication context')
        if kind == 'feather' and expected_pass == 'corrected':
            contract = proof(invocation['expected_batch_contract'])
            require(original / 'manifests' in contract.parents and cli_value(argv, '--expected-batch-contract') == str(contract),
                    'Corrected Feather validator source contract is missing/foreign')
        elif expected_pass != 'corrected':
            require(cli_value(argv, '--expected-batch-contract') is None and cli_value(argv, '--batch-contract') is None,
                    'Unexpected correction contract in benchmark/uncorrected validation')
        if kind == 'rds':
            inputs_file = proof(invocation['input_list'])
            require(cli_value(argv, '--artifact-list') == str(inputs_file)
                    and all(cli_value(argv, name) is None for name in ('--artifact', '--root', '--selection', '--dataset', '--view', '--method'))
                    and '--metadata' not in argv, 'Actual RDS validator input list differs')
            inputs_rows = rows(inputs_file)
            require(inputs_rows and all(len(r) == 5 and all(r) and r[4] in {'0', '1'} for r in inputs_rows)
                    and len({tuple(r) for r in inputs_rows}) == len(inputs_rows), 'Malformed actual RDS validator input list')
            validator_inputs.update(r[0] for r in inputs_rows)
            rds_inputs.update(tuple(r) for r in inputs_rows)
        else:
            require(argv.count('--artifact') == 1, 'Feather validator must name its exact artifact')
            validator_inputs.add(argv[argv.index('--artifact')+1])
    require(observed_validators == needed_validators, 'Production validator kinds differ from exact selected payloads')
    require(validator_inputs == {r['source'] for r in inventory if r['kind'] != 'timing'}, 'Actual production validator input union differs from selected science inventory')
    dispatch = {tuple(r[:3]) for r in rows(evidence['dispatch_selection'])}
    source_ids = read_json(evidence['source_identity'])['entries']
    identities = {(r['dataset'], r['view']): r for r in source_ids}
    references = rows(evidence['input_artifacts']); input_map = {(r[0], r[1]): r for r in references if len(r) == 6}
    require(len(identities) == len(source_ids) and len(input_map) == len(references), 'Duplicate/malformed original input identities')
    allowed_arrays = {r[1] for r in child_ids if r[0] == 'ARRAY'}
    checked_by_path = {}
    covered = set(); expected_payloads = set(); expected_timings = set(); expected_rds_inputs = set(); validated = validation['rows']
    require(catalog['rows'] and len(validated) == len(catalog['rows']), 'Missing selected scientific validation rows')
    for selected, checked in zip(catalog['rows'], validated):
        ds, view, method, parameter = (selected[k] for k in ('dataset', 'view', 'method', 'parameter'))
        key = (ds, view, method, parameter)
        require(key not in covered and (ds, view, method) in dispatch, 'Duplicate/foreign selected scientific row')
        expected_view = 'benchmark_analysis' if context == 'benchmark' else 'batch_effect_' + ('corrected' if context.endswith('/corrected_final') else 'uncorrected')
        entry = configuration.get(ds, {})
        require(not ds.startswith('_') and view == expected_view and entry.get('use_for_benchmark' if context == 'benchmark' else 'use_for_batch_effect') is True
                and view in entry.get('views', {}) and 'columns' not in entry['views'][view], 'Selected cohort/view is not production-enabled')
        covered.add(key)
        require(selected['array_id'] in allowed_arrays and type(selected['task_index']) is int and selected['task_index'] > 0,
                'Selected worker is not in the original array scope')
        # JobIDRaw can be numeric for physical array tasks. The original
        # three-field query cannot reconstruct that mapping retrospectively.
        # Publish validated artifacts, never invent a zero worker exit code.
        require(checked['completion_criterion'] == 'VALIDATED_PUBLISHED_ARTIFACTS'
                and checked['worker_exit_not_established'] is True, 'Unproven scientific worker success claim')
        manifest = proof(selected['manifest']); physical = rows(manifest)
        require(manifest.parent == original / 'manifests' and manifest.name.startswith('matrix_')
                and fields(str(manifest) + '.md5') == dict(PATH=str(manifest), MD5=digest(manifest, 'md5'), SIZE=str(manifest.stat().st_size))
                and selected['task_index'] <= len(physical)
                and physical[selected['task_index']-1] == list(key[:3]) + ([parameter] if parameter else []), 'Selected row/task/manifest differs')
        for name in ('stdout', 'stderr'):
            log = proof(selected[name])
            require(original / 'logs' in log.parents
                    and f"_{selected['array_id']}_{selected['task_index']}." in log.name, 'Worker evidence outside exact original logical task logs')
        identity = identities[(ds, view)]; reference = input_map[(ds, view)]
        from derivative_artifacts import configured_files
        require(Path(reference[2]).name in configured_files(configuration, ds, view, 'stage5'), 'Original input filename differs from registry')
        require(identity['path'] == reference[2] and identity['md5'] == reference[3] and str(identity['size']) == reference[4], 'Current source identity/input reference differs')
        checksum = fields(reference[5])
        require(checksum == dict(PATH=reference[2], MD5=reference[3], SIZE=reference[4]), 'Original input checksum proof changed')
        require(checked['row'] == list(key) and checked['input_md5'] == identity['md5']
                and checked['artifact_contract_valid'] is True and checked['current_labels_valid'] is True
                and checked['sample_order_valid'] is True and checked['timing_valid'] is True,
                'Selected scientific row lacks current artifact/label/order/timing validation')
        retained = checked['sample_ids']; dropped = checked['dropped_sample_ids']
        require(retained and len(set(retained)) == len(retained) and len(set(dropped)) == len(dropped)
                and not set(retained) & set(dropped) and set(retained) | set(dropped) == set(identity['sample_ids'])
                and retained == [sample for sample in identity['sample_ids'] if sample not in set(dropped)]
                and (method == 'scitd' or not dropped), 'Unapproved sample drop/order (only scITD may drop)')
        require(len(checked['labels']) == len(retained) and all(isinstance(v, str) and v for v in checked['labels']), 'Current label vector missing')
        expected, expected_inputs = main_v2_expected_artifacts(original, source, config, binding, context, selected)
        expected_rds_inputs.update(expected_inputs)
        require(expected and set(checked['expected_artifacts']) == expected, 'Incomplete per-row production artifact contract')
        for expected_path in expected:
            checked_by_path[expected_path] = checked
        expected_payloads.update(expected)
        require(checked['timing_artifacts'], 'Selected timing input list missing')
        expected_timings.update(checked['timing_artifacts'])
    require(expected_payloads == {r['source'] for r in inventory if r['kind'] != 'timing'}
            and expected_timings == {r['source'] for r in inventory if r['kind'] == 'timing'}, 'Expected-row contract differs from exact finite payload union')
    require(rds_inputs == expected_rds_inputs, 'Actual RDS validator dataset/view/method/metadata columns differ from production contract')
    preflight_rows = rows(evidence['preflight_manifest'])
    wanted_inputs = {(r[0], r[1]) for r in covered}
    require(len(science['input_statuses']) == len(wanted_inputs), 'Current input preflight witnesses missing')
    for witness in science['input_statuses']:
        key = (witness['dataset'], witness['view']); status_path = proof(witness['status']); status = fields(status_path)
        require(key in wanted_inputs and status_path == original / 'status/h5ad_preflight' / (key[0] + '__' + key[1] + '.status')
                and status['STATE'] == 'OK' and status['RUN_ID'] == original.name
                and status['DATASET'] == key[0] and status['VIEW'] == key[1]
                and status['TASK_ID'].isdigit() and 0 < int(status['TASK_ID']) <= len(preflight_rows)
                and preflight_rows[int(status['TASK_ID'])-1] == [key[0], key[1], input_map[key][2]], 'Original full-input preflight witness differs')
        wanted_inputs.remove(key)
    require(review['status'] == 'APPROVED_SCIENTIFIC_ARTIFACT_REUSE' and review['receipt_sha256'] == control['science_receipt']['sha256']
            and review['catalog_sha256'] == science['catalog']['sha256'] and review['inventory_sha256'] == science['inventory']['sha256']
            and review['context'] == context and review['selected_rows'] == len(covered) and review['payloads'] == len(inventory)
            and review['worker_exit_not_established'] is True
            and review['original_gate_status'] == gate_state and review['publication_completed'] is False, 'Review does not approve this exact finite main-v2 scope')
    snapshots = science['status_snapshots']
    require({str(original / 'status/terminal'), str(original / 'status/aggregate')} <= {x['path'] for x in snapshots}, 'Original terminal/aggregate snapshots missing')
    for snapshot in snapshots:
        path = safe(snapshot['path'], snapshot['exists'])
        require(original / 'status' in path.parents and type(snapshot['exists']) is bool, 'Foreign status snapshot')
        if snapshot['exists']:
            proof(snapshot)
        else:
            require(not path.exists(), 'Originally absent status was created')
    owner_rows = rows(evidence['owners_manifest']); expected_owner_keys = set()
    scope = 'ordinary' if context == 'benchmark' else ('corrected' if context.endswith('/corrected_final') else 'uncorrected')
    for ds, view, method, _ in covered:
        expected_owner_keys.add(f'derivative/{original.name}/{scope}/{ds}/{view}/{method}')
    require(len(science['row_owners']) == len(expected_owner_keys), 'Selected row ownership evidence missing')
    for snapshot in science['row_owners']:
        path = proof(snapshot); value = fields(path)
        require(value == snapshot['fields'] and value['RUN_ID'] == original.name and value['STAGE'] == 'stage5'
                and value['STATE'] in {'OK', 'FAIL', 'ACTIVE'} and value['KEY'] in expected_owner_keys
                and path == original.parent.parent / '_ecoda_owners/stage5' / re.sub(r'[/,:\t |]', '_', value['KEY']) / 'owner'
                and [value['KEY'], str(path.parent)] in owner_rows, 'Original row ownership snapshot changed')
        expected_owner_keys.remove(value['KEY'])
    require(len(science['payloads']) == len(inventory), 'Payload proof coverage differs')
    for row, item in zip(inventory, science['payloads']):
        require(all(item[k] == row[k] for k in keys), 'Approved payload inventory order differs')
        path = safe(row['source']); relative = Path(row['relative'])
        require(not relative.is_absolute() and '..' not in relative.parts and str(relative) == row['relative']
                and relative.parts[0] in {'results', 'embeddings', 'pseudobulks', 'gloscope_dists'}
                and row['kind'] in {'result', 'metadata', 'cache', 'timing'} and path.suffix in {'.rds', '.feather'}, 'Unsafe/unselected payload kind/path')
        expected = original / ('logs/' + path.name if row['kind'] == 'timing' else 'artifacts/' + context + '/' + str(relative))
        require(path == expected and (row['kind'] != 'timing' or str(relative) == 'embeddings/timing_sources/' + path.name), 'Source outside exact original context')
        current = record(path, original, row['original_producer']); proof(item['producer_record']); proof(item['sidecar'])
        require(current['md5'] == row['md5'] and current['size'] == row['size'] and current['record_path'] == item['producer_record']['path']
                and item['sidecar']['path'] == str(path) + '.md5', 'Original payload proof changed')
        if item['owner'] is not None:
            require(owner(path, original, {'OK', 'FAIL', 'ACTIVE'}) == item['owner'], 'Original payload ownership snapshot changed')
        else:
            if row['kind'] == 'timing':
                require(not owner_file(path, original).exists(), 'An unowned timing shard acquired an owner')
            else:
                verify_ownerless_gloscope(row, item, context, original, checked_by_path[str(path)], reboot_coverage)
        for runtime_proof in item['runtime_evidence']:
            proof(runtime_proof)
    adoption = original / 'manifests/prepared_adoption.tsv'
    if adoption.exists() and any(r[4] in expected_payloads for r in rows(adoption)):
        required_lineage = {str(adoption), str(original / 'manifests/prepared_adoption.receipt.json'),
                            str(original / 'manifests/prepared_adoption.json')}
        require({str(proof(item)) for item in science['adoption_lineage']} == required_lineage, 'Selected prepared artifact adoption lineage missing')
        adopted_receipt = original / 'manifests/prepared_adoption.receipt.json'
        adopted_report = original / 'manifests/prepared_adoption.json'
        lineage = read_json(adopted_report)
        require(lineage['operation'] == 'COPY_WITHOUT_REFIT' and lineage['original_gate_state'] == 'FAILED'
                and lineage['receipt_path'] == str(adopted_receipt) and lineage['receipt_sha256'] == digest(adopted_receipt)
                and lineage['binding_sha256'] == digest(evidence['artifact_binding']), 'Prepared adoption provenance differs')
        record(adopted_receipt, original, 'stage5_prepared_adoption_receipt')
        record(adopted_report, original, 'stage5_prepared_adoption')
    destination = safe(str(nas / 'releases' / control['release_id'] / run.name / context), False)
    return control, science, inventory, run, nas, destination


def verified_unc4_scope(control, science, review, run, nas):
    """The fixed0314 four-MrVI fits passed; only their final controller failed."""
    original = safe(science['original_run_root']); context = 'batch_effect/uncorrected_final'
    require(original.parent == run.parent and original.name == UNC4_RUN and original != run
            and science['schema'] == 'ecoda_failed16_unc4_science_receipt_v1'
            and science['status'] == 'NOOP_VALIDATED' and science['science_validated'] is True
            and science['publication_completed'] is False and science['original_gate_status'] == 'FAILED'
            and science['source_commit'] == UNC4_SOURCE and science['worker_exit_not_established'] is True,
            'Wrong fixed uncorrected-four science receipt')
    require(review['status'] == 'APPROVED_SCIENTIFIC_ARTIFACT_REUSE'
            and review['receipt_sha256'] == control['science_receipt']['sha256']
            and review['validation_sha256'] == science['validation']['sha256']
            and review['companions_sha256'] == science['companions']['sha256']
            and review['datasets'] == review['distance_feathers'] == review['runtime_jsons'] == review['timing_feathers'] == 4
            and review['payloads'] == 12 and review['original_gate_status'] == 'FAILED'
            and review['publication_completed'] is False and review['worker_exit_not_established'] is True,
            'Review does not approve the exact four-fit/12-payload release')
    evidence = {key: proof(value) for key, value in science['original_evidence'].items()}
    names = {'source_manifest': 'manifests/source.manifest', 'runtime_identity': 'manifests/runtime.identity',
             'binding': 'manifests/artifact_bindings.json', 'source_identity': 'manifests/source_identity.json',
             'input_artifacts': 'manifests/input_artifacts.tsv', 'selection': 'manifests/selection.tsv',
             'dispatch': 'manifests/dispatch_selection.tsv', 'child_scheduler_ids': 'manifests/scheduler_ids.tsv',
             'metadata': 'metadata', 'owners_manifest': 'manifests/owners.tsv', 'terminal': 'status/terminal', 'aggregate': 'status/aggregate',
             'watchdog': 'status/watchdogs/batch_effect_uncorrected__mrvi.status'}
    require(set(evidence) == set(names) | {'gate_manifest', 'inspect', 'profile', 'matrix', 'parent_scheduler_ids'},
            'Incomplete original fixed-four proof set')
    for name, relative in names.items():
        require(evidence[name] == original / relative, 'Wrong fixed-four original proof namespace')
    source = fields(evidence['source_manifest']); runtime = fields(evidence['runtime_identity'])
    config = safe(source['SOURCE_ROOT'] + '/datasets.json')
    require(source['SOURCE_COMMIT'] == UNC4_SOURCE and source['AUX_ROOT'] == source['SOURCE_ROOT'] + '/aux'
            and digest(config) == source['DATASETS_SHA256'], 'Original four-fit source/config changed')
    require(digest(runtime['RUNTIME_MANIFEST']) == runtime['RUNTIME_MANIFEST_SHA256']
            and fields(runtime['RUNTIME_MANIFEST'])['IMAGE_SHA256'] == runtime['RUNTIME_IMAGE_SHA256']
            and safe(runtime['RUNTIME_IMAGE']).stat().st_size == int(runtime['RUNTIME_IMAGE_SIZE']), 'Original runtime changed')
    binding = read_json(evidence['binding'])
    require(binding['schema'] == 'ecoda_derivative_artifacts_v1' and binding['stage'] == 'stage5'
            and binding['run_id'] == original.name and binding['run_root'] == str(original)
            and binding['selection'] == str(evidence['selection']) and binding['input_artifacts'] == str(evidence['input_artifacts']),
            'Original fixed-four input/output binding differs')
    for path in (config, evidence['source_manifest'], evidence['runtime_identity'], evidence['selection'], evidence['input_artifacts']):
        require(binding['identities'].get(str(path)) == digest(path), 'Original bound proof changed')
    gate = read_json(evidence['gate_manifest']); inspected = read_json(evidence['inspect']); profile = read_json(evidence['profile'])
    require(gate['gate_id'] == 'ecoda_stage5_failed16_20261109_v1' and gate['state'] == inspected['state'] == 'FAILED'
            and inspected['operation'] == 'inspect' and inspected['gate_id'] == gate['gate_id'] and gate['audit'] == inspected['audit']
            and gate['audit_completed'] is True and gate['audit_state'] == 'COMPLETED' and gate['audit']['passed'] is False
            and gate['profile'] == profile['name'] and gate['profile_digest'] == digest(evidence['profile'])
            and gate['command_digest'] == hashlib.sha256(gate['exact_command'].encode()).hexdigest(), 'Original failed16 terminal audit differs')
    terminal = gate['terminal_status']
    require(terminal['schema_version'] == '1.0' and terminal['state'] == 'FAILED' and type(terminal['exit_code']) is int
            and terminal['exit_code'] != 0 and all(terminal[key] == gate[key] for key in ('gate_id', 'exact_command', 'command_digest', 'event_generation'))
            and isinstance(gate['event_generation'], str) and bool(gate['event_generation'])
            and all(isinstance(terminal.get(key), str) and terminal[key] for key in ('started_at', 'finished_at')),
            'Original failed wrapper outcome is not established')
    argv = shlex.split(gate['exact_command'])
    require(cli_value(argv, '--run-id') == 'stage5_failed16_20261109_v1'
            and cli_value(argv, '--source-root') == source['SOURCE_ROOT']
            and safe(cli_value(argv, '--source-manifest')).read_bytes() == evidence['source_manifest'].read_bytes()
            and cli_value(argv, '--scratch-root') == str(original.parent.parent)
            and cli_value(argv, '--runtime-image') == runtime['RUNTIME_IMAGE']
            and cli_value(argv, '--runtime-manifest') == runtime['RUNTIME_MANIFEST']
            and cli_value(argv, '--script') == 'src/5_run_benchmark_methods/stage5_campaign.sh', 'Wrong original failed16 command')
    ids = gate['scheduler_array_ids'] + gate['scheduler_watchdog_ids']; accounting = inspected['audit']['accounting']
    expected_ids = {'4666176', '4666180', '4666183', '4666185', '4666178', '4666187', '4666192', '4666196', '4666184', '4666191'}
    template = 'sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode'
    raw = [line.split('|') for line in accounting['stdout'].splitlines()]
    require(len(ids) == len(set(ids)) == 10 and set(ids) == expected_ids and accounting['scheduler_ids'] == ids
            and profile['policy']['accounting_command'] == template and accounting['command'] == template.replace('{scheduler_ids}', ','.join(ids))
            and accounting['query_count'] == 1 and accounting['attempted'] is True and accounting['actual_exit'] == 0
            and len(raw) == len({row[0] for row in raw}) == 23 and accounting['stdout'].endswith('\n')
            and all(len(row) == 3 and re.fullmatch('[1-9][0-9]*', row[0])
                    and row[1].rstrip('+') in {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED', 'SPECIAL_EXIT'}
                    and re.fullmatch('[0-9]+:[0-9]+', row[2]) for row in raw)
            and expected_ids <= {row[0] for row in raw}
            and accounting['rows'] == [dict(line=i, job_id=row[0], state=row[1], exit_code=row[2]) for i, row in enumerate(raw, 1)]
            and science['original_accounting'] == accounting, 'Original sole10-ID/23-row accounting was changed/incomplete')
    require(all(next(row[1:] for row in raw if row[0] == job) == ['COMPLETED', '0:0'] for job in ('4666187', '4666191', '4666192')),
            'Original logical array/watchdog/aggregate are not completed')
    rows = lambda path: [line.split('\t') for line in path.read_text().splitlines()]
    require(rows(evidence['child_scheduler_ids']) == [['PREFLIGHT', '4666178'], ['ARRAY', '4666187'], ['WATCHDOG', '4666191'], ['AGGREGATE_GATE', '4666192']]
            and fields(str(evidence['child_scheduler_ids']) + '.md5') == dict(PATH=str(evidence['child_scheduler_ids']), MD5=digest(evidence['child_scheduler_ids'], 'md5'), SIZE=str(evidence['child_scheduler_ids'].stat().st_size)),
            'Original child scheduler ledger differs')
    parent = original.parent / 'stage5_failed16_20261109_v1'
    schedulers = rows(evidence['parent_scheduler_ids'])
    require(evidence['parent_scheduler_ids'] == parent / 'manifests/campaign_scheduler_ids.tsv'
            and len(schedulers) == 10 and all(len(row) == 4 for row in schedulers) and {row[3] for row in schedulers} == expected_ids
            and [row[2:] for row in schedulers if row[1] == original.name] == rows(evidence['child_scheduler_ids']), 'Original parent scope differs')
    expected_selection = [[ds, 'batch_effect_uncorrected', 'batch_effect_uncorrected'] for ds in UNC4_DATASETS]
    expected_matrix = [[ds, 'batch_effect_uncorrected', 'mrvi'] for ds in UNC4_DATASETS]
    matrix = evidence['matrix']
    require(rows(evidence['selection']) == expected_selection and rows(evidence['dispatch']) == expected_matrix
            and matrix == original / 'manifests/matrix_batch_effect_uncorrected_mrvi.tsv' and rows(matrix) == expected_matrix
            and fields(str(matrix) + '.md5') == dict(PATH=str(matrix), MD5=digest(matrix, 'md5'), SIZE=str(matrix.stat().st_size)), 'Four physical worker rows differ')
    require(len(science['worker_evidence']) == 4, 'Missing exact logical worker logs')
    for index, item in enumerate(science['worker_evidence'], 1):
        require(item['dataset'] == UNC4_DATASETS[index-1] and item['task_index'] == index, 'Worker logical row mapping differs')
        for key in ('stdout', 'stderr'):
            log = proof(item[key]); suffix = '.log' if key == 'stdout' else '.err'
            require(log == original / 'logs' / (f'5_matrix_batch_effect_uncorrected__mrvi_4666187_{index}' + suffix), 'Wrong original worker log')
    wd = fields(evidence['watchdog']); aggregate_lines = evidence['aggregate'].read_text().splitlines()
    require(fields(evidence['terminal'])['STATE'] == 'FAIL'
            and fields(evidence['terminal'])['REASON'] == 'Stage 5 final matrix artifact validation failed'
            and 'STATE=OK' in aggregate_lines and aggregate_lines.count('WATCHDOG_LABELS=batch_effect_uncorrected__mrvi') == 1
            and [line for line in aggregate_lines if line.startswith('SCHEDULER_ID=')] == ['SCHEDULER_ID=' + job for job in ('4666178', '4666187', '4666191', '4666192')]
            and wd['STATE'] == 'OK' and wd['LABEL'] == 'batch_effect_uncorrected__mrvi' and wd['ARRAY_JOB_ID'] == wd['SCHEDULER_ID'] == '4666187' and wd['RETRY_INDEX'] == '0', 'Original failure/aggregate/WD evidence differs')
    validation = read_json(proof(science['validation'])); companions = read_json(proof(science['companions']))
    require(validation['schema'] == 'ecoda_failed16_unc4_readonly_validation_v1'
            and validation['status'] == 'VALIDATED_CURRENT_FOUR_MRVI_ARTIFACTS' and validation['run_id'] == original.name
            and validation['original_source_commit'] == UNC4_SOURCE
            and all(validation[key] is True for key in ('no_source_writes', 'no_accounting_query', 'no_model_fit', 'no_h5ad_payload_read'))
            and companions['schema'] == 'ecoda_failed16_unc4_companion_audit_v1' and companions['run_id'] == original.name
            and all(companions[key] is True for key in ('no_source_writes', 'no_accounting_query', 'no_scheduler_query', 'no_model_fit', 'no_h5ad_payload_read')),
            'Missing current no-compute science/companion audits')
    validation_source = fields(proof(validation['validation_source']))
    require(validation['validation_source']['commit'] == validation_source['SOURCE_COMMIT'] == '7a6c513f7f3d3bc4b1f0710029e12aef2d26aedd'
            and proof(validation['validator']) == Path(validation_source['SOURCE_ROOT']) / 'src/5_run_benchmark_methods/matrix_artifact_validator.py'
            and validation['validator']['sha256'] == 'c7b9c7f75a9e3bd0de2544d9befc8cf715e9b56a76d2ca00f6d12a036a090d97'
            and proof(validation['config']) == config and proof(validation['binding']) == evidence['binding'], 'Committed patched validator/input binding differs')
    for item in validation['unchanged']:
        proof(item)
    require({item['path'] for item in validation['unchanged']} >= {str(evidence[key]) for key in ('source_identity', 'binding', 'selection', 'metadata', 'terminal', 'aggregate')}, 'Immutable original snapshots missing')
    distances = validation['payloads']; observed = companions['rows']; identities = read_json(evidence['source_identity'])['entries']
    references = rows(evidence['input_artifacts'])
    require([item['dataset'] for item in distances] == [item['dataset'] for item in observed] == list(UNC4_DATASETS)
            and len(identities) == len(references) == 4 and [row[0] for row in references] == list(UNC4_DATASETS)
            and len(science['distance_owners']) == 4, 'Wrong exact dataset/companion/input coverage')
    require(len(science['row_owners']) == 4, 'Missing original method row owners')
    owner_rows = rows(evidence['owners_manifest'])
    for ds, snapshot in zip(UNC4_DATASETS, science['row_owners']):
        path = proof(snapshot); value = fields(path)
        key = f'derivative/{original.name}/uncorrected/{ds}/batch_effect_uncorrected/mrvi'
        require(value == snapshot['fields'] and value['RUN_ID'] == original.name and value['STAGE'] == 'stage5'
                and value['STATE'] in {'OK', 'FAIL', 'ACTIVE'} and value['KEY'] == key
                and path == original.parent.parent / '_ecoda_owners/stage5' / re.sub(r'[/,:\t |]', '_', key) / 'owner'
                and [key, str(path.parent)] in owner_rows, 'Original method owner changed/foreign')
    import pandas as pd
    inventory = []
    for index, (ds, distance, observed_row, identity, reference) in enumerate(zip(UNC4_DATASETS, distances, observed, identities, references)):
        require(identity['dataset'] == ds and identity['view'] == 'batch_effect_uncorrected' and reference[:2] == [ds, identity['view']]
                and reference[2:5] == [identity['path'], identity['md5'], str(identity['size'])]
                and fields(reference[5]) == dict(PATH=identity['path'], MD5=identity['md5'], SIZE=str(identity['size']))
                and distance['input_md5'] == identity['md5'] and distance['n_samples'] == len(identity['sample_ids']) == UNC4_SAMPLES[index]
                and distance['frame_shape'] == [UNC4_SAMPLES[index], UNC4_SAMPLES[index]], 'Current source sample/input universe differs')
        output = original / 'artifacts' / context / 'embeddings' / f'{ds}_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather'
        current = record(output, original, 'stage5_mrvi')
        require(distance['path'] == str(output) and distance['md5'] == current['md5'] and str(distance['size']) == current['size']
                and distance['producer_record'] == fields(current['record_path']), 'Validated original distance proof differs')
        require(owner(output, original, {'OK', 'FAIL', 'ACTIVE'}) == science['distance_owners'][index], 'Original distance owner changed')
        preflight = proof(distance['input_preflight']); status = fields(preflight)
        require(preflight == original / 'status/h5ad_preflight' / (ds + '__batch_effect_uncorrected.status')
                and status['STATE'] == 'OK' and status['RUN_ID'] == original.name and status['DATASET'] == ds, 'Original full-input witness differs')
        runtime_path = Path(str(output) + '.runtime.json'); timing_path = original / 'logs' / f'execution_times_batch_effect_uncorrected_final_mrvi_{ds}.feather'
        payload = read_json(runtime_path); timing = pd.read_feather(timing_path)
        require(set(payload) == {'schema_version', 'artifact_path', 'artifact_md5', 'dataset', 'method', 'time_secs', 'mem_GB'}
                and payload['schema_version'] == 1 and payload['artifact_path'] == str(output) and payload['artifact_md5'] == current['md5']
                and payload['dataset'] == ds and payload['method'] == 'MrVI_hvg2000'
                and observed_row['runtime']['payload'] == payload
                and set(timing.columns) == {'dataset', 'method', 'time_secs', 'mem_GB'} and len(timing) == 1
                and timing.iloc[0]['dataset'] == ds and timing.iloc[0]['method'] == 'MrVI_hvg2000'
                and json.loads(timing.to_json(orient='records')) == observed_row['timing']['rows'], 'Runtime/timing schema or original-fit identity differs')
        for cost in ('time_secs', 'mem_GB'):
            require(type(payload[cost]) in (int, float) and math.isfinite(payload[cost]) and payload[cost] >= 0
                    and math.isfinite(float(timing.iloc[0][cost])) and float(timing.iloc[0][cost]) >= 0
                    and math.isclose(float(timing.iloc[0][cost]), payload[cost], rel_tol=1e-9, abs_tol=1e-8), 'Runtime/timing costs differ')
        for kind, path, producer, observed_item in (('result', output, 'stage5_mrvi', None),
                ('runtime', runtime_path, 'stage5_mrvi', observed_row['runtime']), ('timing', timing_path, 'stage5_execution_log', observed_row['timing'])):
            item = record(path, original, producer)
            if observed_item is not None:
                require(observed_item['exists'] is True and observed_item['path'] == str(path) and observed_item['md5'] == item['md5']
                        and str(observed_item['size']) == item['size'] and proof(observed_item['sidecar']) == Path(str(path) + '.md5')
                        and proof(observed_item['producer_record']) == Path(item['record_path']), 'Original companion proof changed')
                previous = observed_item['owner']
                require(previous['path'] == str(owner_file(path, original)) and type(previous['exists']) is bool, 'Wrong companion owner namespace')
                if previous['exists']:
                    require(owner(path, original, {'OK', 'FAIL', 'ACTIVE'}) == {key: previous[key] for key in ('path', 'sha256', 'fields')}, 'Original companion owner changed')
                else:
                    require(not owner_file(path, original).exists(), 'An originally absent owner was created')
            relative = 'embeddings/' + ('timing_sources/' if kind == 'timing' else '') + path.name
            inventory.append(dict(dataset=ds, kind=kind, source=str(path), relative=relative, md5=item['md5'], size=item['size'], original_producer=producer))
    destination = safe(str(nas / 'releases' / control['release_id'] / run.name / context), False)
    return control, science, inventory, run, nas, destination


def verified_hr22_scope(control, science, review, receipt_path, run, nas):
    from benchmark_request_scope import AUTHOR_HR
    require(review['status'] == 'APPROVED_SCIENTIFIC_ARTIFACT_REUSE'
            and review['receipt_sha256'] == control['science_receipt']['sha256']
            and review['datasets'] == 1 and review['rds_artifacts'] == 23 and review['timing_feathers'] == 1
            and review['original_gate_status'] == 'FAILED' and review['publication_completed'] is False, 'HR22 review does not approve this receipt')
    require(science['schema'] == 'ecoda_bassez_hr22_science_receipt_v1' and science['status'] == 'NOOP_VALIDATED'
            and science['science_validated'] is True and science['publication_completed'] is False
            and science['original_gate_status'] == 'FAILED' and science['source_commit'] == '129c1a46ee7d282773eb090448ed3d061ac8ecf5'
            and (science['n_samples'], science['n_cells'], science['n_cell_types_high_res']) == (26, 74948, 40)
            and science['counts_pseudocount'] == 0.5, 'Wrong fixed HR22 science audit')
    original = safe(science['original_run_root'])
    require(original.name == 'stage5_bassez_author_hr22_20261109_v1' and original.parent == run.parent and original != run
            and run.parent.parent / '_ecoda_inputs' in receipt_path.parents, 'Wrong HR22 original/receipt namespace')
    for key in ('inventory', 'inspection', 'source_manifest', 'runtime_identity', 'source_identity', 'input_record',
                'input_checksum_proof', 'validation_script', 'ari_source'):
        proof(science[key])
    observed = read_json(science['inventory']['path'])
    require(observed['schema'] == 'ecoda_bassez_hr22_failed_run_inventory_v1' and observed['run_id'] == original.name, 'Wrong original HR22 inventory')
    for item in observed['text_evidence']:
        if item.get('exists') and item.get('sha256') and item['path'].startswith(str(original) + '/'):
            require(not item['truncated'], 'Truncated original proof')
            proof(item)
    row_owner = observed['row_owner']; proof(row_owner); row_fields = fields(row_owner['path'])
    require(row_fields['RUN_ID'] == original.name and row_fields['STATE'] == 'OK' and row_fields['STAGE'] == 'stage5'
            and row_fields['KEY'] == f'derivative/{original.name}/ordinary/Bassez/benchmark_analysis/composition', 'Original HR22 row owner differs')
    source = fields(science['source_manifest']['path']); runtime = fields(science['runtime_identity']['path'])
    require(source['SOURCE_COMMIT'] == science['source_commit'] and source['AUX_ROOT'] == source['SOURCE_ROOT'] + '/aux'
            and digest(Path(source['SOURCE_ROOT']) / 'datasets.json') == source['DATASETS_SHA256'], 'HR22 source/config identity differs')
    aux = science['frozen_aux']
    require(aux['path'] == source['AUX_ROOT'] + '/genes.blocklist.rds' and aux['md5'] == '8cacee7893cbb9d170cc9de87e3bd306'
            and digest(aux['path'], 'md5') == aux['md5'], 'HR22 frozen auxiliary identity differs')
    require(digest(runtime['RUNTIME_MANIFEST']) == runtime['RUNTIME_MANIFEST_SHA256']
            and fields(runtime['RUNTIME_MANIFEST'])['IMAGE_SHA256'] == runtime['RUNTIME_IMAGE_SHA256']
            and safe(runtime['RUNTIME_IMAGE']).stat().st_size == int(runtime['RUNTIME_IMAGE_SIZE']), 'Original HR22 runtime differs')
    require(science['validation_execution']['R_version'] == '4.5.2', 'HR22 audit interpreter differs')
    validator = science['production_validator']; proof(validator); proof(validator['log'])
    require(validator['actual_exit'] == 0 and 'benchmark RDS artifact-list contract OK' in Path(validator['log']['path']).read_text().splitlines()
            and validator['path'] == source['SOURCE_ROOT'] + '/src/5_run_benchmark_methods/validate_benchmark_rds_contract.R', 'HR22 strict CLI evidence differs')
    inspection = read_json(science['inspection']['path']); raw = science['original_accounting']
    require(inspection['state'] == 'FAILED' and inspection['audit']['accounting'] == raw
            and raw['attempted'] is True and raw['query_count'] == 1 and raw['actual_exit'] == 0
            and set(raw['scheduler_ids']) == {'4663753', '4663754', '4663757', '4663759', '4663760'}, 'HR22 failed-controller accounting differs')
    parsed = [line.split('|') for line in raw['stdout'].splitlines()]
    require(raw['stdout'].endswith('\n') and len(parsed) == 5 and len({r[0] for r in parsed}) == 5
            and {r[0] for r in parsed} == set(raw['scheduler_ids']) and all(r[1:] == ['COMPLETED', '0:0'] for r in parsed)
            and raw['rows'] == [dict(line=i, job_id=r[0], state=r[1], exit_code=r[2]) for i, r in enumerate(parsed, 1)], 'All five original completed rows must remain visible')
    canonical = {'anosim_score', 'mod_knn3_score', 'cluster_score', 'sil_score', 'lisi_score'}
    require(set(science['canonical_score_keys']) == canonical and len(science['scores']) == 22
            and {row['method'] for row in science['scores']} == set(AUTHOR_HR)
            and all(all(isinstance(row[k], (int, float)) and math.isfinite(row[k]) for k in canonical) for row in science['scores']), 'Exact five-score HR22 screen differs')
    requested = Path(original / 'manifests/requested_selection.tsv').read_text().splitlines()
    require(len(requested) == 22 and set(requested) == {f'Bassez\tbenchmark_analysis\tcomposition\t{key}' for key in AUTHOR_HR}, 'HR22 original request differs')
    expected_names = ['Bassez_metadata.rds', *(f'Bassez_{key}.rds' for key in AUTHOR_HR)]
    require(len(science['artifacts']) == 23 and [Path(item['path']).name for item in science['artifacts']] == expected_names, 'HR22 artifact coverage/order differs')
    require(len(science['timing']['rows']) == 22 and {row['method'] for row in science['timing']['rows']} == set(AUTHOR_HR), 'HR22 timing coverage differs')
    inventory = []
    items = [(item, 'metadata' if i == 0 else 'result', 'stage5_composition') for i, item in enumerate(science['artifacts'])]
    items.append((science['timing']['proof'], 'timing', 'stage5_execution_log'))
    for item, kind, producer in items:
        path = safe(item['path']); relative = 'embeddings/timing_sources/' + path.name if kind == 'timing' else 'results/' + path.name
        expected = original / ('logs/execution_times_composition_Bassez_selected.feather' if kind == 'timing' else 'artifacts/benchmark/' + relative)
        require(path == expected, 'HR22 source payload outside fixed scope')
        current = record(path, original, producer); proof(item['producer_record']); proof(item['sidecar'])
        require(current['md5'] == item['md5'] and current['size'] == item['size'] and current['record_path'] == item['producer_record']['path'], 'HR22 approved payload changed')
        if kind != 'timing':
            current_owner = owner(path, original, {'OK'}); proof(item['owner'])
            require(current_owner['path'] == item['owner']['path'] and current_owner['sha256'] == item['owner']['sha256'], 'Original HR22 output owner changed')
        inventory.append(dict(dataset='Bassez', kind=kind, source=str(path), relative=relative, md5=item['md5'], size=item['size'], original_producer=producer))
    destination = safe(str(nas / 'releases' / control['release_id'] / run.name / 'benchmark'), False)
    return control, science, inventory, run, nas, destination


def tsv(inventory):
    return ''.join('\t'.join(str(row[k]) for k in ('dataset', 'kind', 'source', 'relative', 'md5', 'size', 'original_producer')) + '\n' for row in inventory)


def write_receipt(destination, value):
    destination = safe(str(destination), False)
    temporary = destination.with_name('.bl-receipt-' + secrets.token_hex(16))
    created = False
    try:
        with temporary.open('x') as stream:
            created = True
            json.dump(value, stream, indent=2); stream.write('\n')
            stream.flush(); os.fsync(stream.fileno())
        os.link(temporary, destination)
    finally:
        if created:
            temporary.unlink()


def copy_payload(args):
    """Copy to a new inode; atomic link installs only the new temporary inode."""
    source = safe(args.source); destination = safe(args.destination, False); boundary = safe(args.boundary)
    if source.name.endswith('.runtime.json'):
        control_path = safe(args.control)
        require(digest(control_path, 'md5') == args.control_md5, 'Runtime-copy control digest differs')
        control = read_json(control_path); selected_run = control_path.parent.parent
        require(control['schema'] == UNC4_SCHEMA and control_path == selected_run / 'manifests/control.json'
                and selected_run.parent.name == '_ecoda_runs' and control['run_id'] == selected_run.name, 'Runtime-copy control is not the fixed run-owned scope')
        allowed = {root / 'artifacts/batch_effect/uncorrected_final/embeddings' /
                   f'{ds}_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json'
                   for root in (selected_run, selected_run.parent / UNC4_RUN) for ds in UNC4_DATASETS}
        require(source in allowed, 'Runtime copy is outside the fixed four-fit scope')
    else:
        require(source.suffix in {'.rds', '.feather'} or source.name == 'nas_checksums.md5', 'Copy source is outside finite BL payload kinds')
    require(boundary in destination.parents and not destination.exists(), 'Copy target exists/outside boundary')
    require(digest(source, 'md5') == args.md5 and source.stat().st_size == int(args.size), 'Copy source differs')
    # Pin the parent chain. No symlink can redirect a temporary or publication.
    fds = [os.open(boundary, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)]
    temporary = '.bl-release-' + secrets.token_hex(16)
    created = False
    try:
        for part in destination.parent.relative_to(boundary).parts:
            fds.append(os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fds[-1]))
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600, dir_fd=fds[-1])
        created = True
        with source.open('rb') as reader, os.fdopen(fd, 'wb') as writer:
            shutil.copyfileobj(reader, writer, 1024 * 1024); writer.flush(); os.fsync(writer.fileno())
        require(digest(source, 'md5') == args.md5, 'Copy source changed')
        observed = safe(str(destination.parent)).stat(); pinned = os.fstat(fds[-1])
        require((observed.st_dev, observed.st_ino) == (pinned.st_dev, pinned.st_ino), 'Copy parent changed')
        os.link(temporary, destination.name, src_dir_fd=fds[-1], dst_dir_fd=fds[-1], follow_symlinks=False)
        os.unlink(temporary, dir_fd=fds[-1]); os.fsync(fds[-1])
        require(digest(destination, 'md5') == args.md5 and destination.stat().st_size == int(args.size) and destination.stat().st_nlink == 1
                and (source.stat().st_dev, source.stat().st_ino) != (destination.stat().st_dev, destination.stat().st_ino), 'Copied bytes/inode differ')
    finally:
        try:
            if created:
                os.unlink(temporary, dir_fd=fds[-1])
        except FileNotFoundError:
            pass
        for fd in reversed(fds):
            os.close(fd)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('plan', 'verify', 'finish', 'validate', 'copy'))
    for key in ('control', 'control-md5', 'run-root', 'nas-root', 'inventory', 'report', 'source', 'destination', 'boundary', 'md5', 'size'):
        parser.add_argument('--' + key)
    args = parser.parse_args()
    try:
        if args.mode == 'copy':
            copy_payload(args); return
        control, science, inventory, run, nas, destination = verified_scope(args)
        copy_producer, report_name, report_schema, report_producer = campaign(control)
        listing = safe(args.inventory, args.mode != 'plan')
        require(listing == run / 'manifests/publication_inventory.tsv', 'Inventory must be run-owned')
        if args.mode == 'plan':
            require(not destination.exists(), 'New NAS release already exists')
            with listing.open('x') as stream:
                stream.write(tsv(inventory))
            print(destination); return
        require(listing.read_text() == tsv(inventory), 'Publication inventory differs')
        if args.mode == 'verify':
            return
        require(args.report == str(run / 'artifacts/evidence' / report_name), 'Report must be run-owned')
        mappings = []
        checksums = []
        for row in inventory:
            copied = run / 'artifacts' / analysis_context(control) / row['relative']; published = destination / row['relative']
            for path in (copied, published):
                item = record(path, run, copy_producer); owned = owner(path, run, {'OK'})
                require(item['md5'] == row['md5'] and item['size'] == row['size'], 'Release payload differs from approved science')
                if path == published:
                    require(path.stat().st_mode & 0o7777 == 0o666 and Path(str(path) + '.md5').stat().st_mode & 0o7777 == 0o666, 'NAS file mode differs')
                    parent = path.parent
                    while parent != nas:
                        require(safe(str(parent)).stat().st_mode & 0o7777 == 0o777, 'NAS ancestor mode differs')
                        parent = parent.parent
                else:
                    require(not path.stat().st_mode & 0o222, 'Scratch release copy is writable')
            checksums.extend((f"{row['md5']}  {row['relative']}\n", f"{digest(str(published)+'.md5', 'md5')}  {row['relative']}.md5\n"))
            mappings.append(dict(**row, scratch=str(copied), nas=str(published), scratch_record=record(copied, run, copy_producer), nas_record=record(published, run, copy_producer),
                                 scratch_owner=owner(copied, run, {'OK'}), nas_owner=owner(published, run, {'OK'})))
        checksum_path = safe(str(destination / 'checksums.md5'))
        require(checksum_path.read_text() == ''.join(checksums) and checksum_path.stat().st_mode & 0o7777 == 0o666, 'Global selected checksum manifest differs')
        allowed = {r['relative'] + suffix for r in inventory for suffix in ('', '.md5')} | {'checksums.md5'}
        actual = set()
        for directory, dirs, files in os.walk(destination, followlinks=False):
            for name in dirs + files:
                safe(str(Path(directory) / name))
            actual.update(str((Path(directory) / name).relative_to(destination)) for name in files)
        require(actual == allowed, 'NAS release contains unselected/missing files')
        source = run / 'manifests/source.manifest'; runtime = run / 'manifests/runtime.identity'
        mounts = [run / 'manifests' / name for name in ('nas_preflight.json', 'nas_postflight.json')]
        mounted = [read_json(path) for path in mounts]
        require(mounted[0]['mount_source'] == mounted[1]['mount_source']
                and all(item['nas_root'] == str(nas) for item in mounted), 'NAS mount evidence differs')
        result = dict(schema=report_schema, status='PUBLISHED_VALIDATED', publication_completed=True,
                      new_model_fit=False, original_gate_status=science['original_gate_status'], original_run_root=science['original_run_root'], run_id=run.name,
                      analysis_root=str(run / 'artifacts' / analysis_context(control)), destination_root=str(destination),
                      control=dict(path=args.control, md5=args.control_md5), science_receipt=control['science_receipt'], review=control['review'],
                      publication_source=dict(path=str(source), sha256=digest(source)), publication_runtime=dict(path=str(runtime), sha256=digest(runtime)),
                      nas_mount_checks=[dict(path=str(path), sha256=digest(path)) for path in mounts],
                      original_accounting=science['original_accounting'], mappings=mappings,
                      global_checksums=dict(path=str(checksum_path), sha256=digest(checksum_path)),
                      completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      policy='New release copies of reviewed bytes; original producers, owners, statuses and failed gate remain unchanged. No scheduler query, model fit, or H5AD read.')
        if control['schema'] in MAIN_SCHEMAS:
            result.update(worker_exit_not_established=True, completion_criterion='VALIDATED_PUBLISHED_ARTIFACTS',
                          publication_scope='SELECTED_METHOD_ARTIFACTS_ONLY', sample_metadata_published=False,
                          sample_metadata_policy='Use separately validated metadata/adoption in the analysis release; no export/TASK provenance is invented.',
                          required_sample_metadata_sources=[] if control['context'] == 'benchmark' else sorted({
                              str(Path(science['original_run_root']) / 'artifacts' / control['context'] / 'metadata' / (row['dataset'] + '_sample_metadata.feather'))
                              for row in inventory}))
        if control['schema'] == UNC4_SCHEMA:
            result.update(worker_exit_not_established=True, completion_criterion='VALIDATED_PUBLISHED_ARTIFACTS', publication_scope='FOUR_MRVI_FITS_ONLY', sample_metadata_published=False,
                          runtime_companion_policy='BYTE_IDENTICAL_ORIGINAL_FIT_EVIDENCE_NOT_REBASED',
                          sample_metadata_policy='Separate validated v1 selected metadata-adoption lineage; no exports or TASK provenance are created.')
        if control['schema'] == REBOOT_SCHEMA:
            result.update(incident=control['incident'], orchestration_outcome='LOGIN_REBOOT_LOSS', wrapper_exit_code=None,
                          policy='New reviewed artifact copies; original PRELAUNCH_STOP, absent wrapper outcome, statuses and owners remain unchanged. No scheduler query, fit, or H5AD read.')
        if args.mode == 'validate':
            existing = read_json(args.report)
            require(set(existing) == set(result) and all(existing[k] == value for k, value in result.items() if k != 'completed_at'), 'Publication receipt differs from current verified release')
            completed = datetime.datetime.fromisoformat(existing['completed_at'])
            require(completed.tzinfo is not None and completed <= datetime.datetime.now(datetime.timezone.utc), 'Invalid completion timestamp')
            record(Path(args.report), run, report_producer); owner(Path(args.report), run, {'OK'})
            require(fields(run / 'status/terminal')['STATE'] == 'NOOP_VALIDATED', 'Publication run is not terminal validated')
            print('FAILED16_UNC4_PUBLICATION_CONTRACT_OK' if control['schema'] == UNC4_SCHEMA else 'MAIN_V2_SCIENCE_PUBLICATION_CONTRACT_OK' if control['schema'] in MAIN_SCHEMAS else 'BL_PUBLICATION_CONTRACT_OK'); return
        write_receipt(args.report, result)
    except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
        parser.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
