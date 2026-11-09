"""Fixed BL10 and Bassez HR22 releases; never read H5ADs or fit models.

Control schema ecoda_bl_publication_control_v1 contains run_id, release_id,
nas_root, science_receipt and review (the latter two are {path, sha256}).
The review is APPROVED_SCIENTIFIC_ARTIFACT_REUSE, binds receipt_sha256 and
declares datasets=10, rds_artifacts=20, timing_feathers=10,
original_gate_status=FAILED, publication_completed=false.
The second fixed schema is ecoda_bassez_hr22_publication_control_v1, with
the same keys and review fields but datasets=1, rds_artifacts=23, timings=1.
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
import shutil
import sys

STAT_KEYS = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
ORIGINAL = 'stage5_benchmark_bl10_20261109_v2'
PRODUCER = 'stage5_bl_science_release_copy'
HR_SCHEMA = 'ecoda_bassez_hr22_publication_control_v1'


def campaign(control):
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


def owner(p, run, states):
    key = re.sub(r'[/,:\t |]', '_', p.name) + '_' + hashlib.sha256(str(p).encode()).hexdigest()[:32]
    if len(key.encode()) > 100:
        key = 'artifact_' + hashlib.sha256(str(p).encode()).hexdigest()[:32]
    path = run.parent.parent / '_ecoda_owners/artifact' / key / 'owner'
    value = fields(path)
    require(value.get('RUN_ID') == run.name and value.get('STAGE') == 'stage5' and value.get('PATH') == str(p)
            and value.get('STATE') in states, 'Wrong artifact owner')
    return dict(path=str(path), sha256=digest(path), fields=value)


def verified_scope(args):
    control_path = safe(args.control)
    require(digest(control_path, 'md5') == args.control_md5, 'Control MD5 differs')
    control = read_json(control_path)
    require(set(control) == {'schema', 'run_id', 'release_id', 'nas_root', 'science_receipt', 'review'}
            and control['schema'] in {'ecoda_bl_publication_control_v1', HR_SCHEMA}, 'Unknown publication control')
    run = safe(args.run_root); nas = safe(args.nas_root)
    require(run.is_dir() and run.parent.name == '_ecoda_runs' and control['run_id'] == run.name
            and control['nas_root'] == str(nas) and nas.is_dir()
            and re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*', control['release_id']), 'Publication roots differ')
    require(control_path == run / 'manifests/control.json', 'Publication control must be run-owned')
    receipt_path = proof(control['science_receipt']); review = read_json(proof(control['review'])); science = read_json(receipt_path)
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
            copied = run / 'artifacts/benchmark' / row['relative']; published = destination / row['relative']
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
                      new_model_fit=False, original_gate_status='FAILED', original_run_root=science['original_run_root'], run_id=run.name,
                      analysis_root=str(run / 'artifacts/benchmark'), destination_root=str(destination),
                      control=dict(path=args.control, md5=args.control_md5), science_receipt=control['science_receipt'], review=control['review'],
                      publication_source=dict(path=str(source), sha256=digest(source)), publication_runtime=dict(path=str(runtime), sha256=digest(runtime)),
                      nas_mount_checks=[dict(path=str(path), sha256=digest(path)) for path in mounts],
                      original_accounting=science['original_accounting'], mappings=mappings,
                      global_checksums=dict(path=str(checksum_path), sha256=digest(checksum_path)),
                      completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      policy='New release copies of reviewed bytes; original producers, owners, statuses and failed gate remain unchanged. No scheduler query, model fit, or H5AD read.')
        if args.mode == 'validate':
            existing = read_json(args.report)
            require(set(existing) == set(result) and all(existing[k] == value for k, value in result.items() if k != 'completed_at'), 'Publication receipt differs from current verified release')
            completed = datetime.datetime.fromisoformat(existing['completed_at'])
            require(completed.tzinfo is not None and completed <= datetime.datetime.now(datetime.timezone.utc), 'Invalid completion timestamp')
            record(Path(args.report), run, report_producer); owner(Path(args.report), run, {'OK'})
            require(fields(run / 'status/terminal')['STATE'] == 'NOOP_VALIDATED', 'Publication run is not terminal validated')
            print('BL_PUBLICATION_CONTRACT_OK'); return
        write_receipt(args.report, result)
    except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
        parser.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
