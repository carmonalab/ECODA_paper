"""Finite BL10 release: verify reviewed bytes; never read H5ADs or fit models.

Control schema ecoda_bl_publication_control_v1 contains run_id, release_id,
nas_root, science_receipt and review (the latter two are {path, sha256}).
The review is APPROVED_SCIENTIFIC_ARTIFACT_REUSE, binds receipt_sha256 and
declares datasets=10, rds_artifacts=20, timing_feathers=10,
original_gate_status=FAILED, publication_completed=false.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import secrets
import shutil
import sys

STAT_KEYS = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
ORIGINAL = 'stage5_benchmark_bl10_20261109_v2'
PRODUCER = 'stage5_bl_science_release_copy'


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
            and control['schema'] == 'ecoda_bl_publication_control_v1', 'Unknown publication control')
    run = safe(args.run_root); nas = safe(args.nas_root)
    require(run.is_dir() and run.parent.name == '_ecoda_runs' and control['run_id'] == run.name
            and control['nas_root'] == str(nas) and nas.is_dir()
            and re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*', control['release_id']), 'Publication roots differ')
    require(control_path == run / 'manifests/control.json', 'Publication control must be run-owned')
    receipt_path = proof(control['science_receipt']); review = read_json(proof(control['review'])); science = read_json(receipt_path)
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
        require(args.report == str(run / 'artifacts/evidence/bl_publication.json'), 'Report must be run-owned')
        mappings = []
        checksums = []
        for row in inventory:
            copied = run / 'artifacts/benchmark' / row['relative']; published = destination / row['relative']
            for path in (copied, published):
                item = record(path, run, PRODUCER); owned = owner(path, run, {'OK'})
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
            mappings.append(dict(**row, scratch=str(copied), nas=str(published), scratch_record=record(copied, run, PRODUCER), nas_record=record(published, run, PRODUCER),
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
        result = dict(schema='ecoda_bl_publication_receipt_v1', status='PUBLISHED_VALIDATED', publication_completed=True,
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
            record(Path(args.report), run, 'bl_science_publication'); owner(Path(args.report), run, {'OK'})
            require(fields(run / 'status/terminal')['STATE'] == 'NOOP_VALIDATED', 'Publication run is not terminal validated')
            print('BL_PUBLICATION_CONTRACT_OK'); return
        write_receipt(args.report, result)
    except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
        parser.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
