#!/bin/bash
# Operational fixtures use real owner/record/checksum helpers and real Feather
# validation. No scheduler or scientific worker is allowed in maintenance.
# --sealed-import-only narrows to the sealed-source tamper regression.
# --foreign-promotion-only skips interruptions and checks fresh-run reuse.
# --terminal-status-only checks authoritative producer terminal-file inventory.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON_BIN:-${ROOT}/.pixi/envs/py-cuda13/bin/python}"
[[ -x "$PYTHON" ]] || { echo 'Configured Pixi Python is required' >&2; exit 1; }
export PYTHONDONTWRITEBYTECODE=1
"$PYTHON" - "$ROOT" "$@" <<'PY'
import ast
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

import pandas as pd

repo = Path(sys.argv[1]).resolve()
entrypoint = repo / 'src/5_run_benchmark_methods/finalize_corrected_outputs.sh'
sealed_import_only = sys.argv[2:] == ['--sealed-import-only']
foreign_promotion_only = sys.argv[2:] == ['--foreign-promotion-only']
terminal_status_only = sys.argv[2:] == ['--terminal-status-only']
if sys.argv[2:] and not (sealed_import_only or foreign_promotion_only or terminal_status_only):
    raise SystemExit('usage: test_corrected_output_promotion.sh [--sealed-import-only|--foreign-promotion-only|--terminal-status-only]')
# Exercise the embedded finite I/O module without executing its CLI main.
body = entrypoint.read_text().split("<<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
ns = {'__name__': 'promotion_fixture'}
exec(compile(body, str(entrypoint), 'exec'), ns)
main_ast = next(node for node in ast.parse(body).body if isinstance(node, ast.FunctionDef) and node.name == 'main')
for name in ('receipt_binds', 'promotion_sources', 'owner', 'original_evidence', 'legacy_runtime_evidence', 'collision_evidence', 'scientific', 'inspect'):
    node = next(node for node in main_ast.body if isinstance(node, ast.FunctionDef) and node.name == name)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])), str(entrypoint), 'exec'), ns)
sys.path.insert(0, str(repo / 'src/utils/py'))
validator = ns['module']('promotion_real_validator', repo / 'src/5_run_benchmark_methods/matrix_artifact_validator.py')
identity_module = ns['module']('promotion_real_identity', repo / 'src/utils/py/h5ad_source_identity.py')
from batch_contract import build_batch_contract_identity


def rejects(operation, message):
    try:
        operation()
    except (ValueError, OSError, subprocess.CalledProcessError):
        return
    raise AssertionError(message)


with tempfile.TemporaryDirectory(prefix='ecoda-promotion-') as temporary:
    base = Path(temporary).resolve()
    metadata_probe = base / 'metadata'
    metadata_probe.write_text('RUN_ID=producer\nFORCE_REASON=\n\nTARGET_METHODS=\n')
    assert ns['fields'](metadata_probe, metadata=True) == {
        'RUN_ID': 'producer', 'FORCE_REASON': '', 'TARGET_METHODS': ''}
    rejects(lambda: ns['fields'](metadata_probe), 'empty checksum/owner evidence accepted')
    for malformed in ('RUN_ID=producer\nRUN_ID=other\n', 'RUN_ID=producer\nmissing_equals\n', '=value\n'):
        metadata_probe.write_text(malformed)
        rejects(lambda: ns['fields'](metadata_probe, metadata=True), 'malformed/duplicate metadata accepted')
    scratch = base / 'scratch'
    runs = scratch / '_ecoda_runs'
    owners = scratch / '_ecoda_owners'
    canonical = scratch / 'batch_effect/corrected_final'
    recovery = canonical / 'recovery_35row'
    promotion = runs / 'cutover_promotion'
    manifests = promotion / 'manifests'
    incoming = base / 'incoming'
    for path in (runs, owners, canonical / 'embeddings', recovery / 'embeddings', manifests, incoming):
        path.mkdir(parents=True, exist_ok=True)
    os.environ.update(HPC_SCRATCH_DIR=str(scratch), ECODA_RUNS_ROOT=str(runs), ECODA_OWNERS_ROOT=str(owners),
                      ECODA_RUN_ROOT=str(promotion), ECODA_RUN_ID='cutover_promotion', PROJECT_ROOT=str(repo),
                      ANALYSIS_ROOT=str(canonical), NAS_TARGET_DIR=str(base / 'nas'),
                      ANALYSIS_NAS_ROOT=str(base / 'nas/batch_effect/corrected_final'),
                      ANALYSIS_VARIANT='corrected_final', ANALYSIS_PASS='corrected', PASS_ARG='corrected')
    stub = base / 'bin'
    stub.mkdir()
    (stub / 'sbatch').write_text('#!/bin/bash\necho UNEXPECTED_MAINTENANCE_SUBMIT >&2\nexit 93\n')
    (stub / 'sbatch').chmod(0o755)
    os.environ['PATH'] = str(stub) + os.pathsep + os.environ['PATH']
    def local_shell(code, *args, capture=True, env=None):
        prelude = 'set -euo pipefail\nsource ' + shlex.quote(str(repo / 'src/utils/bash/ecoda_run_common.sh')) + '\nsource ' + shlex.quote(str(repo / 'src/utils/bash/ecoda_stage5_policy.sh')) + '\n'
        return subprocess.run(['bash', '-c', prelude + code, 'fixture', *map(str, args)], check=True,
                              text=True, stdout=subprocess.PIPE if capture else None, env=env).stdout
    ns['shell'] = local_shell
    ns.update(root=repo, runs=runs, scratch=scratch, manifests=manifests, incoming=incoming,
              run_id='cutover_promotion', source_run='stage5_corrected_all9_ygg_retry5_20260916T044420Z',
              validator=validator, identity_module=identity_module, evidence=set(), transfer_rows={})
    source_run = ns['source_run']
    old = runs / source_run
    (old / 'manifests').mkdir(parents=True)
    (old / 'status/watchdogs').mkdir(parents=True)
    (old / 'metadata').write_text(f'STAGE=stage5\nRUN_ID={source_run}\nSTATE=ACTIVE\nANALYSES=\nFORCE_REASON=\n\nTARGET_METHODS=\n')
    (old / 'status/terminal').write_text(f'STATE=FAIL\nRUN_ID={source_run}\nREASON=fixture\nTIME=fixture\n')
    original_status = (old / 'metadata').read_bytes()
    source_record = {'dataset': 'Breast_cancer', 'view': ns['VIEW'], 'path': str(scratch / 'input.h5ad'),
                     'size': 10, 'md5': 'a' * 32, 'sample_ids': ['s1', 's2']}
    identity_payload = {'schema': identity_module.IDENTITY_SCHEMA, 'entries': [source_record]}
    old_identity = old / 'manifests/source_identity.json'
    old_identity.write_text(json.dumps(identity_payload))
    ns['checksum'](old_identity)
    source_manifest = old / 'manifests/source.manifest'
    frozen = scratch / 'snapshots/original'
    (frozen / 'identity').mkdir(parents=True)
    original_archive = frozen / 'identity/source.tar'
    import tarfile
    with tarfile.open(original_archive, 'w') as archive:
        archive.add(repo / 'datasets.json', arcname='datasets.json')
    (frozen / 'COMPLETE').write_text('COMPLETE\n')
    source_manifest.write_text(''.join(f'{key}={ns["digest"](repo / name, "sha256")}\n' for name, key in
        (('datasets.json', 'DATASETS_SHA256'), ('pixi.toml', 'PIXI_TOML_SHA256'), ('pixi.lock', 'PIXI_LOCK_SHA256'))))
    with source_manifest.open('a') as stream:
        stream.write(f'SOURCE_ARCHIVE_PATH={original_archive}\nSOURCE_ARCHIVE_SHA256={ns["digest"](original_archive, "sha256")}\n')
    original_image = frozen / 'runtime.sif'
    original_image.write_bytes(b'finite runtime identity fixture')
    original_runtime_manifest = frozen / 'runtime.manifest'
    original_runtime_manifest.write_text('IMAGE_SHA256=' + ns['digest'](original_image, 'sha256') + '\n')
    runtime_identity = (f'RUNTIME_IMAGE={original_image}\nRUNTIME_MANIFEST={original_runtime_manifest}\n'
                        f'RUNTIME_IMAGE_SHA256={ns["digest"](original_image, "sha256")}\n'
                        f'RUNTIME_MANIFEST_SHA256={ns["digest"](original_runtime_manifest, "sha256")}\n')
    (old / 'manifests/runtime.identity').write_text(runtime_identity)
    (manifests / 'runtime.identity').write_text(runtime_identity)
    ns['checked_runtimes'] = set()
    matrix = old / 'manifests/matrix_batch_effect_corrected_mrvi_default_gpu.tsv'
    matrix.write_text('Breast_cancer\tbatch_effect_corrected\tmrvi\n')
    ns['checksum'](matrix)
    watchdog = old / 'status/watchdogs/batch_effect_corrected__mrvi__default_gpu.status'
    watchdog.write_text('STATE=OK\nLABEL=batch_effect_corrected__mrvi__default_gpu\nSCHEDULER_ID=1\nSCHEDULER_ID=2\n')
    ns['current_records'] = {('Breast_cancer', ns['VIEW']): source_record}
    ns['identity'] = old_identity
    identity = build_batch_contract_identity(['assay', 'suspension_dissociation_time'], sample_column='Sample',
                    method_id='MrVI', model_id='mrvi_composite_v1')
    frame = pd.DataFrame([[0.0, 1.0], [1.0, 0.0]], index=['s1', 's2'], columns=['s1', 's2'])
    def artifact(path, state='OK', values=frame):
        values.to_feather(path)
        ns['checksum'](path)
        runtime = Path(str(path) + '.runtime.json')
        runtime.write_text(json.dumps({'schema_version': 1, 'artifact_path': str(path),
            'artifact_md5': ns['digest'](path), 'dataset': 'Breast_cancer', 'method': 'MrVI_hvg2000',
            'time_secs': 12.5, 'mem_GB': 1.5, 'batch_contract': identity}))
        ns['checksum'](runtime)
        for item in (path, runtime):
            local_shell('ecoda_artifact_owner_acquire "$1" stage5 "$2"; ecoda_write_artifact_record "$1" stage5_mrvi "$2"; ecoda_artifact_owner_release "$1" "$3" fixture', item, source_run, state)
        return runtime
    name = 'Breast_cancer_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather'
    canonical_path = canonical / 'embeddings' / name
    candidate = recovery / 'embeddings' / name
    canonical_runtime = artifact(canonical_path)
    candidate_runtime = artifact(candidate, 'FAIL', frame * 2)
    candidate_before = candidate.read_bytes()
    evidence_before = {path: path.read_bytes() for path in old.rglob('*') if path.is_file()}
    # Valid canonical wins independently; different recovery bytes are irrelevant.
    run, states, _ = ns['inspect']([canonical_path], 'Breast_cancer', 'mrvi')
    assert run == source_run and states == ['OK', 'OK']
    # Older transferred producer matrices have no original sidecar. Accept only
    # their exact pinned transfer bytes, never fabricate evidence or bypass a bad sidecar.
    matrix_sidecar = Path(str(matrix) + '.md5')
    matrix_checksum = matrix_sidecar.read_bytes()
    matrix_bytes = matrix.read_bytes()
    matrix_relative = str(matrix.relative_to(scratch))
    ns['transfer_rows'][matrix_relative] = (matrix.stat().st_size, ns['digest'](matrix))
    matrix_sidecar.unlink()
    ns['inspect']([canonical_path], 'Breast_cancer', 'mrvi')
    assert not matrix_sidecar.exists(), 'historical checksum was fabricated'
    matrix.write_bytes(matrix_bytes + b'other\tbatch_effect_corrected\tmrvi\n')
    rejects(lambda: ns['inspect']([canonical_path], 'Breast_cancer', 'mrvi'), 'changed unchecksummed producer matrix accepted')
    matrix.write_bytes(matrix_bytes)
    del ns['transfer_rows'][matrix_relative]
    rejects(lambda: ns['inspect']([canonical_path], 'Breast_cancer', 'mrvi'), 'unbound producer matrix accepted')
    ns['transfer_rows'][matrix_relative] = (matrix.stat().st_size, ns['digest'](matrix))
    matrix_sidecar.write_bytes(matrix_checksum.replace(b'MD5=', b'MD5=bad'))
    rejects(lambda: ns['inspect']([canonical_path], 'Breast_cancer', 'mrvi'), 'bad existing sidecar bypassed via transfer evidence')
    matrix_sidecar.write_bytes(matrix_checksum)
    runtime_owner = Path(local_shell('ecoda_artifact_owner_dir "$1"', canonical_runtime).strip())
    owner_bytes = (runtime_owner / 'owner').read_bytes()
    (runtime_owner / 'owner').unlink()
    runtime_owner.rmdir()
    runtime_record = Path(local_shell('ecoda_artifact_record_path "$1" "$2"', canonical_runtime, source_run).strip())
    runtime_proof = [canonical_runtime, Path(str(canonical_runtime) + '.md5'), runtime_record]
    rejects(lambda: ns['collision_evidence'](canonical_runtime, 'Breast_cancer', 'mrvi'), 'unbound unowned runtime archived')
    for path in runtime_proof:
        ns['transfer_rows'][str(path.relative_to(scratch))] = (path.stat().st_size, ns['digest'](path))
    archived = ns['collision_evidence'](canonical_runtime, 'Breast_cancer', 'mrvi')
    assert set(runtime_proof).issubset(archived) and not runtime_owner.exists()
    rejects(lambda: ns['inspect']([canonical_path], 'Breast_cancer', 'mrvi'), 'archive-only proof authorized runtime reuse')
    runtime_bytes = canonical_runtime.read_bytes()
    ns['atomic'](canonical_runtime, runtime_bytes + b' ')
    rejects(lambda: ns['collision_evidence'](canonical_runtime, 'Breast_cancer', 'mrvi'), 'changed unowned runtime archived')
    ns['atomic'](canonical_runtime, runtime_bytes)
    payload_owner = Path(local_shell('ecoda_artifact_owner_dir "$1"', canonical_path).strip()) / 'owner'
    payload_owner_bytes = payload_owner.read_bytes()
    payload_owner.write_bytes(payload_owner_bytes.replace(b'STATE=OK', b'STATE=ACTIVE'))
    rejects(lambda: ns['collision_evidence'](canonical_runtime, 'Breast_cancer', 'mrvi'), 'active paired payload authorized archival replacement')
    payload_owner.write_bytes(payload_owner_bytes)
    runtime_owner.mkdir()
    (runtime_owner / 'owner').write_bytes(owner_bytes)
    ns['inspect']([candidate], 'Breast_cancer', 'mrvi', True)
    recovery_runtime_owner = Path(local_shell('ecoda_artifact_owner_dir "$1"', candidate_runtime).strip())
    recovery_runtime_owner_bytes = (recovery_runtime_owner / 'owner').read_bytes()
    (recovery_runtime_owner / 'owner').unlink()
    recovery_runtime_owner.rmdir()
    recovery_record = Path(local_shell('ecoda_artifact_record_path "$1" "$2"', candidate_runtime, source_run).strip())
    for path in (candidate_runtime, Path(str(candidate_runtime) + '.md5'), recovery_record):
        ns['transfer_rows'][str(path.relative_to(scratch))] = (path.stat().st_size, ns['digest'](path))
    recovered_run, recovered_states, recovered_proof = ns['inspect']([candidate], 'Breast_cancer', 'mrvi', True)
    assert recovered_run == source_run and recovered_states == ['FAIL']
    assert recovery_record in recovered_proof and not recovery_runtime_owner.exists()
    rejects(lambda: ns['inspect']([candidate], 'Breast_cancer', 'mrvi'), 'legacy recovery proof weakened ordinary runtime ownership')
    del ns['transfer_rows'][str(recovery_record.relative_to(scratch))]
    rejects(lambda: ns['inspect']([candidate], 'Breast_cancer', 'mrvi', True), 'unbound legacy recovery runtime accepted')
    ns['transfer_rows'][str(recovery_record.relative_to(scratch))] = (recovery_record.stat().st_size, ns['digest'](recovery_record))
    watchdog_before = watchdog.read_bytes()
    watchdog.write_bytes(watchdog_before.replace(b'STATE=OK', b'STATE=FAIL'))
    rejects(lambda: ns['inspect']([candidate], 'Breast_cancer', 'mrvi', True), 'unowned recovery runtime bypassed failed producer watchdog')
    watchdog.write_bytes(watchdog_before)
    recovery_runtime_owner.mkdir()
    (recovery_runtime_owner / 'owner').write_bytes(recovery_runtime_owner_bytes)
    assert candidate.read_bytes() == candidate_before
    assert (old / 'metadata').read_bytes() == original_status
    watchdog_before = watchdog.read_bytes()
    watchdog.write_bytes(watchdog_before.replace(b'STATE=OK', b'STATE=FAIL'))
    rejects(lambda: ns['inspect']([candidate], 'Breast_cancer', 'mrvi', True), 'failed producer watchdog authorized reuse')
    watchdog.write_bytes(watchdog_before)
    # Same-byte destination and path-only JSON rebasing, strict new sidecars.
    imported = canonical / 'embeddings/imported.feather'
    ns['copy_bytes'](candidate, imported)
    ns['checksum'](imported)
    imported_runtime = Path(str(imported) + '.runtime.json')
    ns['rebase_runtime'](candidate_runtime, imported_runtime, imported)
    ns['checksum'](imported_runtime)
    before = json.loads(candidate_runtime.read_text())
    after = json.loads(imported_runtime.read_text())
    assert imported.read_bytes() == candidate_before
    assert set(before) == set(after)
    assert {key for key in before if before[key] != after[key]} == {'artifact_path'}
    local_shell('ecoda_validate_checksum "$1"; ecoda_validate_checksum "$2"', imported, imported_runtime)
    # Imported source-owner evidence is not installed into global ownership.
    source_owner_dir = Path(local_shell('ecoda_artifact_owner_dir "$1"', candidate).strip())
    original_owner = source_owner_dir / 'owner'
    relative = str(original_owner.relative_to(scratch))
    returned_owner = incoming / relative
    returned_owner.parent.mkdir(parents=True)
    returned_owner.write_bytes(original_owner.read_bytes())
    ns['transfer_rows'][relative] = (returned_owner.stat().st_size, ns['digest'](returned_owner))
    original_owner.unlink()
    source_owner_dir.rmdir()
    imported_owner_dir, data = ns['owner'](candidate, True)
    assert imported_owner_dir == returned_owner.parent and data['STATE'] == 'FAIL'
    assert not source_owner_dir.exists()
    rejects(lambda: ns['owner'](candidate, False), 'canonical missing provenance accepted')
    returned_owner.write_text(returned_owner.read_text().replace('STATE=FAIL', 'STATE=ACTIVE'))
    ns['transfer_rows'][relative] = (returned_owner.stat().st_size, ns['digest'](returned_owner))
    rejects(lambda: ns['owner'](candidate, True), 'ACTIVE source ownership accepted')
    # A malformed or wrong-runtime JSON blocks despite valid distances/checksum.
    canonical_runtime.chmod(0o644)
    original_runtime = canonical_runtime.read_bytes()
    payload = json.loads(original_runtime)
    payload['batch_contract']['batch_keys'] = ['sequencing_platform']
    canonical_runtime.write_text(json.dumps(payload))
    ns['checksum'](canonical_runtime)
    rejects(lambda: ns['scientific']([canonical_path], 'Breast_cancer', 'mrvi', source_run), 'stale Breast keys accepted')
    canonical_runtime.write_bytes(original_runtime)
    ns['checksum'](canonical_runtime)
    # Scientific fixture genuinely decodes Feather and rejects wrong sample order.
    rejects(lambda: validator.require_nonempty([canonical_path], 'order fixture', ['s2', 's1'],
        expected_batch_contract=identity, require_runtime_batch_contract=True, require_corrected_summary=False),
        'sample order mismatch accepted')
    # Archive-before-reclaim: exact historical bytes survive terminal transition.
    old_owner = Path(local_shell('ecoda_artifact_owner_dir "$1"', canonical_path).strip()) / 'owner'
    canonical_before = canonical_path.read_bytes()
    archive = base / 'nas/before/owner'
    ns['copy_bytes'](old_owner, archive)
    archived_owner = archive.read_bytes()
    local_shell('ecoda_artifact_owner_acquire "$1" stage5 cutover_promotion 0 1 1; ecoda_artifact_owner_release "$1" FAIL interrupted', canonical_path)
    assert archive.read_bytes() == archived_owner
    local_shell('ecoda_artifact_owner_acquire "$1" stage5 cutover_promotion 0 1 1; ecoda_artifact_owner_release "$1" OK resumed', canonical_path)
    assert canonical_path.read_bytes() == canonical_before
    # Path attacks, malformed evidence, atomic replacement and repeat copies.
    alias = base / 'alias'
    alias.symlink_to(scratch, target_is_directory=True)
    rejects(lambda: ns['safe_path'](str(alias / 'payload')), 'symlink ancestor accepted')
    rejects(lambda: ns['safe_path'](str(scratch) + '/../escape'), 'path escape accepted')
    malformed = base / 'malformed'
    malformed.write_text('STATE=OK\nSTATE=FAIL\n')
    rejects(lambda: ns['fields'](malformed), 'duplicate metadata accepted')
    target = base / 'atomic'
    ns['atomic'](target, b'first')
    ns['atomic'](target, b'second')
    assert target.read_bytes() == b'second'
    ns['copy_bytes'](candidate, imported)
    assert imported.read_bytes() == candidate_before
    for path, content in evidence_before.items():
        assert path.read_bytes() == content, f'historical producer evidence changed: {path}'
    # Pure classifier computes dependency closure and hard-blocks any Python row.
    dispatch = manifests / 'full63.tsv'
    dispatch.write_text(''.join(f'{ds}\t{ns["VIEW"]}\t{method}\n' for ds in ns['DATASETS'] for method in ns['METHODS']))
    pending = manifests / 'pending.tsv'
    classify = 'source ' + shlex.quote(str(repo / 'src/5_run_benchmark_methods/stage5_pending.sh')) + '''
stage5_method_is_forced() { return 1; }
benchmark_selected_artifacts_valid() { case "$3" in mrvi|pilot|qot) return 0;; *) return 1;; esac; }
stage5_prepare_pseudobulk_valid() { [[ "$1" != Alzheimer ]]; }
benchmark_artifacts_for() { ARTIFACT_PATHS=(); }
stage5_classify_pending_rows "$1" "$2"
'''
    local_shell(classify, dispatch, pending)
    rows = [line.split('\t') for line in pending.read_text().splitlines()]
    assert len(rows) == 28 and rows[6] == ['Alzheimer', ns['VIEW'], 'prepare_pseudobulk']
    assert not any(row[2] in ns['PYTHON_METHODS'] for row in rows)
    rejects(lambda: local_shell(classify.replace('mrvi|pilot|qot', 'mrvi|pilot'), dispatch, pending), 'excluded QOT compute accepted')
    # CLI rejects invalid maintenance invocations before sourcing configuration.
    for arguments in ([], ['--inventory-only', '--execute'], ['--force'], ['--execute', '--source-run', 'other']):
        result = subprocess.run(['bash', str(entrypoint), *arguments], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert result.returncode != 0
# The operational seam is the complete Bash CLI, not an extracted function.
# Only the external runtime and NAS boundary are replaced in the private tree.
# The dispatcher, policy, owners, checksums, records and Python consumers remain
# real. The dispatcher is stopped at its runtime boundary before any sbatch.
import shutil
import h5py
import numpy as np

metadata_reader = ns['module']('promotion_fixture_metadata', repo / 'src/utils/py/export_h5ad_sample_metadata.py')
inventory_fixture = repo / '.gate/corrected_return_20261021T130330Z/inventory.tsv'
assert inventory_fixture.is_file(), 'the pinned finite return inventory is required'
assert hashlib.sha256(inventory_fixture.read_bytes()).hexdigest() == '2ec2cb4becda72b07453816acf3e55d804b631bfc4661217c580723c50bc6ea1'

for zero_residual in ((False,) if sealed_import_only or foreign_promotion_only or terminal_status_only else (False, True)):
    with tempfile.TemporaryDirectory(prefix='ecoda-promotion-cli-') as temporary:
        base = Path(temporary).resolve()
        tree = base / 'source'
        shutil.copytree(repo / 'src', tree / 'src')
        for name in ('datasets.json', 'pixi.toml', 'pixi.lock'):
            shutil.copyfile(repo / name, tree / name)
        scratch = base / 'ECODA_paper'
        runs = scratch / '_ecoda_runs'
        canonical = scratch / 'batch_effect/corrected_final'
        recovery = canonical / 'recovery_35row'
        nas = base / 'nas'
        compute = 'fixture_zero' if zero_residual else 'fixture_pending'
        promotion = runs / (compute + '_promotion')
        old = runs / source_run
        incoming = base / '_ecoda_imports/corrected_return_20261021T130330Z'
        for path in (runs, nas, incoming, old / 'manifests', old / 'logs', old / 'status/watchdogs'):
            path.mkdir(parents=True, exist_ok=True)
        selection = base / 'selection.tsv'
        selection.write_text(''.join(f'{ds}\t{ns["VIEW"]}\t{ns["VIEW"]}\n' for ds in ns['DATASETS']))
        matrix = base / 'matrix.tsv'
        matrix.write_text(''.join(f'{ds}\t{ns["VIEW"]}\t{method}\n' for ds in ns['DATASETS'] for method in ns['METHODS']))
        fixture_env = os.environ.copy()
        fixture_env.update(PROJECT_ROOT=str(tree), ECODA_SOURCE_ROOT=str(tree),
            HPC_SCRATCH_DIR=str(scratch), ECODA_RUNS_ROOT=str(runs),
            ECODA_OWNERS_ROOT=str(scratch / '_ecoda_owners'), ECODA_RUN_ID=compute,
            NAS_TARGET_DIR=str(nas), ECODA_NAS_TARGET_DIR=str(nas), ECODA_RUNTIME_MODE='apptainer',
            ECODA_SOURCE_SNAPSHOT_REQUIRED='1', FINALIZE_SCRIPT_DIR=str(tree / 'src/5_run_benchmark_methods'),
            FIXTURE_COMPUTE_ID=compute, FIXTURE_BASE=str(base), PYTHON_BIN=sys.executable)
        fixture_env.pop('ECODA_RUN_ROOT', None)
        fixture_env.update(ANALYSIS_ROOT=str(canonical), ANALYSIS_NAS_ROOT=str(nas / 'batch_effect/corrected_final'),
            ANALYSIS_VARIANT='corrected_final', ANALYSIS_PASS='corrected', PASS_ARG='corrected')
        def fixture_shell(code, *args):
            return subprocess.run(['bash', '-c', 'set -euo pipefail\nsource ' +
                shlex.quote(str(tree / 'src/utils/bash/ecoda_run_common.sh')) + '\nsource ' +
                shlex.quote(str(tree / 'src/utils/bash/ecoda_stage5_policy.sh')) + '\n' + code,
                'fixture', *map(str, args)], env=fixture_env, check=True, text=True, stdout=subprocess.PIPE).stdout
        def fixture_checksum(path):
            fixture_shell('ecoda_write_checksum "$1"', path)
        for path in (selection, matrix):
            fixture_checksum(path)
        for ds in ns['DATASETS']:
            entry, spec = metadata_reader._config_entry(tree / 'datasets.json', ds, ns['VIEW'])
            sample, required, columns = metadata_reader.requested_columns(tree / 'datasets.json', ds, entry, spec, batch_view=True)
            label = metadata_reader._dataset_columns(entry)['label']
            keys, policy = metadata_reader._resolved_batch_configuration(entry, spec, label, dataset=ds)
            obs = {column: ['1', '1', '2', '2'] for column in set(required) | set(keys)}
            obs['Sample'] = ['s1', 's1', 's2', 's2']
            source = identity_module.resolve_h5ad_path(scratch, tree / 'datasets.json', ds, ns['VIEW'])
            source.parent.mkdir(parents=True)
            with h5py.File(source, 'w') as handle:
                group = handle.create_group('obs')
                group.attrs['encoding-type'] = 'dataframe'
                group.attrs['_index'] = '_index'
                group.create_dataset('_index', data=np.asarray(['c1', 'c2', 'c3', 'c4'], dtype=h5py.string_dtype()))
                for column, values in obs.items():
                    group.create_dataset(column, data=np.asarray(values, dtype=h5py.string_dtype()))
            fixture_checksum(source)
            frame = metadata_reader.read_obs_metadata(source, sample, required, columns,
                metadata_reader.READ_CHUNK_SIZE, batch_keys=keys, batch_metadata_policy=policy, biological_column=label)
            frame = frame.loc[:, metadata_reader._ordered_metadata_columns(tree / 'datasets.json', ds, entry, spec, frame.columns)]
            target = recovery / 'metadata' / (ds + '_sample_metadata.feather')
            target.parent.mkdir(parents=True, exist_ok=True)
            frame.to_feather(target)
            fixture_checksum(target)
        input_sources = [identity_module.resolve_h5ad_path(scratch, tree / 'datasets.json', ds, ns['VIEW'])
            for ds in ns['DATASETS']]
        input_sources = [item for path in input_sources for item in (path, Path(str(path) + '.md5'))]
        (base / 'canonical_inputs.txt').write_text(''.join(str(path.relative_to(scratch)) + '\n' for path in input_sources))
        # One reused input, one conflicting predecessor, and unrelated older
        # derivatives exercise normal NAS discovery without broad mutation.
        matching_targets = []
        for path in input_sources[:2]:
            target = nas / path.relative_to(scratch)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
            matching_targets.append(target)
        matching_stats = {path: (path.stat().st_ino, path.stat().st_mtime_ns) for path in matching_targets}
        # Equivalent strict fields do not imply identical predecessor bytes:
        # publish-inputs must archive this reordered sidecar before replacement.
        sidecar_payload = nas / input_sources[4].relative_to(scratch)
        sidecar_payload.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(input_sources[4], sidecar_payload)
        sidecar_conflict = nas / input_sources[5].relative_to(scratch)
        sidecar_conflict.write_text('\n'.join(reversed(input_sources[5].read_text().splitlines())) + '\n')
        assert ns['fields'](sidecar_conflict) == ns['fields'](input_sources[5])
        assert sidecar_conflict.read_bytes() != input_sources[5].read_bytes()
        matching_targets.append(sidecar_payload)
        matching_stats[sidecar_payload] = (sidecar_payload.stat().st_ino, sidecar_payload.stat().st_mtime_ns)
        input_conflict = nas / input_sources[2].relative_to(scratch)
        input_conflict.parent.mkdir(parents=True, exist_ok=True)
        input_conflict.write_bytes(b'historical input predecessor')
        unrelated = input_conflict.parent / 'older_derivative.h5ad'
        unrelated.write_bytes(b'preserve unrelated NAS derivative')
        unrelated.chmod(0o640)
        unrelated_before = (unrelated.read_bytes(), unrelated.stat().st_mode, unrelated.stat().st_ino, unrelated.stat().st_mtime_ns)
        historical_identity = old / 'manifests/source_identity.json'
        historical_identity.write_text(json.dumps(identity_module.build_source_identity(selection, scratch, tree / 'datasets.json')))
        fixture_checksum(historical_identity)
        frozen = base / 'snapshot'
        (frozen / 'identity').mkdir(parents=True)
        (frozen / 'COMPLETE').write_text('COMPLETE\n')
        archive = frozen / 'identity/source.tar'
        with tarfile.open(archive, 'w') as stream:
            stream.add(tree / 'datasets.json', arcname='datasets.json')
        image = frozen / 'runtime.sif'
        image.write_bytes(b'fixture external runtime')
        runtime_manifest = frozen / 'runtime.manifest'
        runtime_manifest.write_text('IMAGE_SHA256=' + ns['digest'](image, 'sha256') + '\n')
        runtime_identity = base / 'runtime.identity'
        runtime_identity.write_text(f'RUNTIME_IMAGE={image}\nRUNTIME_MANIFEST={runtime_manifest}\n'
            f'RUNTIME_IMAGE_SHA256={ns["digest"](image, "sha256")}\nRUNTIME_MANIFEST_SHA256={ns["digest"](runtime_manifest, "sha256")}\n')
        shutil.copyfile(runtime_identity, old / 'manifests/runtime.identity')
        source_manifest = frozen / 'identity/source.manifest'
        source_manifest.write_text(f'SOURCE_COMMIT=7c3eb18ce258fa02552f85b617ce2e22c40a8b41\nSOURCE_ROOT={tree}\nFORMAT=1\n'
            f'SOURCE_ARCHIVE_PATH={archive}\nSOURCE_ARCHIVE_SHA256={ns["digest"](archive, "sha256")}\n' +
            ''.join(f'{key}={ns["digest"](tree / name, "sha256")}\n' for name, key in
            (('datasets.json', 'DATASETS_SHA256'), ('pixi.toml', 'PIXI_TOML_SHA256'), ('pixi.lock', 'PIXI_LOCK_SHA256'))))
        shutil.copyfile(source_manifest, old / 'manifests/source.manifest')
        (old / 'metadata').write_text(f'RUN_ID={source_run}\nSTAGE=stage5\nSTATE=ACTIVE\nANALYSES=\nFORCE_REASON=\n\nTARGET_METHODS=\n')
        (old / 'status/terminal').write_text(f'STATE=FAIL\nRUN_ID={source_run}\nREASON=fixture\nTIME=fixture\n')
        fixture_env['ECODA_SOURCE_MANIFEST'] = str(source_manifest)
        evidence_manifest = base / 'inventory.tsv'
        shutil.copyfile(inventory_fixture, evidence_manifest)
        evidence_receipt = base / 'verified.json'
        evidence_receipt.write_text(json.dumps(dict(completed_at='2026-10-21T13:10:00Z',
            incoming_root=str(incoming), transfer_id='corrected_return_20261021T130330Z',
            selected_count=1253, selected_bytes=216609447099,
            inventory_sha256='2ec2cb4becda72b07453816acf3e55d804b631bfc4661217c580723c50bc6ea1',
            files_sha256='1d9d82c6f6f1e5b741fd1b8be24afe7c5c34d9569664a9fee3d2349dbdbe2dbf',
            transfer_script_sha256='972594ac79a113c13235445d9c5b50041b7e34374339721928156fb81320d28e',
            verifier_sha256='429831ce5f739e5d420c7ed0188575945af15366a2ba26285a5f088a07aef4d2',
            verification='all selected regular files matched inventory size and MD5')))
        fixture_env.update(ECODA_CORRECTED_SOURCE_EVIDENCE_MANIFEST=str(evidence_manifest),
            ECODA_CORRECTED_SOURCE_EVIDENCE_RECEIPT=str(evidence_receipt))
        original_payloads = {}
        collisions = {}
        for method in ns['METHODS']:
            watchdog_matrix = old / 'manifests' / f'matrix_{ns["VIEW"]}_{method}.tsv'
            watchdog_matrix.write_text(''.join(f'{ds}\t{ns["VIEW"]}\t{method}\n' for ds in ns['DATASETS']))
            fixture_checksum(watchdog_matrix)
            label = ns['VIEW'] + '__' + method
            (old / 'status/watchdogs' / (label + '.status')).write_text(f'STATE=OK\nLABEL={label}\n')
            for ds in ns['DATASETS']:
                output = fixture_shell('_ecoda_stage5_artifacts_for "$1" "$2" "$3"; printf "%s\\n" "${ECODA_BENCHMARK_ARTIFACTS[@]}"', ds, ns['VIEW'], method)
                paths = [Path(value) for value in output.splitlines()]
                if method not in ns['PYTHON_METHODS'] and not zero_residual:
                    continue
                for destination in paths:
                    source = recovery / destination.relative_to(canonical) if method in ns['PYTHON_METHODS'] else destination
                    source.parent.mkdir(parents=True, exist_ok=True)
                    if method in ns['PYTHON_METHODS']:
                        pd.DataFrame([[0., 1.], [1., 0.]], index=['s1', 's2'], columns=['s1', 's2']).to_feather(source)
                        batch = validator._build_corrected_batch_contract(
                            validator.read_datasets_json(tree / 'datasets.json', view=ns['VIEW']), ds, ns['VIEW'], method)
                        runtime = Path(str(source) + '.runtime.json')
                        runtime.write_text(json.dumps(dict(schema_version=1, artifact_path=str(source), artifact_md5=ns['digest'](source),
                            dataset=ds, method={'mrvi': 'MrVI_hvg2000', 'pilot': 'PILOT_hvg2000_highres', 'qot': 'QOT_hvg2000_highres'}[method],
                            time_secs=1., mem_GB=1., batch_contract=batch)))
                        companions = [source, runtime]
                        original_payloads[destination] = source.read_bytes()
                    else:
                        source.write_bytes(b'fixture external R runtime artifact')
                        companions = [source]
                    for item in companions:
                        fixture_checksum(item)
                        fixture_shell('ecoda_artifact_owner_acquire "$1" stage5 "$2"; ecoda_write_artifact_record "$1" "$3" "$2"; ecoda_artifact_owner_release "$1" OK fixture',
                            item, source_run, ns['producer'](method))
                log = old / 'logs' / f'execution_times_batch_effect_corrected_final_{method}_{ds}.feather'
                pd.DataFrame(dict(dataset=[ds], method=[method], time_secs=[1.], mem_GB=[1.])).to_feather(log)
                fixture_checksum(log)
                fixture_shell('ecoda_write_artifact_record "$1" stage5_execution_log "$2"', log, source_run)
        # A genuine valid-checksum/record collision fails scientific validation;
        # its exact owner and payload must be archived before replacement.
        collision = next(iter(original_payloads))
        collision.parent.mkdir(parents=True, exist_ok=True)
        collision.write_bytes(b'invalid feather predecessor')
        fixture_checksum(collision)
        fixture_shell('ecoda_artifact_owner_acquire "$1" stage5 "$2"; ecoda_write_artifact_record "$1" stage5_mrvi "$2"; ecoda_artifact_owner_release "$1" FAIL fixture', collision, source_run)
        collision_owner = Path(fixture_shell('ecoda_artifact_owner_dir "$1"', collision).strip()) / 'owner'
        collisions = {collision: collision.read_bytes(), collision_owner: collision_owner.read_bytes()}
        if not zero_residual:
            cache = Path(fixture_shell('_ecoda_stage5_artifacts_for Joanito "$1" prepare_pseudobulk; printf "%s" "${ECODA_BENCHMARK_ARTIFACTS[0]}"', ns['VIEW']))
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_bytes(b'outdated composite-model pseudobulk cache')
            fixture_checksum(cache)
            fixture_shell('ecoda_artifact_owner_acquire "$1" stage5 "$2"; ecoda_write_artifact_record "$1" stage5_prepare_pseudobulk_hvg2000 "$2"; ecoda_artifact_owner_release "$1" OK fixture', cache, source_run)
            (old / 'status/watchdogs' / (ns['VIEW'] + '__prepare_pseudobulk.status')).write_text(
                'STATE=FAIL\nLABEL=' + ns['VIEW'] + '__prepare_pseudobulk\n')
            collisions[cache] = cache.read_bytes()
        collisions[input_conflict] = input_conflict.read_bytes()
        collisions[sidecar_conflict] = sidecar_conflict.read_bytes()
        original_evidence = {path: path.read_bytes() for path in old.rglob('*') if path.is_file()}
        rscript = base / 'Rscript'
        rscript.write_text('#!/bin/bash\nset -euo pipefail\nwhile [[ $# -gt 0 ]]; do\n'
            'case "$1" in --artifact) [[ "$(cat "$2")" == "fixture external R runtime artifact" ]] || exit 4; shift 2;;'
            ' --output) printf "{}\\n" > "$2"; shift 2;; *) shift;; esac\ndone\n')
        rscript.chmod(0o755)
        host_prefix = base / 'runtime'
        (host_prefix / 'bin').mkdir(parents=True)
        (host_prefix / 'bin/python').symlink_to(sys.executable)
        (host_prefix / 'bin/Rscript').symlink_to(rscript)
        fixture_env['ECODA_HOST_ENV_PREFIX'] = str(host_prefix)
        with (tree / 'src/slurm_config.sh').open('a') as stream:
            stream.write('\nexport PYTHON_BIN=' + shlex.quote(sys.executable) +
                '\nexport PIXI_RSCRIPT=' + shlex.quote(str(rscript)) + '\n')
        real_runtime = tree / 'src/utils/bash/ecoda_runtime_real.sh'
        shutil.copyfile(tree / 'src/utils/bash/ecoda_runtime.sh', real_runtime)
        (tree / 'src/utils/bash/ecoda_runtime.sh').write_text('source ' + shlex.quote(str(real_runtime)) + '\n' + '''#!/bin/bash
if [[ "${1:-}" == "$HPC_SCRATCH_DIR/batch_effect/corrected_final/metadata/Joanito_sample_metadata.feather" &&
      -f "$1" && ! -f "$1.md5" && ! -e "$FIXTURE_BASE/metadata_interrupted" ]]; then
  touch "$FIXTURE_BASE/metadata_interrupted"
  return 74
fi
if [[ "${1:-}" == "$HPC_SCRATCH_DIR/batch_effect/corrected_final/embeddings/Joanito_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather" &&
      ! -e "$FIXTURE_BASE/owner_finalization_interrupted" ]]; then
  fixture_payload_record="$(ecoda_artifact_record_path "$1" "$ECODA_RUN_ID")"
  fixture_runtime_record="$(ecoda_artifact_record_path "$1.runtime.json" "$ECODA_RUN_ID")"
  if [[ -f "$fixture_payload_record" && -f "$fixture_runtime_record" ]]; then
    touch "$FIXTURE_BASE/owner_finalization_interrupted"
    return 75
  fi
fi
ecoda_runtime_validate_submission() {
  if [[ "$ECODA_RUN_ID" == "$FIXTURE_COMPUTE_ID" ]]; then
    while IFS= read -r relative; do
      cmp "$HPC_SCRATCH_DIR/$relative" "$NAS_TARGET_DIR/$relative" || return 87
    done < "$FIXTURE_BASE/canonical_inputs.txt"
    printf 'DISPATCH_RUNTIME_BOUNDARY\\n' >> "$FIXTURE_BASE/events"
    return 86
  fi
}
ecoda_runtime_validate_bound_run() {
  cp "$FIXTURE_BASE/runtime.identity" "$ECODA_RUN_ROOT/manifests/runtime.identity"
}
''')
        # The NAS stand-in copies and verifies the actual finite inventories;
        # it interrupts once only AFTER sealing the composite receipt.
        nas_stub = base / 'nas_boundary.py'
        nas_stub.write_text('''import hashlib, json, os, pathlib, shutil, sys
P = pathlib.Path
def md5(p): return hashlib.md5(p.read_bytes()).hexdigest()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def checksum(p):
    P(str(p)+'.md5').write_text(f'MD5={md5(p)}\\nSIZE={p.stat().st_size}\\nPATH={p}\\n')
command = sys.argv[1]
args = dict(zip(sys.argv[2::2], sys.argv[3::2]))
base = P(os.environ['FIXTURE_BASE'])
if command == 'preflight':
    P(args['--report']).write_text('{}')
elif command in ('copy', 'publish-inputs'):
    source, destination = P(args['--source-root']), P(args['--destination-root'])
    inv = P(args['--inventory'])
    if command == 'publish-inputs':
        assert destination == P(args['--nas-root'])
        expected = set((base / 'canonical_inputs.txt').read_text().splitlines())
        assert {line.split('\\t')[0] for line in inv.read_text().splitlines()} == expected
        assert len(expected) == 18
        archive = json.loads(P(args['--archive-receipt']).read_text())
        import csv
        mapping = {row['original_path']: row for row in csv.DictReader(P(args['--archive-mapping']).open(), delimiter='\\t')}
        for line in inv.read_text().splitlines():
            relative, digest, size = line.split('\\t')
            dst = destination / relative
            if dst.exists() and (md5(dst) != digest or dst.stat().st_size != int(size)):
                row = mapping[str(dst)]
                archived = P(row['archive_path'])
                assert dst.stat().st_size == archived.stat().st_size == int(row['size'])
                assert md5(dst) == md5(archived) == row['md5']
                archived_rel = str(archived.relative_to(archive['destination_root']))
                assert f"{archived_rel}\\t{row['md5']}\\t{row['size']}" in P(archive['inventory_path']).read_text().splitlines()
    else:
        assert destination.parts[-2] == os.environ['ECODA_RUN_ID']
    for line in inv.read_text().splitlines():
        relative, digest, size = line.split('\\t')
        src, dst = source / relative, destination / relative
        assert md5(src) == digest and src.stat().st_size == int(size)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists() or md5(dst) != digest or dst.stat().st_size != int(size):
            shutil.copyfile(src, dst)
    receipt = P(args['--receipt'])
    data = dict(schema_version=1, state='VERIFIED', run_id=os.environ['ECODA_RUN_ID'], host='fixture',
        nas_root=args['--nas-root'], mount_source='fixture', source_root=str(source), destination_root=str(destination),
        selection_path=args['--selection-file'], selection_sha256=sha(P(args['--selection-file'])),
        inventory_path=str(inv), inventory_sha256=sha(inv), source_manifest_path=args['--source-manifest'],
        source_manifest_sha256=sha(P(args['--source-manifest'])), completed_at='fixture')
    receipt.write_text(json.dumps(data, sort_keys=True))
    checksum(receipt)
    with (base / 'events').open('a') as stream: stream.write('SEALED ' + receipt.name + '\\n')
    if receipt.name == 'precompute.receipt.json' and not (base / 'interrupted').exists():
        (base / 'interrupted').touch()
        sys.exit(71)
elif command == 'verify':
    receipt = P(args['--receipt'])
    data = json.loads(receipt.read_text())
    inventory = P(data['inventory_path'])
    assert sha(inventory) == data['inventory_sha256']
    for line in inventory.read_text().splitlines():
        relative, digest, size = line.split('\\t')
        target = P(data['destination_root']) / relative
        assert md5(target) == digest and target.stat().st_size == int(size)
else:
    raise AssertionError(command)
''')
        (tree / 'src/utils/bash/ecoda_nas_backup.sh').write_text('#!/bin/bash\nexec ' +
            shlex.quote(sys.executable) + ' ' + shlex.quote(str(nas_stub)) + ' "$@"\n')
        bin_dir = base / 'bin'
        bin_dir.mkdir()
        (bin_dir / 'sbatch').write_text('#!/bin/bash\necho UNEXPECTED_SBATCH >> "$FIXTURE_BASE/events"\nexit 93\n')
        (bin_dir / 'sbatch').chmod(0o755)
        fixture_env['PATH'] = str(bin_dir) + os.pathsep + fixture_env['PATH']
        cli = ['bash', str(tree / 'src/5_run_benchmark_methods/finalize_corrected_outputs.sh'),
            '--selection-file', str(selection), '--method-matrix', str(matrix), '--source-run', source_run, '--nas-root', str(nas)]
        def invoke(arguments):
            return subprocess.run([*cli, *arguments], env=fixture_env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        def replace_fixture_bytes(path, payload):
            temporary = path.with_name(path.name + '.fixture-replacement')
            temporary.write_bytes(payload)
            temporary.chmod(path.stat().st_mode & 0o777)
            os.replace(temporary, path)
        def live_snapshot(paths):
            return {path: (path.read_bytes(), path.stat().st_ino, path.stat().st_mtime_ns)
                for path in paths}
        inventory_report = runs / (compute + '_inventory')
        observed = invoke(['--inventory-only', '--report-root', str(inventory_report)])
        assert observed.returncode == 0, observed.stderr
        assert not (base / 'events').exists(), 'inventory crossed an external mutation boundary'
        expected = inventory_report / 'manifests/expected_pending.tsv'
        expected_text = '' if zero_residual else ''.join(
            f'{ds}\t{ns["VIEW"]}\t{method}\n'
            for ds in ns['DATASETS'] for method in ns['METHODS'][:4])
        assert expected.read_text() == expected_text
        assert all(path.read_bytes() == content for path, content in collisions.items())
        if terminal_status_only:
            rows = list(__import__('csv').DictReader(
                (inventory_report / 'manifests/artifact_decisions.tsv').open(), delimiter='\t'))
            python_rows = [row for row in rows if row['method'] in ns['PYTHON_METHODS']]
            assert len(python_rows) == 27 and all(row['action'] in ('REUSE', 'IMPORT') for row in python_rows)
            assert ns['fields'](old / 'metadata', metadata=True)['STATE'] == 'ACTIVE'
            terminal = old / 'status/terminal'
            protected_paths = [path for directory in (canonical, old, scratch / '_ecoda_owners')
                for path in directory.rglob('*') if path.is_file() and path != terminal]
            protected = live_snapshot(protected_paths)
            terminal_before = live_snapshot([terminal])
            saved_terminal = terminal.with_name('terminal.fixture-original')
            for case in ('missing', 'mismatched'):
                os.replace(terminal, saved_terminal)
                try:
                    if case == 'mismatched':
                        terminal.write_text('STATE=FAIL\nRUN_ID=wrong_producer\nREASON=fixture\nTIME=fixture\n')
                    fixture_env['ECODA_RUN_ID'] = compute + '_' + case
                    rejected_report = runs / (fixture_env['ECODA_RUN_ID'] + '_inventory')
                    rejected_terminal = invoke(['--inventory-only', '--report-root', str(rejected_report)])
                    assert rejected_terminal.returncode != 0 and (
                        'source run lacks matching terminal evidence' in rejected_terminal.stderr
                    ), rejected_terminal.stderr
                    assert live_snapshot(protected) == protected, 'terminal rejection changed artifacts or owners'
                    assert {path for directory in (canonical, old, scratch / '_ecoda_owners')
                        for path in directory.rglob('*') if path.is_file()
                        and path not in (terminal, saved_terminal)} == set(protected_paths)
                    assert not (base / 'events').exists(), 'terminal rejection crossed an external boundary'
                finally:
                    os.replace(saved_terminal, terminal)
                assert live_snapshot([terminal]) == terminal_before
            continue
        if foreign_promotion_only:
            # Skip already-proven interruption injections, not production guards.
            for marker in ('metadata_interrupted', 'owner_finalization_interrupted', 'interrupted'):
                (base / marker).touch()
            completed = invoke(['--execute', '--report-root', str(promotion), '--expected-pending', str(expected)])
            assert completed.returncode != 0, 'fixture must stop at external dispatcher runtime'
            assert 'DISPATCH_RUNTIME_BOUNDARY' in (base / 'events').read_text(), completed.stderr
            assert 'UNEXPECTED_SBATCH' not in (base / 'events').read_text()
            assert ns['fields'](promotion / 'metadata', metadata=True)['STATE'] == 'OK', completed.stderr
            assert (promotion / 'manifests/precompute_backup.json').is_file()
            canonical_paths = [path for path in canonical.rglob('*') if path.is_file()]
            owner_paths = list((scratch / '_ecoda_owners').rglob('owner'))
            prior_paths = [path for path in promotion.rglob('*') if path.is_file()]
            backup_paths = [path for path in nas.rglob('*') if path.is_file()]
            protected = live_snapshot([*canonical_paths, *owner_paths, *prior_paths, *backup_paths])
            fixture_env['ECODA_RUN_ID'] = compute + '_fresh'
            fresh_report = runs / (fixture_env['ECODA_RUN_ID'] + '_inventory')
            fresh = invoke(['--inventory-only', '--report-root', str(fresh_report)])
            assert fresh.returncode == 0, fresh.stderr
            fresh_rows = list(__import__('csv').DictReader(
                (fresh_report / 'manifests/artifact_decisions.tsv').open(), delimiter='\t'))
            python_rows = [row for row in fresh_rows if row['method'] in ns['PYTHON_METHODS']]
            assert len(python_rows) == 27 and all(row['action'] == 'REUSE' for row in python_rows)
            assert (fresh_report / 'manifests/expected_pending.tsv').read_text() == expected_text
            assert live_snapshot(protected) == protected, 'fresh inventory changed prior promotion or canonical evidence'
            assert set((scratch / '_ecoda_owners').rglob('owner')) == set(owner_paths)
            prior_lineage = promotion / 'manifests/promotion_lineage.tsv'
            prior_sidecar = Path(str(prior_lineage) + '.md5')
            restore = {path: path.read_bytes() for path in (prior_lineage, prior_sidecar)}
            try:
                replace_fixture_bytes(prior_lineage, restore[prior_lineage].replace(
                    source_run.encode(), promotion.name.encode()))
                fixture_checksum(prior_lineage)
                fixture_env['ECODA_RUN_ID'] = compute + '_tampered'
                tampered_report = runs / (fixture_env['ECODA_RUN_ID'] + '_inventory')
                rejected_lineage = invoke(['--inventory-only', '--report-root', str(tampered_report)])
                assert rejected_lineage.returncode != 0 and (
                    f'promotion receipt does not bind file: {prior_lineage}' in rejected_lineage.stderr
                ), rejected_lineage.stderr
                unchanged = {path: value for path, value in protected.items() if path not in restore}
                assert live_snapshot(unchanged) == unchanged, 'lineage rejection mutated protected state'
                assert set((scratch / '_ecoda_owners').rglob('owner')) == set(owner_paths)
            finally:
                for path, payload in restore.items():
                    replace_fixture_bytes(path, payload)
            assert all(path.read_bytes() == value[0] for path, value in protected.items())
            continue
        wrong_expected = base / 'wrong_pending.tsv'
        wrong_expected.write_text(f'Joanito\t{ns["VIEW"]}\tmrvi\n')
        fixture_checksum(wrong_expected)
        rejected = invoke(['--execute', '--report-root', str(promotion), '--expected-pending', str(wrong_expected)])
        assert rejected.returncode != 0 and 'expected pending changed' in rejected.stderr, rejected.stderr
        assert not (base / 'events').exists(), 'pending mismatch crossed an external mutation boundary'
        assert all(path.read_bytes() == content for path, content in collisions.items())
        execute = ['--execute', '--report-root', str(promotion), '--expected-pending', str(expected)]
        metadata_interrupted = invoke(execute)
        assert metadata_interrupted.returncode != 0 and (base / 'metadata_interrupted').exists(), metadata_interrupted.stderr
        copied_metadata = canonical / 'metadata/Joanito_sample_metadata.feather'
        assert copied_metadata.read_bytes() == (recovery / 'metadata' / copied_metadata.name).read_bytes()
        assert not Path(str(copied_metadata) + '.md5').exists()
        metadata_stat = (copied_metadata.stat().st_ino, copied_metadata.stat().st_mtime_ns)
        owner_interrupted = invoke(execute)
        assert owner_interrupted.returncode != 0 and (base / 'owner_finalization_interrupted').exists(), owner_interrupted.stderr
        interrupted_payload = next(iter(original_payloads))
        interrupted_owner = Path(fixture_shell('ecoda_artifact_owner_dir "$1"', interrupted_payload).strip()) / 'owner'
        assert ns['fields'](interrupted_owner)['STATE'] == 'FAIL', 'caught interruption did not terminate ownership'
        interrupted = invoke(execute)
        assert interrupted.returncode != 0 and (base / 'interrupted').exists(), interrupted.stderr
        assert (copied_metadata.stat().st_ino, copied_metadata.stat().st_mtime_ns) == metadata_stat
        fixture_shell('ecoda_validate_checksum "$1"', copied_metadata)
        manifests = promotion / 'manifests'
        sealed = {path: path.read_bytes() for name in ('archive', 'canonical', 'inputs', 'precompute')
            for path in (manifests / (name + '.receipt.json'), manifests / (name + '.inventory.tsv'))}
        input_receipt = json.loads((manifests / 'inputs.receipt.json').read_text())
        assert input_receipt['destination_root'] == str(nas)
        assert len((manifests / 'inputs.inventory.tsv').read_text().splitlines()) == 18
        assert 'input_comparisons.tsv' not in (manifests / 'inputs.inventory.tsv').read_text()
        assert not (nas / '_ecoda_backups/corrected_final' / promotion.name / 'inputs').exists()
        assert all((nas / path.relative_to(scratch)).read_bytes() == path.read_bytes() for path in input_sources)
        assert {path: (path.stat().st_ino, path.stat().st_mtime_ns) for path in matching_targets} == matching_stats
        assert (unrelated.read_bytes(), unrelated.stat().st_mode, unrelated.stat().st_ino, unrelated.stat().st_mtime_ns) == unrelated_before
        # A sealed receipt must never release compute when canonical input
        # bytes disappear/change, even though all isolated archives still match.
        canonical_input = nas / input_sources[-2].relative_to(scratch)
        saved_input = canonical_input.read_bytes()
        for case in ('missing', 'changed'):
            if case == 'missing':
                canonical_input.unlink()
            else:
                canonical_input.write_bytes(b'changed canonical input')
            events_before = (base / 'events').read_bytes()
            rejected_input = invoke(execute)
            assert rejected_input.returncode != 0, rejected_input.stderr
            assert (base / 'events').read_bytes() == events_before, 'invalid canonical input reached dispatcher'
            canonical_input.write_bytes(saved_input)
        decisions_before = (manifests / 'artifact_decisions.tsv').read_bytes()
        lineage_before = (manifests / 'promotion_lineage.tsv').read_bytes()
        for destination, content in original_payloads.items():
            assert destination.read_bytes() == content, destination
        assert len(original_payloads) == 27
        mapping = list(__import__('csv').DictReader((manifests / 'archive_mapping.tsv').open(), delimiter='\t'))
        assert str(sidecar_payload) not in {row['original_path'] for row in mapping}, 'matching H5AD redundantly archived for sidecar-only repair'
        for path, content in collisions.items():
            archived = next(row for row in mapping if row['original_path'] == str(path))
            assert Path(archived['archive_path']).read_bytes() == content
        archived_paths = {row['original_path'] for row in mapping}
        for destination in original_payloads:
            historical = recovery / destination.relative_to(canonical)
            for item in (historical, Path(str(historical) + '.runtime.json')):
                assert str(item) in archived_paths, f'historical payload missing from archive: {item}'
                assert str(Path(str(item) + '.md5')) in archived_paths
                owner_path = Path(fixture_shell('ecoda_artifact_owner_dir "$1"', item).strip()) / 'owner'
                record_path = fixture_shell('ecoda_artifact_record_path "$1" "$2"', item, source_run).strip()
                assert str(owner_path) in archived_paths
                assert record_path in archived_paths
        assert '/recovery_35row/' not in (manifests / 'canonical.inventory.tsv').read_text()
        def payload_stats():
            return {item: (item.stat().st_ino, item.stat().st_mtime_ns)
                for destination in original_payloads for item in
                (destination, Path(str(destination) + '.runtime.json'), Path(str(destination) + '.md5'),
                 Path(str(destination) + '.runtime.json.md5'))}
        completed_stats = payload_stats()
        first_destination = next(iter(original_payloads))
        # A self-consistent live recovery runtime must not redefine the
        # immutable IMPORT transition after archive and composite sealing.
        historical_runtime = Path(str(recovery / first_destination.relative_to(canonical)) + '.runtime.json')
        historical_sidecar = Path(str(historical_runtime) + '.md5')
        historical_record = Path(fixture_shell('ecoda_artifact_record_path "$1" "$2"',
            historical_runtime, source_run).strip())
        recovery_before = {path: path.read_bytes()
            for path in (historical_runtime, historical_sidecar, historical_record)}
        owner_paths = list((scratch / '_ecoda_owners').rglob('owner'))
        destination_paths = [path for path in canonical.rglob('*')
            if path.is_file() and not path.is_relative_to(recovery)]
        owners_before = live_snapshot(owner_paths)
        destinations_before = live_snapshot(destination_paths)
        sealed_before = live_snapshot([*sealed,
            *(Path(str(path) + '.md5') for path in sealed),
            manifests / 'archive_mapping.tsv', manifests / 'archive_mapping.tsv.md5'])
        try:
            changed_runtime = json.loads(recovery_before[historical_runtime])
            changed_runtime['time_secs'] = 2.
            replace_fixture_bytes(historical_runtime, json.dumps(changed_runtime).encode())
            changed_md5 = ns['digest'](historical_runtime)
            changed_size = historical_runtime.stat().st_size
            replace_fixture_bytes(historical_sidecar,
                f'MD5={changed_md5}\nSIZE={changed_size}\nPATH={historical_runtime}\n'.encode())
            record_values = {'MD5': changed_md5, 'SIZE': str(changed_size)}
            changed_record = ''.join(
                key + '=' + record_values.get(key, value) + '\n'
                for key, value in (line.split('=', 1)
                    for line in recovery_before[historical_record].decode().splitlines()))
            replace_fixture_bytes(historical_record, changed_record.encode())
            fixture_shell('ecoda_validate_checksum "$1"; ecoda_validate_artifact_record "$1" stage5_mrvi "$2"',
                historical_runtime, source_run)
            rejected_runtime = invoke(execute)
            assert rejected_runtime.returncode != 0 and (
                f'live predecessor changed after archival: {historical_runtime}' in rejected_runtime.stderr
            ), rejected_runtime.stderr
            assert live_snapshot(owner_paths) == owners_before, 'runtime rejection changed owners'
            assert set((scratch / '_ecoda_owners').rglob('owner')) == set(owner_paths)
            assert live_snapshot(destination_paths) == destinations_before, 'runtime rejection changed destinations'
            assert {path for path in canonical.rglob('*') if path.is_file()
                and not path.is_relative_to(recovery)} == set(destination_paths)
            assert live_snapshot(sealed_before) == sealed_before, 'runtime rejection changed sealed evidence'
            assert (manifests / 'artifact_decisions.tsv').read_bytes() == decisions_before
            assert (manifests / 'promotion_lineage.tsv').read_bytes() == lineage_before
        finally:
            for path, payload in recovery_before.items():
                replace_fixture_bytes(path, payload)
        if sealed_import_only:
            continue
        replace_fixture_bytes(first_destination, b'not an original or promoted state')
        tampered = invoke(execute)
        assert tampered.returncode != 0 and 'not an original planned state' in tampered.stderr, tampered.stderr
        assert all(path.read_bytes() == content for path, content in sealed.items())
        replace_fixture_bytes(first_destination, original_payloads[first_destination])
        # The deliberate atomic tamper/restore changes only this payload inode.
        completed_stats[first_destination] = (first_destination.stat().st_ino, first_destination.stat().st_mtime_ns)
        resumed = invoke(execute)
        assert resumed.returncode != 0, 'fixture must stop at external dispatcher runtime'
        assert 'DISPATCH_RUNTIME_BOUNDARY' in (base / 'events').read_text(), resumed.stderr
        assert 'UNEXPECTED_SBATCH' not in (base / 'events').read_text()
        assert payload_stats() == completed_stats, 'completed imports or companions were rewritten on resume'
        assert (copied_metadata.stat().st_ino, copied_metadata.stat().st_mtime_ns) == metadata_stat
        for name in ('archive', 'canonical', 'inputs', 'precompute'):
            assert (base / 'events').read_text().splitlines().count('SEALED ' + name + '.receipt.json') == 1
        assert all(path.read_bytes() == content for path, content in sealed.items()), 'resume changed a sealed backup'
        assert (manifests / 'artifact_decisions.tsv').read_bytes() == decisions_before
        assert (manifests / 'promotion_lineage.tsv').read_bytes() == lineage_before
        assert (manifests / 'pending_after_promotion.tsv').read_text() == expected_text
        assert all(path.read_bytes() == content for path, content in original_evidence.items())
print('PASS corrected promotion producer terminal-status CLI regression' if terminal_status_only else
      'PASS corrected promotion foreign-promotion fresh inventory CLI regression' if foreign_promotion_only else
      'PASS corrected promotion sealed IMPORT source CLI regression' if sealed_import_only else
      'PASS corrected promotion finite I/O, real scientific/schema/sample validation and ownership fixtures')
PY
