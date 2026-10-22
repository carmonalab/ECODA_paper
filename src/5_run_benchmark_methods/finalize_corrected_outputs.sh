#!/bin/bash
# Finite corrected-final maintenance. Scientific validation stays in the existing
# Python/R validators; the public dispatcher remains the only scheduler owner.
set -euo pipefail
fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
mode="" selection="" matrix="" source_run="" nas="" report="" expected=""
seen=" "
while [[ $# -gt 0 ]]; do
  key="$1"
  case "$key" in
    --inventory-only|--execute)
      [[ -z "$mode" ]] || fail 'exactly one maintenance mode is required'
      mode="$key"; shift ;;
    --selection-file|--method-matrix|--source-run|--nas-root|--report-root|--expected-pending)
      case "$seen" in *" $key "*) fail "duplicate option: $key";; esac
      seen="$seen$key "; [[ $# -ge 2 && -n "$2" ]] || fail "missing value: $key"
      case "$key" in
        --selection-file) selection="$2";; --method-matrix) matrix="$2";;
        --source-run) source_run="$2";; --nas-root) nas="$2";;
        --report-root) report="$2";; --expected-pending) expected="$2";;
      esac
      shift 2 ;;
    *) fail "unknown option: $key" ;;
  esac
done
[[ -n "$mode" && -n "$selection" && -n "$matrix" && -n "$nas" && -n "$report" ]] || fail 'missing required argument'
[[ "$source_run" == stage5_corrected_all9_ygg_retry5_20260916T044420Z ]] || fail 'unrecognized source run'
[[ "$mode" != --execute || -n "$expected" ]] || fail '--execute requires --expected-pending'
# Check inherited settings BEFORE the configuration can overwrite them.
[[ -z "${NAS_TARGET_DIR:-}" || "$NAS_TARGET_DIR" == "$nas" ]] || fail 'conflicting inherited NAS_TARGET_DIR'
[[ -z "${ECODA_NAS_TARGET_DIR:-}" || "$ECODA_NAS_TARGET_DIR" == "$nas" ]] || fail 'conflicting inherited ECODA_NAS_TARGET_DIR'
for path in "$selection" "$matrix" "$nas" "$report" ${expected:+"$expected"}; do
  [[ "$path" = /* && "$path" != *$'\n'* && "$path" != *$'\t'* && "$path" != *$'\r'* ]] || fail 'paths must be absolute and delimiter-free'
done
export NAS_TARGET_DIR="$nas" ECODA_NAS_TARGET_DIR="$nas" ECODA_HPC_CLUSTER=bamboo
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../slurm_config.sh"
[[ "$NAS_TARGET_DIR" == "$nas" && "$ECODA_NAS_TARGET_DIR" == "$nas" && "$ECODA_HPC_CLUSTER" == bamboo ]] || fail 'configured Bamboo/NAS identity disagrees'
source "$SCRIPT_DIR/../utils/bash/ecoda_run_common.sh"
source "$SCRIPT_DIR/../utils/bash/ecoda_stage5_policy.sh"
source "$SCRIPT_DIR/../utils/bash/ecoda_runtime.sh"
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && "${ECODA_SOURCE_ROOT:-}" == "$PROJECT_ROOT" ]] || fail 'immutable source executor required'
[[ "$SCRIPT_DIR" == "$ECODA_SOURCE_ROOT/src/5_run_benchmark_methods" ]] || fail 'wrapper is outside immutable source'
[[ "${ECODA_RUNTIME_MODE:-}" == apptainer ]] || fail 'pinned immutable runtime required'
[[ -n "${ECODA_RUN_ID:-}" ]] || fail 'executor compute identity is required'
for stale in ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION ANALYSIS_ROOT_VERSION ECODA_STAGE5_CORRECTED_FINAL_ROOT_IDENTITY; do
  eval 'value=${'"$stale"':-}'
  [[ -z "$value" ]] || fail "stale root identity: $stale"
done
export PYTHONDONTWRITEBYTECODE=1 ECODA_GATE_STAGE=stage5
export ECODA_RUNS_ROOT ECODA_OWNERS_ROOT
export ANALYSIS_VARIANT=corrected_final ANALYSIS_PASS=corrected PASS_ARG=corrected
export ANALYSIS_ROOT="$HPC_SCRATCH_DIR/batch_effect/corrected_final"
export ANALYSIS_NAS_ROOT="$nas/batch_effect/corrected_final"
export FINALIZE_SCRIPT_DIR="$SCRIPT_DIR"
"$PYTHON_BIN" - "$mode" "$selection" "$matrix" "$source_run" "$nas" "$report" "$expected" <<'PY'
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

DATASETS = ('Joanito', 'Stephenson', 'Alzheimer', 'Breast_cancer', 'Covid19_PBMC',
            'Kidney_KPMP_full', 'Diabetes', 'Lupus_PBMC', 'Lung')
METHODS = ('prepare_pseudobulk', 'pseudobulk', 'gloscope', 'composition', 'mrvi', 'pilot', 'qot')
PYTHON_METHODS = ('mrvi', 'pilot', 'qot')
VIEW = 'batch_effect_corrected'
HEADER = ('dataset', 'view', 'method', 'action', 'source_path', 'destination_path', 'source_run', 'md5', 'size', 'reason')
ACQUIRED_OWNERS = []


def safe_path(value, exists=False):
    path = Path(value)
    if not path.is_absolute() or str(path) != value or any(c in value for c in '\n\r\t'):
        raise ValueError(f'noncanonical absolute path: {value!r}')
    for part in (path, *path.parents):
        if part.is_symlink():
            raise ValueError(f'symlink path: {part}')
    if path.resolve() != path or (exists and not path.exists()):
        raise ValueError(f'unsafe or missing path: {path}')
    return path


def digest(path, algorithm='md5'):
    h = hashlib.new(algorithm)
    with safe_path(str(path), True).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fields(path, *, metadata=False):
    result = {}
    for line in safe_path(str(path), True).read_text().splitlines():
        if metadata and not line:
            continue
        key, sep, value = line.partition('=')
        if not sep or not key or (not value and not metadata) or key in result:
            raise ValueError(f'malformed evidence: {path}')
        result[key] = value
    return result


def atomic(path, data):
    path = safe_path(str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + '.build.', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(data.encode() if isinstance(data, str) else data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def shell(code, *args, capture=True, env=None):
    # Only constant internal commands enter Bash; paths are positional arguments.
    prelude = '''set -euo pipefail
source "$FINALIZE_SCRIPT_DIR/../slurm_config.sh"
source "$FINALIZE_SCRIPT_DIR/../utils/bash/ecoda_run_common.sh"
source "$FINALIZE_SCRIPT_DIR/../utils/bash/ecoda_stage5_policy.sh"
source "$FINALIZE_SCRIPT_DIR/../utils/bash/ecoda_runtime.sh"
'''
    return subprocess.run(['bash', '-c', prelude + code, 'finalize', *map(str, args)],
                          check=True, text=True, stdout=subprocess.PIPE if capture else None,
                          env=env).stdout


def checksum(path):
    if path.stat().st_size:
        shell('ecoda_write_checksum "$1"', path)
    else:
        shell('source "$FINALIZE_SCRIPT_DIR/stage5_pending.sh"; stage5_write_pending_checksum "$1"', path)


def write_report(path, data):
    atomic(path, data)
    checksum(path)


def copy_bytes(source, destination):
    shell('ecoda_atomic_copy "$1" "$2"', source, destination)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def producer(method):
    return 'stage5_prepare_pseudobulk_hvg2000' if method == 'prepare_pseudobulk' else 'stage5_' + method


def rebase_runtime(source, destination, artifact):
    # The validator has already rejected duplicate/nonfinite keys. Only this
    # path field may differ; no scientific/history metadata is synthesized.
    data = json.loads(source.read_text())
    data['artifact_path'] = str(artifact)
    atomic(destination, json.dumps(data, indent=2, allow_nan=False) + '\n')


def main(argv):
    mode, selection_arg, matrix_arg, source_run, nas_arg, report_arg, expected_arg = argv
    selection = safe_path(selection_arg, True)
    matrix = safe_path(matrix_arg, True)
    nas = safe_path(nas_arg, True)
    report = safe_path(report_arg)
    root = safe_path(os.environ['ECODA_SOURCE_ROOT'], True)
    scratch = safe_path(os.environ['HPC_SCRATCH_DIR'], True)
    runs = safe_path(os.environ['ECODA_RUNS_ROOT'], True)
    canonical = scratch / 'batch_effect/corrected_final'
    recovery = canonical / 'recovery_35row'
    compute_id = os.environ['ECODA_RUN_ID']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', compute_id):
        raise ValueError('unsafe compute identity')
    run_id = compute_id + ('_inventory' if mode == '--inventory-only' else '_promotion')
    if report != runs / run_id or (runs / compute_id).exists():
        raise ValueError('report must be distinct configured run root; compute root must not exist')
    required_selection = ''.join(f'{ds}\t{VIEW}\t{VIEW}\n' for ds in DATASETS)
    required_matrix = ''.join(f'{ds}\t{VIEW}\t{method}\n' for ds in DATASETS for method in METHODS)
    if selection.read_text() != required_selection or matrix.read_text() != required_matrix:
        raise ValueError('selection/matrix must be the fixed ordered nine-dataset/full63 coverage')
    shell('ecoda_validate_checksum "$1"; ecoda_validate_checksum "$2"', selection, matrix)
    expected = safe_path(expected_arg, True) if expected_arg else None
    if expected:
        shell('source "$FINALIZE_SCRIPT_DIR/stage5_pending.sh"; stage5_validate_pending_checksum "$1"', expected)
    source_manifest = safe_path(os.environ['ECODA_SOURCE_MANIFEST'], True)
    # Runtime helper validates immutable source/config/aux plus pinned image.
    # Inventory may write only its own identity reports, never owners.
    if not report.exists():
        shell('ecoda_init_run stage5 "$1"', run_id)
    else:
        shell('ecoda_open_run "$1" stage5', run_id)
        prior = fields(report / 'metadata', metadata=True)
        if prior.get('CUTOVER_ID') != compute_id:
            raise ValueError('report root is not a matching interrupted maintenance run')
    os.environ.update(ECODA_RUN_ID=run_id, ECODA_RUN_ROOT=str(report))
    manifests = report / 'manifests'
    for src, dst in ((selection, manifests / 'selection.tsv'), (matrix, manifests / 'method_matrix.tsv'),
                     (source_manifest, manifests / 'source.manifest')):
        if dst.exists() and dst.read_bytes() != src.read_bytes():
            raise ValueError(f'changed immutable maintenance input: {dst}')
        if not dst.exists():
            copy_bytes(src, dst)
            checksum(dst)
    shell('ecoda_runtime_validate_submission apptainer')
    shell('ecoda_runtime_validate_bound_run')
    metadata = fields(report / 'metadata', metadata=True)
    metadata.update(CUTOVER_ID=compute_id, MAINTENANCE_MODE=mode[2:], NAS_ROOT=str(nas),
                    SOURCE_RUN=source_run, SOURCE_MANIFEST_SHA256=digest(source_manifest, 'sha256'))
    atomic(report / 'metadata', ''.join(f'{k}={v}\n' for k, v in metadata.items()))
    sys.path.insert(0, str(root / 'src/utils/py'))
    identity_module = module('cutover_source_identity', root / 'src/utils/py/h5ad_source_identity.py')
    validator = module('cutover_matrix_validator', root / 'src/5_run_benchmark_methods/matrix_artifact_validator.py')
    identity = manifests / 'source_identity.json'
    current = identity_module.build_source_identity(selection, scratch, root / 'datasets.json')
    write_report(identity, json.dumps(current, indent=2) + '\n')
    current_records = identity_module.load_source_identity(identity)
    old_root = safe_path(str(runs / source_run), True)
    old_metadata = fields(old_root / 'metadata', metadata=True)
    old_terminal_path = old_root / 'status/terminal'
    old_terminal = fields(old_terminal_path) if old_terminal_path.is_file() else {}
    if (old_metadata.get('RUN_ID') != source_run
            or old_terminal.get('RUN_ID') != source_run or old_terminal.get('STATE') not in ('FAIL', 'OK')):
        raise ValueError('source run lacks matching terminal evidence')
    historical_identity = old_root / 'manifests/source_identity.json'
    historical = identity_module.load_source_identity(historical_identity)
    for key, value in current_records.items():
        if historical.get(key) != value:
            raise ValueError(f'BLOCKED: current input differs from original producer: {key}')
    historical_source = fields(old_root / 'manifests/source.manifest')
    if historical_source.get('SOURCE_COMMIT') != '7c3eb18ce258fa02552f85b617ce2e22c40a8b41':
        raise ValueError('BLOCKED: unexpected original producer source')
    for filename, key in (('datasets.json', 'DATASETS_SHA256'), ('pixi.toml', 'PIXI_TOML_SHA256'), ('pixi.lock', 'PIXI_LOCK_SHA256')):
        if historical_source.get(key) != digest(root / filename, 'sha256'):
            raise ValueError('BLOCKED: original scientific configuration is not current')
    shell('ecoda_validate_checksum "$1"', historical_identity)
    old_runtime = fields(old_root / 'manifests/runtime.identity')
    if not old_runtime.get('RUNTIME_IMAGE_SHA256') or not old_runtime.get('RUNTIME_MANIFEST_SHA256'):
        raise ValueError('BLOCKED: missing original producer runtime provenance')
    # The verified return is evidence only, never an overlay of Bamboo owners.
    transfer_manifest = safe_path(os.environ['ECODA_CORRECTED_SOURCE_EVIDENCE_MANIFEST'], True)
    transfer_receipt = safe_path(os.environ['ECODA_CORRECTED_SOURCE_EVIDENCE_RECEIPT'], True)
    returned = json.loads(transfer_receipt.read_text(), object_pairs_hook=validator._reject_duplicate_json_keys)
    transfer_fields = {'completed_at', 'files_sha256', 'incoming_root', 'inventory_sha256',
                       'selected_bytes', 'selected_count', 'transfer_id', 'transfer_script_sha256', 'verifier_sha256',
                       'verification'}
    if (set(returned) != transfer_fields or returned['transfer_id'] != 'corrected_return_20261021T130330Z'
            or returned['verification'] != 'all selected regular files matched inventory size and MD5'):
        raise ValueError('BLOCKED: unrecognized verified return receipt')
    pins = {
        'inventory_sha256': '2ec2cb4becda72b07453816acf3e55d804b631bfc4661217c580723c50bc6ea1',
        'files_sha256': '1d9d82c6f6f1e5b741fd1b8be24afe7c5c34d9569664a9fee3d2349dbdbe2dbf',
        'transfer_script_sha256': '972594ac79a113c13235445d9c5b50041b7e34374339721928156fb81320d28e',
        'verifier_sha256': '429831ce5f739e5d420c7ed0188575945af15366a2ba26285a5f088a07aef4d2'}
    if any(returned[key] != value for key, value in pins.items()) or digest(transfer_manifest, 'sha256') != pins['inventory_sha256']:
        raise ValueError('BLOCKED: verified return identity changed')
    if returned['selected_count'] != 1253 or returned['selected_bytes'] != 216609447099 or not returned['completed_at']:
        raise ValueError('BLOCKED: incomplete verified return')
    incoming = safe_path(returned['incoming_root'], True)
    if incoming != scratch.parent / '_ecoda_imports' / returned['transfer_id']:
        raise ValueError('BLOCKED: unexpected incoming root')
    transfer_rows = {}
    for line in transfer_manifest.read_text().splitlines():
        relative, size, md5 = line.split('\t')
        if relative in transfer_rows or Path(relative).is_absolute() or '..' in Path(relative).parts:
            raise ValueError('BLOCKED: malformed finite return inventory')
        transfer_rows[relative] = (int(size), md5)
    if len(transfer_rows) != returned['selected_count'] or sum(size for size, _ in transfer_rows.values()) != returned['selected_bytes']:
        raise ValueError('BLOCKED: returned scope mismatch')
    for src, basename in ((transfer_manifest, 'source_evidence.inventory.tsv'), (transfer_receipt, 'source_evidence.receipt.json')):
        dst = manifests / basename
        if dst.exists() and dst.read_bytes() != src.read_bytes():
            raise ValueError('BLOCKED: changed bound source evidence')
        if not dst.exists():
            copy_bytes(src, dst)
            checksum(dst)
    metadata.update(SOURCE_EVIDENCE_MANIFEST_SHA256=digest(transfer_manifest, 'sha256'),
                    SOURCE_EVIDENCE_RECEIPT_SHA256=digest(transfer_receipt, 'sha256'),
                    SOURCE_TRANSFER_ID=returned['transfer_id'])
    atomic(report / 'metadata', ''.join(f'{k}={v}\n' for k, v in metadata.items()))

    evidence = set()
    archives = set()
    decisions = []
    logical = {}
    selected_files = set()
    evidence.update([manifests / 'source_evidence.inventory.tsv', manifests / 'source_evidence.receipt.json'])
    row_paths = {}
    decision_path = manifests / 'artifact_decisions.tsv'
    planned = None
    if decision_path.exists():
        shell('ecoda_validate_checksum "$1"', decision_path)
        with decision_path.open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            if tuple(reader.fieldnames or ()) != HEADER:
                raise ValueError('BLOCKED: malformed original promotion plan')
            planned = [[row[key] for key in HEADER] for row in reader]
        if not planned:
            raise ValueError('BLOCKED: empty original promotion plan')
    checked_runtimes = set()
    completed_promotions = {}

    def receipt_binds(receipt, path):
        inventory = safe_path(receipt['inventory_path'], True)
        relative = str(path.relative_to(safe_path(receipt['source_root'], True)))
        matches = [line.split('\t') for line in inventory.read_text().splitlines()
                   if line.split('\t', 1)[0] == relative]
        if (len(matches) != 1 or len(matches[0]) != 3
                or matches[0][1] != digest(path) or int(matches[0][2]) != path.stat().st_size):
            raise ValueError(f'BLOCKED: promotion receipt does not bind file: {path}')

    def promotion_sources(run, ds, method):
        run_root = safe_path(str(runs / run), True)
        meta = fields(run_root / 'metadata', metadata=True)
        if run != run_id and meta.get('MAINTENANCE_MODE') != 'execute':
            return None
        foreign = run != run_id
        if foreign and run not in completed_promotions:
            if (meta.get('RUN_ID') != run or meta.get('STATE') != 'OK'
                    or run != meta.get('CUTOVER_ID', '') + '_promotion'
                    or meta.get('NAS_ROOT') != str(nas)):
                raise ValueError('BLOCKED: prior promotion is not completed')
            terminal = fields(run_root / 'status/terminal')
            if terminal.get('STATE') != 'OK' or terminal.get('RUN_ID') != run:
                raise ValueError('BLOCKED: prior promotion lacks terminal completion')
            prior = run_root / 'manifests'
            utility = root / 'src/utils/bash/ecoda_nas_backup.sh'
            receipt_path = prior / 'precompute_backup.json'
            subprocess.run(['bash', str(utility), 'verify', '--receipt', str(receipt_path)], check=True)
            composite = json.loads(receipt_path.read_text(), object_pairs_hook=validator._reject_duplicate_json_keys)
            proof = {run_root / 'metadata', run_root / 'status/terminal',
                     receipt_path, Path(str(receipt_path) + '.md5'),
                     Path(composite['inventory_path'])}
            if (composite['run_id'] != run or composite['source_root'] != str(run_root)
                    or composite['nas_root'] != str(nas)
                    or composite['source_manifest_sha256'] != digest(prior / 'source.manifest', 'sha256')):
                raise ValueError('BLOCKED: prior promotion receipt identity mismatch')
            children = {}
            for name in ('archive', 'canonical', 'inputs'):
                child_path = prior / (name + '.receipt.json')
                receipt_binds(composite, child_path)
                child = json.loads(child_path.read_text(), object_pairs_hook=validator._reject_duplicate_json_keys)
                if any(child[key] != composite[key] for key in
                       ('run_id', 'nas_root', 'source_manifest_sha256')):
                    raise ValueError('BLOCKED: prior promotion child identity mismatch')
                proof.update([child_path, Path(str(child_path) + '.md5')])
                # Input receipt bytes are sealed by the composite; lineage does
                # not need another full read of historical H5AD backup payloads.
                if name != 'inputs':
                    receipt_binds(composite, Path(child['inventory_path']))
                    subprocess.run(['bash', str(utility), 'verify', '--receipt', str(child_path)], check=True)
                    proof.add(Path(child['inventory_path']))
                children[name] = child
            lineage_path = prior / 'promotion_lineage.tsv'
            for path in (lineage_path, prior / 'artifact_decisions.tsv'):
                receipt_binds(composite, path)
                proof.update([path, Path(str(path) + '.md5')])
            if lineage_path.read_bytes() != (prior / 'artifact_decisions.tsv').read_bytes():
                raise ValueError('BLOCKED: prior promotion lineage differs from decisions')
            with lineage_path.open() as stream:
                reader = csv.DictReader(stream, delimiter='\t')
                if tuple(reader.fieldnames or ()) != HEADER:
                    raise ValueError('BLOCKED: malformed prior promotion lineage')
                rows = list(reader)
            completed_promotions[run] = (rows, children['canonical'], proof)
        if foreign:
            rows, canonical_receipt, proof = completed_promotions[run]
            evidence.update(proof)
        else:
            lineage_path = manifests / 'promotion_lineage.tsv'
            shell('ecoda_validate_checksum "$1"', lineage_path)
            rows = list(csv.DictReader(lineage_path.open(), delimiter='\t'))
            evidence.update([lineage_path, Path(str(lineage_path) + '.md5')])
        matching = [row for row in rows if row['dataset'] == ds and row['method'] == method]
        if (not matching or any(row['action'] not in ('IMPORT', 'ADOPT') or row['view'] != VIEW
                                for row in matching)
                or [row['destination_path'] for row in matching] !=
                   [str(path) for path in row_paths[(ds, method)]]):
            raise ValueError('BLOCKED: incomplete promotion lineage')
        if foreign:
            for row in matching:
                destination = safe_path(row['destination_path'], True)
                items = [destination]
                if method in PYTHON_METHODS:
                    items.append(Path(str(destination) + '.runtime.json'))
                for item in items:
                    receipt_binds(canonical_receipt, item)
                    receipt_binds(canonical_receipt, Path(str(item) + '.md5'))
        sources = {row['source_run'] for row in matching}
        if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', value) for value in sources):
            raise ValueError('BLOCKED: malformed promotion source identity')
        return sources

    def execution_producer(run, ds, method, trail=()):
        if run in trail:
            raise ValueError('BLOCKED: cyclic promotion lineage')
        sources = promotion_sources(run, ds, method)
        if sources is None:
            return run
        resolved = {execution_producer(value, ds, method, (*trail, run)) for value in sources}
        if len(resolved) != 1:
            raise ValueError('BLOCKED: ambiguous promotion execution producer')
        return resolved.pop()

    def owner(path, recovery_candidate=False):
        owner_dir = Path(shell('ecoda_artifact_owner_dir "$1"', path).strip())
        imported = False
        if not owner_dir.exists():
            if not recovery_candidate:
                raise ValueError(f'BLOCKED: missing owner: {path}')
            relative = str((owner_dir / 'owner').relative_to(scratch))
            if relative not in transfer_rows:
                raise ValueError('BLOCKED: source owner absent from verified finite return')
            imported_owner = safe_path(str(incoming / relative), True)
            size, md5 = transfer_rows[relative]
            if imported_owner.stat().st_size != size or digest(imported_owner) != md5:
                raise ValueError('BLOCKED: returned source owner bytes changed')
            owner_dir = imported_owner.parent
            imported = True
        if not imported:
            shell('_ecoda_artifact_owner_validate_dir "$1" "$2"', owner_dir, path)
        data = fields(owner_dir / 'owner')
        if imported and (not {'RUN_ID', 'STAGE', 'PATH', 'STATE', 'PID'}.issubset(data)
                         or set(data) - {'RUN_ID', 'STAGE', 'PATH', 'STATE', 'PID', 'REASON', 'TIME'}
                         or data['PATH'] != str(path) or data['RUN_ID'] != source_run
                         or not re.fullmatch(r'[1-9][0-9]*', data['PID'])):
            raise ValueError('BLOCKED: malformed imported original owner')
        if data.get('STATE') not in ('OK', 'FAIL') or data.get('STAGE') != 'stage5':
            raise ValueError(f'BLOCKED: active/malformed artifact owner: {path}')
        return owner_dir, data

    def original_evidence(run, ds, method, require_success=True, trail=()):
        if run in trail:
            raise ValueError('BLOCKED: cyclic promotion lineage')
        sources = promotion_sources(run, ds, method)
        if sources is not None:
            paths = set()
            for original_run in sources:
                paths.update(original_evidence(original_run, ds, method, require_success, (*trail, run)))
            return paths
        run_root = safe_path(str(runs / run), True)
        meta = fields(run_root / 'metadata', metadata=True)
        terminal_path = run_root / 'status/terminal'
        terminal = fields(terminal_path) if terminal_path.is_file() else {}
        if (meta.get('RUN_ID') != run or terminal.get('RUN_ID') != run
                or terminal.get('STATE') not in ('OK', 'FAIL')):
            raise ValueError('BLOCKED: producer run lacks matching terminal evidence')
        paths = {run_root / 'metadata', terminal_path, run_root / 'manifests/source.manifest',
                 run_root / 'manifests/runtime.identity', run_root / 'manifests/source_identity.json',
                 run_root / 'manifests/source_identity.json.md5'}
        source = fields(run_root / 'manifests/source.manifest')
        archive_path = safe_path(source.get('SOURCE_ARCHIVE_PATH', ''), True)
        expected_archive = source.get('SOURCE_ARCHIVE_SHA256', '')
        marker = archive_path.parent.parent / 'COMPLETE'
        if not re.fullmatch(r'[0-9a-f]{64}', expected_archive) or digest(archive_path, 'sha256') != expected_archive:
            raise ValueError('BLOCKED: original producer source archive identity mismatch')
        if safe_path(str(marker), True).read_text().strip() != 'COMPLETE':
            raise ValueError('BLOCKED: original producer snapshot is incomplete')
        runtime = fields(run_root / 'manifests/runtime.identity')
        for key in ('RUNTIME_IMAGE_SHA256', 'RUNTIME_MANIFEST_SHA256'):
            if not re.fullmatch(r'[0-9a-f]{64}', runtime.get(key, '')):
                raise ValueError('BLOCKED: missing original runtime provenance')
        if run not in checked_runtimes:
            image_path = safe_path(runtime.get('RUNTIME_IMAGE', ''), True)
            runtime_manifest = safe_path(runtime.get('RUNTIME_MANIFEST', ''), True)
            bound_runtime = fields(manifests / 'runtime.identity')
            image_sha = (bound_runtime['RUNTIME_IMAGE_SHA256'] if str(image_path) == bound_runtime.get('RUNTIME_IMAGE')
                         else digest(image_path, 'sha256'))
            if image_sha != runtime['RUNTIME_IMAGE_SHA256'] or digest(runtime_manifest, 'sha256') != runtime['RUNTIME_MANIFEST_SHA256']:
                raise ValueError('BLOCKED: original producer runtime bytes differ from its identity')
            checked_runtimes.add(run)
            paths.add(runtime_manifest)
        paths.update([archive_path, marker])
        for filename, key in (('datasets.json', 'DATASETS_SHA256'), ('pixi.toml', 'PIXI_TOML_SHA256'), ('pixi.lock', 'PIXI_LOCK_SHA256')):
            if require_success and source.get(key) != digest(root / filename, 'sha256'):
                raise ValueError('BLOCKED: producer scientific configuration mismatch')
        source_records = identity_module.load_source_identity(run_root / 'manifests/source_identity.json')
        if require_success and source_records.get((ds, VIEW)) != current_records[(ds, VIEW)]:
            raise ValueError('BLOCKED: producer current source/sample identity mismatch')
        matrices = sorted((run_root / 'manifests').glob('matrix_*.tsv'))
        matches = []
        for manifest in matrices:
            lines = [line.split('\t') for line in manifest.read_text().splitlines()]
            if any(row[:3] == [ds, VIEW, method] for row in lines):
                sidecar = safe_path(str(manifest) + '.md5')
                if sidecar.exists():
                    shell('ecoda_validate_checksum "$1"', manifest)
                else:
                    # Legacy producer matrices predate sidecar sealing. Their
                    # exact bytes must instead match the pinned return inventory.
                    relative = str(manifest.relative_to(scratch))
                    if transfer_rows.get(relative) != (manifest.stat().st_size, digest(manifest)):
                        raise ValueError(f'BLOCKED: unbound historical producer matrix: {manifest}')
                label = manifest.stem[len('matrix_'):]
                prefix = VIEW + '_' + method
                if label == prefix:
                    label = VIEW + '__' + method
                elif label.startswith(prefix + '_') and label[len(prefix) + 1:] in ('cpu', 'default_gpu', 'any_gpu'):
                    label = VIEW + '__' + method + '__' + label[len(prefix) + 1:]
                else:
                    raise ValueError('BLOCKED: unknown producer watchdog matrix suffix')
                status = run_root / 'status/watchdogs' / (label + '.status')
                status_fields = {}
                for line in safe_path(str(status), True).read_text().splitlines():
                    key, sep, value = line.partition('=')
                    if not sep or (key in status_fields and key != 'SCHEDULER_ID'):
                        raise ValueError('BLOCKED: malformed watchdog evidence')
                    status_fields[key] = value
                if status_fields.get('LABEL') != label or (require_success and status_fields.get('STATE') != 'OK'):
                    raise ValueError('BLOCKED: producer watchdog is not exact terminal OK')
                matches.extend([manifest, status])
                if sidecar.exists():
                    matches.append(sidecar)
        if not matches:
            raise ValueError(f'BLOCKED: missing exact method watchdog: {ds}/{method}')
        paths.update(matches)
        for path in paths:
            safe_path(str(path), True)
        evidence.update(paths)
        return paths

    def collision_evidence(item, ds, method):
        owner_path = safe_path(shell('ecoda_artifact_owner_dir "$1"', item).strip()) / 'owner'
        if owner_path.parent.exists():
            owner_dir, data = owner(item)
            owners = {owner_dir / 'owner'}
        else:
            # Archive-only reconciliation of legacy runtime metadata. Never use
            # a payload's ownership to authorize runtime reuse or invent an owner.
            if method not in PYTHON_METHODS or not str(item).endswith('.runtime.json'):
                raise ValueError(f'BLOCKED: missing collision owner: {item}')
            payload = Path(str(item)[:-len('.runtime.json')])
            owner_dir, data = owner(payload)
            owners = {owner_dir / 'owner'}
        if data['RUN_ID'] == run_id:
            raise ValueError('BLOCKED: promoted owner without a sealed predecessor archive')
        record = Path(shell('ecoda_artifact_record_path "$1" "$2"', item, data['RUN_ID']).strip())
        paths = {item, Path(str(item) + '.md5'), record}
        if not owner_path.parent.exists():
            for path in paths:
                safe_path(str(path), True)
                if transfer_rows.get(str(path.relative_to(scratch))) != (path.stat().st_size, digest(path)):
                    raise ValueError(f'BLOCKED: unbound legacy runtime collision: {path}')
        shell('ecoda_validate_checksum "$1"; ecoda_validate_artifact_record "$1" "$2" "$3"',
              item, producer(method), data['RUN_ID'])
        original_evidence(data['RUN_ID'], ds, method, False)
        return paths | owners


    def scientific(paths, ds, method, run):
        if method in PYTHON_METHODS:
            batch = validator._build_corrected_batch_contract(
                validator.read_datasets_json(root / 'datasets.json', view=VIEW), ds, VIEW, method)
            validator.require_nonempty(paths, f'{ds}/{method}', current_records[(ds, VIEW)]['sample_ids'],
                producer=producer(method), producer_run_id=run, expected_batch_contract=batch,
                require_runtime_batch_contract=True, require_corrected_summary=False)
            for path in paths:
                runtime = json.loads(Path(str(path) + '.runtime.json').read_text())
                expected_method = {'mrvi': 'MrVI_hvg2000', 'pilot': 'PILOT_hvg2000_highres', 'qot': 'QOT_hvg2000_highres'}[method]
                if runtime['dataset'] != ds or runtime['method'] != expected_method:
                    raise ValueError('wrong runtime dataset/method')
        else:
            for path in paths:
                args = shlex.split(os.environ['PIXI_RSCRIPT']) + [str(root / 'src/5_run_benchmark_methods/validate_benchmark_rds_contract.R'),
                    '--artifact', str(path), '--method', method, '--dataset', ds, '--view', VIEW,
                    '--input-root', str(scratch), '--config', str(root / 'datasets.json'),
                    '--batch-pass', 'corrected', '--analysis-variant', 'corrected_final',
                    '--source-identity', str(identity), '--source-identity-verified']
                if path.name.endswith('_metadata.rds'):
                    args.append('--metadata')
                validation_env = os.environ.copy()
                if path.is_relative_to(recovery):
                    # Existing R validator retains read-only historical root
                    # discovery; never export that root to workers/dispatcher.
                    for key in ('ANALYSIS_ROOT', 'ANALYSIS_NAS_ROOT', 'ANALYSIS_ROOT_IDENTITY'):
                        validation_env.pop(key, None)
                subprocess.run(args, check=True, env=validation_env)

    def inspect(paths, ds, method, recovery_candidate=False):
        all_evidence = set()
        states = []
        producer_runs = set()
        for path in paths:
            safe_path(str(path), True)
            companion = [path, Path(str(path) + '.runtime.json')] if method in PYTHON_METHODS else [path]
            for item in companion:
                shell('ecoda_validate_checksum "$1"', item)
                owner_dir, data = owner(item, recovery_candidate)
                run = data['RUN_ID']
                if recovery_candidate and run != source_run:
                    raise ValueError('BLOCKED: recovery producer is not approved source run')
                record = Path(shell('ecoda_validate_artifact_record "$1" "$2" "$3"', item, producer(method), run).strip())
                all_evidence.update([item, Path(str(item) + '.md5'), owner_dir / 'owner', record])
                states.append(data['STATE'])
                producer_runs.add(run)
        run = owner(paths[0], recovery_candidate)[1]['RUN_ID']
        scientific(paths, ds, method, run)
        # A scientifically invalid R cache is a repair candidate even if the
        # historical wave failed. Collision archival still verifies provenance.
        for producer_run in producer_runs:
            original_evidence(producer_run, ds, method)
        return run, states, all_evidence

    for ds in DATASETS:
        for method in METHODS:
            output = shell('_ecoda_stage5_artifacts_for "$1" "$2" "$3"; printf "%s\\n" "${ECODA_BENCHMARK_ARTIFACTS[@]}"', ds, VIEW, method)
            paths = [safe_path(value) for value in output.splitlines()]
            row_paths[(ds, method)] = paths
            original_group = [row for row in planned or [] if row[0] == ds and row[2] == method]
            if original_group:
                if [row[5] for row in original_group] != [str(path) for path in paths]:
                    raise ValueError('BLOCKED: original promotion destinations changed')
                action = original_group[0][3]
                if any(row[3] != action for row in original_group):
                    raise ValueError('BLOCKED: inconsistent original logical transition')
                if action in ('IMPORT', 'ADOPT', 'RECOMPUTE') and not (manifests / 'archive.receipt.json').exists():
                    # A reviewed pending mismatch may stop after the original
                    # plan is written but before its archive exists. Recollect
                    # predecessor bytes, not a new classification, on resume.
                    for path in paths:
                        companions = [path, Path(str(path) + '.runtime.json')] if method in PYTHON_METHODS else [path]
                        for item in companions:
                            if not item.exists():
                                continue
                            archives.update(collision_evidence(item, ds, method))
                logical[(ds, method)] = action
                decisions.extend(original_group)
                continue
            canonical_error = None
            action = None
            origin_paths = paths
            run = ''
            row_evidence = set()
            try:
                run, states, row_evidence = inspect(paths, ds, method)
                action = 'REUSE' if all(state == 'OK' for state in states) else 'ADOPT'
            except (ValueError, OSError, subprocess.CalledProcessError) as exc:
                canonical_error = str(exc)
                # Existing owners are checked independently: a failed scientific
                # validator never grants permission to overwrite an active owner.
                for path in paths:
                    for item in (path, Path(str(path) + '.runtime.json')):
                        if item.exists():
                            archives.update(collision_evidence(item, ds, method))
                candidates = [recovery / path.relative_to(canonical) for path in paths]
                try:
                    run, states, row_evidence = inspect(candidates, ds, method, True)
                    action = 'IMPORT'
                    origin_paths = candidates
                except (ValueError, OSError, subprocess.CalledProcessError) as candidate_error:
                    if method in PYTHON_METHODS:
                        action = 'BLOCKED'
                    elif any(path.exists() for path in candidates):
                        try:
                            for candidate in candidates:
                                if not candidate.exists():
                                    continue
                                _, candidate_owner = owner(candidate, True)
                                candidate_run = candidate_owner['RUN_ID']
                                shell('ecoda_validate_checksum "$1"; ecoda_validate_artifact_record "$1" "$2" "$3"', candidate, producer(method), candidate_run)
                                original_evidence(candidate_run, ds, method)
                            action = 'RECOMPUTE'
                        except (ValueError, OSError, subprocess.CalledProcessError):
                            action = 'BLOCKED'
                    elif canonical_error and ('BLOCKED:' in canonical_error):
                        action = 'BLOCKED'
                    else:
                        action = 'RECOMPUTE'
                    canonical_error = f'{canonical_error}; recovery: {candidate_error}'
            if action == 'ADOPT':
                archives.update(row_evidence)
            logical[(ds, method)] = action
            evidence.update(row_evidence)
            for source, destination in zip(origin_paths, paths):
                decisions.append([ds, VIEW, method, action, str(source), str(destination), run,
                                  digest(source) if source.is_file() else '', source.stat().st_size if source.is_file() else 0,
                                  ('recovery validated; invalid canonical: ' + str(canonical_error)) if action == 'IMPORT'
                                  else 'terminal FAIL independently validated' if action == 'ADOPT'
                                  else 'current independently validated contract' if action == 'REUSE' else canonical_error])
            if action in ('REUSE', 'ADOPT', 'IMPORT'):
                for path in paths:
                    selected_files.update([path, Path(str(path) + '.md5')])
                    if method in PYTHON_METHODS:
                        selected_files.update([Path(str(path) + '.runtime.json'), Path(str(path) + '.runtime.json.md5')])
    # Read the current metadata-only source through the existing exporter. This
    # compares all columns/order/values, including the intended Breast keys,
    # without materializing expression matrices or inventing a second policy.
    metadata_module = module('cutover_metadata', root / 'src/utils/py/export_h5ad_sample_metadata.py')
    metadata_imports = []
    metadata_plan_path = manifests / 'metadata_transitions.tsv'
    metadata_plan = None
    if metadata_plan_path.exists():
        shell('ecoda_validate_checksum "$1"', metadata_plan_path)
        with metadata_plan_path.open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            if reader.fieldnames != ['dataset', 'source_path', 'destination_path', 'md5', 'size']:
                raise ValueError('BLOCKED: malformed original metadata transitions')
            metadata_plan = list(reader)
        if [row['dataset'] for row in metadata_plan] != list(DATASETS):
            raise ValueError('BLOCKED: original metadata scope changed')
    metadata_transitions = []
    for ds in DATASETS:
        entry, view_spec = metadata_module._config_entry(root / 'datasets.json', ds, VIEW)
        sample_col, required, candidates = metadata_module.requested_columns(
            root / 'datasets.json', ds, entry, view_spec, batch_view=True)
        label_col = metadata_module._dataset_columns(entry)['label']
        keys, policy = metadata_module._resolved_batch_configuration(entry, view_spec, label_col, dataset=ds)
        frame = metadata_module.read_obs_metadata(Path(current_records[(ds, VIEW)]['path']),
            sample_col, required, candidates, metadata_module.READ_CHUNK_SIZE,
            batch_keys=keys, batch_metadata_policy=policy, biological_column=label_col)
        frame = frame.loc[:, metadata_module._ordered_metadata_columns(
            root / 'datasets.json', ds, entry, view_spec, frame.columns)]
        destination = canonical / 'metadata' / (ds + '_sample_metadata.feather')
        chosen = None
        original_metadata = metadata_plan[DATASETS.index(ds)] if metadata_plan is not None else None
        if original_metadata is not None and original_metadata['destination_path'] != str(destination):
            raise ValueError('BLOCKED: original metadata destination changed')
        candidates_to_validate = (safe_path(original_metadata['source_path'], True),) if original_metadata else (
            destination, recovery / 'metadata' / destination.name)
        for candidate in candidates_to_validate:
            try:
                safe_path(str(candidate), True)
                metadata_module.validate_output(candidate, required)
                observed = metadata_module.pd.read_feather(candidate)
                metadata_module.pd.testing.assert_frame_equal(frame, observed, check_dtype=False, check_categorical=False)
                chosen = candidate
                break
            except (ValueError, OSError, AssertionError):
                continue
        if chosen is None:
            raise ValueError(f'BLOCKED: no current-contract sample metadata: {ds}')
        if original_metadata and (digest(chosen) != original_metadata['md5']
                                 or chosen.stat().st_size != int(original_metadata['size'])):
            raise ValueError('BLOCKED: original metadata source changed')
        metadata_transitions.append(f'{ds}\t{chosen}\t{destination}\t{digest(chosen)}\t{chosen.stat().st_size}\n')
        if chosen != destination:
            for old in (destination, Path(str(destination) + '.md5')):
                if old.exists():
                    archives.add(old)
            metadata_imports.append((chosen, destination))
        selected_files.update([destination, Path(str(destination) + '.md5')])
    if metadata_plan is None:
        write_report(metadata_plan_path, 'dataset\tsource_path\tdestination_path\tmd5\tsize\n' + ''.join(metadata_transitions))
    evidence.update([metadata_plan_path, Path(str(metadata_plan_path) + '.md5')])
    log_module = module('cutover_execution_logs', root / 'src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py')
    source_shards = set()
    for (ds, method), action in logical.items():
        if action not in ('REUSE', 'IMPORT', 'ADOPT'):
            continue
        row = next(row for row in decisions if row[0] == ds and row[2] == method)
        shard_run = execution_producer(row[6], ds, method)
        shards = sorted((runs / shard_run / 'logs').glob(f'execution_times_batch_effect_corrected_final_{method}_{ds}*.feather'))
        if not shards:
            raise ValueError(f'BLOCKED: missing original execution shard: {ds}/{method}')
        for path in shards:
            frame = log_module._read_recorded_execution_log(path, producer='stage5_execution_log', run_id=shard_run, require_record=True)
            if frame is None or set(frame['dataset'].astype(str)) != {ds}:
                raise ValueError('BLOCKED: execution shard dataset mismatch')
            source_shards.add(path)
            evidence.update([path, Path(str(path) + '.md5'), Path(shell('ecoda_artifact_record_path "$1" "$2"', path, shard_run).strip())])
    if planned is not None:
        if [(row[0], row[1], row[2], row[5]) for row in decisions] != [
                (row[0], row[1], row[2], row[5]) for row in planned]:
            raise ValueError('BLOCKED: original promotion destinations changed')
        # Current classification is observational only. The immutable original
        # transitions control scope, lineage and every resumed mutation.
        for row in planned:
            ds, _, method, action, source, destination, source_run_id, md5, size, _ = row
            if action not in ('REUSE', 'ADOPT', 'IMPORT', 'RECOMPUTE', 'BLOCKED'):
                raise ValueError('BLOCKED: unknown original promotion transition')
            if action == 'IMPORT':
                candidate = safe_path(source, True)
                if digest(candidate) != md5 or candidate.stat().st_size != int(size):
                    raise ValueError('BLOCKED: original import source changed')
                evidence.update(inspect([candidate], ds, method, True)[2])
            elif action == 'REUSE':
                evidence.update(inspect(row_paths[(ds, method)], ds, method)[2])
                candidate = safe_path(destination, True)
                if digest(candidate) != md5 or candidate.stat().st_size != int(size):
                    raise ValueError('BLOCKED: original reused payload changed')
            elif action == 'ADOPT':
                candidate = safe_path(destination, True)
                if digest(candidate) != md5 or candidate.stat().st_size != int(size):
                    raise ValueError('BLOCKED: original adoption payload changed')
                original_evidence(source_run_id, ds, method)
        decisions = planned
        logical = {(row[0], row[2]): row[3] for row in decisions}
        selected_files = set()
        for (ds, method), action in logical.items():
            if action in ('REUSE', 'ADOPT', 'IMPORT'):
                for path in row_paths[(ds, method)]:
                    selected_files.update([path, Path(str(path) + '.md5')])
                    if method in PYTHON_METHODS:
                        selected_files.update([Path(str(path) + '.runtime.json'), Path(str(path) + '.runtime.json.md5')])
        for ds in DATASETS:
            destination = canonical / 'metadata' / (ds + '_sample_metadata.feather')
            selected_files.update([destination, Path(str(destination) + '.md5')])
    text = '\t'.join(HEADER) + '\n' + ''.join('\t'.join(map(str, row)).replace('\r', ' ').replace('\n', ' ') + '\n' for row in decisions)
    if planned is None:
        write_report(decision_path, text)
    if any(action == 'BLOCKED' for action in logical.values()):
        raise ValueError('BLOCKED: see artifact_decisions.tsv; no compute or promotion authorized')

    def classify(destination):
        reusable = manifests / 'classified_reuse.tsv'
        write_report(reusable, ''.join(f'{ds}\t{VIEW}\t{method}\n' for (ds, method), action in logical.items() if action in ('REUSE', 'IMPORT', 'ADOPT')))
        shell('''source "$FINALIZE_SCRIPT_DIR/stage5_pending.sh"
benchmark_selected_artifacts_valid() { local row; row="$(printf '%s\\t%s\\t%s' "$1" "$2" "$3")"; while IFS= read -r value; do [[ "$value" != "$row" ]] || return 0; done < "$3_REUSE"; return 1; }
'''.replace('"$3_REUSE"', shlex.quote(str(reusable))) + '''
stage5_method_is_forced() { return 1; }
stage5_prepare_pseudobulk_valid() { benchmark_selected_artifacts_valid "$1" "$2" prepare_pseudobulk; }
benchmark_artifacts_for() { ARTIFACT_PATHS=(); }
stage5_classify_pending_rows "$1" "$2"
''', manifests / 'method_matrix.tsv', destination)
        checksum(destination)
    pending = manifests / 'expected_pending.tsv'
    classify(pending)
    if any(line.split('\t')[-1] in PYTHON_METHODS for line in pending.read_text().splitlines()):
        raise ValueError('BLOCKED: excluded Python row in pending classification')
    if expected and pending.read_bytes() != expected.read_bytes():
        raise ValueError('expected pending changed; fresh reviewed inventory required')
    # Inspect only the exact selected method-owner keys; do not scan run roots.
    for line in pending.read_text().splitlines():
        ds, view, method = line.split('\t')
        method_owner = Path(shell('ecoda_owner_dir stage5 "$1"', f'corrected/{ds}/{view}/{method}').strip()) / 'owner'
        if method_owner.exists():
            owner_fields = fields(method_owner)
            if owner_fields.get('STATE') not in ('OK', 'FAIL') or not owner_fields.get('RUN_ID'):
                raise ValueError('BLOCKED: pending method owner is active/malformed')
            archives.add(method_owner)
            producer_root = runs / owner_fields['RUN_ID']
            for path in (producer_root / 'metadata', producer_root / 'manifests/source.manifest',
                         producer_root / 'manifests/runtime.identity'):
                safe_path(str(path), True)
                archives.add(path)
    sync_owner = Path(shell('ecoda_owner_dir stage5 "sync/$1"', canonical).strip()) / 'owner'
    sync_record = 'path\tmd5\tsize\trun_id\tstate\tarchive_path\n'
    if sync_owner.exists():
        data = fields(sync_owner)
        if (data.get('STATE') not in ('OK', 'FAIL') or data.get('STAGE') != 'stage5'
                or data.get('KEY') != 'sync/' + str(canonical)
                or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', data.get('RUN_ID', ''))):
            raise ValueError('BLOCKED: global synchronization owner active/malformed')
        archives.add(sync_owner)
        prior_run = runs / data['RUN_ID']
        for path in (prior_run / 'metadata', prior_run / 'manifests/source.manifest', prior_run / 'manifests/runtime.identity'):
            safe_path(str(path), True)
            archives.add(path)
        archive_owner = nas / '_ecoda_backups/corrected_final' / run_id / 'before' / str(sync_owner).lstrip('/')
        sync_record += f'{sync_owner}\t{digest(sync_owner)}\t{sync_owner.stat().st_size}\t{data["RUN_ID"]}\t{data["STATE"]}\t{archive_owner}\n'
    write_report(manifests / 'global_sync_owner.tsv', sync_record)
    # Compare only configured canonical NAS inputs; archive conflicts before
    # publishing current H5ADs and their original PATH-bound sidecars.
    comparisons = []
    input_files = []
    for record in current_records.values():
        path = Path(record['path'])
        relative = path.relative_to(scratch)
        target = safe_path(str(nas / relative))
        payload_matches = target.is_file() and target.stat().st_size == record['size'] and digest(target) == record['md5']
        target_sidecar = safe_path(str(target) + '.md5')
        sidecar_matches = (target_sidecar.is_file() and
            target_sidecar.read_bytes() == Path(str(path) + '.md5').read_bytes())
        matches = payload_matches and sidecar_matches
        comparisons.append(f'{relative}\t{record["md5"]}\t{record["size"]}\t{"MATCH" if matches else "MISSING_OR_MISMATCH"}\n')
        input_files.extend([path, Path(str(path) + '.md5')])
        for old, identical in ((target, payload_matches), (target_sidecar, sidecar_matches)):
            if old.is_file() and not identical:
                archives.add(old)
    write_report(manifests / 'input_comparisons.tsv', 'path\tmd5\tsize\tcomparison\n' + ''.join(comparisons))
    if mode == '--inventory-only':
        atomic(report / 'status/maintenance', 'STATE=INVENTORIED\nNO_COMPUTE=1\n')
        return

    utility = root / 'src/utils/bash/ecoda_nas_backup.sh'
    backup_root = nas / '_ecoda_backups/corrected_final' / run_id
    def backup(name, source_root, files, destination, command='copy', extra=()):
        files = sorted(set(safe_path(str(path), True) for path in files))
        relative = [str(path.relative_to(source_root)) for path in files]
        file_list = manifests / (name + '.files.txt')
        inventory = manifests / (name + '.inventory.tsv')
        receipt = manifests / (name + '.receipt.json')
        inventory_text = ''.join(f'{rel}\t{digest(path)}\t{path.stat().st_size}\n' for rel, path in zip(relative, files))
        if receipt.exists():
            subprocess.run(['bash', str(utility), 'verify', '--receipt', str(receipt)], check=True)
            if inventory.read_text() != inventory_text:
                raise ValueError(f'changed sealed backup scope: {name}')
            if command == 'publish-inputs' and json.loads(receipt.read_text())['destination_root'] != str(nas):
                raise ValueError('isolated historical inputs cannot authorize new compute')
            return receipt
        write_report(file_list, ''.join(value + '\n' for value in relative))
        write_report(inventory, inventory_text)
        subprocess.run(['bash', str(utility), command, '--source-root', str(source_root), '--nas-root', str(nas),
            '--destination-root', str(destination), '--files-from', str(file_list), '--inventory', str(inventory),
            '--receipt', str(receipt), '--selection-file', str(manifests / 'selection.tsv'),
            '--source-manifest', str(source_manifest), *map(str, extra)], check=True)
        subprocess.run(['bash', str(utility), 'verify', '--receipt', str(receipt)], check=True)
        return receipt

    subprocess.run(['bash', str(utility), 'preflight', '--nas-root', str(nas), '--report', str(manifests / 'nas_preflight.json')], check=True)
    archive_receipt = manifests / 'archive.receipt.json'
    if archive_receipt.exists():
        subprocess.run(['bash', str(utility), 'verify', '--receipt', str(archive_receipt)], check=True)
        shell('ecoda_validate_checksum "$1"', manifests / 'archive_mapping.tsv')
        mapping_rows = list(csv.DictReader((manifests / 'archive_mapping.tsv').open(), delimiter='\t'))
        archives = {Path(row['original_path']) for row in mapping_rows}
        for row in mapping_rows:
            archived = safe_path(row['archive_path'], True)
            if digest(archived) != row['md5'] or archived.stat().st_size != int(row['size']):
                raise ValueError('changed sealed collision archive')
    else:
        archives.update(evidence)
        archives.update([manifests / 'artifact_decisions.tsv', manifests / 'global_sync_owner.tsv'])
        # Only selected NAS predecessors enter the immutable staged archive.
        publication_paths = set(selected_files)
        for paths in row_paths.values():
            for path in paths:
                publication_paths.update([path, Path(str(path) + '.md5'),
                    Path(str(path) + '.runtime.json'), Path(str(path) + '.runtime.json.md5')])
        for path in publication_paths:
            if path.is_relative_to(canonical):
                target = nas / path.relative_to(scratch)
                if target.is_file():
                    archives.add(target)
        for relative in ('checksums.md5', 'embeddings/execution_times.feather', 'embeddings/execution_times.feather.md5'):
            target = nas / 'batch_effect/corrected_final' / relative
            if target.is_file():
                archives.add(target)
        staging = report / 'archive_staging'
        staged = []
        mapping = 'original_path\tarchive_path\tmd5\tsize\n'
        for path in sorted(archives):
            target = staging / str(path).lstrip('/')
            copy_bytes(path, target)
            if digest(path) != digest(target):
                raise ValueError('archive staging bytes changed')
            staged.append(target)
            mapping += f'{path}\t{backup_root / "before" / target.relative_to(staging)}\t{digest(path)}\t{path.stat().st_size}\n'
        write_report(manifests / 'archive_mapping.tsv', mapping)
        backup('archive', staging, staged, backup_root / 'before')
    lineage = manifests / 'promotion_lineage.tsv'
    if not lineage.exists():
        copy_bytes(decision_path, lineage)
        checksum(lineage)
    elif lineage.read_bytes() != decision_path.read_bytes():
        raise ValueError('BLOCKED: original promotion lineage changed')
    archive_rows = {row['original_path']: row for row in csv.DictReader(
        (manifests / 'archive_mapping.tsv').open(), delimiter='\t')}
    def assert_archived_live(path):
        row = archive_rows.get(str(path))
        if row is None or path.stat().st_size != int(row['size']) or digest(path) != row['md5']:
            raise ValueError(f'BLOCKED: live predecessor changed after archival: {path}')
    # Pin original IMPORT bytes before completed-row reuse or any live mutation;
    # valid replacement runtime records cannot redefine the sealed transition.
    for row in decisions:
        if row[3] != 'IMPORT':
            continue
        source = Path(row[4])
        originals = [source]
        if row[2] in PYTHON_METHODS:
            originals.append(Path(str(source) + '.runtime.json'))
        for item in originals:
            assert_archived_live(item)
            assert_archived_live(Path(str(item) + '.md5'))
    for source, destination in metadata_imports:
        if destination.is_file() and digest(destination) == digest(source):
            try:
                shell('ecoda_validate_checksum "$1"', destination)
            except subprocess.CalledProcessError:
                # The immutable transition authorizes this exact copied payload;
                # an interrupted checksum write does not authorize recopying it.
                sidecar = Path(str(destination) + '.md5')
                if sidecar.exists():
                    assert_archived_live(sidecar)
                checksum(destination)
            continue
        for old in (destination, Path(str(destination) + '.md5')):
            if old.exists():
                assert_archived_live(old)
        copy_bytes(source, destination)
        checksum(destination)
    # The verified receipt is the first mutation boundary. Sidecars and original
    # records are archived as bytes; original producer evidence is never edited.
    for ds in DATASETS:
        for method in METHODS:
            group = [row for row in decisions if row[0] == ds and row[2] == method]
            action = group[0][3]
            if action == 'REUSE':
                continue
            if action in ('IMPORT', 'ADOPT'):
                try:
                    current_run, states, _ = inspect(row_paths[(ds, method)], ds, method)
                    completed = current_run == run_id and all(state == 'OK' for state in states)
                    if completed:
                        for row in group:
                            destination = Path(row[5])
                            if digest(destination) != row[7] or destination.stat().st_size != int(row[8]):
                                raise ValueError('BLOCKED: completed payload differs from original promotion')
                            if action == 'IMPORT' and method in PYTHON_METHODS:
                                expected_runtime = json.loads(Path(row[4] + '.runtime.json').read_text())
                                expected_runtime['artifact_path'] = str(destination)
                                expected_bytes = (json.dumps(expected_runtime, indent=2, allow_nan=False) + '\n').encode()
                                if Path(str(destination) + '.runtime.json').read_bytes() != expected_bytes:
                                    raise ValueError('BLOCKED: completed runtime differs from original promotion')
                        logical[(ds, method)] = 'REUSE'
                        continue
                except (ValueError, OSError, subprocess.CalledProcessError):
                    # Partial ownership/companions still follow the original
                    # guarded transition below; no row is skipped on failure.
                    pass
            acquired = []
            for row in group:
                source, destination = Path(row[4]), Path(row[5])
                companions = [destination, Path(str(destination) + '.runtime.json')] if method in PYTHON_METHODS else [destination]
                if action == 'RECOMPUTE':
                    for item in companions:
                        if item.exists():
                            for old in (item, Path(str(item) + '.md5'), owner(item)[0] / 'owner'):
                                assert_archived_live(old)
                            shell('owner="$(ecoda_artifact_owner_dir "$1")"; ecoda_owner_reclaim_terminal "$owner" 0 0', item)
                            item.unlink()
                            sidecar = Path(str(item) + '.md5')
                            if sidecar.exists():
                                sidecar.unlink()
                    continue
                for item in companions:
                    if action == 'ADOPT' and owner(item)[1]['STATE'] == 'OK':
                        continue
                    owner_path = safe_path(shell('ecoda_artifact_owner_dir "$1"', item).strip())
                    if item.exists() and not owner_path.exists():
                        if action != 'IMPORT':
                            raise ValueError(f'BLOCKED: unowned destination is not an approved import: {item}')
                        # Original collision proof is sealed before mutation;
                        # the payload may now legitimately belong to this run.
                        for old in (item, Path(str(item) + '.md5')):
                            assert_archived_live(old)
                    elif item.exists() and owner(item)[1]['RUN_ID'] != run_id:
                        for old in (item, Path(str(item) + '.md5'), owner(item)[0] / 'owner'):
                            assert_archived_live(old)
                    elif item.exists() and action == 'IMPORT':
                        if item == destination:
                            promoted_md5 = digest(source)
                        else:
                            data = json.loads(Path(str(source) + '.runtime.json').read_text())
                            data['artifact_path'] = str(destination)
                            promoted_md5 = hashlib.md5((json.dumps(data, indent=2, allow_nan=False) + '\n').encode()).hexdigest()
                        predecessor = archive_rows.get(str(item))
                        if digest(item) != promoted_md5 and (
                                predecessor is None or digest(item) != predecessor['md5']
                                or item.stat().st_size != int(predecessor['size'])):
                            raise ValueError(f'BLOCKED: resumed destination is not an original planned state: {item}')
                    shell('ecoda_artifact_owner_acquire "$1" stage5 "$2" 0 1 1', item, run_id)
                    acquired.append(item)
                    ACQUIRED_OWNERS.append(item)
                if action == 'IMPORT':
                    copy_bytes(source, destination)
                    if method in PYTHON_METHODS:
                        rebase_runtime(Path(str(source) + '.runtime.json'), Path(str(destination) + '.runtime.json'), destination)
                for item in companions:
                    if item in acquired:
                        checksum(item)
                        shell('ecoda_write_artifact_record "$1" "$2" "$3"', item, producer(method), run_id)
            if action != 'RECOMPUTE':
                first = row_paths[(ds, method)][0]
                validation_run = run_id if first in acquired else owner(first)[1]['RUN_ID']
                scientific(row_paths[(ds, method)], ds, method, validation_run)
                for item in acquired:
                    shell('ecoda_artifact_owner_release "$1" OK "verified corrected final promotion"', item)
                    ACQUIRED_OWNERS.remove(item)
                logical[(ds, method)] = 'REUSE'
    # Revalidate every destination as a complete logical row before backing up.
    for (ds, method), action in logical.items():
        if action in ('REUSE', 'IMPORT', 'ADOPT'):
            inspect(row_paths[(ds, method)], ds, method)
    # Original execution-time shards and metadata are explicitly selected from
    # the identified producer manifests, never by a whole-run recursive copy.
    copied_logs = []
    for path in sorted(source_shards):
        target = report / 'logs' / path.name
        if target.exists() and digest(target) != digest(path):
            raise ValueError('changed promotion execution shard')
        if not target.exists():
            copy_bytes(path, target)
            checksum(target)
        copied_logs.append((path, target))
        selected_files.update([target, Path(str(target) + '.md5')])
    write_report(manifests / 'reused_execution_times.tsv', ''.join(f'{src}\t{dst}\t{digest(dst)}\t{dst.stat().st_size}\n' for src, dst in copied_logs))
    # Metadata is validated by the existing consumer contract against current
    # H5ADs and keys, not merely copied because a filename exists.
    subprocess.run(shlex.split(os.environ['PIXI_RSCRIPT']) + [str(root / 'src/5_run_benchmark_methods/validate_corrected_final_consumer_contracts.R'),
        '--selection', str(manifests / 'selection.tsv'), '--method-matrix', str(manifests / 'method_matrix.tsv'),
        '--config', str(root / 'datasets.json'), '--input-root', str(scratch), '--analysis-root', str(canonical),
        '--output', str(manifests / 'consumer_validation.json')], check=True)
    # Historical recovery payloads and original owner bytes live only in the
    # before archive. Validation accumulates evidence as rows are promoted;
    # that changing set must never define a sealed canonical backup scope.
    for basename in ('artifact_decisions.tsv', 'reused_execution_times.tsv', 'metadata_transitions.tsv',
                     'source_evidence.inventory.tsv', 'source_evidence.receipt.json'):
        selected_files.update([manifests / basename, Path(str(manifests / basename) + '.md5')])
    backup('canonical', scratch, selected_files, backup_root / 'canonical')
    # Seal the complete canonical set, including already matching inputs.
    # Comparison and archive mapping remain separate composite-chain evidence.
    backup('inputs', scratch, input_files, nas, command='publish-inputs',
        extra=('--archive-receipt', archive_receipt, '--archive-mapping', manifests / 'archive_mapping.tsv'))
    classify(manifests / 'pending_after_promotion.tsv')
    if pending.read_bytes() != (manifests / 'pending_after_promotion.tsv').read_bytes():
        raise ValueError('pending scope changed during promotion')
    chain = [manifests / 'input_comparisons.tsv', manifests / 'archive_mapping.tsv', manifests / 'global_sync_owner.tsv',
             manifests / 'artifact_decisions.tsv', manifests / 'promotion_lineage.tsv', metadata_plan_path,
             manifests / 'reused_execution_times.tsv', pending]
    for name in ('inputs', 'canonical', 'archive'):
        chain.extend([manifests / (name + '.receipt.json'), manifests / (name + '.receipt.json.md5'), manifests / (name + '.inventory.tsv'), manifests / (name + '.inventory.tsv.md5')])
    composite = backup('precompute', report, chain, backup_root / 'chain')
    final_receipt = manifests / 'precompute_backup.json'
    copy_bytes(composite, final_receipt)
    checksum(final_receipt)
    atomic(report / 'status/maintenance', 'STATE=BACKUP_VERIFIED\nNO_COMPUTE=1\n')
    shell('ecoda_set_run_state OK "corrected promotion backups verified; compute outcome is independent"')
    completed_metadata = fields(report / 'metadata', metadata=True)
    completed_metadata['STATE'] = 'OK'
    atomic(report / 'metadata', ''.join(f'{key}={value}\n' for key, value in completed_metadata.items()))
    env = os.environ.copy()
    env['ECODA_RUN_ID'] = compute_id
    env.pop('ECODA_RUN_ROOT', None)
    # One public launch, never another source executor. Empty pending is handled
    # by its existing finite validated publication path without fake jobs.
    subprocess.run(['bash', str(root / 'src/5_run_benchmark_methods/1_submit_hpc_array.sh'),
        '--selection-file', str(manifests / 'selection.tsv'), '--pass', 'corrected', '--analysis-variant', 'corrected_final',
        '--methods', ','.join(METHODS), '--method-matrix', str(manifests / 'method_matrix.tsv'),
        '--precompute-backup-receipt', str(final_receipt), '--expected-pending', str(pending)], env=env, check=True)


if __name__ == '__main__':
    try:
        main(sys.argv[1:])
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        for path in ACQUIRED_OWNERS:
            try:
                shell('ecoda_artifact_owner_validate "$1" "$ECODA_RUN_ID"; [[ "$ECODA_ARTIFACT_OWNER_STATE" != ACTIVE ]] || ecoda_artifact_owner_release "$1" FAIL "interrupted corrected promotion"', path)
            except (ValueError, OSError, subprocess.CalledProcessError):
                pass
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
PY
