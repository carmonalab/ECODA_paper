#!/bin/bash
# Explicit, run-bound reuse of an already completed scientific validation.
# Only derived outputs and small proof files are hashed before submission;
# immutable input record + exact inode/timestamps bind the large input bytes.
stage5_validate_recovery_manifest() {
  local phase="${1:-submit}"
  [[ -n "${ECODA_STAGE5_RECOVERY_MANIFEST:-}" && -n "${ECODA_STAGE5_RECOVERY_MANIFEST_SHA256:-}" ]] || return 1
  DATASETS_JSON_FILE="${DATASETS_JSON_FILE}" NAS_TARGET_DIR="${NAS_TARGET_DIR}" \
    "${PYTHON_BIN:-python3}" - "${ECODA_STAGE5_RECOVERY_MANIFEST}" "${ECODA_STAGE5_RECOVERY_MANIFEST_SHA256}" "${phase}" <<'PY'
import csv, hashlib, json, os, re, sys
from pathlib import Path

def fail(message):
    raise ValueError('BLOCKED: validated recovery: ' + message)
def safe(value):
    p = Path(value)
    if not p.is_absolute() or p.is_symlink() or not p.is_file() or str(p.resolve()) != str(p):
        fail('unsafe proof path: ' + str(p))
    return p
def digest(p, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with safe(str(p)).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()
def proof(item):
    p = safe(item['path'])
    if not re.fullmatch('[0-9a-f]{64}', item['sha256']) or digest(p) != item['sha256']:
        fail('proof bytes changed: ' + str(p))
    return p
def fields(p):
    result = {}
    for line in safe(str(p)).read_text().splitlines():
        k, sep, v = line.partition('=')
        if not sep or k in result:
            fail('malformed record: ' + str(p))
        result[k] = v
    return result
def rows(p, allow_empty=False):
    result = [tuple(row) for row in csv.reader(safe(str(p)).open(), delimiter='\t')]
    if (not result and not allow_empty) or any(len(r) != 3 or not all(r) for r in result) or len(set(result)) != len(result):
        fail('invalid selection')
    return result
def unique(pairs):
    result = {}
    for k, v in pairs:
        if k in result:
            fail('duplicate JSON key: ' + k)
        result[k] = v
    return result
try:
    manifest, expected, phase = sys.argv[1:]
    if digest(manifest) != expected:
        fail('manifest digest changed')
    m = json.loads(safe(manifest).read_text(), object_pairs_hook=unique)
    run_root = os.environ['ECODA_RUN_ROOT']
    run_id = os.environ.get('ECODA_RUN_ID') or os.environ.get('RUN_ID')
    if m['schema'] != 1 or m['run_root'] != run_root or m['run_id'] != run_id:
        fail('run binding mismatch')
    if str(Path(manifest).parent) != run_root + '/manifests':
        fail('manifest is not run-owned')
    if digest(run_root + '/manifests/source.manifest') != m['source_manifest_sha256']:
        fail('new source binding mismatch')
    if m['promotion_root'] != os.environ['ECODA_STAGE5_RECOVERY_PROMOTION_ROOT']:
        fail('promotion root mismatch')
    control_root = m['promotion_root'] + '/manifests/'
    required = ['artifact_decisions.tsv', 'metadata_transitions.tsv', 'source_identity.json',
                'expected_pending.tsv', 'archive_mapping.tsv', 'global_sync_owner.tsv']
    if any(control_root + name not in m['controls'] for name in required):
        fail('missing promotion control')
    config = os.environ['DATASETS_JSON_FILE']
    if config not in m['controls']:
        fail('scientific configuration not pinned')
    for path, sha in m['controls'].items():
        proof({'path': path, 'sha256': sha})
    archive_path = proof(m['archive_receipt'])
    if m['controls'].get(str(archive_path)) != m['archive_receipt']['sha256']:
        fail('archive receipt is not pinned')
    actual_runtime = fields(run_root + '/manifests/runtime.identity')
    old_runtime = fields(proof(m['preflights']['runtime_identity']))
    for key in ('RUNTIME_IMAGE_SHA256', 'RUNTIME_MANIFEST_SHA256'):
        if not re.fullmatch('[0-9a-f]{64}', m['runtime'][key]) or actual_runtime.get(key) != m['runtime'][key] or old_runtime.get(key) != m['runtime'][key]:
            fail('runtime mismatch')
    prior_source = fields(proof(m['preflights']['source_manifest']))
    if prior_source.get('DATASETS_SHA256') != digest(config):
        fail('scientific configuration differs from validated input preflight')
    r_log = proof(m['preflights']['r_log'])
    accounting = proof(m['preflights']['r_accounting'])
    scheduler_id = str(m['preflights']['r_scheduler_id'])
    if not scheduler_id.isdigit() or not r_log.name.endswith('_' + scheduler_id + '.log'):
        fail('R preflight scheduler identity mismatch')
    if 'benchmark R environment preflight OK' not in r_log.read_text().splitlines():
        fail('R environment success log missing')
    settled = [line.split('|') for line in accounting.read_text().splitlines()]
    if sum(len(r) >= 3 and r[:3] == [scheduler_id, 'COMPLETED', '0:0'] for r in settled) != 1:
        fail('R environment accounting is not completed successfully')
    h5ad_rows = rows(proof(m['preflights']['h5ad_manifest']))
    selection = rows(proof(m['selection']))
    pending = rows(proof(m['pending']), allow_empty=True)
    if m['pending']['sha256'] != m['controls'][control_root + 'expected_pending.tsv']:
        fail('pending differs from validated promotion')
    if len(selection) != 63 or not set(pending) < set(selection):
        fail('recovery requires full 63 coverage and an exact pending subset')
    if (not pending) != (m.get('release_only') is True):
        fail('empty pending requires explicit publication-only scope')
    if any(r[2] in ('mrvi', 'pilot', 'qot') for r in pending):
        fail('protected Python method pending')
    for actual, expected_rows in ((os.environ.get('ECODA_RECOVERY_DISPATCH'), selection),
                                  (os.environ.get('ECODA_RECOVERY_PENDING'), pending)):
        if actual and safe(actual).read_bytes() != safe(m['selection' if expected_rows is selection else 'pending']['path']).read_bytes():
            fail('actual dispatch/pending bytes differ')
    identity = json.loads(proof(m['source_identity']).read_text(), object_pairs_hook=unique)
    if str(m['source_identity']['path']) != control_root + 'source_identity.json':
        fail('source identity is not promotion identity')
    entries = {(e['dataset'], e['view']): e for e in identity['entries']}
    selected_inputs = {(r[0], r[1]) for r in selection}
    if set(entries) != selected_inputs or len(m['inputs']) != len(selected_inputs):
        fail('input coverage mismatch')
    seen_inputs = set()
    for item in m['inputs']:
        key = (item['dataset'], item['view'])
        if key in seen_inputs or key not in entries:
            fail('duplicate or foreign input')
        seen_inputs.add(key)
        if m.get('release_only') is True:
            # A result-only publication does not consume the H5AD inputs.
            e = entries[key]
            if item['path'] != e['path'] or item['md5'] != e['md5'] or int(item['size']) != int(e['size']):
                fail('recorded input provenance mismatch')
            continue
        p = safe(item['path']); st = p.stat(); e = entries[key]
        current = dict(dev=st.st_dev, ino=st.st_ino, size=st.st_size, mtime_ns=st.st_mtime_ns, ctime_ns=st.st_ctime_ns)
        if current != item['stat'] or e['path'] != str(p) or int(e['size']) != st.st_size or e['md5'] != item['md5'] or int(item['size']) != st.st_size:
            fail('input identity changed: ' + str(p))
        sidecar = fields(proof(item['sidecar']))
        record_path = proof(item['record']); record = fields(record_path)
        if st.st_mode & 0o222:
            fail('input payload is writable')
        for data in (sidecar, record):
            if data.get('PATH') != str(p) or data.get('MD5') != item['md5'] or data.get('SIZE') != str(st.st_size):
                fail('source record/sidecar mismatch')
        if record.get('STATE') != 'PUBLISHED' or record.get('PRODUCER') != 'stage5_preflight' or record.get('RUN_ID') != item['preflight_run_id']:
            fail('source record producer mismatch')
        status = fields(proof(item['preflight_status']))
        for k, v in dict(STATE='OK', RUN_ID=item['preflight_run_id'], DATASET=key[0], VIEW=key[1]).items():
            if status.get(k) != v:
                fail('source scientific preflight status mismatch')
        if (key[0], key[1], str(p)) not in h5ad_rows:
            fail('source absent from prior preflight')
        if phase == 'final' and digest(p, 'md5') != item['md5']:
            fail('source payload changed')
    if phase == 'submit':
        decisions = list(csv.DictReader(safe(control_root + 'artifact_decisions.tsv').open(), delimiter='\t'))
        reusable = set(selection) - set(pending)
        covered = set(); prep = set(); metadata = set(); paths = set()
        for path, sha in m['pending_owners'].items():
            proof({'path': path, 'sha256': sha})
        for item in m['artifacts']:
            p = proof(item)
            sidecar = fields(proof(item['sidecar']))
            if item['sidecar']['path'] != str(p) + '.md5' or sidecar != dict(PATH=str(p), MD5=digest(p, 'md5'), SIZE=str(p.stat().st_size)):
                fail('artifact strict sidecar mismatch')
            if str(p) in paths:
                fail('duplicate artifact')
            paths.add(str(p))
            key = (item['dataset'], item['view'], item['method'])
            kind = item['kind']
            if kind == 'metadata':
                metadata.add(item['dataset'])
                transitions = list(csv.DictReader(safe(control_root + 'metadata_transitions.tsv').open(), delimiter='\t'))
                if not any(r['dataset'] == item['dataset'] and r['destination_path'] == str(p) and r['md5'] == digest(p, 'md5') and int(r['size']) == p.stat().st_size for r in transitions):
                    fail('metadata not covered by validated transition')
                continue
            record = fields(proof(item['record'])); owner = fields(proof(item['owner']))
            if record.get('STATE') != 'PUBLISHED':
                fail('artifact record is not published')
            if owner.get('STATE') != 'OK' or owner.get('STAGE') != 'stage5' or owner.get('RUN_ID') != record.get('RUN_ID'):
                fail('artifact owner mismatch')
            if record.get('PATH') != str(p) or record.get('MD5') != digest(p, 'md5') or record.get('SIZE') != str(p.stat().st_size):
                fail('artifact record mismatch')
            if key[2] in ('mrvi', 'pilot', 'qot'):
                runtime = item['runtime']
                runtime_path = proof(runtime)
                runtime_sidecar = fields(proof(runtime['sidecar']))
                runtime_record = fields(proof(runtime['record']))
                if str(runtime_path) != str(p) + '.runtime.json':
                    fail('runtime proof belongs to another artifact')
                for data in (runtime_sidecar, runtime_record):
                    if data.get('PATH') != str(runtime_path) or data.get('MD5') != digest(runtime_path, 'md5') or data.get('SIZE') != str(runtime_path.stat().st_size):
                        fail('runtime record/sidecar mismatch')
                if runtime_record.get('STATE') != 'PUBLISHED' or runtime_record.get('RUN_ID') != record.get('RUN_ID') or runtime_record.get('PRODUCER') != record.get('PRODUCER'):
                    fail('runtime producer proof mismatch')
            if not any((r['dataset'], r['view'], r['method']) == key and r['action'] in ('REUSE', 'ADOPT') and r['destination_path'] == str(p) and r['md5'] == record['MD5'] and int(r['size']) == p.stat().st_size for r in decisions):
                fail('artifact missing validated decision')
            if kind == 'prep':
                if key[2] != 'prepare_pseudobulk': fail('invalid prep proof')
                prep.add(key[:2]); continue
            if kind != 'artifact' or key not in reusable:
                fail('artifact outside reuse scope')
            covered.add(key)
        if covered != {r for r in reusable if r[2] != 'prepare_pseudobulk'} or prep != {r[:2] for r in reusable if r[2] == 'prepare_pseudobulk'} or metadata != {r[0] for r in selection}:
            fail('incomplete reusable artifact/prep/metadata coverage')
    if len(m['nas_receipts']) < 2:
        fail('missing existing input/archive NAS receipts')
    input_receipts = []
    archive_receipts = []
    for receipt in m['nas_receipts']:
        receipt_path = proof(receipt)
        data = json.loads(receipt_path.read_text(), object_pairs_hook=unique)
        if data.get('state') != 'OK' or data.get('nas_root') != os.environ['NAS_TARGET_DIR']:
            fail('NAS receipt state/target mismatch')
        inventory = proof(dict(path=data['inventory_path'], sha256=data['inventory_sha256']))
        if receipt_path.name == 'inputs.receipt.json':
            input_receipts.append((data, inventory))
        if receipt_path == archive_path:
            archive_receipts.append(data)
    if len(input_receipts) != 1 or len(archive_receipts) != 1:
        fail('exact existing input/archive receipts required')
    data, inventory = input_receipts[0]
    inventory_rows = [line.split('\t') for line in inventory.read_text().splitlines()]
    for item in m['inputs']:
        relative = str(Path(item['path']).relative_to(data['source_root']))
        if sum(r == [relative, item['md5'], str(item['size'])] for r in inventory_rows) != 1:
            fail('NAS input backup does not cover current input')
except (OSError, ValueError, KeyError, TypeError) as exc:
    print(str(exc), file=sys.stderr)
    sys.exit(1)
PY
}

# Reuse the already-bound configuration identities in one process instead of
# restarting Python for every dataset/method pair before scientific submission.
stage5_recovery_install_batch_contracts() {
  local manifest="${ECODA_RUN_ROOT}/manifests/batch_contract.tsv"
  "${PYTHON_BIN}" - "${ECODA_STAGE5_RECOVERY_MANIFEST}" \
    "${ECODA_STAGE5_RECOVERY_MANIFEST_SHA256}" "${METHOD_MATRIX}" "${manifest}" <<'PY'
import csv, hashlib, json, os, re, sys
from pathlib import Path

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def publish(path, data):
    temp = Path(str(path) + '.build.' + str(os.getpid()))
    with temp.open('xb') as stream:
        stream.write(data)
    os.link(temp, path)
    temp.unlink()
try:
    source, expected, matrix, destination = sys.argv[1:]
    if sha(source) != expected:
        raise ValueError('recovery manifest changed before contract reuse')
    recovery = json.loads(Path(source).read_text())
    proof = recovery['preflights']['batch_contract_manifest']
    if sha(proof['path']) != proof['sha256']:
        raise ValueError('prior batch-contract manifest changed')
    selected = {tuple(row) for row in csv.reader(Path(matrix).open(), delimiter='\t')}
    required = selected | {(row[0], row[1], 'preprocess') for row in selected}
    output = Path(destination)
    directory = output.parent / 'batch_contracts'
    directory.mkdir(exist_ok=True)
    seen = set()
    rows = []
    for row in csv.reader(Path(proof['path']).open(), delimiter='\t'):
        if len(row) != 6:
            raise ValueError('malformed prior batch-contract row')
        dataset, view, method, path, md5, size = row
        key = (dataset, view, method)
        if key not in required or key in seen or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', value) for value in key):
            raise ValueError('prior batch-contract scope mismatch')
        original = Path(path)
        if original.is_symlink() or str(original.resolve()) != path:
            raise ValueError('unsafe prior batch-contract identity')
        data = original.read_bytes()
        if hashlib.md5(data).hexdigest() != md5 or len(data) != int(size):
            raise ValueError('prior batch-contract identity changed')
        target = directory / ('__'.join(key) + '.json')
        publish(target, data)
        publish(Path(str(target) + '.md5'), f'MD5={md5}\nSIZE={size}\nPATH={target}\n'.encode())
        rows.append('\t'.join((*key, str(target), md5, size)) + '\n')
        seen.add(key)
    if seen != required:
        raise ValueError('incomplete prior batch-contract coverage')
    publish(output, ''.join(rows).encode())
except (OSError, ValueError, KeyError) as error:
    sys.exit(str(error))
PY
  [[ $? -eq 0 ]] || return 1
  ecoda_write_checksum "${manifest}" || return 1
  ECODA_BATCH_CONTRACT_MANIFEST="${manifest}"
  ECODA_BATCH_CONTRACT_MANIFEST_MD5="${ECODA_CHECKSUM_MD5}"
  ECODA_BATCH_CONTRACT_MANIFEST_SIZE="${ECODA_CHECKSUM_SIZE}"
  ECODA_BATCH_CONTRACT_MANIFEST_SHA256="$(ecoda_sha256_file "${manifest}")"
}

stage5_verify_recovery_final_storage() {
  if jq -e '.release_only == true' "${ECODA_STAGE5_RECOVERY_MANIFEST}" >/dev/null; then
    # Verify only the already validated result bytes being transferred.
    # Do not rescan expression inputs or unrelated input-backup payloads.
    stage5_validate_recovery_manifest submit
    return $?
  fi
  stage5_validate_recovery_manifest final || return 1
  local receipt utility
  utility="${ECODA_SOURCE_ROOT:-${PROJECT_ROOT}}/src/utils/bash/ecoda_nas_backup.sh"
  while IFS= read -r receipt; do
    bash "${utility}" verify --receipt "${receipt}" || return 1
  done < <(jq -er '.nas_receipts[].path' "${ECODA_STAGE5_RECOVERY_MANIFEST}")
}

stage5_recovery_row_valid() {
  local ds="$1" view="$2" method="$3" path
  benchmark_artifacts_for "${ds}" "${view}" "${method}" || return 1
  for path in "${ARTIFACT_PATHS[@]}"; do
    jq -e --arg ds "${ds}" --arg view "${view}" --arg method "${method}" --arg path "${path}" \
      'any(.artifacts[]; .dataset==$ds and .view==$view and .method==$method and .path==$path)' \
      "${ECODA_STAGE5_RECOVERY_MANIFEST}" >/dev/null || return 1
  done
}
