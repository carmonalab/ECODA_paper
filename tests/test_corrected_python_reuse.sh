#!/bin/bash
# Exercise the dispatcher's real reuse predicate and real Feather validator.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
PYTHON="${PYTHON_BIN:-${ROOT}/.pixi/envs/default/bin/python}"
"${PYTHON}" - "${ROOT}" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

import h5py
import numpy as np
import pandas as pd

root = Path(sys.argv[1])
sys.path.insert(0, str(root / 'src/utils/py'))
from batch_contract import serialize_batch_contract_identity, validate_batch_metadata
from h5ad_source_identity import build_source_identity

source = (root / 'src/5_run_benchmark_methods/stage5_dispatcher.sh').read_text()
functions = []
for name in ('stage5_validate_path_component', 'stage5_recorded_producer',
             'stage5_allowed_producers_for', 'stage5_producer_allowed',
             'stage5_validate_reusable_artifact', 'stage5_batch_contract_identity_path',
             'benchmark_artifacts_for', 'benchmark_selected_artifacts_valid'):
    match = re.search(r'^' + name + r'\(\) \{\n.*?^\}', source, re.M | re.S)
    assert match, name
    functions.append(match.group())

with tempfile.TemporaryDirectory(prefix='ecoda-python-reuse-') as temporary:
    scratch = Path(temporary).resolve()
    runs = scratch / '_ecoda_runs'
    producer, consumer = 'original_producer', 'new_consumer'
    for run in (producer, consumer):
        (runs / run / 'manifests').mkdir(parents=True)
    analysis = scratch / 'batch_effect/corrected_final'
    (analysis / 'embeddings').mkdir(parents=True)
    nas = scratch / 'nas'
    nas.mkdir()
    env = dict(os.environ, HPC_SCRATCH_DIR=str(scratch), ECODA_RUNS_ROOT=str(runs),
        ECODA_OWNERS_ROOT=str(scratch / '_ecoda_owners'), ECODA_RUN_ROOT=str(runs / consumer),
        ECODA_RUN_ID=consumer, RUN_ID=consumer, PASS_ARG='corrected',
        ANALYSIS_VARIANT='corrected_final', ANALYSIS_PASS='corrected',
        ANALYSIS_ROOT=str(analysis), ANALYSIS_NAS_ROOT=str(nas / 'batch_effect/corrected_final'),
        NAS_TARGET_DIR=str(nas), ANALYSIS_LOG_PREFIX='execution_times_batch_effect_corrected_final_',
        SCRIPT_DIR=str(root / 'src/5_run_benchmark_methods'), DATASETS_JSON_FILE=str(root / 'datasets.json'),
        PYTHON_BIN=sys.executable, BENCHMARK_MATRIX_TEST='0')
    for key in ('ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION', 'ANALYSIS_ROOT_VERSION',
                'ECODA_STAGE5_CORRECTED_FINAL_ROOT_IDENTITY', 'ANALYSIS_ROOT_IDENTITY'):
        env.pop(key, None)
    common = root / 'src/utils/bash/ecoda_run_common.sh'
    prelude = 'set -euo pipefail\nsource ' + shlex.quote(str(common)) + '\n'
    def shell(code, *args):
        return subprocess.run(['bash', '-c', prelude + code, 'fixture', *map(str, args)],
                              env=env, text=True, capture_output=True)
    def checksum(path):
        checked = shell('ecoda_write_checksum "$1"', path)
        assert checked.returncode == 0, checked.stderr

    samples = [f's{i}' for i in range(1, 7)]
    obs = pd.DataFrame({'Sample': np.repeat(samples, 3), 'seqtec': ['A'] * 9 + ['B'] * 9})
    config = json.loads((root / 'datasets.json').read_text())
    h5ad = scratch / 'Joanito/output' / config['Joanito']['views']['batch_effect_corrected']['output_file_name']
    h5ad.parent.mkdir(parents=True)
    with h5py.File(h5ad, 'w') as handle:
        group = handle.create_group('obs')
        group.attrs['encoding-type'] = 'dataframe'
        group.attrs['encoding-version'] = '0.2.0'
        group.create_dataset('_index', data=np.asarray([f'c{i}' for i in range(18)], dtype='S3'))
        group.create_dataset('Sample', data=np.asarray(obs['Sample'], dtype='S2'))
    checksum(h5ad)
    selection = runs / consumer / 'manifests/selection.tsv'
    selection.write_text('Joanito\tbatch_effect_corrected\tbatch_effect_corrected\n')
    identity = runs / consumer / 'manifests/source_identity.json'
    identity.write_text(json.dumps(build_source_identity(selection, scratch, root / 'datasets.json')))
    checksum(identity)
    env['SOURCE_IDENTITY'] = str(identity)
    contract = serialize_batch_contract_identity(validate_batch_metadata(obs, ['seqtec']),
        method_id='MrVI', model_id='mrvi_composite_v1')
    batch = runs / consumer / 'manifests/mrvi_batch.json'
    batch.write_text(json.dumps(contract))
    batch_manifest = runs / consumer / 'manifests/batch_contract.tsv'
    batch_manifest.write_text(f'Joanito\tbatch_effect_corrected\tmrvi\t{batch}\t{hashlib.md5(batch.read_bytes()).hexdigest()}\t{batch.stat().st_size}\n')
    env['ECODA_BATCH_CONTRACT_MANIFEST'] = str(batch_manifest)
    payload = analysis / 'embeddings/Joanito_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather'
    matrix = np.abs(np.subtract.outer(np.arange(6), np.arange(6))).astype(float)
    pd.DataFrame(matrix, index=samples, columns=samples).to_feather(payload)
    runtime = Path(str(payload) + '.runtime.json')
    runtime.write_text(json.dumps(dict(schema_version=1, dataset='Joanito', method='MrVI_hvg2000',
        artifact_path=str(payload), artifact_md5=hashlib.md5(payload.read_bytes()).hexdigest(),
        time_secs=1.0, mem_GB=1.0, batch_contract=contract)))
    for item in (payload, runtime):
        checksum(item)
        created = shell('ecoda_artifact_owner_acquire "$1" stage5 "$2"; '
            'ecoda_write_artifact_record "$1" stage5_mrvi "$2"; '
            'ecoda_artifact_owner_release "$1" OK fixture', item, producer)
        assert created.returncode == 0, created.stderr
    protected = {p: (p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns)
        for p in scratch.rglob('*') if p.is_file() and p not in (identity, Path(str(identity) + '.md5'))}
    definition = '\n'.join(functions)
    policy = 'source ' + shlex.quote(str(root / 'src/utils/bash/ecoda_stage5_policy.sh')) + '\n'
    command = policy + definition + '\nbenchmark_selected_artifacts_valid Joanito batch_effect_corrected mrvi\n'
    reusable = shell(command)
    assert reusable.returncode == 0, reusable.stdout + reusable.stderr + 'valid cached MrVI result rejected'
    changed_identity = json.loads(identity.read_text())
    changed_identity['entries'][0]['sample_ids'].reverse()
    identity.write_text(json.dumps(changed_identity))
    checksum(identity)
    rejected = shell(command)
    assert rejected.returncode != 0, 'cached result with wrong source sample order accepted'
    assert all((p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns) == saved
        for p, saved in protected.items()), 'reuse validation changed artifacts or producer evidence'
    assert not list((runs / consumer / 'manifests').glob('artifact_check_*'))
print('PASS real corrected cached Python reuse, source-order rejection and unchanged artifacts')
PY
