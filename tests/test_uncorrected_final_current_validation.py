"""Validate the real current four-MrVI final scope without fitting any method."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind, context_root, load
from src.utils.py.h5ad_source_identity import build_source_identity

spec = importlib.util.spec_from_file_location('validator', ROOT / 'src/5_run_benchmark_methods/matrix_artifact_validator.py')
validator = importlib.util.module_from_spec(spec); spec.loader.exec_module(validator)
DATASETS = ('Covid19_PBMC', 'Kidney_KPMP_full', 'Lupus_PBMC', 'Lung')
VIEW = 'batch_effect_uncorrected'


def seal(path):
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    Path(str(path) + '.md5').write_text(f'MD5={md5}\nSIZE={path.stat().st_size}\nPATH={path}\n')
    return md5


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-unc-final-current-') as tmp:
        base = Path(tmp).resolve(); nas = base / 'nas'; nas.mkdir()
        config = base / 'datasets.json'; config.write_text((ROOT / 'datasets.json').read_text())
        registry = json.loads(config.read_text())
        for dataset in ('_diagnostic', 'Disabled', 'Viewless', 'NoOutput'):
            registry[dataset] = json.loads(json.dumps(registry['Lupus_PBMC']))
        registry['Disabled']['use_for_batch_effect'] = False
        registry['Viewless']['views'].pop(VIEW)
        registry['NoOutput']['views'][VIEW] = {}
        config.write_text(json.dumps(registry))
        run_root = base / 'scratch/_ecoda_runs/current-four'
        manifests = run_root / 'manifests'; manifests.mkdir(parents=True)
        for name in ('source.manifest', 'runtime.identity'):
            (manifests / name).write_text('immutable fixture identity\n')
        selection = manifests / 'selection.tsv'
        selection.write_text(''.join(f'{ds}\t{VIEW}\t{VIEW}\n' for ds in DATASETS)); seal(selection)
        references = base / 'references.tsv'; refs = []
        for ds in DATASETS:
            source = base / 'sources' / ds / registry[ds]['views'][VIEW]['output_file_name']; source.parent.mkdir(parents=True)
            with h5py.File(source, 'w') as handle:
                obs = handle.create_group('obs'); obs.attrs['encoding-type'] = 'dataframe'; obs.attrs['encoding-version'] = '0.2.0'
                obs.create_dataset('_index', data=np.asarray(['c1', 'c2', 'c3'], dtype='S2'))
                for column in {registry[ds]['columns']['sample'], 'Sample'}:
                    obs.create_dataset(column, data=np.asarray(['s1', 's1', 's2'], dtype='S2'))
            md5 = seal(source); refs.append(f'{ds}\t{VIEW}\t{source}\t{md5}\t{source.stat().st_size}\t{source}.md5\n')
        references.write_text(''.join(refs))
        binding = bind(run_root=run_root, release_id='test', nas_root=nas, stage='stage5', config_path=config, selection=selection, input_artifacts=references)
        value, _ = load(binding, config); root = context_root(value, 'batch_effect/uncorrected_final')
        env = dict(os.environ, PROJECT_ROOT=str(ROOT), PYTHON_BIN=sys.executable, HPC_SCRATCH_DIR=str(base / 'scratch'),
                   DATASETS_JSON_FILE=str(config), ECODA_RUN_ROOT=str(run_root), ECODA_RUN_ID=run_root.name,
                   ECODA_ARTIFACT_BINDING=str(binding), NAS_TARGET_DIR=str(nas), ANALYSIS_ROOT=str(root),
                   ANALYSIS_NAS_ROOT=str(context_root(value, 'batch_effect/uncorrected_final', nas=True)))
        with patch.dict(os.environ, env, clear=True):
            identity = manifests / 'source_identity.json'
            identity.write_text(json.dumps(build_source_identity(selection, base / 'scratch', config, validated_sidecars=True)))
            seal(identity)
            for ds in DATASETS:
                artifact = root / 'embeddings' / f'{ds}_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather'
                artifact.parent.mkdir(parents=True, exist_ok=True)
                pd.DataFrame([[0., 1.], [1., 0.]], index=['s1', 's2'], columns=['s1', 's2']).to_feather(artifact)
                md5 = seal(artifact)
                record = manifests / 'artifacts' / (hashlib.sha256(str(artifact).encode()).hexdigest()[:32] + '.record'); record.parent.mkdir(exist_ok=True)
                record.write_text(f'PATH={artifact}\nSIZE={artifact.stat().st_size}\nMD5={md5}\nRUN_ID={run_root.name}\nPRODUCER=stage5_mrvi\nSTATE=PUBLISHED\n')
            kwargs = dict(batch=True, batch_pass='uncorrected', analysis_variant='final', input_root=base / 'scratch',
                          config_path=config, source_identity=identity, source_identity_verified=True,
                          producer='stage5_mrvi', producer_run_id=run_root.name)
            validator.validate(root, selection, ['mrvi'], **kwargs)
            def rejected(path, message, **overrides):
                try:
                    validator.validate(root, path, ['mrvi'], **dict(kwargs, **overrides))
                except ValueError as error:
                    assert message in str(error), str(error)
                else:
                    raise AssertionError('invalid validation accepted: ' + message)
            probe = manifests / 'validation_probe.tsv'
            for dataset in ('Diabetes', '_debug', '_diagnostic', 'Disabled', 'Viewless', 'NoOutput', 'Unknown'):
                probe.write_text(f'{dataset}\t{VIEW}\t{VIEW}\n'); seal(probe)
                rejected(probe, 'unapproved dataset')
            probe.write_text('Lupus_PBMC\tbenchmark_analysis\tbenchmark_analysis\n'); seal(probe)
            rejected(probe, 'view mismatch')
            foreign_config = base / 'foreign.json'; foreign_config.write_bytes(config.read_bytes())
            rejected(selection, 'immutable --config', config_path=foreign_config)
            artifact = root / 'embeddings/Lupus_PBMC_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather'
            original = artifact.read_bytes(); artifact.write_bytes(original[:-1] + bytes([original[-1] ^ 1]))
            rejected(selection, 'missing/invalid')
            artifact.write_bytes(original)
            record = manifests / 'artifacts' / (hashlib.sha256(str(artifact).encode()).hexdigest()[:32] + '.record')
            original_record = record.read_text(); record.write_text(original_record.replace('PRODUCER=stage5_mrvi', 'PRODUCER=foreign'))
            rejected(selection, 'missing/invalid'); record.write_text(original_record)
            probe.write_text(f'Lupus_PBMC\t{VIEW}\t{VIEW}\n'); seal(probe)
            with patch.dict(os.environ, {}, clear=False):
                os.environ.pop('ECODA_ARTIFACT_BINDING')
                rejected(probe, 'unapproved dataset')
            with patch.dict(os.environ, {'ANALYSIS_ROOT': str(base / 'foreign')}):
                rejected(selection, 'immutable derivative binding')
            validator.validate(root, selection, ['mrvi'], **kwargs)
    print('PASS current bound four-MrVI uncorrected-final validation')


if __name__ == '__main__':
    run()
