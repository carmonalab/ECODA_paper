"""Fixed HR22 publication with real copy/checksum/owner helpers and NAS stubs."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/utils/py'))
from benchmark_request_scope import AUTHOR_HR


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-hr22-publication-') as temporary:
        base = Path(temporary).resolve(); scratch = base / 'scratch'; nas = base / 'nas'
        tree = base / 'snapshot/tree'; old_source = base / 'old-source'
        for p in (scratch / '_ecoda_inputs', scratch / '_ecoda_runs', nas, tree / 'src/utils/bash', tree / 'src/utils/py', base / 'bin'):
            p.mkdir(parents=True)

        def write(p, text='fixture\n'):
            p.parent.mkdir(parents=True, exist_ok=True); p.write_text(text); return p

        def sha(p):
            return hashlib.sha256(p.read_bytes()).hexdigest()

        def proof(p):
            return dict(path=str(p), sha256=sha(p))

        def sealed(p, producer):
            md5 = hashlib.md5(p.read_bytes()).hexdigest(); size = str(p.stat().st_size)
            side = write(Path(str(p) + '.md5'), f'MD5={md5}\nSIZE={size}\nPATH={p}\n')
            key = hashlib.sha256(str(p).encode()).hexdigest()[:32]
            rec = write(old / 'manifests/artifacts' / (key + '.record'),
                        f'PATH={p}\nMD5={md5}\nSIZE={size}\nRUN_ID={old.name}\nPRODUCER={producer}\nSTATE=PUBLISHED\n')
            return dict(path=str(p), md5=md5, size=size, producer_record=proof(rec), sidecar=proof(side))

        def owned(p):
            key = p.name + '_' + hashlib.sha256(str(p).encode()).hexdigest()[:32]
            if len(key.encode()) > 100:
                key = 'artifact_' + hashlib.sha256(str(p).encode()).hexdigest()[:32]
            owner = write(scratch / '_ecoda_owners/artifact' / key / 'owner',
                          f'RUN_ID={old.name}\nSTAGE=stage5\nPATH={p}\nSTATE=OK\nPID=123\n')
            return proof(owner)

        for name in ('publish_bl_science.sh', 'ecoda_run_common.sh'):
            shutil.copyfile(ROOT / 'src/utils/bash' / name, tree / 'src/utils/bash' / name)
        for name in ('bl_science_publication.py', 'benchmark_request_scope.py', 'derivative_artifacts.py'):
            shutil.copyfile(ROOT / 'src/utils/py' / name, tree / 'src/utils/py' / name)
        write(tree / 'source.manifest', f'SOURCE_ROOT={tree}\n')
        write(tree / 'src/slurm_config.sh', '''export PROJECT_ROOT="$TEST_TREE" HPC_SCRATCH_DIR="$TEST_SCRATCH" NAS_TARGET_DIR="$TEST_NAS"
export PYTHON_BIN="$TEST_PYTHON" PYTHONDONTWRITEBYTECODE=1
''')
        write(tree / 'src/utils/bash/ecoda_runtime.sh', '''ecoda_runtime_validate_submission() {
 [[ "$1" == apptainer && "$ECODA_RUNTIME_PROFILE" == stage5 ]]
 printf 'RUNTIME_IMAGE_SHA256=new-runtime\n' > "$ECODA_RUN_ROOT/manifests/runtime.identity"
}
ecoda_runtime_validate_bound_run() {
 grep -q '^RUNTIME_IMAGE_SHA256=new-runtime$' "$ECODA_RUN_ROOT/manifests/runtime.identity"
 cmp -s "$ECODA_SOURCE_MANIFEST" "$ECODA_RUN_ROOT/manifests/source.manifest"
}
''')
        write(tree / 'src/utils/bash/ecoda_nas_backup.sh', '''#!/bin/bash
set -euo pipefail
[[ "$1" == preflight && "$3" == "$TEST_NAS" ]]
[[ "${TEST_MISSING_MOUNT:-0}" != 1 ]] || exit 7
printf '{"mount_source":"fixture-NAS","nas_root":"%s"}\n' "$TEST_NAS" > "$5"
''')
        chmod = write(base / 'bin/chmod', '''#!/bin/bash
if [[ "${TEST_MODE_FAILURE:-0}" == 1 && "$1" == 666 && "$2" == "$TEST_NAS/"* ]]; then exit 0; fi
exec /bin/chmod "$@"
'''); chmod.chmod(0o755)
        for command in ('sbatch', 'sacct', 'squeue', 'Rscript'):
            p = write(base / 'bin' / command, f'#!/bin/bash\ntouch "{base}/forbidden"\nexit 99\n'); p.chmod(0o755)

        old = scratch / '_ecoda_runs/stage5_bassez_author_hr22_20261109_v1'
        old.mkdir(parents=True)
        names = ['Bassez_metadata.rds', *(f'Bassez_{key}.rds' for key in AUTHOR_HR)]
        artifacts = []
        for name in names:
            p = write(old / 'artifacts/benchmark/results' / name, 'opaque approved synthetic result\n')
            item = sealed(p, 'stage5_composition'); item['owner'] = owned(p); artifacts.append(item)
        timing = sealed(write(old / 'logs/execution_times_composition_Bassez_selected.feather'), 'stage5_execution_log')
        requested = write(old / 'manifests/requested_selection.tsv', ''.join(f'Bassez\tbenchmark_analysis\tcomposition\t{key}\n' for key in AUTHOR_HR))
        image = write(base / 'image.sif'); image_manifest = write(base / 'image.sif.manifest', f'IMAGE_SHA256={sha(image)}\n')
        config = write(old_source / 'datasets.json', '{}\n')
        aux = old_source / 'aux/genes.blocklist.rds'; aux.parent.mkdir(); shutil.copyfile(ROOT / 'aux/genes.blocklist.rds', aux)
        source_manifest = write(old / 'manifests/source.manifest', f'SOURCE_ROOT={old_source}\nSOURCE_COMMIT=129c1a46ee7d282773eb090448ed3d061ac8ecf5\nAUX_ROOT={old_source}/aux\nDATASETS_SHA256={sha(config)}\n')
        runtime = write(old / 'manifests/runtime.identity', f'RUNTIME_IMAGE={image}\nRUNTIME_MANIFEST={image_manifest}\nRUNTIME_MANIFEST_SHA256={sha(image_manifest)}\nRUNTIME_IMAGE_SHA256={sha(image)}\nRUNTIME_IMAGE_SIZE={image.stat().st_size}\n')
        source_identity = write(old / 'manifests/source_identity.json', '{}\n')
        ids = ['4663753', '4663754', '4663757', '4663759', '4663760']
        raw = ''.join(f'{i}|COMPLETED|0:0\n' for i in ids)
        accounting = dict(attempted=True, query_count=1, actual_exit=0, scheduler_ids=ids, stdout=raw,
                          rows=[dict(line=i, job_id=job, state='COMPLETED', exit_code='0:0') for i, job in enumerate(ids, 1)])
        inspection = write(base / 'inspection.json', json.dumps(dict(state='FAILED', audit=dict(accounting=accounting))))
        row_owner = write(scratch / '_ecoda_owners/stage5/old-row/owner', f'RUN_ID={old.name}\nSTAGE=stage5\nSTATE=OK\nKEY=derivative/{old.name}/ordinary/Bassez/benchmark_analysis/composition\n')
        inventory = write(base / 'inventory.json', json.dumps(dict(schema='ecoda_bassez_hr22_failed_run_inventory_v1', run_id=old.name,
            row_owner=proof(row_owner), text_evidence=[dict(**proof(p), exists=True, truncated=False) for p in (source_manifest, runtime, source_identity, requested)])))
        cli = write(old_source / 'src/5_run_benchmark_methods/validate_benchmark_rds_contract.R')
        keys = ['anosim_score', 'mod_knn3_score', 'cluster_score', 'sil_score', 'lisi_score']
        science = dict(schema='ecoda_bassez_hr22_science_receipt_v1', status='NOOP_VALIDATED', science_validated=True,
            publication_completed=False, original_gate_status='FAILED', original_run_root=str(old), source_commit='129c1a46ee7d282773eb090448ed3d061ac8ecf5',
            n_samples=26, n_cells=74948, n_cell_types_high_res=40, counts_pseudocount=0.5, canonical_score_keys=keys,
            scores=[dict(method=method, **dict.fromkeys(keys, 0.5)) for method in AUTHOR_HR],
            inventory=proof(inventory), inspection=proof(inspection), original_accounting=accounting,
            source_manifest=proof(source_manifest), runtime_identity=proof(runtime), source_identity=proof(source_identity),
            input_record=proof(write(old / 'manifests/input.record')), input_checksum_proof=proof(write(base / 'input.expected.md5')),
            validation_script=proof(write(base / 'assert.R')), validation_execution=dict(R_version='4.5.2'),
            ari_source=proof(write(old_source / 'src/utils/scoring_metrics.R')), frozen_aux=dict(path=str(aux), md5=hashlib.md5(aux.read_bytes()).hexdigest()),
            production_validator=dict(**proof(cli), actual_exit=0, log=proof(write(base / 'cli.log', 'benchmark RDS artifact-list contract OK\n'))),
            artifacts=artifacts, timing=dict(proof=timing, rows=[dict(method=m) for m in AUTHOR_HR]))
        receipt = write(scratch / '_ecoda_inputs/hr22_science.json', json.dumps(science))
        review = write(scratch / '_ecoda_inputs/review.json', json.dumps(dict(status='APPROVED_SCIENTIFIC_ARTIFACT_REUSE', receipt_sha256=sha(receipt),
            datasets=1, rds_artifacts=23, timing_feathers=1, original_gate_status='FAILED', publication_completed=False)))
        original_paths = [p for top in (old, old_source, scratch / '_ecoda_owners') for p in top.rglob('*') if p.is_file()]
        snapshot = {p: (sha(p), p.stat().st_mode) for p in original_paths}
        environment = dict(os.environ, TEST_TREE=str(tree), TEST_SCRATCH=str(scratch), TEST_NAS=str(nas),
            TEST_PYTHON=sys.executable, ECODA_SOURCE_ROOT=str(tree), ECODA_SOURCE_MANIFEST=str(tree / 'source.manifest'),
            ECODA_SOURCE_SNAPSHOT_REQUIRED='1', PYTHONDONTWRITEBYTECODE='1', PATH=str(base / 'bin') + ':' + os.environ['PATH'])

        def invoke(name, **overrides):
            control = write(scratch / '_ecoda_inputs' / (name + '.json'), json.dumps(dict(schema='ecoda_bassez_hr22_publication_control_v1',
                run_id=name, release_id='release', nas_root=str(nas), science_receipt=proof(receipt), review=proof(review))))
            env = dict(environment, ECODA_RUN_ID=name, ECODA_RUN_ROOT=str(scratch / '_ecoda_runs' / name), **overrides)
            return subprocess.run(['bash', str(tree / 'src/utils/bash/publish_bl_science.sh'), '--control', str(control),
                '--control-md5', hashlib.md5(control.read_bytes()).hexdigest()], env=env, capture_output=True, text=True)

        result = invoke('good')
        assert result.returncode == 0, result.stdout + result.stderr
        run_root = scratch / '_ecoda_runs/good'; report = run_root / 'artifacts/evidence/bassez_hr22_publication.json'
        published = json.loads(report.read_text())
        assert published['schema'] == 'ecoda_bassez_hr22_publication_receipt_v1' and len(published['mappings']) == 24
        assert published['original_gate_status'] == 'FAILED' and published['original_accounting'] == accounting and not published['new_model_fit']
        assert len((Path(published['destination_root']) / 'checksums.md5').read_text().splitlines()) == 48
        assert all('PRODUCER=stage5_bassez_hr22_release_copy\n' in Path(row['scratch_record']['record_path']).read_text() for row in published['mappings'])
        control = run_root / 'manifests/control.json'
        checked = subprocess.run([sys.executable, str(tree / 'src/utils/py/bl_science_publication.py'), 'validate', '--control', str(control),
            '--control-md5', hashlib.md5(control.read_bytes()).hexdigest(), '--run-root', str(run_root), '--nas-root', str(nas),
            '--inventory', str(run_root / 'manifests/publication_inventory.tsv'), '--report', str(report)], env=environment, capture_output=True, text=True)
        assert checked.returncode == 0, checked.stderr
        assert invoke('good').returncode != 0
        assert invoke('mount_failure', TEST_MISSING_MOUNT='1').returncode != 0
        assert not (nas / 'releases/release/mount_failure').exists()
        assert invoke('mode_failure', TEST_MODE_FAILURE='1').returncode != 0
        owner_path = Path(artifacts[0]['owner']['path']); saved = owner_path.read_bytes()
        owner_path.write_text(saved.decode().replace('STATE=OK', 'STATE=ACTIVE'))
        assert invoke('active_old_owner').returncode != 0; owner_path.write_bytes(saved)
        saved = review.read_bytes(); review.write_text('{}')
        assert invoke('bad_review').returncode != 0; review.write_bytes(saved)
        saved = receipt.read_bytes(); changed = json.loads(saved); changed['scores'] = changed['scores'][:-1]
        receipt.write_text(json.dumps(changed)); review_saved = review.read_bytes(); new_review = json.loads(review_saved)
        new_review['receipt_sha256'] = sha(receipt); review.write_text(json.dumps(new_review))
        assert invoke('missing_method').returncode != 0
        receipt.write_bytes(saved); review.write_bytes(review_saved)
        saved = inspection.read_bytes(); inspection.write_text(saved.decode().replace('COMPLETED', 'FAILED'))
        assert invoke('accounting_changed').returncode != 0; inspection.write_bytes(saved)
        for path, before in snapshot.items():
            assert (sha(path), path.stat().st_mode) == before, path
        assert not (base / 'forbidden').exists()
        print('PASS: fixed HR22 24payload/48checksum release, truthful new records, original FAILED/all5 evidence, strict rootreview/source/owner/method/accounting/mount/mode guards')


if __name__ == '__main__':
    run()
