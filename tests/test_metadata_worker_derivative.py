"""Actual metadata worker + immutable runtime guards + obs-only Feather exports."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile

import h5py
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def proof(path):
    Path(f'{path}.md5').write_text(f'MD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\nPATH={path}\n')


def fields(path, values):
    path.write_text(''.join(f'{key}={value}\n' for key,value in values.items()))


def obs_only(path, tag='current'):
    path.parent.mkdir(parents=True, exist_ok=True)
    values = {'Sample': ['s1','s1','s2','s2','s3','s3'], 'specimen': ['p1']*2+['p2']*2+['p3']*2,
        'condition.l1': ['control']*2+['CKD']*4, 'tissue': ['kidney']*2+['cortex of kidney']*2+['renal medulla']*2,
        'suspension_type': ['cell']*2+['nucleus']*4, 'subclass.l1': ['Immune']*6,
        'subclass.l3': ['T cell','B cell']*3}
    if tag == 'legacy': values['Sample'] = ['old1','old1','old2','old2','old3','old3']
    with h5py.File(path,'w') as handle:
        obs=handle.create_group('obs'); obs.attrs['encoding-type']='dataframe'
        obs.create_dataset('_index',data=[f'c{i}' for i in range(6)],dtype=h5py.string_dtype())
        for key,value in values.items(): obs.create_dataset(key,data=value,dtype=h5py.string_dtype())
    proof(path)


def fixture(base, *, override=False):
    scratch=base/'scratch'; nas=base/'nas';nas.mkdir(parents=True)
    snapshot=base/'snapshot'/('1'*40); source=snapshot/'tree'; identity=snapshot/'identity'
    identity.mkdir(parents=True); (source/'aux').mkdir(parents=True)
    for relative in ('src/slurm_config.sh','src/utils/bash/ecoda_run_common.sh','src/utils/bash/ecoda_runtime.sh',
        'src/utils/bash/ecoda_derivative_artifacts.sh','src/utils/bash/h5ad_obs_audit_worker.sh',
        'src/utils/py/derivative_artifacts.py','src/utils/py/h5ad_source_identity.py','src/utils/py/datasets_io.py','src/utils/py/batch_contract.py',
        'src/utils/py/export_h5ad_sample_metadata.py','config_helper.R','pixi.toml','pixi.lock',
        'aux/scGateDB.rds','aux/genes.blocklist.rds','aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz'):
        target=source/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/relative,target)
    registry=json.loads((ROOT/'datasets.json').read_text()); config={'Kidney_KPMP_full':registry['Kidney_KPMP_full']}
    if override: config['Kidney_KPMP_full']['views']['batch_effect_uncorrected']['columns']={'label':'tissue'}
    (source/'datasets.json').write_text(json.dumps(config))
    archive=identity/'source.tar'
    with tarfile.open(archive,'w') as handle:handle.add(source,arcname='.')
    source_manifest=identity/'source.manifest'
    fields(source_manifest,{'FORMAT':1,'SOURCE_ROOT':source,'SOURCE_COMMIT':snapshot.name,
        'SOURCE_ARCHIVE_PATH':archive,'SOURCE_ARCHIVE_SHA256':digest(archive),
        'CONFIG_HELPER_SHA256':digest(source/'config_helper.R'),'DATASETS_SHA256':digest(source/'datasets.json'),
        'PIXI_TOML_SHA256':digest(source/'pixi.toml'),'PIXI_LOCK_SHA256':digest(source/'pixi.lock'),
        'AUX_ROOT':source/'aux','SCGATE_DB_BRANCH':'41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4'})
    (snapshot/'COMPLETE').touch()
    host=base/'host/.pixi/envs/py-cuda13'; (host/'bin').mkdir(parents=True);(host/'lib').mkdir()
    (host/'bin/python').write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n');(host/'bin/python').chmod(0o755)
    (host/'bin/Rscript').write_text('#!/bin/bash\nexit 0\n');(host/'bin/Rscript').chmod(0o755)
    runtime=scratch/'_ecoda_runtime/fixture';runtime.mkdir(parents=True)
    image=runtime/'ecoda.sif';image.write_bytes(b'immutable fixture runtime identity')
    runtime_manifest=runtime/'ecoda.sif.manifest'
    fields(runtime_manifest,{'FORMAT':2,'IMAGE_PATH':image,'IMAGE_SHA256':digest(image),'RUNTIME_ENV':'py-cuda13',
        'RUNTIME_LAYOUT':'relocated','CONTAINER_ENV_PREFIX':'/opt/ecoda/py-cuda13','BASE_IMAGE':'rockylinux:9',
        'PIXITAINER_VERSION':'0.8.3','PIXI_VERSION':'0.49.0','APPTAINER_VERSION':'1.3.2',
        'IMAGE_BUILD_GIT_REVISION':'fixture','IMAGE_PIXI_TOML_SHA256':digest(source/'pixi.toml'),
        'IMAGE_PIXI_LOCK_SHA256':digest(source/'pixi.lock')})
    for path in (image,runtime_manifest):path.chmod(0o444)
    runtime.chmod(0o555)
    (base/'logs').mkdir(); inputs={}; legacy={}
    for view in ('batch_effect_uncorrected','batch_effect_corrected'):
        name=config['Kidney_KPMP_full']['views'][view]['output_file_name']
        path=scratch/'_ecoda_runs/stage3_source/artifacts/datasets/Kidney_KPMP_full/output'/name
        obs_only(path);inputs[view]=path
        old=scratch/'Kidney_KPMP_full/output'/name;obs_only(old,'legacy');legacy[view]=old
    environment=dict(os.environ,ECODA_HPC_CLUSTER='bamboo',ECODA_SOURCE_ROOT=str(source),ECODA_SOURCE_MANIFEST=str(source_manifest),
        ECODA_SOURCE_SNAPSHOT_REQUIRED='1',ECODA_HOST_ENV_PREFIX=str(host),ECODA_RUNTIME_MODE='host',ECODA_RUNTIME_IN_CONTAINER='0',
        ECODA_RUNTIME_IMAGE=str(image),ECODA_RUNTIME_MANIFEST=str(runtime_manifest),HPC_SCRATCH_DIR=str(scratch),
        ECODA_SCRATCH_ROOT=str(scratch),NAS_TARGET_DIR=str(nas),ECODA_NAS_TARGET_DIR=str(nas),ECODA_LOGS_DIR=str(base/'logs'),
        PROJECT_ROOT=str(source),DATASETS_JSON_FILE=str(source/'datasets.json'),PYTHONDONTWRITEBYTECODE='1',
        ECODA_HOST_PYTHON_SHA256=digest(host/'bin/python'),ECODA_HOST_RSCRIPT_SHA256=digest(host/'bin/Rscript'),USER_EMAIL='test@example.invalid')
    # Use genuine global owners; metadata consumer must reject a live writer.
    helper='source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"\n'
    for path in inputs.values():
        result=subprocess.run(['bash','-c',helper+'ecoda_artifact_owner_acquire "$1" stage3 stage3_source 0 0 0 >/dev/null\necoda_artifact_owner_set_state "$1" OK fixture >/dev/null', 'fixture',str(path)],env=environment,capture_output=True,text=True)
        assert result.returncode==0,result.stderr
    contexts=[]
    for pass_name,view in (('uncorrected','batch_effect_uncorrected'),('corrected','batch_effect_corrected')):
        run=scratch/f'_ecoda_runs/export_{pass_name}'; manifests=run/'manifests';manifests.mkdir(parents=True)
        (run/'status/metadata_export').mkdir(parents=True)
        fields(run/'metadata',{'STAGE':'stage5','RUN_ID':run.name,'STATE':'ACTIVE'})
        shutil.copyfile(source_manifest,manifests/'source.manifest')
        fields(manifests/'runtime.identity',{'RUNTIME_IMAGE':image,'RUNTIME_MANIFEST':runtime_manifest,
            'RUNTIME_IMAGE_SHA256':digest(image),'RUNTIME_MANIFEST_SHA256':digest(runtime_manifest),
            'RUNTIME_IMAGE_SIZE':image.stat().st_size,'RUNTIME_MANIFEST_SIZE':runtime_manifest.stat().st_size,
            'IMAGE_PIXI_TOML_SHA256':digest(source/'pixi.toml'),'IMAGE_PIXI_LOCK_SHA256':digest(source/'pixi.lock')})
        selection=manifests/'selection.tsv';selection.write_text(f'Kidney_KPMP_full\t{view}\t{view}\n')
        path=inputs[view]; refs=base/f'{pass_name}.inputs.tsv'
        refs.write_text(f'Kidney_KPMP_full\t{view}\t{path}\t{hashlib.md5(path.read_bytes()).hexdigest()}\t{path.stat().st_size}\t{path}.md5\n')
        binding=bind(run_root=run,release_id='fixture-release',nas_root=nas,stage='stage5',config_path=source/'datasets.json',selection=selection,input_artifacts=refs)
        output=run/f'artifacts/batch_effect/{pass_name}_final/metadata/Kidney_KPMP_full_sample_metadata.feather'
        manifest=manifests/'metadata_export_pending.tsv';manifest.write_text(f'Kidney_KPMP_full\t{view}\t{path}\t{output}\n');proof(manifest)
        env=dict(environment,ECODA_RUN_ID=run.name,ECODA_RUN_ROOT=str(run),ECODA_ARTIFACT_BINDING=str(binding),
            H5AD_OBS_AUDIT_MODE='metadata',H5AD_METADATA_EXPORT_MANIFEST=str(manifest),
            H5AD_METADATA_EXPORT_STATUS_DIR=str(run/'status/metadata_export'),H5AD_METADATA_EXPORT_TASK_ID='1',
            ANALYSIS_VARIANT='final' if pass_name=='uncorrected' else 'corrected_final',ANALYSIS_PASS=pass_name,
            ANALYSIS_ROOT=str(output.parent.parent),ANALYSIS_NAS_ROOT=str(nas/f'releases/fixture-release/{run.name}/batch_effect/{pass_name}_final'))
        contexts.append((view,env,output,manifest))
    for path in source.rglob('*'):
        path.chmod(0o555 if path.is_dir() else 0o444)
    source.chmod(0o555)
    for path in identity.iterdir():path.chmod(0o444)
    protected={path:(path.read_bytes(),path.stat().st_mode) for path in [*inputs.values(),*legacy.values(),*(Path(str(x)+'.md5') for x in [*inputs.values(),*legacy.values()])]}
    return source,contexts,inputs,legacy,protected


def worker(source,env,*,valid=True,match=None):
    completed=subprocess.run(['bash',str(source/'src/utils/bash/h5ad_obs_audit_worker.sh')],env=env,text=True,capture_output=True)
    assert (completed.returncode==0)==valid,(completed.returncode,completed.stdout,completed.stderr)
    if match:assert match in completed.stderr,completed.stderr
    return completed


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-metadata-worker-') as temporary:
        base=Path(temporary).resolve()
        try:
            source,contexts,inputs,legacy,protected=fixture(base/'valid')
            for view,env,output,manifest in contexts:
                # Compute nodes have the bound NAS identity but no mount.
                nas=Path(env['NAS_TARGET_DIR']);hidden=nas.with_name('nas-unmounted')
                nas.rename(hidden)
                try:
                    worker(source,env)
                    assert not nas.exists() # no synthetic mount/fallback
                finally:hidden.rename(nas)
                frame=pd.read_feather(output)
                assert frame.Sample.tolist()==['s1','s2','s3']
                assert frame.columns[:3].tolist()==['Sample','condition.l1','tissue']
                assert frame.tissue.tolist()==['kidney','cortex of kidney','renal medulla']
                before=(output.read_bytes(),Path(str(output)+'.md5').read_bytes())
                worker(source,env) # actual checksum/schema --check reuse
                assert before==(output.read_bytes(),Path(str(output)+'.md5').read_bytes())
                worker(source,dict(env,H5AD_METADATA_EXPORT_ROOT=str(base/'foreign')),valid=False,match='variant/root identity mismatch')
                worker(source,dict(env,H5AD_METADATA_EXPORT_NAS_ROOT=str(base/'foreign')),valid=False,match='NAS root does not match')
                missing=dict(env);missing.pop('ECODA_ARTIFACT_BINDING')
                worker(source,missing,valid=False,match='variant/root identity mismatch')
                original=manifest.read_bytes()
                manifest.write_bytes(original.replace(str(inputs[view]).encode(),str(legacy[view]).encode()));proof(manifest)
                worker(source,env,valid=False,match='not the configured')
                manifest.write_bytes(original);proof(manifest)
                bad=dict(env,ECODA_ARTIFACT_BINDING=str(base/'foreign.json'))
                worker(source,bad,valid=False,match='binding is invalid')
                source_bytes=inputs[view].read_bytes()
                inputs[view].write_bytes(b'X'+source_bytes[1:]) # same size, before --check reuse
                worker(source,env,valid=False,match='checksum is invalid')
                inputs[view].write_bytes(source_bytes)
                helper='source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"\n'
                for state in ('ACTIVE','OK'):
                    result=subprocess.run(['bash','-c',helper+'ecoda_artifact_owner_set_state "$1" "$2" fixture >/dev/null',
                        'fixture',str(inputs[view]),state],env=env,capture_output=True,text=True)
                    assert result.returncode==0,result.stderr
                    if state=='ACTIVE':worker(source,env,valid=False,match='invalid or active owner')
            # The absent-binding canonical lane keeps its established behavior.
            for view,env,_,manifest in contexts:
                legacy_env=dict(env);legacy_env.pop('ECODA_ARTIFACT_BINDING')
                pass_name=env['ANALYSIS_PASS']
                root=Path(env['HPC_SCRATCH_DIR'])/f'batch_effect/{pass_name}_final'
                output=root/'metadata/Kidney_KPMP_full_sample_metadata.feather'
                legacy_env.update(ANALYSIS_ROOT=str(root),ANALYSIS_NAS_ROOT=str(Path(env['NAS_TARGET_DIR'])/f'batch_effect/{pass_name}_final'))
                manifest.write_text(f'Kidney_KPMP_full\t{view}\t{legacy[view]}\t{output}\n');proof(manifest)
                worker(source,legacy_env)
                assert pd.read_feather(output).Sample.tolist()==['old1','old2','old3']
            for path,(data,mode) in protected.items(): assert path.read_bytes()==data and path.stat().st_mode==mode
            # Dataset-level authority rejects view columns before a binding can
            # even be made; do not fabricate a weaker binding for the worker.
            try:fixture(base/'override',override=True)
            except ValueError as error:assert 'invalid dataset/view' in str(error)
            else:raise AssertionError('view-level columns obtained a derivative binding')
            print('PASS: actual immutable metadata worker exports both bound views/axes; foreign roots/inputs/bindings/column overrides reject; source bytes/modes preserved')
        finally:
            for path in base.rglob('*'):
                if not path.is_symlink():path.chmod(0o700 if path.is_dir() else 0o600)


if __name__=='__main__':run()
