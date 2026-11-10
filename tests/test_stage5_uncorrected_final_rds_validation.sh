#!/bin/bash
# Actual derivative binding and RDS CLI: no H5AD reads, model fits or Slurm.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
"${ROOT}/.pixi/envs/default/bin/python" - "${ROOT}" <<'PY'
import hashlib,importlib.util,json,os,pathlib,subprocess,sys,tempfile
root=pathlib.Path(sys.argv[1]);sys.path.insert(0,str(root))
from src.utils.py.derivative_artifacts import bind,load,context_root

def proof(p):pathlib.Path(str(p)+'.md5').write_text(f'MD5={hashlib.md5(p.read_bytes()).hexdigest()}\nSIZE={p.stat().st_size}\nPATH={p}\n')
with tempfile.TemporaryDirectory(prefix='ecoda-unc-final-rds-') as temp:
 b=pathlib.Path(temp).resolve();run=b/'scratch/_ecoda_runs/current';m=run/'manifests';m.mkdir(parents=True);nas=b/'nas';nas.mkdir()
 config=b/'datasets.json';registry=json.loads((root/'datasets.json').read_text());config.write_text(json.dumps(registry))
 for name in ['source.manifest','runtime.identity']:(m/name).write_text('immutable test identity\n')
 selection=m/'selection.tsv';selection.write_text('Joanito\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n');proof(selection)
 source=b/'sources/Joanito'/registry['Joanito']['views']['batch_effect_uncorrected']['output_file_name'];source.parent.mkdir(parents=True);source.write_bytes(b'NOT HDF5; must never open H5AD');proof(source)
 refs=b/'refs.tsv';refs.write_text(f'Joanito\tbatch_effect_uncorrected\t{source}\t{hashlib.md5(source.read_bytes()).hexdigest()}\t{source.stat().st_size}\t{source}.md5\n')
 binding=bind(run_root=run,release_id='fixture',nas_root=nas,stage='stage5',config_path=config,selection=selection,input_artifacts=refs)
 value,_=load(binding,config);analysis=context_root(value,'batch_effect/uncorrected_final');results=analysis/'results';results.mkdir(parents=True)
 ids=[f's{i}' for i in range(1,13)];identity=m/'source_identity.json';identity.write_text(json.dumps(dict(schema=1,entries=[dict(dataset='Joanito',view='batch_effect_uncorrected',path=str(source),size=source.stat().st_size,md5=hashlib.md5(source.read_bytes()).hexdigest(),sample_ids=ids)])));proof(identity)
 artifact=results/'Joanito_batch_effect_uncorrected_final_pseudobulk.rds';R=root/'.pixi/envs/default/bin/Rscript'
 code="""a<-commandArgs(TRUE);ids<-paste0('s',1:12);f<-matrix(as.double(1:36),12,dimnames=list(ids,c('a','b','c')));b<-list(feat_mat=f,dist_mat=dist(f),labels=setNames(factor(rep(c('A','B'),each=6)),ids),scores=as.list(setNames(rep(.5,5),c('anosim_score','mod_knn3_score','cluster_score','sil_score','lisi_score'))),exec_time=1,mem_GB=.1);saveRDS(list(Pseudobulk_hvg2000=b),a[1])"""
 subprocess.run([str(R),'-',str(artifact)],input=code,text=True,check=True,capture_output=True);proof(artifact)
 env={k:v for k,v in os.environ.items() if not k.startswith(('ECODA_','ANALYSIS_'))};env.update(PROJECT_ROOT=str(root),PYTHON_BIN=sys.executable,RETICULATE_PYTHON=sys.executable,DATASETS_JSON_FILE=str(config),ECODA_RUN_ID='current',ECODA_RUN_ROOT=str(run),ECODA_RUNS_ROOT=str(run.parent),ECODA_ARTIFACT_BINDING=str(binding),NAS_TARGET_DIR=str(nas),HPC_SCRATCH_DIR=str(b/'scratch'),ANALYSIS_ROOT=str(analysis),ANALYSIS_NAS_ROOT=str(context_root(value,'batch_effect/uncorrected_final',nas=True)))
 lines=(root/'src/5_run_benchmark_methods/stage5_dispatcher.sh').read_text().splitlines();start=lines.index('stage5_validate_final_rds_rows() {');end=next(i for i in range(start+1,len(lines)) if lines[i]=='}');function='\n'.join(lines[start:end+1]);(b/'function.sh').write_text(function+'\n')
 # Preserve a reproduction of the original omitted-flag call, using exactly the
 # same real artifact/source identity and native CLI as the repaired function.
 direct=[str(R),str(root/'src/5_run_benchmark_methods/validate_benchmark_rds_contract.R'),'--artifact',str(artifact),'--dataset','Joanito','--view','batch_effect_uncorrected','--method','pseudobulk','--config',str(config),'--input-root',str(b/'scratch'),'--batch-pass','uncorrected','--source-identity',str(identity),'--source-identity-verified']
 old=subprocess.run(direct,env=env,text=True,capture_output=True);assert old.returncode!=0;assert 'explicit final batch context' in old.stderr,old.stderr
 new=subprocess.run(direct+['--analysis-variant','final'],env=env,text=True,capture_output=True);assert new.returncode==0,new.stdout+new.stderr
 script=b/'run.sh';script.write_text('''#!/bin/bash
set -euo pipefail
source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"
source "$FIXTURE/function.sh"
ANALYSIS_VARIANT=final;RDS_LABELS=(pseudobulk);MANIFEST="$ECODA_RUN_ROOT/manifests/selection.tsv"
SOURCE_IDENTITY="$ECODA_RUN_ROOT/manifests/source_identity.json";SCRIPT_DIR="$PROJECT_ROOT/src/5_run_benchmark_methods"
PIXI_RSCRIPT="$REAL_R"
benchmark_artifacts_for() { [[ "$1" == Joanito && "$2" == batch_effect_uncorrected && "$3" == pseudobulk ]];ARTIFACT_PATHS=("$REAL_ARTIFACT"); }
stage5_validate_final_rds_rows
''')
 env.update(FIXTURE=str(b),REAL_R=str(R),REAL_ARTIFACT=str(artifact));before={p:p.read_bytes() for p in [artifact,pathlib.Path(str(artifact)+'.md5'),source,pathlib.Path(str(source)+'.md5'),identity,binding]}
 good=subprocess.run(['bash',str(script)],env=env,text=True,capture_output=True);assert good.returncode==0,good.stdout+good.stderr
 path=m/'rds_validation_final_pseudobulk.tsv';stdout=path.with_suffix('.stdout');stderr=path.with_suffix('.stderr');assert not path.exists();assert stdout.is_file() and stderr.is_file();assert 'contract OK' in stdout.read_text()
 # A real checksum failure retains exact selected list, proof, stdout and stderr.
 original=artifact.read_bytes();artifact.write_bytes(original+b'bad')
 bad=subprocess.run(['bash',str(script)],env=env,text=True,capture_output=True);assert bad.returncode!=0;assert path.is_file() and pathlib.Path(str(path)+'.md5').is_file() and stdout.is_file() and stderr.is_file();assert str(path) in bad.stderr and 'exit 1' in bad.stderr
 assert 'checksum' in stderr.read_text().lower(),stderr.read_text();assert path.read_text()==f'{artifact}\tpseudobulk\tJoanito\tbatch_effect_uncorrected\t0\n'
 artifact.write_bytes(original);assert all(p.read_bytes()==v for p,v in before.items());assert source.read_bytes().startswith(b'NOT HDF5')
 print('PASS actual native derivative RDS CLI: missing-final red, explicit-final green, real selected checker and preserved failure list/stdout/stderr; no H5AD/model/Slurm')
PY
