"""Real bound ownership for physical parameter arrays, without bypass flags."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.utils.py.derivative_artifacts import bind
from src.utils.py.benchmark_request_scope import read_request,normalized_rows


def proof(path):
    Path(str(path)+'.md5').write_text(f'MD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\nPATH={path}\n')


def run():
    with tempfile.TemporaryDirectory(prefix='ecoda-parameter-owners-') as temporary:
        base=Path(temporary).resolve();scratch=base/'scratch';nas=base/'nas';nas.mkdir()
        run=scratch/'_ecoda_runs/parameter-wave';manifests=run/'manifests';manifests.mkdir(parents=True)
        (run/'metadata').write_text('STAGE=stage5\nRUN_ID=parameter-wave\nSTATE=ACTIVE\n')
        config=base/'datasets.json'
        config.write_text(json.dumps({dataset:{'columns':{'sample':'Sample','label':'biology'},
            'views':{'benchmark_analysis':{'output_file_name':'source.h5ad'}}} for dataset in ('Stephenson','Bassez')}))
        for name in ('source.manifest','runtime.identity'):(manifests/name).write_text('frozen fixture identity\n')
        selected=('hvg2000_highres_dims15','hvg1000_highres_dims15','hvg2000_highres_dims2')
        requested=manifests/'requested_selection.tsv'
        requested.write_text('Stephenson\tbenchmark_analysis\tgloscope\nStephenson\tbenchmark_analysis\tscpoli\n'+
            ''.join(f'Bassez\tbenchmark_analysis\tscpoli\t{token}\n' for token in selected))
        rows=normalized_rows(read_request(requested,json.loads(config.read_text())))
        selection=manifests/'selection.tsv';selection.write_text(''.join('\t'.join(row)+'\n' for row in rows))
        source=base/'source.h5ad';source.write_bytes(b'input identity metadata only');proof(source)
        refs=base/'inputs.tsv';refs.write_text(''.join(f'{ds}\tbenchmark_analysis\t{source}\t{hashlib.md5(source.read_bytes()).hexdigest()}\t{source.stat().st_size}\t{source}.md5\n' for ds in ('Stephenson','Bassez')))
        binding=bind(run_root=run,release_id='test-release',nas_root=nas,stage='stage5',config_path=config,selection=selection,input_artifacts=refs)
        functions=[]
        for relative,names in (
            ('src/5_run_benchmark_methods/stage5_dispatcher.sh',('stage5_validate_output_ownership','benchmark_method_shard_rows')),
            ('src/5_run_benchmark_methods/matrix_watchdog.sh',('matrix_validate_output_ownership',))):
            lines=(ROOT/relative).read_text().splitlines()
            for name in names:
                start=lines.index(name+'() {');end=next(i for i in range(start+1,len(lines)) if lines[i]=='}')
                functions.extend(lines[start:end+1])
        (base/'functions.sh').write_text('\n'.join(functions)+'\n')
        (base/'exercise.sh').write_text(r'''#!/bin/bash
set -euo pipefail
source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"
source "$PROJECT_ROOT/src/utils/bash/ecoda_stage5_policy.sh"
source "$FIXTURE/functions.sh"
RUN_ID="$ECODA_RUN_ID";RUN_ROOT="$ECODA_RUN_ROOT";PASS_ARG="";ANALYSIS_PASS="";ANALYSIS_VARIANT=""
export ANALYSIS_ROOT="$ECODA_RUN_ROOT/artifacts/benchmark"
export ANALYSIS_NAS_ROOT="$NAS_TARGET_DIR/releases/test-release/$RUN_ID/benchmark"
mkdir -p "$ANALYSIS_ROOT" "$ANALYSIS_NAS_ROOT"
# Exact real dispatcher default physical rows:5 distinct GloScope parameters.
gloscope="$ECODA_RUN_ROOT/manifests/gloscope.tsv"
while IFS=$'\t' read -r token resource; do
  printf 'Stephenson\tbenchmark_analysis\tgloscope\t%s\n' "$token"
done < <(benchmark_method_shard_rows gloscope benchmark_analysis) > "$gloscope"
ecoda_write_checksum "$gloscope"
original="$(ecoda_sha256_file "$gloscope")"
[[ "$(wc -l < "$gloscope" | tr -d ' ')" == 5 ]]
stage5_validate_output_ownership "$gloscope"
[[ ${#ECODA_OUTPUT_PATHS[@]} == 2 && ${#ECODA_OUTPUT_OWNER_DIRS[@]} == 2 ]]
[[ "${ECODA_OUTPUT_PATHS[0]}" == "$ANALYSIS_ROOT/results/Stephenson_gloscope.rds" ]]
for path in "${ECODA_OUTPUT_PATHS[@]}"; do
  ecoda_artifact_owner_validate "$path" "$RUN_ID" >/dev/null
  [[ "$ECODA_ARTIFACT_OWNER_STATE" == ACTIVE && "$ECODA_ARTIFACT_OWNER_STAGE" == stage5 ]]
done
# The real compute retry guard uses scratch only, with no mounted NAS.
mv "$NAS_TARGET_DIR" "$NAS_TARGET_DIR.hidden"
matrix_validate_output_ownership "$gloscope"
[[ ${#ECODA_OUTPUT_PATHS[@]} == 1 ]]
[[ ! -e "$NAS_TARGET_DIR" ]]
mv "$NAS_TARGET_DIR.hidden" "$NAS_TARGET_DIR"
[[ "$(ecoda_sha256_file "$gloscope")" == "$original" ]]
ecoda_validate_manifest "$gloscope" 4;ecoda_validate_checksum "$gloscope"
# Full Stephenson + selected Bassez tokens expand their bound family unions once.
scpoli="$ECODA_RUN_ROOT/manifests/scpoli.tsv"
while IFS=$'\t' read -r token resource; do printf 'Stephenson\tbenchmark_analysis\tscpoli\t%s\n' "$token"; done < <(benchmark_method_shard_rows scpoli benchmark_analysis) > "$scpoli"
printf 'Bassez\tbenchmark_analysis\tscpoli\thvg2000_highres_dims15\nBassez\tbenchmark_analysis\tscpoli\thvg1000_highres_dims15\nBassez\tbenchmark_analysis\tscpoli\thvg2000_highres_dims2\n' >> "$scpoli"
ecoda_write_checksum "$scpoli";original_scpoli="$(ecoda_sha256_file "$scpoli")"
stage5_validate_output_ownership "$scpoli"
[[ ${#ECODA_OUTPUT_PATHS[@]} == 22 && ${#ECODA_OUTPUT_OWNER_DIRS[@]} == 22 ]]
for path in "${ECODA_OUTPUT_PATHS[@]}"; do
  ecoda_artifact_owner_validate "$path" "$RUN_ID" >/dev/null
  [[ "$ECODA_ARTIFACT_OWNER_STATE" == ACTIVE ]]
  case "$path" in *Bassez*lowres*) echo 'unrequested Bassez LR owner' >&2;exit 1;; esac
done
matrix_validate_output_ownership "$scpoli"
[[ ${#ECODA_OUTPUT_PATHS[@]} == 11 ]]
[[ "$(ecoda_sha256_file "$scpoli")" == "$original_scpoli" ]]
ecoda_validate_manifest "$scpoli" 4;ecoda_validate_checksum "$scpoli"
reject() { if "$@" > "$FIXTURE/reject.log" 2>&1; then echo "unexpected acceptance: $*" >&2;exit 1;fi; }
# Both entrypoints retain strict physical-row validation and3-column conflicts.
for function in stage5_validate_output_ownership matrix_validate_output_ownership; do
  bad="$ECODA_RUN_ROOT/manifests/bad.tsv"
  cat "$gloscope" > "$bad";head -1 "$gloscope" >> "$bad";reject "$function" "$bad"
  printf 'Stephenson\tbenchmark_analysis\tgloscope\nStephenson\tbenchmark_analysis\tgloscope\n' > "$bad";reject "$function" "$bad"
  grep -q 'duplicate artifact path' "$FIXTURE/reject.log"
  printf 'Stephenson\tbenchmark_analysis\tgloscope\t\n' > "$bad";reject "$function" "$bad"
  printf 'Stephenson\tbenchmark_analysis\tgloscope\tx\textra\n' > "$bad";reject "$function" "$bad"
  head -1 "$gloscope" > "$bad";printf 'Stephenson\tbenchmark_analysis\tgloscope\n' >> "$bad";reject "$function" "$bad"
  [[ ! -e "$bad.ownership.$$" ]]
done
# Foreign current writer still blocks the exact real artifact owner path.
path="$ANALYSIS_ROOT/results/Stephenson_gloscope.rds"
owner="$(ecoda_artifact_owner_dir "$path")/owner";cp "$owner" "$FIXTURE/owner.saved"
sed 's/RUN_ID=parameter-wave/RUN_ID=foreign_writer/' "$FIXTURE/owner.saved" > "$owner"
reject stage5_validate_output_ownership "$gloscope"
reject matrix_validate_output_ownership "$gloscope"
cp "$FIXTURE/owner.saved" "$owner"
stage5_validate_output_ownership "$gloscope"
[[ "$(ecoda_sha256_file "$gloscope")" == "$original" ]]
printf 'PASS: actual bound5-GloScope/11-scPoli physical rows, unique family owners, unchanged manifest checksums; duplicate/malformed/foreign-owner rejection\n'
''')
        env=dict(os.environ,PROJECT_ROOT=str(ROOT),DATASETS_JSON_FILE=str(config),PYTHON_BIN=sys.executable,
            HPC_SCRATCH_DIR=str(scratch),NAS_TARGET_DIR=str(nas),ECODA_RUNS_ROOT=str(scratch/'_ecoda_runs'),
            ECODA_OWNERS_ROOT=str(scratch/'_ecoda_owners'),ECODA_RUN_ROOT=str(run),ECODA_RUN_ID=run.name,
            ECODA_ARTIFACT_BINDING=str(binding),FIXTURE=str(base))
        env.pop('BENCHMARK_MATRIX_TEST',None)
        result=subprocess.run(['bash',str(base/'exercise.sh')],env=env,capture_output=True,text=True)
        assert result.returncode==0,(result.returncode,result.stdout,result.stderr)
        print(result.stdout.strip())

if __name__=='__main__':run()
