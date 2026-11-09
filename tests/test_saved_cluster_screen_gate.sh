#!/bin/bash
# Actual obs-only R runner/common/artifact validators; Slurm/runtime and scores stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-saved-cluster.XXXXXX")";TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
TREE="${TMP}/snapshot/tree";SCRATCH="${TMP}/scratch"
mkdir -p "${TREE}" "${SCRATCH}/_ecoda_inputs/fixture" "${TMP}/bin" "${TMP}/nas"
cp -R "${ROOT}/src" "${TREE}/src"
cp "${ROOT}/datasets.json" "${TREE}/datasets.json"
printf 'immutable source fixture\n' > "${TREE}/source.manifest"
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_TMP="${TMP}" TEST_NAS="${TMP}/nas"
export TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_REAL_R="${ROOT}/.pixi/envs/default/bin/Rscript"
export TEST_RSCRIPT="${TMP}/bin/Rscript"
export ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_MODE=host PATH="${TMP}/bin:${PATH}"
cat > "${TREE}/src/slurm_config.sh" <<'SH'
export PROJECT_ROOT="${TEST_TREE}" DATASETS_JSON_FILE="${TEST_TREE}/datasets.json"
export HPC_SCRATCH_DIR="${TEST_SCRATCH}" NAS_TARGET_DIR="${TEST_NAS}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners"
export PYTHON_BIN="${TEST_PYTHON}" PIXI_RSCRIPT="${TEST_RSCRIPT}"
SLURM_PARTITION=shared-cpu
SH
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'SH'
ecoda_runtime_validate_submission() { [[ "$1" == apptainer ]]; printf 'runtime fixture\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run() { [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]] && cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"; }
ecoda_runtime_export_csv() { [[ "$1" == stage5 && "$2" == 0 ]]; printf 'ECODA_RUNTIME_PROFILE=stage5'; }
ecoda_runtime_reexec_worker() { [[ "$1" == stage5 && "$3" == --worker && "$2" == "${TEST_TREE}/src/utils/bash/saved_cluster_screen.sh" ]]; ecoda_runtime_validate_bound_run; }
SH
cat > "${TMP}/bin/sbatch" <<'SH'
#!/bin/bash
set -euo pipefail
[[ "$*" == *'--parsable --wait'* && "$*" == *'--cpus-per-task=1'* && "$*" == *'--mem=64G'* && "$*" == *'--time=02:00:00'* && "$*" != *'--gpus'* && "$*" != *'--array'* ]]
printf '%s\n' "$*" >> "${TEST_SCRATCH}/scheduler.calls"
args=("$@");spool="${ECODA_RUN_ROOT}/spool";mkdir -p "${spool}"
cp "${args[${#args[@]}-2]}" "${spool}/slurm_script"
# Compute node has no NAS mount. Its bound identity must not create a mount.
mv "${TEST_NAS}" "${TEST_NAS}.hidden"
set +e
SLURM_JOB_ID=654321 bash "${spool}/slurm_script" --worker > "${ECODA_RUN_ROOT}/logs/saved_cluster_654321.log" 2> "${ECODA_RUN_ROOT}/logs/saved_cluster_654321.err"
rc=$?
set -e
[[ ! -e "${TEST_NAS}" ]];mv "${TEST_NAS}.hidden" "${TEST_NAS}"
if [[ ${rc} -ne 0 ]]; then cat "${ECODA_RUN_ROOT}/logs/saved_cluster_654321.err" >&2; fi
printf '654321;fixture\n'
exit "${rc}"
SH
cat > "${TMP}/bin/Rscript" <<'SH'
#!/bin/bash
set -euo pipefail
if [[ "${2:-}" == *'/run_local_ecoda_derived.R' ]]; then
  original="$2";shift 2
  exec "${TEST_REAL_R}" --vanilla "${TEST_TMP}/runner_fixture.R" "${original}" "$@"
fi
exec "${TEST_REAL_R}" "$@"
SH
cat > "${TMP}/runner_fixture.R" <<'R'
a<-commandArgs(TRUE);original<-a[1];source(original)
.ecoda_runner_path<-function()normalizePath(original,mustWork=TRUE)
source<-function(file,...) {
  result<-base::source(file,...)
  if(basename(file)=='scoring_metrics.R')assign('calc_sep_score',function(...)list(
    anosim_score=.5,mod_knn3_score=.5,cluster_score=.5,sil_score=.5,lisi_score=.5),envir=.GlobalEnv)
  result
}
args<-.ecoda_parse_flags(a[-1])
if(Sys.getenv('TEST_FAIL_HARMONY')=='1' && args$analysis=='harmony')stop('injected scientific worker failure')
.ecoda_run(args)
R
chmod +x "${TMP}/bin/sbatch" "${TMP}/bin/Rscript"
for name in sacct squeue rsync; do printf '#!/bin/bash\necho forbidden >> "%s"\nexit 99\n' "${TMP}/forbidden" > "${TMP}/bin/${name}";chmod +x "${TMP}/bin/${name}";done
mkdir -p "${TMP}/native9ca/src/utils/bash"
for name in h5ad_preflight_worker.sh ecoda_run_common.sh; do
  git -C "${ROOT}" show "9ca2453eea7320a721ef6dfd458f9db108abecac:src/utils/bash/${name}" > "${TMP}/native9ca/src/utils/bash/${name}"
done
"${TEST_PYTHON}" - <<'PY'
import hashlib,json,os,pathlib,copy
import h5py,numpy as np
tmp=pathlib.Path(os.environ['TEST_TMP']);scratch=pathlib.Path(os.environ['TEST_SCRATCH']);tree=pathlib.Path(os.environ['TEST_TREE']);folder=scratch/'_ecoda_inputs/fixture'
config=json.loads((tree/'datasets.json').read_text());cohorts=[ds for ds,row in config.items() if not ds.startswith('_') and row.get('use_for_benchmark') is True and 'benchmark_analysis' in row.get('views',{})]
assert len(cohorts)==11
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def stat(path):
 s=path.stat();return {k:getattr(s,k) for k in ('st_dev','st_ino','st_size','st_mtime_ns','st_ctime_ns','st_mode')}
def save(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value));return dict(path=str(path),sha256=sha(path))
def evidence(path,fields):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(''.join(k+'='+v+'\n' for k,v in fields.items()));return dict(path=str(path),sha256=sha(path),stat=stat(path))
parent=scratch/'_ecoda_runs/prior_preflight';witness=[];references=[];preflight=[]
for i,ds in enumerate(cohorts,1):
 columns=config[ds]['columns'];path=scratch/ds/'output'/config[ds]['views']['benchmark_analysis']['output_file_name'];path.parent.mkdir(parents=True)
 n=53 if ds=='Stephenson' else 26 if ds=='Bassez' else 6;ids=[ds+str(k) for k in range(n)];labels=['A' if k%2==0 else 'B' for k in range(n)]
 with h5py.File(path,'w') as handle:
  group=handle.create_group('obs')
  for name,values in [(columns['sample'],ids),(columns['label'],labels),('leiden_res_50_benchmark_analysis_hvg2000',[str(k%3) for k in range(n)]),('leiden_res_2_benchmark_analysis_hvg2000_harmony',[str(k%2) for k in range(n)])]:
   group.create_dataset(name,data=np.asarray(values,dtype=h5py.string_dtype()))
 md5=hashlib.md5(path.read_bytes()).hexdigest();size=str(path.stat().st_size)
 proof=evidence(pathlib.Path(str(path)+'.md5'),dict(MD5=md5,SIZE=size,PATH=str(path)))
 record=evidence(parent/'manifests/artifacts'/(hashlib.sha256(str(path).encode()).hexdigest()[:32]+'.record'),dict(PATH=str(path),MD5=md5,SIZE=size,RUN_ID=parent.name,PRODUCER='stage5_preflight',STATE='PUBLISHED'))
 status=evidence(parent/'status/h5ad_preflight'/(ds+'__benchmark_analysis.status'),dict(STATE='OK',RUN_ID=parent.name,DATASET=ds,VIEW='benchmark_analysis',TASK_ID=str(i)))
 # Actual current-label witness carries record stats, but only proof/status SHA.
 witness.append(dict(dataset=ds,view='benchmark_analysis',path=str(path),expected_source_md5=md5,size=int(size),stat_before=stat(path),stat_after=stat(path),sample_ids=ids,n_samples=n,named_labels=dict(zip(ids,labels)),source_identity_order_matches_observed=True,configured_sample_column=columns['sample'],label_column=columns['label'],proof={k:proof[k] for k in ('path','sha256')},preflight_record=record,preflight_status={k:status[k] for k in ('path','sha256')}))
 references.append('\t'.join([ds,'benchmark_analysis',str(path),md5,size,proof['path']])+'\n');preflight.append('\t'.join([ds,'benchmark_analysis',str(path)])+'\n')
refs=folder/'references.tsv';refs.write_text(''.join(references));pre=parent/'manifests/h5ad_preflight.tsv';pre.write_text(''.join(preflight))
w=dict(schema='ecoda_old56_current_benchmark_labels_observation_v1',original_source_commit='9ca2453eea7320a721ef6dfd458f9db108abecac',scientific_reuse_approved=False,expression_values_read=False,full_h5ad_hashes=0,fits=0,scores_computed=0,normalization_policy='exact_text_storage_v2_no_numeric_or_logical_conversion',collector_sha256='a'*64,controls=[dict(path=str(tree/'datasets.json'),sha256=sha(tree/'datasets.json'))],rows=witness)
wref=save(folder/'witness.json',w);lr=save(folder/'label_review.json',dict(schema='ecoda_current_primary_label_witness_review_v1',status='APPROVED_CURRENT_PRIMARY_LABEL_WITNESS',observation_sha256=wref['sha256'],scientific_result_reuse_approved=False,collector_sha256=w['collector_sha256'],datasets=11))
r=dict(schema='ecoda_saved_cluster_input_review_v1',status='APPROVED_CURRENT_INPUT_CONTENT_WITNESS_ONLY',review_complete=True,reviewer='synthetic fixture input-only reviewer',scientifically_approved_rows=0,formal_pipeline_dependency=False,parent_gate_approved=False,scheduler_exit_established=False,parent_gate_state='PRELAUNCH_STOP',wrapper_exit_code=None,input_references_sha256=sha(refs),observation=wref,label_review=lr,preflight_manifest=dict(path=str(pre),sha256=sha(pre)))
for key,name in [('native_preflight_worker','h5ad_preflight_worker.sh'),('native_checksum_helper','ecoda_run_common.sh')]:
 path=tmp/'native9ca/src/utils/bash'/name;r[key]=dict(path=str(path),sha256=sha(path))
save(folder/'review.json',r);(folder/'datasets.csv').write_text(','.join(cohorts))
for name,key,value in [('pending','status','PENDING'),('science','scientifically_approved_rows',1),('dependency','formal_pipeline_dependency',True)]:
 bad=copy.deepcopy(r);bad[key]=value;save(folder/(name+'.json'),bad)
PY
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
INPUTS="${SCRATCH}/_ecoda_inputs/fixture"
DATASETS="$(cat "${INPUTS}/datasets.csv")"
invoke() {
  local id="$1" review="$2"
  ECODA_RUN_ID="${id}" ECODA_RUN_ROOT="${SCRATCH}/_ecoda_runs/${id}" bash "${TREE}/src/utils/bash/saved_cluster_screen.sh" \
    --input-references "${INPUTS}/references.tsv" --input-references-sha256 "$(ecoda_sha256_file "${INPUTS}/references.tsv")" \
    --input-review "${review}" --input-review-sha256 "$(ecoda_sha256_file "${review}")" --datasets "${DATASETS}"
}
reject() { if "$@" > "${TMP}/negative.log" 2>&1; then echo 'Unexpected success' >&2;exit 1;fi; }
invoke good "${INPUTS}/review.json" > "${TMP}/good.log" 2>&1 || { tail -40 "${TMP}/good.log" >&2;exit 1; }
RUN="${SCRATCH}/_ecoda_runs/good"
[[ "$(wc -l < "${RUN}/manifests/selection.tsv" | tr -d ' ')" == 22 && "$(wc -l < "${RUN}/manifests/outputs.tsv" | tr -d ' ')" == 30 ]]
ECODA_RUN_ID=good ECODA_RUN_ROOT="${RUN}" bash "${TREE}/src/utils/bash/saved_cluster_screen.sh" --validate > "${TMP}/validated.log"
grep -q '^SAVED_CLUSTER_SCREEN_VALIDATED=good$' "${TMP}/validated.log"
# Real new purpose-handler branch, actual child validator and receipt closure.
eval "$(sed -n '/^_audit_analysis_release() {/,/^}/p' "${ROOT}/src/utils/bash/ecoda_run_audit.sh")"
_audit_metadata_value() { sed -n "s/^$2=//p" "$1"; }
_audit_metadata_count() { awk -F= -v k="$2" '$1==k{n++}END{print n+0}' "$1"; }
_audit_regular_file() { [[ -f "$1" && ! -L "$1" ]]; }
_audit_die() { echo "$*" >&2;return 1; }
_audit_scheduler_ids() { [[ "$(cat "${RUN_ROOT_REAL}/manifests/scheduler_ids.tsv")" == "$(printf 'ARRAY\t654321')" ]]; }
export ECODA_RUN_ID=good ECODA_RUN_ROOT="${RUN}"
RUN_ROOT_REAL="${RUN}" RUN_ID=good SELECTION_ARG="${RUN}/manifests/selection.tsv" AUDIT_TERMINAL_STATE=OK
AUDIT_SOURCE_ROOT="${TREE}" SOURCE_MANIFEST_ARG="${TREE}/source.manifest"
_audit_analysis_release > "${TMP}/generic-audit.log"
grep -q '^SAVED_CLUSTER_SCREEN_VALIDATED=good$' "${TMP}/generic-audit.log"
# A syntactically valid owner from another stage is not this Stage 5 output.
path="$(head -1 "${RUN}/manifests/outputs.tsv")";owner="$(ecoda_artifact_owner_dir "${path}")/owner"
cp "${owner}" "${TMP}/owner.original"
sed 's/^STAGE=stage5$/STAGE=stage3/' "${owner}" > "${TMP}/owner.wrong-stage"
cat "${TMP}/owner.wrong-stage" > "${owner}"
reject bash "${TREE}/src/utils/bash/saved_cluster_screen.sh" --validate
cmp -s "${owner}" "${TMP}/owner.wrong-stage"
cat "${TMP}/owner.original" > "${owner}"
before="$(ecoda_sha256_file "${RUN}/status/terminal")"
printf 'STATE=FAIL\nRUN_ID=good\n' > "${RUN}/status/terminal"
failed_before="$(ecoda_sha256_file "${RUN}/status/terminal")"
reject bash "${TREE}/src/utils/bash/saved_cluster_screen.sh" --validate
[[ "$(ecoda_sha256_file "${RUN}/status/terminal")" == "${failed_before}" ]]
printf 'STATE=OK\nRUN_ID=good\n' > "${RUN}/status/terminal"
[[ "$(ecoda_sha256_file "${RUN}/status/terminal")" == "${before}" ]]
for name in pending science dependency; do reject invoke "bad_${name}" "${INPUTS}/${name}.json";done
( DATASETS="${DATASETS#*,}"; reject invoke bad_subset "${INPUTS}/review.json" )
reject invoke good "${INPUTS}/review.json"
TEST_FAIL_HARMONY=1 invoke partial "${INPUTS}/review.json" > "${TMP}/partial.log" 2>&1 && { echo 'Expected worker failure' >&2;exit 1; }
[[ -s "${SCRATCH}/_ecoda_runs/partial/manifests/scheduler_ids.tsv" ]]
grep -q '^STATE=FAIL$' "${SCRATCH}/_ecoda_runs/partial/status/terminal"
path="${SCRATCH}/_ecoda_runs/partial/artifacts/saved_labels/res50/Adams_ECODA_seuratres_50.rds"
owner="$(ecoda_artifact_owner_dir "${path}")";grep -q '^STATE=OK$' "${owner}/owner"
[[ "$(wc -l < "${SCRATCH}/scheduler.calls" | tr -d ' ')" == 2 && ! -e "${TMP}/forbidden" ]]
# A changed source stat invalidates the prior full-content witness; no H5 hash or job.
path="$(head -1 "${INPUTS}/references.tsv" | cut -f3)";touch "${path}"
reject invoke changed_input "${INPUTS}/review.json"
[[ "$(wc -l < "${SCRATCH}/scheduler.calls" | tr -d ' ')" == 2 ]]
echo 'PASS saved22 exact11 refs/two native obs producers/one scalar/hiddenNAS/ownership/terminal/purpose/failure preservation'
