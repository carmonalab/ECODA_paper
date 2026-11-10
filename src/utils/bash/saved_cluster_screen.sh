#!/bin/bash
# One obs-only saved-label job:11 cohorts x(res50,Harmony2), no model refits.
set -euo pipefail
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && "${ECODA_SOURCE_ROOT:-}" == /* ]] || exit 1
SCRIPT_DIR="${ECODA_SOURCE_ROOT%/}/src/utils/bash"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
SELF="${SCRIPT_DIR}/saved_cluster_screen.sh"
RUNNER="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_local_ecoda_derived.R"
PURPOSE=saved_cluster_counts_screen
PRODUCER=saved_cluster_screen
MODE=controller
if [[ "${1:-}" == --worker || "${1:-}" == --validate ]]; then MODE="${1#--}"; shift; fi
REFERENCES=""; REFERENCES_SHA=""; REVIEW=""; REVIEW_SHA=""; DATASETS=""
while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || exit 1
  case "$1" in
    --input-references) [[ -z "${REFERENCES}" ]] || exit 1; REFERENCES="$2" ;;
    --input-references-sha256) [[ -z "${REFERENCES_SHA}" ]] || exit 1; REFERENCES_SHA="$2" ;;
    --input-review) [[ -z "${REVIEW}" ]] || exit 1; REVIEW="$2" ;;
    --input-review-sha256) [[ -z "${REVIEW_SHA}" ]] || exit 1; REVIEW_SHA="$2" ;;
    --datasets) [[ -z "${DATASETS}" ]] || exit 1; DATASETS="$2" ;;
    *) echo 'ERROR: unknown saved-cluster argument' >&2; exit 1 ;;
  esac
  shift 2
done
ecoda_validate_run_id "${ECODA_RUN_ID:-}" || exit 1
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
ecoda_require_source_script_path "${RUNNER}" "${ECODA_SOURCE_ROOT}" >/dev/null
SCRATCH="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
RUN="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${RUN}" ]] || exit 1
if [[ "${MODE}" == controller ]]; then
  [[ ! -e "${RUN}" && ! -L "${RUN}" && "${REFERENCES}" == "${SCRATCH}/_ecoda_inputs/"* &&
     "${REVIEW}" == "${SCRATCH}/_ecoda_inputs/"* && "${REFERENCES_SHA}" =~ ^[a-f0-9]{64}$ &&
     "${REVIEW_SHA}" =~ ^[a-f0-9]{64}$ && -n "${DATASETS}" ]] || exit 1
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  ecoda_atomic_write "${RUN}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
else
  [[ -z "${REFERENCES}${REFERENCES_SHA}${REVIEW}${REVIEW_SHA}${DATASETS}" ]] || exit 1
  ecoda_open_run "${ECODA_RUN_ID}" stage5
fi
export ECODA_RUN_ROOT="${RUN}" ECODA_RUNTIME_PROFILE=stage5 RETICULATE_PYTHON="${PYTHON_BIN}"
M="${RUN}/manifests"; OUTPUT="${RUN}/artifacts/saved_labels"
STATUS="${RUN}/status/saved_cluster.status"
SCOPE_FILES=(source.manifest runtime.identity operation selection.tsv inputs.tsv input_references.tsv input_review.json input_witness.json label_review.json native_worker.sh native_common.sh preflight.tsv scope.json)
failed() {
  local path owner
  ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\n" || true
  if [[ "${MODE}" == controller ]]; then ecoda_set_run_state FAIL saved_cluster_failed || true; fi
  if [[ -f "${M}/outputs.tsv" ]]; then
    while IFS= read -r path; do
      owner="$(ecoda_artifact_owner_dir "${path}")" || continue
      if _ecoda_artifact_owner_validate_dir "${owner}" "${path}" &&
         [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]; then
        ecoda_artifact_owner_set_state "${path}" FAIL saved_cluster_failed >/dev/null || true
      fi
    done < "${M}/outputs.tsv"
  fi
}
if [[ "${MODE}" != validate ]]; then trap failed EXIT; fi

metadata() {
  "${PYTHON_BIN}" - "$1" "${RUN}" "${PROJECT_ROOT}" "${REFERENCES}" "${REFERENCES_SHA}" "${REVIEW}" "${REVIEW_SHA}" "${DATASETS}" <<'PY'
import hashlib,json,os,pathlib,re,sys,tempfile
mode,root,project,refs,refs_sha,review_path,review_sha,datasets=sys.argv[1:]
root=pathlib.Path(root);m=root/'manifests';project=pathlib.Path(project)
observed_stats={}
def need(ok,text):
    if not ok:raise SystemExit('ERROR: '+text)
def digest(raw):return hashlib.sha256(raw).hexdigest()
def regular(path):
    p=pathlib.Path(path);need(p.is_absolute() and p.is_file() and not p.is_symlink() and p.resolve()==p,'noncanonical metadata/input');return p
def small(path):
    p=regular(path);stat(p);need(p.suffix.lower() not in ('.h5','.h5ad','.rds','.sif') and p.stat().st_size<=2*1024*1024,'bounded metadata only');raw=p.read_bytes();stat(p);return raw
def read_ref(ref):
    need(set(ref)=={'path','sha256'} and re.fullmatch('[a-f0-9]{64}',ref['sha256']),'invalid evidence ref')
    raw=small(ref['path']);need(digest(raw)==ref['sha256'],'evidence changed');return raw
def put(name,raw):
    fd,tmp=tempfile.mkstemp(prefix='.saved-scope-',dir=m)
    try:
        with os.fdopen(fd,'wb') as stream:stream.write(raw);stream.flush();os.fsync(stream.fileno())
        os.link(tmp,m/name)
    finally:os.unlink(tmp)
def fields(raw):
    pairs=[line.split('=',1) for line in raw.decode().splitlines()]
    need(pairs and all(len(r)==2 for r in pairs) and len({r[0] for r in pairs})==len(pairs),'duplicate/malformed native fields');return dict(pairs)
def stat(path):
    p=regular(path);s=p.stat();value={k:getattr(s,k) for k in ('st_dev','st_ino','st_size','st_mtime_ns','st_ctime_ns','st_mode')}
    prior=observed_stats.setdefault(str(p),value);diff={k:dict(expected=prior[k],actual=value[k]) for k in value if prior[k]!=value[k]}
    need(not diff,'in-process stat changed path='+str(p)+' fields='+json.dumps(diff,sort_keys=True));return value
def selected_beegfs(path):
    scratch=root.parent.parent;p=pathlib.Path(path);need(p.is_relative_to(scratch),'selected path escapes scratch: '+str(p))
    mounts=[]
    for line in pathlib.Path('/proc/self/mountinfo').read_text().splitlines():
        fields=line.split();separator=fields.index('-');point=pathlib.Path(fields[4].replace('\\040',' '))
        mounts.append((point,fields[separator+1]))
    def covering(p):return max((r for r in mounts if p==r[0] or r[0] in p.parents),key=lambda r:len(str(r[0])))
    mount=covering(p);need(mount==covering(scratch) and mount[1]=='beegfs','selected path is not on expected BeeGFS mount: '+str(p))
def portable_stat(actual,expected,path):
    # Device numbers are kernel-local. Keep full facts and strict in-process
    # checks, but compare content witnesses portably across login/worker nodes.
    selected_beegfs(path)
    diff={k:dict(expected=expected.get(k),actual=actual.get(k)) for k in set(actual)|set(expected) if k!='st_dev' and actual.get(k)!=expected.get(k)}
    need(set(actual)==set(expected)=={'st_dev','st_ino','st_size','st_mode','st_mtime_ns','st_ctime_ns'} and all(type(v) is int for v in expected.values()) and not diff,'selected stat changed path='+str(path)+' fields='+json.dumps(diff,sort_keys=True));return True
if mode=='bind':
    raw=small(refs);need(digest(raw)==refs_sha,'input references SHA differs');put('input_references.tsv',raw)
    raw=small(review_path);need(digest(raw)==review_sha,'input review SHA differs');put('input_review.json',raw)
review=json.loads(small(m/'input_review.json'))
need(review['schema']=='ecoda_saved_cluster_input_review_v1' and review['status']=='APPROVED_CURRENT_INPUT_CONTENT_WITNESS_ONLY'
     and review['review_complete'] is True and isinstance(review['reviewer'],str) and review['reviewer'].strip()
     and review['scientifically_approved_rows']==0 and review['formal_pipeline_dependency'] is False
     and review['parent_gate_approved'] is False and review['scheduler_exit_established'] is False
     and review['parent_gate_state']=='PRELAUNCH_STOP' and review['wrapper_exit_code'] is None,
     'unapproved input-content-only review')
for key,name in (('observation','input_witness.json'),('label_review','label_review.json'),('native_preflight_worker','native_worker.sh'),('native_checksum_helper','native_common.sh'),('preflight_manifest','preflight.tsv')):
    raw=read_ref(review[key])
    if mode=='bind':put(name,raw)
    else:need(raw==small(m/name),'frozen evidence differs')
w=json.loads(small(m/'input_witness.json'));label_review=json.loads(small(m/'label_review.json'))
need(w['schema']=='ecoda_old56_current_benchmark_labels_observation_v1' and w['scientific_reuse_approved'] is False
     and w['original_source_commit']=='9ca2453eea7320a721ef6dfd458f9db108abecac'
     and w['expression_values_read'] is False and w['full_h5ad_hashes']==w['fits']==w['scores_computed']==0
     and w['normalization_policy']=='exact_text_storage_v2_no_numeric_or_logical_conversion','unsupported prior fullhash/label witness')
need(label_review['schema']=='ecoda_current_primary_label_witness_review_v1' and label_review['status']=='APPROVED_CURRENT_PRIMARY_LABEL_WITNESS'
     and label_review['observation_sha256']==review['observation']['sha256'] and label_review['scientific_result_reuse_approved'] is False
     and label_review['collector_sha256']==w['collector_sha256'] and label_review['datasets']==11,'unreviewed current label witness')
need(review['native_preflight_worker']['sha256']=='352a601a435136101b07f4161111543885cb89c10851a5e4d0f25249796c147f'
     and review['native_checksum_helper']['sha256']=='0ad37f4feaeca5d561d47da60bfd85db38c3bcb9d33c6f02d9284b61ec844cb1','unreviewed prior fullhash producer')
config=json.loads(small(project/'datasets.json'))
eligible={ds for ds,row in config.items() if not ds.startswith('_') and row.get('use_for_benchmark') is True and 'benchmark_analysis' in row.get('views',{})}
rows=[line.split('\t') for line in small(m/'input_references.tsv').decode().splitlines()]
need(len(rows)==11 and all(len(r)==6 and r[1]=='benchmark_analysis' and re.fullmatch('[a-f0-9]{32}',r[3]) and re.fullmatch('[1-9][0-9]*',r[4]) for r in rows)
     and len({r[0] for r in rows})==11 and {r[0] for r in rows}==eligible,'exact enabled eleven benchmark references required')
need(review['input_references_sha256']==digest(small(m/'input_references.tsv')),'review/input references differ')
selected=datasets.split(',') if mode=='bind' else json.loads(small(m/'scope.json'))['datasets']
need(len(selected)==len(set(selected))==11 and set(selected)==eligible,'explicit exact eleven dataset selection required')
witness={r['dataset']:r for r in w['rows']};need(len(witness)==len(w['rows'])==11 and set(witness)==eligible,'all eleven input witnesses required')
preflight=[line.split('\t') for line in small(m/'preflight.tsv').decode().splitlines()]
need(len(preflight)==11 and all(len(row)==3 and row[1]=='benchmark_analysis' for row in preflight)
     and len({row[0] for row in preflight})==11 and {row[0] for row in preflight}==eligible,'exact original eleven-row physical preflight manifest required')
paths=[]
for ref in w['controls']:
    read_ref(ref)
    if ref['path'].endswith('/datasets.json'):need(ref['sha256']==digest(small(project/'datasets.json')),'prior/current config differs')
for ds,view,path,md5,size,proof in rows:
    entry=config[ds];label=witness[ds]
    need('columns' not in entry['views'][view] and label['configured_sample_column']==entry['columns']['sample']
         and label['label_column']==entry['columns']['label'],'top-level column authority differs')
    need(pathlib.Path(path).name==entry['views'][view]['output_file_name'] and label['view']==view and label['path']==path
         and label['expected_source_md5']==md5 and str(label['size'])==size,'current input identity differs')
    need(label['stat_before']==label['stat_after'] and portable_stat(stat(path),label['stat_after'],path) and stat(path)['st_size']==int(size),'current input stat changed; no H5AD hash attempted')
    need(label['source_identity_order_matches_observed'] is True and label['n_samples']==len(label['sample_ids'])
         and len(set(label['sample_ids']))==label['n_samples'] and (ds!='Stephenson' or label['n_samples']==53), 'current sample identity differs')
    need(label['proof']['path']==proof and fields(read_ref({k:label['proof'][k] for k in ('path','sha256')}))==dict(PATH=path,MD5=md5,SIZE=size),'input checksum proof differs')
    record=label['preflight_record'];status=label['preflight_status']
    record_fields=fields(read_ref({k:record[k] for k in ('path','sha256')}));status_fields=fields(read_ref({k:status[k] for k in ('path','sha256')}))
    parent=pathlib.Path(record['path']).parents[2]
    task=next(i for i,row in enumerate(preflight,1) if row==[ds,view,path])
    need(record['path']==str(parent/'manifests/artifacts'/(digest(path.encode())[:32]+'.record'))
         and record_fields==dict(PATH=path,MD5=md5,SIZE=size,RUN_ID=parent.name,PRODUCER='stage5_preflight',STATE='PUBLISHED')
         and status['path']==str(parent/'status/h5ad_preflight'/(ds+'__'+view+'.status'))
         and review['preflight_manifest']['path']==str(parent/'manifests/h5ad_preflight.tsv')
         and status_fields==dict(STATE='OK',RUN_ID=parent.name,DATASET=ds,VIEW=view,TASK_ID=str(task)),
         'prior fullhash published/native task proof differs')
    # The historical observation stores record stats; proof/status bytes are
    # SHA-bound above and their current stats are frozen in this run's scope.
    need(portable_stat(stat(record['path']),record['stat'],record['path']),'native record stat changed')
    need(stat(path)['st_mtime_ns']<=record['stat']['st_mtime_ns'] and stat(path)['st_ctime_ns']<=record['stat']['st_mtime_ns'],'source postdates fullhash record')
    paths.extend([path,proof,record['path'],status['path']])
selection=''.join(ds+'\tbenchmark_analysis\t'+method+'\n' for analysis,method in [('res50','ECODA_seuratres_50'),('harmony','ECODA_seuratres_2_harmony')] for ds in selected).encode()
outputs=''.join(str(root/'artifacts/saved_labels'/analysis/name)+'\n' for analysis,method in [('res50','ECODA_seuratres_50'),('harmony','ECODA_seuratres_2_harmony')]
                for name in [*(ds+'_'+method+'.rds' for ds in selected),'derived_source_manifest.json','derived_run_manifest.json','status.json','.ecoda_run_owner']).encode()
paths=list(dict.fromkeys(paths));inputs=(''.join(p+'\n' for p in paths)).encode()
if mode=='bind':
    put('operation',b'saved-clusters\n');put('selection.tsv',selection);put('outputs.tsv',outputs);put('inputs.tsv',inputs)
    identity=dict(schema=1,purpose='saved_cluster_counts_screen',run_id=root.name,datasets=selected,
        config_sha256=digest(small(project/'datasets.json')),metadata_sha256={p.name:digest(small(p)) for p in m.iterdir() if p.is_file()},
        input_stats={p:stat(p) for p in paths})
    put('scope.json',(json.dumps(identity,sort_keys=True)+'\n').encode())
else:
    identity=json.loads(small(m/'scope.json'))
    need(identity['schema']==1 and identity['purpose']=='saved_cluster_counts_screen' and identity['run_id']==root.name
         and identity['config_sha256']==digest(small(project/'datasets.json')) and small(m/'selection.tsv')==selection
         and small(m/'outputs.tsv')==outputs and small(m/'inputs.tsv')==inputs,'frozen scope differs')
    for name,sha in identity['metadata_sha256'].items():need(digest(small(m/name))==sha,'run-owned metadata changed')
    need(set(identity['input_stats'])==set(paths),'selected stat inventory differs')
    for p in paths:portable_stat(stat(p),identity['input_stats'][p],p)
for path in list(observed_stats):stat(path)
PY
}
if [[ "${MODE}" == controller ]]; then
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest"
  ecoda_runtime_validate_submission apptainer
  metadata bind
  for NAME in "${SCOPE_FILES[@]}" outputs.tsv; do
    ecoda_write_checksum "${M}/${NAME}"
    ecoda_write_artifact_record "${M}/${NAME}" saved_cluster_scope "${ECODA_RUN_ID}" >/dev/null
  done
fi
validate_scope() {
  local name path
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" ]] || return 1
  ecoda_runtime_validate_bound_run || return 1
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${M}/source.manifest" || return 1
  for name in "${SCOPE_FILES[@]}" outputs.tsv; do
    ecoda_validate_artifact_record "${M}/${name}" saved_cluster_scope "${ECODA_RUN_ID}" >/dev/null || return 1
  done
  while IFS= read -r path; do ecoda_require_input_ownership "${path}" "${ECODA_RUN_ID}" || return 1; done < "${M}/inputs.tsv"
  metadata check
}
validate_analysis() {
  ${PIXI_RSCRIPT} --vanilla - "${PROJECT_ROOT}" "${RUN}" "$1" <<'R'
args<-commandArgs(TRUE);source(file.path(args[1],'src/5_run_benchmark_methods/run_local_ecoda_derived.R'))
analysis<-args[3];stopifnot(analysis%in%c('res50','harmony'))
scope<-jsonlite::fromJSON(file.path(args[2],'manifests/scope.json'),simplifyVector=FALSE)
input_witness<-jsonlite::fromJSON(file.path(args[2],'manifests/input_witness.json'),simplifyVector=FALSE)
ids<-unlist(scope$datasets,use.names=FALSE);root<-file.path(args[2],'artifacts/saved_labels',analysis)
method<-if(analysis=='res50')'ECODA_seuratres_50' else 'ECODA_seuratres_2_harmony'
expected<-c(paste0(ids,'_',method,'.rds'),'derived_source_manifest.json','derived_run_manifest.json','status.json','.ecoda_run_owner')
actual<-list.files(root,all.files=TRUE,no..=TRUE);stopifnot(setequal(actual,c(expected,paste0(expected,'.md5'))))
run<-jsonlite::fromJSON(file.path(root,'derived_run_manifest.json'),simplifyVector=FALSE)
native<-jsonlite::fromJSON(file.path(root,'derived_source_manifest.json'),simplifyVector=FALSE)
for(x in list(run,native))ecoda_validate_derived_manifest(x)
stopifnot(run$status=='COMPLETED',run$analysis==analysis,run$run_id==paste0(basename(args[2]),'_',analysis),
 run$source_mode=='bound_saved_cluster_obs',run$input_references_path==file.path(args[2],'manifests/input_references.tsv'),
 run$input_references_sha256==digest::digest(file=file.path(args[2],'manifests/input_references.tsv'),algo='sha256'),
 run$config_path==file.path(args[1],'datasets.json'),run$config_md5==unname(tools::md5sum(run$config_path)),
 run$runner_path==file.path(args[1],'src/5_run_benchmark_methods/run_local_ecoda_derived.R'),
 run$runner_md5==unname(tools::md5sum(run$runner_path)),
 run$source_manifest_path==file.path(root,'derived_source_manifest.json'),
 run$source_manifest_md5==unname(tools::md5sum(run$source_manifest_path)),
 identical(native$sources,run$sources),native$run_id==run$run_id,native$analysis==analysis,
 identical(vapply(run$sources,`[[`,character(1),'dataset'),ids),length(run$artifacts)==11L,
 .ecoda_run_owner_matches(root,run$run_id))
status<-jsonlite::fromJSON(file.path(root,'status.json'),simplifyVector=FALSE)
stopifnot(status$status=='COMPLETED',status$run_id==run$run_id)
for(i in seq_along(ids)) {
  artifact<-run$artifacts[[i]];path<-file.path(root,paste0(ids[i],'_',method,'.rds'))
  proof<-.ecoda_validate_output_sidecar(path);stopifnot(artifact$path==path,artifact$md5==proof$MD5,artifact$method==method)
  observed<-Filter(function(x)identical(x$dataset,ids[i]),input_witness$rows)
  stopifnot(length(observed)==1L)
  expected_ids<-unlist(observed[[1]]$sample_ids,use.names=FALSE)
  truth<-unlist(observed[[1]]$named_labels,use.names=TRUE)[expected_ids]
  value<-readRDS(path);.ecoda_validate_bundle(value,expected_ids,method)
  stopifnot(all(.ecoda_screen_metrics()%in%names(value$scores)),all(is.finite(unlist(value$scores[.ecoda_screen_metrics()]))),
    value$scores$lisi_score>=0,value$scores$lisi_score<=1,identical(value$scoring_policy,'calc_sep_score_five_v1_lisi_separation'),
    identical(as.character(value$labels),unname(truth)),
    value$source_identity$dataset==ids[i],value$source_identity$analysis==analysis,
    value$source_identity$h5ad_md5==run$sources[[i]]$h5ad_md5,
    run$sources[[i]]$h5ad_md5==observed[[1]]$expected_source_md5,
    run$sources[[i]]$h5ad_path==observed[[1]]$path)
  rm(value);gc(verbose=FALSE)
}
R
}
validate_report() {
  local path owner
  validate_scope || return 1
  [[ "$(sed -n 's/^STATE=//p' "${STATUS}")" == OK && "$(sed -n 's/^RUN_ID=//p' "${STATUS}")" == "${ECODA_RUN_ID}" ]] || return 1
  while IFS= read -r path; do
    ecoda_validate_artifact_record "${path}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
    owner="$(ecoda_artifact_owner_dir "${path}")" || return 1
    _ecoda_artifact_owner_validate_dir "${owner}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 &&
       "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  done < "${M}/outputs.tsv"
  validate_analysis res50 && validate_analysis harmony
}
if [[ "${MODE}" == validate ]]; then
  validate_report || exit 1
  [[ "$(sed -n 's/^STATE=//p' "${RUN}/status/terminal")" == OK && "$(sed -n 's/^RUN_ID=//p' "${RUN}/status/terminal")" == "${ECODA_RUN_ID}" ]] || exit 1
  for NAME in scheduler_ids.tsv sbatch_submission.txt; do ecoda_validate_artifact_record "${M}/${NAME}" saved_cluster_submission "${ECODA_RUN_ID}" >/dev/null; done
  echo "SAVED_CLUSTER_SCREEN_VALIDATED=${ECODA_RUN_ID}"
  exit 0
fi
if [[ "${MODE}" == worker ]]; then
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker
  cd "${PROJECT_ROOT}";validate_scope || exit 1
  [[ ! -e "${OUTPUT}" && ! -L "${OUTPUT}" ]] || exit 1
  export ECODA_DERIVED_ROOT="${OUTPUT}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
  BOUND_CSV="$(${PYTHON_BIN} -c 'import json,sys;print(",".join(json.load(open(sys.argv[1]))["datasets"]))' "${M}/scope.json")"
  for ANALYSIS in res50 harmony; do
    ${PIXI_RSCRIPT} --vanilla "${RUNNER}" --analysis "${ANALYSIS}" --scope benchmark_union --config_path "${DATASETS_JSON_FILE}" \
      --input_references "${M}/input_references.tsv" --input_references_sha256 "$(ecoda_sha256_file "${M}/input_references.tsv")" \
      --datasets "${BOUND_CSV}" --output_dir "${OUTPUT}/${ANALYSIS}" --run_id "${ECODA_RUN_ID}_${ANALYSIS}"
    validate_scope || exit 1
    while IFS= read -r PATH_IN; do
      case "${PATH_IN}" in "${OUTPUT}/${ANALYSIS}/"*)
        case "${PATH_IN}" in *.rds) ecoda_validate_checksum "${PATH_IN}" ;; *) ecoda_write_checksum "${PATH_IN}" ;; esac
        ecoda_write_artifact_record "${PATH_IN}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
        ;; esac
    done < "${M}/outputs.tsv"
    validate_analysis "${ANALYSIS}"
    while IFS= read -r PATH_IN; do
      case "${PATH_IN}" in "${OUTPUT}/${ANALYSIS}/"*) ecoda_artifact_owner_set_state "${PATH_IN}" OK saved_cluster_analysis_complete >/dev/null ;; esac
    done < "${M}/outputs.tsv"
  done
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\n"
  validate_report || exit 1
  trap - EXIT
  exit 0
fi
validate_scope || exit 1
mkdir -p "${OUTPUT%/*}"
while IFS= read -r PATH_IN; do ecoda_artifact_owner_acquire "${PATH_IN}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null; done < "${M}/outputs.tsv"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)"
set +e
JOB="$(sbatch --parsable --wait --partition="${SLURM_PARTITION}" --ntasks=1 --cpus-per-task=1 --mem=64G --time=02:00:00 \
  --output="${RUN}/logs/saved_cluster_%j.log" --error="${RUN}/logs/saved_cluster_%j.err" \
  --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" "${SELF}" --worker)"
RC=$?
set -e
ecoda_atomic_write "${M}/sbatch_submission.txt" "${JOB}\n"
JOB="${JOB%%;*}";[[ "${JOB}" =~ ^[0-9]+$ ]] || { echo 'ERROR: ambiguous saved-cluster submission; preserve run' >&2;exit 1; }
echo "SAVED_CLUSTER_SCREEN_JOB_ID=${JOB}"
ecoda_atomic_write "${M}/scheduler_ids.tsv" "ARRAY\t${JOB}\n"
for NAME in scheduler_ids.tsv sbatch_submission.txt; do ecoda_write_checksum "${M}/${NAME}";ecoda_write_artifact_record "${M}/${NAME}" saved_cluster_submission "${ECODA_RUN_ID}" >/dev/null;done
[[ ${RC} -eq 0 ]] || exit "${RC}"
validate_report || exit 1
ecoda_set_run_state OK saved_cluster_complete
trap - EXIT
