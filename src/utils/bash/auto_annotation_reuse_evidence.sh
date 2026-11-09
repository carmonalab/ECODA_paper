#!/bin/bash
# Finite metadata-only evidence; this entry never runs benchmark methods/sync.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && -n "${ECODA_SOURCE_ROOT:-}" ]]; then
  SCRIPT_DIR="${ECODA_SOURCE_ROOT:?}/src/utils/bash"
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
source "${SCRIPT_DIR}/ecoda_derivative_artifacts.sh"
SELF="${PROJECT_ROOT}/src/utils/bash/auto_annotation_reuse_evidence.sh"
PURPOSE=metadata_only_auto_annotation_reuse_evidence
PRODUCER=auto_annotation_reuse_evidence
MODE=controller
if [[ "${1:-}" == --worker || "${1:-}" == --validate ]]; then MODE="${1#--}"; shift; fi
SELECTION=""; REFERENCES=""; RESULTS=""; RELEASE=""; LEGACY_ARTIFACTS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --selection-file) SELECTION="${2:-}"; shift 2 ;;
    --input-artifacts) REFERENCES="${2:-}"; shift 2 ;;
    --results-root) RESULTS="${2:-}"; shift 2 ;;
    --derivative-release) RELEASE="${2:-}"; shift 2 ;;
    --legacy-artifacts) LEGACY_ARTIFACTS="${2:-}"; [[ -n "${LEGACY_ARTIFACTS}" ]] || exit 1; shift 2 ;;
    *) echo "ERROR: unknown evidence argument: $1" >&2; exit 1 ;;
  esac
done
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && -n "${ECODA_RUN_ID:-}" ]] || {
  echo 'ERROR: metadata evidence requires the immutable snapshot executor/run identity' >&2; exit 1;
}
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
if [[ "${MODE}" == controller ]]; then
  [[ -n "${SELECTION}" && -n "${REFERENCES}" && -n "${RELEASE}" &&
     "${RESULTS}" == "${HPC_SCRATCH_DIR}/benchmark/results" ]] || {
    echo 'ERROR: exact selection, references, canonical legacy results and release are required' >&2; exit 1;
  }
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  ecoda_atomic_write "${ECODA_RUN_ROOT}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
  trap 'ecoda_set_run_state FAIL metadata_evidence_failed || true' EXIT
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
  export ECODA_RUNTIME_PROFILE=stage5
  ecoda_runtime_validate_submission apptainer
  ecoda_atomic_copy "${SELECTION}" "${ECODA_RUN_ROOT}/manifests/selection.tsv"
  ecoda_derivative_bind stage5 "${RELEASE}" "${ECODA_RUN_ROOT}/manifests/selection.tsv" "${REFERENCES}"
else
  [[ -z "${SELECTION}${REFERENCES}${RESULTS}${RELEASE}${LEGACY_ARTIFACTS}" ]] || exit 1
  ecoda_open_run "${ECODA_RUN_ID}" stage5
  export ECODA_ARTIFACT_BINDING="${ECODA_RUN_ROOT}/manifests/artifact_bindings.json"
fi
RUN="${ECODA_RUN_ROOT}"
export ECODA_RUN_ROOT ECODA_RUN_ID
MANIFESTS="${RUN}/manifests"
REPORT="${RUN}/artifacts/evidence/auto_annotation_reuse.json"
STATUS="${RUN}/status/auto_annotation_reuse.status"
export ECODA_ARTIFACT_BINDING

# The controller freezes discovery once; workers only consume this finite list.
if [[ "${MODE}" == controller ]]; then
  ecoda_atomic_copy "${RESULTS%/results}/checksums.md5" "${MANIFESTS}/legacy_checksums.md5"
  if [[ -n "${LEGACY_ARTIFACTS}" ]]; then
    [[ "${LEGACY_ARTIFACTS}" == /* && -f "${LEGACY_ARTIFACTS}" && ! -L "${LEGACY_ARTIFACTS}" ]] || {
      echo 'ERROR: --legacy-artifacts requires an absolute regular TSV' >&2; exit 1;
    }
    ecoda_atomic_copy "${LEGACY_ARTIFACTS}" "${MANIFESTS}/legacy_artifacts.requested.tsv"
  fi
  "${PYTHON_BIN}" - "${RUN}" "${RESULTS}" "${PURPOSE}" "${DATASETS_JSON_FILE}" <<'PY'
import hashlib,json,pathlib,re,sys
run,results=map(pathlib.Path,sys.argv[1:3]);purpose=sys.argv[3];m=run/'manifests'
config=json.loads(pathlib.Path(sys.argv[4]).read_text())
known_auto={'layer1','layer2','layer3','scATOMIC_pred'}|{'layer_'+str(i) for i in range(1,7)}
expected={(d,'benchmark_analysis','auto_annotation_reuse_evidence') for d in ('Bassez','Zhang','Pelka')}
rows=[tuple(x.split('\t')) for x in (m/'selection.tsv').read_text().splitlines()]
if len(rows)!=3 or set(rows)!=expected: raise SystemExit('ERROR: evidence selection must name exactly the three approved dataset/views')
binding=json.loads((m/'artifact_bindings.json').read_text());refs=[x.split('\t') for x in (m/'input_artifacts.tsv').read_text().splitlines()]
if len(refs)!=3 or {(x[0],x[1]) for x in refs}!={x[:2] for x in expected}: raise SystemExit('ERROR: evidence input references must match exact three-cohort selection')
checks={}
for line in (m/'legacy_checksums.md5').read_text().splitlines():
    parts=line.split()
    if len(parts)!=2 or not re.fullmatch('[a-fA-F0-9]{32}',parts[0]): raise SystemExit('ERROR: malformed legacy checksum table')
    name=pathlib.Path(parts[1].lstrip('*')).name
    if name in checks: raise SystemExit('ERROR: duplicate checksum basename')
    checks[name]=parts[0].lower()
request=m/'legacy_artifacts.requested.tsv'
explicit=None
if request.exists():
    explicit=[line.split('\t') for line in request.read_text().splitlines()]
    if not explicit or any(len(row)!=2 or row[0] not in {'Bassez','Zhang','Pelka'} or not pathlib.Path(row[1]).is_absolute() for row in explicit): raise SystemExit('ERROR: malformed explicit legacy artifact scope')
    if len({tuple(row) for row in explicit})!=len(explicit) or len({row[1] for row in explicit})!=len(explicit) or {row[0] for row in explicit}!={'Bassez','Zhang','Pelka'}: raise SystemExit('ERROR: duplicate or incomplete explicit legacy artifact scope')
inventory=[]
for dataset in ('Bassez','Zhang','Pelka'):
    columns=config[dataset]['columns']
    author_keys=[key for key,field in [('ECODA_authors_LR','cell_type_low_res'),('ECODA_authors_HR','cell_type_high_res')] if columns.get(field) in known_auto]
    alternatives=['composition','ECODA_HiTME_.*','ECODA_scATOMIC.*']+[re.escape(key) for key in author_keys]
    pattern=re.compile(re.escape(dataset)+r'_('+ '|'.join(alternatives)+r')\.rds$')
    files=sorted(pathlib.Path(row[1]) for row in explicit if row[0]==dataset) if explicit is not None else sorted(p for p in results.iterdir() if pattern.fullmatch(p.name))
    if not files: raise SystemExit('ERROR: no auto-composition evidence artifacts for '+dataset)
    for p in files:
        if p.parent!=results or not pattern.fullmatch(p.name) or not re.fullmatch(r'[A-Za-z0-9_.-]+',p.name) or p.is_symlink() or not p.is_file() or p.resolve()!=p or p.stat().st_size<=0 or p.name not in checks: raise SystemExit('ERROR: unsafe/unrecorded selected RDS: '+str(p))
        if explicit is not None:
            expected=checks[p.name]; digest=hashlib.md5()
            before=p.stat()
            with p.open('rb') as stream:
                for block in iter(lambda:stream.read(8*1024*1024),b''): digest.update(block)
            actual=digest.hexdigest();after=p.stat()
            if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns) or actual!=expected: raise SystemExit(f'ERROR: selected legacy checksum mismatch: {p}; expected={expected}; actual={actual}')
        inventory.append(f'{dataset}\t{p}\t{checks[p.name]}\t{p.stat().st_size}\n')
(m/'legacy_artifacts.tsv').write_text(''.join(inventory))
files=['source.manifest','runtime.identity','selection.tsv','input_artifacts.tsv','artifact_bindings.json','legacy_checksums.md5','legacy_artifacts.tsv']
if explicit is not None: files.append('legacy_artifacts.requested.tsv')
paths=[pathlib.Path(r[2]) for r in refs]+[pathlib.Path(r.split('\t')[1]) for r in inventory]
identity={'schema':1,'purpose':purpose,'run_id':run.name,'results_root':str(results),'legacy_inventory_mode':'explicit' if explicit is not None else 'discovery','input_stats':{str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in paths},'identities':{n:hashlib.sha256((m/n).read_bytes()).hexdigest() for n in files}}
(m/'evidence_identity.json').write_text(json.dumps(identity,sort_keys=True)+'\n')
PY
  for name in source.manifest runtime.identity selection.tsv input_artifacts.tsv artifact_bindings.json legacy_checksums.md5 legacy_artifacts.tsv evidence_identity.json; do
    ecoda_write_checksum "${MANIFESTS}/${name}"
    ecoda_write_artifact_record "${MANIFESTS}/${name}" auto_annotation_reuse_scope "${ECODA_RUN_ID}" >/dev/null
  done
  if [[ -n "${LEGACY_ARTIFACTS}" ]]; then
    ecoda_write_checksum "${MANIFESTS}/legacy_artifacts.requested.tsv"
    ecoda_write_artifact_record "${MANIFESTS}/legacy_artifacts.requested.tsv" auto_annotation_reuse_scope "${ECODA_RUN_ID}" >/dev/null
  fi
  mkdir -p "${REPORT%/*}"
  ecoda_artifact_owner_acquire "${REPORT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
fi

validate_scope() {
  ecoda_derivative_load stage5 || return 1
  local dataset view path md5 size proof owner parent record producer
  while IFS=$'\t' read -r dataset view path md5 size proof; do
    ecoda_require_input_ownership "${path}" "${ECODA_RUN_ID}" || return 1
    case "${path}" in "${ECODA_RUNS_ROOT}/"*)
      owner="$(ecoda_artifact_owner_dir "${path}")" || return 1
      _ecoda_artifact_owner_validate_dir "${owner}" "${path}" || return 1
      [[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
      parent="${path#${ECODA_RUNS_ROOT}/}"; parent="${parent%%/*}"
      [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${parent}" ]] || return 1
      record="$(ecoda_artifact_record_path "${path}" "${parent}")" || return 1
      producer="$(sed -n 's/^PRODUCER=//p' "${record}")" || return 1
      ecoda_validate_artifact_record "${path}" "${producer}" "${parent}" "${proof}" >/dev/null || return 1
      ;;
    esac
  done < "${MANIFESTS}/input_artifacts.tsv"
  while IFS=$'\t' read -r dataset path md5 size; do
    ecoda_require_input_ownership "${path}" "${ECODA_RUN_ID}" || return 1
  done < "${MANIFESTS}/legacy_artifacts.tsv"
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" ]] || return 1
  ecoda_runtime_validate_bound_run || return 1
  ecoda_validate_checksum "${MANIFESTS}/evidence_identity.json" || return 1
  ecoda_validate_artifact_record "${MANIFESTS}/evidence_identity.json" auto_annotation_reuse_scope "${ECODA_RUN_ID}" >/dev/null || return 1
  "${PYTHON_BIN}" - "${RUN}" "${PURPOSE}" <<'PY'
import hashlib,json,pathlib,sys
run=pathlib.Path(sys.argv[1]);m=run/'manifests';identity=json.loads((m/'evidence_identity.json').read_text())
expected={'source.manifest','runtime.identity','selection.tsv','input_artifacts.tsv','artifact_bindings.json','legacy_checksums.md5','legacy_artifacts.tsv'}
if identity.get('legacy_inventory_mode','discovery')=='explicit': expected.add('legacy_artifacts.requested.tsv')
elif identity.get('legacy_inventory_mode','discovery')!='discovery': raise SystemExit('ERROR: unknown legacy inventory mode')
if identity['schema']!=1 or identity['purpose']!=sys.argv[2] or identity['run_id']!=run.name or set(identity['identities'])!=expected: raise SystemExit('ERROR: evidence identity scope differs')
for filename,recorded in identity['input_stats'].items():
    p=pathlib.Path(filename)
    if p.is_symlink() or [p.stat().st_size,p.stat().st_mtime_ns]!=recorded: raise SystemExit('ERROR: evidence input changed since inventory freeze: '+filename)
for name,sha in identity['identities'].items():
    p=m/name
    if p.is_symlink() or hashlib.sha256(p.read_bytes()).hexdigest()!=sha: raise SystemExit('ERROR: frozen evidence identity changed: '+name)
PY
}
validate_report() {
  validate_scope || return 1
  [[ -s "${STATUS}" && "$(sed -n 's/^STATE=//p' "${STATUS}")" == OK &&
     "$(sed -n 's/^RUN_ID=//p' "${STATUS}")" == "${ECODA_RUN_ID}" &&
     "$(sed -n 's/^PURPOSE=//p' "${STATUS}")" == "${PURPOSE}" ]] || return 1
  ecoda_validate_checksum "${REPORT}" || return 1
  ecoda_validate_artifact_record "${REPORT}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
  owner="$(ecoda_artifact_owner_dir "${REPORT}")" || return 1
  _ecoda_artifact_owner_validate_dir "${owner}" "${REPORT}" || return 1
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  "${PYTHON_BIN}" - "${RUN}" "${REPORT}" "${PURPOSE}" "${DATASETS_JSON_FILE}" <<'PY'
import json,pathlib,sys
run,report=map(pathlib.Path,sys.argv[1:3]);x=json.loads(report.read_text());scope={'Bassez','Zhang','Pelka'}
config=json.loads(pathlib.Path(sys.argv[4]).read_text())
known_auto={'layer1','layer2','layer3','scATOMIC_pred'}|{'layer_'+str(i) for i in range(1,7)}
if x['schema']!=1 or x['purpose']!=sys.argv[3] or x['run_id']!=run.name or set(x['scope'])!=scope or set(x['cohorts'])!=scope: raise SystemExit('ERROR: evidence report identity differs')
inventory=[r.split('\t') for r in (run/'manifests/legacy_artifacts.tsv').read_text().splitlines()]
refs={r[0]:r for r in [r.split('\t') for r in (run/'manifests/input_artifacts.tsv').read_text().splitlines()]}
for ds in scope:
    c=x['cohorts'][ds];expected={r[1]:r[2] for r in inventory if r[0]==ds}
    actual=c['artifact_inventory']
    if isinstance(actual,dict): actual=[actual]
    actual={r['V2']:r['V3'] for r in actual}
    if actual!=expected or c['input']!=refs[ds][2] or c['expected_h5ad_md5']!=refs[ds][3] or not c['obs_sha256'] or c['n_cells']<=0: raise SystemExit('ERROR: report input/artifact coverage differs')
    if set(c['artifact_keys'])!=set(expected): raise SystemExit('ERROR: report omitted an artifact key inventory')
    expected_pairs={(p,key) for p,keys in c['artifact_keys'].items() for key in keys}
    pairs=[(r['artifact'],r['method']) for r in c['rows']]
    if len(pairs)!=len(set(pairs)) or set(pairs)!=expected_pairs: raise SystemExit('ERROR: report omitted or duplicated automatic rows')
    fields={f:config[ds]['columns'].get(f) for f in ['cell_type_low_res','cell_type_high_res']}
    overlap={k:v for k,v in fields.items() if v in known_auto}
    if c['configured_ct_fields']!=fields or (c['configured_ct_autofield_overlap'] or {})!=overlap: raise SystemExit('ERROR: configured annotation overlap differs')
    for row in c['rows']:
        dependency=row['method'] in {'ECODA_authors_LR','ECODA_authors_HR'}
        if dependency and fields['cell_type_low_res' if row['method'].endswith('_LR') else 'cell_type_high_res'] not in known_auto: raise SystemExit('ERROR: unselected author annotation scope')
        if row['evidence_role']!=('dependency_evidence' if dependency else 'automatic_composition_evidence'): raise SystemExit('ERROR: evidence role differs')
        if row['artifact'] not in expected or row['artifact_md5']!=expected[row['artifact']] or row['status'] not in {'BLOCKED','INVALID_CURRENT_COUNTS','COUNT_MATCH_PROVENANCE_PENDING'}: raise SystemExit('ERROR: invalid evidence row')
    for field,status in [('invalid_methods','INVALID_CURRENT_COUNTS'),('blocked_methods','BLOCKED')]:
        declared=c[field] if isinstance(c[field],list) else [c[field]]
        if set(declared)!={r['method'] for r in c['rows'] if r['status']==status}: raise SystemExit('ERROR: report summary differs from classified rows')
PY
}
if [[ "${MODE}" == validate ]]; then validate_report || exit 1; echo "AUTO_ANNOTATION_REUSE_EVIDENCE_VALIDATED=${ECODA_RUN_ID}"; exit 0; fi
if [[ "${MODE}" == worker ]]; then
  failed() {
    ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\nPURPOSE=${PURPOSE}\n" || true
    ecoda_artifact_owner_set_state "${REPORT}" FAIL evidence_worker_failed >/dev/null || true
  }
  trap failed EXIT
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker
  cd "${PROJECT_ROOT}"
  validate_scope || exit 1
  temporary="${REPORT}.tmp.${SLURM_JOB_ID:-$$}"
  [[ ! -e "${REPORT}" && ! -L "${REPORT}" && ! -e "${temporary}" && ! -L "${temporary}" ]]
  export RETICULATE_PYTHON="${PYTHON_BIN}"
  ${PIXI_RSCRIPT} --vanilla "${PROJECT_ROOT}/src/utils/r/audit_auto_annotation_reuse.R" \
    "${PROJECT_ROOT}" "${HPC_SCRATCH_DIR}" "${HPC_SCRATCH_DIR}/benchmark/results" \
    "${MANIFESTS}/input_artifacts.tsv" "${MANIFESTS}/legacy_artifacts.tsv" \
    "${MANIFESTS}/legacy_checksums.md5" "${ECODA_OWNERS_ROOT}" "${ECODA_RUN_ID}" "${temporary}"
  validate_scope || exit 1
  ln "${temporary}" "${REPORT}"
  rm -f "${temporary}"
  ecoda_write_checksum "${REPORT}"
  ecoda_write_artifact_record "${REPORT}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
  ecoda_artifact_owner_set_state "${REPORT}" OK evidence_worker_complete >/dev/null
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\nPURPOSE=${PURPOSE}\n"
  validate_report || exit 1
  trap - EXIT
  exit 0
fi
controller_failed() {
  ecoda_set_run_state FAIL metadata_evidence_failed || true
  ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\nPURPOSE=${PURPOSE}\n" || true
  local owner
  owner="$(ecoda_artifact_owner_dir "${REPORT}" 2>/dev/null)" || return 0
  if _ecoda_artifact_owner_validate_dir "${owner}" "${REPORT}" &&
     [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]; then
    ecoda_artifact_owner_set_state "${REPORT}" FAIL evidence_controller_failed >/dev/null || true
  fi
}
trap controller_failed EXIT
validate_scope || exit 1
runtime_export="$(ecoda_runtime_export_csv stage5 0)"
set +e
job="$(sbatch --parsable --wait --partition="${SLURM_PARTITION}" --ntasks=1 --cpus-per-task=1 \
  --mem=64G --time=01:00:00 --output="${RUN}/logs/evidence_%j.log" --error="${RUN}/logs/evidence_%j.err" \
  --export="ALL,${runtime_export},ECODA_RUN_ROOT=${RUN},ECODA_ARTIFACT_BINDING=${ECODA_ARTIFACT_BINDING}" "${SELF}" --worker)"
rc=$?
set -e
job="${job%%;*}"
[[ "${job}" =~ ^[0-9]+$ ]] || exit 1
printf 'AUTO_ANNOTATION_REUSE_EVIDENCE_JOB_ID=%s\n' "${job}"
ecoda_atomic_write "${MANIFESTS}/scheduler_id" "${job}\n"
[[ ${rc} -eq 0 ]] || exit "${rc}"
validate_report || exit 1
ecoda_set_run_state OK metadata_evidence_complete
trap - EXIT
