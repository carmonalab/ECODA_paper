#!/bin/bash
# Actual new entry/common/binding helpers; only scheduler/runtime/R execution stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-auto-evidence.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
TREE="${TMP}/snapshot/0000000000000000000000000000000000000002/tree"
SCRATCH="${TMP}/scratch"
mkdir -p "${TREE}/src/utils/bash" "${TREE}/src/utils/py" "${TREE}/src/utils/r" "${SCRATCH}/benchmark/results" "${TMP}/nas" "${TMP}/bin"
for f in ecoda_run_common.sh ecoda_derivative_artifacts.sh auto_annotation_reuse_evidence.sh; do cp "${ROOT}/src/utils/bash/${f}" "${TREE}/src/utils/bash/${f}"; done
cp "${ROOT}/src/utils/py/derivative_artifacts.py" "${TREE}/src/utils/py/derivative_artifacts.py"
printf '# fixture identity\n' > "${TREE}/source.manifest"
printf '# actual R is covered independently by test_auto_annotation_reuse.R\n' > "${TREE}/src/utils/r/audit_auto_annotation_reuse.R"
cat > "${TREE}/src/slurm_config.sh" <<'EOF'
export PROJECT_ROOT="${TEST_TREE}" DATASETS_JSON_FILE="${TEST_TREE}/datasets.json"
export HPC_SCRATCH_DIR="${TEST_SCRATCH}" NAS_TARGET_DIR="${TEST_NAS}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners"
export PYTHON_BIN="${TEST_PYTHON}" PIXI_RSCRIPT="${TEST_RSCRIPT}"
SLURM_PARTITION=shared-cpu
EOF
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'EOF'
ecoda_runtime_validate_submission() { printf 'immutable fixture runtime\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run() { [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]] && cmp -s "${ECODA_RUN_ROOT}/manifests/source.manifest" "${ECODA_SOURCE_MANIFEST}"; }
ecoda_runtime_export_csv() { printf ECODA_RUNTIME_PROFILE=stage5; }
ecoda_runtime_reexec_worker() { [[ "$1" == stage5 && "$3" == --worker ]]; ecoda_runtime_validate_bound_run; }
EOF
cat > "${TMP}/bin/Rscript" <<'EOF'
#!/bin/bash
set -euo pipefail
[[ "$1" == --vanilla && "$2" == */src/utils/r/audit_auto_annotation_reuse.R ]]
[[ "${EVIDENCE_FIXTURE_FAIL:-}" != R ]] || exit 7
shift 2
"${TEST_PYTHON}" - "$@" <<'PY'
import hashlib,json,pathlib,sys
source,scratch,results,refs,inventory,legacy,owners,run=sys.argv[1:]
rows=[x.split('\t') for x in pathlib.Path(inventory).read_text().splitlines()]
config=json.loads((pathlib.Path(source)/'datasets.json').read_text())
references={x[0]:x for x in (x.split('\t') for x in pathlib.Path(refs).read_text().splitlines())}
x={'schema':1,'purpose':'metadata_only_auto_annotation_reuse_evidence','run_id':run,'scope':['Bassez','Zhang','Pelka'],'cohorts':{}}
for d in x['scope']:
    selected=[r for r in rows if r[0]==d]
    # Deliberately mutate the live directory after freeze: unselected files stay out.
    pathlib.Path(results,d+'_ECODA_HiTME_NEW.rds').write_text('unselected after freeze')
    keys={r[1]:['ECODA_HiTME_HR_layer2' if pathlib.Path(r[1]).stem==d+'_composition' else pathlib.Path(r[1]).stem[len(d)+1:]] for r in selected}
    fields={f:config[d]['columns'].get(f) for f in ['cell_type_low_res','cell_type_high_res']}
    overlap={k:v for k,v in fields.items() if v in {'layer1','layer2','layer3','scATOMIC_pred'}|{'layer_'+str(i) for i in range(1,7)}}
    x['cohorts'][d]={'configured_ct_fields':fields,'configured_ct_autofield_overlap':overlap,'input':references[d][2],'expected_h5ad_md5':references[d][3],'obs_sha256':'f'*64,'n_cells':4,'artifact_inventory':[dict(zip(['V1','V2','V3','V4'],r)) for r in selected],'artifact_keys':keys,'invalid_methods':[k[0] for k in keys.values()],'blocked_methods':[],'rows':[{'method':keys[r[1]][0],'evidence_role':'dependency_evidence' if keys[r[1]][0].startswith('ECODA_authors_') else 'automatic_composition_evidence','artifact':r[1],'artifact_md5':r[2],'status':'INVALID_CURRENT_COUNTS'} for r in selected]}
print(json.dumps(x))
PY
EOF
cat > "${TMP}/bin/sbatch" <<'EOF'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${TEST_SBATCH_CAPTURE}"
[[ "$*" == *'--mem=64G'* && "$*" == *'--cpus-per-task=1'* && "$*" == *'--time=01:00:00'* && "$*" != *'--array'* ]]
if [[ "${EVIDENCE_FIXTURE_FAIL:-}" == OOM ]]; then printf '812345\n'; exit 137; fi
args=("$@")
set +e
mv "${TEST_NAS}" "${TEST_NAS}.hidden"
bash "${args[${#args[@]}-2]}" "${args[${#args[@]}-1]}" >&2
rc=$?
mv "${TEST_NAS}.hidden" "${TEST_NAS}"
set -e
printf '812345\n'
exit "${rc}"
EOF
chmod +x "${TMP}/bin/Rscript" "${TMP}/bin/sbatch"
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_NAS="${TMP}/nas"
export TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_RSCRIPT="${TMP}/bin/Rscript" TEST_SBATCH_CAPTURE="${TMP}/sbatch.calls"
export ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export PATH="${TMP}/bin:${PATH}" ECODA_RUNTIME_MODE=host
"${TEST_PYTHON}" - "${TMP}" "${TREE}" "${SCRATCH}" <<'PY'
import hashlib,json,pathlib,sys
base,tree,scratch=map(pathlib.Path,sys.argv[1:]);refs=[];checks=[];config={}
for d in ('Bassez','Zhang','Pelka'):
    config[d]={'columns':{'sample':'Sample','label':'biology'},'views':{'benchmark_analysis':{'output_file_name':d+'.h5ad'}}}
    if d=='Zhang':config[d]['columns'].update(cell_type_low_res='layer1',cell_type_high_res='layer2')
    else:config[d]['columns'].update(cell_type_low_res='author_lr',cell_type_high_res='author_hr')
    p=scratch/d/'output'/f'{d}.h5ad';p.parent.mkdir(parents=True);p.write_text('small metadata-bound input fixture')
    md5=hashlib.md5(p.read_bytes()).hexdigest();proof=pathlib.Path(str(p)+'.md5');proof.write_text(f'MD5={md5}\nSIZE={p.stat().st_size}\nPATH={p}\n')
    refs.append(f'{d}\tbenchmark_analysis\t{p}\t{md5}\t{p.stat().st_size}\t{proof}\n')
    r=scratch/'benchmark/results'/f'{d}_composition.rds';r.write_text('RDS deserialize stub; real RDS covered by R test')
    checks.append(hashlib.md5(r.read_bytes()).hexdigest()+'  results/'+r.name+'\n')
    for key in ('ECODA_authors_LR','ECODA_authors_HR'):
        r=scratch/'benchmark/results'/f'{d}_{key}.rds';r.write_text('compact dependency evidence fixture')
        checks.append(hashlib.md5(r.read_bytes()).hexdigest()+'  results/'+r.name+'\n')
(tree/'datasets.json').write_text(json.dumps(config));(base/'refs.tsv').write_text(''.join(refs))
(base/'selection.tsv').write_text(''.join(d+'\tbenchmark_analysis\tauto_annotation_reuse_evidence\n' for d in config))
(scratch/'benchmark/checksums.md5').write_text(''.join(checks))
PY
entry="${TREE}/src/utils/bash/auto_annotation_reuse_evidence.sh"
# The repaired Bassez route has a real prior-run owner and producer receipt.
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
ecoda_init_run stage2 bassez_producer >/dev/null
parent="${ECODA_RUN_ROOT}/artifacts/datasets/Bassez/output/Bassez.h5ad"
mkdir -p "${parent%/*}"
cp "${SCRATCH}/Bassez/output/Bassez.h5ad" "${parent}"
ecoda_write_checksum "${parent}"
ecoda_write_artifact_record "${parent}" bassez_h5ad_metadata bassez_producer >/dev/null
ecoda_artifact_owner_acquire "${parent}" stage2 bassez_producer 0 0 0 >/dev/null
ecoda_artifact_owner_set_state "${parent}" OK fixture >/dev/null
"${TEST_PYTHON}" - "${TMP}/refs.tsv" "${parent}" <<'PY_PARENT'
import hashlib,pathlib,sys
f,p=map(pathlib.Path,sys.argv[1:]);rows=[r.split('\t') for r in f.read_text().splitlines()]
for r in rows:
    if r[0]=='Bassez':r[2:]=[str(p),hashlib.md5(p.read_bytes()).hexdigest(),str(p.stat().st_size),str(p)+'.md5']
f.write_text(''.join('\t'.join(r)+'\n' for r in rows))
PY_PARENT
controller() { ECODA_RUN_ID="$1" bash "${entry}" --selection-file "${TMP}/selection.tsv" --input-artifacts "${TMP}/refs.tsv" --results-root "${SCRATCH}/benchmark/results" --derivative-release test-release; }
# Preserve every scientific input/receipt/result byte and mode across the audit.
"${TEST_PYTHON}" - "${SCRATCH}" "${TMP}/before.json" "${parent}" <<'PY'
import hashlib,json,pathlib,stat,sys
root=pathlib.Path(sys.argv[1]);selected=list(root.glob('*/output/*'))+list((root/'benchmark/results').glob('*'))+[root/'benchmark/checksums.md5',pathlib.Path(sys.argv[3]),pathlib.Path(sys.argv[3]+'.md5')]
pathlib.Path(sys.argv[2]).write_text(json.dumps({str(p):[hashlib.sha256(p.read_bytes()).hexdigest(),stat.S_IMODE(p.stat().st_mode)] for p in selected}))
PY
controller good > "${TMP}/good.log" 2>&1
ECODA_RUN_ID=good bash "${entry}" --validate
[[ "$(wc -l < "${TEST_SBATCH_CAPTURE}" | tr -d ' ')" == 1 ]]
grep -q '^AUTO_ANNOTATION_REUSE_EVIDENCE_JOB_ID=812345$' "${TMP}/good.log"
[[ -e "${SCRATCH}/_ecoda_runs/good/artifacts/evidence/auto_annotation_reuse.json" ]]
[[ ! -e "${TMP}/nas/releases/test-release/good" ]]
[[ "$(wc -l < "${SCRATCH}/_ecoda_runs/good/manifests/legacy_artifacts.tsv" | tr -d ' ')" == 5 ]]
! grep -q 'Bassez_ECODA_authors_\|Pelka_ECODA_authors_' "${SCRATCH}/_ecoda_runs/good/manifests/legacy_artifacts.tsv"
"${TEST_PYTHON}" - "${TMP}/before.json" <<'PY'
import hashlib,json,pathlib,stat,sys
for p,(sha,mode) in json.loads(pathlib.Path(sys.argv[1]).read_text()).items():
    f=pathlib.Path(p);assert hashlib.sha256(f.read_bytes()).hexdigest()==sha and stat.S_IMODE(f.stat().st_mode)==mode
PY
# New unbound result files cannot broaden a fresh or frozen audit silently.
if controller unrecorded >"${TMP}/unrecorded.log" 2>&1; then echo 'unrecorded RDS accepted' >&2; exit 1; fi
for ds in Bassez Zhang Pelka; do rm "${SCRATCH}/benchmark/results/${ds}_ECODA_HiTME_NEW.rds"; done
# A published derivative input cannot be consumed with active/failed/missing producer proof.
for state in ACTIVE FAIL; do
  ecoda_artifact_owner_set_state "${parent}" "${state}" fixture >/dev/null
  if controller "bad_owner_${state}" >"${TMP}/bad_owner_${state}.log" 2>&1; then cat "${TMP}/bad_owner_${state}.log" >&2; echo 'invalid input owner accepted' >&2; exit 1; fi
done
ecoda_artifact_owner_set_state "${parent}" OK fixture >/dev/null
record="$(ecoda_artifact_record_path "${parent}" bassez_producer)"
mv "${record}" "${record}.saved"
if controller missing_producer >"${TMP}/missing_producer.log" 2>&1; then echo 'missing producer receipt accepted' >&2; exit 1; fi
mv "${record}.saved" "${record}"
for failure in R OOM; do
  if EVIDENCE_FIXTURE_FAIL="${failure}" controller "failed_${failure}" >"${TMP}/failure_${failure}.log" 2>&1; then echo 'failed evidence job accepted' >&2; exit 1; fi
  grep -q 'STATE=FAIL' "${SCRATCH}/_ecoda_runs/failed_${failure}/status/terminal"
  grep -q 'STATE=FAIL' "${SCRATCH}/_ecoda_runs/failed_${failure}/status/auto_annotation_reuse.status"
done
# Frozen scope/report tampering is rejected without scheduler access.
report="${SCRATCH}/_ecoda_runs/good/artifacts/evidence/auto_annotation_reuse.json"
chmod u+w "${report}"
printf X >> "${report}"
if ECODA_RUN_ID=good bash "${entry}" --validate >"${TMP}/tamper.log" 2>&1; then echo 'corrupt report accepted' >&2; exit 1; fi
[[ "$(wc -l < "${TEST_SBATCH_CAPTURE}" | tr -d ' ')" == 3 ]]
printf 'PASS: one finite metadata evidence job, frozen artifact inventory, actual checksum/owner contract, strict failures, no NAS/scientific writes\n'
