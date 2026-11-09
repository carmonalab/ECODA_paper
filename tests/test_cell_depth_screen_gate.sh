#!/bin/bash
# Real compact RDS preflight/common ownership; only runtime, scoring and Slurm stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-cell-depth-gate.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
cleanup() {
  local RC=$?
  if [[ ${RC} -ne 0 ]]; then
    for LOG in success.log validated.log partial.log failure.log; do
      if [[ -f "${TMP}/${LOG}" ]]; then echo "${LOG}:" >&2; tail -30 "${TMP}/${LOG}" >&2; fi
    done
  fi
  chmod -R u+w "${TMP}"; rm -rf "${TMP}"
}
trap cleanup EXIT
TREE="${TMP}/snapshot/0000000000000000000000000000000000000003/tree"
SCRATCH="${TMP}/scratch"
mkdir -p "${TREE}/src/utils/bash" "${TREE}/src/5_run_benchmark_methods" "${TMP}/bin" "${TMP}/nas" "${SCRATCH}/inputs"
for NAME in ecoda_run_common.sh cell_depth_screen.sh; do cp "${ROOT}/src/utils/bash/${NAME}" "${TREE}/src/utils/bash/${NAME}"; done
cp "${ROOT}/src/utils/analysis_release.R" "${TREE}/src/utils/analysis_release.R"
cp "${ROOT}/src/5_run_benchmark_methods/run_local_ecoda_derived.R" "${TREE}/src/5_run_benchmark_methods/run_local_ecoda_derived.R"
cp "${ROOT}/datasets.json" "${TREE}/datasets.json"
printf 'immutable test source\n' > "${TREE}/source.manifest"
cat > "${TREE}/src/slurm_config.sh" <<'SH'
export PROJECT_ROOT="${TEST_TREE}" DATASETS_JSON_FILE="${TEST_TREE}/datasets.json"
export HPC_SCRATCH_DIR="${TEST_SCRATCH}" NAS_TARGET_DIR="${TEST_NAS}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners"
export PYTHON_BIN="${TEST_PYTHON}" PIXI_RSCRIPT="${TEST_RSCRIPT}"
SLURM_PARTITION=shared-cpu
SH
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'SH'
ecoda_runtime_validate_submission() { [[ "$1" == apptainer ]]; printf 'immutable runtime fixture\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run() { [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]] && cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"; }
ecoda_runtime_export_csv() { [[ "$1" == stage5 && "$2" == 0 ]]; printf 'ECODA_RUNTIME_PROFILE=stage5'; }
ecoda_runtime_reexec_worker() { [[ "$1" == stage5 && "$3" == --worker && "$2" == "${TEST_TREE}/src/utils/bash/cell_depth_screen.sh" ]]; ecoda_runtime_validate_bound_run; }
SH
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_NAS="${TMP}/nas"
export TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_REAL_R="${ROOT}/.pixi/envs/default/bin/Rscript"
export TEST_RSCRIPT="${TMP}/bin/Rscript" TEST_SBATCH_CAPTURE="${TMP}/sbatch.calls" TEST_R_CAPTURE="${TMP}/R.calls"
export ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_MODE=host PATH="${TMP}/bin:${PATH}"
"${TEST_REAL_R}" --vanilla - "${SCRATCH}/inputs" <<'R'
root <- commandArgs(TRUE)[1]; md5 <- function(p) unname(tools::md5sum(p))
evidence <- file.path(root,"review.txt"); writeLines("reviewed current compact fixture",evidence)
ledger <- scope <- list(); contracts <- list(); k <- 0L
for (dataset in c("Bassez","Stephenson")) {
  n <- if(dataset=="Bassez")26L else 53L; types <- if(dataset=="Bassez")40L else 4L
  total <- if(dataset=="Bassez")74948 else 10600
  ids <- paste0(dataset,seq_len(n)); counts <- matrix(2,n,types,dimnames=list(ids,paste0("CT",seq_len(types))))
  counts[1,1] <- counts[1,1]+total-sum(counts)
  labels <- setNames(factor(rep(c("a","b"),length.out=n)),ids)
  hr <- file.path(root,paste0(dataset,"_HR.rds")); meta <- file.path(root,paste0(dataset,"_metadata.rds"))
  saveRDS(list(counts=counts+.5,labels=labels),hr)
  saveRDS(list(labels=labels,cells_per_sample=setNames(rowSums(counts),ids),n_samples=n,n_cells=total,n_cell_types_high_res=types),meta)
  contracts[[dataset]] <- list(sample_ids=ids,n_samples=n,n_cell_types=types,total_cells=total,annotated_cells=total)
  for (method in c("ECODA_authors_HR","__dataset_metadata__")) {
    k <- k+1L; path <- if(method=="ECODA_authors_HR")hr else meta
    row <- data.frame(dataset=dataset,view="benchmark_analysis",method_key=method,kind="standalone_result",
      artifact_path=path,bundle_key="",run_root=root,md5=md5(path),receipt_path=file.path(root,paste0(k,".receipt.rds")),receipt_md5="")
    receipt <- c(as.list(row[c("dataset","view","method_key","kind","artifact_path","bundle_key","run_root","md5")]),
      list(schema="ecoda_analysis_release_receipt_v1",status="reviewed",origin="new",input_md5=strrep("a",32),
        source_policy="current-counts-v1",evidence="reviewed fixture",evidence_path=evidence,evidence_md5=md5(evidence)))
    saveRDS(receipt,row$receipt_path); row$receipt_md5 <- md5(row$receipt_path)
    expected <- row[c("dataset","view","method_key")]; expected$require_new <- TRUE
    expected$current_input_md5 <- strrep("a",32);expected$source_policy <- "current-counts-v1";expected$status <- "approved"
    ledger[[k]] <- row;scope[[k]] <- expected
  }
}
saveRDS(list(schema="ecoda_cell_depth_inputs_v1",status="approved",ledger=do.call(rbind,ledger),
  expected_scope=do.call(rbind,scope),cohort_contract=contracts),file.path(root,"control.rds"))
R
cat > "${TMP}/bin/Rscript" <<'SH'
#!/bin/bash
set -euo pipefail
[[ "$1" == --vanilla ]]; shift
if [[ "$1" == - ]]; then exec "${TEST_REAL_R}" --vanilla "$@"; fi
[[ "$1" == "${TEST_TREE}/src/5_run_benchmark_methods/run_local_ecoda_derived.R" ]]; shift
printf '%s\n' "$*" >> "${TEST_R_CAPTURE}"
"${TEST_PYTHON}" - "$@" <<'PY'
import hashlib,json,os,pathlib,sys
args=sys.argv[1:];a=dict(zip(args[::2],args[1::2]))
assert len(args)%2==0 and len(a)*2==len(args)
assert a['--scope']=='benchmark_union' and a['--datasets']==os.environ.get('TEST_SELECTED','Bassez,Stephenson')
assert '--composition_dir' not in a and '--h5ad' not in a and '--historical_inputs' not in a
assert a['--config_path']==os.environ['TEST_TREE']+'/datasets.json'
control=pathlib.Path(a['--analysis_release']);assert control.parent.name=='manifests'
assert hashlib.md5(control.read_bytes()).hexdigest()==a['--analysis_release_md5']
if os.environ.get('TEST_ON_COMPUTE')=='1':
    assert not pathlib.Path(os.environ['TEST_NAS']).exists(), 'worker must be testable without NAS'
assert os.environ['OMP_NUM_THREADS']==os.environ['OPENBLAS_NUM_THREADS']==os.environ['NUMBA_NUM_THREADS']=='1'
fail=os.environ.get('CELL_DEPTH_FIXTURE_FAIL','')
if '--analysis' in a:
    assert a['--analysis']=='cell_subsetting' and a['--cores']=='1'
    p=pathlib.Path(a['--output_dir']);assert not p.exists();p.mkdir()
    (p/'partial.txt').write_text('preserved only on failure')
    if fail=='R':raise SystemExit(7)
    (p/'partial.txt').unlink()
    for name in ['ECODA_authors_HR_cell_subsetting.rds','ECODA_authors_HR_cell_subsetting.csv','derived_source_manifest.json','derived_run_manifest.json','status.json','.ecoda_run_owner']:
        f=p/name;f.write_text('fixture '+name+' '+a['--run_id']+'\n')
        if name.endswith(('.rds','.csv')):
            pathlib.Path(str(f)+'.md5').write_text(f'MD5={hashlib.md5(f.read_bytes()).hexdigest()}\nSIZE={f.stat().st_size}\nPATH={f}\n')
else:
    assert set(a)=={'--validate_output','--scope','--config_path','--analysis_release','--analysis_release_md5','--datasets'}
    if fail=='VALIDATE':raise SystemExit(8)
    p=pathlib.Path(a['--validate_output'])
    assert (p/'ECODA_authors_HR_cell_subsetting.rds').read_text().startswith('fixture ')
PY
SH
cat > "${TMP}/bin/sbatch" <<'SH'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${TEST_SBATCH_CAPTURE}"
[[ "$*" == *'--parsable --wait'* && "$*" == *'--cpus-per-task=1'* && "$*" == *'--mem=64G'* && "$*" == *'--time=02:00:00'* && "$*" != *'--array'* ]]
if [[ "${CELL_DEPTH_FIXTURE_FAIL:-}" == OOM ]]; then printf '812345\n'; exit 137; fi
if [[ "${CELL_DEPTH_FIXTURE_FAIL:-}" == TAMPER ]]; then printf 'changed after bind\n' >> "${TEST_SCRATCH}/inputs/review.txt"; fi
if [[ "${CELL_DEPTH_FIXTURE_FAIL:-}" == REUSE ]]; then mkdir -p "${ECODA_RUN_ROOT}/artifacts/cell_depth/analysis"; printf keep > "${ECODA_RUN_ROOT}/artifacts/cell_depth/analysis/keep"; fi
ARGS=("$@")
mv "${TEST_NAS}" "${TEST_NAS}.hidden"
set +e
TEST_ON_COMPUTE=1 bash "${ARGS[${#ARGS[@]}-2]}" "${ARGS[${#ARGS[@]}-1]}" >&2
RC=$?
set -e
mv "${TEST_NAS}.hidden" "${TEST_NAS}"
printf '812345;fixture\n'
exit "${RC}"
SH
for NAME in squeue sacct rsync; do
  printf '#!/bin/bash\nprintf unexpected >> "${TEST_SCRATCH}/forbidden.calls"\nexit 99\n' > "${TMP}/bin/${NAME}"
done
chmod +x "${TMP}/bin/"*
ENTRY="${TREE}/src/utils/bash/cell_depth_screen.sh"
CONTROL="${SCRATCH}/inputs/control.rds"
MD5="$("${TEST_PYTHON}" -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "${CONTROL}")"
run_screen() { bash "${ENTRY}" --analysis_release "${CONTROL}" --analysis_release_md5 "${MD5}" --datasets "${TEST_SELECTED:-Bassez,Stephenson}"; }
expect_fail() { if "$@" > "${TMP}/failure.log" 2>&1; then echo 'ERROR: expected failure' >&2; cat "${TMP}/failure.log" >&2; exit 1; fi; }
export ECODA_RUN_ID=screen_ok
if ! run_screen > "${TMP}/success.log" 2>&1; then cat "${TMP}/success.log" >&2; exit 1; fi
grep -q '^CELL_DEPTH_SCREEN_JOB_ID=812345$' "${TMP}/success.log"
RUN="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
[[ "$(cat "${RUN}/manifests/scheduler_ids.tsv")" == "$(printf 'worker\t812345')" ]]
# Validation is read-only; allow NAS to remain hidden during terminal verification too.
mv "${TEST_NAS}" "${TEST_NAS}.hidden"
bash "${ENTRY}" --validate > "${TMP}/validated.log" 2>&1
mv "${TEST_NAS}.hidden" "${TEST_NAS}"
grep -q '^CELL_DEPTH_SCREEN_VALIDATED=screen_ok$' "${TMP}/validated.log"
[[ "$(wc -l < "${TEST_SBATCH_CAPTURE}" | tr -d ' ')" == 1 ]]
[[ "$(find "${RUN}/artifacts/cell_depth/analysis" -type f | wc -l | tr -d ' ')" == 12 ]]
expect_fail run_screen
[[ "$(wc -l < "${TEST_SBATCH_CAPTURE}" | tr -d ' ')" == 1 ]]
# Every output gets a published producer record/owner, not only the score pair.
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
for NAME in ECODA_authors_HR_cell_subsetting.rds ECODA_authors_HR_cell_subsetting.csv derived_source_manifest.json derived_run_manifest.json status.json .ecoda_run_owner; do
  P="${RUN}/artifacts/cell_depth/analysis/${NAME}"
  ecoda_validate_artifact_record "${P}" cell_depth_screen screen_ok >/dev/null
  OWNER="$(ecoda_artifact_owner_dir "${P}")"; _ecoda_artifact_owner_validate_dir "${OWNER}" "${P}"
  [[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK && "${ECODA_ARTIFACT_OWNER_RUN}" == screen_ok ]]
done
# Tamper after success: terminal validation rejects both scope and payload changes.
cp "${SCRATCH}/inputs/review.txt" "${TMP}/review.original"
printf tampered >> "${SCRATCH}/inputs/review.txt"
expect_fail bash "${ENTRY}" --validate
cp "${TMP}/review.original" "${SCRATCH}/inputs/review.txt"
P="${RUN}/artifacts/cell_depth/analysis/ECODA_authors_HR_cell_subsetting.rds"
cp "${P}" "${TMP}/output.original"; chmod u+w "${P}"; printf tampered >> "${P}"
expect_fail bash "${ENTRY}" --validate
cp "${TMP}/output.original" "${P}"; chmod a-w "${P}"
# Scheduler identities cannot be rewritten together to claim another job.
cp "${RUN}/manifests/scheduler_ids.tsv" "${TMP}/scheduler.original"
chmod u+w "${RUN}/manifests/scheduler_ids.tsv"
printf 'worker\t999999\n' > "${RUN}/manifests/scheduler_ids.tsv"
expect_fail bash "${ENTRY}" --validate
cp "${TMP}/scheduler.original" "${RUN}/manifests/scheduler_ids.tsv"; chmod a-w "${RUN}/manifests/scheduler_ids.tsv"
# Exact partial recovery stays one selected dataset; no implicit cohort expansion.
export ECODA_RUN_ID=screen_partial TEST_SELECTED=Stephenson
run_screen > "${TMP}/partial.log" 2>&1
unset TEST_SELECTED
# Controller failures must occur before sbatch.
COUNT="$(wc -l < "${TEST_SBATCH_CAPTURE}" | tr -d ' ')"
export ECODA_RUN_ID=screen_bad_csv TEST_SELECTED=Bassez,Bassez
expect_fail run_screen
unset TEST_SELECTED
export ECODA_RUN_ID=screen_bad_md5
expect_fail bash "${ENTRY}" --analysis_release "${CONTROL}" --analysis_release_md5 00000000000000000000000000000000 --datasets Bassez
export ECODA_RUN_ID=screen_unpinned ECODA_SOURCE_SNAPSHOT_REQUIRED=0
expect_fail run_screen
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
# A valid-looking control outside scratch is not compute-visible authority.
cp "${CONTROL}" "${TMP}/outside.rds"
export ECODA_RUN_ID=screen_outside
expect_fail bash "${ENTRY}" --analysis_release "${TMP}/outside.rds" --analysis_release_md5 "${MD5}" --datasets Bassez
# A run-owned source requires its actual published producer, not just a receipt.
mkdir -p "${SCRATCH}/_ecoda_runs/unpublished"
cp "${CONTROL}" "${SCRATCH}/_ecoda_runs/unpublished/control.rds"
export ECODA_RUN_ID=screen_unpublished
expect_fail bash "${ENTRY}" --analysis_release "${SCRATCH}/_ecoda_runs/unpublished/control.rds" --analysis_release_md5 "${MD5}" --datasets Bassez
# A conflicting active source owner is a prelaunch failure, even with valid bytes.
ecoda_artifact_owner_acquire "${SCRATCH}/inputs/review.txt" stage5 another_run 0 0 0 >/dev/null
export ECODA_RUN_ID=screen_active_input
expect_fail run_screen
ecoda_artifact_owner_set_state "${SCRATCH}/inputs/review.txt" OK test_release >/dev/null
[[ "$(wc -l < "${TEST_SBATCH_CAPTURE}" | tr -d ' ')" == "${COUNT}" ]]
# Runtime worker failures preserve artifacts/evidence and emit exactly one known ID.
for FAILURE in R VALIDATE OOM TAMPER REUSE; do
  export ECODA_RUN_ID="screen_fail_${FAILURE}" CELL_DEPTH_FIXTURE_FAIL="${FAILURE}"
  expect_fail run_screen
  F="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
  grep -q '^STATE=FAIL$' "${F}/status/terminal"
  [[ "$(cat "${F}/manifests/scheduler_ids.tsv")" == "$(printf 'worker\t812345')" ]]
  for NAME in ECODA_authors_HR_cell_subsetting.rds ECODA_authors_HR_cell_subsetting.csv derived_source_manifest.json derived_run_manifest.json status.json .ecoda_run_owner; do
    P="${F}/artifacts/cell_depth/analysis/${NAME}"
    OWNER="$(ecoda_artifact_owner_dir "${P}")"; _ecoda_artifact_owner_validate_dir "${OWNER}" "${P}"
    [[ "${ECODA_ARTIFACT_OWNER_STATE}" == FAIL && "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" ]]
  done
  if [[ "${FAILURE}" == R ]]; then [[ -s "${F}/artifacts/cell_depth/analysis/partial.txt" ]]; fi
  if [[ "${FAILURE}" == REUSE ]]; then [[ "$(cat "${F}/artifacts/cell_depth/analysis/keep")" == keep ]]; fi
  if [[ "${FAILURE}" == TAMPER ]]; then cp "${TMP}/review.original" "${SCRATCH}/inputs/review.txt"; fi
done
unset CELL_DEPTH_FIXTURE_FAIL
[[ ! -e "${SCRATCH}/forbidden.calls" ]]
echo 'cell-depth screen gate tests passed'
