#!/bin/bash
# Real R projection and common ownership/record helpers; runtime is a stub.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-analysis-projection.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
TREE="${TMP}/snapshot/tree"; SCRATCH="${TMP}/scratch"
mkdir -p "${TREE}/src/utils/bash" "${TREE}/src/notebooks" "${TREE}/src/5_run_benchmark_methods" "${SCRATCH}/_ecoda_inputs" "${TMP}/bin"
cp "${ROOT}/src/utils/bash/project_analysis_release.sh" "${ROOT}/src/utils/bash/ecoda_run_common.sh" "${TREE}/src/utils/bash/"
cp "${ROOT}/src/utils/analysis_release.R" "${ROOT}/src/utils/evaluation_scores.R" "${TREE}/src/utils/"
cp "${ROOT}/src/5_run_benchmark_methods/benchmark_pipeline.R" "${TREE}/src/5_run_benchmark_methods/"
for name in datasets_io.R math_utils.R scoring_metrics.R plotting.R; do cp "${ROOT}/src/utils/${name}" "${TREE}/src/utils/"; done
cp "${ROOT}/src/notebooks/benchmark_analysis.R" "${TREE}/src/notebooks/"
cp "${ROOT}/datasets.json" "${TREE}/datasets.json"
printf 'SOURCE_ROOT=%s\n' "${TREE}" > "${TREE}/source.manifest"
cat > "${TREE}/src/slurm_config.sh" <<'SH'
export PROJECT_ROOT="${TEST_TREE}" HPC_SCRATCH_DIR="${TEST_SCRATCH}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners"
export PYTHON_BIN="${TEST_PYTHON}" PIXI_RSCRIPT="${TEST_RSCRIPT}"
SLURM_PARTITION=shared-cpu
SH
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'SH'
ecoda_runtime_validate_submission() {
  [[ "$1" == apptainer && "${ECODA_RUNTIME_PROFILE}" == stage5 ]]
  printf 'fixture pinned runtime\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"
}
ecoda_runtime_validate_bound_run() {
  [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]]
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
}
ecoda_runtime_export_csv() { printf ECODA_RUNTIME_PROFILE=stage5; }
ecoda_runtime_reexec_worker() {
  [[ "$1" == stage5 && "$2" == "${TEST_TREE}/src/utils/bash/project_analysis_release.sh" && "$3" == --worker ]]
  printf 'PINNED_RUNTIME_WORKER\n'
  [[ "${PROJECTION_TEST_FAIL:-0}" != 1 ]]
}
SH
cat > "${TMP}/bin/sbatch" <<'SH'
#!/bin/bash
set -euo pipefail
[[ "$*" == *'--parsable --wait'* && "$*" == *'--cpus-per-task=1'* && "$*" == *'--mem=64G'* && "$*" == *'--time=02:00:00'* && "$*" != *'--gpus'* ]]
printf '%s\n' "$*" >> "${TEST_SCRATCH}/sbatch.calls"
args=("$@")
set +e
SLURM_JOB_ID=765432 bash "${args[${#args[@]}-2]}" "${args[${#args[@]}-1]}" \
  > "${ECODA_RUN_ROOT}/logs/analysis_projection_765432.log" 2> "${ECODA_RUN_ROOT}/logs/analysis_projection_765432.err"
rc=$?
set -e
printf '765432;fixture\n'
exit "${rc}"
SH
chmod +x "${TMP}/bin/sbatch"
for command in sacct squeue; do
  printf '#!/bin/bash\nprintf forbidden >> "%s"\nexit 99\n' "${TMP}/scheduler.called" > "${TMP}/bin/${command}"
  chmod +x "${TMP}/bin/${command}"
done
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}"
export TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_RSCRIPT="${ROOT}/.pixi/envs/default/bin/Rscript"
export ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export PATH="${TMP}/bin:${PATH}"
CONTROL="${SCRATCH}/_ecoda_inputs/control.rds"
"${TEST_RSCRIPT}" --vanilla - "${ROOT}" "${SCRATCH}/_ecoda_inputs" "${CONTROL}" <<'R'
a<-commandArgs(TRUE);source(file.path(a[1],"src/utils/analysis_release.R"));root<-a[2]
ids<-paste0("s",1:6);x<-matrix(1:120,6,dimnames=list(ids,paste0("RNA_gene",1:20)))
labels<-setNames(rep(c("a","b"),each=3),ids)
path<-file.path(root,"RNA.rds");saveRDS(list(feat_mat=x,counts=x,RNA=x,labels=labels,dist_mat=dist(x),scores=list(anosim_score=.5)),path)
evidence<-file.path(root,"review.txt");writeLines("Approved synthetic legacy result",evidence)
md5<-function(p)unname(tools::md5sum(p))
ledger<-data.frame(dataset="D",view="benchmark_analysis",method_key="Pseudobulk_hvg2000",kind="standalone_result",
 artifact_path=path,bundle_key="",run_root=root,md5=md5(path),receipt_path=file.path(root,"receipt.rds"),receipt_md5="")
receipt<-c(as.list(ledger[c("dataset","view","method_key","kind","artifact_path","bundle_key","run_root","md5")]),
 list(schema="ecoda_analysis_release_receipt_v1",status="reviewed",origin="legacy_reuse",input_md5=strrep("a",32),source_policy="fixture",
 evidence="reviewed",evidence_path=evidence,evidence_md5=md5(evidence)))
saveRDS(receipt,ledger$receipt_path);ledger$receipt_md5<-md5(ledger$receipt_path)
scope<-ledger[c("dataset","view","method_key")];scope$require_new<-FALSE;scope$current_input_md5<-strrep("a",32);scope$source_policy<-"fixture";scope$status<-"approved"
saveRDS(list(ledger=ledger,expected_scope=scope),a[3])
R
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
MD5="$(ecoda_md5_file "${CONTROL}")"; SOURCE_BEFORE="$(ecoda_md5_file "${SCRATCH}/_ecoda_inputs/RNA.rds")"
invoke() {
  ECODA_RUN_ID="$1" ECODA_RUN_ROOT="${SCRATCH}/_ecoda_runs/$1" \
    bash "${TREE}/src/utils/bash/project_analysis_release.sh" --analysis-release "${CONTROL}" --analysis-release-md5 "${MD5}"
}
reject() { if "$@" > "${TMP}/failure.log" 2>&1; then echo 'expected projection failure' >&2; exit 1; fi; }
invoke good > "${TMP}/success.log"
run="${SCRATCH}/_ecoda_runs/good"; output="${run}/artifacts/analysis_projection"
grep -q '^STATE=OK$' "${run}/status/terminal"
grep -q '^PURPOSE=analysis_metadata_projection$' "${run}/metadata"
[[ "$(cat "${run}/manifests/selection.tsv")" == $'D\tbenchmark_analysis\tPseudobulk_hvg2000' ]]
ecoda_validate_artifact_record "${run}/manifests/selection.tsv" analysis_projection_scope good >/dev/null
grep -q '^PINNED_RUNTIME_WORKER$' "${run}/logs/analysis_projection_765432.log"
[[ "$(cat "${run}/manifests/scheduler_ids.tsv")" == $'ARRAY\t765432' ]]
ecoda_validate_artifact_record "${run}/manifests/scheduler_ids.tsv" analysis_projection_submission good >/dev/null
while IFS=$'\t' read -r path hash; do
  [[ "$(ecoda_md5_file "${path}")" == "${hash}" ]]
  ecoda_validate_artifact_record "${path}" analysis_metadata_projection good >/dev/null
  _ecoda_artifact_owner_validate_dir "$(ecoda_artifact_owner_dir "${path}")" "${path}"
  [[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]]
done < "${output}/files.tsv"
"${TEST_RSCRIPT}" --vanilla - "${TREE}" "${output}/analysis_release.rds" <<'R'
a<-commandArgs(TRUE);source(file.path(a[1],"src/utils/analysis_release.R"));x<-readRDS(a[2]);r<-load_analysis_release(x$ledger,x$expected_scope,require_projection=TRUE)[[1]]
stopifnot(is.null(r$result$feat_mat),is.null(r$result$counts),is.null(r$result$RNA),r$result$n_features==20L,r$receipt$origin=="legacy_reuse")
R
reject invoke good
[[ "$(ecoda_md5_file "${SCRATCH}/_ecoda_inputs/RNA.rds")" == "${SOURCE_BEFORE}" ]]
ECODA_RUN_ID=good ECODA_RUN_ROOT="${run}" bash "${TREE}/src/utils/bash/project_analysis_release.sh" --validate > "${TMP}/validate.log"
grep -q '^ANALYSIS_PROJECTION_VALIDATED=good$' "${TMP}/validate.log"
# Read-only terminal validation rejects failed/foreign inventory ownership
# without changing completed status, terminal evidence or any owner bytes.
inventory_owner="$(ecoda_artifact_owner_dir "${output}/files.tsv")/owner"
cp "${inventory_owner}" "${TMP}/inventory.owner.saved"
sed 's/STATE=OK/STATE=FAIL/' "${TMP}/inventory.owner.saved" > "${inventory_owner}"
status_hash="$(ecoda_md5_file "${run}/status/analysis_projection.status")"
terminal_hash="$(ecoda_md5_file "${run}/status/terminal")"
owner_hash="$(ecoda_md5_file "${inventory_owner}")"
directory_owner="$(ecoda_artifact_owner_dir "${output}")/owner"
directory_hash="$(ecoda_md5_file "${directory_owner}")"
ECODA_RUN_ID=good ECODA_RUN_ROOT="${run}" reject bash "${TREE}/src/utils/bash/project_analysis_release.sh" --validate
[[ "$(ecoda_md5_file "${run}/status/analysis_projection.status")" == "${status_hash}" &&
   "$(ecoda_md5_file "${run}/status/terminal")" == "${terminal_hash}" &&
   "$(ecoda_md5_file "${inventory_owner}")" == "${owner_hash}" &&
   "$(ecoda_md5_file "${directory_owner}")" == "${directory_hash}" ]]
cp "${TMP}/inventory.owner.saved" "${inventory_owner}"
cp "${run}/status/terminal" "${TMP}/terminal.saved"
for terminal_fault in failed wrong_run; do
  if [[ "${terminal_fault}" == failed ]]; then
    sed 's/STATE=OK/STATE=FAIL/' "${TMP}/terminal.saved" > "${run}/status/terminal"
  else
    sed 's/RUN_ID=good/RUN_ID=foreign/' "${TMP}/terminal.saved" > "${run}/status/terminal"
  fi
  terminal_hash="$(ecoda_md5_file "${run}/status/terminal")"
  ECODA_RUN_ID=good ECODA_RUN_ROOT="${run}" reject bash "${TREE}/src/utils/bash/project_analysis_release.sh" --validate
  [[ "$(ecoda_md5_file "${run}/status/analysis_projection.status")" == "${status_hash}" &&
     "$(ecoda_md5_file "${run}/status/terminal")" == "${terminal_hash}" &&
     "$(ecoda_md5_file "${inventory_owner}")" == "$(ecoda_md5_file "${TMP}/inventory.owner.saved")" &&
     "$(ecoda_md5_file "${directory_owner}")" == "${directory_hash}" ]]
done
cp "${TMP}/terminal.saved" "${run}/status/terminal"
PROJECTION_TEST_FAIL=1 reject invoke failed_worker
failed_run="${SCRATCH}/_ecoda_runs/failed_worker"
grep -q '^STATE=FAIL$' "${failed_run}/status/terminal"
[[ "$(cat "${failed_run}/manifests/scheduler_ids.tsv")" == $'ARRAY\t765432' ]]
ecoda_validate_artifact_record "${failed_run}/manifests/scheduler_ids.tsv" analysis_projection_submission failed_worker >/dev/null
[[ ! -e "${failed_run}/artifacts/analysis_projection" ]]
PCA_RENDER_FIXTURE_ONLY=1 PCA_RENDER_FIXTURE_ROOT="${SCRATCH}/_ecoda_inputs/render_data" \
  "${TEST_RSCRIPT}" --vanilla "${ROOT}/tests/test_benchmark_pca_render.R" > "${TMP}/render-fixture.log" 2>&1
RENDER_CONTROL="${SCRATCH}/_ecoda_inputs/render_data/control.rds"
ECODA_RUN_ID=render ECODA_RUN_ROOT="${SCRATCH}/_ecoda_runs/render" \
  bash "${TREE}/src/utils/bash/project_analysis_release.sh" --analysis-release "${RENDER_CONTROL}" \
    --analysis-release-md5 "$(ecoda_md5_file "${RENDER_CONTROL}")" --operation render-pca > "${TMP}/render.log"
render="${SCRATCH}/_ecoda_runs/render"
grep -q '^STATE=OK$' "${render}/status/terminal"
grep -q '^PURPOSE=analysis_pca_render$' "${render}/metadata"
[[ "$(wc -l < "${render}/manifests/selection.tsv" | tr -d ' ')" == 34 ]]
[[ "$(cat "${render}/manifests/operation")" == render-pca ]]
ecoda_validate_artifact_record "${render}/artifacts/benchmark_pca/render_manifest.rds" analysis_pca_render render >/dev/null
[[ "$(wc -l < "${render}/artifacts/benchmark_pca/files.tsv" | tr -d ' ')" == 20 ]]
ECODA_RUN_ID=render ECODA_RUN_ROOT="${render}" bash "${TREE}/src/utils/bash/project_analysis_release.sh" --validate > "${TMP}/render-validation.log"
ecoda_artifact_owner_acquire "${SCRATCH}/_ecoda_inputs/RNA.rds" stage5 foreign_writer 0 0 0 >/dev/null
reject invoke blocked_writer
grep -q '^STATE=FAIL$' "${SCRATCH}/_ecoda_runs/blocked_writer/status/terminal"
[[ ! -e "${SCRATCH}/_ecoda_runs/blocked_writer/artifacts/analysis_projection" ]]
[[ "$(ecoda_md5_file "${SCRATCH}/_ecoda_inputs/RNA.rds")" == "${SOURCE_BEFORE}" && ! -e "${TMP}/scheduler.called" ]]
[[ "$(wc -l < "${SCRATCH}/sbatch.calls" | tr -d ' ')" == 3 ]]
echo 'PASS: one pinned CPU projection worker, ID-bound logs/records including failed worker, real projection/owners/origin, input-writer guard and read-only terminal validation'
