#!/bin/bash
# Real compact-parent binding/common proofs; only scoring, runtime and Slurm stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-zeroimp-gate.XXXXXX")";TMP="$(cd "${TMP}" && pwd -P)"
cleanup() { local RC=$?; if [[ ${RC} -ne 0 ]]; then for L in success failure; do [[ ! -f "${TMP}/${L}.log" ]] || tail -35 "${TMP}/${L}.log" >&2; done; fi; chmod -R u+w "${TMP}"; rm -rf "${TMP}"; }
trap cleanup EXIT
TREE="${TMP}/snapshot/0000000000000000000000000000000000000004/tree";SCRATCH="${TMP}/scratch"
mkdir -p "${TREE}/src/utils/bash" "${TREE}/src/notebooks" "${TREE}/src/5_run_benchmark_methods" "${SCRATCH}/inputs" "${TMP}/bin" "${TMP}/nas"
for NAME in ecoda_run_common.sh zeroimp_count_screen.sh; do cp "${ROOT}/src/utils/bash/${NAME}" "${TREE}/src/utils/bash/${NAME}"; done
cp "${ROOT}/src/notebooks/run_zeroimp_count_screen.R" "${TREE}/src/notebooks/"
cp "${ROOT}/src/utils/analysis_release.R" "${TREE}/src/utils/"
for NAME in run_local_ecoda_derived.R benchmark_pipeline.R; do cp "${ROOT}/src/5_run_benchmark_methods/${NAME}" "${TREE}/src/5_run_benchmark_methods/"; done
cp "${ROOT}/datasets.json" "${TREE}/datasets.json"
printf 'source fixture\n' > "${TREE}/source.manifest"
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_NAS="${TMP}/nas" TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_REAL_R="${ROOT}/.pixi/envs/default/bin/Rscript"
export TEST_RSCRIPT="${TMP}/bin/Rscript" TEST_CALLS="${TMP}/sbatch.calls" ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export PATH="${TMP}/bin:${PATH}"
ZEROIMP_TEST_ROOT="${SCRATCH}/inputs" ZEROIMP_TEST_BUILD_ONLY=1 "${TEST_REAL_R}" --vanilla "${ROOT}/tests/test_zeroimp_count_screen.R" > "${TMP}/fixture.log" 2>&1
cat > "${TREE}/src/slurm_config.sh" <<'STUB'
export PROJECT_ROOT="${TEST_TREE}" DATASETS_JSON_FILE="${TEST_TREE}/datasets.json" HPC_SCRATCH_DIR="${TEST_SCRATCH}" NAS_TARGET_DIR="${TEST_NAS}"
export ECODA_RUNS_ROOT="${TEST_SCRATCH}/_ecoda_runs" ECODA_OWNERS_ROOT="${TEST_SCRATCH}/_ecoda_owners" PYTHON_BIN="${TEST_PYTHON}" PIXI_RSCRIPT="${TEST_RSCRIPT}"
SLURM_PARTITION=shared-cpu
STUB
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'STUB'
ecoda_runtime_validate_submission() { [[ "$1" == apptainer ]]; printf 'runtime fixture\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"; }
ecoda_runtime_validate_bound_run() { [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]] && cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"; }
ecoda_runtime_export_csv() { printf 'ECODA_RUNTIME_PROFILE=stage5'; }
ecoda_runtime_reexec_worker() { [[ "$1" == stage5 && "$2" == "${TEST_TREE}/src/utils/bash/zeroimp_count_screen.sh" && "$3" == --worker ]]; ecoda_runtime_validate_bound_run; }
STUB
cat > "${TMP}/bin/Rscript" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ "$1" == --vanilla ]];shift
[[ "$1" == "${TEST_TREE}/src/notebooks/run_zeroimp_count_screen.R" ]]
if [[ "$*" == *' --bind '* ]]; then exec "${TEST_REAL_R}" --vanilla "$@"; fi
"${TEST_PYTHON}" - "$@" <<'PY'
import hashlib,os,pathlib,sys
args=sys.argv[2:];a=dict(zip(args[::2],args[1::2]));assert len(a)*2==len(args)
assert a['--project_root']==os.environ['TEST_TREE'] and a['--config_path']==os.environ['TEST_TREE']+'/datasets.json'
assert a['--datasets']==os.environ.get('TEST_DATASETS','Adams') and a['--run_id']==os.environ['ECODA_RUN_ID']
assert not any(k in a for k in ['--h5ad','--input_root','--force','--historical_inputs'])
p=pathlib.Path(a.get('--produce',a.get('--validate')))
if '--produce' in a:
 assert not p.exists();p.mkdir();(p/'partial').write_text('preserve failure')
 if os.environ.get('TEST_FAIL')=='R':raise SystemExit(7)
 (p/'partial').unlink()
 for n in ['screen_manifest.json','status.json']+[ds+suffix for ds in a['--datasets'].split(',') for suffix in ['_zeroimp_counts.rds','_zeroimp_counts.csv']]:
  f=p/n;f.write_text('fixture '+n+'\n')
  if f.suffix in ['.rds','.csv']:pathlib.Path(str(f)+'.md5').write_text(f'MD5={hashlib.md5(f.read_bytes()).hexdigest()}\nSIZE={f.stat().st_size}\nPATH={f}\n')
else:
 assert (p/'Adams_zeroimp_counts.rds').exists()
 if os.environ.get('TEST_FAIL')=='VALIDATE':raise SystemExit(8)
PY
STUB
cat > "${TMP}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${TEST_CALLS}"
[[ "$*" == *'--parsable --wait'* && "$*" == *'--cpus-per-task=1'* && "$*" == *'--mem=64G'* && "$*" == *'--time=02:00:00'* && "$*" != *'--array'* ]]
if [[ "${TEST_FAIL:-}" == OOM ]]; then printf '712345\n';exit 137;fi
if [[ "${TEST_FAIL:-}" == TAMPER ]]; then printf changed >> "${TEST_SCRATCH}/inputs/review.txt";fi
ARGS=("$@");cp "${TEST_TREE}/src/utils/bash/zeroimp_count_screen.sh" "${TEST_SCRATCH}/spool-entry.sh"
set +e
bash "${TEST_SCRATCH}/spool-entry.sh" "${ARGS[${#ARGS[@]}-1]}" >&2
RC=$?
set -e
printf '712345;fixture\n';exit "${RC}"
STUB
for NAME in squeue sacct rsync; do printf '#!/bin/bash\nprintf forbidden >> "${TEST_SCRATCH}/forbidden.calls"\nexit 99\n' > "${TMP}/bin/${NAME}"; done
chmod +x "${TMP}/bin/"*
ENTRY="${TREE}/src/utils/bash/zeroimp_count_screen.sh";CONTROL="${SCRATCH}/inputs/release.rds"
MD5="$("${TEST_PYTHON}" -c 'import hashlib,sys;print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest())' "${CONTROL}")"
run_screen() { bash "${ENTRY}" --analysis_release "${CONTROL}" --analysis_release_md5 "${MD5}" --datasets "${TEST_DATASETS:-Adams}"; }
expect_fail() { if "$@" > "${TMP}/failure.log" 2>&1; then echo 'expected failure' >&2;exit 1;fi; }
export ECODA_RUN_ID=zeroimp_ok TEST_DATASETS=Adams,Bassez,Gongsharma_cmv_young_males,Kfoury,Kim,Lee,Pelka,Smillie,Stephenson,Wu,Zhang
run_screen > "${TMP}/success.log" 2>&1
RUN="${SCRATCH}/_ecoda_runs/zeroimp_ok"
grep -q '^ZEROIMP_COUNT_SCREEN_JOB_ID=712345$' "${TMP}/success.log"
[[ "$(cat "${RUN}/manifests/operation")" == zeroimp-count-screen && "$(wc -l < "${RUN}/manifests/selection.tsv" | tr -d ' ')" == 11 && "$(cut -f1 "${RUN}/manifests/selection.tsv" | paste -sd, -)" == "${TEST_DATASETS}" ]]
bash "${ENTRY}" --validate > "${TMP}/validate.log" 2>&1
source "${TREE}/src/slurm_config.sh";source "${TREE}/src/utils/bash/ecoda_run_common.sh"
for N in operation selection.tsv; do ecoda_validate_artifact_record "${RUN}/manifests/${N}" zeroimp_scope zeroimp_ok >/dev/null; done
for N in screen_manifest.json status.json $(cut -f1 "${RUN}/manifests/selection.tsv" | sed 's/$/_zeroimp_counts.rds/') $(cut -f1 "${RUN}/manifests/selection.tsv" | sed 's/$/_zeroimp_counts.csv/'); do
 P="${RUN}/artifacts/zeroimp_count_screen/analysis/${N}";ecoda_validate_artifact_record "${P}" zeroimp_count_screen zeroimp_ok >/dev/null
 OWNER="$(ecoda_artifact_owner_dir "${P}")";_ecoda_artifact_owner_validate_dir "${OWNER}" "${P}";[[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]]
done
# Wrong-stage output ownership is foreign even with the same RUN_ID/OK.
P="${RUN}/artifacts/zeroimp_count_screen/analysis/Adams_zeroimp_counts.rds"
OWNER_FILE="$(ecoda_artifact_owner_dir "${P}")/owner"
cp "${OWNER_FILE}" "${TMP}/owner.original"
sed 's/^STAGE=stage5$/STAGE=stage3/' "${OWNER_FILE}" > "${TMP}/owner.wrong-stage"
cat "${TMP}/owner.wrong-stage" > "${OWNER_FILE}"
OWNER_SHA="$(ecoda_sha256_file "${OWNER_FILE}")"
expect_fail bash "${ENTRY}" --validate
[[ "$(ecoda_sha256_file "${OWNER_FILE}")" == "${OWNER_SHA}" ]]
cat "${TMP}/owner.original" > "${OWNER_FILE}"
# Actual generic purpose route and scheduler accounting parser, consuming a
# checksummed fixture capture (never querying scheduler).
eval "$(sed -n '/^_audit_analysis_release() {/,/^}/p' "${ROOT}/src/utils/bash/ecoda_run_audit.sh")"
eval "$(sed -n '/^_audit_scheduler_ids() {/,/^}/p' "${ROOT}/src/utils/bash/ecoda_run_audit.sh")"
_audit_metadata_value() { sed -n "s/^$2=//p" "$1"; }
_audit_metadata_count() { awk -F= -v k="$2" '$1==k{n++}END{print n+0}' "$1"; }
_audit_regular_file() { [[ -f "$1" && ! -L "$1" ]]; }
_audit_die() { echo "$*" >&2;return 1; }
export ECODA_RUN_ID=zeroimp_ok ECODA_RUN_ROOT="${RUN}"
RUN_ROOT_REAL="${RUN}" RUN_ID=zeroimp_ok SELECTION_ARG="${RUN}/manifests/selection.tsv" AUDIT_TERMINAL_STATE=OK
AUDIT_SOURCE_ROOT="${TREE}" SOURCE_MANIFEST_ARG="${TREE}/source.manifest"
SCHEDULER_ACCOUNTING_ARG="${RUN}/manifests/scheduler_accounting.tsv"
printf '712345|COMPLETED|0:0\n' > "${SCHEDULER_ACCOUNTING_ARG}"
ecoda_write_checksum "${SCHEDULER_ACCOUNTING_ARG}"
_audit_analysis_release > "${TMP}/generic.log"
grep -q '^ZEROIMP_COUNT_SCREEN_VALIDATED=zeroimp_ok$' "${TMP}/generic.log"
cp "${RUN}/manifests/operation" "${TMP}/operation.original"
chmod u+w "${RUN}/manifests/operation";printf 'projection\n' > "${RUN}/manifests/operation"
expect_fail _audit_analysis_release
cp "${TMP}/operation.original" "${RUN}/manifests/operation";chmod a-w "${RUN}/manifests/operation"
expect_fail run_screen
[[ "$(wc -l < "${TEST_CALLS}" | tr -d ' ')" == 1 ]]
export ECODA_RUN_ID=zeroimp_bad_scope TEST_DATASETS=Other
expect_fail run_screen
unset TEST_DATASETS
export ECODA_RUN_ID=zeroimp_bad_md5
expect_fail bash "${ENTRY}" --analysis_release "${CONTROL}" --analysis_release_md5 00000000000000000000000000000000 --datasets Adams
[[ "$(wc -l < "${TEST_CALLS}" | tr -d ' ')" == 1 ]]
ecoda_artifact_owner_acquire "${SCRATCH}/inputs/Adams_composition.rds" stage5 foreign_writer 0 0 0 >/dev/null
export ECODA_RUN_ID=zeroimp_active_source
expect_fail run_screen
ecoda_artifact_owner_set_state "${SCRATCH}/inputs/Adams_composition.rds" OK fixture_release >/dev/null
[[ "$(wc -l < "${TEST_CALLS}" | tr -d ' ')" == 1 ]]
cp "${SCRATCH}/inputs/review.txt" "${TMP}/review.original"
for FAILURE in R VALIDATE OOM TAMPER; do
 export ECODA_RUN_ID="zeroimp_fail_${FAILURE}" TEST_FAIL="${FAILURE}"
 expect_fail run_screen
 F="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
 grep -q '^STATE=FAIL$' "${F}/status/terminal"
 [[ "$(cat "${F}/manifests/scheduler_ids.tsv")" == "$(printf 'ARRAY\t712345')" ]]
 if [[ "${FAILURE}" == R ]];then [[ -s "${F}/artifacts/zeroimp_count_screen/analysis/partial" ]];fi
 if [[ "${FAILURE}" == TAMPER ]];then cp "${TMP}/review.original" "${SCRATCH}/inputs/review.txt";fi
done
unset TEST_FAIL
[[ ! -f "${SCRATCH}/forbidden.calls" ]]
echo 'zeroimp count-screen gate: scalar scope/proofs/spool, create-only, failed worker/OOM/source tamper: PASS'
