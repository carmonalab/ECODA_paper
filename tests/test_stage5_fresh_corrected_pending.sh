#!/bin/bash
# Actual immutable matrix fixture/dispatcher/binding; Slurm and fits are stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Reuse the maintained fixture's setup, stopping before its first dispatch.
SETUP="$(awk '/^if \[\[ "\$\{ECODA_MATRIX_REQUEST_ONLY:-0\}" == 1/ {found=1; exit} {print} END {if(!found) exit 1}' "${ROOT}/tests/test_benchmark_matrix_submitter.sh")"
eval "${SETUP}"
# Execute the new pure scope check with configured Python rather than stubbing it.
cat > "${HOST_PREFIX}/bin/python" <<'SH'
#!/bin/bash
if [[ "${1:-}" == */derivative_artifacts.py || "${1:-}" == */benchmark_request_scope.py ||
      ( "${1:-}" == - && "${2:-}" == */manifests/artifact_bindings.json ) ]]; then
  exec "${TEST_REAL_PYTHON}" "$@"
fi
exit 0
SH
chmod +x "${HOST_PREFIX}/bin/python"
NAS="${TMP_DIR}/nas"
mkdir -p "${NAS}" "${TMP_DIR}/inputs"
"${TEST_REAL_PYTHON}" - "${SOURCE_TREE}/datasets.json" "${TMP_DIR}" <<'PY'
import hashlib,json,pathlib,sys
config=json.loads(pathlib.Path(sys.argv[1]).read_text());root=pathlib.Path(sys.argv[2])
datasets=[d for d,v in config.items() if not d.startswith('_') and v.get('use_for_batch_effect') is True]
assert len(datasets)==8
refs=[];selection=[]
for ds in datasets:
    view='batch_effect_corrected';p=root/'inputs'/config[ds]['views'][view]['output_file_name']
    p.write_text('small fake H5AD; scientific preflight stubbed\n');md5=hashlib.md5(p.read_bytes()).hexdigest()
    proof=pathlib.Path(str(p)+'.md5');proof.write_text(f'MD5={md5}\nSIZE={p.stat().st_size}\nPATH={p}\n')
    refs.append(f'{ds}\t{view}\t{p}\t{md5}\t{p.stat().st_size}\t{proof}\n')
    selection.append(f'{ds}\t{view}\t{view}\n')
(root/'refs.tsv').write_text(''.join(refs));(root/'selection.tsv').write_text(''.join(selection))
(root/'small.refs.tsv').write_text(refs[0]);(root/'small.selection.tsv').write_text(selection[0])
PY
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
FIRST_DS="$(head -1 "${TMP_DIR}/selection.tsv" | cut -f1)"
LEGACY_OWNER="$(ecoda_owner_acquire stage5 "corrected/${FIRST_DS}/batch_effect_corrected/mrvi" legacy_run 0 0)"
ecoda_owner_set_state "${LEGACY_OWNER}" OK legacy_fixture
cp "${LEGACY_OWNER}/owner" "${TMP_DIR}/legacy-owner.saved"
dispatch_wave() {
  local ID="$1" PREFIX="$2"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
    BENCHMARK_MATRIX_TEST=1 ECODA_RUN_ID="${ID}" NAS_TARGET_DIR="${NAS}" \
    bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
      --selection-file "${TMP_DIR}/${PREFIX}selection.tsv" --pass corrected --analysis-variant corrected_final \
      --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot \
      --derivative-release fresh-corrected-test --input-artifacts "${TMP_DIR}/${PREFIX}refs.tsv"
}
if ! dispatch_wave fresh_corrected '' > "${TMP_DIR}/full.log" 2>&1; then cat "${TMP_DIR}/full.log" >&2; exit 1; fi
cmp "${LEGACY_OWNER}/owner" "${TMP_DIR}/legacy-owner.saved"
RUN="${HPC_ROOT}/_ecoda_runs/fresh_corrected"
"${TEST_REAL_PYTHON}" - "${RUN}" "${CAPTURE}" <<'PY'
import pathlib,shlex,sys
run=pathlib.Path(sys.argv[1]);rows=(run/'manifests/dispatch_selection.tsv').read_text().splitlines()
assert len(rows)==len(set(rows))==56
pending=(run/'manifests/classified_pending.tsv').read_text().splitlines();assert set(pending)==set(rows) and len(pending)==56
jobs=[]
for line in pathlib.Path(sys.argv[2]).read_text().splitlines():
    args=shlex.split(line.split(' ARGS ',1)[1]);i=next(i for i,a in enumerate(args) if not a.startswith('-'))
    opts=dict(a[2:].split('=',1) for a in args[:i] if a.startswith('--') and '=' in a)
    jobs.append((args[i],opts))
arrays=[j for j in jobs if 'array' in j[1]]
assert len(arrays)==7 and sum(int(j[1]['array'].split('%')[0].split('-')[-1]) for j in arrays)==56
assert len([j for j in jobs if j[0].endswith('/matrix_watchdog.sh')])==7
assert len([j for j in jobs if j[0].endswith('/matrix_gate.sh')])==1
assert len(jobs)==15
methods={dict(x.split('=',1) for x in o['export'].split(',') if '=' in x)['METHOD'] for _,o in arrays}
assert methods=={'prepare_pseudobulk','pseudobulk','gloscope','composition','mrvi','pilot','qot'}
assert len((run/'manifests/scheduler_ids.tsv').read_text().splitlines())==15
PY
# The binding permits an explicit smaller fresh cohort selection as well.
: > "${CAPTURE}"
if ! dispatch_wave small_corrected small. > "${TMP_DIR}/small.log" 2>&1; then cat "${TMP_DIR}/small.log" >&2; exit 1; fi
[[ "$(wc -l < "${HPC_ROOT}/_ecoda_runs/small_corrected/manifests/classified_pending.tsv" | tr -d ' ')" == 7 ]]
cmp "${LEGACY_OWNER}/owner" "${TMP_DIR}/legacy-owner.saved"
# A foreign/active owner in the actual derivative namespace still blocks work.
COLLISION_OWNER="$(ecoda_owner_acquire stage5 "derivative/collision_corrected/corrected/${FIRST_DS}/batch_effect_corrected/mrvi" foreign_run 0 0)"
cp "${COLLISION_OWNER}/owner" "${TMP_DIR}/collision-owner.saved"
: > "${CAPTURE}"
if dispatch_wave collision_corrected small. > "${TMP_DIR}/collision.log" 2>&1; then
  echo 'ERROR: derivative owner collision was allowed' >&2; exit 1
fi
[[ ! -s "${CAPTURE}" ]]
cmp "${COLLISION_OWNER}/owner" "${TMP_DIR}/collision-owner.saved"

# Pure classification remains read-only and never broadens repair/reuse scope.
(
  export PROJECT_ROOT="${SOURCE_TREE}" DATASETS_JSON_FILE="${SOURCE_TREE}/datasets.json" NAS_TARGET_DIR="${NAS}"
  export PYTHON_BIN="${HOST_PREFIX}/bin/python" ECODA_RUN_ROOT="${RUN}" ECODA_RUN_ID=fresh_corrected
  export ECODA_ARTIFACT_BINDING="${RUN}/manifests/artifact_bindings.json"
  export ANALYSIS_ROOT="${RUN}/artifacts/batch_effect/corrected_final"
  export ANALYSIS_NAS_ROOT="${NAS}/releases/fresh-corrected-test/fresh_corrected/batch_effect/corrected_final"
  export PASS_ARG=corrected ANALYSIS_VARIANT=corrected_final
  METHODS_ARG=prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
  ECODA_RECOVERY_ACTIVE=0 FORCE_ARG=0 FORCE_TARGETED_ARG=0 TARGET_METHODS_SET=0 METHOD_MATRIX_MODE=0
  SYNC_ONLY_RUN="" EXPECTED_PENDING_ARG="" ECODA_EXPECTED_PENDING="" PRECOMPUTE_BACKUP_RECEIPT_ARG="" ECODA_PRECOMPUTE_BACKUP_RECEIPT=""
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  source "${ROOT}/src/utils/bash/ecoda_derivative_artifacts.sh"
  source "${ROOT}/src/5_run_benchmark_methods/stage5_pending.sh"
  stage5_method_is_forced() { return 1; }
  benchmark_selected_artifacts_valid() { [[ "${VALID_PYTHON:-0}" == 1 && "$3" == mrvi ]]; }
  stage5_prepare_pseudobulk_valid() { return 1; }
  benchmark_artifacts_for() { ARTIFACT_PATHS=("${ANALYSIS_ROOT}/$1-$3.rds"); }
  stage5_pending_owner_valid() { return 0; }
  stage5_pending_target_valid() { return 0; }
  DISPATCH="${RUN}/manifests/dispatch_selection.tsv"
  REPORT="${TMP_DIR}/classified.tsv"
  fail() { if "$@" > "${TMP_DIR}/negative.log" 2>&1; then echo 'ERROR: invalid fresh/legacy classification passed' >&2; exit 1; fi; }
  stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"
  [[ "$(wc -l < "${REPORT}" | tr -d ' ')" == 56 ]]
  VALID_PYTHON=1
  stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"
  [[ "$(wc -l < "${REPORT}" | tr -d ' ')" == 48 ]]
  VALID_PYTHON=0
  BINDING="${ECODA_ARTIFACT_BINDING}"
  ECODA_ARTIFACT_BINDING=""; fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"
  ECODA_ARTIFACT_BINDING="${HPC_ROOT}/_ecoda_runs/small_corrected/manifests/artifact_bindings.json"
  fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"
  ECODA_ARTIFACT_BINDING="${BINDING}"
  ECODA_RECOVERY_ACTIVE=1; fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"; ECODA_RECOVERY_ACTIVE=0
  EXPECTED_PENDING_ARG=historical-residual; fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"; EXPECTED_PENDING_ARG=""
  SYNC_ONLY_RUN=historical; fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"; SYNC_ONLY_RUN=""
  FORCE_ARG=1; fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"; FORCE_ARG=0
  METHOD_MATRIX_MODE=1; fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"; METHOD_MATRIX_MODE=0
  cp "${DISPATCH}" "${TMP_DIR}/dispatch.saved"; chmod u+w "${DISPATCH}"
  printf 'Diabetes\tbatch_effect_corrected\tmrvi\n' >> "${DISPATCH}"
  fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"
  cp "${TMP_DIR}/dispatch.saved" "${DISPATCH}"; chmod a-w "${DISPATCH}"
  # Existing invalid targets are still protected; no unlink/invalidation occurs.
  DS="$(head -1 "${DISPATCH}" | cut -f1)"
  mkdir -p "${ANALYSIS_ROOT}"; printf keep > "${ANALYSIS_ROOT}/${DS}-mrvi.rds"
  fail stage5_classify_pending_rows "${DISPATCH}" "${REPORT}"
  [[ "$(cat "${ANALYSIS_ROOT}/${DS}-mrvi.rds")" == keep ]]
)
echo 'fresh corrected pending: 56 tasks, seven arrays/watchdogs, one aggregate; strict legacy/target protections passed'
