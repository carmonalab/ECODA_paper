#!/bin/bash
# Actual campaign + dispatchers + bindings/ownership; only Slurm and fits stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SETUP="$(awk '/^if \[\[ "\$\{ECODA_MATRIX_REQUEST_ONLY:-0\}" == 1/ {found=1; exit} {print} END {if(!found) exit 1}' "${ROOT}/tests/test_benchmark_matrix_submitter.sh")"
eval "${SETUP}"
cat > "${HOST_PREFIX}/bin/python" <<'SH'
#!/bin/bash
if [[ "${1:-}" == */derivative_artifacts.py || "${1:-}" == */benchmark_request_scope.py ||
      ( "${1:-}" == - && ( "${2:-}" == */manifests/artifact_bindings.json || "${2:-}" == */_ecoda_runs/* ) ) ]]; then
  exec "${TEST_REAL_PYTHON}" "$@"
fi
exit 0
SH
chmod +x "${HOST_PREFIX}/bin/python"
# A campaign submits concurrently: allocate unique stub IDs under one lock.
cat > "${TMP_DIR}/bin/sbatch" <<'SH'
#!/bin/bash
exec "${TEST_REAL_PYTHON}" - "$@" <<'PY'
import fcntl,os,pathlib,shlex,sys
target=pathlib.Path(os.environ['CAPTURE'])
with pathlib.Path(str(target)+'.lock').open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX)
    previous=target.read_text() if target.exists() else ''
    identifier=80000+len(previous.splitlines())+1
    keys=('ECODA_SOURCE_ROOT','ECODA_SOURCE_MANIFEST','ECODA_RUNTIME_IMAGE','ECODA_RUNTIME_MANIFEST','ECODA_RUN_ROOT','ECODA_RUN_ID')
    with target.open('a') as out:
        out.write('ENV '+' '.join(k+'='+os.environ.get(k,'') for k in keys)+' ARGS '+shlex.join(sys.argv[1:])+'\n')
print(identifier)
PY
SH
chmod +x "${TMP_DIR}/bin/sbatch"
NAS="${TMP_DIR}/nas"; mkdir -p "${NAS}" "${TMP_DIR}/inputs"
"${TEST_REAL_PYTHON}" - "${SOURCE_TREE}/datasets.json" "${TMP_DIR}" <<'PY'
import hashlib,json,pathlib,sys
config=json.loads(pathlib.Path(sys.argv[1]).read_text());root=pathlib.Path(sys.argv[2])
datasets=[ds for ds,v in config.items() if not ds.startswith('_') and v.get('use_for_batch_effect') is True]
assert len(datasets)==8
refs=[];selected=[];matrix=[]
for ds,view in [(ds,'batch_effect_corrected') for ds in datasets]+[('Stephenson','benchmark_analysis')]:
    path=root/'inputs'/config[ds]['views'][view]['output_file_name']
    path.write_text('tiny fixture; scientific preflight is stubbed\n')
    md5=hashlib.md5(path.read_bytes()).hexdigest();size=path.stat().st_size
    proof=pathlib.Path(str(path)+'.md5');proof.write_text(f'MD5={md5}\nSIZE={size}\nPATH={path}\n')
    ref=f'{ds}\t{view}\t{path}\t{md5}\t{size}\t{proof}\n'
    if view=='benchmark_analysis':
        (root/'benchmark.inputs.tsv').write_text(ref)
        (root/'benchmark.selection.tsv').write_text(f'{ds}\t{view}\tgloscope\thvg2000_pcadims30\n')
    else:
        refs.append(ref);selected.append(f'{ds}\t{view}\t{view}\n');matrix.append(f'{ds}\t{view}\tqot\n')
(root/'corrected.inputs.tsv').write_text(''.join(refs))
(root/'corrected.selection.tsv').write_text(''.join(selected));(root/'corrected.matrix.tsv').write_text(''.join(matrix))
PY
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
FIRST_DS="$(head -1 "${TMP_DIR}/corrected.selection.tsv" | cut -f1)"
OLD_OWNER="$(ecoda_owner_acquire stage5 "corrected/${FIRST_DS}/batch_effect_corrected/qot" old_valid_run 0 0)"
ecoda_owner_set_state "${OLD_OWNER}" OK old_valid_output
cp "${OLD_OWNER}/owner" "${TMP_DIR}/old-owner.saved"
contexts() {
  printf 'benchmark\t%s_benchmark\t%s/benchmark.selection.tsv\t%s/benchmark.inputs.tsv\n' "$1" "${TMP_DIR}" "${TMP_DIR}"
  printf 'corrected\t%s_corrected\t%s/corrected.selection.tsv\t%s/corrected.inputs.tsv\t%s/corrected.matrix.tsv\n' "$1" "${TMP_DIR}" "${TMP_DIR}" "${TMP_DIR}"
}
contexts recovery > "${TMP_DIR}/contexts.tsv"
if ! HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
    BENCHMARK_MATRIX_TEST=1 ECODA_RUN_ID=partial_campaign NAS_TARGET_DIR="${NAS}" \
    bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
      --contexts "${TMP_DIR}/contexts.tsv" --derivative-release partial-recovery > "${TMP_DIR}/campaign.log" 2>&1; then
  cat "${TMP_DIR}/campaign.log" >&2
  for path in "${HPC_ROOT}/_ecoda_runs/partial_campaign/logs/"*.log; do cat "${path}" >&2; done
  exit 1
fi
cmp "${OLD_OWNER}/owner" "${TMP_DIR}/old-owner.saved"
"${TEST_REAL_PYTHON}" - "${HPC_ROOT}" "${CAPTURE}" <<'PY'
import pathlib,shlex,sys
root=pathlib.Path(sys.argv[1]);runs=root/'_ecoda_runs'
assert not (runs/'recovery_uncorrected').exists()
expected={'recovery_benchmark':('gloscope',1),'recovery_corrected':('qot',8)}
arrays=[];jobs=[]
for line in pathlib.Path(sys.argv[2]).read_text().splitlines():
    args=shlex.split(line.split(' ARGS ',1)[1]);i=next(i for i,a in enumerate(args) if not a.startswith('-'))
    options=dict(a[2:].split('=',1) for a in args[:i] if a.startswith('--') and '=' in a)
    env=dict(a.split('=',1) for a in line.split(' ARGS ',1)[0].split()[1:])
    assert env['ECODA_RUN_ID'] in expected
    jobs.append(args[i])
    if 'array' in options:
        exported=dict(a.split('=',1) for a in options['export'].split(',') if '=' in a)
        method,count=expected[env['ECODA_RUN_ID']]
        assert exported['METHOD']==method
        assert int(options['array'].split('%')[0].split('-')[-1])==count
        arrays.append((method,count))
assert sorted(arrays)==[('gloscope',1),('qot',8)] and sum(n for _,n in arrays)==9
assert len(jobs)==6 and sum(p.endswith('/matrix_watchdog.sh') for p in jobs)==2
assert sum(p.endswith('/matrix_gate.sh') for p in jobs)==2
for child,(method,count) in expected.items():
    folder=runs/child/'manifests'
    dispatch=[line.split('\t') for line in (folder/'dispatch_selection.tsv').read_text().splitlines()]
    assert len(dispatch)==count and all(row[2]==method for row in dispatch)
    assert not any(row[2]=='prepare_pseudobulk' for row in dispatch)
assert len((runs/'partial_campaign/manifests/campaign_scheduler_ids.tsv').read_text().splitlines())==6
PY
# The fresh exception is not a public bypass: exact binding, sealed matrix,
# ordered dispatch, scope/count identity and all legacy-mode guards are required.
(
  export PROJECT_ROOT="${SOURCE_TREE}" DATASETS_JSON_FILE="${SOURCE_TREE}/datasets.json" NAS_TARGET_DIR="${NAS}"
  export PYTHON_BIN="${HOST_PREFIX}/bin/python" ECODA_RUN_ROOT="${HPC_ROOT}/_ecoda_runs/recovery_corrected" ECODA_RUN_ID=recovery_corrected
  export ECODA_ARTIFACT_BINDING="${ECODA_RUN_ROOT}/manifests/artifact_bindings.json"
  export ANALYSIS_ROOT="${ECODA_RUN_ROOT}/artifacts/batch_effect/corrected_final"
  export ANALYSIS_NAS_ROOT="${NAS}/releases/partial-recovery/recovery_corrected/batch_effect/corrected_final"
  export PASS_ARG=corrected ANALYSIS_VARIANT=corrected_final
  METHODS_ARG=prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
  ECODA_RECOVERY_ACTIVE=0 FORCE_ARG=0 FORCE_TARGETED_ARG=0 TARGET_METHODS_SET=0 METHOD_MATRIX_MODE=1
  SYNC_ONLY_RUN="" EXPECTED_PENDING_ARG="" ECODA_EXPECTED_PENDING="" PRECOMPUTE_BACKUP_RECEIPT_ARG="" ECODA_PRECOMPUTE_BACKUP_RECEIPT=""
  source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
  source "${ROOT}/src/utils/bash/ecoda_derivative_artifacts.sh"
  source "${ROOT}/src/5_run_benchmark_methods/stage5_pending.sh"
  METHOD_MATRIX="${ECODA_RUN_ROOT}/manifests/method_matrix.tsv"
  export ECODA_STAGE5_METHOD_MATRIX="${METHOD_MATRIX}"
  METHOD_MATRIX_SHA256="$(ecoda_sha256_file "${METHOD_MATRIX}")"; METHOD_MATRIX_IDENTITY="${METHOD_MATRIX_SHA256}"
  ecoda_validate_checksum "${METHOD_MATRIX}"
  METHOD_MATRIX_MD5="${ECODA_CHECKSUM_MD5}"; METHOD_MATRIX_SIZE="${ECODA_CHECKSUM_SIZE}"; METHOD_MATRIX_DECLARED_COUNT=8
  DISPATCH="${ECODA_RUN_ROOT}/manifests/dispatch_selection.tsv"
  fail() { if stage5_pending_fresh_corrected_derivative_valid "${DISPATCH}" > "${TMP_DIR}/predicate-negative.log" 2>&1; then echo 'unsafe partial scope accepted' >&2; exit 1; fi; }
  stage5_pending_fresh_corrected_derivative_valid "${DISPATCH}"
  METHOD_MATRIX_DECLARED_COUNT=7; fail; METHOD_MATRIX_DECLARED_COUNT=8
  METHOD_MATRIX_MD5=00000000000000000000000000000000; fail
  ecoda_validate_checksum "${METHOD_MATRIX}"; METHOD_MATRIX_MD5="${ECODA_CHECKSUM_MD5}"
  TARGET_METHODS_SET=1; fail; TARGET_METHODS_SET=0
  ECODA_RECOVERY_ACTIVE=1; fail; ECODA_RECOVERY_ACTIVE=0
  FORCE_ARG=1; fail; FORCE_ARG=0
  SYNC_ONLY_RUN=historical; fail; SYNC_ONLY_RUN=""
  EXPECTED_PENDING_ARG=legacy; fail; EXPECTED_PENDING_ARG=""
  BINDING="${ECODA_ARTIFACT_BINDING}"; ECODA_ARTIFACT_BINDING=""; fail; ECODA_ARTIFACT_BINDING="${BINDING}"
  cp "${DISPATCH}" "${TMP_DIR}/dispatch.saved"; chmod u+w "${DISPATCH}"
  printf 'Diabetes\tbatch_effect_corrected\tqot\n' >> "${DISPATCH}"; fail
  cp "${TMP_DIR}/dispatch.saved" "${DISPATCH}"; chmod a-w "${DISPATCH}"
  cp "${METHOD_MATRIX}" "${TMP_DIR}/matrix.saved"; chmod u+w "${METHOD_MATRIX}"
  printf 'Diabetes\tbatch_effect_corrected\tqot\n' >> "${METHOD_MATRIX}"; fail
  cp "${TMP_DIR}/matrix.saved" "${METHOD_MATRIX}"; chmod a-w "${METHOD_MATRIX}"
  stage5_pending_fresh_corrected_derivative_valid "${DISPATCH}"
)
# A foreign owner at the new derivative row namespace rejects the child before
# any fit; the selected healthy peer still completes independently.
FOREIGN="$(ecoda_owner_acquire stage5 "derivative/collision_corrected/corrected/${FIRST_DS}/batch_effect_corrected/qot" foreign_run 0 0)"
cp "${FOREIGN}/owner" "${TMP_DIR}/foreign-owner.saved"
contexts collision > "${TMP_DIR}/collision.tsv"; : > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
    BENCHMARK_MATRIX_TEST=1 ECODA_RUN_ID=collision_campaign NAS_TARGET_DIR="${NAS}" \
    bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
      --contexts "${TMP_DIR}/collision.tsv" --derivative-release partial-recovery > "${TMP_DIR}/collision.log" 2>&1; then
  echo 'foreign derivative owner accepted' >&2; exit 1
fi
cmp "${FOREIGN}/owner" "${TMP_DIR}/foreign-owner.saved"
cmp "${OLD_OWNER}/owner" "${TMP_DIR}/old-owner.saved"
grep -q '^STATE=OK$' "${HPC_ROOT}/_ecoda_runs/collision_benchmark/status/terminal"
[[ $(wc -l < "${CAPTURE}" | tr -d ' ') == 3 ]]
! grep -q 'ECODA_RUN_ID=collision_corrected ARGS' "${CAPTURE}"
echo 'PASS: exact selected recovery has nine scientific tasks/two arrays, no healthy context or PB prep; foreign owner fails in isolation'
