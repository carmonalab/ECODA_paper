#!/bin/bash
# Actual immutable campaign/dispatcher/owners; only Slurm and fits are stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Reuse the existing pinned fixture and concurrency-safe Slurm recorder.
SETUP="$(awk '/^NAS=/ {found=1; exit} {print} END {if(!found) exit 1}' "${ROOT}/tests/test_stage5_partial_recovery.sh")"
eval "${SETUP}"
NAS="${TMP_DIR}/nas"; mkdir -p "${NAS}" "${TMP_DIR}/inputs"
"${TEST_REAL_PYTHON}" - "${SOURCE_TREE}/datasets.json" "${TMP_DIR}" <<'PY'
import hashlib,json,pathlib,sys
config=json.loads(pathlib.Path(sys.argv[1]).read_text());root=pathlib.Path(sys.argv[2])
selection=[];refs=[];matrix=[]
for ds in ('Covid19_PBMC','Kidney_KPMP_full','Lupus_PBMC','Lung'):
    view='batch_effect_uncorrected'
    assert config[ds]['use_for_batch_effect'] is True
    p=root/'inputs'/config[ds]['views'][view]['output_file_name'];p.write_text('tiny fake input; fits/preflight are stubbed\n')
    md5=hashlib.md5(p.read_bytes()).hexdigest();size=p.stat().st_size
    proof=pathlib.Path(str(p)+'.md5');proof.write_text(f'MD5={md5}\nSIZE={size}\nPATH={p}\n')
    refs.append(f'{ds}\t{view}\t{p}\t{md5}\t{size}\t{proof}\n')
    selection.append(f'{ds}\t{view}\t{view}\n');matrix.append(f'{ds}\t{view}\tmrvi\n')
(root/'selection.tsv').write_text(''.join(selection));(root/'refs.tsv').write_text(''.join(refs));(root/'scope.tsv').write_text(''.join(matrix))
PY
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
DS=Covid19_PBMC
OLD="$(ecoda_owner_acquire stage5 "uncorrected/${DS}/batch_effect_uncorrected/mrvi" old_valid_run 0 0)"
ecoda_owner_set_state "${OLD}" OK old_science
cp "${OLD}/owner" "${TMP_DIR}/old.saved"
campaign() {
  local parent="$1" child="$2" matrix="$3" gpu_policy="${4:-auto}"
  printf 'uncorrected\t%s\t%s/selection.tsv\t%s/refs.tsv\t%s\n' "${child}" "${TMP_DIR}" "${TMP_DIR}" "${matrix}" > "${TMP_DIR}/${parent}.tsv"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
    BENCHMARK_MATRIX_TEST=1 ECODA_RUN_ID="${parent}" NAS_TARGET_DIR="${NAS}" \
    bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
      --contexts "${TMP_DIR}/${parent}.tsv" --derivative-release rectangle --batch-gpu-policy "${gpu_policy}"
}
if ! campaign rectangle_parent rectangle_uncorrected "${TMP_DIR}/scope.tsv" any > "${TMP_DIR}/positive.log" 2>&1; then
  cat "${TMP_DIR}/positive.log" >&2
  cat "${HPC_ROOT}/_ecoda_runs/rectangle_parent/logs/uncorrected.log" >&2
  exit 1
fi
cmp "${OLD}/owner" "${TMP_DIR}/old.saved"
"${TEST_REAL_PYTHON}" - "${HPC_ROOT}" "${CAPTURE}" "${TMP_DIR}/scope.tsv" <<'PY'
import pathlib,shlex,sys
runs=pathlib.Path(sys.argv[1])/'_ecoda_runs';child=runs/'rectangle_uncorrected';scope=pathlib.Path(sys.argv[3]).read_bytes()
assert (child/'manifests/dispatch_selection.tsv').read_bytes()==scope
assert (child/'manifests/pending_selection.tsv').read_bytes()==scope
assert not (child/'manifests/method_matrix.tsv').exists()
metadata=dict(line.split('=',1) for line in (child/'metadata').read_text().splitlines() if line)
assert metadata['TARGET_METHODS']==metadata['METHODS']=='mrvi' and metadata['PASS']=='uncorrected' and metadata['ANALYSIS_VARIANT']=='final'
assert not (runs/'rectangle_corrected').exists() and not (runs/'rectangle_benchmark').exists()
jobs=[]
for line in pathlib.Path(sys.argv[2]).read_text().splitlines():
    args=shlex.split(line.split(' ARGS ',1)[1]);i=next(i for i,a in enumerate(args) if not a.startswith('-'))
    options=dict(a[2:].split('=',1) for a in args[:i] if a.startswith('--') and '=' in a)
    jobs.append((args[i],options))
arrays=[j for j in jobs if 'array' in j[1]]
assert len(jobs)==3 and len(arrays)==1 and arrays[0][1]['array'].split('%')[0]=='1-4'
assert 'METHOD=mrvi' in arrays[0][1]['export']
assert 'METHOD_GPU_POLICY=any' in arrays[0][1]['export']
assert 'gpus' in arrays[0][1] or 'gres' in arrays[0][1]
assert 'constraint' not in arrays[0][1]
assert (runs/'rectangle_parent/manifests/campaign_batch_gpu_policy.tsv').read_text()=='any\n'
assert sum(script.endswith('/matrix_watchdog.sh') for script,_ in jobs)==1
assert sum(script.endswith('/matrix_gate.sh') for script,_ in jobs)==1
assert 'prepare_pseudobulk' not in pathlib.Path(sys.argv[2]).read_text()
PY
# No heterogeneous rectangle, implicit PB dependency, wrong pass, unsupported
# method, duplicate or altered order may submit a child or scheduler work.
for fault in heterogeneous implicit_prep wrong_view unknown duplicate reversed; do
  "${TEST_REAL_PYTHON}" - "${TMP_DIR}/scope.tsv" "${TMP_DIR}/bad.tsv" "${fault}" <<'PY'
import pathlib,sys
rows=[line.split('\t') for line in pathlib.Path(sys.argv[1]).read_text().splitlines()];fault=sys.argv[3]
if fault=='heterogeneous':rows[-1][2]='qot'
elif fault=='implicit_prep':
    for row in rows:row[2]='composition'
elif fault=='wrong_view':rows[0][1]='batch_effect_corrected'
elif fault=='unknown':rows[0][2]='unsupported'
elif fault=='duplicate':rows.append(rows[0])
else:rows.reverse()
pathlib.Path(sys.argv[2]).write_text(''.join('\t'.join(row)+'\n' for row in rows))
PY
  : > "${CAPTURE}"
  if campaign "bad_${fault}" "badchild_${fault}" "${TMP_DIR}/bad.tsv" > "${TMP_DIR}/bad.log" 2>&1; then
    echo "invalid uncorrected scope accepted: ${fault}" >&2; exit 1
  fi
  [[ ! -s "${CAPTURE}" && ! -e "${HPC_ROOT}/_ecoda_runs/badchild_${fault}" ]]
done
: > "${CAPTURE}"
if campaign bad_gpu_policy bad_gpu_child "${TMP_DIR}/scope.tsv" cpu > "${TMP_DIR}/gpu-policy.log" 2>&1; then
  echo 'invalid batch GPU policy accepted' >&2; exit 1
fi
[[ ! -s "${CAPTURE}" && ! -e "${HPC_ROOT}/_ecoda_runs/bad_gpu_child" ]]
FOREIGN="$(ecoda_owner_acquire stage5 "derivative/collision_uncorrected/uncorrected/${DS}/batch_effect_uncorrected/mrvi" foreign_run 0 0)"
cp "${FOREIGN}/owner" "${TMP_DIR}/foreign.saved"; : > "${CAPTURE}"
if campaign collision_parent collision_uncorrected "${TMP_DIR}/scope.tsv" > "${TMP_DIR}/collision.log" 2>&1; then
  echo 'foreign derivative owner accepted' >&2; exit 1
fi
[[ ! -s "${CAPTURE}" ]]
cmp "${FOREIGN}/owner" "${TMP_DIR}/foreign.saved"
cmp "${OLD}/owner" "${TMP_DIR}/old.saved"
echo 'PASS: fresh uncorrected rectangle is exactly4GPU-any MrVI, no PB/healthy contexts; policy/scope/dependency/foreign-owner guards fail closed'
