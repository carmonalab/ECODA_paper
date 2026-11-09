#!/bin/bash
# Reuse the matrix test's pinned source/runtime fixture; replace only its child
# dispatcher before the source archive is sealed. No Slurm or NAS operations.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE="${ROOT}/tests/test_benchmark_matrix_submitter.sh"
export TEST_EXPECTED_BATCH_METHODS="$(sed -n 's/^EXPECTED_BATCH_METHODS="\(.*\)"$/\1/p' "${ROOT}/src/5_run_benchmark_methods/stage5_dispatcher.sh")"
eval "$(awk '/^ROOT=/ {next} /^tar -cf / {exit} {print}' "${FIXTURE}")"
printf '#!/bin/bash\nexec "${TEST_REAL_PYTHON}" "$@"\n' > "${HOST_PREFIX}/bin/python"
chmod +x "${HOST_PREFIX}/bin/python"
cat > "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_dispatcher.sh" <<'CHILD'
#!/bin/bash
set -euo pipefail
[[ -z "${ECODA_RUN_ROOT:-}${ECODA_ARTIFACT_BINDING:-}${ANALYSIS_ROOT:-}${ANALYSIS_PASS:-}${ECODA_SELECTION_MANIFEST:-}" ]]
source "${ECODA_SOURCE_ROOT}/src/slurm_config.sh"
source "${ECODA_SOURCE_ROOT}/src/utils/bash/ecoda_run_common.sh"
child="${ECODA_RUN_ID}"; context="${child#*_}"
selection=""; methods_seen=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --selection-file) selection="$2" ;;
    --pass) [[ "$2" == "${context}" ]] ;;
    --analysis-variant) case "${context}:$2" in uncorrected:final|corrected:corrected_final) ;; *) exit 1;; esac ;;
    --methods) [[ "$2" == "${TEST_EXPECTED_BATCH_METHODS}" ]]; methods_seen=1 ;;
    --derivative-release|--input-artifacts) ;;
    *) exit 1 ;;
  esac
  shift 2
done
[[ "${selection}" == */manifests/campaign_"${context}".selection.tsv ]]
[[ "${context}" == benchmark || ${methods_seen} -eq 1 ]]
if [[ "${context}" == benchmark ]]; then
  grep -q $'Stephenson\tbenchmark_analysis\ttrans' "${selection}"
  grep -q $'Stephenson\tbenchmark_analysis\tzeroimp' "${selection}"
fi
ecoda_init_run stage5 "${child}" >/dev/null
case "${context}" in benchmark) scheduler=101;; uncorrected) scheduler=102;; corrected) scheduler=103;; esac
printf 'ARRAY\t%s\n' "${scheduler}" > "${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
mkdir -p "${CAMPAIGN_BARRIER}"
touch "${CAMPAIGN_BARRIER}/${context}.started"
ready=0
for attempt in {1..100}; do
  if [[ $(find "${CAMPAIGN_BARRIER}" -name '*.started' | wc -l | tr -d ' ') == 3 ]]; then ready=1; break; fi
  sleep 0.02
done
[[ ${ready} -eq 1 ]]
if [[ "${CAMPAIGN_FAIL_CONTEXT:-}" == "${context}" ]]; then
  ecoda_set_run_state FAIL 'deliberate isolated child failure'; exit 7
fi
[[ "${CAMPAIGN_MISSING_TERMINAL:-}" != "${context}" ]] || exit 0
ecoda_set_run_state OK 'fixture completed with all peers started'
printf 'BENCHMARK_SCHEDULER_ID=%s:%s\n' "${context}" "${scheduler}"
CHILD
eval "$(awk '/^tar -cf / {copy=1} /^if \[\[ "\$\{ECODA_MATRIX_REQUEST_ONLY/ {exit} copy {print}' "${FIXTURE}")"
NAS_BASE="${TMP_DIR}/campaign-nas"; mkdir -p "${NAS_BASE}"
export CAMPAIGN_BARRIER="${TMP_DIR}/barrier"
"${TEST_REAL_PYTHON}" - "${SOURCE_TREE}/datasets.json" "${TMP_DIR}" <<'PY'
import hashlib,json,pathlib,sys
config=json.loads(pathlib.Path(sys.argv[1]).read_text()); root=pathlib.Path(sys.argv[2])
for context,view,lines in (
    ('benchmark','benchmark_analysis','Stephenson\tbenchmark_analysis\tall\nStephenson\tbenchmark_analysis\ttrans\nStephenson\tbenchmark_analysis\tzeroimp\n'),
    ('uncorrected','batch_effect_uncorrected','Stephenson\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n'),
    ('corrected','batch_effect_corrected','Stephenson\tbatch_effect_corrected\tbatch_effect_corrected\n')):
    source=root/config['Stephenson']['views'][view]['output_file_name']
    source.write_bytes(context.encode()); md5=hashlib.md5(source.read_bytes()).hexdigest(); size=source.stat().st_size
    proof=pathlib.Path(f'{source}.md5'); proof.write_text(f'MD5={md5}\nSIZE={size}\nPATH={source}\n')
    (root/f'{context}.selection.tsv').write_text(lines)
    (root/f'{context}.inputs.tsv').write_text(f'Stephenson\t{view}\t{source}\t{md5}\t{size}\t{proof}\n')
PY
make_contexts() {
  local prefix="$1" context
  for context in benchmark uncorrected corrected; do
    printf '%s\t%s_%s\t%s/%s.selection.tsv\t%s/%s.inputs.tsv\n' \
      "${context}" "${prefix}" "${context}" "${TMP_DIR}" "${context}" "${TMP_DIR}" "${context}"
  done
}
make_contexts child > "${TMP_DIR}/contexts.tsv"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
  ECODA_RUN_ID=campaign_success ANALYSIS_ROOT=/unrelated/parent ANALYSIS_PASS=parent \
  ECODA_ARTIFACT_BINDING=/unrelated/parent/binding ECODA_SELECTION_MANIFEST=/unrelated/parent/selection \
  bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
  --contexts "${TMP_DIR}/contexts.tsv" --derivative-release test-release > "${TMP_DIR}/success.log"
parent="${HPC_ROOT}/_ecoda_runs/campaign_success"
grep -q '^STATE=OK$' "${parent}/status/terminal"
[[ $(wc -l < "${parent}/manifests/campaign_scheduler_ids.tsv" | tr -d ' ') == 3 ]]
cmp "${TMP_DIR}/benchmark.selection.tsv" "${parent}/manifests/campaign_benchmark.selection.tsv"
[[ $(find "${CAMPAIGN_BARRIER}" -name '*.started' | wc -l | tr -d ' ') == 3 ]]
"${TEST_REAL_PYTHON}" - "${parent}" <<'PY'
import hashlib,pathlib,sys
root=pathlib.Path(sys.argv[1]); manifest=root/'manifests'
identities=[line.split('\t') for line in (manifest/'campaign_identities.tsv').read_text().splitlines()]
assert len(identities)==6
for name,expected in identities:
    target=pathlib.Path(name)
    assert target.parent==manifest and hashlib.sha256(target.read_bytes()).hexdigest()==expected
for name in ('source.manifest','runtime.identity','campaign_contexts.tsv'):
    assert (manifest/f'{name}.md5').is_file()
PY

export CAMPAIGN_BARRIER="${TMP_DIR}/failure-barrier"
make_contexts retry > "${TMP_DIR}/retry.tsv"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
    ECODA_RUN_ID=campaign_failure CAMPAIGN_FAIL_CONTEXT=benchmark \
    bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
    --contexts "${TMP_DIR}/retry.tsv" --derivative-release test-release > "${TMP_DIR}/failure.log" 2>&1; then
  echo 'campaign accepted a failed child' >&2; exit 1
fi
grep -q '^STATE=FAIL$' "${HPC_ROOT}/_ecoda_runs/campaign_failure/status/terminal"
for context in uncorrected corrected; do
  grep -q '^STATE=OK$' "${HPC_ROOT}/_ecoda_runs/retry_${context}/status/terminal"
done
[[ $(wc -l < "${HPC_ROOT}/_ecoda_runs/campaign_failure/manifests/campaign_scheduler_ids.tsv" | tr -d ' ') == 3 ]]

export CAMPAIGN_BARRIER="${TMP_DIR}/missing-status-barrier"
make_contexts missing > "${TMP_DIR}/missing.tsv"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
    ECODA_RUN_ID=campaign_missing_status CAMPAIGN_MISSING_TERMINAL=benchmark \
    bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
    --contexts "${TMP_DIR}/missing.tsv" --derivative-release test-release > "${TMP_DIR}/missing.log" 2>&1; then
  echo 'campaign accepted a zero-exit child without terminal evidence' >&2; exit 1
fi
grep -q '^STATE=OK$' "${HPC_ROOT}/_ecoda_runs/missing_corrected/status/terminal"

make_contexts bad > "${TMP_DIR}/bad.tsv"
sed 's/^corrected/uncorrected/' "${TMP_DIR}/bad.tsv" > "${TMP_DIR}/duplicate.tsv"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
    ECODA_RUN_ID=campaign_invalid bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
    --contexts "${TMP_DIR}/duplicate.tsv" --derivative-release test-release > "${TMP_DIR}/invalid.log" 2>&1; then
  echo 'campaign accepted duplicate contexts' >&2; exit 1
fi
for context in benchmark uncorrected corrected; do [[ ! -e "${HPC_ROOT}/_ecoda_runs/bad_${context}" ]]; done
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${TMP_DIR}/absent-nas" \
    ECODA_RUN_ID=campaign_unmounted bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
    --contexts "${TMP_DIR}/bad.tsv" --derivative-release test-release > "${TMP_DIR}/unmounted.log" 2>&1; then
  echo 'campaign controller accepted an absent NAS mount' >&2; exit 1
fi
[[ ! -e "${HPC_ROOT}/_ecoda_runs/campaign_unmounted" && ! -e "${TMP_DIR}/absent-nas" ]]
echo 'PASS: one pinned parent, concurrent isolated contexts, frozen exact scopes, child-failure isolation and collected scheduler IDs'
