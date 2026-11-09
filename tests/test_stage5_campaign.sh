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
# Exercise the production scheduler recorder, rather than pre-sealing a stub
# manifest that hides dispatcher failure before the normal submission tail.
awk '/^stage5_record_scheduler\(\)/ {copy=1} copy {print} copy && /^}/ {exit}' \
  "${ROOT}/src/5_run_benchmark_methods/stage5_dispatcher.sh" \
  > "${SOURCE_TREE}/src/5_run_benchmark_methods/scheduler_recorder_fixture.sh"
grep -q '^stage5_record_scheduler()' "${SOURCE_TREE}/src/5_run_benchmark_methods/scheduler_recorder_fixture.sh"
cat > "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_dispatcher.sh" <<'CHILD'
#!/bin/bash
set -euo pipefail
[[ -z "${ECODA_RUN_ROOT:-}${ECODA_ARTIFACT_BINDING:-}${ANALYSIS_ROOT:-}${ANALYSIS_PASS:-}${ECODA_SELECTION_MANIFEST:-}" ]]
source "${ECODA_SOURCE_ROOT}/src/slurm_config.sh"
source "${ECODA_SOURCE_ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ECODA_SOURCE_ROOT}/src/5_run_benchmark_methods/scheduler_recorder_fixture.sh"
child="${ECODA_RUN_ID}"; context="${child#*_}"
selection=""; methods_seen=0; prepared=""; prepared_sha=""; method_matrix=""; gpu_policy=""; excluded_nodes=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --selection-file) selection="$2" ;;
    --pass) [[ "$2" == "${context}" ]] ;;
    --analysis-variant) case "${context}:$2" in uncorrected:final|corrected:corrected_final) ;; *) exit 1;; esac ;;
    --methods) [[ "$2" == "${TEST_EXPECTED_BATCH_METHODS}" ]]; methods_seen=1 ;;
    --derivative-release|--input-artifacts) ;;
    --prepared-adoption) prepared="$2" ;;
    --prepared-adoption-sha256) prepared_sha="$2" ;;
    --method-matrix) method_matrix="$2" ;;
    --gpu-policy) [[ "${context}" != benchmark ]]; gpu_policy="$2" ;;
    --exclude-nodes) excluded_nodes="$2" ;;
    *) exit 1 ;;
  esac
  shift 2
done
[[ "${excluded_nodes}" == "${CAMPAIGN_EXPECT_EXCLUDED_NODES:-}" ]]
[[ "${selection}" == */manifests/campaign_"${context}".selection.tsv ]]
[[ "${context}" == benchmark || ${methods_seen} -eq 1 ]]
if [[ "${context}" == benchmark || "${CAMPAIGN_EXPECT_BATCH_GPU_POLICY:-auto}" == auto ]]; then
  [[ -z "${gpu_policy}" ]]
else
  [[ "${gpu_policy}" == "${CAMPAIGN_EXPECT_BATCH_GPU_POLICY}" ]]
  mkdir -p "${CAMPAIGN_BARRIER}"; printf '%s\n' "${gpu_policy}" > "${CAMPAIGN_BARRIER}/${context}.gpu_policy"
fi
if [[ "${CAMPAIGN_EXPECT_ADOPTION:-0}" == 1 &&
      ",${CAMPAIGN_EXPECT_ADOPTION_CONTEXTS:-benchmark,uncorrected,corrected}," == *,"${context}",* ]]; then
  [[ "${prepared}" == */manifests/campaign_prepared_adoption.json &&
     "$(ecoda_sha256_file "${prepared}")" == "${prepared_sha}" ]]
  ecoda_validate_checksum "${prepared}"
  mkdir -p "${CAMPAIGN_BARRIER}"
  printf '%s\t%s\t%s\n' "${child}" "${prepared}" "${prepared_sha}" > "${CAMPAIGN_BARRIER}/${context}.adoption"
else
  [[ -z "${prepared}${prepared_sha}" ]]
fi
if [[ -n "${method_matrix}" ]]; then
  [[ "${context}" == corrected && "${method_matrix}" == */manifests/campaign_corrected.method_matrix.tsv ]]
  ecoda_validate_checksum "${method_matrix}"
  printf '%s\n' "${method_matrix}" > "${CAMPAIGN_BARRIER}.matrix"
fi
if [[ "${context}" == benchmark ]]; then
  grep -q $'Stephenson\tbenchmark_analysis\ttrans' "${selection}"
  grep -q $'Stephenson\tbenchmark_analysis\tzeroimp' "${selection}"
fi
ecoda_init_run stage5 "${child}" >/dev/null
case "${context}" in benchmark) scheduler=101;; uncorrected) scheduler=102;; corrected) scheduler=103;; esac
SCHEDULER_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
stage5_record_scheduler ARRAY "${scheduler}"
ecoda_validate_checksum "${SCHEDULER_FILE}"
mkdir -p "${CAMPAIGN_BARRIER}"
touch "${CAMPAIGN_BARRIER}/${context}.started"
ready=0
for attempt in {1..100}; do
  if [[ $(find "${CAMPAIGN_BARRIER}" -name '*.started' | wc -l | tr -d ' ') == "${CAMPAIGN_EXPECT_CHILDREN:-3}" ]]; then ready=1; break; fi
  sleep 0.02
done
[[ ${ready} -eq 1 ]]
if [[ "${CAMPAIGN_FAIL_CONTEXT:-}" == "${context}" ]]; then
  for record in METADATA_EXPORT:201 PREFLIGHT:301 ARRAY:401; do
    stage5_record_scheduler "${record%%:*}" "${record#*:}"
    ecoda_validate_checksum "${SCHEDULER_FILE}"
  done
  ecoda_set_run_state FAIL 'deliberate isolated child failure'; exit 7
fi
if [[ "${CAMPAIGN_BAD_PROOF_CONTEXT:-}" == "${context}" ]]; then
  printf 'invalid checksum evidence\n' > "${SCHEDULER_FILE}.md5"
  before_manifest="$(ecoda_sha256_file "${SCHEDULER_FILE}")"
  before_proof="$(ecoda_sha256_file "${SCHEDULER_FILE}.md5")"
  if stage5_record_scheduler ARRAY "${scheduler}"; then
    echo 'duplicate scheduler record accepted corrupt proof' >&2; exit 1
  fi
  if stage5_record_scheduler PREFLIGHT 999; then
    echo 'new scheduler record repaired corrupt proof' >&2; exit 1
  fi
  [[ "$(ecoda_sha256_file "${SCHEDULER_FILE}")" == "${before_manifest}" &&
     "$(ecoda_sha256_file "${SCHEDULER_FILE}.md5")" == "${before_proof}" ]]
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
  shift
  [[ $# -gt 0 ]] || set -- benchmark uncorrected corrected
  for context in "$@"; do
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
[[ ! -e "${parent}/manifests/campaign_excluded_nodes.tsv" ]]
make_contexts excluded > "${TMP_DIR}/excluded.tsv"
export CAMPAIGN_BARRIER="${TMP_DIR}/excluded-barrier"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
  ECODA_RUN_ID=campaign_excluded CAMPAIGN_EXPECT_EXCLUDED_NODES=cpu046,cpu051,cpu052 \
  bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
  --contexts "${TMP_DIR}/excluded.tsv" --derivative-release test-release \
  --exclude-nodes cpu046,cpu051,cpu052 > "${TMP_DIR}/excluded.log"
excluded_parent="${HPC_ROOT}/_ecoda_runs/campaign_excluded"
[[ "$(cat "${excluded_parent}/manifests/campaign_excluded_nodes.tsv")" == cpu046,cpu051,cpu052 ]]
grep -q 'campaign_excluded_nodes.tsv' "${excluded_parent}/manifests/campaign_identities.tsv"
grep -q '^STATE=OK$' "${excluded_parent}/status/terminal"
for bad in '' 'cpu046,cpu046' 'cpu[046-052]' 'cpu046, cpu051' '--help' 'cpu046,'; do
  if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
      ECODA_RUN_ID=campaign_invalid_exclusion bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
      --contexts "${TMP_DIR}/excluded.tsv" --derivative-release test-release --exclude-nodes "${bad}" >/dev/null 2>&1; then
    echo 'campaign accepted an unsafe or duplicate excluded-node list' >&2; exit 1
  fi
done
[[ ! -e "${HPC_ROOT}/_ecoda_runs/campaign_invalid_exclusion" ]]
export CAMPAIGN_BARRIER="${TMP_DIR}/barrier"
for spec in 'pair benchmark corrected' 'single corrected'; do
  set -- ${spec}; prefix="$1"; shift; count=$#
  policy=auto; policy_args=()
  make_contexts "${prefix}" "$@" > "${TMP_DIR}/${prefix}.tsv"
  if [[ "${prefix}" == pair ]]; then
    printf 'Stephenson\tbatch_effect_corrected\tqot\n' > "${TMP_DIR}/corrected.matrix.tsv"
    awk -v matrix="${TMP_DIR}/corrected.matrix.tsv" 'BEGIN{FS=OFS="\t"} $1=="corrected" {$5=matrix} {print}' \
      "${TMP_DIR}/${prefix}.tsv" > "${TMP_DIR}/${prefix}.matrix.tsv"
    cp "${TMP_DIR}/${prefix}.matrix.tsv" "${TMP_DIR}/${prefix}.tsv"
    policy=any; policy_args=(--batch-gpu-policy any)
  fi
  export CAMPAIGN_BARRIER="${TMP_DIR}/${prefix}-barrier"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
    ECODA_RUN_ID="campaign_${prefix}" CAMPAIGN_EXPECT_CHILDREN="${count}" CAMPAIGN_EXPECT_BATCH_GPU_POLICY="${policy}" \
    bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
      --contexts "${TMP_DIR}/${prefix}.tsv" --derivative-release test-release ${policy_args[@]+"${policy_args[@]}"} > "${TMP_DIR}/${prefix}.log"
  [[ $(wc -l < "${HPC_ROOT}/_ecoda_runs/campaign_${prefix}/manifests/campaign_scheduler_ids.tsv" | tr -d ' ') == "${count}" ]]
  [[ ! -e "${HPC_ROOT}/_ecoda_runs/${prefix}_uncorrected" ]]
done
cmp "${TMP_DIR}/corrected.matrix.tsv" "${HPC_ROOT}/_ecoda_runs/campaign_pair/manifests/campaign_corrected.method_matrix.tsv"
[[ -f "${TMP_DIR}/pair-barrier/corrected.gpu_policy" && ! -e "${TMP_DIR}/pair-barrier/benchmark.gpu_policy" ]]
[[ $(wc -l < "${parent}/manifests/campaign_scheduler_ids.tsv" | tr -d ' ') == 3 ]]
cmp "${TMP_DIR}/benchmark.selection.tsv" "${parent}/manifests/campaign_benchmark.selection.tsv"
[[ $(find "${TMP_DIR}/barrier" -name '*.started' | wc -l | tr -d ' ') == 3 ]]
"${TEST_REAL_PYTHON}" - "${parent}" <<'PY'
import hashlib,pathlib,sys
root=pathlib.Path(sys.argv[1]); manifest=root/'manifests'
identities=[line.split('\t') for line in (manifest/'campaign_identities.tsv').read_text().splitlines()]
assert len(identities)==7
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
[[ $(wc -l < "${HPC_ROOT}/_ecoda_runs/campaign_failure/manifests/campaign_scheduler_ids.tsv" | tr -d ' ') == 6 ]]
for record in METADATA_EXPORT:201 PREFLIGHT:301 ARRAY:401; do
  grep -qx $'benchmark\tretry_benchmark\t'"${record%%:*}"$'\t'"${record#*:}" \
    "${HPC_ROOT}/_ecoda_runs/campaign_failure/manifests/campaign_scheduler_ids.tsv"
done

export CAMPAIGN_BARRIER="${TMP_DIR}/bad-proof-barrier"
make_contexts proof > "${TMP_DIR}/proof.tsv"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
    ECODA_RUN_ID=campaign_bad_proof CAMPAIGN_BAD_PROOF_CONTEXT=benchmark \
    bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
    --contexts "${TMP_DIR}/proof.tsv" --derivative-release test-release > "${TMP_DIR}/proof.log" 2>&1; then
  echo 'campaign accepted invalid scheduler proof' >&2; exit 1
fi
grep -q 'ERROR: campaign scheduler manifest is missing or checksum-invalid: benchmark:proof_benchmark:' "${TMP_DIR}/proof.log"
grep -q '^STATE=FAIL$' "${HPC_ROOT}/_ecoda_runs/campaign_bad_proof/status/terminal"
[[ $(wc -l < "${HPC_ROOT}/_ecoda_runs/campaign_bad_proof/manifests/campaign_scheduler_ids.tsv" | tr -d ' ') == 2 ]]

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
# Optional prepared receipt is frozen once and forwarded unchanged to every
# independent child; the ordinary four-column contexts remain unchanged.
make_contexts adopted > "${TMP_DIR}/adopted.tsv"
"${TEST_REAL_PYTHON}" - "${SOURCE_TREE}/datasets.json" "${TMP_DIR}/prepared.json" <<'PYADOPT'
import hashlib,json,pathlib,sys
config,path=map(pathlib.Path,sys.argv[1:])
path.write_text(json.dumps(dict(schema='ecoda_stage5_prepared_adoption_v1',state='APPROVED_COPY_ONLY',original_gate_state='FAILED',
 config_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),target_run_ids={'benchmark':'adopted_benchmark',
 'batch_effect/uncorrected_final':'adopted_uncorrected','batch_effect/corrected_final':'adopted_corrected'})))
PYADOPT
prepared_sha="$(shasum -a 256 "${TMP_DIR}/prepared.json" | awk '{print $1}')"
export CAMPAIGN_BARRIER="${TMP_DIR}/adopted-barrier"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
  ECODA_RUN_ID=campaign_adoption CAMPAIGN_EXPECT_ADOPTION=1 \
  bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
  --contexts "${TMP_DIR}/adopted.tsv" --derivative-release test-release \
  --prepared-adoption "${TMP_DIR}/prepared.json" --prepared-adoption-sha256 "${prepared_sha}" > "${TMP_DIR}/adopted.log"
[[ $(find "${CAMPAIGN_BARRIER}" -name '*.adoption' | wc -l | tr -d ' ') == 3 ]]
parent="${HPC_ROOT}/_ecoda_runs/campaign_adoption"
cmp "${TMP_DIR}/prepared.json" "${parent}/manifests/campaign_prepared_adoption.json"
grep -q '^STATE=OK$' "${parent}/status/terminal"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
    ECODA_RUN_ID=campaign_bad_adoption bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
    --contexts "${TMP_DIR}/adopted.tsv" --derivative-release test-release \
    --prepared-adoption "${TMP_DIR}/prepared.json" --prepared-adoption-sha256 "$(printf '0%.0s' {1..64})" > "${TMP_DIR}/bad-adoption.log" 2>&1; then
  echo 'campaign accepted altered prepared receipt' >&2; exit 1
fi
[[ ! -e "${HPC_ROOT}/_ecoda_runs/campaign_bad_adoption" ]]

# A selected metadata-only receipt is forwarded only to its target context.
make_contexts meta benchmark corrected > "${TMP_DIR}/metadata-contexts.tsv"
"${TEST_REAL_PYTHON}" - "${SOURCE_TREE}/datasets.json" "${TMP_DIR}/metadata-copy.json" <<'PY'
import hashlib,json,pathlib,sys
config,path=map(pathlib.Path,sys.argv[1:]);ctx='batch_effect/corrected_final'
path.write_text(json.dumps(dict(schema='ecoda_stage5_metadata_adoption_v1',state='APPROVED_COPY_ONLY',original_gate_state='FAILED',
 config_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),target_run_ids={ctx:'meta_corrected'},sources={ctx:{}},
 selected_metadata={ctx:[dict(dataset='Stephenson',view='batch_effect_corrected')]},
 artifacts=[dict(context=ctx,kind='metadata',dataset='Stephenson',view='batch_effect_corrected')])))
PY
metadata_sha="$(shasum -a 256 "${TMP_DIR}/metadata-copy.json" | awk '{print $1}')"
export CAMPAIGN_BARRIER="${TMP_DIR}/metadata-barrier"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
  ECODA_RUN_ID=campaign_metadata CAMPAIGN_EXPECT_CHILDREN=2 CAMPAIGN_EXPECT_ADOPTION=1 CAMPAIGN_EXPECT_ADOPTION_CONTEXTS=corrected \
  bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" --contexts "${TMP_DIR}/metadata-contexts.tsv" \
    --derivative-release test-release --prepared-adoption "${TMP_DIR}/metadata-copy.json" \
    --prepared-adoption-sha256 "${metadata_sha}" > "${TMP_DIR}/metadata.log"
[[ -f "${CAMPAIGN_BARRIER}/corrected.adoption" && ! -e "${CAMPAIGN_BARRIER}/benchmark.adoption" && ! -e "${HPC_ROOT}/_ecoda_runs/meta_uncorrected" ]]

make_contexts invalid benchmark corrected > "${TMP_DIR}/invalid-base.tsv"
for variant in empty unknown extra wrong_matrix duplicate_child matrix_duplicate matrix_foreign matrix_wrongview; do
  case "${variant}" in
    empty) : > "${TMP_DIR}/invalid-${variant}.tsv" ;;
    unknown) sed 's/^benchmark/foreign/' "${TMP_DIR}/invalid-base.tsv" > "${TMP_DIR}/invalid-${variant}.tsv" ;;
    extra) cat "${TMP_DIR}/invalid-base.tsv" "${TMP_DIR}/invalid-base.tsv" > "${TMP_DIR}/invalid-${variant}.tsv" ;;
    wrong_matrix) awk -v matrix="${TMP_DIR}/corrected.matrix.tsv" 'BEGIN{FS=OFS="\t"} {$5=matrix;print}' "${TMP_DIR}/invalid-base.tsv" > "${TMP_DIR}/invalid-${variant}.tsv" ;;
    duplicate_child) sed 's/invalid_corrected/invalid_benchmark/' "${TMP_DIR}/invalid-base.tsv" > "${TMP_DIR}/invalid-${variant}.tsv" ;;
    matrix_*)
      case "${variant}" in
        matrix_duplicate) cat "${TMP_DIR}/corrected.matrix.tsv" "${TMP_DIR}/corrected.matrix.tsv" > "${TMP_DIR}/invalid.matrix.tsv" ;;
        matrix_foreign) printf 'Diabetes\tbatch_effect_corrected\tqot\n' > "${TMP_DIR}/invalid.matrix.tsv" ;;
        matrix_wrongview) printf 'Stephenson\tbatch_effect_uncorrected\tqot\n' > "${TMP_DIR}/invalid.matrix.tsv" ;;
      esac
      awk -v matrix="${TMP_DIR}/invalid.matrix.tsv" 'BEGIN{FS=OFS="\t"} $1=="corrected" {$5=matrix} {print}' \
        "${TMP_DIR}/invalid-base.tsv" > "${TMP_DIR}/invalid-${variant}.tsv" ;;
  esac
  if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" NAS_TARGET_DIR="${NAS_BASE}" \
      ECODA_RUN_ID="campaign_invalid_${variant}" bash "${SOURCE_TREE}/src/5_run_benchmark_methods/stage5_campaign.sh" \
      --contexts "${TMP_DIR}/invalid-${variant}.tsv" --derivative-release test-release > "${TMP_DIR}/invalid-${variant}.log" 2>&1; then
    echo "invalid campaign accepted: ${variant}" >&2; exit 1
  fi
  [[ ! -e "${HPC_ROOT}/_ecoda_runs/invalid_benchmark" && ! -e "${HPC_ROOT}/_ecoda_runs/invalid_corrected" ]]
done
echo 'PASS: one pinned parent, concurrent isolated contexts, frozen exact scopes, child-failure isolation, collected scheduler IDs and prepared receipt forwarding'
