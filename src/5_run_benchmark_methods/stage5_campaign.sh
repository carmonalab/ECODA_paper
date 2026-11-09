#!/bin/bash
# One durable-gate wrapper; one to three selected independent dispatchers.
# Run once through the immutable snapshot executor.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"

CONTEXTS=""; RELEASE=""; PREPARED_ADOPTION=""; PREPARED_ADOPTION_SHA=""; BATCH_GPU_POLICY=auto; BATCH_GPU_POLICY_SET=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --contexts) CONTEXTS="${2:-}"; shift 2 ;;
    --derivative-release) RELEASE="${2:-}"; shift 2 ;;
    --prepared-adoption) PREPARED_ADOPTION="${2:?receipt required}"; shift 2 ;;
    --prepared-adoption-sha256) PREPARED_ADOPTION_SHA="${2:?receipt SHA256 required}"; shift 2 ;;
    --batch-gpu-policy) [[ ${BATCH_GPU_POLICY_SET} -eq 0 ]] || exit 1; BATCH_GPU_POLICY_SET=1; BATCH_GPU_POLICY="${2:-}"; shift 2 ;;
    *) echo 'Usage: stage5_campaign.sh --contexts TSV --derivative-release ID [--batch-gpu-policy auto|default|any] [--prepared-adoption JSON --prepared-adoption-sha256 SHA256]' >&2; exit 1 ;;
  esac
done
case "${BATCH_GPU_POLICY}" in auto|default|any) ;; *) echo 'ERROR: invalid batch GPU policy' >&2; exit 1 ;; esac
[[ -n "${CONTEXTS}" && "${RELEASE}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
   "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 ]] || {
  echo 'ERROR: campaign requires an explicit scope, release and immutable executor' >&2; exit 1;
}
if [[ -n "${PREPARED_ADOPTION}" || -n "${PREPARED_ADOPTION_SHA}" ]]; then
  [[ "${PREPARED_ADOPTION}" == /* && ! -L "${PREPARED_ADOPTION}" &&
     "${PREPARED_ADOPTION_SHA}" =~ ^[0-9a-f]{64}$ &&
     "$(ecoda_realpath_existing "${PREPARED_ADOPTION}")" == "${PREPARED_ADOPTION}" &&
     "$(ecoda_sha256_file "${PREPARED_ADOPTION}")" == "${PREPARED_ADOPTION_SHA}" ]] || {
    echo 'ERROR: campaign prepared-adoption receipt is missing, unsafe or checksum-invalid' >&2; exit 1;
  }
fi
PARENT_ID="${ECODA_RUN_ID:?immutable executor must supply a parent run ID}"
ecoda_validate_run_id "${PARENT_ID}"
ecoda_require_source_script_path "${SCRIPT_DIR}/stage5_campaign.sh" "${ECODA_SOURCE_ROOT}"
ecoda_require_source_script_path "${SCRIPT_DIR}/stage5_dispatcher.sh" "${ECODA_SOURCE_ROOT}"
[[ -d "${NAS_TARGET_DIR}" && ! -L "${NAS_TARGET_DIR}" ]] || {
  echo 'ERROR: campaign controller requires its configured NAS mount' >&2; exit 1;
}
ecoda_init_run stage5 "${PARENT_ID}" >/dev/null
PARENT_ROOT="${ECODA_RUN_ROOT}"
campaign_failed() {
  ecoda_set_run_state FAIL 'campaign did not complete every selected context' || true
}
trap campaign_failed EXIT
ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${PARENT_ROOT}/manifests/source.manifest"
export ECODA_RUNTIME_PROFILE=stage5
ecoda_runtime_validate_submission apptainer
ecoda_atomic_copy "${CONTEXTS}" "${PARENT_ROOT}/manifests/campaign_requested.tsv"
if [[ -n "${PREPARED_ADOPTION}" ]]; then
  ecoda_atomic_copy "${PREPARED_ADOPTION}" "${PARENT_ROOT}/manifests/campaign_prepared_adoption.json"
  PREPARED_ADOPTION="${PARENT_ROOT}/manifests/campaign_prepared_adoption.json"
  [[ "$(ecoda_sha256_file "${PREPARED_ADOPTION}")" == "${PREPARED_ADOPTION_SHA}" ]] || exit 1
fi

# Freeze every scientific request before starting any dispatcher. The metadata
# check reads only sizes/proofs; each existing child performs its full preflight.
"${PYTHON_BIN}" - "${PARENT_ROOT}" "${PROJECT_ROOT}" "${ECODA_RUNS_ROOT}" "${BATCH_GPU_POLICY}" <<'PY'
import hashlib, json, os, pathlib, sys
parent, project, runs = map(pathlib.Path, sys.argv[1:4])
sys.path.insert(0, str(project / 'src/utils/py'))
from benchmark_request_scope import read_request, normalized_rows
from derivative_artifacts import absolute, inputs
if not absolute(os.environ['NAS_TARGET_DIR']).is_dir():
    raise SystemExit('ERROR: campaign controller NAS base is not a mounted directory')
config = json.loads((project / 'datasets.json').read_text())
manifest = parent / 'manifests'
rows = [line.split('\t') for line in (manifest / 'campaign_requested.tsv').read_text().splitlines()]
expected = {'benchmark': 'benchmark_analysis', 'uncorrected': 'batch_effect_uncorrected', 'corrected': 'batch_effect_corrected'}
if not 1 <= len(rows) <= 3 or any(len(row) not in (4, 5) or any(not value for value in row) for row in rows):
    raise SystemExit('ERROR: campaign requires one to three contexts, with an optional batch method scope')
if (len({row[0] for row in rows}) != len(rows) or not {row[0] for row in rows} <= set(expected)
        or len({row[1] for row in rows}) != len(rows)):
    raise SystemExit('ERROR: campaign contexts and child run IDs must be unique')
if any(len(row) == 5 and row[0] == 'benchmark' for row in rows):
    raise SystemExit('ERROR: optional method scopes are batch-context only')
adoption = manifest / 'campaign_prepared_adoption.json'
receipt = None
context_keys = dict(benchmark='benchmark', uncorrected='batch_effect/uncorrected_final', corrected='batch_effect/corrected_final')
if adoption.exists():
    receipt = json.loads(adoption.read_text())
    targets = {context_keys[row[0]]: row[1] for row in rows}
    if (receipt.get('state') != 'APPROVED_COPY_ONLY'
            or (receipt.get('schema') != 'ecoda_stage5_batch_prepared_adoption_v1' and receipt.get('original_gate_state') != 'FAILED')
            or receipt.get('config_sha256') != hashlib.sha256((project / 'datasets.json').read_bytes()).hexdigest()):
        raise SystemExit('ERROR: prepared-adoption receipt differs from the exact campaign scope')
    if receipt.get('schema') == 'ecoda_stage5_prepared_adoption_v1':
        if len(rows) != 3 or receipt.get('target_run_ids') != targets:
            raise SystemExit('ERROR: original prepared adoption requires all three exact contexts')
    elif receipt.get('schema') == 'ecoda_stage5_batch_prepared_adoption_v1':
        from stage5_prepared_adoption import batch_review
        batch_review(receipt, {key for key, value in config.items() if not key.startswith('_') and value.get('use_for_batch_effect') is True})
        if any(targets.get(context) != child for context, child in receipt['target_run_ids'].items()):
            raise SystemExit('ERROR: batch-prepared adoption differs from exact target batch contexts')
    elif receipt.get('schema') == 'ecoda_stage5_metadata_adoption_v1':
        copies = receipt.get('target_run_ids', {})
        if (not copies or not set(copies) <= {'batch_effect/uncorrected_final', 'batch_effect/corrected_final'}
                or set(receipt.get('sources', {})) != set(copies)
                or set(receipt.get('selected_metadata', {})) != set(copies)
                or any(targets.get(context) != child for context, child in copies.items())):
            raise SystemExit('ERROR: metadata-only adoption differs from selected batch contexts')
    else:
        raise SystemExit('ERROR: unsupported campaign adoption schema')
staged, identities = [], []
policy = manifest / 'campaign_batch_gpu_policy.tsv'
with policy.open('x') as handle:
    handle.write(sys.argv[4] + '\n')
policy.chmod(0o444)
identities.append(f'{policy}\t{hashlib.sha256(policy.read_bytes()).hexdigest()}\n')
metadata_union = set()
for row in rows:
    context, run_id, selection_text, references_text = row[:4]
    import re
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', run_id) or run_id == parent.name:
        raise SystemExit('ERROR: invalid or parent-equal child run ID')
    if (runs / run_id).exists() or (runs / run_id).is_symlink():
        raise SystemExit('ERROR: child run root already exists; preserve it and select a new run ID')
    destinations = []
    sources = [(selection_text, 'selection.tsv'), (references_text, 'inputs.tsv')]
    if len(row) == 5:
        sources.append((row[4], 'method_matrix.tsv'))
    for source_text, suffix in sources:
        source = absolute(source_text)
        target = manifest / f'campaign_{context}.{suffix}'
        with target.open('xb') as handle:
            handle.write(source.read_bytes())
        target.chmod(0o444)
        identities.append(f'{target}\t{hashlib.sha256(target.read_bytes()).hexdigest()}\n')
        destinations.append(str(target))
    selection, references = map(pathlib.Path, destinations[:2])
    lines = [line.split('\t') for line in selection.read_text().splitlines()]
    if not lines or any(len(row) not in (3, 4) or row[1] != expected[context] for row in lines):
        raise SystemExit('ERROR: selection is outside its campaign context')
    if context == 'benchmark':
        declared = normalized_rows(read_request(selection, config))
    else:
        if any(len(row) != 3 or row[2] != expected[context] for row in lines) or len({tuple(row) for row in lines}) != len(lines):
            raise SystemExit('ERROR: batch campaign rows must name unique dataset/view/pass triples')
        declared = lines
        if len(row) == 5:
            matrix = [tuple(line.split('\t')) for line in pathlib.Path(destinations[2]).read_text().splitlines()]
            methods = {'prepare_pseudobulk', 'pseudobulk', 'gloscope', 'composition', 'mrvi', 'pilot', 'qot'}
            if (not matrix or len(set(matrix)) != len(matrix)
                    or any(len(item) != 3 or not all(item) or item[1] != expected[context] or item[2] not in methods for item in matrix)
                    or {item[:2] for item in matrix} != {tuple(item[:2]) for item in declared}):
                raise SystemExit('ERROR: batch method scope differs from exact selected dataset/view scope')
            if context == 'uncorrected':
                selected_methods = list(dict.fromkeys(item[2] for item in matrix))
                rectangle = [(item[0], item[1], method) for method in selected_methods for item in declared]
                if (matrix != rectangle or pathlib.Path(destinations[2]).read_bytes() != ''.join('\t'.join(item) + '\n' for item in rectangle).encode()
                        or any(config.get(item[0], {}).get('use_for_batch_effect') is not True or item[0].startswith('_') for item in declared)
                        or (set(selected_methods) & {'pseudobulk', 'composition'} and 'prepare_pseudobulk' not in selected_methods)):
                    raise SystemExit('ERROR: uncorrected scope must be an exact ordered rectangle with explicit preparation dependencies')
    referenced = inputs(references, config, 'stage5')
    if {key[:2] for key in referenced} != {tuple(row[:2]) for row in declared}:
        raise SystemExit('ERROR: input references differ from selected dataset/views')
    if receipt and receipt['schema'] in ('ecoda_stage5_metadata_adoption_v1', 'ecoda_stage5_batch_prepared_adoption_v1') and context_keys[context] in receipt['target_run_ids']:
        requested = (receipt['selected_metadata'][context_keys[context]] if receipt['schema'] == 'ecoda_stage5_metadata_adoption_v1' else
                     [{key: item[key] for key in ('dataset', 'view')} for item in receipt['artifacts'] if item['context'] == context_keys[context] and item['kind'] == 'metadata'])
        pairs = [(item.get('dataset'), item.get('view')) for item in requested if isinstance(item, dict) and set(item) == {'dataset', 'view'}]
        if len(pairs) != len(requested) or len(set(pairs)) != len(pairs) or set(pairs) != {tuple(item[:2]) for item in declared}:
            raise SystemExit('ERROR: metadata-only adoption must cover the exact selected dataset/views')
        metadata_union.update((context_keys[context], dataset, view) for dataset, view in pairs)
    staged.append('\t'.join((context, run_id, *destinations)) + '\n')
if receipt and receipt['schema'] == 'ecoda_stage5_metadata_adoption_v1':
    artifacts = receipt.get('artifacts', [])
    actual = [(item.get('context'), item.get('dataset'), item.get('view')) for item in artifacts if item.get('kind') == 'metadata']
    if len(actual) != len(artifacts) or len(set(actual)) != len(actual) or set(actual) != metadata_union:
        raise SystemExit('ERROR: metadata-only adoption has missing, duplicate or foreign payloads')
for name, content in (('campaign_contexts.tsv', ''.join(staged)), ('campaign_identities.tsv', ''.join(identities))):
    target = manifest / name
    with target.open('x') as handle:
        handle.write(content)
    target.chmod(0o444)
PY

for path in "${PARENT_ROOT}"/manifests/campaign_*; do
  ecoda_write_checksum "${path}"
  ecoda_write_artifact_record "${path}" stage5_campaign_scope "${PARENT_ID}" >/dev/null
done
for path in "${PARENT_ROOT}/manifests/source.manifest" "${PARENT_ROOT}/manifests/runtime.identity"; do
  ecoda_write_checksum "${path}"
  ecoda_write_artifact_record "${path}" stage5_campaign_identity "${PARENT_ID}" >/dev/null
done
while IFS=$'\t' read -r path sha; do
  [[ "$(ecoda_sha256_file "${path}")" == "${sha}" ]] || {
    echo 'ERROR: frozen campaign identity changed before launch' >&2; exit 1;
  }
done < "${PARENT_ROOT}/manifests/campaign_identities.tsv"
ecoda_atomic_write "${PARENT_ROOT}/metadata" \
  "STAGE=stage5\nRUN_ID=${PARENT_ID}\nSTATE=ACTIVE\nCAMPAIGN_RELEASE=${RELEASE}\nCAMPAIGN_CONTEXTS=${PARENT_ROOT}/manifests/campaign_contexts.tsv\nBATCH_GPU_POLICY=${BATCH_GPU_POLICY}\nBATCH_GPU_POLICY_MANIFEST=${PARENT_ROOT}/manifests/campaign_batch_gpu_policy.tsv\nSOURCE_MANIFEST=${PARENT_ROOT}/manifests/source.manifest\nRUNTIME_IDENTITY=${PARENT_ROOT}/manifests/runtime.identity\n"

CHILD_PIDS=(); CHILD_IDS=(); CHILD_CONTEXTS=()
campaign_interrupt() {
  local pid
  trap - INT TERM
  for pid in "${CHILD_PIDS[@]:-}"; do [[ -n "${pid}" ]] && kill -TERM "${pid}" 2>/dev/null || true; done
  for pid in "${CHILD_PIDS[@]:-}"; do [[ -n "${pid}" ]] && wait "${pid}" 2>/dev/null || true; done
  exit 1
}
trap campaign_interrupt INT TERM
while IFS=$'\t' read -r context child selection references method_matrix; do
  (
    unset ECODA_RUN_ROOT ECODA_ARTIFACT_BINDING ECODA_SELECTION_MANIFEST ECODA_STAGE5_METHOD_MATRIX ECODA_STAGE5_METHOD_MATRIX_MODE
    unset ECODA_STAGE5_RECOVERY_MANIFEST ECODA_PRECOMPUTE_BACKUP_RECEIPT ECODA_EXPECTED_PENDING
    unset ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_PASS ANALYSIS_VARIANT BENCHMARK_MANIFEST
    unset ANALYSIS_MANIFEST EXECUTION_LOG_DIR FORCE_BENCHMARK
    export ECODA_RUN_ID="${child}" ECODA_LOGS_DIR="${PARENT_ROOT}/logs/${child}"
    args=(--selection-file "${selection}" --derivative-release "${RELEASE}" --input-artifacts "${references}")
    adoption_context=benchmark
    case "${context}" in uncorrected) adoption_context=batch_effect/uncorrected_final ;; corrected) adoption_context=batch_effect/corrected_final ;; esac
    if [[ -n "${PREPARED_ADOPTION}" ]] &&
       jq -e --arg c "${adoption_context}" '.target_run_ids | has($c)' "${PREPARED_ADOPTION}" >/dev/null; then
      args+=(--prepared-adoption "${PREPARED_ADOPTION}" --prepared-adoption-sha256 "${PREPARED_ADOPTION_SHA}")
    fi
    [[ "${context}" == benchmark || "${BATCH_GPU_POLICY}" == auto ]] || args+=(--gpu-policy "${BATCH_GPU_POLICY}")
    case "${context}" in
      uncorrected)
        args+=(--pass uncorrected --analysis-variant final)
        if [[ -n "${method_matrix}" ]]; then
          ecoda_validate_checksum "${method_matrix}" || exit 1
          target_methods="$(awk -F '\t' '!seen[$3]++ {if(n++) printf ","; printf "%s", $3} END {print ""}' "${method_matrix}")"
          args+=(--target-methods "${target_methods}")
        else
          args+=(--methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot)
        fi ;;
      corrected) args+=(--pass corrected --analysis-variant corrected_final --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot) ;;
    esac
    [[ -z "${method_matrix}" || "${context}" != corrected ]] || args+=(--method-matrix "${method_matrix}")
    exec bash "${SCRIPT_DIR}/stage5_dispatcher.sh" "${args[@]}"
  ) > "${PARENT_ROOT}/logs/${context}.log" 2>&1 &
  CHILD_PIDS+=("$!"); CHILD_IDS+=("${child}"); CHILD_CONTEXTS+=("${context}")
  printf 'CAMPAIGN_CHILD_RUN_ID=%s:%s\n' "${context}" "${child}"
done < "${PARENT_ROOT}/manifests/campaign_contexts.tsv"

failed=0
: > "${PARENT_ROOT}/manifests/campaign_scheduler_ids.build.$$"
for idx in "${!CHILD_PIDS[@]}"; do
  child="${CHILD_IDS[${idx}]}"; context="${CHILD_CONTEXTS[${idx}]}"
  child_rc=0; wait "${CHILD_PIDS[${idx}]}" || child_rc=$?
  CHILD_PIDS[${idx}]=""
  child_root="${ECODA_RUNS_ROOT}/${child}"
  terminal="${child_root}/status/terminal"
  if [[ ! -f "${terminal}" || -L "${terminal}" ]] ||
     [[ "$(sed -n 's/^RUN_ID=//p' "${terminal}")" != "${child}" ]] ||
     [[ "$(grep -c '^STATE=' "${terminal}")" != 1 ]] ||
     ! grep -Eq '^STATE=(OK|NOOP_VALIDATED)$' "${terminal}"; then
    child_rc=1
  fi
  scheduler="${child_root}/manifests/scheduler_ids.tsv"
  if [[ -f "${scheduler}" ]] && ecoda_validate_checksum "${scheduler}"; then
    while IFS=$'\t' read -r kind scheduler_id extra; do
      [[ "${kind}" =~ ^(ARRAY|WATCHDOG|STATUS|AGGREGATE_GATE|PREFLIGHT|METADATA_EXPORT)$ &&
         "${scheduler_id}" =~ ^[0-9]+$ && -z "${extra}" ]] || {
        child_rc=1; continue;
      }
      printf '%s\t%s\t%s\t%s\n' "${context}" "${child}" "${kind}" "${scheduler_id}" \
        >> "${PARENT_ROOT}/manifests/campaign_scheduler_ids.build.$$"
    done < "${scheduler}"
  else
    echo "ERROR: campaign scheduler manifest is missing or checksum-invalid: ${context}:${child}:${scheduler}; preserve child records/logs for terminal accounting" >&2
    child_rc=1
  fi
  [[ ${child_rc} -eq 0 ]] || failed=1
  ecoda_atomic_write "${PARENT_ROOT}/status/${context}" \
    "CONTEXT=${context}\nRUN_ID=${child}\nEXIT_CODE=${child_rc}\nSCHEDULER_MANIFEST=${ECODA_RUNS_ROOT}/${child}/manifests/scheduler_ids.tsv\n"
  cat "${PARENT_ROOT}/logs/${context}.log"
done
ecoda_atomic_install_manifest "${PARENT_ROOT}/manifests/campaign_scheduler_ids.build.$$" \
  "${PARENT_ROOT}/manifests/campaign_scheduler_ids.tsv" 4
rm -f "${PARENT_ROOT}/manifests/campaign_scheduler_ids.build.$$"
ecoda_write_checksum "${PARENT_ROOT}/manifests/campaign_scheduler_ids.tsv"
if [[ ${failed} -ne 0 ]]; then exit 1; fi
ecoda_set_run_state OK 'all selected campaign contexts completed'
trap - EXIT
printf 'CAMPAIGN_RUN_ID=%s\n' "${PARENT_ID}"
