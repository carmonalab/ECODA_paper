#!/bin/bash
# One durable-gate wrapper; three independent existing Stage 5 dispatchers.
# Run once through the immutable snapshot executor, never through three executors.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"

CONTEXTS=""; RELEASE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --contexts) CONTEXTS="${2:-}"; shift 2 ;;
    --derivative-release) RELEASE="${2:-}"; shift 2 ;;
    *) echo 'Usage: stage5_campaign.sh --contexts TSV --derivative-release ID' >&2; exit 1 ;;
  esac
done
[[ -n "${CONTEXTS}" && "${RELEASE}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
   "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 ]] || {
  echo 'ERROR: campaign requires an explicit scope, release and immutable executor' >&2; exit 1;
}
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

# Freeze every scientific request before starting any dispatcher. The metadata
# check reads only sizes/proofs; each existing child performs its full preflight.
"${PYTHON_BIN}" - "${PARENT_ROOT}" "${PROJECT_ROOT}" "${ECODA_RUNS_ROOT}" <<'PY'
import hashlib, json, os, pathlib, sys
parent, project, runs = map(pathlib.Path, sys.argv[1:])
sys.path.insert(0, str(project / 'src/utils/py'))
from benchmark_request_scope import read_request, normalized_rows
from derivative_artifacts import absolute, inputs
if not absolute(os.environ['NAS_TARGET_DIR']).is_dir():
    raise SystemExit('ERROR: campaign controller NAS base is not a mounted directory')
config = json.loads((project / 'datasets.json').read_text())
manifest = parent / 'manifests'
rows = [line.split('\t') for line in (manifest / 'campaign_requested.tsv').read_text().splitlines()]
expected = {'benchmark': 'benchmark_analysis', 'uncorrected': 'batch_effect_uncorrected', 'corrected': 'batch_effect_corrected'}
if len(rows) != 3 or any(len(row) != 4 or any(not value for value in row) for row in rows):
    raise SystemExit('ERROR: campaign requires exactly three four-field context rows')
if {row[0] for row in rows} != set(expected) or len({row[1] for row in rows}) != 3:
    raise SystemExit('ERROR: campaign contexts and child run IDs must be unique')
staged, identities = [], []
for context, run_id, selection_text, references_text in rows:
    import re
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', run_id) or run_id == parent.name:
        raise SystemExit('ERROR: invalid or parent-equal child run ID')
    if (runs / run_id).exists() or (runs / run_id).is_symlink():
        raise SystemExit('ERROR: child run root already exists; preserve it and select a new run ID')
    destinations = []
    for source_text, suffix in ((selection_text, 'selection.tsv'), (references_text, 'inputs.tsv')):
        source = absolute(source_text)
        target = manifest / f'campaign_{context}.{suffix}'
        with target.open('xb') as handle:
            handle.write(source.read_bytes())
        target.chmod(0o444)
        identities.append(f'{target}\t{hashlib.sha256(target.read_bytes()).hexdigest()}\n')
        destinations.append(str(target))
    selection, references = map(pathlib.Path, destinations)
    lines = [line.split('\t') for line in selection.read_text().splitlines()]
    if not lines or any(len(row) not in (3, 4) or row[1] != expected[context] for row in lines):
        raise SystemExit('ERROR: selection is outside its campaign context')
    if context == 'benchmark':
        declared = normalized_rows(read_request(selection, config))
    else:
        if any(len(row) != 3 or row[2] != expected[context] for row in lines) or len({tuple(row) for row in lines}) != len(lines):
            raise SystemExit('ERROR: batch campaign rows must name unique dataset/view/pass triples')
        declared = lines
    referenced = inputs(references, config, 'stage5')
    if {key[:2] for key in referenced} != {tuple(row[:2]) for row in declared}:
        raise SystemExit('ERROR: input references differ from selected dataset/views')
    staged.append('\t'.join((context, run_id, *destinations)) + '\n')
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
  "STAGE=stage5\nRUN_ID=${PARENT_ID}\nSTATE=ACTIVE\nCAMPAIGN_RELEASE=${RELEASE}\nCAMPAIGN_CONTEXTS=${PARENT_ROOT}/manifests/campaign_contexts.tsv\nSOURCE_MANIFEST=${PARENT_ROOT}/manifests/source.manifest\nRUNTIME_IDENTITY=${PARENT_ROOT}/manifests/runtime.identity\n"

CHILD_PIDS=(); CHILD_IDS=(); CHILD_CONTEXTS=()
campaign_interrupt() {
  local pid
  trap - INT TERM
  for pid in "${CHILD_PIDS[@]:-}"; do [[ -n "${pid}" ]] && kill -TERM "${pid}" 2>/dev/null || true; done
  for pid in "${CHILD_PIDS[@]:-}"; do [[ -n "${pid}" ]] && wait "${pid}" 2>/dev/null || true; done
  exit 1
}
trap campaign_interrupt INT TERM
while IFS=$'\t' read -r context child selection references; do
  (
    unset ECODA_RUN_ROOT ECODA_ARTIFACT_BINDING ECODA_SELECTION_MANIFEST ECODA_STAGE5_METHOD_MATRIX
    unset ECODA_STAGE5_RECOVERY_MANIFEST ECODA_PRECOMPUTE_BACKUP_RECEIPT ECODA_EXPECTED_PENDING
    unset ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_PASS ANALYSIS_VARIANT BENCHMARK_MANIFEST
    unset ANALYSIS_MANIFEST EXECUTION_LOG_DIR FORCE_BENCHMARK
    export ECODA_RUN_ID="${child}" ECODA_LOGS_DIR="${PARENT_ROOT}/logs/${child}"
    args=(--selection-file "${selection}" --derivative-release "${RELEASE}" --input-artifacts "${references}")
    case "${context}" in
      uncorrected) args+=(--pass uncorrected --analysis-variant final --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot) ;;
      corrected) args+=(--pass corrected --analysis-variant corrected_final --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot) ;;
    esac
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
ecoda_set_run_state OK 'all three selected campaign contexts completed'
trap - EXIT
printf 'CAMPAIGN_RUN_ID=%s\n' "${PARENT_ID}"
