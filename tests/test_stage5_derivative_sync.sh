#!/bin/bash
# Publication boundaries only: tiny text payloads, stubbed NAS mount probe.
set -euo pipefail
umask 077
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-derivative-sync.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
export PROJECT_ROOT="${ROOT}" PYTHON_BIN="${ROOT}/.pixi/envs/default/bin/python"
export HPC_SCRATCH_DIR="${TMP}/scratch" NAS_TARGET_DIR="${TMP}/nas"
export DATASETS_JSON_FILE="${TMP}/datasets.json"
export ECODA_RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/derivative-sync" ECODA_RUN_ID=derivative-sync
mkdir -p "${ECODA_RUN_ROOT}/manifests" "${NAS_TARGET_DIR}/batch_effect/corrected_final" "${TMP}/canonical"
printf 'old canonical result\n' > "${NAS_TARGET_DIR}/batch_effect/corrected_final/old.rds"
mkdir -p "${NAS_TARGET_DIR}/releases/test-release/unselected"
printf 'unrelated metadata\n' > "${NAS_TARGET_DIR}/releases/test-release/unselected/metadata.json"
printf '{"Toy":{"columns":{"label":"biology"},"views":{"batch_effect_corrected":{"output_file_name":"toy.h5ad"}}}}\n' > "${DATASETS_JSON_FILE}"
printf 'source identity fixture\n' > "${ECODA_RUN_ROOT}/manifests/source.manifest"
printf 'runtime identity fixture\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"
export ECODA_SELECTION_MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'Toy\tbatch_effect_corrected\tmrvi\n' > "${ECODA_SELECTION_MANIFEST}"
"${PYTHON_BIN}" - "${TMP}" <<'PY'
import hashlib, pathlib, sys
root=pathlib.Path(sys.argv[1]); source=root/'canonical/toy.h5ad'; source.write_bytes(b'canonical immutable input')
md5=hashlib.md5(source.read_bytes()).hexdigest(); size=source.stat().st_size; proof=pathlib.Path(str(source)+'.md5')
proof.write_text(f'MD5={md5}\nSIZE={size}\nPATH={source}\n')
(root/'inputs.tsv').write_text(f'Toy\tbatch_effect_corrected\t{source}\t{md5}\t{size}\t{proof}\n')
PY
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/utils/bash/ecoda_derivative_artifacts.sh"
source "${ROOT}/src/utils/bash/ecoda_stage5_policy.sh"
ecoda_derivative_bind stage5 test-release "${ECODA_SELECTION_MANIFEST}" "${TMP}/inputs.tsv"
export ANALYSIS_PASS=corrected PASS_ARG=corrected ANALYSIS_VARIANT=corrected_final
export ANALYSIS_ROOT="$(ecoda_derivative_root batch_effect/corrected_final)"
export ANALYSIS_NAS_ROOT="$(ecoda_derivative_root batch_effect/corrected_final nas)"
export ANALYSIS_LOG_PREFIX=execution_times_batch_effect_corrected_final_
export EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs"
mkdir -p "${ANALYSIS_ROOT}/results" "${ANALYSIS_ROOT}/metadata" "${EXECUTION_LOG_DIR}"
printf 'derived scientific payload\n' > "${ANALYSIS_ROOT}/results/toy.rds"
printf 'Sample\nA\n' > "${ANALYSIS_ROOT}/metadata/Toy_sample_metadata.feather"
printf 'preserved worker timing\n' > "${EXECUTION_LOG_DIR}/execution_times_fixture.feather"
ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}" >/dev/null
ecoda_write_checksum "${ANALYSIS_ROOT}/results/toy.rds" >/dev/null
ecoda_write_checksum "${ANALYSIS_ROOT}/metadata/Toy_sample_metadata.feather" >/dev/null
printf 'Toy\tbatch_effect_corrected\t%s\t%s\n' "${TMP}/canonical/toy.h5ad" \
  "${ANALYSIS_ROOT}/metadata/Toy_sample_metadata.feather" > "${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
source "${ROOT}/src/5_run_benchmark_methods/benchmark_submit_common.sh"
printf '#!/bin/bash\nexit 0\n' > "${TMP}/nas-probe.sh"
benchmark_source_script_path() { printf '%s/nas-probe.sh' "${TMP}"; }
benchmark_require_source_script_path() { [[ -r "$1" ]]; }
benchmark_sync_artifacts_for() { SYNC_ARTIFACTS=("${ANALYSIS_ROOT}/results/toy.rds"); }
notify_sync_status() { return 0; }
benchmark_job_durations_block() { return 0; }
DATASET_NAMES=(Toy)
nas_root_mode="$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}")"
scratch_root_mode="$(_ecoda_artifact_mode_bits "${HPC_SCRATCH_DIR}")"
if [[ "${STAGE5_DERIVATIVE_SYNC_TRANSFER_FAULT:-0}" == 1 ]]; then
  rsync() { return 23; }
  if analysis_merge_sync_cleanup mrvi > "${TMP}/transfer-fault.log" 2>&1; then
    echo 'failed selected transfer reported sync success' >&2; exit 1
  fi
  grep -q 'selected benchmark rsync failed' "${TMP}/transfer-fault.log"
  [[ "$(ecoda_owner_state "$(ecoda_owner_dir stage5 "sync/${ANALYSIS_ROOT}")")" == FAIL ]]
  [[ ! -e "${ANALYSIS_NAS_ROOT}/results/toy.rds" && ! -e "${ANALYSIS_NAS_ROOT}/checksums.md5" ]]
  printf 'PASS: fresh derivative selected transfer failure fails closed\n'
  exit 0
fi
if [[ -n "${STAGE5_DERIVATIVE_SYNC_CHMOD_FAULT:-}" ]]; then
  fault_directory="${NAS_TARGET_DIR}/releases/test-release/derivative-sync"
  mkdir -p "${fault_directory}"
  chmod 700 "${fault_directory}"
  chmod() {
    if [[ "$1" == 777 && "${2:-}" == "${fault_directory}" ]]; then
      [[ "${STAGE5_DERIVATIVE_SYNC_CHMOD_FAULT}" == ineffective ]] && return 0
      return 1
    fi
    command chmod "$@"
  }
  if analysis_merge_sync_cleanup mrvi > "${TMP}/fault.log" 2>&1; then
    echo 'incorrect selected ancestor permissions reported sync success' >&2; exit 1
  fi
  grep -q 'selected NAS ancestor' "${TMP}/fault.log"
  [[ "$(_ecoda_artifact_mode_bits "${fault_directory}")" == 700 ]]
  [[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}")" == "${nas_root_mode}" ]]
  printf 'PASS: Stage5 selected ancestor chmod failure/ineffectiveness fails closed\n'
  exit 0
fi
analysis_merge_sync_cleanup mrvi
ancestor="${ANALYSIS_NAS_ROOT}"
while [[ "${ancestor}" != "${NAS_TARGET_DIR}" ]]; do
  benchmark_require_nas_mode "${ancestor}" 777
  ancestor="$(dirname "${ancestor}")"
done
[[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}")" == "${nas_root_mode}" ]]
[[ "$(_ecoda_artifact_mode_bits "${HPC_SCRATCH_DIR}")" == "${scratch_root_mode}" ]]
[[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}/releases/test-release/unselected")" == 700 ]]
[[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}/releases/test-release/unselected/metadata.json")" == 600 ]]
[[ "$(cat "${NAS_TARGET_DIR}/releases/test-release/unselected/metadata.json")" == 'unrelated metadata' ]]
cmp "${ANALYSIS_ROOT}/results/toy.rds" "${ANALYSIS_NAS_ROOT}/results/toy.rds"
source "${ROOT}/src/utils/bash/ecoda_run_audit_stage5.sh"
_audit_die() { echo "$*" >&2; return 1; }
for path in "${ANALYSIS_NAS_ROOT}/results/toy.rds" \
    "${ANALYSIS_NAS_ROOT}/metadata/Toy_sample_metadata.feather"; do
  _audit_derivative_stage5_nas_artifact "${path}"
  printf 'changed mirror\n' >> "${path}"
  if _audit_derivative_stage5_nas_artifact "${path}" > "${TMP}/bad-nas.log" 2>&1; then
    echo 'terminal audit accepted changed NAS bytes' >&2; exit 1
  fi
  # A newly consistent NAS checksum still cannot substitute for scratch equality.
  ecoda_write_checksum "${path}"; chmod 666 "${path}.md5"
  if _audit_derivative_stage5_nas_artifact "${path}" > "${TMP}/bad-nas-own-checksum.log" 2>&1; then
    echo 'terminal audit accepted independently checksummed different NAS bytes' >&2; exit 1
  fi
  cp "${ANALYSIS_ROOT}/${path#${ANALYSIS_NAS_ROOT}/}" "${path}"
  cp "${ANALYSIS_ROOT}/${path#${ANALYSIS_NAS_ROOT}/}.md5" "${path}.md5"
  chmod 666 "${path}" "${path}.md5"
  _audit_derivative_stage5_nas_artifact "${path}"
  for target in "${path}" "${path}.md5"; do
    chmod 640 "${target}"
    if _audit_derivative_stage5_nas_artifact "${path}" > "${TMP}/bad-nas-mode.log" 2>&1; then
      echo 'terminal audit accepted incorrect NAS file mode' >&2; exit 1
    fi
    chmod 666 "${target}"
  done
  directory="$(dirname "${path}")"
  chmod 750 "${directory}"
  if _audit_derivative_stage5_nas_artifact "${path}" > "${TMP}/bad-nas-dir.log" 2>&1; then
    echo 'terminal audit accepted incorrect selected NAS directory mode' >&2; exit 1
  fi
  chmod 777 "${directory}"
done
# An abort after an actual validated publication can preserve that NAS success;
# a method watchdog alone cannot manufacture the missing mirror/publication.
for function_name in stage5_finalize_tracked_owner_states stage5_derivative_nas_publication_valid; do
  eval "$(awk -v wanted="${function_name}()" '$0==wanted" {" {copy=1} copy {print} copy && /^}$/ {exit}' \
    "${ROOT}/src/5_run_benchmark_methods/stage5_dispatcher.sh")"
done
stage5_watchdog_state_for_artifact_path() { STAGE5_WATCHDOG_STATE=OK; }
RUN_ID="${ECODA_RUN_ID}"
ecoda_owner_clear_tracked
for path in "${ANALYSIS_ROOT}/results/toy.rds" "${ANALYSIS_NAS_ROOT}/results/toy.rds"; do
  ecoda_artifact_owner_acquire "${path}" stage5 "${RUN_ID}" 1 0 0 >/dev/null
  ecoda_owner_track "${ECODA_ARTIFACT_OWNER_DIR}"
done
stage5_finalize_tracked_owner_states 'fixture abort after successful synchronization' > "${TMP}/published-finalization.log"
ecoda_artifact_owner_validate "${ANALYSIS_NAS_ROOT}/results/toy.rds" "${RUN_ID}"
[[ "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]]
[[ "$(cat "${NAS_TARGET_DIR}/batch_effect/corrected_final/old.rds")" == 'old canonical result' ]]
[[ -s "${EXECUTION_LOG_DIR}/execution_times_fixture.feather" ]]
# A repeated identical publication validates; a different existing derivative
# payload fails closed and retains its bytes, without replacement receipts.
analysis_merge_sync_cleanup mrvi
# Identical publication validates ancestor modes without repairing old evidence.
ancestor="${NAS_TARGET_DIR}/releases/test-release/derivative-sync"
chmod 700 "${ancestor}"
if analysis_merge_sync_cleanup mrvi > "${TMP}/identical-bad-ancestor.log" 2>&1; then
  echo 'identical reuse accepted an untraversable ancestor' >&2; exit 1
fi
[[ "$(_ecoda_artifact_mode_bits "${ancestor}")" == 700 ]]
chmod 777 "${ancestor}"
chmod u+w "${ANALYSIS_NAS_ROOT}/results/toy.rds"
printf 'external conflicting bytes\n' > "${ANALYSIS_NAS_ROOT}/results/toy.rds"
if analysis_merge_sync_cleanup mrvi > "${TMP}/conflict.log" 2>&1; then
  echo 'existing derivative payload was overwritten' >&2; exit 1
fi
[[ "$(cat "${ANALYSIS_NAS_ROOT}/results/toy.rds")" == 'external conflicting bytes' ]]
printf 'PASS: derivative NAS publication, identical reuse, canonical preservation and conflict rejection\n'
