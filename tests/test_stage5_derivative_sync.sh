#!/bin/bash
# Publication boundaries only: tiny text payloads, stubbed NAS mount probe.
set -euo pipefail
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
analysis_merge_sync_cleanup mrvi
cmp "${ANALYSIS_ROOT}/results/toy.rds" "${ANALYSIS_NAS_ROOT}/results/toy.rds"
[[ "$(cat "${NAS_TARGET_DIR}/batch_effect/corrected_final/old.rds")" == 'old canonical result' ]]
[[ -s "${EXECUTION_LOG_DIR}/execution_times_fixture.feather" ]]
# A repeated identical publication validates; a different existing derivative
# payload fails closed and retains its bytes, without replacement receipts.
analysis_merge_sync_cleanup mrvi
chmod u+w "${ANALYSIS_NAS_ROOT}/results/toy.rds"
printf 'external conflicting bytes\n' > "${ANALYSIS_NAS_ROOT}/results/toy.rds"
if analysis_merge_sync_cleanup mrvi > "${TMP}/conflict.log" 2>&1; then
  echo 'existing derivative payload was overwritten' >&2; exit 1
fi
[[ "$(cat "${ANALYSIS_NAS_ROOT}/results/toy.rds")" == 'external conflicting bytes' ]]
printf 'PASS: derivative NAS publication, identical reuse, canonical preservation and conflict rejection\n'
