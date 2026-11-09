#!/bin/bash
# Real selected Stage3 synchronization; only data transfer/schema are stubbed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage3-nas-modes.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'command chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
umask 077
export HPC_SCRATCH_DIR="${TMP}/scratch" NAS_TARGET_DIR="${TMP}/nas"
export ECODA_RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/fixture" ECODA_RUN_ID=fixture
export ECODA_ARTIFACT_BINDING="${ECODA_RUN_ROOT}/manifests/artifact_bindings.json"
LOCAL_BASE="${ECODA_RUN_ROOT}/artifacts/datasets"
NAS_BASE="${NAS_TARGET_DIR}/releases/release-v1/fixture/datasets"
mkdir -p "${LOCAL_BASE}/Stephenson/output" "${NAS_TARGET_DIR}/releases/release-v1/unselected"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
eval "$(awk '/^sync_selected\(\)/ {copy=1} copy {print} copy && /^}$/ {exit}' "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh")"
VALIDATED_SYNC_REPORT_SET=0
SYNC_LOCK_ROOT="${TMP}/sync-control"
SOURCE="${LOCAL_BASE}/Stephenson/output/fixture.h5ad"
printf 'tiny payload for selected permission contract\n' > "${SOURCE}"
ecoda_write_checksum "${SOURCE}"
command chmod 444 "${SOURCE}" "${SOURCE}.md5"
printf 'keep unrelated metadata\n' > "${NAS_TARGET_DIR}/releases/release-v1/unselected/metadata.json"
command chmod 600 "${NAS_TARGET_DIR}/releases/release-v1/unselected/metadata.json"
printf 'Stephenson\tbenchmark_analysis\n' > "${TMP}/selection.tsv"
source_md5="$(ecoda_md5_file "${SOURCE}")"
source_mode="$(_ecoda_artifact_mode_bits "${SOURCE}")"
source_proof_md5="$(ecoda_md5_file "${SOURCE}.md5")"
root_mode="$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}")"
scratch_mode="$(_ecoda_artifact_mode_bits "${HPC_SCRATCH_DIR}")"
ecoda_derivative_root() { if [[ "${2:-}" == nas ]]; then printf '%s' "${NAS_BASE}"; else printf '%s' "${LOCAL_BASE}"; fi; }
output_path_for() { printf '%s' "${SOURCE}"; }
validate_h5ad() { [[ "$1" == Stephenson && "$2" == benchmark_analysis && -s "$3" ]]; }
# Keep the actual checksum comparison/receipt producer, using tiny local bytes.
rsync() {
  local files="" argument src dst row
  for argument in "$@"; do case "${argument}" in --files-from=*) files="${argument#*=}" ;; esac; done
  local args=("$@") n=$#
  src="${args[$((n-2))]}"; dst="${args[$((n-1))]}"
  while IFS= read -r row; do mkdir -p "$(dirname "${dst}/${row}")"; cp "${src}/${row}" "${dst}/${row}"; done < "${files}"
}
sync_selected "${TMP}/selection.tsv"
TARGET="${NAS_BASE}/Stephenson/output/fixture.h5ad"
[[ "$(ecoda_md5_file "${TARGET}")" == "${source_md5}" ]]
ecoda_validate_checksum "${TARGET}"
[[ "$(_ecoda_artifact_mode_bits "${TARGET}")" == 666 && "$(_ecoda_artifact_mode_bits "${TARGET}.md5")" == 666 ]]
directory="${TARGET%/*}"
while [[ "${directory}" != "${NAS_TARGET_DIR}" ]]; do
  [[ "$(_ecoda_artifact_mode_bits "${directory}")" == 777 ]]
  directory="${directory%/*}"
done
# No metadata, unselected run, NAS-base, scratch, or source permission changes.
[[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}")" == "${root_mode}" ]]
[[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}/releases/release-v1/unselected")" == 700 ]]
[[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}/releases/release-v1/unselected/metadata.json")" == 600 ]]
[[ "$(_ecoda_artifact_mode_bits "${HPC_SCRATCH_DIR}")" == "${scratch_mode}" ]]
[[ "$(_ecoda_artifact_mode_bits "${SOURCE}")" == "${source_mode}" && "$(ecoda_md5_file "${SOURCE}")" == "${source_md5}" ]]
[[ "$(ecoda_md5_file "${SOURCE}.md5")" == "${source_proof_md5}" ]]
# Fail closed even when chmod reports success but the NAS mode does not stick.
BAD_DIRECTORY="${NAS_TARGET_DIR}/releases/release-v1/fixture"
command chmod 700 "${BAD_DIRECTORY}"
chmod() { if [[ "$1" == 777 && "${2:-}" == "${BAD_DIRECTORY}" ]]; then return 0; fi; command chmod "$@"; }
if sync_selected "${TMP}/selection.tsv"; then echo 'untraversable selected run reported sync success' >&2; exit 1; fi
[[ "$(_ecoda_artifact_mode_bits "${BAD_DIRECTORY}")" == 700 ]]
[[ ! -e "${SYNC_LOCK_ROOT}/sync.lock" ]]
# A failed chmod likewise cannot report synchronization success.
chmod() { if [[ "$1" == 777 && "${2:-}" == "${BAD_DIRECTORY}" ]]; then return 1; fi; command chmod "$@"; }
if sync_selected "${TMP}/selection.tsv"; then echo 'failed selected chmod reported sync success' >&2; exit 1; fi
[[ ! -e "${SYNC_LOCK_ROOT}/sync.lock" ]]
printf 'PASS: selected Stage3 derivative NAS ancestors777/files666, bounded permissions and fail-closed mode verification\n'
