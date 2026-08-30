#!/bin/bash
#SBATCH --job-name=annotation_merge
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage4 \
  "${SCRIPT_DIR}/3.2_merge_worker.sh" || exit 1
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"

MANIFEST="${ANNOTATION_MERGE_MANIFEST:-}"
[[ -r "${MANIFEST}" ]] || { echo "ERROR: merge manifest is unreadable: ${MANIFEST}" >&2; exit 1; }
RUN_ID="${ANNOTATION_RUN_ID:-}"
ecoda_validate_run_id "${RUN_ID}" || { echo "ERROR: ANNOTATION_RUN_ID is invalid or missing." >&2; exit 1; }
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
ecoda_validate_run_owned_path "${MANIFEST}" "${RUN_ROOT}" ||
  { echo "ERROR: merge manifest is outside the Stage 4 run root." >&2; exit 1; }
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
[[ -n "${line}" ]] || { echo "ERROR: no merge row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME VIEWS MANIFEST_RUN_ROOT <<< "${line}"
[[ -n "${DS_NAME}" && -n "${VIEWS}" && -n "${MANIFEST_RUN_ROOT}" ]] || { echo "ERROR: malformed merge manifest row" >&2; exit 1; }
run_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
manifest_root_real="$(ecoda_realpath_existing "${MANIFEST_RUN_ROOT}" 2>/dev/null || true)"
[[ -n "${run_root_real}" && "${run_root_real}" == "${manifest_root_real}" ]] ||
  { echo "ERROR: merge manifest row run root is not canonical." >&2; exit 1; }
RUN_ROOT="${run_root_real}"
[[ -r "${RUN_ROOT}/metadata" ]] || { echo "ERROR: Stage 4 run metadata is missing." >&2; exit 1; }
metadata_stage="$(sed -n 's/^STAGE=//p' "${RUN_ROOT}/metadata" | head -1)"
metadata_run="$(sed -n 's/^RUN_ID=//p' "${RUN_ROOT}/metadata" | head -1)"
[[ "${metadata_stage}" == "stage4" && "${metadata_run}" == "${RUN_ID}" ]] ||
  { echo "ERROR: Stage 4 run metadata does not match merge run." >&2; exit 1; }
OWNER_DIR="$(ecoda_owner_dir stage4 "${DS_NAME}")"
owner_state="$(ecoda_owner_state "${OWNER_DIR}" 2>/dev/null || true)"
owner_run="$(ecoda_owner_run "${OWNER_DIR}" 2>/dev/null || true)"
owner_stage="$(ecoda_owner_field "${OWNER_DIR}" STAGE 2>/dev/null || true)"
owner_key="$(ecoda_owner_field "${OWNER_DIR}" KEY 2>/dev/null || true)"
[[ "${owner_run}" == "${RUN_ID}" && "${owner_stage}" == "stage4" &&
   "${owner_key}" == "${DS_NAME}" &&
   ( "${owner_state}" == "ACTIVE" || "${owner_state}" == "OK" ) ]] ||
  { echo "ERROR: Stage 4 dataset owner is missing or foreign." >&2; exit 1; }
ANNOT_DIR="${RUN_ROOT}/datasets/${DS_NAME}/annotations"
[[ -d "${ANNOT_DIR}" ]] || { echo "ERROR: annotation feather directory missing: ${ANNOT_DIR}" >&2; exit 1; }
IFS=',' read -r -a VIEW_LIST <<< "${VIEWS}"
[[ ${#VIEW_LIST[@]} -gt 0 ]] || { echo "ERROR: merge view list is empty" >&2; exit 1; }
UNION_PATH="${RUN_ROOT}/datasets/${DS_NAME}/union/union.h5ad"
ecoda_validate_checksum "${UNION_PATH}" || { echo "ERROR: annotation union checksum failed: ${UNION_PATH}" >&2; exit 1; }
SOURCE_PATHS=""
SOURCE_RECORDS=""
for view in "${VIEW_LIST[@]}"; do
  output_name="$(jq -r --arg ds "${DS_NAME}" --arg view "${view}" '.[$ds].views[$view].output_file_name // .[$ds].views[$view].output_file // empty' "${DATASETS_JSON_FILE}")"
  [[ -n "${output_name}" ]] || { echo "ERROR: missing output contract for ${DS_NAME}/${view}" >&2; exit 1; }
  h5ad_path="${HPC_SCRATCH_DIR}/${DS_NAME}/output/${output_name}"
  [[ -s "${h5ad_path}" ]] || { echo "ERROR: h5ad missing before merge: ${h5ad_path}" >&2; exit 1; }
  ecoda_validate_checksum "${h5ad_path}" || { echo "ERROR: h5ad checksum failed before merge: ${h5ad_path}" >&2; exit 1; }
  "${PYTHON_BIN}" "${SCRIPT_DIR}/3.1_merge_annotations.py" \
    --h5ad-path "${h5ad_path}" --annot-dir "${ANNOT_DIR}" --union-path "${UNION_PATH}" \
    --output-path "${h5ad_path}"
  "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
    --h5ad "${h5ad_path}" --require-sidecar >/dev/null
  ecoda_write_checksum "${h5ad_path}"
  source_record="${h5ad_path}|$(ecoda_md5_file "${h5ad_path}")|$(wc -c < "${h5ad_path}" | tr -d '[:space:]')"
  [[ -z "${SOURCE_PATHS}" ]] && SOURCE_PATHS="${h5ad_path}" ||
    SOURCE_PATHS="${SOURCE_PATHS};${h5ad_path}"
  [[ -z "${SOURCE_RECORDS}" ]] && SOURCE_RECORDS="${source_record}" ||
    SOURCE_RECORDS="${SOURCE_RECORDS};${source_record}"
done
union_md5="$(ecoda_md5_file "${UNION_PATH}")"
union_size="$(wc -c < "${UNION_PATH}" | tr -d '[:space:]')"
marker="${RUN_ROOT}/datasets/${DS_NAME}/merge.ok"
tmp="${marker}.tmp.$$"
printf 'STATE=OK\nDATASET=%s\nVIEWS=%s\nSOURCE_H5ADS=%s\nSOURCE_RECORDS=%s\nUNION_PATH=%s\nUNION_MD5=%s\nUNION_SIZE=%s\n' \
  "${DS_NAME}" "${VIEWS}" "${SOURCE_PATHS}" "${SOURCE_RECORDS}" \
  "${UNION_PATH}" "${union_md5}" "${union_size}" > "${tmp}"
mv -f "${tmp}" "${marker}"
