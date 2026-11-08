#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-bassez-stage2.XXXXXX")"
FIXTURE="$(cd "${FIXTURE}" && pwd -P)"
cleanup() { chmod -R u+w "${FIXTURE}" 2>/dev/null || true; rm -rf "${FIXTURE}"; }
trap cleanup EXIT
export HPC_SCRATCH_DIR="${FIXTURE}/scratch" NAS_TARGET_DIR="${FIXTURE}/nas"
export ECODA_SOURCE_ROOT="${FIXTURE}/tree" PROJECT_ROOT="${FIXTURE}/tree"
export DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
export REAL_PYTHON="${ROOT}/.pixi/envs/default/bin/python"
export CALLS="${FIXTURE}/python.calls" EVENTS="${FIXTURE}/ownership.calls" RSYNC_CALLS="${FIXTURE}/rsync.calls"
mkdir -p "${HPC_SCRATCH_DIR}" "${NAS_TARGET_DIR}" "${FIXTURE}/bin" \
  "${PROJECT_ROOT}/src/2_dataset_specific_preprocessing" "${PROJECT_ROOT}/src/utils/py" "${PROJECT_ROOT}/src/utils/bash"
cp "${ROOT}/datasets.json" "${DATASETS_JSON_FILE}"
for file in 1_submit_hpc.sh 1.submit.sh bassez_metadata_step.sh bassez_h5ad_metadata.py; do
  cp "${ROOT}/src/2_dataset_specific_preprocessing/${file}" "${PROJECT_ROOT}/src/2_dataset_specific_preprocessing/${file}"
done
cp "${ROOT}/src/utils/py/h5ad_source_identity.py" "${PROJECT_ROOT}/src/utils/py/"
for file in ecoda_run_common.sh ecoda_runtime.sh sync_status_email.sh; do
  cp "${ROOT}/src/utils/bash/${file}" "${PROJECT_ROOT}/src/utils/bash/${file}"
done
cat > "${PROJECT_ROOT}/src/slurm_config.sh" <<'CONFIG'
export PROJECT_ROOT="${ECODA_SOURCE_ROOT}"
export DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
export SLURM_PARTITION=fixture
export MAX_NUM_CHUNKS_PARALLEL=1
CONFIG
cat > "${FIXTURE}/bin/python" <<'PYTHON'
#!/bin/bash
printf '%s\n' "$*" >> "${CALLS}"
exec "${REAL_PYTHON}" "$@"
PYTHON
cat > "${FIXTURE}/bin/rsync" <<'RSYNC'
#!/bin/bash
set -euo pipefail
case "$*" in *'--chmod=ugo+rwX'*) ;; *) exit 90 ;; esac
case "$*" in *'--ignore-existing'*) ;; *) exit 91 ;; esac
args=("$@")
count=${#args[@]}
source="${args[$((count - 2))]}"
target="${args[$((count - 1))]}"
printf '%s\t%s\n' "${source}" "${target}" >> "${RSYNC_CALLS}"
[[ ! -e "${target}" ]] || exit 92
cp "${source}" "${target}"
chmod 666 "${target}"
RSYNC
chmod +x "${FIXTURE}/bin/python" "${FIXTURE}/bin/rsync"
export PATH="${FIXTURE}/bin:${PATH}" PYTHON_BIN="${FIXTURE}/bin/python"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/2_dataset_specific_preprocessing/bassez_metadata_step.sh"
ecoda_init_run stage2 bassez-wrapper >/dev/null
export ECODA_RUN_ID ECODA_RUN_ROOT
printf 'fixture source manifest\n' > "${ECODA_RUN_ROOT}/manifests/source.manifest"
printf 'fixture runtime identity\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"
INPUT="${HPC_SCRATCH_DIR}/Bassez/output/$(ecoda_view_output_name Bassez benchmark_analysis)"
mkdir -p "${INPUT%/*}"
"${REAL_PYTHON}" - "${INPUT}" <<'PY'
import sys
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
x = ad.AnnData(X=np.ones((4, 3), dtype=np.float32),
               obs=pd.DataFrame({"Sample":["A","A","B","B"],
                                 "cellSubType":pd.Categorical(["CD4",None,"","Treg"]),
                                 "cellType":pd.Categorical(["T","B","Cancer","T"])},
                                index=["c1","c2","c3","c4"]))
x.layers["counts"] = sp.csr_matrix(np.ones((4,3), dtype=np.int32))
x.uns["upstream"] = "preserve"
x.write_h5ad(sys.argv[1])
PY
ecoda_write_checksum "${INPUT}"
INPUT_MD5="${ECODA_CHECKSUM_MD5}"
REFERENCE="${FIXTURE}/input.tsv"
printf '%s\t%s\tlegacy\t%s.md5\t-\t-\t%s\n' \
  "${INPUT}" "${INPUT_MD5}" "${INPUT}" "$(_stage2_bassez_sha "${INPUT}.md5")" > "${REFERENCE}"

# Only scheduler/ownership effects are stubbed; the scientific metadata worker
# and checksum/manifest validations execute against tiny real H5AD fixtures.
ecoda_require_input_ownership() {
  printf 'ownership:%s:%s\n' "$1" "$2" >> "${EVENTS}"
  [[ "${BLOCK_INPUT:-0}" == 0 ]]
}
ecoda_validate_input_artifact() {
  printf 'producer:%s:%s:%s\n' "$1" "$2" "$3" >> "${EVENTS}"
  [[ "${BLOCK_PRODUCER:-0}" == 0 ]]
}
stage2_bassez_metadata_bind "${REFERENCE}" fixture-release
stage2_bassez_metadata_load
OUTPUT="$(stage2_bassez_metadata_output)"
NAS_OUTPUT="$(stage2_bassez_metadata_nas_output)"
[[ "${OUTPUT}" == "${ECODA_RUN_ROOT}/artifacts/datasets/Bassez/output/${INPUT##*/}" ]]
[[ "${NAS_OUTPUT}" == "${NAS_TARGET_DIR}/releases/fixture-release/${ECODA_RUN_ID}/datasets/Bassez/output/${INPUT##*/}" ]]
ROW="$(bash "${PROJECT_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" --step-config bassez_h5ad_metadata)"
case "${ROW}" in bassez_h5ad_metadata\|*\|python_bassez_metadata\|bassez_h5ad_metadata.py\|"${OUTPUT};${OUTPUT}.repair.json"\|-) ;; *) echo "wrong exact step row: ${ROW}" >&2; exit 1 ;; esac
printf 'bassez_h5ad_metadata\t%s\t%s\t-\t%s\n' "$(stage2_bassez_metadata_script)" \
  "$(stage2_bassez_metadata_outputs)" "$(ecoda_owner_dir stage2 bassez_h5ad_metadata)" > "${ECODA_RUN_ROOT}/manifests/steps.tsv"
_ecoda_expand_output_selection stage2 "${ECODA_RUN_ROOT}/manifests/steps.tsv"
[[ ${#ECODA_OUTPUT_PATHS[@]} -eq 4 ]]
for expected in "${OUTPUT}" "${OUTPUT}.repair.json" "${NAS_OUTPUT}" "${NAS_OUTPUT}.repair.json"; do
  case " ${ECODA_OUTPUT_PATHS[*]} " in *" ${expected} "*) ;; *) exit 1 ;; esac
done
: > "${EVENTS}"
stage2_bassez_metadata_run create > "${FIXTURE}/create.json"
[[ "$(wc -l < "${EVENTS}" | tr -d '[:space:]')" == 2 ]]
[[ -s "${OUTPUT}" && -s "${OUTPUT}.repair.json" ]]
ecoda_write_checksum "${OUTPUT}.repair.json"
stage2_bassez_metadata_validate "${OUTPUT}" > "${FIXTURE}/validate.json"
BEFORE="$(wc -l < "${CALLS}" | tr -d '[:space:]')"
if BLOCK_INPUT=1 stage2_bassez_metadata_run validate > /dev/null 2>&1; then
  echo "active input conflict was ignored" >&2; exit 1
fi
[[ "$(wc -l < "${CALLS}" | tr -d '[:space:]')" == "${BEFORE}" ]]
printf '#!/bin/bash\nexit 17\n' > "${FIXTURE}/bin/failure"
chmod +x "${FIXTURE}/bin/failure"
: > "${EVENTS}"
if PYTHON_BIN="${FIXTURE}/bin/failure" stage2_bassez_metadata_run validate > /dev/null 2>&1; then
  echo "failed metadata worker reported success" >&2; exit 1
fi
[[ "$(wc -l < "${EVENTS}" | tr -d '[:space:]')" == 2 ]]

mkdir -p "${NAS_TARGET_DIR}/unrelated"
printf 'untouched\n' > "${NAS_TARGET_DIR}/unrelated/artifact"
chmod 600 "${NAS_TARGET_DIR}/unrelated/artifact"
INPUT_MODE="$(_ecoda_artifact_mode_bits "${INPUT}")"
OUTPUT_MODE="$(_ecoda_artifact_mode_bits "${OUTPUT}")"
stage2_bassez_metadata_sync > "${FIXTURE}/sync.log"
[[ "$(wc -l < "${RSYNC_CALLS}" | tr -d '[:space:]')" == 2 ]]
[[ "$(_ecoda_artifact_mode_bits "${NAS_OUTPUT}")" == 666 && "$(_ecoda_artifact_mode_bits "${NAS_OUTPUT%/*}")" == 777 ]]
[[ "$(_ecoda_artifact_mode_bits "${INPUT}")" == "${INPUT_MODE}" && "$(_ecoda_artifact_mode_bits "${OUTPUT}")" == "${OUTPUT_MODE}" ]]
[[ "$(_ecoda_artifact_mode_bits "${NAS_TARGET_DIR}/unrelated/artifact")" == 600 ]]
stage2_bassez_metadata_validate "${NAS_OUTPUT}" > "${FIXTURE}/nas-validate.log"
stage2_bassez_metadata_sync > "${FIXTURE}/sync-noop.log"
[[ "$(wc -l < "${RSYNC_CALLS}" | tr -d '[:space:]')" == 2 ]]
[[ "$(ecoda_md5_file "${INPUT}")" == "${INPUT_MD5}" ]]

# Exercise the production terminal NAS guard after successful publication.
# Source only this function: sourcing the full audit CLI would execute it.
eval "$(sed -n '/^_audit_bassez_nas_artifact() {$/,/^}$/p' "${ROOT}/src/utils/bash/ecoda_run_audit.sh")"
declare -F _audit_bassez_nas_artifact >/dev/null
_audit_die() { echo "$*" >&2; return 1; }
for path in "${NAS_OUTPUT}" "${NAS_OUTPUT}.repair.json"; do
  _audit_bassez_nas_artifact "${path}" > "${FIXTURE}/nas-terminal-valid.log"
  cp "${path}" "${FIXTURE}/nas-payload.backup"
  printf 'corrupted after successful sync\n' >> "${path}"
  if _audit_bassez_nas_artifact "${path}" > "${FIXTURE}/nas-terminal-corrupt.log" 2>&1; then
    echo "terminal audit accepted changed NAS payload: ${path}" >&2; exit 1
  fi
  cp "${FIXTURE}/nas-payload.backup" "${path}"
  _audit_bassez_nas_artifact "${path}" > "${FIXTURE}/nas-terminal-restored.log"
done
for target in "${NAS_OUTPUT}" "${NAS_OUTPUT}.md5" "${NAS_OUTPUT%/*}"; do
  if [[ -d "${target}" ]]; then bad=755; good=777; else bad=444; good=666; fi
  [[ "${target}" != "${NAS_OUTPUT}" ]] || bad=644
  chmod "${bad}" "${target}"
  if _audit_bassez_nas_artifact "${NAS_OUTPUT}" > "${FIXTURE}/nas-terminal-mode.log" 2>&1; then
    echo "terminal audit accepted incorrect NAS mode: ${target}" >&2; exit 1
  fi
  chmod "${good}" "${target}"
done
_audit_bassez_nas_artifact "${NAS_OUTPUT}" > "${FIXTURE}/nas-terminal-final.log"
printf 'corrupted\n' >> "${NAS_OUTPUT}.repair.json"
CORRUPTED_MD5="$(ecoda_md5_file "${NAS_OUTPUT}.repair.json")"
if stage2_bassez_metadata_sync > "${FIXTURE}/sync-corrupted.log" 2>&1; then
  echo "invalid existing NAS artifact was overwritten or accepted" >&2; exit 1
fi
[[ "$(ecoda_md5_file "${NAS_OUTPUT}.repair.json")" == "${CORRUPTED_MD5}" ]]
[[ "$(wc -l < "${RSYNC_CALLS}" | tr -d '[:space:]')" == 2 ]]

# Corrupt bound metadata fails closed before the interpreter executes.
BINDING="${ECODA_RUN_ROOT}/manifests/bassez_metadata.json"
cp "${BINDING}" "${FIXTURE}/binding.backup"
printf 'corrupted\n' >> "${BINDING}"
if stage2_bassez_metadata_load; then echo "changed binding was accepted" >&2; exit 1; fi
cp "${FIXTURE}/binding.backup" "${BINDING}"

# Strict producer mode is explicit and invokes the stronger guard twice.
PRODUCER_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/parent-stage3"
KEY="$(printf '%s' "${INPUT}" | shasum -a 256 | awk '{print substr($1,1,32)}')"
RECORD="${PRODUCER_ROOT}/manifests/artifacts/${KEY}.record"
mkdir -p "${RECORD%/*}"
printf 'PATH=%s\nSIZE=%s\nMD5=%s\nRUN_ID=parent-stage3\nPRODUCER=stage3\nSTATE=PUBLISHED\n' \
  "${INPUT}" "$(wc -c < "${INPUT}" | tr -d '[:space:]')" "${INPUT_MD5}" > "${RECORD}"
ecoda_init_run stage2 bassez-strict >/dev/null
export ECODA_RUN_ID ECODA_RUN_ROOT
printf 'fixture source\n' > "${ECODA_RUN_ROOT}/manifests/source.manifest"
printf 'fixture runtime\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"
printf '%s\t%s\tproducer\t%s\tparent-stage3\tstage3\t%s\n' \
  "${INPUT}" "${INPUT_MD5}" "${RECORD}" "$(_stage2_bassez_sha "${RECORD}")" > "${REFERENCE}"
stage2_bassez_metadata_bind "${REFERENCE}" fixture-release
: > "${EVENTS}"
stage2_bassez_metadata_run create > "${FIXTURE}/strict-create.json"
[[ "$(grep -c '^ownership:' "${EVENTS}")" == 2 && "$(grep -c '^producer:' "${EVENTS}")" == 2 ]]
BEFORE="$(wc -l < "${CALLS}" | tr -d '[:space:]')"
if BLOCK_PRODUCER=1 stage2_bassez_metadata_run validate > /dev/null 2>&1; then
  echo "invalid producer ownership was ignored" >&2; exit 1
fi
[[ "$(wc -l < "${CALLS}" | tr -d '[:space:]')" == "${BEFORE}" ]]

for forbidden in '--force' '--ds_name Bassez' '--datasets Bassez,Adams'; do
  if bash "${PROJECT_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --datasets Bassez --steps bassez_h5ad_metadata --bassez-reference "${REFERENCE}" \
    --derivative-release fixture-release ${forbidden} > "${FIXTURE}/forbidden.log" 2>&1; then
    echo "unscoped metadata selection accepted: ${forbidden}" >&2; exit 1
  fi
done
echo "PASS: exact Stage2 metadata row, 4 output owners, input/producer guards, bound references, create-only roots and scoped NAS permissions"
