#!/bin/bash
# Focused contract for the compute-node H5AD preflight worker.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-h5ad-preflight.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${TMP_DIR}/nas" "${TMP_DIR}/logs"

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}
md5_file() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | cut -d' ' -f1
  else
    md5 -q "$1"
  fi
}
write_sidecar() {
  local path="$1" digest
  digest="$(md5_file "${path}")"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' \
    "${digest}" "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
}

# Use a complete temporary snapshot rather than the mutable checkout.  The
# worker itself is executed from this tree and all identity paths are retained.
SNAPSHOT_COMMIT="0000000000000000000000000000000000000002"
SNAPSHOT_ROOT="${TMP_DIR}/snapshots/${SNAPSHOT_COMMIT}"
SOURCE_ROOT="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY}/source.tar"
mkdir -p "${SOURCE_ROOT}/src/utils/bash" "${SOURCE_ROOT}/src/utils/py" \
  "${SOURCE_ROOT}/aux" "${SOURCE_IDENTITY}"
cp "${ROOT}/src/slurm_config.sh" "${SOURCE_ROOT}/src/slurm_config.sh"
for source_file in ecoda_run_common.sh ecoda_runtime.sh h5ad_preflight_worker.sh; do
  cp "${ROOT}/src/utils/bash/${source_file}" \
    "${SOURCE_ROOT}/src/utils/bash/${source_file}"
done
cp "${ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
  "${SOURCE_ROOT}/src/utils/py/benchmark_h5ad_contract.py"
for source_file in config_helper.R datasets.json pixi.toml pixi.lock; do
  cp "${ROOT}/${source_file}" "${SOURCE_ROOT}/${source_file}"
done
for source_file in scGateDB.rds genes.blocklist.rds \
  EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
  cp "${ROOT}/aux/${source_file}" "${SOURCE_ROOT}/aux/${source_file}"
done
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_ROOT}" .
SOURCE_ARCHIVE_SHA256="$(sha256_file "${SOURCE_ARCHIVE}")"
SOURCE_CONFIG_SHA256="$(sha256_file "${SOURCE_ROOT}/config_helper.R")"
SOURCE_DATASETS_SHA256="$(sha256_file "${SOURCE_ROOT}/datasets.json")"
SOURCE_TOML_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.toml")"
SOURCE_LOCK_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.lock")"
chmod -R a-w "${SOURCE_ROOT}"
printf '%s\n' \
  'FORMAT=1' \
  "SOURCE_ROOT=${SOURCE_ROOT}" \
  "SOURCE_COMMIT=${SNAPSHOT_COMMIT}" \
  "SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}" \
  "SOURCE_ARCHIVE_SHA256=${SOURCE_ARCHIVE_SHA256}" \
  "CONFIG_HELPER_SHA256=${SOURCE_CONFIG_SHA256}" \
  "DATASETS_SHA256=${SOURCE_DATASETS_SHA256}" \
  "PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}" \
  "PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}" \
  "AUX_ROOT=${SOURCE_ROOT}/aux" \
  'SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4' > "${SOURCE_MANIFEST}"
touch "${SNAPSHOT_ROOT}/COMPLETE"
chmod a-w "${SOURCE_MANIFEST}" "${SOURCE_ARCHIVE}" "${SNAPSHOT_ROOT}/COMPLETE"

HOST_ENV_PREFIX="${TMP_DIR}/host/.pixi/envs/py-cuda13"
mkdir -p "${HOST_ENV_PREFIX}/bin" "${HOST_ENV_PREFIX}/lib" \
  "${TMP_DIR}/home/scratch/ECODA_paper"
cat > "${HOST_ENV_PREFIX}/bin/python" <<'STUB'
#!/bin/bash
exit 0
STUB
cat > "${HOST_ENV_PREFIX}/bin/Rscript" <<'STUB'
#!/bin/bash
exit 0
STUB
chmod +x "${HOST_ENV_PREFIX}/bin/python" "${HOST_ENV_PREFIX}/bin/Rscript"
HOST_PYTHON_SHA256="$(sha256_file "${HOST_ENV_PREFIX}/bin/python")"
HOST_RSCRIPT_SHA256="$(sha256_file "${HOST_ENV_PREFIX}/bin/Rscript")"
export ECODA_HOST_PYTHON_SHA256="${HOST_PYTHON_SHA256}"
export ECODA_HOST_RSCRIPT_SHA256="${HOST_RSCRIPT_SHA256}"

APPTAINER_STUB="${TMP_DIR}/bin/apptainer"
cat > "${APPTAINER_STUB}" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ "${1:-}" == "inspect" ]] && exit 0
exit 1
STUB
chmod +x "${APPTAINER_STUB}"

RUNTIME_ID="h5ad-test-runtime"
RUNTIME_DIR="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runtime/${RUNTIME_ID}"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
mkdir -p "${RUNTIME_DIR}"
printf 'format2-runtime-fixture\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA256="$(sha256_file "${RUNTIME_IMAGE}")"
printf '%s\n' \
  'FORMAT=2' \
  "IMAGE_PATH=${RUNTIME_IMAGE}" \
  "IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}" \
  'RUNTIME_ENV=py-cuda13' \
  'RUNTIME_LAYOUT=relocated' \
  'CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13' \
  'BASE_IMAGE=rockylinux:9' \
  'PIXITAINER_VERSION=0.8.3' \
  'PIXI_VERSION=0.49.0' \
  'APPTAINER_VERSION=1.3.2' \
  'IMAGE_BUILD_GIT_REVISION=immutable-runtime-build' \
  "IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}" \
  "IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}" > "${RUNTIME_MANIFEST}"
chmod 444 "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}"
chmod 555 "${RUNTIME_DIR}"

export PATH="${TMP_DIR}/bin:${PATH}"
export HOME="${TMP_DIR}/home"
export HPC_SCRATCH_DIR="${TMP_DIR}/home/scratch/ECODA_paper"
export NAS_TARGET_DIR="${TMP_DIR}/nas"
export ECODA_LOGS_DIR="${TMP_DIR}/logs"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}"
export ECODA_RUNTIME_MODE=host
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_PROFILE=stage3
export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export APPTAINER_BIN="${APPTAINER_STUB}"
export USER_EMAIL="test@example.invalid"

source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
export PROJECT_ROOT="${ROOT}" DATASETS_JSON_FILE="${ROOT}/datasets.json"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/utils/bash/ecoda_runtime.sh"

RUN_ID="run"
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
mkdir -p "${RUN_ROOT}/manifests" "${RUN_ROOT}/status" "${RUN_ROOT}/logs"
printf 'STAGE=stage3\nRUN_ID=%s\nSTATE=ACTIVE\n' "${RUN_ID}" > "${RUN_ROOT}/metadata"
cp "${SOURCE_MANIFEST}" "${RUN_ROOT}/manifests/source.manifest"
RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")"
printf 'RUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA256}" \
  "${RUNTIME_MANIFEST_SHA256}" "$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')" \
  "$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')" \
  "${SOURCE_TOML_SHA256}" "${SOURCE_LOCK_SHA256}" > "${RUN_ROOT}/manifests/runtime.identity"
chmod 600 "${RUN_ROOT}/manifests/runtime.identity"
export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT="${RUN_ROOT}"

SOURCE="${TMP_DIR}/invalid.h5ad"
printf 'not an h5ad\n' > "${SOURCE}"
write_sidecar "${SOURCE}"
MANIFEST="${RUN_ROOT}/manifests/preflight.tsv"
STATUS_DIR="${RUN_ROOT}/status/preflight"
mkdir -p "${STATUS_DIR}"
printf 'Adams\tbenchmark_analysis\t%s\n' "${SOURCE}" > "${MANIFEST}"

run_worker() {
  local mode="$1"
  ECODA_SOURCE_ROOT="${SOURCE_ROOT}" \
  ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
  ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
  ECODA_RUNTIME_MODE=host ECODA_RUNTIME_IN_CONTAINER=0 \
  ECODA_RUNTIME_PROFILE=stage3 ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" \
  ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" \
  ECODA_SCRATCH_ROOT="${HPC_SCRATCH_DIR}" ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" \
  ECODA_AUX_ROOT="${SOURCE_ROOT}/aux" ECODA_RUN_ROOT="${RUN_ROOT}" \
  ECODA_RUN_ID="${RUN_ID}" H5AD_PREFLIGHT_MANIFEST="${MANIFEST}" \
  H5AD_PREFLIGHT_STATUS_DIR="${STATUS_DIR}" H5AD_PREFLIGHT_RUN_ROOT="${RUN_ROOT}" \
  H5AD_PREFLIGHT_MODE="${mode}" \
  H5AD_PREFLIGHT_PYTHON_BIN="${ROOT}/.pixi/envs/default/bin/python" \
  H5AD_PREFLIGHT_TASK_ID=1 SLURM_SUBMIT_DIR="${ROOT}" \
  bash "${SOURCE_ROOT}/src/utils/bash/h5ad_preflight_worker.sh"
}

run_worker classify
STATUS_FILE="${STATUS_DIR}/Adams__benchmark_analysis.status"
[[ "$(sed -n 's/^STATE=//p' "${STATUS_FILE}")" == "REBUILD" ]]
[[ "$(sed -n 's/^RUN_ID=//p' "${STATUS_FILE}")" == "run" ]]

if run_worker require; then
  echo "require-mode accepted malformed H5AD" >&2
  exit 1
fi
[[ "$(sed -n 's/^STATE=//p' "${STATUS_FILE}")" == "FAIL" ]]
[[ "$(sed -n 's/^RUN_ID=//p' "${STATUS_FILE}")" == "run" ]]

# A valid H5AD is strictly checked and published before its run-owned record is
# reused.  The fixture retains the existing benchmark schema/semantic check.
VALID_SOURCE="${TMP_DIR}/valid.h5ad"
pixi run python -c 'import anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,sys; n=3000; x=sp.csr_matrix(np.ones((2,n),dtype="float32")); a=ad.AnnData(X=x,obs=pd.DataFrame({"Sample":["s1","s2"]},index=["c1","c2"]),var=pd.DataFrame({"hvg_rank":np.arange(1,n+1,dtype=float)},index=[f"g{i}" for i in range(n)])); a.layers["counts"]=x.copy(); [a.obsm.__setitem__(k,np.ones((2,2),dtype="float32")) for k in ["X_pca_benchmark_analysis_hvg1000","X_pca_benchmark_analysis_hvg2000","X_pca_benchmark_analysis_hvg3000","X_pca_harmony_benchmark_analysis_hvg2000"]]; a.write_h5ad(sys.argv[1])' "${VALID_SOURCE}"
write_sidecar "${VALID_SOURCE}"
MANIFEST="${RUN_ROOT}/manifests/preflight-valid.tsv"
STATUS_DIR="${RUN_ROOT}/status/preflight-valid"
mkdir -p "${STATUS_DIR}"
printf 'Adams\tbenchmark_analysis\t%s\n' "${VALID_SOURCE}" > "${MANIFEST}"
run_worker require
STATUS_FILE="${STATUS_DIR}/Adams__benchmark_analysis.status"
[[ "$(sed -n 's/^STATE=//p' "${STATUS_FILE}")" == "OK" ]]
RECORD="$(ecoda_artifact_record_path "${VALID_SOURCE}" "${RUN_ID}")"
[[ -s "${RECORD}" ]]
[[ "$(sed -n 's/^PATH=//p' "${RECORD}")" == "${VALID_SOURCE}" ]]
[[ "$(sed -n 's/^PRODUCER=//p' "${RECORD}")" == stage3_preflight ]]
[[ "$(sed -n 's/^STATE=//p' "${RECORD}")" == PUBLISHED ]]

# A same-size mutation is caught by the strict preflight, not hidden by the
# previously published record.
cp "${VALID_SOURCE}" "${VALID_SOURCE}.saved"
chmod u+w "${VALID_SOURCE}"
printf 'X' | dd of="${VALID_SOURCE}" bs=1 count=1 conv=notrunc >/dev/null 2>&1
if run_worker require; then
  echo "strict H5AD preflight accepted a same-size mutation" >&2
  exit 1
fi
[[ "$(sed -n 's/^STATE=//p' "${STATUS_FILE}")" == "FAIL" ]]
cp "${VALID_SOURCE}.saved" "${VALID_SOURCE}"
chmod a-w "${VALID_SOURCE}"

# Record-only reuse verifies sidecar/size/path metadata and does not invoke a
# second full md5 hash after strict publication.
COUNT_BIN="${TMP_DIR}/count-bin"
MD5_COUNT="${TMP_DIR}/md5.calls"
mkdir -p "${COUNT_BIN}"
cat > "${COUNT_BIN}/md5sum" <<'STUB'
#!/bin/bash
printf '%s\n' "$*" >> "${MD5_COUNT:?}"
exit 1
STUB
chmod +x "${COUNT_BIN}/md5sum"
ORIGINAL_PATH="${PATH}"
export PATH="${COUNT_BIN}:${PATH}"
ecoda_validate_artifact_record "${VALID_SOURCE}" stage3_preflight "${RUN_ID}" >/dev/null
export PATH="${ORIGINAL_PATH}"
[[ ! -s "${MD5_COUNT}" ]]

CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
printf '812345\n'
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
export PATH="${TMP_DIR}/bin:${PATH}"
source "${ROOT}/src/utils/bash/h5ad_preflight_submit.sh"
ecoda_wait_h5ad_preflight_status_files "${MANIFEST}" "${STATUS_DIR}"
runtime_export="$(ecoda_runtime_export_csv stage3 0)"
IFS=, read -r -a runtime_fields <<< "${runtime_export}"
[[ "${#runtime_fields[@]}" -eq 20 ]] || {
  echo "FORMAT=2 runtime export field count mismatch" >&2
  exit 1
}
expected_runtime_fields=(
  "ECODA_RUNTIME_MODE=host"
  "ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}"
  "ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}"
  "ECODA_RUNTIME_PROFILE=stage3"
  "ECODA_APPTAINER_NV=0"
  "ECODA_SOURCE_ROOT=${SOURCE_ROOT}"
  "ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}"
  "ECODA_SOURCE_SNAPSHOT_REQUIRED=1"
  "ECODA_RUN_ID=${RUN_ID}"
)
for runtime_index in "${!expected_runtime_fields[@]}"; do
  [[ "${runtime_fields[runtime_index]}" == "${expected_runtime_fields[runtime_index]}" ]] || {
    echo "FORMAT=2 runtime export field ${runtime_index} mismatch" >&2
    exit 1
  }
done
expected_runtime_keys=(
  ECODA_RUNTIME_MODE
  ECODA_RUNTIME_IMAGE
  ECODA_RUNTIME_MANIFEST
  ECODA_RUNTIME_PROFILE
  ECODA_APPTAINER_NV
  ECODA_SOURCE_ROOT
  ECODA_SOURCE_MANIFEST
  ECODA_SOURCE_SNAPSHOT_REQUIRED
  ECODA_RUN_ID
  ECODA_RUNTIME_IMAGE_SHA256
  ECODA_RUNTIME_MANIFEST_SHA256
  ECODA_RUNTIME_IMAGE_SIZE
  ECODA_RUNTIME_MANIFEST_SIZE
  ECODA_IMAGE_PIXI_TOML_SHA256
  ECODA_IMAGE_PIXI_LOCK_SHA256
  ECODA_RUNTIME_IMAGE_READONLY
  ECODA_RUNTIME_MANIFEST_READONLY
  ECODA_RUNTIME_PARENT_READONLY
  ECODA_HOST_PYTHON_SHA256
  ECODA_HOST_RSCRIPT_SHA256
)
for runtime_index in "${!expected_runtime_keys[@]}"; do
  [[ "${runtime_fields[runtime_index]%%=*}" == "${expected_runtime_keys[runtime_index]}" ]] || {
    echo "FORMAT=2 runtime export field name ${runtime_index} mismatch" >&2
    exit 1
  }
done
[[ "${runtime_fields[18]}" == "ECODA_HOST_PYTHON_SHA256=${HOST_PYTHON_SHA256}" ]] || {
  echo "FORMAT=2 host Python digest mismatch" >&2
  exit 1
}
[[ "${runtime_fields[19]}" == "ECODA_HOST_RSCRIPT_SHA256=${HOST_RSCRIPT_SHA256}" ]] || {
  echo "FORMAT=2 host Rscript digest mismatch" >&2
  exit 1
}


if ECODA_RUN_ID="${RUN_ID},comma" \
  ecoda_runtime_export_csv stage3 0 >/dev/null 2>&1; then
  echo "FORMAT=2 runtime export accepted a comma in ECODA_RUN_ID" >&2
  exit 1
fi
if ECODA_RUN_ID="${RUN_ID}"$'\nnewline' \
  ecoda_runtime_export_csv stage3 0 >/dev/null 2>&1; then
  echo "FORMAT=2 runtime export accepted a newline in ECODA_RUN_ID" >&2
  exit 1
fi

preflight_id="$(ecoda_submit_h5ad_preflight \
  "${MANIFEST}" "${STATUS_DIR}" "${RUN_ROOT}" require shared-cpu 1G 1 \
  "${RUN_ROOT}/logs" test "${SOURCE_ROOT}/src/utils/bash/h5ad_preflight_worker.sh" \
  "${runtime_export}")"
[[ "${preflight_id}" == "812345" ]]
CALLS="$(cat "${CAPTURE}")"
case "${CALLS}" in *"--wait"*"--array=1-1%1"*) ;; *) echo "preflight scheduler array contract missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"H5AD_PREFLIGHT_RUN_ID=run"*) ;; *) echo "preflight run ID was not exported" >&2; exit 1 ;; esac
case "${CALLS}" in *"H5AD_PREFLIGHT_MANIFEST=${MANIFEST}"*) ;; *) echo "preflight manifest path was not exported" >&2; exit 1 ;; esac
case "${CALLS}" in *"H5AD_PREFLIGHT_STATUS_DIR=${STATUS_DIR}"*) ;; *) echo "preflight status path was not exported" >&2; exit 1 ;; esac
case "${CALLS}" in *"H5AD_PREFLIGHT_RUN_ROOT=${RUN_ROOT}"*) ;; *) echo "preflight run root was not exported" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_ROOT=${SOURCE_ROOT}"*"ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}"*) ;; *) echo "preflight snapshot source identity missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}"*"ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}"*) ;; *) echo "preflight versioned runtime identity missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"${SOURCE_ROOT}/src/utils/bash/h5ad_preflight_worker.sh"*) ;; *) echo "preflight escaped immutable worker script" >&2; exit 1 ;; esac

# Source-script containment is checked before sbatch, and a legacy run without
# either run-bound manifest cannot reach the scheduler boundary.
OUTSIDE_SCRIPT="${TMP_DIR}/outside-worker.sh"
printf '#!/bin/bash\n' > "${OUTSIDE_SCRIPT}"
chmod +x "${OUTSIDE_SCRIPT}"
: > "${CAPTURE}"
if ecoda_submit_h5ad_preflight \
  "${MANIFEST}" "${STATUS_DIR}" "${RUN_ROOT}" require shared-cpu 1G 1 \
  "${RUN_ROOT}/logs" test "${OUTSIDE_SCRIPT}" "${runtime_export}" >/dev/null 2>&1; then
  echo "preflight accepted a worker outside the source snapshot" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
for missing_manifest in source.manifest runtime.identity; do
  mv "${RUN_ROOT}/manifests/${missing_manifest}" \
    "${RUN_ROOT}/manifests/${missing_manifest}.saved"
  : > "${CAPTURE}"
  set +e
  ecoda_submit_h5ad_preflight \
    "${MANIFEST}" "${STATUS_DIR}" "${RUN_ROOT}" require shared-cpu 1G 1 \
    "${RUN_ROOT}/logs" test "${SOURCE_ROOT}/src/utils/bash/h5ad_preflight_worker.sh" \
    "${runtime_export}" >/dev/null 2>&1
  missing_rc=$?
  set -e
  [[ ${missing_rc} -ne 0 ]]
  [[ ! -s "${CAPTURE}" ]]
  mv "${RUN_ROOT}/manifests/${missing_manifest}.saved" \
    "${RUN_ROOT}/manifests/${missing_manifest}"
done

# An sbatch failure still returns its numeric scheduler ID and nonzero status.
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
printf '812346\n'
exit 1
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
set +e
failed_id="$(ecoda_submit_h5ad_preflight \
  "${MANIFEST}" "${STATUS_DIR}" "${RUN_ROOT}" require shared-cpu 1G 1 \
  "${RUN_ROOT}/logs" test "${SOURCE_ROOT}/src/utils/bash/h5ad_preflight_worker.sh" \
  "${runtime_export}")"
failed_rc=$?
set -e
[[ "${failed_id}" == "812346" && ${failed_rc} -ne 0 ]]

echo "h5ad preflight worker: OK"
