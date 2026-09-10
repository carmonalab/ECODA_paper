#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-watchdog-ids.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home/scratch/ECODA_paper" \
  "${TMP_DIR}/home/reference_atlases/sketched_200ct" \
  "${TMP_DIR}/nas" "${TMP_DIR}/logs" "${TMP_DIR}/tmp"
cat > "${TMP_DIR}/bin/sacct" <<'STUB'
#!/bin/bash
set -euo pipefail
job=""
while [[ $# -gt 0 ]]; do
  case "$1" in -j) job="$2"; shift 2 ;; *) shift ;; esac
done
case "${job}" in
  3001) printf '3001|COMPLETED|0:0\n3001_1|OUT_OF_MEMORY|0:0\n' ;;
  3002|4002|5002|6002) printf '%s|COMPLETED|0:0\n%s_1|COMPLETED|0:0\n' "${job}" "${job}" ;;
  4001) printf '4001|COMPLETED|0:0\n4001_1|OUT_OF_MEMORY|0:0\n' ;;
  5001) printf '5001|COMPLETED|0:0\n5001_1|OUT_OF_MEMORY|0:0\n' ;;
  6001) printf '6001|COMPLETED|0:0\n6001_1|OUT_OF_MEMORY|0:0\n' ;;
  *) printf '%s|COMPLETED|0:0\n' "${job}" ;;
esac
STUB
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "${SBATCH_ID:?}"
STUB
chmod +x "${TMP_DIR}/bin/sacct" "${TMP_DIR}/bin/sbatch"
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

# Watchdogs now validate the same commit-keyed source/runtime identity that
# production workers receive.  Keep the fixture complete but small by copying
# only the repository files required by the immutable snapshot contract.
SNAPSHOT_COMMIT="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SNAPSHOT_ROOT="${TMP_DIR}/snapshots/${SNAPSHOT_COMMIT}"
SOURCE_ROOT="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY}/source.tar"
mkdir -p "${SOURCE_ROOT}" "${SOURCE_IDENTITY}"
cp -R "${ROOT}/src" "${SOURCE_ROOT}/src"
cp -R "${ROOT}/aux" "${SOURCE_ROOT}/aux"
for source_file in config_helper.R datasets.json pixi.toml pixi.lock; do
  cp "${ROOT}/${source_file}" "${SOURCE_ROOT}/${source_file}"
done
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_ROOT}" .
SOURCE_ARCHIVE_SHA256="$(sha256_file "${SOURCE_ARCHIVE}")"
SOURCE_CONFIG_SHA256="$(sha256_file "${SOURCE_ROOT}/config_helper.R")"
SOURCE_DATASETS_SHA256="$(sha256_file "${SOURCE_ROOT}/datasets.json")"
SOURCE_TOML_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.toml")"
SOURCE_LOCK_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.lock")"
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
printf 'COMPLETE\n' > "${SNAPSHOT_ROOT}/COMPLETE"
chmod -R a-w "${SOURCE_ROOT}" "${SOURCE_IDENTITY}" "${SNAPSHOT_ROOT}/COMPLETE"

HOST_ENV_PREFIX="${TMP_DIR}/host/.pixi/envs/py-cuda13"
mkdir -p "${HOST_ENV_PREFIX}/bin" "${HOST_ENV_PREFIX}/lib"
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

RUNTIME_DIR="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runtime/watchdog-test"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
mkdir -p "${RUNTIME_DIR}"
printf 'watchdog runtime fixture\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA256="$(sha256_file "${RUNTIME_IMAGE}")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_BUILD_GIT_REVISION=watchdog-runtime-build
IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}
IMAGE_PATH=${RUNTIME_IMAGE}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.49.0
APPTAINER_VERSION=1.3.2
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}
EOF
chmod 444 "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}"
chmod 555 "${RUNTIME_DIR}"

export HPC_SCRATCH_DIR="${TMP_DIR}/home/scratch/ECODA_paper"
export NAS_PREFIX="${TMP_DIR}/nas"
export NAS_BASE_DIR="${TMP_DIR}/nas/DataCollections"
export NAS_TARGET_DIR="${TMP_DIR}/nas"
export ECODA_LOGS_DIR="${TMP_DIR}/logs"
export TMPDIR="${TMP_DIR}/tmp"
write_run_identity() {
  local run="$1" stage="$2" run_id="${1##*/}"
  local source_copy="${run}/manifests/source.manifest"
  local runtime_identity="${run}/manifests/runtime.identity"
  mkdir -p "${run}/manifests" "${run}/status" "${run}/logs"
  cp "${SOURCE_MANIFEST}" "${source_copy}"
  chmod 600 "${source_copy}"
  printf 'RUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
    "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA256}" \
    "${ECODA_RUNTIME_MANIFEST_SHA256}" "${ECODA_RUNTIME_IMAGE_SIZE}" \
    "${ECODA_RUNTIME_MANIFEST_SIZE}" "${SOURCE_TOML_SHA256}" \
    "${SOURCE_LOCK_SHA256}" > "${runtime_identity}"
  chmod 600 "${runtime_identity}"
  if [[ "${stage}" == "stage3" ]]; then
    printf 'STAGE=stage3\nRUN_ID=%s\nSTATE=ACTIVE\nSOURCE_MANIFEST=%s\nSOURCE_MANIFEST_COPY=%s\nSOURCE_ROOT=%s\nSOURCE_COMMIT=%s\nSOURCE_ARCHIVE_PATH=%s\nSOURCE_ARCHIVE_SHA256=%s\nSOURCE_CONFIG_HELPER_SHA256=%s\nSOURCE_DATASETS_SHA256=%s\nSOURCE_PIXI_TOML_SHA256=%s\nSOURCE_PIXI_LOCK_SHA256=%s\nSOURCE_AUX_ROOT=%s\nSOURCE_SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4\nRUNTIME_IDENTITY=%s\nRUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
      "${run_id}" "${SOURCE_MANIFEST}" "${source_copy}" "${SOURCE_ROOT}" \
      "${SNAPSHOT_COMMIT}" "${SOURCE_ARCHIVE}" "${SOURCE_ARCHIVE_SHA256}" \
      "${SOURCE_CONFIG_SHA256}" "${SOURCE_DATASETS_SHA256}" "${SOURCE_TOML_SHA256}" \
      "${SOURCE_LOCK_SHA256}" "${SOURCE_ROOT}/aux" "${runtime_identity}" \
      "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA256}" \
      "${ECODA_RUNTIME_MANIFEST_SHA256}" "${ECODA_RUNTIME_IMAGE_SIZE}" \
      "${ECODA_RUNTIME_MANIFEST_SIZE}" "${SOURCE_TOML_SHA256}" \
      "${SOURCE_LOCK_SHA256}" > "${run}/metadata"
  else
    printf 'STAGE=stage4\nRUN_ID=%s\nSTATE=ACTIVE\nSOURCE_MANIFEST=%s\nSOURCE_MANIFEST_FORMAT=1\nSOURCE_ROOT=%s\nSOURCE_COMMIT=%s\nSOURCE_ARCHIVE_PATH=%s\nSOURCE_ARCHIVE_SHA256=%s\nCONFIG_HELPER_SHA256=%s\nDATASETS_SHA256=%s\nPIXI_TOML_SHA256=%s\nPIXI_LOCK_SHA256=%s\nAUX_ROOT=%s\nSCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4\nRUNTIME_IDENTITY=%s\nRUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
      "${run_id}" "${source_copy}" "${SOURCE_ROOT}" "${SNAPSHOT_COMMIT}" \
      "${SOURCE_ARCHIVE}" "${SOURCE_ARCHIVE_SHA256}" "${SOURCE_CONFIG_SHA256}" \
      "${SOURCE_DATASETS_SHA256}" "${SOURCE_TOML_SHA256}" "${SOURCE_LOCK_SHA256}" \
      "${SOURCE_ROOT}/aux" "${runtime_identity}" "${RUNTIME_IMAGE}" \
      "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA256}" \
      "${ECODA_RUNTIME_MANIFEST_SHA256}" "${ECODA_RUNTIME_IMAGE_SIZE}" \
      "${ECODA_RUNTIME_MANIFEST_SIZE}" "${SOURCE_TOML_SHA256}" \
      "${SOURCE_LOCK_SHA256}" > "${run}/metadata"
  fi
  chmod 600 "${run}/metadata"
}
export HOME_REF_DIR="${TMP_DIR}/home/reference_atlases/sketched_200ct"
export ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}"
export ECODA_HOST_PYTHON_SHA256="${HOST_PYTHON_SHA256}"
export ECODA_HOST_RSCRIPT_SHA256="${HOST_RSCRIPT_SHA256}"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_MODE=host
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export ECODA_RUNTIME_IMAGE_SHA256="${RUNTIME_IMAGE_SHA256}"
export ECODA_RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")"
export ECODA_RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
export ECODA_RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
export USER_EMAIL="test@example.invalid"


write_artifact_record() {
  local path="$1" run="$2" producer="$3" digest size key canonical
  canonical="$(realpath "${path}")"
  digest="$(md5_file "${path}")"
  size="$(wc -c < "${path}" | tr -d '[:space:]')"
  key="$(printf '%s' "${canonical}" | sha256sum | cut -d' ' -f1)"
  mkdir -p "${run}/manifests/artifacts"
  printf 'PATH=%s\nSIZE=%s\nMD5=%s\nRUN_ID=%s\nPRODUCER=%s\nSTATE=PUBLISHED\n' \
    "${canonical}" "${size}" "${digest}" "${run##*/}" "${producer}" \
    > "${run}/manifests/artifacts/${key:0:32}.record"
}
write_checksum() {
  local path="$1" digest
  digest="$(md5_file "${path}")"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' \
    "${digest}" "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
}
OUTPUT_NAME="$(jq -r '.Stephenson.views.benchmark_analysis.output_file_name' "${ROOT}/datasets.json")"
SOURCE_H5AD="${TMP_DIR}/home/scratch/ECODA_paper/Stephenson/output/${OUTPUT_NAME}"
mkdir -p "$(dirname "${SOURCE_H5AD}")"
pixi run python -c 'import anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,sys; x=sp.csr_matrix(np.ones((2,2000),dtype="float32")); a=ad.AnnData(X=x,obs=pd.DataFrame({"Sample":["s1","s2"]},index=["c1","c2"]),var=pd.DataFrame({"hvg_rank":np.arange(2000,dtype=float)},index=[f"g{i}" for i in range(2000)])); a.layers["counts"]=x.copy(); a.obsm["X_pca_batch_effect_uncorrected_hvg2000"]=np.ones((2,2),dtype="float32"); a.write_h5ad(sys.argv[1])' "${SOURCE_H5AD}"
write_checksum "${SOURCE_H5AD}"
OUTPUT_BATCH_NAME="$(jq -r '.Stephenson.views.batch_effect_uncorrected.output_file_name' "${ROOT}/datasets.json")"
BATCH_H5AD="${TMP_DIR}/home/scratch/ECODA_paper/Stephenson/output/${OUTPUT_BATCH_NAME}"
cp "${SOURCE_H5AD}" "${BATCH_H5AD}"
write_checksum "${BATCH_H5AD}"
run_stage3() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage3"
  local owner_dir="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage3/Stephenson_batch_effect_uncorrected"
  mkdir -p "${owner_dir}"
  write_run_identity "${run}" stage3
  printf 'Stephenson\tbatch_effect_uncorrected\n' > "${run}/manifests/selection.tsv"
  printf 'Stephenson\tbatch_effect_uncorrected\n' > "${run}/manifests/pending.tsv"
  printf 'Stephenson\t%s\n' "${owner_dir}" > "${run}/manifests/owners.tsv"
  write_checksum "${run}/manifests/selection.tsv"
  write_checksum "${run}/manifests/pending.tsv"
  printf 'RUN_ID=stage3\nSTATE=ACTIVE\nSTAGE=stage3\nKEY=Stephenson/batch_effect_uncorrected\n' > "${owner_dir}/owner"
  write_artifact_record "${BATCH_H5AD}" "${run}" stage3
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" NAS_TARGET_DIR="${NAS_TARGET_DIR}" \
    ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" TMPDIR="${TMPDIR}" \
    ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
    ECODA_SOURCE_ROOT="${SOURCE_ROOT}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
    ECODA_SOURCE_MANIFEST_RUN="${run}/manifests/source.manifest" \
    ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_RUNTIME_MODE=host \
    ECODA_RUNTIME_PROFILE=stage3 ECODA_RUNTIME_IDENTITY="${run}/manifests/runtime.identity" \
    ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
    ECODA_RUNTIME_IMAGE_SHA256="${RUNTIME_IMAGE_SHA256}" \
    ECODA_RUNTIME_MANIFEST_SHA256="${ECODA_RUNTIME_MANIFEST_SHA256}" \
    ECODA_RUNTIME_IMAGE_SIZE="${ECODA_RUNTIME_IMAGE_SIZE}" \
    ECODA_RUNTIME_MANIFEST_SIZE="${ECODA_RUNTIME_MANIFEST_SIZE}" \
    PREPROCESS_RUN_ROOT="${run}" PREPROCESS_PENDING_MANIFEST="${run}/manifests/pending.tsv" \
    SBATCH_ID=3002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${ROOT}/src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh" \
      stage3 "${run}/manifests/selection.tsv" 3001 128G 256G shared-cpu 1000
  [[ "$(grep '^STATE=' "${run}/status/watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/watchdog")" == 2 ]]
  [[ "$(grep -c '^SCHEDULER_ID=3001$' "${run}/status/watchdog")" == 1 ]]
  [[ "$(grep -c '^SCHEDULER_ID=3002$' "${run}/status/watchdog")" == 1 ]]
}
run_stage4_prepare() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage4"
  local union="${run}/datasets/Stephenson/union/union.h5ad" chunk
  local owner_dir="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage4/Stephenson"
  write_run_identity "${run}" stage4
  mkdir -p "$(dirname "${union}")" "${run}/datasets/Stephenson/chunks" "${owner_dir}"
  printf 'Stephenson\tbenchmark_analysis\n' > "${run}/manifests/selection.tsv"
  cp "${run}/manifests/selection.tsv" "${run}/manifests/runnable_selection.tsv"
  write_checksum "${run}/manifests/selection.tsv"
  write_checksum "${run}/manifests/runnable_selection.tsv"
  printf 'RUN_ID=stage4\nSTATE=ACTIVE\nSTAGE=stage4\nKEY=Stephenson\n' > "${owner_dir}/owner"
  pixi run python -c 'import importlib.util,sys; from pathlib import Path; s=importlib.util.spec_from_file_location("p","src/4_cell_type_annotation/1.1_prepare_chunks.py"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); m.build_union([Path(sys.argv[1])],Path(sys.argv[2]),"Sample")' "${SOURCE_H5AD}" "${union}"
  source_md5="$(md5_file "${SOURCE_H5AD}")"
  source_size="$(wc -c < "${SOURCE_H5AD}" | tr -d '[:space:]')"
  printf '[{"md5":"%s","path":"%s","size":%s}]\n' "${source_md5}" "${SOURCE_H5AD}" "${source_size}" > "${run}/datasets/Stephenson/source_artifacts.json"
  chunk="${run}/datasets/Stephenson/chunks/chunk_1.txt"
  printf '%s\ns1\ns2\n' "${union}" > "${chunk}"
  printf 'Stephenson\tbenchmark_analysis\t%s\n' "${run}" > "${run}/manifests/preparation.tsv"
  write_checksum "${run}/manifests/preparation.tsv"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" NAS_TARGET_DIR="${NAS_TARGET_DIR}" \
    ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" TMPDIR="${TMPDIR}" \
    ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
    ECODA_SOURCE_ROOT="${SOURCE_ROOT}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
    ECODA_RUNTIME_MODE=host ECODA_RUNTIME_PROFILE=stage4 \
    ECODA_RUNTIME_IDENTITY="${run}/manifests/runtime.identity" \
    ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
    ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_ROOT}/aux" \
    SCGATE_DB_PATH="${SOURCE_ROOT}/aux/scGateDB.rds" \
    SBATCH_ID=4002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${SOURCE_ROOT}/src/4_cell_type_annotation/1.3_prepare_chunks_watchdog.sh" \
      stage4 "${run}/manifests/preparation.tsv" 4001 32G 64G shared-cpu 1000
  [[ "$(grep '^STATE=' "${run}/status/preparation_watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/preparation_watchdog")" == 2 ]]
}
run_stage4_annotation() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage4"
  local union="${run}/datasets/Stephenson/union/union.h5ad"
  local chunk="${run}/datasets/Stephenson/chunks/chunk_1.txt"
  write_run_identity "${run}" stage4
  mkdir -p "${run}/datasets/Stephenson/annotations" "$(dirname "${chunk}")" "$(dirname "${union}")"
  pixi run python -c 'import importlib.util,sys; from pathlib import Path; s=importlib.util.spec_from_file_location("p","src/4_cell_type_annotation/1.1_prepare_chunks.py"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); m.build_union([Path(sys.argv[1])],Path(sys.argv[2]),"Sample")' "${SOURCE_H5AD}" "${union}"
  write_checksum "${union}"
  printf '%s\ns1\n' "${union}" > "${chunk}"
  printf 'Stephenson\t%s\t%s\n' "${chunk}" "${run}/datasets/Stephenson/annotations" > "${run}/manifests/chunks.tsv"
  write_checksum "${run}/manifests/chunks.tsv"
  pixi run python -c 'import hashlib,pandas as pd,sys; from pathlib import Path; p=Path(sys.argv[1]); d={"Sample":["s1"],"cell_barcode":["c1"],"layer1":["T"],"layer2":["T"],"layer3":["T"],"layer_1":["1"],"layer_2":["2"],"layer_3":["3"],"layer_4":["4"],"layer_5":["5"],"layer_6":["6"],"scATOMIC_pred":["T"],"classification_confidence":[.9],"S.Score":[.1],"G2M.Score":[.2],"Phase":["G1"]}; pd.DataFrame(d).to_feather(p); p.with_name(p.name+".md5").write_text(f"MD5={hashlib.md5(p.read_bytes()).hexdigest()}\nSIZE={p.stat().st_size}\nPATH={p}\n")' "${run}/datasets/Stephenson/annotations/annotations_chunk_1.feather"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" NAS_TARGET_DIR="${NAS_TARGET_DIR}" \
    ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" TMPDIR="${TMPDIR}" \
    ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
    ECODA_SOURCE_ROOT="${SOURCE_ROOT}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
    ECODA_RUNTIME_MODE=host ECODA_RUNTIME_PROFILE=stage4 \
    ECODA_RUNTIME_IDENTITY="${run}/manifests/runtime.identity" \
    ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
    ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_ROOT}/aux" \
    SCGATE_DB_PATH="${SOURCE_ROOT}/aux/scGateDB.rds" \
    SBATCH_ID=5002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${SOURCE_ROOT}/src/4_cell_type_annotation/1.2_annotation_watchdog.sh" \
      stage4 "${run}/manifests/chunks.tsv" 5001 32G 64G shared-cpu 1000
  [[ "$(grep '^STATE=' "${run}/status/annotation_watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/annotation_watchdog")" == 2 ]]
  [[ "$(grep -c '^SCHEDULER_ID=5001$' "${run}/status/annotation_watchdog")" == 1 ]]
  [[ "$(grep -c '^SCHEDULER_ID=5002$' "${run}/status/annotation_watchdog")" == 1 ]]
}
run_stage4_merge() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage4"
  local union="${run}/datasets/Stephenson/union/union.h5ad"
  write_run_identity "${run}" stage4
  mkdir -p "${run}/datasets/Stephenson" "${TMP_DIR}/home/scratch/ECODA_paper/Stephenson/output"
  printf 'Stephenson\tbenchmark_analysis\t%s\n' "${run}" > "${run}/manifests/merge.tsv"
  write_checksum "${run}/manifests/merge.tsv"
  mkdir -p "$(dirname "${union}")"
  pixi run python -c 'import importlib.util,sys; from pathlib import Path; s=importlib.util.spec_from_file_location("p","src/4_cell_type_annotation/1.1_prepare_chunks.py"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); m.build_union([Path(sys.argv[1])],Path(sys.argv[2]),"Sample")' "${SOURCE_H5AD}" "${union}"
  pixi run python -c 'import anndata as ad,sys; p=sys.argv[1]; a=ad.read_h5ad(p); a.obs["layer1"]=["T","B"]; a.obs["layer2"]=["T","B"]; a.obs["layer3"]=["T","B"]; a.obs["layer_1"]=["1","1"]; a.obs["layer_2"]=["2","2"]; a.obs["layer_3"]=["3","3"]; a.obs["layer_4"]=["4","4"]; a.obs["layer_5"]=["5","5"]; a.obs["layer_6"]=["6","6"]; a.obs["scATOMIC_pred"]=["T","B"]; a.obs["classification_confidence"]=[.9,.8]; a.obs["S.Score"]=[.1,.2]; a.obs["G2M.Score"]=[.2,.3]; a.obs["Phase"]=["G1","G2"]; a.write_h5ad(p)' "${SOURCE_H5AD}"
  write_checksum "${SOURCE_H5AD}"
  write_checksum "${union}"
  source_record="${SOURCE_H5AD}|$(md5_file "${SOURCE_H5AD}")|$(wc -c < "${SOURCE_H5AD}" | tr -d '[:space:]')"
  union_md5="$(md5_file "${union}")"
  union_size="$(wc -c < "${union}" | tr -d '[:space:]')"
  printf 'STATE=OK\nDATASET=Stephenson\nVIEWS=benchmark_analysis\nSOURCE_H5ADS=%s\nSOURCE_RECORDS=%s\nUNION_PATH=%s\nUNION_MD5=%s\nUNION_SIZE=%s\n' \
    "${SOURCE_H5AD}" "${source_record}" "${union}" "${union_md5}" "${union_size}" > "${run}/datasets/Stephenson/merge.ok"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" NAS_TARGET_DIR="${NAS_TARGET_DIR}" \
    ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" TMPDIR="${TMPDIR}" \
    ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
    ECODA_SOURCE_ROOT="${SOURCE_ROOT}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
    ECODA_RUNTIME_MODE=host ECODA_RUNTIME_PROFILE=stage4 \
    ECODA_RUNTIME_IDENTITY="${run}/manifests/runtime.identity" \
    ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
    ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_ROOT}/aux" \
    SCGATE_DB_PATH="${SOURCE_ROOT}/aux/scGateDB.rds" \
    SBATCH_ID=6002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${SOURCE_ROOT}/src/4_cell_type_annotation/3.3_merge_watchdog.sh" \
      stage4 "${run}/manifests/merge.tsv" 6001 32G 64G shared-cpu 1000 >/dev/null
  [[ "$(grep '^STATE=' "${run}/status/merge_watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/merge_watchdog")" == 2 ]]
}
run_stage3
run_stage4_prepare
run_stage4_annotation
run_stage4_merge
echo "pipeline watchdog scheduler IDs: OK"
