#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-annotation-reuse.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

SNAPSHOT_ID="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SNAPSHOT_ROOT="${TMP_DIR}/snapshots/${SNAPSHOT_ID}"
SOURCE_TREE="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY_DIR="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY_DIR}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY_DIR}/source.tar"
HOST_ENV_PREFIX="${TMP_DIR}/host/.pixi/envs/py-cuda13"
RUNTIME_DIR="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runtime/reuse-test"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
HOME_REF_DIR="${TMP_DIR}/home/reference_atlases/sketched_200ct"
NAS_ROOT="${TMP_DIR}/nas"
LOGS_DIR="${TMP_DIR}/logs"
RUN_TMP_DIR="${TMP_DIR}/tmp"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/manifests" \
  "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/status" \
  "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/chunks" \
  "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/annotations" \
  "${TMP_DIR}/home/scratch/ECODA_paper/Adams/output" \
  "${SOURCE_TREE}" "${SOURCE_IDENTITY_DIR}" "${HOST_ENV_PREFIX}/bin" \
  "${RUNTIME_DIR}" "${HOME_REF_DIR}" "${NAS_ROOT}" "${LOGS_DIR}" "${RUN_TMP_DIR}"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '88000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

cat > "${HOST_ENV_PREFIX}/bin/python" <<STUB
#!/bin/bash
set -euo pipefail
exec "${ROOT}/.pixi/envs/py-cuda13/bin/python" "\$@"
STUB
chmod +x "${HOST_ENV_PREFIX}/bin/python"
printf '#!/bin/bash\nexit 0\n' > "${HOST_ENV_PREFIX}/bin/Rscript"
chmod +x "${HOST_ENV_PREFIX}/bin/Rscript"

OUTPUT_NAME="$(jq -r '.Adams.views.benchmark_analysis.output_file_name' "${ROOT}/datasets.json")"
SOURCE_H5AD="${TMP_DIR}/home/scratch/ECODA_paper/Adams/output/${OUTPUT_NAME}"
UNION_H5AD="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/union/union.h5ad"
mkdir -p "$(dirname "${UNION_H5AD}")"
pixi run python -c 'import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, sys; a=ad.AnnData(X=sp.csr_matrix(np.ones((2,2),dtype="float32")),obs=pd.DataFrame({"Sample":["s1","s2"]},index=["c1","c2"]),var=pd.DataFrame(index=["g1","g2"])); a.layers["counts"]=a.X.copy(); a.write_h5ad(sys.argv[1])' "${SOURCE_H5AD}"
pixi run python -c 'import importlib.util,sys; s=importlib.util.spec_from_file_location("p","src/4_cell_type_annotation/1.1_prepare_chunks.py"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); m.build_union([__import__("pathlib").Path(sys.argv[1])],__import__("pathlib").Path(sys.argv[2]),"Sample")' "${SOURCE_H5AD}" "${UNION_H5AD}"
CHUNK="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/chunks/chunk_1.txt"
printf '%s\ns1\ns2\n' "${UNION_H5AD}" > "${CHUNK}"
pixi run python -c 'import hashlib,pandas as pd,sys; p=sys.argv[1]; d={"Sample":["s1","s2"],"cell_barcode":["c1","c2"],"layer1":["T","B"],"layer2":["T","B"],"layer3":["T","B"],"layer_1":["1","1"],"layer_2":["2","2"],"layer_3":["3","3"],"layer_4":["4","4"],"layer_5":["5","5"],"layer_6":["6","6"],"scATOMIC_pred":["T","B"],"classification_confidence":[.9,.8],"S.Score":[.1,.2],"G2M.Score":[.2,.3],"Phase":["G1","G2"]}; df=pd.DataFrame(d); df.to_feather(p); pth=__import__("pathlib").Path(p); pth.with_name(pth.name+".md5").write_text(f"MD5={hashlib.md5(pth.read_bytes()).hexdigest()}\nSIZE={pth.stat().st_size}\nPATH={pth}\n")' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/annotations/annotations_chunk_1.feather"
pixi run python src/utils/py/annotation_contract.py --path "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/annotations/annotations_chunk_1.feather" --require-sidecar
# Stage 4 reuse must be bound to a complete immutable source snapshot rather
# than the mutable checkout that created the synthetic artifacts.
cp -R "${ROOT}/src" "${SOURCE_TREE}/src"
cp -R "${ROOT}/aux" "${SOURCE_TREE}/aux"
cp "${ROOT}/datasets.json" "${SOURCE_TREE}/datasets.json"
cp "${ROOT}/config_helper.R" "${SOURCE_TREE}/config_helper.R"
cp "${ROOT}/pixi.toml" "${SOURCE_TREE}/pixi.toml"
cp "${ROOT}/pixi.lock" "${SOURCE_TREE}/pixi.lock"
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_TREE}" .
SOURCE_ARCHIVE_SHA="$(sha256_of "${SOURCE_ARCHIVE}")"
SOURCE_CONFIG_SHA="$(sha256_of "${SOURCE_TREE}/config_helper.R")"
SOURCE_DATASETS_SHA="$(sha256_of "${SOURCE_TREE}/datasets.json")"
SOURCE_TOML_SHA="$(sha256_of "${SOURCE_TREE}/pixi.toml")"
SOURCE_LOCK_SHA="$(sha256_of "${SOURCE_TREE}/pixi.lock")"
cat > "${SOURCE_MANIFEST}" <<EOF
FORMAT=1
SOURCE_ROOT=${SOURCE_TREE}
SOURCE_COMMIT=${SNAPSHOT_ID}
SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}
SOURCE_ARCHIVE_SHA256=${SOURCE_ARCHIVE_SHA}
CONFIG_HELPER_SHA256=${SOURCE_CONFIG_SHA}
DATASETS_SHA256=${SOURCE_DATASETS_SHA}
PIXI_TOML_SHA256=${SOURCE_TOML_SHA}
PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA}
AUX_ROOT=${SOURCE_TREE}/aux
SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4
EOF
printf 'COMPLETE\n' > "${SNAPSHOT_ROOT}/COMPLETE"

printf 'synthetic format-2 runtime image\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA="$(sha256_of "${RUNTIME_IMAGE}")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_BUILD_GIT_REVISION=synthetic-runtime-build
IMAGE_SHA256=${RUNTIME_IMAGE_SHA}
IMAGE_PATH=${RUNTIME_IMAGE}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.49.0
APPTAINER_VERSION=1.3.2
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA}
EOF
RUNTIME_MANIFEST_SHA="$(sha256_of "${RUNTIME_MANIFEST}")"
RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
chmod -R a-w "${SNAPSHOT_ROOT}" "${RUNTIME_DIR}"

RUN_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse"
cp "${SOURCE_MANIFEST}" "${RUN_ROOT}/manifests/source.manifest"
cat > "${RUN_ROOT}/manifests/runtime.identity" <<EOF
RUNTIME_IMAGE=${RUNTIME_IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
RUNTIME_IMAGE_SHA256=${RUNTIME_IMAGE_SHA}
RUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA}
RUNTIME_IMAGE_SIZE=${RUNTIME_IMAGE_SIZE}
RUNTIME_MANIFEST_SIZE=${RUNTIME_MANIFEST_SIZE}
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA}
EOF
chmod a-w "${RUN_ROOT}/manifests/source.manifest" "${RUN_ROOT}/manifests/runtime.identity"
cat > "${RUN_ROOT}/metadata" <<EOF
STAGE=stage4
RUN_ID=reuse
STATE=ACTIVE
SOURCE_MANIFEST=${RUN_ROOT}/manifests/source.manifest
SOURCE_MANIFEST_FORMAT=1
SOURCE_ROOT=${SOURCE_TREE}
SOURCE_COMMIT=${SNAPSHOT_ID}
SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}
SOURCE_ARCHIVE_SHA256=${SOURCE_ARCHIVE_SHA}
CONFIG_HELPER_SHA256=${SOURCE_CONFIG_SHA}
DATASETS_SHA256=${SOURCE_DATASETS_SHA}
PIXI_TOML_SHA256=${SOURCE_TOML_SHA}
PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA}
AUX_ROOT=${SOURCE_TREE}/aux
SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4
RUNTIME_IMAGE=${RUNTIME_IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
RUNTIME_IMAGE_SHA256=${RUNTIME_IMAGE_SHA}
RUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA}
RUNTIME_IMAGE_SIZE=${RUNTIME_IMAGE_SIZE}
RUNTIME_MANIFEST_SIZE=${RUNTIME_MANIFEST_SIZE}
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA}
EOF
printf 'Adams\tbenchmark_analysis\n' > "${RUN_ROOT}/manifests/selection.tsv"
selection="${RUN_ROOT}/manifests/selection.tsv"
digest="$(md5sum "${selection}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${selection}" | tr -d '[:space:]')" "${selection}" > "${selection}.md5"
SOURCE_MD5="$(md5sum "${SOURCE_H5AD}" | cut -d' ' -f1)"
SOURCE_SIZE="$(wc -c < "${SOURCE_H5AD}" | tr -d '[:space:]')"
printf '[{"md5":"%s","path":"%s","size":%s}]\n' \
  "${SOURCE_MD5}" "${SOURCE_H5AD}" "${SOURCE_SIZE}" \
  > "${RUN_ROOT}/datasets/Adams/source_artifacts.json"
printf 'Adams\tbenchmark_analysis\t%s\n' "${RUN_ROOT}" > "${RUN_ROOT}/manifests/preparation.tsv"
printf 'Adams\t%s\t%s\n' "${CHUNK}" "${RUN_ROOT}/datasets/Adams/annotations" > "${RUN_ROOT}/manifests/chunks.tsv"
printf 'Adams\tbenchmark_analysis\t%s\n' "${RUN_ROOT}" > "${RUN_ROOT}/manifests/merge.tsv"
OWNER_DIR="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage4/Adams"
mkdir -p "${OWNER_DIR}"
printf 'RUN_ID=reuse\nSTATE=OK\nSTAGE=stage4\nKEY=Adams\n' > "${OWNER_DIR}/owner"
printf 'Adams\t%s\n' "${OWNER_DIR}" > "${RUN_ROOT}/manifests/owners.tsv"
HOME="${TMP_DIR}/home" \
PATH="${TMP_DIR}/bin:${PATH}" \
HPC_SCRATCH_DIR="${TMP_DIR}/home/scratch/ECODA_paper" \
ECODA_LOGS_DIR="${LOGS_DIR}" \
TMPDIR="${RUN_TMP_DIR}" \
HOME_REF_DIR="${HOME_REF_DIR}" \
ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
NAS_TARGET_DIR="${NAS_ROOT}" \
ANNOTATION_SUBMITTER_TEST=1 USER_EMAIL=test@example.invalid \
  bash "${SOURCE_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
    --skip-prepare --reuse-run reuse --views benchmark_analysis >/dev/null
runs=("${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs"/*)
[[ ${#runs[@]} -eq 1 && "${runs[0]}" == "${RUN_ROOT}" ]]
[[ "$(grep -c '1.2_prepare_chunks_worker.sh' "${CAPTURE}" || true)" == 0 ]]
[[ "$(grep -c '2.1_run_worker.sh' "${CAPTURE}")" == 1 ]]
[[ "$(grep -c '3.2_merge_worker.sh' "${CAPTURE}")" == 1 ]]
set +e
NO_REUSE="$(
  HOME="${TMP_DIR}/home" \
  PATH="${TMP_DIR}/bin:${PATH}" \
  HPC_SCRATCH_DIR="${TMP_DIR}/home/scratch/ECODA_paper" \
  NAS_TARGET_DIR="${NAS_ROOT}" \
  bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
    --skip-prepare --datasets Adams 2>&1
)"
RC=$?
set -e
for legacy in \
  "${ROOT}/src/4_cell_type_annotation/1_prepare_chunks.sh" \
  "${ROOT}/src/4_cell_type_annotation/2_submit_hpc_array.sh" \
  "${ROOT}/src/4_cell_type_annotation/3_submit_merge.sh"; do
  set +e
  legacy_output="$(bash "${legacy}" 2>&1)"
  legacy_rc=$?
  set -e
  [[ ${legacy_rc} -eq 64 ]]
  case "${legacy_output}" in *"legacy Stage 4"* ) ;; *) exit 1 ;; esac
done
[[ ${RC} -ne 0 ]]
case "${NO_REUSE}" in *"--skip-prepare requires --reuse-run"*) ;; *) exit 1 ;; esac
echo "annotation reuse and legacy cutover: OK"
