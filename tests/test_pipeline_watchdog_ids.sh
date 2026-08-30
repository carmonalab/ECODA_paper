#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-watchdog-ids.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home/scratch/ECODA_paper"
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
OUTPUT_NAME="$(jq -r '.Stephenson.views.benchmark_analysis.output_file_name' "${ROOT}/datasets.json")"
SOURCE_H5AD="${TMP_DIR}/home/scratch/ECODA_paper/Stephenson/output/${OUTPUT_NAME}"
mkdir -p "$(dirname "${SOURCE_H5AD}")"
pixi run python -c 'import anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,sys; x=sp.csr_matrix(np.ones((2,2000),dtype="float32")); a=ad.AnnData(X=x,obs=pd.DataFrame({"Sample":["s1","s2"]},index=["c1","c2"]),var=pd.DataFrame({"hvg_rank":np.arange(2000,dtype=float)},index=[f"g{i}" for i in range(2000)])); a.layers["counts"]=x.copy(); a.obsm["X_pca_batch_effect_uncorrected_hvg2000"]=np.ones((2,2),dtype="float32"); a.write_h5ad(sys.argv[1])' "${SOURCE_H5AD}"
write_checksum() {
  local path="$1" digest
  digest="$(md5sum "${path}" | cut -d' ' -f1)"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
}
write_checksum "${SOURCE_H5AD}"
OUTPUT_BATCH_NAME="$(jq -r '.Stephenson.views.batch_effect_uncorrected.output_file_name' "${ROOT}/datasets.json")"
BATCH_H5AD="${TMP_DIR}/home/scratch/ECODA_paper/Stephenson/output/${OUTPUT_BATCH_NAME}"
cp "${SOURCE_H5AD}" "${BATCH_H5AD}"
write_checksum "${BATCH_H5AD}"
run_stage3() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage3"
  mkdir -p "${run}/manifests" "${run}/status" "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage3/Stephenson_batch_effect_uncorrected"
  printf 'Stephenson\tbatch_effect_uncorrected\n' > "${run}/manifests/selection.tsv"
  printf 'Stephenson\tbatch_effect_uncorrected\n' > "${run}/manifests/pending.tsv"
  printf 'Stephenson\t%s\n' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage3/Stephenson_batch_effect_uncorrected" > "${run}/manifests/owners.tsv"
  printf 'RUN_ID=stage3\nSTATE=ACTIVE\nSTAGE=stage3\nKEY=Stephenson/batch_effect_uncorrected\n' > "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage3/Stephenson_batch_effect_uncorrected/owner"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" SBATCH_ID=3002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${ROOT}/src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh" stage3 "${run}/manifests/selection.tsv" 3001 128G 256G shared-cpu 1000
  [[ "$(grep '^STATE=' "${run}/status/watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/watchdog")" == 2 ]]
  [[ "$(grep -c '^SCHEDULER_ID=3001$' "${run}/status/watchdog")" == 1 ]]
  [[ "$(grep -c '^SCHEDULER_ID=3002$' "${run}/status/watchdog")" == 1 ]]
}
run_stage4_prepare() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage4prep"
  local union="${run}/datasets/Stephenson/union/union.h5ad" chunk
  mkdir -p "${run}/manifests" "${run}/status" "$(dirname "${union}")" "${run}/datasets/Stephenson/chunks"
  printf 'STAGE=stage4\nRUN_ID=stage4prep\nSTATE=ACTIVE\n' > "${run}/metadata"
  owner_dir="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage4/Stephenson"
  mkdir -p "${owner_dir}"
  printf 'RUN_ID=stage4prep\nSTATE=ACTIVE\nSTAGE=stage4\nKEY=Stephenson\n' > "${owner_dir}/owner"
  pixi run python -c 'import importlib.util,sys; from pathlib import Path; s=importlib.util.spec_from_file_location("p","src/4_cell_type_annotation/1.1_prepare_chunks.py"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); m.build_union([Path(sys.argv[1])],Path(sys.argv[2]),"Sample")' "${SOURCE_H5AD}" "${union}"
  source_md5="$(md5sum "${SOURCE_H5AD}" | cut -d' ' -f1)"
  source_size="$(wc -c < "${SOURCE_H5AD}" | tr -d '[:space:]')"
  printf '[{"md5":"%s","path":"%s","size":%s}]\n' "${source_md5}" "${SOURCE_H5AD}" "${source_size}" > "${run}/datasets/Stephenson/source_artifacts.json"
  chunk="${run}/datasets/Stephenson/chunks/chunk_1.txt"
  printf '%s\ns1\ns2\n' "${union}" > "${chunk}"
  printf 'Stephenson\tbenchmark_analysis\t%s\n' "${run}" > "${run}/manifests/preparation.tsv"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" SBATCH_ID=4002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${ROOT}/src/4_cell_type_annotation/1.3_prepare_chunks_watchdog.sh" stage4prep "${run}/manifests/preparation.tsv" 4001 32G 64G shared-cpu 1000
  [[ "$(grep '^STATE=' "${run}/status/preparation_watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/preparation_watchdog")" == 2 ]]
}
run_stage4_annotation() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage4annot"
  local union="${run}/datasets/Stephenson/union/union.h5ad"
  local chunk="${run}/datasets/Stephenson/chunks/chunk_1.txt"
  mkdir -p "${run}/manifests" "${run}/status" "${run}/datasets/Stephenson/annotations" "$(dirname "${chunk}")" "$(dirname "${union}")"
  pixi run python -c 'import importlib.util,sys; from pathlib import Path; s=importlib.util.spec_from_file_location("p","src/4_cell_type_annotation/1.1_prepare_chunks.py"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); m.build_union([Path(sys.argv[1])],Path(sys.argv[2]),"Sample")' "${SOURCE_H5AD}" "${union}"
  write_checksum "${union}"
  printf '%s\ns1\n' "${union}" > "${chunk}"
  printf 'Stephenson\t%s\t%s\n' "${chunk}" "${run}/datasets/Stephenson/annotations" > "${run}/manifests/chunks.tsv"
  pixi run python -c 'import hashlib,pandas as pd,sys; from pathlib import Path; p=Path(sys.argv[1]); d={"Sample":["s1"],"cell_barcode":["c1"],"layer1":["T"],"layer2":["T"],"layer3":["T"],"layer_1":["1"],"layer_2":["2"],"layer_3":["3"],"layer_4":["4"],"layer_5":["5"],"layer_6":["6"],"scATOMIC_pred":["T"],"classification_confidence":[.9],"S.Score":[.1],"G2M.Score":[.2],"Phase":["G1"]}; pd.DataFrame(d).to_feather(p); p.with_name(p.name+".md5").write_text(f"MD5={hashlib.md5(p.read_bytes()).hexdigest()}\nSIZE={p.stat().st_size}\nPATH={p}\n")' "${run}/datasets/Stephenson/annotations/annotations_chunk_1.feather"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" SBATCH_ID=5002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${ROOT}/src/4_cell_type_annotation/1.2_annotation_watchdog.sh" stage4annot "${run}/manifests/chunks.tsv" 5001 32G 64G shared-cpu 1000
  [[ "$(grep '^STATE=' "${run}/status/annotation_watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/annotation_watchdog")" == 2 ]]
  [[ "$(grep -c '^SCHEDULER_ID=5001$' "${run}/status/annotation_watchdog")" == 1 ]]
  [[ "$(grep -c '^SCHEDULER_ID=5002$' "${run}/status/annotation_watchdog")" == 1 ]]
}
run_stage4_merge() {
  local run="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/stage4merge"
  mkdir -p "${run}/manifests" "${run}/status" "${TMP_DIR}/home/scratch/ECODA_paper/Stephenson/output"
  printf 'Stephenson\tbenchmark_analysis\t%s\n' "${run}" > "${run}/manifests/merge.tsv"
  union="${run}/datasets/Stephenson/union/union.h5ad"
  mkdir -p "$(dirname "${union}")"
  pixi run python -c 'import importlib.util,sys; from pathlib import Path; s=importlib.util.spec_from_file_location("p","src/4_cell_type_annotation/1.1_prepare_chunks.py"); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); m.build_union([Path(sys.argv[1])],Path(sys.argv[2]),"Sample")' "${SOURCE_H5AD}" "${union}"
  mkdir -p "${run}/datasets/Stephenson"
  pixi run python -c 'import anndata as ad,sys; p=sys.argv[1]; a=ad.read_h5ad(p); a.obs["layer1"]=["T","B"]; a.obs["layer2"]=["T","B"]; a.obs["layer3"]=["T","B"]; a.obs["layer_1"]=["1","1"]; a.obs["layer_2"]=["2","2"]; a.obs["layer_3"]=["3","3"]; a.obs["layer_4"]=["4","4"]; a.obs["layer_5"]=["5","5"]; a.obs["layer_6"]=["6","6"]; a.obs["scATOMIC_pred"]=["T","B"]; a.obs["classification_confidence"]=[.9,.8]; a.obs["S.Score"]=[.1,.2]; a.obs["G2M.Score"]=[.2,.3]; a.obs["Phase"]=["G1","G2"]; a.write_h5ad(p)' "${SOURCE_H5AD}"
  write_checksum "${SOURCE_H5AD}"
  write_checksum "${union}"
  source_record="${SOURCE_H5AD}|$(md5sum "${SOURCE_H5AD}" | cut -d' ' -f1)|$(wc -c < "${SOURCE_H5AD}" | tr -d '[:space:]')"
  union_md5="$(md5sum "${union}" | cut -d' ' -f1)"
  union_size="$(wc -c < "${union}" | tr -d '[:space:]')"
  printf 'STATE=OK\nDATASET=Stephenson\nVIEWS=benchmark_analysis\nSOURCE_H5ADS=%s\nSOURCE_RECORDS=%s\nUNION_PATH=%s\nUNION_MD5=%s\nUNION_SIZE=%s\n' \
    "${SOURCE_H5AD}" "${source_record}" "${union}" "${union_md5}" "${union_size}" > "${run}/datasets/Stephenson/merge.ok"
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" SBATCH_ID=6002 ECODA_ACCOUNTING_EMPTY_GRACE=2 \
    bash "${ROOT}/src/4_cell_type_annotation/3.3_merge_watchdog.sh" stage4merge "${run}/manifests/merge.tsv" 6001 32G 64G shared-cpu 1000 >/dev/null
  [[ "$(grep '^STATE=' "${run}/status/merge_watchdog")" == "STATE=OK" ]]
  [[ "$(grep -c '^SCHEDULER_ID=' "${run}/status/merge_watchdog")" == 2 ]]
}
run_stage3
run_stage4_prepare
run_stage4_annotation
run_stage4_merge
echo "pipeline watchdog scheduler IDs: OK"
