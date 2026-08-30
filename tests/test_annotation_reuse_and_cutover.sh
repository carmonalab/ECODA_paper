#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-annotation-reuse.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/manifests" \
  "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/status" \
  "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/chunks" \
  "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse/datasets/Adams/annotations" \
  "${TMP_DIR}/home/scratch/ECODA_paper/Adams/output"
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
RUN_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/reuse"
printf 'STAGE=stage4\nRUN_ID=reuse\nSTATE=ACTIVE\n' > "${RUN_ROOT}/metadata"
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
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" ANNOTATION_SUBMITTER_TEST=1 USER_EMAIL=test@example.invalid \
  bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
    --skip-prepare --reuse-run reuse --views benchmark_analysis >/dev/null
runs=("${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs"/*)
[[ ${#runs[@]} -eq 1 && "${runs[0]}" == "${RUN_ROOT}" ]]
[[ "$(grep -c '1.2_prepare_chunks_worker.sh' "${CAPTURE}" || true)" == 0 ]]
[[ "$(grep -c '2.1_run_worker.sh' "${CAPTURE}")" == 1 ]]
[[ "$(grep -c '3.2_merge_worker.sh' "${CAPTURE}")" == 1 ]]
set +e
NO_REUSE="$(HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" bash "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" --skip-prepare --datasets Adams 2>&1)"
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
