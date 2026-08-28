#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-sync.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${TMP_DIR}/scratch/benchmark/embeddings" "${TMP_DIR}/nas/project"
cat > "${TMP_DIR}/bin/rsync" <<'STUB'
#!/bin/bash
set -euo pipefail
files_from=""
for arg in "$@"; do
  case "${arg}" in
    --files-from=*) files_from="${arg#--files-from=}" ;;
  esac
done
n="$#"
eval "src=\${$((n - 1))}"
eval "dst=\${$n}"
mkdir -p "${dst}"
if [[ -n "${files_from}" ]]; then
  while IFS= read -r rel || [[ -n "${rel}" ]]; do
    [[ -n "${rel}" ]] || continue
    mkdir -p "${dst}/$(dirname "${rel}")"
    cp "${src}/${rel}" "${dst}/${rel}"
  done < "${files_from}"
else
  cp -R "${src}". "${dst}"
fi
STUB
chmod +x "${TMP_DIR}/bin/rsync"
export PATH="${TMP_DIR}/bin:${PATH}"
export HPC_SCRATCH_DIR="${TMP_DIR}/scratch"
export LOGS_DIR="${TMP_DIR}/logs"
export NAS_TARGET_DIR="${TMP_DIR}/nas/project"
export ANALYSIS_ROOT="${TMP_DIR}/scratch/benchmark"
export ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/benchmark"
export USER_EMAIL="test@example.invalid"
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
export HPC_SCRATCH_DIR="${TMP_DIR}/scratch" NAS_TARGET_DIR="${TMP_DIR}/nas/project" \
  ANALYSIS_ROOT="${TMP_DIR}/scratch/benchmark" ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/benchmark"
mkdir -p "${NAS_TARGET_DIR}" "${ANALYSIS_NAS_ROOT}"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/5_run_benchmark_methods/benchmark_submit_common.sh"
ECODA_RUN_ID="sync_run"
ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
export ECODA_RUN_ID ECODA_RUN_ROOT
printf 'Adams\tbenchmark_analysis\tbenchmark_analysis\n' > "${ECODA_RUN_ROOT}/manifests/selection.tsv"
export ECODA_SELECTION_MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv" ECODA_EXACT_SELECTION=0
export ANALYSIS_ROOT ANALYSIS_NAS_ROOT EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs" \
  ANALYSIS_LOG_PREFIX="execution_times_"
DATASET_NAMES=(Adams)
LABELS=(mrvi)
export DATASET_NAMES LABELS
for n in 1000 2000 3000; do
  path="${ANALYSIS_ROOT}/embeddings/Adams_hvg${n}_mrvi_dists.feather"
  pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"s1":[1.0,0.0],"s2":[0.0,1.0]},index=["s1","s2"]).to_feather(sys.argv[1])' "${path}"
  ecoda_write_checksum "${path}"
done
mkdir -p "${ECODA_RUN_ROOT}/logs"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[1.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
ecoda_write_checksum "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
printf 'unrelated remote\n' > "${ANALYSIS_NAS_ROOT}/unrelated.txt"
analysis_merge_sync_cleanup "${LABELS[@]}"
[[ -s "${ANALYSIS_NAS_ROOT}/checksums.md5" ]]
[[ -s "${ANALYSIS_NAS_ROOT}/embeddings/execution_times.feather" ]]
[[ -s "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather" ]]
[[ -s "${ANALYSIS_NAS_ROOT}/unrelated.txt" ]]
[[ ! -e "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather" ]]
[[ ! -d "${ANALYSIS_ROOT}/.sync.lock" ]]
[[ "$(ecoda_owner_state "$(ecoda_owner_dir stage5 "sync/${ANALYSIS_ROOT}")")" == "OK" ]]
[[ "$(grep -c '^checksums.md5$' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv")" == 1 ]]
[[ "$(grep -c 'unrelated' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv" || true)" == 0 ]]
[[ "$(grep -c '^embeddings/Adams_hvg1000_mrvi_dists.feather$' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv")" == 1 ]]
(
  cd "${ANALYSIS_NAS_ROOT}"
  md5sum unrelated.txt >> checksums.md5
)
grep -q 'unrelated.txt' "${ANALYSIS_NAS_ROOT}/checksums.md5"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[2.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
ecoda_write_checksum "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
analysis_merge_sync_cleanup "${LABELS[@]}"
[[ -s "${ANALYSIS_NAS_ROOT}/unrelated.txt" ]]
[[ ! -e "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather" ]]
# A selected checksum failure must stop before merge cleanup and preserve run logs.
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[3.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
ecoda_write_checksum "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
rm -f "${ANALYSIS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather.md5"
set +e
analysis_merge_sync_cleanup "${LABELS[@]}"
rc=$?
set -e
[[ ${rc} -ne 0 ]]
[[ -s "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather" ]]
echo "benchmark sync owner: OK"
