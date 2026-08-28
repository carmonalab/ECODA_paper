#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-matrix-watchdog.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home/scratch/ECODA_paper/run/manifests" "${TMP_DIR}/home/scratch/ECODA_paper/run/status"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sacct" <<'STUB'
#!/bin/bash
case "$*" in
  *"-j 1001"*) printf '1001|COMPLETED|0:0\n1001_1|COMPLETED|0:0\n1001_2|OUT_OF_MEMORY|0:0\n' ;;
  *"-j 1002"*) printf '1002|COMPLETED|0:0\n1002_1|COMPLETED|0:0\n' ;;
  *) printf 'COMPLETED\n' ;;
esac
STUB
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
printf '%s\n' "$*" >> "${CAPTURE}"
printf '1002\n'
STUB
chmod +x "${TMP_DIR}/bin/sacct" "${TMP_DIR}/bin/sbatch"
RUN_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/run"
MANIFEST="${RUN_ROOT}/manifests/mrvi.tsv"
printf 'Adams\tbatch_effect_uncorrected\tmrvi\nBassez\tbatch_effect_uncorrected\tmrvi\n' > "${MANIFEST}"
WATCHDOG_OUTPUT="$(HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" ANALYSIS_PASS=uncorrected MATRIX_WATCHDOG_MAX_POLLS=1 \
  bash "${ROOT}/src/5_run_benchmark_methods/matrix_watchdog.sh" "${RUN_ROOT}" mrvi "${MANIFEST}" 1001 128G 256G shared-gpu 4 "${ROOT}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh" --gpus=1)"
STATUS="${RUN_ROOT}/status/watchdogs/mrvi.status"
[[ "$(grep '^STATE=' "${STATUS}")" == "STATE=OK" ]]
RETRY="${RUN_ROOT}/manifests/mrvi.retry_1.tsv"
[[ "$(cat "${RETRY}")" == $'Bassez\tbatch_effect_uncorrected\tmrvi' ]]
[[ "$(grep -c '^SCHEDULER_ID=' "${STATUS}")" == 2 ]]
[[ "$(grep -c '^SCHEDULER_ID=1001$' "${STATUS}")" == 1 ]]
[[ "$(grep -c '^SCHEDULER_ID=1002$' "${STATUS}")" == 1 ]]
case "${WATCHDOG_OUTPUT}" in *"BATCH_EFFECT_RETRY_ARRAY_JOB_ID=1002"*) ;; *) echo "batch retry marker missing" >&2; exit 1 ;; esac
if grep -q 'BENCHMARK_MANIFEST=' "${CAPTURE}"; then
  echo "batch retry exported BENCHMARK_MANIFEST" >&2
  exit 1
fi
SCHEDULER_MANIFEST="${RUN_ROOT}/manifests/scheduler_ids.tsv"
printf 'ARRAY\t1001\nWATCHDOG\t1003\n' > "${SCHEDULER_MANIFEST}"
SLURM_JOB_ID=1004 HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  bash "${ROOT}/src/5_run_benchmark_methods/matrix_gate.sh" "${RUN_ROOT}" mrvi "${SCHEDULER_MANIFEST}"
[[ "$(grep '^STATE=' "${RUN_ROOT}/status/aggregate")" == "STATE=OK" ]]
[[ "$(grep -c '^SCHEDULER_ID=' "${RUN_ROOT}/status/aggregate")" == 4 ]]
[[ "$(grep -c '^SCHEDULER_ID=1002$' "${RUN_ROOT}/status/aggregate")" == 1 ]]
echo "benchmark matrix watchdog: OK"
