#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-selection.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home"
CAPTURE="${TMP_DIR}/calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '79000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
printf 'Adams\tbenchmark_analysis\tmrvi\nBassez\tbenchmark_analysis\tgloscope\nKfoury\ttrans\ttrans\nKim\tzeroimp\tzeroimp\n' > "${TMP_DIR}/selection.tsv"
OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" BENCHMARK_MATRIX_TEST=1 USER_EMAIL=test@example.invalid \
  bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" --selection-file "${TMP_DIR}/selection.tsv"
)"
RUN_ID="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^BENCHMARK_RUN_ID=//p')"
[[ -n "${RUN_ID}" ]]
RUN_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}"
SELECTION="${RUN_ROOT}/manifests/selection.tsv"
[[ "$(sed -n '3p' "${SELECTION}")" == $'Kfoury\tbenchmark_analysis\ttrans' ]]
[[ "$(sed -n '4p' "${SELECTION}")" == $'Kim\tbenchmark_analysis\tzeroimp' ]]
for row in \
  "benchmark_analysis mrvi Adams" \
  "benchmark_analysis gloscope Bassez" \
  "benchmark_analysis trans Kfoury" \
  "benchmark_analysis zeroimp Kim"; do
  view="${row%% *}"
  rest="${row#* }"
  label="${rest%% *}"
  ds="${rest#* }"
  manifest="${RUN_ROOT}/manifests/matrix_${view}_${label}.tsv"
  [[ "$(wc -l < "${manifest}" | tr -d '[:space:]')" == 1 ]]
  [[ "$(cat "${manifest}")" == $"${ds}\t${view}\t${label}" ]]
done
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == 9 ]]
[[ "$(grep -c -- '--array=1-1' "${CAPTURE}")" == 8 ]]
[[ "$(grep -c 'matrix_gate.sh' "${CAPTURE}")" == 1 ]]
echo "exact benchmark selection and aliases: OK"
