#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-benchmark-deps.XXXXXX")"
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
printf '78000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
printf '_debug\tbenchmark_analysis\tselected\n_debug\tbatch_effect_uncorrected\tselected\nAdams\tbenchmark_analysis\tselected\nBassez\tbenchmark_analysis\tselected\n' > "${TMP_DIR}/selection.tsv"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid BENCHMARK_MATRIX_TEST=1 \
  bash "${ROOT}/src/5_run_benchmark_methods/1_submit_hpc_array.sh" \
  --selection-file "${TMP_DIR}/selection.tsv" --methods mofa,gloscope,mrvi >/dev/null
CALLS="${CAPTURE}"
[[ "$(wc -l < "${CALLS}" | tr -d '[:space:]')" == 19 ]]
RUN_ROOTS=("${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/"*)
[[ ${#RUN_ROOTS[@]} -eq 1 ]]
RUN_ROOT="${RUN_ROOTS[0]}"
MANIFEST_NAMES=(
  matrix_benchmark_analysis_prepare_pseudobulk.tsv
  matrix_batch_effect_uncorrected_prepare_pseudobulk.tsv
  matrix_benchmark_analysis_mofa.tsv
  matrix_batch_effect_uncorrected_mofa.tsv
  matrix_benchmark_analysis_gloscope_cpu.tsv
  matrix_batch_effect_uncorrected_gloscope.tsv
  matrix_benchmark_analysis_mrvi_default_gpu.tsv
  matrix_benchmark_analysis_mrvi_cpu.tsv
  matrix_batch_effect_uncorrected_mrvi.tsv
)
MANIFEST_COUNTS=(2 1 2 1 15 1 3 6 1)
MANIFEST_COLUMNS=(3 3 3 3 4 3 4 4 3)
for manifest_idx in "${!MANIFEST_NAMES[@]}"; do
  manifest="${RUN_ROOT}/manifests/${MANIFEST_NAMES[${manifest_idx}]}"
  [[ -s "${manifest}" ]]
  [[ "$(wc -l < "${manifest}" | tr -d '[:space:]')" == "${MANIFEST_COUNTS[${manifest_idx}]}" ]]
  awk -F '\t' -v expected="${MANIFEST_COLUMNS[${manifest_idx}]}" \
    'NF != expected {exit 1}' "${manifest}"
done
MOFA_BENCH="$(sed -n '5p' "${CALLS}")"
MOFA_BATCH="$(sed -n '7p' "${CALLS}")"
case "${MOFA_BENCH}" in *"--dependency=afterok:780002"*) ;; *) echo "benchmark MOFA dependency missing" >&2; exit 1 ;; esac
case "${MOFA_BATCH}" in *"--dependency=afterok:780004"*) ;; *) echo "batch-view MOFA dependency missing" >&2; exit 1 ;; esac
for line in 9 11 13 15 17; do
  if sed -n "${line}p" "${CALLS}" | grep -q -- '--dependency=afterok:'; then
    echo "independent Stage 5 method carried pseudobulk dependency" >&2
    exit 1
  fi
done
if sed -n '1,18p' "${CALLS}" | grep -q -- '--wait'; then
  echo "an independent array waited before aggregate submission" >&2
  exit 1
fi
case "$(sed -n '19p' "${CALLS}")" in *"matrix_gate.sh"*) ;; *) echo "aggregate gate was not last" >&2; exit 1 ;; esac
echo "benchmark dependency edges: OK"
