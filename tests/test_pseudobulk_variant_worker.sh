#!/bin/bash
# Execute the real shell worker against a temporary runtime/Slurm boundary.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-pb-worker.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'rm -rf "${TMP}"' EXIT
TREE="${TMP}/source"
RUN="${TMP}/scratch/_ecoda_runs/pb-worker-test"
WORKER_REL="src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
mkdir -p "${TREE}/src/5_run_benchmark_methods/run_r_sample_embedding_methods" \
  "${TREE}/src/utils/bash" "${RUN}/manifests"
cp "${ROOT}/${WORKER_REL}" "${TREE}/${WORKER_REL}"
cat > "${TREE}/src/slurm_config.sh" <<'EOF'
export PROJECT_ROOT="${TEST_PROJECT_ROOT}"
export DATASETS_JSON_FILE="${TEST_PROJECT_ROOT}/datasets.json"
export HPC_SCRATCH_DIR="${TEST_SCRATCH_ROOT}"
export NAS_TARGET_DIR="${TEST_SCRATCH_ROOT}/nas"
export LOGS_DIR="${TEST_SCRATCH_ROOT}/logs"
export PIXI_RSCRIPT="${TEST_RSCRIPT}"
export PYTHON_BIN="${TEST_PYTHON}"
EOF
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'EOF'
ecoda_runtime_reexec_worker() { return 0; }
ecoda_runtime_validate_bound_run() { return 0; }
EOF
cat > "${TREE}/src/utils/bash/ecoda_run_common.sh" <<'EOF'
ecoda_validate_run_id() { [[ "$1" == pb-worker-test ]]; }
ecoda_validate_run_owned_path() { case "$1" in "$2"/*) return 0;; *) return 1;; esac; }
EOF
: > "${TREE}/src/utils/bash/ecoda_stage5_policy.sh"
cat > "${TREE}/src/utils/bash/worker_retry.sh" <<'EOF'
worker_clear_retry_count() { return 0; }
worker_requeue_if_transient() { return 1; }
EOF
cat > "${TMP}/Rscript" <<'EOF'
#!/bin/bash
printf '%s\n' "$@" > "${TEST_CAPTURE}"
EOF
chmod +x "${TMP}/Rscript"
printf 'immutable source fixture\n' > "${TMP}/source.manifest"
cp "${TMP}/source.manifest" "${RUN}/manifests/source.manifest"
: > "${RUN}/manifests/runtime.identity"
run_worker() {
  env -i PATH="${PATH}" TEST_PROJECT_ROOT="${TREE}" \
    TEST_SCRATCH_ROOT="${TMP}/scratch" TEST_RSCRIPT="${TMP}/Rscript" \
    TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" \
    TEST_CAPTURE="${TMP}/args" ECODA_SOURCE_SNAPSHOT_REQUIRED=1 \
    ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TMP}/source.manifest" \
    ECODA_RUNTIME_IMAGE="${TMP}/runtime.sif" ECODA_RUNTIME_MANIFEST="${TMP}/runtime.manifest" \
    ECODA_RUNTIME_MODE=host ECODA_RUNTIME_IN_CONTAINER=1 \
    ECODA_RUN_ID=pb-worker-test ECODA_RUN_ROOT="${RUN}" \
    METHOD="$1" ANALYSIS_PASS="${2:-}" \
    ECODA_ARTIFACT_BINDING="${TEST_ARTIFACT_BINDING:-}" \
    ANALYSIS_ROOT="${TEST_ANALYSIS_ROOT:-}" ANALYSIS_NAS_ROOT="${TEST_ANALYSIS_NAS_ROOT:-}" \
    ANALYSIS_MANIFEST="${RUN}/manifests/worker.tsv" \
    EXECUTION_LOG_DIR="${RUN}/logs" SLURM_ARRAY_TASK_ID=1 SLURM_ARRAY_JOB_ID=101 \
    bash "${TREE}/${WORKER_REL}"
}
for method in prepare_pseudobulk pseudobulk; do
  printf 'Toy\tbenchmark_analysis\t%s\thvg2000_bl\n' "${method}" > "${RUN}/manifests/worker.tsv"
  run_worker "${method}"
  grep -qx -- --pseudobulk_variants "${TMP}/args"
  grep -qx hvg2000_bl "${TMP}/args"
  ! grep -qx -- --combo "${TMP}/args"
  ! grep -qx -- --force "${TMP}/args"
  if [[ "${method}" == prepare_pseudobulk ]]; then
    grep -q '/1.1.1_prepare_pseudobulk.R$' "${TMP}/args"
  fi
  printf 'Toy\tbenchmark_analysis\t%s\n' "${method}" > "${RUN}/manifests/worker.tsv"
  run_worker "${method}"
  ! grep -qx -- --pseudobulk_variants "${TMP}/args"
done
for token in invalid hvg2000_bl,hvg500; do
  printf 'Toy\tbenchmark_analysis\tpseudobulk\t%s\n' "${token}" > "${RUN}/manifests/worker.tsv"
  if run_worker pseudobulk >/dev/null 2>&1; then
    echo "ERROR: invalid shell shard token was accepted: ${token}" >&2; exit 1
  fi
done
printf 'Toy\tbenchmark_analysis\tcomposition\thvg2000_bl\n' > "${RUN}/manifests/worker.tsv"
if run_worker composition >/dev/null 2>&1; then
  echo "ERROR: wrong-method selector was accepted" >&2; exit 1
fi
printf 'Toy\tbatch_effect_uncorrected\tpseudobulk\thvg2000_bl\n' > "${RUN}/manifests/worker.tsv"
if run_worker pseudobulk uncorrected >/dev/null 2>&1; then
  echo "ERROR: batch selector was accepted" >&2; exit 1
fi
# The real worker resolves an explicit canonical input reference while writing
# into its exact bound derivative root; arbitrary roots fail before R starts.
mkdir -p "${TREE}/src/utils/py" "${TMP}/scratch/nas" "${TMP}/canonical"
cp "${ROOT}/src/utils/py/derivative_artifacts.py" "${TREE}/src/utils/py/derivative_artifacts.py"
cp "${ROOT}/src/utils/bash/ecoda_derivative_artifacts.sh" "${TREE}/src/utils/bash/ecoda_derivative_artifacts.sh"
printf '{"Toy":{"columns":{"label":"biology"},"views":{"benchmark_analysis":{"output_file_name":"toy.h5ad"}}}}\n' > "${TREE}/datasets.json"
printf 'runtime fixture\n' > "${RUN}/manifests/runtime.identity"
printf 'Toy\tbenchmark_analysis\n' > "${RUN}/manifests/selection.tsv"
printf 'Toy\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n' > "${RUN}/manifests/worker.tsv"
"${ROOT}/.pixi/envs/default/bin/python" - "${TMP}" <<'PY'
import hashlib, pathlib, sys
root = pathlib.Path(sys.argv[1]); source = root / 'canonical/toy.h5ad'
source.write_bytes(b'canonical unchanged fixture')
md5 = hashlib.md5(source.read_bytes()).hexdigest(); size = source.stat().st_size
proof = pathlib.Path(str(source) + '.md5')
proof.write_text(f'MD5={md5}\nSIZE={size}\nPATH={source}\n')
(root/'inputs.tsv').write_text(f'Toy\tbenchmark_analysis\t{source}\t{md5}\t{size}\t{proof}\n')
PY
TEST_ARTIFACT_BINDING="$("${ROOT}/.pixi/envs/default/bin/python" \
  "${TREE}/src/utils/py/derivative_artifacts.py" bind --stage stage5 \
  --run-root "${RUN}" --release-id test-release --nas-root "${TMP}/scratch/nas" \
  --config "${TREE}/datasets.json" --selection "${RUN}/manifests/selection.tsv" \
  --input-artifacts "${TMP}/inputs.tsv")"
TEST_ANALYSIS_ROOT="${RUN}/artifacts/benchmark"
TEST_ANALYSIS_NAS_ROOT="${TMP}/scratch/nas/releases/test-release/pb-worker-test/benchmark"
run_worker pseudobulk
grep -qx "${TMP}/canonical" "${TMP}/args"
grep -qx "${RUN}/artifacts/benchmark/results" "${TMP}/args"
TEST_ANALYSIS_ROOT="${TMP}/scratch/benchmark"
if run_worker pseudobulk >/dev/null 2>&1; then
  echo 'ERROR: derivative worker accepted canonical output root' >&2; exit 1
fi
printf 'PASS: R worker passes exact PB shard selection, binds derivative roots, and rejects contradictory rows\n'
