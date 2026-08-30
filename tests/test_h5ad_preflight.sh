#!/bin/bash
# Focused contract for the compute-node H5AD preflight worker.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-h5ad-preflight.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT

HOME="${TMP_DIR}/home"
export HOME
RUN_ROOT="${HOME}/scratch/ECODA_paper/_ecoda_runs/run"
MANIFEST="${RUN_ROOT}/manifests/preflight.tsv"
STATUS_DIR="${RUN_ROOT}/status/preflight"
SOURCE="${TMP_DIR}/invalid.h5ad"
mkdir -p "$(dirname "${MANIFEST}")"
printf 'not an h5ad\n' > "${SOURCE}"
digest="$(md5sum "${SOURCE}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${SOURCE}" | tr -d '[:space:]')" "${SOURCE}" > "${SOURCE}.md5"
printf 'Adams\tbenchmark_analysis\t%s\n' "${SOURCE}" > "${MANIFEST}"

run_worker() {
  local mode="$1"
  H5AD_PREFLIGHT_MANIFEST="${MANIFEST}" \
  H5AD_PREFLIGHT_STATUS_DIR="${STATUS_DIR}" \
  H5AD_PREFLIGHT_RUN_ROOT="${RUN_ROOT}" \
  H5AD_PREFLIGHT_MODE="${mode}" \
  H5AD_PREFLIGHT_PYTHON_BIN="${ROOT}/.pixi/envs/default/bin/python" \
  H5AD_PREFLIGHT_TASK_ID=1 \
  SLURM_JOB_ID=999999 \
  SLURM_SUBMIT_DIR="${ROOT}" \
  bash "${ROOT}/src/utils/bash/h5ad_preflight_worker.sh"
}

run_worker classify
STATUS_FILE="${STATUS_DIR}/Adams__benchmark_analysis.status"
[[ "$(sed -n 's/^STATE=//p' "${STATUS_FILE}")" == "REBUILD" ]]

if run_worker require; then
  echo "require-mode accepted malformed H5AD" >&2
  exit 1
fi
[[ "$(sed -n 's/^STATE=//p' "${STATUS_FILE}")" == "FAIL" ]]

CAPTURE="${TMP_DIR}/sbatch.calls"
cat > "${TMP_DIR}/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
printf '812345\n'
STUB
chmod +x "${TMP_DIR}/sbatch"
export CAPTURE USER_EMAIL="test@example.invalid"
export ECODA_RUNTIME_MODE=host ECODA_RUNTIME_PROFILE=stage3
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/utils/bash/ecoda_runtime.sh"
source "${ROOT}/src/utils/bash/h5ad_preflight_submit.sh"
export PATH="${TMP_DIR}:${PATH}"
runtime_export="$(ecoda_runtime_export_csv stage3 0)"
preflight_id="$(ecoda_submit_h5ad_preflight \
  "${MANIFEST}" "${STATUS_DIR}" "${RUN_ROOT}" classify shared-cpu 1G 1 \
  "${RUN_ROOT}/logs" test "${ROOT}/src/utils/bash/h5ad_preflight_worker.sh" \
  "${runtime_export}")"
[[ "${preflight_id}" == "812345" ]]
case "$(cat "${CAPTURE}")" in
  *"--wait"*"--array=1-1%1"*) ;;
  *) echo "preflight scheduler array contract missing" >&2; exit 1 ;;
esac
case "$(cat "${CAPTURE}")" in
  *"H5AD_PREFLIGHT_MODE=classify"*) ;;
  *) echo "preflight mode was not exported" >&2; exit 1 ;;
esac

echo "h5ad preflight worker: OK"
