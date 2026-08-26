#!/bin/bash
# Focused contract test for the durable stage-wise preprocessing submitter.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-preprocess-stage.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE

cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
CALL_COUNT="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
case "${CALL_COUNT}" in
  1) printf '600001\n' ;;
  2) printf '600002\n' ;;
  *) printf 'unexpected sbatch call count\n' >&2; exit 1 ;;
esac
STUB
chmod +x "${TMP_DIR}/bin/sbatch"

OUTPUT="$(
  HOME="${TMP_DIR}/home" \
  PATH="${TMP_DIR}/bin:${PATH}" \
  USER_EMAIL="test@example.invalid" \
  PREPROCESS_SUBMITTER_TEST=1 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_batch_effect_stage.sh"
)"

case "${OUTPUT}" in
  *"PREPROCESS_ARRAY_JOB_ID=600001"*) : ;;
  *) echo "missing array job marker" >&2; exit 1 ;;
esac
case "${OUTPUT}" in
  *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) : ;;
  *) echo "missing watchdog job marker" >&2; exit 1 ;;
esac

MANIFEST="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
[[ -s "${MANIFEST}" ]]
[[ "$(wc -l < "${MANIFEST}" | tr -d '[:space:]')" == "11" ]]
[[ "$(sed -n '1p' "${MANIFEST}")" == "Alzheimer" ]]
[[ "$(sed -n '9p' "${MANIFEST}")" == "Parkinson" ]]
[[ "$(sed -n '10p' "${MANIFEST}")" == "Joanito" ]]
[[ "$(sed -n '11p' "${MANIFEST}")" == "Stephenson" ]]

CALLS="$(cat "${CAPTURE}")"
case "${CALLS}" in
  *"--array=1-11%1000"*) : ;;
  *) echo "array was not submitted with all onboarding rows" >&2; exit 1 ;;
esac
case "${CALLS}" in
  *"--dependency=afterany:600001"*) : ;;
  *) echo "watchdog dependency was not bound to the array" >&2; exit 1 ;;
esac
case "${CALLS}" in
  *"PREPROCESS_DATASETS_FILE=${MANIFEST}"*) : ;;
  *) echo "dataset manifest was not exported to the array" >&2; exit 1 ;;
esac

echo "preprocessing stage submitter: OK"
