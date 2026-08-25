#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SUBMITTER="${ROOT_DIR}/src/5_run_benchmark_methods/run_batch_effect_methods/1_submit_hpc_array.sh"

bash -n "${SUBMITTER}"
bash -n "${ROOT_DIR}/src/5_run_benchmark_methods/benchmark_submit_common.sh"
bash -n "${ROOT_DIR}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
bash -n "${ROOT_DIR}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"

set +e
CORRECTED_OUTPUT="$(bash "${SUBMITTER}" --pass corrected --ds_name Alzheimer 2>&1)"
CORRECTED_RC=$?
set -e
if [[ ${CORRECTED_RC} -eq 0 ]]; then
  echo "corrected null-batch guard unexpectedly submitted" >&2
  exit 1
fi
case "${CORRECTED_OUTPUT}" in
  *"corrected batch-effect view requires a confirmed columns.batch"*) ;;
  *)
    echo "unexpected corrected guard output: ${CORRECTED_OUTPUT}" >&2
    exit 1
    ;;
esac

set +e
ALLOWLIST_OUTPUT="$(bash "${SUBMITTER}" --pass uncorrected --ds_name _debug --methods scpoli 2>&1)"
ALLOWLIST_RC=$?
set -e
if [[ ${ALLOWLIST_RC} -eq 0 ]]; then
  echo "disallowed method unexpectedly submitted" >&2
  exit 1
fi
case "${ALLOWLIST_OUTPUT}" in
  *"Unknown batch-effect method 'scpoli'"*) ;;
  *)
    echo "unexpected allow-list output: ${ALLOWLIST_OUTPUT}" >&2
    exit 1
    ;;
esac

printf '%s\n' "batch-effect submitter guards OK"
