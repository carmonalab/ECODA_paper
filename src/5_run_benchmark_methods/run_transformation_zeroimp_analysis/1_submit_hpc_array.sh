#!/bin/bash
# Compatibility entrypoint. Transformation and zero-imputation arrays are
# coordinated with benchmark methods by the root Pipeline 5 wrapper.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARGS=()
ANALYSIS_PRESENT=0
EXPECT_ANALYSIS_VALUE=0
for arg in "$@"; do
  if [[ ${EXPECT_ANALYSIS_VALUE} -eq 1 ]]; then
    ARGS+=("${arg}")
    EXPECT_ANALYSIS_VALUE=0
    continue
  fi
  case "${arg}" in
    --analysis)
      ARGS+=(--analyses)
      ANALYSIS_PRESENT=1
      EXPECT_ANALYSIS_VALUE=1
      ;;
    --analysis=*)
      ARGS+=(--analyses "${arg#*=}")
      ANALYSIS_PRESENT=1
      ;;
    --analyses|--analyses=*)
      ARGS+=("${arg}")
      ANALYSIS_PRESENT=1
      [[ "${arg}" == "--analyses" ]] && EXPECT_ANALYSIS_VALUE=1
      ;;
    *)
      ARGS+=("${arg}")
      ;;
  esac
done
if [[ ${ANALYSIS_PRESENT} -eq 0 ]]; then
  ARGS+=(--analyses trans,zeroimp)
fi
exec "${SCRIPT_DIR}/../1_submit_hpc_array.sh" "${ARGS[@]}"
