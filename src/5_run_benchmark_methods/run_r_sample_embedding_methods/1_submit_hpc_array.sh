#!/bin/bash
# Compatibility entrypoint. The coordinated root wrapper owns all R/Python
# arrays, documented pseudobulk dependencies, and one final synchronization.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METHODS_PRESENT=0
for arg in "$@"; do
  case "${arg}" in --methods|--methods=*) METHODS_PRESENT=1 ;; esac
done
ARGS=("$@")
if [[ ${METHODS_PRESENT} -eq 0 ]]; then
  ARGS+=(--methods gloscope,mofa,pseudobulk,scitd)
fi
exec "${SCRIPT_DIR}/../1_submit_hpc_array.sh" "${ARGS[@]}"
