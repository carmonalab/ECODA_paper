#!/bin/bash
# Compatibility entrypoint for explicit batch-effect passes. The canonical
# root wrapper keeps pass-qualified roots and one final sync owner.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PASS=""
METHODS_PRESENT=0
ARGS=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pass) PASS="${2:-}"; shift 2 ;;
    --pass=*) PASS="${1#*=}"; shift ;;
    --methods|--methods=*) METHODS_PRESENT=1; shift ;;
    *) shift ;;
  esac
done
[[ "${PASS}" == uncorrected || "${PASS}" == corrected ]] || { echo "ERROR: --pass must be uncorrected or corrected" >&2; exit 1; }
for arg in "${ARGS[@]}"; do
  case "${arg}" in *scpoli*) echo "Unknown batch-effect method 'scpoli'" >&2; exit 1 ;; esac
done
if [[ ${METHODS_PRESENT} -eq 0 ]]; then
  ARGS+=(--methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot)
fi
exec "${SCRIPT_DIR}/../1_submit_hpc_array.sh" "${ARGS[@]}"
