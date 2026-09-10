#!/bin/bash
# Execute the adjacent scGate R script without shell-reparsing its path or arguments.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
  RUNTIME_PREFIX="${ECODA_RUNTIME_PREFIX:-}"
else
  RUNTIME_PREFIX="${ECODA_HOST_ENV_PREFIX:-}"
fi
R_BIN="${RUNTIME_PREFIX%/}/bin/Rscript"
R_SCRIPT="${SCRIPT_DIR}/2.0_create_scgate_db.R"
[[ "${RUNTIME_PREFIX}" = /* && -f "${R_BIN}" && ! -L "${R_BIN}" && -x "${R_BIN}" ]] || {
  echo "ERROR: scGate Rscript is missing or unsafe: ${R_BIN}" >&2
  exit 1
}
[[ -f "${R_SCRIPT}" && ! -L "${R_SCRIPT}" && -r "${R_SCRIPT}" ]] || {
  echo "ERROR: scGate R script is missing or unsafe: ${R_SCRIPT}" >&2
  exit 1
}
exec "${R_BIN}" --vanilla "${R_SCRIPT}" "$@"
