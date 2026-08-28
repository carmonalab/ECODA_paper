#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-scgate.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
trap 'rm -rf "${TMP_DIR}"' EXIT
VALID="${TMP_DIR}/valid.rds"
INVALID="${TMP_DIR}/invalid.rds"
pixi run Rscript --vanilla -e 'model <- data.frame(levels="T", use_as="PBMC", name="model", signature="T"); saveRDS(list(human=list(PBMC=list(model=model),HiTME=list(model=model))),commandArgs(TRUE)[1]); saveRDS(list(human=list(PBMC=list(model=model))),commandArgs(TRUE)[2])' "${VALID}" "${INVALID}"
PROJECT_ROOT="${ROOT}" SCGATE_DB_PATH="${VALID}" \
  pixi run Rscript --vanilla "${ROOT}/src/4_cell_type_annotation/2.0_create_scgate_db.R" --validate-only >/dev/null
set +e
INVALID_OUTPUT="$(PROJECT_ROOT="${ROOT}" SCGATE_DB_PATH="${INVALID}" pixi run Rscript --vanilla "${ROOT}/src/4_cell_type_annotation/2.0_create_scgate_db.R" --validate-only 2>&1)"
RC=$?
set -e
[[ ${RC} -ne 0 ]]
case "${INVALID_OUTPUT}" in *"requires non-empty human\$PBMC and human\$HiTME"*) ;; *) exit 1 ;; esac
echo "scGate cache contract: OK"
