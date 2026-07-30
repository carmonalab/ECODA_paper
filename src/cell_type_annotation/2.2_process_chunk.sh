#!/bin/bash
#
# 2.2_process_chunk.sh — Thin wrapper that calls 2.2_process_chunk.R via pixi
#

set -euo pipefail

TASK_ID="$1"
export PROJECT_ROOT="$2"
CHUNK_FILE="$3"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"

log_msg() {
  local msg="$1"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] [TASK ${TASK_ID}] ${msg}"
}

log_msg "============================================"
log_msg "Processing Chunk ID: ${TASK_ID}"
log_msg "============================================"

if [[ ! -f "${CHUNK_FILE}" ]]; then
  log_msg "ERROR: Chunk file ${CHUNK_FILE} not found."
  exit 1
fi

cd "${PROJECT_ROOT}"

"${HOME}/.pixi/bin/pixi" run Rscript --vanilla \
  "${SCRIPT_DIR}/2.2_process_chunk.R" \
  "chunk_file__${CHUNK_FILE}"

log_msg "Chunk ${TASK_ID} processing complete."
