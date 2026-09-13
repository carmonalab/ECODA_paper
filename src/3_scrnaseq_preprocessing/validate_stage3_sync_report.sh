#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
exec "${PYTHON_BIN}" "${SCRIPT_DIR}/validate_stage3_sync_report.py" "$@"
