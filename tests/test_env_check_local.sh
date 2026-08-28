#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "local default dependency check: skipped outside macOS"
  exit 0
fi
cd "${ROOT}"
pixi run check-r-deps
echo "local default dependency check: OK"
