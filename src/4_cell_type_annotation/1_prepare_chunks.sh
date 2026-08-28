#!/bin/bash
set -euo pipefail
printf '%s\n' \
  "ERROR: legacy Stage 4 preparation entrypoint is retired; use" \
  "       1_submit_onboarding_stage.sh." >&2
exit 64
