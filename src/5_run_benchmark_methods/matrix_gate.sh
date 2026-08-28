#!/bin/bash
# Aggregate benchmark gate: one terminal status after every child watchdog.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
[[ $# -eq 3 ]] || { echo "Usage: matrix_gate.sh RUN_ROOT WATCHDOG_LABELS_CSV SCHEDULER_IDS_MANIFEST" >&2; exit 2; }
RUN_ROOT="$1"
WATCHDOG_LABELS_CSV="$2"
SCHEDULER_IDS_MANIFEST="$3"
STATUS_FILE="${RUN_ROOT}/status/aggregate"
mkdir -p "$(dirname "${STATUS_FILE}")"

fail() {
  local tmp="${STATUS_FILE}.tmp.$$"
  printf 'STATE=FAIL\nREASON=%s\nWATCHDOG_LABELS=%s\n' "$1" "${WATCHDOG_LABELS_CSV}" > "${tmp}"
  mv -f "${tmp}" "${STATUS_FILE}"
  exit 1
}

[[ -r "${SCHEDULER_IDS_MANIFEST}" ]] || fail "missing scheduler ID manifest"
ecoda_validate_run_owned_path "${SCHEDULER_IDS_MANIFEST}" "${RUN_ROOT}" ||
  fail "scheduler ID manifest is outside run root"
ecoda_validate_manifest "${SCHEDULER_IDS_MANIFEST}" 2 || fail "malformed scheduler ID manifest"
IFS=',' read -r -a labels <<< "${WATCHDOG_LABELS_CSV}"
SCHEDULER_IDS=()
while IFS=$'\t' read -r kind scheduler_id; do
  [[ -n "${scheduler_id}" ]] || fail "blank scheduler ID in manifest"
  [[ "${scheduler_id}" =~ ^[0-9]+$ ]] || fail "malformed scheduler ID ${scheduler_id}"
  SCHEDULER_IDS+=("${scheduler_id}")
done < "${SCHEDULER_IDS_MANIFEST}"
for label in "${labels[@]}"; do
  [[ -n "${label}" ]] || continue
  safe="$(printf '%s' "${label}" | tr '/:,\t ' '_____')"
  status="${RUN_ROOT}/status/watchdogs/${safe}.status"
  [[ -s "${status}" ]] || fail "missing watchdog status ${label}"
  grep -q '^STATE=OK$' "${status}" || fail "watchdog ${label} did not report OK"
  while IFS= read -r line; do
    case "${line}" in
      SCHEDULER_ID=*)
        scheduler_id="${line#*=}"
        [[ "${scheduler_id}" =~ ^[0-9]+$ ]] || fail "malformed watchdog scheduler ID ${scheduler_id}"
        case " ${SCHEDULER_IDS[*]} " in
          *" ${scheduler_id} "*) ;;
          *) SCHEDULER_IDS+=("${scheduler_id}") ;;
        esac
        ;;
    esac
  done < "${status}"
done
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  [[ "${SLURM_JOB_ID}" =~ ^[0-9]+$ ]] ||
    fail "malformed aggregate gate scheduler ID ${SLURM_JOB_ID}"
  case " ${SCHEDULER_IDS[*]} " in
    *" ${SLURM_JOB_ID} "*) ;;
    *) SCHEDULER_IDS+=("${SLURM_JOB_ID}") ;;
  esac
fi
if ! partials="$(find "${RUN_ROOT}" -type f \( -name '*.tmp.*' -o -name '*.partial' -o -name '*.build.*' \) -print)"; then
  fail "unable to inspect benchmark partial artifacts"
fi
[[ -z "${partials}" ]] || fail "partial benchmark artifacts remain under ${RUN_ROOT}"
tmp="${STATUS_FILE}.tmp.$$"
{
  printf 'STATE=OK\nWATCHDOG_LABELS=%s\n' "${WATCHDOG_LABELS_CSV}"
  for scheduler_id in "${SCHEDULER_IDS[@]}"; do
    printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"
  done
} > "${tmp}"
mv -f "${tmp}" "${STATUS_FILE}"
