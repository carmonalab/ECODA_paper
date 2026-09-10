#!/bin/bash
# Aggregate benchmark gate: one terminal status after every child watchdog.
set -euo pipefail

[[ $# -eq 3 ]] || {
  echo "Usage: matrix_gate.sh RUN_ROOT WATCHDOG_LABELS_CSV SCHEDULER_IDS_MANIFEST" >&2
  exit 2
}
RUN_ROOT="$1"
WATCHDOG_LABELS_CSV="$2"
SCHEDULER_IDS_MANIFEST="$3"
SNAPSHOT_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"

case "${SNAPSHOT_REQUIRED}" in
  1)
    SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
    SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
    BOUND_RUN_ROOT="${ECODA_RUN_ROOT:-}"
    BOUND_RUN_ID="${ECODA_RUN_ID:-}"
    [[ "${SOURCE_ROOT}" = /* && -d "${SOURCE_ROOT}" && ! -L "${SOURCE_ROOT}" ]] || {
      echo "ERROR: Stage 5 snapshot aggregate gate requires an absolute immutable source root." >&2
      exit 1
    }
    [[ "${SOURCE_MANIFEST}" = /* && -f "${SOURCE_MANIFEST}" &&
       ! -L "${SOURCE_MANIFEST}" && -r "${SOURCE_MANIFEST}" && -s "${SOURCE_MANIFEST}" ]] || {
      echo "ERROR: Stage 5 snapshot aggregate gate requires an absolute immutable source manifest." >&2
      exit 1
    }
    GATE_SCRIPT="${SOURCE_ROOT%/}/src/5_run_benchmark_methods/matrix_gate.sh"
    [[ -f "${GATE_SCRIPT}" && ! -L "${GATE_SCRIPT}" && -r "${GATE_SCRIPT}" ]] || {
      echo "ERROR: immutable Stage 5 aggregate gate script is missing or unsafe." >&2
      exit 1
    }
    for helper in \
      "${SOURCE_ROOT%/}/src/slurm_config.sh" \
      "${SOURCE_ROOT%/}/src/utils/bash/ecoda_runtime.sh" \
      "${SOURCE_ROOT%/}/src/utils/bash/ecoda_run_common.sh"; do
      [[ -f "${helper}" && ! -L "${helper}" && -r "${helper}" ]] || {
        echo "ERROR: immutable Stage 5 shared helper is missing or unsafe: ${helper}" >&2
        exit 1
      }
    done
    export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
    export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
    export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
    source "${SOURCE_ROOT%/}/src/slurm_config.sh"
    source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_runtime.sh"
    source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_run_common.sh"
    export ECODA_RUN_ROOT="${BOUND_RUN_ROOT}"
    export ECODA_RUN_ID="${BOUND_RUN_ID}"
    ;;
  0|"")
    echo "ERROR: legacy_source_unpinned: Stage 5 aggregate gate requires an immutable source snapshot." >&2
    exit 1
    ;;
  *)
    echo "ERROR: ECODA_SOURCE_SNAPSHOT_REQUIRED must be 1 for snapshot-backed execution." >&2
    exit 1
    ;;
esac

STATUS_FILE="${RUN_ROOT}/status/aggregate"
if [[ "${SNAPSHOT_REQUIRED}" == "1" ]]; then
  [[ "${RUN_ROOT}" = /* && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] || {
    echo "ERROR: Stage 5 snapshot aggregate gate run root is missing or unsafe." >&2
    exit 1
  }
  [[ "${ECODA_RUN_ROOT:-}" = "${RUN_ROOT}" && "${ECODA_RUN_ROOT}" = /* ]] || {
    echo "ERROR: Stage 5 aggregate gate run root does not match ECODA_RUN_ROOT." >&2
    exit 1
  }
fi
mkdir -p "$(dirname "${STATUS_FILE}")"

fail() {
  local tmp="${STATUS_FILE}.tmp.$$"
  printf 'STATE=FAIL\nREASON=%s\nWATCHDOG_LABELS=%s\n' "$1" "${WATCHDOG_LABELS_CSV}" > "${tmp}"
  mv -f "${tmp}" "${STATUS_FILE}"
  exit 1
}

if [[ "${SNAPSHOT_REQUIRED}" == "1" ]]; then
  [[ "${ECODA_RUN_ID:-}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] ||
    fail "Stage 5 aggregate gate requires an exact ECODA_RUN_ID"
  [[ "${ECODA_RUN_ROOT##*/}" == "${ECODA_RUN_ID}" ]] ||
    fail "Stage 5 aggregate gate run root does not match ECODA_RUN_ID"
  ecoda_validate_run_id "${ECODA_RUN_ID}" ||
    fail "Stage 5 aggregate gate run ID is invalid"
  EXPECTED_RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${ECODA_RUN_ID}"
  [[ "${ECODA_RUN_ROOT}" == "${EXPECTED_RUN_ROOT}" ]] ||
    fail "Stage 5 aggregate gate run root is not the exact bound run root"

  RUN_SOURCE_MANIFEST="${ECODA_RUN_ROOT}/manifests/source.manifest"
  RUN_RUNTIME_IDENTITY="${ECODA_RUN_ROOT}/manifests/runtime.identity"
  [[ -f "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" &&
     -r "${RUN_SOURCE_MANIFEST}" && -s "${RUN_SOURCE_MANIFEST}" ]] ||
    fail "legacy_source_unpinned: Stage 5 run-bound source.manifest is missing or unsafe"
  [[ -f "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" &&
     -r "${RUN_RUNTIME_IDENTITY}" && -s "${RUN_RUNTIME_IDENTITY}" ]] ||
    fail "legacy_source_unpinned: Stage 5 run-bound runtime.identity is missing or unsafe"
  ecoda_validate_run_owned_path "${RUN_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}" ||
    fail "Stage 5 run-bound source.manifest escaped the run root"
  ecoda_validate_run_owned_path "${RUN_RUNTIME_IDENTITY}" "${ECODA_RUN_ROOT}" ||
    fail "Stage 5 run-bound runtime.identity escaped the run root"
  cmp -s "${RUN_SOURCE_MANIFEST}" "${ECODA_SOURCE_MANIFEST}" ||
    fail "Stage 5 run source manifest differs from the immutable source manifest"
  if [[ -n "${ECODA_RUNTIME_IDENTITY:-}" &&
        "${ECODA_RUNTIME_IDENTITY}" != "${RUN_RUNTIME_IDENTITY}" ]]; then
    fail "Stage 5 runtime identity does not match the run-owned runtime.identity"
  fi
  export ECODA_RUN_ROOT ECODA_RUN_ID
  export ECODA_RUNTIME_IDENTITY="${RUN_RUNTIME_IDENTITY}"
  export ECODA_RUNTIME_PROFILE=stage5
  LOGS_DIR="${ECODA_RUN_ROOT}/logs"
  export LOGS_DIR ECODA_LOGS_DIR="${LOGS_DIR}"
  ecoda_runtime_validate_bound_run ||
    fail "Stage 5 run-bound runtime validation failed before aggregate status"
fi

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
