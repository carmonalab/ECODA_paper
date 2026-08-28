#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda_env_lock.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT

ENV_LOCK_FILE="${TMP_DIR}/env_refresh.lock"
source "${ROOT}/src/utils/bash/env_mutation_lock.sh"

acquire_env_lock
[[ -d "${ENV_LOCK_FILE}" && -f "${ENV_LOCK_FILE}/owner" ]]
if (
  ENV_LOCK_FILE="${ENV_LOCK_FILE}"
  source "${ROOT}/src/utils/bash/env_mutation_lock.sh"
  acquire_env_lock
) >"${TMP_DIR}/contender.out" 2>"${TMP_DIR}/contender.err"; then
  echo "ERROR: concurrent lock acquisition unexpectedly succeeded." >&2
  exit 1
fi
release_env_lock
[[ ! -e "${ENV_LOCK_FILE}" ]]

mkdir "${ENV_LOCK_FILE}"
printf '999999 1\n' > "${ENV_LOCK_FILE}/owner"
if acquire_env_lock; then
  echo "ERROR: stale lock was reclaimed automatically." >&2
  exit 1
fi
rm -f "${ENV_LOCK_FILE}/owner"
rmdir "${ENV_LOCK_FILE}"

printf '%s %s\n' "$$" "$(date +%s)" > "${ENV_LOCK_FILE}"
if (
  ENV_LOCK_FILE="${ENV_LOCK_FILE}"
  source "${ROOT}/src/utils/bash/env_mutation_lock.sh"
  acquire_env_lock
) >"${TMP_DIR}/legacy.out" 2>"${TMP_DIR}/legacy.err"; then
  echo "ERROR: legacy active lock acquisition unexpectedly succeeded." >&2
  exit 1
fi
rm -f "${ENV_LOCK_FILE}"

RACE_LOCK="${TMP_DIR}/race.lock"
SUCCESS_FILE="${TMP_DIR}/race.success"
for _ in 1 2 3 4 5 6 7 8; do
  (
    ENV_LOCK_FILE="${RACE_LOCK}"
    source "${ROOT}/src/utils/bash/env_mutation_lock.sh"
    if acquire_env_lock; then
      printf 'winner\n' >> "${SUCCESS_FILE}"
      sleep 1
      release_env_lock
    fi
  ) &
done
wait
[[ "$(wc -l < "${SUCCESS_FILE}" | tr -d '[:space:]')" == 1 ]]
[[ ! -e "${RACE_LOCK}" ]]

echo "environment mutation lock: OK"
