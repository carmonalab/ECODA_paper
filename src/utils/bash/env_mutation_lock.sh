#!/bin/bash
# Atomic lock for the shared Pixi/R environment.
# Source after slurm_config.sh with ENV_LOCK_FILE set.
# The lock is a directory so mkdir is the ownership operation.  A lock left
# by SIGKILL or a lost node is intentionally not reclaimed automatically:
# callers must verify that no writer is active before removing it manually.

acquire_env_lock() {
  if [[ -e "${ENV_LOCK_FILE}" ]]; then
    echo "ERROR: environment lock already exists: ${ENV_LOCK_FILE}" >&2
    echo "       refusing to reclaim it automatically; verify no environment writer is active, then remove the lock directory." >&2
    return 1
  fi
  if ! mkdir "${ENV_LOCK_FILE}" 2>/dev/null; then
    echo "ERROR: environment lock was claimed concurrently: ${ENV_LOCK_FILE}" >&2
    return 1
  fi
  if ! printf '%s %s\n' "$$" "$(date +%s)" > "${ENV_LOCK_FILE}/owner"; then
    rm -f "${ENV_LOCK_FILE}/owner"
    rmdir "${ENV_LOCK_FILE}" 2>/dev/null || true
    echo "ERROR: could not write environment lock metadata: ${ENV_LOCK_FILE}" >&2
    return 1
  fi
  ENV_LOCK_BACKEND="mkdir"
  trap 'release_env_lock' EXIT
}

release_env_lock() {
  local lock_pid lock_ts
  [[ "${ENV_LOCK_BACKEND:-}" == mkdir ]] || return 0
  [[ -f "${ENV_LOCK_FILE}/owner" ]] || return 0
  read -r lock_pid lock_ts < "${ENV_LOCK_FILE}/owner" || return 0
  [[ "${lock_pid}" == "$$" ]] || return 0
  rm -f "${ENV_LOCK_FILE}/owner"
  rmdir "${ENV_LOCK_FILE}" 2>/dev/null || true
}
