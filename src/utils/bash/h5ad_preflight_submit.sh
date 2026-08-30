#!/bin/bash
# Shared scheduler submission for H5AD contract preflight arrays.
# Call after slurm_config.sh and ecoda_run_common.sh are sourced.

ecoda_submit_h5ad_preflight() {
  local manifest="$1"
  local status_dir="$2"
  local run_root="$3"
  local mode="$4"
  local partition="$5"
  local memory="$6"
  local throttle="$7"
  local logs_dir="$8"
  local label="$9"
  local worker_script="${10}"
  local runtime_export="${11:-}"
  local count msg rc scheduler_id

  [[ -r "${manifest}" && -s "${manifest}" ]] || return 1
  [[ -n "${runtime_export}" ]] || return 1
  [[ "${mode}" == "require" || "${mode}" == "classify" ]] || return 1
  count="$(wc -l < "${manifest}" | tr -d '[:space:]')"
  [[ "${count}" =~ ^[0-9]+$ && ${count} -gt 0 ]] || return 1
  mkdir -p "${status_dir}" "${logs_dir}" || return 1

  if msg="$(sbatch --parsable --wait --array="1-${count}%${throttle}" \
      --partition="${partition}" --ntasks=1 --cpus-per-task=1 --mem="${memory}" \
      --time="${WATCHDOG_TIME_LIMIT:-12:00:00}" \
      --output="${logs_dir}/h5ad_preflight_${label}_%A_%a.log" \
      --error="${logs_dir}/h5ad_preflight_${label}_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      --export="ALL,H5AD_PREFLIGHT_MANIFEST=${manifest},H5AD_PREFLIGHT_STATUS_DIR=${status_dir},H5AD_PREFLIGHT_RUN_ROOT=${run_root},H5AD_PREFLIGHT_MODE=${mode},${runtime_export}" \
      "${worker_script}")"; then
    rc=0
  else
    rc=$?
  fi
  scheduler_id="${msg%%;*}"
  if [[ "${scheduler_id}" =~ ^[0-9]+$ ]]; then
    printf '%s\n' "${scheduler_id}"
  else
    return 1
  fi
  return "${rc}"
}
