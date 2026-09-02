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
  local count msg rc scheduler_id run_id

  [[ -r "${manifest}" && -s "${manifest}" ]] || return 1
  [[ -n "${runtime_export}" ]] || return 1
  [[ "${mode}" == "require" || "${mode}" == "classify" ]] || return 1
  run_id="${run_root##*/}"
  [[ "${run_id}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || return 1
  count="$(wc -l < "${manifest}" | tr -d '[:space:]')"
  [[ "${count}" =~ ^[0-9]+$ && ${count} -gt 0 ]] || return 1
  mkdir -p "${status_dir}" "${logs_dir}" || return 1

  if msg="$(sbatch --parsable --wait --array="1-${count}%${throttle}" \
      --partition="${partition}" --ntasks=1 --cpus-per-task=1 --mem="${memory}" \
      --time="${WATCHDOG_TIME_LIMIT:-12:00:00}" \
      --output="${logs_dir}/h5ad_preflight_${label}_%A_%a.log" \
      --error="${logs_dir}/h5ad_preflight_${label}_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      --export="ALL,H5AD_PREFLIGHT_RUN_ID=${run_id},H5AD_PREFLIGHT_MANIFEST=${manifest},H5AD_PREFLIGHT_STATUS_DIR=${status_dir},H5AD_PREFLIGHT_RUN_ROOT=${run_root},H5AD_PREFLIGHT_MODE=${mode},${runtime_export}" \
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

# A completed Slurm array can publish its run-owned status files a short time
# after sbatch --wait returns on the shared filesystem.  Wait only on those
# local files; this never polls scheduler state and never treats missing or
# malformed status as success.
ecoda_wait_h5ad_preflight_status_files() {
  local manifest="$1" status_dir="$2"
  local max_wait="${H5AD_PREFLIGHT_STATUS_GRACE_SECONDS:-60}"
  local elapsed=0 missing ds view path safe status
  [[ "${max_wait}" =~ ^[0-9]+$ ]] || return 1
  while :; do
    missing=0
    while IFS=$'\t' read -r ds view path; do
      safe="$(_ecoda_safe_component "${ds}__${view}")"
      status="${status_dir}/${safe}.status"
      [[ -s "${status}" ]] || missing=1
    done < "${manifest}"
    [[ ${missing} -eq 0 ]] && return 0
    [[ ${elapsed} -lt ${max_wait} ]] || return 1
    sleep 1
    elapsed=$((elapsed + 1))
  done
}
