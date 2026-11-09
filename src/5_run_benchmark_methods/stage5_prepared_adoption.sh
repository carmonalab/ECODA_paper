#!/bin/bash
# Copy only the finite, reviewed preparation set; ordinary dispatch/reuse checks
# remain authoritative. Source after the common helpers and dispatcher setup.
stage5_adopt_prepared_artifacts() {
  local receipt="$1" receipt_sha="$2" context="$3" helper plan bound_receipt schema
  local kind ds view original target producer md5 size owner path
  local cache_list="${ECODA_RUN_ROOT}/manifests/prepared_adoption_caches.tsv"
  local -a adopted_paths=()
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" && "${RUN_ID}" == "${ECODA_RUN_ID:-}" &&
     "${RUN_ID}" == "${ECODA_RUN_ROOT##*/}" && -z "${SYNC_ONLY_RUN:-}" &&
     "${FORCE_ARG:-0}" == 0 && "${FORCE_TARGETED_ARG:-0}" == 0 &&
     "${ECODA_RECOVERY_ACTIVE:-0}" == 0 && "${receipt}" == /* &&
     "${receipt_sha}" =~ ^[0-9a-f]{64}$ ]] || return 1
  helper="${PROJECT_ROOT}/src/utils/py/stage5_prepared_adoption.py"
  plan="${ECODA_RUN_ROOT}/manifests/prepared_adoption.tsv"
  bound_receipt="${ECODA_RUN_ROOT}/manifests/prepared_adoption.receipt.json"
  [[ ! -e "${plan}" && ! -L "${plan}" &&
     ! -e "${bound_receipt}" && ! -L "${bound_receipt}" ]] || return 1
  [[ "$(ecoda_sha256_file "${receipt}")" == "${receipt_sha}" ]] || return 1
  ecoda_atomic_copy "${receipt}" "${bound_receipt}" || return 1
  [[ "$(ecoda_sha256_file "${bound_receipt}")" == "${receipt_sha}" ]] || return 1
  schema="$("${PYTHON_BIN}" -c 'import json,sys; print(json.load(open(sys.argv[1]))["schema"])' "${bound_receipt}")" || return 1
  ecoda_write_checksum "${bound_receipt}" || return 1
  ecoda_write_artifact_record "${bound_receipt}" stage5_prepared_adoption_receipt "${RUN_ID}" >/dev/null || return 1
  "${PYTHON_BIN}" "${helper}" plan --receipt "${bound_receipt}" --sha256 "${receipt_sha}" \
    --config "${DATASETS_JSON_FILE}" --context "${context}" > "${plan}.build.$$" || return 1
  ecoda_atomic_install_manifest "${plan}.build.$$" "${plan}" 8 || return 1
  rm -f "${plan}.build.$$"
  ecoda_write_checksum "${plan}" || return 1
  # Reserve all fresh destinations before copying any payload. No existing
  # destination or terminal owner is reclaimed by this copy-only operation.
  while IFS=$'\t' read -r kind ds view original target producer md5 size; do
    owner="$(ecoda_artifact_owner_dir "${target}")" || return 1
    [[ ! -e "${owner}" && ! -L "${owner}" ]] || return 1
    owner="$(ecoda_artifact_owner_acquire "${target}" stage5 "${RUN_ID}" 0 0 0)" || return 1
    ecoda_owner_track "${owner}" || return 1
    adopted_paths+=("${target}")
  done < "${plan}"
  "${PYTHON_BIN}" "${helper}" copy --receipt "${bound_receipt}" --sha256 "${receipt_sha}" \
    --config "${DATASETS_JSON_FILE}" --context "${context}" || return 1
  : > "${cache_list}.build.$$" || return 1
  while IFS=$'\t' read -r kind ds view original target producer md5 size; do
    ecoda_validate_checksum "${target}" || return 1
    [[ "${ECODA_CHECKSUM_MD5}" == "${md5}" && "${ECODA_CHECKSUM_SIZE}" == "${size}" ]] || return 1
    ecoda_write_artifact_record "${target}" "${producer}" "${RUN_ID}" >/dev/null || return 1
    ecoda_validate_artifact_record "${target}" "${producer}" "${RUN_ID}" >/dev/null || return 1
    case "${kind}" in
      pseudobulk)
        printf '%s\tprepare_pseudobulk\t%s\t%s\t0\n' "${target}" "${ds}" "${view}" >> "${cache_list}.build.$$" || return 1
        ;;
      metadata)
        original="$(ecoda_derivative_input "${ds}" "${view}")" || return 1
        "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/export_h5ad_sample_metadata.py" \
          --config "${DATASETS_JSON_FILE}" --dataset "${ds}" --view "${view}" \
          --input-file "${original}" --output "${target}" --check || return 1
        ;;
      timing) ;;
      *) return 1 ;;
    esac
  done < "${plan}"
  if [[ -s "${cache_list}.build.$$" ]]; then
    ecoda_atomic_install_manifest "${cache_list}.build.$$" "${cache_list}" 5 || return 1
    rm -f "${cache_list}.build.$$"
    ecoda_write_checksum "${cache_list}" || return 1
    if [[ "${schema}" == ecoda_stage5_batch_prepared_adoption_v1 ]]; then
      [[ "${context}" == batch_effect/uncorrected_final ]] || return 1
      ${PIXI_RSCRIPT} "${PROJECT_ROOT}/src/5_run_benchmark_methods/validate_benchmark_rds_contract.R" \
        --artifact-list "${cache_list}" --config "${DATASETS_JSON_FILE}" \
        --input-root "${HPC_SCRATCH_DIR}" --source-identity "${SOURCE_IDENTITY}" \
        --source-identity-verified --batch-pass uncorrected --analysis-variant final || return 1
    else
    ${PIXI_RSCRIPT} "${PROJECT_ROOT}/src/5_run_benchmark_methods/validate_benchmark_rds_contract.R" \
      --artifact-list "${cache_list}" --config "${DATASETS_JSON_FILE}" \
      --input-root "${HPC_SCRATCH_DIR}" --source-identity "${SOURCE_IDENTITY}" \
      --source-identity-verified || return 1
    fi
  fi
  rm -f "${cache_list}.build.$$"
  path="${ECODA_RUN_ROOT}/manifests/prepared_adoption.json"
  ecoda_write_checksum "${path}" || return 1
  ecoda_write_artifact_record "${path}" stage5_prepared_adoption "${RUN_ID}" >/dev/null || return 1
  # Only a fully validated copy set becomes reusable. Its immutable report
  # preserves original fit/timing identity and each original source gate state.
  for path in "${adopted_paths[@]}"; do
    ecoda_artifact_owner_set_state "${path}" OK 'validated prepared artifact byte-copy; no scientific refit' || return 1
  done
}
