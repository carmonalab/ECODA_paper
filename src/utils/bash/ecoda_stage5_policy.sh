#!/bin/bash
# Stage 5 recovery, root, method, and artifact-path policy.
# Source after ecoda_run_common.sh.


ecoda_stage5_validate_identity() {
  local requested_pass="${1:-}" pass variant="${2:-${ANALYSIS_VARIANT:-}}"
  local scratch_root nas_root expected_root expected_nas expected_suffix
  [[ -n "${requested_pass}" ]] || requested_pass="${PASS_ARG:-${ANALYSIS_PASS:-}}"
  pass="${requested_pass}"
  if [[ -n "${PASS_ARG:-}" && "${PASS_ARG}" != "${pass}" ]]; then
    _ecoda_die "Stage 5 pass identity disagrees with PASS_ARG"
    return 1
  fi
  if [[ -n "${ANALYSIS_PASS:-}" && "${ANALYSIS_PASS}" != "${pass}" ]]; then
    _ecoda_die "Stage 5 pass identity disagrees with ANALYSIS_PASS"
    return 1
  fi
  case "${variant}" in
    "")
      return 0
      ;;
    final)
      expected_suffix="uncorrected_final"
      [[ "${pass}" == "uncorrected" ]] || {
        _ecoda_die "final Stage 5 analysis requires the uncorrected pass"
        return 1
      }
      ;;
    corrected_final)
      expected_suffix="corrected_final/recovery_35row"
      [[ "${pass}" == "corrected" ]] || {
        _ecoda_die "corrected_final Stage 5 analysis requires the corrected pass"
        return 1
      }
      ;;
    *)
      _ecoda_die "unsupported Stage 5 analysis variant: ${variant}"
      return 1
      ;;
  esac
  scratch_root="${HPC_SCRATCH_DIR:-}"
  nas_root="${NAS_TARGET_DIR:-}"
  [[ "${scratch_root}" = /* && "${scratch_root}" != *$'\n'* &&
     "${scratch_root}" != *$'\t'* &&
     "${nas_root}" = /* && "${nas_root}" != *$'\n'* &&
     "${nas_root}" != *$'\t'* ]] || {
    _ecoda_die "final Stage 5 analysis requires absolute scratch and NAS roots"
    return 1
  }
  scratch_root="${scratch_root%/}"
  nas_root="${nas_root%/}"
  [[ -n "${scratch_root}" ]] || scratch_root="/"
  [[ -n "${nas_root}" ]] || nas_root="/"
  if [[ "${scratch_root}" == "/" ]]; then
    expected_root="/batch_effect/${expected_suffix}"
  else
    expected_root="${scratch_root}/batch_effect/${expected_suffix}"
  fi
  if [[ "${nas_root}" == "/" ]]; then
    expected_nas="/batch_effect/${expected_suffix}"
  else
    expected_nas="${nas_root}/batch_effect/${expected_suffix}"
  fi
  [[ "${ANALYSIS_ROOT:-}" == "${expected_root}" ]] || {
    _ecoda_die "Stage 5 analysis root must be ${expected_root}"
    return 1
  }
  [[ "${ANALYSIS_NAS_ROOT:-}" == "${expected_nas}" ]] || {
    _ecoda_die "Stage 5 NAS root must be ${expected_nas}"
    return 1
  }
}
ecoda_stage5_validate_final_method() {
  local method="${1:-}"
  [[ -n "${ANALYSIS_VARIANT:-}" ]] || return 0
  case "${method}" in
    prepare_pseudobulk|pseudobulk|gloscope|composition|mrvi|pilot|qot) ;;
    *)
      _ecoda_die "method is not permitted in the final Stage 5 suite: ${method}"
      return 1
      ;;
  esac
}

ecoda_stage5_method_matrix_allows() {
  local ds="${1:-}" view="${2:-}" method="${3:-}"
  local matrix="${ECODA_STAGE5_METHOD_MATRIX:-${METHOD_MATRIX:-}}"
  local row_ds row_view row_method extra
  # Non-matrix runs retain unrestricted method expansion.  A bound matrix is
  # authoritative only when explicitly set.
  [[ -f "${matrix}" && ! -L "${matrix}" && -r "${matrix}" ]] || return 1
  if [[ -n "${ECODA_RUN_ROOT:-}" ]]; then
    ecoda_validate_run_owned_path "${matrix}" "${ECODA_RUN_ROOT}" || return 1
  fi
  while IFS=$'\t' read -r row_ds row_view row_method extra; do
    [[ -z "${extra}" ]] || return 1
    if [[ "${row_ds}" == "${ds}" && "${row_view}" == "${view}" &&
          "${row_method}" == "${method}" ]]; then
      return 0
    fi
  done < "${matrix}"
  return 1
}

ecoda_stage5_method_matrix_path() {
  local matrix="${ECODA_STAGE5_METHOD_MATRIX:-${METHOD_MATRIX:-}}"
  [[ -n "${matrix}" ]] || return 1
  printf '%s' "${matrix}"
}

ecoda_stage5_analysis_root_suffix() {
  local variant="${1:-${ANALYSIS_VARIANT:-}}"
  case "${variant}" in
    final)
      printf 'uncorrected_final'
      ;;
    corrected_final)
      printf 'corrected_final/recovery_35row'
      ;;
    *)
      return 1
      ;;
  esac
}

ecoda_stage5_validate_artifact_path() {
  local path="${1:-}" canonical scratch_root nas_root expected_root expected_nas
  local variant="${ANALYSIS_VARIANT:-}" expected_suffix
  [[ -n "${variant}" ]] || return 0
  ecoda_stage5_validate_identity || return 1
  expected_suffix="$(ecoda_stage5_analysis_root_suffix "${variant}")" || return 1
  [[ -n "${path}" ]] || {
    _ecoda_die "final Stage 5 artifact path is empty"
    return 1
  }
  canonical="$(_ecoda_canonical_path "${path}")" || {
    _ecoda_die "final Stage 5 artifact path cannot be canonicalized: ${path}"
    return 1
  }
  scratch_root="${HPC_SCRATCH_DIR%/}"
  nas_root="${NAS_TARGET_DIR%/}"
  [[ -n "${scratch_root}" ]] || scratch_root="/"
  [[ -n "${nas_root}" ]] || nas_root="/"
  if [[ "${scratch_root}" == "/" ]]; then
    expected_root="/batch_effect/${expected_suffix}"
  else
    expected_root="${scratch_root}/batch_effect/${expected_suffix}"
  fi
  if [[ "${nas_root}" == "/" ]]; then
    expected_nas="/batch_effect/${expected_suffix}"
  else
    expected_nas="${nas_root}/batch_effect/${expected_suffix}"
  fi
  case "${canonical}" in
    "${expected_root}"/*|"${expected_nas}"/*) ;;
    *)
      _ecoda_die "Stage 5 artifact escapes ${expected_suffix} roots: ${path}"
      return 1
      ;;
  esac
}

ecoda_stage5_batch_stem() {
  local ds="${1:-}" pass="${2:-${PASS_ARG:-${ANALYSIS_PASS:-}}}"
  local variant="${3:-${ANALYSIS_VARIANT:-}}" stem
  [[ -n "${ds}" && "${ds}" != *$'\n'* && "${ds}" != *$'\t'* ]] || return 1
  if [[ -n "${pass}" ]]; then
    stem="${ds}_batch_effect_${pass}"
    case "${variant}" in
      "") ;;
      final)
        [[ "${pass}" == "uncorrected" ]] || return 1
        [[ -z "${PASS_ARG:-}" || "${PASS_ARG}" == "${pass}" ]] || return 1
        [[ -z "${ANALYSIS_PASS:-}" || "${ANALYSIS_PASS}" == "${pass}" ]] || return 1
        stem="${stem}_final"
        ;;
      corrected_final)
        [[ "${pass}" == "corrected" ]] || return 1
        [[ -z "${PASS_ARG:-}" || "${PASS_ARG}" == "${pass}" ]] || return 1
        [[ -z "${ANALYSIS_PASS:-}" || "${ANALYSIS_PASS}" == "${pass}" ]] || return 1
        stem="${stem}_final"
        ;;
      *) return 1 ;;
    esac
  else
    [[ -z "${variant}" ]] || return 1
    stem="${ds}"
  fi
  printf '%s' "${stem}"
}

_ecoda_stage5_artifacts_for() {
  local ds="$1" view="$2" label="$3" pass="${PASS_ARG:-${ANALYSIS_PASS:-}}"
  local root nas_root stem batch_stem suffix n
  # A missing pass is ordinary mode even when a caller supplies a batch view.
  # Batch callers bind PASS_ARG/ANALYSIS_PASS before ownership expansion.
  ecoda_stage5_validate_identity "${pass}" || return 1
  if [[ -n "${ECODA_STAGE5_METHOD_MATRIX:-${METHOD_MATRIX:-}}" ]]; then
    ecoda_stage5_method_matrix_allows "${ds}" "${view}" "${label}" || {
      _ecoda_die "Stage 5 method is not authorized by the bound method matrix: ${ds}/${view}/${label}"
      return 1
    }
  fi
  case "${ANALYSIS_VARIANT:-}" in
    final)
      [[ "${pass}" == "uncorrected" &&
         "${view}" == "batch_effect_uncorrected" ]] || {
        _ecoda_die "final Stage 5 artifacts require the uncorrected batch-effect view"
        return 1
      }
      ;;
    corrected_final)
      [[ "${pass}" == "corrected" &&
         "${view}" == "batch_effect_corrected" ]] || {
        _ecoda_die "corrected_final Stage 5 artifacts require the corrected batch-effect view"
        return 1
      }
      ;;
  esac
  ecoda_stage5_validate_final_method "${label}" || return 1
  if [[ -n "${pass}" ]]; then
    batch_stem="$(ecoda_stage5_batch_stem "${ds}" "${pass}" "${ANALYSIS_VARIANT:-}")" ||
      return 1
  else
    batch_stem="${ds}"
  fi
  ECODA_BENCHMARK_ARTIFACTS=()
  if [[ -n "${ANALYSIS_ROOT:-}" ]]; then
    root="${ANALYSIS_ROOT}"
  elif [[ -n "${pass}" ]]; then
    if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
      root="${HPC_SCRATCH_DIR}/batch_effect/$(ecoda_stage5_analysis_root_suffix)"
    else
      root="${HPC_SCRATCH_DIR}/batch_effect/${pass}"
    fi
  else
    root="${HPC_SCRATCH_DIR}/benchmark"
  fi
  if [[ -n "${ANALYSIS_NAS_ROOT:-}" ]]; then
    nas_root="${ANALYSIS_NAS_ROOT}"
  elif [[ -n "${pass}" && -n "${NAS_TARGET_DIR:-}" ]]; then
    if [[ -n "${ANALYSIS_VARIANT:-}" ]]; then
      nas_root="${NAS_TARGET_DIR}/batch_effect/$(ecoda_stage5_analysis_root_suffix)"
    else
      nas_root="${NAS_TARGET_DIR}/batch_effect/${pass}"
    fi
  elif [[ -n "${NAS_TARGET_DIR:-}" ]]; then
    nas_root="${NAS_TARGET_DIR}/benchmark"
  else
    nas_root=""
  fi
  if [[ "${label}" == "prepare_pseudobulk" ]]; then
    if [[ -n "${pass}" ]]; then
      ECODA_BENCHMARK_ARTIFACTS+=("${root}/pseudobulks/${batch_stem}_pseudobulk_hvg2000.rds")
    else
      for stem in schvg2000 hvg2000 hvg500 hvg2000_bl hvg1000 hvg3000; do
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/pseudobulks/${ds}_pseudobulk_${stem}.rds")
      done
    fi
  else
    case "${label}" in
      mrvi)
        if [[ -n "${pass}" ]]; then
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${batch_stem}_hvg2000_highres_mrvi_dists.feather")
        else
          for n in 1000 2000 3000; do
            ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_mrvi_dists.feather")
          done
        fi
        ;;
      scpoli)
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_lowres_scpoli_dims15_embs.feather")
        for n in 1000 3000; do
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_highres_scpoli_dims15_embs.feather")
        done
        for n in 2 3 5 10 15; do
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_highres_scpoli_dims${n}_embs.feather")
        done
        ;;
      pilot|qot)
        suffix="${label}"
        if [[ -n "${pass}" ]]; then
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${batch_stem}_hvg2000_highres_${suffix}_dists.feather")
        else
          ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_lowres_${suffix}_dists.feather")
          for n in 1000 2000 3000; do
            ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg${n}_highres_${suffix}_dists.feather")
          done
        fi
        ;;
      pilotgm)
        [[ -z "${pass}" ]] || return 1
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/embeddings/${ds}_hvg2000_highres_pilotgm_dists.feather")
        ;;
      trans|zeroimp)
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${ds}_${label}.rds")
        ;;
      gloscope|mofa|pseudobulk|scitd)
        stem="${batch_stem}"
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_${label}.rds")
        ;;
      composition)
        stem="${batch_stem}"
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_composition.rds")
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_metadata.rds")
        ;;
      *)
        # Explicit post-baseline methods use one method-specific result key.
        # The submitter's method registry remains authoritative for whether
        # the row is runnable; this fallback gives the ownership layer a
        # deterministic path without broad filesystem discovery.
        stem="${batch_stem}"
        ECODA_BENCHMARK_ARTIFACTS+=("${root}/results/${stem}_${label}.rds")
        ;;
    esac
  fi
  ECODA_BENCHMARK_ARTIFACT_NAS=()
  if [[ -n "${nas_root}" ]]; then
    for stem in "${ECODA_BENCHMARK_ARTIFACTS[@]}"; do
      case "${stem}" in
        "${root}"/*) ECODA_BENCHMARK_ARTIFACT_NAS+=("${nas_root}/${stem#${root}/}") ;;
        *) return 1 ;;
      esac
    done
  fi
}
ecoda_validate_corrected_batch_columns() {
  local config_path="${1:-${DATASETS_JSON_FILE:-}}"
  local dataset="${2:-}"
  local view="${3:-batch_effect_corrected}"
  [[ -r "${config_path}" ]] || {
    _ecoda_die "corrected batch configuration is unreadable: ${config_path}"
    return 1
  }
  [[ -n "${dataset}" && -n "${view}" ]] || {
    _ecoda_die "corrected batch validation requires a dataset and view"
    return 1
  }
  if ! jq -e --arg ds "${dataset}" --arg view "${view}" '
    def nonblank_string:
      if type != "string" then false
      elif length == 0 then false
      else test("[^[:space:]]")
      end;

    if type != "object" then
      error("datasets configuration must be a JSON object")
    elif (has($ds) | not) then
      error("dataset is not configured")
    else
      .[$ds] as $entry |
      if ($entry | type) != "object" then
        error("dataset entry must be a JSON object")
      elif (($entry.views // null) | type) != "object" then
        error("dataset views must be a JSON object")
      elif (($entry.views | has($view)) | not) then
        error("corrected view is not configured")
      else
        $entry.views[$view] as $view_spec |
        if ($view_spec | type) != "object" then
          error("corrected view must be a JSON object")
        elif (($entry.columns // null) | type) != "object" then
          error("dataset columns must be a JSON object")
        else
          ($entry.columns) as $columns |
          ($columns.batch) as $raw_batch |
          (if ($raw_batch | type) == "string" then
             [$raw_batch]
           elif ($raw_batch | type) == "array" then
             $raw_batch
           else
             error("columns.batch must be a string or nonempty array of strings")
           end) as $keys |
          if ($keys | length) == 0 then
            error("columns.batch must not be empty")
          elif any($keys[]; (nonblank_string | not)) then
            error("columns.batch contains an empty, blank, or non-string key")
          elif (($keys | unique | length) != ($keys | length)) then
            error("columns.batch contains duplicate keys")
          elif any($keys[];
                   . == "Sample" or . == "__ecoda_batch_combined_v1") then
            error("columns.batch contains a reserved key")
          elif (($columns.label? != null) and
                (($columns.label | type) == "string") and
                any($keys[]; . == $columns.label)) then
            error("columns.batch overlaps the configured label column")
          elif (($columns.sample? != null) and
                (($columns.sample | type) == "string") and
                any($keys[]; . == $columns.sample)) then
            error("columns.batch overlaps the configured sample column")
          else
            true
          end
        end
      end
    end
  ' "${config_path}" >/dev/null; then
    _ecoda_die "invalid corrected columns.batch for dataset ${dataset}"
    return 1
  fi

}
ecoda_corrected_batch_method_policy() {
  local method="${1:-}"
  ECODA_CORRECTED_BATCH_METHOD_ID=""
  ECODA_CORRECTED_BATCH_MODEL_ID=""
  case "${method}" in
    preprocess)
      ECODA_CORRECTED_BATCH_METHOD_ID="preprocess"
      ECODA_CORRECTED_BATCH_MODEL_ID="hvg_composite_v1"
      ;;
    prepare_pseudobulk|pseudobulk)
      ECODA_CORRECTED_BATCH_METHOD_ID="Pseudobulk"
      ECODA_CORRECTED_BATCH_MODEL_ID="pseudobulk_limma_fixed_effects_v1"
      ;;
    gloscope)
      ECODA_CORRECTED_BATCH_METHOD_ID="GloScope"
      ECODA_CORRECTED_BATCH_MODEL_ID="embedding_consumer_harmony_v1"
      ;;
    pilot)
      ECODA_CORRECTED_BATCH_METHOD_ID="PILOT"
      ECODA_CORRECTED_BATCH_MODEL_ID="embedding_consumer_harmony_v1"
      ;;
    qot)
      ECODA_CORRECTED_BATCH_METHOD_ID="QOT"
      ECODA_CORRECTED_BATCH_MODEL_ID="embedding_consumer_harmony_v1"
      ;;
    mrvi)
      ECODA_CORRECTED_BATCH_METHOD_ID="MrVI"
      ECODA_CORRECTED_BATCH_MODEL_ID="mrvi_composite_v1"
      ;;
    composition)
      ECODA_CORRECTED_BATCH_METHOD_ID="ECODA_authors_HR"
      ECODA_CORRECTED_BATCH_MODEL_ID="limma_fixed_effects_v1"
      ;;
    *)
      _ecoda_die "unsupported corrected batch method policy: ${method}"
      return 1
      ;;
  esac
}

_ecoda_stage5_expand_output_selection() {
  local selection="${1:-}" ds view label extra path nas_path artifact_index
  _ecoda_init_output_arrays
  [[ -r "${selection}" && ! -L "${selection}" ]] || {
    _ecoda_die "Stage 5 selection manifest is missing or unreadable: ${selection}"
    return 1
  }
  [[ -r "${DATASETS_JSON_FILE:-}" ]] || {
    _ecoda_die "configured datasets.json is missing for Stage 5 ownership expansion"
    return 1
  }
  ecoda_validate_manifest "${selection}" 3 || return 1
  ECODA_OUTPUT_PATHS=()
  ECODA_OUTPUT_WRITE_FLAGS=()
  ECODA_OUTPUT_OWNER_DIRS=()
  while IFS=$'\t' read -r ds view label extra; do
    [[ -n "${ds}" && -n "${view}" && -n "${label}" && -z "${extra}" &&
       "${label}" =~ ^[A-Za-z0-9_.-]+$ ]] || return 1
    ecoda_dataset_exists "${ds}" && ecoda_view_exists "${ds}" "${view}" || return 1
    _ecoda_stage5_artifacts_for "${ds}" "${view}" "${label}" || return 1
    artifact_index=0
    for path in "${ECODA_BENCHMARK_ARTIFACTS[@]}"; do
      nas_path=""
      if [[ ${#ECODA_BENCHMARK_ARTIFACT_NAS[@]} -gt 0 ]]; then
        nas_path="${ECODA_BENCHMARK_ARTIFACT_NAS[${artifact_index}]}"
      fi
      artifact_index=$((artifact_index + 1))
      if [[ -n "${nas_path}" ]]; then
        _ecoda_output_add_scratch_nas_pair "${path}" "${nas_path}" || return 1
      else
        _ecoda_output_add_path "${path}" || return 1
      fi
    done
  done < "${selection}"
  [[ ${#ECODA_OUTPUT_PATHS[@]} -gt 0 ]] || {
    _ecoda_die "Stage 5 selection expands to no output artifacts: ${selection}"
    return 1
  }
}

ecoda_stage5_validate_output_ownership() {
  local selection="${1:-}" run_id="${2:-}" reclaim_terminal="${3:-0}"
  [[ "${selection}" = /* ]] || {
    _ecoda_die "Stage 5 ownership selection must be an absolute path"
    return 1
  }
  _ecoda_stage5_expand_output_selection "${selection}" || return 1
  _ecoda_validate_expanded_output_ownership \
    stage5 "${run_id}" "${reclaim_terminal}"
}

ecoda_stage5_artifact_owner_acquire() {
  ecoda_stage5_validate_artifact_path "${1:-}" || return 1
  ecoda_artifact_owner_acquire "$@"
}

ecoda_stage5_artifact_owner_validate() {
  ecoda_stage5_validate_artifact_path "${1:-}" || return 1
  ecoda_artifact_owner_validate "$@"
}