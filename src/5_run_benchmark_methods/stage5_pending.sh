#!/bin/bash
# Pure residual selection: callbacks validate artifacts; only the report is written.
stage5_classify_pending_rows() {
  local dispatch="$1" output="$2" ds view method extra path key pending="" needed="" seen=""
  local corrected=0
  [[ "${ANALYSIS_VARIANT:-${ANALYSIS_VARIANT_ARG:-}}" == corrected_final ]] && corrected=1
  [[ -f "${dispatch}" && ! -L "${dispatch}" && -r "${dispatch}" ]] || return 1
  while IFS=$'\t' read -r ds view method extra; do
    [[ -n "${ds}" && -n "${view}" && -n "${method}" && -z "${extra}" ]] || return 1
    key="${ds}/${view}/${method}"
    case " ${seen} " in *" ${key} "*) return 1 ;; esac
    seen="${seen} ${key}"
    [[ "${method}" == prepare_pseudobulk && ${corrected} -eq 1 ]] && continue
    if ! stage5_method_is_forced "${method}" &&
       benchmark_selected_artifacts_valid "${ds}" "${view}" "${method}"; then
      continue
    fi
    if [[ ${corrected} -eq 1 ]]; then
      if command -v stage5_pending_owner_valid >/dev/null 2>&1; then
        stage5_pending_owner_valid "${ds}" "${view}" "${method}" || return 1
      fi
      case "${method}" in
        mrvi|pilot|qot) echo "BLOCKED: excluded Python row ${key} is not strictly reusable" >&2; return 1 ;;
      esac
      benchmark_artifacts_for "${ds}" "${view}" "${method}" || return 1
      for path in ${ARTIFACT_PATHS[@]+"${ARTIFACT_PATHS[@]}"}; do
        if command -v stage5_pending_target_valid >/dev/null 2>&1; then
          stage5_pending_target_valid "${path}" || return 1
        fi
        [[ ! -e "${path}" && ! -L "${path}" && ! -e "${path}.md5" && ! -L "${path}.md5" &&
           ! -e "${path}.runtime.json" && ! -L "${path}.runtime.json" &&
           ! -e "${path}.runtime.json.md5" && ! -L "${path}.runtime.json.md5" ]] || {
          echo "BLOCKED: existing invalid or unproven target ${path}" >&2; return 1;
        }
      done
    fi
    pending="${pending} ${key}"
    case "${method}" in
      pseudobulk|composition)
        if [[ ${corrected} -eq 1 ]] && ! stage5_prepare_pseudobulk_valid "${ds}" "${view}"; then
          needed="${needed} ${ds}/${view}/prepare_pseudobulk"
        fi ;;
    esac
  done < "${dispatch}"
  if [[ ${corrected} -eq 1 ]]; then
    for key in ${needed}; do
      case " ${seen} " in *" ${key} "*) ;; *) echo "BLOCKED: missing preparation dependency ${key}" >&2; return 1 ;; esac
      pending="${pending} ${key}"
    done
  fi
  : > "${output}" || return 1
  while IFS=$'\t' read -r ds view method extra; do
    key="${ds}/${view}/${method}"
    case " ${pending} " in
      *" ${key} "*)
        if command -v stage5_pending_owner_valid >/dev/null 2>&1; then
          stage5_pending_owner_valid "${ds}" "${view}" "${method}" || return 1
        fi
        if [[ ${corrected} -eq 1 && "${method}" == prepare_pseudobulk ]]; then
          benchmark_artifacts_for "${ds}" "${view}" "${method}" || return 1
          for path in ${ARTIFACT_PATHS[@]+"${ARTIFACT_PATHS[@]}"}; do
            if command -v stage5_pending_target_valid >/dev/null 2>&1; then
              stage5_pending_target_valid "${path}" || return 1
            fi
            [[ ! -e "${path}" && ! -L "${path}" && ! -e "${path}.md5" && ! -L "${path}.md5" ]] || {
              echo "BLOCKED: preparation dependency requires archived retirement ${path}" >&2; return 1;
            }
          done
        fi
        printf '%s\t%s\t%s\n' "${ds}" "${view}" "${method}" >> "${output}" || return 1 ;;
    esac
  done < "${dispatch}"
}

stage5_write_pending_checksum() {
  local path="$1"
  [[ -f "${path}" && ! -L "${path}" ]] || return 1
  if [[ -s "${path}" ]]; then ecoda_write_checksum "${path}"; return $?; fi
  ecoda_atomic_write "${path}.md5" "MD5=d41d8cd98f00b204e9800998ecf8427e\nSIZE=0\nPATH=${path}\n"
}

stage5_validate_pending_checksum() {
  local path="$1"
  [[ -f "${path}" && ! -L "${path}" ]] || return 1
  if [[ -s "${path}" ]]; then ecoda_validate_checksum "${path}"; return $?; fi
  [[ -f "${path}.md5" && ! -L "${path}.md5" ]] || return 1
  awk -F= -v p="${path}" '
    $1=="MD5" {m++; if($2!="d41d8cd98f00b204e9800998ecf8427e") bad=1; next}
    $1=="SIZE" {s++; if($2!="0") bad=1; next}
    $1=="PATH" {n++; if(substr($0,6)!=p) bad=1; next}
    {bad=1}
    END {exit(m==1 && s==1 && n==1 && !bad ? 0 : 1)}
  ' "${path}.md5"
}

stage5_validate_pending_binding() {
  local expected="$1" classified="$2"
  ecoda_validate_run_owned_path "${expected}" "${ECODA_RUN_ROOT}" || return 1
  stage5_validate_pending_checksum "${expected}" || return 1
  [[ ! -s "${expected}" ]] || ecoda_validate_manifest "${expected}" 3 || return 1
  cmp -s "${expected}" "${classified}" || {
    echo "BLOCKED: residual selection differs from expected pending rows" >&2; return 1;
  }
}
