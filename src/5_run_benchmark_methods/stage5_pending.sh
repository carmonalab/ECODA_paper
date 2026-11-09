#!/bin/bash
# Pure residual selection: callbacks validate artifacts; only the report is written.
stage5_pending_fresh_corrected_derivative_valid() {
  local dispatch="$1" matrix="" actual_sha count
  # The historical repair lane never authorizes new excluded Python fits. A
  # fresh derivative wave has separate, immutable roots and exact scope.
  [[ -n "${ECODA_ARTIFACT_BINDING:-}" && "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 &&
     "${RUN_ID:-${ECODA_RUN_ID:-}}" == "${ECODA_RUN_ID:-}" &&
     "${ECODA_RECOVERY_ACTIVE:-0}" == 0 && -z "${SYNC_ONLY_RUN:-}" &&
     "${FORCE_ARG:-0}" == 0 && "${FORCE_TARGETED_ARG:-0}" == 0 &&
     "${TARGET_METHODS_SET:-0}" == 0 &&
     -z "${EXPECTED_PENDING_ARG:-}${ECODA_EXPECTED_PENDING:-}${PRECOMPUTE_BACKUP_RECEIPT_ARG:-}${ECODA_PRECOMPUTE_BACKUP_RECEIPT:-}" &&
     "${PASS_ARG:-}" == corrected && "${ANALYSIS_VARIANT:-${ANALYSIS_VARIANT_ARG:-}}" == corrected_final &&
     "${METHODS_ARG:-}" == prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot ]] || return 1
  command -v ecoda_derivative_stage5_validate >/dev/null 2>&1 || return 1
  ecoda_derivative_stage5_validate corrected corrected_final || return 1
  [[ "${ANALYSIS_ROOT}" == "${ECODA_RUN_ROOT}/artifacts/batch_effect/corrected_final" &&
     "${dispatch}" == "${ECODA_RUN_ROOT}/manifests/dispatch_selection.tsv" ]] || return 1
  case "${METHOD_MATRIX_MODE:-0}" in
    0) ;;
    1)
      matrix="${ECODA_RUN_ROOT}/manifests/method_matrix.tsv"
      [[ "${METHOD_MATRIX:-}" == "${matrix}" && "${ECODA_STAGE5_METHOD_MATRIX:-}" == "${matrix}" &&
         "${METHOD_MATRIX_SHA256:-}" =~ ^[0-9a-f]{64}$ &&
         "${METHOD_MATRIX_IDENTITY:-}" == "${METHOD_MATRIX_SHA256}" ]] || return 1
      ecoda_validate_run_owned_path "${matrix}" "${ECODA_RUN_ROOT}" &&
        ecoda_validate_manifest "${matrix}" 3 && ecoda_validate_checksum "${matrix}" || return 1
      actual_sha="$(ecoda_sha256_file "${matrix}")" || return 1
      count="$(awk 'END {print NR}' "${matrix}")" || return 1
      [[ "${actual_sha}" == "${METHOD_MATRIX_SHA256}" &&
         "${ECODA_CHECKSUM_MD5}" == "${METHOD_MATRIX_MD5:-}" &&
         "${ECODA_CHECKSUM_SIZE}" == "${METHOD_MATRIX_SIZE:-}" &&
         "${count}" == "${METHOD_MATRIX_DECLARED_COUNT:-}" && "${count}" -gt 0 ]] || return 1
      ;;
    *) return 1 ;;
  esac
  "${PYTHON_BIN}" - "${ECODA_ARTIFACT_BINDING}" "${dispatch}" "${DATASETS_JSON_FILE}" "${matrix}" <<'PY'
import itertools,json,pathlib,sys
binding=json.loads(pathlib.Path(sys.argv[1]).read_text());run=pathlib.Path(binding['run_root'])
config=json.loads(pathlib.Path(sys.argv[3]).read_text())
view='batch_effect_corrected'
permitted={d for d,v in config.items() if not d.startswith('_') and v.get('use_for_batch_effect') is True and view in v.get('views',{})}
methods=('prepare_pseudobulk','pseudobulk','gloscope','composition','mrvi','pilot','qot')
if binding.get('requested_selection') or binding['selection']!=str(run/'manifests/selection.tsv'):
    raise SystemExit('BLOCKED: fresh corrected derivative selection identity differs')
selection=[tuple(line.split('\t')) for line in pathlib.Path(binding['selection']).read_text().splitlines()]
datasets={row[0] for row in selection if len(row)==3}
rows=[tuple(line.split('\t')) for line in pathlib.Path(sys.argv[2]).read_text().splitlines()]
expected_rows={(d,view,m) for d,m in itertools.product(datasets,methods)}
if sys.argv[4]:
    matrix=[tuple(line.split('\t')) for line in pathlib.Path(sys.argv[4]).read_text().splitlines()]
    if (not matrix or len(set(matrix))!=len(matrix)
            or any(len(r)!=3 or not all(r) or r[1]!=view or r[2] not in methods for r in matrix)
            or {r[0] for r in matrix}!=datasets or rows!=matrix):
        raise SystemExit('BLOCKED: fresh corrected matrix dispatch differs from its exact declared scope')
    expected_rows=set(matrix)
if not datasets or not datasets<=permitted or len(selection)!=len(datasets) or any(len(r)!=3 or not all(r) or r[1]!=view or r[2]!=view for r in selection) or len(rows)!=len(expected_rows) or set(rows)!=expected_rows:
    raise SystemExit('BLOCKED: fresh corrected derivative dispatch differs from its exact selected cohort/family scope')
PY
}

stage5_classify_pending_rows() {
  local dispatch="$1" output="$2" ds view method extra path key pending="" needed="" seen=""
  local corrected=0
  local fresh_derivative=0
  [[ "${ANALYSIS_VARIANT:-${ANALYSIS_VARIANT_ARG:-}}" == corrected_final ]] && corrected=1
  [[ -f "${dispatch}" && ! -L "${dispatch}" && -r "${dispatch}" ]] || return 1
  if [[ ${corrected} -eq 1 ]] && stage5_pending_fresh_corrected_derivative_valid "${dispatch}"; then
    fresh_derivative=1
  fi
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
        stage5_pending_owner_valid "${ds}" "${view}" "${method}" "${fresh_derivative}" || return 1
      fi
      case "${method}" in
        mrvi|pilot|qot)
          [[ ${fresh_derivative} -eq 1 ]] || {
            echo "BLOCKED: excluded Python row ${key} is not strictly reusable" >&2; return 1;
          } ;;
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
          stage5_pending_owner_valid "${ds}" "${view}" "${method}" "${fresh_derivative}" || return 1
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
