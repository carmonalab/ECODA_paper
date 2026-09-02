#!/bin/bash
# Compute-node watchdog for selected Stage 2 hooks. It owns terminal
# accounting, OOM-only retries, output validation, and owner finalization.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
cd "${PROJECT_ROOT}"

if [[ $# -ne 7 ]]; then
  echo "Usage: stage2_watchdog.sh RUN_ID STEPS_MANIFEST JOBS_FILE MEM MAX_MEM PARTITION THROTTLE" >&2
  exit 2
fi
RUN_ID="$1"
MANIFEST="$2"
JOB_FILE="$3"
CURRENT_MEMORY="$4"
MAX_MEMORY="$5"
PARTITION="$6"
THROTTLE="$7"
ecoda_validate_run_id "${RUN_ID}" || exit 1
ecoda_open_run "${RUN_ID}" || exit 1
RUN_ROOT="${ECODA_RUN_ROOT}"
STATUS_FILE="${RUN_ROOT}/status/watchdog"
RETRY_INDEX=0
SCHEDULER_IDS=()
RUNTIME_EXPORT=""

stage2_step_script() {
  case "$1" in
    gongsharma_cap) printf '%s/1.1_submit_gongsharma.sh' "${SCRIPT_DIR}" ;;
    combinedpbmc) printf '%s/1.2_submit_combinedpbmc.sh' "${SCRIPT_DIR}" ;;
    joanito) printf '%s/1.3_submit_joanito.sh' "${SCRIPT_DIR}" ;;
    kfoury_lowres_ct) printf '%s/1.4_submit_kfoury_lowres_ct.sh' "${SCRIPT_DIR}" ;;
    myocardial_counts) printf '%s/1.5_submit_myocardial.sh' "${SCRIPT_DIR}" ;;
    bassez_cellsubtype) printf '%s/1.6_submit_bassez.sh' "${SCRIPT_DIR}" ;;
    *) return 1 ;;
  esac
}

stage2_step_outputs() {
  case "$1" in
    gongsharma_cap) printf '%s;%s' \
      "${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data/SoundLife_YoungAdult_Male_CMVneg.h5ad" \
      "${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data/SoundLife_YoungAdult_Male_CMVpos.h5ad" ;;
    combinedpbmc) printf '%s/CombinedPBMC/data/combined_pbmc.h5ad' "${HPC_SCRATCH_DIR}" ;;
    joanito) printf '%s/Joanito/data/%s;%s/_debug/data/JoaI_2022_35773407_debug_5samples.h5ad' \
      "${HPC_SCRATCH_DIR}" "$(ecoda_view_input_name Joanito batch_effect_uncorrected)" "${HPC_SCRATCH_DIR}" ;;
    kfoury_lowres_ct) printf '%s/Kfoury/data/Kfoury_2021_34719426.rds' "${HPC_SCRATCH_DIR}" ;;
    myocardial_counts) printf '%s/Myocardial_infarction/data/Myocardial_Infarc_2.h5ad' "${HPC_SCRATCH_DIR}" ;;
    bassez_cellsubtype) printf '%s/Bassez/data/BassezA_2021_33958794whole.rds' "${HPC_SCRATCH_DIR}" ;;
    *) return 1 ;;
  esac
}

atomic_status() {
  local state="$1" reason="${2:-}"
  local tmp="${STATUS_FILE}.tmp.$$"
  mkdir -p "$(dirname "${STATUS_FILE}")"
  {
    printf 'STATE=%s\nRUN_ID=%s\nREASON=%s\n' "${state}" "${RUN_ID}" "${reason}"
    printf 'RETRY_INDEX=%s\n' "${RETRY_INDEX}"
    if [[ -r "${JOB_FILE}" ]]; then
      while IFS=$'\t' read -r step job; do
        [[ -n "${step}" ]] && printf 'JOB=%s:%s\n' "${step}" "${job}"
      done < "${JOB_FILE}"
    fi
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
      printf 'SCHEDULER_ID=%s\n' "${SLURM_JOB_ID}"
    fi
    if [[ ${#SCHEDULER_IDS[@]} -gt 0 ]]; then
      local scheduler_id
      for scheduler_id in "${SCHEDULER_IDS[@]}"; do
        printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"
      done
    fi
  } > "${tmp}" || return 1
  mv -f "${tmp}" "${STATUS_FILE}"
}

fail() {
  local reason="$1"
  local owner_rc=0
  if [[ -r "${MANIFEST}" ]]; then
    while IFS=$'\t' read -r step script outputs dependency owner; do
      [[ -n "${owner}" && "${owner}" != "-" ]] || continue
      if ! ecoda_owner_set_state "${owner}" FAIL "${reason}"; then
        owner_rc=1
      fi
    done < "${MANIFEST}"
  fi
  if ! atomic_status FAIL "${reason}"; then
    owner_rc=1
  fi
  exit 1
}
export ECODA_RUNTIME_PROFILE=stage2
ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE}" || \
  fail "Stage 2 immutable runtime validation failed before retry handling"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage2 0)" || \
  fail "Stage 2 runtime export construction failed"

bump_memory() {
  local value="$1" number suffix
  if [[ "${value}" =~ ^([0-9]+)([GT])$ ]]; then
    number="${BASH_REMATCH[1]}"
    suffix="${BASH_REMATCH[2]}"
    printf '%s%s' "$((number * 2))" "${suffix}"
  else
    return 1
  fi
}

mem_ge() {
  local a="$1" b="$2" an as bn bs
  [[ "${a}" =~ ^([0-9]+)([GT])$ ]] || return 1
  an="${BASH_REMATCH[1]}"; as="${BASH_REMATCH[2]}"
  [[ "${b}" =~ ^([0-9]+)([GT])$ ]] || return 1
  bn="${BASH_REMATCH[1]}"; bs="${BASH_REMATCH[2]}"
  if [[ "${as}" == "T" ]]; then an=$((an * 1024)); fi
  if [[ "${bs}" == "T" ]]; then bn=$((bn * 1024)); fi
  (( an >= bn ))
}

wait_job_terminal() {
  local job="$1"
  ecoda_wait_scalar_accounting "${job}" "${STAGE2_WATCHDOG_POLL_SECONDS:-30}" || return 1
  printf '%s' "${ECODA_ACCOUNTING_STATE}"
}

validate_manifests() {
  ecoda_validate_run_owned_path "${MANIFEST}" "${RUN_ROOT}" || return 1
  ecoda_validate_manifest "${MANIFEST}" 5 || return 1
  ecoda_validate_run_owned_path "${JOB_FILE}" "${RUN_ROOT}" || return 1
  ecoda_validate_manifest "${JOB_FILE}" 2 || return 1
  local seen_steps="" step script outputs dependency owner expected_script expected_outputs
  while IFS=$'\t' read -r step script outputs dependency owner; do
    expected_script="$(stage2_step_script "${step}" 2>/dev/null || true)"
    expected_outputs="$(stage2_step_outputs "${step}" 2>/dev/null || true)"
    [[ -n "${expected_script}" && "${script}" == "${expected_script}" ]] || return 1
    [[ -n "${expected_outputs}" && "${outputs}" == "${expected_outputs}" ]] || return 1
    expected_dependency="-"
    [[ "${step}" == "combinedpbmc" ]] && expected_dependency="gongsharma_cap"
    [[ "${dependency}" == "${expected_dependency}" ]] || return 1
    expected_owner="$(ecoda_owner_dir stage2 "${step}")"
    [[ "${owner}" == "${expected_owner}" ]] || return 1
    case " ${seen_steps} " in *" ${step} "*) return 1 ;; esac
    seen_steps="${seen_steps} ${step}"
  done < "${MANIFEST}"
  local seen_jobs="" job_step job_id known_owner
  while IFS=$'\t' read -r job_step job_id; do
    [[ "${job_id}" =~ ^[0-9]+$ ]] || return 1
    case " ${seen_jobs} " in *" ${job_step} "*) return 1 ;; esac
    known_owner="$(awk -F '\t' -v step="${job_step}" '$1 == step {print $5}' "${MANIFEST}")"
    [[ -n "${known_owner}" && "${known_owner}" != "-" ]] || return 1
    seen_jobs="${seen_jobs} ${job_step}"
  done < "${JOB_FILE}"
  [[ -n "${seen_jobs}" ]] || return 1
}

validate_outputs() {
  local path step script outputs dependency owner old_ifs sidecar sidecar_present
  while IFS=$'\t' read -r step script outputs dependency owner; do
    [[ -n "${step}" ]] || return 1
    old_ifs="${IFS}"
    IFS=';'
    read -r -a paths <<< "${outputs}"
    IFS="${old_ifs}"
    for path in "${paths[@]}"; do
      [[ -s "${path}" ]] || {
        echo "Missing/empty Stage 2 output: ${path}" >&2
        return 1
      }
      sidecar="${path}.md5"
      sidecar_present=0
      if [[ -e "${sidecar}" || -L "${sidecar}" ]]; then
        # Existing checksums are immutable evidence; an invalid sidecar is
        # never replaced by a semantically plausible output.
        ecoda_validate_checksum "${path}" || return 1
        sidecar_present=1
      fi
      ecoda_validate_stage2_output "${step}" "${path}" || {
        echo "Stage 2 semantic output validation failed: ${path}" >&2
        return 1
      }
      case "${path}" in
        *.h5ad)
          "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/artifact_contract.py" \
            --path "${path}" --kind h5ad >/dev/null 2>&1 || return 1
          ;;
      esac
      if [[ ${sidecar_present} -eq 0 ]]; then
        # A genuinely new worker output has no prior checksum to validate.
        # Its semantic contracts must pass before the first sidecar is written.
        ecoda_write_checksum "${path}" ||
          return 1
        ecoda_validate_checksum_record "${path}" "${ECODA_CHECKSUM_MD5}" \
          "${ECODA_CHECKSUM_SIZE}" || return 1
      fi
      if [[ "${step}" == "combinedpbmc" ]]; then
        rm -f "${HPC_SCRATCH_DIR}/CombinedPBMC/data/combined_pbmc_batch_effect_analysis.h5ad" \
          "${HPC_SCRATCH_DIR}/CombinedPBMC/data/combined_pbmc_batch_effect_analysis.h5ad.md5"
      fi
    done
  done < "${MANIFEST}"
}

validate_manifests || fail "Stage 2 manifest or job contract failed"
while IFS=$'\t' read -r initial_step initial_job; do
  SCHEDULER_IDS+=("${initial_job}")
done < "${JOB_FILE}"

while :; do
  OOM_STEPS=()
  FAILED_STEPS=()
  CURRENT_JOBS=()
  while IFS=$'\t' read -r step job; do
    [[ -n "${step}" ]] || continue
    CURRENT_JOBS+=("${step}:${job}")
    state="$(wait_job_terminal "${job}")" || fail "unable to obtain terminal state for ${step}/${job}"
    case "${state}" in
      COMPLETED) ;;
      OUT_OF_MEMORY) OOM_STEPS+=("${step}") ;;
      *) FAILED_STEPS+=("${step}:${state}") ;;
    esac
  done < "${JOB_FILE}"

  if [[ ${#FAILED_STEPS[@]} -gt 0 ]]; then
    # A dependent CombinedPBMC job can terminate as dependency-never-satisfied
    # when the cap OOMs. It is retried with the cap, not treated as a
    # non-OOM scientific failure.
    CAP_OOM=0
    for step in "${OOM_STEPS[@]}"; do [[ "${step}" == gongsharma_cap ]] && CAP_OOM=1; done
    if [[ ${CAP_OOM} -eq 1 ]]; then
      RETAINED_FAILURES=()
      for failure in "${FAILED_STEPS[@]}"; do
        failed_step="${failure%%:*}"
        failed_state="${failure#*:}"
        if [[ "${failed_step}" == combinedpbmc && ( "${failed_state}" == DEPENDENCY* || "${failed_state}" == CANCELLED ) ]]; then
          OOM_STEPS+=(combinedpbmc)
        else
          RETAINED_FAILURES+=("${failure}")
        fi
      done
      FAILED_STEPS=("${RETAINED_FAILURES[@]}")
    fi
  fi
  if [[ ${#FAILED_STEPS[@]} -gt 0 ]]; then
    fail "non-OOM Stage 2 hook failure: ${FAILED_STEPS[*]}"
  fi
  if [[ ${#OOM_STEPS[@]} -eq 0 ]]; then
    break
  fi
  if mem_ge "${CURRENT_MEMORY}" "${MAX_MEMORY}"; then
    fail "OUT_OF_MEMORY Stage 2 hooks at ${MAX_MEMORY} ceiling: ${OOM_STEPS[*]}"
  fi
  NEXT_MEMORY="$(bump_memory "${CURRENT_MEMORY}")" || fail "unparseable Stage 2 memory '${CURRENT_MEMORY}'"
  if mem_ge "${NEXT_MEMORY}" "${MAX_MEMORY}"; then NEXT_MEMORY="${MAX_MEMORY}"; fi
  RETRY_INDEX=$((RETRY_INDEX + 1))
  [[ ${RETRY_INDEX} -le 4 ]] || fail "exceeded Stage 2 OOM retry attempts"

  RETRY_JOB_FILE="${RUN_ROOT}/manifests/jobs.retry_${RETRY_INDEX}.tsv"
  RETRY_JOB_TMP="${RETRY_JOB_FILE}.build.$$"
  : > "${RETRY_JOB_TMP}"
  RETRY_CAP_JOB=""
  for step in "${OOM_STEPS[@]}"; do
    script=""; dependency=""
    while IFS=$'\t' read -r mstep mscript outputs mdependency owner; do
      if [[ "${mstep}" == "${step}" ]]; then script="${mscript}"; dependency="${mdependency}"; break; fi
    done < "${MANIFEST}"
    [[ -n "${script}" ]] || fail "OOM step missing from manifest: ${step}"
    retry_export="ALL,STAGE2_RUN_ROOT=${RUN_ROOT},FORCE_PREPROCESS=${STAGE2_FORCE:-0},${RUNTIME_EXPORT}"
    retry_args=(--parsable --partition="${PARTITION}" --mem="${NEXT_MEMORY}" \
      --output="${LOGS_DIR}/stage2_${step}_retry${RETRY_INDEX}_%j.log" \
      --error="${LOGS_DIR}/stage2_${step}_retry${RETRY_INDEX}_%j.err" \
      --mail-user="${USER_EMAIL}" --export="${retry_export}")
    if [[ "${step}" == "combinedpbmc" && -n "${dependency}" ]]; then
      dep_job="${RETRY_CAP_JOB}"
      if [[ -z "${dep_job}" ]]; then
        while IFS=$'\t' read -r dep_step dep_id; do
          [[ "${dep_step}" == "${dependency}" ]] && dep_job="${dep_id}"
        done < "${JOB_FILE}"
      fi
      [[ -n "${dep_job}" ]] && retry_args+=(--dependency="afterok:${dep_job}")
    fi
    set +e
    retry_output="$(sbatch "${retry_args[@]}" "${script}")"
    retry_rc=$?
    set -e
    [[ ${retry_rc} -eq 0 ]] || fail "sbatch rejected Stage 2 OOM retry for ${step}"
    retry_id="${retry_output%%;*}"
    [[ "${retry_id}" =~ ^[0-9]+$ ]] || fail "invalid retry id for ${step}: ${retry_id}"
    printf '%s\t%s\n' "${step}" "${retry_id}" >> "${RETRY_JOB_TMP}"
    SCHEDULER_IDS+=("${retry_id}")
    [[ "${step}" == gongsharma_cap ]] && RETRY_CAP_JOB="${retry_id}"
    echo "STAGE2_RETRY_JOB_ID=${step}:${retry_id}"
  done
  if ! ecoda_atomic_install_manifest "${RETRY_JOB_TMP}" "${RETRY_JOB_FILE}" 2; then
    fail "failed to install Stage 2 retry scheduler manifest atomically"
  fi
  rm -f "${RETRY_JOB_TMP}"
  ecoda_validate_run_owned_path "${RETRY_JOB_FILE}" "${RUN_ROOT}" ||
    fail "Stage 2 retry manifest escaped the run root"
  ecoda_validate_manifest "${RETRY_JOB_FILE}" 2 ||
    fail "Stage 2 retry scheduler manifest is invalid"
  JOB_FILE="${RETRY_JOB_FILE}"
  CURRENT_MEMORY="${NEXT_MEMORY}"
  if ! ecoda_atomic_install_manifest "${JOB_FILE}" "${RUN_ROOT}/manifests/jobs.tsv" 2; then
    fail "failed to install Stage 2 scheduler manifest atomically"
  fi
done

validate_outputs || fail "Stage 2 output contract/checksum validation failed"
owner_failure=0
while IFS=$'\t' read -r step script outputs dependency owner; do
  [[ -n "${owner}" && "${owner}" != "-" ]] || fail "Stage 2 owner is missing for ${step}"
  if ! ecoda_owner_set_state "${owner}" OK "validated by Stage 2 watchdog"; then
    owner_failure=1
  fi
done < "${MANIFEST}"
[[ ${owner_failure} -eq 0 ]] || fail "failed to finalize Stage 2 owners"
if ! atomic_status OK "all selected hooks completed and outputs validated"; then
  fail "failed to write Stage 2 watchdog success status"
fi
printf 'Stage 2 watchdog completed for run %s\n' "${RUN_ID}"
