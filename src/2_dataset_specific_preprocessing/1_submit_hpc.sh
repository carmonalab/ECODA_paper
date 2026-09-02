#!/bin/bash
# Dataset-specific preprocessing gate. Independent hooks are submitted in one
# wave; only GongSharma cap -> CombinedPBMC is serialized.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
export ECODA_GATE_STAGE=stage2
source "${SCRIPT_DIR}/../utils/bash/sync_status_email.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"

DATASETS_ARG=""
DATASETS_SET=0
STEPS_ARG=""
STEPS_SET=0
FORCE_ARG=0
SYNC_ONLY_RUN=""
SYNC_ONLY_SET=0
MEMORY="${STAGE2_MEM:-128G}"
MAX_MEMORY="${STAGE2_MEM_MAX:-500G}"
PARTITION="${SLURM_PARTITION}"
RUNTIME_EXPORT=""
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"

usage() {
  cat <<'EOF'
Usage: 1_submit_hpc.sh [--datasets LIST] [--steps LIST] [--force]
       [--sync-only RUN_ID] [--mem VALUE] [--max-mem VALUE]
       [--partition NAME] [--throttle N]

Steps: gongsharma_cap, combinedpbmc, joanito, kfoury_lowres_ct,
       myocardial_counts, bassez_cellsubtype
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets) DATASETS_ARG="${2:-}"; DATASETS_SET=1; shift 2 ;;
    --datasets=*) DATASETS_ARG="${1#*=}"; DATASETS_SET=1; shift ;;
    --ds_name) DATASETS_ARG="${2:-}"; DATASETS_SET=1; shift 2 ;;
    --ds_name=*) DATASETS_ARG="${1#*=}"; DATASETS_SET=1; shift ;;
    --steps) STEPS_ARG="${2:-}"; STEPS_SET=1; shift 2 ;;
    --steps=*) STEPS_ARG="${1#*=}"; STEPS_SET=1; shift ;;
    --force) FORCE_ARG=1; shift ;;
    --sync-only) SYNC_ONLY_RUN="${2:-}"; SYNC_ONLY_SET=1; shift 2 ;;
    --sync-only=*) SYNC_ONLY_RUN="${1#*=}"; SYNC_ONLY_SET=1; shift ;;
    --mem) MEMORY="${2:-}"; shift 2 ;;
    --mem=*) MEMORY="${1#*=}"; shift ;;
    --max-mem) MAX_MEMORY="${2:-}"; shift 2 ;;
    --max-mem=*) MAX_MEMORY="${1#*=}"; shift ;;
    --partition) PARTITION="${2:-}"; shift 2 ;;
    --partition=*) PARTITION="${1#*=}"; shift ;;
    --throttle) THROTTLE="${2:-}"; shift 2 ;;
    --throttle=*) THROTTLE="${1#*=}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -n "${SYNC_ONLY_RUN}" && ${FORCE_ARG} -eq 1 ]]; then
  echo "ERROR: --sync-only cannot be combined with --force." >&2
  exit 1
fi
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required for Stage 2 selection." >&2
  exit 1
fi

if [[ ${DATASETS_SET} -eq 1 && -z "${DATASETS_ARG}" ]]; then
  echo "ERROR: --datasets must not be empty." >&2
  exit 1
fi
if [[ ${STEPS_SET} -eq 1 && -z "${STEPS_ARG}" ]]; then
  echo "ERROR: --steps must not be empty." >&2
  exit 1
fi
if [[ ${SYNC_ONLY_SET} -eq 1 && -z "${SYNC_ONLY_RUN}" ]]; then
  echo "ERROR: --sync-only requires a run ID." >&2
  exit 1
fi

stage2_abort() {
  local reason="$1"
  local cleanup_rc=0
  if [[ -n "${ECODA_RUN_ROOT:-}" ]]; then
    ecoda_owner_finalize_tracked FAIL "${reason}" || cleanup_rc=1
    ecoda_set_run_state FAIL "${reason}" || cleanup_rc=1
  fi
  echo "ERROR: ${reason}" >&2
  exit 1
}

step_script() {
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

step_outputs() {
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

if [[ -n "${SYNC_ONLY_RUN}" ]]; then
  ecoda_open_run "${SYNC_ONLY_RUN}" || exit 1
  MANIFEST="${ECODA_RUN_ROOT}/manifests/steps.tsv"
  JOB_FILE="${ECODA_RUN_ROOT}/manifests/jobs.tsv"
  ecoda_validate_run_owned_path "${MANIFEST}" "${ECODA_RUN_ROOT}" ||
    stage2_abort "Stage 2 steps manifest is not run-owned"
  ecoda_validate_manifest "${MANIFEST}" 5 ||
    stage2_abort "Stage 2 steps manifest is invalid"
  seen_steps=""
  pending_count=0
  failed=0
  while IFS=$'\t' read -r step script outputs dependency owner; do
    expected_script="$(step_script "${step}" 2>/dev/null || true)"
    expected_outputs="$(step_outputs "${step}" 2>/dev/null || true)"
    case " ${seen_steps} " in
      *" ${step} "*) echo "ERROR: duplicate Stage 2 step: ${step}" >&2; failed=1 ;;
    esac
    [[ -n "${expected_script}" && "${script}" == "${expected_script}" ]] || {
      echo "ERROR: Stage 2 step script mismatch: ${step}" >&2
      failed=1
    }
    [[ -n "${expected_outputs}" && "${outputs}" == "${expected_outputs}" ]] || {
      echo "ERROR: Stage 2 step output mismatch: ${step}" >&2
      failed=1
    }
    expected_dependency="-"
    [[ "${step}" == "combinedpbmc" ]] && expected_dependency="gongsharma_cap"
    [[ "${dependency}" == "${expected_dependency}" ]] || {
      echo "ERROR: Stage 2 dependency mismatch: ${step}" >&2
      failed=1
    }
    expected_owner="-"
    [[ "${owner}" == "-" ]] || expected_owner="$(ecoda_owner_dir stage2 "${step}")"
    [[ "${owner}" == "${expected_owner}" ]] || {
      echo "ERROR: Stage 2 owner mismatch: ${step}" >&2
      failed=1
    }
    if [[ "${owner}" != "-" ]]; then
      owner_state="$(ecoda_owner_state "${owner}" 2>/dev/null || true)"
      [[ "${owner_state}" == "OK" ]] || {
        echo "ERROR: Stage 2 owner is not terminal OK: ${owner}" >&2
        failed=1
      }
      pending_count=$((pending_count + 1))
    fi
    old_ifs="${IFS}"
    IFS=';'
    read -r -a output_paths <<< "${outputs}"
    IFS="${old_ifs}"
    for path in "${output_paths[@]}"; do
      if ! ecoda_validate_checksum "${path}" ||
         ! ecoda_validate_stage2_output "${step}" "${path}"; then
        echo "ERROR: Stage 2 sync-only artifact failed semantic integrity: ${path}" >&2
        failed=1
      fi
    done
    seen_steps="${seen_steps} ${step}"
  done < "${MANIFEST}"

  if [[ -e "${JOB_FILE}" ]]; then
    ecoda_validate_run_owned_path "${JOB_FILE}" "${ECODA_RUN_ROOT}" ||
      stage2_abort "Stage 2 jobs manifest is not run-owned"
    if [[ ${pending_count} -gt 0 ]]; then
      ecoda_validate_manifest "${JOB_FILE}" 2 || {
        echo "ERROR: Stage 2 jobs manifest is invalid." >&2
        failed=1
      }
      job_count=0
      job_seen=""
      while IFS=$'\t' read -r job_step job_id; do
        job_count=$((job_count + 1))
        expected_owner="$(awk -F '\t' -v step="${job_step}" '$1 == step {print $5}' "${MANIFEST}")"
        [[ -n "${job_step}" && "${job_id}" =~ ^[0-9]+$ &&
           -n "${expected_owner}" && "${expected_owner}" != "-" ]] || failed=1
        case " ${job_seen} " in
          *" ${job_step} "*) failed=1 ;;
        esac
        job_seen="${job_seen} ${job_step}"
      done < "${JOB_FILE}"
      [[ ${job_count} -eq ${pending_count} ]] || failed=1
    elif [[ -s "${JOB_FILE}" ]]; then
      echo "ERROR: all-skipped Stage 2 run has nonempty jobs manifest." >&2
      failed=1
    fi
  elif [[ ${pending_count} -gt 0 ]]; then
    echo "ERROR: Stage 2 jobs manifest is missing for pending work." >&2
    failed=1
  fi
  if [[ ${failed} -ne 0 ]]; then
    stage2_abort "Stage 2 sync-only validation failed"
  fi
  ecoda_set_run_state OK "sync-only artifact and manifest validation passed" ||
    stage2_abort "failed to write Stage 2 terminal OK state"
  exit 0
fi

DATASET_NAMES=()
if [[ -n "${DATASETS_ARG}" ]]; then
  ecoda_split_csv "${DATASETS_ARG}"
  DATASET_NAMES=("${ECODA_ARRAY[@]}")
  ecoda_assert_unique_items "${DATASET_NAMES[@]}"
  for ds in "${DATASET_NAMES[@]}"; do
    ecoda_dataset_exists "${ds}" || { echo "ERROR: unknown dataset '${ds}'." >&2; exit 1; }
  done
else
  while IFS= read -r ds; do DATASET_NAMES+=("${ds}"); done < <(jq -r 'keys[] | select(startswith("_") | not)' "${DATASETS_JSON_FILE}")
fi

ALL_STEPS=(gongsharma_cap combinedpbmc joanito kfoury_lowres_ct myocardial_counts bassez_cellsubtype)
SELECTED_STEPS=()
if [[ -n "${STEPS_ARG}" ]]; then
  ecoda_split_csv "${STEPS_ARG}"
  SELECTED_STEPS=("${ECODA_ARRAY[@]}")
  ecoda_assert_unique_items "${SELECTED_STEPS[@]}"
else
  SELECTED_STEPS=("${ALL_STEPS[@]}")
fi
for step in "${SELECTED_STEPS[@]}"; do
  case " ${ALL_STEPS[*]} " in
    *" ${step} "*) ;;
    *) echo "ERROR: unknown Stage 2 step '${step}'." >&2; exit 1 ;;
  esac
done

# Dataset selection narrows fixed hooks without changing their scientific code.
# Generated datasets pull in their real prerequisites automatically.
if [[ -n "${DATASETS_ARG}" ]]; then
  FILTERED_STEPS=()
  for step in "${SELECTED_STEPS[@]}"; do
    keep=0
    for ds in "${DATASET_NAMES[@]}"; do
      case "${step}:${ds}" in
        gongsharma_cap:Gongsharma_cmv_young_males|gongsharma_cap:CombinedPBMC|combinedpbmc:CombinedPBMC|joanito:Joanito|joanito:_debug|kfoury_lowres_ct:Kfoury|myocardial_counts:Myocardial_infarction|bassez_cellsubtype:Bassez) keep=1 ;;
      esac
    done
    [[ ${keep} -eq 1 ]] && FILTERED_STEPS+=("${step}")
  done
  SELECTED_STEPS=("${FILTERED_STEPS[@]}")
  if [[ " ${DATASET_NAMES[*]} " == *" CombinedPBMC "* ]]; then
    case " ${SELECTED_STEPS[*]} " in *" gongsharma_cap "*) ;; *) SELECTED_STEPS=(gongsharma_cap "${SELECTED_STEPS[@]}") ;; esac
  fi
  if [[ ${#SELECTED_STEPS[@]} -eq 0 ]]; then
    echo "ERROR: selected datasets have no Stage 2 preprocessing hook." >&2
    exit 1
  fi
fi
if [[ " ${SELECTED_STEPS[*]} " == *" combinedpbmc "* ]]; then
  case " ${SELECTED_STEPS[*]} " in
    *" gongsharma_cap "*) ;;
    *) SELECTED_STEPS=(gongsharma_cap "${SELECTED_STEPS[@]}") ;;
  esac
fi
ORDERED_STEPS=()
for known_step in "${ALL_STEPS[@]}"; do
  for selected_step in "${SELECTED_STEPS[@]}"; do
    [[ "${selected_step}" == "${known_step}" ]] && ORDERED_STEPS+=("${known_step}")
  done
done
SELECTED_STEPS=("${ORDERED_STEPS[@]}")

RUN_ID="${ECODA_RUN_ID:-$(ecoda_new_run_id stage2)}"
ecoda_init_run stage2 "${RUN_ID}" >/dev/null
mkdir -p "${LOGS_DIR}" || stage2_abort "failed to create Stage 2 log directory"
MANIFEST="${ECODA_RUN_ROOT}/manifests/steps.tsv"
JOB_FILE="${ECODA_RUN_ROOT}/manifests/jobs.tsv"
MANIFEST_TMP="${MANIFEST}.build.$$"
JOB_FILE_TMP="${JOB_FILE}.build.$$"
: > "${MANIFEST_TMP}"
: > "${JOB_FILE_TMP}"


combinedpbmc_raw_migrate() {
  local data_dir="${HPC_SCRATCH_DIR}/CombinedPBMC/data"
  local old_path="${data_dir}/combined_pbmc_batch_effect_analysis.h5ad"
  local new_path="${data_dir}/combined_pbmc.h5ad"

  if [[ -f "${old_path}" && ! -e "${new_path}" ]]; then
    if ecoda_validate_checksum "${old_path}" &&
       ecoda_validate_stage2_output combinedpbmc "${old_path}"; then
      if ! cp "${old_path}" "${new_path}"; then
        return 1
      fi
      if ! ecoda_write_checksum "${new_path}" ||
         ! ecoda_validate_checksum_record "${new_path}" "${ECODA_CHECKSUM_MD5}" \
           "${ECODA_CHECKSUM_SIZE}" ||
         ! ecoda_validate_stage2_output combinedpbmc "${new_path}"; then
        rm -f "${new_path}" "${new_path}.md5"
        return 1
      fi
      if ! rm -f "${old_path}" "${old_path}.md5"; then
        rm -f "${new_path}" "${new_path}.md5"
        return 1
      fi
      echo "Migrated validated CombinedPBMC raw input to ${new_path}."
    else
      echo "Existing legacy CombinedPBMC raw input failed migration validation; regeneration required." >&2
    fi
  fi
  if [[ -f "${new_path}" && -f "${old_path}" ]] &&
     ecoda_validate_checksum "${new_path}" &&
     ecoda_validate_stage2_output combinedpbmc "${new_path}"; then
    rm -f "${old_path}" "${old_path}.md5"
    echo "Removed duplicate legacy CombinedPBMC raw input after canonical validation."
  fi
}
if [[ " ${SELECTED_STEPS[*]} " == *" combinedpbmc "* ]]; then
  combinedpbmc_raw_migrate || stage2_abort "CombinedPBMC raw input migration failed"
fi
MEMORY_CURRENT="${MEMORY}"
CAP_JOB_ID=""
PENDING_STEPS=()
OWNER_DIRS=()
ecoda_owner_clear_tracked
for step in "${SELECTED_STEPS[@]}"; do
  script="$(step_script "${step}")"
  [[ -f "${script}" ]] || stage2_abort "missing hook script: ${script}"
  outputs="$(step_outputs "${step}")"
  old_ifs="${IFS}"
  IFS=';'
  read -r -a output_paths <<< "${outputs}"
  IFS="${old_ifs}"
  valid=1
  for path in "${output_paths[@]}"; do
    if [[ ${FORCE_ARG} -eq 1 ]] || ! ecoda_validate_checksum "${path}" ||
       ! ecoda_validate_stage2_output "${step}" "${path}"; then
      valid=0
    fi
  done
  if [[ ${FORCE_ARG} -eq 0 && ${valid} -eq 1 ]]; then
    printf '%s\t%s\t%s\t%s\t-\n' "${step}" "${script}" "${outputs}" "-" >> "${MANIFEST_TMP}"
    echo "Skipping validated Stage 2 step ${step}."
    continue
  fi
  set +e
  owner_dir="$(ecoda_owner_acquire stage2 "${step}" "${RUN_ID}" "${FORCE_ARG}" 0)"
  owner_rc=$?
  set -e
  if [[ ${owner_rc} -ne 0 ]]; then
    [[ ${owner_rc} -eq 2 ]] && echo "ERROR: valid terminal owner exists for ${step}; selected artifact unexpectedly passed ownership gate." >&2
    stage2_abort "ownership conflict for ${step}"
  fi
  ecoda_owner_track "${owner_dir}" || stage2_abort "failed to track Stage 2 owner for ${step}"
  if [[ ${FORCE_ARG} -eq 1 ]]; then
    for path in "${output_paths[@]}"; do
      ecoda_invalidate_artifact "${path}" ||
        stage2_abort "failed to invalidate Stage 2 artifact: ${path}"
    done
  fi
  dependency="-"
  [[ "${step}" == "combinedpbmc" ]] && dependency="gongsharma_cap"
  printf '%s\t%s\t%s\t%s\t%s\n' "${step}" "${script}" "${outputs}" "${dependency}" "${owner_dir}" >> "${MANIFEST_TMP}"
  PENDING_STEPS+=("${step}")
  OWNER_DIRS+=("${owner_dir}")
done
if ! ecoda_atomic_install_manifest "${MANIFEST_TMP}" "${MANIFEST}" 5; then
  stage2_abort "failed to install Stage 2 manifest atomically"
fi
rm -f "${MANIFEST_TMP}"

if [[ ${#PENDING_STEPS[@]} -eq 0 ]]; then
  ecoda_atomic_write "${JOB_FILE}" "" || stage2_abort "failed to create empty Stage 2 jobs manifest"
  ecoda_atomic_write "${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv" "" ||
    stage2_abort "failed to create empty Stage 2 scheduler manifest"
  ecoda_set_run_state OK "all selected Stage 2 artifacts already validated" ||
    stage2_abort "failed to write Stage 2 terminal OK state"
  echo "STAGE2_RUN_ID=${RUN_ID}"
  exit 0
fi
export ECODA_RUNTIME_PROFILE=stage2
ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE}" ||
  stage2_abort "Stage 2 immutable runtime validation failed"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage2 0)" ||
  stage2_abort "Stage 2 runtime export construction failed"

# Submit every independent hook immediately. CombinedPBMC is the only explicit
# dependency edge; dependency submission does not serialize unrelated hooks.
for step in "${PENDING_STEPS[@]}"; do
  script="$(step_script "${step}")"
  SBATCH_ARGS=(--parsable --partition="${PARTITION}" --mem="${MEMORY_CURRENT}" \
    --output="${LOGS_DIR}/stage2_${step}_%j.log" \
    --error="${LOGS_DIR}/stage2_${step}_%j.err" --mail-user="${USER_EMAIL}" \
    --export="ALL,STAGE2_RUN_ROOT=${ECODA_RUN_ROOT},FORCE_PREPROCESS=${FORCE_ARG},${RUNTIME_EXPORT}")
  if [[ "${step}" == "combinedpbmc" && -n "${CAP_JOB_ID}" ]]; then
    SBATCH_ARGS+=(--dependency="afterok:${CAP_JOB_ID}")
  fi
  set +e
  job_output="$(sbatch "${SBATCH_ARGS[@]}" "${script}")"
  submit_rc=$?
  set -e
  if [[ ${submit_rc} -ne 0 ]]; then
    stage2_abort "sbatch rejected Stage 2 step ${step}"
  fi
  job_id="${job_output%%;*}"
  [[ "${job_id}" =~ ^[0-9]+$ ]] || stage2_abort "invalid sbatch id for ${step}: ${job_id}"
  printf '%s\t%s\n' "${step}" "${job_id}" >> "${JOB_FILE_TMP}"
  if ! ecoda_atomic_install_manifest "${JOB_FILE_TMP}" "${JOB_FILE}" 2; then
    stage2_abort "failed to install Stage 2 scheduler manifest after ${step}"
  fi
  [[ "${step}" == "gongsharma_cap" ]] && CAP_JOB_ID="${job_id}"
  echo "STAGE2_STEP_JOB_ID=${step}:${job_id}"
done
rm -f "${JOB_FILE_TMP}"
SCHEDULER_IDS_FILE="${ECODA_RUN_ROOT}/manifests/scheduler_ids.tsv"
SCHEDULER_IDS_TMP="${SCHEDULER_IDS_FILE}.build.$$"
: > "${SCHEDULER_IDS_TMP}"
while IFS=$'\t' read -r step scheduler_id; do
  [[ -n "${step}" && "${scheduler_id}" =~ ^[0-9]+$ ]] ||
    stage2_abort "invalid Stage 2 scheduler manifest row"
  printf 'ARRAY\t%s\n' "${scheduler_id}" >> "${SCHEDULER_IDS_TMP}"
done < "${JOB_FILE}"
if ! ecoda_atomic_install_manifest "${SCHEDULER_IDS_TMP}" "${SCHEDULER_IDS_FILE}" 2; then
  stage2_abort "failed to persist Stage 2 array scheduler IDs"
fi
rm -f "${SCHEDULER_IDS_TMP}"
JOB_IDS="$(cut -f2 "${JOB_FILE}" | paste -sd: -)"
set +e
# The watchdog deserializes the full Joanito RDS for mandatory semantic validation.
watchdog_output="$(sbatch --parsable --wait --dependency="afterany:${JOB_IDS}" \
  --partition="${PARTITION}" --ntasks=1 --cpus-per-task=1 --mem="${STAGE2_WATCHDOG_MEM:-64G}" \
  --time="${STAGE2_WATCHDOG_TIME_LIMIT:-12:00:00}" \
  --output="${LOGS_DIR}/stage2_watchdog_%j.log" \
  --error="${LOGS_DIR}/stage2_watchdog_%j.err" --mail-user="${USER_EMAIL}" \
  --export="ALL,STAGE2_RUN_ROOT=${ECODA_RUN_ROOT},STAGE2_FORCE=${FORCE_ARG},${RUNTIME_EXPORT}" \
  "${SCRIPT_DIR}/stage2_watchdog.sh" "${RUN_ID}" "${MANIFEST}" "${JOB_FILE}" \
  "${MEMORY}" "${MAX_MEMORY}" "${PARTITION}" "${THROTTLE}")"
watchdog_rc=$?
set -e
WATCHDOG_ID="${watchdog_output%%;*}"
[[ "${WATCHDOG_ID}" =~ ^[0-9]+$ ]] || stage2_abort "invalid Stage 2 watchdog id: ${WATCHDOG_ID}"
echo "STAGE2_RUN_ID=${RUN_ID}"
echo "STAGE2_WATCHDOG_JOB_ID=${WATCHDOG_ID}"
SCHEDULER_IDS_TMP="${SCHEDULER_IDS_FILE}.watchdog.build.$$"
cp "${SCHEDULER_IDS_FILE}" "${SCHEDULER_IDS_TMP}" ||
  stage2_abort "failed to copy Stage 2 scheduler manifest"
SCHEDULER_ID_SEEN=""
while IFS=$'\t' read -r scheduler_kind scheduler_id; do
  SCHEDULER_ID_SEEN="${SCHEDULER_ID_SEEN} ${scheduler_id}"
done < "${SCHEDULER_IDS_FILE}"
printf 'WATCHDOG\t%s\n' "${WATCHDOG_ID}" >> "${SCHEDULER_IDS_TMP}"
SCHEDULER_ID_SEEN="${SCHEDULER_ID_SEEN} ${WATCHDOG_ID}"
if [[ -s "${ECODA_RUN_ROOT}/status/watchdog" ]]; then
  while IFS= read -r status_line; do
    case "${status_line}" in
      SCHEDULER_ID=*)
        scheduler_id="${status_line#*=}"
        [[ "${scheduler_id}" =~ ^[0-9]+$ ]] ||
          stage2_abort "invalid Stage 2 watchdog scheduler ID"
        case " ${SCHEDULER_ID_SEEN} " in
          *" ${scheduler_id} "*) ;;
          *)
            printf 'STATUS\t%s\n' "${scheduler_id}" >> "${SCHEDULER_IDS_TMP}"
            SCHEDULER_ID_SEEN="${SCHEDULER_ID_SEEN} ${scheduler_id}"
            ;;
        esac
        ;;
    esac
  done < "${ECODA_RUN_ROOT}/status/watchdog"
fi
if ! ecoda_atomic_install_manifest "${SCHEDULER_IDS_TMP}" "${SCHEDULER_IDS_FILE}" 2; then
  stage2_abort "failed to install Stage 2 scheduler ID manifest atomically"
fi
rm -f "${SCHEDULER_IDS_TMP}"
if [[ ${watchdog_rc} -ne 0 ]]; then
  stage2_abort "Stage 2 watchdog submission or execution failed"
fi
if [[ "${STAGE2_SUBMITTER_TEST:-0}" == "1" ]]; then
  ecoda_set_run_state OK "submitter test mode; scheduler calls validated" ||
    stage2_abort "failed to write Stage 2 terminal OK state"
  exit 0
fi
if [[ ! -s "${ECODA_RUN_ROOT}/status/watchdog" ]] ||
   ! grep -q '^STATE=OK$' "${ECODA_RUN_ROOT}/status/watchdog"; then
  stage2_abort "Stage 2 watchdog did not report OK"
fi
while IFS= read -r status_line; do
  case "${status_line}" in
    SCHEDULER_ID=*)
      printf 'STAGE2_SCHEDULER_ID=%s:%s\n' "${RUN_ID}" "${status_line#*=}"
      ;;
  esac
done < "${ECODA_RUN_ROOT}/status/watchdog"
printf 'STAGE2_SCHEDULER_ID=%s:%s\n' "${RUN_ID}" "${WATCHDOG_ID}"
COMPLETE_SCHEDULER_TMP="${SCHEDULER_IDS_FILE}.complete.build.$$"
cp "${SCHEDULER_IDS_FILE}" "${COMPLETE_SCHEDULER_TMP}" ||
  stage2_abort "failed to copy complete Stage 2 scheduler manifest"
SCHEDULER_ID_SEEN=""
while IFS=$'\t' read -r scheduler_kind scheduler_id; do
  SCHEDULER_ID_SEEN="${SCHEDULER_ID_SEEN} ${scheduler_id}"
done < "${SCHEDULER_IDS_FILE}"
while IFS= read -r status_line; do
  case "${status_line}" in
    SCHEDULER_ID=*)
      scheduler_id="${status_line#*=}"
      case " ${SCHEDULER_ID_SEEN} " in
        *" ${scheduler_id} "*) ;;
        *)
          printf 'STATUS\t%s\n' "${scheduler_id}" >> "${COMPLETE_SCHEDULER_TMP}"
          SCHEDULER_ID_SEEN="${SCHEDULER_ID_SEEN} ${scheduler_id}"
          ;;
      esac
      ;;
  esac
done < "${ECODA_RUN_ROOT}/status/watchdog"
if ! ecoda_atomic_install_manifest "${COMPLETE_SCHEDULER_TMP}" "${SCHEDULER_IDS_FILE}" 2; then
  stage2_abort "failed to install complete Stage 2 scheduler ID manifest"
fi
rm -f "${COMPLETE_SCHEDULER_TMP}"
ecoda_set_run_state OK "Stage 2 watchdog completed and scheduler manifest is complete" ||
  stage2_abort "failed to write Stage 2 terminal OK state"
