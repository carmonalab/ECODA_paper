#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-annotation-stage.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

ORIGINAL_PATH="${PATH}"
BIN="${TMP_DIR}/bin"
HOME_ROOT="${TMP_DIR}/home"
REFERENCE_ROOT="${TMP_DIR}/reference"
HOST_ENV_PREFIX="${TMP_DIR}/host-env/.pixi/envs/py-cuda13"
SCGATE_INJECTION_MARKER="${TMP_DIR}/scgate-path-injected"
# Keep the immutable snapshot under a path that would execute a command if
# passed through a shell-reparsed sbatch --wrap string.
SNAPSHOT_PARENT="${TMP_DIR}/snapshots;touch ${SCGATE_INJECTION_MARKER};#"
RUNTIME_PARENT="${TMP_DIR}/runtime/_ecoda_runtime"
PRODUCER_RUN_ID="stage3_fixture"
CAPTURE="${TMP_DIR}/sbatch.calls"
R_CAPTURE="${TMP_DIR}/rscript.calls"
PY_CAPTURE="${TMP_DIR}/python.calls"
PY_AUX_CAPTURE="${TMP_DIR}/python.aux"
mkdir -p "${BIN}" "${HOME_ROOT}" "${REFERENCE_ROOT}" "${HOST_ENV_PREFIX}/bin" \
  "${SNAPSHOT_PARENT}" "${RUNTIME_PARENT}"
export CAPTURE R_CAPTURE PY_CAPTURE PY_AUX_CAPTURE

fail() {
  echo "annotation stage submitter: $*" >&2
  exit 1
}

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | cut -d' ' -f1
  else
    md5 -q "$1"
  fi
}

write_checksum() {
  local path="$1" digest size
  digest="$(md5_of "${path}")"
  size="$(wc -c < "${path}" | tr -d '[:space:]')"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "${size}" "${path}" > "${path}.md5"
}

# The reference-map staging helper has immutable, published MD5 values.  Keep
# the files in a temporary reference root and make only those four lookups
# return the published values; every other checksum remains real.
REAL_MD5SUM="$(command -v md5sum 2>/dev/null || true)"
REAL_MD5="$(command -v md5 2>/dev/null || true)"
export REFERENCE_ROOT REAL_MD5SUM REAL_MD5
cat > "${BIN}/md5sum" <<'STUB'
#!/bin/bash
set -euo pipefail
case "${1:-}" in
  "${REFERENCE_ROOT}/sketched_CD8T_human_ref_v1.rds")
    printf 'be86058ddafdd0154faf0485286b86e7  %s\n' "$1" ;;
  "${REFERENCE_ROOT}/sketched_CD4T_human_ref_v2.rds")
    printf '5540a0ee287e291528c96d476794b194  %s\n' "$1" ;;
  "${REFERENCE_ROOT}/sketched_DC_human_ref_v2.rds")
    printf '033d491ba7ca9bbf0badcae828e55b2c  %s\n' "$1" ;;
  "${REFERENCE_ROOT}/sketched_MoMac_human_v1.rds")
    printf '3043cd9058a8746d972c7be195b18e36  %s\n' "$1" ;;
  *)
    if [[ -n "${REAL_MD5SUM}" ]]; then
      exec "${REAL_MD5SUM}" "$@"
    fi
    [[ -n "${REAL_MD5}" ]] || exit 1
    printf '%s  %s\n' "$("${REAL_MD5}" -q "$1")" "$1"
    ;;
esac
STUB
chmod +x "${BIN}/md5sum"

for reference_name in \
  sketched_CD8T_human_ref_v1.rds \
  sketched_CD4T_human_ref_v2.rds \
  sketched_DC_human_ref_v2.rds \
  sketched_MoMac_human_v1.rds; do
  printf 'temporary reference fixture\n' > "${REFERENCE_ROOT}/${reference_name}"
done

# The scheduler stub is intentionally strict: a snapshot-backed Stage 4 run
# has six array/watchdog submissions and must never submit the frozen scGate
# writer via --wrap.  It also materializes deterministic worker outputs so the
# submitter's later chunk/merge boundaries expose their artifact records.
cat > "${BIN}/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
for arg in "$@"; do
  case "${arg}" in
    --wrap=*)
      echo "scGate writer was submitted through sbatch --wrap" >&2
      exit 97
      ;;
  esac
done
if [[ "$*" == *"1.2_prepare_chunks_worker.sh"* ]]; then
  "${SBATCH_HELPER}" prepare "${ECODA_RUN_ROOT}" >/dev/null
elif [[ "$*" == *"2.1_run_worker.sh"* ]]; then
  "${SBATCH_HELPER}" annotation "${ECODA_RUN_ROOT}" >/dev/null
elif [[ "$*" == *"3.2_merge_worker.sh"* ]]; then
  "${SBATCH_HELPER}" merge "${ECODA_RUN_ROOT}" >/dev/null
fi
if [[ "$*" == *"1.3_prepare_chunks_watchdog.sh"* ]]; then
  printf 'STATE=OK\nSCHEDULER_ID=%s\n' "72000$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" \
    > "${ECODA_RUN_ROOT}/status/preparation_watchdog"
elif [[ "$*" == *"1.2_annotation_watchdog.sh"* ]]; then
  printf 'STATE=OK\nSCHEDULER_ID=%s\n' "72000$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" \
    > "${ECODA_RUN_ROOT}/status/annotation_watchdog"
elif [[ "$*" == *"3.3_merge_watchdog.sh"* ]]; then
  printf 'STATE=OK\nSCHEDULER_ID=%s\n' "72000$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" \
    > "${ECODA_RUN_ROOT}/status/merge_watchdog"
fi
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '72000%s\n' "${N}"
STUB
chmod +x "${BIN}/sbatch"

cat > "${BIN}/rsync" <<'STUB'
#!/bin/bash
set -euo pipefail
FILES_FROM=""
SOURCE=""
DEST=""
for arg in "$@"; do
  case "${arg}" in
    --files-from=*) FILES_FROM="${arg#*=}" ;;
    -*) ;;
    *)
      if [[ -z "${SOURCE}" ]]; then SOURCE="${arg}"
      elif [[ -z "${DEST}" ]]; then DEST="${arg}"
      fi
      ;;
  esac
done
[[ -n "${FILES_FROM}" && -n "${SOURCE}" && -n "${DEST}" ]] || exit 2
while IFS= read -r relative; do
  [[ -n "${relative}" ]] || continue
  mkdir -p "${DEST%/}/$(dirname "${relative}")"
  cp "${SOURCE%/}/${relative}" "${DEST%/}/${relative}"
done < "${FILES_FROM}"
STUB
chmod +x "${BIN}/rsync"

cat > "${BIN}/apptainer" <<'STUB'
#!/bin/bash
set -euo pipefail
case "${1:-}" in
  inspect) exit 0 ;;
  *) exit 0 ;;
esac
STUB
chmod +x "${BIN}/apptainer"

cat > "${HOST_ENV_PREFIX}/bin/python" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ "${ECODA_AUX_ROOT:-}" == "${ECODA_SOURCE_ROOT}/aux" ]] || exit 71
map_path="${ECODA_AUX_ROOT:-}/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
printf '%s\n' "${map_path}" >> "${PY_AUX_CAPTURE}"
printf '%s\n' "$*" >> "${PY_CAPTURE}"
exit 0
STUB
chmod +x "${HOST_ENV_PREFIX}/bin/python"
cat > "${HOST_ENV_PREFIX}/bin/Rscript" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ "${ECODA_RUNTIME_MODE:-}" == "apptainer" ]] || exit 80
printf '%s\n' "$*" >> "${R_CAPTURE}"
case " $* " in
  *' --wrap '*|*' --force '*)
    echo "snapshot scGate validation received a writer flag" >&2
    exit 73
    ;;
esac
script=""
for arg in "$@"; do
  case "${arg}" in
    *.R) script="${arg}" ;;
  esac
done
[[ "${script}" == "${ECODA_SOURCE_ROOT}"/* ]] || exit 74
[[ "${script}" == */2.0_create_scgate_db.R ]] || exit 75
[[ "${SCGATE_DB_PATH:-}" == "${ECODA_AUX_ROOT}/scGateDB.rds" ]] || exit 76
[[ -r "${SCGATE_DB_PATH}" ]] || exit 77
case "$(cat "${SCGATE_DB_PATH}")" in
  *INVALID_FROZEN_SC_GATE*) exit 78 ;;
esac
case " $* " in
  *' --validate-only '*) exit 0 ;;
  *) exit 79 ;;
esac
STUB
chmod +x "${HOST_ENV_PREFIX}/bin/Rscript"

# This helper stands in for the actual worker process behind each captured
# scheduler command.  It publishes the same run-owned artifact records as the
# real preparation/annotation/merge workers, allowing the test to validate the
# records through the shared public validator rather than by inspecting source.
SBATCH_HELPER="${TMP_DIR}/materialize_stage4_outputs.sh"
export SBATCH_HELPER
cat > "${SBATCH_HELPER}" <<'HELPER'
#!/bin/bash
set -euo pipefail
kind="$1"
run_root="$2"
run_id="${run_root##*/}"
export ECODA_RUN_ROOT="${run_root}" ECODA_RUN_ID="${run_id}"
export ECODA_RUNS_ROOT="${run_root%/*}"
source "${ECODA_SOURCE_ROOT}/src/slurm_config.sh"
source "${ECODA_SOURCE_ROOT}/src/utils/bash/ecoda_run_common.sh"
case "${kind}" in
  prepare)
    while IFS=$'\t' read -r ds views owner_root; do
      union="${run_root}/datasets/${ds}/union/union.h5ad"
      chunk_dir="${run_root}/datasets/${ds}/chunks"
      mkdir -p "${chunk_dir}" "$(dirname "${union}")"
      printf 'synthetic union for %s\n' "${ds}" > "${union}"
      ecoda_write_checksum "${union}" >/dev/null
      ecoda_write_artifact_record "${union}" stage4_preparation "${run_id}" >/dev/null
      chunk="${chunk_dir}/chunk_1.txt"
      printf '%s\nfake_sample\n' "${union}" > "${chunk}"
    done < "${run_root}/manifests/preparation.tsv"
    ;;
  annotation)
    while IFS=$'\t' read -r ds chunk feather_dir; do
      chunk_num="${chunk##*/chunk_}"
      chunk_num="${chunk_num%.txt}"
      feather="${feather_dir}/annotations_chunk_${chunk_num}.feather"
      mkdir -p "${feather_dir}"
      printf 'synthetic annotation for %s\n' "${ds}" > "${feather}"
      ecoda_write_checksum "${feather}" >/dev/null
      ecoda_write_artifact_record "${feather}" stage4_annotation "${run_id}" >/dev/null
    done < "${run_root}/manifests/chunks.tsv"
    ;;
  merge)
    while IFS=$'\t' read -r ds views owner_root; do
      IFS=',' read -r -a view_list <<< "${views}"
      for view in "${view_list[@]}"; do
        output_name="$(jq -r --arg ds "${ds}" --arg view "${view}" \
          '.[$ds].views[$view].output_file_name // .[$ds].views[$view].output_file // empty' \
          "${DATASETS_JSON_FILE}")"
        path="${HPC_SCRATCH_DIR}/${ds}/output/${output_name}"
        mkdir -p "$(dirname "${path}")"
        [[ -s "${path}" ]] || printf 'synthetic merged output\n' > "${path}"
        [[ -s "${path}.md5" ]] || ecoda_write_checksum "${path}" >/dev/null
        ecoda_write_artifact_record "${path}" stage4_merge "${run_id}" >/dev/null
      done
    done < "${run_root}/manifests/merge.tsv"
    ;;
  *) exit 2 ;;
esac
HELPER
chmod +x "${SBATCH_HELPER}"

# Build a complete commit-keyed snapshot fixture.  The archive is made before
# the tree is made read-only and its digest is recorded in the manifest, so the
# runtime source validator checks the same archive/tree contract as production.
make_snapshot() {
  local snapshot_id="$1" db_mode="$2"
  local snapshot_root="${SNAPSHOT_PARENT}/${snapshot_id}"
  local tree="${snapshot_root}/tree"
  local identity="${snapshot_root}/identity"
  local archive_sha config_sha datasets_sha toml_sha lock_sha
  mkdir -p "${tree}" "${identity}"
  cp -R "${ROOT}/src" "${tree}/src"
  cp -R "${ROOT}/aux" "${tree}/aux"
  cp "${ROOT}/datasets.json" "${tree}/datasets.json"
  cp "${ROOT}/config_helper.R" "${tree}/config_helper.R"
  cp "${ROOT}/pixi.toml" "${tree}/pixi.toml"
  cp "${ROOT}/pixi.lock" "${tree}/pixi.lock"
  case "${db_mode}" in
    valid) ;;
    invalid) printf 'INVALID_FROZEN_SC_GATE\n' > "${tree}/aux/scGateDB.rds" ;;
    missing) rm -f "${tree}/aux/scGateDB.rds" ;;
    *) fail "unknown snapshot DB mode: ${db_mode}" ;;
  esac
  tar -cf "${identity}/source.tar" -C "${tree}" .
  archive_sha="$(sha256_of "${identity}/source.tar")"
  config_sha="$(sha256_of "${tree}/config_helper.R")"
  datasets_sha="$(sha256_of "${tree}/datasets.json")"
  toml_sha="$(sha256_of "${tree}/pixi.toml")"
  lock_sha="$(sha256_of "${tree}/pixi.lock")"
  cat > "${identity}/source.manifest" <<EOF
FORMAT=1
SOURCE_ROOT=${tree}
SOURCE_COMMIT=${snapshot_id}
SOURCE_ARCHIVE_PATH=${identity}/source.tar
SOURCE_ARCHIVE_SHA256=${archive_sha}
CONFIG_HELPER_SHA256=${config_sha}
DATASETS_SHA256=${datasets_sha}
PIXI_TOML_SHA256=${toml_sha}
PIXI_LOCK_SHA256=${lock_sha}
AUX_ROOT=${tree}/aux
SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4
EOF
  printf 'COMPLETE\n' > "${snapshot_root}/COMPLETE"
  chmod -R a-w "${snapshot_root}"
  printf '%s\n' "${snapshot_root}"
}

VALID_SNAPSHOT="$(make_snapshot aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa valid)"
INVALID_SNAPSHOT="$(make_snapshot bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb invalid)"
MISSING_SNAPSHOT="$(make_snapshot cccccccccccccccccccccccccccccccccccccccc missing)"
VALID_TREE="${VALID_SNAPSHOT}/tree"
VALID_SOURCE_MANIFEST="${VALID_SNAPSHOT}/identity/source.manifest"
INVALID_TREE="${INVALID_SNAPSHOT}/tree"
MISSING_TREE="${MISSING_SNAPSHOT}/tree"
SCGATE_EXPECTED_ROOT="${VALID_TREE}"
PY_EXPECTED_AUX_ROOT="${VALID_TREE}/aux"
export SCGATE_EXPECTED_ROOT PY_EXPECTED_AUX_ROOT

RUNTIME_IMAGE="${RUNTIME_PARENT}/runtime-test/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
mkdir -p "$(dirname "${RUNTIME_IMAGE}")"
printf 'deterministic runtime image fixture\n' > "${RUNTIME_IMAGE}"
IMAGE_SHA="$(sha256_of "${RUNTIME_IMAGE}")"
TOML_SHA="$(sha256_of "${VALID_TREE}/pixi.toml")"
LOCK_SHA="$(sha256_of "${VALID_TREE}/pixi.lock")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_BUILD_GIT_REVISION=build-revision-a
IMAGE_SHA256=${IMAGE_SHA}
IMAGE_PATH=${RUNTIME_IMAGE}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.27.1
APPTAINER_VERSION=1.3.2
IMAGE_PIXI_TOML_SHA256=${TOML_SHA}
IMAGE_PIXI_LOCK_SHA256=${LOCK_SHA}
EOF
chmod -R a-w "${RUNTIME_PARENT}"

SELECTION="${TMP_DIR}/selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\nKidney_KPMP_full\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\n' > "${SELECTION}"

make_base() {
  local base="$1"
  mkdir -p "${base}/home" "${base}/scratch/ECODA_paper" "${base}/logs" \
    "${base}/nas" "${base}/tmp"
}

prepare_inputs() {
  local base="$1" selection="$2" source_root="${3:-${VALID_TREE}}"
  local ds view output_name path unsuitable
  mkdir -p "${base}/scratch/ECODA_paper/_ecoda_runs/${PRODUCER_RUN_ID}/manifests/artifacts"
  while IFS=$'\t' read -r ds view; do
    unsuitable="$(jq -r --arg ds "${ds}" \
      '.[$ds].not_suitable_for_auto_annotation // [] |
       (index("hitme") != null and index("scatomic") != null)' "${ROOT}/datasets.json")"
    [[ "${unsuitable}" == "true" ]] && continue
    output_name="$(jq -r --arg ds "${ds}" --arg view "${view}" \
      '.[$ds].views[$view].output_file_name // .[$ds].views[$view].output_file // empty' \
      "${ROOT}/datasets.json")"
    [[ -n "${output_name}" ]] || fail "missing fixture output name for ${ds}/${view}"
    path="${base}/scratch/ECODA_paper/${ds}/output/${output_name}"
    mkdir -p "$(dirname "${path}")"
    printf 'synthetic Stage 4 input for %s/%s\n' "${ds}" "${view}" > "${path}"
    write_checksum "${path}"
    (
      export HOME="${base}/home" PATH="${BIN}:${ORIGINAL_PATH}" \
        HPC_SCRATCH_DIR="${base}/scratch/ECODA_paper" \
        ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
        ECODA_SOURCE_ROOT="${source_root}" ECODA_AUX_ROOT="${source_root}/aux" \
        ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_RUNTIME_IN_CONTAINER=0 \
        USER_EMAIL="test@example.invalid"
      source "${source_root}/src/slurm_config.sh"
      source "${source_root}/src/utils/bash/ecoda_run_common.sh"
      ecoda_artifact_owner_acquire "${path}" stage3 "${PRODUCER_RUN_ID}" 0 0 0 >/dev/null
      ecoda_write_artifact_record "${path}" stage3 "${PRODUCER_RUN_ID}" >/dev/null
      ecoda_artifact_owner_set_state "${path}" OK "fixture producer published" >/dev/null
    )
  done < "${selection}"
}

run_stage4() {
  local base="$1" snapshot="$2" test_mode="$3" script="$4"
  local runtime_image="$5" runtime_manifest="$6"
  shift 6
  local tree="${snapshot}/tree"
  local source_manifest="${snapshot}/identity/source.manifest"
  env -u ECODA_RUNTIME_MODE \
  HOME="${base}/home" \
  PATH="${BIN}:${ORIGINAL_PATH}" \
  HPC_SCRATCH_DIR="${base}/scratch/ECODA_paper" \
  ECODA_LOGS_DIR="${base}/logs" \
  TMPDIR="${base}/tmp" \
  NAS_TARGET_DIR="${base}/nas" \
  NAS_REF_DIR="${REFERENCE_ROOT}" \
  HOME_REF_DIR="${REFERENCE_ROOT}" \
  USER_EMAIL="test@example.invalid" \
  ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" \
  ECODA_SOURCE_ROOT="${tree}" \
  ECODA_SOURCE_MANIFEST="${source_manifest}" \
  ECODA_SOURCE_SNAPSHOT_REQUIRED=1 \
  ECODA_AUX_ROOT="${tree}/aux" \
  SCGATE_DB_BRANCH="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4" \
  ECODA_RUNTIME_IN_CONTAINER=0 \
  ECODA_RUNTIME_IMAGE="${runtime_image}" \
  ECODA_RUNTIME_MANIFEST="${runtime_manifest}" \
  ECODA_RUNTIME_PROFILE=stage4 \
  ECODA_RUNTIME_BUILD_VALIDATION=0 \
  ECODA_APPTAINER_NV=0 \
  APPTAINER_BIN="${BIN}/apptainer" \
  PYTHONDONTWRITEBYTECODE=1 \
  ANNOTATION_SUBMITTER_TEST="${test_mode}" \
  PREPROCESS_RUN_ID="" STAGE3_RUN_ID="${PRODUCER_RUN_ID}" INPUT_PRODUCER_RUN_ID="" \
  bash "${script}" "$@"
}

# A full snapshot-backed run exercises reference staging, frozen scGate
# validation, all preparation/annotation/merge scheduler boundaries, and the
# final selected sync without using an HPC scheduler.
SUCCESS_BASE="${TMP_DIR}/success"
make_base "${SUCCESS_BASE}"
prepare_inputs "${SUCCESS_BASE}" "${SELECTION}" "${VALID_TREE}"
: > "${CAPTURE}"
: > "${R_CAPTURE}"
: > "${PY_CAPTURE}"
: > "${PY_AUX_CAPTURE}"
SCGATE_SHA_BEFORE="$(sha256_of "${VALID_TREE}/aux/scGateDB.rds")"
SUCCESS_OUTPUT="$(run_stage4 "${SUCCESS_BASE}" "${VALID_SNAPSHOT}" 0 \
  "${VALID_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
  --selection-file "${SELECTION}" --exact-batch-selection --force)"
RUN_ID="$(printf '%s\n' "${SUCCESS_OUTPUT}" | sed -n 's/^ANNOTATION_RUN_ID=//p')"
[[ -n "${RUN_ID}" ]] || fail "successful snapshot run did not report a run ID"
RUN_ROOT="${SUCCESS_BASE}/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}"
[[ "$(cat "${RUN_ROOT}/status/skipped")" == $'Alzheimer\tSKIP_NOT_SUITABLE\nDiabetes\tSKIP_NOT_SUITABLE\nParkinson\tSKIP_NOT_SUITABLE' ]] ||
  fail "exact batch exemption set changed"
[[ "$(wc -l < "${RUN_ROOT}/manifests/runnable_selection.tsv" | tr -d '[:space:]')" == "9" ]] ||
  fail "exact batch selection did not preserve nine runnable rows"
SCHEDULER_MANIFEST="${RUN_ROOT}/manifests/scheduler_ids.tsv"
SCHEDULER_UNIQUE_IDS="$(
  awk -F '\t' '
    NF != 2 || $1 !~ /^(ARRAY|WATCHDOG|STATUS)$/ || $2 !~ /^[0-9]+$/ {
      invalid = 1
      next
    }
    {
      seen[$2] = 1
      if ($1 == "ARRAY" || $1 == "WATCHDOG") {
        job_rows++
        if (job_seen[$2]++) duplicate_job_id = 1
      }
    }
    END {
      if (invalid || job_rows != 6 || duplicate_job_id) exit 2
      count = 0
      for (scheduler_id in seen) count++
      print count
    }
  ' "${SCHEDULER_MANIFEST}"
)" || fail "scheduler ID manifest had malformed or duplicate job rows"
[[ "${SCHEDULER_UNIQUE_IDS}" == "6" ]] ||
  fail "scheduler ID manifest did not contain six unique scheduler IDs"
[[ "$(sed -n '4p' "${RUN_ROOT}/manifests/selection.tsv")" == $'Kidney_KPMP_full\tbatch_effect_uncorrected' ]] ||
  fail "run selection manifest did not retain Kidney_KPMP_full"
[[ -s "${RUN_ROOT}/manifests/source.manifest" ]] || fail "run source manifest was not copied"
[[ -s "${RUN_ROOT}/manifests/runtime.identity" ]] || fail "run runtime identity was not written"
[[ "$(wc -l < "${RUN_ROOT}/manifests/runtime.identity" | tr -d '[:space:]')" == "8" ]] ||
  fail "FORMAT=2 runtime identity did not carry all run-bound fields"
cmp -s "${RUN_ROOT}/manifests/source.manifest" "${VALID_SOURCE_MANIFEST}" ||
  fail "run source manifest differs from immutable snapshot manifest"
case "$(cat "${RUN_ROOT}/metadata")" in
  *"SOURCE_ROOT=${VALID_TREE}"*"RUNTIME_IMAGE=${RUNTIME_IMAGE}"*) ;;
  *) fail "run metadata omitted source/runtime identity" ;;
esac

[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "6" ]] ||
  fail "Stage 4 snapshot run did not submit exactly six array/watchdog jobs"
[[ "$(wc -l < "${R_CAPTURE}" | tr -d '[:space:]')" == "1" ]] ||
  fail "snapshot mode invoked scGate more than once"
R_CALL="$(cat "${R_CAPTURE}")"
case " ${R_CALL} " in *' --validate-only '*) ;; *) fail "snapshot scGate validation did not use --validate-only" ;; esac
case " ${R_CALL} " in *' --wrap '*|*' --force '*) fail "snapshot scGate validation received a writer flag" ;; esac
[[ ! -e "${SCGATE_INJECTION_MARKER}" ]] ||
  fail "hostile snapshot path triggered shell command injection"
[[ "$(sha256_of "${VALID_TREE}/aux/scGateDB.rds")" == "${SCGATE_SHA_BEFORE}" ]] ||
  fail "snapshot scGate validation changed frozen aux/scGateDB.rds"

PREP_CALL=""
ANNOT_CALL=""
MERGE_CALL=""
while IFS= read -r call; do
  for token in \
    "ECODA_SOURCE_ROOT=${VALID_TREE}" \
    "ECODA_SOURCE_MANIFEST=${VALID_SOURCE_MANIFEST}" \
    "ECODA_SOURCE_SNAPSHOT_REQUIRED=1" \
    "ECODA_RUN_ID=${RUN_ID}" \
    "ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}" \
    "ECODA_RUNTIME_MODE=apptainer" \
    "ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}" \
    "ECODA_RUNTIME_IMAGE_SHA256=" \
    "ECODA_RUNTIME_MANIFEST_SHA256=" \
    "ECODA_RUNTIME_IMAGE_SIZE=" \
    "ECODA_RUNTIME_MANIFEST_SIZE="; do
    case "${call}" in *"${token}"*) ;; *) fail "scheduler command omitted ${token}" ;; esac
  done
  case "${call}" in
    *"${ROOT}/src/4_cell_type_annotation/"*) fail "scheduler command used the mutable checkout" ;;
  esac
  case "${call}" in
    *"${VALID_TREE}/src/4_cell_type_annotation/"*) ;;
    *) fail "scheduler command was not rooted in the immutable source tree" ;;
  esac
  case "${call}" in
    *"1.2_prepare_chunks_worker.sh"*) PREP_CALL="${call}" ;;
    *"1.3_prepare_chunks_watchdog.sh"*) ;;
    *"2.1_run_worker.sh"*) ANNOT_CALL="${call}" ;;
    *"1.2_annotation_watchdog.sh"*) ;;
    *"3.2_merge_worker.sh"*) MERGE_CALL="${call}" ;;
    *"3.3_merge_watchdog.sh"*) ;;
    *) fail "captured an unknown Stage 4 scheduler script" ;;
  esac
done < "${CAPTURE}"
[[ -n "${PREP_CALL}" && -n "${ANNOT_CALL}" && -n "${MERGE_CALL}" ]] ||
  fail "preparation, annotation, and merge scheduler commands were not all captured"
case "${PREP_CALL}" in *"ANNOTATION_PREP_MANIFEST=${RUN_ROOT}/manifests/preparation.tsv"*) ;; *) fail "preparation command omitted its run-owned manifest" ;; esac
case "${ANNOT_CALL}" in *"CHUNKS_MANIFEST=${RUN_ROOT}/manifests/chunks.tsv"*) ;; *) fail "annotation command omitted its run-owned chunk manifest" ;; esac
case "${MERGE_CALL}" in *"ANNOTATION_MERGE_MANIFEST=${RUN_ROOT}/manifests/merge.tsv"*) ;; *) fail "merge command omitted its run-owned manifest" ;; esac

# The materialized worker fixtures publish real records.  Validate both
# producers through the shared record API, then prove a missing record and a
# producer mismatch are rejected rather than silently accepted.
(
  export HOME="${SUCCESS_BASE}/home" PATH="${BIN}:${ORIGINAL_PATH}" \
    HPC_SCRATCH_DIR="${SUCCESS_BASE}/scratch/ECODA_paper" \
    ECODA_LOGS_DIR="${SUCCESS_BASE}/logs" TMPDIR="${SUCCESS_BASE}/tmp" \
    ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" ECODA_SOURCE_ROOT="${VALID_TREE}" \
    ECODA_SOURCE_MANIFEST="${VALID_SOURCE_MANIFEST}" ECODA_SOURCE_SNAPSHOT_REQUIRED=1 \
    ECODA_AUX_ROOT="${VALID_TREE}/aux" ECODA_RUNTIME_MODE=apptainer \
    ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
    APPTAINER_BIN="${BIN}/apptainer" ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID="${RUN_ID}"
  source "${VALID_TREE}/src/slurm_config.sh"
  source "${VALID_TREE}/src/utils/bash/ecoda_run_common.sh"
  first_feather=""
  while IFS=$'\t' read -r ds chunk feather_dir; do
    chunk_num="${chunk##*/chunk_}"
    chunk_num="${chunk_num%.txt}"
    first_feather="${feather_dir}/annotations_chunk_${chunk_num}.feather"
    break
  done < "${RUN_ROOT}/manifests/chunks.tsv"
  [[ -s "${first_feather}" ]] || exit 1
  annotation_record="$(ecoda_artifact_record_path "${first_feather}" "${RUN_ID}")"
  [[ -s "${annotation_record}" ]] || exit 1
  ecoda_validate_artifact_record "${first_feather}" stage4_annotation "${RUN_ID}" >/dev/null
  merge_name="$(jq -r '.Breast_cancer.views.batch_effect_uncorrected.output_file_name // .Breast_cancer.views.batch_effect_uncorrected.output_file' "${VALID_TREE}/datasets.json")"
  merge_path="${SUCCESS_BASE}/scratch/ECODA_paper/Breast_cancer/output/${merge_name}"
  merge_record="$(ecoda_artifact_record_path "${merge_path}" "${RUN_ID}")"
  [[ -s "${merge_record}" ]] || exit 1
  ecoda_validate_artifact_record "${merge_path}" stage4_merge "${RUN_ID}" >/dev/null
  if ecoda_validate_artifact_record "${first_feather}" stage4_merge "${RUN_ID}" >/dev/null 2>&1; then
    exit 1
  fi
  rm -f "${annotation_record}"
  if ecoda_validate_artifact_record "${first_feather}" stage4_annotation "${RUN_ID}" >/dev/null 2>&1; then
    exit 1
  fi
  ecoda_write_artifact_record "${first_feather}" stage4_annotation "${RUN_ID}" >/dev/null
  ecoda_validate_artifact_record "${first_feather}" stage4_annotation "${RUN_ID}" >/dev/null
)
[[ -s "${PY_AUX_CAPTURE}" ]] || fail "preprocessing fixture did not consult the Ensembl map"
while IFS= read -r map_path; do
  [[ "${map_path}" == "${VALID_TREE}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz" ]] ||
    fail "preprocessing fixture selected an Ensembl map outside the source snapshot"
done < "${PY_AUX_CAPTURE}"

# A script rooted in the mutable checkout is rejected before the first sbatch.
OUTSIDE_BASE="${TMP_DIR}/outside-root"
make_base "${OUTSIDE_BASE}"
: > "${CAPTURE}"
set +e
run_stage4 "${OUTSIDE_BASE}" "${VALID_SNAPSHOT}" 1 \
  "${ROOT}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
  --datasets _debug --views benchmark_analysis >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]] || fail "mutable-checkout scheduler script was accepted"
[[ ! -s "${CAPTURE}" ]] || fail "script-outside-root failure happened after sbatch"

# Reserve one exact output path for a foreign active run.  Ownership must fail
# before any scheduler command, even though the source/runtime identities are
# otherwise valid.
OWNER_BASE="${TMP_DIR}/ownership-conflict"
make_base "${OWNER_BASE}"
OWNER_OUTPUT_NAME="$(jq -r '._debug.views.benchmark_analysis.output_file_name // ._debug.views.benchmark_analysis.output_file' "${VALID_TREE}/datasets.json")"
OWNER_OUTPUT="${OWNER_BASE}/scratch/ECODA_paper/_debug/output/${OWNER_OUTPUT_NAME}"
mkdir -p "$(dirname "${OWNER_OUTPUT}")"
(
  export HOME="${OWNER_BASE}/home" PATH="${BIN}:${ORIGINAL_PATH}" \
    HPC_SCRATCH_DIR="${OWNER_BASE}/scratch/ECODA_paper" \
    ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}" ECODA_AUX_ROOT="${VALID_TREE}/aux"
  source "${VALID_TREE}/src/slurm_config.sh"
  source "${VALID_TREE}/src/utils/bash/ecoda_run_common.sh"
  ecoda_artifact_owner_acquire "${OWNER_OUTPUT}" stage4 foreign_run 0 0 >/dev/null
)
: > "${CAPTURE}"
set +e
run_stage4 "${OWNER_BASE}" "${VALID_SNAPSHOT}" 1 \
  "${VALID_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
  --datasets _debug --views benchmark_analysis >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]] || fail "foreign active output owner was accepted"
[[ ! -s "${CAPTURE}" ]] || fail "ownership failure happened after sbatch"

# Snapshot compute must fail closed when its upstream producer identity is
# absent, before any scheduler boundary.
NO_PRODUCER_BASE="${TMP_DIR}/missing-producer"
make_base "${NO_PRODUCER_BASE}"
printf '_debug\tbenchmark_analysis\n' > "${TMP_DIR}/debug-selection.tsv"
prepare_inputs "${NO_PRODUCER_BASE}" "${TMP_DIR}/debug-selection.tsv" "${VALID_TREE}"
: > "${CAPTURE}"
SAVED_PRODUCER_RUN_ID="${PRODUCER_RUN_ID}"
PRODUCER_RUN_ID=""
set +e
run_stage4 "${NO_PRODUCER_BASE}" "${VALID_SNAPSHOT}" 0 \
  "${VALID_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
  --selection-file "${TMP_DIR}/debug-selection.tsv" >/dev/null 2>&1
RC=$?
set -e
PRODUCER_RUN_ID="${SAVED_PRODUCER_RUN_ID}"
[[ ${RC} -ne 0 ]] || fail "snapshot compute without producer identity was accepted"
[[ ! -s "${CAPTURE}" ]] || fail "missing producer identity reached sbatch"

# Snapshot mode validates a frozen DB and never falls back to a writer.  An
# invalid DB is rejected by --validate-only, and a missing DB is rejected by
# snapshot input validation, both before scheduler submission.
INVALID_BASE="${TMP_DIR}/invalid-scgate"
make_base "${INVALID_BASE}"
prepare_inputs "${INVALID_BASE}" "${TMP_DIR}/debug-selection.tsv" "${INVALID_TREE}"
: > "${CAPTURE}"
: > "${R_CAPTURE}"
INVALID_DB_SHA="$(sha256_of "${INVALID_TREE}/aux/scGateDB.rds")"
set +e
run_stage4 "${INVALID_BASE}" "${INVALID_SNAPSHOT}" 0 \
  "${INVALID_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
  --selection-file "${TMP_DIR}/debug-selection.tsv" --force >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]] || fail "invalid frozen scGate DB was accepted"
[[ ! -s "${CAPTURE}" ]] || fail "invalid frozen scGate DB reached sbatch"
[[ "$(sha256_of "${INVALID_TREE}/aux/scGateDB.rds")" == "${INVALID_DB_SHA}" ]] ||
  fail "invalid frozen scGate DB was modified"
case " $(cat "${R_CAPTURE}" 2>/dev/null || true) " in *' --validate-only '*) ;; *) fail "invalid frozen DB did not use --validate-only" ;; esac

MISSING_BASE="${TMP_DIR}/missing-scgate"
make_base "${MISSING_BASE}"
: > "${CAPTURE}"
: > "${R_CAPTURE}"
set +e
run_stage4 "${MISSING_BASE}" "${MISSING_SNAPSHOT}" 1 \
  "${MISSING_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
  --datasets _debug --views benchmark_analysis >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]] || fail "missing frozen scGate DB was accepted"
[[ ! -s "${CAPTURE}" ]] || fail "missing frozen scGate DB reached sbatch"
[[ ! -s "${R_CAPTURE}" ]] || fail "missing frozen scGate DB reached R validation"

# Preserve the immutable-runtime preflight contract: an Apptainer submission
# without its versioned image/manifest fails before the scheduler.
MISSING_RUNTIME_BASE="${TMP_DIR}/missing-runtime"
make_base "${MISSING_RUNTIME_BASE}"
: > "${CAPTURE}"
set +e
run_stage4 "${MISSING_RUNTIME_BASE}" "${VALID_SNAPSHOT}" 1 \
  "${VALID_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${MISSING_RUNTIME_BASE}/missing.sif" "${MISSING_RUNTIME_BASE}/missing.sif.manifest" \
  --datasets _debug --views benchmark_analysis >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]] || fail "missing immutable runtime image was accepted"
[[ ! -s "${CAPTURE}" ]] || fail "missing runtime failure happened after sbatch"

# Exact batch selection remains order-sensitive and uses the current full
# Kidney key.  The old key, old view, and missing-row fixtures all fail closed.
for bad_kind in old_key legacy missing; do
  BAD="${TMP_DIR}/bad-${bad_kind}.tsv"
  case "${bad_kind}" in
    old_key) sed '4s/Kidney_KPMP_full/Kidney_KPMP/' "${SELECTION}" > "${BAD}" ;;
    legacy) sed '1s/batch_effect_uncorrected/batch_effect_analysis/' "${SELECTION}" > "${BAD}" ;;
    missing) sed '12d' "${SELECTION}" > "${BAD}" ;;
  esac
  BAD_BASE="${TMP_DIR}/bad-${bad_kind}"
  make_base "${BAD_BASE}"
  : > "${CAPTURE}"
  set +e
  run_stage4 "${BAD_BASE}" "${VALID_SNAPSHOT}" 1 \
    "${VALID_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
    "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
    --selection-file "${BAD}" --exact-batch-selection >/dev/null 2>&1
  RC=$?
  set -e
  [[ ${RC} -ne 0 ]] || fail "bad exact selection was accepted: ${bad_kind}"
  [[ ! -s "${CAPTURE}" ]] || fail "bad exact selection reached sbatch: ${bad_kind}"
done

# Keep the dataset/view manifest mode assertion independent of exact batch
# selection, as a compatibility check for callers that do not use the twelve
# row fixture.
DATASET_BASE="${TMP_DIR}/dataset-mode"
make_base "${DATASET_BASE}"
: > "${CAPTURE}"
DATASET_OUTPUT="$(run_stage4 "${DATASET_BASE}" "${VALID_SNAPSHOT}" 1 \
  "${VALID_TREE}/src/4_cell_type_annotation/1_submit_onboarding_stage.sh" \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" \
  --datasets _debug --views benchmark_analysis,batch_effect_uncorrected)"
DATASET_MANIFEST="$(printf '%s\n' "${DATASET_OUTPUT}" | sed -n 's/^ANNOTATION_SELECTION_MANIFEST=//p')"
[[ "$(wc -l < "${DATASET_MANIFEST}" | tr -d '[:space:]')" == "2" ]] || fail "dataset mode did not preserve both selected views"
[[ "$(sed -n '1p' "${DATASET_MANIFEST}")" == $'_debug\tbenchmark_analysis' ]] || fail "dataset mode changed first view"
[[ "$(sed -n '2p' "${DATASET_MANIFEST}")" == $'_debug\tbatch_effect_uncorrected' ]] || fail "dataset mode changed second view"

echo "annotation stage submitter: OK"
