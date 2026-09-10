#!/bin/bash
# Focused contract test for the canonical manifest-driven preprocessing gate.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-preprocess-stage.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${TMP_DIR}/nas" "${TMP_DIR}/logs"

CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}
md5_file() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | cut -d' ' -f1
  else
    md5 -q "$1"
  fi
}
write_sidecar() {
  local path="$1" digest
  digest="$(md5_file "${path}")"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' \
    "${digest}" "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
}

# Build a small but complete immutable source snapshot.  The source archive,
# manifest, and read-only tree are all temporary fixtures owned by this test.
SNAPSHOT_COMMIT="0000000000000000000000000000000000000001"
SNAPSHOT_ROOT="${TMP_DIR}/snapshots/${SNAPSHOT_COMMIT}"
SOURCE_ROOT="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY}/source.tar"
mkdir -p \
  "${SOURCE_ROOT}/src/3_scrnaseq_preprocessing" \
  "${SOURCE_ROOT}/src/utils/bash" \
  "${SOURCE_ROOT}/src/utils/py" \
  "${SOURCE_ROOT}/aux" \
  "${SOURCE_IDENTITY}"
cp "${ROOT}/src/slurm_config.sh" "${SOURCE_ROOT}/src/slurm_config.sh"
for source_file in \
  ecoda_run_common.sh ecoda_runtime.sh h5ad_preflight_submit.sh \
  h5ad_preflight_worker.sh worker_retry.sh; do
  cp "${ROOT}/src/utils/bash/${source_file}" \
    "${SOURCE_ROOT}/src/utils/bash/${source_file}"
done
for source_file in 1.1_run_worker.sh 1.2_preprocess_watchdog.sh; do
  cp "${ROOT}/src/3_scrnaseq_preprocessing/${source_file}" \
    "${SOURCE_ROOT}/src/3_scrnaseq_preprocessing/${source_file}"
done
cp "${ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
  "${SOURCE_ROOT}/src/utils/py/benchmark_h5ad_contract.py"
for source_file in config_helper.R datasets.json pixi.toml pixi.lock; do
  cp "${ROOT}/${source_file}" "${SOURCE_ROOT}/${source_file}"
done
for source_file in scGateDB.rds genes.blocklist.rds \
  EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
  cp "${ROOT}/aux/${source_file}" "${SOURCE_ROOT}/aux/${source_file}"
done
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_ROOT}" .
SOURCE_ARCHIVE_SHA256="$(sha256_file "${SOURCE_ARCHIVE}")"
SOURCE_CONFIG_SHA256="$(sha256_file "${SOURCE_ROOT}/config_helper.R")"
SOURCE_DATASETS_SHA256="$(sha256_file "${SOURCE_ROOT}/datasets.json")"
SOURCE_TOML_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.toml")"
SOURCE_LOCK_SHA256="$(sha256_file "${SOURCE_ROOT}/pixi.lock")"
chmod -R a-w "${SOURCE_ROOT}"
printf '%s\n' \
  'FORMAT=1' \
  "SOURCE_ROOT=${SOURCE_ROOT}" \
  "SOURCE_COMMIT=${SNAPSHOT_COMMIT}" \
  "SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}" \
  "SOURCE_ARCHIVE_SHA256=${SOURCE_ARCHIVE_SHA256}" \
  "CONFIG_HELPER_SHA256=${SOURCE_CONFIG_SHA256}" \
  "DATASETS_SHA256=${SOURCE_DATASETS_SHA256}" \
  "PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}" \
  "PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}" \
  "AUX_ROOT=${SOURCE_ROOT}/aux" \
  'SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4' > "${SOURCE_MANIFEST}"
touch "${SNAPSHOT_ROOT}/COMPLETE"
chmod a-w "${SOURCE_MANIFEST}" "${SOURCE_ARCHIVE}" "${SNAPSHOT_ROOT}/COMPLETE"

# The runtime fixture is FORMAT=2 and deliberately lives below a versioned
# _ecoda_runtime/<id>/ directory.  The fake Apptainer only services inspect;
# no container or scheduler is launched by this test.
HOST_ENV_PREFIX="${TMP_DIR}/host/.pixi/envs/py-cuda13"
mkdir -p "${HOST_ENV_PREFIX}/bin" "${HOST_ENV_PREFIX}/lib" \
  "${TMP_DIR}/home/scratch/ECODA_paper" "${TMP_DIR}/logs"
cat > "${HOST_ENV_PREFIX}/bin/python" <<'STUB'
#!/bin/bash
exit 0
STUB
cat > "${HOST_ENV_PREFIX}/bin/Rscript" <<'STUB'
#!/bin/bash
exit 0
STUB
chmod +x "${HOST_ENV_PREFIX}/bin/python" "${HOST_ENV_PREFIX}/bin/Rscript"
HOST_PYTHON_SHA256="$(sha256_file "${HOST_ENV_PREFIX}/bin/python")"
HOST_RSCRIPT_SHA256="$(sha256_file "${HOST_ENV_PREFIX}/bin/Rscript")"

APPTAINER_STUB="${TMP_DIR}/bin/apptainer"
cat > "${APPTAINER_STUB}" <<'STUB'
#!/bin/bash
set -euo pipefail
if [[ "${1:-}" == "inspect" ]]; then
  exit 0
fi
if [[ "${1:-}" == "exec" ]]; then
  shift
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --env|--bind|--no-mount) shift 2 ;;
      --containall|--no-home|--nv) shift ;;
      *)
        image="$1"
        shift
        command="$1"
        shift
        script="$1"
        shift
        exec "${command}" "${script}" "$@"
        ;;
    esac
  done
fi
exit 1
STUB
chmod +x "${APPTAINER_STUB}"

RUNTIME_ID="stage3-test-runtime"
RUNTIME_DIR="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runtime/${RUNTIME_ID}"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
mkdir -p "${RUNTIME_DIR}"
printf 'format2-runtime-fixture\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA256="$(sha256_file "${RUNTIME_IMAGE}")"
printf '%s\n' \
  'FORMAT=2' \
  "IMAGE_PATH=${RUNTIME_IMAGE}" \
  "IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}" \
  'RUNTIME_ENV=py-cuda13' \
  'RUNTIME_LAYOUT=relocated' \
  'CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13' \
  'BASE_IMAGE=rockylinux:9' \
  'PIXITAINER_VERSION=0.8.3' \
  'PIXI_VERSION=0.49.0' \
  'APPTAINER_VERSION=1.3.2' \
  'IMAGE_BUILD_GIT_REVISION=immutable-runtime-build' \
  "IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}" \
  "IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}" > "${RUNTIME_MANIFEST}"
chmod 444 "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}"
chmod 555 "${RUNTIME_DIR}"

export PATH="${TMP_DIR}/bin:${PATH}"
export HOME="${TMP_DIR}/home"
export HPC_SCRATCH_DIR="${TMP_DIR}/home/scratch/ECODA_paper"
export NAS_TARGET_DIR="${TMP_DIR}/nas"
export ECODA_LOGS_DIR="${TMP_DIR}/logs"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_HOST_ENV_PREFIX="${HOST_ENV_PREFIX}"
export ECODA_HOST_PYTHON_SHA256="${HOST_PYTHON_SHA256}"
export ECODA_HOST_RSCRIPT_SHA256="${HOST_RSCRIPT_SHA256}"
export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export APPTAINER_BIN="${APPTAINER_STUB}"
export USER_EMAIL="test@example.invalid"

# The scheduler stub records every boundary.  NOOP_PREFLIGHT mode executes the
# immutable preflight worker locally so a valid existing output can complete
# the validator-only path without a real scheduler.
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
CALL_COUNT="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
if [[ "${PREPROCESS_NOOP_PREFLIGHT:-0}" == "1" ]]; then
  export_arg=""
  worker_script=""
  for arg in "$@"; do
    case "${arg}" in
      --export=*) export_arg="${arg#--export=}" ;;
      *.sh) worker_script="${arg}" ;;
    esac
  done
  export_field() {
    local key="$1" item
    local fields=()
    IFS=',' read -r -a fields <<< "${export_arg#ALL,}"
    for item in "${fields[@]}"; do
      case "${item}" in
        "${key}"=*) printf '%s' "${item#*=}"; return 0 ;;
      esac
    done
    return 1
  }
  preflight_manifest="$(export_field H5AD_PREFLIGHT_MANIFEST)"
  preflight_status_dir="$(export_field H5AD_PREFLIGHT_STATUS_DIR)"
  preflight_run_root="$(export_field H5AD_PREFLIGHT_RUN_ROOT)"
  preflight_run_id="$(export_field H5AD_PREFLIGHT_RUN_ID)"
  source_root="$(export_field ECODA_SOURCE_ROOT)"
  source_manifest="$(export_field ECODA_SOURCE_MANIFEST)"
  host_prefix="$(export_field ECODA_HOST_ENV_PREFIX)"
  host_python_sha="$(export_field ECODA_HOST_PYTHON_SHA256)"
  host_rscript_sha="$(export_field ECODA_HOST_RSCRIPT_SHA256)"
  runtime_image="$(export_field ECODA_RUNTIME_IMAGE)"
  runtime_manifest="$(export_field ECODA_RUNTIME_MANIFEST)"
  scratch_root="$(export_field HPC_SCRATCH_DIR)"
  logs_root="$(export_field ECODA_LOGS_DIR)"
  ECODA_RUNTIME_MODE=host ECODA_RUNTIME_IN_CONTAINER=0 \
  ECODA_RUNTIME_PROFILE=stage3 ECODA_SOURCE_ROOT="${source_root}" \
  ECODA_SOURCE_MANIFEST="${source_manifest}" ECODA_SOURCE_SNAPSHOT_REQUIRED=1 \
  ECODA_HOST_ENV_PREFIX="${host_prefix}" ECODA_RUNTIME_IMAGE="${runtime_image}" \
  ECODA_HOST_PYTHON_SHA256="${host_python_sha}" \
  ECODA_HOST_RSCRIPT_SHA256="${host_rscript_sha}" \
  ECODA_RUNTIME_MANIFEST="${runtime_manifest}" HPC_SCRATCH_DIR="${scratch_root}" \
  ECODA_SCRATCH_ROOT="${scratch_root}" ECODA_LOGS_DIR="${logs_root}" \
  LOGS_DIR="${logs_root}" ECODA_AUX_ROOT="${source_root}/aux" \
  H5AD_PREFLIGHT_MANIFEST="${preflight_manifest}" \
  H5AD_PREFLIGHT_STATUS_DIR="${preflight_status_dir}" \
  H5AD_PREFLIGHT_RUN_ROOT="${preflight_run_root}" \
  H5AD_PREFLIGHT_RUN_ID="${preflight_run_id}" H5AD_PREFLIGHT_MODE=classify \
  H5AD_PREFLIGHT_PYTHON_BIN="${host_prefix}/bin/python" \
  H5AD_PREFLIGHT_TASK_ID=1 bash "${worker_script}"
  printf '600003\n'
  exit 0
fi
case "${CALL_COUNT}" in
  1) printf '600001\n' ;;
  2) printf '600002\n' ;;
  *) printf 'unexpected sbatch call count\n' >&2; exit 1 ;;
esac
STUB
chmod +x "${TMP_DIR}/bin/sbatch"

run_stage3() {
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
    USER_EMAIL="test@example.invalid" PREPROCESS_SUBMITTER_TEST=1 \
    bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" "$@"
}

SELECTION="${TMP_DIR}/selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\nBreast_cancer\tbatch_effect_uncorrected\nCovid19_PBMC\tbatch_effect_uncorrected\nKidney_KPMP_full\tbatch_effect_uncorrected\nMyocardial_infarction\tbatch_effect_uncorrected\nDiabetes\tbatch_effect_uncorrected\nLupus_PBMC\tbatch_effect_uncorrected\nLung\tbatch_effect_uncorrected\nParkinson\tbatch_effect_uncorrected\nJoanito\tbatch_effect_uncorrected\nStephenson\tbatch_effect_uncorrected\nCombinedPBMC\tbatch_effect_uncorrected\n' > "${SELECTION}"

OUTPUT="$(run_stage3 --selection-file "${SELECTION}" --exact-batch-selection)"
case "${OUTPUT}" in *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;; *) echo "missing array marker" >&2; exit 1 ;; esac
case "${OUTPUT}" in *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) ;; *) echo "missing watchdog marker" >&2; exit 1 ;; esac
MANIFEST="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
[[ -s "${MANIFEST}" ]]
[[ "$(wc -l < "${MANIFEST}" | tr -d '[:space:]')" == 12 ]]
[[ "$(sed -n '1p' "${MANIFEST}")" == $'Alzheimer\tbatch_effect_uncorrected' ]]
[[ "$(sed -n '4p' "${MANIFEST}")" == $'Kidney_KPMP_full\tbatch_effect_uncorrected' ]]
[[ "$(sed -n '12p' "${MANIFEST}")" == $'CombinedPBMC\tbatch_effect_uncorrected' ]]
RUN_ROOT="$(dirname "${MANIFEST}")/.."
RUN_ROOT="$(cd "${RUN_ROOT}" && pwd)"
SCHEDULER_MANIFEST="${RUN_ROOT}/manifests/scheduler_ids.tsv"
[[ "$(wc -l < "${SCHEDULER_MANIFEST}" | tr -d '[:space:]')" == 2 ]]
[[ -s "${RUN_ROOT}/manifests/source.manifest" && -s "${RUN_ROOT}/manifests/runtime.identity" ]]
cmp -s "${SOURCE_MANIFEST}" "${RUN_ROOT}/manifests/source.manifest"
CALLS="$(cat "${CAPTURE}")"
case "${CALLS}" in *"--array=1-12%1000"*) ;; *) echo "array was not submitted with all selected rows" >&2; exit 1 ;; esac
case "${CALLS}" in *"--dependency=afterany:600001"*) ;; *) echo "watchdog dependency missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"PREPROCESS_SELECTION_FILE=${RUN_ROOT}/manifests/pending.tsv"*) ;; *) echo "pending manifest was not exported" >&2; exit 1 ;; esac
case "${CALLS}" in *"PREPROCESS_RUN_ROOT=${RUN_ROOT}"*) ;; *) echo "run root was not exported at scheduler boundary" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_ROOT=${SOURCE_ROOT}"*) ;; *) echo "snapshot source root export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}"*) ;; *) echo "snapshot source manifest export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_MANIFEST_RUN=${RUN_ROOT}/manifests/source.manifest"*) ;; *) echo "run-bound source manifest export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}"*) ;; *) echo "versioned runtime image export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}"*) ;; *) echo "versioned runtime manifest export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_IDENTITY=${RUN_ROOT}/manifests/runtime.identity"*) ;; *) echo "run-bound runtime identity export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"${SOURCE_ROOT}/src/3_scrnaseq_preprocessing/1.1_run_worker.sh"*) ;; *) echo "array did not use immutable worker script" >&2; exit 1 ;; esac
case "${CALLS}" in *"${SOURCE_ROOT}/src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh"*) ;; *) echo "watchdog did not use immutable watchdog script" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_MODE=apptainer"*"ECODA_RUNTIME_PROFILE=stage3"*) ;; *) echo "Stage 3 runtime export missing" >&2; exit 1 ;; esac

# Missing image identity fails before the first scheduler boundary.
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  ECODA_RUNTIME_IMAGE="${TMP_DIR}/missing.sif" \
  ECODA_RUNTIME_MANIFEST="${TMP_DIR}/missing.sif.manifest" \
  PREPROCESS_SUBMITTER_TEST=1 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
  --selection-file "${SELECTION}" --exact-batch-selection >/dev/null 2>&1; then
  echo "Stage 3 accepted missing immutable runtime image" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]

# Selection identity, schema, and semantic checks remain fail-closed.
BAD="${TMP_DIR}/bad.tsv"
for bad_kind in legacy corrected missing; do
  if [[ "${bad_kind}" == missing ]]; then
    sed '12d' "${SELECTION}" > "${BAD}"
  else
    bad_view="batch_effect_corrected"
    [[ "${bad_kind}" == legacy ]] && bad_view="batch_effect_analysis"
    sed "1s/batch_effect_uncorrected/${bad_view}/" "${SELECTION}" > "${BAD}"
  fi
  RUNS_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs"
  BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
  : > "${CAPTURE}"
  set +e
  run_stage3 --selection-file "${BAD}" --exact-batch-selection >/dev/null 2>&1
  RC=$?
  set -e
  [[ ${RC} -ne 0 ]]
  [[ ! -s "${CAPTURE}" ]]
  [[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]
done
BAD_NONEXACT="${TMP_DIR}/bad-nonexact.tsv"
printf 'Adams\tmissing_view\n' > "${BAD_NONEXACT}"
: > "${CAPTURE}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
set +e
run_stage3 --selection-file "${BAD_NONEXACT}" >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]

# A legacy run cannot enter sync/recovery without both run-bound manifests.
LEGACY_RUN="legacy-unpinned"
mkdir -p "${RUNS_ROOT}/${LEGACY_RUN}/manifests" "${RUNS_ROOT}/${LEGACY_RUN}/status" "${RUNS_ROOT}/${LEGACY_RUN}/logs"
printf 'STAGE=stage3\nRUN_ID=%s\nSTATE=ACTIVE\n' "${LEGACY_RUN}" > "${RUNS_ROOT}/${LEGACY_RUN}/metadata"
: > "${CAPTURE}"
set +e
run_stage3 --sync-only "${LEGACY_RUN}" >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
[[ ! -s "${CAPTURE}" ]]
[[ "$(sed -n 's/^STATE=//p' "${RUNS_ROOT}/${LEGACY_RUN}/status/terminal")" == FAIL ]]

# The shared source-script guard rejects an outside script before any sbatch.
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
export PROJECT_ROOT="${ROOT}" DATASETS_JSON_FILE="${ROOT}/datasets.json"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
OUTSIDE_SCRIPT="${TMP_DIR}/outside-worker.sh"
printf '#!/bin/bash\n' > "${OUTSIDE_SCRIPT}"
chmod +x "${OUTSIDE_SCRIPT}"
if ecoda_require_source_script_path "${OUTSIDE_SCRIPT}" "${SOURCE_ROOT}" >/dev/null 2>&1; then
  echo "outside Stage 3 script was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]

# Output ownership is path based: a same-run OOM retry revalidates an ACTIVE
# owner, while another run is rejected without reclaiming it.
rm -rf "${ECODA_OWNERS_ROOT}"
mkdir -p "${HPC_SCRATCH_DIR}/Adams/output"
OWNER_SELECTION="${TMP_DIR}/owner-selection.tsv"
printf 'Adams\tbenchmark_analysis\n' > "${OWNER_SELECTION}"
ecoda_init_run stage3 retry-owner >/dev/null
OWNER_RUN_ROOT="${ECODA_RUN_ROOT}"
ecoda_validate_output_ownership stage3 "${OWNER_SELECTION}" retry-owner
OWNER_PATH="${ECODA_OUTPUT_OWNER_DIRS[0]}"
ecoda_validate_output_ownership stage3 "${OWNER_SELECTION}" retry-owner
if ecoda_validate_output_ownership stage3 "${OWNER_SELECTION}" other-owner; then
  echo "another run bypassed an active Stage 3 output owner" >&2
  exit 1
fi
[[ "$(sed -n 's/^RUN_ID=//p' "${OWNER_PATH}/owner")" == retry-owner ]]
rm -rf "${ECODA_OWNERS_ROOT}"

# Exercise the watchdog's OOM retry boundary with the immutable source/runtime
# identity and a real same-run ACTIVE global output owner.
WD_RUN_ID="watchdog-retry"
ecoda_init_run stage3 "${WD_RUN_ID}" >/dev/null
WD_ROOT="${ECODA_RUN_ROOT}"
WD_OUTPUT_NAME="$(jq -r '.Adams.views.benchmark_analysis.output_file_name' "${ROOT}/datasets.json")"
WD_OUTPUT="${HPC_SCRATCH_DIR}/Adams/output/${WD_OUTPUT_NAME}"
printf 'watchdog-output\n' > "${WD_OUTPUT}"
write_sidecar "${WD_OUTPUT}"
ecoda_validate_checksum "${WD_OUTPUT}" >/dev/null
ecoda_write_artifact_record "${WD_OUTPUT}" stage3 "${WD_RUN_ID}" >/dev/null
printf 'Adams\tbenchmark_analysis\n' > "${WD_ROOT}/manifests/selection.tsv"
cp "${WD_ROOT}/manifests/selection.tsv" "${WD_ROOT}/manifests/pending.tsv"
ecoda_validate_output_ownership stage3 "${WD_ROOT}/manifests/pending.tsv" "${WD_RUN_ID}"
WD_STAGE_OWNER="$(ecoda_owner_acquire stage3 Adams/benchmark_analysis "${WD_RUN_ID}" 0 0)"
printf 'Adams/benchmark_analysis\t%s\n' "${WD_STAGE_OWNER}" > "${WD_ROOT}/manifests/owners.tsv"
printf 'RUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nRUNTIME_IMAGE_SHA256=%s\nRUNTIME_MANIFEST_SHA256=%s\nRUNTIME_IMAGE_SIZE=%s\nRUNTIME_MANIFEST_SIZE=%s\nIMAGE_PIXI_TOML_SHA256=%s\nIMAGE_PIXI_LOCK_SHA256=%s\n' \
  "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" "${RUNTIME_IMAGE_SHA256}" \
  "$(sha256_file "${RUNTIME_MANIFEST}")" "$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')" \
  "$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')" "${SOURCE_TOML_SHA256}" \
  "${SOURCE_LOCK_SHA256}" > "${WD_ROOT}/manifests/runtime.identity"
chmod 600 "${WD_ROOT}/manifests/runtime.identity"
cp "${SOURCE_MANIFEST}" "${WD_ROOT}/manifests/source.manifest"
WATCHDOG_BIN="${TMP_DIR}/watchdog-bin"
WATCHDOG_CAPTURE="${TMP_DIR}/watchdog.sbatch.calls"
mkdir -p "${WATCHDOG_BIN}"
export WATCHDOG_CAPTURE
cat > "${WATCHDOG_BIN}/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${WATCHDOG_CAPTURE}"
printf '930002\n'
STUB
cat > "${WATCHDOG_BIN}/sacct" <<'STUB'
#!/bin/bash
set -euo pipefail
case "$*" in
  *930001*) printf '930001_1|OUT_OF_MEMORY|0:0\n' ;;
  *930002*) printf '930002_1|COMPLETED|0:0\n' ;;
  *930003*) printf '930003_1|FAILED|1:0\n' ;;
  *) exit 1 ;;
esac
STUB
chmod +x "${WATCHDOG_BIN}/sbatch" "${WATCHDOG_BIN}/sacct"
WATCHDOG_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${WATCHDOG_BIN}:${TMP_DIR}/bin:${PATH}" \
  USER_EMAIL="test@example.invalid" HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" \
  NAS_TARGET_DIR="${NAS_TARGET_DIR}" ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" \
  ECODA_RUNTIME_MODE=host ECODA_RUNTIME_PROFILE=stage3 \
  ECODA_RUNTIME_IMAGE_SHA256="${RUNTIME_IMAGE_SHA256}" \
  ECODA_RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")" \
  ECODA_RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')" \
  ECODA_RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')" \
  PREPROCESS_RUN_ROOT="${WD_ROOT}" PREPROCESS_PENDING_MANIFEST="${WD_ROOT}/manifests/pending.tsv" \
  STAGE3_WATCHDOG_POLL_SECONDS=0 ECODA_ACCOUNTING_EMPTY_GRACE=1 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh" \
    "${WD_RUN_ID}" "${WD_ROOT}/manifests/selection.tsv" 930001 1G 4G shared-cpu 1
)"
case "${WATCHDOG_OUTPUT}" in *"Stage 3 watchdog completed for run ${WD_RUN_ID}"*) ;; *) echo "watchdog retry did not complete" >&2; exit 1 ;; esac
WATCHDOG_CALL="$(cat "${WATCHDOG_CAPTURE}")"
case "${WATCHDOG_CALL}" in *"${SOURCE_ROOT}/src/3_scrnaseq_preprocessing/1.1_run_worker.sh"*) ;; *) echo "OOM retry escaped immutable worker script" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"PREPROCESS_SELECTION_FILE=${WD_ROOT}/manifests/selection.retry_1.tsv"*) ;; *) echo "OOM retry manifest export missing" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"PREPROCESS_RUN_ROOT=${WD_ROOT}"*) ;; *) echo "OOM retry run root export missing" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"ECODA_SOURCE_ROOT=${SOURCE_ROOT}"*) ;; *) echo "OOM retry omitted source root" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}"*) ;; *) echo "OOM retry omitted source manifest" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"ECODA_SOURCE_SNAPSHOT_REQUIRED=1"*) ;; *) echo "OOM retry omitted snapshot requirement" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"ECODA_RUNTIME_IMAGE=${RUNTIME_IMAGE}"*) ;; *) echo "OOM retry omitted runtime image" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}"*) ;; *) echo "OOM retry omitted runtime manifest" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"ECODA_RUNTIME_IDENTITY=${WD_ROOT}/manifests/runtime.identity"*) ;; *) echo "OOM retry omitted runtime identity" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"ECODA_RUN_ID=${WD_RUN_ID}"*) ;; *) echo "OOM retry omitted run ID" >&2; exit 1 ;; esac
[[ "$(sed -n 's/^STATE=//p' "${WD_ROOT}/status/watchdog")" == OK ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_STAGE_OWNER}/owner")" == OK ]]

WD_SCRATCH_OWNER="$(ecoda_artifact_owner_dir "${WD_OUTPUT}")"
WD_NAS_OUTPUT="${NAS_TARGET_DIR}/Adams/output/${WD_OUTPUT_NAME}"
WD_NAS_OWNER="$(ecoda_artifact_owner_dir "${WD_NAS_OUTPUT}")"
[[ "$(sed -n 's/^STATE=//p' "${WD_SCRATCH_OWNER}/owner")" == OK ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_NAS_OWNER}/owner")" == OK ]]

# A submitter-preacquired same-run owner is also finalized FAIL when the
# watchdog reaches a terminal non-OOM scheduler failure.
FAIL_RUN_ID="watchdog-failure"
ecoda_init_run stage3 "${FAIL_RUN_ID}" >/dev/null
FAIL_ROOT="${ECODA_RUN_ROOT}"
FAIL_OUTPUT_NAME="$(jq -r '.Breast_cancer.views.batch_effect_uncorrected.output_file_name' "${ROOT}/datasets.json")"
FAIL_OUTPUT="${HPC_SCRATCH_DIR}/Breast_cancer/output/${FAIL_OUTPUT_NAME}"
mkdir -p "$(dirname "${FAIL_OUTPUT}")"
printf 'watchdog-failure-output\n' > "${FAIL_OUTPUT}"
printf 'Breast_cancer\tbatch_effect_uncorrected\n' > "${FAIL_ROOT}/manifests/selection.tsv"
cp "${FAIL_ROOT}/manifests/selection.tsv" "${FAIL_ROOT}/manifests/pending.tsv"
ecoda_validate_output_ownership stage3 "${FAIL_ROOT}/manifests/pending.tsv" "${FAIL_RUN_ID}"
FAIL_STAGE_OWNER="$(ecoda_owner_acquire \
  stage3 Breast_cancer/batch_effect_uncorrected "${FAIL_RUN_ID}" 0 0)"
printf 'Breast_cancer/batch_effect_uncorrected\t%s\n' "${FAIL_STAGE_OWNER}" \
  > "${FAIL_ROOT}/manifests/owners.tsv"
cp "${SOURCE_MANIFEST}" "${FAIL_ROOT}/manifests/source.manifest"
cp "${WD_ROOT}/manifests/runtime.identity" "${FAIL_ROOT}/manifests/runtime.identity"
chmod 600 "${FAIL_ROOT}/manifests/source.manifest" \
  "${FAIL_ROOT}/manifests/runtime.identity"
FAIL_SCRATCH_OWNER="$(ecoda_artifact_owner_dir "${FAIL_OUTPUT}")"
FAIL_NAS_OUTPUT="${NAS_TARGET_DIR}/Breast_cancer/output/${FAIL_OUTPUT_NAME}"
FAIL_NAS_OWNER="$(ecoda_artifact_owner_dir "${FAIL_NAS_OUTPUT}")"
set +e
FAIL_WATCHDOG_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${WATCHDOG_BIN}:${TMP_DIR}/bin:${PATH}" \
  USER_EMAIL="test@example.invalid" HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" \
  NAS_TARGET_DIR="${NAS_TARGET_DIR}" ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" \
  ECODA_RUNTIME_MODE=host ECODA_RUNTIME_PROFILE=stage3 \
  ECODA_RUNTIME_IMAGE_SHA256="${RUNTIME_IMAGE_SHA256}" \
  ECODA_RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")" \
  ECODA_RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')" \
  ECODA_RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')" \
  PREPROCESS_RUN_ROOT="${FAIL_ROOT}" \
  PREPROCESS_PENDING_MANIFEST="${FAIL_ROOT}/manifests/pending.tsv" \
  STAGE3_WATCHDOG_POLL_SECONDS=0 ECODA_ACCOUNTING_EMPTY_GRACE=1 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh" \
    "${FAIL_RUN_ID}" "${FAIL_ROOT}/manifests/selection.tsv" \
    930003 1G 4G shared-cpu 1
)"
FAIL_RC=$?
set -e
[[ ${FAIL_RC} -ne 0 ]]
[[ "$(sed -n 's/^STATE=//p' "${FAIL_ROOT}/status/watchdog")" == FAIL ]]
[[ "$(sed -n 's/^STATE=//p' "${FAIL_SCRATCH_OWNER}/owner")" == FAIL ]]
[[ "$(sed -n 's/^STATE=//p' "${FAIL_NAS_OWNER}/owner")" == FAIL ]]
[[ "$(sed -n 's/^STATE=//p' "${FAIL_STAGE_OWNER}/owner")" == FAIL ]]


# A missing existing output remains a compute selection and does not trigger a
# redundant H5AD preflight; the test-mode watchdog then fails closed without a
# fabricated terminal status.
MISSING_HOME="${TMP_DIR}/missing-home"
MISSING_HPC="${MISSING_HOME}/scratch/ECODA_paper"
mkdir -p "${MISSING_HPC}" "${MISSING_HOME}/logs"
: > "${CAPTURE}"
set +e
HOME="${MISSING_HOME}" PATH="${TMP_DIR}/bin:${PATH}" \
  HPC_SCRATCH_DIR="${MISSING_HPC}" ECODA_LOGS_DIR="${MISSING_HOME}/logs" \
  NAS_TARGET_DIR="${NAS_TARGET_DIR}" USER_EMAIL="test@example.invalid" \
  PREPROCESS_SUBMITTER_TEST=0 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
  --datasets Adams --views benchmark_analysis >/dev/null 2>&1
RC=$?
set -e
[[ ${RC} -ne 0 ]]
MISSING_CALLS="$(cat "${CAPTURE}")"
case "${MISSING_CALLS}" in *"--array=1-1%1000"*) ;; *) echo "missing-output path did not submit preprocessing work" >&2; exit 1 ;; esac
case "${MISSING_CALLS}" in *"h5ad_preflight"*) echo "missing-output path submitted an unnecessary preflight" >&2; exit 1 ;; esac

# Validator-only no-op still performs its strict preflight, then emits only the
# run-owned NOOP report and never submits a compute array/watchdog.
rm -rf "${ECODA_OWNERS_ROOT}"
NOOP_NAME="$(jq -r '.Adams.views.benchmark_analysis.output_file_name' "${SOURCE_ROOT}/datasets.json")"
NOOP_OUTPUT_PATH="${HPC_SCRATCH_DIR}/Adams/output/${NOOP_NAME}"
mkdir -p "$(dirname "${NOOP_OUTPUT_PATH}")"
chmod u+w "${NOOP_OUTPUT_PATH}"
printf 'already-validated\n' > "${NOOP_OUTPUT_PATH}"
write_sidecar "${NOOP_OUTPUT_PATH}"
NOOP_SELECTION="${TMP_DIR}/noop-selection.tsv"
printf 'Adams\tbenchmark_analysis\n' > "${NOOP_SELECTION}"
: > "${CAPTURE}"
NOOP_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" \
  USER_EMAIL="test@example.invalid" PREPROCESS_NOOP_PREFLIGHT=1 \
  PREPROCESS_SUBMITTER_TEST=0 ECODA_RUN_ID=noop-validated \
  ECODA_RUN_REPORT="${HPC_SCRATCH_DIR}/_ecoda_runs/noop-validated/status/noop-report" \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
    --selection-file "${NOOP_SELECTION}"
)"
chmod a-w "${NOOP_OUTPUT_PATH}"
NOOP_RUN_ID="$(printf '%s\n' "${NOOP_OUTPUT}" | sed -n 's/^PREPROCESS_RUN_ID=//p' | tail -1)"
[[ -n "${NOOP_RUN_ID}" ]]
NOOP_OWNER="$(ecoda_artifact_owner_dir "${NOOP_OUTPUT_PATH}")"
NOOP_NAS_OUTPUT="${NAS_TARGET_DIR}/Adams/output/${NOOP_NAME}"
NOOP_NAS_OWNER="$(ecoda_artifact_owner_dir "${NOOP_NAS_OUTPUT}")"
[[ "$(sed -n 's/^STATE=//p' "${NOOP_OWNER}/owner")" == OK ]]
[[ "$(sed -n 's/^RUN_ID=//p' "${NOOP_OWNER}/owner")" == "${NOOP_RUN_ID}" ]]
[[ "$(sed -n 's/^STATE=//p' "${NOOP_NAS_OWNER}/owner")" == OK ]]
[[ "$(sed -n 's/^RUN_ID=//p' "${NOOP_NAS_OWNER}/owner")" == "${NOOP_RUN_ID}" ]]
NOOP_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${NOOP_RUN_ID}"
case "${NOOP_OUTPUT}" in *"NOOP_VALIDATED=${NOOP_RUN_ID}"*) ;; *) echo "validator-only no-op marker missing" >&2; exit 1 ;; esac
[[ "$(sed -n 's/^STATE=//p' "${NOOP_ROOT}/status/noop-report")" == NOOP_VALIDATED ]]
NOOP_CALLS="$(cat "${CAPTURE}")"
[[ "$(printf '%s\n' "${NOOP_CALLS}" | wc -l | tr -d '[:space:]')" == 1 ]]
case "${NOOP_CALLS}" in *"h5ad_preflight_worker.sh"*) ;; *) echo "no-op did not run validator preflight" >&2; exit 1 ;; esac
case "${NOOP_CALLS}" in *"1.1_run_worker.sh"*|*"1.2_preprocess_watchdog.sh"*) echo "no-op submitted compute work" >&2; exit 1 ;; esac

echo "preprocessing stage submitter: OK"
