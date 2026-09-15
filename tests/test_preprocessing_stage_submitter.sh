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
  h5ad_preflight_worker.sh h5ad_obs_audit_worker.sh worker_retry.sh; do
  cp "${ROOT}/src/utils/bash/${source_file}" \
    "${SOURCE_ROOT}/src/utils/bash/${source_file}"
done
for source_file in 1.1_run_worker.sh 1.2_preprocess_watchdog.sh; do
  cp "${ROOT}/src/3_scrnaseq_preprocessing/${source_file}" \
    "${SOURCE_ROOT}/src/3_scrnaseq_preprocessing/${source_file}"
done
cp "${ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
  "${SOURCE_ROOT}/src/utils/py/benchmark_h5ad_contract.py"
cp "${ROOT}/src/utils/py/artifact_contract.py" \
  "${SOURCE_ROOT}/src/utils/py/artifact_contract.py"
cp "${ROOT}/src/utils/py/batch_contract.py" \
  "${SOURCE_ROOT}/src/utils/py/batch_contract.py"
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
export ECODA_SCRATCH_ROOT="${HPC_SCRATCH_DIR}"
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

# The scheduler stub records every boundary.  The generic preflight branch is
# retained only to catch an accidental validator submission in valid-only runs.
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail

worker_script=""
export_arg=""
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

# A batch-effect selection with pending Covid work submits the read-only obs
# audit. Keep that boundary scheduler-free while emitting the same run-owned
# evidence that the submitter validates.
if [[ "${worker_script}" == *h5ad_obs_audit_worker.sh ]]; then
  preflight_manifest="$(export_field H5AD_PREFLIGHT_MANIFEST)"
  preflight_status_dir="$(export_field H5AD_PREFLIGHT_STATUS_DIR)"
  preflight_run_root="$(export_field H5AD_PREFLIGHT_RUN_ROOT)"
  preflight_run_id="$(export_field H5AD_PREFLIGHT_RUN_ID)"
  source_root="$(export_field ECODA_SOURCE_ROOT)"
  source_manifest="$(export_field ECODA_SOURCE_MANIFEST)"
  runtime_identity="$(export_field ECODA_RUNTIME_IDENTITY || true)"
  runtime_manifest="$(export_field ECODA_RUNTIME_MANIFEST)"
  runtime_image="$(export_field ECODA_RUNTIME_IMAGE)"
  scratch_root="$(export_field HPC_SCRATCH_DIR)"
  input_path="${scratch_root}/Covid19_PBMC/data/Covid19_Ren2021.h5ad"
  input_dir="$(cd "$(dirname "${input_path}")" && pwd -P)" || exit 1
  input_path="${input_dir}/$(basename "${input_path}")"
  [[ -n "${input_path}" ]] || exit 1
  runtime_identity="${runtime_identity:-${preflight_run_root}/manifests/runtime.identity}"
  mkdir -p "${preflight_status_dir}" "${preflight_run_root}/preflight"
  actual_md5="$(
    if command -v md5sum >/dev/null 2>&1; then
      md5sum "${input_path}" | cut -d' ' -f1
    else
      md5 -q "${input_path}"
    fi
  )"
  actual_size="$(wc -c < "${input_path}" | tr -d '[:space:]')"
  task_id=0
  while IFS=$'\t' read -r dataset view manifest_input; do
    task_id=$((task_id + 1))
    report="${preflight_run_root}/preflight/Covid19_PBMC_${view}.json"
    subset_rule="$(jq -cS --arg view "${view}" \
      '.Covid19_PBMC.views[$view].subset_vars' \
      "${source_root}/datasets.json")"
    jq -n \
      --arg dataset "${dataset}" --arg view "${view}" \
      --arg input "${input_path}" --arg md5 "${actual_md5}" \
      --argjson size "${actual_size}" --argjson subset_rule "${subset_rule}" \
      --arg source_root "${source_root}" \
      --arg source_manifest "${source_manifest}" \
      --arg runtime_identity "${runtime_identity}" \
      --arg runtime_manifest "${runtime_manifest}" \
      --arg runtime_image "${runtime_image}" '
      {
        obs_only: true,
        view: $view,
        datasets: [{
          dataset: $dataset,
          view: $view,
          sample_column: "sampleID",
          input_identity: {
            path: $input,
            md5: $md5,
            size: $size
          },
          split_sample_count: 0,
          configured_cardinalities: {sampleID: 1, PatientID: 1},
          sampling_day_audit: {
            column: "Sampling day (Days after symptom onset)",
            raw_unique_values: ["29", "30", "30.5", "control", "unknown", "malformed"]
          },
          subset_vars: $subset_rule
        }],
        provenance: {
          source_root: $source_root,
          source_manifest: {path: $source_manifest},
          runtime_identity: {path: $runtime_identity},
          runtime_manifest: {path: $runtime_manifest},
          runtime_image: {path: $runtime_image}
        }
      }' > "${report}"
    report_md5="$(
      if command -v md5sum >/dev/null 2>&1; then
        md5sum "${report}" | cut -d' ' -f1
      else
        md5 -q "${report}"
      fi
    )"
    printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' \
      "${report_md5}" "$(wc -c < "${report}" | tr -d '[:space:]')" "${report}" \
      > "${report}.md5"
    printf 'STATE=OK\nRUN_ID=%s\nDATASET=%s\nVIEW=%s\nTASK_ID=%s\nINPUT_FILE=%s\nREPORT=%s\n' \
      "${preflight_run_id}" "${dataset}" "${view}" "${task_id}" \
      "${manifest_input}" "${report}" \
      > "${preflight_status_dir}/${dataset}__${view}.status"
  done < "${preflight_manifest}"
  printf '600003\n'
  exit 0
fi

printf '%s\n' "$*" >> "${CAPTURE}"
CALL_COUNT="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
if [[ "${PREPROCESS_NOOP_PREFLIGHT:-0}" == "1" &&
      "${worker_script}" == *h5ad_preflight_worker.sh ]]; then
  preflight_manifest=""
  preflight_status_dir=""
  preflight_run_root=""
  preflight_run_id=""
  source_root=""
  source_manifest=""
  host_prefix=""
  host_python_sha=""
  host_rscript_sha=""
  runtime_image=""
  runtime_manifest=""
  scratch_root=""
  logs_root=""
  for arg in "$@"; do
    case "${arg}" in
      --export=*) export_arg="${arg#--export=}" ;;
    esac
  done
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

# The exact historical run leaves its fixture owners active; discard only
# those temporary owners before the independent batch-effect arrays.
ECODA_OWNERS_ROOT="${HPC_SCRATCH_DIR}/_ecoda_owners"
rm -rf "${ECODA_OWNERS_ROOT}"

# The approved batch-effect release uses two independent arrays.  Each array
# runs the direct Covid obs-only preflight when its Covid row needs compute.
RUNS_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs"
COVID_INPUT="${HPC_SCRATCH_DIR}/Covid19_PBMC/data/Covid19_Ren2021.h5ad"
mkdir -p "$(dirname "${COVID_INPUT}")"
printf 'stub-direct-h5ad\n' > "${COVID_INPUT}"
COVID_INPUT_DIR="$(cd "$(dirname "${COVID_INPUT}")" && pwd -P)"
COVID_INPUT_CANONICAL="${COVID_INPUT_DIR}/$(basename "${COVID_INPUT}")"

UNCORRECTED_SELECTION="${TMP_DIR}/uncorrected-selection.tsv"
printf '%s\n' \
  'Covid19_PBMC	batch_effect_uncorrected' \
  'Diabetes	batch_effect_uncorrected' \
  'Joanito	batch_effect_uncorrected' \
  'Lung	batch_effect_uncorrected' > "${UNCORRECTED_SELECTION}"
: > "${CAPTURE}"
UNCORRECTED_OUTPUT="$(
  run_stage3 --selection-file "${UNCORRECTED_SELECTION}"
)"
case "${UNCORRECTED_OUTPUT}" in
  *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;;
  *) echo "uncorrected selection did not submit its own array" >&2; exit 1 ;;
esac
case "${UNCORRECTED_OUTPUT}" in
  *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) ;;
  *) echo "uncorrected selection did not submit its own watchdog" >&2; exit 1 ;;
esac
UNCORRECTED_MANIFEST="$(printf '%s\n' "${UNCORRECTED_OUTPUT}" |
  sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
[[ -s "${UNCORRECTED_MANIFEST}" ]]
UNCORRECTED_RUN_ROOT="$(dirname "${UNCORRECTED_MANIFEST}")/.."
UNCORRECTED_RUN_ROOT="$(cd "${UNCORRECTED_RUN_ROOT}" && pwd)"
cmp -s "${UNCORRECTED_SELECTION}" "${UNCORRECTED_MANIFEST}"
cmp -s "${UNCORRECTED_SELECTION}" \
  "${UNCORRECTED_RUN_ROOT}/manifests/pending.tsv"
[[ "$(wc -l < "${UNCORRECTED_RUN_ROOT}/manifests/output_ownership.tsv" |
  tr -d '[:space:]')" == 4 ]]
case "$(cat "${CAPTURE}")" in
  *"--array=1-4%1000"*) ;;
  *) echo "uncorrected array escaped its four-row scope" >&2; exit 1 ;;
esac
PREFLIGHT_RECORDS="$(awk -F '	' '$1 == "PREFLIGHT" && $2 == "600003" { count++ } END { print count + 0 }' \
  "${UNCORRECTED_RUN_ROOT}/manifests/scheduler_ids.tsv")"
[[ "${PREFLIGHT_RECORDS}" == 1 ]]
for preflight_view in batch_effect_uncorrected batch_effect_corrected; do
  PREFLIGHT_REPORT="${UNCORRECTED_RUN_ROOT}/preflight/Covid19_PBMC_${preflight_view}.json"
  PREFLIGHT_STATUS="${UNCORRECTED_RUN_ROOT}/status/h5ad_obs_audit/Covid19_PBMC__${preflight_view}.status"
  [[ -s "${PREFLIGHT_REPORT}" && -s "${PREFLIGHT_REPORT}.md5" ]]
  [[ "$(sed -n 's/^MD5=//p' "${PREFLIGHT_REPORT}.md5" | sed -n '1p')" == "$(md5_file "${PREFLIGHT_REPORT}")" ]]
  [[ "$(sed -n 's/^SIZE=//p' "${PREFLIGHT_REPORT}.md5" | sed -n '1p')" == "$(wc -c < "${PREFLIGHT_REPORT}" | tr -d '[:space:]')" ]]
  [[ "$(sed -n 's/^PATH=//p' "${PREFLIGHT_REPORT}.md5" | sed -n '1p')" == "${PREFLIGHT_REPORT}" ]]
  [[ -s "${PREFLIGHT_STATUS}" ]]
  [[ "$(sed -n 's/^STATE=//p' "${PREFLIGHT_STATUS}" | sed -n '1p')" == OK ]]
  jq -e --arg view "${preflight_view}" \
    --arg source_root "${SOURCE_ROOT}" \
    --arg source_manifest "${SOURCE_MANIFEST}" \
    --arg input "${COVID_INPUT_CANONICAL}" \
    '.obs_only == true and
     .view == $view and
     (.datasets | length) == 1 and
     .datasets[0].dataset == "Covid19_PBMC" and
     .datasets[0].view == $view and
     .datasets[0].sample_column == "sampleID" and
     .datasets[0].input_identity.path == $input and
     .datasets[0].split_sample_count == 0 and
     .datasets[0].sampling_day_audit.column == "Sampling day (Days after symptom onset)" and
     (.datasets[0].sampling_day_audit.raw_unique_values | type) == "array" and
     (.provenance.source_root == $source_root) and
     (.provenance.source_manifest.path == $source_manifest) and
     (.provenance.runtime_identity.path | type) == "string" and
     (.provenance.runtime_manifest.path | type) == "string" and
     (.provenance.runtime_image.path | type) == "string"' \
    "${PREFLIGHT_REPORT}" >/dev/null
  expected_subset_rule="$(jq -cS --arg view "${preflight_view}" \
    '.Covid19_PBMC.views[$view].subset_vars' "${SOURCE_ROOT}/datasets.json")"
  jq -e --argjson expected_rule "${expected_subset_rule}" \
    '.datasets[0].subset_vars == $expected_rule' "${PREFLIGHT_REPORT}" >/dev/null
done
[[ ! -e "${UNCORRECTED_RUN_ROOT}/manifests/corrected_batch_preflight.tsv" ]]
[[ ! -e "${UNCORRECTED_RUN_ROOT}/manifests/corrected_batch_preflight" ]]
[[ ! -e "${UNCORRECTED_RUN_ROOT}/status/corrected_batch_preflight" ]]

# The corrected recovery is explicit and excludes Alzheimer.  Its pending
# Covid row requires the direct obs-only reports, never retired metadata/RDS
# preflight state.
rm -rf "${ECODA_OWNERS_ROOT}"
CORRECTED_DATASETS=(
  Joanito
  Stephenson
  Breast_cancer
  Covid19_PBMC
  Kidney_KPMP_full
  Diabetes
  Lupus_PBMC
  Lung
)
CORRECTED_SELECTION="${TMP_DIR}/corrected-selection.tsv"
for corrected_dataset in "${CORRECTED_DATASETS[@]}"; do
  printf '%s	batch_effect_corrected\n' "${corrected_dataset}"
done > "${CORRECTED_SELECTION}"
: > "${CAPTURE}"
CORRECTED_OUTPUT="$(
  run_stage3 --selection-file "${CORRECTED_SELECTION}" --corrected-recovery
)"
case "${CORRECTED_OUTPUT}" in
  *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;;
  *) echo "corrected selection did not submit its own array" >&2; exit 1 ;;
esac
case "${CORRECTED_OUTPUT}" in
  *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) ;;
  *) echo "corrected selection did not submit its own watchdog" >&2; exit 1 ;;
esac
CORRECTED_MANIFEST="$(printf '%s\n' "${CORRECTED_OUTPUT}" |
  sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
[[ -s "${CORRECTED_MANIFEST}" ]]
[[ "$(wc -l < "${CORRECTED_MANIFEST}" | tr -d '[:space:]')" == 8 ]]
CORRECTED_RUN_ROOT="$(dirname "${CORRECTED_MANIFEST}")/.."
CORRECTED_RUN_ROOT="$(cd "${CORRECTED_RUN_ROOT}" && pwd)"
cmp -s "${CORRECTED_SELECTION}" "${CORRECTED_MANIFEST}"
cmp -s "${CORRECTED_SELECTION}" \
  "${CORRECTED_RUN_ROOT}/manifests/pending.tsv"
[[ "$(wc -l < "${CORRECTED_RUN_ROOT}/manifests/output_ownership.tsv" |
  tr -d '[:space:]')" == 8 ]]
[[ "$(wc -l < "${CORRECTED_RUN_ROOT}/manifests/scheduler_ids.tsv" |
  tr -d '[:space:]')" == 3 ]]
CORRECTED_PREFLIGHT_RECORDS="$(awk -F '	' '$1 == "PREFLIGHT" && $2 == "600003" { count++ } END { print count + 0 }' \
  "${CORRECTED_RUN_ROOT}/manifests/scheduler_ids.tsv")"
[[ "${CORRECTED_PREFLIGHT_RECORDS}" == 1 ]]
for preflight_view in batch_effect_uncorrected batch_effect_corrected; do
  PREFLIGHT_REPORT="${CORRECTED_RUN_ROOT}/preflight/Covid19_PBMC_${preflight_view}.json"
  PREFLIGHT_STATUS="${CORRECTED_RUN_ROOT}/status/h5ad_obs_audit/Covid19_PBMC__${preflight_view}.status"
  [[ -s "${PREFLIGHT_REPORT}" && -s "${PREFLIGHT_REPORT}.md5" ]]
  [[ "$(sed -n 's/^MD5=//p' "${PREFLIGHT_REPORT}.md5" | sed -n '1p')" == "$(md5_file "${PREFLIGHT_REPORT}")" ]]
  [[ -s "${PREFLIGHT_STATUS}" ]]
  [[ "$(sed -n 's/^STATE=//p' "${PREFLIGHT_STATUS}" | sed -n '1p')" == OK ]]
done
[[ ! -e "${CORRECTED_RUN_ROOT}/manifests/corrected_batch_preflight.tsv" ]]
[[ ! -e "${CORRECTED_RUN_ROOT}/manifests/corrected_batch_preflight" ]]
[[ ! -e "${CORRECTED_RUN_ROOT}/status/corrected_batch_preflight" ]]
CORRECTED_CALLS="$(cat "${CAPTURE}")"
case "${CORRECTED_CALLS}" in
  *"--array=1-8%1000"*) ;;
  *) echo "corrected array escaped its eight-row scope" >&2; exit 1 ;;
esac
case "${CORRECTED_CALLS}" in
  *"audit_corrected_source_metadata"*|*"audit_corrected_h5ad_source"*|*"corrected_batch_preflight"*)
    echo "corrected selection invoked retired source metadata preflight" >&2
    exit 1
    ;;
esac
# The corrected-recovery scope also admits exactly one targeted Breast row.
# It keeps corrected processing state distinct from the historical eight-row
# recovery and does not trigger a Covid preflight for this one-row input.
rm -rf "${ECODA_OWNERS_ROOT}"
BREAST_TARGET_SELECTION="${TMP_DIR}/breast-target-selection.tsv"
printf 'Breast_cancer\tbatch_effect_corrected\n' > "${BREAST_TARGET_SELECTION}"
: > "${CAPTURE}"
BREAST_TARGET_OUTPUT="$(
  run_stage3 --selection-file "${BREAST_TARGET_SELECTION}" --corrected-recovery
)"
case "${BREAST_TARGET_OUTPUT}" in
  *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;;
  *) echo "targeted Breast corrected selection did not submit its own array" >&2; exit 1 ;;
esac
case "${BREAST_TARGET_OUTPUT}" in
  *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) ;;
  *) echo "targeted Breast corrected selection did not submit its own watchdog" >&2; exit 1 ;;
esac
BREAST_TARGET_MANIFEST="$(printf '%s\n' "${BREAST_TARGET_OUTPUT}" |
  sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
[[ -s "${BREAST_TARGET_MANIFEST}" ]]
[[ "$(wc -l < "${BREAST_TARGET_MANIFEST}" | tr -d '[:space:]')" == 1 ]]
[[ "$(cat "${BREAST_TARGET_MANIFEST}")" == $'Breast_cancer\tbatch_effect_corrected' ]]
BREAST_TARGET_RUN_ROOT="$(dirname "${BREAST_TARGET_MANIFEST}")/.."
BREAST_TARGET_RUN_ROOT="$(cd "${BREAST_TARGET_RUN_ROOT}" && pwd)"
cmp -s "${BREAST_TARGET_SELECTION}" "${BREAST_TARGET_MANIFEST}"
cmp -s "${BREAST_TARGET_SELECTION}" \
  "${BREAST_TARGET_RUN_ROOT}/manifests/pending.tsv"
[[ "$(wc -l < "${BREAST_TARGET_RUN_ROOT}/manifests/output_ownership.tsv" |
  tr -d '[:space:]')" == 1 ]]
[[ "$(sed -n 's/^SELECTION_CLASSIFICATION=//p' \
  "${BREAST_TARGET_RUN_ROOT}/metadata")" == "corrected_breast_targeted" ]]
BREAST_TARGET_CALLS="$(cat "${CAPTURE}")"
case "${BREAST_TARGET_CALLS}" in
  *"--array=1-1%1000"*) ;;
  *) echo "targeted Breast corrected array escaped its one-row scope" >&2; exit 1 ;;
esac
case "${BREAST_TARGET_CALLS}" in
  *"h5ad_obs_audit_worker.sh"*)
    echo "targeted Breast corrected selection unexpectedly submitted Covid obs preflight" >&2
    exit 1
    ;;
esac

# A one-row corrected selection without the explicit corrected-recovery scope
# remains rejected before run initialization or scheduler submission.
: > "${CAPTURE}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
if run_stage3 --selection-file "${BREAST_TARGET_SELECTION}" >/dev/null 2>&1; then
  echo "accidental one-row Breast corrected selection was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]

# The scoped selector also rejects a one-row shape that is not the approved
# corrected Breast row.
BREAST_BAD_SCOPE_SELECTION="${TMP_DIR}/breast-bad-scope-selection.tsv"
printf 'Breast_cancer\tbatch_effect_uncorrected\n' > "${BREAST_BAD_SCOPE_SELECTION}"
: > "${CAPTURE}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
if run_stage3 --selection-file "${BREAST_BAD_SCOPE_SELECTION}" \
    --corrected-recovery >/dev/null 2>&1; then
  echo "non-approved corrected-recovery shape was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]

# A combined launch is retired rather than silently broadening either array.
: > "${CAPTURE}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
if run_stage3 --selection-file "${UNCORRECTED_SELECTION}" \
    --combined-batch-selection >/dev/null 2>&1; then
  echo "retired combined Stage 3 selection was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]
CORRECTED_DATASETS_CSV="$(IFS=,; printf '%s' "${CORRECTED_DATASETS[*]}")"
: > "${CAPTURE}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
if run_stage3 --datasets "${CORRECTED_DATASETS_CSV}" \
    --views batch_effect_corrected >/dev/null 2>&1; then
  echo "generated eight-row corrected recovery bypassed explicit selection-file guard" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]

: > "${CAPTURE}"
BEFORE_RUNS="$(printf '%s\n' "${RUNS_ROOT}"/*)"
if run_stage3 >/dev/null 2>&1; then
  echo "default broad Stage 3 selection was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
[[ "${BEFORE_RUNS}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]


# Validator-only idempotency checks use the same run/ownership primitives as
# the submitter.  Source the immutable fixture copy before initializing runs.
source "${SOURCE_ROOT}/src/slurm_config.sh"
source "${SOURCE_ROOT}/src/utils/bash/ecoda_run_common.sh"
# Alzheimer follow-up first proves that a reviewed producer cannot authorize
# the raw donor-only input.  The producer/output fixture is intentionally
# complete; only the bound config is wrong.
ALZ_RAW_PRODUCER_RUN_ID="alzheimer-stage2-raw-fixture"
ALZ_RAW_INPUT_NAME="SEAAD_Alzheimer.h5ad"
ALZ_RAW_INPUT="${HPC_SCRATCH_DIR}/Alzheimer/data/${ALZ_RAW_INPUT_NAME}"
mkdir -p "$(dirname "${ALZ_RAW_INPUT}")"
printf 'stub-alzheimer-raw-h5ad\n' > "${ALZ_RAW_INPUT}"
(
  set -euo pipefail
  source "${SOURCE_ROOT}/src/slurm_config.sh"
  source "${SOURCE_ROOT}/src/utils/bash/ecoda_run_common.sh"
  ecoda_init_run stage2 "${ALZ_RAW_PRODUCER_RUN_ID}" >/dev/null
  ALZ_PRODUCER_ROOT="${ECODA_RUN_ROOT}"
  printf 'alzheimer_donor_assay\tfixture\t%s\t-\t-\n' "${ALZ_RAW_INPUT}" \
    > "${ALZ_PRODUCER_ROOT}/manifests/steps.tsv"
  ecoda_write_checksum "${ALZ_RAW_INPUT}" >/dev/null
  ecoda_artifact_owner_acquire "${ALZ_RAW_INPUT}" stage2 "${ALZ_RAW_PRODUCER_RUN_ID}" 0 0 0 >/dev/null
  ecoda_write_artifact_record "${ALZ_RAW_INPUT}" alzheimer_donor_assay \
    "${ALZ_RAW_PRODUCER_RUN_ID}" >/dev/null
  ecoda_artifact_owner_set_state "${ALZ_RAW_INPUT}" OK "fixture raw producer published" >/dev/null
  ecoda_set_run_state OK "fixture raw producer validated" >/dev/null
)
export STAGE3_INPUT_PRODUCER_RUN_ID="${ALZ_RAW_PRODUCER_RUN_ID}"
ALZ_UNCORRECTED_SELECTION="${TMP_DIR}/alzheimer-uncorrected-selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\n' > "${ALZ_UNCORRECTED_SELECTION}"
: > "${CAPTURE}"
if run_stage3 --selection-file "${ALZ_UNCORRECTED_SELECTION}" \
    --alzheimer-followup >/dev/null 2>&1; then
  echo "raw donor-only Alzheimer follow-up was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]

# Bind the acceptance cases to a separate immutable post-derivative snapshot.
DERIV_SNAPSHOT_COMMIT="0000000000000000000000000000000000000002"
DERIV_SNAPSHOT_ROOT="${TMP_DIR}/snapshots/${DERIV_SNAPSHOT_COMMIT}"
DERIV_SOURCE_ROOT="${DERIV_SNAPSHOT_ROOT}/tree"
DERIV_SOURCE_IDENTITY="${DERIV_SNAPSHOT_ROOT}/identity"
DERIV_SOURCE_MANIFEST="${DERIV_SOURCE_IDENTITY}/source.manifest"
DERIV_SOURCE_ARCHIVE="${DERIV_SOURCE_IDENTITY}/source.tar"
mkdir -p "${DERIV_SNAPSHOT_ROOT}" "${DERIV_SOURCE_IDENTITY}"
cp -R "${SOURCE_ROOT}" "${DERIV_SOURCE_ROOT}"
chmod -R u+w "${DERIV_SOURCE_ROOT}"
jq '
  .Alzheimer.columns.sample = "donor_id_assay" |
  .Alzheimer.views.batch_effect_uncorrected.input_file_name = "SEAAD_Alzheimer_donor_assay.h5ad" |
  .Alzheimer.views.batch_effect_corrected.input_file_name = "SEAAD_Alzheimer_donor_assay.h5ad"
' "${DERIV_SOURCE_ROOT}/datasets.json" > "${DERIV_SOURCE_ROOT}/datasets.json.build"
mv "${DERIV_SOURCE_ROOT}/datasets.json.build" "${DERIV_SOURCE_ROOT}/datasets.json"
tar -cf "${DERIV_SOURCE_ARCHIVE}" -C "${DERIV_SOURCE_ROOT}" .
DERIV_SOURCE_ARCHIVE_SHA256="$(sha256_file "${DERIV_SOURCE_ARCHIVE}")"
DERIV_SOURCE_CONFIG_SHA256="$(sha256_file "${DERIV_SOURCE_ROOT}/config_helper.R")"
DERIV_SOURCE_DATASETS_SHA256="$(sha256_file "${DERIV_SOURCE_ROOT}/datasets.json")"
DERIV_SOURCE_TOML_SHA256="$(sha256_file "${DERIV_SOURCE_ROOT}/pixi.toml")"
DERIV_SOURCE_LOCK_SHA256="$(sha256_file "${DERIV_SOURCE_ROOT}/pixi.lock")"
chmod -R a-w "${DERIV_SOURCE_ROOT}"
printf '%s\n' \
  'FORMAT=1' \
  "SOURCE_ROOT=${DERIV_SOURCE_ROOT}" \
  "SOURCE_COMMIT=${DERIV_SNAPSHOT_COMMIT}" \
  "SOURCE_ARCHIVE_PATH=${DERIV_SOURCE_ARCHIVE}" \
  "SOURCE_ARCHIVE_SHA256=${DERIV_SOURCE_ARCHIVE_SHA256}" \
  "CONFIG_HELPER_SHA256=${DERIV_SOURCE_CONFIG_SHA256}" \
  "DATASETS_SHA256=${DERIV_SOURCE_DATASETS_SHA256}" \
  "PIXI_TOML_SHA256=${DERIV_SOURCE_TOML_SHA256}" \
  "PIXI_LOCK_SHA256=${DERIV_SOURCE_LOCK_SHA256}" \
  "AUX_ROOT=${DERIV_SOURCE_ROOT}/aux" \
  'SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4' > "${DERIV_SOURCE_MANIFEST}"
touch "${DERIV_SNAPSHOT_ROOT}/COMPLETE"
chmod a-w "${DERIV_SOURCE_MANIFEST}" "${DERIV_SOURCE_ARCHIVE}" "${DERIV_SNAPSHOT_ROOT}/COMPLETE"
SOURCE_ROOT="${DERIV_SOURCE_ROOT}"
SOURCE_MANIFEST="${DERIV_SOURCE_MANIFEST}"
PROJECT_ROOT="${SOURCE_ROOT}"
DATASETS_JSON_FILE="${SOURCE_ROOT}/datasets.json"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
  PROJECT_ROOT DATASETS_JSON_FILE

# The derivative producer/output and both one-row follow-ups must now pass.
source "${SOURCE_ROOT}/src/slurm_config.sh"
source "${SOURCE_ROOT}/src/utils/bash/ecoda_run_common.sh"
ALZ_PRODUCER_RUN_ID="alzheimer-stage2-derivative-fixture"
ALZ_INPUT_NAME="$(jq -r '.Alzheimer.views.batch_effect_uncorrected.input_file_name' \
  "${SOURCE_ROOT}/datasets.json")"
ALZ_INPUT="${HPC_SCRATCH_DIR}/Alzheimer/data/${ALZ_INPUT_NAME}"
mkdir -p "$(dirname "${ALZ_INPUT}")"
printf 'stub-alzheimer-derivative-h5ad\n' > "${ALZ_INPUT}"
(
  set -euo pipefail
  source "${SOURCE_ROOT}/src/slurm_config.sh"
  source "${SOURCE_ROOT}/src/utils/bash/ecoda_run_common.sh"
  ecoda_init_run stage2 "${ALZ_PRODUCER_RUN_ID}" >/dev/null
  ALZ_PRODUCER_ROOT="${ECODA_RUN_ROOT}"
  printf 'alzheimer_donor_assay\tfixture\t%s\t-\t-\n' "${ALZ_INPUT}" \
    > "${ALZ_PRODUCER_ROOT}/manifests/steps.tsv"
  ecoda_write_checksum "${ALZ_INPUT}" >/dev/null
  ecoda_artifact_owner_acquire "${ALZ_INPUT}" stage2 "${ALZ_PRODUCER_RUN_ID}" 0 0 0 >/dev/null
  ecoda_write_artifact_record "${ALZ_INPUT}" alzheimer_donor_assay \
    "${ALZ_PRODUCER_RUN_ID}" >/dev/null
  ecoda_artifact_owner_set_state "${ALZ_INPUT}" OK "fixture derivative published" >/dev/null
  ecoda_set_run_state OK "fixture derivative validated" >/dev/null
)
export STAGE3_INPUT_PRODUCER_RUN_ID="${ALZ_PRODUCER_RUN_ID}"
: > "${CAPTURE}"
ALZ_UNCORRECTED_OUTPUT="$(
  run_stage3 --selection-file "${ALZ_UNCORRECTED_SELECTION}" --alzheimer-followup
)"
case "${ALZ_UNCORRECTED_OUTPUT}" in
  *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;;
  *) echo "Alzheimer uncorrected follow-up did not submit its array" >&2; exit 1 ;;
esac
case "${ALZ_UNCORRECTED_OUTPUT}" in
  *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) ;;
  *) echo "Alzheimer uncorrected follow-up did not submit its watchdog" >&2; exit 1 ;;
esac
ALZ_UNCORRECTED_MANIFEST="$(printf '%s\n' "${ALZ_UNCORRECTED_OUTPUT}" |
  sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
ALZ_UNCORRECTED_ROOT="$(dirname "${ALZ_UNCORRECTED_MANIFEST}")/.."
ALZ_UNCORRECTED_ROOT="$(cd "${ALZ_UNCORRECTED_ROOT}" && pwd)"
[[ "$(wc -l < "${ALZ_UNCORRECTED_MANIFEST}" | tr -d '[:space:]')" == 1 ]]
[[ "$(cat "${ALZ_UNCORRECTED_MANIFEST}")" == $'Alzheimer\tbatch_effect_uncorrected' ]]
[[ "$(wc -l < "${ALZ_UNCORRECTED_ROOT}/manifests/output_ownership.tsv" |
  tr -d '[:space:]')" == 1 ]]
ALZ_INPUT_RECORD="$(printf 'Alzheimer\tbatch_effect_uncorrected\t%s\tdonor_id_assay\t%s\n' \
  "${ALZ_INPUT}" "${ALZ_PRODUCER_RUN_ID}")"
[[ "$(cat "${ALZ_UNCORRECTED_ROOT}/manifests/input_ownership.tsv")" == "${ALZ_INPUT_RECORD}" ]]
[[ "$(sed -n 's/^INPUT_PATH=//p' "${ALZ_UNCORRECTED_ROOT}/metadata")" == "${ALZ_INPUT}" ]]
[[ "$(sed -n 's/^INPUT_SAMPLE_COLUMN=//p' "${ALZ_UNCORRECTED_ROOT}/metadata")" == "donor_id_assay" ]]
[[ "$(sed -n 's/^SELECTION_CLASSIFICATION=//p' \
  "${ALZ_UNCORRECTED_ROOT}/metadata")" == "alzheimer_followup" ]]
[[ "$(sed -n 's/^INPUT_PRODUCER_RUN_ID=//p' \
  "${ALZ_UNCORRECTED_ROOT}/metadata")" == "${ALZ_PRODUCER_RUN_ID}" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == 2 ]]
if awk '$0 ~ /h5ad_obs_audit_worker/ {found=1} END {exit found ? 0 : 1}' \
    "${CAPTURE}"; then
  echo "Alzheimer follow-up unexpectedly submitted Covid obs preflight" >&2
  exit 1
fi

ALZ_CORRECTED_SELECTION="${TMP_DIR}/alzheimer-corrected-selection.tsv"
printf 'Alzheimer\tbatch_effect_corrected\n' > "${ALZ_CORRECTED_SELECTION}"
: > "${CAPTURE}"
ALZ_CORRECTED_OUTPUT="$(
  run_stage3 --selection-file "${ALZ_CORRECTED_SELECTION}" \
    --alzheimer-follow-up-selection
)"
case "${ALZ_CORRECTED_OUTPUT}" in
  *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;;
  *) echo "Alzheimer corrected follow-up did not submit its array" >&2; exit 1 ;;
esac
case "${ALZ_CORRECTED_OUTPUT}" in
  *"PREPROCESS_WATCHDOG_JOB_ID=600002"*) ;;
  *) echo "Alzheimer corrected follow-up did not submit its watchdog" >&2; exit 1 ;;
esac
ALZ_CORRECTED_MANIFEST="$(printf '%s\n' "${ALZ_CORRECTED_OUTPUT}" |
  sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
ALZ_CORRECTED_ROOT="$(dirname "${ALZ_CORRECTED_MANIFEST}")/.."
ALZ_CORRECTED_ROOT="$(cd "${ALZ_CORRECTED_ROOT}" && pwd)"
[[ "$(cat "${ALZ_CORRECTED_MANIFEST}")" == $'Alzheimer\tbatch_effect_corrected' ]]
[[ "$(wc -l < "${ALZ_CORRECTED_ROOT}/manifests/output_ownership.tsv" |
  tr -d '[:space:]')" == 1 ]]
ALZ_CORRECTED_INPUT_RECORD="$(printf 'Alzheimer\tbatch_effect_corrected\t%s\tdonor_id_assay\t%s\n' \
  "${ALZ_INPUT}" "${ALZ_PRODUCER_RUN_ID}")"
[[ "$(cat "${ALZ_CORRECTED_ROOT}/manifests/input_ownership.tsv")" == "${ALZ_CORRECTED_INPUT_RECORD}" ]]
[[ "$(sed -n 's/^INPUT_PATH=//p' "${ALZ_CORRECTED_ROOT}/metadata")" == "${ALZ_INPUT}" ]]
[[ "$(sed -n 's/^INPUT_SAMPLE_COLUMN=//p' "${ALZ_CORRECTED_ROOT}/metadata")" == "donor_id_assay" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == 2 ]]

ALZ_COMBINED_SELECTION="${TMP_DIR}/alzheimer-combined-selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\nAlzheimer\tbatch_effect_corrected\n' \
  > "${ALZ_COMBINED_SELECTION}"
: > "${CAPTURE}"
if run_stage3 --selection-file "${ALZ_COMBINED_SELECTION}" >/dev/null 2>&1; then
  echo "combined Alzheimer follow-up was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]

: > "${CAPTURE}"
if ECODA_SOURCE_SNAPSHOT_REQUIRED=0 \
  run_stage3 --selection-file "${ALZ_UNCORRECTED_SELECTION}" \
    --alzheimer-followup >/dev/null 2>&1; then
  echo "snapshot-unbound Alzheimer follow-up was accepted" >&2
  exit 1
fi
[[ ! -s "${CAPTURE}" ]]
unset STAGE3_INPUT_PRODUCER_RUN_ID


# Valid rows are published into the new run's artifact manifest and skipped;
# a mixed selection submits only its missing/invalid row.
rm -rf "${ECODA_OWNERS_ROOT}"
VALID_RUN_ID="prior-stage3-valid"
ecoda_init_run stage3 "${VALID_RUN_ID}" >/dev/null
VALID_NAME="$(jq -r '.Adams.views.benchmark_analysis.output_file_name' \
  "${SOURCE_ROOT}/datasets.json")"
VALID_SCRATCH="${HPC_SCRATCH_DIR}/Adams/output/${VALID_NAME}"
VALID_NAS="${NAS_TARGET_DIR}/Adams/output/${VALID_NAME}"
mkdir -p "$(dirname "${VALID_SCRATCH}")" "$(dirname "${VALID_NAS}")"
printf 'validated-output\n' > "${VALID_SCRATCH}"
printf 'validated-output\n' > "${VALID_NAS}"
write_sidecar "${VALID_SCRATCH}"
write_sidecar "${VALID_NAS}"
ecoda_write_artifact_record "${VALID_SCRATCH}" stage3 "${VALID_RUN_ID}" >/dev/null
ecoda_write_artifact_record "${VALID_NAS}" stage3 "${VALID_RUN_ID}" >/dev/null
ecoda_artifact_owner_acquire "${VALID_SCRATCH}" stage3 "${VALID_RUN_ID}" 0 1 1 >/dev/null
ecoda_artifact_owner_set_state "${VALID_SCRATCH}" OK "fixture validated" >/dev/null
ecoda_artifact_owner_acquire "${VALID_NAS}" stage3 "${VALID_RUN_ID}" 0 1 1 >/dev/null
ecoda_artifact_owner_set_state "${VALID_NAS}" OK "fixture validated" >/dev/null

NOOP_SELECTION="${TMP_DIR}/noop-selection.tsv"
printf 'Adams	benchmark_analysis\n' > "${NOOP_SELECTION}"
: > "${CAPTURE}"
NOOP_OUTPUT="$(run_stage3 --selection-file "${NOOP_SELECTION}")"
NOOP_RUN_ID="$(printf '%s\n' "${NOOP_OUTPUT}" |
  sed -n 's/^PREPROCESS_RUN_ID=//p' | tail -1)"
[[ -n "${NOOP_RUN_ID}" ]]
case "${NOOP_OUTPUT}" in
  *"NOOP_VALIDATED=${NOOP_RUN_ID}"*) ;;
  *) echo "valid Stage 3 row was not validator-only skipped" >&2; exit 1 ;;
esac
[[ ! -s "${CAPTURE}" ]]
NOOP_ROOT="${RUNS_ROOT}/${NOOP_RUN_ID}"
[[ "$(sed -n 's/^STATE=//p' "${NOOP_ROOT}/status/terminal" | sed -n '1p')" == NOOP_VALIDATED ]]
[[ -s "${NOOP_ROOT}/manifests/artifacts/"*.record ]]
NOOP_SCRATCH_OWNER="$(ecoda_artifact_owner_dir "${VALID_SCRATCH}")"
NOOP_NAS_OWNER="$(ecoda_artifact_owner_dir "${VALID_NAS}")"
[[ "$(sed -n 's/^STATE=//p' "${NOOP_SCRATCH_OWNER}/owner" | sed -n '1p')" == OK ]]
[[ "$(sed -n 's/^RUN_ID=//p' "${NOOP_SCRATCH_OWNER}/owner" | sed -n '1p')" == "${NOOP_RUN_ID}" ]]
[[ "$(sed -n 's/^STATE=//p' "${NOOP_NAS_OWNER}/owner" | sed -n '1p')" == OK ]]
[[ "$(sed -n 's/^RUN_ID=//p' "${NOOP_NAS_OWNER}/owner" | sed -n '1p')" == "${NOOP_RUN_ID}" ]]

rm -rf "${ECODA_OWNERS_ROOT}"
ecoda_artifact_owner_acquire "${VALID_SCRATCH}" stage3 "${VALID_RUN_ID}" 0 1 1 >/dev/null
ecoda_artifact_owner_set_state "${VALID_SCRATCH}" OK "fixture validated" >/dev/null
ecoda_artifact_owner_acquire "${VALID_NAS}" stage3 "${VALID_RUN_ID}" 0 1 1 >/dev/null
ecoda_artifact_owner_set_state "${VALID_NAS}" OK "fixture validated" >/dev/null
MIXED_SELECTION="${TMP_DIR}/mixed-selection.tsv"
printf '%s\n' \
  'Adams	benchmark_analysis' \
  'Bassez	benchmark_analysis' > "${MIXED_SELECTION}"
: > "${CAPTURE}"
MIXED_OUTPUT="$(run_stage3 --selection-file "${MIXED_SELECTION}")"
MIXED_MANIFEST="$(printf '%s\n' "${MIXED_OUTPUT}" |
  sed -n 's/^PREPROCESS_DATASET_MANIFEST=//p')"
MIXED_ROOT="$(dirname "${MIXED_MANIFEST}")/.."
MIXED_ROOT="$(cd "${MIXED_ROOT}" && pwd)"
[[ "$(wc -l < "${MIXED_ROOT}/manifests/pending.tsv" |
  tr -d '[:space:]')" == 1 ]]
[[ "$(sed -n '1p' "${MIXED_ROOT}/manifests/pending.tsv")" == $'Bassez\tbenchmark_analysis' ]]
case "${MIXED_OUTPUT}" in
  *"PREPROCESS_ARRAY_JOB_ID=600001"*) ;;
  *) echo "mixed Stage 3 selection did not submit missing work" >&2; exit 1 ;;
esac
case "$(cat "${CAPTURE}")" in
  *"--array=1-1%1000"*) ;;
  *) echo "mixed Stage 3 array was not narrowed to invalid rows" >&2; exit 1 ;;
esac
[[ "$(printf '%s\n' "$(cat "${CAPTURE}")" | wc -l |
  tr -d '[:space:]')" == 2 ]]

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
CORRECTED_BATCH_BASE="${TMP_DIR}/corrected-batch-base.json"
cat > "${CORRECTED_BATCH_BASE}" <<'JSON'
{"Fixture":{"columns":{"sample":"sample_id","label":"label","batch":"batch_a"},"views":{"batch_effect_corrected":{}}}}
JSON
CORRECTED_BATCH_STATE="${TMP_DIR}/corrected-batch-state"
CORRECTED_BATCH_RUNS="${CORRECTED_BATCH_STATE}/_ecoda_runs"
mkdir -p "${CORRECTED_BATCH_RUNS}"
: > "${CAPTURE}"
for corrected_case in null object empty duplicate blank label_overlap sample_overlap reserved_sample; do
  case "${corrected_case}" in
    null) raw_batch='null' ;;
    object) raw_batch='{"name":"batch_a"}' ;;
    empty) raw_batch='[]' ;;
    duplicate) raw_batch='["batch_a","batch_a"]' ;;
    blank) raw_batch='["   "]' ;;
    label_overlap) raw_batch='["label"]' ;;
    sample_overlap) raw_batch='["sample_id"]' ;;
    reserved_sample) raw_batch='["Sample"]' ;;
  esac
  corrected_fixture="${TMP_DIR}/corrected-${corrected_case}.json"
  jq --argjson batch "${raw_batch}" \
    '.Fixture.columns.batch = $batch' "${CORRECTED_BATCH_BASE}" \
    > "${corrected_fixture}"
  corrected_run_root="${CORRECTED_BATCH_RUNS}/${corrected_case}"
  corrected_runs_before="$(printf '%s\n' "${CORRECTED_BATCH_RUNS}"/*)"
  stage3_runs_before="$(printf '%s\n' "${RUNS_ROOT}"/*)"
  owners_before="$(printf '%s\n' "${ECODA_OWNERS_ROOT}"/*)"
  if ecoda_validate_corrected_batch_columns \
      "${corrected_fixture}" Fixture batch_effect_corrected >/dev/null 2>&1; then
    RC=0
  else
    RC=$?
  fi
  [[ ${RC} -ne 0 ]]
  [[ ! -s "${CAPTURE}" ]]
  [[ "${corrected_runs_before}" == "$(printf '%s\n' "${CORRECTED_BATCH_RUNS}"/*)" ]]
  [[ "${stage3_runs_before}" == "$(printf '%s\n' "${RUNS_ROOT}"/*)" ]]
  [[ "${owners_before}" == "$(printf '%s\n' "${ECODA_OWNERS_ROOT}"/*)" ]]
  [[ ! -e "${corrected_run_root}" ]]
  [[ ! -e "${corrected_run_root}/manifests/scheduler_ids.tsv" ]]
  [[ ! -e "${corrected_run_root}/manifests/selection.tsv" ]]
  [[ ! -e "${corrected_run_root}/manifests/pending.tsv" ]]
  [[ ! -e "${corrected_run_root}/status" ]]
done


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
chmod u+w "${VALID_SCRATCH}" "${VALID_NAS}"

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
WD_SCRATCH_OWNER="$(ecoda_artifact_owner_dir "${WD_OUTPUT}")"
WD_NAS_OUTPUT="${NAS_TARGET_DIR}/Adams/output/${WD_OUTPUT_NAME}"
WD_NAS_OWNER="$(ecoda_artifact_owner_dir "${WD_NAS_OUTPUT}")"
[[ -d "${WD_SCRATCH_OWNER}" && -d "${WD_NAS_OWNER}" ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_SCRATCH_OWNER}/owner" | sed -n '1p')" == ACTIVE ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_NAS_OWNER}/owner" | sed -n '1p')" == ACTIVE ]]
WD_STAGE_OWNER="$(ecoda_owner_acquire stage3 Adams/benchmark_analysis "${WD_RUN_ID}" 0 0)"
printf 'Adams/benchmark_analysis\t%s\n' "${WD_STAGE_OWNER}" > "${WD_ROOT}/manifests/owners.tsv"
ecoda_write_checksum "${WD_ROOT}/manifests/selection.tsv" >/dev/null
ecoda_write_checksum "${WD_ROOT}/manifests/pending.tsv" >/dev/null
printf 'ARRAY\t930001\nWATCHDOG\t930002\n' > "${WD_ROOT}/manifests/scheduler_ids.tsv"
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
# The watchdog only validates scratch outputs.  Owners remain ACTIVE until
# the submitter's existing verified sync/finalization path runs.
[[ -d "${WD_SCRATCH_OWNER}" && -d "${WD_NAS_OWNER}" ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_SCRATCH_OWNER}/owner" | sed -n '1p')" == ACTIVE ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_NAS_OWNER}/owner" | sed -n '1p')" == ACTIVE ]]
WD_SYNC_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${WATCHDOG_BIN}:${TMP_DIR}/bin:${PATH}" \
  USER_EMAIL="test@example.invalid" HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" \
  NAS_TARGET_DIR="${NAS_TARGET_DIR}" ECODA_LOGS_DIR="${ECODA_LOGS_DIR}" \
  ECODA_RUNTIME_MODE=host ECODA_RUNTIME_PROFILE=stage3 \
  ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
  ECODA_RUNTIME_IMAGE_SHA256="${RUNTIME_IMAGE_SHA256}" \
  ECODA_RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")" \
  ECODA_RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')" \
  ECODA_RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')" \
  PREPROCESS_SUBMITTER_TEST=0 \
  bash "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
    --sync-only "${WD_RUN_ID}"
)"
case "${WD_SYNC_OUTPUT}" in
  *"PREPROCESS_RUN_ID=${WD_RUN_ID}"*) ;;
  *) echo "verified Stage 3 sync-only path did not complete" >&2; exit 1 ;;
esac
[[ -d "${WD_SCRATCH_OWNER}" && -d "${WD_NAS_OWNER}" ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_SCRATCH_OWNER}/owner" | sed -n '1p')" == OK ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_NAS_OWNER}/owner" | sed -n '1p')" == OK ]]
[[ "$(sed -n 's/^STATE=//p' "${WD_STAGE_OWNER}/owner" | sed -n '1p')" == OK ]]

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
[[ -d "${FAIL_SCRATCH_OWNER}" && -d "${FAIL_NAS_OWNER}" ]]
[[ "$(sed -n 's/^STATE=//p' "${FAIL_SCRATCH_OWNER}/owner" | sed -n '1p')" == ACTIVE ]]
[[ "$(sed -n 's/^STATE=//p' "${FAIL_NAS_OWNER}/owner" | sed -n '1p')" == ACTIVE ]]
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
[[ -d "${FAIL_SCRATCH_OWNER}" && -d "${FAIL_NAS_OWNER}" ]]
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
  HPC_SCRATCH_DIR="${MISSING_HPC}" ECODA_SCRATCH_ROOT="${MISSING_HPC}" \
  ECODA_LOGS_DIR="${MISSING_HOME}/logs" \
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


echo "preprocessing stage submitter: OK"
