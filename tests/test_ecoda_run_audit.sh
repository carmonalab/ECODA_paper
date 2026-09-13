#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-run-audit.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  if [[ -d "${TMP_DIR}/snapshot" ]]; then
    chmod -R u+w "${TMP_DIR}/snapshot" >/dev/null 2>&1 || true
  fi
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PIXI_PYTHON="${ROOT}/.pixi/envs/default/bin/python"
export PIXI_PYTHON

SCRATCH="${TMP_DIR}/scratch"
NAS="${TMP_DIR}/nas"
AUDIT_TMP="${TMP_DIR}/tmp"
AUDIT_LOG_DIR="${TMP_DIR}/audit-logs"
mkdir -p "${SCRATCH}/_ecoda_runs" "${SCRATCH}/_ecoda_owners" "${NAS}" \
  "${AUDIT_TMP}" "${AUDIT_LOG_DIR}"
export HPC_SCRATCH_DIR="${SCRATCH}" NAS_TARGET_DIR="${NAS}"
export ECODA_OWNERS_ROOT="${SCRATCH}/_ecoda_owners"
export ECODA_RUNS_ROOT="${SCRATCH}/_ecoda_runs"
export TMPDIR="${AUDIT_TMP}" ECODA_LOGS_DIR="${AUDIT_LOG_DIR}"
unset ANALYSIS_VARIANT ANALYSIS_PASS ANALYSIS_ROOT ANALYSIS_NAS_ROOT \
  ANALYSIS_LOG_PREFIX PASS_ARG ECODA_STAGE5_LEGACY_SYNC_SKIP \
  ECODA_STAGE5_LEGACY_SYNC_SKIP_METHODS BENCHMARK_MATRIX_TEST \
  ECODA_SOURCE_ROOT ECODA_SOURCE_MANIFEST ECODA_SOURCE_SNAPSHOT_REQUIRED \
  ECODA_RUNTIME_IDENTITY ECODA_RUNTIME_IMAGE ECODA_RUNTIME_MANIFEST
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

make_artifact() {
  local path="$1" record="$2" owner
  mkdir -p "$(dirname "${path}")"
  printf 'fixture artifact\n' > "${path}"
  ecoda_write_checksum "${path}" >/dev/null
  owner="$(ecoda_artifact_owner_acquire "${path}" stage5 audit-producer 0 1 0)"
  ecoda_artifact_owner_set_state "${path}" OK "audit fixture" >/dev/null
  if [[ "${record}" == 1 ]]; then
    ecoda_write_artifact_record "${path}" mrvi audit-producer >/dev/null
  fi
}

make_run_identities() {
  local run_root="$1" source_manifest="$2" runtime_identity="$3"
  mkdir -p "${run_root}/manifests"
  cp "${source_manifest}" "${run_root}/manifests/source.manifest"
  cp "${runtime_identity}" "${run_root}/manifests/runtime.identity"
  chmod a-w "${run_root}/manifests/source.manifest" \
    "${run_root}/manifests/runtime.identity"
}

make_source_fixture() {
  local snapshot commit tree aux identity source_tar source_manifest
  local source_archive_sha config_sha helper_sha toml_sha lock_sha
  commit=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
  snapshot="${TMP_DIR}/snapshot/${commit}"
  tree="${snapshot}/tree"
  aux="${tree}/aux"
  identity="${snapshot}/identity"
  source_tar="${identity}/source.tar"
  source_manifest="${identity}/source.manifest"
  mkdir -p "${tree}/src/utils/py" "${tree}/src/5_run_benchmark_methods" \
    "${aux}" "${identity}"
  printf 'fixture helper\n' > "${tree}/config_helper.R"
  cat > "${tree}/datasets.json" <<'JSON'
{
  "Fixture": {
    "columns": {"sample": "Sample", "batch": "batch"},
    "views": {
      "batch_effect_corrected": {
        "output_file_name": "Fixture_batch_effect_analysis_corrected.h5ad",
        "columns": {"batch": "batch"}
      },
      "batch_effect_uncorrected": {
        "output_file_name": "Fixture_batch_effect_analysis_uncorrected.h5ad"
      }
    }
  }
}
JSON
  printf 'fixture pixi\n' > "${tree}/pixi.toml"
  printf 'fixture lock\n' > "${tree}/pixi.lock"
  cp "${ROOT}/src/slurm_config.sh" "${tree}/src/slurm_config.sh"
  cp "${ROOT}/src/utils/py/batch_contract.py" \
    "${tree}/src/utils/py/batch_contract.py"
  printf '# fixture matrix validator\n' > \
    "${tree}/src/5_run_benchmark_methods/matrix_artifact_validator.py"
  printf '# fixture h5ad validator\n' > \
    "${tree}/src/utils/py/benchmark_h5ad_contract.py"
  printf 'COMPLETE\n' > "${snapshot}/COMPLETE"
  for aux_file in scGateDB.rds genes.blocklist.rds EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
    printf 'fixture aux\n' > "${aux}/${aux_file}"
  done
  (cd "${tree}" && tar -cf "${source_tar}" .)
  source_archive_sha="$(sha256_file "${source_tar}")"
  config_sha="$(sha256_file "${tree}/datasets.json")"
  helper_sha="$(sha256_file "${tree}/config_helper.R")"
  toml_sha="$(sha256_file "${tree}/pixi.toml")"
  lock_sha="$(sha256_file "${tree}/pixi.lock")"
  cat > "${source_manifest}" <<EOF
FORMAT=1
SOURCE_ROOT=${tree}
SOURCE_COMMIT=${commit}
SOURCE_ARCHIVE_PATH=${source_tar}
SOURCE_ARCHIVE_SHA256=${source_archive_sha}
CONFIG_HELPER_SHA256=${helper_sha}
DATASETS_SHA256=${config_sha}
PIXI_TOML_SHA256=${toml_sha}
PIXI_LOCK_SHA256=${lock_sha}
AUX_ROOT=${aux}
SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4
EOF
  chmod -R a-w "${snapshot}"
  SOURCE_CONFIG="${tree}/datasets.json"
  SOURCE_MANIFEST="${source_manifest}"
  SOURCE_TREE="${tree}"
  export SOURCE_CONFIG SOURCE_MANIFEST SOURCE_TREE
}

make_runtime_fixture() {
  local runtime_dir="${TMP_DIR}/_ecoda_runtime/run-audit" image manifest identity
  local image_sha manifest_sha image_size manifest_size toml_sha lock_sha
  mkdir -p "${runtime_dir}"
  image="${runtime_dir}/ecoda-py-cuda13.sif"
  manifest="${image}.manifest"
  identity="${runtime_dir}/runtime.identity"
  printf 'runtime image\n' > "${image}"
  image_sha="$(sha256_file "${image}")"
  toml_sha="$(sha256_file "${SOURCE_TREE}/pixi.toml")"
  lock_sha="$(sha256_file "${SOURCE_TREE}/pixi.lock")"
  cat > "${manifest}" <<EOF
FORMAT=2
IMAGE_BUILD_GIT_REVISION=fixture
IMAGE_PATH=${image}
IMAGE_SHA256=${image_sha}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.49.0
APPTAINER_VERSION=1.3.2
IMAGE_PIXI_TOML_SHA256=${toml_sha}
IMAGE_PIXI_LOCK_SHA256=${lock_sha}
EOF
  manifest_sha="$(sha256_file "${manifest}")"
  image_size="$(wc -c < "${image}" | tr -d '[:space:]')"
  manifest_size="$(wc -c < "${manifest}" | tr -d '[:space:]')"
  cat > "${identity}" <<EOF
RUNTIME_IMAGE=${image}
RUNTIME_MANIFEST=${manifest}
RUNTIME_IMAGE_SHA256=${image_sha}
RUNTIME_MANIFEST_SHA256=${manifest_sha}
RUNTIME_IMAGE_SIZE=${image_size}
RUNTIME_MANIFEST_SIZE=${manifest_size}
IMAGE_PIXI_TOML_SHA256=${toml_sha}
IMAGE_PIXI_LOCK_SHA256=${lock_sha}
EOF
  chmod -R a-w "${runtime_dir}"
  RUNTIME_IMAGE="${image}"
  RUNTIME_MANIFEST="${manifest}"
  RUNTIME_IDENTITY="${identity}"
  export RUNTIME_IMAGE RUNTIME_MANIFEST RUNTIME_IDENTITY
}

make_validator_stub() {
  VALIDATOR_LOG="${TMP_DIR}/validator.log"
  VALIDATOR_BIN="${TMP_DIR}/validator-python"
  cat > "${VALIDATOR_BIN}" <<'EOF'
#!/bin/bash
set -euo pipefail
printf 'validator:%s\n' "$*" >> "${AUDIT_VALIDATOR_LOG}"
if [[ "${1:-}" == - ]]; then
  exec "${PIXI_PYTHON}" "$@"
fi
exit 0
EOF
  chmod +x "${VALIDATOR_BIN}"
  export AUDIT_VALIDATOR_LOG="${VALIDATOR_LOG}" PYTHON_BIN="${VALIDATOR_BIN}"
}

prepare_artifact_producer_run() {
  local run_root="${SCRATCH}/_ecoda_runs/audit-producer"
  mkdir -p "${run_root}/manifests/artifacts" "${run_root}/status"
  make_run_identities "${run_root}" "${SOURCE_MANIFEST}" "${RUNTIME_IDENTITY}"
  cat > "${run_root}/metadata" <<EOF
STAGE=stage5
SOURCE_MANIFEST=${run_root}/manifests/source.manifest
RUNTIME_IDENTITY=${run_root}/manifests/runtime.identity
RUN_ID=audit-producer
STATE=ACTIVE
PID=$$
CREATED=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
  cat > "${run_root}/status/terminal" <<EOF
STATE=OK
RUN_ID=audit-producer
REASON=fixture producer context
TIME=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
}


make_source_fixture
make_runtime_fixture
make_validator_stub

prepare_artifact_producer_run
prepare_contract_identity() {
  local run_root="$1" file_method="$2" model_id="$3"
  local identity_method="${4:-${file_method}}" path md5 size identity
  path="${run_root}/manifests/batch_contracts/Fixture__batch_effect_corrected__${file_method}.json"
  mkdir -p "$(dirname "${path}")"
  identity="$(
    cd "${ROOT}"
    ECODA_SOURCE_ROOT="${SOURCE_TREE}" PROJECT_ROOT="${ROOT}" \
      DATASETS_JSON_FILE="${SOURCE_CONFIG}" \
      BENCHMARK_MATRIX_TEST=0 PYTHON_BIN="${PIXI_PYTHON}" \
      ecoda_batch_contract_identity "${SOURCE_CONFIG}" Fixture \
      batch_effect_corrected "${identity_method}" "${model_id}"
  )"
  printf '%s\n' "${identity}" > "${path}"
  ecoda_write_checksum "${path}" >/dev/null
  ecoda_validate_checksum "${path}" >/dev/null
  md5="${ECODA_CHECKSUM_MD5}"
  size="${ECODA_CHECKSUM_SIZE}"
  CONTRACT_IDENTITY_PATH="${path}"
  CONTRACT_IDENTITY_MD5="${md5}"
  CONTRACT_IDENTITY_SIZE="${size}"
}

prepare_corrected_final_run() {
  local run_root="${SCRATCH}/_ecoda_runs/corrected_final_fixture"
  local analysis_root="${SCRATCH}/batch_effect/corrected_final"
  local analysis_nas_root="${NAS}/batch_effect/corrected_final"
  local input="${SCRATCH}/Fixture/output/Fixture_batch_effect_analysis_corrected.h5ad"
  local metadata_output="${analysis_root}/metadata/Fixture_sample_metadata.feather"
  local selection pending contract_manifest contract_md5 contract_size contract_sha
  local preprocess_path preprocess_md5 preprocess_size mrvi_path mrvi_md5 mrvi_size
  local mrvi_artifact mrvi_nas_artifact status_report
  local pending_md5 pending_size
  local preprocess_identity mrvi_identity
  mkdir -p "${run_root}/manifests" "${run_root}/status" "${run_root}/status/metadata_export"
  make_run_identities "${run_root}" "${SOURCE_MANIFEST}" "${RUNTIME_IDENTITY}"
  selection="${run_root}/manifests/selection.tsv"
  printf 'Fixture\tbatch_effect_corrected\tbatch_effect_corrected\n' > "${selection}"
  ecoda_write_checksum "${selection}" >/dev/null
  pending="${run_root}/manifests/pending_selection.tsv"
  printf 'Fixture\tbatch_effect_corrected\tmrvi\n' > "${pending}"
  ecoda_write_checksum "${pending}" >/dev/null
  pending_md5="${ECODA_CHECKSUM_MD5}"
  pending_size="${ECODA_CHECKSUM_SIZE}"

  make_artifact "${input}" 0
  mkdir -p "$(dirname "${metadata_output}")"
  printf 'sample metadata fixture\n' > "${metadata_output}"
  ecoda_write_checksum "${metadata_output}" >/dev/null

  preprocess_identity="$(
    cd "${ROOT}"
    ECODA_SOURCE_ROOT="${SOURCE_TREE}" PROJECT_ROOT="${ROOT}" \
      DATASETS_JSON_FILE="${SOURCE_CONFIG}" \
      BENCHMARK_MATRIX_TEST=0 PYTHON_BIN="${PIXI_PYTHON}" \
      ecoda_batch_contract_identity "${SOURCE_CONFIG}" Fixture \
      batch_effect_corrected preprocess hvg_composite_v1
  )"
  prepare_contract_identity "${run_root}" preprocess hvg_composite_v1
  preprocess_path="${CONTRACT_IDENTITY_PATH}"
  preprocess_md5="${CONTRACT_IDENTITY_MD5}"
  preprocess_size="${CONTRACT_IDENTITY_SIZE}"
  mrvi_identity="$(
    cd "${ROOT}"
    ECODA_SOURCE_ROOT="${SOURCE_TREE}" PROJECT_ROOT="${ROOT}" \
      DATASETS_JSON_FILE="${SOURCE_CONFIG}" \
      BENCHMARK_MATRIX_TEST=0 PYTHON_BIN="${PIXI_PYTHON}" \
      ecoda_batch_contract_identity "${SOURCE_CONFIG}" Fixture \
      batch_effect_corrected MrVI mrvi_composite_v1
  )"
  prepare_contract_identity "${run_root}" mrvi mrvi_composite_v1 MrVI
  mrvi_path="${CONTRACT_IDENTITY_PATH}"
  mrvi_md5="${CONTRACT_IDENTITY_MD5}"
  mrvi_size="${CONTRACT_IDENTITY_SIZE}"
  # Keep the computed identities in the files expected by the audit.  The
  # helper above derives their safe path from the method id.
  printf '%s\n' "${preprocess_identity}" > "${preprocess_path}"
  ecoda_write_checksum "${preprocess_path}" >/dev/null
  ecoda_validate_checksum "${preprocess_path}" >/dev/null
  preprocess_md5="${ECODA_CHECKSUM_MD5}"
  preprocess_size="${ECODA_CHECKSUM_SIZE}"
  printf '%s\n' "${mrvi_identity}" > "${mrvi_path}"
  ecoda_write_checksum "${mrvi_path}" >/dev/null
  ecoda_validate_checksum "${mrvi_path}" >/dev/null
  mrvi_md5="${ECODA_CHECKSUM_MD5}"
  mrvi_size="${ECODA_CHECKSUM_SIZE}"
  contract_manifest="${run_root}/manifests/batch_contract.tsv"
  printf 'Fixture\tbatch_effect_corrected\tpreprocess\t%s\t%s\t%s\n' \
    "${preprocess_path}" "${preprocess_md5}" "${preprocess_size}" > "${contract_manifest}"
  printf 'Fixture\tbatch_effect_corrected\tmrvi\t%s\t%s\t%s\n' \
    "${mrvi_path}" "${mrvi_md5}" "${mrvi_size}" >> "${contract_manifest}"
  ecoda_write_checksum "${contract_manifest}" >/dev/null
  ecoda_validate_checksum "${contract_manifest}" >/dev/null
  contract_md5="${ECODA_CHECKSUM_MD5}"
  contract_size="${ECODA_CHECKSUM_SIZE}"
  contract_sha="$(sha256_file "${contract_manifest}")"

  mrvi_artifact="${analysis_root}/embeddings/Fixture_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather"
  mrvi_nas_artifact="${analysis_nas_root}/embeddings/Fixture_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather"
  make_artifact "${mrvi_artifact}" 1
  make_artifact "${mrvi_nas_artifact}" 0
  printf 'Fixture\tbatch_effect_corrected\t%s\t%s\n' \
    "${input}" "${metadata_output}" > "${run_root}/manifests/metadata_export.tsv"
  ecoda_write_checksum "${run_root}/manifests/metadata_export.tsv" >/dev/null
  status_report="${run_root}/status/metadata_export.report"
  cat > "${status_report}" <<EOF
STATE=OK
ANALYSIS_VARIANT=corrected_final
ANALYSIS_ROOT=${analysis_root}
ANALYSIS_NAS_ROOT=${analysis_nas_root}
ANALYSIS_PASS=corrected
ANALYSIS_LOG_PREFIX=execution_times_batch_effect_corrected_final_
RUN_ID=corrected_final_fixture
MANIFEST=${run_root}/manifests/metadata_export.tsv
COUNT=1
PENDING=0
EOF
  ecoda_write_checksum "${status_report}" >/dev/null
  cat > "${run_root}/metadata" <<EOF
STAGE=stage5
RUN_ID=corrected_final_fixture
STATE=ACTIVE
SOURCE_MANIFEST=${run_root}/manifests/source.manifest
RUNTIME_IDENTITY=${run_root}/manifests/runtime.identity
METHODS=mrvi
PASS=corrected
ROOT=${analysis_root}
ANALYSIS_VARIANT=corrected_final
ANALYSIS_ROOT=${analysis_root}
ANALYSIS_NAS_ROOT=${analysis_nas_root}
ANALYSIS_PASS=corrected
ANALYSIS_LOG_PREFIX=execution_times_batch_effect_corrected_final_
METADATA_EXPORT_MANIFEST=${run_root}/manifests/metadata_export.tsv
METADATA_EXPORT_STATUS=${status_report}
PENDING_SELECTION=${pending}
PENDING_SELECTION_MD5=${pending_md5}
PENDING_SELECTION_SIZE=${pending_size}
BATCH_CONTRACT_MANIFEST=${contract_manifest}
BATCH_CONTRACT_MANIFEST_MD5=${contract_md5}
BATCH_CONTRACT_MANIFEST_SIZE=${contract_size}
BATCH_CONTRACT_MANIFEST_SHA256=${contract_sha}
EOF
  cat > "${run_root}/status/terminal" <<EOF
STATE=OK
RUN_ID=corrected_final_fixture
EOF
  : > "${run_root}/manifests/scheduler_ids.tsv"
  CORRECTED_RUN_ROOT="${run_root}"
  CORRECTED_ROOT="${analysis_root}"
  CORRECTED_NAS_ROOT="${analysis_nas_root}"
  CORRECTED_METADATA_MANIFEST="${run_root}/manifests/metadata_export.tsv"
  export CORRECTED_RUN_ROOT CORRECTED_ROOT CORRECTED_NAS_ROOT CORRECTED_METADATA_MANIFEST
}

run_audit() {
  local run_root="$1" stage="$2" selection="$3" expected_rc="$4"
  local rc audit_log
  audit_log="${AUDIT_LOG_DIR}/audit-${run_root##*/}.log"
  : > "${audit_log}"
  set +e
  (
    cd "${ROOT}"
    bash "${ROOT}/src/utils/bash/ecoda_run_audit.sh" \
      --run-root "${run_root}" --stage "${stage}" --selection "${selection}" \
      --source-manifest "${SOURCE_MANIFEST}" \
      --runtime-identity "${RUNTIME_IDENTITY}"
  ) >"${audit_log}" 2>&1
  rc=$?
  set -e
  [[ ${rc} -eq ${expected_rc} ]] || {
    echo "unexpected run audit status: expected ${expected_rc}, got ${rc}" >&2
    cat "${audit_log}" >&2
    return 1
  }
}

prepare_corrected_final_run
CORRECTED_SELECTION="${CORRECTED_RUN_ROOT}/manifests/selection.tsv"
run_audit "${CORRECTED_RUN_ROOT}" stage5 "${CORRECTED_SELECTION}" 0

grep -F -- "${CORRECTED_ROOT}/embeddings/Fixture_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather" \
  "${VALIDATOR_LOG}" >/dev/null
grep -F -- '--analysis-variant corrected_final' "${VALIDATOR_LOG}" >/dev/null
grep -F -- "${CORRECTED_ROOT}/metadata/Fixture_sample_metadata.feather" \
  "${CORRECTED_METADATA_MANIFEST}" >/dev/null
! grep -F -- "${SCRATCH}/batch_effect/corrected/" "${VALIDATOR_LOG}" >/dev/null
! grep -F -- "${SCRATCH}/batch_effect/uncorrected/" "${CORRECTED_METADATA_MANIFEST}" >/dev/null

# A corrected-final run cannot be made to consume a legacy root by changing
# metadata after the producer finished: the run-scoped identity must fail
# before selected artifact expansion.
sed "s#^ANALYSIS_ROOT=.*#ANALYSIS_ROOT=${SCRATCH}/batch_effect/corrected#" \
  "${CORRECTED_RUN_ROOT}/metadata" > "${CORRECTED_RUN_ROOT}/metadata.bad"
mv "${CORRECTED_RUN_ROOT}/metadata.bad" "${CORRECTED_RUN_ROOT}/metadata"
run_audit "${CORRECTED_RUN_ROOT}" stage5 "${CORRECTED_SELECTION}" 1

# Restore the accepted corrected-final identity and prove that the metadata
# exporter path is independently bound to the same variant root.
sed "s#^ANALYSIS_ROOT=.*#ANALYSIS_ROOT=${CORRECTED_ROOT}#" \
  "${CORRECTED_RUN_ROOT}/metadata" > "${CORRECTED_RUN_ROOT}/metadata.bad"
mv "${CORRECTED_RUN_ROOT}/metadata.bad" "${CORRECTED_RUN_ROOT}/metadata"
sed "s#${CORRECTED_ROOT}/metadata/#${SCRATCH}/batch_effect/corrected/metadata/#" \
  "${CORRECTED_METADATA_MANIFEST}" > "${CORRECTED_METADATA_MANIFEST}.bad"
mv "${CORRECTED_METADATA_MANIFEST}.bad" "${CORRECTED_METADATA_MANIFEST}"
ecoda_write_checksum "${CORRECTED_METADATA_MANIFEST}" >/dev/null
run_audit "${CORRECTED_RUN_ROOT}" stage5 "${CORRECTED_SELECTION}" 1

# No ANALYSIS_VARIANT metadata retains the legacy pass root and unqualified
# stem; this run has no final metadata-export contract.
LEGACY_RUN_ROOT="${SCRATCH}/_ecoda_runs/legacy_fixture"
LEGACY_ROOT="${SCRATCH}/batch_effect/uncorrected"
LEGACY_NAS_ROOT="${NAS}/batch_effect/uncorrected"
mkdir -p "${LEGACY_RUN_ROOT}/manifests" "${LEGACY_RUN_ROOT}/status"
make_run_identities "${LEGACY_RUN_ROOT}" "${SOURCE_MANIFEST}" "${RUNTIME_IDENTITY}"
LEGACY_SELECTION="${LEGACY_RUN_ROOT}/manifests/selection.tsv"
printf 'Fixture\tbatch_effect_uncorrected\tmrvi\n' > "${LEGACY_SELECTION}"
ecoda_write_checksum "${LEGACY_SELECTION}" >/dev/null
LEGACY_ARTIFACT="${LEGACY_ROOT}/embeddings/Fixture_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"
LEGACY_NAS_ARTIFACT="${LEGACY_NAS_ROOT}/embeddings/Fixture_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"
make_artifact "${LEGACY_ARTIFACT}" 1
make_artifact "${LEGACY_NAS_ARTIFACT}" 0
cat > "${LEGACY_RUN_ROOT}/metadata" <<EOF
STAGE=stage5
RUN_ID=legacy_fixture
STATE=ACTIVE
SOURCE_MANIFEST=${LEGACY_RUN_ROOT}/manifests/source.manifest
RUNTIME_IDENTITY=${LEGACY_RUN_ROOT}/manifests/runtime.identity
METHODS=mrvi
PASS=uncorrected
ROOT=${LEGACY_ROOT}
EOF
cat > "${LEGACY_RUN_ROOT}/status/terminal" <<EOF
STATE=OK
RUN_ID=legacy_fixture
EOF
: > "${LEGACY_RUN_ROOT}/manifests/scheduler_ids.tsv"
run_audit "${LEGACY_RUN_ROOT}" stage5 "${LEGACY_SELECTION}" 0
grep -F -- "${LEGACY_ARTIFACT}" "${VALIDATOR_LOG}" >/dev/null
! grep -F -- '_batch_effect_uncorrected_final_' "${VALIDATOR_LOG}" >/dev/null

echo "ecoda run audit variant identity: OK"
