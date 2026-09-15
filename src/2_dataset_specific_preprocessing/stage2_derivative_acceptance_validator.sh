#!/usr/bin/env bash
# Validator-only acceptance for the completed Alzheimer donor-by-assay Stage 2
# derivative.  This command consumes one explicit, already-published artifact;
# it never submits work, changes the producer run, or publishes an artifact
# record.  The only files it creates are evidence files below a fresh run root.
set -euo pipefail

SCRIPT_NAME="$(basename "$0")"

SOURCE_RUN_ID=""
SOURCE_RUN_ROOT=""
SOURCE_ROOT=""
SOURCE_MANIFEST_ARG=""
RUNTIME_IDENTITY_ARG=""
SOURCE_TERMINAL_ARG=""
SOURCE_OWNER_ARG=""
ARTIFACT_RECORD_ARG=""
SIDECAR_ARG=""
STEPS_MANIFEST_ARG=""
OWNERSHIP_MANIFEST_ARG=""
JOBS_MANIFEST_ARG=""
SCHEDULER_MANIFEST_ARG=""
WATCHDOG_STATUS_ARG=""
ARTIFACT_ARG=""
NEW_RUN_ID=""
NEW_RUN_ROOT=""
CONFIG_ARG=""
PYTHON_BIN_ARG=""
VALIDATOR_SCRIPT_ARG=""
SCRATCH_ROOT_ARG=""
ARRAY_ID=""
WATCHDOG_ID=""
PRIOR_EVIDENCE_ARG=""
EXPECTED_CELLS=1395601
EXPECTED_SAMPLES=104
EXPECTED_DONORS=83
EXPECTED_ASSAY_COUNTS="10x3v3=83,10xmultiome=21"
EXPECTED_SEX_COUNTS="female=59,male=45"
REQUIRE_EXAMPLES=1

SOURCE_MANIFEST_PATH=""
RUN_SOURCE_MANIFEST_PATH=""
RUNTIME_IDENTITY_PATH=""
RUN_RUNTIME_IDENTITY_PATH=""
SOURCE_METADATA_PATH=""
SOURCE_TERMINAL_PATH=""
SOURCE_OWNER_PATH=""
SOURCE_OWNER_DIR=""
STEPS_MANIFEST_PATH=""
OWNERSHIP_MANIFEST_PATH=""
JOBS_MANIFEST_PATH=""
SCHEDULER_MANIFEST_PATH=""
WATCHDOG_STATUS_PATH=""
ARTIFACT_RECORD_PATH=""
ARTIFACT_SIDECAR_PATH=""
ARTIFACT_PATH=""
RAW_INPUT_PATH=""
CONFIG_PATH=""
PYTHON_BIN=""
VALIDATOR_PATH=""
SCRATCH_ROOT=""
ARTIFACT_MD5=""
ARTIFACT_SIZE=""
OBSERVED_CELLS=""
OBSERVED_SAMPLES=""
SOURCE_MANIFEST_SHA256=""
SOURCE_COMMIT=""
SOURCE_SNAPSHOT_ROOT=""
RUNTIME_IMAGE=""
RUNTIME_MANIFEST=""
RUNTIME_IMAGE_SHA256=""
RUNTIME_MANIFEST_SHA256=""
RUNTIME_IMAGE_SIZE=""
RUNTIME_MANIFEST_SIZE=""
RUNTIME_FORMAT=""
PYTHON_SHA256=""
RAW_INPUT_NAME="SEAAD_Alzheimer.h5ad"
DERIVATIVE_NAME="SEAAD_Alzheimer_donor_assay.h5ad"

usage() {
  cat <<EOF
Usage: ${SCRIPT_NAME} \
  --source-run-id ID --source-run-root PATH \
  --artifact-path PATH --source-terminal PATH --source-owner PATH \
  --artifact-record PATH [--sidecar PATH] \
  --source-root PATH [--source-manifest PATH] [--runtime-identity PATH] \
  --config PATH --python-bin PATH [--validator-script PATH] \
  --run-id ID --run-root PATH \
  [--scratch-root PATH] --scheduler-array-id 4407671 \
  --scheduler-watchdog-id 4407672

Validate one existing Alzheimer donor-by-assay Stage 2 derivative and seal
only source evidence plus validator results into a fresh acceptance run root.
The accepted scheduler IDs are evidence bindings for the old Yggdrasil run;
no scheduler query or submission is performed.
EOF
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

safe_value() {
  local value="${1:-}"
  [[ -n "${value}" && "${value}" != *$'\n'* &&
     "${value}" != *$'\r'* && "${value}" != *$'\t'* &&
     "${value}" != *'='* ]] || return 1
}

require_absolute() {
  local value="${1:-}" label="${2:-path}"
  safe_value "${value}" || die "${label} is empty or contains a record delimiter"
  [[ "${value}" = /* ]] || die "${label} must be absolute: ${value}"
  case "${value}" in
    */../*|*/..|../*|..) die "${label} contains a parent-directory escape: ${value}" ;;
  esac
}

require_regular_file() {
  local path="$1" label="$2"
  require_absolute "${path}" "${label}"
  [[ -f "${path}" && ! -L "${path}" && -s "${path}" ]] ||
    die "${label} must be a non-empty regular file: ${path}"
}

require_regular_dir() {
  local path="$1" label="$2"
  require_absolute "${path}" "${label}"
  [[ -d "${path}" && ! -L "${path}" ]] ||
    die "${label} must be a regular directory: ${path}"
}

canonical_existing_file() {
  local path="$1" label="$2" resolved
  require_regular_file "${path}" "${label}"
  command -v realpath >/dev/null 2>&1 || die "realpath is required for ${label}"
  resolved="$(realpath "${path}" 2>/dev/null || true)"
  [[ -n "${resolved}" && "${resolved}" == "${path}" ]] ||
    die "${label} is not canonical: ${path}"
  printf '%s\n' "${resolved}"
}

canonical_existing_dir() {
  local path="$1" label="$2" resolved
  require_regular_dir "${path}" "${label}"
  command -v realpath >/dev/null 2>&1 || die "realpath is required for ${label}"
  resolved="$(realpath "${path}" 2>/dev/null || true)"
  [[ -n "${resolved}" && "${resolved}" == "${path}" ]] ||
    die "${label} is not canonical: ${path}"
  printf '%s\n' "${resolved}"
}

file_size() {
  wc -c < "$1" | tr -d '[:space:]'
}

sha256_file() {
  local path="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${path}" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${path}" | awk '{print $1}'
  else
    die "sha256sum or shasum is required"
  fi
}

md5_file() {
  local path="$1"
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "${path}" | awk '{print $1}'
  elif command -v md5 >/dev/null 2>&1; then
    md5 -q "${path}"
  else
    die "md5sum or md5 is required to checksum acceptance evidence"
  fi
}

mode_is_nonwritable() {
  local path="$1" mode
  mode="$(stat -c '%a' "${path}" 2>/dev/null || stat -f '%Lp' "${path}" 2>/dev/null || true)"
  [[ "${mode}" =~ ^[0-7]{3,4}$ ]] || return 1
  [[ "${mode: -3}" != *[2367]* ]]
}

require_nonwritable() {
  local path="$1" label="$2"
  mode_is_nonwritable "${path}" || die "${label} is writable or its mode is unavailable: ${path}"
}

path_within() {
  local candidate="$1" root="$2"
  [[ "${candidate}" == "${root}" || "${candidate}" == "${root}/"* ]]
}

field_value() {
  local path="$1" key="$2" value
  value="$(awk -v wanted="${key}" '
    index($0, wanted "=") == 1 {
      count++
      value = substr($0, length(wanted) + 2)
    }
    END {
      if (count != 1 || value == "") exit 1
      print value
    }
  ' "${path}")" || return 1
  safe_value "${value}" || return 1
  printf '%s\n' "${value}"
}

json_quote() {
  local value="${1:-}"
  value="${value//\\/\\\\}"
  value="${value//\"/\\\"}"
  value="${value//$'\n'/\\n}"
  value="${value//$'\r'/\\r}"
  value="${value//$'\t'/\\t}"
  printf '"%s"' "${value}"
}

write_atomic() {
  local destination="$1" content="$2" parent temporary
  require_absolute "${destination}" "atomic destination"
  parent="$(dirname "${destination}")"
  [[ -d "${parent}" && ! -L "${parent}" ]] ||
    die "atomic destination parent is missing or symlinked: ${parent}"
  [[ ! -L "${destination}" ]] || die "atomic destination is a symlink: ${destination}"
  temporary="${destination}.tmp.$$"
  [[ ! -e "${temporary}" && ! -L "${temporary}" ]] ||
    die "atomic temporary destination already exists: ${temporary}"
  umask 077
  printf '%b' "${content}" > "${temporary}" || {
    rm -f "${temporary}"
    die "could not write atomic temporary file: ${destination}"
  }
  mv -f "${temporary}" "${destination}" || {
    rm -f "${temporary}"
    die "could not install atomic file: ${destination}"
  }
}

write_checked_file() {
  local destination="$1" content="$2" digest size
  write_atomic "${destination}" "${content}"
  digest="$(md5_file "${destination}")" || die "could not checksum acceptance evidence: ${destination}"
  size="$(file_size "${destination}")"
  [[ "${digest}" =~ ^[[:xdigit:]]{32}$ && "${size}" =~ ^[1-9][0-9]*$ ]] ||
    die "acceptance evidence checksum construction failed: ${destination}"
  write_atomic "${destination}.md5" \
    "MD5=${digest}\nSIZE=${size}\nPATH=${destination}\n"
  chmod a-w "${destination}" "${destination}.md5" ||
    die "could not seal acceptance evidence: ${destination}"
}

copy_sealed() {
  local source="$1" destination="$2"
  require_regular_file "${source}" "source evidence"
  [[ ! -L "${destination}" ]] || die "evidence destination is a symlink: ${destination}"
  cp -p "${source}" "${destination}" || die "could not copy source evidence: ${source}"
  chmod a-w "${destination}" || die "could not seal copied evidence: ${destination}"
}

validate_manifest_file() {
  local path="$1" label="$2" expected_count="$3"
  local line key value index=0
  local -a seen=()
  require_regular_file "${path}" "${label}"
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    [[ "${line}" == *=* ]] || die "${label} has malformed line ${index}"
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" =~ ^[A-Z][A-Z0-9_]*$ ]] || die "${label} has an unsafe key"
    safe_value "${value}" || die "${label} has an unsafe value for ${key}"
    local prior
    for prior in "${seen[@]:-}"; do
      [[ "${prior}" != "${key}" ]] || die "${label} repeats ${key}"
    done
    seen+=("${key}")
  done < "${path}"
  if [[ "${expected_count}" != 0 ]]; then
    [[ "${index}" == "${expected_count}" ]] ||
      die "${label} has ${index} fields; expected ${expected_count}"
  else
    [[ "${index}" -gt 0 ]] || die "${label} is empty"
  fi
}

validate_source_identity() {
  local line key value expected_key index=0
  local source_snapshot="" source_snapshot_sha256="" complete_marker actual expected
  local required_file source_file
  local -a keys=(FORMAT SOURCE_ROOT SOURCE_COMMIT SOURCE_ARCHIVE_PATH
    SOURCE_ARCHIVE_SHA256 CONFIG_HELPER_SHA256 DATASETS_SHA256
    PIXI_TOML_SHA256 PIXI_LOCK_SHA256 AUX_ROOT SCGATE_DB_BRANCH)
  SOURCE_MANIFEST_PATH="$(canonical_existing_file "${SOURCE_MANIFEST_ARG}" \
    "immutable snapshot source manifest")"
  RUN_SOURCE_MANIFEST_PATH="$(canonical_existing_file \
    "${SOURCE_RUN_ROOT}/manifests/source.manifest" \
    "source-run source manifest")"
  cmp -s "${RUN_SOURCE_MANIFEST_PATH}" "${SOURCE_MANIFEST_PATH}" ||
    die "source-run source manifest differs from the current immutable source identity"
  [[ "${SOURCE_MANIFEST_PATH}" == "${SOURCE_ROOT%/tree}/identity/source.manifest" ]] ||
    die "source manifest is not bound to the immutable source tree"
  validate_manifest_file "${SOURCE_MANIFEST_PATH}" \
    "immutable snapshot source manifest" 11
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    expected_key="${keys[$((index - 1))]}"
    [[ "${line}" == "${expected_key}="* ]] ||
      die "source identity field ${index} must be ${expected_key}"
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" ]] ||
      die "source identity field ${index} is malformed"
    safe_value "${value}" || die "source identity ${key} is unsafe"
    case "${key}" in
      FORMAT) [[ "${value}" == 1 ]] || die "source identity FORMAT must be 1" ;;
      SOURCE_ROOT) [[ "${value}" == "${SOURCE_ROOT}" ]] ||
        die "source identity SOURCE_ROOT mismatch" ;;
      SOURCE_COMMIT) SOURCE_COMMIT="${value}" ;;
      SOURCE_ARCHIVE_PATH) source_snapshot="${value}" ;;
      SOURCE_ARCHIVE_SHA256)
        source_snapshot_sha256="${value}"
        [[ "${value}" =~ ^[[:xdigit:]]{64}$ ]] ||
          die "source archive SHA-256 is malformed"
        ;;
      CONFIG_HELPER_SHA256|DATASETS_SHA256|PIXI_TOML_SHA256|PIXI_LOCK_SHA256)
        [[ "${value}" =~ ^[[:xdigit:]]{64}$ ]] ||
          die "source identity digest is malformed: ${key}"
        case "${key}" in
          CONFIG_HELPER_SHA256) SOURCE_CONFIG_HELPER_SHA256="${value}" ;;
          DATASETS_SHA256) SOURCE_DATASETS_SHA256="${value}" ;;
          PIXI_TOML_SHA256) SOURCE_PIXI_TOML_SHA256="${value}" ;;
          PIXI_LOCK_SHA256) SOURCE_PIXI_LOCK_SHA256="${value}" ;;
        esac
        ;;
      AUX_ROOT) [[ "${value}" == "${SOURCE_ROOT}/aux" ]] ||
        die "source identity AUX_ROOT mismatch" ;;
      SCGATE_DB_BRANCH) [[ -n "${value}" ]] ||
        die "source identity SCGATE_DB_BRANCH is empty" ;;
    esac
  done < "${SOURCE_MANIFEST_PATH}"
  [[ "${SOURCE_COMMIT}" =~ ^[[:xdigit:]]{40}$ ]] ||
    die "source identity commit is not a full commit"
  SOURCE_SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
  [[ "$(basename "${SOURCE_SNAPSHOT_ROOT}")" == "${SOURCE_COMMIT}" ]] ||
    die "source identity commit does not match the snapshot directory"
  complete_marker="${SOURCE_SNAPSHOT_ROOT}/COMPLETE"
  [[ "${source_snapshot}" == "${SOURCE_SNAPSHOT_ROOT}/identity/source.tar" ]] ||
    die "source identity archive is not bound to the snapshot identity directory"
  require_regular_file "${complete_marker}" "immutable snapshot COMPLETE marker"
  [[ "$(cat "${complete_marker}")" == COMPLETE ]] ||
    die "snapshot COMPLETE marker is invalid"
  require_nonwritable "${SOURCE_ROOT}" "immutable source tree"
  require_nonwritable "${SOURCE_MANIFEST_PATH}" "immutable snapshot source manifest"
  require_regular_file "${source_snapshot}" "immutable source archive"
  require_nonwritable "${source_snapshot}" "immutable source archive"
  actual="$(sha256_file "${source_snapshot}")"
  expected="$(printf '%s' "${source_snapshot_sha256}" | tr '[:upper:]' '[:lower:]')"
  [[ -n "${expected}" && "${actual}" == "${expected}" ]] ||
    die "immutable source archive digest does not match source identity"
  for required_file in \
    "${SOURCE_ROOT}/config_helper.R" "${SOURCE_ROOT}/datasets.json" \
    "${SOURCE_ROOT}/pixi.toml" "${SOURCE_ROOT}/pixi.lock" \
    "${SOURCE_ROOT}/aux/scGateDB.rds" \
    "${SOURCE_ROOT}/aux/genes.blocklist.rds" \
    "${SOURCE_ROOT}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"; do
    require_regular_file "${required_file}" "immutable snapshot source file"
    require_nonwritable "${required_file}" "immutable snapshot source file"
  done
  for source_file in \
    "${SOURCE_ROOT}/config_helper.R" "${SOURCE_ROOT}/datasets.json" \
    "${SOURCE_ROOT}/pixi.toml" "${SOURCE_ROOT}/pixi.lock"; do
    case "${source_file}" in
      */config_helper.R) expected="${SOURCE_CONFIG_HELPER_SHA256}" ;;
      */datasets.json) expected="${SOURCE_DATASETS_SHA256}" ;;
      */pixi.toml) expected="${SOURCE_PIXI_TOML_SHA256}" ;;
      */pixi.lock) expected="${SOURCE_PIXI_LOCK_SHA256}" ;;
    esac
    actual="$(sha256_file "${source_file}")"
    [[ "${actual}" == "$(printf '%s' "${expected}" | tr '[:upper:]' '[:lower:]')" ]] ||
      die "immutable snapshot source-file digest does not match identity: ${source_file}"
  done
  require_nonwritable "${complete_marker}" "immutable snapshot COMPLETE marker"
  SOURCE_MANIFEST_SHA256="$(sha256_file "${SOURCE_MANIFEST_PATH}")"
}

validate_runtime_identity() {
  local image_path image_sha manifest_format actual actual_image actual_manifest identity_count
  local -a base_keys=(RUNTIME_IMAGE RUNTIME_MANIFEST RUNTIME_IMAGE_SHA256
    RUNTIME_MANIFEST_SHA256 RUNTIME_IMAGE_SIZE RUNTIME_MANIFEST_SIZE)
  local -a dependency_keys=(IMAGE_PIXI_TOML_SHA256 IMAGE_PIXI_LOCK_SHA256)
  RUNTIME_IDENTITY_PATH="$(canonical_existing_file "${RUNTIME_IDENTITY_ARG}" \
    "current runtime identity manifest")"
  RUN_RUNTIME_IDENTITY_PATH="$(canonical_existing_file \
    "${SOURCE_RUN_ROOT}/manifests/runtime.identity" \
    "source-run runtime identity manifest")"
  require_nonwritable "${RUNTIME_IDENTITY_PATH}" "current runtime identity manifest"
  if [[ -n "${ECODA_RUNTIME_IDENTITY:-}" &&
        "${ECODA_RUNTIME_IDENTITY}" != "${RUNTIME_IDENTITY_PATH}" &&
        "${ECODA_RUNTIME_IDENTITY}" != "${RUN_RUNTIME_IDENTITY_PATH}" ]]; then
    die "exported runtime identity does not match the explicit current identity"
  fi
  cmp -s "${RUNTIME_IDENTITY_PATH}" "${RUN_RUNTIME_IDENTITY_PATH}" ||
    die "source-run runtime identity differs from the current runtime identity"
  identity_count="$(wc -l < "${RUNTIME_IDENTITY_PATH}" | tr -d '[:space:]')"
  [[ "${identity_count}" == 6 || "${identity_count}" == 8 ]] ||
    die "runtime identity has an unexpected field count"
  index=0
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    if [[ ${index} -le 6 ]]; then
      expected_key="${base_keys[$((index - 1))]}"
    else
      expected_key="${dependency_keys[$((index - 7))]}"
    fi
    [[ "${line}" == "${expected_key}="* ]] ||
      die "runtime identity field ${index} must be ${expected_key}"
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" ]] || die "runtime identity field ${index} is malformed"
    safe_value "${value}" || die "runtime identity ${key} is unsafe"
    case "${key}" in
      RUNTIME_IMAGE) RUNTIME_IMAGE="${value}" ;;
      RUNTIME_MANIFEST) RUNTIME_MANIFEST="${value}" ;;
      RUNTIME_IMAGE_SHA256) RUNTIME_IMAGE_SHA256="${value}" ;;
      RUNTIME_MANIFEST_SHA256) RUNTIME_MANIFEST_SHA256="${value}" ;;
      RUNTIME_IMAGE_SIZE) RUNTIME_IMAGE_SIZE="${value}" ;;
      RUNTIME_MANIFEST_SIZE) RUNTIME_MANIFEST_SIZE="${value}" ;;
      IMAGE_PIXI_TOML_SHA256) RUNTIME_IMAGE_PIXI_SHA256="${value}" ;;
      IMAGE_PIXI_LOCK_SHA256) RUNTIME_IMAGE_LOCK_SHA256="${value}" ;;
    esac
  done < "${RUNTIME_IDENTITY_PATH}"
  [[ "${RUNTIME_IMAGE}" = /* && "${RUNTIME_MANIFEST}" = /* ]] ||
    die "runtime identity paths must be absolute"
  [[ "${RUNTIME_IMAGE}" == */_ecoda_runtime/*/*.sif ]] ||
    die "runtime image is not a versioned immutable runtime image"
  [[ "${RUNTIME_MANIFEST}" == "${RUNTIME_IMAGE}.manifest" ]] ||
    die "runtime manifest is not bound beside the runtime image"
  [[ "${RUNTIME_IMAGE_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${RUNTIME_MANIFEST_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${RUNTIME_IMAGE_SIZE}" =~ ^[1-9][0-9]*$ &&
     "${RUNTIME_MANIFEST_SIZE}" =~ ^[1-9][0-9]*$ ]] ||
    die "runtime identity digest or size is malformed"
  RUNTIME_IMAGE="$(canonical_existing_file "${RUNTIME_IMAGE}" "runtime image")"
  RUNTIME_MANIFEST="$(canonical_existing_file "${RUNTIME_MANIFEST}" "runtime manifest")"
  require_nonwritable "${RUNTIME_IMAGE}" "runtime image"
  require_nonwritable "${RUNTIME_MANIFEST}" "runtime manifest"
  [[ "$(file_size "${RUNTIME_IMAGE}")" == "${RUNTIME_IMAGE_SIZE}" &&
     "$(file_size "${RUNTIME_MANIFEST}")" == "${RUNTIME_MANIFEST_SIZE}" ]] ||
    die "runtime identity size does not match the immutable runtime"
  actual_image="$(sha256_file "${RUNTIME_IMAGE}")"
  [[ "${actual_image}" == "$(printf '%s' "${RUNTIME_IMAGE_SHA256}" | tr '[:upper:]' '[:lower:]')" ]] ||
    die "runtime image digest does not match runtime identity"
  actual_manifest="$(sha256_file "${RUNTIME_MANIFEST}")"
  [[ "${actual_manifest}" == "$(printf '%s' "${RUNTIME_MANIFEST_SHA256}" | tr '[:upper:]' '[:lower:]')" ]] ||
    die "runtime manifest digest does not match runtime identity"
  image_path="$(field_value "${RUNTIME_MANIFEST}" IMAGE_PATH 2>/dev/null || true)"
  image_sha="$(field_value "${RUNTIME_MANIFEST}" IMAGE_SHA256 2>/dev/null || true)"
  [[ "${image_path}" == "${RUNTIME_IMAGE}" && "${image_sha}" == "${actual_image}" ]] ||
    die "runtime manifest image path/digest does not match runtime identity"
  validate_manifest_file "${RUNTIME_MANIFEST}" "runtime manifest" 0
  manifest_format="$(field_value "${RUNTIME_MANIFEST}" FORMAT)" || die "runtime manifest FORMAT is missing"
  RUNTIME_FORMAT="${manifest_format}"
  if [[ "${manifest_format}" == 2 ]]; then
    [[ "${identity_count}" == 8 ]] || die "FORMAT=2 runtime requires dependency identity fields"
    local required_key
    for required_key in IMAGE_BUILD_GIT_REVISION RUNTIME_ENV RUNTIME_LAYOUT \
      CONTAINER_ENV_PREFIX BASE_IMAGE PIXITAINER_VERSION PIXI_VERSION \
      APPTAINER_VERSION IMAGE_PIXI_TOML_SHA256 IMAGE_PIXI_LOCK_SHA256; do
      value="$(field_value "${RUNTIME_MANIFEST}" "${required_key}" 2>/dev/null || true)"
      [[ -n "${value}" ]] || die "runtime manifest is missing ${required_key}"
    done
    [[ "$(field_value "${RUNTIME_MANIFEST}" RUNTIME_ENV)" == py-cuda13 &&
       ( "$(field_value "${RUNTIME_MANIFEST}" RUNTIME_LAYOUT)" == relocated ||
         "$(field_value "${RUNTIME_MANIFEST}" RUNTIME_LAYOUT)" == path-preserving ) ]] ||
      die "runtime manifest toolchain layout is invalid"
    [[ "$(field_value "${RUNTIME_MANIFEST}" BASE_IMAGE)" == rockylinux:9 &&
       "$(field_value "${RUNTIME_MANIFEST}" PIXITAINER_VERSION)" == 0.8.3 ]] ||
      die "runtime manifest toolchain identity is invalid"
    [[ "${RUNTIME_IMAGE_PIXI_SHA256}" == "${SOURCE_PIXI_TOML_SHA256}" &&
       "${RUNTIME_IMAGE_LOCK_SHA256}" == "${SOURCE_PIXI_LOCK_SHA256}" ]] ||
      die "runtime dependency identity does not match the source snapshot"
  elif [[ "${manifest_format}" != 1 ]]; then
    die "runtime manifest FORMAT is unsupported"
  fi
}

validate_run_identity() {
  local metadata_state
  SOURCE_METADATA_PATH="$(canonical_existing_file "${SOURCE_RUN_ROOT}/metadata" \
    "source Stage 2 metadata")"
  [[ "$(field_value "${SOURCE_METADATA_PATH}" STAGE 2>/dev/null || true)" == stage2 ]] ||
    die "source run metadata is not Stage 2"
  [[ "$(field_value "${SOURCE_METADATA_PATH}" RUN_ID 2>/dev/null || true)" == "${SOURCE_RUN_ID}" ]] ||
    die "source run metadata RUN_ID mismatch"
  metadata_state="$(field_value "${SOURCE_METADATA_PATH}" STATE 2>/dev/null || true)"
  case "${metadata_state}" in ACTIVE|OK|NOOP_VALIDATED) ;; *)
    die "source run metadata has an invalid state: ${metadata_state}" ;;
  esac
  SOURCE_TERMINAL_PATH="$(canonical_existing_file "${SOURCE_TERMINAL_ARG}" \
    "source Stage 2 terminal status")"
  path_within "${SOURCE_TERMINAL_PATH}" "${SOURCE_RUN_ROOT}" ||
    die "source terminal status escapes the source run root"
  [[ "$(field_value "${SOURCE_TERMINAL_PATH}" STATE 2>/dev/null || true)" == OK ]] ||
    die "source Stage 2 terminal status is not STATE=OK"
  [[ "$(field_value "${SOURCE_TERMINAL_PATH}" RUN_ID 2>/dev/null || true)" == "${SOURCE_RUN_ID}" ]] ||
    die "source terminal RUN_ID mismatch"
}

validate_scheduler_binding() {
  local path="$1" line key value
  [[ -r "${path}" ]] || return 0
  while IFS= read -r line || [[ -n "${line}" ]]; do
    key="${line%%=*}"
    value="${line#*=}"
    case "${key}" in
      ARRAY_JOB_ID|SCHEDULER_ARRAY_ID)
        [[ "${value}" == "${ARRAY_ID}" ]] || die "scheduler array binding mismatch in ${path}" ;;
      WATCHDOG_JOB_ID|SCHEDULER_WATCHDOG_ID)
        [[ "${value}" == "${WATCHDOG_ID}" ]] || die "scheduler watchdog binding mismatch in ${path}" ;;
      SCHEDULER_ID)
        [[ "${value}" == "${ARRAY_ID}" || "${value}" == "${WATCHDOG_ID}" ]] ||
          die "source evidence contains an unbound scheduler ID: ${value}" ;;
    esac
  done < "${path}"
}

validate_sidecar() {
  local sidecar="$1" line1 line2 line3 md5 size recorded_path
  ARTIFACT_SIDECAR_PATH="$(canonical_existing_file "${sidecar}" "artifact checksum sidecar")"
  [[ "${ARTIFACT_SIDECAR_PATH}" == "${ARTIFACT_PATH}.md5" ]] ||
    die "artifact sidecar is not the canonical sibling sidecar"
  line1="$(sed -n '1p' "${ARTIFACT_SIDECAR_PATH}")"
  line2="$(sed -n '2p' "${ARTIFACT_SIDECAR_PATH}")"
  line3="$(sed -n '3p' "${ARTIFACT_SIDECAR_PATH}")"
  [[ "$(wc -l < "${ARTIFACT_SIDECAR_PATH}" | tr -d '[:space:]')" == 3 &&
     "${line1}" == MD5=* && "${line2}" == SIZE=* && "${line3}" == PATH=* ]] ||
    die "artifact checksum sidecar schema is invalid"
  md5="${line1#MD5=}"
  size="${line2#SIZE=}"
  recorded_path="${line3#PATH=}"
  [[ "${md5}" =~ ^[[:xdigit:]]{32}$ ]] || die "artifact sidecar MD5 is malformed"
  [[ "${size}" =~ ^[1-9][0-9]*$ ]] || die "artifact sidecar SIZE is malformed"
  [[ "${recorded_path}" == "${ARTIFACT_PATH}" ]] || die "artifact sidecar PATH mismatch"
  ARTIFACT_SIZE="$(file_size "${ARTIFACT_PATH}")"
  [[ "${size}" == "${ARTIFACT_SIZE}" ]] || die "artifact sidecar SIZE does not match artifact"
  ARTIFACT_MD5="$(printf '%s' "${md5}" | tr '[:upper:]' '[:lower:]')"
}

validate_artifact_record() {
  local line1 line2 line3 line4 line5 line6 path size md5 run producer state
  ARTIFACT_RECORD_PATH="$(canonical_existing_file "${ARTIFACT_RECORD_ARG}" \
    "source artifact record")"
  path_within "${ARTIFACT_RECORD_PATH}" "${SOURCE_RUN_ROOT}" ||
    die "source artifact record escapes the source run root"
  [[ "$(wc -l < "${ARTIFACT_RECORD_PATH}" | tr -d '[:space:]')" == 6 ]] ||
    die "source artifact record schema is invalid"
  line1="$(sed -n '1p' "${ARTIFACT_RECORD_PATH}")"
  line2="$(sed -n '2p' "${ARTIFACT_RECORD_PATH}")"
  line3="$(sed -n '3p' "${ARTIFACT_RECORD_PATH}")"
  line4="$(sed -n '4p' "${ARTIFACT_RECORD_PATH}")"
  line5="$(sed -n '5p' "${ARTIFACT_RECORD_PATH}")"
  line6="$(sed -n '6p' "${ARTIFACT_RECORD_PATH}")"
  [[ "${line1}" == PATH=* && "${line2}" == SIZE=* && "${line3}" == MD5=* &&
     "${line4}" == RUN_ID=* && "${line5}" == PRODUCER=* &&
     "${line6}" == STATE=* ]] || die "source artifact record fields are invalid"
  path="${line1#PATH=}"
  size="${line2#SIZE=}"
  md5="${line3#MD5=}"
  run="${line4#RUN_ID=}"
  producer="${line5#PRODUCER=}"
  state="${line6#STATE=}"
  [[ "${path}" == "${ARTIFACT_PATH}" ]] || die "source artifact record PATH mismatch"
  [[ "${size}" == "${ARTIFACT_SIZE}" ]] || die "source artifact record SIZE mismatch"
  [[ "$(printf '%s' "${md5}" | tr '[:upper:]' '[:lower:]')" == "${ARTIFACT_MD5}" ]] ||
    die "source artifact record MD5 does not match the immutable sidecar"
  [[ "${run}" == "${SOURCE_RUN_ID}" ]] || die "source artifact record RUN_ID mismatch"
  [[ "${producer}" == alzheimer_donor_assay ]] ||
    die "source artifact record producer is not alzheimer_donor_assay"
  [[ "${state}" == PUBLISHED ]] || die "source artifact record is not STATE=PUBLISHED"
}

validate_owner() {
  local owner_arg="$1" owner_run owner_stage owner_state owner_path owner_key
  if [[ -d "${owner_arg}" && ! -L "${owner_arg}" ]]; then
    SOURCE_OWNER_DIR="$(canonical_existing_dir "${owner_arg}" "source artifact owner directory")"
    SOURCE_OWNER_PATH="$(canonical_existing_file "${SOURCE_OWNER_DIR}/owner" \
      "source artifact owner record")"
  else
    SOURCE_OWNER_PATH="$(canonical_existing_file "${owner_arg}" "source artifact owner record")"
    SOURCE_OWNER_DIR="$(dirname "${SOURCE_OWNER_PATH}")"
  fi
  owner_run="$(field_value "${SOURCE_OWNER_PATH}" RUN_ID 2>/dev/null || true)"
  owner_stage="$(field_value "${SOURCE_OWNER_PATH}" STAGE 2>/dev/null || true)"
  owner_state="$(field_value "${SOURCE_OWNER_PATH}" STATE 2>/dev/null || true)"
  [[ "${owner_run}" == "${SOURCE_RUN_ID}" ]] || die "source owner RUN_ID mismatch"
  [[ "${owner_stage}" == stage2 ]] || die "source owner STAGE is not stage2"
  [[ "${owner_state}" == OK ]] || die "source owner is not terminal STATE=OK"
  owner_path="$(field_value "${SOURCE_OWNER_PATH}" PATH 2>/dev/null || true)"
  owner_key="$(field_value "${SOURCE_OWNER_PATH}" KEY 2>/dev/null || true)"
  if [[ -n "${owner_path}" ]]; then
    [[ "${owner_path}" == "${ARTIFACT_PATH}" ]] || die "source owner PATH mismatch"
  elif [[ "${owner_key}" == alzheimer_donor_assay ]]; then
    :
  else
    die "source owner is not bound to the Alzheimer donor-assay step"
  fi
}
validate_stage2_manifests() {
  local line step script outputs dependency owner extra
  local job_step job_id scheduler_kind scheduler_id
  local expected_script expected_owner
  local steps_path ownership_path jobs_path scheduler_path watchdog_path
  local scheduler_count array_count watchdog_count
  steps_path="${STEPS_MANIFEST_ARG:-${SOURCE_RUN_ROOT}/manifests/steps.tsv}"
  ownership_path="${OWNERSHIP_MANIFEST_ARG:-${SOURCE_RUN_ROOT}/manifests/ownership.tsv}"
  jobs_path="${JOBS_MANIFEST_ARG:-${SOURCE_RUN_ROOT}/manifests/jobs.tsv}"
  scheduler_path="${SCHEDULER_MANIFEST_ARG:-${SOURCE_RUN_ROOT}/manifests/scheduler_ids.tsv}"
  watchdog_path="${WATCHDOG_STATUS_ARG:-${SOURCE_RUN_ROOT}/status/watchdog}"
  STEPS_MANIFEST_PATH="$(canonical_existing_file "${steps_path}" \
    "source Stage 2 steps manifest")"
  OWNERSHIP_MANIFEST_PATH="$(canonical_existing_file "${ownership_path}" \
    "source Stage 2 ownership manifest")"
  JOBS_MANIFEST_PATH="$(canonical_existing_file "${jobs_path}" \
    "source Stage 2 jobs manifest")"
  SCHEDULER_MANIFEST_PATH="$(canonical_existing_file "${scheduler_path}" \
    "source Stage 2 scheduler manifest")"
  WATCHDOG_STATUS_PATH="$(canonical_existing_file "${watchdog_path}" \
    "source Stage 2 watchdog status")"
  path_within "${STEPS_MANIFEST_PATH}" "${SOURCE_RUN_ROOT}" ||
    die "source Stage 2 steps manifest escapes the source run root"
  path_within "${OWNERSHIP_MANIFEST_PATH}" "${SOURCE_RUN_ROOT}" ||
    die "source Stage 2 ownership manifest escapes the source run root"
  path_within "${JOBS_MANIFEST_PATH}" "${SOURCE_RUN_ROOT}" ||
    die "source Stage 2 jobs manifest escapes the source run root"
  path_within "${SCHEDULER_MANIFEST_PATH}" "${SOURCE_RUN_ROOT}" ||
    die "source Stage 2 scheduler manifest escapes the source run root"
  path_within "${WATCHDOG_STATUS_PATH}" "${SOURCE_RUN_ROOT}" ||
    die "source Stage 2 watchdog status escapes the source run root"
  [[ "$(wc -l < "${STEPS_MANIFEST_PATH}" | tr -d '[:space:]')" == 1 ]] ||
    die "source Stage 2 steps manifest must contain exactly one row"
  IFS=$'\t' read -r step script outputs dependency owner extra \
    < "${STEPS_MANIFEST_PATH}"
  [[ "${step}" == alzheimer_donor_assay && -n "${script}" &&
     -n "${outputs}" && "${dependency}" == "-" && -n "${owner}" &&
     -z "${extra}" ]] || die "source Stage 2 steps row is not the exact donor-assay row"
  expected_script="${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.7_submit_alzheimer_donor_assay.sh"
  expected_script="$(canonical_existing_file "${expected_script}" \
    "immutable Alzheimer Stage 2 hook")"
  [[ "${script}" == "${expected_script}" ]] ||
    die "source Stage 2 script is not the immutable Alzheimer hook"
  [[ "${outputs}" == "${ARTIFACT_PATH}" ]] ||
    die "source Stage 2 output is not the selected derivative artifact"
  expected_owner="${SOURCE_OWNER_DIR}"
  [[ "${owner}" == "${expected_owner}" ]] ||
    die "source Stage 2 owner differs from the selected owner evidence"
  [[ "$(wc -l < "${OWNERSHIP_MANIFEST_PATH}" | tr -d '[:space:]')" == 1 ]] ||
    die "source Stage 2 ownership manifest must contain exactly one row"
  cmp -s "${STEPS_MANIFEST_PATH}" "${OWNERSHIP_MANIFEST_PATH}" ||
    die "source Stage 2 ownership manifest does not match the one-row steps manifest"
  [[ "$(wc -l < "${JOBS_MANIFEST_PATH}" | tr -d '[:space:]')" == 1 ]] ||
    die "source Stage 2 jobs manifest must contain exactly one row"
  IFS=$'\t' read -r job_step job_id extra < "${JOBS_MANIFEST_PATH}"
  [[ "${job_step}" == alzheimer_donor_assay &&
     "${job_id}" == "${ARRAY_ID}" && -z "${extra}" ]] ||
    die "source Stage 2 jobs manifest is not bound to the accepted array ID"
  [[ "$(wc -l < "${SCHEDULER_MANIFEST_PATH}" | tr -d '[:space:]')" == 2 ]] ||
    die "source Stage 2 scheduler manifest must contain exactly array and watchdog IDs"
  IFS=$'\t' read -r scheduler_kind scheduler_id extra < "${SCHEDULER_MANIFEST_PATH}"
  [[ "${scheduler_kind}" == ARRAY && "${scheduler_id}" == "${ARRAY_ID}" &&
     -z "${extra}" ]] || die "source Stage 2 scheduler ARRAY row is invalid"
  IFS=$'\t' read -r scheduler_kind scheduler_id extra \
    < <(sed -n '2p' "${SCHEDULER_MANIFEST_PATH}")
  [[ "${scheduler_kind}" == WATCHDOG && "${scheduler_id}" == "${WATCHDOG_ID}" &&
     -z "${extra}" ]] || die "source Stage 2 scheduler WATCHDOG row is invalid"
  [[ "$(field_value "${WATCHDOG_STATUS_PATH}" STATE 2>/dev/null || true)" == OK ]] ||
    die "source Stage 2 watchdog status is not STATE=OK"
  [[ "$(field_value "${WATCHDOG_STATUS_PATH}" RUN_ID 2>/dev/null || true)" == "${SOURCE_RUN_ID}" ]] ||
    die "source Stage 2 watchdog RUN_ID mismatch"
  scheduler_count=0
  array_count=0
  watchdog_count=0
  while IFS= read -r line || [[ -n "${line}" ]]; do
    case "${line}" in
      SCHEDULER_ID=*)
        scheduler_id="${line#*=}"
        [[ "${scheduler_id}" == "${ARRAY_ID}" ||
           "${scheduler_id}" == "${WATCHDOG_ID}" ]] ||
          die "source watchdog contains an unbound scheduler ID: ${scheduler_id}"
        if [[ "${scheduler_id}" == "${ARRAY_ID}" ]]; then
          array_count=$((array_count + 1))
        else
          watchdog_count=$((watchdog_count + 1))
        fi
        scheduler_count=$((scheduler_count + 1))
        ;;
    esac
  done < "${WATCHDOG_STATUS_PATH}"
  [[ "${scheduler_count}" == 2 && "${array_count}" == 1 &&
     "${watchdog_count}" == 1 ]] ||
    die "source watchdog must record exactly one accepted array and watchdog ID"
  validate_scheduler_binding "${WATCHDOG_STATUS_PATH}"
}

validate_artifact_binding() {
  local expected_artifact parent
  require_nonwritable "${ARTIFACT_PATH}" "immutable derivative H5AD"
  [[ "${ARTIFACT_PATH}" == *.h5ad ]] || die "derivative artifact is not an H5AD"
  [[ "$(basename "${ARTIFACT_PATH}")" == "${DERIVATIVE_NAME}" ]] ||
    die "derivative artifact filename is not the expected donor-by-assay filename"
  expected_artifact="${SCRATCH_ROOT%/}/Alzheimer/data/${DERIVATIVE_NAME}"
  [[ "${ARTIFACT_PATH}" == "${expected_artifact}" ]] ||
    die "derivative artifact is not the configured Yggdrasil canonical path"
  parent="$(dirname "${ARTIFACT_PATH}")"
  [[ -d "${parent}" && ! -L "${parent}" ]] || die "derivative artifact parent is unsafe"
  [[ "${RAW_INPUT_PATH}" != "${ARTIFACT_PATH}" ]] || die "raw and derivative paths are identical"
}

validate_config_and_inputs() {
  local configured_input configured_corrected configured_raw
  CONFIG_PATH="$(canonical_existing_file "${CONFIG_ARG}" "immutable datasets.json config")"
  path_within "${CONFIG_PATH}" "${SOURCE_ROOT}" ||
    die "datasets.json config escapes the immutable source tree"
  [[ "$(sha256_file "${CONFIG_PATH}")" == "$(printf '%s' "${SOURCE_DATASETS_SHA256}" | tr '[:upper:]' '[:lower:]')" ]] ||
    die "datasets.json does not match the current source identity"
  command -v jq >/dev/null 2>&1 || die "jq is required to resolve Alzheimer configuration"
  configured_input="$(jq -er '.Alzheimer.views.batch_effect_uncorrected.input_file_name // empty' "${CONFIG_PATH}")" ||
    die "Alzheimer uncorrected input_file_name is missing from config"
  [[ "${configured_input}" == "${DERIVATIVE_NAME}" ]] ||
    die "configured Alzheimer derivative filename is not the expected donor-by-assay filename"
  configured_corrected="$(jq -er '.Alzheimer.views.batch_effect_corrected.input_file_name // empty' "${CONFIG_PATH}" 2>/dev/null || true)"
  [[ -z "${configured_corrected}" || "${configured_corrected}" == "${DERIVATIVE_NAME}" ]] ||
    die "configured Alzheimer corrected input filename is not the expected donor-by-assay filename"
  configured_raw="$(jq -er '.Alzheimer.file_names | if type == "string" then . elif type == "array" and length == 1 then .[0] else empty end' "${CONFIG_PATH}" 2>/dev/null || true)"
  [[ -n "${configured_raw}" ]] && RAW_INPUT_NAME="${configured_raw}"
  [[ "${RAW_INPUT_NAME}" =~ ^[A-Za-z0-9_.-]+\.h5ad$ ]] || die "configured raw Alzheimer filename is unsafe"
  RAW_INPUT_PATH="${SCRATCH_ROOT%/}/Alzheimer/data/${RAW_INPUT_NAME}"
  require_regular_file "${RAW_INPUT_PATH}" "authoritative Alzheimer raw H5AD"
}

validate_semantic_arguments() {
  [[ "${EXPECTED_CELLS}" == 1395601 ]] || die "expected cell count must be 1395601"
  [[ "${EXPECTED_SAMPLES}" == 104 ]] || die "expected sample count must be 104"
  [[ "${EXPECTED_DONORS}" == 83 ]] || die "expected donor count must be 83"
  [[ "${EXPECTED_ASSAY_COUNTS}" == "10x3v3=83,10xmultiome=21" ]] ||
    die "expected assay counts do not match the approved contract"
  [[ "${EXPECTED_SEX_COUNTS}" == "female=59,male=45" ]] ||
    die "expected sex counts do not match the approved contract"
  [[ "${REQUIRE_EXAMPLES}" == 1 ]] || die "deterministic example IDs are required"
  VALIDATOR_PATH="$(canonical_existing_file "${VALIDATOR_SCRIPT_ARG}" \
    "semantic Alzheimer donor-assay validator")"
  [[ "${VALIDATOR_PATH}" == "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py" ]] ||
    die "semantic validator is not the checked-in immutable Alzheimer validator"
  PYTHON_BIN="$(canonical_existing_file "${PYTHON_BIN_ARG}" \
    "explicit snapshot Python interpreter")"
  [[ -x "${PYTHON_BIN}" ]] || die "explicit snapshot Python interpreter is not executable"
  PYTHON_SHA256="$(sha256_file "${PYTHON_BIN}")"
}

run_semantic_validation() {
  local semantic_log="$1" line success_line="" success_count=0
  local -a args=(
    --input-file "${RAW_INPUT_PATH}"
    --output-file "${ARTIFACT_PATH}"
    --expected-samples "${EXPECTED_SAMPLES}"
    --expected-donors "${EXPECTED_DONORS}"
    --expected-assay-counts "${EXPECTED_ASSAY_COUNTS}"
    --expected-sex-counts "${EXPECTED_SEX_COUNTS}"
    --require-example-ids
    --validate-only
    --skip-checksum
  )
  (cd "${SOURCE_ROOT}" && "${PYTHON_BIN}" "${VALIDATOR_PATH}" "${args[@]}") >"${semantic_log}" 2>&1 ||
    die "immutable Alzheimer donor-assay semantic validation failed"

  # The checked-in validator reports the observed source contract in its
  # bounded success line.  Treat the line as evidence, rather than copying
  # the approved constants into the acceptance report without checking them.
  while IFS= read -r line || [[ -n "${line}" ]]; do
    if [[ "${line}" == ALZHEIMER_DONOR_ASSAY_VALIDATED=1* ]]; then
      success_count=$((success_count + 1))
      success_line="${line}"
    fi
  done < "${semantic_log}"
  [[ "${success_count}" == 1 ]] ||
    die "semantic validator must emit exactly one bounded success line"
  if [[ "${success_line}" =~ ^ALZHEIMER_DONOR_ASSAY_VALIDATED=1[[:space:]]cells=([0-9]+)[[:space:]]samples=([0-9]+)$ ]]; then
    OBSERVED_CELLS="${BASH_REMATCH[1]}"
    OBSERVED_SAMPLES="${BASH_REMATCH[2]}"
  else
    die "semantic validator success line has an invalid bounded format"
  fi
  [[ "${OBSERVED_CELLS}" == "${EXPECTED_CELLS}" ]] ||
    die "semantic validator observed ${OBSERVED_CELLS} cells; expected ${EXPECTED_CELLS}"
  [[ "${OBSERVED_SAMPLES}" == "${EXPECTED_SAMPLES}" ]] ||
    die "semantic validator observed ${OBSERVED_SAMPLES} samples; expected ${EXPECTED_SAMPLES}"
}

create_acceptance_root() {
  local root_parent report metadata terminal noop
  root_parent="$(dirname "${NEW_RUN_ROOT}")"
  [[ -d "${root_parent}" && ! -L "${root_parent}" ]] ||
    die "fresh acceptance run-root parent is missing or symlinked: ${root_parent}"
  root_parent="$(realpath "${root_parent}" 2>/dev/null || true)"
  [[ -n "${root_parent}" ]] || die "fresh acceptance run-root parent cannot be canonicalized"
  [[ "${root_parent}" == "$(dirname "${NEW_RUN_ROOT}")" ]] ||
    die "fresh acceptance run-root parent is not canonical"
  [[ ! -e "${NEW_RUN_ROOT}" && ! -L "${NEW_RUN_ROOT}" ]] ||
    die "fresh acceptance run root already exists: ${NEW_RUN_ROOT}"
  case "${NEW_RUN_ROOT}" in
    "${SOURCE_RUN_ROOT}"/*) die "fresh acceptance run root is nested under the source run" ;;
  esac
  mkdir "${NEW_RUN_ROOT}" || die "could not create fresh acceptance run root"
  NEW_ROOT_CREATED=1
  mkdir "${NEW_RUN_ROOT}/manifests" "${NEW_RUN_ROOT}/reports" \
    "${NEW_RUN_ROOT}/status" || die "could not create acceptance run subdirectories"

  copy_sealed "${SOURCE_METADATA_PATH}" "${NEW_RUN_ROOT}/manifests/source.metadata"
  copy_sealed "${RUN_SOURCE_MANIFEST_PATH}" "${NEW_RUN_ROOT}/manifests/source.manifest"
  copy_sealed "${SOURCE_MANIFEST_PATH}" "${NEW_RUN_ROOT}/manifests/source.manifest.current"
  copy_sealed "${RUN_RUNTIME_IDENTITY_PATH}" "${NEW_RUN_ROOT}/manifests/runtime.identity"
  copy_sealed "${RUNTIME_IDENTITY_PATH}" "${NEW_RUN_ROOT}/manifests/runtime.identity.current"
  copy_sealed "${SOURCE_TERMINAL_PATH}" "${NEW_RUN_ROOT}/manifests/source_terminal"
  copy_sealed "${WATCHDOG_STATUS_PATH}" "${NEW_RUN_ROOT}/manifests/source_watchdog"
  copy_sealed "${STEPS_MANIFEST_PATH}" "${NEW_RUN_ROOT}/manifests/steps.tsv"
  copy_sealed "${OWNERSHIP_MANIFEST_PATH}" "${NEW_RUN_ROOT}/manifests/ownership.tsv"
  copy_sealed "${JOBS_MANIFEST_PATH}" "${NEW_RUN_ROOT}/manifests/jobs.tsv"
  copy_sealed "${SCHEDULER_MANIFEST_PATH}" "${NEW_RUN_ROOT}/manifests/scheduler_ids.tsv"
  copy_sealed "${SOURCE_OWNER_PATH}" "${NEW_RUN_ROOT}/manifests/source_owner"
  copy_sealed "${ARTIFACT_RECORD_PATH}" "${NEW_RUN_ROOT}/manifests/artifact.record"
  copy_sealed "${ARTIFACT_SIDECAR_PATH}" "${NEW_RUN_ROOT}/manifests/artifact.md5"
  if [[ -n "${PRIOR_EVIDENCE_ARG}" ]]; then
    PRIOR_EVIDENCE_ARG="$(canonical_existing_file "${PRIOR_EVIDENCE_ARG}" "prior durable-gate evidence")"
    copy_sealed "${PRIOR_EVIDENCE_ARG}" "${NEW_RUN_ROOT}/manifests/prior_evidence"
  fi
  write_checked_file "${NEW_RUN_ROOT}/manifests/semantic_validator.args" \
    "PYTHON_BIN=${PYTHON_BIN}\nVALIDATOR=${VALIDATOR_PATH}\nINPUT_FILE=${RAW_INPUT_PATH}\nOUTPUT_FILE=${ARTIFACT_PATH}\nEXPECTED_CELLS=${EXPECTED_CELLS}\nEXPECTED_SAMPLES=${EXPECTED_SAMPLES}\nEXPECTED_DONORS=${EXPECTED_DONORS}\nEXPECTED_ASSAY_COUNTS=${EXPECTED_ASSAY_COUNTS}\nEXPECTED_SEX_COUNTS=${EXPECTED_SEX_COUNTS}\nREQUIRE_EXAMPLE_IDS=1\nVALIDATE_ONLY=1\nSKIP_CHECKSUM=1\nOBSERVED_CELLS=${OBSERVED_CELLS}\nOBSERVED_SAMPLES=${OBSERVED_SAMPLES}\n"
  if [[ -s "${SEMANTIC_LOG}" ]]; then
    copy_sealed "${SEMANTIC_LOG}" "${NEW_RUN_ROOT}/reports/semantic_validation.log"
  fi

  report="${NEW_RUN_ROOT}/reports/stage2_derivative_acceptance.json"
  report_content="{\n"
  report_content+="  \"schema\":\"stage2_derivative_acceptance_v1\",\n"
  report_content+="  \"stage\":\"stage2\",\n"
  report_content+="  \"state\":\"OK\",\n"
  report_content+="  \"status\":\"NOOP_VALIDATED\",\n"
  report_content+="  \"dataset\":\"Alzheimer\",\n"
  report_content+="  \"source_run\":{\"id\":$(json_quote "${SOURCE_RUN_ID}"),\"root\":$(json_quote "${SOURCE_RUN_ROOT}"),\"metadata\":$(json_quote "${SOURCE_METADATA_PATH}"),\"terminal\":$(json_quote "${SOURCE_TERMINAL_PATH}")},\n"
  report_content+="  \"source_evidence\":{\"steps\":$(json_quote "${STEPS_MANIFEST_PATH}"),\"ownership\":$(json_quote "${OWNERSHIP_MANIFEST_PATH}"),\"jobs\":$(json_quote "${JOBS_MANIFEST_PATH}"),\"scheduler\":$(json_quote "${SCHEDULER_MANIFEST_PATH}"),\"watchdog\":$(json_quote "${WATCHDOG_STATUS_PATH}"),\"owner\":$(json_quote "${SOURCE_OWNER_PATH}"),\"artifact_record\":$(json_quote "${ARTIFACT_RECORD_PATH}"),\"artifact_sidecar\":$(json_quote "${ARTIFACT_SIDECAR_PATH}")},\n"
  report_content+="  \"source_identity\":{\"root\":$(json_quote "${SOURCE_ROOT}"),\"snapshot_manifest\":$(json_quote "${SOURCE_MANIFEST_PATH}"),\"run_manifest\":$(json_quote "${RUN_SOURCE_MANIFEST_PATH}"),\"commit\":$(json_quote "${SOURCE_COMMIT}"),\"manifest_sha256\":$(json_quote "${SOURCE_MANIFEST_SHA256}")},\n"
  report_content+="  \"runtime_identity\":{\"path\":$(json_quote "${RUNTIME_IDENTITY_PATH}"),\"image\":$(json_quote "${RUNTIME_IMAGE}"),\"manifest\":$(json_quote "${RUNTIME_MANIFEST}"),\"format\":$(json_quote "${RUNTIME_FORMAT}"),\"image_sha256\":$(json_quote "${RUNTIME_IMAGE_SHA256}"),\"manifest_sha256\":$(json_quote "${RUNTIME_MANIFEST_SHA256}"),\"python_bin\":$(json_quote "${PYTHON_BIN}"),\"python_sha256\":$(json_quote "${PYTHON_SHA256}")},\n"
  report_content+="  \"artifact\":{\"path\":$(json_quote "${ARTIFACT_PATH}"),\"filename\":$(json_quote "${DERIVATIVE_NAME}"),\"size\":${ARTIFACT_SIZE},\"md5\":$(json_quote "${ARTIFACT_MD5}"),\"sidecar\":$(json_quote "${ARTIFACT_SIDECAR_PATH}"),\"record\":$(json_quote "${ARTIFACT_RECORD_PATH}"),\"owner\":$(json_quote "${SOURCE_OWNER_PATH}")},\n"
  report_content+="  \"contract\":{\"cells\":${OBSERVED_CELLS},\"samples\":${OBSERVED_SAMPLES},\"donors\":${EXPECTED_DONORS},\"assay_counts\":{\"10x3v3\":83,\"10xmultiome\":21},\"sex_sample_counts\":{\"female\":59,\"male\":45},\"example_ids\":[\"H20.33.001_10x3v3\",\"H20.33.001_10xmultiome\"]},\n"
  report_content+="  \"scheduler\":{\"array_id\":$(json_quote "${ARRAY_ID}"),\"watchdog_id\":$(json_quote "${WATCHDOG_ID}"),\"binding\":\"prebound-existing-Yggdrasil-IDs\"},\n"
  report_content+="  \"semantic_validator\":{\"script\":$(json_quote "${VALIDATOR_PATH}"),\"python_bin\":$(json_quote "${PYTHON_BIN}"),\"validate_only\":true,\"skip_checksum\":true,\"validated\":true,\"observed_cells\":${OBSERVED_CELLS},\"observed_samples\":${OBSERVED_SAMPLES}}\n"
  report_content+="}\n"
  write_checked_file "${report}" "${report_content}"

  metadata="STAGE=stage2\nRUN_ID=${NEW_RUN_ID}\nSTATE=OK\nSTATUS=NOOP_VALIDATED\nRUN_KIND=derivative_acceptance\nSOURCE_RUN_ID=${SOURCE_RUN_ID}\nSOURCE_RUN_ROOT=${SOURCE_RUN_ROOT}\nSOURCE_ROOT=${SOURCE_ROOT}\nSOURCE_MANIFEST=${SOURCE_MANIFEST_PATH}\nRUNTIME_IDENTITY=${RUNTIME_IDENTITY_PATH}\nCONFIG=${CONFIG_PATH}\nPYTHON_BIN=${PYTHON_BIN}\nVALIDATOR_SCRIPT=${VALIDATOR_PATH}\nSEMANTIC_SKIP_CHECKSUM=1\nSOURCE_TERMINAL=${SOURCE_TERMINAL_PATH}\nSOURCE_OWNER=${SOURCE_OWNER_PATH}\nSTEPS_MANIFEST=${STEPS_MANIFEST_PATH}\nOWNERSHIP_MANIFEST=${OWNERSHIP_MANIFEST_PATH}\nJOBS_MANIFEST=${JOBS_MANIFEST_PATH}\nSCHEDULER_MANIFEST=${SCHEDULER_MANIFEST_PATH}\nWATCHDOG_STATUS=${WATCHDOG_STATUS_PATH}\nARTIFACT_RECORD=${ARTIFACT_RECORD_PATH}\nARTIFACT_SIDECAR=${ARTIFACT_SIDECAR_PATH}\nARTIFACT_PATH=${ARTIFACT_PATH}\nARTIFACT_MD5=${ARTIFACT_MD5}\nARTIFACT_SIZE=${ARTIFACT_SIZE}\nEXPECTED_CELLS=${EXPECTED_CELLS}\nEXPECTED_SAMPLES=${EXPECTED_SAMPLES}\nEXPECTED_DONORS=${EXPECTED_DONORS}\nEXPECTED_ASSAY_COUNTS=${EXPECTED_ASSAY_COUNTS}\nEXPECTED_SEX_COUNTS=${EXPECTED_SEX_COUNTS}\nREQUIRE_EXAMPLE_IDS=1\nSCHEDULER_ARRAY_ID=${ARRAY_ID}\nSCHEDULER_WATCHDOG_ID=${WATCHDOG_ID}\nREPORT=${report}\n"
  write_checked_file "${NEW_RUN_ROOT}/metadata" "${metadata}"

  terminal="STATE=OK\nRUN_ID=${NEW_RUN_ID}\nSTATUS=NOOP_VALIDATED\nREASON=validator-only Stage 2 Alzheimer donor-by-assay derivative acceptance\nSOURCE_RUN_ID=${SOURCE_RUN_ID}\nARTIFACT_PATH=${ARTIFACT_PATH}\nREPORT=${report}\n"
  write_checked_file "${NEW_RUN_ROOT}/status/terminal" "${terminal}"
  noop="STATE=NOOP_VALIDATED\nRUN_ID=${NEW_RUN_ID}\nREASON=immutable derivative already published and semantically validated\nREPORT=${report}\n"
  write_checked_file "${NEW_RUN_ROOT}/status/noop" "${noop}"

  printf 'NOOP_VALIDATED=1\n'
  printf 'STAGE2_DERIVATIVE_ACCEPTANCE_RUN_ID=%s\n' "${NEW_RUN_ID}"
  printf 'STAGE2_DERIVATIVE_ACCEPTANCE_REPORT=%s\n' "${report}"
}

# All validation happens before the fresh run root is created.  The trap only
# removes a root claimed by this invocation if terminal publication did not
# complete; source evidence is never removed or changed.
NEW_ROOT_CREATED=0
cleanup() {
  rm -f "${SEMANTIC_LOG:-}" >/dev/null 2>&1 || true
  if [[ "${NEW_ROOT_CREATED:-0}" == 1 &&
        -n "${NEW_RUN_ROOT:-}" && ! -f "${NEW_RUN_ROOT}/status/terminal" ]]; then
    rm -rf "${NEW_RUN_ROOT}" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT

require_value_after_option() {
  [[ $# -ge 2 && -n "${2:-}" ]] || die "option ${1} requires a value"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-run-id) require_value_after_option "$@"; SOURCE_RUN_ID="$2"; shift 2 ;;
    --source-run-id=*) SOURCE_RUN_ID="${1#*=}"; shift ;;
    --source-run-root|--source-root-run) require_value_after_option "$@"; SOURCE_RUN_ROOT="$2"; shift 2 ;;
    --source-run-root=*|--source-root-run=*) SOURCE_RUN_ROOT="${1#*=}"; shift ;;
    --source-root|--snapshot-root) require_value_after_option "$@"; SOURCE_ROOT="$2"; shift 2 ;;
    --source-root=*|--snapshot-root=*) SOURCE_ROOT="${1#*=}"; shift ;;
    --source-manifest|--snapshot-manifest|--source-identity) require_value_after_option "$@"; SOURCE_MANIFEST_ARG="$2"; shift 2 ;;
    --source-manifest=*|--snapshot-manifest=*|--source-identity=*) SOURCE_MANIFEST_ARG="${1#*=}"; shift ;;
    --runtime-identity|--source-runtime-identity) require_value_after_option "$@"; RUNTIME_IDENTITY_ARG="$2"; shift 2 ;;
    --runtime-identity=*|--source-runtime-identity=*) RUNTIME_IDENTITY_ARG="${1#*=}"; shift ;;
    --source-terminal|--source-terminal-path|--source-run-terminal|--prior-terminal|--terminal) require_value_after_option "$@"; SOURCE_TERMINAL_ARG="$2"; shift 2 ;;
    --source-terminal=*|--source-terminal-path=*|--source-run-terminal=*|--prior-terminal=*|--terminal=*) SOURCE_TERMINAL_ARG="${1#*=}"; shift ;;
    --source-owner|--source-owner-path|--source-run-owner|--prior-owner|--owner) require_value_after_option "$@"; SOURCE_OWNER_ARG="$2"; shift 2 ;;
    --source-owner=*|--source-owner-path=*|--source-run-owner=*|--prior-owner=*|--owner=*) SOURCE_OWNER_ARG="${1#*=}"; shift ;;
    --steps-manifest|--source-steps-manifest) require_value_after_option "$@"; STEPS_MANIFEST_ARG="$2"; shift 2 ;;
    --steps-manifest=*|--source-steps-manifest=*) STEPS_MANIFEST_ARG="${1#*=}"; shift ;;
    --ownership-manifest|--source-ownership-manifest) require_value_after_option "$@"; OWNERSHIP_MANIFEST_ARG="$2"; shift 2 ;;
    --ownership-manifest=*|--source-ownership-manifest=*) OWNERSHIP_MANIFEST_ARG="${1#*=}"; shift ;;
    --jobs-manifest|--source-jobs-manifest) require_value_after_option "$@"; JOBS_MANIFEST_ARG="$2"; shift 2 ;;
    --jobs-manifest=*|--source-jobs-manifest=*) JOBS_MANIFEST_ARG="${1#*=}"; shift ;;
    --scheduler-manifest|--source-scheduler-manifest) require_value_after_option "$@"; SCHEDULER_MANIFEST_ARG="$2"; shift 2 ;;
    --scheduler-manifest=*|--source-scheduler-manifest=*) SCHEDULER_MANIFEST_ARG="${1#*=}"; shift ;;
    --watchdog-status|--source-watchdog-status) require_value_after_option "$@"; WATCHDOG_STATUS_ARG="$2"; shift 2 ;;
    --watchdog-status=*|--source-watchdog-status=*) WATCHDOG_STATUS_ARG="${1#*=}"; shift ;;
    --artifact-record|--artifact-record-path|--source-artifact-record|--prior-artifact-record|--record) require_value_after_option "$@"; ARTIFACT_RECORD_ARG="$2"; shift 2 ;;
    --artifact-record=*|--artifact-record-path=*|--source-artifact-record=*|--prior-artifact-record=*|--record=*) ARTIFACT_RECORD_ARG="${1#*=}"; shift ;;
    --sidecar|--source-sidecar|--artifact-sidecar|--checksum-sidecar) require_value_after_option "$@"; SIDECAR_ARG="$2"; shift 2 ;;
    --sidecar=*|--source-sidecar=*|--artifact-sidecar=*|--checksum-sidecar=*) SIDECAR_ARG="${1#*=}"; shift ;;
    --artifact-path|--output-path|--derivative-path|--artifact) require_value_after_option "$@"; ARTIFACT_ARG="$2"; shift 2 ;;
    --artifact-path=*|--output-path=*|--derivative-path=*|--artifact=*) ARTIFACT_ARG="${1#*=}"; shift ;;
    --config|--datasets-json) require_value_after_option "$@"; CONFIG_ARG="$2"; shift 2 ;;
    --config=*|--datasets-json=*) CONFIG_ARG="${1#*=}"; shift ;;
    --python-bin|--interpreter) require_value_after_option "$@"; PYTHON_BIN_ARG="$2"; shift 2 ;;
    --python-bin=*|--interpreter=*) PYTHON_BIN_ARG="${1#*=}"; shift ;;
    --validator-script|--semantic-validator|--h5ad-validator) require_value_after_option "$@"; VALIDATOR_SCRIPT_ARG="$2"; shift 2 ;;
    --validator-script=*|--semantic-validator=*|--h5ad-validator=*) VALIDATOR_SCRIPT_ARG="${1#*=}"; shift ;;
    --run-id|--acceptance-run-id|--fresh-run-id|--new-run-id) require_value_after_option "$@"; NEW_RUN_ID="$2"; shift 2 ;;
    --run-id=*|--acceptance-run-id=*|--fresh-run-id=*|--new-run-id=*) NEW_RUN_ID="${1#*=}"; shift ;;
    --run-root|--output-root|--acceptance-root|--fresh-run-root|--new-run-root) require_value_after_option "$@"; NEW_RUN_ROOT="$2"; shift 2 ;;
    --run-root=*|--output-root=*|--acceptance-root=*|--fresh-run-root=*|--new-run-root=*) NEW_RUN_ROOT="${1#*=}"; shift ;;
    --scratch-root|--hpc-scratch-dir|--ygg-root) require_value_after_option "$@"; SCRATCH_ROOT_ARG="$2"; shift 2 ;;
    --scratch-root=*|--hpc-scratch-dir=*|--ygg-root=*) SCRATCH_ROOT_ARG="${1#*=}"; shift ;;
    --scheduler-array-id|--accepted-array-id|--array-id) require_value_after_option "$@"; ARRAY_ID="$2"; shift 2 ;;
    --scheduler-array-id=*|--accepted-array-id=*|--array-id=*) ARRAY_ID="${1#*=}"; shift ;;
    --scheduler-watchdog-id|--accepted-watchdog-id|--watchdog-id) require_value_after_option "$@"; WATCHDOG_ID="$2"; shift 2 ;;
    --scheduler-watchdog-id=*|--accepted-watchdog-id=*|--watchdog-id=*) WATCHDOG_ID="${1#*=}"; shift ;;
    --scheduler-ids|--accepted-scheduler-ids) require_value_after_option "$@"; IFS=',' read -r ARRAY_ID WATCHDOG_ID <<< "$2"; shift 2 ;;
    --scheduler-ids=*|--accepted-scheduler-ids=*) IFS=',' read -r ARRAY_ID WATCHDOG_ID <<< "${1#*=}"; shift ;;
    --prior-evidence|--durable-evidence) require_value_after_option "$@"; PRIOR_EVIDENCE_ARG="$2"; shift 2 ;;
    --prior-evidence=*|--durable-evidence=*) PRIOR_EVIDENCE_ARG="${1#*=}"; shift ;;
    --expected-cells) require_value_after_option "$@"; EXPECTED_CELLS="$2"; shift 2 ;;
    --expected-cells=*) EXPECTED_CELLS="${1#*=}"; shift ;;
    --expected-samples) require_value_after_option "$@"; EXPECTED_SAMPLES="$2"; shift 2 ;;
    --expected-samples=*) EXPECTED_SAMPLES="${1#*=}"; shift ;;
    --expected-donors) require_value_after_option "$@"; EXPECTED_DONORS="$2"; shift 2 ;;
    --expected-donors=*) EXPECTED_DONORS="${1#*=}"; shift ;;
    --expected-assay-counts) require_value_after_option "$@"; EXPECTED_ASSAY_COUNTS="$2"; shift 2 ;;
    --expected-assay-counts=*) EXPECTED_ASSAY_COUNTS="${1#*=}"; shift ;;
    --expected-sex-counts) require_value_after_option "$@"; EXPECTED_SEX_COUNTS="$2"; shift 2 ;;
    --expected-sex-counts=*) EXPECTED_SEX_COUNTS="${1#*=}"; shift ;;
    --require-example-ids) REQUIRE_EXAMPLES=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown argument: $1" ;;
  esac
done

[[ -n "${SOURCE_RUN_ID}" ]] || die "source run ID is required"
[[ "${SOURCE_RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die "source run ID is unsafe"
require_absolute "${SOURCE_RUN_ROOT}" "source run root"
require_absolute "${SOURCE_ROOT}" "immutable source root"
require_absolute "${SOURCE_TERMINAL_ARG}" "source terminal status"
require_absolute "${SOURCE_OWNER_ARG}" "source owner evidence"
require_absolute "${ARTIFACT_RECORD_ARG}" "source artifact record"
require_absolute "${ARTIFACT_ARG}" "derivative artifact"
require_absolute "${CONFIG_ARG}" "datasets.json config"
require_absolute "${PYTHON_BIN_ARG}" "explicit snapshot Python interpreter"
[[ -n "${VALIDATOR_SCRIPT_ARG}" ]] ||
  VALIDATOR_SCRIPT_ARG="${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py"
if [[ -z "${SOURCE_MANIFEST_ARG}" ]]; then
  SOURCE_MANIFEST_ARG="${SOURCE_ROOT%/tree}/identity/source.manifest"
fi
if [[ -z "${RUNTIME_IDENTITY_ARG}" ]]; then
  RUNTIME_IDENTITY_ARG="${SOURCE_RUN_ROOT}/manifests/runtime.identity"
fi
if [[ -z "${SIDECAR_ARG}" ]]; then
  SIDECAR_ARG="${ARTIFACT_ARG}.md5"
fi
SCRATCH_ROOT_ARG="${SCRATCH_ROOT_ARG:-${HPC_SCRATCH_DIR:-}}"
require_absolute "${SCRATCH_ROOT_ARG}" "Yggdrasil scratch root"

[[ "${ARRAY_ID}" =~ ^[1-9][0-9]*$ && "${WATCHDOG_ID}" =~ ^[1-9][0-9]*$ ]] ||
  die "both prebound scheduler IDs are required"
[[ "${ARRAY_ID}" == 4407671 && "${WATCHDOG_ID}" == 4407672 ]] ||
  die "only the existing Yggdrasil scheduler IDs 4407671/4407672 are accepted"

if [[ -z "${NEW_RUN_ROOT}" ]]; then
  [[ -n "${NEW_RUN_ID}" ]] || die "fresh acceptance run ID or run root is required"
  NEW_RUN_ROOT="${SCRATCH_ROOT_ARG%/}/_ecoda_runs/${NEW_RUN_ID}"
fi
if [[ -z "${NEW_RUN_ID}" ]]; then
  NEW_RUN_ID="$(basename "${NEW_RUN_ROOT}")"
fi
[[ "${NEW_RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die "fresh acceptance run ID is unsafe"
require_absolute "${NEW_RUN_ROOT}" "fresh acceptance run root"
[[ "$(basename "${NEW_RUN_ROOT}")" == "${NEW_RUN_ID}" ]] ||
  die "fresh acceptance run root basename does not match its ID"
SOURCE_RUN_ROOT="$(canonical_existing_dir "${SOURCE_RUN_ROOT}" "source run root")"
[[ "$(basename "${SOURCE_RUN_ROOT}")" == "${SOURCE_RUN_ID}" ]] ||
  die "source run root basename does not match source run ID"
SOURCE_ROOT="$(canonical_existing_dir "${SOURCE_ROOT}" "immutable source tree")"
[[ "${SOURCE_ROOT##*/}" == tree ]] || die "immutable source root is not a snapshot tree"
SCRATCH_ROOT="$(canonical_existing_dir "${SCRATCH_ROOT_ARG}" "Yggdrasil scratch root")"
ARTIFACT_PATH="$(canonical_existing_file "${ARTIFACT_ARG}" "current Yggdrasil derivative artifact")"

if [[ -n "${ECODA_RUN_ROOT:-}" && "${ECODA_RUN_ROOT}" != "${NEW_RUN_ROOT}" ]]; then
  die "inherited ECODA_RUN_ROOT does not match the fresh acceptance root"
fi
if [[ -n "${ECODA_RUN_ID:-}" && "${ECODA_RUN_ID}" != "${NEW_RUN_ID}" ]]; then
  die "inherited ECODA_RUN_ID does not match the fresh acceptance run"
fi

validate_source_identity
validate_runtime_identity
validate_run_identity
validate_scheduler_binding "${SOURCE_TERMINAL_PATH}"
validate_config_and_inputs
validate_artifact_binding
validate_sidecar "${SIDECAR_ARG}"
validate_artifact_record
validate_owner "${SOURCE_OWNER_ARG}"
validate_stage2_manifests
validate_semantic_arguments

SEMANTIC_LOG="$(mktemp "${TMPDIR:-/tmp}/ecoda-stage2-derivative-acceptance.XXXXXX")" ||
  die "could not create semantic validation log"
run_semantic_validation "${SEMANTIC_LOG}"

create_acceptance_root
exit 0
