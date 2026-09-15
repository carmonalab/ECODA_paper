#!/bin/bash
# Validator-only release gate for one completed Stage 3 retry-2 output.
#
# This script deliberately does not source a submitter, acquire an owner, or
# publish an artifact record.  Every input belonging to the failed source run
# is read-only evidence.  Only the new run root is created.
set -euo pipefail

SCRIPT_NAME="$(basename "$0")"

SOURCE_RUN_ID=""
SOURCE_RUN_ROOT=""
SELECTION_PATH=""
PRIOR_INSPECT_PATH=""
ACCEPTED_ARRAY_ID=""
ACCEPTED_WATCHDOG_ID=""
SUPERSEDED_IDS_CSV=""
NEW_RUN_ID=""
NEW_RUN_ROOT=""
OUTPUT_PATH_ARG=""
WATCHDOG_STATUS_ARG=""
OUTPUT_OWNERSHIP_ARG=""
OWNERS_MANIFEST_ARG=""
SCHEDULER_MANIFEST_ARG=""
ARTIFACT_RECORD_ARG=""
VALIDATOR_SCRIPT_ARG=""
PYTHON_BIN_ARG=""
CONFIG_ARG=""
EXPECTED_BATCH_CONTRACT_ARG=""
SCRATCH_ROOT_ARG=""
NAS_ROOT_ARG=""
PRIOR_INSPECT_SHA256_ARG=""
PRIOR_GATE_ID=""
PRIOR_COMMAND_DIGEST=""
PRIOR_EVENT_GENERATION=""
PRIOR_SOURCE_SELECTION_BINDING=""
PRIOR_SOURCE_MANIFEST_BINDING=""
PRIOR_INSPECT_SHA256=""
SOURCE_ROOT=""
SOURCE_MANIFEST_RUN=""
SOURCE_MANIFEST_ORIGINAL=""
SOURCE_MANIFEST_SHA256=""
SOURCE_SNAPSHOT_ROOT=""
SOURCE_COMMIT=""
SOURCE_ARCHIVE=""
RUNTIME_IDENTITY_PATH=""
RUNTIME_IMAGE=""
RUNTIME_MANIFEST=""
RUNTIME_FORMAT=""
RUNTIME_IMAGE_SHA256=""
RUNTIME_MANIFEST_SHA256=""
RUNTIME_IMAGE_SIZE=""
RUNTIME_MANIFEST_SIZE=""
RUNTIME_IDENTITY_SHA256=""
RUNTIME_IDENTITY_SIZE=""
SOURCE_MANIFEST_SIZE=""
SOURCE_FORMAT=""
SOURCE_ARCHIVE_SHA256=""
SOURCE_CONFIG_HELPER_SHA256=""
SOURCE_DATASETS_SHA256=""
SOURCE_PIXI_TOML_SHA256=""
SOURCE_PIXI_LOCK_SHA256=""
SOURCE_AUX_ROOT=""
SOURCE_SCGATE_DB_BRANCH=""
RUNTIME_IMAGE_PIXI_SHA256=""
RUNTIME_IMAGE_LOCK_SHA256=""
VALIDATED_SIDECAR_MD5=""
VALIDATED_SIDECAR_SIZE=""
VALIDATED_SIDECAR_PATH=""
OUTPUT_SIDECAR_MD5=""
OUTPUT_SIDECAR_SIZE=""
OUTPUT_SIDECAR_PATH=""
VALIDATOR_PATH=""
PYTHON_INTERPRETER=""
CONFIG_PATH=""
OUTPUT_FILE_NAME=""
OUTPUT_PATH=""
NAS_PATH=""
OWNER_DIR=""
WATCHDOG_STATUS_IDS=""
SCHEDULER_MANIFEST_IDS=""
WATCHDOG_STATUS_PATH=""
OUTPUT_OWNERSHIP_PATH=""
OWNERS_MANIFEST_PATH=""
SCHEDULER_MANIFEST_PATH=""
SOURCE_TERMINAL_PATH=""
ARTIFACT_RECORD_PATH=""
ARTIFACT_RECORD_MD5=""
ARTIFACT_RECORD_SIZE=""
EXPECTED_BATCH_CONTRACT_PATH=""
ACCEPTED_ACCOUNTING_QUERY=""
ACCEPTED_ACCOUNTING_ROWS=""

usage() {
  cat <<EOF
Usage: ${SCRIPT_NAME} --source-run-id ID --source-run-root PATH \
       --selection PATH --prior-inspect PATH \
       --accepted-array-id ID --accepted-watchdog-id ID \
       --superseded-attempt-ids ID[,ID...] \
       [--run-id ID] --run-root PATH

Validate one explicit Stage 3 selection whose retry-2 H5AD already exists.
The source run and its artifacts are read-only.  The validator makes one
sacct query for the accepted array/watchdog IDs and never invokes sbatch.

Required acceptance validator:
  --validator-script PATH      explicit semantic H5AD validator (required;
                              artifact_contract.py is only a minimal fixture/
                              watchdog check and is not acceptance evidence)
  --python-bin PATH            configured Python interpreter (required)
Optional paths:
  --watchdog-status PATH       source status/watchdog (default: source root)
  --output-ownership PATH      source manifests/output_ownership.tsv
  --owners-manifest PATH       source manifests/owners.tsv
  --scheduler-manifest PATH    source manifests/scheduler_ids.tsv
  --output-path PATH           equality check against configured H5AD path
  --artifact-record PATH       explicit source artifact record
  --config PATH                immutable snapshot datasets.json
  --scratch-root PATH          configured scratch root (default: HPC_SCRATCH_DIR)
  --nas-root PATH              configured NAS root (default: NAS_TARGET_DIR)
  --expected-batch-contract JSON-or-PATH
                              explicit corrected H5AD contract identity

Aliases accepted for recovery wrappers include --output-root for --run-root,
--accepted-retry2-array-id/--accepted-retry2-watchdog-id,
--prior-inspect-evidence, and --superseded-ids.
EOF
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

# Keep all values that become records free of line/tab delimiters.  This also
# prevents an evidence file from being used as a record-injection primitive.
require_safe_value() {
  local value="$1" label="$2"
  [[ -n "${value}" && "${value}" != *$'\n'* && "${value}" != *$'\r'* &&
     "${value}" != *$'\t'* ]] || die "${label} is empty or contains a record delimiter"
}

require_absolute() {
  local value="$1" label="$2"
  require_safe_value "${value}" "${label}"
  [[ "${value}" = /* ]] || die "${label} must be absolute: ${value}"
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
  local path="$1" label="$2" dir base canonical
  require_regular_file "${path}" "${label}"
  dir="$(cd "$(dirname "${path}")" 2>/dev/null && pwd -P)" ||
    die "cannot canonicalize ${label}: ${path}"
  base="$(basename "${path}")"
  canonical="${dir}/${base}"
  [[ "${canonical}" == "${path}" ]] ||
    die "${label} is not canonical: ${path} (expected ${canonical})"
  printf '%s\n' "${canonical}"
}

canonical_existing_dir() {
  local path="$1" label="$2" canonical
  require_regular_dir "${path}" "${label}"
  canonical="$(cd "${path}" 2>/dev/null && pwd -P)" ||
    die "cannot canonicalize ${label}: ${path}"
  [[ "${canonical}" == "${path}" ]] ||
    die "${label} is not canonical: ${path} (expected ${canonical})"
  printf '%s\n' "${canonical}"
}

field_value() {
  # Read one unique KEY=value field.  The function intentionally returns a
  # non-zero status for duplicate or missing fields.
  local path="$1" key="$2" value
  value="$(awk -v wanted="${key}" '
    index($0, wanted "=") == 1 { count++; value=substr($0, length(wanted) + 2) }
    END { if (count != 1 || value == "") exit 1; print value }
  ' "${path}")" || return 1
  printf '%s\n' "${value}"
}
sha256_file() {
  local path="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${path}" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${path}" | awk '{print $1}'
  else
    return 1
  fi
}


file_size() {
  wc -c < "$1" | tr -d '[:space:]'
}

require_nonwritable_file() {
  local path="$1" label="$2" mode=""
  mode="$(stat -c '%a' "${path}" 2>/dev/null || stat -f '%Lp' "${path}" 2>/dev/null || true)"
  [[ "${mode}" =~ ^[0-7]{3,4}$ ]] || die "cannot inspect ${label} permissions: ${path}"
  [[ "${mode}" != *[2367]* ]] || die "${label} is writable: ${path}"
}
# Identity manifests are line-oriented records.  Keep their values at least
# as strict as the snapshot producer so a copied manifest cannot smuggle
# whitespace or an additional '=' into a field.
require_manifest_value() {
  local value="$1" label="$2"
  require_safe_value "${value}" "${label}"
  [[ "${value}" != *[[:space:]]* && "${value}" != *"="* ]] ||
    die "${label} contains whitespace or '='"
}

require_nonwritable_path() {
  local path="$1" label="$2" mode=""
  [[ -e "${path}" && ! -L "${path}" ]] ||
    die "${label} is missing or symlinked: ${path}"
  mode="$(stat -c '%a' "${path}" 2>/dev/null || stat -f '%Lp' "${path}" 2>/dev/null || true)"
  [[ "${mode}" =~ ^[0-7]{3,4}$ ]] || die "cannot inspect ${label} permissions: ${path}"
  [[ "${mode}" != *[2367]* ]] || die "${label} is writable: ${path}"
}


# Validate the exact three-line sidecar used by every Stage 3 H5AD output.
validate_sidecar() {
  local path="$1" label="$2" sidecar md5 size recorded_path actual_size
  require_regular_file "${path}" "${label}"
  require_nonwritable_file "${path}" "${label}"
  sidecar="${path}.md5"
  require_regular_file "${sidecar}" "${label} checksum sidecar"
  [[ "$(wc -l < "${sidecar}" | tr -d '[:space:]')" == "3" ]] ||
    die "${label} checksum sidecar must contain exactly three rows: ${sidecar}"
  [[ "$(sed -n '1p' "${sidecar}")" == MD5=* &&
     "$(sed -n '2p' "${sidecar}")" == SIZE=* &&
     "$(sed -n '3p' "${sidecar}")" == PATH=* ]] ||
    die "${label} checksum sidecar fields are not MD5/SIZE/PATH: ${sidecar}"
  md5="$(sed -n 's/^MD5=//p' "${sidecar}")"
  size="$(sed -n 's/^SIZE=//p' "${sidecar}")"
  recorded_path="$(sed -n 's/^PATH=//p' "${sidecar}")"
  [[ "${md5}" =~ ^[[:xdigit:]]{32}$ ]] ||
    die "${label} checksum MD5 is malformed: ${sidecar}"
  [[ "${size}" =~ ^[1-9][0-9]*$ ]] ||
    die "${label} checksum SIZE is malformed: ${sidecar}"
  [[ "${recorded_path}" == "${path}" ]] ||
    die "${label} checksum PATH does not match its artifact: ${sidecar}"
  actual_size="$(file_size "${path}")"
  [[ "${size}" == "${actual_size}" ]] ||
    die "${label} checksum SIZE does not match its artifact: ${path}"
  VALIDATED_SIDECAR_MD5="${md5}"
  VALIDATED_SIDECAR_SIZE="${size}"
  VALIDATED_SIDECAR_PATH="${recorded_path}"
}

# Validate the run-bound source identity against the immutable snapshot from
# which it was copied.  This is intentionally local to this validator: a
# recovery gate must not source a mutable checkout or call a writer helper.
validate_source_snapshot() {
  local source_manifest="${SOURCE_RUN_ROOT}/manifests/source.manifest"
  local manifest_source_root="" snapshot_root="" original_manifest="" complete_marker=""
  local source_manifest_from_tree=""
  local line="" key="" value="" expected_key="" index=0
  local path="" required="" actual="" expected=""
  local -a keys=(
    FORMAT SOURCE_ROOT SOURCE_COMMIT SOURCE_ARCHIVE_PATH
    SOURCE_ARCHIVE_SHA256 CONFIG_HELPER_SHA256 DATASETS_SHA256
    PIXI_TOML_SHA256 PIXI_LOCK_SHA256 AUX_ROOT SCGATE_DB_BRANCH
  )

  SOURCE_MANIFEST_RUN="$(canonical_existing_file "${source_manifest}" \
    "source identity manifest")"
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    [[ ${index} -le ${#keys[@]} ]] ||
      die "source identity manifest has extra fields"
    expected_key="${keys[$((index - 1))]}"
    [[ "${line}" == "${expected_key}="* ]] ||
      die "source identity manifest field ${index} must be ${expected_key}"
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" ]] ||
      die "source identity manifest field ${index} is malformed"
    require_manifest_value "${value}" "source identity ${key}"
    case "${key}" in
      FORMAT) SOURCE_FORMAT="${value}" ;;
      SOURCE_ROOT) manifest_source_root="${value}" ;;
      SOURCE_COMMIT) SOURCE_COMMIT="${value}" ;;
      SOURCE_ARCHIVE_PATH) SOURCE_ARCHIVE="${value}" ;;
      SOURCE_ARCHIVE_SHA256) SOURCE_ARCHIVE_SHA256="${value}" ;;
      CONFIG_HELPER_SHA256) SOURCE_CONFIG_HELPER_SHA256="${value}" ;;
      DATASETS_SHA256) SOURCE_DATASETS_SHA256="${value}" ;;
      PIXI_TOML_SHA256) SOURCE_PIXI_TOML_SHA256="${value}" ;;
      PIXI_LOCK_SHA256) SOURCE_PIXI_LOCK_SHA256="${value}" ;;
      AUX_ROOT) SOURCE_AUX_ROOT="${value}" ;;
      SCGATE_DB_BRANCH) SOURCE_SCGATE_DB_BRANCH="${value}" ;;
    esac
  done < "${SOURCE_MANIFEST_RUN}"
  [[ ${index} -eq ${#keys[@]} ]] ||
    die "source identity manifest must contain exactly ${#keys[@]} fields"
  [[ "${SOURCE_FORMAT}" == "1" ]] ||
    die "source identity manifest FORMAT must be 1"
  require_absolute "${manifest_source_root}" "immutable source root"
  [[ "${manifest_source_root}" == */tree ]] ||
    die "immutable source root is not a snapshot tree"
  SOURCE_ROOT="$(canonical_existing_dir "${manifest_source_root}" \
    "immutable source root")"
  [[ "${SOURCE_ROOT}" == "${manifest_source_root}" ]] ||
    die "source identity root is not canonical"
  snapshot_root="${SOURCE_ROOT%/tree}"
  [[ -n "${snapshot_root}" && "${snapshot_root}" != "/" ]] ||
    die "immutable source snapshot parent is missing"
  SOURCE_SNAPSHOT_ROOT="$(canonical_existing_dir "${snapshot_root}" \
    "immutable source snapshot")"
  [[ "$(basename "${SOURCE_SNAPSHOT_ROOT}")" =~ ^[[:xdigit:]]{40}$ ]] ||
    die "immutable source snapshot is not keyed by a full commit"
  [[ "${SOURCE_COMMIT}" == "$(basename "${SOURCE_SNAPSHOT_ROOT}")" ]] ||
    die "source identity commit does not match the snapshot key"
  [[ "${SOURCE_COMMIT}" =~ ^[[:xdigit:]]{40}$ ]] ||
    die "source identity SOURCE_COMMIT is not a full commit"
  [[ "${SOURCE_ARCHIVE_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${SOURCE_CONFIG_HELPER_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${SOURCE_DATASETS_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${SOURCE_PIXI_TOML_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${SOURCE_PIXI_LOCK_SHA256}" =~ ^[[:xdigit:]]{64}$ ]] ||
    die "source identity contains a malformed SHA-256 digest"
  [[ "${SOURCE_AUX_ROOT}" == "${SOURCE_ROOT}/aux" ]] ||
    die "source identity AUX_ROOT is not the snapshot aux directory"
  require_absolute "${SOURCE_ARCHIVE}" "source snapshot archive"
  # Bind the run-owned copy to the canonical form of
  # `${SOURCE_ROOT}/../identity/source.manifest`, never to a mutable checkout.
  source_manifest_from_tree="$(cd "$(dirname "${SOURCE_ROOT}/../identity/source.manifest")" \
    2>/dev/null && pwd -P)/source.manifest" ||
    die "cannot canonicalize the snapshot source manifest path"
  original_manifest="${source_manifest_from_tree}"
  complete_marker="${SOURCE_SNAPSHOT_ROOT}/COMPLETE"
  [[ "${SOURCE_ARCHIVE}" == "${SOURCE_SNAPSHOT_ROOT}/identity/source.tar" ]] ||
    die "source identity archive is outside the snapshot identity directory"
  SOURCE_MANIFEST_ORIGINAL="$(canonical_existing_file "${original_manifest}" \
    "immutable snapshot source manifest")"
  [[ "${SOURCE_MANIFEST_ORIGINAL}" == "${SOURCE_SNAPSHOT_ROOT}/identity/source.manifest" ]] ||
    die "immutable snapshot source manifest is not bound to its source tree"
  cmp -s "${SOURCE_MANIFEST_RUN}" "${SOURCE_MANIFEST_ORIGINAL}" ||
    die "source-run manifest differs from immutable snapshot source.manifest"
  require_nonwritable_file "${SOURCE_MANIFEST_ORIGINAL}" \
    "immutable snapshot source manifest"
  require_regular_file "${complete_marker}" "snapshot COMPLETE marker"
  require_nonwritable_file "${complete_marker}" "snapshot COMPLETE marker"
  [[ "$(cat "${complete_marker}")" == "COMPLETE" ]] ||
    die "snapshot COMPLETE marker is invalid"
  require_regular_file "${SOURCE_ARCHIVE}" "immutable source archive"
  require_nonwritable_file "${SOURCE_ARCHIVE}" "immutable source archive"
  actual="$(sha256_file "${SOURCE_ARCHIVE}")" ||
    die "cannot hash immutable source archive"
  expected="$(printf '%s' "${SOURCE_ARCHIVE_SHA256}" | tr '[:upper:]' '[:lower:]')"
  [[ "${actual}" == "${expected}" ]] ||
    die "immutable source archive digest does not match source identity"
  [[ ! -e "${SOURCE_ROOT}/.git" && -d "${SOURCE_ROOT}/aux" ]] ||
    die "immutable source tree has an invalid Git/aux layout"
  require_regular_dir "${SOURCE_ROOT}/src" "immutable snapshot source tree"
  require_nonwritable_path "${SOURCE_ROOT}/src" "immutable snapshot source tree"
  for required in \
    "${SOURCE_ROOT}/config_helper.R" "${SOURCE_ROOT}/datasets.json" \
    "${SOURCE_ROOT}/pixi.toml" "${SOURCE_ROOT}/pixi.lock" \
    "${SOURCE_ROOT}/aux/scGateDB.rds" \
    "${SOURCE_ROOT}/aux/genes.blocklist.rds" \
    "${SOURCE_ROOT}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"; do
    require_regular_file "${required}" "immutable snapshot source file"
    require_nonwritable_file "${required}" "immutable snapshot source file"
  done
  for path in \
    "${SOURCE_ROOT}/config_helper.R" "${SOURCE_ROOT}/datasets.json" \
    "${SOURCE_ROOT}/pixi.toml" "${SOURCE_ROOT}/pixi.lock"; do
    actual="$(sha256_file "${path}")" ||
      die "cannot hash immutable snapshot source file: ${path}"
    case "${path}" in
      */config_helper.R) expected="${SOURCE_CONFIG_HELPER_SHA256}" ;;
      */datasets.json) expected="${SOURCE_DATASETS_SHA256}" ;;
      */pixi.toml) expected="${SOURCE_PIXI_TOML_SHA256}" ;;
      */pixi.lock) expected="${SOURCE_PIXI_LOCK_SHA256}" ;;
    esac
    expected="$(printf '%s' "${expected}" | tr '[:upper:]' '[:lower:]')"
    [[ "${actual}" == "${expected}" ]] ||
      die "immutable snapshot source-file digest does not match source identity: ${path}"
  done
  while IFS= read -r path || [[ -n "${path}" ]]; do
    require_safe_value "${path}" "immutable snapshot path"
    [[ ! -L "${path}" ]] ||
      die "immutable source snapshot contains a symlink: ${path}"
    require_nonwritable_path "${path}" "immutable source snapshot path"
  done < <(find -P "${SOURCE_SNAPSHOT_ROOT}" -print)
  validate_source_archive_tree
  SOURCE_MANIFEST_SHA256="$(sha256_file "${SOURCE_MANIFEST_RUN}")" ||
    die "cannot hash source-run identity manifest"
  SOURCE_MANIFEST_SIZE="$(file_size "${SOURCE_MANIFEST_RUN}")"
}

validate_source_archive_tree() {
  local extraction="" listing="" archive_entry="" archive_symlink=""
  extraction="$(mktemp -d "${TMP_DIR}/source-archive.XXXXXX")" ||
    die "could not create source archive extraction directory"
  listing="${TMP_DIR}/source-archive.list"
  tar -tf "${SOURCE_ARCHIVE}" > "${listing}" 2>/dev/null ||
    die "could not list immutable source archive"
  while IFS= read -r archive_entry || [[ -n "${archive_entry}" ]]; do
    [[ -n "${archive_entry}" ]] || continue
    case "${archive_entry}" in
      /*|..|../*|*/../*|*/..)
        die "immutable source archive contains a path escape: ${archive_entry}"
        ;;
    esac
  done < "${listing}"
  tar -xf "${SOURCE_ARCHIVE}" -C "${extraction}" 2>/dev/null ||
    die "could not extract immutable source archive"
  [[ ! -e "${extraction}/.git" && ! -L "${extraction}/.git" ]] ||
    die "immutable source archive contains a Git directory"
  archive_symlink="$(find -P "${extraction}" -type l -print -quit 2>/dev/null || true)"
  [[ -z "${archive_symlink}" ]] ||
    die "immutable source archive contains a symlink: ${archive_symlink}"
  diff -qr "${SOURCE_ROOT}" "${extraction}" >/dev/null 2>&1 ||
    die "immutable source archive does not match the frozen source tree"
  rm -rf "${extraction}" "${listing}" || die "could not remove source archive extraction"
}

validate_runtime_binding() {
  local identity="${SOURCE_RUN_ROOT}/manifests/runtime.identity"
  local line="" key="" value="" expected_key="" index=0 identity_count=""
  local runtime_format="" manifest_format="" image_path="" image_sha=""
  local runtime_env="" runtime_layout="" runtime_prefix="" runtime_base="" runtime_pixitainer=""
  local runtime_project_root=""
  local actual_image_sha="" actual_manifest_sha="" actual_image_size="" actual_manifest_size=""
  local runtime_toml="" runtime_lock="" runtime_toml_lower="" runtime_lock_lower=""
  local source_toml_lower="" source_lock_lower="" identity_toml_lower="" identity_lock_lower=""
  local -a base_keys=(
    RUNTIME_IMAGE RUNTIME_MANIFEST RUNTIME_IMAGE_SHA256
    RUNTIME_MANIFEST_SHA256 RUNTIME_IMAGE_SIZE RUNTIME_MANIFEST_SIZE
  )
  local -a dependency_keys=(IMAGE_PIXI_TOML_SHA256 IMAGE_PIXI_LOCK_SHA256)

  RUNTIME_IDENTITY_PATH="$(canonical_existing_file "${identity}" \
    "runtime identity manifest")"
  identity_count="$(wc -l < "${RUNTIME_IDENTITY_PATH}" | tr -d '[:space:]')"
  [[ "${identity_count}" == 6 || "${identity_count}" == 8 ]] ||
    die "runtime identity has an unexpected field count"
  while IFS= read -r line || [[ -n "${line}" ]]; do
    index=$((index + 1))
    if [[ ${index} -le ${#base_keys[@]} ]]; then
      expected_key="${base_keys[$((index - 1))]}"
    else
      expected_key="${dependency_keys[$((index - ${#base_keys[@]} - 1))]}"
    fi
    [[ "${line}" == "${expected_key}="* ]] ||
      die "runtime identity field ${index} must be ${expected_key}"
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" == "${expected_key}" ]] ||
      die "runtime identity field ${index} is malformed"
    require_manifest_value "${value}" "runtime identity ${key}"
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
  [[ ${index} -eq ${identity_count} ]] ||
    die "runtime identity contains malformed trailing fields"
  require_absolute "${RUNTIME_IMAGE}" "runtime image"
  require_absolute "${RUNTIME_MANIFEST}" "runtime manifest"
  [[ "${RUNTIME_IMAGE}" == */_ecoda_runtime/*/*.sif ]] ||
    die "runtime image is not a versioned runtime image"
  [[ "${RUNTIME_IMAGE_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${RUNTIME_MANIFEST_SHA256}" =~ ^[[:xdigit:]]{64}$ &&
     "${RUNTIME_IMAGE_SIZE}" =~ ^[1-9][0-9]*$ &&
     "${RUNTIME_MANIFEST_SIZE}" =~ ^[1-9][0-9]*$ ]] ||
    die "runtime identity has a malformed digest or size"
  RUNTIME_IMAGE="$(canonical_existing_file "${RUNTIME_IMAGE}" \
    "runtime image")"
  RUNTIME_MANIFEST="$(canonical_existing_file "${RUNTIME_MANIFEST}" \
    "runtime manifest")"
  require_nonwritable_file "${RUNTIME_IMAGE}" "runtime image"
  require_nonwritable_file "${RUNTIME_MANIFEST}" "runtime manifest"
  require_nonwritable_path "$(dirname "${RUNTIME_IMAGE}")" "runtime image directory"
  [[ "${RUNTIME_MANIFEST}" == "${RUNTIME_IMAGE}.manifest" ]] ||
    die "runtime manifest is not bound beside the runtime image"
  actual_image_size="$(file_size "${RUNTIME_IMAGE}")"
  actual_manifest_size="$(file_size "${RUNTIME_MANIFEST}")"
  [[ "${actual_image_size}" == "${RUNTIME_IMAGE_SIZE}" &&
     "${actual_manifest_size}" == "${RUNTIME_MANIFEST_SIZE}" ]] ||
    die "runtime identity size does not match the immutable runtime"
  actual_image_sha="$(sha256_file "${RUNTIME_IMAGE}")" ||
    die "cannot hash immutable runtime image"
  actual_manifest_sha="$(sha256_file "${RUNTIME_MANIFEST}")" ||
    die "cannot hash immutable runtime manifest"
  [[ "${actual_image_sha}" == "$(printf '%s' "${RUNTIME_IMAGE_SHA256}" | tr '[:upper:]' '[:lower:]')" &&
     "${actual_manifest_sha}" == "$(printf '%s' "${RUNTIME_MANIFEST_SHA256}" | tr '[:upper:]' '[:lower:]')" ]] ||
    die "runtime identity digest does not match the immutable runtime"
  validate_runtime_manifest_shape "${RUNTIME_MANIFEST}"
  manifest_format="$(field_value "${RUNTIME_MANIFEST}" FORMAT)" ||
    die "runtime manifest FORMAT is missing or duplicated"
  runtime_format="${manifest_format}"
  [[ "${runtime_format}" == 1 || "${runtime_format}" == 2 ]] ||
    die "runtime manifest FORMAT is unsupported"
  [[ "${identity_count}" == 6 && "${runtime_format}" == 1 ||
     "${identity_count}" == 8 && "${runtime_format}" == 2 ]] ||
    die "runtime identity fields do not match runtime manifest FORMAT"
  image_path="$(field_value "${RUNTIME_MANIFEST}" IMAGE_PATH)" ||
    die "runtime manifest IMAGE_PATH is missing or duplicated"
  image_sha="$(field_value "${RUNTIME_MANIFEST}" IMAGE_SHA256)" ||
    die "runtime manifest IMAGE_SHA256 is missing or duplicated"
  require_absolute "${image_path}" "runtime manifest IMAGE_PATH"
  [[ "${image_path}" == "${RUNTIME_IMAGE}" &&
     "${image_sha}" =~ ^[[:xdigit:]]{64}$ &&
     "$(printf '%s' "${image_sha}" | tr '[:upper:]' '[:lower:]')" == "${actual_image_sha}" ]] ||
    die "runtime manifest image binding/digest is invalid"
  if [[ "${runtime_format}" == 2 ]]; then
    for key in IMAGE_BUILD_GIT_REVISION RUNTIME_ENV RUNTIME_LAYOUT \
      CONTAINER_ENV_PREFIX BASE_IMAGE PIXITAINER_VERSION PIXI_VERSION \
      APPTAINER_VERSION IMAGE_PIXI_TOML_SHA256 IMAGE_PIXI_LOCK_SHA256; do
      value="$(field_value "${RUNTIME_MANIFEST}" "${key}")" ||
        die "runtime manifest ${key} is missing or duplicated"
      require_manifest_value "${value}" "runtime manifest ${key}"
    done
    runtime_env="$(field_value "${RUNTIME_MANIFEST}" RUNTIME_ENV)"
    runtime_layout="$(field_value "${RUNTIME_MANIFEST}" RUNTIME_LAYOUT)"
    runtime_prefix="$(field_value "${RUNTIME_MANIFEST}" CONTAINER_ENV_PREFIX)"
    runtime_base="$(field_value "${RUNTIME_MANIFEST}" BASE_IMAGE)"
    runtime_pixitainer="$(field_value "${RUNTIME_MANIFEST}" PIXITAINER_VERSION)"
    [[ "${runtime_env}" == "py-cuda13" &&
       "${runtime_base}" == "rockylinux:9" &&
       "${runtime_pixitainer}" == "0.8.3" ]] ||
      die "runtime manifest toolchain identity is invalid"
    require_absolute "${runtime_prefix}" "runtime container prefix"
    case "${runtime_layout}" in
      relocated)
        [[ "${runtime_prefix}" == "/opt/ecoda/py-cuda13" ]] ||
          die "relocated runtime prefix is invalid"
        ;;
      path-preserving)
        runtime_project_root="$(field_value "${RUNTIME_MANIFEST}" CONTAINER_PROJECT_ROOT)" ||
          die "path-preserving runtime source root is missing"
        require_absolute "${runtime_project_root}" "path-preserving runtime source root"
        [[ "${runtime_project_root}" == "${SOURCE_ROOT}" ]] ||
          die "path-preserving runtime source root is not the snapshot root"
        ;;
      *) die "runtime manifest layout is invalid" ;;
    esac
    runtime_toml="$(field_value "${RUNTIME_MANIFEST}" IMAGE_PIXI_TOML_SHA256)"
    runtime_lock="$(field_value "${RUNTIME_MANIFEST}" IMAGE_PIXI_LOCK_SHA256)"
    runtime_toml_lower="$(printf '%s' "${runtime_toml}" | tr '[:upper:]' '[:lower:]')"
    runtime_lock_lower="$(printf '%s' "${runtime_lock}" | tr '[:upper:]' '[:lower:]')"
    source_toml_lower="$(printf '%s' "${SOURCE_PIXI_TOML_SHA256}" | tr '[:upper:]' '[:lower:]')"
    source_lock_lower="$(printf '%s' "${SOURCE_PIXI_LOCK_SHA256}" | tr '[:upper:]' '[:lower:]')"
    identity_toml_lower="$(printf '%s' "${RUNTIME_IMAGE_PIXI_SHA256}" | tr '[:upper:]' '[:lower:]')"
    identity_lock_lower="$(printf '%s' "${RUNTIME_IMAGE_LOCK_SHA256}" | tr '[:upper:]' '[:lower:]')"
    [[ "${runtime_toml}" =~ ^[[:xdigit:]]{64}$ &&
       "${runtime_lock}" =~ ^[[:xdigit:]]{64}$ &&
       "${runtime_toml_lower}" == "${source_toml_lower}" &&
       "${runtime_lock_lower}" == "${source_lock_lower}" &&
       "${identity_toml_lower}" == "${runtime_toml_lower}" &&
       "${identity_lock_lower}" == "${runtime_lock_lower}" ]] ||
      die "runtime dependency identity does not match the source snapshot"
  fi
  RUNTIME_FORMAT="${runtime_format}"
  RUNTIME_IDENTITY_SHA256="$(sha256_file "${RUNTIME_IDENTITY_PATH}")" ||
    die "cannot hash runtime identity manifest"
  RUNTIME_IDENTITY_SIZE="$(file_size "${RUNTIME_IDENTITY_PATH}")"
}

validate_runtime_manifest_shape() {
  local path="$1" line="" key="" value="" seen="" prior=""
  local -a seen_keys=()
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ "${line}" == *=* ]] || die "runtime manifest has malformed syntax"
    key="${line%%=*}"
    value="${line#*=}"
    [[ "${key}" =~ ^[A-Z][A-Z0-9_]*$ ]] ||
      die "runtime manifest has an unsafe key"
    require_manifest_value "${value}" "runtime manifest ${key}"
    if [[ ${#seen_keys[@]} -gt 0 ]]; then
      for prior in "${seen_keys[@]}"; do
        [[ "${prior}" != "${key}" ]] ||
          die "runtime manifest repeats ${key}"
      done
    fi
    seen_keys+=("${key}")
  done < "${path}"
  [[ ${#seen_keys[@]} -gt 0 ]] ||
    die "runtime manifest is empty"
}



# Validate a source-run artifact record without calling any shared writer or
# owner mutation helper.  The record is expected at the canonical digest name
# used by ecoda_run_common.sh unless an explicit path was supplied.
validate_artifact_record() {
  local path="$1" record_arg="$2" record candidate expected_path
  local actual_size recorded_run producer state sidecar_md5 sidecar_size found=0
  if [[ -n "${record_arg}" ]]; then
    record="${record_arg}"
  else
    # Locate the immutable record by its persisted PATH field.  This avoids
    # re-hashing an already published H5AD merely to reconstruct a filename.
    for candidate in "${SOURCE_RUN_ROOT}/manifests/artifacts/"*.record; do
      [[ -f "${candidate}" && ! -L "${candidate}" ]] || continue
      expected_path="$(field_value "${candidate}" PATH 2>/dev/null || true)"
      if [[ "${expected_path}" == "${path}" ]]; then
        found=$((found + 1))
        record="${candidate}"
      fi
    done
    [[ ${found} -eq 1 ]] || die "source artifact record for the selected H5AD is missing or ambiguous"
  fi
  record="$(canonical_existing_file "${record}" "source artifact record")"
  case "${record}" in
    "${SOURCE_RUN_ROOT}/manifests/artifacts/"*) ;;
    *) die "source artifact record escapes the source run root: ${record}" ;;
  esac
  [[ "$(wc -l < "${record}" | tr -d '[:space:]')" == "6" ]] ||
    die "source artifact record schema is invalid: ${record}"
  expected_path="$(field_value "${record}" PATH)" || die "source artifact record PATH is invalid"
  expected_size="$(field_value "${record}" SIZE)" || die "source artifact record SIZE is invalid"
  expected_md5="$(field_value "${record}" MD5)" || die "source artifact record MD5 is invalid"
  recorded_run="$(field_value "${record}" RUN_ID)" || die "source artifact record RUN_ID is invalid"
  producer="$(field_value "${record}" PRODUCER)" || die "source artifact record PRODUCER is invalid"
  state="$(field_value "${record}" STATE)" || die "source artifact record STATE is invalid"
  [[ "${expected_path}" == "${path}" ]] || die "source artifact record PATH mismatch: ${record}"
  [[ "${expected_size}" =~ ^[1-9][0-9]*$ ]] || die "source artifact record SIZE is malformed: ${record}"
  [[ "${expected_md5}" =~ ^[[:xdigit:]]{32}$ ]] || die "source artifact record MD5 is malformed: ${record}"
  [[ "${recorded_run}" == "${SOURCE_RUN_ID}" ]] || die "source artifact record RUN_ID mismatch: ${record}"
  [[ "${producer}" == stage3 || "${producer}" == stage3_preflight ]] ||
    die "source artifact record producer is not Stage 3: ${record}"
  [[ "${state}" == PUBLISHED ]] || die "source artifact record is not PUBLISHED: ${record}"
  actual_size="$(file_size "${path}")"
  sidecar_md5="${OUTPUT_SIDECAR_MD5:-}"
  sidecar_size="${OUTPUT_SIDECAR_SIZE:-}"
  [[ "${expected_size}" == "${actual_size}" && "${expected_size}" == "${sidecar_size}" ]] ||
    die "source artifact record SIZE does not match the immutable sidecar: ${record}"
  expected_md5_lower="$(printf '%s' "${expected_md5}" | tr '[:upper:]' '[:lower:]')"
  sidecar_md5_lower="$(printf '%s' "${sidecar_md5}" | tr '[:upper:]' '[:lower:]')"
  [[ "${expected_md5_lower}" == "${sidecar_md5_lower}" ]] ||
    die "source artifact record MD5 does not match the immutable sidecar: ${record}"
  ARTIFACT_RECORD_PATH="${record}"
  ARTIFACT_RECORD_MD5="${expected_md5}"
  ARTIFACT_RECORD_SIZE="${actual_size}"
}

# Parse the old watchdog status and bind it to the accepted retry-2 array.
validate_watchdog_status() {
  local path="$1" state run_id retry_index latest_array scheduler_id scheduler_line
  path="$(canonical_existing_file "${path}" "source watchdog status")"
  state="$(field_value "${path}" STATE)" || die "source watchdog STATE is missing or duplicated"
  run_id="$(field_value "${path}" RUN_ID)" || die "source watchdog RUN_ID is missing or duplicated"
  retry_index="$(field_value "${path}" RETRY_INDEX)" || die "source watchdog RETRY_INDEX is missing or duplicated"
  latest_array="$(field_value "${path}" ARRAY_JOB_ID)" || die "source watchdog ARRAY_JOB_ID is missing or duplicated"
  [[ "${state}" == OK ]] || die "source watchdog is not STATE=OK"
  [[ "${run_id}" == "${SOURCE_RUN_ID}" ]] || die "source watchdog RUN_ID does not match source run"
  [[ "${retry_index}" == 2 ]] || die "source watchdog is not the retry-2 watchdog"
  [[ "${latest_array}" == "${ACCEPTED_ARRAY_ID}" ]] ||
    die "source watchdog does not name the accepted retry-2 array"
  scheduler_id_count=0
  scheduler_id_seen=""
  while IFS= read -r scheduler_line || [[ -n "${scheduler_line}" ]]; do
    case "${scheduler_line}" in
      SCHEDULER_ID=*)
        scheduler_id="${scheduler_line#*=}"
        [[ "${scheduler_id}" =~ ^[0-9]+$ ]] || die "source watchdog has malformed SCHEDULER_ID"
        case " ${scheduler_id_seen} " in *" ${scheduler_id} "*) die "source watchdog repeats SCHEDULER_ID" ;; esac
        scheduler_id_seen="${scheduler_id_seen} ${scheduler_id}"
        scheduler_id_count=$((scheduler_id_count + 1))
        ;;
    esac
  done < "${path}"
  [[ ${scheduler_id_count} -gt 0 ]] || die "source watchdog has no scheduler IDs"
  case " ${scheduler_id_seen} " in
    *" ${ACCEPTED_WATCHDOG_ID} "*) ;;
    *) die "source watchdog does not record the accepted watchdog ID" ;;
  esac
  WATCHDOG_STATUS_IDS="${scheduler_id_seen}"
  WATCHDOG_STATUS_PATH="${path}"
}

validate_selection() {
  local path="$1" row ds view extra count=0 old_selection
  path="$(canonical_existing_file "${path}" "Stage 3 selection")"
  while IFS=$'\t' read -r ds view extra || [[ -n "${ds}${view}${extra}" ]]; do
    [[ -n "${ds}" && -n "${view}" && -z "${extra}" ]] || die "Stage 3 selection row must have exactly two fields"
    [[ "${ds}" =~ ^[A-Za-z0-9_.-]+$ ]] || die "Stage 3 selection dataset is unsafe: ${ds}"
    [[ "${view}" == batch_effect_uncorrected || "${view}" == batch_effect_corrected ]] ||
      die "Stage 3 selection view is not a batch-effect view: ${view}"
    count=$((count + 1))
    [[ ${count} -eq 1 ]] || die "Stage 3 retry acceptance requires exactly one selection row"
    SELECTION_DATASET="${ds}"
    SELECTION_VIEW="${view}"
  done < "${path}"
  [[ ${count} -eq 1 ]] || die "Stage 3 selection must contain exactly one row"
  SELECTION_PATH="${path}"
  old_selection="${SOURCE_RUN_ROOT}/manifests/selection.tsv"
  old_selection="$(canonical_existing_file "${old_selection}" "source run selection")"
  cmp -s "${path}" "${old_selection}" || die "explicit selection differs from source run selection"
}

resolve_configured_output_paths() {
  local config_path scratch_root nas_root output_name
  config_path="${CONFIG_ARG:-${SOURCE_ROOT}/datasets.json}"
  config_path="$(canonical_existing_file "${config_path}" "immutable snapshot datasets.json")"
  case "${config_path}" in
    "${SOURCE_ROOT}"/*) ;;
    *) die "configured datasets.json is outside the immutable source root: ${config_path}" ;;
  esac
  command -v jq >/dev/null 2>&1 || die "jq is required to resolve configured Stage 3 output paths"
  output_name="$(jq -er --arg dataset "${SELECTION_DATASET}" --arg view "${SELECTION_VIEW}" '
    .[$dataset].views[$view] as $view_config |
    if ($view_config | type) != "object" then error("view is not configured")
    else ($view_config.output_file_name // $view_config.output_file // empty)
    end
  ' "${config_path}")" || die "configured Stage 3 output filename is missing"
  [[ "${output_name}" =~ ^[A-Za-z0-9_.-]+\.h5ad$ ]] || die "configured Stage 3 output filename is unsafe: ${output_name}"
  scratch_root="${SCRATCH_ROOT_ARG:-${HPC_SCRATCH_DIR:-}}"
  nas_root="${NAS_ROOT_ARG:-${NAS_TARGET_DIR:-}}"
  require_absolute "${scratch_root}" "configured scratch root"
  require_absolute "${nas_root}" "configured NAS root"
  require_nonwritable_file "${config_path}" "immutable snapshot datasets.json"
  scratch_root="$(canonical_existing_dir "${scratch_root}" "configured scratch root")"
  nas_root="$(canonical_existing_dir "${nas_root}" "configured NAS root")"
  EXPECTED_OUTPUT_PATH="${scratch_root}/${SELECTION_DATASET}/output/${output_name}"
  EXPECTED_NAS_PATH="${nas_root}/${SELECTION_DATASET}/output/${output_name}"
  require_absolute "${EXPECTED_OUTPUT_PATH}" "configured scratch H5AD path"
  require_absolute "${EXPECTED_NAS_PATH}" "configured NAS H5AD path"
  CONFIG_PATH="${config_path}"
  OUTPUT_FILE_NAME="${output_name}"
}
validate_scheduler_manifest() {
  local path="$1" kind scheduler_id extra seen="" manifest_order=""
  local expected_order="" row_count=0 id expected_kind index=0
  [[ ${#SUPERSEDED_IDS[@]} -eq 2 ]] || die "retry acceptance requires initial and retry-1 superseded attempts"
  expected_order="${SUPERSEDED_IDS[0]} ${ACCEPTED_WATCHDOG_ID} ${SUPERSEDED_IDS[1]} ${ACCEPTED_ARRAY_ID}"
  path="$(canonical_existing_file "${path}" "source scheduler manifest")"
  while IFS=$'\t' read -r kind scheduler_id extra || [[ -n "${kind}${scheduler_id}${extra}" ]]; do
    [[ -n "${kind}" && "${scheduler_id}" =~ ^[0-9]+$ && -z "${extra}" ]] ||
      die "source scheduler manifest row is malformed"
    case " ${seen} " in *" ${scheduler_id} "*) die "source scheduler manifest repeats scheduler ID: ${scheduler_id}" ;; esac
    seen="${seen} ${scheduler_id}"
    manifest_order="${manifest_order} ${scheduler_id}"
    row_count=$((row_count + 1))
    case "${kind}" in
      ARRAY) expected_kind=ARRAY ;;
      WATCHDOG) expected_kind=WATCHDOG ;;
      STATUS) expected_kind=STATUS ;;
      *) die "source scheduler manifest has unknown kind: ${kind}" ;;
    esac
    case "${scheduler_id}" in
      "${SUPERSEDED_IDS[0]}") [[ "${kind}" == ARRAY ]] || die "initial superseded attempt is not ARRAY" ;;
      "${ACCEPTED_WATCHDOG_ID}") [[ "${kind}" == WATCHDOG ]] || die "accepted watchdog is not WATCHDOG" ;;
      "${SUPERSEDED_IDS[1]}"|"${ACCEPTED_ARRAY_ID}") [[ "${kind}" == STATUS ]] || die "retry attempt is not STATUS" ;;
      *) die "source scheduler manifest contains an unbound scheduler ID: ${scheduler_id}" ;;
    esac
  done < "${path}"
  [[ ${row_count} -eq 4 ]] || die "source scheduler manifest must contain exactly four attempt/watchdog rows"
  [[ "${manifest_order}" == " ${expected_order}" ]] || die "source scheduler manifest order does not match the full retry chain"
  for id in ${WATCHDOG_STATUS_IDS}; do
    case " ${seen} " in *" ${id} "*) ;; *) die "watchdog SCHEDULER_ID is absent from source scheduler manifest: ${id}" ;; esac
  done
  for id in ${seen}; do
    case " ${WATCHDOG_STATUS_IDS} " in *" ${id} "*) ;; *) die "source scheduler ID is absent from watchdog SCHEDULER_ID set: ${id}" ;; esac
  done
  SCHEDULER_MANIFEST_PATH="${path}"
  SCHEDULER_MANIFEST_IDS="${seen}"
}

validate_output_ownership() {
  local path="$1" ds view scratch nas owner extra count=0 owner_row
  path="$(canonical_existing_file "${path}" "source output ownership manifest")"
  while IFS=$'\t' read -r ds view scratch nas owner extra || [[ -n "${ds}${view}${scratch}${nas}${owner}${extra}" ]]; do
    [[ -n "${ds}" && -n "${view}" && -n "${scratch}" && -n "${nas}" && -n "${owner}" && -z "${extra}" ]] ||
      die "source output ownership row is malformed"
    [[ "${ds}" == "${SELECTION_DATASET}" && "${view}" == "${SELECTION_VIEW}" ]] ||
      die "source output ownership contains a row outside the one-row selection"
    count=$((count + 1))
    [[ ${count} -eq 1 ]] || die "source output ownership must contain exactly one row"
    require_absolute "${scratch}" "source scratch output path"
    require_absolute "${nas}" "source NAS output path"
    [[ "${owner}" = /* ]] || die "source output ownership owner must be an absolute path"
    OUTPUT_PATH="$(canonical_existing_file "${scratch}" "source scratch H5AD")"
    NAS_PATH="$(canonical_existing_file "${nas}" "source NAS H5AD")"
    [[ "${OUTPUT_PATH}" == "${EXPECTED_OUTPUT_PATH}" ]] ||
      die "source scratch ownership does not match configured output path"
    [[ "${NAS_PATH}" == "${EXPECTED_NAS_PATH}" ]] ||
      die "source NAS ownership does not match configured output path"
    [[ "${OUTPUT_PATH}" == *.h5ad && "${NAS_PATH}" == *.h5ad ]] ||
      die "source ownership paths must be H5AD files"
    OWNER_DIR="$(canonical_existing_dir "${owner}" "source output owner")"
  done < "${path}"
  [[ ${count} -eq 1 ]] || die "source output ownership must contain exactly one row"
  if [[ -n "${OUTPUT_PATH_ARG}" ]]; then
    local requested_output
    requested_output="$(canonical_existing_file "${OUTPUT_PATH_ARG}" "explicit scratch output")"
    [[ "${requested_output}" == "${OUTPUT_PATH}" ]] || die "explicit output path disagrees with source ownership"
  fi
  OUTPUT_OWNERSHIP_PATH="${path}"
}

validate_owner() {
  local path="$1" key owner extra count=0 owner_path owner_run owner_stage owner_state
  path="$(canonical_existing_file "${path}" "source owners manifest")"
  while IFS=$'\t' read -r key owner extra || [[ -n "${key}${owner}${extra}" ]]; do
    [[ -n "${key}" && -n "${owner}" && -z "${extra}" ]] || die "source owners manifest row is malformed"
    if [[ "${key}" == "${SELECTION_DATASET}/${SELECTION_VIEW}" ]]; then
      count=$((count + 1))
      [[ ${count} -eq 1 ]] || die "source owners manifest repeats selected output owner"
      [[ "${owner}" == "${OWNER_DIR}" ]] || die "source owners owner differs from output ownership"
    elif [[ "${key}" == ARTIFACT ]]; then
      canonical_existing_dir "${owner}" "source artifact owner evidence" >/dev/null
    else
      die "source owners manifest contains an unknown owner row: ${key}"
    fi
  done < "${path}"
  [[ ${count} -eq 1 ]] || die "source owners manifest lacks selected output owner"
  owner_path="$(field_value "${OWNER_DIR}/owner" PATH)" || die "source owner PATH is invalid"
  owner_run="$(field_value "${OWNER_DIR}/owner" RUN_ID)" || die "source owner RUN_ID is invalid"
  owner_stage="$(field_value "${OWNER_DIR}/owner" STAGE)" || die "source owner STAGE is invalid"
  owner_state="$(field_value "${OWNER_DIR}/owner" STATE)" || die "source owner STATE is invalid"
  [[ "${owner_path}" == "${OUTPUT_PATH}" ]] || die "source owner PATH mismatch"
  [[ "${owner_run}" == "${SOURCE_RUN_ID}" ]] || die "source owner RUN_ID mismatch"
  [[ "${owner_stage}" == stage3 ]] || die "source owner STAGE is not stage3"
  [[ "${owner_state}" == OK ]] || die "source output owner is not terminal OK"
  OWNERS_MANIFEST_PATH="${path}"
}

# Prior inspect evidence may be the durable gate's JSON inspect output or the
# plain accounting capture used by small recovery wrappers.  Normalize either
# into ID<TAB>STATE<TAB>EXIT_CODE before checking all attempts.
parse_prior_inspect() {
  local path="$1" first="" json_state="" release="" audit_passed="" accounting_ok="" query_count=""
  path="$(canonical_existing_file "${path}" "prior inspect evidence")"
  first="$(sed -n '/[^[:space:]]/ { s/^[[:space:]]*//; p; q; }' "${path}")"
  : > "${TMP_DIR}/prior.rows"
  : > "${TMP_DIR}/prior.meta"
  if [[ "${first}" == \{* ]]; then
    command -v jq >/dev/null 2>&1 || die "jq is required to parse JSON prior inspect evidence"
    jq -er '
      (.audit.accounting.rows // .accounting.rows // .rows) as $rows |
      if ($rows | type) != "array" then error("accounting rows are missing")
      else $rows[] |
        [(.scheduler_id // .job_id // .id // empty),
         (.state // empty), (.exit_code // .exit // empty)] | @tsv
      end
    ' "${path}" > "${TMP_DIR}/prior.rows" || die "prior inspect accounting rows are malformed"
    json_state="$(jq -r '.state // empty' "${path}")"
    release="$(jq -r '.release_eligible // empty' "${path}")"
    audit_passed="$(jq -r '.audit.passed // empty' "${path}")"
    accounting_ok="$(jq -r '.audit.accounting.ok // .accounting.ok // empty' "${path}")"
    query_count="$(jq -r '.audit.accounting.query_count // .accounting.query_count // empty' "${path}")"
    [[ "${json_state}" == FAILED ]] || die "prior inspect evidence is not failed evidence"
    [[ -z "${release}" || "${release}" == false ]] || die "prior inspect evidence is release-eligible"
    [[ "${audit_passed}" == false ]] || die "prior inspect audit did not fail closed"
    [[ -z "${accounting_ok}" || "${accounting_ok}" == false ]] || die "prior inspect accounting unexpectedly passed"
    [[ "${query_count}" == 1 ]] || die "prior inspect must contain one accounting query"
  else
    while IFS= read -r prior_line || [[ -n "${prior_line}" ]]; do
      [[ -n "${prior_line}" ]] || continue
      case "${prior_line}" in
        GATE_ID=*|COMMAND_DIGEST=*|EVENT_GENERATION=*|SOURCE_RUN_ID=*|SELECTION_PATH=*|SELECTION_SHA256=*|SOURCE_MANIFEST=*|SOURCE_MANIFEST_SHA256=*|STATE=*|RELEASE_ELIGIBLE=*|AUDIT_PASSED=*|ACCOUNTING_OK=*|QUERY_COUNT=*)
          printf '%s\n' "${prior_line}" >> "${TMP_DIR}/prior.meta"
          continue
          ;;
      esac
      printf '%s\n' "${prior_line}" | awk -F '|' 'NF == 3 {print $1 "\t" $2 "\t" $3; next} {exit 1}' >> "${TMP_DIR}/prior.rows" ||
        die "prior inspect accounting row is malformed: ${prior_line}"
    done < "${path}"
    [[ -s "${TMP_DIR}/prior.rows" ]] || die "prior inspect evidence has no accounting rows"
    json_state="$(sed -n 's/^STATE=//p' "${TMP_DIR}/prior.meta" | sed -n '1p')"
    release="$(sed -n 's/^RELEASE_ELIGIBLE=//p' "${TMP_DIR}/prior.meta" | sed -n '1p')"
    audit_passed="$(sed -n 's/^AUDIT_PASSED=//p' "${TMP_DIR}/prior.meta" | sed -n '1p')"
    accounting_ok="$(sed -n 's/^ACCOUNTING_OK=//p' "${TMP_DIR}/prior.meta" | sed -n '1p')"
    query_count="$(sed -n 's/^QUERY_COUNT=//p' "${TMP_DIR}/prior.meta" | sed -n '1p')"
    [[ "${json_state}" == FAILED ]] || die "prior inspect evidence is not failed evidence"
    [[ -z "${release}" || "${release}" == false ]] || die "prior inspect evidence is release-eligible"
    [[ "${audit_passed}" == false ]] || die "prior inspect audit did not fail closed"
    [[ -z "${accounting_ok}" || "${accounting_ok}" == false ]] || die "prior inspect accounting unexpectedly passed"
    [[ "${query_count}" == 1 ]] || die "prior inspect must contain one accounting query"
  fi
  PRIOR_INSPECT_PATH="${path}"
}
optional_field() {
  local path="$1" key="$2" value=""
  value="$(awk -v wanted="${key}" '
    index($0, wanted "=") == 1 { count++; value=substr($0, length(wanted) + 2) }
    END { if (count > 1 || (count == 1 && value == "")) exit 1; if (count == 1) print value }
  ' "${path}")" || return 1
  printf '%s\n' "${value}"
}

validate_prior_inspect_identity() {
  local first="" gate="" command_digest="" event_generation="" source_id=""
  local selection_binding="" selection_digest="" source_root_binding=""
  local source_manifest_binding="" actual_selection_digest=""
  local selection_digest_lower="" actual_selection_digest_lower=""
  local expected_prior_sha_lower="" actual_prior_sha_lower=""
  first="$(sed -n '/[^[:space:]]/ { s/^[[:space:]]*//; p; q; }' "${PRIOR_INSPECT_PATH}")"
  if [[ "${first}" == \{* ]]; then
    command -v jq >/dev/null 2>&1 || die "jq is required to bind JSON prior inspect evidence"
    gate="$(jq -r '.gate_id // empty' "${PRIOR_INSPECT_PATH}")"
    command_digest="$(jq -r '.command_digest // empty' "${PRIOR_INSPECT_PATH}")"
    event_generation="$(jq -r '.event_generation // empty' "${PRIOR_INSPECT_PATH}")"
    source_id="$(jq -r '.source_run_id // .run_id // empty' "${PRIOR_INSPECT_PATH}")"
    selection_binding="$(jq -r '(.selection.path // .selection.file // .selection_manifest // .selection_file // .manifests.selection // empty)' "${PRIOR_INSPECT_PATH}")"
    selection_digest="$(jq -r '(.selection.sha256 // .selection_sha256 // .selection.digest // empty)' "${PRIOR_INSPECT_PATH}")"
    source_root_binding="$(jq -r '(.source_run_root // .source_root // empty)' "${PRIOR_INSPECT_PATH}")"
    source_manifest_binding="$(jq -r '(.source_manifest // .source_manifest_path // .manifests.source_manifest // empty)' "${PRIOR_INSPECT_PATH}")"
  else
    gate="$(optional_field "${PRIOR_INSPECT_PATH}" GATE_ID)" || die "prior inspect GATE_ID is missing or duplicated"
    command_digest="$(optional_field "${PRIOR_INSPECT_PATH}" COMMAND_DIGEST)" || die "prior inspect COMMAND_DIGEST is missing or duplicated"
    event_generation="$(optional_field "${PRIOR_INSPECT_PATH}" EVENT_GENERATION)" || die "prior inspect EVENT_GENERATION is missing or duplicated"
    source_id="$(optional_field "${PRIOR_INSPECT_PATH}" SOURCE_RUN_ID)" ||
      die "prior inspect SOURCE_RUN_ID is duplicated"
    selection_binding="$(optional_field "${PRIOR_INSPECT_PATH}" SELECTION_PATH)" ||
      die "prior inspect SELECTION_PATH is duplicated"
    selection_digest="$(optional_field "${PRIOR_INSPECT_PATH}" SELECTION_SHA256)" ||
      die "prior inspect SELECTION_SHA256 is duplicated"
    source_root_binding=""
    source_manifest_binding="$(optional_field "${PRIOR_INSPECT_PATH}" SOURCE_MANIFEST)" ||
      die "prior inspect SOURCE_MANIFEST is duplicated"
  fi
  require_safe_value "${gate}" "prior inspect gate_id"
  require_safe_value "${command_digest}" "prior inspect command_digest"
  require_safe_value "${event_generation}" "prior inspect event_generation"
  [[ "${gate}" == "${SOURCE_RUN_ID}" ]] || die "prior inspect gate_id does not match source run ID"
  [[ "${command_digest}" =~ ^[A-Za-z0-9._:-]+$ ]] || die "prior inspect command_digest is malformed"
  [[ "${event_generation}" =~ ^[A-Za-z0-9._:-]+$ ]] || die "prior inspect event_generation is malformed"
  if [[ -n "${source_id}" ]]; then
    [[ "${source_id}" == "${SOURCE_RUN_ID}" ]] || die "prior inspect source run ID does not match"
  fi
  if [[ -n "${selection_binding}" ]]; then
    require_absolute "${selection_binding}" "prior inspect selection binding"
    selection_binding="$(canonical_existing_file "${selection_binding}" "prior inspect selection binding")"
    [[ "${selection_binding}" == "${SELECTION_PATH}" ]] ||
      die "prior inspect selection binding does not match the source selection"
    PRIOR_SOURCE_SELECTION_BINDING="${selection_binding}"
  fi
  if [[ -n "${selection_digest}" ]]; then
    [[ "${selection_digest}" =~ ^[0-9a-fA-F]{64}$ ]] || die "prior inspect selection digest is malformed"
    actual_selection_digest="$(sha256_file "${SELECTION_PATH}")" || die "cannot hash source selection for prior binding"
    selection_digest_lower="$(printf '%s' "${selection_digest}" | tr '[:upper:]' '[:lower:]')"
    actual_selection_digest_lower="$(printf '%s' "${actual_selection_digest}" | tr '[:upper:]' '[:lower:]')"
    [[ "${selection_digest_lower}" == "${actual_selection_digest_lower}" ]] ||
      die "prior inspect selection digest does not match the source selection"
  fi
  if [[ -n "${source_root_binding}" ]]; then
    require_absolute "${source_root_binding}" "prior inspect source-run-root binding"
    source_root_binding="$(canonical_existing_dir "${source_root_binding}" "prior inspect source-run-root binding")"
    [[ "${source_root_binding}" == "${SOURCE_RUN_ROOT}" ]] ||
      die "prior inspect source-run-root binding does not match"
  fi
  if [[ -n "${source_manifest_binding}" ]]; then
    require_absolute "${source_manifest_binding}" "prior inspect source-manifest binding"
    source_manifest_binding="$(canonical_existing_file "${source_manifest_binding}" "prior inspect source-manifest binding")"
    [[ "${source_manifest_binding}" == "${SOURCE_RUN_ROOT}/manifests/source.manifest" ]] ||
      die "prior inspect source-manifest binding does not match"
    PRIOR_SOURCE_MANIFEST_BINDING="${source_manifest_binding}"
  fi
  PRIOR_GATE_ID="${gate}"
  PRIOR_COMMAND_DIGEST="${command_digest}"
  PRIOR_EVENT_GENERATION="${event_generation}"
  PRIOR_INSPECT_SHA256="$(sha256_file "${PRIOR_INSPECT_PATH}")" || die "cannot hash prior inspect evidence"
  if [[ -n "${PRIOR_INSPECT_SHA256_ARG}" ]]; then
    [[ "${PRIOR_INSPECT_SHA256_ARG}" =~ ^[0-9a-fA-F]{64}$ ]] || die "expected prior inspect SHA-256 is malformed"
    expected_prior_sha_lower="$(printf '%s' "${PRIOR_INSPECT_SHA256_ARG}" | tr '[:upper:]' '[:lower:]')"
    actual_prior_sha_lower="$(printf '%s' "${PRIOR_INSPECT_SHA256}" | tr '[:upper:]' '[:lower:]')"
    [[ "${expected_prior_sha_lower}" == "${actual_prior_sha_lower}" ]] ||
      die "prior inspect SHA-256 does not match evidence"
  fi
}
validate_prior_attempts() {
  local id state exit_code extra prior_id seen="" role array_order="" expected_array_order=""
  : > "${TMP_DIR}/attempts.tsv"
  while IFS=$'\t' read -r id state exit_code extra || [[ -n "${id}${state}${exit_code}${extra}" ]]; do
    [[ -n "${id}" && -n "${state}" && -n "${exit_code}" && -z "${extra}" ]] || die "prior inspect accounting row is malformed"
    [[ "${id}" =~ ^[0-9]+$ ]] || die "prior inspect scheduler ID is malformed: ${id}"
    case " ${seen} " in *" ${id} "*) die "prior inspect repeats scheduler ID: ${id}" ;; esac
    seen="${seen} ${id}"
    role=""
    for prior_id in "${SUPERSEDED_IDS[@]}"; do
      if [[ "${id}" == "${prior_id}" ]]; then role=superseded; fi
    done
    if [[ "${id}" == "${ACCEPTED_ARRAY_ID}" ]]; then role=accepted_array; fi
    if [[ "${id}" == "${ACCEPTED_WATCHDOG_ID}" ]]; then role=accepted_watchdog; fi
    [[ -n "${role}" ]] || die "prior inspect contains an unbound scheduler ID: ${id}"
    state="${state%%+*}"
    case "${role}" in
      superseded)
        [[ "${state}" == OUT_OF_MEMORY && "${exit_code}" == 0:125 ]] ||
          die "superseded attempt ${id} is not OUT_OF_MEMORY|0:125"
        array_order="${array_order} ${id}"
        ;;
      accepted_array)
        [[ "${state}" == COMPLETED && "${exit_code}" == 0:0* ]] ||
          die "accepted retry-2 array ${id} is not COMPLETED|0:0"
        array_order="${array_order} ${id}"
        ;;
      accepted_watchdog)
        [[ "${state}" == COMPLETED && "${exit_code}" == 0:0* ]] ||
          die "accepted retry-2 watchdog ${id} is not COMPLETED|0:0"
        ;;
    esac
    printf '%s\t%s\t%s\t%s\n' "${id}" "${role}" "${state}" "${exit_code}" >> "${TMP_DIR}/attempts.tsv"
  done < "${TMP_DIR}/prior.rows"
  for id in "${SUPERSEDED_IDS[@]}" "${ACCEPTED_ARRAY_ID}" "${ACCEPTED_WATCHDOG_ID}"; do
    case " ${seen} " in *" ${id} "*) ;; *) die "prior inspect is missing scheduler ID: ${id}" ;; esac
  done
  for id in "${SUPERSEDED_IDS[@]}" "${ACCEPTED_ARRAY_ID}"; do expected_array_order="${expected_array_order} ${id}"; done
  [[ "${array_order}" == "${expected_array_order}" ]] ||
    die "prior inspect array attempts are not initial OOM, retry-1 OOM, retry-2 completed in order"
  [[ "$(wc -l < "${TMP_DIR}/attempts.tsv" | tr -d '[:space:]')" == "$(( ${#SUPERSEDED_IDS[@]} + 2 ))" ]] ||
    die "prior inspect attempt table is incomplete"
}

validate_accepted_accounting() {
  local rows row f1 f2 f3 f4 f5 id state exit_code seen="" root
  local sacct_bin="${SACCT_BIN:-sacct}"
  command -v "${sacct_bin}" >/dev/null 2>&1 || die "sacct is unavailable"
  # This is intentionally the sole accounting query.  Superseded OOM IDs are
  # accepted only from prior inspect evidence and are never queried again.
  rows="$("${sacct_bin}" -n -P -X -j "${ACCEPTED_ARRAY_ID},${ACCEPTED_WATCHDOG_ID}" \
    --format=JobIDRaw,State,ExitCode 2>/dev/null)" || die "sacct accounting query failed"
  [[ -n "${rows//[[:space:]]/}" ]] || die "sacct returned no accepted-attempt rows"
  : > "${TMP_DIR}/accepted.rows"
  while IFS= read -r row || [[ -n "${row}" ]]; do
    [[ -n "${row}" ]] || die "sacct returned a blank row"
    IFS='|' read -r f1 f2 f3 f4 f5 <<< "${row}"
    [[ -z "${f5}" ]] || die "sacct returned an extra field: ${row}"
    if [[ -n "${f4}" ]]; then
      # Four-column JobID|JobIDRaw|State|ExitCode compatibility output.
      id="${f2}"; state="${f3}"; exit_code="${f4}"
    else
      id="${f1}"; state="${f2}"; exit_code="${f3}"
    fi
    [[ "${id}" =~ ^[0-9]+$ || "${id}" =~ ^[0-9]+_[0-9]+$ ]] || die "sacct returned malformed job ID: ${id}"
    if [[ "${id}" == "${ACCEPTED_ARRAY_ID}" || "${id}" == "${ACCEPTED_WATCHDOG_ID}" ]]; then
      case " ${seen} " in *" ${id} "*) die "sacct repeated accepted root ID: ${id}" ;; esac
      seen="${seen} ${id}"
      state="${state%%+*}"
      [[ "${state}" == COMPLETED && "${exit_code}" == 0:0* ]] ||
        die "accepted accounting row is not COMPLETED|0:0: ${row}"
      printf '%s\t%s\t%s\n' "${id}" "${state}" "${exit_code}" >> "${TMP_DIR}/accepted.rows"
    elif [[ "${id}" == "${ACCEPTED_ARRAY_ID}"_* ]]; then
      # Slurm may return array children even with -X in local stubs.  Root
      # validation below remains mandatory; children are not new attempts.
      continue
    else
      die "sacct returned an unexpected accepted-attempt row: ${row}"
    fi
  done <<< "${rows}"
  case " ${seen} " in *" ${ACCEPTED_ARRAY_ID} "*) ;; *) die "sacct is missing accepted array root" ;; esac
  case " ${seen} " in *" ${ACCEPTED_WATCHDOG_ID} "*) ;; *) die "sacct is missing accepted watchdog root" ;; esac
  ACCEPTED_ACCOUNTING_QUERY="sacct -n -P -X -j ${ACCEPTED_ARRAY_ID},${ACCEPTED_WATCHDOG_ID} --format=JobIDRaw,State,ExitCode"
  ACCEPTED_ACCOUNTING_ROWS="$(cat "${TMP_DIR}/accepted.rows")"
}

build_expected_batch_contract() {
  local python_bin="$1" config_path contract_script output
  config_path="${CONFIG_ARG:-${SOURCE_ROOT}/datasets.json}"
  config_path="$(canonical_existing_file "${config_path}" "corrected validation config")"
  contract_script="${SOURCE_ROOT}/src/utils/py/batch_contract.py"
  contract_script="$(canonical_existing_file "${contract_script}" "batch contract builder")"
  output="${TMP_DIR}/expected_batch_contract.json"
  if ! PYTHONDONTWRITEBYTECODE=1 "${python_bin}" - "${contract_script}" "${config_path}" "${SELECTION_DATASET}" > "${output}" <<'PY'
import importlib.util
import json
import sys

script_path, config_path, dataset = sys.argv[1:]
spec = importlib.util.spec_from_file_location("ecoda_retry_batch_contract", script_path)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load batch contract builder")
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
with open(config_path, "r", encoding="utf-8") as handle:
    config = json.load(handle)
entry = config.get(dataset)
if not isinstance(entry, dict):
    raise ValueError("dataset is not configured")
columns = entry.get("columns")
if not isinstance(columns, dict):
    raise ValueError("dataset columns are not configured")
identity = module.build_batch_contract_identity(
    columns.get("batch"),
    sample_column="Sample",
    method_id="preprocess",
    model_id="hvg_composite_v1",
)
json.dump(identity, sys.stdout, sort_keys=True, separators=(",", ":"))
sys.stdout.write("\n")
PY
  then
    die "could not construct corrected H5AD contract identity"
  fi
  [[ -s "${output}" ]] || die "corrected H5AD contract identity is empty"
  EXPECTED_BATCH_CONTRACT_PATH="${output}"
}

validate_python_h5ad() {
  local python_bin="${PYTHON_BIN_ARG:-${PYTHON_BIN:-}}"
  local validator="${VALIDATOR_SCRIPT_ARG}"
  [[ -n "${VALIDATOR_SCRIPT_ARG}" ]] ||
    die "explicit semantic H5AD validator is required; artifact_contract.py is not acceptance evidence"
  [[ -n "${python_bin}" ]] || die "configured Python interpreter is required (use --python-bin or PYTHON_BIN)"
  if [[ "${python_bin}" = */* ]]; then
    [[ -x "${python_bin}" ]] || die "configured Python interpreter is not executable: ${python_bin}"
  else
    command -v "${python_bin}" >/dev/null 2>&1 || die "configured Python interpreter is unavailable: ${python_bin}"
  fi
  validator="$(canonical_existing_file "${validator}" "configured Python H5AD validator")"
  [[ "$(basename "${validator}")" == benchmark_h5ad_contract.py ]] ||
    die "explicit semantic validator must be benchmark_h5ad_contract.py; minimal validators are not acceptance evidence"
  if [[ "${SELECTION_VIEW}" == batch_effect_corrected ]]; then
    if [[ -n "${EXPECTED_BATCH_CONTRACT_ARG}" ]]; then
      EXPECTED_BATCH_CONTRACT_PATH="${EXPECTED_BATCH_CONTRACT_ARG}"
      if [[ "${EXPECTED_BATCH_CONTRACT_PATH}" = /* &&
            -f "${EXPECTED_BATCH_CONTRACT_PATH}" ]]; then
        EXPECTED_BATCH_CONTRACT_PATH="$(canonical_existing_file \
          "${EXPECTED_BATCH_CONTRACT_PATH}" "expected corrected H5AD contract")"
      else
        require_safe_value "${EXPECTED_BATCH_CONTRACT_PATH}" "expected corrected H5AD contract"
      fi
    else
      build_expected_batch_contract "${python_bin}"
    fi
    PYTHONDONTWRITEBYTECODE=1 "${python_bin}" "${validator}" --path "${OUTPUT_PATH}" --view "${SELECTION_VIEW}" \
      --method "Stage 3 preprocessing" \
      --expected-batch-contract "${EXPECTED_BATCH_CONTRACT_PATH}" \
      --allow-missing-corrected-summary >/dev/null 2>&1 ||
      die "configured Python H5AD validator rejected the corrected output"
  else
    PYTHONDONTWRITEBYTECODE=1 "${python_bin}" "${validator}" --path "${OUTPUT_PATH}" --view "${SELECTION_VIEW}" \
      --method "Stage 3 preprocessing" >/dev/null 2>&1 ||
      die "configured Python H5AD validator rejected the output"
  fi
  VALIDATOR_PATH="${validator}"
  PYTHON_INTERPRETER="${python_bin}"
}

write_atomic() {
  local destination="$1" content="$2" temporary
  temporary="${destination}.tmp.$$"
  printf '%b' "${content}" > "${temporary}" || return 1
  mv -f "${temporary}" "${destination}"
}


json_quote() {
  # All report values are constrained by require_safe_value.  Escape the two
  # JSON metacharacters that can still occur in safe absolute paths/tokens.
  printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'
}

create_new_run() {
  local root_parent source_copy runtime_copy selection_copy prior_copy attempts_copy
  local watchdog_copy ownership_copy owners_copy scheduler_copy source_scheduler_copy source_terminal_copy
  local report metadata terminal report_content attempts_json line id role state exit_code comma=0
  root_parent="$(dirname "${NEW_RUN_ROOT}")"
  [[ -d "${root_parent}" && ! -L "${root_parent}" ]] || die "new run root parent is missing or symlinked: ${root_parent}"
  [[ ! -e "${NEW_RUN_ROOT}" && ! -L "${NEW_RUN_ROOT}" ]] || die "new run root already exists: ${NEW_RUN_ROOT}"
  mkdir "${NEW_RUN_ROOT}" || die "could not create fresh run root: ${NEW_RUN_ROOT}"
  NEW_ROOT_CREATED=1
  mkdir "${NEW_RUN_ROOT}/manifests" "${NEW_RUN_ROOT}/status" "${NEW_RUN_ROOT}/reports" || die "could not create fresh run subdirectories"
  source_copy="${NEW_RUN_ROOT}/manifests/source.manifest"
  runtime_copy="${NEW_RUN_ROOT}/manifests/runtime.identity"
  selection_copy="${NEW_RUN_ROOT}/manifests/selection.tsv"
  prior_copy="${NEW_RUN_ROOT}/manifests/prior_inspect"
  watchdog_copy="${NEW_RUN_ROOT}/manifests/watchdog.status"
  ownership_copy="${NEW_RUN_ROOT}/manifests/output_ownership.tsv"
  owners_copy="${NEW_RUN_ROOT}/manifests/owners.tsv"
  scheduler_copy="${NEW_RUN_ROOT}/manifests/scheduler_ids.tsv"
  source_scheduler_copy="${NEW_RUN_ROOT}/manifests/source_scheduler_ids.tsv"
  source_terminal_copy="${NEW_RUN_ROOT}/manifests/source_terminal"
  attempts_copy="${NEW_RUN_ROOT}/manifests/accepted_attempts.tsv"
  cp -p "${SOURCE_RUN_ROOT}/manifests/source.manifest" "${source_copy}" || die "could not copy source identity"
  cp -p "${SOURCE_RUN_ROOT}/manifests/runtime.identity" "${runtime_copy}" || die "could not copy runtime identity"
  cp -p "${SELECTION_PATH}" "${selection_copy}" || die "could not copy selection evidence"
  cp -p "${PRIOR_INSPECT_PATH}" "${prior_copy}" || die "could not copy prior inspect evidence"
  cp -p "${WATCHDOG_STATUS_PATH}" "${watchdog_copy}" || die "could not copy watchdog evidence"
  cp -p "${OUTPUT_OWNERSHIP_PATH}" "${ownership_copy}" || die "could not copy output ownership evidence"
  cp -p "${OWNERS_MANIFEST_PATH}" "${owners_copy}" || die "could not copy owners evidence"
  cp -p "${SCHEDULER_MANIFEST_PATH}" "${source_scheduler_copy}" || die "could not copy scheduler evidence"
  write_atomic "${scheduler_copy}" \
    "ARRAY\t${ACCEPTED_ARRAY_ID}\nWATCHDOG\t${ACCEPTED_WATCHDOG_ID}\n" ||
    die "could not write accepted scheduler manifest"
  cp -p "${SOURCE_TERMINAL_PATH}" "${source_terminal_copy}" || die "could not copy source terminal evidence"
  {
    printf 'ATTEMPT_ID\tROLE\tSTATE\tEXIT_CODE\n'
    cat "${TMP_DIR}/attempts.tsv"
  } > "${attempts_copy}" || die "could not write accepted attempt table"
  chmod a-w "${source_copy}" "${runtime_copy}" "${selection_copy}" "${prior_copy}" "${attempts_copy}" \
    "${watchdog_copy}" "${ownership_copy}" "${owners_copy}" "${scheduler_copy}" \
    "${source_scheduler_copy}" "${source_terminal_copy}" ||
    die "could not seal copied acceptance evidence"
  report="${NEW_RUN_ROOT}/reports/stage3_retry_acceptance.json"
  attempts_json=""
  while IFS=$'\t' read -r id role state exit_code; do
    [[ -n "${id}" ]] || continue
    [[ ${comma} -eq 0 ]] || attempts_json+=","
    attempts_json+="{\"attempt_id\":\"$(json_quote "${id}")\",\"role\":\"$(json_quote "${role}")\",\"state\":\"$(json_quote "${state}")\",\"exit_code\":\"$(json_quote "${exit_code}")\"}"
    comma=1
  done < "${TMP_DIR}/attempts.tsv"
  report_content="{\n  \"schema\":\"stage3_retry_acceptance_v1\",\n  \"stage\":\"stage3\",\n  \"state\":\"OK\",\n  \"source_metadata_state\":\"$(json_quote "${SOURCE_METADATA_STATE}")\",\n  \"source_metadata_control_plane_note\":\"ACTIVE is permitted as preserved source-run control-plane evidence; watchdog/artifact contracts are authoritative\",\n  \"source_run_id\":\"$(json_quote "${SOURCE_RUN_ID}")\",\n  \"source_run_root\":\"$(json_quote "${SOURCE_RUN_ROOT}")\",\n  \"source_identity\":{\"run_manifest\":\"$(json_quote "${SOURCE_MANIFEST_RUN}")\",\"snapshot_manifest\":\"$(json_quote "${SOURCE_MANIFEST_ORIGINAL}")\",\"snapshot_root\":\"$(json_quote "${SOURCE_SNAPSHOT_ROOT}")\",\"commit\":\"$(json_quote "${SOURCE_COMMIT}")\",\"sha256\":\"$(json_quote "${SOURCE_MANIFEST_SHA256}")\",\"size\":\"$(json_quote "${SOURCE_MANIFEST_SIZE}")\"},\n  \"runtime_identity\":{\"path\":\"$(json_quote "${RUNTIME_IDENTITY_PATH}")\",\"sha256\":\"$(json_quote "${RUNTIME_IDENTITY_SHA256}")\",\"size\":\"$(json_quote "${RUNTIME_IDENTITY_SIZE}")\",\"format\":\"$(json_quote "${RUNTIME_FORMAT}")\",\"image\":\"$(json_quote "${RUNTIME_IMAGE}")\",\"manifest\":\"$(json_quote "${RUNTIME_MANIFEST}")\",\"image_sha256\":\"$(json_quote "${RUNTIME_IMAGE_SHA256}")\",\"manifest_sha256\":\"$(json_quote "${RUNTIME_MANIFEST_SHA256}")\"},\n  \"validator\":{\"path\":\"$(json_quote "${VALIDATOR_PATH}")\",\"semantic\":true},\n  \"run_id\":\"$(json_quote "${NEW_RUN_ID}")\",\n  \"run_root\":\"$(json_quote "${NEW_RUN_ROOT}")\",\n  \"selection\":{\"dataset\":\"$(json_quote "${SELECTION_DATASET}")\",\"view\":\"$(json_quote "${SELECTION_VIEW}")\",\"path\":\"$(json_quote "${selection_copy}")\",\"config\":\"$(json_quote "${CONFIG_PATH}")\",\"output_file_name\":\"$(json_quote "${OUTPUT_FILE_NAME}")\"},\n  \"output\":{\"scratch_path\":\"$(json_quote "${OUTPUT_PATH}")\",\"nas_path\":\"$(json_quote "${NAS_PATH}")\",\"md5\":\"$(json_quote "${ARTIFACT_RECORD_MD5}")\",\"size\":\"$(json_quote "${ARTIFACT_RECORD_SIZE}")\",\"artifact_record\":\"$(json_quote "${ARTIFACT_RECORD_PATH}")\"},\n  \"accepted_retry2\":{\"array_id\":\"$(json_quote "${ACCEPTED_ARRAY_ID}")\",\"watchdog_id\":\"$(json_quote "${ACCEPTED_WATCHDOG_ID}")\",\"accounting_query\":\"$(json_quote "${ACCEPTED_ACCOUNTING_QUERY}")\"},\n  \"superseded_attempt_ids\":[$(printf '\"%s\",' "${SUPERSEDED_IDS[@]}" | sed 's/,$//')],\n  \"attempts\":[${attempts_json}],\n  \"prior_inspect_identity\":{\"path\":\"$(json_quote "${prior_copy}")\",\"sha256\":\"$(json_quote "${PRIOR_INSPECT_SHA256}")\",\"gate_id\":\"$(json_quote "${PRIOR_GATE_ID}")\",\"command_digest\":\"$(json_quote "${PRIOR_COMMAND_DIGEST}")\",\"event_generation\":\"$(json_quote "${PRIOR_EVENT_GENERATION}")\",\"selection_binding\":\"$(json_quote "${PRIOR_SOURCE_SELECTION_BINDING:-}")\",\"source_manifest_binding\":\"$(json_quote "${PRIOR_SOURCE_MANIFEST_BINDING:-}")\"},\n  \"evidence\":{\"watchdog_status\":\"$(json_quote "${watchdog_copy}")\",\"output_ownership\":\"$(json_quote "${ownership_copy}")\",\"owners_manifest\":\"$(json_quote "${owners_copy}")\",\"scheduler_manifest\":\"$(json_quote "${scheduler_copy}")\",\"source_manifest\":\"$(json_quote "${source_copy}")\",\"runtime_identity\":\"$(json_quote "${runtime_copy}")\"}\n}\n"
  write_atomic "${report}" "${report_content}" || die "could not write acceptance report"
  metadata="${NEW_RUN_ROOT}/metadata"
  write_atomic "${metadata}" "STAGE=stage3\nRUN_ID=${NEW_RUN_ID}\nSTATE=OK\nRUN_KIND=retry_acceptance\nSOURCE_METADATA_STATE=${SOURCE_METADATA_STATE}\nSOURCE_RUN_ID=${SOURCE_RUN_ID}\nSOURCE_RUN_ROOT=${SOURCE_RUN_ROOT}\nSOURCE_MANIFEST_RUN=${SOURCE_MANIFEST_RUN}\nSOURCE_MANIFEST_ORIGINAL=${SOURCE_MANIFEST_ORIGINAL}\nSOURCE_SNAPSHOT_ROOT=${SOURCE_SNAPSHOT_ROOT}\nSOURCE_COMMIT=${SOURCE_COMMIT}\nSOURCE_MANIFEST_SHA256=${SOURCE_MANIFEST_SHA256}\nSOURCE_MANIFEST_SIZE=${SOURCE_MANIFEST_SIZE}\nRUNTIME_IDENTITY=${RUNTIME_IDENTITY_PATH}\nRUNTIME_IDENTITY_SHA256=${RUNTIME_IDENTITY_SHA256}\nRUNTIME_IDENTITY_SIZE=${RUNTIME_IDENTITY_SIZE}\nRUNTIME_FORMAT=${RUNTIME_FORMAT}\nRUNTIME_IMAGE=${RUNTIME_IMAGE}\nRUNTIME_MANIFEST=${RUNTIME_MANIFEST}\nRUNTIME_IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}\nRUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA256}\nRUNTIME_IMAGE_SIZE=${RUNTIME_IMAGE_SIZE}\nRUNTIME_MANIFEST_SIZE=${RUNTIME_MANIFEST_SIZE}\nVALIDATOR=${VALIDATOR_PATH}\nSELECTION=${selection_copy}\nSELECTION_DATASET=${SELECTION_DATASET}\nSELECTION_VIEW=${SELECTION_VIEW}\nOUTPUT_PATH=${OUTPUT_PATH}\nNAS_PATH=${NAS_PATH}\nOUTPUT_MD5=${ARTIFACT_RECORD_MD5}\nOUTPUT_SIZE=${ARTIFACT_RECORD_SIZE}\nACCEPTED_ARRAY_ID=${ACCEPTED_ARRAY_ID}\nACCEPTED_WATCHDOG_ID=${ACCEPTED_WATCHDOG_ID}\nSUPERSEDED_ATTEMPT_IDS=${SUPERSEDED_IDS_CSV}\nPRIOR_INSPECT=${prior_copy}\nPRIOR_INSPECT_SHA256=${PRIOR_INSPECT_SHA256}\nPRIOR_GATE_ID=${PRIOR_GATE_ID}\nPRIOR_COMMAND_DIGEST=${PRIOR_COMMAND_DIGEST}\nPRIOR_EVENT_GENERATION=${PRIOR_EVENT_GENERATION}\nWATCHDOG_STATUS=${WATCHDOG_STATUS_PATH}\nACCEPTANCE_REPORT=${report}\n"
  terminal="${NEW_RUN_ROOT}/status/terminal"
  write_atomic "${terminal}" "STATE=OK\nRUN_ID=${NEW_RUN_ID}\nREASON=validator-only Stage 3 retry-2 acceptance\nREPORT=${report}\nSOURCE_RUN_ID=${SOURCE_RUN_ID}\nACCEPTED_ARRAY_ID=${ACCEPTED_ARRAY_ID}\nACCEPTED_WATCHDOG_ID=${ACCEPTED_WATCHDOG_ID}\n"
  printf 'STAGE3_RETRY_ACCEPTANCE_RUN_ID=%s\nSTAGE3_RETRY_ACCEPTANCE_REPORT=%s\n' "${NEW_RUN_ID}" "${report}"
}

# Parse arguments before creating anything.  A malformed attempt therefore
# cannot even create an empty recovery run and, importantly, cannot reach sacct.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-run-id) SOURCE_RUN_ID="${2:-}"; shift 2 ;;
    --source-run-id=*) SOURCE_RUN_ID="${1#*=}"; shift ;;
    --source-run-root|--source-root) SOURCE_RUN_ROOT="${2:-}"; shift 2 ;;
    --source-run-root=*|--source-root=*) SOURCE_RUN_ROOT="${1#*=}"; shift ;;
    --selection|--selection-file|--source-selection) SELECTION_PATH="${2:-}"; shift 2 ;;
    --selection=*|--selection-file=*|--source-selection=*) SELECTION_PATH="${1#*=}"; shift ;;
    --prior-inspect|--prior-inspect-evidence|--inspect-evidence|--inspect) PRIOR_INSPECT_PATH="${2:-}"; shift 2 ;;
    --prior-inspect=*|--prior-inspect-evidence=*|--inspect-evidence=*|--inspect=*) PRIOR_INSPECT_PATH="${1#*=}"; shift ;;
    --accepted-array-id|--accepted-retry2-array-id|--accepted-retry2-array) ACCEPTED_ARRAY_ID="${2:-}"; shift 2 ;;
    --accepted-array-id=*|--accepted-retry2-array-id=*|--accepted-retry2-array=*) ACCEPTED_ARRAY_ID="${1#*=}"; shift ;;
    --accepted-watchdog-id|--accepted-retry2-watchdog-id|--accepted-retry2-watchdog) ACCEPTED_WATCHDOG_ID="${2:-}"; shift 2 ;;
    --accepted-watchdog-id=*|--accepted-retry2-watchdog-id=*|--accepted-retry2-watchdog=*) ACCEPTED_WATCHDOG_ID="${1#*=}"; shift ;;
    --accepted-retry2-ids)
      IFS=',' read -r ACCEPTED_ARRAY_ID ACCEPTED_WATCHDOG_ID <<< "${2:-}"; shift 2 ;;
    --accepted-retry2-ids=*) IFS=',' read -r ACCEPTED_ARRAY_ID ACCEPTED_WATCHDOG_ID <<< "${1#*=}"; shift ;;
    --superseded-attempt-id) if [[ -n "${SUPERSEDED_IDS_CSV}" ]]; then SUPERSEDED_IDS_CSV+=",${2:-}"; else SUPERSEDED_IDS_CSV="${2:-}"; fi; shift 2 ;;
    --superseded-attempt-id=*) if [[ -n "${SUPERSEDED_IDS_CSV}" ]]; then SUPERSEDED_IDS_CSV+=",${1#*=}"; else SUPERSEDED_IDS_CSV="${1#*=}"; fi; shift ;;
    --superseded-attempt-ids|--superseded-ids) SUPERSEDED_IDS_CSV="${2:-}"; shift 2 ;;
    --superseded-attempt-ids=*|--superseded-ids=*) SUPERSEDED_IDS_CSV="${1#*=}"; shift ;;
    --run-id|--recovery-run-id|--new-run-id) NEW_RUN_ID="${2:-}"; shift 2 ;;
    --run-id=*|--recovery-run-id=*|--new-run-id=*) NEW_RUN_ID="${1#*=}"; shift ;;
    --run-root|--output-root|--recovery-run-root|--new-run-root) NEW_RUN_ROOT="${2:-}"; shift 2 ;;
    --run-root=*|--output-root=*|--recovery-run-root=*|--new-run-root=*) NEW_RUN_ROOT="${1#*=}"; shift ;;
    --output-path|--h5ad) OUTPUT_PATH_ARG="${2:-}"; shift 2 ;;
    --output-path=*|--h5ad=*) OUTPUT_PATH_ARG="${1#*=}"; shift ;;
    --watchdog-status) WATCHDOG_STATUS_ARG="${2:-}"; shift 2 ;;
    --watchdog-status=*) WATCHDOG_STATUS_ARG="${1#*=}"; shift ;;
    --output-ownership) OUTPUT_OWNERSHIP_ARG="${2:-}"; shift 2 ;;
    --output-ownership=*) OUTPUT_OWNERSHIP_ARG="${1#*=}"; shift ;;
    --owners-manifest) OWNERS_MANIFEST_ARG="${2:-}"; shift 2 ;;
    --owners-manifest=*) OWNERS_MANIFEST_ARG="${1#*=}"; shift ;;
    --scheduler-manifest) SCHEDULER_MANIFEST_ARG="${2:-}"; shift 2 ;;
    --scheduler-manifest=*) SCHEDULER_MANIFEST_ARG="${1#*=}"; shift ;;
    --artifact-record) ARTIFACT_RECORD_ARG="${2:-}"; shift 2 ;;
    --artifact-record=*) ARTIFACT_RECORD_ARG="${1#*=}"; shift ;;
    --validator-script|--h5ad-validator) VALIDATOR_SCRIPT_ARG="${2:-}"; shift 2 ;;
    --validator-script=*|--h5ad-validator=*) VALIDATOR_SCRIPT_ARG="${1#*=}"; shift ;;
    --python-bin) PYTHON_BIN_ARG="${2:-}"; shift 2 ;;
    --python-bin=*) PYTHON_BIN_ARG="${1#*=}"; shift ;;
    --config) CONFIG_ARG="${2:-}"; shift 2 ;;
    --config=*) CONFIG_ARG="${1#*=}"; shift ;;
    --expected-batch-contract) EXPECTED_BATCH_CONTRACT_ARG="${2:-}"; shift 2 ;;
    --expected-batch-contract=*) EXPECTED_BATCH_CONTRACT_ARG="${1#*=}"; shift ;;
    --scratch-root) SCRATCH_ROOT_ARG="${2:-}"; shift 2 ;;
    --scratch-root=*) SCRATCH_ROOT_ARG="${1#*=}"; shift ;;
    --nas-root) NAS_ROOT_ARG="${2:-}"; shift 2 ;;
    --nas-root=*) NAS_ROOT_ARG="${1#*=}"; shift ;;
    --prior-inspect-sha256) PRIOR_INSPECT_SHA256_ARG="${2:-}"; shift 2 ;;
    --prior-inspect-sha256=*) PRIOR_INSPECT_SHA256_ARG="${1#*=}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown argument: $1" ;;
  esac
done

require_safe_value "${SOURCE_RUN_ID}" "source run ID"
[[ "${SOURCE_RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die "source run ID is unsafe: ${SOURCE_RUN_ID}"
require_absolute "${SOURCE_RUN_ROOT}" "source run root"
require_absolute "${SELECTION_PATH}" "selection"
require_absolute "${PRIOR_INSPECT_PATH}" "prior inspect evidence"
[[ -n "${VALIDATOR_SCRIPT_ARG}" ]] ||
  die "explicit semantic H5AD validator is required; artifact_contract.py is not acceptance evidence"
[[ "$(basename "${VALIDATOR_SCRIPT_ARG}")" == benchmark_h5ad_contract.py ]] ||
  die "explicit semantic validator must be benchmark_h5ad_contract.py; minimal validators are not acceptance evidence"
[[ "${ACCEPTED_ARRAY_ID}" =~ ^[1-9][0-9]*$ ]] || die "accepted array ID is invalid"
[[ "${ACCEPTED_WATCHDOG_ID}" =~ ^[1-9][0-9]*$ ]] || die "accepted watchdog ID is invalid"
[[ "${ACCEPTED_ARRAY_ID}" != "${ACCEPTED_WATCHDOG_ID}" ]] || die "accepted array/watchdog IDs must differ"
require_safe_value "${SUPERSEDED_IDS_CSV}" "superseded attempt IDs"
[[ -n "${SUPERSEDED_IDS_CSV}" ]] || die "superseded attempt IDs must not be empty"
IFS=',' read -r -a SUPERSEDED_IDS <<< "${SUPERSEDED_IDS_CSV}"
[[ ${#SUPERSEDED_IDS[@]} -gt 0 ]] || die "superseded attempt IDs must not be empty"
for superseded_id in "${SUPERSEDED_IDS[@]}"; do
  [[ "${superseded_id}" =~ ^[1-9][0-9]*$ ]] || die "superseded attempt ID is invalid: ${superseded_id}"
  [[ "${superseded_id}" != "${ACCEPTED_ARRAY_ID}" && "${superseded_id}" != "${ACCEPTED_WATCHDOG_ID}" ]] || die "accepted ID is also superseded"
done
# Reject duplicate superseded IDs before any evidence query.
for ((superseded_i = 0; superseded_i < ${#SUPERSEDED_IDS[@]}; superseded_i++)); do
  for ((superseded_j = superseded_i + 1; superseded_j < ${#SUPERSEDED_IDS[@]}; superseded_j++)); do
    [[ "${SUPERSEDED_IDS[$superseded_i]}" != "${SUPERSEDED_IDS[$superseded_j]}" ]] || die "superseded attempt IDs contain a duplicate"
  done
done

SOURCE_RUN_ROOT="$(canonical_existing_dir "${SOURCE_RUN_ROOT}" "source run root")"
[[ "$(basename "${SOURCE_RUN_ROOT}")" == "${SOURCE_RUN_ID}" ]] || die "source run root basename does not match source run ID"
require_regular_file "${SOURCE_RUN_ROOT}/metadata" "source run metadata"
require_regular_file "${SOURCE_RUN_ROOT}/manifests/source.manifest" "source identity manifest"
require_regular_file "${SOURCE_RUN_ROOT}/manifests/runtime.identity" "runtime identity manifest"
[[ "$(field_value "${SOURCE_RUN_ROOT}/metadata" STAGE 2>/dev/null || true)" == stage3 ]] || die "source run is not a Stage 3 run"
[[ "$(field_value "${SOURCE_RUN_ROOT}/metadata" RUN_ID 2>/dev/null || true)" == "${SOURCE_RUN_ID}" ]] || die "source run metadata RUN_ID mismatch"
SOURCE_METADATA_STATE="$(field_value "${SOURCE_RUN_ROOT}/metadata" STATE 2>/dev/null || true)"
case "${SOURCE_METADATA_STATE}" in
  ACTIVE|OK) ;;
  *) die "source run metadata has an invalid terminal/control-plane STATE: ${SOURCE_METADATA_STATE}" ;;
esac
SOURCE_TERMINAL_PATH="${SOURCE_RUN_ROOT}/status/terminal"
require_regular_file "${SOURCE_TERMINAL_PATH}" "source terminal status"
[[ "$(field_value "${SOURCE_TERMINAL_PATH}" STATE 2>/dev/null || true)" == OK ]] ||
  die "source terminal status is not STATE=OK"
[[ "$(field_value "${SOURCE_TERMINAL_PATH}" RUN_ID 2>/dev/null || true)" == "${SOURCE_RUN_ID}" ]] ||
  die "source terminal status RUN_ID mismatch"
# If no explicit fresh ID is supplied, the run-root basename is the ID.  The
# explicit forms are both accepted so durable wrappers can bind both values.
if [[ -z "${NEW_RUN_ROOT}" ]]; then
  [[ -n "${NEW_RUN_ID}" ]] || die "--run-root or --run-id is required"
  local_default_root="${ECODA_RUNS_ROOT:-${HPC_SCRATCH_DIR:-${HOME:-/tmp}/scratch/ECODA_paper}/_ecoda_runs}/${NEW_RUN_ID}"
  NEW_RUN_ROOT="${local_default_root}"
fi
require_absolute "${NEW_RUN_ROOT}" "new run root"
if [[ -z "${NEW_RUN_ID}" ]]; then NEW_RUN_ID="$(basename "${NEW_RUN_ROOT}")"; fi
[[ "${NEW_RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die "new run ID is unsafe"
[[ "$(basename "${NEW_RUN_ROOT}")" == "${NEW_RUN_ID}" ]] || die "new run root basename does not match new run ID"
[[ "${NEW_RUN_ROOT}" != "${SOURCE_RUN_ROOT}" ]] || die "new run root must differ from source run root"
NEW_RUN_ROOT_PARENT="$(dirname "${NEW_RUN_ROOT}")"
case "${NEW_RUN_ROOT}" in
  "${SOURCE_RUN_ROOT}"/*) die "new run root must not be nested under source run root" ;;
esac
require_absolute "${NEW_RUN_ROOT_PARENT}" "new run root parent"
[[ -d "${NEW_RUN_ROOT_PARENT}" && ! -L "${NEW_RUN_ROOT_PARENT}" ]] || die "new run root parent must already exist"

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage3-retry-accept.XXXXXX")" || die "could not create temporary evidence directory"
NEW_ROOT_CREATED=0
cleanup() {
  rm -rf "${TMP_DIR:-}" >/dev/null 2>&1 || true
  if [[ "${NEW_ROOT_CREATED:-0}" == 1 && ! -f "${NEW_RUN_ROOT:-}/status/terminal" ]]; then
    # Only remove a root created by this invocation when publication did not
    # reach terminal OK; never touch the source run or any pre-existing root.
    rm -rf "${NEW_RUN_ROOT}" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT
validate_source_snapshot
validate_runtime_binding

validate_selection "${SELECTION_PATH}"
resolve_configured_output_paths
validate_watchdog_status "${WATCHDOG_STATUS_ARG:-${SOURCE_RUN_ROOT}/status/watchdog}"
validate_scheduler_manifest "${SCHEDULER_MANIFEST_ARG:-${SOURCE_RUN_ROOT}/manifests/scheduler_ids.tsv}"
validate_output_ownership "${OUTPUT_OWNERSHIP_ARG:-${SOURCE_RUN_ROOT}/manifests/output_ownership.tsv}"
require_regular_file "${OWNER_DIR}/owner" "source output owner record"
validate_owner "${OWNERS_MANIFEST_ARG:-${SOURCE_RUN_ROOT}/manifests/owners.tsv}"
parse_prior_inspect "${PRIOR_INSPECT_PATH}"
validate_prior_inspect_identity
validate_prior_attempts
validate_accepted_accounting
validate_sidecar "${OUTPUT_PATH}" "source scratch H5AD"
OUTPUT_SIDECAR_MD5="${VALIDATED_SIDECAR_MD5}"
OUTPUT_SIDECAR_SIZE="${VALIDATED_SIDECAR_SIZE}"
validate_sidecar "${NAS_PATH}" "source NAS H5AD"
validate_artifact_record "${OUTPUT_PATH}" "${ARTIFACT_RECORD_ARG}"
validate_python_h5ad

create_new_run
exit 0
