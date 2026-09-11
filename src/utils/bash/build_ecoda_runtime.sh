#!/bin/bash
# Build the immutable ECODA py-cuda13 Apptainer image on a Bamboo compute node.
# This is a build-time operation only; production workers never invoke Pixi or
# Pixitainer.  Bash 3.2-compatible and intentionally fail-closed.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
source "${SCRIPT_DIR}/env_mutation_lock.sh"

PIXITAINER_VERSION="0.8.3"
BASE_IMAGE="rockylinux:9"
LAYOUT=""
OUTPUT=""
RUNTIME_ID=""
RUNTIME_ID_SET=0
FORCE=0

_builder_die() {
  echo "ERROR: $*" >&2
  exit 1
}

# The configured scratch root may intentionally be a site-provided symlink
# (for example, $HOME/scratch).  Only runtime-owned components below it may
# influence publication, so reject symlinks in that appended ancestry.
_builder_no_symlink_components_below() {
  local root="${1:-}"
  local path="${2:-}"
  local suffix component prefix
  local old_ifs="${IFS}"
  local parts=()
  [[ "${root}" = /* && "${path}" = /* ]] || return 1
  case "${root}" in
    /)
      suffix="${path#/}"
      prefix="/"
      ;;
    *)
      root="${root%/}"
      case "${path}" in
        "${root}"/*)
          suffix="${path#${root}/}"
          prefix="${root}"
          ;;
        *) return 1 ;;
      esac
      ;;
  esac
  IFS='/' read -r -a parts <<< "${suffix}"
  IFS="${old_ifs}"
  for component in "${parts[@]}"; do
    [[ -n "${component}" && "${component}" != "." ]] || continue
    [[ "${component}" != ".." ]] || return 1
    prefix="${prefix%/}/${component}"
    [[ ! -L "${prefix}" ]] || return 1
  done
  return 0
}

_builder_usage() {
  cat >&2 <<'USAGE'
Usage: build_ecoda_runtime.sh --layout relocated|path-preserving --output ABSOLUTE_SIF [--runtime-id RUNTIME_ID] [--force]
USAGE
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --layout)
      [[ $# -ge 2 ]] || _builder_usage
      LAYOUT="$2"
      shift 2
      ;;
    --output)
      [[ $# -ge 2 ]] || _builder_usage
      OUTPUT="$2"
      shift 2
      ;;
    --runtime-id)
      [[ $# -ge 2 ]] || _builder_usage
      RUNTIME_ID="$2"
      RUNTIME_ID_SET=1
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --help|-h)
      _builder_usage
      ;;
    *)
      _builder_usage
      ;;
  esac
done

case "${LAYOUT}" in
  relocated|path-preserving) ;;
  *) _builder_die "--layout must be relocated or path-preserving" ;;
esac

RUNTIME_FORMAT=1
RUNTIME_ROOT=""
RUNTIME_DIR=""
if [[ "${RUNTIME_ID_SET}" == 1 ]]; then
  RUNTIME_FORMAT=2
  [[ -n "${RUNTIME_ID}" ]] || _builder_die "--runtime-id requires a nonempty value"
  [[ "${RUNTIME_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || \
    _builder_die "--runtime-id must be a single safe path component"
  [[ -n "${HPC_SCRATCH_DIR:-}" && "${HPC_SCRATCH_DIR}" = /* ]] || \
    _builder_die "format-2 runtime builds require an absolute HPC_SCRATCH_DIR"
  RUNTIME_ROOT="${HPC_SCRATCH_DIR%/}/_ecoda_runtime"
  RUNTIME_DIR="${RUNTIME_ROOT}/${RUNTIME_ID}"
  if [[ -z "${OUTPUT}" ]]; then
    OUTPUT="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
  fi
  _builder_no_symlink_components_below "${HPC_SCRATCH_DIR}" "${RUNTIME_DIR}" || \
    _builder_die "format-2 runtime path contains a symlinked component: ${RUNTIME_DIR}"
  [[ "$(dirname "${OUTPUT}")" == "${RUNTIME_DIR}" ]] || \
    _builder_die "format-2 output must be directly under ${RUNTIME_DIR}"
  [[ ! -e "${RUNTIME_DIR}" && ! -L "${RUNTIME_DIR}" ]] || \
    _builder_die "format-2 runtime output already exists; choose a new --runtime-id: ${RUNTIME_DIR}"
else
  [[ -n "${OUTPUT}" ]] || _builder_die "--output is required for a legacy FORMAT=1 runtime build"
fi

[[ "${OUTPUT}" = /* ]] || _builder_die "--output must be an absolute SIF path"
[[ "${OUTPUT}" != */ ]] || _builder_die "--output must name a SIF file"
[[ "${OUTPUT}" == *.sif ]] || _builder_die "--output must use the .sif suffix"
if [[ "${RUNTIME_FORMAT}" == 1 ]]; then
  [[ "${FORCE}" == 1 || ! -e "${OUTPUT}" ]] || _builder_die "output already exists; use --force to replace it: ${OUTPUT}"
fi

if [[ "$(uname -s)" != "Linux" ]]; then
  _builder_die "runtime image construction is supported only on a Linux Bamboo compute allocation"
fi
[[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_JOB_NODELIST:-}" ]] || \
  _builder_die "runtime image construction requires a Slurm compute allocation"
[[ -n "${SLURM_JOB_PARTITION:-}" ]] || \
  _builder_die "runtime image construction requires SLURM_JOB_PARTITION"

command -v scontrol >/dev/null 2>&1 || _builder_die "scontrol is required to validate the build allocation"
command -v squeue >/dev/null 2>&1 || _builder_die "squeue is required to validate the build allocation"
job_info="$(scontrol show job "${SLURM_JOB_ID}" -o 2>/dev/null)" || \
  _builder_die "could not query the current Slurm allocation"
[[ -n "${job_info}" ]] || _builder_die "current Slurm allocation query returned no data"
[[ "${job_info}" == *"JobState=RUNNING"* ]] || \
  _builder_die "runtime image build allocation is not RUNNING"
host_short="$(hostname -s 2>/dev/null || hostname)"
case "${host_short}" in
  *login*|*Login*|*LOGIN*) _builder_die "runtime image construction must run on a compute node, not ${host_short}" ;;
esac

[[ -f "${PROJECT_ROOT}/pixi.toml" && -r "${PROJECT_ROOT}/pixi.toml" ]] || \
  _builder_die "pixi.toml is missing or unreadable: ${PROJECT_ROOT}/pixi.toml"
[[ -f "${PROJECT_ROOT}/pixi.lock" && -r "${PROJECT_ROOT}/pixi.lock" ]] || \
  _builder_die "pixi.lock is missing or unreadable: ${PROJECT_ROOT}/pixi.lock"
[[ -d "${PROJECT_ROOT}/.pixi/envs/py-cuda13" ]] || \
  _builder_die "realized py-cuda13 environment is missing: ${PROJECT_ROOT}/.pixi/envs/py-cuda13"
[[ -r "${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python" ]] || \
  _builder_die "realized py-cuda13 Python is missing"
[[ -r "${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/Rscript" ]] || \
  _builder_die "realized py-cuda13 Rscript is missing"
[[ ! -e "${LOGS_DIR}/env_refresh.lock" ]] || \
  _builder_die "environment mutation lock is active: ${LOGS_DIR}/env_refresh.lock"
ENV_LOCK_FILE="${LOGS_DIR}/env_refresh.lock"
ecoda_require_no_active_jobs "${SLURM_JOB_ID}" || exit 1

PIXI_BIN="${PIXI_BIN:-$(command -v pixi || true)}"
[[ -n "${PIXI_BIN}" ]] || _builder_die "pixi is unavailable on the build allocation"
[[ -x "${PIXI_BIN}" ]] || _builder_die "configured pixi binary is not executable: ${PIXI_BIN}"
PIXI_VERSION="$(${PIXI_BIN} -V 2>/dev/null | awk 'NF {print $NF; exit}')" || true
[[ -n "${PIXI_VERSION}" ]] || _builder_die "could not determine the build Pixi version"
"${PIXI_BIN}" containerize --help >/dev/null 2>&1 || \
  _builder_die "pinned Pixitainer extension is unavailable as 'pixi containerize'"

APPTAINER_BIN="${APPTAINER_BIN:-apptainer}"
command -v "${APPTAINER_BIN}" >/dev/null 2>&1 || _builder_die "apptainer is unavailable on the build allocation"
APPTAINER_VERSION="$(${APPTAINER_BIN} --version 2>/dev/null | awk 'NF {print $NF; exit}')" || true
[[ -n "${APPTAINER_VERSION}" ]] || _builder_die "could not determine the Apptainer version"

realized_env="$(_ecoda_runtime_realpath_existing "${PROJECT_ROOT}/.pixi/envs/py-cuda13")" || exit 1
toml_sha="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.toml")" || exit 1
lock_sha="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.lock")" || exit 1
git_revision="$(git -C "${PROJECT_ROOT}" rev-parse HEAD 2>/dev/null)" || \
  _builder_die "could not determine the source Git revision"
[[ -n "${git_revision}" ]] || _builder_die "source Git revision is empty"

case "${LAYOUT}" in
  relocated)
    add_destination="/opt/ecoda/py-cuda13"
    container_prefix="/opt/ecoda/py-cuda13"
    ;;
  path-preserving)
    add_destination="${PROJECT_ROOT}/.pixi/envs/py-cuda13"
    container_prefix="${PROJECT_ROOT}/.pixi/envs/py-cuda13"
    ;;
esac

output_parent="$(dirname "${OUTPUT}")"
mkdir -p "${output_parent}"
if [[ "${RUNTIME_FORMAT}" == 2 ]]; then
  _builder_no_symlink_components_below "${HPC_SCRATCH_DIR}" "${output_parent}" || \
    _builder_die "format-2 output parent contains a symlinked component: ${output_parent}"
fi
_ecoda_runtime_realpath_existing "${output_parent}" >/dev/null || exit 1
if [[ "${RUNTIME_FORMAT}" == 2 ]]; then
  canonical_runtime_root="$(_ecoda_runtime_realpath_existing "${RUNTIME_ROOT}")" || exit 1
  canonical_runtime_dir="$(_ecoda_runtime_realpath_existing "${RUNTIME_DIR}")" || exit 1
  canonical_output_parent="$(_ecoda_runtime_realpath_existing "${output_parent}")" || exit 1
  output_name="$(basename "${OUTPUT}")"
  canonical_output="${canonical_output_parent}/${output_name}"
  [[ "${canonical_runtime_dir}" == "${canonical_runtime_root}/${RUNTIME_ID}" ]] || \
    _builder_die "format-2 runtime directory escaped the configured runtime root: ${RUNTIME_DIR}"
  [[ "${canonical_output_parent}" == "${canonical_runtime_dir}" ]] || \
    _builder_die "format-2 output escaped the configured runtime root: ${OUTPUT}"
  case "${canonical_output}" in
    "${canonical_runtime_root}/${RUNTIME_ID}/"*) ;;
    *) _builder_die "format-2 canonical output escaped the configured runtime root: ${canonical_output}" ;;
  esac
  output_parent="${canonical_output_parent}"
  OUTPUT="${canonical_output}"
fi
# Preserve the caller's absolute spelling for legacy FORMAT=1 so the manifest
# path remains identical to slurm_config.sh's runtime default.
temporary_output="${OUTPUT}.partial.$$"
dryrun_def="${OUTPUT}.dryrun.def"
dryrun_stderr="${OUTPUT}.dryrun.stderr"
build_log="${OUTPUT}.build.log"
manifest="${OUTPUT}.manifest"
temporary_manifest="${manifest}.tmp.$$"
kept_def="${temporary_output%.*}.def"

rm -f "${temporary_output}" "${temporary_manifest}" "${dryrun_def}" "${dryrun_stderr}"
trap 'rm -f "${temporary_output}" "${temporary_manifest}"' EXIT

export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${TMPDIR:-/tmp}/ecoda-apptainer-tmp-${USER:-unknown}}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${TMPDIR:-/tmp}/ecoda-apptainer-cache-${USER:-unknown}}"
mkdir -p "${APPTAINER_TMPDIR}" "${APPTAINER_CACHEDIR}"

cd "${PROJECT_ROOT}"
RUNTIME_COMPAT_ENV_ROOT="${PROJECT_ROOT}/.pixi/envs"
RUNTIME_SYSTEM_POST_COMMAND="dnf install -y which jq diffutils && mkdir -p \"${RUNTIME_COMPAT_ENV_ROOT}\" && ln -s /opt/ecoda/py-cuda13 \"${RUNTIME_COMPAT_ENV_ROOT}/py-cuda13\""
dryrun_args=(
  containerize
  --manual
  --no-install
  --env py-cuda13
  --base-image rockylinux:9
  --pixi-version "${PIXI_VERSION}"
  --add-file "${realized_env}:${add_destination}"
  --post-command "${RUNTIME_SYSTEM_POST_COMMAND}"
  --keep-def
  --dry-run
  --quiet
  --output "${temporary_output}"
)
if ! "${PIXI_BIN}" "${dryrun_args[@]}" > "${dryrun_def}" 2> "${dryrun_stderr}"; then
  _builder_die "Pixitainer dry-run failed; inspect ${dryrun_stderr}"
fi
[[ -s "${dryrun_def}" ]] || _builder_die "Pixitainer dry-run produced no definition: ${dryrun_def}"
grep -Fq "From: ${BASE_IMAGE}" "${dryrun_def}" || \
  _builder_die "dry-run definition does not pin ${BASE_IMAGE}"
grep -Fq "${realized_env}" "${dryrun_def}" || \
  _builder_die "dry-run definition omits the realized py-cuda13 source"
grep -Fq "${add_destination}" "${dryrun_def}" || \
  _builder_die "dry-run definition omits the requested environment destination"
grep -Fq "${RUNTIME_SYSTEM_POST_COMMAND}" "${dryrun_def}" || \
  _builder_die "dry-run definition omits required in-image system utilities"
grep -Fq 'exec "$@"' "${dryrun_def}" || \
  _builder_die "dry-run definition is not a manual direct shell entrypoint"
if grep -Fq 'pixi install' "${dryrun_def}"; then
  _builder_die "dry-run definition attempts a fresh Pixi installation despite --no-install"
fi

build_args=(
  containerize
  --manual
  --no-install
  --env py-cuda13
  --base-image rockylinux:9
  --pixi-version "${PIXI_VERSION}"
  --add-file "${realized_env}:${add_destination}"
  --post-command "${RUNTIME_SYSTEM_POST_COMMAND}"
  --keep-def
  --output "${temporary_output}"
)
if ! "${PIXI_BIN}" "${build_args[@]}" > "${build_log}" 2>&1; then
  cat "${build_log}" >&2 || true
  _builder_die "Pixitainer/Apptainer image build failed; inspect ${build_log}"
fi
[[ -s "${temporary_output}" ]] || _builder_die "image build produced an empty SIF: ${temporary_output}"
new_toml_sha="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.toml")" || exit 1
new_lock_sha="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.lock")" || exit 1
[[ "${new_toml_sha}" == "${toml_sha}" ]] || _builder_die "pixi.toml changed during image construction"
[[ "${new_lock_sha}" == "${lock_sha}" ]] || _builder_die "pixi.lock changed during image construction"

image_sha="$(_ecoda_runtime_sha256 "${temporary_output}")" || exit 1
umask 077
if [[ "${RUNTIME_FORMAT}" == 2 ]]; then
  {
    printf '%s\n' \
      'FORMAT=2' \
      "IMAGE_PATH=${OUTPUT}" \
      "IMAGE_SHA256=${image_sha}" \
      'RUNTIME_ENV=py-cuda13' \
      "RUNTIME_LAYOUT=${LAYOUT}" \
      "CONTAINER_ENV_PREFIX=${container_prefix}" \
      "BASE_IMAGE=${BASE_IMAGE}" \
      "PIXITAINER_VERSION=${PIXITAINER_VERSION}" \
      "PIXI_VERSION=${PIXI_VERSION}" \
      "APPTAINER_VERSION=${APPTAINER_VERSION}" \
      "IMAGE_BUILD_GIT_REVISION=${git_revision}" \
      "IMAGE_PIXI_TOML_SHA256=${toml_sha}" \
      "IMAGE_PIXI_LOCK_SHA256=${lock_sha}"
    if [[ "${LAYOUT}" == "path-preserving" ]]; then
      printf 'CONTAINER_PROJECT_ROOT=%s\n' "${PROJECT_ROOT}"
    fi
  } > "${temporary_manifest}"
else
  {
    printf '%s\n' \
      'FORMAT=1' \
      "IMAGE_PATH=${OUTPUT}" \
      "IMAGE_SHA256=${image_sha}" \
      'RUNTIME_ENV=py-cuda13' \
      "RUNTIME_LAYOUT=${LAYOUT}" \
      "CONTAINER_ENV_PREFIX=${container_prefix}" \
      "BASE_IMAGE=${BASE_IMAGE}" \
      "PIXITAINER_VERSION=${PIXITAINER_VERSION}" \
      "PIXI_VERSION=${PIXI_VERSION}" \
      "APPTAINER_VERSION=${APPTAINER_VERSION}" \
      "GIT_REVISION=${git_revision}" \
      "PIXI_LOCK_SHA256=${lock_sha}"
    if [[ "${LAYOUT}" == "path-preserving" ]]; then
      printf 'CONTAINER_PROJECT_ROOT=%s\n' "${PROJECT_ROOT}"
    fi
  } > "${temporary_manifest}"
fi

mv -f "${temporary_output}" "${OUTPUT}"
mv -f "${temporary_manifest}" "${manifest}"
trap - EXIT

export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IMAGE="${OUTPUT}"
export ECODA_RUNTIME_MANIFEST="${manifest}"
export ECODA_RUNTIME_PROFILE=default
export ECODA_APPTAINER_NV=0
if [[ "${RUNTIME_FORMAT}" == 2 ]]; then
  export ECODA_RUNTIME_BUILD_VALIDATION=1
fi
if ! ecoda_runtime_validate_submission apptainer; then
  rm -f "${manifest}"
  _builder_die "published image failed its immutable runtime contract validation"
fi
unset ECODA_RUNTIME_BUILD_VALIDATION

if [[ "${RUNTIME_FORMAT}" == 2 ]]; then
  chmod a-w "${OUTPUT}" "${manifest}" "${output_parent}" || \
    _builder_die "failed to publish format-2 runtime read-only"
  _ecoda_runtime_require_nonwritable "${OUTPUT}" || exit 1
  _ecoda_runtime_require_nonwritable "${manifest}" || exit 1
  _ecoda_runtime_require_nonwritable "${output_parent}" || exit 1
fi

[[ -f "${kept_def}" ]] || _builder_die "Pixitainer --keep-def did not retain the generated definition: ${kept_def}"
[[ -s "${OUTPUT}" && -s "${manifest}" ]] || _builder_die "published SIF/manifest pair is incomplete"
echo "Immutable ECODA runtime built: ${OUTPUT}"
echo "Runtime manifest: ${manifest}"
echo "Runtime layout: ${LAYOUT}"
echo "Runtime image SHA-256: ${image_sha}"
echo "Dry-run definition: ${dryrun_def}"
echo "Build log: ${build_log}"
