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
FORCE=0

_builder_die() {
  echo "ERROR: $*" >&2
  exit 1
}

_builder_usage() {
  cat >&2 <<'USAGE'
Usage: build_ecoda_runtime.sh --layout relocated|path-preserving --output ABSOLUTE_SIF [--force]
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
[[ "${OUTPUT}" = /* ]] || _builder_die "--output must be an absolute SIF path"
[[ "${OUTPUT}" != */ ]] || _builder_die "--output must name a SIF file"
[[ "${OUTPUT}" == *.sif ]] || _builder_die "--output must use the .sif suffix"
[[ "${FORCE}" == 1 || ! -e "${OUTPUT}" ]] || _builder_die "output already exists; use --force to replace it: ${OUTPUT}"

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
_ecoda_runtime_realpath_existing "${output_parent}" >/dev/null || exit 1
# Preserve the caller's absolute spelling (notably $HOME/scratch symlinks) so
# the manifest path is identical to slurm_config.sh's runtime default.
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
RUNTIME_SYSTEM_POST_COMMAND="dnf install -y which jq"
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
"${APPTAINER_BIN}" inspect "${temporary_output}" >/dev/null 2>&1 || \
  _builder_die "apptainer inspect failed for the built SIF"

new_lock_sha="$(_ecoda_runtime_sha256 "${PROJECT_ROOT}/pixi.lock")" || exit 1
[[ "${new_lock_sha}" == "${lock_sha}" ]] || _builder_die "pixi.lock changed during image construction"

image_sha="$(_ecoda_runtime_sha256 "${temporary_output}")" || exit 1
umask 077
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

mv -f "${temporary_output}" "${OUTPUT}"
mv -f "${temporary_manifest}" "${manifest}"
trap - EXIT

export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IMAGE="${OUTPUT}"
export ECODA_RUNTIME_MANIFEST="${manifest}"
export ECODA_RUNTIME_PROFILE=default
export ECODA_APPTAINER_NV=0
if ! ecoda_runtime_validate_submission apptainer; then
  rm -f "${manifest}"
  _builder_die "published image failed its immutable runtime contract validation"
fi

[[ -f "${kept_def}" ]] || _builder_die "Pixitainer --keep-def did not retain the generated definition: ${kept_def}"
[[ -s "${OUTPUT}" && -s "${manifest}" ]] || _builder_die "published SIF/manifest pair is incomplete"
echo "Immutable ECODA runtime built: ${OUTPUT}"
echo "Runtime manifest: ${manifest}"
echo "Runtime layout: ${LAYOUT}"
echo "Runtime image SHA-256: ${image_sha}"
echo "Dry-run definition: ${dryrun_def}"
echo "Build log: ${build_log}"
