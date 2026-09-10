#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
unset HPC_SCRATCH_DIR HOME_REF_DIR ECODA_RUNTIME_IMAGE ECODA_RUNTIME_MANIFEST

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

assert_eq() {
  local expected="$1"
  local actual="$2"
  local label="${3:-value}"
  [[ "${expected}" == "${actual}" ]] || fail "${label}: expected '${expected}', got '${actual}'"
}

assert_contains() {
  local haystack="$1"
  local needle="$2"
  local label="${3:-value}"
  [[ "${haystack}" == *"${needle}"* ]] || fail "${label}: missing '${needle}'"
}

assert_not_contains() {
  local haystack="$1"
  local needle="$2"
  local label="${3:-value}"
  [[ "${haystack}" != *"${needle}"* ]] || fail "${label}: unexpectedly contains '${needle}'"
}

expect_fail() {
  if "$@"; then
    fail "expected failure: $*"
  fi
}

# ECODA_RUNTIME_MODE=apptainer must not change login-side interpreter paths.
(
  export ECODA_RUNTIME_MODE=apptainer
  export ECODA_RUNTIME_IN_CONTAINER=0
  source "${ROOT}/src/slurm_config.sh"
  assert_eq "${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python" "${PYTHON_BIN}" "host Python path in apptainer mode"
  assert_eq "${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/Rscript --vanilla" "${PIXI_RSCRIPT}" "host R path in apptainer mode"
  assert_not_contains "${PYTHON_BIN}" "/opt/" "host Python path"
)

TEST_TMP_BASE="${TMPDIR:-/tmp}"
TEST_TMP_BASE="${TEST_TMP_BASE%/}"
TEST_ROOT="$(mktemp -d "${TEST_TMP_BASE}/ecoda-runtime-test.XXXXXX")"
TEST_ROOT="$(cd "${TEST_ROOT}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TEST_ROOT}" >/dev/null 2>&1 || true
  rm -rf "${TEST_ROOT}"
}
trap cleanup EXIT

TEST_HOME="${TEST_ROOT}/home"
PROJECT="${TEST_ROOT}/project"
SCRATCH_REAL="${TEST_ROOT}/scratch-real"
SCRATCH_LINK="${TEST_HOME}/scratch/ECODA_paper"
NODE_TMP="${TEST_ROOT}/node-tmp"
REFERENCE="${TEST_HOME}/reference_atlases/sketched_200ct"
FAKE_BIN="${TEST_ROOT}/fake-bin"
IMAGE="${TEST_ROOT}/ecoda-runtime.sif"
MANIFEST="${IMAGE}.manifest"
FAKE_APPTAINER_LOG="${TEST_ROOT}/apptainer.log"
FAKE_SBATCH_LOG="${TEST_ROOT}/sbatch.log"
MARKER="${TEST_ROOT}/worker.marker"
mkdir -p \
  "${PROJECT}/.pixi/envs/py-cuda13/bin" \
  "${PROJECT}/.pixi/envs/py-cuda13/lib/R" \
  "${PROJECT}/src/utils/bash" \
  "${PROJECT}/aux" \
  "${PROJECT}/logs" \
  "${SCRATCH_REAL}" \
  "${TEST_HOME}/scratch" \
  "${NODE_TMP}" \
  "${REFERENCE}" \
  "${FAKE_BIN}"
ln -s "${SCRATCH_REAL}" "${SCRATCH_LINK}"
printf 'lock\n' > "${PROJECT}/pixi.lock"
printf '{}\n' > "${PROJECT}/datasets.json"
printf 'config\n' > "${PROJECT}/config_helper.R"
printf 'python\n' > "${PROJECT}/.pixi/envs/py-cuda13/bin/python"
printf '#!/bin/bash\n' > "${PROJECT}/.pixi/envs/py-cuda13/bin/Rscript"
printf 'sif-bytes\n' > "${IMAGE}"
chmod +x "${PROJECT}/.pixi/envs/py-cuda13/bin/python" "${PROJECT}/.pixi/envs/py-cuda13/bin/Rscript"
cp "${ROOT}/src/slurm_config.sh" "${PROJECT}/src/slurm_config.sh"
cp "${ROOT}/src/utils/bash/ecoda_runtime.sh" "${PROJECT}/src/utils/bash/ecoda_runtime.sh"

cat > "${FAKE_BIN}/apptainer" <<'FAKE_APPTAINER'
#!/bin/bash
set -euo pipefail
: "${FAKE_APPTAINER_LOG:?}"
printf 'argv:' >> "${FAKE_APPTAINER_LOG}"
for arg in "$@"; do printf ' <%s>' "${arg}" >> "${FAKE_APPTAINER_LOG}"; done
printf '\n' >> "${FAKE_APPTAINER_LOG}"
if [[ "${1:-}" == "--version" ]]; then
  echo 'apptainer version 1.3.2'
  exit 0
fi
if [[ "${1:-}" == "inspect" ]]; then
  if [[ "${FAKE_APPTAINER_INSPECT_FAIL:-0}" == 1 ]]; then exit 1; fi
  exit 0
fi
[[ "${1:-}" == "exec" ]] || exit 2
shift
image=''
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env)
      export "$2"
      shift 2
      ;;
    --bind)
      printf 'bind=%s\n' "$2" >> "${FAKE_APPTAINER_LOG}"
      shift 2
      ;;
    --no-mount)
      shift 2
      ;;
    --containall|--no-home|--nv)
      if [[ "$1" == "--nv" ]]; then printf 'nv=1\n' >> "${FAKE_APPTAINER_LOG}"; fi
      shift
      ;;
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
printf 'missing exec command for image %s\n' "${image}" >&2
exit 2
FAKE_APPTAINER
chmod +x "${FAKE_BIN}/apptainer"

cat > "${FAKE_BIN}/git" <<'FAKE_GIT'
#!/bin/bash
set -euo pipefail
if [[ "${1:-}" == "-C" && "${3:-}" == "rev-parse" ]]; then
  printf 'runtime-test-rev\n'
  exit 0
fi
exec /usr/bin/git "$@"
FAKE_GIT
chmod +x "${FAKE_BIN}/git"

cat > "${FAKE_BIN}/sbatch" <<'FAKE_SBATCH'
#!/bin/bash
printf 'sbatch called\n' >> "${FAKE_SBATCH_LOG:?}"
exit 0
FAKE_SBATCH
chmod +x "${FAKE_BIN}/sbatch"

export PATH="${FAKE_BIN}:${PATH}"
export FAKE_APPTAINER_LOG FAKE_SBATCH_LOG
export PROJECT_ROOT="${PROJECT}"
export HPC_SCRATCH_DIR="${SCRATCH_LINK}"
export LOGS_DIR="${PROJECT}/logs"
export HOME_REF_DIR="${REFERENCE}"
export TMPDIR="${NODE_TMP}"
export ECODA_RUNTIME_IMAGE="${IMAGE}"
export ECODA_RUNTIME_MANIFEST="${MANIFEST}"
export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_PROFILE=default
export ECODA_APPTAINER_NV=0
export APPTAINER_BIN="${FAKE_BIN}/apptainer"
export PATH="${FAKE_BIN}:${PROJECT}/.pixi/envs/py-cuda13/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="${PROJECT}/.pixi/envs/py-cuda13/lib:/usr/lib"
export PYTHONHOME="${TEST_ROOT}/host-python"
export PYTHONPATH="${TEST_ROOT}/host-python-path"
export R_LIBS_USER="${TEST_ROOT}/host-r-user"
export R_LIBS_SITE="${TEST_ROOT}/host-r-site"
export R_ENVIRON_USER="${TEST_ROOT}/host-r-env"
export R_PROFILE_USER="${TEST_ROOT}/host-r-profile"
export HOME="${TEST_HOME}"
source "${ROOT}/src/slurm_config.sh"
# slurm_config.sh recomputes PROJECT_ROOT; restore the test contract after sourcing.
export PROJECT_ROOT="${PROJECT}"
export HPC_SCRATCH_DIR="${SCRATCH_LINK}"
export LOGS_DIR="${PROJECT}/logs"
export HOME_REF_DIR="${REFERENCE}"
export TMPDIR="${NODE_TMP}"
export ECODA_RUNTIME_IMAGE="${IMAGE}"
export ECODA_RUNTIME_MANIFEST="${MANIFEST}"
export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_PROFILE=default
export ECODA_APPTAINER_NV=0
export APPTAINER_BIN="${FAKE_BIN}/apptainer"
export PATH="${FAKE_BIN}:${PROJECT}/.pixi/envs/py-cuda13/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="${PROJECT}/.pixi/envs/py-cuda13/lib:/usr/lib"
source "${ROOT}/src/utils/bash/ecoda_runtime.sh"

# FORMAT=2 must reject an existing symlink in the requested runtime ancestry
# before any builder preflight or image publication can follow it.
BUILDER_SYMLINK_PARENT="${TEST_ROOT}/builder-scratch"
BUILDER_ESCAPE="${TEST_ROOT}/builder-escape"
BUILDER_FAKE_BIN="${TEST_ROOT}/builder-fake-bin"
BUILDER_ERR="${TEST_ROOT}/builder-symlink.stderr"
BUILDER_PIXI_LOG="${TEST_ROOT}/builder-pixi.log"
mkdir -p "${BUILDER_SYMLINK_PARENT}" "${BUILDER_ESCAPE}" "${BUILDER_FAKE_BIN}"
ln -s "${BUILDER_ESCAPE}" "${BUILDER_SYMLINK_PARENT}/_ecoda_runtime"
printf 'builder toml\n' > "${PROJECT}/pixi.toml"
cp "${ROOT}/src/utils/bash/build_ecoda_runtime.sh" \
  "${PROJECT}/src/utils/bash/build_ecoda_runtime.sh"
cp "${ROOT}/src/utils/bash/env_mutation_lock.sh" \
  "${PROJECT}/src/utils/bash/env_mutation_lock.sh"

cat > "${BUILDER_FAKE_BIN}/uname" <<'BUILDER_UNAME'
#!/bin/bash
printf 'Linux\n'
BUILDER_UNAME
cat > "${BUILDER_FAKE_BIN}/scontrol" <<'BUILDER_SCONTROL'
#!/bin/bash
printf 'JobState=RUNNING\n'
BUILDER_SCONTROL
cat > "${BUILDER_FAKE_BIN}/squeue" <<'BUILDER_SQUEUE'
#!/bin/bash
exit 0
BUILDER_SQUEUE
cat > "${BUILDER_FAKE_BIN}/hostname" <<'BUILDER_HOSTNAME'
#!/bin/bash
printf 'builder-node\n'
BUILDER_HOSTNAME
cat > "${BUILDER_FAKE_BIN}/git" <<'BUILDER_GIT'
#!/bin/bash
printf 'builder-test-revision\n'
BUILDER_GIT
cat > "${BUILDER_FAKE_BIN}/apptainer" <<'BUILDER_APPTAINER'
#!/bin/bash
printf 'apptainer version 1.3.2\n'
BUILDER_APPTAINER
cat > "${BUILDER_FAKE_BIN}/pixi" <<'BUILDER_PIXI'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${BUILDER_PIXI_LOG}"
if [[ "${1:-}" == "-V" ]]; then
  printf 'pixi 0.49.0\n'
  exit 0
fi
if [[ "${1:-}" == "containerize" && "${2:-}" == "--help" ]]; then
  exit 0
fi
output=""
dryrun=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output)
      output="$2"
      shift 2
      ;;
    --dry-run)
      dryrun=1
      shift
      ;;
    *)
      shift
      ;;
  esac
done
[[ -n "${output}" ]] || exit 1
if [[ "${dryrun}" == 1 ]]; then
  printf '%s\n' \
    'From: rockylinux:9' \
    "$(realpath "${BUILDER_PROJECT}/.pixi/envs/py-cuda13")" \
    '/opt/ecoda/py-cuda13' \
    'dnf install -y which jq' \
    'exec "$@"'
else
  printf 'sif-bytes\n' > "${output}"
  printf 'kept-definition\n' > "${output%.*}.def"
fi
BUILDER_PIXI
chmod +x \
  "${BUILDER_FAKE_BIN}/uname" \
  "${BUILDER_FAKE_BIN}/scontrol" \
  "${BUILDER_FAKE_BIN}/squeue" \
  "${BUILDER_FAKE_BIN}/hostname" \
  "${BUILDER_FAKE_BIN}/git" \
  "${BUILDER_FAKE_BIN}/apptainer" \
  "${BUILDER_FAKE_BIN}/pixi"
export BUILDER_PROJECT="${PROJECT}"
export BUILDER_PIXI_LOG
if (
  export HPC_SCRATCH_DIR="${BUILDER_SYMLINK_PARENT}"
  export PROJECT_ROOT="${PROJECT}"
  export LOGS_DIR="${PROJECT}/logs"
  export TMPDIR="${NODE_TMP}"
  export SLURM_JOB_ID=builder-test-job
  export SLURM_JOB_NODELIST=builder-test-node
  export SLURM_JOB_PARTITION=shared-cpu
  export PIXI_BIN="${BUILDER_FAKE_BIN}/pixi"
  export APPTAINER_BIN="${BUILDER_FAKE_BIN}/apptainer"
  export PATH="${BUILDER_FAKE_BIN}:${PATH}"
  bash "${PROJECT}/src/utils/bash/build_ecoda_runtime.sh" \
    --layout relocated --runtime-id hostile
) > /dev/null 2> "${BUILDER_ERR}"; then
  fail "builder accepted a symlinked FORMAT=2 runtime parent"
fi
assert_contains "$(cat "${BUILDER_ERR}")" \
  'format-2 runtime path contains a symlinked component' \
  'symlinked FORMAT=2 runtime parent rejection'
[[ ! -e "${BUILDER_ESCAPE}/hostile" ]] || \
  fail "builder created a runtime directory outside the requested root"
[[ ! -e "${BUILDER_ESCAPE}/ecoda-py-cuda13.sif" ]] || \
  fail "builder published an image outside the requested root"
[[ ! -e "${BUILDER_PIXI_LOG}" ]] || \
  fail "builder reached Pixi before rejecting the symlinked runtime parent"

# A missing FORMAT=2 runtime root remains creatable when its ancestry has no
# symlinks below the trusted configured scratch alias; the fake build
# exercises canonical publication and immutable parent setup.
BUILDER_VALID_ROOT="${SCRATCH_LINK}"
BUILDER_VALID_IMAGE="${SCRATCH_REAL}/_ecoda_runtime/valid/ecoda-py-cuda13.sif"
if (
  export HPC_SCRATCH_DIR="${BUILDER_VALID_ROOT}"
  export PROJECT_ROOT="${PROJECT}"
  export LOGS_DIR="${PROJECT}/logs"
  export TMPDIR="${NODE_TMP}"
  export SLURM_JOB_ID=builder-valid-job
  export SLURM_JOB_NODELIST=builder-valid-node
  export SLURM_JOB_PARTITION=shared-cpu
  export PIXI_BIN="${BUILDER_FAKE_BIN}/pixi"
  export APPTAINER_BIN="${BUILDER_FAKE_BIN}/apptainer"
  export PATH="${BUILDER_FAKE_BIN}:${PATH}"
  bash "${PROJECT}/src/utils/bash/build_ecoda_runtime.sh" \
    --layout relocated --runtime-id valid
) > /dev/null 2> "${TEST_ROOT}/builder-valid.stderr"; then
  :
else
  cat "${TEST_ROOT}/builder-valid.stderr" >&2
  fail "builder rejected a valid missing FORMAT=2 runtime root"
fi
[[ -s "${BUILDER_VALID_IMAGE}" ]] || \
  fail "builder did not publish a valid FORMAT=2 image"
assert_contains "$(cat "${BUILDER_VALID_IMAGE}.manifest")" \
  'FORMAT=2' 'valid FORMAT=2 manifest'
[[ ! -w "${BUILDER_VALID_ROOT}/_ecoda_runtime/valid" ]] || \
  fail "builder did not publish an immutable FORMAT=2 parent"

# Sourcing the config after the Apptainer boundary must preserve inherited
# absolute destinations rather than recomputing them from the container HOME.
container_paths="$(
  HOME="${TEST_HOME}/container-home" PATH="/usr/bin:/bin" \
  ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_IN_CONTAINER=1 \
  ECODA_RUNTIME_PREFIX="${PROJECT}/.pixi/envs/py-cuda13" \
  HPC_SCRATCH_DIR="/host/scratch/ECODA_paper" \
  HOME_REF_DIR="/host/reference_atlases/sketched_200ct" \
  ECODA_RUNTIME_IMAGE="/host/scratch/ECODA_paper/_ecoda_runtime/host-image.sif" \
  bash -c 'source "$1/src/slurm_config.sh"; printf "%s\n%s\n%s" \
    "$HPC_SCRATCH_DIR" "$HOME_REF_DIR" "$ECODA_RUNTIME_IMAGE"' \
    _ "${PROJECT}"
)"
expected_container_paths="$(printf '%s\n' \
  "/host/scratch/ECODA_paper" \
  "/host/reference_atlases/sketched_200ct" \
  "/host/scratch/ECODA_paper/_ecoda_runtime/host-image.sif")"
assert_eq "${expected_container_paths}" "${container_paths}" \
  "inherited container path contracts"
[[ -d "${PROJECT}" ]] || fail "test project vanished before runtime validation: ${PROJECT}"
realpath "${PROJECT}" >/dev/null || fail "external realpath failed: ${PROJECT}"

image_sha="$(_ecoda_runtime_sha256 "${IMAGE}")"
lock_sha="$(_ecoda_runtime_sha256 "${PROJECT}/pixi.lock")"

write_manifest() {
  local layout="$1"
  local prefix="$2"
  {
    printf '%s\n' \
      'FORMAT=1' \
      "IMAGE_PATH=${IMAGE}" \
      "IMAGE_SHA256=${image_sha}" \
      'RUNTIME_ENV=py-cuda13' \
      "RUNTIME_LAYOUT=${layout}" \
      "CONTAINER_ENV_PREFIX=${prefix}" \
      'BASE_IMAGE=rockylinux:9' \
      'PIXITAINER_VERSION=0.8.3' \
      'PIXI_VERSION=0.49.0' \
      'APPTAINER_VERSION=1.3.2' \
      'GIT_REVISION=runtime-test-rev' \
      "PIXI_LOCK_SHA256=${lock_sha}"
    if [[ "${layout}" == path-preserving ]]; then
      printf 'CONTAINER_PROJECT_ROOT=%s\n' "${PROJECT}"
    fi
  } > "${MANIFEST}"
}

# Relocated layout validates and keeps the host-side paths untouched.
write_manifest relocated /opt/ecoda/py-cuda13
ecoda_runtime_validate_submission apptainer
assert_eq relocated "${ECODA_RUNTIME_LAYOUT}" "relocated layout"
assert_eq /opt/ecoda/py-cuda13 "${ECODA_RUNTIME_CONTAINER_PREFIX}" "relocated prefix"
relocated_binds="$(printf '%s\n' "${ECODA_RUNTIME_BIND_ARGS[@]}")"
project_real="$(realpath "${PROJECT}")"
scratch_real="$(realpath "${SCRATCH_LINK}")"
logs_real="$(realpath "${PROJECT}/logs")"
tmp_real="$(realpath "${NODE_TMP}")"
reference_real="$(realpath "${REFERENCE}")"
assert_contains "${relocated_binds}" "${project_real}:${PROJECT}:ro" "relocated project bind"
assert_contains "${relocated_binds}" "${scratch_real}:${SCRATCH_LINK}:rw" "resolved scratch bind"
assert_contains "${relocated_binds}" "${logs_real}:${PROJECT}/logs:rw" "logs bind"
assert_contains "${relocated_binds}" "${tmp_real}:${NODE_TMP}:rw" "node-local bind"

runtime_export="$(ecoda_runtime_export_csv stage4 1)"
assert_eq "ECODA_RUNTIME_MODE=apptainer,ECODA_RUNTIME_IMAGE=${IMAGE},ECODA_RUNTIME_MANIFEST=${MANIFEST},ECODA_RUNTIME_PROFILE=stage4,ECODA_APPTAINER_NV=1" "${runtime_export}" "runtime export"
expect_fail ecoda_runtime_export_csv stage4 2

# Path-preserving layout requires the exact source root and explicit mounts.
write_manifest path-preserving "${PROJECT}/.pixi/envs/py-cuda13"
export ECODA_RUNTIME_PROFILE=stage4
ecoda_runtime_validate_submission apptainer
path_binds="$(printf '%s\n' "${ECODA_RUNTIME_BIND_ARGS[@]}")"
assert_contains "${path_binds}" "${project_real}/src:${PROJECT}/src:ro" "source bind"
assert_contains "${path_binds}" "${project_real}/datasets.json:${PROJECT}/datasets.json:ro" "datasets bind"
assert_contains "${path_binds}" "${project_real}/config_helper.R:${PROJECT}/config_helper.R:ro" "config bind"
assert_contains "${path_binds}" "${project_real}/aux:${PROJECT}/aux:ro" "aux bind"
assert_contains "${path_binds}" "${reference_real}:${REFERENCE}:ro" "reference bind"
assert_not_contains "${path_binds}" "${PROJECT}:${PROJECT}:ro" "path-preserving project root bind"

# Invalid image, manifest, identity, inspect, and bind inputs fail before any worker marker or scheduler submission.
cp "${IMAGE}" "${IMAGE}.saved"
: > "${IMAGE}"
expect_fail ecoda_runtime_validate_submission apptainer
cp "${IMAGE}.saved" "${IMAGE}"
rm -f "${MANIFEST}"
expect_fail ecoda_runtime_validate_submission apptainer
write_manifest path-preserving "${PROJECT}/.pixi/envs/py-cuda13"
printf 'FORMAT=1\n' > "${MANIFEST}"
expect_fail ecoda_runtime_validate_submission apptainer
write_manifest path-preserving "${PROJECT}/.pixi/envs/py-cuda13"
awk 'BEGIN { done=0 } /^IMAGE_SHA256=/ && !done { print "IMAGE_SHA256=" sprintf("%064d", 0); done=1; next } { print }' "${MANIFEST}" > "${MANIFEST}.bad"
mv -f "${MANIFEST}.bad" "${MANIFEST}"
expect_fail ecoda_runtime_validate_submission apptainer
write_manifest path-preserving "${PROJECT}/.pixi/envs/py-cuda13"
export FAKE_APPTAINER_INSPECT_FAIL=1
expect_fail ecoda_runtime_validate_submission apptainer
unset FAKE_APPTAINER_INSPECT_FAIL
write_manifest path-preserving "${PROJECT}/.pixi/envs/py-cuda13"

ECODA_RUNTIME_BIND_ARGS=()
ECODA_RUNTIME_BIND_DESTS=""
ECODA_RUNTIME_CONTAINER_PREFIX=/opt/ecoda/py-cuda13
export ECODA_RUNTIME_CONTAINER_PREFIX
expect_fail _ecoda_runtime_add_bind "${PROJECT}" "/opt/ecoda/py-cuda13/hidden" ro
TMPDIR="${TEST_ROOT}/missing-tmp"
expect_fail ecoda_runtime_build_bind_args stage4
export TMPDIR="${NODE_TMP}"
[[ ! -e "${MARKER}" ]] || fail "invalid cases created a worker marker"
[[ ! -e "${FAKE_SBATCH_LOG}" ]] || fail "invalid cases submitted a scheduler job"

# A scientific worker crosses exactly one boundary; the inner worker sees direct
# image paths and the existing R-to-Python contract without nested Apptainer.
cat > "${PROJECT}/src/runtime_boundary_worker.sh" <<'WORKER'
#!/bin/bash
set -euo pipefail
source "${PROJECT_ROOT}/src/slurm_config.sh"
source "${PROJECT_ROOT}/src/utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage4 "${PROJECT_ROOT}/src/runtime_boundary_worker.sh"
[[ "${ECODA_RUNTIME_IN_CONTAINER}" == 1 ]] || { echo 'inner runtime guard missing' >&2; exit 1; }
[[ "${PYTHON_BIN}" == "${ECODA_RUNTIME_PREFIX}/bin/python" ]] || { echo "PYTHON_BIN=${PYTHON_BIN} prefix=${ECODA_RUNTIME_PREFIX}" >&2; exit 1; }
[[ "${PIXI_RSCRIPT}" == "${ECODA_RUNTIME_PREFIX}/bin/Rscript --vanilla" ]] || { echo "PIXI_RSCRIPT=${PIXI_RSCRIPT} prefix=${ECODA_RUNTIME_PREFIX}" >&2; exit 1; }
[[ "${R_HOME}" == "${ECODA_RUNTIME_PREFIX}/lib/R" ]] || { echo "R_HOME=${R_HOME} prefix=${ECODA_RUNTIME_PREFIX}" >&2; exit 1; }
[[ "${RETICULATE_PYTHON}" == "${ECODA_RUNTIME_PREFIX}/bin/python" ]] || { echo "RETICULATE_PYTHON=${RETICULATE_PYTHON} prefix=${ECODA_RUNTIME_PREFIX}" >&2; exit 1; }
[[ "${PYTHONNOUSERSITE}" == 1 ]] || { echo "PYTHONNOUSERSITE=${PYTHONNOUSERSITE}" >&2; exit 1; }
[[ "${SCGATE_MODEL_CACHE_DIR}" == "${HOME_REF_DIR}/scGate_models" ]] || { echo "SCGATE_MODEL_CACHE_DIR=${SCGATE_MODEL_CACHE_DIR}" >&2; exit 1; }
[[ "${SCGATE_ONTOLOGY_BRANCH}" == "master" ]] || { echo "SCGATE_ONTOLOGY_BRANCH=${SCGATE_ONTOLOGY_BRANCH}" >&2; exit 1; }
[[ -z "${PYTHONHOME:-}" && -z "${PYTHONPATH:-}" ]] || { echo "host Python import variables survived" >&2; exit 1; }
[[ -z "${R_LIBS_USER:-}" && -z "${R_LIBS_SITE:-}" ]] || { echo "host R library variables survived" >&2; exit 1; }
[[ "${PATH}" == "${ECODA_RUNTIME_PREFIX}/bin:"* ]] || { echo "PATH=${PATH}" >&2; exit 1; }
[[ "${LD_LIBRARY_PATH}" == "${ECODA_RUNTIME_PREFIX}/lib:"* ]] || { echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >&2; exit 1; }
if [[ "${ECODA_RUNTIME_PREFIX}" == /opt/* ]]; then
  [[ "${PATH}" != *"${ECODA_HOST_ENV_PREFIX}"* ]] || { echo "host runtime path survived PATH=${PATH}" >&2; exit 1; }
  [[ "${LD_LIBRARY_PATH}" != *"${ECODA_HOST_ENV_PREFIX}"* ]] || { echo "host runtime path survived LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >&2; exit 1; }
fi
touch "${WORKER_MARKER}"
WORKER
chmod +x "${PROJECT}/src/runtime_boundary_worker.sh"
export WORKER_MARKER="${MARKER}"
: > "${FAKE_APPTAINER_LOG}"
export ECODA_RUNTIME_PROFILE=stage4
export ECODA_APPTAINER_NV=1
bash "${PROJECT}/src/runtime_boundary_worker.sh"
exec_count="$(awk '/argv: <exec>/{count++} END {print count + 0}' "${FAKE_APPTAINER_LOG}")"
assert_eq 1 "${exec_count}" "container boundary count with GPU passthrough"
assert_contains "$(cat "${FAKE_APPTAINER_LOG}")" 'nv=1' 'GPU passthrough'
[[ -f "${MARKER}" ]] || fail "inner worker did not run"
rm -f "${MARKER}"
: > "${FAKE_APPTAINER_LOG}"
export ECODA_APPTAINER_NV=0
bash "${PROJECT}/src/runtime_boundary_worker.sh"
assert_eq 1 "$(awk '/argv: <exec>/{count++} END {print count + 0}' "${FAKE_APPTAINER_LOG}")" "container boundary count without GPU passthrough"
assert_not_contains "$(cat "${FAKE_APPTAINER_LOG}")" 'nv=1' 'CPU passthrough'
[[ -f "${MARKER}" ]] || fail "inner CPU worker did not run"

# FORMAT=2 uses a commit-keyed immutable source snapshot while the runtime image
# records the build checkout independently.  Keep this fixture separate from
# the legacy path-preserving checks above so both contracts remain executable.
F2_SNAPSHOT_PARENT="${SCRATCH_REAL}/_ecoda_source_snapshots"
F2_SOURCE_COMMIT="bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
F2_SNAPSHOT_ROOT="${F2_SNAPSHOT_PARENT}/${F2_SOURCE_COMMIT}"
F2_SOURCE_ROOT="${F2_SNAPSHOT_ROOT}/tree"
F2_SOURCE_IDENTITY="${F2_SNAPSHOT_ROOT}/identity"
F2_SOURCE_MANIFEST="${F2_SOURCE_IDENTITY}/source.manifest"
F2_SOURCE_ARCHIVE="${F2_SOURCE_IDENTITY}/source.tar"
F2_RUNTIME_COLLECTION="${SCRATCH_REAL}/_ecoda_runtime"
F2_RUNTIME_DIR="${F2_RUNTIME_COLLECTION}/runtime-A"
F2_IMAGE="${F2_RUNTIME_DIR}/ecoda-py-cuda13.sif"
F2_RUNTIME_MANIFEST="${F2_IMAGE}.manifest"
F2_RUN_ROOT="${TEST_ROOT}/f2-run"
F2_LOGS="${TEST_ROOT}/f2-logs"
F2_SCRATCH="${TEST_ROOT}/f2-scratch"
F2_EXEC_MARKER="${TEST_ROOT}/f2-snapshot.marker"
F2_HASH_LOG="${TEST_ROOT}/f2-image-hash.log"

mkdir -p \
  "${F2_SOURCE_ROOT}/src/utils/bash" \
  "${F2_SOURCE_ROOT}/aux" \
  "${F2_SOURCE_IDENTITY}" \
  "${F2_RUNTIME_DIR}" \
  "${F2_RUN_ROOT}/manifests" \
  "${F2_LOGS}" \
  "${F2_SCRATCH}"
cp "${ROOT}/src/slurm_config.sh" "${F2_SOURCE_ROOT}/src/slurm_config.sh"
cp "${ROOT}/src/utils/bash/ecoda_runtime.sh" \
  "${F2_SOURCE_ROOT}/src/utils/bash/ecoda_runtime.sh"
printf 'format-2 lock\n' > "${F2_SOURCE_ROOT}/pixi.lock"
printf 'format-2 toml\n' > "${F2_SOURCE_ROOT}/pixi.toml"
printf '{"format":2}\n' > "${F2_SOURCE_ROOT}/datasets.json"
printf 'format-2 config\n' > "${F2_SOURCE_ROOT}/config_helper.R"
printf 'f2 scGateDB\n' > "${F2_SOURCE_ROOT}/aux/scGateDB.rds"
printf 'f2 blocklist\n' > "${F2_SOURCE_ROOT}/aux/genes.blocklist.rds"
printf 'f2 ensembl map\n' > \
  "${F2_SOURCE_ROOT}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
cat > "${F2_SOURCE_ROOT}/src/snapshot_worker.sh" <<F2_SNAPSHOT_WORKER
#!/bin/bash
set -euo pipefail
[[ "\${ECODA_SOURCE_ROOT}" == "${F2_SOURCE_ROOT}" ]] || exit 31
[[ "\${ECODA_SOURCE_MANIFEST}" == "${F2_SOURCE_MANIFEST}" ]] || exit 32
[[ "\${ECODA_SOURCE_SNAPSHOT_REQUIRED}" == 1 ]] || exit 33
[[ "\${ECODA_AUX_ROOT}" == "${F2_SOURCE_ROOT}/aux" ]] || exit 34
[[ "\${ECODA_LOGS_DIR}" == "${F2_LOGS}" ]] || exit 35
[[ "\${PYTHONDONTWRITEBYTECODE}" == 1 ]] || exit 36
[[ -f "\${ECODA_AUX_ROOT}/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz" ]] || exit 37
if mkdir "\${ECODA_SOURCE_ROOT}/__pycache__" 2>/dev/null; then
  exit 38
fi
[[ ! -e "\${ECODA_SOURCE_ROOT}/__pycache__" ]] || exit 39
printf 'snapshot-A\n' > "${F2_EXEC_MARKER}"
F2_SNAPSHOT_WORKER
chmod +x "${F2_SOURCE_ROOT}/src/snapshot_worker.sh"
tar -cf "${F2_SOURCE_ARCHIVE}" -C "${F2_SOURCE_ROOT}" .
F2_ARCHIVE_SHA="$(_ecoda_runtime_sha256 "${F2_SOURCE_ARCHIVE}")"
F2_CONFIG_SHA="$(_ecoda_runtime_sha256 "${F2_SOURCE_ROOT}/config_helper.R")"
F2_DATASETS_SHA="$(_ecoda_runtime_sha256 "${F2_SOURCE_ROOT}/datasets.json")"
F2_TOML_SHA="$(_ecoda_runtime_sha256 "${F2_SOURCE_ROOT}/pixi.toml")"
F2_LOCK_SHA="$(_ecoda_runtime_sha256 "${F2_SOURCE_ROOT}/pixi.lock")"
F2_SCGATE_BRANCH="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4"
{
  printf '%s\n' \
    'FORMAT=1' \
    "SOURCE_ROOT=${F2_SOURCE_ROOT}" \
    "SOURCE_COMMIT=${F2_SOURCE_COMMIT}" \
    "SOURCE_ARCHIVE_PATH=${F2_SOURCE_ARCHIVE}" \
    "SOURCE_ARCHIVE_SHA256=${F2_ARCHIVE_SHA}" \
    "CONFIG_HELPER_SHA256=${F2_CONFIG_SHA}" \
    "DATASETS_SHA256=${F2_DATASETS_SHA}" \
    "PIXI_TOML_SHA256=${F2_TOML_SHA}" \
    "PIXI_LOCK_SHA256=${F2_LOCK_SHA}" \
    "AUX_ROOT=${F2_SOURCE_ROOT}/aux" \
    "SCGATE_DB_BRANCH=${F2_SCGATE_BRANCH}"
} > "${F2_SOURCE_MANIFEST}"
printf 'COMPLETE\n' > "${F2_SNAPSHOT_ROOT}/COMPLETE"
chmod -R a-w "${F2_SNAPSHOT_ROOT}"
cp "${F2_SOURCE_MANIFEST}" "${F2_RUN_ROOT}/manifests/source.manifest"

printf 'format-2 image bytes\n' > "${F2_IMAGE}"
F2_IMAGE_SHA="$(_ecoda_runtime_sha256 "${F2_IMAGE}")"
{
  printf '%s\n' \
    'FORMAT=2' \
    "IMAGE_PATH=${F2_IMAGE}" \
    "IMAGE_SHA256=${F2_IMAGE_SHA}" \
    'RUNTIME_ENV=py-cuda13' \
    'RUNTIME_LAYOUT=relocated' \
    'CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13' \
    'BASE_IMAGE=rockylinux:9' \
    'PIXITAINER_VERSION=0.8.3' \
    'PIXI_VERSION=0.49.0' \
    'APPTAINER_VERSION=1.3.2' \
    'IMAGE_BUILD_GIT_REVISION=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa' \
    "IMAGE_PIXI_TOML_SHA256=${F2_TOML_SHA}" \
    "IMAGE_PIXI_LOCK_SHA256=${F2_LOCK_SHA}"
} > "${F2_RUNTIME_MANIFEST}"
chmod -R a-w "${F2_RUNTIME_DIR}"

export PROJECT_ROOT="${PROJECT}"
export ECODA_SOURCE_ROOT="${F2_SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${F2_SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_AUX_ROOT="${F2_SOURCE_ROOT}/aux"
export ECODA_RUN_ROOT="${F2_RUN_ROOT}"
export ECODA_RUN_ID=f2-runtime-A
export HPC_SCRATCH_DIR="${SCRATCH_REAL}"
export ECODA_SCRATCH_ROOT="${F2_SCRATCH}"
export ECODA_LOGS_DIR="${F2_LOGS}"
export HOME_REF_DIR="${REFERENCE}"
export ECODA_HOST_ENV_PREFIX="${PROJECT}/.pixi/envs/py-cuda13"
export ECODA_HOST_PYTHON_BIN="${ECODA_HOST_ENV_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${ECODA_HOST_ENV_PREFIX}/bin/Rscript --vanilla"
export ECODA_RUNTIME_IMAGE="${F2_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${F2_RUNTIME_MANIFEST}"
export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_PROFILE=stage4
export ECODA_APPTAINER_NV=0
unset ECODA_RUNTIME_BUILD_VALIDATION || true

# The image build revision (A) intentionally differs from source commit (B).
# Matching dependency fields make this a valid independent source/runtime pair.
ecoda_runtime_validate_submission apptainer
assert_eq 2 "${ECODA_RUNTIME_FORMAT}" "format-2 runtime format"
assert_eq relocated "${ECODA_RUNTIME_LAYOUT}" "format-2 relocated layout"
assert_contains "$(cat "${F2_RUN_ROOT}/manifests/runtime.identity")" \
  "RUNTIME_IMAGE=${F2_IMAGE}" "run-bound image identity"
assert_contains "$(cat "${F2_RUN_ROOT}/manifests/runtime.identity")" \
  "IMAGE_PIXI_LOCK_SHA256=${F2_LOCK_SHA}" "run-bound dependency identity"
assert_not_contains "$(cat "${F2_RUN_ROOT}/manifests/runtime.identity")" \
  "st_dev" "run identity device independence"
f2_binds="$(printf '%s\n' "${ECODA_RUNTIME_BIND_ARGS[@]}")"
assert_contains "${f2_binds}" "${F2_SOURCE_ROOT}:${F2_SOURCE_ROOT}:ro" \
  "snapshot source bind"
assert_contains "${f2_binds}" "${F2_SOURCE_ROOT}/aux:${F2_SOURCE_ROOT}/aux:ro" \
  "snapshot auxiliary bind"
assert_contains "${f2_binds}" "${F2_LOGS}:${F2_LOGS}:rw" \
  "separate run logs bind"
assert_contains "${f2_binds}" "${F2_SCRATCH}:${F2_SCRATCH}:rw" \
  "separate scratch bind"
assert_not_contains "${f2_binds}" "${PROJECT}/src:${PROJECT}/src:ro" \
  "mutable checkout source bind"
F2_HOST_PYTHON_SHA="$(_ecoda_runtime_sha256 "${ECODA_HOST_ENV_PREFIX}/bin/python")"
F2_HOST_RSCRIPT_SHA="$(_ecoda_runtime_sha256 "${ECODA_HOST_ENV_PREFIX}/bin/Rscript")"
f2_export="$(ecoda_runtime_export_csv stage4 0)"
assert_contains "${f2_export}" "ECODA_SOURCE_ROOT=${F2_SOURCE_ROOT}" \
  "format-2 source export"
assert_contains "${f2_export}" "ECODA_RUNTIME_IMAGE_SHA256=${F2_IMAGE_SHA}" \
  "format-2 image digest export"
assert_contains "${f2_export}" "ECODA_HOST_PYTHON_SHA256=${F2_HOST_PYTHON_SHA}" \
  "format-2 host Python digest export"
assert_contains "${f2_export}" "ECODA_HOST_RSCRIPT_SHA256=${F2_HOST_RSCRIPT_SHA}" \
  "format-2 host Rscript digest export"
export ECODA_HOST_PYTHON_SHA256="${F2_HOST_PYTHON_SHA}"
export ECODA_HOST_RSCRIPT_SHA256="${F2_HOST_RSCRIPT_SHA}"

# Bound validation rechecks small identities and source/archive content, but it
# must not hash the large image a second time.
: > "${F2_HASH_LOG}"
(
  source "${ROOT}/src/utils/bash/ecoda_runtime.sh"
  _ecoda_runtime_sha256() {
    if [[ "${1:-}" == "${F2_IMAGE}" ]]; then
      printf 'large-image-hash\n' >> "${F2_HASH_LOG}"
    fi
    if command -v sha256sum >/dev/null 2>&1; then
      sha256sum "${1}" | awk '{print $1}'
    else
      shasum -a 256 "${1}" | awk '{print $1}'
    fi
  }
  ecoda_runtime_validate_bound_run
)
assert_eq 0 "$(wc -l < "${F2_HASH_LOG}" | tr -d '[:space:]')" \
  "bound validation does not rehash image"

f2_replace_key() {
  local file="$1"
  local key="$2"
  local value="$3"
  local temporary="${file}.tmp.$$"
  awk -v wanted="${key}" -v replacement="${value}" '
    BEGIN { replaced = 0 }
    index($0, wanted "=") == 1 && replaced == 0 {
      print wanted "=" replacement
      replaced = 1
      next
    }
    { print }
    END { exit(replaced == 1 ? 0 : 1) }
  ' "${file}" > "${temporary}"
  mv -f "${temporary}" "${file}"
}

# Runtime dependency identity must match the frozen source lock/config.
cp "${F2_RUNTIME_MANIFEST}" "${TEST_ROOT}/f2-runtime.manifest.saved"
chmod u+w "${F2_RUNTIME_DIR}"
f2_replace_key "${F2_RUNTIME_MANIFEST}" IMAGE_PIXI_LOCK_SHA256 \
  "0000000000000000000000000000000000000000000000000000000000000000"
chmod a-w "${F2_RUNTIME_MANIFEST}" "${F2_RUNTIME_DIR}"
expect_fail ecoda_runtime_validate_submission apptainer
chmod u+w "${F2_RUNTIME_DIR}" "${F2_RUNTIME_MANIFEST}"
cp "${TEST_ROOT}/f2-runtime.manifest.saved" "${F2_RUNTIME_MANIFEST}"
chmod a-w "${F2_RUNTIME_MANIFEST}" "${F2_RUNTIME_DIR}"

# A source lock/config digest mutation is rejected even when the run copy is
# updated consistently; the retained archive and source tree remain trusted
# only when their recorded digests agree.
cp "${F2_SOURCE_MANIFEST}" "${TEST_ROOT}/f2-source.manifest.saved"
cp "${F2_RUN_ROOT}/manifests/source.manifest" \
  "${TEST_ROOT}/f2-run-source.manifest.saved"
for f2_source_key in CONFIG_HELPER_SHA256 PIXI_LOCK_SHA256; do
  chmod u+w "${F2_SOURCE_IDENTITY}"
  f2_replace_key "${F2_SOURCE_MANIFEST}" "${f2_source_key}" \
    "0000000000000000000000000000000000000000000000000000000000000000"
  chmod u+w "${F2_RUN_ROOT}/manifests"
  f2_replace_key "${F2_RUN_ROOT}/manifests/source.manifest" "${f2_source_key}" \
    "0000000000000000000000000000000000000000000000000000000000000000"
  chmod a-w "${F2_SOURCE_MANIFEST}" "${F2_SOURCE_IDENTITY}"
  chmod a-w "${F2_RUN_ROOT}/manifests/source.manifest"
  expect_fail ecoda_runtime_validate_submission apptainer
  chmod u+w "${F2_SOURCE_IDENTITY}" "${F2_SOURCE_MANIFEST}"
  cp "${TEST_ROOT}/f2-source.manifest.saved" "${F2_SOURCE_MANIFEST}"
  chmod a-w "${F2_SOURCE_MANIFEST}" "${F2_SOURCE_IDENTITY}"
  chmod u+w "${F2_RUN_ROOT}/manifests" \
    "${F2_RUN_ROOT}/manifests/source.manifest"
  cp "${TEST_ROOT}/f2-run-source.manifest.saved" \
    "${F2_RUN_ROOT}/manifests/source.manifest"
  chmod a-w "${F2_RUN_ROOT}/manifests/source.manifest" \
    "${F2_RUN_ROOT}/manifests"
done

# Mutating a frozen tree file, retained archive, or source manifest fails
# before the bound worker can execute.
cp "${F2_SOURCE_ROOT}/config_helper.R" "${TEST_ROOT}/f2-config.saved"
chmod -R u+w "${F2_SNAPSHOT_ROOT}"
printf 'tree mutation\n' >> "${F2_SOURCE_ROOT}/config_helper.R"
chmod -R a-w "${F2_SNAPSHOT_ROOT}"
expect_fail ecoda_runtime_validate_bound_run
chmod -R u+w "${F2_SNAPSHOT_ROOT}"
cp "${TEST_ROOT}/f2-config.saved" "${F2_SOURCE_ROOT}/config_helper.R"
chmod -R a-w "${F2_SNAPSHOT_ROOT}"

cp "${F2_SOURCE_ARCHIVE}" "${TEST_ROOT}/f2-source.tar.saved"
chmod u+w "${F2_SOURCE_IDENTITY}" "${F2_SOURCE_ARCHIVE}"
printf 'archive mutation\n' >> "${F2_SOURCE_ARCHIVE}"
chmod a-w "${F2_SOURCE_ARCHIVE}" "${F2_SOURCE_IDENTITY}"
expect_fail ecoda_runtime_validate_bound_run
chmod u+w "${F2_SOURCE_IDENTITY}" "${F2_SOURCE_ARCHIVE}"
cp "${TEST_ROOT}/f2-source.tar.saved" "${F2_SOURCE_ARCHIVE}"
chmod a-w "${F2_SOURCE_ARCHIVE}" "${F2_SOURCE_IDENTITY}"

cp "${F2_SOURCE_MANIFEST}" "${TEST_ROOT}/f2-source.manifest.before-mutation"
chmod u+w "${F2_SOURCE_IDENTITY}"
f2_replace_key "${F2_SOURCE_MANIFEST}" SOURCE_COMMIT \
  "cccccccccccccccccccccccccccccccccccccccc"
chmod a-w "${F2_SOURCE_MANIFEST}" "${F2_SOURCE_IDENTITY}"
expect_fail ecoda_runtime_validate_bound_run
chmod u+w "${F2_SOURCE_IDENTITY}" "${F2_SOURCE_MANIFEST}"
cp "${TEST_ROOT}/f2-source.manifest.before-mutation" "${F2_SOURCE_MANIFEST}"
chmod a-w "${F2_SOURCE_MANIFEST}" "${F2_SOURCE_IDENTITY}"

# Replacing a versioned image/manifest is rejected through bound path, size,
# digest, and read-only checks; no device-ID comparison is part of the proof.
cp "${F2_IMAGE}" "${TEST_ROOT}/f2-image.saved"
chmod u+w "${F2_RUNTIME_DIR}" "${F2_IMAGE}"
printf 'replacement image with a different size\n' > "${F2_IMAGE}"
chmod a-w "${F2_IMAGE}" "${F2_RUNTIME_DIR}"
expect_fail ecoda_runtime_validate_bound_run
chmod u+w "${F2_RUNTIME_DIR}" "${F2_IMAGE}"
cp "${TEST_ROOT}/f2-image.saved" "${F2_IMAGE}"
chmod a-w "${F2_IMAGE}" "${F2_RUNTIME_DIR}"

cp "${F2_RUNTIME_MANIFEST}" "${TEST_ROOT}/f2-runtime.manifest.before-mutation"
chmod u+w "${F2_RUNTIME_DIR}" "${F2_RUNTIME_MANIFEST}"
printf 'manifest replacement\n' >> "${F2_RUNTIME_MANIFEST}"
chmod a-w "${F2_RUNTIME_MANIFEST}" "${F2_RUNTIME_DIR}"
expect_fail ecoda_runtime_validate_bound_run
chmod u+w "${F2_RUNTIME_DIR}" "${F2_RUNTIME_MANIFEST}"
cp "${TEST_ROOT}/f2-runtime.manifest.before-mutation" "${F2_RUNTIME_MANIFEST}"
chmod a-w "${F2_RUNTIME_MANIFEST}" "${F2_RUNTIME_DIR}"

chmod u+w "${F2_RUNTIME_DIR}" "${F2_IMAGE}"
expect_fail ecoda_runtime_validate_bound_run
chmod a-w "${F2_IMAGE}" "${F2_RUNTIME_DIR}"

# A relocated Apptainer worker receives the immutable snapshot and separate
# writable operational paths.  The canonical checkout has a different worker
# at the same relative path, so executing it would produce the wrong marker.
cat > "${PROJECT}/src/snapshot_worker.sh" <<MUTABLE_WORKER
#!/bin/bash
set -euo pipefail
printf 'mutable-checkout\n' > "${F2_EXEC_MARKER}"
MUTABLE_WORKER
chmod +x "${PROJECT}/src/snapshot_worker.sh"
cat > "${PROJECT}/src/f2_outer_worker.sh" <<F2_OUTER_WORKER
#!/bin/bash
set -euo pipefail
source "${ROOT}/src/utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage4 "${F2_SOURCE_ROOT}/src/snapshot_worker.sh"
F2_OUTER_WORKER
chmod +x "${PROJECT}/src/f2_outer_worker.sh"
: > "${FAKE_APPTAINER_LOG}"
rm -f "${F2_EXEC_MARKER}"
bash "${PROJECT}/src/f2_outer_worker.sh"
assert_eq snapshot-A "$(cat "${F2_EXEC_MARKER}")" \
  "relocated worker executes snapshot source"
f2_apptainer_log="$(cat "${FAKE_APPTAINER_LOG}")"
assert_contains "${f2_apptainer_log}" "${F2_SOURCE_ROOT}/src/snapshot_worker.sh" \
  "snapshot worker script path"
assert_contains "${f2_apptainer_log}" \
  "bind=${F2_SOURCE_ROOT}/aux:${F2_SOURCE_ROOT}/aux:ro" \
  "snapshot auxiliary runtime bind"
assert_contains "${f2_apptainer_log}" "bind=${F2_LOGS}:${F2_LOGS}:rw" \
  "run logs runtime bind"
[[ ! -e "${F2_SOURCE_ROOT}/__pycache__" ]] || fail "snapshot worker created __pycache__"

# Snapshot-backed host workers must validate the run-bound runtime identity
# before their worker body is allowed to execute.
cat > "${PROJECT}/src/f2_host_worker.sh" <<F2_HOST_WORKER
#!/bin/bash
set -euo pipefail
source "${ROOT}/src/utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage4 "${F2_SOURCE_ROOT}/src/snapshot_worker.sh"
printf 'host-A\n' > "${F2_EXEC_MARKER}"
F2_HOST_WORKER
chmod +x "${PROJECT}/src/f2_host_worker.sh"
export ECODA_RUNTIME_MODE=host
rm -f "${F2_EXEC_MARKER}"
bash "${PROJECT}/src/f2_host_worker.sh"
assert_eq host-A "$(cat "${F2_EXEC_MARKER}")" \
  "host worker executes after bound runtime validation"
cp "${F2_RUN_ROOT}/manifests/runtime.identity" \
  "${TEST_ROOT}/f2-runtime.identity.host.saved"
chmod u+w "${F2_RUN_ROOT}/manifests" \
  "${F2_RUN_ROOT}/manifests/runtime.identity"
f2_replace_key "${F2_RUN_ROOT}/manifests/runtime.identity" \
  RUNTIME_MANIFEST_SHA256 \
  "0000000000000000000000000000000000000000000000000000000000000000"
chmod 600 "${F2_RUN_ROOT}/manifests/runtime.identity"
chmod a-w "${F2_RUN_ROOT}/manifests"
rm -f "${F2_EXEC_MARKER}"
expect_fail bash "${PROJECT}/src/f2_host_worker.sh"
[[ ! -e "${F2_EXEC_MARKER}" ]] || \
  fail "host worker executed after bound runtime identity mutation"
chmod u+w "${F2_RUN_ROOT}/manifests"
cp "${TEST_ROOT}/f2-runtime.identity.host.saved" \
  "${F2_RUN_ROOT}/manifests/runtime.identity"
chmod 600 "${F2_RUN_ROOT}/manifests/runtime.identity"
chmod a-w "${F2_RUN_ROOT}/manifests"
# A change to a recorded host interpreter is rejected before the host worker
# body executes, without changing the immutable runtime.identity schema.
cp "${ECODA_HOST_ENV_PREFIX}/bin/python" \
  "${TEST_ROOT}/f2-host-python.saved"
printf 'mutated host Python\n' > "${ECODA_HOST_ENV_PREFIX}/bin/python"
chmod +x "${ECODA_HOST_ENV_PREFIX}/bin/python"
rm -f "${F2_EXEC_MARKER}"
expect_fail bash "${PROJECT}/src/f2_host_worker.sh"
[[ ! -e "${F2_EXEC_MARKER}" ]] || \
  fail "host worker executed after host Python mutation"
cp "${TEST_ROOT}/f2-host-python.saved" \
  "${ECODA_HOST_ENV_PREFIX}/bin/python"
chmod +x "${ECODA_HOST_ENV_PREFIX}/bin/python"


# Direct symlinks for either run-bound identity are rejected by both bound
# validation and runtime export before any worker marker can be created.
chmod u+w "${F2_RUN_ROOT}/manifests"
mv "${F2_RUN_ROOT}/manifests/runtime.identity" \
  "${TEST_ROOT}/f2-runtime.identity.symlink-target"
ln -s "${TEST_ROOT}/f2-runtime.identity.symlink-target" \
  "${F2_RUN_ROOT}/manifests/runtime.identity"
expect_fail ecoda_runtime_validate_bound_run
expect_fail ecoda_runtime_export_csv stage4 0
rm -f "${F2_RUN_ROOT}/manifests/runtime.identity"
mv "${TEST_ROOT}/f2-runtime.identity.symlink-target" \
  "${F2_RUN_ROOT}/manifests/runtime.identity"
chmod 600 "${F2_RUN_ROOT}/manifests/runtime.identity"

mv "${F2_RUN_ROOT}/manifests/source.manifest" \
  "${TEST_ROOT}/f2-source.manifest.symlink-target"
ln -s "${TEST_ROOT}/f2-source.manifest.symlink-target" \
  "${F2_RUN_ROOT}/manifests/source.manifest"
expect_fail ecoda_runtime_validate_bound_run
rm -f "${F2_RUN_ROOT}/manifests/source.manifest"
mv "${TEST_ROOT}/f2-source.manifest.symlink-target" \
  "${F2_RUN_ROOT}/manifests/source.manifest"
chmod 600 "${F2_RUN_ROOT}/manifests/source.manifest"
chmod a-w "${F2_RUN_ROOT}/manifests"

# Required-snapshot host mode cannot execute a live checkout or a mutable
# nonmatching environment prefix.
export ECODA_RUNTIME_MODE=host
expect_fail ecoda_runtime_reexec_worker stage4 "${PROJECT}/src/snapshot_worker.sh"
export ECODA_HOST_ENV_PREFIX="${TEST_ROOT}/mutable-host-prefix"
export ECODA_HOST_PYTHON_BIN="${ECODA_HOST_ENV_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${ECODA_HOST_ENV_PREFIX}/bin/Rscript --vanilla"
expect_fail ecoda_runtime_reexec_worker stage4 \
  "${F2_SOURCE_ROOT}/src/snapshot_worker.sh"
export ECODA_HOST_ENV_PREFIX="${PROJECT}/.pixi/envs/py-cuda13"
export ECODA_HOST_PYTHON_BIN="${ECODA_HOST_ENV_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${ECODA_HOST_ENV_PREFIX}/bin/Rscript --vanilla"
export ECODA_RUNTIME_MODE=apptainer
[[ ! -e "${FAKE_SBATCH_LOG}" ]] || fail "runtime validation submitted a scheduler job"

echo "Immutable runtime contract: OK"
