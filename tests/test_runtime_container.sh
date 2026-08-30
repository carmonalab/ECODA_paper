#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

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
cleanup() {
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

echo "Immutable runtime contract: OK"
