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

assert_file() {
  local path="$1"
  local label="${2:-file}"
  [[ -f "${path}" && ! -L "${path}" ]] || fail "${label} is missing or unsafe: ${path}"
}

assert_dir() {
  local path="$1"
  local label="${2:-directory}"
  [[ -d "${path}" && ! -L "${path}" ]] || fail "${label} is missing or unsafe: ${path}"
}

expect_fail() {
  if "$@"; then
    fail "expected failure: $*"
  fi
}

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

manifest_value() {
  awk -F= -v key="$2" '$1 == key { print substr($0, length(key) + 2); found++ } END { exit(found == 1 ? 0 : 1) }' "$1"
}

TEST_TMP_BASE="${TMPDIR:-/tmp}"
TEST_TMP_BASE="${TEST_TMP_BASE%/}"
TEST_ROOT="$(mktemp -d "${TEST_TMP_BASE}/ecoda-source-snapshot.XXXXXX")"
TEST_ROOT="$(cd "${TEST_ROOT}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TEST_ROOT}" >/dev/null 2>&1 || true
  rm -rf "${TEST_ROOT}"
}
trap cleanup EXIT

SOURCE_ROOT="${TEST_ROOT}/source"
SNAPSHOT_PARENT="${TEST_ROOT}/snapshots"
HOST_ENV_PREFIX="${TEST_ROOT}/host-env/.pixi/envs/py-cuda13"
SCRATCH_ROOT="${TEST_ROOT}/scratch"
LOGS_PARENT="${TEST_ROOT}/logs"
LOGS_ROOT="${LOGS_PARENT}/run-A"
RUN_ID="source-snapshot-run-A"
MARKER="${TEST_ROOT}/worker.marker"
HOOK_MARKER="${TEST_ROOT}/hook.marker"
HOOK_FILE="${TEST_ROOT}/bash-hook.sh"
RUNTIME_DIR="${SCRATCH_ROOT}/_ecoda_runtime/runtime-A"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
RUN_ROOT="${SCRATCH_ROOT}/_ecoda_runs/${RUN_ID}"

mkdir -p \
  "${SOURCE_ROOT}/src/utils/bash" \
  "${SOURCE_ROOT}/aux" \
  "${SNAPSHOT_PARENT}" \
  "${HOST_ENV_PREFIX}/bin" \
  "${SCRATCH_ROOT}" \
  "${LOGS_PARENT}"

# Keep the fixture deliberately small while retaining the source files required
# by the production snapshot contract.  Every auxiliary input is tracked in
# exactly one root: source/aux.
cp "${ROOT}/src/utils/bash/ecoda_source_snapshot.sh" \
  "${SOURCE_ROOT}/src/utils/bash/ecoda_source_snapshot.sh"
cp "${ROOT}/src/slurm_config.sh" "${SOURCE_ROOT}/src/slurm_config.sh"
chmod +x "${SOURCE_ROOT}/src/utils/bash/ecoda_source_snapshot.sh"
printf 'config-A\n' > "${SOURCE_ROOT}/config_helper.R"
printf '{"fixture":"A"}\n' > "${SOURCE_ROOT}/datasets.json"
printf 'pixi-toml-A\n' > "${SOURCE_ROOT}/pixi.toml"
printf 'pixi-lock-A\n' > "${SOURCE_ROOT}/pixi.lock"
printf 'scgate-db-A\n' > "${SOURCE_ROOT}/aux/scGateDB.rds"
printf 'blocklist-A\n' > "${SOURCE_ROOT}/aux/genes.blocklist.rds"
printf 'ensembl-map-A\n' > "${SOURCE_ROOT}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
cat > "${SOURCE_ROOT}/src/snapshot_worker.sh" <<'WORKER_A'
#!/bin/bash
set -euo pipefail
: "${SNAPSHOT_MARKER:?}"
: "${EXPECTED_SNAPSHOT_ROOT:?}"
[[ "${ECODA_SOURCE_ROOT}" == "${EXPECTED_SNAPSHOT_ROOT}" ]] || exit 31
[[ "${ECODA_SOURCE_MANIFEST}" == "${EXPECTED_SNAPSHOT_ROOT%/tree}/identity/source.manifest" ]] || exit 32
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED}" == 1 ]] || exit 33
[[ "${ECODA_AUX_ROOT}" == "${EXPECTED_SNAPSHOT_ROOT}/aux" ]] || exit 34
[[ -f "${ECODA_AUX_ROOT}/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz" ]] || exit 35
[[ "${PYTHONDONTWRITEBYTECODE}" == 1 ]] || exit 36
[[ "${ECODA_RUNTIME_MODE}" == apptainer ]] || exit 41
[[ -z "${BASH_ENV:-}" && -z "${ENV:-}" ]] || exit 37
if mkdir "${ECODA_SOURCE_ROOT}/__pycache__" 2>/dev/null; then
  exit 38
fi
[[ ! -e "${ECODA_SOURCE_ROOT}/__pycache__" ]] || exit 39
if [[ -n "${REPLACEMENT_TARGET:-}" ]]; then
  if mv "${EXPECTED_SNAPSHOT_ROOT%/tree}" "${REPLACEMENT_TARGET}" 2>/dev/null; then
    exit 40
  fi
fi
{
  printf 'A\n'
  printf '%s\n' "${ECODA_SOURCE_ROOT}"
  printf '%s\n' "${ECODA_AUX_ROOT}"
  printf '%s\n' "${ECODA_LOGS_DIR}"
} > "${SNAPSHOT_MARKER}"
WORKER_A
chmod +x "${SOURCE_ROOT}/src/snapshot_worker.sh"
cp "${SOURCE_ROOT}/src/snapshot_worker.sh" "${TEST_ROOT}/worker-A.sh"

# Commit A is the immutable identity used for the first snapshot.
git -C "${SOURCE_ROOT}" init -q
git -C "${SOURCE_ROOT}" config user.email ecoda-test@example.invalid
git -C "${SOURCE_ROOT}" config user.name ecoda-test

git -C "${SOURCE_ROOT}" add .
git -C "${SOURCE_ROOT}" commit -qm 'fixture commit A'
COMMIT_A="$(git -C "${SOURCE_ROOT}" rev-parse HEAD)"
[[ "${COMMIT_A}" =~ ^[0-9a-f]{40}$ ]] || fail "commit A is not a full Git identity"

SNAPSHOT_EXECUTOR="${SOURCE_ROOT}/src/utils/bash/ecoda_source_snapshot.sh"
"${SNAPSHOT_EXECUTOR}" create \
  --source-root "${SOURCE_ROOT}" \
  --snapshot-parent "${SNAPSHOT_PARENT}" \
  --commit "${COMMIT_A}" >/dev/null

SNAPSHOT_A="${SNAPSHOT_PARENT}/${COMMIT_A}"
SNAPSHOT_A_TREE="${SNAPSHOT_A}/tree"
SNAPSHOT_A_IDENTITY="${SNAPSHOT_A}/identity"
SNAPSHOT_A_MANIFEST="${SNAPSHOT_A_IDENTITY}/source.manifest"
SNAPSHOT_A_ARCHIVE="${SNAPSHOT_A_IDENTITY}/source.tar"
assert_dir "${SNAPSHOT_A_TREE}" "snapshot A tree"
assert_file "${SNAPSHOT_A_MANIFEST}" "snapshot A manifest"
assert_file "${SNAPSHOT_A_ARCHIVE}" "snapshot A archive"
assert_eq "COMPLETE" "$(tr -d '\n' < "${SNAPSHOT_A}/COMPLETE")" "snapshot completion marker"
assert_eq "${SNAPSHOT_A_TREE}" "$(manifest_value "${SNAPSHOT_A_MANIFEST}" SOURCE_ROOT)" "manifest source root"
assert_eq "${SNAPSHOT_A_TREE}/aux" "$(manifest_value "${SNAPSHOT_A_MANIFEST}" AUX_ROOT)" "manifest auxiliary root"
assert_eq "${COMMIT_A}" "$(manifest_value "${SNAPSHOT_A_MANIFEST}" SOURCE_COMMIT)" "manifest source commit"

for aux_file in \
  aux/scGateDB.rds \
  aux/genes.blocklist.rds \
  aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
  assert_file "${SNAPSHOT_A_TREE}/${aux_file}" "tracked auxiliary file ${aux_file}"
done
[[ ! -e "${SNAPSHOT_A}/aux" ]] || fail "snapshot published a sibling auxiliary root"
assert_eq "$(cat "${TEST_ROOT}/worker-A.sh")" \
  "$(cat "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh")" \
  "snapshot A worker content"

# A later checkout commit must not alter the already published commit-keyed
# tree or archive.
cat > "${SOURCE_ROOT}/src/snapshot_worker.sh" <<'WORKER_B'
#!/bin/bash
set -euo pipefail
printf 'B\n' > "${SNAPSHOT_MARKER:?}"
WORKER_B
chmod +x "${SOURCE_ROOT}/src/snapshot_worker.sh"
git -C "${SOURCE_ROOT}" add src/snapshot_worker.sh
git -C "${SOURCE_ROOT}" commit -qm 'fixture commit B'
COMMIT_B="$(git -C "${SOURCE_ROOT}" rev-parse HEAD)"
[[ "${COMMIT_B}" != "${COMMIT_A}" ]] || fail "commit B did not advance the fixture"
assert_eq "$(cat "${TEST_ROOT}/worker-A.sh")" \
  "$(cat "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh")" \
  "snapshot A remained unchanged after commit B"

# Both dirty tracked state and untracked files are rejected before publication.
printf 'dirty checkout\n' >> "${SOURCE_ROOT}/config_helper.R"
expect_fail "${SNAPSHOT_EXECUTOR}" create \
  --source-root "${SOURCE_ROOT}" \
  --snapshot-parent "${SNAPSHOT_PARENT}" \
  --commit "${COMMIT_B}" >/dev/null 2>&1
cp "${SNAPSHOT_A_TREE}/config_helper.R" "${SOURCE_ROOT}/config_helper.R"
touch "${SOURCE_ROOT}/untracked.fixture"
expect_fail "${SNAPSHOT_EXECUTOR}" create \
  --source-root "${SOURCE_ROOT}" \
  --snapshot-parent "${SNAPSHOT_PARENT}" \
  --commit "${COMMIT_B}" >/dev/null 2>&1
rm -f "${SOURCE_ROOT}/untracked.fixture"

# A path already occupied by a non-directory cannot be adopted as a commit
# snapshot, even when the source checkout itself is valid and clean.
CONFLICT_PARENT="${TEST_ROOT}/conflict-parent"
mkdir -p "${CONFLICT_PARENT}"
printf 'identity conflict\n' > "${CONFLICT_PARENT}/${COMMIT_B}"
expect_fail "${SNAPSHOT_EXECUTOR}" create \
  --source-root "${SOURCE_ROOT}" \
  --snapshot-parent "${CONFLICT_PARENT}" \
  --commit "${COMMIT_B}" >/dev/null 2>&1
rm -f "${CONFLICT_PARENT}/${COMMIT_B}"

# A clean later commit can still receive its own immutable snapshot without
# replacing snapshot A.
"${SNAPSHOT_EXECUTOR}" create \
  --source-root "${SOURCE_ROOT}" \
  --snapshot-parent "${SNAPSHOT_PARENT}" \
  --commit "${COMMIT_B}" >/dev/null
SNAPSHOT_B="${SNAPSHOT_PARENT}/${COMMIT_B}"
assert_file "${SNAPSHOT_B}/tree/src/snapshot_worker.sh" "snapshot B worker"
assert_eq "${COMMIT_A}" "$(manifest_value "${SNAPSHOT_A_MANIFEST}" SOURCE_COMMIT)" "snapshot A identity after snapshot B"

# Build the minimal format-2 runtime identity consumed by exec.  This is a
# fake SIF: the source executor validates its identity and never invokes HPC.
mkdir -p "${RUNTIME_DIR}"
printf 'fake-runtime-image-A\n' > "${RUNTIME_IMAGE}"
IMAGE_SHA256="$(sha256_file "${RUNTIME_IMAGE}")"
IMAGE_TOML_SHA256="$(manifest_value "${SNAPSHOT_A_MANIFEST}" PIXI_TOML_SHA256)"
IMAGE_LOCK_SHA256="$(manifest_value "${SNAPSHOT_A_MANIFEST}" PIXI_LOCK_SHA256)"
{
  printf '%s\n' \
    'FORMAT=2' \
    "IMAGE_PATH=${RUNTIME_IMAGE}" \
    "IMAGE_SHA256=${IMAGE_SHA256}" \
    'RUNTIME_ENV=py-cuda13' \
    'RUNTIME_LAYOUT=relocated' \
    'CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13' \
    'BASE_IMAGE=rockylinux:9' \
    'PIXITAINER_VERSION=0.8.3' \
    'PIXI_VERSION=0.49.0' \
    'APPTAINER_VERSION=1.3.2' \
    'IMAGE_BUILD_GIT_REVISION=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa' \
    "IMAGE_PIXI_TOML_SHA256=${IMAGE_TOML_SHA256}" \
    "IMAGE_PIXI_LOCK_SHA256=${IMAGE_LOCK_SHA256}"
} > "${RUNTIME_MANIFEST}"
chmod -R a-w "${RUNTIME_DIR}"

# Invoke the bootstrap copied into snapshot A itself.  BASH_ENV and ENV are
# installed only after that shell is running, so a valid executor must unset
# them before starting the worker.  The worker also proves that ECODA_AUX_ROOT
# points at the frozen map rather than at the mutable checkout.
cat > "${HOOK_FILE}" <<'HOOK'
printf 'startup hook injected code\n' > "${HOOK_MARKER:?}"
export ECODA_AUX_ROOT=/tmp/injected-aux
HOOK
chmod 600 "${HOOK_FILE}"

run_snapshot_exec() {
  env -u BASH_ENV -u ENV \
    SNAPSHOT_EXECUTOR="${SNAPSHOT_A}/tree/src/utils/bash/ecoda_source_snapshot.sh" \
    SNAPSHOT_MARKER="${MARKER}" \
    EXPECTED_SNAPSHOT_ROOT="${SNAPSHOT_A_TREE}" \
    HOOK_MARKER="${HOOK_MARKER}" \
    HOOK_FILE="${HOOK_FILE}" \
    REPLACEMENT_TARGET="${REPLACEMENT_TARGET:-}" \
    bash -c '
      set +e
      set --
      source "${SNAPSHOT_EXECUTOR}" >/dev/null 2>&1
      set -e
      BASH_ENV="${HOOK_FILE}"
      ENV="${HOOK_FILE}"
      export BASH_ENV ENV
      exec_snapshot \
        "${EXPECTED_SNAPSHOT_ROOT}" \
        "${EXPECTED_SNAPSHOT_ROOT%/tree}/identity/source.manifest" \
        "'"${HOST_ENV_PREFIX}"'" \
        "'"${RUNTIME_IMAGE}"'" \
        "'"${RUNTIME_MANIFEST}"'" \
        "'"${RUN_ID}"'" \
        "'"${SCRATCH_ROOT}"'" \
        "'"${LOGS_ROOT}"'" \
        src/snapshot_worker.sh
    '
}

run_snapshot_exec
assert_eq A "$(sed -n '1p' "${MARKER}")" "snapshot executor worker identity"
assert_eq "${SNAPSHOT_A_TREE}" "$(sed -n '2p' "${MARKER}")" "snapshot executor source root"
assert_eq "${SNAPSHOT_A_TREE}/aux" "$(sed -n '3p' "${MARKER}")" "ECODA_AUX_ROOT resolution"
assert_eq "${LOGS_ROOT}" "$(sed -n '4p' "${MARKER}")" "separate executor logs root"
[[ ! -e "${HOOK_MARKER}" ]] || fail "BASH_ENV/ENV startup hook injected worker code"
assert_dir "${LOGS_ROOT}" "executor-created logs root"
[[ ! -e "${RUN_ROOT}" ]] || fail "executor created a stage run root"
[[ ! -e "${SNAPSHOT_A_TREE}/__pycache__" ]] || fail "snapshot worker created __pycache__"

# The executor must make the shared parent non-replaceable only for the
# worker's lifetime.  The worker attempts the same-parent replacement; a
# successful rename would exit before writing the marker.
REPLACEMENT_PATH="${SNAPSHOT_PARENT}/replacement-${COMMIT_A}"
REPLACEMENT_TARGET="${REPLACEMENT_PATH}"
run_snapshot_exec
unset REPLACEMENT_TARGET
assert_eq A "$(sed -n '1p' "${MARKER}")" "snapshot survived replacement attempt"
[[ ! -e "${REPLACEMENT_PATH}" ]] || fail "worker replaced the commit snapshot"
[[ -w "${SNAPSHOT_PARENT}" ]] || fail "snapshot parent remained locked after execution"
[[ ! -e "${SNAPSHOT_PARENT}/.ecoda-exec-lock" ]] || fail "snapshot execution lock was not released"
"${SNAPSHOT_EXECUTOR}" create \
  --source-root "${SOURCE_ROOT}" \
  --snapshot-parent "${SNAPSHOT_PARENT}" \
  --commit "${COMMIT_B}" >/dev/null

# Symlinked snapshot ancestry is rejected before the worker starts.
SNAPSHOT_PARENT_REAL="${TEST_ROOT}/snapshots-real"
mv "${SNAPSHOT_PARENT}" "${SNAPSHOT_PARENT_REAL}"
ln -s "${SNAPSHOT_PARENT_REAL}" "${SNAPSHOT_PARENT}"
expect_fail run_snapshot_exec
assert_eq A "$(sed -n '1p' "${MARKER}")" "symlinked ancestry blocked worker"
rm "${SNAPSHOT_PARENT}"
mv "${SNAPSHOT_PARENT_REAL}" "${SNAPSHOT_PARENT}"

# Mutating either the published tree or retained archive is detected before
# the snapshot worker can execute.  Restore the fixture and its read-only
# contract after each negative case.
cp "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh" "${TEST_ROOT}/worker-A-saved.sh"
chmod -R u+w "${SNAPSHOT_A}"
printf 'mutated snapshot tree\n' > "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh"
chmod +x "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh"
chmod -R a-w "${SNAPSHOT_A}"
expect_fail run_snapshot_exec
chmod -R u+w "${SNAPSHOT_A}"
cp "${TEST_ROOT}/worker-A-saved.sh" "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh"
chmod +x "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh"
chmod -R a-w "${SNAPSHOT_A}"

cp "${SNAPSHOT_A_ARCHIVE}" "${TEST_ROOT}/source-A.tar.saved"
chmod u+w "${SNAPSHOT_A_IDENTITY}" "${SNAPSHOT_A_ARCHIVE}"
printf 'archive mutation\n' >> "${SNAPSHOT_A_ARCHIVE}"
chmod a-w "${SNAPSHOT_A_ARCHIVE}" "${SNAPSHOT_A_IDENTITY}"
expect_fail run_snapshot_exec
chmod u+w "${SNAPSHOT_A_IDENTITY}" "${SNAPSHOT_A_ARCHIVE}"
cp "${TEST_ROOT}/source-A.tar.saved" "${SNAPSHOT_A_ARCHIVE}"
chmod a-w "${SNAPSHOT_A_ARCHIVE}" "${SNAPSHOT_A_IDENTITY}"

assert_eq "$(cat "${TEST_ROOT}/worker-A.sh")" \
  "$(cat "${SNAPSHOT_A_TREE}/src/snapshot_worker.sh")" \
  "snapshot A restored after mutation checks"

printf 'Source snapshot isolation: OK\n'
