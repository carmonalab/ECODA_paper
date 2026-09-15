#!/bin/bash
# Focused contract test for validator-only Stage 3 retry-2 acceptance.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage3-retry-accept.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

SOURCE_COMMIT="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SNAPSHOT_ROOT="${TMP_DIR}/source-snapshots/${SOURCE_COMMIT}"
SOURCE_TREE="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY}/source.tar"
RUNTIME_DIR="${TMP_DIR}/runtime/_ecoda_runtime/fixture-runtime"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-runtime.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/runs" "${TMP_DIR}/hpc/Alzheimer/output" \
  "${TMP_DIR}/nas/Alzheimer/output" "${TMP_DIR}/owners/artifact" \
  "${SOURCE_TREE}/src" "${SOURCE_TREE}/aux" "${SOURCE_IDENTITY}" "${RUNTIME_DIR}"

md5_file() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | cut -d' ' -f1
  else
    md5 -q "$1"
  fi
}
sha256_text() {
  if command -v sha256sum >/dev/null 2>&1; then
    printf '%s' "$1" | sha256sum | cut -d' ' -f1
  else
    printf '%s' "$1" | shasum -a 256 | cut -d' ' -f1
  fi
}
sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}
write_sidecar() {
  local path="$1"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "$(md5_file "${path}")" \
    "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
}

# The accounting and Python boundaries are deterministic stubs.  The sbatch
# stub is a tripwire: a valid acceptance must never reach a scheduler submit.
SACCT_CALLS="${TMP_DIR}/sacct.calls"
SBATCH_CALLED="${TMP_DIR}/sbatch.called"
MD5_CALLED="${TMP_DIR}/md5.called"
PYTHON_CALLS="${TMP_DIR}/python.calls"
export SACCT_CALLS SBATCH_CALLED MD5_CALLED PYTHON_CALLS
cat > "${TMP_DIR}/bin/sacct" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${SACCT_CALLS:?}"
case "$*" in
  *"-j 9003,9004"*)
    printf '9003|COMPLETED|0:0\n9004|COMPLETED|0:0\n'
    ;;
  *)
    echo "unexpected sacct query: $*" >&2
    exit 1
    ;;
esac
STUB
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
touch "${SBATCH_CALLED:?}"
echo 'sbatch must not be called by retry acceptance' >&2
exit 99
STUB
cat > "${TMP_DIR}/bin/benchmark_h5ad_contract.py" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${PYTHON_CALLS:?}"
# This file stands in for the configured Python interpreter's explicit
# benchmark_h5ad_contract.py target.  Check argv positions rather than
# relying on a substring match over "$*", so a path/name mismatch cannot be
# mistaken for a valid semantic invocation.
[[ "$#" -eq 7 ]] || {
  echo "semantic validator received an unexpected argument count" >&2
  exit 78
}
[[ "$(basename "$1")" == benchmark_h5ad_contract.py ]] || {
  echo "semantic validator target path was not propagated" >&2
  exit 78
}
shift
[[ "$1" == --path && "${2}" == "${EXPECTED_SEMANTIC_PATH:?}" &&
   "$3" == --view && "$4" == batch_effect_uncorrected &&
   "$5" == --method && "$6" == "Stage 3 preprocessing" &&
   -z "${7:-}" ]] || {
  echo "semantic validator arguments were not explicit" >&2
  exit 78
}
STUB
cp "${TMP_DIR}/bin/benchmark_h5ad_contract.py" "${TMP_DIR}/bin/python-validator"
chmod +x "${TMP_DIR}/bin/sacct" "${TMP_DIR}/bin/sbatch" \
  "${TMP_DIR}/bin/benchmark_h5ad_contract.py" "${TMP_DIR}/bin/python-validator"
export PATH="${TMP_DIR}/bin:${PATH}"

SOURCE_RUN_ID="stage3-alzheimer-failed"
SOURCE_RUN_ROOT="${TMP_DIR}/runs/${SOURCE_RUN_ID}"
NEW_RUN_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery"
CONFIG_PATH="${SOURCE_TREE}/datasets.json"
printf 'fixture helper\n' > "${SOURCE_TREE}/config_helper.R"
printf 'fixture pixi\n' > "${SOURCE_TREE}/pixi.toml"
printf 'fixture lock\n' > "${SOURCE_TREE}/pixi.lock"
printf 'fixture source\n' > "${SOURCE_TREE}/src/fixture.py"
for aux_file in scGateDB.rds genes.blocklist.rds EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
  printf 'fixture aux\n' > "${SOURCE_TREE}/aux/${aux_file}"
done
printf '%s\n' '{"Alzheimer":{"views":{"batch_effect_uncorrected":{"output_file_name":"SEAAD_Alzheimer_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad"}}}}' > "${CONFIG_PATH}"
printf 'COMPLETE\n' > "${SNAPSHOT_ROOT}/COMPLETE"
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_TREE}" .
cat > "${SOURCE_MANIFEST}" <<EOF
FORMAT=1
SOURCE_ROOT=${SOURCE_TREE}
SOURCE_COMMIT=${SOURCE_COMMIT}
SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}
SOURCE_ARCHIVE_SHA256=$(sha256_file "${SOURCE_ARCHIVE}")
CONFIG_HELPER_SHA256=$(sha256_file "${SOURCE_TREE}/config_helper.R")
DATASETS_SHA256=$(sha256_file "${SOURCE_TREE}/datasets.json")
PIXI_TOML_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.toml")
PIXI_LOCK_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.lock")
AUX_ROOT=${SOURCE_TREE}/aux
SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4
EOF
chmod -R a-w "${SNAPSHOT_ROOT}"

RUNTIME_IMAGE_SHA256=""
printf 'synthetic format-2 runtime image\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA256="$(sha256_file "${RUNTIME_IMAGE}")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_BUILD_GIT_REVISION=fixture-runtime-build
IMAGE_PATH=${RUNTIME_IMAGE}
IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.49.0
APPTAINER_VERSION=1.3.2
IMAGE_PIXI_TOML_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.toml")
IMAGE_PIXI_LOCK_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.lock")
EOF
RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")"
RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
RUNTIME_IDENTITY_FIXTURE="${TMP_DIR}/runtime.identity.fixture"
cat > "${RUNTIME_IDENTITY_FIXTURE}" <<EOF
RUNTIME_IMAGE=${RUNTIME_IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
RUNTIME_IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}
RUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA256}
RUNTIME_IMAGE_SIZE=${RUNTIME_IMAGE_SIZE}
RUNTIME_MANIFEST_SIZE=${RUNTIME_MANIFEST_SIZE}
IMAGE_PIXI_TOML_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.toml")
IMAGE_PIXI_LOCK_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.lock")
EOF
chmod -R a-w "${RUNTIME_DIR}"
mkdir -p "${SOURCE_RUN_ROOT}/manifests/artifacts" "${SOURCE_RUN_ROOT}/status"
cp -p "${SOURCE_MANIFEST}" "${SOURCE_RUN_ROOT}/manifests/source.manifest"
cp -p "${RUNTIME_IDENTITY_FIXTURE}" "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
printf 'STAGE=stage3\nRUN_ID=%s\nSTATE=ACTIVE\n' "${SOURCE_RUN_ID}" > "${SOURCE_RUN_ROOT}/metadata"
export HPC_SCRATCH_DIR="${TMP_DIR}/hpc" NAS_TARGET_DIR="${TMP_DIR}/nas"

SELECTION="${SOURCE_RUN_ROOT}/manifests/selection.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\n' > "${SELECTION}"

SCRATCH_H5AD="${TMP_DIR}/hpc/Alzheimer/output/SEAAD_Alzheimer_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad"
NAS_H5AD="${TMP_DIR}/nas/Alzheimer/output/SEAAD_Alzheimer_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad"
printf 'retry2-processed-h5ad\n' > "${SCRATCH_H5AD}"
printf 'retry2-processed-h5ad\n' > "${NAS_H5AD}"
write_sidecar "${SCRATCH_H5AD}"
write_sidecar "${NAS_H5AD}"
chmod a-w "${SCRATCH_H5AD}" "${NAS_H5AD}"
export EXPECTED_SEMANTIC_PATH="${SCRATCH_H5AD}"

OWNER_DIR="${TMP_DIR}/owners/artifact/SEAAD_Alzheimer_batch_effect_analysis_uncorrected_ECODAprocessed"
ARTIFACT_OWNER_ONE="${TMP_DIR}/owners/artifact/array-initial"
ARTIFACT_OWNER_TWO="${TMP_DIR}/owners/artifact/retry1"
mkdir -p "${OWNER_DIR}" "${ARTIFACT_OWNER_ONE}" "${ARTIFACT_OWNER_TWO}"
printf 'RUN_ID=%s\nSTAGE=stage3\nPATH=%s\nSTATE=OK\nPID=1\nREASON=retry2 published\n' \
  "${SOURCE_RUN_ID}" "${SCRATCH_H5AD}" > "${OWNER_DIR}/owner"
printf 'Alzheimer/batch_effect_uncorrected\t%s\nARTIFACT\t%s\nARTIFACT\t%s\n' \
  "${OWNER_DIR}" "${ARTIFACT_OWNER_ONE}" "${ARTIFACT_OWNER_TWO}" \
  > "${SOURCE_RUN_ROOT}/manifests/owners.tsv"
printf 'Alzheimer\tbatch_effect_uncorrected\t%s\t%s\t%s\n' \
  "${SCRATCH_H5AD}" "${NAS_H5AD}" "${OWNER_DIR}" \
  > "${SOURCE_RUN_ROOT}/manifests/output_ownership.tsv"
printf 'ARRAY\t9001\nWATCHDOG\t9004\nSTATUS\t9002\nSTATUS\t9003\n' \
  > "${SOURCE_RUN_ROOT}/manifests/scheduler_ids.tsv"
printf 'STATE=OK\nRUN_ID=%s\nRETRY_INDEX=2\nARRAY_JOB_ID=9003\nSCHEDULER_ID=9004\nSCHEDULER_ID=9003\nSCHEDULER_ID=9001\nSCHEDULER_ID=9002\n' \
  "${SOURCE_RUN_ID}" > "${SOURCE_RUN_ROOT}/status/watchdog"
printf 'STATE=OK\nRUN_ID=%s\nREASON=watchdog accepted\n' "${SOURCE_RUN_ID}" \
  > "${SOURCE_RUN_ROOT}/status/terminal"

# The old durable inspect is intentionally failed because its prior attempts
# are OOM; this plain accounting form is also accepted by the generic parser.
PRIOR_INSPECT="${TMP_DIR}/prior.inspect"
{
  printf 'GATE_ID=%s\nCOMMAND_DIGEST=fixture-command-digest\nEVENT_GENERATION=fixture-event-generation\nSOURCE_RUN_ID=%s\nSELECTION_PATH=%s\nSOURCE_MANIFEST=%s\nSTATE=FAILED\nAUDIT_PASSED=false\nACCOUNTING_OK=false\nQUERY_COUNT=1\n' \
    "${SOURCE_RUN_ID}" "${SOURCE_RUN_ID}" "${SELECTION}" \
    "${SOURCE_RUN_ROOT}/manifests/source.manifest"
  printf '9001|OUT_OF_MEMORY|0:125\n9002|OUT_OF_MEMORY|0:125\n9003|COMPLETED|0:0\n9004|COMPLETED|0:0\n'
} > "${PRIOR_INSPECT}"

ARTIFACT_DIGEST="$(sha256_text "${SCRATCH_H5AD}")"
printf 'PATH=%s\nSIZE=%s\nMD5=%s\nRUN_ID=%s\nPRODUCER=stage3\nSTATE=PUBLISHED\n' \
  "${SCRATCH_H5AD}" "$(wc -c < "${SCRATCH_H5AD}" | tr -d '[:space:]')" \
  "$(md5_file "${SCRATCH_H5AD}")" "${SOURCE_RUN_ID}" \
  > "${SOURCE_RUN_ROOT}/manifests/artifacts/${ARTIFACT_DIGEST:0:32}.record"
cat > "${TMP_DIR}/bin/md5sum" <<'STUB'
#!/bin/bash
set -euo pipefail
touch "${MD5_CALLED:?}"
echo 'unexpected MD5 rehash' >&2
exit 99
STUB
cat > "${TMP_DIR}/bin/md5" <<'STUB'
#!/bin/bash
set -euo pipefail
touch "${MD5_CALLED:?}"
echo 'unexpected MD5 rehash' >&2
exit 99
STUB
chmod +x "${TMP_DIR}/bin/md5sum" "${TMP_DIR}/bin/md5"

# Preserve source evidence byte-for-byte to prove the validator does not edit
# the failed run while it creates only the fresh recovery root.
for source_evidence in \
  "${SOURCE_RUN_ROOT}/metadata" \
  "${SOURCE_RUN_ROOT}/manifests/source.manifest" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity" \
  "${SOURCE_RUN_ROOT}/manifests/selection.tsv" \
  "${SOURCE_RUN_ROOT}/manifests/output_ownership.tsv" \
  "${SOURCE_RUN_ROOT}/manifests/owners.tsv" \
  "${SOURCE_RUN_ROOT}/manifests/scheduler_ids.tsv" \
  "${SOURCE_RUN_ROOT}/status/watchdog" \
  "${SOURCE_RUN_ROOT}/status/terminal"; do
  cp "${source_evidence}" "${source_evidence}.before"
done
SOURCE_MANIFEST_SNAPSHOT_BEFORE="${TMP_DIR}/source.manifest.snapshot.before"
RUNTIME_MANIFEST_BEFORE="${TMP_DIR}/runtime.manifest.before"
RUNTIME_IMAGE_BEFORE="${TMP_DIR}/runtime.image.before"
cp "${SOURCE_MANIFEST}" "${SOURCE_MANIFEST_SNAPSHOT_BEFORE}"
cp "${RUNTIME_MANIFEST}" "${RUNTIME_MANIFEST_BEFORE}"
cp "${RUNTIME_IMAGE}" "${RUNTIME_IMAGE_BEFORE}"

"${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
  --source-run-id "${SOURCE_RUN_ID}" \
  --source-run-root "${SOURCE_RUN_ROOT}" \
  --selection "${SELECTION}" \
  --prior-inspect "${PRIOR_INSPECT}" \
  --accepted-array-id 9003 \
  --accepted-watchdog-id 9004 \
  --superseded-attempt-ids 9001,9002 \
  --run-id stage3-alzheimer-recovery \
  --run-root "${NEW_RUN_ROOT}" \
  --config "${CONFIG_PATH}" \
  --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
  --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null

[[ -s "${NEW_RUN_ROOT}/status/terminal" ]]
[[ "$(sed -n 's/^STATE=//p' "${NEW_RUN_ROOT}/status/terminal")" == OK ]]
[[ -s "${NEW_RUN_ROOT}/reports/stage3_retry_acceptance.json" ]]
[[ -s "${NEW_RUN_ROOT}/metadata" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/source.manifest" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/runtime.identity" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/prior_inspect" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/watchdog.status" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/output_ownership.tsv" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/owners.tsv" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/scheduler_ids.tsv" ]]
[[ -s "${NEW_RUN_ROOT}/manifests/source_terminal" ]]
cmp -s "${NEW_RUN_ROOT}/manifests/source.manifest" "${SOURCE_MANIFEST}"
cmp -s "${NEW_RUN_ROOT}/manifests/runtime.identity" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
[[ ! -d "${NEW_RUN_ROOT}/manifests/artifacts" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == 1 ]]
[[ "$(cat "${SACCT_CALLS}")" == *"-j 9003,9004"* ]]
[[ "$(cat "${SACCT_CALLS}")" != *"9001"* && "$(cat "${SACCT_CALLS}")" != *"9002"* ]]
[[ ! -e "${SBATCH_CALLED}" ]]
[[ ! -e "${MD5_CALLED}" ]]
[[ "$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')" == 1 ]]

ATTEMPTS="${NEW_RUN_ROOT}/manifests/accepted_attempts.tsv"
[[ "$(cat "${ATTEMPTS}")" == $'ATTEMPT_ID\tROLE\tSTATE\tEXIT_CODE
9001\tsuperseded\tOUT_OF_MEMORY\t0:125
9002\tsuperseded\tOUT_OF_MEMORY\t0:125
9003\taccepted_array\tCOMPLETED\t0:0
9004\taccepted_watchdog\tCOMPLETED\t0:0' ]]
REPORT_CONTENT="$(cat "${NEW_RUN_ROOT}/reports/stage3_retry_acceptance.json")"
case "${REPORT_CONTENT}" in
  *'"state":"OK"'*) ;;
  *) echo 'acceptance report did not record terminal OK' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"attempt_id":"9001"'*'"state":"OUT_OF_MEMORY"'*'"attempt_id":"9003"'*'"state":"COMPLETED"'*) ;;
  *) echo 'acceptance report did not preserve the full attempt chain' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"source_metadata_state":"ACTIVE"'*'"prior_inspect_identity"'*'"gate_id":"stage3-alzheimer-failed"'*'"command_digest":"fixture-command-digest"'*) ;;
  *) echo 'acceptance report did not record bound source/prior identities' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"source_identity"'*'"snapshot_manifest"'*'"commit":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"'*) ;;
  *) echo 'acceptance report did not record immutable source identity' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"runtime_identity"'*) ;;
  *) echo 'acceptance report did not record runtime identity' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"format":"2"'*) ;;
  *) echo 'acceptance report did not record runtime format 2' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"image":"'"${RUNTIME_IMAGE}"'"'*) ;;
  *) echo 'acceptance report did not record the bound runtime image' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"manifest":"'"${RUNTIME_MANIFEST}"'"'*) ;;
  *) echo 'acceptance report did not record the bound runtime manifest' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"image_sha256":"'"${RUNTIME_IMAGE_SHA256}"'"'*) ;;
  *) echo 'acceptance report did not record the runtime image digest' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"manifest_sha256":"'"${RUNTIME_MANIFEST_SHA256}"'"'*) ;;
  *) echo 'acceptance report did not record the runtime manifest digest' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"validator"'*'"path":"'"${TMP_DIR}/bin/benchmark_h5ad_contract.py"'"'*) ;;
  *) echo 'acceptance report did not record the semantic validator path' >&2; exit 1 ;;
esac
case "${REPORT_CONTENT}" in
  *'"validator"'*'"semantic":true'*) ;;
  *) echo 'acceptance report did not record semantic validation' >&2; exit 1 ;;
esac

for source_evidence in \
  "${SOURCE_RUN_ROOT}/metadata" \
  "${SOURCE_RUN_ROOT}/manifests/source.manifest" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity" \
  "${SOURCE_RUN_ROOT}/manifests/selection.tsv" \
  "${SOURCE_RUN_ROOT}/manifests/output_ownership.tsv" \
  "${SOURCE_RUN_ROOT}/manifests/owners.tsv" \
  "${SOURCE_RUN_ROOT}/manifests/scheduler_ids.tsv" \
  "${SOURCE_RUN_ROOT}/status/watchdog" \
  "${SOURCE_RUN_ROOT}/status/terminal"; do
  cmp -s "${source_evidence}" "${source_evidence}.before"
done
cmp -s "${SOURCE_MANIFEST}" "${SOURCE_MANIFEST_SNAPSHOT_BEFORE}"
cmp -s "${RUNTIME_MANIFEST}" "${RUNTIME_MANIFEST_BEFORE}"
cmp -s "${RUNTIME_IMAGE}" "${RUNTIME_IMAGE_BEFORE}"

# A source-run manifest that is altered after the snapshot was copied must
# fail before accounting or creation of a recovery root.
BAD_SOURCE_BIND_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-source-bind"
SACCT_LINES_BEFORE="$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')"
chmod u+w "${SOURCE_RUN_ROOT}/manifests/source.manifest"
printf 'TAMPERED=source-run-manifest\n' >> "${SOURCE_RUN_ROOT}/manifests/source.manifest"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/source.manifest"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${PRIOR_INSPECT}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-source-bind \
    --run-root "${BAD_SOURCE_BIND_ROOT}" \
    --config "${CONFIG_PATH}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null 2>&1; then
  echo 'altered source-run manifest was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_SOURCE_BIND_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "${SACCT_LINES_BEFORE}" ]]
chmod u+w "${SOURCE_RUN_ROOT}/manifests/source.manifest"
cp "${SOURCE_RUN_ROOT}/manifests/source.manifest.before" \
  "${SOURCE_RUN_ROOT}/manifests/source.manifest"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/source.manifest"

BAD_SNAPSHOT_KEY_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-snapshot-key"
awk '{ if ($0 ~ /^SOURCE_COMMIT=/) print "SOURCE_COMMIT=fixture"; else print }' \
  "${SOURCE_RUN_ROOT}/manifests/source.manifest.before" \
  > "${TMP_DIR}/source.manifest.non-snapshot"
chmod u+w "${SOURCE_RUN_ROOT}/manifests/source.manifest"
cp "${TMP_DIR}/source.manifest.non-snapshot" \
  "${SOURCE_RUN_ROOT}/manifests/source.manifest"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/source.manifest"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${PRIOR_INSPECT}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-snapshot-key \
    --run-root "${BAD_SNAPSHOT_KEY_ROOT}" \
    --config "${CONFIG_PATH}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null 2>&1; then
  echo 'non-snapshot source commit key was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_SNAPSHOT_KEY_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "${SACCT_LINES_BEFORE}" ]]
chmod u+w "${SOURCE_RUN_ROOT}/manifests/source.manifest"
cp "${SOURCE_RUN_ROOT}/manifests/source.manifest.before" \
  "${SOURCE_RUN_ROOT}/manifests/source.manifest"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/source.manifest"

# Runtime manifest and identity mutations must fail closed on their recorded
# path/size/digest bindings, without touching the old run or rehashing H5AD.
BAD_RUNTIME_BIND_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-runtime-bind"
awk '!/^IMAGE_SHA256=/' "${RUNTIME_MANIFEST_BEFORE}" \
  > "${TMP_DIR}/runtime.manifest.tampered"
printf 'IMAGE_SHA256=%064d\n' 0 >> "${TMP_DIR}/runtime.manifest.tampered"
chmod u+w "${RUNTIME_MANIFEST}"
cp "${TMP_DIR}/runtime.manifest.tampered" "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_MANIFEST}"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${PRIOR_INSPECT}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-runtime-bind \
    --run-root "${BAD_RUNTIME_BIND_ROOT}" \
    --config "${CONFIG_PATH}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null 2>&1; then
  echo 'altered runtime manifest was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_RUNTIME_BIND_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "${SACCT_LINES_BEFORE}" ]]
chmod u+w "${RUNTIME_MANIFEST}"
cp "${RUNTIME_MANIFEST_BEFORE}" "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_MANIFEST}"
# A duplicate runtime-manifest key must fail closed after the immutable
# runtime bindings are updated to the duplicate file itself.  This reaches
# the shape validator rather than merely failing on a stale digest.
BAD_RUNTIME_DUPLICATE_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-runtime-duplicate"
RUNTIME_MANIFEST_DUPLICATE="${TMP_DIR}/runtime.manifest.duplicate"
cp "${RUNTIME_MANIFEST_BEFORE}" "${RUNTIME_MANIFEST_DUPLICATE}"
chmod u+w "${RUNTIME_MANIFEST_DUPLICATE}"
printf 'EXTRA_FIELD=one\nEXTRA_FIELD=two\n' >> "${RUNTIME_MANIFEST_DUPLICATE}"
DUPLICATE_RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST_DUPLICATE}")"
DUPLICATE_RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST_DUPLICATE}" | tr -d '[:space:]')"
awk -v sha="${DUPLICATE_RUNTIME_MANIFEST_SHA256}" \
  -v size="${DUPLICATE_RUNTIME_MANIFEST_SIZE}" '
  /^RUNTIME_MANIFEST_SHA256=/ {
    print "RUNTIME_MANIFEST_SHA256=" sha
    next
  }
  /^RUNTIME_MANIFEST_SIZE=/ {
    print "RUNTIME_MANIFEST_SIZE=" size
    next
  }
  { print }
' "${SOURCE_RUN_ROOT}/manifests/runtime.identity.before" \
  > "${TMP_DIR}/runtime.identity.duplicate"
chmod u+w "${RUNTIME_MANIFEST}"
cp "${RUNTIME_MANIFEST_DUPLICATE}" "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_MANIFEST}"
chmod u+w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
cp "${TMP_DIR}/runtime.identity.duplicate" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${PRIOR_INSPECT}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-runtime-duplicate \
    --run-root "${BAD_RUNTIME_DUPLICATE_ROOT}" \
    --config "${CONFIG_PATH}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null 2>&1; then
  echo 'duplicate runtime manifest key was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_RUNTIME_DUPLICATE_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "${SACCT_LINES_BEFORE}" ]]
[[ "$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')" == 1 ]]
[[ ! -e "${MD5_CALLED}" ]]
chmod u+w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
cp "${SOURCE_RUN_ROOT}/manifests/runtime.identity.before" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
chmod u+w "${RUNTIME_MANIFEST}"
cp "${RUNTIME_MANIFEST_BEFORE}" "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_MANIFEST}"


BAD_RUNTIME_IDENTITY_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-runtime-identity"
awk '!/^RUNTIME_IMAGE_SHA256=/' \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity" \
  > "${TMP_DIR}/runtime.identity.tampered"
printf 'RUNTIME_IMAGE_SHA256=%064d\n' 0 >> "${TMP_DIR}/runtime.identity.tampered"
chmod u+w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
cp "${TMP_DIR}/runtime.identity.tampered" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${PRIOR_INSPECT}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-runtime-identity \
    --run-root "${BAD_RUNTIME_IDENTITY_ROOT}" \
    --config "${CONFIG_PATH}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null 2>&1; then
  echo 'altered runtime identity was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_RUNTIME_IDENTITY_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "${SACCT_LINES_BEFORE}" ]]
chmod u+w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
cp "${SOURCE_RUN_ROOT}/manifests/runtime.identity.before" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity"
chmod a-w "${SOURCE_RUN_ROOT}/manifests/runtime.identity"

# An un-hashed source-tree file is still protected by the immutable archive
# comparison rather than being silently accepted after re-sealing.
BAD_SOURCE_TREE_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-source-tree"
SOURCE_TREE_FILE="${SOURCE_TREE}/src/fixture.py"
cp "${SOURCE_TREE_FILE}" "${TMP_DIR}/source-tree-file.before"
chmod u+w "${SOURCE_TREE_FILE}"
printf 'tampered source tree\n' > "${SOURCE_TREE_FILE}"
chmod a-w "${SOURCE_TREE_FILE}"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${PRIOR_INSPECT}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-source-tree \
    --run-root "${BAD_SOURCE_TREE_ROOT}" \
    --config "${CONFIG_PATH}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null 2>&1; then
  echo 'tampered un-hashed source-tree file was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_SOURCE_TREE_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "${SACCT_LINES_BEFORE}" ]]
chmod u+w "${SOURCE_TREE_FILE}"
cp "${TMP_DIR}/source-tree-file.before" "${SOURCE_TREE_FILE}"
chmod a-w "${SOURCE_TREE_FILE}"

# A malformed superseded attempt must fail before the single accepted-ID sacct
# query and before any fresh run root or scheduler boundary is touched.
BAD_PRIOR="${TMP_DIR}/prior.bad.inspect"
sed 's/^9002|OUT_OF_MEMORY|0:125$/9002|FAILED|1:0/' "${PRIOR_INSPECT}" > "${BAD_PRIOR}"
BAD_RUN_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-bad"
SACCT_LINES_BEFORE="$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${BAD_PRIOR}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-bad \
    --run-root "${BAD_RUN_ROOT}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null 2>&1; then
  echo 'malformed superseded attempt was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_RUN_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "${SACCT_LINES_BEFORE}" ]]
[[ ! -e "${SBATCH_CALLED}" ]]

# Optional prior-inspect bindings are genuinely optional.  Omitting them must
# not trip set -u or manufacture an identity; the required gate identity and
# accounting chronology remain enforced.
OPTIONAL_PRIOR="${TMP_DIR}/prior.optional.inspect"
sed -e '/^SOURCE_RUN_ID=/d' -e '/^SELECTION_PATH=/d' -e '/^SOURCE_MANIFEST=/d' \
  "${PRIOR_INSPECT}" > "${OPTIONAL_PRIOR}"
OPTIONAL_RUN_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-optional"
"${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
  --source-run-id "${SOURCE_RUN_ID}" \
  --source-run-root "${SOURCE_RUN_ROOT}" \
  --selection "${SELECTION}" \
  --prior-inspect "${OPTIONAL_PRIOR}" \
  --accepted-array-id 9003 \
  --accepted-watchdog-id 9004 \
  --superseded-attempt-ids 9001,9002 \
  --run-id stage3-alzheimer-recovery-optional \
  --run-root "${OPTIONAL_RUN_ROOT}" \
  --config "${CONFIG_PATH}" \
  --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
  --python-bin "${TMP_DIR}/bin/python-validator" >/dev/null
OPTIONAL_REPORT="$(cat "${OPTIONAL_RUN_ROOT}/reports/stage3_retry_acceptance.json")"
case "${OPTIONAL_REPORT}" in
  *'"selection_binding":""'*'"source_manifest_binding":""'*) ;;
  *) echo 'optional prior-inspect bindings were not initialized safely' >&2; exit 1 ;;
esac
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == 2 ]]
[[ "$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')" == 2 ]]
[[ ! -e "${SBATCH_CALLED}" ]]
[[ ! -e "${MD5_CALLED}" ]]
cmp -s "${SOURCE_RUN_ROOT}/manifests/source.manifest" \
  "${SOURCE_RUN_ROOT}/manifests/source.manifest.before"
cmp -s "${SOURCE_RUN_ROOT}/manifests/runtime.identity" \
  "${SOURCE_RUN_ROOT}/manifests/runtime.identity.before"
cmp -s "${SOURCE_MANIFEST}" "${SOURCE_MANIFEST_SNAPSHOT_BEFORE}"
cmp -s "${RUNTIME_MANIFEST}" "${RUNTIME_MANIFEST_BEFORE}"
cmp -s "${RUNTIME_IMAGE}" "${RUNTIME_IMAGE_BEFORE}"

# A configured interpreter that drops the required method argument must fail
# closed: semantic validation is a release prerequisite, not a best-effort
# check, and no acceptance root may be published.
BAD_SEMANTIC_BIN="${TMP_DIR}/bin/python-validator-malformed"
cat > "${BAD_SEMANTIC_BIN}" <<'STUB'
#!/bin/bash
set -euo pipefail
script="$1"
output="$3"
view="$5"
exec "${script}" --path "${output}" --view "${view}"
STUB
chmod +x "${BAD_SEMANTIC_BIN}"
BAD_SEMANTIC_ROOT="${TMP_DIR}/runs/stage3-alzheimer-recovery-semantic-bad"
SACCT_LINES_BEFORE="$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')"
PYTHON_LINES_BEFORE="$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')"
if "${ROOT}/src/3_scrnaseq_preprocessing/stage3_retry_acceptance_validator.sh" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --selection "${SELECTION}" \
    --prior-inspect "${PRIOR_INSPECT}" \
    --accepted-array-id 9003 \
    --accepted-watchdog-id 9004 \
    --superseded-attempt-ids 9001,9002 \
    --run-id stage3-alzheimer-recovery-semantic-bad \
    --run-root "${BAD_SEMANTIC_ROOT}" \
    --config "${CONFIG_PATH}" \
    --validator-script "${TMP_DIR}/bin/benchmark_h5ad_contract.py" \
    --python-bin "${BAD_SEMANTIC_BIN}" >/dev/null 2>&1; then
  echo 'malformed semantic invocation was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_SEMANTIC_ROOT}" ]]
[[ "$(wc -l < "${SACCT_CALLS}" | tr -d '[:space:]')" == "$((SACCT_LINES_BEFORE + 1))" ]]
[[ "$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')" == "$((PYTHON_LINES_BEFORE + 1))" ]]
[[ ! -e "${SBATCH_CALLED}" ]]
[[ ! -e "${MD5_CALLED}" ]]

printf 'stage3 retry acceptance validator: OK\n'
