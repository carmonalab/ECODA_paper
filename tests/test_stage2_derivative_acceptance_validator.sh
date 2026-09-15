#!/usr/bin/env bash
# Deterministic contract test for the validator-only Stage 2 Alzheimer
# derivative acceptance wrapper.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage2-derivative-accept.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

VALIDATOR="${ROOT}/src/2_dataset_specific_preprocessing/stage2_derivative_acceptance_validator.sh"
REAL_MD5SUM="$(command -v md5sum || true)"
REAL_SHA256SUM="$(command -v sha256sum || true)"
REAL_PIXI="$(command -v pixi || true)"
[[ -n "${REAL_MD5SUM}" && -n "${REAL_SHA256SUM}" && -n "${REAL_PIXI}" ]] || {
  echo "test requires md5sum, sha256sum, and pixi" >&2
  exit 2
}

SOURCE_COMMIT="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SNAPSHOT_ROOT="${TMP_DIR}/source-snapshots/${SOURCE_COMMIT}"
SOURCE_TREE="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY}/source.tar"
RUNTIME_DIR="${TMP_DIR}/runtime/_ecoda_runtime/alzheimer"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
RUNTIME_IDENTITY="${TMP_DIR}/runtime.identity"
HPC_ROOT="${TMP_DIR}/hpc"
SCRATCH_ROOT="${HPC_ROOT}"
RUNS_ROOT="${TMP_DIR}/runs"
SOURCE_RUN_ID="stage2_alzheimer_donor_assay_20260914T172421Z"
SOURCE_RUN_ROOT="${RUNS_ROOT}/${SOURCE_RUN_ID}"
NEW_RUN_ROOT="${RUNS_ROOT}/stage2_alzheimer_donor_assay_acceptance"
ARTIFACT="${SCRATCH_ROOT}/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad"
RAW_INPUT="${SCRATCH_ROOT}/Alzheimer/data/SEAAD_Alzheimer.h5ad"
ARTIFACT_RECORD="${SOURCE_RUN_ROOT}/manifests/artifacts/donor-assay.record"
OWNER_DIR="${TMP_DIR}/owners/stage2/alzheimer_donor_assay"
OWNER_FILE="${OWNER_DIR}/owner"
SOURCE_TERMINAL="${SOURCE_RUN_ROOT}/status/terminal"
WATCHDOG_STATUS="${SOURCE_RUN_ROOT}/status/watchdog"
STEPS_MANIFEST="${SOURCE_RUN_ROOT}/manifests/steps.tsv"
OWNERSHIP_MANIFEST="${SOURCE_RUN_ROOT}/manifests/ownership.tsv"
JOBS_MANIFEST="${SOURCE_RUN_ROOT}/manifests/jobs.tsv"
SCHEDULER_MANIFEST="${SOURCE_RUN_ROOT}/manifests/scheduler_ids.tsv"
CONFIG="${SOURCE_TREE}/datasets.json"
STAGE2_HOOK="${SOURCE_TREE}/src/2_dataset_specific_preprocessing/1.7_submit_alzheimer_donor_assay.sh"
PYTHON_STUB="${TMP_DIR}/bin/python-validator"
SBATCH_STUB="${TMP_DIR}/bin/sbatch"
WRITER_STUB="${TMP_DIR}/bin/ecoda_write_artifact_record"
PYTHON_CALLS="${TMP_DIR}/python.calls"
SBATCH_CALLED="${TMP_DIR}/sbatch.called"
WRITER_CALLED="${TMP_DIR}/writer.called"
REHASH_CALLED="${TMP_DIR}/artifact.rehashed"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/source-snapshots" \
  "${SOURCE_TREE}/src/2_dataset_specific_preprocessing" \
  "${SOURCE_TREE}/src/utils/py" "${SOURCE_TREE}/aux" "${SOURCE_IDENTITY}" \
  "${RUNTIME_DIR}" "${SCRATCH_ROOT}/Alzheimer/data" \
  "${SOURCE_RUN_ROOT}/manifests/artifacts" "${SOURCE_RUN_ROOT}/manifests" \
  "${SOURCE_RUN_ROOT}/status" "${RUNS_ROOT}" "${OWNER_DIR}"

cp "${ROOT}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py" \
  "${SOURCE_TREE}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py"
cat > "${SOURCE_TREE}/src/utils/py/h5ad_source_identity.py" <<'PY'
import numpy as np


def _decode_string(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return value.tobytes().decode("utf-8")
    return value


def _decode_array(values):
    return np.asarray([_decode_string(value) for value in values], dtype=object)


def read_obs_column_values(obs_group, column):
    node = obs_group[column]
    encoding = node.attrs.get("encoding-type")
    if isinstance(encoding, bytes):
        encoding = encoding.decode("utf-8")
    if encoding in (None, "", "array", "string-array"):
        return _decode_array(node[:])
    if encoding == "categorical":
        categories = _decode_array(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        return np.where(codes >= 0, categories[np.clip(codes, 0, None)], "nan")
    raise RuntimeError(f"unsupported fixture obs encoding: {encoding!r}")
PY
cp "${ROOT}/src/2_dataset_specific_preprocessing/1.7_submit_alzheimer_donor_assay.sh" \
  "${STAGE2_HOOK}"
printf 'fixture config helper\n' > "${SOURCE_TREE}/config_helper.R"
printf 'fixture pixi\n' > "${SOURCE_TREE}/pixi.toml"
printf 'fixture lock\n' > "${SOURCE_TREE}/pixi.lock"
for aux_file in scGateDB.rds genes.blocklist.rds EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; do
  printf 'fixture aux\n' > "${SOURCE_TREE}/aux/${aux_file}"
done
cat > "${CONFIG}" <<'JSON'
{
  "Alzheimer": {
    "file_names": "SEAAD_Alzheimer.h5ad",
    "views": {
      "batch_effect_uncorrected": {
        "input_file_name": "SEAAD_Alzheimer_donor_assay.h5ad"
      },
      "batch_effect_corrected": {
        "input_file_name": "SEAAD_Alzheimer_donor_assay.h5ad"
      }
    }
  }
}
JSON
printf 'COMPLETE\n' > "${SNAPSHOT_ROOT}/COMPLETE"
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_TREE}" .
sha256_file() { "${REAL_SHA256SUM}" "$1" | cut -d' ' -f1; }
cat > "${SOURCE_MANIFEST}" <<EOF
FORMAT=1
SOURCE_ROOT=${SOURCE_TREE}
SOURCE_COMMIT=${SOURCE_COMMIT}
SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}
SOURCE_ARCHIVE_SHA256=$(sha256_file "${SOURCE_ARCHIVE}")
CONFIG_HELPER_SHA256=$(sha256_file "${SOURCE_TREE}/config_helper.R")
DATASETS_SHA256=$(sha256_file "${CONFIG}")
PIXI_TOML_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.toml")
PIXI_LOCK_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.lock")
AUX_ROOT=${SOURCE_TREE}/aux
SCGATE_DB_BRANCH=41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4
EOF

printf 'synthetic format-2 runtime image\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA="$(sha256_file "${RUNTIME_IMAGE}")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_BUILD_GIT_REVISION=fixture-runtime-build
IMAGE_PATH=${RUNTIME_IMAGE}
IMAGE_SHA256=${RUNTIME_IMAGE_SHA}
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
RUNTIME_MANIFEST_SHA="$(sha256_file "${RUNTIME_MANIFEST}")"
cat > "${RUNTIME_IDENTITY}" <<EOF
RUNTIME_IMAGE=${RUNTIME_IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
RUNTIME_IMAGE_SHA256=${RUNTIME_IMAGE_SHA}
RUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA}
RUNTIME_IMAGE_SIZE=$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')
RUNTIME_MANIFEST_SIZE=$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')
IMAGE_PIXI_TOML_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.toml")
IMAGE_PIXI_LOCK_SHA256=$(sha256_file "${SOURCE_TREE}/pixi.lock")
EOF

FIXTURE_GENERATOR="${TMP_DIR}/make_h5ad.py"
cat > "${FIXTURE_GENERATOR}" <<'PY'
import sys
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

output = sys.argv[1]
n = 1_395_601
pattern_donors = ["H20.33.001", "H20.33.001"]
pattern_assays = ["10x 3' v3", "10x multiome"]
pattern_sexes = ["female", "female"]
for index in range(1, 58):
    pattern_donors.append(f"D{index:02d}")
    pattern_assays.append("10x 3' v3")
    pattern_sexes.append("female")
for index in range(58, 83):
    pattern_donors.append(f"D{index:02d}")
    pattern_assays.append("10x 3' v3")
    pattern_sexes.append("male")
for index in range(58, 78):
    pattern_donors.append(f"D{index:02d}")
    pattern_assays.append("10x multiome")
    pattern_sexes.append("male")

def repeat_pattern(values):
    return np.resize(np.asarray(values, dtype=object), n)

obs = pd.DataFrame(
    {
        "donor_id": pd.Categorical(repeat_pattern(pattern_donors)),
        "assay": pd.Categorical(repeat_pattern(pattern_assays)),
        "sex": pd.Categorical(repeat_pattern(pattern_sexes)),
        "Cognitive status": pd.Categorical(np.full(n, "Reference", dtype=object)),
    },
    index=pd.RangeIndex(n, name="_index"),
)
# Keep expression payloads genuinely sparse: only the persisted row pointers
# scale with the approved observation count; no dense n-by-gene array exists.
matrix = sp.csr_matrix((n, 2), dtype=np.float32)
data = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=["g1", "g2"]))
data.layers["counts"] = matrix
data.write_h5ad(output)
PY
"${REAL_PIXI}" run -e default python "${FIXTURE_GENERATOR}" "${RAW_INPUT}"
"${REAL_PIXI}" run -e default python \
  "${ROOT}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py" \
  --input-file "${RAW_INPUT}" --output-file "${ARTIFACT}" \
  --expected-samples 104 --expected-donors 83 \
  --expected-assay-counts "10x3v3=83,10xmultiome=21" \
  --expected-sex-counts "female=59,male=45" --require-example-ids
md5_file() { "${REAL_MD5SUM}" "$1" | cut -d' ' -f1; }
write_sidecar() {
  local path="$1"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "$(md5_file "${path}")" \
    "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
}
write_sidecar "${ARTIFACT}"
ARTIFACT_MD5="$(md5_file "${ARTIFACT}")"
ARTIFACT_SIZE="$(wc -c < "${ARTIFACT}" | tr -d '[:space:]')"
printf 'PATH=%s\nSIZE=%s\nMD5=%s\nRUN_ID=%s\nPRODUCER=alzheimer_donor_assay\nSTATE=PUBLISHED\n' \
  "${ARTIFACT}" "${ARTIFACT_SIZE}" "${ARTIFACT_MD5}" "${SOURCE_RUN_ID}" \
  > "${ARTIFACT_RECORD}"
printf 'RUN_ID=%s\nSTATE=OK\nSTAGE=stage2\nKEY=alzheimer_donor_assay\nPID=4242\nREASON=published\n' \
  "${SOURCE_RUN_ID}" > "${OWNER_FILE}"
printf 'STAGE=stage2\nRUN_ID=%s\nSTATE=ACTIVE\n' "${SOURCE_RUN_ID}" \
  > "${SOURCE_RUN_ROOT}/metadata"
printf 'STATE=OK\nRUN_ID=%s\nARRAY_JOB_ID=4407671\nWATCHDOG_JOB_ID=4407672\nREASON=completed\n' \
  "${SOURCE_RUN_ID}" > "${SOURCE_TERMINAL}"
printf 'STATE=OK\nRUN_ID=%s\nRETRY_INDEX=0\nSCHEDULER_ID=4407671\nSCHEDULER_ID=4407672\n' \
  "${SOURCE_RUN_ID}" > "${WATCHDOG_STATUS}"
printf 'alzheimer_donor_assay\t%s\t%s\t-\t%s\n' \
  "${STAGE2_HOOK}" "${ARTIFACT}" "${OWNER_DIR}" > "${STEPS_MANIFEST}"
printf 'alzheimer_donor_assay\t%s\t%s\t-\t%s\n' \
  "${STAGE2_HOOK}" "${ARTIFACT}" "${OWNER_DIR}" > "${OWNERSHIP_MANIFEST}"
printf 'alzheimer_donor_assay\t4407671\n' > "${JOBS_MANIFEST}"
printf 'ARRAY\t4407671\nWATCHDOG\t4407672\n' > "${SCHEDULER_MANIFEST}"
printf 'STATE=FAILED\nRUN_ID=%s\n' "${SOURCE_RUN_ID}" > "${SOURCE_RUN_ROOT}/status/bad-terminal"
cp "${SOURCE_MANIFEST}" "${SOURCE_RUN_ROOT}/manifests/source.manifest"
cp "${RUNTIME_IDENTITY}" "${SOURCE_RUN_ROOT}/manifests/runtime.identity"

# Keep the old run, source snapshot, runtime, and immutable H5AD read-only.
chmod -R a-w "${SNAPSHOT_ROOT}" "${RUNTIME_DIR}" "${SOURCE_RUN_ROOT}" \
  "${RUNTIME_IDENTITY}" "${RAW_INPUT}" "${ARTIFACT}" "${ARTIFACT}.md5"

cat > "${PYTHON_STUB}" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "${PYTHON_CALLS:?}"
[[ "${SEMANTIC_FAIL:-0}" != 1 ]] || exit 77
if [[ "$1" == *"1.7.1_create_alzheimer_donor_assay.py"* ]]; then
  [[ "$*" == *"--input-file ${RAW_INPUT}"* &&
     "$*" == *"--output-file ${ARTIFACT}"* &&
     "$*" == *"--validate-only"* &&
     "$*" == *"--skip-checksum"* &&
     "$*" == *"--expected-samples 104"* &&
     "$*" == *"--expected-donors 83"* &&
     "$*" == *"--expected-assay-counts 10x3v3=83,10xmultiome=21"* &&
     "$*" == *"--expected-sex-counts female=59,male=45"* &&
     "$*" == *"--require-example-ids"* ]] || {
    echo "semantic validator arguments were not explicit" >&2
    exit 78
  }
  if [[ "${SEMANTIC_WRONG_CELLS:-0}" == 1 ]]; then
    printf 'ALZHEIMER_DONOR_ASSAY_VALIDATED=1 cells=104 samples=104\n'
    exit 0
  fi
  cd "${ROOT:?}"
  exec "${REAL_PIXI:?}" run -e default python "$@"
fi
echo "unexpected Python target" >&2
exit 79
STUB
cat > "${SBATCH_STUB}" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
touch "${SBATCH_CALLED:?}"
echo 'scheduler submission is forbidden' >&2
exit 99
STUB
cat > "${WRITER_STUB}" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
touch "${WRITER_CALLED:?}"
echo 'artifact-record writer is forbidden' >&2
exit 99
STUB
cat > "${TMP_DIR}/bin/md5sum" <<STUB
#!/usr/bin/env bash
set -euo pipefail
for argument in "\$@"; do
  if [[ "\${argument}" == "${ARTIFACT}" ]]; then
    touch "${REHASH_CALLED}"
    echo 'immutable derivative was rehashed' >&2
    exit 98
  fi
done
exec "${REAL_MD5SUM}" "\$@"
STUB
chmod +x "${PYTHON_STUB}" "${SBATCH_STUB}" "${WRITER_STUB}" "${TMP_DIR}/bin/md5sum"
: > "${PYTHON_CALLS}"
export HPC_SCRATCH_DIR="${SCRATCH_ROOT}" PATH="${TMP_DIR}/bin:${PATH}"
export REAL_PIXI ROOT RAW_INPUT ARTIFACT PYTHON_CALLS SBATCH_CALLED \
  WRITER_CALLED REHASH_CALLED
export PYTHONDONTWRITEBYTECODE=1
unset ECODA_RUN_ROOT ECODA_RUN_ID

BEFORE_DIR="${TMP_DIR}/before"
mkdir -p "${BEFORE_DIR}"
save_before() {
  cp "$1" "${BEFORE_DIR}/$2"
}
save_before "${SOURCE_RUN_ROOT}/metadata" source.metadata
save_before "${SOURCE_RUN_ROOT}/manifests/source.manifest" source.manifest
save_before "${SOURCE_RUN_ROOT}/manifests/runtime.identity" runtime.identity
save_before "${SOURCE_RUN_ROOT}/manifests/steps.tsv" steps.tsv
save_before "${SOURCE_RUN_ROOT}/manifests/ownership.tsv" ownership.tsv
save_before "${SOURCE_RUN_ROOT}/manifests/jobs.tsv" jobs.tsv
save_before "${SOURCE_RUN_ROOT}/manifests/scheduler_ids.tsv" scheduler_ids.tsv
save_before "${SOURCE_RUN_ROOT}/status/watchdog" source_watchdog
save_before "${SOURCE_RUN_ROOT}/status/terminal" source_terminal
save_before "${ARTIFACT_RECORD}" artifact.record
save_before "${ARTIFACT}.md5" artifact.md5
save_before "${OWNER_FILE}" source_owner

run_acceptance() {
  local output_root="$1"
  shift
  "${VALIDATOR}" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --source-root "${SOURCE_TREE}" \
    --source-manifest "${SOURCE_MANIFEST}" \
    --runtime-identity "${RUNTIME_IDENTITY}" \
    --source-terminal "${SOURCE_TERMINAL}" \
    --source-owner "${OWNER_DIR}" \
    --artifact-record "${ARTIFACT_RECORD}" \
    --sidecar "${ARTIFACT}.md5" \
    --artifact-path "${ARTIFACT}" \
    --scratch-root "${SCRATCH_ROOT}" \
    --config "${CONFIG}" \
    --python-bin "${PYTHON_STUB}" \
    --validator-script "${SOURCE_TREE}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py" \
    --run-id "$(basename "${output_root}")" \
    --run-root "${output_root}" \
    --scheduler-array-id 4407671 \
    --scheduler-watchdog-id 4407672 \
    "$@"
}

run_acceptance "${NEW_RUN_ROOT}" > "${TMP_DIR}/acceptance.out"
[[ "$(sed -n '1p' "${TMP_DIR}/acceptance.out")" == NOOP_VALIDATED=1 ]]
[[ "$(sed -n '2p' "${TMP_DIR}/acceptance.out")" == STAGE2_DERIVATIVE_ACCEPTANCE_RUN_ID=* ]]
[[ "$(sed -n '3p' "${TMP_DIR}/acceptance.out")" == STAGE2_DERIVATIVE_ACCEPTANCE_REPORT=* ]]
REPORT="${NEW_RUN_ROOT}/reports/stage2_derivative_acceptance.json"
[[ -s "${REPORT}" && -s "${REPORT}.md5" ]]
[[ -s "${NEW_RUN_ROOT}/metadata" && -s "${NEW_RUN_ROOT}/metadata.md5" ]]
[[ -s "${NEW_RUN_ROOT}/status/terminal" && -s "${NEW_RUN_ROOT}/status/terminal.md5" ]]
[[ "$(sed -n 's/^STATE=//p' "${NEW_RUN_ROOT}/status/terminal")" == OK ]]
grep -q '^STATUS=NOOP_VALIDATED$' "${NEW_RUN_ROOT}/status/terminal"
grep -q '"state":"OK"' "${REPORT}"
grep -q '"status":"NOOP_VALIDATED"' "${REPORT}"
grep -q '"cells":1395601' "${REPORT}"
grep -q '"observed_cells":1395601' "${REPORT}"
grep -q '^ALZHEIMER_DONOR_ASSAY_VALIDATED=1 cells=1395601 samples=104$' \
  "${NEW_RUN_ROOT}/reports/semantic_validation.log"
grep -q '"samples":104' "${REPORT}"
grep -q '"donors":83' "${REPORT}"
grep -q '"array_id":"4407671"' "${REPORT}"
grep -q '"watchdog_id":"4407672"' "${REPORT}"
[[ ! -e "${NEW_RUN_ROOT}/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad" ]]
[[ "$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')" == 1 ]]
[[ ! -e "${SBATCH_CALLED}" && ! -e "${WRITER_CALLED}" && ! -e "${REHASH_CALLED}" ]]
for copied in source.metadata source.manifest runtime.identity steps.tsv \
  ownership.tsv jobs.tsv scheduler_ids.tsv source_watchdog source_terminal \
  source_owner artifact.record artifact.md5; do
  [[ -s "${NEW_RUN_ROOT}/manifests/${copied}" ]]
done
for pair in \
  "source.metadata:${SOURCE_RUN_ROOT}/metadata" \
  "source.manifest:${SOURCE_RUN_ROOT}/manifests/source.manifest" \
  "runtime.identity:${SOURCE_RUN_ROOT}/manifests/runtime.identity" \
  "steps.tsv:${STEPS_MANIFEST}" \
  "ownership.tsv:${OWNERSHIP_MANIFEST}" \
  "jobs.tsv:${JOBS_MANIFEST}" \
  "scheduler_ids.tsv:${SCHEDULER_MANIFEST}" \
  "source_watchdog:${WATCHDOG_STATUS}" \
  "source_terminal:${SOURCE_TERMINAL}" \
  "artifact.record:${ARTIFACT_RECORD}" \
  "artifact.md5:${ARTIFACT}.md5" \
  "source_owner:${OWNER_FILE}"; do
  copied="${pair%%:*}"
  original="${pair#*:}"
  cmp -s "${original}" "${BEFORE_DIR}/${copied}"
done

# Wrong prebound ID is rejected before the semantic interpreter or any fresh
# root is touched.
BAD_ID_ROOT="${RUNS_ROOT}/bad-id"
if run_acceptance "${BAD_ID_ROOT}" --scheduler-array-id 4407670 > "${TMP_DIR}/bad-id.log" 2>&1; then
  echo 'wrong scheduler ID was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_ID_ROOT}" ]]
[[ "$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')" == 1 ]]

# A terminal marker with a failed state is malformed source evidence and cannot
# create a partial acceptance root.
BAD_TERMINAL="${SOURCE_RUN_ROOT}/status/bad-terminal"
BAD_TERMINAL_ROOT="${RUNS_ROOT}/bad-terminal"
if "${VALIDATOR}" \
    --source-run-id "${SOURCE_RUN_ID}" \
    --source-run-root "${SOURCE_RUN_ROOT}" \
    --source-root "${SOURCE_TREE}" \
    --source-manifest "${SOURCE_MANIFEST}" \
    --runtime-identity "${RUNTIME_IDENTITY}" \
    --source-terminal "${BAD_TERMINAL}" \
    --source-owner "${OWNER_DIR}" \
    --artifact-record "${ARTIFACT_RECORD}" \
    --sidecar "${ARTIFACT}.md5" \
    --artifact-path "${ARTIFACT}" \
    --scratch-root "${SCRATCH_ROOT}" \
    --config "${CONFIG}" \
    --python-bin "${PYTHON_STUB}" \
    --validator-script "${SOURCE_TREE}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py" \
    --run-id bad-terminal --run-root "${BAD_TERMINAL_ROOT}" \
    --scheduler-array-id 4407671 --scheduler-watchdog-id 4407672 \
    > "${TMP_DIR}/bad-terminal.log" 2>&1; then
  echo 'malformed source terminal was accepted' >&2
  exit 1
fi
[[ ! -e "${BAD_TERMINAL_ROOT}" ]]

# A wrong artifact path is rejected against the configured Yggdrasil path,
# rather than being treated as a second derivative or copied into evidence.
BAD_ARTIFACT_ROOT="${RUNS_ROOT}/bad-artifact"
if run_acceptance "${BAD_ARTIFACT_ROOT}" --artifact-path "${RAW_INPUT}" \
    > "${TMP_DIR}/bad-artifact.log" 2>&1; then
  echo 'raw input path was accepted as the derivative' >&2
  exit 1
fi
[[ ! -e "${BAD_ARTIFACT_ROOT}" ]]

# A semantic failure remains validator-only: no writer, scheduler, checksum
# rehash, source mutation, or fresh terminal status is allowed.
SEMANTIC_FAIL=1
export SEMANTIC_FAIL
BAD_SEMANTIC_ROOT="${RUNS_ROOT}/bad-semantic"
if run_acceptance "${BAD_SEMANTIC_ROOT}" > "${TMP_DIR}/bad-semantic.log" 2>&1; then
  echo 'semantic mismatch was accepted' >&2
  exit 1
fi
unset SEMANTIC_FAIL
[[ ! -e "${BAD_SEMANTIC_ROOT}" ]]
[[ ! -e "${SBATCH_CALLED}" && ! -e "${WRITER_CALLED}" && ! -e "${REHASH_CALLED}" ]]
# A successful but false bounded count must not be accepted or recorded.  This
# exercises the wrapper's parsing of the checked-in validator output rather
# than relying on the approved constant in the report.
SEMANTIC_WRONG_CELLS=1
export SEMANTIC_WRONG_CELLS
BAD_COUNT_ROOT="${RUNS_ROOT}/bad-count"
if run_acceptance "${BAD_COUNT_ROOT}" > "${TMP_DIR}/bad-count.log" 2>&1; then
  echo 'wrong observed cell count was accepted' >&2
  exit 1
fi
unset SEMANTIC_WRONG_CELLS
grep -q 'observed 104 cells; expected 1395601' "${TMP_DIR}/bad-count.log"
[[ ! -e "${BAD_COUNT_ROOT}" ]]
[[ "$(wc -l < "${PYTHON_CALLS}" | tr -d '[:space:]')" == 3 ]]
[[ ! -e "${SBATCH_CALLED}" && ! -e "${WRITER_CALLED}" && ! -e "${REHASH_CALLED}" ]]

for pair in \
  "source.metadata:${SOURCE_RUN_ROOT}/metadata" \
  "source.manifest:${SOURCE_RUN_ROOT}/manifests/source.manifest" \
  "runtime.identity:${SOURCE_RUN_ROOT}/manifests/runtime.identity" \
  "steps.tsv:${STEPS_MANIFEST}" \
  "ownership.tsv:${OWNERSHIP_MANIFEST}" \
  "jobs.tsv:${JOBS_MANIFEST}" \
  "scheduler_ids.tsv:${SCHEDULER_MANIFEST}" \
  "source_watchdog:${WATCHDOG_STATUS}" \
  "source_terminal:${SOURCE_TERMINAL}" \
  "artifact.record:${ARTIFACT_RECORD}" \
  "artifact.md5:${ARTIFACT}.md5" \
  "source_owner:${OWNER_FILE}"; do
  copied="${pair%%:*}"
  original="${pair#*:}"
  cmp -s "${original}" "${BEFORE_DIR}/${copied}"
done

echo 'stage2 derivative acceptance validator: OK'
