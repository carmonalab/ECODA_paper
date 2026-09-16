#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage2-watchdog.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
unset ECODA_SCRATCH_ROOT ECODA_LOGS_DIR ECODA_RUN_ROOT ECODA_RUN_ID LOGS_DIR
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${TMP_DIR}/tmp"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sacct" <<'STUB'
#!/bin/bash
set -euo pipefail
if [[ "${WATCHDOG_MUTATE_MANIFEST:-0}" == "1" &&
      "$*" == *1001* && -n "${WATCHDOG_MUTATION_MARKER:-}" &&
      ! -e "${WATCHDOG_MUTATION_MARKER}" ]]; then
  awk -F '\t' -v outside="${WATCHDOG_OUTSIDE_SCRIPT}" \
    'BEGIN { OFS = "\t" } NR == 1 { $2 = outside } { print }' \
    "${WATCHDOG_MANIFEST}" > "${WATCHDOG_MANIFEST}.mutated"
  mv -f "${WATCHDOG_MANIFEST}.mutated" "${WATCHDOG_MANIFEST}"
  touch "${WATCHDOG_MUTATION_MARKER}"
fi
case "$*" in
  *1001*) printf 'OUT_OF_MEMORY\n' ;;
  *1002*) printf 'COMPLETED\n' ;;
  *) printf 'COMPLETED\n' ;;
esac
STUB
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
if [[ -n "${WATCHDOG_EXPECT_RUN_ROOT:-}" ]]; then
  [[ -f "${WATCHDOG_EXPECT_RUN_ROOT}/manifests/source.manifest" &&
     -f "${WATCHDOG_EXPECT_RUN_ROOT}/manifests/runtime.identity" ]] || {
    echo "watchdog submitted before bound identity manifests" >&2
    exit 97
  }
  case "$*" in
    *"${WATCHDOG_EXPECT_SOURCE_SCRIPT}"*) ;;
    *) echo "watchdog retry used a non-snapshot script" >&2; exit 97 ;;
  esac
  case "$*" in
    *"ECODA_SOURCE_ROOT=${WATCHDOG_EXPECT_SOURCE_ROOT}"*) ;;
    *) echo "watchdog retry omitted source root" >&2; exit 97 ;;
  esac
  case "$*" in
    *"ECODA_RUNTIME_IMAGE=${WATCHDOG_EXPECT_RUNTIME_IMAGE}"*) ;;
    *) echo "watchdog retry omitted runtime image" >&2; exit 97 ;;
  esac
  case "$*" in
    *"ECODA_RUN_ID=${WATCHDOG_EXPECT_RUN_ID}"*) ;;
    *) echo "watchdog retry omitted run identity" >&2; exit 97 ;;
  esac
fi
printf '%s\n' "$*" >> "${CAPTURE}"
printf '1002\n'
STUB
chmod +x "${TMP_DIR}/bin/sacct" "${TMP_DIR}/bin/sbatch"
cat > "${TMP_DIR}/bash_env" <<'BASH_ENV'
wait_job_terminal() {
  local job="${1:-}" state
  state="$(sacct -j "${job}" -X -n -P --format=State 2>/dev/null)" || return 1
  case "${state}" in
    OUT_OF_MEMORY|COMPLETED|FAILED|CANCELLED|TIMEOUT|DEPENDENCY*)
      printf '%s\n' "${state}"
      ;;
    *) return 1 ;;
  esac
}
BASH_ENV
export BASH_ENV="${TMP_DIR}/bash_env"
export HPC_SCRATCH_DIR="${TMP_DIR}/home/scratch/ECODA_paper"
mkdir -p "${HPC_SCRATCH_DIR}" "${HPC_SCRATCH_DIR}/logs" \
  "${TMP_DIR}/source-checkout/src/2_dataset_specific_preprocessing" \
  "${TMP_DIR}/source-checkout/src/utils/bash" \
  "${TMP_DIR}/source-checkout/src/utils/py" "${TMP_DIR}/source-checkout/aux" \
  "${TMP_DIR}/source-snapshots" "${TMP_DIR}/runtime/_ecoda_runtime"
HOST_ENV="${TMP_DIR}/host/.pixi/envs/py-cuda13"
mkdir -p "${HOST_ENV}/bin" "${HOST_ENV}/lib/R"
cat > "${HOST_ENV}/bin/python" <<'STUB'
#!/bin/bash
set -euo pipefail
if [[ "${1:-}" == *"1.7.1_create_alzheimer_donor_assay.py"* ]]; then
  printf '%s\n' "$*" >> "${ALZHEIMER_VALIDATION_CAPTURE:-/dev/null}"
  raw="${HPC_SCRATCH_DIR:?}/Alzheimer/data/SEAAD_Alzheimer.h5ad"
  output="${HPC_SCRATCH_DIR:?}/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad"
  [[ "$(cat "${raw}")" == "immutable Alzheimer raw fixture" ]] || exit 50
  if grep -q "semantically mismatched" "${output}"; then exit 42; fi
  [[ "$*" == *"--input-file ${raw}"* ]] || exit 43
  [[ "$*" == *"--output-file ${output}"* ]] || exit 44
  [[ "$*" == *"--expected-samples 104"* ]] || exit 45
  [[ "$*" == *"--expected-donors 83"* ]] || exit 46
  [[ "$*" == *"--expected-assay-counts 10x3v3=83,10xmultiome=21"* ]] || exit 47
  [[ "$*" == *"--expected-sex-counts female=59,male=45"* ]] || exit 48
  [[ "$*" == *"--require-example-ids"* ]] || exit 49
  [[ "$*" == *"--validate-only"* ]] || exit 50
  printf 'ALZHEIMER_DONOR_ASSAY_VALIDATED=1 cells=1395601 samples=104\n'
fi
exit 0
STUB
printf '#!/bin/bash\nexit 0\n' > "${HOST_ENV}/bin/Rscript"
chmod +x "${HOST_ENV}/bin/python" "${HOST_ENV}/bin/Rscript"

STAGE2_SOURCE_FILES=(
  src/2_dataset_specific_preprocessing/1_submit_hpc.sh
  src/2_dataset_specific_preprocessing/1.submit.sh
  src/2_dataset_specific_preprocessing/1.1.1_subset_gongsharma.py
  src/2_dataset_specific_preprocessing/1.2.1_create_combinedpbmc_dataset.py
  src/2_dataset_specific_preprocessing/1.3.1_prepare_joanito.R
  src/2_dataset_specific_preprocessing/1.4.1_create_kfoury_lowres_ct.R
  src/2_dataset_specific_preprocessing/1.5.1_reconstruct_myocardial_counts.py
  src/2_dataset_specific_preprocessing/1.6.1_fill_bassez_cellsubtype.R
  src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py
  src/utils/bash/sync_status_email.sh
  src/utils/py/artifact_contract.py
  src/utils/py/derived_prerequisite_contract.py
  src/utils/py/h5ad_source_identity.py
  src/2_dataset_specific_preprocessing/stage2_watchdog.sh
)
for source_file in "${STAGE2_SOURCE_FILES[@]}"; do
  cp "${ROOT}/${source_file}" "${TMP_DIR}/source-checkout/${source_file}"
done
cat > "${TMP_DIR}/source-checkout/src/slurm_config.sh" <<'SOURCE_CONFIG'
#!/bin/bash
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
export HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR:?}"
export PYTHONDONTWRITEBYTECODE=1
export LOGS_DIR="${ECODA_LOGS_DIR:-${LOGS_DIR:-${HPC_SCRATCH_DIR}/logs}}"
export ECODA_HOST_ENV_PREFIX="${ECODA_HOST_ENV_PREFIX:?}"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
  export PYTHON_BIN="${ECODA_RUNTIME_PREFIX:?}/bin/python"
  export PIXI_RSCRIPT="${ECODA_RUNTIME_PREFIX}/bin/Rscript --vanilla"
else
  export PYTHON_BIN="${ECODA_HOST_ENV_PREFIX}/bin/python"
  export PIXI_RSCRIPT="${ECODA_HOST_ENV_PREFIX}/bin/Rscript --vanilla"
fi
export PATH="${ECODA_HOST_ENV_PREFIX}/bin:${PATH}"
export LD_LIBRARY_PATH="${ECODA_HOST_ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export RETICULATE_PYTHON="${PYTHON_BIN}"
export ECODA_AUX_ROOT="${ECODA_SOURCE_ROOT%/}/aux"
export SCGATE_DB_BRANCH="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4"
SOURCE_CONFIG
cp "${ROOT}/src/utils/bash/ecoda_runtime.sh" \
  "${TMP_DIR}/source-checkout/src/utils/bash/ecoda_runtime.sh"
cp "${ROOT}/src/utils/bash/ecoda_run_common.sh" \
  "${TMP_DIR}/source-checkout/src/utils/bash/ecoda_run_common.sh"
printf 'config fixture\n' > "${TMP_DIR}/source-checkout/config_helper.R"
printf '{}\n' > "${TMP_DIR}/source-checkout/datasets.json"
printf 'pixi toml fixture\n' > "${TMP_DIR}/source-checkout/pixi.toml"
printf 'pixi lock fixture\n' > "${TMP_DIR}/source-checkout/pixi.lock"
printf 'scGate fixture\n' > "${TMP_DIR}/source-checkout/aux/scGateDB.rds"
printf 'gene blocklist fixture\n' \
  > "${TMP_DIR}/source-checkout/aux/genes.blocklist.rds"
printf 'gene map fixture\n' \
  > "${TMP_DIR}/source-checkout/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
git -C "${TMP_DIR}/source-checkout" init -q
git -C "${TMP_DIR}/source-checkout" config user.email test@example.invalid
git -C "${TMP_DIR}/source-checkout" config user.name stage2-watchdog-test
git -C "${TMP_DIR}/source-checkout" add .
git -C "${TMP_DIR}/source-checkout" commit -qm "stage2 watchdog immutable fixture"
SOURCE_COMMIT="$(git -C "${TMP_DIR}/source-checkout" rev-parse HEAD)"
bash "${ROOT}/src/utils/bash/ecoda_source_snapshot.sh" create \
  --source-root "${TMP_DIR}/source-checkout" \
  --snapshot-parent "${TMP_DIR}/source-snapshots" \
  --commit "${SOURCE_COMMIT}" >/dev/null
SNAPSHOT_ROOT="${TMP_DIR}/source-snapshots/${SOURCE_COMMIT}"
SOURCE_ROOT="${SNAPSHOT_ROOT}/tree"
SOURCE_MANIFEST="${SNAPSHOT_ROOT}/identity/source.manifest"
SOURCE_WATCHDOG="${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/stage2_watchdog.sh"
GENERIC_WORKER="${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.submit.sh"

RUNTIME_ID_DIR="${TMP_DIR}/runtime/_ecoda_runtime/stage2-test"
mkdir -p "${RUNTIME_ID_DIR}"
IMAGE="${RUNTIME_ID_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${IMAGE}.manifest"
printf 'format-2 runtime fixture\n' > "${IMAGE}"
sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}
IMAGE_SHA="$(sha256_file "${IMAGE}")"
TOML_SHA="$(sha256_file "${SOURCE_ROOT}/pixi.toml")"
LOCK_SHA="$(sha256_file "${SOURCE_ROOT}/pixi.lock")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_PATH=${IMAGE}
IMAGE_SHA256=${IMAGE_SHA}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.49.0
APPTAINER_VERSION=1.3.2
IMAGE_BUILD_GIT_REVISION=stage2-watchdog-build-a
IMAGE_PIXI_TOML_SHA256=${TOML_SHA}
IMAGE_PIXI_LOCK_SHA256=${LOCK_SHA}
EOF
MANIFEST_SHA="$(sha256_file "${RUNTIME_MANIFEST}")"
IMAGE_SIZE="$(wc -c < "${IMAGE}" | tr -d '[:space:]')"
MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
chmod 444 "${IMAGE}" "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_ID_DIR}"
cat > "${TMP_DIR}/bin/apptainer" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ "${1:-}" == "inspect" ]] && exit 0
exit 0
STUB
chmod +x "${TMP_DIR}/bin/apptainer"
RUNTIME_IDENTITY_TEMPLATE="${TMP_DIR}/runtime.identity"
cat > "${RUNTIME_IDENTITY_TEMPLATE}" <<EOF
RUNTIME_IMAGE=${IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
RUNTIME_IMAGE_SHA256=${IMAGE_SHA}
RUNTIME_MANIFEST_SHA256=${MANIFEST_SHA}
RUNTIME_IMAGE_SIZE=${IMAGE_SIZE}
RUNTIME_MANIFEST_SIZE=${MANIFEST_SIZE}
IMAGE_PIXI_TOML_SHA256=${TOML_SHA}
IMAGE_PIXI_LOCK_SHA256=${LOCK_SHA}
EOF

export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_PROFILE=stage2
export ECODA_RUNTIME_IMAGE="${IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export ECODA_HOST_ENV_PREFIX="${HOST_ENV}"
export APPTAINER_BIN="${TMP_DIR}/bin/apptainer"
export TMPDIR="${TMP_DIR}/tmp"
RUNS_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs"

write_md5() {
  local path="$1" digest
  digest="$(md5sum "${path}" | cut -d' ' -f1)"
  printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" \
    "$(wc -c < "${path}" | tr -d '[:space:]')" "${path}" > "${path}.md5"
}
make_bound_run() {
  local id="$1" root="${RUNS_ROOT}/${1}"
  mkdir -p "${root}/manifests" "${root}/status" "${root}/logs"
  printf 'STAGE=stage2\nRUN_ID=%s\nSTATE=ACTIVE\n' "${id}" > "${root}/metadata"
  cp "${SOURCE_MANIFEST}" "${root}/manifests/source.manifest"
  cp "${RUNTIME_IDENTITY_TEMPLATE}" "${root}/manifests/runtime.identity"
}
global_owner_dir() {
  HOME="${TMP_DIR}/home" HPC_SCRATCH_DIR="${HPC_SCRATCH_DIR}" \
    ECODA_HOST_ENV_PREFIX="${HOST_ENV}" PATH="${TMP_DIR}/bin:${PATH}" \
    bash -c 'source "$1/src/slurm_config.sh" >/dev/null 2>&1; source "$1/src/utils/bash/ecoda_run_common.sh"; ecoda_artifact_owner_dir "$2"' \
    _ "${ROOT}" "$1"
}
RUN_ID="run"
RUN_ROOT="${RUNS_ROOT}/${RUN_ID}"
KFOURY_OUTPUT="${HPC_SCRATCH_DIR}/Kfoury/data/Kfoury_2021_34719426.rds"
mkdir -p "$(dirname "${KFOURY_OUTPUT}")"
printf 'rds payload\n' > "${KFOURY_OUTPUT}"
make_bound_run "${RUN_ID}"
OWNER_DIR="${HPC_SCRATCH_DIR}/_ecoda_owners/stage2/kfoury_lowres_ct"
mkdir -p "${OWNER_DIR}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=kfoury_lowres_ct\n' \
  "${RUN_ID}" > "${OWNER_DIR}/owner"
MANIFEST="${RUN_ROOT}/manifests/steps.tsv"
JOB_FILE="${RUN_ROOT}/manifests/jobs.tsv"
SOURCE_SCRIPT="${GENERIC_WORKER}"
printf 'kfoury_lowres_ct\t%s\t%s\t-\t%s\n' "${SOURCE_SCRIPT}" \
  "${KFOURY_OUTPUT}" "${OWNER_DIR}" > "${MANIFEST}"
write_md5 "${MANIFEST}"
printf 'kfoury_lowres_ct\t1001\n' > "${JOB_FILE}"
cp "${MANIFEST}" "${RUN_ROOT}/manifests/ownership.tsv"
export WATCHDOG_EXPECT_RUN_ROOT="${RUN_ROOT}"
export WATCHDOG_EXPECT_SOURCE_ROOT="${SOURCE_ROOT}"
export WATCHDOG_EXPECT_SOURCE_SCRIPT="${SOURCE_SCRIPT}"
export WATCHDOG_EXPECT_RUNTIME_IMAGE="${IMAGE}"
export WATCHDOG_EXPECT_RUN_ID="${RUN_ID}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  STAGE2_FORCE=1 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${RUN_ID}" "${MANIFEST}" "${JOB_FILE}" \
  128G 256G shared-cpu 1000
[[ "$(grep '^STATE=' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == "STATE=OK" ]]
[[ -s "${TMP_DIR}/home/scratch/ECODA_paper/Kfoury/data/Kfoury_2021_34719426.rds.md5" ]]
[[ "$(grep '^STATE=' "${OWNER_DIR}/owner")" == "STATE=OK" ]]
[[ "$(grep -c '^SCHEDULER_ID=' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == 2 ]]
[[ "$(grep -c '^SCHEDULER_ID=1001$' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == 1 ]]
[[ "$(grep -c '^SCHEDULER_ID=1002$' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == 1 ]]
case "$(cat "${CAPTURE}")" in *"FORCE_PREPROCESS=1"*) ;; *) echo "OOM retry dropped FORCE_PREPROCESS=1" >&2; exit 1 ;; esac
GLOBAL_OWNER="$(global_owner_dir "${KFOURY_OUTPUT}")"
[[ -f "${GLOBAL_OWNER}/owner" ]]
grep -q "^RUN_ID=${RUN_ID}$" "${GLOBAL_OWNER}/owner"
grep -q '^STATE=OK$' "${GLOBAL_OWNER}/owner"
case "$(cat "${CAPTURE}")" in *"${SOURCE_SCRIPT}"*) ;; *) echo "OOM retry did not use snapshot script" >&2; exit 1 ;; esac
case "$(cat "${CAPTURE}")" in *"ECODA_SOURCE_ROOT=${SOURCE_ROOT}"*) ;; *) echo "OOM retry omitted source root" >&2; exit 1 ;; esac
case "$(cat "${CAPTURE}")" in *"ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}"*) ;; *) echo "OOM retry omitted runtime manifest" >&2; exit 1 ;; esac
case "$(cat "${CAPTURE}")" in *"ECODA_SOURCE_SNAPSHOT_REQUIRED=1"*) ;; *) echo "OOM retry omitted snapshot requirement" >&2; exit 1 ;; esac
case "$(cat "${CAPTURE}")" in *"ECODA_RUNTIME_IDENTITY=${RUN_ROOT}/manifests/runtime.identity"*) ;; *) echo "OOM retry omitted runtime identity" >&2; exit 1 ;; esac
case "$(cat "${CAPTURE}")" in *"ECODA_RUN_ID=${RUN_ID}"*) ;; *) echo "OOM retry omitted run ID" >&2; exit 1 ;; esac
rm -rf "${HPC_SCRATCH_DIR}/_ecoda_owners/artifact"
OUTSIDE_RETRY_RUN_ID="outside_retry"
OUTSIDE_RETRY_ROOT="${RUNS_ROOT}/${OUTSIDE_RETRY_RUN_ID}"
make_bound_run "${OUTSIDE_RETRY_RUN_ID}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=kfoury_lowres_ct\n' \
  "${OUTSIDE_RETRY_RUN_ID}" > "${OWNER_DIR}/owner"
OUTSIDE_RETRY_MANIFEST="${OUTSIDE_RETRY_ROOT}/manifests/steps.tsv"
printf 'kfoury_lowres_ct\t%s\t%s\t-\t%s\n' "${SOURCE_SCRIPT}" \
  "${KFOURY_OUTPUT}" "${OWNER_DIR}" > "${OUTSIDE_RETRY_MANIFEST}"
write_md5 "${OUTSIDE_RETRY_MANIFEST}"
printf 'kfoury_lowres_ct\t1001\n' > "${OUTSIDE_RETRY_ROOT}/manifests/jobs.tsv"
cp "${OUTSIDE_RETRY_MANIFEST}" \
  "${OUTSIDE_RETRY_ROOT}/manifests/ownership.tsv"
OUTSIDE_SCRIPT="${TMP_DIR}/outside-stage2.sh"
printf '#!/bin/bash\nexit 0\n' > "${OUTSIDE_SCRIPT}"
chmod +x "${OUTSIDE_SCRIPT}"
OUTSIDE_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
WATCHDOG_MUTATION_MARKER="${TMP_DIR}/manifest.mutated"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  WATCHDOG_EXPECT_RUN_ROOT="" WATCHDOG_MUTATE_MANIFEST=1 \
  WATCHDOG_MANIFEST="${OUTSIDE_RETRY_MANIFEST}" \
  WATCHDOG_OUTSIDE_SCRIPT="${OUTSIDE_SCRIPT}" \
  WATCHDOG_MUTATION_MARKER="${WATCHDOG_MUTATION_MARKER}" \
  STAGE2_FORCE=1 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${OUTSIDE_RETRY_RUN_ID}" \
  "${OUTSIDE_RETRY_MANIFEST}" "${OUTSIDE_RETRY_ROOT}/manifests/jobs.tsv" \
  128G 256G shared-cpu 1000 > "${TMP_DIR}/outside-retry.watchdog.log" 2>&1; then
  echo "outside-root OOM retry unexpectedly succeeded" >&2
  exit 1
fi
grep -Eq '^ERROR:[[:space:]]*[^[:space:]]' \
  "${TMP_DIR}/outside-retry.watchdog.log"
[[ "$(grep '^STATE=' "${OUTSIDE_RETRY_ROOT}/status/watchdog")" == "STATE=FAIL" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${OUTSIDE_CALLS_BEFORE}" ]]

OTHER_RUN_ID="other_active_run"
OTHER_ROOT="${RUNS_ROOT}/${OTHER_RUN_ID}"
make_bound_run "${OTHER_RUN_ID}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=kfoury_lowres_ct\n' \
  "${OTHER_RUN_ID}" > "${OWNER_DIR}/owner"
OTHER_MANIFEST="${OTHER_ROOT}/manifests/steps.tsv"
printf 'kfoury_lowres_ct\t%s\t%s\t-\t%s\n' "${SOURCE_SCRIPT}" \
  "${KFOURY_OUTPUT}" "${OWNER_DIR}" > "${OTHER_MANIFEST}"
write_md5 "${OTHER_MANIFEST}"
printf 'kfoury_lowres_ct\t1001\n' > "${OTHER_ROOT}/manifests/jobs.tsv"
cp "${OTHER_MANIFEST}" "${OTHER_ROOT}/manifests/ownership.tsv"
printf 'RUN_ID=unrelated_active_run\nSTAGE=stage2\nPATH=%s\nSTATE=ACTIVE\nPID=1\n' \
  "${KFOURY_OUTPUT}" > "${GLOBAL_OWNER}/owner"
OTHER_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  WATCHDOG_EXPECT_RUN_ROOT="" STAGE2_FORCE=1 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${OTHER_RUN_ID}" "${OTHER_MANIFEST}" \
  "${OTHER_ROOT}/manifests/jobs.tsv" 128G 256G shared-cpu 1000 \
  > "${TMP_DIR}/other-active.watchdog.log" 2>&1; then
  echo "another run's active global owner was bypassed" >&2
  exit 1
fi
grep -q "active global artifact owner" \
  "${TMP_DIR}/other-active.watchdog.log"
[[ "$(grep '^STATE=' "${OTHER_ROOT}/status/watchdog")" == "STATE=FAIL" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${OTHER_CALLS_BEFORE}" ]]
BASSEZ_RUN_ID="bassez_run"
BASSEZ_ROOT="${RUNS_ROOT}/${BASSEZ_RUN_ID}"
BASSEZ_OUTPUT="${HPC_SCRATCH_DIR}/Bassez/data/BassezA_2021_33958794whole.rds"
BASSEZ_OWNER="${HPC_SCRATCH_DIR}/_ecoda_owners/stage2/bassez_cellsubtype"
mkdir -p "$(dirname "${BASSEZ_OUTPUT}")" "${BASSEZ_OWNER}"
make_bound_run "${BASSEZ_RUN_ID}"
printf 'rds payload\n' > "${BASSEZ_OUTPUT}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=bassez_cellsubtype\n' \
  "${BASSEZ_RUN_ID}" > "${BASSEZ_OWNER}/owner"
BASSEZ_MANIFEST="${BASSEZ_ROOT}/manifests/steps.tsv"
BASSEZ_JOB_FILE="${BASSEZ_ROOT}/manifests/jobs.tsv"
printf 'bassez_cellsubtype\t%s\t%s\t-\t%s\n' \
  "${GENERIC_WORKER}" \
  "${BASSEZ_OUTPUT}" "${BASSEZ_OWNER}" > "${BASSEZ_MANIFEST}"
write_md5 "${BASSEZ_MANIFEST}"
printf 'bassez_cellsubtype\t2001\n' > "${BASSEZ_JOB_FILE}"
cp "${BASSEZ_MANIFEST}" "${BASSEZ_ROOT}/manifests/ownership.tsv"
BASSEZ_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid STAGE2_FORCE=0 \
  STAGE2_WATCHDOG_MAX_POLLS=1 bash "${SOURCE_WATCHDOG}" \
  "${BASSEZ_RUN_ID}" "${BASSEZ_MANIFEST}" "${BASSEZ_JOB_FILE}" 128G 256G shared-cpu 1000
[[ "$(grep '^STATE=' "${BASSEZ_ROOT}/status/watchdog")" == "STATE=OK" ]]
[[ "$(grep '^STATE=' "${BASSEZ_OWNER}/owner")" == "STATE=OK" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${BASSEZ_CALLS_BEFORE}" ]]
rm -rf "${HPC_SCRATCH_DIR}/_ecoda_owners/artifact"
BAD_BASSEZ_DIGEST="00000000000000000000000000000000"
BASSEZ_SIZE="$(wc -c < "${BASSEZ_OUTPUT}" | tr -d '[:space:]')"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${BAD_BASSEZ_DIGEST}" \
  "${BASSEZ_SIZE}" "${BASSEZ_OUTPUT}" > "${BASSEZ_OUTPUT}.md5"
BASSEZ_BAD_SIDECAR="$(cat "${BASSEZ_OUTPUT}.md5")"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid STAGE2_FORCE=0 \
  STAGE2_WATCHDOG_MAX_POLLS=1 bash "${SOURCE_WATCHDOG}" \
  "${BASSEZ_RUN_ID}" "${BASSEZ_MANIFEST}" "${BASSEZ_JOB_FILE}" 128G 256G shared-cpu 1000 \
  > "${TMP_DIR}/invalid-checksum.watchdog.log" 2>&1; then
  echo "Stage 2 rewrote an invalid existing checksum" >&2
  exit 1
fi
[[ "$(cat "${BASSEZ_OUTPUT}.md5")" == "${BASSEZ_BAD_SIDECAR}" ]]
rm -rf "${HPC_SCRATCH_DIR}/_ecoda_owners/artifact"

INVALID_RUN_ID="invalid_bassez_run"
INVALID_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${INVALID_RUN_ID}"
mkdir -p "${INVALID_ROOT}/status" "${INVALID_ROOT}/manifests"
printf 'STAGE=stage2\nRUN_ID=%s\nSTATE=ACTIVE\n' "${INVALID_RUN_ID}" > "${INVALID_ROOT}/metadata"
cp "${SOURCE_MANIFEST}" "${INVALID_ROOT}/manifests/source.manifest"
cp "${RUNTIME_IDENTITY_TEMPLATE}" \
  "${INVALID_ROOT}/manifests/runtime.identity"
INVALID_STEP="basse""x_cellsubtype"
INVALID_MANIFEST="${INVALID_ROOT}/manifests/steps.tsv"
INVALID_JOB_FILE="${INVALID_ROOT}/manifests/jobs.tsv"
sed "s/^bassez_cellsubtype/${INVALID_STEP}/" "${BASSEZ_MANIFEST}" > "${INVALID_MANIFEST}"
sed "s/^bassez_cellsubtype/${INVALID_STEP}/" "${BASSEZ_JOB_FILE}" > "${INVALID_JOB_FILE}"
INVALID_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
write_md5 "${INVALID_MANIFEST}"
cp "${INVALID_MANIFEST}" "${INVALID_ROOT}/manifests/ownership.tsv"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid STAGE2_FORCE=0 \
  STAGE2_WATCHDOG_MAX_POLLS=1 bash "${SOURCE_WATCHDOG}" \
  "${INVALID_RUN_ID}" "${INVALID_MANIFEST}" "${INVALID_JOB_FILE}" 128G 256G shared-cpu 1000 \
  > "${TMP_DIR}/invalid.watchdog.log" 2>&1; then
  echo "legacy Bassez manifest spelling unexpectedly accepted" >&2
  exit 1
fi
[[ "$(grep '^STATE=' "${INVALID_ROOT}/status/watchdog")" == "STATE=FAIL" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${INVALID_CALLS_BEFORE}" ]]
LEGACY_RUN_ID="legacy_unpinned_retry"
LEGACY_ROOT="${RUNS_ROOT}/${LEGACY_RUN_ID}"
mkdir -p "${LEGACY_ROOT}/manifests" "${LEGACY_ROOT}/status"
printf 'STAGE=stage2\nRUN_ID=%s\nSTATE=ACTIVE\n' "${LEGACY_RUN_ID}" \
  > "${LEGACY_ROOT}/metadata"
LEGACY_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  STAGE2_FORCE=1 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${LEGACY_RUN_ID}" \
  "${LEGACY_ROOT}/manifests/steps.tsv" \
  "${LEGACY_ROOT}/manifests/jobs.tsv" 128G 256G shared-cpu 1000 \
  > "${TMP_DIR}/legacy-retry.watchdog.log" 2>&1; then
  echo "legacy unpinned retry unexpectedly succeeded" >&2
  exit 1
fi
grep -q '^REASON=legacy_source_unpinned$' "${LEGACY_ROOT}/status/watchdog"
[[ "$(grep '^STATE=' "${LEGACY_ROOT}/status/watchdog")" == "STATE=FAIL" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${LEGACY_CALLS_BEFORE}" ]]
[[ ! -e "${INVALID_ROOT}/manifests/jobs.retry_1.tsv" ]]

MYOCARDIAL_RUN_ID="myocardial_watchdog"
MYOCARDIAL_ROOT="${RUNS_ROOT}/${MYOCARDIAL_RUN_ID}"
MYOCARDIAL_OUTPUT="${HPC_SCRATCH_DIR}/Myocardial_infarction/data/Myocardial_Infarc_2.h5ad"
mkdir -p "$(dirname "${MYOCARDIAL_OUTPUT}")"
printf 'myocardial H5AD fixture\n' > "${MYOCARDIAL_OUTPUT}"
make_bound_run "${MYOCARDIAL_RUN_ID}"
MYOCARDIAL_OWNER="${HPC_SCRATCH_DIR}/_ecoda_owners/stage2/myocardial_counts"
mkdir -p "${MYOCARDIAL_OWNER}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=myocardial_counts\n' \
  "${MYOCARDIAL_RUN_ID}" > "${MYOCARDIAL_OWNER}/owner"
MYOCARDIAL_MANIFEST="${MYOCARDIAL_ROOT}/manifests/steps.tsv"
MYOCARDIAL_JOB_FILE="${MYOCARDIAL_ROOT}/manifests/jobs.tsv"
printf 'myocardial_counts\t%s\t%s\t-\t%s\n' \
  "${GENERIC_WORKER}" \
  "${MYOCARDIAL_OUTPUT}" "${MYOCARDIAL_OWNER}" > "${MYOCARDIAL_MANIFEST}"
write_md5 "${MYOCARDIAL_MANIFEST}"
printf 'myocardial_counts\t2501\n' > "${MYOCARDIAL_JOB_FILE}"
cp "${MYOCARDIAL_MANIFEST}" "${MYOCARDIAL_ROOT}/manifests/ownership.tsv"
export WATCHDOG_EXPECT_RUN_ROOT="${MYOCARDIAL_ROOT}"
export WATCHDOG_EXPECT_SOURCE_SCRIPT="${GENERIC_WORKER}"
export WATCHDOG_EXPECT_RUN_ID="${MYOCARDIAL_RUN_ID}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  STAGE2_FORCE=0 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${MYOCARDIAL_RUN_ID}" \
  "${MYOCARDIAL_MANIFEST}" "${MYOCARDIAL_JOB_FILE}" 128G 256G shared-cpu 1000
[[ "$(grep '^STATE=' "${MYOCARDIAL_ROOT}/status/watchdog")" == "STATE=OK" ]]
[[ "$(grep '^STATE=' "${MYOCARDIAL_OWNER}/owner")" == "STATE=OK" ]]
[[ -s "${MYOCARDIAL_OUTPUT}.md5" ]]
rm -rf "${HPC_SCRATCH_DIR}/_ecoda_owners/artifact"

ALZHEIMER_RAW="${HPC_SCRATCH_DIR}/Alzheimer/data/SEAAD_Alzheimer.h5ad"
ALZHEIMER_OUTPUT="${HPC_SCRATCH_DIR}/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad"
mkdir -p "$(dirname "${ALZHEIMER_RAW}")"
printf 'immutable Alzheimer raw fixture\n' > "${ALZHEIMER_RAW}"
printf 'valid Alzheimer derivative fixture\n' > "${ALZHEIMER_OUTPUT}"
export ALZHEIMER_VALIDATION_CAPTURE="${TMP_DIR}/alzheimer-validation.calls"
: > "${ALZHEIMER_VALIDATION_CAPTURE}"

ALZHEIMER_RUN_ID="alzheimer_valid"
ALZHEIMER_ROOT="${RUNS_ROOT}/${ALZHEIMER_RUN_ID}"
make_bound_run "${ALZHEIMER_RUN_ID}"
ALZHEIMER_OWNER="${HPC_SCRATCH_DIR}/_ecoda_owners/stage2/alzheimer_donor_assay"
mkdir -p "${ALZHEIMER_OWNER}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=alzheimer_donor_assay\n' \
  "${ALZHEIMER_RUN_ID}" > "${ALZHEIMER_OWNER}/owner"
ALZHEIMER_MANIFEST="${ALZHEIMER_ROOT}/manifests/steps.tsv"
ALZHEIMER_JOB_FILE="${ALZHEIMER_ROOT}/manifests/jobs.tsv"
printf 'alzheimer_donor_assay\t%s\t%s\t-\t%s\n' \
  "${GENERIC_WORKER}" \
  "${ALZHEIMER_OUTPUT}" "${ALZHEIMER_OWNER}" > "${ALZHEIMER_MANIFEST}"
write_md5 "${ALZHEIMER_MANIFEST}"
printf 'alzheimer_donor_assay\t3001\n' > "${ALZHEIMER_JOB_FILE}"
cp "${ALZHEIMER_MANIFEST}" "${ALZHEIMER_ROOT}/manifests/ownership.tsv"
export WATCHDOG_EXPECT_RUN_ROOT="${ALZHEIMER_ROOT}"
export WATCHDOG_EXPECT_SOURCE_ROOT="${SOURCE_ROOT}"
export WATCHDOG_EXPECT_SOURCE_SCRIPT="${GENERIC_WORKER}"
export WATCHDOG_EXPECT_RUN_ID="${ALZHEIMER_RUN_ID}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  STAGE2_FORCE=0 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${ALZHEIMER_RUN_ID}" \
  "${ALZHEIMER_MANIFEST}" "${ALZHEIMER_JOB_FILE}" 128G 256G shared-cpu 1000
[[ "$(grep '^STATE=' "${ALZHEIMER_ROOT}/status/watchdog")" == "STATE=OK" ]]
[[ "$(grep '^STATE=' "${ALZHEIMER_OWNER}/owner")" == "STATE=OK" ]]
[[ -s "${ALZHEIMER_OUTPUT}.md5" ]]
[[ "$(wc -l < "${ALZHEIMER_VALIDATION_CAPTURE}" | tr -d '[:space:]')" == 1 ]]
ALZHEIMER_VALIDATION_CALL="$(cat "${ALZHEIMER_VALIDATION_CAPTURE}")"
case "${ALZHEIMER_VALIDATION_CALL}" in
  *"${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py"*) ;;
  *) echo "watchdog skipped the immutable Alzheimer worker" >&2; exit 1 ;;
esac
case "${ALZHEIMER_VALIDATION_CALL}" in
  *"--input-file ${ALZHEIMER_RAW}"*"--output-file ${ALZHEIMER_OUTPUT}"*) ;;
  *) echo "watchdog passed the wrong Alzheimer source/output paths" >&2; exit 1 ;;
esac
case "${ALZHEIMER_VALIDATION_CALL}" in
  *"--expected-samples 104"*"--expected-donors 83"*"--expected-assay-counts 10x3v3=83,10xmultiome=21"*"--expected-sex-counts female=59,male=45"*"--require-example-ids"*"--validate-only"*) ;;
  *) echo "watchdog omitted Alzheimer semantic expectations" >&2; exit 1 ;;
esac


# A broad step spelling must fail in the owner/step contract before any
# scheduler action, even when its row otherwise resembles the explicit hook.
ALZHEIMER_BROAD_RUN_ID="alzheimer_broad"
ALZHEIMER_BROAD_ROOT="${RUNS_ROOT}/${ALZHEIMER_BROAD_RUN_ID}"
make_bound_run "${ALZHEIMER_BROAD_RUN_ID}"
ALZHEIMER_BROAD_STEP="alzheimer_donor_assay,joanito"
ALZHEIMER_BROAD_OWNER="${HPC_SCRATCH_DIR}/_ecoda_owners/stage2/alzheimer_donor_assay_joanito"
mkdir -p "${ALZHEIMER_BROAD_OWNER}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=%s\n' \
  "${ALZHEIMER_BROAD_RUN_ID}" "${ALZHEIMER_BROAD_STEP}" \
  > "${ALZHEIMER_BROAD_OWNER}/owner"
ALZHEIMER_BROAD_MANIFEST="${ALZHEIMER_BROAD_ROOT}/manifests/steps.tsv"
ALZHEIMER_BROAD_JOB_FILE="${ALZHEIMER_BROAD_ROOT}/manifests/jobs.tsv"
printf '%s\t%s\t%s\t-\t%s\n' "${ALZHEIMER_BROAD_STEP}" \
  "${WATCHDOG_EXPECT_SOURCE_SCRIPT}" "${ALZHEIMER_OUTPUT}" \
  "${ALZHEIMER_BROAD_OWNER}" > "${ALZHEIMER_BROAD_MANIFEST}"
write_md5 "${ALZHEIMER_BROAD_MANIFEST}"
printf '%s\t3002\n' "${ALZHEIMER_BROAD_STEP}" > "${ALZHEIMER_BROAD_JOB_FILE}"
cp "${ALZHEIMER_BROAD_MANIFEST}" \
  "${ALZHEIMER_BROAD_ROOT}/manifests/ownership.tsv"
BROAD_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  STAGE2_FORCE=0 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${ALZHEIMER_BROAD_RUN_ID}" \
  "${ALZHEIMER_BROAD_MANIFEST}" "${ALZHEIMER_BROAD_JOB_FILE}" \
  128G 256G shared-cpu 1000 > "${TMP_DIR}/alzheimer-broad.watchdog.log" 2>&1; then
  echo "broad Alzheimer watchdog step unexpectedly succeeded" >&2
  exit 1
fi
[[ "$(grep '^STATE=' "${ALZHEIMER_BROAD_ROOT}/status/watchdog")" == "STATE=FAIL" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${BROAD_CALLS_BEFORE}" ]]

# Reuse the canonical output path with a fresh run, but make the immutable
# worker report a semantic mismatch.  The watchdog must not write a checksum
# or reach terminal OK after that failure.
rm -rf "${HPC_SCRATCH_DIR}/_ecoda_owners/artifact" "${ALZHEIMER_OWNER}"
chmod u+w "${ALZHEIMER_OUTPUT}"
rm -f "${ALZHEIMER_OUTPUT}.md5"
printf 'semantically mismatched Alzheimer derivative fixture\n' > "${ALZHEIMER_OUTPUT}"
ALZHEIMER_BAD_RUN_ID="alzheimer_semantic_failure"
ALZHEIMER_BAD_ROOT="${RUNS_ROOT}/${ALZHEIMER_BAD_RUN_ID}"
make_bound_run "${ALZHEIMER_BAD_RUN_ID}"
mkdir -p "${ALZHEIMER_OWNER}"
printf 'RUN_ID=%s\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=alzheimer_donor_assay\n' \
  "${ALZHEIMER_BAD_RUN_ID}" > "${ALZHEIMER_OWNER}/owner"
ALZHEIMER_BAD_MANIFEST="${ALZHEIMER_BAD_ROOT}/manifests/steps.tsv"
ALZHEIMER_BAD_JOB_FILE="${ALZHEIMER_BAD_ROOT}/manifests/jobs.tsv"
printf 'alzheimer_donor_assay\t%s\t%s\t-\t%s\n' \
  "${WATCHDOG_EXPECT_SOURCE_SCRIPT}" "${ALZHEIMER_OUTPUT}" \
  "${ALZHEIMER_OWNER}" > "${ALZHEIMER_BAD_MANIFEST}"
write_md5 "${ALZHEIMER_BAD_MANIFEST}"
printf 'alzheimer_donor_assay\t3003\n' > "${ALZHEIMER_BAD_JOB_FILE}"
cp "${ALZHEIMER_BAD_MANIFEST}" \
  "${ALZHEIMER_BAD_ROOT}/manifests/ownership.tsv"
export WATCHDOG_EXPECT_RUN_ROOT="${ALZHEIMER_BAD_ROOT}"
export WATCHDOG_EXPECT_RUN_ID="${ALZHEIMER_BAD_RUN_ID}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid \
  STAGE2_FORCE=0 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${SOURCE_WATCHDOG}" "${ALZHEIMER_BAD_RUN_ID}" \
  "${ALZHEIMER_BAD_MANIFEST}" "${ALZHEIMER_BAD_JOB_FILE}" \
  128G 256G shared-cpu 1000 > "${TMP_DIR}/alzheimer-semantic.watchdog.log" 2>&1; then
  echo "semantically invalid Alzheimer derivative unexpectedly succeeded" >&2
  exit 1
fi
[[ "$(grep '^STATE=' "${ALZHEIMER_BAD_ROOT}/status/watchdog")" == "STATE=FAIL" ]]
[[ "$(grep '^STATE=' "${ALZHEIMER_OWNER}/owner")" == "STATE=FAIL" ]]
[[ ! -e "${ALZHEIMER_OUTPUT}.md5" ]]
[[ "$(wc -l < "${ALZHEIMER_VALIDATION_CAPTURE}" | tr -d '[:space:]')" == 2 ]]

echo "stage2 watchdog: OK"
