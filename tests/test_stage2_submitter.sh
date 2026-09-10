#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage2.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
unset BASH_ENV ECODA_SCRATCH_ROOT ECODA_LOGS_DIR ECODA_RUN_ROOT ECODA_RUN_ID LOGS_DIR
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${TMP_DIR}/tmp"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
if [[ -n "${EXPECT_RUN_ROOT:-}" ]]; then
  [[ -f "${EXPECT_RUN_ROOT}/manifests/source.manifest" &&
     ! -L "${EXPECT_RUN_ROOT}/manifests/source.manifest" ]] || {
    echo "run source manifest was not installed before sbatch" >&2
    exit 97
  }
  [[ -f "${EXPECT_RUN_ROOT}/manifests/runtime.identity" &&
     ! -L "${EXPECT_RUN_ROOT}/manifests/runtime.identity" ]] || {
    echo "run runtime.identity was not installed before sbatch" >&2
    exit 97
  }
  case "$*" in
    *"ECODA_SOURCE_ROOT=${EXPECT_SOURCE_ROOT}"*) ;;
    *) echo "scheduler command omitted immutable source root" >&2; exit 97 ;;
  esac
  case "$*" in
    *"ECODA_RUNTIME_IMAGE=${EXPECT_RUNTIME_IMAGE}"*) ;;
    *) echo "scheduler command omitted immutable runtime image" >&2; exit 97 ;;
  esac
  case "$*" in
    *"ECODA_RUN_ID=${EXPECT_RUN_ID}"*) ;;
    *) echo "scheduler command omitted run identity" >&2; exit 97 ;;
  esac
fi
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '71000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
HOST_ENV="${TMP_DIR}/host/.pixi/envs/py-cuda13"
REAL_PIXI="$(command -v pixi || true)"
export REAL_PIXI
mkdir -p "${HOST_ENV}/bin" "${HOST_ENV}/lib" \
  "${TMP_DIR}/source-checkout/src/2_dataset_specific_preprocessing" \
  "${TMP_DIR}/source-checkout/src/utils/bash" \
  "${TMP_DIR}/source-checkout/src/utils/py" \
  "${TMP_DIR}/source-checkout/aux" "${TMP_DIR}/source-snapshots" \
  "${TMP_DIR}/runtime/_ecoda_runtime"
cat > "${HOST_ENV}/bin/python" <<EOF
#!/bin/bash
if [[ "\${1:-}" == *derived_prerequisite_contract.py && -n "\${REAL_PIXI:-}" ]]; then
  cd "${ROOT}"
  exec "\${REAL_PIXI}" run python "\$@"
fi
target="\${HOOK_FORCE_FILE:-${TMP_DIR}/force.args}"
printf '%s\n' "\$*" > "\${target}"
EOF
printf '#!/bin/bash\nexit 0\n' > "${HOST_ENV}/bin/Rscript"
chmod +x "${HOST_ENV}/bin/python" "${HOST_ENV}/bin/Rscript"

STAGE2_SOURCE_FILES=(
  src/2_dataset_specific_preprocessing/1_submit_hpc.sh
  src/2_dataset_specific_preprocessing/1.1_submit_gongsharma.sh
  src/2_dataset_specific_preprocessing/1.2_submit_combinedpbmc.sh
  src/2_dataset_specific_preprocessing/1.3_submit_joanito.sh
  src/2_dataset_specific_preprocessing/1.4_submit_kfoury_lowres_ct.sh
  src/2_dataset_specific_preprocessing/1.5_submit_myocardial.sh
  src/2_dataset_specific_preprocessing/1.6_submit_bassez.sh
  src/utils/py/derived_prerequisite_contract.py
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
export SLURM_PARTITION=shared-cpu
export MAX_NUM_CHUNKS_PARALLEL=1000
SOURCE_CONFIG
cp "${ROOT}/src/utils/bash/ecoda_runtime.sh" \
  "${TMP_DIR}/source-checkout/src/utils/bash/ecoda_runtime.sh"
cp "${ROOT}/src/utils/bash/ecoda_run_common.sh" \
  "${TMP_DIR}/source-checkout/src/utils/bash/ecoda_run_common.sh"
cp "${ROOT}/src/utils/bash/sync_status_email.sh" \
  "${TMP_DIR}/source-checkout/src/utils/bash/sync_status_email.sh"
printf 'config fixture\n' > "${TMP_DIR}/source-checkout/config_helper.R"
cp "${ROOT}/datasets.json" "${TMP_DIR}/source-checkout/datasets.json"
printf 'pixi toml fixture\n' > "${TMP_DIR}/source-checkout/pixi.toml"
printf 'pixi lock fixture\n' > "${TMP_DIR}/source-checkout/pixi.lock"
printf 'scGate fixture\n' > "${TMP_DIR}/source-checkout/aux/scGateDB.rds"
printf 'gene blocklist fixture\n' > \
  "${TMP_DIR}/source-checkout/aux/genes.blocklist.rds"
printf 'gene map fixture\n' > \
  "${TMP_DIR}/source-checkout/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
git -C "${TMP_DIR}/source-checkout" init -q
git -C "${TMP_DIR}/source-checkout" config user.email test@example.invalid
git -C "${TMP_DIR}/source-checkout" config user.name stage2-test
git -C "${TMP_DIR}/source-checkout" add .
git -C "${TMP_DIR}/source-checkout" commit -qm "stage2 immutable fixture"
SOURCE_COMMIT="$(git -C "${TMP_DIR}/source-checkout" rev-parse HEAD)"
bash "${ROOT}/src/utils/bash/ecoda_source_snapshot.sh" create \
  --source-root "${TMP_DIR}/source-checkout" \
  --snapshot-parent "${TMP_DIR}/source-snapshots" \
  --commit "${SOURCE_COMMIT}" >/dev/null
SNAPSHOT_ROOT="${TMP_DIR}/source-snapshots/${SOURCE_COMMIT}"
SOURCE_ROOT="${SNAPSHOT_ROOT}/tree"
SOURCE_MANIFEST="${SNAPSHOT_ROOT}/identity/source.manifest"

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
HOST_PYTHON_SHA256="$(sha256_file "${HOST_ENV}/bin/python")"
HOST_RSCRIPT_SHA256="$(sha256_file "${HOST_ENV}/bin/Rscript")"
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
IMAGE_BUILD_GIT_REVISION=stage2-build-a
IMAGE_PIXI_TOML_SHA256=${TOML_SHA}
IMAGE_PIXI_LOCK_SHA256=${LOCK_SHA}
EOF
IMAGE_SIZE="$(wc -c < "${IMAGE}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SHA="$(sha256_file "${RUNTIME_MANIFEST}")"
chmod 444 "${IMAGE}" "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_ID_DIR}"
cat > "${TMP_DIR}/bin/apptainer" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ "${1:-}" == "inspect" ]] && exit 0
exit 0
STUB
chmod +x "${TMP_DIR}/bin/apptainer"
REAL_PIXI="$(command -v pixi || true)"
if [[ -n "${REAL_PIXI}" ]]; then
  cat > "${TMP_DIR}/bin/pixi" <<EOF
#!/bin/bash
set -euo pipefail
cd "${ROOT}"
exec "${REAL_PIXI}" "\$@"
EOF
  chmod +x "${TMP_DIR}/bin/pixi"
fi

export HPC_SCRATCH_DIR="${TMP_DIR}/home/scratch/ECODA_paper"
mkdir -p "${HPC_SCRATCH_DIR}/logs"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_MODE=apptainer
export ECODA_RUNTIME_IN_CONTAINER=0
export ECODA_RUNTIME_PROFILE=stage2
export ECODA_RUNTIME_IMAGE="${IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export ECODA_HOST_ENV_PREFIX="${HOST_ENV}"
export ECODA_HOST_PYTHON_SHA256="${HOST_PYTHON_SHA256}"
export ECODA_HOST_RSCRIPT_SHA256="${HOST_RSCRIPT_SHA256}"
export APPTAINER_BIN="${TMP_DIR}/bin/apptainer"
export TMPDIR="${TMP_DIR}/tmp"
export EXPECT_SOURCE_ROOT="${SOURCE_ROOT}"
export EXPECT_RUNTIME_IMAGE="${IMAGE}"
RUNS_ROOT="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs"
INVALID_STEP="basse""x_cellsubtype"
: > "${CAPTURE}"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --datasets Bassez --steps "${INVALID_STEP}" > "${TMP_DIR}/invalid.submitter.log" 2>&1; then
  echo "legacy Bassez spelling unexpectedly accepted" >&2
  exit 1
fi
grep -q "ERROR: unknown Stage 2 step '${INVALID_STEP}'." "${TMP_DIR}/invalid.submitter.log"
[[ ! -s "${CAPTURE}" ]]
[[ ! -e "${RUNS_ROOT}" ]]

SUBMIT_RUN_ID="stage2_submitter_combined"
export ECODA_RUN_ID="${SUBMIT_RUN_ID}"
export EXPECT_RUN_ID="${SUBMIT_RUN_ID}"
export EXPECT_RUN_ROOT="${RUNS_ROOT}/${SUBMIT_RUN_ID}"
OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  STAGE2_SUBMITTER_TEST=1 bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --datasets CombinedPBMC,_debug --steps combinedpbmc,joanito --force
)"
RUN_ID="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^STAGE2_RUN_ID=//p')"
[[ -n "${RUN_ID}" ]]
MANIFEST="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}/manifests/steps.tsv"
[[ "$(wc -l < "${MANIFEST}" | tr -d '[:space:]')" == 3 ]]
SCHEDULER_MANIFEST="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}/manifests/scheduler_ids.tsv"
[[ "$(wc -l < "${SCHEDULER_MANIFEST}" | tr -d '[:space:]')" == 4 ]]
SUBMIT_ROOT="${RUNS_ROOT}/${RUN_ID}"
[[ "${RUN_ID}" == "${SUBMIT_RUN_ID}" ]]
[[ -f "${SUBMIT_ROOT}/manifests/source.manifest" &&
   -f "${SUBMIT_ROOT}/manifests/runtime.identity" ]]
cmp -s "${SOURCE_MANIFEST}" "${SUBMIT_ROOT}/manifests/source.manifest"
[[ "$(wc -l < "${SUBMIT_ROOT}/manifests/runtime.identity" | tr -d '[:space:]')" == 8 ]]
grep -q "^RUNTIME_IMAGE=${IMAGE}$" \
  "${SUBMIT_ROOT}/manifests/runtime.identity"
grep -q "^RUNTIME_MANIFEST=${RUNTIME_MANIFEST}$" \
  "${SUBMIT_ROOT}/manifests/runtime.identity"
grep -q "^SOURCE_ROOT=${SOURCE_ROOT}$" "${SUBMIT_ROOT}/metadata"
grep -q "^RUNTIME_IDENTITY=${SUBMIT_ROOT}/manifests/runtime.identity$" \
  "${SUBMIT_ROOT}/metadata"
case "${OUTPUT}" in *"STAGE2_RUN_ID=${SUBMIT_RUN_ID}"*) ;; *) exit 1 ;; esac
CALLS="$(cat "${CAPTURE}")"
case "${CALLS}" in *"ECODA_RUNTIME_MODE=apptainer"*"ECODA_RUNTIME_PROFILE=stage2"*) ;; *) echo "Stage 2 runtime export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_ROOT=${SOURCE_ROOT}"*) ;; *) echo "snapshot source root missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_MANIFEST=${SOURCE_MANIFEST}"*) ;; *) echo "snapshot source manifest missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_SNAPSHOT_REQUIRED=1"*) ;; *) echo "snapshot-required export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_SOURCE_MANIFEST_RUN=${SUBMIT_ROOT}/manifests/source.manifest"*) ;; *) echo "run source manifest missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_IMAGE=${IMAGE}"*"ECODA_RUNTIME_MANIFEST=${RUNTIME_MANIFEST}"*) ;; *) echo "runtime paths missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_IMAGE_SHA256=${IMAGE_SHA}"*"ECODA_RUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA}"*) ;; *) echo "runtime digest identity missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUNTIME_IDENTITY=${SUBMIT_ROOT}/manifests/runtime.identity"*) ;; *) echo "runtime identity path missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"ECODA_RUN_ID=${SUBMIT_RUN_ID}"*) ;; *) echo "run ID missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.1_submit_gongsharma.sh"*) ;; *) echo "snapshot GongSharma script missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.2_submit_combinedpbmc.sh"*) ;; *) echo "snapshot CombinedPBMC script missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"${ROOT}/src/2_dataset_specific_preprocessing/1.1_submit_gongsharma.sh"*) echo "mutable GongSharma script submitted" >&2; exit 1 ;; esac
case "${CALLS}" in *"--dependency=afterok:710001"*) ;; *) echo "CombinedPBMC cap dependency missing" >&2; exit 1 ;; esac
JOANITO_CALL="$(sed -n '3p' "${CAPTURE}")"
case "${JOANITO_CALL}" in *"--dependency="*) echo "Joanito was artificially serialized" >&2; exit 1 ;; esac
WATCHDOG_CALL="$(sed -n '4p' "${CAPTURE}")"
case "${WATCHDOG_CALL}" in *"--dependency=afterany:710001:710002:710003"*) ;; *) echo "aggregate watchdog dependency missing" >&2; exit 1 ;; esac
case "${WATCHDOG_CALL}" in *"--mem=64G"*) ;; *) echo "watchdog default memory missing" >&2; exit 1 ;; esac
BASSEZ_SUBMIT_RUN_ID="stage2_submitter_bassez"
export ECODA_RUN_ID="${BASSEZ_SUBMIT_RUN_ID}"
export EXPECT_RUN_ID="${BASSEZ_SUBMIT_RUN_ID}"
export EXPECT_RUN_ROOT="${RUNS_ROOT}/${BASSEZ_SUBMIT_RUN_ID}"
BASSEZ_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  STAGE2_SUBMITTER_TEST=1 bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --datasets Bassez --steps bassez_cellsubtype
)"
BASSEZ_RUN_ID="$(printf '%s\n' "${BASSEZ_OUTPUT}" | sed -n 's/^STAGE2_RUN_ID=//p')"
[[ -n "${BASSEZ_RUN_ID}" ]]
BASSEZ_MANIFEST="${RUNS_ROOT}/${BASSEZ_RUN_ID}/manifests/steps.tsv"
IFS=$'\t' read -r BASSEZ_STEP BASSEZ_SCRIPT BASSEZ_OUTPUTS BASSEZ_DEPENDENCY BASSEZ_OWNER \
  < "${BASSEZ_MANIFEST}"
[[ "${BASSEZ_STEP}" == "bassez_cellsubtype" ]]
[[ "${BASSEZ_SCRIPT}" == "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.6_submit_bassez.sh" ]]
case "${BASSEZ_SCRIPT}" in "${ROOT}"/*) echo "Bassez manifest used mutable script" >&2; exit 1 ;; esac
[[ "${BASSEZ_OUTPUTS}" == "${TMP_DIR}/home/scratch/ECODA_paper/Bassez/data/BassezA_2021_33958794whole.rds" ]]
[[ "${BASSEZ_DEPENDENCY}" == "-" && "${BASSEZ_OWNER}" != "-" ]]
BASSEZ_OUTPUT="${TMP_DIR}/home/scratch/ECODA_paper/Bassez/data/BassezA_2021_33958794whole.rds"
mkdir -p "$(dirname "${BASSEZ_OUTPUT}")"
printf 'valid Bassez payload\n' > "${BASSEZ_OUTPUT}"
BASSEZ_DIGEST="$(md5sum "${BASSEZ_OUTPUT}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${BASSEZ_DIGEST}" \
  "$(wc -c < "${BASSEZ_OUTPUT}" | tr -d '[:space:]')" "${BASSEZ_OUTPUT}" \
  > "${BASSEZ_OUTPUT}.md5"
BASSEZ_SKIP_RUN_ID="stage2_submitter_bassez_noop"
export ECODA_RUN_ID="${BASSEZ_SKIP_RUN_ID}"
export EXPECT_RUN_ID="${BASSEZ_SKIP_RUN_ID}"
export EXPECT_RUN_ROOT="${RUNS_ROOT}/${BASSEZ_SKIP_RUN_ID}"
BASSEZ_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
BASSEZ_SKIP_OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --datasets Bassez --steps bassez_cellsubtype
)"
case "${BASSEZ_SKIP_OUTPUT}" in *"NOOP_VALIDATED=${BASSEZ_SKIP_RUN_ID}"*) ;; *) echo "valid row was not validator-only skipped" >&2; exit 1 ;; esac
BASSEZ_SKIP_ROOT="${RUNS_ROOT}/${BASSEZ_SKIP_RUN_ID}"
[[ "$(grep '^STATE=' "${BASSEZ_SKIP_ROOT}/status/terminal")" == "STATE=NOOP_VALIDATED" ]]
[[ -s "${BASSEZ_SKIP_ROOT}/status/noop" ]]
[[ "$(wc -l < "${BASSEZ_SKIP_ROOT}/manifests/scheduler_ids.tsv" | tr -d '[:space:]')" == "0" ]]
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${BASSEZ_CALLS_BEFORE}" ]]
LEGACY_SYNC_RUN_ID="stage2_legacy_sync"
LEGACY_SYNC_ROOT="${RUNS_ROOT}/${LEGACY_SYNC_RUN_ID}"
mkdir -p "${LEGACY_SYNC_ROOT}/manifests" "${LEGACY_SYNC_ROOT}/status"
printf 'STAGE=stage2\nRUN_ID=%s\nSTATE=ACTIVE\n' "${LEGACY_SYNC_RUN_ID}" \
  > "${LEGACY_SYNC_ROOT}/metadata"
LEGACY_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --sync-only "${LEGACY_SYNC_RUN_ID}" > "${TMP_DIR}/legacy.sync.log" 2>&1; then
  echo "legacy unpinned sync unexpectedly succeeded" >&2
  exit 1
fi
grep -q "legacy_source_unpinned" "${TMP_DIR}/legacy.sync.log"
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${LEGACY_CALLS_BEFORE}" ]]

OUTSIDE_RUN_ID="stage2_outside_script"
OUTSIDE_ROOT="${RUNS_ROOT}/${OUTSIDE_RUN_ID}"
mkdir -p "${OUTSIDE_ROOT}/manifests" "${OUTSIDE_ROOT}/status"
printf 'STAGE=stage2\nRUN_ID=%s\nSTATE=ACTIVE\n' "${OUTSIDE_RUN_ID}" \
  > "${OUTSIDE_ROOT}/metadata"
cp "${SOURCE_MANIFEST}" "${OUTSIDE_ROOT}/manifests/source.manifest"
cp "${SUBMIT_ROOT}/manifests/runtime.identity" \
  "${OUTSIDE_ROOT}/manifests/runtime.identity"
OUTSIDE_SCRIPT="${TMP_DIR}/outside-stage2.sh"
printf '#!/bin/bash\nexit 0\n' > "${OUTSIDE_SCRIPT}"
chmod +x "${OUTSIDE_SCRIPT}"
OUTSIDE_MANIFEST="${OUTSIDE_ROOT}/manifests/steps.tsv"
printf 'bassez_cellsubtype\t%s\t%s\t-\t-\n' "${OUTSIDE_SCRIPT}" \
  "${BASSEZ_OUTPUT}" > "${OUTSIDE_MANIFEST}"
OUTSIDE_DIGEST="$(md5sum "${OUTSIDE_MANIFEST}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${OUTSIDE_DIGEST}" \
  "$(wc -c < "${OUTSIDE_MANIFEST}" | tr -d '[:space:]')" "${OUTSIDE_MANIFEST}" \
  > "${OUTSIDE_MANIFEST}.md5"
: > "${OUTSIDE_ROOT}/manifests/jobs.tsv"
OUTSIDE_CALLS_BEFORE="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
if HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --sync-only "${OUTSIDE_RUN_ID}" > "${TMP_DIR}/outside.sync.log" 2>&1; then
  echo "outside-root Stage 2 script unexpectedly accepted" >&2
  exit 1
fi
grep -q "script mismatch" "${TMP_DIR}/outside.sync.log"
[[ "$(wc -l < "${CAPTURE}" | tr -d '[:space:]')" == "${OUTSIDE_CALLS_BEFORE}" ]]


# Guarded CombinedPBMC legacy-raw migration: valid content and sidecar move to
# the canonical raw basename, with PATH rewritten and no duplicate left.
COMBINED_DIR="${TMP_DIR}/home/scratch/ECODA_paper/CombinedPBMC/data"
mkdir -p "${COMBINED_DIR}"
OLD="${COMBINED_DIR}/combined_pbmc_batch_effect_analysis.h5ad"
NEW="${COMBINED_DIR}/combined_pbmc.h5ad"
pixi run python -c 'import anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,sys; a=ad.AnnData(X=sp.csr_matrix([[1,0],[0,2]],dtype="float32"),obs=pd.DataFrame({"Sample":["s1","s2"],"cond":["Healthy","Healthy"],"batch":["A","B"]},index=["c1","c2"]),var=pd.DataFrame(index=["g1","g2"])); a.write_h5ad(sys.argv[1])' "${OLD}"
digest="$(md5sum "${OLD}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${OLD}" | tr -d '[:space:]')" "${OLD}" > "${OLD}.md5"
rm -rf "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners"
: > "${CAPTURE}"
MIGRATION_RUN_ID="stage2_submitter_migration"
export ECODA_RUN_ID="${MIGRATION_RUN_ID}"
export EXPECT_RUN_ID="${MIGRATION_RUN_ID}"
export EXPECT_RUN_ROOT="${RUNS_ROOT}/${MIGRATION_RUN_ID}"
export STAGE2_RUN_ROOT="${RUNS_ROOT}/${MIGRATION_RUN_ID}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  STAGE2_SUBMITTER_TEST=1 bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
  --datasets CombinedPBMC --steps combinedpbmc >/dev/null
[[ -s "${NEW}" && -s "${NEW}.md5" ]]
[[ ! -e "${OLD}" && ! -e "${OLD}.md5" ]]
grep -q "^PATH=${NEW}$" "${NEW}.md5"
MIGRATION_ROOT="${RUNS_ROOT}/${MIGRATION_RUN_ID}"
export ECODA_LOGS_DIR="${MIGRATION_ROOT}/logs"
export ECODA_RUNTIME_IMAGE_SHA256="${IMAGE_SHA}"
export ECODA_RUNTIME_MANIFEST_SHA256="${RUNTIME_MANIFEST_SHA}"
export ECODA_RUNTIME_IMAGE_SIZE="${IMAGE_SIZE}"
export ECODA_RUNTIME_MANIFEST_SIZE="${RUNTIME_MANIFEST_SIZE}"

# Hook-level force propagation with a temporary slurm_config/python stub.
HOOK_ROOT="${TMP_DIR}/hook-project"
mkdir -p "${HOOK_ROOT}/src/2_dataset_specific_preprocessing" "${HOOK_ROOT}/src/utils/bash" "${HOOK_ROOT}/bin"
cp "${ROOT}/src/2_dataset_specific_preprocessing/1.5_submit_myocardial.sh" "${HOOK_ROOT}/src/2_dataset_specific_preprocessing/"
cp "${ROOT}/src/utils/bash/ecoda_runtime.sh" "${HOOK_ROOT}/src/utils/bash/"
cp "${ROOT}/src/2_dataset_specific_preprocessing/1.2_submit_combinedpbmc.sh" "${HOOK_ROOT}/src/2_dataset_specific_preprocessing/"
printf 'PROJECT_ROOT="%s"\nPYTHON_BIN="%s/bin/python"\n' "${HOOK_ROOT}" "${HOOK_ROOT}" > "${HOOK_ROOT}/src/slurm_config.sh"
printf '#!/bin/bash\nprintf "%%s\\n" "$*" > "%s/force.args"\n' "${HOOK_ROOT}" > "${HOOK_ROOT}/bin/python"
printf '#!/bin/bash\n: > "%s/module.called"\n' "${HOOK_ROOT}" > "${HOOK_ROOT}/bin/module"
chmod +x "${HOOK_ROOT}/bin/python" "${HOOK_ROOT}/bin/module"
HOME="${TMP_DIR}/home" FORCE_PREPROCESS=1 HOOK_FORCE_FILE="${HOOK_ROOT}/force.args" \
  ECODA_RUNTIME_MODE=host ECODA_RUNTIME_IN_CONTAINER=0 \
  bash "${HOOK_ROOT}/src/2_dataset_specific_preprocessing/1.5_submit_myocardial.sh"
case "$(cat "${HOOK_ROOT}/force.args")" in *"--force"*) ;; *) echo "myocardial hook dropped --force" >&2; exit 1 ;; esac
ECODA_RUNTIME_IN_CONTAINER=1 ECODA_RUNTIME_MODE=apptainer \
  ECODA_RUNTIME_PREFIX="${HOST_ENV}" SLURM_JOB_ID=999999 \
  HOME="${TMP_DIR}/home" HOOK_FORCE_FILE="${HOOK_ROOT}/force.args" \
  PATH="${HOOK_ROOT}/bin:${TMP_DIR}/bin:${PATH}" \
  bash "${HOOK_ROOT}/src/2_dataset_specific_preprocessing/1.2_submit_combinedpbmc.sh"
[[ ! -e "${HOOK_ROOT}/module.called" ]] || { echo "CombinedPBMC loaded host module inside container" >&2; exit 1; }
echo "stage2 submitter: OK"
