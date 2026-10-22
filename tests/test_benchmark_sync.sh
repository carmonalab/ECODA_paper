#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
unset ANALYSIS_VARIANT ANALYSIS_PASS ANALYSIS_LOG_PREFIX
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-sync.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
cleanup() {
  chmod -R u+w "${TMP_DIR}" >/dev/null 2>&1 || true
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT
mode_of() {
  if stat -c '%a' "$1" >/dev/null 2>&1; then
    stat -c '%a' "$1"
  else
    stat -f '%Lp' "$1"
  fi
}
export REAL_RSYNC="$(command -v rsync)"
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home" "${TMP_DIR}/scratch/benchmark/embeddings" "${TMP_DIR}/nas/project"
cat > "${TMP_DIR}/bin/rsync" <<'STUB'
#!/bin/bash
set -euo pipefail
[[ -z "${RSYNC_CAPTURE:-}" ]] || printf '%s\n' "$*" >> "${RSYNC_CAPTURE}"
[[ "${FAIL_RSYNC:-0}" != 1 ]] || exit 23
if [[ "${USE_REAL_RSYNC:-0}" == 1 ]]; then
  exec "${REAL_RSYNC}" "$@"
fi
files_from=""
chmod_spec=""
for arg in "$@"; do
  case "${arg}" in
    --files-from=*) files_from="${arg#--files-from=}" ;;
    --chmod=*) chmod_spec="${arg#--chmod=}" ;;
  esac
done
n="$#"
eval "src=\${$((n - 1))}"
eval "dst=\${$n}"
src="${src%/}"
dst="${dst%/}"
mkdir -p "${dst}"
if [[ -n "${files_from}" ]]; then
  while IFS= read -r rel || [[ -n "${rel}" ]]; do
    [[ -n "${rel}" ]] || continue
    mkdir -p "${dst}/$(dirname "${rel}")"
    cp -f "${src}/${rel}" "${dst}/${rel}"
    if [[ "${chmod_spec}" == "ugo+rwX" ]]; then
      chmod 666 "${dst}/${rel}"
      dir="${dst}/$(dirname "${rel}")"
      while [[ "${dir}" != "${dst}" ]]; do
        chmod 777 "${dir}"
        next_dir="$(dirname "${dir}")"
        [[ "${next_dir}" != "${dir}" ]] || break
        dir="${next_dir}"
      done
      chmod 777 "${dst}"
    fi
    if [[ "${FAIL_RSYNC_AFTER_METADATA:-0}" == 1 && "${rel}" == metadata/* && "${rel}" != *.md5 ]]; then
      exit 23
    fi
    if [[ "${FAIL_RSYNC_AFTER_LOG:-0}" == 1 && "${rel}" == embeddings/execution_times.feather ]]; then
      exit 23
    fi
  done < "${files_from}"
else
  cp -R "${src}". "${dst}"
fi
STUB
chmod +x "${TMP_DIR}/bin/rsync"
SCHEDULER_CAPTURE="${TMP_DIR}/scheduler.calls"
export SCHEDULER_CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${SCHEDULER_CAPTURE}"
exit 0
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
export PATH="${TMP_DIR}/bin:${PATH}"
export RSYNC_CAPTURE="${TMP_DIR}/rsync.calls"
export HPC_SCRATCH_DIR="${TMP_DIR}/scratch"
export LOGS_DIR="${TMP_DIR}/logs"
export NAS_TARGET_DIR="${TMP_DIR}/nas/project"
export ANALYSIS_ROOT="${TMP_DIR}/scratch/benchmark"
export ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/benchmark"
export USER_EMAIL="test@example.invalid"
source "${ROOT}/src/slurm_config.sh" >/dev/null 2>&1 || true
export HPC_SCRATCH_DIR="${TMP_DIR}/scratch" NAS_TARGET_DIR="${TMP_DIR}/nas/project" \
  ANALYSIS_ROOT="${TMP_DIR}/scratch/benchmark" ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/benchmark"
mkdir -p "${NAS_TARGET_DIR}" "${ANALYSIS_NAS_ROOT}"
mkdir -p "${TMP_DIR}/tmp"
export TMPDIR="${TMP_DIR}/tmp"
source "${ROOT}/src/utils/bash/ecoda_run_common.sh"
source "${ROOT}/src/5_run_benchmark_methods/benchmark_submit_common.sh"
source "${ROOT}/src/5_run_benchmark_methods/stage5_pending.sh"
source "${ROOT}/src/utils/bash/ecoda_runtime.sh"
unset ECODA_PRODUCER_RUN_ID ECODA_ARTIFACT_PRODUCER
export ECODA_EXECUTION_LOG_PRODUCER=stage5_execution_log
export PYTHONDONTWRITEBYTECODE=1
# Simulate only the finite NAS boundary; publication/owner behavior stays real.
NAS_UTILITY_STUB="${TMP_DIR}/bin/nas-fixture"
cat > "${NAS_UTILITY_STUB}" <<'STUB'
#!/bin/bash
set -euo pipefail
if [[ "$1" == preflight ]]; then
  printf '{"state":"OK"}\n' > "$5"
  exit 0
fi
[[ "$1" == verify && "$2" == --receipt ]]
pixi run python - "$3" <<'PY'
import hashlib, json, pathlib, sys
r = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert r["state"] == "OK"
for key in ("selection", "inventory", "source_manifest"):
    assert hashlib.sha256(pathlib.Path(r[key + "_path"]).read_bytes()).hexdigest() == r[key + "_sha256"]
for line in pathlib.Path(r["inventory_path"]).read_text().splitlines():
    rel, md5, size = line.split("\t")
    for root in ("destination_root",):
        data = (pathlib.Path(r[root]) / rel).read_bytes()
        assert len(data) == int(size) and hashlib.md5(data).hexdigest() == md5
PY
STUB
chmod +x "${NAS_UTILITY_STUB}"
benchmark_source_script_path() {
  if [[ "$1" == src/utils/bash/ecoda_nas_backup.sh ]]; then
    printf '%s\n' "${NAS_UTILITY_STUB}"
  else
    printf '%s/%s\n' "${ECODA_SOURCE_ROOT:-${ROOT}}" "$1"
  fi
}

seal_publication_receipt() {
  local promotion="${ECODA_RUN_ROOT}/promotion" owner
  owner="$(ecoda_owner_dir stage5 "sync/${ANALYSIS_ROOT}")/owner"
  mkdir -p "${promotion}/manifests"
  cp "${ECODA_SELECTION_MANIFEST}" "${promotion}/manifests/selection.tsv"
  export ECODA_EXPECTED_PENDING="${ECODA_RUN_ROOT}/manifests/expected_pending.tsv"
  : > "${ECODA_EXPECTED_PENDING}"
  stage5_write_pending_checksum "${ECODA_EXPECTED_PENDING}"
  cp "${ECODA_EXPECTED_PENDING}" "${promotion}/manifests/expected_pending.tsv"
  pixi run python - "${promotion}" "${TMP_DIR}" "${NAS_TARGET_DIR}" \
    "${ANALYSIS_NAS_ROOT}" "${owner}" "${SOURCE_MANIFEST}" "${ECODA_RUN_ID}" <<'PY'
import csv, hashlib, json, pathlib, shutil, sys
p, base, nas, remote, owner, source = map(pathlib.Path, sys.argv[1:7])
run_id = sys.argv[7]
def digest(path, algorithm="sha256"):
    return hashlib.new(algorithm, path.read_bytes()).hexdigest()
def receipt(kind, root, destination, files):
    inv = p / "manifests" / (kind + ".inventory.tsv")
    inv.write_text("".join(f"{f.relative_to(root)}\t{digest(f, 'md5')}\t{f.stat().st_size}\n" for f in files))
    for f in files:
        target = destination / f.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f, target)
    selection = p / "manifests/selection.tsv"
    r = dict(schema_version=1, state="OK", run_id=run_id, host="fixture",
             nas_root=str(nas), mount_source="fixture:/NAS", source_root=str(root),
             destination_root=str(destination), selection_path=str(selection),
             selection_sha256=digest(selection), inventory_path=str(inv),
             inventory_sha256=digest(inv), source_manifest_path=str(source),
             source_manifest_sha256=digest(source), completed_at="2026-10-21T00:00:00Z")
    target = p / "manifests" / (kind + ".receipt.json")
    target.write_text(json.dumps(r) + "\n")
    return target
archive_destination = nas / "_ecoda_backups" / run_id / "before"
archived = [f for f in remote.rglob("*") if f.is_file()]
if owner.exists():
    archived.append(owner)
receipt("archive", base, archive_destination, archived)
mapping = p / "manifests/archive_mapping.tsv"
with mapping.open("w") as stream:
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(["original_path", "archive_path", "md5", "size"])
    for f in archived:
        writer.writerow([f, archive_destination / f.relative_to(base),
                         digest(f, "md5"), f.stat().st_size])
receipt("inputs", p, nas / "_ecoda_backups" / run_id / "inputs", [])
receipt("canonical", p, nas / "_ecoda_backups" / run_id / "canonical", [])
auth = p / "manifests/global_sync_owner.tsv"
with auth.open("w") as stream:
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(["path", "md5", "size", "run_id", "state", "archive_path"])
    if owner.exists():
        fields = dict(line.split("=", 1) for line in owner.read_text().splitlines())
        writer.writerow([owner, digest(owner, "md5"), owner.stat().st_size,
                         fields["RUN_ID"], fields["STATE"], archive_destination / owner.relative_to(base)])
files = [auth, mapping, p / "manifests/expected_pending.tsv"]
for kind in ("inputs", "canonical", "archive"):
    files.extend([p / "manifests" / (kind + ".receipt.json"),
                  p / "manifests" / (kind + ".inventory.tsv")])
receipt("precompute_backup", p, nas / "_ecoda_backups" / run_id / "composite", files)
PY
  export ECODA_PRECOMPUTE_BACKUP_RECEIPT="${promotion}/manifests/precompute_backup.receipt.json"
  export ECODA_PRECOMPUTE_BACKUP_RECEIPT_SHA256="$(sha256_file "${ECODA_PRECOMPUTE_BACKUP_RECEIPT}")"
  export ECODA_EXPECTED_PENDING_SHA256="$(sha256_file "${ECODA_EXPECTED_PENDING}")"
}

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

SNAPSHOT_ID="$(printf 'a%.0s' {1..40})"
SNAPSHOT_ROOT="${TMP_DIR}/source_snapshots/${SNAPSHOT_ID}"
SOURCE_TREE="${SNAPSHOT_ROOT}/tree"
SOURCE_IDENTITY_DIR="${SNAPSHOT_ROOT}/identity"
SOURCE_MANIFEST="${SOURCE_IDENTITY_DIR}/source.manifest"
SOURCE_ARCHIVE="${SOURCE_IDENTITY_DIR}/source.tar"
mkdir -p "${SOURCE_TREE}/src" "${SOURCE_TREE}/aux" "${SOURCE_IDENTITY_DIR}"
printf '#!/bin/bash\nexit 0\n' > "${SOURCE_TREE}/src/snapshot-worker.sh"
chmod +x "${SOURCE_TREE}/src/snapshot-worker.sh"
printf 'config\n' > "${SOURCE_TREE}/config_helper.R"
printf '{}\n' > "${SOURCE_TREE}/datasets.json"
printf '[project]\n' > "${SOURCE_TREE}/pixi.toml"
printf 'lock\n' > "${SOURCE_TREE}/pixi.lock"
printf 'scgate\n' > "${SOURCE_TREE}/aux/scGateDB.rds"
printf 'blocklist\n' > "${SOURCE_TREE}/aux/genes.blocklist.rds"
printf 'ensembl\n' > "${SOURCE_TREE}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
tar -cf "${SOURCE_ARCHIVE}" -C "${SOURCE_TREE}" .
SOURCE_ARCHIVE_SHA256="$(sha256_file "${SOURCE_ARCHIVE}")"
SOURCE_CONFIG_SHA256="$(sha256_file "${SOURCE_TREE}/config_helper.R")"
SOURCE_DATASETS_SHA256="$(sha256_file "${SOURCE_TREE}/datasets.json")"
SOURCE_TOML_SHA256="$(sha256_file "${SOURCE_TREE}/pixi.toml")"
SOURCE_LOCK_SHA256="$(sha256_file "${SOURCE_TREE}/pixi.lock")"
SOURCE_SC_GATE_BRANCH="${SCGATE_DB_BRANCH:-main}"
cat > "${SOURCE_MANIFEST}" <<EOF
FORMAT=1
SOURCE_ROOT=${SOURCE_TREE}
SOURCE_COMMIT=${SNAPSHOT_ID}
SOURCE_ARCHIVE_PATH=${SOURCE_ARCHIVE}
SOURCE_ARCHIVE_SHA256=${SOURCE_ARCHIVE_SHA256}
CONFIG_HELPER_SHA256=${SOURCE_CONFIG_SHA256}
DATASETS_SHA256=${SOURCE_DATASETS_SHA256}
PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}
PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}
AUX_ROOT=${SOURCE_TREE}/aux
SCGATE_DB_BRANCH=${SOURCE_SC_GATE_BRANCH}
EOF
printf 'complete\n' > "${SNAPSHOT_ROOT}/COMPLETE"
chmod -R a-w "${SOURCE_TREE}" "${SOURCE_IDENTITY_DIR}" "${SNAPSHOT_ROOT}/COMPLETE"

RUNTIME_DIR="${TMP_DIR}/runtime/_ecoda_runtime/runtime-a"
RUNTIME_IMAGE="${RUNTIME_DIR}/ecoda-py-cuda13.sif"
RUNTIME_MANIFEST="${RUNTIME_IMAGE}.manifest"
mkdir -p "${RUNTIME_DIR}"
printf 'runtime-image\n' > "${RUNTIME_IMAGE}"
RUNTIME_IMAGE_SHA256="$(sha256_file "${RUNTIME_IMAGE}")"
cat > "${RUNTIME_MANIFEST}" <<EOF
FORMAT=2
IMAGE_PATH=${RUNTIME_IMAGE}
IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}
RUNTIME_ENV=py-cuda13
RUNTIME_LAYOUT=relocated
CONTAINER_ENV_PREFIX=/opt/ecoda/py-cuda13
BASE_IMAGE=rockylinux:9
PIXITAINER_VERSION=0.8.3
PIXI_VERSION=0.49.0
APPTAINER_VERSION=1.3.2
IMAGE_BUILD_GIT_REVISION=build-a
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}
EOF
RUNTIME_MANIFEST_SHA256="$(sha256_file "${RUNTIME_MANIFEST}")"
RUNTIME_IMAGE_SIZE="$(wc -c < "${RUNTIME_IMAGE}" | tr -d '[:space:]')"
RUNTIME_MANIFEST_SIZE="$(wc -c < "${RUNTIME_MANIFEST}" | tr -d '[:space:]')"
chmod 444 "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}"
chmod 555 "${RUNTIME_DIR}"

prepare_run_source_write() {
  local run_root="$1"
  local source_copy="${run_root}/manifests/source.manifest"
  chmod u+w "${run_root}/manifests"
  if [[ -e "${source_copy}" ]]; then
    chmod u+w "${source_copy}"
  fi
}

write_run_identity() {
  local run_root="$1"
  local source_copy="${run_root}/manifests/source.manifest"
  local runtime_identity="${run_root}/manifests/runtime.identity"
  prepare_run_source_write "${run_root}"
  if [[ -e "${runtime_identity}" ]]; then
    chmod u+w "${runtime_identity}"
  fi
  cp "${SOURCE_MANIFEST}" "${source_copy}"
  cat > "${runtime_identity}" <<EOF
RUNTIME_IMAGE=${RUNTIME_IMAGE}
RUNTIME_MANIFEST=${RUNTIME_MANIFEST}
RUNTIME_IMAGE_SHA256=${RUNTIME_IMAGE_SHA256}
RUNTIME_MANIFEST_SHA256=${RUNTIME_MANIFEST_SHA256}
RUNTIME_IMAGE_SIZE=${RUNTIME_IMAGE_SIZE}
RUNTIME_MANIFEST_SIZE=${RUNTIME_MANIFEST_SIZE}
IMAGE_PIXI_TOML_SHA256=${SOURCE_TOML_SHA256}
IMAGE_PIXI_LOCK_SHA256=${SOURCE_LOCK_SHA256}
EOF
  printf 'SOURCE_MANIFEST=%s\nRUNTIME_IDENTITY=%s\nSOURCE_ROOT=%s\nRUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\n' \
    "${source_copy}" \
    "${runtime_identity}" \
    "${SOURCE_TREE}" "${RUNTIME_IMAGE}" "${RUNTIME_MANIFEST}" >> "${run_root}/metadata"
  chmod a-w "${source_copy}" "${runtime_identity}"
}

write_execution_log() {
  local seconds="$1"
  local log="${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
  mkdir -p "${ECODA_RUN_ROOT}/logs"
  pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[float(sys.argv[2])],"mem_GB":[1.0]}).to_feather(sys.argv[1])' \
    "${log}" "${seconds}"
  ecoda_write_checksum "${log}"
  ecoda_write_artifact_record "${log}" stage5_execution_log "${ECODA_RUN_ID}" >/dev/null
}

prepare_sync_run() {
  local run_id="$1" seconds="$2"
  ecoda_init_run stage5 "${run_id}" >/dev/null
  export ECODA_RUN_ID ECODA_RUN_ROOT
  write_run_identity "${ECODA_RUN_ROOT}"
  export EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs"
  write_execution_log "${seconds}"
  export ECODA_SELECTION_MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
  printf 'Adams\tbenchmark_analysis\tbenchmark_analysis\n' > "${ECODA_SELECTION_MANIFEST}"
  ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}"
}

expect_sync_failure() {
  local label="$1"
  local owner before
  owner="$(ecoda_owner_dir stage5 "sync/${ANALYSIS_ROOT}")/owner"
  before="$(sha256_file "${owner}")"
  if analysis_merge_sync_cleanup "${LABELS[@]}" >/dev/null 2>&1; then
    echo "expected benchmark sync failure: ${label}" >&2
    exit 1
  fi
  [[ -s "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather" ]] ||
    { echo "failed sync cleaned run logs: ${label}" >&2; exit 1; }
  [[ "$(sha256_file "${owner}")" == "${before}" ]] ||
    { echo "validation failure changed terminal owner: ${label}" >&2; exit 1; }
}

ECODA_RUN_ID="sync_run"
ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
export ECODA_RUN_ID ECODA_RUN_ROOT
printf 'Adams\tbenchmark_analysis\tbenchmark_analysis\n' > "${ECODA_RUN_ROOT}/manifests/selection.tsv"
export ECODA_SELECTION_MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv" ECODA_EXACT_SELECTION=0
ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}"
write_run_identity "${ECODA_RUN_ROOT}"
export ANALYSIS_ROOT ANALYSIS_NAS_ROOT EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs" \
  ANALYSIS_LOG_PREFIX="execution_times_"
DATASET_NAMES=(Adams)
LABELS=(mrvi)
export DATASET_NAMES LABELS
for n in 1000 2000 3000; do
  path="${ANALYSIS_ROOT}/embeddings/Adams_hvg${n}_mrvi_dists.feather"
  pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"s1":[1.0,0.0],"s2":[0.0,1.0]},index=["s1","s2"]).to_feather(sys.argv[1])' "${path}"
  ecoda_write_checksum "${path}"
  ecoda_write_artifact_record "${path}" mrvi "${ECODA_RUN_ID}" >/dev/null
  pixi run python -c 'import hashlib,json,pathlib,sys; artifact=pathlib.Path(sys.argv[1]); metadata=pathlib.Path(str(artifact)+".runtime.json"); metadata.write_text(json.dumps({"schema_version":1,"artifact_path":str(artifact),"artifact_md5":hashlib.md5(artifact.read_bytes()).hexdigest(),"dataset":"Adams","method":f"MrVI_hvg{sys.argv[2]}","time_secs":1.0,"mem_GB":1.0},separators=(",",":"))+"\n")' "${path}" "${n}"
  ecoda_write_checksum "${path}.runtime.json"
  ecoda_write_artifact_record "${path}.runtime.json" mrvi "${ECODA_RUN_ID}" >/dev/null
done
for n in 1000 2000 3000; do
  path="${ANALYSIS_ROOT}/embeddings/Adams_hvg${n}_mrvi_dists.feather"
  [[ -s "${path}.runtime.json" ]]
  [[ -s "${path}.runtime.json.md5" ]]
done
mkdir -p "${ECODA_RUN_ROOT}/logs"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[1.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
ecoda_write_checksum "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather"
ecoda_write_artifact_record "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather" \
  stage5_execution_log "${ECODA_RUN_ID}" >/dev/null
for combo in hvg1000 hvg3000; do
  shard_log="${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams_${combo}.feather"
  pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":[sys.argv[2]],"time_secs":[1.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' \
    "${shard_log}" "MrVI_${combo}"
  ecoda_write_checksum "${shard_log}"
  ecoda_write_artifact_record "${shard_log}" stage5_execution_log "${ECODA_RUN_ID}" >/dev/null
done
export DATASETS_JSON_FILE="${ROOT}/datasets.json"
OWNER_SELECTION="${TMP_DIR}/owner-selection.tsv"
printf 'Adams\tbenchmark_analysis\tmrvi\n' > "${OWNER_SELECTION}"
ecoda_stage5_validate_output_ownership "${OWNER_SELECTION}" "${ECODA_RUN_ID}" >/dev/null
ecoda_stage5_validate_output_ownership "${OWNER_SELECTION}" "${ECODA_RUN_ID}" >/dev/null
if ecoda_stage5_validate_output_ownership "${OWNER_SELECTION}" other_sync_run \
    >/dev/null 2>&1; then
  sbatch --wrap=unexpected-owner-retry
  echo "other-run active owner was accepted before retry boundary" >&2
  exit 1
fi
[[ ! -e "${SCHEDULER_CAPTURE}" ]]
ecoda_owner_finalize_tracked OK "benchmark sync test publication" >/dev/null
printf 'unrelated remote\n' > "${ANALYSIS_NAS_ROOT}/unrelated.txt"
UNRELATED_MODE="$(mode_of "${ANALYSIS_NAS_ROOT}/unrelated.txt")"
analysis_merge_sync_cleanup "${LABELS[@]}"
[[ "$(pixi run python -c 'import pandas as pd,sys; print(len(pd.read_feather(sys.argv[1])))' "${ANALYSIS_NAS_ROOT}/embeddings/execution_times.feather")" == 3 ]]
[[ -s "${ANALYSIS_NAS_ROOT}/checksums.md5" ]]
[[ -s "${ANALYSIS_NAS_ROOT}/embeddings/execution_times.feather" ]]
[[ -s "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather" ]]
[[ -s "${ANALYSIS_NAS_ROOT}/unrelated.txt" ]]
[[ "$(mode_of "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather")" == "666" ]]
[[ "$(mode_of "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather.md5")" == "666" ]]
[[ "$(mode_of "${ANALYSIS_NAS_ROOT}/embeddings")" == "777" ]]
[[ "$(mode_of "${ANALYSIS_NAS_ROOT}/checksums.md5")" == "666" ]]
[[ "$(mode_of "${ANALYSIS_NAS_ROOT}/unrelated.txt")" == "${UNRELATED_MODE}" ]]
[[ ! -e "${ECODA_RUN_ROOT}/logs/execution_times_mrvi_Adams.feather" ]]
[[ ! -d "${ANALYSIS_ROOT}/.sync.lock" ]]
[[ "$(ecoda_owner_state "$(ecoda_owner_dir stage5 "sync/${ANALYSIS_ROOT}")")" == "OK" ]]
[[ "$(grep -c '^checksums.md5$' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv")" == 1 ]]
[[ "$(grep -c 'unrelated' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv" || true)" == 0 ]]
[[ "$(grep -c '^embeddings/Adams_hvg1000_mrvi_dists.feather$' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv")" == 1 ]]
[[ "$(grep -c '^embeddings/Adams_hvg1000_mrvi_dists.feather.runtime.json$' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv")" == 1 ]]
[[ "$(grep -c '^embeddings/Adams_hvg1000_mrvi_dists.feather.runtime.json.md5$' "${ECODA_RUN_ROOT}/manifests/sync_files.tsv")" == 1 ]]
[[ "$(grep -c 'embeddings/Adams_hvg1000_mrvi_dists.feather.runtime.json' "${ANALYSIS_NAS_ROOT}/checksums.md5")" == 1 ]]
[[ -s "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather.runtime.json" ]]
[[ -s "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather.runtime.json.md5" ]]
cmp "${SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
[[ "$(grep -c '^RUNTIME_IMAGE=' "${ECODA_RUN_ROOT}/manifests/runtime.identity")" == 1 ]]
[[ "$(grep -c "^RUNTIME_IMAGE=${RUNTIME_IMAGE}$" "${ECODA_RUN_ROOT}/manifests/runtime.identity")" == 1 ]]
[[ "$(grep -c "^RUNTIME_MANIFEST=${RUNTIME_MANIFEST}$" "${ECODA_RUN_ROOT}/manifests/runtime.identity")" == 1 ]]
grep -q "^SOURCE_MANIFEST=${ECODA_RUN_ROOT}/manifests/source.manifest$" \
  "${ECODA_RUN_ROOT}/metadata"
grep -q "^RUNTIME_IDENTITY=${ECODA_RUN_ROOT}/manifests/runtime.identity$" \
  "${ECODA_RUN_ROOT}/metadata"
for n in 1000 2000 3000; do
  path="${ANALYSIS_ROOT}/embeddings/Adams_hvg${n}_mrvi_dists.feather"
  ecoda_validate_artifact_record "${path}" mrvi sync_run >/dev/null
  ecoda_artifact_owner_validate "${path}" sync_run >/dev/null
  [[ "${ECODA_ARTIFACT_OWNER_STATE}" == "OK" ]]
done
MERGED_LOG="${ANALYSIS_ROOT}/embeddings/execution_times.feather"
ecoda_validate_artifact_record "${MERGED_LOG}" stage5_execution_log sync_run >/dev/null
(
  cd "${ANALYSIS_NAS_ROOT}"
  md5sum unrelated.txt >> checksums.md5
)
grep -q 'unrelated.txt' "${ANALYSIS_NAS_ROOT}/checksums.md5"
PAYLOAD="${ANALYSIS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather"
PAYLOAD_RECORD="$(ecoda_artifact_record_path "${PAYLOAD}" sync_run)"
cp "${PAYLOAD_RECORD}" "${PAYLOAD_RECORD}.saved"
printf 'PATH=%s\n' "${PAYLOAD}" > "${PAYLOAD_RECORD}"
if ecoda_validate_artifact_record "${PAYLOAD}" mrvi sync_run >/dev/null 2>&1; then
  echo "malformed artifact record was accepted" >&2
  exit 1
fi
rm -f "${PAYLOAD_RECORD}"
if ecoda_validate_artifact_record "${PAYLOAD}" mrvi sync_run >/dev/null 2>&1; then
  echo "missing artifact record was accepted" >&2
  exit 1
fi
mv -f "${PAYLOAD_RECORD}.saved" "${PAYLOAD_RECORD}"
ecoda_validate_artifact_record "${PAYLOAD}" mrvi sync_run >/dev/null

export ECODA_RUN_ROOT="${ECODA_RUNS_ROOT}/sync_run" ECODA_RUN_ID=sync_run \
  ECODA_SOURCE_ROOT="${SOURCE_TREE}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}" \
  ECODA_SOURCE_SNAPSHOT_REQUIRED=1 ECODA_AUX_ROOT="${SOURCE_TREE}/aux" \
  ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}" ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}" \
  ECODA_RUNTIME_MODE=apptainer ECODA_RUNTIME_IN_CONTAINER=1 ECODA_RUNTIME_PROFILE=stage5 \
  ECODA_SCRATCH_ROOT="${HPC_SCRATCH_DIR}" ECODA_LOGS_DIR="${ECODA_RUN_ROOT}/logs"
expect_runtime_failure() {
  local label="$1"
  if ecoda_runtime_validate_bound_run >/dev/null 2>&1; then
    echo "expected bound-runtime failure: ${label}" >&2
    exit 1
  fi
}
ecoda_runtime_validate_bound_run >/dev/null

mv "${ECODA_RUN_ROOT}/manifests/runtime.identity" \
  "${ECODA_RUN_ROOT}/manifests/runtime.identity.saved"
expect_runtime_failure missing-runtime-identity
mv "${ECODA_RUN_ROOT}/manifests/runtime.identity.saved" \
  "${ECODA_RUN_ROOT}/manifests/runtime.identity"
chmod u+w "${ECODA_RUN_ROOT}/manifests" \
  "${ECODA_RUN_ROOT}/manifests/runtime.identity"
printf 'RUNTIME_IMAGE=%s\n' "${RUNTIME_IMAGE}" > \
  "${ECODA_RUN_ROOT}/manifests/runtime.identity"
chmod a-w "${ECODA_RUN_ROOT}/manifests/runtime.identity"
expect_runtime_failure malformed-runtime-identity
write_run_identity "${ECODA_RUN_ROOT}"

mv "${ECODA_RUN_ROOT}/manifests/source.manifest" \
  "${ECODA_RUN_ROOT}/manifests/source.manifest.saved"
expect_runtime_failure missing-source-identity
mv "${ECODA_RUN_ROOT}/manifests/source.manifest.saved" \
  "${ECODA_RUN_ROOT}/manifests/source.manifest"
prepare_run_source_write "${ECODA_RUN_ROOT}"
printf 'FORMAT=1\n' > "${ECODA_RUN_ROOT}/manifests/source.manifest"
chmod a-w "${ECODA_RUN_ROOT}/manifests/source.manifest"
expect_runtime_failure malformed-source-identity
prepare_run_source_write "${ECODA_RUN_ROOT}"
cp "${SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
chmod a-w "${ECODA_RUN_ROOT}/manifests/source.manifest"

cp "${SOURCE_MANIFEST}" "${TMP_DIR}/source.manifest.saved"
chmod u+w "${SOURCE_IDENTITY_DIR}" "${SOURCE_MANIFEST}"
printf 'FORMAT=1\n' > "${SOURCE_MANIFEST}"
chmod a-w "${SOURCE_MANIFEST}" "${SOURCE_IDENTITY_DIR}"
expect_runtime_failure malformed-source-snapshot-record
chmod u+w "${SOURCE_IDENTITY_DIR}" "${SOURCE_MANIFEST}"
mv -f "${TMP_DIR}/source.manifest.saved" "${SOURCE_MANIFEST}"
chmod a-w "${SOURCE_MANIFEST}" "${SOURCE_IDENTITY_DIR}"

cp "${RUNTIME_MANIFEST}" "${TMP_DIR}/runtime.manifest.saved"
chmod u+w "${RUNTIME_DIR}" "${RUNTIME_MANIFEST}"
printf 'FORMAT=2\n' > "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_MANIFEST}" "${RUNTIME_DIR}"
expect_runtime_failure malformed-runtime-manifest-record
chmod u+w "${RUNTIME_DIR}" "${RUNTIME_MANIFEST}"
mv -f "${TMP_DIR}/runtime.manifest.saved" "${RUNTIME_MANIFEST}"
chmod a-w "${RUNTIME_MANIFEST}" "${RUNTIME_DIR}"

INSIDE_SCRIPT="${SOURCE_TREE}/src/snapshot-worker.sh"
OUTSIDE_SCRIPT="${TMP_DIR}/outside-worker.sh"
printf '#!/bin/bash\nexit 0\n' > "${OUTSIDE_SCRIPT}"
chmod +x "${OUTSIDE_SCRIPT}"
[[ "$(ecoda_require_source_script_path "${INSIDE_SCRIPT}" "${SOURCE_TREE}")" == "${INSIDE_SCRIPT}" ]]
if ecoda_require_source_script_path "${OUTSIDE_SCRIPT}" "${SOURCE_TREE}" \
    >/dev/null 2>&1; then
  sbatch --wrap=unexpected-outside-source-retry
  echo "outside source script reached retry boundary" >&2
  exit 1
fi
[[ ! -e "${SCHEDULER_CAPTURE}" ]]
unset ECODA_SOURCE_ROOT ECODA_SOURCE_MANIFEST ECODA_SOURCE_SNAPSHOT_REQUIRED \
  ECODA_AUX_ROOT ECODA_RUNTIME_IMAGE ECODA_RUNTIME_MANIFEST ECODA_RUNTIME_MODE \
  ECODA_RUNTIME_IN_CONTAINER ECODA_RUNTIME_PROFILE ECODA_SCRATCH_ROOT ECODA_LOGS_DIR

prepare_sync_run sync_remote_corrupt 2.0
chmod u+w "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather"
printf 'remote mutation\n' >> \
  "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather"
expect_sync_failure remote-md5-boundary
chmod u+w "${ANALYSIS_NAS_ROOT}/embeddings"
cp -f "${ANALYSIS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather" \
  "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather"
chmod a-w "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather"
cp -f "${ANALYSIS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather.md5" \
  "${ANALYSIS_NAS_ROOT}/embeddings/Adams_hvg1000_mrvi_dists.feather.md5"
chmod a-w "${ANALYSIS_NAS_ROOT}/embeddings"

prepare_sync_run sync_remote_valid 2.0
ecoda_owner_set_state "$(ecoda_owner_dir stage5 "sync/${ANALYSIS_ROOT}")" FAIL \
  "ordinary incomplete publication fixture"
chmod u+w "${ANALYSIS_NAS_ROOT}/embeddings"
analysis_merge_sync_cleanup "${LABELS[@]}"
chmod a-w "${ANALYSIS_NAS_ROOT}/embeddings"

prepare_sync_run sync_missing_checksum 3.0
rm -f "${PAYLOAD}.md5"
expect_sync_failure missing-local-checksum
ecoda_write_checksum "${PAYLOAD}"

prepare_sync_run sync_malformed_checksum 4.0
cp "${PAYLOAD}.md5" "${PAYLOAD}.md5.saved"
printf 'MD5=not-a-digest\nSIZE=1\nPATH=%s\n' "${PAYLOAD}" > "${PAYLOAD}.md5"
expect_sync_failure malformed-local-checksum
mv -f "${PAYLOAD}.md5.saved" "${PAYLOAD}.md5"

RUNTIME_OUTPUT="${PAYLOAD}.runtime.json"
prepare_sync_run sync_malformed_runtime 5.0
chmod u+w "${RUNTIME_OUTPUT}"
cp "${RUNTIME_OUTPUT}" "${RUNTIME_OUTPUT}.saved"
printf '{"schema_version":1}\n' > "${RUNTIME_OUTPUT}"
ecoda_write_checksum "${RUNTIME_OUTPUT}"
expect_sync_failure malformed-runtime-schema
mv -f "${RUNTIME_OUTPUT}.saved" "${RUNTIME_OUTPUT}"
ecoda_write_checksum "${RUNTIME_OUTPUT}"
chmod a-w "${RUNTIME_OUTPUT}"

prepare_sync_run sync_missing_runtime 6.0
rm -f "${RUNTIME_OUTPUT}" "${RUNTIME_OUTPUT}.md5"
expect_sync_failure missing-runtime-record

echo "benchmark sync owner: OK"

LEGACY_ANALYSIS_ROOT="${TMP_DIR}/scratch/batch_effect/uncorrected"
LEGACY_ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/batch_effect/uncorrected"
mkdir -p "${LEGACY_ANALYSIS_ROOT}" "${LEGACY_ANALYSIS_NAS_ROOT}"
unset ANALYSIS_VARIANT
export ANALYSIS_PASS=uncorrected \
  ANALYSIS_ROOT="${LEGACY_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${LEGACY_ANALYSIS_NAS_ROOT}" \
  ANALYSIS_LOG_PREFIX="execution_times_batch_effect_uncorrected_"
benchmark_sync_artifacts_for Adams pseudobulk
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${LEGACY_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_pseudobulk.rds" ]]
benchmark_sync_artifacts_for Adams composition
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${LEGACY_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_composition.rds" ]]
[[ "${SYNC_ARTIFACTS[1]}" == \
  "${LEGACY_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_metadata.rds" ]]
_ecoda_stage5_artifacts_for Adams batch_effect_uncorrected pseudobulk
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${LEGACY_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_pseudobulk.rds" ]]

echo "legacy benchmark root/stem compatibility: OK"

FINAL_ANALYSIS_ROOT="${TMP_DIR}/scratch/batch_effect/uncorrected_final"
FINAL_ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/batch_effect/uncorrected_final"
mkdir -p "${FINAL_ANALYSIS_ROOT}" "${FINAL_ANALYSIS_NAS_ROOT}"
export ANALYSIS_VARIANT=final ANALYSIS_PASS=uncorrected \
  ANALYSIS_ROOT="${FINAL_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${FINAL_ANALYSIS_NAS_ROOT}" \
  ANALYSIS_LOG_PREFIX="execution_times_batch_effect_uncorrected_final_"

benchmark_sync_artifacts_for Adams prepare_pseudobulk
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/pseudobulks/Adams_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds" ]]
benchmark_sync_artifacts_for Adams pseudobulk
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_pseudobulk.rds" ]]
benchmark_sync_artifacts_for Adams gloscope
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_gloscope.rds" ]]
benchmark_sync_artifacts_for Adams composition
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_composition.rds" ]]
[[ "${SYNC_ARTIFACTS[1]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_metadata.rds" ]]
benchmark_sync_artifacts_for Adams mrvi
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather" ]]
benchmark_sync_artifacts_for Adams pilot
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_uncorrected_final_hvg2000_highres_pilot_dists.feather" ]]
benchmark_sync_artifacts_for Adams qot
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_uncorrected_final_hvg2000_highres_qot_dists.feather" ]]
for final_path in "${SYNC_ARTIFACTS[@]}"; do
  [[ "${final_path}" == *"_batch_effect_uncorrected_final_"* ]]
done

_ecoda_stage5_artifacts_for Adams batch_effect_uncorrected pseudobulk
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_pseudobulk.rds" ]]
_ecoda_stage5_artifacts_for Adams batch_effect_uncorrected composition
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_composition.rds" ]]
[[ "${ECODA_BENCHMARK_ARTIFACTS[1]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_metadata.rds" ]]

HISTORICAL_CORRECTED_ROOT="${TMP_DIR}/scratch/batch_effect/corrected_final/recovery_35row"
HISTORICAL_CORRECTED_NAS="${TMP_DIR}/nas/project/batch_effect/corrected_final/recovery_35row"
mkdir -p "${HISTORICAL_CORRECTED_ROOT}" "${HISTORICAL_CORRECTED_NAS}"
printf 'immutable historical evidence\n' > "${HISTORICAL_CORRECTED_NAS}/preserved"
export ANALYSIS_VARIANT=corrected_final ANALYSIS_PASS=corrected \
  ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION=recovery_35row \
  ANALYSIS_ROOT="${HISTORICAL_CORRECTED_ROOT}" \
  ANALYSIS_NAS_ROOT="${HISTORICAL_CORRECTED_NAS}" \
  ANALYSIS_LOG_PREFIX="execution_times_batch_effect_corrected_final_"
HISTORICAL_RSYNC_BEFORE="$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')"
if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "publisher accepted a historical corrected-final alias" >&2; exit 1
fi
[[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${HISTORICAL_RSYNC_BEFORE}" ]]
[[ "$(cat "${HISTORICAL_CORRECTED_NAS}/preserved")" == "immutable historical evidence" ]]
[[ ! -e "${HISTORICAL_CORRECTED_NAS}/checksums.md5" ]]
unset ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION
CORRECTED_FINAL_ANALYSIS_ROOT="${TMP_DIR}/scratch/batch_effect/corrected_final"
CORRECTED_FINAL_ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/batch_effect/corrected_final"
mkdir -p "${CORRECTED_FINAL_ANALYSIS_ROOT}" "${CORRECTED_FINAL_ANALYSIS_NAS_ROOT}"
export ANALYSIS_VARIANT=corrected_final ANALYSIS_PASS=corrected \
  ANALYSIS_ROOT="${CORRECTED_FINAL_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${CORRECTED_FINAL_ANALYSIS_NAS_ROOT}" \
  ANALYSIS_LOG_PREFIX="execution_times_batch_effect_corrected_final_"
benchmark_sync_artifacts_for Adams prepare_pseudobulk
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/pseudobulks/Adams_batch_effect_corrected_final_pseudobulk_hvg2000.rds" ]]
benchmark_sync_artifacts_for Adams pseudobulk
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_pseudobulk.rds" ]]
benchmark_sync_artifacts_for Adams gloscope
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_gloscope.rds" ]]
benchmark_sync_artifacts_for Adams composition
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_composition.rds" ]]
[[ "${SYNC_ARTIFACTS[1]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_metadata.rds" ]]
benchmark_sync_artifacts_for Adams mrvi
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather" ]]
benchmark_sync_artifacts_for Adams pilot
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather" ]]
benchmark_sync_artifacts_for Adams qot
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_qot_dists.feather" ]]
for corrected_path in "${SYNC_ARTIFACTS[@]}"; do
  [[ "${corrected_path}" == *"_batch_effect_corrected_final_"* ]]
done
_ecoda_stage5_artifacts_for Adams batch_effect_corrected pseudobulk
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_pseudobulk.rds" ]]
_ecoda_stage5_artifacts_for Adams batch_effect_corrected composition
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_composition.rds" ]]
[[ "${ECODA_BENCHMARK_ARTIFACTS[1]}" == \
  "${CORRECTED_FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_metadata.rds" ]]
export ANALYSIS_VARIANT=corrected_final ANALYSIS_PASS=uncorrected \
  ANALYSIS_ROOT="${CORRECTED_FINAL_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${CORRECTED_FINAL_ANALYSIS_NAS_ROOT}"
if benchmark_sync_artifacts_for Adams batch_effect_uncorrected pseudobulk >/dev/null 2>&1; then
  echo "corrected-final helper accepted an uncorrected root/pass" >&2
  exit 1
fi
export ANALYSIS_VARIANT=final ANALYSIS_PASS=corrected \
  ANALYSIS_ROOT="${FINAL_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${FINAL_ANALYSIS_NAS_ROOT}"
if benchmark_sync_artifacts_for Adams batch_effect_corrected pseudobulk >/dev/null 2>&1; then
  echo "final helper accepted corrected pass" >&2
  exit 1
fi
export ANALYSIS_VARIANT=final ANALYSIS_PASS=uncorrected \
  ANALYSIS_ROOT="${FINAL_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${FINAL_ANALYSIS_NAS_ROOT}" \
  ANALYSIS_LOG_PREFIX="execution_times_batch_effect_uncorrected_final_"

_ecoda_stage5_artifacts_for Adams batch_effect_uncorrected pseudobulk
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_pseudobulk.rds" ]]
_ecoda_stage5_artifacts_for Adams batch_effect_uncorrected composition
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_composition.rds" ]]
[[ "${ECODA_BENCHMARK_ARTIFACTS[1]}" == \
  "${FINAL_ANALYSIS_ROOT}/results/Adams_batch_effect_uncorrected_final_metadata.rds" ]]

echo "final benchmark variant paths: OK"
export ANALYSIS_VARIANT=corrected_final ANALYSIS_PASS=corrected
MATRIX_CORRECTED_ANALYSIS_ROOT="${TMP_DIR}/scratch/batch_effect/corrected_final"
MATRIX_CORRECTED_ANALYSIS_NAS_ROOT="${TMP_DIR}/nas/project/batch_effect/corrected_final"
mkdir -p "${MATRIX_CORRECTED_ANALYSIS_ROOT}" \
  "${MATRIX_CORRECTED_ANALYSIS_NAS_ROOT}"
unset ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION
export ANALYSIS_ROOT="${MATRIX_CORRECTED_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${MATRIX_CORRECTED_ANALYSIS_NAS_ROOT}"
benchmark_sync_artifacts_for Adams pseudobulk
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${MATRIX_CORRECTED_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_pseudobulk.rds" ]]
_ecoda_stage5_artifacts_for Adams batch_effect_corrected composition
[[ "${ECODA_BENCHMARK_ARTIFACTS[0]}" == \
  "${MATRIX_CORRECTED_ANALYSIS_ROOT}/results/Adams_batch_effect_corrected_final_composition.rds" ]]
MATRIX_METHOD_FILTER_SOURCE="${TMP_DIR}/matrix-method-filter.tsv"
MATRIX_METHOD_FILTER="${ECODA_RUN_ROOT}/manifests/method_matrix.tsv"
printf 'Adams\tbatch_effect_corrected\tmrvi\n' > "${MATRIX_METHOD_FILTER_SOURCE}"
ecoda_atomic_install_manifest "${MATRIX_METHOD_FILTER_SOURCE}" \
  "${MATRIX_METHOD_FILTER}" 3
ecoda_write_checksum "${MATRIX_METHOD_FILTER}"
MATRIX_METHOD_FILTER_MD5="${ECODA_CHECKSUM_MD5}"
MATRIX_METHOD_FILTER_SIZE="${ECODA_CHECKSUM_SIZE}"
MATRIX_METHOD_FILTER_SHA256="$(ecoda_sha256_file "${MATRIX_METHOD_FILTER}")"
ecoda_validate_manifest "${MATRIX_METHOD_FILTER}" 3
ecoda_validate_checksum "${MATRIX_METHOD_FILTER}"
export ECODA_STAGE5_METHOD_MATRIX="${MATRIX_METHOD_FILTER}" \
  METHOD_MATRIX="${MATRIX_METHOD_FILTER}" \
  METHOD_MATRIX_SOURCE_PATH="${MATRIX_METHOD_FILTER_SOURCE}" \
  METHOD_MATRIX_MD5="${MATRIX_METHOD_FILTER_MD5}" \
  METHOD_MATRIX_SIZE="${MATRIX_METHOD_FILTER_SIZE}" \
  METHOD_MATRIX_SHA256="${MATRIX_METHOD_FILTER_SHA256}" \
  METHOD_MATRIX_IDENTITY="${MATRIX_METHOD_FILTER_SHA256}" \
  METHOD_MATRIX_DECLARED_COUNT=1 ECODA_STAGE5_METHOD_MATRIX_MODE=1
benchmark_sync_artifacts_for Adams mrvi
[[ "${SYNC_ARTIFACTS[0]}" == \
  "${MATRIX_CORRECTED_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather" ]]
benchmark_sync_artifacts_for Adams pilot
[[ "${#SYNC_ARTIFACTS[@]}" == 0 ]]
unset ECODA_STAGE5_METHOD_MATRIX METHOD_MATRIX METHOD_MATRIX_SOURCE_PATH \
  METHOD_MATRIX_MD5 METHOD_MATRIX_SIZE METHOD_MATRIX_SHA256 METHOD_MATRIX_IDENTITY \
  METHOD_MATRIX_DECLARED_COUNT ECODA_STAGE5_METHOD_MATRIX_MODE
unset ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION
export ANALYSIS_VARIANT=final ANALYSIS_PASS=uncorrected \
  ANALYSIS_ROOT="${FINAL_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${FINAL_ANALYSIS_NAS_ROOT}"

FINAL_PAYLOAD="${FINAL_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather"
mkdir -p "$(dirname "${FINAL_PAYLOAD}")"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"s1":[1.0,0.0],"s2":[0.0,1.0]},index=["s1","s2"]).to_feather(sys.argv[1])' \
  "${FINAL_PAYLOAD}"
ecoda_write_checksum "${FINAL_PAYLOAD}"
FINAL_PAYLOAD_MD5="${ECODA_CHECKSUM_MD5}"
FINAL_RUNTIME="${FINAL_PAYLOAD}.runtime.json"
printf '{"schema_version":1,"artifact_path":"%s","artifact_md5":"%s","dataset":"Adams","method":"MrVI_hvg2000","time_secs":1.0,"mem_GB":1.0}\n' \
  "${FINAL_PAYLOAD}" "${FINAL_PAYLOAD_MD5}" > "${FINAL_RUNTIME}"
ecoda_write_checksum "${FINAL_RUNTIME}"

FINAL_METADATA_OUTPUT="${FINAL_ANALYSIS_ROOT}/metadata/Adams_sample_metadata.feather"
mkdir -p "$(dirname "${FINAL_METADATA_OUTPUT}")"
printf 'Sample\ns1\ns2\n' > "${FINAL_METADATA_OUTPUT}"
ecoda_write_checksum "${FINAL_METADATA_OUTPUT}"

FINAL_RUN_ID="final_sync_run"
ecoda_init_run stage5 "${FINAL_RUN_ID}" >/dev/null
export ECODA_RUN_ID ECODA_RUN_ROOT
FINAL_SELECTION_MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'Adams\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n' \
  > "${FINAL_SELECTION_MANIFEST}"
export ECODA_SELECTION_MANIFEST="${FINAL_SELECTION_MANIFEST}" ECODA_EXACT_SELECTION=0
ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}"
DATASET_NAMES=(Adams)
LABELS=(mrvi)
export DATASET_NAMES LABELS
export EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs"
mkdir -p "${EXECUTION_LOG_DIR}"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[1.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_uncorrected_final_mrvi_Adams.feather"
ecoda_write_checksum \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_uncorrected_final_mrvi_Adams.feather"
ecoda_write_artifact_record \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_uncorrected_final_mrvi_Adams.feather" \
  stage5_execution_log "${ECODA_RUN_ID}" >/dev/null
METADATA_EXPORT_MANIFEST="${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
printf 'Adams\tbatch_effect_uncorrected\t/immutable/Adams.h5ad\t%s\n' \
  "${FINAL_METADATA_OUTPUT}" > "${METADATA_EXPORT_MANIFEST}"
grep -q $'^Adams\tbatch_effect_uncorrected\t/immutable/Adams.h5ad\t'"${FINAL_METADATA_OUTPUT}"'$' \
  "${METADATA_EXPORT_MANIFEST}"

analysis_merge_sync_cleanup mrvi
FINAL_SYNC_FILES="${ECODA_RUN_ROOT}/manifests/sync_files.tsv"
grep -q '^embeddings/Adams_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather$' \
  "${FINAL_SYNC_FILES}"
grep -q '^embeddings/Adams_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather.md5$' \
  "${FINAL_SYNC_FILES}"
grep -q '^metadata/Adams_sample_metadata.feather$' "${FINAL_SYNC_FILES}"
grep -q '^metadata/Adams_sample_metadata.feather.md5$' "${FINAL_SYNC_FILES}"
if grep -Eq 'uncorrected_(hvg|pseudobulk|gloscope|composition)' \
    "${FINAL_SYNC_FILES}"; then
  echo "final sync expanded a legacy analysis path" >&2
  exit 1
fi
[[ -s "${FINAL_ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather" ]]
[[ -s "${FINAL_ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather.md5" ]]
[[ -s "${FINAL_ANALYSIS_NAS_ROOT}/embeddings/Adams_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather" ]]
[[ ! -e "${FINAL_ANALYSIS_NAS_ROOT}/pseudobulks/Adams_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds" ]]
echo "final benchmark sync variant: OK"
export ANALYSIS_VARIANT=corrected_final ANALYSIS_PASS=corrected \
  ANALYSIS_ROOT="${MATRIX_CORRECTED_ANALYSIS_ROOT}" \
  ANALYSIS_NAS_ROOT="${MATRIX_CORRECTED_ANALYSIS_NAS_ROOT}" \
  ANALYSIS_LOG_PREFIX="execution_times_batch_effect_corrected_final_"
CORRECTED_PAYLOAD="${MATRIX_CORRECTED_ANALYSIS_ROOT}/embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather"
mkdir -p "$(dirname "${CORRECTED_PAYLOAD}")"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"s1":[1.0,0.0],"s2":[0.0,1.0]},index=["s1","s2"]).to_feather(sys.argv[1])' \
  "${CORRECTED_PAYLOAD}"
ecoda_write_checksum "${CORRECTED_PAYLOAD}"
CORRECTED_PAYLOAD_MD5="${ECODA_CHECKSUM_MD5}"
CORRECTED_RUNTIME="${CORRECTED_PAYLOAD}.runtime.json"
printf '{"schema_version":1,"artifact_path":"%s","artifact_md5":"%s","dataset":"Adams","method":"MrVI_hvg2000","time_secs":1.0,"mem_GB":1.0}\n' \
  "${CORRECTED_PAYLOAD}" "${CORRECTED_PAYLOAD_MD5}" > "${CORRECTED_RUNTIME}"
ecoda_write_checksum "${CORRECTED_RUNTIME}"
CORRECTED_METADATA_OUTPUT="${MATRIX_CORRECTED_ANALYSIS_ROOT}/metadata/Adams_sample_metadata.feather"
mkdir -p "$(dirname "${CORRECTED_METADATA_OUTPUT}")"
printf 'Sample\ns1\ns2\n' > "${CORRECTED_METADATA_OUTPUT}"
ecoda_write_checksum "${CORRECTED_METADATA_OUTPUT}"

# A failed direct-root synchronization leaves a FAIL owner that a later
# receipt-authorized publisher can reclaim without changing the payload.
CORRECTED_FAILED_RUN_ID="corrected_final_sync_failed"
ecoda_init_run stage5 "${CORRECTED_FAILED_RUN_ID}" >/dev/null
export ECODA_RUN_ID="${CORRECTED_FAILED_RUN_ID}" ECODA_RUN_ROOT
CORRECTED_FAILED_SELECTION="${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'Adams\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${CORRECTED_FAILED_SELECTION}"
export ECODA_SELECTION_MANIFEST="${CORRECTED_FAILED_SELECTION}" \
  ECODA_EXACT_SELECTION=0
ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}"
DATASET_NAMES=(Adams)
LABELS=(mrvi)
export DATASET_NAMES LABELS
EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs"
export EXECUTION_LOG_DIR
CORRECTED_FAILED_METADATA_EXPORT="${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
printf 'Adams\tbatch_effect_corrected\t/immutable/Adams-corrected.h5ad\t%s\n' \
  "${CORRECTED_METADATA_OUTPUT}" > "${CORRECTED_FAILED_METADATA_EXPORT}"
seal_publication_receipt
CORRECTED_PAYLOAD_MD5_SAVED="${CORRECTED_PAYLOAD}.md5.failed-recovery"
cp "${CORRECTED_PAYLOAD}.md5" "${CORRECTED_PAYLOAD_MD5_SAVED}"
rm -f "${CORRECTED_PAYLOAD}.md5"
if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "failed corrected-final synchronization unexpectedly succeeded" >&2
  exit 1
fi
CORRECTED_FAILED_SYNC_OWNER="$(ecoda_owner_dir stage5 \
  "sync/${MATRIX_CORRECTED_ANALYSIS_ROOT}")"
[[ ! -e "${CORRECTED_FAILED_SYNC_OWNER}/owner" ]]
mv -f "${CORRECTED_PAYLOAD_MD5_SAVED}" "${CORRECTED_PAYLOAD}.md5"

CORRECTED_RUN_ID="corrected_final_sync_run"
ecoda_init_run stage5 "${CORRECTED_RUN_ID}" >/dev/null
export ECODA_RUN_ID ECODA_RUN_ROOT
CORRECTED_SELECTION_MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'Adams\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${CORRECTED_SELECTION_MANIFEST}"
export ECODA_SELECTION_MANIFEST="${CORRECTED_SELECTION_MANIFEST}" ECODA_EXACT_SELECTION=0
ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}"
DATASET_NAMES=(Adams)
LABELS=(mrvi)
export DATASET_NAMES LABELS
export EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs"
mkdir -p "${EXECUTION_LOG_DIR}"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[1.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_corrected_final_mrvi_Adams.feather"
ecoda_write_checksum \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_corrected_final_mrvi_Adams.feather"
ecoda_write_artifact_record \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_corrected_final_mrvi_Adams.feather" \
  stage5_execution_log "${ECODA_RUN_ID}" >/dev/null
CORRECTED_METADATA_EXPORT_MANIFEST="${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
printf 'Adams\tbatch_effect_corrected\t/immutable/Adams-corrected.h5ad\t%s\n' \
  "${CORRECTED_METADATA_OUTPUT}" > "${CORRECTED_METADATA_EXPORT_MANIFEST}"
seal_publication_receipt

# Replacing a child with a newly self-checksummed receipt cannot change the
# child bytes authorized by the composite (whose NAS copy remains sealed).
BOUND_PROMOTION="$(jq -er '.source_root' "${ECODA_PRECOMPUTE_BACKUP_RECEIPT}")"
BOUND_CHILD="${BOUND_PROMOTION}/manifests/inputs.receipt.json"
cp "${BOUND_CHILD}" "${TMP_DIR}/bound-child.saved"
BOUND_OWNER="$(ecoda_owner_dir stage5 "sync/${ANALYSIS_ROOT}")/owner"
[[ ! -e "${BOUND_OWNER}" && ! -L "${BOUND_OWNER}" ]]
BOUND_RSYNC_BEFORE="$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')"
pixi run python - "${BOUND_CHILD}" <<'PY'
import json, pathlib, sys
p = pathlib.Path(sys.argv[1])
r = json.loads(p.read_text())
r["completed_at"] = "2026-10-21T01:00:00Z"
p.write_text(json.dumps(r) + "\n")
PY
ecoda_write_checksum "${BOUND_CHILD}"
if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "publication accepted replaced self-checksummed child receipt" >&2; exit 1
fi
[[ ! -e "${BOUND_OWNER}" && ! -L "${BOUND_OWNER}" ]]
[[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${BOUND_RSYNC_BEFORE}" ]]
cp "${TMP_DIR}/bound-child.saved" "${BOUND_CHILD}"
ecoda_write_checksum "${BOUND_CHILD}"

# A new self-sidecar and fresh in-process SHA cannot authorize either a changed
# compute copy or a replaced promotion source outside the composite inventory.
for changed_expected in compute promotion; do
  printf 'Adams\tbatch_effect_corrected\tgloscope\n' > "${ECODA_EXPECTED_PENDING}"
  ecoda_write_checksum "${ECODA_EXPECTED_PENDING}"
  export ECODA_EXPECTED_PENDING_SHA256="$(sha256_file "${ECODA_EXPECTED_PENDING}")"
  if [[ "${changed_expected}" == promotion ]]; then
    cp "${ECODA_EXPECTED_PENDING}" "${BOUND_PROMOTION}/manifests/expected_pending.tsv"
    ecoda_write_checksum "${BOUND_PROMOTION}/manifests/expected_pending.tsv"
  fi
  if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
    echo "publication accepted changed ${changed_expected} expected pending with fresh sidecar" >&2; exit 1
  fi
  [[ ! -e "${BOUND_OWNER}" && ! -L "${BOUND_OWNER}" ]]
  [[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${BOUND_RSYNC_BEFORE}" ]]
done
: > "${ECODA_EXPECTED_PENDING}"
: > "${BOUND_PROMOTION}/manifests/expected_pending.tsv"
stage5_write_pending_checksum "${ECODA_EXPECTED_PENDING}"
stage5_write_pending_checksum "${BOUND_PROMOTION}/manifests/expected_pending.tsv"
export ECODA_EXPECTED_PENDING_SHA256="$(sha256_file "${ECODA_EXPECTED_PENDING}")"

analysis_merge_sync_cleanup mrvi
CORRECTED_SYNC_FILES="${ECODA_RUN_ROOT}/manifests/sync_files.tsv"
grep -q '^embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather$' \
  "${CORRECTED_SYNC_FILES}"
grep -q '^embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather.md5$' \
  "${CORRECTED_SYNC_FILES}"
grep -q '^metadata/Adams_sample_metadata.feather$' "${CORRECTED_SYNC_FILES}"
grep -q '^metadata/Adams_sample_metadata.feather.md5$' "${CORRECTED_SYNC_FILES}"
[[ -s "${MATRIX_CORRECTED_ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather" ]]
[[ -s "${MATRIX_CORRECTED_ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather.md5" ]]
[[ -s "${MATRIX_CORRECTED_ANALYSIS_NAS_ROOT}/embeddings/Adams_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather" ]]
[[ ! -e "${MATRIX_CORRECTED_ANALYSIS_NAS_ROOT}/pseudobulks/Adams_batch_effect_corrected_final_pseudobulk_hvg2000.rds" ]]
echo "corrected-final benchmark sync variant: OK"
CORRECTED_SYNC_OWNER="$(ecoda_owner_dir stage5 \
  "sync/${MATRIX_CORRECTED_ANALYSIS_ROOT}")"
[[ "$(ecoda_owner_state "${CORRECTED_SYNC_OWNER}")" == "OK" ]]
[[ "$(ecoda_owner_run "${CORRECTED_SYNC_OWNER}")" == \
  "corrected_final_sync_run" ]]
CORRECTED_SYNC_OWNER_SNAPSHOT="${TMP_DIR}/corrected_sync_owner.before_conflict"
cp "${CORRECTED_SYNC_OWNER}/owner" "${CORRECTED_SYNC_OWNER_SNAPSHOT}"
CORRECTED_SYNC_OWNER_STATE="$(ecoda_owner_state "${CORRECTED_SYNC_OWNER}")"
CORRECTED_SYNC_OWNER_RUN="$(ecoda_owner_run "${CORRECTED_SYNC_OWNER}")"
CORRECTED_SYNC_OWNER_BYTES="$(wc -c < "${CORRECTED_SYNC_OWNER}/owner" | tr -d '[:space:]')"


# Another run may validate identical publication without rewriting owner/NAS.
CORRECTED_CONFLICT_RUN_ID="corrected_final_sync_conflict"
ecoda_init_run stage5 "${CORRECTED_CONFLICT_RUN_ID}" >/dev/null
export ECODA_RUN_ID="${CORRECTED_CONFLICT_RUN_ID}" ECODA_RUN_ROOT
CORRECTED_CONFLICT_SELECTION="${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'Adams\tbatch_effect_corrected\tbatch_effect_corrected\n' \
  > "${CORRECTED_CONFLICT_SELECTION}"
export ECODA_SELECTION_MANIFEST="${CORRECTED_CONFLICT_SELECTION}" \
  ECODA_EXACT_SELECTION=0
ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}"
DATASET_NAMES=(Adams)
LABELS=(mrvi)
export DATASET_NAMES LABELS
export EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs"
mkdir -p "${EXECUTION_LOG_DIR}"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000"],"time_secs":[1.0],"mem_GB":[1.0]}).to_feather(sys.argv[1])' \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_corrected_final_mrvi_Adams.feather"
ecoda_write_checksum \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_corrected_final_mrvi_Adams.feather"
ecoda_write_artifact_record \
  "${EXECUTION_LOG_DIR}/execution_times_batch_effect_corrected_final_mrvi_Adams.feather" \
  stage5_execution_log "${ECODA_RUN_ID}" >/dev/null
CORRECTED_CONFLICT_METADATA_EXPORT="${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
printf 'Adams\tbatch_effect_corrected\t/immutable/Adams-corrected.h5ad\t%s\n' \
  "${CORRECTED_METADATA_OUTPUT}" > "${CORRECTED_CONFLICT_METADATA_EXPORT}"
seal_publication_receipt
REMOTE_CHECKSUM_SNAPSHOT="${TMP_DIR}/checksums.before_reuse"
cp "${ANALYSIS_NAS_ROOT}/checksums.md5" "${REMOTE_CHECKSUM_SNAPSHOT}"
RSYNC_BEFORE="$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')"
analysis_merge_sync_cleanup mrvi
cmp "${CORRECTED_SYNC_OWNER}/owner" "${CORRECTED_SYNC_OWNER_SNAPSHOT}"
cmp "${ANALYSIS_NAS_ROOT}/checksums.md5" "${REMOTE_CHECKSUM_SNAPSHOT}"
[[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${RSYNC_BEFORE}" ]]

# Missing selection must fail closed despite populated Cartesian globals.
SAVED_SELECTION="${ECODA_SELECTION_MANIFEST}"
unset ECODA_SELECTION_MANIFEST
DATASET_NAMES=(Adams Breast_cancer)
LABELS=(mrvi pilot qot)
if analysis_merge_sync_cleanup "${LABELS[@]}" >/dev/null 2>&1; then
  echo "missing sync selection widened publication" >&2; exit 1
fi
export ECODA_SELECTION_MANIFEST="${SAVED_SELECTION}"
DATASET_NAMES=(Adams)
LABELS=(mrvi)
cmp "${CORRECTED_SYNC_OWNER}/owner" "${CORRECTED_SYNC_OWNER_SNAPSHOT}"
[[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${RSYNC_BEFORE}" ]]

# A changed publication cannot reclaim terminal OK without the sealed archive.
SAVED_RECEIPT="${ECODA_PRECOMPUTE_BACKUP_RECEIPT}"
# Deliberately retain size and mtime while changing content: rsync quick-check
# must not skip this selected replacement.
pixi run python - "${CORRECTED_METADATA_OUTPUT}" \
  "${ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather" <<'PY'
import os, pathlib, sys
local, remote = map(pathlib.Path, sys.argv[1:])
old = remote.read_bytes()
new = bytes([old[0] ^ 1]) + old[1:]
local.write_bytes(new)
stat = remote.stat()
os.utime(local, ns=(stat.st_atime_ns, stat.st_mtime_ns))
assert len(new) == len(old) and new != old
assert local.stat().st_mtime_ns == remote.stat().st_mtime_ns
PY
ecoda_write_checksum "${CORRECTED_METADATA_OUTPUT}"
unset ECODA_PRECOMPUTE_BACKUP_RECEIPT
if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "changed publication reclaimed terminal OK without receipt" >&2; exit 1
fi
cmp "${CORRECTED_SYNC_OWNER}/owner" "${CORRECTED_SYNC_OWNER_SNAPSHOT}"
cmp "${ANALYSIS_NAS_ROOT}/checksums.md5" "${REMOTE_CHECKSUM_SNAPSHOT}"
[[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${RSYNC_BEFORE}" ]]
export ECODA_PRECOMPUTE_BACKUP_RECEIPT="${SAVED_RECEIPT}"
# Retry the same successful compute run with the same sealed archive after a
# failed first publication; no recomputation or new authorization is needed.
if FAIL_RSYNC=1 analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "injected NAS transfer failure unexpectedly succeeded" >&2; exit 1
fi
[[ "$(ecoda_owner_state "${CORRECTED_SYNC_OWNER}")" == FAIL ]]
[[ "$(ecoda_owner_run "${CORRECTED_SYNC_OWNER}")" == "${CORRECTED_CONFLICT_RUN_ID}" ]]
cmp "${ANALYSIS_NAS_ROOT}/checksums.md5" "${REMOTE_CHECKSUM_SNAPSHOT}"
USE_REAL_RSYNC=1 analysis_merge_sync_cleanup mrvi
[[ "$(ecoda_owner_state "${CORRECTED_SYNC_OWNER}")" == OK ]]
[[ "$(ecoda_owner_run "${CORRECTED_SYNC_OWNER}")" == "${CORRECTED_CONFLICT_RUN_ID}" ]]
cmp "${CORRECTED_METADATA_OUTPUT}" "${ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather"
ARCHIVED_OWNER="$(pixi run python -c 'import csv,sys; print(next(csv.DictReader(open(sys.argv[1]), delimiter="\t"))["archive_path"])' \
  "${ECODA_RUN_ROOT}/promotion/manifests/global_sync_owner.tsv")"
cmp "${ARCHIVED_OWNER}" "${CORRECTED_SYNC_OWNER_SNAPSHOT}"

# A later publication fails after replacing a selected payload but before its
# sidecar/checksum manifest. Retry uses the same seal and original continuity.
PARTIAL_SELECTION="${ECODA_SELECTION_MANIFEST}"
ecoda_init_run stage5 corrected_final_sync_partial >/dev/null
export ECODA_RUN_ID ECODA_RUN_ROOT
export ECODA_SELECTION_MANIFEST="${ECODA_RUN_ROOT}/manifests/selection.tsv"
cp "${PARTIAL_SELECTION}" "${ECODA_SELECTION_MANIFEST}"
ecoda_write_checksum "${ECODA_SELECTION_MANIFEST}"
export EXECUTION_LOG_DIR="${ECODA_RUN_ROOT}/logs"
mkdir -p "${EXECUTION_LOG_DIR}"
PARTIAL_LOG="${EXECUTION_LOG_DIR}/execution_times_batch_effect_corrected_final_mrvi_Adams.feather"
pixi run python -c 'import pandas as pd,sys; pd.DataFrame({"dataset":["Adams"],"method":["MrVI_hvg2000_partial"],"time_secs":[7.0],"mem_GB":[2.0]}).to_feather(sys.argv[1])' \
  "${PARTIAL_LOG}"
ecoda_write_checksum "${PARTIAL_LOG}"
ecoda_write_artifact_record "${PARTIAL_LOG}" stage5_execution_log "${ECODA_RUN_ID}" >/dev/null
printf 'Adams\tbatch_effect_corrected\t/immutable/Adams-corrected.h5ad\t%s\n' \
  "${CORRECTED_METADATA_OUTPUT}" > "${ECODA_RUN_ROOT}/manifests/metadata_export.tsv"
printf 'unselected continuity\n' > "${ANALYSIS_NAS_ROOT}/unselected.txt"
chmod 666 "${ANALYSIS_NAS_ROOT}/unselected.txt"
printf '%s  unselected.txt\n' "$(ecoda_md5_file "${ANALYSIS_NAS_ROOT}/unselected.txt")" \
  >> "${ANALYSIS_NAS_ROOT}/checksums.md5"
seal_publication_receipt
PARTIAL_RECEIPT_SHA="${ECODA_PRECOMPUTE_BACKUP_RECEIPT_SHA256}"
printf 'Sample\npartial-retry\n' > "${CORRECTED_METADATA_OUTPUT}"
ecoda_write_checksum "${CORRECTED_METADATA_OUTPUT}"
if FAIL_RSYNC_AFTER_METADATA=1 analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "mid-transfer failure unexpectedly succeeded" >&2; exit 1
fi
[[ "$(ecoda_owner_state "${CORRECTED_SYNC_OWNER}")" == FAIL ]]
cmp "${CORRECTED_METADATA_OUTPUT}" "${ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather"
if (cd "${ANALYSIS_NAS_ROOT}" && md5sum -c checksums.md5 >/dev/null 2>&1); then
  echo "partial-transfer fixture did not invalidate prior canonical checksum" >&2; exit 1
fi
# Advance the same transition through the aggregate payload but not its sidecar.
if FAIL_RSYNC_AFTER_LOG=1 analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "post-log transfer failure unexpectedly succeeded" >&2; exit 1
fi
[[ "$(ecoda_owner_state "${CORRECTED_SYNC_OWNER}")" == FAIL ]]
if ecoda_validate_checksum_remote \
    "${ANALYSIS_NAS_ROOT}/embeddings/execution_times.feather" \
    "${ANALYSIS_NAS_ROOT}/embeddings/execution_times.feather.md5" >/dev/null 2>&1; then
  echo "partial aggregate log fixture retained a matching sidecar" >&2; exit 1
fi
PARTIAL_RSYNC_BEFORE="$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')"
cp "${CORRECTED_SYNC_OWNER}/owner" "${TMP_DIR}/partial-owner.saved"
printf 'unexpected unselected bytes\n' > "${ANALYSIS_NAS_ROOT}/unselected.txt"
if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "partial retry accepted changed unselected canonical file" >&2; exit 1
fi
printf 'unselected continuity\n' > "${ANALYSIS_NAS_ROOT}/unselected.txt"
printf 'neither original nor intended\n' > "${ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather"
if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "partial retry accepted unrecognized selected bytes" >&2; exit 1
fi
cp "${CORRECTED_METADATA_OUTPUT}" "${ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather"
cmp "${CORRECTED_SYNC_OWNER}/owner" "${TMP_DIR}/partial-owner.saved"
[[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${PARTIAL_RSYNC_BEFORE}" ]]
analysis_merge_sync_cleanup mrvi
[[ "$(ecoda_owner_state "${CORRECTED_SYNC_OWNER}")" == OK ]]
[[ "$(ecoda_owner_run "${CORRECTED_SYNC_OWNER}")" == "${ECODA_RUN_ID}" ]]
[[ "$(sha256_file "${ECODA_PRECOMPUTE_BACKUP_RECEIPT}")" == "${PARTIAL_RECEIPT_SHA}" ]]
cmp "${CORRECTED_METADATA_OUTPUT}" "${ANALYSIS_NAS_ROOT}/metadata/Adams_sample_metadata.feather"
pixi run python - "${ANALYSIS_NAS_ROOT}/embeddings/execution_times.feather" <<'PY'
import pandas as pd, sys
frame = pd.read_feather(sys.argv[1])
assert {"MrVI_hvg2000", "MrVI_hvg2000_partial"} <= set(frame["method"])
assert len(frame[frame["method"] == "MrVI_hvg2000_partial"]) == 1
PY

# Even identical publication never bypasses an ACTIVE global writer.
ecoda_owner_set_state "${CORRECTED_SYNC_OWNER}" ACTIVE "active writer fixture"
seal_publication_receipt
cp "${CORRECTED_SYNC_OWNER}/owner" "${TMP_DIR}/active-owner.saved"
RSYNC_BEFORE="$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')"
if analysis_merge_sync_cleanup mrvi >/dev/null 2>&1; then
  echo "ACTIVE global sync owner was admitted" >&2; exit 1
fi
cmp "${CORRECTED_SYNC_OWNER}/owner" "${TMP_DIR}/active-owner.saved"
[[ "$(wc -l < "${RSYNC_CAPTURE}" | tr -d '[:space:]')" == "${RSYNC_BEFORE}" ]]
echo "corrected-final receipt reclaim, identical reuse and ACTIVE exclusion: OK"
