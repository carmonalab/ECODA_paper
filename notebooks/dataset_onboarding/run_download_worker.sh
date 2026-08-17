#!/bin/bash
#SBATCH --job-name=onboard_download
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
# ==============================================================================
# run_download_worker.sh -- one Phase-5 dataset download task (SLURM array)
# ==============================================================================
# Downloads one PILOT-GM-VAE study file into ${HPC_SCRATCH_DIR}/_downloads/
# (BeeGFS scratch via the $HOME/scratch symlink: 1.1 PB, no per-user size
# quota, NOT $HOME itself -- SSD, 1 TB quota, backed up). Never touches the
# NAS: the login-node submitter (download_datasets_hpc.sh) rsyncs after all
# tasks pass.
#
# * One task per key (8 keys). The submitter exports DOWNLOAD_KEYS
#   (space-separated, array-task 1-8 indexed); the key comes from
#   SLURM_ARRAY_TASK_ID or --key <key> for direct login-node fallback runs
#   (nice + DOWNLOAD_LIMIT_RATE).
# * The URL/md5/tar-pattern catalog lives in download_sources.sh (single
#   source of truth shared with the Mac->NAS fallback download_datasets.sh).
#   "verify" md5s (CellxGene) fetch the `.h5ad.md5` sidecar and ENFORCE it.
# * `curl -L -C -` is resumable -- a killed/requeued task resumes from the
#   partial file on re-submission; already-verified files are skipped.
# * Tar entries do SELECTIVE extraction of only the needed h5ads (excluded
#   members -- MI(1)/PDAC/Kidney cancer/follicular lymphoma -- are never
#   extracted), then delete the archive. Member matching greps the listing
#   already written to _status/<key>.tar_listing (no repeated multi-GB
#   tar -tzf decompressions).
# * Per-task status (OK/FAIL + md5s + size) is written atomically to
#   ${HPC_SCRATCH_DIR}/_downloads/_status/<key>.status for the submitter tail;
#   tar-key statuses also carry the extracted files' md5s (MD5_<file>=) so the
#   tail verifies the NAS copies against worker-verified values.
# * A script-level EXIT trap removes any leftover _tar_tmp.* extraction dirs
#   (success or failure -- set -e exits skip the in-function rm otherwise).
#
# Usage:
#   sbatch ... run_download_worker.sh          # key from SLURM_ARRAY_TASK_ID + DOWNLOAD_KEYS
#   bash run_download_worker.sh --key breast   # direct (login-node fallback)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch the script is copied to /var/spool/slurmd/job<id>/slurm_script,
# so BASH_SOURCE no longer points at the repo. Recover the real path from the
# job record (Command= field); BASH_SOURCE fallback for login-node execution.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../../src/slurm_config.sh"
cd "${PROJECT_ROOT}"

DOWNLOAD_DIR="${HPC_SCRATCH_DIR}/_downloads"
STATUS_DIR="${DOWNLOAD_DIR}/_status"
mkdir -p "${DOWNLOAD_DIR}" "${STATUS_DIR}"

# Clean any leftover extraction temp dirs on exit (also covers _tar_tmp.* leaks
# from previously failed tasks).
trap 'rm -rf "${DOWNLOAD_DIR}"/_tar_tmp.* 2>/dev/null || true' EXIT

# --- key dispatch --------------------------------------------------------------
KEY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --key) KEY="${2:-}"; shift 2 ;;
    *) echo "ERROR: unknown argument $1 (only --key <key> is accepted)." >&2; exit 1 ;;
  esac
done
if [[ -z "${KEY}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  read -r -a KEYS <<<"${DOWNLOAD_KEYS:-}"
  KEY="${KEYS[$((SLURM_ARRAY_TASK_ID - 1))]:-}"
fi
if [[ -z "${KEY}" ]]; then
  echo "ERROR: cannot resolve download key (set --key <key> or export DOWNLOAD_KEYS for array runs)." >&2
  exit 1
fi

source "${SCRIPT_DIR}/download_sources.sh"
if ! source_download_key "${KEY}"; then
  echo "ERROR: unknown --only key '${KEY}'." >&2
  exit 1
fi

LIMIT_RATE="${DOWNLOAD_LIMIT_RATE:-}"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${KEY}: $*"; }

# --- helpers -------------------------------------------------------------------
md5_of() {
  md5sum "$1" | awk '{print $1}'
}

# write_status <status> [extra_kv...] -- atomic per-key status file; the
# submitter tail reads these to build the download_log.md report.
write_status() {
  local st="$1"
  shift
  {
    echo "STATUS=${st}"
    echo "KEY=${KEY}"
    echo "UPDATED_AT=$(date -Is)"
    [[ -n "${EXPECTED_MD5:-}" ]] && echo "MD5_EXPECTED=${EXPECTED_MD5}"
    [[ -n "${ACTUAL_MD5:-}" ]] && echo "MD5_ACTUAL=${ACTUAL_MD5}"
    if [[ -n "${TARGET_FILE:-}" ]]; then
      echo "FILE=${TARGET_FILE}"
      if [[ -f "${DOWNLOAD_DIR}/${TARGET_FILE}" ]]; then
        echo "SIZE=$(stat -c %s "${DOWNLOAD_DIR}/${TARGET_FILE}" 2>/dev/null || echo '?')"
      fi
    fi
    for kv in "$@"; do echo "${kv}"; done
  } >"${STATUS_DIR}/${KEY}.status.tmp"
  mv "${STATUS_DIR}/${KEY}.status.tmp" "${STATUS_DIR}/${KEY}.status"
}

# fetch_and_verify <local_file> -- downloads SRC_URL into DOWNLOAD_DIR with the
# catalog's expected md5 (or the fetched sidecar for "verify" keys).
fetch_and_verify() {
  local local_file="$1"
  local url="${SRC_URL}" expected="${SRC_MD5}"
  TARGET_FILE="${local_file}"
  local dest="${DOWNLOAD_DIR}/${local_file}"

  if [[ "${expected}" == "verify" ]]; then
    if [[ -z "${SRC_SIDECAR}" ]]; then
      echo "ERROR: 'verify' mode needs a sidecar URL." >&2
      write_status FAIL "ERROR=no-sidecar-url"
      exit 1
    fi
    log "fetching md5 sidecar ${SRC_SIDECAR}"
    local sidecar_txt
    sidecar_txt="$(curl -L --fail --silent --show-error "${SRC_SIDECAR}")"
    EXPECTED_MD5="$(printf '%s' "${sidecar_txt}" | awk 'NR==1{print $1}')"
    if [[ ! "${EXPECTED_MD5}" =~ ^[0-9a-f]{32}$ ]]; then
      echo "ERROR: unparsable md5 sidecar: '${sidecar_txt}'" >&2
      write_status FAIL "ERROR=bad-md5-sidecar"
      exit 1
    fi
    expected="${EXPECTED_MD5}"
  else
    EXPECTED_MD5="${expected}"
  fi

  if [[ -f "${dest}" ]]; then
    local got
    got="$(md5_of "${dest}")"
    if [[ "${got}" == "${expected}" ]]; then
      ACTUAL_MD5="${got}"
      log "${local_file} already downloaded and verified (md5 ${got}); skipping."
      write_status OK "SKIPPED=already-verified"
      return 0
    fi
    log "${local_file} exists but md5 differs (expected ${expected}, got ${got}); resuming/re-downloading."
  fi

  local curl_args=(-L -C - --retry 5 --retry-delay 10 --fail --silent --show-error)
  if [[ -n "${LIMIT_RATE}" ]]; then
    curl_args+=(--limit-rate "${LIMIT_RATE}")
  fi
  curl_args+=(-o "${dest}" "${url}")
  log "downloading ${url} -> ${local_file} (resumable; Ctrl-C/re-queue safe)"
  set +e
  curl "${curl_args[@]}"
  local RC=$?
  set -e
  if [[ ${RC} -ne 0 ]]; then
    echo "ERROR: curl failed (rc=${RC}) for ${url}" >&2
    write_status FAIL "ERROR=curl-rc-${RC}"
    exit "${RC}"
  fi

  ACTUAL_MD5="$(md5_of "${dest}")"
  if [[ "${ACTUAL_MD5}" != "${expected}" ]]; then
    echo "ERROR: md5 mismatch for ${local_file}: expected ${expected}, got ${ACTUAL_MD5}" >&2
    write_status FAIL "ERROR=md5-mismatch"
    return 1
  fi
  log "md5 verified: ${ACTUAL_MD5}"
  write_status OK
}

# extract_tar_members <tar> -- selective extraction driven by SRC_TAR_PATTERNS
# ("egrep_pattern:canonical" pairs; empty pattern = the only .h5ad member),
# then delete the tar. Member matching greps the already-written
# _status/<key>.tar_listing (no repeated multi-GB tar -tzf decompressions).
extract_tar_members() {
  local tar="$1"
  log "listing ${tar} members -> _status/${KEY}.tar_listing"
  tar -tzf "${tar}" >"${STATUS_DIR}/${KEY}.tar_listing"
  [[ -n "${SRC_TAR_PATTERNS}" ]] || return 0

  local pair pattern canonical members n tmp extracted
  for pair in ${SRC_TAR_PATTERNS}; do
    pattern="${pair%%:*}"
    canonical="${pair#*:}"
    if [[ -n "${pattern}" ]]; then
      members="$(grep -iE '\.h5ad$' "${STATUS_DIR}/${KEY}.tar_listing" | grep -iE "${pattern}" || true)"
    else
      members="$(grep -iE '\.h5ad$' "${STATUS_DIR}/${KEY}.tar_listing" || true)"
    fi
    n="$(printf '%s\n' "${members}" | grep -c . || true)"
    if [[ "${n}" -ne 1 ]]; then
      echo "ERROR: ${n} h5ad member(s) match '${pattern}' (${canonical}) in ${tar}:" >&2
      printf '%s\n' "${members}" >&2
      write_status FAIL "ERROR=tar-member-ambiguous-${canonical}"
      return 1
    fi
    tmp="$(mktemp -d "${DOWNLOAD_DIR}/_tar_tmp.XXXXXX")"
    tar -xzf "${tar}" -C "${tmp}" "${members}"
    extracted="$(find "${tmp}" -type f -name '*.h5ad' | head -1)"
    mv "${extracted}" "${DOWNLOAD_DIR}/${canonical}"
    rm -rf "${tmp}"
    log "extracted + canonicalized ${members} -> ${canonical}"
  done
  rm -f "${tar}"
  log "${tar} deleted after selective extraction (excluded members left unextracted)."
}

# record_extracted <canonical files...> -- final OK status carrying the
# extracted files' md5s (MD5_<file>=) so the submitter tail can verify the NAS
# copies against worker-verified values.
record_extracted() {
  local kvs=("FILES=$(printf '%s,' "$@" | sed 's/,$//')")
  local nf
  for nf in "$@"; do
    kvs+=("MD5_${nf}=$(md5_of "${DOWNLOAD_DIR}/${nf}")")
  done
  TARGET_FILE=""
  write_status OK "${kvs[@]}"
}

# --- per-key dispatch (catalog-driven) -------------------------------------------
if [[ -n "${SRC_TAR_PATTERNS}" ]]; then
  # Already-extracted check: all canonical outputs present + tar gone -> done.
  nf_all=1
  nf=""
  for nf in $(printf '%s\n' "${SRC_TAR_PATTERNS}" | cut -d: -f2); do
    [[ -f "${DOWNLOAD_DIR}/${nf}" ]] || nf_all=0
  done
  tar="${DOWNLOAD_DIR}/${SRC_FILE}"
  if [[ ${nf_all} -eq 1 && ! -f "${tar}" ]]; then
    log "extracted files already present and ${SRC_FILE} gone; nothing to do."
    record_extracted $(printf '%s\n' "${SRC_TAR_PATTERNS}" | cut -d: -f2)
    exit 0
  fi
  fetch_and_verify "${SRC_FILE}"
  extract_tar_members "${tar}"
  record_extracted $(printf '%s\n' "${SRC_TAR_PATTERNS}" | cut -d: -f2)
else
  fetch_and_verify "${SRC_FILE}"
fi

log "done."
exit 0