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
#   "verify" md5s (CellxGene) are SIZE-verified via HEAD content-length (no
#   `.h5ad.md5` sidecar exists -- 403 -- and the S3 ETag is a multipart
#   digest); the computed md5 is recorded as informational.
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
# catalog's expected md5; for "verify" keys (CellxGene) no md5 sidecar exists
# (the `.h5ad.md5` URL returns 403 and the S3 ETag is a multipart digest, so
# neither is usable) -- the final file SIZE is checked against a HEAD
# content-length instead, and the computed md5 is recorded as informational.
fetch_and_verify() {
  local local_file="$1"
  local url="${SRC_URL}" expected="${SRC_MD5}"
  TARGET_FILE="${local_file}"
  local dest="${DOWNLOAD_DIR}/${local_file}"
  EXPECTED_MD5=""
  EXPECTED_SIZE=""

  if [[ "${expected}" == "verify" ]]; then
    log "resolving expected size via HEAD ${url}"
    local hdr
    hdr="$(curl -sIL --fail --silent --show-error "${url}" | tr -d '\r' | awk -F': ' 'tolower($1)=="content-length"{len=$2} END{print len}')"
    EXPECTED_SIZE="$(printf '%s' "${hdr}" | tr -d ' ')"
    if [[ ! "${EXPECTED_SIZE}" =~ ^[0-9]+$ ]]; then
      echo "ERROR: could not resolve content-length for ${url} ('${hdr}')" >&2
      write_status FAIL "ERROR=no-content-length"
      exit 1
    fi
    expected=""
  fi

  if [[ -f "${dest}" ]]; then
    local got
    if [[ -n "${EXPECTED_SIZE}" ]]; then
      got="$(stat -c %s "${dest}")"
      if [[ "${got}" == "${EXPECTED_SIZE}" ]]; then
        ACTUAL_MD5="$(md5_of "${dest}")"
        log "${local_file} already downloaded and size-verified (${got} bytes); skipping."
        write_status OK "SKIPPED=already-verified" "VERIFY=size-${EXPECTED_SIZE}" "MD5_RECORDED=${ACTUAL_MD5}"
        return 0
      fi
      log "${local_file} exists but size differs (expected ${EXPECTED_SIZE}, got ${got}); resuming/re-downloading."
    else
      got="$(md5_of "${dest}")"
      if [[ "${got}" == "${expected}" ]]; then
        ACTUAL_MD5="${got}"
        log "${local_file} already downloaded and verified (md5 ${got}); skipping."
        write_status OK "SKIPPED=already-verified"
        return 0
      fi
      log "${local_file} exists but md5 differs (expected ${expected}, got ${got}); resuming/re-downloading."
    fi
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

  if [[ -n "${EXPECTED_SIZE}" ]]; then
    local got_size
    got_size="$(stat -c %s "${dest}")"
    if [[ "${got_size}" != "${EXPECTED_SIZE}" ]]; then
      echo "ERROR: size mismatch for ${local_file}: expected ${EXPECTED_SIZE}, got ${got_size}" >&2
      write_status FAIL "ERROR=size-mismatch"
      return 1
    fi
    ACTUAL_MD5="$(md5_of "${dest}")"
    log "size verified: ${got_size} bytes; md5 (informational): ${ACTUAL_MD5}"
    write_status OK "VERIFY=size-${EXPECTED_SIZE}" "MD5_RECORDED=${ACTUAL_MD5}"
  else
    ACTUAL_MD5="$(md5_of "${dest}")"
    if [[ "${ACTUAL_MD5}" != "${expected}" ]]; then
      echo "ERROR: md5 mismatch for ${local_file}: expected ${expected}, got ${ACTUAL_MD5}" >&2
      write_status FAIL "ERROR=md5-mismatch"
      return 1
    fi
    log "md5 verified: ${ACTUAL_MD5}"
    write_status OK
  fi
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

# tar_canonicals -- canonical filenames for the tar key, one per line.
# SRC_TAR_PATTERNS is a space-separated list of "pattern:canonical" pairs and
# CAN hold multiple pairs (e.g. covid_lupus_tar). A single `cut -d: -f2` over
# the whole string mis-fields multi-pair keys (the second pair's PATTERN lands
# in the middle of field 2, corrupting the arg list for record_extracted ->
# md5sum on a bogus "pattern" filename -> task FAIL exit 1). Derive each
# canonical per-pair with ${pair#*:} instead.
tar_canonicals() {
  local pair
  for pair in ${SRC_TAR_PATTERNS}; do
    printf '%s\n' "${pair#*:}"
  done
}

# --- per-key dispatch (catalog-driven) ---
if [[ -n "${SRC_TAR_PATTERNS}" ]]; then
  # Already-extracted check: all canonical outputs present + tar gone -> done.
  nf_all=1
  nf=""
  for nf in $(tar_canonicals); do
    [[ -f "${DOWNLOAD_DIR}/${nf}" ]] || nf_all=0
  done
  tar="${DOWNLOAD_DIR}/${SRC_FILE}"
  if [[ ${nf_all} -eq 1 && ! -f "${tar}" ]]; then
    log "extracted files already present and ${SRC_FILE} gone; nothing to do."
    record_extracted $(tar_canonicals)
    exit 0
  fi
  fetch_and_verify "${SRC_FILE}"
  extract_tar_members "${tar}"
  record_extracted $(tar_canonicals)
else
  fetch_and_verify "${SRC_FILE}"
fi

log "done."
exit 0