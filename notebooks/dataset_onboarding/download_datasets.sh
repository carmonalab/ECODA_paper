#!/bin/bash
# ==============================================================================
# download_datasets.sh -- download the 9 Phase-5 datasets to the NAS
# ==============================================================================
# FALLBACK route ONLY (NAS-stable): downloads the PILOT-GM-VAE study (Joodaki
# et al. 2025, BIB bbaf547, PMID 41097818) author-provided files into the NAS
# folder
#   JooM_2025_41097818/output/
# i.e. the datasets.json `folder_name` target used by 1_stage_data.sh:
#   ${NAS_SC_DIR}/JooM_2025_41097818/output/<file>
#
# PRIMARY ROUTE (user decision 2026-08-17): the HPC downloader
# `download_datasets_hpc.sh` (compute-node array or login-node fallback). Use
# this Mac->NAS script only when the NAS mount is stable and the HPC route is
# unavailable.
#
# The download catalog (URLs, md5s, CellxGene sidecar handling, tar member
# patterns) lives in download_sources.sh -- the single source of truth shared
# with the HPC worker. "verify" md5s fetch the `.h5ad.md5` sidecar and ENFORCE
# the checksum (same semantics as the HPC worker).
#
# RULES:
#   * Sequential ONLY -- one download at a time, each md5-verified before the
#     next (bandwidth + NAS SMB). Do NOT run two copies of this script.
#   * `curl -L -C -` resumable; re-running the script continues where it left
#     off (already-verified files are skipped; partial files resume).
#   * Selective tar extraction: only the needed h5ads are pulled out of
#     Datasets.tar.gz / lungatlas.h5ad.tar.gz, then the tars are DELETED (the
#     Mac has only ~108 GB free).
#   * Excluded datasets are never extracted/registered: Kidney cancer
#     (kidney_cancer_processed.h5ad in record 14615923), PDAC, MI(1),
#     follicular lymphoma (follicular_lymphoma.h5ad in record 8370081).
#
# Usage (run from the repo root; NAS must be mounted):
#   ./notebooks/dataset_onboarding/download_datasets.sh            # all, sequentially
#   ./notebooks/dataset_onboarding/download_datasets.sh --only breast   # one entry
#   --only accepts: alzheimer|breast|covid_lupus_tar|diabetes|kidney|lung_tar|myocardial|parkinson
#
# Zenodo note: the paper's Data Availability cites records 7956950/7957118/
# 14615923, but record 8370081 is the latest version of the same concept
# record (7435911 -> 7956950 -> 8370081) and its files are byte-identical
# (same md5). We download from 8370081 (latest).
# ==============================================================================
set -euo pipefail

ONBOARD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ONBOARD_DIR}/download_sources.sh"

NAS_JOO_OUTPUT_DIR="${NAS_JOO_OUTPUT_DIR:-/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets/JooM_2025_41097818/output}"
LOG_FILE="${ONBOARD_DIR}/download_log.md"

mkdir -p "${NAS_JOO_OUTPUT_DIR}"

MD5_BIN=""
if command -v md5 >/dev/null 2>&1; then
  MD5_BIN="md5 -q"
elif command -v md5sum >/dev/null 2>&1; then
  MD5_BIN="md5sum"
else
  echo "ERROR: no md5/md5sum available." >&2
  exit 1
fi

md5_of() {
  if [ "$MD5_BIN" = "md5sum" ]; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

log_entry() {
  # log_entry <status> <key> <file> <size> <md5_expected> <md5_actual> <extra>
  local ts
  ts="$(date '+%Y-%m-%d %H:%M')"
  printf '\n### %s\n\n- **Status:** %s\n- **Key:** %s\n- **File:** %s\n- **Size:** %s\n- **md5 expected:** %s\n- **md5 actual:** %s\n- **Extra:** %s\n' \
    "${ts}" "$1" "$2" "$3" "$4" "$5" "$6" "$7" >> "${LOG_FILE}"
}

fetch_and_verify() {
  # fetch_and_verify <key> -- downloads SRC_URL (from download_sources.sh) into
  # the NAS output dir, enforcing SRC_MD5 ("verify" = fetch the .h5ad.md5
  # sidecar and enforce it).
  local key="$1" url="${SRC_URL}" expected="${SRC_MD5}" local_file="${SRC_FILE}"
  local dest="${NAS_JOO_OUTPUT_DIR}/${local_file}"

  if [[ "${expected}" == "verify" ]]; then
    echo "  [sidecar] fetching ${SRC_SIDECAR}"
    local sidecar_txt got_side
    sidecar_txt="$(curl -L --fail --silent --show-error "${SRC_SIDECAR}")"
    expected="$(printf '%s' "${sidecar_txt}" | awk 'NR==1{print $1}')"
    if [[ ! "${expected}" =~ ^[0-9a-f]{32}$ ]]; then
      echo "  [FAIL] unparsable md5 sidecar: '${sidecar_txt}'" >&2
      log_entry "FAIL (sidecar)" "${key}" "${local_file}" "" "${expected}" "" ""
      return 1
    fi
  fi

  if [ -f "${dest}" ]; then
    local got
    got="$(md5_of "${dest}")"
    if [ "${got}" = "${expected}" ]; then
      echo "  [skip] ${local_file} already downloaded and verified (md5 ${got})."
      log_entry "OK (cached)" "${key}" "${local_file}" "$(du -h "${dest}" | cut -f1)" "${expected}" "${got}" ""
      return 0
    fi
    echo "  [warn] ${local_file} exists but md5 differs (expected ${expected}, got ${got}); re-downloading."
  fi

  echo "  [downloading] ${url}"
  echo "  -> ${local_file}  (resumable; Ctrl-C safe, re-run to continue)"
  curl -L -C - --retry 5 --retry-delay 10 --fail \
    -o "${dest}" "${url}"

  out="$(ls -l "${dest}" 2>/dev/null || echo '')"
  [ -n "${out}" ] && echo "  [done] size: $(du -h "${dest}" | cut -f1)"

  local got="$(md5_of "${dest}")"
  if [ "${got}" != "${expected}" ]; then
    echo "  [FAIL] md5 mismatch for ${local_file}: expected ${expected}, got ${got}" >&2
    log_entry "FAIL (md5)" "${key}" "${local_file}" "$(du -h "${dest}" | cut -f1)" "${expected}" "${got}" ""
    return 1
  fi
  echo "  [OK] md5 verified: ${got}"
  log_entry "OK" "${key}" "${local_file}" "$(du -h "${dest}" | cut -f1)" "${expected}" "${got}" ""
}

# ---------------------------------------------------------------------------
# Selective tar extraction driven by SRC_TAR_PATTERNS ("egrep_pattern:canonical"
# pairs; empty pattern = the only .h5ad member). Only the needed h5ads are
# extracted; excluded members stay in the archive, then the tar is DELETED.
# ---------------------------------------------------------------------------
extract_tar_members() {
  local key="$1" tar="${NAS_JOO_OUTPUT_DIR}/${SRC_FILE}"
  echo "  [tar] listing members of ${SRC_FILE} (printing to download_log.md)..."
  tar -tzf "${tar}" | tee -a "${LOG_FILE}" > /dev/null
  [[ -n "${SRC_TAR_PATTERNS}" ]] || return 0

  local pair pattern canonical members n tmp extracted
  for pair in ${SRC_TAR_PATTERNS}; do
    pattern="${pair%%:*}"
    canonical="${pair#*:}"
    if [[ -n "${pattern}" ]]; then
      members="$(tar -tzf "${tar}" | grep -iE '\.h5ad$' | grep -iE "${pattern}" || true)"
    else
      members="$(tar -tzf "${tar}" | grep -iE '\.h5ad$' || true)"
    fi
    n="$(printf '%s\n' "${members}" | grep -c . || true)"
    if [[ "${n}" -ne 1 ]]; then
      echo "  [FAIL] ${n} h5ad member(s) matching '${pattern}' (${canonical}) in ${SRC_FILE}:" >&2
      printf '%s\n' "${members}" | sed 's/^/    /' >&2
      exit 1
    fi
    tmp="$(mktemp -d /tmp/onboard_tar.XXXXXX)"
    echo "  [tar] extracting '${members}' -> ${canonical}..."
    tar -xzf "${tar}" -C "${tmp}" "${members}"
    extracted="$(find "${tmp}" -type f -name '*.h5ad' | head -1)"
    mv "${extracted}" "${NAS_JOO_OUTPUT_DIR}/${canonical}"
    echo "  [OK] extracted + canonicalized: $(basename "${extracted}") -> ${canonical} ($(du -h "${NAS_JOO_OUTPUT_DIR}/${canonical}" | cut -f1))"
    rm -rf "${tmp}"
  done
  echo "  [tar] deleting ${SRC_FILE} (disk space)..."
  rm -f "${tar}"
  echo "  [tar] done. Excluded members (MI(1)/PDAC/Kidney cancer/follicular lymphoma) left unextracted."
}

# ---------------------------------------------------------------------------
run_entry() {
  local key="$1"
  if ! source_download_key "${key}"; then
    echo "ERROR: unknown --only key '${key}'." >&2
    exit 1
  fi
  if [[ -n "${SRC_TAR_PATTERNS}" ]]; then
    fetch_and_verify "${key}"
    extract_tar_members "${key}"
  else
    fetch_and_verify "${key}"
  fi
}

# ---------------------------------------------------------------------------
ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only)
      ONLY="${2:-}"
      shift 2
      ;;
    *)
      echo "ERROR: unknown argument $1 (only --only <key> is accepted)." >&2
      exit 1
      ;;
  esac
done

if [ ! -d "$(dirname "${NAS_JOO_OUTPUT_DIR}")" ]; then
  echo "ERROR: NAS not mounted? $(dirname "${NAS_JOO_OUTPUT_DIR}") is missing." >&2
  exit 1
fi

if [ ! -f "${LOG_FILE}" ]; then
  cat > "${LOG_FILE}" <<'EOF'
# Dataset download log (Phase 5 onboarding -- JooM_2025_41097818)

Source records (see .kilo/plans/1786899069337-onboard-new-datasets-phase5.md):
- Zenodo 8370081 (Part 1, latest version of 7435911 -> 7956950 -> 8370081):
  Datasets.tar.gz, diabetes.h5ad
- Zenodo 7957118 (Part 2): lungatlas.h5ad.tar.gz
- Zenodo 14615923 (Part 3): BreastCncr_processed.h5ad, Kidney_KPMP.h5ad,
  Myocardial_Infarc_2.h5ad (also holds kidney_cancer_processed.h5ad --
  EXCLUDED dataset, not used)
- CellxGene: SEA-AD collection 1ca90a2d-2943-483d-b678-b809bf464c30
  (dataset_version_id c2b49431-9288-4d94-8ca5-f6723b72217e, 1,395,601 cells
  -- exact match with the paper) and Parkinson collection
  d5d0df8f-4eee-49d8-a221-a288f50a1590 (dataset_version_id
  0270e5e5-ce1d-4165-828e-699210189a92, 2,096,155 cells -- exact match).

Recorded 2026-08-17 (planning): Zenodo md5 values confirmed via the API.
EOF
fi

ALL_KEYS="alzheimer breast covid_lupus_tar diabetes kidney lung_tar myocardial parkinson"
if [ -n "${ONLY}" ]; then
  if ! printf '%s\n' ${ALL_KEYS} | grep -qx "${ONLY}"; then
    echo "ERROR: --only '${ONLY}' invalid. Valid: ${ALL_KEYS}" >&2
    exit 1
  fi
  echo "=== Downloading: ${ONLY} (sequential, md5-verified) ==="
  run_entry "${ONLY}"
else
  echo "=== Downloading all datasets (SEQUENTIAL, one at a time) ==="
  for k in ${ALL_KEYS}; do
    echo ""
    echo ">>> ${k}"
    run_entry "${k}"
  done
fi

echo ""
echo "All requested downloads complete. Verify the destination:"
echo "  ls -lh ${NAS_JOO_OUTPUT_DIR}"
echo "Log: ${LOG_FILE}"