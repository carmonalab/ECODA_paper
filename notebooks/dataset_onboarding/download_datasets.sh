#!/bin/bash
# ==============================================================================
# download_datasets.sh -- download the 9 Phase-5 datasets to the NAS
# ==============================================================================
# Downloads the PILOT-GM-VAE study (Joodaki et al. 2025, BIB bbaf547,
# PMID 41097818) author-provided files into the NAS folder
#   JooM_2025_41097818/output/
# i.e. the datasets.json `folder_name` target used by 1_stage_data.sh:
#   ${NAS_SC_DIR}/JooM_2025_41097818/output/<file>
#
# RULES (per the Phase-5 plan, decided with the user):
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

NAS_JOO_OUTPUT_DIR="${NAS_JOO_OUTPUT_DIR:-/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets/JooM_2025_41097818/output}"
LOG_FILE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/download_log.md"
ZENODO_PART1="https://zenodo.org/api/records/8370081/files"
ZENODO_PART2="https://zenodo.org/api/records/7957118/files"
ZENODO_PART3="https://zenodo.org/api/records/14615923/files"
CELLXGENE="https://datasets.cellxgene.cziscience.com"

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
  # fetch_and_verify <key> <url> <local_file> <expected_md5|verify>
  local key="$1" url="$2" local_file="$3" expected_md5="$4"
  local dest="${NAS_JOO_OUTPUT_DIR}/${local_file}"

  if [ -f "${dest}" ]; then
    local got
    got="$(md5_of "${dest}")"
    if [ "${expected_md5}" = "verify" ] || [ "${got}" = "${expected_md5}" ]; then
      echo "  [skip] ${local_file} already downloaded and verified (md5 ${got})."
      log_entry "OK (cached)" "${key}" "${local_file}" "$(du -h "${dest}" | cut -f1)" "${expected_md5}" "${got}" ""
      return 0
    fi
    echo "  [warn] ${local_file} exists but md5 differs (expected ${expected_md5}, got ${got}); re-downloading."
  fi

  echo "  [downloading] ${url}"
  echo "  -> ${local_file}  (resumable; Ctrl-C safe, re-run to continue)"
  curl -L -C - --retry 5 --retry-delay 10 --fail \
    -o "${dest}" "${url}"

  out="$(ls -l "${dest}" 2>/dev/null || echo '')"
  [ -n "${out}" ] && echo "  [done] size: $(du -h "${dest}" | cut -f1)"

  local got="$(md5_of "${dest}")"
  if [ "${expected_md5}" != "verify" ] && [ "${got}" != "${expected_md5}" ]; then
    echo "  [FAIL] md5 mismatch for ${local_file}: expected ${expected_md5}, got ${got}" >&2
    log_entry "FAIL (md5)" "${key}" "${local_file}" "$(du -h "${dest}" | cut -f1)" "${expected_md5}" "${got}" ""
    return 1
  fi
  echo "  [OK] md5 verified: ${got}"
  log_entry "OK" "${key}" "${local_file}" "$(du -h "${dest}" | cut -f1)" "${expected_md5}" "${got}" ""
}

# ---------------------------------------------------------------------------
# Datasets.tar.gz selective extraction: pick Covid-19 (Ren 2021) + Lupus
# (Perez 2022) h5ad members by filename pattern, canonicalize their names.
# ---------------------------------------------------------------------------
extract_covid_lupus() {
  local tar="${NAS_JOO_OUTPUT_DIR}/Datasets.tar.gz"
  echo "  [tar] listing members of Datasets.tar.gz (printing to download_log.md)..."
  tar -tzf "${tar}" | tee -a "${LOG_FILE}" > /dev/null

  local tmp
  extract_single_member() {
    # extract_single_member <pattern_egrep> <canonical_name> <context_label>
    local pattern="$1" canonical="$2" context="$3"
    local members
    members="$(tar -tzf "${tar}" | grep -iE '\.h5ad$' | grep -iE "${pattern}" || true)"
    local n
    n="$(printf '%s\n' "${members}" | grep -c . || true)"
    if [ "${n}" -eq 0 ]; then
      echo "  [FAIL] no h5ad member matching '${pattern}' (${context}) in Datasets.tar.gz." >&2
      echo "  All h5ad members in the tar:" >&2
      tar -tzf "${tar}" | grep -iE '\.h5ad$' | sed 's/^/    /' >&2
      exit 1
    fi
    if [ "${n}" -gt 1 ]; then
      echo "  [FAIL] ${n} members match '${pattern}' (${context}); canonical name ambiguous:" >&2
      printf '%s\n' "${members}" | sed 's/^/    /' >&2
      exit 1
    fi
    tmp="$(mktemp -d /tmp/onboard_tar.XXXXXX)"
    echo "  [tar] extracting '${members}' (${context})..."
    tar -xzf "${tar}" -C "${tmp}" "${members}"
    local extracted
    extracted="$(find "${tmp}" -type f -name '*.h5ad' | head -1)"
    mv "${extracted}" "${NAS_JOO_OUTPUT_DIR}/${canonical}"
    echo "  [OK] extracted + canonicalized: $(basename "${extracted}") -> ${canonical} ($(du -h "${NAS_JOO_OUTPUT_DIR}/${canonical}" | cut -f1))"
    rm -rf "${tmp}"
  }

  extract_single_member 'covid|ren_|ren-|GSE158055|coron' 'Covid19_Ren2021.h5ad' 'Covid-19'
  extract_single_member 'lupus|perez|GSE174188|SLE' 'Lupus_Perez2022.h5ad' 'Lupus'
  echo "  [tar] deleting Datasets.tar.gz (disk space)..."
  rm -f "${tar}"
  echo "  [tar] done. MI(1)/PDAC/Kidney cancer/follicular lymphoma members left unextracted (excluded datasets)."
}

extract_lung() {
  local tar="${NAS_JOO_OUTPUT_DIR}/lungatlas.h5ad.tar.gz"
  echo "  [tar] listing members of lungatlas.h5ad.tar.gz (printing to download_log.md)..."
  tar -tzf "${tar}" | tee -a "${LOG_FILE}" > /dev/null
  local tmp
  tmp="$(mktemp -d /tmp/onboard_lung.XXXXXX)"
  echo "  [tar] extracting lungatlas.h5ad..."
  tar -xzf "${tar}" -C "${tmp}"
  local found
  found="$(find "${tmp}" -type f -name '*.h5ad' | head -20 || true)"
  local n
  n="$(printf '%s\n' "${found}" | grep -c . || true)"
  if [ "${n}" -eq 0 ]; then
    echo "  [FAIL] no .h5ad member in lungatlas.h5ad.tar.gz" >&2
    tar -tzf "${tar}" | sed 's/^/    /' >&2
    exit 1
  fi
  if [ "${n}" -gt 1 ]; then
    echo "  [FAIL] multiple .h5ad members in lungatlas.h5ad.tar.gz:" >&2
    printf '%s\n' "${found}" | sed 's/^/    /' >&2
    exit 1
  fi
  mv "${found}" "${NAS_JOO_OUTPUT_DIR}/lungatlas.h5ad"
  echo "  [OK] extracted: $(basename "${found}") -> lungatlas.h5ad ($(du -h "${NAS_JOO_OUTPUT_DIR}/lungatlas.h5ad" | cut -f1))"
  rm -rf "${tmp}"
  echo "  [tar] deleting lungatlas.h5ad.tar.gz (disk space)..."
  rm -f "${tar}"
}

# ---------------------------------------------------------------------------
run_entry() {
  case "$1" in
    alzheimer)
      fetch_and_verify alzheimer \
        "${CELLXGENE}/c2b49431-9288-4d94-8ca5-f6723b72217e.h5ad" \
        SEAAD_Alzheimer.h5ad verify
      ;;
    breast)
      fetch_and_verify breast \
        "${ZENODO_PART3}/BreastCncr_processed.h5ad/content" \
        BreastCncr_processed.h5ad 8b28a349c2c3638ddbfb3946a32d12ba
      ;;
    covid_lupus_tar)
      fetch_and_verify covid_lupus_tar \
        "${ZENODO_PART1}/Datasets.tar.gz/content" \
        Datasets.tar.gz d105b52dbba38ac49c2ffe8b3cf34e24
      extract_covid_lupus
      ;;
    diabetes)
      fetch_and_verify diabetes \
        "${ZENODO_PART1}/diabetes.h5ad/content" \
        diabetes.h5ad 38189a381bad630fa39ce2d7ad3a0855
      ;;
    kidney)
      fetch_and_verify kidney \
        "${ZENODO_PART3}/Kidney_KPMP.h5ad/content" \
        Kidney_KPMP.h5ad 36ceb02ba23c559f80625ec7bef6884f
      ;;
    lung_tar)
      fetch_and_verify lung_tar \
        "${ZENODO_PART2}/lungatlas.h5ad.tar.gz/content" \
        lungatlas.h5ad.tar.gz 0d0c97924f1b7a405b6ec3b55da02882
      extract_lung
      ;;
    myocardial)
      fetch_and_verify myocardial \
        "${ZENODO_PART3}/Myocardial_Infarc_2.h5ad/content" \
        Myocardial_Infarc_2.h5ad 7431ae99250c99f11bf63e3034798af4
      ;;
    parkinson)
      fetch_and_verify parkinson \
        "${CELLXGENE}/0270e5e5-ce1d-4165-828e-699210189a92.h5ad" \
        Parkinson.h5ad verify
      ;;
    *)
      echo "ERROR: unknown --only key '$1'." >&2
      exit 1
      ;;
  esac
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