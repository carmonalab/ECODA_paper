#!/bin/bash
# ==============================================================================
# download_sources.sh -- single source of truth for the Phase-5 download catalog
# ==============================================================================
# Shared by the HPC download worker (run_download_worker.sh) and the Mac->NAS
# fallback (download_datasets.sh). Keep the catalog HERE only -- never hardcode
# URLs/checksums/tar patterns in the two scripts.
#
# API: source_download_key <key> sets, for the 8 keys
# (alzheimer breast covid_lupus_tar diabetes kidney lung_tar myocardial
# parkinson):
#   SRC_URL          full download URL
#   SRC_FILE         canonical filename in the download dir / on the NAS
#   SRC_MD5          expected md5, or the literal token "verify" for CellxGene
#                    files whose checksum must be fetched from the
#                    `<url>.md5` sidecar (SRC_SIDECAR) -- "verify" ALWAYS
#                    means "fetch the sidecar and enforce the checksum".
#   SRC_TAR_PATTERNS for tar keys: space-separated "egrep_pattern:canonical"
#                    pairs for selective extraction (empty pattern = the only
#                    .h5ad member); empty for direct-file keys.
# Returns 1 for unknown keys (sets nothing).
#
# Canonical filenames (identical on scratch and NAS, per the onboarding plan):
#   SEAAD_Alzheimer.h5ad  BreastCncr_processed.h5ad  Covid19_Ren2021.h5ad
#   diabetes.h5ad  Kidney_KPMP.h5ad  lungatlas.h5ad  Myocardial_Infarc_2.h5ad
#   Parkinson.h5ad
# Zenodo records: 8370081 (Part 1, latest of 7435911 -> 7956950 -> 8370081),
# 7957118 (Part 2), 14615923 (Part 3). CellxGene: SEA-AD dataset_version_id
# c2b49431-9288-4d94-8ca5-f6723b72217e, Parkinson 0270e5e5-ce1d-4165-828e-699210189a92.
# ==============================================================================

CELLXGENE_URL="https://datasets.cellxgene.cziscience.com"
ZENODO_PART1="https://zenodo.org/api/records/8370081/files"
ZENODO_PART2="https://zenodo.org/api/records/7957118/files"
ZENODO_PART3="https://zenodo.org/api/records/14615923/files"

source_download_key() {
  local key="$1"
  SRC_URL=""
  SRC_FILE=""
  SRC_MD5=""
  SRC_SIDECAR=""
  SRC_TAR_PATTERNS=""
  case "${key}" in
    alzheimer)
      SRC_URL="${CELLXGENE_URL}/c2b49431-9288-4d94-8ca5-f6723b72217e.h5ad"
      SRC_FILE="SEAAD_Alzheimer.h5ad"
      SRC_MD5="verify"
      ;;
    breast)
      SRC_URL="${ZENODO_PART3}/BreastCncr_processed.h5ad/content"
      SRC_FILE="BreastCncr_processed.h5ad"
      SRC_MD5="8b28a349c2c3638ddbfb3946a32d12ba"
      ;;
    covid_lupus_tar)
      SRC_URL="${ZENODO_PART1}/Datasets.tar.gz/content"
      SRC_FILE="Datasets.tar.gz"
      SRC_MD5="d105b52dbba38ac49c2ffe8b3cf34e24"
      SRC_TAR_PATTERNS="covid|ren_|ren-|GSE158055|coron:Covid19_Ren2021.h5ad lupus|perez|GSE174188|SLE:Lupus_Perez2022.h5ad"
      ;;
    diabetes)
      SRC_URL="${ZENODO_PART1}/diabetes.h5ad/content"
      SRC_FILE="diabetes.h5ad"
      SRC_MD5="38189a381bad630fa39ce2d7ad3a0855"
      ;;
    kidney)
      SRC_URL="${ZENODO_PART3}/Kidney_KPMP.h5ad/content"
      SRC_FILE="Kidney_KPMP.h5ad"
      SRC_MD5="36ceb02ba23c559f80625ec7bef6884f"
      ;;
    lung_tar)
      SRC_URL="${ZENODO_PART2}/lungatlas.h5ad.tar.gz/content"
      SRC_FILE="lungatlas.h5ad.tar.gz"
      SRC_MD5="0d0c97924f1b7a405b6ec3b55da02882"
      SRC_TAR_PATTERNS=":lungatlas.h5ad"
      ;;
    myocardial)
      SRC_URL="${ZENODO_PART3}/Myocardial_Infarc_2.h5ad/content"
      SRC_FILE="Myocardial_Infarc_2.h5ad"
      SRC_MD5="7431ae99250c99f11bf63e3034798af4"
      ;;
    parkinson)
      SRC_URL="${CELLXGENE_URL}/0270e5e5-ce1d-4165-828e-699210189a92.h5ad"
      SRC_FILE="Parkinson.h5ad"
      SRC_MD5="verify"
      ;;
    *)
      return 1
      ;;
  esac
  if [[ "${SRC_MD5}" == "verify" ]]; then
    SRC_SIDECAR="${SRC_URL}.md5"
  fi
}
