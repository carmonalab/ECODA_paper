#!/bin/bash
# Stage and checksum the four annotation reference maps exactly once.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

REF_MAP_NAMES=(
  sketched_CD8T_human_ref_v1.rds
  sketched_CD4T_human_ref_v2.rds
  sketched_DC_human_ref_v2.rds
  sketched_MoMac_human_v1.rds
)
ref_map_id() {
  case "$1" in
    sketched_CD8T_human_ref_v1.rds) printf '47714158' ;;
    sketched_CD4T_human_ref_v2.rds) printf '47714155' ;;
    sketched_DC_human_ref_v2.rds) printf '47714161' ;;
    sketched_MoMac_human_v1.rds) printf '47714164' ;;
    *) return 1 ;;
  esac
}
ref_map_md5() {
  case "$1" in
    sketched_CD8T_human_ref_v1.rds) printf 'be86058ddafdd0154faf0485286b86e7' ;;
    sketched_CD4T_human_ref_v2.rds) printf '5540a0ee287e291528c96d476794b194' ;;
    sketched_DC_human_ref_v2.rds) printf '033d491ba7ca9bbf0badcae828e55b2c' ;;
    sketched_MoMac_human_v1.rds) printf '3043cd9058a8746d972c7be195b18e36' ;;
    *) return 1 ;;
  esac
}
md5_of() {
  if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi
}
mkdir -p "${HOME_REF_DIR}"
for name in "${REF_MAP_NAMES[@]}"; do
  destination="${HOME_REF_DIR}/${name}"
  if [[ -s "${destination}" && "$(md5_of "${destination}")" == "$(ref_map_md5 "${name}")" ]]; then
    continue
  fi
  rm -f "${destination}"
  if [[ -r "${NAS_REF_DIR}/${name}" ]]; then
    tmp="${destination}.tmp.$$"
    cp "${NAS_REF_DIR}/${name}" "${tmp}"
    if [[ "$(md5_of "${tmp}")" == "$(ref_map_md5 "${name}")" ]]; then mv -f "${tmp}" "${destination}"; else rm -f "${tmp}"; fi
  fi
  if [[ ! -s "${destination}" ]]; then
    url="https://ndownloader.figshare.com/files/$(ref_map_id "${name}")"
    tmp="${destination}.tmp.$$"
    curl -f -L --retry 3 -o "${tmp}" "${url}"
    [[ "$(md5_of "${tmp}")" == "$(ref_map_md5 "${name}")" ]] || { rm -f "${tmp}"; echo "ERROR: MD5 mismatch for ${name}" >&2; exit 1; }
    mv -f "${tmp}" "${destination}"
  fi
done
for name in "${REF_MAP_NAMES[@]}"; do
  [[ -s "${HOME_REF_DIR}/${name}" ]] || { echo "ERROR: reference map missing: ${name}" >&2; exit 1; }
  [[ "$(md5_of "${HOME_REF_DIR}/${name}")" == "$(ref_map_md5 "${name}")" ]] || { echo "ERROR: reference map checksum failed: ${name}" >&2; exit 1; }
done
echo "Annotation reference maps staged and MD5-verified in ${HOME_REF_DIR}."
