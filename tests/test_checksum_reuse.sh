#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-checksum-reuse.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT

export HPC_SCRATCH_DIR="${TMP_DIR}/scratch"
FILE="${TMP_DIR}/artifact.bin"
printf 'checksum reuse fixture\n' > "${FILE}"
digest="$(md5sum "${FILE}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${FILE}" | tr -d '[:space:]')" "${FILE}" > "${FILE}.md5"

source "${ROOT}/src/utils/bash/ecoda_run_common.sh"

md5_calls=0
ecoda_md5_file() {
  md5_calls=$((md5_calls + 1))
  command md5sum "$1" | cut -d' ' -f1
}

ecoda_validate_checksum "${FILE}"
[[ "${ECODA_CHECKSUM_PATH}" == "${FILE}" ]]
[[ "${ECODA_CHECKSUM_MD5}" == "${digest}" ]]
validated_calls="${md5_calls}"
ecoda_validate_checksum_record "${FILE}" "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}"
[[ "${md5_calls}" == "${validated_calls}" ]]

printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${FILE}" | tr -d '[:space:]')" "${TMP_DIR}/foreign.bin" > "${FILE}.md5"
if ecoda_validate_checksum_record "${FILE}" "${digest}" "$(wc -c < "${FILE}" | tr -d '[:space:]')"; then
  echo "foreign checksum PATH was accepted" >&2
  exit 1
fi

printf 'replacement content\n' > "${FILE}"
ecoda_write_checksum "${FILE}"
written_calls="${md5_calls}"
ecoda_validate_checksum_record "${FILE}" "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}"
[[ "${md5_calls}" == "${written_calls}" ]]

REMOTE="${TMP_DIR}/remote.bin"
cp "${FILE}" "${REMOTE}"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${ECODA_CHECKSUM_MD5}" \
  "${ECODA_CHECKSUM_SIZE}" "${REMOTE}" > "${REMOTE}.md5"
md5_calls=0
ecoda_compare_checksum_remote "${FILE}" "${REMOTE}" "${REMOTE}.md5" \
  "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}"
[[ "${md5_calls}" == "1" ]]

echo "checksum reuse: OK"
