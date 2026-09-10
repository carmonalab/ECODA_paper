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

MD5_CALLS_FILE="${TMP_DIR}/md5-calls"
printf '0\n' > "${MD5_CALLS_FILE}"
md5_calls=0
ecoda_md5_file() {
  local count
  count="$(cat "${MD5_CALLS_FILE}")"
  count=$((count + 1))
  printf '%s\n' "${count}" > "${MD5_CALLS_FILE}"
  command md5sum "$1" | cut -d' ' -f1
}
md5_call_count() {
  tr -d '[:space:]' < "${MD5_CALLS_FILE}"
}

RUN_A_ID="checksum-reuse-a"
RUN_B_ID="checksum-reuse-b"
ecoda_init_run stage5 "${RUN_A_ID}" >/dev/null
if ecoda_init_run stage5 "${RUN_A_ID}" >/dev/null 2>&1; then
  echo "duplicate run ID was accepted" >&2
  exit 1
fi
FOREIGN_RUN="${TMP_DIR}/foreign-run"
mkdir -p "${FOREIGN_RUN}"
ln -s "${FOREIGN_RUN}" "${HPC_SCRATCH_DIR}/_ecoda_runs/symlink-run"
if ecoda_open_run symlink-run >/dev/null 2>&1; then
  echo "symlinked run root was accepted" >&2
  exit 1
fi
REDIRECT_TARGET="${TMP_DIR}/redirect-target"
mkdir -p "${REDIRECT_TARGET}"
ln -s "${REDIRECT_TARGET}" "${HPC_SCRATCH_DIR}/redirect"
if _ecoda_output_add_path "${HPC_SCRATCH_DIR}/redirect/new.bin" >/dev/null 2>&1; then
  echo "symlinked output parent was accepted" >&2
  exit 1
fi
FUTURE_PATH="${HPC_SCRATCH_DIR}/future-dir/artifact.bin"
FUTURE_ROOT="$(realpath "${HPC_SCRATCH_DIR}")"
[[ "$(ecoda_canonical_path "${FUTURE_PATH}")" == "${FUTURE_ROOT}/future-dir/artifact.bin" ]]

ecoda_validate_checksum "${FILE}"
[[ "${ECODA_CHECKSUM_PATH}" == "${FILE}" ]]
[[ "${ECODA_CHECKSUM_MD5}" == "${digest}" ]]
md5_calls="$(md5_call_count)"
[[ "${md5_calls}" == "1" ]]

# Publication reuses the immediately preceding strict digest.  Validation
# calls after publication must continue to use the recorded digest/size only.
RECORD="$(ecoda_write_artifact_record "${FILE}" checksum-producer "${RUN_A_ID}")"
[[ -s "${RECORD}" ]]
_ecoda_artifact_is_nonwritable "${FILE}"
validated_calls="${md5_calls}"
for repeat in 1 2 3; do
  ecoda_validate_checksum_record "${FILE}" "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}"
  ecoda_validate_artifact_record "${FILE}" checksum-producer "${RUN_A_ID}" >/dev/null
done
[[ "$(md5_call_count)" == "${validated_calls}" ]]

# A same-size mutation cannot use the published record.  Normal users are
# stopped by the immutable mode; if a privileged write succeeds, writable
# artifacts take the fresh checksum path and are still rejected.
cp "${FILE}" "${TMP_DIR}/artifact.original"
chmod u+w "${TMP_DIR}/artifact.original"
same_size_mutated=0
if printf 'X' | dd of="${FILE}" bs=1 count=1 conv=notrunc \
    >/dev/null 2>&1; then
  same_size_mutated=1
  if ecoda_validate_artifact_record "${FILE}" checksum-producer "${RUN_A_ID}" \
      >/dev/null 2>&1; then
    echo "same-size artifact mutation was accepted" >&2
    exit 1
  fi
  chmod u+w "${FILE}"
  cp "${TMP_DIR}/artifact.original" "${FILE}"
  chmod a-w "${FILE}"
fi
_ecoda_artifact_is_nonwritable "${FILE}"
validated_calls="$(md5_call_count)"

# A changed artifact size is rejected without trusting the old record.
cp "${FILE}.md5" "${TMP_DIR}/artifact.md5.original"
chmod u+w "${FILE}"
printf 'changed artifact size\n' >> "${FILE}"
if ecoda_validate_artifact_record "${FILE}" checksum-producer "${RUN_A_ID}" >/dev/null 2>&1; then
  echo "changed artifact size was accepted" >&2
  exit 1
fi
cp "${TMP_DIR}/artifact.original" "${FILE}"
chmod a-w "${FILE}"
validated_calls="$(md5_call_count)"
# A changed sidecar digest is rejected even when the artifact record is valid.
sed 's/^MD5=.*/MD5=00000000000000000000000000000000/' \
  "${FILE}.md5" > "${TMP_DIR}/artifact.bad.md5"
mv "${TMP_DIR}/artifact.bad.md5" "${FILE}.md5"
if ecoda_validate_artifact_record "${FILE}" checksum-producer "${RUN_A_ID}" >/dev/null 2>&1; then
  echo "changed sidecar digest was accepted" >&2
  exit 1
fi
cp "${TMP_DIR}/artifact.md5.original" "${FILE}.md5"

# Producer identity is part of the published record contract.
if ecoda_validate_artifact_record "${FILE}" wrong-producer "${RUN_A_ID}" >/dev/null 2>&1; then
  echo "wrong artifact producer was accepted" >&2
  exit 1
fi

ecoda_init_run stage5 "${RUN_B_ID}" >/dev/null
RUN_B_RECORD="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_B_ID}/manifests/artifacts/$(basename "${RECORD}")"
mkdir -p "$(dirname "${RUN_B_RECORD}")"
cp "${RECORD}" "${RUN_B_RECORD}"
if ecoda_validate_artifact_record "${FILE}" checksum-producer "${RUN_B_ID}" >/dev/null 2>&1; then
  echo "wrong artifact run ID was accepted" >&2
  exit 1
fi
[[ "$(md5_call_count)" == "${validated_calls}" ]]

# Preserve the existing strict sidecar PATH rejection.
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${FILE}" | tr -d '[:space:]')" \
  "${TMP_DIR}/foreign.bin" > "${FILE}.md5"
if ecoda_validate_checksum_record "${FILE}" "${digest}" "$(wc -c < "${FILE}" | tr -d '[:space:]')"; then
  echo "foreign checksum PATH was accepted" >&2
  exit 1
fi
cp "${TMP_DIR}/artifact.md5.original" "${FILE}.md5"

ecoda_invalidate_artifact "${FILE}"
[[ -w "${FILE}" ]]
printf 'replacement content\n' > "${FILE}"
ecoda_write_checksum "${FILE}"
chmod a-w "${FILE}"
written_calls="$(md5_call_count)"
ecoda_validate_checksum_record "${FILE}" "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}"
[[ "$(md5_call_count)" == "${written_calls}" ]]

REMOTE="${TMP_DIR}/remote.bin"
cp "${FILE}" "${REMOTE}"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${ECODA_CHECKSUM_MD5}" \
  "${ECODA_CHECKSUM_SIZE}" "${REMOTE}" > "${REMOTE}.md5"
printf '0\n' > "${MD5_CALLS_FILE}"
md5_calls=0
ecoda_compare_checksum_remote "${FILE}" "${REMOTE}" "${REMOTE}.md5" \
  "${ECODA_CHECKSUM_MD5}" "${ECODA_CHECKSUM_SIZE}"
md5_calls="$(md5_call_count)"
[[ "${md5_calls}" == "1" ]]

echo "checksum reuse: OK"
