#!/bin/bash
# Simulated findmnt/capacity boundaries with real descriptor-relative copies.
# This is not live NAS proof. The independent Bamboo probe remains mandatory.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd -P)"
: "${PYTHON_BIN:?set PYTHON_BIN to the configured interpreter}"
real_python="${PYTHON_BIN}"
real_chmod="$(command -v chmod)"
work="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-nas-test.XXXXXX")"
work="$(cd "${work}" && pwd -P)"
trap 'chmod -R u+w "${work}"; rm -rf "${work}"' EXIT
export HPC_SCRATCH_DIR="${work}/scratch"
export ECODA_RUN_ID=fixture_promotion
export ECODA_RUN_ROOT="${work}/${ECODA_RUN_ID}"
export FIXTURE_ROOT="${work}"
export MOUNT_KIND=valid
mkdir -p "${work}/bin" "${work}/python" "${work}/nas,comma" "${work}/source" \
  "${work}/snapshot/tree/aux" "${work}/snapshot/identity" "${ECODA_RUN_ROOT}/manifests"
source "${repo}/src/utils/bash/ecoda_run_common.sh"
cat > "${work}/bin/findmnt" <<'SH'
#!/bin/bash
case "${MOUNT_KIND}" in
  absent) exit 1 ;;
  wrong) source='othernas:/wrong'; type=nfs ;;
  scratch) source='/dev/beegfs'; type=beegfs ;;
  autofs) source='systemd-1'; type=autofs ;;
  *) source='smednas515:/volume1/Shared'; type=nfs ;;
esac
target=''
while (( $# )); do
  if [[ "$1" == --target ]]; then target="$2"; break; fi
  shift
done
if [[ -n "${WRONG_MOUNT_TARGET:-}" &&
      ( "${target}" == "${WRONG_MOUNT_TARGET}" || "${target}" == "${WRONG_MOUNT_TARGET}/"* ) ]]; then
  source='/dev/nested'; type=ext4
fi
printf '{"filesystems":[{"source":"%s","target":"%s/nas,comma","fstype":"%s"}]}\n' \
  "$source" "$FIXTURE_ROOT" "$type"
SH
cat > "${work}/bin/rsync" <<'SH'
#!/bin/bash
printf 'copy\n' >> "${FIXTURE_ROOT}/rsync.calls"
exit 99
SH
cat > "${work}/bin/sbatch" <<'SH'
#!/bin/bash
printf 'FORBIDDEN\n' >> "${FIXTURE_ROOT}/sbatch.calls"
exit 99
SH
# Capacity, ineffective fchmod, and exact descriptor race timing are injected here.
cat > "${work}/python/sitecustomize.py" <<'PY'
import os
if os.environ.get('NO_CAPACITY') == '1':
    original = os.statvfs
    def statvfs(path):
        value = original(path)
        return os.statvfs_result((value.f_bsize, value.f_frsize, value.f_blocks,
                                 value.f_bfree, 0, value.f_files, value.f_ffree,
                                 value.f_favail, value.f_flag, value.f_namemax))
    os.statvfs = statvfs
original_fchmod = os.fchmod
def fchmod(fd, mode):
    attack = os.environ.get('SWAP_PATH')
    if attack and not os.path.exists(attack + '.saved'):
        current, pinned = os.stat(attack, follow_symlinks=False), os.fstat(fd)
        if (current.st_dev, current.st_ino) == (pinned.st_dev, pinned.st_ino):
            os.rename(attack, attack + '.saved')
            os.symlink(os.environ['SWAP_TARGET'], attack)
            with open(os.environ['FIXTURE_ROOT'] + '/swap.triggered', 'w') as stream:
                stream.write(attack)
    if os.environ.get('INEFFECTIVE_CHMOD') != '1':
        original_fchmod(fd, mode)
os.fchmod = fchmod
original_open = os.open
def descriptor_open(path, flags, mode=0o777, *, dir_fd=None):
    if path == 'created' and dir_fd is not None and os.environ.get('PROBE_SWAP_TARGET'):
        pinned = os.fstat(dir_fd)
        with os.scandir(os.environ['FIXTURE_ROOT'] + '/nas,comma') as entries:
            for entry in entries:
                if not entry.name.startswith('.ecoda-nas-probe-'):
                    continue
                current = entry.stat(follow_symlinks=False)
                if (current.st_dev, current.st_ino) == (pinned.st_dev, pinned.st_ino):
                    os.rename(entry.path, entry.path + '.saved')
                    os.symlink(os.environ['PROBE_SWAP_TARGET'], entry.path)
                    with open(os.environ['FIXTURE_ROOT'] + '/probe-swap.triggered', 'w') as marker:
                        marker.write(entry.path)
                    break
    if str(path).startswith('.ecoda-copy-') and dir_fd is not None:
        with open(os.environ['FIXTURE_ROOT'] + '/copy.calls', 'a') as calls:
            calls.write('copy\n')
        attack = os.environ.get('COPY_SWAP_PATH')
        if attack and not os.path.exists(attack + '.saved'):
            current, pinned = os.stat(attack, follow_symlinks=False), os.fstat(dir_fd)
            if (current.st_dev, current.st_ino) == (pinned.st_dev, pinned.st_ino):
                os.rename(attack, attack + '.saved')
                os.symlink(os.environ['COPY_SWAP_TARGET'], attack)
                with open(os.environ['FIXTURE_ROOT'] + '/copy-swap.triggered', 'w') as marker:
                    marker.write(attack)
        for key in ('MUTATE_SELECTION', 'MUTATE_SOURCE'):
            if os.environ.get(key):
                with open(os.environ[key], 'a') as stream:
                    stream.write('changed\n')
    return original_open(path, flags, mode, dir_fd=dir_fd)
os.open = descriptor_open
PY
"${real_chmod}" +x "${work}/bin/"*
export PATH="${work}/bin:${PATH}" PYTHONPATH="${work}/python${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHON_BIN="${real_python}"
printf 'dataset\tcorrected\n' > "${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'real fixture payload\n' > "${work}/source/input.h5ad"
printf 'MD5=original-producer\nSIZE=21\nPATH=/original/source/input.h5ad\n' > "${work}/source/input.h5ad.md5"
mkdir "${work}/source/nested"
printf 'nested fixture payload\n' > "${work}/source/nested/payload"
printf 'input.h5ad\ninput.h5ad.md5\nnested/payload\n' > "${work}/files.txt"
for name in input.h5ad input.h5ad.md5 nested/payload; do
  printf '%s\t%s\t%s\n' "${name}" "$(ecoda_md5_file "${work}/source/${name}")" \
    "$(wc -c < "${work}/source/${name}" | tr -d '[:space:]')"
done > "${ECODA_RUN_ROOT}/manifests/inventory.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/inventory.tsv"
printf 'archive bytes\n' > "${work}/snapshot/identity/source.tar"
printf 'complete\n' > "${work}/snapshot/COMPLETE"
for name in config_helper.R datasets.json pixi.toml pixi.lock; do
  printf 'immutable fixture %s\n' "${name}" > "${work}/snapshot/tree/${name}"
done
manifest="${work}/snapshot/identity/source.manifest"
printf 'FORMAT=1\nSOURCE_ROOT=%s/snapshot/tree\nSOURCE_COMMIT=1111111111111111111111111111111111111111\nSOURCE_ARCHIVE_PATH=%s/snapshot/identity/source.tar\nSOURCE_ARCHIVE_SHA256=%s\nCONFIG_HELPER_SHA256=%s\nDATASETS_SHA256=%s\nPIXI_TOML_SHA256=%s\nPIXI_LOCK_SHA256=%s\nAUX_ROOT=%s/snapshot/tree/aux\nSCGATE_DB_BRANCH=fixture\n' \
  "${work}" "${work}" "$(ecoda_sha256_file "${work}/snapshot/identity/source.tar")" \
  "$(ecoda_sha256_file "${work}/snapshot/tree/config_helper.R")" \
  "$(ecoda_sha256_file "${work}/snapshot/tree/datasets.json")" \
  "$(ecoda_sha256_file "${work}/snapshot/tree/pixi.toml")" \
  "$(ecoda_sha256_file "${work}/snapshot/tree/pixi.lock")" "${work}" > "${manifest}"
"${real_chmod}" 444 "${manifest}"
utility="${repo}/src/utils/bash/ecoda_nas_backup.sh"
nas="${work}/nas,comma"
destination="${nas}/_ecoda_backups/corrected_final/${ECODA_RUN_ID}/data"
receipt="${ECODA_RUN_ROOT}/manifests/backup.json"
copy() {
  bash "${utility}" copy --source-root "${work}/source" --nas-root "${nas}" \
    --destination-root "${destination}" --files-from "${work}/files.txt" \
    --inventory "${ECODA_RUN_ROOT}/manifests/inventory.tsv" --receipt "${receipt}" \
    --selection-file "${ECODA_RUN_ROOT}/manifests/selection.tsv" --source-manifest "${manifest}"
}
verify() { bash "${utility}" verify --receipt "${receipt}"; }
: > "${work}/rsync.calls"
: > "${work}/copy.calls"
reject() { if "$@"; then echo "unexpected success: $*" >&2; exit 1; fi; }
for MOUNT_KIND in absent wrong scratch autofs; do
  export MOUNT_KIND
  reject bash "${utility}" preflight --nas-root "${nas}" --report "${work}/preflight.json"
  [[ ! -e "${work}/preflight.json" ]]
done
export MOUNT_KIND=valid
bash "${utility}" preflight --nas-root "${nas}" --report "${work}/preflight.json"
"${PYTHON_BIN}" - "${work}/preflight.json" "${nas}" <<'PY'
import json, sys
value = json.load(open(sys.argv[1]))
assert value['nas_root'] == sys.argv[2]
assert value['mount_source'] == 'smednas515:/volume1/Shared'
assert value['available_bytes'] > 0
PY
# Swap the probe after its directory-mode/chain checks, immediately before create.
# Both preflight and copy must reject without touching scratch/home or publishing.
for victim in "${work}/scratch/probe-victim" "${work}/home-probe-victim"; do
  mkdir -p "${victim}"
  printf 'private sentinel\n' > "${victim}/renamed"
  "${real_chmod}" 700 "${victim}"
  "${real_chmod}" 600 "${victim}/renamed"
  for probe_command in preflight copy; do
    if [[ "${probe_command}" == preflight ]]; then
      PROBE_SWAP_TARGET="${victim}" reject bash "${utility}" preflight \
        --nas-root "${nas}" --report "${work}/probe-race.json"
    else
      PROBE_SWAP_TARGET="${victim}" reject copy
    fi
    [[ -e "${work}/probe-swap.triggered" ]]
    [[ ! -e "${work}/probe-race.json" && ! -e "${receipt}" ]]
    "${real_python}" - "${victim}" <<'PY'
import os, stat, sys
from pathlib import Path
victim = Path(sys.argv[1])
assert os.listdir(victim) == ['renamed']
assert (victim / 'renamed').read_text() == 'private sentinel\n'
assert stat.S_IMODE(victim.stat().st_mode) == 0o700
assert stat.S_IMODE((victim / 'renamed').stat().st_mode) == 0o600
PY
    probe_path="$(cat "${work}/probe-swap.triggered")"
    rm "${probe_path}" "${work}/probe-swap.triggered"
    rmdir "${probe_path}.saved"
  done
done
NO_CAPACITY=1 reject copy
[[ ! -e "${receipt}" ]]
INEFFECTIVE_CHMOD=1 reject copy
[[ ! -e "${receipt}" ]]
# The actual selected ancestor is a different mount, despite a valid NAS root.
mkdir -p "${destination}"
count="$(wc -l < "${work}/copy.calls")"
WRONG_MOUNT_TARGET="${nas}/_ecoda_backups/corrected_final" reject copy
[[ ! -e "${receipt}" ]]
[[ "$(wc -l < "${work}/copy.calls")" == "${count}" ]]
# Swap an opened selected directory to a scratch link exactly at fchmod.
mkdir -p "${work}/scratch/permission-victim"
"${real_chmod}" 700 "${work}/scratch/permission-victim"
SWAP_PATH="${destination}" SWAP_TARGET="${work}/scratch/permission-victim" reject copy
[[ -e "${work}/swap.triggered" && ! -e "${receipt}" ]]
[[ "$("${real_python}" -c 'import os,stat,sys; print(oct(stat.S_IMODE(os.stat(sys.argv[1]).st_mode)))' \
  "${work}/scratch/permission-victim")" == 0o700 ]]
rm "${destination}" "${work}/swap.triggered"
mv "${destination}.saved" "${destination}"
# Wrapper-timed swaps immediately before descriptor-relative copy creation.
# Neither a root swap nor a selected nested-directory swap can redirect writes.
mkdir -p "${destination}/nested" "${work}/home-copy"
"${real_chmod}" 700 "${work}/home-copy"
for copy_path in "${destination}" "${destination}/nested"; do
  COPY_SWAP_PATH="${copy_path}" COPY_SWAP_TARGET="${work}/home-copy" reject copy
  [[ -e "${work}/copy-swap.triggered" && ! -e "${receipt}" ]]
  "${real_python}" - "${work}/home-copy" <<'PY'
import os, stat, sys
assert os.listdir(sys.argv[1]) == []
assert stat.S_IMODE(os.stat(sys.argv[1]).st_mode) == 0o700
PY
  rm "${copy_path}" "${work}/copy-swap.triggered"
  mv "${copy_path}.saved" "${copy_path}"
done
cp "${work}/files.txt" "${work}/files.saved"
printf '../escape\n' > "${work}/files.txt"
reject copy
mv "${work}/files.saved" "${work}/files.txt"
mv "${work}/source/input.h5ad" "${work}/source/payload.saved"
ln -s payload.saved "${work}/source/input.h5ad"
reject copy
rm "${work}/source/input.h5ad"
mv "${work}/source/payload.saved" "${work}/source/input.h5ad"
copy
verify
cmp "${work}/source/input.h5ad.md5" "${destination}/input.h5ad.md5"
[[ ! -e "${work}/sbatch.calls" ]]
# A fully matching inventory means no streaming copy, even with a new receipt.
count="$(wc -l < "${work}/copy.calls")"
receipt="${ECODA_RUN_ROOT}/manifests/reused.json"
copy
[[ "$(wc -l < "${work}/copy.calls")" == "${count}" ]]
verify
# Receipt verification must inspect the destination mount, not only nas_root.
WRONG_MOUNT_TARGET="${destination}" reject verify
mv "${destination}" "${destination}.saved"
ln -s "${work}/scratch/permission-victim" "${destination}"
reject verify
rm "${destination}"
mv "${destination}.saved" "${destination}"
# Even byte-identical hard-linked selected files are not safe permission targets.
ln "${destination}/input.h5ad" "${work}/external-hardlink"
reject verify
rm "${work}/external-hardlink"
# A file swapped to a home-like target at descriptor chmod cannot change that target.
mkdir "${work}/home"
printf 'private\n' > "${work}/home/private"
"${real_chmod}" 600 "${work}/home/private"
receipt="${ECODA_RUN_ROOT}/manifests/raced-file.json"
SWAP_PATH="${destination}/input.h5ad" SWAP_TARGET="${work}/home/private" reject copy
[[ -e "${work}/swap.triggered" && ! -e "${receipt}" ]]
[[ "$("${real_python}" -c 'import os,stat,sys; print(oct(stat.S_IMODE(os.stat(sys.argv[1]).st_mode)))' \
  "${work}/home/private")" == 0o600 ]]
rm "${destination}/input.h5ad" "${work}/swap.triggered"
mv "${destination}/input.h5ad.saved" "${destination}/input.h5ad"
receipt="${ECODA_RUN_ROOT}/manifests/reused.json"
verify
# Receipt corruption and omitted archive content both fail closed.
printf 'corruption\n' >> "${receipt}"
reject verify
receipt="${ECODA_RUN_ROOT}/manifests/backup.json"
rm "${destination}/input.h5ad.md5"
reject verify
receipt="${ECODA_RUN_ROOT}/manifests/restored.json"
copy
verify
"${real_chmod}" 600 "${destination}/input.h5ad"
reject verify
"${real_chmod}" 666 "${destination}/input.h5ad"
printf 'unlisted\n' > "${destination}/unlisted"
reject verify
receipt="${ECODA_RUN_ROOT}/manifests/unlisted.json"
reject copy
rm "${destination}/unlisted"
mkdir "${destination}/hidden-subtree"
printf 'unlisted\n' > "${destination}/hidden-subtree/payload"
"${real_chmod}" 000 "${destination}/hidden-subtree"
receipt="${ECODA_RUN_ROOT}/manifests/restored.json"
reject verify
receipt="${ECODA_RUN_ROOT}/manifests/hidden-subtree.json"
reject copy
"${real_chmod}" 700 "${destination}/hidden-subtree"
rm "${destination}/hidden-subtree/payload"
rmdir "${destination}/hidden-subtree"
printf 'corrupt\n' > "${destination}/input.h5ad"
receipt="${ECODA_RUN_ROOT}/manifests/changed-selection.json"
MUTATE_SELECTION="${ECODA_RUN_ROOT}/manifests/selection.tsv" reject copy
[[ ! -e "${receipt}" ]]
# Restore the fixture's selection identity, then prove mid-copy source change fails.
printf 'dataset\tcorrected\n' > "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'corrupt again\n' > "${destination}/input.h5ad"
receipt="${ECODA_RUN_ROOT}/manifests/changed-source.json"
MUTATE_SOURCE="${work}/source/input.h5ad" reject copy
[[ ! -e "${receipt}" ]]
# Changed, even correctly rechecksummed, selection invalidates an older receipt.
printf 'other\tcorrected\n' > "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
receipt="${ECODA_RUN_ROOT}/manifests/restored.json"
reject verify
[[ ! -e "${work}/sbatch.calls" ]]
[[ ! -s "${work}/rsync.calls" ]]
echo 'PASS: finite NAS backup fixtures (not live NAS proof)'
