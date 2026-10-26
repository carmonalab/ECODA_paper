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
import io
# Count bytes actually returned by binary content reads, including the Python
# interpreter used by the real hashing PATH adapter below and nested verify.
def record_read(path, count):
    targets = os.environ.get('READ_BUDGET_PATHS', '').splitlines()
    if str(path) in targets and count:
        with open(os.environ['READ_BUDGET_LOG'], 'a') as log:
            log.write(str(path) + '\t' + str(count) + '\n')
    trigger = os.environ.get('MUTATE_AFTER_READ_PATH')
    marker = os.environ['FIXTURE_ROOT'] + '/read-mutation.triggered'
    if count == 0 and str(path) == trigger and not os.path.exists(marker):
        with open(os.environ['MUTATE_AFTER_READ_TARGET'], 'r+b') as stream:
            stream.write(b'!')
        with open(marker, 'w') as stream:
            stream.write(str(path))
original_io_open, original_read = io.open, os.read
class CountedReader:
    def __init__(self, stream, path):
        self.stream, self.path = stream, path
    def __enter__(self):
        self.stream.__enter__()
        return self
    def __exit__(self, *args):
        return self.stream.__exit__(*args)
    def __getattr__(self, name):
        return getattr(self.stream, name)
    def read(self, *args):
        data = self.stream.read(*args)
        record_read(self.path, len(data))
        return data
def counted_open(file, mode='r', *args, **kwargs):
    stream = original_io_open(file, mode, *args, **kwargs)
    return CountedReader(stream, file) if mode == 'rb' else stream
def counted_read(fd, size):
    data = original_read(fd, size)
    if os.environ.get('READ_BUDGET_PATHS') or os.environ.get('MUTATE_AFTER_READ_PATH'):
        info = os.fstat(fd)
        targets = set(os.environ.get('READ_BUDGET_PATHS', '').splitlines())
        targets.update(os.environ.get('MUTATE_AFTER_READ_PATH', '').splitlines())
        for target in targets:
            current = os.stat(target)
            if (info.st_dev, info.st_ino) == (current.st_dev, current.st_ino):
                record_read(target, len(data))
    return data
io.open, os.read = counted_open, counted_read
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
            if os.environ.get('SWAP_HARDLINK') == '1':
                os.link(os.environ['SWAP_TARGET'], attack)
            else:
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
original_rename, original_link = os.rename, os.link
def inject_publication_race():
    target = os.environ.get('PUBLICATION_RACE_TARGET')
    marker = os.environ['FIXTURE_ROOT'] + '/publication-race.triggered'
    if target and not os.path.exists(marker):
        with open(target + '.racing', 'w') as stream:
            stream.write('unexpected concurrent canonical bytes\n')
        os.replace(target + '.racing', target)
        with open(marker, 'w') as stream:
            stream.write(target)
def race_rename(src, dst, *, src_dir_fd=None, dst_dir_fd=None):
    if dst == 'predecessor' and os.environ.get('PUBLICATION_RACE_KIND') == 'conflict':
        inject_publication_race()
    return original_rename(src, dst, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd)
def race_link(src, dst, *, src_dir_fd=None, dst_dir_fd=None, follow_symlinks=True):
    late = os.environ.get('LATE_RECEIPT_TARGET')
    if late and dst == os.path.basename(late) and str(src).startswith('.nas-receipt-'):
        with open(late, 'x') as stream:
            stream.write('late bad control file\n')
    if str(src).startswith('.ecoda-copy-') and os.environ.get('PUBLICATION_RACE_KIND') == 'absent':
        inject_publication_race()
    return original_link(src, dst, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd,
                         follow_symlinks=follow_symlinks)
os.rename, os.link = race_rename, race_link
PY
# Real streaming MD5 (not an echoed checksum): makes shell-validator reads
# observable with the same cross-process byte counter on either macOS or Linux.
cat > "${work}/bin/md5sum" <<'SH'
#!/bin/bash
exec "${PYTHON_BIN}" - "$@" <<'PY'
import hashlib, sys
from pathlib import Path
identity = hashlib.md5()
with Path(sys.argv[1]).open('rb') as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        identity.update(block)
print(identity.hexdigest() + '  ' + sys.argv[1])
PY
SH
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
printf '{"fixture":{"folder_name":"fixture-data","views":{"batch_effect_corrected":{"output_file_name":"corrected.h5ad"}}}}\n' \
  > "${work}/snapshot/tree/datasets.json"
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
# Standalone strict receipt verification hashes its content exactly once.
READ_BUDGET_PATHS="${receipt}" READ_BUDGET_LOG="${work}/strict-reads" verify
"${real_python}" - "${receipt}" "${work}/strict-reads" <<'PY'
import sys
from pathlib import Path
receipt, log = map(Path, sys.argv[1:])
assert sum(int(line.split('\t')[1]) for line in log.read_text().splitlines()) == receipt.stat().st_size
PY
# A strict record changed immediately after hashing still fails confirmation.
cp "${receipt}.md5" "${work}/strict-sidecar.saved"
MUTATE_AFTER_READ_PATH="${receipt}" MUTATE_AFTER_READ_TARGET="${receipt}.md5" reject verify
[[ -e "${work}/read-mutation.triggered" ]]
cp "${work}/strict-sidecar.saved" "${receipt}.md5"
rm "${work}/read-mutation.triggered"
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
# Canonical inputs publish only configured selected companions, never walk NAS.
printf 'fixture\tbatch_effect_corrected\tbatch_effect_corrected\n' > "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
input_name='fixture/output/corrected.h5ad'
mkdir -p "${work}/source/fixture/output" "${nas}/fixture/output" "${nas}/unrelated"
printf 'corrected canonical bytes\n' > "${work}/source/${input_name}"
ecoda_write_checksum "${work}/source/${input_name}"
printf 'unrelated private bytes\n' > "${nas}/unrelated/private"
"${real_chmod}" 700 "${nas}/unrelated"
"${real_chmod}" 600 "${nas}/unrelated/private"
printf 'old canonical bytes\n' > "${nas}/${input_name}"
printf 'old sidecar bytes\n' > "${nas}/${input_name}.md5"
mkdir -p "${work}/archive-source"
cp "${nas}/${input_name}" "${work}/archive-source/old"
cp "${nas}/${input_name}.md5" "${work}/archive-source/old.md5"
archive_destination="${nas}/_ecoda_backups/corrected_final/${ECODA_RUN_ID}/before"
archive_receipt="${ECODA_RUN_ROOT}/manifests/archive.receipt.json"
mapping="${ECODA_RUN_ROOT}/manifests/archive_mapping.tsv"
printf 'old\nold.md5\n' > "${work}/archive-files"
for name in old old.md5; do
  printf '%s\t%s\t%s\n' "${name}" "$(ecoda_md5_file "${work}/archive-source/${name}")" \
    "$(wc -c < "${work}/archive-source/${name}" | tr -d '[:space:]')"
done > "${ECODA_RUN_ROOT}/manifests/archive_inventory.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/archive_inventory.tsv"
bash "${utility}" copy --source-root "${work}/archive-source" --nas-root "${nas}" \
  --destination-root "${archive_destination}" --files-from "${work}/archive-files" \
  --inventory "${ECODA_RUN_ROOT}/manifests/archive_inventory.tsv" --receipt "${archive_receipt}" \
  --selection-file "${ECODA_RUN_ROOT}/manifests/selection.tsv" --source-manifest "${manifest}"
printf 'original_path\tarchive_path\tmd5\tsize\n' > "${mapping}"
for suffix in '' .md5; do
  printf '%s\t%s\t%s\t%s\n' "${nas}/${input_name}${suffix}" "${archive_destination}/old${suffix}" \
    "$(ecoda_md5_file "${nas}/${input_name}${suffix}")" \
    "$(wc -c < "${nas}/${input_name}${suffix}" | tr -d '[:space:]')" >> "${mapping}"
done
ecoda_write_checksum "${mapping}"
printf '%s\n%s.md5\n' "${input_name}" "${input_name}" > "${work}/input-files"
input_inventory="${ECODA_RUN_ROOT}/manifests/input_inventory.tsv"
make_input_inventory() {
  while IFS= read -r name; do
    printf '%s\t%s\t%s\n' "${name}" "$(ecoda_md5_file "${work}/source/${name}")" \
      "$(wc -c < "${work}/source/${name}" | tr -d '[:space:]')"
  done < "${work}/input-files" > "${input_inventory}"
  ecoda_write_checksum "${input_inventory}"
}
make_input_inventory
receipt="${ECODA_RUN_ROOT}/manifests/inputs.receipt.json"
publish() {
  bash "${utility}" publish-inputs --source-root "${work}/source" --nas-root "${nas}" \
    --destination-root "${nas}" --files-from "${work}/input-files" --inventory "${input_inventory}" \
    --receipt "${receipt}" --selection-file "${ECODA_RUN_ROOT}/manifests/selection.tsv" \
    --source-manifest "${manifest}" --archive-receipt "${archive_receipt}" --archive-mapping "${mapping}"
}
# Mutations at the actual publication syscall, after all pathname prechecks,
# cannot overwrite unarchived bytes or create a success receipt.
quarantine="$("${real_python}" -c \
  'import hashlib,sys; from pathlib import Path; p=Path(sys.argv[1]); print(p.parent / (".ecoda-displaced-" + hashlib.sha256(p.name.encode()).hexdigest()))' \
  "${nas}/${input_name}")"
PUBLICATION_RACE_TARGET="${nas}/${input_name}" PUBLICATION_RACE_KIND=conflict reject publish
[[ ! -e "${receipt}" && -e "${work}/publication-race.triggered" ]]
[[ "$(cat "${quarantine}/predecessor")" == 'unexpected concurrent canonical bytes' ]]
[[ ! -e "${nas}/${input_name}" ]]
# An interrupted/failed displacement remains a fail-closed barrier on retry.
reject publish
[[ "$(cat "${quarantine}/predecessor")" == 'unexpected concurrent canonical bytes' ]]
rm "${quarantine}/predecessor" "${work}/publication-race.triggered"
rmdir "${quarantine}"
PUBLICATION_RACE_TARGET="${nas}/${input_name}" PUBLICATION_RACE_KIND=absent reject publish
[[ ! -e "${receipt}" && -e "${work}/publication-race.triggered" ]]
[[ "$(cat "${nas}/${input_name}")" == 'unexpected concurrent canonical bytes' ]]
reject publish
[[ "$(cat "${nas}/${input_name}")" == 'unexpected concurrent canonical bytes' ]]
rm "${work}/publication-race.triggered"
cp "${work}/archive-source/old" "${nas}/${input_name}"
# A changed predecessor is not authorization to overwrite even with sealed archive.
printf 'not archived\n' > "${nas}/${input_name}"
reject publish
[[ ! -e "${receipt}" ]]
cp "${work}/archive-source/old" "${nas}/${input_name}"
# Missing proof, bad scope, omitted companion, duplicate selection all fail closed.
mv "${mapping}.md5" "${mapping}.saved"
reject publish
mv "${mapping}.saved" "${mapping}.md5"
cp "${work}/input-files" "${work}/input-files.saved"
printf 'nested/payload\n' >> "${work}/input-files"
make_input_inventory
reject publish
cp "${work}/input-files.saved" "${work}/input-files"
printf '%s\n' "${input_name}" > "${work}/input-files"
make_input_inventory
reject publish
cp "${work}/input-files.saved" "${work}/input-files"
make_input_inventory
cp "${ECODA_RUN_ROOT}/manifests/selection.tsv" "${work}/canonical-selection"
cat "${work}/canonical-selection" >> "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
reject publish
cp "${work}/canonical-selection" "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
printf 'fixture\tbatch_effect_uncorrected\tbatch_effect_corrected\n' > "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
reject publish
cp "${work}/canonical-selection" "${ECODA_RUN_ROOT}/manifests/selection.tsv"
ecoda_write_checksum "${ECODA_RUN_ROOT}/manifests/selection.tsv"
WRONG_MOUNT_TARGET="${nas}/fixture" reject publish
ln "${nas}/${input_name}" "${work}/canonical-hardlink"
reject publish
rm "${work}/canonical-hardlink"
mv "${nas}/fixture/output" "${nas}/fixture/output.saved"
ln -s "${work}/home-copy" "${nas}/fixture/output"
reject publish
rm "${nas}/fixture/output"
mv "${nas}/fixture/output.saved" "${nas}/fixture/output"
# Directory swap while creating the copy cannot redirect canonical publication.
COPY_SWAP_PATH="${nas}/fixture/output" COPY_SWAP_TARGET="${work}/home-copy" reject publish
[[ ! -e "${receipt}" && -e "${work}/copy-swap.triggered" ]]
rm "${nas}/fixture/output" "${work}/copy-swap.triggered"
mv "${nas}/fixture/output.saved" "${nas}/fixture/output"
# A byte-identical source hardlink and malformed producer sidecar fail closed.
ln "${work}/source/${input_name}" "${work}/source-hardlink"
reject publish
rm "${work}/source-hardlink"
cp "${work}/source/${input_name}.md5" "${work}/input-sidecar.saved"
printf 'MD5=wrong\n' > "${work}/source/${input_name}.md5"
make_input_inventory
reject publish
cp "${work}/input-sidecar.saved" "${work}/source/${input_name}.md5"
make_input_inventory
publish
# Reused matching bytes never copy, and a descriptor-time hardlink swap cannot
# chmod unrelated private files even when a selected mode needs repair.
receipt="${ECODA_RUN_ROOT}/manifests/hardlink-race-input.receipt.json"
"${real_chmod}" 600 "${nas}/${input_name}"
SWAP_PATH="${nas}/${input_name}" SWAP_TARGET="${work}/home/private" SWAP_HARDLINK=1 reject publish
[[ ! -e "${receipt}" && -e "${work}/swap.triggered" ]]
rm "${nas}/${input_name}" "${work}/swap.triggered"
mv "${nas}/${input_name}.saved" "${nas}/${input_name}"
"${real_chmod}" 666 "${nas}/${input_name}"
publish
verify
cmp "${work}/source/${input_name}" "${nas}/${input_name}"
cmp "${work}/source/${input_name}.md5" "${nas}/${input_name}.md5"
# A protected receipt cannot be overwritten, including malformed prior receipts.
reject publish
receipt="${ECODA_RUN_ROOT}/manifests/bad-input.receipt.json"
printf 'bad existing receipt\n' > "${receipt}"
reject publish
[[ "$(cat "${receipt}")" == 'bad existing receipt' ]]
# A non-cooperating late writer cannot have its receipt or sidecar overwritten.
receipt="${ECODA_RUN_ROOT}/manifests/late-input.receipt.json"
LATE_RECEIPT_TARGET="${receipt}" reject publish
[[ "$(cat "${receipt}")" == 'late bad control file' ]]
reject verify
receipt="${ECODA_RUN_ROOT}/manifests/late-sidecar-input.receipt.json"
LATE_RECEIPT_TARGET="${receipt}.md5" reject publish
[[ ! -e "${receipt}" && "$(cat "${receipt}.md5")" == 'late bad control file' ]]
# Matching targets are reused; an absent selected companion can be resumed.
receipt="${ECODA_RUN_ROOT}/manifests/matching-input.receipt.json"
count="$(wc -l < "${work}/copy.calls")"
READ_BUDGET_PATHS="$(printf '%s\n%s\n' "${work}/source/${input_name}" "${nas}/${input_name}")" \
  READ_BUDGET_LOG="${work}/publisher-reads" publish
"${real_python}" - "${work}/source/${input_name}" "${nas}/${input_name}" \
  "${work}/publisher-reads" <<'PY'
import sys
from pathlib import Path
source, nas, log = map(Path, sys.argv[1:])
totals = {str(source): 0, str(nas): 0}
for line in log.read_text().splitlines():
    name, count = line.split('\t')
    totals[name] += int(count)
assert totals == {str(source): 2 * source.stat().st_size, str(nas): 2 * nas.stat().st_size}, totals
PY
[[ "$(wc -l < "${work}/copy.calls")" == "${count}" ]]
verify
receipt="${ECODA_RUN_ROOT}/manifests/chmod-input.receipt.json"
"${real_chmod}" 600 "${nas}/${input_name}"
INEFFECTIVE_CHMOD=1 reject publish
[[ ! -e "${receipt}" ]]
"${real_chmod}" 666 "${nas}/${input_name}"
# Even a no-copy publication must re-read mutable source and NAS bytes before
# sealing success. Mutate just after the initial NAS content read.
for target in "${work}/source/${input_name}" "${nas}/${input_name}"; do
  cp "${target}" "${work}/read-mutation.saved"
  receipt="${ECODA_RUN_ROOT}/manifests/read-mutation-input.receipt.json"
  MUTATE_AFTER_READ_PATH="${nas}/${input_name}" MUTATE_AFTER_READ_TARGET="${target}" reject publish
  [[ -e "${work}/read-mutation.triggered" && ! -e "${receipt}" ]]
  cp "${work}/read-mutation.saved" "${target}"
  rm "${work}/read-mutation.triggered"
done
receipt="${ECODA_RUN_ROOT}/manifests/matching-input.receipt.json"
rm "${nas}/${input_name}.md5"
reject verify
receipt="${ECODA_RUN_ROOT}/manifests/resumed-input.receipt.json"
publish
verify
printf 'changed canonical bytes\n' > "${nas}/${input_name}"
reject verify
"${real_python}" - "${nas}" "${work}/home/private" <<'PY'
import os, stat, sys
from pathlib import Path
root = Path(sys.argv[1])
assert (root / 'unrelated/private').read_text() == 'unrelated private bytes\n'
assert stat.S_IMODE((root / 'unrelated').stat().st_mode) == 0o700
assert stat.S_IMODE((root / 'unrelated/private').stat().st_mode) == 0o600
assert not list(root.glob('.ecoda-nas-probe-*'))
assert stat.S_IMODE(Path(sys.argv[2]).stat().st_mode) == 0o600
assert Path(sys.argv[2]).read_text() == 'private\n'
PY
[[ ! -e "${work}/sbatch.calls" ]]
[[ ! -s "${work}/rsync.calls" ]]
echo 'PASS: finite NAS backup fixtures (not live NAS proof)'
