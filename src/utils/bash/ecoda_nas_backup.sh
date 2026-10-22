#!/bin/bash
# Finite, non-publishing backups. Login-side only; never mounts or submits work.
set -euo pipefail
utility_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
: "${HPC_SCRATCH_DIR:=${TMPDIR:-/tmp}}"
source "${utility_dir}/ecoda_run_common.sh"
: "${PYTHON_BIN:?set PYTHON_BIN to the configured interpreter}"
export ECODA_RUN_ROOT ECODA_RUN_ID
"${PYTHON_BIN}" - "${utility_dir}/ecoda_run_common.sh" "$@" <<'PY'
import argparse
from contextlib import contextmanager
import secrets
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import stat
import subprocess
import sys
import tempfile

common = sys.argv[1]
EXPECTED_SOURCE = 'smednas515:/volume1/Shared'
FIELDS = {'schema_version', 'state', 'run_id', 'host', 'nas_root', 'mount_source',
          'source_root', 'destination_root', 'selection_path', 'selection_sha256',
          'inventory_path', 'inventory_sha256', 'source_manifest_path',
          'source_manifest_sha256', 'completed_at'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def path(value):
    p = Path(value)
    require(p.is_absolute() and str(p) == value and not any(c in value for c in '\n\r\t'),
            'noncanonical absolute path: ' + value)
    require('..' not in p.parts and '.' not in p.parts, 'unsafe path: ' + value)
    for part in [p, *p.parents]:
        require(not part.is_symlink(), 'symlink path: ' + str(part))
    require(str(p.resolve()) == value, 'redirected path: ' + value)
    return p


def digest(p, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def shared(function, *args):
    subprocess.run(['bash', '-c', 'source "$1"; shift; "$@"', 'nas-backup',
                    common, function, *map(str, args)], check=True)


def strict_checksum(p):
    path(str(p) + '.md5')
    require(p.is_file(), 'missing control file: ' + str(p))
    lines = Path(str(p) + '.md5').read_text().splitlines()
    require(len(lines) == 3 and lines[0] == 'MD5=' + digest(p, 'md5') and
            lines[1] == 'SIZE=' + str(p.stat().st_size) and
            lines[2] == 'PATH=' + str(p), 'invalid strict sidecar: ' + str(p))
    shared('ecoda_validate_checksum', p)


def atomic_json(destination, data, checksum=False):
    path(str(destination))
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.nas-receipt-', dir=str(destination.parent))
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(data, stream, sort_keys=True)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        shared('ecoda_atomic_copy', temporary, destination)
        if checksum:
            shared('ecoda_write_checksum', destination)
    finally:
        os.unlink(temporary)


@contextmanager
def opened(p):
    """Pin every component; never follow a replacement or chmod through a link."""
    descriptors = [os.open('/', os.O_RDONLY | os.O_DIRECTORY)]
    chain = []
    try:
        for index, name in enumerate(p.parts[1:]):
            flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK
            if index < len(p.parts) - 2:
                flags |= os.O_DIRECTORY
            parent = descriptors[-1]
            fd = os.open(name, flags, dir_fd=parent)
            descriptors.append(fd)
            chain.append((parent, name, fd))

        def unchanged():
            for parent, name, fd in chain:
                current = os.stat(name, dir_fd=parent, follow_symlinks=False)
                pinned = os.fstat(fd)
                require((current.st_dev, current.st_ino, current.st_mode) ==
                        (pinned.st_dev, pinned.st_ino, pinned.st_mode),
                        'replaced NAS path chain: ' + str(p))
            info = os.fstat(descriptors[-1])
            require(stat.S_ISDIR(info.st_mode) or
                    (stat.S_ISREG(info.st_mode) and info.st_nlink == 1),
                    'unsafe NAS file type or hard link: ' + str(p))

        unchanged()
        yield descriptors[-1], unchanged
        unchanged()
    finally:
        for fd in reversed(descriptors):
            os.close(fd)


def mode(p, expected, change=False):
    with opened(p) as (fd, unchanged):
        info = os.fstat(fd)
        require(stat.S_ISDIR(info.st_mode) if expected == 0o777 else stat.S_ISREG(info.st_mode),
                'unexpected NAS file type: ' + str(p))
        mount(p)
        unchanged()
        if change:
            os.fchmod(fd, expected)
        unchanged()
        require(stat.S_IMODE(os.fstat(fd).st_mode) == expected,
                'ineffective NAS mode: ' + str(p))


def mount(nas):
    # Access first: an idle autofs map is not the mounted NFS filesystem.
    require(nas.exists(), 'NAS target missing')
    nas.stat()
    result = subprocess.run(['findmnt', '--json', '--target', str(nas),
                             '--output', 'SOURCE,TARGET,FSTYPE'], check=True,
                            text=True, stdout=subprocess.PIPE)
    entries = json.loads(result.stdout)['filesystems']
    require(len(entries) == 1, 'ambiguous NAS mount')
    entry = entries[0]
    require(entry['source'] == EXPECTED_SOURCE and
            entry['fstype'] in ('nfs', 'nfs4'), 'wrong NAS mount identity')
    # SOURCE is authoritative; the same genuine share may mount elsewhere.
    target = path(entry['target'])
    require(nas == target or target in nas.parents, 'NAS root outside mounted share')
    return entry['source']


def destination_mount(nas, destination, rows):
    """Check the existing write ancestors and selected targets, not just the share root."""
    source = mount(nas)
    targets = {destination, *(destination / name for name, _, _ in rows)}
    ancestors = {nas}
    for target in targets:
        ancestors.update(p for p in (target, *target.parents)
                         if p == nas or nas in p.parents)
    for target in sorted(ancestors, key=lambda p: len(p.parts)):
        if target.exists():
            with opened(target) as (_, unchanged):
                require(mount(target) == source, 'destination mount identity changed')
                unchanged()
    return source


def destination_directory(directory):
    if not directory.exists():
        destination_directory(directory.parent)
        with opened(directory.parent) as (parent, unchanged):
            mount(directory.parent)
            unchanged()
            try:
                os.mkdir(directory.name, mode=0o777, dir_fd=parent)
            except FileExistsError:
                pass
            unchanged()
    with opened(directory) as (fd, unchanged):
        require(stat.S_ISDIR(os.fstat(fd).st_mode), 'non-directory backup ancestor')
        mount(directory)
        unchanged()


def copy_file(source, destination, md5, size):
    destination_directory(destination.parent)
    with opened(destination.parent) as (parent, unchanged):
        mount(destination.parent)
        unchanged()
        temporary = '.ecoda-copy-' + secrets.token_hex(16)
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                     0o600, dir_fd=parent)
        try:
            with os.fdopen(fd, 'wb') as output, opened(source) as (source_fd, source_unchanged):
                require(stat.S_ISREG(os.fstat(source_fd).st_mode), 'nonregular copy source')
                identity, copied = hashlib.md5(), 0
                for block in iter(lambda: os.read(source_fd, 1024 * 1024), b''):
                    output.write(block)
                    identity.update(block)
                    copied += len(block)
                source_unchanged()
                require(copied == size and identity.hexdigest() == md5,
                        'source changed during copy: ' + str(source))
                output.flush()
                require(os.fstat(output.fileno()).st_nlink == 1, 'hard-linked copy temporary')
                os.fchmod(output.fileno(), 0o666)
                info = os.fstat(output.fileno())
                require(stat.S_IMODE(info.st_mode) == 0o666 and info.st_nlink == 1,
                        'ineffective or unsafe NAS copy mode')
                os.fsync(output.fileno())
            unchanged()
            os.replace(temporary, destination.name, src_dir_fd=parent, dst_dir_fd=parent)
            os.fsync(parent)
            unchanged()
        finally:
            try:
                os.unlink(temporary, dir_fd=parent)
            except FileNotFoundError:
                pass


def available(nas):
    info = os.statvfs(nas)
    return info.f_bavail * info.f_frsize


def probe(nas):
    with opened(nas) as (parent, parent_unchanged):
        mount(nas)
        parent_unchanged()
        name = '.ecoda-nas-probe-' + secrets.token_hex(16)
        os.mkdir(name, mode=0o700, dir_fd=parent)
        directory = os.open(name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent)

        def unchanged():
            parent_unchanged()
            current = os.stat(name, dir_fd=parent, follow_symlinks=False)
            pinned = os.fstat(directory)
            require((current.st_dev, current.st_ino) == (pinned.st_dev, pinned.st_ino) and
                    stat.S_ISDIR(current.st_mode), 'replaced NAS probe directory')

        try:
            unchanged()
            os.fchmod(directory, 0o777)
            require(stat.S_IMODE(os.fstat(directory).st_mode) == 0o777,
                    'ineffective NAS probe directory mode')
            unchanged()
            payload = os.urandom(128)
            fd = os.open('created', os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                         0o600, dir_fd=directory)
            with os.fdopen(fd, 'w+b') as stream:
                unchanged()
                stream.write(payload)
                stream.flush()
                os.fsync(stream.fileno())
                info = os.fstat(stream.fileno())
                require(stat.S_ISREG(info.st_mode) and info.st_nlink == 1,
                        'unsafe NAS probe file')
                os.fchmod(stream.fileno(), 0o666)
                require(stat.S_IMODE(os.fstat(stream.fileno()).st_mode) == 0o666,
                        'ineffective NAS probe file mode')
                os.rename('created', 'renamed', src_dir_fd=directory, dst_dir_fd=directory)
                renamed = os.stat('renamed', dir_fd=directory, follow_symlinks=False)
                info = os.fstat(stream.fileno())
                require((renamed.st_dev, renamed.st_ino) == (info.st_dev, info.st_ino) and
                        stat.S_ISREG(renamed.st_mode) and info.st_nlink == 1,
                        'replaced NAS probe file')
                stream.seek(0)
                require(stream.read() == payload, 'NAS probe readback failed')
            unchanged()
            require(stat.S_IMODE(os.fstat(directory).st_mode) == 0o777,
                    'changed NAS probe directory mode')
            os.fsync(directory)
        finally:
            try:
                for entry in ('created', 'renamed'):
                    try:
                        os.unlink(entry, dir_fd=directory)
                    except FileNotFoundError:
                        pass
                try:
                    current = os.stat(name, dir_fd=parent, follow_symlinks=False)
                    pinned = os.fstat(directory)
                    if (current.st_dev, current.st_ino) == (pinned.st_dev, pinned.st_ino):
                        os.rmdir(name, dir_fd=parent)
                except FileNotFoundError:
                    pass
            finally:
                os.close(directory)


def inventory(p):
    result = []
    for line in p.read_text().splitlines():
        parts = line.split('\t')
        require(len(parts) == 3, 'inventory requires relative_path, MD5, SIZE')
        name, md5, size = parts
        require(name and not name.startswith('/') and not any(c in name for c in '\r\n\t\\') and
                all(piece not in ('', '.', '..') for piece in name.split('/')),
                'unsafe inventory path')
        require(re.fullmatch('[0-9a-f]{32}', md5) and re.fullmatch('0|[1-9][0-9]*', size),
                'invalid inventory identity')
        result.append((name, md5, int(size)))
    names = [row[0] for row in result]
    require(result and len(names) == len(set(names)), 'empty or duplicate inventory')
    require(not any(parent in set(names) for name in names
                    for parent in map(str, Path(name).parents) if parent != '.'),
            'overlapping inventory entries')
    return result


def matches(p, md5, size):
    path(str(p))
    return p.is_file() and p.stat().st_size == size and digest(p, 'md5') == md5


def bindings(selection, inv, manifest, run_id, run_root):
    require(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*', run_id), 'invalid promotion run ID')
    require(run_root.name == run_id, 'promotion root/run ID mismatch')
    for p in (selection, inv):
        shared('ecoda_validate_run_owned_path', p, run_root)
        require(run_root in p.parents, 'control file outside promotion run')
        strict_checksum(p)
    # The original immutable manifest is bound, not rewritten for backup paths.
    require(manifest.is_file() and not (manifest.stat().st_mode & 0o222),
            'source manifest must be immutable')
    source = dict(line.split('=', 1) for line in manifest.read_text().splitlines())
    keys = {'FORMAT', 'SOURCE_ROOT', 'SOURCE_COMMIT', 'SOURCE_ARCHIVE_PATH',
            'SOURCE_ARCHIVE_SHA256', 'CONFIG_HELPER_SHA256', 'DATASETS_SHA256',
            'PIXI_TOML_SHA256', 'PIXI_LOCK_SHA256', 'AUX_ROOT', 'SCGATE_DB_BRANCH'}
    require(set(source) == keys and len(manifest.read_text().splitlines()) == 11 and
            source['FORMAT'] == '1' and re.fullmatch('[0-9a-f]{40}', source['SOURCE_COMMIT']),
            'invalid immutable source manifest')
    root = path(source['SOURCE_ROOT'])
    snapshot = root.parent
    require(root.name == 'tree' and manifest == snapshot / 'identity/source.manifest' and
            path(source['SOURCE_ARCHIVE_PATH']) == snapshot / 'identity/source.tar' and
            (snapshot / 'COMPLETE').is_file() and (root / 'aux').is_dir() and
            not (root / '.git').exists(), 'source manifest is not a completed immutable snapshot')
    require(source['AUX_ROOT'] == str(root / 'aux'), 'wrong frozen auxiliary root')
    checks = [(path(source['SOURCE_ARCHIVE_PATH']), source['SOURCE_ARCHIVE_SHA256']),
              (root / 'config_helper.R', source['CONFIG_HELPER_SHA256']),
              (root / 'datasets.json', source['DATASETS_SHA256']),
              (root / 'pixi.toml', source['PIXI_TOML_SHA256']),
              (root / 'pixi.lock', source['PIXI_LOCK_SHA256'])]
    for p, expected in checks:
        path(str(p))
        require(p.is_file() and re.fullmatch('[0-9a-f]{64}', expected) and digest(p) == expected,
                'immutable source identity mismatch: ' + str(p))
    return {key: digest(p) for key, p in
            [('selection_sha256', selection), ('inventory_sha256', inv),
             ('source_manifest_sha256', manifest)]}


def destination_scope(nas, destination, run_id, source):
    boundary = nas / '_ecoda_backups' / 'corrected_final' / run_id
    require(destination == boundary or boundary in destination.parents,
            'destination is not an isolated run-owned NAS backup')
    require(source != destination and source not in destination.parents and
            destination not in source.parents, 'source/destination overlap')


def traversal_error(error):
    raise error


def inspect_destination(destination, rows, change=False):
    names = {row[0] for row in rows}
    allowed_dirs = {str(parent) for name in names for parent in Path(name).parents}
    if destination.exists():
        for directory, dirs, files in os.walk(destination, followlinks=False, onerror=traversal_error):
            for name in dirs + files:
                p = path(str(Path(directory) / name))
                require(not p.is_symlink(), 'symlink in backup')
                if p.is_dir():
                    require(str(p.relative_to(destination)) in allowed_dirs,
                            'unlisted directory in backup: ' + str(p))
                if not p.is_dir():
                    require(p.is_file() and str(p.relative_to(destination)) in names,
                            'unlisted file in backup: ' + str(p))
    directories = {destination}
    for name, md5, size in rows:
        p = destination / name
        require(matches(p, md5, size), 'incomplete or corrupt backup: ' + str(p))
        mode(p, 0o666, change)
        directories.update(parent for parent in p.parents if parent == destination or destination in parent.parents)
    for directory in directories:
        mode(directory, 0o777, change)


parser = argparse.ArgumentParser(description=__doc__)
commands = parser.add_subparsers(dest='command', required=True)
p = commands.add_parser('preflight')
p.add_argument('--nas-root', required=True)
p.add_argument('--report', required=True)
p = commands.add_parser('copy')
for name in ('source-root', 'nas-root', 'destination-root', 'files-from', 'inventory',
             'receipt', 'selection-file', 'source-manifest'):
    p.add_argument('--' + name, required=True)
p = commands.add_parser('verify')
p.add_argument('--receipt', required=True)
args = parser.parse_args(sys.argv[2:])
try:
    if args.command == 'preflight':
        nas, report = path(args.nas_root), path(args.report)
        mount_source = mount(nas)
        probe(nas)
        atomic_json(report, {'mount_source': mount_source, 'nas_root': str(nas),
                             'available_bytes': available(nas)})
    elif args.command == 'copy':
        source, nas, destination = map(path, (args.source_root, args.nas_root, args.destination_root))
        selection, inv, manifest, receipt, files = map(path, (args.selection_file, args.inventory,
                                                            args.source_manifest, args.receipt, args.files_from))
        run_id = os.environ.get('ECODA_RUN_ID', '')
        run_root = path(os.environ.get('ECODA_RUN_ROOT', ''))
        require(run_root in receipt.parents, 'receipt outside promotion run')
        destination_scope(nas, destination, run_id, source)
        pins = bindings(selection, inv, manifest, run_id, run_root)
        rows = inventory(inv)
        require(files.read_text().splitlines() == [row[0] for row in rows],
                'files-from differs from exact inventory')
        files_digest = digest(files)
        mount_source = destination_mount(nas, destination, rows)
        probe(nas)
        # Reject unlisted existing content before writing anything to the backup.
        names = {row[0] for row in rows}
        allowed_dirs = {str(parent) for name in names for parent in Path(name).parents}
        if destination.exists():
            for directory, dirs, entries in os.walk(destination, followlinks=False, onerror=traversal_error):
                for name in dirs + entries:
                    p = path(str(Path(directory) / name))
                    relative = str(p.relative_to(destination))
                    require((p.is_dir() and relative in allowed_dirs) or
                            (p.is_file() and relative in names), 'unlisted backup content')
        missing = []
        for name, md5, size in rows:
            require(matches(source / name, md5, size), 'source changed: ' + name)
            if not matches(destination / name, md5, size):
                missing.append((name, size))
        require(available(nas) >= sum(size for _, size in missing), 'insufficient NAS capacity')
        # Old successful receipts must not survive an interrupted replacement attempt.
        require(not receipt.exists() and not Path(str(receipt) + '.md5').exists(),
                'receipt already exists; verify it or use a new run-owned receipt')
        destination_directory(destination)
        run_boundary = nas / '_ecoda_backups' / 'corrected_final' / run_id
        current = destination
        while current == run_boundary or run_boundary in current.parents:
            mode(current, 0o777, True)
            current = current.parent
        # Shared directory names may be replaced at any time. Never hand rsync
        # a NAS pathname; stream and publish each selected file through its dirfd.
        missing_names = {name for name, _ in missing}
        for name, md5, size in rows:
            if name in missing_names:
                copy_file(source / name, destination / name, md5, size)
        for name, md5, size in rows:
            require(matches(source / name, md5, size), 'source changed during copy: ' + name)
        require(digest(files) == files_digest, 'files-from changed during copy')
        require(bindings(selection, inv, manifest, run_id, run_root) == pins,
                'control identity changed during copy')
        require(destination_mount(nas, destination, rows) == mount_source,
                'mount changed during copy')
        inspect_destination(destination, rows, True)
        data = dict(schema_version=1, state='OK', run_id=run_id, host=socket.gethostname(),
                    nas_root=str(nas), mount_source=mount_source, source_root=str(source),
                    destination_root=str(destination), selection_path=str(selection),
                    inventory_path=str(inv), source_manifest_path=str(manifest),
                    completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), **pins)
        atomic_json(receipt, data, True)
    else:
        receipt = path(args.receipt)
        strict_checksum(receipt)
        data = json.loads(receipt.read_text())
        require(set(data) == FIELDS and type(data['schema_version']) is int and
                data['schema_version'] == 1 and data['state'] == 'OK', 'invalid receipt schema/state')
        require(isinstance(data['host'], str) and data['host'], 'missing receipt host')
        completed = datetime.datetime.fromisoformat(data['completed_at'])
        require(completed.tzinfo is not None and completed <= datetime.datetime.now(datetime.timezone.utc),
                'invalid receipt completion time')
        selection, inv, manifest = map(path, (data['selection_path'], data['inventory_path'],
                                             data['source_manifest_path']))
        # Resolve ownership from the standard run/manifests layout, not caller environment.
        candidates = [p for p in selection.parents if p.name == data['run_id']]
        require(len(candidates) == 1, 'ambiguous promotion run root')
        run_root = candidates[0]
        require(run_root in receipt.parents, 'receipt outside promotion run')
        pins = bindings(selection, inv, manifest, data['run_id'], run_root)
        require(all(data[k] == v for k, v in pins.items()), 'stale receipt bindings')
        nas, destination, source = map(path, (data['nas_root'], data['destination_root'], data['source_root']))
        destination_scope(nas, destination, data['run_id'], source)
        rows = inventory(inv)
        require(destination_mount(nas, destination, rows) == data['mount_source'],
                'receipt mount identity changed')
        inspect_destination(destination, rows)
        run_boundary = nas / '_ecoda_backups' / 'corrected_final' / data['run_id']
        current = destination
        while current == run_boundary or run_boundary in current.parents:
            mode(current, 0o777)
            current = current.parent
except (ValueError, OSError, KeyError, TypeError, subprocess.CalledProcessError) as error:
    print('ERROR: NAS backup: ' + str(error), file=sys.stderr)
    sys.exit(1)
PY
