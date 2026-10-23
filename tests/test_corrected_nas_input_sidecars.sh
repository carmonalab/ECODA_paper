#!/bin/bash
# The corrected-final publisher preserves the source PATH in NAS sidecars.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
PYTHON="${PYTHON_BIN:-${ROOT}/.pixi/envs/default/bin/python}"
"${PYTHON}" - "${ROOT}" <<'PY'
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

root = Path(sys.argv[1])
source = (root / 'src/5_run_benchmark_methods/stage5_dispatcher.sh').read_text()
function = re.search(r'^stage5_repair_missing_h5ad_sidecars\(\) \{\n.*?^\}', source, re.M | re.S)
assert function
with tempfile.TemporaryDirectory(prefix='ecoda-nas-sidecars-') as temporary:
    base = Path(temporary).resolve()
    scratch, nas = base / 'scratch', base / 'nas'
    scratch.mkdir()
    nas.mkdir()
    original, mirror = scratch / 'input.h5ad', nas / 'input.h5ad'
    payload = b'checksum contract fixture\n'
    original.write_bytes(payload)
    mirror.write_bytes(payload)
    manifest = base / 'selection.tsv'
    manifest.write_text('Joanito\tbatch_effect_corrected\tbatch_effect_corrected\n')
    run = scratch / '_ecoda_runs/sidecar_consumer'
    (run / 'manifests').mkdir(parents=True)
    env = dict(os.environ, HPC_SCRATCH_DIR=str(scratch), NAS_TARGET_DIR=str(nas),
        ECODA_RUN_ROOT=str(run), RUN_ID='sidecar_consumer', PASS_ARG='corrected',
        ANALYSIS_VARIANT='corrected_final', MANIFEST=str(manifest), BENCHMARK_MATRIX_TEST='0')
    prelude = 'set -euo pipefail\nsource ' + shlex.quote(str(root / 'src/utils/bash/ecoda_run_common.sh')) + '\n'
    def shell(code):
        return subprocess.run(['/bin/bash', '-c', prelude + code], env=env, text=True, capture_output=True)
    checked = shell('ecoda_write_checksum ' + shlex.quote(str(original)))
    assert checked.returncode == 0, checked.stderr
    original_sidecar = Path(str(original) + '.md5')
    mirror_sidecar = Path(str(mirror) + '.md5')
    source_record = original_sidecar.read_bytes()
    mirror_sidecar.write_bytes(source_record)
    setup = ('stage5_input_path() { printf "%s\\n" ' + shlex.quote(str(original)) + '; }\n'
        'stage5_output_path_for_root() { printf "%s\\n" ' + shlex.quote(str(mirror)) + '; }\n')
    def validate(variant='corrected_final'):
        return shell('ANALYSIS_VARIANT=' + shlex.quote(variant) + '\n' + setup + function.group() + '\nstage5_repair_missing_h5ad_sidecars\n')
    before = {p: (p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns)
              for p in (original, mirror, original_sidecar, mirror_sidecar)}
    accepted = validate()
    assert accepted.returncode == 0, accepted.stderr
    assert all((p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns) == state for p, state in before.items())
    assert validate('other').returncode != 0, 'source-PATH exception leaked outside corrected-final'
    mirror.write_bytes(b'X' + payload[1:])
    assert validate().returncode != 0, 'same-size NAS corruption accepted'
    mirror.write_bytes(payload)
    mirror_sidecar.write_bytes(source_record.replace(str(original).encode(), str(base / 'foreign.h5ad').encode()))
    assert validate().returncode != 0, 'unrelated provenance PATH accepted'
    mirror_sidecar.write_bytes(source_record.replace(b'SIZE=' + str(len(payload)).encode(), b'SIZE=1'))
    assert validate().returncode != 0, 'incorrect NAS size record accepted'
    checked = shell('ecoda_write_checksum ' + shlex.quote(str(mirror)))
    assert checked.returncode == 0, checked.stderr
    assert validate().returncode == 0, 'native NAS-PATH sidecar rejected'
    print('PASS corrected-final copied and native NAS sidecars; corrupt bytes, foreign PATH and wrong size rejected; valid inputs unchanged')
PY
