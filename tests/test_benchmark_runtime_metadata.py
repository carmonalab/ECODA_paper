#!/usr/bin/env python3
"""The sync guard accepts producer batch annotations without weakening binding."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    with tempfile.TemporaryDirectory(prefix="ecoda-runtime-schema-") as temporary:
        directory = Path(temporary)
        output = directory / "scratch/batch_effect/corrected_final/embedding.feather"
        output.parent.mkdir(parents=True)
        output.write_bytes(b"runtime artifact binding fixture\n")
        digest = hashlib.md5(output.read_bytes()).hexdigest()
        metadata = Path(str(output) + ".runtime.json")
        environment = dict(os.environ, PYTHON_BIN=sys.executable,
                           HPC_SCRATCH_DIR=str(directory / "scratch"),
                           NAS_TARGET_DIR=str(directory / "nas"),
                           ANALYSIS_ROOT=str(output.parent),
                           ANALYSIS_NAS_ROOT=str(directory / "nas/batch_effect/corrected_final"),
                           ANALYSIS_VARIANT="corrected_final", ANALYSIS_PASS="corrected",
                           PASS_ARG="corrected")
        payload = dict(schema_version=1, artifact_path=str(output), artifact_md5=digest,
                       dataset="Synthetic", method="MrVI_hvg2000", time_secs=1.0, mem_GB=None)

        def accepted(value):
            metadata.write_text(json.dumps(value) + "\n")
            result = subprocess.run(
                ["/bin/bash", "-c", 'source "$1"; benchmark_validate_runtime_metadata "$2" "$3" "$4"',
                 "runtime-schema", str(ROOT / "src/5_run_benchmark_methods/benchmark_submit_common.sh"),
                 str(metadata), str(output), digest],
                env=environment, capture_output=True, text=True, check=False,
            )
            return result.returncode == 0

        assert accepted(payload), "legacy runtime metadata rejected"
        # Scientific identity validation is separate; the runtime schema must
        # preserve the producer's optional mapping rather than reject its field.
        corrected = dict(payload, batch_contract={"ordered_source_keys": ["technology"]})
        assert accepted(corrected), "producer batch_contract annotation rejected"
        for invalid in (None, [], "technology"):
            assert not accepted(dict(corrected, batch_contract=invalid)), "non-mapping contract accepted"
        assert not accepted(dict(corrected, unexpected=True)), "unknown runtime field accepted"
        missing = dict(corrected)
        del missing["dataset"]
        assert not accepted(missing), "missing required runtime field accepted"
        assert not accepted(dict(corrected, artifact_md5="0" * 32)), "wrong output digest accepted"
        assert not accepted(dict(corrected, artifact_path=str(output) + ".other")), "wrong output path accepted"
    print("PASS legacy/corrected runtime schema, strict fields and output binding")


if __name__ == "__main__":
    main()
