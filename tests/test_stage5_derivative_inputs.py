"""Bound Stage5 inputs/proofs replace no legacy input or NAS artifact."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from unittest.mock import patch

import h5py

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind
from src.utils.py.h5ad_source_identity import (
    _sidecar_fields, _source_sidecar_fields, build_source_identity, verify_source_identity,
)


def proof(path, target=None):
    destination = Path(target or f"{path}.md5")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(f"MD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\nPATH={path}\n")
    return destination


def shell(code, env, *, valid=True):
    prefix = 'set -euo pipefail\nsource "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"\n' \
        'source "$PROJECT_ROOT/src/utils/bash/ecoda_derivative_artifacts.sh"\n'
    result = subprocess.run(["bash", "-c", prefix + code], env=env, text=True, capture_output=True)
    assert (result.returncode == 0) == valid, (result.returncode, result.stdout, result.stderr)
    return result.stdout.strip()


def rejects(call, message):
    try:
        call()
    except (ValueError, OSError) as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError(f"expected rejection: {message}")


def h5ad(path, tag):
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as handle:
        handle.attrs["fixture"] = tag
        obs = handle.create_group("obs"); obs.attrs["encoding-type"] = "dataframe"
        obs.create_dataset("_index", data=[b"c1", b"c2", b"c3"])
        obs.create_dataset("Sample", data=[b"s1", b"s1", b"s2"])


def run():
    with tempfile.TemporaryDirectory(prefix="ecoda-bound-inputs-") as temporary:
        base = Path(temporary).resolve(); scratch = base / "scratch"; nas = base / "nas"; nas.mkdir()
        runs = scratch / "_ecoda_runs"
        producer = runs / "bassez-producer"; consumer = runs / "bl-consumer"
        for run_root, stage in ((producer, "stage2"), (consumer, "stage5")):
            (run_root / "manifests").mkdir(parents=True)
            (run_root / "metadata").write_text(f"STAGE={stage}\nRUN_ID={run_root.name}\nSTATE=OK\n")
        repaired = producer / "artifacts/datasets/Bassez/output/bassez.h5ad"
        legacy = scratch / "Adams/output/adams.h5ad"
        old_bassez = scratch / "Bassez/output/bassez.h5ad"
        h5ad(repaired, "HR-only repair"); h5ad(legacy, "unchanged Adams"); h5ad(old_bassez, "legacy Bassez")
        old_nas = nas / "Bassez/output/bassez.h5ad"; old_nas.parent.mkdir(parents=True)
        old_nas.write_bytes(old_bassez.read_bytes()); proof(old_nas)
        release_nas = nas / "releases/test-release/bassez-producer/datasets/Bassez/output/bassez.h5ad"
        release_nas.parent.mkdir(parents=True); release_nas.write_bytes(repaired.read_bytes()); proof(release_nas)
        assert release_nas.read_bytes() != old_nas.read_bytes()
        alias = base / "home-scratch"; alias.symlink_to(scratch, target_is_directory=True)
        legacy_sidecar = proof(legacy)
        legacy_sidecar.write_text(legacy_sidecar.read_text().replace(str(scratch), str(alias)))
        expected = proof(legacy, base / "receipts/Adams.expected.md5")
        proof(old_bassez)
        config = base / "datasets.json"
        config.write_text(json.dumps({dataset: {"columns": {"sample": "Sample", "label": "biology"},
            "views": {"benchmark_analysis": {"output_file_name": filename}}}
            for dataset, filename in (("Adams", "adams.h5ad"), ("Bassez", "bassez.h5ad"))}))
        env = dict(os.environ, PROJECT_ROOT=str(ROOT), DATASETS_JSON_FILE=str(config),
            PYTHON_BIN=sys.executable, HPC_SCRATCH_DIR=str(scratch), NAS_TARGET_DIR=str(nas),
            ECODA_RUNS_ROOT=str(runs), ECODA_OWNERS_ROOT=str(scratch / "_ecoda_owners"),
            ECODA_RUN_ROOT=str(consumer), ECODA_RUN_ID=consumer.name, REPAIRED=str(repaired))
        shell('ecoda_write_checksum "$REPAIRED"\n'
            'ecoda_write_artifact_record "$REPAIRED" bassez_h5ad_metadata bassez-producer >/dev/null\n'
            'ecoda_artifact_owner_acquire "$REPAIRED" stage2 bassez-producer 0 0 0 >/dev/null\n'
            'ecoda_artifact_owner_set_state "$REPAIRED" OK fixture >/dev/null\n', env)
        protected = {path: (path.read_bytes(), path.stat().st_mode) for path in
            (legacy, legacy_sidecar, expected, old_bassez, Path(f"{old_bassez}.md5"), old_nas, Path(f"{old_nas}.md5"),
             repaired, Path(f"{repaired}.md5"), release_nas, Path(f"{release_nas}.md5"))}
        manifests = consumer / "manifests"
        (manifests / "source.manifest").write_text("pinned source fixture\n")
        (manifests / "runtime.identity").write_text("pinned runtime fixture\n")
        selection = manifests / "selection.tsv"
        selection.write_text("Adams\tbenchmark_analysis\tpseudobulk\nBassez\tbenchmark_analysis\tpseudobulk\n")
        references = base / "input_artifacts.tsv"
        references.write_text("".join(f"{dataset}\tbenchmark_analysis\t{source}\t{hashlib.md5(source.read_bytes()).hexdigest()}\t{source.stat().st_size}\t{receipt}\n"
            for dataset, source, receipt in (("Adams", legacy, expected), ("Bassez", repaired, Path(f"{repaired}.md5")))))
        binding = bind(run_root=consumer, release_id="test-release", nas_root=nas, stage="stage5",
            config_path=config, selection=selection, input_artifacts=references)
        env["ECODA_ARTIFACT_BINDING"] = str(binding)
        lines = (ROOT / "src/5_run_benchmark_methods/stage5_dispatcher.sh").read_text().splitlines()
        functions = []
        for name in ("stage5_recorded_producer", "stage5_validate_derivative_input_row", "stage5_repair_missing_h5ad_sidecars"):
            start = lines.index(f"{name}() {{"); end = next(i for i in range(start + 1, len(lines)) if lines[i] == "}")
            functions.extend(lines[start:end + 1])
        preflight = "\n".join(functions) + '\nRUN_ID="$ECODA_RUN_ID"\nMANIFEST="$ECODA_RUN_ROOT/manifests/selection.tsv"\nstage5_repair_missing_h5ad_sidecars\n'
        shell(preflight, env)
        # The actual old adjacent alias remains rejected by the default strict
        # parser; only the explicitly bound canonical projection is selected.
        rejects(lambda: _sidecar_fields(legacy), "PATH mismatch")
        with patch.dict(os.environ, env):
            payload = build_source_identity(selection, scratch, config)
            identity = manifests / "source_identity.json"; identity.write_text(json.dumps(payload)); proof(identity)
            verify_source_identity(identity, selection, scratch, config)
            rejects(lambda: _source_sidecar_fields(repaired, config, "Adams", "benchmark_analysis", verify_content=True), "outside the bound")
            saved = expected.read_bytes(); expected.write_text(expected.read_text().replace("MD5=", "MD5=0", 1))
            rejects(lambda: build_source_identity(selection, scratch, config), "checksum proof differs")
            expected.write_bytes(saved)
            original = legacy.read_bytes(); legacy.write_bytes(original[:-1] + bytes([original[-1] ^ 1]))
            rejects(lambda: build_source_identity(selection, scratch, config), "MD5 mismatch")
            shell(preflight, env, valid=False)
            legacy.write_bytes(original)
        # Required run-owned producer evidence cannot be replaced by a checksum.
        record = shell('ecoda_artifact_record_path "$REPAIRED" bassez-producer', env)
        record_path = Path(record); saved_record = record_path.with_suffix(".saved")
        record_path.rename(saved_record); shell(preflight, env, valid=False); saved_record.rename(record_path)
        for state in ("ACTIVE", "FAIL"):
            shell(f'ecoda_artifact_owner_set_state "$REPAIRED" {state} fixture >/dev/null', env)
            shell(preflight, env, valid=False)
        shell('ecoda_artifact_owner_set_state "$REPAIRED" OK fixture >/dev/null', env)
        # Missing adjacent receipts are diagnosed without a repair/write.
        missing = legacy_sidecar.with_suffix(".saved")
        legacy_sidecar.rename(missing); shell(preflight, env, valid=False); missing.rename(legacy_sidecar)
        shell(preflight, env)
        for path, (content, mode) in protected.items():
            assert path.read_bytes() == content and path.stat().st_mode == mode, path
    print("PASS: bound repaired/legacy inputs, canonical proof projection, full hashing, producer evidence, foreign/tampered proof rejection, and no old input/NAS writes")


if __name__ == "__main__":
    run()
