"""Standalone bounded tests for derivative-root/input binding, using text fixtures."""
from pathlib import Path
import hashlib
import json
import sys
import subprocess
import os
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind, load, resolve_input, context_root


def rejects(call, text):
    try:
        call()
    except (ValueError, OSError) as exc:
        assert text in str(exc), str(exc)
    else:
        raise AssertionError(f"expected failure: {text}")


def run():
    with tempfile.TemporaryDirectory(prefix="ecoda-derivative-") as temporary:
        root = Path(temporary).resolve()
        nas = root / "nas"
        nas.mkdir()
        source_dir = root / "canonical"
        source_dir.mkdir()
        source = source_dir / "toy.h5ad"
        original = b"existing canonical source; never replace"
        source.write_bytes(original)
        md5 = hashlib.md5(original).hexdigest()
        proof = source.with_suffix(".h5ad.md5")
        proof.write_text(f"MD5={md5}\nSIZE={len(original)}\nPATH={source}\n")
        manifest = root / "inputs.tsv"
        manifest.write_text(f"Toy\tbenchmark_analysis\t{source}\t{md5}\t{len(original)}\t{proof}\n")
        config = root / "datasets.json"
        config.write_text(json.dumps({"Toy": {"columns": {"label": "biology"},
            "views": {"benchmark_analysis": {"input_file_name": "toy.h5ad",
                                             "output_file_name": "toy.h5ad"}}}}))
        for stage in ("stage3", "stage5"):
            run_root = root / "scratch/_ecoda_runs" / f"{stage}-test"
            control = run_root / "manifests"
            control.mkdir(parents=True)
            (control / "source.manifest").write_text("source snapshot identity\n")
            (control / "runtime.identity").write_text("runtime identity\n")
            selection = control / "selection.tsv"
            selection.write_text("Toy\tbenchmark_analysis\n")
            binding = bind(run_root=run_root, release_id="20261108-confounding-subsets-v1",
                nas_root=nas, stage=stage, config_path=config, selection=selection,
                input_artifacts=manifest)
            value, _ = load(binding, config, verify_content=True)
            assert resolve_input(binding, config, "Toy", "benchmark_analysis") == source
            assert resolve_input(binding, config, "Toy", "benchmark_analysis", directory=True) == source_dir
            context = "datasets" if stage == "stage3" else "benchmark"
            assert context_root(value, context) == run_root / "artifacts" / context
            assert context_root(value, context, nas=True) == nas / "releases" / value["release_id"] / run_root.name / context
            assert not (run_root / "artifacts").exists(), "binding must not create scientific outputs"
            shell_env = dict(os.environ, PROJECT_ROOT=str(ROOT), DATASETS_JSON_FILE=str(config),
                PYTHON_BIN=sys.executable, HPC_SCRATCH_DIR=str(root / "scratch"),
                ECODA_RUN_ROOT=str(run_root), NAS_TARGET_DIR=str(nas),
                ECODA_ARTIFACT_BINDING=str(binding))
            shell = 'source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"\n' \
                'source "$PROJECT_ROOT/src/utils/bash/ecoda_derivative_artifacts.sh"\n' \
                f'ecoda_derivative_load {stage}\n' \
                f'ecoda_derivative_root {context}\n' \
                'ecoda_derivative_input Toy benchmark_analysis file content\n'
            result = subprocess.run(["bash", "-euo", "pipefail", "-c", shell],
                env=shell_env, text=True, capture_output=True)
            assert result.returncode == 0, result.stderr
            assert result.stdout.splitlines() == [str(run_root / "artifacts" / context), str(source)]
            bad_env = dict(shell_env, ECODA_RUN_ROOT=str(root))
            assert subprocess.run(["bash", "-euo", "pipefail", "-c", shell],
                env=bad_env, capture_output=True).returncode != 0
            rejects(lambda: resolve_input(binding, config, "Other", "benchmark_analysis"), "outside derivative selection")
            rejects(lambda: context_root(value, "../canonical"), "invalid derivative output")
            wrong = run_root / "artifacts"
            wrong.symlink_to(source_dir, target_is_directory=True)
            rejects(lambda: context_root(value, context), "not canonical")
            wrong.unlink()
            selection.write_text("Other\tbenchmark_analysis\n")
            rejects(lambda: load(binding, config), "identity changed")
            selection.write_text("Toy\tbenchmark_analysis\n")
            source.write_bytes(b"x" * len(original))
            rejects(lambda: resolve_input(binding, config, "Toy", "benchmark_analysis", verify_content=True), "checksum changed")
            source.write_bytes(original)
            assert source.read_bytes() == original
            rejects(lambda: bind(run_root=run_root, release_id=value["release_id"],
                nas_root=nas, stage=stage, config_path=config, selection=selection,
                input_artifacts=manifest), "File exists")
        extra = manifest.read_text() + manifest.read_text()
        manifest.write_text(extra)
        from src.utils.py.derivative_artifacts import inputs
        rejects(lambda: inputs(manifest, json.loads(config.read_text()), "stage3"), "duplicate input")
    print("PASS: exact derivative roots, immutable canonical inputs, scope/proof checks, and tamper rejection")


if __name__ == "__main__":
    run()
