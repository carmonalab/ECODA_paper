"""Exact partial-family scope, independent BL preparation, and request tampering."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.benchmark_request_scope import (AUTHOR_HR, active_scope, normalized_rows,
    read_request, selected_files, tokens_for)
from src.utils.py.derivative_artifacts import bind


def rejects(call, message):
    try:
        call()
    except (ValueError, OSError) as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError(f"expected rejection: {message}")


def run():
    with tempfile.TemporaryDirectory(prefix="ecoda-request-") as temporary:
        base = Path(temporary).resolve()
        config = {dataset: {"columns": {"label": "biology"},
            "views": {"benchmark_analysis": {"output_file_name": "source.h5ad"}}}
            for dataset in ("Stephenson", "Adams", "Bassez")}
        request = base / "request.tsv"
        original = "Stephenson\tbenchmark_analysis\tall\nAdams\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n" \
            "Bassez\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n" \
            "Bassez\tbenchmark_analysis\tpseudobulk\tCT_HR_hvg2000\n" \
            "Bassez\tbenchmark_analysis\tpseudobulk\tCT_HR_hvg500\n" + \
            "".join(f"Bassez\tbenchmark_analysis\tcomposition\t{token}\n" for token in AUTHOR_HR) + \
            "Bassez\tbenchmark_analysis\tscpoli\thvg2000_highres_dims15\n"
        request.write_text(original)
        scope = read_request(request, config)
        assert len(AUTHOR_HR) == 22
        assert tokens_for(scope, "Stephenson", "benchmark_analysis", "pseudobulk") is None
        assert tokens_for(scope, "Stephenson", "benchmark_analysis", "prepare_pseudobulk") is None
        assert tokens_for(scope, "Adams", "benchmark_analysis", "prepare_pseudobulk") == "ABSENT"
        assert tokens_for(scope, "Bassez", "benchmark_analysis", "prepare_pseudobulk") == "ABSENT"
        assert tokens_for(scope, "Bassez", "benchmark_analysis", "pseudobulk") == (
            "hvg2000_bl", "CT_HR_hvg2000", "CT_HR_hvg500")
        assert tokens_for(scope, "Bassez", "benchmark_analysis", "mofa") == "ABSENT"
        files = selected_files("Adams", "pseudobulk", ("hvg2000_bl",))
        assert files == ["pseudobulks/Adams_pseudobulk_hvg2000_bl.rds", "results/Adams_Pseudobulk_hvg2000_bl.rds"]
        hr_files = selected_files("Bassez", "composition", AUTHOR_HR)
        assert len(hr_files) == 23 and not any("authors_LR" in path or "seuratres" in path for path in hr_files)
        assert selected_files("Bassez", "scpoli", ("hvg2000_highres_dims15",)) == [
            "embeddings/Bassez_hvg2000_highres_scpoli_dims15_embs.feather"]
        rows = normalized_rows(scope)
        assert rows.count(("Bassez", "benchmark_analysis", "composition")) == 1
        assert ("Adams", "benchmark_analysis", "prepare_pseudobulk") not in rows
        for bad, error in (
            (original + "Adams\tbenchmark_analysis\tpseudobulk\n", "overlapping whole-family"),
            ("Adams\tbenchmark_analysis\tpseudobulk\nAdams\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n", "whole-family/parameter overlap"),
            (original + "Adams\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n", "duplicate parameter"),
            ("Adams\tbenchmark_analysis\tpseudobulk\tCT_LR_hvg2000\n", "unknown parameter"),
            ("Adams\tbenchmark_analysis\tcomposition\tECODA_seuratres_2\n", "unknown parameter"),
            ("Adams\tbenchmark_analysis\tmofa\nAdams\tbenchmark_analysis\tpseudobulk\thvg2000\n", "cache overlap"),
            ("Adams\tbenchmark_analysis\tprepare_pseudobulk\nAdams\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n", "cache overlap")):
            request.write_text(bad)
            rejects(lambda: read_request(request, config), error)
        request.write_text("Adams\tbenchmark_analysis\tmofa\nAdams\tbenchmark_analysis\tpseudobulk\thvg2000_bl\n")
        assert tokens_for(read_request(request, config), "Adams", "benchmark_analysis", "prepare_pseudobulk") == (
            "hvg2000", "hvg1000", "hvg3000")
        pilotgm = tuple(f"hvg{n}_highres" for n in (1000, 2000, 3000))
        request.write_text("".join(f"Bassez\tbenchmark_analysis\tpilotgm\t{token}\n" for token in pilotgm))
        assert tokens_for(read_request(request, config), "Bassez", "benchmark_analysis", "pilotgm") == pilotgm
        assert selected_files("Bassez", "pilotgm", pilotgm) == [
            f"embeddings/Bassez_{token}_pilotgm_dists.feather" for token in pilotgm]
        step_gm = ("hvg1000_highres", "hvg3000_highres", "hvg2000_lowres")
        request.write_text("".join(f"Stephenson\tbenchmark_analysis\tpilotgm\t{token}\n" for token in step_gm))
        assert tokens_for(read_request(request, config), "Stephenson", "benchmark_analysis", "pilotgm") == step_gm
        assert selected_files("Stephenson", "pilotgm", step_gm) == [
            f"embeddings/Stephenson_{token}_pilotgm_dists.feather" for token in step_gm]
        for token in ("hvg1000_lowres", "hvg3000_lowres"):
            request.write_text(f"Bassez\tbenchmark_analysis\tpilotgm\t{token}\n")
            rejects(lambda: read_request(request, config), "unknown parameter token")
        request.write_text(original)

        config_path = base / "datasets.json"; config_path.write_text(json.dumps(config))
        nas = base / "nas"; nas.mkdir()
        source = base / "source.h5ad"; source.write_bytes(b"immutable input metadata fixture")
        md5 = hashlib.md5(source.read_bytes()).hexdigest(); size = source.stat().st_size
        proof = Path(f"{source}.md5"); proof.write_text(f"MD5={md5}\nSIZE={size}\nPATH={source}\n")
        refs = base / "inputs.tsv"
        refs.write_text("".join(f"{dataset}\tbenchmark_analysis\t{source}\t{md5}\t{size}\t{proof}\n" for dataset in config))
        run_root = base / "scratch/_ecoda_runs/request-scope"
        manifests = run_root / "manifests"; manifests.mkdir(parents=True)
        (manifests / "source.manifest").write_text("source identity\n")
        (manifests / "runtime.identity").write_text("runtime identity\n")
        exact = manifests / "requested_selection.tsv"; exact.write_text(original)
        selection = manifests / "selection.tsv"; selection.write_text("".join("\t".join(row) + "\n" for row in rows))
        binding = bind(run_root=run_root, release_id="test-release", nas_root=nas, stage="stage5",
            config_path=config_path, selection=selection, input_artifacts=refs)
        environment = dict(ECODA_RUN_ROOT=str(run_root), ECODA_ARTIFACT_BINDING=str(binding), NAS_TARGET_DIR=str(nas))
        with patch.dict(os.environ, environment):
            assert active_scope(config_path) == scope
            spec = importlib.util.spec_from_file_location("matrix_artifact_validator",
                ROOT / "src/5_run_benchmark_methods/matrix_artifact_validator.py")
            validator = importlib.util.module_from_spec(spec); spec.loader.exec_module(validator)
            with patch.dict(os.environ, ANALYSIS_ROOT=str(run_root / "artifacts/benchmark")):
                selected = validator.expected_artifacts(run_root / "artifacts/benchmark", "Bassez", "scpoli",
                    False, None, config_path=config_path)
                assert [path.name for path in selected] == ["Bassez_hvg2000_highres_scpoli_dims15_embs.feather"]
                whole = validator.expected_artifacts(run_root / "artifacts/benchmark", "Stephenson", "scpoli",
                    False, None, config_path=config_path)
                assert len(whole) == 8 and any("lowres" in path.name for path in whole)
                # The explicit LR token does not widen a whole-family GM request.
                gm_whole = validator.expected_artifacts(run_root / "artifacts/benchmark", "Stephenson", "pilotgm",
                    False, None, config_path=config_path)
                assert [path.name for path in gm_whole] == ["Stephenson_hvg2000_highres_pilotgm_dists.feather"]
            exact.write_text(original.replace("hvg2000_bl", "hvg500"))
            rejects(lambda: active_scope(config_path), "identity changed")
            exact.write_text(original)
            assert active_scope(config_path) == scope
    print("PASS: mixed full/partial scope, self-contained BL caches, HR-only outputs, overlap/tamper rejection")


if __name__ == "__main__":
    run()
