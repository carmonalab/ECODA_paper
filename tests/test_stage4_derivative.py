"""Real tiny Stage 4 copy, preparation and merge; no annotator or scheduler."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind, stage4_input_copy
from src.utils.py.annotation_contract import validate_h5ad


def checksum(path):
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    Path(f"{path}.md5").write_text(f"MD5={md5}\nSIZE={path.stat().st_size}\nPATH={path}\n")
    return md5


def run():
    with tempfile.TemporaryDirectory(prefix="ecoda-stage4-derivative-") as temporary:
        base = Path(temporary).resolve()
        scratch = base / "scratch"
        nas = base / "nas"; nas.mkdir()
        run_root = scratch / "_ecoda_runs" / "stage4-fixture"
        manifests = run_root / "manifests"; manifests.mkdir(parents=True)
        (run_root / "metadata").write_text("STAGE=stage4\nRUN_ID=stage4-fixture\n")
        for name in ("source.manifest", "runtime.identity"):
            (manifests / name).write_text(f"fixture {name}\n")
        config = json.loads((ROOT / "datasets.json").read_text())
        name = config["Stephenson"]["views"]["benchmark_analysis"]["output_file_name"]
        source = scratch / "_ecoda_runs/stage3-parent/artifacts/datasets/Stephenson/output" / name
        source.parent.mkdir(parents=True)
        counts = sp.csr_matrix(np.arange(24, dtype=np.int32).reshape(6, 4))
        obs = pd.DataFrame({"Sample": ["s1"] * 3 + ["s2"] * 3,
                            "Status": ["Healthy"] * 3 + ["Covid"] * 3,
                            "initial_clustering": ["author LR"] * 6,
                            "full_clustering": ["author HR"] * 6,
                            "layer2": ["legacy label"] * 6}, index=[f"c{x}" for x in range(6)])
        data = ad.AnnData(X=counts.astype(np.float32) / 10, obs=obs,
                          var=pd.DataFrame({"hvg_rank": np.arange(4.)}, index=[f"g{x}" for x in range(4)]))
        data.layers["counts"] = counts.copy(); data.layers["other"] = counts * 2
        data.raw = ad.AnnData(X=counts.copy(), obs=obs.copy(), var=data.var.copy())
        data.obsm["X_pca"] = np.arange(12.).reshape(6, 2)
        data.obsp["connectivities"] = sp.eye(6, format="csr")
        data.uns["upstream"] = {"identity": "preserved"}
        data.write_h5ad(source)
        parent_md5 = checksum(source)
        selection = manifests / "selection.tsv"
        selection.write_text("Stephenson\tbenchmark_analysis\n")
        refs = base / "inputs.tsv"
        refs.write_text(f"Stephenson\tbenchmark_analysis\t{source}\t{parent_md5}\t{source.stat().st_size}\t{source}.md5\n")
        binding = bind(run_root=run_root, release_id="test-release", nas_root=nas, stage="stage4",
                       config_path=ROOT / "datasets.json", selection=selection, input_artifacts=refs,
                       producer_run="stage3-parent")
        owner_dir = scratch / "_ecoda_owners/stage4/derivative_stage4-fixture_Stephenson"
        owner_dir.mkdir(parents=True)
        (owner_dir / "owner").write_text("RUN_ID=stage4-fixture\nSTAGE=stage4\nKEY=derivative/stage4-fixture/Stephenson\nSTATE=ACTIVE\n")
        environment = dict(PROJECT_ROOT=str(ROOT), DATASETS_JSON_FILE=str(ROOT / "datasets.json"),
                           HPC_SCRATCH_DIR=str(scratch), NAS_TARGET_DIR=str(nas), ECODA_ARTIFACT_BINDING=str(binding),
                           ECODA_RUN_ROOT=str(run_root), ECODA_RUN_ID="stage4-fixture", ANNOTATION_RUN_ID="stage4-fixture",
                           DS_NAME="Stephenson", NUMBA_CACHE_DIR=str(base / "numba"), MPLCONFIGDIR=str(base / "mpl"))
        with patch.dict(os.environ, environment):
            nas.rename(base / "nas.hidden")
            output = stage4_input_copy(binding, ROOT / "datasets.json", "Stephenson", "benchmark_analysis")
            assert output != source and output.read_bytes() == source.read_bytes()
            assert stage4_input_copy(binding, ROOT / "datasets.json", "Stephenson", "benchmark_analysis") == output
            subprocess.run([sys.executable, str(ROOT / "src/4_cell_type_annotation/1.1_prepare_chunks.py"),
                            "--views", "benchmark_analysis", "--run-root", str(run_root)], check=True)
            union = run_root / "datasets/Stephenson/union/union.h5ad"
            chunks = list((run_root / "datasets/Stephenson/chunks").glob("chunk_*.txt"))
            assert len(chunks) == 1 and chunks[0].read_text().splitlines() == [str(union), "s1", "s2"]
            assert hashlib.md5(source.read_bytes()).hexdigest() == parent_md5
            # A synthetic complete dual-method result stands in only for fitting.
            frame = pd.DataFrame({"Sample": obs.Sample.values, "cell_barcode": obs.index.values})
            for column in ("layer1", "layer2", "layer3", *(f"layer_{n}" for n in range(1, 7)), "scATOMIC_pred"):
                frame[column] = ["T", "B"] * 3
            frame["classification_confidence"] = "confident"
            frame["S.Score"] = 0.1; frame["G2M.Score"] = 0.2; frame["Phase"] = "G1"
            annot_dir = run_root / "datasets/Stephenson/annotations"; annot_dir.mkdir(exist_ok=True)
            feather = annot_dir / "annotations_chunk_1.feather"; frame.to_feather(feather); checksum(feather)
            spec = importlib.util.spec_from_file_location("stage4_merge", ROOT / "src/4_cell_type_annotation/3.1_merge_annotations.py")
            merge = importlib.util.module_from_spec(spec); spec.loader.exec_module(merge)
            merge.merge_annotations(str(output), str(annot_dir), str(output), union_path=str(union))
            validate_h5ad(output)
            actual = ad.read_h5ad(output)
            assert actual.obs_names.tolist() == data.obs_names.tolist()
            assert actual.var_names.tolist() == data.var_names.tolist()
            assert np.array_equal(actual.var["hvg_rank"].values, data.var["hvg_rank"].values)
            for column in ("Sample", "Status", "initial_clustering", "full_clustering"):
                assert actual.obs[column].tolist() == data.obs[column].tolist()
            assert np.array_equal(actual.X.toarray(), data.X.toarray())
            assert np.array_equal(actual.raw.X.toarray(), data.raw.X.toarray())
            for layer in data.layers:
                assert np.array_equal(actual.layers[layer].toarray(), data.layers[layer].toarray())
            assert np.array_equal(actual.obsm["X_pca"], data.obsm["X_pca"])
            assert np.array_equal(actual.obsp["connectivities"].toarray(), data.obsp["connectivities"].toarray())
            assert actual.uns["upstream"]["identity"] == "preserved"
            assert hashlib.md5(source.read_bytes()).hexdigest() == parent_md5
            final_md5 = hashlib.md5(output.read_bytes()).hexdigest()
            try:
                stage4_input_copy(binding, ROOT / "datasets.json", "Stephenson", "benchmark_analysis")
            except ValueError:
                pass
            else:
                raise AssertionError("annotated output was replaced by its original input copy")
            assert hashlib.md5(output.read_bytes()).hexdigest() == final_md5
    print("PASS: create-only Stage4 copy, hidden-NAS preparation, complete dual merge and immutable parent/arrays")


if __name__ == "__main__":
    run()
