#!/usr/bin/env python3
"""Regression test for mixed annotation labels and numeric scores."""

import hashlib
import importlib.util
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src/4_cell_type_annotation/3.1_merge_annotations.py"
spec = importlib.util.spec_from_file_location("merge_annotations", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

obs = pd.DataFrame(
    {"layer2": ["T", None, 2], "S.Score": [0.1, np.nan, 0.3]},
    index=["c1", "c2", "c3"],
)
coerced = module._coerce_annotation_columns(obs, ["layer2", "S.Score"])
assert pd.api.types.is_string_dtype(coerced["layer2"])
assert pd.isna(coerced.loc["c2", "layer2"])
assert coerced.loc["c3", "layer2"] == "2"
assert pd.api.types.is_numeric_dtype(coerced["S.Score"])


def write_sidecar(path: Path) -> None:
    path.with_name(f"{path.name}.md5").write_text(
        f"MD5={hashlib.md5(path.read_bytes()).hexdigest()}\n"
        f"SIZE={path.stat().st_size}\nPATH={path}\n"
    )


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    input_path = root / "mixed_input.h5ad"
    output_path = root / "mixed_output.h5ad"
    annot_dir = root / "annotations"
    annot_dir.mkdir()
    input_obs = pd.DataFrame({"Sample": ["s1", "s1", "s1"]}, index=["c1", "c2", "c3"])
    ad.AnnData(X=np.ones((3, 1), dtype=np.float32), obs=input_obs).write_h5ad(input_path)
    write_sidecar(input_path)
    annotation_frame = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s1"],
            "cell_barcode": ["c1", "c2", "c3"],
            "layer1": ["T", "B", "T"],
            "layer2": ["T", None, 2],
            "layer3": ["T", "B", "T"],
            "layer_1": ["1", "1", "1"], "layer_2": ["2", "2", "2"],
            "layer_3": ["3", "3", "3"], "layer_4": ["4", "4", "4"],
            "layer_5": ["5", "5", "5"], "layer_6": ["6", "6", "6"],
            "scATOMIC_pred": ["T", "B", "T"],
            "classification_confidence": [0.9, 0.8, 0.7],
            "S.Score": [0.1, np.nan, 0.3],
            "G2M.Score": [0.2, 0.3, 0.4],
            "Phase": ["G1", "G2", "G1"],
        }
    )
    annotation_frame = module._coerce_annotation_columns(
        annotation_frame, module.ANNOT_OUTPUT_COLS
    )
    annotation_path = annot_dir / "annotations_chunk_1.feather"
    annotation_frame.to_feather(annotation_path)
    write_sidecar(annotation_path)
    module.merge_annotations(str(input_path), str(annot_dir), str(output_path))
    written = ad.read_h5ad(output_path)
    assert pd.api.types.is_string_dtype(written.obs["layer2"])
    assert written.obs.loc["c1", "layer2"] == "T"
    assert pd.isna(written.obs.loc["c2", "layer2"])
    assert written.obs.loc["c3", "layer2"] == "2"
    assert pd.api.types.is_numeric_dtype(written.obs["S.Score"])
print("annotation merge safety contracts OK")
