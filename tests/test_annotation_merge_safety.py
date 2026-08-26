#!/usr/bin/env python3
"""Regression test for mixed annotation labels in h5ad writes."""

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
spec.loader.exec_module(module)

obs = pd.DataFrame(
    {
        "layer2": ["T", None, 2],
        "S.Score": [0.1, np.nan, 0.3],
    },
    index=["c1", "c2", "c3"],
)
coerced = module._coerce_annotation_columns(obs, ["layer2", "S.Score"])
assert str(coerced["layer2"].dtype) == "string"
assert pd.isna(coerced.loc["c2", "layer2"])
assert coerced.loc["c3", "layer2"] == "2"
assert pd.api.types.is_numeric_dtype(coerced["S.Score"])

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "mixed_annotation.h5ad"
    ad.AnnData(X=np.ones((3, 1), dtype=np.float32), obs=coerced).write_h5ad(path)

print("annotation merge safety contracts OK")
