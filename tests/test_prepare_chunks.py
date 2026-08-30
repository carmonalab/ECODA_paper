#!/usr/bin/env python3
"""Exercise run-owned multi-view annotation chunk preparation."""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "ecoda_prepare_chunks", ROOT / "src/4_cell_type_annotation/1.1_prepare_chunks.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-chunks-") as raw:
        root = Path(raw)
        scratch_target = root / "scratch_target"
        scratch_target.mkdir()
        scratch = root / "scratch"
        scratch.symlink_to(scratch_target, target_is_directory=True)
        output = scratch / "Synthetic" / "output"
        output.mkdir(parents=True)
        config = {
            "Synthetic": {
                "views": {
                    "view_a": {"input_file_name": "raw_a.h5ad", "output_file_name": "view_a.h5ad"},
                    "view_b": {"input_file_name": "raw_b.h5ad", "output_file_name": "view_b.h5ad"},
                }
            }
        }
        config_path = root / "datasets.json"
        config_path.write_text(json.dumps(config))
        for name in ("view_a.h5ad", "view_b.h5ad"):
            obs = pd.DataFrame({"Sample": ["s1", "s1", "s2", "s2"]}, index=["c1", "c2", "c3", "c4"])
            var = pd.DataFrame(index=["g1", "g2"])
            data = sparse.csr_matrix(np.ones((4, 2), dtype="int32"))
            obj = ad.AnnData(X=data.astype("float32"), obs=obs, var=var)
            obj.layers["counts"] = data
            obj.write_h5ad(output / name)
        run_root = scratch / "_ecoda_runs" / "run"
        run_root.mkdir(parents=True)
        (run_root / "metadata").write_text("STAGE=stage4\nRUN_ID=run\nSTATE=ACTIVE\n")
        owner_dir = scratch / "_ecoda_owners" / "stage4" / "Synthetic"
        owner_dir.mkdir(parents=True)
        (owner_dir / "owner").write_text(
            "RUN_ID=run\nSTATE=ACTIVE\nSTAGE=stage4\nKEY=Synthetic\n"
        )
        os.environ.update(
            {
                "PROJECT_ROOT": str(ROOT),
                "HPC_SCRATCH_DIR": str(scratch),
                "DATASETS_JSON_FILE": str(config_path),
                "DS_NAME": "Synthetic",
                "SAMPLE_COLNAME": "Sample",
                "ANNOTATION_RUN_ID": "run",
                "ANNOTATION_RUN_ROOT": str(run_root),
            }
        )
        sys.argv = ["1.1_prepare_chunks.py", "--views", "view_a,view_b", "--run-root", str(run_root)]
        MODULE.main()
        records = json.loads(
            (run_root / "datasets" / "Synthetic" / "source_artifacts.json").read_text()
        )
        assert records[0]["path"] == os.path.abspath(os.fspath(output / "view_a.h5ad"))
        assert str(scratch_target) not in records[0]["path"]
        sys.argv = [
            "1.1_prepare_chunks.py",
            "--views",
            "view_a,view_b",
            "--run-root",
            str(run_root),
            "--validate-only",
        ]
        MODULE.main()
        chunk_dir = run_root / "datasets" / "Synthetic" / "chunks"
        union = run_root / "datasets" / "Synthetic" / "union" / "union.h5ad"
        chunks = sorted(chunk_dir.glob("chunk_*.txt"))
        assert chunks and union.is_file()
        assert chunks[0].read_text().splitlines()[0] == str(union.resolve())
        assert not (output / "chunks").exists()
    print("run-owned chunk preparation: OK")


if __name__ == "__main__":
    main()
