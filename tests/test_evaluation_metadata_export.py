"""Real obs-only exports for Kidney's two axes and scalar-label compatibility."""
from argparse import Namespace
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.export_h5ad_sample_metadata import export, write_metadata
from src.utils.py.datasets_io import read_datasets_json


def rejects(call, message):
    try:
        call()
    except ValueError as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError(f"expected failure: {message}")


def write_h5ad(path, obs):
    data = ad.AnnData(X=np.arange(len(obs) * 2, dtype=np.float32).reshape(-1, 2),
                     obs=obs.copy(), var=pd.DataFrame(index=["g1", "g2"]))
    data.layers["counts"] = data.X.astype(np.int32)
    data.write_h5ad(path)


def run():
    with tempfile.TemporaryDirectory(prefix="ecoda-evaluation-metadata-") as temporary:
        base = Path(temporary).resolve()
        registry = json.loads((ROOT / "datasets.json").read_text())
        config = {name: registry[name] for name in ("Kidney_KPMP_full", "Adams")}
        config_path = base / "datasets.json"; config_path.write_text(json.dumps(config))
        tissues = ["kidney"] * 45 + ["cortex of kidney"] * 23 + ["renal medulla"] * 12
        obs = pd.DataFrame([{
            "Sample": f"S{sample:02d}", "specimen": f"source{sample:02d}",
            "condition.l1": "control" if sample < 40 else "CKD",
            "tissue": tissue, "suspension_type": "cell" if tissue == "kidney" else "nucleus",
            "subclass.l1": "Immune", "subclass.l3": "T cell" if cell == 0 else "B cell",
        } for sample, tissue in enumerate(tissues) for cell in range(2)])
        obs.index = [f"cell{i}" for i in range(len(obs))]
        source = base / "kidney.h5ad"; write_h5ad(source, obs)
        source_md5 = hashlib.md5(source.read_bytes()).hexdigest()
        args = Namespace(config=config_path, dataset="Kidney_KPMP_full",
            view="batch_effect_uncorrected", input_file=source, output=base / "kidney_uncorrected.feather",
            chunk_size=3, analysis_variant=None, check=False)
        environment = {key: value for key, value in os.environ.items()
                       if key not in {"ECODA_ARTIFACT_BINDING", "ANALYSIS_VARIANT", "ANALYSIS_ROOT", "ANALYSIS_NAS_ROOT"}}
        with patch.dict(os.environ, environment, clear=True):
            for view in ("batch_effect_uncorrected", "batch_effect_corrected"):
                args.view = view; args.output = base / f"{view}.feather"
                export(args)
                frame = pd.read_feather(args.output)
                assert len(frame) == 80 and frame["Sample"].tolist() == [f"S{i:02d}" for i in range(80)]
                assert frame["tissue"].value_counts().to_dict() == {"kidney": 45, "cortex of kidney": 23, "renal medulla": 12}
                assert frame.columns[:3].tolist() == ["Sample", "condition.l1", "tissue"]
                assert read_datasets_json(config_path)["Kidney_KPMP_full"]["label_col"] == "condition.l1"
                args.check = True; export(args); args.check = False
            original_args = copy.copy(args)
            for name, frame, message in (
                ("missing", obs.drop(columns="tissue"), "missing configured obs columns"),
                ("mixed", obs.assign(tissue=obs.tissue.mask(obs.index == "cell1", "renal medulla")), "varies within Sample"),
                ("blank", obs.assign(tissue=obs.tissue.mask(obs.index == "cell1", "")), "missing/blank"),
                ("mixed_condition", obs.assign(**{"condition.l1": obs["condition.l1"].mask(obs.index == "cell1", "CKD")}), "varies within Sample"),
            ):
                path = base / f"{name}.h5ad"; write_h5ad(path, frame)
                args.input_file = path; args.output = base / f"{name}.feather"
                rejects(lambda: export(args), message)
                assert not args.output.exists()
            args = copy.copy(original_args)
            invalid = copy.deepcopy(config)
            invalid["Kidney_KPMP_full"]["columns"]["evaluation_labels"]["tissue"] = "unknown_column"
            invalid_path = base / "unknown.json"; invalid_path.write_text(json.dumps(invalid))
            args.config = invalid_path; args.output = base / "unknown.feather"
            rejects(lambda: export(args), "missing configured obs columns")
            invalid["Kidney_KPMP_full"]["columns"]["evaluation_labels"]["tissue"] = ""
            invalid_path.write_text(json.dumps(invalid))
            rejects(lambda: export(args), "Invalid dataset-level evaluation_labels")
            invalid["Kidney_KPMP_full"]["columns"]["evaluation_labels"]["tissue"] = "suspension_type"
            invalid_path.write_text(json.dumps(invalid))
            rejects(lambda: export(args), "cannot be technical batch keys")

            # A checksum-valid old export without the newly declared axis is not reusable.
            args = copy.copy(original_args); args.check = True
            frame = pd.read_feather(args.output)
            write_metadata(frame.drop(columns="tissue"), args.output)
            rejects(lambda: export(args), "missing configured columns")
            frame.loc[0, "tissue"] = None
            write_metadata(frame, args.output)
            rejects(lambda: export(args), "missing/blank")

            # An ordinary dataset keeps its scalar label and existing column order.
            ordinary_obs = pd.DataFrame({"Sample": ["A", "A", "B", "B"],
                "Disease_Identity": ["healthy", "disease", "disease", "disease"],
                "OriginalAnnotationLevel1": ["T", "B", "T", "B"],
                "OriginalAnnotationLevel2": ["T1", "B1", "T1", "B1"]},
                index=["a", "b", "c", "d"])
            ordinary = base / "adams.h5ad"; write_h5ad(ordinary, ordinary_obs)
            args = Namespace(config=config_path, dataset="Adams", view="benchmark_analysis",
                input_file=ordinary, output=base / "adams.feather", chunk_size=1, analysis_variant=None, check=False)
            export(args)
            frame = pd.read_feather(args.output)
            assert frame.columns.tolist() == ["Sample", "Disease_Identity", "OriginalAnnotationLevel1", "OriginalAnnotationLevel2"]
            assert frame.Sample.tolist() == ["A", "B"]
            assert frame.Disease_Identity.tolist() == ["healthy", "disease"]
        assert hashlib.md5(source.read_bytes()).hexdigest() == source_md5
    print("PASS: two Kidney evaluation axes/all80 samples, stable biological labels, missing/malformed-axis rejection, scalar-label compatibility, and unchanged H5AD")


if __name__ == "__main__":
    run()
