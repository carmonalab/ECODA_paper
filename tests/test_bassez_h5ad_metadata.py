#!/usr/bin/env python3
"""Tiny H5AD fixtures for source-bound, HR-only Bassez metadata derivatives."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
WORKER = ROOT / "src/2_dataset_specific_preprocessing/bassez_h5ad_metadata.py"
spec = importlib.util.spec_from_file_location("bassez_metadata", WORKER)
worker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(worker)


def expect_error(callback, fragment):
    try:
        callback()
    except (ValueError, OSError, RuntimeError) as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected failure containing {fragment!r}")


def fixture(path, *, invalid_fallback=False, complete=False):
    high = ["CD4", "", None, " UnKnOwN ", "  Myeloid  ", "NA", "<NA>", "none"]
    low = ["T", "B", "Cancer", "Fibro", "unknown", "Endo", "irrelevant", "Mast"]
    if invalid_fallback:
        low[2] = " NaN "
    if complete:
        high = [f"valid-{index}" for index in range(8)]
    data = ad.AnnData(
        X=np.arange(48, dtype=np.float32).reshape(8, 6),
        obs=pd.DataFrame({"Sample": ["S2", "S2", "S1", "S1", "S2", "S1", "S1", "S2"],
                          "sample": ["raw2", "raw2", "raw1", "raw1", "raw2", "raw1", "raw1", "raw2"],
                          "cellSubType": pd.Categorical(high), "cellType": pd.Categorical(low),
                          "expansion": ["post", "post", "pre", "pre", "post", "pre", "pre", "post"]},
                         index=["cell7", "cell2", "cell8", "cell1", "cell6", "cell3", "cell5", "cell4"]),
        var=pd.DataFrame({"hvg_rank": np.arange(6, dtype=float)}, index=[f"gene{i}" for i in range(6)]),
    )
    data.layers["counts"] = sp.csr_matrix(np.arange(48, dtype=np.int32).reshape(8, 6))
    data.raw = data.copy()
    data.obsm["X_pca"] = np.arange(24, dtype=np.float64).reshape(8, 3)
    data.obsm["X_pca_harmony"] = data.obsm["X_pca"] * 0.5
    data.uns["cohort_selection"] = {"upstream_identity": "preserve-exactly"}
    data.uns["batch_contract"] = {"fingerprint": "original-contract", "keys": ["batch"]}
    data.uns["optional_none"] = None
    data.uns["structured_metadata"] = np.asarray(
        [("geneA", 1.25), ("geneB", np.nan)], dtype=[("name", object), ("score", np.float64)]
    )
    data.write_h5ad(path, compression="gzip")
    with h5py.File(path, "r+") as handle:
        handle.attrs["fixture_attribute"] = np.asarray([1, 3, 5])
        handle["layers/counts/data"].attrs["array_attribute"] = "keep this"


def source_arguments(source, scratch, run_root):
    md5 = hashlib.md5(source.read_bytes()).hexdigest()
    key = hashlib.sha256(str(source).encode()).hexdigest()[:32]
    record = scratch / "_ecoda_runs/parent-stage3/manifests/artifacts" / f"{key}.record"
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(f"PATH={source}\nSIZE={source.stat().st_size}\nMD5={md5}\nRUN_ID=parent-stage3\nPRODUCER=stage3\nSTATE=PUBLISHED\n")
    sidecar = Path(f"{source}.md5")
    sidecar.write_text(f"MD5={md5}\nSIZE={source.stat().st_size}\nPATH={source}\n")
    return dict(input_file=source, expected_md5=md5, parent_record=record,
                parent_run_id="parent-stage3", scratch_root=scratch, run_root=run_root,
                output_file=run_root / "artifacts/datasets/Bassez/output" / source.name)


def main():
    assert all(worker._missing(value) for value in (None, "", " NA\t", "NaN", "nOnE", " Unknown\r\n"))
    # Match R trimws()/the exact listed sentinels: do not broaden the policy.
    assert not worker._missing("<NA>")
    assert not worker._missing("\u00a0NA\u00a0")
    assert not worker._missing("NA-like")
    for selection in worker._slices((7, 19, 23), 4, buffer_bytes=64):
        assert np.prod([part.stop - part.start for part in selection]) * 4 <= 64

    with tempfile.TemporaryDirectory(prefix="ecoda-bassez-metadata-") as temporary:
        scratch = Path(temporary).resolve() / "scratch"
        scratch.mkdir()
        run_root = scratch / "_ecoda_runs/bassez-repair"
        manifests = run_root / "manifests"
        manifests.mkdir(parents=True)
        (manifests / "source.manifest").write_text("fixture pinned source identity\n")
        (manifests / "runtime.identity").write_text("fixture pinned runtime identity\n")
        source = scratch / "Bassez.h5ad"
        fixture(source)
        args = source_arguments(source, scratch, run_root)
        original = source.read_bytes()
        original_sidecar = Path(f"{source}.md5").read_bytes()
        original_record = args["parent_record"].read_bytes()
        expected, before = worker.annotations(source)
        assert expected == ["CD4", "B", "Cancer", "Fibro", "  Myeloid  ", "Endo", "<NA>", "Mast"]
        assert before["n_replaced"] == 5 and before["before_valid_categories"] == 3

        report = worker.create_derivative(**args)
        output = args["output_file"]
        assert report["identity"]["coverage"]["n_cells"] == 8
        assert report["identity"]["coverage"]["n_samples"] == 2
        assert report["identity"]["coverage"]["after_categories"] == 8
        assert worker.annotations(output)[1]["n_replaced"] == 0
        assert worker.protected_tree(source) == worker.protected_tree(output)
        assert worker.protected_tree(source, buffer_bytes=64) == worker.protected_tree(output, buffer_bytes=64)
        assert worker.validate_derivative(**args) == report
        reopened = ad.read_h5ad(output)
        assert reopened.obs.cellSubType.astype(str).tolist() == expected
        assert reopened.obs.index.tolist() == ["cell7", "cell2", "cell8", "cell1", "cell6", "cell3", "cell5", "cell4"]
        assert np.array_equal(reopened.layers["counts"].toarray(), np.arange(48, dtype=np.int32).reshape(8, 6))
        assert np.array_equal(reopened.raw.X, np.arange(48, dtype=np.float32).reshape(8, 6))
        with h5py.File(output, "r") as handle:
            assert handle["uns/cohort_selection/upstream_identity"].asstr()[()] == "preserve-exactly"
            assert handle["uns/batch_contract/fingerprint"].asstr()[()] == "original-contract"
            assert np.array_equal(handle["X"][:], np.arange(48, dtype=np.float32).reshape(8, 6))
            assert np.array_equal(handle["obsm/X_pca"][:], np.arange(24, dtype=np.float64).reshape(8, 3))
        assert source.read_bytes() == original
        assert Path(f"{source}.md5").read_bytes() == original_sidecar
        assert args["parent_record"].read_bytes() == original_record
        changed_attribute = scratch / "changed-attribute.h5ad"
        shutil.copyfile(source, changed_attribute)
        with h5py.File(changed_attribute, "r+") as handle:
            handle["layers/counts/data"].attrs["array_attribute"] = "changed"
        assert worker.protected_tree(changed_attribute) != worker.protected_tree(source)
        ragged = scratch / "numeric-ragged.h5ad"
        shutil.copyfile(source, ragged)
        with h5py.File(ragged, "r+") as handle:
            values = handle["uns"].create_dataset("numeric_ragged", shape=(1,), dtype=h5py.vlen_dtype(np.int32))
            values[0] = np.arange(3, dtype=np.int32)
        expect_error(lambda: worker.protected_tree(ragged), "variable-length numeric")
        old_outputs = {path: path.read_bytes() for path in (output, Path(f"{output}.md5"), Path(f"{output}.repair.json"))}
        expect_error(lambda: worker.create_derivative(**args), "existing repair outputs are preserved")
        assert all(path.read_bytes() == value for path, value in old_outputs.items())

        command = [sys.executable, str(WORKER), "audit", "--input-file", str(source),
                   "--expected-md5", args["expected_md5"], "--parent-record", str(args["parent_record"]),
                   "--parent-run-id", args["parent_run_id"], "--scratch-root", str(scratch)]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        assert json.loads(result.stdout)["n_replaced"] == 5
        incomplete = subprocess.run([*command, "--require-complete"], capture_output=True, text=True)
        assert incomplete.returncode != 0 and "5 missing HR labels" in incomplete.stderr

        legacy = scratch / "legacy.h5ad"
        fixture(legacy)
        legacy_md5 = hashlib.md5(legacy.read_bytes()).hexdigest()
        legacy_proof = Path(f"{legacy}.md5")
        legacy_proof.write_text(f"MD5={legacy_md5}\nSIZE={legacy.stat().st_size}\nPATH={legacy}\n")
        legacy_bytes, legacy_proof_bytes = legacy.read_bytes(), legacy_proof.read_bytes()
        legacy_args = {"input_file": legacy, "expected_md5": legacy_md5,
                       "legacy_input_proof": legacy_proof, "scratch_root": scratch,
                       "run_root": run_root, "output_file": output.with_name("legacy.h5ad")}
        legacy_report = worker.create_derivative(**legacy_args)
        assert worker.validate_derivative(**legacy_args) == legacy_report
        parent = legacy_report["identity"]["parent"]
        assert parent["provenance_kind"] == "legacy_validated_input"
        assert "producer" not in parent and "producer_run_id" not in parent
        assert parent["checksum_proof_sha256"] == worker._digest(legacy_proof)
        assert legacy.read_bytes() == legacy_bytes and legacy_proof.read_bytes() == legacy_proof_bytes
        legacy_command = [sys.executable, str(WORKER), "create", "--input-file", str(legacy),
                          "--expected-md5", legacy_md5, "--legacy-input-proof", str(legacy_proof),
                          "--scratch-root", str(scratch), "--run-root", str(run_root),
                          "--output-file", str(output.with_name("legacy-cli.h5ad"))]
        created = subprocess.run(legacy_command, check=True, capture_output=True, text=True)
        assert json.loads(created.stdout)["identity"]["parent"]["provenance_kind"] == "legacy_validated_input"
        validated_command = [*legacy_command]
        validated_command[2] = "validate"
        subprocess.run(validated_command, check=True, capture_output=True, text=True)
        expect_error(lambda: worker.create_derivative(**{**legacy_args, "output_file": scratch / "unsafe.h5ad"}), "run-owned")
        wrong_proof = scratch / "wrong-proof.md5"
        wrong_proof.write_text(f"MD5={legacy_md5}\nSIZE={legacy.stat().st_size}\nPATH={source}\n")
        wrong_args = {**legacy_args, "legacy_input_proof": wrong_proof, "output_file": output.with_name("wrong-proof.h5ad")}
        expect_error(lambda: worker.create_derivative(**wrong_args), "provenance differs")
        wrong_proof.write_text(f"SIZE={legacy.stat().st_size}\nMD5={legacy_md5}\nPATH={legacy}\n")
        expect_error(lambda: worker.create_derivative(**wrong_args), "invalid MD5")
        no_fallback = {**args, "legacy_input_proof": Path(f"{source}.md5"), "output_file": output.with_name("no-fallback.h5ad")}
        expect_error(lambda: worker.create_derivative(**no_fallback), "exactly one")
        noncanonical = scratch / "unowned-producer.record"
        noncanonical.write_bytes(original_record)
        expect_error(lambda: worker.create_derivative(**{**args, "parent_record": noncanonical, "output_file": output.with_name("unowned-record.h5ad")}), "canonical run-owned")
        bad_record = scratch / "_ecoda_runs/invalid-parent/manifests/artifacts" / args["parent_record"].name
        bad_record.parent.mkdir(parents=True)
        bad_record.write_text(original_record.decode().replace("RUN_ID=parent-stage3", "RUN_ID=invalid-parent").replace("STATE=PUBLISHED", "STATE=FAILED"))
        expect_error(lambda: worker.create_derivative(**{**args, "parent_record": bad_record, "parent_run_id": "invalid-parent", "output_file": output.with_name("bad-record.h5ad")}), "provenance differs")
        assert not wrong_args["output_file"].exists()
        assert not no_fallback["output_file"].exists()

        invalid = scratch / "invalid.h5ad"
        fixture(invalid, invalid_fallback=True)
        invalid_args = source_arguments(invalid, scratch, run_root)
        expect_error(lambda: worker.create_derivative(**invalid_args), "invalid cellType fallback")
        assert not invalid_args["output_file"].exists()

        complete = scratch / "complete.h5ad"
        fixture(complete, complete=True)
        complete_args = source_arguments(complete, scratch, run_root)
        expect_error(lambda: worker.create_derivative(**complete_args), "already complete")
        assert not complete_args["output_file"].exists()

        # Failure removes only the newly allocated temporary files.
        failed_args = {**args, "output_file": output.with_name("failure.h5ad")}
        original_writer = worker._write_repaired_column

        def fail_after_edit(handle, filled):
            original_writer(handle, filled)
            raise RuntimeError("injected repair failure")

        worker._write_repaired_column = fail_after_edit
        expect_error(lambda: worker.create_derivative(**failed_args), "injected repair failure")
        worker._write_repaired_column = original_writer
        assert not failed_args["output_file"].exists()
        assert not list(output.parent.glob(".failure.h5ad.repair-*"))
        assert all(path.read_bytes() == value for path, value in old_outputs.items())
        assert source.read_bytes() == original

        # Even self-consistent replacement sidecars cannot hide changed data.
        with h5py.File(output, "r+") as handle:
            handle["layers/counts/data"][0] += 1
        report["output"].update(md5=worker._digest(output, "md5"), size=output.stat().st_size)
        Path(f"{output}.repair.json").write_text(worker._json(report) + "\n")
        Path(f"{output}.md5").write_text(f"MD5={report['output']['md5']}\nSIZE={output.stat().st_size}\nPATH={output}\n")
        expect_error(lambda: worker.validate_derivative(**args), "protected H5AD datasets/attributes differ")
    print("PASS: exact Bassez fallback, all other H5AD data/attributes unchanged, source proof, create-only publication, bounded hashes, cleanup")


if __name__ == "__main__":
    main()
