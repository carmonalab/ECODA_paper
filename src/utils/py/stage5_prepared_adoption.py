"""Copy reviewed preparation artifacts into exact fresh derivative namespaces.

The original contract requires all 23 payloads. A separate metadata-only
contract copies precisely the selected batch dataset/view union, without fits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

import derivative_artifacts as da

SCHEMA = "ecoda_stage5_prepared_adoption_v1"
METADATA_SCHEMA = "ecoda_stage5_metadata_adoption_v1"
CONTEXTS = ("benchmark", "batch_effect/uncorrected_final", "batch_effect/corrected_final")
VARIANTS = ("schvg2000", "hvg2000", "hvg500", "hvg2000_bl", "hvg1000", "hvg3000")


def require(condition, message):
    if not condition:
        raise ValueError("prepared adoption: " + message)


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with da.absolute(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fields(path):
    result = {}
    for line in da.absolute(path).read_text().splitlines():
        key, sep, value = line.partition("=")
        require(sep and key not in result, "malformed record " + str(path))
        result[key] = value
    return result


def checked(proof):
    require(set(proof) == {"path", "sha256"} and re.fullmatch("[0-9a-f]{64}", proof["sha256"]), "invalid proof")
    path = da.absolute(proof["path"])
    require(digest(path) == proof["sha256"], "proof changed " + str(path))
    return path


def input_rows(binding):
    result = {}
    for line in Path(binding["input_artifacts"]).read_text().splitlines():
        dataset, view, path, md5, size, _ = line.split("\t")
        key = (dataset, view)
        require(key not in result, "expected one processed H5AD per dataset/view")
        result[key] = (path, md5, size)
    return result


def owner_path(source, scratch):
    hashed = hashlib.sha256(str(source).encode()).hexdigest()[:32]
    name = re.sub(r"[/:,\t |]", "_", source.name) + "_" + hashed
    if len(name.encode()) > 100:
        name = "artifact_" + hashed
    return scratch / "_ecoda_owners/artifact" / name / "owner"


def prepare(receipt_path, sha256, config_path, context):
    receipt_path = da.absolute(receipt_path)
    require(re.fullmatch("[0-9a-f]{64}", sha256) and digest(receipt_path) == sha256, "receipt checksum differs")
    receipt = json.loads(receipt_path.read_text())
    require(receipt["schema"] in (SCHEMA, METADATA_SCHEMA) and receipt["state"] == "APPROVED_COPY_ONLY"
            and receipt["original_gate_state"] == "FAILED", "receipt purpose/state differs")
    metadata_only = receipt["schema"] == METADATA_SCHEMA
    config_path = da.absolute(config_path)
    require(receipt["config_sha256"] == digest(config_path), "configuration differs")
    current, _ = da.active_binding(config_path)
    contexts = set(receipt["target_run_ids"])
    require(context in contexts and current["run_id"] == receipt["target_run_ids"][context], "target run differs")
    require(set(receipt["sources"]) == contexts and
            (bool(contexts) and contexts <= set(CONTEXTS[1:]) if metadata_only else contexts == set(CONTEXTS)),
            "context coverage differs")
    require(len(set(receipt["target_run_ids"].values())) == len(contexts), "target runs must be distinct")
    require(len(receipt["evidence"]) == 3 and len({item["path"] for item in receipt["evidence"]}) == 3,
            "missing/duplicate reviewed scientific/accounting evidence")
    for evidence in receipt["evidence"]:
        checked(evidence)
    config = json.loads(config_path.read_text())
    batch = {key for key, value in config.items() if not key.startswith("_") and value.get("use_for_batch_effect") is True}
    require(len(batch) == 8, "this finite recovery requires exactly eight configured batch cohorts")
    if metadata_only:
        require(set(receipt["selected_metadata"]) == contexts, "metadata selection contexts differ")
        expected = set()
        for ctx, rows in receipt["selected_metadata"].items():
            view = "batch_effect_" + ("uncorrected" if ctx.endswith("uncorrected_final") else "corrected")
            require(isinstance(rows, list) and rows and all(set(row) == {"dataset", "view"} and
                    row["dataset"] in batch and row["view"] == view and
                    view in config[row["dataset"]].get("views", {}) for row in rows), "invalid selected batch metadata")
            keys = {(ctx, "metadata", row["dataset"], row["view"], "") for row in rows}
            require(len(keys) == len(rows), "duplicate selected metadata")
            expected.update(keys)
        selected_scope = {(row["dataset"], row["view"]) for row in receipt["selected_metadata"][context]}
        binding_scope = {tuple(line.split("\t")[:2]) for line in Path(current["selection"]).read_text().splitlines()}
        require(binding_scope == selected_scope, "metadata scope differs from selected target dataset/view union")
    else:
        expected = {("benchmark", "pseudobulk", "Stephenson", "benchmark_analysis", variant) for variant in VARIANTS}
        expected.add(("benchmark", "timing", "Stephenson", "benchmark_analysis", ""))
        for ctx in CONTEXTS[1:]:
            view = "batch_effect_" + ("uncorrected" if ctx.endswith("uncorrected_final") else "corrected")
            expected.update((ctx, "metadata", dataset, view, "") for dataset in batch)
    scratch = Path(current["run_root"]).parent.parent
    current_inputs = input_rows(current)
    current_runtime = fields(Path(current["run_root"]) / "manifests/runtime.identity")
    old = {}
    metadata_order = {}
    for ctx, source in receipt["sources"].items():
        old_root = da.absolute(source["run_root"])
        require(old_root.parent == scratch / "_ecoda_runs" and old_root.name != receipt["target_run_ids"][ctx], "source run root differs")
        binding_path = checked(source["binding"])
        old_config = checked(source["config"])
        require(binding_path == old_root / "manifests/artifact_bindings.json" and digest(old_config) == digest(config_path), "source binding/config differs")
        binding, _ = da.load(binding_path, old_config)
        require(binding["run_root"] == str(old_root) and binding["stage"] == "stage5", "source binding stage differs")
        for key, name in (("source_manifest", "source.manifest"), ("runtime_identity", "runtime.identity")):
            require(checked(source[key]) == old_root / "manifests" / name, "source control path differs")
        runtime = fields(old_root / "manifests/runtime.identity")
        require(all(re.fullmatch("[0-9a-f]{64}", runtime.get(key, "")) and runtime[key] == current_runtime.get(key)
                    for key in ("RUNTIME_IMAGE_SHA256", "RUNTIME_MANIFEST_SHA256")), "runtime image/manifest differs")
        old[ctx] = (old_root, input_rows(binding))
        metadata_order[ctx] = list(dict.fromkeys(tuple(line.split("\t")[:2])
            for line in Path(binding["selection"]).read_text().splitlines()))
    selected = []
    observed = set()
    for item in receipt["artifacts"]:
        ctx, kind, dataset, view = (item[key] for key in ("context", "kind", "dataset", "view"))
        source = da.absolute(item["source"])
        old_root, old_inputs = old[ctx]
        variant = ""
        if kind == "pseudobulk":
            variant = source.name.removeprefix("Stephenson_pseudobulk_").removesuffix(".rds")
            relative = Path("pseudobulks") / f"Stephenson_pseudobulk_{variant}.rds"
            producer = "stage5_prepare_pseudobulk_" + variant
        elif kind == "metadata":
            relative = Path("metadata") / f"{dataset}_sample_metadata.feather"
            producer = "stage5_metadata_adopted_copy"
        elif kind == "timing":
            relative = Path("execution_times_prepare_pseudobulk_Stephenson.feather")
            producer = "stage5_execution_log"
        else:
            raise ValueError("prepared adoption: unknown payload kind")
        key = (ctx, kind, dataset, view, variant)
        require(key in expected and key not in observed, "foreign/duplicate payload scope")
        observed.add(key)
        expected_source = old_root / (Path("logs") / relative if kind == "timing" else Path("artifacts") / ctx / relative)
        require(source == expected_source, "source path differs from exact producer namespace")
        require(re.fullmatch("[0-9a-f]{32}", item["md5"]) and type(item["size"]) is int
                and item["size"] > 0 and source.stat().st_size == item["size"] and digest(source, "md5") == item["md5"], "source content differs")
        sidecar = checked(dict(path=str(source) + ".md5", sha256=item["sidecar_sha256"]))
        da.proof(sidecar, source, item["md5"], item["size"])
        old_input = old_inputs[(dataset, view)]
        if ctx == context:
            require(current_inputs.get((dataset, view)) == old_input, "selected H5AD identity differs")
        record = old_root / "manifests/artifacts" / (hashlib.sha256(str(source).encode()).hexdigest()[:32] + ".record")
        if item["record_sha256"] is not None:
            checked(dict(path=str(record), sha256=item["record_sha256"]))
            require(fields(record) == dict(PATH=str(source), MD5=item["md5"], SIZE=str(item["size"]),
                                          RUN_ID=old_root.name, PRODUCER=producer, STATE="PUBLISHED"), "source producer record differs")
        else:
            require(kind == "metadata" and not record.exists(), "missing/undeclared source producer record")
        owner = owner_path(source, scratch)
        if item["owner_sha256"] is not None:
            checked(dict(path=str(owner), sha256=item["owner_sha256"]))
            values = fields(owner)
            require(values.get("PATH") == str(source) and values.get("RUN_ID") == old_root.name
                    and values.get("STAGE") == "stage5" and values.get("STATE") in ("OK", "FAIL"), "source has active/foreign owner")
        else:
            require(kind != "pseudobulk" and not owner.exists(), "missing/undeclared source owner")
        if kind == "metadata":
            status_path = checked(item["status"])
            require(status_path == old_root / "status/metadata_export" / f"{dataset}__{view}.status", "metadata status namespace differs")
            status = fields(status_path)
            require(all(status.get(k) == v for k, v in dict(STATE="OK", RUN_ID=old_root.name, DATASET=dataset,
                    VIEW=view, INPUT_FILE=old_input[0], OUTPUT_FILE=str(source),
                    ANALYSIS_ROOT=str(old_root / "artifacts" / ctx),
                    ANALYSIS_PASS=view.removeprefix("batch_effect_"),
                    ANALYSIS_VARIANT="final" if ctx.endswith("uncorrected_final") else "corrected_final").items())
                    and status.get("TASK_ID") == str(metadata_order[ctx].index((dataset, view)) + 1), "metadata status differs")
            if metadata_only:
                require(status.get("STATUS") == "EXPORTED", "metadata-only adoption requires an original exported task")
        else:
            require(item["status"] is None, "unexpected metadata status")
        if ctx == context:
            target = Path(current["run_root"]) / "logs" / relative if kind == "timing" else da.context_root(current, ctx) / relative
            selected.append(dict(kind=kind, dataset=dataset, view=view, source=str(source), target=str(target),
                                 producer=producer, md5=item["md5"], size=item["size"]))
    require(observed == expected and len(receipt["artifacts"]) == len(expected), "finite payload coverage differs")
    return receipt, current, selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "copy"))
    for name in ("receipt", "sha256", "config", "context"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    receipt, binding, selected = prepare(args.receipt, args.sha256, args.config, args.context)
    report_path = Path(binding["run_root"]) / "manifests/prepared_adoption.json"
    if args.command == "plan":
        require(not report_path.exists(), "adoption report already exists; preserved")
        for item in selected:
            target = da.absolute(item["target"], exists=False)
            require(all(not p.exists() and not p.is_symlink() for p in (target, Path(str(target) + ".md5"))), "target already exists; preserved")
            print("\t".join(str(item[key]) for key in ("kind", "dataset", "view", "source", "target", "producer", "md5", "size")))
        return
    # Validate all old proofs/content again before installation. Never overwrite.
    for item in selected:
        target = da.absolute(item["target"], exists=False)
        require(not target.exists() and not Path(str(target) + ".md5").exists(), "target already exists; preserved")
    for item in selected:
        target = Path(item["target"])
        target.parent.mkdir(parents=True, exist_ok=True)
        da.absolute(target.parent)
        fd, name = tempfile.mkstemp(prefix=".prepared-copy-", dir=target.parent)
        os.close(fd)
        temporary = Path(name)
        try:
            shutil.copyfile(item["source"], temporary)
            require(digest(temporary, "md5") == item["md5"] and temporary.stat().st_size == item["size"], "copied content differs")
            os.link(temporary, target)
            da.create_only(Path(str(target) + ".md5"), f"MD5={item['md5']}\nSIZE={item['size']}\nPATH={target}\n".encode())
        finally:
            temporary.unlink(missing_ok=True)
    report = dict(schema=receipt["schema"], operation="COPY_WITHOUT_REFIT", original_gate_state="FAILED",
                  receipt_path=args.receipt, receipt_sha256=args.sha256, binding_sha256=digest(Path(binding["run_root"]) / "manifests/artifact_bindings.json"),
                  fit_timing=("metadata byte-copy only; no cache, timing payload or scientific refit" if receipt["schema"] == METADATA_SCHEMA else
                              "original cache/timing bytes and timing IDs preserved; copy cost is not a fit"), artifacts=selected)
    da.create_only(report_path, (json.dumps(report, sort_keys=True, indent=2) + "\n").encode())


if __name__ == "__main__":
    main()
