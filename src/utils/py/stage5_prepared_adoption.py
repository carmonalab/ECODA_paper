"""Copy reviewed preparation artifacts into exact fresh derivative namespaces.

The original23 and selected metadata-only contracts remain distinct. The fixed
batch32 contract retains v2 STOP preparation and v1 FAILED metadata source roles.
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
BATCH_SCHEMA = "ecoda_stage5_batch_prepared_adoption_v1"
BATCH_TARGETS = {"batch_effect/uncorrected_final": "stage5_remaining95_uncorrected_20261109_v1",
                 "batch_effect/corrected_final": "stage5_remaining95_corrected_20261109_v1"}
BATCH_ROLES = {
    "uncorrected_preparation": ("batch_effect/uncorrected_final", "stage5_main_uncorrected_20261109_v2", "PRELAUNCH_STOP", "9ca2453eea7320a721ef6dfd458f9db108abecac", "ecoda_stage5_main_confounds_20261109_v2"),
    "uncorrected_metadata": ("batch_effect/uncorrected_final", "stage5_main_uncorrected_20261109_v1", "FAILED", "c534cdfd28fc9a77c009f4e25fa4c4ddb7ee4400", "ecoda_stage5_main_confounds_20261109_v1"),
    "corrected_metadata": ("batch_effect/corrected_final", "stage5_main_corrected_20261109_v1", "FAILED", "c534cdfd28fc9a77c009f4e25fa4c4ddb7ee4400", "ecoda_stage5_main_confounds_20261109_v1"),
}


def scope_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def batch_review(receipt, batch=None):
    require("original_gate_state" not in receipt and receipt["target_run_ids"] == BATCH_TARGETS
            and set(receipt["sources"]) == set(BATCH_ROLES), "fixed batch source-role/target scope differs")
    batch = set(batch) if batch is not None else {item["dataset"] for item in receipt["artifacts"] if item["kind"] == "pseudobulk"}
    expected = {(ctx, "metadata", dataset, "batch_effect_" + ("uncorrected" if ctx.endswith("uncorrected_final") else "corrected"),
                 "uncorrected_metadata" if ctx.endswith("uncorrected_final") else "corrected_metadata")
                for ctx in BATCH_TARGETS for dataset in batch}
    expected.update(("batch_effect/uncorrected_final", kind, dataset, "batch_effect_uncorrected", "uncorrected_preparation")
                    for dataset in batch for kind in ("pseudobulk", "timing"))
    actual = [(item["context"], item["kind"], item["dataset"], item["view"], item["source_role"]) for item in receipt["artifacts"]]
    require(len(batch) == 8 and len(actual) == len(set(actual)) == 32 and set(actual) == expected, "batch fixed32 payload/source-role coverage differs")
    review = json.loads(checked(receipt["reuse_review"]).read_text())
    require(len(receipt["evidence"]) == 3, "batch requires cache/metadata/scoped-writer evidence")
    cache = json.loads(checked(receipt["evidence"][0]).read_text())
    require(cache["schema"] == "ecoda_unc8_cache_semantics_audit_v1"
            and cache["cache_semantics_validated"] is True and cache["timing_semantics_validated"] is True
            and type(cache["actual_cli"]["status"]) is int and cache["actual_cli"]["status"] == 0
            and cache["reuse_approved"] is False and cache["original_gate_state"] == "PRELAUNCH_STOP"
            and cache["original_wrapper_exit_code"] is None
            and Path(cache["original_run_root"]).name == BATCH_ROLES["uncorrected_preparation"][1], "batch cache semantics are not established")
    metadata = json.loads(checked(cache["metadata"]).read_text())
    prepared = [item for item in receipt["artifacts"] if item["kind"] in ("pseudobulk", "timing")]
    observed = {(row["dataset"], kind): row["cache" if kind == "pseudobulk" else "timing"]
                for row in metadata["rows"] for kind in ("pseudobulk", "timing")}
    require(len(metadata["rows"]) == len(cache["rows"]) == 8
            and len(observed) == 16 and {row["dataset"] for row in cache["rows"]} == {row["dataset"] for row in metadata["rows"]}
            and all(all(row["identity"][key] == observed[(row["dataset"], "pseudobulk")][key]
                        and row["timing_file"]["identity"][key] == observed[(row["dataset"], "timing")][key]
                        for key in ("md5", "size")) for row in cache["rows"])
            and all((item["dataset"], item["kind"]) in observed
                    and all(observed[(item["dataset"], item["kind"])][key] == item["source" if key == "path" else key]
                            for key in ("path", "md5", "size")) for item in prepared), "cache/timing science evidence differs from selected bytes")
    writers = json.loads(checked(receipt["evidence"][2]).read_text())
    require(writers["schema"] == "ecoda_stage5_batch_prepared_writer_review_v1"
            and writers["status"] == "APPROVED_SCOPED_NO_STALE_WRITERS"
            and writers["source_run_root"] == cache["original_run_root"]
            and writers["original_gate_state"] == "PRELAUNCH_STOP" and writers["wrapper_exit_code"] is None
            and writers["scoped_no_stale_writers"] is True and writers["artifacts_sha256"] == scope_digest(prepared)
            and writers["evidence"], "scoped preparation-writer review is missing or unapproved")
    observations = []
    for item in writers["evidence"]:
        path = checked(item)
        if path.suffix == ".json":
            value = json.loads(path.read_text())
            if isinstance(value, dict) and value.get("schema") == "ecoda_unc8_scoped_writer_observation_v1":
                observations.append(value)
    canonical = sorted((dict(dataset=row["dataset"], kind=kind,
                             path=observed[(row["dataset"], kind)]["path"],
                             **{key: identity[key] for key in ("md5", "sha256", "size")})
                        for row in cache["rows"] for kind, identity in
                        (("pseudobulk", row["identity"]), ("timing", row["timing_file"]["identity"]))),
                       key=lambda item: (item["dataset"], item["kind"], item["path"]))
    require(len(observations) == 1, "actual scoped writer observation is missing/ambiguous")
    observation = observations[0]
    require(observation["status"] == "READONLY_CANDIDATE" and observation["artifact_identity_stable"] is True
            and observation["blockers"] == [] and observation["scoped_no_stale_writers_approved"] is False
            and observation["scientific_reuse_approved"] is False
            and observation["source_run_root"] == cache["original_run_root"]
            and observation["original_gate_state"] == "PRELAUNCH_STOP" and observation["wrapper_exit_code"] is None
            and observation["canonical_artifacts"] == canonical and observation["canonical_artifacts_sha256"] == scope_digest(canonical)
            and cache["metadata"] in observation["control_proofs"] and receipt["evidence"][0] in observation["control_proofs"],
            "scoped writer observation does not bind these exact16 preparation artifacts")
    for item in observation["control_proofs"]:
        checked(item)

    require(review["schema"] == "ecoda_stage5_batch_prepared_review_v1" and review["status"] == "APPROVED_COPY_ONLY"
            and review["config_sha256"] == receipt["config_sha256"]
            and review["sources_sha256"] == scope_digest(receipt["sources"])
            and review["artifacts_sha256"] == scope_digest(receipt["artifacts"])
            and review["evidence"] == receipt["evidence"]
            and review["source_gate_states"] == {role: values[2] for role, values in BATCH_ROLES.items()}
            and (review["cache_artifacts"], review["timing_artifacts"], review["metadata_artifacts"]) == (8, 8, 16)
            and all(review.get(key) is True for key in ("cache_semantics_validated", "timing_semantics_validated",
                    "metadata_semantics_validated", "scoped_no_stale_writers")), "batch preparation/scoped-writer review is missing or unapproved")

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
    batch_prepared = receipt["schema"] == BATCH_SCHEMA
    require(receipt["schema"] in (SCHEMA, METADATA_SCHEMA, BATCH_SCHEMA) and receipt["state"] == "APPROVED_COPY_ONLY"
            and (batch_prepared or receipt["original_gate_state"] == "FAILED"), "receipt purpose/state differs")
    if batch_prepared:
        batch_review(receipt)
    metadata_only = receipt["schema"] == METADATA_SCHEMA
    config_path = da.absolute(config_path)
    require(receipt["config_sha256"] == digest(config_path), "configuration differs")
    current, _ = da.active_binding(config_path)
    contexts = set(receipt["target_run_ids"])
    require(context in contexts and current["run_id"] == receipt["target_run_ids"][context], "target run differs")
    require((contexts == set(BATCH_TARGETS) and set(receipt["sources"]) == set(BATCH_ROLES)) if batch_prepared else
            (set(receipt["sources"]) == contexts and
             (bool(contexts) and contexts <= set(CONTEXTS[1:]) if metadata_only else contexts == set(CONTEXTS))),
            "context coverage differs")
    require(len(set(receipt["target_run_ids"].values())) == len(contexts), "target runs must be distinct")
    require(len(receipt["evidence"]) == 3 and len({item["path"] for item in receipt["evidence"]}) == 3,
            "missing/duplicate reviewed scientific/accounting evidence")
    for evidence in receipt["evidence"]:
        checked(evidence)
    config = json.loads(config_path.read_text())
    batch = {key for key, value in config.items() if not key.startswith("_") and value.get("use_for_batch_effect") is True}
    require(len(batch) == 8, "this finite recovery requires exactly eight configured batch cohorts")
    if batch_prepared:
        expected = set()
        for ctx in BATCH_TARGETS:
            view = "batch_effect_" + ("uncorrected" if ctx.endswith("uncorrected_final") else "corrected")
            expected.update((ctx, "metadata", dataset, view, "") for dataset in batch)
        expected.update((CONTEXTS[1], kind, dataset, "batch_effect_uncorrected", "hvg2000" if kind == "pseudobulk" else "")
                        for dataset in batch for kind in ("pseudobulk", "timing"))
        view = "batch_effect_" + ("uncorrected" if context.endswith("uncorrected_final") else "corrected")
        require({tuple(line.split("\t")[:2]) for line in Path(current["selection"]).read_text().splitlines()} == {(ds, view) for ds in batch},
                "batch prepared metadata must cover the exact eight target cohorts")
    elif metadata_only:
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
    for role, source in receipt["sources"].items():
        ctx = BATCH_ROLES[role][0] if batch_prepared else role
        old_root = da.absolute(source["run_root"])
        require(old_root.parent == scratch / "_ecoda_runs" and old_root.name != receipt["target_run_ids"][ctx], "source run root differs")
        binding_path = checked(source["binding"])
        old_config = checked(source["config"])
        require(binding_path == old_root / "manifests/artifact_bindings.json" and digest(old_config) == digest(config_path), "source binding/config differs")
        binding, _ = da.load(binding_path, old_config)
        require(binding["run_root"] == str(old_root) and binding["stage"] == "stage5", "source binding stage differs")
        for key, name in (("source_manifest", "source.manifest"), ("runtime_identity", "runtime.identity")):
            require(checked(source[key]) == old_root / "manifests" / name, "source control path differs")
        if batch_prepared:
            _, source_run, gate_state, source_commit, gate_id = BATCH_ROLES[role]
            gate = json.loads(checked(source["gate_manifest"]).read_text())
            manifest = fields(old_root / "manifests/source.manifest")
            require(source["context"] == ctx and source["gate_state"] == gate_state and old_root.name == source_run
                    and manifest["SOURCE_COMMIT"] == source_commit and gate["gate_id"] == gate_id and gate["state"] == gate_state,
                    "batch source role/run/commit/gate state differs")
            if role == "uncorrected_preparation":
                require(gate.get("terminal_status") is None and gate.get("remote_exit_code") is None,
                        "stopped preparation must retain its unknown wrapper outcome")
        runtime = fields(old_root / "manifests/runtime.identity")
        require(all(re.fullmatch("[0-9a-f]{64}", runtime.get(key, "")) and runtime[key] == current_runtime.get(key)
                    for key in ("RUNTIME_IMAGE_SHA256", "RUNTIME_MANIFEST_SHA256")), "runtime image/manifest differs")
        old[role] = (old_root, input_rows(binding))
        metadata_order[role] = list(dict.fromkeys(tuple(line.split("\t")[:2])
            for line in Path(binding["selection"]).read_text().splitlines()))
    selected = []
    observed = set()
    for item in receipt["artifacts"]:
        ctx, kind, dataset, view = (item[key] for key in ("context", "kind", "dataset", "view"))
        source = da.absolute(item["source"])
        role = item["source_role"] if batch_prepared else ctx
        if batch_prepared:
            expected_role = "uncorrected_preparation" if kind in ("pseudobulk", "timing") else ("uncorrected_metadata" if ctx == CONTEXTS[1] else "corrected_metadata")
            require(role == expected_role and role in BATCH_ROLES and ctx == BATCH_ROLES[role][0], "payload source role differs from kind/context")
        old_root, old_inputs = old[role]
        variant = ""
        if kind == "pseudobulk" and batch_prepared:
            variant = "hvg2000"
            relative = Path("pseudobulks") / f"{dataset}_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds"
            producer = "stage5_prepare_pseudobulk_hvg2000"
        elif kind == "pseudobulk":
            variant = source.name.removeprefix("Stephenson_pseudobulk_").removesuffix(".rds")
            relative = Path("pseudobulks") / f"Stephenson_pseudobulk_{variant}.rds"
            producer = "stage5_prepare_pseudobulk_" + variant
        elif kind == "metadata":
            relative = Path("metadata") / f"{dataset}_sample_metadata.feather"
            producer = "stage5_metadata_adopted_copy"
        elif kind == "timing":
            relative = Path(f"execution_times_batch_effect_uncorrected_final_prepare_pseudobulk_{dataset}.feather") if batch_prepared else Path("execution_times_prepare_pseudobulk_Stephenson.feather")
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
            require(kind != "pseudobulk" and not owner.exists() and not owner.is_symlink(), "missing/undeclared source owner")
        if batch_prepared and kind == "pseudobulk":
            require(item["owner_sha256"] is not None and fields(owner).get("STATE") == "FAIL", "batch cache source owner must remain exactly FAIL")
        if batch_prepared and kind == "timing":
            require(item["owner_sha256"] is None and not owner.exists() and not owner.is_symlink(), "original preparation timing owner must remain absent")
        if kind == "metadata":
            status_path = checked(item["status"])
            require(status_path == old_root / "status/metadata_export" / f"{dataset}__{view}.status", "metadata status namespace differs")
            status = fields(status_path)
            require(all(status.get(k) == v for k, v in dict(STATE="OK", RUN_ID=old_root.name, DATASET=dataset,
                    VIEW=view, INPUT_FILE=old_input[0], OUTPUT_FILE=str(source),
                    ANALYSIS_ROOT=str(old_root / "artifacts" / ctx),
                    ANALYSIS_PASS=view.removeprefix("batch_effect_"),
                    ANALYSIS_VARIANT="final" if ctx.endswith("uncorrected_final") else "corrected_final").items())
                    and status.get("TASK_ID") == str(metadata_order[role].index((dataset, view)) + 1), "metadata status differs")
            if metadata_only or batch_prepared:
                require(status.get("STATUS") == "EXPORTED", "metadata-only adoption requires an original exported task")
        else:
            require(item["status"] is None, "unexpected metadata status")
        if ctx == context:
            target = Path(current["run_root"]) / "logs" / relative if kind == "timing" else da.context_root(current, ctx) / relative
            selected.append(dict(kind=kind, dataset=dataset, view=view, source=str(source), target=str(target),
                                 producer=producer, md5=item["md5"], size=item["size"],
                                 **({"source_role": role} if batch_prepared else {})))
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
    report = dict(schema=receipt["schema"], operation="COPY_WITHOUT_REFIT",
                  receipt_path=args.receipt, receipt_sha256=args.sha256, binding_sha256=digest(Path(binding["run_root"]) / "manifests/artifact_bindings.json"),
                  fit_timing=("metadata byte-copy only; no cache, timing payload or scientific refit" if receipt["schema"] == METADATA_SCHEMA else
                              "original cache/timing bytes and timing IDs preserved; copy cost is not a fit"), artifacts=selected)
    if receipt["schema"] == BATCH_SCHEMA:
        report["source_gate_states"] = {role: values[2] for role, values in BATCH_ROLES.items()}
        report["reuse_review"] = receipt["reuse_review"]
        report["fit_timing"] = "byte-copy only; canonical cache class retained with source_role lineage; original timing IDs/costs unchanged, no new fit"
    else:
        report["original_gate_state"] = "FAILED"
    da.create_only(report_path, (json.dumps(report, sort_keys=True, indent=2) + "\n").encode())


if __name__ == "__main__":
    main()
