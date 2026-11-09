#!/usr/bin/env python3
"""Finite Stage 3 artifact reuse evidence; never releases or changes a gate."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import stat
import sys
import tempfile


def require(value, message):
    if not value:
        raise ValueError(message)


def regular(path):
    p = Path(path)
    require(p.is_absolute() and p.is_file() and not p.is_symlink() and p.resolve() == p,
            f"Expected canonical regular file: {p}")
    return p


def identity(path):
    p = regular(path)
    return {"path": str(p), "size": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def kv(path):
    result = {}
    for line in regular(path).read_text().splitlines():
        key, sep, value = line.partition("=")
        require(sep and key not in result, f"Malformed/duplicate field: {path}/{key}")
        result[key] = value
    return result


def rows(path, width):
    result = [tuple(line.split("\t")) for line in regular(path).read_text().splitlines()]
    require(result and all(len(r) == width and all(r) for r in result) and len(set(result)) == len(result),
            f"Invalid finite manifest: {path}")
    return result


def option(tokens, name):
    require(tokens.count(name) == 1, f"Expected exact wrapper option {name}")
    return tokens[tokens.index(name) + 1]


def history(args):
    run = Path(args.run_root)
    gate = json.loads(regular(args.failed_gate_manifest).read_text())
    inspect = json.loads(regular(args.terminal_inspect).read_text())
    require(gate["state"] == inspect["state"] == "FAILED" and gate["gate_id"] == inspect["gate_id"], "Original failed gate identity differs")
    require(gate["audit_completed"] and gate["audit_state"] == "COMPLETED" and gate["audit"] == inspect["audit"], "Original terminal audit differs/incomplete")
    terminal = gate["terminal_status"]
    require(terminal["state"] == "COMPLETED" and terminal["exit_code"] == 0 and
            terminal["gate_id"] == gate["gate_id"] and terminal["exact_command"] == gate["exact_command"] and
            terminal["command_digest"] == gate["command_digest"] == hashlib.sha256(gate["exact_command"].encode()).hexdigest() and
            terminal["event_generation"] == gate["event_generation"], "Original wrapper did not complete with exact identity")
    audit = gate["audit"]
    require(all(x["ok"] for k in ("commands", "artifact_contracts", "immutable_fingerprints") for x in audit[k]), "Original audit has a non-accounting failure")
    require(gate["discrepancies"] and all(x["kind"] == "terminal_audit_failed" for x in gate["discrepancies"]), "Unexpected original discrepancy")
    tokens = shlex.split(gate["exact_command"])
    require(option(tokens, "--run-id") == run.name and option(tokens, "--scratch-root") == str(run.parent.parent) and
            option(tokens, "--script") == "src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" and
            regular(option(tokens, "--source-manifest")).read_bytes() == regular(args.source_manifest).read_bytes(), "Wrapper producer/source scope differs")
    require(option(tokens, "--mem") == "128G" and option(tokens, "--max-mem") == "500G", "Unexpected bounded retry memory policy")
    source_config = json.loads(regular(Path(kv(args.source_manifest)["SOURCE_ROOT"]) / "datasets.json").read_text())
    expected_original = {(d, v) for d, entry in source_config.items() if not d.startswith("_") and entry.get("use_for_batch_effect") for v in ("batch_effect_uncorrected", "batch_effect_corrected")} | {("Stephenson", "benchmark_analysis")}
    original = rows(args.selection, 2)
    selected = rows(args.reuse_selection, 2)
    require(set(selected) <= set(original), "Reuse scope is not an exact original selection subset")
    require(regular(args.selection) == run / "manifests/selection.tsv" and len(original) == 17 and set(original) == expected_original,
            "Reuse mode requires this original seventeen-row Stage 3 cohort run")
    require(regular(option(tokens, "--selection-file")).read_bytes() == regular(args.selection).read_bytes(), "Wrapper selection differs from original run")
    scheduler = rows(run / "manifests/scheduler_ids.tsv", 2)
    require([r[0] for r in scheduler] == ["PREFLIGHT", "ARRAY", "WATCHDOG", "STATUS", "STATUS"] and
            all(re.fullmatch(r"[0-9]+", r[1]) for r in scheduler) and len({r[1] for r in scheduler}) == 5, "Unexpected retry scheduler lineage")
    ids = [r[1] for r in scheduler]
    require(set(ids) == set(gate["scheduler_array_ids"] + gate["scheduler_watchdog_ids"]), "Gate/run scheduler identities differ")
    accounting = audit["accounting"]
    require(accounting["query_count"] == 1 and accounting["actual_exit"] == 0 and accounting["attempted"] and
            not accounting["malformed_rows"] and not accounting["missing_scheduler_ids"] and set(accounting["scheduler_ids"]) == set(ids), "Original single accounting query incomplete")
    query = shlex.split(accounting["command"])
    require(query[0] == "sacct" and set(option(query, "-j").split(",")) == set(ids) and
            "--format=JobIDRaw,State,ExitCode" in query, "Original accounting query scope differs")
    raw = [line.split("|") for line in accounting["stdout"].splitlines()]
    require(len(raw) == 24 and all(len(r) == 3 and re.fullmatch(r"[0-9]+", r[0]) for r in raw) and
            len({r[0] for r in raw}) == 24, "Missing/extra/malformed original accounting rows")
    require([(r["job_id"], r["state"], r["exit_code"]) for r in accounting["rows"]] == [tuple(r) for r in raw], "Raw/parsed accounting rows differ")
    states = {r[0]: r[1:] for r in raw}
    require(all(r[1:] in (["COMPLETED", "0:0"], ["OUT_OF_MEMORY", "0:125"]) for r in raw), "Non-OOM/unresolved accounting failure")
    require(all(states.get(ids[i]) == ["COMPLETED", "0:0"] for i in (0, 1, 2, 4)), "Preflight/initial/final/watchdog root not successful")
    retry_paths = [run / f"manifests/selection.retry_{i}.tsv" for i in (1, 2)]
    retries = [rows(p, 2) for p in retry_paths]
    expected = {("Alzheimer", "batch_effect_uncorrected"), ("Alzheimer", "batch_effect_corrected")}
    require(all(len(r) == 2 and set(r) == expected for r in retries) and expected <= set(original), "Unexpected retry cohort scope")
    preflight = rows(run / "manifests/h5ad_obs_audit.tsv", 3)
    require(len(preflight) == 2 and {(r[0], r[1]) for r in preflight} == {("Covid19_PBMC", "batch_effect_uncorrected"), ("Covid19_PBMC", "batch_effect_corrected")}, "Unexpected preflight task scope")
    log_root = Path(option(tokens, "--logs-root"))
    require(log_root == run.parent.parent / "_ecoda_run_logs" / run.name and log_root.resolve() == log_root, "Unexpected run log root")
    witnesses = []
    for attempt, previous, current, array in ((0, original, retries[0], ids[1]), (1, retries[0], retries[1], ids[3])):
        for pair in current:
            task = previous.index(pair) + 1
            prefix = "3_scrnaseq_preprocessing_" if attempt == 0 else "3_scrnaseq_preprocessing_retry1_"
            p = log_root / f"{prefix}{array}_{task}.err"
            require(regular(p).stat().st_size <= 1024 * 1024, "OOM witness exceeds bounded text scope")
            text = p.read_text()
            found = set(re.findall(r"StepId=([0-9]+)\.batch", text))
            require(len(found) == 1 and re.search(r"oom|out.of.memory", text, re.I), "Missing/ambiguous exact OOM StepId witness")
            raw_id = found.pop()
            require(states.get(raw_id) == ["OUT_OF_MEMORY", "0:125"], "OOM witness does not match original accounting")
            witnesses.append(dict(identity(p), raw_job_id=raw_id, array_id=array, task_id=task, dataset=pair[0], view=pair[1]))
    require(len({w["raw_job_id"] for w in witnesses}) == 4 and
            {w["raw_job_id"] for w in witnesses} == {r[0] for r in raw if r[1] == "OUT_OF_MEMORY"}, "Unwitnessed/extra historical OOM")
    watchdog = regular(run / "status/watchdog").read_text().splitlines()
    fields = [line.split("=", 1) for line in watchdog]
    for key, value in (("STATE", "OK"), ("RUN_ID", run.name), ("RETRY_INDEX", "2"), ("ARRAY_JOB_ID", ids[4])):
        require([v for k, v in fields if k == key] == [value], "Watchdog final retry proof differs: " + key)
    require({v for k, v in fields if k == "SCHEDULER_ID"} == set(ids[1:]), "Watchdog scheduler lineage differs")
    terminal_status = kv(run / "status/terminal")
    require(terminal_status["STATE"] == "OK" and terminal_status["RUN_ID"] == run.name, "Original stage terminal is not OK")
    proof_paths = [args.failed_gate_manifest, args.terminal_inspect, args.selection, args.reuse_selection,
                   run / "manifests/scheduler_ids.tsv", run / "status/watchdog", run / "status/terminal", *retry_paths]
    return {"gate_id": gate["gate_id"], "gate_state": "FAILED", "wrapper_terminal": terminal,
            "original_accounting_stdout": accounting["stdout"], "accounting_rows": raw,
            "oom_witnesses": witnesses, "source_defined_memory_progression": ["128G", "256G", "500G"], "proofs": [identity(p) for p in proof_paths], "selected": selected}


def md5(path):
    digest = hashlib.md5()
    with regular(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checksum(path):
    p = regular(path)
    proof = kv(str(p) + ".md5")
    require(set(proof) == {"MD5", "SIZE", "PATH"} and proof["PATH"] == str(p) and
            proof["SIZE"] == str(p.stat().st_size) and proof["MD5"] == md5(p), "Read-only checksum mismatch: " + str(p))
    return {"path": str(p), "md5": proof["MD5"], "size": p.stat().st_size, "proof": identity(str(p) + ".md5")}


def ownership(path, scratch, run_id, record):
    p = regular(path); digest = hashlib.sha256(str(p).encode()).hexdigest()[:32]
    safe = re.sub(r"[/ :,\t|]", "_", p.name)
    key = f"{safe}_{digest}" if len(f"{safe}_{digest}".encode()) <= 100 else f"artifact_{digest}"
    owner_path = scratch / "_ecoda_owners/artifact" / key / "owner"
    owner = kv(owner_path)
    require(owner["PATH"] == str(p) and owner["RUN_ID"] == run_id and owner["STATE"] == "OK", "Artifact owner differs: " + str(p))
    result = {"owner": identity(owner_path), "owner_fields": owner}
    if record:
        r = scratch / "_ecoda_runs" / run_id / "manifests/artifacts" / f"{digest}.record"
        value = kv(r)
        require(set(value) == {"PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE"} and value["PATH"] == str(p) and value["RUN_ID"] == run_id and value["PRODUCER"] == "stage3" and value["STATE"] == "PUBLISHED" and value["SIZE"] == str(p.stat().st_size) and value["MD5"] == kv(str(p) + ".md5")["MD5"], "Stage 3 producer differs")
        result.update(record=identity(r), record_fields=value)
    return result


def write_report(args, evidence):
    import h5py  # Metadata only; no X/count/layer arrays are opened.
    run = Path(args.run_root); scratch = run.parent.parent
    binding_path = run / "manifests/artifact_bindings.json"
    binding = json.loads(regular(binding_path).read_text())
    source = kv(args.source_manifest)["SOURCE_ROOT"]
    config_path = Path(source) / "datasets.json"
    config = json.loads(regular(config_path).read_text())
    refs = {r[:2]: r for r in rows(run / "manifests/input_artifacts.tsv", 6)}
    sys.path.insert(0, str(Path(source) / "src/utils/py"))
    from datasets_io import read_datasets_json
    from cohort_selection import cohort_configuration, validate_cohort_h5ad_path
    registry = read_datasets_json(config_path)
    artifacts = []
    for dataset, view in evidence["selected"]:
        name = config[dataset]["views"][view]["output_file_name"]
        output = Path(binding["artifact_root"]) / "datasets" / dataset / "output" / name
        nas = Path(binding["nas_root"]) / "datasets" / dataset / "output" / name
        item = {"dataset": dataset, "view": view, "scratch": checksum(output), "nas": checksum(nas)}
        require(item["scratch"]["md5"] == item["nas"]["md5"] and item["scratch"]["size"] == item["nas"]["size"], "NAS/scratch checksum identity differs")
        item["scratch"].update(ownership(output, scratch, run.name, True))
        item["nas"].update(ownership(nas, scratch, run.name, False))
        require(stat.S_IMODE(nas.stat().st_mode) == stat.S_IMODE(Path(str(nas) + ".md5").stat().st_mode) == 0o666, "Selected NAS file/proof mode differs")
        directory = nas.parent; directories = []
        while directory != Path(binding["nas_base"]):
            require(directory.is_relative_to(Path(binding["nas_base"])) and not directory.is_symlink() and stat.S_IMODE(directory.stat().st_mode) == 0o777, "Selected NAS directory mode differs")
            directories.append(str(directory)); directory = directory.parent
        item["nas"]["directories_mode_777"] = directories
        ref = refs[(dataset, view)]
        source_path = regular(ref[2]); before = source_path.stat(); proof_identity = identity(ref[5])
        proof = kv(ref[5])
        require(proof == {"MD5": ref[3], "SIZE": ref[4], "PATH": ref[2]} and before.st_size == int(ref[4]), "Selected production input reference/proof differs")
        entry = registry[dataset]
        configuration = cohort_configuration(dataset=dataset, view=view, sample_column=entry["sample_col"],
            label_column=entry["label_col"], batch_keys=entry["batch_col"], subset_vars=entry["views"][view]["subset_vars"],
            batch_metadata_policy=entry["batch_metadata_policy"], evaluation_labels=entry["evaluation_labels"])
        validate_cohort_h5ad_path(output, configuration, [{"path": ref[2], "size": int(ref[4]), "md5": ref[3]}])
        after = source_path.stat()
        require((before.st_size, before.st_mtime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ino) and
                identity(ref[5]) == proof_identity, "Historical raw source/proof changed during validation")
        item["input"] = {"path": ref[2], "md5": ref[3], "size": int(ref[4]), "proof": proof_identity,
                         "verification": "original production source MD5/reference reused; no fresh raw content hash",
                         "stat_before_after": {"size": before.st_size, "mtime_ns": before.st_mtime_ns, "inode": before.st_ino}}
        with h5py.File(output, "r") as handle:
            record = handle["uns/cohort_selection"]
            def scalar(key):
                x = record[key][()]
                return x.decode() if isinstance(x, bytes) else x.item() if hasattr(x, "item") else x
            item["cohort"] = {k: scalar(k) for k in ("schema_version", "configuration_sha256", "membership_sha256", "output_observations_sha256")}
            item["cohort"]["configuration"] = json.loads(scalar("configuration_json"))
            item["cohort"]["membership"] = json.loads(scalar("membership_json"))
            item["cohort"]["source_identities"] = json.loads(scalar("source_identity_json"))
            item["n_samples"] = len(item["cohort"]["membership"]["retained_after_cellcount_ids"])
        artifacts.append(item)
    report = {"schema": "ecoda_stage3_subset_reuse_v1", "status": "NOOP_VALIDATED", "purpose": "individually_validated_stage3_artifact_reuse", "run_id": run.name,
              "failed_gate_released": False, "scheduler_queries": 0, "scheduler_submissions": 0, "history": evidence,
              "identities": [identity(p) for p in (args.source_manifest, args.runtime_identity, binding_path, config_path, __file__, args.audit_entrypoint)],
              "artifacts": artifacts}
    destination = Path(args.stage3_reuse_report)
    require(destination.is_absolute() and destination.parent.resolve() == destination.parent and destination.parent.is_dir() and not destination.is_relative_to(run), "Receipt must be outside original run in a canonical existing directory")
    payload = (json.dumps(report, sort_keys=True) + "\n").encode()
    publish(destination, payload)
    proof = f"MD5={hashlib.md5(payload).hexdigest()}\nSIZE={len(payload)}\nPATH={destination}\n".encode()
    publish(Path(str(destination) + ".md5"), proof)


def publish(destination, payload):
    descriptor, temporary = tempfile.mkstemp(prefix="." + destination.name + ".", dir=destination.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload); stream.flush(); os.fsync(stream.fileno())
        os.chmod(temporary, 0o444)
        os.link(temporary, destination)  # Create-only; preserve any earlier evidence.
    finally:
        os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run-root", "selection", "reuse-selection", "source-manifest", "runtime-identity", "terminal-inspect", "failed-gate-manifest"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--check-history", action="store_true")
    parser.add_argument("--stage3-reuse-report")
    parser.add_argument("--audit-entrypoint")
    args = parser.parse_args()
    try:
        evidence = history(args)
        if not args.check_history:
            require(args.stage3_reuse_report and args.audit_entrypoint, "Explicit report and audit entrypoint required")
            write_report(args, evidence)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
