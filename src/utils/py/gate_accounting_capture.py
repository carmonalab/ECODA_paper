"""Bridge one completed gate inspect to the run audit without another query.

The raw JobIDRaw output remains unchanged. The audit input is explicitly a
derived projection of actual requested-root rows, after *every* raw row passes.
Single-child mode requires exact gate/child ID equality. Campaign mode validates
the sealed three-context receipts and writes one projection per child audit.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import sys
import tempfile


ACCOUNTING_COMMAND = "sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode"
ID = re.compile(r"[1-9][0-9]*")


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _file(path, *, allow_empty=False):
    path = Path(path)
    if not path.is_absolute() or path.is_symlink() or not path.is_file() or path.resolve() != path:
        raise ValueError(f"input must be a canonical regular absolute file: {path}")
    data = path.read_bytes()
    if not data and not allow_empty:
        raise ValueError(f"empty evidence file: {path}")
    return path, data


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key in evidence")
        result[key] = value
    return result


def _option(argv, name):
    values = [argv[index + 1] for index, arg in enumerate(argv[:-1]) if arg == name]
    values += [arg[len(name) + 1:] for arg in argv if arg.startswith(name + "=")]
    if len(values) != 1:
        raise ValueError(f"wrapper must bind exactly one {name}")
    return values[0]


def _ids(values):
    if not isinstance(values, list) or any(not isinstance(value, str) or not ID.fullmatch(value) for value in values):
        raise ValueError("invalid recorded scheduler IDs")
    if len(values) != len(set(values)):
        raise ValueError("duplicate recorded scheduler IDs")
    return values


def _gate_evidence(*, gate_manifest, inspect, profile, run_root):
    root = Path(run_root)
    if (not root.is_absolute() or not root.is_dir() or root.resolve() != root
            or root.parent.name != "_ecoda_runs" or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", root.name)):
        raise ValueError("run root must be canonical scratch/_ecoda_runs/RUN_ID")
    paths = {}
    for name, path in (("gate_manifest", gate_manifest), ("inspect", inspect), ("profile", profile)):
        paths[name] = _file(path)
    gate, inspected, configured = [json.loads(paths[name][1], object_pairs_hook=_object)
                                 for name in ("gate_manifest", "inspect", "profile")]
    if (gate["state"] != "COMPLETED" or gate["audit_state"] != "COMPLETED"
            or gate["audit_completed"] is not True or gate["audit"]["passed"] is not True
            or inspected["operation"] != "inspect" or inspected["state"] != "COMPLETED"
            or inspected["gate_id"] != gate["gate_id"] or inspected["audit"] != gate["audit"]):
        raise ValueError("gate and first inspect do not bind the same completed audit")
    if (configured["name"] != gate["profile"] or _sha(paths["profile"][1]) != gate["profile_digest"]
            or configured["policy"]["accounting_command"] != ACCOUNTING_COMMAND):
        raise ValueError("profile identity/digest or accounting command differs from this bridge")
    if gate["command_digest"] != _sha(gate["exact_command"].encode("utf-8")):
        raise ValueError("wrapper command differs from its recorded digest")
    wrapper = shlex.split(gate["exact_command"])
    if (_option(wrapper, "--run-id") != root.name
            or Path(_option(wrapper, "--scratch-root")) != root.parent.parent):
        raise ValueError("wrapper does not bind this exact child run root")
    arrays = _ids(gate["scheduler_array_ids"])
    watchdogs = _ids(gate["scheduler_watchdog_ids"])
    requested = list(dict.fromkeys(arrays + watchdogs))
    if not requested or len(requested) != len(arrays) + len(watchdogs):
        raise ValueError("missing or multiply classified gate scheduler IDs")
    accounting = inspected["audit"]["accounting"]
    expected_command = ACCOUNTING_COMMAND.replace("{scheduler_ids}", shlex.quote(",".join(requested)))
    if (accounting["attempted"] is not True or type(accounting["query_count"]) is not int
            or accounting["query_count"] != 1 or type(accounting["actual_exit"]) is not int
            or accounting["actual_exit"] != 0 or accounting["scheduler_ids"] != requested
            or accounting["command"] != expected_command):
        raise ValueError("accounting is not exactly one successful query of the recorded scope")
    stdout = accounting.get("stdout")
    if (not isinstance(stdout, str) or not stdout or len(stdout) > 100_000
            or "[output truncated]" in stdout):
        raise ValueError("complete original accounting stdout is required")
    rows, parsed = {}, []
    for number, line in enumerate(stdout.splitlines(keepends=True), 1):
        if not line.endswith("\n"):
            raise ValueError("unterminated accounting row")
        fields = line[:-1].split("|")
        if len(fields) != 3 or not ID.fullmatch(fields[0]) or fields[0] in rows:
            raise ValueError("malformed or duplicate raw accounting row")
        if fields[1:] != ["COMPLETED", "0:0"]:
            raise ValueError("every raw accounting row, including unmatched children, must succeed")
        rows[fields[0]] = line
        parsed.append({"line": number, "job_id": fields[0], "state": fields[1], "exit_code": fields[2]})
    if not set(requested).issubset(rows):
        raise ValueError("accounting is missing a requested root row")
    if accounting["rows"] != parsed:
        raise ValueError("raw accounting differs from inspect parsed evidence")
    return root, paths, gate, requested, rows, stdout.encode("utf-8"), expected_command


def _scheduler_rows(data):
    rows = []
    for line in data.decode("utf-8").splitlines():
        fields = line.split("\t")
        if (len(fields) != 2 or not re.fullmatch(r"[A-Za-z0-9_.-]+", fields[0])
                or not ID.fullmatch(fields[1]) or fields[1] in [row[1] for row in rows]):
            raise ValueError("malformed or duplicate child scheduler ID")
        rows.append(fields)
    return rows


def _projection_proof(root, paths, gate, requested, rows, raw, expected_command, *, raw_path=None):
    projection = "".join(line for identifier, line in rows.items() if identifier in requested).encode("utf-8")
    folder = root / "manifests"
    if folder.resolve() != folder:
        raise ValueError("unsafe manifest directory")
    raw_path = raw_path or folder / "scheduler_accounting.raw.tsv"
    projection_path = folder / "scheduler_accounting.tsv"
    proof_path = folder / "scheduler_accounting.proof.json"
    proof = {"schema": "ecoda_gate_accounting_projection_v1", "role": "DERIVED_ROOT_PROJECTION",
             "gate_id": gate["gate_id"], "run_root": str(root), "query_command": expected_command,
             "query_count": 1, "requested_ids": requested, "all_raw_rows_success": True,
             "raw_row_count": len(rows), "projection_row_count": len(requested),
             "inputs": {name: {"path": str(path), "sha256": _sha(data)} for name, (path, data) in paths.items()},
             "raw": {"path": str(raw_path), "sha256": _sha(raw)},
             "projection": {"path": str(projection_path), "sha256": _sha(projection)},
             "rule": "Only unchanged actual requested-root rows; every original raw row validated before projection"}
    return proof, [(proof_path, _json_bytes(proof)), (projection_path, projection)]


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def _publish(payloads, paths):
    outputs = []
    for path, data in payloads:
        sidecar = f"MD5={hashlib.md5(data).hexdigest()}\nSIZE={len(data)}\nPATH={path}\n".encode()
        outputs.extend([(Path(f"{path}.md5"), sidecar), (path, data)])
    if any(path.exists() or path.is_symlink() for path, _ in outputs):
        raise ValueError("accounting evidence already exists; preserve it without replacement")
    if any(path.read_bytes() != data for path, data in paths.values()):
        raise ValueError("input evidence changed during validation")
    # Publish each complete file create-only; projection appears last. Preserve
    # any published partial evidence on failure, removing only our temp files.
    for path, data in outputs:
        fd, temporary = tempfile.mkstemp(prefix=".accounting-capture-", dir=path.parent)
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            os.link(temporary, path)
        finally:
            os.unlink(temporary)


def create_capture(*, gate_manifest, inspect, profile, scheduler_ids, run_root):
    root, paths, gate, requested, rows, raw, command = _gate_evidence(
        gate_manifest=gate_manifest, inspect=inspect, profile=profile, run_root=run_root)
    paths["scheduler_ids"] = _file(scheduler_ids)
    if paths["scheduler_ids"][0] != root / "manifests/scheduler_ids.tsv":
        raise ValueError("scheduler IDs must be the exact child run-owned manifest")
    if {row[1] for row in _scheduler_rows(paths["scheduler_ids"][1])} != set(requested):
        raise ValueError("child scheduler IDs differ from the exact gate union")
    proof, projections = _projection_proof(root, paths, gate, requested, rows, raw, command)
    _publish([(root / "manifests/scheduler_accounting.raw.tsv", raw)] + projections, paths)
    return proof


def _fields(data):
    result = {}
    for line in data.decode("utf-8").splitlines():
        key, separator, value = line.partition("=")
        if not separator or not key or key in result:
            raise ValueError("malformed or duplicate evidence field")
        result[key] = value
    return result


def _receipt(paths, path, *, producer=None, allow_empty=False):
    path = Path(path)
    data = _remember(paths, path, allow_empty=allow_empty)
    checksum = _remember(paths, f"{path}.md5")
    expected = {"PATH": str(path), "SIZE": str(len(data)), "MD5": hashlib.md5(data).hexdigest()}
    if _fields(checksum) != expected:
        raise ValueError(f"receipt checksum differs: {path}")
    if producer:
        root = path.parent.parent
        record = root / "manifests/artifacts" / f"{_sha(str(path).encode())[:32]}.record"
        contents = _remember(paths, record)
        if (_fields(contents) != dict(expected, RUN_ID=root.name, PRODUCER=producer, STATE="PUBLISHED")
                or path.stat().st_mode & 0o222):
            raise ValueError(f"sealed producer receipt differs: {path}")
    return data


def _remember(paths, path, *, allow_empty=False):
    path, data = _file(path, allow_empty=allow_empty)
    if str(path) in paths and paths[str(path)][1] != data:
        raise ValueError("input evidence changed between reads")
    paths[str(path)] = (path, data)
    return data


def create_campaign_capture(*, gate_manifest, inspect, profile, run_root):
    parent, paths, gate, requested, rows, raw, command = _gate_evidence(
        gate_manifest=gate_manifest, inspect=inspect, profile=profile, run_root=run_root)
    folder = parent / "manifests"
    wrapper = shlex.split(gate["exact_command"])
    if _option(wrapper, "--script") != "src/5_run_benchmark_methods/stage5_campaign.sh":
        raise ValueError("campaign mode requires the exact campaign wrapper")
    release = _option(wrapper, "--derivative-release")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", release) or release in (".", ".."):
        raise ValueError("invalid campaign release")
    for name in ("source.manifest", "runtime.identity"):
        _receipt(paths, folder / name, producer="stage5_campaign_identity")
    source = _fields(paths[str(folder / "source.manifest")][1])
    runtime = _fields(paths[str(folder / "runtime.identity")][1])
    if (source["SOURCE_ROOT"] != _option(wrapper, "--source-root")
            or source["AUX_ROOT"] != str(Path(source["SOURCE_ROOT"]) / "aux")
            or _remember(paths, _option(wrapper, "--source-manifest")) != paths[str(folder / "source.manifest")][1]
            or runtime["RUNTIME_IMAGE"] != _option(wrapper, "--runtime-image")
            or runtime["RUNTIME_MANIFEST"] != _option(wrapper, "--runtime-manifest")):
        raise ValueError("campaign source/runtime/aux identities differ from the exact wrapper")
    configuration = Path(source["SOURCE_ROOT"]) / "datasets.json"
    config_bytes = _remember(paths, configuration)
    if _sha(config_bytes) != source["DATASETS_SHA256"]:
        raise ValueError("campaign source configuration differs")
    contexts = _receipt(paths, folder / "campaign_contexts.tsv", producer="stage5_campaign_scope")
    identities = _receipt(paths, folder / "campaign_identities.tsv", producer="stage5_campaign_scope")
    original = _receipt(paths, folder / "campaign_requested.tsv", producer="stage5_campaign_scope")
    if _remember(paths, _option(wrapper, "--contexts")) != original:
        raise ValueError("campaign request differs from the exact wrapper input")
    context_rows = [line.split("\t") for line in contexts.decode().splitlines()]
    original_rows = [line.split("\t") for line in original.decode().splitlines()]
    expected_views = {"benchmark": "benchmark_analysis", "uncorrected": "batch_effect_uncorrected",
                      "corrected": "batch_effect_corrected"}
    if (len(context_rows) != 3 or any(len(row) != 4 for row in context_rows)
            or {row[0] for row in context_rows} != set(expected_views)
            or len({row[1] for row in context_rows}) != 3 or len(original_rows) != 3
            or any(len(row) != 4 for row in original_rows)
            or [(row[0], row[1]) for row in original_rows] != [(row[0], row[1]) for row in context_rows]):
        raise ValueError("campaign requires the three exact unique context/child rows")
    sealed_hashes = {}
    for line in identities.decode().splitlines():
        fields = line.split("\t")
        if len(fields) != 2 or fields[0] in sealed_hashes or not re.fullmatch(r"[0-9a-f]{64}", fields[1]):
            raise ValueError("invalid campaign scope identities")
        sealed_hashes[fields[0]] = fields[1]
    scheduler_capture = _receipt(paths, folder / "campaign_scheduler_ids.tsv")
    declared = [line.split("\t") for line in scheduler_capture.decode().splitlines()]
    actual, children, scope_paths = [], [], set()
    parent_terminal = _fields(_remember(paths, parent / "status/terminal"))
    if parent_terminal.get("RUN_ID") != parent.name or parent_terminal.get("STATE") != "OK":
        raise ValueError("campaign parent terminal receipt is not successful")
    for index, (context, child, selection_text, references_text) in enumerate(context_rows):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", child) or child == parent.name:
            raise ValueError("unsafe campaign child ID")
        root = parent.parent / child
        if not root.is_dir() or root.resolve() != root:
            raise ValueError("campaign child root is not canonical")
        selection_path = folder / f"campaign_{context}.selection.tsv"
        references_path = folder / f"campaign_{context}.inputs.tsv"
        if selection_text != str(selection_path) or references_text != str(references_path):
            raise ValueError("campaign scopes must use their exact sealed parent paths")
        selections = _receipt(paths, selection_path, producer="stage5_campaign_scope")
        references = _receipt(paths, references_path, producer="stage5_campaign_scope")
        for path, data, original_text in ((selection_path, selections, original_rows[index][2]),
                                          (references_path, references, original_rows[index][3])):
            scope_paths.add(str(path))
            if sealed_hashes.get(str(path)) != _sha(data) or _remember(paths, original_text) != data:
                raise ValueError("campaign sealed scope/input identity differs")
        selection_rows = [line.split("\t") for line in selections.decode().splitlines()]
        if not selection_rows or any(len(row) not in (3, 4) or row[1] != expected_views[context] for row in selection_rows):
            raise ValueError("campaign selection has the wrong context view")
        child_folder = root / "manifests"
        for name in ("source.manifest", "runtime.identity"):
            if _remember(paths, child_folder / name) != paths[str(folder / name)][1]:
                raise ValueError("child source/runtime identity differs from campaign")
        child_selection = child_folder / ("requested_selection.tsv" if context == "benchmark" else "selection.tsv")
        child_scope = _receipt(paths, child_selection)
        if (child_scope != selections or
                _remember(paths, child_folder / "input_artifacts.tsv") != references):
            raise ValueError("child scientific request/input references differ from sealed campaign scope")
        binding_bytes = _receipt(paths, child_folder / "artifact_bindings.json")
        binding = json.loads(binding_bytes, object_pairs_hook=_object)
        if (binding["schema"] != "ecoda_derivative_artifacts_v1" or binding["stage"] != "stage5"
                or binding["run_id"] != child or binding["run_root"] != str(root)
                or binding["release_id"] != release or binding["input_artifacts"] != str(child_folder / "input_artifacts.tsv")
                or binding["selection"] != str(child_folder / "selection.tsv")):
            raise ValueError("child derivative binding has a different campaign identity")
        for name in ("source.manifest", "runtime.identity", "input_artifacts.tsv", "selection.tsv"):
            path = child_folder / name
            if binding["identities"].get(str(path)) != _sha(_remember(paths, path)):
                raise ValueError("child derivative binding input identity differs")
        if binding["identities"].get(str(configuration)) != _sha(config_bytes):
            raise ValueError("child derivative configuration differs")
        if context == "benchmark" and (binding.get("requested_selection") != str(child_selection)
                or binding["identities"].get(str(child_selection)) != _sha(child_scope)):
            raise ValueError("child benchmark requested-scope binding differs")
        terminal = _fields(_remember(paths, root / "status/terminal"))
        receipt = _fields(_remember(paths, parent / "status" / context))
        scheduler_path = child_folder / "scheduler_ids.tsv"
        if (terminal.get("RUN_ID") != child or terminal.get("STATE") not in ("OK", "NOOP_VALIDATED")
                or receipt != {"CONTEXT": context, "RUN_ID": child, "EXIT_CODE": "0", "SCHEDULER_MANIFEST": str(scheduler_path)}):
            raise ValueError("child/campaign terminal receipt is not successful")
        child_schedulers = _scheduler_rows(_receipt(paths, scheduler_path, allow_empty=True))
        if not child_schedulers and terminal["STATE"] != "NOOP_VALIDATED":
            raise ValueError("only an explicit NOOP child may lack scheduler IDs")
        actual.extend([[context, child, kind, identifier] for kind, identifier in child_schedulers])
        children.append((context, root, [row[1] for row in child_schedulers]))
    union = [row[3] for row in actual]
    if (set(sealed_hashes) != scope_paths or declared != actual or len(union) != len(set(union))
            or set(union) != set(requested)):
        raise ValueError("campaign scheduler union or sealed identity scope differs")
    raw_path = folder / "scheduler_accounting.raw.tsv"
    payloads = [(raw_path, raw)]
    summaries = []
    for context, root, child_ids in children:
        proof, outputs = _projection_proof(root, paths, gate, child_ids, rows, raw, command, raw_path=raw_path)
        proof.update(campaign_run_root=str(parent), context=context, gate_requested_ids=requested)
        outputs[0] = (outputs[0][0], _json_bytes(proof))
        summaries.append({"context": context, "run_root": str(root), "requested_ids": child_ids,
                          "proof": str(outputs[0][0]), "projection": proof["projection"]})
        payloads.extend(outputs)
    proof = {"schema": "ecoda_campaign_accounting_projection_v1", "role": "DERIVED_CAMPAIGN_CHILD_PROJECTIONS",
             "gate_id": gate["gate_id"], "run_root": str(parent), "query_command": command, "query_count": 1,
             "requested_ids": requested, "all_raw_rows_success": True, "raw_row_count": len(rows),
             "raw": {"path": str(raw_path), "sha256": _sha(raw)}, "children": summaries,
             "inputs": {name: {"path": str(path), "sha256": _sha(data)} for name, (path, data) in paths.items()},
             "audit_rule": "Run all three child audits with their projections; never an ordinary parent audit"}
    payloads.append((folder / "scheduler_accounting.proof.json", _json_bytes(proof)))
    _publish(payloads, paths)
    return proof


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gate-manifest", "inspect", "profile", "run-root"):
        parser.add_argument(f"--{name}", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--scheduler-ids")
    mode.add_argument("--campaign", action="store_true")
    try:
        args = vars(parser.parse_args())
        campaign = args.pop("campaign")
        if campaign:
            args.pop("scheduler_ids")
            proof = create_campaign_capture(**args)
        else:
            proof = create_capture(**args)
    except (ValueError, KeyError, TypeError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(json.dumps(proof, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
