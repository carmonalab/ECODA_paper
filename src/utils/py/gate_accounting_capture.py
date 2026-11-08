"""Bridge one completed gate inspect to the run audit without another query.

The raw JobIDRaw output remains unchanged. The audit input is explicitly a
derived projection of actual requested-root rows, after *every* raw row passes.
This single-child bridge requires exact equality of gate and child scheduler IDs.
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


def _file(path):
    path = Path(path)
    if not path.is_absolute() or path.is_symlink() or not path.is_file() or path.resolve() != path:
        raise ValueError(f"input must be a canonical regular absolute file: {path}")
    data = path.read_bytes()
    if not data:
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


def create_capture(*, gate_manifest, inspect, profile, scheduler_ids, run_root):
    root = Path(run_root)
    if (not root.is_absolute() or not root.is_dir() or root.resolve() != root
            or root.parent.name != "_ecoda_runs" or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", root.name)):
        raise ValueError("run root must be canonical scratch/_ecoda_runs/RUN_ID")
    paths = {}
    for name, path in (("gate_manifest", gate_manifest), ("inspect", inspect),
                       ("profile", profile), ("scheduler_ids", scheduler_ids)):
        paths[name] = _file(path)
    if paths["scheduler_ids"][0] != root / "manifests/scheduler_ids.tsv":
        raise ValueError("scheduler IDs must be the exact child run-owned manifest")
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
    child_ids = []
    for line in paths["scheduler_ids"][1].decode("utf-8").splitlines():
        fields = line.split("\t")
        if (len(fields) != 2 or not re.fullmatch(r"[A-Za-z0-9_.-]+", fields[0])
                or not ID.fullmatch(fields[1]) or fields[1] in child_ids):
            raise ValueError("malformed or duplicate child scheduler ID")
        child_ids.append(fields[1])
    if set(child_ids) != set(requested):
        raise ValueError("child scheduler IDs differ from the exact gate union")
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
    rows, roots, parsed = {}, [], []
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
        if fields[0] in requested:
            roots.append(line)
    if not set(requested).issubset(rows):
        raise ValueError("accounting is missing a requested root row")
    if accounting["rows"] != parsed:
        raise ValueError("raw accounting differs from inspect parsed evidence")
    raw, projection = stdout.encode("utf-8"), "".join(roots).encode("utf-8")
    folder = root / "manifests"
    if folder.resolve() != folder:
        raise ValueError("unsafe manifest directory")
    raw_path = folder / "scheduler_accounting.raw.tsv"
    projection_path = folder / "scheduler_accounting.tsv"
    proof_path = folder / "scheduler_accounting.proof.json"
    proof = {"schema": "ecoda_gate_accounting_projection_v1", "role": "DERIVED_ROOT_PROJECTION",
             "gate_id": gate["gate_id"], "run_root": str(root), "query_command": expected_command,
             "query_count": 1, "requested_ids": requested, "all_raw_rows_success": True,
             "raw_row_count": len(rows), "projection_row_count": len(roots),
             "inputs": {name: {"path": str(path), "sha256": _sha(data)} for name, (path, data) in paths.items()},
             "raw": {"path": str(raw_path), "sha256": _sha(raw)},
             "projection": {"path": str(projection_path), "sha256": _sha(projection)},
             "rule": "Only unchanged actual requested-root rows; every original raw row validated before projection"}
    payloads = [(raw_path, raw), (proof_path, (json.dumps(proof, sort_keys=True, indent=2) + "\n").encode()),
                (projection_path, projection)]
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
        fd, temporary = tempfile.mkstemp(prefix=".accounting-capture-", dir=folder)
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            os.link(temporary, path)
        finally:
            os.unlink(temporary)
    return proof


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gate-manifest", "inspect", "profile", "scheduler-ids", "run-root"):
        parser.add_argument(f"--{name}", required=True)
    try:
        proof = create_capture(**vars(parser.parse_args()))
    except (ValueError, KeyError, TypeError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(json.dumps(proof, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
