"""Run-bound derivative destinations and read-only input artifact references.

Input TSV: dataset, view, absolute input path, MD5, size, checksum-proof path.
Proofs use the existing MD5/SIZE/PATH record format. Binding checks metadata;
compute-node consumers can request content verification before reading inputs.
Stage 4 may seed an explicit create-only working copy; parents stay read-only.
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

SCHEMA = "ecoda_derivative_artifacts_v1"
COMPONENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
CONTEXTS = {"datasets", "benchmark", "batch_effect/uncorrected_final", "batch_effect/corrected_final"}


def digest(path: Path, algorithm: str = "sha256") -> str:
    result = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def component(value: str) -> str:
    if not COMPONENT.fullmatch(value) or value in {".", ".."}:
        raise ValueError(f"invalid path component: {value!r}")
    return value


def absolute(path: str | Path, *, exists: bool = True) -> Path:
    text = str(path)
    result = Path(text)
    if not result.is_absolute() or any(c in text for c in "\t\r\n") or os.path.normpath(text) != text:
        raise ValueError(f"path must be normalized and absolute: {text!r}")
    if result.resolve(strict=exists) != result or result.is_symlink():
        raise ValueError(f"path is not canonical: {text}")
    return result


def proof(path: Path, payload: Path, md5: str, size: int) -> None:
    fields = {}
    for line in path.read_text().splitlines():
        key, separator, value = line.partition("=")
        if separator and key in {"MD5", "SIZE", "PATH"}:
            if key in fields:
                raise ValueError(f"duplicate checksum-proof field: {path}")
            fields[key] = value
    if fields != {"MD5": md5, "SIZE": str(size), "PATH": str(payload)}:
        raise ValueError(f"input checksum proof differs: {path}")


def configured_files(config: dict, dataset: str, view: str, stage: str) -> list[str]:
    component(dataset); component(view)
    entry = config.get(dataset, {})
    spec = entry.get("views", {}).get(view)
    if not isinstance(spec, dict) or "columns" in spec:
        raise ValueError(f"unknown or invalid dataset/view: {dataset}/{view}")
    value = spec.get("output_file_name") if stage in {"stage4", "stage5"} else spec.get("input_file_name")
    if stage == "stage4" and (view != "benchmark_analysis" or entry.get("use_for_benchmark") is not True):
        raise ValueError("Stage 4 derivatives require an eligible benchmark view")
    names = value if isinstance(value, list) else [value]
    if not names or any(not isinstance(name, str) or component(name) != name for name in names):
        raise ValueError(f"missing configured input filenames: {dataset}/{view}")
    return names


def inputs(path: Path, config: dict, stage: str, *, verify_content: bool = False,
           verify_scope: tuple[str, str] | None = None) -> dict:
    rows = {}
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) != 6 or any(not value for value in parts):
            raise ValueError("input artifact rows require exactly six nonempty TSV fields")
        dataset, view, filename, md5, size_text, proof_text = parts
        payload = absolute(filename)
        checksum_proof = absolute(proof_text)
        names = configured_files(config, dataset, view, stage)
        if payload.name not in names or not payload.is_file() or payload.stat().st_size <= 0:
            raise ValueError(f"input filename or type disagrees with registry: {payload}")
        if not re.fullmatch(r"[0-9a-f]{32}", md5) or not size_text.isdigit() or int(size_text) <= 0:
            raise ValueError("invalid input MD5/size")
        size = int(size_text)
        if str(size) != size_text or payload.stat().st_size != size:
            raise ValueError(f"input artifact size changed: {payload}")
        proof(checksum_proof, payload, md5, size)
        if verify_content and (verify_scope is None or verify_scope == (dataset, view)) and digest(payload, "md5") != md5:
            raise ValueError(f"input artifact checksum changed: {payload}")
        key = (dataset, view, payload.name)
        if key in rows:
            raise ValueError(f"duplicate input artifact: {key}")
        rows[key] = payload
    if not rows:
        raise ValueError("input artifact manifest is empty")
    for dataset, view in {(key[0], key[1]) for key in rows}:
        expected = set(configured_files(config, dataset, view, stage))
        selected = [value for key, value in rows.items() if key[:2] == (dataset, view)]
        if {value.name for value in selected} != expected or len({value.parent for value in selected}) != 1:
            raise ValueError(f"input files require exact configured names and one directory: {dataset}/{view}")
    return rows


def create_only(path: Path, payload: bytes) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("xb") as handle:
            handle.write(payload)
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def bind(*, run_root: Path, release_id: str, nas_root: Path, stage: str,
         config_path: Path, selection: Path, input_artifacts: Path, producer_run: str | None = None) -> Path:
    run_root = absolute(run_root)
    nas_root = absolute(nas_root)
    release_id = component(release_id)
    run_id = component(run_root.name)
    if run_root.parent.name != "_ecoda_runs" or not run_root.is_dir() or not nas_root.is_dir():
        raise ValueError("derivative binding requires an existing exact run root and NAS base")
    if stage not in {"stage2", "stage3", "stage4", "stage5"}:
        raise ValueError("unsupported derivative stage")
    if (stage == "stage4") != (producer_run is not None):
        raise ValueError("Stage 4 requires an explicit upstream producer run only")
    if producer_run is not None:
        component(producer_run)
    config = json.loads(config_path.read_text())
    source = run_root / "manifests/source.manifest"
    runtime = run_root / "manifests/runtime.identity"
    for path in (source, runtime, selection, input_artifacts, config_path):
        absolute(path)
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"missing binding identity: {path}")
    if selection.parent != source.parent:
        raise ValueError("selection must be run-owned")
    selected = set()
    for line in selection.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) not in (2, 3, 4) or any(not field for field in fields):
            raise ValueError("invalid derivative selection row")
        configured_files(config, fields[0], fields[1], stage)
        selected.add(tuple(fields[:2]))
    rows = inputs(input_artifacts, config, stage)
    if not selected or selected != {key[:2] for key in rows}:
        raise ValueError("input artifacts must cover exactly the selected dataset/views")
    copied_inputs = source.parent / "input_artifacts.tsv"
    create_only(copied_inputs, input_artifacts.read_bytes())
    binding = source.parent / "artifact_bindings.json"
    value = dict(schema=SCHEMA, stage=stage, run_id=run_id, run_root=str(run_root),
                 release_id=release_id, nas_base=str(nas_root),
                 artifact_root=str(run_root / "artifacts"),
                 nas_root=str(nas_root / "releases" / release_id / run_id),
                 selection=str(selection), input_artifacts=str(copied_inputs),
                 identities={str(path): digest(path) for path in
                             (source, runtime, selection, copied_inputs, config_path)})
    requested = source.parent / "requested_selection.tsv"
    if stage == "stage4":
        value["input_producer_run"] = producer_run
    if stage == "stage5" and requested.exists():
        absolute(requested)
        if not requested.is_file() or requested.stat().st_size == 0:
            raise ValueError("requested benchmark selection is empty or invalid")
        value["requested_selection"] = str(requested)
        value["identities"][str(requested)] = digest(requested)
    create_only(binding, (json.dumps(value, sort_keys=True) + "\n").encode())
    create_only(Path(f"{binding}.md5"), (
        f"MD5={digest(binding, 'md5')}\nSIZE={binding.stat().st_size}\nPATH={binding}\n"
    ).encode())
    return binding


def load(binding: str | Path, config_path: str | Path, *, verify_content: bool = False) -> tuple[dict, dict]:
    binding = absolute(binding)
    proof(absolute(f"{binding}.md5"), binding, digest(binding, "md5"), binding.stat().st_size)
    value = json.loads(binding.read_text())
    run_root = absolute(value["run_root"])
    # The immutable bind proves the mounted NAS destination on the controller.
    # Compute nodes resolve this bound identity without mounting that filesystem.
    nas_base = absolute(value["nas_base"], exists=False)
    if (value["schema"] != SCHEMA or binding != run_root / "manifests/artifact_bindings.json"
        or run_root.parent.name != "_ecoda_runs" or value["run_id"] != run_root.name
        or value["stage"] not in {"stage2", "stage3", "stage4", "stage5"}
        or value["artifact_root"] != str(run_root / "artifacts")
        or value["nas_root"] != str(nas_base / "releases" / component(value["release_id"]) / run_root.name)):
        raise ValueError("derivative binding root/run identity differs")
    if value["stage"] == "stage4":
        component(value["input_producer_run"])
    expected = {str(run_root / "manifests" / name) for name in
                ("source.manifest", "runtime.identity", "input_artifacts.tsv")}
    expected.update((value["selection"], str(absolute(config_path))))
    if "requested_selection" in value:
        requested = binding.parent / "requested_selection.tsv"
        if value["stage"] != "stage5" or value["requested_selection"] != str(requested):
            raise ValueError("requested benchmark selection is not run-owned")
        expected.add(str(requested))
    if set(value["identities"]) != expected or Path(value["selection"]).parent != binding.parent:
        raise ValueError("derivative binding identity scope differs")
    if value["input_artifacts"] != str(binding.parent / "input_artifacts.tsv"):
        raise ValueError("derivative input manifest is not run-owned")
    for filename, sha256 in value["identities"].items():
        path = absolute(filename)
        if not path.is_file() or digest(path) != sha256:
            raise ValueError(f"derivative binding identity changed: {path}")
    config = json.loads(Path(config_path).read_text())
    rows = inputs(Path(value["input_artifacts"]), config, value["stage"], verify_content=verify_content)
    return value, rows


def resolve_input(binding: str | Path, config_path: str | Path, dataset: str, view: str,
                  *, directory: bool = False, verify_content: bool = False) -> Path:
    value, rows = load(binding, config_path)
    selected = [path for key, path in rows.items() if key[:2] == (dataset, view)]
    if not selected:
        raise ValueError(f"input reference is outside derivative selection: {dataset}/{view}")
    if verify_content:
        inputs(Path(value["input_artifacts"]), json.loads(Path(config_path).read_text()),
               value["stage"], verify_content=True, verify_scope=(dataset, view))
    if directory:
        return selected[0].parent
    if len(selected) != 1:
        raise ValueError("file resolution requires exactly one input artifact")
    return selected[0]


def context_root(value: dict, context: str, *, nas: bool = False) -> Path:
    if context not in CONTEXTS:
        raise ValueError(f"invalid derivative output context: {context}")
    if (value["stage"] in {"stage2", "stage3", "stage4"}) != (context == "datasets"):
        raise ValueError("derivative context disagrees with the selected stage")
    root = Path(value["nas_root"] if nas else value["artifact_root"]) / context
    return absolute(root, exists=False)


def active_binding(config_path: str | Path, *, stage: str = "stage5") -> tuple[dict, Path]:
    """Validate the environment's reference against the exact active run."""
    run_root = absolute(os.environ.get("ECODA_RUN_ROOT", ""))
    binding = run_root / "manifests/artifact_bindings.json"
    if os.environ.get("ECODA_ARTIFACT_BINDING") != str(binding):
        raise ValueError("derivative binding is not owned by the active run")
    value, _ = load(binding, config_path)
    if value["run_root"] != str(run_root) or value["stage"] != stage:
        raise ValueError("derivative binding disagrees with the active run/stage")
    nas = os.environ.get("NAS_TARGET_DIR")
    if nas and value["nas_base"] != str(absolute(nas, exists=False)):
        raise ValueError("derivative binding disagrees with the active NAS base")
    return value, binding


def stage4_input_copy(binding: Path, config_path: Path, dataset: str, view: str) -> Path:
    """Seed a new annotation working copy; the immutable parent stays read-only."""
    value, _ = active_binding(config_path, stage="stage4")
    if binding != Path(value["run_root"]) / "manifests/artifact_bindings.json":
        raise ValueError("annotation copy requires the active binding")
    source = resolve_input(binding, config_path, dataset, view, verify_content=True)
    root = context_root(value, "datasets")
    output = absolute(root / component(dataset) / "output" / source.name, exists=False)
    receipt = Path(value["run_root"]) / "manifests" / f"stage4_input_{component(dataset)}_{component(view)}.json"
    copy_identity = dict(schema="ecoda_stage4_input_copy_v1", binding_sha256=digest(binding),
                         source=str(source), source_md5=digest(source, "md5"), size=source.stat().st_size,
                         output=str(output), input_producer_run=value["input_producer_run"])
    sidecar = Path(f"{output}.md5")
    if any(path.exists() or path.is_symlink() for path in (output, sidecar, receipt)):
        if (not output.is_file() or output.is_symlink() or not receipt.is_file()
                or receipt.is_symlink() or json.loads(receipt.read_text()) != copy_identity
                or digest(output, "md5") != copy_identity["source_md5"]):
            raise ValueError("existing annotation working copy/proof differs; preserved")
        proof(sidecar, output, copy_identity["source_md5"], copy_identity["size"])
        return output
    output.parent.mkdir(parents=True, exist_ok=True)
    absolute(output.parent)
    descriptor, filename = tempfile.mkstemp(prefix=f".{output.name}.copy-", dir=output.parent)
    os.close(descriptor)
    temporary = Path(filename)
    try:
        shutil.copyfile(source, temporary)
        if digest(temporary, "md5") != copy_identity["source_md5"]:
            raise ValueError("annotation input copy checksum differs")
        resolve_input(binding, config_path, dataset, view, verify_content=True)
        create_only(receipt, (json.dumps(copy_identity, sort_keys=True) + "\n").encode())
        create_only(sidecar, (f"MD5={copy_identity['source_md5']}\nSIZE={copy_identity['size']}\nPATH={output}\n").encode())
        os.link(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("bind", "validate", "input", "root", "copy-stage4"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--binding", type=Path)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--release-id")
    parser.add_argument("--nas-root", type=Path)
    parser.add_argument("--stage")
    parser.add_argument("--producer-run")
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--input-artifacts", type=Path)
    parser.add_argument("--dataset")
    parser.add_argument("--view")
    parser.add_argument("--directory", action="store_true")
    parser.add_argument("--verify-content", action="store_true")
    parser.add_argument("--context")
    parser.add_argument("--nas", action="store_true")
    args = parser.parse_args()
    if args.command == "bind":
        print(bind(run_root=args.run_root, release_id=args.release_id, nas_root=args.nas_root,
                   stage=args.stage, config_path=args.config, selection=args.selection,
                   input_artifacts=args.input_artifacts, producer_run=args.producer_run))
    elif args.command == "copy-stage4":
        print(stage4_input_copy(args.binding, args.config, args.dataset, args.view))
    else:
        value, _ = load(args.binding, args.config,
                        verify_content=args.verify_content and args.command != "input")
        if ((args.run_root is not None and value["run_root"] != str(absolute(args.run_root))) or
            (args.nas_root is not None and value["nas_base"] != str(absolute(args.nas_root, exists=False))) or
            (args.stage is not None and value["stage"] != args.stage)):
            raise ValueError("derivative binding disagrees with the active run/stage/NAS base")
        if args.command == "input":
            print(resolve_input(args.binding, args.config, args.dataset, args.view,
                                directory=args.directory, verify_content=args.verify_content))
        else:
            print(context_root(value, args.context, nas=args.nas) if args.command == "root" else "DERIVATIVE_BINDING_VALIDATED")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError) as exc:
        raise SystemExit(f"ERROR: {exc}")
