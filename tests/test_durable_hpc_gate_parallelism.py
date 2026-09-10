#!/usr/bin/env python3
"""Standalone synthetic regressions for durable gate dependency parallelism."""
from __future__ import annotations

import importlib.util
import json
import os
import shlex
import signal
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable


GATE_SCRIPT = (
    Path.home()
    / ".agents"
    / "skills"
    / "durable-hpc-gate"
    / "scripts"
    / "durable_hpc_gate.py"
)


def _load_gate_module() -> Any:
    spec = importlib.util.spec_from_file_location("durable_hpc_gate_parallelism", GATE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {GATE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_profile(path: Path) -> None:
    profile = {
        "name": "synthetic-parallelism",
        "execution": {
            "mode": "synthetic",
            "runner_shell": "/bin/bash",
            "ssh_binary": None,
            "tmux_binary": None,
        },
        "policy": {
            "invariants": [
                {"name": "workdir-exists", "command": "test -d \"$PWD\"", "expected_exit": 0}
            ],
            "audit_commands": [
                {"name": "synthetic-audit", "command": "true", "expected_exit": 0}
            ],
            "artifact_contracts": [
                {"name": "synthetic-artifact", "command": "true", "expected_exit": 0}
            ],
            "immutable_fingerprints": [],
            "accounting_command": None,
            "require_scheduler_ids": False,
            "reviewer_required": True,
            "manifest_constraints": {},
            "known_discrepancies": [],
        },
    }
    path.write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")


def _invoke(module: Any, output_dir: Path, operation: str, *arguments: str) -> tuple[int, dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{operation}-{time.monotonic_ns()}.json"
    rc = module.main([operation, *arguments, "--output", str(output)])
    if not output.exists():
        raise AssertionError(f"{operation} did not write structured output (rc={rc})")
    return rc, json.loads(output.read_text(encoding="utf-8"))


def _prepare(
    module: Any,
    root: Path,
    profile: Path,
    gate_id: str,
    *,
    project: str = "ecoda",
    serialization_group: str,
    exact_command: str,
    dependencies: Iterable[Path] = (),
) -> Path:
    gate_root = root / gate_id
    workdir = gate_root / "work"
    workdir.mkdir(parents=True, exist_ok=True)
    manifest_root = root / "manifests"
    manifest_root.mkdir(parents=True, exist_ok=True)
    manifest = manifest_root / f"{gate_id}.json"
    remote = gate_root / "remote"
    args = [
        "--manifest",
        str(manifest),
        "--profile",
        str(profile),
        "--project",
        project,
        "--gate-id",
        gate_id,
        "--remote-host",
        "synthetic-host",
        "--remote-workdir",
        str(workdir),
        "--exact-command",
        exact_command,
        "--serialization-group",
        serialization_group,
        "--tmux-session",
        f"{gate_id}-tmux",
        "--completion-channel",
        f"file:{remote / 'completion.marker'}",
        "--remote-manifest",
        str(remote / "manifest.json"),
        "--remote-runner",
        str(remote / "runner.sh"),
        "--remote-log",
        str(remote / "runner.log"),
        "--remote-status",
        str(remote / "status.json"),
    ]
    for dependency in dependencies:
        args.extend(("--dependency-manifest", str(dependency)))
    rc, payload = _invoke(module, root / "results", "prepare", *args)
    assert rc == 0, payload
    assert payload["ok"] is True
    assert payload["state"] == "PREPARED"
    return manifest


def _launch(module: Any, root: Path, profile: Path, manifest: Path) -> tuple[int, dict[str, Any]]:
    return _invoke(
        module,
        root / "results",
        "launch",
        "--manifest",
        str(manifest),
        "--profile",
        str(profile),
    )


def _wait(module: Any, root: Path, manifest: Path) -> tuple[int, dict[str, Any]]:
    return _invoke(
        module,
        root / "results",
        "wait",
        "--manifest",
        str(manifest),
    )


def _inspect(
    module: Any,
    root: Path,
    profile: Path,
    manifest: Path,
    *extra: str,
) -> tuple[int, dict[str, Any]]:
    return _invoke(
        module,
        root / "results",
        "inspect",
        "--manifest",
        str(manifest),
        "--profile",
        str(profile),
        *extra,
    )


def _wait_for_lines(path: Path, expected: list[str]) -> None:
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline:
        try:
            lines = path.read_text(encoding="utf-8").splitlines()
        except FileNotFoundError:
            lines = []
        if lines == expected:
            return
        time.sleep(0.02)
    actual = path.read_text(encoding="utf-8").splitlines() if path.exists() else []
    raise AssertionError(f"{path} did not reach expected lines {expected!r}: {actual!r}")

def _read_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _complete_and_review(module: Any, root: Path, profile: Path, manifest: Path) -> None:
    rc, payload = _launch(module, root, profile, manifest)
    assert rc == 0, payload
    assert payload["state"] == "RUNNING"
    rc, payload = _wait(module, root, manifest)
    assert rc == 0, payload
    assert payload["state"] == "COMPLETED"
    rc, payload = _inspect(module, root, profile, manifest)
    assert rc == 0, payload
    assert payload["audit"]["passed"] is True
    rc, payload = _inspect(module, root, profile, manifest, "--approve-reviewer", "reviewer")
    assert rc == 0, payload
    completed = _read_manifest(manifest)
    assert completed["state"] == "COMPLETED"
    assert completed["audit_state"] == "COMPLETED"
    assert completed["audit_completed"] is True
    assert completed["audit"]["passed"] is True
    assert completed["immutable_fingerprints"] == []
    assert completed["audit"]["immutable_fingerprints"] == []
    assert completed["reviewer"]["approved"] is True


def _runner_pids(manifests: Iterable[Path]) -> list[int]:
    pids: list[int] = []
    for manifest in manifests:
        try:
            value = _read_manifest(manifest)
        except (OSError, json.JSONDecodeError):
            continue
        for evidence in reversed(value.get("launcher_evidence", [])):
            if not isinstance(evidence, dict) or evidence.get("kind") != "tmux_runner":
                continue
            pid = evidence.get("pid")
            if isinstance(pid, int) and not isinstance(pid, bool) and pid > 0:
                pids.append(pid)
            break
    return pids


def _cleanup_synthetic_processes(manifests: Iterable[Path]) -> None:
    for pid in _runner_pids(manifests):
        try:
            os.killpg(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            pass
        deadline = time.monotonic() + 2.0
        while time.monotonic() < deadline:
            try:
                waited, _status = os.waitpid(pid, os.WNOHANG)
            except (ChildProcessError, ProcessLookupError):
                break
            if waited == pid:
                break
            time.sleep(0.02)
        else:
            try:
                os.killpg(pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            try:
                os.waitpid(pid, 0)
            except (ChildProcessError, ProcessLookupError):
                pass


def _assert_dependency_rejected(
    module: Any,
    root: Path,
    profile: Path,
    manifest: Path,
    *,
    code: str,
) -> None:
    rc, payload = _launch(module, root, profile, manifest)
    assert rc == 1, payload
    assert payload["ok"] is False
    assert payload["code"] == code
    rejected = _read_manifest(manifest)
    assert rejected["state"] == "PRELAUNCH_STOP"
    assert any(item["kind"] == "dependency_blocked" for item in rejected["discrepancies"])
    assert any(
        history["state"] == "PRELAUNCH_STOP" for history in rejected["state_history"]
    )


def _write_overlap_fixture(root: Path) -> dict[str, Path]:
    """Create two distinct Stage 3 selections that share one output path."""
    scratch = root / "ownership-scratch"
    nas = root / "ownership-nas"
    datasets = root / "ownership-datasets.json"
    selection_one = root / "selection-one.tsv"
    selection_two = root / "selection-two.tsv"
    sbatch_log = root / "fake-sbatch.calls"
    fake_sbatch = root / "bin" / "sbatch"
    fake_sbatch.parent.mkdir(parents=True, exist_ok=True)
    for base in (scratch, nas):
        for dataset in ("dataset-one", "dataset-two"):
            (base / dataset / "output").mkdir(parents=True, exist_ok=True)
    datasets.write_text(
        json.dumps(
            {
                "dataset-one": {
                    "views": {
                        "benchmark_analysis": {"output_file_name": "shared.h5ad"}
                    }
                },
                "dataset-two": {
                    "views": {
                        "benchmark_analysis": {"output_file_name": "other.h5ad"}
                    }
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    selection_one.write_text("dataset-one\tbenchmark_analysis\n", encoding="utf-8")
    selection_two.write_text(
        "dataset-one\tbenchmark_analysis\n"
        "dataset-two\tbenchmark_analysis\n",
        encoding="utf-8",
    )
    fake_sbatch.write_text(
        "#!/bin/sh\n"
        'printf "%s\\n" "$1" >> "${ECODA_FAKE_SBATCH_LOG:?}"\n',
        encoding="utf-8",
    )
    fake_sbatch.chmod(0o755)
    return {
        "scratch": scratch,
        "nas": nas,
        "datasets": datasets,
        "selection_one": selection_one,
        "selection_two": selection_two,
        "sbatch_log": sbatch_log,
        "fake_sbatch": fake_sbatch,
    }


def _ownership_command(
    fixture: dict[str, Path], selection: Path, run_id: str
) -> str:
    common = Path(__file__).resolve().parents[1] / "src/utils/bash/ecoda_run_common.sh"
    exports = " ".join(
        (
            f"HPC_SCRATCH_DIR={shlex.quote(str(fixture['scratch']))}",
            f"NAS_TARGET_DIR={shlex.quote(str(fixture['nas']))}",
            f"DATASETS_JSON_FILE={shlex.quote(str(fixture['datasets']))}",
            f"ECODA_FAKE_SBATCH_LOG={shlex.quote(str(fixture['sbatch_log']))}",
        )
    )
    return (
        f"export {exports}; "
        f"source {shlex.quote(str(common))} && "
        f"ecoda_validate_output_ownership stage3 {shlex.quote(str(selection))} "
        f"{shlex.quote(run_id)} && "
        f"{shlex.quote(str(fixture['fake_sbatch']))} {shlex.quote(run_id)} && "
        "sleep 60"
    )


def test_durable_gate_parallelism_contract() -> None:
    module = _load_gate_module()
    manifests: list[Path] = []
    with tempfile.TemporaryDirectory(prefix="ecoda-durable-parallelism-") as raw:
        root = Path(raw)
        profile = root / "synthetic-profile.json"
        _write_profile(profile)
        counter = root / "shared-invocations.log"
        sleep_command = "sleep 60"
        try:
            predecessor = _prepare(
                module,
                root,
                profile,
                "reviewed-predecessor",
                serialization_group="predecessor-wave",
                exact_command="true",
            )
            manifests.append(predecessor)
            _complete_and_review(module, root, profile, predecessor)
            predecessor_value = _read_manifest(predecessor)
            assert predecessor_value["serialization_group"] == "predecessor-wave"

            b4 = _prepare(
                module,
                root,
                profile,
                "B4",
                serialization_group="b4-wave",
                exact_command=sleep_command,
                dependencies=(predecessor,),
            )
            benchmark_p3 = _prepare(
                module,
                root,
                profile,
                "benchmark-P3",
                serialization_group="benchmark-p3-wave",
                exact_command=sleep_command,
                dependencies=(predecessor,),
            )
            manifests.extend((b4, benchmark_p3))
            assert _read_manifest(b4)["dependency_manifests"] == [str(predecessor.resolve())]
            assert _read_manifest(benchmark_p3)["dependency_manifests"] == [str(predecessor.resolve())]
            assert _read_manifest(b4)["serialization_group"] != _read_manifest(benchmark_p3)["serialization_group"]

            rc, payload = _launch(module, root, profile, b4)
            assert rc == 0, payload
            assert payload["ok"] is True
            assert payload["state"] == "RUNNING"
            rc, payload = _launch(module, root, profile, benchmark_p3)
            assert rc == 0, payload
            assert payload["ok"] is True
            assert payload["state"] == "RUNNING"

            overlap = _write_overlap_fixture(root)
            assert overlap["selection_one"].read_text(encoding="utf-8") != overlap[
                "selection_two"
            ].read_text(encoding="utf-8")
            overlap_first = _prepare(
                module,
                root,
                profile,
                "overlap-first",
                serialization_group="overlap-wave-one",
                exact_command=_ownership_command(
                    overlap, overlap["selection_one"], "overlap-first"
                ),
            )
            overlap_second = _prepare(
                module,
                root,
                profile,
                "overlap-second",
                serialization_group="overlap-wave-two",
                exact_command=_ownership_command(
                    overlap, overlap["selection_two"], "overlap-second"
                ),
            )
            manifests.extend((overlap_first, overlap_second))
            assert _read_manifest(overlap_first)["serialization_group"] != _read_manifest(
                overlap_second
            )["serialization_group"]

            rc, payload = _launch(module, root, profile, overlap_first)
            assert rc == 0, payload
            assert payload["state"] == "RUNNING"
            _wait_for_lines(overlap["sbatch_log"], ["overlap-first"])

            # Distinct durable serialization groups do not bypass concrete
            # artifact ownership.  The second selection overlaps the first
            # expanded scratch/NAS paths and must fail before fake sbatch.
            rc, payload = _launch(module, root, profile, overlap_second)
            assert rc == 0, payload
            assert payload["state"] == "RUNNING"
            rc, payload = _wait(module, root, overlap_second)
            assert rc == 0, payload
            assert payload["state"] == "FAILED"
            overlap_second_value = _read_manifest(overlap_second)
            assert overlap_second_value["state"] == "FAILED"
            assert overlap["sbatch_log"].read_text(encoding="utf-8").splitlines() == [
                "overlap-first"
            ]
            overlap_second_log = Path(overlap_second_value["remote_log"]).read_text(
                encoding="utf-8"
            )
            assert "active global artifact owner" in overlap_second_log
            assert str(
                (overlap["scratch"] / "dataset-one" / "output" / "shared.h5ad").resolve()
            ) in overlap_second_log

            malformed_predecessor = root / "malformed-predecessor.json"
            malformed_predecessor.write_text("{not-json\n", encoding="utf-8")
            malformed_child = _prepare(
                module,
                root,
                profile,
                "malformed-child",
                serialization_group="malformed-wave",
                exact_command=sleep_command,
                dependencies=(malformed_predecessor,),
            )
            manifests.append(malformed_child)
            _assert_dependency_rejected(
                module,
                root,
                profile,
                malformed_child,
                code="invalid_json",
            )

            mismatched_predecessor = _prepare(
                module,
                root,
                profile,
                "mismatched-predecessor",
                project="other-project",
                serialization_group="other-wave",
                exact_command="true",
            )
            manifests.append(mismatched_predecessor)
            _complete_and_review(module, root, profile, mismatched_predecessor)
            mismatched_child = _prepare(
                module,
                root,
                profile,
                "mismatched-child",
                serialization_group="mismatched-wave",
                exact_command=sleep_command,
                dependencies=(mismatched_predecessor,),
            )
            manifests.append(mismatched_child)
            _assert_dependency_rejected(
                module,
                root,
                profile,
                mismatched_child,
                code="dependency_blocked",
            )

            conflict_counter = shlex.quote(str(counter))
            first = _prepare(
                module,
                root,
                profile,
                "shared-first",
                serialization_group="shared-resource-wave",
                exact_command=f"printf 'first\\n' >> {conflict_counter}; sleep 60",
            )
            second = _prepare(
                module,
                root,
                profile,
                "shared-second",
                serialization_group="shared-resource-wave",
                exact_command=f"printf 'second\\n' >> {conflict_counter}; sleep 60",
            )
            manifests.extend((first, second))
            rc, payload = _launch(module, root, profile, first)
            assert rc == 0, payload
            assert payload["ok"] is True
            assert payload["state"] == "RUNNING"
            _wait_for_lines(counter, ["first"])

            rc, payload = _launch(module, root, profile, second)
            assert rc == 1, payload
            assert payload["ok"] is False
            assert payload["code"] == "resource_lock_conflict"
            rejected = _read_manifest(second)
            assert rejected["state"] == "PRELAUNCH_STOP"
            assert not any(
                evidence.get("kind") == "tmux_runner"
                for evidence in rejected["launcher_evidence"]
            )
            assert any(item["kind"] == "resource_lock_conflict" for item in rejected["discrepancies"])
            assert any(
                history["state"] == "PRELAUNCH_STOP" for history in rejected["state_history"]
            )
            assert counter.read_text(encoding="utf-8").splitlines() == ["first"]
        finally:
            _cleanup_synthetic_processes(manifests)


def main() -> None:
    test_durable_gate_parallelism_contract()
    print("durable gate parallelism: OK")


if __name__ == "__main__":
    main()
