"""Selected original batch exports copy without caches, fits or export jobs."""
from pathlib import Path
import copy
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("metadata_fixture", ROOT / "tests/test_stage5_prepared_metadata.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
CONTEXT = "batch_effect/corrected_final"
SCHEMA = "ecoda_stage5_metadata_adoption_v1"


def selected_fixture(base):
    f = m.build_fixture(base)
    v = f["value"]
    v["schema"] = SCHEMA
    v["target_run_ids"] = {CONTEXT: v["target_run_ids"][CONTEXT]}
    v["sources"] = {CONTEXT: v["sources"][CONTEXT]}
    v["artifacts"] = [a for a in v["artifacts"] if a["context"] == CONTEXT]
    v["selected_metadata"] = {CONTEXT: [{"dataset": a["dataset"], "view": a["view"]} for a in v["artifacts"]]}
    f["artifacts"] = v["artifacts"]
    m.save_receipt(f)
    return f


def plan(f, valid):
    m.save_receipt(f)
    result = subprocess.run([sys.executable, str(ROOT / "src/utils/py/stage5_prepared_adoption.py"),
        "plan", "--receipt", str(f["receipt"]), "--sha256", m.digest(f["receipt"]),
        "--config", str(f["config"]), "--context", CONTEXT],
        env=f["contexts"][CONTEXT]["env"], text=True, capture_output=True)
    assert (result.returncode == 0) == valid, (result.stdout, result.stderr)
    return result


def main():
    with tempfile.TemporaryDirectory(prefix="ecoda-selected-metadata-") as directory:
        root = Path(directory).resolve()
        f = selected_fixture(root / "valid")
        protected = m.protected_sources(f)
        original = copy.deepcopy(f["value"])
        result = plan(f, True)
        assert len(result.stdout.splitlines()) == 8
        assert all(line.startswith("metadata\t") for line in result.stdout.splitlines())
        # Receipt scope cannot widen or shrink the exact current binding union.
        for fault in ("missing", "extra", "foreign_context", "selection", "hash", "old23"):
            f["value"] = copy.deepcopy(original)
            if fault == "missing":
                f["value"]["artifacts"].pop()
            elif fault == "extra":
                f["value"]["artifacts"].append(copy.deepcopy(f["value"]["artifacts"][0]))
            elif fault == "foreign_context":
                f["value"]["target_run_ids"]["benchmark"] = "unused_benchmark"
            elif fault == "selection":
                f["value"]["selected_metadata"][CONTEXT].pop()
                f["value"]["artifacts"].pop()
            elif fault == "hash":
                f["value"]["artifacts"][0]["md5"] = "0" * 32
            else:
                f["value"]["schema"] = "ecoda_stage5_prepared_adoption_v1"
            plan(f, False)
        f["value"] = copy.deepcopy(original)
        # A correctly bound alternative H5AD is still the wrong scientific input.
        current = f["contexts"][CONTEXT]
        binding_path = current["new_binding"]
        binding_bytes = binding_path.read_bytes()
        binding = json.loads(binding_bytes)
        inputs = Path(binding["input_artifacts"])
        input_bytes = inputs.read_bytes()
        lines = inputs.read_text().splitlines()
        row = lines[0].split("\t")
        alternate = root / "alternate" / Path(row[2]).name
        alternate.parent.mkdir(); alternate.write_bytes(Path(row[2]).read_bytes()); m.proof(alternate)
        row[2], row[5] = str(alternate), str(alternate) + ".md5"
        lines[0] = "\t".join(row)
        inputs.chmod(0o600); inputs.write_text("\n".join(lines) + "\n")
        binding["identities"][str(inputs)] = m.digest(inputs)
        binding_path.chmod(0o600); binding_path.write_text(json.dumps(binding) + "\n")
        m.proof(binding_path)
        result = plan(f, False)
        assert "selected H5AD identity differs" in result.stderr, result.stderr
        inputs.write_bytes(input_bytes); binding_path.write_bytes(binding_bytes); m.proof(binding_path)
        # Source TASK_ID=0/NOOP is not an original exported array task.
        status = Path(f["value"]["artifacts"][0]["status"]["path"])
        status_bytes = status.read_bytes()
        status.write_text(status.read_text().replace("TASK_ID=1", "TASK_ID=0").replace("STATUS=EXPORTED", "STATUS=NOOP_VALIDATED"))
        f["value"]["artifacts"][0]["status"]["sha256"] = m.digest(status)
        plan(f, False)
        status.write_bytes(status_bytes)
        f["value"] = copy.deepcopy(original); m.save_receipt(f)
        # No R/cache preparation path is permitted for metadata-only adoption.
        forbidden = root / "forbidden-r"
        forbidden.write_text("#!/bin/bash\necho FORBIDDEN_R_OR_REFIT >&2\nexit 99\n"); forbidden.chmod(0o700)
        current["env"]["PIXI_RSCRIPT"] = str(forbidden)
        result = m.adopt(f, CONTEXT, validate_coordinator=True)
        assert result.returncode == 0, (result.stdout, result.stderr)
        new = current["new_run"]
        report = json.loads((new / "manifests/prepared_adoption.json").read_text())
        assert report["schema"] == SCHEMA and report["operation"] == "COPY_WITHOUT_REFIT"
        assert len(report["artifacts"]) == 8 and {a["kind"] for a in report["artifacts"]} == {"metadata"}
        for a in report["artifacts"]:
            assert Path(a["source"]).read_bytes() == Path(a["target"]).read_bytes()
            assert a["producer"] == "stage5_metadata_adopted_copy"
            status = new / "status/metadata_export" / f"{a['dataset']}__{a['view']}.status"
            assert "STATUS=NOOP_VALIDATED\n" in status.read_text() and "TASK_ID=0\n" in status.read_text()
            m.common(current["env"], 'ecoda_validate_artifact_record "$1" stage5_metadata_adopted_copy "$2" >/dev/null\n'
                'owner="$(ecoda_artifact_owner_dir "$1")"\n_ecoda_artifact_owner_validate_dir "$owner" "$1"\n'
                '[[ "$ECODA_ARTIFACT_OWNER_STATE" == OK ]] || exit 1\n', a["target"], new.name)
        for context in ("benchmark", "batch_effect/uncorrected_final"):
            assert not (f["contexts"][context]["new_run"] / "manifests/prepared_adoption.json").exists()
            assert not (f["contexts"][context]["new_run"] / "artifacts").exists()
        for path, expected in protected.items():
            assert (path.read_bytes(), path.stat().st_mode) == expected, path
        assert m.adopt(f, CONTEXT).returncode != 0

        # Fresh target collisions and foreign ACTIVE ownership are never reclaimed.
        for fault in ("target", "owner"):
            f = selected_fixture(root / fault)
            a = f["artifacts"][0]; current = f["contexts"][CONTEXT]
            target = current["new_run"] / "artifacts" / CONTEXT / "metadata" / Path(a["source"]).name
            if fault == "target":
                target.parent.mkdir(parents=True); target.write_bytes(b"preserve existing target")
            else:
                m.common(current["env"], 'ecoda_artifact_owner_acquire "$1" stage5 foreign_run 0 0 0 >/dev/null\n', target)
            protected = m.protected_sources(f)
            result = m.adopt(f, CONTEXT)
            assert result.returncode != 0, fault
            for path, expected in protected.items():
                assert (path.read_bytes(), path.stat().st_mode) == expected, (fault, path)
            if fault == "target": assert target.read_bytes() == b"preserve existing target"
            else: assert not target.exists()
    print("Selected original batch metadata adoption: PASS; exact8 union, no other contexts/caches/refits/jobs, strict proofs/input/status/owners")


if __name__ == "__main__":
    main()
