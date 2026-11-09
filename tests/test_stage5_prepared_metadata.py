"""Tiny, real Feather fixtures for finite prepared-artifact adoption."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile

import h5py
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind
from src.utils.py.export_h5ad_sample_metadata import requested_columns, write_metadata

BATCH_DATASETS = ("Joanito", "Stephenson", "Alzheimer", "Breast_cancer",
                  "Covid19_PBMC", "Kidney_KPMP_full", "Lupus_PBMC", "Lung")
VARIANTS = ("schvg2000", "hvg2000", "hvg500", "hvg2000_bl", "hvg1000", "hvg3000")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def proof(path):
    path = Path(path)
    Path(f"{path}.md5").write_text(
        f"MD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\nPATH={path}\n")


def fields(path, values):
    Path(path).write_text("".join(f"{key}={value}\n" for key, value in values.items()))


def common(environment, script, *arguments):
    result = subprocess.run(["bash", "-euo", "pipefail", "-c",
        'source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"\n' + script,
        "fixture", *map(str, arguments)], env=environment, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    return result


def build_fixture(base, *, run_ids=None, old_source_commit="1" * 40):
    """Return three bound old/new contexts and exactly 23 declared artifacts.

    Cache/timing bytes are opaque placeholders for the peer's RDS fixture to
    replace; metadata are actual Feather tables and inputs contain obs only.
    """
    base = Path(base).resolve()
    scratch = base / "scratch"
    nas = base / "nas"
    nas.mkdir(parents=True)
    source = base / "source"
    source.mkdir()
    config = source / "datasets.json"
    config.write_bytes((ROOT / "datasets.json").read_bytes())
    registry = json.loads(config.read_text())
    runtime = scratch / "runtime"
    runtime.mkdir(parents=True)
    image = runtime / "fixture.sif"
    image.write_bytes(b"unchanged immutable runtime fixture")
    image_manifest = runtime / "fixture.sif.manifest"
    fields(image_manifest, {"FORMAT": 2, "IMAGE_PATH": image, "IMAGE_SHA256": digest(image)})
    runtime_identity = {"RUNTIME_IMAGE": image, "RUNTIME_MANIFEST": image_manifest,
        "RUNTIME_IMAGE_SHA256": digest(image), "RUNTIME_MANIFEST_SHA256": digest(image_manifest),
        "RUNTIME_IMAGE_SIZE": image.stat().st_size,
        "RUNTIME_MANIFEST_SIZE": image_manifest.stat().st_size}
    contexts = {}
    artifacts = []
    sources = {}
    for context in ("benchmark", "batch_effect/uncorrected_final", "batch_effect/corrected_final"):
        lane = "benchmark" if context == "benchmark" else context.split("/")[1].removesuffix("_final")
        view = "benchmark_analysis" if lane == "benchmark" else "batch_effect_" + lane
        datasets = ("Stephenson",) if lane == "benchmark" else BATCH_DATASETS
        old_id, new_id = (run_ids or {}).get(context, ("old_" + lane, "new_" + lane))
        old_run = scratch / "_ecoda_runs" / old_id
        new_run = scratch / "_ecoda_runs" / new_id
        inputs = base / (lane + ".inputs.tsv")
        input_rows = []
        for dataset in datasets:
            path = scratch / "inputs" / dataset / registry[dataset]["views"][view]["output_file_name"]
            path.parent.mkdir(parents=True, exist_ok=True)
            entry = registry[dataset]
            _, required, _ = requested_columns(config, dataset, entry,
                entry["views"][view], batch_view=lane != "benchmark")
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                obs.attrs["encoding-type"] = "dataframe"
                obs.create_dataset("_index", data=["c1", "c2", "c3", "c4"], dtype=h5py.string_dtype())
                for column in required:
                    values = ([dataset + "_s1"] * 2 + [dataset + "_s2"] * 2
                              if column == "Sample" else ["present"] * 4)
                    obs.create_dataset(column, data=values, dtype=h5py.string_dtype())
            proof(path)
            input_rows.append(f"{dataset}\t{view}\t{path}\t{hashlib.md5(path.read_bytes()).hexdigest()}\t{path.stat().st_size}\t{path}.md5\n")
        inputs.write_text("".join(input_rows))
        bindings = []
        for run in (old_run, new_run):
            manifests = run / "manifests"
            manifests.mkdir(parents=True)
            fields(run / "metadata", {"STAGE": "stage5", "RUN_ID": run.name,
                                      "STATE": "FAIL" if run == old_run else "ACTIVE"})
            fields(manifests / "source.manifest", {"FORMAT": 1, "SOURCE_ROOT": source,
                "SOURCE_COMMIT": old_source_commit if run == old_run else "2" * 40,
                "DATASETS_SHA256": digest(config)})
            fields(manifests / "runtime.identity", runtime_identity)
            selection = manifests / "selection.tsv"
            selection.write_text("".join(f"{dataset}\t{view}\tall\n" for dataset in datasets))
            bindings.append(bind(run_root=run, release_id="fixture-release", nas_root=nas,
                stage="stage5", config_path=config, selection=selection, input_artifacts=inputs))
        environment = dict(os.environ, PROJECT_ROOT=str(ROOT), DATASETS_JSON_FILE=str(config),
            PYTHON_BIN=sys.executable, HPC_SCRATCH_DIR=str(scratch), NAS_TARGET_DIR=str(nas),
            ECODA_NAS_TARGET_DIR=str(nas), ECODA_RUN_ROOT=str(new_run), ECODA_RUN_ID=new_run.name,
            RUN_ID=new_run.name, PIXI_RSCRIPT=str(ROOT / ".pixi/envs/default/bin/Rscript"),
            ECODA_ARTIFACT_BINDING=str(bindings[1]), ANALYSIS_ROOT=str(new_run / "artifacts" / context),
            ANALYSIS_NAS_ROOT=str(nas / "releases/fixture-release" / new_run.name / context))
        if lane != "benchmark":
            environment.update(ANALYSIS_PASS=lane, PASS_ARG=lane,
                ANALYSIS_VARIANT="final" if lane == "uncorrected" else "corrected_final")
        contexts[context] = {"old_run": old_run, "new_run": new_run,
            "old_binding": bindings[0], "new_binding": bindings[1], "env": environment,
            "environment": environment, "inputs": inputs, "view": view}
        sources[context] = {"run_root": str(old_run)}
        for key, path in (("binding", bindings[0]), ("config", config),
            ("source_manifest", old_run / "manifests/source.manifest"),
            ("runtime_identity", old_run / "manifests/runtime.identity")):
            sources[context][key] = {"path": str(path), "sha256": digest(path)}
        names = [("pseudobulk", f"Stephenson_pseudobulk_{variant}.rds", variant)
                 for variant in VARIANTS] if lane == "benchmark" else [
                     ("metadata", f"{dataset}_sample_metadata.feather", dataset) for dataset in datasets]
        if lane == "benchmark":
            names.append(("timing", "execution_times_prepare_pseudobulk_Stephenson.feather", "Stephenson"))
        for index, (kind, filename, item) in enumerate(names, 1):
            dataset = item if kind == "metadata" else "Stephenson"
            path = (old_run / "logs" / filename if kind == "timing" else
                    old_run / "artifacts" / context / ("pseudobulks" if kind == "pseudobulk" else "metadata") / filename)
            path.parent.mkdir(parents=True, exist_ok=True)
            status = None
            record_sha = owner_sha = None
            if kind == "metadata":
                entry = registry[dataset]
                _, required, _ = requested_columns(config, dataset, entry, entry["views"][view], batch_view=True)
                frame = pd.DataFrame({column: [dataset + "_s1", dataset + "_s2"]
                    if column == "Sample" else ["present", "present"] for column in required})
                write_metadata(frame, path)
                status_path = old_run / "status/metadata_export" / f"{dataset}__{view}.status"
                status_path.parent.mkdir(parents=True, exist_ok=True)
                input_path = input_rows[index - 1].split("\t")[2]
                fields(status_path, {"STATE": "OK", "STATUS": "EXPORTED", "RUN_ID": old_run.name,
                    "ANALYSIS_VARIANT": environment["ANALYSIS_VARIANT"], "ANALYSIS_PASS": lane,
                    "ANALYSIS_ROOT": old_run / "artifacts" / context, "DATASET": dataset,
                    "VIEW": view, "TASK_ID": index, "INPUT_FILE": input_path, "OUTPUT_FILE": path})
                status = {"path": str(status_path), "sha256": digest(status_path)}
            else:
                path.write_bytes(("opaque fixture " + filename).encode())
                proof(path)
                producer = "stage5_execution_log" if kind == "timing" else "stage5_prepare_pseudobulk_" + item
                old_environment = dict(environment, ECODA_RUN_ROOT=str(old_run), ECODA_RUN_ID=old_run.name)
                result = common(old_environment,
                    'ecoda_write_artifact_record "$1" "$2" "$3" >/dev/null\n'
                    'ecoda_artifact_record_path "$1" "$3"\n', path, producer, old_run.name)
                record_sha = digest(Path(result.stdout.strip()))
                if kind == "pseudobulk":
                    result = common(old_environment,
                        'ecoda_artifact_owner_acquire "$1" stage5 "$2" 0 0 0 >/dev/null\n'
                        'ecoda_artifact_owner_set_state "$1" OK fixture >/dev/null\n'
                        'ecoda_artifact_owner_dir "$1"\n', path, old_run.name)
                    owner_sha = digest(Path(result.stdout.strip()) / "owner")
            artifacts.append({"context": context, "kind": kind, "dataset": dataset,
                "view": view, "source": str(path), "md5": hashlib.md5(path.read_bytes()).hexdigest(),
                "size": path.stat().st_size, "sidecar_sha256": digest(Path(str(path) + ".md5")),
                "record_sha256": record_sha, "owner_sha256": owner_sha, "status": status})
    evidence = []
    for name in ("cache-audit", "metadata-audit", "terminal-inspect"):
        path = base / (name + ".json")
        path.write_text(json.dumps({"scope": name, "original_gate_state": "FAILED"}) + "\n")
        evidence.append({"path": str(path), "sha256": digest(path)})
    value = {"schema": "ecoda_stage5_prepared_adoption_v1", "state": "APPROVED_COPY_ONLY",
        "original_gate_state": "FAILED", "config_sha256": digest(config),
        "target_run_ids": {key: row["new_run"].name for key, row in contexts.items()},
        "sources": sources, "evidence": evidence, "artifacts": artifacts}
    receipt = base / "prepared.json"
    receipt.write_text(json.dumps(value, sort_keys=True) + "\n")
    return {"base": base, "scratch": scratch, "nas": nas, "config": config,
            "receipt": receipt, "contexts": contexts, "artifacts": artifacts, "value": value}


def save_receipt(fixture):
    fixture["receipt"].write_text(json.dumps(fixture["value"], sort_keys=True) + "\n")


def adopt(fixture, context, *, validate_coordinator=False):
    script = '''source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"
source "$PROJECT_ROOT/src/utils/bash/ecoda_derivative_artifacts.sh"
source "$PROJECT_ROOT/src/5_run_benchmark_methods/stage5_prepared_adoption.sh"
ecoda_derivative_load stage5
stage5_adopt_prepared_artifacts "$1" "$2" "$3"
'''
    if validate_coordinator:
        script += '''
eval "$(awk '/^stage5_export_final_sample_metadata\\(\\)/ {copy=1} copy {print} copy && /^}$/ {exit}' "$PROJECT_ROOT/src/5_run_benchmark_methods/stage5_dispatcher.sh")"
source "$PROJECT_ROOT/src/utils/bash/ecoda_stage5_policy.sh"
ECODA_RECOVERY_ACTIVE=0; SYNC_ONLY_RUN=''; BENCHMARK_MATRIX_TEST=0
MANIFEST="$ECODA_RUN_ROOT/manifests/selection.tsv"; ANALYSIS_LOG_PREFIX=fixture
stage5_source_script() { printf '%s/%s\\n' "$PROJECT_ROOT" "$1"; }
stage5_require_source_script() { [[ -f "$1" ]]; }
stage5_input_path() { ecoda_derivative_input "$1" "$2"; }
sbatch() { echo 'UNEXPECTED_EXPORT_JOB' >&2; return 99; }
stage5_export_final_sample_metadata 0
grep -q '^PENDING=0$' "$ECODA_RUN_ROOT/status/metadata_export.report"
'''
    return subprocess.run(["bash", "-euo", "pipefail", "-c", script, "fixture",
        str(fixture["receipt"]), digest(fixture["receipt"]), context],
        env=fixture["contexts"][context]["env"], text=True, capture_output=True)


def protected_sources(fixture):
    paths = []
    for row in fixture["contexts"].values():
        paths.extend(path for path in row["old_run"].rglob("*") if path.is_file())
    paths.extend(path for path in (fixture["scratch"] / "inputs").rglob("*") if path.is_file())
    paths.extend(path for path in (fixture["scratch"] / "_ecoda_owners").rglob("*") if path.is_file())
    return {path: (path.read_bytes(), path.stat().st_mode) for path in paths}


def run():
    with tempfile.TemporaryDirectory(prefix="ecoda-prepared-metadata-") as temporary:
        base = Path(temporary).resolve()
        fixture = build_fixture(base / "valid")
        assert len(fixture["artifacts"]) == 23
        protected = protected_sources(fixture)
        for context in ("batch_effect/uncorrected_final", "batch_effect/corrected_final"):
            result = adopt(fixture, context, validate_coordinator=True)
            assert result.returncode == 0, result.stderr
            current = fixture["contexts"][context]
            report = json.loads((current["new_run"] / "manifests/prepared_adoption.json").read_text())
            assert report["operation"] == "COPY_WITHOUT_REFIT" and report["original_gate_state"] == "FAILED"
            assert len(report["artifacts"]) == 8
            for row in report["artifacts"]:
                original, target = Path(row["source"]), Path(row["target"])
                assert target.read_bytes() == original.read_bytes()
                assert target.stat().st_ino != original.stat().st_ino
                assert f"PATH={target}\n" in Path(str(target) + ".md5").read_text()
                assert row["producer"] == "stage5_metadata_adopted_copy"
                common(current["env"],
                    'ecoda_validate_artifact_record "$1" stage5_metadata_adopted_copy "$2" >/dev/null\n'
                    'owner="$(ecoda_artifact_owner_dir "$1")"\n'
                    '_ecoda_artifact_owner_validate_dir "$owner" "$1"\n'
                    '[[ "$ECODA_ARTIFACT_OWNER_STATE" == OK ]]\n', target, current["new_run"].name)
            result = adopt(fixture, context)
            assert result.returncode != 0, "existing adoption must be preserved, not repeated"
        assert protected_sources(fixture).items() >= protected.items()
        for path, expected in protected.items():
            assert (path.read_bytes(), path.stat().st_mode) == expected, path

        context = "batch_effect/uncorrected_final"
        for fault in ("source", "proof", "status", "evidence", "input", "target", "schema"):
            fixture = build_fixture(base / fault)
            row = next(item for item in fixture["artifacts"] if item["context"] == context)
            target = fixture["contexts"][context]["new_run"] / "artifacts" / context / "metadata" / Path(row["source"]).name
            if fault == "source":
                path = Path(row["source"])
                payload = path.read_bytes()
                path.write_bytes(bytes([payload[0] ^ 1]) + payload[1:])
            elif fault == "proof":
                Path(row["source"] + ".md5").write_text("unreviewed proof\n")
            elif fault == "status":
                path = Path(row["status"]["path"])
                path.write_text(path.read_text().replace("STATE=OK", "STATE=FAIL"))
                row["status"]["sha256"] = digest(path)
                save_receipt(fixture)
            elif fault == "evidence":
                Path(fixture["value"]["evidence"][0]["path"]).write_text("changed review\n")
            elif fault == "input":
                # Rebind the new child to a different source with valid content
                # and proof. It must not inherit metadata for the old identity.
                current = fixture["contexts"][context]
                binding = json.loads(Path(current["new_binding"]).read_text())
                input_manifest = Path(binding["input_artifacts"])
                lines = input_manifest.read_text().splitlines()
                values = lines[0].split("\t")
                alternate = base / fault / "other" / Path(values[2]).name
                alternate.parent.mkdir(parents=True)
                alternate.write_bytes(Path(values[2]).read_bytes())
                proof(alternate)
                values[2], values[5] = str(alternate), str(alternate) + ".md5"
                lines[0] = "\t".join(values)
                input_manifest.chmod(0o600)
                input_manifest.write_text("\n".join(lines) + "\n")
                binding["identities"][str(input_manifest)] = digest(input_manifest)
                Path(current["new_binding"]).chmod(0o600)
                Path(current["new_binding"]).write_text(json.dumps(binding) + "\n")
            elif fault == "target":
                target.parent.mkdir(parents=True)
                target.write_bytes(b"preexisting target must survive")
            else:
                frame = pd.read_feather(row["source"]).drop(columns=["Sample"])
                # Bypass the scientific writer only to inject a checksum-valid
                # malformed source for the adoption validator to reject.
                frame.to_feather(row["source"])
                proof(Path(row["source"]))
                row["md5"] = hashlib.md5(Path(row["source"]).read_bytes()).hexdigest()
                row["size"] = Path(row["source"]).stat().st_size
                row["sidecar_sha256"] = digest(Path(row["source"] + ".md5"))
                save_receipt(fixture)
            protected = protected_sources(fixture)
            result = adopt(fixture, context)
            assert result.returncode != 0, fault
            assert "UNEXPECTED_EXPORT_JOB" not in result.stderr
            for path, expected in protected.items():
                assert (path.read_bytes(), path.stat().st_mode) == expected, (fault, path)
            if fault == "target":
                assert target.read_bytes() == b"preexisting target must survive"
            elif fault != "schema":
                assert not target.exists(), (fault, result.stderr)
            else:
                common(fixture["contexts"][context]["env"],
                    'owner="$(ecoda_artifact_owner_dir "$1")"\n'
                    '_ecoda_artifact_owner_validate_dir "$owner" "$1"\n'
                    '[[ "$ECODA_ARTIFACT_OWNER_STATE" != OK ]]\n', target)
        print("Prepared metadata adoption: PASS")


if __name__ == "__main__":
    run()
