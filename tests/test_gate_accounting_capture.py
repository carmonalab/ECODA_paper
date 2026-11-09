"""Offline terminal-accounting bridge regression; no Slurm/SSH invocation."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shlex
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("gate_capture", ROOT / "src/utils/py/gate_accounting_capture.py")
capture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(capture)


class CaptureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name).resolve()
        self.run = self.base / "scratch/_ecoda_runs/stage3-fixture"
        self.folder = self.run / "manifests"
        self.folder.mkdir(parents=True)
        self.profile = {"name": "fixture", "policy": {"accounting_command": capture.ACCOUNTING_COMMAND}}
        self.accounting = {"attempted": True, "query_count": 1, "actual_exit": 0,
                           "scheduler_ids": ["123", "124"],
                           "command": "sacct -n -P -X -j 123,124 --format=JobIDRaw,State,ExitCode"}
        self.raw("123|COMPLETED|0:0\n125|COMPLETED|0:0\n124|COMPLETED|0:0\n")
        self.gate = {"gate_id": "fixture", "profile": "fixture", "state": "COMPLETED",
                     "audit_state": "COMPLETED", "audit_completed": True,
                     "scheduler_array_ids": ["123"], "scheduler_watchdog_ids": ["124"],
                     "exact_command": shlex.join(["executor", "--run-id", self.run.name,
                                                  "--scratch-root", str(self.base / "scratch")])}
        self.gate["command_digest"] = hashlib.sha256(self.gate["exact_command"].encode()).hexdigest()
        self.child = "array\t123\nwatchdog\t124\n"

    def raw(self, stdout):
        self.accounting["stdout"] = stdout
        self.accounting["rows"] = []
        for index, line in enumerate(stdout.splitlines(), 1):
            fields = line.split("|")
            if len(fields) == 3:
                self.accounting["rows"].append(dict(line=index, job_id=fields[0], state=fields[1], exit_code=fields[2]))

    def arguments(self):
        profile = json.dumps(self.profile).encode()
        self.gate["profile_digest"] = hashlib.sha256(profile).hexdigest()
        self.gate["audit"] = {"passed": True, "accounting": copy.deepcopy(self.accounting)}
        inspected = {"operation": "inspect", "gate_id": "fixture", "state": "COMPLETED",
                     "audit": copy.deepcopy(self.gate["audit"])}
        files = {"gate_manifest": ("gate.json", json.dumps(self.gate).encode()),
                 "inspect": ("inspect.json", json.dumps(inspected).encode()),
                 "profile": ("profile.json", profile),
                 "scheduler_ids": ("scheduler_ids.tsv", self.child.encode())}
        arguments = {"run_root": str(self.run)}
        for name, (filename, data) in files.items():
            path = self.folder / filename
            path.write_bytes(data)
            arguments[name] = str(path)
        return arguments

    def test_preserves_raw_and_projects_only_actual_successful_roots(self):
        arguments = self.arguments()
        proof = capture.create_capture(**arguments)
        self.assertEqual(proof["role"], "DERIVED_ROOT_PROJECTION")
        self.assertEqual(proof["raw_row_count"], 3)
        self.assertEqual((self.folder / "scheduler_accounting.raw.tsv").read_text(), self.accounting["stdout"])
        self.assertEqual((self.folder / "scheduler_accounting.tsv").read_text(), "123|COMPLETED|0:0\n124|COMPLETED|0:0\n")
        for filename in ("scheduler_accounting.raw.tsv", "scheduler_accounting.tsv", "scheduler_accounting.proof.json"):
            path = self.folder / filename
            sidecar = dict(line.split("=", 1) for line in Path(f"{path}.md5").read_text().splitlines())
            self.assertEqual(sidecar, {"PATH": str(path), "SIZE": str(path.stat().st_size),
                                       "MD5": hashlib.md5(path.read_bytes()).hexdigest()})

    def test_all_raw_rows_must_succeed_even_unmatched_children(self):
        for state, code in (("FAILED", "1:0"), ("OUT_OF_MEMORY", "0:125"), ("RUNNING", "0:0"), ("COMPLETED", "0:1")):
            self.raw(f"123|COMPLETED|0:0\n125|{state}|{code}\n124|COMPLETED|0:0\n")
            with self.subTest(state=state), self.assertRaisesRegex(ValueError, "every raw"):
                capture.create_capture(**self.arguments())

    def test_missing_root_duplicate_malformed_and_truncated_rejected(self):
        for stdout in ("123|COMPLETED|0:0\n", "123|COMPLETED|0:0\n123|COMPLETED|0:0\n124|COMPLETED|0:0\n",
                       "123_1|COMPLETED|0:0\n124|COMPLETED|0:0\n", "123|COMPLETED|0:0\n\n124|COMPLETED|0:0\n",
                       "123|123|COMPLETED|0:0\n124|124|COMPLETED|0:0\n", "123|COMPLETED|0:0\n124|COMPLETED|0:0",
                       "123|COMPLETED|0:0\n[output truncated]"):
            self.raw(stdout)
            with self.subTest(stdout=stdout), self.assertRaises(ValueError):
                capture.create_capture(**self.arguments())

    def test_scope_query_and_wrapper_binding_rejected(self):
        for key, value in (("command", "sacct -n -P -X -j 123,124,999 --format=JobIDRaw,State,ExitCode"),
                           ("query_count", 2), ("query_count", True), ("actual_exit", 1),
                           ("scheduler_ids", ["123"]), ("attempted", False)):
            previous = self.accounting[key]
            self.accounting[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                capture.create_capture(**self.arguments())
            self.accounting[key] = previous
        self.child += "extra\t999\n"
        with self.assertRaisesRegex(ValueError, "gate union"):
            capture.create_capture(**self.arguments())
        self.child = "array\t123\nwatchdog\t124\n"
        self.gate["exact_command"] += " --run-id other"
        with self.assertRaises(ValueError):
            capture.create_capture(**self.arguments())

    def test_tampered_evidence_and_profile_fail(self):
        arguments = self.arguments()
        path = Path(arguments["profile"])
        path.write_bytes(path.read_bytes() + b" ")
        with self.assertRaisesRegex(ValueError, "profile identity"):
            capture.create_capture(**arguments)
        arguments = self.arguments()
        path = Path(arguments["inspect"])
        inspected = json.loads(path.read_text())
        inspected["audit"]["accounting"]["stdout"] += "999|COMPLETED|0:0\n"
        path.write_text(json.dumps(inspected))
        with self.assertRaisesRegex(ValueError, "same completed audit"):
            capture.create_capture(**arguments)
        self.accounting["rows"][0]["job_id"] = "999"
        with self.assertRaisesRegex(ValueError, "parsed evidence"):
            capture.create_capture(**self.arguments())

    def test_command_tamper_with_unchanged_digest_fails(self):
        self.gate["exact_command"] += " --unexpected changed"
        with self.assertRaisesRegex(ValueError, "recorded digest"):
            capture.create_capture(**self.arguments())

    def test_existing_projection_or_partial_sidecar_preserved(self):
        for name in ("scheduler_accounting.tsv", "scheduler_accounting.raw.tsv.md5"):
            path = self.folder / name
            path.write_bytes(b"existing evidence")
            with self.assertRaisesRegex(ValueError, "already exists"):
                capture.create_capture(**self.arguments())
            self.assertEqual(path.read_bytes(), b"existing evidence")
            path.unlink()


class CampaignTests(unittest.TestCase):
    """Real receipt layouts from stage5_campaign.sh; only scheduler data is fake."""
    raw = CaptureTests.raw
    arguments = CaptureTests.arguments

    def setUp(self):
        CaptureTests.setUp(self)
        self.project = self.base / "snapshot/tree"
        self.project.mkdir(parents=True)
        (self.project / "datasets.json").write_text('{"Toy": {}}\n')
        self.source = (f"SOURCE_ROOT={self.project}\nAUX_ROOT={self.project}/aux\n"
                       f"DATASETS_SHA256={hashlib.sha256((self.project / 'datasets.json').read_bytes()).hexdigest()}\n").encode()
        self.runtime = b"RUNTIME_IMAGE=/pinned/image.sif\nRUNTIME_MANIFEST=/pinned/image.sif.manifest\n"
        original_source = self.base / "source.manifest"
        original_source.write_bytes(self.source)
        self.context_request = self.base / "contexts.tsv"
        self.gate["exact_command"] = shlex.join(["executor", "--run-id", self.run.name,
            "--scratch-root", str(self.base / "scratch"), "--source-root", str(self.project),
            "--source-manifest", str(original_source), "--runtime-image", "/pinned/image.sif",
            "--runtime-manifest", "/pinned/image.sif.manifest", "--script",
            "src/5_run_benchmark_methods/stage5_campaign.sh", "--contexts", str(self.context_request),
            "--derivative-release", "fixture-release"])
        self.gate["command_digest"] = hashlib.sha256(self.gate["exact_command"].encode()).hexdigest()
        self.gate["scheduler_array_ids"] = ["123", "124", "126"]
        self.gate["scheduler_watchdog_ids"] = []
        self.accounting["scheduler_ids"] = ["123", "124", "126"]
        self.accounting["command"] = capture.ACCOUNTING_COMMAND.replace("{scheduler_ids}", "123,124,126")
        self.raw("123|COMPLETED|0:0\n125|COMPLETED|0:0\n124|COMPLETED|0:0\n126|COMPLETED|0:0\n")
        self.seal(self.folder / "source.manifest", self.source, "stage5_campaign_identity")
        self.seal(self.folder / "runtime.identity", self.runtime, "stage5_campaign_identity")
        self.children = {}
        contexts, identities, original, schedulers = [], [], [], []
        for context, view, identifier in (("benchmark", "benchmark_analysis", "123"),
                                          ("uncorrected", "batch_effect_uncorrected", "124"),
                                          ("corrected", "batch_effect_corrected", "126")):
            child = f"child_{context}"
            root = self.run.parent / child
            folder = root / "manifests"
            folder.mkdir(parents=True)
            (root / "status").mkdir()
            self.children[context] = root
            selection = f"Toy\t{view}\t{'pseudobulk' if context == 'benchmark' else view}\n".encode()
            references = f"Toy\t{view}\t/input\tmd5\t1\t/proof\n".encode()
            scope = self.folder / f"campaign_{context}.selection.tsv"
            inputs = self.folder / f"campaign_{context}.inputs.tsv"
            self.seal(scope, selection, "stage5_campaign_scope")
            self.seal(inputs, references, "stage5_campaign_scope")
            for path in (scope, inputs):
                identities.append(f"{path}\t{hashlib.sha256(path.read_bytes()).hexdigest()}\n")
            original.extend([f"{context}\t{child}\t{scope}\t{inputs}\n"])
            contexts.extend(original[-1:])
            schedulers.append(f"{context}\t{child}\tARRAY\t{identifier}\n")
            (folder / "source.manifest").write_bytes(self.source)
            (folder / "runtime.identity").write_bytes(self.runtime)
            self.seal(folder / "selection.tsv", selection)
            (folder / "input_artifacts.tsv").write_bytes(references)
            binding = dict(schema="ecoda_derivative_artifacts_v1", stage="stage5", run_id=child,
                run_root=str(root), release_id="fixture-release", selection=str(folder / "selection.tsv"),
                input_artifacts=str(folder / "input_artifacts.tsv"), identities={
                    str(folder / name): hashlib.sha256((folder / name).read_bytes()).hexdigest()
                    for name in ("source.manifest", "runtime.identity", "selection.tsv", "input_artifacts.tsv")})
            binding["identities"][str(self.project / "datasets.json")] = hashlib.sha256(
                (self.project / "datasets.json").read_bytes()).hexdigest()
            if context == "benchmark":
                request = folder / "requested_selection.tsv"
                self.seal(request, selection)
                binding["requested_selection"] = str(request)
                binding["identities"][str(request)] = hashlib.sha256(selection).hexdigest()
            self.seal(folder / "artifact_bindings.json", json.dumps(binding).encode())
            self.seal(folder / "scheduler_ids.tsv", f"ARRAY\t{identifier}\n".encode())
            (root / "status/terminal").write_text(f"RUN_ID={child}\nSTATE=OK\n")
            (self.run / "status").mkdir(exist_ok=True)
            (self.run / "status" / context).write_text(
                f"CONTEXT={context}\nRUN_ID={child}\nEXIT_CODE=0\nSCHEDULER_MANIFEST={folder}/scheduler_ids.tsv\n")
        self.seal(self.folder / "campaign_contexts.tsv", "".join(contexts).encode(), "stage5_campaign_scope")
        self.seal(self.folder / "campaign_identities.tsv", "".join(identities).encode(), "stage5_campaign_scope")
        self.seal(self.folder / "campaign_requested.tsv", "".join(original).encode(), "stage5_campaign_scope")
        self.context_request.write_text("".join(original))
        self.seal(self.folder / "campaign_scheduler_ids.tsv", "".join(schedulers).encode())
        (self.run / "status/terminal").write_text(f"RUN_ID={self.run.name}\nSTATE=OK\n")

    def seal(self, path, data, producer=None):
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            path.chmod(0o600)
        path.write_bytes(data)
        sidecar = f"MD5={hashlib.md5(data).hexdigest()}\nSIZE={len(data)}\nPATH={path}\n"
        Path(f"{path}.md5").write_text(sidecar)
        if producer:
            record = path.parent / "artifacts" / f"{hashlib.sha256(str(path).encode()).hexdigest()[:32]}.record"
            record.parent.mkdir(exist_ok=True)
            record.write_text(sidecar + f"RUN_ID={self.run.name}\nPRODUCER={producer}\nSTATE=PUBLISHED\n")
            path.chmod(0o444)

    def campaign_arguments(self):
        arguments = self.arguments()
        arguments.pop("scheduler_ids")
        return arguments

    def test_campaign_success_preserves_one_raw_and_three_actual_root_projections(self):
        result = capture.create_campaign_capture(**self.campaign_arguments())
        self.assertEqual(len(result["children"]), 3)
        self.assertFalse((self.folder / "scheduler_accounting.tsv").exists())
        self.assertEqual((self.folder / "scheduler_accounting.raw.tsv").read_text(), self.accounting["stdout"])
        for context, identifier in (("benchmark", "123"), ("uncorrected", "124"), ("corrected", "126")):
            folder = self.children[context] / "manifests"
            self.assertEqual((folder / "scheduler_accounting.tsv").read_text(), f"{identifier}|COMPLETED|0:0\n")
            proof = json.loads((folder / "scheduler_accounting.proof.json").read_text())
            self.assertEqual(proof["context"], context)
            self.assertEqual(proof["raw"]["path"], str(self.folder / "scheduler_accounting.raw.tsv"))
            self.assertEqual(proof["gate_requested_ids"], ["123", "124", "126"])

    def test_campaign_failure_in_unmatched_child_blocks_all_publication(self):
        self.raw("123|COMPLETED|0:0\n125|OUT_OF_MEMORY|0:125\n124|COMPLETED|0:0\n126|COMPLETED|0:0\n")
        with self.assertRaisesRegex(ValueError, "every raw"):
            capture.create_campaign_capture(**self.campaign_arguments())
        self.assertFalse((self.folder / "scheduler_accounting.raw.tsv").exists())

    def test_campaign_union_scope_source_runtime_and_receipt_tamper(self):
        targets = [self.folder / "campaign_scheduler_ids.tsv",
                   self.children["benchmark"] / "manifests/requested_selection.tsv",
                   self.children["corrected"] / "manifests/input_artifacts.tsv",
                   self.children["uncorrected"] / "manifests/source.manifest",
                   self.children["corrected"] / "manifests/runtime.identity",
                   self.run / "status/benchmark"]
        for path in targets:
            original = path.read_bytes()
            path.write_bytes(original + b"tamper\n")
            with self.subTest(path=path), self.assertRaises(ValueError):
                capture.create_campaign_capture(**self.campaign_arguments())
            path.write_bytes(original)

    def test_campaign_union_is_checked_even_with_valid_changed_checksum(self):
        path = self.folder / "campaign_scheduler_ids.tsv"
        self.seal(path, path.read_bytes().replace(b"\t126\n", b"\t999\n"))
        with self.assertRaisesRegex(ValueError, "scheduler union"):
            capture.create_campaign_capture(**self.campaign_arguments())

    def test_campaign_requires_sealed_contexts_and_no_existing_child_projection(self):
        path = self.folder / "campaign_contexts.tsv"
        path.chmod(0o644)
        with self.assertRaisesRegex(ValueError, "sealed producer"):
            capture.create_campaign_capture(**self.campaign_arguments())
        path.chmod(0o444)
        projection = self.children["corrected"] / "manifests/scheduler_accounting.tsv"
        projection.write_text("existing immutable evidence")
        with self.assertRaisesRegex(ValueError, "already exists"):
            capture.create_campaign_capture(**self.campaign_arguments())
        self.assertEqual(projection.read_text(), "existing immutable evidence")
        self.assertFalse((self.folder / "scheduler_accounting.raw.tsv").exists())




if __name__ == "__main__":
    unittest.main()
