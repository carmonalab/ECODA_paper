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


if __name__ == "__main__":
    unittest.main()
