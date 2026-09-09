from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/christianhalter/Desktop/ECODA_paper")
EMBEDDINGS = ROOT / "data/benchmark/embeddings"
CANDIDATE = EMBEDDINGS / "execution_times.candidate_before_mem_repair_1788900485.feather"
REPAIRED_BASELINE = EMBEDDINGS / "execution_times.repaired_baseline_before_restore_1788900485.feather"
CANONICAL = EMBEDDINGS / "execution_times.feather"
EVIDENCE = EMBEDDINGS / "execution_times_candidate_restore_1788900485.json"

EXPECTED_CANDIDATE_MD5 = "3d9cdb84633e02e52606feb3b0ffccf4"
EXPECTED_REPAIRED_BASELINE_MD5 = "bc1fc17936ebd6eaba61f2c7f9058ebf"
COLUMNS = ["dataset", "method", "time_secs", "mem_GB"]


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    if not path.is_file() or path.stat().st_size <= 0:
        raise RuntimeError(f"missing/empty source: {path}")
    actual = md5(path)
    if actual != expected:
        raise RuntimeError(f"source hash mismatch for {path}: {actual}")


require_hash(CANDIDATE, EXPECTED_CANDIDATE_MD5)
require_hash(REPAIRED_BASELINE, EXPECTED_REPAIRED_BASELINE_MD5)
if EVIDENCE.exists():
    raise RuntimeError(f"refusing to overwrite evidence: {EVIDENCE}")

candidate = pd.read_feather(CANDIDATE)
repaired_baseline = pd.read_feather(REPAIRED_BASELINE)
for name, frame, expected_rows in (
    ("candidate", candidate, 1403),
    ("repaired baseline", repaired_baseline, 1376),
):
    if list(frame.columns) != COLUMNS:
        raise RuntimeError(f"{name} schema mismatch: {list(frame.columns)}")
    if len(frame) != expected_rows:
        raise RuntimeError(f"{name} row count mismatch: {len(frame)}")

key_columns = ["dataset", "method"]
candidate_keys = list(zip(candidate.dataset, candidate.method))
baseline_keys = list(zip(repaired_baseline.dataset, repaired_baseline.method))
if len(set(candidate_keys)) != len(candidate_keys):
    raise RuntimeError("candidate has duplicate dataset/method keys")
if len(set(baseline_keys)) != len(baseline_keys):
    raise RuntimeError("repaired baseline has duplicate dataset/method keys")

baseline_by_key = repaired_baseline.set_index(key_columns)
keys = pd.MultiIndex.from_arrays([candidate.dataset, candidate.method])
source_mem = baseline_by_key.reindex(keys)["mem_GB"].to_numpy()
fillable = candidate["mem_GB"].isna().to_numpy() & pd.notna(source_mem)
if int(fillable.sum()) != 84:
    raise RuntimeError(f"unexpected fill count: {int(fillable.sum())}")

restored = candidate.copy()
restored.loc[fillable, "mem_GB"] = source_mem[fillable]

# The candidate is the authoritative row universe for this restoration: preserve
# all identifiers, ordering, runtimes, and already-populated memory values.
if not restored[key_columns].equals(candidate[key_columns]):
    raise RuntimeError("restore changed identifiers/order")
if not restored["time_secs"].equals(candidate["time_secs"]):
    raise RuntimeError("restore changed time_secs")
nonfilled = ~fillable
if not restored.loc[nonfilled, "mem_GB"].equals(candidate.loc[nonfilled, "mem_GB"]):
    raise RuntimeError("restore changed a non-target mem_GB value")
if restored[key_columns].duplicated().any():
    raise RuntimeError("restored duplicate dataset/method keys")
if restored["mem_GB"].notna().any():
    values = pd.to_numeric(restored.loc[restored["mem_GB"].notna(), "mem_GB"])
    if not values.map(lambda value: pd.notna(value) and value >= 0 and value != float("inf")).all():
        raise RuntimeError("restored invalid mem_GB value")

# Install the restored candidate atomically after write/read-back validation.
tmp = EMBEDDINGS / f".execution_times.restore.tmp.{os.getpid()}"
restored.to_feather(tmp)
try:
    read_back = pd.read_feather(tmp)
    if not read_back.equals(restored):
        raise RuntimeError("restored Feather failed read-back equality")
    if md5(tmp) == "":
        raise RuntimeError("restored Feather hash is empty")
    os.replace(tmp, CANONICAL)
finally:
    if tmp.exists():
        tmp.unlink()

output_md5 = md5(CANONICAL)
output_size = CANONICAL.stat().st_size
changed = candidate.loc[fillable, key_columns].copy()
changed["mem_GB"] = restored.loc[fillable, "mem_GB"].to_numpy()
changed_records = changed.to_dict(orient="records")

payload = {
    "schema_version": 1,
    "status": "RESTORED_CANDIDATE_WITH_TARGETED_MEMORY_FILL",
    "candidate_source": {
        "path": str(CANDIDATE),
        "md5": EXPECTED_CANDIDATE_MD5,
        "size": int(CANDIDATE.stat().st_size),
        "rows": int(len(candidate)),
    },
    "repaired_baseline_source": {
        "path": str(REPAIRED_BASELINE),
        "md5": EXPECTED_REPAIRED_BASELINE_MD5,
        "size": int(REPAIRED_BASELINE.stat().st_size),
        "rows": int(len(repaired_baseline)),
    },
    "canonical_output": {
        "path": str(CANONICAL),
        "md5": output_md5,
        "size": int(output_size),
        "rows": int(len(restored)),
    },
    "changed_row_count": int(fillable.sum()),
    "changed_rows": changed_records,
    "candidate_only_rows_preserved": int(len(set(candidate_keys) - set(baseline_keys))),
    "unresolved_candidate_na_count": int(restored["mem_GB"].isna().sum()),
}
EVIDENCE.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps({
    "canonical_md5": output_md5,
    "canonical_size": output_size,
    "rows": len(restored),
    "changed_rows": int(fillable.sum()),
    "candidate_only_preserved": payload["candidate_only_rows_preserved"],
    "unresolved_na": payload["unresolved_candidate_na_count"],
}, sort_keys=True))
