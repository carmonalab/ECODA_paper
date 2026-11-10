"""Freeze three compact application inputs; never read H5AD expression data.

Run locally with Pixi. The selected obs-only query runs through Bamboo's
configured interpreter via SSH stdin; it installs nothing on the cluster.
Yan's historical import/archive is documented separately in the application plan.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


def md5(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "md5").hexdigest()


def collect(args):
    import h5py
    import numpy as np

    remote_repo = Path(args.remote_repo)
    sys.path.insert(0, str(remote_repo / "src/utils/py"))
    from h5ad_source_identity import read_obs_column_values, _sidecar_fields

    entry = json.loads(args.entry_json)
    path = Path(args.remote_root) / args.dataset / "output" / entry["output_file"]
    nas_path = Path(args.remote_nas) / args.dataset / "output" / entry["output_file"]
    source_sidecar = _sidecar_fields(path, verify_content=False)
    nas_sidecar = _sidecar_fields(nas_path, verify_content=False)
    digest = source_sidecar["MD5"]
    if digest != nas_sidecar["MD5"]:
        raise ValueError("Scratch/NAS source checksum identities disagree")
    if path.stat().st_size <= 0 or path.stat().st_size != nas_path.stat().st_size:
        raise ValueError("Scratch/NAS source sizes disagree")

    clinical = {
        "Bassez": ["expansion", "BC_type", "timepoint", "cohort"],
        "Zhang": ["Clinical.efficacy.", "Tissue", "Treatment", "Treatment_type"],
        "Pelka": ["MMRStatus", "Tissue"],
    }[args.dataset]
    sample_col = entry["sample_col"]
    columns = list(dict.fromkeys([sample_col, "Sample", "layer2"] + clinical
                                + list(entry["subset_vars"])))
    before = path.stat()
    with h5py.File(path, "r") as handle:
        if "ecoda_annotation_version" not in handle.get("uns", {}):
            raise ValueError("Source lacks automatic-annotation provenance marker")
        obs = handle["obs"]
        data = {col: read_obs_column_values(obs, col) for col in columns}
    if (before.st_size, before.st_mtime_ns) != (path.stat().st_size, path.stat().st_mtime_ns):
        raise ValueError("Source changed during metadata extraction")

    def clean(value):
        if value is None or (isinstance(value, (float, np.floating)) and np.isnan(value)):
            return None
        return str(value)

    data = {key: [clean(value) for value in values] for key, values in data.items()}
    n_cells = len(data[sample_col])
    if not n_cells or any(len(values) != n_cells for values in data.values()):
        raise ValueError("Empty or inconsistent obs columns")
    if any(value is None or not value for value in data[sample_col]):
        raise ValueError("Missing configured sample identity")
    # Canonical preprocessing introduces Sample; prove the configured ID mapping
    # is bijective rather than silently treating it as a column override.
    aliases = {}
    reverse = {}
    for source_id, canonical_id in zip(data[sample_col], data["Sample"]):
        if source_id in aliases and aliases[source_id] != canonical_id:
            raise ValueError("Configured sample splits across canonical identities")
        if canonical_id in reverse and reverse[canonical_id] != source_id:
            raise ValueError("Canonical sample merges configured identities")
        aliases[source_id] = canonical_id
        reverse[canonical_id] = source_id
    for column, rule in entry["subset_vars"].items():
        values = set(map(str, rule["values"]))
        if rule["op"] not in ("in", "notin"):
            raise ValueError("Unsupported application source restriction")
        matches = [value in values for value in data[column]]
        if not all(matches if rule["op"] == "in" else [not x for x in matches]):
            raise ValueError(f"Source violates configured restriction: {column}")

    samples = list(dict.fromkeys(data[sample_col]))
    metadata = []
    for sample in samples:
        rows = [i for i, value in enumerate(data[sample_col]) if value == sample]
        record = {"sample": sample, "canonical_sample": aliases[sample]}
        for column in clinical:
            values = {data[column][i] for i in rows if data[column][i] is not None}
            if len(values) > 1:
                raise ValueError(f"Conflicting sample metadata: {sample}/{column}")
            record[column] = next(iter(values)) if values else None
        metadata.append(record)
    missing = {None, "", "NA", "nan", "None", "Unknown"}
    types = list(dict.fromkeys(value for value in data["layer2"] if value not in missing))
    if not types:
        raise ValueError("No valid HiTME layer-2 annotations")
    sample_idx = {value: i for i, value in enumerate(samples)}
    type_idx = {value: i for i, value in enumerate(types)}
    counts = [[0] * len(types) for _ in samples]
    for sample, label in zip(data[sample_col], data["layer2"]):
        if label not in missing:
            counts[sample_idx[sample]][type_idx[label]] += 1
    obs_hash = hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
    return {
        "schema_version": 1, "dataset": args.dataset,
        "samples": samples, "cell_types": types, "counts": counts,
        "metadata": metadata,
        "provenance": {
            "host": "bamboo", "source": str(path), "nas_source": str(nas_path),
            "source_md5_identity": digest, "source_size": before.st_size,
            "source_uid": before.st_uid, "source_mtime_ns": before.st_mtime_ns,
            "checksum_scope": "Recorded existing H5AD sidecars; matrix bytes not rehashed. Selected obs columns SHA256 verified by snapshot.",
            "obs_sha256": obs_hash, "obs_columns": columns, "n_cells": n_cells,
            "sample_column": sample_col, "annotation_column": "layer2",
            "decoder_md5": md5(remote_repo / "src/utils/py/h5ad_source_identity.py"),
            "python": sys.version, "h5py": h5py.__version__, "numpy": np.__version__,
            "dataset_contract": entry,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="data/applications/cancer_clustering/20261108-v1")
    parser.add_argument("--nas-root", default="/Volumes/Shared/Projects/ECODA_paper")
    parser.add_argument("--host", default="bamboo", choices=["bamboo"])
    parser.add_argument("--remote-repo", default="/home/users/h/halterc/ECODA_paper")
    parser.add_argument("--remote-root", default="/home/users/h/halterc/scratch/ECODA_paper")
    parser.add_argument("--remote-nas", default="/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper")
    parser.add_argument("--remote-python", default="/home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13/bin/python")
    parser.add_argument("--collect", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--dataset", help=argparse.SUPPRESS)
    parser.add_argument("--entry-json", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.collect:
        print(json.dumps(collect(args), allow_nan=False))
        return

    repo = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo / "src/utils/py"))
    from datasets_io import read_datasets_json
    config = read_datasets_json(repo / "datasets.json", view="benchmark_analysis")
    output = Path(args.output_root).resolve()
    benchmark = Path(args.nas_root) / "benchmark"
    checksum_text = (benchmark / "checksums.md5").read_text()
    recorded = {row.split(maxsplit=1)[1].lstrip("*"): row.split()[0]
                for row in checksum_text.splitlines() if row.strip()}
    for dataset in ("Bassez", "Zhang", "Pelka"):
        entry = config[dataset]
        if not entry["use_for_benchmark"] or entry["view_name"] != "benchmark_analysis":
            raise ValueError(f"Dataset not enabled: {dataset}")
        selected = output / dataset
        snapshot_path = selected / "obs_snapshot.json"
        bundle_name = f"{dataset}_ECODA_HiTME_HR_layer2.rds"
        source_bundle = benchmark / "results" / bundle_name
        expected = recorded[f"results/{bundle_name}"]
        selected.mkdir(parents=True, exist_ok=True)
        if md5(source_bundle) != expected:
            raise ValueError(f"Source bundle checksum mismatch: {dataset}")
        bundle = selected / bundle_name
        if bundle.exists() and md5(bundle) != expected:
            raise ValueError(f"Existing imported bundle differs: {bundle}")
        if not bundle.exists():
            temporary = selected / ("." + bundle_name + ".tmp")
            subprocess.run(["rsync", "--checksum", "--", str(source_bundle), str(temporary)], check=True)
            if md5(temporary) != expected:
                raise ValueError("Transferred bundle checksum mismatch")
            os.replace(temporary, bundle)
        bundle_sidecar = selected / (bundle_name + ".md5")
        expected_sidecar = f"{expected}  {bundle_name}\n"
        if bundle_sidecar.exists() and bundle_sidecar.read_text() != expected_sidecar:
            raise ValueError("Existing imported sidecar differs")
        if not bundle_sidecar.exists():
            bundle_sidecar.write_text(expected_sidecar)

        command = shlex.join([
            args.remote_python, "-", "--collect", "--dataset", dataset,
            "--entry-json", json.dumps(entry), "--remote-repo", args.remote_repo,
            "--remote-root", args.remote_root, "--remote-nas", args.remote_nas,
        ])
        result = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
                                 args.host, command], input=Path(__file__).read_text(),
                                text=True, capture_output=True)
        if result.returncode:
            raise RuntimeError(f"{dataset}: read-only Bamboo query failed:\n{result.stderr}")
        snapshot = json.loads(result.stdout)
        snapshot["provenance"].update(
            bundle_md5=expected, checksum_manifest_md5=md5(benchmark / "checksums.md5"),
            datasets_json_md5=md5(repo / "datasets.json"), importer_md5=md5(__file__),
        )
        payload = json.dumps(snapshot, indent=2, allow_nan=False) + "\n"
        if snapshot_path.exists() and snapshot_path.read_text() != payload:
            raise ValueError(f"Existing snapshot differs; select a new run root: {snapshot_path}")
        if not snapshot_path.exists():
            with tempfile.NamedTemporaryFile(mode="w", dir=selected, delete=False) as handle:
                handle.write(payload)
                temporary = Path(handle.name)
            os.replace(temporary, snapshot_path)
        digest_path = Path(str(snapshot_path) + ".md5")
        line = f"{md5(snapshot_path)}  {snapshot_path.name}\n"
        if digest_path.exists() and digest_path.read_text() != line:
            raise ValueError("Existing snapshot sidecar differs")
        if not digest_path.exists():
            digest_path.write_text(line)
        print(f"{dataset}: froze {len(snapshot['samples'])} samples, "
              f"{len(snapshot['cell_types'])} types; {snapshot_path}", flush=True)


if __name__ == "__main__":
    main()
