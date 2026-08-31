"""Small fail-closed artifact validator used by scheduler watchdogs.

The validator intentionally has no project-specific scientific assumptions: for
h5ad it checks the minimum AnnData groups needed to prove that an atomic write
finished; RDS/other binary artifacts are validated by non-empty size and their
checksum sidecar at the caller.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def validate_artifact(path: str | Path, kind: str) -> None:
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size == 0:
        raise ValueError(f"artifact is missing or empty: {artifact}")
    if kind != "h5ad":
        return
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - worker environment issue
        raise RuntimeError("h5py is required to validate h5ad artifacts") from exc
    with h5py.File(artifact, "r") as handle:
        missing = [group for group in ("X", "obs", "var") if group not in handle]
        if missing:
            raise ValueError(f"h5ad missing required groups {missing}: {artifact}")
        obs = handle["obs"]
        if "_index" in obs.attrs:
            index_name = obs.attrs["_index"]
            if isinstance(index_name, bytes):
                try:
                    index_name = index_name.decode("utf-8")
                except UnicodeDecodeError:
                    index_name = None
            if (
                not isinstance(index_name, str)
                or not index_name.strip()
                or index_name not in obs
            ):
                raise ValueError(f"h5ad obs index is missing: {artifact}")
        elif "_index" not in obs:
            raise ValueError(f"h5ad obs index is missing: {artifact}")
        shape = handle.attrs.get("shape")
        if shape is not None and any(int(value) <= 0 for value in shape):
            raise ValueError(f"h5ad shape is empty: {artifact}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", required=True)
    parser.add_argument("--kind", choices=("h5ad", "binary"), default="binary")
    args = parser.parse_args()
    validate_artifact(args.path, args.kind)
    print(f"artifact contract OK: {args.path}")


if __name__ == "__main__":
    main()
