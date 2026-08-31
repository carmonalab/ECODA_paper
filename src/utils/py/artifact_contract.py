"""Small fail-closed artifact validator used by scheduler watchdogs.

The validator intentionally has no project-specific scientific assumptions: for
h5ad it checks the minimum AnnData groups needed to prove that an atomic write
finished; RDS/other binary artifacts are validated by non-empty size and their
checksum sidecar at the caller.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def _shape_is_positive(shape: object) -> bool:
    try:
        dimensions = tuple(shape)
    except TypeError:
        return False
    if len(dimensions) != 2:
        return False
    for dimension in dimensions:
        try:
            integer = int(dimension)
        except (TypeError, ValueError, OverflowError):
            return False
        if integer <= 0 or integer != dimension:
            return False
    return True


def _group_has_storage(group: object, h5py: object) -> bool:
    for child_name in group:
        child = group[child_name]
        if isinstance(child, h5py.Dataset):
            if child.size > 0:
                return True
        elif isinstance(child, h5py.Group) and _group_has_storage(child, h5py):
            return True
    return False


def _matrix_node_has_storage(node: object, h5py: object) -> bool:
    if isinstance(node, h5py.Dataset):
        return node.ndim == 2 and node.size > 0
    if not isinstance(node, h5py.Group):
        return False
    shape = node.attrs.get("shape")
    if shape is not None and not _shape_is_positive(shape):
        return False
    return _group_has_storage(node, h5py)


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
        missing = [group for group in ("obs", "var") if group not in handle]
        if missing:
            raise ValueError(f"h5ad missing required groups {missing}: {artifact}")
        invalid = [
            group for group in ("obs", "var") if not isinstance(handle[group], h5py.Group)
        ]
        if invalid:
            raise ValueError(f"h5ad required groups are not groups {invalid}: {artifact}")

        matrix_candidates = ("X", "layers/counts")
        matrix_name = next(
            (
                candidate
                for candidate in matrix_candidates
                if candidate in handle and _matrix_node_has_storage(handle[candidate], h5py)
            ),
            None,
        )
        if matrix_name is None:
            if any(candidate in handle for candidate in matrix_candidates):
                raise ValueError(f"h5ad matrix storage is empty: {artifact}")
            raise ValueError(
                f"h5ad missing matrix storage (expected X or layers/counts): {artifact}"
            )

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
        if shape is not None and not _shape_is_positive(shape):
            raise ValueError(f"h5ad shape is empty or invalid: {artifact}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", required=True)
    parser.add_argument("--kind", choices=("h5ad", "binary"), default="binary")
    args = parser.parse_args()
    validate_artifact(args.path, args.kind)
    print(f"artifact contract OK: {args.path}")


if __name__ == "__main__":
    main()
