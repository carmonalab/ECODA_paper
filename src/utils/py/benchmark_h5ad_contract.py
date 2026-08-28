"""Content contract for benchmark-analysis AnnData artifacts."""

import os
import numpy as np


REQUIRED_OBSM = {
    "benchmark_analysis": {
        "X_pca_benchmark_analysis_hvg1000",
        "X_pca_benchmark_analysis_hvg2000",
        "X_pca_benchmark_analysis_hvg3000",
        "X_pca_harmony_benchmark_analysis_hvg2000",
    },
    "batch_effect_uncorrected": {
        "X_pca_batch_effect_uncorrected_hvg2000",
    },
    "batch_effect_corrected": {
        "X_pca_batch_effect_corrected_hvg2000",
        "X_pca_harmony_batch_effect_corrected_hvg2000",
    },
}


def _h5_shape(node):
    """Return a node's persisted shape, including group-backed values."""
    if node is None:
        return ()
    shape = getattr(node, "shape", None)
    if shape is None:
        shape = node.attrs.get("shape")
    if shape is None:
        return ()
    try:
        return tuple(int(value) for value in shape)
    except (TypeError, ValueError):
        return ()


def _h5_nonempty_node(node):
    shape = _h5_shape(node)
    return bool(shape) and all(value > 0 for value in shape)


def _h5_vector_non_na(node):
    """Return the number of non-NA values in a scalar/vector HDF5 field."""
    value_node = node
    if hasattr(node, "keys"):
        for name in ("values", "data", "_values", "codes"):
            if name in node:
                value_node = node[name]
                break
    try:
        values = np.asarray(value_node[()]).reshape(-1)
    except (TypeError, ValueError, OSError):
        return 0
    if "mask" in getattr(node, "keys", lambda: ())():
        try:
            mask = np.asarray(node["mask"][()]).reshape(-1).astype(bool)
            return int((~mask).sum())
        except (TypeError, ValueError, OSError):
            return 0
    if getattr(value_node, "name", "").endswith("/codes"):
        try:
            return int((values >= 0).sum())
        except TypeError:
            return 0
    if np.issubdtype(values.dtype, np.number):
        return int((~np.isnan(values)).sum())
    return int(values.size)


def _h5_index_name(group):
    index_name = group.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    return str(index_name)


def _validate_count_values(values, label):
    values = np.asarray(values)
    if values.size == 0:
        raise ValueError(f"{label} is empty")
    try:
        finite = np.isfinite(values).all()
        nonnegative = (values >= 0).all()
        integral = np.equal(values, np.floor(values)).all()
    except (TypeError, ValueError):
        raise ValueError(f"{label} is not numeric") from None
    if not finite:
        raise ValueError(f"{label} contains nonfinite values")
    if not nonnegative or not integral:
        raise ValueError(f"{label} must be finite, nonnegative, integer-valued counts")


def _validate_counts_layer(layer, shape, label="layers['counts']"):
    import scipy.sparse as sparse

    layer_shape = _h5_shape(layer)
    if layer_shape and tuple(shape) != layer_shape:
        raise ValueError(f"{label} shape {layer_shape} does not match {tuple(shape)}")
    if sparse.issparse(layer):
        values = layer.data
    else:
        values = np.asarray(layer)
    _validate_count_values(values, label)
def _validate_h5_count_node(node, label, max_values=1024 * 1024):
    shape = _h5_shape(node)
    if not shape:
        raise ValueError(f"{label} has no persisted shape")
    if hasattr(node, "shape"):
        if len(shape) == 1:
            for start in range(0, shape[0], max_values):
                _validate_count_values(
                    node[start:start + max_values], label
                )
        else:
            rows_per = max(1, max_values // max(1, int(np.prod(shape[1:])))
                           )
            for start in range(0, shape[0], rows_per):
                _validate_count_values(node[start:start + rows_per], label)
    else:
        raise ValueError(f"{label} is not a readable HDF5 dataset")


def _validate_h5_counts_layer(layer, expected_shape, label="layers['counts']"):
    shape = _h5_shape(layer)
    if tuple(shape) != tuple(expected_shape):
        raise ValueError(
            f"{label} shape {shape} does not match {tuple(expected_shape)}"
        )
    if hasattr(layer, "keys"):
        encoding = layer.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding != "csr_matrix":
            raise ValueError(f"{label} has unsupported encoding {encoding!r}")
        if "data" not in layer:
            raise ValueError(f"{label} has no data values")
        _validate_h5_count_node(layer["data"], label)
    else:
        _validate_h5_count_node(layer, label)


def validate_benchmark_h5ad_contract(adata, view, method):
    """Reject incomplete AnnData artifacts before benchmark computation."""
    if view not in REQUIRED_OBSM:
        raise ValueError(f"Unknown preprocessing view for h5ad contract: {view}")

    missing = []
    try:
        x_present = adata.X is not None
    except Exception:
        x_present = False
    if not x_present:
        missing.append("X")
    if "counts" not in adata.layers:
        missing.append("layers['counts']")
    if "counts" in adata.layers:
        _validate_counts_layer(
            adata.layers["counts"], (adata.n_obs, adata.n_vars)
        )

    obsm = getattr(adata, "obsm", {})
    missing_obsm = sorted(REQUIRED_OBSM[view] - set(obsm.keys()))
    missing.extend(f"obsm['{key}']" for key in missing_obsm)

    if "hvg_rank" not in adata.var.columns:
        missing.append("var['hvg_rank']")
    else:
        required_hvg = 3000 if view == "benchmark_analysis" else 2000
        n_ranked = int(adata.var["hvg_rank"].notna().sum())
        if n_ranked < required_hvg:
            missing.append(
                f"var['hvg_rank'] with at least {required_hvg} non-NA ranks "
                f"(found {n_ranked})"
            )

    if "Sample" not in adata.obs.columns:
        missing.append("obs['Sample']")
    if getattr(adata, "n_obs", 0) <= 0 or getattr(adata, "n_vars", 0) <= 0:
        missing.append("non-empty obs/var")

    if missing:
        raise ValueError(
            f"h5ad content contract failed for {method} ({view}): "
            f"missing or invalid {', '.join(missing)}. "
            "Re-run 1.1.1_preprocess.py with --force and use the "
            "authoritative processed h5ad."
        )


def validate_benchmark_h5ad_path(path, view, method):
    """Validate persisted h5ad structure without materializing AnnData."""
    import h5py

    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError(f"h5ad path is missing or empty: {path}")
    if view not in REQUIRED_OBSM:
        raise ValueError(f"Unknown preprocessing view for h5ad contract: {view}")

    missing = []
    with h5py.File(path, "r") as handle:
        x_shape = _h5_shape(handle["X"]) if "X" in handle else ()
        if "X" not in handle:
            missing.append("X")
        elif len(x_shape) != 2 or any(value <= 0 for value in x_shape):
            missing.append("X with a non-empty persisted shape")
        if "layers" not in handle or "counts" not in handle["layers"]:
            missing.append("layers['counts']")
        else:
            counts = handle["layers"]["counts"]
            counts_shape = _h5_shape(counts)
            expected_shape = x_shape if len(x_shape) == 2 else counts_shape
            try:
                _validate_h5_counts_layer(counts, expected_shape)
            except (OSError, TypeError, ValueError) as exc:
                missing.append(str(exc))
        if "obsm" not in handle:
            missing.extend(f"obsm['{key}']" for key in sorted(REQUIRED_OBSM[view]))
        else:
            for key in sorted(REQUIRED_OBSM[view]):
                if key not in handle["obsm"]:
                    missing.append(f"obsm['{key}']")

        required_hvg = 3000 if view == "benchmark_analysis" else 2000
        if "var" not in handle or "hvg_rank" not in handle["var"]:
            missing.append("var['hvg_rank']")
        else:
            hvg_rank = handle["var"]["hvg_rank"]
            if _h5_vector_non_na(hvg_rank) < required_hvg:
                missing.append(
                    f"var['hvg_rank'] with at least {required_hvg} non-NA ranks"
                )

        if "obs" not in handle or "Sample" not in handle["obs"]:
            missing.append("obs['Sample']")
        if "obs" not in handle:
            missing.append("non-empty obs index")
        else:
            obs = handle["obs"]
            index_name = _h5_index_name(obs)
            if index_name not in obs or not _h5_nonempty_node(obs[index_name]):
                missing.append("non-empty obs index")
        if "var" not in handle:
            missing.append("non-empty var index")
        else:
            var = handle["var"]
            index_name = _h5_index_name(var)
            if index_name not in var or not _h5_nonempty_node(var[index_name]):
                missing.append("non-empty var index")

    if missing:
        raise ValueError(
            f"h5ad content contract failed for {method} ({view}): "
            f"missing or invalid {', '.join(dict.fromkeys(missing))}. "
            "Re-run 1.1.1_preprocess.py with --force and use the "
            "authoritative processed h5ad."
        )



def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", required=True)
    parser.add_argument("--view", required=True)
    parser.add_argument("--method", required=True)
    args = parser.parse_args()
    validate_benchmark_h5ad_path(args.path, args.view, args.method)
    print(f"h5ad contract OK: {args.path}")


if __name__ == "__main__":
    main()
