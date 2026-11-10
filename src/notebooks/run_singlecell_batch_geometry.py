"""Fixed Breast UNC pilot: selected PCs/graph only, never parent expression.

The Alzheimer mode inspects headers only. Bulk scheduling and other representations
are deliberately outside this pilot. Source artifacts and configured columns stay intact.
"""
import argparse
import contextlib
import gc
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import threading
import time
import tempfile

DATASET = "Breast_cancer"
PILOT_CELLS = 114579
VIEW = "batch_effect_uncorrected"
PC = "X_pca_batch_effect_uncorrected_hvg2000"
NEIGHBORS = "neighbors_batch_effect_uncorrected_hvg2000"
LEIDEN = "leiden_res_2_batch_effect_uncorrected_hvg2000"
RECEIPT_SHA = "c8e1d9b5b0a4b485ecdb228c706c6f8a3de65f0861fc57b148de68cd39bb05b4"
REVIEW_SHA = "5c05ac542f7aabc4924749a07159092184cec83a58c53deb59fdd68c7db2444c"
STAT = ("st_ino", "st_size", "st_mode", "st_mtime_ns", "st_ctime_ns")
OUTPUTS = ("embedding.h5ad", "metrics.json", "display.feather", "pilot_pca.pdf", "pilot_umap.pdf", "phases.json")


def need(ok, text):
    if not ok:
        raise ValueError(text)


def regular(path):
    p = Path(path)
    need(p.is_absolute() and p.resolve(strict=True) == p and p.is_file() and not p.is_symlink(), "Noncanonical file: " + str(p))
    return p


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with regular(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@contextlib.contextmanager
def atomic_file(path):
    fd, name = tempfile.mkstemp(prefix=".pilot-", suffix=path.suffix, dir=path.parent)
    os.close(fd)
    try:
        yield Path(name)
        # Atomic create-only installation; never overwrite an existing derivative.
        os.link(name, path)
    finally:
        os.unlink(name)


def ids_digest(ids):
    return hashlib.sha256(json.dumps(list(ids), ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def small(ref):
    p = regular(ref["path"])
    need(p.stat().st_size <= 2 * 1024 * 1024, "Oversized native metadata")
    raw = p.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == ref["sha256"], "Changed metadata: " + str(p))
    return raw


def fields(raw):
    rows = [line.split("=", 1) for line in raw.decode().splitlines()]
    need(all(len(r) == 2 for r in rows) and len(dict(rows)) == len(rows), "Duplicate/malformed native fields")
    return dict(rows)


def source_guard(c, dataset=DATASET):
    need(c["schema"] == "ecoda_singlecell_breast_unc_pilot_v1", "Wrong pilot control")
    need(c["stage3_receipt"]["sha256"] == RECEIPT_SHA and c["stage3_review"]["sha256"] == REVIEW_SHA, "Unreviewed Stage3 source")
    receipt = json.loads(small(c["stage3_receipt"]))
    review = json.loads(small(c["stage3_review"]))
    config = json.loads(regular(Path(os.environ["ECODA_SOURCE_ROOT"]) / "datasets.json").read_text())
    row = config[DATASET]
    need(row.get("use_for_batch_effect") is True and VIEW in row["views"] and "columns" not in row["views"][VIEW]
         and row["columns"]["batch"] == ["assay", "suspension_dissociation_time"]
         and row["columns"]["sample"] == "sample_id" and row["columns"]["cell_type_low_res"] == "broad_cell_type"
         and row["columns"]["label"] == "disease", "Top-level Breast column authority differs")
    need(review["status"] == "APPROVED_SELECTED_ARTIFACT_REUSE" and review["receipt_sha256"] == RECEIPT_SHA
         and review["original_gate_status"] == "FAILED" and review["failed_gate_released"] is False, "Selected reuse review differs")
    entry = next(a for a in receipt["artifacts"] if a["dataset"] == dataset and a["view"] == VIEW)
    w = entry["scratch"]; p = regular(w["path"])
    need(p.stat().st_size == w["size"], "Parent size differs")
    proof = fields(small(w["proof"])); record = fields(small(w["record"])); owner = fields(small(w["owner"]))
    need(proof == dict(PATH=str(p), MD5=w["md5"], SIZE=str(w["size"])), "Parent native proof differs")
    need(record == w["record_fields"] and owner == w["owner_fields"] and record["STATE"] == "PUBLISHED"
         and owner["STATE"] == "OK" and owner["STAGE"] == "stage3" and record["RUN_ID"] == owner["RUN_ID"] == "stage3_confounds_20261108_v3", "Parent owner/publication differs")
    # st_dev is kernel-local; cross-node identity is native bytes plus all other
    # witness fields and the actual BeeGFS mount. In-process stats include st_dev.
    expected = c["parent_stats"][dataset]
    actual = p.stat()
    need(all(getattr(actual, k) == expected[k] for k in STAT), "Parent stat differs")
    if sys.platform.startswith("linux"):
        mounts = []
        for line in Path("/proc/self/mountinfo").read_text().splitlines():
            cols = line.split(); i = cols.index("-")
            mount = Path(cols[4].replace("\\040", " "))
            if p == mount or mount in p.parents:
                mounts.append((len(str(mount)), cols[i + 1]))
        need(max(mounts)[1] == "beegfs", "Selected parent is not on BeeGFS")
    snapshot = {str(q): tuple(getattr(regular(q).stat(), k) for k in ("st_dev", *STAT))
                for q in [p, w["proof"]["path"], w["record"]["path"], w["owner"]["path"]]}
    return entry, snapshot


def load_selected(path, expected_cells=None):
    import anndata as ad
    from anndata.io import read_elem
    import h5py
    import numpy as np
    import pandas as pd
    if expected_cells is None: expected_cells = PILOT_CELLS
    columns = ["Sample", "sample_id", "broad_cell_type", "disease", "assay", "suspension_dissociation_time", LEIDEN]
    with h5py.File(path, "r") as h:
        obs = h["obs"]; key = obs.attrs["_index"]
        if isinstance(key, bytes): key = key.decode()
        ids = pd.Index(read_elem(obs[key]))
        need(len(ids) == expected_cells and ids.is_unique and not ids.hasnans and all(isinstance(x, str) and x for x in ids), "Invalid ordered cell IDs")
        frame = pd.DataFrame({k: read_elem(obs[k]) for k in columns}, index=ids)
        pcs = read_elem(h["obsm"][PC]); graph = read_elem(h["obsp"][NEIGHBORS + "_connectivities"])
        neighbors = read_elem(h["uns"][NEIGHBORS])
    need(pcs.shape == (expected_cells, 50) and np.isfinite(pcs).all(), "PC50 invalid")
    need(graph.shape == (expected_cells, expected_cells) and graph.format == "csr" and np.isfinite(graph.data).all()
         and (graph.data >= 0).all(), "Cached connectivity invalid")
    params = neighbors["params"]
    need(params["use_rep"] == PC and params["n_pcs"] == 50 and params["n_neighbors"] == 15 and params["metric"] == "euclidean", "Cached graph provenance differs")
    need(frame["Sample"].notna().all() and frame["sample_id"].notna().all()
         and frame.groupby("sample_id", observed=True)["Sample"].nunique().eq(1).all()
         and frame.groupby("Sample", observed=True)["sample_id"].nunique().eq(1).all(), "Canonical Sample mapping is not bijective")
    need(frame.Sample.nunique() == 33, "Breast sample universe differs")
    need(set(frame["assay"].dropna()) == {"10x 3' v3"}, "Breast assay is not the reviewed constant")
    value = ad.AnnData(X=None, obs=frame, obsm={PC: pcs}, obsp={neighbors["connectivities_key"]: graph}, uns={NEIGHBORS: neighbors})
    need(value.X is None and value.raw is None and not value.layers, "Expression-free construction failed")
    return value


def ann90(x, threads):
    import numpy as np
    from pynndescent import NNDescent
    need(len(x) > 90, "LISI population must exceed90 for perplexity30")
    index = NNDescent(x, n_neighbors=91, metric="euclidean", random_state=0, n_jobs=threads)
    ids, distances = index.neighbor_graph
    out_ids = np.empty((len(x), 90), dtype=np.int32); out_d = np.empty((len(x), 90), dtype=np.float64)
    for i in range(len(x)):
        valid = (ids[i] >= 0) & (ids[i] < len(x)) & (ids[i] != i)
        j = np.argsort(distances[i][valid], kind="stable")[:90]
        selected = ids[i][valid][j]; d = distances[i][valid][j]
        need(len(selected) == len(set(selected)) == 90 and np.isfinite(d).all() and (d >= 0).all(), "ANN90 incomplete/nonunique")
        out_ids[i] = selected; out_d[i] = d
    del index; gc.collect()
    return out_ids, out_d


def lisi_neighbors(indices, distances, codes, nlabels):
    """Entropy-calibrated inverse Simpson on90 nonself Euclidean neighbors.

    The Gaussian bandwidth search targets log(perplexity30), matching native
    LISI. Cache/JIT output is disabled so immutable source trees stay untouched.
    """
    import numpy as np
    from numba import njit
    @njit(cache=False)
    def kernel(ix, ds, labels, categories):
        out = np.empty(len(ix))
        for i in range(len(ix)):
            d = ds[i] - np.min(ds[i]); beta = 1.0; lo = 0.0; hi = np.inf
            for _ in range(50):
                w = np.exp(-beta * d); total = np.sum(w)
                entropy = np.log(total) + beta * np.sum(d * w) / total
                delta = entropy - np.log(30.0)
                if abs(delta) < 1e-5: break
                if delta > 0:
                    lo = beta; beta = beta * 2 if np.isinf(hi) else (beta + hi) / 2
                else:
                    hi = beta; beta = (beta + lo) / 2
            bucket = np.zeros(categories)
            for j in range(ix.shape[1]): bucket[labels[ix[i, j]]] += w[j] / total
            out[i] = 1.0 / np.sum(bucket * bucket)
        return out
    return kernel(indices, distances, codes, nlabels)


def summarize(raw, scaled, samples, valid, levels):
    import numpy as np
    import pandas as pd
    out = dict(evaluable_cells=int(valid.sum()), total_cells=len(valid), coverage=float(valid.mean()), actual_label_count=levels,
               raw_cellwise_mean=float(np.mean(raw[valid])), scaled_cellwise_mean=None if levels < 2 else float(np.mean(scaled[valid])))
    table = pd.DataFrame(dict(sample=np.asarray(samples)[valid], raw=raw[valid], scaled=scaled[valid]))
    means = table.groupby("sample", observed=True)[["raw", "scaled"]].mean()
    out.update(equal_sample_raw_mean=float(means.raw.mean()), equal_sample_scaled_mean=None if levels < 2 else float(means.scaled.mean()), evaluable_samples=len(means))
    return out


def score_population(pcs, labels, samples, direction, threads, shared=None):
    import numpy as np
    import pandas as pd
    valid = np.asarray(pd.notna(labels)) & np.asarray([str(x).strip() != "" for x in labels])
    need(valid.sum() > 90, "Insufficient valid labels for ANN90/perplexity30")
    codes, uniques = pd.factorize(np.asarray(labels)[valid], sort=True); levels = len(uniques)
    nn = shared if valid.all() and shared is not None else ann90(pcs[valid], threads)
    raw = np.full(len(pcs), np.nan); scaled = raw.copy()
    raw[valid] = lisi_neighbors(*nn, codes, levels)
    need(np.isfinite(raw[valid]).all() and (raw[valid] >= 1 - 1e-8).all() and (raw[valid] <= levels + 1e-8).all(), "Invalid LISI bounds")
    if levels > 1:
        scaled[valid] = np.clip((raw[valid] - 1) / (levels - 1) if direction == "mixing" else (levels - raw[valid]) / (levels - 1), 0, 1)
    return raw, scaled, summarize(raw, scaled, samples, valid, levels)


@contextlib.contextmanager
def phase(phases, name):
    def rss():
        if Path("/proc/self/statm").exists(): return int(Path("/proc/self/statm").read_text().split()[1]) * os.sysconf("SC_PAGE_SIZE")
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    start = time.perf_counter(); peak = [rss()]; stop = threading.Event()
    def monitor():
        while not stop.wait(.1): peak[0] = max(peak[0], rss())
    thread = threading.Thread(target=monitor, daemon=True); thread.start()
    try: yield
    finally:
        stop.set(); thread.join(); peak[0] = max(peak[0], rss())
        phases[name] = dict(wall_seconds=time.perf_counter() - start, peak_RSS_bytes=peak[0], rss_sampling_seconds=.1)


def render(display, metrics, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    columns = ["Sample", "broad_cell_type", "disease", "suspension_dissociation_time"]
    annotation = "50D normalized cellwise means: iLISI mixing={}; cLISI preservation={}".format(
        metrics["batch"]["scaled_cellwise_mean"], metrics["cell_type"]["scaled_cellwise_mean"])
    for kind, axes, labels in [("pca", ["PC1", "PC2"], columns), ("umap", ["UMAP1", "UMAP2"], columns[:2] + [LEIDEN] + columns[2:])]:
        fig, panels = plt.subplots(1, len(labels), figsize=(5 * len(labels), 5), squeeze=False)
        for ax, column in zip(panels[0], labels):
            values = pd.Categorical(display[column]); colors = plt.get_cmap("turbo")(np.linspace(0, 1, max(1, len(values.categories))))
            for code, label in enumerate(values.categories):
                keep = values.codes == code; ax.scatter(display.loc[keep, axes[0]], display.loc[keep, axes[1]], s=.5, color=colors[code], label=str(label), rasterized=True)
            missing = values.codes < 0
            if missing.any(): ax.scatter(display.loc[missing, axes[0]], display.loc[missing, axes[1]], s=.5, color="#bdbdbd", label="Missing", rasterized=True)
            ax.set(xlabel=axes[0], ylabel=axes[1], title=column); ax.legend(fontsize=4, markerscale=4, loc="best", ncol=2)
        fig.suptitle("Breast UNC pilot; " + annotation, fontsize=8); fig.tight_layout()
        with atomic_file(output / ("pilot_" + kind + ".pdf")) as path: fig.savefig(path)
        plt.close(fig)


def produce(c, output, run_id):
    phases = {}
    with phase(phases, "package_imports"):
        import numpy as np
        import pandas as pd
        import scanpy as sc
    need(not output.exists() and not output.is_symlink(), "Pilot output is create-only")
    entry, before = source_guard(c); threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "1")); need(1 <= threads <= 8, "Pilot CPU allocation outside1..8")
    with phase(phases, "selected_input_load"): value = load_selected(entry["scratch"]["path"])
    with phase(phases, "UMAP"):
        sc.tl.umap(value, neighbors_key=NEIGHBORS, random_state=0)
        need(value.obsm["X_umap"].shape == (len(value), 2) and np.isfinite(value.obsm["X_umap"]).all(), "UMAP invalid")
    with phase(phases, "ANN90"):
        nn = ann90(value.obsm[PC], threads)
    with phase(phases, "LISI"):
        summaries = {}
        for key, column, direction in [("batch", "suspension_dissociation_time", "mixing"), ("cell_type", "broad_cell_type", "preservation")]:
            raw, scaled, summaries[key] = score_population(value.obsm[PC], value.obs[column], value.obs.Sample, direction, threads, nn)
            value.obs[key + "_lisi_raw"] = raw; value.obs[key + "_lisi_scaled"] = scaled
    del nn; gc.collect(); output.mkdir()
    metrics = dict(schema="ecoda_singlecell_pilot_metrics_v1", dataset=DATASET, view=VIEW, representation=PC, run_id=run_id,
                   cells=len(value), dimensions=50, n_neighbors=90, perplexity=30, seed=0, threads=threads,
                   batch=summaries["batch"], cell_type=summaries["cell_type"], source_md5=entry["scratch"]["md5"],
                   original_global_gate_state="FAILED", stored_scores_approved=False, source_expression_read=False,
                   missing_label_policy="NA_OR_EMPTY_ONLY_NO_RECODING", UMAP_fixed_seed_optimization_may_be_serial=True,
                   ordered_cell_id_sha256=ids_digest(value.obs_names), display_policy="ALL_BREAST_CELLS_LE_200000")
    # All Breast cells are displayed by the approved <=200000 rule; no sampling RNG is used.
    display = value.obs.iloc[:, :7].copy(); display.insert(0, "cell_id", display.index)
    display[["PC1", "PC2"]] = value.obsm[PC][:, :2]; display[["UMAP1", "UMAP2"]] = value.obsm["X_umap"]
    with phase(phases, "rendering"): render(display, summaries, output)
    with phase(phases, "serialization"):
        # Scalar JSON retains explicit null/degenerate-score fields without
        # relying on HDF5's heterogeneous nested-dictionary serialization.
        value.uns["singlecell_pilot_json"] = json.dumps(metrics, allow_nan=False)
        with atomic_file(output / "embedding.h5ad") as path: value.write_h5ad(path)
        with atomic_file(output / "display.feather") as path: display.reset_index(drop=True).to_feather(path)
        need((output / "display.feather").stat().st_size <= 32 * 1024 * 1024, "Display extract exceeds32MiB; do not silently thin")
        metrics["display_contract"] = dict(rows=len(display), columns=list(display.columns), ordered_cell_id_sha256=ids_digest(display.cell_id),
                                            md5=digest(output / "display.feather", "md5"), max_bytes=32 * 1024 * 1024)
        with atomic_file(output / "metrics.json") as path: path.write_text(json.dumps(metrics, allow_nan=False, indent=2) + "\n")
    need(source_guard(c)[1] == before, "Parent/proofs changed during pilot")
    with atomic_file(output / "phases.json") as path: path.write_text(json.dumps(phases, indent=2) + "\n")
    with atomic_file(output / "files.tsv") as path:
        with path.open("w") as f:
            for name in OUTPUTS: f.write(str(output / name) + "\t" + digest(output / name, "md5") + "\n")


def memory_preflight(c):
    import h5py
    entry, before = source_guard(c, "Alzheimer")
    with h5py.File(entry["scratch"]["path"], "r") as h:
        pc = h["obsm"][PC]; graph = h["obsp"][NEIGHBORS + "_connectivities"]
        need(pc.shape == (1305022, 50), "Alzheimer PC header differs")
        arrays = {k: dict(shape=list(graph[k].shape), dtype=str(graph[k].dtype), uncompressed_bytes=int(graph[k].size * graph[k].dtype.itemsize),
                         allocated_HDF5_bytes=int(graph[k].id.get_storage_size())) for k in ("data", "indices", "indptr")}
        pc_dtype = str(pc.dtype); pc_storage = int(pc.id.get_storage_size())
        payload = int(pc.size * pc.dtype.itemsize) + sum(x["uncompressed_bytes"] for x in arrays.values())
        n = pc.shape[0]; ann = n * 90 * (8 + 4)
    need(source_guard(c, "Alzheimer")[1] == before, "Alzheimer metadata changed")
    return dict(schema="ecoda_alzheimer_geometry_header_memory_v1", cells=n, pc_shape=[n, 50], pc_dtype=pc_dtype,
                PC_allocated_HDF5_bytes=pc_storage, parent_H5AD_bytes=entry["scratch"]["size"], graph_arrays=arrays,
                PC_and_graph_payload_bytes=payload, ANN90_index64_distance32_bytes=ann,
                conservative_provisional_allocation_bytes=max(128 * 1024**3, 8 * (payload + ann)),
                overhead_factor_not_measured=8, full_Alzheimer_compute_approved=False, coordinate_values_read=0, expression_objects_opened=0)


def validate_output(c, output, run_id):
    import h5py
    import numpy as np
    import pandas as pd
    from anndata.io import read_elem
    entry, before = source_guard(c)
    rows = [line.split("\t") for line in regular(output / "files.tsv").read_text().splitlines()]
    need(len(rows) == len(OUTPUTS) and [r[0] for r in rows] == [str(output / n) for n in OUTPUTS], "Unexpected pilot output inventory")
    for path, md5 in rows: need(digest(path, "md5") == md5, "Changed pilot payload")
    metrics = json.loads((output / "metrics.json").read_text()); phases = json.loads((output / "phases.json").read_text())
    need(metrics["schema"] == "ecoda_singlecell_pilot_metrics_v1" and metrics["run_id"] == run_id and metrics["cells"] == PILOT_CELLS
         and metrics["source_md5"] == entry["scratch"]["md5"] and metrics["dimensions"] == 50
         and metrics["n_neighbors"] == 90 and metrics["perplexity"] == 30 and metrics["seed"] == 0, "Pilot report scope differs")
    need(set(phases) == {"package_imports", "selected_input_load", "UMAP", "ANN90", "LISI", "rendering", "serialization"}
         and all(r["wall_seconds"] >= 0 and r["peak_RSS_bytes"] > 0 for r in phases.values()), "Missing phase measurements")
    with h5py.File(output / "embedding.h5ad", "r") as h:
        need(("X" not in h or h["X"].shape is None) and "raw" not in h and len(h["layers"]) == 0, "Unexpected expression in derivative")
        need(h["obsm"][PC].shape == (PILOT_CELLS, 50) and h["obsm"]["X_umap"].shape == (PILOT_CELLS, 2), "Derivative geometry shape differs")
        need(all(k in h["obs"] for k in ["batch_lisi_raw", "batch_lisi_scaled", "cell_type_lisi_raw", "cell_type_lisi_scaled"]), "Missing per-cell metrics")
        obs = h["obs"]; index = obs.attrs["_index"]; index = index.decode() if isinstance(index, bytes) else index
        ids = read_elem(obs[index]); need(ids_digest(ids) == metrics["ordered_cell_id_sha256"], "Derivative ordered IDs differ")
        display = pd.read_feather(output / "display.feather"); contract = metrics["display_contract"]
        need(len(display) == contract["rows"] == PILOT_CELLS and list(display.columns) == contract["columns"]
             and display.cell_id.is_unique and ids_digest(display.cell_id) == contract["ordered_cell_id_sha256"] == ids_digest(ids)
             and (output / "display.feather").stat().st_size <= contract["max_bytes"] == 32 * 1024 * 1024
             and digest(output / "display.feather", "md5") == contract["md5"], "Display scope/ID/provenance differs")
        for columns, key in [(["PC1", "PC2"], PC), (["UMAP1", "UMAP2"], "X_umap")]:
            coordinates = h["obsm"][key][:, :2]
            need(np.isfinite(coordinates).all() and np.array_equal(display[columns].to_numpy(), coordinates), "Display coordinates misaligned")
        for column in ["Sample", "sample_id", "broad_cell_type", "disease", "assay", "suspension_dissociation_time", LEIDEN]:
            need(pd.Series(np.asarray(read_elem(obs[column]))).equals(pd.Series(display[column].to_numpy())), "Display annotation misaligned: " + column)
    need(source_guard(c)[1] == before, "Parent changed during validation")


def main():
    p = argparse.ArgumentParser(); p.add_argument("--control", required=True); p.add_argument("--control-sha256", required=True)
    p.add_argument("--output"); p.add_argument("--run-id"); mode = p.add_mutually_exclusive_group()
    mode.add_argument("--memory-preflight", action="store_true"); mode.add_argument("--check-control", action="store_true")
    mode.add_argument("--validate-output", action="store_true"); a = p.parse_args()
    c = json.loads(small(dict(path=a.control, sha256=a.control_sha256)))
    if a.memory_preflight: print(json.dumps(memory_preflight(c), indent=2)); return
    if a.check_control: source_guard(c); print("BREAST_UNC_PILOT_CONTROL_OK_NO_VALUES_READ"); return
    need(a.output and a.run_id, "Explicit run/output required")
    if a.validate_output: validate_output(c, Path(a.output), a.run_id)
    else: produce(c, Path(a.output), a.run_id)


if __name__ == "__main__": main()
