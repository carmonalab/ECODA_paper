"""Reviewed six-cohort UNC/Harmony geometry; selected arrays only, never expression.

The original Breast pilot CLI remains supported. Local rendering reads only bounded
display extracts; it cannot fit or score. Alzheimer memory mode reads headers only.
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
BULK_OUTPUTS = ("embedding.h5ad", "metrics.json", "display.feather", "phases.json")
DATASETS = ("Alzheimer", "Breast_cancer", "Covid19_PBMC", "Joanito", "Lupus_PBMC", "Stephenson")
COHORT_SIZE = {"Alzheimer": (1305022, 101), "Breast_cancer": (114579, 33), "Covid19_PBMC": (254174, 44),
               "Joanito": (252519, 124), "Lupus_PBMC": (458785, 94), "Stephenson": (543080, 110)}
REPRESENTATIONS = ("unc_raw", "cor_harmony")


def row_settings(dataset, representation):
    need(dataset in DATASETS and representation in REPRESENTATIONS, "Unsupported main representation")
    config = json.loads(regular(Path(os.environ["ECODA_SOURCE_ROOT"]) / "datasets.json").read_text())
    row = config[dataset]; corrected = representation == "cor_harmony"
    view = "batch_effect_corrected" if corrected else VIEW
    need(row.get("use_for_batch_effect") is True and view in row["views"] and "columns" not in row["views"][view], "Dataset flag/view/column authority differs")
    columns = row["columns"]; batch = columns["batch"]
    batch = [batch] if isinstance(batch, str) else batch
    if dataset == DATASET: need(batch == ["assay", "suspension_dissociation_time"], "Breast correction authority differs")
    else: need(len(batch) == 1, "Unreviewed multi-key display batch")
    suffix = view + "_hvg2000"; neighbor_suffix = suffix + ("_harmony" if corrected else "")
    return dict(dataset=dataset, representation=representation, view=view, sample=columns["sample"], biology=columns["label"],
                author="louvain" if dataset == "Lupus_PBMC" else columns["cell_type_low_res"],
                batch="suspension_dissociation_time" if dataset == DATASET else batch[0], correction_batch_keys=batch,
                pc=("X_pca_harmony_" if corrected else "X_pca_") + suffix,
                neighbors="neighbors_" + neighbor_suffix, leiden="leiden_res_2_" + neighbor_suffix,
                cells=COHORT_SIZE[dataset][0], samples=COHORT_SIZE[dataset][1])


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


def source_guard(c, dataset=DATASET, representation="unc_raw"):
    bulk = c["schema"] == "ecoda_singlecell_batch_main_v2"
    need(bulk or c["schema"] == "ecoda_singlecell_breast_unc_pilot_v1", "Wrong geometry control")
    settings = row_settings(dataset, representation)
    need(c["stage3_receipt"]["sha256"] == RECEIPT_SHA and c["stage3_review"]["sha256"] == REVIEW_SHA, "Unreviewed Stage3 source")
    receipt = json.loads(small(c["stage3_receipt"]))
    review = json.loads(small(c["stage3_review"]))
    need(review["status"] == "APPROVED_SELECTED_ARTIFACT_REUSE" and review["receipt_sha256"] == RECEIPT_SHA
         and review["original_gate_status"] == "FAILED" and review["failed_gate_released"] is False, "Selected reuse review differs")
    entry = next(a for a in receipt["artifacts"] if a["dataset"] == dataset and a["view"] == settings["view"])
    w = entry["scratch"]; p = regular(w["path"])
    need(p.stat().st_size == w["size"], "Parent size differs")
    proof = fields(small(w["proof"])); record = fields(small(w["record"])); owner = fields(small(w["owner"]))
    need(proof == dict(PATH=str(p), MD5=w["md5"], SIZE=str(w["size"])), "Parent native proof differs")
    need(record == w["record_fields"] and owner == w["owner_fields"] and record["STATE"] == "PUBLISHED"
         and owner["STATE"] == "OK" and owner["STAGE"] == "stage3" and record["RUN_ID"] == owner["RUN_ID"] == "stage3_confounds_20261108_v3", "Parent owner/publication differs")
    # st_dev is kernel-local; cross-node identity is native bytes plus all other
    # witness fields and the actual BeeGFS mount. In-process stats include st_dev.
    expected = c["parent_stats"][dataset][settings["view"]] if bulk else c["parent_stats"][dataset]
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


def guard_control(c):
    if c["schema"] != "ecoda_singlecell_batch_main_v2": return {(DATASET, "unc_raw"): source_guard(c)}
    need(set(c["parent_stats"]) == set(DATASETS) and all(set(c["parent_stats"][d]) == {VIEW, "batch_effect_corrected"} for d in DATASETS), "Exact twelve parent scope required")
    return {(d, r): source_guard(c, d, r) for d in DATASETS for r in REPRESENTATIONS}


def load_selected(path, expected_cells=None, settings=None):
    import anndata as ad
    from anndata.io import read_elem
    import h5py
    import numpy as np
    import pandas as pd
    if expected_cells is None: expected_cells = PILOT_CELLS
    bulk = settings is not None
    if settings is None: settings = dict(pc=PC, neighbors=NEIGHBORS, leiden=LEIDEN, sample="sample_id", author="broad_cell_type", biology="disease", batch="suspension_dissociation_time", dataset=DATASET, samples=33)
    columns = list(dict.fromkeys(["Sample", settings["sample"], settings["author"], settings["biology"], settings["batch"], settings["leiden"]] + (["assay"] if settings["dataset"] == DATASET else [])))
    pc_key = settings["pc"]; neighbor_key = settings["neighbors"]
    with h5py.File(path, "r") as h:
        obs = h["obs"]; key = obs.attrs["_index"]
        if isinstance(key, bytes): key = key.decode()
        ids = pd.Index(read_elem(obs[key]))
        need(len(ids) == expected_cells and ids.is_unique and not ids.hasnans and all(isinstance(x, str) and x for x in ids), "Invalid ordered cell IDs")
        frame = pd.DataFrame({k: read_elem(obs[k]) for k in columns}, index=ids)
        pcs = read_elem(h["obsm"][pc_key]); graph = read_elem(h["obsp"][neighbor_key + "_connectivities"])
        neighbors = read_elem(h["uns"][neighbor_key])
    need(pcs.shape == (expected_cells, 50) and np.isfinite(pcs).all(), "PC50 invalid")
    need(graph.shape == (expected_cells, expected_cells) and graph.format == "csr" and np.isfinite(graph.data).all()
         and (graph.data >= 0).all(), "Cached connectivity invalid")
    params = neighbors["params"]
    need(params["use_rep"] == pc_key and params["n_pcs"] == 50 and params["n_neighbors"] == 15 and params["metric"] == "euclidean", "Cached graph provenance differs")
    need(frame["Sample"].notna().all() and frame[settings["sample"]].notna().all()
         and frame.groupby(settings["sample"], observed=True)["Sample"].nunique().eq(1).all()
         and frame.groupby("Sample", observed=True)[settings["sample"]].nunique().eq(1).all(), "Canonical Sample mapping is not bijective")
    need(frame.Sample.nunique() == settings["samples"], "Sample universe differs")
    if settings["dataset"] == DATASET: need(set(frame["assay"].dropna()) == {"10x 3' v3"}, "Breast assay is not the reviewed constant")
    # H5AD obs retains original names; display aliases stay in the typed Feather,
    # avoiding collision with Breast's distinct source HR author_cell_type.
    value = ad.AnnData(X=None, obs=frame, obsm={pc_key: pcs}, obsp={neighbors["connectivities_key"]: graph}, uns={neighbor_key: neighbors})
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


def display_indices(ids, batch, limit=100000, full_threshold=200000):
    """Equal batch quotas with rare groups retained, using cell-ID hash ranks.

    This is display only. Neither UMAP nor LISI sees this selection. Missing
    technical labels form a display stratum, never a new scoring category.
    """
    import numpy as np
    import pandas as pd
    need(len(ids) == len(batch) and len(set(ids)) == len(ids) and 0 < limit <= full_threshold, "Display identity/budget differs")
    if len(ids) <= full_threshold: return np.arange(len(ids)), {"policy": "ALL_CELLS_LE_200000", "rows": len(ids)}
    groups = {}
    for i, label in enumerate(batch):
        key = "MISSING" if pd.isna(label) else json.dumps([type(label).__name__, str(label)], ensure_ascii=False)
        groups.setdefault(key, []).append(i)
    keys = sorted(groups); capacity = np.array([len(groups[k]) for k in keys]); quota = np.zeros(len(keys), dtype=int); left = limit
    while left:
        active = np.flatnonzero(quota < capacity); base, remainder = divmod(left, len(active))
        if not base:
            quota[active[:remainder]] += 1; break
        take = np.minimum(capacity[active] - quota[active], base); quota[active] += take; left -= int(take.sum())
    chosen = []
    for key, n in zip(keys, quota):
        ranked = sorted(groups[key], key=lambda i: (hashlib.sha256(("display-seed0\0" + ids[i]).encode()).digest(), ids[i]))
        chosen.extend(ranked[:n])
    chosen.sort(key=lambda i: ids[i])
    need(len(chosen) == len(set(chosen)) == limit, "Display budget is not exact/no-replacement")
    return np.asarray(chosen), dict(policy="BATCH_BALANCED_100000_NO_REPLACEMENT", seed=0, rows=len(chosen),
        strata=[dict(source_label=k, population=int(n), selected=int(q)) for k, n, q in zip(keys, capacity, quota)])


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


def produce(c, output, run_id, dataset=None, representation=None):
    phases = {}
    with phase(phases, "package_imports"):
        import numpy as np
        import pandas as pd
        import scanpy as sc
    bulk = dataset is not None or representation is not None
    need(output.is_absolute() and output.resolve() == output and not output.exists() and not output.is_symlink(), "Output is canonical and create-only")
    if bulk:
        need(c["schema"] == "ecoda_singlecell_batch_main_v2" and (dataset, representation) != (DATASET, "unc_raw"), "Published Breast UNC is excluded from bulk compute")
        settings = row_settings(dataset, representation); before = guard_control(c); entry = before[(dataset, representation)][0]
    else:
        need(c["schema"] == "ecoda_singlecell_breast_unc_pilot_v1", "Bulk control requires explicit row selectors")
        settings = None; entry, before = source_guard(c)
    pc_key = settings["pc"] if bulk else PC; neighbor_key = settings["neighbors"] if bulk else NEIGHBORS
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "1")); need(1 <= threads <= 8, "CPU allocation outside1..8")
    with phase(phases, "selected_input_load"): value = load_selected(entry["scratch"]["path"], settings["cells"] if bulk else None, settings)
    with phase(phases, "UMAP"):
        sc.tl.umap(value, neighbors_key=neighbor_key, random_state=0)
        need(value.obsm["X_umap"].shape == (len(value), 2) and np.isfinite(value.obsm["X_umap"]).all(), "UMAP invalid")
    with phase(phases, "ANN90"):
        nn = ann90(value.obsm[pc_key], threads)
    with phase(phases, "LISI"):
        summaries = {}
        for key, column, direction in [("batch", settings["batch"] if bulk else "suspension_dissociation_time", "mixing"), ("cell_type", settings["author"] if bulk else "broad_cell_type", "preservation")]:
            raw, scaled, summaries[key] = score_population(value.obsm[pc_key], value.obs[column], value.obs.Sample, direction, threads, nn)
            value.obs[key + "_lisi_raw"] = raw; value.obs[key + "_lisi_scaled"] = scaled
    del nn; gc.collect(); output.mkdir()
    metrics = dict(schema="ecoda_singlecell_main_metrics_v2" if bulk else "ecoda_singlecell_pilot_metrics_v1", dataset=dataset if bulk else DATASET,
                   view=settings["view"] if bulk else VIEW, representation=pc_key, run_id=run_id,
                   cells=len(value), dimensions=50, n_neighbors=90, perplexity=30, seed=0, threads=threads,
                   batch=summaries["batch"], cell_type=summaries["cell_type"], source_md5=entry["scratch"]["md5"],
                   original_global_gate_state="FAILED", stored_scores_approved=False, source_expression_read=False,
                   missing_label_policy="NA_OR_EMPTY_ONLY_NO_RECODING", UMAP_fixed_seed_optimization_may_be_serial=True,
                   ordered_cell_id_sha256=ids_digest(value.obs_names), display_policy="ALL_BREAST_CELLS_LE_200000")
    if bulk:
        metrics.update(representation_id=representation, source_columns=settings, HPC_rendering_performed=False,
            annotation_semantics="SOURCE_LOUVAIN_CLUSTER_PRESERVATION" if dataset == "Lupus_PBMC" else "SOURCE_BROAD_ANNOTATION_PRESERVATION")
        with phase(phases, "display_selection"):
            selected, policy = display_indices(list(value.obs_names), value.obs[settings["batch"]])
        metrics["display_policy"] = policy
        display = pd.DataFrame({"Sample": value.obs.Sample, "author_cell_type": value.obs[settings["author"]],
            "biological_label": value.obs[settings["biology"]], "batch_label": value.obs[settings["batch"]], "leiden2": value.obs[settings["leiden"]]}).iloc[selected].copy()
    else:
        selected = np.arange(len(value)); display = value.obs.iloc[:, :7].copy()
    display.insert(0, "cell_id", display.index)
    display[["PC1", "PC2"]] = value.obsm[pc_key][selected, :2]; display[["UMAP1", "UMAP2"]] = value.obsm["X_umap"][selected]
    if not bulk:
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
    need((guard_control(c) if bulk else source_guard(c)[1]) == before, "Parent/proofs changed during geometry computation")
    with atomic_file(output / "phases.json") as path: path.write_text(json.dumps(phases, indent=2) + "\n")
    with atomic_file(output / "files.tsv") as path:
        with path.open("w") as f:
            for name in (BULK_OUTPUTS if bulk else OUTPUTS): f.write(str(output / name) + "\t" + digest(output / name, "md5") + "\n")


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


def validate_output(c, output, run_id, dataset=None, representation=None):
    import h5py
    import numpy as np
    import pandas as pd
    from anndata.io import read_elem
    bulk = dataset is not None or representation is not None
    if bulk:
        need(c["schema"] == "ecoda_singlecell_batch_main_v2" and (dataset, representation) != (DATASET, "unc_raw"), "Bulk row scope differs")
        settings = row_settings(dataset, representation); before = guard_control(c); entry = before[(dataset, representation)][0]
    else:
        settings = None; entry, before = source_guard(c)
    n = settings["cells"] if bulk else PILOT_CELLS; pc_key = settings["pc"] if bulk else PC
    names = BULK_OUTPUTS if bulk else OUTPUTS
    rows = [line.split("\t") for line in regular(output / "files.tsv").read_text().splitlines()]
    need(len(rows) == len(names) and [r[0] for r in rows] == [str(output / name) for name in names], "Unexpected output inventory")
    for path, md5 in rows: need(digest(path, "md5") == md5, "Changed pilot payload")
    metrics = json.loads((output / "metrics.json").read_text()); phases = json.loads((output / "phases.json").read_text())
    need(metrics["schema"] == ("ecoda_singlecell_main_metrics_v2" if bulk else "ecoda_singlecell_pilot_metrics_v1") and metrics["run_id"] == run_id and metrics["cells"] == n
         and metrics["source_md5"] == entry["scratch"]["md5"] and metrics["dimensions"] == 50
         and metrics["n_neighbors"] == 90 and metrics["perplexity"] == 30 and metrics["seed"] == 0
         and metrics["representation"] == pc_key, "Report scope differs")
    if bulk: need(metrics["source_columns"] == settings and metrics["representation_id"] == representation and metrics["HPC_rendering_performed"] is False, "Bulk column/representation scope differs")
    need(set(phases) == {"package_imports", "selected_input_load", "UMAP", "ANN90", "LISI", "display_selection" if bulk else "rendering", "serialization"}
         and all(r["wall_seconds"] >= 0 and r["peak_RSS_bytes"] > 0 for r in phases.values()), "Missing phase measurements")
    with h5py.File(output / "embedding.h5ad", "r") as h:
        need(("X" not in h or h["X"].shape is None) and "raw" not in h and len(h["layers"]) == 0, "Unexpected expression in derivative")
        need(h["obsm"][pc_key].shape == (n, 50) and h["obsm"]["X_umap"].shape == (n, 2), "Derivative geometry shape differs")
        need(all(k in h["obs"] for k in ["batch_lisi_raw", "batch_lisi_scaled", "cell_type_lisi_raw", "cell_type_lisi_scaled"]), "Missing per-cell metrics")
        obs = h["obs"]; index = obs.attrs["_index"]; index = index.decode() if isinstance(index, bytes) else index
        ids = read_elem(obs[index]); need(len(ids) == len(set(ids)) == n and ids_digest(ids) == metrics["ordered_cell_id_sha256"], "Derivative ordered IDs differ")
        display = pd.read_feather(output / "display.feather"); contract = metrics["display_contract"]
        display_n = n if n <= 200000 else 100000
        if bulk: selected, display_policy = display_indices(list(ids), read_elem(obs[settings["batch"]]))
        else: selected = np.arange(n)
        need(len(display) == contract["rows"] == display_n and list(display.columns) == contract["columns"]
             and display.cell_id.is_unique and ids_digest(display.cell_id) == contract["ordered_cell_id_sha256"] == ids_digest(np.asarray(ids)[selected])
             and (output / "display.feather").stat().st_size <= contract["max_bytes"] == 32 * 1024 * 1024
             and digest(output / "display.feather", "md5") == contract["md5"], "Display scope/ID/provenance differs")
        if bulk:
            need(list(display.columns) == ["cell_id", "Sample", "author_cell_type", "biological_label", "batch_label", "leiden2", "PC1", "PC2", "UMAP1", "UMAP2"]
                 and metrics["display_policy"] == display_policy, "Unreviewed display columns/selection policy")
        for columns, key in [(["PC1", "PC2"], pc_key), (["UMAP1", "UMAP2"], "X_umap")]:
            coordinates = h["obsm"][key][:, :2][selected]
            need(np.isfinite(coordinates).all() and np.array_equal(display[columns].to_numpy(), coordinates), "Display coordinates misaligned")
        columns = {"Sample": "Sample", "author_cell_type": settings["author"], "biological_label": settings["biology"], "batch_label": settings["batch"], "leiden2": settings["leiden"]} if bulk else {k: k for k in ["Sample", "sample_id", "broad_cell_type", "disease", "assay", "suspension_dissociation_time", LEIDEN]}
        for column, original in columns.items():
            need(pd.Series(np.asarray(read_elem(obs[original]))[selected]).equals(pd.Series(display[column].to_numpy())), "Display annotation misaligned: " + column)
        if bulk:
            samples = read_elem(obs["Sample"])
            for key, column, direction in [("batch", settings["batch"], "mixing"), ("cell_type", settings["author"], "preservation")]:
                labels = np.asarray(read_elem(obs[column])); valid = pd.notna(labels) & np.asarray([str(x).strip() != "" for x in labels])
                raw = np.asarray(read_elem(obs[key + "_lisi_raw"])); scaled = np.asarray(read_elem(obs[key + "_lisi_scaled"]))
                levels = len(pd.unique(labels[valid])); expected = np.full(n, np.nan)
                if levels > 1: expected[valid] = np.clip((raw[valid] - 1) / (levels - 1) if direction == "mixing" else (levels - raw[valid]) / (levels - 1), 0, 1)
                need(np.isfinite(raw[valid]).all() and np.isnan(raw[~valid]).all() and (raw[valid] >= 1 - 1e-8).all()
                     and (raw[valid] <= levels + 1e-8).all() and np.allclose(scaled, expected, equal_nan=True, atol=1e-12), "Per-cell LISI/coverage/normalization differs")
                need(summarize(raw, scaled, samples, valid, levels) == metrics[key], "Full-population mean/coverage differs")
    need((guard_control(c) if bulk else source_guard(c)[1]) == before, "Parent changed during validation")


def load_display(row):
    """Only Root-reviewed, bounded local display+summary bytes; never a parent H5."""
    import pandas as pd
    import numpy as np
    review = json.loads(small(row["review"]))
    need(review["schema"] == "ecoda_singlecell_geometry_root_review_v1" and review["status"].startswith("APPROVED_")
         and review["reviewer"] and review["reviewed_at"] and review["source_expression_read"] is False, "Display source review missing")
    reviewed = {x["path"]: x for x in review["output_artifacts"]}
    for name, limit in [("display", 32 * 1024 * 1024), ("metrics", 2 * 1024 * 1024)]:
        ref = row[name]; p = regular(ref["path"])
        need(p.stat().st_size <= limit and reviewed[ref["source_path"]]["md5"] == ref["md5"] == digest(p, "md5"), "Unreviewed display/metrics bytes")
    metrics = json.loads(Path(row["metrics"]["path"]).read_text())
    need(row["dataset"] in DATASETS and row["representation"] in REPRESENTATIONS, "Unexpected display row")
    view = "batch_effect_corrected" if row["representation"] == "cor_harmony" else VIEW
    expected_pc = ("X_pca_harmony_" if row["representation"] == "cor_harmony" else "X_pca_") + view + "_hvg2000"
    need(metrics["dataset"] == row["dataset"] and metrics["run_id"] == review["run_id"] and metrics["dimensions"] == 50
         and metrics["view"] == view and metrics["representation"] == expected_pc and metrics["cells"] == COHORT_SIZE[row["dataset"]][0]
         and metrics["n_neighbors"] == 90 and metrics["perplexity"] == 30 and metrics["source_expression_read"] is False, "Display metric scope differs")
    frame = pd.read_feather(row["display"]["path"]); contract = metrics["display_contract"]
    need(len(frame) == contract["rows"] and frame.cell_id.is_unique and not frame.cell_id.isna().any()
         and ids_digest(frame.cell_id) == contract["ordered_cell_id_sha256"] and digest(row["display"]["path"], "md5") == contract["md5"], "Display ID contract differs")
    need(len(frame) == (metrics["cells"] if metrics["cells"] <= 200000 else 100000), "Wrong capped display population")
    if metrics["schema"] == "ecoda_singlecell_pilot_metrics_v1":
        need((row["dataset"], row["representation"]) == (DATASET, "unc_raw") and metrics["representation"] == PC
             and metrics["cells"] == PILOT_CELLS and len(frame) == PILOT_CELLS, "Only the approved Breast UNC pilot may use V1 display schema")
        mapping = dict(author_cell_type="broad_cell_type", biological_label="disease", batch_label="suspension_dissociation_time", leiden2=LEIDEN)
        frame = frame.rename(columns={v: k for k, v in mapping.items()})
        metrics["source_columns"] = dict(author="broad_cell_type", biology="disease", batch="suspension_dissociation_time", leiden=LEIDEN, sample="sample_id")
    else:
        need(metrics["schema"] == "ecoda_singlecell_main_metrics_v2" and metrics["representation_id"] == row["representation"], "Main representation differs")
    need(np.isfinite(frame[["PC1", "PC2", "UMAP1", "UMAP2"]].to_numpy()).all(), "Display coordinates not finite")
    for name in ["display", "metrics"]: need(digest(row[name]["path"], "md5") == row[name]["md5"], "Local input changed while reading")
    small(row["review"])
    return frame, metrics


def render_main(c, output):
    """Six pairs, one cohort at a time. No fitting, scoring or whole-H5 reads."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    import pandas as pd
    need(c["schema"] == "ecoda_singlecell_main_display_control_v1" and c["status"] == "APPROVED_LOCAL_MAIN_DISPLAY_EXTRACTS", "Local display control not approved")
    expected = [(d, r) for d in DATASETS for r in REPRESENTATIONS]
    need([(r["dataset"], r["representation"]) for r in c["rows"]] == expected and not output.exists() and output.resolve() == output, "Exact twelve display rows/fresh output required")
    output.mkdir(); manifest = []
    for dataset in DATASETS:
        refs = [r for r in c["rows"] if r["dataset"] == dataset]
        loaded = [load_display(r) for r in refs]; frames = [v[0] for v in loaded]; summaries = [v[1] for v in loaded]
        ids = sorted(frames[0].cell_id); need(set(ids) == set(frames[1].cell_id), "Representation display ID sets differ")
        frames = [f.set_index("cell_id").loc[ids].reset_index() for f in frames]
        for column in ["Sample", "author_cell_type", "biological_label", "batch_label"]:
            need(pd.Series(frames[0][column].to_numpy()).equals(pd.Series(frames[1][column].to_numpy())), "Cross-view source metadata differs: " + column)
        # One globally shuffled cell-ID order removes last-category draw bias and
        # is identical in every panel, row and geometry of this cohort.
        order = sorted(range(len(ids)), key=lambda i: (hashlib.sha256(("draw-seed0\0" + ids[i]).encode()).digest(), ids[i]))
        palettes = {}
        for column in ["Sample", "author_cell_type", "leiden2", "biological_label", "batch_label"]:
            levels = sorted({str(x) for f in frames for x in f[column] if pd.notna(x) and str(x).strip()})
            palettes[column] = {label: matplotlib.colors.hsv_to_rgb([int(hashlib.sha256((column + "\0" + label).encode()).hexdigest()[:8], 16) / 2**32, .65, .8]).tolist() for label in levels}
        for kind, coords, columns in [("PCA", ["PC1", "PC2"], ["Sample", "author_cell_type", "biological_label", "batch_label"]),
                                      ("UMAP", ["UMAP1", "UMAP2"], ["Sample", "author_cell_type", "leiden2", "biological_label", "batch_label"])]:
            fig, axes = plt.subplots(2, len(columns), figsize=(5 * len(columns), 16))
            for i, (frame, metrics) in enumerate(zip(frames, summaries)):
                number = lambda v: "NA (single label)" if v is None else f"{v:.4f}"
                score = "50D mean iLISI mixing=" + number(metrics["batch"]["scaled_cellwise_mean"]) + "\nmean cLISI preservation=" + number(metrics["cell_type"]["scaled_cellwise_mean"])
                for j, column in enumerate(columns):
                    ax = axes[i, j]; palette = palettes[column]
                    colors = [palette.get(str(x), [.74, .74, .74]) if pd.notna(x) and str(x).strip() else [.74, .74, .74] for x in frame[column]]
                    ax.scatter(frame[coords[0]].to_numpy()[order], frame[coords[1]].to_numpy()[order], c=np.asarray(colors)[order], s=.5, rasterized=True)
                    title = {"Sample": "Sample", "author_cell_type": "Source clusters" if dataset == "Lupus_PBMC" else "Source cell type", "leiden2": "Leiden res. 2", "biological_label": "Biology", "batch_label": "Batch: dissociation time" if dataset == DATASET else "Batch"}[column]
                    ax.set_title(title, fontsize=9, pad=27); ax.text(.5, 1.02, score, transform=ax.transAxes, fontsize=6, ha="center", va="bottom")
                    ax.set_xlabel(("PC 1" if i == 0 else "Harmony 1") if kind == "PCA" else "UMAP 1")
                    ax.set_ylabel(("UNC raw" if i == 0 else "COR Harmony") + "\n" + (coords[1] if kind == "UMAP" else ("PC 2" if i == 0 else "Harmony 2")))
                    handles = [Line2D([], [], linestyle="", marker="o", color=color, markersize=3, label=label) for label, color in palette.items()]
                    if frame[column].isna().any() or any(not str(x).strip() for x in frame[column] if pd.notna(x)):
                        handles.append(Line2D([], [], linestyle="", marker="o", color="#bdbdbd", markersize=3, label="Missing"))
                    ax.legend(handles=handles, fontsize=4 if column == "Sample" else 5, ncol=4 if column == "Sample" else 3,
                              loc="upper center", bbox_to_anchor=(.5, -.15), frameon=False)
            display_policy = "all cells" if len(ids) == summaries[0]["cells"] else "batch-balanced; not original cell abundances"
            fig.suptitle(dataset + ": actual UNC raw versus COR Harmony\n"
                         + f"Display: {len(ids):,}/{summaries[0]['cells']:,} cells per row ({display_policy}); scores: full-population 50D cellwise means", fontsize=11)
            if kind == "UMAP": fig.text(.5, .015, "Leiden IDs are row-local; equal numbers do not establish matching clusters.", ha="center", fontsize=7)
            fig.subplots_adjust(hspace=1.0, bottom=.23, top=.90, wspace=.35)
            path = output / (dataset + "_" + kind + ".pdf")
            with atomic_file(path) as tmp: fig.savefig(tmp)
            plt.close(fig)
            manifest.append(dict(path=str(path), md5=digest(path, "md5"), dataset=dataset, geometry=kind,
                display_rows=len(ids), ordered_display_id_sha256=ids_digest(ids), draw_order_sha256=ids_digest([ids[k] for k in order]),
                palettes=palettes, source_columns=[m["source_columns"] for m in summaries], inputs=refs, score_computations=0))
        del frames, loaded, summaries, frame, metrics, fig, axes, ax, ids, order, colors, handles, palette, palettes, refs
        gc.collect()
    with atomic_file(output / "render_manifest.json") as path: path.write_text(json.dumps(dict(schema="ecoda_singlecell_main_render_v1", figures=manifest, model_fits=0, score_calls=0), indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(); p.add_argument("--control", required=True); p.add_argument("--control-sha256", required=True)
    p.add_argument("--output"); p.add_argument("--run-id"); mode = p.add_mutually_exclusive_group()
    p.add_argument("--dataset"); p.add_argument("--representation", choices=REPRESENTATIONS)
    mode.add_argument("--memory-preflight", action="store_true"); mode.add_argument("--check-control", action="store_true")
    mode.add_argument("--render-main", action="store_true")
    mode.add_argument("--validate-output", action="store_true"); a = p.parse_args()
    c = json.loads(small(dict(path=a.control, sha256=a.control_sha256)))
    if a.render_main:
        need(a.output and not a.dataset and not a.representation and not a.run_id, "Local rendering has no compute-row selectors")
        render_main(c, Path(a.output)); return
    if a.memory_preflight: print(json.dumps(memory_preflight(c), indent=2)); return
    if a.check_control:
        need(not a.dataset and not a.representation, "Control check has no row selectors"); guard_control(c)
        print("SINGLECELL_MAIN12_CONTROL_OK_NO_VALUES_READ" if c["schema"] == "ecoda_singlecell_batch_main_v2" else "BREAST_UNC_PILOT_CONTROL_OK_NO_VALUES_READ"); return
    need(a.output and a.run_id, "Explicit run/output required")
    need(bool(a.dataset) == bool(a.representation), "Both row selectors are required")
    if a.validate_output: validate_output(c, Path(a.output), a.run_id, a.dataset, a.representation)
    else: produce(c, Path(a.output), a.run_id, a.dataset, a.representation)


if __name__ == "__main__": main()
