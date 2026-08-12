#!/usr/bin/env python3
#
# 1.1_prepare_chunks.py — Build per-dataset sample chunk files (production or test mode)
#
# Native Python replacement for the former reticulate-based 1.1_prepare_chunks.r.
#
# One chunk set per DATASET (not per view): all preprocessed view .h5ad files in
# ${HPC_SCRATCH_DIR}/${DS_NAME}/output are merged into a single union h5ad
# (deduplicated on (sample, barcode)), written OUTSIDE the output dir to
# ${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union/union.h5ad so the NAS sync
# glob over output/*.h5ad stays clean. Samples are grouped into consecutive
# chunks of 2 (or 1 in --test mode) and chunk_<N>.txt files are written under
# ${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks. Line 1 of each chunk file is the
# absolute union .h5ad path, subsequent lines are sample IDs.
#
# The union is ALWAYS written as a minimal h5ad (even for single-view datasets,
# whose preprocessed view keeps log-norm X + layers["counts"] and must not be
# chunked directly): X = raw counts (float64 CSR), obs + var only, no
# layers/obsp/obsm/uns/varm/varp. anndata's read_h5ad(backed="r") eagerly
# materializes the entire `layers` group into RAM at open (verified in anndata
# 0.12.10 and the HPC-pinned 0.12.19), so a union that kept counts in layers
# would add a full-matrix memory floor per annotation worker (Kfoury:
# +1.31 GB; a GongSharma-class union: 26-40 GB — OOM at any chunk size). With
# counts in X the backed open costs ~10 MB; X stays lazy and per-sample backed
# subsetting stays selective (CSR-on-disk). The annotation worker
# (2.1.1_process_chunk.R) treats the counts-absent -> X fallback as the
# designed primary path for union files (message, not warning).
#
# Union construction is memory-lean for the common case:
#   - exactly one view  -> streaming h5py copy (obs/var groups + counts as X)
#   - multiple views, one of which already contains the full union (e.g.
#     Stephenson: benchmark view is a cell subset of the batch-effect view)
#     -> same streaming copy of that view (h5py-only obs key-set reads — no
#     anndata open, since backed open eagerly materializes the whole layers
#     group into RAM on anndata 0.12.x (~6-10 GB for a full-cohort view, vs
#     ~100-200 MB for the key-set read); keeps the 4G srun allocation valid)
#   - otherwise (partial overlaps between views) -> full in-memory concat +
#     dedup; this path needs a large srun allocation and prints a warning.

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pyarrow as pa
import pyarrow.feather as pafeather
import scipy.sparse as sp
import scanpy as sc


# Legacy annotation column classification for the clean-entry check (--check-
# clean). Tier 1: unconditional pattern matches (any scGate_*,
# functional.cluster*, *_UCell, scATOMIC*, layer<digit>*/layer_<digit>*
# column — a plain `layer` column survives). Tier 2: exact names only, so
# author columns like a plain `cellCycle` or `Phase_variant` never collide.
# Must mirror 3.1_merge_annotations.py.
LEGACY_ANNOT_TIER1 = [r"^scGate", r"^functional\.cluster", r"_UCell$",
                      r"^scATOMIC", r"^layer_?\d"]
LEGACY_ANNOT_TIER2 = ["S.Score", "G2M.Score", "Phase", "classification_confidence",
                      "cellCycle.G1S_UCell", "cellCycle.G2M_UCell"]


def cell_keys(adata, sample_col):
    """(sample, barcode) composite keys for the obs of `adata`."""
    samples = adata.obs[sample_col].astype(str).values
    return np.char.add(np.char.add(samples, "_"), adata.obs_names.astype(str).values)


def _read_str_dataset(ds):
    """
    Read a 1-D string h5py dataset into a numpy str array. Covers the
    encodings anndata writes for obs string columns (fixed-length bytes
    'S', variable-length h5py string dtype 'O', and numeric columns, which
    mirror pandas `.astype(str)` semantics). Fails closed on anything else.
    """
    arr = ds[:]
    if arr.dtype.kind == "S":
        return np.char.decode(arr, "utf-8")
    if arr.dtype.kind == "O":
        out = np.empty(arr.shape, dtype=object)
        for i, x in enumerate(arr):
            if isinstance(x, bytes):
                out[i] = x.decode("utf-8")
            else:
                out[i] = str(x)
        return out.astype(str)
    if arr.dtype.kind in "iuf":
        return arr.astype(str)
    raise RuntimeError(
        f"Unsupported dtype '{arr.dtype}' for h5py dataset '{ds.name}'."
    )


def _read_obs_column(obs_grp, col):
    """
    Read one obs column as a numpy str array, mirroring the semantics of
    `adata.obs[col].astype(str)` for the encodings anndata 0.12 writes
    (plain string/numeric datasets and 'categorical' groups). A missing
    column raises KeyError like pandas; unknown encodings raise
    RuntimeError (fail closed, never silently mis-decode).
    """
    ds = obs_grp[col]
    enc = ds.attrs.get("encoding-type")
    if isinstance(enc, bytes):
        enc = enc.decode("utf-8")
    if enc is None or enc in ("", "array", "string-array"):
        return _read_str_dataset(ds)
    if enc == "categorical":
        categories = _read_str_dataset(ds["categories"])
        codes = ds["codes"][:].astype(np.int64)
        # -1 = missing category; pandas .astype(str) renders NaN as "nan".
        # np.clip keeps -1 from ever indexing categories.
        return np.where(codes >= 0, categories[np.clip(codes, 0, None)], "nan")
    raise RuntimeError(
        f"Unsupported encoding-type '{enc}' for obs column '{col}'."
    )


def read_obs_keys_h5py(path, sample_col):
    """
    Memory-lean (sample, barcode) key-set reader for one preprocessed view,
    using h5py only — no anndata open, so the backed-open eager `layers`
    materialization (anndata 0.12.x; ~6-10 GB for a full-cohort view) never
    happens and the key-set comparison stays at ~100-200 MB. Reproduces
    `cell_keys()` semantics exactly (str sample column, "_" separator, str
    obs index). Fails closed on layout/encoding surprises and on
    obs<->matrix row-count mismatches.
    """
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        obs_enc = obs.attrs.get("encoding-type")
        if isinstance(obs_enc, bytes):
            obs_enc = obs_enc.decode("utf-8")
        if str(obs_enc) != "dataframe":
            raise RuntimeError(
                f"{path}: obs encoding-type is '{obs_enc}', expected "
                "'dataframe' (anndata 0.12 h5ad layout)."
            )
        barcodes = _read_str_dataset(obs["_index"])
        samples = _read_obs_column(obs, sample_col)
        # Row-count guard against the counts source, mirroring
        # write_annotation_union's fallback (layers/counts else X) and CSR
        # checks: a corrupt obs<->matrix mismatch must fail loudly.
        if "layers/counts" in f:
            counts = f["layers/counts"]
        else:
            counts = f["X"]
        counts_enc = counts.attrs.get("encoding-type")
        if isinstance(counts_enc, bytes):
            counts_enc = counts_enc.decode("utf-8")
        if str(counts_enc) != "csr_matrix":
            raise RuntimeError(
                f"{path}: counts source encoding-type is '{counts_enc}', "
                "expected 'csr_matrix'."
            )
        indptr = counts["indptr"][:]
        n_rows = len(indptr) - 1
        if n_rows != len(barcodes):
            raise RuntimeError(
                f"{path}: obs/matrix row mismatch ({len(barcodes)} obs rows "
                f"vs {n_rows} counts rows)."
            )
    return {f"{s}_{b}" for s, b in zip(samples, barcodes)}


def read_obs_colnames_h5py(path):
    """
    Memory-lean obs column-name reader (h5py only, no anndata open) for the
    clean-entry check. Column names are the keys of the obs dataframe group
    minus the row-index key `_index`. Fails closed on layout surprises,
    mirroring read_obs_keys_h5py.
    """
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        obs_enc = obs.attrs.get("encoding-type")
        if isinstance(obs_enc, bytes):
            obs_enc = obs_enc.decode("utf-8")
        if str(obs_enc) != "dataframe":
            raise RuntimeError(
                f"{path}: obs encoding-type is '{obs_enc}', expected "
                "'dataframe' (anndata 0.12 h5ad layout)."
            )
        return {str(k) for k in obs.keys() if k != "_index"}


def read_feather_colnames(path):
    """
    Column names of one annotations feather (schema-only read). The R arrow
    workers write feather v2 (= Arrow IPC file format), so a RecordBatchFile
    reader exposes the schema without loading data; falls back to a full read
    for any other layout (feathers are small).
    """
    try:
        reader = pa.ipc.open_file(path)
        names = list(reader.schema.names)
        reader.close()
        return names
    except Exception:
        return list(pafeather.read_table(path).column_names)


def check_clean(path_data):
    """
    Clean-entry check for one dataset: a dataset entering annotation is
    "clean" iff every tier-matching obs column of every view is also a column
    of its own annotation feathers (i.e. produced by THIS pipeline); legacy
    columns are tier matches NOT present in the feathers. Returns an exit
    code: 0 = clean / no views, 2 = legacy found, 1 = error (raised).
    """
    h5ad_files = sorted(
        f for f in path_data.glob("*.h5ad") if not f.name.endswith("_raw.h5ad")
    )
    feather_files = sorted(path_data.glob("annotations_chunk_*.feather"))
    feather_cols = set()
    for ff in feather_files:
        feather_cols |= set(read_feather_colnames(ff))

    if not h5ad_files:
        print(f"No preprocessed .h5ad views in {path_data} (and "
              f"{len(feather_files)} feather files); nothing to check — clean.")
        return 0

    found_legacy = False
    for h5ad_file in h5ad_files:
        cols = read_obs_colnames_h5py(h5ad_file)
        tier1_matches = sorted(c for c in cols if any(re.search(p, c) for p in LEGACY_ANNOT_TIER1))
        tier2_matches = sorted(set(LEGACY_ANNOT_TIER2) & cols)
        legacy = sorted((set(tier1_matches) | set(tier2_matches)) - feather_cols)
        if not legacy:
            print(f"CLEAN: {h5ad_file.name} — {len(cols)} obs columns, "
                  f"tier matches all present in this run's feathers "
                  f"(T1: {len(tier1_matches)}, T2: {len(tier2_matches)})")
            continue
        found_legacy = True
        prev_annot = sorted(set(tier1_matches) - feather_cols)
        author_meta = sorted(set(tier2_matches) - feather_cols - set(tier1_matches))
        print(f"NOT CLEAN: {h5ad_file.name} — legacy annotation column(s): "
              f"{legacy}")
        print(f"  of which previous-annotation-like (Tier-1): {prev_annot}")
        print(f"  possibly author metadata (Tier-2 only): {author_meta}")
    if found_legacy:
        print(f"CLEAN-ENTRY CHECK FAILED for {path_data}: {len(feather_files)} "
              "feather file(s) on record; legacy columns will leak into "
              "annotation unless chunks are rebuilt and the dataset is "
              "re-annotated (worker wipe + merge tiered drop then scrub them).")
        return 2
    return 0


def verify_annotation_union(union_path, expected_nnz, sample_col):
    """
    Backed-mode self-check of a written annotation union: X must be a lazy
    backed CSR dataset (not eagerly materialized), its nnz must equal the
    source counts nnz, the file must carry no layers/obsp/obsm data, and the
    sample column must be present. Raises on any mismatch.
    """
    adata = ad.read_h5ad(union_path, backed="r")
    try:
        if type(adata.X).__name__ != "_CSRDataset":
            raise RuntimeError(
                f"Union X is '{type(adata.X).__name__}' (expected lazy backed "
                "CSR) — on-disk layout is not the minimal union."
            )
        if len(adata.layers) or len(adata.obsp) or len(adata.obsm):
            raise RuntimeError("Union unexpectedly contains layers/obsp/obsm data.")
        if sample_col not in adata.obs.columns:
            raise KeyError(f"Sample column '{sample_col}' not found in union obs.")
        with h5py.File(union_path, "r") as f:
            nnz = int(f["X/indptr"][-1])
        if nnz != expected_nnz:
            raise RuntimeError(
                f"Union X nnz ({nnz}) != source counts nnz ({expected_nnz})."
            )
    finally:
        adata.file.close()
    print(f"  Union self-check OK: lazy X ({nnz} nnz), no layers/obsp/obsm, "
          f"'{sample_col}' present.")


def write_annotation_union(source_path, union_path, sample_col):
    """
    Streaming, memory-lean annotation-union writer: copies the obs/var groups
    (datasets + attrs) from `source_path` and re-writes the raw counts as X
    (float64 CSR, same dtype as today's layers["counts"]), so the worker's
    astype("float64") path is a no-op. No layers/obsp/obsm/uns/varm/varp are
    copied. Writes atomically (tmp + os.replace) and self-verifies the result.
    """
    union_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = union_path.with_name(union_path.name + ".tmp")
    try:
        with h5py.File(source_path, "r") as src, h5py.File(tmp_path, "w") as dst:
            dst.attrs["encoding-type"] = src.attrs.get("encoding-type", "anndata")
            dst.attrs["encoding-version"] = src.attrs.get("encoding-version", "0.1.0")
            for grp in ("obs", "var"):
                src.copy(src[grp], dst, name=grp)
            if "layers/counts" in src:
                counts = src["layers/counts"]
            else:
                print(f"  No layers/counts in {source_path.name}; using its X as "
                      "the counts source (preprocessed views should always "
                      "carry a counts layer).")
                counts = src["X"]
            # Fail closed on non-CSR sources: preprocessed views are CSR by
            # construction (1.1.1_preprocess.py), so anything else means the
            # file was replaced/regenerated elsewhere. Copying a CSC layout
            # under the csr_matrix encoding-type below would silently corrupt
            # the union (the nnz/shape self-check cannot catch it).
            counts_enc = str(counts.attrs.get("encoding-type"))
            if counts_enc != "csr_matrix":
                raise RuntimeError(
                    f"{source_path.name}: counts source encoding-type is "
                    f"'{counts_enc}', expected 'csr_matrix'."
                )
            x = dst.create_group("X")
            x.attrs["encoding-type"] = "csr_matrix"
            x.attrs["encoding-version"] = "0.1.0"
            x.attrs["shape"] = counts.attrs["shape"]
            for ds_name in ("data", "indices", "indptr"):
                src.copy(counts[ds_name], x, name=ds_name)
            expected_nnz = int(counts["indptr"][-1])
        os.replace(tmp_path, union_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    verify_annotation_union(union_path, expected_nnz, sample_col)


def build_union(h5ad_files, union_path, sample_col):
    """
    Build the per-dataset union h5ad (dedup on (sample, barcode)).
    Always writes the minimal annotation union and returns `union_path`.
    """
    if len(h5ad_files) == 1:
        print(f"Single view: writing annotation union from {h5ad_files[0].name} "
              "(streaming copy; the view itself keeps its full layout).")
        write_annotation_union(h5ad_files[0], union_path, sample_col)
        return union_path

    print(f"Building union from {len(h5ad_files)} views: "
          + ", ".join(f.name for f in h5ad_files))

    # 1. Compare (sample, barcode) key sets using h5py-only obs reads (no
    #    anndata open: on anndata 0.12.x, backed open eagerly materializes
    #    the whole layers group into RAM, which OOMs the 4G srun on
    #    full-cohort views like Stephenson's batch-effect view). If one view
    #    already equals the union, stream it.
    key_sets = []
    for f in h5ad_files:
        keys = read_obs_keys_h5py(f, sample_col)
        key_sets.append(keys)
        print(f"  {f.name}: {len(keys)} cells")

    union_keys = set().union(*key_sets)
    largest_idx = int(np.argmax([len(k) for k in key_sets]))
    print(f"Union: {len(union_keys)} cells")

    if len(key_sets[largest_idx]) == len(union_keys):
        print(f"  Union == view {h5ad_files[largest_idx].name}; streaming copy "
              "(no in-memory concat).")
        write_annotation_union(h5ad_files[largest_idx], union_path, sample_col)
        return union_path

    # 2. Partial-overlap views: concat in memory and dedup. Memory-heavy; the
    #    srun allocation for this dataset must be raised accordingly.
    print("WARNING: views partially overlap; concatenating in memory. "
          "If this OOMs, raise the srun --mem for this dataset in "
          "1_prepare_chunks.sh.")
    adatas = [ad.read_h5ad(f) for f in h5ad_files]
    adata = sc.concat(adatas, join="outer", index_unique=None)
    del adatas
    keys = cell_keys(adata, sample_col)
    # keep first occurrence of each (sample, barcode) key; fancy indexing
    # already allocates fresh matrices, so no extra .copy() is needed
    _, first_idx = np.unique(keys, return_index=True)
    adata = adata[np.sort(first_idx)]

    # Minimal union layout: X = raw counts (CSR), obs/var only — mirrors the
    # streaming writer. anndata writes empty group stubs for emptied mappings,
    # so strip them on disk afterwards.
    counts = adata.layers.get("counts")
    if counts is None:
        counts = adata.X
        if counts is None:
            raise RuntimeError("No counts layer or X to use as the union's X.")
    if not sp.isspmatrix_csr(counts):
        counts = sp.csr_matrix(counts)
    adata.X = counts
    adata.layers = {}
    adata.obsm = {}
    adata.obsp = {}
    adata.uns = {}
    adata.varm = {}
    adata.varp = {}
    expected_nnz = int(counts.nnz)

    union_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(union_path))
    del adata
    with h5py.File(union_path, "r+") as f:
        for grp in ("layers", "obsm", "obsp", "uns", "varm", "varp"):
            if grp in f:
                del f[grp]
    verify_annotation_union(union_path, expected_nnz, sample_col)
    return union_path


def main():
    parser = argparse.ArgumentParser(description="Build per-dataset sample chunk files.")
    parser.add_argument(
        "--test",
        action="store_true",
        default=False,
        help="Test mode: 1 sample per chunk (production default: 2).",
    )
    parser.add_argument(
        "--check-clean",
        action="store_true",
        default=False,
        help="Clean-entry check only: report tier-matching obs columns per "
             "view that are NOT produced by this pipeline's feathers; exit 0 "
             "= clean, 2 = legacy found, 1 = error. No chunk/union writes.",
    )
    args = parser.parse_args()

    project_root = os.environ.get("PROJECT_ROOT")
    if not project_root:
        sys.exit("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")

    ds_name = os.environ.get("DS_NAME")
    if not ds_name:
        sys.exit("CRITICAL Error: DS_NAME not set. Ensure it is exported before calling this script.")

    sample_col = os.environ.get("SAMPLE_COLNAME")
    if not sample_col:
        sys.exit("CRITICAL Error: SAMPLE_COLNAME not set. Ensure it is exported before calling this script.")

    hpc_scratch_dir = os.environ.get("HPC_SCRATCH_DIR")
    if not hpc_scratch_dir:
        sys.exit("CRITICAL Error: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling this script.")

    path_data = Path(hpc_scratch_dir) / ds_name / "output"
    path_union_dir = Path(hpc_scratch_dir) / ds_name / "annotation_union"
    path_output_chunks = path_data / "chunks"

    print(f"Path is: {path_data}")

    if args.check_clean:
        try:
            code = check_clean(path_data)
        except Exception as e:
            print(f"ERROR: clean-entry check failed for {ds_name}: {e}")
            sys.exit(1)
        sys.exit(code)

    # Skip rds->h5ad conversion caches written into the output dir by
    # src.utils.py.preprocess_utils.load_single_input() (they lack the standardized sample
    # column and are not preprocessed views). The union h5ad lives in
    # annotation_union/, outside this glob.
    h5ad_files = sorted(
        f for f in path_data.glob("*.h5ad") if not f.name.endswith("_raw.h5ad")
    )
    print(f"Files found: {', '.join(f.name for f in h5ad_files)}")
    if not h5ad_files:
        sys.exit(f"CRITICAL Error: No preprocessed .h5ad files found in {path_data} "
                 "(run the preprocess array first).")

    # Defensive fail-closed completeness check (mirrors the bash guard in
    # 1_prepare_chunks.sh, mirroring 1.1.1_preprocess.py skip semantics):
    # every view output the preprocess array is expected to produce must
    # already exist. Catches bypasses/drift of the bash check — the dataset
    # lands in FAILED_DATASETS instead of building a partial union.
    datasets_json = os.environ.get("DATASETS_JSON_FILE")
    if not datasets_json:
        print("WARNING: DATASETS_JSON_FILE not set; skipping the expected-view "
              "check (source slurm_config.sh on HPC).")
    else:
        with open(datasets_json) as f:
            ds_entry = json.load(f).get(ds_name, {})
        expected = {
            v.get("output_file_name")
            for v in ds_entry.get("views", {}).values()
            if v.get("input_file_name") is not None and v.get("output_file_name")
        }
        missing = expected - {f.name for f in h5ad_files}
        if missing:
            sys.exit(
                f"CRITICAL Error: preprocessing incomplete for {ds_name}: missing "
                f"expected view file(s) {sorted(missing)} in {path_data} "
                "(run the preprocess array first)."
            )

    # Delete chunk file folder recursively to ensure a perfectly clean start
    shutil.rmtree(path_output_chunks, ignore_errors=True)
    path_output_chunks.mkdir(parents=True, exist_ok=True)

    # Delete any stale union from a previous run (e.g. after --force
    # preprocess re-runs regenerated the views).
    shutil.rmtree(path_union_dir, ignore_errors=True)
    path_union_dir.mkdir(parents=True, exist_ok=True)

    # Number of samples per chunk
    chunk_size = 1 if args.test else 2

    # Build the per-dataset minimal annotation union
    union_file = build_union(h5ad_files, path_union_dir / "union.h5ad", sample_col)

    # Read the union file in backed mode and extract its unique samples
    adata = ad.read_h5ad(union_file, backed="r")
    if sample_col not in adata.obs.columns:
        raise KeyError(f"Sample column '{sample_col}' not found in obs of {union_file}")
    file_samples = adata.obs[sample_col].astype(str).unique()
    adata.file.close()
    print(f"Union file {union_file.name} has {len(file_samples)} samples.")

    # Split the union's samples into consecutive groups of chunk_size
    sample_groups = [
        file_samples[i : i + chunk_size] for i in range(0, len(file_samples), chunk_size)
    ]

    # Write each group to a unique chunk file (global counter within the dataset)
    global_chunk_counter = 1
    for group in sample_groups:
        # We write the source file as the VERY FIRST line, followed by the sample IDs
        chunk_path = path_output_chunks / f"chunk_{global_chunk_counter}.txt"
        chunk_path.write_text(
            "\n".join([str(union_file.resolve())] + [str(s) for s in group]) + "\n"
        )

        global_chunk_counter += 1

    print(f"Successfully generated {global_chunk_counter - 1} total chunk files.")

    # Delete stale per-chunk annotation feathers from earlier runs (chunk
    # numbering is per-dataset and changes with chunk size / sample set, so a
    # rerun must not merge leftovers from a previous numbering). Only after
    # all chunks were generated successfully, and never in --test mode
    # (test runs must not destroy production annotation results).
    if not args.test:
        for stale in path_data.glob("annotations_chunk_*.feather"):
            stale.unlink()
            print(f"Removed stale annotations file: {stale.name}")

        # Also delete stale per-sample checkpoint intermediates
        # (output/annotation_tmp/): chunk numbering/sample sets change on
        # rebuild, so old intermediates must not be resumed by 2.1.1
        # (which maps sample_<NN>.feather to positions in the chunk file).
        stale_tmp = path_data / "annotation_tmp"
        if stale_tmp.exists():
            shutil.rmtree(stale_tmp, ignore_errors=True)
            print(f"Removed stale annotation checkpoints: {stale_tmp}")


if __name__ == "__main__":
    main()
