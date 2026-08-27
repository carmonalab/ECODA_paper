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
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


import anndata as ad
import h5py
import numpy as np
import pyarrow as pa
import pyarrow.feather as pafeather
import scipy.sparse as sp
import scanpy as sc
from src.utils.py.h5ad_source_identity import (
    read_obs_column as _read_obs_column,
    read_str_dataset as _read_str_dataset,
)


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
        if expected_nnz is not None and nnz != expected_nnz:
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
    tmp_path = union_path.with_name(f".{union_path.name}.tmp.{os.getpid()}")
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
    write_sidecar(union_path)


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
          "If this OOMs, raise the srun --mem for this dataset in the "
          "canonical Stage 4 submitter.")
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
    tmp_path = union_path.with_name(union_path.name + f".tmp.{os.getpid()}")
    try:
        adata.write_h5ad(str(tmp_path))
        if not tmp_path.is_file() or tmp_path.stat().st_size == 0:
            raise RuntimeError(f"empty annotation union: {tmp_path}")
        with h5py.File(tmp_path, "r+") as f:
            for grp in ("layers", "obsm", "obsp", "uns", "varm", "varp"):
                if grp in f:
                    del f[grp]
        os.replace(tmp_path, union_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    del adata
    verify_annotation_union(union_path, expected_nnz, sample_col)
    write_sidecar(union_path)
    return union_path


def parse_view_names(raw):
    if not raw:
        return []
    names = [part.strip() for part in raw.split(",")]
    if not names or any(not name for name in names):
        raise ValueError("annotation view selection contains an empty item")
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate annotation view: {raw}")
    return names


def selected_h5ad_files(ds_name, path_data, datasets_json, raw_views):
    with open(datasets_json) as handle:
        entry = json.load(handle).get(ds_name)
    if entry is None:
        raise KeyError(f"dataset {ds_name!r} is absent from datasets.json")
    declared = entry.get("views") or {}
    view_names = parse_view_names(raw_views)
    if not view_names:
        view_names = list(declared.keys())
    if not view_names:
        raise ValueError(f"dataset {ds_name} declares no annotation views")
    files = []
    for view_name in view_names:
        view = declared.get(view_name)
        if not view:
            raise ValueError(f"{ds_name}/{view_name} is not declared in datasets.json")
        output_name = view.get("output_file_name") or view.get("output_file")
        input_name = view.get("input_file_name") or view.get("input_file")
        if not input_name or not output_name:
            raise ValueError(f"{ds_name}/{view_name} has no input/output contract")
        path = path_data / output_name
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"selected preprocessed view is missing: {path}")
        files.append(path)
    return files


def validate_prepared(chunk_dir, union_path, sample_col, source_records=None):
    if not union_path.is_file() or union_path.stat().st_size == 0:
        raise ValueError(f"annotation union is missing/empty: {union_path}")
    validate_sidecar_record(union_path)
    verify_annotation_union(union_path, None, sample_col)
    if source_records is not None:
        metadata_path = union_path.parent.parent / "source_artifacts.json"
        if not metadata_path.is_file() or metadata_path.stat().st_size == 0:
            raise ValueError(f"source artifact metadata is missing: {metadata_path}")
        try:
            recorded = json.loads(metadata_path.read_text())
        except (OSError, ValueError) as exc:
            raise ValueError(f"source artifact metadata is invalid: {metadata_path}") from exc
        if recorded != source_records:
            raise ValueError(f"selected preprocessed source contents changed: {metadata_path}")
    union = ad.read_h5ad(union_path, backed="r")
    try:
        if sample_col not in union.obs.columns or union.n_obs <= 0:
            raise ValueError(f"annotation union has no non-empty {sample_col}: {union_path}")
        union_samples = set(union.obs[sample_col].astype(str).unique())
    finally:
        union.file.close()
    chunks = sorted(chunk_dir.glob("chunk_*.txt"))
    if not chunks:
        raise ValueError(f"no annotation chunks in {chunk_dir}")
    seen = []
    for chunk in chunks:
        lines = chunk.read_text().splitlines()
        if len(lines) < 2 or Path(lines[0]).resolve() != union_path.resolve():
            raise ValueError(f"chunk does not point at the current union: {chunk}")
        if any(not sample.strip() for sample in lines[1:]):
            raise ValueError(f"chunk contains a blank sample ID: {chunk}")
        seen.extend(lines[1:])
    if not seen or set(seen) != union_samples or len(seen) != len(set(seen)):
        raise ValueError(
            f"chunk sample coverage mismatch in {chunk_dir}: "
            f"chunks={len(set(seen))}, union={len(union_samples)}"
        )
    return len(chunks), len(union_samples)



def file_md5(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def artifact_record(path):
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "md5": file_md5(resolved),
        "size": resolved.stat().st_size,
    }
def validate_run_context(run_root, run_id, ds_name, scratch_dir):
    """Require a genuine Stage 4 run and its dataset owner before any access."""
    expected_root = (Path(scratch_dir) / "_ecoda_runs" / run_id).resolve()
    if run_root.resolve() != expected_root:
        raise ValueError(f"annotation run root is not canonical for {run_id}: {run_root}")
    metadata_path = run_root / "metadata"
    if not metadata_path.is_file() or metadata_path.stat().st_size == 0:
        raise ValueError(f"Stage 4 run metadata is missing: {metadata_path}")
    metadata = {}
    for line in metadata_path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            metadata[key] = value
    if metadata.get("STAGE") != "stage4" or metadata.get("RUN_ID") != run_id:
        raise ValueError(f"Stage 4 run metadata does not match {run_id}: {metadata_path}")
    safe_dataset = re.sub(r"[/,:	 |]", "_", ds_name)
    owner_dir = Path(scratch_dir) / "_ecoda_owners" / "stage4" / safe_dataset
    owner_path = owner_dir / "owner"
    if not owner_path.is_file() or owner_path.stat().st_size == 0:
        raise ValueError(f"Stage 4 dataset owner is missing: {owner_path}")
    owner = {}
    for line in owner_path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            owner[key] = value
    if (
        owner.get("RUN_ID") != run_id
        or owner.get("STAGE") != "stage4"
        or owner.get("KEY") != ds_name
        or owner.get("STATE") not in {"ACTIVE", "OK"}
    ):
        raise ValueError(f"Stage 4 dataset owner does not match {ds_name}: {owner_path}")



def validate_sidecar_record(path):
    sidecar = Path(f"{path}.md5")
    if not sidecar.is_file() or sidecar.stat().st_size == 0:
        raise ValueError(f"artifact checksum sidecar is missing: {sidecar}")
    fields = {}
    for line in sidecar.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            fields[key] = value
    if fields.get("MD5", "").lower() != file_md5(path):
        raise ValueError(f"artifact MD5 mismatch: {path}")
    if fields.get("SIZE") != str(path.stat().st_size):
        raise ValueError(f"artifact SIZE mismatch: {path}")
    try:
        recorded_path = Path(fields.get("PATH", "")).resolve()
    except (OSError, RuntimeError):
        recorded_path = None
    if recorded_path != path.resolve():
        raise ValueError(f"artifact PATH mismatch: {path}")


def write_sidecar(path):
    sidecar = Path(f"{path}.md5")
    atomic_text(
        sidecar,
        f"MD5={file_md5(path)}\nSIZE={path.stat().st_size}\nPATH={path}\n",
    )
def atomic_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        tmp.write_text(content)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def main():
    parser = argparse.ArgumentParser(description="Build per-dataset sample chunk files.")
    parser.add_argument("--test", action="store_true",
                        help="Test mode: 1 sample per chunk (production default: 2).")
    parser.add_argument("--check-clean", action="store_true",
                        help="Run the clean-entry check without writing chunks.")
    parser.add_argument("--validate-only", action="store_true",
                        help="Validate an existing run-owned union and chunk set.")
    parser.add_argument("--view", default=None,
                        help="Only use one declared output view.")
    parser.add_argument("--views", default=None,
                        help="Comma-separated declared output views.")
    parser.add_argument("--run-root", default=None,
                        help="Run-owned annotation root; avoids fixed-path races.")
    parser.add_argument("--force", action="store_true",
                        help="Rebuild only the selected run-owned dataset artifacts.")
    args = parser.parse_args()

    project_root = os.environ.get("PROJECT_ROOT")
    if not project_root:
        sys.exit("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")
    ds_name = os.environ.get("DS_NAME")
    if not ds_name:
        sys.exit("CRITICAL Error: DS_NAME not set. Ensure it is exported before calling this script.")
    sample_col = os.environ.get("SAMPLE_COLNAME") or "Sample"
    hpc_scratch_dir = os.environ.get("HPC_SCRATCH_DIR")
    if not hpc_scratch_dir:
        sys.exit("CRITICAL Error: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling this script.")
    datasets_json = os.environ.get("DATASETS_JSON_FILE") or str(Path(project_root) / "datasets.json")
    raw_views = args.views or args.view or os.environ.get("ANNOTATION_VIEWS") or ""
    if args.view and args.views:
        raise ValueError("use exactly one of --view and --views")
    path_data = Path(hpc_scratch_dir) / ds_name / "output"
    requested_run_root = args.run_root or os.environ.get("ANNOTATION_RUN_ROOT")
    if requested_run_root:
        run_id = os.environ.get("ANNOTATION_RUN_ID") or os.environ.get("RUN_ID")
        if not run_id:
            raise ValueError("run-owned annotation preparation requires ANNOTATION_RUN_ID")
        expected_root = (Path(hpc_scratch_dir) / "_ecoda_runs" / run_id).resolve()
        run_root = Path(requested_run_root).resolve()
        if run_root != expected_root:
            raise ValueError(f"annotation run root is not canonical for {run_id}: {run_root}")
        if not run_root.is_dir():
            raise ValueError(f"annotation run root is missing: {run_root}")
        dataset_root = run_root / "datasets" / ds_name
        path_union_dir = dataset_root / "union"
        path_output_chunks = dataset_root / "chunks"
        annotation_output_dir = dataset_root / "annotations"
    else:
        run_root = None
        path_union_dir = Path(hpc_scratch_dir) / ds_name / "annotation_union"
        path_output_chunks = path_data / "chunks"
        annotation_output_dir = path_data
    union_file = path_union_dir / "union.h5ad"
    if run_root is not None and not args.check_clean:
        validate_run_context(run_root, run_id, ds_name, hpc_scratch_dir)

    try:
        h5ad_files = selected_h5ad_files(ds_name, path_data, datasets_json, raw_views)
        source_records = [artifact_record(path) for path in h5ad_files]
        if args.check_clean:
            sys.exit(check_clean(path_data))
        if args.validate_only:
            n_chunks, n_samples = validate_prepared(
                path_output_chunks, union_file, sample_col,
                source_records if run_root is not None else None,
            )
            print(f"Prepared annotation artifacts valid: {n_chunks} chunks, {n_samples} samples")
            return

        path_output_chunks.parent.mkdir(parents=True, exist_ok=True)
        if args.force or not path_output_chunks.exists():
            shutil.rmtree(path_output_chunks, ignore_errors=True)
        path_output_chunks.mkdir(parents=True, exist_ok=True)
        if args.force or not path_union_dir.exists():
            shutil.rmtree(path_union_dir, ignore_errors=True)
        path_union_dir.mkdir(parents=True, exist_ok=True)
        if args.force and annotation_output_dir != path_data:
            shutil.rmtree(annotation_output_dir, ignore_errors=True)
        annotation_output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Path is: {path_data}")
        print(f"Selected views: {', '.join(f.name for f in h5ad_files)}")
        union_file = build_union(h5ad_files, union_file, sample_col)
        if run_root is not None:
            atomic_text(
                dataset_root / "source_artifacts.json",
                json.dumps(source_records, sort_keys=True, indent=2) + "\n",
            )
        adata = ad.read_h5ad(union_file, backed="r")
        if sample_col not in adata.obs.columns:
            adata.file.close()
            raise KeyError(f"Sample column '{sample_col}' not found in h5ad union")
        file_samples = adata.obs[sample_col].astype(str).unique()
        adata.file.close()
        chunk_size = 1 if args.test else 2
        for chunk_number, start in enumerate(range(0, len(file_samples), chunk_size), start=1):
            group = file_samples[start:start + chunk_size]
            chunk_path = path_output_chunks / f"chunk_{chunk_number}.txt"
            atomic_text(
                chunk_path,
                "\n".join([str(union_file.resolve())] + [str(sample) for sample in group]) + "\n",
            )
        if not args.test:
            # Production rebuilds own the run directory, never another run's
            # fixed scratch artifacts.
            for stale in annotation_output_dir.glob("annotations_chunk_*.feather"):
                stale.unlink()
            stale_tmp = annotation_output_dir / "annotation_tmp"
            if stale_tmp.exists():
                shutil.rmtree(stale_tmp, ignore_errors=True)
        n_chunks, n_samples = validate_prepared(
            path_output_chunks, union_file, sample_col,
            source_records if run_root is not None else None,
        )
        print(f"Successfully generated {n_chunks} chunks for {n_samples} samples.")
        print(f"ANNOTATION_CHUNK_DIR={path_output_chunks}")
        print(f"ANNOTATION_UNION={union_file}")
        print(f"ANNOTATION_FEATHER_DIR={annotation_output_dir}")
    except SystemExit:
        raise
    except Exception as exc:
        print(f"ERROR: chunk preparation failed for {ds_name}: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
