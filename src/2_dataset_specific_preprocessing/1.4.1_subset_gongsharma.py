"""
Per-sample 5000-cell cap of the staged GongSharma SoundLife h5ads (in place).

Fixes the preprocess-array OOM for `Gongsharma_cmv_young_males` (task 4, 128G
worker): the preprocess worker held both files + `sc.concat` (2.92M cells) and
densified `sc.pp.scale` matrices per HVG pass (~70 GB for 3000 top genes) —
over 128G. After this step, every sample (`columns.sample` =
`specimen.specimenGuid`) has at most 5000 cells.

Files overwritten IN PLACE (atomic: temp file + `os.replace`) in
${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data:
- SoundLife_YoungAdult_Male_CMVneg.h5ad  (1,712,244 -> 531,291 cells)
- SoundLife_YoungAdult_Male_CMVpos.h5ad   (1,206,761 -> 365,000 cells)
Total 896,291 cells / 180 samples / max 5000 per sample — matching the
historical NAS artifact Gongsharma_cmv_young_males.h5ad.

Sampling mirrors the historical `downsample_by_group` (git 3a4711e,
src/py/preprocess_gongsharma.qmd): a seed-42 `RandomState` created per file,
groups iterated in pandas `unique()` order of first appearance,
`RandomState.choice(pos, 5000, replace=False)` per group. A boolean keep-mask
(not the historical index-list subset, which reorders rows) preserves the
original row order while selecting the identical cell set.

Semantics / caveats:
- Idempotent: re-running on an already-capped file keeps every cell (every
  group is <= 5000), so the result is unchanged. No skip logic needed.
- RE-STAGING: `src/1_stage_data/1_stage_data.sh` rsyncs the NAS originals back
  over the capped files — re-run `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`
  (which runs this step again) before the preprocess array after any re-stage.
- ORDERING (race hazard): `1_submit_hpc.sh` submits this step FIRST and gates
  the CombinedPBMC step (`1.1_submit_combinedpbmc.sh`) behind it via
  `--dependency=afterok`, because that step reads the SAME staged files in
  backed mode — an in-place overwrite racing that read would nondeterminize
  the CombinedPBMC dataset.

Usage (HPC, via 1.4_submit_gongsharma.sh; pure Python, no R interop):
    ${PYTHON_BIN} 1.4.1_subset_gongsharma.py [--config_path ...] [--data_dir ...]
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.utils.py.datasets_io import read_datasets_json

DS_NAME = "Gongsharma_cmv_young_males"
EXPECTED_FILE_NAMES = [
    "SoundLife_YoungAdult_Male_CMVneg.h5ad",
    "SoundLife_YoungAdult_Male_CMVpos.h5ad",
]
# Cells expected after the cap, keyed by file name (NAS historical artifact
# Gongsharma_cmv_young_males.h5ad: 896,291 cells / 180 samples / max 5000).
EXPECTED_CELLS_AFTER_CAP = {
    "SoundLife_YoungAdult_Male_CMVneg.h5ad": 531291,
    "SoundLife_YoungAdult_Male_CMVpos.h5ad": 365000,
}
EXPECTED_TOTAL_AFTER_CAP = 896291
MAX_CELLS_PER_SAMPLE = 5000
SEED = 42


def build_keep_mask(adata, group_col, max_cells=MAX_CELLS_PER_SAMPLE, seed=SEED):
    """Boolean keep-mask mirroring the historical `downsample_by_group`.

    One seed-42 RandomState per file; groups in pandas `unique()` order of
    first appearance; `RandomState.choice` per group — identical selected
    cell SET as the historical index-list code, but order-preserving.
    """
    keep = np.zeros(adata.n_obs, dtype=bool)
    rng = np.random.RandomState(seed)
    for group in adata.obs[group_col].unique():
        pos = np.where(adata.obs[group_col] == group)[0]
        if len(pos) > max_cells:
            keep[rng.choice(pos, size=max_cells, replace=False)] = True
        else:
            keep[pos] = True
    return keep


def cap_file(path, group_col):
    print(f"=== Capping {path.name} ===", flush=True)
    adata = sc.read_h5ad(str(path), backed="r")
    n_before = adata.n_obs
    n_samples = adata.obs[group_col].nunique()
    max_before = int(adata.obs[group_col].value_counts().max())
    try:
        mask = build_keep_mask(adata, group_col)
        capped = adata[mask].to_memory()
    finally:
        try:
            adata.file.close()
        except Exception:
            pass
    n_after = capped.n_obs
    max_after = int(capped.obs[group_col].value_counts().max())

    tmp_path = path.with_name(path.name + ".capped_tmp.h5ad")
    capped.write_h5ad(str(tmp_path))
    os.replace(tmp_path, path)

    print(
        f"  {path.name}: {n_before:,} -> {n_after:,} cells "
        f"({n_samples} samples, max {max_before:,} -> {max_after:,}/sample)",
        flush=True,
    )
    expected = EXPECTED_CELLS_AFTER_CAP.get(path.name)
    if expected is not None and n_after != expected:
        print(
            f"  WARNING: expected {expected:,} cells after cap — check against "
            "the NAS historical artifact Gongsharma_cmv_young_males.h5ad.",
            flush=True,
        )
    return n_after


def main():
    default_data_dir = os.environ.get("HPC_SCRATCH_DIR")
    parser = argparse.ArgumentParser(
        description=(
            "Cap each sample of the staged GongSharma SoundLife h5ads at "
            "5000 cells (seed 42), overwriting them in place."
        )
    )
    parser.add_argument(
        "--config_path",
        default=os.environ.get("DATASETS_JSON_FILE") or "datasets.json",
        help="Path to datasets.json (default: ${DATASETS_JSON_FILE} or 'datasets.json' in CWD).",
    )
    parser.add_argument(
        "--data_dir",
        default=None if not default_data_dir else str(Path(default_data_dir) / DS_NAME / "data"),
        help=(
            f"Directory holding the staged input files "
            f"(default: ${{HPC_SCRATCH_DIR}}/{DS_NAME}/data)."
        ),
    )
    args = parser.parse_args()

    config_path = Path(args.config_path)
    if not config_path.is_absolute():
        config_path = Path.cwd() / config_path
    config = read_datasets_json(str(config_path))

    entry = config.get(DS_NAME)
    if entry is None:
        raise KeyError(f"Dataset {DS_NAME!r} not found in {config_path}.")

    file_names = entry.get("file_names")
    if not isinstance(file_names, list) or sorted(file_names) != sorted(EXPECTED_FILE_NAMES):
        raise ValueError(
            f"Unexpected file_names for {DS_NAME}: {file_names!r} — expected exactly "
            f"{EXPECTED_FILE_NAMES!r}. The staged inputs drifted; re-check datasets.json "
            "before capping."
        )
    group_col = entry.get("sample_col")
    if not group_col:
        raise ValueError(f"No columns.sample for {DS_NAME} in {config_path}.")

    if args.data_dir is None:
        parser.error(
            "--data_dir required (HPC_SCRATCH_DIR unset): pass the staged-data dir explicitly."
        )
    data_dir = Path(args.data_dir)
    missing = [n for n in file_names if not (data_dir / n).exists()]
    if missing:
        raise FileNotFoundError(
            f"Staged input(s) missing in {data_dir}: {missing} — "
            "run src/1_stage_data/1_stage_data.sh first."
        )

    total_after = 0
    for name in file_names:
        total_after += cap_file(data_dir / name, group_col)

    if total_after != EXPECTED_TOTAL_AFTER_CAP:
        print(
            f"WARNING: total after cap = {total_after:,} cells, expected "
            f"{EXPECTED_TOTAL_AFTER_CAP:,} (NAS historical artifact) — investigate "
            "before re-running the preprocess array.",
            flush=True,
        )
    else:
        print(
            f"Total after cap: {total_after:,} cells "
            f"(expected {EXPECTED_TOTAL_AFTER_CAP:,}) — OK.",
            flush=True,
        )


if __name__ == "__main__":
    main()
