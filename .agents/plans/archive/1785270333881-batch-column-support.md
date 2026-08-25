# Plan: Step 3d — Batch Column Support & Combined PBMC Dataset

## Goal

Enable `preprocess.py` to properly handle batch columns for `batch_effect_analysis` views, create the Joanito batch column, overhaul cell type annotation output, create the Combined PBMC dataset (with per-source gene standardization before merging), and update `datasets.json`.

## Prerequisites / Dependencies

- Step 1 (directory restructure + SLURM centralization) — **done**
- Step 3a (fix datasets.json reading in preprocess.py) is **not done** and must be included — preprocess.py currently crashes on the per-view datasets.json format (`json.load(f)["datasets"]` + `ds_info.get("file_name")`)

## Key Design Decisions

### 1. Joanito batch column: pre-create via one-time R script
`seqtec` is derived from the `dataset` column (CRC-SG1/KUL5 → "5' seq", rest → "3' seq"). Rather than adding derivation logic to preprocess.py, write `_create_joanito_batch_col.R` that modifies the .rds in-place.

### 2. Combined PBMC: standardize genes per-source THEN merge
The combine script loads each source dataset's raw data, applies `standardize_gene_symbols()` per-source (same logic as preprocess.py), THEN finds the common gene set, subsets, and concatenates. This avoids merging non-standardized gene names and makes the combined .h5ad immediately ready for preprocess.py.

### 3. Cell type annotation overhaul (prerequisite for unified layer2)
The cell type annotation pipeline is restructured:
- **Remove** the redundant `STACAS::StandardizeGeneSymbols()` call from `2.2_process_chunk.sh` (preprocess.py already does this)
- **Change output** from per-sample .rds files to per-chunk CSV files (cell_barcode, layer1, layer2, layer3)
- **Add merge-back script** that reads all chunk annotation CSVs and writes layer1/layer2/layer3 into the input .h5ad's `obs`
- This benefits both per-dataset analysis and the combined PBMC

### 4. Combined PBMC gets layer2 via annotation pipeline (not in combine script)
The combine script does NOT attempt to map cell types. Instead:
1. Combine script creates raw combined .h5ad (standardized genes, batch column, no cell types)
2. preprocess.py processes it (normalization, HVG, PCA, Harmony with batch_key="batch")
3. Overhauled cell type annotation runs on the combined preprocessed .h5ad → writes layer2 into it
4. ECODA analysis reads the annotated combined .h5ad

### 5. Zhu removed from datasets.json
Zhu is referenced nowhere except `batch_effect_analysis.rmd` (old-pipeline RDS path). Once CombinedPBMC exists, standalone Zhu is dead weight.

---

## Implementation Tasks

### Task A — Fix preprocess.py config reading (Step 3a)

**File**: `src/preprocess/1.2_preprocess.py`

Changes to `main()`:
1. Remove `["datasets"]` wrapper: `json.load(f)["datasets"]` → `json.load(f)`
2. Iterate per-view: move `input_file_name`/`output_file_name` reading inside the view loop using `view_info.get("input_file_name")` and `view_info.get("output_file_name")`
3. Skip if `input_file_name` is missing (print message, continue)
4. Use `view_info["output_file_name"]` as `processed_file_path` directly (instead of constructing from `file_name_no_ext`)
5. Handle array `input_file_name`: write helper `_load_input(input_file_name, base_path)` that:
   - If `str`: convert .rds→.h5ad (via R interop, with cached raw.h5ad) or load .h5ad directly
   - If `list`: for each entry, convert/load, then `sc.concat()` with `index_unique` to avoid obs_name collisions

Add `_load_input()` function near the top of the file (after the R interop block).

### Task B — Batch column reading in preprocess.py (Step 3d)

**File**: `src/preprocess/1.2_preprocess.py`

1. Read `batch_col` once per dataset from `ds_info.get("columns", {}).get("batch", sample_col)` (already at line 260 — unchanged)
2. Keep current logic (lines 290-293):
   ```python
   is_batch_view = view_name == "batch_effect_analysis" and ds_info.get("use_for_batch_effect", False)
   batch_key = batch_col if is_batch_view else None
   ```
3. Add validation: raise `ValueError` if `is_batch_view` is True but `batch_col` is not in `adata_view.obs.columns` (catches missing columns like an un-created Joanito `seqtec`)
4. **No change** to `process_view()` logic

### Task C — Create `_create_joanito_batch_col.R`

**File**: `src/preprocess/_create_joanito_batch_col.R`

```r
library(Seurat)

input <- "data/JoaI_2022_35773407_Nofilt_whole.rds"
seurat <- readRDS(input)

seurat$seqtec <- ifelse(
  seurat$dataset %in% c("CRC-SG1", "KUL5"),
  "5' seq",
  "3' seq"
)

saveRDS(seurat, input)
```

Run once before preprocess.py processes Joanito's batch_effect_analysis view.

### Task D — Overhaul cell type annotation output

The goal: eliminate per-sample .rds files, remove redundant STACAS gene standardization, write annotations back into the input .h5ad via a merge step.

#### D1. Remove STACAS gene standardization from `2.2_process_chunk.sh`
- Delete lines 203-207 (`STACAS::StandardizeGeneSymbols(...)`)
- The `geneRef` loading (line 166) can also be removed if not needed elsewhere in the script
- **Verify**: HiTME's `Run.HiTME()` does not implicitly depend on STACAS-standardized symbols

#### D2. Change chunk output from .rds to annotation feather files
In `2.2_process_chunk.sh`, replace lines 228 (per-sample `saveRDS`) and 228-270 (per-sample ECODA file saving) with:

After HiTME runs and `layer1`/`layer2`/`layer3` are in `seurat_obj@meta.data`:
```r
library(arrow)

# Extract annotations for this sample
annotations <- seurat_obj@meta.data[, c("layer1", "layer2", "layer3"), drop = FALSE]
annotations$cell_barcode <- rownames(annotations)

# Write per-chunk annotation file (one per chunk, feather format)
chunk_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
annot_file <- file.path(paths$path_output, paste0("annotations_chunk_", chunk_id, ".feather"))
write_feather(annotations, annot_file)
```

- Each chunk processes 5 samples sequentially and writes a single `.feather` file per chunk
- Remove `saveRDS(seurat_obj, file = processed_file_path)` (line 228) and the ECODA section (lines 230-270)
- `get_celltype_comp` parameter and `path_output_ecoda` can be deprecated

#### D3. Add merge-back script

**File**: `src/cell_type_annotation/3_merge_annotations.py`

Logic (uses feather for fast I/O):
```python
import pandas as pd
import anndata as ad
import glob

h5ad_path = ...  # same input as cell type annotation
output_path = ...  # same path, overwrite with annotations merged in

# Read all chunk annotation feather files
annot_files = sorted(glob.glob(f"{annot_dir}/annotations_chunk_*.feather"))
annotations = pd.concat(
    [pd.read_feather(f) for f in annot_files], ignore_index=True
)
annotations = annotations.set_index("cell_barcode")

# Load h5ad
adata = ad.read_h5ad(h5ad_path)

# Merge annotations into obs
adata.obs = adata.obs.join(annotations, how="left")

adata.write_h5ad(output_path)
```

This script runs on the HPC login node after all annotation chunks complete.

### Task E — Create `_create_combinedpbmc_dataset.py`

**File**: `src/preprocess/_create_combinedpbmc_dataset.py`

**Key change from original plan**: standardize gene symbols per-source BEFORE finding common genes and merging.

Logic:

```python
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys

# Add parent dir to path so we can import from preprocess.py
sys.path.insert(0, str(Path(__file__).parent))
from preprocess_1_2 import standardize_gene_symbols, _load_ensembl105_map

def load_and_prepare_source(ds_name, view_info, base_path):
    """Load raw data, apply filters, add metadata, standardize genes."""
    # 1. Load via _load_input logic (handles .rds and array .h5ad)
    ...

    # 2. Apply subset_vars
    ...

    # 3. Standardize gene symbols (key refinement)
    standardize_gene_symbols(adata)

    # 4. Add batch and condition metadata
    ...

    return adata

# Main flow:
# 1. Read datasets.json for Stephenson, GongSharma, Zhu configs
# 2. Load + prepare each source
# 3. Find common genes from standardized var_names
# 4. Subset each to common genes
# 5. Ensure unique sample IDs (prefix with dataset key)
# 6. sc.concat()
# 7. Save to data/combined_pbmc_batch_effect_analysis.h5ad
```

Handling per-source:
- **Stephenson**: Load from .rds (R interop), apply batch_effect_analysis `subset_vars`, standardize genes, add `batch`="Stephenson", `cond`=`Status`, keep raw counts
- **GongSharma**: Load from array .h5ad files, concat, `np.random.seed(123)` → pick 15 random samples, standardize genes, add `batch`="GongSharma", `cond`="Healthy"
- **Zhu**: Load from .rds (R interop), standardize genes, add `batch`="Zhu", `cond`="Healthy"
  - **New**: Do NOT attempt cell type column mapping. Cell types will come from annotation pipeline (Task D) running on the combined preprocessed .h5ad later.
  - Only keep obs columns: `Sample`, `batch`, `cond` (no cell type column in combine script)

Clean only the three columns above from obs, drop the rest (to minimize file size).

### Task F — Update `datasets.json`

Changes:
1. **Stephenson**: Add `"batch": "Site"` to `columns`
2. **Joanito**: Add `"batch": "seqtec"` to `columns`
3. **Add CombinedPBMC** entry (no cell_type columns since they come from annotation pipeline):
   ```json
   "CombinedPBMC": {
     "display_name": "Combined PBMC (Stephenson, GongSharma, Zhu)",
     "folder_name": null,
     "tissue": "Blood",
     "normal_tissue": true,
     "use_for_benchmark": false,
     "use_for_batch_effect": true,
     "columns": {
       "sample": "Sample",
       "label": "cond",
       "batch": "batch",
       "cell_type_low_res": null,
       "cell_type_high_res": null
     },
     "views": {
       "batch_effect_analysis": {
         "input_file_name": "combined_pbmc_batch_effect_analysis.h5ad",
         "output_file_name": "combined_pbmc_batch_effect_analysis_batch_effect_analysis_ECODAprocessed.h5ad",
         "subset_vars": {}
       }
     }
   }
   ```
   Note: `cell_type_high_res` is null because cell types will be added by Task D (annotation pipeline) post-preprocessing.
4. **Remove Zhu** entry entirely

### Task G — R-side: Update `read_datasets_json()` if needed

**File**: `src/utils/datasets_io.R`

Add `batch_col` to the returned entry list:
```r
entry <- list(
  ...
  batch_col = ds[["columns"]][["batch"]]
)
```

---

## Execution Order

```
Task A+B (preprocess.py)         ← must go first (currently crashes)
    │
    ├── Task C (Joanito batch col)     ← can run anytime before Joanito is preprocessed
    ├── Task D (cell type annotation)  ← independent; needed before combined PBMC analysis but not before creation
    ├── Task E (combine script)        ← depends on A+B (uses standardize_gene_symbols import)
    └── Task F (datasets.json)         ← after E (file must exist for CombinedPBMC entry)
         │
         └── Task G (read_datasets_json) ← independent, can be anytime after A+B
```

Post-combined-PBMC flow:
1. `_create_combinedpbmc_dataset.py` → data/combined_pbmc_batch_effect_analysis.h5ad
2. `preprocess.py` processes CombinedPBMC → combined preprocessed .h5ad (with Harmony by batch)
3. Cell type annotation (overhauled) on combined preprocessed .h5ad → writes layer2 into obs
4. R notebooks load annotated combined .h5ad for ECODA, GloScope, MrVI, PILOT-GM-VAE

---

## Validation

| Task | Validation |
|---|---|
| A | `python src/preprocess/1.2_preprocess.py` processes a dataset (e.g., Adams, Lee) without crashing |
| A | GongSharma array input loads, concatenates correctly |
| B | Stephenson batch_effect_analysis: `batch_key="Site"`, Harmony runs, outputs `X_pca_harmony_...` |
| B | Joanito batch_effect_analysis: errors with clear message if `seqtec` missing; passes after Task C |
| C | `readRDS(...)$seqtec` has only values "5' seq" / "3' seq" |
| D | `2.2_process_chunk.sh` runs without STACAS call; output is CSV not .rds |
| D | `3_merge_annotations.py` adds layer1/layer2/layer3 columns to .h5ad obs |
| E | Combined .h5ad loads, `var_names` are standardized (ensembl105), shape is correct |
| E | GongSharma contributes exactly 15 samples |
| E | No cell type column in obs (comes from annotation pipeline later) |
| F | preprocess.py processes CombinedPBMC with `batch_key="batch"`, produces Harmony output |
| F | `'Zhu' not in json.load(open('datasets.json'))` |
| G | `read_datasets_json("datasets.json", "batch_effect_analysis")` includes `batch_col` |

---

## Risk Register

| Risk | Mitigation |
|---|---|
| **Combine script imports from preprocess.py** — preprocess.py imports rpy2 which may not be available in all contexts | Move `_load_ensembl105_map()` and `standardize_gene_symbols()` to a shared module `src/gene_utils.py`, import from both scripts |
| **HiTME depends on STACAS-standardized symbols** (despite preprocess.py already standardizing) | Test on one dataset: run preprocess.py then annotation without STACAS; compare layer2 output with STACAS version. If HiTME needs STACAS, add a flag to skip only the redundant standardization |
| **Combined PBMC gene intersection too small** | If < 5000 common genes after standardization, use union instead; preprocess.py's `filter_genes` will drop dataset-specific genes with low expression |
| **Zhu raw data completely unannotated** | Not a problem — cell types come from annotation pipeline (Task D) running on the combined preprocessed .h5ad. No cell type column needed in combine script |
| **Feather write must not conflict across chunks** | Each chunk writes its OWN `.feather` file (`annotations_chunk_<task_id>.feather`); no concurrent writes to the same file. Chunk task IDs are unique by design |
| **Merge-back script runs on potentially large .h5ad** | Must load full .h5ad into memory. If too large, use `backed=True` and write in chunks; but preprocessed .h5ad sizes are typically < 1GB |
