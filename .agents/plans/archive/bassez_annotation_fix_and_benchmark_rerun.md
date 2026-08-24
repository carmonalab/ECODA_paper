# Execution & Parallelization Plan: Bassez Cell Type Fix + Commit 76be859 Rollout

> **Status update (2026-08-24):** This original proposal is retained for context
> but is superseded for execution by
> `.kilo/plans/1787522239980-bassez-annotation-benchmark-rerun.md`. Use the
> corrected plan for serialization, canonical paths, atomic writes, checksum
> validation, and render safeguards.

## Summary of Workflow Architecture

1. **Worker-Node In-Place Patching for Bassez (Prerequisite):**
   - **NO calculation runs on the login node.**
   - Staged Seurat RDS: Submitted to a worker node via `sbatch src/2_dataset_specific_preprocessing/1.6_submit_bassez.sh` (64GB RAM, 4 CPUs).
   - Preprocessed `.h5ad`: Executed on a compute node via `srun --mem=64G` inline command.
   - Bypasses re-running Stage 3 (preprocessing/HVG/PCA/Harmony) and Stage 4 (chunking/automated annotation arrays), saving hours of computation while preserving all embeddings and automated annotations.

2. **Bassez Method Scope in Stage 5:**
   - Because `cell_type_low_res` (`cellType`) did not change, methods that only use low-res annotations (`mrvi`, `scitd`) or whole-sample embeddings (`gloscope`, `mofa`) do **not** need to be re-run for Bassez.
   - Only methods consuming `cell_type_high_res` (`cellSubType`) need re-running for Bassez:
     - **Python (GPU/CPU):** `scpoli`, `pilot`, `qot`, `pilotgm` (`--methods scpoli,pilot,qot,pilotgm`).
     - **R (CPU):** `composition` (recomputes `ECODA_authors_HR*`, `GloProp`, `Freq_highres`, etc.) and `pseudobulk` (`Pseudobulk_CT_HR_*`).
     - **Trans & Zeroimp:** `trans` and `zeroimp` sensitivity analyses.

3. **Commit `76be859` Scope (`valid_mask` NA / Unknown Filtering across ALL Datasets):**
   - Commit `76be859cb3cc98cf3c5feb1f472a70a276bce4f8` introduced a `valid_mask` (`!is.na(cts) & !cts %in% c("NA", "nan", "None", "Unknown") & cts != ""`) in:
     - `get_ct_comp_df()` & `get_ct_comp_df_seurat()` (`src/utils/seurat_utils.R`)
     - `process_gloprop_fig()` (`src/5_run_benchmark_methods/benchmark_methods_r.R`)
   - Prevents unassigned cells (`"NA"`, `"Unknown"`, `"None"`, `"nan"`) from polluting sample-by-celltype composition tables across **all** benchmark datasets (author annotations with "Unknown" in `Kim`, `Wu`, `Pelka`, `Adams`, etc., and automated annotations from `scATOMIC`/`HiTME`).
   - Requires re-running `composition`, `trans`, and `zeroimp` arrays across **ALL benchmark datasets**.

---

## Execution Structure: Sequential vs. Parallel Dependency Graph

```mermaid
flowchart TD
    subgraph Phase 1: Prerequisite [Phase 1: Compute Node Prerequisite (Sequential)]
        P1A["1. sbatch 1.6_submit_bassez.sh\n(Patches staged Seurat .rds on 64GB node)"]
        P1B["2. srun --mem=64G python -c '...'\n(Patches preprocessed .h5ad on 64GB node)"]
        P1A --> P1B
    end

    subgraph Phase 2: Parallel [Phase 2: Independent Benchmark Arrays (Run in Parallel Shells)]
        direction LR
        J_PY["Terminal A: Python Array (Bassez)\nscpoli, pilot, qot, pilotgm"]
        J_COMP["Terminal B: R Composition Array (All Datasets)\ncomposition"]
        J_PB["Terminal C: R Pseudobulk Array (Bassez)\npseudobulk (CT-HR)"]
        J_TZ["Terminal D: Trans & Zeroimp Array (All Datasets)\ntrans, zeroimp"]
    end

    subgraph Phase 3: Validation [Phase 3: Downstream Validation & Render (Sequential)]
        P3["Local macOS Terminal:\npixi run rmarkdown::render('notebooks/benchmark_analysis.rmd')"]
    end

    P1B --> J_PY
    P1B --> J_COMP
    P1B --> J_PB
    P1B --> J_TZ
    J_PY --> P3
    J_COMP --> P3
    J_PB --> P3
    J_TZ --> P3
```

---

## Step-by-Step Execution Guide

### Phase 1: Prerequisite In-Place Patch on Compute Node (Sequential)

> [!IMPORTANT]
> Run these on `bamboo`. Do **not** proceed to Phase 2 until both commands finish with success.

```bash
# On bamboo terminal:
source src/slurm_config.sh
cd "${PROJECT_ROOT}"

# 1. Update the staged Seurat RDS object on worker node (64GB RAM):
sbatch src/2_dataset_specific_preprocessing/1.6_submit_bassez.sh

# (Wait for job completion via: squeue -u $USER)

# 2. Update the preprocessed .h5ad directly on compute node (64GB RAM):
srun --partition="${SLURM_PARTITION}" --mem=64G --cpus-per-task=4 ${PYTHON_BIN} -c '
import anndata as ad
import numpy as np
import os, sys

scratch_dir = os.environ.get("HPC_SCRATCH_DIR", "")
h5ad_path = f"{scratch_dir}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad"
if not os.path.exists(h5ad_path):
    h5ad_path = f"{scratch_dir}/BassezA_2021_33958794/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad"

print(f"Reading {h5ad_path} ...")
adata = ad.read_h5ad(h5ad_path)

raw_subtype = adata.obs["cellSubType"].astype(str)
is_na = adata.obs["cellSubType"].isna() | raw_subtype.isin(["NA", "nan", "None", "Unknown", ""])
n_na = is_na.sum()
print(f"Filling {n_na}/{len(adata)} missing cellSubType entries with cellType...")

adata.obs["cellSubType"] = raw_subtype
adata.obs.loc[is_na, "cellSubType"] = adata.obs.loc[is_na, "cellType"].astype(str)
adata.obs["cellSubType"] = adata.obs["cellSubType"].astype("category")

adata.write_h5ad(h5ad_path)
print("Successfully updated Bassez .h5ad in-place on compute node!")
'
```

---

### Phase 2: Independent Benchmark Arrays (Run in Parallel Across 4 Shells)

Once Phase 1 completes, the following 4 commands are completely independent (distinct output files, distinct log files, separate array IDs). They can be submitted **concurrently in 4 separate terminal sessions** on `bamboo`.

#### Terminal A (Python GPU/CPU Array — Bassez High-Res Methods)
Submits `scpoli`, `pilot`, `qot`, and `pilotgm` on GPU/CPU for Bassez:
```bash
# Terminal A (bamboo):
source src/slurm_config.sh
cd "${PROJECT_ROOT}"

./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name Bassez \
  --methods scpoli,pilot,qot,pilotgm \
  --force
```

#### Terminal B (R Composition Array — ALL Benchmark Datasets)
Applies the `76be859` NA-filtering fix across all datasets (including updated Bassez):
```bash
# Terminal B (bamboo):
source src/slurm_config.sh
cd "${PROJECT_ROOT}"

./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --methods composition \
  --force
```

#### Terminal C (R Pseudobulk Array — Bassez CT-HR Combos)
Recomputes cell-type specific pseudobulks (`Pseudobulk_CT_HR_*`) for Bassez:
```bash
# Terminal C (bamboo):
source src/slurm_config.sh
cd "${PROJECT_ROOT}"

./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name Bassez \
  --methods pseudobulk \
  --force
```

#### Terminal D (Trans & Zero-Imputation Arrays — ALL Benchmark Datasets)
Recomputes transformation and zero-imputation robustness checks across all datasets:
```bash
# Terminal D (bamboo):
source src/slurm_config.sh
cd "${PROJECT_ROOT}"

./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh \
  --force
```

---

### Phase 3: Validation & Publication Figure Rendering (Sequential)

Once all arrays from Phase 2 have completed and synchronized back to the NAS:

#### 1. Sanity Check Bassez Clean Annotations
```bash
# On bamboo:
source src/slurm_config.sh
${PYTHON_BIN} -c '
import anndata as ad
import os
scratch = os.environ.get("HPC_SCRATCH_DIR", "data")
adata = ad.read_h5ad(f"{scratch}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad", backed="r")
obs = adata.obs
print("Total cells:", len(obs))
print("NA/Unknown in cellSubType:", (obs["cellSubType"].isna() | obs["cellSubType"].astype(str).isin(["NA", "Unknown", "None", ""])).sum())
print("Unique cellSubTypes:", obs["cellSubType"].nunique())
'
```

#### 2. Render Benchmark Analysis Notebook (Local macOS)
```bash
# On local macOS terminal:
pixi run Rscript -e "rmarkdown::render('notebooks/benchmark_analysis.rmd')"
```
Inspect updated publication figures: `Figure 2A`, `Supp fig 2`, `Supp fig 15`, and the transformation/zero-imputation figures.

---

## Current Implementation Status

The source-hardening portion is implemented in the files listed by the current
working-tree diff, with focused regression coverage in
`tests/test_bassez_and_benchmark_regressions.R`. Local R, Pixi R, parse, Bash
syntax, and whitespace checks pass. HPC patching, benchmark arrays, NAS sync,
local mirror refresh, notebook rendering, commit, and push have not been run;
follow the corrected `.kilo` plan before starting those steps.

> **Execution status (2026-08-24 01:10):** Source hardening was committed as
> `7a2d591` and pushed to `origin/master`; the bamboo clone was fast-forwarded
> to that commit and preflight gates passed. The initial debug-cpu submission
> (`4334701`) was pending with `PartitionTimeLimit` and was cancelled before
> starting. The replacement shared-cpu job (`4334702`) completed successfully
> in 8m30s. Its log confirms 162,753 of 226,454 missing `cellSubType` values
> were filled, subtype cardinality changed from 33 to 40, and the validated
> RDS was atomically installed at the canonical scratch path. The processed
> `.h5ad`, NAS synchronization, benchmark reruns, local mirror refresh, and
> notebook render remain pending; do not begin them until explicitly resumed.
