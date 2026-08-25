# TODO.md Cleanup + Pipeline Completion Ordering

## Goal

Rewrite `TODO.md`: remove completed/stale sections, merge duplicates, and re-order the
remaining work as:

1. **Phase 1 — Complete pipelines src/1–4** (agent, no HPC needed): finish debug/test
   mode, switch to per-dataset annotation + integrate `3_merge_annotations.py`, audit
   legacy `Preprocess_datasets.Rmd`, audit `.Rprofile`/`.Renviron`.
2. **Phase 2 — Debug run on HPC** (explicitly requires user: connect NAS + log in to HPC):
   staged run of src/1 → src/4 on the `_debug` dataset.
3. **Phase 3 — src/5_run_benchmark_methods** (implementation; see `notebooks/benchmark_analysis.rmd`).
4. Phase 4 (batch effect analysis) + human-managed tasks as later sections.

## Verified current state (facts checked against code)

- `3_merge_annotations.py` exists (orphan — no bash wrapper; manual CLI only); joins
  annotation feathers to h5ad on `(sample, barcode)`, so it already works per-view when
  the feather set is per-dataset.
- No `_create_debug_dataset.R`, no `data/debug/`, no `_debug` in datasets.json (14 keys).
- `1_stage_data.sh` has NO dataset filter; `1_submit_hpc_array.sh` (preprocess) has NO
  single-dataset mode. `1_prepare_chunks.sh` supports `test <DS>` + `--test`;
  `2_submit_hpc_array.sh` supports single-dataset `<DS>`.
- `1.1.1_preprocess.py` has the "Already processed" skip (no `--force`).
- Views per dataset: only Stephenson has 2 views (benchmark ⊂ batch superset); Joanito +
  CombinedPBMC batch-only; Zhu view-less. Future multi-view datasets (PILOT-GM-VAE set)
  have unknown overlap relationships.
- **Kfoury gap**: `datasets.json` declares `columns.cell_type_low_res: "cells_lowres"`, but NO code in `src/` creates `cells_lowres` (legacy-only: `Preprocess_datasets.Rmd:444`). → must be ported.
- `Preprocess_datasets.Rmd` is at REPO ROOT (905 lines; legacy Seurat pipeline).
- `.Rprofile` (621 B): sets `RETICULATE_PYTHON` to pixi python when `PIXI_PROJECT_ROOT`
  is set, IDE fallbacks otherwise. **Influences pipelines**: `2.1.1_process_chunk.R:28-30`
  uses `reticulate` + `anndata`; `src/utils/imports.R:27` requires reticulate. Sourced
  because all HPC R runs start from `${PROJECT_ROOT}` (AGENTS.md convention).
- `.Renviron` (17 B): `R_MAX_VSIZE=200Gb` — unix virtual-memory cap, affects all R
  pipeline processes (annotation workers, pseudobulk, DESeq2).
- `1.2_benchmark_methods_py.qmd` is the only file in `run_python_sample_embedding_methods/`;
  `benchmark_analysis.rmd` still reads `_ECODAprocessed.rds` inputs.

## Decisions

- **USER-APPROVED**: register `_debug` in `datasets.json` (views: `benchmark_analysis` +
  `batch_effect_analysis`). Convention: default-all script loops skip keys starting with
  `_` unless explicitly requested via `--ds_name`.
- **USER-DECIDED**: annotation per-dataset, not per-view (see Phase 1.2).
- **USER-DECIDED**: `_create_debug_dataset.R` → decimal naming
  `src/2_dataset_specific_preprocessing/1.3.1_create_debug_dataset.R` (local script, NOT
  part of the `1_submit_hpc.sh` dispatcher glob — no `1.3_submit_*.sh` wrapper; document).
- **Recommended defaults (implementation agent follows unless user objects)**:
  - Per-dataset annotation approach: build a per-dataset union h5ad (concat all views,
    dedup on `(sample, barcode)`), write it OUTSIDE the synced output dir
    (`${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union/union.h5ad`) → no glob/NAS-sync
    interference; chunk + annotate it once; merge feathers into EACH view h5ad (existing
    join logic); delete union file after merge. Robust for future multi-view datasets
    with unknown overlap. (Fallback: chunk the largest view if union is infeasible.)
  - Merge integration: dedicated `3.1_submit_merge.sh` wrapper (per-dataset `<DS>`,
    srun, `PYTHON_BIN`, after annotation array + NAS sync; loops over view h5ads; rsyncs
    merged h5ads back to NAS).
  - Item E: add `--force` flag to `1.1.1_preprocess.py` (needed for debug re-runs).
  - Completed-history: remove all DONE sections from TODO.md entirely (git history
    preserves); keep a one-line changelog pointer.
  - `.Rprofile`: KEEP (needed by reticulate in `2.1.1_process_chunk.R`); verify
    `Sys.which("python")` under `pixi run` resolves to `.pixi/envs/default/bin/python`;
    document in AGENTS.md/README.
  - `.Renviron`: KEEP unless evidence of redundancy; verify under pixi R (pixi sets no
    R_MAX_VSIZE) whether any R step hit allocation errors; update value to match largest
    worker allocation + headroom, or delete if confirmed unnecessary (ask user).

## Phase 1 — Complete pipelines (src/1–4) [agent, no HPC needed]

### 1.1 Debug/test execution mode
- [ ] `src/2_dataset_specific_preprocessing/1.3.1_create_debug_dataset.R`: read Joanito
  input (per datasets.json `input_file_name`), select 5 samples covering both biological
  conditions and batches, subset 500 cells/sample (random), minimal obs cols (incl.
  `seqtec`, `Site`, sample/patient cols), write `data/debug/*.rds` + `*.h5ad`.
- [ ] datasets.json: add `_debug` entry (USER APPROVED; batch col `Site`, label cols
  following Joanito conventions, views → debug input/output files).
- [ ] `src/1_stage_data/1_stage_data.sh`: `--ds_name` filter; default-all skips `_*`
  keys unless explicitly requested; document that the debug file must exist on NAS
  (user places it) or be staged manually.
- [ ] `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`: single-dataset mode
  (`--ds_name` arg → map DS→jq index, 1-task array); default-all skips `_*`.
- [ ] `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`: `--force` flag to bypass the
  "Already processed" skip.

### 1.2 Per-dataset annotation + merge-back (replaces per-view; USER-DECIDED)
- [ ] `src/4_cell_type_annotation/1.1_prepare_chunks.py`: one chunk set per DATASET:
  build union h5ad (concat all view h5ads in `${HPC_SCRATCH_DIR}/${DS_NAME}/output`,
  dedup `(sample, barcode)`) → `${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union/union.h5ad`
  (outside output dir; not synced); chunk the union (line 1 of chunk file pins the union
  path); keep `--test` mode (1 sample/chunk).
- [ ] `3_merge_annotations.py`: polish CLI (argparse `--h5ad-path`, `--annot-dir`
  default = h5ad parent, `--output-path`); logic unchanged (SAMPLE_COLNAME env-driven
  `(sample, barcode)` join + whitelist obs subsetting); runs once per view h5ad of the
  dataset.
- [ ] New `3.1_submit_merge.sh` (per-dataset `<DS>`, srun, `PYTHON_BIN`): loop over the
  dataset's view h5ads, merge feathers, delete `annotation_union/union.h5ad`, rsync
  merged h5ads to `${NAS_TARGET_DIR}/${DS_NAME}/output/`.
- [ ] `2_submit_hpc_array.sh` / `2.1_run_worker.sh`: unchanged (chunk-file driven).
  Verify chunk manifest/feather naming still consistent (global counter across datasets).
- [ ] Docs: ARCHITECTURE.md (union + merge as wrapped stages; NAS sync unchanged),
  AGENTS.md repo-structure note.

### 1.3 Legacy `Preprocess_datasets.Rmd` audit (repo root; delete when done)
- [ ] Enumerate legacy steps; map to new pipeline (draft checklist for TODO.md):
  - scGate models + ProjecTILs loading (lines 408–427) → scGate: ported (scGate DB
    cache in `4_cell_type_annotation`); ProjecTILs: VERIFY still used by
    `2.1.1_process_chunk.R` or obsolete.
  - "Create low res cell types for Kfoury" (line 444) → **MISSING**: `cells_lowres`
    (author-CT collapse → Tcells/NKcells/Bcells/MoMac/DCcells) declared in datasets.json
    but never created in the new pipeline. **CHECK WITH USER**: port as dataset-specific
    step (`src/2_dataset_specific_preprocessing/`, e.g. `1.4.1_create_kfoury_lowres_ct.R`,
    run via `1_submit_hpc.sh`) vs fold into `1.1.1_preprocess.py`.
  - gene symbols (line 483) → ported (`src/gene_utils.py`); preprocess loop (line 489)
    → ported (`1.1.1_preprocess.py`); HiTME/scATOMIC whitelist (lines 501–518) → ported
    (`3_merge_annotations.py`).
  - "Create clearcut age for GongSharma" (line 826) → handled differently: current
    GongSharma condition (cmv/young/males) fixed via file selection; age split belongs
    to future other-subsetting conditions (draft `preprocess_gongsharma.qmd`).
    **CHECK WITH USER**: drop vs defer.
  - "Export datasets without author annotation" (line 868; Lee/Zhang for scPoli) →
    likely obsolete (new pipeline annotates all datasets; scPoli consumes preprocessed
    h5ad). **CHECK WITH USER**: drop vs adapt.
- [ ] Implement missing items (after user confirmation); keep `Figure_workflow_schematic.Rmd`
  untouched (out of scope).
- [ ] When audit complete + user OK: delete `Preprocess_datasets.Rmd`; update docs
  (AGENTS.md, ARCHITECTURE.md) if they reference it.

### 1.4 `.Rprofile` / `.Renviron` audit
- [ ] `.Rprofile` (621 B): KEEP (reticulate needed by `2.1.1_process_chunk.R` +
  `imports.R`). Verify under `pixi run Rscript`: `Sys.which("python")` →
  `.pixi/envs/default/bin/python`; confirm fallback branches still correct
  (Windows branch irrelevant on HPC/macOS — simplify if desired); document the wiring.
- [ ] `.Renviron` (`R_MAX_VSIZE=200Gb`): verify whether any R pipeline step needs it
  (annotation workers, pseudobulk/DESeq2 on cluster allocations); if needed keep (update
  value/comment to match allocations + headroom), else delete (ask user).

## Phase 2 — Debug run on HPC [REQUIRES USER: connect NAS + log in to HPC]

Prereqs (explicitly user): mount NAS
(`/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets`), place
`data/debug/` files on NAS, `ssh [REDACTED_HOST]` (password entry).
All sbatch/srun execution on compute nodes (never login node).

- [ ] Stage: `src/1_stage_data/1_stage_data.sh --ds_name _debug`.
- [ ] Dataset-specific preprocessing: Joanito `seqtec` step — decide skip vs run full
  (debug subset already carries `seqtec` from `1.3.1_create_debug_dataset.R`).
- [ ] Preprocess: `1_submit_hpc_array.sh --ds_name _debug`; validation: h5ad loads,
  `X_pca`/`X_pca_harmony` present, ~2500 cells, runtime < 30s.
- [ ] Chunks: `1_prepare_chunks.sh test _debug` (per-dataset union, 1 sample/chunk).
- [ ] Annotation: `2_submit_hpc_array.sh _debug`.
- [ ] Merge: `3.1_submit_merge.sh _debug`; validation: each view h5ad obs gains
  layer1–3 / scATOMIC cols (NA where absent), counts layer intact, annotated h5ad loads
  in R/Python.
- [ ] Cluster verify items: CombinedPBMC 64G baseline; preprocess 16G (GongSharma);
  annotation 2h/16G vs 5×2 retries; `aux/scGateDB.rds` committed-cache note in
  `2_submit_hpc_array.sh` comment.
- [ ] After debug passes: run one real dataset (e.g. Kim) before full rollout.

## Phase 3 — src/5_run_benchmark_methods [agent implements; HPC runs]

- [ ] **3.1 Python methods**: rename + convert `1.2_benchmark_methods_py.qmd` →
  `1.1.1_benchmark_methods_py.py` (strip quarto; CLI dataset/view/method; consume
  preprocessed obsm/layers/counts; output `.feather`); add `1_submit_hpc_array.sh`
  (per method+dataset array, `SLURM_PARTITION` override incl. GPU) + `1.1_run_worker.sh`
  (`PYTHON_BIN`); align exec-time/memory logging to R `exec_time()` format.
- [ ] **3.2 R methods**: run `benchmark_methods_r.R` + `benchmark_pipeline.R` on HPC via
  single-worker bash script; NAS targets `benchmark/{embeddings,plots}/`.
- [ ] **3.3 New methods**: PILOT-GM-VAE (add to py script + `constants.R` + R ingest;
  Harmony `X_pca_harmony` input for batch-effect views); QOT (feasibility test, deps
  from `QOT_PDAC_Example.ipynb`, no package); PULSAR (requirements test: UCE input,
  GPU/VRAM — may not be runnable).
- [ ] **3.4 Notebook adaptation**: `benchmark_analysis.rmd` + `batch_effect_analysis.rmd`
  read preprocessed h5ad (Step 4a approach: `ReadH5AD`/reticulate — benchmark on debug
  dataset), paths from datasets.json view outputs, ingest `.feather` from NAS; strip
  data-processing steps moved to HPC scripts.
- [ ] **3.5 SLURM config**: per-pipeline partition strategy (GPU for some methods;
  runtime check) in `slurm_config.sh`.
- [ ] **3.6 Docs**: README usage/workflow, ARCHITECTURE.md, AGENTS.md.
- [ ] Validation: `bash -n`/`py_compile`/R parse; debug-dataset run on HPC once implemented.

## Phase 4 — Batch effect analysis (later)

- Methods: ECODA batch-associated CT removal (t-test/Wilcoxon, ANOVA/Kruskal-Wallis);
  Pseudobulk DESeq2+limma with `batch_col`; MrVI native `batch_key`; GloScope on
  `X_pca_harmony`; PILOT-GM-VAE on `X_pca_harmony`; CombinedPBMC (Stephenson,
  GongSharma, Zhu) dataset handling; `columns.batch` in datasets.json (Joanito `seqtec`
  DONE via `1.2.1_create_joanito_batch_col.R`).

## Human-managed tasks (not agent)

- Batch effect datasets: whole Stephenson (n=143), KPMP Kidney (n=45), breast cancer
  (n=126), Covid-19 PBMC (n=151), diabetes (n=52), possibly Sikkema Lung.
- Benchmark datasets: Alzheimer (n=83), Lupus PBMC (n=261), myocardial infarction
  (n=23), possibly KPMP (n=45); GongSharma other subsetting conditions.
- Place debug dataset on NAS before Phase 2.

## Ideas for later

- GloScope on HPC; MOFAcellular; cell/sample/annotation counts from h5ad without full
  load.

## Keep-draft notes

- `preprocess_gongsharma.qmd` + `TODO_STUMP_preprocess_sikkema.qmd` are intentional
  drafts — do NOT delete (mark as such in TODO.md).

## Sections to remove from current TODO.md (all DONE/stale; git history preserves)

- "Open reviewer points" + "Suggested order / priority — DONE"
- "HPC pipeline review — Phase 1 fixes", "Round 2 Phase 1 fixes", "Completed (HPC layout…)"
- Phase-2 Item B (done); Steps 1, 2 (stale), 3 (done), 7 (folded); "Other major goals"
- "Risk Register", "New Methods to Be Added" (folded into Phase 3.3), "New Datasets to Be
  Added" (→ Human-managed tasks), "Completed" (→ one-line changelog pointer)

## Validation (no pipeline runs — per AGENTS.md)

- `bash -n` on touched `.sh`; `python -m py_compile` on touched `.py`; R parse on
  touched `.R` (incl. `1.3.1_create_debug_dataset.R`); `jq` smoke tests (staging filter
  emits `_debug` only when requested; default-all excludes `_*`).
- Greps: no references to removed TODO sections; no `*.h5ad` glob interference from
  `annotation_union/`; `3.1_submit_merge.sh` documented in ARCHITECTURE.md; draft qmds
  still present and marked.
- After Phase 1.3: `cells_lowres` creation exists for Kfoury (per user decision);
  `Preprocess_datasets.Rmd` deleted only after user confirmation.
- After Phase 1.4: `.Rprofile` consumers verified (reticulate); `.Renviron` decision
  documented (keep/update/delete).
