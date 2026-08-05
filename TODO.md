# ECODA_paper — TODO

Implementation plan for the remaining pipeline work. Phases are ordered:
**Phase 1** (agent, done) → **Phase 2** (debug run on HPC, requires user) →
**Phase 3** (benchmark methods, agent + HPC) → **Phase 4** (batch effect
analysis) → human-managed tasks. Completed history is preserved in git; see
`git log` and the changelog at the bottom.

## Phase 1 — Complete pipelines src/1–4 (DONE, agent)

### 1.1 Debug/test execution mode
- Joanito dataset-specific step (renamed `1.2_submit_joanito.sh` →
  `1.2.1_prepare_joanito.R`): computes the `seqtec` batch column in place
  (single source of truth for the batch definition) and, from the SAME
  in-memory object, builds the `_debug` 5-sample subset into
  `${HPC_SCRATCH_DIR}/_debug/data/JoaI_2022_35773407_debug_5samples.h5ad` —
  5 samples covering both biological conditions (sample.origin), batches
  (seqtec) and sites, 500 cells/sample (seeded random, seed 321), minimal obs
  cols (incl. `seqtec`, `Site`, sample/patient cols); h5ad via anndataR:
  X=None + layers["counts"], handled by the preprocess X-promotion. No local
  script, no manual rsync, no NAS fallback — the raw subset only ever lives on
  HPC scratch (built by the step; never staged).
- datasets.json: `_debug` entry (batch col `Site`, label/sample cols
  following Joanito conventions, `use_for_benchmark` + `use_for_batch_effect`,
  both views → the debug h5ad). `_debug.folder_name` is `null` (CombinedPBMC
  precedent) so `1_stage_data.sh` skips staging the raw subset; `_debug`
  *outputs* (preprocessed/annotated h5ads) sync to NAS as usual via the
  preprocess-array and merge syncs. Convention: default-all script loops skip
  `_*` keys unless explicitly requested via `--ds_name`.
- `src/1_stage_data/1_stage_data.sh`: `--ds_name <DS>` filter; default-all skips
  `_*` keys; `folder_name: null` datasets (incl. `_debug`) are skipped.
- `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`: `--ds_name` single-dataset
  mode (1-task array at the dataset's position in the sorted jq key list, so
  `1.1_run_worker.sh`'s `sed -n ${SLURM_ARRAY_TASK_ID}p` mapping needs no change);
  `--force` recomputes existing outputs (forwarded via `FORCE_PREPROCESS`);
  default-all skips `_*`; single-dataset NAS sync syncs only the requested
  dataset; sacct fail-closed gate holds for 1-task arrays.
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`: `--force` flag bypassing the
  "Already processed" skip; additionally promotes `X=None` inputs with a `counts`
  layer to `X` (anndataR-written files, e.g. the debug h5ad) before any
  X-dependent step.

### 1.2 Per-dataset annotation + merge-back (replaces per-view)
- `src/4_cell_type_annotation/1.1_prepare_chunks.py`: one chunk set per DATASET —
  builds a union h5ad (concat all view h5ads, dedup `(sample, barcode)`) at
  `${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union/union.h5ad` (outside the synced
  output dir). Memory-lean: single view → chunked directly; a view equal to the
  union (e.g. Stephenson benchmark ⊂ batch-effect) → hardlink (copy fallback);
  partial-overlap views only → in-memory concat + dedup (warning). `--test` mode
  kept (1 sample/chunk).
- `src/4_cell_type_annotation/3_merge_annotations.py`: argparse CLI
  (`--h5ad-path` required, `--annot-dir` defaults to the h5ad's parent,
  `--output-path` defaults to in-place); idempotent re-merge (drops existing
  annotation columns before the join); atomic write (temp file + `os.replace`).
- `src/4_cell_type_annotation/3.1_submit_merge.sh` (new, per-dataset `<DS>`):
  coverage gate against the chunk manifest (fails on partial annotation arrays),
  NAS reachability check BEFORE any destructive step, merges feathers into each
  view h5ad via `srun` (64G baseline), deletes `annotation_union/` + stale
  `output/chunks/`, rsyncs annotated h5ads to `${NAS_TARGET_DIR}/${DS}/output/`.
- `2_submit_hpc_array.sh` / `2.1_run_worker.sh`: unchanged (chunk-file driven);
  chunk manifest + feather naming verified consistent (global counter across
  datasets, feather names derive from chunk file names).
- Docs: ARCHITECTURE.md (union + `3.1_submit_merge.sh` as wrapped stages; HPC
  layout incl. `annotation_union/`), AGENTS.md (repo-structure notes, `_debug`
  convention).

### 1.3 Legacy `Preprocess_datasets.Rmd` audit — DONE, file deleted
- `Preprocess_datasets.Rmd` (repo root, 905 lines) audited, ported, and
  **deleted** (git history preserves it). Mapping of its steps:
  - scGate models + ProjecTILs loading → ported (scGate DB cache
    `aux/scGateDB.rds` in `4_cell_type_annotation`; ProjecTILs ref maps via
    `HOME_REF_DIR` in `2.1.1_process_chunk.R`).
  - "Create low res cell types for Kfoury" → ported:
    `src/2_dataset_specific_preprocessing/1.4.1_create_kfoury_lowres_ct.R` +
    `1.4_submit_kfoury_lowres_ct.sh` (dataset-specific step via `1_submit_hpc.sh`,
    in-place idempotent saveRDS; collapses author `cells` → `cells_lowres`
    Tcells/NKcells/Bcells/MoMac/DCcells).
  - gene symbols → `src/gene_utils.py`; preprocess loop → `1.1.1_preprocess.py`;
    HiTME/scATOMIC whitelist → `3_merge_annotations.py` + `2.1.1_process_chunk.R`.
  - "Create clearcut age for GongSharma" → dropped (user decision): deferred to
    the kept draft `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd`
    (future other-subsetting conditions).
  - "Export datasets without author annotation" (Lee/Zhang for scPoli) → dropped
    as a pipeline step (user decision); the legacy context is noted in
    ARCHITECTURE.md ("Legacy pipeline notes") since it explains decisions in
    `notebooks/benchmark_analysis.rmd` and the paper.
- `Figure_workflow_schematic.Rmd` untouched (out of scope).

### 1.4 `.Rprofile` / `.Renviron` audit
- `.Rprofile`: KEPT — fixed a real bug: the file had NO trailing newline, so R's
  startup profile loader silently dropped the final line and the profile was
  never applied under `pixi run Rscript` (RETICULATE_PYTHON stayed unset).
  Verified after the fix: under `pixi run Rscript`, `Sys.which("python")` → pixi
  env python and `RETICULATE_PYTHON` gets set.
- `src/slurm_config.sh`: exports `RETICULATE_PYTHON="${PROJECT_ROOT}/.pixi/envs/default/bin/python"`
  — under `PIXI_RSCRIPT` (`--vanilla`) the `.Rprofile` is NOT read, and
  reticulate's own discovery was observed to pick a stray
  `~/.virtualenvs/r-reticulate` (local Mac artifact; anndata import failed). The
  export pins the pixi python deterministically on HPC workers.
- `.Renviron` (`R_MAX_VSIZE=200Gb`): KEPT + commented (user decision). Verified
  that `--vanilla` (--no-environ) means it does NOT apply to HPC pipeline
  workers — it only affects interactive/notebook R sessions.

## Phase 2 — Debug run on HPC [REQUIRES USER: connect NAS + log in to HPC]

Prereqs (explicitly user):
1. Connect NAS + log in to HPC:
   `ssh [REDACTED_HOST]` (password entry). All sbatch/srun
   execution on compute nodes (never login node).
2. Stage the full datasets (incl. the Joanito raw .rds):
   `src/1_stage_data/1_stage_data.sh` on the login node (`_debug` is skipped —
   its raw subset never lives on the NAS; the Joanito step builds it).

- [ ] Dataset-specific preprocessing: `1_submit_hpc.sh` — the Joanito step
      (`1.2_submit_joanito.sh` → `1.2.1_prepare_joanito.R`) computes `seqtec`
      AND builds the `_debug` subset into `${HPC_SCRATCH_DIR}/_debug/data/`;
      verify the h5ad exists after the step. Run `1_submit_hpc.sh` only if the
      full datasets are staged.
- [ ] Preprocess: `1_submit_hpc_array.sh --ds_name _debug`; validation: h5ad loads, `X_pca`/`X_pca_harmony` present, ~2500 cells, runtime < 30s.
- [ ] Chunks: `1_prepare_chunks.sh test _debug` (per-dataset union, 1 sample/chunk).
- [ ] Annotation: `2_submit_hpc_array.sh _debug`.
- [ ] Merge: `3.1_submit_merge.sh _debug`; validation: each view h5ad obs gains layer1–3 / scATOMIC cols (NA where absent), counts layer intact, annotated h5ad loads in R/Python.
- [ ] Cluster verify items: CombinedPBMC 64G baseline; preprocess 16G (GongSharma); annotation 2h/16G vs 5×2 retries; `aux/scGateDB.rds` committed-cache note in `2_submit_hpc_array.sh` comment.
- [ ] After debug passes: run one real dataset (e.g. Kim) before full rollout.

## Phase 3 — src/5_run_benchmark_methods [agent implements; HPC runs]

- [ ] **3.1 Python methods**: rename + convert `1.2_benchmark_methods_py.qmd` →
  `1.1.1_benchmark_methods_py.py` (strip quarto; CLI dataset/view/method; consume
  preprocessed obsm/layers/counts; output `.feather`); add `1_submit_hpc_array.sh`
  (per method+dataset array, `SLURM_PARTITION` override incl. GPU) + `1.1_run_worker.sh`
  (`PYTHON_BIN`); align exec-time/memory logging to R `exec_time()` format.
- [ ] **3.2 R methods**: run `benchmark_methods_r.R` + `benchmark_pipeline.R` on HPC via
  a single-worker bash script; NAS targets `benchmark/{embeddings,plots}/`.
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
  DONE via `1.2.1_prepare_joanito.R`; Kfoury `cells_lowres` DONE via
  `1.4.1_create_kfoury_lowres_ct.R`).

## Human-managed tasks (not agent)

- Batch effect datasets: whole Stephenson (n=143), KPMP Kidney (n=45), breast cancer
  (n=126), Covid-19 PBMC (n=151), diabetes (n=52), possibly Sikkema Lung.
- Benchmark datasets: Alzheimer (n=83), Lupus PBMC (n=261), myocardial infarction
  (n=23), possibly KPMP (n=45); GongSharma other subsetting conditions.

## Ideas for later

- GloScope on HPC; MOFAcellular; cell/sample/annotation counts from h5ad without full
  load.

## Keep-draft notes

- `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd` (GongSharma other-subsetting
  conditions) + `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd` (future
  Sikkema Lung dataset) are intentional drafts — do NOT delete.

## Changelog

- Phase 1 (debug mode, per-dataset annotation + merge, legacy audit, .Rprofile/.Renviron
  audit) implemented — details preserved in git history (see previous TODO.md versions
  in `git log`).
- `Preprocess_datasets.Rmd` deleted after the legacy audit (user-confirmed); debug
  dataset location decision: gitignored project-root `data/debug/` instead of NAS
  `Standardized_SingleCell_Datasets/debug/output/`.
- Debug subset workflow simplified (user-approved): the local
  `1.3.1_create_debug_dataset.R` + manual rsync / NAS fallback are replaced by the
  merged Joanito step `1.2_submit_joanito.sh` → `1.2.1_prepare_joanito.R`, which
  computes `seqtec` (single source of truth) and builds the `_debug` 5-sample subset
  into `${HPC_SCRATCH_DIR}/_debug/data/` from the same in-memory object;
  `_debug.folder_name` is now `null` (raw subset never on NAS; `_debug` outputs sync
  to NAS as usual). `data/debug/` removed from `.gitignore`.
