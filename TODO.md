# ECODA_paper — TODO

Implementation plan for the remaining pipeline work. Phases are ordered:
**Phase 1** (agent, done) → **Phase 2** (debug run on HPC, done) →
**Phase 3** (benchmark methods + full dataset rollout, agent + HPC) →
**Phase 4** (batch effect analysis) → human-managed tasks. Completed history
is preserved in git; see `git log`.

## Phase 3 — src/5_run_benchmark_methods [agent implements; HPC runs]

- [ ] **3.1 Run pipelines for all remaining datasets** [IN PROGRESS — USER
      RUNNING ON HPC]: (2026-08-11) `_debug` + Kfoury validated end-to-end
      (Stage 1 → Stage 5, incl. all benchmark methods); implementation +
      syntax/parse checks + HPC debug validation DONE (former 3.1/3.2).
      Remaining commands (preprocess rollout currently running; then):
      - `./src/4_cell_type_annotation/1_prepare_chunks.sh production`
      - `./src/4_cell_type_annotation/2_submit_hpc_array.sh`
      - `./src/4_cell_type_annotation/3_submit_merge.sh`
      - `./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh`
      - `./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh`
      - `./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh`
      Dataset coverage — benchmark: Adams, Bassez, Gongsharma_cmv_young_males,
      Kim, Lee, Pelka, Smillie, Stephenson (benchmark view), Wu, Zhang;
      batch-effect: Joanito, Stephenson (batch_effect view), CombinedPBMC.
      Zhu: no views (feeds only the CombinedPBMC stage-2 step 1.2) — stage 2
      only, confirm participation.
      GongSharma cap validation (cap log: 531,291 + 365,000 = 896,291 cells,
      max 5000 per sample) checked when its preprocess task reaches the NAS
      sync gate.
      After all datasets complete: verify NAS outputs (preprocessed h5ads +
      benchmark bundles), then resume 3.2 (new methods), 3.3 (notebook
      adaptation), 3.4 (docs), 3.5 (SLURM config cleanup), and Phase 4.
- [ ] **3.2 New methods**: PILOT-GM-VAE (add to py script + `constants.R` + R ingest;
  Harmony `X_pca_harmony` input for batch-effect views); QOT (feasibility test, deps
  from `QOT_PDAC_Example.ipynb`, no package); PULSAR (requirements test: UCE input,
  GPU/VRAM — may not be runnable).
- [x] **3.3 Notebook adaptation**: `benchmark_analysis.rmd` + `batch_effect_analysis.rmd`
  read preprocessed h5ad (Step 4a approach: `ReadH5AD`/reticulate — benchmark on debug
  dataset), paths from datasets.json view outputs, ingest `.feather` from NAS; strip
  data-processing steps moved to HPC scripts.
  - `benchmark_analysis.rmd` DONE (2026-08-12): backed h5ad obs-only reads,
    NAS benchmark/pseudobulks/embeddings paths, unified exec times (NAS
    feather + bundle-derived rows), RAM plot (Supp fig 14B), funky-heatmap
    refactor (`build_funky_heatmap` in `src/utils/plotting.R`) with the
    `benchmark_metrics` notebook parameter, zeroimp flattening fix +
    underscore method-key rename in `run_zeroimp_analysis` (breaking change:
    zeroimp bundles must be re-run with `--force`, see user to-dos).
    `batch_effect_analysis.rmd` still pending (Phase 4).
- [ ] **3.4 Docs**: README usage/workflow, ARCHITECTURE.md, AGENTS.md.
- [ ] **3.5 SLURM config cleanup**: resolve or drop the leftover
      `# TODO: Adapt for specific pipelines` comment on `SLURM_PARTITION`
      (`src/slurm_config.sh:114`) — decide per-stage partitions for stages 2–4
      or remove the comment.

## Phase 4 — Batch effect analysis (later)

- Methods: ECODA batch-associated CT removal (t-test/Wilcoxon, ANOVA/Kruskal-Wallis);
  Pseudobulk DESeq2+limma with `batch_col` (DONE for the batch-effect notebook:
  `DESeq2.normalize()`/`get_pb_deseq2()` take `batch_col`/`blind`/`correct_batch`,
  wired batch-only at all 4 `batch_effect_analysis.rmd` call sites; benchmark
  defaults unchanged); MrVI native `batch_key`; GloScope on
  `X_pca_harmony`; PILOT-GM-VAE on `X_pca_harmony`; CombinedPBMC (Stephenson,
  GongSharma, Zhu) dataset handling; `columns.batch` in datasets.json (Joanito `seqtec`
  DONE via `1.3.1_prepare_joanito.R`; Kfoury `cells_lowres` DONE via
  `1.4.1_create_kfoury_lowres_ct.R`).

## Human-managed tasks (not agent)

- Batch effect datasets: whole Stephenson (n=143), KPMP Kidney (n=45), breast cancer
  (n=126), Covid-19 PBMC (n=151), diabetes (n=52), possibly Sikkema Lung.
- Benchmark datasets: Alzheimer (n=83), Lupus PBMC (n=261), myocardial infarction
  (n=23), possibly KPMP (n=45); GongSharma other subsetting conditions.

## Ideas for later

- R-method peak RAM: backfill sacct `MaxRSS` into `execution_times.feather`
  (R rows currently have `mem_GB = NA`; the notebook's RAM figure
  `Supp_fig_14B_benchmark_mem_GB.pdf` shows python methods only).
- ECODA+Pseudobulk distance combos (`ECODA_PB_combo_*`): legacy, disabled in
  `run_benchmark_analysis`, kept commented-out in `benchmark_analysis.rmd` for
  internal testing only — NOT shown in publication figures.
- Optional: HPC workers dumping per-dataset `obs.rds` + PCA-mean matrix if the
  notebook should ever avoid h5ad reads entirely (currently not needed —
  obs/obsm backed reads are light).
- MOFAcellular; cell/sample/annotation counts from h5ad without full load.
- Gene blacklist before HVG selection: dump `aux/genes.blocklist.rds` (STACAS
  default_black_list) to a text file (one gene per line; `full` and `no_sex`
  variants), add `load_blacklist(path, exclude_sex=True)` to
  `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`, apply before HVG selection
  (`adata = adata[:, ~adata.var_names.isin(blacklist)].copy()`).
- Batch effect analysis: decide whether to run with and/or without batch
  correction — more important to only do WITH batch correction; non-corrected
  results possibly in the paper appendix.
- Phase 4 details: `DESeq2.normalize()` `batch_col` is now correctly implemented
  and wired (batch-only, never `"Sample"` as batch column; no-leakage — see
  AGENTS.md); ECODA batch-associated
  CT removal should print a warning naming the significant cell types
  (t-test/Wilcoxon for 2 batches, ANOVA/Kruskal-Wallis for >2, p < 0.05); test
  each cell type separately vs. checking global variance of cell type
  composition across batches.

## Keep-draft notes

- `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd` (GongSharma other-subsetting
  conditions) + `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd` (future
  Sikkema Lung dataset) are intentional drafts — do NOT delete.
