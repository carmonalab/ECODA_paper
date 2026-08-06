# ECODA_paper — TODO

Implementation plan for the remaining pipeline work. Phases are ordered:
**Phase 1** (agent, done) → **Phase 2** (debug run on HPC, requires user) →
**Phase 3** (benchmark methods, agent + HPC) → **Phase 4** (batch effect
analysis) → human-managed tasks. Completed history is preserved in git; see
`git log` and the changelog at the bottom.

## Phase 1 — Complete pipelines src/1–4 (DONE, agent)

Implemented: debug/test execution mode (`--ds_name`/`--force`/`--test`, `_debug`
5-sample subset built by the Joanito step `1.2.1_prepare_joanito.R` into
`${HPC_SCRATCH_DIR}/_debug/data/`), per-dataset annotation + merge-back (union
h5ad, `3_submit_merge.sh` coverage gate), legacy `Preprocess_datasets.Rmd`
audit (deleted), `.Rprofile`/`.Renviron`/`RETICULATE_PYTHON` fixes. Details are
preserved in git history (see previous TODO.md versions in `git log`).

## Phase 2 — Debug run on HPC [REQUIRES USER: connect NAS + log in to HPC]

Prereqs (explicitly user):
1. Connect NAS + log in to HPC:
   `ssh halterc@login1.bamboo.hpc.unige.ch` (password entry). All sbatch/srun
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
- [ ] Merge: `3_submit_merge.sh _debug`; validation: each view h5ad obs gains layer1–3 / scATOMIC cols (NA where absent), counts layer intact, annotated h5ad loads in R/Python.
- [ ] Cluster verify items: CombinedPBMC 64G baseline; preprocess 16G (GongSharma); annotation 2h/16G vs 5×2 retries; `aux/scGateDB.rds` committed-cache note in `2_submit_hpc_array.sh` comment.
- [ ] After debug passes: run one real dataset (e.g. Kim) before full rollout.

## Phase 3 — src/5_run_benchmark_methods [agent implements; HPC runs]

- [x] **3.1 Python methods**: rename + convert `1.2_benchmark_methods_py.qmd` →
  `1.1.1_benchmark_methods_py.py` (strip quarto; CLI dataset/view/method; consume
  preprocessed obsm/layers/counts; output `.feather`); add `1_submit_hpc_array.sh`
  (per method+dataset array, `SLURM_PARTITION` override incl. GPU) + `1.1_run_worker.sh`
  (`PYTHON_BIN`); align exec-time/memory logging to R `exec_time()` format.
  CODE COMPLETE (slurm_config.sh benchmark vars + ARCHITECTURE.md updated);
  HPC debug validation PENDING — smoke test:
  `./1_submit_hpc_array.sh --ds_name _debug --methods mrvi` (then scpoli, pilot);
  check `benchmark/embeddings/` feathers + `execution_times.feather`.
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
  `1.3.1_create_kfoury_lowres_ct.R`).

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

- Documentation cleanup (README / ARCHITECTURE / AGENTS / TODO): single-source-of-truth
  per file type — README (usage + reference data), ARCHITECTURE (pipeline details, the
  authority), AGENTS (short, pointer-heavy), TODO (task tracking). Fixed: stale rpy2
  claim in README, hardcoded cohort count (now points to datasets.json), broken
  `#cell-type-annotation-pipeline-*` anchor, empty "Batch Effect Analysis" section and
  duplicated gene-reference provenance in ARCHITECTURE; AGENTS pipeline overview +
  HPC layout compressed into pointers; stale changelog bullets dropped (git history
  preserves them); TODO Phase 1 DONE section compressed.
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
