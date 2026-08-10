# ECODA_paper — TODO

Implementation plan for the remaining pipeline work. Phases are ordered:
**Phase 1** (agent, done) → **Phase 2** (debug run on HPC, requires user) →
**Phase 3** (benchmark methods, agent + HPC) → **Phase 4** (batch effect
analysis) → **Phase 5** (pipeline robustness fixes, agent) → human-managed
tasks. Completed history is preserved in git; see `git log` and the changelog
at the bottom.

## Phase 2 — Debug run on HPC [REQUIRES USER: connect NAS + log in to HPC]

Prereqs (explicitly user):
1. Connect NAS + log in to HPC:
   `ssh [REDACTED_HOST]` (password entry). All sbatch/srun
   execution on compute nodes (never login node).
2. [X] Stage the full datasets (incl. the Joanito raw .rds):
   `src/1_stage_data/1_stage_data.sh` on the login node (`_debug` is skipped —
   its raw subset never lives on the NAS; the Joanito step builds it).

- [ ] Dataset-specific preprocessing: `1_submit_hpc.sh` — the Joanito step
      (`1.2_submit_joanito.sh` → `1.2.1_prepare_joanito.R`) computes `seqtec`
      AND builds the `_debug` subset into `${HPC_SCRATCH_DIR}/_debug/data/`;
      verify the h5ad exists after the step. Run `1_submit_hpc.sh` only if the
      full datasets are staged. (Joanito preprocessing needs to be run again, see Phase 6)
- [X] Preprocess: `1_submit_hpc_array.sh --ds_name _debug`; validation: h5ad loads, `X_pca`/`X_pca_harmony` present, ~2500 cells, runtime < 30s.
- [X] Chunks: `1_prepare_chunks.sh test _debug` (per-dataset union, 1 sample/chunk).
- [X] Annotation: `2_submit_hpc_array.sh _debug`.
- [X] Merge: `3_submit_merge.sh _debug`; validation: each view h5ad obs gains layer1–3 / scATOMIC cols (NA where absent), counts layer intact, annotated h5ad loads in R/Python.
- [ ] Cluster verify items: CombinedPBMC parallelized (3 in-job fork workers, backed GongSharma read, `_intermediates/` per-source outputs; expect ~20-25 min first run vs ~40-60 min serial, faster on reruns via the `_raw.h5ad` cache; 128G/16 cpus sbatch); preprocess 16G (GongSharma); annotation 2h/16G vs 5×2 retries; `aux/scGateDB.rds` committed-cache note in `2_submit_hpc_array.sh` comment.
- [ ] After debug passes: run one real dataset (e.g. Kfoury) before full rollout. (partly done, src/5_run_benchmark_methods/run_python_sample_embedding_methods is waiting in the HPC queue)(started running the full pipeline on all datasets, see phase 6, still some debugging necessary, currenltly running src/3_scrnaseq_preprocessing and debugging)

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
  Combos run defaults-first (MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
  PILOT_hvg2000_highres) so the main methods' `mem_GB`/`time_secs` are measured
  before in-process memory bloat (peak RSS is monotonic within a process;
  non-default combos keep their relative order via stable sort).
- [ ] **3.2 R methods**: run `benchmark_methods_r.R` + `benchmark_pipeline.R` on HPC via
  a single-worker bash script; NAS targets `benchmark/{embeddings,plots}/`.
  CODE COMPLETE (R benchmark pipeline `run_r_sample_embedding_methods/` — GloScope,
  MOFA, Pseudobulk, scITD + `prepare_pseudobulk` prep step — plus the
  transformation/zero-imputation pipeline `run_transformation_zeroimp_analysis/`;
  `benchmark_pipeline.R` split into `run_benchmark_analysis` (fast,
  composition-based, notebook-side) + `load_hpc_benchmark_results` + HPC drivers;
  GloScope `_sqrtmat` variant merged; shared exec-log schema + NAS `benchmark/`
  target; see ARCHITECTURE.md). HPC debug validation PENDING — smoke test:
  Pipeline B first (`--ds_name _debug --analysis trans,zeroimp`), then Pipeline A
  `--ds_name _debug --methods prepare_pseudobulk,pseudobulk`; check
  `benchmark/results/`, `benchmark/pseudobulks/`, `execution_times.feather` on NAS.
- [ ] **3.3 New methods**: PILOT-GM-VAE (add to py script + `constants.R` + R ingest;
  Harmony `X_pca_harmony` input for batch-effect views); QOT (feasibility test, deps
  from `QOT_PDAC_Example.ipynb`, no package); PULSAR (requirements test: UCE input,
  GPU/VRAM — may not be runnable).
- [ ] **3.4 Notebook adaptation**: `benchmark_analysis.rmd` + `batch_effect_analysis.rmd`
  read preprocessed h5ad (Step 4a approach: `ReadH5AD`/reticulate — benchmark on debug
  dataset), paths from datasets.json view outputs, ingest `.feather` from NAS; strip
  data-processing steps moved to HPC scripts.
- [ ] **3.6 Docs**: README usage/workflow, ARCHITECTURE.md, AGENTS.md.
- [ ] **3.7 SLURM config cleanup**: resolve or drop the leftover
      `# TODO: Adapt for specific pipelines` comment on `SLURM_PARTITION`
      (`src/slurm_config.sh:114`) — decide per-stage partitions for stages 2–4
      or remove the comment.
- [ ] Validation: `bash -n`/`py_compile`/R parse; debug-dataset run on HPC once implemented.

## Phase 4 — Batch effect analysis (later)

- Methods: ECODA batch-associated CT removal (t-test/Wilcoxon, ANOVA/Kruskal-Wallis);
  Pseudobulk DESeq2+limma with `batch_col`; MrVI native `batch_key`; GloScope on
  `X_pca_harmony`; PILOT-GM-VAE on `X_pca_harmony`; CombinedPBMC (Stephenson,
  GongSharma, Zhu) dataset handling; `columns.batch` in datasets.json (Joanito `seqtec`
  DONE via `1.2.1_prepare_joanito.R`; Kfoury `cells_lowres` DONE via
  `1.3.1_create_kfoury_lowres_ct.R`).

## Phase 5 — Annotation completeness guard [agent]

- [x] `1_prepare_chunks.sh` / `1.1_prepare_chunks.py`: guard against building the
      annotation union from a PARTIALLY preprocessed dataset. The current skip
      predicate only checks for ≥1 `output/*.h5ad`; for multi-view datasets
      (Stephenson, `_debug`) a still-running preprocess task can
      have written only one of its views, so the union/chunks are built from an
      incomplete view set — and the dataset is then marked "already annotated"
      and skipped forever by later runs (silent incomplete annotations; only
      `--force` repairs). Done via (b): `1_prepare_chunks.sh` verifies every
      expected view from `datasets.json` exists (missing → WARNING-skip into a
      `SKIPPED_INCOMPLETE` summary bucket, exit 0, picked up on re-run;
      no-views datasets like Zhu keep the plain h5ad check), and
      `1.1_prepare_chunks.py` repeats the check fail-closed (CRITICAL exit →
      `FAILED_DATASETS`) as a bypass/drift guard.
- [x] Optional: `2_submit_hpc_array.sh` per-chunk skip (feather already present)
      so re-runs on annotated-but-not-yet-merged datasets don't redo work (no
      skip logic today; redundancy is only prevented by 1_prepare's skip + the
      merge deleting `output/chunks/`). Done in `2.1_run_worker.sh`:
      a chunk whose `annotations_chunk_<N>.feather` already exists exits 0
      (task COMPLETED in sacct); the chunk manifest stays unfiltered so the
      `3_submit_merge.sh` coverage gate is unchanged.

## Phase 6 — Preprocess array 4294824 recovery [agent code fixes + user HPC steps]

Context: 9/14 tasks COMPLETED, 5 failed; gate failed closed → nothing synced to
NAS. Full diagnosis + implementation detail: `.kilo/plans/1786319218000-preprocess-array-4294824-recovery.md`.

Agent code fixes:
- [X] `select_hvgs_ranked` (`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`):
      catch loess `ValueError` (Smillie/Zhang) and retry once with a
      deterministic jittered counts layer (seed 42, 1e-6 on CSR `.data`,
      `check_values=False`); delete the jitter layer afterwards; re-raise the
      original error if the retry fails. (Validated locally with pinned
      scanpy 1.12.2.)
- [X] `convert_rds_to_raw_h5ad` (`src/utils/py/preprocess_utils.py`): sanitize
      character meta.data columns to ASCII (`iconv(from="latin1", to="ASCII",
      sub=" ")`) before `write_h5ad` — anndataR writes `encoding='ascii'` even
      for non-ASCII bytes (Joanito `Stage.TNM` NBSP → UnicodeDecodeError).
- [X] `convert_rds_to_raw_h5ad` (same file): defensively reindex Assay5 layer
      rownames to the assay features before writing (Wu `validate_aligned_mapping`
      mismatch), with a diagnostic message naming the differing genes. Note:
      layers without rownames (fresh Assay5 objects, the healthy path) are
      skipped — only layers that actually carry rownames are aligned; a
      rowname set that is not the same gene set stops loudly.

HPC manual steps [REQUIRES USER — agents cannot run HPC]:
- [ ] Repair py-cuda13 R env (login node, installs allowed):
      `~/.pixi/bin/pixi run -e py-cuda13 Rscript -e 'install.packages("mime", repos = "https://cloud.r-project.org")'`,
      then `~/.pixi/bin/pixi run -e py-cuda13 setup` (idempotent). Fixes Zhu
      (missing `mime` lazyload DB); also protects annotation/benchmark workers.
- [X] Delete broken Joanito conversion cache:
      `rm "${HOME}/scratch/ECODA_paper/Joanito/output/JoaI_2022_35773407_Nofilt_whole_raw.h5ad"`
      (regenerates sanitized after the converter fix).
- [X] Check Stephenson task 11 (stale RUNNING at gate):
      `sacct -j 4294824_11 --format=JobID,State,Elapsed,End -X` +
      `ls "${HOME}/scratch/ECODA_paper/Stephenson/output/"`; re-run
      `--ds_name Stephenson` if outputs missing.
- [ ] Re-run failed datasets individually (each run syncs its own outputs):
      `--ds_name` Joanito, Smillie, Wu, Zhang, Stephenson.
- [ ] Sync the 9 COMPLETED-but-unsynced datasets by re-running their `--ds_name`
      (already-processed outputs are skipped, then synced): Adams, Bassez,
      CombinedPBMC, Kfoury, Kim, Lee, Pelka.
- [ ] GongSharma OOM fix (task 4, 128G): implement
      `src/2_dataset_specific_preprocessing/1.4_submit_gongsharma.sh` +
      `1.4.1_subset_gongsharma.py` (new step, auto-discovered by
      `1_submit_hpc.sh`): per-sample cap of 5000 cells
      (`specimen.specimenGuid`, `np.random.RandomState(42)`) — historical
      `downsample_by_group` strategy (git 3a4711e, `src/py/preprocess_gongsharma.qmd`).
      Decide: (a) in-place overwrite of the two staged SoundLife h5ads (no
      datasets.json change; re-run the step after any re-staging), or (b) write
      a single capped h5ad + update `datasets.json` (needs explicit approval —
      ground truth). Then re-run:
      `./src/2_dataset_specific_preprocessing/1_submit_hpc.sh` followed by
      `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name Gongsharma_cmv_young_males`.

## Human-managed tasks (not agent)

- Batch effect datasets: whole Stephenson (n=143), KPMP Kidney (n=45), breast cancer
  (n=126), Covid-19 PBMC (n=151), diabetes (n=52), possibly Sikkema Lung.
- Benchmark datasets: Alzheimer (n=83), Lupus PBMC (n=261), myocardial infarction
  (n=23), possibly KPMP (n=45); GongSharma other subsetting conditions.

## Ideas for later

- MOFAcellular; cell/sample/annotation counts from h5ad without full load.
- Gene blacklist before HVG selection: dump `aux/genes.blocklist.rds` (STACAS
  default_black_list) to a text file (one gene per line; `full` and `no_sex`
  variants), add `load_blacklist(path, exclude_sex=True)` to
  `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`, apply before HVG selection
  (`adata = adata[:, ~adata.var_names.isin(blacklist)].copy()`).
- Batch effect analysis: decide whether to run with and/or without batch
  correction — more important to only do WITH batch correction; non-corrected
  results possibly in the paper appendix.
- Phase 4 details: verify `DESeq2.normalize()` `batch_col` is correctly
  implemented and does not get `"Sample"` as batch column; ECODA batch-associated
  CT removal should print a warning naming the significant cell types
  (t-test/Wilcoxon for 2 batches, ANOVA/Kruskal-Wallis for >2, p < 0.05); test
  each cell type separately vs. checking global variance of cell type
  composition across batches.

## Keep-draft notes

- `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd` (GongSharma other-subsetting
  conditions) + `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd` (future
  Sikkema Lung dataset) are intentional drafts — do NOT delete.
