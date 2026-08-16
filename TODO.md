# ECODA_paper — TODO

Implementation plan for the remaining pipeline work. Phases are ordered:
Phase 3 (benchmark rollout) → Phase 4 (batch effect analysis) → Phase 5 (new
datasets) → Phase 6 (additional reviewer analyses) → human-managed tasks.
Completed history is preserved in git; see `git log`.
Merged 2026-08-13 with TODO_revision_plan.md (now deleted) and
new_datasets_to_implement.md (kept as appendix, linked from Phase 5).

## Priority overview (from the 2026-07-29 revision plan, merged)

- Prio 1: multi-batch benchmark with batch-mixing metrics
  (4.2); Ecotypes TNBC patient clustering (6.1); MrVI-vs-ECODA signal
  attribution (6.2); batch-correction impact on unsupervised annotation (4.5);
  clustering-resolution impact (6.5); Figure 3B marker-gene heatmap (6.3);
  zero-handling range extension (6.6); new datasets (Phase 5).
- Prio 2: downsampling robustness (6.4); PULSAR (3.2, later step) /
  MOFA cellular (3.2); batch-mixing quantification supp fig (4.2).
- Response text only (no code): circularity per-dataset table, MrVI/scPoli
  objective mismatch, discovery-of-subgroups scope, translational-claim
  tone-down (6.12).

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
      benchmark bundles), then resume 3.2 (new methods), 3.4 (docs), 3.5
      (SLURM config cleanup), and Phase 4.
- [x] **3.2 New methods** (DONE 2026-08-14, benchmark view; see
      `.kilo/plans/1786651957910-pilotgm-qot-benchmark-implementation.md`):
      - PILOT-GM-VAE (Prio 1) + QOT (Prio 2) implemented in
        `1.1.1_benchmark_methods_py.py` (`--method qot|pilotgm`), vendored
        `qot_utils_re.py` (PennShenLab/QOT @ 28cd529880c1, two hotfixes —
        traceability in `docs/qot_hotfixes.md`), R ingest
        (`process_qot_fig`/`process_pilotgm_fig`, keys `QOT_hvg{n}` /
        `PILOT-GM-VAE_hvg{n}`), `constants.R` labels, HPC submitter arrays
        (qot → CPU, pilotgm → GPU). NOTE: the current `PILOT` method key is
        the plain-EMD variant (`pl.tl.wasserstein_distance`, no GM-VAE) —
        `PILOT-GM-VAE` is a distinct key. Validation: `_debug` e2e on NAS
        (both methods) + R ingest smoke test done locally; HPC validation
        pending (see `docs/qot_hotfixes.md` / plan §6).
      - Batch-effect view (`X_pca_harmony` input for PILOT-GM-VAE) stays
        Phase 4 (`4.4` below).
      - **PILOT NaN fix (commit `c547613`, `fill_unknown_ct` in
        `1.1.1_benchmark_methods_py.py:326`)**: HiTME `layer2` annotations on
        Lee/Zhang contain NaN for unclassified cells, which silently produced
        all-zero PILOT EMD distances. The fix is committed but the **Lee/Zhang
        PILOT bundles on NAS are STALE (all-zero EMD) — re-run pending**:
        `run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Lee --methods pilot --force`
        and same for `Zhang` (parallel-friendly, one terminal each).
      - **QOT + PILOT-GM-VAE `_debug` HPC validation**: `_debug` e2e run
        submitted 2026-08-16; once it passes, roll both methods out to all
        benchmark datasets (`--methods qot,pilotgm` on the python submitter;
        qot → CPU, pilotgm → GPU arrays). The notebook already ingests their
        feathers (skip-with-message while absent).
      - PULSAR (+UCE): LATER step, NOT included — foundation-model scale (UCE
        1280-dim embeddings, multi-GB pretrained weights, PBMC-specific
        pretrained checkpoints), unclear GPU/VRAM, out-of-domain for the
        multi-tissue benchmark. Feasibility check needed; candidate for a
        dedicated plan. See also 6.9.
      - MOFA cellular (Prio 2, feasibility) still pending.
      - Method-annotation columns for figures: "originally designed for sample
        representation", "provides batch correction".
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
- [x] Metrics in main benchmark: Silhouette + LISI included in the metrics
      vector (`benchmark_analysis.rmd:2557`); Modularity with multiple KNN
      variants (sqrt(n)/3/6/9) in `src/utils/scoring_metrics.R`
      (`calc_modularity`). DONE.
- [x] **3.6 Composition methods moved to HPC** (DONE 2026-08-16; plan
      `.kilo/plans/1786711296469-move-composition-methods-to-hpc.md`): the
      notebook-local composition methods (ECODA_* family, GloProp,
      Freq_highres, Avg_PCA_embedding, ECODA_deconv) now run as one
      `composition` array task per dataset (obs-only worker: backed h5ad obs
      + hvg2000 obsm + precomputed hvg2000 pseudobulk; set.seed(123) per
      dataset; HiTME/scATOMIC combos guarded on ct-column presence). The
      notebook reads ZERO h5ad files: bundles +
      `<ds>_metadata.rds` (labels, cell/sample counts) + python-method
      feathers; `result_list.rds` checkpointing dropped (fresh load every
      knit). R workers log peak RSS (`peak_rss_gb()`, `/proc/self/status`
      VmHWM) into the exec-log + each bundle (`mem_GB`) — R rows no longer
      NA. Adams pseudobulk sample-name fix (hyphen/underscore mismatch, the
      notebook silhouette crash) in the same plan (Phase 0).
- [ ] **3.4 Docs**: README usage/workflow, ARCHITECTURE.md, AGENTS.md (also
      update the AGENTS.md reference to the TODO phases after this restructure).
- [ ] **3.5 SLURM config cleanup**: resolve or drop the leftover
      `# TODO: Adapt for specific pipelines` comment on `SLURM_PARTITION`
      (`src/slurm_config.sh:114`) — decide per-stage partitions for stages 2–4
      or remove the comment.

## Phase 4 — Batch effect analysis [agent + HPC]

- [x] Pseudobulk DESeq2+limma with batch-only correction (Joanito + Stephenson
      wired via `DESeq2.normalize()`/`get_pb_deseq2()` `batch_col`/`blind`/
      `correct_batch` at all 4 `batch_effect_analysis.rmd` call sites; benchmark
      defaults unchanged).
- [x] CombinedPBMC (Stephenson+GongSharma+Zhu) handling in
      `batch_effect_analysis.rmd` (analyzed at line 594+).
- [x] `columns.batch` in datasets.json: Joanito `seqtec`
      (`1.3.1_prepare_joanito.R`), Kfoury `cells_lowres`
      (`1.4.1_create_kfoury_lowres_ct.R`).
- [ ] **4.1 ECODA generic `batch` argument**: generalize the hardcoded
      Stephenson `ct_remove` subset (`batch_effect_analysis.rmd:244-248`) —
      per-CT t-test/Wilcoxon (2 batches) or ANOVA/Kruskal-Wallis (>2), p<0.05;
      print a warning naming significant cell types. No-leakage: batch input only.
- [ ] **4.2 Multi-batch benchmark (Prio 1)**: datasets = KPMP Kidney, whole
      Stephenson (by center, n=143), Joanito, CombinedPBMC, full GongSharma
      (3 covariates), Covid-19 PBMC, Breast cancer, Diabetes. Separation with
      AND without batch correction; when >1 major batch, process one
      batch/covariate at a time (all methods except scPoli). Metrics: bio
      separation + batch-mixing (silhouette, LISI; cLISI vs the other metrics
      correlation on one dataset) → Supp fig.
- [ ] **4.3 MrVI native batch_key** in the batch-effect notebook (MrVI
      currently runs only at `_lowres` in `1.2_benchmark_methods_py.qmd`).
- [ ] **4.4 GloScope + PILOT-GM-VAE on `X_pca_harmony`** in the batch-effect
      notebook.
- [ ] **4.5 Impact of batch correction on unsupervised cell type annotation for
      ECODA (Prio 1)**.
- [ ] **4.6 GongSharma**: run Harmony on all samples (batch = sample), Leiden
      clustering in corrected space, compare labels uncorrected vs corrected,
      add to Figure 3.

## Phase 5 — New datasets (from the BIB benchmark study) [review task first]

Reviewer asked to add datasets from the study
https://academic.oup.com/bib/article/26/5/bbaf547/8287234 (BIB 2025). Prefer
downloading from where that study's authors provide the data, with no or minor
changes; whether a dataset shows batch effects decides benchmark vs batch-effect
usage in this repo. Full extraction (counts, classes, feasibility colors,
comments) is in `new_datasets_to_implement.md` (Excel source:
`/Users/christianhalter/Desktop/ECODA_PAPER_DATASETS.xlsx`).
- [ ] **5.1 Check the BIB study (bbaf547) for author-provided data locations**;
      verify availability (cellxgene/GEO per the table) and the amount of
      changes needed. [separate review task — not started]
- [ ] **5.2 Benchmark candidates (green)**: Alzheimer (n=83), Lupus PBMC
      (n=261), Myocardial infarction (n=23; has only clustering, no
      high-granularity CTs).
- [ ] **5.3 Batch-effect candidates (green/yellow)**: whole Stephenson by
      center (n=143; extension of the existing batch-effect view), KPMP Kidney
      (n=45; check batch effects first), Breast cancer (n=126; no cancer cells
      — contralateral unaffected samples), Covid-19 PBMC (n=151; GEO), Diabetes
      (n=52; yellow — 9 sub-datasets, many conditions).
- [ ] **5.4 Check-before-use (orange)**: Lung (n=165; too many technical
      conditions), Parkinson (n=97; composition dominated by brain region —
      possible negative control), Kidney cancer (n=17; too few samples),
      Pancreas PDAC (n=35; separates by "Ductal 2" cells only).
- [ ] **5.5 Sikkema Lung HCA** (reviewer-mentioned, data readily available;
      draft `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd`
      exists) — test sample separation by tissue.
- [ ] **5.6 Register approved datasets in datasets.json (ASK THE USER FIRST —
      datasets.json must not be changed without asking)**, with the appropriate
      view(s); then stage → preprocess → annotate → benchmark/batch-effect arrays.

## Phase 6 — Additional reviewer analyses [agent, Prio 1 first]

- [ ] 6.1 Ecotypes: fully unsupervised ECODA patient clustering on a large
      pre-treatment TNBC cohort (n≈100) with known chemotherapy response;
      TME-type clusters vs clinical response. Addresses Rev#1 circularity
      (also see 6.12). (skip this for now (but keep it here as open point). User will come back to this. it's crucial but user will let you know when and how to implement this)
- [ ] 6.2 MrVI vs ECODA (Adams + Stephenson): identify the main gene program
      driving MrVI separation that ECODA missed; PCA/UMAP colored by
      gene-program UCell scores. (scPoli batch-effect observation — adjR2
      quote — is discussion context for 6.2/6.12.)
- [ ] 6.3 Marker-cell-type stability: per-CT DE genes (padj/FC) or top-10-50
      markers; clustered Jaccard-overlap heatmap authors HR vs HiTME, then
      authors HR vs Leiden_res_5. Plus Figure 3B shared-marker dotplot/heatmap. (prio 2)
- [ ] 6.4 Downsampling sensitivity (Prio 2): separation score vs # cells per
      cell type/cluster on two datasets; release all min.cell thresholds. (prio 2, skip for now, user will come back to this)
- [ ] 6.5 Leiden resolution scan: extend the existing ECODA_seuratres_* range
      (currently 0.1-20) until separation drops; no min.cell filtering. (already done. extended range to max 50) (just leave this point here)
- [X] 6.6 Zero-handling range extension (current range too narrow). (already implented previously by user. range already extended)
- [ ] 6.7 Runtime table: shuffled vs ECODA/GloProp comparison. (very low prio)
- [ ] 6.8 LASSO-penalized classification comparison (does variance-based
      selection approximate supervised selection?). (very low prio)
- [ ] 6.9 Foundation models (PULSAR) + large-scale/multi-study scenarios
      (OneK1K, HLCA) — discuss, benchmark where runnable. (prio 2-3)
- [ ] 6.10 Supp fig: separation on females<40 and males>40 (same pattern as
      males<40 — no cherry-picking). [needs user guidance on cohorts] (prio 2 but very easy to implement. just tap h5ad files on nas (already there), just eeds to get cell type composition. minimal computational cost. just need to add supplementary figure, showing the different combinations, all same results, whole dataset already shown in figure 1 (i think))
- [ ] 6.11 Check how expert annotations were generated (manual vs automated
      classifiers). [needs user knowledge]
- [ ] 6.12 Response text (no code): circularity table per dataset (condition
      defined by clinical tests, not cell composition — Adams: HRCT/PFTs, etc.);
      MrVI/scPoli different objectives than unsupervised patient clustering;
      discovery of subgroups is out of scope (validated in 6.1); tone down
      translational marker-claim (multi-parameter panels not routine in
      clinical practice). (user did that, will be done by the user)

## Human-managed tasks (not agent)

- HPC rollout of Phase 5 datasets (stage/preprocess/annotate/benchmark arrays),
  guided by the BIB-study check (5.1).
- GongSharma other-subsetting conditions (draft `preprocess_gongsharma.qmd`).
- Where user input is required: 6.10 (cohort definition), 6.11 (annotation
  provenance), 6.12 (response text).

## Ideas for later

- CPU benchmark array throttling (`BENCHMARK_CPU_ARRAY_THROTTLE`, e.g. 4) for
  the R/transzeroimp/PILOT submitters — deferred from the 2026-08-13
  direct-env imports work: if the slim imports still show startup metadata
  storms on the CPU arrays, throttle concurrency as a separate design
  (`MAX_NUM_CHUNKS_PARALLEL` stays untouched). Long-term alternative:
  node-shared `/srv/share/users/...` staging (documented, not implemented).
- Optional (backfill): sacct `MaxRSS` backfill of `execution_times.feather`
  `mem_GB` for R-method bundles computed BEFORE 2026-08-16 (workers only
  started logging peak RSS then; legacy bundles keep `mem_GB = NA` until
  their next `--force` run).
- ECODA+Pseudobulk distance combos (`ECODA_PB_combo_*`): legacy, disabled in
  `run_benchmark_analysis`, kept commented-out in `benchmark_analysis.rmd` for
  internal testing only — NOT shown in publication figures.
- Cell/sample/annotation counts from h5ad without full load (MOFAcellular
  moved to 3.2; the benchmark notebook now gets them from
  `<ds>_metadata.rds` — the metadata-bundle pattern could be extended to the
  batch-effect notebook).
- Gene blacklist before HVG selection: dump `aux/genes.blocklist.rds` (STACAS
  default_black_list) to a text file (one gene per line; `full` and `no_sex`
  variants), add `load_blacklist(path, exclude_sex=True)` to
  `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`, apply before HVG selection
  (`adata = adata[:, ~adata.var_names.isin(blacklist)].copy()`).
- Batch effect analysis: decide whether to run with and/or without batch
  correction — more important to only do WITH batch correction; non-corrected
  results possibly in the paper appendix (4.2 runs both).
- Phase 4 details: `DESeq2.normalize()` `batch_col` is now correctly implemented
  and wired (batch-only, never `"Sample"` as batch column; no-leakage — see
  AGENTS.md); the ECODA batch-associated CT-removal warning detail (test each
  cell type separately vs. checking global variance, p < 0.05) moved to 4.1.
- on CUDA, torch's non-deterministic kernels (atomic reductions, cudnn autotuning) can still produce tiny run-to-run differences even with the same seed. Full GPU determinism would require torch.use_deterministic_algorithms(True) + CUBLAS_WORKSPACE_CONFIG, which scvi-tools doesn't enable.

## Keep-draft notes

- `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd` (GongSharma other-subsetting
  conditions) + `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd` (future
  Sikkema Lung dataset) are intentional drafts — do NOT delete.
