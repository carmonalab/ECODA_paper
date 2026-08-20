# ECODA_paper — Project Roadmap & Task Tracking

Implementation plan and task tracking for manuscript revisions, benchmark extensions, and pipeline work.

---

## Priority Overview

- **Priority 1:**
  - Multi-batch benchmark with batch-mixing metrics (4.2).
  - Ecotypes TNBC patient clustering (6.1).
  - MrVI vs. ECODA signal attribution (6.2).
  - Batch-correction impact on unsupervised cell-type annotation (4.5).
  - Clustering-resolution scan (6.5).
  - New dataset onboarding (Phase 5).
- **Priority 2:**
  - Downsampling robustness (6.4).
  - Marker-cell-type stability and Figure 3B heatmap (6.3).
  - Foundation models (PULSAR) feasibility (6.9).
- **Manuscript Text Only (Author Managed):**
  - Clinical circularity table per dataset, objective alignment (MrVI/scPoli vs. unsupervised clustering), and translational claim adjustments (6.12).

---

## Phase 3 — Benchmark Pipeline Rollout (`src/5_run_benchmark_methods/`)

- [x] **3.2 New Benchmark Methods Integration:**
  - Implemented PILOT-GM-VAE and QOT in `1.1.1_benchmark_methods_py.py` (`--methods qot,pilotgm`).
  - Vendored and hotfixed `qot_utils_re.py` (documented in [`docs/qot_hotfixes.md`](docs/qot_hotfixes.md)).
  - R ingest functions (`process_qot_fig`, `process_pilotgm_fig`) and labels wired in `benchmark_analysis.rmd`.
- [x] **3.3 Notebook Decoupling & Interactive Analysis:**
  - `benchmark_analysis.rmd` reads zero raw `.h5ad` files directly; ingests precomputed `.feather` matrices and `.rds` bundles via `load_hpc_benchmark_results()`.
  - Added Silhouette and LISI metrics into the core evaluation vector.
- [x] **3.6 HPC Composition Methods Array:**
  - Migrated notebook-local composition methods (ECODA family, GloProp, EPIC, Avg_PCA, Freq_highres) into an HPC array task (`composition` method).
  - Added peak RSS memory logging (`peak_rss_gb()`, VmHWM) to all R workers.
- [x] **3.4 Documentation Overhaul:**
  - Reorganized and decluttered `README.md`, `AGENTS.md`, and `docs/ARCHITECTURE.md` with strict separation of concerns.
- [ ] **3.1 Complete HPC Execution Across All Datasets [User Action]:**
  - Run full benchmark arrays on HPC for all registered cohorts.
  - Re-run stale PILOT bundles for Lee and Zhang (`1_submit_hpc_array.sh --methods pilot --force`).
  - Execute QOT and PILOT-GM-VAE across full benchmark cohorts.
- [ ] **3.5 SLURM Partition Config Cleanup:**
  - Review partition defaults in `src/slurm_config.sh`.

---

## Phase 4 — Batch Effect Analysis

- [x] **Pseudobulk DESeq2 + limma:** Batch-only correction implemented in `DESeq2.normalize()` / `get_pb_deseq2()`.
- [x] **CombinedPBMC Cohort:** Staging and combination of Stephenson + GongSharma + Zhu.
- [x] **Batch Metadata Tracking:** Integrated `columns.batch` in `datasets.json` (Joanito `seqtec`, Kfoury `cells_lowres`).
- [ ] **4.1 Generic ECODA Batch-Associated Cell Type Removal:**
  - Implement statistical test (t-test/Wilcoxon for 2 batches, ANOVA/Kruskal-Wallis for >2 batches, $p < 0.05$) to identify and optionally exclude batch-confounded cell types.
- [ ] **4.2 Multi-Batch Benchmark Suite (Priority 1):**
  - Evaluate datasets with defined technical batch structures (Stephenson, Joanito, CombinedPBMC, KPMP Kidney, Breast Cancer, Covid-19 PBMC).
  - Quantify biological separation vs. batch-mixing metrics (Silhouette, LISI) across all methods with and without batch correction.
- [ ] **4.3 MrVI Native `batch_key` Integration:**
  - Wire native `batch_key` parameter in batch-effect benchmark evaluations.
- [ ] **4.4 GloScope & PILOT-GM-VAE on Corrected Embeddings:**
  - Benchmark on `X_pca_harmony` embeddings.
- [ ] **4.5 Batch Correction Impact on Unsupervised Annotation:**
  - Test the effect of batch correction prior to unsupervised Leiden clustering and ECODA stratification.

---

## Phase 5 — New Dataset Onboarding (BIB Study Cohorts)

Reference plan: [`.agents/plans/dataset_onboarding_and_debug_overhaul.md`](.agents/plans/dataset_onboarding_and_debug_overhaul.md). Feasibility details: [`new_datasets_to_implement.md`](new_datasets_to_implement.md).

- [x] **5.1 Data Sources & Catalog:**
  - Identified Zenodo and CellxGene source locations for 9 target cohorts (Alzheimer, Breast Cancer, Covid-19 PBMC, Diabetes, Kidney KPMP, Lupus PBMC, Lung, Myocardial Infarction, Parkinson).
- [x] **5.2 HPC Download Pipeline (T1 & T1.1):**
  - Implemented `notebooks/dataset_onboarding/download_datasets_hpc.sh` and `run_download_worker.sh` (resumable `curl -C -`, Zenodo MD5 verification, CellxGene size checks, and automated NAS synchronization).
- [x] **5.3 Diagnostic Onboarding Tooling (T3.1):**
  - Scaffolded 9 per-dataset Quarto check notebooks (`dataset_check_<Name>.qmd`).
  - Implemented `onboarding_utils.py` (sample-first subsetting, count validation, crosstabs) and standalone `onboarding_metrics.R` (cell-level LISI on unintegrated PCA).
- [x] **5.4 Annotation Suitability Guardrails:**
  - Added `not_suitable_for_auto_annotation` flag handling in `benchmark_pipeline.R` and Figure 3 / Supp Fig 19.
- [ ] **5.5 Execute HPC Downloads [User Action]:**
  - Run `download_datasets_hpc.sh` on HPC to transfer datasets to scratch and sync to NAS.
- [ ] **5.6 Run Diagnostic Check Notebooks & Count Validation:**
  - Render onboarding notebooks, verify raw integer count matrices, and inspect unintegrated PCA/UMAP batch separation.
- [ ] **5.7 Dataset Review & Registration Checkpoint [User Decision]:**
  - Decide cohort categorization (benchmark vs. batch-effect vs. negative control vs. exclude).
  - Register approved cohorts in `datasets.json` (**ask user before editing**).
- [ ] **5.8 Diabetes Mouse-Gene Handling (Conditional):**
  - Add mouse gene ortholog mapping to `gene_utils.py` if Diabetes is approved for benchmark inclusion.
- [ ] **5.9 Pipeline Rollout across Approved Cohorts:**
  - Execute stage -> preprocess -> annotate -> benchmark array pipeline for new datasets.

---

## Phase 6 — Additional Reviewer Analyses

- [ ] **6.1 Ecotypes TNBC Patient Clustering (Priority 1):**
  - Unsupervised ECODA stratification on pre-treatment TNBC cohort ($n \approx 100$) evaluated against known chemotherapy response.
- [ ] **6.2 MrVI vs. ECODA Signal Attribution:**
  - Identify gene programs driving MrVI patient separation on Adams and Stephenson; score via UCell and compare with ECODA HVCs.
- [ ] **6.3 Marker-Cell-Type Stability & Figure 3B:**
  - Jaccard overlap of top marker genes across author annotations, HiTME, and multi-resolution Leiden clusters.
- [ ] **6.4 Downsampling Robustness (Priority 2):**
  - Measure separation score stability across varying cell counts per sample.
- [x] **6.5 Leiden Resolution Scan:**
  - Extended clustering resolution range up to resolution 50.
- [x] **6.6 Zero-Imputation Range Extension:**
  - Extended parameter evaluation range for zero-handling strategies.
- [ ] **6.7 Shuffled Baseline Runtime Comparison:**
  - Benchmark execution time table comparing ECODA and GloProp against randomized matrices.
- [ ] **6.8 LASSO-Penalized Classification Comparison:**
  - Compare unsupervised variance-based HVC selection against supervised sparse feature selection.
- [ ] **6.9 Foundation Models (PULSAR / UCE) Feasibility Check:**
  - Assess PBMC pretrained foundation model feasibility and requirements.
- [ ] **6.10 Demographics Supplementary Figure:**
  - Evaluate stratification consistency across age and sex stratifications.
- [x] **6.12 Author Response Text:**
  - Circularity discussion, objective alignment notes, and manuscript revisions handled by the author.

---

## Ideas & Technical Notes for Later

- **Python `.feather` Atomic Writes:** Python workers write `.feather` files directly; consider writing to `.tmp` followed by `os.replace` if interrupted tasks ever leave partial feathers.
- **PILOT-GM-VAE Runtime Documentation:** On large cohorts, PILOT-GM-VAE execution time is an intrinsic algorithmic property (several hours per combo on GPU); default combo (`hvg2000_highres`) serves as the primary benchmark comparison.
- **Watchdog Stale-RUNNING Recovery:** If SlurmDBD accounting lag exceeds the 20-minute poll cap and causes a false-fail status on a completed array, recovery is achieved via `1_submit_hpc_array.sh --sync-only <array_id>` (bypassing the watchdog's `STATE=FAIL` cache file).
