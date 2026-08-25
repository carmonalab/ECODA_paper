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
  - Implement statistical test (t-test/Wilcoxon for 2 batches, ANOVA/Kruskal-Wallis for >2 batches, $p < 0.05$) to identify and optionally exclude batch-confounded cell types. (will only be shown for information only but is not used for analysis as it is superseded by a more sophisticated method using lme4 as explained in NOTES.md).
- [ ] **4.2 Multi-Batch Benchmark Suite (Priority 1 — `batch_effect_analysis.rmd`):**
  - Implement two-pass benchmark evaluation:
    1. **Pass 1 (Uncorrected):** Run all benchmark methods on raw counts/embeddings. Evaluate metadata collinearity (NMI) and distance matrix PERMANOVA. Conduct cross-modality confounding checks (expression vs. composition) on partially confounded cohorts (e.g. `Site`/`Institution`). Finalize modality-appropriate batch variables based on uncorrected results.
    2. **Pass 2 (Corrected):** Run methods with modality-appropriate batch correction (`limma` for ECODA/pseudobulk, `Harmony` on PCA for PILOT/GloScope, native multi-batch for MrVI) adhering strictly to the No-Leakage Principle.
  - **Comparative Quantitative Metrics:**
    - PERMANOVA distance variance partitioning ($R^2_{\text{Bio}}$, $R^2_{\text{Batch}}$, $R^2_{\text{Shared}}$, $R^2_{\text{Residual}}$).
    - $\text{Bio} / \text{Batch}$ signal-to-noise ratio in uncorrected vs. corrected runs ($\text{Ratio}_{\text{Raw}}$ vs. $\text{Ratio}_{\text{Corr}}$).
    - Shift in biological retention ($\Delta R^2_{\text{Bio}}$) and batch suppression ($\Delta R^2_{\text{Batch}}$).
    - **Ratio-of-Ratios (Correction Benefit Index):** $\text{RoR} = \text{Ratio}_{\text{Corr}} / \text{Ratio}_{\text{Raw}}$.
  - **Visualization Decision [Open]:** Grouped PERMANOVA variance bar plot per method/covariate vs. FunkyHeatmap summary.
- [ ] **4.3 MrVI Native `batch_key` Integration:**
  - Wire native `batch_key` parameter in batch-effect benchmark evaluations.
- [ ] **4.4 GloScope & PILOT-GM-VAE on Corrected Embeddings:**
  - Benchmark on `X_pca_harmony` embeddings.
- [ ] **4.5 Batch Correction Impact on Unsupervised Annotation:**
  - Test the effect of batch correction prior to unsupervised Leiden clustering and ECODA stratification.

---

## Phase 5 — New Dataset Onboarding (BIB Study Cohorts)

Reference plan: [`.agents/plans/archive/implementation_plan_new_datasets_json-superseded-20260824.md`](.agents/plans/archive/implementation_plan_new_datasets_json-superseded-20260824.md). The checked-in authorities are `datasets.json`, `notebooks/dataset_onboarding/dataset_specs.py`, and regenerated audit JSON.

- [x] **5.1 Data sources and catalog:** Identify the nine Joodaki et al. cohorts and canonical source files.
- [x] **5.2 HPC download pipeline:** Keep resumable download, checksum, and NAS synchronization gates.
- [x] **5.3 Declared sample and annotation roles:** Record the user-confirmed sample, label, low/high tiers, annotation source, decision notes, and audit-warning policy in `dataset_specs.py`.
- [x] **5.4 Derived annotation gate:** Keep HiTME/Leiden roles at `PASS_PENDING_DERIVED_ANNOTATION` until processed h5ad evidence proves presence, nonzero coverage, cardinality ≥2, and the selected hierarchy.
- [x] **5.5 Lung registry subset:** Apply exact `platform == "10x"` filtering before audits; persist pre/post cell/sample counts, mask expression, assay-mask comparison, and the filtered `ann_coarse → ann_fine` hierarchy.
- [ ] **5.6 Regenerate full-file audits:** Run all nine cohorts on Bamboo, pull metadata, and require exact agreement with `datasets.json`; no stale heuristic choice may activate a registry entry.
- [x] **5.7 Two-pass registry contract:** Declare `batch_effect_uncorrected` and `batch_effect_corrected` with pass-qualified outputs. New cohorts remain `use_for_benchmark=false`, `use_for_batch_effect=true`, and `columns.batch=null`; corrected invocation fails closed.
- [x] **5.8 Pass-aware preprocessing and runners:** Route exact views/passes, isolate scratch/NAS artifacts, enforce the fixed method allow-list, and never reuse benchmark-named batch artifacts.
- [ ] **5.9 Uncorrected evidence checkpoint:** Run the ordered cohort suite, generate nine candidate CSVs plus `batch_candidate_review.csv`, and stop for one explicit technical-batch confirmation per new cohort.
- [ ] **5.10 Corrected execution after confirmation:** Validate the confirmed columns, run the paired corrected suite, verify pass-separated checksums/keys, update the mapping, and remove only the five authorized obsolete files.

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

- **Dataset Onboarding Batch Effect & Collinearity Synthesis:**
  - After re-running diagnostic reports with balanced subsets ($\min(N_{\text{samples}}, 25)$ samples, 1,000 cells/sample for new datasets; 500 cells/sample for `_debug`):
    - Document concise summary findings per dataset in `notebooks/dataset_onboarding/README.md` (identified batch effect variables, degree of collinearity with biological conditions, LISI within-cell-type separation).
    - Post-onboarding extension (full cohort): generate cross-metadata collinearity matrices (Cramér's V / normalized mutual information) and atlas-wide LISI heatmaps to comprehensively capture all technical vs biological variance sources.
    - Finalize "Batch condition sequencing/sample prep (ECODA)" metadata assignments in `datasets.json` / `README.md` once all full-cohort samples are processed and evaluated.
- **Python `.feather` Atomic Writes:** Python workers write `.feather` files directly; consider writing to `.tmp` followed by `os.replace` if interrupted tasks ever leave partial feathers.
- **PILOT-GM-VAE Runtime Documentation:** On large cohorts, PILOT-GM-VAE execution time is an intrinsic algorithmic property (several hours per combo on GPU); default combo (`hvg2000_highres`) serves as the primary benchmark comparison.
- **Watchdog Stale-RUNNING Recovery:** If SlurmDBD accounting lag exceeds the 20-minute poll cap and causes a false-fail status on a completed array, recovery is achieved via `1_submit_hpc_array.sh --sync-only <array_id>` (bypassing the watchdog's `STATE=FAIL` cache file).
