# Plan: Integrate TODO_revision_plan.md + new_datasets_to_implement.md into TODO.md

## Goal

Merge `TODO_revision_plan.md` (reviewer-driven priority tasks) and
`new_datasets_to_implement.md` (new-dataset info from the BIB benchmark study)
into the single source of truth `TODO.md`, marking what is already implemented,
preserving the revision plan's prioritization, and adding the new-datasets
workstream (incl. a check against the BIB study from which the datasets come).

Decisions already made with the user:
- After the merge, **delete `TODO_revision_plan.md`** (fully absorbed).
- **Keep `new_datasets_to_implement.md`** as an appendix (it holds the full
  Excel dump) and link it from the new datasets section in TODO.md.

## Repo facts (checked 2026-08-13)

- `TODO.md` currently: Phase 3 (benchmark, agent+HPC) / Phase 4 (batch effect,
  later) / Human-managed tasks / Ideas for later / Keep-draft notes. Phase 3.1
  is IN PROGRESS (user running HPC rollout).
- `AGENTS.md` line "Pipeline/code work ... see TODO.md (Phase 3 + Phase 4)" —
  must be updated after restructuring.
- `datasets.json` must not be changed without asking (AGENTS.md rule) — the new
  dataset registration step must carry that rule.
- `TODO_revision_plan.md` and `new_datasets_to_implement.md` are untracked
  (git status shows `??`).

## Implementation status of revision-plan items (evidence)

### Already implemented — do NOT re-add as open items; reference as done
- **Silhouette + LISI in main benchmark**: `notebooks/benchmark_analysis.rmd:2557`
  (`metrics <- c("ANOSIM", "Modularity", "ARI", "Silhouette", "LISI")`).
- **Leiden clustering resolution benchmark** (`ECODA_seuratres_*` 0.1/0.4/2/5/20):
  `benchmark_analysis.rmd:1033-1037` (extended figure; the *broader scan until
  separation drops* is still open → 6.5).
- **Modularity with multiple KNN variants**: `src/utils/scoring_metrics.R`
  (`calc_modularity`, KNN sqrt(n)/3/6/9).
- **Zero-imputation analysis pipeline** + real-dataset plots:
  `benchmark_analysis.rmd:3048+` (the *range too narrow* extension is open → 6.6).
- **Pseudobulk DESeq2 + limma batch-only correction** (`batch_col`/`blind`/
  `correct_batch`): `src/utils/pseudobulk.R:7-93`, wired at all
  `batch_effect_analysis.rmd` call sites.
- **ECODA batch-associated CT removal — PARTIAL**: hardcoded `ct_remove`
  subset for Stephenson only (`batch_effect_analysis.rmd:244-248`). The generic
  ECODA `batch` argument + per-CT statistical test + warning is NOT done → 4.x.
- **CombinedPBMC dataset** (Stephenson+GongSharma+Zhu): registered in
  `datasets.json`, analyzed in `batch_effect_analysis.rmd:594+`.
- **MrVI, scPoli, GloScope, PILOT in main benchmark**:
  `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd`
  (PILOT is the plain EMD variant, lines 219-239 — GM-VAE NOT done → 3.2).
- **`columns.batch`**: Joanito `seqtec` (`1.3.1_prepare_joanito.R`), Kfoury
  `cells_lowres` (`1.4.1_create_kfoury_lowres_ct.R`).
- **batch_effect_analysis.rmd partial**: Joanito (ECODA+PB), Stephenson
  (ECODA+PB+CT-removal subset), CombinedPBMC (combined pseudobulk). **MrVI and
  GloScope on harmony NOT yet in this notebook** → Phase 4.

### Not implemented (goes into the new TODO.md)
PILOT-GM-VAE; QOT; PULSAR; MOFA cellular; new datasets (Alzheimer, Lupus,
MI, KPMP, Breast cancer, Covid-19 PBMC, Diabetes, Lung, Sikkema, whole
Stephenson by center); Ecotypes/TNBC unsupervised patient clustering;
MrVI-vs-ECODA gene-program attribution; marker-gene stability/Jaccard heatmap
(Fig 3B); downsampling sensitivity; extended resolution scan; zero-handling
range extension; shuffled-vs-GloProp runtime table; LASSO comparison;
foundation-model (PULSAR) + large-scale (OneK1K/HLCA) discussion; batch-mixing
quantification supp fig; cLISI-vs-other-metrics correlation; impact of batch
correction on unsupervised annotation; full GongSharma (3 covariates);
GongSharma harmony-corrected Leiden (Figure 3); females<40/males>40 supp fig;
expert-annotation generation check; circularity/response text items.

## Target TODO.md structure (drafted content)

Replace `TODO.md` with the following. Keep existing Phase 1/2 history note
unchanged. Headers/priorities below are the merge result; adapt wording as
needed while preserving every item from the mapping table (appendix).

```markdown
# ECODA_paper — TODO

Implementation plan for the remaining pipeline work. Phases are ordered:
Phase 3 (benchmark rollout) → Phase 4 (batch effect analysis) → Phase 5 (new
datasets) → Phase 6 (additional reviewer analyses) → human-managed tasks.
Completed history is preserved in git; see `git log`.
Merged 2026-08-13 with TODO_revision_plan.md (now deleted) and
new_datasets_to_implement.md (kept as appendix, linked from Phase 5).

## Priority overview (from the 2026-07-29 revision plan, merged)
- Prio 1: PILOT-GM-VAE (3.2); multi-batch benchmark with batch-mixing metrics
  (4.2); Ecotypes TNBC patient clustering (6.1); MrVI-vs-ECODA signal
  attribution (6.2); batch-correction impact on unsupervised annotation (4.5);
  clustering-resolution impact (6.5); Figure 3B marker-gene heatmap (6.3);
  zero-handling range extension (6.6); new datasets (Phase 5).
- Prio 2: downsampling robustness (6.4); PULSAR/QOT/MOFA cellular (3.2);
  batch-mixing quantification supp fig (4.2).
- Response text only (no code): circularity per-dataset table, MrVI/scPoli
  objective mismatch, discovery-of-subgroups scope, translational-claim
  tone-down (6.12).

## Phase 3 — src/5_run_benchmark_methods [agent implements; HPC runs]
- [ ] 3.1 (unchanged, IN PROGRESS — user running on HPC; see current text)
- [ ] 3.2 New methods:
      - PILOT-GM-VAE (Prio 1): add to `1.2_benchmark_methods_py.qmd` +
        `constants.R` + R ingest; `X_pca_harmony` input for batch-effect views.
      - QOT (Prio 2; feasibility, deps from `QOT_PDAC_Example.ipynb`, no package).
      - PULSAR (Prio 2; requirements test: UCE input, GPU/VRAM — may not be runnable).
      - MOFA cellular (Prio 2, feasibility).
      - Method-annotation columns for figures: "originally designed for sample
        representation", "provides batch correction".
- [x] 3.3 Notebook adaptation — benchmark_analysis.rmd DONE (2026-08-12);
      batch_effect_analysis.rmd → Phase 4.
- [ ] 3.4 Docs (README/ARCHITECTURE/AGENTS) — also update the AGENTS.md
      reference to TODO phases after this restructure.
- [ ] 3.5 SLURM config cleanup (unchanged).

## Phase 4 — Batch effect analysis [agent + HPC]
- [x] Pseudobulk DESeq2+limma with batch-only correction (Joanito + Stephenson
      wired; benchmark defaults unchanged).
- [x] CombinedPBMC (Stephenson+GongSharma+Zhu) handling in batch_effect_analysis.rmd.
- [x] columns.batch in datasets.json: Joanito `seqtec`, Kfoury `cells_lowres`.
- [ ] 4.1 ECODA generic `batch` argument: generalize the hardcoded Stephenson
      ct_remove subset (batch_effect_analysis.rmd:244-248) — per-CT
      t-test/Wilcoxon (2 batches) or ANOVA/Kruskal-Wallis (>2), p<0.05; print
      a warning naming significant cell types. No-leakage: batch input only.
- [ ] 4.2 Multi-batch benchmark (Prio 1): datasets = KPMP Kidney, whole
      Stephenson (by center, n=143), Joanito, CombinedPBMC, full GongSharma
      (3 covariates), Covid-19 PBMC, Breast cancer, Diabetes. Separation with
      AND without batch correction; when >1 major batch, process one
      batch/covariate at a time (all methods except scPoli). Metrics: bio
      separation + batch-mixing (silhouette, LISI; cLISI vs the other metrics
      correlation on one dataset) → Supp fig.
- [ ] 4.3 MrVI native batch_key in the batch-effect notebook.
- [ ] 4.4 GloScope + PILOT-GM-VAE on `X_pca_harmony` in the batch-effect notebook.
- [ ] 4.5 Impact of batch correction on unsupervised cell type annotation for
      ECODA (Prio 1).
- [ ] 4.6 GongSharma: run Harmony on all samples (batch = sample), Leiden
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
- [ ] 5.1 Check the BIB study (bbaf547) for author-provided data locations;
      verify availability (cellxgene/GEO per the table) and the amount of
      changes needed. [separate review task — not started]
- [ ] 5.2 Benchmark candidates (green): Alzheimer (n=83), Lupus PBMC (n=261),
      Myocardial infarction (n=23; has only clustering, no high-granularity CTs).
- [ ] 5.3 Batch-effect candidates (green/yellow): whole Stephenson by center
      (n=143; extension of the existing batch-effect view), KPMP Kidney (n=45;
      check batch effects first), Breast cancer (n=126; no cancer cells —
      contralateral unaffected samples), Covid-19 PBMC (n=151; GEO), Diabetes
      (n=52; yellow — 9 sub-datasets, many conditions).
- [ ] 5.4 Check-before-use (orange): Lung (n=165; too many technical
      conditions), Parkinson (n=97; composition dominated by brain region —
      possible negative control), Kidney cancer (n=17; too few samples),
      Pancreas PDAC (n=35; separates by "Ductal 2" cells only).
- [ ] 5.5 Sikkema Lung HCA (reviewer-mentioned, data readily available; draft
      `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd` exists) —
      test sample separation by tissue.
- [ ] 5.6 Register approved datasets in datasets.json (ASK THE USER FIRST —
      datasets.json must not be changed without asking), with the appropriate
      view(s); then stage → preprocess → annotate → benchmark/batch-effect arrays.

## Phase 6 — Additional reviewer analyses [agent, Prio 1 first]
- [ ] 6.1 Ecotypes: fully unsupervised ECODA patient clustering on a large
      pre-treatment TNBC cohort (n≈100) with known chemotherapy response;
      TME-type clusters vs clinical response. Addresses Rev#1 circularity
      (also see 6.12).
- [ ] 6.2 MrVI vs ECODA (Adams + Stephenson): identify the main gene program
      driving MrVI separation that ECODA missed; PCA/UMAP colored by
      gene-program UCell scores.
- [ ] 6.3 Marker-cell-type stability: per-CT DE genes (padj/FC) or top-10-50
      markers; clustered Jaccard-overlap heatmap authors HR vs HiTME, then
      authors HR vs Leiden_res_5. Plus Figure 3B shared-marker dotplot/heatmap.
- [ ] 6.4 Downsampling sensitivity (Prio 2): separation score vs # cells per
      cell type/cluster on two datasets; release all min.cell thresholds.
- [ ] 6.5 Leiden resolution scan: extend the existing ECODA_seuratres_* range
      (currently 0.1-20) until separation drops; no min.cell filtering.
- [ ] 6.6 Zero-handling range extension (current range too narrow).
- [ ] 6.7 Runtime table: shuffled vs ECODA/GloProp comparison.
- [ ] 6.8 LASSO-penalized classification comparison (does variance-based
      selection approximate supervised selection?).
- [ ] 6.9 Foundation models (PULSAR) + large-scale/multi-study scenarios
      (OneK1K, HLCA) — discuss, benchmark where runnable.
- [ ] 6.10 Supp fig: separation on females<40 and males>40 (same pattern as
      males<40 — no cherry-picking). [needs user guidance on cohorts]
- [ ] 6.11 Check how expert annotations were generated (manual vs automated
      classifiers). [needs user knowledge]
- [ ] 6.12 Response text (no code): circularity table per dataset (condition
      defined by clinical tests, not cell composition — Adams: HRCT/PFTs, etc.);
      MrVI/scPoli different objectives than unsupervised patient clustering;
      discovery of subgroups is out of scope (validated in 6.1); tone down
      translational marker-claim (multi-parameter panels not routine in
      clinical practice).

## Human-managed tasks (not agent)
- HPC rollout of Phase 5 datasets (stage/preprocess/annotate/benchmark arrays),
  guided by the BIB-study check (5.1).
- Where user input is required: 6.10 (cohort definition), 6.11 (annotation
  provenance), 6.12 (response text).

## Ideas for later (unchanged)
- (keep existing bullets: CPU array throttle, R-method peak RAM backfill,
  ECODA_PB_combo legacy, obs.rds dump, MOFAcellular → now 3.2, gene blacklist,
  with/without batch correction appendix decision, Phase 4 details notes —
  the ECODA CT-removal warning detail moves to 4.1)

## Keep-draft notes (unchanged)
- preprocess_gongsharma.qmd + TODO_STUMP_preprocess_sikkema.qmd drafts — do NOT delete.
```

## Tasks for the implementing agent

1. Rewrite `TODO.md` with the structure above, folding in the existing
   Phase 3.1/3.3/3.4/3.5, Ideas-for-later and Keep-draft content verbatim where
   marked "unchanged" (read the current file first).
2. Verify no revision-plan item is lost: check each bullet of
   `TODO_revision_plan.md` against the mapping table below.
3. Update the AGENTS.md reference line to the new phase list
   ("see TODO.md (Phase 3 + Phase 4)" → "Phase 3-6" or per final headers).
4. Delete `TODO_revision_plan.md`; keep `new_datasets_to_implement.md` (linked
   from Phase 5). Do NOT touch `datasets.json`.
5. Commit (follow repo Task Completion Workflow: move this plan to
   `.kilo/plans/archive/`, `git add .`, commit, push).

## Validation

- Every bullet of `TODO_revision_plan.md` appears in the new TODO.md (mapping
  table below) — including the "1st/2nd/Prio 1" priority markers.
- Implemented items are marked [x] with the evidence references above; no
  done item is re-opened.
- `new_datasets_to_implement.md` linked from Phase 5, incl. the BIB-study
  download/check instruction and datasets.json ask-first rule.
- `git status` clean after commit; AGENTS.md reference consistent.

## Mapping table (revision-plan item → TODO.md destination)

| Revision-plan item | Destination |
| --- | --- |
| New benchmark, datasets with 1 batch/covariate (KPMP, whole Stephenson, COVID, Joanito, full G&S); sep. ± batch correction; one batch at a time, all methods except scPoli | 4.2 |
| Characterize newly discovered patient subgroups (breast cancer cohort) | 6.1 |
| Identify signals driving MrVI/PB separation vs ECODA | 6.2 |
| Impact of batch correction on unsupervised annotation | 4.5 |
| Impact of clustering resolution (high side) | 6.5 |
| Figure 3B shared-marker heatmap/dotplot | 6.3 |
| Range of zero handling | 6.6 |
| Add PILOT-GM-VAE | 3.2 |
| Downsampling robustness (release min.cell thresholds) | 6.4 |
| PULSAR, QOT, MOFA cellular | 3.2 |
| Supp fig: batch mixing quantification (silhouette, LISI) | 4.2 |
| Stratification in presence of batch effect: ECODA batch arg + warning | 4.1 |
| PB DESeq2 batch correction | [x] Phase 4 header |
| MrVI native batch correction | 4.3 |
| Harmony-corrected GloScope, PILOT-GM-VAE | 4.4 |
| Main benchmark: +Alzheimer, Lupus, MI (+KPMP, Lung check) | 5.2, 5.4 |
| Method columns: sample-representation design / batch correction | 3.2 |
| Fully unsupervised ECODA clustering (breast cancer n=100) | 6.1 |
| Circularity table per dataset (Adams...Zhang) | 6.12 |
| Sikkema lung HCA | 5.5 |
| GongSharma labels uncorrected vs harmony; Leiden corrected → Fig 3 | 4.6 |
| Quantify separation in full G&S (3 covariates) | 4.2 |
| Include Silhouette in benchmark | [x] (benchmark_analysis.rmd:2557) |
| Ecotypes dataset analysis (Rev#1) | 6.1, 6.12 |
| Runtime table shuffled vs ECODA GloProp | 6.7 |
| Leiden clustering benchmark (corrected latent space, G&S, resolution guidance) | 4.6, 6.5 |
| LASSO comparison | 6.8 |
| Translational marker claim tone-down | 6.12 |
| Foundation models (PULSAR), OneK1K/HLCA | 6.9 |
| Modularity NN / cLISI correlation | [x] modularity (scoring_metrics.R); cLISI → 4.2 |
| Expert-annotation generation check | 6.11 |
| Females<40 / males>40 supp fig | 6.10 |
| Marker cell-type stability (jaccard HR vs HiTME, vs Leiden_res_5) | 6.3 |
| scPoli batch-effect observation (adjR2 quote) | discussion context → 6.2/6.12 |
