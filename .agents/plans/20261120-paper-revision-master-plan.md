# ECODA Paper Revision Master Plan (RC-2026-03528)

- **Created:** 2026-11-20
- **Owner:** Orchestrator (`Orchestrator (Neon)`) — **SINGLE-WRITER POLICY:** only the orchestrator edits this plan file. Subagents/workers read this file for context and report findings/status via board/messages to prevent parallel write conflicts.
- **Target Manuscript:** *Cell type composition enables effective patient stratification in single-cell RNA-seq cohorts* (Review Commons `RC-2026-03528`; PLOS Computational Biology transfer)
- **Status:** `✅ ROUND 2 REVISIONS COMPLETE & VERIFIED — 2026-11-26 (6-Cohort Batch / CT-drop, sLISI, 16-Row Classic SF15, PLOS DOCX)`

---

## 1. Goal Summary

Complete all scientific, analytical, figure, table, documentation, and manuscript updates required to address the official journal revision commitments (`revision-plan.docx`) and all point-by-point reviewer comments (`Review_RC-2026-03528_1781259937_1_CH_2.docx`) consistently across:
1. Repository documentation & final publication/supplementary figures.
2. The revised manuscript (`Halter_etal_bioRxiv_ReviewCommons.docx`) — all edits as native Word track changes authored by **`Christian Halter`**.
3. The reviewer response document (`Review_RC-2026-03528_1781259937_1_CH_2.docx`) — green tracked-change responses below each reviewer section, zero edits to existing black/red text or existing Word comments.
4. Progress annotations in `revision-plan.docx` — concise Word comments authored by **`Christian Halter`** only, zero body edits.

---

## 2. Current Progress Dashboard

| Phase | Scope | Status | Gate |
| :--- | :--- | :--- | :--- |
| **Phase 0** | Read-only audit of all 5 `.docx` files, `datasets.json`, `NOTES.md`, `ECODA_PAPER_DATASETS.xlsx`, notebooks & `plots/` + `/grill-me` user alignment (Q1–Q11) | ✅ Done (2026-11-20) | Plan saved |
| **Phase 1** | Read-only Bamboo metadata & Adams marker query + `Breast_cancer` bio/batch cell-type diagnostic + `Lupus_PBMC` strict filtering design | ✅ Done (Gate 1 & Confirmations 1–3 locked 2026-11-21) | Gate 1 approved: keep `Breast_cancer` (7 batch cohorts), `Lupus_PBMC` Option 2 (`majority_only`), `Figure 3C` `Leiden_res_5`, `Supp_fig_16A` `counts_all__0.5`, git push & Bamboo Stage 3 |
| **Phase 2** | `Lupus_PBMC` Option 2 Bamboo rerun (`durable-hpc-gate-ecoda`) + centralize `NOTES.md` + recompute `Supp_fig_16A` (`counts_all_0.5`) + build/polish new & updated figures (`Figure 3C/3`, `Figure 5`, `Supp_fig_X_*`) | ✅ Done & Approved (2026-11-21) | **User Gate 2 approved** (method parenthetical labels + larger BatchMDS text updated; proceeded to Phases 3 & 4) |
| **Phase 3** | Update manuscript (`revision_docs/Halter_etal_bioRxiv_ReviewCommons_revised_20261121.docx`; native `Christian Halter` track changes, dynamic `SEQ`/`REF` for Main Figs 1–7, Supp Figs 1–36 & Supp Tables 1–6, `[ZOTERO: ...]` placeholders, OOXML/`xmllint`/reject-all QA) | ✅ Done (2026-11-21) | **User Review** (`Cmd+A -> F9` field refresh & link 17 `[ZOTERO: ...]` placeholders) |
| **Phase 4** | Populate green `#007A33` tracked-change responses in `revision_docs/Review_RC-2026-03528_1781259937_1_CH_2_revised_20261121.docx` + 13 Word progress comments in `revision_docs/revision-plan_commented_20261121.docx` + OOXML/`xmllint` QA | ✅ Done (2026-11-21) | **Final Release Sign-off** |
| **Round 2** | 16-row Bamboo classic `SF15` profiling (`EPYC 7742` / `H200 NVL`), `Alzheimer` excluded (6 batch / 17 unique cohorts), `CT-drop (|ΔCLR| ≥ 2)`, `sLISI`, `Supp_fig_2` `benchmark_access1133_v1`, figure polish & `*_revised_20261126.docx` | ✅ Done (2026-11-26) | Verified & committed |

---

## 3. Authoritative File & Document Registry

### 3.1 Revision Documents (`revision_docs/`)

| File | Role & Authority | Strict Edit Rules |
| :--- | :--- | :--- |
| `revision_docs/Halter_etal_bioRxiv_ReviewCommons.docx` | Working manuscript draft (currently 5 Main Figs, 22 Supp Figs, 4 Supp Tables, live Zotero CSL fields, Word `SEQ`/`REF` fields on Main Figs & Supp Tables). | **EDITABLE (Phase 3).** Every single edit MUST use native Word track changes (`<w:ins>`/`<w:del>`) authored by **`Christian Halter`**. Preserve all Zotero CSL fields; insert new citations as tracked `[ZOTERO: Author Year]` placeholders. Use dynamic Word `SEQ`/`_Ref*`/`REF` fields for all Main Figs, Supp Figs, and Supp Tables. Validate with `validate.py --original ... --author "Christian Halter"` + LibreOffice PDF render. |
| `revision_docs/revision-plan.docx` | Official journal revision plan (`RC-2026-03528`). Ground truth for 10 in-scope points + 2 out-of-scope points + already-incorporated title change. | **NO BODY EDITS (Phase 4).** Add only very short, concise Word comments (`author="Christian Halter"`) tracking completion & exact manuscript section/figure pointers. |
| `revision_docs/Review_RC-2026-03528_1781259937_1_CH_2.docx` | Raw reviewer comments (black text) with author highlights, red inline notes, and 22 unresolved Word comments. | **APPEND-ONLY RESPONSES (Phase 4).** **NEVER edit** existing black text, red text, highlights, or existing Word comments, and **NEVER add Word comments**. Add responses below each reviewer section/paragraph in **green text with track changes** (`author="Christian Halter"`). Put any internal clarifications/notes in `[square brackets]` inside the green text. May include compact reviewer-only tables/mini-figures where warranted (e.g., R1.4 runtime table). |
| `revision_docs/Thesis.docx` | PhD thesis (`Ch. 2.2`: 3 breast + 1 CRC cancer cohorts, 10 figures, `+1` pseudocount). | **READ-ONLY. NEVER EDIT.** Source narrative only; all quantitative numbers & figures in the paper MUST come from the ported repo `counts_all_0.5` pipeline. |
| `revision_docs/draftECODA_revision_plan_1.docx` | Older, incomplete internal discussion notes (some items marked green). | **READ-ONLY. NEVER EDIT.** Do not treat as ground truth; cross-checked during Phase 0 only. |

> **Pre-Edit Safety Rule:** Before Phases 3–4 begin, confirm no Word lock files (`~$*.docx`) exist in `revision_docs/` and save timestamped untouched backup copies outside the working targets.

### 3.2 Key Repository & External Reference Files

- **Config & metadata ground truth:** `datasets.json` (**do not modify without explicit user confirmation**), `NOTES.md` (central documentation target), `notebooks/dataset_onboarding/README.md` (reconcile conflicts vs `NOTES.md`), `new_datasets_to_implement.md`, `/Users/christianhalter/Desktop/ECODA_PAPER_DATASETS.xlsx`.
- **Analysis notebooks:**
  - `notebooks/benchmark_analysis.rmd` (main benchmark, parameter screen, zero imputation, cell subsetting, modularity kNN, Adams correspondence, Stephenson UCell/CD14)
  - `notebooks/batch_effect_analysis_corrected.rmd` (multi-method batch correction, funkyheatmaps, arrows, MDS)
  - `notebooks/batch_effect_ECODA_becm_comparison.rmd` (ECODA BECM comparison, `sigdiff_ctdrop_<t>` screening)
  - `notebooks/ecoda_application_cancer_clustering.rmd` + `src/notebooks/ecoda_application_cancer_clustering.R` (4-cohort cancer application, `counts_all_0.5`)
- **Derived local data (verified 2026-11-08 to 2026-11-20):**
  - `data/derived/benchmark_reporting_20261110_v2/` (`dataset_statistics.csv`, `core_scores.csv`, `main_ranking.csv`)
  - `data/derived/access_v2/` (includes integer `composition_counts` for all 11 benchmark cohorts + `__trans__` files)
  - `data/derived/ECODA_becm_comp_20261116_v1/comparison_cache.rds` & `data/batch_effect/batch_effect_corrected/p6_batch_majority_current_20261110_v3/bio_vs_batch/`
  - `data/applications/cancer_clustering/20261108-v1/` (frozen H5AD `obs` counts + `Yan/ecoda_data.rds` + `derived/analysis.rds`)
  - `plots/batch_figure_confounding_audit_20261111_v1/` (`Lupus_PBMC` strict vs majority audit)
  - `notebooks/Stephenson_MrVI_research_2026-10-23.md`, `notebooks/Bassez_composition_research_2026-10-23.md`
- **HPC Policy:** Active host is **Bamboo** (`~/ECODA_paper`, `~/scratch/ECODA_paper`). Never copy full H5ADs locally. Phase 1 uses **read-only metadata/marker queries only**. No SLURM/HPC compute jobs may be launched without explicit user approval and the checked-in `durable-hpc-gate-ecoda` workflow.

---

## 4. Settled User Decisions Ledger (Rounds 1–2, 2026-11-20)

1. **Title & Causal Framing:** Current title (*"Cell type composition enables effective patient stratification in single-cell RNA-seq cohorts"*) is finalized (0 `"drives"` remaining); leave the 6 `"driven"` / 2 `"driver(s)"` occurrences untouched.
2. **Scoring Metrics & LISI Definition:**
   - Main benchmark uses **5 metrics**: `ANOSIM`, `Modularity (kNN=3)`, `ARI`, `Silhouette`, and **`LISI`** (all higher = better separation; min–max scaled per dataset/metric before composite mean).
   - **Methods must explicitly define LISI** as normalized inverted cLISI (`LISI_sep = mean(clamp_[0,1]((N - cLISI) / (N - 1)))`, dynamic perplexity `max(2, min(30, floor(0.10 * n_samples)))` via `thisutils::compute_lisi`, `N` = number of biological groups; `1` = perfect separation, `0` = complete mixing). Elsewhere in prose/figures, simplified `"LISI"` is approved.
3. **Zero-Imputation Consistency (`counts_all_0.5`):**
   - Verified across all 4 production pipelines (main benchmark, batch-effect/BECM, cell-subsetting `Supp_fig_X`, cancer clustering): all already add **`+0.5` to ALL counts** (`counts_all_0.5`).
   - **Action:** Fix Methods `P092` (currently mistakenly says "added to zero values") to `counts_all_0.5`, AND recompute **`Supp_fig_16A`** (transformation comparison, currently pinned to `counts_all__1`) locally in pure R with **`counts_all__0.5`** so every ECODA figure is 100% consistent. Leave `EPIC` deconvolution (`(2/3) * min` on zeros) unchanged. **Zero HPC recomputation needed for `counts_all_0.5`.**
4. **Main Benchmark Scope & Synchronization (`Figure 2–4`, `Supp Figs 1–21`):**
   - Synchronize ALL main & supplementary benchmark figures to the refreshed **Nov 10** plots (`11` cohorts, `11` real methods + `Shuffled_baseline` = `12` rows in `Figure 2A`, `18` methods in `Supp Fig 15`, 5 metrics, `QOT` + `PILOT-GM-VAE`, extended Leiden up to `res = 50` + `Leiden_res_2_harmony` in Main `Figure 3`, expanded zero-handling 24 recipes in `Supp Fig 16B`).
   - Fix `Supp Fig 7` (Kim) & `Supp Fig 8` (Lee) alphabetical order and `Supp Fig 18` vs `Fig 3A` dataset count mismatch.
   - Update `Supp Table 2` (`Bassez` 33→40 CTs, `Kfoury` 76,709→76,661 cells, `Stephenson` 55/418,768→53/415,482, `Pelka` `Only immune cells?` Yes→No, verify `Kim`/`Lee`/`Zhang` 1-feature metadata diffs, updated ANOSIM/HVCs, label-independence note / Bassez exception) and `Supp Table 3` (fix `GloProp` typo, add `QOT` & `PILOT-GM-VAE`, add 2 columns: `Originally designed for sample representation` and `Provides batch correction`).
5. **Batch-Effect Cohort Selection & Narrative (`Figure 5`, `Supp Table 4–5`):**
   - **Finalized active set = ALL 7 cohorts** from `ECODA_becm_comp_ANOSIM.pdf`: `Alzheimer` (101), **`Breast_cancer` (33 — RETAINED per Gate 1 decision!)**, `Covid19_PBMC` (44), `Joanito` (124), `Kidney_KPMP_full` (80 — **`condition.l1` evaluation only**, excluding `tissue` evaluation), `Lupus_PBMC` (94 — **Option 2 `cell_filter: "majority_only"`: 374,642 cells, 94/94 samples, min 1,189, 0 multi-batch samples**, Stage 3 + Stage 5 rerun on Bamboo), `Stephenson` (110).
   - **`Breast_cancer` honest limitation framing (Gate 1 decision):** Retain `Breast_cancer` and present results neutrally/transparently, confirming Sikkema et al. HLCA dissociation-protocol compositional batch effects: explicitly note in funkyheatmap discussion that `ECODA_authors` (`limma`) falls below `Shuffled_baseline`, and in the ANOSIM arrow plot explain that two-sided bio/batch imbalance (`12h` 70.6% vs `15m` 43.8% cancer) + overlapping dissociation-sensitive cell types (`pericytes`, `Fibro-SFRP4`, `vsmc`, `Fibro-major`, `Vas-capillary`, `Vas-arterial`, `Lumsec-basal`, etc.) causes bio-blind batch correction to remove both batch and biological separation (inevitable without bio-informed modeling, contrasting with `Stephenson` where remaining unaffected cell types retain signal).
   - **Exclusions:**
     - `Lung`: excluded (post-confounding filtering left only cancer subtypes with ~0 separation across all methods).
     - Legacy `CombinedPBMC` (old `Supp Fig 22`): **dropped** from manuscript and **omitted from `Supp Table 5`**; document in `NOTES.md` and explain transparently in the reviewer response why artificial merging was superseded by real multi-batch cohorts.
     - Other PILOT-GM-VAE / candidate cohorts (`Diabetes`, `Myocardial_infarction`, `Parkinson`, `Lung`, `Kidney cancer`, `Pancreas PDAC`): document centralized exclusion & confounding pre-filtering rationale in **`NOTES.md`** AND distill into **`Supp Table 5`**.
   - **Core scientific message for Batch Section:** Pure sequencing batches (e.g., 3' vs 5') can strongly perturb gene expression while minimally affecting cell-type composition. In contrast, **tissue dissociation / sample preparation protocols** (`Kidney_KPMP_full` & `Breast_cancer`, citing Sikkema et al. HLCA *Nat Med* 2023 Fig 4a, `10.1038/s41591-023-02327-2`) and multi-site processing (`Stephenson`, noting MrVI bioRxiv reviewer discussion) can introduce genuine compositional batch effects, which batch-aware ECODA / correction mitigates when biological signal is orthogonal (`Stephenson`, `Kidney`), whereas severe bio-batch confounding (`Breast_cancer`) limits bio-blind correction; note `Alzheimer` is primarily resolved by MrVI.
   - **Batch figures:** Drop `(A)` uncorrected-only and `(D)` internal single-cell batch geometry from paper. **Main `Figure 5`** = combined 7-cohort batch-effect corrected funkyheatmap + cleaned ANOSIM before/after correction arrows (without `Lung` & without `Kidney` tissue). **Supplementary** = cleaned `ECODA_becm_comp_ANOSIM` (title/subtitle moved to legend, polished design) + one **4-column MDS figure per retained batch cohort (7 figures)** (`Batch Uncorrected | Batch Corrected` ∣ `Bio Uncorrected | Bio Corrected`, vertical middle divider, harmonized group/subplot titles).
   - **Point 10 ("≥5 new datasets / ≥16 cohorts total"):** Satisfied cleanly via `11` clean benchmark cohorts + `6` new non-Stephenson batch cohorts + `Yan` cancer cohort (**18 unique cohorts total**).
6. **Cancer Application (`Thesis.docx` Ch. 2.2 → Paper, Point 2):**
   - **Placement:** New **final Results/Application section immediately AFTER `scECODA`** (and before `Discussion`).
   - **Data & figures:** Strictly use the **ported `counts_all_0.5` / current-obs pipeline** (`notebooks/ecoda_application_cancer_clustering.rmd`, 34 shared HiTME subsets across `Yan` N=99, `Zhang` tumor N=13, `Bassez` N=26, `Pelka` tumor N=63).
     - **`Figure 6` (Main):** `breast_cross_cohort_pca_and_stratification.pdf` (3-breast PCA biplots, response boxplots, consensus loadings).
     - **`Figure 7` (Main):** `bassez_hitme_shared_heatmap.pdf` (Bassez clustered immune-module heatmap showing lymphoid-high / activation-low non-expanders).
     - **Supplementary:** `appendix_yan_zhang_hr_heatmaps.pdf`, `cross_cohort_all_datasets_pc1_pc2_clustered_correlation_matrix_heatmap.pdf`, and `pelka_hitme_shared_heatmap.pdf`, plus **`Supp Table 6`** (cancer cohorts summary).
     - **Omitted:** Thesis Figs 2.1–2.4, 2.8 (forest plot), 2.9A–D, and 2.2 (Bassez marker dotplot — omitted since Adams diagnostic plots are included).
   - Note concisely in Methods/Data Availability that `Yan` was processed with identical HiTME reference maps/scGate models (`ecoda_data.rds` archived) and `NAC == "NA"` reflects the early-switch `Low` group per author clarification.
7. **Other `revision-plan.docx` Items (Points 3, 4, 5, 7, 9) & Minor Comments:**
   - **Point 3 (Downsampling):** Existing `Supp_fig_X_ECODA_authors_HR_cell_subsetting.pdf` (all 11 cohorts × 5 metrics, `2000 → 50` cells/sample) is sufficient; add short text in `"ECODA is robust to different cell type annotation strategies"` (or check fit in HVC section) and place figure in Supplementary.
   - **Point 4 (MrVI / transcriptional complementarity):** Add concise text in Benchmark Results + Discussion explaining MrVI captures global IFN-response transcriptional programs in Stephenson (citing Stephenson et al. & noting MrVI authors also evaluated Stephenson); refresh **both** `CD14_LR_pseudobulk_hvg2000_PCA` and `UCell_IFN_UCell_authors_highres` on current **53-sample** Stephenson cohort as a Supplementary Figure.
   - **Point 5(iii) (Adams marker consistency — Gate 1 & Confirmation 1 Updated):** Switch Main **`Figure 3C`** Adams unsupervised PCA panel from `Leiden_res_0.4` to **`Leiden_res_5`** (consistent with `Leiden_res_5` `ANOSIM 0.648` vs `0.445` at `res 0.4`); export full **`Authors_HR × Leiden_res_5` (`38 × 92`) correspondence matrix as a Supplementary Excel table (`.xlsx`)** (dropping the arbitrary `res 0.4` appendix matrix); and match **`Leiden_res_5`** clusters with canonical Adams et al. (2020) markers and `Authors_HR`/`HiTME`/`Authors_LR` in the Supplementary canonical marker-gene dotplot.
   - **Point 7 (PULSAR):** Benchmark includes `PILOT-GM-VAE` + `QOT`; exclude `PULSAR` and explain primarily via **PBMC-trained domain limitation** (secondarily computational cost) in reviewer response & Discussion.
   - **Point 9 & R1.12 (Complete Gong & Sharma):** Point to existing **Figure 1B–C** (already all 868 samples with CMV, age, sex); no new Stage 5 run.
   - **R1.4 (Shuffled vs ECODA runtime):** Explain sub-second timing parity (`GloProp` 0.581s, `Shuffled_baseline` 0.638s, `ECODA` 0.753s mean) with a compact 3-row table in the green reviewer response.
8. **Dynamic Numbering & Citations:**
   - Convert **ALL Supplementary Figures** (existing & new) to **dynamic Word `SEQ`/`_Ref*`/`REF` fields** matching Main Figures & Supp Tables.
   - Use **temporary IDs (`Supp_fig_TMP_*`)** during drafting and perform a **final sequential renumbering pass at the end of Phase 3** strictly ordered by first mention in the main text.
   - New citations: tracked `[ZOTERO: Author Year]` placeholders + DOI table (Section 8) for user linking in Word/Zotero.

---

## 5. Ground-Truth Traceability Matrix (`revision-plan.docx` & Reviewers)

| Plan Item | Reviewer Ref(s) | Target Manuscript Location | Target Figure / Table | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Title** (*"enables effective"*) | R2 #6 | Title (`P000`) | — | ✅ Already in `.docx` |
| **Point 1** (Label circularity) | R1.1, R1.5, R1.9 | Benchmark Results + Discussion | `Supp Table 2` (label provenance / Bassez exception) | ⬜ Phase 3 |
| **Point 2** (De novo discovery / Cancer) | R1.2, R1.9, R1.10 | NEW Final Results section (post-`scECODA`) + Methods | `Figure 6`, `Figure 7`, `Supp_TMP_Cancer_*`, `Supp Table 6` | ⬜ Phase 3 (PDFs ready) |
| **Point 3** (Small-sample / cell downsampling) | R1.7, R3.2 | Annotation Results (or HVC) + Methods | `Supp_TMP_CellSub` (`Supp_fig_X_..._cell_subsetting.pdf`) | ⬜ Phase 3 (PDF ready) |
| **Point 4** (MrVI / transcriptional complementarity) | R2.1, R3.11 | Benchmark Results + Discussion | `Supp_TMP_StephIFN` (53-sample refresh) | ⬜ Phase 3 (PDF ready) |
| **Point 5(i–iii)** (Leiden Harmony, extended res 50, Adams markers) | R2.2–R2.5 | Annotation Results + Methods | `Figure 3` (with `2_harmony`, res 50, & `Leiden_res_5` PCA), `Supp Fig 18`, `Supp_TMP_AdamsConsist` (`.xlsx` + dotplot) | ⬜ Phase 3 (PDFs/table ready) |
| **Point 6** (Expanded zero-handling range) | R2.6 | Benchmark Results + Methods + Supp Note | `Supp Fig 16B` (24 recipes) + `Supp Fig 16A` (`counts_all_0.5`) | ⬜ Phase 3 (PDFs ready) |
| **Point 7** (PILOT-GM-VAE, QOT; PULSAR rationale) | R2.11, R3.5 | Benchmark Results + Discussion + Methods | `Figure 2A`, `Supp Figs 1, 2, 15, 17`, `Supp Table 3` | ⬜ Phase 3 (PDFs ready) |
| **Point 8** (Multi-batch cohorts + Silhouette & LISI) | R2.12, R2.14, R3.1, R3.7, R3.9, R3.10 | Batch Results + Methods + Metrics | `Figure 5`, `Supp_TMP_BECM`, `Supp_TMP_BatchMDS_*`, `Supp Table 4–5` | ⬜ Phase 3 (PDFs/tables ready) |
| **Point 9** (Unfiltered cohorts / full G&S) | R1.12, R3.8, R3.10 | Batch Results + Figure 1B–C reference | `Figure 1B–C` (868 samples), `Figure 5` | ⬜ Phase 3–4 |
| **Point 10** (≥5 new datasets / ≥16 cohorts) | R2.12, R3.6, R3.10 | Results + Batch + Cancer + Methods | 11 benchmark + 5–6 batch + `Yan` = 17–18 unique cohorts (`Supp Tables 2, 4, 5, 6`) | ⬜ Phase 3 |
| **OOS-1** (Supervised LASSO vs HVC) | R2.7 | Reviewer doc + brief Discussion note | `Supp Fig 21` | ⬜ Phase 4 |
| **OOS-2** (Compositional methods on same HVC) | R3.4 | HVC Results clarification + Reviewer doc | Clarify all compositional methods used all CTs | ⬜ Phase 3–4 |
| **Minor items** | R1.3, R1.4, R1.6, R1.8, R1.10, R1.11, R1.13, R2.8, R2.9, R2.13, R3.3, R3.7 | Throughout text + Methods + Discussion + Reviewer doc | `Supp_TMP_ModKNN` (`Supp_fig_X2.pdf`), R1.4 3-row runtime table | ⬜ Phase 3–4 |

---

## 6. Temporary Figure & Table Registry (Final Sequential Numbering Locked at End of Phase 3)

### 6.1 Main Figures (Dynamic `SEQ Figure`)
- **`Figure 1`**: Concept + Gong & Sharma 868-sample cohort *(unchanged)*
- **`Figure 2`**: Benchmark ranking (11 cohorts × 11 methods + `Shuffled_baseline`, 5 metrics incl. Silhouette & LISI) + Bassez MDS *(refresh from Nov 10 `plots/Figure_2.pdf`)*
- **`Figure 3`**: Annotation granularity & method (incl. `Leiden_res_50` + `Leiden_res_2_harmony`) + Adams PCA (**updated in Phase 2A to `Leiden_res_5` for consistency**)
- **`Figure 4`**: HVC variance explained + Adams/Kfoury 2-HVC PCA & log-ratios *(refresh from Nov 10 `plots/Figure_4.pdf`)*
- **`Figure 5` (REDESIGNED)**: Multi-cohort batch effects — **Panel A**: 7-cohort batch-effect corrected funkyheatmap (excl. `Lung` & `Kidney` tissue) + **Panel B**: cleaned 7-cohort ANOSIM before/after correction arrows *(built in Phase 2B)*
- **`Figure 6` (NEW)**: Cancer application — 3-breast-cohort cross-cohort 34-immune PCA, response boxplots, and consensus loadings (`plots/ecoda_application_cancer_clustering/20261108-v1/breast_cross_cohort_pca_and_stratification.pdf`)
- **`Figure 7` (NEW)**: Cancer application — Bassez clustered 34-immune heatmap revealing lymphoid-high / activation-low non-expanders (`plots/ecoda_application_cancer_clustering/20261108-v1/bassez_hitme_shared_heatmap.pdf`)

### 6.2 Supplementary Figures (Convert ALL to Dynamic `SEQ Supplementary_Figure`; Final Numbers Assigned in First-Mention Order in Phase 3)

| Temp ID | Content / Source | Planned First Mention in Main Text |
| :--- | :--- | :--- |
| `Supp_TMP_Metrics` (was `1`) | Metric consistency across 5 metrics (`plots/Supp_fig_1.pdf`) | Benchmark Results (metrics intro) |
| `Supp_TMP_ModKNN` *(NEW)* | Modularity kNN sensitivity: `k = 3, 6, 9, √n` (`plots/Supp_fig_X2.pdf`) | Benchmark Results (immediately after `Supp_TMP_Metrics`) |
| `Supp_TMP_ParamScreen` (was `2`) | Parameter screening (`plots/Supp_fig_2_methodparameter_screening.pdf`, GM33-refreshed) | Benchmark Results (parameter sensitivity) |
| `Supp_TMP_MDS_Adams..Zhang` (was `3–13`) | Per-cohort MDS across all 11 benchmark cohorts (`plots/Supp_fig_3..13_*.pdf`, alphabetical: **Kim before Lee** fixed) | Benchmark Results |
| `Supp_TMP_Runtime` (was `14`) | Execution time & memory (`plots/Supp_fig_14.pdf`) | Benchmark Results |
| `Supp_TMP_Extended` (was `15`) | Extended 18-method benchmark (`plots/Supp_fig_15_benchmark_extended.pdf`) | Benchmark Results |
| `Supp_TMP_TransZero` (was `16`) | Transformations (**A recomputed with `counts_all_0.5`**) + expanded 24-recipe zero-handling (**B** `Supp_fig_16B_*`) | Benchmark Results + Supplemental Note |
| `Supp_TMP_Unscaled` (was `17`) | Unscaled 5 metrics per dataset (`plots/Supp_fig_17_*.pdf`) | Benchmark Results |
| `Supp_TMP_StephIFN` *(NEW)* | Stephenson 53-sample CD14-monocyte pseudobulk PCA + IFN UCell PCA (refreshed in Phase 2) | Benchmark Results (short sentence on MrVI / Stephenson IFN complementarity) |
| `Supp_TMP_AnnotMetrics` (was `18`) | Annotation strategies (`mod/ari/sil/lisi`, extended Leiden res `50` + `2_harmony`) (`plots/Supp_fig_18_*.pdf`) | Annotation Results |
| `Supp_TMP_AnnotAutoVsLeiden` (was `19`) | Bassez, Lee, Stephenson auto vs Leiden MDS (`plots/Supp_fig_19_Bassez_Lee_Stephenson.pdf`) | Annotation Results |
| `Supp_TMP_AdamsConsist` *(NEW)* | Adams `Leiden_res_5` canonical marker-gene dotplot + **Supplementary `.xlsx` `Authors_HR × Leiden_res_5` (`38 × 92`) correspondence matrix** (Phase 2A) | Annotation Results |
| `Supp_TMP_CellSub` *(NEW)* | Cell-per-sample downsampling `2000 → 50` across 11 cohorts × 5 metrics (`plots/Supp_fig_X_ECODA_authors_HR_cell_subsetting.pdf`) | Annotation Results (or check fit in HVC Results) |
| `Supp_TMP_HVC_MeanVar` (was `20`) | Mean-variance CLR plots highlighting HVCs (`plots/Supp_fig_20_HVCs_per_dataset.pdf`) | HVC Results |
| `Supp_TMP_HVC_Metrics` (was `21`) | HVC `mod/ari/sil/lisi` (`plots/Supp_fig_21_*.pdf`) | HVC Results |
| `Supp_TMP_BECM` *(NEW)* | Cleaned ECODA BECM comparison (`ECODA_becm_comp_ANOSIM`, title/subtitle moved to legend, polished in Phase 2) | Batch Results |
| `Supp_TMP_BatchMDS_1..7` *(NEW, 1 per retained cohort)* | 4-column per-cohort batch MDS (`Batch Uncorr | Batch Corr` ∣ `Bio Uncorr | Bio Corr`, center divider, harmonized headers) | Batch Results |
| `Supp_TMP_Cancer_YanZhang` *(NEW)* | Yan & Zhang clustered 34-immune heatmaps (`appendix_yan_zhang_hr_heatmaps.pdf`) | Cancer Application Results |
| `Supp_TMP_Cancer_CRC` *(NEW)* | Cross-cancer 8-axis loading correlation + Pelka CRC clustered heatmap (`cross_cohort_..._heatmap.pdf` + `pelka_hitme_shared_heatmap.pdf`) | Cancer Application Results |

*(Legacy `Supp Fig 22` `CombinedPBMC` is removed. Thesis Fig 2.2 Bassez marker dotplot is omitted.)*

### 6.3 Supplementary Tables (Dynamic `SEQ Supplementary_table`)
- **`Supp Table 1`**: Gong & Sharma 868-sample cohort (verify values).
- **`Supp Table 2`**: 11 Benchmark cohorts (sync with `dataset_statistics.csv`; fix `Bassez` 40 CTs, `Kfoury` 76,661 cells, `Stephenson` 53 samples / 415,482 cells, `Pelka` `Only immune cells?` → `No`; verify `Kim`/`Lee`/`Zhang` 1-feature diffs; updated ANOSIM/HVCs; add label-independence note / Bassez exception).
- **`Supp Table 3`**: Benchmarked methods overview (fix `GloProp` typo; add `QOT` & `PILOT-GM-VAE`; add columns `Originally designed for sample representation` & `Provides batch correction`).
- **`Supp Table 4`**: Retained multi-batch cohorts (7 cohorts; Bamboo-verified `#samples`, `#cells`, `mean ± SD`, `#CTs`, batch type [`Sequencing` / `Tissue dissociation` / `Multi-site`], pre-filtering rules).
- **`Supp Table 5` (NEW)**: Candidate multi-batch cohort screening & exclusion ledger (`Diabetes`, `Myocardial_infarction`, `Parkinson`, `Lung`, `Kidney cancer`, `Pancreas PDAC`; **omit `CombinedPBMC`**).
- **`Supp Table 6` (NEW)**: Cancer application cohorts (`Yan` N=99, `Zhang` tumor N=13, `Bassez` N=26, `Pelka` tumor N=63; 34 shared HiTME subsets, clinical endpoints).

---

## 7. Detailed 4-Phase Execution Plan

### Phase 1 — Data Verification & Decision Gates *(Read-only Bamboo + Local Diagnostics; Zero `.docx` Edits; Zero HPC Compute)*

- [x] **Task 1.1 — Lightweight Read-Only Bamboo Metadata & Adams Marker Query:**
  - Exact post-Stage-3 stats extracted from Bamboo across all 8 batch views (`data/derived/phase1_metadata_audit_20261120/batch_cohort_metadata.{json,csv}`), verifying identical `uncorrected` vs `corrected` counts (`Alzheimer` 101 / 1,305,022; `Breast_cancer` 33 / 114,579; `Covid19_PBMC` 44 / 254,174; `Joanito` 124 / 252,519; `Kidney_KPMP_full` 80 / 283,290; `Lupus_PBMC` 94 / 458,785 [`layer1` 0 non-NA]; `Lung` 80 / 244,454; `Stephenson` 110 / 543,080).
  - Resolved `Kim`/`Lee`/`Zhang` 1-CT diff (`data/derived/phase1_metadata_audit_20261120/benchmark_metadata_audit.json`): metadata counts include `"NA"` / unassigned cells (`Kim` literal `"NA"` category: 50 vs 49; `Lee` 37,829 NA codes: 43 vs 42; `Zhang` 24,548 NA codes: 42 vs 41), whereas `ECODA_authors_HR` retains 49 / 42 / 41 valid cell types. Updated `Supp Table 2` sync values recorded for `Bassez` (40 CTs, ANOSIM `0.461`), `Kfoury` (76,661 cells, ANOSIM `0.425`), `Lee` (`0.256`), `Stephenson` (53 samples, 415,482 cells, ANOSIM `0.313`, top2 `Mono_prolif, CD16_mono`), and `Zhang` (ANOSIM `0.756`, top2 `Mono_CD14, CD8.MAIT`).
  - Extracted Adams H5AD markers on Bamboo via chunked CSR scan (`data/derived/adams_marker_dotplot_20261120/`): 62/63 canonical markers matched (`CD68` absent; `C1QA`/`MARCO`/`FABP4`/`MSR1` present), and top enriched markers verified for `Leiden_res_0.4` clusters `1, 6, 11, 12, 15, 18, 20` matching `adams_correspondence_top_clusters.csv` (Cluster 15 = `VE_Peribronchial`/`VE_Capillary_B`/`VE_Arterial`; Cluster 12 = `ATII`/`Club`/`Basal`/`Aberrant_Basaloid`; Cluster 11 = `Myofibroblast`/`Fibroblast`; Cluster 18 = `SMC`/`Pericyte`; Cluster 6 = `Macrophage_Alveolar`; Cluster 1 = `Macrophage`).
- [x] **Task 1.2 — `Breast_cancer` Bio-vs-Batch Cell-Type Diagnostic (`notebooks/Breast_cancer_batch_bio_diagnostic_2026-11-20.md`):**
  - Audited all 58 `Authors_HR` cell types across the 33-sample cohort (`12h`: 12 cancer / 5 normal [70.6%]; `15m`: 7 cancer / 9 normal [43.8%]) in `data/derived/phase1_metadata_audit_20261120/breast_cancer_bio_vs_batch_celltypes.csv`.
  - Found uncorrected `batch ANOSIM = 0.379` (>2× `bio ANOSIM = 0.183`), driven by 11 BH-significant dissociation types (`pericytes`, `Fibro-SFRP4`, `vsmc`, `Vas-capillary`, `Fibro-major`, `Vas-arterial`, `CD4-Tem`, `Lumsec-basal`, `Lumsec-prol`, `cDC2`, `bmem_unswitched`), 7 of which overlap with marginal cancer-vs-normal differences; within-`12h` stromal/mural differences disappear or reverse, while `limma` collapses `bio ANOSIM` (`0.183 → 0.094`, `Harmony_limma` `0.251 → 0.024`, `Pseudobulk` `0.150 → 0.073`) and `droplow50pct` drops `#1` biological hit `CD4-activated` (`p_bh = 0.00092`). Recommended **dropping `Breast_cancer`** (6 clean cohorts keep ≥17 unique cohorts / ≥6 new cohorts).
- [x] **Task 1.3 — `Lupus_PBMC` Strict Cell-Level Filtering Design:**
  - User confirmed switching `Lupus_PBMC` to strict 7-batch cell filtering within the exact 94-sample allowlist (`384,718` cells, 94/94 samples, min `1,189` ≥ 500 cells) — requiring targeted Bamboo rerun.
  - Audited exact implementation touchpoints (`datasets.json`, `src/utils/py/cohort_selection.py`, `src/utils/py/export_h5ad_sample_metadata.py`, `src/utils/py/audit_corrected_h5ad_source.py`, `tests/test_batch_majority_contract.py`, `tests/test_cohort_selection.py`): scoped dataset policy flag (`cell_filter: "subset_majority_keys_v1"`, approved only for `Lupus_PBMC`) selects the exact 94 samples via `majority_v1` first (preventing re-entry of excluded donors), then intersects cell rows with `evaluate_subset_mask` on `batch_cov`, yielding 100% single-batch samples (`n_samples_with_multiple_raw_batches = 0`) and clean Bamboo Stage 3 (2 views) + Stage 5 (`Lupus_PBMC` batch only) rerun into a new derivative root under `durable-hpc-gate-ecoda` (Stage 4 skipped: `not_suitable_for_auto_annotation`).

> **✅ USER DECISION GATE 1 RESOLVED (2026-11-21):**
> 1. **`Breast_cancer`**: **RETAINED** (7 batch cohorts) with neutral/honest limitation framing in `NOTES.md`, `Figure 5` funkyheatmap/arrow discussion, and text.
> 2. **`Lupus_PBMC` Option 2 (`cell_filter: "majority_only"`)**: Approved (374,642 cells, 94/94 samples, min 1,189, **0 multi-batch samples**), including scoped commit/push + `durable-hpc-gate-ecoda` Stage 3 & Stage 5 rerun on Bamboo.
> 3. **Main `Figure 3C` & Adams consistency (`Leiden_res_5`)**: Approved switching `Figure 3C` to `Leiden_res_5`, exporting `38 × 92` correspondence matrix as `.xlsx`, and matching `Leiden_res_5` clusters in the canonical marker dotplot.
> 4. **`Supp_fig_16A` (`counts_all__0.5`)**: Update canonical `transformation_analysis_policy()` in `benchmark_pipeline.R` to `"counts_all__0.5"` and regenerate locally.

---

### Phase 2 — Centralized Docs & Final Figure Generation

- [x] **Task 2.1 — Centralize Cohort Documentation in `NOTES.md`:**
  - Added canonical section `## Canonical Paper Revision Cohort Registry & Decisions — 2026-11-21` in `NOTES.md` documenting: (a) final 7 batch-effect cohorts (`Breast_cancer` retained honestly with bio-blind dissociation limitation), (b) confounding pre-filtering rules per cohort (`Lupus_PBMC` strict `majority_only`), (c) `Kidney_KPMP_full` tissue-dissociation classification & HLCA Fig 4a citation, (d) `Breast_cancer` cell-type bio/batch diagnostic outcome, (e) complete exclusion rationale for candidate cohorts (`Diabetes`, `Myocardial_infarction`, `Parkinson`, `Lung`, `Kidney cancer`, `Pancreas PDAC`) and legacy `CombinedPBMC`, and (f) reconciliations with `notebooks/dataset_onboarding/README.md`.
- [x] **Task 2.2 — Recompute `Supp_fig_16A` with `counts_all__0.5`:**
  - Updated `transformation_analysis_policy()` in `src/5_run_benchmark_methods/benchmark_pipeline.R` (`4a531ff`) and recomputed via `src/notebooks/recompute_supp_fig_16A_counts_all_0.5.R` (bundles & summary in `data/derived/transformation_counts_all_0.5_20261121/`).
  - Regenerated & visually verified `plots/Supp_fig_16A_transformations_anosim.pdf`, `plots/Supp_fig_16A_transformations_metric_companions.pdf`, and `plots/Supp_fig_16AB_trans_zero_anosim.pdf` (plus individual metric panels; zero star clipping).
- [x] **Task 2.3 — Generate & Polish New / Updated Local Figures (Visual QA Passed):**
  - **`Figure 3` (`Leiden_res_5` Adams consistency package):** Updated `notebooks/benchmark_analysis.rmd` Figure 3C chunk to `ECODA_seuratres_5` (`Leiden_res_5`, top-5 PCA loading clusters `40, 60, 44, 47, 53`), regenerated & visually verified `plots/Figure_3.pdf` via `src/notebooks/regenerate_figure3_leiden5.R`; exported 4-sheet correspondence matrix `plots/Supp_Table_Adams_AuthorsHR_vs_LeidenRes5.xlsx` (`src/notebooks/build_adams_leiden5_xlsx.R`, `305,701` cells, `92` clusters, `100` samples); built & visually verified `plots/Supp_fig_TMP_Adams_marker_dotplot_leiden5.pdf` (`src/notebooks/build_adams_marker_dotplot_leiden5.R`).
  - **`Supp_TMP_StephIFN`:** Refreshed both label-blind `CD14` pseudobulk PCA (ANOSIM `0.489`) and `IFN_UCell` `authors_highres` PCA (ANOSIM `0.218`, PC1 `91%`) on the current **53-sample** Stephenson cohort (`plots/Supp_fig_TMP_Stephenson_IFN.pdf`, `src/notebooks/build_stephenson_53_complementarity.R`).
  - **`Figure 5`, `Supp_TMP_BECM`, & `Supp_TMP_BatchMDS_1..7`:** Built & visually verified all 9 7-cohort batch figures (`src/notebooks/build_batch_7cohort_figures_20261121.R`, `data/derived/batch_7cohort_release_20261121/`, `data/derived/ECODA_becm_comp_20261121_v2/`) incorporating the strict `majority_only` `Lupus_PBMC` Stage 5 results: `plots/Figure_5.pdf` (+ standalone `Figure_5_A_*.pdf` & `Figure_5_B_*.pdf` + 5 companion arrow plots under `plots/batch_effect_corrected_7cohorts_20261121/`), `plots/Supp_fig_TMP_ECODA_becm_comp_ANOSIM.pdf` (cleaned 2×4 layout, title/subtitle moved to caption, non-overlapping `ggrepel` cutoffs), and 7 4-column MDS files `plots/Supp_fig_TMP_BatchMDS_{Alzheimer,Breast_cancer,Covid19_PBMC,Joanito,Kidney_KPMP_full,Lupus_PBMC,Stephenson}.pdf`.
  - **7-cohort corrected mean-scaled ranking (`batch_7cohort_ranking_corrected.csv`):** `MrVI` (`0.838`), `GloScope` (`0.719`), `ECODA_Harmony_Leiden_CLR` (`0.640`), `ECODA_authors` (`0.610`), `ECODA_Harmony_Leiden_CLR_limma` (`0.603`), `Pseudobulk` (`0.543`), `QOT` (`0.493`), `PILOT` (`0.299`), `Shuffled_baseline` (`0.123`).
- [x] **Task 2.4 — `Lupus_PBMC` Strict `majority_only` Bamboo Rerun (`durable-hpc-gate-ecoda`):**
  - **Scoped commits pushed to `origin/master`:** `4a531ffcbe8d53d5989f9bd34acdf237f9aab801` (`Lupus_PBMC.batch_metadata_policy.cell_filter = "majority_only"`, `counts_all__0.5`, tests) and `16bfb302cf10c6403aaaffc2c5023cda9cc2d216` (`resolve_cohort_selection()` post-filter sample-order fix + regression tests).
  - **Approved Stage 3 gate (`ecoda_stage3_lupus_strict_20261121_v2`):** Slurm jobs `4802040/4802041/4802042` (`COMPLETED|0:0`), `ecoda_run_audit.sh` (`ECODA_RUN_AUDIT_OK`), reviewer approved (`release_eligible: true`), NAS `666` & MD5 verified — uncorrected `e10d61f033a09032b2b31ffdeb29b1b6`, corrected `8fc265128af35ac572b81154df8f0757`; both `374,642` cells, `94` samples (`46 Case / 48 Healthy`), min `1,189`, `max_batch_per_sample = 1`, `0` multi-batch samples.
  - **Approved Stage 5 batch campaign (`ecoda_stage5_lupus_strict_batch_20261121_v1`):** Children `stage5_lupus_strict_uncorrected_20261121_v1` & `stage5_lupus_strict_corrected_20261121_v1`, derivative release `20261121-lupus-strict-v2`, `36/36` scheduler jobs `COMPLETED|0:0`, both `ECODA_RUN_AUDIT_OK`, NAS `checksums.md5` 100% OK + `666` permissions, reviewer approved (`release_eligible: true`).
  - **Local scoped sync:** `4.3 MB` sample-level results & metadata only (no H5ADs/pseudobulks) into `data/derived/lupus_strict_stage5_20261121/`.
  - **Key strict vs old `Lupus_PBMC` shifts (`lupus_old_vs_strict_comparison.csv`):** `ECODA_authors` corrected ANOSIM `0.348 -> 0.350` (batch ANOSIM `0.094 -> -0.053`), `ECODA_Harmony_Leiden_CLR_limma` Mod k=3 `0.340 -> 0.513`, `MrVI` ANOSIM `0.354 -> 0.311` & ARI `0.313 -> 0.022`, `Pseudobulk` ARI `0.179 -> 0.307`.

> **✅ USER DECISION GATE 2 APPROVED (2026-11-21):** User inspected all Gate 2 figures, requested parenthetical batch-method labels plus larger `BatchMDS` typography (`src/notebooks/build_batch_7cohort_figures_20261121.R`, regenerated & verified), and authorized proceeding directly through Phases 3 and 4.

---

### Phase 3 — Manuscript Revision (`revision_docs/Halter_etal_bioRxiv_ReviewCommons_revised_20261121.docx`)

- [x] **Task 3.0 — Pre-flight Backup & Safe Output Policy:** Original `Halter_etal_bioRxiv_ReviewCommons.docx` held an active Word lock (`~$lter_etal_bioRxiv_ReviewCommons.docx`). Saved untouched MD5-verified backups under scratch (`docx_backups_20261121/original_md5s.txt`) and wrote all three deliverables as non-conflicting dated drafts (`*_revised_20261121.docx` / `*_commented_20261121.docx`).
- [x] **Task 3.1 — Tracked-Change Prose Updates (`author="Christian Halter"`):**
  - *Abstract & Introduction*: Updated to 11 clean cohorts (`695` samples), 7 multi-batch cohorts (`586` samples), 4 solid-tumor discovery cohorts (`201` samples; `18` unique cohorts total), 5 separation metrics, `ECODA` vs `scECODA`, `QOT`, `PILOT-GM-VAE`, and *de novo* unsupervised stratification.
  - *Results §1–2 (Benchmark)*: Verified protocol label independence in `10/11` clean cohorts (`Bassez` T-cell expansion exception; `ST2`), `11` methods + `Shuffled_baseline` (`0.051`), `5` metrics (`ANOSIM`, `Modularity kNN=3`, `ARI`, `Silhouette`, `LISI`), `SF01`, `SF02` (`kNN` sensitivity `r = 0.93–0.98`), `SF03` parameters, `SF15` runtime (`0.581 s` GloProp, `0.638 s` Shuffled, `0.753 s` ECODA, `192.5 s` QOT, `2968.4 s` / `3.89 h` max PILOT-GM-VAE) & peak RAM, and `SF18` Stephenson (`n=53`) MrVI–ECODA interferon complementarity (`CD14+` PB ANOSIM `0.489`, `IFN_UCell` PC1 `91%`).
  - *Results §3 (Annotation strategies)*: Extended Leiden up to `res 50` (`845.5` clusters, ANOSIM `0.416`, `p < 0.05`) + `Leiden_res_2_harmony` (`0.412`), documented Gong & Sharma Leiden nuance (`0.477/0.299/0.300` vs `0.695` Authors_HR), `Figure 3C` + `SF22` + `Supplementary Data 1` (`38 × 92`) Adams `Leiden_res_5` marker convergence (`40/44/60` profibrotic vs `47/53` alveolar macrophages), and `SF23` cell-per-sample downsampling (`2,000 -> 50`, `90.8%` ANOSIM retention at `500` cells).
  - *Results §4 (HVCs)*: Synced `Figure 4`, `SF24–25` (all 5 metrics), clarified OOS-2 (all competitors used `100%` of cell types), softened clinical flow-cytometry claim (`Discussion` + `Results`), and added rare/condition-specific HVC caveat.
  - *Results §5 (Multi-cohort Batch Evaluation)*: Retired legacy `CombinedPBMC`, rewrote around 7 multi-batch cohorts (`586` samples, `3.13M` cells; `Figure 5`, `ST4`, `ST5`, `SF26–SF33`), distinguishing sequencing-chemistry invariance (`Joanito`) from sc/sn & dissociation shifts (`Kidney KPMP`, `Breast cancer`, `Stephenson`), honest bio-blind limitation on `Breast cancer` (`ECODA_authors (limma)` below `Shuffled_baseline`, `SF28`), strict `majority_only` `Lupus PBMC` (`n=94`, 0 mixed samples, `SF32`), and MrVI complementarity on `Alzheimer` (`SF27`).
  - *NEW Final Results Section (after `scECODA`, before `Discussion`)*: *"De novo exploratory patient stratification across solid tumor cohorts"* (`Figure 6`, `Figure 7`, `ST6`, `SF34–SF36`) across 201 primary tumors (`Yan n=99`, `Zhang n=13`, `Bassez n=26`, `Pelka CRC n=63`; `34` shared HiTME subsets; `counts_all_0.5`), showing conserved PC1/PC2 axes (`r = 0.507–0.833`) and unsupervised clinical/MMRd stratification.
  - *Discussion*: Added transcriptional vs compositional complementarity, Bassez exception & de novo cancer validation, annotation subjectivity caveat, rare/low-variance HVC caveat, softened clinical biomarker framing + CD8/response citations (`Tumeh et al. 2014`, `Yan et al. 2026`), and dedicated `PULSAR` discussion (`Pang et al. 2025`, PBMC-domain limit & foundation-model compute).
  - *Methods & Supplemental Note*: Fixed `P092` (`counts_all_0.5`), added exact `Silhouette` (`cluster::silhouette`) & normalized inverted `cLISI` (`LISI_sep`), `QOT`, `PILOT-GM-VAE`, `Leiden_res_50` & `Leiden_res_2_harmony`, Gong & Sharma full `n=868` `Figure 1B–C`/`ST1` cross-ref, HVC vs HVG rationale, downsampling, batch correction methods & contingency screens, cancer methods (`Yan` `ecoda_data.rds`), and updated `SF17` (`counts_all_0.5` + expanded `0.5–200` & `multLN 0.001%–5%` range).
- [x] **Task 3.2 — Supplementary Tables (`ST1–ST6`) & External `Supplementary Data 1`:** Updated `ST2` (11 cohorts), `ST3` (`11` methods + GloProp typo fix + populated cols 6–7), `ST4` (7 retained batch cohorts; legacy `Gong & Sharma`/`Zhu` track-deleted), added landscape `ST5` (6 excluded cohorts, explicitly omitting legacy `CombinedPBMC`), added landscape `ST6` (4 cancer cohorts), and documented external `Supp_Table_Adams_AuthorsHR_vs_LeidenRes5.xlsx` as `Supplementary Data 1`.
- [x] **Task 3.3 — Dynamic Figures, `SEQ`/`REF` Fields & Monotonic Numbering:** Embedded `Fig2..Fig7` and `SF01..SF36` (240/180 DPI PNGs); converted all captions and in-text citations to dynamic Word `SEQ`/`_Ref*`/`REF _Ref* \h` fields; fixed Kim (`SF08`) vs Lee (`SF09`) order; set `<w:updateFields w:val="true"/>`; zero `Supp_TMP_` tokens remain.
- [x] **Task 3.4 — Automated OOXML Validation:** Verified via `scratch/docx_pipeline/validate_revision_docx.py` (`ALL_DOCX_VALIDATION_CHECKS_PASSED`): `xmllint --noout` across all `24` XML/RELS parts; `new_ins=676`, `new_del=244` (all `Christian Halter`); `71` unique `wp:docPr` IDs & all media present; `reject_all` preserves all `73` Zotero fields & `22` `<m:oMath>` equations; `SEQ Figure=7`, `Supplementary_table=6`, `Supplementary_Figure=36`; `104` resolved `REF _Ref*` bookmarks; strictly monotonic `1..36` first-mention SF sequence in main text. *(Note: `soffice`/`libreoffice` is not installed locally, so full XML/XSD/reject-all/PIL figure QA ran cleanly; open in Microsoft Word and press `Cmd+A -> F9` to refresh cached field display numbers.)*

---

### Final Numbering Crosswalk (Phase 3 / Phase 4 Locked)

- **Main Figures:** `Figure 1` (workflow & Gong & Sharma n=868) | `Figure 2` (11-cohort 5-metric benchmark & Bassez MDS) | `Figure 3` (annotation counts, strategies incl `Leiden_res_50` & `Leiden_res_2_harmony`, Adams `Leiden_res_5` PCA) | `Figure 4` (HVC variance selection & CD4 driver ratios) | `Figure 5` (7-cohort batch benchmark ranking & bio/batch arrow plot) | `Figure 6` (*de novo* breast cross-cohort PCA & clinical stratification) | `Figure 7` (Bassez hierarchical immune-module heatmap).
- **Supplementary Tables:** `ST1` (Gong & Sharma `n=868`) | `ST2` (11 benchmark cohorts & label origins) | `ST3` (sample representation methods) | `ST4` (7 retained batch cohorts) | `ST5` (6 excluded candidate cohorts; no `CombinedPBMC`) | `ST6` (4 cancer discovery cohorts) | **`Supplementary Data 1`** (`plots/Supp_Table_Adams_AuthorsHR_vs_LeidenRes5.xlsx`).
- **Supplementary Figures (`SF01–SF36` in strict first-mention order):**
  `SF01` Metric consistency (5 metrics) -> `SF02` Modularity kNN sensitivity (`k=3,6,9,√n`) -> `SF03` Parameter screening -> `SF04–SF14` Clean benchmark MDS (`Adams, Bassez, Gong & Sharma, Kfoury, Kim [alphabetical fix], Lee, Pelka, Smillie, Stephenson, Wu, Zhang`) -> `SF15` Method runtime & peak RAM -> `SF16` Extended 18-method benchmark -> `SF17` Transformations (`counts_all_0.5`) & expanded zero-handling (`0.5–200`, `multLN`) -> `SF18` Stephenson (`n=53`) CD14+ PB & `IFN_UCell` -> `SF19` Unscaled 5 metrics -> `SF20` Annotation Mod/ARI/Sil/LISI (`11` cohorts, `res 50`, `2_harmony`) -> `SF21` Bassez/Lee/Stephenson auto vs Leiden MDS -> `SF22` Adams `Leiden_res_5` canonical marker dotplot -> `SF23` Cell-per-sample downsampling (`2000->50`) -> `SF24` HVC mean-variance plots -> `SF25` HVC Mod/ARI/Sil/LISI curves -> `SF26` ECODA BECM batch mitigation summary -> `SF27–SF33` 4-col paired BatchMDS (`Alzheimer, Breast cancer, COVID-19 PBMC, Joanito, Kidney KPMP, Lupus PBMC, Stephenson`) -> `SF34` Yan & Zhang tumor heatmaps -> `SF35` Cross-cancer PC1/PC2 eigenvector correlation matrix -> `SF36` Pelka CRC heatmap.

---

### Phase 4 — Reviewer Response & Revision Plan Tracking

- [x] **Task 4.1 — Populated `revision_docs/Review_RC-2026-03528_1781259937_1_CH_2_revised_20261121.docx`:**
  - Inserted `504` green (`#007A33`) tracked-change (`<w:ins w:author="Christian Halter">`) point-by-point responses below every section/point across Reviewers #1, #2, and #3 plus compact 3-tier runtime table (`Fast 0.58–9.89 s`, `Medium 28.5 s–17.8 min`, `Heavy 49.5 min / 3.89 h max`). Zero `[ZOTERO:` tokens remain in reviewer doc; internal author checks placed in `[square brackets]`.
- [x] **Task 4.2 — Commented `revision_docs/revision-plan_commented_20261121.docx`:**
  - **Zero body edits** (100% character parity with original); anchored `13` concise `"Christian Halter"` (`CH`) Word comments across `P10..P49` (Points 1–10, Title, OOS-1, OOS-2) with modern Word `commentsExtended`/`commentsIds`/`commentsExtensible`/`people.xml` parts.
- [x] **Task 4.3 — Final Verification (`ALL_DOCX_VALIDATION_CHECKS_PASSED`):**
  - Reviewer response: `xmllint --noout` (`16` XML/RELS parts), `word/comments.xml` byte-identical SHA-256 (`22` existing comments, `0` added), `reject_all` plain text `100%` identical to backup, all inserted runs/cells `#007A33`.
  - Revision plan: `xmllint --noout` (`21` XML/RELS parts), body plain text `100%` identical to backup, `13` balanced comments.
  - All `6` original files in `revision_docs/` (`original_md5s.txt`) + `datasets.json` untouched.

---

## 8. New Zotero Citations Checklist (17 Unique Placeholders / 33 Occurrences in Manuscript)

Replace `[ZOTERO: ...]` tracked placeholders in `revision_docs/Halter_etal_bioRxiv_ReviewCommons_revised_20261121.docx` via Word's Zotero plugin:

| Placeholder Token in Revised Manuscript | Author / Subject | Reference DOI / Source |
| :--- | :--- | :--- |
| `[ZOTERO: Gabitto et al. 2024]` | SEA-AD Alzheimer brain MTG atlas | `10.1038/s41593-024-01774-5` |
| `[ZOTERO: Hrovatin et al. 2023]` | Mouse pancreatic islet atlas (excluded candidate) | `10.1038/s42255-023-00876-x` |
| `[ZOTERO: Joodaki et al. 2025]` | PILOT-GM-VAE (`Briefings in Bioinformatics`) | `10.1093/bib/bbae713` |
| `[ZOTERO: Kamath et al. 2022]` | Human substantia nigra / Parkinson (excluded candidate) | `10.1038/s41593-022-01061-1` |
| `[ZOTERO: Korsunsky et al. 2019]` | Harmony integration & LISI | `10.1038/s41592-019-0619-0` |
| `[ZOTERO: Kumar et al. 2023]` | Human Breast Cell Atlas (`Breast_cancer` batch cohort) | `10.1038/s41586-023-06252-9` |
| `[ZOTERO: Kuppe et al. 2022]` | Human myocardial infarction (excluded candidate) | `10.1038/s41586-022-05060-x` |
| `[ZOTERO: Lake et al. 2023]` | KPMP Kidney Atlas (`Kidney_KPMP_full` batch cohort) | `10.1038/s41586-023-05769-3` |
| `[ZOTERO: Pang et al. 2025]` | PULSAR multi-scale foundation model | Verify Pang et al. bioRxiv/publication entry in Zotero |
| `[ZOTERO: Perez et al. 2022]` | Lupus PBMC multiplexed cohort (`Lupus_PBMC`) | `10.1126/science.abf1970` |
| `[ZOTERO: Prashant et al. 2024]` | Parkinson multi-cohort snRNA-seq study (excluded candidate) | Verify exact Prashant et al. entry in Zotero |
| `[ZOTERO: Ren et al. 2021]` | COVID-19 PBMC large multi-study atlas (`Covid19_PBMC`) | `10.1016/j.cell.2021.01.053` |
| `[ZOTERO: Rousseeuw 1987]` | Silhouettes graphical validation (`J Comput Appl Math`) | `10.1016/0377-0427(87)90125-7` |
| `[ZOTERO: Sikkema et al. 2023]` | Human Lung Cell Atlas (HLCA, Fig 4a dissociation & exclusion) | `10.1038/s41591-023-02327-2` |
| `[ZOTERO: Tumeh et al. 2014]` | PD-1 blockade induces responses via pre-existing CD8+ T cells | `10.1038/nature13954` |
| `[ZOTERO: Wang Z et al. 2024]` | QOT & GloProp sample distance / divergence methods | Verify against existing Ref in `Bibliography` / `10.1093/bioinformatics/btae597` |
| `[ZOTERO: Yan et al. 2026]` | Early-stage TNBC neoadjuvant chemotherapy response cohort | `10.1038/s41586-026-10469-9` |

---

## 9. Round 2 Key Updates & Decisions (2026-11-26)

1. **Batch cohort exclusion (`Alzheimer`):** Moved `Alzheimer` (`n = 101`) from `Supp Table 4` to `Supp Table 5` due to low biological separation across all methods (uncorrected ANOSIM `0.059–0.110`, corrected `0.018–0.143`), leaving **6 retained multi-batch cohorts** (`485` samples; **17 unique cohorts total**) and `SF27–SF32` (`BatchMDS` without score subtitles; internal scored copies archived under `plots/batch_mds_internal_with_scores_20261126/`).
2. **ECODA batch mitigation & combined ranking (`Figure 5`, `SF26`):** Replaced ECODA `limma` and `droplow50pct` in `Figure 5` and `BatchMDS` with **`CT-drop (|ΔCLR| ≥ 2)`** (`ECODA_authors` `0.737`, rank 2 behind `MrVI` `0.760` under joint 16-row min-max scaling; `ECODA_Leiden` `Harmony` `0.666`; `ECODA_Leiden` `Harmony + CT-drop` `0.637`), while retaining `limma` alongside `|ΔCLR| ≥ 0..4` in `SF26`.
3. **Metric & framing updates:** Standardized **`sLISI`** everywhere (native Word equation in Methods, `cLISI` restricted to Methods only); replaced Bassez "sole exception" framing with longitudinal post-treatment outcome wording; simplified `SF18` (`Supp_fig_TMP_Stephenson_IFN.pdf`) to single-panel `CD14+` monocyte pseudobulk PCA (ANOSIM-only comparison); removed `"ns"` across `Figure 3`, `SF17`, `SF20` and added Wilcoxon asterisks (`hide.ns = TRUE`) to `SF19` and `SF23`.
4. **`Supp_fig_2` & `SF15` profiling harmonization:** Switched `Supp_fig_2` to `data/derived/benchmark_access1133_v1/analysis_release.rds` to include all 4 `PILOT-GM-VAE` parameter-screening variants; ran scoped 16-row classic-full-H5AD Bamboo profiling (`13 AMD EPYC 7742 CPU + 3 NVIDIA H200 NVL GPU`, `16/16 COMPLETED|0:0`, `data/derived/sf15_classic_profile_20261126/`) resolving `Stephenson`/`Bassez` `SF15` outliers, and documented Bamboo CPU/GPU hardware plus `PILOT-GM-VAE` CPU routing (`batch_size = 32` + CPU `ot.emd2`) in Methods and `NOTES.md`.
5. **Deliverables & final repository sync:** Generated `revision_docs/Halter_etal_bioRxiv_ReviewCommons_revised_20261126.docx` and `revision_docs/Review_RC-2026-03528_1781259937_1_CH_2_revised_20261126.docx` (native `Christian Halter` track changes, PLOS-style Arial typography, zero em dashes in tracked insertions, `SF01–SF35` monotonic numbering), validated OOXML/reject-all/R contract tests, and committed and pushed all changes to `origin/master`.
