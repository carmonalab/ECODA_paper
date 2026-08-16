# Onboard 9 new datasets from the PILOT-GM-VAE study (Phase 5)

## Context & goal

Reviewer asked for additional datasets. The chosen source study is PILOT-GM-VAE
(Joodaki et al. 2025, *Brief Bioinform* 26(5):bbaf547, DOI 10.1093/bib/bbaf547,
PMID 41097818 — matches the NAS folder name `JooM_2025_41097818`), whose authors
provide preprocessed, annotated AnnData objects on Zenodo + CellxGene. Full
feasibility table (colors, comments): `new_datasets_to_implement.md` (Excel
source `/Users/christianhalter/Desktop/ECODA_PAPER_DATASETS.xlsx`); the paper's
own dataset description: Zenodo record 15575593 (`datasets.pdf`, saved during
planning).

**Scope**: download 9 of the 12 paper datasets to the NAS, run per-dataset
checks (raw-count sanity, metadata exploration, UMAP batch check), document
each study, then (after user review) register them in datasets.json and roll
out the pipeline. Excluded: Kidney cancer [29] (too few samples), Pancreas
PDAC [31], Myocardial infarction (1) [6] (duplicate of MI-2).

**Decisions already made with the user**:
- Downloads go **directly to NAS** `JooM_2025_41097818/output/` (Mac has only
  ~108 GB free vs ~100–150 GB of data), **sequentially**, one at a time,
  checksum-verified.
- Check notebooks live in **`notebooks/dataset_onboarding/`** (one
  `dataset_check_<Name>.qmd` per dataset + `README.md` summary table).
- The repo **annotation pipeline (scGate/HiTME/scATOMIC) RUNS** for the new
  datasets, with new edge-case safety checks (no-crash when a method annotates
  0/too-few cells), a new datasets.json suitability flag, and Figure-3
  handling in `benchmark_analysis.rmd`.
- `batch_effect_analysis.rmd` funkyheatmap (separation by bio_col vs
  batch_col(s)) is **out of scope** — only TODO.md items get added.
- Diabetes (mouse) and Lung are downloaded but conditional: Diabetes needs a
  mouse-gene decision (T9); Lung is checked by UMAP first (known batch effects).
- **Metrics**: LISI-based separation instead of plain silhouette — the repo's
  min/max-bounded `calc_lisi` (`src/utils/scoring_metrics.R:168`, transform
  `(N-LISI)/(N-1)` → [0,1], 1 = perfect separation, 0 = perfect mixing), which
  already feeds the paper's `calc_sep_score` (`scoring_metrics.R:6-34`).
  Computed at the **cell level** (X_pca/UMAP are cell embeddings; no
  sample-level centroids needed) **per cell type** (author CT column, all CTs
  with ≥ min cells, subsampled ≤1–2k cells per CT for speed/balance) on the
  **unintegrated PCA embedding (top 30–50 PCs)** — metrics on UMAP coords
  distorts distances; NEVER on `X_pca_harmony`/integrated space — for the
  bio label and each batch candidate. Only question answered: "is there a
  batch effect?" — a cell type with high batch separation = batch effect
  signal; low bio separation per CT is expected and fine (cells can be
  identical across conditions while composition differs). scIntegrationMetrics (carmonalab, R/Seurat-based;
  Andreatta et al. 2024, Nat Commun s41467-024-45240-z) and scib (Python) are
  noted as richer options but **deferred**: scIntegrationMetrics is the
  recommended metric suite for the Phase-4 batch-effect analysis (needs a pixi
  setup/install addition → user approval); scib/pylisi are not installed and
  would add a pixi.toml dependency.
- Onboarding notebooks use the **pixi default Python env (scanpy 1.12.2)** and
  do NOT use the R notebook loader (`load_all_functions.R` / the downsized
  20-package `imports.R` from commit `30c4caa` — that loader serves only
  `benchmark_analysis.rmd` / `batch_effect_analysis.rmd`). The LISI helper is a
  standalone Rscript. Verified: `thisutils` 0.4.7 present in the default env
  (`pixi run -e default Rscript`), sklearn 1.9.0 present; `scib`/`lisi`/`pylisi`
  NOT installed (avoid pixi.toml changes unless approved).
- **Metric validation on `_debug` first**: the per-cell-type LISI helper is
  validated against the `_debug` dataset (Joanito 5-sample subset with bio
  `sample.origin` + batch, e.g. `seqtec`) before rolling out to the 9 new
  datasets. Preprocessed `_debug` views available on the NAS at
  `/Volumes/Shared/Projects/ECODA_paper/_debug/output/`
  (`JoaI_2022_35773407_debug_5samples_{benchmark,batch_effect}_analysis_ECODAprocessed.h5ad`).

## Dataset inventory & sources (9 + 0 exclusions)

Author-provided files preferred (they match the paper's Table 1 numbers).
Alzheimer + Parkinson are NOT on Zenodo — the paper obtained them from
CellxGene (must be downloaded from there).

| # | Dataset (short key) | Study (PMID) | Paper: cells / samples / CTs | Author-provided file | Source | Size (GB) | md5 |
|---|---|---|---|---|---|---|---|
| 1 | Alzheimer (SEA-AD) | Gabitto 2024 Nat Neurosci (42486312) | 1,395,601 / 83 / 18 | snRNA h5ad | CellxGene collection `1ca90a2d-2943-483d-b678-b809bf464c30` (≈50 datasets: per-region × assay; pick the 10x 3' v3 nucleus one(s) matching paper counts) | ~15–25 | verify on download |
| 2 | Breast cancer | Kumar 2023 Nature (37380767) | 714,331 / 126 / 10 | `BreastCncr_processed.h5ad` | zenodo 14615923 | 28.9 | `8b28a349c2c3638ddbfb3946a32d12ba` |
| 3 | Covid-19 PBMC | Ren 2021 Cell (34767776) | 993,171 / 151 / 10 | inside `Datasets.tar.gz` | zenodo 8370081 | 36.35 (tar; selective-extract) | `d105b52dbba38ac49c2ffe8b3cf34e24` |
| 4 | Diabetes | Hrovatin 2023 Nat Metab (37697055) — **mouse** | 264,235 / 52 / 13 | `diabetes.h5ad` | zenodo 8370081 | 4.1 | `38189a381bad630fa39ce2d7ad3a0855` |
| 5 | Kidney (KPMP) | Lake 2023 Nature (41648348) | 104,314 / 45 / 14 | `Kidney_KPMP.h5ad` | zenodo 14615923 | 2.75 | `36ceb02ba23c559f80625ec7bef6884f` |
| 6 | Lupus PBMC | Perez 2022 Science (42115607) | 1,263,676 / 261 / 11 | inside `Datasets.tar.gz` | zenodo 8370081 | (in tar) | (in tar) |
| 7 | Lung | Sikkema 2023 Nat Med (42362693) | 941,504 / 165 / 12 | `lungatlas.h5ad` (inside tar.gz) | zenodo 7957118 | 17.2 (tar) | `0d0c97924f1b7a405b6ec3b55da02882` |
| 8 | Myocardial infarction (2) | Kuppe 2022 Nature (41937210) | 132,888 / 23 / 11 | `Myocardial_Infarc_2.h5ad` | zenodo 14615923 | 3.6 | `7431ae99250c99f11bf63e3034798af4` |
| 9 | Parkinson | Prashant 2024 Sci Data (39580497) | 2,096,155 / 97 / 11 | h5ad | CellxGene collection `d5d0df8f-4eee-49d8-a221-a288f50a1590` (single dataset) | ~25–40 | verify on download |

Excluded (NOT downloaded/registered): Kidney cancer (GEO GSE242299), PDAC
(GSA CRA001160), MI(1) (inside the same tar — leave unextracted or move to
`_excluded/`, not registered).

Fallback sources if the author file lacks raw counts (count check FAIL):
GEO GSE158055 (Covid), GSE174188 (Lupus), GSE195665 (Breast), GSE211799
(Diabetes); CellxGene collections `4195ab4c-…` (Breast), `bcb61471-…` (KPMP),
`8191c283-…` (MI); Zenodo 7227571 (Lung original); CellxGene for Alzheimer /
Parkinson (same collections).

## Tasks

### T1 — Download to NAS (user executes prepared script; agent prepares it)
- Write `notebooks/dataset_onboarding/download_datasets.sh` (bash, `curl -L -C -`
  resumable; md5 verification against the table above after each file; log +
  append results to `notebooks/dataset_onboarding/download_log.md`). Run
  **sequentially** (bandwidth + NAS SMB), each verified before the next.
- `Datasets.tar.gz` and `lungatlas.h5ad.tar.gz`: **selective extraction**
  (`tar -xzf … <members>`) — write only the needed h5ads to
  `JooM_2025_41097818/output/`, then delete the tars to save disk.
- CellxGene files: resolve dataset h5ad URLs via the curation API
  (`https://datasets.cellxgene.cziscience.com/<dataset_version_id>.h5ad`);
  for SEA-AD pick the snRNA dataset(s) matching the paper (1,395,601 nuclei,
  83 donors, 18 CTs); verify cell counts and document any delta vs the paper
  (collection was updated 2026-06-10). Keep the `datasets.pdf` description
  file in the onboarding folder for reference.
- Save a record in the download log: source URL, date, size, md5, extraction
  status (feeds the AGENTS.md documentation requirement "when and where it
  was downloaded/retrieved from").

### T2 — Count sanity check (immediately after each download; agent runs)
- Python (`.pixi/envs/default/bin/python`, scanpy 1.12.2 available; backed read
  `sc.read_h5ad(backed='r')` for big files): locate counts
  (`layers['counts']` → `X` → `raw.X`), report slot used; verify sparse matrix
  and **integer VALUES** — do NOT require integer dtype: float-encoded CSR
  (1.0/2.0…) is valid raw counts; check via `np.all(np.mod(sample, 1) == 0)`
  (or `np.allclose(v, np.round(v))`) on a value sample (e.g. 100k cells) with
  epsilon tolerance; **non-negative**; no NaN/Inf; also flag if X looks
  log-normalized (max >> integer range / non-integer fractions in X while
  counts live elsewhere). Report shape, sparsity, n genes, integer-value
  verdict, slot used.
- Verdict PASS/FAIL logged in the onboarding `README.md` table + per-dataset
  qmd. On FAIL → fallback source from the table above (notify user first).

### T3 — Onboarding notebooks (`notebooks/dataset_onboarding/`)
- **T3.0 — Validate the metric + embedding approach on `_debug` first**
  (before writing the 9 notebooks): run `embed_and_umap` + `onboarding_metrics.R`
  (per-cell-type LISI separation on the **unintegrated `X_pca`** — never
  `X_pca_harmony` — bio label `sample.origin` vs batch `seqtec`; per-CT
  subsample caps + confounded-CT guards in effect) on
  `/Volumes/Shared/Projects/ECODA_paper/_debug/output/JoaI_2022_35773407_debug_5samples_batch_effect_analysis_ECODAprocessed.h5ad`.
  Checks: helper runs end-to-end, `thisutils`/`vegan`/`cluster`/`mclust`
  resolve when sourcing `scoring_metrics.R` standalone, per-CT table + heatmap
  produced, PCA input path exercised (incl. the precomputed-`X_pca` case), and
  the result is sensible (batch separation visible within CTs; bio/batch
  expectations from the known Joanito batch structure). Fix the helper until
  it passes, then replicate for the 9 new datasets. Plots for `_debug` also
  serve as a sanity baseline for the UMAP + sample-bar sections.
- Shared helper `onboarding_utils.py` (backed-mode friendly): `locate_counts`,
  `count_sanity_check` (**integer-VALUE check with epsilon tolerance, not
  strict dtype — float-encoded CSR like 1.0/2.0 is valid raw counts**),
  `obs_summary`, `candidate_col_detection` (heuristics for sample /
  biological-label / batch / CT columns), `paper_table_compare`,
  `cells_per_sample_stats`, `embed_and_umap` (subsampled scanpy recipe on
  RAW counts — unintegrated only + precomputed-obsm detection),
  `confounding_crosstab` (bio × batch contingency table + collinearity
  warning), `save_png` — plus a standalone R helper `onboarding_metrics.R`
  (cell-level per-cell-type `calc_lisi` on unintegrated PCA, see section 5;
  **validated on `_debug` before rollout**, see Decisions). One
  `dataset_check_<Name>.qmd` per dataset (9 files; names: `Alzheimer`,
  `Breast_cancer`, `Covid19_PBMC`, `Diabetes`, `Kidney_KPMP`, `Lupus_PBMC`,
  `Lung`, `Myocardial_infarction`, `Parkinson`).
- Notebook sections:
  0. **Study summary card** (agent-filled, from the original paper): title,
     authors, journal, PMID, DOI/link, download source + date + checksum;
     biological groups (with n per group); batch-candidate variables;
     notes (e.g. mouse — Diabetes; brain-region dominance — Parkinson).
  1. **File structure**: n cells × n genes, X / layers / raw summary, obsm
     keys (precomputed `X_umap`/`X_pca`?), obs columns + dtypes.
  2. **Count sanity check** (re-runs T2 with verdict).
  3. **Metadata exploration**: value counts per categorical obs column
     (bounded, top-20); **sorted bar plot of cells per sample**; n samples,
     total n cells, cells-per-sample min/median/mean±sd/max, n samples < 50
     cells; **n unique cell types per CT-annotation column**; **NaN check in
     CT columns** (precedent: Lee/Zhang HiTME NaNs broke PILOT);
     gene-symbol sanity (mouse vs human — Diabetes); **comparison vs
     PILOT-GM-VAE Table 1** (n cells / samples / CTs); **confounding /
     design-matrix check**: cross-tab `table(bio_label, batch_candidate)`
     per batch candidate — if a batch is (near-)perfectly collinear with the
     bio label (e.g. all Disease at Center 1, all Control at Center 2),
     print a loud warning: batch and bio are statistically indistinguishable
     in that dataset, and metrics/UMAP interpretation must account for it.
  4. **UMAP** (unintegrated space ONLY — never Harmony/integrated): use
     precomputed `X_umap` if present (subsample ≤200k for plotting), else
     compute on a subsample ≤100k cells from RAW counts (normalize_total →
     log1p → HVG 2000 → PCA → neighbors → UMAP). Panels: biological label,
     each batch candidate, CT low-res, CT high-res.
  5. **Quick separation metrics** (auxiliary, seconds): cell-level LISI-based
     separation via the repo's own functions — a standalone R helper
     `notebooks/dataset_onboarding/onboarding_metrics.R` sources
     `src/utils/scoring_metrics.R` directly (NOT the notebook loader — it is
     not sourced by these notebooks) and runs `calc_lisi()` on the
     **unintegrated PCA embedding (top 30–50 PCs)** — NOT UMAP coords
     (2D UMAP distorts distances; PCA is the standard input for
     neighborhood metrics per scIB/scIntegrationMetrics) and NEVER
     `X_pca_harmony`/integrated space. Use the global PCA (all subsampled
     cells in the same space); UMAP coords only as documented fallback when
     PCA is unavailable. For the whole subsample AND for **each cell type**
     (author CT column; rand seed fixed): **subsample up to N≈1000–2000
     cells per cell type** (balanced power across abundant and rare CTs,
     keeps kNN/memory trivial regardless of dataset size); separation score
     on (a) the biological label and (b) each batch candidate within the CT.
     Guards: CTs with < min_cells (e.g. 50) skipped (listed); CTs with
     >90% cells from a single batch or single bio group are NOT informative
     → tagged "confounded/uninformative" instead of a score (avoids
     false-positive batch warnings, e.g. rare CT present in one donor);
     single unique label → calc_lisi returns NA (existing behavior, note in
     output). Output: cell-type × label separation table (heatmap) + verdict:
     "batch effect detected in cell types X, Y" (high batch score within a
     CT = batch effect; interpret together with the confounding crosstab —
     collinear designs make bio vs batch indistinguishable). Low bio score
     within a CT is expected and NOT a problem (cells can be identical across
     conditions while composition differs — only "is there a batch effect?"
     matters here). Plots remain the primary deliverable for the user's
     visual inspection. The richer metric suite (CiLISI/norm_cLISI/
     celltype_ASW of carmonalab `scIntegrationMetrics`, and scib in Python)
     is documented as a Phase-4 option, not implemented here (would need pixi
     env additions → user approval).
  6. **Summary + recommendation** (agent): benchmark / batch-effect /
     negative control / exclude, with rationale incl. annotation-suitability
     note (e.g. brain dataset: no/minimal immune cells → HiTME/scATOMIC
     likely unsuitable).
- Plots → `data/new_dataset_checks/<ds>/*.png` (gitignored); optional NAS
  mirror to `JooM_2025_41097818/plots/`.
- `README.md` in the folder: summary table (dataset, study, PMID, file, size,
  checksum, count-check verdict, n cells/samples/CTs, batch candidates, UMAP
  verdict, recommended use, status) + links to `new_datasets_to_implement.md`
  and this plan.

### T4 — Study summaries
- Agent reads each original paper (abstract + methods) for T3 section 0 and
  the metadata-column descriptions; confirms candidate columns against the
  actual obs columns.

### T5 — Docs
- **AGENTS.md**: new section "Onboarding new datasets" — the procedure
  (download author-provided data → count sanity check → onboarding notebook +
  UMAP check → user review → datasets.json registration → pipeline rollout),
  what to document per dataset (study, PMID, groups, batch variables,
  metadata-column description, download source/date), pointers to
  `notebooks/dataset_onboarding/`, the `JooM_2025_41097818` NAS folder, and
  the appendix `new_datasets_to_implement.md`. Write this section against the
  **current** AGENTS.md wording (recent commit `30c4caa` downsized the R
  notebook loader to ~20 pkgs and added `src/utils/env_check.R`; the
  onboarding section must state that these Python check notebooks use the
  pixi default env and do NOT source the R loader).
- **TODO.md**: mark 5.1–5.4 progress; add items:
  - `batch_effect_analysis.rmd`: funkyheatmap like Figure 2a showing separation
    by bio_col and batch_col(s) (multiple batch columns possible) — out of
    scope for this plan;
  - Diabetes mouse-gene preprocessing support (T9);
  - annotation edge-case safety checks + `not_suitable_for_auto_annotation`
    flag + Figure 3 handling (T6);
  - per-dataset annotation-rate documentation (%-annotated cells per method).

### T6 — Annotation safety checks + suitability flag
- `src/4_cell_type_annotation/2.1.1_process_chunk.R`: after each annotation
  method (scGate, HiTME layer2/3, scATOMIC), guard on results — 0 annotated
  cells, <2 unique cell types, or all-NA → **warn + write the column as
  NA/unclassified (no crash)** + record per-method stats (n cells, n types).
  Keep the existing `cutoff.scATOMIC::em` patch and wall-time guards intact.
- `datasets.json` new **optional** field (proposal — flat, like existing
  fields): `"not_suitable_for_auto_annotation": ["hitme", "scatomic"]`
  (absent/empty = suitable for all methods). Consumers:
  - `benchmark_pipeline.R` HiTME/scATOMIC combo guards (currently
    `benchmark_pipeline.R:1398-1416`; scATOMIC combos likewise) — skip the
    flagged methods per dataset;
  - `benchmark_analysis.rmd` Figure 3 (ECODA_automated bars): exclude flagged
    datasets from the HiTME/scATOMIC bars with a documented note (dataset
    stays in Figure 2a where author CTs drive the methods; Leiden-based
    rescue possible later, see below).
- Annotation-rate documentation: after the first HPC annotation run of the
  new datasets, compute %-annotated cells + n unique types per method from the
  merged view h5ads → `notebooks/dataset_onboarding/annotation_summary.json` +
  update the per-dataset qmd section 6 with a short rationale paragraph on why
  HiTME/scATOMIC may be inappropriate (e.g. no/minimal immune cells).
- Fallback for datasets with unusable annotations (future, only if needed):
  unsupervised Leiden-derived CT column via a dataset-specific step like
  `1.4.1_create_kfoury_lowres_ct.R`.

### T7 — Usage-decision checkpoint (user)
- User reviews the UMAP plots + summaries; decides per dataset: benchmark /
  batch-effect / negative control (Parkinson?) / exclude (Lung?, Diabetes?).
  Expected mapping (per the feasibility table): benchmark → Alzheimer, Lupus,
  MI(2), Kidney; batch-effect → Breast, Covid-19, Diabetes (conditional),
  Lung (conditional); Parkinson → decide after check.

### T8 — datasets.json registration (ASK THE USER FIRST — AGENTS.md rule)
- Entries for the approved datasets: `folder_name: "JooM_2025_41097818"`,
  `file_names` (as stored on NAS), `columns` (sample / label / batch / CT
  columns confirmed in T3), `meta_cols_keep`, `views` with `subset_vars`
  (e.g. Diabetes minus embryos E12–E15; Covid-19 PBMC-only if not already
  subset; Kidney already scRNA-only in the author file), `use_for_benchmark`
  / `use_for_batch_effect` per T7, `tissue` / `normal_tissue`, and the
  annotation flag per T6. MI(1), PDAC, Kidney cancer NOT registered.
- Note for views: benchmark views use the **author CT columns** (same as
  existing datasets); batch-effect views use `columns.batch` (one batch col
  per entry; multi-batch handling belongs to the out-of-scope
  `batch_effect_analysis.rmd` funkyheatmap task).

### T9 — Conditional: Diabetes (mouse) pipeline support
- Only if Diabetes passes T3 and the user wants it: add mouse handling to
  `standardize_gene_symbols` (`src/utils/py/gene_utils.py`) + `1.1.1_preprocess.py`
  — detect mouse symbols, map to human orthologs (verify STACAS mouse gene
  conversion availability, e.g. `STACAS::convert_mouse_genes`; else MGI/biomaRt
  orthology table; pixi.toml untouched unless a package is really needed),
  adapt mitochondrial-gene detection (`mt-` prefix). **Requires user sign-off**
  before touching the pipeline.

### T10 — Pipeline rollout (HPC; user runs; agent provides commands)
- Validation first with the two smallest green datasets (Kidney ~2.75 GB /
  104k cells, MI(2) ~3.6 GB / 133k cells), end-to-end:
  `./src/1_stage_data/1_stage_data.sh --ds_name <ds>` →
  `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name <ds>` →
  `./src/4_cell_type_annotation/1_prepare_chunks.sh production` (per-dataset
  chunking; then `2_submit_hpc_array.sh` + `3_submit_merge.sh`) → benchmark
  submitters (`run_r_sample_embedding_methods`, `run_python_sample_embedding_methods`,
  `run_transformation_zeroimp_analysis`) / batch-effect submitters as
  applicable. Verify NAS outputs, then roll out to the remaining datasets.
- Big datasets (Breast 28.9 GB file, Alzheimer/Parkinson 1.4–2.1 M cells):
  preprocess tasks need high-mem nodes — check/raise `--mem` in the
  preprocess submitter before rollout.

## Risks & fallbacks
- **Recent-commit impact (checked 2026-08-16)**: commit `30c4caa` (plan
  `1786902582051-downsize-notebook-imports-loader.md`) downsized the R
  notebook loader (`imports.R` → 20 pkgs) and added `src/utils/env_check.R`.
  No structural impact on this plan: the onboarding notebooks are Python
  (pixi default env) and do not source the R loader; the LISI helper is a
  standalone Rscript sourcing `src/utils/scoring_metrics.R` directly (needs
  `thisutils`/`vegan`/`cluster`/`mclust` — verify at implementation; default
  env currently has `thisutils` 0.4.7). The AGENTS.md section (T5) is written
  against the current post-downsize wording.
- **SEA-AD collection updated** (2026-06-10): counts may differ from the
  paper; verify against Table 1 and document the delta; old dataset versions
  may be retrievable via `dataset_version_id` if a closer match is needed.
- **Author files may lack raw counts** → T2 FAIL → fallback source table
  (above); inform the user before substituting.
- **Huge files** (29 GB h5ad, 2 M cells): backed reads + subsampling in
  notebooks; high-mem nodes on HPC.
- **Annotation failure on new tissues** (brain/heart/kidney/pancreas): safety
  checks + flag (T6); tissue values in datasets.json must map to the worker's
  supported tissue dispatch — validate on Kidney/MI(2) first (T10).
- **NaN in author CT columns** (Lee/Zhang precedent): checked in T3, guarded
  downstream.
- **Diabetes mouse genes** break `standardize_gene_symbols` (human Ensembl105
  map) → T9 conditional; until then Diabetes stays unregistered/benchmark-
  disabled.
- **Disk space**: sequential downloads, selective tar extraction, delete tars
  after verification.

## Validation
- md5 checksums verified per download (Zenodo API values in the table).
- Count sanity check PASS/FAIL logged per dataset.
- Notebooks render (quarto) with plots + summary tables; `README.md` table
  complete.
- Paper-table comparison (cells/samples/CTs) within tolerance.
- **Metric + embedding helpers validated end-to-end on `_debug` first** (T3.0;
  NAS `_debug` batch-effect view; bio `sample.origin` vs batch `seqtec`) —
  helper passes before any of the 9 datasets is processed; per-dataset
  cell-type × label separation tables produced via
  `pixi run -e default Rscript` `onboarding_metrics.R` (thisutils 0.4.7
  verified present).
- User visually reviews UMAPs → usage decisions (T7) → datasets.json approved
  (T8).
- HPC: Kidney + MI(2) end-to-end validation before full rollout (T10).

## Out of scope
- `batch_effect_analysis.rmd` funkyheatmap (bio vs batch_col(s)) — TODO.md
  item only.
- PILOT-GM-VAE/QOT method rollout (separate plan; TODO 3.2).
- Phase 6 items (6.1–6.12), whole-Stephenson-by-center (4.2) extension.
- Excluded datasets (Kidney cancer, PDAC, MI(1)) — not downloaded/registered.
