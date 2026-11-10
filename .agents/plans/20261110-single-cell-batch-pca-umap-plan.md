# Single-cell/nucleus batch PCA and UMAP — approved pilot first

## Current authority — Q14–15 supersede the earlier expanded draft

**APPROVED_PILOT_FIRST,2026-11-10.** The user approved implementation, with a timed first task and an Alzheimer memory check before bulk expansion. Q15 reduces the initial scope to **12 main figures and12uniqueUMAP fits**: one2×4PCA and one2×5UMAP per cohort, comparing actual production UNC PCA with COR Harmony. Use broad annotations for five cohorts, source `louvain` for Lupus, and each row's own cached Leiden2 labels. Matched-HVG, fine/source-named and fixed-anchor variants are **deferred**, not scheduled or scored. The72-figure/18-fit proposal below is retained only as design history and deferred options.

- First task: Breast UNC, all114,579cells, time UMAP, separate ANN90/50D LISI and rendering individually. No bulk launch until measured costs are reviewed.
- Alzheimer: inspect selected obs/PC/graph shapes/dtypes and estimate memory including UMAP/ANN overhead before allocation. Load only selected arrays into expression-free AnnData; never load the whole RNA/count matrix or whole parent AnnData.
- Q12: display all cells for cohorts≤200,000; otherwise at most100,000, balanced across actual per-cell batches without replacement. Retain small batches in full and redistribute unused quota. Identical selected IDs across both representations/panels. All cells remain in fitting/scoring.
- Q13 clarifies1,000cells/sample was the unused alternative, not1million.
- Accepted technical contract: standard Scanpy UMAP settings/seed0; validated approximate90-neighbor LISI at perplexity30; normalized cellwise means on both figure types with raw/scaled/equal-Sample summaries saved; gray missing labels and valid-label-only cLISI with explicit coverage; separate expression-free H5ADs on HPC/NAS; bounded display extracts for local PDFs.
- Finish GM33 data release before diagnostic compute. Gong's sole waiter68872 completed18:45:07UTC and its terminal/native/scientific/NAS/official reviews passed; no more model fitting is needed there. New source still requires normal committed snapshot/runtime binding and durable execution. Do not overwrite existing parents, modify datasets.json or add runtime packages.

The remaining text records the availability audit, scientific contracts and earlier expanded proposal. Where it says Q14 pending,72figures,18fits/rows,three representations, or future approval, this current section overrides it: approval is **pilot first**, initial scope is **two representations/cohort**, and additional variants remain deferred. Implement only the smaller current scope; do not build a generalized variant scheduler.

### Pilot implementation checkpoint

- Implementer `docs_checkpoint_audit` completed the fixed BreastUNC Python worker, scalar shell entry point and narrow native-audit PURPOSE branch. No bulk array, package change, parent write or datasets.json change. Six atomic pilot products contain embedding-only H5AD, full-population metrics, bounded display metadata, two one-row diagnostic PDFs and phase wall-time/RSS measurements. These diagnostic PDFs measure rendering and are not the final two-row cohort deliverables.
- Source paths: `src/notebooks/run_singlecell_batch_geometry.py`, `src/utils/bash/singlecell_batch_geometry.sh`, `src/utils/bash/ecoda_run_audit.sh`, `tests/test_singlecell_batch_geometry.py`.
- Author's four focused tests pass29.109s; independent reviewer `benchmark_dependency_audit` accepts the exact frozen files, with two independent checks passing27.112s. Coverage includes native Harmonypy LISI parity, missing/degenerate labels, real expression-free UMAP with forbidden expression/PCA/graph refits, display-ID/coordinate alignment, atomic outputs, Alzheimer header-only reads, copied Slurm-spool execution, native ownership and failure evidence.
- Fixed pilot resources are8CPUs/64G/2h/shared-cpu, with no CPU-model/GPU restriction. Stage3 selected-artifact receipt/review and exact observed parent stats are bound in controlSHAac7bf25e30a0ad652f4695b4ff9a3a9c7113bde595a0a2c6ea788843cfd7c8e1. Root must commit/push only task-owned files, create a clean commit-keyed snapshot, stage the three small control files, run actual control/header preflight, then use one durable gate/wait/terminal review. Existing downstream dirty files remain excluded.
- Implementation is source-reviewed; no pilot measurement or bulk-compute success is claimed at this checkpoint. NAS publication uses the existing finite backup helper after actual output review.

The availability audit and earlier expanded options below support the approved pilot. They do not enlarge the current12-figure scope.

## Responsibility and current boundary

The Supervisor assigns implementers and reviewers. Root is the orchestrator/operator. Anyone with an unclear role asks the Supervisor before acting. This draft does not transfer ownership of existing analysis sessions or downstream figure work.

Gong fitting/review is complete. Its remaining GM33 metadata projection/release retains a separate run ledger and sole waiter; this task adds no second monitor or scheduler query to it.

## Requested figures

For each of six curated batch cohorts:

| Figure | Rows | Columns, left to right |
|---|---|---|
| PCA, 2 × 4 | Uncorrected/pre-Harmony; Harmony | Sample; author cell type; biological label; batch |
| UMAP, 2 × 5 | Uncorrected/pre-Harmony; Harmony | Sample; author cell type; Leiden resolution 2; biological label; batch |

These are cell/nucleus-level points, not sample embeddings. “Raw PCA” means PCA of normalized/log-transformed/scaled expression before Harmony, not PCA of raw counts. Annotation variants follow the user decisions below; Lupus source-cluster and source-named-label semantics must remain explicit.

Round1 selects **two comparisons**: primary actual production UNC raw versus COR Harmony; and a separate matched-HVG comparison of COR pre-Harmony versus COR Harmony. Reuse the **three unique PC representations and three UMAP fits per cohort**, rather than fitting the shared COR Harmony UMAP twice.

For UMAP, show each row's own resolution-2 Leiden labels in the primary display and fixed uncorrected-anchor labels on both rows in a separate alternative display. Reuse identical coordinates between label variants; valid stored Leiden assignments require no refit. Record the anchor identity and align labels by cell ID. Each annotation variant has two PCA figures plus four UMAP figures per cohort. The user now selects **two annotation variants for all six cohorts:72 combined figures**, reusing the same **18 unique UMAP fits**.

Planning sizes are **506 samples and 2,928,159 cells/nuclei**. The completed UNC6/COR6 header checks confirm the corresponding coordinate row counts and50-column PC shapes; ordered-ID/value-level validation remains separate.

| Cohort | Samples | Reported cells/nuclei |
|---|---:|---:|
| Alzheimer | 101 | 1,305,022 |
| Breast_cancer | 33 | 114,579 |
| Covid19_PBMC | 44 | 254,174 |
| Joanito | 124 | 252,519 |
| Lupus_PBMC | 94 | 458,785 |
| Stephenson | 110 | 543,080 |

Stephenson's **110-sample batch cohort is distinct from its 53-sample benchmark cohort**. All six have `use_for_batch_effect=true` and both batch views; only Stephenson also participates in the benchmark. Top-level `datasets.json` columns remain authoritative; no configuration change is planned.

| Cohort | Configured sample | Author LR / HR | Biology, display only | Configured batch |
|---|---|---|---|---|
| Alzheimer | `donor_id_assay` | `Subclass` / `Supertype` | `Cognitive status` | `assay` |
| Breast_cancer | `sample_id` | `broad_cell_type` / `author_cell_type` | `disease` | `assay`, `suspension_dissociation_time` |
| Covid19_PBMC | `sampleID` | `majorType` / `celltype` | `CoVID-19 severity` | `datasets` |
| Joanito | `sample.ID` | `cell.type` / `cell.type_new` | `sample.origin` | `seqtec` |
| Lupus_PBMC | `sampleID` | `layer1` / `louvain` | `Status` | `batch_cov` |
| Stephenson | `Sample` | `initial_clustering` / `full_clustering` | `Status` | `Site` |

Batch outputs also contain canonical `Sample`. Its mapping must remain tied to the configured raw sample column and source witness, rather than being standardized again. Both UNC/COR header queries establish that Lupus `layer1` is **absent** and `louvain` has24 numeric declared categories. The user has explicitly selected Lupus `louvain`/`cell_types` as two diagnostic annotation variants; this is not a production configuration substitution or change.

## Established provenance and source behavior

The maintained Stage3 selection contains **17 rows: eight UNC, eight COR, and Stephenson benchmark**. Producer run: `stage3_confounds_20261108_v3`; immutable source: `7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3`. Source manifest SHA256: `4d9a8d9a75ecbbe9926d130a88f3c40e90d83627816770790baa91e80d843c69`; runtime identity SHA256: `1288531120165abdfe1b50edf5632a8e61cddf7f2656260335539f1201d58231`.

The selected16 batch artifacts were independently validated and approved in `.gate/stage3_noop_controls_20261109_v1/batch_views.json` (SHA256 `c8e1d9b5b0a4b485ecdb228c706c6f8a3de65f0861fc57b148de68cd39bb05b4`) and `batch_review.json`. The original global Stage3 gate remains **FAILED**; selected-artifact approval does not rewrite that outcome. Scratch/NAS identities and selected OK ownership were checked in that receipt.

For all six requested cohorts, UNC/COR have identical ordered cell-ID/canonical-Sample digests and retained-sample order. Their membership hashes differ because view configuration is part of that hash. Authoritative COR parents are under:

`/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage3_confounds_20261108_v3/artifacts/datasets/<dataset>/output/<configured corrected-view filename>`

| COR parent | Actual recorded MD5 |
|---|---|
| Alzheimer | `11387ee22bd059c3b15f73fa1c9b250f` |
| Breast_cancer | `ace473009bce4d61aefa70acc2f44938` |
| Covid19_PBMC | `30491bb345424a1e0f627096a9835761` |
| Joanito | `81293c230eff77615bca97ad6220c1d4` |
| Lupus_PBMC | `f26d8170cd10d2afd2928cfefb1fb06e` |
| Stephenson | `e697cdcf2bd859d4673ee7e8e9b148e3` |

Pinned `1.1.1_preprocess.py` (SHA256 `97db667758b730b72ce63f5cca86a9b91096443a8d30371df5b1c2eece902d1e`) establishes:

- A COR file stores `X_pca_batch_effect_corrected_hvg2000` and applies Harmony to that same PCA, storing `X_pca_harmony_batch_effect_corrected_hvg2000`. Thus one COR parent can isolate Harmony on the **same cells and same HVG2000 feature selection**.
- Separate UNC preprocessing ranks HVGs by `Sample`; COR ranks by configured batch, using an ephemeral tuple composite for multiple batch columns. UNC PCA versus COR Harmony would compare different HVG-selection regimes as well as correction.
- PCA requests 50 components, capped by `n_vars - 1` and `n_obs - 1`; clustering uses all stored PCs. HVG2000 is the clustering pass. Expected resolutions are 0.1, 0.4, 2, 5, 20, 50, with igraph Leiden, two iterations and an undirected graph.
- Expected COR resolution-2 keys are `leiden_res_2_batch_effect_corrected_hvg2000` and `leiden_res_2_batch_effect_corrected_hvg2000_harmony`, with corresponding raw/Harmony neighbors graphs.
- The maintained producer does **not** compute UMAP. Input embeddings may have been inherited; their presence would not establish current subset/graph provenance. The selected comparisons therefore retain the plan for18fresh, unique UMAP fits, with implementation/launch approval still pending.

The initial COR6 query (`.gate/curated6_singlecell_schema_query_20261110_v2/schema.observed.v1.json`, SHA256 `9306ce94ab18874e267027eb46ebb8cce0b1e23528921d84664d27ec9067b45c`) was followed by successful paired structural queries:

- `.gate/curated6_singlecell_schema_query_20261110_v3/schema.COR6.observed.v1.json`, SHA256 `58be55584a99e4392b5d181fd08d5503bfc7d4be813a487d9edbd629bb6cdb8f`.
- `.gate/curated6_singlecell_schema_query_20261110_v3/schema.UNC6.observed.v1.json`, SHA256 `3d4d2296052a275f72c26cd0a9e0353c0f844a68f375c6f76f3cf2eda962b20b`.

All18 unique PC matrices are50-dimensional. Each cohort has all three corresponding graphs and resolution-2 Leiden keys. The production UNC keys are `X_pca_batch_effect_uncorrected_hvg2000`, `neighbors_batch_effect_uncorrected_hvg2000`, and `leiden_res_2_batch_effect_uncorrected_hvg2000`; the two COR representations/keys are listed above. Neighbor metadata records15neighbors,50PCs, Euclidean distance and seed0. Its `method='umap'` describes neighbor-graph construction, not proof of a maintained fitted UMAP. These are observed structural facts, **not** coordinate-finiteness, cellwise missingness, ordered-ID, or categorical-code validation.

Only Alzheimer, Breast and Lupus contain parent `X_umap` keys; they are inherited and lack maintained UMAP provenance. They are not substitutes for the selected new fits. The other cohorts do not establish a parent UMAP option in this query.

Declared author-category counts from the header query are shown below. These are **dictionary levels, not counts of used/nonmissing cell labels**.

| Cohort | Declared HR levels | Declared LR levels |
|---|---:|---:|
| Alzheimer | 131 | 24 |
| Breast_cancer | 58 | 10 |
| Covid19_PBMC | 54 | 10 |
| Joanito | 12 | 11 |
| Lupus_PBMC | 24 numeric `louvain` levels | `layer1` absent |
| Stephenson | 50 | 18 |

### Lupus annotation provenance — not resolved by category names

Both selected views contain `cell_types` with11 declared named levels: `B`, `NK`, `PB`, `Progen`, `Prolif`, `T4`, `T8`, `cDC`, `cM`, `ncM`, `pDC`. `ct_cov` has14 named lymphoid/progenitor levels; its vocabulary lacks monocyte/DC labels, but cellwise coverage/missingness has **not** been read. `pop_cov` is ancestry (`African American`, `Asian`, `European`, `Hispanic`), not cell type. `L3` has binary `0.0`/`1.0` levels, not an established cell-type annotation. Numeric `louvain` clusters are not interchangeable with named author cell types or the new resolution-2 Leiden labels.

The local provenance trace establishes **source-provided columns**, not an ECODA annotation mapping:

- `notebooks/dataset_onboarding/download_sources.sh:56–60` selects Lupus/Perez/GSE174188/SLE from the standardized Zenodo8370081 `Datasets.tar.gz` archive (recorded archive MD5 `d105b52dbba38ac49c2ffe8b3cf34e24`) as `Lupus_Perez2022.h5ad`; Stage1 copies that standardized file from the JooM collection. The Stage3 input witness records MD5 `001658910686c61a5010da95b7b14a15`.
- The maintained Stage2 step table has no Lupus-specific hook. Pinned Stage3 and its H5AD loader preserve source `obs` columns while subsetting/filtering and adding canonical Sample/new representation keys; neither defines a `cell_types`/`ct_cov` remapping or automated annotation step. Current Stage3/loader bytes match producer7e67 (`97db6677...` / `da1d3a44...`).
- `notebooks/dataset_onboarding/dataset_specs.py:211–234` explicitly calls `cell_types`/`ct_cov` “raw author” annotations, distinct from an old onboarding intention for HiTME `layer1`/`layer2`. That intention is not evidence that Pipeline4 ran on these batch views. The batch exception remains in force.

The repository does not carry the original Perez annotation dictionary or the standardized archive's label-creation mapping. Consequently the defensible current description is **source-provided named `cell_types`**, with original-author versus upstream-curator provenance unestablished. Do not claim original-Perez author ground truth solely from names or silently replace configured `layer1`/`louvain`. The user has chosen `louvain` for the main diagnostic and `cell_types` for the alternative; original-Perez genealogy is **not** an additional approval barrier, while cellwise coverage still needs validation. This new-figure selection authorizes neither a production `datasets.json` edit nor reprocessing of already-reviewed batch methods.

## Round1–2 decisions and Round3 frontier

1. **Q1 answered:** primary actual production UNC raw versus COR Harmony, plus separate matched-HVG COR pre-Harmony versus COR Harmony. Reuse the common Harmony representation/UMAP.
2. **Q2 answered:** both Leiden displays: each row's own labels in the primary version, fixed uncorrected-anchor labels across both rows in a separate alternative. Coordinates are shared; independently generated cluster numbers do not establish cross-row identity.
3. **Q3 answered:** Breast plots and batch iLISI use **dissociation time only**. The actual H5 header has one declared assay category, consistent with constant `10x 3' v3`. Keep both stored/configured correction columns unchanged; this display/metric choice is not a `datasets.json` change. `sequencing_platform` remains outside correction covariates.
4. **Q4 answered:** create a **new embedding-only H5AD on HPC/NAS**, without copying expression `X`, counts, or raw matrices. The storage choice is approved; its exact schema/provenance and implementation still need agreement.

The user also requests **batch iLISI and cLISI against the chosen annotation**. Raw iLISI increases with batch mixing; raw cLISI near1 indicates label purity. These are cell-level integration diagnostics, separate from the benchmark's `lisi_score` and its normalized inverse biological-label definition. Compute/report once per distinct representation/annotation in50D, not once per color panel or repeated display; batch iLISI is reused between annotation variants. The five named broad/fine cohorts use author annotations. Lupus main cLISI measures **source-louvain-cluster preservation**, while its alternative measures preservation of source-provided named `cell_types`; neither interpretation is silently upgraded to proven original-author biology. The new resolution-2 Leiden labels are not cLISI truth.

5. **Q5 answered:** use **all cells for fitting and scoring**. A cap applies to display only, now settled by Q12. No display subset becomes the fitted or scored population.
6. **Q6 answered:** compute LISI in the corresponding **50-dimensional PC/Harmony representation only**, and annotate both PCA and UMAP figures with those results. Labels/captions must explicitly state50D; compute no UMAP2D LISI. The shared Harmony score is reused across its repeated comparison/color panels.
7. **Q7 answered:** display normalized scores with **iLISI mixing higher-is-better** and **cLISI preservation higher-is-better**. Retain raw and scaled per-cell values and summaries. Normalization uses label-count bounds, not within-dataset minimum/maximum scaling. Report arithmetic means, not the SciB median convention; Q11 now selects cellwise means for panels, while both mean summaries are saved.

**Round3 annotation and weighting answers:**

8. **Q8 answered:** broad main figures plus fine alternatives for Alzheimer, Breast, Covid, Joanito and Stephenson. Declared broad/fine counts are24/131,10/58,10/54,11/12 and18/50 respectively. Counts do not yet establish used levels/missingness.
9. **Q9 answered:** Lupus `louvain`24 MAIN and source-provided `cell_types`11 ALTERNATIVE. Main semantics are source-cluster preservation, not named-cell-type ground truth. `ct_cov` is not selected.
11. **Q11 answered:** cellwise arithmetic means annotate panels; equal-Sample arithmetic means are also saved. Preserve both raw/scaled summaries and valid-label coverage.

**Display selection is now settled:**

10. **Q10 answered through Q12:** save one deterministic selected cell-ID set for reuse across all representations/comparisons/colors/annotation variants. Full-cell fitting/scoring is unchanged.
12. **Q12 answered:** display **all cells when the cohort has at most200,000**; otherwise select **100,000 batch-balanced cells**, retaining rare batches in full and redistributing unused quota among the remaining batches, **without replacement**. Use actual original per-cell batch labels, not sample-majority recoding. Breast's reported114,579cells are all displayed; the other five cohorts each have100,000 display cells under the current confirmed counts.
13. **Q13 clarified:** the unused alternative meant **1,000 per sample**, not1million. It is not the selected display method.

**Q14 answered:** the user approved pilot-first implementation, with measured timing and Alzheimer memory checks before expansion. **Q15 answered:** twelve main figures/twelveUMAP fits initially; all extra variants deferred. The accepted scientific/storage/sampling contract is summarized at the top.

Rasterized point layers, readable legends and local PDF output are accepted. Used-label coverage and numerical validation are implementation checks, not additional approval flows.

### LISI implementation facts and proposed implementation path

Root's installed-source inspection establishes `harmonypy`2.0.0's native LISI path: exact KDTree queries for90nonself neighbors at perplexity30. Its per-label loop repeats neighbor search for the two label vectors, which may be expensive at2.93million cells/50D. The existing15-neighbor graphs do **not** provide the90-neighbor support required by that path. Do not run perplexity30 against those15-neighbor graphs or overwrite UMAP graphs to obtain LISI neighbors.

Reuse the cached15-neighbor graphs for UMAP. A separate ANN90 neighborhood computation **for LISI only**, shared where populations permit, is proposed. A source-pinned, tested existing scorer or minimal native/Numba scoring adapter is needed: this installed Python API cannot directly consume precomputed neighbor indices/distances; the runtime R `lisi` path is unverified, and no local `lisi` package is installed. No new packages, runtime changes or environment refresh are planned. Scorer/adapter selection is a routine implementation decision within final implementation approval, not a separate user approval flow.

The normalized directions use actual evaluable label-count bounds (batch categories for mixing; chosen source-cell-type categories for preservation), with degenerate single-label cases handled explicitly. Preserve raw values, scaled values, exact label counts, population/coverage and parameter provenance. Do not substitute benchmark normalization or dataset min/max scaling.

Missing-label handling remains a proposal: keep all cells in geometry and use gray for missing display labels; evaluate cLISI only where the chosen fixed source annotation is valid, report coverage, and build a separate evaluation neighborhood if required by that valid-label population. Do not silently introduce an `UNKNOWN` biological/cell-type class or fit geometry to phenotype labels. Numerical arrays/codes have not yet been read to establish actual missingness.

## Scientific and storage boundaries

Biology is coloring/ground truth only: it must not enter HVG selection, PCA, Harmony, neighbors, Leiden or UMAP. Selected source annotations, including the explicitly labeled Lupus source-cluster variant, may evaluate cLISI after representations are fixed, but must not train them. Per-cell batch colors/iLISI use actual cell metadata, **not** the sample-majority `e75...` reporting derivative. Constant Breast assay remains recorded, while the agreed Breast display/iLISI uses dissociation time. No Pipeline4 refresh is implied by these batch-only figures.

**Lupus whole-sample curation is intentional.** The approved sample-level selection retains `dmx_AbFlare-4` and the20/21/22pool groups recorded in the original plan (selection table/clarification around lines294/307), using the established majority policy where needed. Retaining whole samples can preserve minority **cell-level** batch labels outside those seven dominant groups. Root's existing UNC source observation already contains14labels in the inspected slice; that is not a full census. Neither “seven active cell batches” nor “all17declared levels are active” is established. The future full-cell job must count actually used original per-cell `batch_cov` labels for display/iLISI and normalization bounds, without an `e75` recode or pruning minority cells. This is not evidence of a curation-scope bug and triggers no current-method/Stage3 reprocessing.

PCA variance percentages require actual matching stored provenance; none is established here. Harmony dimensions must not be labeled with PCA variance explained. Raw and Harmony UMAPs use their respective graphs; independently fitted axes are not a physical alignment or point-displacement map.

RNA matrices and H5AD parents remain on authoritative HPC/NAS, immutable. New derived outputs use fresh roots and contain coordinates/selected annotations/provenance only; no parent overwrite or `X` copy. Local transfer is limited to explicitly approved small metadata/rendered outputs, not full H5ADs or whole expression/count matrices. Million-cell coordinate exports are not automatically small local metadata.

### Bounded local rendering proposal — not approved or implemented

After Q14/final implementation agreement, each representation worker could additionally emit one capped **display-only Feather**, containing only selected cell IDs, the explicitly chosen obs fields, PC1/PC2 and UMAP1/UMAP2, plus a small JSON of validated **full-population** metric summaries/provenance. Proposed hard budget: at most32MiB per Feather and the exact Q12-derived row count; JSON has a separate small bound fixed in the final contract. All full-cell50D PCs, UMAP arrays, raw/scaled per-cell LISI, NN90 caches and H5ADs remain on HPC/NAS.

Use a dedicated display-selection RNG independent of UMAP/LISI. The same canonical ordered cell IDs, seed and selected sampling rule must yield one recorded cohort display-ID set reused across all three representations and every panel/annotation variant. Check its ordered hash across the extracts; never independently resample after fitting. Exceeding a row/byte budget fails extraction rather than silently thinning or omitting fields; it does not reactivate the unused per-sample alternative.

This can satisfy the repository's small-metadata exception **only** when the final confirmed contract names the finite files, row/byte budgets and allowed columns. A local renderer would load only the three checked extracts and summary JSONs for one cohort, render from stored coordinates/validated full-population means, release them, then process the next cohort. It performs no fitting, scoring, neighbor construction, H5AD loading,50D-data transfer or mean recomputation from the display subset. A local `--render` mode or one small dedicated plotter is sufficient; existing benchmark notebook/helpers need no changes.

With that verified contract, downstream rendering can avoid six additional HPC render jobs. If the extracts cannot meet the agreed small-artifact boundary or coordinate/summary identity checks, retain HPC rendering instead; no silent scope or sampling change is allowed. Optional preservation of the original15-neighbor graph inside each expression-free HPC H5AD remains compatible with future Scanpy use; NN90 reuse/checkpoint products stay remote-only.

## Authorized availability check and future execution boundary

The user's question about data availability already authorizes read-only schema fact-finding; it does **not** require another permission request or completion of the remaining scientific decisions. Fixed-six UNC and COR header probes and the annotation-schema follow-up are complete as recorded above; annotation-source investigation remains availability work. Structural metadata includes parent proof/stat/source identities, column/key names, dtypes, declared category dictionaries, array shapes, graph/Leiden/UMAP key presence and available provenance. It reads no expression `X`, coordinate values, categorical codes, obs rows, or full-file hashes.

After shared understanding is recorded, future compute validation may check ordered cell IDs, actual coordinate finiteness, used annotation categories/missingness and numerical graph/embedding consistency. Those value-level reads are distinct from the header-only availability query and still exclude expression/count arrays. Pending Q14 blocks implementation, new UMAP computation, metric execution and derivative writes, not the authorized availability investigation.

Only after that schema review and explicit implementation/launch approval should an assigned implementer prepare new coordinate-only work. New compute must use the existing checked-in durable ECODA profile and a committed immutable snapshot containing the reviewed new source. The existing runtime image SHA `8fcd00b7a02592f82e8ad96ec9dffdcfb40214b3c47218f770dfb8e828e05e34` is proposed for reuse; no package/environment change is planned. Resources, UMAP parameters and the tested LISI adapter still need specification.

### Execution seam — pilot first; bulk scheduling deferred

Implement the smallest fixed-purpose scalar pilot first: Breast UNC,8CPUs/64G/2h on any CPU model, one task performing UMAP then ANN90/LISI then one-row diagnostic rendering, with separate phase wall-time/RSS. A separate metadata-only Alzheimer header preflight estimates selected-array/working memory;128G is a provisional conservative floor, not a measured peak or authorization to launch Alzheimer. Reuse existing scalar/native audit/finite publication patterns. Do not implement a bulk array or general variant scheduler before the pilot's measurements justify it. The earlier array discussion below records a possible later extension only; if adopted, its current selection is twelve representations, not eighteen.

**Prefer one18-row representation array:** six cohorts × `unc_raw`, `cor_pre`, `cor_harmony`, each task doing UMAP then NN90/LISI sequentially. This makes the three long Alzheimer representations independently schedulable while retaining one native manifest and watchdog. It preserves18unique UMAP fits. A36-operation UMAP/LISI split is not proposed: the existing code supplies no simpler phase-dependency/retry mapping, and it would add scheduling state.

Initial scheduler shape would be one ARRAY root with18tasks and one compute-node watchdog with `afterany` on that array. The controller waits for the watchdog, which owns any bounded OOM retry arrays, then performs final validation and finite NAS publication. Record all actual IDs immediately; no separate aggregate job is needed for this proposed single-array scope.

- **Small source footprint, if assigned:** one new `src/utils/bash/singlecell_batch_geometry.sh` with controller/worker/watchdog/read-only `--validate` modes, plus `src/notebooks/run_singlecell_batch_geometry.py` for the data operation. Reuse `ecoda_run_common.sh` atomic manifests/checksums/PUBLISHED records/owners and compute-node accounting helpers, `ecoda_runtime.sh`, and the immutable source executor. The wrapper follows the source-root/spool bootstrap and frozen-scope patterns of `project_analysis_release.sh`/`zeroimp_count_screen.sh`; it does not invoke benchmark dispatch or pretend these are method fits.
- **Existing watchdogs are patterns, not drop-ins.** `matrix_watchdog.sh` is coupled to Stage5 method-family ownership/runtime contracts; the Stage3 watchdog hardcodes the preprocessing worker and Stage3 output contract. A narrow purpose-specific watchdog mode must reuse their OOM-only classification, bound memory escalation, exact retry-subset manifest/readback and complete ID recording, without importing those unrelated output validators. Non-OOM failures/exhausted ceilings fail closed. Agent-side polling remains prohibited.
- **One final H5AD per representation:18 disjoint outputs.** Proposed scratch namespace is `_ecoda_runs/<new-run>/artifacts/singlecell_batch_geometry/<dataset>/<representation>/`, with separate run logs and a fresh NAS release prefix. No benchmark `ANALYSIS_ROOT`, benchmark sync destination, parent H5AD, or existing method file is written. Store selected obs/50D PCs/UMAP/Leiden/metrics/provenance only; no `X`, counts or raw matrix. After data review, prefer the bounded local-display-extract route above if its confirmed contract passes; otherwise six HPC cohort-render tasks can each read their three derived files. Both routes reuse cell IDs/colors/anchor labels. A merged cohort H5AD is unnecessary unless requested; downstream rendering remains separately assigned.
- **Bounded phase checkpoints:** atomically publish validated UMAP and NN90/LISI phase products with source/runtime/input/representation/ordered-ID/graph/parameter/label-policy fingerprints and native proofs/records. Within the same run, an OOM retry may reuse only matching current-run PUBLISHED/OK checkpoints and compute the incomplete phase. Final embedding H5AD is a fresh assembly after both phases validate. Failure cleanup marks incomplete ACTIVE phase owners, preserving already-OK checkpoints; the generic tracked-owner finalizer would otherwise mark every tracked owner FAIL. No generic failed/foreign-input ownership guard is relaxed, and no cross-run checkpoint adoption is implied.
- **Narrow generic audit addition is necessary.** `_audit_analysis_release()` currently requires exactly one scalar scheduler-ledger row (`ecoda_run_audit.sh:581–587`). Proposed exact PURPOSE `singlecell_batch_geometry` / operation `singlecell-batch-geometry` must have a separate array-ledger check, then reuse common source/runtime/terminal/selection/accounting checks and delegate to source-pinned `singlecell_batch_geometry.sh --validate`. Scope/submission producers can be `singlecell_batch_scope` / `singlecell_batch_submission`. Preserve all existing scalar cases and unknown-purpose rejection. Validate the exact18 selection and initial array +watchdog +every actual retry ID; no new general dispatcher or audit framework.
- **Accounting remains truthful.** Compute-node watchdogs may manage bounded OOM retries; the durable sole terminal inspection still retains every emitted attempt. The current profile requires all captured rows `COMPLETED|0:0`, so historical OOM rows may leave the original gate failed even when selected final products validate. Preserve that outcome; do not promise successful global approval or hide attempt IDs. Any later artifact-only release would require its actual reviewed evidence, not a fabricated green gate.

After complete selected-output validation, the existing `ecoda_nas_backup.sh` finite inventory/copy/verify path is a publication candidate for the new NAS prefix. Its source-manifest/checksum/mount/selected666-file/777-directory checks stay intact; it performs no scratch chmod. Exact final/checkpoint/evidence inventory must be specified during implementation, with no directory-wide discovery or invented backup receipt.

Source-level feasibility of expression-free UMAP is established in installed AnnData/Scanpy: AnnData accepts `X=None` with shape inferred from obs/obsm, and Scanpy's UMAP consumes the selected stored graph plus the explicit `obsm` representation named by neighbors `use_rep`/`n_pcs`. The implementation must supply that explicit50D key and cached15-neighbor graph, avoiding an `X`/automatic-PCA fallback. This is source inspection, not an executed runtime smoke test or authorization to read coordinate arrays now.

**Serialization proposal: queue, not bypass.** The checked-in durable profile accepts an opaque group but defines no dedicated evidence group. Existing `ecoda-annotation` and scoped preprocessing groups represent other real stage scopes; their names do not establish permission for this new derived Stage5 purpose. Use the existing `ecoda-benchmark` group and queue after Gong and the33-projection terminal review. Private artifact/NAS roots do not alone justify relabeling the gate to evade the active mutex. No user operational exception is needed for this queued path.

Scheduling must respect the existing serialization/sync-owner contract and the ongoing Gong terminal wait/inspect/review. Resolve that ordering with the Supervisor; create no alternate group or bypass. Exact selection, fresh run/NAS roots, one durable waiter, sole complete terminal inspection, native validation, selected NAS checksums/modes and reviewer approval remain required. This draft supplies no runnable command and authorizes no jobs or data writes.
