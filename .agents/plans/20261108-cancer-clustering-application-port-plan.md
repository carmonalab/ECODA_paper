# Cancer clustering application port

Status: **IMPLEMENTED / VALIDATED; available for user review**. User confirmed the design and authorized implementation on 2026-11-08 (Europe/Zurich). No commits or pushes. Keep this progress record at the requested active path through user review.

## Authority and boundaries

- **DO NOT COMMIT OR PUSH.** This applies to every part of this task. Another agent owns the concurrent confounding/reprocessing plan and pipeline changes; preserve all unrelated working-tree changes. That agent may independently include this work in its own commits.
- This task owns the new application notebook, its focused helpers/tests, compact application inputs, application figures, Yan NAS archive, and this plan. Do not edit the other agent's plan, notebook guidance, pipeline files, datasets.json, pixi.toml, pixi.lock, or vendored QOT source.
- No cohort preprocessing, annotation reruns, benchmark reruns, scheduler submissions, environment refresh, or production artifact replacement. Full expression objects remain on NAS/HPC. Read selected metadata only, cohort by cohort.

## Agreed implementation

Create `notebooks/ecoda_application_cancer_clustering.rmd`, using the current `/Users/christianhalter/Desktop/analysis_Yan_2026_42129561/Analysis_unsup_strat.rmd` as code reference and historical PDFs as visual references. Put focused reusable functions in `src/notebooks/ecoda_application_cancer_clustering.R`.

Outputs:

1. Complete breast composite: Yan/Zhang/Bassez HR immune PCA, response boxplots, and breast consensus PC1/PC2 stacked loadings.
2. Bassez shared-cell clustered heatmap.
3. Yan/Zhang shared-cell heatmaps and their appendix composite.
4. Pelka shared-cell heatmap.
5. Four-cohort clustered PC1/PC2 loading-correlation heatmap.

Use descriptive filenames in run-owned application directories under data/ and plots/. Do not renumber publication figures. Omit unrelated LR, prediction, projection, mutation-track and meta-analysis work.

### Inputs and archive

- Current validated benchmark-analysis Bassez/Zhang/Pelka inputs, resolved through dataset/view helpers and reconciled with `.agents/plans/20261108-confounding-subsets-benchmark-reprocessing-plan.md`. Unchanged HiTME outputs may be reused with dependency evidence; do not wait for unrelated methods.
- Scope compact export to those datasets. Use verified HiTME layer-2 composition artifacts and obs-only sample metadata. Recover raw counts only after confirming pseudocount convention, nonnegative integer counts and matching IDs. Reject invalid provenance or mismatch; no automatic repair compute.
- Additional application restrictions: Zhang Tissue == tumor; Pelka Tissue == T; preserve configured upstream cohort restrictions.
- Import existing Yan ecoda_data.rds (33,676 bytes; observed MD5 35cc7ebec2a54ed17052756d3a65c42c): 99 samples, layer 1–3 HiTME counts and metadata.
- Archive existing Yan annotated per-sample Seurat outputs, composition summaries and relevant legacy provenance directly within NAS into a new ECODA application-owned destination. Verify checksums; normalize only selected destination files to 666/directories to 777. No full H5AD/Seurat transfer to workstation.
- Yan remains application-only, outside datasets.json. Document legacy annotation and missing historical reference/runtime identities without asserting current Pipeline 4 compliance.

### Scientific contract

- Copy exact source immune exclusions/regex, including neutrophils. Retain samples with >=10 eligible immune cells before sparsity/intersection. Drop types with zero fraction >0.9 independently per cohort; exactly 0.9 stays.
- Intersect retained cell types across Yan, Zhang, Bassez AND Pelka; no second post-intersection sample cutoff.
- Add 0.5 to EVERY retained count, then natural-log CLR. Use one common CLR calculation for PCA and heatmaps. This intentionally changes source PCA (+1).
- Independent centered, unscaled PCAs; preserve lineage-based PC1 orientation and PC2 alignment to Yan. Biological labels are annotation/evaluation only.
- Breast consensus bars average Yan/Zhang/Bassez only.
- Preserve response mappings: Yan pCR=High, RD=Medium, literal NA=Low; actual missing stays missing. **User explicitly confirmed the NA group comprises early dropouts with worst/Low response.** Supplementary table itself only says NA means not available; distinguish user interpretation from table evidence.
- Heatmaps: feature-centered CLR, no variance scaling, Euclidean distances, Ward.D2. Cuts: Bassez 3/3, Yan 4/3, Zhang 3/3, Pelka 2/2. Pelka annotates MMRStatus; others Response.
- Correlation plot compares eight loading vectors (four cohorts x PC1/PC2), Pearson correlations; Ward.D2 on 1-correlation, lower triangle, original raw unadjusted Pearson significance stars.

### Runtime and traceability

- Use installed Pixi R 4.5.2/Python 3.13 and existing packages; no external sibling-library loading, old library paths or namespace monkey patches.
- Notebook consumes frozen compact inputs; import/archive is a separate explicit preparation operation.
- Atomic writes/checksums; retain input identities, source hashes, package versions, retained/excluded sample/type lists and reasons, shared ordering, CLR matrices, scores/loadings/variance, correlations and cluster assignments.
- One persistent R session for figure iteration. Review PDF legibility, clipping, legend completeness and appropriate numeric formatting.

## Validation and acceptance

- One focused regression script covering cutoff boundaries, four-way intersection, CLR consistency, invalid counts/IDs and independence from response labels.
- Compare scoped source calculations on identical current inputs with +0.5 against port at tolerance 1e-8 with consistent signs. Historical PDFs are visual references, not numerical equality requirements.
- Execute relevant notebook chunks; inspect every requested PDF. Verify archive/source checksums and permission scope. Record unavailable/invalid input dependencies without launching jobs.
- All requested code, compact inputs, archive evidence and figures available for review; no commits or pushes.

## Progress

| Work package | Status | Evidence / next action |
|---|---|---|
| Save approved plan | DONE | This file created before implementation |
| Reconcile inputs/provenance | DONE | Current obs aggregates frozen with sidecars; stale legacy bundles retained only for comparison |
| Notebook/helper implementation | DONE | Notebook, analysis/plot helper and explicit preparation script written; all notebook chunks execute successfully |
| Yan NAS archive | DONE | 216 selected source files / 4,297,021,736 bytes; all hashes verified; final 220 files/7 directories, Bamboo mode audit PASS |
| Regression/numerical validation | DONE | Focused R regression PASS; source numerical parity max2.176037e-14; final read-only review accepted |
| Figure execution/visual QA | DONE | Seven PDFs rendered/inspected; source defaults and footer layout corrected; persistent R session22115 retained for iteration |

## Findings and execution log

- 2026-11-08 planning: source Rmd currently emits figure_2.3/2.5/2.8, whereas some requested examples use older names. Current Rmd is the implementation reference.
- Source Bassez/Zhang/Pelka inputs came from sibling scECODA example_data.rda, not directly from current canonical H5ADs. Current inputs plus +0.5 can change historical numbers.
- Source four-way intersection includes Pelka even in breast plots; preserved by user decision. Source >=50 comment contradicts actual >=10 constant; preserve actual behavior and correct documentation.
- Legacy Yan uses author H5AD X, STACAS gene-symbol standardization and HiTME with tumor scGate models; reference filenames match ECODA but historical map checksums are absent. Legacy scGate dev DB was force-refreshed. No current scATOMIC claim.
- Planning runtime check: installed scECODA 1.1.6 RemoteSha 055bc5efb51cac7c3b6e0ea7eae84420f0a8d21a; ggrepel, patchwork, ggpubr, cowplot, pheatmap and readxl available. No package/environment changes needed.
- Implementation start: shared working tree has unrelated pipeline/test/plan edits. They remain outside this task's ownership. Current concurrent plan records Stage 3 READY_TO_PREPARE and no launch; recheck relevant dependency evidence before selecting inputs.
- Bamboo read-only preflight: HOME=/home/users/h/halterc, host login1.bamboo, scratch resolves /srv/beegfs/scratch/users/h/halterc/ECODA_paper, configured Python 3.13.14. Selected scratch/NAS processed files exist with checksum sidecars and automatic-annotation marker. No matrices read. Sizes Bassez 3,301,058,992; Zhang 6,657,128,030; Pelka 18,636,721,109 bytes.
- Existing benchmark/checksums.md5 records HiTME layer-2 RDS identities: Bassez 6d80bf0ebe957350b174132bf232997e; Zhang 063b4946da8c2f32495f1c79c81f0173; Pelka 2930f77463b61d975a62d11ba21fae18. Reuse remains subject to obs-count reconciliation.
- Implementation addition: a separate explicit preparation entry point, src/notebooks/prepare_cancer_clustering_inputs.py, freezes only the three selected compact inputs. It sends a read-only obs query through the configured Bamboo interpreter using SSH stdin; it installs no cluster code and submits no work. This keeps the notebook offline/reproducible and avoids executing the existing all-dataset exporter.
- Cross-chat coordination: user identified session 01a1bc2a-e71d-77e3-a838-e32048cd72cc. Its orchestrator independently confirmed ownership separation, will exclude our changes from its commits, and permits validated existing HiTME reuse without claiming new campaign provenance. No request to interrupt or change its work was made.
- Archive finding: macOS SMB reports mode 700 even after destination-only chmod. Keep archive unpublished until required destination modes can be checked from Bamboo's authoritative NAS mount; source permissions remain untouched.
- **Input reconciliation finding and scoped adjustment:** legacy HiTME bundles pass recorded MD5 but do not match current obs counts after subtracting 0.5. Sample sets agree. Bassez current obs adds one AS-DC cell (BIOKEY_19_Pre), absent from bundle; maximum common-type count differences are Bassez 25, Zhang 144, Pelka 44. These bundles are not eligible for current-count reuse. Retain them only as comparison evidence. Use fresh raw composition aggregates from the already-collected, checksummed current obs snapshots instead; no preprocessing, annotation, or benchmark fit is required. This follows the user's current-input decision and changes no production artifact. Record exact comparison tables and source identities with outputs. Other orchestrator informed with no action requested.
- Initial importer incident: H5AD sidecars use repository MD5=/SIZE=/PATH= fields, not md5sum syntax. Corrected importer to call the existing validated sidecar parser. Scratch/NAS MD5/size/PATH identities agree. Full expression bytes are not rehashed by this obs-only audit; record the metadata scope explicitly.
- Current aggregate inputs: Bassez 26 samples/42 types; Zhang 31 samples/41 types (before application tumor restriction); Pelka 98 samples/46 types (before tumor restriction); Yan 99 samples. Application retains Yan 99, Zhang 13, Bassez 26, Pelka 63, with 34 common immune types.
- Focused R regression PASS: >=10 boundary, exactly0.9 sparsity boundary, Pelka-dependent intersection, +0.5 CLR/PCA reconstruction, heatmap same-CLR centering, response-label permutation independence, invalid counts and IDs. Command: pixi run -e default Rscript tests/test_cancer_clustering_application.R.
- Yan archive DONE: /Volumes/Shared/Projects/ECODA_paper/applications/cancer_clustering/imports/Yan/20261108-v1. Manifest MD5 6683d38b61426d850557afb44bc09319. Includes 99 annotated sample RDS, 99 summaries, aggregate RDS, 15 historical code/reference files and two spreadsheets. Bamboo final audit verified 220 files666 / seven directories777, owner halterc:hpc_users. Local import RDS hash35cc7ebec2a54ed17052756d3a65c42c; provenance hash17da971efb380fe6d9e7d8ae0f990546. Only exact interrupted unpublished Mac staging removed after successful replacement; no source permissions changed or staging left.
- Source numerical parity PASS: selectively parsed six original pure helpers from frozen archived functions.R (never evaluated its top-level namespace patches), changed only original CLR +1 to +0.5, and compared current-input raw counts/shared order/CLR, PC1-PC2 scores and loadings, full explained-variance spectrum and breast consensus means/SEM. Maximum difference 2.176037e-14 <1e-8. Evidence and reproducible validation script saved under derived/source_parity.rds and derived/source_parity_validation.R with MD5 sidecars.
- Visual QA: rendered and inspected all seven standalone PDFs. Fixed heatmap footer clipping by reserving Samples-label height within the original canvas rather than extending 1npc. Restored source ggpubr Wilcoxon labeling defaults after its asymptotic significance screen. Deterministic horizontal-only point jitter is a display adaptation. Complete figures are under plots/ecoda_application_cancer_clustering/20261108-v1/.
- Notebook execution incident: all analysis/figure chunks completed, but knitr's per-chunk working-directory restoration made relative source hashes NA. Fixed with opts_knit(root.dir=project_root), absolute source paths, and an explicit no-NA hash assertion. Preserved the initial unreleased invalid identity RDS/sidecar under derived/quarantine/initial-rootdir/; scientific matrices were valid and unchanged. Corrected complete notebook execution PASS without the working-directory warning; existing numerical RDS evidence was verified/reused identically.

## Final evidence and findings

All paths below are rooted at /Users/christianhalter/Desktop/ECODA_paper unless explicitly NAS/Bamboo paths.

- Notebook: notebooks/ecoda_application_cancer_clustering.rmd (MD5 9865a40c51776f30f6068ee63ca87746).
- Helper: src/notebooks/ecoda_application_cancer_clustering.R (MD5 293ae4fc814eedf5f39362bba46b6b23).
- Explicit no-compute importer: src/notebooks/prepare_cancer_clustering_inputs.py (MD5 40e7297de0de5e9d71a1e3d5ad07f6a8).
- Focused regression: tests/test_cancer_clustering_application.R; passed after final helper changes. Python importer syntax compilation passed; live selected metadata extraction completed. No broad test/pipeline run.
- Complete notebook execution: knitr::knit(input='notebooks/ecoda_application_cancer_clustering.rmd', output='/tmp/ecoda_cancer_application.md', envir=.GlobalEnv), in persistent Pixi R4.5.2 session22115. All chunks execute, including evidence idempotency. The session remains open for user-directed figure iteration; do not close without user confirmation per notebooks/AGENTS.md.
- Compact input root: data/applications/cancer_clustering/20261108-v1/. Each Bassez/Zhang/Pelka directory contains a checksummed obs_snapshot.json with exact configured/canonical sample mapping, selected cell-type counts, clinical metadata, source sidecar identities, selected-obs SHA256 and runtime/decoder identities. Each also contains the checksummed legacy HiTME bundle strictly as comparison evidence. Yan contains its verified RDS and archive provenance.

| Dataset | Current H5AD recorded MD5 | Legacy HiTME bundle MD5 | Current samples before/after application selection |
|---|---|---|---|
| Bassez | 22be4eb1f41433b4a8a33dd8e1a3b7d2 | 6d80bf0ebe957350b174132bf232997e | 26 / 26 |
| Zhang | 7c0854d54bbfdc6b8eb971def8aa7f84 | 063b4946da8c2f32495f1c79c81f0173 | 31 / 13 |
| Pelka | f3f2fc8fff93d32a805ccef15d67041b | 2930f77463b61d975a62d11ba21fae18 | 98 / 63 |

- H5AD paths are captured exactly in each obs_snapshot.json provenance.source and provenance.nas_source. Scratch/NAS existing MD5/size/PATH sidecars agree. This is explicitly an obs-only audit, not a fresh full-expression checksum validation. Legacy/current count discrepancies were not explained by a documented historical filtering convention in the inspected code; no equivalence claim is made.
- Scientific evidence: data/applications/cancer_clustering/20261108-v1/derived/{input_and_source_identity.rds,analysis.rds,legacy_bundle_reconciliation.rds,figure_numerical_evidence.rds,session_info.rds,source_parity.rds,source_parity_validation.R}, all with MD5 sidecars. Final audit verifies every selected evidence RDS and figure sidecar, and the recorded source hashes equal current notebook/helper/importer files with no missing values.
- Final composition: 34 common types; Yan99/Zhang13/Bassez26/Pelka63 samples. Yan labels retained as user confirmed:44 High/37 Medium/18 Low. PC1/PC2 percent variance (two significant digits): Yan27/14, Zhang30/22, Bassez27/16, Pelka18/16.
- Figure root: plots/ecoda_application_cancer_clustering/20261108-v1/. Complete breast composite, four individual cohort heatmaps, Yan/Zhang appendix composite and four-cohort loading correlation PDF; all have MD5 sidecars. Visual QA accepted final breast composite after restoring source Wilcoxon labeling defaults, all heatmaps after Samples-footer fix, and correlation plot. Notebook preview PNGs are under the same root/notebook/.
- Independent bounded read-only review of notebook/helper/importer/plan found no material correctness issue. Review confirmed scientific order, current-obs source, shared CLR, label independence and provenance limitations.
- NAS archive ownership/permission/source checksum audit is complete as recorded above. No full H5AD/Seurat file or expression matrix was copied locally. No canonical data/result artifact, dataset configuration, pinned environment, other agent file or publication figure was modified by this task. No commits, pushes, scheduler jobs, environment refreshes or cluster-checkout edits were performed.
