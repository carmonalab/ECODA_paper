# Confounding subsets, benchmark repairs, and Bamboo reprocessing

Status: **IMPLEMENTING. User confirmed shared understanding and authorized implementation, checkpoint/push, documentation harmonization and parallel work.**
Created: 2026-11-08, Europe/Zurich. Owner: the orchestrator of this chat.

## 1. Roles, authority, and update contract

**You are an orchestrator, implementor, or reviewer. Your assignment must name
your role and scope. If either is unclear, ask your supervisor before continuing.**

- **Orchestrator:** maintains this plan, resolves decisions with the user,
  assigns bounded work, records dependencies/selections and accepts evidence.
  Does not implement pipeline changes. Authoring this plan is orchestration.
- **Implementor:** changes only assigned files/behaviors, runs focused checks,
  and reports commits, artifacts, limitations, and proposed status updates.
- **Reviewer:** independently checks the assigned diff or run evidence against
  this plan. Reports actionable findings and explicit acceptance; does not
  expand the compute selection. The operating agent may perform required
  durable-gate review; a separate reviewer is useful for implementation work.
- The orchestrator alone changes authoritative decisions and the launch ledger.
  Other agents return updates; they do not concurrently rewrite the plan.
- **Concurrent external agents:** the user is assigning separate agents to
  `notebooks/` and other files outside this task, mainly downstream analysis.
  The user will give them this plan and instruct them not to commit or push.
  Their working-tree changes are expected, are not evidence of corruption,
  and remain owned by those assignments. Preserve them and exclude them from
  this task's commits/pushes unless explicit shared ownership is agreed.
  Before touching an overlapping file, coordinate scope rather than replacing,
  stashing, reverting, or indiscriminately staging their work. Use exact file/
  hunk staging and verify ownership immediately before each commit.
- P6 notebook/plot implementation must coordinate with those external agents;
  the scientific reporting requirements here remain the agreed contract, but
  downstream execution is not automatically claimed by this orchestrator.
  The current P1 documentation assignment includes only `notebooks/AGENTS.md`
  and `notebooks/dataset_onboarding/README.md` inside the notebook tree; flag
  any overlap there explicitly. Notebook source/helper edits outside those
  assigned documentation files are not part of the current P1 package.
- Scientific workers execute the exact committed immutable snapshot. Never
  pull uncommitted notebook changes into run sources, refresh a running
  snapshot, or broaden a commit to make the shared working tree appear clean.
- Keep each work package marked `PLANNED`, `IN_PROGRESS`, `BLOCKED`, `REVIEW`, or
  `DONE`, with commit/evidence links and an explicit completion criterion.
  Configuration recorded, derivative generated, validated, synchronized, and
  reviewer-approved are distinct milestones.

This plan governs the work agreed in this chat. Existing run evidence and
immutable manifests retain their original meaning. Link older overlapping
plans to this checkpoint during documentation harmonization; preserve their
historical records and unrelated proposals.

Authorization: the user explicitly confirmed shared understanding and asked
to start implementation with logical checkpoint commits pushed to master,
then documentation harmonization and parallel independent work. P0 is now
complete. Later compute remains subject to the recorded integrity, immutable
source/runtime, exact selection, and durable-gate prerequisites. Scheduler
submissions and their outcomes are recorded in the live checkpoint and launch ledger below.

Provenance: repository planning anchor
`5302671ad94556edcf9acccf372d2dc34121d714`; inspected HEAD
`70b6033`. Record subsequent checkpoint, documentation, implementation, and
execution commits separately. A dirty configuration is not proof that its
decisions have been applied to existing artifacts.

## Live checkpoint — 2026-11-09 13:20Z, live MrVI progress confirmed; recovery pending terminal evidence

- P0/P1 DONE: logical checkpoints and harmonized documentation are pushed.
  Active controller/reuse source is9ca2453, committed/pushed and
  deployed as an immutable Bamboo snapshot. Later reviewed diagnostics6bcf26b
  and compact-analysis code14e6749 are pushed but do not alter active workers. External application work remains excluded.
- **Main wave progress confirmed at13:20Z:** source9ca, launched08:18:46Z;32method
  arrays and32watchdogs submitted across benchmark/uncorrected/corrected.
  Existing Positron SSH terminal confirms uncorrected MrVI task4665560_3
  around epoch44/50 and corrected task4665623_3 around24/50, with recent
  log writes. User queue at13:09Z shows both RUNNING and benchmark GPU
  tasks4665568_1..2 waiting for unavailable/reserved hardware. No context
  terminal marker or terminal accounting inspection has been observed.
- **Hardware policy clarified/applied:** benchmarktimings retainEPYC-7742/
  configuredGPU model; bothbatch views permit anyhardwaremodel. Exact18
  pendingbatch tasks had onlySlurmFeatures cleared with before/afterproof;
  zero benchmark/runningjobs changed, zero jobs submitted/cancelled.
- **Failures being catalogued:** priorStephensonGloScope parameter4665516_2
  and now8correctedQOT tasks4665627_1..8. QOT has no stdout/stderr files;
  task1 failed in1sec oncpu046 with0:53. StartupI/O is a supported hypothesis,
  not a provedcause or QOTmodel failure. Preservevalidwork; recovery scope
  must come from complete terminalevidence, not a broad rerun.
- **Changed inputs ready:** approved Bassez metadata repair (74948cells,
  26samples,40HRtypes); all17 Stage3 outputs produced/NAS-synced;53-sample
  Stephenson benchmark input individually NOOP-reviewed. Original Stage3 gate
  remains FAILED with its historical OOM evidence preserved.
- **All16 batch inputs reviewed for reuse:** exact expected counts in both
  views101/33/44/124/80/94/80/110.52 batch tables pass except the documented
  two Ncl exceptions. Receipt SHAc8e1d9b5...; no old failed-gate approval.
- **Kidney two-axis metadata validated:** both views have the same80 samples,
  condition21AKI/26CKD/33Ref and tissue45kidney/23cortex/12medulla. Tissue is
  fully confounded with cell/nucleus suspension and remains descriptive.
  No extra H5AD, feature fit or cohort filtering is needed for this evaluation.
- **Stephenson Stage4 COMPLETED and formally approved:** wrapper0 at04:50:11Z,
  all32 raw scheduler rows0 across6roots; accounting bridge and pinned
  protected-copy/annotation/source/runtime/NAS audit pass. New input MD5
  `7d92fb6e0d3d7dd5c3500212fd951493`,size13,357,978,168.
- **BL10 DONE and published/reviewed:**20RDS+10timings copied into a fresh
  validated NAS release, source3794e3e, publication finished05:03:40Z; no fit
  repeated. Original failed BL gate remains failed. Stephenson BL stays in
  its current whole-method run, once only; worker watchdogOK, final
  source/artifact/terminal review pending. Aux remains unchanged.
- **Bassez authorHR22 science and publication approved:** original gate
  failed only in controller ownership validation after all5jobs completed0.
  Its23RDS+one22-row timing file pass exact CLI,5finite metrics, current
  metadata/NULL/count/source/proof checks with no fit/H5AD read. Original
  abort cleanup already finalized ownersOK; no old record was edited.
  New science receipt SHA78523c4...; new24-payload publication completed
  05:48:23Z and passed terminal validation/root review without refitting.
- **Historical mainv1 setup failure preserved:** sourcec534cdf,06:17–06:46Z,
  all6actual scheduler roots/30rawrowsCOMPLETED0. DuplicateGloScope family
  ownership and scheduler-proof controller issues were repaired/reviewed/
  committed; source9ca finite23prepared copies reused validated outputs in
  currentv2. The oldfailed gate was not approved or relabelled. No BL10,
  authorHR22 or completed preprocessing was recomputed.
- **User HiTME clarification: NO targeted refresh added.** Plain layer2
  counts match in all3audited cohorts and maintained final selectors use that
  plain key. Disputed variance entries are unused in maintained figures;
  their historical experimental provenance is unproven. The earlier100-row
  refresh proposal is withdrawn from this implementation scope. The user
  will decide any further action after clarification; do not infer approval.
- P6 five-metric/Kidney reporting implementation is committed. Final figures
  and remaining method/screen waves are not claimed complete or released.

## 2. Confirmed scientific decisions

### Participation: flags and views must both permit a production row

View existence alone never enables a dataset. Resolve participation through
`use_for_benchmark` / `use_for_batch_effect`, then the relevant configured
view. Underscore-prefixed diagnostic datasets are excluded from production.

| Participation | Dataset IDs |
|---|---|
| Benchmark only, 10 | `Adams`, `Bassez`, `Gongsharma_cmv_young_males`, `Kfoury`, `Kim`, `Lee`, `Pelka`, `Smillie`, `Wu`, `Zhang` |
| Benchmark and both batch views, 1 | `Stephenson` |
| Both batch views only, 7 | `Alzheimer`, `Breast_cancer`, `Covid19_PBMC`, `Joanito`, `Kidney_KPMP_full`, `Lupus_PBMC`, `Lung` |
| Disabled | `Zhu`, `CombinedPBMC`, `Kidney_KPMP`, `Myocardial_infarction`, `Diabetes`, `Parkinson` |
| Diagnostic only | `_debug` |

Thus benchmark blacklist repair covers **11 production datasets**; batch
preprocessing covers **8 datasets x 2 views = 16 rows**. Alzheimer, Joanito,
and Kidney full do not receive benchmark runs. Stephenson is the sole
production overlap. Add this distinction to `AGENTS.md` and operator docs.

### Fixed cohort selections and balance validation

The named selections below are authoritative decisions, not examples for an
automatic pruning algorithm. Validate them and return discrepancies to the
user; do not silently derive different cohorts.

- Apply approved exclusions and existing compatible view restrictions first.
- For each retained batch, require at least two represented biological groups.
  Each represented group needs at least two samples and at least **20% for two
  represented groups**, or **10% for three represented groups**. Equality passes.
- Absent global groups remain reported structural zeros; they do not by
  themselves invalidate a batch. A batch with only one represented group fails.
  Unexpected four-or-more-group cases require a decision; no threshold is invented.
- COVID additionally requires healthy controls in each retained batch; d07 is
  explicitly excluded. Missing disease categories are allowed.
- **Approved exception:** Stephenson Ncl remains after quality exclusion of
  two Healthy outliers: 43 Covid / 10 Healthy, Healthy **18.87%**. Record the
  exception explicitly; do not downsample Covid to manufacture 20% balance.
- Balance eligibility is assessed against the configured primary biological
  label. Kidney's extra tissue evaluation is descriptive and does not trigger
  another cohort filter.

| Dataset | Approved batch-analysis selection | Expected samples per view |
|---|---|---:|
| Alzheimer | Globally exclude `Cognitive status == "Reference"` (3 samples); retain Dementia / No dementia across assay batches; no sex correction | 101 |
| Breast_cancer | Retain only `10x 3' v3` with dissociation `12 hour` or `15 minute` | 33 |
| Covid19_PBMC | Retain only d02 (22), d09 (12), d13 (10), AND the existing `control` or finite numeric day <=30 rule; exclude d07 and all other batches | 44 |
| Joanito | Globally exclude `Tumor-2`, retain the existing LymphNode exclusion | 124 |
| Kidney_KPMP_full | Exclude Papilla; correct by `suspension_type` only, never sex | 80 |
| Lupus_PBMC | `dmx_AbFlare-4` (11), `dmx_YS-JY-20_pool3/4`, `dmx_YS-JY-21_pool1/2`, `dmx_YS-JY-22_pool5/6` (83 combined) | 94 |
| Lung | Six exact batches listed below; preserve existing 10x/lung/COPD restrictions. Normal disappears with the excluded batches | 80 |
| Stephenson | Exclude Sanger, plus all existing sample exclusions including the two new outliers | 110 |

Lung's exact retained batch values are:

1. `Goveia_Carmeliet_2020` (13).
2. `Lambrechts_Thienpont_2018_6149v2` (12).
3. `Lambrechts_Thienpont_2018_6653` (11).
4. `Leader_Merad_2021_10x_3p_v2_beads_cite` (18).
5. `Leader_Merad_2021_10x_3p_v2_digest-deadcell_cite` (22; ASCII hyphen in metadata).
6. `Leader_Merad_2021_10x_5p_v1_CD2` (4).

The original Lupus note omitted batch 22; the user confirmed 20+21+22 is
correct. Breast v2 was an exclusion example; the retained chemistry is v3.
Diabetes remains entirely disabled, rather than retaining its old subset.

These counts are now confirmed by the sixteen-view Bamboo NOOP receipt and
independent review recorded below. Raw/current sample masks and >=500-cell
retention agree in both views; no selection discrepancy was found. Future
unexpected discrepancies must return to the user, never change the selection
automatically.

### Other benchmark changes

- **Stephenson:** exclude `MH8919226` and `MH8919227` from benchmark and both
  batch views; preserve previous exclusions `BGCV10_CV0198` and `MH8919230`.
  Both new exclusions are Healthy/Ncl samples with no annotated T cells.
  Benchmark remains Ncl-only, expected 53 samples. Recompute all ten baseline
  method families with maintained eligible variants, plus `trans` and
  `zeroimp` screens. Sample/rank eligibility still applies.
- **Bassez:** valid `cellSubType` remains HR; otherwise use configured broad
  `cellType`. Preserve valid labels and all cells. Existing implementation
  exists, but the current published NAS processed H5AD still lacks the complete
  fallback (verified metadata-only on 2026-11-08): 74,948 cells, 24,139 valid
  HR cells, 50,809 missing HR cells, 26 samples. Valid HR has 33 categories;
  the proposed mixed HR/LR labels have 40, with zero invalid broad fallbacks.
  This independently reproduces the October diagnostic. Bamboo's current
  producer remains to be reconciled before repair/compute.
  **User checkpoint resolved:** after the NAS audit and method-impact report,
  the user selected a validated HR-only metadata derivative and targeted
  dependent reruns. Verify Bamboo's producer and unchanged inputs first.
  A broader rebuild or expanded scientific invalidation requires a new scope
  decision if the audit contradicts the metadata-only premise.
- **Pseudobulk blacklist:** use only the existing frozen human aux groups G1/S,
  G2/M, Mito, Ribo, and TCR, **leaving the aux file unchanged**. Treat sex genes
  exactly as the existing file does: add neither a sex-gene exclusion nor a
  special sex-gene exemption. The user's clarification supersedes the earlier
  all-X/Y proposal; no chromosome annotation is needed.
  Exclude the five-group union before normalization / fitting; select the top
  2,000 eligible genes afterward. This needs a separate BL fit/cache identity.
  Retain ordinary non-BL fits when validated. Genes are not excluded merely
  for being sex-linked; membership in one of the five chosen groups still applies.

### Corrected batch views: separate unsupervised Harmony and Harmony + limma

The user requested two explicitly named unsupervised ECODA results for every
active **batch-corrected** dataset, distinguishing cell-level integration from
sample-composition correction:

| Result identity | Cell annotation | Sample composition | Sample-level limma |
|---|---|---|---|
| `ecoda_unsupervised_harmony` | Leiden on Harmony-corrected cell PCA | Existing count + 0.5 CLR transform | **None** |
| `ecoda_unsupervised_harmony_limma` | The exact same Leiden assignments | The same initial CLR matrix | Apply the established technical-only limma correction **once** |

- Use the current batch unsupervised default: HVG2000 and Leiden resolution 2.
  Both branches have identical cells, sample order, cluster identities,
  composition counts, zero handling and pre-limma CLR features. Build that
  common input once; derive two separate feature/distance/score bundles.
- Reuse the corrected view's verified Harmony embedding and Leiden labels.
  No second H5AD, preprocessing run, Harmony fit, Leiden fit, or Pipeline 4
  annotation run is required solely to add the second composition branch.
- Limma in the second branch uses only the dataset's configured technical
  batch keys, with intercept-only preservation and no biological covariates.
  Keep author-annotation ECODA correction behavior as already agreed.
- Give the two results distinct artifact/cache identities and method-specific
  provenance. Record the Harmony source/Leiden obs key, HVG size/resolution,
  technical batch keys, and whether sample CLR limma was applied. The
  Harmony-only branch must never claim a limma correction model merely because
  its input view is named corrected.
- Update R/Python correction-contract allowlists/validators, bundle expansion,
  selection and sync/audit expectations, method labels, timing records and
  analysis readers for both identities. Preserve strict validation for other
  methods; do not remove corrected-result provenance checks to permit the new arm.
- In corrected comparisons, replace the ambiguous single unsupervised entry
  with these two named entries, rather than counting the old unqualified
  `ECODA_seuratres_2` as an additional third arm. Preserve historical bundles
  and map their meaning only with explicit source/correction evidence.
- Report both arms separately in corrected funkyheatmaps, MDS panels, score
  tables and applicable method comparisons, using all five confirmed metrics.
  Include both Kidney evaluations for each arm. They are separate method
  comparisons, not separate datasets or additional biological evaluation axes.
- Scope is the eight corrected batch views. The existing uncorrected-view
  comparator remains non-Harmony PCA -> Leiden -> CLR without limma; do not
  label it Harmony. Benchmark unsupervised method definitions remain unchanged
  by this specific addition.

### Kidney: exactly two biological evaluations, one computation

Keep `columns.label = "condition.l1"` and `columns.batch = "suspension_type"`.
Add optional named evaluation metadata under the **dataset-level** columns:

```json
"evaluation_labels": {
  "condition": "condition.l1",
  "tissue": "tissue"
}
```

- Preserve scalar `label` behavior. Extend both R/Python loaders to expose the
  optional mapping; datasets without it retain one primary evaluation.
- Kidney gets two evaluation entries in both batch analyses, funkyheatmaps,
  MDS panels and score exports. Keep the real dataset ID; do not create aliases,
  another H5AD, another preprocessing view, or extra embedding jobs.
- Reuse each validated method distance matrix and the same MDS coordinates;
  calculate ANOSIM, adjusted kNN-3 modularity, mean Ward/PAM ARI, silhouette,
  and LISI separation separately for each axis.
- Normalize by evaluation and metric separately. **Initially count both Kidney
  evaluations as independent, equally weighted entries** in overall/per-metric
  summaries. The user may revise reporting/weighting after seeing results.
- Tissue categories after Papilla removal are literal `kidney` (45),
  `cortex of kidney` (23), and `renal medulla` (12): both axes score all 80.
  `region.l1 == "NA"` is a literal source value mapping to `tissue == "kidney"`.
  Never infer Cortex/Medulla for generic biopsies. Region remains diagnostic,
  not a third scored evaluation.
- Existing NAS obs audit proves one-to-one region/tissue partitions and
  per-specimen consistency. Repeat against the selected Bamboo producer.
  Existing sample Feather exports lack `tissue`; export it explicitly.
- The generic kidney group is perfectly confounded with cell suspension;
  Cortex/Medulla are nuclei. Describe this as recorded tissue/sampling category
  recovery, not unconfounded anatomical recovery. It is not a correction input
  or a reason to add new exclusions.
- For `ECODA_authors_HR_NULL`, deterministically shuffle each axis's own labels
  using the established seed 123. Scoring null distances against true tissue
  labels would incorrectly duplicate the ordinary HR result.
  For an existing null artifact, preserve its recorded shuffled label vector;
  do not replace it by a new permutation to complete missing scores.

### Five-metric results, figures, and rankings

The user explicitly expanded the default reporting contract from three to
**five metrics**: ANOSIM, adjusted kNN-3 modularity, mean Ward/PAM ARI,
silhouette, and LISI separation. This supersedes earlier three-metric defaults
in this plan, `notebooks/AGENTS.md`, and existing figure selectors.

- Include all five in refreshed result/score exports, main and extended
  funkyheatmaps, presentation summaries, MDS separation summaries, and
  applicable per-metric supplementary figures.
- Min-max scale each metric within dataset/evaluation across the compared
  methods, consistently for the two added metrics and the existing three.
  All five contribute **equally** to overall/per-evaluation summary scores
  and resulting method rankings. Preserve explicitly rank-based diagnostic
  presentations as rank views over the agreed five-metric inputs.
- Add Silhouette and LISI separation summary columns on the right of
  funkyheatmaps, after ANOSIM, Modularity and ARI. Extend plot dimensions and
  expected result cardinalities accordingly.
- Existing `calc_sep_score()` already computes these metrics. LISI is already
  transformed to a higher-is-better separation score, `(N - LISI)/(N - 1)`,
  bounded to [0,1]; do not invert it again. Preserve the current scorer and
  its distance/perplexity behavior rather than changing its method here.
- Reuse validated persisted scores first. Complete only missing metrics from
  existing distances and exact named evaluation labels, in new derived outputs.
  Adding metrics alone does not justify preprocessing or embedding reruns.
- Reject incomplete/nonfinite metric grids before aggregation; do not silently
  drop metrics/rows or average a mixture of three- and five-metric results.
  Respect scITD's documented dropped-sample universe and NULL label provenance.
- Refresh **all affected maintained figures**, including screens and per-metric
  supplements. Supp fig 1/2 already select five metrics in part; audit/reuse
  that implementation and refresh inputs rather than duplicating the work.
  An explicitly ANOSIM-only plot can remain ANOSIM-only with refreshed values.

## 3. Evidence and implementation boundaries already established

- Config already contains the new Stephenson exclusions, Papilla exclusions
  in both full-kidney batch views, Alzheimer/KPMP sex removal, and Diabetes
  disablement. Other requested batch selections remain to be encoded.
- Stage 3 applies `subset_vars` before filtering/normalization/HVG/PCA/Harmony.
  Reuse that boundary; a general new Stage 2 subsetting pipeline is unnecessary.
- Sample-level `majority_v1` technical winners underlie the confounding tables.
  Raw cell metadata may vary within a sample for Alzheimer assay, Breast
  dissociation, and Lupus batch. Direct cell-level predicates can split samples.
  Materialize audited whole-sample membership from the approved sample policy,
  then retain all cells for selected raw sample IDs. Keep column authority and
  scientific selection in the registry; bind resolved IDs/checksums in run manifests.
- Legacy Lupus uncorrected versus corrected exports differ in 33 batch
  assignments. Resolve both new views from the same authoritative source and
  aggregation policy, rather than copying mismatched historical lane metadata.
- **COVID temporal verification completed:** the exact Feather underlying
  `plots/batch_effect_bio_vs_batch/Covid19_PBMC_bio_vs_batch_datasets.pdf` has
  MD5 `93085b5c8fd83d6eb8c0da10ca1c0a46`, matching the plot provenance. Its
  111 samples are 20 literal controls plus 91 finite numeric days 1-30; none
  exceed 30 or have missing/malformed/nonfinite days. Its batch x biology
  crosstab exactly equals the saved plotted CSV. Summary RDS MD5
  `57c9b8fbb56427e67a1c2c96d5791f54` and sample order/labels also match.
  A metadata-only read of the existing corrected-final H5AD confirms the
  same 111 samples / 702,665 cells and sample-consistent day/batch/labels.
  Both predicates give d02=5/13/4, d09=2/7/3, d13=3/7/0 in
  control/mild/severe order: 44 samples total. No replacement temporal-filter
  plot or preprocessing is needed for this clarification. Older uncorrected
  metadata (164 samples, including 52 above day 30) is a different historical lane.
- Stage 3 checks integrity before reuse but does not adequately prove the
  current subset/sample membership matches an old same-name output. Add
  cohort identity verification and use new derivative destinations.
- Bassez fallback helper and Stage 2 dispatch exist. The historical worker
  replaces staged RDS in place; adapt any required repair to write a derivative.
  Stage 4 merges automated labels; it is not the author-HR fallback fix.
- Canonical BL selection explicitly preserves the old no-op; legacy code also
  filters on a literal flag instead of the resolved genes. Fix active paths and
  the regression expectation. Loading `default_black_list` from pinned STACAS
  is also wrong; package-name changes alone do not fix the method.
- Frozen `aux/genes.blocklist.rds` Hs union is 1,928 unique symbols; TCR has
  228 entries. SHA256:
  `0bbf15296350940082b843f6cc8505f24a14f3365af569435aa2c7598d6a3113`.
  Its broad Ribo/TCR definitions are intentionally chosen over HiTME's smaller sets.
- Audit exact stable-ID/symbol matching for the frozen five-group union,
  including standardized symbols and duplicate suffixes. Record matching
  coverage and per-group/union exclusions; resolve ambiguity without guessing.
  A chromosome annotation is no longer required: sex-gene exclusion was withdrawn.
- Preserve public `Pseudobulk_hvg2000_bl` / `hvg2000_bl` identifiers for readers,
  but give the corrected policy a new provenance version. Replace misleading
  internal `default_without_sex_genes` usage with the explicit resolved policy.
- Existing heterogeneous benchmark selection files can union methods and
  deduplicate Stephenson/Bassez pseudobulk. R workers currently lack sufficiently
  narrow BL/HR variant selection; add the minimum selector required to leave
  valid variants outside recomputation.
- Corrected unsupervised two-stage correction is confirmed in current code:
  Stage 3 stores Leiden assignments from Harmony PCA; `rename_leiden_cols()`
  selects the corrected/Harmony key for `RNA_snn_res.2`; batch composition
  builds `ECODA_seuratres_2` and runs it through the same `correct_result()`
  limma wrapper as author-HR composition. Separate both branches at this shared
  composition boundary, before limma, using the new identities above. Current
  corrected R/Python contracts assume the old unsupervised method has the
  limma model; extend that identity contract rather than reusing it incorrectly.
- The dispatcher owns one analysis context per invocation. Three contexts need
  separate invocations. A shared durable gate may envelope a lightweight launcher;
  it must not make scientific job survival depend on sibling success.
- Existing IDs isolate logs/manifests, **not scientific outputs**. Stage 3/5 and
  sync still target canonical paths. New derivative output/input binding is
  required before running, not merely a new RUN_ID.

### Bassez current-state audit and conditional method impact

On 2026-11-08 a bounded read of only `Sample`, `sample`, `cellType`, and
`cellSubType` in the mounted NAS file
`/Volumes/Shared/Projects/ECODA_paper/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad`
confirmed the missing-HR counts above. No X/raw/count matrices were read or
copied; no full H5AD checksum was recomputed. This proves the inspected NAS
metadata state, not complete input integrity or the current Bamboo copy.

The Bamboo SSH attempt failed hostname resolution inside the sandbox; its
permission-enabled retry was aborted before connecting. No remote read or
mutation succeeded. Before any reuse or repair, verify the exact selected
Bamboo producer, sidecars/ownership and dependency identity. Do not infer that
the existing helper was applied from its presence or historical plan claims.
Trace and repair source routing/validation if needed so future production
outputs apply the same fallback, rather than fixing only a one-off analysis.

For an **HR-only metadata change**, with cells, samples, counts, features,
embeddings, LR and clinical labels proven unchanged, the dependency audit is:

| Method/result scope | Required response if repaired HR labels change |
|---|---|
| Author-HR composition | Regenerate `ECODA_authors_HR`, `_NULL`, `GloProp`, `Freq_highres` |
| HR composition variants | Regenerate top-variance thresholds 0 through 0.9, 3-most/2-least/3-least-variable-type variants, and PCA dimensions 2/3/5/10/15 |
| HR CT pseudobulk | Regenerate `Pseudobulk_CT_HR_hvg2000` and `Pseudobulk_CT_HR_hvg500`; LR CT remains reusable |
| scPoli | Regenerate seven HR combinations: HVG1000/dim15, HVG2000/dims2/3/5/10/15, HVG3000/dim15; preserve the LR combination |
| PILOT and QOT | Regenerate each method's HR HVG1000/2000/3000 combinations; preserve LR |
| PILOT-GM | HR category count controls model component count: regenerate HR combinations when that effective count changes; preserve LR |
| Transformation and zero-imputation | Both consume configured HR composition; regenerate `trans_analysis` and `zeroimp_analysis` |
| Other composition | Author LR, averaged PCA, deconvolution, Leiden and unchanged automated-annotation results are reusable |
| Other baseline families | MrVI, GloScope, MOFA and LR scITD are reusable |
| Whole-sample pseudobulk | Ordinary/PCA variants remain reusable; the BL variant changes independently under the approved global blacklist repair |

Update HR category-count summaries and all affected figures. Reuse is
conditional on exact input-dependency proof, not solely scientific intuition.
Changing the whole-file checksum does not by itself mean every method's
scientific input changed. Conversely, a Stage 3 rebuild that changes the
features or cell universe invalidates this narrow reuse matrix and requires
a revised scope decision. **The user selected this targeted metadata repair
and dependent reruns after reviewing the report (Q22).** Existing valid work
remains outside recomputation selections. The separately approved global BL
repair includes Bassez once; it is not duplicated by the HR selection.

## 4. Bamboo, preservation, and execution contract

- The user explicitly returned work to **Bamboo** after maintenance. Bamboo
  holds the most up-to-date synchronized data; Yggdrasil is stale/unmaintained.
  Record this as user-reported transfer completion, then perform scoped integrity
  reconciliation on Bamboo. Do not automatically resync from Yggdrasil.
- The checked-in durable profile already pins Bamboo; harmonize contradictory
  AGENTS/skill/active-plan prose and verify profile/preflight against actual
  storage/runtime. Do not mass-replace historical host evidence.
- Author locally, commit/push, and pull the exact commit on Bamboo. Resolve
  remote home/scratch/NAS/runtime paths rather than assuming old mount state.
- Execute full-cohort work only through the checked-in durable profile and
  commit-keyed immutable snapshot/versioned runtime. Preserve Pixi constraints,
  lockfile, no-label-leakage rules, counts/sample-order contracts and QOT source.
- Original H5ADs, RDS/cache bundles, sidecars, manifests, and gates stay immutable.
  Add explicit derivative output roots and an input-artifact manifest keyed by
  dataset/view, binding path, checksum, provenance and cohort identity. Unchanged
  benchmark H5ADs remain read-only references rather than copied full files.
- Mirror only the selected new artifacts to a new finite NAS campaign root.
  Bind roots in ownership, worker, exporter, sync and auditor contracts; changing
  a shell variable alone is insufficient. Readers consume the approved release
  manifest, which can reference validated reused and new artifacts.
- Canonical overwrite/promotion is outside this initial plan. Preserve selected
  NAS file/directory permission normalization and checksum verification.
- No full-cohort data or expression matrices on the workstation. Read metadata
  on HPC/NAS, process one dataset at a time, and transfer only small explicit
  debug/metadata/derived artifacts.

### Parallel execution and failure isolation

Use existing Slurm arrays, watchdogs and explicit selections. Add only a small
launcher/status summary for three concurrent contexts: benchmark, batch
uncorrected, batch corrected. Preserve their separate IDs, logs, selections,
statuses, scientific output roots and one synchronization owner per root.

- Launch independent contexts without fail-fast chaining. A failed child must
  not cancel sibling arrays/watchdogs or prevent another child from starting.
- Once submitted, Slurm jobs survive launcher/SSH disconnection. Persist every
  submission ID before advancing; reconcile partial/ambiguous submission, never
  replay the whole launcher blindly.
- Narrow preparation dependencies to the relevant dataset and consumers.
  Current whole-view preparation watchdog dependencies can block unrelated
  datasets after one failure. Preserve one explicit selection wave and use
  dataset-scoped completion evidence/dependencies within it.
- Only real descendants of failed prerequisites are blocked. Independent rows
  continue; OOM-only retries remain row-scoped and bounded.
- Wait for every launched context to settle, even when another fails. Preserve
  successful outputs and record per-context failure. A summary failure never
  means all work must be rerun.
- Synchronize successful contexts independently with their own root owners;
  a failed context cannot report complete publication. Overall completion
  requires all selected work and reviews. Never mark a failed durable manifest
  complete or rewrite it to hide partial failure.
- Keep the existing `ecoda-benchmark` serialization group; do not invent groups
  or run separate same-root gates concurrently. There is one durable envelope,
  exactly one unbounded durable wait, all emitted scheduler/watchdog IDs, and
  one terminal accounting inspection, explicit child run audits, and reviewer approval. A failed
  campaign preserves success evidence for later targeted recovery/release;
  it does not authorize consumption solely from a successful status string.
- Do not introduce a general workflow engine, distributed queue, or replacement
  orchestration framework. Shared publication is distinct from compute survival.

## 5. Ordered work packages and acceptance

| Package | Status | Work and completion criterion |
|---|---|---|
| P0: safety checkpoint | DONE | Six logical commits pushed to origin/master; remote head verified as `6a1726d4a79f74998c54e968bb045b452b9c93bd` |
| P1: documentation | DONE | NOTES, onboarding README, both AGENTS guides, Bamboo durable-skill examples and historical-plan pointer harmonized; scoped whitespace and skill validation passed; Git publication recorded below |
| P2: scoped contracts | IN_PROGRESS | Cohort/BL core, derivative roots and two corrected Harmony arms committed; exact rerun selectors, failure isolation, Bassez Stage 2 wrapper and gate evidence integration continue |
| P3: source reconciliation | IN_PROGRESS | Orchestrator handles bounded read-only Bamboo source/runtime/NAS preflight; targeted Bassez scope is already user-approved, subject to unchanged-input proof |
| P4: changed inputs | DONE | All16 batch derivatives and53-sample Stephenson input individually reviewed; Bassez metadata repair and Stephenson dual annotation validated/synchronized/approved |
| P5: method waves | IN_PROGRESS | Run three independent concurrent contexts with exact deduplicated scope, durable evidence, focused recovery and root-specific sync |
| P6: evaluation/reporting | IN_PROGRESS | Complete five-metric scores, score both Kidney axes, refresh all affected maintained reports/figures from approved results, verify plots and document limitations/provenance |

### P0: first implementation action, before harmonization or pipeline changes

At inspection: `master...origin/master`, HEAD `70b6033`, 12 modified tracked
files and 24 untracked files; no staging. Recheck before acting. Existing
changes belong to the user and must be preserved, including work outside
this plan. Review the push destination `origin` (carmonalab/ECODA_paper), fetch
and reconcile remote state without force-pushing. Do not push to `upstream`.

Suggested logical checkpoint packages:

1. Existing registry decisions and existing onboarding notes.
2. Batch diagnostics/notebooks/helpers/tests and completed diagnostic plans.
3. Local composition subsetting work and its plan updates.
4. Benchmark/Stephenson exploratory utilities, plotting changes, tests and research.
5. Kidney investigations and discussion-only feature-coverage proposal.
6. Deferred architecture/concurrency drafts and this consolidated plan.

Prefer file-level packages when intertwined hunks make finer splitting risky.
Preserve the original work before later fixes. Record validation limitations:
the initial `git diff --check` found three whitespace lines in the architecture
draft. Check JSON and bounded sensitive-content/file-size concerns; do not
describe backup commits as completed scientific validation.

Leave incidental `.vscode/settings.json` (empty local settings) and `Rplots.pdf`
outside source commits, without deleting them. Source commits do not back up
ignored `data/` or `plots/`; authoritative existing HPC/NAS artifacts remain
separately preserved. No broad `git add .`.

### P1: documentation ownership

- This plan owns current decisions, task scope, exceptions, progress and run IDs.
- `NOTES.md` explains scientific/method changes, evidence and limitations.
- Onboarding README explains inclusion, covariates, group balance and tissue
  interpretation, linking this plan for execution details. Fix its stale nine-
  cohort/Diabetes text and old Alzheimer/Kidney/Breast covariate table.
- `AGENTS.md` states Bamboo policy, flags-plus-view participation and durable
  operational invariants; link task-specific selections here.
- Preserve earlier run histories and unrelated discussion-only proposals.
  Kidney feature-coverage investigation does not authorize new masks,
  post-hoc batch labels, exclusions, or additional remediation runs here.

### P2-P4: input and repair details

- Keep selections in both batch views equivalent; resolve sample membership
  from common source evidence, before any normalization/HVG/embedding.
- For Bassez, first inspect current source/staged/processed obs and provenance
  on Bamboo. The user's selected path is validated reuse if already correct,
  or an atomic HR metadata derivative preserving cell/sample/count/features/
  embeddings if only HR metadata is wrong, followed by the targeted reruns
  in section 3. A broader rebuild needs evidence of an additional invalid
  dependency, an updated rerun matrix, and a new user scope decision.
- If benchmark Stage 3 changes, satisfy HiTME/scATOMIC requirements. Validate
  existing per-sample checkpoints for exact reuse where supported; otherwise
  recompute only invalid/missing annotation work. Batch views skip Pipeline 4.
- Stage 3 changed-selection reuse must fail closed. Existing same-name valid
  MD5s alone cannot establish the new sample universe.
- BL-only repair for unchanged benchmark inputs must not recompute unrelated
  ordinary/CT variants. Shared aggregation may be reused only with verified
  counts/input identity; BL fitting and timing have distinct provenance.

### P5: rerun matrix and collision prevention

| Input/change | Required downstream scope |
|---|---|
| All 16 changed batch-view cohorts | Full supported batch suite for each dataset/view |
| Eight corrected-view unsupervised compositions | Two named arms sharing Harmony/Leiden/counts/initial CLR; Harmony-only and Harmony + one limma correction |
| Stephenson benchmark, changed sample universe | Ten baseline families and maintained eligible variants, plus trans and zeroimp |
| Bassez, HR-only metadata repair | HR composition/HR parameter and sensitivity results, HR CT pseudobulk, HR scPoli/PILOT/QOT/PILOT-GM inputs; exact scope follows dependency audit |
| All 11 enabled benchmark datasets | Corrected whole-sample `hvg2000_bl` variant |
| Kidney extra evaluation | Posthoc scoring/plotting only, using the same newly validated distances |

For a Bassez HR-only repair, whole-sample non-BL pseudobulk, MrVI, GloScope,
MOFA, and unaffected LR methods are scientifically reusable if provenance
proves their inputs unchanged. If cells/features change, expand only the
demonstrably invalid rows and record why.

Baseline families: `gloscope`, `mofa`, `pseudobulk`, `composition`, `scitd`,
`mrvi`, `scpoli`, `pilot`, `qot`, `pilotgm`. Batch dispatcher suite:
`prepare_pseudobulk` plus `pseudobulk`, `gloscope`, `composition`, `mrvi`,
`pilot`, `qot`: **112 declared dataset/view/method rows** for 16 views before
reuse/sharding. The corrected composition row now produces both named
unsupervised arms; this does not add a preprocessing or dispatcher method row.
Expand and record result-artifact/score counts for the extra corrected arm
separately. These are selection rows, not scheduler-task counts.

Build the benchmark selection as a union, then deduplicate by exact artifact
scope/variant: Stephenson all-method work includes its BL result once; Bassez
HR work and its whole-sample BL result likewise cannot collide. Screens may
need distinct existing analysis invocations; inventory their artifact paths
and place genuinely shared writes under the same synchronization owner.

Record exact wrappers, selected rows/variants, expected counts and concrete
destinations before prepare/launch. If emitted scope exceeds approval, apply
the repository's cancellation/FAILED inspection rule for that launch. That
integrity incident is distinct from an ordinary isolated worker failure.

### P6: evaluation and presentation

Extend the existing long score table with `evaluation_id`/label column and
source/analysis provenance. Unique keys include dataset, view, evaluation,
method/variant and metric. Preserve original bundle labels/scores; derived
evaluation outputs have their own ownership/checksums. Adding an evaluation
specification alone never triggers preprocessing or retraining.

Use all five confirmed metrics in default aggregate figures and results.
Complete missing scores from existing validated distances rather than
recomputing methods merely to add silhouette/LISI. Include both right-hand
metric summary columns and refresh the applicable per-metric supplements.
Report both Kidney axes and the initial two-entry weighting explicitly.
Refresh comparison metadata/confounding tables so they describe the new
cohorts, rather than silently relabeling historical matrices.

Publication figures are preserved; no renumbering. Figure 2A retains default
settings, Supp fig 15 extended methods, Supp fig 2 parameter screening;
exclude legacy `ECODA_PB_combo_*`. Execute affected plotting chunks, not a full
knitr render. Review labels, clipping, numeric precision, scales and layout.

## 6. Focused validation and release criteria

- R/Python loaders retain scalar-label compatibility and expose only declared
  evaluation mappings; no view-level column overrides.
- Flags exclude inactive datasets even when views exist; production selections
  are exactly 11 benchmark datasets and 16 batch dataset/view rows.
- Whole-sample masks, raw/canonical IDs, majority-policy winners, min-count and
  exact 20%/10% boundaries; structural zeros and Ncl exception are explicit.
- Corrected/uncorrected sample universes and exported evaluation labels match.
  Old subset identity cannot pass changed-cohort reuse.
- BL actually removes the chosen five groups before fitting, without adding a
  sex/chromosome mask; correct
  identifier namespace, no silent missing-list fallback, no legacy no-op;
  ordinary variants and validated source artifacts remain unchanged.
- Bassez valid HR is retained and fallback covers every eligible cell; exact
  input-dependent rerun selection, no use of biological outcomes as annotation.
- Two Kidney scores use identical distances/coordinates, sample-aligned metadata,
  separate normalization/keys and axis-specific null labels. Changing only the
  evaluation specification submits no Stage 3/5 method compute.
- Five-metric coverage, direction and scaling are consistent in benchmark and
  batch outputs, including both Kidney evaluations. LISI is not inverted twice;
  all five enter means/rankings equally; absent/nonfinite metrics fail rather
  than silently changing denominators. Existing five-metric screens remain valid.
- Corrected unsupervised arms use exactly the same Harmony Leiden assignments,
  sample/cluster counts and initial CLR. The Harmony-only branch calls no
  sample-level limma; the second applies technical-only limma exactly once.
  Method contracts, caches and readers distinguish both correction histories,
  reject mislabeled legacy reuse and expose both results in five-metric plots.
  Adding the branch submits no extra preprocessing/Harmony/Leiden work.
- Inject one row failure, one child-wave failure, partial submission failure,
  and launcher/transport loss: unrelated submitted work survives; IDs/statuses
  remain recoverable; only actual dependents block; recovery selects failed rows.
- Stub Slurm/NAS to test root isolation, single sync owner, checksums, permission
  handling, OOM-only retries, manifest scope and fail-closed terminal review.
- Run narrow existing regressions and `_debug`; no full cohorts for routine QA.
  All required source/runtime/artifact audits pass before release eligibility.

## 7. Decision closure and implementation prerequisites

**No unresolved scientific scope preferences remain from the clarification
rounds.** The user selected targeted Bassez HR metadata repair and dependent
reruns after the impact report. The user subsequently confirmed shared
understanding and explicitly authorized implementation; the grilling
confirmation requirement is satisfied. P0 has completed.

Bamboo producer reconciliation, exact input identities and focused contract
checks are implementation/prelaunch prerequisites; they are not unanswered
questions for the user. If evidence contradicts an agreed cohort or the
Bassez metadata-only premise, report the concrete discrepancy and return only
that decision to the user. Do not silently broaden reruns or change selections.

Resolved latest replies: COVID retains both temporal and batch predicates;
the existing plot already satisfies the temporal restriction. Refresh all
affected maintained figures and results, including silhouette/LISI. All five
metrics contribute equally after min-max scaling to aggregate means/rankings.
Corrected unsupervised ECODA now has separate Harmony-only and Harmony + limma
arms, sharing the same cell annotation and pre-limma composition inputs.

Discoverable prerequisites stay with agents: current Bamboo source/runtime/
storage state, exact retained raw IDs, blacklist identifier matching, Bassez coverage,
screen output paths, and finalized per-variant job counts. Report scientific
contradictions to the user; do not turn them into silent defaults.

## 8. Progress and launch ledger

| Date | Event | Evidence/status |
|---|---|---|
| 2026-11-08 | Read-only repository/metadata audits and decision rounds | Participation/filter tables checked; Lupus typo corrected, COVID d07 dropped, Ncl exception approved; no source implementation or compute |
| 2026-11-08 | Kidney two-axis evaluation | User chose condition.l1 + literal tissue, all 80 samples, initially two independent evaluation weights |
| 2026-11-08 | BL policy, latest instruction | User chose existing aux five groups including TCR, pre-fit filtering; subsequently withdrew sex-gene exclusion. No added X/Y mask or chromosome-reference prerequisite |
| 2026-11-08 | Orchestration revision | User accepted independent concurrent waves, narrow dependencies, retained successful work and a small launcher |
| 2026-11-08 | This plan consolidated | Documentation-only change; remaining decisions above pending |
| 2026-11-08 | COVID day-cutoff clarification | User retained both filters; checksum/metadata/obs audit proves original matrix already day-filtered and selected batches retain 44 samples |
| 2026-11-08 | Bassez decision checkpoint | User requested an exact annotation/method impact report first and a new scope question afterward; repair remains unapproved |
| 2026-11-08 | Five-metric expansion | Refresh all affected maintained figures/results; add scaled silhouette/LISI columns and per-metric panels; user confirmed equal five-metric means/rankings |
| 2026-11-08 | Bassez bounded NAS obs audit and dependency report | Current published H5AD has 50,809 missing HR labels; all have valid LR fallbacks; mixed labels give 40 categories. Exact conditional affected/reusable scopes above. Bamboo SSH read did not complete; no repair or compute performed |
| 2026-11-08 | Bassez repair scope selected (Q22) | User chose targeted HR metadata repair and dependent reruns, with unchanged-input proof and Bamboo producer verification; no broader preprocessing rebuild authorized |
| 2026-11-08 | Corrected unsupervised ECODA split | User requested `ecoda_unsupervised_harmony` and `ecoda_unsupervised_harmony_limma`; code path confirms both correction stages currently occur. Add two composition branches from the same Harmony/Leiden inputs; no extra preprocessing |
| 2026-11-08 | Implementation authorized | User confirmed shared understanding, requested checkpoint commits/push to master, documentation harmonization, subagents, and parallel subset/Bassez/BL work |
| 2026-11-08 | P0 complete | `a325e82` registry decisions; `3039302` batch diagnostics; `b15288d` local composition subsetting; `ae75c1c` benchmark/Stephenson diagnostics; `0d89c6c` Kidney research/limitation; `6a1726d` proposals and this plan. Push succeeded; `git ls-remote` independently verified remote master `6a1726d4a79f74998c54e968bb045b452b9c93bd` |
| 2026-11-08 | Bamboo connectivity and bounded metadata preflight | Scoped permission-enabled SSH succeeded. Home `/home/users/h/halterc`, host `login1.bamboo`, scratch exists, NAS mounted with about 90T free; checkout clean at `b1854c6e786c998f3a786e6885e3119c11d772e2` (must pull exact new source later). Configured py-cuda13 Python is 3.13.14. Obs-only Bassez read independently confirms 74,948 cells / 26 samples / 24,139 valid HR / 50,809 missing HR / zero invalid fallback / 40 mixed categories. No matrix read or scheduler submission |
| 2026-11-08 | Runtime selection prerequisite | Generic `ecoda-py-cuda13.sif` and its manifest are absent on Bamboo; finite runtime-directory inspection confirms versioned images/manifests exist. Select and verify an explicit compatible runtime identity before prepare; no implicit alias fallback or environment refresh |
| 2026-11-08 | BL timing implementation choice | Keep ordinary shared-fit timing unchanged. BL is a self-contained preparation unit with its own timing ID, zero shared-time charge, and measured aggregation + BL fit + selection in its local cost. Mixed ordinary/BL requests run the BL preparation unit separately; this bounded extra aggregation avoids charging an unrelated ordinary fit or inventing a new timing-accounting framework |
| 2026-11-08 | Concurrent external downstream agents | User reported separate agents will edit notebooks and unrelated downstream-analysis files, will receive this plan, and will not commit/push. Preserve and exclude their uncommitted work from this task's Git packages; coordinate overlapping files, particularly the two P1 notebook documentation files and future P6 work |
| 2026-11-08 | P1 documentation QA | Harmonized active host, participation, subsets, pending/applied status, Bassez repair, unchanged aux BL, two unsupervised corrected arms and five-metric reporting. Preserved historical evidence with a supersession banner. Scoped `git diff --check` passed; skill-creator `quick_validate.py` reports `Skill is valid!` |
| 2026-11-08 | BL core and worker review | New BL pre-fit regression failed before implementation, then passed; independent reviewer accepted core and three new BL/worker checks. Existing Bassez benchmark and batch correction regressions also pass. R workers support exact standalone variant outputs; top-level dispatcher/release integration remains in progress |
| 2026-11-08 | Cohort implementation returned | Approved registry filters, evaluation mapping loaders, whole-sample selection/source/config/membership proof, create-only publication and isolated RDS conversion implemented. Six focused checks pass, including evaluation-only mapping reuse; independent review in progress. Shell/sync root wiring remains in progress, so no production launch yet |
| 2026-11-08 | Runtime candidate fully verified | Selected existing format-2 image `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif` and adjacent `.manifest`. Full SHA256 verified as `8fcd00b7a02592f82e8ad96ec9dffdcfb40214b3c47218f770dfb8e828e05e34`; image size 8,863,477,760 bytes; image/manifest/parent have no write bits; manifest dependency hashes match unchanged local Pixi files. Runtime tooling present; this verification submitted no scheduler work |
| 2026-11-08 | Cohort package independently accepted | Reviewer found no material issues and independently passed `test_cohort_selection.py` and `test_datasets_evaluation_columns.R`; existing registry, subset, corrected-source and atomic-preservation regressions also pass. Release remains gated on shell/sync/audit integration and real-source campaign preflight |
| 2026-11-08 | First implementation packages published | Pushed `dfad68d` documentation harmonization, `feb3550` pre-fit BL core and exact worker variants, and `1b6969e` approved whole-sample cohort selection, evaluation mappings and create-only Stage 3 Python publication to origin/master. These source packages passed their scoped tests and independent review; no scientific derivatives or scheduler jobs have run yet |
| 2026-11-08 | Runtime executable probes | The verified relocated image reports Python 3.13.14 and R 4.5.2 under `apptainer exec --cleanenv --containall`; dependency files remain unchanged |
| 2026-11-08 | Bassez truthful legacy input provenance | Full Bamboo source MD5 `22be4eb1f41433b4a8a33dd8e1a3b7d2`, size 3,301,058,992 bytes; NAS and scratch checksum sidecars agree. Neither exact source has a global artifact owner. The worker-specific assumption of a six-field producer record was withdrawn: an explicit legacy-input mode will verify the existing three-field proof and actual content, check ownership conflicts, and record `legacy_validated_input`; it must never fabricate a published producer/run or fall back from an invalid supplied producer record |
| 2026-11-08 | Stage 3 integration review in progress | Derivative input references and scratch/NAS roots are implemented and narrow helper tests pass. Review identified optional Alzheimer producer validation and first-manifest-row view selection; both were fixed and the new focused boundary regression passes |
| 2026-11-08 | Stage 3 isolation package published | `5d37aaabaf0d01172652f2a9d64bc40b7fc669cb` pushed to origin/master. Independent helper, submitter, terminal-audit and Alzheimer boundary tests passed. Bamboo clean detached checkout advanced only after ancestry verification; verified immutable snapshot created under the release source-snapshot parent |
| 2026-11-08 | Input checksum staging complete | Eight unique raw inputs fully MD5-verified, existing checksum proofs compared where present, no active writer before/after; 17-row canonical input references and canonical proof receipts written create-only under `_ecoda_inputs/20261108-confounding-subsets-v1/stage3-inputs-v2`. Original files/sidecars unchanged. Attempt v1 stopped before hashing because the read-only ownership helper required a scope ID; its marker is preserved, v2 completed |
| 2026-11-08 | Bassez worker independently accepted | Metadata-only create/validate/audit worker and focused tests pass independently. Strict parent mode requires the exact canonical run-owned record and wrapper producer validation; explicit legacy mode is separate. Streamed non-target HDF5 data/attribute preservation, valid HR retention, fallback completeness, and create-only outputs covered. Stage 2 shell integration remains in progress |
| 2026-11-08 | Harmony backend independently checked | Both arms have distinct method/model/history contracts; actual counts/CLR and technical limma tests prove one shared source and zero versus one sample-level correction. Four focused/broad regressions pass; standalone artifact identity and shell expansion integration continue |
| 2026-11-08 | External application workstream identified | User named `.agents/plans/20261108-cancer-clustering-application-port-plan.md`, session `01a1bd25-d0d9-7bf2-a34c-31a27de20974`, and explicitly authorized two-way coordination. Read its boundaries and sent a coordination message: application notebook/helper/tests/inputs/figures and Yan archive belong there, with no commits/pushes, scheduler work or pipeline/environment edits. Preserve and exclude its changes from this plan's Git packages. Validated existing HiTME compositions may remain usable independently of author-HR repair; record exact input identities |
| 2026-11-08 | Bassez core checkpoint pushed | `48a6503` contains only the independently accepted metadata worker and preservation tests; pushed to origin/master. No repair jobs submitted; Stage 2 integration continues |
| 2026-11-08 | Stage 3 gate prepared; reconciliation diagnostic | Prepare succeeded for `ecoda_stage3_confounds_20261108_v1`. Read-only reconciliation is stuck on `ssh bamboo ps -axo pid=,command=` before launch. A separate bounded remote enumeration returned timeout 124. Reviewer is diagnosing process enumeration/transport and the global helper's failure handling; no launch retry, no scheduler jobs and no weakened runner identity check |
| 2026-11-08 | Stage 5/Harmony checkpoint pushed | `3388a45` contains reviewed derivative input/output routing, exact scratch/NAS validation, both named corrected Harmony arms and their tests. Independent three-context, scientific-sync and real counts/CLR/limma tests pass. Requested variant propagation and failure-isolation launcher are separate remaining work |
| 2026-11-08 | Exact Bassez dependent scope | 22 HR composition outputs inventoried; two HR CT pseudobulk variants; seven HR scPoli combinations; three HR PILOT and QOT combinations; affected HR PILOT-GM. All maintained trans/zeroimp settings use the HR table, so their complete rerun is within approved scope. Preserve unrelated LR/automated/Leiden/deconvolution/average-PCA outputs |
| 2026-11-08 | Independent BL execution choice | Selected BL tasks use the existing self-contained R path and own their exact BL cache/result/sidecars/timing together, avoiding the shared per-view preparation watchdog. A failed Stephenson full-preparation task must not block unrelated selected BL tasks. Full-family Stephenson defaults retain their existing preparation dependency |
| 2026-11-08 | Process-evidence diagnosis and repair | SSH works; ps all-process and current-user enumeration timed out. Pinned-Python current-UID procfs scan completed in 0.44 s (9 own process command lines, 4,163 visible PIDs). The stalled unlaunched reconcile was stopped; a subagent permission review also stalled 805 s and was canceled. Root now owns permission-requiring operations. Global helper updated atomically after backup; old SHA256 `c2a41dc887bd6c13cfec48c275bbd1ecf6f05f791bdea37aa410f0ccf1318146`, new `1a80a80ec3f79a787788e1e33b7c576aa252a2e7327d9c8248ca8bc0922eb2c9`. Eight focused process-evidence tests and two lifecycle/reconnect tests pass; exact argv/PID/decoy/tmux checks remain. Explicit Bamboo procfs profile changes its digest, so prepare fresh v2 and never launch superseded v1 |
| 2026-11-08 | Current benchmark source inventory | All 11 configured benchmark H5AD paths and MD5/size sidecars exist on Bamboo. Nine unchanged BL candidates still need content verification before use; Bassez and Stephenson must respect their changed-input dependencies. Existing sidecars alone do not prove a source/result is reusable |
| 2026-11-08 | External application evidence handoff | Application agent completed its isolated task without commits/jobs/canonical mutations. Current obs versus legacy HiTME count differences are preserved in `data/applications/cancer_clustering/20261108-v1/derived/legacy_bundle_reconciliation.rds` plus checksum; exact H5AD/bundle identities and scope are in its own plan. No documented historical-filter explanation was found. It used compact current-obs counts, not those bundles. This is an evidence flag for reuse review, not authorization for additional benchmark fits |
| 2026-11-08 | Bassez parent receipt ready | Canonical legacy receipt under `_ecoda_inputs/20261108-confounding-subsets-v1/bassez-hr-input-v1/Bassez_benchmark_analysis_input.md5`; full parent MD5/size and ownership conflicts rechecked; receipt SHA256 `8106926f4ee19edd6cc5264dd8ded7a5fb1c413601a44e06d831f5fb31e2d9cd` |
| 2026-11-08 | Fresh Stage 3 v2 gate reconciled | Profile/tool checkpoint `0e97438` pushed. Gate `ecoda_stage3_confounds_20261108_v2` is PREPARED; reconciliation proves ABSENT with seven passing invariants and no discrepancies. Profile digest `765697b9df3337f5795d3e73e753abf75d196d8d264ae6d1d5513263257e9e1c`; exact wrapper/source/runtime/17-row scope unchanged. Ready for its first launch; v1 is superseded/unlaunched |
| 2026-11-08 | Stage 3 v2 early controller failure; scoped recovery | Controller ran 22:36:24–22:36:44 UTC and exited1: existing parser mishandles `--throttle 17`. No child run root or scheduler IDs created. Single waiter recorded FAILED; terminal inspect failed closed because no IDs exist (`v2.inspect.json`), so no successful audit/release is claimed. Help-only validation confirms all recovery arguments with `--throttle=17` parse successfully. Fresh gate v3 / child run `stage3_confounds_20261108_v2` uses identical17row/source/runtime scope and separate output root; no force or prior-gate relaunch. Source parser regression fix delegated independently |
| 2026-11-08 | Stage 3 v3 controller launched | Fresh gate `ecoda_stage3_confounds_20261108_v3`, child `stage3_confounds_20261108_v2`, controller RUNNING at `2026-11-08T22:41:28Z`; wrapper digest `b8ac3b2e7e92f5304b003757374ed5fae8adc48d2c36130e2513aa806443ce81`. All invariants passed before launch. Exactly one unbounded waiter armed (`v3.wait.json`, local exec session95573). Scheduler IDs not yet observed. Continue independent local work; do not relaunch or arm another waiter |
| 2026-11-09 | Stage 3 v3 failed; exact recovery scope proved | Wait completed FAILED. IDs preflight4663460, array4663462, watchdog4663463; one terminal inspect queried exactly those IDs once and failed as expected (2COMPLETED metadata tasks,18FAILED rows including watchdog). Every expected new H5AD path was checked:17missing,0existing. Worker/watchdog errors occur before preprocessing because compute nodes lack `/srv/smednas...`; preserve all evidence and target only the17missing rows in a fresh source/run |
| 2026-11-09 | NAS boundary repaired and independently accepted | `d4fb10e` pushed: bound NAS identity may be read without a compute mount; fresh bind/publication still require actual NAS. Explicit scratch-only Stage2/3 ownership preserves NAS identity. Bassez NAS sync, producer records, permissions and owner completion moved to login controller. Independent real-H5AD NAS-hidden, restored/missing mount, four-record, corruption, mode and retained-success tests pass; broad existing watchdog/submitter regressions pass |
| 2026-11-09 | Single-query accounting bridge accepted | `6551371` pushed; seven focused tests pass. Preserve complete original accounting output, validate every returned row and exact recorded scope, then emit a new labeled root projection/proof. Run first inspect once, bridge, pinned source audit with the projection, selected NAS integrity checks, then reviewer approval only if all pass. This corrects contradictory old prose ordering without a second query, changed running profile, invented child IDs or ignored failures |
| 2026-11-09 | Exact method scope and coordinator review | R22HR/2CT package, immutable public3/4-field selection, exact7/3/3HR Python rows, valid-owned-prep consumer recovery and three-context concurrent controller accepted independently. Within a failed context, successful scratch artifacts remain preserved but no automatic partial NAS publication is claimed. New Stage5 OOM/NAS ownership and terminal/failure-state findings are being fixed before launch |

### Repaired-source gates — launch identities and current status

Both source snapshots use full commit `7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3`
and source archive SHA256
`aa79db2073ed5aaf138b1c2c88a2021de110e79b1a1e61aff4f859fe637159f9`.
They have **separate snapshot-parent locks**. Runtime remains the verified
relocated format-2 image/manifest recorded above; Pixi and aux inputs unchanged.
Fresh run/log/gate roots were checked absent before preparation. Earlier
failed/prepared gate manifests are not changed or reused. No force flags.

#### Stage 3 missing-only recovery

- Status RUNNING; launched `2026-11-08T23:44:50Z`; gate `ecoda_stage3_confounds_20261108_v4`; child `stage3_confounds_20261108_v3`. Preflight `4663496`, array `4663499`; one waiter local exec `69025`.
- Exact scope: 17 missing H5AD dataset/view rows (16 batch views + Stephenson benchmark); initial128G, OOM-only ceiling500G, throttle17,16CPUs/fourhours perrow; two routine COVID metadata preflight rows.
- Serialization group `ecoda-preprocess-20261108-confounding-v1`; source manifest
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-recovery-v4/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/identity/source.manifest`.
- Local evidence prefix `.gate/ecoda_stage3_confounds_20261108_v4`; remote prefix
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/gates/ecoda_stage3_confounds_20261108_v4`.
- Wrapper SHA256 `0d3833eeaf94c8ed02c43c2b6d8f4330979229854fb8cb899acc051c1b42de97`. Output scratch root
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage3_confounds_20261108_v3/artifacts/datasets`; NAS release root
  `/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper/releases/20261108-confounding-subsets-v1/stage3_confounds_20261108_v3/datasets`.

Exact wrapper:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper STAGE3_INPUT_PRODUCER_RUN_ID=stage2_alzheimer_donor_assay_20260914T172421Z /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-recovery-v4/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-recovery-v4/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-recovery-v4/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage3_confounds_20261108_v3 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage3_confounds_20261108_v3 --script src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-inputs-v2/selection.tsv --derivative-release 20261108-confounding-subsets-v1 --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-inputs-v2/input_artifacts.tsv --throttle=17 --mem 128G --max-mem 500G
```

#### Bassez metadata-only repair

- Status DONE: completed, source/artifact/NAS audited and reviewer-approved; gate `ecoda_stage2_bassez_hr_20261108_v1`; child `stage2_bassez_hr_20261108_v1`. Worker `4663493`, watchdog `4663494`, both `COMPLETED|0:0`.
- Exact scope: one explicit bassez_h5ad_metadata step; one worker plus watchdog;1CPU,16G initial,64G ceiling,fourhours,throttle1; target74948cells/26samples/50809fallbacks/40HRcategories; full unchanged non-target-data proof.
- Serialization group `ecoda-bassez-metadata-20261108-v1`; source manifest
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-bassez-hr-metadata-v1/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/identity/source.manifest`.
- Local evidence prefix `.gate/ecoda_stage2_bassez_hr_20261108_v1`; remote prefix
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/gates/ecoda_stage2_bassez_hr_20261108_v1`.
- Wrapper SHA256 `38e4c92550ac3eeaa72d78d1121058ecff7fb35c5bc1805b035e22093ef8ebb9`. Output scratch root
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage2_bassez_hr_20261108_v1/artifacts/datasets`; NAS release root
  `/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper/releases/20261108-confounding-subsets-v1/stage2_bassez_hr_20261108_v1/datasets`.

Exact wrapper:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-bassez-hr-metadata-v1/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-bassez-hr-metadata-v1/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-bassez-hr-metadata-v1/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage2_bassez_hr_20261108_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage2_bassez_hr_20261108_v1 --script src/2_dataset_specific_preprocessing/1_submit_hpc.sh -- --datasets Bassez --steps bassez_h5ad_metadata --bassez-reference /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/bassez-hr-input-v1/reference.tsv --derivative-release 20261108-confounding-subsets-v1 --mem=16G --max-mem=64G --throttle=1
```

Bassez uses the one-row seven-field legacy reference whose SHA256 is
`43e64b9c5e298c27162b544090a08578d937c6ddbfd24ad098e3959ddc620a10`.
It makes no fabricated original-producer claim. The older48dbf0e proposal was
never prepared/launched and is preserved as a pre-topology proposal copy.
Stage3 reuses unchanged input-reference/selection identities recorded above;
all17 scientific outputs are missing, so all17 are legitimate recovery rows.

### Active assignments

- Orchestrator `/root`: this plan, Git integration/push, bounded Bamboo
  preflight and run/release coordination; no pipeline implementation.
- Implementor/reviewer `docs_checkpoint_audit`: reviewed topology/R selectors
  and coordinator; authored accepted single-child accounting bridge. Next
  reviews and the explicit multi-child accounting extension stay bounded.
- Implementor/reviewer `batch_subset_audit`: cohort/Bassez core, topology and
  exact R selectors returned/accepted; independently reviewed requested scope.
- Implementor `benchmark_dependency_audit`: accepted public request/coordinator
  packages and Stage 5 OOM/NAS followup accepted. Now prepares the exact
  early ten-dataset BL-only launch contract locally. No remote operations.
- Subagents work locally without escalation. Return network or protected
  operations to the orchestrator. Use `followup_task` for a new assignment
  when an agent is completed; `send_message` alone does not resume it.

Checkpoint exceptions: `.vscode/settings.json`, `Rplots.pdf` and the unrelated
`Halter_etal_bioRxiv_ReviewCommons.docx` remain untouched and uncommitted. The architecture draft's three existing whitespace
lines were preserved in its backup commit. Backup commits are not scientific
release validation. No full H5ADs or expression matrices were copied locally.

### First Stage 3 gate: launch specification

- Status: v1 was superseded before launch. v2 controller launched and FAILED
  before creating a run root or submitting scheduler work (spaced throttle
  parser bug). v3 is FAILED after terminal wait/inspection. It submitted preflight4663460,
  array4663462 and watchdog4663463; NAS-dependent compute validation stopped
  every preprocessing row. Do not relaunch prior gates. Recovery needs the
  repaired source and a fresh run/gate.
- Failed gate `ecoda_stage3_confounds_20261108_v3`; child run
  `stage3_confounds_20261108_v2`; host Bamboo.
- Serialization group `ecoda-preprocess-20261108-confounding-v1` protects
  this disjoint derivative preprocessing destination. Stage 5 retains its
  existing `ecoda-benchmark` synchronization group.
- Source commit `5d37aaabaf0d01172652f2a9d64bc40b7fc669cb`; immutable
  source archive SHA256 `a1172287c5a7f65e55ce0a34846ec7f30c192576198aed678e554e1b0574eeb9`;
  snapshot `tree/aux` is the frozen auxiliary root. Reuse the fully verified
  relocated format-2 runtime recorded above; no environment mutation.
- Exact selection: the eight approved datasets, each corrected and
  uncorrected batch view, plus Stephenson benchmark: **17 rows**.
  The submitter also performs the two selected COVID obs-only preflight rows;
  these are validation work, not additional cohort preprocessing. Initial
  memory 128G, OOM-only ceiling 500G, throttle 17 (scheduler availability
  governs actual concurrency), standard 16 CPUs / four hours per worker.
- Control root on Bamboo:
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-inputs-v2`.
  `selection.tsv` SHA256
  `e91880b664163ca6480549f2c5ad29839ce7b95956aef9aa1de08bc48dabe54d`;
  `input_artifacts.tsv` SHA256
  `ed7464e40cdeaa4a6921119bf5859bc857b1cf42865442470876318b0d0dab20`.
- Scratch outputs: the child run's `artifacts/datasets/<dataset>/output`;
  NAS outputs: `releases/20261108-confounding-subsets-v1/stage3_confounds_20261108_v2/datasets/<dataset>/output`.
  No canonical historical output is replaced.
- Local gate/spec/evidence prefix:
  `.gate/ecoda_stage3_confounds_20261108_v3`; remote gate prefix:
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/gates/ecoda_stage3_confounds_20261108_v3`.
- Before release, wait once, audit the exact child run/selection/source/runtime,
  collect every emitted scheduler/watchdog ID, terminal-inspect and review.
  Stephenson benchmark additionally needs a validated Pipeline 4 derivative;
  existing annotation checkpoints cannot be blindly adopted by changed source
  identity. Reuse requires exact retained cells/counts/gene inputs and pinned
  annotator/reference provenance; otherwise recompute the necessary annotations.

Exact immutable wrapper:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper STAGE3_INPUT_PRODUCER_RUN_ID=stage2_alzheimer_donor_assay_20260914T172421Z /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-v1/5d37aaabaf0d01172652f2a9d64bc40b7fc669cb/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-v1/5d37aaabaf0d01172652f2a9d64bc40b7fc669cb/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-v1/5d37aaabaf0d01172652f2a9d64bc40b7fc669cb/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage3_confounds_20261108_v2 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage3_confounds_20261108_v2 --script src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-inputs-v2/selection.tsv --derivative-release 20261108-confounding-subsets-v1 --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-inputs-v2/input_artifacts.tsv --throttle=17 --mem 128G --max-mem 500G
```

Stage 3 v1 was never launched; v2 controller failed before scheduler
submission; v3 failed on compute NAS validation and was terminal-inspected.
No new H5ADs exist and no successful release is claimed. For every
future launch append:
package, source/runtime/aux identities, host and resolved roots, dependency
reviews, exact wrapper/selection and expected row count, gate/child run IDs,
scheduler/watchdog IDs, wait/inspect/reviewer evidence, incident/recovery scope,
and approved output-manifest location.


### 2026-11-09 terminal clearance — Bassez metadata-only derivative

- Exact wrapper/source/runtime and input identities are the Stage 2 launch
  specification above. Source commit `7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3`;
  source archive SHA256 `aa79db2073ed5aaf138b1c2c88a2021de110e79b1a1e61aff4f859fe637159f9`.
- Completion `2026-11-08T23:49:57Z`; one terminal accounting query recorded
  worker `4663493` and watchdog `4663494`, both `COMPLETED|0:0`.
  First inspection, raw capture/projection, pinned run audit and approval are
  separate preserved evidence. No repeated accounting query was issued.
- Pinned audit: child `status/terminal_audit.log` ends
  `ECODA_RUN_AUDIT_OK=stage2_bassez_hr_20261108_v1`. It validates the four
  scratch/NAS payload/repair records, ownership, full checksums and selected
  NAS modes. Bridge proofs live in child `manifests/`; original gate/inspect/
  profile copies live in `manifests/gate_terminal/`.
- Output: `artifacts/datasets/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad`
  under the recorded child run, mirrored to its separate NAS release root;
  size 3,301,360,648 bytes, MD5 `678656c738597471d884f88e345e6970`.
  Legacy input MD5 `22be4eb1f41433b4a8a33dd8e1a3b7d2` is unchanged.
- Only `cellSubType` and explicit repair provenance changed: 50,809 fallback
  fills, 74,948 cells, 26 samples, 33 to 40 categories. The unchanged protected
  H5AD tree covers 178 datasets / 259 objects, SHA256
  `171ead0cf9820c03600dc584099b463f2e2bb78c66f4bcfc532beeca38c85845`.
- Independent reviewer `batch_subset_audit` accepted the evidence and no
  further scientific compute check was required. The operating agent recorded
  formal reviewer approval through the durable helper. This clears the input
  for targeted HR reruns and the independently required BL rerun; those method
  results have not yet been computed.
- Early BL wave scope: exactly ten `benchmark_analysis / pseudobulk /
  hvg2000_bl` requests for Adams, Bassez, Gongsharma_cmv_young_males, Kfoury,
  Kim, Lee, Pelka, Smillie, Wu and Zhang. BL preparation is self-contained
  within those workers. No ordinary pseudobulk preparation row, LR variant,
  Stephenson row or broad force is authorized by this early wave. Exact input
  references, wrapper, roots and resource settings must be appended before
  durable prepare/launch.


### Early BL wave preparation — pending final derivative preflight fix

- Ten-row selection and six-field input references are frozen on Bamboo at
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/benchmark-bl10-inputs-v1`.
  Selection SHA256 `202dfbe2aaaa433ab7abe1b641735fa6c5dd022647931eb6fcbd73070fb5e371`;
  input references SHA256 `69839b146cd86b8032b53117d8e61700d47e40b0a3d5682ce3a6835946dc0dd3`.
- Bassez references its approved new derivative; the other nine reference
  existing canonical benchmark H5ADs. Separate `<dataset>.expected.md5`
  receipts project the original expected MD5/size to the canonical absolute
  payload path. Original payloads and sidecars are untouched. This is an
  expected-checksum projection, **not a new full-content verification claim**;
  durable preflight and workers must verify complete payload MD5 before use.
- `expected_checksum_provenance.json` preserves original receipt bytes and
  SHA256 identities, the explicit pending content check, and the distinction
  between the approved Bassez derivative and legacy inputs. Its SHA256 is
  `c39f404f25becc43037bc49840881be96fd443b327d80feff9206429d736e49c`.
- Read-only launch review found one remaining dispatcher input assumption:
  the legacy sidecar repair preflight compares against canonical NAS even
  when the input binding deliberately selects a derivative release. The
  implementor is adding a narrow derivative-mode read-only validation branch
  with a regression fixture; legacy behavior remains unchanged. This fix is
  now independently accepted, root-tested, committed/pushed as `9b06d15` and
  deployed in the new snapshot. Bound source identity uses the exact selected
  canonical checksum proof; original adjacent receipts remain unchanged.
  Other active preprocessing jobs are unaffected.


### Early BL10 gate — exact immutable launch contract

- Status **FAILED**, terminal-inspected once after preflight failure; launched
  `2026-11-09T00:17:52Z`, finished `2026-11-09T00:22:47Z`; host Bamboo; gate
  `ecoda_stage5_benchmark_bl10_20261109_v1`, child
  `stage5_benchmark_bl10_20261109_v1`; serialization group `ecoda-benchmark`.
  This gate is independent of the active Stage 3 gate. Review it before a
  later same-group Stage 5 gate; no second synchronization group.
- Source full commit `9b06d15cf6388225b2523aeb2677eadf68338937`, archive SHA256
  `b0809b3ffeced3fd6845880786650a0219edaaba0eafdf99b8f2ea760b783b86`.
  Snapshot has its own parent `20261109-benchmark-bl10-v1`; source manifest:
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v1/9b06d15cf6388225b2523aeb2677eadf68338937/identity/source.manifest`.
- Same verified relocated format-2 runtime and manifest as previous gates;
  image SHA256 `8fcd00b7a02592f82e8ad96ec9dffdcfb40214b3c47218f770dfb8e828e05e34`,
  manifest SHA256 `1c4c2cc340353fe86f1b89235f1f3e5c7b3545f7f303d2b85612cfa8cdedfd50`.
  Auxiliary root is the snapshot's `tree/aux`; unchanged blacklist pickle
  SHA256 `0bbf15296350940082b843f6cc8505f24a14f3365af569435aa2c7598d6a3113`.
  Five configured gene groups are excluded before fitting; no sex mask.
- Exact 10 selection rows and their input/provenance digests are recorded
  above. One standalone `hvg2000_bl` result and one BL preparation/cache per
  dataset: **20 RDS outputs**, no ordinary/PCA/CT/aggregate method results.
  No `prepare_pseudobulk` task, broad force, screen, GPU method or Stephenson.
- Expected initial scheduling: H5AD preflight array (10 tasks), R environment
  preflight (one task), pseudobulk array (10 tasks), its watchdog (one task),
  aggregate gate (one task): **five scheduler IDs / 23 scheduled tasks**.
  Zero sample metadata export tasks in this benchmark-only wave. Extra work
  is limited to documented OOM-only retries for failed selected PB rows.
- CPU workers: existing shared-cpu/EPYC-7742/16 CPU contract, four hours;
  initial 128G, OOM-only ceiling 500G. Explicit preflight throttle 10; the
  method array's default throttle is 1000 but its ten-row scope permits only
  ten concurrent tasks. No GPU requests or environment refresh.
- New scratch root `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage5_benchmark_bl10_20261109_v1`;
  outputs under `artifacts/benchmark`, separate NAS release under
  `releases/20261108-confounding-subsets-v1/stage5_benchmark_bl10_20261109_v1`.
  Fresh run/log/NAS destinations checked absent. Existing artifacts immutable.
- Local gate evidence prefix `.gate/ecoda_stage5_benchmark_bl10_20261109_v1`;
  wrapper SHA256 `396a9b39b247a5f2d37a2117368888c8d3b05091578f1ac22472746bd201178b`.
  Prepare, reconcile, launch exactly once and arm one unbounded waiter.
  Before release use one accounting inspection, single-child bridge, pinned
  Stage 5 audit, selected NAS integrity checks, and reviewer approval.

Exact wrapper (no additional broad selectors):

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper BENCHMARK_CPU_TIME_LIMIT=04:00:00 /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v1/9b06d15cf6388225b2523aeb2677eadf68338937/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v1/9b06d15cf6388225b2523aeb2677eadf68338937/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v1/9b06d15cf6388225b2523aeb2677eadf68338937/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_benchmark_bl10_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_benchmark_bl10_20261109_v1 --script src/5_run_benchmark_methods/1_submit_hpc_array.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/benchmark-bl10-inputs-v1/selection.tsv --derivative-release 20261108-confounding-subsets-v1 --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/benchmark-bl10-inputs-v1/input_artifacts.tsv --mem=128G --max-mem=500G --throttle=10
```

### Stephenson benchmark annotation dependency

Pipeline 4's existing canonical entry point cannot yet safely target the new
Stage 3 derivative: it hardcodes canonical input/output paths and modifies
the H5AD during merge. A bounded implementation is underway to reuse the
existing annotation engine with a strict Stage 4 input binding and a new
run-owned H5AD copy. It must preserve the Stage 3 producer, use fresh owned
chunks/checkpoints, dual HiTME/scATOMIC annotation, scratch-only compute
ownership and controller NAS publication. Batch-only views remain exempt.
The intended scope is one Stephenson benchmark input, one preparation task,
27 two-sample annotation tasks for 53 samples, and one merge task, with the
existing phase watchdogs. Final exact resource/reference/wrapper contract
must be recorded after review and before launch.

Audit additionally found that the existing annotation worker exports
`normal_tissue` but otherwise defaults HiTME tissue to Tumor. The bounded new
Stephenson mode must use its already-declared registry `tissue=Blood` and
`normal_tissue=true` explicitly, pinned by the source configuration identity.
This is an additional reason not to adopt unproven legacy model checkpoints.
It does not authorize a broad annotation rerun of other cohorts.


### 2026-11-09 BL10 launch event

Prepare and reconcile passed; reconciliation found ABSENT remote execution
with all invariants green. The controller launched exactly once at
`2026-11-09T00:17:52Z`, wrapper digest
`396a9b39b247a5f2d37a2117368888c8d3b05091578f1ac22472746bd201178b`.
One unbounded durable waiter is armed (local exec session `24353`);
`.gate/ecoda_stage5_benchmark_bl10_20261109_v1.wait.json` is its final output.
Do not relaunch or arm another waiter. Initial log read was empty while the
immutable executor performs its validation; no scheduler IDs observed yet.
The exact expected ten method rows/five initial scheduler IDs remain the scope
contract above. Stage 3's independent waiter `69025` is still active. Method
results are not release-eligible until terminal accounting, pinned audits,
selected NAS integrity checks and reviewer approval pass.


### P6 readiness — metadata export complete; result/evaluation helpers underway

- `3f36464` is pushed: sample metadata exports resolve the dataset-level
  `evaluation_labels` map through the existing config helper. Both Kidney
  batch views export `condition.l1` and `tissue`, retaining 80 samples and
  their 45/23/12 tissue counts. Missing, malformed, blank or sample-inconsistent
  evaluation axes fail closed; technical aggregation cannot absorb biology.
  Existing scalar-only export behavior remains unchanged. New actual-H5AD
  export tests and the existing batch-majority contract independently pass.
- Core `calc_sep_score` already computes all five approved metrics. Reporting
  is incomplete: default benchmark/shared plots and both batch notebooks still
  have three-metric/hardcoded legacy paths; standalone BL/HR outputs are not
  read by the legacy aggregate-only loader. Changes must join explicit approved
  artifact selections, never silently substitute old required BL/HR/Stephenson
  or batch results when a new result is missing.
- Implementor is adding two bounded analysis helpers with focused tests:
  one explicit finite release-artifact ledger and one evaluation/scaling helper.
  Source bundles stay immutable. New derived score tables bind labels,
  distances, input/release identities and scoring policy. Known-compatible
  stored metrics may be reused; missing or incompatible metrics are rescored
  from validated distances, without fitting methods. LISI provenance must
  establish the already-transformed separation definition; do not invert it.
- Preserve primary NULL's actual stored shuffled labels; a scores-only legacy
  artifact cannot supply new metrics without recovering its validated distance
  and label bundle. Do not fabricate a primary shuffle or score true labels.
  Kidney's additional tissue NULL evaluation uses its own explicit seed 124
  and records exact shuffled labels. Default overall ranking is the arithmetic
  mean of the five per-evaluation/per-metric min-max scores, each equally
  weighted. Constant metric ranges use a documented neutral display value
  0.5 in the new aggregation helper, so tied methods remain tied; global math
  utilities and the separate cancer application workflow remain unchanged.
  Existing separately labeled rank diagnostics may retain their definitions,
  but must include all five metrics and not replace this default mean.
- Expected batch report grids: eight physical cohorts per view, nine equally
  weighted evaluation entries (Kidney condition and tissue separate), five
  metrics; uncorrected eight comparators including NULL gives 360 values,
  corrected nine including NULL/two Harmony arms gives 405. Disabled Diabetes
  is absent. These are score entries, not new preprocessing/method tasks.
- External cancer-application agent completed without commits, jobs or shared
  environment changes. Its separate files remain excluded. Its read-only
  reconciliation found legacy HiTME layer-2 composition counts inconsistent
  with current obs for Bassez, Zhang and Pelka (max common-type differences
  25/144/44 cells; Bassez also has one current AS-DC). Evidence is
  `data/applications/cancer_clustering/20261108-v1/derived/legacy_bundle_reconciliation.rds`.
  Flag those exact baseline rows for bounded provenance review before current
  input reuse; this does not authorize broad method fitting or expand the
  approved HR/BL work. The BL10 wave is independent of that reuse question.
- Notebook/plot integration follows helper review; no publication figure is
  refreshed or declared final before its required releases and complete
  finite five-metric grid are validated.


### Bassez conditional PILOT-GM dependency resolved

A bounded metadata-only Bamboo audit read `obs/cellSubType` on the unchanged
parent and approved repaired derivative with the pinned `9b06d15` helper;
no counts, embeddings or method fits were loaded. The current worker's exact
`notna & stripped_nonempty` / `nunique` component rule would request **34**
components on the parent and **40** on the repaired input. The parent has 33
valid HR categories plus one placeholder category accepted by that worker's
older narrow validity rule; this does not change the repair's valid-category
33-to-40 report. We do not independently claim the historical fitted model's
component count.

The approved conditional rerun is therefore required for all three Bassez HR
PILOT-GM variants (HVG 1000/2000/3000). LR remains reusable subject to its
ordinary input provenance. A tiny public token allowlist extension is needed
for the 1000/3000 HR rows; the worker already supports those combinations.
No model/scientific worker or vendored QOT change is required. Metadata audit
evidence: `.gate/bassez_pilotgm_hr_component_count_20261109.json`.


### BL10 v1 terminal failure and narrow recovery dependency

- Single durable wait completed FAILED. First terminal inspection performed
  once for the only emitted scheduler ID, preflight array `4663527`, recorded
  nine failed tasks and one completed Bassez task. No second accounting query,
  bridge projection or success approval is authorized for this failed gate.
- Exact run-owned status inventory: Bassez OK; Adams, Gongsharma_cmv_young_males,
  Kfoury, Kim, Lee, Pelka, Smillie, Wu and Zhang FAIL. Run-owned scientific
  RDS inventory is **zero**. No method fitting, R environment preflight,
  watchdog or aggregate scheduler jobs were reached. Existing data unchanged.
- Compute logs were empty because the worker suppresses contract errors and
  generic checksum validation emits no diagnostic. Pinned Adams metadata/
  schema contract passes on the login host; all ten bound full content hashes
  already passed controller validation. Legacy adjacent receipts name a home
  alias while the runtime uses `--containall --no-home`; Bassez's canonical
  receipt succeeds. The old source-identity adapter did not cover compute-side
  generic checksum/validation-record calls.
- Implementor is wiring an explicit exact-bound checksum proof through all
  compute input consumers and validation-record checks. Unbound checksum
  behavior and ordinary output immutability remain strict. A narrow
  `stage5_preflight` input-validation record mode may preserve input file and
  receipt modes while writing only its consumer-run evidence; it must fully
  hash writable inputs and never rewrite/chmod old inputs or proofs. New
  fixture hides the old home alias and exercises actual worker checks.
  Future worker failures must log phase/return code/path diagnostics.
- After independent review/commit/push, use a fresh run/gate/snapshot parent
  and the same exact ten missing BL result rows, no force and no unrelated
  rows. Preserve the failed gate and all evidence. Stage 3 remains independent
  with its original waiter active; do not disturb it.

### Stage 4 implementation checkpoint (not an annotation run)

`51fc6ce` contains the independently accepted 14-file Stage 4 derivative
package. Root also ran the real create-only copy/preparation/merge fixture;
full controller, existing derivative, worker-policy, watchdog ID, reuse/cutover
and shared-run regressions pass. Parent X/raw/counts/embeddings/graph/order/
author labels are preserved. Copy receipt checks are required at merge,
publication and terminal audit. New lane explicitly uses registry Blood and
normal_tissue, validates pinned references without downloads, retains NAS
identity on compute and publishes only on the login host.

No Stage 4 HPC job is launched yet: it waits for the approved Stage 3
Stephenson benchmark output. Exact future scope is one dataset/view row,
53 samples, 27 two-sample dual annotation tasks, one preparation task and one
merge task, plus three phase watchdogs: **six initial scheduler IDs / 32
scheduled tasks**. Set `STAGE4_INPUT_PRODUCER_RUN_ID=stage3_confounds_20261108_v3`
and use a fresh committed immutable source/runtime and new run-owned output.
The 1000/3000 HR PILOT-GM public tokens are separately checkpointed in `7a474ef`;
focused exact-three-token selection and LR rejection tests pass.


### BL10 v2 exact recovery launch contract

- Status RUNNING; launched `2026-11-09T00:54:50Z`; fresh gate `ecoda_stage5_benchmark_bl10_20261109_v2`,
  child `stage5_benchmark_bl10_20261109_v2`; Bamboo; group `ecoda-benchmark`.
  Failed v1 remains FAILED and terminal-inspected. No scientific RDS exists
  there, so all ten requested BL results remain legitimate missing rows.
- Explicit proof fix `084f512` is reviewed, root-tested and pushed. Actual
  isolated worker fixtures hide the old home alias, retain writable source
  modes, pass the bound proof through compute/downstream checks, and reject
  same-size mutation, foreign proof and unsupported preserve-input exceptions.
  Existing default output freezing and legacy/unbound checks remain unchanged.
  Shared run audit, derivative-input and R worker regressions pass.
- Source full commit `084f512cf89f677fd72b98937980a456f6842b1f`; separate snapshot
  parent `20261109-benchmark-bl10-v2`; archive SHA256
  `29e114df96db92117516890adb75ddee92a85a4db5075e1334a497fa0504b09d`.
  Source manifest `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v2/084f512cf89f677fd72b98937980a456f6842b1f/identity/source.manifest`.
- Runtime/image/manifest/aux and all ten source/checksum-reference identities
  are unchanged from v1; selected blacklist aux remains unchanged, no sex mask.
  Exact selection SHA256 `202dfbe2aaaa433ab7abe1b641735fa6c5dd022647931eb6fcbd73070fb5e371`;
  input reference SHA256 `69839b146cd86b8032b53117d8e61700d47e40b0a3d5682ce3a6835946dc0dd3`.
- Exact scope/resources unchanged: ten self-contained `pseudobulk/hvg2000_bl`
  rows, twenty RDS outputs (cache + standalone each), zero ordinary prep,
  screens, metadata export, CT or GPU tasks; initial five scheduler IDs /
  23 tasks (10 H5AD preflight, one R environment, ten fits, one watchdog,
  one aggregate), plus only documented failed-row OOM retries. Initial128G,
  ceiling500G, CPU16/four hours, explicit throttle10, no force.
- New scratch root `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage5_benchmark_bl10_20261109_v2`, outputs `artifacts/benchmark`;
  NAS release `releases/20261108-confounding-subsets-v1/stage5_benchmark_bl10_20261109_v2`.
  Fresh scratch/log/NAS roots checked absent and no environment refresh active.
  Exact wrapper digest `0cc285d59729b806d066b062db7a8b61e4540d18389dbce8cf800f676245fc40`.
- Local evidence prefix `.gate/ecoda_stage5_benchmark_bl10_20261109_v2`.
  Prepare/reconcile/launch once, arm exactly one durable waiter; terminal
  inspect once with every emitted ID, then bound accounting bridge, pinned
  Stage5/NAS audit and reviewer approval before release.

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper BENCHMARK_CPU_TIME_LIMIT=04:00:00 /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v2/084f512cf89f677fd72b98937980a456f6842b1f/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v2/084f512cf89f677fd72b98937980a456f6842b1f/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-benchmark-bl10-v2/084f512cf89f677fd72b98937980a456f6842b1f/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_benchmark_bl10_20261109_v2 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_benchmark_bl10_20261109_v2 --script src/5_run_benchmark_methods/1_submit_hpc_array.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/benchmark-bl10-inputs-v1/selection.tsv --derivative-release 20261108-confounding-subsets-v1 --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/benchmark-bl10-inputs-v1/input_artifacts.tsv --mem=128G --max-mem=500G --throttle=10
```

### P6 helper implementation accepted

`94cf3e6` checkpoints the four bounded release/evaluation helper and test files.
Both independent R checks pass after review fixes. Explicit canonical paths,
MD5-checked row-specific receipts and evidence, expected/blocked scopes and
new-result requirements prevent discovery/fallback to old results. The
evaluation helper preserves primary NULL's category/count permutation,
records tissue's seed124 labels, binds source/distance/label identities,
rejects invalid LISI separation values and enforces complete finite equally
weighted five-metric/evaluation grids. Notebook integration and production
ledger population are next; no figure refresh is yet claimed.


### BL10 v2 launch event and ongoing assignments

Prepare and reconcile passed: ABSENT remote execution, no discrepancies.
Launched once `2026-11-09T00:54:50Z`, wrapper digest
`0cc285d59729b806d066b062db7a8b61e4540d18389dbce8cf800f676245fc40`.
One unbounded durable waiter is armed (local exec `66518`); do not duplicate
the waiter or relaunch the wrapper. Initial recorded preflight ID `4663549`.
At the one post-launch scope check, the controller had reached strict
downstream source-record validation for Adams and Bassez, demonstrating that
the old home-alias preflight failure is resolved for the representative
legacy source. No final fit/release success is claimed yet.

P6 helpers accepted in `94cf3e6` are now being integrated in parallel:
`batch_subset_audit` owns the batch helper/two maintained batch notebooks;
`docs_checkpoint_audit` owns benchmark notebook/helper plus shared plotting
and constants. They coordinate interfaces and exclude external cancer files.
The release loader may add one narrow `metadata_feather` kind for explicit
checked sample metadata, without a second provenance system. Benchmarks can
use an explicitly selected reviewed non-NULL reference bundle for scalar
primary labels; no eleven new metadata exports are needed just for labels.

`benchmark_dependency_audit` owns a separate two-file/two-test read-only
automatic-annotation reuse evidence entry. The existing three-cohort layer2
differences cannot be generalized to other auto-composition rows without
inspection. A future **metadata-only** durable evidence gate will freeze the
exact three-cohort artifact inventory, input refs and expected checksums,
then use one bounded 64G/1CPU/one-hour job to read obs/compact RDS only and
produce a JSON evidence report. Same `ecoda-benchmark` group after BLv2
terminal review; no new sync group, fits, source mutations or NAS publication.
Unknown annotation mappings/missing producer proof remain blocked/pending;
no additional method rerun is approved from that audit alone. Exact audit
wrapper/source/scope must be recorded after code review and before prepare.


## Execution checkpoint — BL10 terminal failure and screen coverage, 2026-11-09

BL10 v2 evidence: `.gate/ecoda_stage5_benchmark_bl10_20261109_v2.inspect.json`;
first terminal accounting inspection was performed `2026-11-09T02:19:27Z`,
query count exactly one. The 21 successful raw rows cover ten input checks,
one R-environment check and ten fits; the two failed rows are the watchdog
and aggregate only. The existing durable waiter has finished. This failed
gate remains immutable historical evidence and cannot release a dependency.

Root-cause repair is assigned to the reporting/reviewer agent as an urgent,
bounded implementation: `stage5_dispatcher.sh` appends the explicit
`ECODA_ARTIFACT_BINDING` and `ANALYSIS_NAS_ROOT` to the worker runtime export,
but `matrix_watchdog.sh` reconstructs only its base runtime fields before
strict equality. Restore consistent construction without relaxing source,
runtime, input, ownership or retry checks. Test actual runtime helpers and
preserved derivative context. Inspect existing validated-recovery/sync-only
contracts first; propose a validator-only release path without repeating fits
or manually synchronizing partial outputs. No new recovery launch is approved
by this checkpoint alone.

P6 coverage audit additionally found that the maintained cell-depth screen
retained only ANOSIM summaries, and transformation summaries omitted silhouette
and LISI and used a noncanonical modularity setting. Existing summaries cannot
supply the missing scores. The approved five-metric figure scope therefore
requires lightweight composition-only screen refreshes, with the same draw
schedule/seeds and persisted per-draw evidence. Transformation refresh covers
all eleven benchmark datasets; nine additional transformation-only requests
are added to the proposed later wave. Zero-imputation already contains the
five scores and needs only input-dependent Bassez/Stephenson reruns unless a
separate integrity audit proves further invalid rows. None of this authorizes
unrelated RNA-model refits.

The later three-context proposal is currently 63 public benchmark requests,
30 normalized benchmark family rows plus 56 uncorrected and 56 corrected
family rows: **142 family dispatch rows**, not a scheduler-task count.
Physical task/artifact scope must still be frozen before launch. This wave
has not launched. The cell-depth job contract is still under review.


## Three-cohort metadata evidence — exact launch contract, 2026-11-09

State: **FAILED**, finished `2026-11-09T02:42:53Z`; launched once
`2026-11-09T02:41:20Z`. Prepare and reconcile
passed (ABSENT, no discrepancies), and the reviewed Bassez dependency was
verified by the durable launcher. Wrapper digest
`f62dbcafba7437c9ac5cef35ba8391dec629de61f660390fa9564e7070c31cbb`.
The sole durable waiter (local exec `97611`) has finished. Job `4663622`
failed `1:0`; its one terminal accounting inspection is preserved. Do not
relaunch this gate, approve it, or claim a completed comparison report.
This independent evidence gate uses the already reviewer-approved Bassez HR
input and unchanged Zhang/Pelka inputs. BL10 v2's failure is terminal-inspected
and understood; that failed gate remains failed and is not a predecessor.
This audit consumes none of its scientific outputs, writes no shared
`ANALYSIS_ROOT`, and may proceed while BL science validation is repaired.
Retain the existing `ecoda-benchmark` group; do not create another sync group.
Its sole formal data dependency is the reviewed Bassez Stage 2 gate.

- Gate: `ecoda_stage5_auto_annotation_reuse_evidence_20261109_v1`; run: `stage5_auto_annotation_reuse_evidence_20261109_v1`.
- Source: `669532bf01fc23d5602eac41a9f47be77fbca416`, separate snapshot parent
  `20261109-auto-annotation-evidence-v1`; verified archive SHA256
  `8a5c771f006bd4ede3f80e938c47dae31bfa55aa2e9e56cb5c7bf49f7fdbb965`.
- Same pinned format-2 runtime image/manifest and frozen tracked aux as the
  BL10 v2 contract. No environment refresh or model/annotation fitting.
- Exact selection: Bassez, Zhang, Pelka; each `benchmark_analysis` /
  `auto_annotation_reuse_evidence`. Three rows, SHA256
  `477b3c88ebb2164123d081552ba02dc29b2adb7d95106af7ca0b6c42e4828bb7`.
- Exact six-field input refs SHA256
  `315f9ec8b22f9126ef5197f60b255e95f1b08fb0ec553d0256412da5effb66ea`.
  They are a finite subset of the already content-verified BL10 refs.
  Bassez MD5 `678656c738597471d884f88e345e6970`; Zhang
  `7c0854d54bbfdc6b8eb971def8aa7f84`; Pelka
  `f3f2fc8fff93d32a805ccef15d67041b`.
- One ordinary Slurm job, one CPU, 64G, one hour; no array, GPU, watchdog,
  retry, broad discovery, RNA matrix load, fit, or NAS publication.
  Expected exactly **one scheduler ID / one scheduled task**.
- The controller freezes the selected nonrecursive legacy RDS inventory,
  checksum-table bytes and input file stats; the worker consumes only that
  inventory and reads bounded obs/uns plus checked compact RDS artifacts.
  Include author composition only for configured auto-column overlap (Zhang).
- Expected output: new run-owned `artifacts/evidence/auto_annotation_reuse.json`,
  checksum and truthful producer record. Classifications are evidence only;
  `COUNT_MATCH_PROVENANCE_PENDING` does not approve reuse, and no additional
  RNA-method rerun follows automatically from a discrepancy.
- Terminal procedure: one existing durable waiter; one first `inspect` with
  its sole emitted ID; pinned custom `--validate` entry against this run;
  evidence review and second `inspect --approve-reviewer`. No standard fit
  artifact auditor or new scheduler accounting query.

Immutable wrapper:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v1/669532bf01fc23d5602eac41a9f47be77fbca416/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v1/669532bf01fc23d5602eac41a9f47be77fbca416/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v1/669532bf01fc23d5602eac41a9f47be77fbca416/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_auto_annotation_reuse_evidence_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_auto_annotation_reuse_evidence_20261109_v1 --script src/utils/bash/auto_annotation_reuse_evidence.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v1/selection.tsv --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v1/input_artifacts.tsv --results-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/benchmark/results --derivative-release 20261108-confounding-subsets-v1
```

The initial snapshot command stopped before creating a snapshot because the
new parent directory was absent. After explicitly creating that new parent,
source creation and remote control/archive hashes passed. No job had been
submitted during that setup correction.

Reviewed implementation checkpoints `d29551d` (five-metric transformation
screen) and `669532b` (shared derivative runtime export) are pushed to master.
The production BL-cache validator's extra-field mismatch was also reproduced
read-only on the actual Adams output; its isolated compatibility fix is in
progress. No completed scientific fit was repeated.


### Metadata audit integrity stop and narrowed recovery — 2026-11-09

`ecoda_stage5_auto_annotation_reuse_evidence_20261109_v1` stopped at the first
legacy RDS/table checksum disagreement, before deserializing that file.
Its sole job `4663622`, worker stderr, partial 46-byte temporary report,
original frozen inventory, and first terminal inspection remain preserved.
No method fit or source mutation occurred; no valid final audit report exists.

A separate bounded no-compute hash/metadata investigation inspected exactly
that 74-file frozen inventory, without deserializing unverified payloads:
Bassez 24 files, Zhang 26, Pelka 24. Only one file per cohort matches the
legacy checksum table. All 74 per-file sidecars match actual bytes, but their
historical PATHs use the home/scratch alias and no current canonical-path
owner record was found. This is conflicting legacy integrity/provenance
metadata, not permission to rewrite old checksums or approve reuse.
Small diagnostic evidence is in
`.gate/auto_annotation_legacy_checksum_proofs_20261109.json`.

The three table-valid files are the exact `<Dataset>_composition.rds`
aggregates. Checksum-before-read inspection with the configured Bamboo Pixi R
interpreter confirms each has 53 result keys, including all 23 automatic
composition keys; Zhang also contains its two configured author/auto-column
overlap keys. This covers the required **69 automatic bundle comparisons plus
two Zhang dependency comparisons** without reading the 71 conflicting
standalone files. Key/dimension evidence is
`.gate/auto_annotation_verified_aggregate_keys_20261109.json`.

A narrow existing-entry extension is assigned: explicit finite
`--legacy-artifacts TSV` physical scope, with strict dataset/path coverage and
full checksum-table validation before submitting any job. A fresh recovery
scope may select only the three verified aggregate files; preserve the old
default and all failed-run evidence. Do not add a checksum fallback, silently
reconcile sidecars, mutate a legacy table, or infer that excluded standalones
are reusable. The future recovery needs its own committed source, new run,
exact wrapper and recorded scope before launch. No additional RNA-method
rerun is authorized by this integrity finding.

Local implementation checkpoints now include `032fc7a` (cell-depth draw
retention and all affected supplemental five-metric report code) and
`2847032` (strict current BL-cache/result validation). Root focused real CLI,
cell-depth and supplemental report tests pass. A legacy transformation fixture
in the broader RDS contract suite is being updated to the committed v2
contract; its old summary-only success expectation is no longer valid.

One local test initially overwrote the untracked scratch `Rplots.pdf`; the
test now opens an explicit temporary device and verifies the scratch file is
unchanged before/after. No maintained publication figure changed. The prior
scratch bytes could not be restored; the file remains excluded from commits.


## Stage 3 terminal evidence and validator-only reuse — 2026-11-09

The immutable source remains `7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3`,
producer run `stage3_confounds_20261108_v3`. All 17 selected outputs and their
sidecars were synchronized by the original wrapper to its bound NAS release;
total H5AD/sidecar size reported by rsync is 492,399,725,768 bytes. No manual
sync or source overwrite followed. Terminal producer/watchdog status is OK,
retry index 2, final array `4663537`.

Exactly one terminal inspection queried roots `4663496`, `4663499`,
`4663522`, `4663537`, `4663500`; full 24-row capture is in
`.gate/ecoda_stage3_confounds_20261108_v4.inspect.json`. Initial Alzheimer
attempts `4663501`/`4663502` and retry-1 `4663523`/`4663522` report
`OUT_OF_MEMORY|0:125`. The four actual stderr files explicitly identify those
StepIds, including the reversed numeric order in retry 1; do not infer task
mapping from consecutive numeric IDs. Both retry selection TSVs contain only
the same two Alzheimer views. Final retry rows `4663537`/`4663538` and
watchdog `4663500` are `COMPLETED|0:0`; all remaining captured rows complete.

The global failed gate and every OOM row remain immutable evidence. The
checked-in durable ECODA skill's validator-first recovery rule (§6) requires
valid artifact rows to remain outside recomputation: when every selected
artifact validates, write `NOOP_VALIDATED` without gate prepare/launch. Use
that already-authorized no-compute path to review individually produced
artifacts; it does not approve/relabel the old failed gate or turn it into a
passing dependency-manifest. Future jobs bind reviewed exact artifact/proof
identities, never a fictional successful predecessor gate.

Assigned implementation is a narrow Stage3-only reuse-report path in the
existing run auditor, retaining its standard mode unchanged. It consumes the
original failed gate/terminal inspection and complete retry lineage, validates
source/runtime/config/full-run identities, then performs expensive ownership,
checksum, cohort/schema, NAS equality and NAS-mode checks only for an explicit
nonempty unique subset of the original 17 rows. Report creation is separate,
create-only and read-only with respect to producer outputs/status/proofs.
Reject any non-OOM/unresolved failure, missing proof or foreign scope.

First validate the one Stephenson benchmark input (expected 53 samples) so
its independent annotation can proceed after review. Validate the remaining
16 batch outputs separately, avoiding a large NAS hashing prerequisite for
that one input. Eventually reconcile complete 17-row coverage. No additional
preprocessing is selected. Exact validator source/report paths and evidence
will be recorded before execution/release.

## Narrowed metadata evidence v2 — exact launch contract, 2026-11-09

State: **PLANNED, ready for prepare/reconcile/launch after this checkpoint**.
The old v1 evidence gate stays failed; it is not a dependency. The reviewed
Bassez metadata gate is the sole formal source dependency. Same existing
`ecoda-benchmark` group, new root and separate snapshot parent; no shared
analysis output, fits, RNA-array reads, or NAS publication.

- Gate: `ecoda_stage5_auto_annotation_reuse_evidence_20261109_v2`; run: `stage5_auto_annotation_reuse_evidence_20261109_v2`.
- Source `b03cdd8d4eeaeb72a3b58cdde1b2aa5845911b17`, snapshot parent
  `20261109-auto-annotation-evidence-v2`; archive SHA256
  `74b3f99a779ec530b417cf0c1013d7b775d2de04e6b3f92092ee132e24393347`.
- Same immutable format-2 runtime/image/aux identities as v1. Exact three
  cohort selection/input-reference hashes remain
  `477b3c88ebb2164123d081552ba02dc29b2adb7d95106af7ca0b6c42e4828bb7` /
  `315f9ec8b22f9126ef5197f60b255e95f1b08fb0ec553d0256412da5effb66ea`.
- Explicit physical artifact TSV SHA256
  `e4caa7e50e560763af53e740d93708141217d45460f38fcbdbdf549bcc402ec4`:
  exactly Bassez/Zhang/Pelka `<Dataset>_composition.rds` beneath canonical
  scratch `benchmark/results`. All three match the preserved checksum table.
  The 71 conflicting standalone files remain outside this scope and unread.
- Expected 69 automatic composition comparisons plus two Zhang configured
  author/auto-column overlap comparisons, **71 total**, from three checked
  aggregate payloads. Counts/provenance findings do not authorize reruns.
- One Slurm CPU job, 64G, one hour; exactly one ID/task; no array, GPU,
  watchdog or retry. Controller hashes selected payloads before submission;
  scope, source/runtime, checksum table and physical request are frozen.
- Terminal procedure: exactly one durable waiter, one accounting inspection
  with the sole emitted ID, pinned custom `--validate`, evidence review,
  then reviewer approval only if all gate/audit requirements pass.

Exact immutable wrapper:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v2/b03cdd8d4eeaeb72a3b58cdde1b2aa5845911b17/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v2/b03cdd8d4eeaeb72a3b58cdde1b2aa5845911b17/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v2/b03cdd8d4eeaeb72a3b58cdde1b2aa5845911b17/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_auto_annotation_reuse_evidence_20261109_v2 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_auto_annotation_reuse_evidence_20261109_v2 --script src/utils/bash/auto_annotation_reuse_evidence.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v2/selection.tsv --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v2/input_artifacts.tsv --results-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/benchmark/results --derivative-release 20261108-confounding-subsets-v1 --legacy-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v2/legacy_artifacts.tsv
```

Source and all three control hashes were independently checked on Bamboo.
Reviewed local commits are pushed through `b03cdd8`; includes `032fc7a`,
`2847032` and test-only `40fec4f`. The full RDS contract suite now passes its
current transformation fixture, retaining old-summary rejection.

The proposed main3context physical dry-run audit found a separate legacy
corrected repair guard in `stage5_pending.sh` that blocks new derivative
MRVI/PILOT/QOT rows before any submission. A minimal exact-bound fresh-run
exception is assigned, preserving historical residual restrictions. Until
that fix and actual corrected scope test pass, corrected physical counts are
intended estimates only. The main method wave remains unlaunched.


## Validator execution checkpoint — 2026-11-09, fc2949f

All reviewed code is committed and pushed through
`fc2949f7deff0b4a1aea6f7751b0226518643bfa`: includes corrected fresh-wave
`e45cc0c`, Stage3 subset validator `3d3eb2e`, owned BL review `da8081e`,
and explicit evidence JSON-file output `fc2949f`. Root actual-entrypoint,
strict-CLI, ownership/runtime, metadata-channel tests pass; independent review
accepted the corrected-wave guard/owner fix. The corrected dry-run now proves
56 tasks/seven arrays/seven watchdogs/one aggregate, plus smaller fresh cohort
support with canonical old owners preserved. No main method wave launched.

Three separate immutable snapshot parents were created on Bamboo for this
same full commit: `20261109-stage3-noop-validation-v1`,
`20261109-bl-science-validation-v1`, and
`20261109-auto-annotation-evidence-v3`. Each archive SHA256 is
`6cfe3f1399632d5989965062725f4797398020d78ab5bcde254ffbc9d87cad84`.
Separate parents prevent independent executor locks from serializing work.
The canonical Bamboo checkout is clean at this exact revision; existing
producer snapshots/data/results are unchanged.

Metadata evidence v2 launched once `2026-11-09T03:10:46Z`, wrapper digest
`2b8c05072860d64ece93e28be4d56e0f50ba40e794c100914e0b2a4e8735e90d`,
sole waiter exec `6520` finished **FAILED** at `2026-11-09T03:12:42Z`.
Job `4663649` was terminal-inspected exactly once. All three comparisons ran
and a 236,163-byte report plus sidecar exists, but R startup stdout prefixed
`WARNING: ignoring environment value of R_HOME` to its JSON and final parsing
failed. Preserve that output/status/failed gate; do not strip/repair it in
place or claim release. The committed narrow fix writes JSON directly to a
create-only file and keeps stdout in the log. A fresh v3 run may repeat only
this failed small metadata evidence job, never RNA-method fits.

BL10's exact 20 saved cache/result files **pass the corrected strict RDS CLI**
on Bamboo from committed `b03cdd8`. This is a successful preliminary contract
check, not the final coupled source/policy/timing/ownership review or NAS
publication. The checked-in validator-only BL receipt will create a separate
owned run and preserve the original failed gate/all 23 accounting rows.

### Stage3 exact validator-only scopes (no scheduler work)

New control root:
`/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-reuse-evidence-v1`.
Original 17-row selection SHA256
`e91880b664163ca6480549f2c5ad29839ce7b95956aef9aa1de08bc48dabe54d`;
Stephenson benchmark one-row scope
`641deda51c98f11590693ca0f3b552a96f565f35b76d3b6817ad73abd3378bf3`;
remaining 16 batch-view scope
`bd594bb23ed34ce83593d523a8214a584cf7e26b815d8d4476f8f98c2a05ebc0`.
Preserved original failed gate/inspect copies have SHA256
`68d7cd6a90872b0b7869d533c307943656527eda98fd51b346301138bb9b062d` /
`4470c08d9b8dcb685573ea921020481a1c42479fa4a98452143cbc45066ea62b`.
Only these small text controls were transferred locally/remotely.

The new immutable auditor is invoked directly (the documented no-compute
validator path), with original producer7e source/runtime/run identities.
It does not execute a newsource snapshot wrapper that could override original
source environment. Clear inherited run/source/runtime/binding/analysis vars;
set canonical NAS and configured Pixi Python explicitly. Its receipt records
actual validator code hashes separately from the complete old producer
archive/config/runtime identity. No new Slurm query or submission is permitted.

Exact commands are recorded in
`.gate/stage3_noop_controls_20261109_v1/commands.json`; root executes only its
named `stephenson_benchmark` then `batch_views` entries. Reports are create-only
`reports/stephenson_benchmark.json` and `reports/batch_views.json` under the
control root, with strict checksum sidecars. They retain the failed gate,
all 24 accounting rows and four StepId OOM witnesses. They fully hash selected
scratch/NAS outputs and verify schemas, cohort/obs/source identity, owners,
producer records, NAS equality and 666/777 modes. Historical raw production
MD5 references/proofs and stable stats are reused; no fresh raw-input content
hash or model fit is claimed. Root must review actual sample counts against
53 and the approved batch counts before any dependent launch.

First exact command:

```bash
/usr/bin/env -i PATH=/usr/local/bin:/usr/bin:/bin HOME=/home/users/h/halterc USER=halterc LOGNAME=halterc PYTHON_BIN=/home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13/bin/python NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 ECODA_GATE_STAGE=stage3 /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-stage3-noop-validation-v1/fc2949f7deff0b4a1aea6f7751b0226518643bfa/tree/src/utils/bash/ecoda_run_audit.sh --run-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage3_confounds_20261108_v3 --stage stage3 --selection /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage3_confounds_20261108_v3/manifests/selection.tsv --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261108-confounding-subsets-recovery-v4/7e67ef71ed317c4d81f94d6bfc5b0252f2776ef3/identity/source.manifest --runtime-identity /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage3_confounds_20261108_v3/manifests/runtime.identity --stage3-reuse-report /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-reuse-evidence-v1/reports/stephenson_benchmark.json --reuse-selection /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-reuse-evidence-v1/stephenson_benchmark.tsv --terminal-inspect /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-reuse-evidence-v1/original_inspect.json --failed-gate-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage3-reuse-evidence-v1/original_gate.json
```

Command-record SHA256: `f9505470869b0cebc2ba43ce377255bab4c1be3f487e78997af349412d332fe3`.


### BL10 validator-only execution contract — 2026-11-09

State: **READY for synchronous no-compute execution**. The checked-in durable
profile's validator-first NOOP rule applies: no gate prepare/launch, worker,
Slurm query, model fit or publication is performed. New source/runtime/control
and report ownership are bound through standard helpers in a fresh owned run;
the actual execution is controller Pixi R, recorded truthfully in the receipt.
This is independent of the Stage3 artifact/NAS checks and shares no analysis
or sync root with them.

- Validation run: `stage5_bl10_science_validation_20261109_v1`.
- Source `fc2949f7deff0b4a1aea6f7751b0226518643bfa`, separate
  `20261109-bl-science-validation-v1` snapshot; archive SHA256
  `6cfe3f1399632d5989965062725f4797398020d78ab5bcde254ffbc9d87cad84`.
- External control: canonical `_ecoda_inputs/20261108-confounding-subsets-v1/bl-science-validation-v1/control.json`;
  MD5 `f946b3f31a4d00002e031505620fa02a`, SHA256
  `5aa9af1621ea754d3da455971ed9883f4a86b5254391b3fa11b79b46e21d5713`, 19,999 bytes.
  Collector SHA256 `320279c81f3b2c75eb031c9192f9115d4845811d3c5399a29d60eea5d674b3ec`.
- Exact physical ten-row matrix order: Adams, Bassez,
  Gongsharma_cmv_young_males, Kfoury, Kim, Lee, Pelka, Smillie, Wu, Zhang.
  Only benchmark `pseudobulk/hvg2000_bl`; no Stephenson or unrelated variants.
- Exact existing evidence: 20 RDS caches/results, 10 scientific timing
  Feathers, original successful worker logs and completed preflight proofs,
  old source/runtime/aux/config/input identities and current observed stats.
  No fresh whole-input H5AD hash/consumption is claimed.
- All 23 raw accounting rows are preserved: 21 successful rows and original
  watchdog4663573 / aggregate4663574 failures. Original failed gate/profile/
  inspection/raw-capture copies are checksum-bound; no success-only projection
  or new accounting query. Original gates/owners/statuses/bytes stay unchanged.
- Output: `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage5_bl10_science_validation_20261109_v1/artifacts/evidence/bl_science.json`, create-only JSON/checksum/truthful
  producer record/global owner. Receipt status `NOOP_VALIDATED`, scientific
  validation only, original gate FAILED, publication explicitly pending.

Exact immutable execution command:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bl-science-validation-v1/fc2949f7deff0b4a1aea6f7751b0226518643bfa/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bl-science-validation-v1/fc2949f7deff0b4a1aea6f7751b0226518643bfa/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bl-science-validation-v1/fc2949f7deff0b4a1aea6f7751b0226518643bfa/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_bl10_science_validation_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_bl10_science_validation_20261109_v1 --script src/utils/bash/validate_bl_science.sh -- --control /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/bl-science-validation-v1/control.json --control-md5 f946b3f31a4d00002e031505620fa02a --report /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage5_bl10_science_validation_20261109_v1/artifacts/evidence/bl_science.json
```

### Selected NAS permission repair before Stage3 reuse — 2026-11-09

The first one-row validator reached NAS payload/hash checks but stopped before
creating a receipt because two newly created selected NAS namespace directories
had mode0700. This is a producer publication-permission defect, not evidence
that H5AD bytes changed. `datasets/Stephenson` and `output` were777, selected
files/proofs666, release ancestors777. Root verified the immutable binding
and changed **only** the bound Stage3 NAS run directory and its `datasets`
directory from0700 to0777. No recursive chmod, scratch/home permission change,
resync, model fit, or payload/old-proof rewrite occurred.

Create-only before/after repair evidence:
`_ecoda_inputs/20261108-confounding-subsets-v1/stage3-reuse-evidence-v1/reports/nas_directory_permission_repair.json`,
MD5 `7035a30325f95acc7ada564e3a65f00a`.
The same one-row read-only validation was restarted after that repair
(local exec `67613`); the initial failed validation exec `98826` produced no
receipt. A minimal future Stage3 sync ancestor-mode fix is assigned with
focused tests. All old producer and durable accounting statuses remain intact.


## Reviewed artifact readiness and next independent launches — 2026-11-09

**Stephenson Stage3 benchmark input: individually validated and reviewed for
reuse.** Exact 53 retained samples, all Ncl, 43 Covid /10 Healthy; all four
excluded IDs absent. Correct cohort source/config/obs identity, raw-count/HVG/
PCA schema, full scratch/NAS payload hashes/equality, producer/ownership and
NAS modes pass. H5AD MD5 `657c72165dfbc4dd1caf83fc656e9c39`, size
13,333,003,144. NOOP receipt MD5 `0e793f1f1bb1bb60c8257c82a3edf47e`, SHA256
`27266a3613315905d9905a74b5d23fd1eba8908943f8104496c1a3e318b7bd21`;
small local copy `.gate/stage3_noop_controls_20261109_v1/stephenson_benchmark.json`.
This approves this specific artifact; original durable gate remains FAILED.
Remaining16 batch artifacts are being validated independently (exec `52786`).

**BL10 science: validated and reviewed; publication still pending.** The new
validator run `stage5_bl10_science_validation_20261109_v1` finished
`NOOP_VALIDATED`. All20RDS/10timing files, exact features/sample identities,
five-group pre-fit policy/audits/timing, original preflights/source/runtime
and preserved23 accounting rows pass. Every dataset has2000features; sample
counts Adams100/Bassez26/Gongsharma180/Kfoury32/Kim58/Lee27/Pelka98/Smillie64/
Wu26/Zhang31. All use unchanged aux MD5 `8cacee7893cbb9d170cc9de87e3bd306`.
Receipt MD5 `22dae0a1b19a0e0545c309f1308ba2e8`, SHA256
`acaab097e284bea392e0f79afe8156c8eaa1f555753056593a7460e95f29a3ae`.
Root review `.gate/bl_science_controls_20261109_v1/review.json` approves
scientific reuse only: old gate FAILED, NAS publication false. No fit repeated.
Stephenson BL remains solely in its eventual whole-method benchmark wave.

### Stage4 Stephenson annotation — RUNNING, launched 2026-11-09 04:01:47Z

- Gate `ecoda_stage4_stephenson_annotation_20261109_v1`; run `stage4_stephenson_annotation_20261109_v1`;
  group `ecoda-annotation`, independent of metadata evidence/Stage3 checks.
- Exact source `fc2949f7deff0b4a1aea6f7751b0226518643bfa`, separate snapshotparent
  `20261109-stephenson-annotation-v1`; archive SHA256
  `6cfe3f1399632d5989965062725f4797398020d78ab5bcde254ffbc9d87cad84`.
- Exact one-row selection SHA256 `641deda51c98f11590693ca0f3b552a96f565f35b76d3b6817ad73abd3378bf3`;
  six-field source reference SHA256 `550d2791571b00713d90e0b303958fcb9098e4f698de520c185392c3d659ab82`.
  Bound source producer `stage3_confounds_20261108_v3`, reviewed input above.
- No failed predecessor gate is passed as a dependency-manifest. Helper binds
  and verifies the reviewed NOOP receipt; exact producer/input/proof/MD5/size
  are in the immutable wrapper/control and new run-owned source binding.
- Same pinned8fcd... runtime/1c4c...manifest/frozenaux. Four reference maps
  and scGate DB/cache validation pass on the exact source. Initial standalone
  reference command needed parent-shell slurm_config exported for the second
  validator; corrected read-only invocation passed, with no regeneration.
- Fresh scratch/log/NAS/gate namespaces verified. Source author metadata,
  raw/counts/features/embeddings/graphs/order remain protected; annotation
  writes only a fresh run-owned copy. Registry Blood/normal_tissue=true;
  dual HiTME layers1–3 and scATOMIC1–6, predictions/confidence/cell-cycle.
- Expected53samples /two samples per chunk →27 annotation tasks; preparation1,
  merge1, watchdogs3: **6 initial scheduler IDs /32 scheduled tasks**.
  CPU4 workers,32Ginitial/128Gmax, shared-cpu, throttle27; prep1h,
  annotation+merge2h. CPU1 watchdogs2G/2G/32G,12h. OOM-only bounded retries
  may add IDs/tasks and must remain explicit terminal evidence.
- Preserve legacy checkpoints; no reuse from unproven old annotation sources.
  Exactly one durable waiter; one full accounting inspection and pinned
  selected artifact/copy/source/NAS review before downstream use.

Exact immutable wrapper:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper STAGE4_INPUT_PRODUCER_RUN_ID=stage3_confounds_20261108_v3 ANNOTATION_WORKER_TIME_LIMIT=02:00:00 ANNOTATION_WATCHDOG_TIME_LIMIT=12:00:00 ANNOTATION_MERGE_WATCHDOG_MEMORY=32G /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-stephenson-annotation-v1/fc2949f7deff0b4a1aea6f7751b0226518643bfa/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-stephenson-annotation-v1/fc2949f7deff0b4a1aea6f7751b0226518643bfa/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-stephenson-annotation-v1/fc2949f7deff0b4a1aea6f7751b0226518643bfa/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage4_stephenson_annotation_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage4_stephenson_annotation_20261109_v1 --script src/4_cell_type_annotation/1_submit_onboarding_stage.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage4-stephenson-annotation-v1/selection.tsv --derivative-release 20261108-confounding-subsets-v1 --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/stage4-stephenson-annotation-v1/input_artifacts.tsv --mem=32G --max-mem=128G --partition=shared-cpu --throttle=27
```

### Metadata evidence v3 — COMPLETED and reviewed 2026-11-09 04:11:43Z

Gate `ecoda_stage5_auto_annotation_reuse_evidence_20261109_v3`, run `stage5_auto_annotation_reuse_evidence_20261109_v3`.
Same `ecoda-benchmark` group and sole reviewed Bassez data dependency as v2;
no dependency on a failed metadata/BL gate. Sourcefc2949f, existing separate
`20261109-auto-annotation-evidence-v3` snapshot/archive6cfe... .
Exact unchanged selection/input/physical-inventory hashes:
`477b3c88ebb2164123d081552ba02dc29b2adb7d95106af7ca0b6c42e4828bb7` /
`315f9ec8b22f9126ef5197f60b255e95f1b08fb0ec553d0256412da5effb66ea` /
`e4caa7e50e560763af53e740d93708141217d45460f38fcbdbdf549bcc402ec4`.
Three verified aggregate RDS files, Bassez/Zhang/Pelka,71 comparisons; the71
conflicting standalones remain excluded/unread. One CPU64G1h job, exactly1ID/
1task; no fit, array, GPU, watchdog, retry, sharedanalysisroot or NASpublication.
Only failed metadata-evidence generation repeats. R writes JSON directly to
an exclusive file; actual startup-noise/create-only regressions pass.
One durable wait, one accounting inspection, pinned custom --validate, then
reviewer approval if all contracts pass.

Exact immutable wrapper:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v3/fc2949f7deff0b4a1aea6f7751b0226518643bfa/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v3/fc2949f7deff0b4a1aea6f7751b0226518643bfa/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-auto-annotation-evidence-v3/fc2949f7deff0b4a1aea6f7751b0226518643bfa/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_auto_annotation_reuse_evidence_20261109_v3 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_auto_annotation_reuse_evidence_20261109_v3 --script src/utils/bash/auto_annotation_reuse_evidence.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v3/selection.tsv --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v3/input_artifacts.tsv --results-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/benchmark/results --derivative-release 20261108-confounding-subsets-v1 --legacy-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/auto-annotation-evidence-v3/legacy_artifacts.tsv
```

FutureStage3NAS ancestor-permission fix is independently tested/committed
`29d202e`; it normalizes/verifies only selected NAS path ancestors, preserves
NASroot/unrelatedmetadata/scratch, and fails if chmod or verification fails.
The active scientific producer and validator snapshots remain unchanged.


## Superseded live checkpoint preserved — before 2026-11-09 04:23 UTC


- P0/P1 complete. Reviewed code pushed through `084f512`; `0514b48`
  adds the exact HR selectors and `a9237c1` adds explicit variant requests,
  independent campaign contexts and strict derivative NAS/OOM checks.
  External cancer application files remain excluded from these commits.
- Stage 3 v4 wrapper **COMPLETED `0` at `2026-11-09T02:53:27Z`**:
  all 17 new H5ADs were produced and synchronized to NAS. The strict global
  gate is nevertheless **FAILED** after its sole terminal inspection because
  it retains historical OOM attempts. Alzheimer uncorrected/corrected alone
  retried at 256G, then succeeded at 500G; the other 15 rows completed.
  All 24 raw accounting rows and five scheduler roots are preserved, including
  four OOM attempts, final retry success and the successful watchdog.
  No gate approval or downstream artifact release is claimed yet. Validate
  existing outputs for explicit NOOP reuse; do not repeat successful processing.
- Bassez gate `ecoda_stage2_bassez_hr_20261108_v1` is **COMPLETED,
  audited and reviewer-approved**. Worker `4663493` and watchdog `4663494`
  completed `0:0`; the one accounting inspection, accounting bridge, pinned
  Stage 2 source/artifact/NAS audit and independent evidence review passed.
  Formal approval is recorded in `.gate/ecoda_stage2_bassez_hr_20261108_v1.review.json`.
  The new derivative filled 50,809 HR labels (33 to 40 categories), retaining
  74,948 cells and 26 samples. All protected data are unchanged. Its MD5 is
  `678656c738597471d884f88e345e6970`; the legacy input remains unchanged.
- Stage 5 early BL10 v1 is **FAILED**, terminal-inspected once. Its sole
  preflight array `4663527` passed Bassez but failed the nine legacy inputs;
  no R environment/method jobs launched and **zero scientific RDS outputs**
  exist. The controller finished `2026-11-09T00:22:47Z`. Do not relaunch it.
  The repaired Bassez input remains independently approved. A narrow explicit
  bound-proof fix is underway for compute-side legacy checksum/record checks;
  no old input or receipt is rewritten. A fresh BL gate is required after
  review. Stephenson BL remains solely in its later full-method wave.
- BL10 recovery v2 is **FAILED at the watchdog/release layer**, terminal-inspected
  once. It launched `2026-11-09T00:54:50Z` from immutable `084f512` and
  finished `2026-11-09T01:31:21Z`. Both preflights (`4663549`, `4663561`)
  and all ten fit workers in array `4663563` completed `0:0`. The expected
  ten BL pseudobulk caches and ten standalone result RDS files exist.
  Watchdog `4663573` rejected the augmented derivative runtime export;
  aggregate `4663574` consequently failed. First inspection preserves all
  23 raw accounting rows, including both failures; no accounting projection,
  gate approval or NAS publication is claimed. Preserve these scientific
  outputs and validate their exact contracts before choosing no-fit recovery.
  Do not rerun successful fits, edit old statuses or hide failed accounting.
- Reviewed local commits after `084f512`: `94cf3e6` release/evaluation helpers,
  `33fa9c9` five-metric benchmark/batch reporting, `c72eb4d` harmonized operator
  documentation, `7719d4d` metadata-only annotation-reuse evidence, and
  `f91a2fe` focused physical GM dispatch coverage. These have not changed any
  running immutable snapshot. Supplemental screen updates remain in progress.
- Earlier Stage 3 v3 failed before producing new H5ADs because compute nodes
  cannot mount NAS (17 missing / 0 existing outputs). Its evidence is retained;
  the reviewed fix validates scratch on compute and publishes NAS from the
  controller. The current recovery uses fresh roots without force.
- The multi-child accounting adapter is accepted and committed in `eb26bdb`;
  12 tests pass, including seven unchanged single-child cases. It preserves
  one raw capture and derives three exact child projections before three
  pinned audits. No ordinary parent audit or extra accounting query.


## Terminal metadata evidence and HiTME closure finding — 2026-11-09

Metadata v3 exact sourcefc2949f, command digest
`a990f45608a0865271cee3c474195ca37892370cad0e7e34fe00cb6e048ecaf8`.
The sole accounting inspection for job4663686 passed; the pinned custom
validator passed after explicitly binding both scratch-root environment
variables to the canonical scratch path. Reviewer approval is recorded in
`.gate/ecoda_stage5_auto_annotation_reuse_evidence_20261109_v3.review.json`.
Receipt `.gate/auto_annotation_reuse_inputs_20261109_v3/auto_annotation_reuse.json`:
MD5 `679af625528283dcf105e6d15387e163`, SHA256
`6f9adb904a1b3f1550d9441e203f80f750add3d8d0e76c1427e441465ec29219`.
Three verified aggregates: Bassez23 comparisons, Zhang25, Pelka23.
41 COUNT_MATCH_PROVENANCE_PENDING and30 INVALID_CURRENT_COUNTS. Each cohort's
10 invalid rows are exactly `ECODA_HiTME_HR_layer2_top_varexp{0,...,0.9}`.
Plain HiTME L2/L3 and scATOMIC baselines match; Zhang author LR/HR also match.
No inference is made from excluded checksum-conflicted standalone artifacts.

Actual delayed composition producer reproduces the defect: the outer
`hitme_ct` loop advances to layer3 before deferred layer2 methods execute.
The minimal fix captures the layer and forces the variance threshold within
that iteration's callback environment. The focused real counts/CLR/HVC test
checks all10 L2 and10 L3 variants, exact counts/features/distances/labels,
checksum-validated saved RDS, and unchanged author/plain baselines. Independent
review and selected-author-HR regression pass. No RNA/model fit was run.

Proposed additional scope, pending user decision: ten exact L2 variants across
Adams,Bassez,Gongsharma_cmv_young_males,Kfoury,Kim,Lee,Pelka,Smillie,Wu,Zhang
(100 public parameter rows). Stephenson is already covered by its full rerun.
Use existing composition metadata-only workers and selected-variant plumbing;
no preprocessing, embedding or RNA-model recomputation is justified. Keep
proposed scope separate from the already approved Bassez author-HR22 early
wave while its exact physical contract is prepared.


Reviewed HiTME closure implementation checkpoint: `e776105`; independent reviewer
reproduced the old failure and passed the current actual-builder regression.
Future Stage3 NAS ancestor fix: `29d202e`. These commits are source fixes only;
active gates remain pinned to their original immutablefc2949f snapshots.


## Bassez author-HR22 independent early wave — FAILED, terminal-inspected 2026-11-09

Purpose: run the already approved repaired-author-HR composition results now,
while independent Stephenson annotation and Stage3 batch artifact validation
continue. No new HiTME refresh scope is selected. This is one selected
composition worker, not a whole-composition or RNA-method rerun.

- Gate `ecoda_stage5_bassez_author_hr22_20261109_v1`; run
  `stage5_bassez_author_hr22_20261109_v1`; group `ecoda-benchmark`.
  Sole formal dependency: approved Bassez Stage2 gate. All same-group gates
  are terminal; metadata v3 review is approved. Independent Stage4 uses its
  existing separate annotation group and snapshot.
- Source `129c1a46ee7d282773eb090448ed3d061ac8ecf5`, separate snapshotparent
  `20261109-bassez-author-hr22-v2`, archive SHA256
  `1b9a587450628b6082591aa86279b8544a4c12d8604148f232c5a2d685a2b0f8`.
  Earlier cc5e470/v1 source snapshot was created but never executed; do not
  delete/replace it or confuse it with this launch. Existing fc active runs
  stay unchanged. Same pinned format2 runtime and unchanged frozen auxiliary.
- Source input is the approved Bassez derivative, MD5
  `678656c738597471d884f88e345e6970`, size3,301,360,648. Exact six-field
  reference SHA256 `e972a6c154ae37420c15fb4ea26fa20f5d8fbf5af43950faacbcabd18940eb30`.
- Public22-row selection SHA256
  `c6d219bb9c9b99e369467ab79509dbacd9d9cc23a58a4ee2225a5ecb155f23ee`;
  actual normalized one-row composition CSV SHA256
  `6320d84831019444bd909648e9cdb77f5e336352199270daa71d9b9e8693532a`.
  Controls are staged/hashed/sealed under
  `_ecoda_inputs/20261108-confounding-subsets-v1/bassez-author-hr22-v1`.
- Exact22 public methods: authorHR baseline and NULL, GloProp, Freq_highres,
  ten variance thresholds0..0.9, 3most/2least/3least-variable CT settings,
  and composition PCA2/3/5/10/15. No HiTME, pseudobulk preparation, RNA PCA,
  Python/GPU, or unrelated composition variant. Output23 RDS: metadata1+
  standalone22, sidecars and ordinary timing/runtime/owner/log records.
- Physical dispatcher stub passed: method array1row +watchdog1+aggregate1;
  production adds H5AD preflight1 and Renv preflight1. Expected initial
  **5 IDs /5 scheduled task units**. OOM-only bounded retries, if any, remain
  explicit. Method CPU16/shared-cpu/EPYC-7742,128G→500G ceiling,4h,GPU0;
  effective concurrency1 (existing method throttle literal1000 over one row).
- Fresh run/log/NAS/gate namespaces verified before prepare. New Stage5 NAS
  selected-ancestor fix129c1a4 is independently reviewed; normal, failedchmod,
  ineffectivechmod and selected-transfer-failure fixtures all pass. NASroot,
  scratch and unselected metadata are unchanged. A broader legacy canonical
  retry fixture fails identically against the unchanged pre-fix helper
  (`post-log transfer failure unexpectedly succeeded`); this separate known
  failure is preserved in `.gate/stage5_nas_permission_tests_20261109_v1/`.
  No legacy behavior is weakened to make that test pass.
- Dedup is conditional on successful terminal audit/review of this wave.
  Only then remove exactly these22 from the later63-row proposal: remaining
  41 public /29 normalized families, SHA256
  `56119d971ed5d1bb2581c2765a1b88c112c9fb8d732e3a387491f22320cbbc79`.
  Bassez's other20 affected requests remain, as do all Stephenson/all-nine
  additional transformation requests. A pending user decision on100 HiTME
  variants stays separate. Release ledger consumes the early reviewed root;
  do not recompute/copy these22 into the main method root.
- Exactly one durable waiter after launch, one terminal accounting inspection
  with every emitted ID, selected pinned artifact/source/NAS audit, then root
  reviewer approval. No completion/reuse is claimed at readiness.

Exact immutable wrapper (command SHA256
`6fc8dc0e04179be474693dd018eb9e7e650168e83f9de36339907cd50bf3f0a7`):

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper BENCHMARK_CPU_TIME_LIMIT=04:00:00 /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bassez-author-hr22-v2/129c1a46ee7d282773eb090448ed3d061ac8ecf5/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bassez-author-hr22-v2/129c1a46ee7d282773eb090448ed3d061ac8ecf5/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bassez-author-hr22-v2/129c1a46ee7d282773eb090448ed3d061ac8ecf5/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_bassez_author_hr22_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_bassez_author_hr22_20261109_v1 --script src/5_run_benchmark_methods/1_submit_hpc_array.sh -- --selection-file /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/bassez-author-hr22-v1/selection.tsv --derivative-release 20261108-confounding-subsets-v1 --input-artifacts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/bassez-author-hr22-v1/input_artifacts.tsv --mem=128G --max-mem=500G --throttle=1
```

Small Kidney secondary-label audit helper is staged read-only at the Stage3
reuse controlroot, SHA256
`5c6393f9d49997fcad74b3371af60221213f0c9937fc16e6b6e0086c1a3c797a`.
It will run after the batch NOOP receipt exists, uses configured specimen and
its producer-proven Sample alias, reads only relevant obs columns, and checks
80 identical sample universes /condition/tissue homogeneity /45kidney+23cortex+
12medulla /region NA mapping /noPapilla /suspension confounding. No expression
read, full H5AD transfer/rehash or new fit is authorized by this check.


### Bassez HR22 launch recorded — 2026-11-09 04:43:52Z

Prepare passed; reconcile classified ABSENT with all invariants passing and
no discrepancies. Launch returned RUNNING, exact command digest
`6fc8dc0e04179be474693dd018eb9e7e650168e83f9de36339907cd50bf3f0a7`.
Exactly one unbounded durable wait is armed (local exec4481). Initial source/
runtime verification precedes Slurm submission; the first immediate log read
was empty and no relaunch occurred. Subsequent scope check proves actual
requested22 SHA c6d219... and dispatched one-row CSV SHA6320d848... match
this contract; the exact physical matrix contains only the named22 tokens.
Current emitted preflight IDs are4663753 and4663754. Terminal completion,
remaining emitted IDs, output contracts, NAS checks and review are pending.
No main-wave deduplication or completed-result reuse is claimed yet.

Reviewed no-compute BL publication adapter is committed/pushed as `3794e3e`:
three new finite files, original20RDS+10timings preserved, fresh copied release,
original failed gate remains FAILED. Two independent full fixtures, strict
special-bit permission rejection, atomic create-only receipt and transfer/
mount/checksum/owner/review failure checks pass. Control has been staged,
MD5 `fe1ecb6f7d3bd0a8cf81567b8446353b`, SHA256
`5f70cbd8d4e7ce0ffce0c044d816433dea4098df0a9f68c35085603ba713d84e`;
unchanged root-review proof SHA256
`ded1c6ee10b554386208c60378660b5541fc0605af26cd77a92df03c1ba39567`.
Publication execution remains pending until the current ecoda-benchmark gate
has terminal inspection/reviewer approval. No dummy Slurm job, second mutex,
new scientific fit or failed-gate status edit is permitted for this adapter.


## Bassez HR22 terminal failure and independent BL publication — 2026-11-09

HR22 gate ended FAILED1 at04:50:47Z. Exactly one terminal inspection covered
all five emitted IDs4663753/4663754/4663757/4663759/4663760. **Every accounting
row is COMPLETED0:0**: input and Renv preflights, composition worker, watchdog,
and aggregate all succeeded. The controller then rejected
`results/Bassez_metadata.rds` as an ACTIVE input, before the exact-request RDS
CLI ran. This is an evidenced output-vs-input validation ordering defect:
current selected output owners intentionally stay ACTIVE until publication,
but the post-aggregate validator used the reusable-input guard. Generic
foreign/ACTIVE input guards must remain strict. A narrowly scoped owned
current-output validation branch and a real-owner/real-RDS-CLI regression are
assigned. Check existing23RDS and timing/record/owner evidence before deciding
no-fit recovery. Original gate/statuses/owners/manifests remain unchanged;
no reviewer approval, bridge, successful-science reuse or main dedup is claimed.
The conditional41-row main proposal therefore remains conditional, not a launch.

Failure operational review (NOT formal gate/scientific approval):
`.gate/bassez_author_hr22_inputs_20261109_v1/failed_terminal_review.json`, SHA256
`6d891b46f6c93ab0bdc797bc607c9d0f84cf43d4c3e1d1bbb6c0d4f8e16f01aa`.
It binds the failed terminal manifest and sole inspection, keeps all HR22
outputs unapproved, and permits only the independent already-reviewed BL10
publication run to proceed. That run consumes no HR22 data and has a distinct
output/sync root; there is no active ecoda-benchmark gate. The earlier local
operational precondition of a successful HR22 approval was unnecessarily
strong for this unrelated no-compute publication and is superseded by this
explicit terminal-inspected failure review. No failed gate is approved or
passed as a successful dependency. Stephenson annotation and16batch validation
continue independently.

### BL10 publication-only execution contract — READY

- Run `stage5_bl10_publication_20261109_v1`, no scheduler submissions, no H5AD
  reads/transfers, no model fit. Scope remains exactly20 reviewed RDS+
  10 original scientific timing Feathers. New copied bytes/sidecars and
  60-entry NAS checksum table; source producers and failed accounting remain
  explicit. Old BL run/gate/status/owners/checksums stay immutable.
- Source `3794e3e98b6d87bb18164a86747b0b624dc25c55`, snapshotparent
  `20261109-bl10-publication-v1`, archive SHA256
  `932528f9bee6354b190e5254936991eb0ca151e6982b16aed973c1198bedb40a`.
  Same pinned runtime image/manifest and frozen auxiliary identity.
- External control/review are staged and sealed in
  `_ecoda_inputs/20261108-confounding-subsets-v1/bl-science-publication-v1`.
  Control MD5 fe1ecb6f7d3bd0a8cf81567b8446353b, SHA5f70cbd8...;
  original root scientific review SHA ded1c6ee... . Exact full hashes above.
- Fresh run/log/NAS namespace and actual NAS mount are required. Common
  selected sync owner protects the new root. Root operational helper enforces
  the exact terminal-inspected HR22 failure review and unchanged BL controls;
  no new global group/lock/gate or dummy scheduler ID is invented.
- Completion requires the pinned helper's read-only `validate` mode over
  original evidence and all copied mappings/owners/records/sidecars/modes/
  global checksum plus terminal state, followed by explicit root review.
  Adapter writes a fresh publication receipt and truthful NOOP_VALIDATED
  no-fit terminal record. It never changes the original failed BL gate.

Exact immutable command (SHA256
`5c92d791b52f625925476c382a962816965ed8dc0a2e210de6e3bb6db74fb57a`):

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bl10-publication-v1/3794e3e98b6d87bb18164a86747b0b624dc25c55/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bl10-publication-v1/3794e3e98b6d87bb18164a86747b0b624dc25c55/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bl10-publication-v1/3794e3e98b6d87bb18164a86747b0b624dc25c55/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_bl10_publication_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_bl10_publication_20261109_v1 --script src/utils/bash/publish_bl_science.sh -- --control /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/bl-science-publication-v1/control.json --control-md5 fe1ecb6f7d3bd0a8cf81567b8446353b
```


## Reviewed input and scientific-result milestones — 2026-11-09

All16 Stage3 batch-view artifact reuse is approved independently of the old
FAILED global gate. Receipt MD5 `04e7131fe49d0904435c4bee0d6fcea5`, SHA256
`c8e1d9b5b0a4b485ecdb228c706c6f8a3de65f0861fc57b148de68cd39bb05b4`.
Independent review reapplied named masks, min500 and canonical-ID rules,
compared paired memberships/configuration hashes, and recomputed all52tables.
Only the two documented Ncl43Covid/10Healthy18.87% exceptions remain; COVID
counts22/12/10, structural missing severe group allowed. All32scratch/NAS
ownersOK, payload/proof hash/size agreement,16published scratch records and
selectedNAS777 chains are evidenced by the pinned audit. Old24accounting rows
and4OOM witnesses are preserved; no new scheduler queries/submissions.

Kidney obs-only report MD5 `fc3d0b74371feaca6dcc1d1cde115691`, SHA256
`db1e75e47d95eb2a275ab5d65042eae9d96fcbc11084e75f143048169f6df99a`.
Configured specimen IDs and the producer's canonical Sample alias agree;
condition/tissue/region/suspension are homogeneous per sample and match across
bothviews. NoPapilla.45kidney samples are cells;23cortex+12medulla are nuclei.
This confirms the stated anatomy/protocol confounding caveat, not another
filter or correction covariate. Tissue and condition share all feature and
distance outputs; evaluation/reporting alone supplies separate results.

Stephenson Stage4 firstinspect is complete: six roots4663688/4663694/4663744
and4663689/4663721/4663745, all32rowsCOMPLETED0. One accounting operation only.
Bridge preserved raw SHA7f804bda... and derived6-root projection SHA71b3330a...;
no rows were invented or failed rows omitted. Pinned fc2949f
`ecoda_run_audit.sh` passed, including source/runtime/frozenaux, parent/copy
protection, annotation completeness, payload/schema/ownership and actual NAS
hashes/equality/modes. Formal root approval is
`.gate/ecoda_stage4_stephenson_annotation_20261109_v1.review.json`.
The derivative publication helper already normalizes/validates its full
selected NAS ancestor chain; a suspected issue in a separate legacy wrapper
branch was not applicable. No unnecessary chmod/repair/source change occurred.
New annotated input MD5 `7d92fb6e0d3d7dd5c3500212fd951493`,size13,357,978,168.

BL10 publication completed05:03:40Z and pinned terminal validation passed.
New run `stage5_bl10_publication_20261109_v1`: receipt MD5
`3fc1f22c595acea6d071c50aebcd0799`, SHA256
`61ad382358ad576736e8c66336d88aa612d64f3caface951f4f50a7540a836db`.
Review `.gate/bl_publication_inputs_20261109_v1/publication_review.json`, SHA
`c2901b0f57ef8753b3a69919ae612ca5d36faafa84abda46041d5c61c8efa401`.
Exact30payloads/60checksum entries, strictNAS666/777, source/runtime/oldFAILED
accounting preserved; no fits or H5AD reads repeated. New NAS root is
`releases/20261108-confounding-subsets-v1/stage5_bl10_publication_20261109_v1/benchmark`.

HR22 actual read-only science audit passed:23RDS+one22-row timing Feather,
26ordered samples/74948cells/40HRtypes, canonical5finite scores, groundtruth
and storedNULL semantics, baselinecounts and unchanged features/distances,
all producer records/owners/full hashes/input/source/runtime/aux proofs.
Configured R4.5.2; no H5AD consumption/rehash, rescoring or method rerun.
External R assertion SHA `c3de2329c9a3e52a515bc19d667a622ffcdc520dc264c2fcd4329dbabef8ba70`,
exact23-row CLI manifest SHA026c8103486dd05b60120421f4522d0e01df9d57595a3dfbd933e655abbdc398.
Science receipt MD5 `3ccc48e3e352286e500d862da67c9a67`, SHA256
`78523c4ad31a9768085c0be52a6576dd7cb97a511d1fdc54b69619a25ede74c4`;
root scientific review SHA
`9828da9b4c4046275b4774f1f4e93777c994b0928e6be8c724e71bb1ee6b32b1`.
This approves scientific reuse; old FAILED gate/statuses stay unchanged.
Future pending-output lifecycle fix039d79c is independently accepted and
source-frozen for the main wave; generic input guards remain unchanged.

## HiTME clarification from the user — NO additional targeted refresh

The user requested clarification and explicitly withheld any targeted refresh
of100parameter results. None was selected, implemented as a selector, or run.
The prior proposed100-row expansion is withdrawn from the current scope.
Do not present it as a prerequisite or infer approval from elapsed time.

What the audit establishes: plain `ECODA_HiTME_HR_layer2` matches current L2
counts exactly in Bassez,Zhang,Pelka; plain L3 also matches its own metadata.
Separate `ECODA_HiTME_HR_layer2_top_varexp0..0.9` entries in three explicitly
selected checksum-table-valid aggregate RDS contain L3-style state categories
and fail comparison to currentL2. The current deferred variance builder had
a reproducible loop-capture bug; fixe776105 affects only those future variants,
not plain baselines. Historical producer ownership/source provenance for those
old aggregates is absent, so this does not establish whether the mismatched
entries came from intentional L3 tests or that historical bug, nor which RDS
bytes generated a past published figure.

Maintained code/retainedHTML traces support the user's intended distinction:
Fig2A omits HiTME; Supp15 extended/presentation selects exact plain L2
(`notebooks/benchmark_analysis.rmd:2908,3005`); annotation/PCA/MDS comparisons
also select plainL2. Supp2's active filter excludes HiTME and variance entries;
its earlier method list is unused. Fig4/Supp20 HVC analysis uses authorsHR.
The L2variance0.5 reference at2429 is diagnostic feature-name printing.
`src/utils/constants.R` extended normalization pool contains only plainL2;
exact membership filtering prevents excluded variants from affecting scaling.
There is no evidence here of L3 accidentally replacing the final plain L2
benchmark. Further targeted work requires a new user decision after this
clarification. Preserve intentional L3 tests and old artifacts.

## Bassez HR22 no-fit publication contract — READY

New run `stage5_bassez_hr22_publication_20261109_v1`; fixed24payloads=
23reviewedRDS+one original selected timingFeather,48NASchecksum entries.
No scheduler, new fit, rescoring, H5AD read or oldrun write. Reuse the accepted
publisher's copy/owner/NAS/mode/atomicreceipt machinery via the second explicit
fixed scope `ecoda_bassez_hr22_publication_control_v1`; BL10 behavior unchanged.
New copy producer `stage5_bassez_hr22_release_copy`, receipt
`ecoda_bassez_hr22_publication_receipt_v1` /`bassez_hr22_publication.json`.
Author and independent HR fixture plus unchangedBL regression pass.

Source `fc46e0366226a4dcce2c2cb0366678df167d3aab`, separate snapshotparent
`20261109-bassez-hr22-publication-v1`, archive SHA256
`e2e00e96e3bc2293567ff889d5e639919fd149aa602f9fc8716c7488eb2d52e8`.
Same frozen runtime/aux. Control MD5 `214bd3b9054026f4dccc523065fbcf28`, SHA256
`46bfbe58ca8edae16160fc03f7c2e7aec67454aa70629463df7046031516c69c`;
review SHA9828da9b... binds the science78523c4... receipt. Stage both under
`_ecoda_inputs/20261108-confounding-subsets-v1/bassez-hr22-publication-v1`.
Require fresh owned run/log/NAS roots and terminal customvalidate+rootreview.
Historical stdout validation marker is shared withBL; HR schema/report/copy
producer and24mappings determine scope, not that text. Original gate remains
FAILED even after successful new publication.

Exact immutable command (SHA256
`752af53115482b45e9bdd4e4913e74ce8d152db1282aec15c73bb64b38e884de`):

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bassez-hr22-publication-v1/fc46e0366226a4dcce2c2cb0366678df167d3aab/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bassez-hr22-publication-v1/fc46e0366226a4dcce2c2cb0366678df167d3aab/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-bassez-hr22-publication-v1/fc46e0366226a4dcce2c2cb0366678df167d3aab/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_bassez_hr22_publication_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_bassez_hr22_publication_20261109_v1 --script src/utils/bash/publish_bl_science.sh -- --control /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/bassez-hr22-publication-v1/control.json --control-md5 214bd3b9054026f4dccc523065fbcf28
```


## Main parallel method wave — READY for exact-once launch

- Gate `ecoda_stage5_main_confounds_20261109_v1`; parent run
  `stage5_main_confounds_20261109_v1`; existing group `ecoda-benchmark`.
  One existing campaign launches three isolated child dispatchers. Failure of
  one child does not cancel its peers; only documented intra-dataset producer
  dependencies serialize methods. No new framework/group/global sync owner.
- Final source `c534cdfd28fc9a77c009f4e25fa4c4ddb7ee4400`, snapshotparent
  `20261109-main-confounds-v2`; archive SHA256
  `3e70d93de6ae914c5d377bf80502e05befca02d158fb2673e5588ef446e09044`.
  Earlier source039/main-v1 snapshot and proposed wrapper were never executed
  and remain immutable history. Final metadata-mode path fix resolves bound
  scratch/NAS identities and input, preserving all legacy branches/ownership/
  checksum/column guards. Real-worker/Feather tests pass for both views,
  hidden NAS, legacy positives, Kidney2axes and foreign/tampered negatives.
  Author+independent review passed. No scope/reference/task-count change.
- Runtime image8fcd00... and manifest1c4c2c... remain pinned, R4.5.2/Python3.13,
  no environment refresh or dependency/aux change. Aux MD58cacee... unchanged;
  no new sex-gene exclusion/chromosome annotation. Dataset configSHA522342...
  remains the confirmed top-level column/view authority.
- Sole formal predecessors: approved Bassez Stage2 and approved Stephenson
  Stage4 gate manifests. Never pass failed Stage3/BL/earlyHR22 gates as
  successful dependencies. All16 batch artifacts have explicit NOOP reuse
  review5c05ac54...; current Kidney metadata separately passed.
- Bassez HR22 publication completed05:48:23Z,24payloads/48checksum entries.
  Receipt MD5 `1be5722bf1a3e839da7f5e28e9d15635`, SHA256
  `dcbd31bb47403c4ab1181e808d0dbc0d12bbd250617c070bc615ed13f6d7f4df`.
  Root publication review SHA256
  `3d2227f027706bf7738b9d3aba5fda09fa4aa9fe88f6e41ccbbd9c4ad5457b62`.
  The pinned HR terminal validator passed. Original failed gate remainsfailed;
  original23RDS/timing/owners unchanged. Thus the exact22 valid author-HR
  composition rows stay outside recomputation. Source/science approval—not
  oldgate success—is the reuse authority. Other20Bassez requests are retained.
- Benchmark selection41 public rows /29 normalized families: Stephenson12
  whole-family/screen requests; Bassez20 affected HR/screen requests; nine
  additional transformation-only requests. No BassezBL duplication; StephBL
  occurs once within its whole-method row. No extra HiTME100 selector or fit.
- Both batch contexts select exactly8cohorts/56families each. Inputs are27
  actual six-field references:11benchmark,8uncorrected,8corrected. All
  path/size/strict3fieldproof references checked onBamboo without new H5AD
  reads; scheduled preflights validate payloads. External controls13files
  are hashed/sealed under `_ecoda_inputs/20261108-confounding-subsets-v1/main-confounds-v1`.
  Existing039 proposal files remain evidence only; `wrapper.actual.sh` is
  SHA09c8c482b659ea207ebc613cc24fc01a2f90a4ed0d9dbbe0501056b33527b551.
- Exact selection/input hashes: benchmark selection56119d97..., normalized
  038658d1..., inputs3fad7dad...; uncorrected selection1da3b269...,
  inputs7fe1eb3f...; corrected selection4bbc10c8..., inputsb5f506b5...;
  contextsdd5d41be... . Full hashes remain in `.gate/ecoda_stage5_main_confounds_20261109_v1.spec.json`.
- Physical scope: benchmark54method/preptasks in19arrays+19watchdogs+1aggregate;
  eachbatch56tasks in7arrays+7watchdogs+1aggregate. Production preflights:
  H5AD27tasks/3IDs, Renv3/3, metadata exports16/2. **77 initial schedulerIDs /
  248 scheduled task units**, before explicit bounded OOM retries. Total166
  method/prep tasks,33watchdogs,3aggregate. Expected scientific payloads223
  (benchmark63RDS/Feather,uncorrected64,corrected96), plus16sample-metadata
  exports and ordinary sidecars/timing/log/control records.
- Resources:128G initial→500G OOM ceiling; methods12h; CPU16 except trans/
  zeroimp4; GPU8CPU+1GPU,16benchmark GPUtasks(MRVI2k1+scPoli15), batchCPUonly.
  GPU class throttle4, CPU literal1000 bounded by selected rows. Watchdog12h.
  Explicit source-bound wrappers retain OOM-only row retries and fail closed
  for non-OOM/exhausted limits. If physical rows exceed declared selections,
  cancel every emitted schedulerID and runner, preserve evidence and fail.
- Parent/threechild run/log/NAS/gate namespaces verified fresh before prepare.
  Actual NAS preflight passes on `smednas515:/volume1/Shared`; no synthetic
  scratch fallback. Three independent roots/sync owners are:
  `stage5_main_benchmark_20261109_v1/artifacts/benchmark`,
  `stage5_main_uncorrected_20261109_v1/artifacts/batch_effect/uncorrected_final`,
  `stage5_main_corrected_20261109_v1/artifacts/batch_effect/corrected_final`.
  NAS mirrors are under `releases/20261108-confounding-subsets-v1/<child>/...`.
  Child logs are within parentrun/logs/<child>, as the campaign actually sets.
- After prepare/reconcile/launch arm exactly one unbounded durable wait.
  Terminal review uses one union accounting query for every emitted child ID,
  the accepted campaign raw-capture/three-child projection adapter and three
  pinned selected run audits. No ordinary parent-only audit or repeat query.
  Root approval/release follows all selected science/artifact/NAS checks;
  failed child rows are recovered narrowly without recomputing valid rows.

Exact immutable command (SHA256
`950bb25c9e224f1642e2b63929e59cdfb8145c7eda745ddd7103475e9ee3f518`):

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper BENCHMARK_MEM=128G BENCHMARK_MEM_MAX=500G BENCHMARK_CPU_CPUS_PER_TASK=16 BENCHMARK_GPU_CPUS_PER_TASK=8 BENCHMARK_CPU_TIME_LIMIT=12:00:00 BENCHMARK_GPU_DEFAULT_TIME_LIMIT=12:00:00 BENCHMARK_GPU_ANY_TIME_LIMIT=12:00:00 WATCHDOG_TIME_LIMIT=12:00:00 MAX_NUM_CHUNKS_PARALLEL=1000 BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE=4 BENCHMARK_GPU_ANY_ARRAY_THROTTLE=4 /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-main-confounds-v2/c534cdfd28fc9a77c009f4e25fa4c4ddb7ee4400/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-main-confounds-v2/c534cdfd28fc9a77c009f4e25fa4c4ddb7ee4400/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-main-confounds-v2/c534cdfd28fc9a77c009f4e25fa4c4ddb7ee4400/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_main_confounds_20261109_v1 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_main_confounds_20261109_v1 --script src/5_run_benchmark_methods/stage5_campaign.sh -- --contexts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/main-confounds-v1/contexts.tsv --derivative-release 20261108-confounding-subsets-v1
```


## Main wave setup failure — terminal evidence, 2026-11-09

Gate `ecoda_stage5_main_confounds_20261109_v1` launched06:17:39Z from
`c534cdfd28fc9a77c009f4e25fa4c4ddb7ee4400`, command SHA256
`950bb25c9e224f1642e2b63929e59cdfb8145c7eda745ddd7103475e9ee3f518`,
and endedFAILED1 at06:46:38Z. The one unbounded wait has finished; one
terminal inspect already captured all actual scheduler IDs. Preserve this
failed gate, source snapshot, owners, logs, manifests and accounting unchanged.

Actual IDs:4663922 benchmark input preflight(11tasks),4665386 R environment,
4665398 Stephenson pseudobulk preparation,4665399 its watchdog,
4663913 uncorrected metadata(8tasks),4665389 corrected metadata(8tasks).
All30raw accounting rows completed0. Bulk method arrays never launched.

Confirmed benchmark failure: parameter-matrix ownership projection repeated
the same GloScope family, triggering the existing duplicate-artifact guard.
Future repair deduplicates family ownership while preserving distinct physical
parameter rows and the duplicate-path protection. A separate scheduler recorder
checksum lifecycle bug prevented failed children from exposing their IDs to
the parent campaign; root supplied the corroborated6IDs to terminal inspect.

All16 metadata worker logs report successful exports with expected sample
counts; both controllers nevertheless rejected final acceptance. The omitted
explicit analysis-variant flag is NOT a proven cause: the exporter deliberately
uses the exported environment fallback and real16-table fixtures pass both
forms. Diagnose the original unsuppressed host check/status/proof before
patching or rerunning any export. No production metadata fix is yet approved.

Stephenson preparation may already contain all6pseudobulk variants, including
the corrected blacklist cache. Inspect exact caches/proofs/input identity and
timing without compute before deciding reuse. Its scored blacklist result is
still missing; do not count preparation alone as completed pseudobulk_bl.

Recovery must use a fresh run identity and reviewed committed source, adopt
only individually validated prepared artifacts, exclude completed BL10 and
Bassez authorHR22, and record the revised exact selection before launch.
Current work is controller diagnosis/review; no HPC method arrays are running.


## Controller recovery — reviewed fixes and prepared-output audit

User explicitly directed continued controller repair and parallel launch with
completed-output reuse. Independent subagents implement/review bounded source
changes; the orchestrator owns Git, plan, controls and HPC execution.

- `cba58c3702c82c73217f82353048924a7c94fc76` committed/pushed: unique family
  ownership projection preserves every physical parameter row; scheduler IDs
  receive a checksum at every append and existing evidence must validate first.
  Real ownership and watchdog fixtures, campaign failure/accounting fixture,
  syntax and independent cross-review pass.
- `80e1d88` committed/pushed: metadata acceptance now reports the precise
  scheduler/checksum/schema/status/report failure boundary. Both8-row context
  tests and11negative boundaries pass; independent review accepted. This is
  diagnostics, not a claim to have reproduced the historical batch failure.
- Actual original16metadata host checks all pass with exact original outputs,
  status filenames/TASK_ID1..8, MD5/size/PATH. Evidence is
  `.gate/main_confounds_failure_20261109_v1/metadata_host_checks.json`.
- Original Stephenson6caches pass the actual c534 strict RDS CLI with53ordered
  samples, expected feature dimensions, current annotated input identity and
  unchanged content/proofs/records. BL uses unchanged auxMD5
  `8cacee7893cbb9d170cc9de87e3bd306`, five approved groups, before-fit policy,
  and independent370.6966s variant timing. No sex-gene exclusion was added.
  Original7-row preparation timing is retained; no fit or H5AD expression read
  occurred during this audit. Cache audit SHA256
  `5ece0f9a7a5f70bad1509ac731261ce192aed1b7a6f3c55c6bf538d7c691b30c`.
- All6cache owners were finalizedFAIL by the later controller abort. Preserve
  that truth. Their preparation watchdog reportsOK and the original sole
  terminal inspect includes array4665398/watchdog4665399 completed0. Fresh
  copy-only acceptance is based on independent artifact evidence, not an old
  owner/gate promotion. No old status, owner, checksum or manifest was edited.
- Approved exact23-artifact copy-only receipt SHA256
  `be7cfc7ab3dda369c8167c4fdef043f3501d4f97b87e823f63f194f73cfe1800`:
  six Stephenson caches +one original timing Feather +16batch sample metadata
  Feathers. Local controls/review: `.gate/main_recovery_inputs_20261109_v2/`;
  remote control: `_ecoda_inputs/20261108-confounding-subsets-v1/main-confounds-v2/prepared_adoption.json`.
  New target children are `stage5_main_{benchmark,uncorrected,corrected}_20261109_v2`.
- The finite adoption hook is DONE/reviewed in9ca2453, using existing derivative bindings,
  strict proofs and create-only new copies. It will not relax the legacy
  validated-recovery contract or submit work. Existing scientific validators
  must accept the new copies before completed preparation is excluded.
- Remaining41public benchmark requests and both8-dataset batch selections are
  unchanged. No BL10 or Bassez authorHR22 recomputation, no targetedHiTME100.
  Fresh recovery gate is PLANNED, not launched. Final immutable source,
  selected pending row count and exact wrapper will be recorded before launch.


## Recovery parallel method wave v2 — RUNNING since2026-11-09 08:18:46Z

- Gate `ecoda_stage5_main_confounds_20261109_v2`; parent
  `stage5_main_confounds_20261109_v2`; existing group `ecoda-benchmark`.
  Children `stage5_main_benchmark_20261109_v2`,
  `stage5_main_uncorrected_20261109_v2`, `stage5_main_corrected_20261109_v2`.
  The three controllers start independently/concurrently and a peer failure
  does not cancel the other contexts. Only documented within-context
  preparation dependencies serialize tasks.
- Exact source `9ca2453eea7320a721ef6dfd458f9db108abecac`; snapshot
  `_ecoda_source_snapshots/20261109-main-confounds-recovery-v1/9ca2453eea7320a721ef6dfd458f9db108abecac`;
  source archive SHA256
  `39b4d7bedad8806e69d44dcd1c7ef545cf08193d63b0ecac0e36ae4d48af77d0`.
  Existing versioned format2 runtime/aux/config/lock constraints unchanged.
- Authoritative final spec is `.gate/main_recovery_inputs_20261109_v2/spec.ready.json`,
  SHA256 `b687bde28eff9eefda175eb8d61a55620e3460efec583d4354ccc322f1dfd9d0`.
  Wrapper SHA256 `399cbff4e968595749ca4549e430af781d8c621c2feb0fd667f365088cae6bca`.
  The unlaunched `spec.json` proposal is preserved and superseded: review caught
  stale descriptive v1 run/control paths; finalready fields and executable
  command all agree onv2. No gate was prepared/launched from that proposal.
- Public selection remains41benchmark requests→29selected families; validated
  preparation reuse leaves28pending families/53physical benchmark rows, plus
  56uncorrected+56corrected rows =165pending method/preparation tasks.
  Expected initial scheduler scope:32method arrays,32watchdogs,3aggregate gates,
  3H5AD-preflight arrays(27tasks),3R environment jobs =73IDs/230task units.
  Bounded OOM-only retries may add IDs for the same selected rows. No metadata
  export arrays and no Stephenson preparation array/watchdog are required.
- This arithmetic derives from the original77IDs/248units minus one completed
  preparation task/watchdog and two8-task metadata arrays. The full33-matrix
  scope/real-owner audit passed; separate actual-cache dispatcher/dependency
  regression proves that completed preparation is skipped. Pending scientific
  payload scope remains223. No completed BL10/HR22 or HiTME100 expansion.
- Independent review accepted real six-cache/RDS/BL-policy/cache-reader reuse,
  both8-cohort realFeather/exporter/coordinator reuse (`PENDING=0`/nojobs),
  create-only/proof/schema/foreign-owner negatives, and campaign sealed-receipt
  forwarding/concurrency/failure isolation. No fits/count reads in these tests.
- Remote source archive/control hashes and all27input references passed
  preflight; newrun roots absent. NAS `df` proves
  `smednas515:/volume1/Shared` mounted at the configured Bamboo target with90TB
  available. Formal dependencies remain only the reviewed Bassez Stage2 and
  Stephenson Stage4 gates; no failed gate is promoted or passed as a dependency.
- Resources unchanged: CPU16(method screens4), GPU8CPU/oneGPU, GPU throttle4,
  128G initial/500G bounded ceiling,12h method/watchdog limits. Per-row scheduler
  manifests and checksum proofs are persisted at each submission boundary.
- Prepare/reconcile/launch exactly once, then arm one unbounded waiter. Record
  launch outcome immediately. Terminal inspection/review remains required;
  launch does not claim scientific results or figures are complete.

Exact immutable wrapper command:

```bash
env ECODA_HPC_CLUSTER=bamboo ECODA_NAS_TARGET_DIR=/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper BENCHMARK_MEM=128G BENCHMARK_MEM_MAX=500G BENCHMARK_CPU_CPUS_PER_TASK=16 BENCHMARK_GPU_CPUS_PER_TASK=8 BENCHMARK_CPU_TIME_LIMIT=12:00:00 BENCHMARK_GPU_DEFAULT_TIME_LIMIT=12:00:00 BENCHMARK_GPU_ANY_TIME_LIMIT=12:00:00 WATCHDOG_TIME_LIMIT=12:00:00 MAX_NUM_CHUNKS_PARALLEL=1000 BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE=4 BENCHMARK_GPU_ANY_ARRAY_THROTTLE=4 /bin/bash /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-main-confounds-recovery-v1/9ca2453eea7320a721ef6dfd458f9db108abecac/tree/src/utils/bash/ecoda_source_snapshot.sh exec --source-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-main-confounds-recovery-v1/9ca2453eea7320a721ef6dfd458f9db108abecac/tree --source-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/20261109-main-confounds-recovery-v1/9ca2453eea7320a721ef6dfd458f9db108abecac/identity/source.manifest --host-env-prefix /home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13 --runtime-image /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif --runtime-manifest /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-6bbf70b-relocated/ecoda-py-cuda13.sif.manifest --run-id stage5_main_confounds_20261109_v2 --scratch-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper --logs-root /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_run_logs/stage5_main_confounds_20261109_v2 --script src/5_run_benchmark_methods/stage5_campaign.sh -- --contexts /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/main-confounds-v2/contexts.tsv --derivative-release 20261108-confounding-subsets-v1 --prepared-adoption /srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_inputs/20261108-confounding-subsets-v1/main-confounds-v2/prepared_adoption.json --prepared-adoption-sha256 be7cfc7ab3dda369c8167c4fdef043f3501d4f97b87e823f63f194f73cfe1800
```


### Recovery exact-once launch — 2026-11-09 08:18:46Z /09:18CET

Gate `ecoda_stage5_main_confounds_20261109_v2` prepare succeeded; reconcile
classifiedABSENT with no discrepancies and every profile invariant passing.
Launch returnedRUNNING and accepted only the two approved Stage2/Stage4
dependencies. Command SHA256
`ad0a8c35f98860c8e3626e843f2023fa6ef9d876b941ea325e5717666c4463ee`.
Exactly one unbounded durable waiter is armed (local exec36153); no scheduler
polling and no duplicate runner/waiter. Immutable scientific source remains
9ca2453; later plan-only commits do not change the executing snapshot.

Startup may still be in copy-validation/preflight/submission. RUNNING refers
to the durable campaign, not a claim that every method array has already
started or completed. Inspect actual run-owned adoption/classification/ID
records for scope evidence. The formerv1 gate remainsFAILED unchanged.
Terminal accounting and pinned artifact/source/runtime/NAS review remain
required before accepting scientific results or releasing downstream gates.


### User-requested one-off startup check — 2026-11-09,09:36–09:40CET

The user requested a current check, any necessary failure repair with reuse,
then a short status and a pause pending their response. No periodic30-minute
monitor or additional waiter was created; the existing durable waiter remains
the sole completion observer. The server-side tmux runner is independent of
the workstation/SSH connection and can continue while the workstation sleeps.

- Durable status remainsRUNNING/nonterminal. All3named child controllers are
  alive in the same server runner's process tree; no child terminalFAIL exists.
- Each child scheduler_ids.tsv is still empty and no pending_selection.tsv
  has been finalized: **zero method jobs submitted/running/queued at this
  check**. CampaignRUNNING must not be described as method jobs running.
- Benchmark has copied6Stephenson caches+originaltiming; the actual strict
  RDS artifact-list validator reportsOK. Uncorrected copied8sample metadata
  Feathers and its metadata report isSTATEOK,COUNT8,PENDING0: no export rerun.
  Corrected is still before adoption, validating its input contracts.
- One read-only process-tree snapshot (no Slurm polling) identifies active
  content hashing through derivative_artifacts.py --verify-content for
  benchmark/uncorrected inputs, including the large Alzheimer input, and
  benchmark_h5ad_contract.py for corrected Alzheimer. This explains the
  current startup delay; no demonstrated controller/worker failure warrants
  interrupting the run or bypassing validation. No fit or valid output was
  recomputed during this check, and no source/runtime/owner was changed.
- Existing BL10/HR22 and all original prepared artifacts remain preserved.
  Further intervention/monitor setup is deferred pending the user's reply;
  this does not pause or cancel the already launched Bamboo campaign.


### User resumed work; submission progress and P6 preparation — 2026-11-09

The user will leave the workstation awake and explicitly directed continued
work. This supersedes the earlier request to pause after a one-off check.
Keep the existing server runner/single waiter; do not create another launch,
waiter or periodic scheduler-accounting loop.

At08:41:45Z benchmark compute-preflight4665493 had11/11workerOK statuses;
uncorrected4665505 had7/8OK (Alzheimer pending). At08:48:48Z benchmark had
also recorded R-environment preflight4665514 and passed all selected source
contract checks. Both batch metadata reports now exist and areOK,COUNT8,
PENDING0. Therefore all23approved prepared payloads have been copied without
refits/export jobs; current full input validation continues. No worker or
controller failure was observed and no running source was changed.

Important evidence correction: an empty scheduler ledger alone cannot prove
that no scheduler work exists. The shared H5AD helper records its ID only
after sbatch --wait returns; worker filenames/statuses can show in-flight
preflight work sooner. All earlier zero-method statements concern method
arrays, not these compute-preflight arrays. Future source-only improvement
may persist IDs before waiting; the current source is immutable.

Read-only startup/performance audit explains the delay: corrected source
contract checks stream counts on the login controller, bound sidecar checks
hash whole selected H5ADs sequentially, compute preflight repeats mandatory
schema/content checks concurrently, and source artifact checks can repeat
count scans. Source-identity helper itself reads obs/proofs rather than doing
another full hash; immutable runtime validation uses the existing fast path.
Do not interrupt healthy checks or bypass integrity to accelerate this run.
Future-only candidates: timestamped per-input phase logging; persist an
in-flight ID before waiting; reuse exact current-run preflight evidence only
for proven immutable inputs with a focused regression. No such change is
implemented or hot-swapped into the running snapshot.

P6 implementors resumed the previously frozen compact metadata projection,
streaming local intake and exact supplementary-scope audit. No production
figure execution/intake or RNA-bearing local transfers are authorized by a
code fixture passing. The batch projection-required guard is independently
assigned and tests both Kidney axes/360+405values. Publication parameter
scope beyond the231core benchmark rows is being enumerated explicitly;
no targetedHiTME100 or additional model rerun is inferred.

Finite checksum/provenance audit completed outside scientific outputs:
`_ecoda_inputs/20261108-confounding-subsets-v1/p6-legacy-evidence-v1/legacy_inventory.observed.json`,
SHA256 `f7e732aa07a28acce44a6764b64277ecfcf9fae5dfb954a92f4c3ab7c3115e4a`;
local metadata copy `.gate/p6_legacy_evidence_20261109_v1/`. Exact114result
paths plus one timing Feather were stream-hashed, with no deserialization,
H5AD/count/SIF read, scoring, fit, scheduler query or old-file write.
106payloads match the frozen canonical checksum table; nine existing
zero-imputation companion RDS files disagree (Adams,Gongsharma,Kfoury,Kim,
Lee,Pelka,Smillie,Wu,Zhang). This is **not** scientific reuse approval and
does not authorize nine reruns or checksum-table repair.

All115legacy global owners are absent; adjacent receipts retain historical
home-alias PATH values. Missing historical ownership/provenance must remain
truthful and cannot be invented/promoted; it is not by itself evidence that
the stored scientific distances are invalid. Prepare the bounded science
audit against externally frozen exact MD5/size/current-input/key/label scope,
and investigate the nine checksum disagreements using finite independent
backup/receipt evidence before deserializing or selecting any recompute.
Actual full science/projection compute remains gated after the active wave's
terminal review in the existing ecoda-benchmark serialization group.


### Actual method launch and isolated GloScope failure — 2026-11-09

At09:21:02Z benchmark had28pending family rows and4method arrays/4watchdogs;
uncorrected had56pending rows and4method arrays/3watchdogs (another submission
in progress); corrected had8/8compute preflightOK under4665528 and was finishing
source/consumer checks. Benchmark preflight IDs4665493/4665514 and uncorrected
4665505/4665515 are recorded. No source/runtime/config/snapshot was changed.

MOFA array4665525 published its7-combo Stephenson result and its watchdog4665526
reportsOK. GloScope4665516 siblings progressed/completed, but task2 failed
before R with ENOENT for the run-owned manifests/input_artifacts.tsv. Exact
actual matrix row2 is Stephenson/benchmark_analysis/gloscope/hvg2000_pcadims30;
matrix MD5 `9732a9b209be8401c14679c321773e45`,size285. Local read-only forensics:
`.gate/main_recovery_failure_20261109_v2/worker_forensics.json` and
`failed_row_matrix_and_stat.txt`.

The host TSV exists with exact binding SHA3fad7dad..., size3939, inode
4761855699885302067, and unchanged creation-era mtime/ctime. The same failed
worker passed earlier post-container binding validation and emitted its
ECODA_WORKER_DISPATCH banner, then failed during ecoda_derivative_input
directory/content resolution. Other parameter tokens/MOFA use the same mounts.
No deletion or parameter-specific mount path was found; a static mount defect
is not established. An exact tiny resolver fixture reproduces the phase when
the TSV is deliberately removed after successful validation and passes when
restored, proving the error boundary only, not the production cause. The one
retained scontrol record lookup found the old task already purged; no terminal
accounting query was performed. Do not infer missing data or repair sidecars:
the TSV is bound by SHA in artifact_bindings.json and has no separate.md5 by
design.

Independent review accepted failure-only shared-resolver diagnostics plus its
focused and existing derivative regressions: identical successful stdout/args,
original failure code retained, selected context/path stats on stderr only,
no retries, source/proof/owner/permission changes or cause-fix claim. This
applies to a future immutable snapshot only; active9ca workers stay unchanged.

Concrete one-row v3 recovery proposal is prepared but NOT launched:
`/private/tmp/ecoda-single-gloscope-recovery-proposal-v3/proposal.json`, SHA256
`fd618fd792bc03af8fd9bad362b4d7cb1a2c17a5bfe20ecebed2bc6fcc3a11f2`.
Existing requested-scope sharding/actual standalone RDS validator tests pass:
one CPU16/128G/12h gloscope row, one watchdog and one aggregate, plus one
H5AD+oneR environment preflight =5initialIDs/5units. No preparation, other four
GloScope parameters or MOFA fit is selected. Formal dependency would be only
the approved Stage4 manifest, never the failed/current parent. Wait for active
v2 terminal accounting/inspection and the complete failed-row catalog before
sealing a new recovery source/control; retain every valid output. Other
contexts/arrays continue and have not been cancelled.

P6 bounded observations completed meanwhile:
- All9NAS zeroimp counterparts match the frozen canonical table and differ
  from scratch versions; no NAS adjacent sidecars exist. Identity metadata
  SHA76e360294240c9a3593a447e01bc43f1a09c55fb001490050f6d6de1c0550834.
  Nine tiny checksum-backed NAS summaries contain41–43settings each and all
  eight stored scores, including the canonical five, all finite. No distances,
  labels or scoring-policy witness are stored. Schema metadata SHA
  `e1e77131b91e2e3da1f5129eddb4047512cc425a75c9c34832c07a479587f696`.
  Scientific/current-input/LISI-definition compatibility still needs review;
  no nine reruns, old-sidecar/table edits or scientific approval are inferred.
- Eleven current H5ADs were probed for group/column names/shapes only: all22
  Leiden50/Harmony2 source columns and corresponding embeddings exist. No
  dataset values/expression groups/full hashes were read. Observation SHA
  `74c20960df40f10a3625edff15c2bdabab60044b50cb399e802eb006f2d11c77`.
- Numbered-figure scope proposal adds673explicit rows to231benchmark core
  results (904total;1089intake rows including existing metadata/batch/screens).
  Candidate sources:14Bassez publication,69main-wave,568legacy,22derivedTODO.
  No HiTME variance grid/PILOT-GM noncanonical expansion or optional9QC rows.
  Twenty derived rows are reuse/compatibility candidates, not new fits;
  Stephenson's two need current53-compatible metadata results. Old fast-obs
  CSVs exist locally as compact sample×cluster metadata; oldSteph55 is not
  valid for current53 without an exact current count/label comparison. Current
  category-array lengths866/34 versus old active counts881/31 show a schema
  difference, not by themselves proof of active assignment differences. No
  new RNA/PCA/Harmony/Leiden fit is needed to tally saved current labels.
- P6 projection/mirroring/streaming code and focused tests are green but FROZEN
  and not independently released/committed. Exact scientific origins and NULL
  permutations remain preserved, local RNA guards precede reads, and a future
  oneCPU64G2h pinned projection worker is prepared. Urgent failure diagnosis
  takes priority; no production projections/figures or extra gate has run.


### Compact analysis code accepted; all contexts submitted — 2026-11-09

At09:32:56Z benchmark had all18planned remaining method arrays/18watchdogs,
uncorrected all7/7, and corrected had begun its7-group method submission.
At09:34:10Z a SINGLE user-requested live queue snapshot, restricted to the
29then-recorded method-array roots, showed5benchmark tasksRUNNING/28PENDING,
5uncorrectedRUNNING/42PENDING, corrected32PENDING. This is a time-stamped live
snapshot, not a periodic polling/accounting loop or a claim that every task
is running. Metadata evidence is
`.gate/main_recovery_failure_20261109_v2/live_queue_snapshot.json`.
The only observed failure remains GloScope4665516_2; its watchdog reports
non-OOMFAIL, while MOFA7, BassezHR pseudobulk parameters, and the new zeroimp
group watchdogs reportOK. Keep successful work and let all independent jobs
finish before final accounting/targeted recovery.

Independent review accepted the P6 code package (no production run yet):
`analysis_release.R`, benchmark helper, batch helper, projection entry, and
focused projection/ledger/benchmark/supplement/cell-depth/batch tests.
Local intake rejects unprojected RNA-bearing RDS before reading it; thin
projections keep exact distances/labels/NULL/drop IDs/scalar scores/timing and
only the finite CT-feature whitelist. Per-cohort streaming avoids retaining
whole RNA releases. Explicit copied-root→local-root mirroring changes only
new access metadata/receipts/control hashes; raw copies, witnesses and parent
scientific origins are preserved. A fresh projection cannot launder legacy
data into required-new rows.

The projection entry uses the established oneCPU64G2h controller/worker and
pinned immutable runtime pattern; it records a returned scheduler ID even on
worker failure. Read-only --validate installs no mutating failure trap and
requires current-run stage5/OK ownership of root, inventory and each payload.
Negative validation tests preserve old evidence/owners byte-for-byte. All
focused independent tests passed; protectedRplots.pdf hash is unchanged.
No projection/evidence compute gate is launched during the active same-group
wave; no RNA data was transferred locally and no figure is claimed refreshed.

The remaining maintained active unnumbered true-feature PCA panels need a
separate bounded metadata-coordinate extension:11PB_hvg2000 cases plus4Pelka
latent/CT cases, with genuine source-feature PC1/2, explained variances and
selected5(top7Bassez)loading arrows. Full-space scores stay separate; do not
substitute distance-profile PCA for feature PCA. Numbered CT panels are
already covered. Avg_PCA_embedding score comparisons require no coordinate
export. The explicit old matched-allcell/55sample manual diagnostics stay
historical unless reopened. Bassez E/NE display must meanExpander/Non-expander,
not clinical responder categories. No coordinate extension or gene feature
localization has been implemented by this checkpoint.

Zeroimp code-history audit adds compatibility evidence, not new provenance:
canonical LISI normalization/clamping dates toJuly15, before NASAug25–27;
parsed relevant score/zeroimp functions match current definitions, apart from
ANOSIM execution parallelism with the same statistic. All31required and6
multLN settings are present in every NASsummary; missing8settings are only
try-wrapped conditionalmultRepl branches. Actual producer/runtime and current
count/label compatibility remain unproven until the bounded science audit.
No9zeroimp reruns or LISI reinversion are justified by current evidence.
Machine evidence `/private/tmp/ecoda-zeroimp-policy-audit-v1/policy_evidence.json`,
SHA256 `801c827a3d425aa224d802864cf29ee954153ab0146a7a3ae526f13ecd7f154a`.


### Continued execution after user withdrew sleep handoff — 2026-11-09 10:13Z

The user will leave the workstation awake and explicitly authorized continued
work. The existing single durable waiter remains active; no second waiter or
automation was created. Durable status is RUNNING, not terminal or released.

All32 method-array roots and32watchdogs are now recorded: benchmark18+18,
uncorrected7+7, corrected7+7. Corrected roots are4665615/16,4665617/18,
4665619/20,4665621/22,4665623/24,4665625/26,4665627/28 (array/watchdog pairs).
All three contexts passed their H5AD and R-environment preflights and reached
method submission. Corrected worker execution is not yet proven by logs;
submission/queue waiting is distinct from execution. No additional live queue
or accounting query was made. Aggregate --wait IDs may appear only at return;
collect every emitted ID before the sole terminal accounting inspect.

Benchmark watchdogs now reportOK for MOFA, Stephenson ordinary/BL/CT
pseudobulk family, Bassez HR pseudobulk parameters, Stephenson composition,
scITD, transformation and zero-imputation groups. The only observed failed
watchdog remains GloScope parameter task4665516_2. These are worker-level
observations pending final artifact/terminal review, not new publication
approvals. Independent arrays continue; no healthy work was cancelled.

P6 compact projection/mirroring/intake package14e6749 is committed and pushed;
no production projection/figure run occurred. A separate agent is implementing
the bounded existing-feature PCA rendering operation in the existing entry
point and helper. An independent reviewer is checking integration with the
mandatory run-scoped terminal auditor before any launch; the entry point's
own --validate check alone is not sufficient evidence of that integration.
A second implementor is preparing the finite106-file/209-logical-row legacy
science evidence reader from the recorded proposal. It performs no fits or
new scores and cannot approve historical provenance by itself. These code
packages remain separate from the immutable active main wave.


### Reporting renderer and terminal audit reviewed — 2026-11-09

The bounded render-pca implementation uses34 explicitly approved source rows:
authorHR, defaultHVG2000 pseudobulk and dataset metadata for11cohorts, plus
BassezHVG3000 for the existing fixed-marker diagnostic. It generates the
existing11benchmark_pca PDFs, Bassez diagnosticPDF/PNG and six strictly bounded
aggregateCSV schemas on HPC. RNA result features remain on HPC; local intake
receives finished figures/aggregate metadata only. Coverage is validated from
current author counts and metadata, not hardcoded33types or old55sample data.
Bassez E/NE display consistently meansExpander/Non-expander with a Clonal
expansion legend, including activeFigure2B/Supp19/presentation. Full-space
five-metric summaries remain distinct from two-dimensional plot coordinates.
No production figure has been regenerated by these code tests.

Default-summary PB selectors now use the defaultHVG2000 arm; realPB10 remains
in the parameter screen. AvgPCA remains an explicit diagnostic selection.
Historical unnumbered HiTME/HVC/Pelka RNA viewer prints are manual-only.
Maintained numberedSupp21 uses authorHR counts and remains ACTIVE: review
caught and corrected an accidental gating of that figure before commit.
Existing filenames/numbering and the protectedRplots.pdf were preserved.

The mandatory genericrun auditor now recognizes only the exact projection and
PCA-render purposes. Shared immutable source/runtime, terminal, selection and
one-scalar-job accounting checks precede source-pinned read-only operation
validation; unknown/partial purposes fail closed. Ordinary benchmark auditing
is unchanged. Focused real generic-entry tests, the full existing audit suite,
actual synthetic11cohort PCA rendering and independent review pass. Standalone
validation explicitly rejects wrong/failed terminal state or RUN_ID without
changing previous status/terminal/owner bytes. These are future-source changes;
active9ca workers and run manifests remain immutable.

A finite metadata-only header observation resolved the remaining empty Feather
ID-column proposal fields:46files,2,963,340bytes total, pinnedMD5/size/stat
checked before and after; Arrow schemas/pandas index metadata only, ZEROrecord
batches read. Every file explicitly declares__index_level_0__. Evidence
`.gate/p6_legacy_evidence_20261109_v1/feather_headers.observed.json`, SHA256
`68fb54e7f172e34ad987f0e75db42080aa911fdbccef8b972fb6237980694d2f`.
New versioned proposal `/private/tmp/ecoda-p6-legacy-science-controls-v4-s1fv7jt8/`
preserves all earlier controls. Its control SHA256 is
`14b6292938e97491a2c34062c89f2937924d5542d447f9366241e05737e57aee`;
rows209 SHA256`232716792fbbf3153ce6c07a8679fb84ed5dda5f3790eda4c4438395cbe0286d`.
Only46explicit ID fields changed, based on declared metadata;106physical files
and209logical rows remain unchanged. It is NOT a launched or science-approved
control. The separate evidence reader is under review and still needs its
minimal pinned execution/terminal contract and prior full-H5AD preflight
witness binding before launch. No model/preprocessing/score recomputation or
newgate submission occurred during this reporting code work.


### Legacy evidence reader independently accepted — 2026-11-09

The two-file evidence-only R reader and its focused tiny-H5AD/Feather/RDS
fixture passed independent review. Exact scope106physical files/209logical
rows/9authorHRmetadata pairs is fail-closed; declared sample-ID fields and
header witnesses are mandatory. Selected payloads are externally MD5-checked
before/after reads on HPC, one cohort at a time. Existing H5ADs supply only
configured obs columns; serialized resultRNA fields may materialize on HPC
but are neither exported nor retained across cohorts. Stored NULL shuffle
labels and matching author features/distances are preserved; onlyscITD may
drop IDs and every drop requires explicit review. Current author rawcounts
+.5/CLR/Euclidean distances and metadata are checked, with selected timing
limitations separated from unrelated old timing rows. No fits, score
calculations, historical producer/runtime claims or scientific approvals are
created by the reader. Current H5AD proof/size checking is explicitly
conditional on a separately bound previous full-content validation witness.

Focused fixture and independent run bothPASS; actualv4/header metadata
contract alsoPASS without payload reads. The minimal one-job controller and
terminal-audit integration is a separate pending package; controls are still
proposals, and no evidence job was launched. Renderer/audit checkpoint0b8729c
is pushed. Active main9ca source and all completed outputs remain unchanged.


### CPU-node reboot allocation delay and continued progress — 2026-11-09 10:44Z

A bounded run-owned log/status check at10:41:36Z found no additional error
signatures: benchmark7watchdogsOK/1FAIL (sameGloScope4665516_2),44method/
watchdoglogs; uncorrected20logs with actual output progress (e.g.JoanitoMrVI
completed its50epochs and saved distances); corrected0methodlogs despite
all7methodarrays submitted. CorrectedR-environment log4665565 and input
preflight statuses areOK. These are worker-level observations, not terminal
accounting, completed publication or a claim allthree contexts are executing.

To distinguish scheduler waiting from a broken controller/dependency, a SINGLE
targeted scontrol read of corrected4665615,4665619,4665623 was made (no repeat
squeue/sacct poll). Jobs arePENDING with ReqNodeNotAvail/UnavailableNodes;
preparation hasDependency=null. Node inventory shows all configuredEPYC-7742
shared-CPU nodesDRAIN/REBOOT_REQUESTED, reasonRebootASAP since~10:12Z, with
NextStateUNDRAIN oncpu007. Existing allocations may finish before reboot.
cpu046 isIDLE EPYC-7763; cpu047 isDOWN/nonresponding, cpu048DRAINED. This is
a new external allocation constraint, not evidence of another pipeline defect.

The existing benchmarkCPU-model constraintEPYC-7742 is an approved timing
contract. An async user question is pending: temporarily also allowEPYC-7763
for still-pendingCPU work and explicitly record mixed-hardware timings, or
retainEPYC-7742 and wait for rebooted nodes. No resources, constraints, scopes,
IDs, sources, oldmanifests or outputs have changed. Do not treat elapsed time
as approval. Continue independent code/review work and the original single
durable waiter while awaiting the decision. GPU policy is outside this choice.

Reporting checkpoints0b8729c and3477954 are pushed. The thin legacy evidence
controller/terminal audit is being tested independently; a minimal existing
derived-runner obs-only lane for22savedcluster-label results is being
implemented separately. No new fits/PCA/Harmony/Leiden runs are authorized
for those derived results. No additional gate has been launched.

Current input-content candidate witness: ten selected inputs' existing
full-checksum stage5preflight records, source/runtime/inputmanifest, exact
status/log/proof hashes and integer file stats were observed read-only. All
ten currentmtime/ctime precede their immutable successful preflightrecord
mtime; ZEROH5AD payloadbytes read. Small observation
`.gate/p6_legacy_evidence_20261109_v1/input_preflight_witness.observed.json`,
SHA256`5f838063ed15ec2d58dc3b50a0d3de7affcdc7b392c787d332428ee12a639684`.
It remains unapproved and terminal-accounting-review-pending. A future
separate APPROVED_INPUT_CONTENT_WITNESS_ONLY receipt may bind the sole
terminalinspect/accounting and selected10preflightPASS even if the overall
parentFAILED, without approving that gate or making it a formal dependency.


### Saved-cluster derived implementation independently accepted — 2026-11-09

The minimal existing run_local_ecoda_derived.R lane and focused regression are
accepted. Current res50/Harmony2 operations require explicit six-field bound
inputrefs/SHA and dataset selection; top-level sample/label columns are sole
authority (including Bassezsample and GongSharmaspecimen.specimenGuid), with
AND(enabled+benchmarkview) participation and no canonical-path fallback. It
reads only savedobs cluster labels throughh5py; no RNA/countlayers/obsm or
new PCA/Harmony/Leiden/model fit. Active categories only are tallied before
+.5CLR; unused categorical/global-zero padding is removed. Sample order and
production missing-cluster semantics are preserved; wholly missing samples
fail, and currentStephenson must have53samples. Biological labels affect only
evaluation, verified by label-permutation invariance of counts/features/
distances. Existing counts-bundle implementation supplies allfive scores.

Expected production result scope is11res50+11Harmony2 standalone bundles in
two fresh analysis roots. These are new distance/score derivations from saved
labels, not rerun embeddings/clustering or claims of old-derived RDS reuse.
Old20CSV count slices remain immutable compatibility candidates; this lane
does NOT claim they were validated/reused. Oldraw/snapshot modes require
explicit historical_inputs=true. Cell-depth extraction/snapshot schema and
its181draw schedules are unchanged and are not run by this new lane.

Independent actual obs-only fixture, existing cell-depth R regression and
stubbed cell-depth gate allPASS. Current config confirms exactly11enabled
benchmark datasets and top-level column authority. ProtectedRplots.pdf SHA
remains unchanged. No production run occurred: a pinned controller/common
terminal audit and prior full-content input-witness binding remain required
before these22outputs can be generated or approved.


### Legacy audit controller and terminal contract accepted — 2026-11-09

The dedicated legacy_analysis_audit.sh controller/worker and narrow generic
run-audit purpose route are independently accepted. API: --control ABS_JSON
--control-sha256 SHA. Typed requestecoda_legacy_analysis_audit_gate_v1 has
six explicit path/SHA refs: science_control(proposalv4), observation(legacy
inventory), zeroimp_schema, zeroimp_policy, feather_headers andinput_review.
It freezes bounded metadata under a fresh run root and emits exactly one
pinnedCPU64G1h scalar job; no newserializationgroup or NASsync. The accepted
Rreader alone opens selected resultRDS/Feather and configuredH5ADobs.

The separate input_review is mandatory: schema
ecoda_p6_preflight_witness_review_v1; decision
APPROVED_INPUT_CONTENT_WITNESS_ONLY; review_complete=true;
scientifically_approved_rows=0; formal_pipeline_dependency=false; exact
parent_terminal_state FAILEDorCOMPLETED; four path/SHA refs namedobservation
(the10input preflight candidate), gate_manifest, terminal_inspect andaccounting.
selected_preflight_tasks must have exactly10PASS rows withdataset,
fullhash_array_id, integer task_id, current_input_path, current_input_md5 and
current_input_size. It may accept only those currentinput-content witnesses
from a terminalFAILED parent; it does not approve that parent gate, oldlegacy
results or historical producer/runtime identity. Pending candidate alone is
rejected. No such final production review exists yet whilemainv2runs.

The actual physical manifest ish5ad_preflight.tsv with11rows, Stephensonfirst;
selected10 task IDs are matched to their exact dataset/view/path line. Review
caught a mistaken input_artifacts.tsv fixture assumption and the invented
terminal vocabularySUCCEEDED; both were corrected before acceptance using
the real durable1.0 gate/inspect accounting schema. Complete rawone-query
accounting and exactgateIDunion are retained, including unrelatedFAIL rows.
Current integer H5ADstat and original record/proof/status/log/source/runtime
bytes must match the prior witness; mtime/ctime must not postdate successful
fullhash recordcreation. This reuses previously validated content rather
than rehashing largeH5ADs solely for repair/validation.

Both independent focused controller/generic fixturesPASS, including pending
review/changedinput/wrongtask rejection, correctFAILED/COMPLETED cases,
failed-worker emittedID and mismatchJSON preservation, and read-only failed/
foreign-terminal validation. Syntax/whitespacechecksPASS. No scientific
approvals or original file/owner/proof changes are made by this controller.
No production gate has been prepared/launched: first finish and terminal
inspectmainv2, review selectedpreflight evidence, seal/stage exactcontrols
and create the next committed immutable snapshot. Saved-cluster runner
checkpoint781b166 is committed; active9ca workers remain unchanged.


### CPU decision resolved — retain EPYC-7742

The user answered the scheduling question: **Keep EPYC-7742; wait for nodes
to return.** Retain the existing benchmarkCPU model/threads/memory/runtime
contract for queued and recovery method jobs. No EPYC-7763 expansion, job
cancellation, resource adjustment or re-submission is made. Existing running
work and the original single durable waiter continue. Future one-CPU evidence
operations do not supply benchmark fit timings and retain their reviewed
operation-specific resource contracts. Do not ask again or infer a different
CPU preference from continued maintenance delays.

Reviewed checkpoints781b166 (saved-cluster runner) andf1be3fa (finite legacy
evidence controller/auditor) are pushed. All scoped local implementation
packages currently pass focused checks/independent review; remaining work is
execution/terminal review, the targeted actual failed-row recovery, and P6
source compatibility/derived/reporting runs with their recorded prerequisites.
The tracked working tree is otherwise clean; external-agent untracked files
remain preserved and excluded.

The existing durable waitprocess(session36153) has one asynchronous output
consumer(cell977). This is not a second gatewait or scheduler poll. Observe
completion through that consumer; do not attach another consumer directly to
session36153 or start another durablewait command.


### User clarified hardware policy: benchmark-only timings; batch unrestricted

The user explicitly clarified: **only benchmarkdatasets/views are used for
timing; anybatch uncorrected orcorrected analyses can run on anyCPUorGPU.**
This narrows the preceding retainEPYC-7742 decision to benchmarkviews ONLY.
Keep benchmarkCPU/GPU model constraints for timing comparability. For the
exacttwo batch views remove hardwareMODEL constraints; preserve eachmethod's
CPU-vs-GPU dispatch, accelerator allocation/VRAM, threads, memory, time limit,
scientific inputs/methods and pinned runtime. Do not ask again for this
already-authorized batch flexibility. No environment or immutable source is
hotpatched; a minimal dispatcher fix for futurewaves is being implemented
and reviewed separately, including watchdog/OOM flag propagation.

A SINGLE explicit pending-task inventory of the14batch arrays (not aperiodic
scheduler poll) found only18pending CPUtasks at this point; otheroriginal
batch tasks have leftPENDING. Exact selected scope:
4665623_7/8 (correctedMrVI),4665625_1..8 (correctedPILOT),4665627_1..8
(correctedQOT). All werePENDING, EPYC-7742,16CPU,128G in the frozen inventory.
No benchmark IDs, running/completed jobs, newmethod rows or cancellations
are included. Source remains9ca; original runmanifests/proofs are unchanged.

The proposed operational amendment clears onlySlurmFeatures on those
individual pending task IDs. Recheck each exacttask's owner/run/source/log
identity, state and resource contract immediately before updating; if it
has started, leave it untouched and recordSKIPPED_RUNNING. Preserve every
other resource and jobID; capture before/result/after metadata and report any
failures without broadresubmission. Installed scontroldocumentation explicitly
supports blankFeatures= and task-specificArrayJobID_TaskID updates.
Selection `.gate/main_v2_batch_hardware_20261109_v1/selection.proposed.json`,
SHA256`fefa835e72a71334ef0182e077af9ced87f7f06d1f8157d6367aeb6f0f0cc69e`,
18tasks. No update has yet occurred by this checkpoint; record outcome next.
Existing9ca watchdogs retain original submissionflags, so an actualfuture
OOM retry is assessed under this same user policy if it appears; do not
modify activewatchdogs or oldrecords.


### Batch flexibility applied; future source accepted; startup failures separated

The approved18task scheduler amendment completed: all18UPDATED_FEATURES_ONLY,
zero benchmarkjobs changed and zero newjobs/cancellations. Each task was
verifiedPENDING/ownedbyhalterc, exactcorrected run/source/log path,16CPU/128G/
12h/sharedcpu immediately before clearingFeatures; every other recorded
resource and identity field matched afterwards. Current9ca source and old
manifests/proofs remain untouched. Remoteappend-only journal/result are under
`_ecoda_inputs/20261108-confounding-subsets-v1/main-v2-batch-hardware-amendment-v1/`.
Local result `.gate/main_v2_batch_hardware_20261109_v1/result.json`, SHA256
`0ca85a0c907f6fd9dc95491709d941643891640f5e36c83dcb3916e04de9da24`.

The future dispatcher change is independently accepted: onlyexactbatchviews
filter model--constraint flags aftermethod/resource classification. Both
initial sbatch and capturedwatchdog/OOMretry receive identicalfilteredflags.
Benchmarkpins, CPU-vs-GPU policy,explicitVRAM/GPUcount,NV,partition,threads,
memory andtime are unchanged. ActualBash3.2 test covers45view/method/class
cases, initialsubmission and actualretrycommand; existingcorrected/OOM
regressionPASS, syntax/whitespacePASS, protectedRplots.pdf unchanged. Old9ca
watchdogs were not edited and may retain their originalpin on an actualretry;
handle only a realpendingretry under the user's existingbatch policy.

At11:12–11:15Z corrected has actualexecution, includingMrVItraining and
PILOToutputs. CorrectedPILOT4665625_2(Stephenson) and_4(Breast) COMPLETED0 on
cpu054(EPYC9754), with featuresnull, after the amendment. CorrectedMrVI_4
isRUNNING oncpu035/EPYC7742. This proves the thirdcontext passedstartup;
submission alone is no longer the evidence. Finalscience/publication review
remains pending.

A NEW incident is distinct from the earlierGloScope inputpath failure:
correctedQOT watchdog4665628 reports all8tasks4665627 FAILED0:53, nonOOM.
All8actual expectedstdout andstderr paths areabsent. Task1/raw4665756 on
cpu046 ran12:10:51–12:10:52clusterlocal (one second), RaisedSignal53,
featuresnull, norestarts. SuccessfulPILOTcomparison wascpu054, notcpu046.
cpu046 reportsIDLE; a directreadonlySSH mount/path probe timedout onport22,
so no daemonerrno/compute-path visibility proof was obtained. Currentnode
health/loginvisibility cannot provehistorical computeI/O. NoQOT code or
vendoredsource changed, and no QOTrecovery job was submitted.

InstalledSlurm is26.05.4. Exactofficialsource permits a misleading0:53 for
I/Osetup: ESLURMD_IO_ERROR4021 can survive stepstatusaggregation and be
formatted throughsignalbits. This supports—but does not establish—a batch
stdout/stderr startupI/O failure, rather than proving an actualRTsignal
sender or QOTalgorithm defect. Primarysource:
https://github.com/SchedMD/slurm/blob/slurm-26-05-4-1/src/slurmd/slurmstepd/mgr.c#L2028-L2037
andaggregation#L302-L306; formatter
https://github.com/SchedMD/slurm/blob/slurm-26-05-4-1/src/slurmctld/job_mgr.c#L5919-L5927 .
Exactdaemon error/errno is stillneeded forcausalconfirmation.

Small incident metadata `.gate/main_v2_batch_hardware_20261109_v1/qot_startup_incident.json`,
SHA256`2392af9c54fef05c0ae27dc2d442695f504e52b2c48d94fcce707295d2f47a58`.
Continuevalidindependent work and the single durablewait. The oldone-row
GloScope-only proposal is no longer a complete recoveryselection; freeze the
actualfailed/missing/invalid rowcatalog afterterminalinspection before any
newgate. Consider a run-scoped exclusion of suspectcpu046 only in recovery
while investigating; do not permanentlyhardcode nodeexclusions or silently
refit successful rows. Noaccountingquery hasbeen made for activev2.


### Explicit GPU flexibility for future batch MrVI accepted

Follow-up code/history review found that batchMrVI's CPU-only guard was a
previous queue/resource workaround(commit3531e67), not a model/scientific
requirement. Under the user's new anyCPU/GPU batch authorization, the minimal
implementation preserves auto→CPU to avoid unexpectedly moving work intoGPU
queues, but honors existing explicit --gpu-policy any/default for future
batchMrVI. No newCLI or automatic GPU migration is introduced. Benchmark
resource/timing policy and all current runningCPU fits remain unchanged.

The existing CUDAavailability/runtime check remains fail-closed. Scientific
MrVI setup, sample/technicalbatch keys, seed,50epochs and distance logic are
unchanged; only the permitted acceleratorselection changes. ExplicitGPU
count/VRAM,partition,time,threads andNV survive normaldispatch and OOMretry.
IncreasinghostRAM does not solve GPUVRAM exhaustion; no numericalbitidentity
claim is made for CPU-vs-GPU stochastic model output. Actualbatch device/GPU
name is logged. The newprovenance logger is exactbatchview-only, guarded by
CUDAavailability: it performs ZERO CUDAproperty calls on benchmarkviews, so
benchmark context initialization is not moved outside the method timer.

README,NOTES and onboardingREADME now state benchmark-only comparative
timings and unrestricted compatiblebatch hardware models, with CPUdefault
and explicitGPU option. Focused actualshell/AST device+allocation testPASS,
real watchdogVRAM/NV/retry regressionPASS, actualrequested-scope dispatcher
regressionPASS. Root independently reran/reviewed the finalfocused test and
timing-boundary correction; syntax/whitespacechecksPASS, protectedRplotsSHA
unchanged. This changes future committedsources only; no newGPU job, backend
change to running work, or recomputation of valid results occurred.

Inparallel, bounded controller fixes are being tested for exactselected
recoverycontexts/correctedmethodmatrix and eight metadata-only copies. The
current mainwave must still terminalinspect before sealing actualfailed-row
scope and launchingrecovery. Its originalone-rowGloScope proposal remains
incomplete after the correctedQOT startup failures.


### Selected metadata-only adoption independently accepted

New helper receipt schemaecoda_stage5_metadata_adoption_v1 supports exactly
the selected nonempty batchcontext subset and its explicitselected_metadata
dataset/view union. targets,sources andselected_metadata must have identical
context keys; targetbindings determine the exact scientificdataset/view
selection. Artifacts aremetadata ONLY: noPB caches, timing, benchmark copies
or unusedcontexts. The originalecoda_stage5_prepared_adoption_v1 contract
remains exactly23payloads/all3contexts. No fake sourceTASK0 provenance is
accepted: futureselected8 copies use originalv1 EXPORTED task1..8 evidence,
whose bytes already passed the23copy validation.

Current runtime/config/H5ADidentity, payload/proof/owner/status and fresh
target checks are retained. Sourcebytes/modes/statuses are untouched; new
truthful copy records/NOOP statuses/owners reachOK only after actualmetadata
semanticvalidation. The8Feather fixture proves PENDING0 and NO metadataexport
Slurm jobs, with no unusedbenchmark/uncorrected roots. Missing/extra/foreign
contexts/payloads, corruption, alteredinput, TASK0, existingtarget andforeign
ACTIVE owner failclosed. Review corrected one negativefixture to refresh its
mutatedbindingproof and assert the exactH5AD-identity rejection; production
comparison was already intact.

Newselected8 test and independent correctedfixturePASS. Both original23
metadata and actual6-cache/RDS/BL/timing/no-refit regressionsPASS. Twofile
implementation/test package accepted, sourceSHA
`6e31fdf0cd86ced0ead9930c265b755be92f7af65c54ef7edb2bd5923fe8d287`,
testSHA`14dfb78b84f86bb3d45d97b4a653ba6639866d068b0b5df8f8f03a3d6fc96e8d`.
This is code acceptance; no productionreceipt/recoverycontext or newgate is
sealed/launched before the currentwave's terminalcatalog and review.

### SSH observation delay — status remains unconfirmed after lastgoodcheck

The next roughly30minute readonlywatchdog/error collector did not complete.
A separate minimalSSH date command stalled too. A fresh connection (noSSH
multiplexing) establishedTCP, verified the knownhostkey and progressed
through keyexchange/serviceaccept before timingout; localssh-add public-key
listing returns normally. No key/authconfiguration or accesspolicy changed.
This does not prove a compute failure or a stoppedgate; remote jobs were
launched durably and must be reconciled from real evidence whenaccessworks.

Only the two stalled LOCALreadonlydiagnostic SSHclients were terminated,
verified by exactcommand+collectorparent. ZEROHPCworkers or durablewaiter
processes were terminated. Original waitsession36153 and its sole output
consumer977 remain untouched. Do not start anotherwait/relaunchwrapper or
claim a freshRUNNING/COMPLETED state from silence. Noactivev2accountingquery
hasoccurred. Continuecode checkpointing and retry oneboundedtransportcheck
later; hostname/clusterpolicy remainsBamboo, without a Yggdrasilfallback.


### Exact selected-context recovery controller accepted

The six-file campaign/pending/accounting package is independently accepted.
Campaigns may select1..3 nonempty recognized UNIQUEcontexts and freshchildIDs;
all input/config/runtime/source/scope checks remain. Healthyunselected
contexts get NO childroot, jobs, owners or accountingprojection. A fifth
context field may carry an explicitcorrected-only methodmatrix; it is
frozen/hashed with the parent scope and forwarded through the existingCLI.
Newmetadata-only adoption is forwarded only to its declaredselected target
contexts; original23payload/all3context behavior remains strictunchanged.

Fresh-bound corrected partialselection now passes ONLY through its exact
run-owned methodmatrix: strict checksum/MD5/size/SHA/count, uniqueallowed
three-tuples, enabledANDview participation, exactcohortunion and ordered
dispatch equality. Existingfull7family CLI remainsfixed; TARGET_METHODS_SET
stays0. Immutablebinding/disjointfreshroot/no-force/no legacyrecovery/sync/
residualoverride and ownershipguards remain mandatory. Uncorrectedmatrix
support is NOT expanded without actualcatalog need.

Accounting capture/audit is limited to the exactselected children and their
complete schedulerIDunion, with original/staged orderedcontext identity,
sealedmatrix/childdispatch bytes and strictmetadata identities. It does not
invent a childmatrix producerrecord for older schemas that have only their
existing strict sidecar and metadata hashes. The existing single complete
rawaccounting query contract remains; failedoriginalgates are not approved
or used as formal dependencies by this code change.

Actual dispatcher/stub fixture proves exactly9scientifictasks: oneStephenson
GloScope parameter and8correctedQOT rows, with NOpseudobulk preparation,
uncorrected context or healthy method rows. SixstubIDs are2methodarrays+
2watchdogs+2aggregates; productionH5AD/R-environment preflights are separate
and must be counted explicitly in the eventualsealed launchspec. Foreign
derivativeowner failure remains isolated while its selectedpeer completes;
oldcanonical ownerbytes remain unchanged.

Final unchanged nine-taskfixturePASS; old56/full7 andsmall7/legacyprotection
regressionsPASS; campaignfixturePASS; accounting15testsPASS; independent
repeat of allfocusedchecksPASS. Syntax/AST/whitespacePASS and protectedRplots
SHA unchanged. NoHPCoperation happened in these tests. Selectedmetadata
helper checkpointbb16dd0 is committed; actualv3rowcatalog/wrapper/resources
remainUNSEALED untilmainv2terminalinspection and access recovery.


### User-assisted status handoff / possible Codex restart

The user reports the already-established PositronSSH session stillworks and
our queue appearsrunning. A NEWPositronSSH session failskex/resetbypeer,
matching the agent's newconnectiontimeouts. This corroborates a new-login
transport/service problem, not proof of stopped Slurmjobs or a Codexapproval
block. Localssh-agent publickeylistingworks. Freshper-method statuses after
11:15Z remain unverified by the orchestrator; no activev2accounting query.

Nativecomputer access toPositron timedout on itsfirstgetApp call; no terminal
command or jobchange was executed. CUAcompletion was recorded. The user
allowed window/screenshot access and reportsmacOS mayrequireCodexrestart
for the newpermission. They explicitly offered to pastecommands into the
workingSSH terminal. Nextstep is a boundedreadonlythree-context watchdog/
terminalstatus reader (noSlurmaccounting/submission) and pastedoutput. Keep
the existingPositronSSH connectionopen. NoChromeextensionpermission is
required for thiscopy/paste route.

Before anyapprestart: sourcecheckpointsbb16dd0(selectedmetadata) andf2df1c9
(selectedrecoverycampaign) are committed/pushed; scopedtrackedwork isclean.
Externalapplication/Word/editor files remainuntracked and untouched. Remote
Slurmjobs/persistentgate are independentofCodex; localtool sessions/monitor
attachment mayneedrecovery. Onresume firstcheck the originalsinglewait
process/session36153 and consumer977, thenreconcile existinggate identity
whenSSHworks. Neverblindlylaunch anotherwrapper or secondwait. Active
scientificsource remains9ca; no newsource hasbeenpulled/hotpatched remotely.
Actualrecoveryrows/wrapper/run controls remainUNSEALED pendingterminal
inspection andcompletefailed-row catalog. Minimumknownfailures remain
StephensonGloScope30 +8correctedQOTstartup failures.


### Restart recovery and user-assisted evidence — 2026-11-09 13:20Z

The user restarted Codex and supplied a 13:09:20Z queue/status snapshot.
Benchmark reports are eight OK and one FAIL out of eighteen watchdogs;
uncorrected has one OK report out of seven; corrected has one FAIL report
out of seven. The known failures remain Stephenson GloScope parameter task2
and all eight corrected QOT tasks. Missing reports are unresolved, never
inferred successful from absence in the queue. Dependency-pending aggregate
IDs in the supplied queue include4665594,4665612,4665633; include these and
all other actual emitted/retry IDs in the eventual single terminal inspect.

Native Positron access now works after restart. The orchestrator raised the
ECODA workspace and used only the existing Bamboo SSH terminal for bounded
read-only filesystem/log commands. The connection remains open at its prompt;
no submit/cancel/accounting/source mutation was performed. A fresh standalone
SSH attempt still timed out, corroborating the new-connection problem.
At13:20:40Z, uncorrected MrVI task4665560_3 log was around44/50 epochs
(last write13:18:27Z) and corrected task4665623_3 around24/50
(last write13:16:58Z). Both scientific fits are making progress. Watchdog
log discovery found only the same reporting families; other watchdog startup
or completion status still requires the final catalog. No additional method
failure is asserted without evidence. The durable gate status file is absent;
that absence alone does not establish failure or terminal state.

Local exact-process inventory after restart proves the original durable wait
and output-consumer processes are no longer present. Historical claims that
sessions36153/977 remain attached are superseded by this restart observation.
No duplicate waiter or wrapper was started. Reconcile the existing remote gate
identity when fresh SSH is restored, then restore exactly one durable wait;
the active remote Slurm jobs must not be relaunched. Current9ca sources and
all completed artifacts remain unchanged. Recovery selection is still unsealed
pending complete terminal/accounting/artifact evidence. No blanket rerun or
new HPC job was launched during this check.
