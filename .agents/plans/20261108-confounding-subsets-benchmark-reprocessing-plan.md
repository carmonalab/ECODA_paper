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
source/runtime, exact selection, and durable-gate prerequisites. No scheduler
jobs have been launched under this plan.

Provenance: repository planning anchor
`5302671ad94556edcf9acccf372d2dc34121d714`; inspected HEAD
`70b6033`. Record subsequent checkpoint, documentation, implementation, and
execution commits separately. A dirty configuration is not proof that its
decisions have been applied to existing artifacts.

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

These counts are expectations from existing confounding evidence, not new
Bamboo acceptance results. Reconcile raw IDs, current sample-level metadata,
and the existing >=500-cell/QC rules before launch. Any unexpected discrepancy
returns to the user rather than changing the selection automatically.

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
  explicit child run audits followed by terminal inspect/review. A failed
  campaign preserves success evidence for later targeted recovery/release;
  it does not authorize consumption solely from a successful status string.
- Do not introduce a general workflow engine, distributed queue, or replacement
  orchestration framework. Shared publication is distinct from compute survival.

## 5. Ordered work packages and acceptance

| Package | Status | Work and completion criterion |
|---|---|---|
| P0: safety checkpoint | DONE | Six logical commits pushed to origin/master; remote head verified as `6a1726d4a79f74998c54e968bb045b452b9c93bd` |
| P1: documentation | DONE | NOTES, onboarding README, both AGENTS guides, Bamboo durable-skill examples and historical-plan pointer harmonized; scoped whitespace and skill validation passed; Git publication recorded below |
| P2: scoped contracts | IN_PROGRESS | Separate implementors assigned cohort selection and BL core/cache fix; derivative roots, failure isolation and remaining interfaces follow integration |
| P3: source reconciliation | IN_PROGRESS | Orchestrator handles bounded read-only Bamboo source/runtime/NAS preflight; targeted Bassez scope is already user-approved, subject to unchanged-input proof |
| P4: changed inputs | PLANNED | Generate and validate 16 batch-view derivatives, Stephenson benchmark derivative, and only the necessary Bassez repair; satisfy benchmark annotation contracts |
| P5: method waves | PLANNED | Run three independent concurrent contexts with exact deduplicated scope, durable evidence, focused recovery and root-specific sync |
| P6: evaluation/reporting | PLANNED | Complete five-metric scores, score both Kidney axes, refresh all affected maintained reports/figures from approved results, verify plots and document limitations/provenance |

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

### Active assignments

- Orchestrator `/root`: this plan, Git integration/push, bounded Bamboo
  preflight and run/release coordination; no pipeline implementation.
- Implementor `docs_checkpoint_audit`: documentation harmonization, no current
  plan edits or production code/config changes. P1 returned; now independently
  reviewing BL and cohort packages. The orchestrator completed the remaining
  protected coordination-prose edits through scoped approved execution after
  the subagent's editing tool stalled; no pipeline implementation by root.
- Implementor `batch_subset_audit`: cohort/loader/Stage 3 Python package returned
  for review; now owns new Stage 2 Bassez metadata-only derivative worker and
  its fixtures. No shell submitter or Stage 5 edits.
- Implementor `benchmark_dependency_audit`: BL core/worker package returned and
  independently accepted; now owns derivative root/input-reference integration
  in submitter/worker/sync/audit seams, followed by exact variant propagation.
  No downstream notebook changes. Previously completed BL files stay frozen
  until their separate commit is recorded.

Checkpoint exceptions: `.vscode/settings.json` and `Rplots.pdf` remain local
and uncommitted as agreed. The architecture draft's three existing whitespace
lines were preserved in its backup commit. Backup commits are not scientific
release validation. No full H5ADs or expression matrices were copied locally.

No gates prepared or launched under this plan. For every future launch append:
package, source/runtime/aux identities, host and resolved roots, dependency
reviews, exact wrapper/selection and expected row count, gate/child run IDs,
scheduler/watchdog IDs, wait/inspect/reviewer evidence, incident/recovery scope,
and approved output-manifest location.
