# Batch-figure confounding consistency audit — 2026-11-11

Status: READ-ONLY AUDIT COMPLETE; COHORT POLICY DECISION PENDING. User requests a double-check of Lupus and all six PCA/UMAP
cohorts. This is a read-only follow-up to the completed geometry plan. No
cohort, configuration, embedding, score, figure or existing checksum changes
are authorized by this audit. Keep the current outputs intact while assessing
whether the cohort definition matches the user's intended batch exclusions.

## Roles

Root orchestrates and performs operational read-only checks. `batch_subset_audit`
owns Lupus lineage, raw batch counts and a proposed finite metadata probe.
`benchmark_dependency_audit` owns the other five cohorts and figure/input joins.
Any agent with an unclear role asks Root before acting. No new compute jobs,
source edits, commits or pushes are delegated.

## Questions to resolve

1. Are extra Lupus labels real observed batches or unused legend levels?
2. Did excluded samples re-enter, or did selected samples retain cells from
   other raw batches under the documented whole-sample majority policy?
3. Do the other five cohorts satisfy their approved selections, including
   assay/time, COVID day limits, biological exclusions and named outliers?
4. If strict per-cell batch exclusion is intended, what cells and samples would
   remain, and which downstream results would need a targeted update?

## Initial reproduced finding

Both Lupus views contain 17 active raw batch labels. The seven named selection
labels contain 384,718 cells; 10 other labels contain 74,067 cells (16.14% of
458,785). These other labels contribute 54,191 of 100,000 displayed cells
because the display balances across all observed batches. They are real cells,
not unused categories. All 94 approved samples appear in the bounded display,
with no excluded sample IDs. Full population cross-tabs by sample and biology
are being checked before any policy-change recommendation.

The cohort selector evaluates filters on sample-level majority metadata, then
retains every cell of selected samples. The previous plan explicitly chose
whole-sample retention. This can meet the recorded sample policy while failing
a literal requirement that every retained cell belong to the named batches.
The distinction must be made explicit; do not hide it by renaming batch labels.

## Evidence and boundaries

- Final figure set: `plots/singlecell_batch_geometry_20261110_v2/`.
- Approved display control SHA-256: `4133c7af3609a2399e5c3afda3b77ce8209cc41c71084b805a15b6459387d970`.
- Stage 3 parent receipt: `.gate/stage3_noop_controls_20261109_v1/batch_views.json`,
  SHA-256 `c8e1d9b5b0a4b485ecdb228c706c6f8a3de65f0861fc57b148de68cd39bb05b4`.
- Whole-sample rule: `src/utils/py/cohort_selection.py`; original processing
  plan archived under `.agents/plans/archive/`.
- Read only bounded local display/metric artifacts and existing source-bound
  receipts first. Any further remote binary inspection uses the configured
  Bamboo interpreter and reads only required obs metadata, sequentially.
  Never transfer full H5ADs, expression, PC50 or neighbour arrays locally.
- Final report will distinguish exact full-population facts from capped-display
  observations. A change to datasets.json or the cohort selection policy needs
  explicit user confirmation after impact is reported.

## Completed findings

All twelve figures were traced to their approved inputs. The other five
cohorts match their approved sample lists, predicates and active batch
categories in both views:

| Cohort | Samples | Checks |
|---|---:|---|
| Alzheimer | 101 | Three Reference samples excluded; two assay categories. |
| Breast | 33 | v3 only; 12-hour / 15-minute groups; all 114,579 displayed cells checked. |
| COVID | 44 | d02=22, d09=12, d13=10; control-or-finite-day≤30 rule retained. |
| Joanito | 124 | Tumor-2 and LymphNode absent. |
| Stephenson | 110 | Sanger and all four named exclusions absent; Ncl's approved 18.87% Healthy exception remains. |

These checks replayed the full source-bound sample metadata and verified
input MD5s, approved allowlists and actual figure/display bindings. Capped
plots were not treated as independent full-cell censuses. Actual full-metric
batch counts equal observed plotted batch sets. Existing source contracts
provide the sample-consistency checks for strict selection fields.

A single additional read-only Lupus inspection on Bamboo read only the obs
index, Sample, sampleID, batch_cov and Status in each of the two pinned H5ADs.
No expression, embeddings or graphs were accessed; no remote files or jobs
were created. Before/after parent stats and native proof/record/owner bytes
agree. Full cell-ID hashes match the geometry. The views have identical full
sample/batch/biology aggregates: no extra or missing samples, no ID mapping or
biology discrepancies, and all source cell counts agree.

- Current Lupus: 458,785 cells, 94 approved samples, 17 raw batches.
- Outside the seven allowed labels: 74,067 cells (16.1442%), from 20 retained
  samples. Six of the ten additional batches have only one biological group;
  all ten fail the minimum-two-samples-per-group rule in this retained cohort.
- All seven allowed raw batches pass the 20% and minimum-two rules, using
  sample incidence rather than cell weighting. Samples may span raw batches.
- Strict seven-batch filtering WITHIN THE EXISTING 94-SAMPLE ALLOWLIST retains
  384,718 cells and all 94 samples (46 Case / 48 Healthy). The smallest sample
  has 1,189 cells, and none falls below 500. This is a metadata counterfactual,
  not an executed change. Any subsequent preprocessing must validate it anew.
- Display balancing gives the outside-label cells 54.191% of plotted points,
  compared with 16.1442% of the full population. Drawing-order shuffling does
  not fix this cohort-definition difference.

## Report and decision boundary

The report and two CSV tables are under
`plots/batch_figure_confounding_audit_20261111_v1/`:

- README.md: six-cohort findings, cause, implications and source hashes.
- Lupus_raw_batch_counts.csv: 17 batch rows with cell and unique-sample counts,
  biological-group fractions, whitelist status and balance-rule results.
- perSample_strict_impact.csv: 94 raw/canonical samples with original, retained
  and removed cell counts; exactly 20 samples would lose cells.

Actual obs-only observation SHA-256:
`9f3146f089032f299f7029ff81d72ad215714e78fa25d611f369aa0d0bed56f3`.
Root independently validated all table totals, sample IDs, counts and hashes;
review is `.gate/batch_figure_confounding_audit_20261111_v1/ROOT_audit.review.v1.json`.

Root asked the user whether to implement strict cell-level exclusion within
the existing 94 samples or retain the documented whole-sample policy. This is
a scientific cohort-definition decision; AGENTS.md requires explicit user
confirmation before datasets.json changes. No filter, input, numerical result
or existing figure has been modified, and no rerun is authorized yet.

If strict filtering is confirmed, preserve the exact current sample allowlist
AND remove off-list batch cells before canonical preprocessing. Do not filter
all 261 source samples and recompute majority winners first: that could admit
previously excluded samples. Scope the repair to the two Lupus batch views,
their dependent method results, geometry/LISI and downstream reports. Do not
select unaffected cohorts or benchmark results for recomputation. Merely
hiding or recolouring points would leave the fitted inputs unchanged.

The previous geometry review should have highlighted this distinction instead
of leaving the raw 17-batch versus seven-majority-group implication in the
technical ledger. Keep historical artifacts and reviews intact; the new audit
qualifies their applicability to the user's intended strict selection.
