# Batch confounding diagnostic plots

## Goal
Create a standalone notebook that visualizes sample-level representation of biological labels across configured and extended batch covariates for the active uncorrected batch-effect cohorts.

## Scope
- Add `notebooks/batch_effect_analysis_uncorrected_batchconfounding_contingency.rmd`.
- Use the twelve authoritative `BATCH_EFFECT_DATASET_ORDER` datasets, including `Kidney_KPMP_full`.
- Reuse validated `data/batch_effect/uncorrected/metadata/<dataset>_sample_metadata.feather` files and the existing batch metadata loader contract. If the canonical `Kidney_KPMP_full` artifacts are absent, use only the existing `Kidney_KPMP` artifacts through an explicit warning/status-marked fallback.
- Produce one PDF per dataset in `plots/batch_effect_uncorrected_batchconfounding/`:
  - `<dataset>.pdf`: `datasets.json` `columns.batch` covariates.
  - `<dataset>_extended.pdf`: deduplicated onboarding secondary biological plus technical candidates, excluding the primary biological label.
- Normalize each covariate-specific stacked bar to 100% within each biological-label group and label every nonzero segment with its unique-sample count.
- Include missing values as `Missing/unknown`, show covariate-specific grouped stacks, and add conservative potential-confounding warnings when a covariate level occurs in exactly one biological group.
- Write `batch_confounding_summary.csv` beside the PDFs with one row per observed covariate level (or unavailable candidate), including sample totals, represented biological groups, availability, source status, and flags.

## Non-goals
- No changes to `datasets.json`, preprocessing, benchmark artifacts, correction settings, or HPC execution.
- No statistical tests or batch removal.
- No reconstruction from cell-level H5AD matrices.

## Verification
- Validate notebook syntax and helper functions locally.
- Run the notebook plotting path against the available metadata/debug artifacts without full-cohort HPC work.
- Confirm deterministic dataset names, plot files, summary rows, percentage totals, explicit missing levels, and confounding flags.
