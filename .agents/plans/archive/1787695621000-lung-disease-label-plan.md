# Lung atlas disease-label migration

## Scope

Switch the Lung atlas registry label from `origin` to `disease` consistently in the authoritative configuration and active onboarding documentation. Keep `origin` as an observed secondary metadata field because it remains present in the source data and is useful for tumor-versus-normal diagnostics.

## Changes

1. Update the `Lung.columns.label` entry in `datasets.json` to `disease`.
2. Update `notebooks/dataset_onboarding/dataset_specs.py` so `Lung.bio_col` and `Lung.registry_roles.label` are `disease`; retain `origin`/`origin_fine` only as secondary or batch candidates.
3. Update the onboarding role table and Lung diagnostic report to describe `disease` as primary and `origin` as secondary.
4. Check active Lung-specific references and generated audit evidence; do not rewrite historical/superseded plans.

## Verification

- Parse `datasets.json` and load the Python onboarding spec.
- Assert the two authoritative Lung label fields are `disease`.
- Search active Lung onboarding files for stale role declarations.
- Run the focused onboarding metadata validator against the local Lung audit if regeneration is available.

## Notes

The source dataset contains five `disease` values, including `chronic obstructive pulmonary disease`; no source metadata columns are renamed or removed.

## Derived audit follow-up

The checked-in source/config migration is complete. The local ignored
`data/new_dataset_checks/subsets/Lung_meta.json` was generated before this
role change and still records `origin` as the selected registry label; the
metadata validator detects that mismatch. Do not hand-edit generated evidence.
Regenerate the Lung full-file audit through the durable ECODA HPC gate before
using that local audit as registry evidence.
