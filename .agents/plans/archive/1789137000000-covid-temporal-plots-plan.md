# Covid temporal diagnostic plots

## Goal
Refine the Covid19 extended covariate plot and add time-binned ECODA/pseudobulk MDS and ANOSIM diagnostics to the existing batch-confounding notebook.

## Scope
- Exclude `PatientID` from the Covid19 extended covariate plot.
- Render `Sampling day (Days after symptom onset)` as a continuous viridis-colored bar with no one-legend-entry-per-day explosion.
- Use the persisted ECODA and `Pseudobulk_hvg2000` feature matrices, subset rows by disjoint time bins (`control` or day <= 15; day > 15 and <= 30; day > 30), recompute Euclidean distances per subset, and plot MDS colored by the primary Covid severity label.
- Save one two-method MDS panel PDF per time bin and one ANOSIM line PDF under a `covid19_time_bins` output directory.
- Keep artifact checksum/size validation; tolerate only the known local mirror PATH mismatch while requiring MD5 and SIZE matches.

## Non-goals
- No H5AD reads, preprocessing, feature normalization, batch correction, or HPC work.
- No changes to `datasets.json` or existing benchmark artifacts.

## Verification
- Execute only the affected Covid plotting code in the persistent Pixi R session.
- Verify bin counts, feature/sample alignment, non-empty PDFs, checksum-backed ANOSIM table, and visual readability.
