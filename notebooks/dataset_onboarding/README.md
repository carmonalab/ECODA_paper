# Dataset onboarding — PILOT-GM-VAE study (Phase 5)

Reviewer-requested additional datasets from Joodaki et al. 2025, *Brief
Bioinform* 26(5):bbaf547 (PMID 41097818). Author-provided files are downloaded
to the NAS folder `JooM_2025_41097818/output/` and checked with one
`dataset_check_<Name>.qmd` notebook per dataset before `datasets.json`
registration (user approval required) and pipeline rollout.

Sources: Zenodo records 8370081 (Part 1, latest version of 7435911 → 7956950 →
8370081 — byte-identical files), 7957118 (Part 2), 14615923 (Part 3); CellxGene
for Alzheimer (SEA-AD) and Parkinson. Excluded: Kidney cancer, PDAC, MI(1)
(duplicate of MI-2), follicular lymphoma. Full feasibility table (colors,
comments): `new_datasets_to_implement.md` (Excel source
`/Users/christianhalter/Desktop/ECODA_PAPER_DATASETS.xlsx`). Implementation
plan: `.kilo/plans/archive/1786899069337-onboard-new-datasets-phase5.md`.

## Download (sequential, md5-verified — user runs)

```bash
./notebooks/dataset_onboarding/download_datasets.sh            # all datasets
./notebooks/dataset_onboarding/download_datasets.sh --only breast  # one entry
```

Progress + md5s are logged to `notebooks/dataset_onboarding/download_log.md`.
Mac disk is tight (~108 GB free): the script selectively extracts only the
needed h5ads out of the tars and deletes the tars afterwards; run one download
at a time.

## Render a check notebook (needs the downloaded file + NAS mounted)

```bash
# pixi default env on PATH so quarto finds the python3 kernelspec
PATH="$PWD/.pixi/envs/default/bin:$PATH" quarto render notebooks/dataset_onboarding/dataset_check_Alzheimer.qmd
```

Outputs (plots/feathers/csv) go to `data/new_dataset_checks/<Name>/` (gitignored).

## Helpers

- `onboarding_utils.py` — shared Python helpers (Pixi **default** env, scanpy
  1.12.2; NOT the R notebook loader): `locate_counts`, `count_sanity_check`
  (integer-VALUE check with epsilon tolerance — float-encoded CSR is valid),
  `obs_summary`, `candidate_col_detection`, `cells_per_sample_stats`,
  `paper_table_compare`, `confounding_crosstab` (bio × batch collinearity),
  `embed_and_umap_workflow` (precomputed or computed UMAP from RAW counts,
  unintegrated only), `write_metrics_input`.
- `onboarding_metrics.R` — standalone Rscript (/ default env): cell-level,
  per-cell-type `calc_lisi` separation (repo `src/utils/scoring_metrics.R`)
  on the **unintegrated PCA embedding** (never `X_pca_harmony*`), per-CT
  subsample caps + confounded-CT guards; validated on the NAS `_debug` view
  (`_debug_validation.py`, 2026-08-17).
- `_debug_validation.py` — T3.0 helper validation on the Joanito 5-sample
  `_debug` batch-effect view (bio `sample.origin` vs batch `Site`/`seqtec`).

## Summary table (fill as downloads + checks complete)

| # | Dataset | Study (PMID) | File (JooM_2025_41097818/output) | Size | md5 | Count check | n cells / samples / CTs (paper) | Batch candidates | UMAP verdict | Recommended use | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Alzheimer | Gabitto 2024 Nat Neurosci (42486312) | `SEAAD_Alzheimer.h5ad` | ~15–25 GB | verify on download | - | 1,395,601 / 83 / 18 | region / assay (10x 3' v3) | - | benchmark | pending download |
| 2 | Breast cancer | Kumar 2023 Nature (37380767) | `BreastCncr_processed.h5ad` | 28.9 GB | `8b28a349c2c3638ddbfb3946a32d12ba` | - | 714,331 / 126 / 10 | 10x v2/v3, sample prep | - | batch-effect | pending download |
| 3 | Covid-19 PBMC | Ren 2021 Cell (34767776) | `Covid19_Ren2021.h5ad` (from `Datasets.tar.gz`) | (in tar 36.35 GB) | tar `d105b52dbba38ac49c2ffe8b3cf34e24` | - | 993,171 / 151 / 10 | collection site | - | batch-effect | pending download |
| 4 | Diabetes (mouse) | Hrovatin 2023 Nat Metab (37697055) | `diabetes.h5ad` | 4.1 GB | `38189a381bad630fa39ce2d7ad3a0855` | - | 264,235 / 52 / 13 | different studies, 10x v2/v3 | - | batch-effect (conditional) | pending download |
| 5 | Kidney (KPMP) | Lake 2023 Nature (41648348) | `Kidney_KPMP.h5ad` | 2.75 GB | `36ceb02ba23c559f80625ec7bef6884f` | - | 104,314 / 45 / 14 | multi-site | - | benchmark + batch-effect | pending download |
| 6 | Lupus PBMC | Perez 2022 Science (42115607) | `Lupus_Perez2022.h5ad` (from `Datasets.tar.gz`) | (in tar) | (in tar) | - | 1,263,676 / 261 / 11 | - | - | benchmark | pending download |
| 7 | Lung | Sikkema 2023 Nat Med (42362693) | `lungatlas.h5ad` (from `lungatlas.h5ad.tar.gz`) | 17.2 GB (tar) | tar `0d0c97924f1b7a405b6ec3b55da02882` | - | 941,504 / 165 / 12 | study / platform (technical) | - | check-before-use | pending download |
| 8 | Myocardial infarction (2) | Kuppe 2022 Nature (41937210) | `Myocardial_Infarc_2.h5ad` | 3.6 GB | `7431ae99250c99f11bf63e3034798af4` | - | 132,888 / 23 / 11 | - | - | benchmark (note: no HR CTs) | pending download |
| 9 | Parkinson | Prashant 2024 Sci Data (39580497) | `Parkinson.h5ad` | ~25–40 GB | verify on download | - | 2,096,155 / 97 / 11 | brain region | - | negative control candidate | pending download |

`-` = to be filled by the check notebook after download (T3/T7). Paper numbers
from Zenodo record 15575593 (`datasets.pdf`, the study's own dataset
descriptions) + `new_datasets_to_implement.md`.