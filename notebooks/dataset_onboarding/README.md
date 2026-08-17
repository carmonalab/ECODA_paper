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
plan: `.kilo/plans/1786899069337-onboard-new-datasets-phase5.md`.

## Download (HPC route first — user runs; Mac→NAS script is the fallback)

The Mac NAS mount is unstable and the Mac disk is tight (~108 GB free vs
~150–180 GB of files), so downloads go over the **HPC** into BeeGFS scratch
(`${HPC_SCRATCH_DIR}/_downloads/`, 1.1 PB, no per-user size quota) and a
login-node tail rsyncs them to the NAS folder
`JooM_2025_41097818/output/`:

```bash
cd "${HOME}/ECODA_paper" && git pull        # prerequisite: get the new scripts
# from the HPC login node (repo root):
./notebooks/dataset_onboarding/download_datasets_hpc.sh                  # all 8 keys (array job, 3 concurrent)
./notebooks/dataset_onboarding/download_datasets_hpc.sh --only breast    # one key (resumable)
./notebooks/dataset_onboarding/download_datasets_hpc.sh --sync-only <job-id>  # resume: gate + NAS sync tail only
./notebooks/dataset_onboarding/download_datasets_hpc.sh --login-node     # fallback if compute nodes lack egress
```

Compute-node egress is smoke-tested first (debug-cpu curl against real
object URLs — a bucket-root HEAD would 403 and falsely fail); on failure the
submitter automatically runs the login-node path (`nice -n 19` +
`--limit-rate`). Tasks are `curl -L -C -` resumable and verified per key
(Zenodo md5s; CellxGene files size-verified via HEAD content-length — no
`.h5ad.md5` sidecar exists and the S3 ETag is a multipart digest; computed
md5 recorded as informational); tar entries extract only the needed
h5ads and the tars are deleted. Progress + per-key md5s + resume commands are
logged to `notebooks/dataset_onboarding/download_log.md` (commit from the
Mac). The original `download_datasets.sh` (Mac→NAS, sequential) is kept as a
NAS-stable fallback only.

## Render a check notebook (needs the downloaded file + NAS mounted)

```bash
# pixi default env on PATH so quarto finds the python3 kernelspec
PATH="$PWD/.pixi/envs/default/bin:$PATH" quarto render notebooks/dataset_onboarding/dataset_check_Alzheimer.qmd
```

Each notebook opens the h5ad `backed='r'`, runs the full-file count sanity +
obs exploration, then **section 1.5 builds an in-memory sample-first subset**
(`subset_by_samples`, ~10k cells) and UMAP + per-CT LISI metrics run on that
subset — target wall time <2–3 min and ~1–2 GB peak RSS on the 64 GB Mac.

Outputs (plots/feathers/csv) go to `data/new_dataset_checks/<Name>/` (gitignored).

## Helpers

- `onboarding_utils.py` — shared Python helpers (Pixi **default** env, scanpy
  1.12.2; NOT the R notebook loader): `locate_counts`, `count_sanity_check`
  (integer-VALUE check with epsilon tolerance — float-encoded CSR is valid),
  `obs_summary`, `candidate_col_detection`, `cells_per_sample_stats`,
  `paper_table_compare`, `confounding_crosstab` (bio × batch collinearity),
  `subset_by_samples` + `SUBSET_CONFIG` (sample-first, RAM-bounded subsetting
  — select ~10–20 samples stratified by bio condition + round-robin over the
  batch candidates, read only those samples' cells via per-sample
  `.to_memory()` slices, cap per sample / per CT / overall ~10k (max 50k);
  precomputed unintegrated `obsm` arrays are re-sliced to the subset rows so
  UMAP/metrics stay in the original embedding space; diagnostic only — the
  HPC pipeline uses full data), `embed_and_umap_workflow` (precomputed or
  computed UMAP from RAW counts, unintegrated only), `write_metrics_input`.
- `onboarding_metrics.R` — standalone Rscript (/ default env): cell-level,
  per-cell-type `calc_lisi` separation (repo `src/utils/scoring_metrics.R`)
  on the **unintegrated PCA embedding** (never `X_pca_harmony*`), per-CT
  subsample caps + confounded-CT guards; validated on the NAS `_debug` view
  (`_debug_validation.py`, 2026-08-17).
- `_debug_validation.py` — T3.0 helper validation on the Joanito 5-sample
  `_debug` batch-effect view (bio `sample.origin` vs batch `Site`/`seqtec`),
  including the T3.1 `subset_by_samples` path (reduced config: 4 samples,
  per-sample cap 2000, target 5000 — checks per-sample slice reads, budgeted
  wall time/RSS, structural batch signal).
- `run_download_worker.sh` + `download_datasets_hpc.sh` — Phase-5 HPC download
  worker (sbatch array, one task per key: resumable `curl -L -C -`,
  Zenodo md5 + CellxGene size verification via HEAD content-length, selective
  tar extraction) + login-node submitter
  (egress smoke test, sacct gate, NAS rsync + md5 verify + log append).

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