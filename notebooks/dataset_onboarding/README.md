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
plan: `.agents/plans/dataset_onboarding_and_debug_overhaul.md`.

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

## Fast Workflow: Generate Subsets on HPC & Render Locally on Mac

Reading 53 GB `.h5ad` files over a network NAS mount incurs heavy SMB random-seek latency. To achieve fast, sub-second loads and renders on your Mac:

### Step 1: Generate diagnostic subsets on HPC scratch (~1–2 minutes)
On the HPC login node or compute node (where the ~195 GB datasets sit on BeeGFS scratch):

```bash
cd "${HOME}/ECODA_paper" && git pull
./notebooks/dataset_onboarding/run_subset_hpc.sh               # runs across all 9 datasets
# or for a single dataset:
./notebooks/dataset_onboarding/run_subset_hpc.sh --only alzheimer
```

This runs `create_subsets_hpc.py` to perform the full-file `count_sanity_check()`, extract full metadata summaries (`<Name>_meta.json`), and write lightweight diagnostic `.h5ad` subsets (~15–40 MB each, ~200 MB total) into `${HPC_SCRATCH_DIR}/_downloads/subsets/`.

### Step 2: Pull the lightweight subsets to your Mac (~200 MB)
From the repo root on your local Mac:

```bash
mkdir -p data/new_dataset_checks/subsets
rsync -avP bamboo:scratch/ECODA_paper/_downloads/subsets/ data/new_dataset_checks/subsets/
```

### Step 3: Render check notebooks instantly on your Mac
With the local subsets in place, notebooks load in <0.2s and render completely in ~5–10s without touching the network:

```bash
PATH="$PWD/.pixi/envs/default/bin:$PATH" quarto render notebooks/dataset_onboarding/dataset_check_Alzheimer.qmd
```

*(Note: If local subsets are missing, the notebooks automatically fall back to opening the full file from the mounted NAS in backed mode).*

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

## Summary Table & Diagnostic Results

| # | Dataset | Study (PMID) | File (`JooM_2025_41097818/output`) | Size | md5 | Count Check | n cells / samples / CTs (paper) | Bio Condition (PILOT-GM-VAE) | Bio Condition (ECODA) | Batch Condition (PILOT-GM-VAE) | Batch Sequencing (ECODA) | Batch Sample Prep (ECODA) | Technical Batch Candidates | Suitable for Auto-Annotation | Recommended Use | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Alzheimer | Gabitto 2024 Nat Neurosci (42486312) | `SEAAD_Alzheimer.h5ad` | 53.2 GB | `c2ad4c584f31f40e8aae0b32608e8146` | NOTE (X: log-normalized; no raw.X) | 1,395,601 / 83 / 18 | `Cognitive status` | `Cognitive status`, `ADNC`, `Braak stage`, `CERAD score` | None | `assay` | `tissue_type`, `PMI` | `assay`, `tissue_type`, `PMI` | No | `benchmark` | Checked & rendered |
| 2 | Breast cancer | Wu 2021 Nat Genet (34493872) | `BreastCncr_processed.h5ad` | 28.9 GB | `8b28a349c2c3638ddbfb3946a32d12ba` | PASS (raw.X: raw integer counts; X: log1p) | 714,331 / 126 / 10 | `disease` | `disease`, `sample_preservation_method`, `suspension_dissociation_time`, `suspension_dissociation_reagent` | None | `assay`, `sequencing_platform` | `sample_preservation_method`, `suspension_dissociation_time`, `suspension_dissociation_reagent` | `assay`, `sequencing_platform`, `sample_source`, `suspension_dissociation_time` | Yes | `batch-effect` | Checked & rendered |
| 3 | Covid-19 PBMC | Ren 2021 Cell (33657410) | `Covid19_Ren2021.h5ad` | 30.4 GB | `ae2fab89414914b6001879c01f822381` | PASS (raw.X: raw integer counts; X: log1p) | 993,171 / 151 / 10 | `CoVID-19 severity` | `CoVID-19 severity` | `Single cell sequencing platform`, `datasets` | `Single cell sequencing platform` | `Sample type` (TBD) | `Single cell sequencing platform`, `City`, `datasets`, `Sample type` | Yes | `batch-effect` | Checked & rendered |
| 4 | Diabetes (mouse) | Hrovatin 2023 Nat Metab (37697055) | `diabetes.h5ad` | 4.1 GB | `38189a381bad630fa39ce2d7ad3a0855` | PASS (raw.X: raw integer counts; X: log1p) | 264,235 / 52 / 13 | `disease` | `disease` | `dataset` | `assay` | `dataset`, `design` | `dataset`, `design`, `assay` | No | `batch-effect` | Checked & rendered |
| 5 | Kidney (KPMP) | Lake 2023 Nature (41648348) | `Kidney_KPMP.h5ad` | 2.75 GB | `36ceb02ba23c559f80625ec7bef6884f` | PASS (raw.X: raw integer counts; X: log1p) | 104,314 / 45 / 14 | `condition.l1` | `condition.l1`, `condition.l2` | None | `assay`, `library` | `tissue_type`, `region.l1` | `assay`, `tissue_type`, `region.l1`, `library` | Yes | `batch-effect` | Checked & rendered |
| 6 | Lupus PBMC | Perez 2022 Science (35389779) | `Lupus_Perez2022.h5ad` | 24.4 GB | `001658910686c61a5010da95b7b14a15` | PASS (raw.X: raw integer counts; X: scaled) | 1,263,676 / 261 / 11 | `Status` | `Status` | `batch_cov` | `batch_cov` | `Processing_Cohort` | `batch_cov`, `Processing_Cohort` | Yes | `batch-effect` | Checked & rendered |
| 7 | Lung | Sikkema 2023 Nat Med (37291214) | `lungatlas.h5ad` | 17.4 GB | `010cd8b233ac569b711ea0cbd80980be` | PASS (raw.X: raw integer counts; X: log1p) | 941,504 / 165 / 12 | `disease` | `disease`, `origin` | `dataset` | `dataset`, `platform`, `study`, `assay` | `dataset`, `study`, `origin`, `origin_fine`, `tissue sampling method`, `tissue dissociation protocol`, `anatomical region`, `donor status` | `dataset`, `study`, `platform`, `assay` | Yes | `batch-effect` | Checked & rendered |
| 8 | Myocardial infarction | Kuppe 2022 Nature (35948637) | `Myocardial_Infarc_2.h5ad` | 3.6 GB | `7431ae99250c99f11bf63e3034798af4` | NOTE (X: log-normalized; no raw.X) | 132,888 / 23 / 11 | `patient_group` | `patient_group` | None | `batch` | `sampleType` | `batch`, `sampleType` | Yes | `batch-effect` | Checked & rendered |
| 9 | Parkinson | Kamath 2022 Nat Neurosci (35513515) | `Parkinson.h5ad` | 30.5 GB | `f576bcf5eb28366aeaecff01c50fff34` | NOTE (X: log-normalized; no raw.X) | 2,096,155 / 97 / 11 | `disease` | `disease` | None | `assay` | `Brain_bank`, `tissue_type` | `Brain_bank`, `assay`, `tissue_type` | No | `benchmark` | Checked & rendered |
## Key Findings & Guidelines

1. **Experimental Design & Attribution:**
   - PILOT-GM-VAE author selections, ECODA biological conditions, and sequencing vs sample preparation batch covariates are consolidated in the master Summary Table above.
   - For datasets like **Lung Atlas**, PILOT-GM-VAE evaluated `disease`, whereas ECODA additionally evaluates `origin` (tumor vs normal) and extensive consortium-level technical batch structures (`dataset`, `study`, `platform`, `assay`).
2. **Dual Evaluation Architecture:**
   - Every onboarding report evaluates both **expression-level LISI batch mixing** (on unintegrated PCA embeddings) and **cell-type compositional variance partitioning** (quantifying the fraction of total inter-sample variance explained by each biological vs technical covariate across cell types, following Sikkema et al. 2023 Fig 4a).
3. **Count Layer Integrity:**
   - Datasets from CellxGene / CZ CELLxGENE (Breast Cancer, Covid-19 PBMC, Diabetes, Kidney KPMP, Lung, Lupus PBMC) store raw integer count matrices in `adata.raw.X`, while `adata.X` holds log1p-normalized / scaled expression. `onboarding_utils.locate_counts()` automatically detects and routes raw integer counts from `raw.X`.
   - Datasets where the downloaded file only provides log-normalized expression (Alzheimer, Myocardial Infarction, Parkinson) are clearly noted with `NOTE (X: log-normalized)` status.
4. **Assay Covariate & Level Distribution:**
   - In single-platform atlases like **Kidney KPMP** and **Parkinson**, the `assay` metadata column is genuinely single-valued (`10x 3' v3`) across all 104k and 2M cells, as authored by the consortium. Technical batch variance in these cohorts is driven by `region.l1` (Cortex vs Medulla), `library`, and `Brain_bank`.
   - In multi-platform atlases (**Lung**, **Covid-19 PBMC**, **Breast Cancer**, **Diabetes**), `assay` and `platform` contain multiple distinct chemistries (`10x v2/v3`, `BD-Rhapsody`, `Smart-seq2`, `NovaSeq 6000` vs `HiSeq 4000/3000`).
5. **Auto-Annotation Suitability:**
   - **Unsuitable (`No`):** Brain tissue cohorts (**Alzheimer**, **Parkinson**) lack hematopoietic/immune marker lineages for HiTME/scATOMIC, and mouse diabetes (**Diabetes**) requires ortholog symbol conversion.
   - **Suitable (`Yes`):** Solid tumor & PBMC cohorts (**Breast Cancer**, **Covid-19 PBMC**, **Kidney KPMP**, **Lung**, **Lupus PBMC**, **Myocardial Infarction**, **Debug**) have well-defined immune and stromal lineages.