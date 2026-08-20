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

| # | Dataset | Study (PMID) | File (`JooM_2025_41097818/output`) | Size | md5 | Count check | n cells / samples / CTs (paper) | Primary Bio Condition | Technical Batch Candidates | Recommended Cell Type | Recommended Use | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Alzheimer | Gabitto 2024 Nat Neurosci (42486312) | `SEAAD_Alzheimer.h5ad` | 53.2 GB | `c2ad4c584f31f40e8aae0b32608e8146` | PASS (X, raw integers) | 1,395,601 / 83 / 18 | `Cognitive status` / `disease` | `assay`, `tissue_type`, `PMI` | `cell_type` / `Subclass` | `benchmark` (`not_suitable_for_auto_annotation`) | Checked & rendered |
| 2 | Breast cancer | Wu 2021 Nat Genet (34493872) | `BreastCncr_processed.h5ad` | 28.9 GB | `8b28a349c2c3638ddbfb3946a32d12ba` | PASS (raw.X, raw integers) | 714,331 / 126 / 10 | `disease` (normal vs cancer) | `assay`, `sequencing_platform`, `sample_source`, `suspension_dissociation_time` | `broad_cell_type` / `cell_type` | `benchmark` | Checked & rendered |
| 3 | Covid-19 PBMC | Ren 2021 Cell (33657410) | `Covid19_Ren2021.h5ad` | 30.4 GB | `ae2fab89414914b6001879c01f822381` | PASS (X, raw integers) | 993,171 / 151 / 10 | `CoVID-19 severity` | `Single cell sequencing platform`, `City`, `datasets`, `Sample type` | `majorType` / `celltype` | `batch-effect` | Checked & rendered |
| 4 | Diabetes (mouse) | Hrovatin 2023 Nat Metab (37697055) | `diabetes.h5ad` | 4.1 GB | `38189a381bad630fa39ce2d7ad3a0855` | PASS (raw.X, raw integers) | 264,235 / 52 / 13 | `disease` (T2D vs ND) | `dataset`, `design`, `assay` | `cell_type` | `batch-effect` / `benchmark` | Checked & rendered |
| 5 | Kidney (KPMP) | Lake 2023 Nature (41648348) | `Kidney_KPMP.h5ad` | 2.75 GB | `36ceb02ba23c559f80625ec7bef6884f` | PASS (X, raw integers) | 104,314 / 45 / 14 | `condition.l1` (CKD, AKI, Ref) | `assay`, `tissue_type`, `region.l1`, `library` | `subclass.l1` / `cell_type` | `benchmark` | Checked & rendered |
| 6 | Lupus PBMC | Perez 2022 Science (35389779) | `Lupus_Perez2022.h5ad` | 24.4 GB | `001658910686c61a5010da95b7b14a15` | PASS (X, raw integers) | 1,263,676 / 261 / 11 | `Status` (SLE vs Control) | `batch_cov`, `Processing_Cohort` | `cell_types` / `ct_cov` | `batch-effect` | Checked & rendered |
| 7 | Lung | Sikkema 2023 Nat Med (37291214) | `lungatlas.h5ad` | 17.4 GB | `010cd8b233ac569b711ea0cbd80980be` | PASS (raw.X, raw integers) | 941,504 / 165 / 12 | `origin` (tumor vs normal) | `dataset`, `study`, `platform`, `assay` | `ann_fine` / `cell_type_major` | `batch-effect` / `benchmark` | Checked & rendered |
| 8 | Myocardial infarction | Kuppe 2022 Nature (35948637) | `Myocardial_Infarc_2.h5ad` | 3.6 GB | `7431ae99250c99f11bf63e3034798af4` | PASS (X, raw integers) | 132,888 / 23 / 11 | `patient_group` (infarct stage) | `batch` (library batches) | `cell_type` | `batch-effect` | Checked & rendered |
| 9 | Parkinson | Kamath 2022 Nat Neurosci (35513515) | `Parkinson.h5ad` | 30.5 GB | `f576bcf5eb28366aeaecff01c50fff34` | PASS (X, raw integers) | 2,096,155 / 97 / 11 | `disease` (PD vs Control) | `Brain_bank`, `assay`, `tissue_type` | `cell_type` | `benchmark` (`not_suitable_for_auto_annotation`) | Checked & rendered |
| 10 | Debug (_debug) | Joanito 2022 Nat Genet (35773407) | `_debug_5samples.h5ad` | 38 MB | - | PASS (X, raw integers) | 6,000 / 12 / 10 | `sample.origin` | `seqtec`, `Site` | `cell.type` | `reference` | Checked & rendered |

## Key Findings & Guidelines

1. **Dual Evaluation Architecture:**
   - Every onboarding report evaluates both **expression-level LISI batch mixing** (on unintegrated PCA embeddings) and **cell-type compositional shifts** (sample-level CLR abundance distributions with Mann-Whitney/Kruskal-Wallis statistical testing and Benjamini-Hochberg FDR correction).
2. **Cell Type Harmonization:**
   - Evaluated using `ou.cell_type_harmonization_check()`. Datasets compiled from multi-study consortia (Lung, Covid-19 PBMC, Diabetes) exhibit clear distinctions between atlas-wide harmonized annotations (present across $\ge 80\%$ of studies) and study-specific sub-annotations.
3. **Balanced Diagnostic Subsetting:**
   - Fast local Mac rendering is powered by $\min(N_{\text{samples}}, 20)$ sample allocation ($500$ cells/sample target, ~15–40 MB `.h5ad` subsets) generated on HPC via `create_subsets_hpc.py` / `select_balanced_samples()`. This guarantees balanced coverage across biological conditions $\times$ batch candidates without collinearity dropping test rows.
4. **Heatmap & Visual Sizing:**
   - Cell-level LISI separation heatmaps (`ou.plot_separation_heatmap`) and compositional significance matrices dynamically scale figure widths and label margins, ensuring clear numeric annotations and standalone colorbars without horizontal compression.