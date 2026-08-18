# Onboard 9 new datasets from the PILOT-GM-VAE study (Phase 5)

## Implementation status (2026-08-18)

Implemented + committed (2026-08-17 – 2026-08-18):
- **T1 + T1.1: HPC parallel download + verification + NAS sync COMPLETE (2026-08-18)**:
  - All 8 download tasks (`alzheimer`, `breast`, `covid_lupus_tar`, `diabetes`, `kidney`, `lung_tar`, `myocardial`, `parkinson`) completed in BeeGFS scratch (`${HPC_SCRATCH_DIR}/_downloads/`).
  - Bugfixes implemented and committed:
    - Smoke test probed real object URLs (commit `52e785e`).
    - CellxGene size-verification vs HEAD content-length (commit `52e785e`).
    - Multi-pair tar canonical parsing fix (commit `eb28cb8`).
    - Submitter grep no-match fix under `set -e` (commit `f2c54fb`).
    - Rsync SMB permissions flag fix (commit `3ba440d`).
  - All 9 canonical `.h5ad` files (~195 GB) rsynced cleanly to `${NAS_SC_DIR}/JooM_2025_41097818/output/`.
  - NAS MD5 verification passed for all 9 files against expected/worker-recorded digests:
    - `SEAAD_Alzheimer.h5ad` (53,187,995,660 B, md5: `c2ad4c584f31f40e8aae0b32608e8146`)
    - `BreastCncr_processed.h5ad` (28,939,228,608 B, md5: `8b28a349c2c3638ddbfb3946a32d12ba`)
    - `Covid19_Ren2021.h5ad` (30,408,338,999 B, md5: `ae2fab89414914b6001879c01f822381`)
    - `Lupus_Perez2022.h5ad` (24,425,598,482 B, md5: `001658910686c61a5010da95b7b14a15`)
    - `diabetes.h5ad` (4,134,240,780 B, md5: `38189a381bad630fa39ce2d7ad3a0855`)
    - `Kidney_KPMP.h5ad` (2,755,120,874 B, md5: `36ceb02ba23c559f80625ec7bef6884f`)
    - `lungatlas.h5ad` (17,356,989,429 B, md5: `010cd8b233ac569b711ea0cbd80980be`)
    - `Myocardial_Infarc_2.h5ad` (3,605,875,880 B, md5: `7431ae99250c99f11bf63e3034798af4`)
    - `Parkinson.h5ad` (30,547,659,019 B, md5: `f576bcf5eb28366aeaecff01c50fff34`)
  - Intermediate `.tar.gz` archives cleaned from scratch; `download_log.md` appended.
- **T3 partial + T5/T6 (2026-08-17)**:
  - 9 `dataset_check_<Name>.qmd` notebooks, `onboarding_utils.py` with `subset_by_samples()` (T3.1), `onboarding_metrics.R` (T3.0).
  - T3.0 `_debug` validation passed; T6 annotation safety guards implemented in `process_chunk.R`.

**Pending**:
- T2 & T3 execution: Run count sanity check and render onboarding check notebooks locally on Mac across the 9 datasets from the mounted NAS.
- T4: Study summaries verification.
- T7: User usage decisions (benchmark / batch-effect / exclude).
- T8: `datasets.json` registration (pending user approval).
- T9: Diabetes (mouse) pipeline support (conditional).
- T10: Pipeline rollout on HPC.
This plan is NOT complete — it must stay at `.kilo/plans/` (do NOT archive).

## Context & goal

Reviewer asked for additional datasets. The chosen source study is PILOT-GM-VAE
(Joodaki et al. 2025, *Brief Bioinform* 26(5):bbaf547, DOI 10.1093/bib/bbaf547,
PMID 41097818 — matches the NAS folder name `JooM_2025_41097818`), whose authors
provide preprocessed, annotated AnnData objects on Zenodo + CellxGene. Full
feasibility table (colors, comments): `new_datasets_to_implement.md` (Excel
source `/Users/christianhalter/Desktop/ECODA_PAPER_DATASETS.xlsx`); the paper's
own dataset description: Zenodo record 15575593 (`datasets.pdf`, saved during
planning).

**Scope**: download 9 of the 12 paper datasets to the NAS, run per-dataset
checks (raw-count sanity, metadata exploration, UMAP batch check), document
each study, then (after user review) register them in datasets.json and roll
out the pipeline. Excluded: Kidney cancer [29] (too few samples), Pancreas
PDAC [31], Myocardial infarction (1) [6] (duplicate of MI-2).

**Decisions already made with the user**:
- Downloads now go over the **HPC, into BeeGFS scratch FIRST, then rsync to the
  NAS** (user decision 2026-08-17, replacing the original "direct to NAS from
  the Mac" plan): the Mac NAS mount is unstable (disconnected mid-download;
  partial `SEAAD_Alzheimer.h5ad` ≈ 5–10 GB left on the NAS — ignorable,
  re-downloaded on HPC), the Mac has only ~108 GB free vs ~150–180 GB of
  files, and HPC downloads are faster AND leave the data on the HPC for later
  pipeline use. Downloads land in `${HPC_SCRATCH_DIR}/_downloads/` (= BeeGFS
  scratch via the `$HOME/scratch` symlink: 1.1 PB, **no per-user size quota**,
  not backed up, HDD) — intentionally NOT `$HOME` itself (SSD, 1 TB quota,
  backed up, and the repo/scratch already holds pipeline data). After all
  tasks verify, a login-node tail `rsync`s scratch → NAS
  (`${NAS_SC_DIR}/JooM_2025_41097818/output/`, where the onboarding notebooks
  read them from the Mac), keeps the `_downloads` copy for later use, and
  deletes the original `.tar.gz` archives.
- **Parallel HPC downloads** (user request): one SLURM **array job with 8
  download tasks** (one per `--only` key: alzheimer, breast, covid_lupus_tar,
  diabetes, kidney, lung_tar, myocardial, parkinson), concurrent tasks
  (configurable, e.g. 3–4) so Zenodo-S3 and CellxGene-S3 transfers overlap
  and per-server limits are better used. Each task is `curl -L -C -`
  resumable (Ctrl-C/re-queue safe; partial files persist in scratch and are
   resumed on re-run) + verified per key (Zenodo md5s; CellxGene files
   size-verified via HEAD content-length — no `.md5` sidecar exists, see table
   footnote + T1.1) + tar-selective extraction for the two
   archive entries, then archive deletion. **Compute-node internet is
   smoke-tested first** (`srun -p debug-cpu curl -sI https://zenodo.org` +
   a real CellxGene object URL — the original bucket-root probe falsely
   failed, fixed in T1.1);
  HPC forum evidence says compute nodes have internet (baobab incident
  resolved, 2024; login-node transfers remain the documented high-bandwidth
  path). If the smoke test fails → login-node fallback: same script entry
  run locally on the login node under `nice` + `--limit-rate` (documented
  transfer path; heavy login use is otherwise against AGENTS.md).
  The original Mac→NAS script stays as a fallback only (NAS mount stable).
- Check notebooks live in **`notebooks/dataset_onboarding/`** (one
  `dataset_check_<Name>.qmd` per dataset + `README.md` summary table).
- The repo **annotation pipeline (scGate/HiTME/scATOMIC) RUNS** for the new
  datasets, with new edge-case safety checks (no-crash when a method annotates
  0/too-few cells), a new datasets.json suitability flag, and Figure-3
  handling in `benchmark_analysis.rmd`.
- `batch_effect_analysis.rmd` funkyheatmap (separation by bio_col vs
  batch_col(s)) is **out of scope** — only TODO.md items get added.
- Diabetes (mouse) and Lung are downloaded but conditional: Diabetes needs a
  mouse-gene decision (T9); Lung is checked by UMAP first (known batch effects).
- **Metrics**: LISI-based separation instead of plain silhouette — the repo's
  min/max-bounded `calc_lisi` (`src/utils/scoring_metrics.R:168`, transform
  `(N-LISI)/(N-1)` → [0,1], 1 = perfect separation, 0 = perfect mixing), which
  already feeds the paper's `calc_sep_score` (`scoring_metrics.R:6-34`).
  Computed at the **cell level** (X_pca/UMAP are cell embeddings; no
  sample-level centroids needed) **per cell type** (author CT column, all CTs
  with ≥ min cells, subsampled ≤1–2k cells per CT for speed/balance) on the
  **unintegrated PCA embedding (top 30–50 PCs)** — metrics on UMAP coords
  distorts distances; NEVER on `X_pca_harmony`/integrated space — for the
  bio label and each batch candidate. Only question answered: "is there a
  batch effect?" — a cell type with high batch separation = batch effect
  signal; low bio separation per CT is expected and fine (cells can be
  identical across conditions while composition differs). scIntegrationMetrics (carmonalab, R/Seurat-based;
  Andreatta et al. 2024, Nat Commun s41467-024-45240-z) and scib (Python) are
  noted as richer options but **deferred**: scIntegrationMetrics is the
  recommended metric suite for the Phase-4 batch-effect analysis (needs a pixi
  setup/install addition → user approval); scib/pylisi are not installed and
  would add a pixi.toml dependency.
- Onboarding notebooks use the **pixi default Python env (scanpy 1.12.2)** and
  do NOT use the R notebook loader (`load_all_functions.R` / the downsized
  20-package `imports.R` from commit `30c4caa` — that loader serves only
  `benchmark_analysis.rmd` / `batch_effect_analysis.rmd`). The LISI helper is a
  standalone Rscript. Verified: `thisutils` 0.4.7 present in the default env
  (`pixi run -e default Rscript`), sklearn 1.9.0 present; `scib`/`lisi`/`pylisi`
  NOT installed (avoid pixi.toml changes unless approved).
- **Metric validation on `_debug` first**: the per-cell-type LISI helper is
  validated against the `_debug` dataset (Joanito 5-sample subset with bio
  `sample.origin` + batch, e.g. `seqtec`) before rolling out to the 9 new
  datasets. Preprocessed `_debug` views available on the NAS at
  `/Volumes/Shared/Projects/ECODA_paper/_debug/output/`
  (`JoaI_2022_35773407_debug_5samples_{benchmark,batch_effect}_analysis_ECODAprocessed.h5ad`).
- **Local interactive run on the 64 GB Mac — sample-first, RAM-bounded
  subsetting (user requirement, added 2026-08-17)**: all onboarding work runs
  on the user's macOS laptop (64 GB RAM) with **minimal computation and fast,
  interactive turnaround**. No full-matrix loads ever: h5ad files are opened
  `backed='r'`; obs is read first (cheap even for 2M cells); then **10–20
  samples are selected stratified by bio condition and batch candidates**
  (≥1 per bio group; batch values covered when possible), their cells are
  read sample-by-sample (`.to_memory()` per sample → concat — bounded peak
  memory regardless of file size, works for the 29 GB Breast file), and caps
  are applied (configurable; exact numbers not critical): per-sample
  ≤2000 cells (range 500–5000), optional per-CT ≤100 cells to (stratified by
  batch to avoid artificial imbalance), overall target **~10,000 cells
  (max 50,000)**. Count sanity check stays full-file but value-sampled
  (200k values, sub-second); UMAP/PCA + per-CT metrics run on the in-memory
  subset (10–50k cells → seconds). Notebook target: <~2–3 min wall time per
  dataset, ~1–2 GB peak RSS. The full files remain on NAS — subsetting
  happens at read time, no subset downloads. Subsetting is onboarding-
  diagnostics only: the real pipeline (HPC) still uses the full data.

## Dataset inventory & sources (9 + 0 exclusions)

Author-provided files preferred (they match the paper's Table 1 numbers).
Alzheimer + Parkinson are NOT on Zenodo — the paper obtained them from
CellxGene (must be downloaded from there).

| # | Dataset (short key) | Study (PMID) | Paper: cells / samples / CTs | Author-provided file | Source | Size (GB) | md5 |
|---|---|---|---|---|---|---|---|
| 1 | Alzheimer (SEA-AD) | Gabitto 2024 Nat Neurosci (42486312) | 1,395,601 / 83 / 18 | `SEAAD_Alzheimer.h5ad` | CellxGene collection `1ca90a2d-2943-483d-b678-b809bf464c30`, dataset_version_id `c2b49431-9288-4d94-8ca5-f6723b72217e` (count matches paper exactly) | **49.5** (observed 2026-08-17; HEAD content-length 53,187,995,660 B) | size-verified\* |
| 2 | Breast cancer | Kumar 2023 Nature (37380767) | 714,331 / 126 / 10 | `BreastCncr_processed.h5ad` | zenodo 14615923 | 28.9 | `8b28a349c2c3638ddbfb3946a32d12ba` |
| 3 | Covid-19 PBMC | Ren 2021 Cell (34767776) | 993,171 / 151 / 10 | `Covid19_Ren2021.h5ad` (inside `Datasets.tar.gz`) | zenodo 8370081 | 36.35 (tar; selective-extract) | `d105b52dbba38ac49c2ffe8b3cf34e24` |
| 4 | Diabetes | Hrovatin 2023 Nat Metab (37697055) — **mouse** | 264,235 / 52 / 13 | `diabetes.h5ad` | zenodo 8370081 | 4.1 | `38189a381bad630fa39ce2d7ad3a0855` |
| 5 | Kidney (KPMP) | Lake 2023 Nature (41648348) | 104,314 / 45 / 14 | `Kidney_KPMP.h5ad` | zenodo 14615923 | 2.75 | `36ceb02ba23c559f80625ec7bef6884f` |
| 6 | Lupus PBMC | Perez 2022 Science (42115607) | 1,263,676 / 261 / 11 | `Lupus_Perez2022.h5ad` (inside `Datasets.tar.gz`) | zenodo 8370081 | (in tar) | (in tar) |
| 7 | Lung | Sikkema 2023 Nat Med (42362693) | 941,504 / 165 / 12 | `lungatlas.h5ad` (inside `lungatlas.h5ad.tar.gz`) | zenodo 7957118 | 17.2 (tar) | `0d0c97924f1b7a405b6ec3b55da02882` |
| 8 | Myocardial infarction (2) | Kuppe 2022 Nature (41937210) | 132,888 / 23 / 11 | `Myocardial_Infarc_2.h5ad` | zenodo 14615923 | 3.6 | `7431ae99250c99f11bf63e3034798af4` |
| 9 | Parkinson | Prashant 2024 Sci Data (39580497) | 2,096,155 / 97 / 11 | `Parkinson.h5ad` | CellxGene collection `d5d0df8f-4eee-49d8-a221-a288f50a1590`, dataset_version_id `0270e5e5-ce1d-4165-828e-699210189a92` (count matches paper exactly) | ~40–60 (estimate; observed SEAAD 49.5) | size-verified\* |

\* **CellxGene files CANNOT be md5-verified** (discovered 2026-08-17): the
`.h5ad.md5` sidecar URL returns 403 for everyone, the S3 `ETag` is a
multipart digest (e.g. `"72c502238a7281c1517f8b7218637949-6341"` — not the
md5), and the curation API
(`https://api.cellxgene.cziscience.com/curation/v1/collections/<id>`) exposes
only `filesize`/`filetype`/`url` per asset. Verification = final file **SIZE**
checked against a HEAD `content-length` (the object itself is publicly
reachable: HEAD 200) + computed md5 recorded as informational. Implemented in
T1.1.

Excluded (NOT downloaded/registered): Kidney cancer (GEO GSE242299), PDAC
(GSA CRA001160), MI(1) (inside the same tar — leave unextracted or move to
`_excluded/`, not registered).

Fallback sources if the author file lacks raw counts (count check FAIL):
GEO GSE158055 (Covid), GSE174188 (Lupus), GSE195665 (Breast), GSE211799
(Diabetes); CellxGene collections `4195ab4c-…` (Breast), `bcb61471-…` (KPMP),
`8191c283-…` (MI); Zenodo 7227571 (Lung original); CellxGene for Alzheimer /
Parkinson (same collections).

## Tasks

### T1 — HPC downloads into BeeGFS scratch + NAS sync (user executes; agent prepares scripts)
- **New scripts** (kept next to the existing Mac→NAS `download_datasets.sh`, which
  stays as a NAS-stable fallback):
  - `download_datasets_hpc.sh` — login-node **submitter**: sources
    `src/slurm_config.sh` (`source "${PROJECT_ROOT}/src/slurm_config.sh"`),
    `cd "${PROJECT_ROOT}"`, runs the **compute-node connectivity smoke test**
    first (`srun -p debug-cpu --time=00:05:00 curl -sI https://zenodo.org`
    and the CellxGene host; fail-closed → shortcut to login-node mode), then
    submits **one SLURM array job, 8 tasks** (`--array=1-8`, concurrency
    `--array` throttling ≈ 3–4, `--partition="${SLURM_PARTITION_CPU}"` i.e.
    `shared-cpu`, 1 cpu / 2–4 G per task, generous `--time` (e.g. 10–12h,
    resumable anyway)), waits for completion, gates every task state via
    `sacct` (like the repo's submitter tails), then the **login-node NAS sync
    tail**: `rsync` `${HPC_SCRATCH_DIR}/_downloads/` →
    `${NAS_SC_DIR}/JooM_2025_41097818/output/` (canonical file names stay
    identical), keeps the `_downloads` copy AND deletes the tar.gz archives
    in scratch; appends the per-task status + md5s + `--sync-only <job-id>`
    resume to `notebooks/dataset_onboarding/download_log.md` (repo log on the
    login node; user pulls/commits from the Mac).
  - `run_download_worker.sh` — **sbatch worker** (must not resolve
    `slurm_config.sh` from `BASH_SOURCE` — Slurm copies scripts to the spool
    dir; recover `SCRIPT_DIR` via `scontrol show job <SLURM_JOB_ID>`
    `Command=` field with the `BASH_SOURCE` fallback, per AGENTS.md):
    one task per `--only` key, `curl -L -C - --retry 5 --retry-delay 10
    --fail` into `${HPC_SCRATCH_DIR}/_downloads/<file>`, verified per key
    (Zenodo md5s from the table; CellxGene size-verified — the `.h5ad.md5`
    sidecar does NOT exist, see T1.1), tar tasks do **selective extraction**
    then delete the archive.
    Partial files persist in scratch; killed/requeued tasks **resume** on
    re-submission (idempotent, like the repo's other workers).
- **Parallelism**: the array achieves server overlap (Zenodo-S3 vs
  CellxGene-S3) and utilization of the HPC uplink; individual file
  downloads must NOT be run twice (one task per key).
- **Login-node fallback** (only if the smoke test fails): run the same worker
  body locally on the login node under `nice -n 19` + `curl --limit-rate`
  (documented high-bandwidth transfer path; otherwise login-node use is
  against AGENTS.md).
- **Storage**: `_downloads` lives on `$HOME/scratch` (→ `/srv/beegfs/scratch/`,
  1.1 PB, NO per-user size quota, not backed up, HDD) — never `$HOME` itself
  (SSD, 1 TB quota incl. the existing repo + scratch pipeline data; the raw
  download volume ~150–180 GB would risk filling it). Scratch is
  re-downloadable, so a retention purge is harmless; sync to NAS promptly ran
  by the tail.
- **CellxGene resolution**: SEA-AD `dataset_version_id c2b49431-…` and
  Parkinson `0270e5e5-…` already resolved (counts match the paper exactly);
  verify cell counts again after download + document any delta vs the paper's
  Table 1 (SEA-AD collection was updated 2026-06-10).
- **Log entry per file** (in `download_log.md`, also committed): source URL,
  date, size, md5 expected + actual, extraction status — feeds the AGENTS.md
  "when and where it was downloaded" requirement.
- **Prerequisite**: the HPC clone must be current — user runs
  `cd "${HOME}/ECODA_paper" && git pull` (login node, quick) after the new
  scripts are implemented+committed, before first use.
- **Cancel/resume**: cancellation is safe at any point (Ctrl-C / `scancel`);
  re-running the submitter (`--only <key>` or all, or `--sync-only <id>` for
  the tail) resumes partial files via `curl -C -`.

### T1.1 — First-run bugfixes (2026-08-17, detected on the first HPC run) — DO BEFORE RE-RUNNING

Two bugs were exposed by the user's first run and verified from the Mac
(HEAD requests + curation API). Implement these fixes in
`notebooks/dataset_onboarding/`, commit + push, then the user re-runs.

- **Bug 1 — egress smoke test probed the wrong URL (false login-node
  fallback).** `download_datasets_hpc.sh` smoke-tested
  `curl -sI https://datasets.cellxgene.cziscience.com` (bucket root); S3
  answers 403 to a bucket-root HEAD **always**, so the test failed even with
  working egress and the submitter fell into the slow `--limit-rate 2m`
  login-node mode. Verified: HEAD on the real SEA-AD object → 200
  (`content-length: 53187995660`). **Fix**: probe real object URLs —
  `https://zenodo.org` and the SEA-AD `.h5ad` object — with
  `curl -sIL -o /dev/null -w '%{http_code}'`, accept 2xx/3xx, and do NOT
  swallow `srun` stderr (capture `2>&1`, print the failing URL + output on
  failure) so a genuine `srun`/partition problem is visible.
- **Bug 2 — the CellxGene `.h5ad.md5` sidecar does not exist.** The
  planning discovery was wrong: the sidecar URL returns 403 for everyone,
  the S3 ETag is a multipart digest (not the md5), and the curation API
  exposes only `filesize`/`filetype`/`url` (no checksum). The alzheimer task
  failed here (`curl: (22) ... 403`). **Fix** (size-based verification):
  - `run_download_worker.sh` `fetch_and_verify()`: for `verify` keys,
    resolve the expected size via `curl -sIL <url>` → last
    `content-length` (regex-guard `^[0-9]+$`, else FAIL
    `ERROR=no-content-length`); the "already downloaded" skip path compares
    SIZE, not md5; after the download compare the final size
    (`stat -c %s`, mismatch → FAIL `ERROR=size-mismatch`); compute + record
    the md5 as informational (`MD5_RECORDED=<md5>` + `VERIFY=size-<bytes>`
    status kvs). Update the header comment (`"verify"` = size-verified, no
    sidecar).
  - `download_datasets.sh` (Mac fallback): same size-based logic with a
    portable size command (`wc -c < file`, macOS has no `stat -c`); keep
    Zenodo md5 checks unchanged; `log_entry` statuses `OK (cached)` /
    `FAIL (size)`.
  - `download_sources.sh`: update the `SRC_MD5="verify"` comment to the
    size-verification semantics; drop the now-unused `SRC_SIDECAR` lines.
  - `download_log.md` report: the submitter's `append_log` picks up the new
    `VERIFY`/`MD5_RECORDED` status keys (extend the note grep
    `^(SKIPPED|FILES|ERROR)=` → include `VERIFY`), so CellxGene entries show
    "size-verified: N bytes; md5 (informational): …".
- **Docs** (same commit): AGENTS.md onboarding step 1 and
  `notebooks/dataset_onboarding/README.md` — replace all `.h5ad.md5`
  sidecar wording with the size-verification description (this plan's table
  footnote is already updated).
- **Validation**: `bash -n` on the 3 scripts; from the Mac, re-verify
  HEAD-content-length parsing on the SEA-AD object; on the HPC re-run, the
  smoke test must print `OK: compute nodes reach …` for BOTH URLs → array
  mode (if it still fails, the printed `srun` output decides whether egress
  is genuinely missing — then login-node mode is the correct fallback);
  alzheimer/parkinson tasks must log `size verified: <bytes>; md5
  (informational): …` and write `STATUS=OK` status files.
- **User re-run** (after `git pull` on the HPC clone): Ctrl-C of the aborted
  run is safe; `./notebooks/dataset_onboarding/download_datasets_hpc.sh`
  resumes the partial `BreastCncr_processed.h5ad` via `curl -C -`.

### T2 — Count sanity check (immediately after each download; agent runs)
- Python (`.pixi/envs/default/bin/python`, scanpy 1.12.2 available; backed read
  `sc.read_h5ad(backed='r')` for big files): locate counts
  (`layers['counts']` → `X` → `raw.X`), report slot used; verify sparse matrix
  and **integer VALUES** — do NOT require integer dtype: float-encoded CSR
  (1.0/2.0…) is valid raw counts; check via `np.all(np.mod(sample, 1) == 0)`
  (or `np.allclose(v, np.round(v))`) on a value sample (e.g. 100k cells) with
  epsilon tolerance; **non-negative**; no NaN/Inf; also flag if X looks
  log-normalized (max >> integer range / non-integer fractions in X while
  counts live elsewhere). Report shape, sparsity, n genes, integer-value
  verdict, slot used.
- Verdict PASS/FAIL logged in the onboarding `README.md` table + per-dataset
  qmd. On FAIL → fallback source from the table above (notify user first).

### T3 — Onboarding notebooks (`notebooks/dataset_onboarding/`)
- **T3.1 — Implement `subset_by_samples()` + `SUBSET_CONFIG` (NEW, added
  2026-08-17 — the key local-interactive requirement)**: add to
  `onboarding_utils.py` a RAM-bounded, sample-first subsetting routine and a
  config dict with these defaults (overridable per notebook in the header
  block; exact numbers NOT critical):
  - `MAX_SAMPLES = 15` (range 10–20), `N_PER_BIO = 3–5`
  - `MAX_CELLS_PER_SAMPLE = 2000` (range 500–5000)
  - `MAX_CELLS_PER_CT = None` (optional; e.g. 100 — applied after concat,
    stratified by batch to avoid artificial imbalance)
  - `CELLS_TARGET = 10_000`, `CELLS_MAX = 50_000`, `SEED = 0`
  Flow: (1) `sc.read_h5ad(path, backed='r')` + obs-only access; (2) select
  10–20 samples stratified by bio condition and batch candidates (≥1 per bio
  group; distinct batch values covered when possible — round-robin over the
  batch-candidate values within each bio group, else random with seed); small
  datasets (<10 samples) keep all samples; (3) read the selected samples'
  cells **per sample**: boolean mask → slice → `.to_memory()` (bounded peak
  memory regardless of file size; works for CSR- and CSC-on-disk, incl. the
  29 GB Breast file) → concat; (4) apply per-sample cap (random within sample,
  seed), optional per-CT cap, then overall target cap (stratified by sample so
  every selected sample stays represented); (5) return the in-memory AnnData
  + a summary dict (samples per bio group, cells per sample, total cells).
  Notebooks get a new **section 1.5** that calls it and prints the subset
  summary with the "diagnostic subset only — HPC pipeline uses full data"
  note. UMAP (sec. 4) and metrics (sec. 5) then run on this subset; target
  wall time per notebook < ~2–3 min, peak RSS ~1–2 GB on the 64 GB Mac.
- **T3.0 — Validate the metric + embedding approach on `_debug` first**
  (before writing the 9 notebooks): run `embed_and_umap` + `onboarding_metrics.R`
  (per-cell-type LISI separation on the **unintegrated `X_pca`** — never
  `X_pca_harmony` — bio label `sample.origin` vs batch `seqtec`; per-CT
  subsample caps + confounded-CT guards in effect) on
  `/Volumes/Shared/Projects/ECODA_paper/_debug/output/JoaI_2022_35773407_debug_5samples_batch_effect_analysis_ECODAprocessed.h5ad`.
  Checks: helper runs end-to-end, `thisutils`/`vegan`/`cluster`/`mclust`
  resolve when sourcing `scoring_metrics.R` standalone, per-CT table + heatmap
  produced, PCA input path exercised (incl. the precomputed-`X_pca` case), and
  the result is sensible (batch separation visible within CTs; bio/batch
  expectations from the known Joanito batch structure). Fix the helper until
  it passes, then replicate for the 9 new datasets. **Subsetting-path check
  (from T3.1)**: run `subset_by_samples` on the same `_debug` view with a
  reduced config (e.g. `MAX_SAMPLES=4`, `MAX_CELLS_PER_SAMPLE=2000`,
  `CELLS_TARGET=5000`) and confirm the batch signal (Site within CTs) is
  still detected on the subset, per-sample slice reads work, and the measured
  wall time + peak RSS stay in budget — this validates the exact code path
  the 9 notebooks will use on the Mac. Plots for `_debug` also serve as a
  sanity baseline for the UMAP + sample-bar sections.
- Shared helper `onboarding_utils.py` (backed-mode friendly): `locate_counts`,
  `count_sanity_check` (**integer-VALUE check with epsilon tolerance, not
  strict dtype — float-encoded CSR like 1.0/2.0 is valid raw counts**),
  `obs_summary`, `candidate_col_detection` (heuristics for sample /
  biological-label / batch / CT columns), `paper_table_compare`,
  `cells_per_sample_stats`, `embed_and_umap` (unintegrated scanpy recipe on
  RAW counts + precomputed-obsm detection),
  `subset_by_samples` (NEW — sample-first RAM-bounded subsetting, see T3.1),
  `confounding_crosstab` (bio × batch contingency table + collinearity
  warning), `save_png` — plus a standalone R helper `onboarding_metrics.R`
  (cell-level per-cell-type `calc_lisi` on unintegrated PCA, see section 5;
  **validated on `_debug` before rollout**, see Decisions). One
  `dataset_check_<Name>.qmd` per dataset (9 files; names: `Alzheimer`,
  `Breast_cancer`, `Covid19_PBMC`, `Diabetes`, `Kidney_KPMP`, `Lupus_PBMC`,
  `Lung`, `Myocardial_infarction`, `Parkinson`).
- Notebook sections (add a **subsetting section** right after section 1: the
  whole downstream analysis — sections 3–6 — operates on the in-memory
  sample subset from `subset_by_samples`):
  - **1.5. Sample-first subsetting (RAM-bounded, ~10k cells)**: backed
    `obs`-only read → candidate column confirmation → select **10–20 samples
    stratified by bio condition + batch candidates** (≥1 per bio group;
    batch values covered when possible; seed fixed) → read those samples'
    cells per-sample (cap ≤ per-sample cells, `.to_memory()` → concat) →
    optional per-CT cap (stratified by batch) → overall target cap
    (~10k, max 50k) → print the subset summary (n samples selected per bio
    group, n cells per sample, n cells total) with a note that this is a
    diagnostic subset only (HPC pipeline uses full data).
  0. **Study summary card** (agent-filled, from the original paper): title,
     authors, journal, PMID, DOI/link, download source + date + checksum;
     biological groups (with n per group); batch-candidate variables;
     notes (e.g. mouse — Diabetes; brain-region dominance — Parkinson).
  1. **File structure**: n cells × n genes, X / layers / raw summary, obsm
     keys (precomputed `X_umap`/`X_pca`?), obs columns + dtypes.
  2. **Count sanity check** (re-runs T2 with verdict).
  3. **Metadata exploration**: value counts per categorical obs column
     (bounded, top-20); **sorted bar plot of cells per sample**; n samples,
     total n cells, cells-per-sample min/median/mean±sd/max, n samples < 50
     cells; **n unique cell types per CT-annotation column**; **NaN check in
     CT columns** (precedent: Lee/Zhang HiTME NaNs broke PILOT);
     gene-symbol sanity (mouse vs human — Diabetes); **comparison vs
     PILOT-GM-VAE Table 1** (n cells / samples / CTs); **confounding /
     design-matrix check**: cross-tab `table(bio_label, batch_candidate)`
     per batch candidate — if a batch is (near-)perfectly collinear with the
     bio label (e.g. all Disease at Center 1, all Control at Center 2),
     print a loud warning: batch and bio are statistically indistinguishable
     in that dataset, and metrics/UMAP interpretation must account for it.
  4. **UMAP** (unintegrated space ONLY — never Harmony/integrated): computed
     on the **in-memory sample subset (~10–50k cells)** — use precomputed
     `X_umap` if present (sliced to the subset rows), else compute on the
     subset's RAW counts (normalize_total → log1p → HVG 2000 → PCA →
     neighbors → UMAP; seconds at this size). Panels: biological label,
     each batch candidate, CT low-res, CT high-res.
  5. **Quick separation metrics** (auxiliary, seconds): cell-level LISI-based
     separation via the repo's own functions — a standalone R helper
     `notebooks/dataset_onboarding/onboarding_metrics.R` sources
     `src/utils/scoring_metrics.R` directly (NOT the notebook loader — it is
     not sourced by these notebooks) and runs `calc_lisi()` on the
     **unintegrated PCA embedding (top 30–50 PCs)** of the sample subset —
     NOT UMAP coords (2D UMAP distorts distances; PCA is the standard input
     for neighborhood metrics per scIB/scIntegrationMetrics) and NEVER
     `X_pca_harmony`/integrated space. Use the global PCA (all subset cells
     in the same space); UMAP coords only as documented fallback when PCA is
     unavailable. For the whole subset AND for **each cell type** (author CT
     column; rand seed fixed): subsample up to N≈500–1000 cells per CT
     (balanced power across abundant and rare CTs; keeps kNN/memory trivial);
     separation score on (a) the biological label and (b) each batch
     candidate within the CT.
     Guards: CTs with < min_cells (e.g. 50) skipped (listed); CTs with
     >90% cells from a single batch or single bio group are NOT informative
     → tagged "confounded/uninformative" instead of a score (avoids
     false-positive batch warnings, e.g. rare CT present in one donor);
     single unique label → calc_lisi returns NA (existing behavior, note in
     output). Output: cell-type × label separation table (heatmap) + verdict:
     "batch effect detected in cell types X, Y" (high batch score within a
     CT = batch effect; interpret together with the confounding crosstab —
     collinear designs make bio vs batch indistinguishable). Low bio score
     within a CT is expected and NOT a problem (cells can be identical across
     conditions while composition differs — only "is there a batch effect?"
     matters here). Plots remain the primary deliverable for the user's
     visual inspection. The richer metric suite (CiLISI/norm_cLISI/
     celltype_ASW of carmonalab `scIntegrationMetrics`, and scib in Python)
     is documented as a Phase-4 option, not implemented here (would need pixi
     env additions → user approval).
  6. **Summary + recommendation** (agent): benchmark / batch-effect /
     negative control / exclude, with rationale incl. annotation-suitability
     note (e.g. brain dataset: no/minimal immune cells → HiTME/scATOMIC
     likely unsuitable).
- Plots → `data/new_dataset_checks/<ds>/*.png` (gitignored); optional NAS
  mirror to `JooM_2025_41097818/plots/`.
- `README.md` in the folder: summary table (dataset, study, PMID, file, size,
  checksum, count-check verdict, n cells/samples/CTs, batch candidates, UMAP
  verdict, recommended use, status) + links to `new_datasets_to_implement.md`
  and this plan.

### T4 — Study summaries
- Agent reads each original paper (abstract + methods) for T3 section 0 and
  the metadata-column descriptions; confirms candidate columns against the
  actual obs columns.

### T5 — Docs
- **AGENTS.md**: new section "Onboarding new datasets" — the procedure
  (download author-provided data → count sanity check → onboarding notebook +
  UMAP check → user review → datasets.json registration → pipeline rollout),
  what to document per dataset (study, PMID, groups, batch variables,
  metadata-column description, download source/date), pointers to
  `notebooks/dataset_onboarding/`, the `JooM_2025_41097818` NAS folder, and
  the appendix `new_datasets_to_implement.md`. Write this section against the
  **current** AGENTS.md wording (recent commit `30c4caa` downsized the R
  notebook loader to ~20 pkgs and added `src/utils/env_check.R`; the
  onboarding section must state that these Python check notebooks use the
  pixi default env and do NOT source the R loader).
- **TODO.md**: mark 5.1–5.4 progress; add items:
  - `batch_effect_analysis.rmd`: funkyheatmap like Figure 2a showing separation
    by bio_col and batch_col(s) (multiple batch columns possible) — out of
    scope for this plan;
  - Diabetes mouse-gene preprocessing support (T9);
  - annotation edge-case safety checks + `not_suitable_for_auto_annotation`
    flag + Figure 3 handling (T6);
  - per-dataset annotation-rate documentation (%-annotated cells per method).

### T6 — Annotation safety checks + suitability flag
- `src/4_cell_type_annotation/2.1.1_process_chunk.R`: after each annotation
  method (scGate, HiTME layer2/3, scATOMIC), guard on results — 0 annotated
  cells, <2 unique cell types, or all-NA → **warn + write the column as
  NA/unclassified (no crash)** + record per-method stats (n cells, n types).
  Keep the existing `cutoff.scATOMIC::em` patch and wall-time guards intact.
- `datasets.json` new **optional** field (proposal — flat, like existing
  fields): `"not_suitable_for_auto_annotation": ["hitme", "scatomic"]`
  (absent/empty = suitable for all methods). Consumers:
  - `benchmark_pipeline.R` HiTME/scATOMIC combo guards (currently
    `benchmark_pipeline.R:1398-1416`; scATOMIC combos likewise) — skip the
    flagged methods per dataset;
  - `benchmark_analysis.rmd` Figure 3 (ECODA_automated bars): exclude flagged
    datasets from the HiTME/scATOMIC bars with a documented note (dataset
    stays in Figure 2a where author CTs drive the methods; Leiden-based
    rescue possible later, see below).
- Annotation-rate documentation: after the first HPC annotation run of the
  new datasets, compute %-annotated cells + n unique types per method from the
  merged view h5ads → `notebooks/dataset_onboarding/annotation_summary.json` +
  update the per-dataset qmd section 6 with a short rationale paragraph on why
  HiTME/scATOMIC may be inappropriate (e.g. no/minimal immune cells).
- Fallback for datasets with unusable annotations (future, only if needed):
  unsupervised Leiden-derived CT column via a dataset-specific step like
  `1.4.1_create_kfoury_lowres_ct.R`.

### T7 — Usage-decision checkpoint (user)
- User reviews the UMAP plots + summaries; decides per dataset: benchmark /
  batch-effect / negative control (Parkinson?) / exclude (Lung?, Diabetes?).
  Expected mapping (per the feasibility table): benchmark → Alzheimer, Lupus,
  MI(2), Kidney; batch-effect → Breast, Covid-19, Diabetes (conditional),
  Lung (conditional); Parkinson → decide after check.

### T8 — datasets.json registration (ASK THE USER FIRST — AGENTS.md rule)
- Entries for the approved datasets: `folder_name: "JooM_2025_41097818"`,
  `file_names` (as stored on NAS), `columns` (sample / label / batch / CT
  columns confirmed in T3), `meta_cols_keep`, `views` with `subset_vars`
  (e.g. Diabetes minus embryos E12–E15; Covid-19 PBMC-only if not already
  subset; Kidney already scRNA-only in the author file), `use_for_benchmark`
  / `use_for_batch_effect` per T7, `tissue` / `normal_tissue`, and the
  annotation flag per T6. MI(1), PDAC, Kidney cancer NOT registered.
- Note for views: benchmark views use the **author CT columns** (same as
  existing datasets); batch-effect views use `columns.batch` (one batch col
  per entry; multi-batch handling belongs to the out-of-scope
  `batch_effect_analysis.rmd` funkyheatmap task).

### T9 — Conditional: Diabetes (mouse) pipeline support
- Only if Diabetes passes T3 and the user wants it: add mouse handling to
  `standardize_gene_symbols` (`src/utils/py/gene_utils.py`) + `1.1.1_preprocess.py`
  — detect mouse symbols, map to human orthologs (verify STACAS mouse gene
  conversion availability, e.g. `STACAS::convert_mouse_genes`; else MGI/biomaRt
  orthology table; pixi.toml untouched unless a package is really needed),
  adapt mitochondrial-gene detection (`mt-` prefix). **Requires user sign-off**
  before touching the pipeline.

### T10 — Pipeline rollout (HPC; user runs; agent provides commands)
- Validation first with the two smallest green datasets (Kidney ~2.75 GB /
  104k cells, MI(2) ~3.6 GB / 133k cells), end-to-end:
  `./src/1_stage_data/1_stage_data.sh --ds_name <ds>` →
  `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name <ds>` →
  `./src/4_cell_type_annotation/1_prepare_chunks.sh production` (per-dataset
  chunking; then `2_submit_hpc_array.sh` + `3_submit_merge.sh`) → benchmark
  submitters (`run_r_sample_embedding_methods`, `run_python_sample_embedding_methods`,
  `run_transformation_zeroimp_analysis`) / batch-effect submitters as
  applicable. Verify NAS outputs, then roll out to the remaining datasets.
- Big datasets (Breast 28.9 GB file, Alzheimer/Parkinson 1.4–2.1 M cells):
  preprocess tasks need high-mem nodes — check/raise `--mem` in the
  preprocess submitter before rollout.

## Risks & fallbacks
- **Recent-commit impact (checked 2026-08-16)**: commit `30c4caa` (plan
  `1786902582051-downsize-notebook-imports-loader.md`) downsized the R
  notebook loader (`imports.R` → 20 pkgs) and added `src/utils/env_check.R`.
  No structural impact on this plan: the onboarding notebooks are Python
  (pixi default env) and do not source the R loader; the LISI helper is a
  standalone Rscript sourcing `src/utils/scoring_metrics.R` directly (needs
  `thisutils`/`vegan`/`cluster`/`mclust` — verify at implementation; default
  env currently has `thisutils` 0.4.7). The AGENTS.md section (T5) is written
  against the current post-downsize wording.
- **SEA-AD collection updated** (2026-06-10): counts may differ from the
  paper; verify against Table 1 and document the delta; old dataset versions
  may be retrievable via `dataset_version_id` if a closer match is needed.
- **Author files may lack raw counts** → T2 FAIL → fallback source table
  (above); inform the user before substituting.
- **Huge files** (29 GB h5ad, 2 M cells): on the Mac, **never load full
  matrices** — backed reads + `subset_by_samples` (T3.1) bound peak memory to
  ~1–2 GB regardless of file size; per-sample `.to_memory()` slices keep the
  per-read allocation small even for the 29 GB BreastCncr file; count check
  is value-sampled; UMAP never exceeds the ~10–50k-cell subset. High-mem
  nodes on HPC for the real pipeline.
- **Annotation failure on new tissues** (brain/heart/kidney/pancreas): safety
  checks + flag (T6); tissue values in datasets.json must map to the worker's
  supported tissue dispatch — validate on Kidney/MI(2) first (T10).
- **NaN in author CT columns** (Lee/Zhang precedent): checked in T3, guarded
  downstream.
- **Diabetes mouse genes** break `standardize_gene_symbols` (human Ensembl105
  map) → T9 conditional; until then Diabetes stays unregistered/benchmark-
  disabled.
- **Disk space**: sequential downloads, selective tar extraction, delete tars
  after verification.
- **HPC storage quotas/limits**: downloads go to `$HOME/scratch` (BeeGFS
  scratch, no size quota, 10M-file limit, not backed up, re-downloadable) —
  keep them out of `$HOME` itself (SSD, 1 TB quota already loaded with the
  repo + scratch pipeline data). Do NOT run the download tasks concurrently in
  duplicate (one task per key) and confirm free space beforehand
  (`df -h $HOME/scratch`); if the scratch retention policy purges idle files,
  re-download (idempotent).
- **Compute-node internet**: evidence OK (forum 2024 incident resolved) but
  verify with the `debug-cpu` smoke test; if nodes have no egress, fall back
  to quiet login-node downloads (`nice`, `--limit-rate`) — the documented
  high-bandwidth transfer path. First run (2026-08-17) falsely failed the
  smoke test: it HEAD-ed the CellxGene bucket ROOT, which 403s always —
  fixed in T1.1 to probe real object URLs (whether compute nodes genuinely
  have egress is still unverified until the fixed smoke test runs).
- **CellxGene files have no md5** (discovered 2026-08-17): the `.h5ad.md5`
  sidecar 403s for everyone, the S3 ETag is a multipart digest, the curation
  API has no checksum. Verification is size-based (HEAD content-length) +
  informational md5 — a size match plus the onboarding count-sanity check
  (T2, which verifies cell counts vs the paper) together guard integrity.
  If the size check ever fails, re-download (idempotent, `curl -C -`).
- **Mac NAS mount instability** (observed 2026-08-17): the original
  Mac→NAS serial downloads are replaced by the HPC route; the Mac script is
  kept only as a fallback. A partial `SEAAD_Alzheimer.h5ad` (~5–10 GB) may
  remain on the NAS — harmless, overwritten at sync.

## Validation
- Downloads verified per key (Zenodo: md5 against the table values; CellxGene
  Alzheimer/Parkinson: final SIZE vs HEAD content-length + informational md5
  — no md5 sidecar exists, see T1.1) — logged per file in `download_log.md`.
- Compute-node connectivity smoke test (debug-cpu) passes with the FIXED
  URLs (zenodo.org + real CellxGene object) → compute-node downloads; else
  quiet login-node fallback engaged (T1.1).
- All 8 array tasks `COMPLETED` (sacct-gated), NAS sync verified:
  `ls -lh ${NAS_SC_DIR}/JooM_2025_41097818/output/` from the Mac matches the
  `_downloads` listing (Zenodo files: same md5s; CellxGene files: NAS md5 ==
  worker-recorded `MD5_RECORDED`), tar archives gone.
- Count sanity check PASS/FAIL logged per dataset.
- Notebooks render (quarto) with plots + summary tables; `README.md` table
  complete.
- Paper-table comparison (cells/samples/CTs) within tolerance.
- **Metric + embedding helpers validated end-to-end on `_debug` first** (T3.0;
  NAS `_debug` batch-effect view; bio `sample.origin` vs batch `seqtec`) —
  helper passes before any of the 9 datasets is processed; per-dataset
  cell-type × label separation tables produced via
  `pixi run -e default Rscript` `onboarding_metrics.R` (thisutils 0.4.7
  verified present).
- **Sample-first subsetting (T3.1) validated on `_debug` + the first real
  download** (Kidney or MI(2) once on NAS): subset summary printed, batch
  signal still detected at reduced sample counts, wall time < ~2–3 min and
  peak RSS ~1–2 GB on the Mac; `SUBSET_CONFIG` tuned if out of budget.
- User visually reviews UMAPs → usage decisions (T7) → datasets.json approved
  (T8).
- HPC: Kidney + MI(2) end-to-end validation before full rollout (T10).

## Out of scope
- `batch_effect_analysis.rmd` funkyheatmap (bio vs batch_col(s)) — TODO.md
  item only.
- PILOT-GM-VAE/QOT method rollout (separate plan; TODO 3.2).
- Phase 6 items (6.1–6.12), whole-Stephenson-by-center (4.2) extension.
- Excluded datasets (Kidney cancer, PDAC, MI(1)) — not downloaded/registered.
