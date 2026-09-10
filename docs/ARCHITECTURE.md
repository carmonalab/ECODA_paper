# ECODA Pipeline Architecture & System Design

This document serves as the single technical source of truth for the **ECODA** reproducible benchmark suite: pipeline stages, HPC execution layout, data flow, and module reference.

---

## 1. System Overview & Core Design Principles

The ECODA benchmarking pipeline is designed for **cohort-level exploratory analysis of scRNA-seq datasets**. It compares compositional representations against gene-expression embeddings and baselines across standardized cohorts in a strictly unsupervised setting.

```
+---------------------------------------------------------------------------------------------------+
|                                       RAW DATA (NAS Storage)                                      |
+---------------------------------------------------------------------------------------------------+
                                                  |
                                                  v (1_stage_data.sh)
+---------------------------------------------------------------------------------------------------+
|                               STAGE 2: PREPROCESSING & ANNOTATION                                 |
|                                                                                                   |
|  [2_dataset_specific_preprocessing] -> [3_scrnaseq_preprocessing] -> [4_cell_type_annotation]       |
|    - GongSharma sample cap               - Standardized HVG, PCA       - Chunked scATOMIC & HiTME |
|    - CombinedPBMC combine                - Harmony batch correction    - Backed CSR union reads   |
|    - Joanito seqtec & _debug             - Leiden clusterings          - Atomic obs metadata join |
+---------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+---------------------------------------------------------------------------------------------------+
|                             STAGE 5: BENCHMARK & TRANSFORMATION ARRAYS                                |
|                                                                                                   |
|  [Python Benchmark Array]           [R Benchmark Array]            [Transformation & Zero-Imp]    |
|   - MrVI, scPoli, PILOT-GM-VAE       - GloScope, MOFA, Pseudobulk   - 7 CLR transformations       |
|   - PILOT (EMD), QOT                 - scITD, Composition (ECODA)   - 6 Zero-imputation methods   |
|   -> Emits .feather dists/embs       -> Emits .rds result bundles   -> Emits .rds bundles         |
+---------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+---------------------------------------------------------------------------------------------------+
|                             EVALUATION & PUBLICATION NOTEBOOKS                                    |
|                                                                                                   |
|  notebooks/benchmark_analysis.rmd                      notebooks/batch_effect_analysis.rmd        |
|  - Ingests precomputed .rds & .feather results          - Evaluates bio vs batch separation       |
|  - Computes ANOSIM, ARI, Silhouette, Modularity, LISI   - Tests batch-correction robustness       |
|  - Renders Figure 2A, Supp Fig 2, Supp Fig 14-21       - Generates batch evaluation panels        |
+---------------------------------------------------------------------------------------------------+
```

### Key Technical Contracts
1. **Backed HDF5 / CSR on Disk:** Preprocessed `.h5ad` files retain raw counts in `layers["counts"]` and use CSR storage. Annotation unions deliberately contain only raw counts in `X`, so their backed chunk workers can slice samples selectively. Login-side source-ID and annotation-contract validators use h5py-only metadata/`obs` reads; they must not use `anndata.read_h5ad(..., backed="r")`, because anndata 0.12.x may eagerly materialize `layers`.
2. **Run-Owned Minimal Annotation Union:** Each Stage 4 run writes a minimal union at `${HPC_SCRATCH_DIR}/_ecoda_runs/<RUN_ID>/datasets/<DS_NAME>/union/union.h5ad` (`X` = raw counts CSR, `obs` + `var` only, no `layers`/`obsp`/`obsm`). Feathers are validated once against the complete union `(Sample, cell_barcode)` key set, then projected onto each selected view.
3. **Cross-Language Data Exchange (`.feather`):** Python methods serialize embeddings and pairwise distance matrices to Apache Arrow `.feather` files, which R consumers read with zero format drift.
4. **Atomic Writes & MD5 Checksums:** All R result bundles, annotation feathers, unions, and merged h5ads are written atomically (temp file + rename) with MD5/SIZE/PATH sidecars before synchronization to persistent NAS storage.

---

## 2. Dataset Metadata & Schema (`datasets.json`)

All evaluated cohorts and views are configured in [`datasets.json`](../datasets.json).

### Schema Structure
```json
{
  "<DS_NAME>": {
    "display_name": "Human-readable label",
    "file_names": "Raw source file(s) in NAS storage",
    "folder_name": "Source subfolder in NAS_SC_DIR",
    "tissue": "Tissue of origin",
    "normal_tissue": true,
    "use_for_benchmark": true,
    "use_for_batch_effect": false,
    "columns": {
      "sample": "Metadata column holding sample/patient ID",
      "label": "Ground truth biological condition column",
      "batch": "Technical batch variable column (if applicable)",
      "cell_type_low_res": "Broad cell type annotation column",
      "cell_type_high_res": "High-granularity cell type annotation column"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "Input raw file",
        "output_file_name": "<DS>_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": { "column": { "values": ["exclude_val"], "op": "notin" } }
      },
      "batch_effect_uncorrected": {
        "input_file_name": "Input raw file",
        "output_file_name": "<DS>_batch_effect_uncorrected_ECODAprocessed.h5ad",
        "subset_vars": {}
      },
      "batch_effect_corrected": {
        "input_file_name": "Input raw file",
        "output_file_name": "<DS>_batch_effect_corrected_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  }
}
```

### The `_debug` Dataset
- Built by `src/2_dataset_specific_preprocessing/1.3.1_prepare_joanito.R` into `${HPC_SCRATCH_DIR}/_debug/data/`.
- A 5-sample subset (500 cells/sample) covering combinations of biological origin and technical batches.
- Used for fast end-to-end integration tests of stages 2, 3, and 4. Default execution scripts skip `_*` entries unless `--ds_name _debug` is explicitly specified.

---

## 3. HPC & Storage Directory Layout

Execution is split between compute scratch (fast local I/O on worker nodes) and the NAS archive (persistent storage accessible via login nodes).

```
# HPC Scratch ($HOME/scratch/ECODA_paper -> ${HPC_SCRATCH_DIR})
${HPC_SCRATCH_DIR}/
├── <DS_NAME>/data/               # Staged raw inputs per dataset (1_stage_data.sh)
├── <DS_NAME>/output/             # Preprocessed .h5ad, annotation chunks, and merged outputs
├── <DS_NAME>/annotation_union/   # legacy fixed-path location; new runs use _ecoda_runs/<RUN_ID>/datasets/<DS_NAME>/union/
├── benchmark/
│   ├── embeddings/               # Python & R method distance matrices (.feather)
│   ├── results/                  # Evaluated R result bundles (<ds>_<method>.rds)
│   ├── pseudobulks/              # Shared DESeq2 normalized pseudobulk matrices
│   ├── gloscope_dists/           # Cached GloScope distance matrices
│   └── checksums.md5             # MD5 verification sidecar
├── chunks_manifest_<pid>.txt     # Per-submission chunk manifests
├── _worker_retries/              # Worker transient retry counters
└── _benchmark_watchdog/          # SLURM watchdog status logs

# Run ownership / manifests (scratch only)
${HPC_SCRATCH_DIR}/_ecoda_runs/<RUN_ID>/
├── metadata
├── manifests/                 # immutable selected dataset/view/method rows
└── status/                    # atomic watchdog/aggregate/terminal records
${HPC_SCRATCH_DIR}/_ecoda_owners/<stage>/<artifact>/
└── owner                      # ACTIVE/OK/FAIL ownership state

# NAS Storage (carmona_smb / Shared Collections; login node only)
${NAS_SC_DIR}/                    # Standardized raw source datasets
${NAS_REF_DIR}/                   # Carmona Lab reference maps (HiTME)
${NAS_TARGET_DIR}/                # Synced project results
├── <DS_NAME>/output/             # Rsynced preprocessed and annotated .h5ad files
└── benchmark/                    # Rsynced embeddings, results, and execution times
```

### Centralized Environment Variables (`src/slurm_config.sh`)

| Variable | Default Path / Value | Description |
|---|---|---|
| `PROJECT_ROOT` | `$HOME/ECODA_paper` | Root directory of the git clone |
| `HPC_SCRATCH_DIR` | `$HOME/scratch/ECODA_paper` | High-speed scratch directory on compute cluster |
| `NAS_SC_DIR` | `/srv/.../Standardized_SingleCell_Datasets` | Cold storage path for raw datasets |
| `NAS_TARGET_DIR` | `/srv/.../Projects/ECODA_paper` | Target sync directory on persistent storage |
| `HOME_REF_DIR` | `$HOME/reference_atlases/sketched_200ct` | Local sketched reference maps for HiTME |
| `SCGATE_DB_PATH` | `${PROJECT_ROOT}/aux/scGateDB.rds` | Pre-downloaded scGate model database cache |
| `PYTHON_BIN` | `${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python` | Target Python binary on HPC |
| `PIXI_RSCRIPT` | `${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/Rscript --vanilla` | Direct target Rscript command on HPC; avoids the r-base activation hook's shared `R CMD javareconf` writes |
| `SAMPLE_COLNAME` | `"Sample"` | Standardized sample metadata column name |

---

## 4. Pipeline Stages & Execution

### Stage 1 — QC Filtering (`notebooks/QC_filtering/`)
Per-dataset R Markdown notebooks performing study-specific initial quality control (filtering low-quality cells, doublets, and uninformative genes) before dataset standardization.

---

### Stage 2 — Preprocessing & Cell Type Annotation

```
+------------------+      +--------------------------------+      +-------------------------------+      +--------------------------+
| 1_stage_data.sh  | ---> | 2_dataset_specific_preproc     | ---> | 3_scrnaseq_preprocessing      | ---> | 4_cell_type_annotation   |
| (NAS -> Scratch) |      | (GongSharma, Combined, Joanito)|      | (Standardized Scanpy array)   |      | (scATOMIC + HiTME array) |
+------------------+      +--------------------------------+      +-------------------------------+      +--------------------------+
```

#### 1. Data Staging (`src/1_stage_data/`)
- `1_stage_data.sh`: Login-node utility that queries `datasets.json` and rsyncs required raw input files from `${NAS_SC_DIR}` to `${HPC_SCRATCH_DIR}/<DS_NAME>/data/`. Supports `--ds_name <DS>`.

#### 2. Dataset-Specific Preprocessing (`src/2_dataset_specific_preprocessing/`)
- `1_submit_hpc.sh`: selected hook gate. Independent hooks are submitted in
  one wave; only the GongSharma cap -> CombinedPBMC read/write edge uses
  `afterok`. `stage2_watchdog.sh` owns terminal accounting, OOM-only retries,
  semantic prerequisite validation, and atomic checksums. Stage 2 remains
  scratch-only.
- Hook outputs are installed atomically before the next numbered stage reads
  them. The CombinedPBMC legacy raw basename is accepted only for guarded
  content/checksum validation and one-time rename to `combined_pbmc.h5ad`.

#### 3. Standardized Preprocessing (`src/3_scrnaseq_preprocessing/`)
- `1_submit_hpc_array.sh` -> `1.1_run_worker.sh` -> `1.1.1_preprocess.py`
  consumes one immutable `DATASET<TAB>VIEW` row per selected view. The exact
  batch mode requires the canonical twelve-row uncorrected selection and
  rejects legacy/corrected rows before scheduler submission.
- Standardizes gene names using Ensembl 105, applies view-specific subsets,
  preserves raw counts in `layers["counts"]`, normalizes/log-transforms,
  ranks HVGs, and computes PCA/Harmony/Leiden outputs.
- `1.2_preprocess_watchdog.sh` performs OOM-only reduced-row retries and full
  H5AD schema/checksum validation. Existing output candidates are validated by
  the compute-node H5AD preflight array before idempotent skip decisions, so
  the login submitter does not open large H5AD content for the full contract.
  Selected h5ads only are synchronized to NAS after the watchdog succeeds;
  Stage 3 transfers use one per-stage files-from manifest plus per-artifact
  remote checksum comparisons. `batch_effect_analysis` is not an accepted
  logical view; historical filename components are migration/documentation
  text only.

#### 4. Parallel Cell Type Annotation (`src/4_cell_type_annotation/`)
- `1_submit_onboarding_stage.sh` is the only production entrypoint. It stages
  and validates reference maps/scGate once, submits dataset-parallel
  preparation rows, one global annotation chunk array, and a dataset-parallel
  merge array (`3.2_merge_worker.sh`). View updates remain serial inside each
  dataset.
- Run-owned union/chunk/checkpoint/feather paths prevent concurrent runs from
  deleting one another. Watchdogs retry OOM rows only and validate union
  membership, Feather key uniqueness, full `(Sample, cell_barcode)` coverage, and
  atomic h5ad installation.
- The union is the immutable annotation key authority for a run. A target view
  may be a strict subset; merge validates the complete union once, then projects
  only matching keys into each target view. Reuse requires matching source
  h5ad path/MD5/SIZE records and canonical run-owned paths.
- `--skip-prepare` is only a reuse operation and requires
  `--reuse-run RUN_ID`; it never creates a new empty run root. The old
  preparation/annotation/merge submitters fail closed.
- Login-side synchronization runs once per selected dataset in bounded
  parallel after merge validation. Stage 4 batches each dataset's H5AD and
  sidecar transfer with an rsync files-from manifest, then compares every
  remote artifact checksum; failed sync preserves run artifacts.
All runnable datasets must produce the complete HiTME (`layer1`, `layer2`,
`layer3`) and scATOMIC (`layer_1` through `layer_6`, `scATOMIC_pred`,
confidence, and cell-cycle) schema. A sample with zero output is retained as
an all-NA checkpoint and reported in per-sample stats; the dataset-level
`layer1` and `scATOMIC_pred` anchors remain mandatory. The three
`not_suitable_for_auto_annotation` cohorts are skipped a priori and are not
failed annotation results. scATOMIC `breast_mode` is intentionally disabled:
all callers use the upstream default `FALSE` for cross-cohort comparability.

---

### Stage 5 — Benchmark & Method Analyses (`src/5_run_benchmark_methods/`)

`1_submit_hpc_array.sh` is the canonical coordinated wrapper. In ordinary
mode it accepts the benchmark matrix. The default batch mode requires
`--pass uncorrected` and the exact three-column twelve-row matrix selection.
The immutable batch matrix selection is
`DATASET<TAB>VIEW<TAB>SCOPE`; in pass mode `VIEW` and `SCOPE` must both equal
the selected `batch_effect_<pass>` view. `SCOPE` is not a method label; the
default batch methods are the fixed ordered suite
`prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot`.
For fail-closed recovery, an explicit `--target-methods LIST` may be combined
with `--pass`, `--selection-file`, and an explicit `--partition` to select only
the named batch methods for the supplied rows. It cannot be combined with the
exact twelve-row selection or ordinary analyses; the fixed suite remains the
default.

- `matrix_watchdog.sh`: compute-node OOM-only retry and terminal task gate for
  one method matrix. Retries bind both `MATRIX_RETRY_MANIFEST` and
  `ANALYSIS_MANIFEST` to the reduced run-owned manifest; batch retries never
  export `BENCHMARK_MANIFEST`.
- `matrix_gate.sh`: aggregate gate that requires every child watchdog to report
  `STATE=OK`.
- RDS validation requires the ordered source sample universe for every method
  except scITD. scITD may report an ordered subset and its dropped sample IDs;
  no other method receives this exception.
- Each Stage 5 run records a run-owned `manifests/source_identity.json` with
  source path, size, MD5, and ordered sample IDs. The identity sidecar is
  verified against current source contents at stage entry and before final
  validation; matrix/RDS validators consume the verified IDs.
- Source Sample IDs and annotation H5AD contracts use h5py-only `obs` reads.
  Full persisted-count H5AD validation is performed by the compute-node
  preflight array, avoiding anndata 0.12.x backed opens on the login node.
- GloScope, composition, PILOT, QOT, and PILOT-GM-VAE use the h5py-backed
  counts-free loader and receive only required obs columns plus stored
  embeddings. MrVI and scPoli require raw counts, but their loaders stream
  only the stored HVG columns into a minimal AnnData; scPoli densifies only
  that selected subset.
- Canonical Stage 5 pseudobulk uses only the persisted raw CSR counts in
  `layers["counts"]`; normalized/log-transformed `X`, labels, and Seurat
  aggregates are never a source of canonical counts. The Sample-only and
  Sample × cell-type paths below are separate contracts. MOFA consumes
  precomputed pseudobulks and only uses the bounded raw-count path when a
  required cache is missing.
- The canonical path passes raw aggregate matrices directly to the R
  matrix-to-DESeq2 boundary. It does not create a one-sample-per-column
  Seurat object or call `AggregateExpression()` a second time.
- Batch composition requires `ECODA_authors_HR`,
  `ECODA_authors_HR_NULL`, and `ECODA_seuratres_2`. Existing
  `ECODA_HiTME_HR_layer2` and `ECODA_scATOMIC_HR` bundles are recognized
  legacy extras but are not required or regenerated in batch mode; ordinary
  benchmark composition retains its annotation-specific outputs.
- Batch Feather skip checks use one-row dataset manifests, and fully populated
  per-dataset RDS skip checks are grouped into one R validator invocation.
- Pass roots, logs, watchdog status, manifests, and markers are scoped to
  `batch_effect/uncorrected` or `batch_effect/corrected`. Batch markers use the
  `BATCH_EFFECT_*` namespace; ordinary benchmark markers retain
  `BENCHMARK_*`.
- The family submitters are compatibility entrypoints that delegate to the
  canonical wrapper; they do not own independent synchronization.
- `batch_effect_analysis` is never a logical view or loader fallback.

#### Canonical pseudobulk dataflow

Stage 5 has one canonical count source and two deliberately separate
aggregation contracts:

- **Source and leakage boundary.** The only canonical pseudobulk source is
  the H5AD `layers["counts"]` CSR matrix. Values must be finite,
  nonnegative, integer-valued raw counts. Normalized/log-transformed `X`,
  `raw.X`, embeddings, and biological labels are not interchangeable count
  sources; labels are evaluation metadata and never DESeq2 covariates.
- **Sample-only contract.** Group by `Sample` only, retaining the first-seen
  sample order and first-observation metadata. The Python reducer returns
  checked-int64 genes-by-samples counts. The R boundary
  `get_pb_deseq2_from_counts()` consumes that matrix directly and publishes
  canonical samples-by-genes output with the same sample universe and
  metadata alignment. Ordinary pseudobulk uses `~1`, `blind=TRUE`, and no
  batch correction. Corrected batch-effect pseudobulk uses only the configured
  technical batch, `blind=FALSE`, and batch-only correction; biological
  labels remain outside the model.
- **Shared full-gene fit.** `fit_pseudobulk_deseq2()` performs the full-gene
  DESeq2 size-factor/normalization/VST fit once and returns the normalized/VST
  matrix plus its variance ordering. `select_pseudobulk_deseq2()` derives
  `hvg500`, `hvg1000`, `hvg2000`, and `hvg3000` from that fit, so shared work
  is not repeated per variant. `schvg2000` remains a separate restricted
  gene-universe fit. The documented `hvg2000_bl` no-op behavior is preserved.
  `validate_pseudobulk_counts_matrix()` enforces the R-side boundary before
  `DESeqDataSetFromMatrix`; `get_pb()` and `get_pb_deseq2()` are not part of
  this canonical path.
- **Sample × cell-type contract.** Group only present `(Sample, cell_type)`
  combinations; never materialize the Cartesian product. Cell types retain
  first-observation order (missing values excluded), and a group is eligible
  only when it has at least five raw cells. Each cell type is normalized
  independently; failed cell types are isolated, while the final distance
  universe is `sort(unique(Sample))` and successful per-cell-type pairwise
  distances are averaged with the existing denominator and counters. The
  canonical entry point is `process_pseudobulk_ct_h5ad_fig()`, which consumes
  the composite raw-count store rather than a full Seurat object.
- **Checked count transport.** Every persistent accumulator addition is
  checked in int64 before the addition; signed negative CSR values are
  rejected. Any sparse multiplication has a division-safe pre-bound for
  the selected limit. The canonical Stage 5 reducer requests the DESeq2
  ceiling `INT_MAX = 2,147,483,647`, and the R validator repeats finite,
  nonnegative, integer-valued, and `INT_MAX` checks before reticulate data
  reaches DESeq2. Overflow or an unsafe value fails closed; it is never
  wrapped, saturated, or silently coerced.

#### Run-owned cell-type stores and timing

`prepare_h5ad_ct_group_store()` performs a metadata pass followed by one raw
CSR count pass. It writes an exclusive, append-only/chunked HDF5/CSR store
under a unique run-owned path such as
`<run-owned scratch>/pseudobulk_ct/<run_id>/<unique-token>`. The store keeps
contribution `group_ids`, CSR `indptr`, `indices`, `data`, group cell counts,
and first-observation group metadata; it must remain sparse and bounded rather
than retaining a dense all-group array or a persistent
`Sample × cell_type × gene` dictionary. `read_h5ad_ct_group_store()` selects
only requested contribution rows without densifying all groups.

Each store has an ownership manifest containing run ID, PID/scheduler
identity, source H5AD identity/checksum, stage, schema, and creation time.
The canonical CT routine closes handles and cleans its unique store on normal
and ordinary error paths. Before reuse it invokes
`audit_h5ad_ct_group_store()` in no-compute mode: cleanup is allowed only for
a valid matching manifest, a demonstrably dead owner, and an expired
age/lock policy. Active, ambiguous, malformed, or mismatched stores fail
closed and are never removed.

New pseudobulk cache records use timing schema 2:
`$pb`, `$time_secs`, `$mem_GB`, `aggregate_time_secs`,
`shared_fit_time_secs`, `shared_time_secs`, `variant_time_secs`,
`shared_mem_GB`, `timing_id`, and `timing_schema`. Here `$time_secs` and
`$mem_GB` are variant-local; `shared_time_secs` is exactly
`aggregate_time_secs + shared_fit_time_secs`; and
`timing_id` is
`${ECODA_RUN_ID}:${cache_stem}:${view}:${analysis_pass_or_none}`.
The shared aggregate/fit row is logged once per timing ID, never once per
variant. Records without schema 2 remain readable with their old inclusive
timing interpretation. `shared_mem_GB` is a shared-stage peak or upper
bound: RSS/VmHWM is process-cumulative and must not be described as an
isolated component allocation.

#### Backend selection policy

The custom checked HDF5/CSR reducer and run-owned store are the production
baseline. Scanpy and decoupler remain optional candidates only. Neither is
adopted, nor added as an implicit dependency, unless a controlled local
contract/performance benchmark demonstrates exact contract preservation and a
clear wall-time or peak-RSS improvement. Existing production artifacts and
their cache/checksum/producer/force/missing-only semantics are not
invalidated by candidate evaluation.

#### Legacy Seurat boundary and maintained-caller audit

The maintained-caller audit found no maintained canonical Stage 5 caller for
the legacy Seurat pseudobulk APIs. The disposition is:

| API | Observed callers and disposition |
|---|---|
| `get_pb()` | Called internally by legacy `get_pb_deseq2()` and by historical `notebooks/batch_effect_analysis_legacy.rmd` sections. No canonical Stage 5 caller remains; retain as an isolated legacy boundary. |
| `get_pb_deseq2()` | Former canonical calls in `prepare_pseudobulks_hpc()` and `process_pseudobulk_ct_fig()` migrate to the direct matrix APIs. The Seurat branch of deprecated `run_benchmark_analysis()` and the historical batch-effect notebook remain legacy-only callers; retain the wrapper and do not use it for canonical production pseudobulk. |
| `load_h5ad_pseudobulk_seurat()` | Former preparation and missing-cache fallback callsites in the Stage 5 workers migrate to the raw CSR reducer/direct matrix path. The focused adapter fixture may exercise this helper, but it is not a maintained production caller; no canonical path may construct a sample-level Seurat object solely for pseudobulk. |
| `run_benchmark_analysis()` | Explicitly deprecated and notebook-only. `notebooks/benchmark_analysis.rmd` loads HPC result bundles rather than invoking it; its retained Seurat branch is a compatibility boundary and is not a canonical Stage 5 dependency. |

This audit does not edit historical notebooks or delete the legacy
implementations. Cell-level Seurat loaders needed by methods with a genuine
cell-count contract remain distinct from the legacy sample-pseudobulk
adapter.

#### Benchmark Analysis Notebook (`notebooks/benchmark_analysis.rmd`)
- Loads precomputed `.rds` result bundles and `.feather` matrices via
  `load_hpc_benchmark_results()`.
- Does **not** load monolithic raw `.h5ad` files into memory, enabling rapid
  local knitting on macOS.
- Generates publication figures (Figure 2A, Supp Fig 2, Supp Figs 14–21) and
  summary rank heatmaps.
---

### Stage 4 — Batch Effect Analysis (`notebooks/batch_effect_analysis.rmd`)

The notebook retains its historical filename, but registry/submitter identifiers
are the explicit `batch_effect_uncorrected` and `batch_effect_corrected` views.
The uncorrected twelve-cohort pass is the evidence gate; no corrected processing
is launched before its scientific review. Batch formulas consume only configured
high-resolution cell types. Pseudobulk, Harmony, and native MrVI correction
settings are batch-only and never protect biological labels.

---

### Dataset Onboarding Architecture (`notebooks/dataset_onboarding/`)

Standardized protocol for onboarding new external cohorts:
1. **HPC Download Worker (`download_datasets_hpc.sh` / `run_download_worker.sh`):** Array-based resumable `curl` downloading directly into BeeGFS scratch, followed by MD5 verification and NAS sync.
2. **Sample-First Subsetting (`onboarding_utils.py`):** Fast diagnostic inspection reading ~10–20 stratified samples directly from backed `.h5ad` files into memory (<2 GB RAM) for immediate exploratory analysis.
3. **Count Integrity Check:** Identifies raw integer count layers (`layers["counts"]` vs `raw.X` vs `X`).
4. **Standalone LISI Scoring (`onboarding_metrics.R`):** Calculates cell-level biological and batch LISI on unintegrated PCA space.
5. **Onboarding Check Notebooks (`dataset_check_<Name>.qmd`):** Diagnostic Quarto notebooks for review before registering cohorts into `datasets.json`.

The batch-effect onboarding scope is the ordered twelve-dataset selection
maintained in `dataset_specs.py`: nine audit cohorts, Joanito, Stephenson, and
CombinedPBMC. Stage 4 records `not_suitable_for_auto_annotation` as an
a-priori exclusion for non-immune/unsupported cohorts; it is not an annotation
failure. The annotation contract requires all dual-method columns, exact
`(Sample, cell_barcode)` union coverage, checksums, and dataset-level
`layer1`/`scATOMIC_pred` anchors. The batch analysis consumes only configured
high-resolution cell-type columns. scATOMIC `breast_mode` remains at default
`FALSE` and must not be passed by callers.

---

## 5. File & Module Reference Tables

### Core Utilities (`src/utils/`)

| Script / Module | Primary Functions & Purpose |
|---|---|
| `src/utils/math_utils.R` | `clr_transform()`, `zero_imputation()`, ALR/ILR transformations |
| `src/utils/hvcs.R` | `select_hvcs()` — variance-based selection of Highly Variable Cell Types |
| `src/utils/pseudobulk.R` | `validate_pseudobulk_counts_matrix()`, `fit_pseudobulk_deseq2()`, `select_pseudobulk_deseq2()`, `get_pb_deseq2_from_counts()` and `DESeq2.normalize()` — direct raw-matrix normalization; `get_pb()`/`get_pb_deseq2()` remain isolated legacy Seurat boundaries |
| `src/utils/scoring_metrics.R` | `calc_sep_score()` (ANOSIM), `clust_eval()` (ARI), `calc_sil()`, `calc_lisi()`, `calc_modularity()` |
| `src/utils/datasets_io.R` | `read_datasets_json()`, `get_dataset_view_info()` (R parser) |
| `src/utils/py/datasets_io.py` | `read_datasets_json()` (Python parser matching R semantics) |
| `src/utils/py/h5ad_pseudobulk.py` | Checked-int64 Sample-only aggregation plus `prepare_h5ad_ct_group_store()`, `read_h5ad_ct_group_store()`, and `audit_h5ad_ct_group_store()` for bounded run-owned Sample × cell-type CSR stores |
| `src/5_run_benchmark_methods/benchmark_hpc_utils.R` | H5AD source/metadata validation, direct raw aggregate bridge, cache identity, and schema-2 timing integration |
| `src/5_run_benchmark_methods/benchmark_methods_r.R` | `process_pseudobulk_ct_h5ad_fig()` canonical raw-matrix CT processing; `process_pseudobulk_ct_fig()` retained only as a legacy Seurat boundary |
| `src/5_run_benchmark_methods/benchmark_pipeline.R` | Stage 5 orchestration; `run_benchmark_analysis()` is deprecated notebook-only compatibility code, not a canonical pseudobulk entry point |
| `src/utils/py/h5ad_counts_subset.py` | Stored-HVG raw-count loading for MrVI/scPoli without full-gene materialization |
| `src/utils/py/gene_utils.py` | `standardize_gene_symbols()` using Ensembl 105 reference dictionary |
| `src/utils/bash/worker_retry.sh` | Sourced by SLURM workers for automated self-requeue on transient I/O faults |
| `src/utils/bash/sync_status_email.sh` | Best-effort email notification utility with per-task duration reports |
| `src/utils/env_check.R` | Validates package attachment and namespace loads across Pixi environments |

---

## 6. Resilience & Fault-Tolerance Mechanisms
1. **Transient Fault Recovery:** Array workers source `worker_retry.sh`. On non-zero exit codes caused by transient BeeGFS cache misses, workers grep the task `.err` log for known signatures and self-requeue via `scontrol requeue` (capped at `WORKER_MAX_RETRIES=3`).
2. **Deterministic Locking:** Environment refreshes acquire the atomic directory lock `logs/env_refresh.lock`; interrupted locks fail closed and require independent verification before removal. This prevents concurrent Pixi/R writes to the shared library.
3. **Fail-Closed Verification:** All submitter sync tails verify terminal scheduler/accounting state, artifact schemas, and MD5/SIZE/PATH sidecars before initiating NAS synchronization. If any task or transfer fails, synchronization is aborted and an alert email is dispatched.
4. **Run-ID Recovery:** `--sync-only RUN_ID` validates the immutable run manifests, scheduler records, watchdog/aggregate state, fresh artifacts, and owners before completing an interrupted login tail; it never resubmits. Numeric/CSV scheduler-ID recovery remains a separate compatibility path and requires the caller's original dataset/view selection.