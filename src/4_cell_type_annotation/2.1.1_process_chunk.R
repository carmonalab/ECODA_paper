# ==============================================================================
# 2.1.1_process_chunk.R — Process one chunk of samples for cell type annotation
# ==============================================================================
# Called by 2.1_run_worker.sh (pixi run Rscript --vanilla)
# Expects a single argument: <path_to_chunk_txt>
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") stop("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")

source(file.path(project_root, "config_helper.R"))
paths <- get_pipeline_config()

# Sourcing: seurat_utils.R (get_seurat_obj_from_h5ad) — lighter than
# load_all_functions.R (no package re-imports / conflicts with scGate & HiTME).
source(file.path(project_root, "src/utils/seurat_utils.R"))

library(scGate)
library(ProjecTILs)
library(SignatuR)
library(HiTME)
library(Seurat)
library(arrow)
library(scATOMIC)
library(cutoff.scATOMIC)
library(R.utils)

library(reticulate)

ad <- import("anndata", convert = FALSE)
set.seed(123)

# ==============================================================================
# HELPER FUNCTIONS & PARAMETER PARSING
# ==============================================================================

env_sample_col    <- Sys.getenv("SAMPLE_COLNAME")
env_tissue        <- Sys.getenv("TISSUE_TYPE")
env_normal_tissue <- Sys.getenv("NORMAL_TISSUE")

defaults <- list(
  chunk_file            = NULL,
  sample_colname        = if (env_sample_col != "") env_sample_col else "Sample",
  tissue_type           = if (env_tissue != "") env_tissue else "Tumor",
  normal_tissue         = if (env_normal_tissue != "") as.logical(env_normal_tissue) else TRUE
)

hitme_cols_keep <- c("IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                     "cellCycle.G2M_UCell", "layer1", "layer2", "layer3")
scatomic_cols <- c("layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                   "scATOMIC_pred", "S.Score", "G2M.Score", "Phase", "classification_confidence")
annot_cols <- c(hitme_cols_keep, scatomic_cols)

raw_args <- commandArgs(trailingOnly = TRUE)
args <- defaults
if (length(raw_args) > 0) args$chunk_file <- raw_args[1]

if (is.null(args$chunk_file) || !file.exists(args$chunk_file)) {
  stop("Valid 'chunk_file' parameter not parsed from execution context!")
}

# Wall-clock budget: SLURM_TIME_LIMIT is in MINUTES (Slurm sets it as the
# --time value, e.g. "120" for a 2h job) — multiply by 60 for the seconds
# scale. The 7200 s fallback matches the worker's #SBATCH --time=02:00:00.
# SLURM_JOB_START_TIME is not a standard Slurm env var, so proc.time()[3]
# (R launch, i.e. after staging) is the actual elapsed source; it under-counts
# wall time by the staging duration, which is a safe (conservative) direction.
wall_limit_s <- suppressWarnings(as.numeric(Sys.getenv("SLURM_TIME_LIMIT", "")))
if (is.na(wall_limit_s) || wall_limit_s <= 0) {
  wall_limit_s <- 7200
} else {
  wall_limit_s <- wall_limit_s * 60
}
job_start_s <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_START_TIME", "")))
wall_elapsed <- function() if (!is.na(job_start_s)) as.numeric(Sys.time()) - job_start_s else proc.time()[3]
wall_left <- function() wall_limit_s - wall_elapsed()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load ref maps ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Load scGate models ####
# Load from a shared cache (created by 2.0_create_scgate_db.R via
# 2_submit_hpc_array.sh) so array workers do not all download the model DB in
# parallel. If the cache is missing, download once and persist it.
scgate_db_path <- Sys.getenv("SCGATE_DB_PATH", unset = file.path(project_root, "aux", "scGateDB.rds"))
scgate_db_branch <- Sys.getenv("SCGATE_DB_BRANCH", unset = "41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4")
if (scgate_db_path != "" && file.exists(scgate_db_path)) {
  scGate_models_DB <- readRDS(scgate_db_path)
  message("Loaded scGate DB from cache: ", scgate_db_path)
} else {
  message("scGate DB cache not found; downloading models (force_update=TRUE)...")
  scGate_models_DB <- get_scGateDB(branch = scgate_db_branch, force_update = TRUE)
  if (scgate_db_path != "") {
    dir.create(dirname(scgate_db_path), showWarnings = FALSE, recursive = TRUE)
    tmp_path <- paste0(scgate_db_path, ".tmp.", Sys.getpid())
    saveRDS(scGate_models_DB, tmp_path)
    file.rename(tmp_path, scgate_db_path)
    message("Saved scGate DB cache to: ", scgate_db_path)
  }
}
scGate_models_blood <- scGate_models_DB$human$PBMC
scGate_models_blood$MoMac <- scGate_models_blood$Monocyte
scGate_models_blood$Monocyte <- NULL
scGate_models_tumor <- scGate_models_DB$human$HiTME
scGate_models_blood <- c(scGate_models_blood, scGate_models_tumor[!names(scGate_models_tumor) %in% names(scGate_models_blood)])

### Load ProjecTILs ref maps ####
ref.maps_sketched <- list(
  CD8 = load.reference.map(file.path(paths$path_ref, "sketched_CD8T_human_ref_v1.rds")),
  CD4 = load.reference.map(file.path(paths$path_ref, "sketched_CD4T_human_ref_v2.rds")),
  DC = load.reference.map(file.path(paths$path_ref, "sketched_DC_human_ref_v2.rds")),
  MoMac = load.reference.map(file.path(paths$path_ref, "sketched_MoMac_human_v1.rds"))
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Process data ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chunk_lines <- readLines(args$chunk_file)
h5ad_file <- chunk_lines[1]
samples_to_process <- chunk_lines[-1]

### Tap h5ad file ####
adata <- ad$read_h5ad(h5ad_file, backed = "r")
obs <- py_to_r(adata$obs)

# Hardening: warn (not stop) if the on-disk X format is not CSR. anndata only
# overrides selective row-indexing for backed CSR matrices; a backed CSC matrix
# falls back to a full in-memory read per-sample subset (silent OOM).
# Files produced by 1.1.1_preprocess.py are CSR by construction, so a non-CSR
# file here means preprocessing output was replaced/regenerated elsewhere.
x_format <- tryCatch(py_to_r(adata$X$format), error = function(e) NULL)
if (!is.null(x_format) && x_format != "csr") {
  warning("On-disk X format is '", x_format, "' (expected 'csr'); backed per-sample ",
          "subsetting will fully load the matrix into memory. Re-run preprocessing ",
          "(1.1.1_preprocess.py) to force CSR on-disk.")
}

if (!args$sample_colname %in% colnames(obs)) {
  stop(args$sample_colname, " not found in h5ad obs colnames!")
}

# ==============================================================================
# PROCESSING LOOP
# ==============================================================================
message(paste("--- Starting processing for chunk file:", args$chunk_file, "---"))

# Feather name is derived from the chunk file (chunk_<N>.txt ->
# annotations_chunk_<N>.feather), NOT from SLURM_ARRAY_TASK_ID: task IDs are
# global across datasets and renumber on reruns, which would merge stale
# feathers. Chunk numbers are per-dataset and stable as long as the input
# files/samples do not change; 1.1_prepare_chunks.py deletes leftover feathers
# on every rerun.
annot_file <- file.path(
  paths$path_output,
  paste0("annotations_", sub("\\.txt$", ".feather", basename(args$chunk_file)))
)
if (file.exists(annot_file)) {
  message(paste("Chunk already processed. Annotations exist at:", annot_file))
} else {
  # Per-sample checkpoints: each sample is written to
  # output/annotation_tmp/chunk_<N>/sample_<NN>.feather as it completes, so a
  # killed/failed chunk re-runs only the missing sample(s). The final
  # annotations_chunk_<N>.feather is written atomically only when every sample
  # succeeded; annotation_tmp/ is removed on success. The subdirectory name
  # (annotation_tmp/) never matches the merge/coverage glob
  # annotations_chunk_*.feather, and 1.1_prepare_chunks.py deletes stale
  # annotation_tmp/ on every chunk rebuild (production mode).
  tmp_dir <- file.path(paths$path_output, "annotation_tmp",
                       sub("\\.txt$", "", basename(args$chunk_file)))
  failed_samples <- character(0)
  for (i in seq_along(samples_to_process)) {
    target_sample <- samples_to_process[i]
    sample_tmp <- file.path(tmp_dir, sprintf("sample_%02d.feather", i))
    if (file.exists(sample_tmp)) {
      message(paste("resume: sample", target_sample, "already annotated; skipping"))
      next
    }
    message(paste("--- Processing sample:", target_sample, "---"))

    sample_ok <- tryCatch({

    # Python 3 keys() returns a view object (dict_keys/KeysView) that py_to_r()
    # does NOT convert; materialize it to a Python list first so py_to_r returns
    # an R character vector usable with %in%.
    layer_keys <- py_to_r(import_builtins(convert = FALSE)$list(adata$layers$keys()))
    counts_layer <- if ("counts" %in% layer_keys) "counts" else "X"
    if (counts_layer == "X") {
      # Annotation-union files carry the raw counts in X by design (minimal
      # layout: no layers group — see 1.1_prepare_chunks.py), so this is the
      # designed primary path for unions, not a warning case. Preprocessed
      # view files still carry layers["counts"] and take the "counts" branch.
      message("Union carries counts in X by design (no 'counts' layer in ", h5ad_file,
              "); using X as counts input for scATOMIC/HiTME.")
    }

    seurat_obj <- get_seurat_obj_from_h5ad(
      adata, obs, target_sample,
      sample_colname = args$sample_colname,
      counts_layer = counts_layer
    )

    # HiTME (via scGate/UCell/ProjecTILs) requires the log-normalized "data"
    # layer, but CreateSeuratObject only populates "counts" — without it,
    # Run.HiTME fails with "Cannot find layer data in assay RNA" on every
    # sample. The worker object always comes fresh from CreateSeuratObject, so
    # normalize unconditionally (a v5-only Layers() guard would crash on the
    # pinned Seurat 4.4.x build, which does not export it). NormalizeData
    # leaves the counts layer untouched (scATOMIC input).
    message("  Adding log-normalized 'data' layer (NormalizeData) for HiTME...")
    seurat_obj <- NormalizeData(seurat_obj)

    timeout <- max(60, ncol(seurat_obj) / 10000 * 60 * 10)
    # Per-attempt cap (wall-clock policy: 1800 s max, and only when the
    # remaining wall time exceeds the attempt by a 300 s margin). withTimeout
    # is best-effort — it fires at R evaluation points, not inside blocking
    # python (reticulate) calls — so the real budget is enforced here at R
    # level and Slurm's wall limit stays the backstop.
    attempt_timeout <- min(timeout, 1800)

    ### scATOMIC annotation ####
    if (is.null(seurat_obj@meta.data[["layer_1"]])) {
      for (a in 1:5) {
        eff <- min(attempt_timeout, wall_left() - 300)
        if (eff < 60) {
          message(paste("  scATOMIC: insufficient wall time remaining (", round(wall_left()), "s); skipping attempt"))
          break
        }
        message(paste("  scATOMIC attempt", a, "with", round(eff), "s timeout"))
        result <- tryCatch({
          withTimeout({
            sca_preds <- run_scATOMIC(seurat_obj@assays$RNA$counts)
            sca_results <- create_summary_matrix(
              prediction_list = sca_preds,
              raw_counts = seurat_obj@assays$RNA$counts,
              normal_tissue = args$normal_tissue
            )
            "Complete"
          }, timeout = eff)
        }, TimeoutException = function(te) {
          message("  scATOMIC timeout, retrying...")
          NULL
        }, error = function(er) {
          message(paste("  scATOMIC error:", er$message, "- retrying..."))
          NULL
        })
        if (!is.null(result)) {
          sca_cols <- intersect(scatomic_cols, colnames(sca_results))
          seurat_obj <- AddMetaData(seurat_obj, sca_results[, sca_cols, drop = FALSE])
          break
        }
      }
    }

    ### HiTME annotation ####
    for (a in 1:5) {
      eff <- min(attempt_timeout, wall_left() - 300)
      if (eff < 60) {
        message(paste("  HiTME: insufficient wall time remaining (", round(wall_left()), "s); skipping attempt"))
        break
      }
      message(paste("  HiTME attempt", a, "with", round(eff), "s timeout"))
      result <- tryCatch({
        withTimeout({
          if (args$tissue_type == "Blood") {
            seurat_obj <- Run.HiTME(
              object = seurat_obj,
              scGate.model = scGate_models_blood,
              ref.maps = ref.maps_sketched,
              verbose = FALSE,
              ncores = 1
            )
          } else {
            seurat_obj <- Run.HiTME(
              object = seurat_obj,
              scGate.model = scGate_models_tumor,
              ref.maps = ref.maps_sketched,
              verbose = FALSE,
              ncores = 1
            )
          }
          "Complete"
        }, timeout = eff)
      }, TimeoutException = function(te) {
        message("  HiTME timeout, retrying...")
        NULL
      }, error = function(er) {
        message(paste("  HiTME error:", er$message, "- retrying..."))
        NULL
      })
      if (!is.null(result)) break
    }

    ### Extract annotations ####
    meta <- seurat_obj@meta.data
    keep_cols <- intersect(annot_cols, colnames(meta))
    annot <- meta[, keep_cols, drop = FALSE]
    annot$cell_barcode <- rownames(annot)
    annot[[args$sample_colname]] <- target_sample

    rm(seurat_obj)
    gc()
    TRUE
    }, error = function(e) {
      message(paste("SAMPLE FAILED (", target_sample, "): ", conditionMessage(e), " - continuing"))
      FALSE
    })

    if (sample_ok) {
      dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
      sample_tmp_write <- paste0(sample_tmp, ".tmp.", Sys.getpid())
      write_feather(annot, sample_tmp_write)
      file.rename(sample_tmp_write, sample_tmp)
      message(paste("  Checkpoint written:", sample_tmp))
    } else {
      failed_samples <- c(failed_samples, target_sample)
    }
  }

  if (length(failed_samples) > 0) {
    message(paste("chunk INCOMPLETE - failed samples:", paste(failed_samples, collapse = ", ")))
    message(paste("Sample intermediates kept in:", tmp_dir, "- a re-run resumes only the failed sample(s)."))
    quit(status = 1)
  }

  annotations_list <- list()
  for (i in seq_along(samples_to_process)) {
    annotations_list[[samples_to_process[i]]] <- read_feather(
      file.path(tmp_dir, sprintf("sample_%02d.feather", i))
    )
  }
  annotations_df <- do.call(rbind, annotations_list)
  rownames(annotations_df) <- NULL
  annot_file_tmp <- paste0(annot_file, ".tmp.", Sys.getpid())
  write_feather(annotations_df, annot_file_tmp)
  file.rename(annot_file_tmp, annot_file)
  unlink(tmp_dir, recursive = TRUE)
  message(paste("Wrote annotations to:", annot_file))
}

message(paste("---", args$chunk_file, "processing complete! ---"))
