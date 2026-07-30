# ==============================================================================
# 2.2_process_chunk.R — Process one chunk of samples for cell type annotation
# ==============================================================================
# Called by 2.2_process_chunk.sh (pixi run Rscript --vanilla)
# Expects a single argument: chunk_file__<path_to_chunk_txt>
# ==============================================================================

project_root <- getwd()

source(file.path(project_root, "config_helper.R"))
paths <- get_pipeline_config()

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
get_sample_seurat_obj <- function(adata, r_obs, target_sample, sample_colname) {
  sample_indices <- which(r_obs[[sample_colname]] == target_sample) - 1

  subset_py <- adata[as.integer(sample_indices)]
  raw_X_py <- subset_py$X$astype("float64")$tocsc()

  counts_matrix <- py_to_r(raw_X_py)

  counts_matrix <- t(counts_matrix)
  rownames(counts_matrix) <- as.character(py_to_r(subset_py$var_names$values))
  colnames(counts_matrix) <- as.character(py_to_r(subset_py$obs_names$values))

  seurat_obj <- CreateSeuratObject(
    counts = counts_matrix,
    meta.data = as.data.frame(py_to_r(subset_py$obs))
  )
  return(seurat_obj)
}

env_sample_col    <- Sys.getenv("SAMPLE_COLNAME")
env_tissue        <- Sys.getenv("TISSUE_TYPE")
env_auth_annots   <- Sys.getenv("AUTHOR_ANNOT_COLNAMES")
env_normal_tissue <- Sys.getenv("NORMAL_TISSUE")

defaults <- list(
  chunk_file            = NULL,
  sample_colname        = if (env_sample_col != "") env_sample_col else "sample",
  tissue_type           = if (env_tissue != "") env_tissue else "Tumor",
  author_annot_colnames = if (env_auth_annots != "") unlist(strsplit(env_auth_annots, ",")) else character(),
  normal_tissue         = if (env_normal_tissue != "") as.logical(env_normal_tissue) else TRUE
)

hitme_cols_keep <- c("IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                     "cellCycle.G2M_UCell", "layer1", "layer2", "layer3")
scatomic_cols <- c("layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                   "scATOMIC_pred", "S.Score", "G2M.Score", "Phase", "classification_confidence")
annot_cols <- c(hitme_cols_keep, scatomic_cols)

raw_args <- commandArgs(trailingOnly = TRUE)
args <- defaults

if (length(raw_args) > 0) {
  parsed_args_list <- unlist(strsplit(raw_args[1], "__"))
  keys <- parsed_args_list[seq(1, length(parsed_args_list), by = 2)]
  vals <- parsed_args_list[seq(2, length(parsed_args_list), by = 2)]

  args_list <- as.list(vals)
  names(args_list) <- keys
  args <- modifyList(defaults, args_list)
}

if (is.null(args$chunk_file) || !file.exists(args$chunk_file)) {
  stop("Valid 'chunk_file' parameter not parsed from execution context!")
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load ref maps ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Load scGate models ####
scGate_models_DB <- get_scGateDB(branch = "dev", verbose = TRUE, force_update = TRUE)
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

if (!args$sample_colname %in% colnames(obs)) {
  stop(args$sample_colname, " not found in h5ad obs colnames!")
}

# ==============================================================================
# PROCESSING LOOP
# ==============================================================================
message(paste("--- Starting processing for chunk file:", args$chunk_file, "---"))

chunk_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
annot_file <- file.path(paths$path_output, paste0("annotations_chunk_", chunk_id, ".feather"))
if (file.exists(annot_file)) {
  message(paste("Chunk already processed. Annotations exist at:", annot_file))
} else {
  annotations_list <- list()
  for (target_sample in samples_to_process) {
    message(paste("--- Processing sample:", target_sample, "---"))

    seurat_obj <- get_sample_seurat_obj(
      adata, obs, target_sample, args$sample_colname
    )

    timeout <- max(60, ncol(seurat_obj) / 10000 * 60 * 10)

    ### scATOMIC annotation ####
    if (is.null(seurat_obj@meta.data[["layer_1"]])) {
      for (a in 1:5) {
        message(paste("  scATOMIC attempt", a, "with", round(timeout), "s timeout"))
        result <- tryCatch({
          withTimeout({
            sca_preds <- run_scATOMIC(seurat_obj@assays$RNA$counts)
            sca_results <- create_summary_matrix(
              prediction_list = sca_preds,
              raw_counts = seurat_obj@assays$RNA$counts,
              normal_tissue = args$normal_tissue
            )
            "Complete"
          }, timeout = timeout)
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
      message(paste("  HiTME attempt", a, "with", timeout, "s timeout"))
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
        }, timeout = timeout)
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
    annot$sample <- target_sample
    annotations_list[[target_sample]] <- annot

    rm(seurat_obj)
    gc()
  }

  annotations_df <- do.call(rbind, annotations_list)
  rownames(annotations_df) <- NULL
  write_feather(annotations_df, annot_file)
  message(paste("Wrote annotations to:", annot_file))
}

message(paste("---", args$chunk_file, "processing complete! ---"))
