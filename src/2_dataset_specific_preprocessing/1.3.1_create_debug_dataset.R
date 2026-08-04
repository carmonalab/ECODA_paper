# ==============================================================================
# 1.3.1_create_debug_dataset.R — Build the 5-sample Joanito debug subset
# ==============================================================================
# LOCAL script: runs on the user's computer, NOT on the HPC. It is deliberately
# NOT part of the 1_submit_hpc.sh dispatcher glob (which matches 1.*_submit_*.sh
# only), so no 1.3_submit_*.sh wrapper exists for it.
#
# Reads the Joanito raw input declared in datasets.json (input_file_name of the
# batch_effect_analysis view), selects 5 samples covering both biological
# conditions (sample.origin) and batches (seqtec), subsets 500 cells per sample
# (random, seeded), keeps minimal obs columns (incl. seqtec, Site, sample and
# patient columns), and writes:
#   data/debug/JoaI_2022_35773407_debug_5samples.rds   (Seurat, for R workflows)
#   data/debug/JoaI_2022_35773407_debug_5samples.h5ad  (anndata, X=None +
#       layers["counts"] — handled by 1.1.1_preprocess.py's X=None promotion)
#
# Before the HPC debug run (Phase 2), the user must place these files on the
# NAS under Standardized_SingleCell_Datasets/debug/output/ (the _debug entry's
# folder_name), or scp/rsync them to ${HPC_SCRATCH_DIR}/_debug/data/ manually,
# because 1_stage_data.sh reads only from the NAS.
#
# Usage:
#   pixi run Rscript src/2_dataset_specific_preprocessing/1.3.1_create_debug_dataset.R
#   pixi run Rscript src/2_dataset_specific_preprocessing/1.3.1_create_debug_dataset.R \
#       --input /path/to/JoaI_2022_35773407_Nofilt_whole.rds \
#       --output-dir data/debug --n-samples 5 --cells-per-sample 500 --seed 321
# ==============================================================================

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(anndataR))

# -----------------------------------------------------------------------------
# Argument parsing (plain R, no extra deps)
# -----------------------------------------------------------------------------
parse_args <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  args <- list(
    input = NULL,
    output_dir = file.path("data", "debug"),
    n_samples = 5,
    cells_per_sample = 500,
    seed = 321
  )
  i <- 1
  while (i <= length(raw)) {
    key <- raw[i]
    val <- raw[i + 1]
    if (key == "--input") {
      args$input <- val
    } else if (key == "--output-dir") {
      args$output_dir <- val
    } else if (key == "--n-samples") {
      args$n_samples <- as.integer(val)
    } else if (key == "--cells-per-sample") {
      args$cells_per_sample <- as.integer(val)
    } else if (key == "--seed") {
      args$seed <- as.integer(val)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 2
  }
  args
}

args <- parse_args()
set.seed(args$seed)

# -----------------------------------------------------------------------------
# Locate the Joanito raw input (CLI arg > local data/ copy > NAS)
# -----------------------------------------------------------------------------
config <- fromJSON("datasets.json", simplifyVector = FALSE)
joanito_input <- config[["Joanito"]][["views"]][["batch_effect_analysis"]][["input_file_name"]]
if (is.null(joanito_input)) {
  stop("No input_file_name for the Joanito batch_effect_analysis view in datasets.json.")
}

if (!is.null(args$input)) {
  input_path <- args$input
} else {
  local_candidates <- file.path("data", joanito_input)
  nas_candidates <- file.path(
    "/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets",
    config[["Joanito"]][["folder_name"]],
    "output",
    joanito_input
  )
  input_path <- local_candidates
  if (!file.exists(input_path)) {
    message("Local input not found at ", input_path, "; trying NAS.")
    input_path <- nas_candidates
  }
}
if (!file.exists(input_path)) {
  stop(
    "Joanito raw input not found. Expected at:\n  ",
    file.path("data", joanito_input),
    "\nor:\n  ", nas_candidates,
    "\nPass --input <path> to use another location."
  )
}
message("Reading Joanito raw input: ", input_path)

seurat <- readRDS(input_path)
message("Input has ", ncol(seurat), " cells in ", length(unique(seurat$sample.ID)), " samples.")

# -----------------------------------------------------------------------------
# Compute the seqtec batch column. NOTE: this mapping must stay in sync with
# 1.2.1_create_joanito_batch_col.R (single source of truth for the batch
# definition of the full Joanito dataset) — update both if it changes.
# -----------------------------------------------------------------------------
seurat$seqtec <- ifelse(
  seurat$dataset %in% c("CRC-SG1", "KUL5"),
  "5' seq",
  "3' seq"
)

# -----------------------------------------------------------------------------
# Sample selection: 5 samples covering both biological conditions (sample.origin)
# and batches (seqtec), preferring samples with >= min_cells cells
# -----------------------------------------------------------------------------
select_debug_samples <- function(meta, n, min_cells) {
  counts <- sort(table(meta$sample.ID), decreasing = TRUE)
  candidates <- names(counts)[counts >= min_cells]
  if (length(candidates) < n) {
    stop(sprintf(
      "Only %d samples have >= %d cells, but %d are requested. Lower --min-cells or --n-samples.",
      length(candidates), min_cells, n
    ))
  }

  info <- unique(meta[meta$sample.ID %in% candidates,
                      c("sample.ID", "sample.origin", "seqtec", "Site")])
  info <- info[order(-counts[info$sample.ID]), ]
  # Diversity key: one sample per (condition x batch x site) combo first
  info$combo <- paste(info$sample.origin, info$seqtec, info$Site, sep = "|")

  chosen <- character(0)
  for (k in unique(info$combo)) {
    cand <- info$sample.ID[info$combo == k & !info$sample.ID %in% chosen]
    if (length(cand) > 0) chosen <- c(chosen, cand[1])
    if (length(chosen) >= n) break
  }
  if (length(chosen) < n) {
    rest <- info$sample.ID[!info$sample.ID %in% chosen]
    chosen <- c(chosen, rest[seq_len(n - length(chosen))])
  }
  chosen[seq_len(n)]
}

chosen <- select_debug_samples(seurat@meta.data, args$n_samples, 500)
message("Chosen samples: ", paste(chosen, collapse = ", "))

meta <- seurat@meta.data
set.seed(args$seed)
keep_cells <- unlist(lapply(chosen, function(s) {
  s_cells <- rownames(meta)[meta$sample.ID == s]
  sample(s_cells, min(args$cells_per_sample, length(s_cells)))
}))
message("Keeping ", length(keep_cells), " cells (",
        args$cells_per_sample, " per sample x ", length(chosen), " samples).")

seurat <- seurat[, keep_cells]
meta <- seurat@meta.data

# -----------------------------------------------------------------------------
# Keep minimal obs columns (incl. seqtec, Site, sample/patient cols)
# -----------------------------------------------------------------------------
cols_keep <- intersect(
  c("sample.ID", "sample.origin", "patient.ID", "iCMS",
    "dataset", "Gender", "Site", "cell.type", "seqtec"),
  colnames(meta)
)
if (length(cols_keep) < 5) {
  warning("Fewer than 5 of the intended metadata columns were found: ",
          paste(cols_keep, collapse = ", "))
}
seurat@meta.data <- meta[, cols_keep, drop = FALSE]

# Clean v5 object (counts + metadata only), consistent with the pipeline's
# create_clean_seuratv5_object() step for rds inputs.
seurat <- CreateSeuratObject(
  counts = seurat@assays$RNA$counts,
  meta.data = seurat@meta.data
)

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
stem <- "JoaI_2022_35773407_debug_5samples"
rds_path <- file.path(args$output_dir, paste0(stem, ".rds"))
h5ad_path <- file.path(args$output_dir, paste0(stem, ".h5ad"))

saveRDS(seurat, rds_path)
message("Saved: ", rds_path)

# anndataR writes X=None + layers["counts"] for Seurat v5 objects; the
# preprocess pipeline promotes the counts layer to X (see base_preprocessing()).
write_h5ad(seurat, h5ad_path)
message("Saved: ", h5ad_path)

# -----------------------------------------------------------------------------
# Summary for verification
# -----------------------------------------------------------------------------
meta <- seurat@meta.data
summary_table <- unique(meta[, intersect(c("sample.ID", "sample.origin", "seqtec", "Site"), colnames(meta))])
summary_table$n_cells <- as.integer(table(meta$sample.ID)[summary_table$sample.ID])
print(summary_table, row.names = FALSE)

message("Debug dataset created: ", ncol(seurat), " cells across ",
        length(unique(meta$sample.ID)), " samples.")
