# ==============================================================================
# 1.6.1_fill_bassez_cellsubtype.R — Fill missing cellSubType with cellType
# ==============================================================================
# In the Bassez dataset, original author subclustering (cellSubType) was only
# performed for T cells and Myeloid cells. All other cells (Cancer cells,
# Fibroblasts, B cells, Endothelial, Mast, etc. — ~68% of the dataset) have
# cellSubType == NA.
#
# Following the original preprocessing pipeline, missing cellSubType annotations
# are filled with the corresponding broad cellType annotation to preserve
# distinct non-immune and B-cell populations in high-resolution compositions.
#
# Must run AFTER 1_stage_data.sh and BEFORE the preprocess array
# (src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh).
# ==============================================================================

suppressPackageStartupMessages(library(Seurat))

script_file_arg <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_file_arg) != 1) {
  stop("Could not determine the Bassez subtype patch script directory.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg)))
source(file.path(script_dir, "bassez_cellsubtype_utils.R"))

scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR")
if (scratch_dir == "") {
  stop("CRITICAL: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling this script.")
}

input <- file.path(scratch_dir, "Bassez", "data", "BassezA_2021_33958794whole.rds")
if (!file.exists(input)) {
  stop("Bassez input RDS not found at: ", input)
}

seurat <- readRDS(input)
if (!inherits(seurat, "Seurat")) {
  stop("Bassez input is not a Seurat object: ", input)
}

metadata_before <- seurat@meta.data
if (nrow(metadata_before) != ncol(seurat)) {
  stop(
    "Bassez metadata row count (", nrow(metadata_before),
    ") does not match cell count (", ncol(seurat), ")."
  )
}
required_metadata_columns <- c("cellType", "cellSubType")
missing_metadata_columns <- setdiff(required_metadata_columns, colnames(metadata_before))
if (length(missing_metadata_columns) > 0) {
  stop(
    "Bassez metadata is missing required column(s): ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

replacement_rows <- bassez_missing_annotation(metadata_before$cellSubType)
n_replaced <- sum(replacement_rows)
before_valid_subtypes <- as.character(metadata_before$cellSubType)[!replacement_rows]
before_cardinality <- length(unique(before_valid_subtypes))
metadata_columns <- colnames(metadata_before)
metadata_rownames <- rownames(metadata_before)
expected_cell_count <- ncol(seurat)
original_cell_type <- as.character(metadata_before$cellType)

message(
  "Filling ", n_replaced, "/", expected_cell_count,
  " missing cellSubType entries with cellType in Bassez..."
)

seurat@meta.data <- bassez_fill_cell_subtype(metadata_before)

# Write and validate on the same filesystem before replacing the staged source.
tmp_input <- tempfile(
  pattern = paste0(".", basename(input), ".tmp-"),
  tmpdir = dirname(input),
  fileext = ".rds"
)
tmp_created <- TRUE
on.exit({
  if (isTRUE(tmp_created) && file.exists(tmp_input)) unlink(tmp_input)
}, add = TRUE)

saveRDS(seurat, tmp_input)
validated <- readRDS(tmp_input)
if (!inherits(validated, "Seurat")) {
  stop("Temporary Bassez RDS did not reopen as a Seurat object.")
}

validated_metadata <- validated@meta.data
if (ncol(validated) != expected_cell_count || nrow(validated_metadata) != expected_cell_count) {
  stop("Temporary Bassez RDS failed cell-count validation.")
}
if (!identical(colnames(validated_metadata), metadata_columns) ||
    !identical(rownames(validated_metadata), metadata_rownames)) {
  stop("Temporary Bassez RDS failed metadata index/column validation.")
}
if (!identical(as.character(validated_metadata$cellType), original_cell_type)) {
  stop("Temporary Bassez RDS changed cellType metadata unexpectedly.")
}
if (!is.factor(validated_metadata$cellSubType)) {
  stop("Temporary Bassez RDS did not preserve categorical cellSubType metadata.")
}
if (any(bassez_missing_annotation(validated_metadata$cellSubType))) {
  stop("Temporary Bassez RDS still contains missing or sentinel cellSubType values.")
}
after_valid_subtypes <- as.character(validated_metadata$cellSubType)[!replacement_rows]
if (!identical(after_valid_subtypes, before_valid_subtypes)) {
  stop("Temporary Bassez RDS changed valid cellSubType values unexpectedly.")
}

after_cardinality <- length(unique(as.character(validated_metadata$cellSubType)))
if (!file.rename(tmp_input, input)) {
  stop("Could not atomically install validated Bassez RDS over: ", input)
}
tmp_created <- FALSE

message(
  "Installed updated Bassez Seurat object: ", input,
  " (replacements=", n_replaced,
  ", subtype_cardinality=", before_cardinality, " -> ", after_cardinality, ")"
)
