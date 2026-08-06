# ==============================================================================
# 1.2.1_prepare_joanito.R — Joanito seqtec batch column + _debug 5-sample subset
# ==============================================================================
# HPC script (submitted via 1_submit_hpc.sh -> 1.2_submit_joanito.sh). Reads
# the staged Joanito raw .rds into memory once and:
#
#   1. Computes the `seqtec` batch column and saves the full object back in
#      place (only when the column is absent — on re-runs the in-place save is
#      skipped). Must run AFTER 1_stage_data.sh and BEFORE the preprocess array
#      (the batch_effect_analysis view uses seqtec as batch_col).
#
#   2. Derives the _debug 5-sample subset from the SAME in-memory object:
#      5 samples covering (sample.origin x seqtec x Site) combos (candidates
#      with >= 500 cells), 500 cells/sample (random, seeded 321), minimal obs
#      columns (incl. seqtec, Site, sample/patient cols), written as
#      ${HPC_SCRATCH_DIR}/_debug/data/JoaI_2022_35773407_debug_5samples.h5ad
#      (anndataR: X=None + layers["counts"] — handled by the X-promotion in
#      1.1.1_preprocess.py).
#
# NOTE: This script is the SINGLE SOURCE OF TRUTH for the seqtec batch mapping
# (dataset -> "5' seq"/"3' seq"). If the mapping changes, only this file
# changes. The _debug subset inherits the mapping from the in-memory object.
# ==============================================================================

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(anndataR))

scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR")
if (scratch_dir == "") {
  stop("CRITICAL: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling this script.")
}
config_path <- Sys.getenv("DATASETS_JSON_FILE")
if (config_path == "") {
  stop("CRITICAL: DATASETS_JSON_FILE not set. Source slurm_config.sh before calling this script.")
}
config <- fromJSON(config_path, simplifyVector = FALSE)

joanito_input <- config[["Joanito"]][["views"]][["batch_effect_analysis"]][["input_file_name"]]
if (is.null(joanito_input)) {
  stop("No input_file_name for the Joanito batch_effect_analysis view in datasets.json.")
}
input <- file.path(scratch_dir, "Joanito", "data", joanito_input)
message("Reading Joanito raw input: ", input)

seurat <- readRDS(input)
message("Input has ", ncol(seurat), " cells in ",
        length(unique(seurat$sample.ID)), " samples.")

# -----------------------------------------------------------------------------
# 1. seqtec batch column (single source of truth — see header note)
# -----------------------------------------------------------------------------
if (!"seqtec" %in% colnames(seurat@meta.data)) {
  seurat$seqtec <- ifelse(
    seurat$dataset %in% c("CRC-SG1", "KUL5"),
    "5' seq",
    "3' seq"
  )
  # In-place save is idempotent. Must run AFTER 1_stage_data.sh and BEFORE
  # the preprocess array.
  saveRDS(seurat, input)
  message("Saved seqtec back to: ", input)
} else {
  message("seqtec already present — skipping recompute and in-place save.")
}

# -----------------------------------------------------------------------------
# 2. _debug 5-sample subset (same in-memory object, no extra read)
# -----------------------------------------------------------------------------
set.seed(321)

# Sample selection: 5 samples covering both biological conditions (sample.origin)
# and batches (seqtec), preferring samples with >= min_cells cells
select_debug_samples <- function(meta, n, min_cells) {
  counts <- sort(table(meta$sample.ID), decreasing = TRUE)
  candidates <- names(counts)[counts >= min_cells]
  if (length(candidates) < n) {
    stop(sprintf(
      "Only %d samples have >= %d cells, but %d are requested.",
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

chosen <- select_debug_samples(seurat@meta.data, 5, 500)
message("Chosen samples: ", paste(chosen, collapse = ", "))

meta <- seurat@meta.data
set.seed(321)
keep_cells <- unlist(lapply(chosen, function(s) {
  s_cells <- rownames(meta)[meta$sample.ID == s]
  sample(s_cells, min(500, length(s_cells)))
}))
message("Keeping ", length(keep_cells), " cells (500 per sample x ",
        length(chosen), " samples).")

seurat <- seurat[, keep_cells]
meta <- seurat@meta.data

# Keep minimal obs columns (incl. seqtec, Site, sample/patient cols)
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

# anndataR writes X=None + layers["counts"] for Seurat v5 objects; the
# preprocess pipeline promotes the counts layer to X (see base_preprocessing()).
# Hardcoded stem MUST match datasets.json _debug.file_names / view
# input_file_name. Raw subset stays in _debug/data/ only — never copy it into
# _debug/output/ (the merge script treats every output/*.h5ad as a view).
debug_dir <- file.path(scratch_dir, "_debug", "data")
dir.create(debug_dir, showWarnings = FALSE, recursive = TRUE)
h5ad_path <- file.path(debug_dir, "JoaI_2022_35773407_debug_5samples.h5ad")
write_h5ad(seurat, h5ad_path)
message("Saved: ", h5ad_path)

# -----------------------------------------------------------------------------
# Summary for verification
# -----------------------------------------------------------------------------
summary_table <- unique(meta[, intersect(c("sample.ID", "sample.origin", "seqtec", "Site"), colnames(meta))])
summary_table$n_cells <- as.integer(table(meta$sample.ID)[summary_table$sample.ID])
print(summary_table, row.names = FALSE)

message("Debug dataset created: ", ncol(seurat), " cells across ",
        length(unique(meta$sample.ID)), " samples.")
