# ==============================================================================
# 1.3.1_prepare_joanito.R — Joanito seqtec batch column + _debug 5-sample subset
# ==============================================================================
# HPC script (submitted via 1_submit_hpc.sh -> 1.3_submit_joanito.sh). Reads
# the staged Joanito raw .rds into memory once and:
#
#   1. Computes the `seqtec` batch column and saves the full object back in
#      place (only when the column is absent — on re-runs the in-place save is
#      skipped). Must run AFTER 1_stage_data.sh and BEFORE the preprocess array
#      (the batch_effect_uncorrected/corrected views use seqtec as batch_col).
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

joanito_input <- config[["Joanito"]][["views"]][["batch_effect_uncorrected"]][["input_file_name"]]
if (is.null(joanito_input)) {
  stop("No input_file_name for the Joanito batch-effect view in datasets.json.")
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

# Sample selection: Balanced cohort (12 samples, 500 cells each = 6,000 cells):
# 1. seqtec: exactly 50% 3' seq (6 samples) and 50% 5' seq (6 samples).
# 2. sample.origin: Normal and Tumor present in BOTH 3' seq and 5' seq
#    (allowing direct assessment of batch vs bio separation within identical cell types/states).
# 3. Site: matched anatomical locations (Ascending colon, Sigmoid colon, Rectum) across technologies.
select_debug_samples <- function(meta, min_cells = 500) {
  counts <- sort(table(meta$sample.ID), decreasing = TRUE)
  candidates <- names(counts)[counts >= min_cells]

  info <- unique(meta[meta$sample.ID %in% candidates,
                      c("sample.ID", "sample.origin", "seqtec", "Site")])
  info$n_cells <- as.integer(counts[info$sample.ID])
  info <- info[order(-info$n_cells), ]

  pick_samples <- function(st, bo, sites, max_n) {
    chosen_sub <- character(0)
    for (st_loc in sites) {
      if (length(chosen_sub) >= max_n) break
      cand <- info$sample.ID[info$seqtec == st & info$sample.origin == bo & info$Site == st_loc & !info$sample.ID %in% chosen_sub]
      if (length(cand) > 0) chosen_sub <- c(chosen_sub, cand[1])
    }
    if (length(chosen_sub) < max_n) {
      cand <- info$sample.ID[info$seqtec == st & info$sample.origin == bo & !info$sample.ID %in% chosen_sub]
      needed <- max_n - length(chosen_sub)
      chosen_sub <- c(chosen_sub, cand[seq_len(min(needed, length(cand)))])
    }
    return(chosen_sub)
  }

  sites_pref <- c("Ascending colon", "Sigmoid colon", "Rectum", "Upper rectum", "Caecum")

  s_3p_norm  <- pick_samples("3' seq", "Normal", sites_pref, 3)
  s_3p_tumor <- pick_samples("3' seq", "Tumor",  sites_pref, 3)

  s_5p_norm  <- pick_samples("5' seq", "Normal", sites_pref, 3)
  s_5p_tumor <- pick_samples("5' seq", "Tumor",  sites_pref, 2)
  s_5p_ln    <- pick_samples("5' seq", "LymphNode", sites_pref, 1)

  chosen <- c(s_3p_norm, s_3p_tumor, s_5p_norm, s_5p_tumor, s_5p_ln)
  chosen_meta <- info[info$sample.ID %in% chosen, ]

  message("\n=== Selected Debug Samples (", length(chosen), " samples) ===")
  print(table(chosen_meta$seqtec, chosen_meta$sample.origin))
  print(chosen_meta[order(chosen_meta$seqtec, chosen_meta$sample.origin, chosen_meta$Site),
                    c("sample.ID", "seqtec", "sample.origin", "Site", "n_cells")], row.names = FALSE)

  # Strict validation
  if (sum(chosen_meta$seqtec == "3' seq") != sum(chosen_meta$seqtec == "5' seq")) {
    stop("select_debug_samples: seqtec is not 50/50 balanced!")
  }
  if (!all(c("Normal", "Tumor") %in% chosen_meta$sample.origin[chosen_meta$seqtec == "3' seq"]) ||
      !all(c("Normal", "Tumor") %in% chosen_meta$sample.origin[chosen_meta$seqtec == "5' seq"])) {
    stop("select_debug_samples: Normal and Tumor must be present in BOTH 3' and 5' seq!")
  }

  return(chosen)
}

chosen <- select_debug_samples(seurat@meta.data, 500)
message("\nChosen sample IDs (", length(chosen), "): ", paste(chosen, collapse = ", "))

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
write_h5ad(seurat, h5ad_path, mode = "w")
message("Saved: ", h5ad_path)

# -----------------------------------------------------------------------------
# Summary for verification
# -----------------------------------------------------------------------------
summary_table <- unique(meta[, intersect(c("sample.ID", "sample.origin", "seqtec", "Site"), colnames(meta))])
summary_table$n_cells <- as.integer(table(meta$sample.ID)[summary_table$sample.ID])
print(summary_table, row.names = FALSE)

message("Debug dataset created: ", ncol(seurat), " cells across ",
        length(unique(meta$sample.ID)), " samples.")
