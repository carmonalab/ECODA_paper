# ==============================================================================
# 1.3.1_prepare_joanito.R — Joanito seqtec batch column + _debug 5-sample subset
# ==============================================================================
# HPC script (submitted via 1_submit_hpc.sh -> 1.3_submit_joanito.sh). Reads
# the staged Joanito raw .rds into memory once and:
#
#   1. Computes deterministic `seqtec` and idempotent `cell.type_new` metadata
#      from the raw `dataset`, `cell.type`, and `iCMS` columns. The full RDS is
#      atomically rewritten whenever either derived column is missing or stale.
#
#   2. Derives the _debug 5-sample subset from the SAME in-memory object:
#      exactly one sample for each required (seqtec x sample.origin) target,
#      including the preferred 5' seq/LymphNode sample. Each selected sample
#      contributes exactly 500 cells, and the result is atomically written as
#      ${HPC_SCRATCH_DIR}/_debug/data/JoaI_2022_35773407_debug_5samples.h5ad.
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
# 1. Derived metadata (single source of truth — see header note)
# -----------------------------------------------------------------------------
meta_full <- seurat@meta.data
required_raw <- c("dataset", "cell.type", "iCMS")
missing_raw <- setdiff(required_raw, colnames(meta_full))
if (length(missing_raw) > 0) {
  stop("Joanito raw metadata is missing required columns: ",
       paste(missing_raw, collapse = ", "))
}

expected_seqtec <- ifelse(
  as.character(meta_full$dataset) %in% c("CRC-SG1", "KUL5"),
  "5' seq",
  "3' seq"
)
derive_cell_type_new <- function(cell_type, iCMS) {
  cell_type_base <- as.character(cell_type)
  iCMS_present <- !is.na(iCMS)
  expected <- cell_type_base
  expected[iCMS_present] <- paste0(
    cell_type_base[iCMS_present],
    "_",
    ifelse(as.character(iCMS[iCMS_present]) == "Normal",
           "Normal", "Cancer")
  )
  expected
}
expected_cell_type_new <- derive_cell_type_new(
  meta_full$cell.type, meta_full$iCMS
)

column_matches <- function(meta, column, expected) {
  if (!column %in% colnames(meta) || length(meta[[column]]) != length(expected)) {
    return(FALSE)
  }
  actual <- as.character(meta[[column]])
  same_na <- is.na(actual) & is.na(expected)
  same_value <- !is.na(actual) & !is.na(expected) & actual == expected
  all(same_na | same_value)
}

seqtec_current <- column_matches(meta_full, "seqtec", expected_seqtec)
cell_type_new_current <- column_matches(
  meta_full, "cell.type_new", expected_cell_type_new
)
meta_full$seqtec <- expected_seqtec
meta_full$cell.type_new <- expected_cell_type_new
seurat@meta.data <- meta_full

if (!seqtec_current || !cell_type_new_current) {
  # Install the staged RDS atomically. A killed saveRDS must never leave a
  # readable-looking but truncated source for the next numbered stage.
  tmp_rds <- tempfile(pattern = paste0(".", basename(input), ".tmp-"),
                      tmpdir = dirname(input), fileext = ".rds")
  saveRDS(seurat, tmp_rds)
  if (!file.exists(tmp_rds) || file.info(tmp_rds)$size <= 0) {
    unlink(tmp_rds)
    stop("Atomic Joanito RDS write produced an empty file: ", tmp_rds)
  }
  if (!file.rename(tmp_rds, input)) {
    unlink(tmp_rds)
    stop("Could not atomically install Joanito RDS: ", input)
  }
  message("Saved repaired seqtec/cell.type_new metadata back to: ", input)
} else {
  message("seqtec and cell.type_new are current — skipping full-RDS rewrite.")
}

# The _debug subset is always regenerated from this repaired in-memory object.

# -----------------------------------------------------------------------------
# 2. _debug 5-sample subset (same in-memory object, no extra read)
# -----------------------------------------------------------------------------
set.seed(321)

# Sample selection: exactly five samples, each with at least 500 cells:
# one each for the four technology/state combinations and the preferred
# 5' seq/LymphNode sample. Site preference only breaks ties deterministically.
select_debug_samples <- function(meta, min_cells = 500) {
  required <- c("sample.ID", "sample.origin", "seqtec", "Site")
  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) {
    stop("select_debug_samples: missing metadata columns: ",
         paste(missing, collapse = ", "))
  }

  if (any(is.na(meta$sample.ID)) ||
      any(!nzchar(trimws(as.character(meta$sample.ID))))) {
    stop("select_debug_samples: sample.ID contains missing or blank values.")
  }
  sample_ids <- as.character(meta$sample.ID)
  counts <- sort(table(sample_ids), decreasing = TRUE)
  eligible <- names(counts)[as.integer(counts) >= min_cells]
  if (length(eligible) < 5) {
    stop("select_debug_samples: fewer than five samples have at least ",
         min_cells, " cells.")
  }

  for (sample_id in eligible) {
    sample_rows <- meta[sample_ids == sample_id, , drop = FALSE]
    for (column in c("sample.origin", "seqtec", "Site")) {
      values <- as.character(sample_rows[[column]])
      if (any(is.na(values)) || any(!nzchar(trimws(values))) ||
          length(unique(values)) != 1L) {
        stop("select_debug_samples: sample ", sample_id,
             " has inconsistent ", column, " metadata.")
      }
    }
  }

  info <- meta[match(eligible, sample_ids), required, drop = FALSE]
  info$sample.ID <- as.character(info$sample.ID)
  info$sample.origin <- as.character(info$sample.origin)
  info$seqtec <- as.character(info$seqtec)
  info$Site <- as.character(info$Site)
  info$n_cells <- as.integer(counts[info$sample.ID])

  sites_pref <- c(
    "Ascending colon", "Sigmoid colon", "Rectum", "Upper rectum", "Caecum"
  )
  pick_sample <- function(seqtec_value, origin_value) {
    candidates <- info[
      info$seqtec == seqtec_value & info$sample.origin == origin_value,
      , drop = FALSE
    ]
    if (nrow(candidates) == 0) {
      stop("select_debug_samples: no eligible sample for (",
           seqtec_value, ", ", origin_value, ").")
    }
    site_rank <- match(candidates$Site, sites_pref)
    site_rank[is.na(site_rank)] <- length(sites_pref) + 1L
    candidates <- candidates[order(
      site_rank, -candidates$n_cells, candidates$sample.ID
    ), , drop = FALSE]
    candidates$sample.ID[[1]]
  }

  chosen <- c(
    pick_sample("3' seq", "Normal"),
    pick_sample("3' seq", "Tumor"),
    pick_sample("5' seq", "Normal"),
    pick_sample("5' seq", "Tumor"),
    pick_sample("5' seq", "LymphNode")
  )
  if (length(chosen) != 5L || anyDuplicated(chosen)) {
    stop("select_debug_samples: selected sample IDs are not exactly five unique IDs.")
  }

  chosen_meta <- info[match(chosen, info$sample.ID), , drop = FALSE]
  expected_pairs <- c(
    "3' seq\tNormal", "3' seq\tTumor",
    "5' seq\tNormal", "5' seq\tTumor", "5' seq\tLymphNode"
  )
  actual_pairs <- paste(chosen_meta$seqtec, chosen_meta$sample.origin, sep = "\t")
  if (!identical(actual_pairs, expected_pairs)) {
    stop("select_debug_samples: selected technology/state pairs are not exact.")
  }
  if (any(chosen_meta$n_cells < min_cells)) {
    stop("select_debug_samples: selected sample has fewer than min_cells cells.")
  }

  message("\n=== Selected Debug Samples (", length(chosen), " samples) ===")
  print(chosen_meta[, c("sample.ID", "seqtec", "sample.origin", "Site", "n_cells")],
        row.names = FALSE)
  chosen
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

# Keep minimal obs columns (including both configured cell-type roles).
cols_keep <- intersect(
  c("sample.ID", "sample.origin", "patient.ID", "iCMS",
    "dataset", "Gender", "Site", "cell.type", "cell.type_new", "seqtec"),
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
h5ad_tmp <- tempfile(pattern = paste0(".", basename(h5ad_path), ".tmp-"),
                     tmpdir = dirname(h5ad_path), fileext = ".h5ad")
write_h5ad(seurat, h5ad_tmp, mode = "w")
if (!file.exists(h5ad_tmp) || file.info(h5ad_tmp)$size <= 0) {
  unlink(h5ad_tmp)
  stop("Atomic Joanito debug h5ad write produced an empty file: ", h5ad_tmp)
}
if (!file.rename(h5ad_tmp, h5ad_path)) {
  unlink(h5ad_tmp)
  stop("Could not atomically install Joanito debug h5ad: ", h5ad_path)
}
message("Saved: ", h5ad_path)

# -----------------------------------------------------------------------------
# Summary for verification
# -----------------------------------------------------------------------------
summary_table <- unique(meta[, intersect(c("sample.ID", "sample.origin", "seqtec", "Site"), colnames(meta))])
summary_table$n_cells <- as.integer(table(meta$sample.ID)[summary_table$sample.ID])
print(summary_table, row.names = FALSE)

message("Debug dataset created: ", ncol(seurat), " cells across ",
        length(unique(meta$sample.ID)), " samples.")
