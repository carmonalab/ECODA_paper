# ==============================================================================
# onboarding_metrics.R -- quick batch-effect separation check for new datasets
# ==============================================================================
# Cell-level LISI-based separation (repo's own scoring_metrics.R::calc_lisi)
# on the UNINTEGRATED PCA embedding, per cell type and for the whole sample.
#
# Sourced/run with the pixi DEFAULT env ONLY (`pixi run -e default Rscript`).
# It sources src/utils/scoring_metrics.R directly -- NOT the notebook loader
# (load_all_functions.R / imports.R serve only benchmark_analysis.rmd and
# batch_effect_analysis.rmd).
#
# Input : feather written by notebooks/dataset_onboarding/onboarding_utils.py
#         `write_metrics_input`: one row per cell with columns
#           cell_index, <ct_col>, <bio_col>, <batch_cols>..., PC_1..PC_n
# Output: CSV (CT x label separation table) + JSON (verdict + skip/confound
#         bookkeeping), consumed by the dataset_check_*.qmd notebooks.
#
# Guards (per the Phase-5 plan):
#   * per-CT subsample cap (default 2000 cells, fixed seed)
#   * CTs with < min_cells (default 50) are SKIPPED and listed
#   * CTs with >90% cells from a single batch or bio group are tagged
#     "confounded/uninformative" (no score -- avoids false-positive warnings,
#     e.g. a rare CT present in one donor)
#   * single unique label within a CT -> calc_lisi returns NA (existing
#     behavior; recorded as "single_label")
#
# Usage:
#   pixi run -e default Rscript --vanilla onboarding_metrics.R \
#     --input <feather> --ct-col <col> --bio-col <col> \
#     --batch-cols <c1,c2> --out-csv <path> --out-json <path> \
#     [--seed 0] [--per-ct-cap 2000] [--min-cells 50] [--confounded-frac 0.9]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) return(default)
  args[i + 1]
}

this_file <- sub("--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
repo_root <- normalizePath(file.path(dirname(this_file), "..", ".."))
source(file.path(repo_root, "src", "utils", "scoring_metrics.R"))

input_path    <- parse_arg("--input")
ct_col        <- parse_arg("--ct-col")
bio_col       <- parse_arg("--bio-col")
batch_cols_raw <- parse_arg("--batch-cols")
batch_cols    <- if (is.null(batch_cols_raw) || batch_cols_raw == "") character(0) else strsplit(batch_cols_raw, ",", fixed = TRUE)[[1]]
out_csv       <- parse_arg("--out-csv")
out_json      <- parse_arg("--out-json")
seed          <- as.integer(parse_arg("--seed", "0"))
per_ct_cap    <- as.integer(parse_arg("--per-ct-cap", "2000"))
min_cells     <- as.integer(parse_arg("--min-cells", "50"))
confounded_frac <- as.numeric(parse_arg("--confounded-frac", "0.9"))

stopifnot(
  !is.null(input_path), file.exists(input_path),
  !is.null(ct_col), !is.null(bio_col), !is.null(out_csv), !is.null(out_json)
)
stopifnot(requireNamespace("arrow", quietly = TRUE))

set.seed(seed)
dat <- arrow::read_feather(input_path)
dat <- as.data.frame(dat)
pc_cols <- grep("^PC_", colnames(dat), value = TRUE)
if (length(pc_cols) == 0) stop("No PC_* columns found in ", input_path)
if (!ct_col %in% colnames(dat)) stop("ct column '", ct_col, "' not in feather")
if (!bio_col %in% colnames(dat)) stop("bio column '", bio_col, "' not in feather")
missing_batch <- setdiff(batch_cols, colnames(dat))
if (length(missing_batch) > 0) stop("batch column(s) missing: ", paste(missing_batch, collapse = ", "))

label_cols <- c(bio_col, batch_cols)

# Treat literal "<NA>" (python-side missing marker) as missing and drop.
for (lbl in label_cols) dat[[lbl]][dat[[lbl]] == "<NA>"] <- NA_character_
dat <- dat[!is.na(dat[[ct_col]]), , drop = FALSE]
dat[[ct_col]] <- as.character(dat[[ct_col]])
if (nrow(dat) == 0) stop("no cells left after NA ct-column drop")

all_cts <- unique(dat[[ct_col]])

sep_lisi <- function(lbl, cells_df) {
  n <- nrow(cells_df)
  if (n < min_cells) return(list(status = "skipped_small", score = NA_real_))
  vc <- table(cells_df[[lbl]], useNA = "ifany")
  vc <- vc[!is.na(names(vc))]
  dom <- max(vc) / sum(vc)
  if (dom > confounded_frac) {
    return(list(
      status = "confounded",
      score = NA_real_,
      note = sprintf("%.0f%% of cells in a single %s level", 100 * dom, lbl)
    ))
  }
  if (length(unique(droplevels(factor(cells_df[[lbl]])))) <= 1) {
    return(list(status = "single_label", score = NA_real_))
  }
  idx <- seq_len(min(per_ct_cap, n))
  cells_df <- cells_df[idx, , drop = FALSE]
  features <- as.matrix(cells_df[, pc_cols, drop = FALSE])
  rownames(features) <- paste0("cell_", seq_len(nrow(features)))
  score <- tryCatch(
    calc_lisi(labels = as.character(cells_df[[lbl]]), features = features),
    error = function(e) NA_real_
  )
  list(status = "ok", score = score)
}

# rows: each CT + the whole-sample pseudo-CT
cts_to_run <- unique(c(all_cts, "<ALL>"))
res_rows <- list()
verdict_notes <- list()

for (ct in cts_to_run) {
  if (ct == "<ALL>") {
    cells_df <- dat
  } else {
    cells_df <- dat[dat[[ct_col]] == ct, , drop = FALSE]
  }
  n <- nrow(cells_df)
  new_row <- function(ct, n, bio_score, bio_status) {
    row <- data.frame(
      cell_type = ct,
      n_cells = n,
      bio_separation = bio_score,
      bio_status = bio_status,
      stringsAsFactors = FALSE
    )
    for (bcol in batch_cols) {
      row[[paste0("batch_", bcol, "_separation")]] <- NA_real_
      row[[paste0("batch_", bcol, "_status")]] <- "skipped_small"
    }
    row
  }
  if (n < min_cells) {
    res_rows[[ct]] <- new_row(ct, n, NA_real_, "skipped_small")
    verdict_notes[[paste0(ct, "|size")]] <- sprintf("%s: skipped (n=%d < %d)", ct, n, min_cells)
    next
  }
  # balance: fixed-seed per-CT subsample so kNN/thisutils stays cheap
  set.seed(seed)
  keep <- sort(sample.int(n, size = min(per_ct_cap, n)))
  cells_df <- cells_df[keep, , drop = FALSE]

  bio_res <- sep_lisi(bio_col, cells_df)
  row <- new_row(ct, n, bio_res$score, bio_res$status)
  for (bcol in batch_cols) {
    b_res <- sep_lisi(bcol, cells_df)
    row[[paste0("batch_", bcol, "_separation")]] <- b_res$score
    row[[paste0("batch_", bcol, "_status")]] <- b_res$status
    if (!is.null(b_res$note)) verdict_notes[[paste0(ct, "|", bcol)]] <- sprintf(
      "%s / %s: %s", ct, bcol, b_res$note
    )
  }
  res_rows[[ct]] <- row
}

sep_table <- do.call(rbind, res_rows)
rownames(sep_table) <- NULL
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(sep_table, out_csv, row.names = FALSE)

# Verdict: batch effect signal = CT(s) with a batch separation score that is
# both "ok" and above threshold, excluding the whole-sample row.
batch_cols_table <- grep("^batch_.*_separation$", colnames(sep_table), value = TRUE)
verdict_scores <- list()
for (bcol in batch_cols) {
  sc_col <- paste0("batch_", bcol, "_separation")
  st_col <- paste0("batch_", bcol, "_status")
  if (!sc_col %in% colnames(sep_table)) next
  ok <- !is.na(sep_table[[sc_col]]) & sep_table[[st_col]] == "ok" & sep_table$cell_type != "<ALL>"
  above <- sep_table$cell_type[ok & sep_table[[sc_col]] > 0.5]
  verdict_scores[[bcol]] <- list(
    above_threshold_cts = above,
    max_per_ct = setNames(round(sep_table[[sc_col]][ok], 3), sep_table$cell_type[ok])
  )
}

verdict_lines <- c(
  sprintf(
    "LISI separation (1 = perfect separation, 0 = perfect mixing) on the unintegrated PCA embedding (%d PCs), %d cells, seed %d; per-CT cap %d, min cells %d.",
    length(pc_cols), nrow(dat), seed, per_ct_cap, min_cells
  ),
  sprintf(
    "Biological label '%s': overall separation = %s (low per-CT bio separation is expected -- cells can be identical across conditions; composition may still differ).",
    bio_col,
    if (is.na(sep_table$bio_separation[sep_table$cell_type == "<ALL>"])) {
      paste0("NA/", sep_table$bio_status[sep_table$cell_type == "<ALL>"])
    } else {
      sprintf("%.3f", sep_table$bio_separation[sep_table$cell_type == "<ALL>"])
    }
  )
)
for (bcol in batch_cols) {
  vs <- verdict_scores[[bcol]]
  if (length(vs$above_threshold_cts) == 0) {
    verdict_lines <- c(verdict_lines, sprintf(
      "Batch candidate '%s': no cell type with separation > 0.5 -- no clear batch effect signal within the checked cell types.", bcol
    ))
  } else {
    verdict_lines <- c(verdict_lines, sprintf(
      "Batch candidate '%s': separation > 0.5 within cell type(s): %s -- batch effect signal; interpret together with the confounding crosstab (collinear designs make bio vs batch indistinguishable).",
      bcol, paste(vs$above_threshold_cts, collapse = ", ")
    ))
  }
}
if (length(verdict_notes) > 0) {
  verdict_lines <- c(verdict_lines, paste0("Note: ", paste(unlist(verdict_notes), collapse = "; ")))
}

out_json_obj <- list(
  seed = seed, per_ct_cap = per_ct_cap, min_cells = min_cells,
  confounded_frac = confounded_frac,
  n_cells_input = nrow(dat),
  n_pcs = length(pc_cols),
  verdict = verdict_lines,
  batch_verdicts = verdict_scores,
  cts = setNames(
    lapply(all_cts, function(ct) {
      r <- sep_table[sep_table$cell_type == ct, , drop = FALSE]
      as.list(r)
    }),
    all_cts
  )
)
write(jsonlite::toJSON(out_json_obj, auto_unbox = TRUE, digits = 4), out_json)
cat(paste(verdict_lines, collapse = "\n"), "\n")
cat("Wrote separation table:", out_csv, "\n")
cat("Wrote JSON summary:", out_json, "\n")