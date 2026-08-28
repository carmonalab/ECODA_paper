#!/usr/bin/env Rscript
# Evaluate only the pure selector function; do not load Seurat/anndataR or raw data.
expressions <- parse("src/2_dataset_specific_preprocessing/1.3.1_prepare_joanito.R")
selector_assignment <- Filter(function(expr) {
  is.call(expr) && identical(as.character(expr[[1L]]), "<-") &&
    identical(as.character(expr[[2L]]), "select_debug_samples")
}, as.list(expressions))
stopifnot(length(selector_assignment) == 1L)
eval(selector_assignment[[1L]], envir = .GlobalEnv)

pairs <- data.frame(
  seqtec = c("3' seq", "3' seq", "5' seq", "5' seq", "5' seq"),
  sample.origin = c("Normal", "Tumor", "Normal", "Tumor", "LymphNode"),
  Site = c("Ascending colon", "Sigmoid colon", "Rectum", "Upper rectum", "Caecum"),
  stringsAsFactors = FALSE
)
meta <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  data.frame(
    sample.ID = rep(paste0("sample_", i), 500),
    sample.origin = rep(pairs$sample.origin[i], 500),
    seqtec = rep(pairs$seqtec[i], 500),
    Site = rep(pairs$Site[i], 500),
    stringsAsFactors = FALSE
  )
}))
meta$dataset <- "OTHER"
meta$cell.type <- "T"
meta$iCMS <- "Normal"
chosen <- select_debug_samples(meta, min_cells = 500)
stopifnot(length(chosen) == 5L, length(unique(chosen)) == 5L)
chosen_rows <- meta[match(chosen, meta$sample.ID), c("seqtec", "sample.origin")]
stopifnot(identical(
  paste(chosen_rows$seqtec, chosen_rows$sample.origin, sep = "\t"),
  paste(pairs$seqtec, pairs$sample.origin, sep = "\t")
))

derive_assignment <- Filter(function(expr) {
  is.call(expr) && identical(as.character(expr[[1L]]), "<-") &&
    identical(as.character(expr[[2L]]), "derive_cell_type_new")
}, as.list(expressions))
stopifnot(length(derive_assignment) == 1L)
eval(derive_assignment[[1L]], envir = .GlobalEnv)
stopifnot(identical(
  derive_cell_type_new(c("T_Normal", "B"), c("Normal", NA_character_)),
  c("T_Normal_Normal", "B")
))

expect_error <- function(expr, label) {
  failed <- FALSE
  tryCatch(force(expr), error = function(e) failed <<- TRUE)
  if (!failed) stop("expected selector failure: ", label)
}

over_count <- rbind(
  meta,
  data.frame(
    sample.ID = rep("sample_1", 100),
    sample.origin = rep("Normal", 100),
    seqtec = rep("3' seq", 100),
    Site = rep("Ascending colon", 100),
    dataset = rep("OTHER", 100),
    cell.type = rep("T", 100),
    iCMS = rep("Normal", 100),
    stringsAsFactors = FALSE
  )
)
stopifnot("sample_1" %in% select_debug_samples(over_count, min_cells = 500))

expect_error(
  select_debug_samples(meta[meta$sample.ID != "sample_5", , drop = FALSE], 500),
  "missing technology/state pair"
)
under_count <- meta[-which(meta$sample.ID == "sample_1")[1], , drop = FALSE]
expect_error(select_debug_samples(under_count, 500), "sample below minimum count")

mixed <- meta
mixed$Site[which(mixed$sample.ID == "sample_1")[1]] <- "Rectum"
expect_error(select_debug_samples(mixed, 500), "mixed per-sample Site")

blank_id <- meta
blank_id$sample.ID[1] <- " "
expect_error(select_debug_samples(blank_id, 500), "blank sample ID")

duplicate_candidate <- meta
duplicate_candidate$sample.ID[which(duplicate_candidate$sample.ID == "sample_2")] <- "sample_1"
expect_error(select_debug_samples(duplicate_candidate, 500), "duplicate candidate ID")

retained <- meta[match(chosen, meta$sample.ID), c("sample.ID", "sample.origin", "seqtec", "Site")]
stopifnot(
  all(nzchar(retained$sample.ID)),
  all(nzchar(retained$sample.origin)),
  all(nzchar(retained$seqtec)),
  all(nzchar(retained$Site))
)
cat("Joanito five-sample selector: OK\n")
