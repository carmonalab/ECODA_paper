#!/usr/bin/env Rscript
# Build uncorrected-pass technical-batch decision evidence for the exact
# twelve-cohort onboarding selection. Biological labels are evaluation fields;
# they are never supplied to preprocessing or benchmark methods.

suppressPackageStartupMessages({
  library(jsonlite)
  library(reticulate)
  library(arrow)
})

N_PERMUTATIONS <- 999L

parse_arg <- function(args, flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1L]]
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Cannot determine repository root")
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

args <- commandArgs(trailingOnly = TRUE)
config_path <- parse_arg(args, "--config", file.path(repo_root, "datasets.json"))
analysis_root <- parse_arg(args, "--analysis-root")
input_root <- parse_arg(args, "--input-root")
output_dir <- parse_arg(args, "--output-dir", file.path(repo_root, "data", "batch_candidate_evidence"))
selection_file <- parse_arg(args, "--selection-file")
if (is.null(analysis_root) || is.null(input_root) || is.null(selection_file)) {
  stop("Usage: build_batch_candidate_evidence.R --selection-file <TSV> --analysis-root <scratch/batch_effect/uncorrected> --input-root <scratch> [--config <datasets.json>] [--output-dir <dir>]")
}
if (!is.null(parse_arg(args, "--ds_name", NULL))) {
  stop("--ds_name is not supported: evidence scope is exactly twelve cohorts")
}

config <- fromJSON(config_path, simplifyVector = FALSE)
spec_path <- file.path(repo_root, "notebooks", "dataset_onboarding", "dataset_specs.py")
spec_module <- import_from_path("dataset_specs", path = dirname(spec_path), convert = TRUE)
BATCH_EFFECT_DATASET_ORDER <- as.character(spec_module$BATCH_EFFECT_DATASET_ORDER)
BATCH_EFFECT_SPECS <- spec_module$BATCH_EFFECT_SPECS
if (length(BATCH_EFFECT_DATASET_ORDER) != 12L ||
    length(unique(BATCH_EFFECT_DATASET_ORDER)) != 12L ||
    !setequal(names(BATCH_EFFECT_SPECS), BATCH_EFFECT_DATASET_ORDER)) {
  stop("BATCH_EFFECT_SPECS must contain exactly the twelve ordered datasets")
}
for (ds in BATCH_EFFECT_DATASET_ORDER) {
  candidates <- as.character(unlist(BATCH_EFFECT_SPECS[[ds]] %||% character(0)))
  if (length(candidates) == 0L || any(!nzchar(candidates))) {
    stop(ds, ": BATCH_EFFECT_SPECS candidate list is empty")
  }
}

read_exact_selection <- function(path) {
  validate_source_artifact(path, "selection")
  if (!file.exists(path) || file.info(path)$size <= 0) stop("selection file is missing or empty: ", path)
  lines <- readLines(path, warn = FALSE)
  if (length(lines) != 12L || any(!nzchar(lines))) stop("selection file must contain exactly twelve nonblank rows")
  parsed <- lapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]])
  if (any(vapply(parsed, length, integer(1)) != 2L)) stop("selection file must be headerless DATASET<TAB>VIEW")
  ds <- vapply(parsed, `[[`, character(1), 1L)
  view <- vapply(parsed, `[[`, character(1), 2L)
  if (!identical(ds, BATCH_EFFECT_DATASET_ORDER)) stop("selection dataset order/set does not match BATCH_EFFECT_DATASET_ORDER")
  if (any(view != "batch_effect_uncorrected")) stop("selection file must use only batch_effect_uncorrected")
  data.frame(dataset = ds, view = view, stringsAsFactors = FALSE)
}

.sidecar_fields <- function(path) {
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || file.info(sidecar)$size <= 0) stop("missing checksum sidecar: ", sidecar)
  lines <- readLines(sidecar, warn = FALSE)
  fields <- list()
  for (line in lines) {
    parts <- strsplit(line, "=", fixed = TRUE)[[1]]
    if (length(parts) >= 2L) fields[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
  }
  fields
}

validate_source_artifact <- function(path, kind = "artifact") {
  if (!file.exists(path) || file.info(path)$size <= 0) stop(kind, " is missing or empty: ", path)
  fields <- .sidecar_fields(path)
  actual_md5 <- unname(tools::md5sum(path))
  actual_size <- as.character(file.info(path)$size)
  if (!identical(tolower(as.character(fields$MD5 %||% "")), tolower(actual_md5))) stop(kind, " MD5 mismatch: ", path)
  if (!identical(as.character(fields$SIZE %||% ""), actual_size)) stop(kind, " SIZE mismatch: ", path)
  if (!identical(as.character(fields$PATH %||% ""), path)) stop(kind, " PATH mismatch: ", path)
  invisible(list(path = path, md5 = actual_md5, size = as.numeric(actual_size), kind = kind))
}

artifact_checksum_ok <- function(path) {
  isTRUE(tryCatch({ validate_source_artifact(path); TRUE }, error = function(e) FALSE))
}

read_checked_rds <- function(path, kind = "RDS") {
  validate_source_artifact(path, kind)
  tryCatch(readRDS(path), error = function(e) stop(kind, " is malformed: ", path, " (", conditionMessage(e), ")"))
}

read_checked_feather <- function(path, kind = "Feather") {
  validate_source_artifact(path, kind)
  tryCatch(arrow::read_feather(path), error = function(e) stop(kind, " is malformed: ", path, " (", conditionMessage(e), ")"))
}
selection <- read_exact_selection(selection_file)

effective_columns <- function(entry, view) {
  modifyList(entry$columns %||% list(), view$columns %||% list())
}

uncorrected_view <- function(entry, ds) {
  view <- entry$views$batch_effect_uncorrected
  if (is.null(view) || is.null(view$output_file_name) || !identical(view$input_file_name %||% "", entry$file_names %||% view$input_file_name)) {
    stop(ds, ": malformed canonical batch_effect_uncorrected registry view")
  }
  view
}

read_sample_metadata <- function(ds, entry) {
  view <- uncorrected_view(entry, ds)
  cols <- effective_columns(entry, view)
  h5ad_path <- file.path(input_root, ds, "output", view$output_file_name)
  validate_source_artifact(h5ad_path, "uncorrected h5ad")

  ad <- import("anndata", convert = FALSE)
  adata <- tryCatch(ad$read_h5ad(h5ad_path, backed = "r"), error = function(e) stop(ds, ": malformed h5ad: ", conditionMessage(e)))
  on.exit(try(adata$file$close(), silent = TRUE), add = TRUE)
  shape <- as.integer(py_to_r(adata$shape))
  if (length(shape) < 2L || any(shape[seq_len(2L)] <= 0L)) stop(ds, ": uncorrected h5ad is empty")
  obs <- as.data.frame(py_to_r(adata$obs), stringsAsFactors = FALSE)
  if (!"Sample" %in% colnames(obs)) stop(ds, ": standardized Sample column missing")
  label_col <- as.character(cols$label %||% "")
  if (!nzchar(label_col) || !label_col %in% colnames(obs)) stop(ds, ": biological label missing: ", label_col)
  label_values <- as.character(obs[[label_col]])
  if (any(is.na(label_values)) || any(!nzchar(trimws(label_values)))) {
    stop(ds, ": biological label contains missing/blank values")
  }
  sample_ids <- as.character(obs$Sample)
  if (any(is.na(sample_ids)) || any(!nzchar(trimws(sample_ids)))) stop(ds, ": blank standardized sample IDs")
  keep <- !duplicated(sample_ids)
  meta <- obs[keep, , drop = FALSE]
  rownames(meta) <- sample_ids[keep]
  meta$Sample <- sample_ids[keep]
  list(
    meta = meta,
    sample_ids = sample_ids[keep],
    sample_col = "Sample",
    label_col = label_col,
    candidates = as.character(unlist(BATCH_EFFECT_SPECS[[ds]])),
    path = h5ad_path
  )
}

nmi_score <- function(x, y) {
  valid <- !is.na(x) & !is.na(y) & nzchar(trimws(as.character(x))) & nzchar(trimws(as.character(y)))
  x <- as.character(x[valid]); y <- as.character(y[valid])
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  tab <- table(x, y); pxy <- tab / sum(tab); px <- rowSums(pxy); py <- colSums(pxy)
  nz <- pxy > 0; mi <- sum(pxy[nz] * log(pxy[nz] / outer(px, py)[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0])); hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

candidate_summary <- function(meta, candidate, label_col) {
  if (!candidate %in% colnames(meta)) {
    return(list(present = FALSE, completeness = 0, levels = 0L, samples_per_level = "",
                nmi_biology = NA_real_, constant = NA, sample_unique = NA,
                perfect_confounded = NA, warnings = paste0("missing candidate column ", candidate)))
  }
  values <- as.character(meta[[candidate]]); valid <- !is.na(values) & nzchar(trimws(values))
  label <- as.character(meta[[label_col]]); n_samples <- nrow(meta)
  level_counts <- table(values[valid]); n_levels <- length(level_counts)
  warnings <- character(0); constant <- n_levels < 2L
  sample_unique <- n_levels == sum(valid) && sum(valid) == n_samples
  if (constant) warnings <- c(warnings, "constant candidate")
  if (sample_unique) warnings <- c(warnings, "sample-unique candidate")
  if (any(!valid)) warnings <- c(warnings, "incomplete candidate")
  valid_pair <- valid & !is.na(label) & nzchar(trimws(label))
  perfect_confounded <- FALSE
  if (sum(valid_pair) > 0L) {
    tab <- table(label[valid_pair], values[valid_pair])
    perfect_confounded <- all(rowSums(tab > 0) == 1L) && all(colSums(tab > 0) == 1L)
    if (perfect_confounded) warnings <- c(warnings, "perfectly confounded with biology")
  }
  list(present = TRUE, completeness = mean(valid), levels = n_levels,
       samples_per_level = paste(names(level_counts), as.integer(level_counts), sep = ":", collapse = ";"),
       nmi_biology = nmi_score(label, values), constant = constant,
       sample_unique = sample_unique, perfect_confounded = perfect_confounded,
       warnings = paste(warnings, collapse = "; "))
}

read_feather_matrix <- function(path) {
  df <- read_checked_feather(path, "Feather distance artifact")
  if (nrow(df) <= 0L || ncol(df) < 2L) stop("Feather distance file is empty or lacks an index column: ", path)
  ids <- as.character(df[[ncol(df)]])
  matrix_names <- colnames(df)[seq_len(ncol(df) - 1L)]
  if (any(is.na(ids)) || any(!nzchar(trimws(ids))) || anyDuplicated(ids) ||
      !identical(as.character(matrix_names), ids)) {
    stop("Feather distance sample-order mismatch, blank, or duplicated IDs: ", path)
  }
  mat <- as.matrix(df[, seq_len(ncol(df) - 1L), drop = FALSE]); storage.mode(mat) <- "double"
  if (nrow(mat) != ncol(mat) || nrow(mat) != length(ids) || any(!is.finite(mat))) stop("Feather distance matrix is non-square/non-finite: ", path)
  rownames(mat) <- ids; colnames(mat) <- ids
  mat
}

read_r_bundle_matrix <- function(path, combo) {
  bundle <- read_checked_rds(path, "RDS distance bundle")
  result <- bundle[[combo]]
  if (is.null(result) || is.null(result$dist_mat)) stop("RDS bundle lacks applicable matrix ", combo, ": ", path)
  mat <- as.matrix(result$dist_mat); storage.mode(mat) <- "double"
  if (nrow(mat) != ncol(mat) || is.null(rownames(mat)) || is.null(colnames(mat)) ||
      anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat)) || any(!is.finite(mat))) {
    stop("RDS distance matrix is malformed: ", path, " / ", combo)
  }
  if (!identical(as.character(rownames(mat)), as.character(colnames(mat)))) stop("RDS matrix row/column IDs differ: ", path, " / ", combo)
  mat
}

unavailable_annotation_combo <- function(ds, combo, entry) {
  excluded <- as.character(entry$not_suitable_for_auto_annotation %||% character(0))
  has_exclusion <- all(c("hitme", "scatomic") %in% excluded)
  if (has_exclusion && !ds %in% c("Alzheimer", "Diabetes", "Parkinson")) {
    stop(ds, ": unexpected dual-annotation exclusion in batch evidence scope")
  }
  has_exclusion && combo %in% c("ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR")
}
matrix_entry <- function(matrix = NULL, available = FALSE, applicable = TRUE, reason = "", artifact = "") {
  list(matrix = matrix, available = available, applicable = applicable, reason = reason, artifact = artifact)
}

load_method_matrices <- function(ds, result_root, entry) {
  stem <- paste0(ds, "_batch_effect_uncorrected")
  results_dir <- file.path(result_root, "results"); embeddings_dir <- file.path(result_root, "embeddings")
  composition_names <- c("ECODA_authors_HR", "ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR", "ECODA_seuratres_2", "ECODA_authors_HR_NULL")
  composition_file <- file.path(results_dir, paste0(stem, "_composition.rds"))
  if (!file.exists(composition_file)) stop(ds, ": missing applicable composition RDS: ", composition_file)
  composition <- read_checked_rds(composition_file, "composition RDS")
  methods <- list()
  for (name in composition_names) {
    if (unavailable_annotation_combo(ds, name, entry)) {
      methods[[name]] <- matrix_entry(reason = "not_suitable_for_auto_annotation", applicable = FALSE, artifact = composition_file)
    } else if (is.null(composition[[name]]) || is.null(composition[[name]]$dist_mat)) {
      stop(ds, ": applicable composition subresult is missing: ", name)
    } else {
      mat <- as.matrix(composition[[name]]$dist_mat); storage.mode(mat) <- "double"
      if (nrow(mat) != ncol(mat) || is.null(rownames(mat)) || is.null(colnames(mat)) ||
          anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat)) ||
          !identical(as.character(rownames(mat)), as.character(colnames(mat))) ||
          any(!is.finite(mat))) stop(ds, ": malformed composition matrix: ", name)
      methods[[name]] <- matrix_entry(mat, TRUE, TRUE, "", composition_file)
    }
  }

  bundle_specs <- list(
    Pseudobulk_hvg2000 = list(file = file.path(results_dir, paste0(stem, "_pseudobulk.rds")), combo = "Pseudobulk_hvg2000"),
    GloScope_hvg2000_pcadims30 = list(file = file.path(results_dir, paste0(stem, "_gloscope.rds")), combo = "GloScope_hvg2000_pcadims30")
  )
  for (name in names(bundle_specs)) {
    spec <- bundle_specs[[name]]; if (!file.exists(spec$file)) stop(ds, ": missing applicable ", name, " RDS: ", spec$file)
    bundle <- read_checked_rds(spec$file, paste0(name, " RDS")); result <- bundle[[spec$combo]]
    if (is.null(result) || is.null(result$dist_mat)) stop(ds, ": applicable RDS subresult is missing: ", name)
    mat <- as.matrix(result$dist_mat); storage.mode(mat) <- "double"
    if (nrow(mat) != ncol(mat) || is.null(rownames(mat)) || is.null(colnames(mat)) ||
        anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat)) ||
        !identical(as.character(rownames(mat)), as.character(colnames(mat))) ||
        any(!is.finite(mat))) stop(ds, ": malformed RDS distance matrix: ", name)
    methods[[name]] <- matrix_entry(mat, TRUE, TRUE, "", spec$file)
  }

  feather_specs <- c(
    PILOT_hvg2000_highres = "pilot",
    MrVI_hvg2000 = "mrvi",
    QOT_hvg2000_highres = "qot"
  )
  for (name in names(feather_specs)) {
    path <- file.path(embeddings_dir, paste0(stem, "_hvg2000_highres_", feather_specs[[name]], "_dists.feather"))
    if (!file.exists(path)) stop(ds, ": missing applicable Feather: ", path)
    methods[[name]] <- matrix_entry(read_feather_matrix(path), TRUE, TRUE, "", path)
  }
  methods
}

check_matrix_alignment <- function(mat, sample_ids, method, ds) {
  if (is.null(mat)) return(FALSE)
  ids <- rownames(mat)
  if (is.null(ids) || anyDuplicated(ids) || !identical(as.character(ids), as.character(sample_ids))) stop(ds, "/", method, ": strict sample-order mismatch")
  if (nrow(mat) != ncol(mat) || any(!is.finite(mat))) stop(ds, "/", method, ": invalid distance matrix")
  TRUE
}

permanova_stats <- function(mat, meta, label_col, candidate) {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("vegan is required for candidate PERMANOVA evidence")
  if (nrow(meta) < 3L || length(unique(meta[[candidate]])) < 2L) return(c(marginal_r2 = NA_real_, marginal_p = NA_real_, joint_r2 = NA_real_, joint_p = NA_real_))
  valid <- !is.na(meta[[candidate]]) & nzchar(trimws(as.character(meta[[candidate]]))) & !is.na(meta[[label_col]]) & nzchar(trimws(as.character(meta[[label_col]])))
  if (sum(valid) < 3L) return(c(marginal_r2 = NA_real_, marginal_p = NA_real_, joint_r2 = NA_real_, joint_p = NA_real_))
  data <- data.frame(bio = factor(as.character(meta[[label_col]][valid])), candidate = factor(as.character(meta[[candidate]][valid])), row.names = rownames(meta)[valid])
  dist_obj <- stats::as.dist(mat[valid, valid, drop = FALSE])
  marginal <- vegan::adonis2(dist_obj ~ candidate, data = data, permutations = N_PERMUTATIONS, by = "margin")
  joint <- vegan::adonis2(dist_obj ~ bio + candidate, data = data, permutations = N_PERMUTATIONS, by = "margin")
  joint_row <- which(rownames(joint) == "candidate"); if (length(joint_row) != 1L) stop("Unable to locate candidate term in joint PERMANOVA")
  c(marginal_r2 = as.numeric(marginal$R2[1]), marginal_p = as.numeric(marginal$`Pr(>F)`[1]), joint_r2 = as.numeric(joint$R2[joint_row]), joint_p = as.numeric(joint$`Pr(>F)`[joint_row]))
}

build_dataset_evidence <- function(ds, entry, spec = NULL) {
  sample_info <- read_sample_metadata(ds, entry); meta <- sample_info$meta
  method_entries <- load_method_matrices(ds, analysis_root, entry); sample_ids <- sample_info$sample_ids
  rows <- list(); row_idx <- 0L
  for (method in names(method_entries)) {
    entry_method <- method_entries[[method]]; available <- isTRUE(entry_method$available)
    if (available) check_matrix_alignment(entry_method$matrix, sample_ids, method, ds)
    for (candidate in sample_info$candidates) {
      row_idx <- row_idx + 1L; summary <- candidate_summary(meta, candidate, sample_info$label_col)
      stats <- c(marginal_r2 = NA_real_, marginal_p = NA_real_, joint_r2 = NA_real_, joint_p = NA_real_)
      if (available && isTRUE(summary$present) && !isTRUE(summary$constant) && !isTRUE(summary$sample_unique) && !isTRUE(summary$perfect_confounded)) stats <- permanova_stats(entry_method$matrix, meta, sample_info$label_col, candidate)
      warning_text <- summary$warnings
      if (!available) warning_text <- paste(c(warning_text, entry_method$reason), collapse = "; ")
      rows[[row_idx]] <- data.frame(
        dataset = ds, method = method, method_available = available,
        method_applicable = isTRUE(entry_method$applicable), method_reason = entry_method$reason,
        artifact = entry_method$artifact, candidate = candidate,
        present = summary$present, completeness = summary$completeness,
        levels = summary$levels, samples_per_level = summary$samples_per_level,
        nmi_biology = summary$nmi_biology, constant_candidate = summary$constant,
        sample_unique_candidate = summary$sample_unique,
        perfect_confounded = summary$perfect_confounded,
        marginal_r2 = stats[["marginal_r2"]], marginal_p = stats[["marginal_p"]],
        joint_r2 = stats[["joint_r2"]], joint_p = stats[["joint_p"]],
        warnings = warning_text, stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows); result$marginal_p_holm <- NA_real_; result$joint_p_holm <- NA_real_
  for (method in unique(result$method)) {
    idx <- result$method == method
    result$marginal_p_holm[idx] <- p.adjust(result$marginal_p[idx], method = "holm")
    result$joint_p_holm[idx] <- p.adjust(result$joint_p[idx], method = "holm")
  }
  result
}

write_artifact_sidecar <- function(path) {
  digest <- unname(tools::md5sum(path)); size <- file.info(path)$size
  sidecar <- paste0(path, ".md5"); tmp <- paste0(sidecar, ".tmp.", Sys.getpid())
  writeLines(c(paste0("MD5=", digest), paste0("SIZE=", size), paste0("PATH=", path)), tmp)
  if (!file.rename(tmp, sidecar)) {
    unlink(tmp)
    stop("cannot install artifact sidecar: ", sidecar)
  }
}

write_csv_atomic <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  write.csv(data, tmp, row.names = FALSE, na = "")
  if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("empty evidence CSV: ", tmp)
  if (!file.rename(tmp, path)) stop("cannot install evidence CSV: ", path)
  write_artifact_sidecar(path)
  validate_source_artifact(path, "evidence CSV")
}

for (ds in BATCH_EFFECT_DATASET_ORDER) {
  if (is.null(config[[ds]])) stop(ds, ": missing datasets.json entry")
  message("Building uncorrected candidate evidence for ", ds)
  result <- build_dataset_evidence(ds, config[[ds]], BATCH_EFFECT_SPECS[[ds]])
  output_path <- file.path(output_dir, paste0(ds, "_batch_candidate_evidence.csv"))
  write_csv_atomic(result, output_path)
}
combined <- do.call(rbind, lapply(
  BATCH_EFFECT_DATASET_ORDER,
  function(ds) read.csv(
    file.path(output_dir, paste0(ds, "_batch_candidate_evidence.csv")),
    stringsAsFactors = FALSE
  )
))
write_csv_atomic(combined, file.path(output_dir, "batch_candidate_review.csv"))
message("Wrote twelve cohort evidence CSVs and batch_candidate_review.csv to ", output_dir)
