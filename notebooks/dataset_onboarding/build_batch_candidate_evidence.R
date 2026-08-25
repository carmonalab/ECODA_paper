#!/usr/bin/env Rscript
# Build the uncorrected-pass technical-batch decision evidence.
#
# The script is deliberately strict about sample identity. It reads sample-level
# metadata from each uncorrected h5ad, then aligns every available method
# distance matrix before computing candidate summaries and PERMANOVA statistics.
# Biological labels are evaluation metadata only; they are never supplied to a
# method or correction model. Missing method artifacts are reported as
# unavailable, while malformed artifacts and sample-order mismatches fail closed.

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
repo_root <- normalizePath(
  file.path(dirname(sub("^--file=", "", file_arg)), "../..")
)
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x)) y else x
}
args <- commandArgs(trailingOnly = TRUE)
config_path <- parse_arg(args, "--config", file.path(repo_root, "datasets.json"))
analysis_root <- parse_arg(args, "--analysis-root")
input_root <- parse_arg(args, "--input-root")
output_dir <- parse_arg(args, "--output-dir", file.path(repo_root, "data", "batch_candidate_evidence"))
only_ds <- parse_arg(args, "--ds_name", NULL)
if (is.null(analysis_root) || is.null(input_root)) {
  stop("Usage: build_batch_candidate_evidence.R --analysis-root <scratch/batch_effect/uncorrected> --input-root <scratch> [--config <datasets.json>] [--output-dir <dir>] [--ds_name <dataset>]")
}

config <- fromJSON(config_path, simplifyVector = FALSE)
spec_path <- file.path(repo_root, "notebooks", "dataset_onboarding", "dataset_specs.py")
spec_module <- import_from_path("dataset_specs", path = dirname(spec_path), convert = TRUE)
specs <- spec_module$DATASET_SPECS
if (is.null(specs) || length(specs) == 0L) stop("No DATASET_SPECS loaded from ", spec_path)

dataset_names <- names(specs)
if (!is.null(only_ds)) dataset_names <- only_ds
missing_config <- setdiff(dataset_names, names(config))
if (length(missing_config) > 0L) stop("Dataset(s) absent from datasets.json: ", paste(missing_config, collapse = ", "))

uncorrected_view <- function(entry) {
  view <- entry$views$batch_effect_uncorrected
  if (is.null(view) || is.null(view$output_file_name)) {
    stop("Dataset lacks batch_effect_uncorrected output view")
  }
  view
}

effective_columns <- function(entry, view) {
  base <- entry$columns %||% list()
  override <- view$columns %||% list()
  modifyList(base, override)
}

read_sample_metadata <- function(ds, entry, spec) {
  view <- uncorrected_view(entry)
  cols <- effective_columns(entry, view)
  h5ad_path <- file.path(input_root, ds, "output", view$output_file_name)
  if (!file.exists(h5ad_path)) stop(ds, ": uncorrected h5ad not found: ", h5ad_path)

  ad <- import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(h5ad_path, backed = "r")
  on.exit(try(adata$file$close(), silent = TRUE), add = TRUE)
  obs <- py_to_r(adata$obs)
  obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  sample_col <- "Sample"
  if (!sample_col %in% colnames(obs)) stop(ds, ": standardized Sample column missing from ", h5ad_path)
  label_col <- cols$label
  if (is.null(label_col) || !label_col %in% colnames(obs)) stop(ds, ": biological label missing: ", label_col)

  sample_ids <- as.character(obs[[sample_col]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) stop(ds, ": missing standardized sample IDs")
  keep <- !duplicated(sample_ids)
  meta <- obs[keep, , drop = FALSE]
  rownames(meta) <- sample_ids[keep]
  meta[[sample_col]] <- sample_ids[keep]
  list(
    meta = meta,
    sample_col = sample_col,
    label_col = label_col,
    candidates = as.character(unlist(spec$batch_cols %||% character(0))),
    path = h5ad_path
  )
}

nmi_score <- function(x, y) {
  valid <- !is.na(x) & !is.na(y) & nzchar(as.character(x)) & nzchar(as.character(y))
  x <- as.character(x[valid]); y <- as.character(y[valid])
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  tab <- table(x, y)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy); py <- colSums(pxy)
  nz <- pxy > 0
  mi <- sum(pxy[nz] * log(pxy[nz] / outer(px, py)[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

candidate_summary <- function(meta, candidate, label_col) {
  if (!candidate %in% colnames(meta)) {
    return(list(
      present = FALSE, completeness = 0, levels = NA_integer_,
      samples_per_level = NA_character_, nmi_biology = NA_real_,
      constant = NA, sample_unique = NA, perfect_confounded = NA,
      warnings = paste0("missing candidate column ", candidate)
    ))
  }
  values <- as.character(meta[[candidate]])
  valid <- !is.na(values) & nzchar(values)
  n_samples <- nrow(meta)
  level_counts <- table(values[valid])
  n_levels <- length(level_counts)
  label <- as.character(meta[[label_col]])
  warnings <- character(0)
  constant <- n_levels < 2L
  sample_unique <- n_levels == sum(valid) && sum(valid) == n_samples
  perfect_confounded <- FALSE
  if (constant) warnings <- c(warnings, "constant candidate")
  if (sample_unique) warnings <- c(warnings, "sample-unique candidate")
  if (any(!valid)) warnings <- c(warnings, "incomplete candidate")
  if (candidate == label_col) warnings <- c(warnings, "candidate equals biological label")
  valid_pair <- valid & !is.na(label) & nzchar(label)
  if (sum(valid_pair) > 0L) {
    pair <- unique(paste(label[valid_pair], values[valid_pair], sep = "\r"))
    perfect_confounded <- length(pair) == length(unique(label[valid_pair])) &&
      length(pair) == length(unique(values[valid_pair]))
    if (perfect_confounded) warnings <- c(warnings, "perfectly confounded with biology")
  }
  list(
    present = TRUE,
    completeness = mean(valid),
    levels = n_levels,
    samples_per_level = paste(names(level_counts), as.integer(level_counts), sep = ":", collapse = ";"),
    nmi_biology = nmi_score(label, values),
    constant = constant,
    sample_unique = sample_unique,
    perfect_confounded = perfect_confounded,
    warnings = paste(warnings, collapse = "; ")
  )
}

read_feather_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- arrow::read_feather(path)
  if (ncol(df) < 2L) stop("Feather distance file has no index column: ", path)
  ids <- as.character(df[[ncol(df)]])
  mat <- as.matrix(df[, seq_len(ncol(df) - 1L), drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  colnames(mat) <- ids
  mat
}

read_r_bundle_matrix <- function(path, combo) {
  if (!file.exists(path)) return(NULL)
  bundle <- readRDS(path)
  result <- bundle[[combo]]
  if (is.null(result) || is.null(result$dist_mat)) return(NULL)
  mat <- as.matrix(result$dist_mat)
  if (is.null(rownames(mat))) stop("Distance bundle lacks sample rownames: ", path, " / ", combo)
  storage.mode(mat) <- "double"
  mat
}

load_method_matrices <- function(ds, result_root, entry) {
  stem <- paste0(ds, "_batch_effect_uncorrected")
  results_dir <- file.path(result_root, "results")
  embeddings_dir <- file.path(result_root, "embeddings")
  composition_file <- file.path(results_dir, paste0(stem, "_composition.rds"))
  composition <- if (file.exists(composition_file)) readRDS(composition_file) else list()
  methods <- list()
  composition_names <- c(
    "ECODA_authors_HR", "ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR",
    "ECODA_seuratres_2", "ECODA_authors_HR_NULL"
  )
  for (name in composition_names) {
    methods[[name]] <- if (!is.null(composition[[name]])) as.matrix(composition[[name]]$dist_mat) else NULL
  }
  pseudobulk_file <- file.path(results_dir, paste0(stem, "_pseudobulk.rds"))
  pseudobulk <- if (file.exists(pseudobulk_file)) readRDS(pseudobulk_file) else list()
  methods[["Pseudobulk_hvg2000"]] <- if (!is.null(pseudobulk[["Pseudobulk_hvg2000"]])) as.matrix(pseudobulk[["Pseudobulk_hvg2000"]]$dist_mat) else NULL
  gloscope_file <- file.path(results_dir, paste0(stem, "_gloscope.rds"))
  gloscope <- if (file.exists(gloscope_file)) readRDS(gloscope_file) else list()
  methods[["GloScope_hvg2000_pcadims30"]] <- if (!is.null(gloscope[["GloScope_hvg2000_pcadims30"]])) as.matrix(gloscope[["GloScope_hvg2000_pcadims30"]]$dist_mat) else NULL

  feather_specs <- c(
    PILOT_hvg2000_highres = "pilot",
    `PILOT-GM-VAE_hvg2000_highres` = "pilotgm",
    MrVI_hvg2000 = "mrvi",
    QOT_hvg2000_highres = "qot"
  )
  for (name in names(feather_specs)) {
    suffix <- feather_specs[[name]]
    path <- file.path(
      embeddings_dir,
      paste0(stem, "_hvg2000_highres_", suffix, "_dists.feather")
    )
    methods[[name]] <- read_feather_matrix(path)
  }
  methods
}

check_matrix_alignment <- function(mat, sample_ids, method, ds) {
  if (is.null(mat)) return(FALSE)
  ids <- rownames(mat)
  if (is.null(ids) || anyDuplicated(ids) || !identical(as.character(ids), as.character(sample_ids))) {
    stop(ds, "/", method, ": strict sample-order mismatch")
  }
  if (nrow(mat) != ncol(mat) || any(!is.finite(mat))) stop(ds, "/", method, ": invalid distance matrix")
  TRUE
}

permanova_stats <- function(mat, meta, label_col, candidate) {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("vegan is required for candidate PERMANOVA evidence")
  if (length(unique(meta[[candidate]])) < 2L) {
    return(c(marginal_r2 = NA_real_, marginal_p = NA_real_, joint_r2 = NA_real_, joint_p = NA_real_))
  }
  data <- data.frame(
    bio = factor(as.character(meta[[label_col]])),
    candidate = factor(as.character(meta[[candidate]])),
    row.names = rownames(meta)
  )
  dist_obj <- stats::as.dist(mat)
  marginal <- vegan::adonis2(
    dist_obj ~ candidate, data = data, permutations = N_PERMUTATIONS, by = "margin"
  )
  joint <- vegan::adonis2(
    dist_obj ~ bio + candidate, data = data, permutations = N_PERMUTATIONS, by = "margin"
  )
  joint_row <- which(rownames(joint) == "candidate")
  if (length(joint_row) != 1L) stop("Unable to locate candidate term in joint PERMANOVA")
  c(
    marginal_r2 = as.numeric(marginal$R2[1]),
    marginal_p = as.numeric(marginal$`Pr(>F)`[1]),
    joint_r2 = as.numeric(joint$R2[joint_row]),
    joint_p = as.numeric(joint$`Pr(>F)`[joint_row])
  )
}

build_dataset_evidence <- function(ds, entry, spec) {
  sample_info <- read_sample_metadata(ds, entry, spec)
  meta <- sample_info$meta
  candidate_names <- sample_info$candidates
  method_mats <- load_method_matrices(ds, analysis_root, entry)
  sample_ids <- rownames(meta)
  rows <- list(); row_idx <- 0L
  for (method in names(method_mats)) {
    mat <- method_mats[[method]]
    available <- !is.null(mat)
    if (available) check_matrix_alignment(mat, sample_ids, method, ds)
    for (candidate in candidate_names) {
      row_idx <- row_idx + 1L
      summary <- candidate_summary(meta, candidate, sample_info$label_col)
      stats <- c(marginal_r2 = NA_real_, marginal_p = NA_real_, joint_r2 = NA_real_, joint_p = NA_real_)
      if (available && isTRUE(summary$present) && !isTRUE(summary$constant) && !isTRUE(summary$sample_unique)) {
        stats <- permanova_stats(mat, meta, sample_info$label_col, candidate)
      }
      rows[[row_idx]] <- data.frame(
        dataset = ds,
        method = method,
        method_available = available,
        candidate = candidate,
        completeness = summary$completeness,
        levels = summary$levels,
        samples_per_level = summary$samples_per_level,
        nmi_biology = summary$nmi_biology,
        constant_candidate = summary$constant,
        sample_unique_candidate = summary$sample_unique,
        perfect_confounded = summary$perfect_confounded,
        marginal_r2 = stats[["marginal_r2"]],
        marginal_p = stats[["marginal_p"]],
        joint_r2 = stats[["joint_r2"]],
        joint_p = stats[["joint_p"]],
        warnings = paste(c(summary$warnings, if (!available) "method artifact unavailable" else character(0)), collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  result$marginal_p_holm <- NA_real_
  result$joint_p_holm <- NA_real_
  for (method in unique(result$method)) {
    idx <- result$method == method
    result$marginal_p_holm[idx] <- p.adjust(result$marginal_p[idx], method = "holm")
    result$joint_p_holm[idx] <- p.adjust(result$joint_p[idx], method = "holm")
  }
  result
}

mkdir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)
mkdir(output_dir)
all_results <- list()
for (ds in dataset_names) {
  message("Building uncorrected candidate evidence for ", ds)
  result <- build_dataset_evidence(ds, config[[ds]], specs[[ds]])
  all_results[[ds]] <- result
  write.csv(result, file.path(output_dir, paste0(ds, "_batch_candidate_evidence.csv")), row.names = FALSE, na = "")
}
combined <- do.call(rbind, all_results)
write.csv(combined, file.path(output_dir, "batch_candidate_review.csv"), row.names = FALSE, na = "")
message("Wrote ", length(all_results), " cohort evidence CSVs and batch_candidate_review.csv to ", output_dir)
