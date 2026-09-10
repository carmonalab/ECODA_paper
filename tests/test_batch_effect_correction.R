#!/usr/bin/env Rscript

raw_args <- commandArgs(trailingOnly = FALSE)
script_arg <- raw_args[grepl("^--file=", raw_args)][1]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), ".."))
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
source(file.path(root, "src/utils/seurat_utils.R"))
source(file.path(root, "src/utils/pseudobulk.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_methods_r.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_pipeline.R"))

# CT timing method names use injective UTF-8 byte tokens.  In particular,
# punctuation and underscores must not collapse into one four-column log key.
ct_dot_token <- ct_timing_token("cell.type")
ct_underscore_token <- ct_timing_token("cell_type")
stopifnot(
  nzchar(ct_dot_token),
  nzchar(ct_underscore_token),
  !identical(ct_dot_token, ct_underscore_token),
  !identical(
    ct_shared_timing_method("cell.type"),
    ct_shared_timing_method("cell_type")
  ),
  grepl(
    "^prepare_pseudobulk_ct_shared_[A-Za-z0-9_-]+$",
    ct_shared_timing_method("cell.type")
  ),
  identical(
    ct_shared_timing_method_from_timing_id(
      paste(
        "run-1",
        paste0("Synthetic_ct_", ct_dot_token),
        "benchmark_analysis",
        "none",
        sep = ":"
      )
    ),
    ct_shared_timing_method("cell.type")
  )
)

# A schema-2 local/downstream bundle cannot replay a divergent exec_time.
divergent_exec_bundle <- list(
  shared_time_secs = 4,
  variant_time_secs = 2,
  shared_mem_GB = NA_real_,
  timing_id = "run-1:Synthetic:benchmark_analysis:none",
  timing_schema = 2L,
  exec_time = 3
)
divergent_exec_error <- tryCatch(
  validate_hpc_timing_bundle(divergent_exec_bundle),
  error = identity
)
stopifnot(inherits(divergent_exec_error, "error"))
divergent_exec_bundle$exec_time <- divergent_exec_bundle$variant_time_secs
stopifnot(isTRUE(
  validate_hpc_timing_bundle(divergent_exec_bundle)
))

# A configured author annotation may itself be a generated Leiden column.
# Preserve that source while adding the legacy alias consumed by composition.
leiden_source <- "leiden_res_5_batch_effect_uncorrected_hvg2000"
leiden_obs <- data.frame(
  Sample = paste0("s", 1:3),
  check = c("x", "y", "x"),
  stringsAsFactors = FALSE
)
leiden_obs[[leiden_source]] <- c("0", "1", "0")
leiden_mapped <- rename_leiden_cols(
  leiden_obs,
  view = "batch_effect_uncorrected",
  preserve_source = TRUE
)
stopifnot(
  leiden_source %in% colnames(leiden_mapped),
  "RNA_snn_res.5" %in% colnames(leiden_mapped),
  identical(leiden_mapped[[leiden_source]], leiden_mapped[["RNA_snn_res.5"]])
)
leiden_renamed <- rename_leiden_cols(
  leiden_obs,
  view = "batch_effect_uncorrected"
)
stopifnot(
  !leiden_source %in% colnames(leiden_renamed),
  "RNA_snn_res.5" %in% colnames(leiden_renamed)
)

# Batch composition emits only its three required keys. Legacy optional
# HiTME/scATOMIC bundles are not part of the batch execution contract.
original_process_coda_fig <- process_coda_fig
original_save_rds_atomic <- save_rds_atomic
original_artifact_checksum_ok <- artifact_checksum_ok
original_exec_time <- exec_time
original_peak_rss_gb <- peak_rss_gb
original_log_exec_row <- log_exec_row
process_coda_fig <- function(seurat, labels, ...) {
  feat <- matrix(
    1,
    nrow = length(labels),
    ncol = 1,
    dimnames = list(names(labels), "feature")
  )
  list(scores = 1, feat_mat = feat, dist_mat = dist(feat), labels = labels)
}
save_rds_atomic <- function(...) invisible(NULL)
artifact_checksum_ok <- function(...) FALSE
exec_time <- function(expr) {
  force(expr)
  0
}
peak_rss_gb <- function() 0
log_exec_row <- function(...) invisible(NULL)
batch_obs <- data.frame(
  Sample = rep(paste0("s", 1:3), each = 2),
  label = rep(c("A", "B", "A"), each = 2),
  author = rep(c("T", "B", "T"), each = 2),
  RNA_snn_res.2 = rep(c("0", "1", "0"), each = 2),
  layer2 = rep(c("T", "B", "T"), each = 2),
  scATOMIC_pred = rep(c("T", "B", "T"), each = 2)
)
batch_labels <- structure(
  factor(c("A", "B", "A")),
  names = paste0("s", 1:3)
)
batch_metadata <- data.frame(
  Sample = paste0("s", 1:3),
  label = c("A", "B", "A")
)
batch_composition <- run_composition_methods_hpc(
  batch_labels,
  batch_metadata,
  pca_emb = NULL,
  pb_hvg2000 = NULL,
  obs = batch_obs,
  label_col = "label",
  ct_col_high_res = "author",
  results_dir = tempfile("batch-composition-"),
  ds = "Synthetic",
  batch_mode = TRUE,
  result_stem = "Synthetic_batch_effect_uncorrected"
)
stopifnot(identical(
  names(batch_composition),
  c("ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2")
))
process_coda_fig <- original_process_coda_fig
save_rds_atomic <- original_save_rds_atomic
artifact_checksum_ok <- original_artifact_checksum_ok
exec_time <- original_exec_time
peak_rss_gb <- original_peak_rss_gb
log_exec_row <- original_log_exec_row

# Corrected CLR removes the fitted technical effect and preserves the CLR
# invariant exactly after row recentering.
set.seed(11)
n <- 12
feat <- matrix(
  rnorm(n * 4, sd = 0.2),
  nrow = n,
  dimnames = list(paste0("s", seq_len(n)), paste0("ct", seq_len(4)))
)
batch <- rep(c("A", "B"), each = n / 2)
feat[batch == "B", 1:2] <- feat[batch == "B", 1:2] + 4
feat <- feat - rowMeans(feat)
meta <- data.frame(Sample = rownames(feat), tech = factor(batch))
corrected <- correct_clr_batch_lmm(feat, meta, "tech")
stopifnot(max(abs(rowSums(corrected))) < 1e-8)
stopifnot(identical(dimnames(corrected), dimnames(feat)))

bad_order <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]
order_error <- tryCatch(
  correct_clr_batch_lmm(feat, bad_order, "tech"),
  error = identity
)
stopifnot(inherits(order_error, "error"))

missing_batch <- meta
missing_batch$tech[1] <- NA
missing_error <- tryCatch(
  correct_clr_batch_lmm(feat, missing_batch, "tech"),
  error = identity
)
stopifnot(inherits(missing_error, "error"))

# The canonical pseudobulk driver uses one raw H5AD aggregate and forwards the
# two batch modes to the direct DESeq2 fit without exposing biological labels.
captured <- new.env(parent = emptyenv())
fake_counts <- matrix(
  c(1, 2, 3, 4, 5, 6),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(c("g1", "g2"), paste0("s", 1:3))
)
fake_metadata <- data.frame(
  Sample = paste0("s", 1:3),
  tech = factor(c("A", "B", "A")),
  row.names = paste0("s", 1:3),
  stringsAsFactors = FALSE
)
prepare_pseudobulk_env <- environment(prepare_pseudobulks_hpc)
old_h5ad_aggregate <- get(
  "load_h5ad_sample_aggregate",
  envir = prepare_pseudobulk_env,
  inherits = TRUE
)
old_fit_pseudobulk <- get(
  "fit_pseudobulk_deseq2",
  envir = prepare_pseudobulk_env,
  inherits = TRUE
)
old_select_pseudobulk <- get(
  "select_pseudobulk_deseq2",
  envir = prepare_pseudobulk_env,
  inherits = TRUE
)
old_exec_time <- get("exec_time", envir = prepare_pseudobulk_env, inherits = TRUE)
old_peak_rss <- get("peak_rss_gb", envir = prepare_pseudobulk_env, inherits = TRUE)
assign(
  "load_h5ad_sample_aggregate",
  function(
    h5ad_path,
    sample_col = "Sample",
    metadata_columns = character(),
    chunk_size = 4096L,
    max_value = .Machine$integer.max
  ) {
    captured$aggregate <- list(
      path = h5ad_path,
      sample_col = sample_col,
      metadata_columns = metadata_columns,
      chunk_size = chunk_size,
      max_value = max_value
    )
    requested_columns <- unique(c(sample_col, metadata_columns))
    if (any(!requested_columns %in% colnames(fake_metadata))) {
      stop("fixture requested metadata column is unavailable")
    }
    returned_metadata <- fake_metadata[
      , requested_columns,
      drop = FALSE
    ]
    list(
      counts = fake_counts,
      sample_ids = colnames(fake_counts),
      gene_names = rownames(fake_counts),
      metadata = returned_metadata
    )
  },
  envir = prepare_pseudobulk_env
)
assign(
  "fit_pseudobulk_deseq2",
  function(
    counts,
    metadata,
    batch_col = NULL,
    blind = TRUE,
    correct_batch = FALSE
  ) {
    captured$fit <- list(
      counts = counts,
      metadata = metadata,
      batch_col = batch_col,
      blind = blind,
      correct_batch = correct_batch
    )
    variance_order <- rownames(counts)
    list(
      norm_matrix = counts,
      normalized_matrix = counts,
      variance_order = variance_order,
      variance_ordering = variance_order,
      row_variances = setNames(numeric(nrow(counts)), variance_order),
      counts = counts,
      metadata = metadata,
      batch_col = batch_col,
      blind = blind,
      correct_batch = correct_batch
    )
  },
  envir = prepare_pseudobulk_env
)
assign(
  "select_pseudobulk_deseq2",
  function(fit, n_hvg, black_list = "none") {
    captured$select <- list(n_hvg = n_hvg, black_list = black_list)
    fit$norm_matrix[seq_len(min(n_hvg, nrow(fit$norm_matrix))), , drop = FALSE]
  },
  envir = prepare_pseudobulk_env
)
assign("exec_time", function(expr) {
  force(expr)
  0
}, envir = prepare_pseudobulk_env)
assign("peak_rss_gb", function() 0, envir = prepare_pseudobulk_env)
fixture_h5ad <- tempfile("ecoda-pseudobulk-driver-")
writeBin(charToRaw("fixture-only"), fixture_h5ad)
invisible(prepare_pseudobulks_hpc(
  h5ad_path = fixture_h5ad,
  hvg_rank_genes = c("g1", "g2"),
  variants = "hvg2000",
  batch_col = "tech",
  blind = FALSE,
  correct_batch = TRUE,
  cache_stem = "fixture",
  view = "batch_effect_corrected",
  analysis_pass = "corrected",
  run_id = "fixture-run"
))
stopifnot(
  identical(captured$fit$counts, fake_counts),
  identical(captured$fit$metadata, fake_metadata),
  !"label" %in% colnames(captured$fit$metadata),
  identical(captured$fit$batch_col, "tech"),
  identical(captured$fit$blind, FALSE),
  identical(captured$fit$correct_batch, TRUE),
  identical(captured$aggregate$metadata_columns, c("Sample", "tech")),
  !"label" %in% captured$aggregate$metadata_columns
)
invisible(prepare_pseudobulks_hpc(
  h5ad_path = fixture_h5ad,
  hvg_rank_genes = c("g1", "g2"),
  variants = "hvg2000"
))
stopifnot(
  identical(captured$fit$counts, fake_counts),
  identical(
    captured$fit$metadata,
    fake_metadata[, "Sample", drop = FALSE]
  ),
  !"label" %in% colnames(captured$fit$metadata),
  identical(captured$fit$batch_col, NULL),
  identical(captured$fit$blind, TRUE),
  identical(captured$fit$correct_batch, FALSE),
  identical(captured$aggregate$metadata_columns, "Sample"),
  !"label" %in% captured$aggregate$metadata_columns
)
assign(
  "load_h5ad_sample_aggregate",
  old_h5ad_aggregate,
  envir = prepare_pseudobulk_env
)
assign(
  "fit_pseudobulk_deseq2",
  old_fit_pseudobulk,
  envir = prepare_pseudobulk_env
)
assign(
  "select_pseudobulk_deseq2",
  old_select_pseudobulk,
  envir = prepare_pseudobulk_env
)
assign("exec_time", old_exec_time, envir = prepare_pseudobulk_env)
assign("peak_rss_gb", old_peak_rss, envir = prepare_pseudobulk_env)
unlink(fixture_h5ad)

single_batch_error <- tryCatch(
  correct_clr_batch_lmm(feat, transform(meta, tech = factor("A")), "tech"),
  error = identity
)
stopifnot(inherits(single_batch_error, "error"))

# The correction model is batch-only by construction.
correction_body <- paste(deparse(body(correct_clr_batch_lmm)), collapse = " ")
stopifnot(grepl("y ~ 1 \\+ \\(1 \\| batch\\)", correction_body))
stopifnot(!grepl("label|Status|disease", correction_body, ignore.case = TRUE))

# Batch pseudobulk result bundles use a pass-qualified stem and only the
# hvg2000 high-resolution result.
process_pseudobulk_fig <- function(feat_mat, labels, ...) {
  list(feat_mat = feat_mat, labels = labels)
}
tmp_results <- tempfile()
dir.create(tmp_results)
pb <- list(
  hvg2000 = list(
    pb = matrix(1, nrow = 3, ncol = 3,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:3))),
    time_secs = 0
  )
)
batch_pb <- run_pseudobulk_hpc(
  list(),
  labels = structure(factor(c("A", "B", "A")), names = paste0("s", 1:3)),
  pb_variants = pb,
  results_dir = tmp_results,
  ds = "DS",
  batch_mode = TRUE,
  result_stem = "DS_batch_effect_uncorrected"
)

stopifnot(identical(names(batch_pb), "Pseudobulk_hvg2000"))
stopifnot(file.exists(file.path(
  tmp_results, "DS_batch_effect_uncorrected_Pseudobulk_hvg2000.rds"
)))

cat("batch-effect CLR and pseudobulk modes OK\n")
# CT pseudobulk must fail closed when every per-CT normalization fails, and
# successful runs must expose contribution counts in the result bundle.
ct_counts <- matrix(
  1,
  nrow = 4,
  ncol = 6,
  dimnames = list(paste0("g", 1:4), paste0("cell", 1:6))
)
ct_meta <- data.frame(
  Sample = rep(paste0("s", 1:3), each = 2),
  ct = rep(c("A", "B"), 3),
  row.names = colnames(ct_counts)
)
ct_seurat <- Seurat::CreateSeuratObject(
  counts = ct_counts,
  meta.data = ct_meta
)
ct_labels <- structure(
  factor(c("A", "B", "A")),
  names = paste0("s", 1:3)
)
get_pb_deseq2 <- function(...) {
  stop("synthetic per-CT failure")
}
all_failed <- tryCatch(
  process_pseudobulk_ct_fig(
    ct_seurat,
    ct_labels,
    ct_col = "ct",
    sample_col = "Sample",
    hvg = 2,
    min_cells = 1
  ),
  error = identity
)
stopifnot(
  inherits(all_failed, "error"),
  grepl("no successful cell-type pseudobulks", conditionMessage(all_failed))
)

create_result_bundle <- function(feat_mat, labels, dist_mat = NULL, extra = list()) {
  c(list(feat_mat = feat_mat, labels = labels), extra)
}
get_pb_deseq2 <- function(...) {
  matrix(
    seq_len(6),
    nrow = 3,
    ncol = 2,
    dimnames = list(paste0("s", 1:3), c("g1", "g2"))
  )
}
ct_success <- process_pseudobulk_ct_fig(
  ct_seurat,
  ct_labels,
  ct_col = "ct",
  sample_col = "Sample",
  hvg = 2,
  min_cells = 1
)
stopifnot(
  identical(ct_success$n_ct_success, 2L),
  identical(ct_success$n_sample_pairs_contributed, 3L),
  identical(ct_success$n_ct_pair_contributions, 6L),
  identical(ct_success$successful_cell_types, c("A", "B"))
)

# Canonical CT pseudobulk uses a fixture-only H5AD and the one-pass composite
# store.  The fixture deliberately has unsorted samples, an absent combination,
# a below-threshold group, and one normalization failure.
ct_h5ad_root <- tempfile("ecoda-ct-h5ad-")
dir.create(ct_h5ad_root, recursive = TRUE)
ct_h5ad <- file.path(ct_h5ad_root, "fixture.h5ad")
ct_python <- paste(
  "import anndata as ad, numpy as np, pandas as pd, sys;",
  "from scipy import sparse;",
  "samples = (['s2'] * 5 + ['s3'] * 5 + ['s1'] * 5 +",
  " ['s2'] * 5 + ['s1'] * 5 + ['s3'] * 4 +",
  " ['s2'] * 5 + ['s1'] * 5);",
  "cell_types = (['B'] * 15 + ['A'] * 14 + ['C'] * 10);",
  "b_values = {'s1': 1, 's2': 3, 's3': 7};",
  "a_values = {'s1': 30, 's2': 4, 's3': 8};",
  "values = [99 if ct == 'C' else (",
  " a_values if ct == 'A' else b_values)[sample]",
  " for sample, ct in zip(samples, cell_types)];",
  "counts = np.asarray([[value, 1]",
  " for i, value in enumerate(values)], dtype=np.int64);",
  "obs = pd.DataFrame({'Sample': samples, 'ct': cell_types},",
  " index=[f'cell{i}' for i in range(len(samples))]);",
  "adata = ad.AnnData(X=counts.astype(np.float32), obs=obs,",
  " var=pd.DataFrame(index=['g1', 'g2']));",
  "adata.layers['counts'] = sparse.csr_matrix(counts);",
  "adata.write_h5ad(sys.argv[1])"
)
ct_status <- system2(
  "pixi",
  c("run", "python", "-c", shQuote(ct_python), shQuote(ct_h5ad)),
  stdout = FALSE,
  stderr = FALSE
)
if (!identical(ct_status, 0L) || !file.exists(ct_h5ad)) {
  stop("could not create the synthetic CT H5AD fixture")
}
old_project_root <- Sys.getenv("PROJECT_ROOT", unset = NA_character_)
Sys.setenv(PROJECT_ROOT = root)
old_direct_ct_normalizer <- get_pb_deseq2_from_counts
get_pb_deseq2_from_counts <- function(counts, metadata, ...) {
  # The direct CT boundary receives per-group aggregates: A/s1 is 5 * 30 =
  # 150, while C aggregates are 5 * 99 = 495.  Fail only the C-sized
  # aggregates without confusing a valid high-count A group for C.
  if (max(counts) >= 400) stop("synthetic C normalization failure")
  sample_ids <- colnames(counts)
  matrix(
    as.numeric(colSums(counts)),
    nrow = length(sample_ids),
    ncol = 1L,
    dimnames = list(sample_ids, "score")
  )
}
old_ct_bundle <- create_result_bundle
create_result_bundle <- function(feat_mat, labels, dist_mat = NULL, extra = list()) {
  c(list(feat_mat = feat_mat, dist_mat = dist_mat, labels = labels), extra)
}
ct_temp_root <- tempfile("ecoda-ct-store-root-")
ct_labels_h5ad <- structure(
  factor(c("A", "B", "A")),
  names = c("s1", "s2", "s3")
)
ct_result_h5ad <- tryCatch(
  process_pseudobulk_ct_h5ad_fig(
    ct_h5ad,
    ct_labels_h5ad,
    sample_col = "Sample",
    ct_col = "ct",
    hvg = 1L,
    min_cells = 5L,
    chunk_size = 4L,
    run_id = "ct-fixture-run",
    temp_root = ct_temp_root,
    source_identity = "fixture-source"
  ),
  error = function(error) {
    create_result_bundle <<- old_ct_bundle
    get_pb_deseq2_from_counts <<- old_direct_ct_normalizer
    stop(error)
  }
)
create_result_bundle <- old_ct_bundle
get_pb_deseq2_from_counts <- old_direct_ct_normalizer
stopifnot(
  identical(ct_result_h5ad$successful_cell_types, c("B", "A")),
  identical(ct_result_h5ad$n_ct_success, 2L),
  identical(ct_result_h5ad$n_sample_pairs_contributed, 3L),
  identical(ct_result_h5ad$n_ct_pair_contributions, 4L),
  identical(
    rownames(ct_result_h5ad$feat_mat),
    c("s1", "s2", "s3")
  ),
  identical(
    colnames(ct_result_h5ad$feat_mat),
    c("s1", "s2", "s3")
  )
)
ct_expected_distance <- matrix(
  c(
    0, 70, 30,
    70, 0, 20,
    30, 20, 0
  ),
  nrow = 3L,
  byrow = TRUE,
  dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
)
stopifnot(isTRUE(all.equal(
  ct_result_h5ad$feat_mat,
  ct_expected_distance,
  tolerance = 1e-12
)))
ct_distance_matrix <- as.matrix(ct_result_h5ad$dist_mat)
stopifnot(isTRUE(all.equal(
  ct_distance_matrix,
  ct_expected_distance,
  tolerance = 1e-12
)))
ct_store_paths <- if (dir.exists(ct_temp_root)) {
  list.files(ct_temp_root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
} else {
  character()
}
stopifnot(!any(grepl(
  "(groups\\.h5|\\.manifest\\.json|\\.lock)$",
  ct_store_paths,
  perl = TRUE
)))

all_failed_temp_root <- tempfile("ecoda-ct-store-failure-root-")
old_direct_ct_normalizer <- get_pb_deseq2_from_counts
get_pb_deseq2_from_counts <- function(...) {
  stop("synthetic normalization failure")
}
all_failed_h5ad <- tryCatch(
  process_pseudobulk_ct_h5ad_fig(
    ct_h5ad,
    ct_labels_h5ad,
    sample_col = "Sample",
    ct_col = "ct",
    hvg = 1L,
    min_cells = 5L,
    chunk_size = 4L,
    run_id = "ct-failure-run",
    temp_root = all_failed_temp_root,
    source_identity = "fixture-source"
  ),
  error = identity
)
get_pb_deseq2_from_counts <- old_direct_ct_normalizer
stopifnot(
  inherits(all_failed_h5ad, "error"),
  grepl("no successful cell-type pseudobulks", conditionMessage(all_failed_h5ad))
)
failed_store_paths <- if (dir.exists(all_failed_temp_root)) {
  list.files(
    all_failed_temp_root,
    recursive = TRUE,
    all.files = TRUE,
    full.names = TRUE
  )
} else {
  character()
}
stopifnot(!any(grepl(
  "(groups\\.h5|\\.manifest\\.json|\\.lock)$",
  failed_store_paths,
  perl = TRUE
)))

# Exercise the store API directly so the ownership manifest and conservative
# stale-store cleanup remain covered independently of the canonical wrapper.
ct_module <- reticulate::import_from_path(
  "h5ad_pseudobulk",
  path = file.path(root, "src", "utils", "py"),
  convert = TRUE
)
ct_store_root <- tempfile("ecoda-ct-store-api-")
dir.create(ct_store_root, recursive = TRUE)
ct_store_path <- file.path(ct_store_root, "groups.h5")
ct_store_payload <- ct_module$prepare_h5ad_ct_group_store(
  path = ct_h5ad,
  sample_col = "Sample",
  cell_type_col = "ct",
  metadata_columns = as.list(c("Sample", "ct")),
  chunk_size = 4L,
  max_value = as.integer(.Machine$integer.max),
  store_path = ct_store_path,
  run_id = "store-fixture-run",
  source_identity = "fixture-source"
)
ct_manifest_path <- paste0(ct_store_path, ".manifest.json")
ct_manifest <- jsonlite::fromJSON(ct_manifest_path, simplifyVector = FALSE)
write_ct_manifest <- function(value) {
  jsonlite::write_json(
    value,
    ct_manifest_path,
    auto_unbox = TRUE,
    null = "null",
    pretty = FALSE
  )
  manifest_json <- paste(readLines(ct_manifest_path, warn = FALSE), collapse = "")
  manifest_json <- sub(
    '"scheduler_identity":\\[\\]',
    '"scheduler_identity":{}',
    manifest_json,
    fixed = FALSE
  )
  writeLines(manifest_json, ct_manifest_path, useBytes = TRUE)
}
stopifnot(
  is.list(ct_store_payload),
  identical(
    sort(names(ct_store_payload)),
    sort(c(
      "store_path", "group_ids", "sample_ids", "all_sample_ids",
      "cell_type_ids", "group_cell_counts", "group_metadata",
      "gene_names", "n_vars"
    ))
  ),
  identical(
    as.character(ct_store_payload$store_path),
    normalizePath(ct_store_path)
  ),
  identical(
    as.character(ct_store_payload$group_ids),
    c(
      "Sample=s2;cell_type=B",
      "Sample=s3;cell_type=B",
      "Sample=s1;cell_type=B",
      "Sample=s2;cell_type=A",
      "Sample=s1;cell_type=A",
      "Sample=s3;cell_type=A",
      "Sample=s2;cell_type=C",
      "Sample=s1;cell_type=C"
    )
  ),
  identical(
    as.character(ct_store_payload$sample_ids),
    c("s2", "s3", "s1", "s2", "s1", "s3", "s2", "s1")
  ),
  identical(
    as.character(ct_store_payload$all_sample_ids),
    c("s2", "s3", "s1")
  ),
  identical(
    as.character(ct_store_payload$cell_type_ids),
    c("B", "B", "B", "A", "A", "A", "C", "C")
  ),
  identical(
    as.integer(ct_store_payload$group_cell_counts),
    c(5L, 5L, 5L, 5L, 5L, 4L, 5L, 5L)
  ),
  identical(
    colnames(ct_store_payload$group_metadata),
    c("Sample", "ct")
  ),
  identical(
    as.character(ct_store_payload$group_metadata$Sample),
    c("s2", "s3", "s1", "s2", "s1", "s3", "s2", "s1")
  ),
  identical(
    as.character(ct_store_payload$group_metadata$ct),
    c("B", "B", "B", "A", "A", "A", "C", "C")
  ),
  identical(
    as.character(ct_store_payload$gene_names),
    c("g1", "g2")
  ),
  identical(as.integer(ct_store_payload$n_vars), 2L),
  all(c(
    "run_id", "pid", "scheduler_identity", "source_identity",
    "source_checksum", "stage", "schema", "created_at"
  ) %in% names(ct_manifest)),
  identical(as.character(ct_manifest$run_id), "store-fixture-run"),
  identical(as.character(ct_manifest$stage), "pseudobulk_ct"),
  identical(as.integer(ct_manifest$schema), 1L)
)
ct_active_audit <- ct_module$audit_h5ad_ct_group_store(
  ct_store_path,
  expected_run_id = "store-fixture-run",
  max_age_seconds = 0,
  cleanup = TRUE
)
stopifnot(
  isTRUE(ct_active_audit$valid),
  !isTRUE(ct_active_audit$cleanup_performed),
  file.exists(ct_store_path)
)
unlink(ct_manifest_path)
ct_malformed_audit <- ct_module$audit_h5ad_ct_group_store(
  ct_store_path,
  expected_run_id = "store-fixture-run",
  max_age_seconds = 0,
  cleanup = TRUE
)
stopifnot(
  !isTRUE(ct_malformed_audit$valid),
  !isTRUE(ct_malformed_audit$cleanup_performed),
  file.exists(ct_store_path)
)
write_ct_manifest(ct_manifest)
ct_manifest$pid <- 2147483647
ct_manifest$scheduler_identity <- list()
ct_manifest$created_at <- 0
write_ct_manifest(ct_manifest)
ct_stale_audit <- ct_module$audit_h5ad_ct_group_store(
  ct_store_path,
  expected_run_id = "store-fixture-run",
  max_age_seconds = 1,
  cleanup = TRUE
)
stopifnot(
  isTRUE(ct_stale_audit$valid),
  identical(as.character(ct_stale_audit$owner_status), "dead"),
  isTRUE(ct_stale_audit$expired),
  isTRUE(ct_stale_audit$cleanup_performed),
  !file.exists(ct_store_path),
  !file.exists(ct_manifest_path)
)
unlink(
  c(ct_h5ad_root, ct_temp_root, all_failed_temp_root, ct_store_root),
  recursive = TRUE,
  force = TRUE
)
if (is.na(old_project_root)) {
  Sys.unsetenv("PROJECT_ROOT")
} else {
  Sys.setenv(PROJECT_ROOT = old_project_root)
}

cat("batch-effect CLR, pseudobulk modes, and CT contribution guard OK\n")
