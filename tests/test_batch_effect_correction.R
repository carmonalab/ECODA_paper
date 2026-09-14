#!/usr/bin/env Rscript

raw_args <- commandArgs(trailingOnly = FALSE)
script_arg <- raw_args[grepl("^--file=", raw_args)][1]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), ".."))
if (!nzchar(Sys.getenv("PROJECT_ROOT", unset = ""))) {
  Sys.setenv(PROJECT_ROOT = root)
}
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
source(file.path(root, "src/utils/seurat_utils.R"))
source(file.path(root, "src/utils/pseudobulk.R"))
source(file.path(root, "src/utils/scoring_metrics.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_methods_r.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_pipeline.R"))

# Corrected-final H5ADs intentionally carry configuration-only identity. The R
# path validator must pass the explicit summary-free opt-in through to Python,
# even when the consumer method is not a Pipeline 3 method name.
summary_free_h5ad <- tempfile(
  "ecoda-summary-free-h5ad-", fileext = ".h5ad"
)
summary_free_h5ad_python <- paste(
  "import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, sys;",
  "sys.path.insert(0, sys.argv[2]);",
  "from src.utils.py.batch_contract import build_batch_contract_identity;",
  "n_genes = 2000;",
  "obs = pd.DataFrame({'Sample': ['s1', 's2'],",
  "                    'assay': ['a', 'a'], 'sex': ['F', 'M']},",
  "                   index=['c1', 'c2']);",
  "adata = ad.AnnData(",
  "X=np.ones((2, n_genes), dtype=np.float32), obs=obs,",
  "var=pd.DataFrame({'hvg_rank': np.arange(1, n_genes + 1, dtype=float)},",
  "                  index=[f'g{i}' for i in range(n_genes)]));",
  "adata.layers['counts'] = sp.csr_matrix(",
  "np.ones((2, n_genes), dtype=np.int64));",
  "adata.obsm['X_pca_batch_effect_corrected_hvg2000'] =",
  "np.ones((2, 2), dtype=np.float32);",
  "adata.obsm['X_pca_harmony_batch_effect_corrected_hvg2000'] =",
  "np.ones((2, 2), dtype=np.float32);",
  "adata.uns['batch_contract'] = build_batch_contract_identity(",
  "['assay', 'sex'], sample_column='Sample', method_id='preprocess',",
  "model_id='hvg_composite_v1');",
  "adata.write_h5ad(sys.argv[1])"
)
summary_free_h5ad_status <- system2(
  "pixi",
  c(
    "run", "python", "-c", shQuote(summary_free_h5ad_python),
    shQuote(summary_free_h5ad), shQuote(root)
  ),
  stdout = FALSE,
  stderr = FALSE
)
if (!identical(summary_free_h5ad_status, 0L)) {
  stop("could not create the summary-free corrected H5AD fixture")
}
old_analysis_variant <- Sys.getenv("ANALYSIS_VARIANT", unset = "")
old_analysis_pass <- Sys.getenv("ANALYSIS_PASS", unset = "")
Sys.setenv(
  ANALYSIS_VARIANT = "corrected_final",
  ANALYSIS_PASS = "corrected"
)
summary_free_identity <- ecoda_hpc_batch_contract_identity(
  batch_keys = list("assay", "sex"),
  sample_col = "Sample",
  method_id = "preprocess",
  model_id = "hvg_composite_v1"
)
ecoda_hpc_validate_h5ad_path_identity(
  h5ad_path = summary_free_h5ad,
  view = "batch_effect_corrected",
  method = "gloscope",
  expected_batch_contract = summary_free_identity,
  allow_missing_summary = TRUE
)
if (nzchar(old_analysis_variant)) {
  Sys.setenv(ANALYSIS_VARIANT = old_analysis_variant)
} else {
  Sys.unsetenv("ANALYSIS_VARIANT")
}
if (nzchar(old_analysis_pass)) {
  Sys.setenv(ANALYSIS_PASS = old_analysis_pass)
} else {
  Sys.unsetenv("ANALYSIS_PASS")
}
unlink(summary_free_h5ad)

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

# Corrected CLR uses limma's documented fixed-effect boundary: build a
# separate technical-key model matrix, preserve only the intercept in
# removeBatchEffect(), and restore the CLR row-sum invariant.
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
one_key_validation <- ecoda_batch_validate_metadata(
  meta,
  batch_keys = "tech",
  sample_col = "Sample"
)
one_key_design <- ecoda_batch_fixed_effect_design(
  metadata = meta,
  batch_keys = "tech",
  validation = one_key_validation,
  sample_col = "Sample"
)
corrected <- correct_clr_batch_limma(
  feat,
  meta,
  batch_keys = "tech",
  sample_col = "Sample",
  metadata_validation = one_key_validation
)
limma_reference <- limma::removeBatchEffect(
  x = t(feat),
  covariates = one_key_design$design[
    , one_key_design$technical_columns,
    drop = FALSE
  ],
  design = matrix(
    1,
    nrow = nrow(one_key_design$design),
    ncol = 1L,
    dimnames = list(rownames(one_key_design$design), "(Intercept)")
  )
)
limma_reference <- t(limma_reference)
limma_reference <- sweep(
  limma_reference,
  1L,
  rowMeans(limma_reference),
  FUN = "-"
)
limma_reference[, ncol(limma_reference)] <- -rowSums(
  limma_reference[, -ncol(limma_reference), drop = FALSE]
)
dimnames(limma_reference) <- dimnames(feat)
stopifnot(
  identical(one_key_validation$effective_batch_keys, "tech"),
  identical(one_key_validation$non_estimable_batch_keys, character()),
  identical(one_key_design$rank, 2L),
  identical(one_key_design$columns, 2L),
  one_key_design$residual_df > 0L,
  isTRUE(all.equal(corrected, limma_reference, tolerance = 1e-8)),
  all(is.finite(corrected)),
  max(abs(rowSums(corrected))) < 1e-8,
  identical(dimnames(corrected), dimnames(feat)),
  !"__ecoda_batch_combined_v1" %in% colnames(meta)
)

bad_order <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]
order_error <- tryCatch(
  correct_clr_batch_limma(
    feat,
    bad_order,
    batch_keys = "tech",
    sample_col = "Sample"
  ),
  error = identity
)
stopifnot(inherits(order_error, "error"))

missing_batch <- meta
missing_batch$tech[1] <- NA
missing_error <- tryCatch(
  correct_clr_batch_limma(
    feat,
    missing_batch,
    batch_keys = "tech",
    sample_col = "Sample"
  ),
  error = identity
)
stopifnot(inherits(missing_error, "error"))

# The canonical pseudobulk driver uses one raw H5AD aggregate and forwards the
# limma-corrected path through an intercept-only DESeq2 fit without exposing
# biological labels or synthesizing a combined batch column.
captured <- new.env(parent = emptyenv())
fake_counts <- matrix(
  c(1, 2, 3, 4, 5, 6, 7, 8),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(c("g1", "g2"), paste0("s", 1:4))
)
fake_metadata <- data.frame(
  Sample = paste0("s", 1:4),
  tech = factor(c("A", "B", "A", "B")),
  row.names = paste0("s", 1:4),
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
fixture_cell_metadata <- fake_metadata[
  rep(seq_len(nrow(fake_metadata)), each = 2L),
  ,
  drop = FALSE
]
rownames(fixture_cell_metadata) <- paste0(
  "fixture-cell",
  seq_len(nrow(fixture_cell_metadata))
)
fixture_batch_context <- ecoda_hpc_batch_context(
  metadata = fixture_cell_metadata,
  batch_keys = "tech",
  sample_col = "Sample"
)
fixture_h5ad <- tempfile("ecoda-pseudobulk-driver-", fileext = ".h5ad")
fixture_h5ad_python <- paste(
  "import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, sys;",
  "adata = ad.AnnData(",
  "X=np.ones((4, 1), dtype=np.float32),",
  "obs=pd.DataFrame({'Sample': ['s1', 's2', 's3', 's4'],",
  "                  'tech': ['A', 'B', 'A', 'B']},",
  "                 index=['c1', 'c2', 'c3', 'c4']),",
  "var=pd.DataFrame(index=['g1']));",
  "adata.layers['counts'] = sp.csr_matrix(np.ones((4, 1), dtype=np.int64));",
  "adata.write_h5ad(sys.argv[1])"
)
fixture_h5ad_status <- system2(
  "pixi",
  c(
    "run", "python", "-c", shQuote(fixture_h5ad_python),
    shQuote(fixture_h5ad)
  ),
  stdout = FALSE,
  stderr = FALSE
)
if (!identical(fixture_h5ad_status, 0L) || !file.exists(fixture_h5ad)) {
  stop("could not create the synthetic pseudobulk H5AD fixture")
}
corrected_pb <- prepare_pseudobulks_hpc(
  h5ad_path = fixture_h5ad,
  hvg_rank_genes = c("g1", "g2"),
  variants = "hvg2000",
  batch_col = "tech",
  blind = FALSE,
  correct_batch = TRUE,
  cache_stem = "fixture",
  view = "batch_effect_corrected",
  analysis_pass = "corrected",
  run_id = "fixture-run",
  batch_context = fixture_batch_context
)
stopifnot(
  identical(captured$fit$counts, fake_counts),
  identical(captured$fit$metadata, fake_metadata),
  !"label" %in% colnames(captured$fit$metadata),
  identical(captured$fit$batch_col, NULL),
  identical(captured$fit$blind, FALSE),
  identical(captured$fit$correct_batch, FALSE),
  identical(captured$aggregate$metadata_columns, c("Sample", "tech")),
  !"label" %in% captured$aggregate$metadata_columns,
  !"__ecoda_batch_combined_v1" %in% colnames(captured$fit$metadata),
  is.matrix(corrected_pb$hvg2000$pb),
  all(is.finite(corrected_pb$hvg2000$pb)),
  identical(
    dimnames(corrected_pb$hvg2000$pb),
    list(paste0("s", 1:4), c("g1", "g2"))
  ),
  identical(
    corrected_pb$hvg2000$batch_contract$model_id,
    "pseudobulk_limma_fixed_effects_v1"
  ),
  identical(
    corrected_pb$hvg2000$batch_contract$correction_mode,
    "limma_fixed_effects_pseudobulk"
  ),
  grepl(
    "DESeq2 design=~ 1; model.matrix(~ 1 + batch_key_1); limma::removeBatchEffect",
    corrected_pb$hvg2000$batch_contract$correction_formula,
    fixed = TRUE
  ),
  !grepl(
    "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
    corrected_pb$hvg2000$batch_contract$correction_formula,
    perl = TRUE
  )
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

no_correction_cell_metadata <- fixture_cell_metadata
no_correction_cell_metadata$tech <- factor(
  rep("A", nrow(no_correction_cell_metadata))
)
no_correction_validation <- ecoda_batch_validate_metadata(
  no_correction_cell_metadata,
  batch_keys = "tech",
  sample_col = "Sample"
)
no_correction_context <- ecoda_hpc_batch_context(
  metadata = no_correction_cell_metadata,
  batch_keys = "tech",
  sample_col = "Sample"
)
no_correction_pb <- prepare_pseudobulks_hpc(
  h5ad_path = fixture_h5ad,
  hvg_rank_genes = c("g1", "g2"),
  variants = "hvg2000",
  batch_col = "tech",
  blind = FALSE,
  correct_batch = TRUE,
  cache_stem = "fixture-no-correction",
  view = "batch_effect_corrected",
  analysis_pass = "corrected",
  run_id = "fixture-no-correction-run",
  batch_context = no_correction_context,
  batch_contract = ecoda_hpc_batch_contract_identity(
    batch_keys = "tech",
    sample_col = "Sample",
    method_id = "Pseudobulk",
    model_id = "pseudobulk_limma_fixed_effects_v1"
  )
)
stopifnot(
  identical(no_correction_validation$correction_state, "NO_CORRECTION"),
  identical(captured$fit$batch_col, NULL),
  identical(captured$fit$blind, TRUE),
  identical(captured$fit$correct_batch, FALSE),
  identical(
    no_correction_pb$hvg2000$pb,
    t(fake_counts)
  ),
  identical(
    no_correction_pb$hvg2000$batch_contract$correction_state,
    "NO_CORRECTION"
  ),
  identical(
    no_correction_pb$hvg2000$batch_contract$correction_formula,
    "NO_CORRECTION: no estimable technical batch key"
  )
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

single_batch_meta <- transform(meta, tech = factor("A"))
single_batch_validation <- ecoda_batch_validate_metadata(
  single_batch_meta,
  batch_keys = "tech",
  sample_col = "Sample"
)
single_batch_result <- correct_clr_batch_limma(
  feat,
  single_batch_meta,
  batch_keys = "tech",
  sample_col = "Sample",
  metadata_validation = single_batch_validation
)
stopifnot(
  identical(single_batch_validation$correction_state, "NO_CORRECTION"),
  identical(single_batch_validation$effective_batch_keys, character()),
  identical(single_batch_result, feat)
)

# A constant component of a multi-key technical design remains in the
# contract and in the metadata, while the estimable varying component drives
# the separate limma covariate correction.
multikey_meta <- data.frame(
  Sample = rownames(feat),
  assay = factor(rep("10x 3' v3", n)),
  sex = factor(batch),
  stringsAsFactors = FALSE
)
multikey_validation <- ecoda_hpc_sample_metadata_validation(
  multikey_meta,
  batch_keys = list("assay", "sex"),
  sample_col = "Sample"
)
stopifnot(
  identical(multikey_validation$key_level_counts, list(assay = 1L, sex = 2L)),
  identical(multikey_validation$effective_batch_keys, "sex"),
  identical(multikey_validation$non_estimable_batch_keys, "assay"),
  identical(multikey_validation$fixed_effect_aliases, c(sex = "batch_key_2")),
  identical(multikey_validation$correction_state, "BATCH_CORRECTION"),
  identical(multikey_validation$design_rank, 2L),
  identical(multikey_validation$design_columns, 2L),
  multikey_validation$design_residual_df > 0L,
  length(multikey_validation$composite_levels) == 2L,
  all(c("assay", "sex") %in% colnames(multikey_validation$sample_metadata)),
  !"__ecoda_batch_combined_v1" %in% colnames(multikey_validation$sample_metadata)
)
multikey_identity <- ecoda_hpc_batch_contract_identity(
  batch_keys = list("assay", "sex"),
  sample_col = "Sample",
  method_id = "ECODA_authors_HR",
  model_id = "limma_fixed_effects_v1"
)
multikey_contract <- ecoda_hpc_augment_batch_contract(
  identity = multikey_identity,
  validation = multikey_validation,
  method_id = "ECODA_authors_HR",
  batch_keys = list("assay", "sex")
)
multikey_contract_again <- ecoda_hpc_augment_batch_contract(
  identity = multikey_contract,
  validation = multikey_validation,
  method_id = "ECODA_authors_HR",
  batch_keys = list("assay", "sex")
)
stopifnot(
  identical(multikey_contract_again, multikey_contract),
  identical(multikey_contract$effective_batch_keys, "sex"),
  identical(multikey_contract$non_estimable_batch_keys, "assay"),
  identical(multikey_contract$correction_state, "BATCH_CORRECTION"),
  identical(multikey_contract$correction_mode, "limma_fixed_effects"),
  identical(multikey_contract$fixed_effect_aliases, c(sex = "batch_key_2")),
  identical(multikey_contract$design_rank, 2L),
  identical(multikey_contract$design_columns, 2L),
  multikey_contract$design_residual_df > 0L,
  grepl(
    "model.matrix(~ 1 + batch_key_2); limma::removeBatchEffect",
    multikey_contract$correction_formula,
    fixed = TRUE
  ),
  !grepl(
    "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
    multikey_contract$correction_formula,
    perl = TRUE
  )
)
ecoda_hpc_validate_batch_contract(
  multikey_contract,
  multikey_contract,
  label = "strict limma corrected batch contract",
  require_effective_metadata = TRUE
)
strict_source_error <- tryCatch(
  ecoda_hpc_validate_batch_contract(
    multikey_identity,
    multikey_contract,
    label = "strict corrected batch contract",
    require_effective_metadata = TRUE
  ),
  error = identity
)
stopifnot(inherits(strict_source_error, "error"))
multikey_corrected <- correct_clr_batch_limma(
  feat,
  multikey_meta,
  batch_keys = list("assay", "sex"),
  sample_col = "Sample",
  metadata_validation = multikey_validation
)
stopifnot(
  all(is.finite(multikey_corrected)),
  identical(dimnames(multikey_corrected), dimnames(feat)),
  max(abs(rowSums(multikey_corrected))) < 1e-8
)

# A two-level assay/sex design keeps both original categorical columns as
# separate estimable limma covariates.
assay_sex_meta <- data.frame(
  Sample = rownames(feat),
  assay = factor(rep(c("10x 3' v3", "10x multiome"), each = 6)),
  sex = factor(rep(c("F", "M"), times = 6)),
  stringsAsFactors = FALSE
)
assay_sex_validation <- ecoda_batch_validate_metadata(
  assay_sex_meta,
  batch_keys = list("assay", "sex"),
  sample_col = "Sample"
)
assay_sex_design <- ecoda_batch_fixed_effect_design(
  metadata = assay_sex_meta,
  batch_keys = list("assay", "sex"),
  validation = assay_sex_validation,
  sample_col = "Sample"
)
assay_sex_corrected <- correct_clr_batch_limma(
  feat,
  assay_sex_meta,
  batch_keys = list("assay", "sex"),
  sample_col = "Sample",
  metadata_validation = assay_sex_validation
)
assay_sex_contract <- ecoda_hpc_augment_batch_contract(
  identity = ecoda_hpc_batch_contract_identity(
    batch_keys = list("assay", "sex"),
    sample_col = "Sample",
    method_id = "ECODA_authors_HR",
    model_id = "limma_fixed_effects_v1"
  ),
  validation = assay_sex_validation,
  method_id = "ECODA_authors_HR",
  batch_keys = list("assay", "sex")
)
stopifnot(
  identical(assay_sex_validation$effective_batch_keys, c("assay", "sex")),
  identical(assay_sex_validation$non_estimable_batch_keys, character()),
  identical(assay_sex_validation$fixed_effect_aliases, c(
    assay = "batch_key_1",
    sex = "batch_key_2"
  )),
  identical(assay_sex_design$rank, 3L),
  identical(assay_sex_design$columns, 3L),
  assay_sex_design$residual_df > 0L,
  identical(assay_sex_contract$correction_mode, "limma_fixed_effects"),
  identical(assay_sex_contract$design_rank, 3L),
  identical(assay_sex_contract$design_columns, 3L),
  assay_sex_contract$design_residual_df > 0L,
  all(c("assay", "sex") %in% names(assay_sex_contract$validation_summary$per_key_levels)),
  all(is.finite(assay_sex_corrected)),
  identical(dimnames(assay_sex_corrected), dimnames(feat)),
  max(abs(rowSums(assay_sex_corrected))) < 1e-8,
  grepl(
    "model.matrix(~ 1 + batch_key_1 + batch_key_2); limma::removeBatchEffect(covariates=technical_covariates, design=intercept)",
    assay_sex_contract$correction_formula,
    fixed = TRUE
  ),
  !grepl(
    "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
    assay_sex_contract$correction_formula,
    perl = TRUE
  )
)

assay_sex_pb_matrix <- matrix(
  rnorm(3L * n),
  nrow = 3L,
  dimnames = list(paste0("pb-g", 1:3), rownames(feat))
)
assay_sex_pb <- ecoda_hpc_apply_limma_batch_correction(
  fit = list(norm_matrix = assay_sex_pb_matrix),
  metadata = assay_sex_meta,
  batch_keys = list("assay", "sex"),
  sample_col = "Sample"
)
stopifnot(
  isTRUE(assay_sex_pb$correct_batch),
  is.null(assay_sex_pb$batch_col),
  all(is.finite(assay_sex_pb$norm_matrix)),
  identical(dimnames(assay_sex_pb$norm_matrix), dimnames(assay_sex_pb_matrix)),
  identical(assay_sex_pb$batch_correction$correction_mode, "limma_fixed_effects_pseudobulk"),
  identical(
    assay_sex_pb$batch_correction$configured_batch_keys,
    c("assay", "sex")
  ),
  identical(
    assay_sex_pb$batch_correction$effective_batch_keys,
    c("assay", "sex")
  ),
  assay_sex_pb$batch_correction$design_rank == 3L,
  assay_sex_pb$batch_correction$design_columns == 3L,
  assay_sex_pb$batch_correction$design_residual_df > 0L,
  grepl(
    "DESeq2 design=~ 1; model.matrix(~ 1 + batch_key_1 + batch_key_2); limma::removeBatchEffect(covariates=technical_covariates, design=intercept)",
    assay_sex_pb$batch_correction$correction_formula,
    fixed = TRUE
  ),
  !grepl(
    "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
    assay_sex_pb$batch_correction$correction_formula,
    perl = TRUE
  )
)

# A Breast-like three-key design uses the real configured technical names.
# Every key remains an original metadata column; only fixed internal aliases
# enter the additive model.
breast_samples <- paste0("breast-s", seq_len(12))
breast_sample_meta <- data.frame(
  Sample = breast_samples,
  disease = factor(rep(c("tumor", "normal"), each = 6)),
  assay = factor(rep(c("10x 3' v3", "10x multiome"), each = 6)),
  sequencing_platform = factor(
    rep(rep(c("NovaSeq", "Illumina"), each = 3), 2)
  ),
  suspension_dissociation_time = factor(
    rep(c("fresh", "frozen", "ambient"), 4)
  ),
  stringsAsFactors = FALSE
)
breast_keys <- c(
  "assay",
  "sequencing_platform",
  "suspension_dissociation_time"
)
breast_cell_meta <- breast_sample_meta[
  rep(seq_len(nrow(breast_sample_meta)), each = 2L),
  ,
  drop = FALSE
]
rownames(breast_cell_meta) <- paste0("breast-cell-", seq_len(nrow(breast_cell_meta)))
breast_validation <- ecoda_batch_validate_metadata(
  breast_cell_meta,
  batch_keys = as.list(breast_keys),
  sample_col = "Sample",
  biological_label = "disease"
)
breast_context <- ecoda_hpc_batch_context(
  metadata = breast_cell_meta,
  batch_keys = as.list(breast_keys),
  sample_col = "Sample",
  biological_label = "disease"
)
breast_design <- ecoda_batch_fixed_effect_design(
  metadata = breast_context$sample_metadata,
  batch_keys = as.list(breast_keys),
  validation = breast_context$validation,
  sample_col = "Sample"
)
breast_feat <- matrix(
  rnorm(length(breast_samples) * 4L),
  nrow = length(breast_samples),
  dimnames = list(breast_samples, paste0("breast-ct", 1:4))
)
breast_feat <- breast_feat - rowMeans(breast_feat)
breast_corrected <- correct_clr_batch_limma(
  breast_feat,
  breast_sample_meta,
  batch_keys = as.list(breast_keys),
  sample_col = "Sample",
  metadata_validation = breast_validation
)
breast_contract <- ecoda_hpc_augment_batch_contract(
  identity = ecoda_hpc_batch_contract_identity(
    batch_keys = as.list(breast_keys),
    sample_col = "Sample",
    method_id = "ECODA_authors_HR",
    model_id = "limma_fixed_effects_v1"
  ),
  validation = breast_validation,
  method_id = "ECODA_authors_HR",
  batch_keys = as.list(breast_keys)
)
breast_applied_metadata <- ecoda_hpc_apply_batch_context(
  breast_sample_meta[, c("Sample", "disease"), drop = FALSE],
  breast_context,
  sample_col = "Sample"
)
stopifnot(
  identical(breast_validation$effective_batch_keys, breast_keys),
  identical(breast_validation$non_estimable_batch_keys, character()),
  identical(breast_validation$fixed_effect_aliases, c(
    assay = "batch_key_1",
    sequencing_platform = "batch_key_2",
    suspension_dissociation_time = "batch_key_3"
  )),
  identical(breast_design$rank, 5L),
  identical(breast_design$columns, 5L),
  identical(breast_design$residual_df, 7L),
  identical(
    colnames(breast_applied_metadata),
    c("Sample", "disease", breast_keys)
  ),
  identical(
    as.character(breast_applied_metadata$assay),
    as.character(breast_sample_meta$assay)
  ),
  identical(
    as.character(breast_applied_metadata$sequencing_platform),
    as.character(breast_sample_meta$sequencing_platform)
  ),
  identical(
    as.character(breast_applied_metadata$suspension_dissociation_time),
    as.character(breast_sample_meta$suspension_dissociation_time)
  ),
  !("__ecoda_batch_combined_v1" %in% colnames(breast_applied_metadata)),
  all(is.finite(breast_corrected)),
  identical(dimnames(breast_corrected), dimnames(breast_feat)),
  max(abs(rowSums(breast_corrected))) < 1e-8,
  identical(breast_contract$model_id, "limma_fixed_effects_v1"),
  identical(breast_contract$correction_mode, "limma_fixed_effects"),
  identical(breast_contract$design_rank, 5L),
  identical(breast_contract$design_columns, 5L),
  identical(breast_contract$design_residual_df, 7L),
  grepl(
    "model.matrix(~ 1 + batch_key_1 + batch_key_2 + batch_key_3); limma::removeBatchEffect(covariates=technical_covariates, design=intercept); effective_batch_keys=[assay,sequencing_platform,suspension_dissociation_time]",
    breast_contract$correction_formula,
    fixed = TRUE
  ),
  !grepl(
    "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
    breast_contract$correction_formula,
    perl = TRUE
  )
)

breast_pb_matrix <- matrix(
  rnorm(3L * length(breast_samples)),
  nrow = 3L,
  dimnames = list(paste0("breast-g", 1:3), breast_samples)
)
breast_pb <- ecoda_hpc_apply_limma_batch_correction(
  fit = list(norm_matrix = breast_pb_matrix),
  metadata = breast_sample_meta,
  batch_keys = as.list(breast_keys),
  sample_col = "Sample"
)
stopifnot(
  isTRUE(breast_pb$correct_batch),
  is.null(breast_pb$batch_col),
  all(is.finite(breast_pb$norm_matrix)),
  identical(dimnames(breast_pb$norm_matrix), dimnames(breast_pb_matrix)),
  identical(
    breast_pb$batch_correction$configured_batch_keys,
    breast_keys
  ),
  identical(
    breast_pb$batch_correction$effective_batch_keys,
    breast_keys
  ),
  identical(
    breast_pb$batch_correction$correction_mode,
    "limma_fixed_effects_pseudobulk"
  ),
  identical(breast_pb$batch_correction$design_rank, 5L),
  identical(breast_pb$batch_correction$design_columns, 5L),
  identical(breast_pb$batch_correction$design_residual_df, 7L),
  grepl(
    "DESeq2 design=~ 1; model.matrix(~ 1 + batch_key_1 + batch_key_2 + batch_key_3); limma::removeBatchEffect(covariates=technical_covariates, design=intercept); effective_batch_keys=[assay,sequencing_platform,suspension_dissociation_time]",
    breast_pb$batch_correction$correction_formula,
    fixed = TRUE
  ),
  !grepl(
    "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
    breast_pb$batch_correction$correction_formula,
    perl = TRUE
  )
)

# Exercise the composition driver itself for a one-key and a three-key
# corrected run. The fixture process emits deterministic CLR features, so the
# observable contract is the corrected matrix, identifiers, and bundle policy.
run_corrected_composition_fixture <- function(
  sample_metadata,
  cell_metadata,
  batch_keys,
  validation,
  result_stem
) {
  old_process <- process_coda_fig
  old_save <- save_rds_atomic
  old_exec <- exec_time
  old_peak <- peak_rss_gb
  old_log <- log_exec_row
  on.exit({
    process_coda_fig <<- old_process
    save_rds_atomic <<- old_save
    exec_time <<- old_exec
    peak_rss_gb <<- old_peak
    log_exec_row <<- old_log
  }, add = TRUE)

  sample_ids <- as.character(sample_metadata[["Sample"]])
  features <- matrix(
    rnorm(length(sample_ids) * 3L),
    nrow = length(sample_ids),
    dimnames = list(sample_ids, paste0("composition-ct", 1:3))
  )
  features <- features - rowMeans(features)
  process_coda_fig <<- local({
    template <- features
    function(seurat, labels, ...) {
      ids <- names(labels)
      feature <- template[ids, , drop = FALSE]
      list(
        scores = list(sil_score = 0.5),
        feat_mat = feature,
        dist_mat = dist(feature),
        labels = labels
      )
    }
  })
  save_rds_atomic <<- function(...) invisible(NULL)
  exec_time <<- function(expr) {
    force(expr)
    0
  }
  peak_rss_gb <<- function() 0
  log_exec_row <<- function(...) invisible(NULL)

  if (!"label" %in% colnames(sample_metadata)) {
    sample_metadata[["label"]] <- factor(
      rep(c("A", "B"), length.out = nrow(sample_metadata))
    )
  }
  labels <- structure(
    factor(rep(c("A", "B"), length.out = length(sample_ids))),
    names = sample_ids
  )
  cell_metadata[["author_cell_type"]] <- rep(
    c("T", "B"),
    length.out = nrow(cell_metadata)
  )
  cell_metadata[["RNA_snn_res.2"]] <- rep(
    c("0", "1"),
    length.out = nrow(cell_metadata)
  )
  identity <- ecoda_hpc_batch_contract_identity(
    batch_keys = as.list(unname(batch_keys)),
    sample_col = "Sample",
    method_id = "ECODA_authors_HR",
    model_id = "limma_fixed_effects_v1"
  )
  contract <- ecoda_hpc_augment_batch_contract(
    identity = identity,
    validation = validation,
    method_id = "ECODA_authors_HR",
    batch_keys = as.list(unname(batch_keys))
  )
  run_composition_methods_hpc(
    labels = labels,
    metadata = sample_metadata,
    pca_emb = NULL,
    pb_hvg2000 = NULL,
    obs = cell_metadata,
    label_col = "label",
    ct_col_high_res = "author_cell_type",
    sample_col = "Sample",
    results_dir = tempfile("composition-limma-"),
    ds = "Synthetic",
    batch_mode = TRUE,
    result_stem = result_stem,
    corrected = TRUE,
    batch_keys = as.list(unname(batch_keys)),
    metadata_validation = validation,
    batch_contract = contract
  )
}

one_key_composition <- run_corrected_composition_fixture(
  sample_metadata = meta,
  cell_metadata = meta[rep(seq_len(nrow(meta)), each = 2L), , drop = FALSE],
  batch_keys = "tech",
  validation = one_key_validation,
  result_stem = "Synthetic_batch_effect_uncorrected_one_key"
)
stopifnot(
  identical(
    names(one_key_composition),
    c("ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2", "batch_contract")
  ),
  identical(
    one_key_composition$batch_contract$model_id,
    "limma_fixed_effects_v1"
  ),
  identical(
    one_key_composition$batch_contract$effective_batch_keys,
    "tech"
  ),
  all(vapply(
    one_key_composition[c(
      "ECODA_authors_HR",
      "ECODA_authors_HR_NULL",
      "ECODA_seuratres_2"
    )],
    function(result) {
      all(is.finite(result$feat_mat)) &&
        identical(rownames(result$feat_mat), rownames(feat)) &&
        max(abs(rowSums(result$feat_mat))) < 1e-8 &&
        !grepl(
          "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
          result$batch_contract$correction_formula,
          perl = TRUE
        )
    },
    logical(1)
  ))
)

breast_composition <- run_corrected_composition_fixture(
  sample_metadata = breast_sample_meta,
  cell_metadata = breast_cell_meta,
  batch_keys = breast_keys,
  validation = breast_validation,
  result_stem = "Synthetic_breast_batch_effect_uncorrected"
)
stopifnot(
  identical(
    names(breast_composition),
    c("ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2", "batch_contract")
  ),
  identical(
    breast_composition$batch_contract$effective_batch_keys,
    breast_keys
  ),
  all(vapply(
    breast_composition[c(
      "ECODA_authors_HR",
      "ECODA_authors_HR_NULL",
      "ECODA_seuratres_2"
    )],
    function(result) {
      all(is.finite(result$feat_mat)) &&
        identical(rownames(result$feat_mat), breast_samples) &&
        identical(names(result$labels), breast_samples) &&
        max(abs(rowSums(result$feat_mat))) < 1e-8 &&
        grepl(
          "model.matrix(~ 1 + batch_key_1 + batch_key_2 + batch_key_3); limma::removeBatchEffect",
          result$batch_contract$correction_formula,
          fixed = TRUE
        ) &&
        !grepl(
          "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
          result$batch_contract$correction_formula,
          perl = TRUE
        )
    },
    logical(1)
  ))
)
# All configured technical keys may be constant. The full-cell validator still
# validates every row, while both corrected consumers retain an exact no-op.
all_constant_meta <- data.frame(
  Sample = rownames(feat),
  assay = factor(rep("10x 3' v3", n)),
  sex = factor(rep("F", n)),
  stringsAsFactors = FALSE
)
all_constant_cells <- all_constant_meta[
  rep(seq_len(nrow(all_constant_meta)), each = 2L),
  ,
  drop = FALSE
]
rownames(all_constant_cells) <- paste0("constant-cell-", seq_len(nrow(all_constant_cells)))
all_constant_validation <- ecoda_batch_validate_metadata(
  all_constant_cells,
  batch_keys = list("assay", "sex"),
  sample_col = "Sample"
)
all_constant_corrected <- correct_clr_batch_limma(
  feat,
  all_constant_meta,
  batch_keys = list("assay", "sex"),
  sample_col = "Sample",
  metadata_validation = all_constant_validation
)
all_constant_contract <- ecoda_hpc_augment_batch_contract(
  identity = multikey_identity,
  validation = all_constant_validation,
  method_id = "ECODA_authors_HR",
  batch_keys = list("assay", "sex")
)
stopifnot(
  identical(all_constant_validation$correction_state, "NO_CORRECTION"),
  identical(all_constant_validation$effective_batch_keys, character()),
  identical(all_constant_corrected, feat),
  identical(all_constant_contract$effective_batch_keys, character()),
  identical(all_constant_contract$correction_state, "NO_CORRECTION"),
  identical(all_constant_contract$correction_mode, "limma_fixed_effects"),
  identical(all_constant_contract$fixed_effect_aliases, character()),
  identical(all_constant_contract$correction_design_formula, "~1"),
  identical(all_constant_contract$design_rank, 1L),
  identical(all_constant_contract$design_columns, 1L),
  all_constant_contract$design_residual_df > 0L,
  identical(
    all_constant_contract$correction_formula,
    "NO_CORRECTION: no estimable technical batch key"
  ),
  !"__ecoda_batch_combined_v1" %in% colnames(all_constant_validation$sample_metadata)
)
ecoda_hpc_validate_batch_contract(
  all_constant_contract,
  all_constant_contract,
  label = "all-constant limma no-op",
  require_effective_metadata = TRUE
)

# Full-cell strictness remains mandatory and catches a within-Sample
# disagreement or a pre-existing reserved temporary column.
inconsistent_cells <- all_constant_cells
inconsistent_cells$sex <- as.character(inconsistent_cells$sex)
inconsistent_cells$sex[2L] <- "M"
strict_full_cell_error <- tryCatch(
  ecoda_batch_validate_metadata(
    inconsistent_cells,
    batch_keys = list("assay", "sex"),
    sample_col = "Sample"
  ),
  error = identity
)
reserved_cells <- all_constant_cells
reserved_cells[["__ecoda_batch_combined_v1"]] <- "forbidden"
reserved_column_error <- tryCatch(
  ecoda_batch_validate_metadata(
    reserved_cells,
    batch_keys = list("assay", "sex"),
    sample_col = "Sample"
  ),
  error = identity
)
stopifnot(
  inherits(strict_full_cell_error, "error"),
  grepl("disagrees within Sample", conditionMessage(strict_full_cell_error), fixed = TRUE),
  inherits(reserved_column_error, "error"),
  grepl("__ecoda_batch_combined_v1", conditionMessage(reserved_column_error), fixed = TRUE)
)

# A rank-deficient additive design and a saturated design with no residual
# degrees of freedom are both rejected before any limma call.
rank_deficient_meta <- data.frame(
  Sample = paste0("rank-s", 1:4),
  assay = factor(c("rna", "rna", "atac", "atac")),
  sex = factor(c("F", "F", "M", "M")),
  stringsAsFactors = FALSE
)
rank_error <- tryCatch(
  ecoda_batch_validate_metadata(
    rank_deficient_meta,
    batch_keys = list("assay", "sex"),
    sample_col = "Sample"
  ),
  error = identity
)
zero_residual_meta <- data.frame(
  Sample = paste0("df-s", 1:4),
  assay = factor(c("rna", "atac", "rna", "rna")),
  sex = factor(c("F", "F", "M", "F")),
  sequencing_platform = factor(c("NovaSeq", "NovaSeq", "NovaSeq", "Illumina")),
  stringsAsFactors = FALSE
)
residual_df_error <- tryCatch(
  ecoda_batch_validate_metadata(
    zero_residual_meta,
    batch_keys = list("assay", "sex", "sequencing_platform"),
    sample_col = "Sample"
  ),
  error = identity
)
stopifnot(
  inherits(rank_error, "error"),
  grepl("rank deficient|confounded", conditionMessage(rank_error), ignore.case = TRUE),
  inherits(residual_df_error, "error"),
  grepl("residual degrees of freedom", conditionMessage(residual_df_error), fixed = TRUE)
)
# Direct pseudobulk correction uses one, two, or three original technical
# factors. The fixture includes categorical suspension-like values and a
# biological label that must never enter the corrected fit.
utility_pb_samples <- paste0("utility-pb-s", seq_len(12L))
utility_pb_counts <- outer(
  seq_len(8L),
  seq_len(12L),
  FUN = function(gene, sample) as.integer(20L + 3L * gene + sample + (sample %% 3L) * gene)
)
dimnames(utility_pb_counts) <- list(
  paste0("utility-pb-g", seq_len(nrow(utility_pb_counts))),
  utility_pb_samples
)
utility_pb_metadata <- data.frame(
  Sample = utility_pb_samples,
  assay = factor(rep(c("10x 3' v3", "10x multiome"), each = 6L)),
  platform = factor(rep(c("NovaSeq", "NextSeq"), times = 6L)),
  suspension_dissociation_time = factor(
    rep(c("fresh 0 min", "frozen 30 min", "ambient 2 h"), each = 4L)
  ),
  biological_label = factor(rep(c("T cell", "B cell"), times = 6L)),
  row.names = utility_pb_samples,
  stringsAsFactors = FALSE
)
utility_pb_one <- fit_pseudobulk_deseq2(
  utility_pb_counts,
  utility_pb_metadata,
  batch_col = "assay",
  blind = FALSE,
  correct_batch = TRUE
)
utility_pb_two <- fit_pseudobulk_deseq2(
  utility_pb_counts,
  utility_pb_metadata,
  batch_col = c("assay", "platform"),
  blind = FALSE,
  correct_batch = TRUE
)
utility_pb_three <- fit_pseudobulk_deseq2(
  utility_pb_counts,
  utility_pb_metadata,
  batch_col = c("assay", "platform", "suspension_dissociation_time"),
  blind = FALSE,
  correct_batch = TRUE
)
stopifnot(
  identical(
    utility_pb_one$batch_correction$configured_batch_keys,
    "assay"
  ),
  identical(
    utility_pb_two$batch_correction$configured_batch_keys,
    c("assay", "platform")
  ),
  identical(
    utility_pb_three$batch_correction$configured_batch_keys,
    c("assay", "platform", "suspension_dissociation_time")
  ),
  identical(
    utility_pb_three$batch_correction$effective_batch_keys,
    c("assay", "platform", "suspension_dissociation_time")
  ),
  identical(utility_pb_one$batch_correction$design_columns, 2L),
  identical(utility_pb_two$batch_correction$design_columns, 3L),
  identical(utility_pb_three$batch_correction$design_columns, 5L),
  is.factor(utility_pb_three$metadata[["assay"]]),
  is.factor(utility_pb_three$metadata[["platform"]]),
  is.factor(utility_pb_three$metadata[["suspension_dissociation_time"]]),
  !"__ecoda_batch_combined_v1" %in% colnames(utility_pb_three$metadata),
  !"biological_label" %in% colnames(utility_pb_three$metadata),
  all(is.finite(utility_pb_one$norm_matrix)),
  all(is.finite(utility_pb_two$norm_matrix)),
  all(is.finite(utility_pb_three$norm_matrix)),
  identical(dimnames(utility_pb_three$norm_matrix), dimnames(utility_pb_counts)),
  grepl(
    "DESeq2 design=~ 1; model.matrix(~ 1 + batch_key_1 + batch_key_2 + batch_key_3)",
    utility_pb_three$batch_correction$correction_formula,
    fixed = TRUE
  ),
  !grepl(
    "__ecoda_batch_combined_v1|biological_label|label|Status|disease",
    utility_pb_three$batch_correction$correction_formula,
    ignore.case = TRUE,
    perl = TRUE
  )
)
# The utility's output is the documented limma fixed-effect operation on the
# uncorrected intercept-only normalized matrix, with separate dummy columns
# and an intercept-only protected design.
utility_pb_reference <- fit_pseudobulk_deseq2(
  utility_pb_counts,
  utility_pb_metadata,
  blind = FALSE,
  correct_batch = FALSE
)
utility_pb_reference_data <- data.frame(
  batch_key_1 = factor(
    as.character(utility_pb_metadata$assay),
    levels = unique(as.character(utility_pb_metadata$assay))
  ),
  batch_key_2 = factor(
    as.character(utility_pb_metadata$platform),
    levels = unique(as.character(utility_pb_metadata$platform))
  ),
  batch_key_3 = factor(
    as.character(utility_pb_metadata$suspension_dissociation_time),
    levels = unique(as.character(utility_pb_metadata$suspension_dissociation_time))
  ),
  row.names = utility_pb_samples
)
utility_pb_reference_design <- model.matrix(
  ~ 1 + batch_key_1 + batch_key_2 + batch_key_3,
  data = utility_pb_reference_data
)
utility_pb_reference_expected <- limma::removeBatchEffect(
  x = utility_pb_reference$norm_matrix,
  covariates = utility_pb_reference_design[, -1L, drop = FALSE],
  design = matrix(
    1,
    nrow = nrow(utility_pb_reference_design),
    ncol = 1L,
    dimnames = list(utility_pb_samples, "(Intercept)")
  )
)
dimnames(utility_pb_reference_expected) <- dimnames(utility_pb_counts)
stopifnot(isTRUE(all.equal(
  utility_pb_three$norm_matrix,
  utility_pb_reference_expected,
  tolerance = 1e-8
)))


# DESeq2.normalize retains the public scalar batch_col argument while exposing
# the same corrected genes-by-samples orientation and identifiers.
utility_pb_normalized <- DESeq2.normalize(
  utility_pb_counts,
  utility_pb_metadata,
  n_hvg = 4L,
  batch_col = c("assay", "platform", "suspension_dissociation_time"),
  blind = FALSE,
  correct_batch = TRUE
)
stopifnot(
  is.matrix(utility_pb_normalized),
  identical(
    dimnames(utility_pb_normalized),
    list(
      utility_pb_three$variance_order[seq_len(4L)],
      utility_pb_samples
    )
  ),
  all(is.finite(utility_pb_normalized))
)

# A biological-label permutation cannot affect the technical-only corrected
# fit. This also guards against accidental use of all metadata columns.
utility_pb_label_permuted <- utility_pb_metadata
utility_pb_label_permuted$biological_label <- factor(
  rev(as.character(utility_pb_label_permuted$biological_label)),
  levels = levels(utility_pb_metadata$biological_label)
)
utility_pb_three_label_permuted <- fit_pseudobulk_deseq2(
  utility_pb_counts,
  utility_pb_label_permuted,
  batch_col = c("assay", "platform", "suspension_dissociation_time"),
  blind = FALSE,
  correct_batch = TRUE
)
stopifnot(isTRUE(all.equal(
  utility_pb_three$norm_matrix,
  utility_pb_three_label_permuted$norm_matrix,
  tolerance = 1e-8
)))

# All-constant technical keys are an exact no-op after the same intercept-only
# DESeq2 normalization.
utility_pb_constant_metadata <- utility_pb_metadata
utility_pb_constant_metadata$assay <- factor(rep("single assay", 12L))
utility_pb_constant_metadata$platform <- factor(rep("single platform", 12L))
utility_pb_unbatched <- utility_pb_reference
utility_pb_constant <- fit_pseudobulk_deseq2(
  utility_pb_counts,
  utility_pb_constant_metadata,
  batch_col = c("assay", "platform"),
  blind = FALSE,
  correct_batch = TRUE
)
stopifnot(
  identical(
    utility_pb_constant$batch_correction$effective_batch_keys,
    character()
  ),
  isTRUE(all.equal(
    utility_pb_constant$norm_matrix,
    utility_pb_unbatched$norm_matrix,
    tolerance = 1e-8
  )),
  identical(
    dimnames(utility_pb_constant$norm_matrix),
    dimnames(utility_pb_counts)
  )
)

utility_pb_expect_error <- function(value, pattern) {
  captured_error <- tryCatch(value, error = identity)
  stopifnot(
    inherits(captured_error, "error"),
    grepl(pattern, conditionMessage(captured_error), ignore.case = TRUE)
  )
  invisible(TRUE)
}
# Breast's configured suspension-duration key accepts the literal ``unknown``
# category, while the same value remains invalid for unrelated technical keys.
utility_pb_suspension_unknown <- .pseudobulk_factor_column(
  c("fresh 0 min", "unknown", "ambient 2 h"),
  "suspension_dissociation_time"
)
stopifnot(
  is.factor(utility_pb_suspension_unknown),
  identical(
    levels(utility_pb_suspension_unknown),
    c("fresh 0 min", "unknown", "ambient 2 h")
  ),
  identical(
    as.character(utility_pb_suspension_unknown),
    c("fresh 0 min", "unknown", "ambient 2 h")
  )
)
utility_pb_expect_error(
  .pseudobulk_factor_column(c("NovaSeq", "unknown"), "assay"),
  "missing or blank"
)
utility_pb_expect_error(
  .pseudobulk_factor_column(c("NovaSeq", "unknown"), "platform"),
  "missing or blank"
)
utility_pb_expect_error(
  .pseudobulk_factor_column(
    c("fresh 0 min", "Unknown"),
    "suspension_dissociation_time"
  ),
  "missing or blank"
)
utility_pb_expect_error(
  .pseudobulk_factor_column(
    c("fresh 0 min", " unknown "),
    "suspension_dissociation_time"
  ),
  "missing or blank"
)
utility_pb_expect_error(
  .pseudobulk_factor_column(
    c("fresh 0 min", NA_character_),
    "suspension_dissociation_time"
  ),
  "missing"
)
utility_pb_expect_error(
  .pseudobulk_factor_column(
    c(1, Inf),
    "suspension_dissociation_time"
  ),
  "non-finite"
)

utility_pb_rank_deficient_metadata <- utility_pb_metadata
utility_pb_rank_deficient_metadata$platform <- utility_pb_rank_deficient_metadata$assay
utility_pb_expect_error(
  fit_pseudobulk_deseq2(
    utility_pb_counts,
    utility_pb_rank_deficient_metadata,
    batch_col = c("assay", "platform"),
    blind = FALSE,
    correct_batch = TRUE
  ),
  "rank deficient"
)
utility_pb_saturated_metadata <- utility_pb_metadata
utility_pb_saturated_metadata$unique_batch <- factor(utility_pb_samples)
utility_pb_expect_error(
  fit_pseudobulk_deseq2(
    utility_pb_counts,
    utility_pb_saturated_metadata,
    batch_col = "unique_batch",
    blind = FALSE,
    correct_batch = TRUE
  ),
  "residual degrees of freedom"
)
utility_pb_blank_metadata <- utility_pb_metadata
utility_pb_blank_metadata$assay <- as.character(utility_pb_blank_metadata$assay)
utility_pb_blank_metadata$assay[[1L]] <- " "
utility_pb_expect_error(
  fit_pseudobulk_deseq2(
    utility_pb_counts,
    utility_pb_blank_metadata,
    batch_col = "assay",
    blind = FALSE,
    correct_batch = TRUE
  ),
  "blank"
)
utility_pb_expect_error(
  fit_pseudobulk_deseq2(
    utility_pb_counts,
    utility_pb_metadata,
    batch_col = "__ecoda_batch_combined_v1",
    blind = FALSE,
    correct_batch = TRUE
  ),
  "reserved"
)


# The active composition correction is limma-only and never advertises the
# retired lme4/random-intercept or combined-key model.
correction_body <- paste(deparse(body(correct_clr_batch_limma)), collapse = " ")
stopifnot(
  grepl("limma::removeBatchEffect", correction_body, fixed = TRUE),
  !grepl("lme4|lmer|lmerTest|\\(1 \\|", correction_body, perl = TRUE),
  !grepl("__ecoda_batch_combined_v1", correction_body, fixed = TRUE),
  !grepl("label|Status|disease", correction_body, ignore.case = TRUE)
)

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
