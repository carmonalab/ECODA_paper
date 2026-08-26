# Focused regression checks for the Bassez metadata patch and benchmark loading.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])))
project_root <- dirname(script_dir)

assert_error <- function(expr, pattern) {
  condition <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  if (is.null(condition)) {
    stop("Expected an error matching /", pattern, "/.")
  }
  if (!grepl(pattern, condition$message)) {
    stop(
      "Error did not match /", pattern, "/: ",
      condition$message
    )
  }
  invisible(TRUE)
}

utils_env <- new.env(parent = globalenv())
sys.source(
  file.path(
    project_root,
    "src",
    "2_dataset_specific_preprocessing",
    "bassez_cellsubtype_utils.R"
  ),
  envir = utils_env
)

sentinel_values <- c(NA, "", "  ", "NA", " nan ", "None", " Unknown ", "valid")
stopifnot(identical(
  utils_env$bassez_missing_annotation(sentinel_values),
  c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
))
stopifnot(identical(
  utils_env$bassez_missing_annotation(factor(sentinel_values)),
  c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
))

metadata <- data.frame(
  cellType = c("T cell", "B cell", "Myeloid", "Cancer", "Endothelial", "Mast", "NK"),
  cellSubType = factor(c(
    "  preserved subtype  ", NA, "NA", " nan ", "None", "Unknown", "another subtype"
  )),
  row.names = paste0("cell", seq_len(7)),
  stringsAsFactors = FALSE
)
filled <- utils_env$bassez_fill_cell_subtype(metadata)
stopifnot(is.factor(filled$cellSubType))
stopifnot(identical(
  as.character(filled$cellSubType),
  c(
    "  preserved subtype  ", "B cell", "Myeloid", "Cancer",
    "Endothelial", "Mast", "another subtype"
  )
))
stopifnot(identical(as.character(filled$cellType), as.character(metadata$cellType)))

invalid_metadata <- data.frame(
  cellType = factor(c(" Unknown ", "B cell")),
  cellSubType = factor(c(NA, "preserved")),
  row.names = c("invalid_fallback", "valid_fallback"),
  stringsAsFactors = FALSE
)
assert_error(
  utils_env$bassez_fill_cell_subtype(invalid_metadata),
  "invalid cellType fallback"
)
assert_error(
  utils_env$bassez_fill_cell_subtype(data.frame(cellSubType = factor(NA))),
  "missing required column"
)
pseudobulk_env <- new.env(parent = globalenv())
pseudobulk_env$standardize_sample_names <- function(sample_names) {
  gsub("-", "_", as.character(sample_names), fixed = TRUE)
}
sys.source(
  file.path(project_root, "src", "utils", "pseudobulk.R"),
  envir = pseudobulk_env
)
pseudobulk_fixture <- matrix(
  c(1, 2, 3, 4),
  nrow = 2,
  dimnames = list(c("BIOKEY-2-Pre", "exact"), c("gene1", "gene2"))
)
aligned_pseudobulk <- pseudobulk_env$align_pseudobulk_sample_names(
  pseudobulk_fixture,
  c("BIOKEY_2_Pre", "exact")
)
stopifnot(identical(
  rownames(aligned_pseudobulk),
  c("BIOKEY_2_Pre", "exact")
))
assert_error(
  pseudobulk_env$align_pseudobulk_sample_names(
    pseudobulk_fixture,
    c("BIOKEY_2_Pre", "different")
  ),
  "do not match canonical metadata IDs"
)

pipeline_env <- new.env(parent = globalenv())
sys.source(
  file.path(
    project_root,
    "src",
    "5_run_benchmark_methods",
    "benchmark_pipeline.R"
  ),
  envir = pipeline_env
)

hpc_env <- new.env(parent = globalenv())
sys.source(
  file.path(
    project_root,
    "src",
    "5_run_benchmark_methods",
    "benchmark_hpc_utils.R"
  ),
  envir = hpc_env
)

pseudobulk_dir <- tempfile("ecoda_pseudobulk_cache-")
dir.create(pseudobulk_dir, recursive = TRUE)
for (variant in hpc_env$PB_VARIANT_NAMES) {
  saveRDS(
    list(pb = matrix(1, nrow = 1, ncol = 1), time_secs = 0),
    file.path(pseudobulk_dir, paste0("Toy_pseudobulk_", variant, ".rds"))
  )
}
stopifnot(identical(
  hpc_env$pb_variants_missing(pseudobulk_dir, "Toy", force = FALSE),
  character(0)
))
stopifnot(identical(
  hpc_env$pb_variants_missing(pseudobulk_dir, "Toy", force = TRUE),
  hpc_env$PB_VARIANT_NAMES
))
composition_pseudobulks <- hpc_env$load_composition_pb_variants(
  sample_col = "Sample",
  hvg_rank_genes = character(0),
  pseudobulk_dir = pseudobulk_dir,
  ds = "Toy"
)
stopifnot(identical(
  names(composition_pseudobulks),
  hpc_env$PB_VARIANT_NAMES
))
composition_loader_calls <- new.env(parent = emptyenv())
composition_loader <- function(
  seurat,
  sample_col,
  hvg_rank_genes,
  pseudobulk_dir,
  ds,
  force = FALSE,
  log_file = NULL,
  cache_stem = ds,
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE,
  variants = hpc_env$PB_VARIANT_NAMES
) {
  composition_loader_calls$seurat <- seurat
  composition_loader_calls$force <- force
  composition_loader_calls$args <- list(
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    pseudobulk_dir = pseudobulk_dir,
    ds = ds,
    log_file = log_file,
    cache_stem = cache_stem,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch,
    variants = variants
  )
  setNames(
    lapply(variants, function(variant) list(variant = variant)),
    variants
  )
}
injected_composition_pseudobulks <- hpc_env$load_composition_pb_variants(
  sample_col = "Sample",
  hvg_rank_genes = character(0),
  pseudobulk_dir = pseudobulk_dir,
  ds = "Toy",
  log_file = "composition-test.log",
  loader = composition_loader
)
stopifnot(is.null(composition_loader_calls$seurat))
stopifnot(identical(composition_loader_calls$force, FALSE))
stopifnot(identical(
  composition_loader_calls$args,
  list(
    sample_col = "Sample",
    hvg_rank_genes = character(0),
    pseudobulk_dir = pseudobulk_dir,
    ds = "Toy",
    log_file = "composition-test.log",
    cache_stem = "Toy",
    batch_col = NULL,
    blind = TRUE,
    correct_batch = FALSE,
    variants = hpc_env$PB_VARIANT_NAMES
  )
))
stopifnot(identical(
  names(injected_composition_pseudobulks),
  hpc_env$PB_VARIANT_NAMES
))
stopifnot(identical(
  unname(vapply(
    injected_composition_pseudobulks,
    function(variant) variant$variant,
    character(1)
  )),
  hpc_env$PB_VARIANT_NAMES
))
unlink(pseudobulk_dir, recursive = TRUE, force = TRUE)

# `--force` still invalidates composition result bundles, while the
# obs-only pseudobulk loader reuses the prepared cache above.
pipeline_env$peak_rss_gb <- function() NA_real_
pipeline_env$save_rds_atomic <- function(object, file) saveRDS(object, file)
pipeline_env$log_exec_row <- function(...) invisible(NULL)
pipeline_env$process_avg_pca_embedding_fig <- function(...) {
  list(marker = "fresh")
}
pipeline_env$process_deconv_fig <- function(...) {
  list(marker = "fresh")
}
pipeline_env$process_coda_fig <- function(...) {
  list(marker = "fresh")
}


composition_results_dir <- tempfile("ecoda_composition_results-")
dir.create(composition_results_dir, recursive = TRUE)
composition_obs <- data.frame(
  Sample = c("sample_1", "sample_2"),
  stringsAsFactors = FALSE
)
composition_labels <- factor(c("group_1", "group_2"))
names(composition_labels) <- composition_obs$Sample
composition_metadata <- composition_obs
composition_pca <- matrix(1, nrow = 2, ncol = 1)
composition_pb <- list(pb = matrix(1, nrow = 1, ncol = 2))
run_composition <- function(force) {
  pipeline_env$run_composition_methods_hpc(
    labels = composition_labels,
    metadata = composition_metadata,
    pca_emb = composition_pca,
    pb_hvg2000 = composition_pb,
    obs = composition_obs,
    label_col = "label",
    sample_col = "Sample",
    results_dir = composition_results_dir,
    ds = "Toy",
    force = force,
    factors_test = integer(0),
    seurat_res = 0.1,
    ECODA_top_varexp_hvct = numeric(0)
  )
}
run_composition(force = FALSE)
composition_bundle <- file.path(
  composition_results_dir,
  "Toy_Avg_PCA_embedding.rds"
)
saveRDS(
  list(marker = "cached", exec_time = 0, mem_GB = NA_real_),
  composition_bundle
)
cached_composition <- run_composition(force = FALSE)
stopifnot(identical(cached_composition$Avg_PCA_embedding$marker, "cached"))
forced_composition <- run_composition(force = TRUE)
stopifnot(identical(forced_composition$Avg_PCA_embedding$marker, "fresh"))
unlink(composition_results_dir, recursive = TRUE, force = TRUE)



fixture_root <- tempfile("ecoda_checksum_test-")
results_dir <- file.path(fixture_root, "benchmark", "results")
dir.create(results_dir, recursive = TRUE)
on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)

composition_file <- file.path(results_dir, "Toy_composition.rds")
trans_file <- file.path(results_dir, "Toy_trans.rds")
zeroimp_file <- file.path(results_dir, "Toy_zeroimp.rds")
saveRDS(list(example = 42), composition_file)
saveRDS(list(transformed = TRUE), trans_file)
saveRDS(list(imputed = TRUE), zeroimp_file)

hash_for <- function(path) unname(tools::md5sum(path))
checksum_file <- file.path(fixture_root, "benchmark", "checksums.md5")
writeLines(c(
  paste0(hash_for(composition_file), "  results/Toy_composition.rds"),
  paste0(hash_for(trans_file), "  results/Toy_trans.rds"),
  paste0(hash_for(zeroimp_file), "  results/Toy_zeroimp.rds")
), checksum_file)

checksum_lines <- readLines(checksum_file)
checksum_lines[1] <- paste0(strrep("0", 32), "  results/Toy_composition.rds")
writeLines(checksum_lines, checksum_file)
assert_error(
  pipeline_env$load_hpc_benchmark_results(
    list(), "Toy", results_dir, methods = "composition"
  ),
  "Checksum mismatch"
)

writeLines(c(
  paste0(hash_for(composition_file), "  results/Toy_composition.rds"),
  paste0(hash_for(trans_file), "  results/Toy_trans.rds"),
  paste0(hash_for(zeroimp_file), "  results/Toy_zeroimp.rds")
), checksum_file)
loaded <- pipeline_env$load_hpc_benchmark_results(
  list(), "Toy", results_dir, methods = "composition"
)
stopifnot(identical(loaded$bmark$Toy$example, 42))
stopifnot(isTRUE(loaded$trans$Toy$transformed))
stopifnot(isTRUE(loaded$zeroimp$Toy$imputed))

message("Bassez and benchmark regression checks passed.")
