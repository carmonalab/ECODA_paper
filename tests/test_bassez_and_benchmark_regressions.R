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
  c(10, 20, 30, 40),
  nrow = 2,
  dimnames = list(c("exact", "BIOKEY-2-Pre"), c("gene1", "gene2"))
)
aligned_pseudobulk <- pseudobulk_env$align_pseudobulk_sample_names(
  pseudobulk_fixture,
  c("BIOKEY_2_Pre", "exact")
)
stopifnot(identical(
  rownames(aligned_pseudobulk),
  c("BIOKEY_2_Pre", "exact")
))
stopifnot(identical(
  as.numeric(aligned_pseudobulk[1, ]),
  c(20, 40)
))
stopifnot(identical(
  as.numeric(aligned_pseudobulk[2, ]),
  c(10, 30)
))
assert_error(
  pseudobulk_env$align_pseudobulk_sample_names(
    pseudobulk_fixture,
    c("BIOKEY_2_Pre", "different")
  ),
  "do not match canonical metadata IDs"
)

# The direct matrix boundary rejects values that cannot be represented by
# DESeq2 before constructing a DESeqDataSet.
overflow_counts <- matrix(
  as.numeric(.Machine$integer.max) + 1,
  nrow = 1L,
  ncol = 2L,
  dimnames = list("overflow_gene", c("s1", "s2"))
)
overflow_metadata <- data.frame(
  Sample = c("s1", "s2"),
  row.names = c("s1", "s2"),
  stringsAsFactors = FALSE
)
assert_error(
  pseudobulk_env$get_pb_deseq2_from_counts(
    overflow_counts,
    overflow_metadata
  ),
  "exceed"
)

# Full-gene HVG variants are selections from one fit and therefore agree with
# independent legacy normalization for every shared full-gene prefix.
direct_counts <- matrix(
  c(
    30, 10, 40,
    20, 40, 50,
    7, 14, 21,
    1, 2, 3,
    1000, 1, 1
  ),
  nrow = 5L,
  byrow = TRUE,
  dimnames = list(paste0("g", 1:5), c("s2", "s1", "s3"))
)
direct_metadata <- data.frame(
  Sample = c("s1", "s2", "s3"),
  batch = factor(c("A", "B", "A")),
  row.names = c("s1", "s2", "s3"),
  stringsAsFactors = FALSE
)
shared_fit <- pseudobulk_env$fit_pseudobulk_deseq2(
  direct_counts,
  direct_metadata
)
stopifnot(
  identical(rownames(shared_fit$norm_matrix), rownames(direct_counts)),
  identical(colnames(shared_fit$norm_matrix), colnames(direct_counts)),
  length(shared_fit$variance_order) == nrow(direct_counts),
  setequal(shared_fit$variance_order, rownames(direct_counts))
)
assert_matrix_close <- function(actual, expected, tolerance = 1e-7) {
  if (!identical(dim(actual), dim(expected)) ||
      !identical(dimnames(actual), dimnames(expected))) {
    stop("normalized matrix dimensions or identifiers differ")
  }
  delta <- max(abs(as.numeric(actual) - as.numeric(expected)))
  if (!is.finite(delta) || delta > tolerance) {
    stop("normalized matrices differ by ", delta)
  }
  invisible(TRUE)
}
for (n_hvg in c(2L, 3L, 5L)) {
  shared_selected <- pseudobulk_env$select_pseudobulk_deseq2(
    shared_fit,
    n_hvg = n_hvg
  )
  independent_selected <- pseudobulk_env[["DESeq2.normalize"]](
    direct_counts,
    direct_metadata,
    n_hvg = n_hvg
  )
  assert_matrix_close(shared_selected, independent_selected)
}
direct_published <- pseudobulk_env$get_pb_deseq2_from_counts(
  direct_counts,
  direct_metadata,
  n_hvg = 3L
)
expected_published <- t(
  pseudobulk_env$select_pseudobulk_deseq2(shared_fit, n_hvg = 3L)
)
expected_published <- expected_published[
  c("s1", "s2", "s3"),
  ,
  drop = FALSE
]
assert_matrix_close(direct_published, expected_published)
stopifnot(isTRUE(all.equal(
  pseudobulk_env$select_pseudobulk_deseq2(
    shared_fit,
    n_hvg = 5L,
    black_list = "default_without_sex_genes"
  ),
  pseudobulk_env$select_pseudobulk_deseq2(
    shared_fit,
    n_hvg = 5L,
    black_list = "none"
  ),
  tolerance = 1e-7
)))

# schvg2000 is a separate fit on the prefiltered raw gene universe, rather
# than a prefix selected from the full-gene fit.
schvg_genes <- c("g1", "g2", "g3")
schvg_fit <- pseudobulk_env$fit_pseudobulk_deseq2(
  direct_counts[schvg_genes, , drop = FALSE],
  direct_metadata
)
schvg_selected <- pseudobulk_env$select_pseudobulk_deseq2(
  schvg_fit,
  n_hvg = 2000L
)
stopifnot(
  identical(rownames(schvg_fit$norm_matrix), schvg_genes),
  identical(rownames(schvg_selected), schvg_genes)
)
full_shared_for_schvg <- shared_fit$norm_matrix[
  schvg_genes,
  colnames(direct_counts),
  drop = FALSE
]
if (max(abs(schvg_fit$norm_matrix - full_shared_for_schvg)) <= 1e-7) {
  stop("schvg fit unexpectedly reused full-gene normalization")
}

suppressPackageStartupMessages(library(dplyr))
metadata_env <- new.env(parent = globalenv())
sys.source(
  file.path(project_root, "src", "utils", "seurat_utils.R"),
  envir = metadata_env
)
metadata_with_unused_level <- data.frame(
  Sample = factor(
    c("sample_keep", "sample_keep"),
    levels = c("sample_keep", "sample_removed")
  ),
  label = c("case", "case"),
  stringsAsFactors = FALSE
)
collapsed_metadata <- metadata_env$collapse_sample_metadata(
  metadata_with_unused_level
)
stopifnot(nrow(collapsed_metadata) == 1L)
stopifnot(identical(as.character(collapsed_metadata$Sample), "sample_keep"))
ordered_metadata <- metadata_env$collapse_sample_metadata(
  data.frame(
    Sample = factor(c("sample_2", "sample_1", "sample_2", "sample_1")),
    label = c("case", "control", "case", "control"),
    stringsAsFactors = FALSE
  )
)
stopifnot(identical(as.character(ordered_metadata$Sample), c("sample_2", "sample_1")))
assert_error(
  metadata_env$collapse_sample_metadata(
    data.frame(Sample = c("sample_keep", ""), stringsAsFactors = FALSE)
  ),
  "missing or blank sample IDs"
)
composition_obs <- data.frame(
  Sample = factor(
    c("sample_2", "sample_1", "sample_2", "sample_1"),
    levels = c("sample_1", "sample_2")
  ),
  ct = factor(c("B", "A", "A", "B")),
  stringsAsFactors = FALSE
)
composition_counts <- metadata_env$get_ct_comp_df(
  composition_obs, sample_col = "Sample", ct_col = "ct"
)
stopifnot(identical(rownames(composition_counts), c("sample_2", "sample_1")))


methods_env <- new.env(parent = globalenv())
sys.source(
  file.path(project_root, "src", "5_run_benchmark_methods", "benchmark_methods_r.R"),
  envir = methods_env
)
mofa_metadata <- methods_env$prepare_mofa_metadata(
  data.frame(
    Sample = factor(c("sample_1", "sample_2")),
    label = c("case", "control"),
    stringsAsFactors = FALSE
  )
)
stopifnot(identical(rownames(mofa_metadata), c("sample_1", "sample_2")))
stopifnot(identical(as.character(mofa_metadata$sample), c("sample_1", "sample_2")))
assert_error(
  methods_env$prepare_mofa_metadata(
    data.frame(Sample = c("sample_1", "sample_1"), stringsAsFactors = FALSE)
  ),
  "nonmissing, non-empty, and unique"
)
result_features <- matrix(
  c(10, 20, 30, 40),
  nrow = 2,
  dimnames = list(c("sample_2", "sample_1"), c("PC_1", "PC_2"))
)
result_labels <- factor(c("case", "control"))
names(result_labels) <- c("sample_1", "sample_2")
result_dist <- as.dist(matrix(
  c(0, 9, 9, 0),
  nrow = 2,
  dimnames = list(c("sample_2", "sample_1"), c("sample_2", "sample_1"))
))
aligned_result <- methods_env$align_result_samples(
  result_features, result_labels, result_dist
)
stopifnot(identical(
  rownames(aligned_result$feat_mat),
  c("sample_1", "sample_2")
))
stopifnot(identical(
  as.numeric(aligned_result$feat_mat[1, ]),
  c(20, 40)
))
stopifnot(identical(
  names(aligned_result$labels),
  c("sample_1", "sample_2")
))
stopifnot(identical(
  rownames(aligned_result$dist_mat),
  c("sample_1", "sample_2")
))
stopifnot(identical(
  colnames(aligned_result$dist_mat),
  c("sample_1", "sample_2")
))
square_features <- matrix(
  c(1, 2, 3, 4),
  nrow = 2,
  dimnames = list(c("sample_2", "sample_1"), c("sample_2", "sample_1"))
)
square_result <- methods_env$align_result_samples(
  square_features, result_labels
)
stopifnot(identical(
  rownames(square_result$feat_mat),
  c("sample_1", "sample_2")
))
stopifnot(identical(
  colnames(square_result$feat_mat),
  c("sample_1", "sample_2")
))
stopifnot(identical(
  as.numeric(square_result$feat_mat[1, ]),
  c(4, 2)
))



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

write_checked_fixture <- function(path, value) {
  saveRDS(value, path)
  writeLines(
    c(
      paste0("MD5=", unname(tools::md5sum(path))),
      paste0("SIZE=", file.info(path)$size),
      paste0("PATH=", path)
    ),
    paste0(path, ".md5")
  )
}
validator_path <- file.path(
  project_root,
  "src",
  "5_run_benchmark_methods",
  "validate_benchmark_rds_contract.R"
)
validator_status <- function(path) {
  system2(
    "pixi",
    c(
      "run", "Rscript", "--vanilla", validator_path,
      "--artifact", path,
      "--method", "prepare_pseudobulk"
    ),
    stdout = FALSE,
    stderr = FALSE
  )
}
expect_validator_ok <- function(path, label) {
  status <- validator_status(path)
  if (!identical(status, 0L)) {
    stop("expected pseudobulk validator success: ", label)
  }
  invisible(TRUE)
}
expect_validator_failure <- function(path, label) {
  status <- validator_status(path)
  if (identical(status, 0L)) {
    stop("expected pseudobulk validator failure: ", label)
  }
  invisible(TRUE)
}
schema2_record <- function(
  variant,
  timing_id = "run-cache:Toy:benchmark_analysis:none"
) {
  list(
    pb = matrix(
      match(variant, hpc_env$PB_VARIANT_NAMES),
      nrow = 1L,
      dimnames = list(c("sample_1"), c("gene_1"))
    ),
    time_secs = 0.25,
    mem_GB = NA_real_,
    aggregate_time_secs = 1.25,
    shared_fit_time_secs = 2.75,
    shared_time_secs = 4,
    variant_time_secs = 0.25,
    shared_mem_GB = NA_real_,
    timing_id = timing_id,
    timing_schema = 2L
  )
}
schema2_dir <- tempfile("ecoda_pseudobulk_schema2-")
dir.create(schema2_dir, recursive = TRUE)
schema2_paths <- setNames(
  file.path(
    schema2_dir,
    paste0("Toy_pseudobulk_", hpc_env$PB_VARIANT_NAMES, ".rds")
  ),
  hpc_env$PB_VARIANT_NAMES
)
for (variant in names(schema2_paths)) {
  write_checked_fixture(schema2_paths[[variant]], schema2_record(variant))
}
expect_validator_ok(schema2_paths[[1L]], "schema-2 cache")
legacy_matrix_path <- file.path(schema2_dir, "legacy-matrix.rds")
write_checked_fixture(
  legacy_matrix_path,
  matrix(1, nrow = 1L, ncol = 1L,
         dimnames = list("sample_1", "gene_1"))
)
expect_validator_ok(legacy_matrix_path, "legacy matrix cache")
legacy_wrapper_path <- file.path(schema2_dir, "legacy-wrapper.rds")
write_checked_fixture(
  legacy_wrapper_path,
  list(
    pb = matrix(1, nrow = 1L, ncol = 1L,
                dimnames = list("sample_1", "gene_1")),
    time_secs = 3.5,
    mem_GB = NA_real_
  )
)
expect_validator_ok(legacy_wrapper_path, "legacy wrapped cache")
stopifnot(identical(
  hpc_env$pb_variants_missing(schema2_dir, "Toy"),
  character(0)
))
cached_schema2 <- hpc_env$load_pb_variants(
  seurat = NULL,
  sample_col = "Sample",
  hvg_rank_genes = paste0("g", 1:5),
  pseudobulk_dir = schema2_dir,
  ds = "Toy",
  h5ad_path = file.path(schema2_dir, "must-not-be-read.h5ad")
)
stopifnot(identical(names(cached_schema2), hpc_env$PB_VARIANT_NAMES))
missing_schema2 <- schema2_paths[["hvg500"]]
unlink(c(missing_schema2, paste0(missing_schema2, ".md5")))
stopifnot(identical(
  hpc_env$pb_variants_missing(schema2_dir, "Toy"),
  "hvg500"
))
stopifnot(identical(
  hpc_env$pb_variants_missing(schema2_dir, "Toy", force = TRUE),
  hpc_env$PB_VARIANT_NAMES
))
original_prepare_pseudobulks <- hpc_env$prepare_pseudobulks_hpc
prepared_missing <- character()
hpc_env$prepare_pseudobulks_hpc <- function(h5ad_path, variants, ...) {
  prepared_missing <<- as.character(variants)
  setNames(lapply(variants, schema2_record), as.character(variants))
}
repaired_schema2 <- hpc_env$load_pb_variants(
  seurat = NULL,
  sample_col = "Sample",
  hvg_rank_genes = paste0("g", 1:5),
  pseudobulk_dir = schema2_dir,
  ds = "Toy",
  h5ad_path = file.path(schema2_dir, "synthetic.h5ad")
)
hpc_env$prepare_pseudobulks_hpc <- original_prepare_pseudobulks
stopifnot(
  identical(prepared_missing, "hvg500"),
  identical(names(repaired_schema2), hpc_env$PB_VARIANT_NAMES)
)
bad_shared_total_path <- file.path(schema2_dir, "bad-shared-total.rds")
bad_shared_total <- schema2_record("hvg2000")
bad_shared_total$shared_time_secs <- 4.5
write_checked_fixture(bad_shared_total_path, bad_shared_total)
expect_validator_failure(bad_shared_total_path, "inconsistent shared total")
bad_timing_id_path <- file.path(schema2_dir, "bad-timing-id.rds")
bad_timing_id <- schema2_record("hvg2000")
bad_timing_id$timing_id <- ""
write_checked_fixture(bad_timing_id_path, bad_timing_id)
expect_validator_failure(bad_timing_id_path, "blank timing ID")
bad_timing_id$timing_id <- "malformed"
write_checked_fixture(bad_timing_id_path, bad_timing_id)
expect_validator_failure(bad_timing_id_path, "malformed timing ID")
bad_schema_fields_path <- file.path(schema2_dir, "bad-schema-fields.rds")
bad_schema_fields <- schema2_record("hvg2000")
bad_schema_fields$shared_mem_GB <- NULL
write_checked_fixture(bad_schema_fields_path, bad_schema_fields)
expect_validator_failure(bad_schema_fields_path, "missing schema-2 field")
bad_schema_version_path <- file.path(schema2_dir, "bad-schema-version.rds")
bad_schema_version <- schema2_record("hvg2000")
bad_schema_version$timing_schema <- 3L
write_checked_fixture(bad_schema_version_path, bad_schema_version)
expect_validator_failure(bad_schema_version_path, "unsupported schema version")
bad_variant_time_path <- file.path(schema2_dir, "bad-variant-time.rds")
bad_variant_time <- schema2_record("hvg2000")
bad_variant_time$time_secs <- 0.5
write_checked_fixture(bad_variant_time_path, bad_variant_time)
expect_validator_failure(bad_variant_time_path, "inconsistent variant timing")
bad_memory_path <- file.path(schema2_dir, "bad-memory.rds")
bad_memory <- schema2_record("hvg2000")
bad_memory$shared_mem_GB <- Inf
write_checked_fixture(bad_memory_path, bad_memory)
expect_validator_failure(bad_memory_path, "nonfinite shared memory")
timing_rows <- list()
original_hpc_log_exec_row <- hpc_env$log_exec_row
hpc_env$log_exec_row <- function(
  dataset, method, time_secs, log_file, mem_gb = NA_real_, ...
) {
  timing_rows[[length(timing_rows) + 1L]] <<- list(
    dataset = dataset,
    method = method,
    time_secs = as.numeric(time_secs),
    mem_GB = mem_gb
  )
  invisible(NULL)
}
schema2_variants <- setNames(
  lapply(hpc_env$PB_VARIANT_NAMES, schema2_record),
  hpc_env$PB_VARIANT_NAMES
)
hpc_env$emit_pseudobulk_timing_rows(
  schema2_variants,
  ds = "Toy",
  log_file = "fixture-execution-times.feather"
)
shared_rows <- vapply(
  timing_rows,
  function(row) identical(row$method, "prepare_pseudobulk_shared"),
  logical(1L)
)
stopifnot(
  sum(shared_rows) == 1L,
  length(timing_rows) == length(hpc_env$PB_VARIANT_NAMES) + 1L,
  identical(timing_rows[[which(shared_rows)[[1L]]]]$time_secs, 4),
  sum(vapply(timing_rows, function(row) row$time_secs, numeric(1L))) ==
    4 + 0.25 * length(hpc_env$PB_VARIANT_NAMES)
)
different_timing_id <- schema2_variants
different_timing_id[["hvg500"]]$timing_id <- "other-run:Toy:benchmark_analysis:none"
assert_error(
  hpc_env$emit_pseudobulk_timing_rows(
    different_timing_id,
    ds = "Toy",
    log_file = "fixture-execution-times.feather"
  ),
  "timing"
)
blank_timing_id <- schema2_variants
blank_timing_id[["hvg500"]]$timing_id <- ""
assert_error(
  hpc_env$emit_pseudobulk_timing_rows(
    blank_timing_id,
    ds = "Toy",
    log_file = "fixture-execution-times.feather"
  ),
  "timing"
)
hpc_env$log_exec_row <- original_hpc_log_exec_row
pipeline_env$read_rds_checked <- hpc_env$read_rds_checked
pipeline_env$artifact_checksum_ok <- function(file) file.exists(file) && file.info(file)$size > 0

pseudobulk_dir <- tempfile("ecoda_pseudobulk_cache-")
dir.create(pseudobulk_dir, recursive = TRUE)
for (variant in hpc_env$PB_VARIANT_NAMES) {
  path <- file.path(pseudobulk_dir, paste0("Toy_pseudobulk_", variant, ".rds"))
  saveRDS(
    list(pb = matrix(1, nrow = 1, ncol = 1), time_secs = 0),
    path
  )
  writeLines(
    c(
      paste0("MD5=", unname(tools::md5sum(path))),
      paste0("SIZE=", file.info(path)$size),
      paste0("PATH=", path)
    ),
    paste0(path, ".md5")
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
stopifnot(identical(
  hpc_env$PB_VARIANT_PRODUCERS,
  setNames(
    paste0("stage5_prepare_pseudobulk_", hpc_env$PB_VARIANT_NAMES),
    hpc_env$PB_VARIANT_NAMES
  )
))

record_env <- Sys.getenv(
  c("ECODA_RUN_ID", "ECODA_RUNS_ROOT", "ECODA_ARTIFACT_PRODUCER"),
  unset = NA_character_
)
recorded_pseudobulk_dir <- tempfile("ecoda_pseudobulk_records-")
recorded_runs_root <- tempfile("ecoda_pseudobulk_runs-")
dir.create(recorded_pseudobulk_dir, recursive = TRUE)
dir.create(recorded_runs_root, recursive = TRUE)
record_run_id <- "benchmark-cache-records"
Sys.setenv(
  ECODA_RUN_ID = record_run_id,
  ECODA_RUNS_ROOT = normalizePath(recorded_runs_root, mustWork = FALSE),
  ECODA_ARTIFACT_PRODUCER = "mofa"
)
recorded_object <- function(variant) {
  list(
    pb = matrix(match(variant, hpc_env$PB_VARIANT_NAMES), nrow = 1L),
    time_secs = 0
  )
}
for (variant in hpc_env$PB_VARIANT_NAMES) {
  hpc_env$save_rds_atomic(
    recorded_object(variant),
    file.path(
      recorded_pseudobulk_dir,
      paste0("Toy_pseudobulk_", variant, ".rds")
    ),
    producer = hpc_env$PB_VARIANT_PRODUCERS[[variant]],
    run_id = record_run_id
  )
}
stopifnot(identical(
  hpc_env$pb_variants_missing(recorded_pseudobulk_dir, "Toy"),
  character(0)
))
recorded_loaded <- hpc_env$load_pb_variants(
  seurat = NULL,
  sample_col = "Sample",
  hvg_rank_genes = character(0),
  pseudobulk_dir = recorded_pseudobulk_dir,
  ds = "Toy"
)
stopifnot(identical(
  names(recorded_loaded),
  hpc_env$PB_VARIANT_NAMES
))

tamper_path <- file.path(
  recorded_pseudobulk_dir,
  paste0("Toy_pseudobulk_", hpc_env$PB_VARIANT_NAMES[[1L]], ".rds")
)
tampered_bytes <- readBin(
  tamper_path, what = "raw", n = file.info(tamper_path)$size
)
tampered_bytes[[1L]] <- as.raw(bitwXor(as.integer(tampered_bytes[[1L]]), 1L))
writeBin(tampered_bytes, tamper_path)
assert_error(
  hpc_env$pb_variants_missing(recorded_pseudobulk_dir, "Toy"),
  "Artifact checksum validation failed"
)
assert_error(
  hpc_env$load_pb_variants(
    seurat = NULL,
    sample_col = "Sample",
    hvg_rank_genes = character(0),
    pseudobulk_dir = recorded_pseudobulk_dir,
    ds = "Toy",
    variants = hpc_env$PB_VARIANT_NAMES[[1L]]
  ),
  "Artifact checksum validation failed"
)

mismatched_pseudobulk_dir <- tempfile("ecoda_pseudobulk_mismatch-")
dir.create(mismatched_pseudobulk_dir, recursive = TRUE)
for (variant in hpc_env$PB_VARIANT_NAMES) {
  hpc_env$save_rds_atomic(
    recorded_object(variant),
    file.path(
      mismatched_pseudobulk_dir,
      paste0("Toy_pseudobulk_", variant, ".rds")
    ),
    producer = hpc_env$PB_VARIANT_PRODUCERS[[variant]],
    run_id = record_run_id
  )
}
mismatched_variant <- hpc_env$PB_VARIANT_NAMES[[1L]]
hpc_env$save_rds_atomic(
  recorded_object(mismatched_variant),
  file.path(
    mismatched_pseudobulk_dir,
    paste0("Toy_pseudobulk_", mismatched_variant, ".rds")
  ),
  producer = "mofa",
  run_id = record_run_id
)
assert_error(
  hpc_env$pb_variants_missing(mismatched_pseudobulk_dir, "Toy"),
  "Artifact record binding is invalid"
)
assert_error(
  hpc_env$load_pb_variants(
    seurat = NULL,
    sample_col = "Sample",
    hvg_rank_genes = character(0),
    pseudobulk_dir = mismatched_pseudobulk_dir,
    ds = "Toy",
    variants = mismatched_variant
  ),
  "Artifact record binding is invalid"
)
unlink(
  c(recorded_pseudobulk_dir, recorded_runs_root, mismatched_pseudobulk_dir),
  recursive = TRUE,
  force = TRUE
)
for (name in names(record_env)) {
  if (is.na(record_env[[name]])) {
    Sys.unsetenv(name)
  } else {
    do.call(Sys.setenv, setNames(list(record_env[[name]]), name))
  }
}

unlink(pseudobulk_dir, recursive = TRUE, force = TRUE)
unlink(schema2_dir, recursive = TRUE, force = TRUE)

# `--force` still invalidates composition result bundles, while the
# obs-only pseudobulk loader reuses the prepared cache above.
pipeline_env$peak_rss_gb <- function() NA_real_
pipeline_env$save_rds_atomic <- function(object, file) {
  hpc_env$save_rds_atomic(object, file)
}
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
  Sample = factor(
    c("sample_2", "sample_1", "sample_2"),
    levels = c("sample_1", "sample_2")
  ),
  stringsAsFactors = FALSE
)
composition_labels <- factor(c("group_1", "group_2"))
names(composition_labels) <- unique(as.character(composition_obs$Sample))
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
metadata_bundle <- readRDS(file.path(
  composition_results_dir,
  "Toy_metadata.rds"
))
stopifnot(identical(
  names(metadata_bundle$labels),
  c("sample_2", "sample_1")
))
stopifnot(identical(
  names(metadata_bundle$cells_per_sample),
  c("sample_2", "sample_1")
))
stopifnot(identical(
  as.integer(metadata_bundle$cells_per_sample),
  c(2L, 1L)
))
composition_bundle <- file.path(
  composition_results_dir,
  "Toy_Avg_PCA_embedding.rds"
)
hpc_env$save_rds_atomic(
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
