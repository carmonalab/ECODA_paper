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
