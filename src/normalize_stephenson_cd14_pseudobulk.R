#!/usr/bin/env Rscript
# Run-owned CD14 LR normalization; biological labels never enter the fit.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  required <- c("--input-dir", "--output-dir", "--normalization-source")
  if (length(args) != 6L || anyDuplicated(args[c(1L, 3L, 5L)]) ||
      !setequal(args[c(1L, 3L, 5L)], required)) {
    stop("Required: --input-dir ABS --output-dir NEW_ABS --normalization-source ABS")
  }
  options <- stats::setNames(as.list(args[c(2L, 4L, 6L)]), args[c(1L, 3L, 5L)])
  if (any(!startsWith(unlist(options), "/"))) stop("All paths must be absolute")
  input <- options[["--input-dir"]]
  output <- options[["--output-dir"]]
  normalization_source <- options[["--normalization-source"]]
  output_link <- Sys.readlink(output)
  if (file.exists(output) || dir.exists(output) ||
      (!is.na(output_link) && nzchar(output_link))) {
    stop("Output directory already exists or is a symlink")
  }
  hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
  manifest_path <- file.path(input, "aggregation.json")
  aggregation_hash <- hash(manifest_path)
  aggregation <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  expected <- list(schema_version = 1L, dataset = "Stephenson", view = "benchmark_analysis",
                   annotation_column = "initial_clustering", cell_type = "CD14")
  for (key in names(expected)) {
    if (!identical(aggregation[[key]], expected[[key]])) stop("Invalid aggregation field: ", key)
  }
  ids <- unlist(aggregation$sample_ids, use.names = FALSE)
  if (!is.character(ids) || !length(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) ||
      anyDuplicated(ids) || !identical(ids, sort(ids, method = "radix"))) {
    stop("Invalid aggregation sample axis")
  }
  for (name in c("counts.feather", "samples.feather")) {
    declared <- aggregation$files[[name]]
    if (!is.character(declared) || length(declared) != 1L ||
        !grepl("^[0-9a-f]{64}$", declared) ||
        !identical(hash(file.path(input, name)), declared)) {
      stop("Raw artifact SHA256 mismatch: ", name)
    }
  }
  if (is.null(aggregation$source) || !is.list(aggregation$source) ||
      !is.character(aggregation$source$path) ||
      !identical(aggregation$source$full_file_checksum_verified, FALSE)) {
    stop("Missing or invalid source provenance")
  }
  samples <- as.data.frame(arrow::read_feather(file.path(input, "samples.feather")))
  if (!identical(names(samples), c("Sample", "Status", "n_cells")) ||
      !identical(samples$Sample, ids) || anyNA(samples) ||
      !all(samples$Status %in% c("Covid", "Healthy")) ||
      !is.numeric(samples$n_cells) || any(!is.finite(samples$n_cells)) ||
      any(samples$n_cells < 5 | samples$n_cells != floor(samples$n_cells)) ||
      !identical(as.numeric(sum(samples$n_cells)), as.numeric(aggregation$n_selected_cells))) {
    stop("Invalid sample metadata or selected-cell count")
  }
  raw <- as.data.frame(arrow::read_feather(file.path(input, "counts.feather")))
  if (!identical(names(raw), c("gene", ids)) || !is.character(raw$gene) ||
      anyNA(raw$gene) || any(!nzchar(trimws(raw$gene))) || anyDuplicated(raw$gene) ||
      nrow(raw) < 2000L || !identical(as.numeric(nrow(raw)), as.numeric(aggregation$n_genes))) {
    stop("Invalid raw gene/sample axes (at least 2000 genes required)")
  }
  # Explicit per-column conversion also handles Arrow integer64 without bit reinterpretation.
  count_columns <- lapply(raw[ids], function(column) {
    if (!is.numeric(column) && !inherits(column, "integer64")) stop("Nonnumeric counts")
    as.numeric(column)
  })
  counts <- as.matrix(as.data.frame(count_columns, check.names = FALSE))
  dimnames(counts) <- list(raw$gene, ids)
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts)) || any(counts > .Machine$integer.max)) {
    stop("Raw counts must be nonnegative integers within R integer capacity")
  }
  source_hash <- hash(normalization_source)
  normalization_env <- new.env(parent = globalenv())
  sys.source(normalization_source, envir = normalization_env)
  expression <- normalization_env$get_pb_deseq2_from_counts(
    counts = counts, metadata = data.frame(Sample = ids), hvg = NULL,
    n_hvg = 2000, black_list = "none", batch_col = NULL, blind = TRUE,
    correct_batch = FALSE
  )
  if (!is.matrix(expression) || !identical(rownames(expression), ids) ||
      ncol(expression) != 2000L || any(!is.finite(expression)) ||
      is.null(colnames(expression)) || anyDuplicated(colnames(expression)) ||
      any(!nzchar(colnames(expression))) || "Sample" %in% colnames(expression)) {
    stop("Invalid normalized expression matrix")
  }
  if (!identical(hash(manifest_path), aggregation_hash) ||
      !identical(hash(normalization_source), source_hash)) stop("Inputs changed during normalization")
  for (name in c("counts.feather", "samples.feather")) {
    if (!identical(hash(file.path(input, name)), aggregation$files[[name]])) {
      stop("Raw artifact changed during normalization: ", name)
    }
  }
  if (!dir.create(output, recursive = FALSE, showWarnings = FALSE)) stop("Cannot create new output")
  exported <- data.frame(Sample = ids, expression, check.names = FALSE)
  arrow::write_feather(exported, file.path(output, "expression.feather"))
  if (!file.copy(file.path(input, "samples.feather"), file.path(output, "samples.feather"),
                 overwrite = FALSE)) stop("Cannot copy sample metadata")
  result <- list(
    schema_version = 1L, dataset = "Stephenson", view = "benchmark_analysis",
    annotation_column = "initial_clustering", cell_type = "CD14", sample_ids = ids,
    n_selected_cells = aggregation$n_selected_cells, n_features = ncol(expression),
    normalization = list(blind = TRUE, batch_col = NULL, correct_batch = FALSE,
                         n_hvg = 2000L, black_list = "none"),
    source = aggregation$source,
    aggregation = list(path = manifest_path, sha256 = aggregation_hash, manifest = aggregation),
    normalization_code_sha256 = source_hash,
    runtime = list(R = R.version.string, DESeq2 = as.character(utils::packageVersion("DESeq2"))),
    files = stats::setNames(lapply(c("expression.feather", "samples.feather"), function(name) {
      hash(file.path(output, name))
    }), c("expression.feather", "samples.feather"))
  )
  jsonlite::write_json(result, file.path(output, "manifest.json"), auto_unbox = TRUE,
                       pretty = TRUE, null = "null", digits = NA)
}

main()
