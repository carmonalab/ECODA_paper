# ============================================================
# HPC BENCHMARK WORKER UTILITIES
# Shared by the R benchmark workers (Pipeline A:
# run_r_sample_embedding_methods/; Pipeline B:
# run_transformation_zeroimp_analysis/). NOT added to
# load_all_functions.R (HPC-only); scripts source this file
# explicitly after src/utils/load_worker_functions.R.
# ============================================================

# Single source of truth for the shared DESeq2 pseudobulk variants
# (prepare_pseudobulks_hpc spec names; used by the prep worker,
# load_pb_variants and the per-variant cache files).
PB_VARIANT_NAMES <- c(
  "schvg2000", "hvg2000", "hvg500", "hvg2000_bl", "hvg1000", "hvg3000"
)
# Every prepared cache record is published by the Stage 5 producer for its
# variant, regardless of which downstream benchmark method consumes it.
PB_VARIANT_PRODUCERS <- setNames(
  paste0("stage5_prepare_pseudobulk_", PB_VARIANT_NAMES),
  PB_VARIANT_NAMES
)

# Validate a prepared pseudobulk cache against its canonical producer. A
# run-owned record is authoritative: malformed or mismatched records fail
# closed instead of being treated as a cache miss. Artifacts from before
# run-owned records were introduced retain the strict sidecar fallback.
.pb_variant_cache_valid <- function(path, variant) {
  producer <- PB_VARIANT_PRODUCERS[[variant]]
  context <- .artifact_context(producer = producer)
  if (!is.null(context)) {
    record_path <- artifact_record_path(
      path, context$run_id, runs_root = context$runs_root
    )
    if (file.exists(record_path)) {
      record <- artifact_record_for_load(
        path, producer = producer, run_id = context$run_id
      )
      if (is.null(record)) {
        stop("Artifact record validation failed: ", record_path)
      }
      sidecar <- .artifact_sidecar(path, verify = TRUE)
      if (is.null(sidecar) ||
          !identical(sidecar$MD5, record$MD5) ||
          !identical(sidecar$SIZE, record$SIZE)) {
        stop("Artifact checksum validation failed: ", path)
      }
      return(TRUE)
    }
  }
  artifact_checksum_ok(path, producer = producer)
}

# Tiny "--flag value" / "--flag=value" / "--flag" (TRUE) arg parser
parse_flags <- function(raw_args) {
  args <- list()
  i <- 1
  while (i <= length(raw_args)) {
    flag <- raw_args[i]
    if (!startsWith(flag, "--")) {
      stop("Unexpected positional argument: ", flag)
    }
    name <- sub("^--", "", flag)
    if (grepl("=", name, fixed = TRUE)) {
      kv <- strsplit(name, "=", fixed = TRUE)[[1]]
      args[[kv[1]]] <- kv[2]
      i <- i + 1
    } else {
      if (i < length(raw_args) && !startsWith(raw_args[i + 1], "--")) {
        args[[name]] <- raw_args[i + 1]
        i <- i + 2
      } else {
        args[[name]] <- TRUE
        i <- i + 1
      }
    }
  }
  return(args)
}

# Resolve the preprocessed view h5ad path for a dataset. `config` comes from
# read_datasets_json(config_path, view = view) (datasets_io.R), which already
# maps columns.label -> label_col and output_file_name -> output_file.
get_h5ad_path <- function(config, ds, view, input_dir) {
  entry <- config[[ds]]
  if (is.null(entry)) {
    stop("Dataset '", ds, "' not found in datasets.json (view '", view, "').")
  }
  views <- entry[["views"]]
  if (is.null(views) || is.null(views[[view]])) {
    stop("Dataset '", ds, "' has no '", view, "' view in datasets.json.")
  }
  out_file <- views[[view]][["output_file"]]
  if (is.null(out_file) || is.na(out_file) || out_file == "") {
    stop("Dataset '", ds, "' view '", view,
         "' has no output_file_name in datasets.json.")
  }
  return(file.path(input_dir, out_file))
}

# Full ranked gene list from the stored var["hvg_rank"] (set by
# 1.1.1_preprocess.py's select_hvgs_ranked; batch-aware). Mirrors the Python
# top_n_hvg_genes() (dropna + ascending rank sort, names = genes).
get_hvg_rank_genes <- function(adata) {
  var_df <- py_to_r(adata$var)
  if (!"hvg_rank" %in% colnames(var_df)) {
    stop("Column 'hvg_rank' not found in adata.var. Re-run preprocessing ",
         "(1.1.1_preprocess.py).")
  }
  ranks <- var_df[["hvg_rank"]]
  keep <- !is.na(ranks)
  genes <- rownames(var_df)[keep]
  genes[order(ranks[keep])]
}

# Named hvg sets (hvg1000/hvg2000/hvg3000) for the scITD driver; sizes
# clamped to the available ranked genes.
make_hvg_sets <- function(hvg_rank_genes, sizes = c(1000, 2000, 3000)) {
  sets <- lapply(sizes, function(n) {
    hvg_rank_genes[seq_len(min(n, length(hvg_rank_genes)))]
  })
  names(sets) <- paste0("hvg", sizes)
  return(sets)
}

# Variants of the shared DESeq2 pseudobulks missing from pseudobulk_dir
# (all of them under --force).
pb_variants_missing <- function(pseudobulk_dir, ds, force = FALSE,
                                cache_stem = ds) {
  if (force) return(PB_VARIANT_NAMES)
  valid <- vapply(PB_VARIANT_NAMES, function(variant) {
    path <- file.path(
      pseudobulk_dir,
      paste0(cache_stem, "_pseudobulk_", variant, ".rds")
    )
    .pb_variant_cache_valid(path, variant)
  }, logical(1))
  PB_VARIANT_NAMES[!valid]
}

# Validate the on-disk AnnData contract before any benchmark worker builds a
# Seurat object or consumes a precomputed pseudobulk. This deliberately checks
# content, not only file existence: a metadata/PCA-only h5ad can otherwise
# pass the scheduler gate and fail later as a misleading zero result.
validate_benchmark_h5ad_contract <- function(
  adata,
  obs = NULL,
  view = "benchmark_analysis",
  method = NULL
) {
  required_obsm <- switch(
    view,
    benchmark_analysis = c(
      "X_pca_benchmark_analysis_hvg1000",
      "X_pca_benchmark_analysis_hvg2000",
      "X_pca_benchmark_analysis_hvg3000",
      "X_pca_harmony_benchmark_analysis_hvg2000"
    ),
    batch_effect_uncorrected = "X_pca_batch_effect_uncorrected_hvg2000",
    batch_effect_corrected = c(
      "X_pca_batch_effect_corrected_hvg2000",
      "X_pca_harmony_batch_effect_corrected_hvg2000"
    ),
    stop("Unknown preprocessing view for h5ad contract: ", view)
  )
  required_hvg <- if (identical(view, "benchmark_analysis")) 3000 else 2000

  layer_keys <- py_to_r(import_builtins(convert = FALSE)$list(
    adata$layers$keys()
  ))
  missing <- character()
  if (!"counts" %in% layer_keys) {
    missing <- c(missing, "layers['counts']")
  }

  x_present <- tryCatch(
    !is.null(adata$X),
    error = function(e) FALSE
  )
  if (!x_present) {
    missing <- c(missing, "X")
  }
  dimensions <- tryCatch(
    as.integer(py_to_r(adata$shape)),
    error = function(e) integer()
  )
  if (length(dimensions) < 2L || any(dimensions[seq_len(2L)] <= 0L)) {
    missing <- c(missing, "non-empty obs/var")
  }

  obsm_keys <- py_to_r(import_builtins(convert = FALSE)$list(
    adata$obsm$keys()
  ))
  missing_obsm <- setdiff(required_obsm, obsm_keys)
  if (length(missing_obsm) > 0) {
    missing <- c(missing, paste0("obsm['", missing_obsm, "']"))
  }

  var_df <- py_to_r(adata$var)
  if (!"hvg_rank" %in% colnames(var_df)) {
    missing <- c(missing, "var['hvg_rank']")
  } else if (sum(!is.na(var_df[["hvg_rank"]])) < required_hvg) {
    missing <- c(
      missing,
      paste0(
        "var['hvg_rank'] with at least ", required_hvg,
        " non-NA ranks"
      )
    )
  }

  if (is.null(obs)) {
    obs <- py_to_r(adata$obs)
  }
  if (!"Sample" %in% colnames(obs)) {
    missing <- c(missing, "obs['Sample']")
  }

  if (length(missing) > 0) {
    method_suffix <- if (is.null(method)) "" else paste0(" for method ", method)
    stop(
      "h5ad content contract failed", method_suffix, " (view ", view,
      "): missing or invalid ", paste(missing, collapse = ", "),
      ". Re-run src/3_scrnaseq_preprocessing/1.1.1_preprocess.py with ",
      "--force and use the authoritative processed h5ad."
    )
  }
  invisible(TRUE)
}

# Load only the metadata and requested PCA embeddings from an H5AD.  The
# pinned anndata backed reader can materialize layers["counts"] at open, so
# methods whose algorithms consume no counts use this h5py/minimal-AnnData
# path instead.  Count-dependent methods keep load_benchmark_seurat()'s
# explicit counts-layer path.
load_h5ad_counts_free <- function(
  h5ad_path,
  obs_columns,
  embedding_keys,
  obs_prefixes = character(),
  view = NULL,
  method = NULL
) {
  project_root <- Sys.getenv("PROJECT_ROOT")
  if (project_root == "") {
    stop("PROJECT_ROOT not set; cannot load a counts-free H5AD.")
  }
  module_dir <- normalizePath(
    file.path(project_root, "src", "utils", "py"),
    mustWork = TRUE
  )
  python_sys <- reticulate::import("sys", convert = FALSE)
  python_sys$path$insert(0L, module_dir)
  loader <- reticulate::import_from_path(
    "h5ad_counts_free",
    path = module_dir,
    convert = FALSE
  )
  if (!is.null(view) && !is.null(method)) {
    loader$validate_h5ad_counts_free_input(
      h5ad_path,
      as.character(view),
      as.character(method)
    )
  }
  loader$load_h5ad_counts_free(
    h5ad_path,
    as.list(as.character(obs_columns)),
    as.list(as.character(embedding_keys)),
    as.list(as.character(obs_prefixes))
  )
}

# Read only the first metadata row for each sample. This keeps MOFA and
# batch-mode pseudobulk result workers counts-free when their pseudobulks are
# already cached.
load_h5ad_sample_metadata <- function(
  h5ad_path,
  sample_col = "Sample",
  metadata_columns = character(),
  chunk_size = 4096L
) {
  project_root <- Sys.getenv("PROJECT_ROOT")
  if (project_root == "") {
    stop("PROJECT_ROOT not set; cannot load sample metadata.")
  }
  module_dir <- normalizePath(
    file.path(project_root, "src", "utils", "py"),
    mustWork = TRUE
  )
  python_sys <- reticulate::import("sys", convert = FALSE)
  python_sys$path$insert(0L, module_dir)
  loader <- reticulate::import_from_path(
    "h5ad_pseudobulk",
    path = module_dir,
    convert = FALSE
  )
  metadata <- reticulate::py_to_r(loader$read_h5ad_sample_metadata(
    h5ad_path,
    sample_col,
    as.list(as.character(unique(c(sample_col, metadata_columns)))),
    as.integer(chunk_size)
  ))
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  }
  if (!sample_col %in% colnames(metadata) ||
      is.null(rownames(metadata)) ||
      anyNA(rownames(metadata)) ||
      any(!nzchar(rownames(metadata)))) {
    stop("Streaming sample metadata has no valid sample index.")
  }
  metadata
}
# Read the raw Sample aggregate directly from the H5AD CSR layer.  This is
# the canonical Stage 5 path: the returned genes-by-samples matrix crosses
# the reticulate boundary once and is passed directly to the matrix DESeq2
# helpers; no one-sample-per-column Seurat object is created.
load_h5ad_sample_aggregate <- function(
  h5ad_path,
  sample_col = "Sample",
  metadata_columns = character(),
  chunk_size = 4096L,
  max_value = .Machine$integer.max
) {
  project_root <- Sys.getenv("PROJECT_ROOT")
  if (project_root == "") {
    stop("PROJECT_ROOT not set; cannot aggregate an H5AD.")
  }
  module_dir <- normalizePath(
    file.path(project_root, "src", "utils", "py"),
    mustWork = TRUE
  )
  python_sys <- reticulate::import("sys", convert = FALSE)
  python_sys$path$insert(0L, module_dir)
  loader <- reticulate::import_from_path(
    "h5ad_pseudobulk",
    path = module_dir,
    convert = FALSE
  )
  columns <- unique(c(sample_col, metadata_columns))
  if (is.null(max_value)) {
    aggregated <- loader$aggregate_h5ad_counts_by_sample(
      h5ad_path,
      sample_col,
      as.list(as.character(columns)),
      as.integer(chunk_size)
    )
  } else {
    aggregated <- loader$aggregate_h5ad_counts_by_sample(
      h5ad_path,
      sample_col,
      as.list(as.character(columns)),
      as.integer(chunk_size),
      as.numeric(max_value)
    )
  }
  counts <- reticulate::py_to_r(aggregated$counts)
  sample_ids <- as.character(reticulate::py_to_r(aggregated$sample_ids))
  gene_names <- as.character(reticulate::py_to_r(aggregated$gene_names))
  metadata <- reticulate::py_to_r(aggregated$metadata)
  if (!is.matrix(counts)) counts <- as.matrix(counts)
  if (length(dim(counts)) != 2L ||
      nrow(counts) != length(gene_names) ||
      ncol(counts) != length(sample_ids)) {
    stop("H5AD Sample aggregate has inconsistent matrix dimensions.")
  }
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  }
  if (nrow(metadata) != length(sample_ids) ||
      !sample_col %in% colnames(metadata)) {
    stop("H5AD Sample aggregate metadata has inconsistent sample rows.")
  }
  metadata_samples <- as.character(metadata[[sample_col]])
  if (length(metadata_samples) != length(sample_ids) ||
      anyNA(metadata_samples) ||
      !identical(metadata_samples, sample_ids)) {
    stop("H5AD Sample aggregate metadata is not aligned to first-seen samples.")
  }
  if (anyDuplicated(sample_ids) || anyNA(sample_ids) ||
      any(!nzchar(trimws(sample_ids)))) {
    stop("H5AD Sample aggregate IDs are invalid.")
  }
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) {
    stop("H5AD Sample aggregate counts are not finite nonnegative integers.")
  }
  if (!is.null(max_value) && any(counts > max_value)) {
    stop("H5AD Sample aggregate exceeds the requested count ceiling.")
  }
  rownames(metadata) <- sample_ids
  rownames(counts) <- gene_names
  colnames(counts) <- sample_ids
  list(
    counts = counts,
    sample_ids = sample_ids,
    gene_names = gene_names,
    metadata = metadata
  )
}

# Metadata/HVG-only loading used before a canonical pseudobulk cache decides
# whether the raw counts layer is needed.  h5ad_obs_free validates persisted
# shapes and reads only obs; h5ad_counts_subset reads var["hvg_rank"] without
# reading any X/layers values.
load_h5ad_pseudobulk_metadata <- function(
  h5ad_path,
  sample_col = "Sample",
  metadata_columns = character(),
  n_hvg = 3000L,
  required_nonmissing_columns = sample_col
) {
  project_root <- Sys.getenv("PROJECT_ROOT")
  if (project_root == "") {
    stop("PROJECT_ROOT not set; cannot load H5AD metadata.")
  }
  module_dir <- normalizePath(
    file.path(project_root, "src", "utils", "py"),
    mustWork = TRUE
  )
  python_sys <- reticulate::import("sys", convert = FALSE)
  python_sys$path$insert(0L, module_dir)
  obs_loader <- reticulate::import_from_path(
    "h5ad_obs_free", path = module_dir, convert = FALSE
  )
  hvg_loader <- reticulate::import_from_path(
    "h5ad_counts_subset", path = module_dir, convert = FALSE
  )
  columns <- unique(c(sample_col, metadata_columns))
  obs <- reticulate::py_to_r(obs_loader$load_h5ad_obs_free(
    h5ad_path, as.list(as.character(columns))
  ))
  hvg_rank_genes <- as.character(reticulate::py_to_r(
    hvg_loader$read_h5ad_hvg_genes(h5ad_path, as.integer(n_hvg))
  ))
  if (!is.data.frame(obs) || !sample_col %in% colnames(obs)) {
    stop("H5AD metadata has no valid ", sample_col, " column.")
  }
  required_nonmissing_columns <- unique(as.character(
    required_nonmissing_columns
  ))
  missing_required <- setdiff(required_nonmissing_columns, colnames(obs))
  if (length(missing_required) > 0L) {
    stop("H5AD metadata lacks required columns: ",
         paste(missing_required, collapse = ", "))
  }
  for (column in required_nonmissing_columns) {
    values <- as.character(obs[[column]])
    if (anyNA(values) || any(!nzchar(trimws(values))) ||
        any(tolower(values) == "nan")) {
      stop("H5AD metadata column '", column,
           "' contains missing or blank values.")
    }
  }
  if (length(hvg_rank_genes) == 0L || anyNA(hvg_rank_genes) ||
      any(!nzchar(hvg_rank_genes)) || anyDuplicated(hvg_rank_genes)) {
    stop("H5AD hvg_rank metadata is invalid.")
  }
  list(obs = obs, hvg_rank_genes = hvg_rank_genes)
}

# Legacy-only adapter: this retains the historical Seurat boundary for
# maintained callers, but canonical Stage 5 preparation/fallback never calls
# it.  A sample-level Seurat round trip must not be used for AggregateExpression.
load_h5ad_pseudobulk_seurat <- function(
  h5ad_path,
  sample_col = "Sample",
  batch_col = NULL,
  chunk_size = 4096L
) {
  project_root <- Sys.getenv("PROJECT_ROOT")
  if (project_root == "") {
    stop("PROJECT_ROOT not set; cannot aggregate a pseudobulk.")
  }
  module_dir <- normalizePath(
    file.path(project_root, "src", "utils", "py"),
    mustWork = TRUE
  )
  python_sys <- reticulate::import("sys", convert = FALSE)
  python_sys$path$insert(0L, module_dir)
  loader <- reticulate::import_from_path(
    "h5ad_pseudobulk",
    path = module_dir,
    convert = FALSE
  )
  metadata_columns <- unique(c(sample_col, batch_col))
  aggregated <- loader$aggregate_h5ad_counts_by_sample(
    h5ad_path,
    sample_col,
    as.list(as.character(metadata_columns)),
    as.integer(chunk_size)
  )
  counts <- reticulate::py_to_r(aggregated$counts)
  sample_ids <- as.character(reticulate::py_to_r(aggregated$sample_ids))
  gene_names <- as.character(reticulate::py_to_r(aggregated$gene_names))
  metadata <- reticulate::py_to_r(aggregated$metadata)
  if (!is.matrix(counts)) counts <- as.matrix(counts)
  if (length(dim(counts)) != 2L ||
      nrow(counts) != length(gene_names) ||
      ncol(counts) != length(sample_ids)) {
    stop("Streaming pseudobulk counts have inconsistent dimensions.")
  }
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  }
  if (nrow(metadata) != length(sample_ids) ||
      !sample_col %in% colnames(metadata)) {
    stop("Streaming pseudobulk metadata has inconsistent sample rows.")
  }
  if (anyDuplicated(sample_ids) || anyNA(sample_ids) ||
      any(!nzchar(sample_ids))) {
    stop("Streaming pseudobulk sample IDs are invalid.")
  }
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) {
    stop("Streaming pseudobulk counts are not finite nonnegative integers.")
  }
  rownames(metadata) <- sample_ids
  rownames(counts) <- gene_names
  colnames(counts) <- sample_ids
  Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
}


# Build a Seurat object from a validated H5AD. Count materialization is
# opt-in: GloScope and composition consume only obs plus precomputed embeddings,
# while pseudobulk, scITD, and count-backed model paths request the counts layer.
load_benchmark_seurat <- function(
  adata,
  obs,
  sample_col = "Sample",
  fetch_embedding = c(
    "X_pca_benchmark_analysis_hvg1000",
    "X_pca_benchmark_analysis_hvg2000",
    "X_pca_benchmark_analysis_hvg3000"
  ),
  counts_layer = "counts"
) {
  all_samples <- unique(obs[[sample_col]])
  if (!is.null(counts_layer)) {
    layer_keys <- py_to_r(import_builtins(convert = FALSE)$list(
      adata$layers$keys()
    ))
    if (!counts_layer %in% layer_keys) {
      stop(
        "load_benchmark_seurat: validated h5ad is missing layers['",
        counts_layer, "']; refusing to build a count-backed Seurat object."
      )
    }
  }
  seurat <- get_seurat_obj_from_h5ad(
    adata, obs, all_samples,
    sample_colname = sample_col,
    counts_layer = counts_layer,
    fetch_embedding = fetch_embedding
  )
  return(seurat)
}

# The canonical fallback consumes the H5AD directly.  The first argument is
# retained as an inert compatibility slot for old callers; no canonical path
# may pass a Seurat object here.
pb_timing_id <- function(
  cache_stem,
  view = "benchmark_analysis",
  analysis_pass = NULL,
  run_id = NULL
) {
  if (is.null(run_id) || !nzchar(as.character(run_id))) {
    run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
  }
  if (!nzchar(as.character(run_id))) run_id <- "unbound"
  pass <- if (is.null(analysis_pass)) "none" else as.character(analysis_pass)
  paste(
    as.character(run_id), as.character(cache_stem), as.character(view), pass,
    sep = ":"
  )
}

.pb_timing_schema2 <- function(value, label = "Pseudobulk cache") {
  if (!is.list(value)) stop(label, " is not a list.")
  if (!"timing_schema" %in% names(value)) return(FALSE)
  schema <- value[["timing_schema"]]
  if (!is.numeric(schema) || length(schema) != 1L ||
      is.na(schema) || !is.finite(schema) ||
      schema != 2 || schema != floor(schema)) {
    stop(label, " has an unsupported timing_schema.")
  }
  TRUE
}

.pb_validate_nonnegative <- function(
  value, label, allow_na = FALSE
) {
  if (!is.numeric(value) || length(value) != 1L) {
    stop(label, " is invalid.")
  }
  if (is.na(value)) {
    if (allow_na && !is.nan(value)) return(invisible(TRUE))
    stop(label, " is invalid.")
  }
  if (!is.finite(value) || value < 0) stop(label, " is invalid.")
  invisible(TRUE)
}

pb_variant_time_secs <- function(variant) {
  if (!is.list(variant)) stop("Pseudobulk variant record must be a list.")
  schema2 <- .pb_timing_schema2(variant)
  field <- if (isTRUE(schema2)) "variant_time_secs" else "time_secs"
  if (!field %in% names(variant)) {
    stop("Pseudobulk variant timing is missing ", field, ".")
  }
  .pb_validate_nonnegative(
    variant[[field]], paste0("Pseudobulk ", field)
  )
  as.numeric(variant[[field]])
}

pb_shared_time_secs <- function(variant) {
  if (!is.list(variant)) {
    stop("Pseudobulk variant record must be a list.")
  }
  if (!isTRUE(.pb_timing_schema2(variant))) return(0)
  .pb_validate_nonnegative(
    variant[["shared_time_secs"]],
    "Pseudobulk shared_time_secs"
  )
  as.numeric(variant[["shared_time_secs"]])
}

# Validate cache timing fields independently of whether an execution log is
# being emitted.  This keeps malformed schema-2 payloads from being treated as
# cache hits or from reaching a raw-count fallback.
validate_pseudobulk_timing_record <- function(value, variant = NULL) {
  label <- if (is.null(variant)) "Pseudobulk cache" else {
    paste0("Pseudobulk cache ", variant)
  }
  if (!is.list(value)) stop(label, " is not a list.")
  schema2 <- .pb_timing_schema2(value, label)
  if (!isTRUE(schema2)) {
    if (!"time_secs" %in% names(value)) {
      # A legacy matrix wrapper may have no timing metadata at all.  It remains
      # a legacy cache; callers that need a timing row will reject it at the
      # point where a numeric time is required.
      return(invisible(TRUE))
    }
    .pb_validate_nonnegative(value[["time_secs"]],
                             paste0(label, " legacy time_secs"))
    return(invisible(TRUE))
  }

  required <- c(
    "pb", "time_secs", "mem_GB", "aggregate_time_secs",
    "shared_fit_time_secs", "shared_time_secs", "variant_time_secs",
    "shared_mem_GB", "timing_id", "timing_schema"
  )
  missing <- setdiff(required, names(value))
  if (length(missing) > 0L) {
    stop(label, " is missing timing fields: ", paste(missing, collapse = ", "))
  }

  pb <- value[["pb"]]
  if (!is.matrix(pb) || length(dim(pb)) != 2L ||
      any(dim(pb) <= 0L) || !is.numeric(pb) ||
      any(!is.finite(pb))) {
    stop(label, " has a non-matrix or invalid pb payload.")
  }

  numeric_fields <- c(
    "time_secs", "aggregate_time_secs", "shared_fit_time_secs",
    "shared_time_secs", "variant_time_secs"
  )
  for (field in numeric_fields) {
    .pb_validate_nonnegative(value[[field]], paste0(label, " ", field))
  }
  for (field in c("mem_GB", "shared_mem_GB")) {
    .pb_validate_nonnegative(
      value[[field]], paste0(label, " ", field), allow_na = TRUE
    )
  }

  timing_id <- value[["timing_id"]]
  timing_id_parts <- if (is.character(timing_id) && length(timing_id) == 1L &&
                         !is.na(timing_id)) {
    strsplit(timing_id, ":", fixed = TRUE)[[1L]]
  } else {
    character()
  }
  if (!is.character(timing_id) || length(timing_id) != 1L ||
      is.na(timing_id) || !nzchar(trimws(timing_id)) ||
      length(timing_id_parts) != 4L ||
      any(!nzchar(trimws(timing_id_parts))) ||
      grepl("[[:cntrl:]]", timing_id, perl = TRUE)) {
    stop(label, " has an invalid timing_id.")
  }
  if (!isTRUE(all.equal(
    as.numeric(value[["time_secs"]]),
    as.numeric(value[["variant_time_secs"]]),
    tolerance = 0
  ))) {
    stop(label, " time_secs does not equal variant_time_secs.")
  }
  if (!isTRUE(all.equal(
    as.numeric(value[["shared_time_secs"]]),
    as.numeric(value[["aggregate_time_secs"]]) +
      as.numeric(value[["shared_fit_time_secs"]]),
    tolerance = 0
  ))) {
    stop(label, " shared_time_secs is not aggregate + shared_fit.")
  }
  invisible(TRUE)
}


# Emit one shared timing row and one variant-local row per cache record.  The
# execution Feather schema intentionally remains the historical four columns;
# the cache timing_id is the run-scoped identity used to deduplicate the
# shared stage during report/merge.
emit_pseudobulk_timing_rows <- function(
  variants, ds, log_file = NULL
) {
  if (length(variants) == 0L) return(invisible(NULL))
  if (is.null(names(variants)) || any(!nzchar(names(variants)))) {
    stop("Pseudobulk timing records must be named by variant.")
  }
  has_timing <- vapply(
    variants,
    function(value) is.list(value) && (
      "time_secs" %in% names(value) ||
      "timing_schema" %in% names(value)
    ),
    logical(1)
  )
  for (variant_name in names(variants)[has_timing]) {
    validate_pseudobulk_timing_record(
      variants[[variant_name]], variant = variant_name
    )
  }
  if (is.null(log_file) || is.na(log_file) || !nzchar(log_file)) {
    return(invisible(NULL))
  }
  schema2 <- has_timing & vapply(
    variants,
    function(value) is.list(value) && isTRUE(
      .pb_timing_schema2(value)
    ),
    logical(1)
  )
  if (any(schema2)) {
    shared_ids <- unique(vapply(
      variants[schema2],
      function(value) as.character(value[["timing_id"]]),
      character(1)
    ))
    if (length(shared_ids) != 1L || !nzchar(shared_ids[[1L]])) {
      stop("Pseudobulk schema-2 variants do not share one timing_id.")
    }
    shared_times <- unique(vapply(
      variants[schema2],
      function(value) as.numeric(value[["shared_time_secs"]]),
      numeric(1)
    ))
    if (length(shared_times) != 1L || !is.finite(shared_times[[1L]]) ||
        shared_times[[1L]] < 0) {
      stop("Pseudobulk schema-2 shared timing is invalid.")
    }
    shared_mem <- variants[[which(schema2)[[1L]]]][["shared_mem_GB"]]
    if (!all(vapply(
      variants[schema2],
      function(value) identical(value[["shared_mem_GB"]], shared_mem),
      logical(1)
    ))) {
      stop("Pseudobulk schema-2 shared memory differs.")
    }
    log_exec_row(
      ds, "prepare_pseudobulk_shared", shared_times[[1L]], log_file,
      mem_gb = shared_mem
    )
  }
  for (variant_name in names(variants)[has_timing]) {
    value <- variants[[variant_name]]
    local_time <- pb_variant_time_secs(value)
    if (!is.finite(local_time) || local_time < 0) {
      stop("Pseudobulk variant timing is invalid for ", variant_name)
    }
    log_exec_row(
      ds, paste0("prepare_pseudobulk_", variant_name), local_time, log_file,
      mem_gb = if (is.list(value)) value$mem_GB else NA_real_
    )
  }
  invisible(NULL)
}

# Load the precomputed pseudobulk variants from pseudobulks/.  Cache
# validation and deserialization happen before the H5AD path is touched.
# Missing variants are rebuilt with one bounded raw Sample aggregation and
# one shared full-gene fit; a Seurat object is never used by this fallback.
load_pb_variants <- function(
  seurat = NULL,
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
  variants = PB_VARIANT_NAMES,
  h5ad_path = NULL,
  view = "benchmark_analysis",
  analysis_pass = NULL,
  run_id = NULL,
  source_identity = NULL,
  chunk_size = 4096L
) {
  target_variants <- unique(as.character(variants))
  if (!all(target_variants %in% PB_VARIANT_NAMES)) {
    stop("Unknown pseudobulk variant requested: ",
         paste(setdiff(target_variants, PB_VARIANT_NAMES), collapse = ", "))
  }
  missing <- character()
  variants_out <- list()
  computed_out <- list()
  for (v in target_variants) {
    cache_path <- file.path(
      pseudobulk_dir, paste0(cache_stem, "_pseudobulk_", v, ".rds")
    )
    if (force || !.pb_variant_cache_valid(cache_path, v)) {
      missing <- c(missing, v)
      next
    }
    value <- read_rds_checked(
      cache_path,
      producer = PB_VARIANT_PRODUCERS[[v]]
    )
    validate_pseudobulk_timing_record(value, variant = v)
    variants_out[[v]] <- value
  }
  if (length(missing) > 0L) {
    if (is.null(h5ad_path) || !nzchar(as.character(h5ad_path))) {
      stop("Pseudobulk variant(s) missing in ", pseudobulk_dir,
           " (", paste(missing, collapse = ", "),
           ") and no H5AD path was provided for the raw fallback.")
    }
    message("Pseudobulk variant(s) missing in ", pseudobulk_dir,
            ", computing on the fly: ", paste(missing, collapse = ", "))
    computed <- prepare_pseudobulks_hpc(
      h5ad_path = h5ad_path,
      sample_col = sample_col,
      hvg_rank_genes = hvg_rank_genes,
      variants = missing,
      batch_col = batch_col,
      blind = blind,
      correct_batch = correct_batch,
      cache_stem = cache_stem,
      view = view,
      analysis_pass = analysis_pass,
      run_id = run_id,
      source_identity = source_identity,
      chunk_size = chunk_size
    )
    for (v in names(computed)) {
      value <- computed[[v]]
      validate_pseudobulk_timing_record(value, variant = v)
      cache_path <- file.path(
        pseudobulk_dir, paste0(cache_stem, "_pseudobulk_", v, ".rds")
      )
      save_rds_atomic(
        value, cache_path,
        producer = PB_VARIANT_PRODUCERS[[v]]
      )
      variants_out[[v]] <- value
      computed_out[[v]] <- value
    }
  }
  # A partial repair must not mix an older cached timing identity with the
  # newly computed shared stage in the canonical four-column log.  Cache-only
  # loads retain the historical all-cached replay.
  timing_records <- if (length(missing) > 0L) computed_out else variants_out
  emit_pseudobulk_timing_rows(
    timing_records, ds, log_file = log_file
  )
  variants_out[target_variants]
}

# Composition is obs-only and therefore cannot rebuild missing pseudobulks.
load_composition_pb_variants <- function(
  sample_col,
  hvg_rank_genes,
  pseudobulk_dir,
  ds,
  log_file = NULL,
  loader = load_pb_variants,
  cache_stem = ds,
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE,
  variants = PB_VARIANT_NAMES
) {
  loader(
    seurat = NULL,
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    pseudobulk_dir = pseudobulk_dir,
    ds = ds,
    force = FALSE,
    log_file = log_file,
    cache_stem = cache_stem,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch,
    variants = variants
  )
}


# Artifact records are run-owned, immutable publication metadata. A valid
# record lets cache checks reuse the already computed MD5; an absent record
# falls back to the strict sidecar/content check for legacy artifacts.
.artifact_canonical_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) stop("artifact path must be one non-empty string")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.artifact_sha256_text <- function(value) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(tolower(digest::digest(value, algo = "sha256", serialize = FALSE)))
  }
  commands <- c("shasum", "sha256sum")
  executable <- commands[nzchar(Sys.which(commands))][1L]
  if (is.na(executable)) stop("digest package or SHA-256 utility is required")
  temporary <- tempfile("ecoda_sha256_")
  on.exit(unlink(temporary), add = TRUE)
  writeBin(charToRaw(value), temporary)
  arguments <- if (identical(executable, "shasum")) {
    c("-a", "256", temporary)
  } else {
    temporary
  }
  output <- system2(executable, arguments, stdout = TRUE, stderr = TRUE)
  digest <- sub("[[:space:]].*$", "", output[grepl("^[[:xdigit:]]{64}", output)][1L])
  if (length(digest) != 1L || is.na(digest) ||
      !grepl("^[0-9a-fA-F]{64}$", digest)) {
    stop("could not compute SHA-256 for artifact record path")
  }
  tolower(digest)
}

.artifact_context <- function(producer = NULL, run_id = NULL,
                              execution_log = FALSE) {
  if (is.null(run_id)) run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
  if (is.null(producer)) {
    env_name <- if (execution_log) {
      "ECODA_EXECUTION_LOG_PRODUCER"
    } else {
      "ECODA_ARTIFACT_PRODUCER"
    }
    producer <- Sys.getenv(env_name, unset = "")
    if (!nzchar(producer) && !execution_log) {
      producer <- Sys.getenv("METHOD", unset = "")
      if (!nzchar(producer)) producer <- Sys.getenv("ANALYSIS", unset = "")
    }
    if (execution_log && !nzchar(producer)) producer <- "stage5_execution_log"
  }
  if (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) return(NULL)
  if (!is.character(producer) || length(producer) != 1L || is.na(producer) ||
      !nzchar(producer) || grepl("[\r\n]", producer, perl = TRUE)) return(NULL)
  runs_root <- Sys.getenv("ECODA_RUNS_ROOT", unset = "")
  if (!nzchar(runs_root)) {
    scratch <- Sys.getenv("HPC_SCRATCH_DIR", unset = "")
    if (nzchar(scratch)) runs_root <- file.path(scratch, "_ecoda_runs")
  }
  if (!nzchar(runs_root) || !grepl("^/", runs_root)) return(NULL)
  list(
    run_id = run_id,
    producer = producer,
    runs_root = runs_root
  )
}

artifact_record_path <- function(path, run_id, runs_root = NULL) {
  canonical <- .artifact_canonical_path(path)
  if (!is.character(run_id) || length(run_id) != 1L ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) {
    stop("artifact record run ID is invalid")
  }
  if (is.null(runs_root)) {
    context <- .artifact_context(run_id = run_id, producer = "record")
    if (is.null(context)) stop("artifact record root is unavailable")
    runs_root <- context$runs_root
  } else {
    if (!is.character(runs_root) || length(runs_root) != 1L ||
        is.na(runs_root) || !grepl("^/", runs_root)) {
      stop("artifact record root must be absolute")
    }
    runs_root <- as.character(runs_root)
  }
  key <- substr(.artifact_sha256_text(canonical), 1L, 32L)
  file.path(runs_root, run_id, "manifests", "artifacts",
            paste0(key, ".record"))
}

.artifact_sidecar <- function(path, verify = FALSE, allow_canonical = FALSE) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) return(NULL)
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || !isTRUE(file.info(sidecar)$size > 0)) return(NULL)
  lines <- tryCatch(readLines(sidecar, warn = FALSE),
                    error = function(e) character())
  if (length(lines) != 3L ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    return(NULL)
  }
  md5 <- sub("^MD5=", "", lines[[1L]])
  size <- sub("^SIZE=", "", lines[[2L]])
  recorded <- sub("^PATH=", "", lines[[3L]])
  canonical <- .artifact_canonical_path(path)
  valid_path <- identical(recorded, as.character(path)) ||
    (allow_canonical && identical(recorded, canonical))
  if (!valid_path || !grepl("^[0-9a-fA-F]{32}$", md5) ||
      !grepl("^[1-9][0-9]*$", size) ||
      !identical(size, as.character(file.info(path)$size))) return(NULL)
  if (verify) {
    actual <- unname(tools::md5sum(path))
    if (length(actual) != 1L || is.na(actual) ||
        !identical(tolower(md5), tolower(actual))) return(NULL)
  }
  list(MD5 = tolower(md5), SIZE = size, PATH = recorded)
}

# Atomic RDS write plus a sidecar used by idempotency checks.
save_rds_atomic <- function(object, file, producer = NULL, run_id = NULL) {
  validate_rds_object <- function(value) {
    if (is.null(value) || (is.list(value) && length(value) == 0L)) {
      stop("RDS artifact is empty: ", file)
    }
    dimensions <- dim(value)
    if (length(dimensions) > 0L && any(dimensions <= 0L)) {
      stop("RDS artifact has empty dimensions: ", file)
    }
  }
  validate_rds_object(object)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(file, ".tmp.", Sys.getpid())
  checksum_tmp <- paste0(file, ".md5.tmp.", Sys.getpid())
  sidecar <- paste0(file, ".md5")
  backup <- paste0(file, ".previous.", Sys.getpid())
  sidecar_backup <- paste0(sidecar, ".previous.", Sys.getpid())
  had_file <- file.exists(file)
  had_sidecar <- file.exists(sidecar)
  installed <- FALSE
  sidecar_installed <- FALSE
  digest <- NULL
  size <- NULL
  restore <- function() {
    if (installed && file.exists(file)) unlink(file)
    if (had_file && file.exists(backup)) file.rename(backup, file)
    if (sidecar_installed && file.exists(sidecar)) unlink(sidecar)
    if (had_sidecar && file.exists(sidecar_backup)) file.rename(sidecar_backup, sidecar)
    if (!had_file && file.exists(file)) unlink(file)
    if (!had_sidecar && file.exists(sidecar)) unlink(sidecar)
  }
  on.exit({
    for (temporary in c(tmp, checksum_tmp, backup, sidecar_backup)) {
      if (file.exists(temporary)) unlink(temporary)
    }
  }, add = TRUE)
  tryCatch({
    saveRDS(object, tmp)
    if (!file.exists(tmp) || file.info(tmp)$size <= 0) {
      stop("Empty RDS temporary file: ", tmp)
    }
    if (had_file && !isTRUE(file.link(file, backup))) {
      stop("Could not preserve existing RDS artifact: ", file)
    }
    if (had_sidecar && !isTRUE(file.link(sidecar, sidecar_backup))) {
      stop("Could not preserve existing RDS checksum: ", file)
    }
    if (!file.rename(tmp, file)) stop("Could not atomically install RDS: ", file)
    installed <- TRUE
    digest <- tolower(unname(tools::md5sum(file)))
    size <- as.character(file.info(file)$size)
    writeLines(c(
      paste0("MD5=", digest),
      paste0("SIZE=", size),
      paste0("PATH=", file)
    ), checksum_tmp)
    if (!file.rename(checksum_tmp, sidecar)) {
      stop("Could not atomically install RDS checksum: ", file)
    }
    sidecar_installed <- TRUE
  }, error = function(error) {
    restore()
    stop(error)
  })
  artifact_write_record(file, producer = producer, run_id = run_id,
                        md5 = digest, size = size)
  invisible(NULL)
}

.artifact_record_from_file <- function(path, producer = NULL, run_id = NULL) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (is.null(context)) return(NULL)
  record_path <- artifact_record_path(
    path, context$run_id, runs_root = context$runs_root
  )
  if (!file.exists(record_path) || !isTRUE(file.info(record_path)$size > 0)) {
    return(NULL)
  }
  lines <- tryCatch(readLines(record_path, warn = FALSE),
                    error = function(e) character())
  keys <- c("PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE")
  if (length(lines) != length(keys) ||
      !identical(sub("=.*$", "", lines), keys)) {
    stop("Artifact record has the wrong schema: ", record_path)
  }
  values <- sub("^[^=]*=", "", lines)
  record <- as.list(values)
  names(record) <- keys
  canonical <- .artifact_canonical_path(path)
  if (!identical(record$PATH, canonical) ||
      !identical(record$RUN_ID, context$run_id) ||
      !identical(record$PRODUCER, context$producer) ||
      !identical(record$STATE, "PUBLISHED") ||
      !grepl("^[0-9a-f]{32}$", record$MD5) ||
      !grepl("^[1-9][0-9]*$", record$SIZE)) {
    stop("Artifact record binding is invalid: ", record_path)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0 ||
      !identical(record$SIZE, as.character(info$size))) {
    stop("Artifact record SIZE mismatch: ", path)
  }
  sidecar <- .artifact_sidecar(path, allow_canonical = TRUE)
  if (is.null(sidecar) || !identical(sidecar$MD5, record$MD5) ||
      !identical(sidecar$SIZE, record$SIZE)) {
    stop("Artifact record checksum sidecar mismatch: ", path)
  }
  record
}

artifact_validate_record <- function(path, producer = NULL, run_id = NULL) {
  .artifact_record_from_file(path, producer = producer, run_id = run_id)
}

artifact_record_for_load <- function(path, producer = NULL, run_id = NULL) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (!is.null(context)) {
    record_path <- artifact_record_path(
      path, context$run_id, runs_root = context$runs_root
    )
    if (file.exists(record_path)) {
      record <- .artifact_record_from_file(
        path, producer = context$producer, run_id = context$run_id
      )
      if (is.null(record)) {
        stop("Artifact checksum validation failed: ", path)
      }
      sidecar <- .artifact_sidecar(path, verify = TRUE)
      if (is.null(sidecar) ||
          !identical(sidecar$MD5, record$MD5) ||
          !identical(sidecar$SIZE, record$SIZE)) {
        stop("Artifact checksum validation failed: ", path)
      }
      return(record)
    }
  }
  sidecar <- .artifact_sidecar(path, verify = TRUE)
  if (is.null(sidecar)) {
    stop("Artifact checksum validation failed: ", path)
  }
  sidecar
}

artifact_write_record <- function(
  path, producer = NULL, run_id = NULL, md5 = NULL, size = NULL
) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (is.null(context)) return(invisible(NULL))
  canonical <- .artifact_canonical_path(path)
  sidecar <- .artifact_sidecar(path, allow_canonical = TRUE)
  if (is.null(sidecar)) {
    stop("Cannot publish artifact record without a valid sidecar: ", path)
  }
  if (is.null(md5)) md5 <- sidecar$MD5
  if (is.null(size)) size <- sidecar$SIZE
  if (!identical(tolower(as.character(md5)), sidecar$MD5) ||
      !identical(as.character(size), sidecar$SIZE)) {
    stop("Artifact record checksum does not match sidecar: ", path)
  }
  record_path <- artifact_record_path(
    canonical, context$run_id, runs_root = context$runs_root
  )
  dir.create(dirname(record_path), showWarnings = FALSE, recursive = TRUE)
  temporary <- paste0(record_path, ".tmp.", Sys.getpid())
  writeLines(c(
    paste0("PATH=", canonical),
    paste0("SIZE=", sidecar$SIZE),
    paste0("MD5=", sidecar$MD5),
    paste0("RUN_ID=", context$run_id),
    paste0("PRODUCER=", context$producer),
    "STATE=PUBLISHED"
  ), temporary)
  if (!file.rename(temporary, record_path)) {
    if (file.exists(temporary)) unlink(temporary)
    stop("Could not atomically install artifact record: ", record_path)
  }
  invisible(list(
    PATH = canonical, SIZE = sidecar$SIZE, MD5 = sidecar$MD5,
    RUN_ID = context$run_id, PRODUCER = context$producer,
    STATE = "PUBLISHED"
  ))
}

artifact_checksum_ok <- function(file, producer = NULL, run_id = NULL) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (!is.null(context)) {
    record_path <- artifact_record_path(
      file, context$run_id, runs_root = context$runs_root
    )
    if (file.exists(record_path)) {
      return(isTRUE(tryCatch(
        !is.null(.artifact_record_from_file(
          file, producer = context$producer, run_id = context$run_id
        )),
        error = function(e) FALSE
      )))
    }
  }
  isTRUE(!is.null(.artifact_sidecar(file, verify = TRUE)))
}

read_rds_checked <- function(path, producer = NULL, run_id = NULL) {
  artifact_record_for_load(path, producer = producer, run_id = run_id)
  readRDS(path)
}

read_feather_checked <- function(path, producer = NULL, run_id = NULL) {
  artifact_record_for_load(path, producer = producer, run_id = run_id)
  arrow::read_feather(path)
}
# ---------------------------------------------------------------------------
# Runtime metadata for cached benchmark artifacts. The JSON and its checksum
# are published only after the output artifact itself has been installed.
# ---------------------------------------------------------------------------
runtime_metadata_path <- function(artifact_path) {
  paste0(artifact_path, ".runtime.json")
}

runtime_metadata_checksum_path <- function(artifact_path) {
  paste0(runtime_metadata_path(artifact_path), ".md5")
}

runtime_require_string <- function(value, label) {
  if (!is.character(value) || length(value) != 1L ||
      is.na(value) || !nzchar(value)) {
    stop(label, " must be one non-empty string")
  }
  value
}

runtime_require_nonnegative_number <- function(value, label) {
  if (!is.numeric(value) || length(value) != 1L ||
      is.na(value) || !is.finite(value) || value < 0) {
    stop(label, " must be one finite nonnegative number")
  }
  as.numeric(value)
}

runtime_normalize_mem <- function(value, label = "mem_GB") {
  if (is.null(value)) return(NULL)
  if (!is.numeric(value) || length(value) != 1L) {
    stop(label, " must be one finite nonnegative number or null")
  }
  # NA_real_ is the in-memory representation of the JSON null used when
  # /proc is unavailable; NaN is nonfinite and must fail closed.
  if (is.na(value)) {
    if (is.nan(value)) stop(label, " must be finite or null")
    return(NULL)
  }
  if (!is.finite(value) || value < 0) {
    stop(label, " must be one finite nonnegative number or null")
  }
  as.numeric(value)
}

# Strictly validate an MD5/SIZE/PATH sidecar and return its canonical record.
# Runtime metadata uses this for both the published RDS and the JSON sidecar;
# unlike artifact_checksum_ok(), malformed or extra sidecar rows are errors.
runtime_validate_checksum_sidecar <- function(
  file, description = "artifact", expected_md5 = NULL, expected_size = NULL
) {
  file <- runtime_require_string(file, paste0(description, " path"))
  if (!file.exists(file)) stop(description, " is missing: ", file)
  info <- file.info(file)
  if (is.na(info$size) || info$size <= 0) {
    stop(description, " is missing or empty: ", file)
  }
  sidecar <- paste0(file, ".md5")
  if (!file.exists(sidecar)) {
    stop(description, " checksum sidecar is missing: ", sidecar)
  }
  sidecar_info <- file.info(sidecar)
  if (is.na(sidecar_info$size) || sidecar_info$size <= 0) {
    stop(description, " checksum sidecar is missing or empty: ", sidecar)
  }
  lines <- tryCatch(
    readLines(sidecar, warn = FALSE),
    error = function(error) {
      stop("Could not read ", description, " checksum sidecar: ", sidecar)
    }
  )
  if (length(lines) != 3L || any(!nzchar(lines)) ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    stop(description, " checksum sidecar has the wrong schema: ", sidecar)
  }
  recorded_md5 <- sub("^MD5=", "", lines[[1L]])
  recorded_size <- sub("^SIZE=", "", lines[[2L]])
  recorded_path <- sub("^PATH=", "", lines[[3L]])
  if (!grepl("^[0-9a-f]{32}$", recorded_md5, perl = TRUE)) {
    stop(description, " checksum sidecar has a malformed MD5: ", sidecar)
  }
  if (!grepl("^[0-9]+$", recorded_size, perl = TRUE) ||
      !nzchar(recorded_path)) {
    stop(description, " checksum sidecar has malformed SIZE/PATH: ", sidecar)
  }
  if (!identical(recorded_path, file)) {
    stop(description, " checksum sidecar PATH mismatch: ", sidecar)
  }
  if (!identical(recorded_size, as.character(info$size))) {
    stop(description, " checksum sidecar SIZE mismatch: ", file)
  }
  if (!is.null(expected_md5)) {
    if (!identical(tolower(as.character(expected_md5)), recorded_md5)) {
      stop(description, " checksum sidecar MD5 mismatch: ", file)
    }
  } else {
    actual_md5 <- tolower(unname(tools::md5sum(file)))
    if (length(actual_md5) != 1L || is.na(actual_md5) ||
        !identical(recorded_md5, actual_md5)) {
      stop(description, " checksum sidecar MD5 mismatch: ", file)
    }
  }
  if (!is.null(expected_size) &&
      !identical(as.character(expected_size), recorded_size)) {
    stop(description, " checksum sidecar SIZE mismatch: ", file)
  }
  list(MD5 = recorded_md5, SIZE = recorded_size, PATH = recorded_path)
}

validate_runtime_metadata <- function(
  artifact_path,
  dataset,
  method,
  producer = NULL,
  run_id = NULL,
  artifact_record = NULL
) {
  artifact_path <- runtime_require_string(artifact_path, "artifact path")
  dataset <- runtime_require_string(dataset, "dataset")
  method <- runtime_require_string(method, "method")
  if (is.null(artifact_record)) {
    artifact_record <- artifact_record_for_load(
      artifact_path, producer = producer, run_id = run_id
    )
  }
  runtime_validate_checksum_sidecar(
    artifact_path, "RDS artifact",
    expected_md5 = artifact_record[["MD5"]],
    expected_size = artifact_record[["SIZE"]]
  )
  metadata_path <- runtime_metadata_path(artifact_path)
  metadata_record <- artifact_record_for_load(
    metadata_path, producer = producer, run_id = run_id
  )
  runtime_validate_checksum_sidecar(
    metadata_path, "runtime metadata",
    expected_md5 = metadata_record[["MD5"]],
    expected_size = metadata_record[["SIZE"]]
  )
  parse_error <- NULL
  payload <- tryCatch(
    jsonlite::fromJSON(metadata_path, simplifyVector = FALSE),
    error = function(error) {
      parse_error <<- conditionMessage(error)
      NULL
    }
  )
  if (!is.null(parse_error) || !is.list(payload)) {
    stop("Runtime metadata is not a JSON object: ", metadata_path)
  }
  required <- c(
    "schema_version", "artifact_path", "artifact_md5", "dataset",
    "method", "time_secs", "mem_GB"
  )
  if (!identical(sort(names(payload)), sort(required))) {
    stop("Runtime metadata has the wrong fields: ", metadata_path)
  }
  schema_version <- payload[["schema_version"]]
  if (!is.integer(schema_version) || !identical(schema_version, 1L)) {
    stop("Runtime metadata schema_version is invalid: ", metadata_path)
  }
  if (!is.character(payload[["artifact_path"]]) ||
      length(payload[["artifact_path"]]) != 1L ||
      is.na(payload[["artifact_path"]]) ||
      !identical(payload[["artifact_path"]], artifact_path)) {
    stop("Runtime metadata artifact_path mismatch: ", metadata_path)
  }
  artifact_md5 <- payload[["artifact_md5"]]
  if (!is.character(artifact_md5) || length(artifact_md5) != 1L ||
      is.na(artifact_md5) ||
      !grepl("^[0-9a-f]{32}$", artifact_md5, perl = TRUE) ||
      !identical(artifact_md5, artifact_record[["MD5"]])) {
    stop("Runtime metadata artifact_md5 mismatch: ", metadata_path)
  }
  if (!is.character(payload[["dataset"]]) ||
      length(payload[["dataset"]]) != 1L || is.na(payload[["dataset"]]) ||
      !identical(payload[["dataset"]], dataset)) {
    stop("Runtime metadata dataset mismatch: ", metadata_path)
  }
  if (!is.character(payload[["method"]]) ||
      length(payload[["method"]]) != 1L || is.na(payload[["method"]]) ||
      !identical(payload[["method"]], method)) {
    stop("Runtime metadata method mismatch: ", metadata_path)
  }
  time_secs <- runtime_require_nonnegative_number(
    payload[["time_secs"]], "Runtime metadata time_secs"
  )
  mem_gb <- runtime_normalize_mem(payload[["mem_GB"]])
  payload[["schema_version"]] <- 1L
  payload[["time_secs"]] <- time_secs
  payload["mem_GB"] <- list(mem_gb)
  payload
}
write_runtime_metadata <- function(
  artifact_path,
  dataset,
  method,
  time_secs,
  mem_gb = NA_real_,
  producer = NULL,
  run_id = NULL
) {
  artifact_path <- runtime_require_string(artifact_path, "artifact path")
  dataset <- runtime_require_string(dataset, "dataset")
  method <- runtime_require_string(method, "method")
  time_secs <- runtime_require_nonnegative_number(time_secs, "time_secs")
  mem_gb <- runtime_normalize_mem(mem_gb)
  artifact_record <- artifact_record_for_load(
    artifact_path, producer = producer, run_id = run_id
  )
  runtime_validate_checksum_sidecar(
    artifact_path, "RDS artifact",
    expected_md5 = artifact_record[["MD5"]],
    expected_size = artifact_record[["SIZE"]]
  )
  metadata_path <- runtime_metadata_path(artifact_path)
  checksum_path <- runtime_metadata_checksum_path(artifact_path)
  payload <- list(
    schema_version = 1L,
    artifact_path = artifact_path,
    artifact_md5 = artifact_record[["MD5"]],
    dataset = dataset,
    method = method,
    time_secs = time_secs,
    mem_GB = mem_gb
  )
  metadata_json <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = 17
  )
  dir.create(dirname(metadata_path), showWarnings = FALSE, recursive = TRUE)
  metadata_tmp <- paste0(metadata_path, ".tmp.", Sys.getpid())
  checksum_tmp <- paste0(checksum_path, ".tmp.", Sys.getpid())
  metadata_backup <- paste0(metadata_path, ".previous.", Sys.getpid())
  checksum_backup <- paste0(checksum_path, ".previous.", Sys.getpid())
  had_metadata <- file.exists(metadata_path)
  had_checksum <- file.exists(checksum_path)
  metadata_installed <- FALSE
  checksum_installed <- FALSE
  metadata_digest <- NULL
  metadata_size <- NULL
  restore <- function() {
    if (metadata_installed && file.exists(metadata_path)) unlink(metadata_path)
    if (had_metadata && file.exists(metadata_backup)) {
      file.rename(metadata_backup, metadata_path)
    }
    if (checksum_installed && file.exists(checksum_path)) unlink(checksum_path)
    if (had_checksum && file.exists(checksum_backup)) {
      file.rename(checksum_backup, checksum_path)
    }
    if (!had_metadata && file.exists(metadata_path)) unlink(metadata_path)
    if (!had_checksum && file.exists(checksum_path)) unlink(checksum_path)
  }
  on.exit({
    for (temporary in c(
      metadata_tmp, checksum_tmp, metadata_backup, checksum_backup
    )) {
      if (file.exists(temporary)) unlink(temporary)
    }
  }, add = TRUE)
  tryCatch({
    writeLines(as.character(metadata_json), metadata_tmp, useBytes = TRUE)
    if (!file.exists(metadata_tmp) || file.info(metadata_tmp)$size <= 0) {
      stop("Empty runtime metadata temporary file: ", metadata_tmp)
    }
    if (had_metadata && !isTRUE(file.link(metadata_path, metadata_backup))) {
      stop("Could not preserve existing runtime metadata: ", metadata_path)
    }
    if (had_checksum && !isTRUE(file.link(checksum_path, checksum_backup))) {
      stop("Could not preserve existing runtime metadata checksum: ",
           checksum_path)
    }
    if (!file.rename(metadata_tmp, metadata_path)) {
      stop("Could not atomically install runtime metadata: ", metadata_path)
    }
    metadata_installed <- TRUE
    metadata_digest <- tolower(unname(tools::md5sum(metadata_path)))
    metadata_size <- as.character(file.info(metadata_path)$size)
    writeLines(c(
      paste0("MD5=", metadata_digest),
      paste0("SIZE=", metadata_size),
      paste0("PATH=", metadata_path)
    ), checksum_tmp)
    if (!file.rename(checksum_tmp, checksum_path)) {
      stop("Could not atomically install runtime metadata checksum: ",
           checksum_path)
    }
    checksum_installed <- TRUE
  }, error = function(error) {
    restore()
    stop(error)
  })
  artifact_write_record(
    metadata_path, producer = producer, run_id = run_id,
    md5 = metadata_digest, size = metadata_size
  )
  invisible(NULL)
}

replay_runtime_metadata <- function(
  artifact_path,
  dataset,
  method,
  log_file,
  producer = NULL,
  run_id = NULL,
  artifact_record = NULL
) {
  payload <- validate_runtime_metadata(
    artifact_path, dataset, method,
    producer = producer, run_id = run_id,
    artifact_record = artifact_record
  )
  mem_gb <- payload[["mem_GB"]]
  if (is.null(mem_gb)) mem_gb <- NA_real_
  log_exec_row(
    dataset,
    method,
    payload[["time_secs"]],
    log_file,
    mem_gb = mem_gb
  )
  invisible(payload)
}


# Peak resident set size of the current R process in GB, mirroring the
# python worker's peak_rss_gb() (getrusage().ru_maxrss: KB on Linux, bytes
# on macOS). On Linux, VmHWM from /proc/self/status is the process peak-RSS
# equivalent of ru_maxrss (no extra packages, base R only). Off-Linux (and
# when /proc is unavailable) returns NA_real_.
# Same monotonic-cumulative semantics as python: call it at each combo's
# completion; combos running earlier report the least bloated peak.
peak_rss_gb <- function() {
  if (.Platform$OS.type != "unix") return(NA_real_)
  status_file <- "/proc/self/status"
  if (!file.exists(status_file)) return(NA_real_)
  hwm <- grep("^VmHWM:", readLines(status_file, warn = FALSE), value = TRUE)
  if (length(hwm) == 0) return(NA_real_)
  kb <- suppressWarnings(
    as.numeric(sub("^VmHWM:\\s*([0-9]+)\\s*kB\\s*$", "\\1", hwm[1]))
  )
  if (is.na(kb)) return(NA_real_)
  return(kb / 1024^2)
}

# Append/overwrite one (dataset, method) row in the per-task exec log feather.
# The write and its checksum are installed atomically so concurrent reruns
# cannot leave a partial file that a later worker treats as complete.
write_feather_atomic <- function(
  df, path, producer = NULL, run_id = NULL
) {
  if (nrow(df) <= 0L ||
      !all(c("dataset", "method", "time_secs", "mem_GB") %in% names(df))) {
    stop("Execution log has an empty or incomplete schema: ", path)
  }
  if (anyNA(df[["dataset"]]) || anyNA(df[["method"]])) {
    stop("Execution log has NA identifiers: ", path)
  }
  identifiers <- paste(df[["dataset"]], df[["method"]], sep = "\r")
  if (anyDuplicated(identifiers)) {
    stop("Execution log has duplicate identifiers: ", path)
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- file.path(dirname(path), paste0(".", basename(path), ".tmp.", Sys.getpid()))
  checksum_tmp <- paste0(path, ".md5.tmp.", Sys.getpid())
  sidecar <- paste0(path, ".md5")
  backup <- paste0(path, ".previous.", Sys.getpid())
  sidecar_backup <- paste0(sidecar, ".previous.", Sys.getpid())
  had_file <- file.exists(path)
  had_sidecar <- file.exists(sidecar)
  installed <- FALSE
  sidecar_installed <- FALSE
  digest <- NULL
  size <- NULL
  restore <- function() {
    if (installed && file.exists(path)) unlink(path)
    if (had_file && file.exists(backup)) file.rename(backup, path)
    if (sidecar_installed && file.exists(sidecar)) unlink(sidecar)
    if (had_sidecar && file.exists(sidecar_backup)) file.rename(sidecar_backup, sidecar)
    if (!had_file && file.exists(path)) unlink(path)
    if (!had_sidecar && file.exists(sidecar)) unlink(sidecar)
  }
  on.exit({
    for (temporary in c(tmp, checksum_tmp, backup, sidecar_backup)) {
      if (file.exists(temporary)) unlink(temporary)
    }
  }, add = TRUE)
  tryCatch({
    arrow::write_feather(df, tmp)
    if (!file.exists(tmp) || file.info(tmp)$size <= 0) {
      stop("Atomic Feather write produced an empty file: ", tmp)
    }
    if (had_file && !isTRUE(file.link(path, backup))) {
      stop("Could not preserve existing Feather artifact: ", path)
    }
    if (had_sidecar && !isTRUE(file.link(sidecar, sidecar_backup))) {
      stop("Could not preserve existing Feather checksum: ", path)
    }
    if (!file.rename(tmp, path)) stop("Could not atomically install Feather: ", path)
    installed <- TRUE
    digest <- tolower(unname(tools::md5sum(path)))
    size <- as.character(file.info(path)$size)
    writeLines(c(
      paste0("MD5=", digest),
      paste0("SIZE=", size),
      paste0("PATH=", path)
    ), checksum_tmp)
    if (!file.rename(checksum_tmp, sidecar)) {
      stop("Could not atomically install Feather checksum: ", path)
    }
    sidecar_installed <- TRUE
  }, error = function(error) {
    restore()
    stop(error)
  })
  artifact_write_record(path, producer = producer, run_id = run_id,
                        md5 = digest, size = size)
  invisible(NULL)
}

# Append/overwrite one (dataset, method) row in the per-task exec log feather.
log_exec_row <- function(
  dataset, method, time_secs, log_file, mem_gb = NA_real_,
  producer = NULL, run_id = NULL
) {
  if (is.null(log_file) || is.na(log_file) || log_file == "") {
    return(invisible(NULL))
  }
  if (is.null(producer)) {
    producer <- Sys.getenv(
      "ECODA_EXECUTION_LOG_PRODUCER", unset = "stage5_execution_log"
    )
  }
  if (is.null(mem_gb)) mem_gb <- NA_real_
  new_row <- data.frame(
    dataset = as.character(dataset),
    method = as.character(method),
    time_secs = as.numeric(time_secs),
    mem_GB = mem_gb,
    stringsAsFactors = FALSE
  )
  if (file.exists(log_file) && file.info(log_file)$size > 0) {
    if (!artifact_checksum_ok(
      log_file, producer = producer, run_id = run_id
    )) {
      stop("Execution log checksum or artifact record validation failed: ",
           log_file)
    }
    df_existing <- read_feather_checked(
      log_file, producer = producer, run_id = run_id
    )
    df_existing <- df_existing[
      !(df_existing[["dataset"]] == dataset &
        df_existing[["method"]] == method),
      ,
      drop = FALSE
    ]
    df_final <- rbind(df_existing, new_row)
  } else {
    df_final <- new_row
  }
  write_feather_atomic(
    df_final, log_file, producer = producer, run_id = run_id
  )
  return(invisible(NULL))
}

# ============================================================================
# Pipeline B worker driver (shared by 1.1.1_run_transformation_analysis.R and
# 1.1.1_run_zeroimp_analysis.R): parses the CLI args, reads datasets.json and
# the preprocessed view h5ad (obs-only backed read, no counts matrix), builds
# the per-sample cell-type composition table + labels, runs the analysis and
# saves <ds><out_suffix>.rds atomically with one exec-log row per dataset.
# ============================================================================
run_ct_comps_analysis_worker <- function(
  analysis_label,  # "trans" / "zeroimp" (messages)
  run_fun,         # function(ct_comps, labels)
  out_suffix,      # "_trans" / "_zeroimp"
  log_method       # "trans_analysis" / "zeroimp_analysis"
) {
  raw_args <- commandArgs(trailingOnly = TRUE)
  args <- parse_flags(raw_args)

  for (req in c("config_path", "ds_name", "view", "input_dir",
                "output_dir", "log_file")) {
    if (is.null(args[[req]]) || identical(args[[req]], TRUE)) {
      stop("Missing required --", req, " argument")
    }
  }
  force <- isTRUE(args[["force"]]) || identical(args[["force"]], "TRUE")

  artifact_producer <- paste0("stage5_", analysis_label)
  artifact_run_id <- Sys.getenv(
    "ECODA_ARTIFACT_PRODUCER_RUN_ID",
    unset = ""
  )
  if (!nzchar(artifact_run_id)) {
    artifact_run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
  }
  if (!nzchar(artifact_run_id)) artifact_run_id <- NULL

  config <- read_datasets_json(args$config_path, view = args$view)
  ds <- args$ds_name
  entry <- config[[ds]]
  if (is.null(entry)) {
    stop("Dataset '", ds, "' not found in ", args$config_path)
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- file.path(args$output_dir, paste0(ds, out_suffix, ".rds"))
  if (artifact_checksum_ok(
    out_file,
    producer = artifact_producer,
    run_id = artifact_run_id
  ) && !force) {
    replay_runtime_metadata(
      out_file,
      ds,
      log_method,
      args$log_file,
      producer = artifact_producer,
      run_id = artifact_run_id
    )
    message(analysis_label, " results already exist and passed checksum validation: ", out_file)
    return(invisible(NULL))
  }

  h5ad_path <- get_h5ad_path(config, ds, args$view, args$input_dir)
  if (!file.exists(h5ad_path)) {
    stop("Input h5ad not found: ", h5ad_path)
  }

  ad <- import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(h5ad_path, backed = "r")
  obs <- py_to_r(adata$obs)

  sample_col <- "Sample"
  if (!sample_col %in% colnames(obs)) {
    stop(sample_col, " not found in obs columns of ", h5ad_path)
  }
  ct_col <- entry$cell_type_high_res
  if (is.null(ct_col)) {
    stop("Dataset '", ds, "' has no cell_type_high_res column in datasets.json")
  }
  if (!ct_col %in% colnames(obs)) {
    stop("Cell type column '", ct_col, "' not found in obs columns of ",
         h5ad_path)
  }
  label_col <- entry$label_col
  if (!label_col %in% colnames(obs)) {
    stop("Label column '", label_col, "' not found in obs columns of ",
         h5ad_path)
  }

  # Cell-type composition per sample (get_ct_comp_df: rows = samples, cols =
  # cell types; rowSums != 0 filter). obs["Sample"] is already the
  # standardized sample column; a plain data.frame keeps the dplyr verbs in
  # run_transformation_analysis / run_zeroimp_analysis working.
  ct_comps <- get_ct_comp_df(obs, sample_col, ct_col)

  # Labels: per-sample slice(1) of label_col, names = Sample
  # (get_labels-equivalent)
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
  labels <- as.factor(metadata[[label_col]])
  names(labels) <- metadata[[sample_col]]

  time_secs <- exec_time(res <- run_fun(ct_comps, labels))
  time_secs_numeric <- as.numeric(time_secs, units = "secs")
  mem_gb <- peak_rss_gb()
  save_rds_atomic(
    res,
    out_file,
    producer = artifact_producer,
    run_id = artifact_run_id
  )
  write_runtime_metadata(
    out_file,
    ds,
    log_method,
    time_secs_numeric,
    mem_gb = mem_gb,
    producer = artifact_producer,
    run_id = artifact_run_id
  )
  log_exec_row(
    ds,
    log_method,
    time_secs_numeric,
    args$log_file,
    mem_gb = mem_gb,
    run_id = artifact_run_id
  )
  message("Saved: ", out_file, " (",
          round(time_secs_numeric, 1), "s)")
  message("--- ", analysis_label, " analysis for ", ds, " complete ---")
}
