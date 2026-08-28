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
  paths <- file.path(
    pseudobulk_dir,
    paste0(cache_stem, "_pseudobulk_", PB_VARIANT_NAMES, ".rds")
  )
  PB_VARIANT_NAMES[!vapply(paths, artifact_checksum_ok, logical(1))]
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

# Build the full Seurat object from the validated raw counts layer with the
# benchmark obsm PCA embeddings attached as reductions named
# pca_benchmark_analysis_hvg{n} (via get_seurat_obj_from_h5ad). All samples of
# the dataset are included; obs is read once by the caller.
load_benchmark_seurat <- function(
  adata,
  obs,
  sample_col = "Sample",
  fetch_embedding = c(
    "X_pca_benchmark_analysis_hvg1000",
    "X_pca_benchmark_analysis_hvg2000",
    "X_pca_benchmark_analysis_hvg3000"
  )
) {
  all_samples <- unique(obs[[sample_col]])
  layer_keys <- py_to_r(import_builtins(convert = FALSE)$list(
    adata$layers$keys()
  ))
  if (!"counts" %in% layer_keys) {
    stop(
      "load_benchmark_seurat: validated h5ad is missing layers['counts']; ",
      "refusing to use log-normalized X as DESeq2 counts."
    )
  }
  seurat <- get_seurat_obj_from_h5ad(
    adata, obs, all_samples,
    sample_colname = sample_col,
    counts_layer = "counts",
    fetch_embedding = fetch_embedding
  )
  return(seurat)
}

# Load the precomputed pseudobulk variants from pseudobulks/ (produced by the
# prepare_pseudobulk worker); variants missing on disk (e.g. --methods mofa
# without the prep array on a stale scratch) are computed on the fly via
# prepare_pseudobulks_hpc (only the missing subset), saved atomically and
# logged as prepare_pseudobulk_<variant> rows. `seurat` (counts) is required
# only when the on-the-fly fallback triggers; callers that cannot afford the
# full object may pre-check pb_variants_missing() and build it lazily.
# Returns per-variant list(pb, time_secs).
load_pb_variants <- function(
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
  variants = PB_VARIANT_NAMES
) {
  target_variants <- unique(as.character(variants))
  if (!all(target_variants %in% PB_VARIANT_NAMES)) {
    stop("Unknown pseudobulk variant requested: ",
         paste(setdiff(target_variants, PB_VARIANT_NAMES), collapse = ", "))
  }
  cache_paths <- file.path(
    pseudobulk_dir,
    paste0(cache_stem, "_pseudobulk_", target_variants, ".rds")
  )
  missing <- target_variants[
    force | !vapply(cache_paths, artifact_checksum_ok, logical(1))
  ]
  loaded <- target_variants[!target_variants %in% missing]
  variants_out <- list()
  for (v in loaded) {
    variants_out[[v]] <- readRDS(
      file.path(pseudobulk_dir, paste0(cache_stem, "_pseudobulk_", v, ".rds"))
    )
  }
  if (length(missing) > 0) {
    if (is.null(seurat)) {
      stop("Pseudobulk variant(s) missing in ", pseudobulk_dir,
           " (", paste(missing, collapse = ", "),
           ") and no Seurat object was provided for the on-the-fly fallback.")
    }
    message("Pseudobulk variant(s) missing in ", pseudobulk_dir,
            ", computing on the fly: ", paste(missing, collapse = ", "))
    computed <- prepare_pseudobulks_hpc(
      seurat,
      sample_col = sample_col,
      hvg_rank_genes = hvg_rank_genes,
      variants = missing,
      batch_col = batch_col,
      blind = blind,
      correct_batch = correct_batch
    )
    for (v in names(computed)) {
      save_rds_atomic(
        computed[[v]],
        file.path(pseudobulk_dir, paste0(cache_stem, "_pseudobulk_", v, ".rds"))
      )
      log_exec_row(ds, paste0("prepare_pseudobulk_", v),
                   computed[[v]]$time_secs, log_file)
      variants_out[[v]] <- computed[[v]]
    }
  }
  return(variants_out)
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


# Atomic RDS write plus a sidecar used by idempotency checks.
save_rds_atomic <- function(object, file) {
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
    writeLines(c(
      paste0("MD5=", unname(tools::md5sum(file))),
      paste0("SIZE=", file.info(file)$size),
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
  invisible(NULL)
}

artifact_checksum_ok <- function(file) {
  sidecar <- paste0(file, ".md5")
  if (!file.exists(file) || file.info(file)$size <= 0 || !file.exists(sidecar)) return(FALSE)
  lines <- readLines(sidecar, warn = FALSE)
  md5 <- sub("^MD5=", "", lines[grepl("^MD5=", lines)][1])
  size <- sub("^SIZE=", "", lines[grepl("^SIZE=", lines)][1])
  recorded <- sub("^PATH=", "", lines[grepl("^PATH=", lines)][1])
  identical(recorded, file) &&
    identical(md5, unname(tools::md5sum(file))) &&
    identical(size, as.character(file.info(file)$size))
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
write_feather_atomic <- function(df, path) {
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
    digest <- unname(tools::md5sum(path))
    writeLines(c(
      paste0("MD5=", digest),
      paste0("SIZE=", file.info(path)$size),
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
  invisible(NULL)
}

# Append/overwrite one (dataset, method) row in the per-task exec log feather.
log_exec_row <- function(dataset, method, time_secs, log_file,
                         mem_gb = NA_real_) {
  if (is.null(log_file) || is.na(log_file) || log_file == "") {
    return(invisible(NULL))
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
    df_existing <- arrow::read_feather(log_file)
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
  write_feather_atomic(df_final, log_file)
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

  config <- read_datasets_json(args$config_path, view = args$view)
  ds <- args$ds_name
  entry <- config[[ds]]
  if (is.null(entry)) {
    stop("Dataset '", ds, "' not found in ", args$config_path)
  }

  h5ad_path <- get_h5ad_path(config, ds, args$view, args$input_dir)
  if (!file.exists(h5ad_path)) {
    stop("Input h5ad not found: ", h5ad_path)
  }
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- file.path(args$output_dir, paste0(ds, out_suffix, ".rds"))
  if (artifact_checksum_ok(out_file) && !force) {
    message(analysis_label, " results already exist and passed checksum validation: ", out_file)
    quit(save = "no", status = 0)
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
  save_rds_atomic(res, out_file)
  log_exec_row(ds, log_method, as.numeric(time_secs, units = "secs"),
               args$log_file, mem_gb = peak_rss_gb())
  message("Saved: ", out_file, " (",
          round(as.numeric(time_secs, units = "secs"), 1), "s)")
  message("--- ", analysis_label, " analysis for ", ds, " complete ---")
}
