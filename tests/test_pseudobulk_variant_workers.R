# Execute the actual R entrypoints with metadata/count/score boundaries stubbed.
# No packages, H5AD expression matrices, or scheduler access are required.
root <- normalizePath(".")
worker_dir <- file.path(root, "src/5_run_benchmark_methods/run_r_sample_embedding_methods")
assert_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
local({
  prior <- Sys.getenv(c("PROJECT_ROOT", "ANALYSIS_VARIANT", "ANALYSIS_PASS"), unset = NA_character_)
  on.exit(for (name in names(prior)) {
    if (is.na(prior[[name]])) Sys.unsetenv(name) else
      do.call(Sys.setenv, setNames(list(prior[[name]]), name))
  })
  Sys.setenv(PROJECT_ROOT = root)
  Sys.unsetenv(c("ANALYSIS_VARIANT", "ANALYSIS_PASS"))
  tmp <- tempfile("pseudobulk-worker-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  input <- file.path(tmp, "input.h5ad")
  writeLines("metadata reader is stubbed", input)
  execute <- function(preparation, selector = c("--pseudobulk_variants", "hvg2000_bl"),
                      method = "pseudobulk", pass = NULL) {
    env <- new.env(parent = globalenv())
    sys.source(file.path(root, "src/utils/pseudobulk.R"), envir = env)
    sys.source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"), envir = env)
    script <- file.path(worker_dir, if (preparation) "1.1.1_prepare_pseudobulk.R" else
      "1.1.1_run_benchmark_methods_r.R")
    argv <- c("--config_path", "stub.json", "--ds_name", "Toy",
              "--view", "benchmark_analysis", "--method", method,
              "--input_dir", tmp, "--results_dir", file.path(tmp, "results"),
              "--pseudobulk_dir", file.path(tmp, "pseudobulks"),
              "--log_file", file.path(tmp, "timing.feather"), selector,
              if (!is.null(pass)) c("--analysis_pass", pass))
    env$source <- function(...) invisible(NULL)
    env$attached <- character()
    env$library <- function(package) {
      name <- as.character(substitute(package))
      if (!name %in% c("reticulate", "dplyr", "GloScope")) stop("unexpected method-specific package loading")
      env$attached <- c(env$attached, name)
    }
    env$commandArgs <- function(trailingOnly = FALSE) {
      if (trailingOnly) argv else paste0("--file=", script)
    }
    env$read_datasets_json <- function(...) list(Toy = list(
      label_col = "biology", cell_type_low_res = "LR", cell_type_high_res = "HR"
    ))
    env$get_h5ad_path <- function(...) input
    env$metadata_columns <- NULL
    env$load_h5ad_pseudobulk_metadata <- function(..., metadata_columns) {
      env$metadata_columns <- metadata_columns
      list(obs = data.frame(Sample = c("s1", "s2", "s3", "s4"),
                            biology = c("A", "B", "A", "B")),
           hvg_rank_genes = "unused_rank")
    }
    env$load_h5ad_author_hr_metadata <- function(h5ad_path, obs_columns, view) {
      env$metadata_columns <- obs_columns
      data.frame(Sample = paste0("s", 1:4), biology = c("A", "B", "A", "B"), HR = c("T", "B", "T", "M"))
    }
    env$load_h5ad_counts_free <- function(...) stop("unselected embedding loader reached")
    env$py_to_r <- identity
    env$get_hvg_rank_genes <- env$rename_leiden_cols <- env$load_composition_pb_variants <-
      function(...) stop("unselected composition prerequisite reached")
    env$collapse_sample_metadata <- function(obs, ...) obs
    env$artifact_checksum_ok <- function(...) TRUE # Old aggregate must not short-circuit a shard.
    env$read_rds_checked <- function(...) stop("must not read an aggregate or unselected cache")
    env$.pb_variant_cache_valid <- function(...) FALSE
    env$prepared <- env$loaded <- env$scored <- env$published <- list()
    value <- list(pb = matrix(1, nrow = 4L, ncol = 1L), time_secs = 1)
    env$prepare_pseudobulks_hpc <- function(..., variants) {
      env$prepared <- variants
      setNames(rep(list(value), length(variants)), variants)
    }
    env$load_pb_variants <- function(..., variants) {
      env$loaded <- variants
      setNames(rep(list(value), length(variants)), variants)
    }
    env$run_pseudobulk_hpc <- function(..., variants) {
      env$scored <- variants
      setNames(rep(list(value), length(variants)), paste0("Pseudobulk_", variants))
    }
    env$run_composition_methods_hpc <- function(labels, metadata, pca_emb, pb_hvg2000, obs, ..., variants) {
      stopifnot(is.null(pca_emb), is.null(pb_hvg2000), identical(obs$HR, metadata$HR))
      env$scored <- variants
      list()
    }
    env$save_rds_atomic <- function(object, file, ...) {
      env$published[[file]] <- object
    }
    env$emit_pseudobulk_timing_rows <- function(...) invisible(NULL)
    env$log_exec_row <- function(...) invisible(NULL)
    env$import <- env$load_benchmark_seurat <- function(...) stop("heavy boundary reached")
    sys.source(script, envir = env)
    env
  }
  prep <- execute(TRUE)
  stopifnot(identical(prep$prepared, "hvg2000_bl"),
            identical(basename(names(prep$published)), "Toy_pseudobulk_hvg2000_bl.rds"))
  result <- execute(FALSE)
  stopifnot(identical(result$loaded, "hvg2000_bl"),
            identical(result$scored, "hvg2000_bl"),
            !any(c("LR", "HR") %in% result$metadata_columns),
            length(result$published) == 0L) # Standalone pipeline owns publication.
  for (selection in c("CT_HR_hvg2000,CT_HR_hvg500", "hvg2000_bl,CT_HR_hvg2000,CT_HR_hvg500")) {
    ct <- execute(FALSE, c("--pseudobulk_variants", selection))
    stopifnot(identical(ct$scored, strsplit(selection, ",", fixed = TRUE)[[1L]]),
              identical(ct$loaded, if (startsWith(selection, "hvg2000_bl")) "hvg2000_bl" else list()),
              "HR" %in% ct$metadata_columns, !"LR" %in% ct$metadata_columns,
              length(ct$published) == 0L)
    assert_error(execute(TRUE, c("--pseudobulk_variants", selection)), "pseudobulk_variants")
  }
  hr_keys <- result$COMPOSITION_HR_VARIANT_NAMES
  composition <- execute(FALSE, c("--composition_variants", paste(hr_keys, collapse = ",")), "composition")
  stopifnot(length(hr_keys) == 22L, identical(composition$scored, hr_keys),
            identical(composition$metadata_columns, c("Sample", "biology", "HR")),
            identical(composition$attached, c("reticulate", "dplyr", "GloScope")),
            length(composition$loaded) == 0L, length(composition$published) == 0L)
  plain_hr <- execute(FALSE, c("--composition_variants", "ECODA_authors_HR"), "composition")
  stopifnot(identical(plain_hr$attached, c("reticulate", "dplyr")))
  for (selection in c("ECODA_deconv", "ECODA_authors_HR,ECODA_authors_HR", "ECODA_authors_HR,")) {
    assert_error(execute(FALSE, c("--composition_variants", selection), "composition"), "composition_variants")
  }
  assert_error(execute(FALSE, c("--composition_variants", "ECODA_authors_HR")), "ordinary benchmark composition")
  assert_error(execute(FALSE, c("--composition_variants", "ECODA_authors_HR"), "composition", "uncorrected"), "ordinary benchmark composition")
  assert_error(execute(FALSE, c("--composition_variants=ECODA_authors_HR", "--composition_variants=GloProp"), "composition"), "Duplicate")
  full_prep <- execute(TRUE, selector = character())
  stopifnot(identical(full_prep$prepared, full_prep$PB_VARIANT_NAMES))
  for (preparation in c(TRUE, FALSE)) {
    for (selector in list(c("--pseudobulk_variants", "unknown"),
                          c("--pseudobulk_variants", "hvg2000_bl,hvg2000_bl"),
                          c("--pseudobulk_variants", ",hvg2000_bl"),
                          "--pseudobulk_variants")) {
      assert_error(execute(preparation, selector), "pseudobulk_variants")
    }
    assert_error(execute(preparation, c("--pseudobulk_variants=hvg2000_bl",
      "--pseudobulk_variants=hvg2000")), "Duplicate")
    assert_error(execute(preparation, pass = "uncorrected"), "ordinary benchmark")
    assert_error(execute(preparation, c("--pseudobulk-variants", "hvg2000_bl")),
                 "Unknown")
  }
  assert_error(execute(FALSE, method = "composition"), "ordinary benchmark")
})
cat("PASS: exact BL/HR-CT and 22 HR composition entrypoints, defaults and rejected contradictory selectors\n")
