# Exercise actual selected orchestration with only score/store I/O stubbed.
root <- normalizePath(".")
env <- new.env(parent = globalenv())
sys.source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"), envir = env)
sys.source(file.path(root, "src/5_run_benchmark_methods/benchmark_pipeline.R"), envir = env)
assert_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
env$artifact_checksum_ok <- function(...) FALSE
env$read_rds_checked <- function(...) stop("unselected cached result read")
env$exec_time <- function(expr) { force(expr); as.difftime(1, units = "secs") }
env$peak_rss_gb <- function() 0
env$published <- list()
env$logged <- character()
env$save_rds_atomic <- function(object, file, ...) env$published[[basename(file)]] <- object
env$log_exec_row <- function(ds, method, ...) env$logged <- c(env$logged, method)
ids <- paste0("s", 1:6)
metadata <- data.frame(Sample = ids, biology = rep(c("A", "B"), 3))
obs <- metadata[rep(1:6, each = 3), ]
obs$HR <- rep(c("T", "B", "M"), 6)
labels <- setNames(factor(metadata$biology), ids)
bundle <- list(labels = labels, feat_mat = matrix(1, 6, 2, dimnames = list(ids, c("a", "b"))))
env$process_coda_fig <- function(seurat, labels, ..., ct_col, obs) {
  stopifnot(is.null(seurat), identical(ct_col, "HR"), "HR" %in% names(obs))
  bundle
}
env$process_gloprop_fig <- function(seurat, metadata, ct_col, label_col, obs, ...) {
  stopifnot(is.null(seurat), identical(ct_col, "HR"), identical(label_col, "biology"))
  bundle
}
env$process_avg_pca_embedding_fig <- env$process_deconv_fig <-
  function(...) stop("unselected composition method reached")
keys <- env$COMPOSITION_HR_VARIANT_NAMES
result <- env$run_composition_methods_hpc(
  labels, metadata, NULL, NULL, obs, "biology", ct_col_high_res = "HR",
  results_dir = tempdir(), ds = "Toy", variants = keys
)
stopifnot(length(keys) == 22L, identical(names(result), keys),
          setequal(names(env$published), paste0("Toy_", c("metadata", keys), ".rds")),
          !"Toy_composition.rds" %in% names(env$published),
          identical(env$published$Toy_metadata.rds$n_cells, nrow(obs)))
assert_error(env$run_composition_methods_hpc(
  labels, metadata, NULL, NULL, obs, "biology", results_dir = tempdir(), ds = "Toy",
  variants = "ECODA_deconv"), "ordinary author-HR")
assert_error(env$run_composition_methods_hpc(
  labels, metadata, NULL, NULL, obs, "biology", results_dir = tempdir(), ds = "Toy",
  variants = "ECODA_authors_HR"), "high-resolution")

# Both selected HVGs use one HR store/timing identity; no LR, plain or PCA fit.
env$published <- list()
env$logged <- character()
store_calls <- 0L
env$process_pseudobulk_fig <- function(...) stop("ordinary/PCA fit reached")
env$process_pseudobulk_ct_h5ad_variants_fig <- function(
  h5ad_path, labels, sample_col, ct_col, hvgs, ..., timing_id
) {
  store_calls <<- store_calls + 1L
  stopifnot(identical(ct_col, "HR"), identical(hvgs, c(2000L, 500L)))
  value <- c(bundle, list(exec_time = 1, variant_time_secs = 1,
                          shared_time_secs = 2, shared_mem_GB = 0,
                          timing_id = timing_id, timing_schema = 2L))
  setNames(rep(list(value), length(hvgs)), paste0("hvg", hvgs))
}
env$process_pseudobulk_ct_h5ad_fig <- function(..., timing_id) stop("duplicate CT store reached")
ct_keys <- env$PB_CT_HR_VARIANT_NAMES
ct <- env$run_pseudobulk_hpc(
  labels, list(), results_dir = tempdir(), ds = "Toy", h5ad_path = "stub.h5ad",
  ct_col_low_res = "forbidden_LR", ct_col_high_res = "HR",
  variants = ct_keys, run_id = "selected-ct"
)
stopifnot(store_calls == 1L, identical(names(ct), paste0("Pseudobulk_", ct_keys)),
          setequal(names(env$published), paste0("Toy_Pseudobulk_", ct_keys, ".rds")),
          sum(grepl("^prepare_pseudobulk_ct_shared_", env$logged)) == 1L,
          identical(ct[[1L]]$timing_id, ct[[2L]]$timing_id))
assert_error(env$run_pseudobulk_hpc(
  labels, list(), results_dir = tempdir(), ds = "Toy", variants = ct_keys), "high-resolution")
assert_error(env$run_pseudobulk_hpc(
  labels, list(), results_dir = tempdir(), ds = "Toy", variants = "CT_LR_hvg2000"), "ordinary variant")

# Real obs reader accepts a preprocessed H5AD without embeddings or Leiden.
# The existing path-contract boundary is isolated here; its own tests cover
# producer identities. This check exercises the metadata-only read itself.
local({
  prior <- Sys.getenv(c("PROJECT_ROOT", "RETICULATE_PYTHON", "RETICULATE_USE_MANAGED_VENV"), unset = NA_character_)
  on.exit(for (name in names(prior)) {
    if (is.na(prior[[name]])) Sys.unsetenv(name) else
      do.call(Sys.setenv, setNames(list(prior[[name]]), name))
  })
  Sys.setenv(PROJECT_ROOT = root, RETICULATE_PYTHON = file.path(root, ".pixi/envs/default/bin/python"),
             RETICULATE_USE_MANAGED_VENV = "no")
  fixture <- tempfile(fileext = ".h5ad")
  on.exit(unlink(fixture), add = TRUE)
  reticulate::py_run_string(paste0(
    "import anndata as ad, pandas as pd, scipy.sparse as sp\n",
    "data=ad.AnnData(X=sp.csr_matrix((4,2)),obs=pd.DataFrame({'Sample':['s1','s1','s2','s2'],",
    "'biology':['A','A','B','B'],'HR':['T','B','T','M']},index=['c1','c2','c3','c4']))\n",
    "data.layers['counts']=data.X.copy()\ndata.var['hvg_rank']=[0.,1.]\ndata.write_h5ad(",
    deparse(fixture), ")"
  ))
  env$ecoda_hpc_validate_h5ad_path_identity <- function(...) invisible(TRUE)
  actual <- env$load_h5ad_author_hr_metadata(fixture, c("Sample", "biology", "HR"), "benchmark_analysis")
  stopifnot(identical(actual$Sample, c("s1", "s1", "s2", "s2")),
            identical(actual$HR, c("T", "B", "T", "M")))
})
cat("PASS: selected 22 HR composition shards and one-store HR CT pair without unrelated work\n")
