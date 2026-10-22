#!/usr/bin/env Rscript
# Real corrected-final preparation and pseudobulk workers; synthetic data only.
# Run: pixi run -e default Rscript --vanilla tests/test_corrected_worker_recovery.R
# Requires repository .pixi/envs/default Python/Rscript, Python anndata/numpy/
# pandas/pyarrow/scipy/h5py, and the real worker R dependencies listed below.

script_arg <- commandArgs(trailingOnly = FALSE)[grepl(
  "^--file=", commandArgs(trailingOnly = FALSE)
)][1]
stopifnot(length(script_arg) == 1L, !is.na(script_arg))
root <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg)), ".."),
  mustWork = TRUE
)
python <- file.path(root, ".pixi", "envs", "default", "bin", "python")
rscript <- file.path(root, ".pixi", "envs", "default", "bin", "Rscript")
for (executable in c(python, rscript)) {
  if (file.access(executable, mode = 1L) != 0L) {
    stop("Required repository Pixi executable is unavailable: ", executable)
  }
}
required_packages <- c(
  "jsonlite", "Seurat", "reticulate", "dplyr", "arrow", "DESeq2", "limma",
  "Matrix", "MatrixGenerics", "SummarizedExperiment", "vegan", "cluster",
  "igraph", "mclust", "thisutils"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Required real-worker R dependencies are missing: ",
       paste(missing_packages, collapse = ", "))
}

with_temporary <- function(code) {
  directory <- tempfile("ecoda-corrected-worker-")
  dir.create(directory, recursive = TRUE)
  directory <- normalizePath(directory, mustWork = TRUE)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  previous_wd <- getwd()
  setwd(root) # load_worker_functions.R intentionally uses project-relative paths.
  on.exit(setwd(previous_wd), add = TRUE)
  eval(substitute(code), envir = environment())
}

with_temporary({
  input_root <- file.path(directory, "input")
  analysis_root <- file.path(directory, "batch_effect", "corrected_final")
  pseudobulk_dir <- file.path(analysis_root, "pseudobulks")
  results_dir <- file.path(analysis_root, "results")
  log_root <- file.path(directory, "logs")
  runs_root <- file.path(directory, "runs")
  for (path in c(input_root, pseudobulk_dir, results_dir, log_root,
                 runs_root, file.path(analysis_root, "metadata"))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  h5ad_path <- file.path(input_root, "synthetic-corrected.h5ad")
  metadata_path <- file.path(
    analysis_root, "metadata", "Synthetic_sample_metadata.feather"
  )
  config_path <- file.path(directory, "datasets.json")
  jsonlite::write_json(
    list(Synthetic = list(
      use_for_batch_effect = TRUE,
      columns = list(
        sample = "Sample", label = "label", batch = list("assay", "sex"),
        cell_type_low_res = "author_cell_type",
        cell_type_high_res = "author_cell_type"
      ),
      views = list(batch_effect_corrected = list(
        input_file_name = "synthetic.h5ad",
        output_file_name = basename(h5ad_path), subset_vars = list()
      ))
    )),
    config_path, auto_unbox = TRUE, pretty = TRUE
  )
  # Nonlexical encounter order catches accidental sample sorting in either worker.
  sample_ids <- c("s4", "s1", "s6", "s2", "s5", "s3")
  fixture_python <- paste(c(
    "import hashlib, sys",
    "from pathlib import Path",
    "import anndata as ad, numpy as np, pandas as pd, pyarrow.feather as feather",
    "from scipy import sparse",
    "import h5py",
    "sys.path.insert(0, sys.argv[2])",
    "from src.utils.py.batch_contract import build_batch_contract_identity",
    "rng = np.random.default_rng(20261021)",
    "samples = ['s4', 's1', 's6', 's2', 's5', 's3']",
    "metadata = pd.DataFrame({'Sample': samples, 'label': ['A']*3 + ['B']*3,",
    "    'assay': ['a','b','a','b','a','b'], 'sex': ['F','F','M','M','F','M']})",
    "design = np.column_stack([np.ones(6), (metadata.assay == 'b').astype(float), (metadata.sex == 'M').astype(float)])",
    "assert np.linalg.matrix_rank(design) == 3",
    "obs = metadata.loc[metadata.index.repeat(12)].reset_index(drop=True)",
    "obs['author_cell_type'] = np.tile(['T']*6 + ['B']*6, 6)",
    "obs.index = [f'cell{i}' for i in range(len(obs))]",
    "means = rng.lognormal(2.5, 0.8, size=(1, 2000))",
    "sample_effects = np.repeat(rng.lognormal(0, 0.35, size=(6, 2000)), 12, axis=0)",
    "counts = rng.negative_binomial(5, 5 / (5 + means * sample_effects)).astype(np.int64) + 1",
    "assert counts.shape == (72, 2000) and np.all(counts > 0)",
    "adata = ad.AnnData(X=sparse.csr_matrix(np.log1p(counts).astype(np.float32)), obs=obs,",
    "    var=pd.DataFrame({'hvg_rank': np.arange(1, 2001, dtype=float)}, index=[f'g{i}' for i in range(2000)]))",
    "adata.layers['counts'] = sparse.csr_matrix(counts)",
    "adata.obsm['X_pca_batch_effect_corrected_hvg2000'] = rng.normal(size=(72, 5)).astype(np.float32)",
    "adata.obsm['X_pca_harmony_batch_effect_corrected_hvg2000'] = rng.normal(size=(72, 5)).astype(np.float32)",
    "adata.uns['batch_contract'] = build_batch_contract_identity(['assay','sex'], sample_column='Sample', method_id='preprocess', model_id='hvg_composite_v1')",
    "adata.write_h5ad(sys.argv[1])",
    "path = Path(sys.argv[3])",
    "feather.write_feather(metadata, path)",
    "digest = hashlib.md5(path.read_bytes()).hexdigest()",
    "Path(str(path) + '.md5').write_text(f'MD5={digest}\\nSIZE={path.stat().st_size}\\nPATH={path}\\n')"
  ), collapse = "\n")
  environment_values <- c(
    PROJECT_ROOT = root, ECODA_SOURCE_ROOT = root,
    RETICULATE_PYTHON = python, PIXI_PYTHON = python,
    ANALYSIS_ROOT = analysis_root, ANALYSIS_VARIANT = "corrected_final",
    ANALYSIS_PASS = "corrected", ECODA_RUN_ID = "corrected-worker-test",
    ECODA_RUNS_ROOT = runs_root,
    ECODA_RUN_ROOT = file.path(runs_root, "corrected-worker-test"),
    ECODA_SOURCE_IDENTITY = "", ECODA_ARTIFACT_PRODUCER = "",
    METHOD = "pseudobulk", ANALYSIS = "batch_effect_corrected",
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
  )
  worker_env <- paste0(names(environment_values), "=", shQuote(environment_values))
  invoke <- function(executable, arguments, label) {
    output <- suppressWarnings(system2(
      executable, shQuote(arguments), stdout = TRUE, stderr = TRUE,
      env = worker_env
    ))
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
      stop(label, " failed (exit ", status, "):\n",
           paste(output, collapse = "\n"), call. = FALSE)
    }
    invisible(output)
  }
  invoke(
    python, c("-c", fixture_python, h5ad_path, root, metadata_path),
    "Synthetic corrected H5AD/Feather creation"
  )
  worker_dir <- file.path(
    root, "src", "5_run_benchmark_methods", "run_r_sample_embedding_methods"
  )
  common_args <- c(
    "--config_path", config_path, "--ds_name", "Synthetic",
    "--view", "batch_effect_corrected", "--analysis_pass", "corrected",
    "--input_dir", input_root, "--pseudobulk_dir", pseudobulk_dir
  )
  prepare <- function() invoke(
    rscript,
    c("--vanilla", file.path(worker_dir, "1.1.1_prepare_pseudobulk.R"),
      common_args, "--log_file", file.path(log_root, "prepare.csv")),
    "Actual corrected preparation worker"
  )
  benchmark <- function() invoke(
    rscript,
    c("--vanilla", file.path(worker_dir, "1.1.1_run_benchmark_methods_r.R"),
      common_args, "--method", "pseudobulk", "--results_dir", results_dir,
      "--gloscope_cache_dir", file.path(analysis_root, "gloscope_cache"),
      "--log_file", file.path(log_root, "pseudobulk.csv")),
    "Actual corrected pseudobulk benchmark worker"
  )
  # No force flag: this must fail at the unavailable callee before its repair.
  prepare()
  cache_path <- file.path(
    pseudobulk_dir,
    "Synthetic_batch_effect_corrected_final_pseudobulk_hvg2000.rds"
  )
  stopifnot(file.exists(cache_path), file.exists(paste0(cache_path, ".md5")))
  cache_md5 <- unname(tools::md5sum(cache_path))
  cache <- readRDS(cache_path)
  assert_contract <- function(contract) {
    stopifnot(
      identical(unname(unlist(contract$ordered_source_keys)), c("assay", "sex")),
      identical(unname(unlist(contract$required_source_obs_columns)),
                c("Sample", "assay", "sex")),
      identical(contract$method_id, "Pseudobulk"),
      identical(contract$model_id, "pseudobulk_limma_fixed_effects_v1"),
      identical(unname(unlist(contract$effective_batch_keys)), c("assay", "sex")),
      length(unlist(contract$non_estimable_batch_keys)) == 0L,
      identical(contract$correction_state, "BATCH_CORRECTION"),
      identical(contract$correction_mode, "limma_fixed_effects_pseudobulk"),
      identical(contract$correction_formula, paste0(
        "DESeq2 design=~ 1; model.matrix(~ 1 + batch_key_1 + batch_key_2); ",
        "limma::removeBatchEffect(covariates=technical_covariates, ",
        "design=intercept); effective_batch_keys=[assay,sex]"
      )),
      identical(unname(unlist(contract$fixed_effect_aliases)),
                c("batch_key_1", "batch_key_2")),
      contract$design_rank == 3L,
      contract$design_residual_df == 3L
    )
  }
  stopifnot(
    is.matrix(cache$pb), identical(dim(cache$pb), c(6L, 2000L)),
    identical(rownames(cache$pb), sample_ids), all(is.finite(cache$pb))
  )
  assert_contract(cache$batch_contract)
  benchmark()
  result_path <- file.path(
    results_dir, "Synthetic_batch_effect_corrected_final_pseudobulk.rds"
  )
  stopifnot(file.exists(result_path), file.exists(paste0(result_path, ".md5")))
  result <- readRDS(result_path)
  bundle <- result$Pseudobulk_hvg2000
  stopifnot(
    identical(rownames(bundle$feat_mat), sample_ids),
    identical(bundle$feat_mat, cache$pb), all(is.finite(bundle$feat_mat)),
    identical(rownames(as.matrix(bundle$dist_mat)), sample_ids),
    identical(colnames(as.matrix(bundle$dist_mat)), sample_ids),
    all(is.finite(as.matrix(bundle$dist_mat))),
    identical(names(bundle$labels), sample_ids),
    identical(as.character(bundle$labels), c("A", "A", "A", "B", "B", "B"))
  )
  assert_contract(bundle$batch_contract)
  assert_contract(result$batch_contract)
  prepare()
  benchmark()
  stopifnot(
    identical(unname(tools::md5sum(cache_path)), cache_md5),
    identical(readRDS(cache_path), cache)
  )
})

cat("corrected-final real-worker recovery: OK\n")
