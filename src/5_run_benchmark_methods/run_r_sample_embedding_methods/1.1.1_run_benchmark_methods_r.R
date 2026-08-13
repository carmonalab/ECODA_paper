# ==============================================================================
# 1.1.1_run_benchmark_methods_r.R — Run one R benchmark method (gloscope, mofa,
# pseudobulk or scitd) for one dataset (Pipeline A).
#
# Called by 1.1_run_worker.sh via ${PIXI_RSCRIPT} with:
#   --config_path --ds_name --view benchmark_analysis --method {gloscope,mofa,
#   pseudobulk,scitd} --input_dir --results_dir --pseudobulk_dir
#   --gloscope_cache_dir --log_file [--force]
# Loads the preprocessed benchmark view h5ad -> Seurat (raw counts +
# X_pca_benchmark_analysis_hvg{n} obsm embeddings via reticulate), sets
# seurat@misc$cell_type_low_res / label_col from datasets.json, dispatches on
# --method to the T2 driver (benchmark_pipeline.R), writes the per-combo
# bundle files (<ds>_<combo>.rds; combo names are method-prefixed) + the
# method-level RDS (<ds>_<method>.rds, a named list of result bundles) +
# per-combo exec-log rows. Skip-if-exists: the method RDS exists -> re-emit
# its exec-log rows
# (failure-resume must not lose timing from an aborted run) and skip all
# unless --force; otherwise per-combo cache files are reused.
#
# Memory: mofa consumes only the precomputed pseudobulks, so the full Seurat
# object (multi-GB counts matrix) is built lazily only when a pb variant is
# missing on disk; gloscope fetches the embeddings, pseudobulk/scitd the
# counts (fetch_embedding = NULL).
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  stop("PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")
}

source(file.path(project_root, "src/utils/imports_worker_core.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))

raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_flags(raw_args)

for (req in c("config_path", "ds_name", "view", "method", "input_dir",
              "results_dir", "log_file")) {
  if (is.null(args[[req]]) || identical(args[[req]], TRUE)) {
    stop("Missing required --", req, " argument")
  }
}
force <- isTRUE(args[["force"]]) || identical(args[["force"]], "TRUE")

method <- args$method
if (!method %in% c("gloscope", "mofa", "pseudobulk", "scitd")) {
  stop("Unknown method '", method,
       "' (expected gloscope, mofa, pseudobulk or scitd)")
}

# Method-specific attaches: MOFA2/scITD are needed only by their methods
# (bare create_mofa / initialize_params + make_new_container); gloscope needs
# only the installed namespace (GloScope::gloscope is called qualified).
if (method == "mofa") library(MOFA2)
if (method == "scitd") library(scITD)

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
dir.create(args$results_dir, showWarnings = FALSE, recursive = TRUE)

method_rds <- file.path(args$results_dir, paste0(ds, "_", method, ".rds"))
if (file.exists(method_rds) && !force) {
  message("Method results already exist: ", method_rds,
          " (use --force to recompute)")
  cached <- readRDS(method_rds)
  for (nm in names(cached)) {
    if (!is.null(cached[[nm]]$exec_time)) {
      log_exec_row(ds, nm, cached[[nm]]$exec_time, args$log_file)
    }
  }
  quit(save = "no", status = 0)
}

ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(h5ad_path, backed = "r")
obs <- py_to_r(adata$obs)

sample_col <- "Sample"
if (!sample_col %in% colnames(obs)) {
  stop(sample_col, " not found in obs columns of ", h5ad_path)
}
hvg_rank_genes <- get_hvg_rank_genes(adata)

pb_variants <- NULL
seurat <- NULL
metadata <- NULL
labels <- NULL

if (method == "mofa") {
  # MOFA consumes only the precomputed pseudobulks: skip the multi-GB counts
  # materialization unless a pb variant is missing (on-the-fly fallback
  # needs the Seurat). Metadata/labels come straight from obs.
  if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
    stop("Missing required --pseudobulk_dir argument for method mofa")
  }
  dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
  if (length(pb_variants_missing(args$pseudobulk_dir, ds, force)) > 0) {
    message("Building Seurat for on-the-fly pseudobulk variants...")
    seurat <- load_benchmark_seurat(adata, obs, sample_col = sample_col,
                                    fetch_embedding = NULL)
    seurat@meta.data[[sample_col]] <- standardize_sample_names(
      seurat@meta.data[[sample_col]]
    )
  }
  pb_variants <- load_pb_variants(
    seurat, sample_col, hvg_rank_genes,
    pseudobulk_dir = args$pseudobulk_dir, ds = ds,
    force = force, log_file = args$log_file
  )
  obs[[sample_col]] <- standardize_sample_names(obs[[sample_col]])
  metadata <- obs %>%
    dplyr::group_by(!!sym(sample_col)) %>%
    dplyr::slice(1)
  labels <- as.factor(metadata[[entry$label_col]])
  names(labels) <- metadata[[sample_col]]
} else {
  seurat <- load_benchmark_seurat(
    adata, obs, sample_col = sample_col,
    fetch_embedding = if (method == "gloscope") {
      c("X_pca_benchmark_analysis_hvg1000",
        "X_pca_benchmark_analysis_hvg2000",
        "X_pca_benchmark_analysis_hvg3000")
    } else {
      NULL
    }
  )
  seurat@meta.data[[sample_col]] <- standardize_sample_names(
    seurat@meta.data[[sample_col]]
  )
  if (length(hvg_rank_genes) > 0) {
    VariableFeatures(seurat) <- hvg_rank_genes[
      seq_len(min(2000, length(hvg_rank_genes)))
    ]
  }
  seurat@misc$label_col <- entry$label_col
  seurat@misc$cell_type_low_res <- entry$cell_type_low_res
  seurat@misc$cell_type_high_res <- entry$cell_type_high_res

  metadata <- get_metadata(seurat)
  labels <- get_labels(seurat, entry$label_col)

  if (method == "pseudobulk") {
    if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
      stop("Missing required --pseudobulk_dir argument for method pseudobulk")
    }
    dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
    pb_variants <- load_pb_variants(
      seurat, sample_col, hvg_rank_genes,
      pseudobulk_dir = args$pseudobulk_dir, ds = ds,
      force = force, log_file = args$log_file
    )
  }
  if (method == "gloscope") {
    if (is.null(args$gloscope_cache_dir) ||
        identical(args$gloscope_cache_dir, TRUE)) {
      stop("Missing required --gloscope_cache_dir argument for method gloscope")
    }
    dir.create(args$gloscope_cache_dir, showWarnings = FALSE, recursive = TRUE)
  }
}

results <- switch(
  method,
  gloscope = run_gloscope_hpc(
    seurat, metadata, label_col = entry$label_col,
    sample_col = sample_col,
    gloscope_cache_dir = args$gloscope_cache_dir,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file
  ),
  mofa = run_mofa_hpc(
    metadata, labels, pb_variants,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file
  ),
  pseudobulk = run_pseudobulk_hpc(
    seurat, labels, pb_variants,
    sample_col = sample_col,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file
  ),
  scitd = run_scitd_hpc(
    seurat, label_col = entry$label_col,
    hvg_sets = make_hvg_sets(hvg_rank_genes),
    sample_col = sample_col,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file
  )
)

save_rds_atomic(results, method_rds)
message("Saved method results: ", method_rds, " (", length(results), " combos)")
message("--- ", method, " for ", ds, " complete ---")
