#!/usr/bin/env Rscript
# Scoped 16-row classic full-H5AD profiling worker for R benchmark methods.
# Reproduces the pre-streaming (fdf137d^) H5AD loading path in memory and
# records wall-clock runtime (time_secs) and Linux VmHWM peak RSS (mem_GB).

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  stop("PROJECT_ROOT not set. Source src/slurm_config.sh before running.")
}

source(file.path(project_root, "src/utils/imports_worker_core.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))

raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_flags(raw_args)
for (req in c("ds_name", "method", "out_csv", "work_dir")) {
  if (is.null(args[[req]]) || identical(args[[req]], TRUE) || !nzchar(args[[req]])) {
    stop("Missing required --", req, " argument")
  }
}

ds <- args$ds_name
method <- args$method
out_csv <- args$out_csv
work_dir <- args$work_dir
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

allowed_methods <- c(
  "ECODA_authors_HR",
  "GloProp",
  "Pseudobulk_hvg2000",
  "MOFA_hvg2000_factors15",
  "GloScope_hvg2000_pcadims30"
)
if (!method %in% allowed_methods) {
  stop("Unsupported R profiling method: ", method)
}

if (method == "GloProp") library(GloScope)
if (method == "MOFA_hvg2000_factors15") library(MOFA2)

config_path <- file.path(project_root, "datasets.json")
scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR", unset = "")
if (!nzchar(scratch_dir)) {
  stop("HPC_SCRATCH_DIR not set.")
}
input_dir <- file.path(scratch_dir, ds, "output")
config <- read_datasets_json(config_path, view = "benchmark_analysis")
entry <- config[[ds]]
if (is.null(entry)) {
  stop("Dataset not found in datasets.json: ", ds)
}
h5ad_path <- get_h5ad_path(config, ds, "benchmark_analysis", input_dir)
if (!file.exists(h5ad_path)) {
  stop("Input h5ad not found: ", h5ad_path)
}

driver_tag <- switch(
  method,
  ECODA_authors_HR = "composition",
  GloProp = "composition",
  Pseudobulk_hvg2000 = "pseudobulk",
  MOFA_hvg2000_factors15 = "mofa",
  GloScope_hvg2000_pcadims30 = "gloscope"
)

sample_col <- "Sample"
load_start <- proc.time()[["elapsed"]]
ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(h5ad_path, backed = "r")
obs <- py_to_r(adata$obs)
validate_benchmark_h5ad_contract(
  adata,
  obs = obs,
  view = "benchmark_analysis",
  method = driver_tag
)
hvg_rank_genes <- get_hvg_rank_genes(adata)

if (driver_tag == "composition") {
  emb_key <- "X_pca_benchmark_analysis_hvg2000"
  pca_emb <- py_to_r(adata$obsm[[emb_key]])
  if (is.null(rownames(pca_emb))) rownames(pca_emb) <- rownames(obs)
  colnames(pca_emb) <- paste0("PC_", seq_len(ncol(pca_emb)))
  obs <- rename_leiden_cols(obs, view = "benchmark_analysis", preserve_source = TRUE)
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
  labels <- as.factor(metadata[[entry$label_col]])
  names(labels) <- metadata[[sample_col]]
  seurat <- NULL
} else {
  fetch_emb <- if (method == "GloScope_hvg2000_pcadims30") {
    c(
      "X_pca_benchmark_analysis_hvg1000",
      "X_pca_benchmark_analysis_hvg2000",
      "X_pca_benchmark_analysis_hvg3000"
    )
  } else {
    NULL
  }
  seurat <- load_benchmark_seurat(
    adata,
    obs,
    sample_col = sample_col,
    fetch_embedding = fetch_emb,
    counts_layer = "counts"
  )
  if (length(hvg_rank_genes) > 0) {
    VariableFeatures(seurat) <- hvg_rank_genes[seq_len(min(2000, length(hvg_rank_genes)))]
  }
  seurat@misc$label_col <- entry$label_col
  seurat@misc$cell_type_low_res <- entry$cell_type_low_res
  seurat@misc$cell_type_high_res <- entry$cell_type_high_res
  metadata <- get_metadata(seurat)
  labels <- get_labels(seurat, entry$label_col)
}
load_secs <- as.numeric(proc.time()[["elapsed"]] - load_start)

set.seed(123)
if (method == "ECODA_authors_HR") {
  t_run <- exec_time(
    res <- process_coda_fig(NULL, labels, ct_col = entry$cell_type_high_res, obs = obs)
  )
  time_secs <- as.numeric(t_run, units = "secs")
} else if (method == "GloProp") {
  t_run <- exec_time(
    res <- process_gloprop_fig(
      NULL, metadata,
      ct_col = entry$cell_type_high_res,
      label_col = entry$label_col,
      obs = obs
    )
  )
  time_secs <- as.numeric(t_run, units = "secs")
} else if (method == "Pseudobulk_hvg2000") {
  t_pb <- exec_time(
    pb <- get_pb_deseq2(seurat, sample_col = sample_col, n_hvg = 2000, black_list = "none")
  )
  t_proc <- exec_time(
    res <- process_pseudobulk_fig(pb, labels)
  )
  time_secs <- as.numeric(t_pb, units = "secs") + as.numeric(t_proc, units = "secs")
} else if (method == "MOFA_hvg2000_factors15") {
  t_pb <- exec_time(
    pb <- get_pb_deseq2(seurat, sample_col = sample_col, n_hvg = 2000, black_list = "none")
  )
  t_mofa <- exec_time(
    res <- process_mofa_bulk_fig(pb, metadata = metadata, labels, num_factors = 15)
  )
  time_secs <- as.numeric(t_pb, units = "secs") + as.numeric(t_mofa, units = "secs")
} else if (method == "GloScope_hvg2000_pcadims30") {
  dist_file <- file.path(work_dir, paste0(ds, "_gloscope_hvg2000_pcadims30_dists.rds"))
  t_gs <- exec_time(
    res <- process_gloscope_fig(
      embedding_matrix = seurat@reductions[["pca_benchmark_analysis_hvg2000"]]@cell.embeddings,
      sample_ids = seurat@meta.data[[sample_col]],
      metadata = metadata,
      label_col = entry$label_col,
      gloscope_dist_file = dist_file,
      n_pca_dims = 30,
      force = TRUE
    )
  )
  time_secs <- as.numeric(t_gs, units = "secs")
}

mem_gb <- peak_rss_gb()
total_wall_secs <- load_secs + time_secs

row_df <- data.frame(
  dataset = ds,
  method = method,
  lang = "R",
  resource_mode = "cpu",
  time_secs = time_secs,
  load_secs = load_secs,
  total_wall_secs = total_wall_secs,
  mem_GB = mem_gb,
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = ""),
  slurm_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
  node = Sys.getenv("ECODA_PROF_NODE", unset = Sys.info()[["nodename"]]),
  cpu_model = Sys.getenv("ECODA_PROF_CPU_MODEL", unset = ""),
  gpu_model = Sys.getenv("ECODA_PROF_GPU_MODEL", unset = "none"),
  stringsAsFactors = FALSE
)

tmp_csv <- paste0(out_csv, ".tmp.", Sys.getpid())
utils::write.csv(row_df, tmp_csv, row.names = FALSE, quote = TRUE)
if (!file.rename(tmp_csv, out_csv)) {
  stop("Failed to atomically write CSV: ", out_csv)
}
message(sprintf(
  "PROFILE_R_DONE dataset=%s method=%s time_secs=%.3f load_secs=%.3f mem_GB=%.4f",
  ds, method, time_secs, load_secs, mem_gb
))
