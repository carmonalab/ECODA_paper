# ==============================================================================
# 1.1.1_run_zeroimp_analysis.R — ECODA zero-imputation analysis for one
# dataset (Pipeline B, zeroimp array).
#
# Called by 1.1_run_worker.sh via ${PIXI_RSCRIPT} with:
#   --config_path --ds_name --view benchmark_analysis --input_dir
#   --output_dir --log_file [--force]
# obs-only backed read (reticulate, no counts matrix); builds the per-sample
# cell-type composition table + labels from datasets.json's label_col, runs
# run_zeroimp_analysis(), saves results/<ds>_zeroimp.rds atomically with one
# exec-log row per dataset (method "zeroimp_analysis"). Skip-if-exists unless
# --force. The shared driver lives in benchmark_hpc_utils.R
# (run_ct_comps_analysis_worker).
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  stop("PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")
}

source(file.path(project_root, "src/utils/imports_worker_transzeroimp.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))

run_ct_comps_analysis_worker(
  analysis_label = "zeroimp",
  run_fun = run_zeroimp_analysis,
  out_suffix = "_zeroimp",
  log_method = "zeroimp_analysis"
)
