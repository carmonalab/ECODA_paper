get_pipeline_config <- function(
  ds_name = Sys.getenv("DS_NAME")
) {
  if (ds_name == "") {
    stop("CRITICAL: DS_NAME not set. Ensure DS_NAME is set before calling R.")
  }

  hpc_scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR")
  if (hpc_scratch_dir == "") {
    stop("CRITICAL: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling R.")
  }

  home_ref_dir <- Sys.getenv("HOME_REF_DIR")
  if (home_ref_dir == "") {
    stop("CRITICAL: HOME_REF_DIR not set. Source slurm_config.sh before calling R.")
  }

  scratch_output_dir <- file.path(hpc_scratch_dir, ds_name, "output")
  annotation_output_dir <- Sys.getenv("ANNOTATION_OUTPUT_DIR")
  if (annotation_output_dir == "") annotation_output_dir <- scratch_output_dir

  config_data <- list(
    ds_name             = ds_name,
    path_data           = scratch_output_dir,
    path_output         = annotation_output_dir,
    path_ref            = home_ref_dir
  )

  return(config_data)
}
