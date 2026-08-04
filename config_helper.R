get_pipeline_config <- function(
  ds_name = Sys.getenv("DS_NAME"),
  force_overwrite = FALSE
) {
  if (ds_name == "") {
    stop("CRITICAL: DS_NAME not set. Ensure it is exported before calling R.")
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

  config_data <- list(
    ds_name             = ds_name,
    # All per-dataset dirs live under HPC_SCRATCH_DIR/<DS_NAME>/output
    # (preprocessed output = annotation input). Matches 2_submit_hpc_array.sh
    # CHUNKS_DIR.
    path_data           = scratch_output_dir,
    path_output         = scratch_output_dir,
    path_output_chunks  = file.path(scratch_output_dir, "chunks"),
    path_ref            = home_ref_dir
  )

  return(config_data)
}
