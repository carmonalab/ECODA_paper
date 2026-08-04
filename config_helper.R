get_pipeline_config <- function(
  ds_name = Sys.getenv("DS_NAME"),
  force_overwrite = FALSE,
  test_mode = FALSE
) {
  if (ds_name == "") {
    stop("CRITICAL: DS_NAME not set. Ensure it is exported before calling R.")
  }

  hpc_scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR")
  if (hpc_scratch_dir == "") {
    stop("CRITICAL: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling R.")
  }

  scratch_output_dir <- Sys.getenv("SCRATCH_OUTPUT_DIR")
  if (scratch_output_dir == "") {
    scratch_output_dir <- file.path(hpc_scratch_dir, "output")
  }

  home_ref_dir <- Sys.getenv("HOME_REF_DIR")
  if (home_ref_dir == "") {
    stop("CRITICAL: HOME_REF_DIR not set. Source slurm_config.sh before calling R.")
  }

  config_data <- list(
    ds_name             = ds_name,
    test_mode           = test_mode,
    chunk_size          = ifelse(test_mode, 1, 5),
    max_test_array_jobs = 3,
    # All per-dataset dirs live under SCRATCH_OUTPUT_DIR/<DS_NAME> (preprocessed
    # output = annotation input). Matches 2_submit_hpc_array.sh HOME_CHUNKS_DIR.
    path_data           = file.path(scratch_output_dir, ds_name),
    path_output         = file.path(scratch_output_dir, ds_name),
    path_output_chunks  = file.path(scratch_output_dir, ds_name, "chunks"),
    path_ref            = home_ref_dir
  )

  return(config_data)
}
