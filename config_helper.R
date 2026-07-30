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
  gene_ref_file <- Sys.getenv("GENE_REF_FILE")

  config_data <- list(
    ds_name             = ds_name,
    test_mode           = test_mode,
    chunk_size          = ifelse(test_mode, 1, 5),
    max_test_array_jobs = 3,
    path_data           = file.path(hpc_scratch_dir, "data"),
    path_plots          = file.path(hpc_scratch_dir, "plots"),
    path_output         = scratch_output_dir,
    path_output_samples = file.path(scratch_output_dir, "samples"),
    path_output_chunks  = file.path(scratch_output_dir, "chunks"),
    path_output_ecoda   = file.path(scratch_output_dir, "ecoda"),
    path_ref            = if (home_ref_dir != "") home_ref_dir else file.path(Sys.getenv("HOME"), "reference_atlases", "sketched_200ct"),
    gene_ref            = if (gene_ref_file != "") gene_ref_file else file.path(getwd(), "EnsemblGenes105_Hsa_GRCh38.p13.txt.gz")
  )

  return(config_data)
}
