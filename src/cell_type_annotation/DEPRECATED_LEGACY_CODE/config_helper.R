readRenviron("config.env")

# --- ROBUST ENVIRONMENT & LIBRARY RESOLUTION ---
# If running inside an isolated cluster worker, renv's sandbox might mask packages.
# We explicitly ensure the project library and the Pixi system library are exposed.
project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") project_root <- getwd()
r_version_major_minor <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))

# Ensure renv paths are properly exposed on worker nodes
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
local_renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
pixi_lib <- file.path(project_root, ".pixi", "envs", "default", "lib", "R", "library")

.libPaths(unique(c(local_renv_lib, pixi_lib, .libPaths())))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("CRITICAL: 'jsonlite' could not be resolved. Check .libPaths(): ", paste(.libPaths(), collapse = ", "))
}
library(jsonlite)
# -----------------------------------------------

get_pipeline_config <- function(
  ds_name = Sys.getenv("DS_NAME"),
  force_overwrite = FALSE,
  test_mode = FALSE
) {
  config_file <- file.path(project_root, "pipeline_config.json")

  if (file.exists(config_file) && force_overwrite) {
    unlink(config_file)
  }

  if (file.exists(config_file) && !force_overwrite) {
    return(read_json(config_file, simplifyVector = TRUE))
  }


  nas_prefix <- Sys.getenv("NAS_PREFIX")

  path_home <- file.path(Sys.getenv("HOME"), ds_name)

  config_data <- list(
    ds_name             = ds_name,
    test_mode           = test_mode,
    chunk_size          = ifelse(test_mode, 1, 5),
    max_test_array_jobs = 3,
    path_nas            = file.path(nas_prefix, "DataCollections"),
    path_nas_root       = file.path(nas_prefix, "DataCollections", "Standardized_SingleCell_Datasets", ds_name),
    path_root           = path_home,
    path_data           = file.path(path_home, "data"),
    path_plots          = file.path(path_home, "plots"),
    path_output         = file.path(path_home, "output"),
    path_output_samples = file.path(path_home, "output", "samples"),
    path_output_chunks  = file.path(path_home, "output", "chunks"),
    path_output_ecoda   = file.path(path_home, "output", "ecoda"),
    path_ref            = file.path(path_home, "reference_atlases", "sketched_200ct"),
    gene_ref            = file.path(project_root, "EnsemblGenes105_Hsa_GRCh38.p13.txt.gz")
  )

  # Save it so it's available for subsequent scripts
  write_json(config_data, config_file, pretty = TRUE, auto_unbox = TRUE)
  message("Created new pipeline_config.json pointer file.")

  return(config_data)
}
