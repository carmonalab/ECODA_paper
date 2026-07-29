project_root <- getwd() # Since the bash script changes directory to PROJECT_ROOT

# Load central config into R's environment
readRenviron(file.path(project_root, "config.env"))

r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

raw_args <- commandArgs(trailingOnly = TRUE)

# Set up parameter parser defaults
defaults <- list(test = "false")
args <- defaults

if (length(raw_args) > 0) {
  parsed_args_list <- unlist(strsplit(raw_args[1], "__"))
  keys <- parsed_args_list[seq(1, length(parsed_args_list), by = 2)]
  vals <- parsed_args_list[seq(2, length(parsed_args_list), by = 2)]

  args_list <- as.list(vals)
  names(args_list) <- keys
  args <- modifyList(defaults, args_list)
}

# Convert string argument to an R logical boolean flag
RUN_AS_TEST <- as.logical(args$test)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ###### Set paths ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Fetch dataset name dynamically from config.env
ds_name <- Sys.getenv("DS_NAME")
if (ds_name == "") stop("CRITICAL Error: DS_NAME not found in config.env")

source("config_helper.R")
paths <- get_pipeline_config(ds_name, force_overwrite = TRUE, test_mode = RUN_AS_TEST)

message(paste("Path is:", paths$path_data))
message(paste("Files found:", paste(list.files(paths$path_data), collapse = ", ")))

dir.create(paths$path_output, showWarnings = FALSE)
dir.create(paths$path_output_samples, showWarnings = FALSE)
# Delete chunk file folder recursively to ensure a perfectly clean start
if (dir.exists(paths$path_output_chunks)) {
  unlink(paths$path_output_chunks, recursive = TRUE, force = TRUE)
}
dir.create(paths$path_output_chunks, showWarnings = FALSE)
dir.create(paths$path_output_ecoda, showWarnings = FALSE)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ###### Create chunk.txt files ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# They contain the information to run per cpu core (file and samples)

# Import anndata WITHOUT automatic R conversion
library(reticulate)

# # Use pixi python
# pixi_python <- file.path(getwd(), ".pixi", "envs", "default", "bin", "python")
# reticulate::use_python(pixi_python, required = TRUE)
py_require("anndata")
ad <- import("anndata", convert = FALSE)

h5ad_files <- list.files(paths$path_data, pattern = "\\.h5ad$", full.names = TRUE)

# Number of samples per chunk
chunk_size <- if (RUN_AS_TEST) 1 else 5
global_chunk_counter <- 1

sample_col <- Sys.getenv("SAMPLE_COLNAME")

# Loop through each file individually
for (f in h5ad_files) {
  message(paste("Processing file-specific chunks for:", basename(f)))

  # 1. Read the h5ad file in backed mode and extract its unique samples
  adata <- ad$read_h5ad(f, backed = "r")
  # 2. Convert ONLY the obs data frame into a native R data.frame
  r_obs <- py_to_r(adata$obs)
  
  # 3. Use standard R matrix/data.frame extraction rules cleanly
  file_samples <- as.character(unique(r_obs[[sample_col]]))

  # 2. Split ONLY this file's samples into groups of 10
  # ceiling(seq_along(...) / 10) creates grouping factors specific to this pool
  sample_groups <- split(file_samples, ceiling(seq_along(file_samples) / chunk_size))

  # 3. Write each group to a unique chunk file
  for (i in seq_along(sample_groups)) {
    # We write the source file as the VERY FIRST line, followed by the sample IDs
    writeLines(
      text = c(f, sample_groups[[i]]),
      con = file.path(
        paths$path_output_chunks,
        sprintf("chunk_%d.txt", global_chunk_counter)
      )
    )

    global_chunk_counter <- global_chunk_counter + 1
  }
}

message(sprintf("Successfully generated %d total chunk files.", global_chunk_counter - 1))
