# Extract the compact sample-level inputs used by the scECODA example.
# H5AD access is restricted to obs; expression matrices and count layers are
# never opened.

.ecoda_example_source_file <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[[1L]])
    if (!script_path %in% c("-e", "-")) {
      return(normalizePath(script_path, mustWork = TRUE))
    }
  }
  source_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1))
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files)) normalizePath(tail(source_files, 1L), mustWork = TRUE) else ""
})

.ecoda_example_script_path <- function() {
  if (nzchar(.ecoda_example_source_file)) return(.ecoda_example_source_file)
  source_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1))
  source_files <- source_files[nzchar(source_files)]
  normalizePath(tail(source_files, 1L), mustWork = TRUE)
}

.ecoda_example_project_root <- function() {
  dirname(dirname(dirname(.ecoda_example_script_path())))
}

.ecoda_example_python_vector <- function(value) {
  value <- reticulate::py_to_r(value)
  if (is.list(value)) {
    value <- lapply(value, function(item) {
      if (is.null(item) || !length(item)) NA else item
    })
    value <- unlist(value, use.names = FALSE)
  }
  as.vector(value)
}

.ecoda_example_is_missing <- function(value) {
  text <- tolower(trimws(as.character(value)))
  is.na(value) | text %in% c("", "na", "nan", "none", "<na>", "null")
}

.ecoda_example_drop_obs_column <- function(column, source_cell_type_columns) {
  column %in% c(
    source_cell_type_columns,
    "seurat_clusters", "S.Score", "G2M.Score", "Phase",
    "classification_confidence"
  ) ||
    grepl(
      "^(RNA_snn_res\\.|leiden_res_|scGate|functional\\.cluster|scATOMIC|cellCycle\\.|layer_?[1-6]$)|_UCell$",
      column
    )
}

.ecoda_example_collapse_column <- function(values, sample_rows) {
  collapsed <- lapply(sample_rows, function(rows) {
    sample_values <- values[rows]
    sample_values <- sample_values[!.ecoda_example_is_missing(sample_values)]
    distinct_values <- unique(sample_values)
    if (length(distinct_values) > 1L) return(NULL)
    if (length(distinct_values)) distinct_values[[1L]] else NA
  })
  if (any(vapply(collapsed, is.null, logical(1)))) return(NULL)
  do.call(c, collapsed)
}

.ecoda_example_read_h5ad_obs <- function(
  path,
  sample_column,
  source_cell_type_columns,
  include_hitme_lr,
  decoder,
  h5py,
  builtins,
  get_ct_comp_df
) {
  handle <- h5py$File(path, "r")
  on.exit(handle$close(), add = TRUE)
  obs <- handle[["obs"]]
  index_name <- as.character(reticulate::py_to_r(obs$attrs$get("_index")))
  obs_columns <- as.character(reticulate::py_to_r(builtins$list(obs$keys())))
  sample_values <- as.character(.ecoda_example_python_vector(
    decoder$read_obs_column_values(obs, sample_column)
  ))
  sample_ids <- unique(sample_values[!.ecoda_example_is_missing(sample_values)])
  sample_rows <- split(
    seq_along(sample_values),
    factor(match(sample_values, sample_ids), levels = seq_along(sample_ids))
  )

  metadata <- list()
  layer1_counts <- NULL
  for (column in setdiff(obs_columns, index_name)) {
    if (include_hitme_lr && identical(column, "layer1")) {
      layer1_values <- .ecoda_example_python_vector(
        decoder$read_obs_column_values(obs, column)
      )
      layer1_obs <- data.frame(
        sample_values, layer1_values,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      names(layer1_obs) <- c(sample_column, "layer1")
      layer1_counts <- get_ct_comp_df(
        layer1_obs, sample_col = sample_column, ct_col = "layer1"
      )
    }
    if (.ecoda_example_drop_obs_column(column, source_cell_type_columns)) next
    values <- if (identical(column, sample_column)) {
      sample_values
    } else {
      .ecoda_example_python_vector(decoder$read_obs_column_values(obs, column))
    }
    collapsed <- .ecoda_example_collapse_column(values, sample_rows)
    if (!is.null(collapsed)) metadata[[column]] <- collapsed
  }

  metadata <- as.data.frame(
    metadata, check.names = FALSE, stringsAsFactors = FALSE
  )
  rownames(metadata) <- sample_ids
  list(metadata = metadata, layer1_counts = layer1_counts)
}

.ecoda_example_result_counts <- function(results_dir, dataset, method) {
  bundle <- readRDS(file.path(results_dir, paste0(dataset, "_", method, ".rds")))
  t(as.matrix(bundle$counts) - 0.5)
}

.ecoda_example_gongsharma_counts <- function(data, layer, sample_ids) {
  sample_column <- "specimen.specimenGuid"
  cell_type_column <- paste0("AIFI_", layer)
  count_column <- paste0(cell_type_column, "_count")
  cell_types <- gsub(" ", "_", as.character(data[[cell_type_column]]), fixed = TRUE)
  counts_data <- data.frame(
    sample = factor(as.character(data[[sample_column]]), levels = sample_ids),
    cell_type = factor(
      cell_types,
      levels = unique(cell_types[!is.na(cell_types)])
    ),
    count = as.numeric(data[[count_column]])
  )
  counts <- t(xtabs(count ~ sample + cell_type, data = counts_data))
  counts[!rownames(counts) %in% c("Platelet", "Erythrocyte"), sample_ids, drop = FALSE]
}

extract_scecoda_example_data <- function(project_root = .ecoda_example_project_root()) {
  source(file.path(project_root, "src", "utils", "datasets_io.R"), local = environment())
  source(file.path(project_root, "src", "utils", "seurat_utils.R"), local = environment())
  suppressPackageStartupMessages(library(dplyr))

  python_sys <- reticulate::import("sys", convert = FALSE)
  python_sys$dont_write_bytecode <- TRUE
  module_dir <- file.path(project_root, "src", "utils", "py")
  decoder <- reticulate::import_from_path(
    "h5ad_source_identity", path = module_dir, convert = FALSE
  )
  h5py <- reticulate::import("h5py", convert = FALSE)
  builtins <- reticulate::import_builtins(convert = FALSE)

  datasets <- read_datasets_json(
    path = file.path(project_root, "datasets.json"),
    view = "benchmark_analysis"
  )
  results_dir <- Sys.getenv(
    "ECODA_BENCHMARK_RESULTS_DIR",
    unset = file.path(project_root, "data", "benchmark", "results")
  )
  scratch_dir <- Sys.getenv(
    "HPC_SCRATCH_DIR",
    unset = file.path(path.expand("~"), "scratch", "ECODA_paper")
  )
  source_cell_types <- c("cell_type_low_res", "cell_type_high_res")
  benchmark_datasets <- c(
    "Adams", "Bassez", "Gongsharma_cmv_young_males", "Kfoury", "Kim",
    "Lee", "Pelka", "Smillie", "Stephenson", "Wu", "Zhang"
  )
  study_names <- c(
    "Adams", "Bassez", "Gongsharma_cmv_young_males", "Gongsharma_full",
    "Kfoury", "Kim", "Lee", "Pelka", "Smillie", "Stephenson", "Wu", "Zhang"
  )
  example_data <- setNames(vector("list", length(study_names)), study_names)

  for (dataset in benchmark_datasets) {
    config <- datasets[[dataset]]
    sample_column <- config$sample_col
    output_file <- get_view_h5ad_path(datasets, dataset, "benchmark_analysis")
    h5ad_path <- file.path(scratch_dir, dataset, "output", output_file)
    obs <- .ecoda_example_read_h5ad_obs(
      h5ad_path,
      sample_column,
      unlist(config[source_cell_types], use.names = FALSE),
      dataset %in% c("Lee", "Zhang"),
      decoder,
      h5py,
      builtins,
      get_ct_comp_df
    )

    if (dataset %in% c("Lee", "Zhang")) {
      cell_counts <- list(
        HiTME_LR = t(as.matrix(obs$layer1_counts)),
        HiTME_HR = .ecoda_example_result_counts(
          results_dir, dataset, "ECODA_HiTME_HR_layer2"
        ),
        scATOMIC_HR = .ecoda_example_result_counts(
          results_dir, dataset, "ECODA_scATOMIC_HR"
        )
      )
    } else {
      cell_counts <- list(
        authors_LR = .ecoda_example_result_counts(
          results_dir, dataset, "ECODA_authors_LR"
        ),
        authors_HR = .ecoda_example_result_counts(
          results_dir, dataset, "ECODA_authors_HR"
        ),
        HiTME_HR = .ecoda_example_result_counts(
          results_dir, dataset, "ECODA_HiTME_HR_layer2"
        ),
        scATOMIC_HR = .ecoda_example_result_counts(
          results_dir, dataset, "ECODA_scATOMIC_HR"
        )
      )
    }
    sample_ids <- rownames(obs$metadata)
    cell_counts <- lapply(cell_counts, function(counts) {
      counts[, sample_ids, drop = FALSE]
    })
    example_data[[dataset]] <- list(
      cell_counts = cell_counts,
      metadata = obs$metadata,
      biocond_colname = config$label_col
    )
  }

  gongsharma_dir <- Sys.getenv(
    "ECODA_GONGSHARMA_DATA_DIR",
    unset = file.path(
      Sys.getenv("NAS_SC_DIR"),
      "GongSharma_2024_PrePrintTBD", "data"
    )
  )
  gongsharma_data <- list(
    L1 = read.csv(
      file.path(gongsharma_dir, "sound_life_AIFI_L1_frequencies.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    L2 = read.csv(
      file.path(gongsharma_dir, "sound_life_AIFI_L2_frequencies.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    L3 = read.csv(
      file.path(gongsharma_dir, "sound_life_AIFI_L3_frequencies.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  )
  gongsharma_metadata_columns <- c(
    names(gongsharma_data$L3)[seq_len(14)],
    "total_cells", "scrna.lymphocyte_count", "bc.lymphocyte_count", "alc_ratio"
  )
  gongsharma_metadata <- unique(
    gongsharma_data$L3[, gongsharma_metadata_columns, drop = FALSE]
  )
  gongsharma_samples <- as.character(gongsharma_metadata$specimen.specimenGuid)
  gongsharma_metadata$age_group <- factor(
    ifelse(gongsharma_metadata$cohort.cohortGuid == "BR1", "Young", "Old"),
    levels = c("Young", "Old")
  )
  gongsharma_metadata$combined_cmv_age <- factor(
    paste0(
      "CMV",
      ifelse(as.character(gongsharma_metadata$subject.cmv) == "Positive", "+", "-"),
      " ", as.character(gongsharma_metadata$age_group)
    ),
    levels = c("CMV- Old", "CMV+ Old", "CMV- Young", "CMV+ Young")
  )
  rownames(gongsharma_metadata) <- gongsharma_samples
  example_data$Gongsharma_full <- list(
    cell_counts = list(
      authors_LR = .ecoda_example_gongsharma_counts(
        gongsharma_data$L1, "L1", gongsharma_samples
      ),
      authors_MR = .ecoda_example_gongsharma_counts(
        gongsharma_data$L2, "L2", gongsharma_samples
      ),
      authors_HR = .ecoda_example_gongsharma_counts(
        gongsharma_data$L3, "L3", gongsharma_samples
      )
    ),
    metadata = gongsharma_metadata,
    biocond_colname = "combined_cmv_age"
  )

  example_data <- example_data[study_names]
  save(example_data, file = file.path(project_root, "example_data.rda"))
  invisible(example_data)
}

if (sys.nframe() == 0L) extract_scecoda_example_data()
