read_datasets_json <- function(path = "datasets.json", view = NULL) {
  datasets <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  result <- list()

  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    views <- ds[["views"]]
    base_columns <- ds[["columns"]]
    if (is.null(base_columns)) base_columns <- list()

    if (is.null(views)) next

    matched_views <- list()
    for (v_name in names(views)) {
      if (!is.null(view) && v_name != view) next

      v <- views[[v_name]]
      output_file <- v[["output_file_name"]]
      if (is.null(output_file)) next

      view_columns <- v[["columns"]]
      if (is.null(view_columns)) view_columns <- list()
      matched_views[[v_name]] <- list(
        view_name = v_name,
        input_file = v[["input_file_name"]],
        output_file = output_file,
        subset_vars = v[["subset_vars"]],
        columns = modifyList(base_columns, view_columns)
      )
    }

    if (length(matched_views) == 0) next

    first_v_name <- names(matched_views)[1]
    first_v <- matched_views[[first_v_name]]
    columns <- first_v[["columns"]]

    entry <- list(
      # dataset-level fields (order mirrors datasets.json)
      display_name = ds[["display_name"]],
      file_names = ds[["file_names"]],
      folder_name = ds[["folder_name"]],
      tissue = ds[["tissue"]],
      normal_tissue = ds[["normal_tissue"]],
      use_for_benchmark = ds[["use_for_benchmark"]],
      use_for_batch_effect = ds[["use_for_batch_effect"]],
      # columns (order mirrors datasets.json "columns")
      sample_col = columns[["sample"]],
      label_col = columns[["label"]],
      batch_col = columns[["batch"]],
      cell_type_low_res = columns[["cell_type_low_res"]],
      cell_type_high_res = columns[["cell_type_high_res"]],
      # Optional flat flag: annotation methods this dataset is NOT suitable
      # for (values: "hitme", "scatomic"). Absent/empty = suitable for all
      # (see AGENTS.md "Onboarding new datasets" / benchmark_pipeline.R).
      not_suitable_for_auto_annotation = if (is.null(ds[["not_suitable_for_auto_annotation"]])) {
        character(0)
      } else {
        ds[["not_suitable_for_auto_annotation"]]
      },
      # first-matching-view summary (R-compat; mirrors view entry order)
      view_name = first_v_name,
      input_file = first_v[["input_file"]],
      output_file = first_v[["output_file"]],
      subset_vars = first_v[["subset_vars"]],
      views = matched_views
    )

    result[[ds_name]] <- entry
  }

  return(result)
}

# Resolve the preprocessed view h5ad file name for a dataset (mirrors
# get_h5ad_path in benchmark_hpc_utils.R, but returns the bare file name —
# the notebook composes the full path from its own NAS mount). `datasets` is
# the read_datasets_json(view = ...) output; `ds` the datasets.json key.
get_view_h5ad_path <- function(datasets, ds, view) {
  entry <- datasets[[ds]]
  if (is.null(entry)) {
    stop("Dataset '", ds, "' not found in datasets.json (view '", view, "').")
  }
  views <- entry[["views"]]
  if (is.null(views) || is.null(views[[view]])) {
    stop("Dataset '", ds, "' has no '", view, "' view in datasets.json.")
  }
  out_file <- views[[view]][["output_file"]]
  if (is.null(out_file) || is.na(out_file) || out_file == "") {
    stop("Dataset '", ds, "' view '", view,
         "' has no output_file_name in datasets.json.")
  }
  return(out_file)
}
