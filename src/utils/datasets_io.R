read_datasets_json <- function(path = "datasets.json", view = NULL) {
  datasets <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  result <- list()

  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    views <- ds[["views"]]

    if (is.null(views)) next

    matched_views <- list()
    for (v_name in names(views)) {
      if (!is.null(view) && v_name != view) next

      v <- views[[v_name]]
      output_file <- v[["output_file_name"]]
      if (is.null(output_file)) next

      matched_views[[v_name]] <- list(
        view_name = v_name,
        input_file = v[["input_file_name"]],
        output_file = output_file,
        subset_vars = v[["subset_vars"]]
      )
    }

    if (length(matched_views) == 0) next

    first_v_name <- names(matched_views)[1]
    first_v <- matched_views[[first_v_name]]
    columns <- ds[["columns"]]

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
      low_res_ct_col = columns[["cell_type_low_res"]],
      hi_res_ct_col = columns[["cell_type_high_res"]],
      meta_cols_keep = ds[["meta_cols_keep"]],
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
