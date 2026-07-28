read_datasets_json <- function(path = "datasets.json", view = NULL) {
  datasets <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  result <- list()

  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    views <- ds[["views"]]

    if (is.null(views)) next

    for (v_name in names(views)) {
      if (!is.null(view) && v_name != view) next

      v <- views[[v_name]]
      output_file <- v[["output_file_name"]]
      if (is.null(output_file)) next

      entry <- list(
        output_file = output_file,
        label_col = ds[["columns"]][["label"]],
        low_res_ct_col = ds[["columns"]][["cell_type_low_res"]],
        hi_res_ct_col = ds[["columns"]][["cell_type_high_res"]],
        display_name = ds[["display_name"]],
        tissue = ds[["tissue"]]
      )

      result[[ds_name]] <- entry
      break
    }
  }

  return(result)
}
