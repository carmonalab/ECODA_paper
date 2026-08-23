# Shared validation and replacement helpers for the Bassez subtype patch.

BASSEZ_MISSING_ANNOTATIONS <- c("", "na", "nan", "none", "unknown")

bassez_missing_annotation <- function(values) {
  values_chr <- as.character(values)
  normalized <- tolower(trimws(values_chr))
  is.na(values_chr) | normalized %in% BASSEZ_MISSING_ANNOTATIONS
}

bassez_fill_cell_subtype <- function(metadata) {
  if (!is.data.frame(metadata)) {
    stop("Bassez metadata must be a data.frame.")
  }

  required <- c("cellType", "cellSubType")
  missing_columns <- setdiff(required, colnames(metadata))
  if (length(missing_columns) > 0) {
    stop(
      "Bassez metadata is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  replacement_rows <- bassez_missing_annotation(metadata$cellSubType)
  invalid_fallback <- replacement_rows & bassez_missing_annotation(metadata$cellType)
  if (any(invalid_fallback)) {
    bad_rows <- which(invalid_fallback)
    row_labels <- rownames(metadata)[bad_rows]
    if (is.null(row_labels)) row_labels <- as.character(bad_rows)
    stop(
      "Cannot fill Bassez cellSubType: invalid cellType fallback at row(s): ",
      paste(row_labels, collapse = ", ")
    )
  }

  filled_subtype <- as.character(metadata$cellSubType)
  filled_subtype[replacement_rows] <- as.character(metadata$cellType)[replacement_rows]
  metadata$cellSubType <- factor(filled_subtype)
  metadata
}
