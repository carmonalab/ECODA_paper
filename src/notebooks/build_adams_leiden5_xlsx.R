#!/usr/bin/env Rscript
# Build 4-sheet Supplementary Excel correspondence workbook using writexl:
# Authors_HR (38) x Leiden_res_5 (92) on Adams IPF (305,701 cells).

root <- normalizePath(".")
counts_csv <- file.path(root, "data/derived/adams_annotation_correspondence_20261023/counts.csv")
top_csv <- file.path(root, "data/derived/adams_leiden5_20261121/leiden5_top_enriched_markers.csv")
hitme_csv <- file.path(root, "data/derived/adams_leiden5_20261121/leiden5_hitme_correspondence.csv")

df_counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
clusters_sorted <- colnames(df_counts)[order(as.integer(colnames(df_counts)))]
authors_sorted <- sort(rownames(df_counts))
df_counts <- df_counts[authors_sorted, clusters_sorted, drop = FALSE]
mat_counts <- as.matrix(df_counts)
stopifnot(nrow(mat_counts) == 38L, ncol(mat_counts) == 92L, sum(mat_counts) == 305701)

row_totals <- rowSums(mat_counts)
col_totals <- colSums(mat_counts)
mat_row_pct <- round(100 * sweep(mat_counts, 1, row_totals, "/"), 2)
mat_col_pct <- round(100 * sweep(mat_counts, 2, col_totals, "/"), 2)

top_markers <- read.csv(top_csv, stringsAsFactors = FALSE)
hitme_df <- read.csv(hitme_csv, stringsAsFactors = FALSE)
rownames(hitme_df) <- as.character(hitme_df$cluster)

summary_rows <- lapply(clusters_sorted, function(cl) {
  vals <- sort(mat_counts[, cl], decreasing = TRUE)
  tot <- col_totals[[cl]]
  top1_a <- names(vals)[1]
  top1_n <- unname(vals[1])
  top1_pct <- round(100 * top1_n / tot, 2)
  top2_a <- if (length(vals) > 1 && vals[2] > 0) names(vals)[2] else ""
  top2_pct <- if (nzchar(top2_a)) round(100 * unname(vals[2]) / tot, 2) else 0
  h <- hitme_df[cl, , drop = FALSE]
  m <- utils::head(top_markers[as.character(top_markers$cluster) == cl, ], 5)
  data.frame(
    Leiden_res_5_Cluster = paste0("cl_", cl),
    Cluster_ID = as.integer(cl),
    N_Cells = as.integer(tot),
    Top1_Authors_HR = top1_a,
    Top1_Authors_HR_Cells = as.integer(top1_n),
    Top1_Authors_HR_Pct = top1_pct,
    Top2_Authors_HR = top2_a,
    Top2_Authors_HR_Pct = top2_pct,
    Dominant_HiTME_Layer2 = as.character(h$hitme_top1[[1]]),
    Dominant_HiTME_Pct = as.numeric(h$hitme_top1_pct[[1]]),
    Top5_Enriched_Markers = paste(m$gene, collapse = ", "),
    stringsAsFactors = FALSE
  )
})
df_summary <- do.call(rbind, summary_rows)

as_matrix_sheet <- function(mat, add_col_totals = TRUE, add_row_totals = TRUE, total_label = "Total_Cells") {
  out <- as.data.frame(mat, check.names = FALSE)
  colnames(out) <- paste0("cl_", colnames(mat))
  out <- cbind(Authors_HR = rownames(mat), out, stringsAsFactors = FALSE)
  if (add_row_totals) {
    out[[total_label]] <- if (identical(total_label, "Total_Cells")) as.numeric(rowSums(mat)) else round(rowSums(mat), 2)
  }
  if (add_col_totals) {
    tot_row <- as.list(c("Total_Cells", as.numeric(colSums(mat))))
    if (add_row_totals) tot_row[[length(tot_row) + 1L]] <- sum(mat)
    names(tot_row) <- colnames(out)
    out <- rbind(out, as.data.frame(tot_row, check.names = FALSE, stringsAsFactors = FALSE))
  }
  out
}

sheets <- list(
  Leiden5_Cluster_Summary = df_summary,
  Counts_38x92 = as_matrix_sheet(mat_counts, TRUE, TRUE, "Total_Cells"),
  RowPct_Author_in_Cluster = as_matrix_sheet(mat_row_pct, FALSE, TRUE, "Row_Sum_Pct"),
  ColPct_Cluster_Composition = as_matrix_sheet(mat_col_pct, TRUE, FALSE, "Col_Sum_Pct")
)

out_paths <- c(
  file.path(root, "data/derived/adams_leiden5_20261121/Supp_Table_Adams_AuthorsHR_vs_LeidenRes5.xlsx"),
  file.path(root, "plots/Supp_Table_Adams_AuthorsHR_vs_LeidenRes5.xlsx")
)
for (p in out_paths) {
  writexl::write_xlsx(sheets, path = p)
  cat("Wrote:", p, "\n")
}
