#!/usr/bin/env Rscript
# Application-specific scientific boundaries, without cohort files or plotting.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(root, "src/notebooks/ecoda_application_cancer_clustering.R"))

cohorts <- c("Yan", "Zhang", "Bassez", "Pelka")
ids <- paste0("s", seq_len(11L))
types <- c("CD4_Naive", "CD8_Effector", "Macrophages", "Bcell",
           "Rare90", "Rare100", "MissingInPelka", "Neutrophils", "CAFs")
counts <- matrix(0, length(ids), length(types), dimnames = list(ids, types))
counts[, "CD4_Naive"] <- c(4, 5, 8, 12, 6, 2, 9, 11, 3, 7, 9)
counts[, "CD8_Effector"] <- c(3, 8, 2, 7, 11, 5, 4, 6, 13, 9, 0)
counts[, "Macrophages"] <- c(2, 2, 9, 4, 3, 12, 7, 5, 8, 6, 0)
counts[, "Bcell"] <- c(1, 4, 5, 2, 8, 6, 3, 9, 7, 11, 0)
counts["s10", "Rare90"] <- 3
counts[2:10, "MissingInPelka"] <- 2
# Excluded cells must not rescue the sample with only nine immune cells.
counts[, c("Neutrophils", "CAFs")] <- 1000
metadata <- data.frame(
  Sample = ids,
  `NAC response` = rep(c("pCR", "RD", "NA"), length.out = length(ids)),
  Tissue = "tumor", Clinical.efficacy. = rep(c("PR", "SD"), length.out = length(ids)),
  expansion = rep(c("E", "NE"), length.out = length(ids)),
  MMRStatus = rep(c("MMRd", "MMRp"), length.out = length(ids)),
  row.names = ids, check.names = FALSE, stringsAsFactors = FALSE
)
inputs <- setNames(lapply(cohorts, function(cohort) {
  list(counts = counts, metadata = metadata)
}), cohorts)
inputs$Pelka$metadata$Tissue <- "T"
inputs$Pelka$counts[, "MissingInPelka"] <- 0

result <- cancer_prepare(inputs)
expected_types <- c("CD4_Naive", "CD8_Effector", "Macrophages", "Bcell", "Rare90")
stopifnot(
  setequal(result$shared_cells, expected_types),
  all(vapply(result$counts, function(x) identical(rownames(x), ids[1:10]), logical(1))),
  all(vapply(result$counts, function(x) "Rare90" %in% colnames(x), logical(1))),
  all(vapply(result$counts, function(x) !any(c("Rare100", "Neutrophils", "CAFs",
                                            "MissingInPelka") %in% colnames(x)), logical(1)))
)
for (cohort in cohorts) {
  actual_counts <- as.matrix(result$counts[[cohort]])
  expected <- log(actual_counts + 0.5)
  expected <- sweep(expected, 1L, rowMeans(expected), "-")
  stopifnot(
    identical(rownames(result$metadata[[cohort]]), rownames(actual_counts)),
    identical(colnames(result$clr[[cohort]]), colnames(actual_counts)),
    isTRUE(all.equal(as.matrix(result$clr[[cohort]]), expected, tolerance = 1e-12)),
    max(abs(rowSums(as.matrix(result$clr[[cohort]])))) < 1e-12
  )
  # Scores and loadings reconstruct the same centered CLR used for plotting.
  pca <- result$pca[[cohort]]
  reconstructed <- pca$x %*% t(pca$rotation)
  centered <- sweep(expected, 2L, colMeans(expected), "-")
  stopifnot(isTRUE(all.equal(unname(reconstructed), unname(centered), tolerance = 1e-10)))
}

# The heatmap must center these same CLR values, without variance scaling.
grDevices::pdf(NULL)
heatmap <- tryCatch(
  cancer_heatmap(result$clr$Bassez, result$metadata$Bassez, "Bassez"),
  finally = grDevices::dev.off()
)
expected_heatmap <- t(sweep(result$clr$Bassez, 2L, colMeans(result$clr$Bassez), "-"))
stopifnot(isTRUE(all.equal(heatmap$matrix, expected_heatmap, tolerance = 1e-12)))

# Pelka affects the shared feature universe even when breast figures are used.
with_pelka_feature <- inputs
with_pelka_feature$Pelka$counts[2:10, "MissingInPelka"] <- 2
stopifnot("MissingInPelka" %in% cancer_prepare(with_pelka_feature)$shared_cells)

# A response-label permutation must not influence any PCA fit or orientation.
permuted <- inputs
for (cohort in cohorts) {
  for (column in c("NAC response", "Clinical.efficacy.", "expansion", "MMRStatus")) {
    values <- permuted[[cohort]]$metadata[[column]]
    permuted[[cohort]]$metadata[[column]] <- c(values[-1L], values[1L])
  }
}
permuted_result <- cancer_prepare(permuted)
for (cohort in cohorts) {
  stopifnot(
    identical(result$clr[[cohort]], permuted_result$clr[[cohort]]),
    identical(result$pca[[cohort]]$x, permuted_result$pca[[cohort]]$x),
    identical(result$pca[[cohort]]$rotation, permuted_result$pca[[cohort]]$rotation)
  )
}

must_reject <- function(bad, label) {
  error <- tryCatch({ cancer_prepare(bad); NULL }, error = identity)
  if (!inherits(error, "error")) stop("Expected rejection: ", label)
}
bad <- inputs
bad$Yan$counts[1, 1] <- 1.5
must_reject(bad, "noninteger count")
bad <- inputs
bad$Yan$counts[1, 1] <- -1
must_reject(bad, "negative count")
bad <- inputs
bad$Yan$counts[1, 1] <- NA_real_
must_reject(bad, "missing count")
bad <- inputs
rownames(bad$Yan$counts)[2] <- rownames(bad$Yan$counts)[1]
must_reject(bad, "duplicate sample ID")
bad <- inputs
rownames(bad$Yan$counts)[1] <- ""
must_reject(bad, "blank sample ID")
bad <- inputs
rownames(bad$Yan$metadata)[1] <- "unexpected"
must_reject(bad, "unmatched metadata ID")
bad <- inputs
colnames(bad$Yan$counts)[2] <- colnames(bad$Yan$counts)[1]
must_reject(bad, "duplicate cell type")

cat("PASS cancer application filtering, four-cohort intersection, CLR, PCA and input validation\n")
