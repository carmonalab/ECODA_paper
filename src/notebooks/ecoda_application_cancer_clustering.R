# Application-specific port of Analysis_unsup_strat.rmd. These functions read no
# cohort files and never fit preprocessing or annotation models.

cancer_verify_file <- function(path) {
  sidecar <- paste0(path, ".md5")
  if (!file.exists(path) || file.info(path)$size <= 0 || !file.exists(sidecar)) {
    stop("Missing/empty application input or checksum: ", path)
  }
  expected <- strsplit(readLines(sidecar, warn = FALSE)[1], "[[:space:]]+")[[1]][1]
  if (!identical(unname(tools::md5sum(path)), expected)) stop("Checksum mismatch: ", path)
  invisible(expected)
}

cancer_read_inputs <- function(input_dir) {
  inputs <- provenance <- reconciliation <- list()
  for (cohort in c("Bassez", "Zhang", "Pelka")) {
    snapshot_path <- file.path(input_dir, cohort, "obs_snapshot.json")
    bundle_path <- file.path(input_dir, cohort, paste0(cohort, "_ECODA_HiTME_HR_layer2.rds"))
    cancer_verify_file(snapshot_path)
    cancer_verify_file(bundle_path)
    snapshot <- jsonlite::fromJSON(snapshot_path)
    if (!identical(snapshot$schema_version, 1L) || !identical(snapshot$dataset, cohort)) {
      stop("Unexpected compact snapshot schema: ", cohort)
    }
    observed <- as.matrix(snapshot$counts)
    dimnames(observed) <- list(snapshot$samples, snapshot$cell_types)
    meta <- as.data.frame(snapshot$metadata, check.names = FALSE, stringsAsFactors = FALSE)
    cancer_check_ids(meta$sample, paste(cohort, "snapshot sample IDs"))
    cancer_check_ids(meta$canonical_sample, paste(cohort, "canonical sample IDs"))
    rownames(meta) <- meta$sample
    bundle <- readRDS(bundle_path)
    raw <- as.matrix(bundle$counts) - 0.5
    mapped <- match(rownames(raw), meta$canonical_sample)
    if (anyNA(mapped)) stop(cohort, ": benchmark sample lacks configured/canonical ID mapping")
    rownames(raw) <- meta$sample[mapped]
    if (!setequal(rownames(raw), rownames(observed))) stop(cohort, ": legacy/current sample mismatch")
    common <- intersect(colnames(raw), colnames(observed))
    difference <- observed[rownames(raw), common, drop = FALSE] - raw[, common, drop = FALSE]
    reconciliation[[cohort]] <- list(
      legacy_bundle_used_as_input = FALSE,
      types_current_only = setdiff(colnames(observed), colnames(raw)),
      types_legacy_only = setdiff(colnames(raw), colnames(observed)),
      max_common_count_difference = max(abs(difference)),
      common_count_difference = difference
    )
    # Current obs is the explicit source of raw counts. The stale legacy
    # benchmark artifact is comparison evidence only, never a fallback input.
    inputs[[cohort]] <- list(counts = observed, metadata = meta[rownames(observed), , drop = FALSE])
    provenance[[cohort]] <- snapshot$provenance
    provenance[[cohort]]$count_source <- "fresh aggregate of frozen current obs; legacy RDS is comparison only"
  }
  yan_path <- file.path(input_dir, "Yan", "ecoda_data.rds")
  yan_provenance <- file.path(input_dir, "Yan", "provenance.json")
  cancer_verify_file(yan_path)
  cancer_verify_file(yan_provenance)
  yan <- readRDS(yan_path)
  counts <- as.data.frame(yan$cell_counts$HiTME$layer2, check.names = FALSE)
  cancer_check_ids(as.character(counts$sample), "Yan count sample IDs")
  rownames(counts) <- as.character(counts$sample)
  counts$sample <- NULL
  meta <- as.data.frame(yan$metadata, check.names = FALSE)
  cancer_check_ids(as.character(meta$sample), "Yan metadata sample IDs")
  rownames(meta) <- as.character(meta$sample)
  inputs$Yan <- list(counts = as.matrix(counts), metadata = meta)
  provenance$Yan <- jsonlite::fromJSON(yan_provenance)
  list(inputs = inputs[c("Yan", "Zhang", "Bassez", "Pelka")], provenance = provenance,
    legacy_reconciliation = reconciliation)
}

cancer_save_rds <- function(object, path) {
  if (file.exists(path)) {
    cancer_verify_file(path)
    if (!identical(readRDS(path), object)) stop("Existing evidence differs; choose a new run ID: ", path)
    return(invisible(path))
  }
  temporary <- tempfile(pattern = ".cancer-evidence-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Cannot publish evidence: ", path)
  sidecar <- tempfile(pattern = ".cancer-checksum-", tmpdir = dirname(path))
  on.exit(unlink(sidecar), add = TRUE)
  writeLines(paste0(unname(tools::md5sum(path)), "  ", basename(path)), sidecar)
  if (!file.rename(sidecar, paste0(path, ".md5"))) stop("Cannot publish evidence checksum")
  invisible(path)
}

cancer_immune_excluded <- function(cell_types) {
  cell_types %in% c(
    "Erythrocyte", "Neutrophils", "Platelet", "Endothelial", "Epithelial",
    "Fibroblast", "abnormal cell", "fibroblast", "perivascular cell",
    "endothelial cell"
  ) | grepl("__|endo-|fibro-|VSMC-|peri-|CAFs|Unresolved", cell_types)
}

cancer_clr <- function(counts) {
  # Adding to all counts is intentional, including cohorts with no zeros.
  logs <- log(counts + 0.5)
  sweep(logs, 1L, rowMeans(logs), "-")
}

cancer_check_ids <- function(ids, context) {
  if (is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop(context, " must be unique, nonmissing, nonempty identifiers")
  }
}

cancer_align_pca <- function(pca, reference = NULL) {
  anchors <- list(
    c("CD4|CD8|Bcell|Plasma|NK|gdT", "Mono|Mac|DC|panDC|MoMac"),
    c("TEX|CTL_Exh|CXCL9", "Naive")
  )
  for (pc in 1:2) {
    loadings <- pca$rotation[, pc]
    score <- sum(loadings[grepl(anchors[[pc]][1], names(loadings), ignore.case = TRUE)]) -
      sum(loadings[grepl(anchors[[pc]][2], names(loadings), ignore.case = TRUE)])
    if (score < 0) {
      pca$rotation[, pc] <- -pca$rotation[, pc]
      pca$x[, pc] <- -pca$x[, pc]
    }
  }
  if (!is.null(reference) && stats::cor(pca$rotation[, 2], reference$rotation[, 2]) < 0) {
    pca$rotation[, 2] <- -pca$rotation[, 2]
    pca$x[, 2] <- -pca$x[, 2]
  }
  pca
}

cancer_prepare <- function(inputs) {
  cohorts <- c("Yan", "Zhang", "Bassez", "Pelka")
  if (!is.list(inputs) || !setequal(names(inputs), cohorts) || anyDuplicated(names(inputs))) {
    stop("inputs must contain exactly Yan, Zhang, Bassez and Pelka")
  }
  counts <- metadata <- sample_audits <- type_audits <- setNames(vector("list", 4), cohorts)
  for (cohort in cohorts) {
    input <- inputs[[cohort]]
    mat <- input$counts
    meta <- input$metadata
    if (!is.matrix(mat) || !is.numeric(mat) || !length(mat) ||
        any(!is.finite(mat)) || any(mat < 0) || any(abs(mat - round(mat)) > 1e-8)) {
      stop(cohort, ": counts must be a nonempty matrix of finite nonnegative integers")
    }
    cancer_check_ids(rownames(mat), paste(cohort, "count sample IDs"))
    cancer_check_ids(colnames(mat), paste(cohort, "cell types"))
    if (!is.data.frame(meta)) stop(cohort, ": metadata must be a data.frame")
    cancer_check_ids(rownames(meta), paste(cohort, "metadata sample IDs"))
    if (!setequal(rownames(mat), rownames(meta))) stop(cohort, ": metadata/count sample mismatch")
    meta <- meta[rownames(mat), , drop = FALSE]
    tissue_keep <- rep(TRUE, nrow(mat))
    if (cohort %in% c("Zhang", "Pelka")) {
      if (!"Tissue" %in% names(meta)) stop(cohort, ": Tissue is required")
      tissue_keep <- !is.na(meta$Tissue) & meta$Tissue == if (cohort == "Zhang") "tumor" else "T"
    }
    excluded <- cancer_immune_excluded(colnames(mat))
    immune_sums <- rowSums(mat[, !excluded, drop = FALSE])
    keep <- tissue_keep & immune_sums >= 10
    if (sum(keep) < 3) stop(cohort, ": fewer than three samples after tissue/immune filtering")
    immune <- mat[keep, !excluded, drop = FALSE]
    zero_fraction <- colMeans(immune == 0)
    sparse <- zero_fraction > 0.9
    counts[[cohort]] <- immune[, !sparse, drop = FALSE]
    metadata[[cohort]] <- meta[keep, , drop = FALSE]
    sample_audits[[cohort]] <- data.frame(
      cohort = cohort, sample = rownames(mat), input_cells = rowSums(mat),
      eligible_immune_cells = immune_sums, shared_cells = NA_real_, retained = keep,
      reason = ifelse(!tissue_keep, "tissue", ifelse(!keep, "fewer_than_10_immune_cells", "retained")),
      stringsAsFactors = FALSE, row.names = NULL
    )
    type_audits[[cohort]] <- data.frame(
      cohort = cohort, cell_type = colnames(mat), immune_excluded = excluded,
      zero_fraction = unname(zero_fraction[match(colnames(mat), names(zero_fraction))]),
      sparse_excluded = colnames(mat) %in% names(zero_fraction)[sparse],
      in_shared = FALSE, stringsAsFactors = FALSE
    )
  }
  shared <- Reduce(intersect, lapply(counts, colnames))
  if (length(shared) < 3) stop("Fewer than three shared immune cell types")
  clr <- pca <- setNames(vector("list", 4), cohorts)
  for (cohort in cohorts) {
    counts[[cohort]] <- counts[[cohort]][, shared, drop = FALSE]
    clr[[cohort]] <- cancer_clr(counts[[cohort]])
    fit <- stats::prcomp(clr[[cohort]], center = TRUE, scale. = FALSE)
    if (length(fit$sdev) < 2 || !is.finite(fit$sdev[2]) ||
        fit$sdev[2] <= max(fit$sdev[1], 1) * 1e-12) {
      stop(cohort, ": shared CLR matrix has fewer than two nondegenerate axes")
    }
    pca[[cohort]] <- cancer_align_pca(fit, if (cohort == "Yan") NULL else pca$Yan)
    metadata[[cohort]]$PC1 <- pca[[cohort]]$x[, 1]
    metadata[[cohort]]$PC2 <- pca[[cohort]]$x[, 2]
    sample_audits[[cohort]]$shared_cells[match(rownames(counts[[cohort]]), sample_audits[[cohort]]$sample)] <-
      rowSums(counts[[cohort]])
    types <- type_audits[[cohort]]
    types$in_shared <- types$cell_type %in% shared
    types$reason <- ifelse(types$immune_excluded, "immune_exclusion",
      ifelse(types$sparse_excluded, "more_than_90_percent_zeros",
        ifelse(!types$in_shared, "four_cohort_intersection", "retained")))
    type_audits[[cohort]] <- types
  }
  yan_response <- intersect(c("NAC response", "NAC.response"), names(metadata$Yan))[1]
  if (is.na(yan_response)) stop("Yan: NAC response or NAC.response is required")
  if (!"Clinical.efficacy." %in% names(metadata$Zhang)) stop("Zhang: Clinical.efficacy. is required")
  if (!"expansion" %in% names(metadata$Bassez)) stop("Bassez: expansion is required")
  if (!"MMRStatus" %in% names(metadata$Pelka)) stop("Pelka: MMRStatus is required")
  y <- as.character(metadata$Yan[[yan_response]])
  # The literal NA group means early dropout/Low by the user's interpretation;
  # actual missing values remain missing. Neither enters feature construction.
  metadata$Yan$Response <- factor(ifelse(y == "pCR", "High",
    ifelse(y == "RD", "Medium", ifelse(y == "NA", "Low", NA))),
    levels = c("High", "Medium", "Low"))
  z <- as.character(metadata$Zhang$Clinical.efficacy.)
  metadata$Zhang$Response <- factor(ifelse(z %in% c("PR", "CR"), "High",
    ifelse(z == "SD", "Medium", ifelse(z == "PD", "Low", NA))),
    levels = c("High", "Medium", "Low"))
  b <- as.character(metadata$Bassez$expansion)
  metadata$Bassez$Response <- factor(ifelse(b == "E", "High", ifelse(b == "NE", "Low", NA)),
    levels = c("High", "Low"))
  metadata$Pelka$MMRStatus <- factor(metadata$Pelka$MMRStatus, levels = c("MMRd", "MMRp"))
  consensus <- setNames(lapply(1:2, function(pc) {
    values <- vapply(pca[c("Yan", "Zhang", "Bassez")], function(fit) fit$rotation[shared, pc], numeric(length(shared)))
    out <- data.frame(cell_type = shared, values, mean_loading = rowMeans(values),
      sem = apply(values, 1, stats::sd) / sqrt(3), row.names = NULL, check.names = FALSE)
    out[order(out$mean_loading), , drop = FALSE]
  }), c("PC1", "PC2"))
  all_types <- unique(unlist(lapply(type_audits, function(x) x$cell_type)))
  harmonization <- data.frame(cell_type = all_types, stringsAsFactors = FALSE)
  for (cohort in cohorts) {
    types <- type_audits[[cohort]]
    harmonization[[cohort]] <- all_types %in% types$cell_type[!types$immune_excluded & !types$sparse_excluded]
  }
  harmonization$in_shared <- all_types %in% shared
  list(counts = counts, metadata = metadata, clr = clr, pca = pca, shared_cells = shared,
    audit = list(samples = do.call(rbind, sample_audits), cell_types = do.call(rbind, type_audits),
      harmonization = harmonization), consensus = consensus,
    settings = list(min_immune_cells = 10, max_zero_fraction = 0.9, pseudocount = 0.5,
      intersection_cohorts = cohorts, pca_center = TRUE, pca_scale = FALSE))
}

cancer_response_colors <- c(High = "#2ca02c", Medium = "#1f77b4", Low = "#d62728")

cancer_biplot <- function(fit, meta, title, legend) {
  variance <- fit$sdev^2 / sum(fit$sdev^2) * 100
  top <- head(order(rowSums(fit$rotation[, 1:2, drop = FALSE]^2), decreasing = TRUE), 8)
  multiplier <- min(vapply(1:2, function(pc) diff(range(fit$x[, pc])) / diff(range(fit$rotation[, pc])), 0)) * 0.7
  arrows <- data.frame(feature = rownames(fit$rotation)[top],
    x = fit$rotation[top, 1] * multiplier, y = fit$rotation[top, 2] * multiplier)
  ggplot2::ggplot(meta, ggplot2::aes(x = PC1, y = PC2, color = Response)) +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
      inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.2, "cm")), color = "black", linewidth = 0.6) +
    ggrepel::geom_text_repel(data = arrows, ggplot2::aes(x = x, y = y, label = feature),
      inherit.aes = FALSE, color = "black", size = 3.2, seed = 42) +
    ggplot2::scale_color_manual(values = cancer_response_colors, name = legend) +
    ggplot2::theme_classic() + ggplot2::coord_fixed() +
    ggplot2::labs(title = paste0(title, " - Pre-treatment"), subtitle = "Immune cells only - High granularity",
      x = paste0("PC1 (", formatC(variance[1], digits = 2, format = "fg"), "%)"),
      y = paste0("PC2 (", formatC(variance[2], digits = 2, format = "fg"), "%)")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8.5), legend.position = "bottom")
}

cancer_response_box <- function(meta, pc, title, label) {
  valid <- meta[!is.na(meta$Response), , drop = FALSE]
  comparisons <- if (nlevels(meta$Response) > 2) list(c("High", "Medium"), c("High", "Low"), c("Medium", "Low")) else list(c("High", "Low"))
  if (grepl("Zhang", title)) comparisons <- list(c("High", "Medium"))
  significant <- Filter(function(pair) {
    a <- valid[[pc]][valid$Response == pair[1]]
    b <- valid[[pc]][valid$Response == pair[2]]
    length(a) > 0 && length(b) > 0 && isTRUE(suppressWarnings(stats::wilcox.test(a, b, exact = FALSE)$p.value <= 0.05))
  }, comparisons)
  p <- ggplot2::ggplot(valid, ggplot2::aes(x = Response, y = .data[[pc]], fill = Response)) +
    ggplot2::geom_boxplot(width = 0.38, outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 42), size = 1, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = cancer_response_colors) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.22)), labels = function(x) format(x, digits = 2)) +
    ggplot2::theme_classic() + ggplot2::labs(title = paste(title, "-", pc), x = label, y = pc) +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 10, face = "bold"),
      axis.text = ggplot2::element_text(size = 8.5), axis.title = ggplot2::element_text(size = 9.5))
  if (length(significant)) p <- p + ggpubr::stat_compare_means(comparisons = significant,
    method = "wilcox.test", label = "p.signif", hide.ns = TRUE)
  p
}

cancer_consensus_plot <- function(values, pc) {
  cohorts <- c("Yan", "Zhang", "Bassez")
  long <- do.call(rbind, lapply(cohorts, function(cohort) data.frame(
    cell_type = factor(values$cell_type, levels = values$cell_type), cohort = cohort, loading = values[[cohort]])))
  ggplot2::ggplot(long, ggplot2::aes(x = loading, y = cell_type, fill = cohort)) +
    ggplot2::geom_col(position = "stack", color = "black", linewidth = 0.25, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = c(Yan = "gray30", Zhang = "gray60", Bassez = "gray85"), name = "Cohort") +
    ggplot2::theme_classic() + ggplot2::labs(title = if (pc == "PC1") "Lymphoid-Myeloid Axis (PC1)" else "Immune Activation Axis (PC2)",
      x = paste("Cumulative", pc, "Eigenvector Loading"), y = "") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9, face = "bold"))
}

cancer_heatmap <- function(clr, meta, cohort) {
  annotation <- if (cohort == "Pelka") "MMRStatus" else "Response"
  colors <- if (cohort == "Pelka") c(MMRd = "#1b9e77", MMRp = "#d95f02") else cancer_response_colors
  annotations <- droplevels(meta[, annotation, drop = FALSE])
  colors <- colors[names(colors) %in% as.character(annotations[[1]])]
  centered <- t(sweep(clr, 2, colMeans(clr), "-"))
  cuts <- switch(cohort, Yan = c(4L, 3L), Zhang = c(3L, 3L), Bassez = c(3L, 3L), Pelka = c(2L, 2L))
  if (nrow(centered) < cuts[1] || ncol(centered) < cuts[2]) stop(cohort, ": insufficient dimensions for heatmap cuts")
  heat <- pheatmap::pheatmap(centered, annotation_col = annotations,
    annotation_colors = setNames(list(colors), annotation), scale = "none",
    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2", cutree_rows = cuts[1], cutree_cols = cuts[2],
    show_colnames = FALSE, silent = TRUE)
  matrix_cell <- heat$gtable$layout[heat$gtable$layout$name == "matrix", ]
  # Reserve footer space inside the original 1npc canvas; adding absolute
  # height without subtracting it clips both the tree and the Samples label.
  heat$gtable$heights[matrix_cell$t] <- heat$gtable$heights[matrix_cell$t] - grid::unit(18, "pt")
  heat$gtable <- gtable::gtable_add_rows(heat$gtable, grid::unit(18, "pt"))
  heat$gtable <- gtable::gtable_add_grob(heat$gtable,
    grid::textGrob("Samples", gp = grid::gpar(fontface = "bold", fontsize = 10)),
    t = nrow(heat$gtable), l = matrix_cell$l, r = matrix_cell$r, clip = "off", name = "samples_label")
  list(plot = heat, matrix = centered,
    clusters = list(rows = stats::cutree(heat$tree_row, cuts[1]), columns = stats::cutree(heat$tree_col, cuts[2]),
      tree_rows = heat$tree_row, tree_columns = heat$tree_col))
}

cancer_correlation <- function(pca) {
  order <- c("Yan", "Bassez", "Zhang", "Pelka")
  values <- do.call(cbind, lapply(1:2, function(pc) {
    out <- vapply(pca[order], function(fit) fit$rotation[, pc], numeric(nrow(pca$Yan$rotation)))
    colnames(out) <- paste0(order, " (PC", pc, ")")
    out
  }))
  correlations <- stats::cor(values)
  tree <- stats::hclust(stats::as.dist(1 - correlations), method = "ward.D2")
  axes <- colnames(values)[tree$order]
  p_values <- matrix(0, 8, 8, dimnames = dimnames(correlations))
  cells <- list()
  for (j in 1:7) for (i in (j + 1):8) {
    test <- stats::cor.test(values[, axes[j]], values[, axes[i]], method = "pearson")
    p <- test$p.value
    p_values[axes[i], axes[j]] <- p_values[axes[j], axes[i]] <- p
    stars <- if (p <= 0.0001) "****" else if (p <= 0.001) "***" else if (p <= 0.01) "**" else if (p <= 0.05) "*" else ""
    cells[[length(cells) + 1L]] <- data.frame(x = axes[j], y = axes[i], r = round(unname(test$estimate), 2), stars = stars)
  }
  cells <- do.call(rbind, cells)
  cells$x <- factor(cells$x, levels = axes[1:7])
  cells$y <- factor(cells$y, levels = rev(axes[2:8]))
  plot <- ggplot2::ggplot(cells, ggplot2::aes(x = x, y = y, fill = r)) +
    ggplot2::geom_tile(color = "white", linewidth = 1.3, width = 0.95, height = 0.95) +
    ggplot2::geom_text(ggplot2::aes(label = stars), size = 4.2, fontface = "bold", vjust = 0.75) +
    ggplot2::scale_fill_gradient2(low = "#1f77b4", mid = "white", high = "#d62728", midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::theme_classic() + ggplot2::labs(x = "", y = "",
      title = "Clustered Cross-Cohort PC1 & PC2 Loadings Correlation Matrix",
      subtitle = "* p <= 0.05, ** p <= 0.01, *** p <= 0.001, **** p <= 0.0001 (unadjusted)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9.5, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 9.5, face = "bold"), axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 7.5, hjust = 0.5), legend.title = ggplot2::element_text(size = 9))
  list(plot = plot, matrix = correlations, p_values = p_values, order = axes, hclust = tree, loadings = values)
}

cancer_save_pdf <- function(plot, path, width, height) {
  temporary <- tempfile(pattern = ".cancer-plot-", tmpdir = dirname(path), fileext = ".pdf")
  on.exit(unlink(temporary), add = TRUE)
  ggplot2::ggsave(temporary, plot = plot, width = width, height = height, units = "in", limitsize = FALSE)
  if (!file.rename(temporary, path)) stop("Could not publish PDF: ", path)
  checksum <- tempfile(pattern = ".cancer-md5-", tmpdir = dirname(path))
  on.exit(unlink(checksum), add = TRUE)
  writeLines(paste0(unname(tools::md5sum(path)), "  ", basename(path)), checksum)
  if (!file.rename(checksum, paste0(path, ".md5"))) stop("Could not publish PDF checksum: ", path)
  path
}

cancer_plot_all <- function(analysis, output_dir) {
  packages <- c("ggplot2", "ggrepel", "ggpubr", "patchwork", "pheatmap", "gtable")
  for (package in packages) if (!requireNamespace(package, quietly = TRUE)) stop("Missing installed package: ", package)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory")
  titles <- c(Yan = "Yan TNBC", Zhang = "Zhang TNBC", Bassez = "Bassez BC")
  response_labels <- c(Yan = "Chemo resp.", Zhang = "ICB resp.", Bassez = "ICB resp.")
  biplots <- boxes <- list()
  for (cohort in names(titles)) {
    biplots[[cohort]] <- cancer_biplot(analysis$pca[[cohort]], analysis$metadata[[cohort]], titles[[cohort]], response_labels[[cohort]])
    for (pc in c("PC1", "PC2")) boxes[[paste(cohort, pc)]] <- cancer_response_box(analysis$metadata[[cohort]], pc, titles[[cohort]], response_labels[[cohort]])
  }
  row1 <- patchwork::wrap_plots(biplots, ncol = 3)
  row2 <- patchwork::wrap_plots(boxes, ncol = 6)
  bars <- lapply(c("PC1", "PC2"), function(pc) cancer_consensus_plot(analysis$consensus[[pc]], pc))
  row3 <- patchwork::wrap_plots(bars, ncol = 2)
  composite <- patchwork::wrap_plots(
    patchwork::wrap_elements(row1) + ggplot2::labs(tag = "A"),
    patchwork::wrap_elements(row2) + ggplot2::labs(tag = "B"),
    patchwork::wrap_elements(row3) + ggplot2::labs(tag = "C"),
    ncol = 1, heights = c(1, 0.95, 1.8))
  files <- c(breast_composite = cancer_save_pdf(composite,
    file.path(output_dir, "breast_cross_cohort_pca_and_stratification.pdf"), 12, 14))
  heatmaps <- clusters <- centered <- list()
  for (cohort in c("Yan", "Zhang", "Bassez", "Pelka")) {
    heat <- cancer_heatmap(analysis$clr[[cohort]], analysis$metadata[[cohort]], cohort)
    heatmaps[[cohort]] <- heat$plot
    clusters[[cohort]] <- heat$clusters
    centered[[cohort]] <- heat$matrix
    dimensions <- switch(cohort, Yan = c(11, 8.5), Zhang = c(7, 8.5), Bassez = c(11, 9), Pelka = c(10, 8))
    files[paste0(tolower(cohort), "_heatmap")] <- cancer_save_pdf(heat$plot$gtable,
      file.path(output_dir, paste0(tolower(cohort), "_hitme_shared_heatmap.pdf")), dimensions[1], dimensions[2])
  }
  appendix <- patchwork::wrap_plots(patchwork::wrap_elements(heatmaps$Yan$gtable),
    patchwork::wrap_elements(heatmaps$Zhang$gtable), ncol = 2, widths = c(1.35, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
  files["appendix_heatmaps"] <- cancer_save_pdf(appendix,
    file.path(output_dir, "appendix_yan_zhang_hr_heatmaps.pdf"), 18, 9)
  correlation <- cancer_correlation(analysis$pca)
  files["correlation"] <- cancer_save_pdf(correlation$plot,
    file.path(output_dir, "cross_cohort_all_datasets_pc1_pc2_clustered_correlation_matrix_heatmap.pdf"), 6, 5)
  list(plots = list(breast_composite = composite, biplots = biplots, boxplots = boxes,
      consensus = bars, heatmaps = heatmaps, appendix_heatmaps = appendix, correlation = correlation$plot),
    files = files, correlation = correlation[names(correlation) != "plot"], clusters = clusters,
    heatmap_matrices = centered)
}
