# Compact exported UCell artifacts only; no H5AD loading or benchmark mutation.
stephenson_ucell_scores <- c(
  "IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell", "cellCycle.G2M_UCell"
)
stephenson_ucell_annotations <- c(
  "authors_lowres", "authors_highres", "hitme_layer1", "hitme_layer2"
)

load_stephenson_ucell <- function(root) {
  if (!dir.exists(root)) stop("UCell artifact directory does not exist: ", root)
  manifest <- jsonlite::fromJSON(file.path(root, "manifest.json"), simplifyVector = FALSE)
  ids <- unlist(manifest$sample_ids, use.names = FALSE)
  valid_ids <- function(x) is.character(x) && length(x) > 0L &&
    !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
  if (!identical(manifest$schema_version, 1L) ||
      !identical(manifest$dataset, "Stephenson") ||
      !identical(manifest$view, "benchmark_analysis") ||
      !identical(manifest$label_column, "Status") || !valid_ids(ids) ||
      !identical(ids, sort(ids, method = "radix")) ||
      !setequal(names(manifest$annotations), stephenson_ucell_annotations) ||
      !identical(unlist(manifest$scores, use.names = FALSE), stephenson_ucell_scores)) {
    stop("Invalid Stephenson UCell manifest schema/IDs")
  }
  keys <- unlist(lapply(stephenson_ucell_annotations, function(a) {
    paste(a, stephenson_ucell_scores, sep = "__")
  }), use.names = FALSE)
  expected <- c(paste0(keys, ".feather"),
                paste0(stephenson_ucell_annotations, "__n_cells.feather"), "samples.feather")
  if (!setequal(names(manifest$files), expected)) stop("Manifest file inventory mismatch")
  # Authenticate every compact input before parsing any Feather data.
  for (name in expected) {
    checksum <- manifest$files[[name]]
    path <- file.path(root, name)
    if (!is.character(checksum) || length(checksum) != 1L ||
        !grepl("^[0-9a-f]{64}$", checksum) || !file.exists(path) ||
        !identical(digest::digest(path, algo = "sha256", file = TRUE), checksum)) {
      stop("UCell checksum mismatch: ", name)
    }
  }
  read_table <- function(name) as.data.frame(arrow::read_feather(file.path(root, name)))
  samples <- read_table("samples.feather")
  if (!identical(names(samples), c("Sample", "Status")) ||
      !identical(samples$Sample, ids) || anyNA(samples$Status) ||
      any(!nzchar(trimws(as.character(samples$Status))))) {
    stop("Invalid UCell sample labels/IDs")
  }
  matrices <- counts <- list()
  coverage <- list()
  for (annotation in stephenson_ucell_annotations) {
    count_table <- read_table(paste0(annotation, "__n_cells.feather"))
    types <- unlist(manifest$cell_type_columns[[annotation]], use.names = FALSE)
    if (!valid_ids(types) || !identical(types, sort(types, method = "radix"))) {
      stop("Invalid manifest cell-type axes: ", annotation)
    }
    as_feature_matrix <- function(tab, name) {
      if (!identical(names(tab), c("Sample", types)) || !identical(tab$Sample, ids) ||
          !all(vapply(tab[-1L], is.numeric, logical(1)))) {
        stop("UCell matrix axes/types mismatch: ", name)
      }
      mat <- as.matrix(tab[-1L])
      rownames(mat) <- ids
      if (any(is.infinite(mat)) || any(is.nan(mat))) stop("Invalid nonfinite values: ", name)
      mat
    }
    n <- as_feature_matrix(count_table, annotation)
    if (anyNA(n) || any(n < 0 | n != floor(n))) stop("Invalid cell counts: ", annotation)
    counts[[annotation]] <- n
    excluded <- manifest$excluded_annotation_cells[[annotation]]
    if (!is.numeric(excluded) || length(excluded) != 1L || !is.finite(excluded) ||
        excluded < 0 || sum(n) + excluded != manifest$n_cells) {
      stop("Cell counts disagree with manifest: ", annotation)
    }
    for (score in stephenson_ucell_scores) {
      key <- paste(annotation, score, sep = "__")
      mat <- as_feature_matrix(read_table(paste0(key, ".feather")), key)
      if (any(!is.na(mat[n == 0]))) stop("Absent activity groups must be NA: ", key)
      matrices[[key]] <- mat
      selected <- prepare_stephenson_ucell_features(mat, require_pca = FALSE)
      coverage[[key]] <- data.frame(
        annotation = annotation, score = score, n_samples = nrow(mat),
        samples_with_annotation = sum(rowSums(n) > 0),
        observed_types = ncol(mat), retained_types = ncol(selected$features),
        observed_sample_type_groups = sum(n > 0),
        finite_activity_groups = sum(is.finite(mat)),
        excluded_annotation_cells = excluded,
        dropped_incomplete = paste(selected$dropped_incomplete, collapse = "; "),
        dropped_constant = paste(selected$dropped_constant, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }
  }
  labels <- factor(samples$Status, levels = sort(unique(samples$Status)))
  names(labels) <- ids
  list(manifest = manifest, matrices = matrices, counts = counts, samples = samples,
       labels = labels, coverage = do.call(rbind, coverage))
}

prepare_stephenson_ucell_features <- function(mat, require_pca = TRUE) {
  if (!is.matrix(mat) || !is.numeric(mat) || is.null(colnames(mat)) ||
      anyDuplicated(colnames(mat))) stop("Expected a numeric matrix with unique feature names")
  complete <- colSums(!is.finite(mat)) == 0L
  constant <- rep(FALSE, ncol(mat))
  constant[complete] <- vapply(seq_len(ncol(mat))[complete], function(j) {
    length(unique(mat[, j])) < 2L
  }, logical(1))
  features <- mat[, complete & !constant, drop = FALSE]
  if (require_pca && (nrow(features) < 4L || ncol(features) < 2L)) {
    stop("PCA and k=3 comparison require >=4 samples and >=2 complete nonconstant features")
  }
  list(features = features, dropped_incomplete = colnames(mat)[!complete],
       dropped_constant = colnames(mat)[constant])
}

align_stephenson_bundle <- function(bundle, sample_ids, labels) {
  mat <- as.matrix(bundle$feat_mat)
  ids <- rownames(mat)
  label_ids <- names(bundle$labels)
  valid_ids <- function(x) !is.null(x) && !anyNA(x) &&
    all(nzchar(trimws(x))) && !anyDuplicated(x)
  if (!valid_ids(sample_ids) || !valid_ids(ids) || !valid_ids(label_ids) ||
      !setequal(ids, sample_ids) || !setequal(label_ids, sample_ids) ||
      !identical(names(labels), sample_ids)) {
    stop("Benchmark sample IDs must match the complete UCell sample universe exactly")
  }
  aligned_labels <- bundle$labels[match(sample_ids, label_ids)]
  if (anyNA(aligned_labels) || !identical(as.character(aligned_labels), as.character(labels))) {
    stop("Benchmark and UCell labels disagree for aligned sample IDs")
  }
  mat <- mat[match(sample_ids, ids), , drop = FALSE]
  if (!is.numeric(mat) || any(!is.finite(mat))) stop("Benchmark features must be finite numeric values")
  prepare_stephenson_ucell_features(mat)
}

validate_stephenson_n_loadings <- function(n_loadings) {
  if (!is.numeric(n_loadings) || length(n_loadings) != 1L ||
      !is.finite(n_loadings) || n_loadings < 1 || n_loadings != floor(n_loadings)) {
    stop("n_loadings must be a positive integer")
  }
  invisible(n_loadings)
}

stephenson_pca_comparison <- function(mat, labels, title, n_loadings, seed) {
  distances <- stats::dist(mat)
  set.seed(seed)
  anosim <- unname(vegan::anosim(distances, labels, parallel = 1)$statistic)
  set.seed(seed)
  ari <- clust_eval(distances, labels)
  set.seed(seed)
  modularity <- calc_modularity(distances, labels, knn_k = 3)
  # Plot and report the same full-feature statistics, never PC1/2 distances.
  set.seed(seed)
  panel <- plot_pca(
    mat, labels = labels, scale. = FALSE, knn_k = 3, parallel = 1,
    coord_equal = FALSE, title = title, invisible = "quali",
    n_ct_show = min(n_loadings, ncol(mat)), repel = TRUE, labelsize = 3.5,
    cluster_score = FALSE, mod_score = FALSE, anosim_score = TRUE
  ) + ggplot2::scale_color_manual(
    values = c(Covid = "#F8766D", Healthy = "#00BFC4"),
    limits = levels(labels), drop = FALSE
  ) + ggplot2::labs(color = "Status", shape = "Status")
  # fviz_pca fixes variable colors independently of the sample Status scale.
  for (layer in panel$layers) {
    if (inherits(layer$geom, c("GeomSegment", "GeomText", "GeomTextRepel"))) {
      layer$aes_params$colour <- "black"
    }
  }
  list(plot = panel, ANOSIM = anosim, modularity_k3 = modularity, ARI = ari)
}

compare_stephenson_ucell <- function(
  artifacts, ecoda_bundles, score = "IFN_UCell",
  annotations = c("authors_lowres", "authors_highres"), seed = 20261022L,
  n_loadings = 8L
) {
  validate_stephenson_n_loadings(n_loadings)
  methods <- c(authors_lowres = "ECODA_authors_LR", authors_highres = "ECODA_authors_HR")
  annotation_titles <- c(
    authors_lowres = "Authors LR", authors_highres = "Authors HR",
    hitme_layer1 = "HiTME layer 1", hitme_layer2 = "HiTME layer 2"
  )
  if (!score %in% stephenson_ucell_scores || !length(annotations) ||
      anyDuplicated(annotations) || !all(annotations %in% stephenson_ucell_annotations)) {
    stop("Unknown or repeated UCell annotation/score selection")
  }
  ids <- artifacts$samples$Sample
  labels <- artifacts$labels
  panels <- features <- tables <- list()
  for (annotation in annotations) {
    key <- paste(annotation, score, sep = "__")
    method <- unname(methods[annotation])
    selected <- list(UCell = prepare_stephenson_ucell_features(artifacts$matrices[[key]]))
    if (!is.na(method)) {
      if (is.null(ecoda_bundles[[method]])) stop("Missing matching ECODA bundle: ", method)
      selected$ECODA <- align_stephenson_bundle(ecoda_bundles[[method]], ids, labels)
    }
    for (representation in names(selected)) {
      item <- selected[[representation]]
      mat <- item$features
      name <- paste(annotation, representation, sep = "__")
      features[[name]] <- mat
      comparison <- stephenson_pca_comparison(
        mat, labels,
        title = paste0(annotation_titles[[annotation]], ": ",
                       if (representation == "UCell") gsub("_", " ", score, fixed = TRUE) else "ECODA"),
        n_loadings = n_loadings, seed = seed
      )
      panels[[name]] <- comparison$plot
      row <- artifacts$coverage[artifacts$coverage$annotation == annotation &
                                  artifacts$coverage$score == score, , drop = FALSE]
      row$representation <- representation
      row$comparison_samples <- nrow(mat)
      row$comparison_features <- ncol(mat)
      row$comparison_dropped_incomplete <- paste(item$dropped_incomplete, collapse = "; ")
      row$comparison_dropped_constant <- paste(item$dropped_constant, collapse = "; ")
      row$ANOSIM <- comparison$ANOSIM
      row$modularity_k3 <- comparison$modularity_k3
      row$ARI <- comparison$ARI
      tables[[name]] <- row
    }
  }
  list(plot = patchwork::wrap_plots(panels, ncol = 2, guides = "collect"),
       table = do.call(rbind, tables), features = features, labels = labels)
}

compare_stephenson_pseudobulk <- function(
  artifacts, bundle, n_loadings = 15L, seed = 20261022L
) {
  validate_stephenson_n_loadings(n_loadings)
  labels <- artifacts$labels
  selected <- align_stephenson_bundle(bundle, artifacts$samples$Sample, labels)
  excluded_genes <- colnames(selected$features)[grepl("^IG[HKL]", colnames(selected$features))]
  mat <- selected$features[, !colnames(selected$features) %in% excluded_genes, drop = FALSE]
  if (ncol(mat) < 2L) stop("IG exclusion leaves fewer than two nonconstant features for PCA")
  comparison <- stephenson_pca_comparison(
    mat, labels, title = "Pseudobulk HVG2000: IG excluded",
    n_loadings = n_loadings, seed = seed
  )
  list(
    plot = comparison$plot,
    table = data.frame(
      method = "Pseudobulk_hvg2000_IGexcluded", n_samples = nrow(mat), n_features = ncol(mat),
      n_excluded_ig = length(excluded_genes),
      ANOSIM = comparison$ANOSIM, modularity_k3 = comparison$modularity_k3,
      ARI = comparison$ARI
    ),
    features = mat, labels = labels, excluded_genes = excluded_genes
  )
}

# CD14 gene expression is normalized/selected upstream without Status or IG exclusion.
load_stephenson_cd14_pseudobulk <- function(root) {
  if (!dir.exists(root)) stop("CD14 artifact directory does not exist: ", root)
  manifest <- jsonlite::fromJSON(file.path(root, "manifest.json"), simplifyVector = FALSE)
  ids <- unlist(manifest$sample_ids, use.names = FALSE)
  norm <- manifest$normalization
  if (!identical(manifest$schema_version, 1L) ||
      !identical(manifest$dataset, "Stephenson") ||
      !identical(manifest$view, "benchmark_analysis") ||
      !identical(manifest$annotation_column, "initial_clustering") ||
      !identical(manifest$cell_type, "CD14") ||
      !is.character(ids) || length(ids) != 55L || anyNA(ids) ||
      any(!nzchar(trimws(ids))) || anyDuplicated(ids) ||
      !identical(ids, sort(ids, method = "radix")) ||
      !identical(manifest$n_features, 2000L) ||
      !is.list(norm) ||
      !setequal(names(norm), c("blind", "batch_col", "correct_batch", "n_hvg", "black_list")) ||
      !identical(norm$blind, TRUE) || !is.null(norm$batch_col) ||
      !identical(norm$correct_batch, FALSE) || !identical(norm$n_hvg, 2000L) ||
      !identical(norm$black_list, "none")) {
    stop("Invalid Stephenson CD14 manifest identity/normalization/IDs")
  }
  expected <- c("expression.feather", "samples.feather")
  if (!setequal(names(manifest$files), expected)) stop("CD14 manifest file inventory mismatch")
  for (name in expected) {
    checksum <- manifest$files[[name]]
    path <- file.path(root, name)
    if (!is.character(checksum) || length(checksum) != 1L ||
        !grepl("^[0-9a-f]{64}$", checksum) || !file.exists(path) ||
        !identical(digest::digest(path, algo = "sha256", file = TRUE), checksum)) {
      stop("CD14 checksum mismatch: ", name)
    }
  }
  samples <- as.data.frame(arrow::read_feather(file.path(root, "samples.feather")))
  if (!identical(names(samples), c("Sample", "Status", "n_cells")) ||
      !identical(samples$Sample, ids) || !is.character(samples$Status) ||
      anyNA(samples$Status) || !setequal(samples$Status, c("Covid", "Healthy")) ||
      !is.numeric(samples$n_cells) || any(!is.finite(samples$n_cells)) ||
      any(samples$n_cells < 5 | samples$n_cells != floor(samples$n_cells))) {
    stop("Invalid CD14 sample labels/IDs/cell counts")
  }
  total <- manifest$n_selected_cells
  if (!is.numeric(total) || length(total) != 1L || !is.finite(total) ||
      total != sum(samples$n_cells)) stop("CD14 cell counts disagree with manifest")
  expression <- as.data.frame(arrow::read_feather(file.path(root, "expression.feather")))
  genes <- names(expression)[-1L]
  if (!identical(names(expression)[1L], "Sample") ||
      !identical(expression$Sample, ids) || length(genes) != manifest$n_features ||
      anyNA(genes) || any(!nzchar(trimws(genes))) || anyDuplicated(genes) ||
      "Sample" %in% genes || !all(vapply(expression[-1L], is.numeric, logical(1)))) {
    stop("CD14 expression sample/feature axes or types mismatch")
  }
  mat <- as.matrix(expression[-1L])
  rownames(mat) <- ids
  if (any(!is.finite(mat))) stop("CD14 expression must be finite; no zero imputation")
  if (any(vapply(seq_len(ncol(mat)), function(j) length(unique(mat[, j])) < 2L,
                 logical(1)))) stop("CD14 HVGs must be nonconstant")
  labels <- factor(samples$Status, levels = c("Covid", "Healthy"))
  names(labels) <- ids
  list(feat_mat = mat, labels = labels, metadata = samples, manifest = manifest)
}
