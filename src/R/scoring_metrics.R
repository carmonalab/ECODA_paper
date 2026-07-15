# ============================================================
# SCORING FUNCTIONS
# ============================================================

# Calculate separation scores (ANOSIM, modularity, clustering, silhouette)
calc_sep_score <- function(dist_mat, labels, knn_k = NULL) {
  num_labels <- as.numeric(as.factor(labels))

  sil_score <- calc_sil(dist_mat, num_labels)
  mod_knnsqrtn_score <- calc_modularity(dist_mat, num_labels, knn_k)
  mod_knn3_score <- calc_modularity(dist_mat, num_labels, 3)
  mod_knn6_score <- calc_modularity(dist_mat, num_labels, 6)
  mod_knn9_score <- calc_modularity(dist_mat, num_labels, 9)
  anosim_score <- vegan::anosim(
    x = dist_mat,
    grouping = num_labels,
    distance = "euclidean"
  )[["statistic"]]
  cluster_score <- clust_eval(dist_mat, num_labels)
  lisi_score <- calc_lisi(num_labels, dist_mat = dist_mat)

  res <- list(
    sil_score = sil_score,
    mod_knnsqrtn_score = mod_knnsqrtn_score,
    mod_knn3_score = mod_knn3_score,
    mod_knn6_score = mod_knn6_score,
    mod_knn9_score = mod_knn9_score,
    anosim_score = anosim_score,
    cluster_score = cluster_score,
    lisi_score = lisi_score
  )

  return(res)
}


# Calculate average silhouette width
calc_sil <- function(dist_mat, labels) {
  sils <- cluster::silhouette(
    x = as.numeric(factor(labels)),
    dist = dist_mat
  ) %>%
    as.data.frame()
  score <- mean(sils[["sil_width"]])
  return(score)
}


# Calculate modularity score (adjusted for number of groups)
calc_modularity <- function(dist_mat, labels, knn_k = NULL) {
  ngroups <- length(unique(labels))

  if (is.null(knn_k)) {
    knn_k <- max(3, round(sqrt(attr(dist_mat, "Size"))))
  }

  # Create a graph object
  knn <- compute_KNN_from_dist(dist_mat, knn_k)
  g <- compute_snn_graph(knn)
  # Compute modularity
  modularity_score <- igraph::modularity(
    g,
    membership = as.numeric(factor(labels))
  )

  # NOTE:
  # Maximum modularity depends on the number of groups: max(mod) = 1 - 1 / (number of groups)
  # see Brandes, Ulrik, et al. On finding graph clusterings with maximum modularity.

  # Adjust modularity score for number of groups
  maximum_modularity_score <- 1 - (1 / ngroups)
  adjusted_modularity_score <- modularity_score / maximum_modularity_score

  return(adjusted_modularity_score)
}

# Compute SNN graph from KNN matrix
compute_snn_graph <- function(knn) {
  n <- nrow(knn)
  k <- ncol(knn)

  # 1. Create a binary sparse matrix (A)
  i_idx <- rep(1:n, each = k)
  j_idx <- as.vector(t(knn))

  adj_bin <- Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = 1,
    dims = c(n, n)
  )

  # 2. Compute SNN weights for ALL pairs using Matrix Multiplication
  snn_matrix <- adj_bin %*% t(adj_bin)
  Matrix::diag(snn_matrix) <- 0

  # 3. Convert to igraph
  g <- igraph::graph_from_adjacency_matrix(
    snn_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  return(g)
}

# Compute KNN from distance matrix
compute_KNN_from_dist <- function(dist_mat, knn_k) {
  # dist_mat should be a square matrix or 'dist' object
  dist_mat <- as.matrix(dist_mat)
  n <- nrow(dist_mat)
  knn <- matrix(0, nrow = n, ncol = knn_k)

  for (i in 1:n) {
    # Sort distances in row i, get indices of the 2nd to (k+1)th closest
    knn[i, ] <- order(dist_mat[i, ])[2:(knn_k + 1)]
  }
  return(knn)
}


# Cluster samples and compare to original annotation
clust_eval <- function(
  dist_mat,
  labels,
  nclusts = NULL,
  return_mean = TRUE
) {
  results <- list()

  if (is.null(nclusts)) {
    nclusts <- length(unique(labels))
  }

  # Perform hierarchical clustering
  hc <- stats::hclust(dist_mat, method = "ward.D2")
  clust_labels <- stats::cutree(hc, k = nclusts)
  results[["hclust_accuracy"]] <- mclust::adjustedRandIndex(
    as.numeric(as.factor(labels)),
    clust_labels
  )

  # Perform PAM clustering
  clust_labels <- cluster::pam(dist_mat, k = nclusts)$cluster
  results[["pamclust_accuracy"]] <- mclust::adjustedRandIndex(
    as.numeric(as.factor(labels)),
    clust_labels
  )

  if (return_mean) {
    return(mean(unlist(results)))
  } else {
    return(results)
  }
}


calc_lisi <- function(labels, features = NULL, dist_mat = NULL) {
  # 1. Input Validation & Diagnostics
  if (is.null(features) && is.null(dist_mat)) {
    stop("calc_lisi: You must provide either 'features' or 'dist_mat'.")
  }

  N <- length(unique(labels))
  if (N <= 1) {
    warning("calc_lisi: Only one unique label found. Returning NA.")
    return(NA_real_)
  }

  # 2. Determine number of samples (n)
  n_samples <- if (!is.null(features)) {
    nrow(features)
  } else {
    attr(as.dist(dist_mat), "Size")
  }

  # 3. Dynamically set perplexity for small datasets
  # Default is 30. For small n, we cap perplexity at 10% of samples (minimum of 2 to keep it valid)
  dynamic_perplexity <- max(2, min(30, floor(0.10 * n_samples)))

  # 4. Handle Distance Matrix Input (using lossless MDS)
  if (is.null(features)) {
    if (!inherits(dist_mat, "dist")) {
      dist_mat <- as.dist(dist_mat)
    }

    # Perfectly represent the distance matrix using n - 1 dimensions
    k_dims <- n_samples - 1

    # Perform Classical MDS (Lossless)
    features <- stats::cmdscale(dist_mat, k = k_dims)
    rownames(features) <- paste0("Sample_", 1:nrow(features))
  }

  # 5. Format metadata for compute_lisi
  metadata <- data.frame(
    label = as.character(labels),
    row.names = rownames(features)
  )

  # 6. Compute standard LISI
  # We pass our adjusted dynamic_perplexity to prevent solver errors
  standard_lisi <- thisutils::compute_lisi(
    X = features,
    meta_data = metadata,
    label_colnames = "label",
    perplexity = dynamic_perplexity
  )

  # 7. Transform to 0-1 Separation Score: (N - LISI) / (N - 1)
  # 1 = perfect separation, 0 = perfect mixing
  # Transform to 0-1 Separation Score: (N - LISI) / (N - 1)
  raw_scores <- (N - standard_lisi$label) / (N - 1)
  
  # OPTIMIZATION: Enforce physical boundaries [0, 1] to handle numerical float drift
  separation_scores <- pmax(0, pmin(1, raw_scores))

  return(mean(separation_scores, na.rm = TRUE))
}
