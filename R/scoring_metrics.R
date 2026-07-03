# ============================================================
# SCORING FUNCTIONS
# ============================================================

# Calculate separation scores (ANOSIM, modularity, clustering, silhouette)
calc_sep_score <- function(dist_mat, labels, knn_k = NULL) {
  # sil_score <- round(calc_sil(df, labels), 3)
  mod_knnsqrtn_score <- unlist(round(
    calc_modularity(dist_mat, labels, knn_k),
    3
  ))
  mod_knn3_score <- unlist(round(calc_modularity(dist_mat, labels, 3), 3))
  mod_knn6_score <- unlist(round(calc_modularity(dist_mat, labels, 6), 3))
  mod_knn9_score <- unlist(round(calc_modularity(dist_mat, labels, 9), 3))
  anosim_score <- vegan::anosim(
    x = dist_mat,
    grouping = labels,
    distance = "euclidean"
  )[["statistic"]]
  cluster_score <- clust_eval(dist_mat, labels)
  sil_score <- calc_sil(dist_mat, labels)

  res <- list(
    # sil_score = sil_score,
    mod_knnsqrtn_score = mod_knnsqrtn_score,
    mod_knn3_score = mod_knn3_score,
    mod_knn6_score = mod_knn6_score,
    mod_knn9_score = mod_knn9_score,
    anosim_score = anosim_score,
    cluster_score = cluster_score,
    sil_score = sil_score
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
  digits = 3,
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
    return(round(mean(unlist(results)), digits))
  } else {
    results[["hclust_accuracy"]] <- round(results[["hclust_accuracy"]], digits)
    results[["pamclust_accuracy"]] <- round(
      results[["pamclust_accuracy"]],
      digits
    )
    return(results)
  }
}