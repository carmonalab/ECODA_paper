# ============================================================
# MATH UTILITIES
# ============================================================

# Centered log-ratio transformation
clr <- function(df) {
  geometric_mean <- apply(df, 1, function(row) exp(mean(log(row))))
  clr_df <- apply(df, 2, function(row) log(row) - log(geometric_mean)) %>%
    as.data.frame()

  return(clr_df)
}

# Calculate coefficient of variation
cv <- function(x) {
  sd_x <- sd(x) # Standard deviation
  mean_x <- mean(x) # Mean
  cv_value <- abs(sd_x / mean_x) * 100 # Coefficient of variation (%)
  return(cv_value)
}

# Convert count matrix to percentages
calc_perc_df <- function(df) {
  df <- t(apply(df, 1, function(row) (row / sum(row)) * 100)) %>%
    as.data.frame()
  return(df)
}

# Zero imputation for compositional data
impute_zeros <- function(
  df,
  clr_zero_impute_method = c(
    "percentage_zeros",
    "percentage_all",
    "counts_zeros",
    "counts_all"
  ),
  clr_zero_impute_num = 1
) {
  if (
    !clr_zero_impute_method %in%
      c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all")
  ) {
    stop("clr_zero_impute_method not found")
  }

  # Apply specified zero imputation method
  if (clr_zero_impute_method == "percentage_zeros") {
    for (row in 1:nrow(df)) {
      df[row, ][df[row, ] == 0] <- sum(df[row, ]) / 100 * clr_zero_impute_num
    }
  } else if (clr_zero_impute_method == "percentage_all") {
    for (row in 1:nrow(df)) {
      df[row, ] <- df[row, ] + sum(df[row, ]) / 100 * clr_zero_impute_num
    }
  } else if (clr_zero_impute_method == "counts_zeros") {
    # Impute zeros by replacing them with a small non-zero value (1 in this case)
    df[df == 0] <- clr_zero_impute_num
  } else if (clr_zero_impute_method == "counts_all") {
    # Impute by adding a fixed count to all values
    df <- df + clr_zero_impute_num
  }

  return(df)
}

# Min-max-scaling score metrics
min_max <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  if (max(x, na.rm = T) == min(x, na.rm = T)) {
    return(x)
  }
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Z-score transformation for distance matrices
zscore_transform <- function(dist_mat) {
  mu <- mean(dist_mat)
  sigma <- sd(dist_mat)
  z_score_matrix <- (dist_mat - mu) / sigma

  return(z_score_matrix)
}

# Global quantile normalization to Gaussian distribution
global_quantile_norm_gaussian <- function(dist_mat) {
  # 1. Flatten the matrix to a single vector
  v <- as.vector(dist_mat)

  # 2. Calculate ranks (handling ties by averaging)
  r <- rank(v, ties.method = "average")

  # 3. Convert ranks to probabilities (Blom's method to avoid Inf)
  probs <- (r - 0.5) / length(v)

  # 4. Apply Inverse Normal CDF (qnorm)
  v_norm <- qnorm(probs)

  # 5. Reshape back to original matrix dimensions
  mat_norm <- matrix(v_norm, nrow = nrow(dist_mat), ncol = ncol(dist_mat))

  # Restore names if they existed
  dimnames(mat_norm) <- dimnames(dist_mat)

  return(mat_norm)
}