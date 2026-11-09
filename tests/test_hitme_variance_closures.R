# Actual composition builder and CLR/variance/count producer; no RNA/model fits.
root <- normalizePath('.')
Sys.setenv(PROJECT_ROOT = root)
suppressPackageStartupMessages(library(dplyr))
env <- new.env(parent = globalenv())
for (file in c('src/utils/math_utils.R', 'src/utils/hvcs.R', 'src/utils/seurat_utils.R',
               'src/5_run_benchmark_methods/benchmark_methods_r.R',
               'src/5_run_benchmark_methods/benchmark_hpc_utils.R',
               'src/5_run_benchmark_methods/benchmark_pipeline.R')) {
  sys.source(file.path(root, file), envir = env)
}
# Isolate metrics and non-composition methods; real counts/CLR/HVC math is used.
env$create_result_bundle <- function(feat_mat, labels, dist_mat = stats::dist(feat_mat), extra = list()) {
  c(list(feat_mat = feat_mat, labels = labels, dist_mat = dist_mat, scores = list(sil_score = 0)), extra)
}
real_coda <- env$process_coda_fig
env$process_coda_fig <- function(seurat, labels, ..., ct_col, ECODA_top_varexp_hvct = NULL) {
  result <- real_coda(seurat, labels, ..., ct_col = ct_col, ECODA_top_varexp_hvct = ECODA_top_varexp_hvct)
  result$observed_ct_col <- ct_col
  result$observed_variance_threshold <- ECODA_top_varexp_hvct
  result
}
dummy <- list(feat_mat = matrix(1, 6, 2), scores = list(sil_score = 0))
env$process_avg_pca_embedding_fig <- env$process_deconv_fig <- env$process_gloprop_fig <- function(...) dummy
env$exec_time <- function(expr) { force(expr); as.difftime(1, units = 'secs') }
env$peak_rss_gb <- function() 0
env$log_exec_row <- function(...) invisible(NULL)
observations <- lapply(1:6, function(i) {
  fine <- rep(c('T_resting', 'T_cycle', 'B_resting', 'B_cycle', 'M_resting', 'M_cycle'),
              c(5+i, 15+i, 9+2*i, 22-i, 14-i, 15-2*i))
  data.frame(Sample = paste0('s', i), biology = if (i <= 3) 'A' else 'B',
             author_hr = rep(c('author_T', 'author_B', 'author_M', 'author_other'), length.out = length(fine)),
             layer2 = sub('_(resting|cycle)$', '', fine), layer3 = fine)
})
obs <- do.call(rbind, observations)
obs$RNA_snn_res.0.1 <- obs$author_hr
labels <- setNames(factor(rep(c('A', 'B'), each = 3)), paste0('s', 1:6))
metadata <- unique(obs[, c('Sample', 'biology')])
thresholds <- seq(0, 0.9, 0.1)
output <- tempfile('ecoda-hitme-closure-'); dir.create(output)
result <- suppressWarnings(env$run_composition_methods_hpc(
  labels, metadata, NULL, list(pb = matrix(1, 2, 6)), obs, 'biology',
  ct_col_high_res = 'author_hr', results_dir = output, ds = 'Tiny',
  factors_test = 2L, seurat_res = 0.1,
  ECODA_top_varexp_hvct = thresholds, not_suitable_for_auto_annotation = 'scatomic'
))
for (layer in c('layer2', 'layer3')) {
  plain <- result[[paste0('ECODA_HiTME_HR_', layer)]]
  stopifnot(identical(plain$observed_ct_col, layer))
  for (v in thresholds) {
    key <- paste0('ECODA_HiTME_HR_', layer, '_top_varexp', v)
    actual <- result[[key]]
    if (!identical(actual$observed_ct_col, layer)) {
      stop(key, ' used ', actual$observed_ct_col, ' instead of ', layer)
    }
    stopifnot(identical(actual$observed_variance_threshold, v))
    expected <- real_coda(NULL, labels, ct_col = layer, obs = obs, ECODA_top_varexp_hvct = v)
    for (field in c('counts', 'feat_mat', 'dist_mat', 'labels')) stopifnot(identical(actual[[field]], expected[[field]]))
    # Check the real persisted payload, with checksum validation before reading.
    saved <- get('read_rds_checked', envir = env, inherits = TRUE)(file.path(output, paste0('Tiny_', key, '.rds')))
    stopifnot(identical(saved$counts, expected$counts), identical(saved$observed_ct_col, layer))
  }
}
for (v in thresholds) {
  actual <- result[[paste0('ECODA_authors_HR_top_varexp', v)]]
  expected <- real_coda(NULL, labels, ct_col = 'author_hr', obs = obs, ECODA_top_varexp_hvct = v)
  stopifnot(identical(actual$counts, expected$counts), identical(actual$feat_mat, expected$feat_mat))
}
stopifnot(!any(grepl('_(resting|cycle)$', colnames(result$ECODA_HiTME_HR_layer2_top_varexp0$counts))),
          any(grepl('_(resting|cycle)$', colnames(result$ECODA_HiTME_HR_layer3_top_varexp0$counts))))
unlink(output, recursive = TRUE)
cat('PASS: deferred HiTME L2/L3 variance closures retain exact layer and threshold; plain/author results unchanged\n')
