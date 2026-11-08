# Focused regression for the pre-fit human five-group pseudobulk blacklist.
root <- normalizePath(".")
env <- new.env(parent = globalenv())
sys.source(file.path(root, "src/utils/pseudobulk.R"), envir = env)

assert_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

groups <- readRDS(file.path(root, "aux/genes.blocklist.rds"))$Hs
blocked <- unique(vapply(groups, function(x) unname(x[[1L]]), character(1)))
genes <- c(blocked, "XIST", "SRY", "retained_gene")
counts <- matrix(seq_len(length(genes) * 3L), nrow = length(genes),
                 dimnames = list(genes, c("s1", "s2", "s3")))
metadata <- data.frame(Sample = c("s3", "s1", "s2"),
                       row.names = c("s3", "s1", "s2"))
fitted_genes <- NULL
env$.fit_pseudobulk_deseq2_validated <- function(counts, ...) {
  fitted_genes <<- rownames(counts)
  list(norm_matrix = as.matrix(counts), variance_order = rownames(counts))
}
actual <- env$get_pb_deseq2_from_counts(
  counts, metadata, n_hvg = 2000L,
  black_list = "default_without_sex_genes"
)
stopifnot(!any(blocked %in% fitted_genes),
          identical(fitted_genes, c("XIST", "SRY", "retained_gene")),
          identical(rownames(actual), metadata$Sample))
cat("PASS: five-group blacklist is applied before pseudobulk fitting\n")

resolved <- env$pseudobulk_blacklist()
stopifnot(identical(names(resolved$groups), names(groups)),
          identical(resolved$policy$source_md5,
                    unname(tools::md5sum(file.path(root, "aux/genes.blocklist.rds")))))
for (flag in c("human_aux_five_groups", "default", "default_without_sex_genes")) {
  stopifnot(identical(
    env$get_pb_deseq2_from_counts(counts, metadata, black_list = flag), actual
  ))
}
stopifnot(identical(
  env$get_pb_deseq2_from_counts(counts, metadata),
  t(counts[, metadata$Sample, drop = FALSE])
))

# Stable IDs, their versions, and AnnData's duplicate-symbol suffixes retain
# the same membership. Sex-linked genes outside these groups are retained;
# RPS4Y1 is itself in the chosen Ribo group and receives no exemption.
stable <- names(groups[[1L]])[[1L]]
variants <- c(blocked[[1L]], paste0(blocked[[1L]], "-1"), stable,
              paste0(stable, ".7"), paste0(stable, ".7-1"),
              "RPS4Y1", "XIST", "SRY")
namespaced <- matrix(1, nrow = length(variants), ncol = 1L,
                     dimnames = list(variants, "s1"))
filtered <- env$filter_pseudobulk_blacklist(namespaced, "human_aux_five_groups")
stopifnot(identical(rownames(filtered$counts), c("XIST", "SRY")),
          length(filtered$audit$excluded_genes) == 6L,
          filtered$audit$matched_rows_per_group[["cellCycle.G1S"]] >= 5L)
orphan <- namespaced[c(2L, 7L), , drop = FALSE]
assert_error(env$filter_pseudobulk_blacklist(orphan, "human_aux_five_groups"),
             "without an unsuffixed feature")
assert_error(env$select_pseudobulk_deseq2(
  list(norm_matrix = counts, variance_order = rownames(counts)), 2L,
  "human_aux_five_groups"), "pre-fit")
assert_error(env$filter_pseudobulk_blacklist(counts[blocked, , drop = FALSE],
                                           "human_aux_five_groups"), "every input gene")

# Legacy Seurat callers reach the same pre-fit mask without attaching STACAS.
env$get_pb <- function(...) counts
env$get_metadata <- function(...) metadata
env$DESeq2.normalize <- function(count_mat, ...) {
  stopifnot(identical(rownames(count_mat), c("XIST", "SRY", "retained_gene")))
  count_mat
}
stopifnot(identical(env$get_pb_deseq2(NULL, black_list = "human_aux_five_groups"), actual))

local({
  prior <- Sys.getenv("ECODA_AUX_ROOT", unset = NA_character_)
  on.exit(if (is.na(prior)) Sys.unsetenv("ECODA_AUX_ROOT") else
    Sys.setenv(ECODA_AUX_ROOT = prior))
  temporary <- tempfile("pseudobulk-blacklist-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  Sys.setenv(ECODA_AUX_ROOT = temporary)
  assert_error(env$pseudobulk_blacklist(), "missing")
  saveRDS(list(Hs = groups[names(groups) != "TCR"]),
          file.path(temporary, "genes.blocklist.rds"))
  assert_error(env$pseudobulk_blacklist(), "five-group schema")
})

sys.source(file.path(root, "src/5_run_benchmark_methods/benchmark_pipeline.R"), envir = env)
sys.source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"), envir = env)
aggregations <- 0L
fits <- list()
env$load_h5ad_sample_aggregate <- function(..., metadata_columns) {
  stopifnot(identical(metadata_columns, "Sample"))
  aggregations <<- aggregations + 1L
  list(counts = counts, metadata = metadata, sample_ids = metadata$Sample)
}
env$fit_pseudobulk_deseq2 <- function(counts, ...) {
  fits[[length(fits) + 1L]] <<- rownames(counts)
  list(norm_matrix = as.matrix(counts), variance_order = rownames(counts))
}
env$exec_time <- function(expr) { force(expr); as.difftime(1, units = "secs") }
env$peak_rss_gb <- function() 0.25
input <- tempfile("pseudobulk-blacklist-input-")
writeLines("stubbed H5AD read", input)
prepare <- function(variants) env$prepare_pseudobulks_hpc(
  input, hvg_rank_genes = rownames(counts), variants = variants,
  cache_stem = "Toy", run_id = "blacklist-test"
)
ordinary <- prepare(c("hvg500", "hvg2000"))
stopifnot(aggregations == 1L, length(fits) == 1L,
          identical(fits[[1L]], rownames(counts)))
aggregations <- 0L; fits <- list()
blacklist <- prepare("hvg2000_bl")
bl <- blacklist$hvg2000_bl
stopifnot(aggregations == 1L, length(fits) == 1L,
          identical(fits[[1L]], c("XIST", "SRY", "retained_gene")),
          bl$aggregate_time_secs == 0, bl$shared_fit_time_secs == 0,
          bl$shared_time_secs == 0, bl$variant_time_secs == 2,
          bl$blacklist_aggregate_time_secs == 1,
          bl$timing_id != ordinary$hvg2000$timing_id,
          identical(bl$blacklist_policy, resolved$policy))
env$validate_pseudobulk_timing_record(bl, "hvg2000_bl")
aggregations <- 0L; fits <- list()
mixed <- prepare(c("hvg500", "hvg2000_bl", "hvg2000"))
stopifnot(aggregations == 2L, length(fits) == 2L,
          identical(mixed[c("hvg500", "hvg2000")], ordinary),
          identical(mixed$hvg2000_bl, bl))
rows <- list()
env$log_exec_row <- function(ds, method, time_secs, ..., mem_gb) {
  rows[[length(rows) + 1L]] <<- list(method = method, time_secs = time_secs)
}
env$emit_pseudobulk_timing_rows(mixed, "Toy", "stubbed.feather")
stopifnot(identical(vapply(rows, `[[`, "", "method"), c(
  "prepare_pseudobulk_shared", "prepare_pseudobulk_hvg500",
  "prepare_pseudobulk_hvg2000_bl", "prepare_pseudobulk_hvg2000"
)), rows[[3L]]$time_secs == 2)

cache <- tempfile("pseudobulk-blacklist-cache-", fileext = ".rds")
write_cache <- function(value) {
  saveRDS(value, cache)
  writeLines(c(paste0("MD5=", unname(tools::md5sum(cache))),
               paste0("SIZE=", file.info(cache)$size), paste0("PATH=", cache)),
             paste0(cache, ".md5"))
}
write_cache(bl)
stopifnot(env$.pb_variant_cache_valid(cache, "hvg2000_bl"))
legacy <- bl; legacy$blacklist_policy <- NULL
write_cache(legacy)
stopifnot(!env$.pb_variant_cache_valid(cache, "hvg2000_bl"))
stale <- bl; stale$blacklist_policy$source_md5 <- paste(rep("0", 32), collapse = "")
write_cache(stale)
stopifnot(!env$.pb_variant_cache_valid(cache, "hvg2000_bl"))
write_cache(ordinary$hvg2000)
stopifnot(env$.pb_variant_cache_valid(cache, "hvg2000"))

# BL result reuse also binds the policy and preparation identity; a valid
# checksum on a historical no-op result must not bypass the new preparation.
results_dir <- tempfile("pseudobulk-blacklist-results-")
dir.create(results_dir)
all_variants <- prepare(env$PB_VARIANT_NAMES)
score_calls <- 0L
env$process_pseudobulk_fig <- function(pb, labels, ...) {
  score_calls <<- score_calls + 1L
  list(feat_mat = pb, dist_mat = dist(pb), labels = labels)
}
labels <- setNames(factor(c("A", "B", "A")), metadata$Sample)
run <- function() suppressWarnings(env$run_pseudobulk_hpc(
  labels, all_variants, results_dir = results_dir, ds = "Toy"
))
shard_dir <- tempfile("pseudobulk-blacklist-shard-")
dir.create(shard_dir)
shard <- env$run_pseudobulk_hpc(
  labels, blacklist, results_dir = shard_dir, ds = "Toy",
  variants = "hvg2000_bl", ct_col_low_res = "must_not_read_LR",
  ct_col_high_res = "must_not_read_HR"
)
stopifnot(identical(names(shard), "Pseudobulk_hvg2000_bl"),
          identical(sort(list.files(shard_dir)), sort(c(
            "Toy_Pseudobulk_hvg2000_bl.rds", "Toy_Pseudobulk_hvg2000_bl.rds.md5"
          ))))
assert_error(env$run_pseudobulk_hpc(
  labels, blacklist, results_dir = shard_dir, ds = "Toy",
  variants = "hvg2000_bl", batch_mode = TRUE), "ordinary")
unlink(shard_dir, recursive = TRUE)
first <- run()
stopifnot(identical(first$Pseudobulk_hvg2000_bl$blacklist_policy, resolved$policy))
score_calls <- 0L
second <- run()
stopifnot(score_calls == 0L,
          identical(first$Pseudobulk_hvg2000_bl, second$Pseudobulk_hvg2000_bl))
bl_result_path <- file.path(results_dir, "Toy_Pseudobulk_hvg2000_bl.rds")
historical <- first$Pseudobulk_hvg2000_bl
historical$blacklist_policy <- NULL
save_rds_atomic(historical, bl_result_path)
third <- run()
stopifnot(score_calls == 1L,
          identical(third$Pseudobulk_hvg2000_bl$blacklist_policy, resolved$policy))
unlink(results_dir, recursive = TRUE)
unlink(c(input, cache, paste0(cache, ".md5")))
cat("PASS: direct/legacy masks, ordinary equivalence, independent BL cost, and policy-aware reuse\n")
