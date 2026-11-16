#!/usr/bin/env Rscript
# Focused scientific contracts; execute only the notebook's pure core chunk.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg)) setwd(normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..")))

core <- new.env(parent = baseenv())
draft <- Sys.getenv("BECM_CORE_FILE", "")
if (nzchar(draft)) {
  sys.source(draft, envir = core)
} else {
  notebook <- readLines("notebooks/batch_effect_ECODA_becm_comparison.rmd", warn = FALSE)
  start <- grep("^```\\{r becm-core[,}]", notebook)
  stopifnot(length(start) == 1L)
  end <- start + which(grepl("^```[[:space:]]*$", notebook[(start + 1L):length(notebook)]))[1L]
  stopifnot(is.finite(end), end > start + 1L)
  eval(parse(text = notebook[(start + 1L):(end - 1L)]), envir = core)
}
required <- c("becm_screen", "becm_drop", "becm_variant", "becm_build_variants",
  "becm_anosim", "becm_score_variants", "becm_native_input")
stopifnot(all(vapply(required, function(name) is.function(core[[name]]), logical(1))))
assert_error <- function(expr) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
}
ids <- paste0("sample_", c(12, 2, 10, 4, 8, 6, 7, 5, 9, 3, 11, 1))
counts <- cbind(mono = c(1:6, 70:75), T = c(100:105, 10:15),
  B = 40:51, sparse50 = c(rep(0, 6), 1:6),
  sparse42 = c(rep(0, 5), 1:7), stable = rep(40, 12))
rownames(counts) <- ids
clr <- log(counts + 0.5) - rowMeans(log(counts + 0.5))
batch2 <- rep(c("A", "B"), each = 6)
batch3 <- rep(c("A", "B", "C"), each = 4)

# Independent checks of the two tests, full initial BH family and effect size.
screen2 <- core$becm_screen(clr, batch2, counts)
screen3 <- core$becm_screen(clr, batch3, counts)
expected_wilcox <- apply(clr, 2, function(y) stats::wilcox.test(
  y ~ factor(batch2), exact = FALSE, correct = TRUE)$p.value)
expected_kw <- apply(clr, 2, function(y) stats::kruskal.test(y ~ factor(batch3))$p.value)
expected_span <- apply(clr, 2, function(y) diff(range(tapply(y, batch3, stats::median))))
stopifnot(identical(screen2$cell_type, colnames(counts)),
  all(screen2$test == "Wilcoxon rank-sum"), all(screen3$test == "Kruskal-Wallis"),
  isTRUE(all.equal(screen2$p, unname(expected_wilcox))),
  isTRUE(all.equal(screen3$p, unname(expected_kw))),
  isTRUE(all.equal(screen2$p_bh, unname(stats::p.adjust(expected_wilcox, "BH")))),
  isTRUE(all.equal(screen3$p_bh, unname(stats::p.adjust(expected_kw, "BH")))),
  isTRUE(all.equal(screen3$median_clr_range, unname(expected_span))),
  screen2$zero_fraction[screen2$cell_type == "sparse50"] == 0.5,
  screen2$zero_fraction[screen2$cell_type == "sparse42"] == 5 / 12)

# Deterministic selection boundary: significance AND inclusive effect cutoff.
mock <- screen2
mock$p_bh <- mock$q <- c(0.01, 0.049, 0.1, 0.5, 0.03, 1)
mock$median_clr_range <- c(2, 1, 10, 4, 0.5, 0)
dropped <- core$becm_drop(clr, mock, cutoff = 2)
stopifnot(identical(dropped$removed_types, "mono"),
  identical(dropped$n_initial, ncol(clr)), dropped$n_remaining == 5L,
  !"B" %in% dropped$removed_types, # large effect alone is insufficient
  !"T" %in% dropped$removed_types, # significance alone is insufficient
  identical(rownames(dropped$feat_mat), ids))
boundary <- core$becm_drop(clr, mock, cutoff = 1)
stopifnot(setequal(boundary$removed_types, c("mono", "T")))
mock_at_alpha <- mock
mock_at_alpha$p_bh[mock_at_alpha$cell_type == "T"] <- 0.05
stopifnot(!"T" %in% core$becm_drop(clr, mock_at_alpha, cutoff = 1)$removed_types)

# >=50% zeros is an independent UNION, not conditional on significance.
sparse <- core$becm_drop(clr, mock, cutoff = 2, drop_sparse = TRUE)
stopifnot(setequal(sparse$removed_types, c("mono", "sparse50")),
  identical(sparse$sparse_removed_types, "sparse50"),
  !"sparse42" %in% sparse$removed_types)
retained <- counts[, !colnames(counts) %in% sparse$removed_types, drop = FALSE]
expected_clr <- log(retained + 0.5) - rowMeans(log(retained + 0.5))
stopifnot(isTRUE(all.equal(sparse$feat_mat, expected_clr)),
  max(abs(rowSums(sparse$feat_mat))) < 1e-12)

# Invalid low-dimensional outputs carry an explicit status and no fake scores.
too_few <- mock
too_few$p_bh <- too_few$q <- c(rep(0.01, 5), 1)
too_few$median_clr_range <- rep(3, 6)
one <- core$becm_drop(clr, too_few, 2)
v_one <- core$becm_variant(one$feat_mat, "authors_drop", "corrected", 2, one)
stopifnot(one$n_remaining == 1L, grepl("^UNSCORABLE", one$status), is.null(v_one$dist_mat))
zero <- core$becm_drop(clr, transform(too_few, p_bh = rep(0.01, 6)), 2)
stopifnot(zero$n_remaining == 0L, grepl("^UNSCORABLE", zero$status))
identical_rows <- matrix(rep(clr[1, ], each = nrow(clr)), nrow = nrow(clr),
  dimnames = dimnames(clr))
same <- core$becm_variant(identical_rows, "fixture", "corrected")
stopifnot(grepl("^UNSCORABLE", same$status), is.null(same$dist_mat))
metadata <- data.frame(Sample = ids, bio = rep(c("case", "control"), 6), batch = batch2)
unscorable <- core$becm_score_variants(list(one = v_one), metadata, "Fixture")
stopifnot(is.na(unscorable$biological_anosim), is.na(unscorable$batch_anosim),
  grepl("^UNSCORABLE", unscorable$status))

# ANOSIM R is a rank contrast on the exact named sample distance matrix.
stopifnot(requireNamespace("vegan", quietly = TRUE))
distance <- stats::dist(clr)
rank_anosim <- function(distance, groups) {
  rank_matrix <- matrix(NA_real_, length(groups), length(groups))
  lower <- lower.tri(rank_matrix)
  rank_matrix[lower] <- rank(as.numeric(distance))
  same_group <- outer(groups, groups, `==`)[lower]
  ranked <- rank_matrix[lower]
  (mean(ranked[!same_group]) - mean(ranked[same_group])) /
    (length(groups) * (length(groups) - 1) / 4)
}
stopifnot(isTRUE(all.equal(core$becm_anosim(distance, batch2),
  rank_anosim(distance, batch2))))
normal <- core$becm_variant(clr, "fixture", "uncorrected")
scored <- core$becm_score_variants(list(original = normal), metadata, "Fixture")
stopifnot(scored$status == "PASS",
  isTRUE(all.equal(scored$biological_anosim, rank_anosim(distance, metadata$bio))),
  isTRUE(all.equal(scored$batch_anosim, rank_anosim(distance, metadata$batch))))
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
expect_error(core$becm_score_variants(list(original = normal), metadata[rev(seq_len(nrow(metadata))), ], "Fixture"))

# Screen once per initial annotation family, never after a threshold removal.
input <- list(clr = clr, counts = counts)
inputs <- setNames(rep(list(input), 5),
  c("authors", "leiden", "authors_limma", "harmony", "harmony_limma"))
screen_calls <- list()
original_screen <- core$becm_screen
core$becm_screen <- function(clr, batch, counts) {
  screen_calls[[length(screen_calls) + 1L]] <<- colnames(clr)
  original_screen(clr, batch, counts)
}
built <- core$becm_build_variants(inputs, batch2)
stopifnot(length(screen_calls) == 2L,
  all(vapply(screen_calls, identical, logical(1), colnames(clr))),
  length(built$variants) == 24L, !anyDuplicated(names(built$variants)),
  sum(vapply(built$variants, function(v) v$state == "uncorrected", logical(1))) == 2L)
for (family in c("ECODA_authors_sigdiff_ctdrop_", "ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_")) {
  keys <- paste0(family, vapply(seq(0, 4, 0.5), core$becm_threshold_id, character(1)))
  stopifnot(all(keys %in% names(built$variants)))
  for (k in seq_len(length(keys) - 1L)) {
    stopifnot(all(built$variants[[keys[k + 1L]]]$removed_types %in%
      built$variants[[keys[k]]]$removed_types))
  }
}
for (v in built$variants) if (v$status == "PASS") {
  stopifnot(identical(attr(v$dist_mat, "Labels"), ids))
}
core$becm_screen <- original_screen

# Biology is unavailable to the screen/drop/transform APIs.
for (fun in c("becm_screen", "becm_drop", "becm_build_variants")) {
  stopifnot(!any(c("bio", "biology", "labels", "Status") %in% names(formals(core[[fun]]))))
}
core$Status <- rev(metadata$bio)
stopifnot(identical(core$becm_build_variants(inputs, batch2), built))
assert_error(core$becm_build_variants(inputs, batch2, bio = metadata$bio))
assert_error(core$becm_screen(clr, rep("A", 12), counts))
assert_error(core$becm_screen(clr, replace(batch2, 1, NA), counts))
missing_ids <- clr
rownames(missing_ids) <- NULL
missing_count_ids <- counts
rownames(missing_count_ids) <- NULL
assert_error(core$becm_screen(missing_ids, batch2, missing_count_ids))
bad_order <- inputs
bad_order$harmony$clr <- bad_order$harmony$clr[rev(ids), , drop = FALSE]
assert_error(core$becm_build_variants(bad_order, batch2))

# A native bundle's distance axis names are part of its sample identity.
local({
  directory <- tempfile("becm-native-fixture-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  path <- file.path(directory, "composition.rds")
  receipt_path <- file.path(directory, "receipt.rds")
  install <- function(distance) {
    saveRDS(list(ECODA_authors_HR = list(feat_mat = clr, counts = counts + 0.5,
      dist_mat = distance)), path)
    md5 <- unname(tools::md5sum(path))
    saveRDS(list(status = "reviewed", md5 = md5, artifact_path = path,
      bundle_key = "ECODA_authors_HR"), receipt_path)
    list(artifact_path = path, md5 = md5, receipt_path = receipt_path,
      receipt_md5 = unname(tools::md5sum(receipt_path)),
      bundle_key = "ECODA_authors_HR", role = "authors")
  }
  good <- install(stats::dist(clr))
  restored <- core$becm_native_input(good, ids, colnames(clr))
  stopifnot(identical(restored$counts, counts), isTRUE(all.equal(restored$clr, clr)))
  bad_distance <- stats::dist(clr)
  attr(bad_distance, "Labels") <- paste0("unrelated_", seq_along(ids))
  bad <- install(bad_distance)
  assert_error(core$becm_native_input(bad, ids, colnames(clr)))
})

cat("ECODA becm scientific contracts passed: tests/BH, inclusive median cutoff, sparse union, retained CLR, fixed screens, 24 variants, identity and unscorable states.\n")
