source("src/utils/analysis_release.R")
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
tmp <- tempfile("release-"); dir.create(tmp); tmp <- normalizePath(tmp)
md5 <- function(path) unname(tools::md5sum(path))
ids <- paste0("s", 1:6)
features <- matrix(seq_len(12), 6, dimnames = list(ids, c("x", "y")))
labels <- setNames(rep(c("a", "b"), each = 3), ids)
bundle <- list(feat_mat = features, dist_mat = dist(features), labels = labels,
               scores = list(anosim_score = 0.7))
aggregate <- file.path(tmp, "old.rds"); saveRDS(list(ordinary = bundle), aggregate)
standalone <- file.path(tmp, "new_bl.rds"); saveRDS(bundle, standalone)
evidence <- file.path(tmp, "review.txt"); writeLines("Bounded row review; input and artifact verified; no invented producer", evidence)
ledger <- data.frame(dataset = "D", view = "benchmark_analysis", method_key = c("ordinary", "Pseudobulk_hvg2000_bl"),
  kind = c("rds_bundle", "standalone_result"), artifact_path = c(aggregate, standalone),
  bundle_key = c("ordinary", ""), run_root = tmp, md5 = c(md5(aggregate), md5(standalone)),
  receipt_path = file.path(tmp, c("ordinary.receipt.rds", "bl.receipt.rds")), receipt_md5 = "")
scope <- ledger[c("dataset", "view", "method_key")]
scope$require_new <- c(FALSE, TRUE); scope$current_input_md5 <- strrep("a", 32)
scope$source_policy <- "input-v1"; scope$status <- "approved"; scope$blocked_reason <- ""
for (i in 1:2) {
  receipt <- c(as.list(ledger[i, c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")]),
    list(schema = "ecoda_analysis_release_receipt_v1", status = "reviewed", origin = c("legacy_reuse", "new")[i],
         input_md5 = scope$current_input_md5[i], source_policy = "input-v1", evidence = "bounded review",
         evidence_path = evidence, evidence_md5 = md5(evidence)))
  saveRDS(receipt, ledger$receipt_path[i]); ledger$receipt_md5[i] <- md5(ledger$receipt_path[i])
}
loaded <- load_analysis_release(ledger, scope)
stopifnot(length(loaded) == 2, identical(loaded[[1]]$result, bundle), identical(loaded[[2]]$result$labels, labels))
# Declared aliases cannot escape the reviewed root or mask distinct identities.
symlink <- file.path(tmp, "linked.rds")
stopifnot(file.symlink(standalone, symlink))
bad <- ledger; bad$artifact_path[2] <- symlink
fail(load_analysis_release(bad, scope))
bad <- ledger; bad$artifact_path[2] <- file.path(tmp, ".", basename(standalone))
fail(load_analysis_release(bad, scope))
root_alias <- file.path(tmp, "linked-root"); stopifnot(file.symlink(tmp, root_alias))
bad <- ledger; bad$run_root <- root_alias
fail(load_analysis_release(bad, scope))
receipt_alias <- file.path(tmp, "receipt-alias.rds"); stopifnot(file.symlink(ledger$receipt_path[2], receipt_alias))
bad <- ledger; bad$receipt_path[2] <- receipt_alias
fail(load_analysis_release(bad, scope))
evidence_alias <- file.path(tmp, "evidence-alias.txt"); stopifnot(file.symlink(evidence, evidence_alias))
receipt <- readRDS(ledger$receipt_path[2]); receipt$evidence_path <- evidence_alias
alias_receipt <- file.path(tmp, "aliased-evidence.receipt.rds"); saveRDS(receipt, alias_receipt)
bad <- ledger; bad$receipt_path[2] <- alias_receipt; bad$receipt_md5[2] <- md5(alias_receipt)
fail(load_analysis_release(bad, scope))
# A symlink to an existing external artifact fails even with a valid checksum.
outside <- tempfile("outside-release-"); saveRDS(bundle, outside)
escape <- file.path(tmp, "escape.rds"); stopifnot(file.symlink(outside, escape))
bad <- ledger; bad$artifact_path[2] <- escape; bad$md5[2] <- md5(outside)
fail(load_analysis_release(bad, scope)); unlink(outside)
fail(load_analysis_release(ledger[1, ], scope))
fail(load_analysis_release(rbind(ledger, ledger[1, ]), scope))
bad <- scope; bad$status[1] <- "blocked"; bad$blocked_reason[1] <- "HiTME count mismatch pending review"
fail(load_analysis_release(ledger, bad))
bad <- scope; bad$require_new[1] <- TRUE
fail(load_analysis_release(ledger, bad))
bad <- scope; bad$current_input_md5[1] <- strrep("b", 32)
fail(load_analysis_release(ledger, bad))
bad <- scope; bad$source_policy[1] <- "changed-input-policy"
fail(load_analysis_release(ledger, bad))
bad <- ledger; bad$bundle_key[1] <- "missing"
fail(load_analysis_release(bad, scope))
# A missing required standalone never falls back to an old aggregate key.
stopifnot(file.rename(standalone, paste0(standalone, ".saved")))
fail(load_analysis_release(ledger, scope))
stopifnot(file.rename(paste0(standalone, ".saved"), standalone))
writeLines("tampered review", evidence)
fail(load_analysis_release(ledger, scope))
# The same finite ledger reads exact distance/embedding Feather rows, without
# conflating sample identity with a numeric feature column.
writeLines("Bounded row review; input and artifact verified; no invented producer", evidence)
for (kind in c("distance_feather", "embedding_feather")) {
  row <- ledger[2, ]; row$kind <- kind; row$method_key <- "feather_method"
  row$artifact_path <- file.path(tmp, paste0(kind, ".feather")); row$sample_id_col <- "sample_id"
  matrix <- if (kind == "distance_feather") as.matrix(bundle$dist_mat) else features
  arrow::write_feather(cbind(data.frame(sample_id = ids), as.data.frame(matrix)), row$artifact_path)
  row$md5 <- md5(row$artifact_path)
  receipt <- readRDS(ledger$receipt_path[2])
  for (field in c("method_key", "kind", "artifact_path", "md5")) receipt[[field]] <- row[[field]]
  saveRDS(receipt, row$receipt_path); row$receipt_md5 <- md5(row$receipt_path)
  expected <- scope[2, ]; expected$method_key <- row$method_key
  read <- load_analysis_release(row, expected)[[1]]
  stopifnot(isTRUE(all.equal(as.matrix(read$result$dist_mat), as.matrix(bundle$dist_mat))))
  if (kind == "embedding_feather") stopifnot(identical(read$result$feat_mat, features))
}
unlink(tmp, recursive = TRUE)
cat("analysis release ledger: OK\n")
