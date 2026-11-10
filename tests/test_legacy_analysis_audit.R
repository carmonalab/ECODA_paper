# Small synthetic H5AD obs and sample-level bundles only. No production reads.
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/seurat_utils.R")
source("src/notebooks/audit_legacy_analysis_inputs.R")
suppressPackageStartupMessages(library(dplyr))
reject <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
root <- tempfile("legacy-science-audit-"); dir.create(root); root <- normalizePath(root)
md5 <- function(path) unname(tools::md5sum(path))
Sys.setenv(RETICULATE_PYTHON = normalizePath(".pixi/envs/default/bin/python"))
module <- normalizePath("src/utils/py")
ids <- paste0("s", 1:6)
obs <- data.frame(donor = rep(ids, each = 4), condition = rep(rep(c("a", "b"), each = 3), each = 4),
  annotation = rep(c("T", "T", "B", "B"), 6), stringsAsFactors = FALSE)
obs$annotation[1] <- "Unknown"
labels <- setNames(rep(c("a", "b"), each = 3), ids)
entry <- list(sample_col = "donor", label_col = "condition", cell_type_high_res = "annotation", use_for_benchmark = TRUE)
reticulate::py_run_string("import h5py, numpy as np\ndef make_audit_obs(path, donor, condition, annotation):\n    with h5py.File(path, 'w') as h:\n        o=h.create_group('obs')\n        for name, values in [('donor',donor),('condition',condition),('annotation',annotation)]:\n            d=o.create_dataset(name,data=np.array(values,dtype=h5py.string_dtype()))\n            d.attrs['encoding-type']='string-array'\n        h.create_dataset('X',data=np.ones((len(donor),3)))\n        h.create_group('layers').create_dataset('counts',data=np.ones((len(donor),3)))\n        h.create_group('obsm').create_dataset('RNA_embedding',data=np.ones((len(donor),3)))\n_old_audit_getitem=h5py.Dataset.__getitem__\ndef obs_only_audit_getitem(self,key):\n    assert self.name.startswith('/obs/'), self.name\n    return _old_audit_getitem(self,key)\nh5py.Dataset.__getitem__=obs_only_audit_getitem\n")
inputs <- setNames(file.path(root, paste0(c("D", "B"), ".h5ad")), c("D", "B"))
for (p in inputs) {
  reticulate::py$make_audit_obs(p, obs$donor, obs$condition, obs$annotation)
  writeLines(c(paste0("MD5=", md5(p)), paste0("SIZE=", file.info(p)$size), paste0("PATH=", p)), paste0(p, ".md5"))
}
loaded_obs <- .legacy_read_obs(inputs[[1]], names(obs), module)
stopifnot(identical(loaded_obs, obs))
calc_sep_score <- function(...) stop("AUDIT MUST NOT SCORE")
counts <- as.matrix(get_ct_comp_df(obs, "donor", "annotation"))
imputed <- counts + .5; clr <- log(imputed) - rowMeans(log(imputed))
scores <- setNames(as.list(rep(.5, 5)), evaluation_metrics())
author <- list(feat_mat = as.data.frame(clr), counts = as.data.frame(imputed), labels = labels,
  dist_mat = dist(clr), scores = scores, exec_time = 3, timing_id = "original:fit", shared_time_secs = 4)
null <- author; null$labels <- setNames(labels[c(1,4,2,5,3,6)], ids)
rna <- matrix(1:48, 6, dimnames = list(ids, paste0("RNA_GENE_", 1:8)))
pb <- list(feat_mat = rna, counts = rna, RNA_hidden = rna, labels = labels, dist_mat = dist(rna), scores = scores)
scitd <- pb; scitd$feat_mat <- rna[-6, ]; scitd$dist_mat <- dist(scitd$feat_mat); scitd$labels <- labels[-6]
meta <- list(labels = labels, n_samples = 6L, n_cells = 24L,
  cells_per_sample = setNames(rep(4L, 6), ids), n_cell_types_high_res = 3L)
composition <- file.path(root, "D_composition.rds"); saveRDS(list(ECODA_authors_HR = author,
  ECODA_authors_HR_NULL = null, Avg_PCA_embedding = pb, UNSELECTED_HIDDEN_RNA = rna), composition)
pb_path <- file.path(root, "D_pseudobulk.rds"); saveRDS(list(Pseudobulk_hvg2000 = pb, scITD_hvg2000_factors5 = scitd), pb_path)
meta_path <- file.path(root, "D_metadata.rds"); saveRDS(meta, meta_path)
other_path <- file.path(root, "B_composition.rds"); saveRDS(list(Avg_PCA_embedding = pb), other_path)
emb <- file.path(root, "D_emb.feather"); arrow::write_feather(data.frame(`__index_level_0__` = ids, Dim_1 = 1:6, Dim_2 = 3:8, check.names = FALSE), emb)
dist_path <- file.path(root, "D_dist.feather"); arrow::write_feather(cbind(as.data.frame(as.matrix(pb$dist_mat)), data.frame(Sample = ids)), dist_path)
make_row <- function(dataset, method, path, kind = "rds_bundle") data.frame(dataset = dataset, view = "benchmark_analysis", method_key = method,
  kind = kind, artifact_path = path, bundle_key = if (kind == "rds_bundle") method else "",
  sample_id_col = if (kind == "embedding_feather") "__index_level_0__" else if (kind == "distance_feather") "Sample" else "",
  md5 = md5(path), size = as.character(file.info(path)$size),
  current_input_path = inputs[[dataset]], current_input_md5 = md5(inputs[[dataset]]), current_input_size = as.character(file.info(inputs[[dataset]])$size),
  current_input_proof = paste0(inputs[[dataset]], ".md5"), origin = "legacy_reuse", status = "BLOCKED_UNTIL_SCIENTIFIC_REVIEW", stringsAsFactors = FALSE)
rows <- do.call(rbind, list(make_row("D", "ECODA_authors_HR", composition), make_row("D", "ECODA_authors_HR_NULL", composition),
  make_row("D", "Avg_PCA_embedding", composition), make_row("D", "Pseudobulk_hvg2000", pb_path),
  make_row("D", "scITD_hvg2000_factors5", pb_path), make_row("D", "__dataset_metadata__", meta_path, "standalone_result"),
  make_row("D", "scPoli_hvg2000_dims15_highres", emb, "embedding_feather"), make_row("D", "PILOT_hvg2000", dist_path, "distance_feather"),
  make_row("B", "Avg_PCA_embedding", other_path)))
files <- rows[!duplicated(rows$artifact_path), c("dataset", "artifact_path", "md5", "size")]
files$needed_bundle_keys <- vapply(files$artifact_path, function(p) paste(rows$bundle_key[rows$artifact_path == p & rows$kind == "rds_bundle"], collapse = ";"), character(1))
timing <- file.path(root, "execution_times.feather")
timing_frame <- data.frame(dataset = rows$dataset[-6], method = rows$method_key[-6], time_secs = 1, mem_GB = NA_real_)
timing_frame <- rbind(timing_frame, data.frame(dataset = "UNRELATED", method = "irrelevant", time_secs = NA_real_, mem_GB = -1))
arrow::write_feather(timing_frame, timing)
files <- rbind(files, data.frame(dataset = "__all__", artifact_path = timing, md5 = md5(timing), size = as.character(file.info(timing)$size), needed_bundle_keys = ""))
pairs <- data.frame(dataset = "D", author_path = composition, author_md5 = md5(composition), author_bundle_key = "ECODA_authors_HR",
  metadata_path = meta_path, metadata_md5 = md5(meta_path), current_input_path = inputs[[1]], current_input_md5 = md5(inputs[[1]]))
datasets <- list(D = entry, B = entry)
reads <- character(); base_reader <- base::readRDS
readRDS <- function(file, ...) { reads <<- c(reads, file); base_reader(file, ...) }
obs_calls <- list()
obs_reader <- function(path, columns) { obs_calls[[length(obs_calls) + 1L]] <<- columns; .legacy_read_obs(path, columns, module) }
protected <- tools::md5sum(c(files$artifact_path, inputs, paste0(inputs, ".md5")))
set.seed(114); rng_before <- .Random.seed
report <- audit_legacy_analysis_inputs(rows, files, pairs, datasets, obs_reader)
rm(readRDS)
stopifnot(report$compatible, identical(rng_before, .Random.seed), report$scientifically_approved_rows == 0L, report$score_calculations == 0L,
  report$model_fits == 0L, !report$h5ad_expression_values_read, identical(protected, tools::md5sum(names(protected))),
  length(reads) == 4L, !anyDuplicated(reads), identical(obs_calls[[1]], names(obs)), identical(obs_calls[[2]], c("donor", "condition")),
  report$rows$`D/ECODA_authors_HR_NULL`$null_features_and_distances_match_author,
  identical(unlist(report$rows$`D/scITD_hvg2000_factors5`$dropped_sample_ids), "s6"),
  report$rows$`D/scITD_hvg2000_factors5`$dropped_ids_require_explicit_review,
  report$author_hr_pairs$D$valid_annotated_cells == 23L, report$author_hr_pairs$D$n_recorded_categories == 3L,
  identical(report$author_hr_pairs$D$sample_ids, unname(names(meta$cells_per_sample))),
  identical(report$rows$`D/scPoli_hvg2000_dims15_highres`$sample_id_col, "__index_level_0__"),
  report$timing$selected_costs_valid, length(report$timing$invalid_unrelated_rows) == 1L,
  report$rows$`D/ECODA_authors_HR`$original_timing$timing_id == "original:fit")
encoded <- jsonlite::toJSON(report, auto_unbox = TRUE, null = "null")
stopifnot(!grepl("RNA_GENE_|RNA_hidden", encoded), !any(vapply(report$rows, function(x) !is.null(x$feat_mat), logical(1))))
bad <- files; bad$needed_bundle_keys[1] <- "ECODA_authors_HR"; reject(audit_legacy_analysis_inputs(rows, bad, pairs, datasets, obs_reader))
bad <- rows; bad$origin[1] <- "new"; reject(audit_legacy_analysis_inputs(bad, files, pairs, datasets, obs_reader))
bad <- pairs; bad$author_md5 <- strrep("0", 32); reject(audit_legacy_analysis_inputs(rows, files, bad, datasets, obs_reader))
bad <- rows; bad$sample_id_col[7] <- ""; reject(audit_legacy_analysis_inputs(bad, files, pairs, datasets, obs_reader))
bad <- rows; bad$sample_id_col[7] <- "missing_column"
bad_report <- audit_legacy_analysis_inputs(bad, files, pairs, datasets, obs_reader)
stopifnot(!bad_report$compatible, bad_report$rows$`D/scPoli_hvg2000_dims15_highres`$status == "INCOMPATIBLE")
headers <- list(schema = "ecoda_p6_feather_header_observation_v1", rows = lapply(which(rows$kind %in% c("embedding_feather", "distance_feather")), function(i) {
  row <- as.list(rows[i, ]); row$fields <- list(list(name = row$sample_id_col), list(name = "other"))
  row$pandas <- list(index_columns = row$sample_id_col); row$record_batches_read <- 0L; row
}))
.legacy_header_scope(rows, headers)
bad_headers <- headers; bad_headers$rows[[1]]$pandas$index_columns <- "wrong"
reject(.legacy_header_scope(rows, bad_headers))
reject(.legacy_labels(transform(obs, condition = replace(condition, 2, "b")), entry))
bad_null <- null; bad_null$labels[] <- "a"; reject(.legacy_bundle(bad_null, rows[2, ], labels))
bad_pb <- pb; bad_pb$labels[1] <- "b"; reject(.legacy_bundle(bad_pb, rows[4, ], labels))
reject(.legacy_bundle(scitd, rows[4, ], labels))
bad_pb <- pb; bad_pb$dist_mat <- as.matrix(pb$dist_mat); bad_pb$dist_mat[1,2] <- -1; reject(.legacy_bundle(bad_pb, rows[4, ], labels))
unknown <- pb; unknown$scores$lisi_score <- 4
stopifnot(.legacy_bundle(unknown, rows[4, ], labels)$score_disposition == "needs_validated_distance_rescore")
unknown$scoring_policy <- evaluation_scoring_policy(); stopifnot("lisi_score" %in% .legacy_bundle(unknown, rows[4, ], labels)$invalid_core_scores)
unknown$scores$lisi_score <- NaN; stopifnot(.legacy_bundle(unknown, rows[4, ], labels)$score_disposition == "needs_validated_distance_rescore")
bad_author <- author; bad_author$counts[1,1] <- bad_author$counts[1,1] + 1; reject(.legacy_pair(bad_author, meta, obs, entry, labels))
bad_author <- author; bad_author$feat_mat[1,1] <- bad_author$feat_mat[1,1] + 1; reject(.legacy_pair(bad_author, meta, obs, entry, labels))
bad_meta <- meta; bad_meta$n_cells <- 25L; reject(.legacy_pair(author, bad_meta, obs, entry, labels))
# Actual Adams metadata IDs are already the configured Sample authority.
# Reapplying the raw-alias standardizer would corrupt58 validated hyphen IDs.
adams_ids <- strsplit("g001C g002C g003C g010I g021I g022C-a g022C-b g222C g022I g025I g034C g034I-a g040I-b g040I g041I-b g051I-a g051I g052CO-a g053I-d g053I-n g053I g056CO g063I-b g065C g081C g092C g098C-a g098C-b g123I g133C-a g135I-a g135I-b g137C-a g137C-b g137C g138I-a g145I-a g152CO-a g153CO-a g153CO-b g157I-a g157I-b g157I g158I-b g160C g166I-a g174I-a g177I g178CO g179I g184CO-a g184CO-b g186CO-b g192C-a g192CO-a g192CO-b g192CO g193CO-a g193CO-b g194CO-a g207CO g208C g209I-a g210CO g212I-a g214I-a g217CO-a g218C-a g218C-b g219C-a g219C-b g221I-a g222I-a g222I-b g225I-a g225I-b g226C-a g226C-b g228I-a g228I-b g235CO g237CO g238CO g23CO g253C g296C g29I g388C g396C g439C-b g439C g454C g465C g47I-a g47I-b g483C g484C g49I-a g49I-b g59I", " ", fixed = TRUE)[[1]]
adams_obs <- data.frame(Sample = adams_ids, condition = rep(c("a", "b"), 50), stringsAsFactors = FALSE)
adams_entry <- list(sample_col = "Sample", label_col = "condition")
adams_view <- .legacy_canonical_obs(adams_obs, adams_entry)
stopifnot(length(adams_ids) == 100L, sum(grepl("-", adams_ids, fixed = TRUE)) == 58L,
  identical(adams_view$obs, adams_obs), identical(adams_view$aliases$raw, adams_ids),
  identical(adams_view$aliases$canonical, adams_ids),
  identical(names(.legacy_labels(adams_view$obs, adams_entry)), adams_ids))
# Gong-style raw IDs differ from the producer's existing canonical Sample.
# Read both obs columns, prove the exact one-to-one production mapping, and
# change only this in-memory validation view. No result/source IDs are renamed.
canonical_path <- file.path(root, "canonical.h5ad")
raw_samples <- paste0("p-", rep(1:6, each = 4))
canonical_samples <- gsub("-", "_", raw_samples, fixed = TRUE)
reticulate::py$make_audit_obs(canonical_path, raw_samples, obs$condition, obs$annotation)
reticulate::py_run_string("def add_audit_canonical(path, values):\n    with h5py.File(path, 'a') as h:\n        h['obs'].create_dataset('Sample',data=np.array(values,dtype=h5py.string_dtype()))\n")
reticulate::py$add_audit_canonical(canonical_path, canonical_samples)
canonical_md5 <- md5(canonical_path)
raw_obs <- .legacy_read_obs(canonical_path, names(obs), module)
stopifnot(identical(raw_obs$donor, raw_samples), identical(raw_obs$Sample, canonical_samples))
mapped <- .legacy_canonical_obs(raw_obs, entry)
stopifnot(identical(mapped$obs$donor, canonical_samples), identical(raw_obs$donor, raw_samples),
  identical(mapped$aliases$raw, unique(raw_samples)), identical(mapped$aliases$canonical, unique(canonical_samples)))
canonical_ids <- unique(canonical_samples); canonical_labels <- .legacy_labels(mapped$obs, entry)
canonical_author <- author
rownames(canonical_author$counts) <- rownames(canonical_author$feat_mat) <- canonical_ids
attr(canonical_author$dist_mat, "Labels") <- canonical_ids; names(canonical_author$labels) <- canonical_ids
canonical_meta <- meta; names(canonical_meta$labels) <- names(canonical_meta$cells_per_sample) <- canonical_ids
stopifnot(.legacy_pair(canonical_author, canonical_meta, mapped$obs, entry, canonical_labels)$status == "COMPATIBLE_PENDING_REVIEW",
  identical(.legacy_bundle(canonical_author, rows[1, ], canonical_labels)$sample_ids, canonical_ids),
  identical(md5(canonical_path), canonical_md5))
bad_mapping <- raw_obs; bad_mapping$Sample[1] <- "foreign"; reject(.legacy_canonical_obs(bad_mapping, entry))
bad_mapping <- raw_obs; bad_mapping$Sample <- NULL; reject(.legacy_canonical_obs(bad_mapping, entry))
collision <- data.frame(donor = c("p-1", "p_1"), Sample = c("p_1", "p_1"))
reject(.legacy_canonical_obs(collision, entry))
# GloScope carries incidental names on its distance-ID vector. Value/order
# equality remains strict, including actual label and permutation negatives.
glo <- pb; glo$dist_mat <- as.matrix(glo$dist_mat)
names(dimnames(glo$dist_mat)[[1]]) <- ids; names(dimnames(glo$dist_mat)[[2]]) <- ids
glo_before <- serialize(glo, NULL)
stopifnot(identical(.legacy_bundle(glo, rows[4, ], labels)$sample_ids, ids), identical(serialize(glo, NULL), glo_before))
bad_glo <- glo; bad_glo$labels <- bad_glo$labels[6:1]; reject(.legacy_bundle(bad_glo, rows[4, ], labels))
# Kim-style stored count axes may differ from the features' valid axes.
aligned_author <- author
aligned_author$counts <- aligned_author$counts[6:1, ncol(imputed):1, drop = FALSE]
aligned_before <- serialize(aligned_author, NULL)
stopifnot(.legacy_pair(aligned_author, meta, obs, entry, labels)$status == "COMPATIBLE_PENDING_REVIEW",
  identical(serialize(aligned_author, NULL), aligned_before))
bad_axes <- aligned_author; colnames(bad_axes$feat_mat)[1] <- "foreign"
reject(.legacy_pair(bad_axes, meta, obs, entry, labels))
bad_axes <- aligned_author; bad_axes$feat_mat <- as.matrix(bad_axes$feat_mat)
rownames(bad_axes$feat_mat)[1] <- rownames(bad_axes$feat_mat)[2]
reject(.legacy_pair(bad_axes, meta, obs, entry, labels))
bad_axes <- aligned_author; bad_axes$dist_mat <- as.matrix(author$dist_mat)[6:1, 6:1]
reject(.legacy_pair(bad_axes, meta, obs, entry, labels))
alias <- file.path(root, "alias.rds"); stopifnot(file.symlink(composition, alias)); reject(.legacy_checked_read(alias, md5(composition), file.info(composition)$size, readRDS))
backup <- readBin(composition, "raw", n = file.info(composition)$size)
reject(.legacy_checked_read(composition, md5(composition), file.info(composition)$size, function(p) { x <- readRDS(p); saveRDS(list(changed = TRUE), p); x }))
writeBin(backup, composition)
stopifnot(identical(protected, tools::md5sum(names(protected))))
reticulate::py_run_string("h5py.Dataset.__getitem__=_old_audit_getitem")
unlink(root, recursive = TRUE)
cat("legacy analysis science audit: PASS; exact scope, one-read streaming, obs-only configured columns, NULL/scITD, counts/CLR/distances, no scoring/RNA export, immutable proofs\n")
