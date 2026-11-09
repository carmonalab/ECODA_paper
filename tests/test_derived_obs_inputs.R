# Real saved-label counts and five scores; synthetic H5AD obs only.
source("src/5_run_benchmark_methods/run_local_ecoda_derived.R")
reject <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
project <- normalizePath(".")
runner <- file.path(project, "src/5_run_benchmark_methods/run_local_ecoda_derived.R")
.ecoda_runner_path <- function() runner
root <- tempfile("derived-bound-obs-"); dir.create(root); root <- normalizePath(root)
Sys.setenv(PROJECT_ROOT = project, ECODA_DERIVED_ROOT = file.path(root, "derived"),
  RETICULATE_PYTHON = normalizePath(".pixi/envs/default/bin/python"))
dir.create(Sys.getenv("ECODA_DERIVED_ROOT"))
md5 <- function(p) unname(tools::md5sum(p))
sha <- function(p) digest::digest(file = p, algo = "sha256")
sample_ids <- paste0("id", c(10,2,1,12,4,3,11,8,5,9,7,6))
sample <- bio <- cluster50 <- cluster2 <- character()
for (i in seq_along(sample_ids)) {
  count <- c(i + 1L, 3L + i %% 4L, 28L - i - i %% 4L)
  values <- rep(c("0", "1", "2"), count)
  sample <- c(sample, rep(sample_ids[i], length(values)))
  bio <- c(bio, rep(if (i <= 6) "a" else "b", length(values)))
  cluster50 <- c(cluster50, values)
  cluster2 <- c(cluster2, c(rep("0", 3L + i), rep("1", 4L + i %% 3L), rep("2", length(values) - 7L - i - i %% 3L)))
}
cluster50[1] <- NA_character_
reticulate::py_run_string("import h5py, numpy as np\ndef make_derived_obs(path, samplecol, samples, labels, r50, h2):\n    with h5py.File(path,'w') as h:\n        o=h.create_group('obs')\n        for name, values in [(samplecol,samples),('biology',labels)]:\n            d=o.create_dataset(name,data=np.array(values,dtype=h5py.string_dtype()))\n            d.attrs['encoding-type']='string-array'\n        for name, values in [('leiden_res_50_benchmark_analysis_hvg2000',r50),('leiden_res_2_benchmark_analysis_hvg2000_harmony',h2)]:\n            g=o.create_group(name);g.attrs['encoding-type']='categorical'\n            g.create_dataset('categories',data=np.array(['0','1','2','UNUSED'],dtype=h5py.string_dtype()))\n            codes=[-1 if v is None or str(v)=='nan' or str(v)=='NA' else int(v) for v in values]\n            g.create_dataset('codes',data=np.array(codes,dtype='i4'))\n        h.create_dataset('X',data=np.ones((len(samples),3)))\n        h.create_group('layers').create_dataset('counts',data=np.ones((len(samples),3)))\n        h.create_group('obsm').create_dataset('X_pca',data=np.ones((len(samples),3)))\n_old_derived_getitem=h5py.Dataset.__getitem__\ndef derived_obs_only_getitem(self,key):\n    assert self.name.startswith('/obs/'), self.name\n    return _old_derived_getitem(self,key)\nh5py.Dataset.__getitem__=derived_obs_only_getitem\n")
paths <- setNames(file.path(root, c("D.h5ad", "E.h5ad")), c("D", "E"))
for (d in names(paths)) {
  reticulate::py$make_derived_obs(paths[[d]], if (d == "D") "donor" else "specimen.guid", sample, bio,
    as.list(cluster50), as.list(cluster2))
  writeLines(c(paste0("MD5=", md5(paths[[d]])), paste0("SIZE=", file.info(paths[[d]])$size), paste0("PATH=", paths[[d]])), paste0(paths[[d]], ".md5"))
}
config <- list(D = list(use_for_benchmark = TRUE, columns = list(sample = "donor", label = "biology", cell_type_high_res = "DO_NOT_READ"),
  views = list(benchmark_analysis = list(output_file_name = "D.h5ad"))),
  E = list(use_for_benchmark = TRUE, columns = list(sample = "specimen.guid", label = "biology"),
    views = list(benchmark_analysis = list(output_file_name = "E.h5ad"))),
  DISABLED = list(use_for_benchmark = FALSE, columns = list(sample = "Sample", label = "biology"),
    views = list(benchmark_analysis = list(output_file_name = "D.h5ad"))))
config_path <- file.path(root, "datasets.json"); jsonlite::write_json(config, config_path, auto_unbox = TRUE)
refs <- data.frame(dataset = names(paths), view = "benchmark_analysis", path = unname(paths),
  md5 = vapply(paths, md5, character(1)), size = as.character(file.info(paths)$size), proof = paste0(unname(paths), ".md5"))
refs_path <- file.path(root, "inputs.tsv")
utils::write.table(refs, refs_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
args <- list(analysis = "res50", scope = "benchmark_union", datasets = "D,E", config_path = config_path,
  output_dir = file.path(root, "derived/res50"), run_id = "test_res50", input_references = refs_path, input_references_sha256 = sha(refs_path))
protected <- tools::md5sum(c(paths, paste0(paths, ".md5"), refs_path, config_path))
load_h5ad_counts_free <- function(...) stop("NO EXPRESSION/EMBEDDING LOADER")
ecoda_subsetting_plan <- function(...) stop("NO CELL-DEPTH DRAWS")
options(mc.cores = 7L)
first <- .ecoda_run(args)
stopifnot(getOption("mc.cores") == 7L, first$status == "COMPLETED", length(first$artifacts) == 2L,
  first$source_mode == "bound_saved_cluster_obs", !"Seurat" %in% loadedNamespaces())
options(mc.cores = 1L)
res <- readRDS(file.path(args$output_dir, "D_ECODA_seuratres_50.rds"))
stopifnot(identical(rownames(res$counts), sample_ids), identical(names(res$labels), sample_ids),
  identical(colnames(res$counts), c("0", "1", "2")), !"UNUSED" %in% colnames(res$counts),
  res$source_identity$dropped_annotation_cells == 1L, res$source_identity$n_obs == length(sample),
  identical(res$scoring_policy, "calc_sep_score_five_v1_lisi_separation"),
  all(is.finite(unlist(res$scores[.ecoda_screen_metrics()]))), !is.null(res$source_identity$stat$proof_md5))
# The production count/CLR core agrees, with no global padded category.
raw <- res$counts - .5
stopifnot(all(as.matrix(raw) >= 0), all(as.matrix(raw) == floor(as.matrix(raw))))
expected <- .ecoda_counts_bundle(raw, res$labels, sample_ids, "parity")
stopifnot(isTRUE(all.equal(res$feat_mat, expected$feat_mat)), identical(res$dist_mat, expected$dist_mat),
  isTRUE(all.equal(unlist(res$scores[.ecoda_screen_metrics()]), unlist(expected$scores[.ecoda_screen_metrics()]))))
padded <- cbind(raw, GLOBAL_PADDING = 0)
padded <- padded[, colSums(padded) > 0, drop = FALSE]
stopifnot(identical(padded, raw))
# Exact target recovery selects only its declared cohort; other refs are not opened.
second_args <- args; second_args$analysis <- "harmony"; second_args$datasets <- "E"
second_args$output_dir <- file.path(root, "derived/harmony"); second_args$run_id <- "test_harmony"
second <- .ecoda_run(second_args)
stopifnot(length(second$artifacts) == 1L, second$artifacts[[1]]$dataset == "E", second$artifacts[[1]]$method == "ECODA_seuratres_2_harmony")
sources <- .ecoda_bound_sources(args, config, c("D", "E"), "res50")
input <- .ecoda_bound_counts(sources$D, project)
again <- .ecoda_bound_counts(sources$D, project)
stopifnot(identical(input$counts, again$counts), identical(input$count_fingerprint, again$count_fingerprint),
  identical(protected, tools::md5sum(names(protected))))
# Biology changes evaluation only, never the saved-label tally or CLR features.
dir.create(file.path(root, "permuted")); perm_path <- file.path(root, "permuted/D.h5ad")
permuted_bio <- unname(setNames(rep(c("a", "b"), 6), sample_ids)[sample])
reticulate::py$make_derived_obs(perm_path, "donor", sample, permuted_bio, as.list(cluster50), as.list(cluster2))
writeLines(c(paste0("MD5=", md5(perm_path)), paste0("SIZE=", file.info(perm_path)$size), paste0("PATH=", perm_path)), paste0(perm_path, ".md5"))
permuted_source <- sources$D; permuted_source$h5ad_path <- perm_path; permuted_source$h5ad_md5 <- md5(perm_path)
permuted_source$h5ad_size <- as.numeric(file.info(perm_path)$size); permuted_source$input_proof <- paste0(perm_path, ".md5")
permuted_source$stat <- .ecoda_bound_input_stat(permuted_source)
permuted <- .ecoda_bound_counts(permuted_source, project)
permuted_bundle <- .ecoda_counts_bundle(permuted$counts, permuted$labels, permuted$sample_ids, "permuted")
stopifnot(identical(permuted$counts, input$counts), !identical(permuted$labels, input$labels),
  identical(permuted_bundle$feat_mat, res$feat_mat), identical(permuted_bundle$dist_mat, res$dist_mat))
bad_manifest <- first; bad_manifest$sources[[1]]$sample_col <- NULL
reject(ecoda_validate_derived_manifest(bad_manifest))
bad_manifest <- first; bad_manifest$input_references_sha256 <- NULL
reject(ecoda_validate_derived_manifest(bad_manifest))
bad <- args; bad$input_references_sha256 <- strrep("0", 64); reject(.ecoda_bound_sources(bad, config, "D", "res50"))
bad <- args; bad$input_references_sha256 <- NULL; reject(.ecoda_bound_sources(bad, config, "D", "res50"))
foreign <- refs; foreign$dataset[2] <- "DISABLED"
foreign_path <- file.path(root, "foreign.tsv")
utils::write.table(foreign, foreign_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
bad <- args; bad$input_references <- foreign_path; bad$input_references_sha256 <- sha(foreign_path)
reject(.ecoda_bound_sources(bad, config, "D", "res50"))
reject(ecoda_select_screen_datasets(ecoda_raw_config_union(config), "DISABLED"))
bad_config <- config; bad_config$D$views$benchmark_analysis$columns <- list(sample = "WRONG")
reject(.ecoda_bound_sources(args, bad_config, "D", "res50"))
bad_config <- config; bad_config$D$columns$sample <- NULL; reject(.ecoda_bound_sources(args, bad_config, "D", "res50"))
bad_source <- sources$D; bad_source$h5ad_md5 <- strrep("0", 32); reject(.ecoda_bound_input_stat(bad_source))
alias <- file.path(root, "proof.alias"); stopifnot(file.symlink(sources$D$input_proof, alias))
bad_source <- sources$D; bad_source$input_proof <- alias; reject(.ecoda_bound_input_stat(bad_source))
bad_source <- sources$D; bad_source$sample_col <- "Sample"; reject(.ecoda_bound_counts(bad_source, project))
bad_source <- sources$D; bad_source$dataset <- "Stephenson"; reject(.ecoda_bound_counts(bad_source, project))
bad <- args; bad$input_dir <- root; reject(.ecoda_run(bad))
bad <- args; bad$datasets <- NULL; reject(.ecoda_run(bad))
bad <- args; bad$analysis <- "cell_subsetting"; reject(.ecoda_run(bad))
bad <- args; bad$input_references <- bad$input_references_sha256 <- NULL; bad$input_dir <- root
reject(.ecoda_run(bad))
reject(.ecoda_run(args))
stopifnot(identical(protected, tools::md5sum(names(protected))))
reticulate::py_run_string("h5py.Dataset.__getitem__=_old_derived_getitem")
unlink(root, recursive = TRUE)
cat("derived bound obs inputs: PASS; configured columns, active categories, exact refs, five-score parity, subset recovery, no expression/models/depth draws\n")
