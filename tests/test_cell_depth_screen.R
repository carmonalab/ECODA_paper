suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2)})
source("src/utils/constants.R")
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/math_utils.R")
source("src/utils/scoring_metrics.R")
source("src/5_run_benchmark_methods/run_local_ecoda_derived.R")
source("src/notebooks/benchmark_analysis.R")
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
root <- tempfile("cell-depth-"); dir.create(root); root <- normalizePath(root)
md5 <- function(path) unname(tools::md5sum(path))
scratch_plot <- normalizePath("Rplots.pdf", mustWork=FALSE)
scratch_md5 <- if (file.exists(scratch_plot)) md5(scratch_plot) else NULL
evidence <- file.path(root, "review.txt"); writeLines("Synthetic approved compact counts/current input review", evidence)
metrics <- evaluation_metrics(); stopifnot(identical(metrics, .ecoda_screen_metrics()))
raw_config <- jsonlite::fromJSON("datasets.json", simplifyVector = FALSE)
stopifnot(identical(ecoda_select_screen_datasets(c("Bassez", "Stephenson"), "Stephenson"), "Stephenson"))
fail(.ecoda_parse_flags(c("--dataset", "Bassez")))
fail(ecoda_select_screen_datasets(c("Bassez"), "Bassez,"))
fail(ecoda_select_screen_datasets(c("Bassez"), "Bassez,Bassez")); fail(ecoda_select_screen_datasets(c("Bassez"), "Other"))
fixture <- function(dataset, n, types, total, tag = dataset) {
  ids <- paste0(dataset, seq_len(n)); counts <- matrix(2, n, types, dimnames = list(ids, paste0("CT", seq_len(types))))
  counts[1,1] <- counts[1,1] + total - sum(counts)
  labels <- factor(rep(c("a", "b"), length.out = n)); names(labels) <- ids
  hr <- file.path(root, paste0(tag, "_HR.rds")); meta <- file.path(root, paste0(tag, "_metadata.rds"))
  saveRDS(list(counts = counts + .5, labels = labels), hr)
  saveRDS(list(labels = labels, cells_per_sample = setNames(rowSums(counts), ids),
    n_samples = n, n_cells = total, n_cell_types_high_res = types), meta)
  list(hr = hr, metadata = meta, counts = counts, labels = labels,
    contract = list(sample_ids = ids, n_samples = n, n_cell_types = types, total_cells = total, annotated_cells = total))
}
bassez <- fixture("Bassez", 26, 40, 74948)
steph <- fixture("Stephenson", 53, 4, 10600)
make_control <- function(fixtures, tag) {
  ledger <- scope <- list(); index <- 0L
  for (dataset in names(fixtures)) for (method in c("ECODA_authors_HR", "__dataset_metadata__")) {
    index <- index + 1L; f <- fixtures[[dataset]]
    path <- if (method == "ECODA_authors_HR") f$hr else f$metadata
    row <- data.frame(dataset = dataset, view = "benchmark_analysis", method_key = method,
      kind = "standalone_result", artifact_path = path, bundle_key = "", run_root = root,
      md5 = md5(path), receipt_path = file.path(root, paste0(tag, index, ".receipt.rds")), receipt_md5 = "")
    receipt <- c(as.list(row[c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")]),
      list(schema = "ecoda_analysis_release_receipt_v1", status = "reviewed", origin = "new", input_md5 = strrep("a",32),
        source_policy = "current-counts-v1", evidence = "reviewed fixture", evidence_path = evidence, evidence_md5 = md5(evidence)))
    saveRDS(receipt, row$receipt_path); row$receipt_md5 <- md5(row$receipt_path)
    expected <- row[c("dataset", "view", "method_key")]; expected$require_new <- TRUE
    expected$current_input_md5 <- strrep("a",32); expected$source_policy <- "current-counts-v1"; expected$status <- "approved"
    ledger[[index]] <- row; scope[[index]] <- expected
  }
  control <- list(schema = "ecoda_cell_depth_inputs_v1", status = "approved", ledger = do.call(rbind, ledger),
    expected_scope = do.call(rbind, scope), cohort_contract = lapply(fixtures, `[[`, "contract"))
  path <- file.path(root, paste0(tag,".rds")); saveRDS(control,path); path
}
input <- make_control(list(Bassez = bassez, Stephenson = steph), "inputs")
configs <- ecoda_screen_release_sources(input, md5(input), c("Bassez", "Stephenson"), raw_config)
fail(ecoda_screen_release_sources(input, md5(input), "Missing", raw_config))
old <- fixture("Bassez",26,33,24139,"oldBassez")
old_input <- make_control(list(Bassez=old),"old-input")
fail(ecoda_screen_release_sources(old_input,md5(old_input),"Bassez",raw_config))
old <- fixture("Stephenson",55,4,11000,"oldSteph")
old_input <- make_control(list(Stephenson=old),"old-steph-input")
fail(ecoda_screen_release_sources(old_input,md5(old_input),"Stephenson",raw_config))
# Actual five-metric scorer parity for one compact draw; subsequent schedule
# checks stub only scoring, keeping real count sampling and CLR distances.
core_counts <- bassez$counts[1:12,1:6]; core_counts[,2:6] <- rep(1:5,each=12)
core_labels <- bassez$labels[1:12]
core_features <- clr(as.data.frame(core_counts + .5))
core <- calc_sep_score(dist(core_features),core_labels,parallel=1)
stopifnot(all(is.finite(unlist(core[metrics]))))
create_result_bundle <- function(feat_mat,labels,dist_mat,extra) {
  runif(if (as.character(labels[1]) == "a") 3L else 7L) # Deliberately label-dependent scoring RNG.
  c(list(feat_mat=feat_mat,labels=labels,dist_mat=dist_mat, scores=setNames(as.list(seq(.1,.5,.1)),metrics)),extra)
}
original <- tools::md5sum(c(bassez$hr,bassez$metadata,steph$hr,steph$metadata))
result <- .ecoda_process_local_dataset("Bassez",configs$Bassez,101:120)
steph_result <- .ecoda_process_local_dataset("Stephenson",configs$Stephenson,101:120)
stopifnot(nrow(steph_result$results)==181L, all(vapply(steph_result$draws,function(x)nrow(x$counts)==53L,logical(1))))
stopifnot(nrow(result$results)==181L, length(result$draws)==181L,
          identical(result$results$ANOSIM,result$results$anosim_score))
source_row <- c(list(dataset="Bassez"),configs$Bassez)
package <- .ecoda_screen_package(result$results,result$draws,list(source_row),strrep("b",32),strrep("c",32),strrep("d",32),"compact_counts_sampling")
bad <- result$results; bad$sil_score[1] <- NA_real_; fail(.ecoda_validate_cell_subsetting_results(bad))
bad_draws <- result$draws; bad_draws[[1]]$counts[1,1] <- bad_draws[[1]]$counts[1,1]+1
fail(.ecoda_validate_screen_draws(bad_draws,result$results))
# Biological label permutation cannot change draws or the selected features.
permuted <- bassez; labs <- bassez$labels[rev(seq_along(bassez$labels))];names(labs)<-names(bassez$labels)
x<-readRDS(bassez$hr);x$labels<-labs;permuted$hr<-file.path(root,"permutedHR.rds");saveRDS(x,permuted$hr)
x<-readRDS(bassez$metadata);x$labels<-labs;permuted$metadata<-file.path(root,"permutedMeta.rds");saveRDS(x,permuted$metadata)
permuted_input<-make_control(list(Bassez=permuted),"permuted")
permuted_config<-ecoda_screen_release_sources(permuted_input,md5(permuted_input),"Bassez",raw_config)$Bassez
other <- .ecoda_process_local_dataset("Bassez",permuted_config,101:120)
stopifnot(!identical(result$draws[[1]]$labels, other$draws[[1]]$labels))
stopifnot(all(vapply(seq_len(181),function(i) identical(result$draws[[i]]$counts,other$draws[[i]]$counts) &&
  identical(result$draws[[i]]$features_md5,other$draws[[i]]$features_md5),logical(1))))
# Saved draw rescoring cannot invoke or alter the sampler.
sampler <- .ecoda_local_sample_counts
.ecoda_local_sample_counts <- function(...) stop("Sampler must not run during saved-draw rescoring")
rescored <- ecoda_rescore_screen_draws(package)
stopifnot(identical(rescored$draws,package$draws),identical(rescored$results,package$results))
.ecoda_local_sample_counts <- sampler
# Exercise the real entry orchestration with scoring-only imports stubbed;
# one selected dataset runs serially, without Slurm, sockets or model imports.
runner_path <- normalizePath("src/5_run_benchmark_methods/run_local_ecoda_derived.R")
.ecoda_runner_path <- function() runner_path
real_source <- source
source <- function(file, ...) {
  if (!basename(file) %in% c("analysis_release.R", "math_utils.R", "scoring_metrics.R", "benchmark_methods_r.R"))
    stop("Unexpected entrypoint import: ", file)
  invisible(NULL)
}
Sys.setenv(ECODA_DERIVED_ROOT = root)
manifest <- .ecoda_run(list(analysis="cell_subsetting", scope="benchmark_union", datasets="Bassez", cores="1",
  config_path=normalizePath("datasets.json"), analysis_release=input, analysis_release_md5=md5(input),
  output_dir=file.path(root,"entrypoint"), run_id="entrypoint"))
source <- real_source
validated_entry <- ecoda_validate_screen_output(list(validate_output=file.path(root,"entrypoint"),
  config_path=normalizePath("datasets.json"), analysis_release=input, analysis_release_md5=md5(input),
  datasets="Bassez", scope="benchmark_union"))
stopifnot(nrow(validated_entry$results)==181L)
stopifnot(manifest$status == "COMPLETED", identical(manifest$selected_datasets,"Bassez"), length(manifest$artifacts)==2L)
entry <- readRDS(file.path(root,"entrypoint","ECODA_authors_HR_cell_subsetting.rds"))
stopifnot(entry$schema=="ecoda_cell_depth_five_v1", length(entry$draws)==181L,
          identical(entry$current_inputs$Bassez,strrep("a",32)))
# Portable reviewed output, with explicitly relocated sibling CSV identity.
csv<-file.path(root,"screen.csv");write.csv(package$results,csv,row.names=FALSE)
package$csv<-list(basename=basename(csv),md5=md5(csv));package$source_release_md5<-md5(input)
artifact<-file.path(root,"screen.rds");saveRDS(package,artifact)
row<-data.frame(dataset="ALL",view="benchmark_analysis",method_key="ECODA_authors_HR_cell_subsetting",kind="standalone_result",
 artifact_path=artifact,bundle_key="",run_root=root,md5=md5(artifact),receipt_path=file.path(root,"output.receipt.rds"),receipt_md5="")
receipt<-c(as.list(row[c("dataset","view","method_key","kind","artifact_path","bundle_key","run_root","md5")]),
 list(schema="ecoda_analysis_release_receipt_v1",status="reviewed",origin="new",metadata_schema="ecoda_cell_depth_five_v1",input_md5=md5(input),source_policy="screen-v1",
 evidence="reviewed output",evidence_path=evidence,evidence_md5=md5(evidence)))
saveRDS(receipt,row$receipt_path);row$receipt_md5<-md5(row$receipt_path)
scope<-row[c("dataset","view","method_key")];scope$require_new<-TRUE;scope$current_input_md5<-md5(input);scope$source_policy<-"screen-v1";scope$status<-"approved"
output_control<-file.path(root,"output.rds");saveRDS(list(ledger=row,expected_scope=scope),output_control)
read_current<-function() load_derived_cell_subsetting(NULL,NULL,"Bassez",list(Bassez=names(bassez$labels)),
 release_control_path=output_control,release_control_md5=md5(output_control),
 expected_labels_by_dataset=list(Bassez=bassez$labels),expected_current_inputs=list(Bassez=strrep("a",32)))
loaded<-read_current();stopifnot(nrow(loaded$results)==181L, all(metrics%in%names(loaded$results)))
# Execute only the maintained cell-depth plotting chunk, intercepting save;
# no production plot or notebook rendering occurs.
lines<-readLines("notebooks/benchmark_analysis.rmd")
start<-grep("# Supp fig X: ECODA author annotation sensitivity to cell depth",lines,fixed=TRUE)
finish<-grep("# Exploratory Stephenson UCell activity versus ECODA",lines,fixed=TRUE)
section<-lines[start:(finish-1)]
first<-grep("^```\\{r",section)[1];last<-tail(grep("^```$",section),1)
cell_subsetting_requested<-TRUE;cell_subsetting_results<-loaded;path_plots<-root
# Reporting fixture only: all 11 cohorts/10 targets/5 metrics, no new draws or scores.
display_datasets<-names(Filter(function(x)isTRUE(x$use_for_benchmark)&&!is.null(x$views$benchmark_analysis),raw_config))
display_datasets<-display_datasets[!grepl("^_",display_datasets)]
stopifnot(length(display_datasets)==11L)
cell_subsetting_results$results<-do.call(rbind,lapply(seq_along(display_datasets),function(i){
  frame<-loaded$results;frame$dataset<-display_datasets[i]
  for(k in seq_along(metrics))frame[[metrics[k]]]<-i/20+as.integer(factor(frame$target))/100+k/1000+frame$replicate/10000
  frame
}))
real_ggsave<-ggsave;ggsave<-function(filename,plot,...)invisible(filename)
grDevices::pdf(file.path(root,"cell-depth-fixture.pdf"))
tryCatch(eval(parse(text=section[(first+1):(last-1)])),finally={grDevices::dev.off();ggsave<-real_ggsave})
built<-suppressWarnings(ggplot_build(p_cell_subsetting))
stopifnot(length(unique(built$layout$layout$PANEL))==5L)
stopifnot(nrow(cell_subsetting_dataset_means)==11L*10L*5L,nrow(cell_subsetting_aggregate)==50L,
 all(is.finite(cell_subsetting_aggregate$se)),all(cell_subsetting_aggregate$se>0),
 nrow(built$data[[1L]])==50L,all(is.finite(built$data[[1L]]$ymin)),all(is.finite(built$data[[1L]]$ymax)))
display_long<-pivot_longer(cell_subsetting_results$results,all_of(metrics),names_to="metric",values_to="value")
expected_means<-aggregate(value~dataset+target+metric,display_long,mean)
expected_summary<-do.call(rbind,lapply(split(expected_means,paste(expected_means$target,expected_means$metric)),function(x){
  stopifnot(nrow(x)==11L,length(unique(x$dataset))==11L)
  data.frame(target=x$target[1],metric=x$metric[1],expected_value=mean(x$value),expected_se=sd(x$value)/sqrt(11L))
}))
comparison<-merge(as.data.frame(cell_subsetting_aggregate),expected_summary,by=c("target","metric"))
stopifnot(nrow(comparison)==50L,max(abs(comparison$value-comparison$expected_value))<1e-14,
 max(abs(comparison$se-comparison$expected_se))<1e-14)
writeLines("corrupt",csv);fail(read_current())
stopifnot(identical(tools::md5sum(names(original)),original))
# Source publication remains create-only.
.ecoda_save_rds(package,file.path(root,"create-only.rds"));before<-md5(file.path(root,"create-only.rds"))
fail(.ecoda_save_rds(list(corrupt=TRUE),file.path(root,"create-only.rds")));stopifnot(md5(file.path(root,"create-only.rds"))==before)
if (is.null(scratch_md5)) stopifnot(!file.exists(scratch_plot)) else stopifnot(identical(md5(scratch_plot),scratch_md5))
unlink(root,recursive=TRUE)
cat("cell-depth screen: five metrics, immutable compact draws/current inputs and portable release: OK\n")
