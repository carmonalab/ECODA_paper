suppressPackageStartupMessages({library(dplyr);library(survival)})
source('src/utils/analysis_release.R');source('src/5_run_benchmark_methods/run_local_ecoda_derived.R')
source('src/5_run_benchmark_methods/benchmark_pipeline.R');source('src/utils/math_utils.R');source('src/utils/scoring_metrics.R')
source('src/notebooks/run_zeroimp_count_screen.R')
fail <- function(x)stopifnot(inherits(tryCatch(force(x),error=identity),'error'))
root <- Sys.getenv('ZEROIMP_TEST_ROOT',unset=tempfile('zeroimp-count-screen-'));dir.create(root,showWarnings=FALSE);root <- normalizePath(root)
md5 <- function(p)unname(tools::md5sum(p))
ids <- paste0('s',1:12);counts <- matrix(as.double(1:48),12,dimnames=list(ids,paste0('CT',1:4)));counts[c(1,18,27,40)] <- 0
labels <- setNames(factor(rep(c('A','B'),each=6)),ids)
input <- list(counts=counts,labels=labels)
stopifnot(length(zeroimp_menu())==43L,!anyDuplicated(zeroimp_menu()),identical(zeroimp_policy()$metrics,.ecoda_screen_metrics()),zeroimp_policy()$parent_counts_offset==.5)
# Actual ordinary recipe body, including nested multLN and any native dropped
# settings. Scoring is isolated; label-dependent RNG consumption cannot alter a
# subsequent recipe's distance matrix or imputation draw.
metrics <- zeroimp_policy()$metrics;seen <- list()
stub <- function(dist_mat,labels,parallel=1) {
 stopifnot(parallel==1);seen[[length(seen)+1L]] <<- as.matrix(dist_mat)
 runif(if(as.character(labels[1])=='A')3L else 7L)
 setNames(as.list(rep(.5,8)),c(metrics,'mod_knnsqrtn_score','mod_knn6_score','mod_knn9_score'))
}
set.seed(701);before <- .Random.seed;kind <- RNGkind()
a <- zeroimp_run_cohort(input,score_fun=stub);geometry <- seen
stopifnot(identical(before,.Random.seed),identical(kind,RNGkind()),nrow(a$results)==43L,
 identical(a$results$setting,zeroimp_menu()),length(geometry)>0L,a$recipes$counts_all_0.5$status=="PASS",all(a$results$status%in%c('PASS','FAILED_OR_DROPPED_BY_ORDINARY_PRODUCER')))
seen <- list();permuted <- input;permuted$labels <- setNames(rev(labels),ids)
b <- zeroimp_run_cohort(permuted,score_fun=stub)
stopifnot(identical(geometry,seen),identical(a$counts,b$counts),!identical(a$labels,b$labels),identical(a$counts_md5,b$counts_md5),!identical(a$labels_md5,b$labels_md5))
named_input <- input;rownames(named_input$counts) <- setNames(ids,ids)
named_result <- zeroimp_run_cohort(named_input,score_fun=stub)
stopifnot(named_result$recipes$counts_all_0.5$status=='PASS')
# Compare wrapper scores against the actual current scorer, with its native
# default k3/ARI/silhouette/ANOSIM/LISI definitions. No fake metric backfill.
recipe <- function(df,labels)list(counts_all_0.5=calc_sep_score(dist(clr(df+.5)),labels))
RNGkind("L'Ecuyer-CMRG","Inversion","Rejection");set.seed(123)
expected <- calc_sep_score(dist(clr(as.data.frame(counts)+.5)),labels,parallel=1)
actual <- zeroimp_run_cohort(input,recipe_fun=recipe)
stopifnot(isTRUE(all.equal(actual$recipes$counts_all_0.5$scores,expected)),sum(actual$results$status=='PASS')==1L,
 all(is.na(actual$results[actual$results$status!='PASS',metrics])))
fail(zeroimp_run_cohort(list(counts=counts+.5,labels=labels),score_fun=stub))
# A multLN-like sample drop is explicit and never scored/approved. Geometry
# is retained only for PASS settings with the exact original sample universe.
drop_recipe <- function(df,labels)list(counts_all_0.5=calc_sep_score(dist(clr(df[-1,,drop=FALSE]+.5)),labels[-1]))
never <- function(...)stop('partial samples must not reach scoring')
dropped <- zeroimp_run_cohort(input,score_fun=never,recipe_fun=drop_recipe)
stopifnot(dropped$recipes$counts_all_0.5$status=='FAILED_OR_DROPPED_BY_ORDINARY_PRODUCER',identical(dropped$recipes$counts_all_0.5$dropped_ids,ids[1]),is.null(dropped$recipes$counts_all_0.5$dist_mat),is.null(dropped$recipes$counts_all_0.5$scores),identical(a$rng_after,b$rng_after))
# Checked parent release fixture: real integer+.5 recovery, metadata truth,
# receipts and exact eleven inventory and original nine lineage. Legacy origins are never rewritten.
evidence <- file.path(root,'review.txt');writeLines('ROOT fixture: approved legacy209 current HR/metadata counts',evidence)
ledger <- scope <- parents <- list();contracts <- list();n <- 0L
for(ds in zeroimp_datasets()) {
 current_counts <- counts;current_labels <- labels;current_ids <- ids
 if(ds%in%c('Bassez','Stephenson')) {
  n_samples <- if(ds=='Bassez')26L else 53L;types <- if(ds=='Bassez')40L else 4L;total <- if(ds=='Bassez')74948 else 10600
  current_ids <- paste0(ds,seq_len(n_samples));current_counts <- matrix(2,n_samples,types,dimnames=list(current_ids,paste0('CT',seq_len(types))));current_counts[1,1] <- current_counts[1,1]+total-sum(current_counts);current_labels <- setNames(factor(rep(c('A','B'),length.out=n_samples)),current_ids)
 }
 hr <- file.path(root,paste0(ds,'_composition.rds'));meta <- file.path(root,paste0(ds,'_metadata.rds'))
 saveRDS(list(ECODA_authors_HR=list(counts=current_counts+.5,labels=current_labels)),hr)
 saveRDS(list(labels=current_labels,cells_per_sample=setNames(rowSums(current_counts),current_ids)),meta)
 parents[[ds]] <- data.frame(dataset=ds,author_path=hr,author_md5=md5(hr),author_bundle_key='ECODA_authors_HR',metadata_path=meta,metadata_md5=md5(meta),current_input_path=file.path(root,paste0(ds,'.ABSENT.h5ad')),current_input_md5=strrep('a',32),expected_counts_policy='production get_ct_comp_df valid labels; counts_all +0.5; compare raw/imputed/CLR without scoring')
 contracts[[ds]] <- list(sample_ids=current_ids,n_samples=nrow(current_counts),n_cell_types=ncol(current_counts),total_cells=sum(current_counts),annotated_cells=sum(current_counts))
 for(key in c('ECODA_authors_HR','__dataset_metadata__')) {
  n <- n+1L;p <- if(key=='ECODA_authors_HR')hr else meta
  row <- data.frame(dataset=ds,view='benchmark_analysis',method_key=key,kind=if(key=='ECODA_authors_HR')'rds_bundle' else 'standalone_result',artifact_path=p,bundle_key=if(key=='ECODA_authors_HR')key else '',run_root=root,md5=md5(p),receipt_path=file.path(root,paste0(n,'.receipt.rds')),receipt_md5='')
  receipt <- c(as.list(row[c('dataset','view','method_key','kind','artifact_path','bundle_key','run_root','md5')]),list(schema='ecoda_analysis_release_receipt_v1',status='reviewed',origin=if(ds%in%c('Bassez','Stephenson'))'new' else 'legacy_reuse',input_md5=strrep('a',32),source_policy='approved-legacy209-counts',evidence='synthetic approved counts',evidence_path=evidence,evidence_md5=md5(evidence)))
  saveRDS(receipt,row$receipt_path);row$receipt_md5 <- md5(row$receipt_path);ledger[[n]] <- row
  e <- row[c('dataset','view','method_key')];e$require_new <- ds%in%c("Bassez","Stephenson");e$current_input_md5 <- strrep('a',32);e$source_policy <- 'approved-legacy209-counts';e$status <- 'approved';scope[[n]] <- e
 }
}
parent_path <- file.path(root,'author_hr_metadata_pairs.11.tsv');write.table(do.call(rbind,parents),parent_path,sep='\t',quote=FALSE,row.names=FALSE)
legacy_path <- file.path(root,'author_hr_metadata_pairs.9.tsv');write.table(do.call(rbind,parents[setdiff(names(parents),c('Bassez','Stephenson'))]),legacy_path,sep='\t',quote=FALSE,row.names=FALSE)
control <- list(schema='ecoda_cell_depth_inputs_v1',status='approved',ledger=do.call(rbind,ledger),expected_scope=do.call(rbind,scope),cohort_contract=contracts,parent_inventory=list(path=parent_path,md5=md5(parent_path)),legacy_parent_inventory=list(path=legacy_path,md5=md5(legacy_path)))
control_path <- file.path(root,'release.rds');saveRDS(control,control_path)
if(identical(Sys.getenv('ZEROIMP_TEST_BUILD_ONLY'),'1'))quit(status=0L)
observed <- zeroimp_control(control_path,md5(control_path),normalizePath('datasets.json'),'Adams')
parent <- zeroimp_input(control_path,md5(control_path),'Adams',observed$raw)
stopifnot(identical(parent$input$counts,counts),identical(as.character(parent$input$labels),as.character(labels)))
stopifnot(nrow(zeroimp_input(control_path,md5(control_path),'Bassez',observed$raw)$input$counts)==26L,nrow(zeroimp_input(control_path,md5(control_path),'Stephenson',observed$raw)$input$counts)==53L)
fail(zeroimp_control(control_path,md5(control_path),normalizePath('datasets.json'),'Other'))
bad <- control;bad$status <- 'pending';bad_path <- file.path(root,'pending.rds');saveRDS(bad,bad_path)
fail(zeroimp_control(bad_path,md5(bad_path),normalizePath('datasets.json'),'Adams'))
source_row <- c(list(dataset='Adams'),parent$config)
a$dataset <- 'Adams';a$schema <- zeroimp_policy()$schema;a$source <- source_row;a$run_id <- 'fixture';a$config_md5 <- strrep('b',32);a$runner_md5 <- strrep('c',32);a$release_md5 <- md5(control_path)
zeroimp_validate_package(a,parent$input,source_row,'fixture',strrep('b',32),strrep('c',32),md5(control_path))
bad <- a;bad$results$setting[1] <- 'wrong';fail(zeroimp_validate_package(bad,parent$input,source_row,'fixture',strrep('b',32),strrep('c',32),md5(control_path)))
bad <- a;bad$counts[1,1] <- 1;fail(zeroimp_validate_package(bad,parent$input,source_row,'fixture',strrep('b',32),strrep('c',32),md5(control_path)))
bad <- a;passed <- which(bad$results$status=='PASS')[1];bad$recipes[[passed]]$dist_mat[1] <- 0;fail(zeroimp_validate_package(bad,parent$input,source_row,'fixture',strrep('b',32),strrep('c',32),md5(control_path)))
# Validation does not invoke the scorer/recipe, and source parents stay immutable.
original <- tools::md5sum(c(control$ledger$artifact_path,control$ledger$receipt_path,evidence,parent_path))
scorer <- calc_sep_score;calc_sep_score <- function(...)stop('scoring must not run during validation')
zeroimp_validate_package(a,parent$input,source_row,'fixture',strrep('b',32),strrep('c',32),md5(control_path));calc_sep_score <- scorer
bad <- a;bad$recipe_body_md5 <- 'bad';fail(zeroimp_validate_package(bad,parent$input,source_row,'fixture',strrep('b',32),strrep('c',32),md5(control_path)))
# Actual create-only entrypoint and no-score validation of a single selected
# cohort. Only compact category counts are read; current H5AD paths remain absent.
source_manifest <- file.path(root,'source.manifest');writeLines('immutable source fixture',source_manifest)
Sys.setenv(ECODA_SOURCE_ROOT=normalizePath('.'),ECODA_SOURCE_MANIFEST=source_manifest,ECODA_SOURCE_SNAPSHOT_REQUIRED='1',ECODA_RUN_ID='fixture')
base <- c('--project_root',normalizePath('.'),'--analysis_release',control_path,'--analysis_release_md5',md5(control_path),'--datasets','Adams','--config_path',normalizePath('datasets.json'),'--run_id','fixture')
# Controller binding may read only checked control/receipt metadata, never
# count-parent RDS before the shell's native ownership/producer guards.
bind_dir <- file.path(root,'bind');dir.create(bind_dir)
reader <- readRDS;loader <- zeroimp_input
readRDS <- function(file,...) { if(is.character(file)&&file%in%control$ledger$artifact_path)stop('core RDS read before native guards');reader(file,...) }
zeroimp_input <- function(...)stop('core parent loader forbidden during --bind')
tryCatch(zeroimp_main(c(base,'--bind',bind_dir)),finally={readRDS<-reader;zeroimp_input<-loader})
stopifnot(nrow(read.delim(file.path(bind_dir,'selection.tsv'),header=FALSE))==1L)
output <- file.path(root,'produced');zeroimp_main(c(base,'--produce',output))
produced <- readRDS(file.path(output,'Adams_zeroimp_counts.rds'));stopifnot(all(vapply(produced$source$parent_receipts,function(r)identical(r$origin,'legacy_reuse'),logical(1))))
calc_sep_score <- function(...)stop('scoring must not run during --validate')
zeroimp_main(c(base,'--validate',output));calc_sep_score <- scorer
fail(zeroimp_main(c(base,'--produce',output)))
stopifnot(!any(file.exists(vapply(parents,function(x)x$current_input_path,character(1)))))
stopifnot(identical(original,tools::md5sum(names(original))))
unlink(root,recursive=TRUE)
cat('zeroimp counts: current five metrics, raw +0.5 recovery, 43 explicit settings, isolated RNG, checked22 current/legacy parents and immutable9 lineage: OK\n')
