# Synthetic metadata/distances only; no production figures or RNA reads.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
run <- function() {
  device <- tempfile(fileext='.pdf');pdf(device);on.exit({dev.off();unlink(device)},add=TRUE)
  path <- 'notebooks/batch_effect_analysis_uncorrected_batchconfounding_contingency.rmd'
  text <- readLines(path,warn=FALSE)
  starts <- grep('^```\\{r ',text); chunks <- list()
  for(i in starts) {
    label <- sub('^```\\{r ([^,}]+).*','\\1',text[i]);end <- i+which(text[(i+1):length(text)]=='```')[1]
    body <- text[(i+1):(end-1)];parse(text=body);chunks[[label]]<-body
  }
  stopifnot(length(chunks)==6L,!any(grepl('readRDS|feat_mat|legacy_Kidney|input_root',unlist(chunks))))
  e <- new.env(parent=globalenv());sys.source('src/utils/batch_effect_analysis.R',e)
  eval(parse(text=chunks[['metadata-helpers']]),e);eval(parse(text=chunks[['plotting-helpers']]),e)
  wrapped <- e$.blank_temporal_plot('Unavailable', paste(rep('unavailable',30),collapse=' '))$labels$subtitle
  stopifnot(length(strsplit(wrapped,'\n',fixed=TRUE)[[1]])>1L,
    max(nchar(strsplit(wrapped,'\n',fixed=TRUE)[[1]]))<=64L)
  fail <- function(expr)stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
  ids <- paste0('s',1:12);labels <- c('control','control','mild','severe','mild','severe','mild','severe','mild','severe','mild','severe')
  metadata <- data.frame(Sample=ids,raw_sample=paste0('raw-',ids),severity=labels,batch=rep(c('a','b'),6),check.names=FALSE)
  metadata[['Sampling day (Days after symptom onset)']]<-c('control','control',1,5,10,15,20,20,25,25,30,30)
  e$datasets <- list(Covid19_PBMC=list(label_col='severity',batch_col='batch',evaluation_labels=list(primary='severity')))
  e$control <- list(ledger=data.frame(dataset='Covid19_PBMC',method_key='__sample_metadata__'),expected_scope=data.frame(dataset='Covid19_PBMC',method_key='__sample_metadata__'))
  e$load_analysis_release <- function(rows,scope,require_projection) {
    stopifnot(require_projection,nrow(rows)==1L,nrow(scope)==1L,rows$method_key=='__sample_metadata__')
    list(list(result=list(metadata=metadata),source=list(sample_id_col='Sample',artifact_path='/approved/thin/metadata.feather')))
  }
  info <- e$.read_sample_metadata('Covid19_PBMC');stopifnot(identical(info$sample_ids,ids),identical(info$metadata$raw_sample,metadata$raw_sample),info$source_status=='approved_projected_release')
  original <- metadata;metadata$Sample[2]<-metadata$Sample[1];fail(e$.read_sample_metadata('Covid19_PBMC'));metadata<-original
  metadata$batch<-NULL;fail(e$.read_sample_metadata('Covid19_PBMC'));metadata<-original
  features <- matrix(seq_len(48),12,dimnames=list(ids,paste0('g',1:4)));distances<-list(ECODA_authors_HR=list(dist_mat=dist(features)),Pseudobulk_hvg2000=list(dist_mat=dist(features*2)))
  e$load_batch_release_dataset <- function(control,dataset)list(sample_ids=ids,methods=distances)
  loaded <- e$.load_covid_temporal_distances(info);stopifnot(identical(loaded$ECODA,distances$ECODA_authors_HR$dist_mat))
  e$load_batch_release_dataset <- function(control,dataset)list(sample_ids=rev(ids),methods=distances);fail(e$.load_covid_temporal_distances(info))
  e$load_batch_release_dataset <- function(control,dataset)list(sample_ids=ids,methods=distances)
  bins <- e$.covid_time_bins(info);stopifnot(identical(unname(bins$summary$n_samples),c(6L,6L,0L)),identical(bins$summary$status,c('ESTIMABLE','ESTIMABLE','EMPTY_BIN')))
  for(mask in bins$masks[1:2])stopifnot(isTRUE(all.equal(as.matrix(dist(features[mask,,drop=FALSE])),as.matrix(loaded$ECODA)[ids[mask],ids[mask]],check.attributes=FALSE)))
  seen <- list();written <- character();e$plot_root<-tempdir();e$n_permutations<-1L;e$parallel_cores<-1L
  e$compute_batch_anosim <- function(dist_mat,grouping,...) {seen[[length(seen)+1L]]<<-dist_mat;list(statistic=.2,p_value=.5,status='ESTIMABLE',warning='')}
  e$plot_mds <- function(...) {stopifnot(list(...)$digits == 2);ggplot()+theme_void()}
  saved_plots <- list()
  e$.save_plot <- function(plot,path,...) {written<<-c(written,basename(path));saved_plots[[basename(path)]]<<-plot;invisible(path)}
  e$write_batch_analysis_table <- function(...)invisible(NULL)
  result <- e$.build_covid_temporal_plots(info)
  stopifnot(length(seen)==4L,length(written)==4L,nrow(result$anosim)==6L,all(result$anosim$status[5:6]=='EMPTY_BIN'),all(is.na(result$anosim$anosim_statistic[5:6])))
  stopifnot(isTRUE(all.equal(as.matrix(seen[[1]]),as.matrix(loaded$ECODA)[ids[1:6],ids[1:6]],check.attributes=FALSE)))
  summary_plot <- saved_plots[['Covid19_PBMC_time_bins_anosim.pdf']]
  text <- unlist(lapply(ggplot_build(summary_plot)$data,function(x)if('label'%in%names(x))as.character(x$label)else character()))
  stopifnot('Unavailable (n=0)'%in%text,all(c('0.20','NA')%in%text),
    !any(grepl('[0-9]\\.[0-9]{3}',text)),all(result$anosim$anosim_statistic[1:4]==.2))
  missing <- info;missing$metadata[['Sampling day (Days after symptom onset)']]<-NULL;seen<-list();written<-character();result <- e$.build_covid_temporal_plots(missing)
  stopifnot(length(seen)==0L,length(written)==4L,all(result$anosim$status=='MISSING_SAMPLING_DAY'),all(is.na(result$anosim$anosim_statistic)),all(result$bins$n_samples==0L))
  under <- info;under$metadata[['Sampling day (Days after symptom onset)']]<-c(1,2,rep(NA,10));under$metadata$severity<-rep('control',12)
  stopifnot(e$.covid_time_bins(under)$summary$status[1]=='UNDERPOWERED')
  # Execute the actual metadata-only routing loop with tiny eight-cohort tables.
  sys.source('src/utils/datasets_io.R',e)
  configured <- e$read_datasets_json(view='batch_effect_uncorrected')
  e$datasets <- configured[vapply(names(configured),function(ds)!startsWith(ds,'_')&&isTRUE(configured[[ds]]$use_for_batch_effect),logical(1))]
  e$batch_datasets <- names(e$datasets);stopifnot(length(e$batch_datasets)==8L)
  e$.read_sample_metadata <- function(dataset) {
    ds<-e$datasets[[dataset]];columns<-unique(c(unlist(ds$evaluation_labels),unlist(ds$batch_col)))
    m<-data.frame(Sample=ids[1:4],row.names=ids[1:4]);for(column in columns)m[[column]]<-rep(c('a','b'),2)
    list(dataset=dataset,metadata=m,sample_ids=ids[1:4],label_col=ds$label_col)
  }
  e$extended_candidates <- function(...)data.frame(covariate=character(),candidate_class=character())
  e$.build_plot <- e$.build_covid_extended_plot <- function(metadata_info,candidates,plot_type='extended')list(plot=ggplot(),summary=data.frame(dataset=metadata_info$dataset,biological_label=metadata_info$label_col,plot_type=plot_type,sample_count=4L,potential_confounding=FALSE))
  directory<-tempfile('contingency-route-');dir.create(directory);on.exit(unlink(directory,recursive=TRUE),add=TRUE);e$plot_root<-directory;written<-character()
  e$write_batch_analysis_table <- function(value,path)write.csv(value,path,row.names=FALSE)
  eval(parse(text=chunks[['generate-plots']]),e)
  stopifnot(length(written)==18L,all(c('Kidney_KPMP_full.pdf','Kidney_KPMP_full_extended.pdf','Kidney_KPMP_full_tissue.pdf','Kidney_KPMP_full_tissue_extended.pdf')%in%written),
    length(unique(paste(e$batch_confounding_summary$dataset,e$batch_confounding_summary$evaluation_id)))==9L,
    identical(names(e$metadata_results),'Covid19_PBMC'))
  eval(parse(text=chunks[['summary']]),e)
  kidney<-e$summary_overview[e$summary_overview$dataset=='Kidney_KPMP_full',,drop=FALSE]
  stopifnot(nrow(e$summary_overview)==18L,nrow(kidney)==4L,setequal(kidney$evaluation_id,c('condition','tissue')),all(kidney$sample_count==4L))
  cat('PASS: projected metadata only; exact8/nine-axis intake; saved-distance equivalence; empty/missing/underpowered temporal bins retain explicit unavailable outputs\n')
}
run()
