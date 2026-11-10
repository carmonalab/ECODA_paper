# Actual maintained MDS chunks, tiny saved distances; no scoring or production data.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
source('src/utils/constants.R');source('src/utils/plotting.R');source('src/notebooks/benchmark_analysis.R')
run<-function(){
 root<-tempfile('mds-display-');dir.create(root);on.exit(unlink(root,recursive=TRUE))
 prior<-if(file.exists('Rplots.pdf'))unname(tools::md5sum('Rplots.pdf'))else NULL
 pdf(file.path(root,'scratch.pdf'));on.exit({dev.off();stopifnot(if(is.null(prior))!file.exists('Rplots.pdf')else identical(prior,unname(tools::md5sum('Rplots.pdf'))))},add=TRUE)
 lines<-readLines('notebooks/benchmark_analysis.rmd');starts<-grep('^```\\{r',lines);chunk<-function(i){first<-starts[i];last<-first+which(lines[(first+1):length(lines)]=='```')[1];parse(text=lines[(first+1):(last-1)])}
 ids<-paste0('s',1:6);x<-matrix(c(1,2,5,4,3,6,2,5,1,6,3,4),ncol=2,dimnames=list(ids,c('x','y')));d<-dist(x);coordinates<-cmdscale(d,k=2,eig=TRUE)$points
 e<-new.env(parent=.GlobalEnv);e$result_list<-list(bmark=list())
 for(ds in c('Bassez','Adams')){
  labels<-setNames(factor(rep(if(ds=='Bassez')c('E','NE')else c('a','b'),each=3)),ids)
  e$result_list$bmark[[ds]]<-setNames(lapply(unique(c(benchmark_main_method_keys,'Freq_highres')),function(m)list(dist_mat=d,labels=labels)),unique(c(benchmark_main_method_keys,'Freq_highres')))
 }
 # Only the display map is overridden: scientific config and named labels stay literal.
 config_before<-unname(tools::md5sum('datasets.json'))
 e$dataset_display_map<-c(Adams='Adams (pulmonary fibrosis)',Bassez='Bassez (ICB treatment response)')
 map_before<-e$dataset_display_map;source_before<-digest::digest(e$result_list)
 display_expr<-Filter(function(x)identical(paste(deparse(x),collapse=' '),'dataset_display_map[["Bassez"]] <- "Bassez (Clonal expansion)"'),as.list(chunk(1L)))
 stopifnot(length(display_expr)==1L);eval(display_expr[[1L]],e)
 stopifnot(identical(e$dataset_display_map[['Bassez']],'Bassez (Clonal expansion)'),
   identical(e$dataset_display_map['Adams'],map_before['Adams']),identical(source_before,digest::digest(e$result_list)),
   identical(config_before,unname(tools::md5sum('datasets.json'))))
 before<-digest::digest(e$result_list);saved<-list();calls<-0L
 e$methods_in_scope_named<-function(m)m;e$path_plots<-root;e$df_benchmark<-data.frame()
 e$benchmark_main_normalization_methods<-benchmark_main_method_keys;e$benchmark_extended_normalization_methods<-benchmark_extended_method_keys
 e$expected_pool<-benchmark_main_method_keys
 e$benchmark_normalized_mean<-function(score_frame,dataset,method,normalization_methods){stopifnot(identical(normalization_methods,e$expected_pool),method%in%normalization_methods);calls<<-calls+1L;.5}
 e$ggsave<-function(filename,...){saved[[basename(filename)]]<<-e$plot_list;invisible(filename)}
 # Any fresh metric calculation would violate this display-only fixture.
 e$calc_sep_score<-e$clust_eval<-e$calc_modularity<-function(...)stop('SCORING FORBIDDEN')
 eval(chunk(14L),e)
 check<-function(plots,ds){
  for(p in plots){
   stopifnot(identical(p$data$Dim1,unname(coordinates[,1])),identical(p$data$Dim2,unname(coordinates[,2])))
   if(ds=='Bassez')stopifnot(identical(as.character(p$data$Cluster),rep(c('Expander','Non-expander'),each=3)),identical(p$labels$colour,'Clonal expansion'),identical(p$labels$shape,'Clonal expansion'))
   else stopifnot(identical(as.character(p$data$Cluster),rep(c('a','b'),each=3)))
  }
 }
 stopifnot(length(saved)==2L,calls==22L);check(saved[[grep('_Bassez.pdf$',names(saved),value=TRUE)]],'Bassez');check(saved[[grep('_Adams.pdf$',names(saved),value=TRUE)]],'Adams')
 # Figure2B intentionally contains frequency and therefore retains extended21.
 e$expected_pool<-benchmark_extended_method_keys;eval(chunk(42L),e);check(saved[['Figure_2_B_MDS_Bassez.pdf']],'Bassez')
 stopifnot(calls==30L,identical(digest::digest(e$result_list),before))
 cat('PASS: actual SuppMDS expansion labels/legend, main12; Figure2B frequency/extended21; unchanged coordinates/source labels; no scoring\n')
}
run()
