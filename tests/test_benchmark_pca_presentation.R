# Presentation-only synthetic PCA; no release/result/H5AD or scientific scoring.
suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
source('src/utils/plotting.R');source('src/notebooks/benchmark_analysis.R')
run <- function() {
  root<-Sys.getenv('PCA_PRESENTATION_FIXTURE_ROOT',unset='')
  if(!nzchar(root))root<-tempfile('pca-presentation-')
  stopifnot(!file.exists(root),dir.create(root,recursive=TRUE));root<-normalizePath(root)
  prior<-if(file.exists('Rplots.pdf'))unname(tools::md5sum('Rplots.pdf')) else NULL
  pdf(file.path(root,'scratch-device.pdf'));on.exit({dev.off();stopifnot(if(is.null(prior))!file.exists('Rplots.pdf') else identical(prior,unname(tools::md5sum('Rplots.pdf'))))},add=TRUE)
  set.seed(42);ids<-paste0('s',1:32);groups<-setNames(factor(rep(c('Negative','Positive'),each=16)),ids)
  x<-matrix(rnorm(32*8),32,dimnames=list(ids,c('KLRF1+ GZMB+ CD27- memory CD4 T cell','KLRF1+ GZMB+ CD27- EM CD8 T cell','KLRF1+ effector Vd1 gdT','GZMB+ Vd2 gdT','Adaptive NK cell','Naive CD4+ T cell','Follicular B cell','Alveolar Mac')))
  x[,1]<-x[,1]*3;x[,2]<-x[,2]*2
  genes<-x;colnames(genes)<-c('SFTPA1','SFTPA2','SFTPB','SCGB1A1','SFTPC','COL1A1','COL3A1','HBB')
  matrix_hash<-digest::digest(list(x,genes,groups))
  make<-function(value,title)plot_pca(value,groups,cluster_score=FALSE,mod_score=FALSE,anosim_score=FALSE,
    coord_equal=FALSE,n_ct_show=5,invisible='quali',repel=TRUE,title=title)
  original<-make(x,'Cell type composition (ECODA)');point<-which(vapply(original$layers,function(x)inherits(x$geom,'GeomPoint'),logical(1)))
  text<-which(vapply(original$layers,function(x)inherits(x$geom,'GeomTextRepel'),logical(1)))
  arrows<-which(vapply(original$layers,function(x)inherits(x$geom,'GeomSegment'),logical(1)))
  stopifnot(length(text)==1L,length(point)==1L,length(arrows)==1L)
  original_names<-original$layers[[text]]$data$name;original_params<-original$layers[[text]]$geom_params
  presented<-benchmark_pca_presentation(original)
  stopifnot(identical(original$layers[[text]]$data$name,original_names),identical(original$layers[[text]]$geom_params,original_params),
    identical(original$layers[[point]]$data,presented$layers[[point]]$data),identical(original$layers[[arrows]]$data,presented$layers[[arrows]]$data),
    identical(original$layers[[text]]$data[c('x','y')],presented$layers[[text]]$data[c('x','y')]),
    length(presented$layers[[text]]$data$name)==length(original_names),any(grepl('\n',presented$layers[[text]]$data$name)),
    identical(gsub('\n',' ',presented$layers[[text]]$data$name),as.character(original_names)),
    presented$layers[[text]]$geom_params$seed==123,presented$layers[[text]]$geom_params$force==3,
    presented$layers[[text]]$aes_params$size==3)
  built_text<-ggplot_build(presented)$data[[text]]
  stopifnot(identical(as.character(built_text$label),as.character(presented$layers[[text]]$data$name)),
    identical(built_text$x,original$layers[[text]]$data$x),identical(built_text$y,original$layers[[text]]$data$y))
  stopifnot(identical(ggplot_build(original)$data[[point]],ggplot_build(presented)$data[[point]]),
    identical(ggplot_build(original)$data[[arrows]],ggplot_build(presented)$data[[arrows]]))
  precision<-benchmark_pca_presentation(original+labs(x='PCA dim1 (31.4%)',y='PCA dim2 (9.76%)'))
  stopifnot(precision$labels$x=='PCA dim1 (31%)',precision$labels$y=='PCA dim2 (9.8%)')
  geneplot<-make(genes,'Pseudobulk')
  before<-patchwork::wrap_plots(list(geneplot,original),nrow=1,guides='collect')
  after<-patchwork::wrap_plots(list(benchmark_pca_presentation(geneplot),presented),nrow=1,guides='collect')
  ggsave(file.path(root,'presentation_before.pdf'),before,width=8,height=3.5)
  ggsave(file.path(root,'presentation_after.pdf'),after,width=8,height=3.5)
  stopifnot(identical(matrix_hash,digest::digest(list(x,genes,groups))))
  jsonlite::write_json(list(status='SYNTHETIC_PRESENTATION_PASS',feature_sample_arrow_geometry_unchanged=TRUE,
    labels_selected_unchanged=TRUE,source_matrices_labels_unchanged=TRUE,variance_significant_digits=2L,
    wrapped_loading_labels=TRUE,model_fits=0L,scientific_scores_computed=0L,root=root),file.path(root,'fixture.observed.json'),auto_unbox=TRUE,pretty=TRUE)
  cat('PASS renderer-only two-significant-digit axes/wrapped repelled loading labels; unchanged samples/features/arrows/selected labels; debug root=',root,'\n',sep='')
}
run()
