# Real producer assembly and actual artifact CLI; normalization/scoring are
# deliberately stubbed so this is a schema/blacklist check, never a model fit.
root<-normalizePath('.')
Sys.setenv(PROJECT_ROOT=root,ECODA_AUX_ROOT=file.path(root,'aux'))
env<-new.env(parent=globalenv())
for(file in c('src/utils/pseudobulk.R','src/utils/artifact_records.R',
 'src/5_run_benchmark_methods/benchmark_pipeline.R','src/5_run_benchmark_methods/benchmark_hpc_utils.R'))
 sys.source(file.path(root,file),envir=env)
resolved<-env$pseudobulk_blacklist()
blocked<-unique(unlist(lapply(resolved$groups,function(x)head(unname(x),2)),use.names=FALSE))
genes<-c(blocked,'XIST','SRY','retained_A','retained_B','retained_C')
counts<-matrix(seq_len(length(genes)*5),nrow=length(genes),dimnames=list(genes,paste0('s',1:5)))
metadata<-data.frame(Sample=colnames(counts),row.names=colnames(counts))
env$load_h5ad_sample_aggregate<-function(...,metadata_columns) {
 stopifnot(identical(metadata_columns,'Sample'))
 list(counts=counts,metadata=metadata,sample_ids=metadata$Sample)
}
env$fit_pseudobulk_deseq2<-function(counts,...) list(norm_matrix=counts,variance_order=rownames(counts))
env$exec_time<-function(expr) {force(expr);as.difftime(1,units='secs')}
env$peak_rss_gb<-function() 0.25
labels<-stats::setNames(factor(c('a','a','a','b','b')),colnames(counts))
env$process_pseudobulk_fig<-function(pb,labels) list(
 scores=list(sil_score=0.2,mod_knn3_score=0.3,cluster_score=0.4,anosim_score=0.5,lisi_score=0.6),
 feat_mat=pb,dist_mat=dist(pb),labels=labels)
tmp<-tempfile('ecoda-bl-cli-');dir.create(tmp)
input<-file.path(tmp,'stub.h5ad');writeLines('no H5AD read in this producer assembly fixture',input)
bl<-env$prepare_pseudobulks_hpc(input,hvg_rank_genes=rownames(counts),variants='hvg2000_bl',cache_stem='Toy',run_id='bl-cli')$hvg2000_bl
ordinary<-env$prepare_pseudobulks_hpc(input,hvg_rank_genes=rownames(counts),variants='hvg2000',cache_stem='Toy',run_id='bl-cli')$hvg2000
cache_dir<-file.path(tmp,'pseudobulks');result_dir<-file.path(tmp,'results')
dir.create(cache_dir);dir.create(result_dir)
cache<-file.path(cache_dir,'Toy_pseudobulk_hvg2000_bl.rds')
ordinary_file<-file.path(cache_dir,'Toy_pseudobulk_hvg2000.rds')
write_checked<-function(value,path) {
 saveRDS(value,path)
 writeLines(c(paste0('MD5=',unname(tools::md5sum(path))),paste0('SIZE=',file.info(path)$size),paste0('PATH=',path)),paste0(path,'.md5'))
}
write_checked(bl,cache);write_checked(ordinary,ordinary_file)
results<-env$run_pseudobulk_hpc(labels,list(hvg2000_bl=bl),results_dir=result_dir,ds='Toy',variants='hvg2000_bl')
result_file<-file.path(result_dir,'Toy_Pseudobulk_hvg2000_bl.rds')
original_result<-results$Pseudobulk_hvg2000_bl
cli<-function(path,method,valid=TRUE,pattern=NULL) {
 output<-suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',
  'src/5_run_benchmark_methods/validate_benchmark_rds_contract.R','--artifact',shQuote(path),'--method',method),stdout=TRUE,stderr=TRUE))
 status<-attr(output,'status');if(is.null(status)) status<-0L
 if((status==0L)!=valid) stop(paste(c(paste('Unexpected CLI status',status,'for',basename(path)),output),collapse='\n'))
 if(!is.null(pattern)) stopifnot(any(grepl(pattern,output)))
 invisible(output)
}
cli(cache,'prepare_pseudobulk')
cli(result_file,'pseudobulk')
cli(ordinary_file,'prepare_pseudobulk')
aggregate_file<-file.path(result_dir,'Toy_pseudobulk.rds')
write_checked(results,aggregate_file);cli(aggregate_file,'pseudobulk')
# A fresh checksum never makes a stale policy or blocked feature legitimate.
for(kind in c('stale','missing','malformed_audit','blocked','shared_cost','aggregate_cost','extra')) {
 value<-bl
 if(kind=='stale') value$blacklist_policy$source_md5<-paste(rep('0',32),collapse='')
 if(kind=='missing') value$blacklist_policy<-NULL
 if(kind=='malformed_audit') value$blacklist_audit$matched_rows_per_group[[1]]<- -1L
 if(kind=='blocked') colnames(value$pb)[1]<-blocked[[1]]
 if(kind=='shared_cost') {value$aggregate_time_secs<-1;value$shared_time_secs<-1}
 if(kind=='aggregate_cost') value$blacklist_aggregate_time_secs<-value$variant_time_secs+1
 if(kind=='extra') value$unrequested_metadata<-TRUE
 write_checked(value,cache);cli(cache,'prepare_pseudobulk',FALSE)
}
for(kind in c('stale','missing','malformed_audit','blocked','malformed_source')) {
 value<-original_result
 if(kind=='stale') value$blacklist_policy$version<-'historical_noop'
 if(kind=='missing') value$blacklist_policy<-NULL
 if(kind=='malformed_audit') value$blacklist_audit$retained_genes<-0L
 if(kind=='blocked') colnames(value$feat_mat)[1]<-blocked[[1]]
 if(kind=='malformed_source') value$blacklist_source_timing_id<-'run::view:pass'
 write_checked(value,result_file);cli(result_file,'pseudobulk',FALSE)
 write_checked(list(Pseudobulk_hvg2000_bl=value),aggregate_file);cli(aggregate_file,'pseudobulk',FALSE)
}
# Ordinary schema2 stays exact; BL metadata is not a general extra-field escape.
value<-ordinary;value$blacklist_policy<-resolved$policy
write_checked(value,ordinary_file);cli(ordinary_file,'prepare_pseudobulk',FALSE)
write_checked(bl,cache);write_checked(original_result,result_file);write_checked(ordinary,ordinary_file)
cli(cache,'prepare_pseudobulk');cli(result_file,'pseudobulk');cli(ordinary_file,'prepare_pseudobulk')
write_checked(results,aggregate_file);cli(aggregate_file,'pseudobulk')
unlink(tmp,recursive=TRUE)
cat('PASS: actual CLI accepts current BL cache/result and ordinary schema2, rejects stale/malformed/blocked BL and ordinary extras\n')
