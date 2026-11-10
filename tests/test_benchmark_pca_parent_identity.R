# Metadata-only PCA/A join: exact checked parent bytes, not purpose-string aliases.
source('src/utils/analysis_release.R');source('src/utils/evaluation_scores.R');source('src/notebooks/benchmark_analysis.R')
run<-function(){
 root<-tempfile('pca-parent-identity-');dir.create(root);root<-normalizePath(root);on.exit(unlink(root,recursive=TRUE))
 md5<-function(p)unname(tools::md5sum(p));fail<-function(expr)stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
 datasets<-c('Stephenson','Bassez','Adams','Gongsharma_cmv_young_males','Kfoury','Kim','Lee','Pelka','Smillie','Wu','Zhang')
 ledger<-scope<-origins<-list();i<-0L;ids<-paste0('s',1:4);y<-setNames(c('a','a','b','b'),ids)
 x<-matrix(1:12,4,dimnames=list(ids,letters[1:3]));bundle<-list(feat_mat=x,counts=x+.5,dist_mat=dist(x),labels=y)
 metadata<-list(labels=y,n_samples=4L,n_cells=sum(x),cells_per_sample=setNames(rowSums(x),ids),n_cell_types_high_res=3L)
 evidence<-file.path(root,'review.txt');writeLines('Synthetic approved parent identity',evidence)
 for(ds in datasets)for(method in c('ECODA_authors_HR','Pseudobulk_hvg2000','__dataset_metadata__',if(ds=='Bassez')'Pseudobulk_hvg3000')){
  i<-i+1L;p<-file.path(root,paste0(i,'.rds'));saveRDS(if(method=='__dataset_metadata__')metadata else bundle,p)
  r<-data.frame(dataset=ds,view='benchmark_analysis',method_key=method,kind='standalone_result',artifact_path=p,bundle_key='',run_root=root,md5=md5(p),receipt_path=file.path(root,paste0(i,'.receipt.rds')),receipt_md5='')
  policy<-if(ds=='Adams')'legacy_row_current_input_compatibility_v1'else if(ds=='Bassez')'reviewed_published_current_input_v1'else'new_current_input_required_v1'
  receipt<-c(as.list(r[c('dataset','view','method_key','kind','artifact_path','bundle_key','run_root','md5')]),list(schema='ecoda_analysis_release_receipt_v1',status='reviewed',origin=if(ds=='Adams'||ds=='Bassez'&&grepl('^Pseudobulk_',method))'legacy_reuse'else'new',input_md5=strrep('a',32),source_policy=policy,evidence='reviewed',evidence_path=evidence,evidence_md5=md5(evidence)))
  saveRDS(receipt,r$receipt_path);r$receipt_md5<-md5(r$receipt_path);ledger[[i]]<-r;origins[[i]]<-receipt$origin
  s<-r[c('dataset','view','method_key')];s$status<-'approved';s$require_new<-receipt$origin=='new';s$current_input_md5<-receipt$input_md5;s$source_policy<-policy;scope[[i]]<-s
 }
 control<-list(ledger=do.call(rbind,ledger),expected_scope=do.call(rbind,scope));stopifnot(nrow(control$ledger)==34L)
 cp<-file.path(root,'raw-control.rds');saveRDS(control,cp)
 source<-file.path(root,'source.manifest');runtime<-file.path(root,'runtime.identity');writeLines('source identity',source);writeLines('runtime identity',runtime)
 projected<-project_analysis_release(control,file.path(root,'projection'),c(source=source,runtime=runtime,control=cp))
 pc<-file.path(root,'projected-control.rds');saveRDS(projected,pc)
 # Fail immediately if metadata reader tries deserializing either scientific generation.
 native_read<-readRDS;blocked<-c(control$ledger$artifact_path,projected$ledger$artifact_path)
 assign('readRDS',function(file,...){if(as.character(file)%in%blocked)stop('SCIENTIFIC BODY READ FORBIDDEN');native_read(file,...)},envir=.GlobalEnv)
 on.exit(rm('readRDS',envir=.GlobalEnv),add=TRUE)
 current<-read_benchmark_analysis_control(pc,md5(pc));stopifnot(identical(current$expected_scope$source_artifact_md5,control$ledger$md5),
  !identical(current$expected_scope$source_artifact_md5,projected$ledger$md5))
 stopifnot(identical(read_benchmark_analysis_control(cp,md5(cp)),control)) # raw/same-purpose backcompat
 render<-file.path(root,'render');dir.create(render);names<-benchmark_pca_render_files(datasets)
 for(name in names)writeLines('synthetic finished output',file.path(render,name))
 schemas<-benchmark_pca_csv_columns();for(name in names(schemas)){
  table<-as.data.frame(setNames(rep(list(1),length(schemas[[name]])),schemas[[name]]));write.csv(table,file.path(render,paste0('Bassez_',name,'.csv')),row.names=FALSE)
 }
 sources<-control$ledger[c('dataset','method_key','artifact_path')];sources$artifact_md5<-control$ledger$md5;sources$origin<-unlist(origins);sources$input_md5<-strrep('a',32)
 sources$source_policy<-ifelse(sources$method_key%in%c('ECODA_authors_HR','__dataset_metadata__'),'CURRENT_REVIEWED_AUTHOR_HR_INTEGER_COUNTS_PLUS_0.5_AND_DATASET_METADATA_V1','CURRENT_REVIEWED_SAVED_PSEUDOBULK_FEATURES_AND_NAMED_GEOMETRY_V1')
 report<-list(schema='ecoda_saved_pca_render_v1',datasets=datasets,operation='RENDER_SAVED_FEATURES_WITHOUT_MODEL_REFIT',control_md5=strrep('b',32),sources=sources,
  files=data.frame(basename=names,md5=vapply(file.path(render,names),md5,character(1))),csv_bounds=setNames(as.list(rep(1L,6)),names(schemas)))
 manifest<-file.path(render,'render_manifest.rds');saveRDS(report,manifest)
 read<-function(s=current$expected_scope)read_benchmark_pca_render(manifest,md5(manifest),datasets,s)
 stopifnot(nrow(read()$sources)==34L);fail(read(control$expected_scope))
 bad<-current$expected_scope;bad$source_artifact_md5[1]<-strrep('0',32);fail(read(bad))
 bad<-current$expected_scope;bad$current_input_md5[1]<-strrep('0',32);fail(read(bad))
 original<-report;legacy_index<-which(!control$expected_scope$require_new)[1]
 report$sources$origin[legacy_index]<-'new';saveRDS(report,manifest);fail(read())
 report<-original;report$sources$source_policy[1]<-'unrecognized-purpose';saveRDS(report,manifest);fail(read())
 report<-original;report$sources$origin[1]<-'legacy_reuse';saveRDS(report,manifest);fail(read())
 report<-original;report$sources$source_policy<-control$expected_scope$source_policy;saveRDS(report,manifest);stopifnot(nrow(read(control$expected_scope)$sources)==34L)
 # A resealed altered parent hash cannot bypass the frozen source control proof.
 row<-projected$ledger[1,,drop=FALSE];receipt<-native_read(row$receipt_path);proof<-native_read(receipt$evidence_path)
 proof$parent$md5<-strrep('0',32);saveRDS(proof,receipt$evidence_path);receipt$evidence_md5<-md5(receipt$evidence_path);saveRDS(receipt,row$receipt_path)
 projected$ledger$receipt_md5[1]<-md5(row$receipt_path);saveRDS(projected,pc);fail(read_benchmark_analysis_control(pc,md5(pc)))
 cat('PASS:34 checked parent MD5s; known cross-purpose policies only; no core reads; input/origin/hash/proof negatives; same-purpose backcompat\n')
}
run()
