# Real compact RDS copy/rekey path: no obs reads, counts, fitting or scoring.
source('src/5_run_benchmark_methods/run_local_ecoda_derived.R')
source('src/utils/seurat_utils.R')
reject<-function(x)stopifnot(inherits(tryCatch(force(x),error=identity),'error'))
# Configured canonical Sample is already authoritative, including hyphens.
obs<-data.frame(Sample=c('g001C-a','g002C-b'),label=c('A','B'))
stopifnot(identical(.ecoda_saved_sample_obs(obs,'Sample')$obs,obs))
obs<-data.frame(raw=c('PB-1','PB-2'),Sample=c('PB_1','PB_2'),label=c('A','B'))
stopifnot(identical(.ecoda_saved_sample_obs(obs,'raw')$obs$raw,obs$Sample))
x<-obs;x$Sample<-c('PB_1','PB_1');reject(.ecoda_saved_sample_obs(x,'raw'))
x<-obs;x$raw<-c('PB-1','PB_1');reject(.ecoda_saved_sample_obs(x,'raw'))
x<-obs;x$Sample<-NULL;reject(.ecoda_saved_sample_obs(x,'raw'))
root<-tempfile('saved11-reuse-');dir.create(root);root<-normalizePath(root)
project<-file.path(root,'source');dir.create(project);dir.create(file.path(project,'src'),recursive=TRUE)
dir.create(file.path(project,'src/5_run_benchmark_methods'),recursive=TRUE)
file.copy('datasets.json',file.path(project,'datasets.json'))
file.copy('src/5_run_benchmark_methods/run_local_ecoda_derived.R',file.path(project,'src/5_run_benchmark_methods/run_local_ecoda_derived.R'))
config<-jsonlite::fromJSON('datasets.json',simplifyVector=FALSE);selected<-ecoda_raw_config_union(config,'benchmark_union');stopifnot(length(selected)==11L)
runs<-file.path(root,'_ecoda_runs');dir.create(runs);old<-file.path(runs,'stage5_p6_saved22_20261110_v2');new<-file.path(runs,'recovery')
dir.create(file.path(old,'artifacts/saved_labels/res50'),recursive=TRUE);dir.create(file.path(new,'manifests'),recursive=TRUE)
refs<-list();expected<-list();inputs<-list()
for(ds in selected) {
 path<-file.path(root,ds,config[[ds]]$views$benchmark_analysis$output_file_name);dir.create(dirname(path));writeLines('not an H5AD; contents must not be opened by reuse',path)
 .ecoda_write_text_atomic(paste0(path,'.md5'),c(paste0('MD5=',.ecoda_md5(path)),paste0('SIZE=',file.info(path)$size),paste0('PATH=',path)))
 refs[[ds]]<-c(ds,'benchmark_analysis',path,.ecoda_md5(path),as.character(file.info(path)$size),paste0(path,'.md5'))
 ids<-if(ds=='Gongsharma_cmv_young_males')paste0('PB_',1:6) else if(ds=='Adams')paste0('g',1:6,'C-a') else paste0(ds,1:6)
 expected[[ds]]<-list(dataset=ds,sample_ids=as.list(ids),named_labels=as.list(setNames(rep(c('A','B'),3),ids)))
}
refpath<-file.path(new,'manifests/input_references.tsv');writeLines(vapply(refs,paste,collapse='\t',character(1)),refpath)
args<-list(input_references=refpath,input_references_sha256=digest::digest(file=refpath,algo='sha256'),scope='benchmark_union')
sources<-.ecoda_bound_sources(args,config,selected,'res50')
originals<-list();facts<-list()
for(ds in selected) {
 raw<-unlist(expected[[ds]]$sample_ids,use.names=FALSE);if(ds=='Gongsharma_cmv_young_males')raw<-gsub('_','-',raw,fixed=TRUE)
 counts<-matrix(as.double(1:18)+.5,6,dimnames=list(raw,paste0('cluster',1:3)));features<-log(counts)-rowMeans(log(counts))
 bundle<-list(feat_mat=as.data.frame(features),counts=as.data.frame(counts),labels=setNames(factor(rep(c('A','B'),3)),raw),
 dist_mat=stats::dist(features),scores=setNames(as.list(rep(.5,8)),c(.ecoda_screen_metrics(),'mod_knnsqrtn_score','mod_knn6_score','mod_knn9_score')),
 scoring_policy='calc_sep_score_five_v1_lisi_separation',source_identity=sources[[ds]])
 bundle$source_identity$count_fingerprint<-.ecoda_screen_digest(bundle$counts)
 path<-file.path(old,'artifacts/saved_labels/res50',paste0(ds,'_ECODA_seuratres_50.rds'));.ecoda_save_rds(bundle,path)
 originals[[ds]]<-bundle;facts[[path]]<-list(md5=.ecoda_md5(path))
}
control<-list(schema='ecoda_saved_res50_reuse_v1',status='APPROVED_SAVED_RES50_REUSE',review_complete=TRUE,reviewer='synthetic fixture',datasets=as.list(selected),original_run_id=basename(old))
diag<-list(original_run_id=basename(old),original_gate_state='FAILED',native_res50_manifest_state='COMPLETED',source_commit='457fd153c45a1c6b79b0ae281ea784f33019619d',source_file_facts=facts,rows=setNames(lapply(selected,function(x)list(dataset=x)),selected))
diagpath<-file.path(new,'manifests/res50_diagnostic.json');jsonlite::write_json(diag,diagpath,auto_unbox=TRUE)
control$diagnostic<-list(path=diagpath,sha256=digest::digest(file=diagpath,algo='sha256'))
jsonlite::write_json(control,file.path(new,'manifests/res50_reuse.json'),auto_unbox=TRUE)
jsonlite::write_json(list(rows=unname(expected)),file.path(new,'manifests/input_witness.json'),auto_unbox=TRUE)
before<-tools::md5sum(names(facts));Sys.setenv(ECODA_RUN_ROOT=new)
calc_sep_score<-get_ct_comp_df<-.ecoda_bound_counts<-function(...)stop('res50 reuse must not recount/read obs/score')
manifest<-.ecoda_reuse_saved_res50(project,new)
stopifnot(manifest$status=='COMPLETED',length(manifest$artifacts)==11L,!manifest$res50_reuse$counts_recomputed,!manifest$res50_reuse$scores_computed)
for(ds in selected) {
 p<-file.path(new,'artifacts/saved_labels/res50',paste0(ds,'_ECODA_seuratres_50.rds'));value<-readRDS(p);ids<-unlist(expected[[ds]]$sample_ids,use.names=FALSE)
 .ecoda_validate_bundle(value,ids,'ECODA_seuratres_50')
 if(ds=='Gongsharma_cmv_young_males') {
  stopifnot(!identical(.ecoda_md5(p),before[[names(facts)[match(ds,selected)]]]),
   identical(unname(as.matrix(value$feat_mat)),unname(as.matrix(originals[[ds]]$feat_mat))),identical(value$scores,originals[[ds]]$scores),
   identical(value$source_identity$identifier_projection$original_count_fingerprint,originals[[ds]]$source_identity$count_fingerprint),
   identical(value$source_identity$count_fingerprint,.ecoda_screen_digest(value$counts)),
   !identical(value$source_identity$count_fingerprint,originals[[ds]]$source_identity$count_fingerprint))
 } else stopifnot(identical(.ecoda_md5(p),unname(before[[names(facts)[match(ds,selected)]]])))
}
stopifnot(identical(before,tools::md5sum(names(facts))))
reject(.ecoda_reuse_saved_res50(project,new))
x<-originals[['Gongsharma_cmv_young_males']];reject(.ecoda_rekey_saved_res50(x,'Adams',unlist(expected[['Gongsharma_cmv_young_males']]$sample_ids,use.names=FALSE)))
# A late corrupt source is rejected before any RDS read or new output directory.
new2<-file.path(runs,'bad');dir.create(file.path(new2,'manifests'),recursive=TRUE);file.copy(list.files(file.path(new,'manifests'),full.names=TRUE),file.path(new2,'manifests'))
Sys.setenv(ECODA_RUN_ROOT=new2);last<-tail(names(facts),1);cat('corrupt',file=last,append=TRUE)
reader<-readRDS;readRDS<-function(...)stop('binary read before all11 guards')
err<-tryCatch(.ecoda_reuse_saved_res50(project,new2),error=identity);readRDS<-reader
stopifnot(inherits(err,'error'),grepl('checksum|bytes differ',conditionMessage(err)),!dir.exists(file.path(new2,'artifacts/saved_labels/res50')))
unlink(root,recursive=TRUE)
cat('PASS eleven real saved RDS reuse:10bytecopies,Gongidentifier-only,currentlabels,Adams-hyphens,no res50count/score/H5,guard-all-before-read\n')
