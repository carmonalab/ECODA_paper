# Focused metadata/RDS evidence: no cohort source or model fits.
root <- normalizePath('.')
source(file.path(root,'src/utils/r/audit_auto_annotation_reuse.R'))
reject <- function(call) stopifnot(inherits(tryCatch({force(call);NULL},error=function(e)e),'error'))
raw <- matrix(c(1,2,1,0),2,dimnames=list(c('s1','s2'),c('A','B')))
bundle <- function(counts) {
 imputed<-counts+0.5;feat<-log(imputed)-rowMeans(log(imputed))
 list(counts=as.data.frame(imputed),feat_mat=as.data.frame(feat),dist_mat=dist(feat),labels=c(s1='x',s2='y'))
}
b <- bundle(raw)
stopifnot(identical(auto_annotation_bundle_counts(b),raw),
 auto_annotation_compare(b,raw)$status=='COUNT_MATCH_PROVENANCE_PENDING',
 auto_annotation_compare(bundle(raw+matrix(c(0,1,0,0),2)),raw)$status=='INVALID_CURRENT_COUNTS',
 auto_annotation_compare(b,NULL)$status=='BLOCKED',
 auto_annotation_compare(list(),raw)$status=='BLOCKED',
 auto_annotation_column('ECODA_HiTME_HR_layer1')=='layer1',
 auto_annotation_column('ECODA_HiTME_HR_layer3_top_varexp0.5')=='layer3',
 auto_annotation_column('ECODA_scATOMIC_HR_layer_6')=='layer_6',
 auto_annotation_column('ECODA_scATOMIC_HR')=='scATOMIC_pred',
 is.na(auto_annotation_column('ECODA_HiTME_unknown')))
bad<-b;bad$dist_mat<-as.matrix(b$dist_mat)+1
stopifnot(auto_annotation_compare(bad,raw)$status=='BLOCKED')
bad<-b;bad$counts<-as.matrix(b$counts)+0.5
stopifnot(auto_annotation_compare(bad,raw)$status=='BLOCKED')
# Exercise the real obs loader, manifest-bound RDS reads, and complete cohort loop.
tmp<-tempfile('ecoda-auto-reuse-');dir.create(tmp);tmp<-normalizePath(tmp)
Sys.setenv(RETICULATE_PYTHON=file.path(root,'.pixi/envs/default/bin/python'))
suppressPackageStartupMessages(library(reticulate))
py_run_string("import anndata as ad,numpy as np,pandas as pd\ndef fixture(path):\n a=ad.AnnData(X=np.ones((4,2)),obs=pd.DataFrame({'Sample':['s1','s1','s2','s2'],'sample':['s1','s1','s2','s2'],'layer1':['A','B','A','A'],'layer2':['A','B','A','A'],'layer3':['A','B','A','A'],'scATOMIC_pred':['A','B','A','A'],'layer_6':['A','B','A','A']},index=['c1','c2','c3','c4']),var=pd.DataFrame({'hvg_rank':[1.,2.]},index=['g1','g2']))\n a.layers['counts']=a.X.copy();a.uns['ecoda_annotation_version']='fixture';a.write_h5ad(path)\n")
refs<-inventory<-checksum<-character()
results<-file.path(tmp,'benchmark/results');dir.create(results,recursive=TRUE)
for(ds in c('Bassez','Zhang','Pelka')) {
 h<-file.path(tmp,paste0(ds,'.h5ad'));py$fixture(h)
 md5<-unname(tools::md5sum(h));size<-file.info(h)$size;proof<-paste0(h,'.md5')
 writeLines(c(paste0('MD5=',md5),paste0('SIZE=',size),paste0('PATH=',h)),proof)
 refs<-c(refs,paste(ds,'benchmark_analysis',h,md5,size,proof,sep='\t'))
 values<-switch(ds,Bassez=list(ECODA_HiTME_HR_layer2=b,ECODA_authors_HR=b,ECODA_authors_LR=b),
  Zhang=list(ECODA_HiTME_HR_layer1=b,ECODA_scATOMIC_HR=bundle(raw+matrix(c(0,1,0,0),2)),ECODA_authors_LR=b,ECODA_authors_HR=bundle(raw+matrix(c(0,1,0,0),2))),
  Pelka=list(ECODA_HiTME_HR_layer3=bad,ECODA_HiTME_unknown=b,ECODA_authors_HR=b,ECODA_authors_LR=b))
 p<-file.path(results,paste0(ds,'_composition.rds'));saveRDS(values,p)
 inventory<-c(inventory,paste(ds,p,unname(tools::md5sum(p)),file.info(p)$size,sep='\t'))
 checksum<-c(checksum,paste(unname(tools::md5sum(p)),paste0('results/',basename(p))))
}
refpath<-file.path(tmp,'refs.tsv');ipath<-file.path(tmp,'inventory.tsv');cpath<-file.path(tmp,'checksums.md5')
writeLines(refs,refpath);writeLines(inventory,ipath);writeLines(checksum,cpath)
# Create the historical producer receipt through the actual shared API, not a
# reconstructed filename. This locks the report lookup to the real namespace.
producer_path<-file.path(results,'Bassez_composition.rds')
Sys.setenv(PROJECT_ROOT=root,HPC_SCRATCH_DIR=tmp,
 ECODA_RUNS_ROOT=file.path(tmp,'_ecoda_runs'),ECODA_OWNERS_ROOT=file.path(tmp,'_ecoda_owners'),
 EVIDENCE_TEST_RDS=producer_path)
code<-'set -euo pipefail
source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"
ecoda_init_run stage5 historical_producer >/dev/null
ecoda_write_checksum "$EVIDENCE_TEST_RDS"
ecoda_write_artifact_record "$EVIDENCE_TEST_RDS" composition historical_producer >/dev/null
ecoda_artifact_owner_acquire "$EVIDENCE_TEST_RDS" stage5 historical_producer 0 0 0 >/dev/null
ecoda_artifact_owner_set_state "$EVIDENCE_TEST_RDS" OK fixture >/dev/null'
stopifnot(system2('bash',c('-c',shQuote(code)))==0L)
identity<-file.path(tmp,'_ecoda_runs/historical_producer/manifests/source_identity.json')
expected_h5ad<-strsplit(refs[1],'\t')[[1]][4]
jsonlite::write_json(list(schema=1L,entries=list(list(dataset='Bassez',view='benchmark_analysis',path=file.path(tmp,'Bassez.h5ad'),size=file.info(file.path(tmp,'Bassez.h5ad'))$size,md5=expected_h5ad,sample_ids=c('s1','s2')))),identity,auto_unbox=TRUE)
Sys.setenv(EVIDENCE_TEST_IDENTITY=identity)
code<-'set -euo pipefail
source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"
ecoda_write_checksum "$EVIDENCE_TEST_IDENTITY"
ecoda_write_artifact_record "$EVIDENCE_TEST_IDENTITY" stage5_source_identity historical_producer >/dev/null'
stopifnot(system2('bash',c('-c',shQuote(code)))==0L)
verified<-auto_annotation_source_identity(identity,'Bassez',expected_h5ad,dirname(dirname(identity)),'historical_producer')
stopifnot(verified$status=='VERIFIED_METADATA_PROOF',verified$selected_h5ad_md5_matches_current,
 !auto_annotation_source_identity(identity,'Bassez',paste(rep('0',32),collapse=''),dirname(dirname(identity)),'historical_producer')$selected_h5ad_md5_matches_current)

unknown<-file.path(dirname(identity),'unknown_source_identity.json')
jsonlite::write_json(list(schema=99L,other='unknown'),unknown,auto_unbox=TRUE)
writeLines(c(paste0('MD5=',unname(tools::md5sum(unknown))),paste0('SIZE=',file.info(unknown)$size),paste0('PATH=',normalizePath(unknown))),paste0(unknown,'.md5'))
unknown_result<-auto_annotation_source_identity(unknown,'Bassez',expected_h5ad,dirname(dirname(identity)),'historical_producer')
stopifnot(unknown_result$status=='UNKNOWN_SCHEMA','other'%in%unknown_result$schema_keys)
writeLines('invalid checksum metadata',paste0(unknown,'.md5'))
stopifnot(auto_annotation_source_identity(unknown,'Bassez',expected_h5ad,dirname(dirname(identity)),'historical_producer')$status=='UNVERIFIED_OR_CONFLICTING_METADATA_PROOF')
report<-auto_annotation_audit(root,tmp,results,refpath,ipath,cpath,file.path(tmp,'_ecoda_owners'),'fixture')
stopifnot(report$cohorts$Bassez$rows[[1]]$status=='COUNT_MATCH_PROVENANCE_PENDING',
 'ECODA_scATOMIC_HR'%in%report$cohorts$Zhang$invalid_methods,
 setequal(report$cohorts$Pelka$blocked_methods,c('ECODA_HiTME_HR_layer3','ECODA_HiTME_unknown')),
 report$cohorts$Bassez$n_cells==4L,
 length(report$cohorts$Bassez$rows)==1L,length(report$cohorts$Pelka$rows)==2L,
 setequal(report$cohorts$Zhang$dependency_methods,c('ECODA_authors_LR','ECODA_authors_HR')),
 length(report$cohorts$Zhang$rows)==4L,
 all(vapply(Filter(function(x)x$method%in%c('ECODA_authors_LR','ECODA_authors_HR'),report$cohorts$Zhang$rows),function(x)x$evidence_role=='dependency_evidence',logical(1))),
 report$cohorts$Bassez$rows[[1]]$producer$selected_source_identity$selected_h5ad_md5_matches_current,
 report$cohorts$Bassez$rows[[1]]$producer$selected_source_identity$artifact_link_verified,
 report$cohorts$Bassez$rows[[1]]$producer$producer_record$PRODUCER=='composition',
 report$cohorts$Bassez$rows[[1]]$producer$producer_record$RUN_ID=='historical_producer',
 report$cohorts$Bassez$rows[[1]]$producer$producer_record$MD5==unname(tools::md5sum(producer_path)),
 length(report$cohorts$Bassez$obs_sha256)==1)
 # Actual Rscript startup writes its R_HOME warning to stdout before R code.
 # The explicit create-only output remains pure JSON and cannot be replaced.
 cli_output<-file.path(tmp,'cli-evidence.json')
 cli_args<-c('--vanilla',file.path(root,'src/utils/r/audit_auto_annotation_reuse.R'),
   root,tmp,results,refpath,ipath,cpath,file.path(tmp,'_ecoda_owners'),'fixture',cli_output)
 inventory_paths<-vapply(strsplit(inventory,'\t'),`[[`,character(1),2)
 before<-tools::md5sum(c(refpath,ipath,cpath,inventory_paths))
 invoke_cli<-function() suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),
   vapply(cli_args,shQuote,character(1)),stdout=TRUE,stderr=TRUE,
   env=c('R_HOME=/ecoda-deliberately-wrong-r-home')))
 logs<-invoke_cli()
 if(!is.null(attr(logs,'status'))) stop(paste(logs,collapse='\n'))
 stopifnot(is.null(attr(logs,'status')),
   any(grepl('ignoring environment value of R_HOME',logs,fixed=TRUE)),
   jsonlite::fromJSON(cli_output,simplifyVector=FALSE)$cohorts$Bassez$n_cells==4L,
   identical(tools::md5sum(names(before)),before))
 cli_md5<-unname(tools::md5sum(cli_output));logs<-invoke_cli()
 stopifnot(!is.null(attr(logs,'status')),unname(tools::md5sum(cli_output))==cli_md5)
 reject(auto_annotation_cli(c(rep('not-read',8),'relative.json')))
 p<-file.path(results,'Bassez_composition.rds');Sys.chmod(p,'0644');saveRDS(list(unexpected='changed'),p)
failure<-tryCatch(auto_annotation_audit(root,tmp,results,refpath,ipath,cpath,file.path(tmp,'_ecoda_owners'),'fixture'),error=function(e)e)
stopifnot(inherits(failure,'error'),grepl(p,conditionMessage(failure),fixed=TRUE),
 grepl('expected=',conditionMessage(failure),fixed=TRUE),
 grepl(paste0('actual=',unname(tools::md5sum(p))),conditionMessage(failure),fixed=TRUE))
unlink(tmp,recursive=TRUE)
cat('PASS: exact auto row inventory/counts, stored CLR/distance/pseudocount boundary, blocked ambiguity, obs-only loader, and checksum-before-RDS\n')
