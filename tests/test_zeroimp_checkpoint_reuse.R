# Four validated compact bytepairs are copied; only seven cohorts are evaluated.
suppressPackageStartupMessages({library(dplyr);library(survival)})
for(p in c('src/utils/analysis_release.R','src/5_run_benchmark_methods/run_local_ecoda_derived.R','src/5_run_benchmark_methods/benchmark_pipeline.R','src/utils/math_utils.R','src/notebooks/run_zeroimp_count_screen.R'))source(p)
fail<-function(x)stopifnot(inherits(tryCatch(force(x),error=identity),'error'))
root<-tempfile('zeroimp-reuse-');dir.create(root);root<-normalizePath(root)
on.exit<-function(...)NULL
md5<-function(p)unname(tools::md5sum(p));sha<-function(p)digest::digest(file=p,algo='sha256')
ids<-paste0('s',1:12);counts<-matrix(as.double(1:48),12,dimnames=list(ids,paste0('CT',1:4)));counts[c(1,18,27,40)]<-0
input<-list(counts=counts,labels=setNames(factor(rep(c('A','B'),each=6)),ids))
score<-function(...)setNames(as.list(rep(.5,8)),c(zeroimp_policy()$metrics,'mod_knnsqrtn_score','mod_knn6_score','mod_knn9_score'))
template<-suppressWarnings(zeroimp_run_cohort(input,score_fun=score))
old_id<-'stage5_p6_zeroimp11_20261110_v1';old<-file.path(root,'_ecoda_runs',old_id);output<-file.path(old,'artifacts/zeroimp_count_screen/analysis');dir.create(output,recursive=TRUE)
rows<-facts<-list();results<-list()
for(ds in zeroimp_datasets())for(ext in c('rds','csv')){
 path<-file.path(output,paste0(ds,'_zeroimp_counts.',ext));rows[[length(rows)+1L]]<-list(path=path,exists=ds%in%zeroimp_datasets()[1:4])
}
for(ds in zeroimp_datasets()[1:4]){
 p<-template;p$schema<-zeroimp_policy()$schema;p$dataset<-ds;p$source<-list(dataset=ds);p$run_id<-old_id;p$config_md5<-strrep('b',32);p$runner_md5<-strrep('c',32);p$release_md5<-'e14f1b40ba776731c86f2e1c421e8b64'
 csv<-file.path(output,paste0(ds,'_zeroimp_counts.csv'));rds<-file.path(output,paste0(ds,'_zeroimp_counts.rds'))
 .ecoda_save_csv(p$results,csv);p$csv_md5<-md5(csv);.ecoda_save_rds(p,rds)
 results[[ds]]<-list(run_id=old_id,release_md5=p$release_md5,runner_md5=p$runner_md5,full_PASS_geometry_and_five_validated=TRUE,recipes43=lapply(p$recipes,function(x)x[c("setting","status")]))
 for(path in c(rds,csv)){
  i<-which(vapply(rows,function(x)identical(x$path,path),logical(1)));rows[[i]]$md5<-md5(path);rows[[i]]$stat<-list(st_size=file.info(path)$size)
  owner<-file.path(root,paste0(basename(path),'.owner'));writeLines(c(paste0('RUN_ID=',old_id),'STAGE=stage5',paste0('PATH=',path),'STATE=FAIL','PID=1'),owner)
  for(q in c(path,paste0(path,'.md5'),owner))facts[[length(facts)+1L]]<-list(path=q,sha256=sha(q))
 }
}
put<-function(x,name){path<-file.path(root,name);jsonlite::write_json(x,path,auto_unbox=TRUE,null='null',digits=NA);list(path=path,sha256=sha(path))}
# Geometry is validated again from the actual RDS; JSON evidence stores only its completion facts.
science<-put(list(schema='ecoda_zeroimp_checkpoint4_science_observation_v1',status='VALIDATED_CHECKPOINTS_PENDING_ROOT_REUSE_REVIEW',results4=results,checked_before_after=facts,unpublished_checkpoint_records_absent=TRUE,scientific_reuse_approved=FALSE,original_source_commit='3736052f4451ca4dda80684a4af840dac5f37f4b',original_gate='FAILED_PRESERVED',original_owners='FAIL_PRESERVED',model_fits=0,score_calculations=0,H5AD_reads=0),'science.json')
inventory<-put(list(schema='ecoda_zeroimp11_partial_byte_inventory_v1',rows=rows),'inventory.json')
review_value<-list(schema='ecoda_zeroimp_checkpoint4_root_reuse_review_v1',status='APPROVED_FIXED_FOUR_TYPED_CHECKPOINT_REUSE',review_complete=TRUE,typed_checkpoint_science_validated=TRUE,reuse_without_imputation_or_scoring=TRUE,generic_failed_input_bypass_authorized=FALSE,original_gate_approved=FALSE,original_gate_state='FAILED',original_stage5_owners='FAIL_PRESERVED',original_published_records='ABSENT_NOT_FABRICATED',original_run_id=old_id,original_source_commit='3736052f4451ca4dda80684a4af840dac5f37f4b',original_control_md5='e14f1b40ba776731c86f2e1c421e8b64',datasets=zeroimp_datasets()[1:4],recipes_each=43L,science_report_sha256=science$sha256,byte_inventory_sha256=inventory$sha256,new_model_fits=FALSE)
review<-put(review_value,'review.json');control<-list(zeroimp_checkpoint_reuse=list(schema='ecoda_zeroimp_checkpoint_reuse_v1',review=review,science=science,inventory=inventory))
reuse<-zeroimp_reuse_plan(control);zeroimp_check_reuse(reuse)
original<-tools::md5sum(c(vapply(facts,`[[`,character(1),'path'),science$path,review$path,inventory$path))
destination<-file.path(root,'new');dir.create(destination);calls<-0L
run_fun<-function(input){calls<<-calls+1L;template}
for(ds in zeroimp_datasets()){
 rds<-file.path(destination,paste0(ds,'_zeroimp_counts.rds'));csv<-file.path(destination,paste0(ds,'_zeroimp_counts.csv'))
 p<-zeroimp_produce_package(input,list(dataset=ds),ds,rds,csv,list(run_id='new_run',config_md5=strrep('b',32),runner_md5=strrep('d',32),release_md5=strrep('f',32)),reuse,run_fun)
 stopifnot(length(p$recipes)==43L)
 if(ds%in%reuse$datasets)stopifnot(identical(md5(rds),md5(file.path(output,basename(rds)))),identical(md5(csv),md5(file.path(output,basename(csv)))),p$run_id==old_id,p$runner_md5==strrep('c',32),p$release_md5=='e14f1b40ba776731c86f2e1c421e8b64')else stopifnot(p$run_id=='new_run',p$runner_md5==strrep('d',32))
 stopifnot(.ecoda_validate_output_sidecar(rds)$PATH==rds,.ecoda_validate_output_sidecar(csv)$PATH==csv)
}
stopifnot(calls==7L,identical(original,tools::md5sum(names(original))),length(list.files(destination))==44L)
fail(zeroimp_produce_package(input,list(dataset='Adams'),'Adams',file.path(destination,'Adams_zeroimp_counts.rds'),file.path(destination,'Adams_zeroimp_counts.csv'),list(run_id='new',config_md5=strrep('b',32)),reuse,run_fun))
owner<-facts[[3]]$path;lines<-readLines(owner);writeLines(sub('STATE=FAIL','STATE=ACTIVE',lines),owner);fail(zeroimp_check_reuse(reuse));writeLines(lines,owner)
record<-file.path(old,'manifests/artifacts',paste0(substr(digest::digest(file.path(output,'Adams_zeroimp_counts.rds'),algo='sha256',serialize=FALSE),1,32),'.record'));dir.create(dirname(record),recursive=TRUE);writeLines('unexpected',record);fail(zeroimp_check_reuse(reuse));unlink(record)
bad<-review_value;bad$typed_checkpoint_science_validated<-FALSE;control$zeroimp_checkpoint_reuse$review<-put(bad,'bad-review.json');fail(zeroimp_reuse_plan(control))
stopifnot(identical(original,tools::md5sum(names(original))));unlink(root,recursive=TRUE)
cat('PASS: 4 byte-identical checkpoints +7 evaluations,43statuses each, original identities/proofs immutable; invalid review/owner/record/target rejected\n')
