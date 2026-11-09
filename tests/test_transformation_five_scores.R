# Real, tiny composition screen; no RNA/preprocessing/embedding-model fits.
root<-normalizePath('.')
Sys.setenv(PROJECT_ROOT=root,RETICULATE_PYTHON=file.path(root,'.pixi/envs/default/bin/python'))
source('src/utils/imports_worker_transzeroimp.R')
source('src/utils/math_utils.R');source('src/utils/scoring_metrics.R')
source('src/utils/seurat_utils.R');source('src/5_run_benchmark_methods/benchmark_pipeline.R')
source('src/5_run_benchmark_methods/benchmark_hpc_utils.R')
reject<-function(call) stopifnot(inherits(tryCatch({force(call);NULL},error=function(e)e),'error'))
v<-rep(-5:4,2);w<-rep(c(-2,1,3,-1,0,2,-3,4,1,-2),2)
counts<-as.data.frame(cbind(reference=rep(10,20),A=30+v,B=25+w,C=35-v-w))
rownames(counts)<-paste0('s',1:20)
labels<-stats::setNames(factor(rep(c('x','y'),each=10)),rownames(counts))
result<-suppressWarnings(run_transformation_analysis(counts,labels))
validate_transformation_result(result)
reject(validate_transformation_result(result,expected_samples='foreign'))
policy<-transformation_analysis_policy();draws<-attr(result,'transformation_replicates')
stopifnot(nrow(result)==6L,nrow(draws)==120L,policy$modularity_k==3L,
 all(draws$alr_reference[draws$trans_method=='alr_mincvref']=='reference'),
 all(is.na(draws$alr_reference[!grepl('^alr_',draws$trans_method)])))
# Minimum CV is uniquely the constant-frequency reference, not a random draw.
imputed<-as.matrix(counts)+1;freq<-imputed/rowSums(imputed)*100
stopifnot(which.min(apply(freq,2,cv))==1L)
set.seed(policy$seed,kind=policy$rng_kind,normal.kind=policy$rng_normal_kind,sample.kind=policy$rng_sample_kind)
expected_random<-replicate(20,sample(colnames(counts),1))
stopifnot(identical(draws$alr_reference[draws$trans_method=='alr_randref'],expected_random))
# Exact feature math and all five canonical score meanings for one stored draw.
i<-which(draws$trans_method=='alr_mincvref')[1]
expected<-log(sweep(freq[,c('A','B','C')],1,freq[,'reference'],'/'))
stopifnot(max(abs(as.matrix(draws$feat_mat[[i]])-expected))<1e-12)
assign('.Random.seed',draws$rng_score_state[[i]],envir=.GlobalEnv)
score<-suppressWarnings(calc_sep_score(dist(draws$feat_mat[[i]]),labels,parallel=1L))
for(key in names(policy$score_columns)) stopifnot(abs(draws[[policy$score_columns[[key]]]][i]-score[[key]])<1e-10)
# Labels are evaluation-only, including their permutation RNG consumption.
permuted<-labels[c(1:5,11:15,6:10,16:20)];names(permuted)<-names(labels)
changed<-suppressWarnings(datrans(counts,permuted,0,0,reps=20,n_cores=1,rng_seed=policy$seed))
stopifnot(identical(draws$alr_reference,changed$alr_reference),
 identical(draws$feat_mat,changed$feat_mat),identical(draws$rng_feature_state,changed$rng_feature_state),
 any(abs(draws$Silhouette_score-changed$Silhouette_score)>1e-8))
# Validation neither fits models nor silently accepts old three-score summaries.
legacy<-result[,c('trans_method','ANOSIM_score','Modularity_score','Adjusted_Rand_Index')]
attributes(legacy)<-attributes(legacy)[c('names','row.names','class')]
reject(validate_transformation_result(legacy))
bad<-result;attr(bad,'transformation_policy')$modularity_k<-6L
reject(validate_transformation_result(bad))
bad<-result;bad$LISI_score[1]<-bad$LISI_score[1]+0.1
reject(validate_transformation_result(bad))
bad<-result;attr(bad,'transformation_replicates')<-draws[-1,]
reject(validate_transformation_result(bad))
namespace<-new.env(parent=baseenv());loaded<-loadedNamespaces()
sys.source(file.path(root,'src/5_run_benchmark_methods/benchmark_pipeline.R'),envir=namespace)
namespace$validate_transformation_result(result)
stopifnot(!any(c('Seurat','DESeq2')%in%setdiff(loadedNamespaces(),loaded)))
# Exercise the actual strict artifact CLI with checksums before reads.
tmp<-tempfile('ecoda-trans-five-');dir.create(tmp)
file<-file.path(tmp,'Toy_trans.rds')
write_bundle<-function(value) {
 saveRDS(value,file)
 writeLines(c(paste0('MD5=',unname(tools::md5sum(file))),paste0('SIZE=',file.info(file)$size),paste0('PATH=',file)),paste0(file,'.md5'))
}
validate_cli<-function() {
 output<-suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),
  c('--vanilla','src/5_run_benchmark_methods/validate_benchmark_rds_contract.R','--artifact',shQuote(file),'--method','trans'),stdout=TRUE,stderr=TRUE))
 status<-attr(output,'status');if(is.null(status)) status<-0L
 list(status=status,output=output)
}
write_bundle(result);check<-validate_cli();stopifnot(check$status==0L)
write_bundle(legacy);check<-validate_cli();stopifnot(check$status!=0L,any(grepl('version-2',check$output)))
# Exercise the real obs-only driver against a tiny H5AD whose RNA arrays are
# deliberately unreadable external links. Opening AnnData/layers would fail.
reticulate::py_run_string(r"(
import h5py, numpy as np

def metadata_fixture(path, counts, samples, types, labels):
    sample_values, type_values, label_values = [], [], []
    for i, sample in enumerate(samples):
        for j, cell_type in enumerate(types):
            n = int(counts[i,j])
            sample_values.extend([sample]*n)
            type_values.extend([cell_type]*n)
            label_values.extend([labels[i]]*n)
    with h5py.File(path,'w') as h:
        for key in ['X','layers/counts']:
            node=h.create_group(key); node.attrs['shape']=[len(sample_values),2]
            node.attrs['encoding-type']='csr_matrix'
            node['data']=h5py.ExternalLink('/missing-fixture-RNA.h5','/data')
        obs=h.create_group('obs');obs.attrs['encoding-type']='dataframe';obs.attrs['_index']='_index'
        for key,values in {'_index':['c'+str(i) for i in range(len(sample_values))],
                          'Sample':sample_values,'ct':type_values,'biology':label_values}.items():
            obs.create_dataset(key,data=np.asarray(values,dtype='S'))
        var=h.create_group('var');var.attrs['_index']='_index'
        var.create_dataset('_index',data=np.asarray(['g1','g2'],dtype='S'))
        var.create_dataset('hvg_rank',data=[1.,2.])
)")
h5ad<-file.path(tmp,'tiny_metadata_only.h5ad')
reticulate::py$metadata_fixture(h5ad,as.matrix(counts),as.list(rownames(counts)),as.list(colnames(counts)),as.list(as.character(labels)))
config_path<-file.path(tmp,'config.json');writeLines('{}',config_path)
output_dir<-file.path(tmp,'driver');dir.create(output_dir)
flags<-c('--config_path',config_path,'--ds_name','Toy','--view','benchmark_analysis',
 '--input_dir',tmp,'--output_dir',output_dir,'--log_file',file.path(tmp,'runtime.feather'))
commandArgs<-function(trailingOnly=FALSE) flags
read_datasets_json<-function(...) list(Toy=list(sample_col='Sample',label_col='biology',cell_type_high_res='ct'))
get_h5ad_path<-function(...) h5ad
# The real h5py helper remains in use; forbid the old unqualified AnnData import.
import<-function(...) stop('AnnData must not be opened for obs-only transformation input')
artifact_checksum_ok<-function(path,...) file.exists(path)
exec_time<-function(expr) {force(expr);as.difftime(0,units='secs')}
peak_rss_gb<-function() 0
save_rds_atomic<-function(object,file,...) saveRDS(object,file)
write_runtime_metadata<-function(...) invisible(NULL)
log_exec_row<-function(...) invisible(NULL)
replayed<-0L
replay_runtime_metadata<-function(...) {replayed<<-replayed+1L;invisible(NULL)}
fit_calls<-0L
screen<-function(ct,eval_labels) {
 fit_calls<<-fit_calls+1L
 stopifnot(isTRUE(all.equal(as.matrix(ct),as.matrix(counts))),identical(eval_labels,labels))
 result
}
run_ct_comps_analysis_worker('trans',screen,'_trans','trans_analysis')
persisted<-readRDS(file.path(output_dir,'Toy_trans.rds'))
validate_transformation_result(persisted)
stopifnot(fit_calls==1L,attr(persisted,'transformation_input')$cell_type_column=='ct')
run_ct_comps_analysis_worker('trans',screen,'_trans','trans_analysis')
stopifnot(fit_calls==1L,replayed==1L)
# Old score-only cache stays unchanged and cannot trigger silent replacement.
old_path<-file.path(output_dir,'Toy_trans.rds');saveRDS(legacy,old_path)
old_md5<-tools::md5sum(old_path)
reject(run_ct_comps_analysis_worker('trans',screen,'_trans','trans_analysis'))
stopifnot(identical(tools::md5sum(old_path),old_md5),fit_calls==1L)
# Zero-imputation still uses its existing producer/cache policy on the shared
# metadata path; it is not made to satisfy the transformation-v2 contract.
run_ct_comps_analysis_worker('zeroimp',function(ct,eval_labels) {
 stopifnot(isTRUE(all.equal(as.matrix(ct),as.matrix(counts))),identical(eval_labels,labels))
 list(existing_eight_scores=score)
},'_zeroimp','zeroimp_analysis')
stopifnot(length(readRDS(file.path(output_dir,'Toy_zeroimp.rds'))[[1]])==8L)
unlink(tmp,recursive=TRUE)
cat('PASS: five canonical transformation scores, corrected ALR semantics, seeded/label-blind features, retained draws and strict old-summary rejection\n')
