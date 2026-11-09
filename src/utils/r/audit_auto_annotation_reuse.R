# Finite read-only source audit; CLI creates only its new evidence JSON.
# CLI: SOURCE_ROOT SCRATCH_ROOT RESULT_ROOT INPUT_REFS INVENTORY CHECKSUM_TABLE OWNERS RUN_ID OUTPUT_JSON
AUTO_ANNOTATION_FIELDS <- c('layer1','layer2','layer3',paste0('layer_',1:6),'scATOMIC_pred')
auto_annotation_column <- function(key,columns=list()) {
 if(key %in% c('ECODA_authors_LR','ECODA_authors_HR')) {
  field<-if(key=='ECODA_authors_LR')'cell_type_low_res' else 'cell_type_high_res'
  value<-columns[[field]]
  if(length(value)==1L&&value%in%AUTO_ANNOTATION_FIELDS) return(value)
  return(NA_character_)
 }
 if(grepl('^ECODA_HiTME_.*layer[123]($|_)',key)) return(sub('^.*(layer[123])($|_.*)','\\1',key))
 if(identical(key,'ECODA_scATOMIC_HR')) return('scATOMIC_pred')
 if(grepl('^ECODA_scATOMIC_.*layer_?[1-6]($|_)',key)) return(paste0('layer_',sub('^.*layer_?([1-6])($|_.*)','\\1',key)))
 NA_character_
}

# Subtract 0.5 only after the convention is consistent with the actual stored
# CLR and Euclidean distance. This checks the old bundle; it fits no model.
auto_annotation_bundle_counts <- function(bundle) {
 counts <- as.matrix(bundle$counts); features <- as.matrix(bundle$feat_mat)
 distance <- as.matrix(bundle$dist_mat)
 stopifnot(is.numeric(counts),length(counts)>0,all(is.finite(counts)),
   all(counts>=0.5),any(abs(counts-0.5)<1e-8),
   all(abs(counts-0.5-round(counts-0.5))<1e-8),
   !is.null(rownames(counts)),!is.null(colnames(counts)),
   !anyDuplicated(rownames(counts)),!anyDuplicated(colnames(counts)),
   is.numeric(features),all(is.finite(features)),
   setequal(rownames(features),rownames(counts)),
   ncol(features)>0,all(colnames(features)%in%colnames(counts)),
   setequal(rownames(distance),rownames(counts)),setequal(colnames(distance),rownames(counts)))
 imputed <- counts[rownames(features),colnames(features),drop=FALSE]
 expected <- log(imputed)-rowMeans(log(imputed))
 stopifnot(max(abs(expected-features))<1e-8,
   max(abs(as.matrix(stats::dist(features))-distance[rownames(features),rownames(features),drop=FALSE]))<1e-8)
 counts-0.5
}

auto_annotation_compare <- function(bundle,current) {
 raw <- tryCatch(auto_annotation_bundle_counts(bundle),error=function(e)e)
 if(inherits(raw,'error')) return(list(status='BLOCKED',reason='stored counts do not substantiate counts_all +0.5 with stored CLR/distance'))
 if(is.null(current)) return(list(status='BLOCKED',reason='current annotation column is absent'))
 common_samples <- intersect(rownames(raw),rownames(current)); common_types<-intersect(colnames(raw),colnames(current))
 difference<-current[common_samples,common_types,drop=FALSE]-raw[common_samples,common_types,drop=FALSE]
 exact<-setequal(rownames(raw),rownames(current))&&setequal(colnames(raw),colnames(current))&&all(difference==0)
 list(status=if(exact)'COUNT_MATCH_PROVENANCE_PENDING' else 'INVALID_CURRENT_COUNTS',
   reason=if(exact)'counts match; producer/reuse provenance still required' else 'stored raw counts differ from current metadata',
   pseudocount_check='counts_all +0.5 consistent with stored CLR/distance',
   max_abs_difference=if(length(difference))max(abs(difference)) else NA_real_,
   affected_samples=sum(rowSums(abs(difference))>0),n_different_bins=sum(difference!=0),
   current_only_samples=setdiff(rownames(current),rownames(raw)),legacy_only_samples=setdiff(rownames(raw),rownames(current)),
   current_only_types=setdiff(colnames(current),colnames(raw)),legacy_only_types=setdiff(colnames(raw),colnames(current)))
}

# Inspect only this exact producer-owned metadata file. Unknown schemas remain
# evidence, never an inferred source identity or permission to reuse a fit.
auto_annotation_source_identity <- function(path,dataset,expected_md5,run_root,run_id) {
 result<-list(path=path,status='MISSING',current_expected_h5ad_md5=expected_md5)
 if(!file.exists(path)) return(result)
 if(any(nzchar(Sys.readlink(c(path,dirname(path),run_root))))) {
  result$status<-'UNSAFE_METADATA_PATH';return(result)
 }
 if(!startsWith(normalizePath(path),paste0(normalizePath(run_root),'/manifests/'))) {
  result$status<-'UNSAFE_METADATA_PATH';return(result)
 }
 result$md5<-unname(tools::md5sum(path));result$size<-file.info(path)$size
 read_fields<-function(p) {
  if(!file.exists(p)) return(NULL)
  if(nzchar(Sys.readlink(p))) return(NULL)
  lines<-readLines(p,warn=FALSE)
  if(!length(lines)||any(!grepl('=',lines,fixed=TRUE))) return(NULL)
  keys<-sub('=.*$','',lines)
  if(anyDuplicated(keys)) return(NULL)
  as.list(stats::setNames(sub('^[^=]*=','',lines),keys))
 }
 sidecar<-read_fields(paste0(path,'.md5'))
 key<-substr(digest::digest(normalizePath(path),algo='sha256',serialize=FALSE),1,32)
 record_path<-file.path(run_root,'manifests/artifacts',paste0(key,'.record'))
 record<-read_fields(record_path)
 result$checksum_metadata<-sidecar;result$producer_record_path<-record_path;result$producer_record_metadata<-record
 same_path<-function(x) is.character(x)&&length(x)==1L&&file.exists(x)&&normalizePath(x)==normalizePath(path)
 proof_ok<-function(x) !is.null(x)&&identical(x$MD5,result$md5)&&identical(x$SIZE,as.character(result$size))&&same_path(x$PATH)
 sidecar_ok<-proof_ok(sidecar)&&setequal(names(sidecar),c('MD5','SIZE','PATH'))
 record_ok<-proof_ok(record)&&setequal(names(record),c('PATH','SIZE','MD5','RUN_ID','PRODUCER','STATE'))&&identical(record$RUN_ID,run_id)&&identical(record$STATE,'PUBLISHED')&&is.character(record$PRODUCER)&&nzchar(record$PRODUCER)
 # Invalid existing receipts cannot be replaced by a different fallback claim.
 if((file.exists(paste0(path,'.md5'))&&!sidecar_ok)||(file.exists(record_path)&&!record_ok)||(!sidecar_ok&&!record_ok)) {
  result$status<-'UNVERIFIED_OR_CONFLICTING_METADATA_PROOF';return(result)
 }
 result$status<-'VERIFIED_METADATA_PROOF'
 if(result$size>1024^2) {result$status<-'BOUNDED_JSON_TOO_LARGE';return(result)}
 payload<-tryCatch(jsonlite::fromJSON(path,simplifyVector=FALSE),error=function(e)NULL)
 result$schema_keys<-head(names(payload),32L)
 if(is.null(payload)||!identical(payload$schema,1L)||!is.list(payload$entries)) {
  result$status<-'UNKNOWN_SCHEMA';return(result)
 }
 selected<-Filter(function(x)is.list(x)&&identical(x$dataset,dataset)&&identical(x$view,'benchmark_analysis'),payload$entries)
 if(length(selected)!=1L) {result$status<-'MISSING_OR_AMBIGUOUS_SELECTED_ROW';return(result)}
 entry<-selected[[1L]]
 required<-c('dataset','view','path','size','md5','sample_ids')
 if(!setequal(names(entry),required)||!is.character(entry$md5)||length(entry$md5)!=1L||!grepl('^[a-f0-9]{32}$',entry$md5)) {
  result$status<-'UNKNOWN_SELECTED_ROW_SCHEMA';return(result)
 }
 stopifnot(identical(unname(tools::md5sum(path)),result$md5))
 result$selected_row<-entry
 result$selected_h5ad_md5_matches_current<-identical(entry$md5,expected_md5)
 result
}

auto_annotation_audit <- function(source_root,scratch,results,references,inventory,legacy_manifest,owners,run_id) {
source_root <- normalizePath(source_root, mustWork=TRUE)
scratch <- normalizePath(scratch, mustWork=TRUE)
results <- normalizePath(results, mustWork=TRUE)
refs <- read.delim(references,header=FALSE,sep='\t',quote='',colClasses='character',stringsAsFactors=FALSE)
 frozen <- read.delim(inventory,header=FALSE,sep='\t',quote='',colClasses='character',stringsAsFactors=FALSE)
 stopifnot(ncol(frozen)==4,setequal(refs[[1]],c('Bassez','Zhang','Pelka')),nrow(refs)==3,all(refs[[2]]=='benchmark_analysis'),all(frozen[[1]]%in%refs[[1]]),!anyDuplicated(frozen[[2]]))
stopifnot(ncol(refs)==6)
stopifnot(requireNamespace('reticulate',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE),requireNamespace('digest',quietly=TRUE))
source(file.path(source_root,'src/utils/seurat_utils.R'))
registry<-jsonlite::fromJSON(file.path(source_root,'datasets.json'),simplifyVector=FALSE)
suppressPackageStartupMessages(library(dplyr))
# Selected obs and tiny uns values only. The shared loader never reads X/counts.
reticulate::py_run_string(paste0("import sys\nsys.path.insert(0, ",
 jsonlite::toJSON(file.path(source_root,'src/utils/py'),auto_unbox=TRUE),")"))
reticulate::py_run_string(r"(
import hashlib, h5py, numpy as np
from h5ad_obs_free import load_h5ad_obs_free

def metadata_only(path, sample_col):
    wanted = list(dict.fromkeys(['Sample', sample_col, 'layer1', 'layer2', 'layer3']
                  + ['layer_' + str(i) for i in range(1, 7)] + ['scATOMIC_pred']))
    with h5py.File(path, 'r') as h:
        present = [x for x in wanted if x in h['obs']]
        def small(node, depth=0):
            if isinstance(node, h5py.Group):
                return ({k: small(node[k], depth+1) for k in node.keys()}
                        if depth < 2 else {'keys': list(node.keys())})
            if node.size > 64:
                return {'shape': list(node.shape), 'dtype': str(node.dtype)}
            value = node[()]
            if isinstance(value, bytes): return value.decode()
            if isinstance(value, np.ndarray): return value.astype(str).tolist()
            if isinstance(value, np.generic): return value.item()
            return str(value)
        uns = {k: small(h['uns'][k]) for k in h.get('uns', {})
               if any(t in k.lower() for t in
                      ['annot', 'hitme', 'scatomic', 'cohort_selection', 'bassez'])}
    obs = load_h5ad_obs_free(path, present)
    sha = hashlib.sha256(obs.to_csv(index=True).encode()).hexdigest()
    return {'obs': obs, 'obs_sha256': sha, 'uns': uns}
)")
read_kv <- function(p) {
 if(!file.exists(p)) return(NULL)
 z<-readLines(p,warn=FALSE); ok<-grepl('=',z,fixed=TRUE)
 stats::setNames(sub('^[^=]*=','',z[ok]),sub('=.*$','',z[ok]))
}
producer_evidence <- function(p,dataset,current_expected_md5,artifact_md5) {
 sha<-digest::digest(normalizePath(p),algo='sha256',serialize=FALSE)
 key<-paste0(chartr('/:,\t |','______',basename(p)), '_',substr(sha,1,32))
 if(nchar(key,type='bytes')>100) key<-paste0('artifact_',substr(sha,1,32))
 owner<-file.path(owners,'artifact',key,'owner')
 # Metadata fields vary historically; report exact existing owner files only.
 directory<-dirname(owner)
 evidence<-list(owner_directory=directory,owner_files=if(dir.exists(directory)) list.files(directory) else character())
 if(file.exists(owner)) {
  fields<-read_kv(owner)
  evidence$owner_metadata<-as.list(fields)
  run<-unname(fields['RUN_ID'])
  if(length(run)==1&&!is.na(run)&&grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$',run)) {
   root<-file.path(scratch,'_ecoda_runs',run)
   evidence$producer_record<-as.list(read_kv(file.path(root,'manifests/artifacts',paste0(substr(sha,1,32),'.record'))))
   rec<-evidence$producer_record
   evidence$producer_record_matches_artifact<-length(rec)==6L&&setequal(names(rec),c('PATH','SIZE','MD5','RUN_ID','PRODUCER','STATE'))&&is.character(rec$PRODUCER)&&nzchar(rec$PRODUCER)&&identical(unname(fields['STATE']),'OK')&&identical(unname(fields['PATH']),normalizePath(p))&&identical(rec$PATH,normalizePath(p))&&identical(rec$MD5,artifact_md5)&&identical(rec$SIZE,as.character(file.info(p)$size))&&identical(rec$RUN_ID,run)&&identical(rec$STATE,'PUBLISHED')
   evidence$selected_source_identity<-auto_annotation_source_identity(file.path(root,'manifests/source_identity.json'),dataset,current_expected_md5,root,run)
   evidence$selected_source_identity$artifact_link_verified<-evidence$producer_record_matches_artifact
   evidence$run_metadata<-as.list(read_kv(file.path(root,'metadata')))
   evidence$source_manifest<-as.list(read_kv(file.path(root,'manifests/source.manifest')))
   evidence$runtime_identity<-as.list(read_kv(file.path(root,'manifests/runtime.identity')))
   evidence$available_manifests<-list.files(file.path(root,'manifests'))
   ip<-file.path(root,'manifests/input_artifacts.tsv')
   if(file.exists(ip)) {
    lines<-readLines(ip,warn=FALSE)
    evidence$selected_input_rows<-lines[startsWith(lines,paste0(dataset,'\t'))]
   }
  }
 }
 evidence
}
checksum_manifest <- legacy_manifest
stopifnot(file.exists(checksum_manifest))
manifest_lines <- readLines(checksum_manifest,warn=FALSE)
verify_rds <- function(p) {
 row<-frozen[frozen[[2]]==p,,drop=FALSE]
 stopifnot(nrow(row)==1)
 expected<-row[[3]]
 actual<-unname(tools::md5sum(p))
 if(!identical(actual,expected)||file.info(p)$size!=as.numeric(row[[4]])) stop('Legacy artifact checksum mismatch: ',p,'; expected=',expected,'; actual=',actual,'; expected_size=',row[[4]],'; actual_size=',file.info(p)$size)
 list(value=readRDS(p),md5=actual)
}
output<-list(schema=1L,purpose='metadata_only_auto_annotation_reuse_evidence',run_id=run_id,audit='auto_annotation_composition_reuse',checksum_manifest=checksum_manifest,checksum_manifest_md5=unname(tools::md5sum(checksum_manifest)),scope=c('Bassez','Zhang','Pelka'),cohorts=list())
for(ds in output$scope) {
 ref<-refs[refs[[1]]==ds & refs[[2]]=='benchmark_analysis',,drop=FALSE]
 stopifnot(nrow(ref)==1,file.exists(ref[[3]]),file.info(ref[[3]])$size==ref[[5]])
 proof<-read_kv(ref[[6]])
 stopifnot(!is.null(proof),identical(unname(proof['MD5']),ref[[4]]),identical(unname(proof['SIZE']),as.character(ref[[5]])),normalizePath(proof['PATH'])==normalizePath(ref[[3]]))
 sample_col<-registry[[ds]]$columns$sample
 stopifnot(is.character(sample_col),length(sample_col)==1L)
 meta<-reticulate::py$metadata_only(ref[[3]],sample_col)
 obs<-meta$obs
 stopifnot(all(c('Sample',sample_col)%in%names(obs)),!anyNA(obs$Sample),!anyNA(obs[[sample_col]]))
 alias<-unique(obs[,unique(c(sample_col,'Sample')),drop=FALSE])
 stopifnot(nrow(alias)==length(unique(obs$Sample)),nrow(alias)==length(unique(obs[[sample_col]])))
 auto_columns<-intersect(c('layer1','layer2','layer3',paste0('layer_',1:6),'scATOMIC_pred'),names(obs))
 current<-lapply(auto_columns,function(ct) as.matrix(get_ct_comp_df(obs,'Sample',ct)));names(current)<-auto_columns
 configured_ct_fields<-registry[[ds]]$columns[c('cell_type_low_res','cell_type_high_res')]
 overlap<-Filter(function(x)length(x)==1L&&x%in%AUTO_ANNOTATION_FIELDS,configured_ct_fields)
 dependency_keys<-c(ECODA_authors_LR='cell_type_low_res',ECODA_authors_HR='cell_type_high_res')
 dependency_keys<-names(dependency_keys)[dependency_keys%in%names(overlap)]
 rows<-list(); keys_by_artifact<-list(); files<-frozen[[2]][frozen[[1]]==ds]
 stopifnot(length(files)>0)
 for(p in files) {
  checked<-verify_rds(p); value<-checked$value
  if(!is.null(value$counts)) {
   key<-sub(paste0('^',ds,'_'),'',sub('\\.rds$','',basename(p)))
   stopifnot(grepl('HiTME|scATOMIC',key)||key%in%dependency_keys)
   bundles<-setNames(list(value),key)
  } else {
   keys<-names(value)[grepl('HiTME|scATOMIC',names(value))|names(value)%in%dependency_keys];bundles<-value[keys]
  }
  stopifnot(!anyDuplicated(names(bundles)))
  keys_by_artifact[[p]]<-as.list(names(bundles))
  for(key in names(bundles)) {
   b<-bundles[[key]];ct<-auto_annotation_column(key,registry[[ds]]$columns)
   row<-list(method=key,artifact=p,artifact_md5=checked$md5,annotation_column=ct,evidence_role=if(key%in%dependency_keys)'dependency_evidence' else 'automatic_composition_evidence',producer=producer_evidence(p,ds,ref[[4]],checked$md5),status='BLOCKED',reason='unknown annotation mapping or missing stored counts/current column')
   if(!is.na(ct)) row<-c(row[c('method','artifact','artifact_md5','annotation_column','evidence_role','producer')],auto_annotation_compare(b,current[[ct]]))
   rows[[length(rows)+1]]<-row
  }
  rm(value,checked,bundles);gc(verbose=FALSE)
 }
 output$cohorts[[ds]]<-list(input=ref[[3]],input_producer=producer_evidence(ref[[3]],ds,ref[[4]],ref[[4]]),expected_h5ad_md5=ref[[4]],expected_size=ref[[5]],proof=ref[[6]],checksum_scope='bound expected H5AD checksum/size metadata; no new whole-H5AD hash; current obs digest is fresh',obs_sha256=meta$obs_sha256,n_cells=nrow(obs),configured_sample_column=sample_col,canonical_sample_mapping=alias,configured_ct_fields=configured_ct_fields,configured_ct_autofield_overlap=overlap,dependency_methods=dependency_keys,annotation_columns_present=auto_columns,annotation_columns_absent=setdiff(c('layer1','layer2','layer3',paste0('layer_',1:6),'scATOMIC_pred'),auto_columns),annotation_uns=meta$uns,artifact_inventory=frozen[frozen[[1]]==ds,,drop=FALSE],artifact_keys=keys_by_artifact,invalid_methods=unique(vapply(Filter(function(x)x$status=='INVALID_CURRENT_COUNTS',rows),function(x)x$method,character(1))),blocked_methods=unique(vapply(Filter(function(x)x$status=='BLOCKED',rows),function(x)x$method,character(1))),rows=rows)
 rm(meta,obs,current,rows,keys_by_artifact,alias);gc(verbose=FALSE)
 reticulate::py_run_string('import gc; gc.collect()')
}
output
}

auto_annotation_cli <- function(a) {
 stopifnot(length(a)==9L)
 output<-a[[9L]];parent<-dirname(output);link<-Sys.readlink(output)
 if(!startsWith(output,'/')||!dir.exists(parent)||
    !identical(normalizePath(parent,mustWork=TRUE),parent)||
    !identical(file.path(parent,basename(output)),output)||
    file.exists(output)||(!is.na(link)&&nzchar(link)))
  stop('Evidence output must be an absent canonical absolute file in an existing directory')
 result<-auto_annotation_audit(a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
 temporary<-tempfile('.annotation-evidence-',tmpdir=parent)
 on.exit(unlink(temporary),add=TRUE)
 jsonlite::write_json(result,temporary,auto_unbox=TRUE,null='null',na='null',pretty=TRUE,digits=16)
 if(!file.link(temporary,output)) stop('Evidence output already exists or cannot be created exclusively')
 invisible(output)
}

if(sys.nframe()==0L) auto_annotation_cli(commandArgs(TRUE))
