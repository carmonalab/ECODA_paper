# Tiny synthetic metadata only; never opens production RNA/H5AD artifacts.
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = base::identity), "error"))
root <- tempfile("analysis-projection-"); dir.create(root); root <- normalizePath(root)
md5 <- function(path) unname(tools::md5sum(path))
ids <- paste0("s", 1:8)
labels <- setNames(factor(rep(c("a", "b"), each = 4)), ids)
ct <- matrix(seq_len(24), 8, dimnames = list(ids, c("T", "B", "myeloid")))
rna <- matrix(seq_len(160), 8, dimnames = list(ids, paste0("gene", 1:20)))
eight <- c(evaluation_metrics(), "mod_knnsqrtn_score", "mod_knn6_score", "mod_knn9_score")
scores <- setNames(as.list(seq_along(eight) / 10), eight)
bundle <- function(features) list(feat_mat = features, dist_mat = dist(features), labels = labels,
  scores = scores, counts = features + .5, RNA = list(matrix = rna), hidden_cache = rna,
  exec_time = 3, mem_GB = .2, timing_schema = 2L, timing_id = "original:shared",
  shared_time_secs = 10, shared_mem_GB = 1, variant_time_secs = 3,
  aggregate_time_secs = 4, shared_fit_time_secs = 6)
author <- bundle(ct); null <- author; null$labels <- labels[c(1,5,2,6,3,7,4,8)]; names(null$labels) <- ids
pb <- bundle(rna)
scitd <- pb; scitd$feat_mat <- rna[-8,,drop=FALSE]; scitd$dist_mat <- dist(scitd$feat_mat); scitd$labels <- labels[-8]
# Real legacy Glo IDs may carry incidental vector names through matrix dimnames.
# Validate an unname()d view only; retain source distance/labels/fingerprints.
named_pb<-pb;attr(named_pb$dist_mat,"Labels")<-setNames(ids,paste0("incidental",seq_along(ids)))
stopifnot(!identical(names(named_pb$labels),rownames(as.matrix(named_pb$dist_mat))))
item<-list(result=named_pb,source=list(method_key="GloScope_hvg2000_pcadims30",view="benchmark_analysis"))
original_fp<-c(distance=.analysis_fingerprint(named_pb$dist_mat),labels=.analysis_fingerprint(named_pb$labels))
thin_named<-.analysis_thin_result(item)
stopifnot(identical(thin_named$dist_mat,named_pb$dist_mat),identical(thin_named$labels,named_pb$labels),
 identical(original_fp,c(distance=.analysis_fingerprint(thin_named$dist_mat),labels=.analysis_fingerprint(thin_named$labels))))
named_ct<-author;dimnames(named_ct$feat_mat)[[1]]<-setNames(ids,paste0("featureAttr",seq_along(ids)))
stopifnot(!identical(rownames(named_ct$feat_mat),rownames(as.matrix(named_ct$dist_mat))))
# Kim's already-reviewed count rows may have their own exact name order.
named_ct$counts<-named_ct$counts[rev(ids),,drop=FALSE]
item_ct<-list(result=named_ct,source=list(method_key="ECODA_authors_HR",view="benchmark_analysis"))
thin_ct<-.analysis_thin_result(item_ct)
stopifnot(identical(thin_ct$feat_mat,named_ct$feat_mat),identical(thin_ct$counts,named_ct$counts),
 identical(thin_ct$labels,named_ct$labels),identical(thin_ct$dist_mat,named_ct$dist_mat))
bad<-item;attr(bad$result$dist_mat,"Labels")<-rev(ids);fail(.analysis_thin_result(bad))
bad<-item;names(bad$result$labels)<-rev(ids);fail(.analysis_thin_result(bad))
bad<-item;attr(bad$result$dist_mat,"Labels")<-c("wrong_sample",ids[-1]);fail(.analysis_thin_result(bad))
bad<-item;bad$result$labels[1]<-NA;fail(.analysis_thin_result(bad))
bad<-item_ct;rownames(bad$result$feat_mat)<-rev(ids);fail(.analysis_thin_result(bad))
bad<-item_ct;rownames(bad$result$counts)[1]<-"wrong_sample";fail(.analysis_thin_result(bad))
metadata <- list(labels=labels,n_samples=8L,n_cells=80L,cells_per_sample=setNames(rep(10,8),ids),n_cell_types_high_res=3L,RNA=rna)
aggregate <- file.path(root,"source.rds")
keys <- c("ECODA_authors_HR","ECODA_authors_HR_NULL","Pseudobulk_hvg2000","scITD_hvg2000_factors5","__dataset_metadata__")
saveRDS(setNames(list(author,null,pb,scitd,metadata),keys),aggregate)
evidence <- file.path(root,"science.txt"); writeLines("Explicit synthetic scientific row review; original origins unchanged",evidence)
ledger <- data.frame(dataset="D",view="benchmark_analysis",method_key=keys,kind="rds_bundle",artifact_path=aggregate,
  bundle_key=keys,run_root=root,md5=md5(aggregate),receipt_path=file.path(root,paste0(1:5,".rds")),receipt_md5="")
scope <- ledger[c("dataset","view","method_key")];scope$require_new<-FALSE;scope$current_input_md5<-strrep("a",32)
scope$source_policy<-"synthetic-current-input";scope$status<-"approved"
for(i in seq_len(nrow(ledger))) {
  receipt <- c(as.list(ledger[i,c("dataset","view","method_key","kind","artifact_path","bundle_key","run_root","md5")]),
    list(schema="ecoda_analysis_release_receipt_v1",status="reviewed",origin="legacy_reuse",input_md5=scope$current_input_md5[i],
      source_policy=scope$source_policy[i],evidence="reviewed",evidence_path=evidence,evidence_md5=md5(evidence),
      scoring_policy=evaluation_scoring_policy(),allowed_dropped_sample_ids=if(i==4L) "s8" else character()))
  saveRDS(receipt,ledger$receipt_path[i]);ledger$receipt_md5[i]<-md5(ledger$receipt_path[i])
}
control <- list(ledger=ledger,expected_scope=scope,reference_method_key=c(D="ECODA_authors_HR"))
timing_path<-file.path(root,"timing.feather")
arrow::write_feather(data.frame(dataset="D",method="Pseudobulk_hvg2000",time_secs=3,mem_GB=.2),timing_path)
control$execution_times<-list(artifact_path=timing_path,md5=md5(timing_path))
control_path<-file.path(root,"control.rds");saveRDS(control,control_path)
source_path<-file.path(root,"source.manifest");runtime_path<-file.path(root,"runtime.identity")
writeLines("pinned source fixture",source_path);writeLines("pinned runtime fixture",runtime_path)
validation_identity<-c(source=source_path,runtime=runtime_path,control=control_path)
# The early guard must not deserialize a source aggregate, including hidden keys.
readRDS_base <- base::readRDS; reads <- character()
readRDS <- function(file,...) {reads<<-c(reads,as.character(file));readRDS_base(file,...)}
fail(load_analysis_release(ledger,scope,require_projection=TRUE))
stopifnot(!aggregate %in% reads)
rm(readRDS)
protected <- tools::md5sum(c(aggregate,ledger$receipt_path,evidence))
projected <- project_analysis_release(control,file.path(root,"projected"),validation_identity)
stopifnot(identical(protected,tools::md5sum(names(protected))))
items <- load_analysis_release(projected$ledger,projected$expected_scope,require_projection=TRUE)
stopifnot(length(items)==5L,identical(items[[1]]$result$feat_mat,ct),identical(items[[1]]$result$counts,author$counts),
  identical(items[[2]]$result$labels,null$labels),identical(items[[3]]$result$dist_mat,pb$dist_mat),
  identical(items[[3]]$result$labels,labels),identical(items[[3]]$result$scores,scores),
  identical(items[[3]]$result$n_features,20L),is.null(items[[3]]$result$feat_mat),is.null(items[[3]]$result$counts),
  identical(items[[4]]$receipt$allowed_dropped_sample_ids,"s8"),
  identical(items[[4]]$result$labels,scitd$labels),is.null(items[[5]]$result$RNA),
  all(vapply(items,function(x) identical(x$receipt$origin,"legacy_reuse"),logical(1))))
for(item in items[1:4]) {
  stopifnot(is.null(item$result$RNA),is.null(item$result$hidden_cache),identical(item$result$timing_id,"original:shared"))
}
# A metadata-only transfer may land at an unrelated canonical local root.
# The raw copied producer tree and scientific witnesses remain unchanged.
producer_root<-file.path(root,"projected");raw_copy<-file.path(root,"copied-thin-root");dir.create(raw_copy)
manifest<-read.delim(file.path(producer_root,"files.tsv"),header=FALSE,col.names=c("path","md5"))
for(i in seq_len(nrow(manifest))) {
  target<-file.path(raw_copy,substring(manifest$path[i],nchar(producer_root)+2L))
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  stopifnot(file.copy(manifest$path[i],target))
}
stopifnot(file.copy(file.path(producer_root,"files.tsv"),file.path(raw_copy,"files.tsv")))
raw_files<-c(file.path(raw_copy,substring(manifest$path,nchar(producer_root)+2L)),file.path(raw_copy,"files.tsv"))
raw_hashes<-tools::md5sum(raw_files)
stopifnot(file.rename(producer_root,file.path(root,"unavailable-HPC-projection")))
mirrored<-mirror_analysis_projection(raw_copy,file.path(root,"local-thin-root"),md5(file.path(raw_copy,"analysis_release.rds")))
local_items<-load_analysis_release(mirrored$ledger,mirrored$expected_scope,require_projection=TRUE)
stopifnot(identical(local_items[[2]]$result$labels,null$labels),is.null(local_items[[3]]$result$feat_mat),
 identical(local_items[[3]]$result$dist_mat,pb$dist_mat),identical(raw_hashes,tools::md5sum(raw_files)),
 identical(md5(mirrored$execution_times$artifact_path),md5(timing_path)),
 identical(md5(file.path(raw_copy,"parents/validation.control")),md5(file.path(root,"local-thin-root/parents/validation.control"))))
wrong<-mirrored$expected_scope;wrong$require_new[3]<-TRUE
fail(load_analysis_release(mirrored$ledger,wrong,require_projection=TRUE))
fail(mirror_analysis_projection(raw_copy,file.path(root,"local-thin-root"),md5(file.path(raw_copy,"analysis_release.rds"))))
stopifnot(file.rename(file.path(root,"unavailable-HPC-projection"),producer_root))
fail(project_analysis_release(control,file.path(root,"projected"),validation_identity))
blocked<-control;blocked$expected_scope$status[1]<-"blocked";blocked$expected_scope$blocked_reason<-"Unapproved legacy evidence"
fail(project_analysis_release(blocked,file.path(root,"blocked"),validation_identity));stopifnot(!dir.exists(file.path(root,"blocked")))
loaded_control<-control;loaded_control$loaded<-list(RNA=rna)
fail(project_analysis_release(loaded_control,file.path(root,"loaded-control"),validation_identity))
# A new file/run is a projection, never a new scientific result. Rehashing the
# access receipt cannot promote its copied parent receipt from legacy to new.
bad<-projected;r<-readRDS(bad$ledger$receipt_path[3]);r$origin<-"new";saveRDS(r,bad$ledger$receipt_path[3])
bad$ledger$receipt_md5[3]<-md5(bad$ledger$receipt_path[3]);bad$expected_scope$require_new[3]<-TRUE
fail(load_analysis_release(bad$ledger,bad$expected_scope,require_projection=TRUE))
# Even rewriting a fresh projection witness and its copied parent receipt
# cannot change the scientific receipt frozen by the approved parent control.
witness_path<-r$evidence_path; original_witness<-readRDS(witness_path)
parent_path<-file.path(root,"projected",original_witness$parent_receipt_copy)
original_parent<-readRDS(parent_path);forged_parent<-original_parent;forged_parent$origin<-"new";saveRDS(forged_parent,parent_path)
forged_witness<-original_witness;forged_witness$origin<-"new";forged_witness$parent$receipt_md5<-md5(parent_path)
saveRDS(forged_witness,witness_path);r$evidence_md5<-md5(witness_path);saveRDS(r,bad$ledger$receipt_path[3])
bad$ledger$receipt_md5[3]<-md5(bad$ledger$receipt_path[3])
fail(load_analysis_release(bad$ledger,bad$expected_scope,require_projection=TRUE))
saveRDS(original_parent,parent_path);saveRDS(original_witness,witness_path);r$evidence_md5<-md5(witness_path)
# Restore the access receipt, then change projected labels with fresh artifact
# and witness digests. The parent distance/label witness must still reject it.
r$origin<-"legacy_reuse";saveRDS(r,projected$ledger$receipt_path[3]);projected$ledger$receipt_md5[3]<-md5(projected$ledger$receipt_path[3])
path<-projected$ledger$artifact_path[3];x<-readRDS(path);x[[1]]$labels[1]<-"b";saveRDS(x,path)
projected$ledger$md5[3]<-md5(path);r$md5<-md5(path)
witness<-readRDS(r$evidence_path);witness$output_md5<-r$md5;saveRDS(witness,r$evidence_path);r$evidence_md5<-md5(r$evidence_path)
saveRDS(r,projected$ledger$receipt_path[3]);projected$ledger$receipt_md5[3]<-md5(projected$ledger$receipt_path[3])
fail(load_analysis_release(projected$ledger,projected$expected_scope,require_projection=TRUE))
inventory <- analysis_release_input_inventory(control)
stopifnot(aggregate %in% inventory$path,all(inventory$md5==vapply(inventory$path,md5,character(1))))
# Native distances are already metadata. Their projected copy preserves the
# matrix, absent stored labels and explicitly unknown source feature count.
native<-ledger[3,,drop=FALSE];native$method_key<-"PILOT_hvg2000";native$kind<-"distance_feather";native$bundle_key<-"";native$sample_id_col<-"Sample"
native$artifact_path<-file.path(root,"native.feather")
arrow::write_feather(cbind(data.frame(Sample=ids),as.data.frame(as.matrix(pb$dist_mat))),native$artifact_path)
native$md5<-md5(native$artifact_path);native$receipt_path<-file.path(root,"native.receipt.rds")
r<-readRDS(ledger$receipt_path[3]);for(field in c("method_key","kind","artifact_path","bundle_key","md5"))r[[field]]<-native[[field]]
saveRDS(r,native$receipt_path);native$receipt_md5<-md5(native$receipt_path)
native_scope<-scope[3,,drop=FALSE];native_scope$method_key<-native$method_key
native_control<-list(ledger=native,expected_scope=native_scope)
native_path<-file.path(root,"native.control.rds");saveRDS(native_control,native_path)
native_projected<-project_analysis_release(native_control,file.path(root,"native.projected"),
 c(source=source_path,runtime=runtime_path,control=native_path))
native_item<-load_analysis_release(native_projected$ledger,native_projected$expected_scope,require_projection=TRUE)[[1]]
stopifnot(is.null(native_item$result$feat_mat),is.null(native_item$result$labels),is.na(native_item$result$n_features),
 identical(as.matrix(native_item$result$dist_mat),as.matrix(pb$dist_mat)))
unlink(root,recursive=TRUE)
cat("analysis projection: PASS; source immutability, pre-read boundary, CT whitelist, distances/NULL/drop/timing identity, origin and tamper controls\n")
