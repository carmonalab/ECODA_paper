# Tiny synthetic saved-result features only; never reads production RNA/H5AD.
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
for(path in c("src/utils/analysis_release.R","src/utils/evaluation_scores.R","src/utils/datasets_io.R",
 "src/utils/math_utils.R","src/utils/scoring_metrics.R","src/utils/plotting.R","src/notebooks/benchmark_analysis.R"))source(path)
fail<-function(expr)stopifnot(inherits(tryCatch(force(expr),error=base::identity),"error"))
root<-Sys.getenv("PCA_RENDER_FIXTURE_ROOT",unset="")
if(!nzchar(root))root<-tempfile("saved-pca-render-")
dir.create(root,recursive=TRUE);root<-normalizePath(root)
dir.create(file.path(root,"results"));md5<-function(path)unname(tools::md5sum(path))
before_plot<-if(file.exists("Rplots.pdf"))md5("Rplots.pdf") else NULL
grDevices::pdf(file.path(root,"scratch-device.pdf"))
datasets<-read_datasets_json("datasets.json",view="benchmark_analysis")
datasets<-datasets[vapply(datasets,function(x)isTRUE(x$use_for_benchmark),logical(1))&!grepl("_debug$",names(datasets))]
stopifnot(length(datasets)==11L)
evidence<-file.path(root,"science.txt");writeLines("Synthetic approved saved-feature result review",evidence)
ledger<-scope<-contracts<-list();index<-0L;set.seed(123)
for(dataset in names(datasets)) {
 ids<-paste0(dataset,"_s",1:8);y<-setNames(factor(rep(if(dataset=="Bassez")c("E","NE") else c("a","b"),each=4)),ids)
 counts<-matrix(sample(3:20,32,replace=TRUE),8,dimnames=list(ids,paste0("CT",1:4)))
 clr<-log(counts+.5)-rowMeans(log(counts+.5))
 genes<-c("IGHA1","IGKC","CD3D","CD8A",paste0("gene",1:16))
 rna<-matrix(rnorm(160),8,dimnames=list(ids,genes))
 bundles<-list(ECODA_authors_HR=list(feat_mat=clr,counts=counts+.5,labels=y,dist_mat=dist(clr)),
   Pseudobulk_hvg2000=list(feat_mat=rna,labels=y,dist_mat=dist(rna)),
   `__dataset_metadata__`=list(labels=y,n_samples=8L,n_cells=sum(counts),cells_per_sample=setNames(rowSums(counts),ids),n_cell_types_high_res=4L))
 if(dataset=="Kim") {
  # Reviewed partial author coverage; stored CLR axes may have another order.
  bundles$`__dataset_metadata__`$n_cells<-sum(counts)+16
  bundles$`__dataset_metadata__`$cells_per_sample<-setNames(rowSums(counts)+2,ids)
  bundles$`__dataset_metadata__`$n_cell_types_high_res<-5L
  bundles$ECODA_authors_HR$feat_mat<-clr[rev(ids),rev(colnames(clr)),drop=FALSE]
  bundles$Pseudobulk_hvg2000$feat_mat<-rna[rev(ids),,drop=FALSE]
 }
 contracts[[dataset]]<-list(sample_ids=ids,n_samples=8L,n_cell_types=4L,
  total_cells=bundles$`__dataset_metadata__`$n_cells,annotated_cells=sum(counts))
 reference<-benchmark_pca_reference(bundles$ECODA_authors_HR,bundles$`__dataset_metadata__`,contracts[[dataset]])
 stopifnot(identical(reference$ids,ids),reference$coverage$annotated_cells==sum(counts))
 if(dataset=="Kim") {
  fail(benchmark_pca_reference(bundles$ECODA_authors_HR,bundles$`__dataset_metadata__`))
  bad<-bundles$ECODA_authors_HR;bad$counts[1,1]<-bad$counts[1,1]+1
  fail(benchmark_pca_reference(bad,bundles$`__dataset_metadata__`,contracts[[dataset]]))
  bad<-bundles$ECODA_authors_HR;bad$labels[1]<-bad$labels[5]
  fail(benchmark_pca_reference(bad,bundles$`__dataset_metadata__`,contracts[[dataset]]))
  bad<-bundles$ECODA_authors_HR;rownames(bad$counts)[1]<-"foreign"
  fail(benchmark_pca_reference(bad,bundles$`__dataset_metadata__`,contracts[[dataset]]))
  bad<-bundles$`__dataset_metadata__`;names(bad$cells_per_sample)[2]<-names(bad$cells_per_sample)[1]
  fail(benchmark_pca_reference(bundles$ECODA_authors_HR,bad,contracts[[dataset]]))
 }
 if(dataset=="Bassez")bundles$Pseudobulk_hvg3000<-bundles$Pseudobulk_hvg2000
 for(method in names(bundles)) {
  index<-index+1L;value<-bundles[[method]]
  if(method!="__dataset_metadata__")value$scores<-setNames(as.list(rep(.5,8)),c(evaluation_metrics(),"mod_knnsqrtn_score","mod_knn6_score","mod_knn9_score"))
  path<-file.path(root,"results",paste0(dataset,"_",method,".rds"));saveRDS(value,path)
  row<-data.frame(dataset=dataset,view="benchmark_analysis",method_key=method,kind="standalone_result",artifact_path=path,bundle_key="",run_root=root,
   md5=md5(path),receipt_path=file.path(root,paste0(index,".receipt.rds")),receipt_md5="")
  receipt<-c(as.list(row[c("dataset","view","method_key","kind","artifact_path","bundle_key","run_root","md5")]),
   list(schema="ecoda_analysis_release_receipt_v1",status="reviewed",origin=if(dataset=="Stephenson"||dataset=="Bassez"&&method=="ECODA_authors_HR")"new" else "legacy_reuse",
    input_md5=strrep("a",32),source_policy="synthetic-current-input",scoring_policy=evaluation_scoring_policy(),evidence="reviewed",evidence_path=evidence,evidence_md5=md5(evidence)))
  saveRDS(receipt,row$receipt_path);row$receipt_md5<-md5(row$receipt_path);ledger[[index]]<-row
  expected<-row[c("dataset","view","method_key")];expected$status<-"approved";expected$require_new<-receipt$origin=="new"
  expected$current_input_md5<-receipt$input_md5;expected$source_policy<-receipt$source_policy;scope[[index]]<-expected
 }
}
control<-list(ledger=do.call(rbind,ledger),expected_scope=do.call(rbind,scope),cohort_contract=contracts);stopifnot(nrow(control$ledger)==34L)
if(Sys.getenv("PCA_RENDER_FIXTURE_ONLY",unset="0")=="1") {
 saveRDS(control,file.path(root,"control.rds"));grDevices::dev.off();quit(save="no",status=0)
}
protected<-tools::md5sum(c(control$ledger$artifact_path,control$ledger$receipt_path,evidence))
# Model fitting and H5AD/count loading are not part of display rendering.
fit_pseudobulk_deseq2<-load_h5ad_sample_aggregate<-function(...)stop("NEW MODEL/H5AD READ FORBIDDEN")
# A foreign feature ID must fail even when its receipt binds the exact bytes.
bad<-control;i<-which(bad$ledger$dataset==names(datasets)[1]&bad$ledger$method_key=="Pseudobulk_hvg2000")
value<-readRDS(bad$ledger$artifact_path[i]);rownames(value$feat_mat)[1]<-"foreign_sample"
bad$ledger$artifact_path[i]<-file.path(root,"foreign-features.rds");saveRDS(value,bad$ledger$artifact_path[i]);bad$ledger$md5[i]<-md5(bad$ledger$artifact_path[i])
receipt<-readRDS(bad$ledger$receipt_path[i]);receipt$artifact_path<-bad$ledger$artifact_path[i];receipt$md5<-bad$ledger$md5[i]
bad$ledger$receipt_path[i]<-file.path(root,"foreign-features.receipt.rds");saveRDS(receipt,bad$ledger$receipt_path[i]);bad$ledger$receipt_md5[i]<-md5(bad$ledger$receipt_path[i])
fail(render_benchmark_pca_release(bad,file.path(root,"foreign-feature-render"),datasets,strrep("b",32)))
value$feat_mat<-as.matrix(value$dist_mat)[,rev(rownames(as.matrix(value$dist_mat))),drop=FALSE]
saveRDS(value,bad$ledger$artifact_path[i]);bad$ledger$md5[i]<-md5(bad$ledger$artifact_path[i]);receipt$md5<-bad$ledger$md5[i]
saveRDS(receipt,bad$ledger$receipt_path[i]);bad$ledger$receipt_md5[i]<-md5(bad$ledger$receipt_path[i])
fail(render_benchmark_pca_release(bad,file.path(root,"distance-as-feature-render"),datasets,strrep("b",32)))
report<-render_benchmark_pca_release(control,file.path(root,"render"),datasets,strrep("b",32))
stopifnot(nrow(report$files)==19L,nrow(report$sources)==34L,report$coverage$Bassez$n_cell_types==4L,
 report$coverage$Kim$unannotated_cells==16,report$coverage$Kim$recorded_cell_types==5L,
 identical(protected,tools::md5sum(names(protected))))
path<-file.path(root,"render/render_manifest.rds")
read<-read_benchmark_pca_render(path,md5(path),names(datasets),control$expected_scope)
stopifnot(nrow(read$tables$methods)==3L,nrow(read$tables$features)==24L,nrow(read$tables$marker_coverage)==36L,
 identical(as.character(benchmark_pca_display_labels(setNames(c("NE","E"),c("x","y")),"Bassez")),c("Non-expander","Expander")))
bad<-control$expected_scope;bad$require_new[which(bad$dataset!="Bassez"&bad$dataset!="Stephenson")[1]]<-TRUE
fail(read_benchmark_pca_render(path,md5(path),names(datasets),bad))
bad_tables<-read$tables;bad_tables$features$sample1<-seq_len(nrow(bad_tables$features))
fail(validate_benchmark_pca_tables(bad_tables,report$csv_bounds))
fail(render_benchmark_pca_release(control,file.path(root,"render"),datasets,strrep("b",32)))
bad<-control;bad$expected_scope$status[1]<-"blocked";bad$expected_scope$blocked_reason<-"Unreviewed source"
fail(render_benchmark_pca_release(bad,file.path(root,"blocked"),datasets,strrep("b",32)));stopifnot(!dir.exists(file.path(root,"blocked")))
writeLines("tampered finished figure",read$files$path[1]);fail(read_benchmark_pca_render(path,md5(path),names(datasets)))
grDevices::dev.off()
stopifnot(if(is.null(before_plot))!file.exists("Rplots.pdf") else identical(before_plot,md5("Rplots.pdf")))
lines<-readLines("notebooks/benchmark_analysis.rmd");active<-FALSE;n<-0L
supp21<-match("## Supp fig 21",lines)
supp21_chunk<-lines[seq.int(supp21+1L,length(lines))][grepl("^```\\{r",lines[seq.int(supp21+1L,length(lines))])][1]
manual_hitme<-match("# Most highly variable cell types per dataset",lines)
manual_hitme_chunk<-lines[seq.int(manual_hitme+1L,length(lines))][grepl("^```\\{r",lines[seq.int(manual_hitme+1L,length(lines))])][1]
stopifnot(!grepl("eval=FALSE",supp21_chunk),grepl("eval=FALSE",manual_hitme_chunk))
for(line in lines) {
 if(startsWith(line,"```{r")){active<-TRUE;code<-character();next}
 if(active&&startsWith(line,"```")){invisible(parse(text=code));n<-n+1L;active<-FALSE;next}
 if(active)code<-c(code,line)
}
stopifnot(!active,n==62L)
unlink(root,recursive=TRUE)
cat("PASS: exact11 saved-feature PCA PDFs, Bassez pair/six bounded CSVs/current coverage, inherited source origins, no model/H5AD reads, immutable/tamper checks;62 chunks parsed\n")
