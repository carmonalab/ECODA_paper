source('src/5_run_benchmark_methods/run_local_ecoda_derived.R')
reject<-function(x)stopifnot(inherits(tryCatch(force(x),error=identity),'error'))
root<-tempfile('local-label-align-');dir.create(root);root<-normalizePath(root)
ids<-c('sample-C','sample-A','sample-B')
counts<-matrix(as.double(c(1,2,3,4,5,6)),3,dimnames=list(ids,c('CT1','CT2')))
labels<-setNames(factor(c('X','Y','X')),ids)
hr<-file.path(root,'HR.rds');meta<-file.path(root,'metadata.rds')
saveRDS(list(counts=counts+.5,labels=labels[c(2,3,1)]),hr)
canonical<-ids[c(2,1,3)]
saveRDS(list(labels=labels[canonical],cells_per_sample=setNames(rowSums(counts)[canonical],canonical)),meta)
config<-list(composition_path=hr,metadata_path=meta)
before<-tools::md5sum(c(hr,meta))
value<-.ecoda_load_local_composition(config)
stopifnot(identical(value$sample_ids,canonical),identical(value$counts,counts[canonical,,drop=FALSE]),
 identical(as.character(value$labels),unname(as.character(labels[canonical]))),identical(names(value$labels),canonical),
 identical(before,tools::md5sum(c(hr,meta))))
# Real label changes are rejected by the independent named metadata truth.
x<-readRDS(hr);x$labels['sample-A']<-'X';saveRDS(x,hr);reject(.ecoda_load_local_composition(config))
x<-list(counts=counts+.5,labels=labels);names(x$labels)[2]<-names(x$labels)[1];saveRDS(x,hr);reject(.ecoda_load_local_composition(config))
x<-list(counts=counts+.5,labels=labels[-1]);saveRDS(x,hr);reject(.ecoda_load_local_composition(config))
# An already aligned checkpoint keeps exactly the same canonical counts/labels.
saveRDS(list(counts=counts+.5,labels=labels),hr);aligned<-.ecoda_load_local_composition(config)
stopifnot(identical(value$counts,aligned$counts),identical(value$labels,aligned$labels))
unlink(root,recursive=TRUE)
cat('PASS exact complete-name HR label alignment, wrong/duplicate/missing labels reject, counts/source unchanged\n')
