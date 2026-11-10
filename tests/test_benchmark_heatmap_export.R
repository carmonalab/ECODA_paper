# Actual export block with several live devices; no production inputs/scoring.
suppressPackageStartupMessages(library(pheatmap))
run<-function(){
 root<-tempfile('heatmap-export-');dir.create(root);on.exit(unlink(root,recursive=TRUE))
 prior<-if(file.exists('Rplots.pdf'))unname(tools::md5sum('Rplots.pdf'))else NULL
 on.exit({for(i in rev(dev.list()))if(i%in%dev.list())dev.off(i);stopifnot(if(is.null(prior))!file.exists('Rplots.pdf')else identical(prior,unname(tools::md5sum('Rplots.pdf'))))},add=TRUE)
 devices<-integer()
 for(i in 1:3){pdf(file.path(root,paste0('scratch',i,'.pdf')),width=6,height=6);devices<-c(devices,dev.cur());grid::grid.newpage();grid::grid.rect()}
 dev.set(devices[1L]);stopifnot(length(dev.list())==3L)
 lines<-readLines('notebooks/benchmark_analysis.rmd');starts<-grep('^```\\{r',lines);first<-starts[18L];last<-first+which(lines[(first+1):length(lines)]=='```')[1];expressions<-as.list(parse(text=lines[(first+1):(last-1)]))
 start<-which(vapply(expressions,function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('bassez_heatmap')),logical(1)))
 finish<-which(vapply(expressions,function(x)is.call(x)&&identical(x[[1]],as.name('tryCatch')),logical(1)))
 stopifnot(length(start)==1L,length(finish)==1L,finish-start==3L)
 block<-expressions[start:finish]
 e<-new.env(parent=.GlobalEnv);e$path_plots<-root;e$ds<-'Bassez'
 e$df_heatmap<-matrix(sin(seq_len(40*26)/5),40,26,dimnames=list(paste0('ct',1:40),paste0('sample',1:26)))
 e$my_sample_col<-data.frame(subtype=factor(rep(c('E','NE'),c(9,17))),row.names=colnames(e$df_heatmap))
 # Internal pheatmap device operations must finish BEFORE the named target opens.
 e$pheatmap<-function(...){stopifnot(!file.exists(file.path(root,'Heatmap_Bassez.pdf')));pheatmap::pheatmap(...)}
 before<-digest::digest(list(e$df_heatmap,e$my_sample_col));eval(block[[1L]],e)
 result<-e$bassez_heatmap;devices_after_build<-dev.list()
 for(x in block[-1L])eval(x,e)
 stopifnot(identical(before,digest::digest(list(e$df_heatmap,e$my_sample_col))),identical(dev.list(),devices_after_build),!e$bassez_heatmap_device%in%dev.list())
 stopifnot(!any(vapply(.Devices,function(x)identical(attr(x,'filepath'),file.path(root,'Heatmap_Bassez.pdf')),logical(1))))
 reference<-pheatmap::pheatmap(e$df_heatmap,annotation_col=e$my_sample_col,scale='none',clustering_method='ward.D2',show_colnames=FALSE,cutree_rows=4,cutree_cols=4,silent=TRUE)
 stopifnot(identical(result$tree_row,reference$tree_row),identical(result$tree_col,reference$tree_col),identical(result$gtable$layout,reference$gtable$layout))
 matrix_grob<-function(z)z$gtable$grobs[[which(z$gtable$layout$name=='matrix')]]$children[[1L]]
 actual_cells<-matrix_grob(result);expected_cells<-matrix_grob(reference)
 for(field in c('x','y','width','height','just','hjust','vjust','gp'))stopifnot(identical(actual_cells[[field]],expected_cells[[field]]))
 row_labels<-function(z)z$gtable$grobs[[which(z$gtable$layout$name=='row_names')]]$label
 stopifnot(identical(row_labels(result),row_labels(reference)))
 file<-file.path(root,'Heatmap_Bassez.pdf');stopifnot(file.exists(file),file.info(file)$size>1000)
 info<-system2('/opt/homebrew/bin/pdfinfo',file,stdout=TRUE,stderr=TRUE)
 pages<-grep('^Pages:',info,value=TRUE)
 stopifnot(is.null(attr(info,'status')),length(pages)==1L,identical(as.integer(sub('^Pages: +','',pages)),1L))
 raster<-file.path(root,'page')
 stopifnot(system2('/opt/homebrew/bin/pdftoppm',c('-f','1','-singlefile','-scale-to','300','-png',file,raster),stdout=FALSE,stderr=FALSE)==0L)
 pixels<-png::readPNG(paste0(raster,'.png'));stopifnot(mean(pixels[,,1:3]<.95)>.1)
 # A drawing failure still closes the exact owned target, leaving other devices alive.
 e$ds<-'Bassez_error';for(x in block[2:3])eval(x,e)
 error_device<-e$bassez_heatmap_device;other_devices<-dev.list()[dev.list()!=error_device]
 bad<-block[[4L]];bad[[2L]][[3L]]<-quote(stop('synthetic drawing failure'))
 failure<-tryCatch({eval(bad,e);NULL},error=identity)
 stopifnot(inherits(failure,'error'),identical(conditionMessage(failure),'synthetic drawing failure'),!error_device%in%dev.list(),identical(dev.list(),other_devices))
 cat('PASS: actual grob-before-device export yields one nonblank PDF page; exact target closes on success/failure; other post-build devices and scientific inputs/cells/labels/clustering unchanged; no scoring\n')
}
run()
