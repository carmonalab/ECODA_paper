# Build only the three actual plotting expressions; no scoring, notebook execution or production output.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(forcats);library(ggplot2);library(patchwork)})
source('src/utils/constants.R');source('src/utils/evaluation_scores.R')
run<-function(){
 lines<-readLines('notebooks/benchmark_analysis.rmd');starts<-grep('^```\\{r',lines)
 chunk<-function(i){first<-starts[i];last<-first+which(lines[(first+1):length(lines)]=='```')[1];parse(text=lines[(first+1):(last-1)])}
 pick<-function(i,needle){expressions<-as.list(chunk(i));match<-vapply(expressions,function(x)grepl(needle,paste(deparse(x),collapse=' '),fixed=TRUE),logical(1));stopifnot(sum(match)==1L);expressions[[which(match)]]}
 datasets<-c('Adams','Bassez','GongSharma','Kfoury','Kim','Lee','Pelka','Smillie','Stephenson','Wu','Zhang')
 frame<-expand.grid(dataset=datasets,method=c('m1','m2'),score=evaluation_metrics(),stringsAsFactors=FALSE)
 frame$value<-seq_len(nrow(frame))/nrow(frame);frame$trans_method<-frame$method;frame$metric<-frame$score
 e<-new.env(parent=.GlobalEnv);e$df_plot<-e$df_normalization<-e$df_zeroimp<-frame;e$pointsize<-1.5
 before<-digest::digest(list(e$df_plot,e$df_normalization,e$df_zeroimp))
 # The native score table has BOTH columns: bare metric is masked by its data column.
 metric<-'sil_score';stopifnot(nrow(filter(frame,score==metric))==nrow(frame),nrow(filter(frame,score==.env$metric))==22L)
 eval(pick(8L,'extra_metric_plots <-'),e)
 plots<-c(e$extra_metric_plots,list(eval(pick(48L,'ggplot(df_normalization, aes(x = trans_method, y = value))'),e),eval(pick(49L,'ggplot(df_zeroimp, aes(x = method, y = value))'),e)))
 prior<-if(file.exists('Rplots.pdf'))unname(tools::md5sum('Rplots.pdf'))else NULL
 scratch<-tempfile(fileext='.pdf');pdf(scratch);on.exit({dev.off();unlink(scratch);stopifnot(if(is.null(prior))!file.exists('Rplots.pdf')else identical(prior,unname(tools::md5sum('Rplots.pdf'))))},add=TRUE)
 for(i in seq_along(plots)){
  p<-plots[[i]];layer<-which(vapply(p$layers,function(x)inherits(x$geom,'GeomPoint'),logical(1)));stopifnot(length(layer)==1L)
  warnings<-character();built<-withCallingHandlers(ggplot_build(p),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  points<-built$data[[layer]];stopifnot(nrow(points)==if(i<=2L)22L else 110L,!anyNA(points$shape),length(unique(points$shape))==11L,
   all(is.finite(points$x)),all(is.finite(points$y)),!any(grepl('shape palette|Removed',warnings)))
  source<-if(i<=2L)filter(frame,score==c('sil_score','lisi_score')[i])else frame
  stopifnot(identical(sort(points$y),sort(source$value)))
  # Red control: removing only the manual scale drops five cohorts via NA shapes.
  bad<-p;bad$scales$scales<-Filter(function(s)!'shape'%in%s$aesthetics,bad$scales$scales)
  red<-suppressWarnings(ggplot_build(bad))$data[[layer]]
  stopifnot(sum(is.na(red$shape))==if(i<=2L)10L else 50L)
 }
 stopifnot(identical(before,digest::digest(list(e$df_plot,e$df_normalization,e$df_zeroimp))))
 # Test the actual three label pipelines, retaining low-resolution and explicit HR aliases.
 label_lines<-lines[grepl('stringr::str_replace("^',lines,fixed=TRUE)]
 stopifnot(length(label_lines)==3L)
 for(i in seq_along(label_lines)){
  name<-c('PILOT','QOT','PILOT-GM-VAE')[i]
  keys<-paste0(name,c('_hvg2000','_hvg2000_lowres','_hvg2000_highres','_hvg2000_HR'))
  text<-gsub(' %>%$','',trimws(label_lines[i]))
  actual<-eval(parse(text=paste0('stringr::str_replace_all(keys,c(highres="HR",lowres="LR")) %>% ',text)))
  stopifnot(identical(actual,paste0(name,c('_hvg2000_HR','_hvg2000_LR','_hvg2000_HR','_hvg2000_HR'))))
 }
 # Actual statistics expressions, with complete paired cohorts and values at/beyond the old limits.
 suppressPackageStartupMessages(library(ggpubr))
 methods<-c('CLR','Frequency(%)','Counts','ALR_minCVref','ALR_randref','Arcsine_sqrt')
 screen<-expand.grid(dataset=datasets,method=methods,score=evaluation_metrics(),stringsAsFactors=FALSE)
 screen$value<- -0.18 + (seq_len(nrow(screen))%%11L)/10*1.18
 screen$trans_method<-factor(screen$method,levels=methods)
 e$df_normalization<-screen;e$df_normalization_anosim<-filter(screen,score=='anosim_score')
 z<-screen;z$method<-factor(ifelse(z$method=='CLR','counts_all_0.5',as.character(z$method)))
 e$df_zeroimp<-z;e$my_comparisons<-list(c('CLR','Frequency(%)'),c('CLR','Counts'))
 stat_plots<-list()
 for(index in c(48L,49L)){
  expressions<-as.list(chunk(index))
  for(x in expressions){
   text<-paste(deparse(x),collapse=' ')
   if(grepl('stat_summary(',text,fixed=TRUE)&&!grepl('group_by(',text,fixed=TRUE)){
    stat_plots[[length(stat_plots)+1L]]<-eval(x,e)
   }
  }
 }
 stopifnot(length(stat_plots)==6L)
 for(p in stat_plots){
  stopifnot(is.null(p$scales$get_scales('y')))
  # Every method/metric statistical comparison receives all eleven paired observations.
  groups<-interaction(p$data$score,p$data[[if('trans_method'%in%names(p$data)&&!is.null(p$mapping$x)&&grepl('trans_method',paste(deparse(p$mapping$x),collapse=' ')))'trans_method'else'method']],drop=TRUE)
  stopifnot(all(table(groups)==11L))
  warnings<-character()
  built<-withCallingHandlers({b<-ggplot_build(p);ggplot_gtable(b);b},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  point<-which(vapply(p$layers,function(x)inherits(x$geom,'GeomPoint'),logical(1)))
  stopifnot(length(point)==1L,nrow(built$data[[point]])==nrow(p$data),all(is.finite(built$data[[point]]$y)),
    identical(sort(built$data[[point]]$y),sort(p$data$value)),!any(grepl('Removed|Computation failed|same length',warnings)))
  bars<-which(vapply(p$layers,function(x)inherits(x$geom,'GeomBar'),logical(1)))
  stopifnot(length(bars)==1L,all(is.finite(built$data[[bars]]$y)))
 }
 # The combined figure preserves the standalone2.5in/7.5in panel ratio and all data.
 e$p1<-stat_plots[[1L]];e$p2<-stat_plots[[3L]];e$path_plots<-tempdir()
 source_panels<-digest::digest(list(e$p1$data,e$p2$data));combined_expressions<-as.list(chunk(50L))
 stopifnot(length(combined_expressions)==2L)
 combined<-eval(combined_expressions[[1L]],e)
 stopifnot(identical(combined$patches$layout$widths,c(1,3)),
   identical(combined[[1L]]$data,e$p1$data),identical(combined[[2L]]$data,e$p2$data))
 capture<-NULL;e$ggsave<-function(filename,width,height){capture<<-list(filename=filename,width=width,height=height)}
 eval(combined_expressions[[2L]],e)
 stopifnot(capture$width==10,capture$height==4.5,basename(capture$filename)=='Supp_fig_16AB_trans_zero_anosim.pdf',
   identical(source_panels,digest::digest(list(e$p1$data,e$p2$data))))
 cat('PASS: actual combinedA/B uses1:3panelwidths and10x4.5canvas with identical sourcepanels\n')
 cat('PASS: metric data-mask regression; bare HR/LR labels; six actual statistics panels retain all paired11 observations and annotations; no scoring\n')
 cat('PASS: four actual panels retain all11 cohorts/finitepoints/unchangedvalues; removed-scale redcontrols lose5cohorts; no scoring\n')
}
run()
