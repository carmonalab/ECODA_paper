# Pure plotting contracts: no models, cohort inputs, or score recomputation.
source('src/utils/evaluation_scores.R')
source('src/notebooks/benchmark_analysis.R')
reject <- function(expr) stopifnot(inherits(tryCatch({ force(expr); NULL }, error=identity), 'error'))
metrics <- evaluation_metrics()
scores <- expand.grid(dataset=c('b','a'),trans_method=c('CLR','ALR_minCVref'),
 score=metrics,stringsAsFactors=FALSE)
scores$value <- match(scores$score,metrics) / 10
scores <- scores[rev(seq_len(nrow(scores))),]
means <- benchmark_unscaled_transformation_means(scores)
stopifnot(nrow(means)==4L,all(abs(means$value-0.3)<1e-12),
 identical(means$dataset,rep(c('a','b'),2)))
reject(benchmark_unscaled_transformation_means(scores[-1,]))
reject(benchmark_unscaled_transformation_means(rbind(scores,scores[1,])))
bad <- scores; bad$value[1] <- Inf
reject(benchmark_unscaled_transformation_means(bad))
# Fits follow observed log-log slopes, including a legitimate negative slope.
x <- data.frame(method=rep(c('up','down'),each=3),n_cells=rep(c(10,100,1000),2),
 time_secs=c(100,10000,1000000,10,1,0.1))
lines <- fit_log_log_lines(x,'time_secs')
stopifnot(nrow(lines)==200L,
 max(abs(lines$fitted[lines$method=='up']-lines$n_cells[lines$method=='up']^2))<1e-8,
 max(abs(lines$fitted[lines$method=='down']-100/lines$n_cells[lines$method=='down']))<1e-10)
cat('PASS publication raw-mean completeness/pairing and observed log-log fits\n')
source('src/utils/helpers.R')
log <- data.frame(dataset=c('a','b'),method='MOFA',time_secs=c(2,3),
 mem_GB=c(NA,7),timing_source='reviewed Feather')
bundle <- data.frame(dataset=c('b','a'),method='MOFA',time_secs=c(30,20),mem_GB=c(70,5))
completed <- benchmark_complete_execution_memory(log,bundle)
stopifnot(identical(completed$time_secs,log$time_secs),identical(completed$mem_GB,c(5,7)),
 identical(completed$memory_source,c('selected result bundle','reviewed Feather')))
foreign <- bundle;foreign$method <- 'other'
stopifnot(is.na(benchmark_complete_execution_memory(log,foreign)$mem_GB[1]))
cat('PASS missing-memory fallback keeps runtime and exact dataset/method identity\n')
# Standalone pseudobulk display cost counts required prep once per method,
# while the original accounting keeps the shared charge once for the run.
make_pb <- function(method,local) execution_time_rows_from_bundle('a',method,list(
 exec_time=local,mem_GB=3,timing_schema=2L,variant_time_secs=local,
 shared_time_secs=10,shared_mem_GB=5,timing_id='run:a:view:none'))
pb_rows <- dplyr::bind_rows(make_pb('Pseudobulk_hvg2000',2),make_pb('Pseudobulk_hvg500',1),
 data.frame(dataset='b',method='Pseudobulk_hvg2000',time_secs=7,mem_GB=4))
original <- deduplicate_exec_times(pb_rows)
original$timing_source <- 'reviewed Feather'
original <- benchmark_complete_execution_memory(original,pb_rows)
before <- summarize_exec_times(original)
costs <- benchmark_method_execution_costs(original,pb_rows)
stopifnot(nrow(costs)==3L,!.execution_is_shared_timing_method(costs$method),
 all(costs$time_secs==c(12,11,7)),all(costs$mem_GB==c(5,5,4)),
 identical(summarize_exec_times(original),before),
 identical(normalize_exec_times(costs)$time_secs,costs$time_secs),
 sum(before$total_time_secs)==20)
cat('PASS standalone PB shared+variant display preserves shared-once global accounting\n')
