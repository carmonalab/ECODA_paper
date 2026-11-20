# Missing historical methods and mixed contracts must not imply paired corrections.
nb <- readLines("notebooks/batch_effect_ECODA_becm_comparison.rmd")
start <- grep("^```\\{r becm-all-available-functions", nb)
stopifnot(length(start)==1L)
end <- start+which(grepl("^```[[:space:]]*$",nb[(start+1):length(nb)]))[1]
e <- new.env(parent=globalenv())
eval(parse(text=nb[(start+1):(end-1)]),e)
panel <- function(id,problem,valid=TRUE) data.frame(
  panel_id=id,dataset=id,status="Historical",problem=problem,n_samples=10L,n_types=4L,
  method_id=c("ECODA_authors_uncorrected","ECODA_authors_limma"),
  biological_anosim=c(.2,.3),batch_anosim=if(valid)c(.8,.1) else c(NA_real_,NA_real_),threshold=NA_real_)
scores <- rbind(panel("paired","Matched"),
  panel("mixed","Annotation basis differs: UNC 4 / COR 3; no limma arrow"),
  panel("missing_batch","No configured batch",FALSE))
pages <- e$becm_all_available_pages(scores,setNames(rep("all",3),unique(scores$panel_id)))
stopifnot(length(pages)==1L,nrow(attr(pages,"panel_index"))==3L,nrow(pages[[1]]$data)==4L)
arrows <- Filter(function(layer)inherits(layer$geom,"GeomSegment"),pages[[1]]$layers)
stopifnot(length(arrows)==1L,nrow(arrows[[1]]$data)==1L,
  grepl("paired",as.character(arrows[[1]]$data$dataset)),
  arrows[[1]]$data$x==.2,arrows[[1]]$data$xend==.3)
cat("Historical figure contracts passed: missing batch facet retained; mixed-contract arrow suppressed.\n")
