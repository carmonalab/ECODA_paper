# Helpers for notebooks/benchmark_analysis.rmd.
# Source into the notebook evaluation environment after load_all_functions.R;
# scope helpers resolve the notebook's selected methods, datasets, and metrics
# when called, not when this file is sourced.

# Analysis scope
# Scope helpers (closed over analysis_methods, resolved after df_benchmark is
# built; only called from figure chunks, which run later).
methods_in_scope <- function(x) {
  requested <- if (exists("analysis_requested_methods", inherits = TRUE)) analysis_requested_methods else analysis_methods
  selected <- if (is.null(requested)) x else intersect(x, requested)
  if (exists("analysis_scope", inherits = TRUE) && identical(analysis_scope, "preprint"))
    selected <- selected[!grepl(extended_method_regex, selected)]
  missing <- setdiff(selected, analysis_methods)
  if (length(missing)) stop("Required figure methods are absent from the reviewed release: ", paste(missing, collapse = ", "))
  selected
}
methods_in_scope_named <- function(m) m[names(m) %in% methods_in_scope(names(m))]

# The reviewed control names every artifact and its expected current-input
# policy. Historical aggregates and new standalone shards use the same loader;
# absent required rows never fall back to a directory or another generation.
load_benchmark_analysis_release <- function(control, datasets, score_fn = calc_sep_score) {
  if (!is.list(control) || !all(c("ledger", "expected_scope", "reference_method_key") %in% names(control)))
    stop("Benchmark release requires ledger, expected_scope and explicit reference_method_key")
  selected <- names(datasets)
  references <- control$reference_method_key
  if (!is.character(references) || is.null(names(references)) || anyDuplicated(names(references)) ||
      !all(selected %in% names(references)) || anyNA(references) || any(!nzchar(references)))
    stop("Each selected dataset requires its explicit reviewed primary reference method")
  scope <- control$expected_scope
  ledger <- control$ledger
  scope <- scope[scope$dataset %in% selected & scope$view == "benchmark_analysis", , drop = FALSE]
  ledger <- ledger[ledger$dataset %in% selected & ledger$view == "benchmark_analysis", , drop = FALSE]
  for (dataset in selected) {
    bl <- scope$dataset == dataset & scope$method_key == "Pseudobulk_hvg2000_bl"
    if (sum(bl) != 1L || !isTRUE(scope$require_new[bl]))
      stop("Every benchmark dataset requires its new pseudobulk BL result: ", dataset)
  }
  validate_analysis_release(ledger, scope, require_projection = TRUE)
  reserved <- c("__dataset_metadata__", "__sample_metadata__", "__trans__", "__zeroimp__")
  result_list <- list(bmark = list(), trans = list(), zeroimp = list())
  scores <- diagnostics <- provenance <- statistics <- list()
  diagnostic_metrics <- c("mod_knnsqrtn_score", "mod_knn6_score", "mod_knn9_score")
  for (dataset in selected) {
    take <- ledger$dataset == dataset & !ledger$method_key %in% c("__trans__", "__zeroimp__")
    rows <- load_analysis_release(ledger[take, , drop = FALSE],
      scope[match(.analysis_key(ledger[take, , drop = FALSE]), .analysis_key(scope)), , drop = FALSE],
      require_projection = TRUE)
    keys <- vapply(rows, function(x) x$source$method_key, character(1))
    reference_key <- references[[dataset]]
    if (!reference_key %in% keys || grepl("NULL", reference_key) || grepl("^scITD", reference_key))
      stop("Primary reference must be a reviewed full-sample non-NULL result: ", dataset)
    reference <- rows[[match(reference_key, keys)]]
    labels <- reference$result$labels
    if (is.null(labels) || is.null(names(labels)) || anyDuplicated(names(labels)) || anyNA(labels) ||
        anyNA(names(labels)) || any(!nzchar(names(labels)))) stop("Invalid reviewed reference labels: ", dataset)
    if (any(vapply(rows, function(x) !identical(x$identity$input_md5, reference$identity$input_md5), logical(1))))
      stop("Benchmark rows disagree on the approved current input: ", dataset)
    axes <- datasets[[dataset]]$evaluation_labels
    if (length(axes) != 1L || !identical(unname(unlist(axes)), datasets[[dataset]]$label_col))
      stop("Benchmark reference labels support exactly the configured scalar primary axis")
    primary <- names(axes)[1L]
    evaluations <- setNames(list(labels), primary)
    metadata_index <- match("__dataset_metadata__", keys)
    if (is.na(metadata_index)) stop("Reviewed dataset metadata is required for counts/runtime figures: ", dataset)
    metadata <- rows[[metadata_index]]$result
    ids <- names(labels)
    if (is.null(metadata$labels) || !setequal(names(metadata$labels), ids) ||
        !identical(as.character(metadata$labels[ids]), as.character(labels)) ||
        !is.numeric(metadata$n_cells) || length(metadata$n_cells) != 1L || !is.finite(metadata$n_cells) || metadata$n_cells <= 0 ||
        !identical(as.integer(metadata$n_samples), as.integer(length(labels))) ||
        !setequal(names(metadata$cells_per_sample), ids) || any(!is.finite(metadata$cells_per_sample)) ||
        any(metadata$cells_per_sample <= 0) || sum(metadata$cells_per_sample) != metadata$n_cells ||
        !is.numeric(metadata$n_cell_types_high_res) || length(metadata$n_cell_types_high_res) != 1L ||
        !is.finite(metadata$n_cell_types_high_res) || metadata$n_cell_types_high_res < 1)
      stop("Reviewed dataset metadata disagrees with the reference sample universe: ", dataset)
    statistics[[dataset]] <- data.frame(dataset = dataset, n_samples = length(ids), n_cells = metadata$n_cells,
      cells_per_sample_mean = mean(metadata$cells_per_sample), cells_per_sample_sd = stats::sd(metadata$cells_per_sample),
      n_cell_types = metadata$n_cell_types_high_res)
    for (i in seq_along(rows)) {
      item <- rows[[i]]
      method <- unname(keys[i])
      if (method %in% reserved) {
        next
      }
      if (grepl("^ECODA_PB_combo_", method)) stop("Legacy ECODA_PB_combo methods are outside publication scope")
      calculated <- NULL
      scorer <- function(distance, labels) {
        calculated <<- score_fn(distance, labels)
        calculated
      }
      is_null <- grepl("(^|_)NULL($|_)", method)
      evaluated <- evaluate_release_result(item$result, dataset, "benchmark_analysis", method,
        evaluations, primary, is_null = is_null, score_fn = scorer,
        stored_policy = item$receipt$scoring_policy, source_identity = item$identity,
        allowed_dropped_sample_ids = if (is.null(item$receipt$allowed_dropped_sample_ids)) character()
                                    else item$receipt$allowed_dropped_sample_ids)
      bundle <- item$result
      core <- setNames(as.list(evaluated$scores$value), evaluated$scores$metric)
      bundle$scores <- utils::modifyList(if (is.null(bundle$scores)) list() else bundle$scores,
                                       if (is.null(calculated)) core else utils::modifyList(calculated, core))
      bundle$labels <- evaluated$labels[[primary]]
      result_list$bmark[[dataset]][[method]] <- bundle
      frame <- evaluated$scores
      frame$method <- frame$method_key
      frame$score <- frame$metric
      frame$n_samples <- length(bundle$labels)
      frame$n_features <- if (is.null(bundle$feat_mat)) {
        if (is.null(bundle$n_features)) NA_integer_ else bundle$n_features
      } else ncol(bundle$feat_mat)
      scores[[length(scores) + 1L]] <- frame
      provenance[[paste(dataset, method, sep = "/")]] <- evaluated$provenance
      # Keep the existing neighborhood diagnostics explicitly separate from
      # the five-metric aggregate. Never let a missing diagnostic alter a mean.
      diagnostic <- intersect(diagnostic_metrics, names(bundle$scores))
      if (length(diagnostic)) {
        extra <- frame[rep(1L, length(diagnostic)), , drop = FALSE]
        extra$metric <- extra$score <- diagnostic
        extra$value <- vapply(bundle$scores[diagnostic], function(x) {
          if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) stop("Nonfinite neighborhood diagnostic")
          x
        }, numeric(1))
        diagnostics[[length(diagnostics) + 1L]] <- extra
      }
    }
    rm(rows, reference, metadata, item, bundle, evaluated)
    invisible(gc(FALSE))
  }
  frame <- do.call(rbind, scores)
  expected_grid <- frame[c("dataset", "view", "evaluation_id", "method_key", "metric")]
  # Validate every five-metric result now; individual figures also verify the
  # rectangular set of methods they compare before calculating means.
  for (method in unique(frame$method_key)) {
    selected_rows <- frame$method_key == method
    aggregate_evaluation_scores(frame[selected_rows, ], expected_grid[selected_rows, ])
  }
  list(results = result_list, scores = rbind(frame, if (length(diagnostics)) do.call(rbind, diagnostics)),
       core_scores = frame, expected_grid = expected_grid, statistics = do.call(rbind, statistics),
       provenance = provenance, control = list(ledger = ledger, expected_scope = scope))
}

read_benchmark_analysis_control <- function(path, md5) {
  path <- .analysis_checked_file(path, md5)
  if (file.info(path)$size > 16 * 1024^2) stop("Analysis control must be small metadata")
  control <- readRDS(path)
  .analysis_checked_file(path, md5)
  # PCA and metric release receipts describe different purposes. Join their
  # original bytes through checked projection metadata, never RNA features.
  ledger <- control$ledger
  selected <- which(ledger$view == "benchmark_analysis" &
    (ledger$method_key %in% c("ECODA_authors_HR", "__dataset_metadata__", "Pseudobulk_hvg2000") |
      (ledger$dataset == "Bassez" & ledger$method_key == "Pseudobulk_hvg3000")))
  parents <- rep(NA_character_, nrow(control$expected_scope))
  origins <- rep(NA_character_, nrow(control$expected_scope))
  for (i in selected) {
    row <- ledger[i, , drop = FALSE]
    receipt_path <- .analysis_checked_file(row$receipt_path, row$receipt_md5)
    if (file.info(receipt_path)$size > 16 * 1024^2) stop("PCA parent receipt must be small metadata")
    receipt <- readRDS(receipt_path)
    .analysis_checked_file(receipt_path, row$receipt_md5)
    if (is.null(receipt$projection_schema)) next
    index <- which(.analysis_key(control$expected_scope) == .analysis_key(row))
    fields <- c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")
    if (length(index) != 1L || control$expected_scope$status[index] != "approved" ||
        !identical(receipt$schema, "ecoda_analysis_release_receipt_v1") || !identical(receipt$status, "reviewed") ||
        !identical(receipt$input_md5, control$expected_scope$current_input_md5[index]) ||
        !identical(receipt$source_policy, control$expected_scope$source_policy[index]) ||
        any(!vapply(fields, function(field) identical(receipt[[field]], as.character(row[[field]])), logical(1))))
      stop("PCA parent projection differs from the reviewed analysis scope")
    proof <- .analysis_projection_proof(row, receipt)
    if (!is.character(proof$parent$md5) || length(proof$parent$md5) != 1L ||
        is.na(proof$parent$md5) || !grepl("^[a-f0-9]{32}$", proof$parent$md5)) stop("Invalid original PCA artifact MD5")
    parents[index] <- proof$parent$md5
    origins[index] <- proof$origin
  }
  if (any(!is.na(parents))) {
    control$expected_scope$source_artifact_md5 <- parents
    control$expected_scope$source_artifact_origin <- origins
  }
  control
}

read_benchmark_execution_times <- function(control) {
  timing <- control$execution_times
  if (!is.list(timing) || !all(c("artifact_path", "md5") %in% names(timing)))
    stop("Explicit reviewed execution_times artifact_path and md5 are required; no legacy timing fallback")
  path <- .analysis_checked_file(timing$artifact_path, timing$md5)
  frame <- as.data.frame(arrow::read_feather(path))
  .analysis_checked_file(path, timing$md5)
  normalize_exec_times(frame)
}

# Finished RNA PCA figures are rendered on HPC from reviewed result features.
# Their local counterpart reads only this finite figure/summary manifest.
benchmark_pca_presentation <- function(plot) {
  percent <- function(label) {
    match <- regexpr("[0-9]+(?:\\.[0-9]+)?(?=%)", label, perl=TRUE)
    if (match[1] < 0L) stop("Missing PCA variance axis label")
    regmatches(label, match) <- format(signif(as.numeric(regmatches(label, match)), 2L),
      trim=TRUE, scientific=FALSE)
    label
  }
  plot <- plot + ggplot2::labs(x=percent(plot$labels$x), y=percent(plot$labels$y))
  for (i in seq_along(plot$layers)) if (inherits(plot$layers[[i]]$geom, "GeomTextRepel")) {
    # Rebuild the native text layer; feature/sample coordinates and arrows stay exact.
    layer <- plot$layers[[i]]; data <- layer$data
    data$name <- vapply(as.character(data$name), function(label)
      paste(strwrap(label, width=22L), collapse="\n"), character(1))
    params <- utils::modifyList(c(layer$geom_params, layer$aes_params), list(size=3,
      box.padding=grid::unit(0.4, "lines"), point.padding=grid::unit(0.2, "lines"),
      force=3, force_pull=0.3, max.time=2, max.iter=20000, min.segment.length=0, seed=123))
    plot$layers[[i]] <- do.call(ggrepel::geom_text_repel,
      c(list(data=data, mapping=layer$mapping, inherit.aes=layer$inherit.aes,
        show.legend=layer$show.legend), params))
  }
  plot
}

benchmark_pca_render_scope <- function(control, datasets) {
  selected <- names(datasets)
  if (length(selected) != 11L || anyDuplicated(selected) || !"Bassez" %in% selected || any(grepl("_debug$", selected)))
    stop("Saved PCA rendering requires the explicit eleven-cohort benchmark scope")
  expected <- expand.grid(dataset=selected, view="benchmark_analysis",
    method_key=c("ECODA_authors_HR", "Pseudobulk_hvg2000", "__dataset_metadata__"), stringsAsFactors=FALSE)
  expected <- rbind(expected, data.frame(dataset="Bassez",view="benchmark_analysis",method_key="Pseudobulk_hvg3000"))
  keys <- .analysis_key(expected)
  if (anyDuplicated(.analysis_key(control$ledger)) || anyDuplicated(.analysis_key(control$expected_scope)))
    stop("Duplicate saved PCA source keys are prohibited")
  ledger <- control$ledger[match(keys,.analysis_key(control$ledger)),,drop=FALSE]
  scope <- control$expected_scope[match(keys,.analysis_key(control$expected_scope)),,drop=FALSE]
  if (anyNA(ledger$method_key) || anyNA(scope$method_key) ||
      any(!ledger$kind %in% c("rds_bundle","standalone_result")) ||
      any(grepl("/pseudobulks/|\\.h5ad$",ledger$artifact_path)))
    stop("Saved PCA inputs require every declared result/metadata row; caches and H5ADs are prohibited")
  validate_analysis_release(ledger,scope)
  list(ledger=ledger,expected_scope=scope)
}

benchmark_pca_csv_columns <- function() list(
  methods=c("method","n_samples","n_features","n_E","n_NE","n_samples_overlap","n_samples_missing","sample_scope","ANOSIM","mod_knn3","ARI","Silhouette","LISI"),
  features=c("method","feature","variance","variance_rank","mean_E","mean_NE","delta_E_minus_NE","group_R2","loading_PC1","loading_PC2","contribution_PC1_percent","contribution_PC2_percent","contribution_PC1_PC2_percent","group_contrast_percent"),
  pcs=c("method","PC","variance_percent","group_R2"),
  proportions=c("method","feature","median_raw_percent_E","median_raw_percent_NE","delta_median_raw_percent_E_minus_NE","prevalence_E","prevalence_NE","zero_fraction_E","zero_fraction_NE","denominator","n_denominator_categories","pseudocount_removed","denominator_min_E","denominator_max_E","denominator_min_NE","denominator_max_NE"),
  marker_coverage=c("method","marker","present","variance","variance_rank","interpretation"),
  sensitivity=c("method","comparison","ANOSIM_R","ANOSIM_R_min","ANOSIM_R_max","n_samples","n_features","n_features_removed","n_evaluable","status","IG_variance_percent","IG_PC1_contribution_percent","IG_PC2_contribution_percent","retained_raw_cells_median_E","retained_raw_cells_median_NE","PC1_log_retained_cells_correlation","rarefaction_depth","seed"))

benchmark_pca_render_files <- function(datasets) c(paste0("benchmark_pca_",datasets,".pdf"),
  "Bassez_feature_PCA_ECODA_HR_PB2000.pdf","Bassez_feature_PCA_ECODA_HR_PB2000.png",
  paste0("Bassez_",names(benchmark_pca_csv_columns()),".csv"))

benchmark_pca_display_labels <- function(labels,dataset) {
  if(dataset=="Bassez") {
    map<-c(E="Expander",NE="Non-expander");values<-unname(map[as.character(labels)])
    if(anyNA(values))stop("Unexpected Bassez expansion label")
    return(setNames(factor(values,levels=unname(map)),names(labels)))
  }
  if(dataset=="Zhang") {labels<-factor(labels);levels(labels)<-tools::toTitleCase(levels(labels))}
  labels
}

validate_benchmark_pca_tables <- function(tables,bounds) {
  schemas<-benchmark_pca_csv_columns()
  if(!setequal(names(tables),names(schemas)))stop("Missing saved PCA summary tables")
  for(name in names(schemas)) {
    frame<-tables[[name]]
    if(!is.data.frame(frame)||!identical(names(frame),schemas[[name]])||!nrow(frame)||
       nrow(frame)>bounds[[name]]||any(!vapply(frame,function(x)is.atomic(x)&&is.null(dim(x)),logical(1))))
      stop("Invalid/broad PCA summary CSV schema: ",name)
  }
  for(column in c("ANOSIM","mod_knn3","ARI","Silhouette","LISI"))
    if(!is.numeric(tables$methods[[column]])||any(!is.finite(tables$methods[[column]])))stop("Missing full-space PCA summary metric")
  invisible(TRUE)
}

benchmark_pca_reference <- function(hr, metadata, contract=NULL) {
  valid_ids<-function(x)!is.null(x)&&length(x)>0L&&!anyNA(x)&&all(nzchar(x))&&!anyDuplicated(x)
  ids<-unname(names(metadata$cells_per_sample))
  if(!valid_ids(ids)||!valid_ids(names(metadata$labels))||!valid_ids(names(hr$labels))||
     !setequal(names(metadata$labels),ids)||!setequal(names(hr$labels),ids)||
     anyNA(metadata$labels)||anyNA(hr$labels)||any(!nzchar(as.character(metadata$labels)))||any(!nzchar(as.character(hr$labels)))||
     !identical(as.character(metadata$labels[ids]),as.character(hr$labels[ids]))||
     !identical(as.integer(metadata$n_samples),as.integer(length(ids))))
    stop("PCA reference/metadata sample identity mismatch")
  counts<-as.matrix(hr$counts)
  if(!is.numeric(counts)||!valid_ids(rownames(counts))||!valid_ids(colnames(counts))||
     !setequal(rownames(counts),ids)||any(!is.finite(counts)))stop("Invalid PCA author count axes")
  counts<-counts[ids,,drop=FALSE]-.5
  totals<-as.numeric(metadata$cells_per_sample[ids]);annotated<-rowSums(counts)
  keep<-colSums(counts)>0
  invalid<-tolower(trimws(colnames(counts)))%in%c("","na","nan","none","<na>","unknown","n/a","null","unassigned")
  if(any(counts<0)||any(counts!=floor(counts))||any(annotated<=0)||
     length(totals)!=length(ids)||any(!is.finite(totals))||any(totals<=0)||any(totals!=floor(totals))||
     any(annotated>totals)||length(metadata$n_cells)!=1L||!is.finite(metadata$n_cells)||
     sum(totals)!=metadata$n_cells||!any(keep)||any(invalid&keep)||
     length(metadata$n_cell_types_high_res)!=1L||!is.finite(metadata$n_cell_types_high_res)||
     metadata$n_cell_types_high_res<sum(keep))stop("PCA author counts disagree with current metadata")
  shape<-c(n_samples=length(ids),n_cell_types=sum(keep),total_cells=sum(totals),annotated_cells=sum(annotated))
  # Missing author annotations are recorded coverage, not missing samples. A
  # partial table needs the exact reviewed count contract, never a guessed total.
  if(is.null(contract)) {
    if(any(annotated!=totals)||metadata$n_cell_types_high_res!=sum(keep))stop("Partial PCA annotation coverage requires its reviewed count contract")
  } else if(!identical(as.character(contract$sample_ids),ids)||
            !identical(as.numeric(unlist(contract[names(shape)])),as.numeric(shape)))
    stop("PCA reviewed count contract mismatch")
  list(ids=ids,labels=hr$labels[ids],counts=counts[,keep,drop=FALSE],
    coverage=list(n_samples=length(ids),n_cells=sum(totals),n_cell_types=sum(keep),
      annotated_cells=sum(annotated),unannotated_cells=sum(totals)-sum(annotated),
      recorded_cell_types=metadata$n_cell_types_high_res,empty_categories=colnames(counts)[!keep]))
}

render_benchmark_pca_release <- function(control,output_root,datasets,control_md5) {
  selected<-benchmark_pca_render_scope(control,datasets)
  if(.analysis_exists(output_root)||!grepl("^/",output_root))stop("PCA render root is create-only")
  if(!dir.create(output_root))stop("Cannot create PCA render root")
  output_root<-normalizePath(output_root,mustWork=TRUE)
  publish<-function(writer,name) {
    path<-file.path(output_root,name);tmp<-tempfile(".pca-render-",tmpdir=output_root,fileext=paste0(".",tools::file_ext(name)))
    on.exit(unlink(tmp));writer(tmp)
    if(!file.link(tmp,path))stop("Existing PCA output must be preserved: ",path)
  }
  source_rows<-list();coverage<-list();bounds<-NULL
  set.seed(123)
  for(dataset in names(datasets)) {
    take<-selected$ledger$dataset==dataset
    rows<-load_analysis_release(selected$ledger[take,,drop=FALSE],selected$expected_scope[take,,drop=FALSE])
    methods<-vapply(rows,function(x)x$source$method_key,character(1));names(rows)<-methods
    hr<-rows$ECODA_authors_HR$result;metadata<-rows$`__dataset_metadata__`$result
    reference<-benchmark_pca_reference(hr,metadata,control$cohort_contract[[dataset]])
    ids<-reference$ids;counts<-reference$counts;hr$labels<-reference$labels
    coverage[[dataset]]<-reference$coverage
    bundles<-list()
    for(method in setdiff(methods,"__dataset_metadata__")) {
      item<-rows[[method]];bundle<-item$result;x<-as.matrix(bundle$feat_mat)
      if(!is.numeric(x)||is.null(rownames(x))||anyNA(rownames(x))||anyDuplicated(rownames(x))||
         !setequal(rownames(x),ids)||!ncol(x)||is.null(colnames(x))||anyNA(colnames(x))||
         any(!nzchar(colnames(x)))||anyDuplicated(colnames(x))||any(!is.finite(x))||
         (nrow(x)==ncol(x)&&setequal(colnames(x),ids))||
         is.null(names(bundle$labels))||anyDuplicated(names(bundle$labels))||!setequal(names(bundle$labels),ids)||
         !identical(as.character(bundle$labels[ids]),as.character(hr$labels)))
        stop("A genuine, sample-aligned feature matrix is required for PCA: ",dataset,"/",method)
      calculated<-NULL;scorer<-function(distance,labels){calculated<<-calc_sep_score(distance,labels,parallel=1L);calculated}
      evaluated<-evaluate_release_result(bundle,dataset,"benchmark_analysis",method,list(primary=hr$labels),"primary",
        source_identity=item$identity,stored_policy=item$receipt$scoring_policy,score_fn=scorer)
      five<-setNames(as.list(evaluated$scores$value),evaluated$scores$metric)
      bundle$scores<-utils::modifyList(if(is.null(bundle$scores))list() else bundle$scores,
        if(is.null(calculated))five else utils::modifyList(calculated,five))
      # Display views align by names after scoring the unchanged source geometry.
      bundle$feat_mat<-x[ids,,drop=FALSE];bundle$labels<-bundle$labels[ids]
      bundle$dist_mat<-stats::as.dist(as.matrix(evaluation_distance(bundle$dist_mat))[ids,ids,drop=FALSE])
      if(!is.null(bundle$counts))bundle$counts<-bundle$counts[ids,,drop=FALSE]
      bundles[[method]]<-bundle
    }
    panels<-lapply(c("Pseudobulk_hvg2000","ECODA_authors_HR"),function(method) {
      p<-plot_pca(bundles[[method]]$feat_mat,benchmark_pca_display_labels(hr$labels,dataset),
        cluster_score=FALSE,mod_score=FALSE,anosim_score=FALSE,coord_equal=FALSE,
        n_ct_show=5,invisible="quali",repel=TRUE,
        title=if(method=="ECODA_authors_HR")"Cell type composition (ECODA)" else "Pseudobulk")
      if(dataset=="Bassez")p<-p+ggplot2::labs(color="Clonal expansion",shape="Clonal expansion")
      if(dataset=="Zhang")p<-p+ggplot2::labs(color="Tissue",shape="Tissue")
      benchmark_pca_presentation(p)
    })
    paired<-patchwork::wrap_plots(panels,nrow=1,guides="collect")
    publish(function(path)ggplot2::ggsave(path,paired,width=8,height=3.5),paste0("benchmark_pca_",dataset,".pdf"))
    if(dataset=="Bassez") {
      summaries<-bassez_feature_diagnostics(bundles)
      metric_columns<-c(Silhouette="sil_score",LISI="lisi_score")
      for(column in names(metric_columns)) {
        field<-metric_columns[[column]]
        summaries$methods[[column]]<-vapply(summaries$methods$method,function(method)bundles[[method]]$scores[[field]],numeric(1))
      }
      summaries$proportions$denominator<-"Recovered raw cells in current approved author-HR categories; total/coverage verified against current metadata"
      tables<-summaries[names(benchmark_pca_csv_columns())]
      bounds<-c(methods=3L,features=sum(vapply(bundles[c("ECODA_authors_HR","Pseudobulk_hvg2000")],function(x)ncol(x$feat_mat),integer(1))),
        pcs=2L*(length(ids)-1L),proportions=ncol(counts),marker_coverage=36L,sensitivity=5L)
      validate_benchmark_pca_tables(tables,bounds)
      for(name in names(tables))publish(function(path)utils::write.csv(tables[[name]],path,row.names=FALSE),paste0("Bassez_",name,".csv"))
      panels<-lapply(c("ECODA_authors_HR","Pseudobulk_hvg2000"),function(method)
        benchmark_pca_presentation(plot_pca(bundles[[method]]$feat_mat,benchmark_pca_display_labels(hr$labels,dataset),scale.=FALSE,
          n_ct_show=7,invisible="quali",repel=TRUE,coord_equal=FALSE,title=if(method=="ECODA_authors_HR")"ECODA authors HR" else "Pseudobulk HVG2000",
          cluster_score=FALSE,mod_score=FALSE,sil_score=FALSE,anosim_score=FALSE) +
          ggplot2::labs(color="Clonal expansion",shape="Clonal expansion")))
      paired<-patchwork::wrap_plots(panels,nrow=1,guides="collect")+patchwork::plot_annotation(title="Bassez: centered, unscaled feature PCA",
        subtitle="E/NE denotes expansion, not clinical response; arrows selected without labels")
      for(ext in c("pdf","png"))publish(function(path)ggplot2::ggsave(path,paired,width=12,height=6,dpi=300),paste0("Bassez_feature_PCA_ECODA_HR_PB2000.",ext))
    }
    for(item in rows)source_rows[[length(source_rows)+1L]]<-data.frame(dataset=dataset,method_key=item$source$method_key,
      artifact_path=item$source$artifact_path,artifact_md5=item$source$md5,receipt_path=item$source$receipt_path,
      receipt_md5=item$source$receipt_md5,origin=item$receipt$origin,input_md5=item$identity$input_md5,source_policy=item$identity$source_policy)
    rm(rows,bundles,panels,paired,hr,metadata,counts,reference,x,bundle,item,evaluated);invisible(gc(FALSE))
  }
  filenames<-benchmark_pca_render_files(names(datasets));paths<-file.path(output_root,filenames)
  report<-list(schema="ecoda_saved_pca_render_v1",control_md5=control_md5,datasets=names(datasets),
    operation="RENDER_SAVED_FEATURES_WITHOUT_MODEL_REFIT",sources=do.call(rbind,source_rows),coverage=coverage,csv_bounds=bounds,
    files=data.frame(basename=filenames,md5=unname(tools::md5sum(paths))))
  .analysis_save_new(report,file.path(output_root,"render_manifest.rds"))
  inventory<-data.frame(path=c(paths,file.path(output_root,"render_manifest.rds")),md5=unname(tools::md5sum(c(paths,file.path(output_root,"render_manifest.rds")))))
  utils::write.table(inventory,file.path(output_root,"files.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  invisible(report)
}

read_benchmark_pca_render <- function(path,md5,datasets,expected_scope=NULL) {
  .analysis_checked_file(path,md5)
  if(file.info(path)$size>16*1024^2)stop("PCA manifest must be small metadata")
  report<-readRDS(path);.analysis_checked_file(path,md5)
  expected<-benchmark_pca_render_files(report$datasets)
  if(!identical(report$schema,"ecoda_saved_pca_render_v1")||!all(datasets%in%report$datasets)||length(report$datasets)!=11L||
     anyDuplicated(report$datasets)||!identical(report$operation,"RENDER_SAVED_FEATURES_WITHOUT_MODEL_REFIT")||
     !grepl("^[a-f0-9]{32}$",report$control_md5)||
     anyDuplicated(report$files$basename)||!setequal(report$files$basename,expected))stop("Invalid finite saved-PCA manifest")
  required_sources<-expand.grid(dataset=report$datasets,method_key=c("ECODA_authors_HR","Pseudobulk_hvg2000","__dataset_metadata__"),stringsAsFactors=FALSE)
  required_sources<-rbind(required_sources,data.frame(dataset="Bassez",method_key="Pseudobulk_hvg3000"))
  source_key<-function(x)paste(x$dataset,x$method_key,sep="\r")
  if(anyDuplicated(source_key(report$sources))||!setequal(source_key(report$sources),source_key(required_sources))||
     any(!report$sources$origin%in%c("new","legacy_reuse")))stop("Saved PCA source coverage/origins differ")
  if(!is.null(expected_scope)) {
    for(i in seq_len(nrow(report$sources))) {
      source<-report$sources[i,];row<-expected_scope[expected_scope$dataset==source$dataset&expected_scope$view=="benchmark_analysis"&expected_scope$method_key==source$method_key,,drop=FALSE]
      parent<-if("source_artifact_md5"%in%names(row))row$source_artifact_md5 else NA_character_
      known_policy<-if(source$method_key%in%c("ECODA_authors_HR","__dataset_metadata__"))
        "CURRENT_REVIEWED_AUTHOR_HR_INTEGER_COUNTS_PLUS_0.5_AND_DATASET_METADATA_V1" else
        "CURRENT_REVIEWED_SAVED_PSEUDOBULK_FEATURES_AND_NAMED_GEOMETRY_V1"
      parent_matches<-length(parent)==1L&&!is.na(parent)&&grepl("^[a-f0-9]{32}$",parent)&&identical(parent,source$artifact_md5)
      if(nrow(row)!=1L||row$status!="approved"||row$current_input_md5!=source$input_md5||
         (!identical(row$source_policy,source$source_policy)&&!(identical(source$source_policy,known_policy)&&parent_matches))||
         (any(!is.na(parent))&&!parent_matches)||
         (any(!is.na(parent))&&!identical(row$source_artifact_origin,source$origin))||
         (row$require_new&&source$origin!="new"))stop("Saved PCA source differs from the current scientific scope")
    }
  }
  report$files$path<-file.path(dirname(path),report$files$basename)
  for(i in seq_len(nrow(report$files))) .analysis_checked_file(report$files$path[i],report$files$md5[i])
  report$tables<-setNames(lapply(names(benchmark_pca_csv_columns()),function(name)
    utils::read.csv(file.path(dirname(path),paste0("Bassez_",name,".csv")),stringsAsFactors=FALSE,check.names=FALSE)),names(benchmark_pca_csv_columns()))
  validate_benchmark_pca_tables(report$tables,report$csv_bounds)
  report
}

# These supplements show raw metrics individually. Validate all five before
# display recoding or plotting; diagnostic modularities cannot replace k=3.
benchmark_validate_screen_scores <- function(frame) {
  fields <- c("dataset", "method", "score", "value")
  if (!is.data.frame(frame) || !all(fields %in% names(frame))) stop("Malformed supplemental score table")
  core <- as.data.frame(frame[fields])
  core$view <- "benchmark_analysis"
  core$evaluation_id <- "primary"
  core$method_key <- as.character(core$method)
  core$metric <- as.character(core$score)
  grid <- merge(unique(core[c("dataset", "view", "evaluation_id")]),
    expand.grid(method_key = unique(core$method_key), metric = evaluation_metrics(), stringsAsFactors = FALSE), by = NULL)
  aggregate_evaluation_scores(core, grid)
  lisi <- core$value[core$metric == "lisi_score"]
  if (any(lisi < 0 | lisi > 1)) stop("Supplemental LISI must use the higher-is-better separation scale")
  frame
}

plot_benchmark_raw_metric <- function(frame, metric, method_order = unique(as.character(frame$method)),
                                      x_label = "", reference = NULL) {
  frame <- benchmark_validate_screen_scores(frame)
  if (!metric %in% evaluation_metrics() || !setequal(method_order, unique(as.character(frame$method))))
    stop("Supplemental metric or method order differs from the complete comparison")
  frame <- frame[frame$score == metric, , drop = FALSE]
  frame$method <- factor(frame$method, levels = method_order)
  panel <- ggplot2::ggplot(frame, ggplot2::aes(x = method, y = value)) +
    ggplot2::stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
    ggplot2::stat_summary(fun = mean, geom = "bar", fill = "white", color = "black") +
    ggplot2::geom_point(ggplot2::aes(shape = dataset, color = dataset), size = 1.5) +
    ggplot2::scale_shape_manual(values = seq_along(unique(frame$dataset))) +
    ggplot2::geom_line(ggplot2::aes(group = dataset, color = dataset)) +
    ggplot2::labs(x = x_label, y = paste("Raw", score_facet_labels[[metric]])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, vjust = 1, hjust = 1))
  if (!is.null(reference)) {
    if (!reference %in% method_order) stop("Supplemental reference method is absent")
    panel <- panel + ggpubr::stat_compare_means(ref.group = reference,
      method = "wilcox.test", paired = TRUE, method.args = list(exact = FALSE),
      ggplot2::aes(label = ggplot2::after_stat(p.signif)), size = 3)
  }
  panel
}

benchmark_reviewed_screen_rows <- function(release, analysis, datasets) {
  if (!analysis %in% c("trans", "zeroimp") || !is.list(release$control) || !length(datasets) || anyDuplicated(datasets))
    stop("Explicit reviewed supplemental release and dataset scope are required")
  key <- paste0("__", analysis, "__")
  ledger <- release$control$ledger
  ledger <- ledger[ledger$method_key == key & ledger$view == "benchmark_analysis" & ledger$dataset %in% datasets, , drop = FALSE]
  if (anyDuplicated(ledger$dataset) || !setequal(ledger$dataset, datasets))
    stop("Every selected dataset requires its explicit reviewed ", key, " ledger row")
  ledger[match(datasets, ledger$dataset), , drop = FALSE]
}

benchmark_load_screen_row <- function(release, row) {
  scope <- release$control$expected_scope
  scope <- scope[match(.analysis_key(row), .analysis_key(scope)), , drop = FALSE]
  item <- load_analysis_release(row, scope, require_projection = TRUE)[[1L]]
  if (row$method_key == "__trans__" && !identical(item$receipt$origin, "new"))
    stop("Transformation screens require the reviewed new producer; no historical summary fallback")
  item
}

benchmark_zeroimp_scores <- function(release, datasets) {
  items <- benchmark_reviewed_screen_rows(release, "zeroimp", datasets)
  frames <- outcomes <- list(); typed <- logical(nrow(items))
  for (i in seq_len(nrow(items))) {
    item <- benchmark_load_screen_row(release, items[i, , drop = FALSE])
    methods <- item$result
    typed[i] <- !is.null(methods$schema)
    if (typed[i]) {
      .analysis_zeroimp_screen(methods)
      if (!identical(item$receipt$origin, "new") || !identical(methods$dataset, item$source$dataset))
        stop("Counts-only screen requires its new reviewed cohort producer")
      table <- methods$results
      outcomes[[i]] <- data.frame(dataset = item$source$dataset, table, check.names = FALSE)
      for (metric in evaluation_metrics()) frames[[length(frames) + 1L]] <- data.frame(
        dataset = item$source$dataset, method = table$setting, score = metric, value = table[[metric]],
        status = table$status, reason = table$reason)
      rm(item, methods)
      next
    }
    if (!is.list(methods) || !length(methods) || is.null(names(methods)) ||
        anyDuplicated(names(methods)) || any(!nzchar(names(methods)))) stop("Malformed reviewed zero-imputation results")
    for (method in names(methods)) {
      values <- methods[[method]]
      if (is.list(values) && length(values) == 1L && is.list(values[[1L]])) values <- values[[1L]]
      metrics <- evaluation_metrics()
      if (!is.list(values) || anyDuplicated(names(values)) || !all(metrics %in% names(values)))
        stop("Zero-imputation method lacks the canonical five metrics: ", method)
      values <- vapply(values[metrics], function(x) {
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) stop("Nonfinite zero-imputation metric")
        x
      }, numeric(1))
      frames[[length(frames) + 1L]] <- data.frame(dataset = item$source$dataset,
        method = method, score = metrics, value = unname(values))
    }
    rm(item, methods)
  }
  frame <- do.call(rbind, frames)
  if (!any(typed)) return(benchmark_validate_screen_scores(frame))
  if (!all(typed)) stop("Counts-only screens cannot fall back to historical cohort summaries")
  statuses <- do.call(rbind, outcomes)
  eligible <- Reduce(intersect, lapply(datasets, function(dataset)
    statuses$setting[statuses$dataset == dataset & statuses$status == "PASS"]))
  comparison <- frame[frame$method %in% eligible, , drop = FALSE]
  benchmark_validate_screen_scores(comparison)
  excluded <- setdiff(unique(statuses$setting), eligible)
  if (length(excluded)) warning(length(excluded), " of 43 zero-imputation recipes lack PASS in every selected cohort; ",
    "comparisons use the complete common PASS menu, with all failures retained in recipe_status/all_recipe_scores")
  attr(comparison, "recipe_status") <- statuses
  attr(comparison, "all_recipe_scores") <- frame
  attr(comparison, "excluded_methods") <- excluded
  attr(comparison, "comparison_policy") <- paste("Common PASS recipes; exactly", length(datasets),
    "cohorts per method/metric; all 43 outcomes per cohort retained; failed scores are NA, never imputed")
  comparison
}

benchmark_transformation_scores <- function(release, datasets) {
  if (!exists("validate_transformation_result", mode = "function") ||
      !exists("transformation_analysis_policy", mode = "function"))
    stop("Source the versioned transformation producer's pure validator before reporting")
  items <- benchmark_reviewed_screen_rows(release, "trans", datasets)
  columns <- transformation_analysis_policy()$score_columns
  frames <- list()
  for (j in seq_len(nrow(items))) {
    item <- benchmark_load_screen_row(release, items[j, , drop = FALSE])
    table <- item$result
    validate_transformation_result(table)
    for (i in seq_len(nrow(table))) {
      frames[[length(frames) + 1L]] <- data.frame(dataset = item$source$dataset,
        method = table$trans_method[i], score = names(columns),
        value = vapply(unname(columns), function(column) table[[column]][i], numeric(1)))
    }
    rm(item, table)
  }
  benchmark_validate_screen_scores(do.call(rbind, frames))
}

# Diagnostic only: average the five original metric values within each
# dataset/transformation, without changing the canonical scaled aggregates.
benchmark_unscaled_transformation_means <- function(scores) {
  metrics <- evaluation_metrics()
  required <- c("dataset", "trans_method", "score", "value")
  if (!all(required %in% names(scores))) stop("Missing transformation mean columns")
  scores <- scores[scores$score %in% metrics, required, drop = FALSE]
  if (!nrow(scores) || anyNA(scores) || any(!is.finite(scores$value)) ||
      anyDuplicated(scores[c("dataset", "trans_method", "score")]))
    stop("Transformation mean requires unique finite five-metric rows")
  datasets <- sort(unique(scores$dataset))
  methods <- unique(as.character(scores$trans_method))
  for (method in methods) for (dataset in datasets) {
    rows <- scores$dataset == dataset & scores$trans_method == method
    if (!setequal(scores$score[rows], metrics))
      stop("Incomplete transformation five-metric grid: ", dataset, "/", method)
  }
  result <- stats::aggregate(value ~ dataset + trans_method, scores, mean)
  # Paired tests must see the same dataset order within every method.
  result[order(match(result$trans_method, methods), result$dataset), , drop = FALSE]
}

# Derived-run loading
# Strict loader for local, run-owned derived analyses.  A derived RDS is never
# trusted from its filename alone: the completed run record, source manifest,
# output record, output checksum, and (for ordinary H5AD mode) source H5AD
# digest are all checked before readRDS() is called.  Snapshot mode verifies
# the compact local composition snapshot and preserves remote H5AD paths only
# as provenance.
.derived_require_scalar <- function(x, label, allow_na = FALSE) {
  if (!is.character(x) || length(x) != 1L ||
      (!allow_na && (is.na(x) || !nzchar(x)))) {
    stop("Invalid ", label, " in derived-run manifest.")
  }
  x
}

.derived_resolve_path <- function(path, base_dir = NULL) {
  .derived_require_scalar(as.character(path), "path")
  if (!is.null(base_dir) &&
      !grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path <- file.path(base_dir, path)
  }
  normalizePath(path, mustWork = FALSE)
}

.derived_path_within <- function(path, root) {
  path <- .derived_resolve_path(path)
  root <- .derived_resolve_path(root)
  boundary <- if (identical(root, .Platform$file.sep)) {
    root
  } else {
    paste0(root, .Platform$file.sep)
  }
  identical(path, root) || startsWith(path, boundary)
}

.derived_run_artifact_path <- function(run_root, relative_path,
                                       label = "derived artifact") {
  root <- .derived_resolve_path(run_root)
  relative_path <- .derived_require_scalar(
    as.character(relative_path),
    paste0(label, " relative path")
  )
  if (!dir.exists(root) ||
      grepl("^(/|[A-Za-z]:[/\\\\])", relative_path, perl = TRUE) ||
      grepl("^[A-Za-z]:", relative_path, perl = TRUE) ||
      relative_path %in% c(".", "..") ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", relative_path, perl = TRUE) ||
      grepl("[\\\\]", relative_path, perl = TRUE) ||
      grepl("[[:cntrl:]]", relative_path, perl = TRUE)) {
    stop(label, " is not a safe path within the requested run root.")
  }
  path <- .derived_resolve_path(file.path(root, relative_path))
  if (!.derived_path_within(path, root) || identical(path, root)) {
    stop(label, " escapes the requested run root: ", path)
  }
  path
}

.derived_dataset_artifact_path <- function(run_root, dataset, suffix,
                                           label = "derived artifact") {
  dataset <- .derived_require_scalar(
    as.character(dataset),
    paste0(label, " dataset")
  )
  suffix <- .derived_require_scalar(
    as.character(suffix),
    paste0(label, " suffix")
  )
  if (dataset %in% c(".", "..") ||
      grepl("^[A-Za-z]:", dataset, perl = TRUE) ||
      grepl("[/\\\\]", dataset, perl = TRUE) ||
      grepl("[[:cntrl:]]", dataset, perl = TRUE) ||
      grepl("^(/|[A-Za-z]:[/\\\\])", suffix, perl = TRUE) ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", suffix, perl = TRUE) ||
      grepl("[\\\\]", suffix, perl = TRUE) ||
      grepl("[[:cntrl:]]", suffix, perl = TRUE)) {
    stop(label, " has an unsafe dataset/path component.")
  }
  .derived_run_artifact_path(
    run_root,
    paste0(dataset, suffix),
    label
  )
}

.derived_md5 <- function(path, label = "artifact") {
  if (!file.exists(path) || !isTRUE(file.info(path)$isdir == FALSE) ||
      !isTRUE(file.info(path)$size > 0)) {
    stop("Missing or empty ", label, ": ", path)
  }
  digest <- tolower(unname(tools::md5sum(path)))
  if (length(digest) != 1L || is.na(digest) ||
      !grepl("^[[:xdigit:]]{32}$", digest)) {
    stop("Could not compute a valid MD5 for ", label, ": ", path)
  }
  digest
}

.derived_read_json <- function(path, label) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    stop("Missing or empty ", label, ": ", path)
  }
  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) {
      stop("Invalid ", label, ": ", path, ": ", conditionMessage(e))
    }
  )
}

.derived_require_fields <- function(x, fields, label) {
  if (!is.list(x) || !all(fields %in% names(x))) {
    stop(label, " is missing required fields: ",
         paste(setdiff(fields, names(x)), collapse = ", "))
  }
  invisible(x)
}

.derived_sidecar <- function(path, require_path = TRUE) {
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || !isTRUE(file.info(sidecar)$size > 0)) {
    stop("Missing derived-artifact checksum sidecar: ", sidecar)
  }
  lines <- readLines(sidecar, warn = FALSE)
  if (length(lines) != 3L || any(!nzchar(lines)) ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    stop("Malformed derived-artifact checksum sidecar: ", sidecar)
  }
  md5 <- sub("^MD5=", "", lines[[1L]])
  size <- sub("^SIZE=", "", lines[[2L]])
  recorded_path <- sub("^PATH=", "", lines[[3L]])
  if (!grepl("^[[:xdigit:]]{32}$", md5) ||
      !grepl("^[0-9]+$", size) || !nzchar(recorded_path)) {
    stop("Malformed derived-artifact checksum sidecar: ", sidecar)
  }
  actual_size <- as.character(file.info(path)$size)
  actual_md5 <- .derived_md5(path)
  if (!identical(tolower(md5), actual_md5) ||
      !identical(size, actual_size)) {
    stop("Derived-artifact checksum sidecar mismatch: ", path)
  }
  if (isTRUE(require_path)) {
    path_matches <- identical(
      .derived_resolve_path(recorded_path),
      .derived_resolve_path(path)
    )
  } else {
    path_matches <- identical(basename(recorded_path), basename(path))
  }
  if (!path_matches) {
    stop("Derived-artifact checksum sidecar PATH mismatch: ", path)
  }
  list(md5 = actual_md5, size = as.numeric(actual_size))
}
.derived_snapshot_context <- function(
  snapshot_path_ref,
  snapshot_md5_ref,
  snapshot_size_ref = NULL,
  run_dir,
  expected_scope = NULL,
  expected_config_md5 = NULL
) {
  snapshot_path <- .derived_resolve_snapshot_file(
    snapshot_path_ref,
    run_dir,
    "composition snapshot"
  )
  snapshot_md5_ref <- as.character(snapshot_md5_ref)
  if (length(snapshot_md5_ref) != 1L ||
      !grepl("^[[:xdigit:]]{32}$", snapshot_md5_ref)) {
    stop("Invalid composition snapshot MD5 in derived manifest.")
  }
  snapshot_sidecar <- .derived_sidecar(snapshot_path, require_path = FALSE)
  if (!identical(snapshot_sidecar$md5, tolower(snapshot_md5_ref))) {
    stop("Composition snapshot checksum mismatch: ", snapshot_path)
  }
  if (!is.null(snapshot_size_ref)) {
    snapshot_size_ref <- suppressWarnings(as.numeric(snapshot_size_ref))
    if (length(snapshot_size_ref) != 1L ||
        !is.finite(snapshot_size_ref) ||
        snapshot_size_ref <= 0 ||
        snapshot_size_ref != snapshot_sidecar$size) {
      stop("Composition snapshot size mismatch: ", snapshot_path)
    }
  }

  snapshot_manifest_path <- .derived_resolve_path(
    file.path(dirname(snapshot_path), "composition_snapshot_manifest.json")
  )
  snapshot <- .derived_read_json(
    snapshot_manifest_path,
    "composition_snapshot_manifest.json"
  )
  .derived_require_fields(
    snapshot,
    c(
      "schema_version", "run_id", "scope", "status", "source_mode",
      "config_path", "config_md5", "extractor_path", "extractor_md5",
      "snapshot_path", "snapshot_md5", "snapshot_size", "sources"
    ),
    "composition_snapshot_manifest.json"
  )
  schema_version <- suppressWarnings(as.numeric(snapshot$schema_version))
  if (length(schema_version) != 1L || !is.finite(schema_version) ||
      schema_version != 1) {
    stop("Composition snapshot schema_version is invalid.")
  }
  .derived_require_scalar(
    as.character(snapshot$run_id),
    "composition snapshot run_id"
  )
  .derived_require_scalar(
    as.character(snapshot$scope),
    "composition snapshot scope"
  )
  if (!identical(as.character(snapshot$status), "COMPLETED")) {
    stop("Composition snapshot manifest is not a completed snapshot.")
  }
  if (!identical(as.character(snapshot$source_mode), "composition_snapshot")) {
    stop("Composition snapshot manifest source_mode is invalid.")
  }
  if (!is.null(expected_scope) &&
      !identical(as.character(snapshot$scope), as.character(expected_scope))) {
    stop("Composition snapshot scope does not match the derived run.")
  }
  for (field in c("config_path", "extractor_path")) {
    .derived_require_scalar(
      as.character(snapshot[[field]]),
      paste0("composition snapshot ", field)
    )
  }
  if (!grepl("^/", as.character(snapshot$extractor_path))) {
    stop("Composition snapshot extractor_path must be absolute.")
  }

  for (field in c("config_md5", "extractor_md5", "snapshot_md5")) {
    if (!grepl("^[[:xdigit:]]{32}$", as.character(snapshot[[field]]))) {
      stop("Invalid composition snapshot ", field, ".")
    }
  }
  if (!is.null(expected_config_md5) &&
      !identical(
        tolower(as.character(snapshot$config_md5)),
        tolower(as.character(expected_config_md5))
      )) {
    stop("Composition snapshot config checksum does not match the derived run.")
  }
  manifest_snapshot_md5 <- tolower(as.character(snapshot$snapshot_md5))
  if (!identical(manifest_snapshot_md5, snapshot_sidecar$md5) ||
      !identical(manifest_snapshot_md5, tolower(snapshot_md5_ref))) {
    stop("Composition snapshot manifest checksum does not match the RDS.")
  }
  manifest_snapshot_size <- suppressWarnings(as.numeric(snapshot$snapshot_size))
  if (length(manifest_snapshot_size) != 1L ||
      !is.finite(manifest_snapshot_size) ||
      manifest_snapshot_size <= 0 ||
      manifest_snapshot_size != snapshot_sidecar$size) {
    stop("Composition snapshot manifest size does not match the RDS.")
  }
  manifest_snapshot_path <- .derived_require_scalar(
    as.character(snapshot$snapshot_path),
    "composition snapshot snapshot_path"
  )
  json_fields <- c(
    "snapshot_json_path", "snapshot_json_md5", "snapshot_json_size"
  )
  json_present <- json_fields %in% names(snapshot)
  snapshot_json_path <- NULL
  if (any(json_present) && !all(json_present)) {
    stop("Composition snapshot JSON identity is incomplete.")
  }
  if (all(json_present)) {
    snapshot_json_path <- .derived_resolve_snapshot_file(
      snapshot$snapshot_json_path,
      dirname(snapshot_path),
      "composition snapshot JSON"
    )
    json_sidecar <- .derived_sidecar(
      snapshot_json_path, require_path = FALSE
    )
    json_md5 <- tolower(as.character(snapshot$snapshot_json_md5))
    json_size <- suppressWarnings(as.numeric(snapshot$snapshot_json_size))
    if (!grepl("^[[:xdigit:]]{32}$", json_md5) ||
        !identical(json_md5, json_sidecar$md5) ||
        length(json_size) != 1L || !is.finite(json_size) ||
        json_size != json_sidecar$size) {
      stop("Composition snapshot JSON identity does not match its sidecar.")
    }
  }


  if (!is.list(snapshot$sources) || length(snapshot$sources) == 0L) {
    stop("Composition snapshot manifest has no source provenance rows.")
  }
  snapshot_sources <- lapply(seq_along(snapshot$sources), function(i) {
    row <- snapshot$sources[[i]]
    .derived_require_fields(
      row,
      c("dataset", "view", "h5ad_path", "h5ad_md5"),
      paste0("composition snapshot source entry ", i)
    )
    dataset <- .derived_require_scalar(
      as.character(row$dataset),
      paste0("composition snapshot source ", i, " dataset")
    )
    view <- .derived_require_scalar(
      as.character(row$view),
      paste0("composition snapshot source ", i, " view")
    )
    if (!identical(view, "benchmark_analysis")) {
      stop("Composition snapshot source entry ", i,
           " is not the benchmark_analysis view.")
    }
    h5ad_path <- as.character(row$h5ad_path)
    if (length(h5ad_path) != 1L || is.na(h5ad_path) ||
        !grepl("^/", h5ad_path) || !nzchar(h5ad_path)) {
      stop("Composition snapshot source entry ", i,
           " must retain an absolute HPC H5AD path.")
    }
    h5ad_path <- .derived_resolve_path(h5ad_path)
    h5ad_md5 <- tolower(as.character(row$h5ad_md5))
    if (length(h5ad_md5) != 1L ||
        !grepl("^[[:xdigit:]]{32}$", h5ad_md5)) {
      stop("Invalid H5AD MD5 in composition snapshot source entry ", i)
    }
    get_source_field <- function(candidates, label) {
      present <- candidates[candidates %in% names(row)]
      if (!length(present)) {
        stop("Composition snapshot source entry ", i,
             " is missing ", label, ".")
      }
      value <- row[[present[[1L]]]]
      if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
        stop("Invalid ", label, " in composition snapshot source entry ", i)
      }
      value
    }
    source_size <- suppressWarnings(as.numeric(
      get_source_field(c("h5ad_size", "size"), "H5AD size")
    ))
    if (length(source_size) != 1L || !is.finite(source_size) ||
        source_size <= 0 || source_size != floor(source_size)) {
      stop("Invalid H5AD size in composition snapshot source entry ", i)
    }
    source_mtime <- get_source_field(c("h5ad_mtime", "mtime"), "H5AD mtime")
    label_col <- as.character(get_source_field(
      c("label_col", "label_column", "label"),
      "biological-label column"
    ))
    high_col <- as.character(get_source_field(
      c(
        "high_res_col", "high_res_column", "cell_type_high_res",
        "cell_type_high_res_column"
      ),
      "high-resolution cell-type column"
    ))
    list(
      dataset = dataset,
      view = view,
      h5ad_path = h5ad_path,
      h5ad_md5 = h5ad_md5,
      size = source_size,
      mtime = as.character(source_mtime),
      label_column = label_col,
      high_res_column = high_col
    )
  })
  source_datasets <- vapply(snapshot_sources, `[[`, character(1), "dataset")
  if (anyDuplicated(source_datasets)) {
    stop("Composition snapshot source provenance contains duplicate datasets.")
  }
  list(
    path = snapshot_path,
    md5 = snapshot_sidecar$md5,
    size = snapshot_sidecar$size,
    manifest_path = snapshot_manifest_path,
    snapshot_config_path = as.character(snapshot$config_path),
    snapshot_config_md5 = tolower(as.character(snapshot$config_md5)),
    snapshot_remote_path = manifest_snapshot_path,
    snapshot_json_path = if (exists("snapshot_json_path")) {
      snapshot_json_path
    } else {
      NULL
    },
    manifest = snapshot,
    sources = snapshot_sources
  )
}
.derived_snapshot_source_field <- function(row, candidates, label) {
  present <- candidates[candidates %in% names(row)]
  if (!length(present)) {
    stop("Composition snapshot source is missing ", label, ".")
  }
  values <- vapply(present, function(name) {
    value <- row[[name]]
    if (length(value) != 1L || is.na(value) ||
        !nzchar(as.character(value))) {
      stop("Invalid ", label, " in composition snapshot source.")
    }
    as.character(value)
  }, character(1))
  if (length(unique(values)) != 1L) {
    stop("Conflicting ", label, " fields in composition snapshot source.")
  }
  values[[1L]]
}

.derived_compare_snapshot_sources <- function(snapshot_sources, source_rows) {
  if (length(snapshot_sources) != length(source_rows)) {
    stop("Derived source rows do not match the composition snapshot.")
  }
  for (i in seq_along(snapshot_sources)) {
    for (field in c(
      "dataset", "view", "h5ad_path", "h5ad_md5", "size", "mtime",
      "label_column", "high_res_column"
    )) {
      if (!field %in% names(source_rows[[i]]) ||
          !identical(
            as.character(snapshot_sources[[i]][[field]]),
            as.character(source_rows[[i]][[field]])
          )) {
        stop("Derived source row ", i,
             " does not match the composition snapshot provenance.")
      }
    }
  }
  invisible(TRUE)
}


.derived_resolve_bound_file <- function(path, run_dir, label) {
  path <- .derived_require_scalar(as.character(path), label)
  candidates <- if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path
  } else {
    c(file.path(run_dir, path), file.path(getwd(), path))
  }
  candidates <- unique(vapply(
    candidates,
    function(x) .derived_resolve_path(x),
    character(1)
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) != 1L) {
    stop("Could not resolve exactly one ", label, ": ", path)
  }
  existing[[1L]]
}

.derived_resolve_snapshot_file <- function(path_ref, run_dir, label) {
  path_ref <- .derived_require_scalar(as.character(path_ref), label)
  candidates <- if (grepl("^(/|[A-Za-z]:[/\\\\])", path_ref)) {
    path_ref
  } else {
    c(file.path(run_dir, path_ref), file.path(getwd(), path_ref))
  }
  candidates <- c(candidates, file.path(run_dir, basename(path_ref)))
  candidates <- unique(vapply(
    candidates,
    function(x) .derived_resolve_path(x),
    character(1)
  ))
  info <- file.info(candidates)
  existing <- candidates[
    file.exists(candidates) & !is.na(info$isdir) & !info$isdir
  ]
  if (length(existing) != 1L) {
    stop("Could not resolve exactly one ", label, ": ", path_ref)
  }
  existing[[1L]]
}

.derived_validate_record <- function(record, run_dir, label) {
  .derived_require_fields(
    record,
    c("path", "size", "md5", "dataset", "method"),
    label
  )
  artifact_path <- .derived_resolve_path(record$path, run_dir)
  if (!.derived_path_within(artifact_path, run_dir) ||
      identical(artifact_path, .derived_resolve_path(run_dir))) {
    stop(label, " escapes the run-owned directory: ", artifact_path)
  }
  if (!grepl("^[[:xdigit:]]{32}$", as.character(record$md5)) ||
      !grepl("^[1-9][0-9]*$", as.character(record$size))) {
    stop("Invalid checksum/size in ", label, ": ", artifact_path)
  }
  actual_size <- as.character(file.info(artifact_path)$size)
  actual_md5 <- .derived_md5(artifact_path, label)
  if (!identical(as.character(record$size), actual_size) ||
      !identical(tolower(as.character(record$md5)), actual_md5)) {
    stop("Derived-run artifact manifest mismatch: ", artifact_path)
  }
  if (!is.character(record$dataset) || length(record$dataset) != 1L ||
      is.na(record$dataset) ||
      !is.character(record$method) || length(record$method) != 1L ||
      is.na(record$method)) {
    stop("Invalid dataset/method provenance in ", label)
  }
  list(
    path = artifact_path,
    size = as.numeric(actual_size),
    md5 = actual_md5,
    dataset = record$dataset,
    method = record$method
  )
}

.derived_run_context <- function(path, expected_analysis = NULL) {
  path <- .derived_resolve_path(path)
  run_dir <- dirname(path)
  run_manifest_path <- file.path(run_dir, "derived_run_manifest.json")
  source_manifest_default <- file.path(run_dir, "derived_source_manifest.json")
  status_path <- file.path(run_dir, "status.json")
  run <- .derived_read_json(run_manifest_path, "derived_run_manifest.json")
  .derived_require_fields(
    run,
    c(
      "run_id", "analysis", "scope", "status", "output_dir", "derived_root",
      "config_path",
      "config_md5", "runner_path", "runner_md5", "source_manifest_path",
      "source_manifest_md5", "artifacts", "sources"
    ),
    "derived_run_manifest.json"
  )
  source_mode <- if (is.null(run$source_mode)) {
    "h5ad"
  } else {
    .derived_require_scalar(
      as.character(run$source_mode),
      "derived_run_manifest.json source_mode"
    )
  }
  if (!source_mode %in%
      c("h5ad", "composition_snapshot", "local_composition")) {
    stop("Unsupported derived source_mode: ", source_mode)
  }
  snapshot_context <- NULL

  run_id <- .derived_require_scalar(
    as.character(run$run_id),
    "derived_run_manifest.json run_id"
  )
  owner_path <- file.path(run_dir, ".ecoda_run_owner")
  if (!file.exists(owner_path) ||
      !isTRUE(file_test("-f", owner_path)) ||
      isTRUE(nzchar(Sys.readlink(owner_path), keepNA = TRUE))) {
    stop("Missing or invalid .ecoda_run_owner marker: ", owner_path)
  }
  owner_lines <- tryCatch(
    readLines(owner_path, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(owner_lines) || length(owner_lines) != 1L ||
      !identical(owner_lines[[1L]], paste0("RUN_ID=", run_id))) {
    stop(".ecoda_run_owner marker does not match derived run: ", owner_path)
  }

  run_derived_root <- .derived_require_scalar(
    run$derived_root,
    "derived_run_manifest.json derived_root"
  )
  if (!identical(run_derived_root, .derived_resolve_path(run_derived_root))) {
    stop("derived_run_manifest.json derived_root must be a canonical path")
  }
  if (!.derived_path_within(run_dir, run_derived_root) ||
      identical(run_dir, run_derived_root)) {
    stop(
      "Derived run directory must be a strict descendant of derived_root: ",
      run_derived_root
    )
  }
  if (!identical(as.character(run$status), "COMPLETED") ||
      (!is.null(expected_analysis) &&
       !identical(as.character(run$analysis), expected_analysis))) {
    stop(
      "Derived run is not a completed ",
      if (is.null(expected_analysis)) "" else expected_analysis,
      " run: ", run_manifest_path
    )
  }
  status <- .derived_read_json(status_path, "status.json")
  .derived_require_fields(
    status,
    c("status", "run_id", "analysis", "scope", "output_dir", "derived_root"),
    "status.json"
  )
  if (!identical(as.character(status$status), "COMPLETED")) {
    stop("Derived run terminal status is not COMPLETED: ", status_path)
  }
  status_derived_root <- .derived_require_scalar(
    status$derived_root,
    "status.json derived_root"
  )
  if (!identical(status_derived_root,
                 .derived_resolve_path(status_derived_root))) {
    stop("status.json derived_root must be a canonical path")
  }
  if (!identical(status_derived_root, run_derived_root)) {
    stop("status.json derived_root does not match derived_run_manifest.json")
  }
  for (field in c("run_id", "analysis", "scope", "output_dir")) {
    .derived_require_scalar(
      as.character(status[[field]]),
      paste0("status.json ", field)
    )
    .derived_require_scalar(
      as.character(run[[field]]),
      paste0("derived_run_manifest.json ", field)
    )
    if (!identical(as.character(status[[field]]),
                   as.character(run[[field]]))) {
      stop("status.json ", field, " does not match derived_run_manifest.json")
    }
  }
  if (!identical(.derived_resolve_path(run$output_dir), run_dir)) {
    stop(
      "derived_run_manifest.json output_dir does not match the run-owned directory"
    )
  }
  for (field in c("config_md5", "runner_md5", "source_manifest_md5")) {
    if (!grepl("^[[:xdigit:]]{32}$", as.character(run[[field]]))) {
      stop("Invalid ", field, " in derived_run_manifest.json")
    }
  }

  config_path <- .derived_resolve_bound_file(
    run$config_path, run_dir, "config_path"
  )
  runner_path <- .derived_resolve_bound_file(
    run$runner_path, run_dir, "runner_path"
  )
  if (!identical(.derived_md5(config_path, "config_path"),
                 tolower(as.character(run$config_md5))) ||
      !identical(.derived_md5(runner_path, "runner_path"),
                 tolower(as.character(run$runner_md5)))) {
    stop("Derived run config/runner checksum mismatch: ", run_manifest_path)
  }
  if (identical(source_mode, "composition_snapshot")) {
    .derived_require_fields(
      run,
      c(
        "snapshot_path", "snapshot_md5",
        "snapshot_config_path", "snapshot_config_md5"
      ),
      "composition_snapshot derived_run_manifest.json"
    )
    run_snapshot_config_path <- .derived_require_scalar(
      as.character(run$snapshot_config_path),
      "derived run snapshot_config_path"
    )
    run_snapshot_config_md5 <- .derived_require_scalar(
      as.character(run$snapshot_config_md5),
      "derived run snapshot_config_md5"
    )
    if (!grepl("^[[:xdigit:]]{32}$", run_snapshot_config_md5)) {
      stop("Invalid snapshot_config_md5 in derived_run_manifest.json")
    }
    snapshot_context <- .derived_snapshot_context(
      snapshot_path_ref = run$snapshot_path,
      snapshot_md5_ref = run$snapshot_md5,
      snapshot_size_ref = if ("snapshot_size" %in% names(run)) {
        run$snapshot_size
      } else {
        NULL
      },
      run_dir = run_dir,
      expected_scope = run$scope,
      expected_config_md5 = run$config_md5
    )
    if (!identical(
      run_snapshot_config_path,
      snapshot_context$snapshot_config_path
    ) ||
        !identical(
          tolower(run_snapshot_config_md5),
          snapshot_context$snapshot_config_md5
        )) {
      stop("Derived run snapshot config provenance does not match the snapshot.")
    }
  } else if (any(c(
    "snapshot_path", "snapshot_md5", "snapshot_size",
    "snapshot_config_path", "snapshot_config_md5"
  ) %in% names(run))) {
    stop("Snapshot identity fields require source_mode=composition_snapshot.")
  }


  source_manifest_ref <- .derived_require_scalar(
    as.character(run$source_manifest_path),
    "derived_run_manifest.json source_manifest_path"
  )
  source_manifest_path <- .derived_resolve_bound_file(
    source_manifest_ref,
    run_dir,
    "source_manifest_path"
  )
  canonical_source_manifest <- .derived_resolve_path(source_manifest_default)

  if (!identical(source_manifest_path, canonical_source_manifest)) {
    stop(
      "derived_run_manifest.json source_manifest_path must resolve to ",
      "the canonical run-root derived_source_manifest.json"
    )
  }

  if (!.derived_path_within(source_manifest_path, run_dir)) {
    stop("source_manifest_path escapes the run-owned directory.")
  }
  if (!identical(.derived_md5(source_manifest_path, "source manifest"),
                 tolower(as.character(run$source_manifest_md5)))) {
    stop("Derived source manifest checksum mismatch: ", source_manifest_path)
  }
  source <- .derived_read_json(source_manifest_path, "derived_source_manifest.json")
  .derived_require_fields(
    source,
    c(
      "run_id", "analysis", "scope", "derived_root", "config_path",
      "config_md5", "runner_path", "runner_md5", "sources"
    ),
    "derived_source_manifest.json"
  )
  source_mode_source <- if (is.null(source$source_mode)) {
    "h5ad"
  } else {
    .derived_require_scalar(
      as.character(source$source_mode),
      "derived_source_manifest.json source_mode"
    )
  }
  if (!source_mode_source %in%
      c("h5ad", "composition_snapshot", "local_composition") ||
      !identical(source_mode_source, source_mode)) {
    stop("Derived source manifest source_mode does not match the run.")
  }
  if (identical(source_mode, "composition_snapshot")) {
    .derived_require_fields(
      source,
      c(
        "snapshot_path", "snapshot_md5",
        "snapshot_config_path", "snapshot_config_md5"
      ),
      "composition_snapshot derived_source_manifest.json"
    )
    source_snapshot_path <- .derived_require_scalar(
      as.character(source$snapshot_path),
      "derived source snapshot_path"
    )
    run_snapshot_path <- .derived_require_scalar(
      as.character(run$snapshot_path),
      "derived run snapshot_path"
    )
    source_snapshot_md5 <- .derived_require_scalar(
      as.character(source$snapshot_md5),
      "derived source snapshot_md5"
    )
    source_snapshot_config_path <- .derived_require_scalar(
      as.character(source$snapshot_config_path),
      "derived source snapshot_config_path"
    )
    source_snapshot_config_md5 <- .derived_require_scalar(
      as.character(source$snapshot_config_md5),
      "derived source snapshot_config_md5"
    )
    if (!grepl("^[[:xdigit:]]{32}$", source_snapshot_md5) ||
        !grepl("^[[:xdigit:]]{32}$", source_snapshot_config_md5)) {
      stop("Invalid snapshot identity in derived_source_manifest.json.")
    }
    if (!identical(
      basename(source_snapshot_path),
      basename(run_snapshot_path)
    ) ||
        !identical(
          tolower(source_snapshot_md5),
          tolower(as.character(run$snapshot_md5))
        )) {
      stop("Derived source/run snapshot identity does not match.")
    }
    if (!identical(
      source_snapshot_config_path,
      as.character(run$snapshot_config_path)
    ) ||
        !identical(
          tolower(source_snapshot_config_md5),
          tolower(as.character(run$snapshot_config_md5))
        ) ||
        !identical(
          source_snapshot_config_path,
          snapshot_context$snapshot_config_path
        ) ||
        !identical(
          tolower(source_snapshot_config_md5),
          snapshot_context$snapshot_config_md5
        )) {
      stop("Derived source/run snapshot config provenance does not match.")
    }
    snapshot_size_matches <- function(value, label) {
      value <- suppressWarnings(as.numeric(value))
      if (length(value) != 1L || !is.finite(value) ||
          value != snapshot_context$size) {
        stop(label, " does not match the local snapshot.")
      }
      invisible(TRUE)
    }
    if ("snapshot_size" %in% names(run)) {
      snapshot_size_matches(run$snapshot_size, "Derived run snapshot size")
    }
    if ("snapshot_size" %in% names(source)) {
      snapshot_size_matches(
        source$snapshot_size,
        "Derived source snapshot size"
      )
    }
  } else if (any(c(
    "snapshot_path", "snapshot_md5", "snapshot_size",
    "snapshot_config_path", "snapshot_config_md5"
  ) %in% names(source))) {
    stop("Snapshot identity fields require source_mode=composition_snapshot.")
  }

  source_derived_root <- .derived_require_scalar(
    source$derived_root,
    "derived_source_manifest.json derived_root"
  )
  if (!identical(source_derived_root,
                 .derived_resolve_path(source_derived_root))) {
    stop("derived_source_manifest.json derived_root must be a canonical path")
  }
  for (field in c("run_id", "analysis", "scope", "config_path", "runner_path")) {
    .derived_require_scalar(
      as.character(source[[field]]),
      paste0("derived source manifest ", field)
    )
  }
  for (field in c("config_md5", "runner_md5")) {
    if (!grepl("^[[:xdigit:]]{32}$", as.character(source[[field]]))) {
      stop("Invalid ", field, " in derived_source_manifest.json")
    }
  }
  if (!identical(as.character(source$run_id), as.character(run$run_id)) ||
      !identical(as.character(source$analysis), as.character(run$analysis)) ||
      !identical(as.character(source$scope), as.character(run$scope)) ||
      !identical(source_derived_root, run_derived_root) ||
      !identical(as.character(source$config_path), as.character(run$config_path)) ||
      !identical(tolower(as.character(source$config_md5)),
                 tolower(as.character(run$config_md5))) ||
      !identical(as.character(source$runner_path), as.character(run$runner_path)) ||
      !identical(tolower(as.character(source$runner_md5)),
                 tolower(as.character(run$runner_md5)))) {
    stop("Derived source manifest run provenance does not match the run.")
  }
  if (!identical(as.character(source$analysis), as.character(run$analysis)) ||
      !is.list(source$sources) || length(source$sources) == 0L) {
    stop("Derived source manifest analysis/sources do not match the run.")
  }
  source_rows <- lapply(seq_along(source$sources), function(i) {
    row <- source$sources[[i]]
    dataset <- .derived_require_scalar(
      as.character(row$dataset),
      paste0("source entry ", i, " dataset")
    )
    if (identical(source_mode, "local_composition")) {
      local_fields <- c(
        "dataset", "view", "composition_path", "composition_md5",
        "composition_size", "composition_mtime", "metadata_path",
        "metadata_md5", "metadata_size", "metadata_mtime", "analysis"
      )
      .derived_require_fields(
        row, local_fields, paste0("local composition source entry ", i)
      )
      if (!identical(as.character(row$view), "local_composition") ||
          !identical(as.character(row$analysis), as.character(run$analysis)) ||
          !grepl("^[[:xdigit:]]{32}$", as.character(row$composition_md5)) ||
          !grepl("^[[:xdigit:]]{32}$", as.character(row$metadata_md5))) {
        stop("Invalid local composition source entry ", i)
      }
      validate_local_input <- function(
        path_ref,
        expected_md5,
        expected_size,
        expected_mtime,
        label
      ) {
        path <- .derived_resolve_path(path_ref)
        info <- file.info(path)
        if (!file.exists(path) || isTRUE(info$isdir) ||
            is.na(info$size) || info$size <= 0 || is.na(info$mtime)) {
          stop("Invalid ", label, ": ", path)
        }
        actual_md5 <- .derived_md5(path, label)
        if (!identical(actual_md5, tolower(as.character(expected_md5))) ||
            !identical(
              as.character(info$size),
              as.character(expected_size)
            ) ||
            !identical(
              as.character(as.numeric(info$mtime)),
              as.character(as.numeric(expected_mtime))
            )) {
          stop(label, " identity mismatch: ", path)
        }
        list(
          path = path,
          md5 = actual_md5,
          size = as.numeric(info$size),
          mtime = as.numeric(info$mtime)
        )
      }
      composition_identity <- validate_local_input(
        row$composition_path,
        row$composition_md5,
        row$composition_size,
        row$composition_mtime,
        "local HR composition bundle"
      )
      metadata_identity <- validate_local_input(
        row$metadata_path,
        row$metadata_md5,
        row$metadata_size,
        row$metadata_mtime,
        "local metadata bundle"
      )
      return(list(
        dataset = dataset,
        view = as.character(row$view),
        composition_path = composition_identity$path,
        composition_md5 = composition_identity$md5,
        composition_size = composition_identity$size,
        composition_mtime = composition_identity$mtime,
        metadata_path = metadata_identity$path,
        metadata_md5 = metadata_identity$md5,
        metadata_size = metadata_identity$size,
        metadata_mtime = metadata_identity$mtime,
        analysis = as.character(row$analysis)
      ))
    }
    .derived_require_fields(
      row, c("dataset", "view", "h5ad_path", "h5ad_md5", "analysis"),
      paste0("derived source entry ", i)
    )
    if (!identical(as.character(row$view), "benchmark_analysis") ||
        !identical(as.character(row$analysis), as.character(run$analysis)) ||
        !grepl("^[[:xdigit:]]{32}$", as.character(row$h5ad_md5))) {
      stop("Invalid derived source entry ", i)
    }
    h5ad_path <- as.character(row$h5ad_path)
    if (length(h5ad_path) != 1L || is.na(h5ad_path) || !nzchar(h5ad_path)) {
      stop("Invalid H5AD path in derived source entry ", i)
    }
    h5ad_path <- .derived_resolve_path(h5ad_path)
    h5ad_md5 <- tolower(as.character(row$h5ad_md5))
    if (identical(source_mode, "composition_snapshot")) {
      if (!grepl("^/", h5ad_path)) {
        stop("Snapshot source H5AD path is not absolute: ", h5ad_path)
      }
    } else {
      .derived_sidecar(h5ad_path)
      actual_h5ad_md5 <- .derived_md5(h5ad_path, "source H5AD")
      if (!identical(actual_h5ad_md5, h5ad_md5)) {
        stop("Source H5AD checksum mismatch: ", h5ad_path)
      }
      h5ad_md5 <- actual_h5ad_md5
    }
    source_row <- list(
      dataset = dataset,
      view = as.character(row$view),
      h5ad_path = h5ad_path,
      h5ad_md5 = h5ad_md5,
      analysis = as.character(row$analysis)
    )
    if (identical(source_mode, "composition_snapshot")) {
      source_size <- suppressWarnings(as.numeric(
        .derived_snapshot_source_field(
          row, c("size", "h5ad_size"), "H5AD size"
        )
      ))
      if (length(source_size) != 1L || !is.finite(source_size) ||
          source_size <= 0 || source_size != floor(source_size)) {
        stop("Invalid H5AD size in derived source entry ", i)
      }
      source_row$size <- source_size
      source_row$mtime <- as.character(.derived_snapshot_source_field(
        row, c("mtime", "h5ad_mtime"), "H5AD mtime"
      ))
      source_row$label_column <- as.character(
        .derived_snapshot_source_field(
          row, c("label_column", "label_col"), "biological-label column"
        )
      )
      source_row$high_res_column <- as.character(
        .derived_snapshot_source_field(
          row,
          c("high_res_column", "high_res_col"),
          "high-resolution cell-type column"
        )
      )
    }
    source_row
  })
  source_datasets <- vapply(source_rows, `[[`, character(1), "dataset")
  if (anyDuplicated(source_datasets)) {
    stop("Derived source manifest contains duplicate datasets.")
  }
  if (!is.list(run$sources) || length(run$sources) != length(source_rows)) {
    stop("derived_run_manifest.json sources do not match source manifest.")
  }
  if (identical(source_mode, "composition_snapshot")) {
    .derived_compare_snapshot_sources(
      snapshot_context$sources,
      source_rows
    )
  }
  for (i in seq_along(source_rows)) {
    run_source <- run$sources[[i]]
    if (identical(source_mode, "local_composition")) {
      local_fields <- c(
        "dataset", "view", "composition_path", "composition_md5",
        "composition_size", "composition_mtime", "metadata_path",
        "metadata_md5", "metadata_size", "metadata_mtime", "analysis"
      )
      .derived_require_fields(
        run_source, local_fields,
        paste0("derived run local source entry ", i)
      )
      for (field in local_fields) {
        if (!identical(
          as.character(run_source[[field]]),
          as.character(source_rows[[i]][[field]])
        )) {
          stop("Derived run/source local provenance mismatch at entry ", i)
        }
      }
    } else {
      .derived_require_fields(
        run_source,
        c("dataset", "view", "h5ad_path", "h5ad_md5", "analysis"),
        paste0("derived run source entry ", i)
      )
      for (field in c(
        "dataset", "view", "h5ad_path", "h5ad_md5", "analysis"
      )) {
        if (!identical(
          as.character(run_source[[field]]),
          as.character(source_rows[[i]][[field]])
        )) {
          stop("derived run/source manifest source mismatch at entry ", i)
        }
      }
      if (identical(source_mode, "composition_snapshot")) {
        provenance_fields <- list(
          size = c("size", "h5ad_size"),
          mtime = c("mtime", "h5ad_mtime"),
          label_column = c("label_column", "label_col"),
          high_res_column = c("high_res_column", "high_res_col")
        )
        for (field in names(provenance_fields)) {
          run_value <- .derived_snapshot_source_field(
            run_source,
            provenance_fields[[field]],
            field
          )
          if (!identical(
            as.character(run_value),
            as.character(source_rows[[i]][[field]])
          )) {
            stop("Derived run/source snapshot provenance mismatch at entry ", i)
          }
        }
      }
    }
  }

  if (!is.list(run$artifacts) || length(run$artifacts) == 0L) {
    stop("derived_run_manifest.json has no output artifact records.")
  }
  artifact_rows <- lapply(seq_along(run$artifacts), function(i) {
    .derived_validate_record(
      run$artifacts[[i]], run_dir,
      paste0("derived output artifact ", i)
    )
  })
  artifact_paths <- vapply(artifact_rows, `[[`, character(1), "path")
  if (anyDuplicated(artifact_paths)) {
    stop("Derived run manifest contains duplicate output paths.")
  }
  artifact_datasets <- vapply(artifact_rows, `[[`, character(1), "dataset")
  artifact_methods <- vapply(artifact_rows, `[[`, character(1), "method")
  if (run$analysis %in% c("res50", "harmony")) {
    expected_method <- if (identical(run$analysis, "res50")) {
      "ECODA_seuratres_50"
    } else {
      "ECODA_seuratres_2_harmony"
    }
    expected_paths <- paste0(
      source_datasets, "_", expected_method, ".rds"
    )
    if (!identical(artifact_datasets, source_datasets) ||
        !identical(
          artifact_methods,
          rep(expected_method, length(source_datasets))
        ) ||
        !identical(basename(artifact_paths), expected_paths)) {
      stop(
        "Derived run manifest does not record the exact ",
        run$analysis, " output set."
      )
    }
  } else if (identical(run$analysis, "cell_subsetting")) {
    expected_paths <- c(
      "ECODA_authors_HR_cell_subsetting.rds",
      "ECODA_authors_HR_cell_subsetting.csv"
    )
    if (!identical(artifact_datasets, c("ALL", "ALL")) ||
        !identical(
          artifact_methods,
          c(
            "ECODA_authors_HR_cell_subsetting",
            "ECODA_authors_HR_cell_subsetting"
          )
        ) ||
        !identical(basename(artifact_paths), expected_paths)) {
      stop("Derived run manifest does not record the exact cell-subsetting output set.")
    }
  }
  list(
    run_dir = run_dir,
    run = run,
    status = status,
    source = source,
    source_mode = source_mode,
    snapshot_path = if (is.null(snapshot_context)) NULL else snapshot_context$path,
    snapshot_md5 = if (is.null(snapshot_context)) NULL else snapshot_context$md5,
    snapshot_size = if (is.null(snapshot_context)) NULL else snapshot_context$size,
    snapshot_config_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_config_path
    },
    snapshot_config_md5 = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_config_md5
    },
    snapshot_remote_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_remote_path
    },
    snapshot_json_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_json_path
    },
    snapshot_manifest_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$manifest_path
    },
    snapshot_manifest = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$manifest
    },
    sources = source_rows,
    artifacts = artifact_rows,
    config_path = config_path,
    runner_path = runner_path,
    source_manifest_path = source_manifest_path,
    derived_root = run_derived_root
  )
}

load_derived_bundle <- function(
  path,
  dataset,
  method,
  expected_samples = NULL,
  allow_dropped = FALSE
) {
  if (isTRUE(allow_dropped) &&
      length(method) == 1L &&
      !is.na(method) &&
      method %in% c(
        "ECODA_seuratres_2_harmony",
        "ECODA_seuratres_50"
      )) {
    stop(
      "allow_dropped=TRUE is not permitted for complete-sample method: ",
      method
    )
  }
  path <- .derived_resolve_path(path)
  expected_analysis <- switch(
    method,
    ECODA_seuratres_2_harmony = "harmony",
    ECODA_seuratres_50 = "res50",
    NULL
  )
  context <- .derived_run_context(path, expected_analysis = expected_analysis)
  records <- Filter(function(x) identical(x$path, path), context$artifacts)
  if (length(records) != 1L) {
    stop("Derived bundle is not uniquely recorded in the run manifest: ", path)
  }
  record <- records[[1L]]
  if (!identical(record$dataset, dataset) || !identical(record$method, method)) {
    stop("Derived bundle provenance mismatch: ", path)
  }
  source_rows <- Filter(function(x) identical(x$dataset, dataset), context$sources)
  if (length(source_rows) != 1L) {
    stop("Derived bundle has no unique source H5AD entry: ", dataset)
  }
  .derived_sidecar(path)
  bundle <- tryCatch(
    readRDS(path),
    error = function(e) stop("Could not read derived bundle ", path, ": ",
                             conditionMessage(e))
  )
  required <- c("feat_mat", "labels", "dist_mat", "scores", "counts")
  if (!is.list(bundle) || !all(required %in% names(bundle))) {
    stop("Derived bundle missing required fields: ", path)
  }
  feat_mat <- bundle$feat_mat
  feat_values <- suppressWarnings(as.numeric(as.matrix(feat_mat)))
  if ((!is.matrix(feat_mat) && !is.data.frame(feat_mat)) ||
      nrow(feat_mat) == 0L || ncol(feat_mat) == 0L ||
      is.null(rownames(feat_mat)) || anyNA(rownames(feat_mat)) ||
      any(!nzchar(rownames(feat_mat))) || anyDuplicated(rownames(feat_mat)) ||
      any(!is.finite(feat_values))) {
    stop("Invalid feat_mat in derived bundle: ", path)
  }
  actual_samples <- rownames(feat_mat)
  labels <- bundle$labels
  if (length(labels) != length(actual_samples) ||
      is.null(names(labels)) || anyDuplicated(names(labels)) ||
      !identical(as.character(names(labels)), as.character(actual_samples))) {
    stop("Derived bundle labels do not match feat_mat sample IDs: ", path)
  }
  distances <- bundle$dist_mat
  distance_matrix <- as.matrix(distances)
  if (!is.matrix(distance_matrix) ||
      nrow(distance_matrix) != length(actual_samples) ||
      ncol(distance_matrix) != length(actual_samples) ||
      is.null(rownames(distance_matrix)) ||
      is.null(colnames(distance_matrix)) ||
      !identical(as.character(rownames(distance_matrix)),
                 as.character(actual_samples)) ||
      !identical(as.character(colnames(distance_matrix)),
                 as.character(actual_samples)) ||
      any(!is.finite(distance_matrix))) {
    stop("Invalid dist_mat in derived bundle: ", path)
  }
  if (!is.list(bundle$scores) && !is.numeric(bundle$scores)) {
    stop("Invalid scores in derived bundle: ", path)
  }
  score_values <- unlist(bundle$scores, use.names = FALSE)
  if (!is.numeric(score_values) ||
      length(score_values) == 0L || any(!is.finite(score_values))) {
    stop("Derived bundle scores are empty or nonfinite: ", path)
  }
  counts <- bundle$counts
  counts_values <- suppressWarnings(as.numeric(as.matrix(counts)))
  if ((!is.matrix(counts) && !is.data.frame(counts)) ||
      nrow(counts) == 0L || ncol(counts) == 0L ||
      is.null(rownames(counts)) || anyDuplicated(rownames(counts)) ||
      !identical(as.character(rownames(counts)),
                 as.character(actual_samples)) ||
      any(!is.finite(counts_values))) {
    stop("Invalid counts/sample IDs in derived bundle: ", path)
  }
  if (!is.null(expected_samples)) {
    expected_samples <- as.character(expected_samples)
    if (length(expected_samples) == 0L || anyNA(expected_samples) ||
        anyDuplicated(expected_samples)) {
      stop("Invalid expected sample IDs for derived bundle: ", path)
    }
    extras <- setdiff(actual_samples, expected_samples)
    missing <- setdiff(expected_samples, actual_samples)
    if (length(extras) > 0L || (!allow_dropped && length(missing) > 0L)) {
      stop(
        "Derived bundle sample policy failed for ", dataset, "/", method,
        "; missing=", paste(missing, collapse = ","),
        "; extras=", paste(extras, collapse = ",")
      )
    }
    if (!allow_dropped && !identical(actual_samples, expected_samples)) {
      stop("Derived bundle sample order differs from expected samples: ", path)
    }
    if (allow_dropped && length(missing) > 0L) {
      attr(bundle, "derived_dropped_samples") <- missing
    }
  }
  bundle
}

load_derived_cell_subsetting <- function(
  rds_path,
  csv_path,
  expected_datasets,
  expected_samples_by_dataset,
  release_control_path = NULL, release_control_md5 = NULL,
  expected_labels_by_dataset = NULL, expected_current_inputs = NULL
) {
  portable <- !is.null(release_control_path)
  if (portable) {
    if (!all(expected_datasets %in% names(expected_current_inputs)) ||
        !all(expected_datasets %in% names(expected_labels_by_dataset)) ||
        any(!vapply(expected_current_inputs[expected_datasets], function(x)
          is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[a-f0-9]{32}$", x), logical(1))))
      stop("Current cell-depth reading requires explicit current input hashes and named labels")
    .analysis_checked_file(release_control_path, release_control_md5)
    control <- readRDS(release_control_path)
    .analysis_checked_file(release_control_path, release_control_md5)
    if (nrow(control$expected_scope) != 1L ||
        !isTRUE(control$expected_scope$require_new) ||
        control$expected_scope$method_key != "ECODA_authors_HR_cell_subsetting")
      stop("Current cell-depth output requires one required-new reviewed artifact")
    entries <- load_analysis_release(control$ledger, control$expected_scope, require_projection = TRUE)
    entry <- entries[[1L]]; raw <- entry$result; rds_path <- entry$source$artifact_path
    if (!identical(raw$schema, "ecoda_cell_depth_five_v1") ||
        !identical(raw$scoring_policy, evaluation_scoring_policy()) ||
        !identical(raw$source_release_md5, entry$identity$input_md5))
      stop("Current cell-depth source/scoring identity mismatch")
    for (dataset in expected_datasets) {
      if (!identical(raw$current_inputs[[dataset]], expected_current_inputs[[dataset]]))
        stop("Cell-depth output does not use the current approved input: ", dataset)
    }
    name <- raw$csv$basename
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !identical(basename(name), name) || name %in% c(".", ".."))
      stop("Cell-depth CSV requires an explicit safe relative basename")
    csv_path <- file.path(dirname(rds_path), name)
    .analysis_checked_file(csv_path, raw$csv$md5)
    context <- list(sources = raw$sources, source_mode = "reviewed_screen", artifact = entry$source)
  } else {
    rds_path <- .derived_resolve_path(rds_path)
    csv_path <- .derived_resolve_path(csv_path)
    context <- .derived_run_context(rds_path, expected_analysis = "cell_subsetting")
  }
  expected_datasets <- as.character(expected_datasets)
  if (length(expected_datasets) == 0L || anyDuplicated(expected_datasets)) {
    stop("Invalid expected datasets for cell-subsetting results.")
  }
  source_datasets <- vapply(context$sources, `[[`, character(1), "dataset")
  if (length(source_datasets) == 0L || anyDuplicated(source_datasets) ||
      !setequal(source_datasets, expected_datasets)) {
    stop("Cell-subsetting source manifest does not cover the expected datasets.")
  }
  if (!portable) {
  rds_records <- Filter(function(x) identical(x$path, rds_path), context$artifacts)
  csv_records <- Filter(function(x) identical(x$path, csv_path), context$artifacts)
  if (length(rds_records) != 1L || length(csv_records) != 1L) {
    stop("Cell-subsetting RDS/CSV are not uniquely recorded in the run manifest.")
  }
  if (!identical(rds_records[[1L]]$method,
                 "ECODA_authors_HR_cell_subsetting") ||
      !identical(csv_records[[1L]]$method,
                 "ECODA_authors_HR_cell_subsetting")) {
    stop("Cell-subsetting artifact method provenance is invalid.")
  }
  .derived_sidecar(rds_path)
  .derived_sidecar(csv_path)
  raw <- tryCatch(
    readRDS(rds_path),
    error = function(e) stop("Could not read cell-subsetting RDS: ",
                             conditionMessage(e))
  )
  }
  results <- if (is.data.frame(raw)) {
    raw
  } else if (is.list(raw) && is.data.frame(raw$results)) {
    raw$results
  } else {
    stop("Cell-subsetting RDS must contain a results data.frame.")
  }
  rds_table <- results
  if (portable) {
    if (!all(c(evaluation_metrics(), "draw_id") %in% names(results)) ||
        !is.list(raw$draws) || length(raw$draws) != nrow(results) || anyDuplicated(results$draw_id))
      stop("Current cell-depth release lacks all five metrics and exact saved draws")
    for (metric in evaluation_metrics()) {
      if (!is.numeric(results[[metric]]) || any(!is.finite(results[[metric]])))
        stop("Cell-depth release has missing/nonfinite metric ", metric)
    }
    if (any(results$lisi_score < 0 | results$lisi_score > 1) ||
        !identical(results$anosim_score, results$ANOSIM)) stop("Cell-depth score policy mismatch")
    for (i in seq_along(raw$draws)) {
      draw <- raw$draws[[i]]; identity <- draw$draw_id; draw$draw_id <- NULL
      if (!identical(identity, results$draw_id[i]) || !identical(identity, .evaluation_fingerprint(draw)))
        stop("Cell-depth saved draw identity mismatch")
      for (field in c("dataset", "target", "replicate", "seed")) {
        if (!identical(as.character(draw[[field]]), as.character(results[[field]][i])))
          stop("Cell-depth saved draw schedule mismatch")
      }
      counts <- draw$counts; ids <- expected_samples_by_dataset[[draw$dataset]]
      source <- Filter(function(x) identical(x$dataset, draw$dataset), raw$sources)
      fields <- c("composition_path", "composition_md5", "composition_bundle_key", "metadata_path",
        "metadata_md5", "metadata_bundle_key", "input_md5", "h5ad_path", "h5ad_md5")
      if (length(source) != 1L || !identical(draw$source_identity,
          .evaluation_fingerprint(c(list(dataset = draw$dataset), source[[1]][intersect(fields, names(source[[1]]))]))))
        stop("Cell-depth draw does not bind its declared source")
      if (!is.matrix(counts) || !is.numeric(counts) || any(!is.finite(counts)) ||
          any(counts < 0) || any(counts != floor(counts)) || any(rowSums(counts) <= 0) ||
          !identical(rownames(counts), ids) || sum(counts) != results$effective_cells[i] ||
          !identical(as.character(draw$labels), as.character(expected_labels_by_dataset[[draw$dataset]][ids])) ||
          !identical(names(draw$labels), ids)) stop("Cell-depth draw counts/labels disagree with current cohort")
      if (draw$target == "all cells" && draw$dataset == "Bassez" &&
          (nrow(counts) != 26L || ncol(counts) != 40L || sum(counts) != 74948L))
        stop("Stale Bassez cell-depth population")
      if (draw$dataset == "Stephenson" && nrow(counts) != 53L) stop("Stale Stephenson cell-depth sample universe")
    }
  }
  required <- c(
    "dataset", "target", "target_cells", "replicate", "seed", "sample_ids",
    "original_cells_per_sample", "effective_cells_per_sample",
    "total_cells", "ANOSIM"
  )
  if (!all(required %in% names(results))) {
    stop("Cell-subsetting results missing fields: ",
         paste(setdiff(required, names(results)), collapse = ", "))
  }
  if (identical(context$source_mode, "local_composition")) {
    methodology_fields <- c(
      "source_mode", "sampling_method", "source_limitation"
    )
    if (!all(methodology_fields %in% names(results)) ||
        any(!nzchar(as.character(results$source_mode))) ||
        any(as.character(results$source_mode) != "local_composition") ||
        any(!nzchar(as.character(results$sampling_method))) ||
        any(!nzchar(as.character(results$source_limitation)))) {
      stop(
        "Local composition subsetting methodology is missing or invalid."
      )
    }
  }
  if (!setequal(unique(as.character(results$dataset)), expected_datasets)) {
    stop("Cell-subsetting results do not cover the expected datasets.")
  }
  target_levels <- c("all cells", "2000", "1000", "500", "400",
                     "300", "200", "150", "100", "50")
  results$dataset <- as.character(results$dataset)
  results$target <- as.character(results$target)
  if (any(!results$target %in% target_levels)) {
    stop("Cell-subsetting results contain an unexpected target.")
  }
  target_cells_raw <- results$target_cells
  results$target_cells <- suppressWarnings(as.numeric(target_cells_raw))
  target_cells_expected <- ifelse(
    results$target == "all cells",
    NA_real_,
    suppressWarnings(as.numeric(results$target))
  )
  if (length(results$target_cells) != nrow(results) ||
      any(is.na(results$target_cells) & !is.na(target_cells_raw)) ||
      any(is.nan(results$target_cells)) ||
      any(is.na(results$target_cells) != is.na(target_cells_expected)) ||
      any(!is.na(results$target_cells) &
          (!is.finite(results$target_cells) |
           results$target_cells < 1 |
           results$target_cells != floor(results$target_cells))) ||
      any(!is.na(results$target_cells) &
          results$target_cells != target_cells_expected)) {
    stop("Cell-subsetting target_cells values are invalid.")
  }
  if (any(!is.finite(as.numeric(results$replicate))) ||
      any(as.integer(results$replicate) != as.numeric(results$replicate)) ||
      any(!is.finite(as.numeric(results$ANOSIM)))) {
    stop("Cell-subsetting replicate/ANOSIM values are invalid.")
  }
  results$replicate <- as.integer(results$replicate)
  seed_raw <- results$seed
  results$seed <- suppressWarnings(as.numeric(seed_raw))
  if (length(results$seed) != nrow(results) ||
      any(is.na(results$seed) & !is.na(seed_raw)) ||
      any(is.nan(results$seed)) ||
      any(!is.na(results$seed) & !is.finite(results$seed))) {
    stop("Cell-subsetting seed values are invalid.")
  }
  results$anosim <- as.numeric(results$ANOSIM)
  as_ids <- function(x) {
    if (is.list(x)) {
      if (length(x) == 1L && is.list(x[[1L]])) x <- x[[1L]]
      else x <- unlist(x, use.names = FALSE)
    }
    if (length(x) == 1L && is.character(x)) {
      x <- strsplit(x, "|", fixed = TRUE)[[1L]]
    }
    trimws(as.character(x))
  }
  as_named_counts <- function(x) {
    if (is.list(x)) {
      if (length(x) == 1L && is.list(x[[1L]])) x <- x[[1L]]
      else x <- unlist(x, use.names = FALSE)
    }
    if (length(x) == 1L && is.character(x)) {
      tokens <- strsplit(x, "|", fixed = TRUE)[[1L]]
      fields <- strsplit(tokens, ":", fixed = TRUE)
      if (any(lengths(fields) != 2L)) return(numeric())
      values <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1), 2L)))
      names(values) <- trimws(vapply(fields, `[[`, character(1), 1L))
      return(values)
    }
    as.numeric(x)
  }
  sample_lists <- lapply(results$sample_ids, as_ids)
  original_count_lists <- lapply(
    results$original_cells_per_sample,
    as_named_counts
  )
  count_lists <- lapply(
    results$effective_cells_per_sample,
    as_named_counts
  )
  for (i in seq_len(nrow(results))) {
    ds <- results$dataset[[i]]
    expected_ids <- as.character(expected_samples_by_dataset[[ds]])
    ids <- sample_lists[[i]]
    original_counts <- original_count_lists[[i]]
    counts <- count_lists[[i]]
    if (length(expected_ids) == 0L || anyDuplicated(expected_ids) ||
        !identical(ids, expected_ids) || anyDuplicated(ids) ||
        length(original_counts) != length(ids) ||
        length(counts) != length(ids) ||
        !identical(names(original_counts), ids) ||
        !identical(names(counts), ids) ||
        any(!is.finite(original_counts)) || any(!is.finite(counts)) ||
        any(original_counts < 1) || any(counts < 1) ||
        any(original_counts != as.integer(original_counts)) ||
        any(counts != as.integer(counts)) ||
        !isTRUE(as.numeric(results$total_cells[[i]]) == sum(counts))) {
      stop("Cell-subsetting sample/count diagnostics are invalid at row ", i)
    }
  }
  for (ds in expected_datasets) {
    ds_rows <- results[results$dataset == ds, , drop = FALSE]
    if (!identical(sort(unique(ds_rows$target)), sort(target_levels))) {
      stop("Cell-subsetting targets are incomplete for dataset ", ds)
    }
    baseline <- ds_rows[ds_rows$target == "all cells", , drop = FALSE]
    if (nrow(baseline) != 1L || baseline$replicate[[1L]] != 0L ||
        !is.na(baseline$seed[[1L]])) {
      stop("Cell-subsetting baseline policy failed for dataset ", ds)
    }
    baseline_i <- which(rownames(results) == rownames(baseline)[[1L]])
    baseline_original <- original_count_lists[[baseline_i]]
    baseline_counts <- count_lists[[baseline_i]]
    if (!identical(baseline_counts, baseline_original)) {
      stop("Cell-subsetting all-cells baseline changed effective counts.")
    }
    for (target in target_levels[-1L]) {
      rows <- ds_rows[ds_rows$target == target, , drop = FALSE]
      if (nrow(rows) != 20L ||
          !setequal(rows$replicate, seq_len(20L)) ||
          !setequal(rows$seed, seq(101, 120))) {
        stop("Cell-subsetting replicate/seed policy failed for ", ds, "/", target)
      }
      target_counts <- as.numeric(target)
      for (row_name in rownames(rows)) {
        row_i <- which(rownames(results) == row_name)
        if (!identical(
          original_count_lists[[row_i]],
          baseline_original
        ) || !identical(
          count_lists[[row_i]],
          pmin(baseline_counts, target_counts)
        )) {
          stop("Cell-subsetting effective counts violate min(original,target).")
        }
      }
    }
  }
  csv <- tryCatch(
    utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Could not read cell-subsetting CSV: ",
                             conditionMessage(e))
  )
  csv_required <- c(
    "dataset", "target", "target_cells", "replicate", "seed", "sample_ids",
    "original_cells_per_sample", "effective_cells_per_sample",
    "total_cells", "ANOSIM"
  )
  if (!identical(names(csv), names(rds_table)) ||
      !all(csv_required %in% names(csv)) ||
      nrow(csv) != nrow(rds_table)) {
    stop("Cell-subsetting CSV does not match the RDS result table.")
  }
  exact_columns <- setdiff(names(rds_table), c("seed", "ANOSIM", if (portable) evaluation_metrics()))
  for (column in exact_columns) {
    if (!identical(as.character(csv[[column]]),
                   as.character(rds_table[[column]]))) {
      stop("Cell-subsetting CSV column does not match the RDS result table: ",
           column)
    }
  }
  csv_seed_text <- as.character(csv$seed)
  csv_seed <- suppressWarnings(as.numeric(csv_seed_text))
  rds_seed <- results$seed
  allowed_blank_seed <- is.na(csv_seed_text) |
    !nzchar(trimws(csv_seed_text)) |
    trimws(csv_seed_text) == "NA"
  if (any(is.na(csv_seed) & !allowed_blank_seed) ||
      any(!is.na(csv_seed) & !is.finite(csv_seed)) ||
      any(is.na(csv_seed) != is.na(rds_seed)) ||
      any(!is.na(csv_seed) & csv_seed != rds_seed)) {
    stop("Cell-subsetting CSV seed diagnostics do not match the RDS result table.")
  }
  csv_anosim <- suppressWarnings(as.numeric(csv$ANOSIM))
  if (any(!is.finite(csv_anosim))) {
    stop("Cell-subsetting CSV ANOSIM diagnostics are invalid.")
  }
  anosim_tolerance <- 1e-8
  if (any(abs(csv_anosim - results$anosim) >
          anosim_tolerance * pmax(1, abs(results$anosim)))) {
    stop("Cell-subsetting CSV ANOSIM diagnostics do not match the RDS result table.")
  }
  if (portable) {
    for (metric in evaluation_metrics()) {
      if (!is.numeric(csv[[metric]]) || any(!is.finite(csv[[metric]])) ||
          any(abs(csv[[metric]] - results[[metric]]) > 1e-8 * pmax(1, abs(results[[metric]]))))
        stop("Cell-depth CSV metric differs from the checked RDS: ", metric)
    }
    .analysis_checked_file(csv_path, raw$csv$md5)
    .analysis_checked_file(rds_path, entry$source$md5)
  }
  list(results = results, raw = raw, context = context)
}

# Optional derived-run roots
.derived_requested_root <- function(env_name, label) {
  value <- Sys.getenv(env_name, unset = NA_character_)
  if (is.na(value)) return(NULL)
  if (!nzchar(value)) {
    stop(env_name, " must be a non-empty directory when set.")
  }
  path <- normalizePath(path.expand(value), mustWork = FALSE)
  if (!dir.exists(path)) {
    stop("Completed ", label, " derived-run directory not found: ", path)
  }
  path
}

# Annotation eligibility
# datasets.json `not_suitable_for_auto_annotation` flag (see AGENTS.md):
# datasets flagged for the annotation method families are dropped from the
# HiTME/scATOMIC annotation-method figures (Figure 3 A + Supp fig 18) with
# this documented message. The dataset stays in Figure 2a, where the author
# CT columns drive the methods. Absent/empty flag = suitable for all.
flag_status_for <- function(ds_key) {
  e <- datasets[[ds_key]]
  if (is.null(e) || is.null(e$not_suitable_for_auto_annotation)) {
    character(0)
  } else {
    e$not_suitable_for_auto_annotation
  }
}

flag_drop_annotation_methods <- function(df) {
  df %>%
    filter(!(dataset %in% hitme_flagged_ds &
               method %in% c("ECODA_HiTME_HR_layer2", "ECODA_HiTME_HR_layer3"))) %>%
    filter(!(dataset %in% scatomic_flagged_ds & method == "ECODA_scATOMIC_HR"))
}

# Annotation cell-type counts
annotation_cell_type_count <- function(bundle, method, dataset) {
  if (!is.null(bundle$cell_type_names)) {
    names <- trimws(as.character(bundle$cell_type_names))
    if (!length(names) || anyNA(names) || any(!nzchar(names)) || anyDuplicated(names) ||
        any(names %in% c("NA", "nan", "None", "Unknown"))) stop("Invalid projected cell-type identities")
    return(length(names))
  }
  if (is.null(bundle) || !"feat_mat" %in% names(bundle)) {
    stop(
      "Missing feature matrix for annotation method '", method,
      "' in dataset '", dataset, "'."
    )
  }

  if (identical(method, "ECODA_deconv")) {
    # ECODA_deconv has inferred EPIC components but no counts field.
    cell_type_names <- colnames(bundle$feat_mat)
    count_source <- "feat_mat component columns"
  } else {
    if (!"counts" %in% names(bundle) || is.null(bundle$counts)) {
      stop(
        "Missing top-level counts field for annotation method '", method,
        "' in dataset '", dataset, "'."
      )
    }
    cell_type_names <- colnames(bundle$counts)
    count_source <- "counts columns"
  }

  if (is.null(cell_type_names) || length(cell_type_names) == 0L) {
    stop(
      "Empty ", count_source, " for annotation method '", method,
      "' in dataset '", dataset, "'."
    )
  }

  cell_type_names <- trimws(as.character(cell_type_names))
  invalid <- is.na(cell_type_names) |
    !nzchar(cell_type_names) |
    cell_type_names %in% c("NA", "nan", "None", "Unknown")
  valid_names <- cell_type_names[!invalid]
  unique_names <- unique(valid_names)
  if (length(valid_names) == 0L || length(unique_names) != length(valid_names)) {
    stop(
      "Invalid or duplicated cell-type identities in ", count_source,
      " for annotation method '", method, "' in dataset '", dataset, "'."
    )
  }
  length(unique_names)
}

# Execution-log keys
# Normalize keys so the merge with df_benchmark (paste0(dataset, method)
# join, helpers.R::merge_exec_times) hits:
# - GongSharma: internal key Gongsharma_cmv_young_males is displayed as
#   GongSharma (same rename as bmark_result_list above);
# - PILOT/QOT/PILOT-GM-VAE: the feather rows are <method>_hvg<n>_highres,
#   the bundles <method>_hvg<n>.
normalize_exec_keys <- function(df) {
  df %>%
    mutate(
      dataset = recode(dataset, "Gongsharma_cmv_young_males" = "GongSharma"),
      method = gsub("^(PILOT_hvg[0-9]+)_highres$", "\\1", method),
      method = gsub("^(QOT_hvg[0-9]+)_highres$", "\\1", method),
      method = gsub("^(PILOT-GM-VAE_hvg[0-9]+)_highres$", "\\1", method)
    )
}

# Feather rows win for ordinary dataset/method records.  Shared preparation
# rows are different when their explicit timing identities differ; an
# anti_join on only (dataset, method) would silently discard those contexts.
execution_union_key <- function(df) {
  key <- paste(
    as.character(df[["dataset"]]),
    as.character(df[["method"]]),
    sep = "\r"
  )
  shared <- .execution_is_shared_timing_method(
    as.character(df[["method"]])
  )
  if ("timing_id" %in% names(df)) {
    identified <- shared & .execution_nonblank(df[["timing_id"]])
    if (any(identified)) {
      key[identified] <- paste(
        key[identified],
        as.character(df[["timing_id"]][identified]),
        sep = "\r"
      )
    }
  }
  key
}

# A reviewed timing log may predate memory logging. Preserve its runtime but
# recover absent memory from the exact selected dataset/method bundle.
benchmark_complete_execution_memory <- function(data, bundle_rows) {
  bundle_rows <- deduplicate_exec_times(bundle_rows)
  index <- match(execution_union_key(data), execution_union_key(bundle_rows))
  memory <- bundle_rows$mem_GB[index]
  fill <- is.na(data$mem_GB) & is.finite(memory) & memory > 0
  data$memory_source <- data$timing_source
  data$mem_GB[fill] <- memory[fill]
  data$memory_source[fill] <- "selected result bundle"
  data
}

# Per-method display costs include the one shared preparation needed to run
# that method alone. This copy must never replace global shared-once totals.
benchmark_method_execution_costs <- function(data, bundle_rows) {
  data <- data[!.execution_is_shared_timing_method(data$method), , drop = FALSE]
  bundle_rows <- deduplicate_exec_times(bundle_rows)
  index <- match(execution_union_key(data), execution_union_key(bundle_rows))
  schema2 <- .execution_schema2_mask(bundle_rows)[index]
  take <- grepl("^Pseudobulk", data$method) & !is.na(schema2) & schema2
  data$cost_accounting <- "recorded method cost"
  if (any(take)) {
    selected <- bundle_rows[index[take], , drop = FALSE]
    data$time_secs[take] <- selected$shared_time_secs + selected$variant_time_secs
    data$cost_accounting[take] <- "standalone shared + variant"
    # Peak RSS for a standalone method includes its required preparation.
    shared_peak <- selected$shared_mem_GB
    previous_peak <- data$mem_GB[take]
    use_shared <- is.finite(shared_peak) & (is.na(previous_peak) | shared_peak > previous_peak)
    positions <- which(take)[use_shared]
    data$mem_GB[positions] <- shared_peak[use_shared]
    data$memory_source[positions] <- "selected shared preparation bundle"
    data$timing_source[take] <- "selected shared + variant bundle"
  }
  # Schema 2 explicitly means variant-local time_secs to generic timing
  # helpers. The display copy has a different, explicitly labeled meaning.
  data$timing_schema <- NULL
  data
}

# Execution-time and memory plots
scientific_A_10n <- function(x) {
  s <- scales::scientific_format(digits = 2)(x) 
  s_converted <- gsub("e\\+", " %*% 10^", s) 
  return(parse(text = s_converted)) 
}
# Fit and draw the regression in log10(number of cells) / log10(response)
# coordinates. The raw axes are log-scaled below, so the displayed line is
# straight in the same log-log space used for the model.
fit_log_log_lines <- function(data, y_col) {
  groups <- split(data, data$method)
  dplyr::bind_rows(lapply(groups, function(df) {
    x <- as.numeric(df$n_cells)
    y <- as.numeric(df[[y_col]])
    keep <- is.finite(x) & x > 0 & is.finite(y) & y > 0
    x <- x[keep]
    y <- y[keep]
    if (length(x) < 2L) return(NULL)
    fit <- stats::lm(log10(y) ~ log10(x))
    x_grid <- 10^seq(
      log10(min(x)),
      log10(max(x)),
      length.out = 100L
    )
    data.frame(
      method = as.character(df$method[[1L]]),
      n_cells = x_grid,
      fitted = 10^stats::predict(
        fit,
        newdata = data.frame(x = x_grid)
      ),
      stringsAsFactors = FALSE
    )
  }))
}
# Place labels in a deterministic right-hand lane. Repulsion in raw
# coordinates is unstable after log transformation; spacing labels in the
# already-log10 y coordinate guarantees no text overlap.
spread_log_label_y <- function(data, min_gap = 0.2, x_offset = 0.18) {
  data <- data[order(data$log_pred_y, decreasing = TRUE), , drop = FALSE]
  data$log_label_y <- data$log_pred_y
  if (nrow(data) > 1L) {
    for (i in 2:nrow(data)) {
      data$log_label_y[i] <- min(
        data$log_label_y[i],
        data$log_label_y[i - 1L] - min_gap
      )
    }
  }
  data$log_label_x <- max(data$log_max_x) + x_offset
  data
}
spread_log_label_y_centered <- function(
  data,
  min_gap = 0.15,
  x_offset = 0.18
) {
  data <- data[order(data$log_pred_y, decreasing = TRUE), , drop = FALSE]
  n_labels <- nrow(data)
  if (n_labels == 0L) return(data)
  center_y <- mean(range(data$log_pred_y))
  span_y <- min_gap * max(n_labels - 1L, 0L)
  data$log_label_y <- seq(
    center_y + span_y / 2,
    center_y - span_y / 2,
    length.out = n_labels
  )
  data$log_label_x <- max(data$log_max_x) + x_offset
  data
}

# Rank scaling
scale_ranks_global <- function(x, n_meth) {
  if(all(is.na(x))) return(x)
  # BEST rank (1) should map to 1.0
  # WORST rank (n_meth) should map to 0.0
  # Formula: (Worst - Current) / (Worst - Best)
  (n_meth - x) / (n_meth - 1)
}

# Means always name their compared methods; raw metric ranges are never mixed.
get_global_avg <- function(score_frame, method, compared_methods) {
  frame <- prepare_evaluation_plot_data(score_frame, compared_methods)
  if (!method %in% compared_methods) stop("Mean requires a selected method")
  mean(frame$scaled_value[frame$method == method])
}

benchmark_normalized_mean <- function(score_frame, dataset, method, normalization_methods) {
  if (identical(dataset, "Gongsharma_cmv_young_males") && !dataset %in% score_frame$dataset)
    dataset <- "GongSharma"
  frame <- score_frame[score_frame$dataset == dataset, , drop = FALSE]
  if (!nrow(frame)) stop("No reviewed score rows for MDS dataset: ", dataset)
  get_global_avg(frame, method, normalization_methods)
}

# Metric axes
# --- 2. REUSABLE AXIS AND PLOTTING HELPERS ---
# Use each variable's observed range rather than one broad range for every
# metric. This keeps narrow metrics such as Silhouette from looking flat while
# retaining a 5% margin around every observed point.
metric_axis_range <- function(x, metric_name) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    stop("Supp fig 1: no finite values available for ", metric_name)
  }

  observed <- range(x)
  span <- diff(observed)
  padding <- if (span > 0) 0.05 * span else max(0.05, 0.05 * abs(observed[1]))
  observed + c(-padding, padding)
}

metric_axis_breaks <- function(x, limits) {
  observed <- x[is.finite(x)]
  if (length(observed) == 0) {
    stop("Supp fig 1: no finite values available for tick breaks")
  }

  breaks <- scales::breaks_extended(n = 4)(limits)

  # Keep a negative boundary tick visible when the observed data cross zero
  # but the extended-break algorithm starts at zero.
  if (min(observed) < 0 && min(breaks) >= 0) {
    negative_break <- floor(min(observed) / 0.1) * 0.1
    breaks <- sort(unique(c(negative_break, breaks)))
  }

  if (length(breaks) < 3L || length(breaks) > 5L) {
    stop("Supp fig 1: expected 3-5 meaningful ticks, got ", length(breaks))
  }
  breaks
}

build_metric_axis_spec <- function(data, variables) {
  metric_values <- setNames(
    lapply(variables, function(metric) data[[metric]]),
    variables
  )
  axis_ranges <- setNames(
    lapply(variables, function(metric) {
      metric_axis_range(metric_values[[metric]], metric)
    }),
    variables
  )
  axis_breaks <- Map(metric_axis_breaks, metric_values, axis_ranges)
  axis_ranges <- Map(
    function(lims, breaks) range(c(lims, breaks)),
    axis_ranges,
    axis_breaks
  )
  stopifnot(
    all(vapply(
      axis_breaks,
      function(breaks) length(breaks) >= 3L && length(breaks) <= 5L,
      logical(1)
    )),
    all(vapply(names(axis_ranges), function(metric) {
      values <- metric_values[[metric]]
      values <- values[is.finite(values)]
      all(
        values >= axis_ranges[[metric]][1] &
          values <= axis_ranges[[metric]][2]
      )
    }, logical(1)))
  )
  list(
    values = metric_values,
    ranges = axis_ranges,
    breaks = axis_breaks
  )
}

format_metric_ticks <- function(x) {
  x[abs(x) < 0.005] <- 0
  sprintf("%.2f", x)
}

# Pairwise metric plots
# Add variable-specific meaningful 3-5 tick scales to each GGally panel. The
# coordinate limits are applied with coord_cartesian() so smoothing/statistical
# layers still see all rows; limits include every finite observation.
add_pairwise_axes <- function(plot, mapping, axis_spec, include_y = TRUE) {
  x_name <- rlang::as_name(mapping$x)
  plot <- plot +
    scale_x_continuous(
      breaks = axis_spec$breaks[[x_name]],
      labels = format_metric_ticks,
      expand = c(0, 0)
    )

  if (include_y) {
    y_name <- rlang::as_name(mapping$y)
    plot <- plot +
      scale_y_continuous(
        breaks = axis_spec$breaks[[y_name]],
        labels = format_metric_ticks,
        expand = c(0, 0)
      ) +
      coord_cartesian(
        xlim = axis_spec$ranges[[x_name]],
        ylim = axis_spec$ranges[[y_name]],
        expand = FALSE
      )
  } else {
    plot <- plot +
      coord_cartesian(
        xlim = axis_spec$ranges[[x_name]],
        expand = FALSE
      )
  }

  plot
}

make_pairwise_upper <- function(digits = 2, text_size = 6.5) {
  function(data, mapping, ...) {
    ggally_cor(
      data,
      mapping,
      ...,
      stars = FALSE,
      digits = digits,
      title_args = list(size = text_size),
      group_args = list(size = text_size)
    )
  }
}

make_pairwise_lower <- function(axis_spec) {
  function(data, mapping, ...) {
    plot <- ggplot(data, mapping) +
      geom_point(alpha = 0.4, color = "#000000") +
      geom_smooth(method = "lm", color = "black", se = TRUE)
    add_pairwise_axes(plot, mapping, axis_spec)
  }
}

make_pairwise_diag <- function(axis_spec) {
  function(data, mapping, ...) {
    add_pairwise_axes(
      ggally_densityDiag(data, mapping, ...),
      mapping,
      axis_spec,
      include_y = FALSE
    )
  }
}


# Metric-versus-mean plots
# --- 3. DYNAMICALLY GENERATE "METRIC VS MEAN" PLOTS ---
# purrr::map iterates over your metrics and creates a list of plots
plot_mean_corr <- function(data, metric_var, axis_spec) {
  ggplot(data, aes(x = Mean_Score, y = .data[[metric_var]])) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.8
    ) +
    geom_point(alpha = 0.4, color = "#000000", size = 1.5) +
    geom_smooth(method = "lm", color = "black", linewidth = 1, se = FALSE) +
    stat_cor(
      aes(label = after_stat(rr.label)),
      label.x.npc = "right",
      label.y.npc = "bottom",
      hjust = 1,
      vjust = 0,
      size = 5
    ) +
    labs(x = NULL, y = metric_var) +
    scale_x_continuous(
      breaks = axis_spec$breaks[["Mean_Score"]],
      labels = format_metric_ticks,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = axis_spec$breaks[[metric_var]],
      labels = format_metric_ticks,
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = axis_spec$ranges[["Mean_Score"]],
      ylim = axis_spec$ranges[[metric_var]],
      expand = FALSE
    ) +
    theme_classic(base_size = 15) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 13),
      aspect.ratio = 1
    )
}


# Bassez compact-bundle diagnostics. E/NE denote Expander/Non-expander, not
# clinical response. All associations are descriptive, not causal/inferential.
# PCA uses no labels; variance is not group association, signs are arbitrary,
# and VST feature values do not establish low expression or expression absence.
bassez_feature_diagnostics <- function(method_results) {
  main <- c("ECODA_authors_HR", "Pseudobulk_hvg2000")
  if (!is.list(method_results) || is.null(names(method_results)) ||
      anyNA(names(method_results)) || any(!nzchar(names(method_results))) ||
      anyDuplicated(names(method_results))) {
    stop("Bassez diagnostics require uniquely named method bundles.")
  }
  if (!all(main %in% names(method_results)) ||
      any(vapply(method_results[main], is.null, logical(1)))) {
    stop("Bassez diagnostics require ECODA_authors_HR and Pseudobulk_hvg2000.")
  }
  # Deliberately exclude the shuffled-label NULL control from aligned summaries.
  selected <- c(
    main, "ECODA_authors_LR", "Freq_highres",
    "ECODA_HiTME_HR_layer2", "ECODA_HiTME_HR_layer3", "ECODA_scATOMIC_HR",
    paste0("ECODA_seuratres_", c("0.1", "0.4", "2", "5", "20")),
    paste0("Pseudobulk_hvg", c(500, 1000, 3000)), "Pseudobulk_schvg2000",
    "MOFA_hvg2000_factors15", "scITD_hvg2000_factors5",
    "MrVI_hvg2000", "scPoli_hvg2000_dims15_highres",
    "PILOT_hvg2000", "PILOT_hvg2000_HR", "PILOT-GM-VAE_hvg2000",
    "QOT_hvg2000", "GloScope_hvg2000_pcadims30", "GloProp",
    "Pseudobulk_CT_HR_hvg2000"
  )
  selected <- intersect(selected, names(method_results))
  selected <- selected[!vapply(method_results[selected], is.null, logical(1))]
  valid_names <- function(x) {
    !is.null(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  }
  finite_matrix <- function(x, context, require_columns = TRUE) {
    x <- as.matrix(x)
    # Some compact artifacts attach names to their dimname vectors; these
    # attributes are not sample IDs and must not affect exact ID comparison.
    dimnames(x) <- lapply(dimnames(x), unname)
    if (!is.numeric(x) || length(dim(x)) != 2L ||
        nrow(x) < 1L || ncol(x) < 1L || any(!is.finite(x)) ||
        !valid_names(rownames(x)) ||
        (require_columns && !valid_names(colnames(x)))) {
      stop(context, " must be finite numeric values with unique sample names",
           if (require_columns) " and unique column names." else ".")
    }
    x
  }
  matrices <- labels <- distances <- list()
  reference_samples <- reference_labels <- NULL
  method_rows <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    method <- selected[[i]]
    bundle <- method_results[[method]]
    if (!is.list(bundle) ||
        !all(c("feat_mat", "labels", "dist_mat", "scores") %in% names(bundle))) {
      stop("Incomplete Bassez bundle: ", method)
    }
    x <- finite_matrix(
      bundle$feat_mat, paste(method, "feat_mat"),
      require_columns = method %in% c(main, "Pseudobulk_hvg3000")
    )
    y <- bundle$labels
    if (!is.null(names(y))) names(y) <- unname(names(y))
    if (!is.atomic(y) || !valid_names(names(y)) ||
        !identical(names(y), rownames(x)) ||
        anyNA(y) || !setequal(as.character(y), c("E", "NE"))) {
      stop(method, ": named labels must match feature rows exactly and contain E/NE.")
    }
    y <- setNames(as.character(y), names(y))
    if (is.null(reference_samples)) {
      reference_samples <- rownames(x)
      reference_labels <- y
    } else if (identical(method, "scITD_hvg2000_factors5")) {
      if (!all(rownames(x) %in% reference_samples) ||
          !identical(y, reference_labels[rownames(x)])) {
        stop(method, ": retained sample IDs or E/NE labels do not match reference.")
      }
    } else if (!identical(rownames(x), reference_samples) ||
               !identical(y, reference_labels)) {
      stop(method, ": sample order or E/NE labels differ across Bassez methods.")
    }
    d <- finite_matrix(bundle$dist_mat, paste(method, "dist_mat"))
    if (!identical(rownames(d), rownames(x)) ||
        !identical(colnames(d), rownames(x)) ||
        any(d < 0) || any(diag(d) != 0) ||
        !isTRUE(all.equal(d, t(d), tolerance = 1e-10))) {
      stop(method, ": invalid or sample-misaligned distance matrix.")
    }
    if (method %in% main && nrow(x) == ncol(x) &&
        setequal(colnames(x), rownames(x))) {
      stop(method, ": distance-as-features matrix is not a gene/composition matrix.")
    }
    scores <- bundle$scores
    stored_score <- function(key) {
      if (!key %in% names(scores)) return(NA_real_)
      value <- scores[[key]]
      if (is.null(value)) return(NA_real_)
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop(method, ": invalid stored score ", key)
      }
      as.numeric(value)
    }
    if ((!is.list(scores) && !is.numeric(scores)) || is.null(names(scores))) {
      stop(method, ": scores must be named stored values.")
    }
    method_rows[[i]] <- data.frame(
      method = method, n_samples = nrow(x), n_features = ncol(x),
      n_E = sum(y == "E"), n_NE = sum(y == "NE"),
      n_samples_overlap = sum(rownames(x) %in% reference_samples),
      n_samples_missing = length(setdiff(reference_samples, rownames(x))),
      sample_scope = if (nrow(x) < length(reference_samples)) {
        "Retained scITD subset; stored scores describe this subset"
      } else {
        "Complete reference sample set"
      },
      ANOSIM = stored_score("anosim_score"),
      mod_knn3 = stored_score("mod_knn3_score"),
      ARI = stored_score("cluster_score"),
      stringsAsFactors = FALSE
    )
    matrices[[method]] <- x
    labels[[method]] <- y
    distances[[method]] <- d
  }
  if (length(reference_samples) < 3L) {
    stop("Bassez diagnostics require at least three samples.")
  }
  group_r2 <- function(x, y) {
    total <- sum((x - mean(x))^2)
    if (total == 0) return(NA_real_)
    fitted <- ifelse(y == "E", mean(x[y == "E"]), mean(x[y == "NE"]))
    sum((fitted - mean(x))^2) / total
  }
  feature_rows <- pc_rows <- fits <- list()
  for (method in main) {
    x <- matrices[[method]]
    y <- labels[[method]]
    variance <- apply(x, 2L, stats::var)
    delta <- colMeans(x[y == "E", , drop = FALSE]) -
      colMeans(x[y == "NE", , drop = FALSE])
    fit <- if (any(variance > 0)) {
      stats::prcomp(x, center = TRUE, scale. = FALSE)
    } else {
      NULL
    }
    nonzero <- if (is.null(fit)) integer() else {
      which(fit$sdev > max(fit$sdev) * sqrt(.Machine$double.eps))
    }
    pc_percent <- if (is.null(fit)) numeric() else {
      100 * fit$sdev^2 / sum(fit$sdev^2)
    }
    loading <- contribution <- matrix(
      NA_real_, nrow = ncol(x), ncol = 2L,
      dimnames = list(colnames(x), c("PC1", "PC2"))
    )
    for (pc in intersect(1:2, nonzero)) {
      loading[, pc] <- fit$rotation[, pc]
      contribution[, pc] <- 100 * fit$rotation[, pc]^2
    }
    first_two <- intersect(1:2, nonzero)
    weighted <- if (length(first_two)) {
      as.vector(contribution[, first_two, drop = FALSE] %*%
                  (pc_percent[first_two] / sum(pc_percent[first_two])))
    } else {
      rep(NA_real_, ncol(x))
    }
    feature_rows[[method]] <- data.frame(
      method = method, feature = colnames(x),
      variance = variance, variance_rank = rank(-variance, ties.method = "min"),
      mean_E = colMeans(x[y == "E", , drop = FALSE]),
      mean_NE = colMeans(x[y == "NE", , drop = FALSE]),
      delta_E_minus_NE = delta,
      group_R2 = apply(x, 2L, group_r2, y = y),
      loading_PC1 = loading[, 1L], loading_PC2 = loading[, 2L],
      contribution_PC1_percent = contribution[, 1L],
      contribution_PC2_percent = contribution[, 2L],
      contribution_PC1_PC2_percent = weighted,
      group_contrast_percent = if (sum(delta^2) > 0) {
        100 * delta^2 / sum(delta^2)
      } else {
        rep(NA_real_, ncol(x))
      },
      stringsAsFactors = FALSE, row.names = NULL
    )
    pc_rows[[method]] <- data.frame(
      method = rep(method, length(nonzero)),
      PC = paste0("PC", nonzero)[seq_along(nonzero)],
      variance_percent = pc_percent[nonzero],
      group_R2 = vapply(nonzero, function(pc) group_r2(fit$x[, pc], y), numeric(1)),
      stringsAsFactors = FALSE
    )
    fits[method] <- list(fit)
  }
  features <- do.call(rbind, feature_rows)
  rownames(features) <- NULL
  hr <- main[[1L]]
  counts <- finite_matrix(method_results[[hr]]$counts, "ECODA HR counts")
  if (!identical(rownames(counts), reference_samples) ||
      !all(colnames(matrices[[hr]]) %in% colnames(counts))) {
    stop("ECODA HR counts must align exactly with samples and cover retained features.")
  }
  # process_coda_fig() stores counts_all +0.5 in every entry. Match the
  # recovery assertion in run_local_ecoda_derived.R before claiming raw counts.
  counts <- counts[, colnames(matrices[[hr]]), drop = FALSE] - 0.5
  if (any(counts < 0) || any(counts != floor(counts))) {
    stop("ECODA HR counts are not recoverable nonnegative integer+0.5 values.")
  }
  totals <- rowSums(counts)
  if (any(totals <= 0)) stop("ECODA HR retained categories have an empty sample.")
  proportions <- 100 * counts / totals
  y <- labels[[hr]]
  med_E <- apply(proportions[y == "E", , drop = FALSE], 2L, stats::median)
  med_NE <- apply(proportions[y == "NE", , drop = FALSE], 2L, stats::median)
  proportions <- data.frame(
    method = hr, feature = colnames(counts),
    median_raw_percent_E = med_E, median_raw_percent_NE = med_NE,
    delta_median_raw_percent_E_minus_NE = med_E - med_NE,
    prevalence_E = colMeans(counts[y == "E", , drop = FALSE] > 0),
    prevalence_NE = colMeans(counts[y == "NE", , drop = FALSE] > 0),
    zero_fraction_E = colMeans(counts[y == "E", , drop = FALSE] == 0),
    zero_fraction_NE = colMeans(counts[y == "NE", , drop = FALSE] == 0),
    denominator = paste(
      "Recovered raw cells in retained author-subtyped HR feature categories",
      "(T/myeloid in this artifact), not all immune or tumor cells"
    ),
    n_denominator_categories = ncol(counts), pseudocount_removed = 0.5,
    denominator_min_E = min(totals[y == "E"]),
    denominator_max_E = max(totals[y == "E"]),
    denominator_min_NE = min(totals[y == "NE"]),
    denominator_max_NE = max(totals[y == "NE"]),
    stringsAsFactors = FALSE, row.names = NULL
  )
  panel <- c(
    "CXCL13", "PDCD1", "CTLA4", "LAG3", "TIGIT", "TOX", "IFNG", "GZMB",
    "TCF7", "CCR7", "CD3D", "CD8A", "CD4", "FOXP3", "MS4A1", "CD79A",
    "MZB1", "JCHAIN"
  )
  marker_rows <- lapply(
    intersect(c("Pseudobulk_hvg2000", "Pseudobulk_hvg3000"), selected),
    function(method) {
      variance <- apply(matrices[[method]], 2L, stats::var)
      index <- match(panel, names(variance))
      data.frame(
        method = method, marker = panel, present = !is.na(index),
        variance = unname(variance[index]),
        variance_rank = unname(rank(-variance, ties.method = "min")[index]),
        interpretation = "Coverage of selected VST features only; absence is not full-expression exclusion",
        stringsAsFactors = FALSE
      )
    }
  )
  pb <- main[[2L]]
  x <- matrices[[pb]]
  y <- labels[[pb]]
  ig <- grepl("^IG[HKL][VDJC]", colnames(x))
  variable <- apply(x, 2L, stats::var) > 0
  vegan_available <- requireNamespace("vegan", quietly = TRUE)
  anosim_r <- function(x, y, distance = FALSE) {
    if (!vegan_available || (!distance && ncol(x) == 0L) ||
        length(unique(y)) < 2L || length(y) < 3L) return(NA_real_)
    d <- if (distance) stats::as.dist(x) else stats::dist(x)
    if (all(as.vector(d) == 0)) return(NA_real_)
    as.numeric(vegan::anosim(d, grouping = y, permutations = 0)$statistic)
  }
  no_ig <- anosim_r(x[, !ig, drop = FALSE], y)
  unit <- if (any(variable)) {
    anosim_r(scale(x[, variable, drop = FALSE], center = TRUE, scale = TRUE), y)
  } else {
    NA_real_
  }
  loo <- vapply(seq_along(reference_samples), function(i) {
    anosim_r(distances[[hr]][-i, -i, drop = FALSE], labels[[hr]][-i], TRUE)
  }, numeric(1))
  valid_loo <- loo[is.finite(loo)]
  pb_features <- feature_rows[[pb]]
  sensitivity <- data.frame(
    method = c(pb, pb, hr),
    comparison = c("Exclude immunoglobulin genes ^IG[HKL][VDJC]",
                   "Unit-variance features; zero-variance features dropped",
                   "Leave one sample out; stored Euclidean HR distances"),
    ANOSIM_R = c(no_ig, unit, NA_real_),
    ANOSIM_R_min = c(NA_real_, NA_real_,
                     if (length(valid_loo)) min(valid_loo) else NA_real_),
    ANOSIM_R_max = c(NA_real_, NA_real_,
                     if (length(valid_loo)) max(valid_loo) else NA_real_),
    n_samples = c(nrow(x), nrow(x), nrow(x) - 1L),
    n_features = c(sum(!ig), sum(variable), ncol(matrices[[hr]])),
    n_features_removed = c(sum(ig), sum(!variable), 0L),
    n_evaluable = c(as.integer(is.finite(no_ig)), as.integer(is.finite(unit)),
                    length(valid_loo)),
    status = if (!vegan_available) {
      rep("Skipped: optional vegan namespace unavailable", 3L)
    } else {
      ifelse(c(is.finite(no_ig), is.finite(unit), length(valid_loo) > 0),
             "Exploratory descriptive statistic; permutations=0", "Not estimable")
    },
    stringsAsFactors = FALSE
  )
  # Shares describe selected PB2000 features, not immune-cell causality/noise.
  sensitivity$IG_variance_percent <- c(
    if (sum(pb_features$variance) > 0) {
      100 * sum(pb_features$variance[ig]) / sum(pb_features$variance)
    } else NA_real_, NA_real_, NA_real_
  )
  sensitivity$IG_PC1_contribution_percent <- c(
    if (all(is.finite(pb_features$contribution_PC1_percent))) {
      sum(pb_features$contribution_PC1_percent[ig])
    } else NA_real_, NA_real_, NA_real_
  )
  sensitivity$IG_PC2_contribution_percent <- c(
    if (all(is.finite(pb_features$contribution_PC2_percent))) {
      sum(pb_features$contribution_PC2_percent[ig])
    } else NA_real_, NA_real_, NA_real_
  )
  hr_labels <- labels[[hr]]
  keep_non_EX <- !grepl("_EX", colnames(matrices[[hr]]))
  without_EX <- anosim_r(
    matrices[[hr]][, keep_non_EX, drop = FALSE], hr_labels
  )
  # Exact multivariate-hypergeometric sampling without replacement, conditional
  # on the retained HR categories only. Neither sampling nor CLR uses labels.
  depth <- min(totals)
  rarefied <- rep(NA_real_, 100L)
  if (vegan_available) {
    old_rng_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      do.call(RNGkind, as.list(old_rng_kind))
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(20261023, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    for (draw in seq_along(rarefied)) {
      sampled <- matrix(0, nrow(counts), ncol(counts), dimnames = dimnames(counts))
      for (i in seq_len(nrow(counts))) {
        remaining_draw <- depth
        remaining_total <- totals[[i]]
        for (j in seq_len(ncol(counts) - 1L)) {
          value <- if (remaining_draw == 0 || counts[i, j] == 0) 0 else {
            stats::rhyper(1L, counts[i, j],
                          remaining_total - counts[i, j], remaining_draw)
          }
          sampled[i, j] <- value
          remaining_draw <- remaining_draw - value
          remaining_total <- remaining_total - counts[i, j]
        }
        sampled[i, ncol(counts)] <- remaining_draw
      }
      log_counts <- log(sampled + 0.5)
      rarefied[[draw]] <- anosim_r(log_counts - rowMeans(log_counts), hr_labels)
    }
  }
  valid_rarefied <- rarefied[is.finite(rarefied)]
  extra <- sensitivity[rep(3L, 2L), , drop = FALSE]
  extra$comparison <- c(
    "Exclude _EX features from existing CLR feature space; no re-CLR",
    "Retained-HR depth-matched rarefaction; 100 draws; CLR(counts+0.5)"
  )
  extra$ANOSIM_R <- c(
    without_EX, if (length(valid_rarefied)) mean(valid_rarefied) else NA_real_
  )
  extra$ANOSIM_R_min <- c(
    NA_real_, if (length(valid_rarefied)) min(valid_rarefied) else NA_real_
  )
  extra$ANOSIM_R_max <- c(
    NA_real_, if (length(valid_rarefied)) max(valid_rarefied) else NA_real_
  )
  extra$n_samples <- nrow(counts)
  extra$n_features <- c(sum(keep_non_EX), ncol(counts))
  extra$n_features_removed <- c(sum(!keep_non_EX), 0L)
  extra$n_evaluable <- c(as.integer(is.finite(without_EX)), length(valid_rarefied))
  extra$status <- if (!vegan_available) {
    rep("Skipped: optional vegan namespace unavailable", 2L)
  } else {
    ifelse(extra$n_evaluable > 0,
           "Exploratory descriptive statistic; permutations=0", "Not estimable")
  }
  sensitivity <- rbind(sensitivity, extra)
  rownames(sensitivity) <- NULL
  sensitivity$retained_raw_cells_median_E <- c(
    NA_real_, NA_real_, rep(stats::median(totals[hr_labels == "E"]), 3L)
  )
  sensitivity$retained_raw_cells_median_NE <- c(
    NA_real_, NA_real_, rep(stats::median(totals[hr_labels == "NE"]), 3L)
  )
  hr_fit <- fits[[hr]]
  pc1_depth_cor <- if (!is.null(hr_fit) && hr_fit$sdev[[1L]] > 0 &&
                       stats::sd(log(totals)) > 0) {
    stats::cor(hr_fit$x[, 1L], log(totals))
  } else NA_real_
  sensitivity$PC1_log_retained_cells_correlation <- c(
    NA_real_, NA_real_, rep(pc1_depth_cor, 3L)
  )
  sensitivity$rarefaction_depth <- c(rep(NA_real_, 4L), depth)
  sensitivity$seed <- c(rep(NA_integer_, 4L), 20261023L)
  list(
    methods = do.call(rbind, method_rows), features = features,
    pcs = do.call(rbind, pc_rows), proportions = proportions,
    marker_coverage = do.call(rbind, marker_rows), sensitivity = sensitivity,
    fits = fits
  )
}

# Matched all-cell comparison from original metadata, never an HR/LR append:
# subtype-missing T/myeloid cells must also enter the broad fallback categories.
bassez_allcell_comparison <- function(annotation_counts, method_results) {
  required <- c("ECODA_authors_HR", "ECODA_authors_LR", "Pseudobulk_hvg2000")
  if (!is.list(method_results) || is.null(names(method_results)) ||
      anyNA(names(method_results)) || anyDuplicated(names(method_results)) ||
      !all(required %in% names(method_results)) ||
      any(vapply(method_results[required], is.null, logical(1)))) {
    stop("All-cell comparison requires unique HR, LR and PB HVG2000 bundles.")
  }
  valid_names <- function(x) {
    !is.null(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  }
  matrix_values <- function(x, context) {
    x <- as.matrix(x)
    dimnames(x) <- lapply(dimnames(x), unname)
    if (!is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L ||
        any(!is.finite(x)) || !valid_names(rownames(x)) ||
        !valid_names(colnames(x))) {
      stop(context, ": expected finite numeric matrix with unique dimnames.")
    }
    x
  }
  hr <- method_results[["ECODA_authors_HR"]]
  samples <- unname(rownames(hr$feat_mat))
  labels <- hr$labels
  if (!is.atomic(labels) || !valid_names(names(labels)) ||
      !identical(unname(names(labels)), samples) || length(samples) != 26L ||
      anyNA(labels) || sum(labels == "E") != 9L || sum(labels == "NE") != 17L) {
    stop("Expected the original 26 Bassez samples with named 9 E / 17 NE labels.")
  }
  labels <- setNames(as.character(labels), samples)
  selected <- c(required, "ECODA_seuratres_2")
  selected <- selected[selected %in% names(method_results)]
  selected <- selected[!vapply(method_results[selected], is.null, logical(1))]
  features <- distances <- list()
  for (method in selected) {
    bundle <- method_results[[method]]
    x <- matrix_values(bundle$feat_mat, paste(method, "features"))
    y <- setNames(as.character(bundle$labels), unname(names(bundle$labels)))
    if (!identical(rownames(x), samples) || !identical(y, labels)) {
      stop(method, ": sample IDs, order or E/NE labels differ from HR.")
    }
    d <- matrix_values(bundle$dist_mat, paste(method, "distances"))
    if (!identical(rownames(d), samples) || !identical(colnames(d), samples) ||
        any(d < 0) || any(diag(d) != 0) ||
        !isTRUE(all.equal(d, t(d), tolerance = 1e-10))) {
      stop(method, ": invalid or misaligned distance matrix.")
    }
    features[[method]] <- x
    distances[[method]] <- stats::as.dist(d)
  }
  recover_counts <- function(method) {
    raw <- matrix_values(method_results[[method]]$counts,
                         paste(method, "counts")) - 0.5
    if (!identical(rownames(raw), samples) ||
        any(raw < 0) || any(raw != floor(raw))) {
      stop(method, ": counts must align with samples and be integer+0.5.")
    }
    raw
  }
  hr_raw <- recover_counts("ECODA_authors_HR")
  lr_raw <- recover_counts("ECODA_authors_LR")
  if (ncol(hr_raw) != 33L || sum(hr_raw) != 24139 ||
      sum(lr_raw) != 74948 ||
      !all(colnames(features[["ECODA_authors_HR"]]) %in% colnames(hr_raw))) {
    stop("Baseline counts do not match original HR33 (24139) / all-cell (74948) scope.")
  }
  columns <- c("Sample", "cellType", "cellSubType", "count")
  if (!is.data.frame(annotation_counts) ||
      !all(columns %in% names(annotation_counts)) || !nrow(annotation_counts)) {
    stop("annotation_counts requires Sample, cellType, cellSubType, count.")
  }
  tab <- annotation_counts[, columns, drop = FALSE]
  tab[1:3] <- lapply(tab[1:3], as.character)
  if (anyNA(tab$Sample) || any(!nzchar(tab$Sample)) ||
      !setequal(tab$Sample, samples) ||
      anyNA(tab$cellType) ||
      any(tab$cellType %in% c("NA", "nan", "None", "Unknown", "")) ||
      !is.numeric(tab$count) || any(!is.finite(tab$count)) ||
      any(tab$count < 0) || any(tab$count != floor(tab$count)) ||
      anyDuplicated(tab[1:3])) {
    stop("Cross-tab must have unique keys, exact sample IDs, valid broad types and nonnegative integer counts.")
  }
  fine <- !is.na(tab$cellSubType) &
    !tab$cellSubType %in% c("NA", "nan", "None", "Unknown", "")
  if (!setequal(tab$cellSubType[fine], colnames(hr_raw)) ||
      !setequal(tab$cellType, colnames(lr_raw))) {
    stop("Cross-tab fine/broad categories differ from original HR/LR counts.")
  }
  # Original metadata can assign the same fine subtype to different broad
  # parents. Preserve the author fine labels; broad parentage is not a filter.
  # Broad and fine pDC intentionally denote the same population: literal
  # coalescing combines its 241 fine-labelled and 10 fallback cells.
  aliases <- intersect(tab$cellType[!fine], tab$cellSubType[fine])
  if (length(setdiff(aliases, "pDC"))) {
    stop("Unsupported fine/fallback category collision beyond source-verified pDC.")
  }
  count_matrix <- function(rows, categories, levels) {
    out <- matrix(0, length(samples), length(levels),
                  dimnames = list(samples, levels))
    for (i in which(rows)) {
      s <- match(tab$Sample[[i]], samples)
      c <- match(categories[[i]], levels)
      out[s, c] <- out[s, c] + tab$count[[i]]
    }
    out
  }
  reconstructed_hr <- count_matrix(fine, tab$cellSubType, colnames(hr_raw))
  reconstructed_lr <- count_matrix(rep(TRUE, nrow(tab)), tab$cellType,
                                   colnames(lr_raw))
  if (!identical(as.vector(reconstructed_hr), as.vector(hr_raw)) ||
      !identical(as.vector(reconstructed_lr), as.vector(lr_raw))) {
    stop("Original metadata must reconstruct every raw HR and LR sample/category exactly.")
  }
  category <- ifelse(fine, tab$cellSubType, tab$cellType)
  raw_counts <- count_matrix(rep(TRUE, nrow(tab)), category, sort(unique(category)))
  if (sum(raw_counts) != 74948 ||
      !identical(unname(rowSums(raw_counts)), unname(rowSums(lr_raw)))) {
    stop("Mixed-resolution counts do not conserve all 74948 original cells.")
  }
  counts <- raw_counts + 0.5
  # Existing row-wise CLR convention; all features retained without E/NE input.
  feat_mat <- as.matrix(clr(counts))
  dimnames(feat_mat) <- dimnames(counts)
  mixed <- "ECODA_authors_mixed_allcells"
  features[[mixed]] <- feat_mat
  distances[[mixed]] <- stats::dist(feat_mat)
  methods <- c("ECODA_authors_HR", mixed, "Pseudobulk_hvg2000",
               intersect("ECODA_seuratres_2", selected))
  n_cells <- c(ECODA_authors_HR = sum(hr_raw),
               ECODA_authors_mixed_allcells = sum(raw_counts),
               Pseudobulk_hvg2000 = sum(lr_raw))
  if ("ECODA_seuratres_2" %in% methods) {
    leiden_raw <- recover_counts("ECODA_seuratres_2")
    if (!identical(unname(rowSums(leiden_raw)), unname(rowSums(lr_raw)))) {
      stop("Leiden2 counts do not describe the same all-cell sample totals.")
    }
    n_cells["ECODA_seuratres_2"] <- sum(leiden_raw)
  }
  comparison <- do.call(rbind, lapply(methods, function(method) {
    d <- distances[[method]]
    data.frame(
      method = method, n_cells = unname(n_cells[[method]]),
      n_samples = nrow(features[[method]]), n_features = ncol(features[[method]]),
      ANOSIM = unname(vegan::anosim(d, labels, permutations = 0,
                                   parallel = 1)[["statistic"]]),
      mod_knn3 = calc_modularity(d, labels, knn_k = 3),
      ARI = clust_eval(d, labels, return_mean = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  fallback <- count_matrix(!fine, tab$cellType, colnames(lr_raw))
  coverage <- data.frame(
    Sample = samples, expansion = unname(labels),
    n_cells = rowSums(raw_counts), n_fine = rowSums(reconstructed_hr),
    n_fallback = rowSums(fallback),
    fallback_fraction = rowSums(fallback) / rowSums(raw_counts),
    row.names = NULL
  )
  fallback_counts <- data.frame(Sample = samples, fallback, check.names = FALSE,
                                row.names = NULL)
  list(methods = comparison, raw_counts = raw_counts, counts = counts,
       feat_mat = feat_mat, labels = labels, dist = distances[[mixed]],
       coverage = coverage, fallback_counts = fallback_counts)
}

# Stephenson metadata-only cell scores: no H5AD/expression matrix is opened.
load_stephenson_ucell_cell_scores <- function(root) {
  root <- .derived_require_scalar(root, "Stephenson cell-score root")
  manifest <- .derived_read_json(file.path(root, "manifest.json"),
                                 "Stephenson cell-score manifest")
  .derived_require_fields(
    manifest,
    c("schema_version", "dataset", "view", "n_cells", "n_samples",
      "sample_ids", "thresholds", "files", "source"),
    "Stephenson cell-score manifest"
  )
  if (!identical(manifest$schema_version, 1L) &&
      !identical(manifest$schema_version, 1)) {
    stop("Unsupported Stephenson cell-score schema_version.")
  }
  if (!identical(manifest$dataset, "Stephenson") ||
      !identical(manifest$view, "benchmark_analysis") ||
      !is.numeric(manifest$n_cells) || length(manifest$n_cells) != 1L ||
      !isTRUE(manifest$n_cells == 418768) ||
      !is.numeric(manifest$n_samples) || length(manifest$n_samples) != 1L ||
      !isTRUE(manifest$n_samples == 55) ||
      !is.list(manifest$source) || !length(manifest$source)) {
    stop("Stephenson cell-score manifest identity or source provenance is invalid.")
  }
  sample_ids <- manifest$sample_ids
  if (!is.list(sample_ids) || length(sample_ids) != 55L ||
      !all(vapply(sample_ids, function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) &&
          nzchar(trimws(x)) && !grepl("[[:cntrl:]]", x)
      }, logical(1)))) {
    stop("Stephenson manifest requires 55 nonempty sample IDs.")
  }
  sample_ids <- unname(unlist(sample_ids, use.names = FALSE))
  if (anyDuplicated(sample_ids) ||
      !identical(sample_ids, sort(sample_ids, method = "radix"))) {
    stop("Stephenson manifest sample IDs must be unique and sorted.")
  }
  scores <- c("cellCycle.G1S_UCell", "cellCycle.G2M_UCell", "IFN_UCell")
  if (!is.list(manifest$thresholds) ||
      !setequal(names(manifest$thresholds), scores) ||
      anyDuplicated(names(manifest$thresholds))) {
    stop("Stephenson manifest must define exactly the G1/S, G2/M and IFN thresholds.")
  }
  thresholds <- do.call(rbind, lapply(scores, function(score) {
    entry <- manifest$thresholds[[score]]
    .derived_require_fields(entry, c("threshold", "comparison"), score)
    value <- entry$threshold
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 1 || !identical(entry$comparison, ">")) {
      stop(score, ": threshold must be numeric in [0, 1] with comparison '>'.")
    }
    data.frame(score = score, threshold = value, comparison = ">",
               stringsAsFactors = FALSE)
  }))
  .derived_require_fields(manifest$files, "cells.feather", "Stephenson files")
  expected_digest <- .derived_require_scalar(
    manifest$files[["cells.feather"]], "cells.feather SHA256"
  )
  if (!grepl("^[[:xdigit:]]{64}$", expected_digest)) {
    stop("Invalid cells.feather SHA256.")
  }
  path <- .derived_run_artifact_path(root, "cells.feather", "Stephenson cells")
  if (!file.exists(path) || !isTRUE(file.info(path)$isdir == FALSE) ||
      !isTRUE(file.info(path)$size > 0) ||
      !identical(digest::digest(file = path, algo = "sha256"),
                 tolower(expected_digest))) {
    stop("Missing, empty or SHA256-mismatched Stephenson cells.feather.")
  }
  # Select only the bounded metadata contract, never expression columns.
  columns <- c("Sample", "Status", "layer2", "layer3", "initial_clustering",
               "full_clustering", "IFN_UCell", "HeatShock_UCell",
               "cellCycle.G1S_UCell", "cellCycle.G2M_UCell")
  cells <- as.data.frame(arrow::read_feather(path, col_select = tidyselect::all_of(columns)))
  if (!identical(names(cells), columns) || nrow(cells) != 418768L) {
    stop("Stephenson cell-score columns or 418768-cell scope do not match.")
  }
  for (column in columns[1:6]) {
    value <- cells[[column]]
    if (!is.character(value) && !is.factor(value)) {
      stop(column, ": expected character metadata.")
    }
    cells[[column]] <- as.character(value)
    present <- !is.na(cells[[column]])
    if (any(!nzchar(trimws(cells[[column]][present]))) ||
        any(grepl("[[:cntrl:]]", cells[[column]][present]))) {
      stop(column, ": empty or control-character metadata labels are invalid.")
    }
  }
  if (anyNA(cells$Sample) || anyNA(cells$Status) ||
      !identical(sort(unique(cells$Sample), method = "radix"), sample_ids) ||
      !all(cells$Status %in% c("Covid", "Healthy"))) {
    stop("Stephenson cells require the manifest sample IDs and Covid/Healthy labels.")
  }
  sample_status <- unique(cells[c("Sample", "Status")])
  if (nrow(sample_status) != 55L || anyDuplicated(sample_status$Sample) ||
      sum(sample_status$Status == "Covid") != 43L ||
      sum(sample_status$Status == "Healthy") != 12L) {
    stop("Expected constant per-sample labels: 43 Covid and 12 Healthy.")
  }
  for (score in columns[7:10]) {
    value <- cells[[score]]
    if (!is.numeric(value) || any(is.nan(value)) ||
        any(!is.na(value) & (!is.finite(value) | value < 0 | value > 1))) {
      stop(score, ": expected numeric scores in [0, 1] or NA, without imputation.")
    }
  }
  list(cells = cells, thresholds = thresholds, manifest = manifest)
}

# Individual-cell distributions pooled across samples, not sample means.
plot_stephenson_ucell_score_violins <- function(artifacts, annotation = "layer2") {
  annotations <- c("layer2", "layer3", "initial_clustering", "full_clustering")
  if (!is.character(annotation) || length(annotation) != 1L ||
      is.na(annotation) || !annotation %in% annotations) {
    stop("annotation must be layer2, layer3, initial_clustering or full_clustering.")
  }
  cells <- artifacts$cells
  thresholds <- artifacts$thresholds
  scores <- c(G1S = "cellCycle.G1S_UCell", G2M = "cellCycle.G2M_UCell",
              IFN = "IFN_UCell")
  if (!is.data.frame(cells) || !is.data.frame(thresholds) ||
      !all(c(annotation, unname(scores)) %in% names(cells)) ||
      !all(c("score", "threshold", "comparison") %in% names(thresholds)) ||
      nrow(thresholds) != 3L ||
      !setequal(thresholds$score, unname(scores)) ||
      anyDuplicated(thresholds$score) ||
      !is.numeric(thresholds$threshold) ||
      any(!is.finite(thresholds$threshold)) ||
      any(thresholds$threshold < 0 | thresholds$threshold > 1) ||
      anyNA(thresholds$comparison) || any(thresholds$comparison != ">")) {
    stop("Expected validated Stephenson cell-score artifacts and strict thresholds.")
  }
  category <- as.character(cells[[annotation]])
  category[is.na(category)] <- "(unclassified)"
  levels <- sort(unique(category), method = "radix")
  groups <- split(seq_len(nrow(cells)), factor(category, levels = levels))
  titles <- c(G1S = "G1/S cell cycle", G2M = "G2/M cell cycle", IFN = "IFN response")
  plots <- summaries <- list()
  for (key in names(scores)) {
    score <- scores[[key]]
    value <- cells[[score]]
    if (!is.numeric(value) || any(is.nan(value)) ||
        any(!is.na(value) & (!is.finite(value) | value < 0 | value > 1))) {
      stop(score, ": expected scores in [0, 1] or NA.")
    }
    cutoff <- thresholds$threshold[match(score, thresholds$score)]
    tab <- do.call(rbind, lapply(names(groups), function(type) {
      values <- value[groups[[type]]]
      finite <- values[is.finite(values)]
      n <- length(finite)
      data.frame(
        score = score, annotation = annotation, cell_type = type,
        threshold = cutoff, comparison = ">", n_cells = length(values),
        n_finite = n, n_missing = sum(is.na(values)),
        n_positive = sum(finite > cutoff),
        pct_positive = if (n) 100 * sum(finite > cutoff) / n else NA_real_,
        median = if (n) stats::median(finite) else NA_real_,
        minimum = if (n) min(finite) else NA_real_,
        maximum = if (n) max(finite) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
    tab$cell_type <- factor(tab$cell_type, levels = levels)
    # Constant/single-cell groups have a median marker instead of a spurious
    # density; entirely missing groups retain their axis position and table row.
    density_types <- as.character(tab$cell_type[
      tab$n_finite >= 2L & !is.na(tab$minimum) & tab$minimum < tab$maximum
    ])
    draw <- data.frame(cell_type = factor(category, levels = levels), score = value)
    draw <- draw[is.finite(draw$score) &
                   as.character(draw$cell_type) %in% density_types, , drop = FALSE]
    missing_types <- as.character(tab$cell_type[tab$n_finite == 0L])
    axis_labels <- setNames(levels, levels)
    axis_labels[missing_types] <- paste0(missing_types, " [no scores]")
    boundary <- paste0("Positive: score > ", format(cutoff, trim = TRUE))
    upper <- min(1, 1.1 * max(0.05, cutoff, value[is.finite(value)]))
    p <- ggplot2::ggplot() +
      ggplot2::geom_violin(
        data = draw, ggplot2::aes(x = cell_type, y = score),
        fill = "#80B1D3", colour = "#4D4D4D", linewidth = 0.25,
        scale = "width", trim = TRUE
      ) +
      ggplot2::geom_hline(yintercept = cutoff, linetype = "dashed",
                          colour = "#B2182B", linewidth = 0.4) +
      ggplot2::geom_point(
        data = tab[tab$n_finite > 0L, , drop = FALSE],
        ggplot2::aes(x = cell_type, y = median), size = 0.9
      ) +
      ggplot2::annotate("text", x = Inf, y = cutoff, label = boundary,
                        hjust = 1.02, vjust = -0.5, size = 3,
                        colour = "#B2182B") +
      ggplot2::scale_x_discrete(limits = levels, labels = axis_labels, drop = FALSE) +
      ggplot2::scale_y_continuous(limits = c(0, upper),
                                breaks = pretty(c(0, upper), n = 5)) +
      ggplot2::labs(
        title = paste("Stephenson:", titles[[key]]),
        subtitle = paste("Individual cells pooled across samples; grouped by", annotation),
        x = annotation, y = "UCell score",
        caption = paste(
          "Dashed line: strict positivity threshold (equality is not positive).",
          "Points: medians. NA scores are not imputed;",
          "[no scores] categories contain no finite scores.",
          "Positive percentages use finite-score cells only."
        )
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1,
                                           size = 8),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0)
      )
    plots[[key]] <- p
    tab$cell_type <- as.character(tab$cell_type)
    summaries[[key]] <- tab
  }
  list(plots = plots, summary = do.call(rbind, unname(summaries)))
}

# Exploratory author-HR split; does not change the cell-score manifest cutoffs.
# Unobserved states are omitted from feature_map (per_type records their zero
# counts), so no globally empty state creates a pseudocount-only CLR feature.
stephenson_ifn_split_composition <- function(cell_artifacts, bundle, threshold = 0.1) {
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 || threshold > 1) {
    stop("IFN threshold must be a finite numeric scalar in [0, 1].")
  }
  cells <- cell_artifacts$cells
  required <- c("Sample", "Status", "full_clustering", "IFN_UCell")
  if (!is.data.frame(cells) || !nrow(cells) ||
      !all(required %in% names(cells))) {
    stop("Stephenson IFN split requires Sample, Status, full_clustering and IFN_UCell.")
  }
  for (column in required[1:3]) {
    value <- cells[[column]]
    if ((!is.character(value) && !is.factor(value)) || anyNA(value) ||
        any(!nzchar(trimws(as.character(value)))) ||
        any(grepl("[[:cntrl:]]", value))) {
      stop(column, ": expected complete, nonempty metadata labels.")
    }
  }
  score <- cells$IFN_UCell
  if (!is.numeric(score) || any(!is.finite(score)) ||
      any(score < 0 | score > 1)) {
    stop("IFN_UCell must contain only finite scores in [0, 1]; missing scores cannot be classified.")
  }
  if (!is.list(bundle) ||
      !all(c("counts", "feat_mat", "labels") %in% names(bundle))) {
    stop("Stephenson IFN split requires an author-HR counts/feat_mat/labels bundle.")
  }
  counts <- as.matrix(bundle$counts)
  valid_names <- function(x) {
    !is.null(x) && length(x) > 0L && !anyNA(x) &&
      all(nzchar(trimws(x))) && !anyDuplicated(x) &&
      !any(grepl("[[:cntrl:]]", x))
  }
  if (!is.numeric(counts) || !nrow(counts) || !ncol(counts) ||
      !valid_names(rownames(counts)) || !valid_names(colnames(counts)) ||
      any(!is.finite(counts)) || any(counts < 0.5) ||
      any(counts - 0.5 != round(counts - 0.5))) {
    stop("Author-HR counts must have unique named axes and equal nonnegative integer raw counts + 0.5.")
  }
  samples <- rownames(counts)
  types <- colnames(counts)
  if (!identical(rownames(bundle$feat_mat), samples) ||
      !identical(colnames(bundle$feat_mat), types)) {
    stop("Author-HR counts and feat_mat sample/category axes must match exactly.")
  }
  labels <- bundle$labels
  if (!is.atomic(labels) || !is.null(dim(labels)) ||
      length(labels) != length(samples) || !valid_names(names(labels)) ||
      !identical(names(labels), samples) || anyNA(labels) ||
      !all(as.character(labels) %in% c("Covid", "Healthy"))) {
    stop("Author-HR labels must be named Covid/Healthy values in exact sample order.")
  }
  sample <- as.character(cells$Sample)
  category <- as.character(cells$full_clustering)
  if (!setequal(unique(sample), samples) || !setequal(unique(category), types)) {
    stop("Cell metadata samples/categories must match the complete author-HR bundle.")
  }
  sample_status <- unique(data.frame(
    Sample = sample, Status = as.character(cells$Status), stringsAsFactors = FALSE
  ))
  if (anyDuplicated(sample_status$Sample) ||
      !identical(sample_status$Status[match(samples, sample_status$Sample)],
                 unname(as.character(labels)))) {
    stop("Cell Status must be constant within each sample and match author-HR labels.")
  }
  original_raw <- counts - 0.5
  recovered <- table(
    factor(sample, levels = samples), factor(category, levels = types)
  )
  if (any(recovered != original_raw)) {
    stop("Cell metadata does not exactly recover the original author-HR raw counts.")
  }
  positive_names <- paste0(types, "_IFN_responding")
  if (any(positive_names %in% types) || anyDuplicated(positive_names)) {
    stop("IFN-responding feature names collide with original author-HR categories.")
  }
  # Labels are checked above only; assignment uses the strict score cutoff alone.
  responding <- score > threshold
  positive <- table(
    factor(sample[responding], levels = samples),
    factor(category[responding], levels = types)
  )
  negative <- original_raw - positive
  keep_negative <- colSums(negative) > 0
  keep_positive <- colSums(positive) > 0
  feature_map <- data.frame(
    feature = c(types[keep_negative], positive_names[keep_positive]),
    cell_type = c(types[keep_negative], types[keep_positive]),
    ifn_responding = c(rep(FALSE, sum(keep_negative)), rep(TRUE, sum(keep_positive))),
    stringsAsFactors = FALSE
  )
  raw_counts <- cbind(
    negative[, keep_negative, drop = FALSE],
    positive[, keep_positive, drop = FALSE]
  )
  dimnames(raw_counts) <- list(samples, feature_map$feature)
  collapsed <- t(rowsum(t(raw_counts), feature_map$cell_type, reorder = FALSE))
  collapsed <- collapsed[, types, drop = FALSE]
  if (any(raw_counts < 0) ||
      any(rowSums(raw_counts) != rowSums(original_raw)) ||
      any(collapsed != original_raw) || sum(raw_counts) != nrow(cells) ||
      sum(raw_counts[, feature_map$ifn_responding, drop = FALSE]) != sum(responding)) {
    stop("IFN split failed sample/category count conservation.")
  }
  per_type <- data.frame(
    cell_type = types, n_total = unname(colSums(original_raw)),
    n_ifn_responding = unname(colSums(positive)),
    n_nonresponding = unname(colSums(negative)), stringsAsFactors = FALSE
  )
  counts <- raw_counts + 0.5
  log_counts <- log(counts)
  list(
    raw_counts = raw_counts, counts = counts,
    feat_mat = log_counts - rowMeans(log_counts), labels = labels,
    per_type = per_type, feature_map = feature_map, threshold = threshold
  )
}
