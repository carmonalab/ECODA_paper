# Counts-only refresh of the eleven cohorts using one current zero-imputation policy.
zeroimp_datasets <- function() c("Adams", "Bassez", "Gongsharma_cmv_young_males", "Kfoury", "Kim", "Lee", "Pelka", "Smillie", "Stephenson", "Wu", "Zhang")
zeroimp_need <- function(ok, message) if (!isTRUE(ok)) stop(message, call.=FALSE)
zeroimp_menu <- function() {
  counts <- c("0.5", "2/3", "1", "5", "10", "20", "50", "100", "200")
  percent <- c("0.001", "0.01", "0.1", "1", "2", "5")
  c(unlist(lapply(counts, function(x) paste0(c("counts_zeros_", "counts_all_"), x)), use.names=FALSE),
    unlist(lapply(percent, function(x) paste0(c("percentage_all_", "percentage_zeros_", "multLN_", "multRepl_"), x, "%")), use.names=FALSE), "asinsqrt")
}
zeroimp_policy <- function() list(schema="ecoda_zeroimp_counts_five_v1", version=1L,
  recipe_menu=zeroimp_menu(), recipe="unchanged run_zeroimp_analysis ordinary count-screen body",
  parent_counts_offset=0.5, recipe_input="recovered nonnegative integer category counts; no offset before recipes",
  metrics=c("anosim_score", "mod_knn3_score", "cluster_score", "sil_score", "lisi_score"),
  scoring_policy="calc_sep_score_five_v1_lisi_separation", seed=123L,
  rng_kind=c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
  rng_policy="reset seed per cohort; restore RNG after each score call, so evaluation labels cannot affect subsequent imputation RNG",
  legacy_source9ca_comparability="NOT_USED_FOR_FINALS: original nine legacy and source9ca Stephenson/Bassez summary-only screens remain evidence; this refresh computes all eleven with one menu/RNG/full-sample policy",
  geometry_policy="PASS retains compact named distance and exact full-universe primary labels; dropped/extra/reordered samples cannot pass",
  failed_settings="explicit FAILED_OR_DROPPED_BY_ORDINARY_PRODUCER with actual scoring IDs when reached; no partial-universe or imputed scores; unhandled producer error fails the job")

zeroimp_control <- function(path, md5, config_path, datasets) {
  .analysis_checked_file(path, md5); control <- readRDS(path); .analysis_checked_file(path, md5)
  zeroimp_need(identical(control$schema,"ecoda_cell_depth_inputs_v1") && identical(control$status,"approved"), "Approved compact HR/metadata input release required")
  raw <- jsonlite::fromJSON(config_path, simplifyVector=FALSE)
  selected <- ecoda_select_screen_datasets(ecoda_raw_config_union(raw,"benchmark_union"),datasets)
  zeroimp_need(length(selected)>0L && all(selected%in%zeroimp_datasets()), "Only explicitly selected members of the eleven approved benchmark cohorts are allowed")
  zeroimp_need(is.list(control$parent_inventory), "Explicit checked author_hr_metadata_pairs.11 inventory required")
  inventory <- .analysis_checked_file(control$parent_inventory$path,control$parent_inventory$md5)
  parents <- read.delim(inventory,check.names=FALSE,stringsAsFactors=FALSE)
  required <- c("dataset","author_path","author_md5","author_bundle_key","metadata_path","metadata_md5","current_input_path","current_input_md5","expected_counts_policy")
  zeroimp_need(identical(names(parents),required) && identical(parents$dataset,zeroimp_datasets()), "Exact eleven-parent inventory schema/order required")
  zeroimp_need(is.list(control$legacy_parent_inventory),"Immutable original nine-parent inventory reference required")
  legacy <- read.delim(.analysis_checked_file(control$legacy_parent_inventory$path,control$legacy_parent_inventory$md5),check.names=FALSE,stringsAsFactors=FALSE)
  legacy_ids <- setdiff(zeroimp_datasets(),c("Bassez","Stephenson"))
  zeroimp_need(identical(names(legacy),required) && identical(legacy$dataset,legacy_ids) && identical(unname(as.matrix(parents[parents$dataset%in%legacy_ids,,drop=FALSE])),unname(as.matrix(legacy))),"Original nine-parent inventory lineage changed")
  expected <- expand.grid(dataset=zeroimp_datasets(),view="benchmark_analysis",method_key=c("ECODA_authors_HR","__dataset_metadata__"),stringsAsFactors=FALSE)
  zeroimp_need(setequal(.analysis_key(control$ledger),.analysis_key(expected)) && setequal(.analysis_key(control$expected_scope),.analysis_key(expected)),"Only the 22 reviewed count/metadata parents are allowed")
  for (i in seq_along(zeroimp_datasets())) {
    ds <- zeroimp_datasets()[i]; rows <- control$ledger[control$ledger$dataset==ds,,drop=FALSE]
    hr <- rows[rows$method_key=="ECODA_authors_HR",,drop=FALSE]; meta <- rows[rows$method_key=="__dataset_metadata__",,drop=FALSE]
    scope <- control$expected_scope[control$expected_scope$dataset==ds,,drop=FALSE]
    zeroimp_need(nrow(hr)==1L && nrow(meta)==1L && identical(hr$artifact_path,parents$author_path[i]) && identical(hr$md5,parents$author_md5[i]) && identical(hr$bundle_key,parents$author_bundle_key[i]) && identical(meta$artifact_path,parents$metadata_path[i]) && identical(meta$md5,parents$metadata_md5[i]) && identical(meta$bundle_key,"") && all(scope$current_input_md5==parents$current_input_md5[i]) && (!ds%in%c("Bassez","Stephenson") || all(scope$require_new)),"Parent inventory differs from reviewed ledger/current input")
  }
  list(control=control,raw=raw,selected=selected,parents=parents)
}
zeroimp_input <- function(path,md5,ds,raw) {
  config <- ecoda_screen_release_sources(path,md5,ds,raw)[[ds]]
  input <- .ecoda_load_local_composition(config)
  metadata <- readRDS(.analysis_checked_file(config$metadata_path,config$metadata_md5))
  truth <- metadata$labels
  zeroimp_need(!is.null(names(truth)) && !anyDuplicated(names(truth)) && setequal(names(truth),input$sample_ids) && identical(unname(as.character(truth[input$sample_ids])),unname(as.character(input$labels))),"Reviewed HR and primary metadata labels differ by Sample")
  list(config=config,input=input)
}
zeroimp_run_cohort <- function(input,score_fun=calc_sep_score,recipe_fun=run_zeroimp_analysis) {
  counts <- input$counts; labels <- input$labels; policy <- zeroimp_policy()
  zeroimp_need(is.matrix(counts) && is.numeric(counts) && all(is.finite(counts)) && all(counts>=0) && all(counts==floor(counts)) && all(rowSums(counts)>0) && identical(unname(names(labels)),unname(rownames(counts))),"Recovered raw integer counts/labels required")
  old_kind <- RNGkind(); old_seed <- if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) get(".Random.seed",.GlobalEnv) else NULL
  on.exit({do.call(RNGkind,as.list(old_kind));if(is.null(old_seed)){if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)}else assign(".Random.seed",old_seed,.GlobalEnv)},add=TRUE)
  do.call(RNGkind,as.list(policy$rng_kind));set.seed(policy$seed);rng_before <- .Random.seed
  # Shallow-copy the existing recipe function; isolate scoring RNG only. The
  # imputation/transformation body and its historical try/drop behavior stay intact.
  env <- new.env(parent=environment(recipe_fun));env$calc_sep_score <- function(dist_mat,labels,...) {
    kind <- RNGkind(); state <- get(".Random.seed",.GlobalEnv)
    on.exit({do.call(RNGkind,as.list(kind));assign(".Random.seed",state,.GlobalEnv)},add=TRUE)
    distance <- as.matrix(dist_mat); ids <- unname(rownames(distance)); expected <- unname(rownames(counts))
    geometry <- list(sample_ids=ids,dropped_ids=setdiff(expected,ids),extra_ids=setdiff(ids,expected),labels=labels)
    full <- identical(ids,expected) && identical(unname(colnames(distance)),ids) &&
      identical(unname(names(labels)),ids) && identical(as.character(labels),as.character(input$labels)) &&
      is.numeric(distance) && all(is.finite(distance)) && all(distance>=-1e-10) &&
      max(abs(distance-t(distance)))<=1e-10 && max(abs(diag(distance)))<=1e-10
    # No sample-drop exception applies to zero-imputation (scITD alone has one).
    if(full) {
      geometry$dist_mat <- stats::as.dist(distance)
      score <- tryCatch(score_fun(dist_mat,labels,...,parallel=1),error=function(e)structure(list(),score_error=conditionMessage(e)))
    } else score <- structure(list(),score_error="Imputer/scorer did not preserve exact full sample universe, primary labels or named distance geometry")
    attr(score,"zeroimp_geometry") <- geometry
    score
  }
  environment(recipe_fun) <- env
  messages <- character(); values <- withCallingHandlers(recipe_fun(as.data.frame(counts),labels),
    warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart("muffleWarning")})
  zeroimp_need(is.list(values) && !is.null(names(values)) && !anyDuplicated(names(values)) && all(names(values)%in%policy$recipe_menu),"Ordinary producer emitted unknown/duplicate settings")
  records <- setNames(vector("list",length(policy$recipe_menu)),policy$recipe_menu)
  rows <- vector("list",length(records))
  for(i in seq_along(records)) {
    key <- names(records)[i]; score <- values[[key]]
    if(startsWith(key,"multLN_") && is.list(score) && identical(names(score),"multLN")) score <- score$multLN
    geometry <- attr(score,"zeroimp_geometry");score_error <- attr(score,"score_error")
    attr(score,"zeroimp_geometry") <- NULL;attr(score,"score_error") <- NULL
    ok <- is.list(geometry) && !is.null(geometry$dist_mat) && is.list(score) && !anyDuplicated(names(score)) && all(policy$metrics%in%names(score)) && all(vapply(score[policy$metrics],function(x)is.numeric(x)&&length(x)==1L&&is.finite(x),logical(1)))
    if(ok) ok <- score$lisi_score>=0 && score$lisi_score<=1
    status <- if(ok) "PASS" else "FAILED_OR_DROPPED_BY_ORDINARY_PRODUCER"
    records[[i]] <- list(setting=key,status=status,scores=if(ok)score else NULL,
      sample_ids=geometry$sample_ids,dropped_ids=geometry$dropped_ids,extra_ids=geometry$extra_ids,
      dist_mat=if(ok)geometry$dist_mat else NULL,labels=if(ok)geometry$labels else NULL,
      distance_md5=if(ok).ecoda_screen_digest(geometry$dist_mat) else NULL,
      reason=if(ok)"" else if(!is.null(score_error))score_error else "No complete finite canonical-five result/geometry returned by ordinary producer; scoring IDs unknown when not reached")
    row <- data.frame(setting=key,status=status,reason=records[[i]]$reason,scoring_sample_ids=as.character(jsonlite::toJSON(geometry$sample_ids,auto_unbox=FALSE,null="null")),dropped_sample_ids=as.character(jsonlite::toJSON(geometry$dropped_ids,auto_unbox=FALSE,null="null")),stringsAsFactors=FALSE)
    for(metric in policy$metrics)row[[metric]] <- if(ok)score[[metric]] else NA_real_
    rows[[i]] <- row
  }
  list(results=do.call(rbind,rows),recipes=records,counts=counts,labels=labels,
    counts_md5=.ecoda_screen_digest(counts),labels_md5=.ecoda_screen_digest(labels),
    policy=policy,policy_md5=.ecoda_screen_digest(policy),rng_before=rng_before,rng_after=.Random.seed,
    recipe_body_md5=.ecoda_screen_digest(body(recipe_fun)),warnings=unique(messages))
}
zeroimp_validate_package <- function(package,input,source,run_id,config_md5,runner_md5,release_md5) {
  p <- zeroimp_policy()
  zeroimp_need(identical(package$schema,p$schema) && identical(package$dataset,source$dataset) && identical(package$policy,p) && identical(package$run_id,run_id) && identical(package$source,source) && identical(package$config_md5,config_md5) && identical(package$runner_md5,runner_md5) && identical(package$release_md5,release_md5),"Zeroimp package source/code/policy identity differs")
  zeroimp_need(identical(package$counts,input$counts) && identical(package$labels,input$labels) && identical(package$counts_md5,.ecoda_screen_digest(input$counts)) && identical(package$labels_md5,.ecoda_screen_digest(input$labels)) && identical(package$policy_md5,.ecoda_screen_digest(p)),"Zeroimp count/label/policy fingerprint differs")
  zeroimp_need(identical(package$results$setting,p$recipe_menu) && identical(names(package$recipes),p$recipe_menu) && all(package$results$status%in%c("PASS","FAILED_OR_DROPPED_BY_ORDINARY_PRODUCER")) && is.integer(package$rng_before) && length(package$rng_before)==7L && is.integer(package$rng_after) && length(package$rng_after)==7L && identical(package$recipe_body_md5,.ecoda_screen_digest(body(run_zeroimp_analysis))),"Complete recipe/RNG/body identities required")
  zeroimp_need(identical(package$recipes$counts_all_0.5$status,"PASS"),"Reference counts_all_0.5 setting must pass; an all-failed screen cannot be published")
  for(i in seq_along(p$recipe_menu)) {
    record <- package$recipes[[i]]; row <- package$results[i,,drop=FALSE]
    zeroimp_need(identical(record$setting,p$recipe_menu[i]) && identical(record$status,row$status) && identical(record$reason,row$reason),"Setting status identity differs")
    zeroimp_need(identical(row$scoring_sample_ids,as.character(jsonlite::toJSON(record$sample_ids,auto_unbox=FALSE,null="null"))) && identical(row$dropped_sample_ids,as.character(jsonlite::toJSON(record$dropped_ids,auto_unbox=FALSE,null="null"))),"Scoring/dropped sample evidence differs")
    if(row$status=="PASS") {
      distance <- as.matrix(record$dist_mat)
      zeroimp_need(inherits(record$dist_mat,"dist") && identical(record$sample_ids,unname(rownames(input$counts))) && !length(record$dropped_ids) && !length(record$extra_ids) && identical(unname(rownames(distance)),unname(rownames(input$counts))) && identical(unname(colnames(distance)),unname(rownames(input$counts))) && all(is.finite(distance)) && all(distance>=-1e-10) && max(abs(distance-t(distance)))<=1e-10 && max(abs(diag(distance)))<=1e-10 && identical(record$labels,input$labels) && identical(record$distance_md5,.ecoda_screen_digest(record$dist_mat)),"PASS distance/labels/full sample universe fingerprint differs")
      zeroimp_need(all(vapply(p$metrics,function(m)is.numeric(record$scores[[m]])&&length(record$scores[[m]])==1L&&is.finite(record$scores[[m]])&&identical(record$scores[[m]],row[[m]]),logical(1))) && row$lisi_score>=0 && row$lisi_score<=1,"PASS requires exact finite canonical five scores")
    } else zeroimp_need(is.null(record$scores) && is.null(record$dist_mat) && is.null(record$labels) && all(is.na(row[p$metrics])) && nzchar(record$reason),"Dropped setting cannot fabricate scores/geometry")
  }
  invisible(TRUE)
}
zeroimp_reuse_plan <- function(control) {
  refs <- control$zeroimp_checkpoint_reuse
  if(is.null(refs))return(NULL)
  zeroimp_need(identical(refs$schema,"ecoda_zeroimp_checkpoint_reuse_v1")&&identical(sort(names(refs)),sort(c("schema","review","science","inventory"))),"Unknown zeroimp checkpoint reuse schema")
  checked <- function(ref) {
    zeroimp_need(is.list(ref)&&identical(sort(names(ref)),c("path","sha256"))&&grepl("^[a-f0-9]{64}$",ref$sha256),"Invalid checkpoint evidence reference")
    zeroimp_need(is.character(ref$path)&&length(ref$path)==1L&&grepl("\\.json$",ref$path)&&file.exists(ref$path)&&!dir.exists(ref$path)&&file.info(ref$path)$size<8*1024^2,"Checkpoint evidence must be bounded JSON metadata")
    path <- .ecoda_validate_local_file(ref$path,"checkpoint metadata")$PATH
    zeroimp_need(file.info(path)$size<8*1024^2&&identical(digest::digest(file=path,algo="sha256"),ref$sha256),"Checkpoint evidence changed")
    jsonlite::fromJSON(path,simplifyVector=FALSE)
  }
  review<-checked(refs$review);science<-checked(refs$science);inventory<-checked(refs$inventory)
  datasets<-zeroimp_datasets()[1:4];old_run<-"stage5_p6_zeroimp11_20261110_v1"
  zeroimp_need(identical(review$schema,"ecoda_zeroimp_checkpoint4_root_reuse_review_v1")&&identical(review$status,"APPROVED_FIXED_FOUR_TYPED_CHECKPOINT_REUSE")&&isTRUE(review$review_complete)&&isTRUE(review$typed_checkpoint_science_validated)&&isTRUE(review$reuse_without_imputation_or_scoring)&&identical(review$generic_failed_input_bypass_authorized,FALSE)&&identical(review$original_gate_approved,FALSE)&&identical(review$original_gate_state,"FAILED")&&identical(review$original_stage5_owners,"FAIL_PRESERVED")&&identical(review$original_published_records,"ABSENT_NOT_FABRICATED")&&identical(review$original_run_id,old_run)&&identical(review$original_source_commit,"3736052f4451ca4dda80684a4af840dac5f37f4b")&&identical(review$original_control_md5,"e14f1b40ba776731c86f2e1c421e8b64")&&identical(unlist(review$datasets,use.names=FALSE),datasets)&&review$recipes_each==43L&&identical(review$science_report_sha256,refs$science$sha256)&&identical(review$byte_inventory_sha256,refs$inventory$sha256),"Only ROOT-reviewed original four zeroimp checkpoints may be reused")
  zeroimp_need(identical(science$schema,"ecoda_zeroimp_checkpoint4_science_observation_v1")&&identical(science$status,"VALIDATED_CHECKPOINTS_PENDING_ROOT_REUSE_REVIEW")&&identical(names(science$results4),datasets)&&isTRUE(science$unpublished_checkpoint_records_absent)&&identical(science$scientific_reuse_approved,FALSE)&&identical(inventory$schema,"ecoda_zeroimp11_partial_byte_inventory_v1"),"Checkpoint scientific evidence differs")
  files<-Filter(function(x)isTRUE(x$exists),inventory$rows)
  zeroimp_need(identical(review$new_model_fits,FALSE)&&identical(science$original_source_commit,review$original_source_commit)&&identical(science$original_gate,"FAILED_PRESERVED")&&identical(science$original_owners,"FAIL_PRESERVED")&&science$model_fits==0&&science$score_calculations==0&&science$H5AD_reads==0,"Checkpoint review cannot promote the original gate or claim new computation")
  zeroimp_need(length(inventory$rows)==22L&&length(files)==8L,"Exactly eight original checkpoint files in the eleven-cohort inventory are required")
  old_root<-dirname(dirname(dirname(files[[1]]$path)))
  # path is <run>/artifacts/zeroimp_count_screen/analysis/<file>.
  old_root<-dirname(old_root)
  zeroimp_need(identical(basename(old_root),old_run)&&identical(basename(dirname(old_root)),"_ecoda_runs"),"Wrong checkpoint producer namespace")
  for(ds in datasets) {
    result<-science$results4[[ds]]
    zeroimp_need(identical(result$run_id,old_run)&&identical(result$release_md5,review$original_control_md5)&&isTRUE(result$full_PASS_geometry_and_five_validated)&&length(result$recipes43)==43L,"Incomplete validated checkpoint cohort")
    for(ext in c("rds","csv"))zeroimp_need(sum(vapply(files,function(x)identical(x$path,file.path(old_root,"artifacts/zeroimp_count_screen/analysis",paste0(ds,"_zeroimp_counts.",ext))),logical(1)))==1L,"Missing/duplicate checkpoint bytepair")
  }
  list(refs=refs,datasets=datasets,run_root=old_root,files=files,science=science)
}
zeroimp_check_reuse <- function(reuse) {
  if(is.null(reuse))return(invisible(TRUE))
  for(x in reuse$science$checked_before_after) {
    zeroimp_need(!grepl("\\.(h5ad|h5|sif)$",x$path,ignore.case=TRUE),"Checkpoint evidence cannot read RNA/runtime images")
    path<-.ecoda_validate_local_file(x$path,"checkpoint source/proof/owner")$PATH
    zeroimp_need(identical(digest::digest(file=path,algo="sha256"),x$sha256),"Checkpoint source/proof/owner changed since scientific review")
  }
  for(x in reuse$files) {
    proof<-.ecoda_validate_output_sidecar(x$path)
    zeroimp_need(identical(proof$MD5,x$md5)&&as.numeric(proof$SIZE)==x$stat$st_size,"Checkpoint byte identity differs")
    key<-substr(digest::digest(x$path,algo="sha256",serialize=FALSE),1,32)
    record<-file.path(reuse$run_root,"manifests/artifacts",paste0(key,".record"))
    zeroimp_need(!file.exists(record)&&!.ecoda_is_symlink(record),"Unexpected checkpoint producer record")
  }
  invisible(TRUE)
}
zeroimp_produce_package <- function(input,source,ds,rds,csv,identity,reuse=NULL,run_fun=zeroimp_run_cohort) {
  reused<-!is.null(reuse)&&ds%in%reuse$datasets
  if(reused) {
    zeroimp_check_reuse(reuse)
    paths<-setNames(vapply(c("rds","csv"),function(ext)file.path(reuse$run_root,"artifacts/zeroimp_count_screen/analysis",paste0(ds,"_zeroimp_counts.",ext)),character(1)),c("rds","csv"))
    original<-reuse$science$results4[[ds]];package<-readRDS(paths[["rds"]])
    identity$run_id<-original$run_id;identity$runner_md5<-original$runner_md5;identity$release_md5<-original$release_md5
  } else {
    package<-run_fun(input);package$schema<-zeroimp_policy()$schema;package$dataset<-ds;package$source<-source
    for(k in names(identity))package[[k]]<-identity[[k]]
    package$origin<-"new_counts_only_refresh_from_reviewed_parents";package$parent_origin<-"retained in checked release receipts; not promoted to new model outputs"
  }
  zeroimp_validate_package(package,input,source,identity$run_id,identity$config_md5,identity$runner_md5,identity$release_md5)
  if(reused) {
    zeroimp_need(identical(package$csv_md5,.ecoda_md5(paths[["csv"]])),"Reused checkpoint CSV identity differs")
    for(ext in c("rds","csv")) {
      target<-if(ext=="rds")rds else csv;temporary<-paste0(target,".tmp.",Sys.getpid())
      zeroimp_need(!file.exists(target)&&!.ecoda_is_symlink(target),"Checkpoint copy destination exists")
      on.exit(unlink(temporary),add=TRUE)
      zeroimp_need(file.copy(paths[[ext]],temporary,overwrite=FALSE)&&identical(.ecoda_md5(temporary),.ecoda_md5(paths[[ext]]))&&file.link(temporary,target),"Create-only checkpoint copy failed")
      unlink(temporary)
      .ecoda_write_text_atomic(paste0(target,".md5"),c(paste0("MD5=",.ecoda_md5(target)),paste0("SIZE=",file.info(target)$size),paste0("PATH=",target)))
    }
  } else {.ecoda_save_csv(package$results,csv);package$csv_md5<-.ecoda_md5(csv);.ecoda_save_rds(package,rds)}
  invisible(package)
}
zeroimp_main <- function(args=commandArgs(TRUE)) {
  zeroimp_need(identical(as.character(getRversion()),"4.5.2"),"Configured R4.5.2 required")
  zeroimp_need(length(args)%%2L==0L,"Arguments require values")
  keys <- args[seq.int(1,length(args),2)];zeroimp_need(!anyDuplicated(keys),"Duplicate flags")
  allowed <- c("--project_root","--analysis_release","--analysis_release_md5","--datasets","--config_path","--bind","--produce","--validate","--run_id")
  zeroimp_need(all(keys%in%allowed),"Unknown zeroimp flag");a <- setNames(as.list(args[seq.int(2,length(args),2)]),sub("^--","",keys))
  zeroimp_need(all(c("project_root","analysis_release","analysis_release_md5","datasets","config_path","run_id")%in%names(a)) && sum(c("bind","produce","validate")%in%names(a))==1L,"Exact bound zeroimp action required")
  project <- normalizePath(a$project_root,mustWork=TRUE)
  zeroimp_need(identical(Sys.getenv("ECODA_SOURCE_SNAPSHOT_REQUIRED"),"1") && identical(Sys.getenv("ECODA_SOURCE_ROOT"),project) && identical(Sys.getenv("ECODA_RUN_ID"),a$run_id),"Immutable run-bound zeroimp executor required")
  source(file.path(project,"src/utils/analysis_release.R"));source(file.path(project,"src/5_run_benchmark_methods/run_local_ecoda_derived.R"));source(file.path(project,"src/5_run_benchmark_methods/benchmark_pipeline.R"))
  control <- zeroimp_control(a$analysis_release,a$analysis_release_md5,a$config_path,a$datasets)
  reuse<-zeroimp_reuse_plan(control$control)
  if(!is.null(reuse))zeroimp_need(identical(control$selected,zeroimp_datasets()),"Checkpoint recovery retains the complete eleven-cohort denominator")
  if(!is.null(a$bind)) {
    paths <- c(a$analysis_release,control$control$parent_inventory$path,control$control$legacy_parent_inventory$path);hashes <- c(a$analysis_release_md5,control$control$parent_inventory$md5,control$control$legacy_parent_inventory$md5)
    for(ds in control$selected) {
      rows <- control$control$ledger[control$control$ledger$dataset==ds,,drop=FALSE]
      paths <- c(paths,rows$artifact_path,rows$receipt_path);hashes <- c(hashes,rows$md5,rows$receipt_md5)
      for(i in seq_len(nrow(rows))) {
        receipt <- readRDS(.analysis_checked_file(rows$receipt_path[i],rows$receipt_md5[i]));paths <- c(paths,receipt$evidence_path);hashes <- c(hashes,receipt$evidence_md5)
      }
    }
    if(!is.null(reuse))for(ref in reuse$refs[c("review","science","inventory")]) {paths<-c(paths,ref$path);hashes<-c(hashes,.ecoda_md5(ref$path))}
    inventory <- unique(data.frame(path=paths,md5=hashes));zeroimp_need(!anyDuplicated(inventory$path),"Conflicting input identities")
    for(i in seq_len(nrow(inventory))).analysis_checked_file(inventory$path[i],inventory$md5[i])
    write.table(inventory,file.path(a$bind,"inputs.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(data.frame(control$selected,"benchmark_analysis","zeroimp_count_screen"),file.path(a$bind,"selection.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    .ecoda_write_json_atomic(list(schema="ecoda_zeroimp_checkpoint_reuse_binding_v1",reuse=if(is.null(reuse))NULL else reuse$refs),file.path(a$bind,"checkpoint_reuse.json"))
    return(invisible(NULL))
  }
  output <- a$produce %||% a$validate
  if(!is.null(a$produce)) {
    zeroimp_need(!file.exists(output) && !.ecoda_is_symlink(output),"Zeroimp output is create-only")
    dir.create(output,recursive=TRUE,showWarnings=FALSE)
    suppressPackageStartupMessages({library(dplyr);library(survival)})
    source(file.path(project,"src/utils/math_utils.R"));source(file.path(project,"src/utils/scoring_metrics.R"))
    options(mc.cores=1L)
  }
  runner_md5 <- .ecoda_md5(file.path(project,"src/notebooks/run_zeroimp_count_screen.R"));config_md5 <- .ecoda_md5(a$config_path)
  zeroimp_check_reuse(reuse)
  for(ds in control$selected) {
    parent <- zeroimp_input(a$analysis_release,a$analysis_release_md5,ds,control$raw);input <- parent$input;source <- c(list(dataset=ds),parent$config)
    parent_rows <- control$control$ledger[control$control$ledger$dataset==ds,,drop=FALSE]
    source$parent_receipts <- lapply(seq_len(nrow(parent_rows)),function(i) {
      row <- parent_rows[i,,drop=FALSE];receipt <- readRDS(.analysis_checked_file(row$receipt_path,row$receipt_md5));.analysis_checked_file(row$receipt_path,row$receipt_md5)
      list(method_key=row$method_key,artifact_path=row$artifact_path,artifact_md5=row$md5,
        receipt_path=row$receipt_path,receipt_md5=row$receipt_md5,origin=receipt$origin,
        input_md5=receipt$input_md5,source_policy=receipt$source_policy)
    })
    rds <- file.path(output,paste0(ds,"_zeroimp_counts.rds"));csv <- file.path(output,paste0(ds,"_zeroimp_counts.csv"))
    if(!is.null(a$produce)) {
      package<-zeroimp_produce_package(input,source,ds,rds,csv,list(run_id=a$run_id,config_md5=config_md5,runner_md5=runner_md5,release_md5=a$analysis_release_md5),reuse)
    } else {
      .ecoda_validate_output_sidecar(rds);.ecoda_validate_output_sidecar(csv);package <- readRDS(rds)
      original<-if(!is.null(reuse)&&ds%in%reuse$datasets)reuse$science$results4[[ds]] else NULL
      zeroimp_validate_package(package,input,source,if(is.null(original))a$run_id else original$run_id,config_md5,if(is.null(original))runner_md5 else original$runner_md5,if(is.null(original))a$analysis_release_md5 else original$release_md5)
      zeroimp_need(identical(package$csv_md5,.ecoda_md5(csv)) && isTRUE(all.equal(read.csv(csv,check.names=FALSE,colClasses=c(setting="character",status="character",reason="character",scoring_sample_ids="character",dropped_sample_ids="character")),package$results,check.attributes=FALSE,tolerance=1e-10)),"Zeroimp CSV/RDS differs")
    }
    rm(input,parent,package);gc()
  }
  zeroimp_check_reuse(reuse)
  if(!is.null(a$produce)) {
    .ecoda_write_json_atomic(list(schema="ecoda_zeroimp_counts_manifest_v1",run_id=a$run_id,datasets=control$selected,recipe_menu=zeroimp_menu(),scoring_policy=zeroimp_policy()$scoring_policy,source_manifest_md5=.ecoda_md5(Sys.getenv("ECODA_SOURCE_MANIFEST")),release_md5=a$analysis_release_md5,legacy_source9ca_comparability=zeroimp_policy()$legacy_source9ca_comparability,checkpoint_reuse=if(is.null(reuse))NULL else list(refs=reuse$refs,datasets=reuse$datasets,action="BYTE_IDENTICAL_COPY_NO_IMPUTATION_OR_SCORING",original_embedded_identity_preserved=TRUE)),file.path(output,"screen_manifest.json"))
    .ecoda_write_json_atomic(list(state="OK",run_id=a$run_id,cohorts=length(control$selected),settings_per_cohort=43L),file.path(output,"status.json"))
  }
  manifest <- jsonlite::fromJSON(file.path(output,"screen_manifest.json"),simplifyVector=FALSE)
  status <- jsonlite::fromJSON(file.path(output,"status.json"),simplifyVector=FALSE)
  lineage<-manifest$checkpoint_reuse
  zeroimp_need(if(is.null(reuse))is.null(lineage) else is.list(lineage)&&identical(lineage$refs,reuse$refs)&&identical(unlist(lineage$datasets,use.names=FALSE),reuse$datasets)&&identical(lineage$action,"BYTE_IDENTICAL_COPY_NO_IMPUTATION_OR_SCORING")&&isTRUE(lineage$original_embedded_identity_preserved),"Checkpoint copy lineage differs")
  zeroimp_need(identical(manifest$schema,"ecoda_zeroimp_counts_manifest_v1") && identical(manifest$run_id,a$run_id) && identical(unlist(manifest$datasets,use.names=FALSE),control$selected) && identical(unlist(manifest$recipe_menu,use.names=FALSE),zeroimp_menu()) && identical(manifest$scoring_policy,zeroimp_policy()$scoring_policy) && identical(manifest$release_md5,a$analysis_release_md5) && identical(manifest$source_manifest_md5,.ecoda_md5(Sys.getenv("ECODA_SOURCE_MANIFEST"))) && identical(manifest$legacy_source9ca_comparability,zeroimp_policy()$legacy_source9ca_comparability) && identical(status$state,"OK") && identical(status$run_id,a$run_id) && status$cohorts==length(control$selected) && status$settings_per_cohort==43L,"Zeroimp manifest/status/complete recipe denominator differs")
  invisible(NULL)
}
if(sys.nframe()==0L)zeroimp_main()
