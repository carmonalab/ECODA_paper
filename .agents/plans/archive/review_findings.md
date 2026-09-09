 Review findings after implementation of refined plan in /Users/christianhalter/.omp/agent/sessions/-Desktop-ECODA_paper/2026-08-27T08-17-45-367Z_01a0424b-e597-7272-8e6a-23ffe1c37aab/local/onboarding-uncorrected-two-stage-plan.md
 
(start review findings by agent)

 Verdict: REQUEST_CHANGES                                                                                                                                                                          
                                                                                                                                                                                                   
 Stage 1 is not safe to execute. The review found multiple P1 blockers in scheduler recovery, artifact validation, annotation reuse, and the required exact Stage 5 path.                          
                                                                                                                                                                                                   
 P1 findings                                                                                                                                                                                       
                                                                                                                                                                                                   
 1. Exact Stage 5 batch runs always fail during final validation                                                                                                                                   
    src/5_run_benchmark_methods/1_submit_hpc_array.sh:774-783                                                                                                                                      
    Exact matrix rows use:                                                                                                                                                                         
    ```text                                                                                                                                                                                        
DATASET<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
    ```                                                                                                                                                                                            
    but LABELS contains method names. Final Python/R validators and the sync path treat the third field as a method label. They reject batch_effect_uncorrected as an unknown artifact label after 
    all arrays complete. The required uncorrected run cannot reach synchronization.                                                                                                                
    Fix: keep the third field as pass/view scope and create a method-scoped manifest for final artifact validation/sync, or make validators consume LABELS while separately validating the pass    
    scope.                                                                                                                                                                                         
 2. Ordinary trans/zeroimp outputs are rejected by the RDS validator                                                                                                                               
    src/5_run_benchmark_methods/validate_benchmark_rds_contract.R:126-138                                                                                                                          
    validate_result_file() requires every bundle element to contain scores, feat_mat, dist_mat, and labels. The actual transformation and zero-imputation workers emit different result shapes.    
    Ordinary analyses therefore fail the final validator despite successful worker output.                                                                                                         
    Fix: add analysis-specific validators or stop routing these outputs through the benchmark-combo validator.                                                                                     
 3. Stage 5 RDS validation is fail-open on scientific integrity                                                                                                                                    
    src/5_run_benchmark_methods/validate_benchmark_rds_contract.R:86-99,123-138,165-172                                                                                                            
    The validator:                                                                                                                                                                                 
   - uses set equality instead of ordered sample comparison;                                                                                                                                       
   - accepts non-finite feature, distance, and score values;                                                                                                                                       
   - does not require method/pass-specific combo keys;                                                                                                                                             
   - misinterprets the exact batch scope as an artifact label.                                                                                                                                     
    A checksum-valid but reordered, malformed, or incomplete bundle can pass or the required exact run can fail.                                                                                   
 4. Matrix validator accepts malformed reusable Feather artifacts                                                                                                                                  
    src/5_run_benchmark_methods/matrix_artifact_validator.py:36-73                                                                                                                                 
    It checks nonempty rows and identifiers but does not enforce square shape, numeric values, finite values, or ordered row/column identity. Rectangular or all-NaN distance outputs can be       
    accepted and synchronized.                                                                                                                                                                     
 5. Python benchmark workers fail after writing method results                                                                                                                                     
    src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py:247-249                                                                                          
    log_execution_time() creates an indexless RangeIndex frame, then sends it through atomic_to_feather(), whose validator requires sample identifiers. Python method output can be produced       
    successfully but the task exits nonzero while writing its execution log.                                                                                                                       
    Fix: use a separate execution-log writer/validator that accepts the documented log schema.                                                                                                     
 6. Stage 5 permits a pass/view mismatch                                                                                                                                                           
    src/5_run_benchmark_methods/1_submit_hpc_array.sh:174-183                                                                                                                                      
    A non-exact pass selection file can contain benchmark_analysis while --pass uncorrected or --pass corrected is active. The R worker accepts the row and runs ordinary embeddings under a       
    pass-qualified output namespace, silently mixing analysis regimes.                                                                                                                             
    Fix: require VIEW == batch_effect_${PASS} for pass selections, regardless of exact mode.                                                                                                       
 7. Stage 5 transient retries use the wrong error-log prefix                                                                                                                                       
    src/5_run_benchmark_methods/1_submit_hpc_array.sh:561-568                                                                                                                                      
    Matrix jobs write errors under 5_matrix_*, but the worker environment omits JOB_LOG_PREFIX. Workers look for their old default prefixes, so transient failures are not detected/requeued.      
 8. Valid Myocardial prerequisites are rejected                                                                                                                                                    
    src/utils/py/derived_prerequisite_contract.py:39-44; called from src/2_dataset_specific_preprocessing/1_submit_hpc.sh:242-245                                                                  
    The semantic validator requires obs["Sample"]. Raw Myocardial data uses the registry sample column orig_ident; the reconstruction helper only adds layers["counts"]. A valid reconstructed     
    artifact therefore never satisfies the prerequisite gate.                                                                                                                                      
    Fix: validate only the required counts layer, or validate against the configured raw sample column.                                                                                            
 9. Joanito semantic validation is broken by shell quoting                                                                                                                                         
    src/utils/bash/ecoda_run_common.sh:495-506                                                                                                                                                     
    The R program is embedded in a shell single-quoted string but contains literals such as "5' seq" and "3' seq". Those apostrophes terminate the shell string and split the command. Correct     
    Joanito metadata will fail semantic validation.                                                                                                                                                
 10. Stage 2 can race CombinedPBMC when steps are supplied out of order                                                                                                                            
     src/2_dataset_specific_preprocessing/1_submit_hpc.sh:287-300                                                                                                                                  
     If --steps combinedpbmc,gongsharma_cap is supplied, CAP_JOB_ID is empty when CombinedPBMC is submitted, so no afterok dependency is attached. CombinedPBMC can read source files while the    
     cap hook overwrites them.                                                                                                                                                                     
     Fix: normalize gongsharma_cap before combinedpbmc, independent of caller order.                                                                                                               
 11. Stage 2 has no terminal success state for real runs                                                                                                                                           
     src/2_dataset_specific_preprocessing/1_submit_hpc.sh:369-374                                                                                                                                  
     After the watchdog and scheduler manifest complete, the normal path does not call ecoda_set_run_state OK. Only test mode writes terminal success. A successful real run remains without       
     status/terminal, and its own --sync-only path rejects it.                                                                                                                                     
 12. Stage 2 scheduler failures leave runs and owners active                                                                                                                                       
     src/2_dataset_specific_preprocessing/1_submit_hpc.sh:296-301,309-317; stage2_watchdog.sh:197-203                                                                                              
     set -e surrounds both hook and retry sbatch calls. A scheduler rejection exits before the explicit failure branches, leaving no durable FAIL status and leaving already-acquired owners in    
     ACTIVE.                                                                                                                                                                                       
 13. Stage 2 scheduler log directories are not created                                                                                                                                             
     src/2_dataset_specific_preprocessing/1_submit_hpc.sh:287-294                                                                                                                                  
     The dispatcher submits --output/--error paths under ${LOGS_DIR} without creating the directory. A clean checkout may have no logs/ directory, causing jobs to fail before execution.          
 14. Stage 3 --sync-only breaks documented recovery                                                                                                                                                
     src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh:121-132                                                                                                                                    
     The previous recovery contract accepts scheduler/array IDs and is intended to bypass a false watchdog failure. The new branch interprets the argument as an ECODA run ID and requires         
     terminal STATE=OK, so recovery from a failed or interrupted run cannot work.                                                                                                                  
 15. Stage 3 synchronization can report success with stale NAS output                                                                                                                              
     src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh:112-114                                                                                                                                    
     rsync failure is not guarded. If NAS already contains an older valid file/sidecar pair, the subsequent remote checksum check validates that stale pair against itself, not against the local  
     source.                                                                                                                                                                                       
 16. Stage 3 OOM retry transient handling is broken                                                                                                                                                
     src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh:141-147                                                                                                                               
     Retry jobs write errors under a retry-specific prefix, but PREPROCESS_ERROR_PREFIX is not exported. Worker transient-error detection looks at a nonexistent path. Additionally, retry sbatch  
     is executed under set -e, so scheduler rejection skips fail(), durable failure status, and owner cleanup.                                                                                     
 17. Forced Stage 4 runs delete the input checksum required by merge                                                                                                                               
     src/4_cell_type_annotation/1_submit_onboarding_stage.sh:377-379                                                                                                                               
     --force calls ecoda_invalidate_artifact() on selected preprocessed h5ads, deleting their .md5 sidecars. 3.2_merge_worker.sh now requires the input h5ad checksum before merging. Every forced 
     Stage 4 run therefore fails at merge.                                                                                                                                                         
 18. Multi-view annotation unions cannot be merged                                                                                                                                                 
     src/4_cell_type_annotation/3.2_merge_worker.sh:30-37; 3.1_merge_annotations.py                                                                                                                
     Preparation builds one annotation feather set over the union of all selected views. Merge then validates that same feather set against each individual view's keys. For a view that is a      
     subset of the union, such as Stephenson's benchmark view, union-only rows are treated as extras and the merge fails.                                                                          
     Fix: validate feathers against the run-owned union once, then project to each target view, or merge the union once and reuse it.                                                              
 19. Stage 4 reuse can annotate excluded/unselected datasets and foreign unions                                                                                                                    
     src/4_cell_type_annotation/1_submit_onboarding_stage.sh:130-152,590-614                                                                                                                       
     Reuse validation does not cross-check chunks.tsv against the current runnable selection or the three a-priori exclusions. It also only checks the chunk's union path for readability; it may  
     point to another run or a global union. The entire reused chunk manifest is then submitted.                                                                                                   
 20. Stage 4 sync-only bypasses ownership and exclusions                                                                                                                                           
     src/4_cell_type_annotation/1_submit_onboarding_stage.sh:201-211                                                                                                                               
     Sync-only derives its dataset list from all twelve selection rows, including the three excluded cohorts, and does not acquire or verify dataset ownership. An older terminal run can          
     overwrite a newer run's NAS output, and a successful exact run can fail while trying to sync unannotated excluded h5ads.                                                                      
 21. Annotation setup errors are converted into successful all-NA samples                                                                                                                          
     src/4_cell_type_annotation/2.1.1_process_chunk.R:522-531                                                                                                                                      
     The sample-wide tryCatch() catches setup and structural failures such as counts extraction, get_seurat_obj_from_h5ad(), normalization, and output materialization. Those errors are converted 
     into all-NA checkpoints. If other samples provide anchors, the dataset can pass final validation while silently losing a failed sample.                                                       
     Only method-level exceptions should use the nonfatal all-NA path.                                                                                                                             
 22. Annotation sidecar validation can OOM on full h5ads                                                                                                                                           
     src/utils/py/annotation_contract.py:44-55                                                                                                                                                     
     validate_sidecar() computes MD5 with artifact.read_bytes(), allocating the complete h5ad in memory before the backed validator opens it. This defeats the intended backed-mode behavior and   
     can exceed worker memory for full cohorts.                                                                                                                                                    
     Fix: hash in bounded chunks.                                                                                                                                                                  
 23. Nonnumeric annotation scores pass the strict contract                                                                                                                                         
     src/utils/py/annotation_contract.py:144-149                                                                                                                                                   
     Required score columns are checked only for existence. A value such as "not-a-number" passes if anchors and checksums are valid. Merge later coerces it to NaN, and reuse/skip validation can 
     bypass merge entirely.                                                                                                                                                                        
 24. Evidence wrapper can publish extra or unsafe artifacts                                                                                                                                        
     notebooks/dataset_onboarding/submit_batch_candidate_evidence.sh:52-53,119-132                                                                                                                 
   - RUN_ID is unsanitized and can contain ../, escaping both scratch and NAS roots.                                                                                                               
   - Wildcard loops include arbitrary extra *_batch_candidate_evidence.csv files.                                                                                                                  
   - Extra files are added to the checksum manifest and synchronized without dataset/schema validation.                                                                                            
   - mkdir/rsync failures under set -e can leave the evidence run ACTIVE with partial NAS output.                                                                                                  
     Validate a safe run-ID component, enumerate exactly twelve cohort files plus the review file, and record terminal failure for every synchronization error.                                    
 25. Durable profile audit masks malformed historical roots                                                                                                                                        
     .agents/skills/durable-hpc-gate-ecoda/references/profile.json:117-119                                                                                                                         
     The terminal run-root loop executes test -d for each root but does not accumulate failures. Bash returns the status of only the final iteration, so a malformed earlier run root can be       
     masked by a later valid one. Scope the check to the current run and aggregate all failures.                                                                                                   
                                                                                                                                                                                                   
 P2 findings                                                                                                                                                                                       
                                                                                                                                                                                                   
 - src/4_cell_type_annotation/2.1.1_process_chunk.R removed the previously commented breast_mode call arguments instead of retaining explicit commented rationale. The policy is documented, but   
   the requested source-level guard should remain adjacent to both upstream calls.                                                                                                                 
 - src/2_dataset_specific_preprocessing/1.3.1_prepare_joanito.R:70-76 strips terminal _Normal/_Cancer from raw cell.type, which can corrupt legitimate source labels.                              
 - src/2_dataset_specific_preprocessing/1.3.1_prepare_joanito.R:141-146 selects metadata from only the first row of each sample. Mixed per-cell sample.origin/Site values can pass pair validation 
   while the selected 500 cells violate the declared design.                                                                                                                                       
 - src/2_dataset_specific_preprocessing/1_submit_hpc.sh:107-132 now uses the canonical bassez_cellsubtype identifier, so explicit Bassez hook selection is preserved.
 - src/4_cell_type_annotation/2.0_create_scgate_db.R:16-23 validates only nonempty top-level payloads. Structurally invalid PBMC/HiTME caches can pass preflight and fail only after annotation    
   workers start.                                                                                                                                                                                  
 - docs/ARCHITECTURE.md:207-208 documents (Sample, barcode) while the implementation contract is (Sample, cell_barcode).                                                                           
 - NOTES.md:145-150 omits the deliberate “historical batch_effect_analysis.rmd unchanged/future extension” boundary and the deferred publication-figure handoff note.                              
                                                                                                                                                                                                   
 Test-quality blockers                                                                                                                                                                             
                                                                                                                                                                                                   
 Several new regressions do not actually protect the required contracts:                                                                                                                           
                                                                                                                                                                                                   
 - tests/test_benchmark_matrix_submitter.sh treats a malformed 1×2 Feather distance matrix as a valid reusable artifact.                                                                           
 - tests/test_benchmark_h5ad_contract.py derives its valid embedding fixtures from REQUIRED_OBSM, so removing a required key from the implementation constant can make the test self-consistent.   
   It also uses floating-point counts for the supposedly valid fixture.                                                                                                                            
 - tests/test_annotation_worker_policy.sh performs static grep and tests generic R coercion; it does not execute stubbed worker calls or verify the emitted NORMAL_TISSUE value.                   
 - tests/test_batch_candidate_evidence.py auto-creates the selection file, so it does not test the required selection-file failure contract.                                                       
 - tests/test_batch_candidate_evidence_wrapper.sh does not create post-job outputs or fake NAS synchronization, so missing-output, stale-output, and unvalidated-sync regressions pass.            
 - Exact-selection tests do not consistently assert no new run root/pending state, full middle-row order, wrong labels, or extra rows.                                                             
                                                                                                                                                                                                   
 No files were modified during this review.

 (end review findings by agent)

 Subsequent user prompt:
 " please write down your findings and proposed solutions, checks and validations in a new plan to fix and adress these findings. also, stage 2 of the initial plan                                  
 /Users/christianhalter/.omp/agent/sessions/-Desktop-ECODA_paper/2026-08-27T08-17-45-367Z_01a0424b-e597-7272-8e6a-23ffe1c37aab/local/onboarding-uncorrected-two-stage-plan.md remains an open todo 
 item; would you suggest to add it to the new plan or leave it in the original plan and just point out in the new plan (to implement the fixes) to see the original old plan and to finish stage 2 
 there?   "