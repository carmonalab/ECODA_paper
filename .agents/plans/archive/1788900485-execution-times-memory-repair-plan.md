# Local execution-time memory repair

## Goal
Fill missing `mem_GB` values for existing benchmark execution-log rows using validated local result bundles. Do not run HPC jobs, workers, submitters, pipelines, SSH, or schedulers.

## Provenance boundary
- Treat `data/benchmark/embeddings/execution_times_OLD.feather` as the authoritative baseline: MD5 `4df267b7492cce9937db26a717d3aa24`, size 26994 bytes, 1376 rows; it matches `data/benchmark/checksums.md5` and the reconciliation manifest.
- Preserve the current `execution_times.feather` unchanged as a candidate backup before installing any repaired artifact. It is MD5 `3d9cdb84633e02e52606feb3b0ffccf4`, size 27538 bytes, 1403 rows, and contains 27 additional PILOT-GM-VAE rows that the reconciliation manifest marks as not installed/blocked.
- Do not silently accept or erase those 27 candidate rows during the memory repair.

## Repair
- Copy the current candidate Feather to a uniquely named local pre-repair backup.
- Start the repaired table from the authoritative OLD baseline, preserving every baseline key, runtime, and non-NA memory value.
- For existing baseline `(dataset, method)` keys only, fill `mem_GB` when it is NA and the matching validated local RDS bundle contains one finite scalar `mem_GB`. Never append duplicate keys, change `time_secs`, fabricate values, or use zero as a substitute.
- Write the repaired Feather atomically and record a JSON repair-evidence file containing source/output hashes, changed keys, and validation results.
- Update only the aggregate checksum entry for `embeddings/execution_times.feather` after successful validation; do not overwrite the reconciliation manifest or claim the unapproved candidate was installed.

## Verification
- Validate exact four-column schema, nonblank identifiers, unique `(dataset, method)` keys, unchanged baseline times/keys, finite-or-NA memory values, and expected MOFA/GloScope default coverage.
- Confirm the repaired default MOFA/GloScope rows have no recoverable `mem_GB` gaps for production datasets.
- Regenerate Supp. Fig. 14B locally, without executing any pipeline code, and verify the output exists and the source data contains the repaired points.
