# Alzheimer donor-by-assay processing and cross-cluster backup

Status: Draft. User-approved scope, backup order, and topology are recorded.
The approved global corrected composition/pseudobulk estimator is limma with
separate original technical covariates; lme4 is removed from active corrected
work. No implementation or HPC launch is authorized until the plan is
integrated, focused-tested, and the maintenance feasibility gate passes.

## Relationship to the current plan

The current batch-effect plan is
`.agents/plans/20260912-final-batch-effect-subset-plan.md`. Its superseding
scope amendment now owns the eight-dataset corrected-final Stage 5 recovery:
Joanito, Stephenson, Breast_cancer, Covid19_PBMC, Kidney_KPMP_full, Diabetes,
Lupus_PBMC, and Lung. This plan owns the excluded Alzheimer dataset and the
backup/transfer work. Historical gate failures and existing artifacts remain
read-only evidence.

The current plan's eight-dataset corrected gate must be snapshotted before
canonical `datasets.json` is changed for Alzheimer. Commit-keyed snapshots
make later canonical edits non-interfering, but the two corrected-final Stage
5 waves still share `batch_effect/corrected_final` and must serialize at the
shared Stage 5 synchronization boundary.

## Goals

1. Create an immutable, checksummed Alzheimer H5AD derivative whose
   standardized sample identity is donor-by-assay.
2. Change only the Alzheimer-dependent configuration and downstream rows after
   the eight-dataset corrected snapshot is sealed.
3. Run Alzheimer Stage 3 uncorrected and corrected preprocessing as explicit
   one-dataset view gates, concurrently only when durable ownership proves
   their output namespaces disjoint.
4. Run Alzheimer Stage 5 uncorrected and corrected-final methods with exact
   one-dataset selections. Corrected-final waits behind the eight-dataset
   corrected gate; uncorrected-final may overlap only after a disjoint-root
   gate review.
5. Resolve the corrected composition estimator before any composition worker
   is submitted.
6. Establish a tiny transfer proof of concept, then clone the repository and
   scratch data to an alternate HPC cluster, while retaining the existing
   explicit result/processed-data synchronization to NAS.

## Non-goals and invariants

- Do not mutate or delete the raw Alzheimer H5AD, existing processed H5ADs,
  RDS bundles, Feather files, checksums, manifests, gates, or logs.
- Do not rerun Pipeline 4 annotation for batch-effect views. Preserve the
  configured author/source cell-type columns.
- Biological labels, including `Cognitive status`, remain evaluation-only and
  never enter filtering, HVG selection, normalization, Harmony, batch
  correction, model covariates, or feature construction.
- Keep configured technical columns `assay` and `sex` as metadata and batch
  candidates. The new sample identity does not replace either technical
  column.
- Use only immutable source snapshots, run-bound source/runtime/auxiliary
  manifests, exact selections, atomic outputs, and checksum sidecars for full
  cohort work.
- Never use `$HOME/scratch` in an executor or backup command on Bamboo because
  it is a symlink. Resolve to `/srv/beegfs/scratch/users/h/halterc` first.
- No broad `--force`, historical matrix selection, or inferred dataset scope.

## Observed Alzheimer source contract

The authoritative raw input is currently:

```text
/srv/beegfs/scratch/users/h/halterc/ECODA_paper/Alzheimer/data/SEAAD_Alzheimer.h5ad
```

Metadata-only inspection found 1,395,601 cells, 83 donors, exactly two assay
values (`10x 3' v3` and `10x multiome`), and 104 unique donor-by-assay pairs.
The derived sample universe is therefore expected to contain 104 samples:
83 `10x3v3` samples and 21 `10xmultiome` samples. Sex is consistent within all
104 pairs and has 59 female and 45 male samples.

Existing configuration currently has `datasets.json:447`:

```json
"sample": "donor_id"
```

and `datasets.json:449`:

```json
"batch": ["assay", "sex"]
```

The existing corrected H5AD/metadata contract collapsed assay within a donor
under the old sample definition. Its Alzheimer row is not reusable after this
sample-universe change.

## Phase 0 — freeze the independent eight-dataset recovery

Before changing Alzheimer configuration:

1. Finish the required source-contract fixes, including corrected-final
   target/sync handling, effective-key propagation, cache-hit validation,
   canonical batch-value handling, and corrected-final regression coverage.
2. Resolve and implement the corrected composition estimator policy below.
3. Create a fresh immutable snapshot from the old configuration and an exact
   eight-row corrected-final selection:

   ```text
   Joanito<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   Stephenson<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   Breast_cancer<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   Covid19_PBMC<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   Kidney_KPMP_full<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   Diabetes<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   Lupus_PBMC<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   Lung<TAB>batch_effect_corrected<TAB>batch_effect_corrected
   ```

4. Target only:

   ```text
   prepare_pseudobulk,pseudobulk,gloscope,composition
   ```

The eight existing corrected prepare caches are checksum-valid but stale under
the approved global policy: they were generated with the prohibited
`__ecoda_batch_combined_v1` correction covariate. Preserve them as immutable
historical artifacts, but do not reuse them as new-policy inputs. The recovery
therefore emits 32 method rows: eight `prepare_pseudobulk`, eight
`pseudobulk`, eight GloScope, and eight composition rows. The validated
corrected MRVI/PILOT/QOT rows remain outside recomputation.
5. Run this gate only if a cutoff feasibility calculation proves that
   preflight, queue margin, workers/watchdogs, terminal audit, review, and
   synchronization finish before the maintenance boundary.
6. After the old-config snapshot is sealed and its gate is launched or
   otherwise terminally resolved, canonical Alzheimer source/configuration
   edits may proceed independently. The old gate cannot see those edits.

The current corrected-final selector accepts neither this exact eight-row
scope nor the later one-row Alzheimer scope. Both selector contracts must be
implemented and focused-tested before any corresponding gate.

## Phase 1 — Stage 2 donor-by-assay derivative

Integrate a new `alzheimer_donor_assay` step into the existing
`src/2_dataset_specific_preprocessing/1_submit_hpc.sh`. Do not create an
unguarded second submitter. The exact logical selector is:

```text
src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets Alzheimer \
  --steps alzheimer_donor_assay
```

Add a dedicated immutable-snapshot worker/hook following the existing Stage 2
runtime, ownership, watchdog, atomic-output, and checksum contracts. The hook
must:

1. read only the authoritative raw Alzheimer H5AD;
2. preserve the raw file unchanged;
3. require nonblank `donor_id` and `assay` values for every cell;
4. accept exactly `10x 3' v3` and `10x multiome`;
5. map them exactly to `10x3v3` and `10xmultiome`;
6. write `donor_id_assay = donor_id + "_" + assay_token`;
7. reject duplicate or colliding derived IDs and mixed metadata within a
   derived sample;
8. preserve `donor_id`, `assay`, `sex`, `Cognitive status`, source cell-type
   columns, and all other source observations;
9. write a new derivative atomically under the Alzheimer data area, with a
   run-owned artifact record and MD5/SIZE/PATH sidecar; and
10. validate the output schema, 104-sample expectation for the current source,
    source row count, exact assay vocabulary, and deterministic example IDs:

    ```text
    H20.33.001 + 10x 3' v3  -> H20.33.001_10x3v3
    H20.33.001 + 10x multiome -> H20.33.001_10xmultiome
    ```

Unexpected assay values, missing values, blank donor IDs, collisions, row
count changes, or inconsistent pair metadata are hard failures. Generic
string sanitization is not permitted.

Use a new explicit derivative filename rather than replacing
`SEAAD_Alzheimer.h5ad`. After the derivative is validated, update the
Alzheimer `file_names`/view input contract to point Stage 3 at that derivative
and set `columns.sample` to `donor_id_assay`. Keep `columns.batch` as
`["assay", "sex"]`; review the existing `majority_keys: ["assay"]` policy
against the now-constant assay-within-sample identity, but do not silently
remove it without a contract decision.

The Stage 2 run must be Alzheimer-only, snapshot-backed, and terminally
reviewed before Stage 3. If a valid derivative already exists, the run must be
`NOOP_VALIDATED` and submit no worker.

## Phase 2 — configuration and source identity

After Stage 2 derivative validation:

- update only the Alzheimer `datasets.json` sample/input contract and the
  corresponding `NOTES.md` explanation;
- preserve every other dataset/view filename and all Pixi constraints;
- record that the new sample universe is donor×assay, not donor-only;
- ensure the changed configuration invalidates only dependent Alzheimer
  Stage 3/Stage 5 rows;
- create a new commit-keyed snapshot containing the validated derivative
  contract and source code; and
- do not reuse pre-change Alzheimer H5AD, prepare, composition, MRVI, PILOT,
  or QOT outputs merely because their checksums pass. Their sample identity is
  stale under the new contract.

The updated H5AD source identity must record `Sample = donor_id_assay` and
retain `assay`/`sex` as technical metadata. All downstream validators must
compare the new 104-sample order exactly.

## Phase 3 — Stage 3 Alzheimer views

The existing Stage 3 submitter rejects a combined uncorrected-plus-corrected
selection. Use two explicit one-row selections:

```text
Alzheimer<TAB>batch_effect_uncorrected
Alzheimer<TAB>batch_effect_corrected
```

Each gate depends on the terminally reviewed Stage 2 derivative. The two
outputs have distinct configured filenames and may be submitted concurrently
only after:

- the durable profile accepts the selected gate serialization topology;
- each run has its own immutable source/runtime identity;
- ownership manifests prove no shared mutable output path;
- no shared metadata/export/checksum owner exists; and
- both exact selections are recorded before launch.

If any condition is not provable, serialize the two Stage 3 gates. Within each
view, the Alzheimer row is the only dataset row. Do not use the full configured
corrected selector, which would accidentally include the other eight cohorts.

Validate both outputs for:

- 104 unique ordered `Sample` IDs;
- preserved raw counts and required observation columns;
- `donor_id_assay` consistency with donor and assay;
- no biological-label covariate use;
- uncorrected semantic representation without Harmony; and
- corrected semantic representation using only configured technical keys.

## Phase 4 — Stage 5 Alzheimer methods

Add explicit one-row Stage 5 selectors for both lanes:

```text
Alzheimer<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Alzheimer<TAB>batch_effect_corrected<TAB>batch_effect_corrected
```

Use the established seven-method baseline:

```text
prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
```

Do not use broad/default selections, `--force`, or the historical exact matrix.
Existing pre-change Alzheimer rows are not reusable because the sample
identity changed; each new lane requires validator-only source/metadata
preflight and then targeted computation.

The uncorrected lane writes the separate `batch_effect/uncorrected_final`
root. The corrected lane writes the separate `batch_effect/corrected_final`
root. The corrected Alzheimer gate must wait behind the eight-dataset
corrected-final gate because both corrected lanes share the same corrected
root and synchronization owner. A second serialization group is not a safe
way around that boundary.

The uncorrected Alzheimer gate can overlap a corrected gate only if the two
roots, NAS destinations, checksum manifests, execution-time files, artifact
owners, and gate serialization locks are proven disjoint. Safety default:
serialize Stage 5 gates unless a focused durable-gate contract test and the
prelaunch manifest establish safe distinct-root concurrency.

Stage 5 corrected consumers must derive effective technical keys from the new
104-sample metadata. A configured key with one level is retained as metadata
but excluded from correction; if no key has at least two levels, correction is
`NO_CORRECTION`. The new Alzheimer source currently has two assay and two sex
levels, so both are expected to be effective after donor×assay grouping.

## Phase 5 — corrected composition estimator decision

### Evidence

`src/5_run_benchmark_methods/benchmark_methods_r.R:171-349` currently fits one
lme4 random-intercept model per CLR feature. A local experiment using the
existing 83-sample Alzheimer CLR composition bundle found:

```text
assay sample levels: 76 / 7
sex sample levels:   48 / 35
additive fits:       111 / 131 singular
zero variance:        99 / 131
lme4 convergence messages: 112 / 131
```

The current `correct_clr_batch_lmm()` fails on the first boundary fit with a
`boundary (singular) fit` message. The lme4 documentation says singular fits
are statistically defined but have higher numerical/inferential risk; it does
not define five or ten levels as a hard prohibition. The observed failure rate
is nevertheless decisive evidence that the current corrected composition path
must not launch unchanged.

Official references:

- <https://lme4.github.io/lme4/reference/isSingular.html>
- <https://lme4.github.io/lme4/reference/lmer.html>
- <https://lme4.github.io/lme4/reference/convergence.html>

### Approved global corrected limma policy

The user approved limma for **all corrected composition and pseudobulk
correction**, using separate original technical covariates. lme4 is removed
from the active corrected implementation. Existing artifacts remain immutable
and are not silently rewritten, but no new corrected output may claim or use
the lme4 random-intercept model.

For effective technical keys:

```text
design = model.matrix(~ 1 + batch_key_1 + batch_key_2 + ...)
```

The implementation must:

1. select only effective technical keys;
2. use fixed internal aliases, never configured names in executable formulas;
3. require full design rank and positive residual degrees of freedom;
4. fit pseudobulk counts with `DESeq2 design=~ 1`, without technical
   concatenation or biological labels;
5. apply limma correction to normalized pseudobulk values using each original
   technical covariate as a separate design term;
6. apply the analogous rank-checked limma fixed-effect correction to CLR
   composition values;
7. remove only non-intercept technical coefficient contributions;
8. preserve sample/gene or sample/cell-type identifiers and restore the CLR
   row zero-sum invariant exactly;
9. record a new explicit limma model identity such as
   `limma_separate_covariates_v1`, with exact correction mode and formula; and
10. retain configured one-level keys as metadata-only
    `non_estimable_batch_keys`; if all keys are constant, emit
    `NO_CORRECTION`, do not call limma, and leave values unchanged.

For Breast_cancer's three-key design, rank deficiency must fail closed rather
than silently dropping aliases. Perfectly confounded or zero-residual-degree
designs are errors. No `__ecoda_batch_combined_v1` value may be passed as a
correction model covariate; the original columns must remain visible in the
design and metadata.

The observed 83-sample Alzheimer bundle provides the decision evidence:
111/131 additive lme4 fits were singular, 99/131 had zero random-effect
variance, and 112/131 emitted convergence messages. The new donor×assay
sample universe still has only two assay and two sex levels, so retaining
lme4 is not an acceptable fallback.

The old lme4 model identity/formula is historical evidence only. Update
`batch_contract.R`, `benchmark_hpc_utils.R`, the composition/pseudobulk
workers, validators, and focused tests so new corrected artifacts cannot
silently mix old random-intercept or concatenated-batch contracts with the
approved limma contract.

Required focused regressions before the Stage 5 composition row:

- observed Alzheimer-like two-level assay/sex design;
- Breast-like three-key design;
- full-rank fixed design with finite output;
- rank-deficient design rejection;
- all-constant `NO_CORRECTION` identity result;
- exact coefficient removal/intercept preservation;
- CLR zero-sum preservation; and
- source inspection proving biological labels are absent.

## Phase 6 — backup and alternate-cluster transfer

### Documented facts

Local HPC documentation states:

- Bamboo, Baobab, and Yggdrasil have separate private networks, storage, and
  login nodes; jobs cannot be submitted across clusters.
- Home is backed up, but scratch is operational storage and is not backed up.
- NASAC/GVfs mounts are per-process/per-node and may require
  `dbus-launch bash` plus `gio mount`; network shares can be unreliable.
- documented cluster transfer uses `rsync -aviuzPrg` and supports `-n` dry runs.
- current local documentation and the newer retention announcement disagree
  on the exact scratch inactivity deletion threshold; confirm with HPC support.

Relevant local documents:

```text
docs/hpc_docs/storage_on_hpc.md
docs/hpc_docs/best_practices.md
docs/hpc_docs/access_the_hpc_clusters.md
docs/hpc_docs/hpc_clusters.md
docs/hpc_docs/data_life_cycle.md
```

Current Bamboo measurements:

```text
/home/users/h/halterc/ECODA_paper                         21 GB
/srv/beegfs/scratch/users/h/halterc/ECODA_paper          2.3 TB
/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper 92 TB free
```

Current noninteractive SSH from Bamboo to `baobab` and `yggdrasil` was not
available. This is an access fact, not proof that the clusters cannot reach
the NAS; each cluster needs a documented access/mount test.

### Backup order

1. **Inventory only.** Resolve canonical paths, cluster hostnames, destination
   capacity, quotas, file counts, symlink targets, and active writers. Do not
   recursively copy an active tree.
2. **Tiny transfer proof of concept.** After access is confirmed, transfer a
   small explicitly created test file from Bamboo to each intended alternate
   cluster and, separately, to the approved NAS destination. Verify source and
   destination SHA-256/MD5, permissions, and a second dry-run with no changes.
   The test must not use a project output path.
3. **Repository clone.** Transfer the complete 21 GB Bamboo repository clone
   first, including `.git`, source, plans, runtime metadata references, and
   hidden files, to a timestamped alternate-cluster destination. Use
   resumable, non-destructive rsync and an explicit sorted path/size/hash
   manifest. Do not use `--delete` or `--inplace` for the first mirror.
4. **Scratch clone.** After a quiescence checkpoint, mirror the complete
   scratch `ECODA_paper` tree to the alternate cluster in explicit top-level
   batches. Preserve raw/processed data, results, logs, gates, manifests, and
   checksums there. Use resumable transfer, per-batch logs, source/destination
   hash manifests, and a second dry run. Do not rely on scratch retention.
5. **NAS synchronization.** Keep the existing Stage 5 behavior: sync only
   explicitly selected validated results and processed-data artifacts to NAS,
   with exact checksums and sidecars. The repository clone, logs, gates, and
   full scratch control plane do not need a second NAS copy if the alternate
   cluster clone is verified.
6. **Restore/review.** Sample-restore small and large files from the alternate
   cluster, compare manifests, verify runtime/source identities, and obtain
   human review before considering the clone complete.

The full clone must be quiescent enough that manifests are not changing during
hash verification. Active Stage 5/Stage 3 writers must never be copied and
then treated as a consistent completed backup.

## Parallel execution matrix

| Work | Can overlap | Required boundary |
|---|---|---|
| Focused source/tests, LME design, and backup documentation | Each other | No HPC or artifact mutation |
| Tiny transfer PoC and local source work | Yes | Separate non-production paths |
| Eight-dataset corrected Stage 5 and Alzheimer Stage 2 | Yes, after old-config snapshot is sealed | Distinct source snapshots and output/owner roots |
| Alzheimer Stage 3 uncorrected and corrected | Yes | Stage 2 terminal review plus disjoint output ownership; otherwise serialize |
| Eight-dataset corrected Stage 5 and Alzheimer corrected Stage 5 | No | Same `batch_effect/corrected_final` root and shared sync owner |
| Alzheimer uncorrected Stage 5 and corrected Stage 5 | Only with explicit disjoint-root proof | Separate roots, checksum manifests, NAS paths, and allowed serialization groups |
| Full scratch backup and active writers | No | Quiescence checkpoint and stable manifests |

Every full-cohort gate remains snapshot-backed and follows exactly one durable
wait, one terminal inspect over every emitted scheduler/watchdog ID, and the
required review before dependent work.

## Maintenance feasibility gate

The maintenance notice's absolute boundary is:

```text
2026-09-15 08:00 +0100
2026-09-15 07:00 UTC
2026-09-15 09:00 Bamboo local time (UTC+0200)
```

At the last check Bamboo was `2026-09-14 12:22 UTC`. A launch decision must
include at least the requested worker/watchdog walltime plus preparation,
source/preflight, queue-start margin, terminal accounting, reviewer, and
synchronization time. The planning default is `8h` submitted walltime plus
`2h` operational margin; a 12-hour default is allowed only if the computed
launch time and all dependencies still end before the absolute cutoff with
additional buffer. If not provable, do not launch before maintenance.

## Verification and release criteria

- Stage 2 derivative test proves exact mapping, 104 samples, no collisions,
  raw preservation, atomic checksum, and invalid-input failures.
- Both Stage 3 one-row outputs validate the new sample identity and exact
  configured technical metadata.
- Corrected composition tests pass the selected estimator policy before any
  composition worker is submitted.
- Stage 5 one-row Alzheimer selectors and eight-row non-Alzheimer selector
  reject broad/legacy scopes and preserve validated rows outside recompute.
- Every new gate has immutable source/runtime/auxiliary identity and exact
  run-owned selections.
- Alternate-cluster clone has matching path/size/hash manifests and a tested
  restore; NAS receives only explicit validated result/processed-data lists.
- No existing artifact is overwritten, invalidated, or deleted.
