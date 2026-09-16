# ECODA final batch-effect subset processing

## Current scope — authoritative (2026-09-15)

This section is the single active contract for this plan. It supersedes every
older active Context, Approach, Verification, Assumptions, and status section
that conflicts with it. The incident appendix is evidence only; it does not
create selectors, authorize a launch, or preserve an obsolete nine-row scope.
This plan now records the user's latest execution override: run only two
independent Stage 5 jobs on Yggdrasil. The corrected job covers all nine
configured batch-effect datasets and all seven corrected methods; the
Alzheimer uncorrected job covers its one dataset and the same seven methods.
No Stage 2 or Stage 3 review, acceptance gate, rerun, or dependency is part
of this execution. Every Stage 5 job still uses the required immutable
snapshot/runtime identity, durable terminal accounting, artifact audit,
synchronization, and Luna Max review.
The former 35-row matrix and separate corrected-Alzheimer follow-up are
superseded for this execution and remain historical evidence only.
### Global dataset-level column authority — current decision

The production correction contract resolves metadata columns only from each
dataset's top-level `columns` object. `views.*.columns` is prohibited: views
may select input/output/subset behavior, but they cannot replace, narrow, or
add correction columns. Stage 5 cannot select a separate correction key, and
no historical cohort-specific override is active.

For `Breast_cancer`, `columns.batch` is exactly
`["assay", "suspension_dissociation_time"]`. Corrected Stage 5 consumes this
dataset-level pair. `sequencing_platform` remains in H5AD `obs` and exported
sample metadata when present, but it is not a Breast correction column.
`disease` remains a biological label and never enters correction.

The previous Breast corrected Stage 3 contract and output history remain
immutable evidence. Under the user's override, no Stage 3 review, acceptance
gate, or rerun is performed. Required source, H5AD, metadata, checksum, and
artifact checks remain inside the two Stage 5 jobs.

### Confirmed two-job Stage 5 execution scope — 2026-09-15

The user explicitly overrides the earlier staged ordering. The only new
processing is two independent Stage 5 jobs on Yggdrasil:

1. One corrected batch-effect Stage 5 job covering all nine configured
   datasets in this exact order:
   `Joanito`, `Stephenson`, `Alzheimer`, `Breast_cancer`, `Covid19_PBMC`,
   `Kidney_KPMP_full`, `Diabetes`, `Lupus_PBMC`, and `Lung`.
2. One Alzheimer uncorrected batch-effect Stage 5 job covering `Alzheimer`
   only.

Both jobs use the complete seven-method suite:
`prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot`.
The corrected job therefore declares 63 logical rows (9 datasets × 7
methods), including the full corrected Alzheimer Stage 5 workload. The
Alzheimer uncorrected job declares 7 logical rows.

There is no separate corrected Alzheimer Stage 5 gate. The two jobs are
launched independently and in parallel. The corrected job uses
`batch_effect/corrected_final`; the Alzheimer uncorrected job uses
`batch_effect/uncorrected_final`.

All historical Stage 2, Stage 3, and corrected Stage 5 gate records are
evidence only. The user explicitly declines Stage 2/3 review or acceptance
work; do not create, review, rerun, or use those records as prerequisites.
Existing valid Stage 5 artifacts remain eligible for ordinary idempotent
skipping; no broad force rerun is implied.

This two-job scope supersedes the former 35-row per-dataset method matrix,
the separate corrected Alzheimer follow-up, the obsolete eight-dataset/32-row
task, and any uncorrected changed-dataset/Kidney-only wave.

The Stage 5 durable workflow remains mandatory for both jobs: immutable
Yggdrasil source snapshot/runtime binding, explicit selection, terminal
accounting, artifact audit, synchronization, and Luna Max review.

This section is the active authorization and ordering contract; later
historical sections cannot expand or replace it.

### Historical implementation and gate evidence

All earlier Stage 2/3 acceptance plans, the former 35-row matrix, the
separate corrected Alzheimer follow-up, the obsolete eight-dataset/32-row
task, and the uncorrected changed-dataset/Kidney-only wave are historical
evidence only. Their files and artifacts remain immutable and do not
authorize new work.

## Release checklist

The active release scope is exactly two Stage 5 jobs:

- corrected batch effect: all nine datasets
  (`Joanito`, `Stephenson`, `Alzheimer`, `Breast_cancer`, `Covid19_PBMC`,
  `Kidney_KPMP_full`, `Diabetes`, `Lupus_PBMC`, `Lung`) × the seven methods,
  for 63 logical method rows;
- Alzheimer uncorrected batch effect: `Alzheimer` × the seven methods,
  for 7 logical method rows.

The jobs use explicit selection files, separate Yggdrasil source snapshots
when run concurrently, and their distinct output roots. Existing valid
artifacts may be skipped individually. No Stage 2 or Stage 3 review,
acceptance, rerun, or dependency is performed.

Each Stage 5 job still requires its own durable wait, one terminal accounting
and artifact inspect, Luna Max review, and manifest-listed synchronization.
The corrected job has no separate Alzheimer follow-up: corrected Alzheimer is
already included in its nine-dataset scope.

### Open items and exact next steps

1. Complete the currently launched corrected all-nine Stage 5 job.
2. Launch the Alzheimer uncorrected Stage 5 job independently and in parallel
   using an isolated Yggdrasil source snapshot.
3. After both jobs reach terminal review, synchronize only their manifest-listed
   Stage 5 artifacts and run the final analysis.
### Alzheimer Stage 3 retry-acceptance decision

The old Alzheimer Stage 3 gates remain `FAILED` and unreleased because their
complete recorded attempt chains include the initial and retry1
`OUT_OF_MEMORY|0:125` attempts. The retry2 H5AD outputs remain preserved and
immutable; they are evidence for recovery, not permission to rewrite the old
run.

A fresh recovery gate may be validator-only and may accept only explicitly
validated retry2 array/watchdog outcomes. It must preserve the full initial,
retry1, and retry2 attempt chain in its new run-owned acceptance report,
including the OOM states. It must not manufacture scheduler IDs for
superseded attempts, submit preprocessing, mutate the failed Stage 3 roots or
manifests, or create an artifact record for an old H5AD. No preprocessing
rerun is allowed merely to erase OOM history. Alzheimer Stage 5 may depend
only on fresh recovery records that pass terminal audit and required review;
the failed original gates and unreviewed retry2 outputs are not sufficient.



### Objective and hard boundaries

- Local selector/source-contract implementation, backup feasibility, and the
  completed transfer-sanity mirror are recorded. The scratch mirror is now
  the active Yggdrasil working tree and the repository is at its canonical
  Yggdrasil home path. Remaining work is scoped execution, not Bamboo-first
  migration.
- The Alzheimer Stage 2 derivative is complete on Bamboo under its exact
  snapshot/run contract. Its validated derivative is an input for the
  Yggdrasil follow-up; it must not be recomputed or modified outside an
  explicitly approved recovery.
- The production source is the current `datasets.json` plus authoritative HPC
  data. Local mirrors and historical JSON reports are diagnostic evidence only.
- Existing H5ADs, RDS bundles, pseudobulk caches, Feather files, manifests,
  checksums, logs, and plots are immutable. Reuse is validator-only and
  artifact-by-artifact; no broad `--force`, historical matrix, inferred scope,
  or overwrite is allowed.
- Batch-effect views do not invoke Pipeline 4 annotation. Preserve configured
  source/author cell-type columns and keep biological labels evaluation-only.
- Every future full-cohort operation uses the checked-in
  `durable-hpc-gate-ecoda` workflow, an exact run-owned selection, immutable
  source/runtime/auxiliary identities, atomic outputs, checksums, one
  unbounded durable wait, one terminal inspection over every emitted ID, and
  the required Luna Max review.
### Yggdrasil host, runtime, and storage contract

Yggdrasil is the default compute, data, results, and backup host for this plan
until the user explicitly directs a return to Bamboo. Bamboo is fallback/source
infrastructure only; do not infer a host change from maintenance dates.
All authoring and pipeline-file changes happen on the local workstation, then
are committed and pushed; Yggdrasil pulls the exact committed revision. Never
edit the Yggdrasil checkout directly.

Operationally, `PORTABILITY_AUDIT=READY_FOR_SCOPED_LANES`: every new gate uses
an exact local commit, immutable source snapshot, runtime identity, explicit
selection, and durable manifest. The active paths are
`~/ECODA_paper` and `~/scratch/ECODA_paper`. The mirrored pinned environment
provides Python `3.13.14` and R `4.5.2`; the scoped batch-effect lanes use CPU
resources. The explicit NAS target is
`~/scratch/ECODA_paper/_nas_mirror/Projects/ECODA_paper` because NASAC is not
mounted.

The durable profile targets `remote_host=yggdrasil`, and its root checks use
the explicit NAS mirror. There is no implicit Bamboo NAS fallback. Before any
new pipeline gate, validate the source snapshot, runtime, auxiliary root,
selection, and configured roots. No broad selection or partial method subset
is permitted.
### Phase order (superseded)

The active phase order is the confirmed execution scope near the top of this
plan. This earlier phase-order narrative is retained only as historical
context and must not derive a selector, method matrix, or scheduler launch.
The confirmed sequence is review-only prerequisites, one explicit 35-row
corrected-final gate, Alzheimer uncorrected Stage 5 on its disjoint root, and
Alzheimer corrected Stage 5 after the shared corrected-final owner is reviewed.
### Explicit execution ordering and parallelism (superseded)

The active execution and parallelism contract is the confirmed section near
the top of this plan. The earlier ordering is retained only as historical
context. In particular, the former eight-dataset/32-row gate is not active;
Breast contributes seven corrected methods, and the non-Breast scope
contributes 28 rows to the single 35-row corrected-final matrix.
### Eight-dataset/32-row corrected-final Stage 5 recovery (superseded)

This former recovery is historical evidence only. Its failed gate,
selection, manifests, and partial artifacts remain immutable and cannot
authorize reuse or launch. The active replacement is the explicit 35-row
per-dataset method matrix documented near the top of this plan.
## Historical release checklist (superseded)

The active release checklist is the confirmed checklist near the top of this
plan. This older checklist is retained only as historical context and must
not authorize the former eight-dataset/32-row task.
### Open items and exact next steps (superseded)

The current status and next steps are recorded in the active confirmed-scope
section near the top of this plan. This older checkpoint is retained only as
historical evidence; it must not be used to launch the superseded 8×32 task.
### Final corrected-method policy

This is a final policy, not an experiment or an alternative:

- **Dataset-level authority is mandatory.** Corrected Stage 3/Harmony and
  corrected Stage 5 composition/pseudobulk resolve correction columns only
  from each dataset's top-level `columns.batch`. Views cannot replace, narrow,
  or add those fields, and Stage 5 cannot select a separate correction key.
  No historical cohort-specific override can change that source of truth.
- **Every corrected composition and corrected pseudobulk mode uses limma
  fixed effects with the configured dataset-level technical covariates as
  separate design columns.** For effective keys, use a fixed internal-alias
  design such as `model.matrix(~ 1 + technical_key_1 + technical_key_2 + ...)`;
  any effective/non-estimable report is derived from those dataset-level
  fields and is not an independent configuration.
- For `Breast_cancer`, `columns.batch` is exactly
  `["assay", "suspension_dissociation_time"]`. `sequencing_platform` remains
  visible metadata in H5AD `obs` and exported sample metadata when present,
  but it is not a configured correction column. The biological `disease`
  label is never a covariate.
- Never construct a combined/artificial batch key such as an interaction or
  concatenated `batch_key` in place of the separate technical columns. Keep
  configured technical fields separately visible in metadata, source identity,
  manifests, and the exact recorded design string.
- Biological labels and sample IDs are never model covariates. A one-level
  dataset-level technical field is retained as metadata but omitted from the
  effective design; if no technical field varies, record `NO_CORRECTION` and
  return the uncorrected object for that method.
- Fail closed on missing metadata, non-finite output, rank deficiency, or
  non-positive residual degrees of freedom. Do not silently drop an aliased
  configured technical column. Remove only the non-intercept technical
  contribution, preserve the intercept, and restore the CLR row-sum invariant
  for composition.
- Record a stable identity such as `ecoda_additive_fixed_effects_v1`, the
  exact formula/design, effective and non-estimable keys, and correction state
  in both composition and pseudobulk results. The limma pseudobulk operation
  remains method-specific; sharing the batch-only fixed-effect family does not
  collapse it into the composition feature method.
- `lme4` random-intercept fitting is removed from all new corrected work. No
  corrected worker, validator, selector, or recovery may call or advertise it
  as a supported path. Existing lme4 payloads are immutable historical files
  only and are not corrected-result reuse candidates under this policy.


### Strict Alzheimer donor-by-assay follow-up

The current raw input evidence is
`/srv/beegfs/scratch/users/h/halterc/ECODA_paper/Alzheimer/data/SEAAD_Alzheimer.h5ad`.
It remains unchanged. The new Stage 2 logical selector is:

```text
src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets Alzheimer \
  --steps alzheimer_donor_assay
```

The dedicated snapshot-backed hook writes a new derivative rather than
replacing `SEAAD_Alzheimer.h5ad`. It must read only the raw input, require
nonblank `donor_id` and `assay` for every cell, accept exactly these assay
values, and apply only these mappings:

```text
10x 3' v3       -> 10x3v3
10x multiome    -> 10xmultiome
```

It writes `donor_id_assay = donor_id + "_" + assay_token`, rejects any
unexpected/missing/blank assay, donor collision, duplicate derived ID, mixed
metadata within a derived sample, or row-count change, and preserves
`donor_id`, `assay`, `sex`, `Cognitive status`, source cell types, and all
other observations. It validates atomically with an artifact record,
checksum, size, and deterministic example IDs. A valid derivative is
`NOOP_VALIDATED`; it is never rebuilt merely because a later gate is created.
After derivative validation, the Alzheimer input contract uses
`columns.sample = donor_id_assay` while retaining `columns.batch = ["assay",
"sex"]` and the configured biological/cell-type metadata.

Metadata-only inspection observed:

- 1,395,601 cells;
- 83 donors;
- 104 unique donor-by-assay samples;
- assay sample counts 83 (`10x 3' v3`) and 21 (`10x multiome`);
- sex sample counts 59 and 45; and
- no donor-by-assay sample has mixed sex.

These values are the acceptance contract for the derivative and downstream
sample order, subject to a fresh source-bound report. Pre-derivative
Alzheimer H5AD, metadata, prepare, composition, pseudobulk, MRVI, PILOT, and
QOT outputs are not reusable merely because their checksums pass: their
sample identity is donor-only.

#### Alzheimer Stage 3

After the reviewed Stage 2 derivative, use two separate one-row,
view-specific selection manifests:

```text
Alzheimer<TAB>batch_effect_uncorrected
Alzheimer<TAB>batch_effect_corrected
```

Each gate must export
`STAGE3_INPUT_PRODUCER_RUN_ID=stage2_alzheimer_donor_assay_20260914T172421Z`.
The submitter materializes and validates its run-owned `input_ownership.tsv`,
binding `Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad` to that validated
Stage 2 producer before selecting the view. The current Stage 3 submitter
rejects a combined Alzheimer selection, so the two one-row gates use separate
serialization groups and separate run roots, output paths, source/runtime
identities, ownership/checksum records, and snapshot parents; these groups
represent disjoint view owners rather than a bypass of a shared owner. Both
outputs must contain the ordered 104 `donor_id_assay` samples, preserve
original technical metadata and raw counts, exclude biological labels from
processing covariates, and use semantic uncorrected/corrected representations
respectively.
#### Alzheimer Stage 5

Use explicit one-row lanes, never a broad/default selection:

```text
Alzheimer<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Alzheimer<TAB>batch_effect_corrected<TAB>batch_effect_corrected
```

The approved one-row baseline method list is:

```text
prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
```
These Alzheimer Stage 5 lanes are downstream of the two Alzheimer Stage 3
gates and are not part of the selected-eight 32-row recovery. They use
explicit one-row lanes and never a broad/default selection.

Valid rows may be skipped only after validator-only source/metadata,
checksum, ownership, and model-contract checks. Pre-derivative Alzheimer
rows are not valid reuse candidates. The uncorrected lane uses
`batch_effect/uncorrected_final`; the corrected lane uses the physical
replacement root `batch_effect/corrected_final/recovery_35row` and the final
limma policy above.

The required Stage 3 predecessor, same-root serialization, and allowed
parallelism rules are defined once in **Explicit execution ordering and
parallelism** above. Do not create a second corrected-final owner.

### Snapshots, locks, and concurrency

- Before each gate, record the exact snapshot commit, runtime identity,
  auxiliary manifests, selector, expected rows, output roots, and dependency
  reason. A running gate sees only its immutable snapshot; later canonical
  configuration edits cannot change it.
- Snapshot executor locks are **parent-wide**: an `.ecoda-exec-lock` held by a
  snapshot parent covers every child snapshot/executor operation below that
  parent. Do not start a concurrent child, launch from a partially created
  parent, remove a lock blindly, or evade it with another serialization-group
  name. Reconcile the owner and terminal state first. Safe parallel work must
  use deliberately disjoint canonical snapshot parents and disjoint run/root
  ownership.
- Corrected-final Stage 5 gates that target the physical replacement root
  `batch_effect/corrected_final/recovery_35row` share the Stage 5
  synchronization owner, checksum merge, execution-time merge, and NAS
  destination. They **must** use the `ecoda-benchmark` policy group and
  serialize. The confirmed 35-row matrix gate (Breast seven methods plus
  seven non-Breast four-method rows) must reach terminal audit/review before
  the one-row Alzheimer corrected gate begins. A different durable group cannot
  make same-root synchronization safe.
- The temporary large-mirror verification is intentionally `TRANSFER_SANITY_ONLY`.
  A full `rsync --checksum` scan is not required: it rereads terabytes and
  hundreds of thousands of files, and the repository scan already exceeded
  the practical time budget. Record `CONTENT_CHECKSUM=DEFERRED`; do not call
  the mirror cryptographically verified.
- The required completion evidence is a quiescent source, rsync exit `0`, no
  failure marker, a destination root, and the transfer success marker.
  A rough `du -sh` or coarse file-count sanity check is optional when cheap;
  exact source/destination equality, per-file hashes, permissions, and a
  second checksum dry run are not required for this explicitly approved
  temporary backup. Temporary rsync metadata/partial directories may differ
  while a transfer is active and are cleaned or recorded before final status.

### Backup and alternate-cluster priority

The local operational references are `AGENTS.md` (durable gate, snapshot,
selection, and same-root rules), `README.md` (repository onboarding),
`NOTES.md` (metadata audit and policy evidence), and the HPC knowledge-base
snapshots:

```text
docs/hpc_docs/storage_on_hpc.md
docs/hpc_docs/best_practices.md
docs/hpc_docs/access_the_hpc_clusters.md
docs/hpc_docs/hpc_clusters.md
docs/hpc_docs/data_life_cycle.md
```

The structural cleanup follow-up remains separate at
`.agents/plans/1789248891506-ecoda-pipeline-structure-plan.md`; it is not an
authorization to refactor this run.

Recorded Bamboo measurements are:

```text
/home/users/h/halterc/ECODA_paper                                      21 GB
/srv/beegfs/scratch/users/h/halterc/ECODA_paper                       2.3 TB
/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper             92 TB free
```


The complete Bamboo→Yggdrasil scratch transfer passed the approved
`TRANSFER_SANITY_ONLY` check and was explicitly reclassified by the user as
the active writable `~/scratch/ECODA_paper` working tree. It is no longer an
independent immutable backup.

The repository mirror move also completed. The canonical Yggdrasil repository
is `~/ECODA_paper`; the active data/results tree is `~/scratch/ECODA_paper`.

Large mirrors use bounded verification only: a quiescent source, rsync exit
`0`, no failure marker, destination presence, and a success marker. Rough
size or file-count checks are optional. Full content checksum scans and exact
source/destination equality are intentionally deferred and recorded as
`CONTENT_CHECKSUM=DEFERRED`.

The Bamboo NAS mount
`/srv/smednas515.unige.ch/carmona_smb` is not present on Yggdrasil. Yggdrasil
resolves `nasac-evs2.unige.ch` and has `gio`/D-Bus clients; a user-scoped
interactive NASAC mount is required before result/NAS synchronization. No
password or private key may be stored.

The migration sequence completed: the repository move, exact local
commit/push and Yggdrasil revision pull, canonical source/runtime/auxiliary
checks, scheduler/profile/path checks, scratch activation, and explicit NAS
mirror routing were verified. The remaining operations are only the exact
approved selections below, each through a new snapshot-backed durable gate.

### Maintenance feasibility gate — superseded for Yggdrasil

The prior Bamboo maintenance boundary of `2026-09-15 07:00 UTC` remains
historical evidence only. It no longer constrains the active run because the
user explicitly selected Yggdrasil as the default host and approved controlled
Wave 1 execution there.

Every new Yggdrasil launch still requires validator-only preflight, an exact
immutable source/runtime/auxiliary identity, explicit selected rows and
expected counts, no active serialization conflict, and queue feasibility for
the declared CPU resources. Record those checks in the run manifest before
launch. This supersession does not authorize a broad selection, a stale
manifest reuse, or a bypass of terminal accounting, artifact audit, or Luna
Max review.
## Historical evidence appendix — non-authoritative

This compact record preserves only major implementation evidence, terminal
gate outcomes, and validator reports. It is not an active scope and does not
authorize reuse or relaunch. The old nine-row corrected launch is historical;
its Alzheimer row used donor-only samples and is superseded by the
donor-by-assay contract. The former eight-dataset corrected-final gate had a
32-row scope; its old combined-key prepare caches and all lme4 payloads remain
immutable historical artifacts, not reuse candidates for the confirmed
35-row replacement.

### Major commits retained as evidence

| Commit | Evidence retained |
|---|---|
| `e9ee50add76c2e7826980d7333e6f9440d5c647b` | Initial subset/preflight, final-variant, metadata-export, and manifest implementation. |
| `986c6c7` | Stabilization wave and Stage 3 regression before later trust-boundary fixes. |
| `332f7c4` | Immutable-source/RDS metadata auditing and focused regression. |
| `b8f7aec15dfc91bc493caeb04f3d2cc066db8c5d` | Lightweight subset and source-bound H5AD audit; verified source snapshot retained. |
| `734b174a0b0b2a9c4e07edbf1e11d03c9fbf8206` | Container source-root bootstrap repair and verified snapshot. |
| `59b781fa31e6d9fb015e1d7911aa283f5131ea6c` | Reviewed Stage 3/Stage 5 runtime and ownership integration; full source snapshot. |
| `837eafb80c5b204ba6999f3c51ec8e4d95b09ddb` | Historical targeted selector; its hvg2000 dependency-reuse validation is superseded for the eight prohibited combined-key prepare caches. |
| `001243f351037590dc8df7db22cb2d34561e1452` | Historical producer-bound pseudobulk-cache validation; not evidence for reuse of the combined-key prepare caches. |
| `18906147debaf88d7c439642ea032212663e9868` | High-resolution cell-type/exporter contract and corrected-final snapshot. |
| `ba8b9dad556cd44e6cfffbb045a2b6180098d757` | H5AD preflight terminal-failure publication and long-checksum grace. |
| `9147dc0f2b0285bf2dc201f388153fb0962f9356` | Corrected preflight array-accounting wait. |
| `7f5604d5891d8f4ef8dcae0b457a2a4ad2162b4` | Positive preflight cardinality and zero-row fail-closed guard. |
| `b1f127e` | Earlier R corrected-final H5AD summary-policy propagation. |
| `6bf9307` | Explicit Python corrected-final consumer validation context. |
| `0d07b1d67e694c638e247ecb608dd2abf453d466` | Historical corrected-final selector/sync-only acceptance; superseded by the mandatory 32-row recovery and final limma identity. |

### Failed or superseded gate evidence

| Gate / IDs | Terminal reason and disposition |
|---|---|
| `stage2_joanito_final_20260912` | Snapshot executor rejected textual `$HOME/scratch` symlink before submission; no scheduler IDs. |
| `stage3_batch_final_20260912b`, `stage3_batch_final_20260912c` / preflight `4403668` | Slurm-spool/source-root bootstrap failure, then Covid obs-only preflight failure; no Stage 3 rows released. |
| `stage3_uncorrected_final_20260912` / `4403790`, with extra `4403791` | Covid preflight resolved the worker's relative source path incorrectly; `4403791` later settled independently. Replacement `...20260912b` completed IDs `4403794/4403795` and was reviewed; its NAS owner discrepancy was later reconciled validator-only. |
| `stage3_corrected_final_20260913a` / `4403877`, `4403887`, `4403903`, `4403989`, `4403888` | Old nine-row corrected run had an OOM retry and aggregate failure; later sync repair was validator-only. Its nine-row scope is historical, not current. |
| `stage5_uncorrected_final_20260913a` | Held by the shared source-snapshot parent lock before scheduler submission; no scheduler/artifact. |
| `stage5_uncorrected_final_20260913b` / aggregate `4403988` | Scheduler rows for GloScope, composition, and MRVI launched and failed; preparation/pseudobulk/PILOT/QOT evidence was retained and valid Kidney rows were skipped individually. The aggregate gate then failed on those row outcomes; this was not a prelaunch dependency block. This was an uncorrected historical lane and does not qualify corrected reuse. |
| `stage5_uncorrected_composition_recovery_20260913a` / `4404199`, `4404204`, `4404406`, `4404408`, `4404410`, `4404412` | Scope mismatch emitted GloScope + composition + MRVI instead of composition-only; arrays were canceled and the gate failed closed. |
| `stage5_corrected_final_20260913b` | Shared corrected-final snapshot parent lock held; no scheduler submission. |
| `stage5_corrected_final_20260913c`, `...d` | Parser contract failure, then Lupus metadata export requested an unavailable low-resolution column; no corrected method result was accepted. |
| `stage5_corrected_final_20260913e` / metadata `4406405`, preflight `4406437` | Long Alzheimer H5AD checksum/preflight status was not available in the grace window; no method array was emitted. |
| `stage5_corrected_final_20260913f`, `...g`, `...h` | Wrong serialization group, stale parent lock, and zero preflight-row counter respectively; all stopped before compute with no reusable gate. |
| `stage5_corrected_final_20260913i` | Scheduler rows launched; Alzheimer prepare and GloScope rows failed, while MRVI/PILOT/QOT completed. Completion transport could not be recovered, so the durable gate stopped in `PRELAUNCH_STOP`; this is not evidence that the launched rows were dependency-blocked. Its Alzheimer row is invalid under donor-by-assay scope. |

| `stage3_alzheimer_wave_ygg_20260915T115500Z`, `stage5_eight_corrected_wave_ygg_20260915T115439Z` | First Yggdrasil Wave 1 attempt failed before scheduler submission: Stage 3 supplied a combined Alzheimer selection, and Stage 5 supplied a four-method option reserved by the submitter for the fixed seven-method suite. No scheduler IDs or reusable artifacts. |
| `stage5_eight_corrected_ygg_retry_20260915T115439Z` | Corrected four-method selector reached source-H5AD sidecar validation and failed because the explicit NAS mirror was empty; no scheduler IDs. Destination-only normalization then completed validator-only for all eight H5ADs. Durable report: `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_logs/20260915T123309Z_stage5_inputs/SIDECAR_VALIDATION.tsv`; marker `SIDECAR_SUCCESS`; report `ROWS=8`, `COMPARISON=SOURCE_DESTINATION_MD5_SIZE_PATH_PASS`. The failed manifest remains unreused. |
| `stage5_eight_corrected_ygg_parallel_20260915T131349Z` / preflights `45687614`, `45687624` | Corrected-final consumer barrier passed seven dataset rows but failed `Breast_cancer` on the rank-deficient `assay`/`sequencing_platform`/`suspension_dissociation_time` design. No target method array was emitted; the first terminal inspect accounted for both preflight IDs once and recorded `FAILED`. This gate is historical evidence only and cannot authorize a retry. |
| `stage3_breast_corrected_assay_dissociation_ygg_20260915T154236Z` | The new global-column Breast Stage 3 one-row attempt failed before scheduler submission because the submitter accepted only the exact eight-row corrected recovery or the one-row Alzheimer follow-up. No scheduler IDs, output, or reviewer evidence. Preserve the manifest as failed evidence; a targeted selector and fresh run are required. |
### Validator-only reports and accepted historical artifacts

- The source-bound H5AD audit from snapshot
  `b8f7aec15dfc91bc493caeb04f3d2cc066db8c5d` passed
  `Covid19_PBMC`, `Kidney_KPMP_full`, `Diabetes`, and `Lung` with zero split
  samples. It failed old corrected-source assumptions for Alzheimer (`assay`
  disagreement within 21 donor samples), Breast_cancer (65,359 literal
  `unknown` dissociation-time cells), and Lupus_PBMC (`batch_cov` disagreement
  within `sampleID`). These are evidence for targeted contracts, not a broad
  rerun authorization.
- The reviewed Stage 2 predecessor `stage2_joanito_final_20260912b` completed
  with scheduler IDs `4403663/4403664`; its watchdog recorded the explicit
  Joanito hook, 373,058 cells, 189 samples, current `seqtec`/`cell.type_new`,
  and the five-sample debug artifact.
- Validator-only corrected Stage 3 sync repair
  `stage3_corrected_sync_repair_20260913a` produced `STATE=OK`, verified all
  nine old corrected H5AD destinations, and retained the original failed
  terminal evidence. It is not a current Alzheimer or eight-dataset Stage 5
  selector.
- The old corrected-final consumer reports
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_runs/stage5_corrected_consumer_preflight_20260914a/manifests/consumer_contract.json`
  (MD5 `08fc88e0c99cb2bc0e2b2a9040f09996`) and `...20260914b/...` (MD5
  `219a9bf92bd6af8d8cf3a544388a774f`) validated eight historical rows under
  then-current assumptions and found the old Alzheimer sample universe had
  one assay level after donor grouping. They are evidence only for artifacts
  that pass the current contract; they do not override the mandatory 32-row
  recovery or authorize reuse of combined-key prepare caches. The old
  Alzheimer row is excluded.
- The old uncorrected selected-sync post-audit
  `stage5_uncorrected_final_20260913b/manifests/selected_sync_16_post_audit.tsv`
  recorded 16 validated method rows and synchronized payloads without copying
  H5ADs locally. This does not authorize a new broad uncorrected wave.
- The local lme4 experiment against the old 83-sample Alzheimer CLR bundle
  found assay levels `76/7`, sex levels `48/35`, `111/131` singular additive
  fits, `99/131` zero random-effect variance, and `112/131` convergence
  messages. This is historical evidence for the final global limma policy
  only; lme4 is removed from all new corrected work and its payloads are not
  reuse candidates.

## Historical execution checkpoint — 2026-09-15 (non-authoritative)
This checkpoint preserves earlier implementation and execution evidence only.
Its former eight-dataset/32-row wording is superseded by the confirmed
35-row matrix and the execution order near the top of this plan. It must not
derive a selector, release a gate, or authorize a launch.

### Completed

- Consolidated this file into the single authoritative plan. The separate
  Alzheimer/backup draft was archived at
  `.agents/plans/archive/1789391587-alzheimer-donor-assay-backup-plan.md`.
- Simplified `AGENTS.md` to durable scientific, artifact, snapshot, ownership,
  and gate rules while preserving the required baseline anchor
  `5302671ad94556edcf9acccf372d2dc34121d714`, the full HiTME/scATOMIC
  annotation contract, the user-authored plan-reference text, Yggdrasil as the
  default host, and the explicit Bamboo fallback policy.
- Updated `docs/ARCHITECTURE.md` to remain a general config-driven overview;
  exact dataset/method/row selections now belong here in the active plan.
  Updated `NOTES.md` with current Alzheimer, global-limma, 32-row, and backup
  evidence while retaining detailed historical notes. Updated the onboarding
  README with source/paper-reported major cell-type counts and the
  CombinedPBMC legacy/confounding rationale.
- Implemented and focused-tested the global separate-covariate limma boundary:
  corrected composition and pseudobulk use categorical factor designs and
  `limma::removeBatchEffect` with an intercept-preservation design; DESeq2
  corrected fitting is `design=~1`, `batch_col=NULL`, followed by limma.
  New identities are `limma_fixed_effects_v1` and
  `pseudobulk_limma_fixed_effects_v1`. lme4 and artificial combined keys are
  prohibited for new corrected artifacts.
- Implemented strict Alzheimer donor-by-assay Stage 2 source/worker/hook
  contracts, including exact assay tokens, 104 samples, 83/21 assay sample
  counts, 59/45 sex sample counts, collision/mixed-metadata rejection, raw
  immutability, and derivative-bound Stage 3 input validation.
- Implemented explicit corrected-final eight-row and one-row Alzheimer
  selector contracts, Stage 3 derivative binding, Stage 2 watchdog/common
  validation, keyed pseudobulk RDS validation, active identity whitelist
  migration, and focused test fixtures.
- Implemented the global dataset-level column contract in commit
  `bfa3ac63d3e9dd27d64093fe70be2b9857f8d431`: Breast uses the exact
  `assay`/`suspension_dissociation_time` pair, the corrected output is
  versioned, Parkinson view overrides are removed, and all production
  callers reject view-level column declarations. The commit was pushed and
  pulled exactly on Yggdrasil; focused configuration, source, and regression
  checks passed.
- Parent verification is green for shell syntax, Python compilation,
  H5AD/matrix/multibatch contracts, Stage 2 submitter/watchdog, Stage 3
  submitter, Stage 5 selection, benchmark matrix submitter/synchronization,
  H5AD preflight, batch registry, corrected limma, corrected consumer, and RDS
  contracts. Expected negative diagnostics and DESeq2/`cmdscale` warnings are
  non-fatal.
- Real-data smoke passed without writing artifacts using existing uncorrected
  Alzheimer Stage 5 files:
  `data/batch_effect/uncorrected/results/Alzheimer_batch_effect_uncorrected_composition.rds`
  (83 samples × 131 composition features) and
  `data/batch_effect/uncorrected/results/Alzheimer_batch_effect_uncorrected_Pseudobulk_hvg2000.rds`
  (83 samples × 2,000 pseudobulk features), with one-key and two-key designs,
  finite outputs, preserved identifiers, exact CLR row sums, and pseudobulk
  design rank 3.
- The remote-only transfer proof and repository backup succeeded without local
  staging. `ssh -A bamboo` reached
  `login1.yggdrasil.hpc.unige.ch`; source/destination POC SHA-256 was
  `ac5944fad030a07ad4257a3d7b7b44a83925c3e0fcba83196f9cbec6d670dcb2`, and the
  second checksum-aware dry run was empty. The repository clone is at
  `yggdrasil:~/scratch/_ecoda_backups/ECODA_paper_repo_20260914`; both clones
  report `751c3f7fd8d9a863d6940bc37b1269fa785c06d4`, and the full checksum
  dry run was empty. No Mac staging was used.

### Open items and exact next steps

1. **Alzheimer Stage 2 is provisionally complete.** The exact one-step gate
   ran worker `4407671` and watchdog `4407672`; both completed with exit
   `0:0`. The derivative validator passed for 1,395,601 cells and 104
   `donor_id_assay` samples. The local gate retains a completion-transport
   `PRELAUNCH_STOP`; its accounting/artifact audit evidence is preserved and
   requires explicit reviewer disposition before formal release.
2. **Yggdrasil scratch working tree is active.** The transfer-sanity-passed
   scratch mirror was moved from its backup path to the active canonical
   `~/scratch/ECODA_paper` path on Yggdrasil. It is now a writable working
   tree, not an independent immutable backup. The transfer used no
   `--delete` or `--inplace`; exact equality and content hashes remain
   intentionally deferred under `TRANSFER_SANITY_ONLY=PASSED`.
3. **Repository backup transfer completed with transfer sanity verification.**
   The timestamped Yggdrasil copy
   `ECODA_paper_repo_20260914T212143Z_69a7443` is approximately `23G` and
   includes source commit
   `69a744344c6ce0cb1a91a3904a07daabc6bb8070`. The subsequent plan-status
   commit(s), including this checkpoint, are not in that mirror; treat it as
   a historical repository backup rather than current compute provenance.
   The full content checksum dry run was stopped after exceeding the
   practical time budget; its repository status scan was not used as a gate.
   Record `CONTENT_CHECKSUM=DEFERRED` and rely on the successful rsync marker,
   destination presence, commit identity, and optional rough size sanity.
4. **Yggdrasil portability checks pass for the scoped CPU lanes.** Canonical
   `~/scratch/ECODA_paper` and `~/ECODA_paper` are present, the mirrored
   `.pixi` environment provides Python `3.13.14` and R `4.5.2`, the FORMAT 2
   runtime smoke passes, and CPU Slurm execution passes on `cpu001`
   (`45683950`). Host-aware defaults select `shared-cpu,shared-gpu` and
   `public-gpu`; local NAS routing and the Ygg profile tests pass. NASAC is
   unmounted, so the explicit local scratch result mirror is required.
   All approved batch-effect lanes run on CPU; GPU compatibility is not a
   blocker because GPU resources were only required for benchmark views, which
   are complete. `PORTABILITY_AUDIT=READY_FOR_SCOPED_LANES`; no broad
   selection or partial method subset is allowed.
5. **Wave 1 execution status.** The first combined Stage 3 and four-method
   Stage 5 attempts are preserved as failed pre-scheduler evidence. The
   producer-bound Alzheimer Stage 3 retries used separate one-row manifests.
   The exact uncorrected inspection for
   `stage3_alzheimer_uncorrected_ygg_retry_20260915T115439Z` completed at
   `2026-09-15T14:37:25Z` with `state=FAILED`, `audit.passed=false`, and
   `release_eligible=false`; its one accounting query returned
   `45686974|OUT_OF_MEMORY|0:125`, `45686976|COMPLETED|0:0` (watchdog),
   `45686982|OUT_OF_MEMORY|0:125`, and
   `45687166|COMPLETED|0:0` (retry2 array task
   `45687166_1`). The exact corrected inspection for
   `stage3_alzheimer_corrected_ygg_retry_20260915T115439Z` completed at
   `2026-09-15T14:37:21Z` with `state=FAILED`, `audit.passed=false`, and
   `release_eligible=false`; its one accounting query returned
   `45686975|OUT_OF_MEMORY|0:125`, `45686977|COMPLETED|0:0` (watchdog),
   `45686981|OUT_OF_MEMORY|0:125`, and
   `45687165|COMPLETED|0:0` (retry2 array task
   `45687165_1`). The corresponding watchdog ownership is
   uncorrected `45686976` → `45687166_1` and corrected `45686977` →
   `45687165_1`.
   The corresponding run-scoped `ecoda_run_audit.sh` attempts were
   validator-only: each used its immutable snapshot `source.manifest`, the
   run-root selection, and the run-root `runtime.identity`; each timed out
   after 300 seconds while processing the large H5AD and produced no passing
   report. Neither audit attempt created scheduler work. The retry2 H5AD
   outputs remain preserved; do not rerun preprocessing merely to erase OOM
   history. A fresh recovery gate may accept only validated retry2
   array/watchdog outcomes while preserving this complete OOM attempt chain.
   The selected-eight corrected-final Stage 5 retry
   `stage5_eight_corrected_ygg_parallel_20260915T131349Z` also remains failed:
   all source checks passed, seven dataset consumer rows passed, and the
   Breast row failed on the rank-deficient three-key design. Preflight IDs
   `45687614` and `45687624` were accounted for exactly once; no target method
   array ran. Its durable sidecar report remains historical evidence at
   `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_logs/20260915T123309Z_stage5_inputs/SIDECAR_VALIDATION.tsv`.
   The current global policy is
   `columns.batch=["assay","suspension_dissociation_time"]`. The first
   replacement Breast Stage 3 gate
   `stage3_breast_corrected_assay_dissociation_ygg_20260915T154236Z` failed
   before scheduler submission on its unsupported one-row classification. It
   emitted no scheduler IDs or H5AD. The targeted `--corrected-recovery`
   Breast-one-row selector is implemented in source commit
   `63c124568aa3f2ccb30e27553a799707b43133da`; its focused submitter test
   passes, and that exact revision was pulled on Yggdrasil. Fresh snapshot
   `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/stage3_breast_corrected_assay_dissociation_retry2/63c124568aa3f2ccb30e27553a799707b43133da`
   and one-row selection are sealed. Gate
   `stage3_breast_corrected_assay_dissociation_retry2_ygg_20260915T161013Z`
   launched once at `2026-09-15T16:11:32Z`; its single durable wait is armed.
   Validate its terminal H5AD before preparing the selected-eight Stage 5 gate
   with the exact 32-row scope. Never reuse the failed manifests.
6. **Alzheimer Stage 5 follow-up.** After the two Alzheimer Stage 3 lanes
   reach terminal accounting, run-scoped audit, synchronization, and required
   review, run the two explicit one-row Alzheimer Stage 5 lanes under the
   ordering and shared-root rules above. Keep all existing validated artifacts
   outside every recomputation selection.
7. **Finalize only after reviewed artifacts.** Synchronize manifest-listed
   outputs and checksums, then execute the final analysis lane.
