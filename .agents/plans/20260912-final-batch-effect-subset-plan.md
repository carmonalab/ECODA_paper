# ECODA final batch-effect subset processing

## Current scope — authoritative (2026-09-14)

This section is the single active contract for this plan. It supersedes every
older active Context, Approach, Verification, Assumptions, and status section
that conflicts with it. The incident appendix is evidence only; it does not
create selectors, authorize a launch, or preserve an obsolete nine-row scope.
This maintenance change is plan-only: it makes no source/configuration edits,
runs no validation command, and launches no HPC, scheduler, backup, or SIF
operation.

### Objective and hard boundaries

- First complete the local selector/source-contract implementation and the
  backup-feasibility work. No compute may start before both are complete and
  reviewed.
- The production source is the current `datasets.json` plus authoritative HPC
  data. Local mirrors and historical JSON reports are diagnostic evidence only.
- Existing H5ADs, RDS bundles, pseudobulk caches, Feather files, manifests,
  checksums, logs, and plots are immutable. Reuse is validator-only and
  artifact-by-artifact; no broad `--force`, historical matrix, inferred scope,
  or overwrite is allowed.
- Batch-effect views do not invoke Pipeline 4 annotation. Preserve configured
  source/author cell-type columns and keep biological labels evaluation-only.
- Bamboo remains the default compute host for all production work. Yggdrasil
  is never an implicit compute host; it is a temporary backup destination
  only during the 2026-09-15–18 Bamboo maintenance window and only when
  explicitly named by the approved backup plan or by the user.
- Every future full-cohort operation uses the checked-in
  `durable-hpc-gate-ecoda` workflow, an exact run-owned selection, immutable
  source/runtime/auxiliary identities, atomic outputs, checksums, one
  unbounded durable wait, one terminal inspection over every emitted ID, and
  the required Luna Max review.

### Phase order

1. **Local selector and source-contract implementation.** Implement and
  focused-test the exact eight-row corrected-final selector and mandatory
  32-row corrected Stage 5 recovery (8 prepare, 8 pseudobulk, 8 GloScope,
  8 composition), the strict Alzheimer donor-by-assay selector/derivative
  contract, corrected method path and sync validation, fixed-effect limma
  consumers, and all no-op/reuse guards. Do not launch compute while these
  contracts are incomplete.
2. **Backup feasibility before compute.** Resolve canonical paths, access,
  quotas, capacities, symlink targets, and active writers; prove a tiny
  non-production transfer; then establish the alternate-cluster clone plan
  below. NAS is not a repository or scratch backup. It receives only the
  explicitly expected processed/results data and its checksums/sidecars.
3. **Eight-dataset corrected-final recovery.** Seal the full-hash
  source/runtime snapshot and exact eight-row/32-method-row manifest, then
  run this gate if and only if the maintenance feasibility envelope below is
  provable.
4. **Alzheimer follow-up.** Only after the eight-row snapshot is sealed may
  the Alzheimer derivative/configuration be changed and snapshotted. Run the
  strict donor-by-assay Stage 2, one-row Stage 3 views, and one-row Stage 5
  lanes under the contracts below.

### Eight-dataset corrected-final Stage 5 recovery

This is the first compute target and is exactly one corrected-final selection,
not the historical matrix and not a dynamically expanded config selection. The
order is contractual:

```text
Joanito
Stephenson
Breast_cancer
Covid19_PBMC
Kidney_KPMP_full
Diabetes
Lupus_PBMC
Lung
```

The headerless run-owned rows are:

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

The only target methods are:

```text
prepare_pseudobulk,pseudobulk,gloscope,composition
```

`prepare_pseudobulk` remains a declared target and dependency, and all eight
old corrected prepare caches use the prohibited
`__ecoda_batch_combined_v1` identity. They are immutable-but-stale historical
artifacts, must be recomputed for this recovery, and are ineligible for reuse
under the separate-covariate limma contract. Therefore the current pending
recovery manifest contains exactly 32 method rows: eight
`prepare_pseudobulk`, eight `pseudobulk`, eight GloScope, and eight
composition rows. All eight prepare rows are mandatory in this recovery; no
validator may remove or reuse one of the stale caches.
Valid existing `MRVI`/`PILOT`/`QOT` rows are explicitly skipped and remain
outside this selection; they are not recomputed to satisfy a historical
seven-method matrix. Any other artifact can be reused individually only after
it passes the current source, configuration, method, root, checksum,
ownership, and model-identity contracts, and never if its correction identity
uses a combined key or lme4. Prior lme4-produced corrected composition or
pseudobulk results do not satisfy the new corrected contract and remain
untouched historical evidence.

The gate uses `--pass corrected --analysis-variant corrected_final` and the
variant-qualified roots:

```text
${HPC_SCRATCH_DIR}/batch_effect/corrected_final
${NAS_TARGET_DIR}/batch_effect/corrected_final
```

All output paths, run metadata, ownership records, watchdogs, validators,
sync lists, and execution logs must use this root and the
`_batch_effect_corrected_final_` stem. No corrected result may be reconstructed
under the legacy pass root. The single 32-row gate owns the shared
corrected-final root and its synchronization boundary; independent dataset
and method rows may dispatch concurrently within this gate subject to
declared dependencies, but no other corrected-final gate may overlap or share
this root.

### Final corrected-method policy

This is a final policy, not an experiment or an alternative:

- **Every corrected composition and corrected pseudobulk mode uses limma
  fixed effects with the original technical covariates as separate design
  columns.** For effective keys, use a fixed internal-alias design such as
  `model.matrix(~ 1 + technical_key_1 + technical_key_2 + ...)`.
- Never construct a combined/artificial batch key such as an interaction or
  concatenated `batch_key` in place of the separate technical columns. Keep
  `assay`, `sex`, and other configured technical fields separately visible in
  metadata, source identity, manifests, and the exact recorded design string.
- Biological labels (including Alzheimer `Cognitive status`) and sample IDs
  are never model covariates. A one-level technical field is retained as
  metadata but omitted from the effective design; if no technical field varies,
  record `NO_CORRECTION` and return the uncorrected object for that method.
- Fail closed on missing metadata, non-finite output, rank deficiency, or
  non-positive residual degrees of freedom. Do not silently drop an aliased
  technical column. Remove only the non-intercept technical contribution,
  preserve the intercept, and restore the CLR row-sum invariant for
  composition.
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
view-specific selections:

```text
Alzheimer<TAB>batch_effect_uncorrected
Alzheimer<TAB>batch_effect_corrected
```

Each has its own run root, output path, source/runtime identity, ownership and
checksums. The two gates may run concurrently only when exact manifests prove
disjoint mutable output paths, metadata/export/checksum owners, serialization
locks, and source snapshots. Otherwise they are serialized. Do not use a
full configured corrected selector that could include the eight other
cohorts. Both outputs must contain the ordered 104 `donor_id_assay` samples,
preserve original technical metadata and raw counts, exclude biological labels
from processing covariates, and use semantic uncorrected/corrected
representations respectively.

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
This is a separate later Alzheimer follow-up, not part of the current
eight-dataset corrected-final recovery or its mandatory 32-row manifest.

Valid rows may be skipped only after validator-only source/metadata,
checksum, ownership, and model-contract checks. Pre-derivative Alzheimer
rows are not valid reuse candidates. The uncorrected lane uses
`batch_effect/uncorrected_final`; the corrected lane uses
`batch_effect/corrected_final` and the final limma policy above. The corrected
Alzheimer lane is mandatory **after** the eight-dataset corrected-final gate
has terminal accounting, run-scoped audit, synchronization, and Luna Max
review. A second durable serialization group is not a safe bypass. The
uncorrected Alzheimer lane may overlap another lane only after a focused
review proves its roots, NAS destinations, checksum/execution-time files,
artifact owners, and locks are disjoint; the default is serialization.

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
- Corrected-final Stage 5 gates that target
  `batch_effect/corrected_final` share the Stage 5 synchronization owner,
  checksum merge, execution-time merge, and NAS destination. They **must** use
  the `ecoda-benchmark` policy group and serialize: the eight-dataset
  32-row gate completes terminal audit/review before the one-row Alzheimer
  corrected gate begins. A different durable group cannot make same-root
  synchronization safe.
- Safe parallelism is limited to genuinely disjoint contracts: local
  selector/limma implementation and backup feasibility; the sealed eight-row,
  32-method-row corrected-final compute and Alzheimer Stage 2 if their
  snapshots, roots, and owners are disjoint; and the two Alzheimer Stage 3
  rows only when the disjointness test passes. The eight-row corrected-final
  gate and Alzheimer corrected-final Stage 5 gate never overlap. Alzheimer
  uncorrected Stage 5 may overlap only with an explicit disjoint-root proof.
- A complete repository or scratch clone must not run against active writers.
  Take a quiescence checkpoint and hold stable path/size/hash manifests before
  declaring that backup complete. A failed row releases only its own target;
  successful rows and validated sync fragments remain immutable.

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

Home is backed up; scratch is operational and not backed up. `$HOME/scratch`
is a symlink, so backup, snapshot, and gate commands must resolve the
canonical `/srv/beegfs/scratch/users/h/halterc` path. Bamboo remains the
default compute host. The user's local setup reaches Yggdrasil with
`ssh yggdrasil`; no local staging is needed.
Yggdrasil is a temporary backup destination only for the 2026-09-15–18 Bamboo
maintenance window when explicitly named by the approved backup plan or user,
never an implicit compute host. The direct Bamboo→Yggdrasil route is proven
conditionally via `ssh -A bamboo` and the documented FQDN
`login1.yggdrasil.hpc.unige.ch` while the originating Mac agent session
remains alive.
The remote-only POC root
`/srv/beegfs/scratch/users/h/halterc/_ecoda_backup_poc_20260914_direct`
transferred a non-production file without local staging: source and destination
SHA-256 both matched
`ac5944fad030a07ad4257a3d7b7b44a83925c3e0fcba83196f9cbec6d670dcb2`, and the
second checksum-aware dry run was empty.
The Bamboo `yggdrasil` alias remains optional/unconfigured and currently fails
DNS; do not treat its local success as evidence for a Bamboo transfer.
Agent-forwarded SSH (for example, `ssh -A`) depends on the originating Mac
session and agent remaining alive: a sleeping or disconnected Mac can break a
multi-hour Bamboo-launched `rsync` even inside `tmux`. Require either an alive
agent session or a separately approved Bamboo key registered centrally; never
copy or store a private key or password.

Execute the backup priority in this order, without copying an active tree:

1. Inventory canonical source/destination paths, cluster access, quotas,
   capacities, file counts, symlink targets, and active writers.
2. After explicit backup authorization, retain the direct POC above as the
   tiny-transfer proof and require an alive forwarded agent session or a
   separately approved Bamboo key registered centrally for any multi-hour
   transfer. Compare SHA-256/MD5, size, permissions, and a second no-change
   dry run. This proof and authorization are required before cloning; no
   repository or scratch clone is feasible or authorized before they succeed.
3. Clone the complete repository (including `.git`, hidden files, source,
   plans, and runtime references) to an explicitly authorized, timestamped
   Yggdrasil backup path by resumable, non-destructive `rsync`; keep a sorted
   path/size/hash manifest and do not use `--delete` or `--inplace` for the
   first mirror.
4. At a quiescence checkpoint, clone the complete scratch `ECODA_paper` tree
   to that explicitly authorized Yggdrasil backup path in top-level batches.
   Include raw, processed, results, logs, gates, manifests, and checksums;
   retain per-batch logs, source/destination hash manifests, and a second dry
   run. Do not rely on scratch retention.
5. Send NAS only the expected validated processed/results artifacts, metadata,
   checksums, sidecars, and explicit manifests required by the active lanes.
   Do not copy the complete repository, control plane, logs/gates, full
   scratch tree, raw H5ADs, or unrelated legacy artifacts to NAS.
6. Restore small and large samples from the Yggdrasil backup clone, compare
   complete manifests and source/runtime identities, and obtain human review
   before declaring the clone complete.

### Maintenance feasibility gate

The hard maintenance boundary is:

```text
2026-09-15 08:00 +0100
= 2026-09-15 07:00 UTC
= 2026-09-15 09:00 Bamboo local time (UTC+0200)
```

The last recorded Bamboo time check was `2026-09-14 12:22 UTC`. A launch is
permitted only when the complete chain—wrapper and source preflight, queue
start, dependency release, worker/watchdog compute, terminal accounting,
artifact validation, review, backup/synchronization, and recovery margin—fits
before that boundary. The planning envelope is **at most 8 hours of submitted
worker/watchdog compute plus 2 hours of operational preparation, queue/review,
and synchronization**. Treat this as a feasibility budget, not merely a
worker walltime request; reject a launch when the bound and a safety margin
cannot be demonstrated.

## Release checklist

Before the first eight-dataset compute gate, the run owner must have:
- completed selector/source-contract implementation and focused contract
  checks, with the exact eight dataset rows, four corrected target methods,
  and mandatory 32-row recovery scope recorded;
- completed backup inventory and tiny transfer proof, with alternate-cluster
  repository/scratch clone capacity and quiescence plan recorded;
- marked all eight old `__ecoda_batch_combined_v1` prepare caches
  immutable-but-stale historical artifacts; each must be recomputed for this
  recovery and none is eligible for reuse. All eight prepare rows remain
  mandatory in the 32-row selection;
- Valid `MRVI`/`PILOT`/`QOT` rows remain outside the eight-row pending selection;
- recorded the fixed-effect limma identity/design for corrected composition and
  pseudobulk, separate original technical covariates, no composite key, and
  no lme4 path;
- sealed a full-hash source/runtime snapshot under a verified parent lock
  and proved the 8h + 2h cutoff envelope; and
- recorded one durable gate command, selector checksum, expected rows, roots,
  and dependency/review boundary.

Before Alzheimer work, additionally require the reviewed eight-dataset gate,
a new full-hash source/runtime snapshot, the strict derivative schema and
104-sample contract, and explicit one-row Stage 3/Stage 5 manifests. The
Alzheimer corrected-final Stage 5 lane is serialized behind the completed
eight-dataset corrected-final root; no stale nine-row artifact or donor-only
metadata may authorize it.

After all approved lanes, synchronize only manifest-listed processed/results
artifacts and metadata/checksums to the workstation/NAS. Final analysis must
load paths from explicit manifests, may read approved legacy result artifacts
read-only where declared, never read legacy H5ADs for convenience, and must
write only its lane-specific output root. Existing legacy artifacts and their
modification times remain untouched.

## Historical evidence appendix — non-authoritative

This compact record preserves only major implementation evidence, terminal
gate outcomes, and validator reports. It is not an active scope and does not
authorize reuse or relaunch. The old nine-row corrected launch is historical;
its Alzheimer row used donor-only samples and is superseded by the
donor-by-assay contract. The current eight-dataset corrected-final gate still
requires all 32 method rows; the eight old combined-key prepare caches and
all lme4 payloads remain immutable historical artifacts, not reuse candidates.

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

## Execution checkpoint — 2026-09-14

### Completed

- Consolidated this file into the single authoritative plan. The separate
  Alzheimer/backup draft was archived at
  `.agents/plans/archive/1789391587-alzheimer-donor-assay-backup-plan.md`.
- Simplified `AGENTS.md` to durable scientific, artifact, snapshot, ownership,
  and gate rules while preserving the required baseline anchor
  `5302671ad94556edcf9acccf372d2dc34121d714`, the full HiTME/scATOMIC
  annotation contract, the user-authored plan-reference text, Bamboo as the
  default host, and the explicitly named Yggdrasil backup exception.
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

1. **Commit/push this implementation checkpoint.** Review the explicit
   allowlist, preserve unrelated user worktree changes, and commit/push the
   source, tests, documentation, and plan changes. `datasets.json`,
   `pixi.toml`, and `pixi.lock` have not been changed by this checkpoint.
2. **Finish backup before production compute.** The repository clone is
   verified; the complete 2.3 TB scratch clone has not started. Use only
   direct Bamboo→Yggdrasil rsync, never Mac staging. Agent forwarding is
   conditional on the originating Mac session remaining alive; a persistent
   centrally registered Bamboo key is preferable for a multi-hour transfer.
   Quiesce writers, split scratch into explicit top-level batches, preserve
   logs/gates/manifests/checksums, and verify source/destination manifests.
3. **Run the eight-dataset corrected-final gate.** After the source checkpoint
   is snapshotted, use the exact eight-row selection and mandatory 32 method
   rows (8 prepare, 8 pseudobulk, 8 GloScope, 8 composition). The old eight
   combined-key prepare caches remain immutable and stale; none may be reused.
   Launch only if the `2026-09-15 07:00 UTC` cutoff envelope is provable.
4. **Run the Alzheimer follow-up from a new configuration snapshot.** Preserve
   raw H5AD, execute the guarded Stage 2 derivative, then the explicit
   derivative-bound Stage 3 views and one-row Stage 5 lanes. Corrected Stage 5
   serializes behind the eight-dataset corrected-final gate because the root
   and synchronization owner are shared.
5. **Finalize analysis only after reviewed artifacts are available.** Sync
   explicit processed/result files and manifests to NAS; do not copy full
   H5ADs to the Mac or run final notebooks from stale local mirrors.

No new full-cohort Stage 2/3/5 job, durable gate, or full scratch transfer is
authorized merely by this checkpoint. The next session must re-read this
checkpoint, verify the working-tree allowlist, and continue with item 1 or the
explicitly authorized backup item.
