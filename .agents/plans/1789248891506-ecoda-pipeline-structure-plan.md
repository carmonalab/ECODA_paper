# ECODA Pipelines 2--5 structural review and simplification

## Status

Proposed follow-up plan. This plan is separate from
`.agents/plans/20260912-final-batch-effect-subset-plan.md` and is not authorized
for implementation in the final batch-effect execution session. The current
session records triage only; no code, artifact, gate, runtime, SIF, or HPC
changes belong here.

## Origin and disposition

A read-only reviewer assessed the current final batch-effect implementation and
Pipelines 2--5. The review correctly identifies architectural growth and
repeated orchestration/contract code, but it also contains stale findings from
a moving worktree and placeholder commit identifiers. Every claimed line
number, checksum, entry point, and dead-code status requires independent
verification before action.

The immediate decision is to separate **trust-boundary correctness** from
**structural simplification**:

- Trust-boundary blockers remain in the final batch-effect plan and must be
  resolved before any resumed runtime or gate work.
- Broad deletion, DRY rewrites, checksum-format changes, and cross-language
  contract retirement are deferred to this plan.
- No recommendation below authorizes deletion or refactoring by itself.

## Non-goals for the final batch-effect execution session

Do not, as part of the current final batch-effect run:

- delete annotation stubs, forwarding shims, diagnostic scripts, or Stage 2
  wrappers;
- rewrite the Stage 2--5 submitters around a new shared orchestration layer;
- replace or normalize all MD5 sidecar implementations;
- retire `src/utils/batch_contract.R` in favor of Python;
- change manifest, ownership, runtime, checksum, or artifact schemas;
- rerun valid cohorts or launch exploratory HPC work;
- rebuild the runtime solely to accommodate a structural refactor.

These changes require a call-graph, contract, compatibility, and migration
review first.

## Review findings to investigate

### 1. Orchestration size and duplicated boundaries

The submitters and shared shell helpers contain large amounts of manifest
parsing, atomic file handling, ownership transitions, scheduler state, retry
logic, and terminal auditing. Similar stage-prefixed functions appear in
Pipelines 2--5.

Investigate whether each repetition is:

- genuinely stage-specific and safer at the boundary;
- a pure helper that can be centralized without changing semantics; or
- obsolete compatibility code with no live callers.

Do not infer that similar names imply interchangeable behavior. Record inputs,
side effects, failure states, and ownership assumptions for every candidate.

### 2. Checksum and sidecar duplication

Bash, Python, R worker code, and the analysis loader each read or write
`MD5`, `SIZE`, and `PATH` sidecars. The review’s count of validation passes is a
risk signal, not proof that all checks should be removed.

First define a compatibility matrix covering:

- accepted sidecar fields and ordering;
- path canonicalization and symlink policy;
- MD5/SHA-256 roles;
- immutable artifact ownership records;
- producer/run identity requirements;
- legacy artifacts with missing ownership records.

Only then consider centralizing pure parsing. Boundary-specific checks must
remain if they protect different trust transitions.

### 3. Dual Python/R batch contracts

`src/utils/py/batch_contract.py` and `src/utils/batch_contract.R` encode related
corrected-mode semantics. Retiring one implementation is high risk because
existing artifacts and workers depend on byte-level identity and cross-language
parity.

Before any change, establish one canonical serialized contract or prove that
both implementations are equivalent across:

- key normalization and ordering;
- UTF-8 sorting and hashing;
- scalar versus composite modes;
- missing/sentinel handling;
- within-sample constancy;
- rank and near-unique checks;
- persisted validation summaries.

Migration must preserve validation of existing artifacts and must not alter the
current run’s source/runtime identity.

### 4. Hardcoded paths and deployment assumptions

Audit institutional defaults such as NAS, scratch, user email, and container
prefixes separately from true correctness defects. The Bamboo `$HOME/scratch`
symlink is an operationally proven issue; other defaults may be intentional
deployment configuration.

Classify each path as:

- required site configuration;
- safe environment override;
- source/runtime identity field; or
- accidental hardcode requiring removal.

Do not replace a trusted path with a weaker fallback merely to shorten code.

### 5. Candidate orphaned scripts

The review names short stubs, forwarding shims, diagnostic scripts, and
single-dataset wrappers as possible dead code. Each candidate requires a
repository-wide reference and operator-entrypoint audit before deletion.

Required evidence for deletion:

- no shell, documentation, scheduler, or user-facing caller;
- no compatibility purpose recorded in repository rules;
- replacement entry point is explicit and documented;
- focused regression coverage proves the supported path;
- no existing gate or artifact manifest refers to it.

## Proposed work sequence

1. **Read-only inventory:** map Pipeline 2--5 entry points, callers, wrappers,
   workers, watchdogs, validators, manifests, ownership records, checksums,
   and runtime boundaries.
2. **Contract map:** document the canonical schemas and state transitions for
   source snapshots, run roots, selection manifests, output ownership,
   artifact records, sidecars, and terminal reviews.
3. **Duplication classification:** mark each repeated helper as intentional,
   safely centralizable, or high-risk; quantify actual semantic differences.
4. **Compatibility design:** choose one small refactor candidate and define
   before/after behavior, migration handling, and failure invariants.
5. **Incremental implementation:** refactor one boundary at a time with
   focused contract tests; no broad rewrite or simultaneous schema change.
6. **Orphan cleanup:** delete only candidates that pass the reference and
   compatibility evidence above.
7. **Post-change audit:** verify all entry points, artifacts, checksums,
   ownership states, legacy paths, and runtime identities remain compatible.

## Acceptance criteria

This plan is ready for implementation only when a future review provides:

- a complete caller/entry-point map for Pipelines 2--5;
- a contract compatibility matrix for manifests, sidecars, owners, and
  runtime/source identities;
- an evidence-backed orphan list with explicit replacements;
- a bounded first refactor with no artifact-schema migration hidden inside it;
- focused tests defending each changed observable boundary;
- an explicit decision that no current validated artifact or scientific
  invariant is invalidated.

Until then, preserve the current pipeline implementation and keep all final
batch-effect execution blockers in the approved execution plan.
