---
name: durable-hpc-gate-ecoda
description: >-
  Applies durable-hpc-gate to ECODA on Bamboo, adding repository invariants,
  authoritative-plan checkpoints, canonical scratch/NAS contracts, clean
  benchmark serialization, and Luna Max task/reviewer requirements.
---

# Durable HPC Gate — ECODA/Bamboo Profile

## Overview

This is the ECODA policy layer for the global `durable-hpc-gate` skill. It
keeps long-running Slurm wrappers on Bamboo durable across SSH/session loss
without copying the global orchestration implementation. The profile is
repository-aware: it checks the ECODA checkout, canonical configuration and
artifact roots, scientific invariants, benchmark wave serialization, and the
terminal accounting/audit contract.

The profile contains no credentials, email addresses, scheduler secrets, or
alternate runtime. Use the global CLI and the checked-in profile at
`references/profile.json`; do not add a second runner or wrapper here.

## Dependencies

- **`durable-hpc-gate` (required runtime dependency):** supplies the
  stdlib-only `uv run` CLI, durable manifest, remote runner, event-driven
  waiter, monotonic lifecycle, and terminal inspection.
- **`workflow-skill-creator` (authoring only):** used to structure this skill;
  it is not needed to run an ECODA gate.

## Quick Start

Run these commands from the ECODA repository. Substitute a new gate ID,
wrapper, and artifact paths for the gate being run. Resolve Bamboo's home
directory remotely; local `$HOME` is not a valid prefix for remote artifact
paths. The remote checkout is `${BAMBOO_HOME}/ECODA_paper`; its scratch data
root is `${BAMBOO_HOME}/scratch/ECODA_paper` (not a git clone).

```bash
GLOBAL="$HOME/.agents/skills/durable-hpc-gate/scripts/durable_hpc_gate.py"
PROFILE="$PWD/.agents/skills/durable-hpc-gate-ecoda/references/profile.json"
GATE_ID="ecoda_<unique-id>"
MANIFEST="$PWD/.gate/${GATE_ID}.json"
BAMBOO_HOME="$(ssh bamboo 'printf %s "$HOME"')"
REMOTE_WORKDIR="${BAMBOO_HOME}/ECODA_paper"
WRAPPER='cd "$HOME/ECODA_paper" && source src/slurm_config.sh && <serialized-wrapper>'
REMOTE_ROOT="${BAMBOO_HOME}/scratch/ECODA_paper/gates/${GATE_ID}"
REMOTE_RUNNER="${REMOTE_ROOT}.runner.sh"

uv run "$GLOBAL" prepare \
  --profile "$PROFILE" --manifest "$MANIFEST" \
  --output "$PWD/.gate/${GATE_ID}.prepare.json" \
  --project ECODA --gate-id "$GATE_ID" --remote-host bamboo \
  --remote-workdir "$REMOTE_WORKDIR" --exact-command "$WRAPPER" \
  --serialization-group ecoda-benchmark \
  --tmux-session "$GATE_ID" --completion-channel "file:${REMOTE_ROOT}.done" \
  --remote-manifest "${REMOTE_ROOT}.manifest.json" \
  --remote-runner "$REMOTE_RUNNER" \
  --remote-log "${REMOTE_ROOT}.log" --remote-status "${REMOTE_ROOT}.status.json"


uv run "$GLOBAL" reconcile --manifest "$MANIFEST" --profile "$PROFILE" \
  --output "$PWD/.gate/${GATE_ID}.reconcile.json"
uv run "$GLOBAL" launch --manifest "$MANIFEST" --profile "$PROFILE" \
  --output "$PWD/.gate/${GATE_ID}.launch.json"
```

`serialization_group` is the required explicit cross-manifest mutex identity.
Launch takes a deterministic lock keyed by project, profile, profile digest,
and group; an active `RUNNING` sibling with that identity records a
`resource_lock_conflict` discrepancy, moves this gate to `PRELAUNCH_STOP`, and
is not launched. Distinct groups may be prepared and launched in parallel,
including gates that share a reviewed predecessor. Reviewed dependency lineage
is independent of `serialization_group`.

After `launch` returns, the harness owns exactly one durable waiter. Start it
as an asynchronous, unbounded call (`async=true`, `timeout=0` in the harness;
do not keep a task agent alive and do not send periodic polling requests):

```bash
uv run "$GLOBAL" wait --manifest "$MANIFEST" \
  --output "$PWD/.gate/${GATE_ID}.wait.json"
```

When the completion event wakes the workflow, a fresh bounded Luna Max task
runs one terminal inspection and records the evidence. Before that invocation,
the task reads the durable wrapper log and extracts every scheduler array or
watchdog ID the wrapper emitted after launch; it passes those values explicitly
on the first inspect and never edits the manifest by hand. A passing audit is
not release-eligible until the Luna Max reviewer explicitly approves it through
the second inspect invocation:

```bash
uv run "$GLOBAL" inspect --manifest "$MANIFEST" --profile "$PROFILE" \
  --scheduler-array-id "$ARRAY_ID" \
  --scheduler-watchdog-id "$WATCHDOG_ID" \
  --output "$PWD/.gate/${GATE_ID}.inspect.json"
uv run "$GLOBAL" inspect --manifest "$MANIFEST" --profile "$PROFILE" \
  --approve-reviewer "Luna Max" \
  --output "$PWD/.gate/${GATE_ID}.review.json"
```

Pass each emitted ID with a repeatable option; omit an option when that class
of ID was not emitted. The first inspect token-validates and atomically
deduplicates supplied IDs before setting `audit_state: IN_PROGRESS`. Supplying
IDs after audit begins or completes is rejected.

Use `status` for a non-mutating local/remote view. It is not a scheduler
poll and cannot replace `inspect`:

```bash
uv run "$GLOBAL" status --manifest "$MANIFEST" \
  --output "$PWD/.gate/${GATE_ID}.status.json"
```

Every subcommand requires `--output`; output files are the structured record,
while stdout stays concise. If a command cannot prove the launch identity or
terminal state, stop and reconcile rather than retrying the wrapper.

## Utility Scripts (if CLI-based)

The local profile has no utility scripts. All six approved subcommands belong
to the global CLI:

- `reconcile`: read-only adoption/recovery of the manifest against remote
  tmux/process/status and profile invariants; never launches or polls Slurm.
- `prepare`: validates the explicit gate specification and writes a
  `PREPARED` manifest atomically.
- `launch`: reconciles first, repairs only monitoring infrastructure when safe,
  and launches the exact command once in persistent tmux.
- `wait`: blocks on the completion event, then treats durable status as the
  authority if event transport reconnects.
- `inspect`: is legal only after terminal durable status; performs the single
  accounting query and profile artifact/checksum audits.
- `status`: reads durable state without mutation or scheduler polling.

Do not add aliases, implicit defaults, a second waiter, or a local copy of the
Python implementation.

## Workflow

### 1. Check the authoritative plan before preparing

- Use the existing authoritative `.kilo` execution-plan checkpoint. Record the
  gate ID, exact immutable wrapper, opaque serialization group, manifest, and
  profile.
- On successful launch set that checkpoint to `RUNNING`; after terminal
  inspection set it to `COMPLETED` or `FAILED` with links to status and audit
  evidence. Preserve the discrepancy history.
- Never create a competing plan, edit the current plan during skill
  installation, or mark a gate complete before `inspect` and review.

### 2. Enforce Bamboo and repository invariants

- Remote host is `bamboo`; the remote working directory is the repository
  clone `$HOME/ECODA_paper`. The clone is not `$HOME/scratch/ECODA_paper`.
- Source `src/slurm_config.sh` in every HPC wrapper. It is the authority for
  `PROJECT_ROOT`, `HPC_SCRATCH_DIR`, `NAS_TARGET_DIR`, `PYTHON_BIN`,
  `PIXI_RSCRIPT`, scheduler resources, and retry ceilings.
- Scratch data belongs under `$HOME/scratch/ECODA_paper`; NAS artifacts belong
  under the configured `NAS_TARGET_DIR`. Do not hard-code a second path or
  treat scratch as a git checkout.
- Workers use the immutable config-provided interpreter commands; never invoke
  bare `python`, `Rscript`, or an ordinary lock-mutating `pixi run` in a job.
  Login nodes are for editing, staging, synchronization, and submission;
  heavy preprocessing/benchmarking runs in Slurm.
- Preserve `pixi.toml` and the resolved `pixi.lock`. Do not mutate
  `datasets.json` without explicit user confirmation. Publication figures are
  fixed rather than removed.

### 3. Coordinate benchmark waves and reviewed lineage

- `serialization_group` is the required explicit cross-manifest mutex identity
  for a matching project, profile, and profile digest. Before any remote
  launch, `launch` takes its deterministic group lock and scans sibling
  authoritative manifests for an active `RUNNING` gate with that identity.
  Such a sibling records a `resource_lock_conflict` discrepancy, enters
  `PRELAUNCH_STOP`, and is not launched. Distinct groups may be prepared and
  launched concurrently, including gates that share a reviewed predecessor.
- A top-level benchmark wrapper declares reviewed predecessor lineage through
  repeatable absolute `--dependency-manifest` paths. Launch checks every
  predecessor for the same project, profile, and profile digest, state
  `COMPLETED`, a passing first audit, required reviewer approval, and a bound
  `local_manifest`, without requiring matching `serialization_group`.
  Dependency lineage is independent of the mutex identity.
- The exact command string in the manifest is immutable. Do not manually split
  a wrapper or rerun a wrapper after an ambiguous launch.
- Existing remote evidence may be adopted only after the complete remote
  manifest validates and binds every required identity/path field, including
  `serialization_group` and `dependency_manifests`. A tmux session must have
  exactly one pane across all windows whose `pane_start_command` is the
  declared runner; process fallback parses exact argv rather than accepting a
  runner-path substring. Any mismatch is a fail-closed stop.

### 4. Launch, wait, and inspect durably

- A bounded Luna Max **task** executes `reconcile → prepare → launch` and
  exits after updating the checkpoint. Main arms exactly one harness-owned
  asynchronous `wait` with no periodic update requests.
- The runner uses `/bin/bash` on Bamboo, starts in the declared
  `remote_workdir`, and changes to that directory before invoking the exact
  wrapper. It writes terminal status atomically and signals completion only
  after a durable primary or fallback failure status exists. Every terminal
  status must carry the current manifest `event_generation`; stale generations
  are rejected. Monitoring files/events may be repaired only when that cannot
  duplicate or disturb the wrapper; stale generation markers are ignored.
- A fresh bounded Luna Max **task** reads the durable wrapper log, extracts
  every scheduler ID emitted after launch, and runs the first `inspect`
  exactly once with repeatable `--scheduler-array-id`/
  `--scheduler-watchdog-id` options. The ECODA profile sets
  `require_scheduler_ids: true`; `inspect` fails closed before any accounting
  or audit command when neither an array nor watchdog ID was recorded or
  supplied. Otherwise it atomically records `audit_state: IN_PROGRESS` before
  issuing any query, then performs the only terminal scheduler accounting
  operation: one strict
  `sacct -n -P -X -j ... --format=JobIDRaw,State,ExitCode` query for all
  recorded scheduler IDs. Every ID needs an exact or array-task
  `JobIDRaw|State|ExitCode` row and every matched row must be
  `COMPLETED|0:0`; missing, partial, pending, failed, malformed, or nonzero
  rows fail the audit. The final evidence write atomically records
  `audit_state: COMPLETED`; an interrupted retry fails closed on
  `IN_PROGRESS` without repeating accounting or audits. Never poll
  `squeue`/`sacct` or use an event as evidence.
- A Luna Max **reviewer** invokes the second `inspect --approve-reviewer`
  phase. That phase records approval without rerunning audit/accounting.
  Only a reviewed `COMPLETED` predecessor releases a dependent gate;
  dependency lineage is independent of `serialization_group`. Before remote
  launch, the deterministic same-group lock rejects an active `RUNNING`
  sibling with the same project/profile/profile-digest/group, records
  `resource_lock_conflict`, and moves this gate to `PRELAUNCH_STOP`. Distinct
  groups may proceed in parallel, including gates sharing a reviewed
  predecessor.

### 5. Preserve scientific and artifact contracts

- Biological labels (`Status`, `sample.origin`, `cond`,
  `Disease_Identity`, and equivalent ground truth) are never preprocessing,
  HVG, normalization, batch-correction, embedding, or model covariates.
- Keep the established DESeq2 modes: default `blind=TRUE`, no batch column,
  no batch correction; batch-effect mode is explicitly batch-only with
  `blind=FALSE`, `batch_col=<batch>`, and `correct_batch=TRUE`.
- Preserve raw counts in `layers["counts"]`, sample identity/order in Feather,
  and features/labels/distances in atomic RDS result bundles. Reject missing,
  NA, or mismatched identifiers. Verify the canonical NAS
  `${NAS_TARGET_DIR}/benchmark/checksums.md5` is nonempty (`test -s`) before
  running `md5sum -c` or consuming RDS.
- Missing status, failed workers, exhausted OOM retry, incomplete accounting,
  checksum mismatch, or malformed wrapper-produced artifacts is a failure, not
  an optional warning.
- Gate-specific artifact shape, identifier, and provenance checks are defined
  by the authoritative execution plan and performed by the completion agent
  against artifacts the wrapper actually produced; do not invent generic
  provenance or source-manifest files.
- Use the `_debug` Joanito five-sample subset for routine verification. Do not
  run full cohorts for minor checks.

### 6. Apply fail-closed recovery

- An ambiguous launch becomes `PRELAUNCH_STOP`: preserve all evidence and
  reconcile before any retry. Never guess that the wrapper did not start.
- A nonzero wrapper, unexpected terminal accounting state, unknown process,
  path/config discrepancy, or audit failure records `FAILED`/`STOP`, starts no
  dependent gate, and leaves the checkpoint unresolved for review.
- Never manually synchronize partial outputs after wrapper failure. Do not
  delete logs, manifests, status, checksums, or discrepancy records. A later
  retry is a new explicitly reconciled gate, not an ad-hoc `rsync`.

## Profile severity and audit contracts

`references/profile.json` is the machine-readable policy. The following
severity rules govern any discrepancy it reports:

- **BLOCK:** unknown discrepancy; exact-command or gate-identity mismatch;
  unexpected tmux/process; missing or non-atomic status; wrapper nonzero;
  non-terminal/non-`COMPLETED` accounting; missing artifact; checksum, path,
  source-level label-leakage, or scientific-invariant failure; a failed
  gate-specific shape/identifier/provenance check defined by the authoritative
  plan and run by the completion agent; manual/partial synchronization. Do
  not launch a dependent gate.
- **APPROVAL_REQUIRED:** a discrepancy whose `kind` is explicitly listed in
  `policy.known_discrepancies`. Preserve the raw evidence and obtain explicit
  user approval before proceeding; approval does not erase or downgrade the
  recorded discrepancy.
- **PASS:** all exact identity, status, accounting, configured-root, checksum,
  immutable-fingerprint, and plan-defined scientific/artifact audits succeed.
  Only a reviewed `COMPLETED` result is releasable.

Audit contracts are evaluated only at the terminal `inspect` step:

1. gate/manifest identity and exact command;
2. atomic terminal status and completion signal;
3. one `sacct` query covering every recorded scheduler ID;
4. repository/configured scratch/NAS roots and wrapper-produced artifact
   presence;
5. nonempty canonical NAS checksums and immutable datasets/pixi/`HEAD`
   fingerprints; and
6. gate-specific artifact shape, identifier, and provenance checks defined by
   the authoritative execution plan and recorded by the completion agent,
   followed by the authoritative plan checkpoint and reviewer evidence.

Gate-specific artifact shape, identifier, and provenance checks are
intentionally not generic profile contracts. The authoritative execution plan
must name the existing wrapper-produced artifacts and expected checks, and the
completion agent must run and record those checks against those artifacts.
Do not require generic `provenance.json` or `source_manifest.json` files; their
absence is not a failure.

The source-level no-label-leakage rule remains mandatory (see Section 5).
This profile does not scan benchmark output strings for biological labels:
ground-truth labels are legitimate result content, and runtime output strings
cannot establish which covariates a model received.

The event is a wake-up mechanism, not proof. If transport is lost, reconnect
and read durable status; do not relaunch or poll.

## Rate Limiting (if applicable)

None. The profile uses SSH, tmux, filesystem state, and one terminal `sacct`
query; it calls no external HTTP API. Scheduler accounting is lifecycle-limited
and must never be converted into a polling loop.

## Common Mistakes

- Running the wrapper from `$HOME/scratch/ECODA_paper`, a login node, or a
  second tmux session instead of the Bamboo git clone and explicit gate session.
- Using bare interpreters, mutating `datasets.json`/the lockfile, passing labels
  into a method, or skipping `src/slurm_config.sh`.
- Relaunching after an SSH drop, treating event delivery as completion, or
  manually rsyncing partial outputs after a failure.
- Polling while a subagent remains alive, running more than one terminal
  `sacct`, releasing a dependent gate before inspect and Luna Max review, or
  treating `serialization_group` as dependency lineage rather than the
  required mutex identity, or ignoring a `resource_lock_conflict`
  `PRELAUNCH_STOP`.
