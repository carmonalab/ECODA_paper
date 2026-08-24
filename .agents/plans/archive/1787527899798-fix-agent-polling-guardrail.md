# Configure OMP Sol/Luna Plan Orchestration

## Goal

Configure oh-my-pi (OMP) as a resumable orchestrator for large ECODA plans:

- GPT-5.6 Sol at `medium` is the default planner/orchestrator.
- GPT-5.6 Luna at `max` implements and reviews one verified step at a time.
- Workers share the current working tree and run serially.
- Commands already approved in an authoritative plan may run on Bamboo without another prompt.
- OMP waits on registered long-running work without model-side polling while its session remains alive.
- Bamboo work survives Mac sleep/restart; one manual `omp --resume` continues from durable state afterward.

Do not migrate or rename existing `.kilo` plans. OMP can read them directly.

## How The User Will Work With It

Run OMP directly, not through a Kilo agent. The recommended interface is the
VS Code integrated terminal opened at the repository root; a standalone
terminal is behaviorally equivalent. Kilo remains available for small IDE
tasks and independent review, but it must not wrap or launch the OMP agent
session.

### Existing Bassez Plan

After the one-time setup and validation in this plan are complete:

```bash
cd /Users/christianhalter/Desktop/ECODA_paper
omp
```

In the OMP TUI, confirm the normal/default model is Sol Medium, then submit:

```text
/skill:execute-plan .kilo/plans/1787522239980-bassez-annotation-benchmark-rerun.md
```

Do not enter OMP Plan mode for this run. The Bassez plan already exists, and
OMP Plan mode deliberately makes task agents read-only. The normal Sol session
is the orchestrator; it delegates writable gates to Luna.

### New Large Work

1. Start `omp` at the repository root.
2. Enter OMP's native Plan mode using the mode control shown by the installed
   TUI. During setup, inspect `/help` and record the exact installed shortcut;
   do not assume a Kilo slash command maps to OMP.
3. Confirm the active `plan` role is `openai-codex/gpt-5.6-sol:medium`.
4. Ask Sol to investigate, interview, and produce an implementation-ready
   plan. Read-only planning subagents may assist.
5. In OMP's plan-review UI, save the approved plan under
   `.agents/plans/<topic>.md` and choose the action that starts a fresh normal
   session.
6. In that fresh session, confirm `default` is Sol Medium and invoke
   `/skill:execute-plan .agents/plans/<topic>.md`.
7. Let Sol coordinate serial Luna implementation/review gates. Use `Alt+A` only
   when human inspection or steering of a child is needed; routine completion
   is delivered to the parent automatically.

Avoid `--plan-yolo` for publication/HPC plans: it bypasses the explicit plan
review boundary and selects an implementation target rather than enforcing the
Sol-orchestrator/Luna-worker loop.

### Resume After Interruption

From the repository root, run `omp --resume`, select the saved session, and
invoke the same execute-plan skill/path if the session is idle. The skill must
reconcile the plan, git state, remote tmux/status logs, and Slurm accounting
before it delegates another gate. A completed remote job is never sufficient
evidence by itself.

### Small Direct Tasks

Use normal OMP with Sol when coordination or analysis matters. For a narrow
task that does not need a planner, launch OMP with the configured `@worker`
role and let Luna work directly. Do not use the execute-plan skill for trivial
single-step edits.

## Decisions And Boundaries

- Use a hybrid configuration: reusable model roles, agents, and execute-plan skill under `~/.omp/agent`; ECODA-specific concurrency and safety settings under project `.omp/`.
- Keep Kilo installed and its config unchanged. OMP is an additional execution harness, not a repository-format migration.
- Use `openai-codex/...` selectors because the installed OMP config already authenticates through that provider. Verify both requested models before changing defaults.
- Do not use OMP read-only Plan mode to execute an existing plan: plan-mode subagents cannot edit. Run the main Sol session normally and invoke the execute-plan skill.
- Use the current tree, not OMP isolation. OMP branch isolation can stash/merge and conflicts with this repository's no-stash convention and sequential shared state.
- Set project task concurrency to one. The Bassez benchmark wrappers share outputs and must remain serialized.
- The Sol parent may update plan status/todos but must not implement source changes. Luna workers edit/run; Luna reviewers are read-only.
- OMP child sessions are headless and force tool approval to `yolo`; enforce destructive-command denials in project config and retain all `AGENTS.md` safety rules.
- Automatic post-reboot agent invocation is out of scope. It would require `launchd`, durable session routing, locking, VPN/NAS readiness, and unattended credentials. Initial recovery is one manual resume command.
- The implementation of this setup plan ends after configuration and smoke
  validation. It must not immediately launch the production Bassez workflow;
  the user starts that workflow from OMP using the command above after reviewing
  the setup validation results.

## 1. Fix The Portable Background-Work Rule

Replace the current `AGENTS.md` strict launch-and-stop block before testing OMP. Keep these semantics:

- Never issue repeated model-side status calls, timer calls, or sleep loops.
- A bounded directly tracked command may run normally and return its result.
- A registered background task may use exactly one documented blocking/event-driven wait call.
- Detached `tmux -d`, `nohup`, and `sbatch` launch success means submission only, not workload completion.
- For external detached work without a managed handle, report host, session/job ID, log/status path, and one bounded inspection command, then yield without promising automatic wake-up.
- Pipeline scripts may perform their existing bounded/fail-closed scheduler monitoring; that is orchestration inside one managed process, not model tool-call polling.

Do not mention a nonexistent global `kilo.terminalTimeout`, and do not use `tail -f` as a one-shot status check.

## 2. Configure Global OMP Model Roles

Patch `~/.omp/agent/config.yml` without overwriting the existing theme, composer, symbol, or setup settings:

```yaml
modelRoles:
  default: openai-codex/gpt-5.6-sol:medium
  plan: openai-codex/gpt-5.6-sol:medium
  worker: openai-codex/gpt-5.6-luna:max
  review: openai-codex/gpt-5.6-luna:max
  task: openai-codex/gpt-5.6-luna:max

task:
  agentModelOverrides:
    task: "@worker"
    reviewer: "@review"

async:
  enabled: true
  pollWaitDuration: smart

bash:
  autoBackground:
    enabled: true
    thresholdMs: 60000

launch:
  enabled: true

skills:
  enableSkillCommands: true
```

Before patching, confirm OMP 18.0.3 or newer and verify that `openai-codex/gpt-5.6-sol` and `openai-codex/gpt-5.6-luna` are visible. If either selector or effort suffix is rejected, stop and use OMP's `/model` Roles UI to select the exact discovered ID; do not silently substitute another model.

## 3. Add Reusable Luna Agents

Create `~/.omp/agent/agents/luna-worker.md`:

- `name: luna-worker`, `model: "@worker"`, `blocking: true`.
- Tools: `read`, `grep`, `glob`, `edit`, `write`, `bash`, `lsp`, `hub`, and `yield`; no child spawning.
- Operate only on the assigned plan step in the current tree.
- Read root `AGENTS.md`, authoritative plan status, and current repo/HPC state before acting.
- Never rerun a completed or ambiguous step; inspect evidence once and return a blocker when state cannot be proven.
- For an approved long command, use managed async Bash with `timeout: 0`, then one `hub wait` for its job ID with `timeoutMs: 0`.
- Return a compact structured result: status, step/gate, files changed, commands/jobs/session IDs, checks/evidence, plan status change, blockers, and next safe state.

Create `~/.omp/agent/agents/luna-reviewer.md`:

- `name: luna-reviewer`, `model: "@review"`, `blocking: true`.
- Read-only tools: `read`, `grep`, `glob`, `lsp`, bounded/non-mutating `bash`, `hub`, and `yield`.
- Verify the assigned step against the plan, diff, checks, and HPC evidence without editing or rerunning heavy pipelines.
- Lead with concrete findings and return `PASS`, `FIX`, or `BLOCKED`, with exact evidence and the next permitted action.

Give both agents an inline JSON output schema in frontmatter so malformed or verbose returns fail validation instead of polluting the Sol parent context.

## 4. Add The Global Execute-Plan Skill

Create `~/.omp/agent/skills/execute-plan/SKILL.md` with an explicit name and description. Its workflow must instruct the Sol parent to:

1. Read the supplied plan path and treat its authoritative status plus verified external state as the source of truth.
2. Inspect repo status and any referenced Bamboo/tmux/Slurm state before creating work; never infer that an old launch completed.
3. Build/update an OMP todo list at gate granularity, not one item per prose bullet.
4. Select only the next dependency-ready gate.
5. Spawn exactly one `luna-worker` with a self-contained assignment and required evidence.
6. Wait for the blocking worker result; then spawn exactly one `luna-reviewer` for that result.
7. On `PASS`, update only the plan's authoritative status/checkpoint and advance. On `FIX`, message the same worker when its context is useful or issue one focused repair task. On `BLOCKED`, stop safely.
8. Keep raw logs/diffs in files or agent artifacts; retain only compact summaries in the parent context.
9. Continue until the plan's completion boundary, then follow its archive/commit/push instructions exactly.
10. On resumed sessions, reconcile the plan, git state, remote status/logs, tmux sessions, and Slurm accounting before deciding the next gate.

The skill must prohibit the Sol parent from editing source code or directly running heavy pipelines. It may edit the plan status and coordinate agents.

## 5. Add ECODA Project Settings

Create `.omp/config.yml` in the repository:

```yaml
task:
  maxConcurrency: 1
  maxRuntimeMs: 0
  isolation:
    mode: none

async:
  enabled: true
  pollWaitDuration: smart

bash:
  autoBackground:
    enabled: true
    thresholdMs: 60000
  patterns:
    - match: "rm -rf *"
      approval: deny
    - match: "rm --recursive *"
      approval: deny
    - match: "git reset --hard*"
      approval: deny
    - match: "git checkout -- *"
      approval: deny
    - match: "git stash*"
      approval: deny
```

Add a final broad rule only if OMP's schema requires one; do not accidentally deny normal commands after the listed destructive patterns. Do not create `.omp/AGENTS.md`: the existing root `AGENTS.md` must remain the sole project instruction source rather than being shadowed by a native OMP context file.

## 6. Define The Bamboo Long-Task Pattern

For approved long login-node wrappers or NAS transfers:

1. Use a unique remote tmux name, log path, and durable status path derived from the plan/run identifier.
2. Run the workload inside remote tmux, but keep an attached SSH client (`ssh -tt ... tmux new-session -A -s ...`) as the process OMP tracks.
3. Launch that SSH command through OMP managed async Bash with `timeout: 0`.
4. Have the Luna worker call `hub wait` once with the returned job ID and `timeoutMs: 0`; do not issue repeated `hub jobs`, `squeue`, `sacct`, or sleep calls.
5. On completion, independently verify the durable status, expected artifacts, checksums, and/or `sacct` terminal state before advancing. Do not trust tmux/SSH exit alone.
6. If the Mac sleeps or restarts, remote tmux/Slurm continues. After restart, resume OMP and inspect the named tmux session, log/status file, and Slurm state once before relaunching anything.
7. If SSH drops while the remote command continues, mark local completion unknown. Reattach or use the plan's documented recovery path (`--sync-only`, watchdog ID, or one-shot artifact gate).

For Slurm-heavy multi-stage work, prefer the repository's existing wrappers and watchdog/dependency logic. OMP should supervise one wrapper process, not reproduce Slurm polling in model turns.

## 7. Validate Before Production HPC Use

Run these gates in order:

1. **Config/model gate:** launch OMP in the repo, inspect `/model` Roles, and require Sol Medium for `default`/`plan` and Luna Max for `worker`/`review`/`task`.
2. **Context gate:** ask OMP to summarize sacred figure, no-leakage, dataset permission, and HPC login-node rules; require answers from root `AGENTS.md`.
3. **Routing gate:** invoke `luna-worker` on a read-only repository inspection and confirm the child reports Luna Max while the parent remains Sol Medium.
4. **Review gate:** invoke `luna-reviewer` on a harmless known diff/status and confirm it cannot edit.
5. **Managed-wait gate:** have a worker run one harmless, bounded Bamboo status command as an explicit async Bash job and use one indefinite `hub wait`; require automatic result delivery without model polling.
6. **Persistence gate:** use a harmless uniquely named remote tmux command that writes a log/status marker, track it through attached SSH, and verify the marker once. Do not use a production pipeline for this test.
7. **Resume gate:** exit and resume the OMP session, invoke the skill again, and confirm it reconciles existing status rather than repeating the test.
8. **Operator UX gate:** use the installed OMP `/help`/mode UI to identify the
   exact Plan-mode control, plan-review save/start-new-session action, model-role
   display, Agent Hub control, and resume picker. Return these exact controls in
   the setup completion summary rather than relying on guessed Kilo-equivalent
   commands.

Stop rollout on wrong model routing, missing `AGENTS.md`, parallel mutating workers, inability to enforce destructive-command denials, or failure to recover a persisted remote status.

## 8. Resume The Existing Bassez Plan

Use `.kilo/plans/1787522239980-bassez-annotation-benchmark-rerun.md` unchanged as the authoritative execution plan.

- Its status says source hardening and jobs `4334702`-`4334704` are complete; do not repeat them.
- The first candidate gate is Section 3's dataset-output NAS persistence, but verify current scratch/NAS hashes and any existing sync session/log before launching it.
- Keep all four benchmark wrappers serialized as required by the plan and project `maxConcurrency: 1`.
- Require worker plus reviewer evidence before updating the authoritative status and advancing to the next wrapper.
- Use the plan's watchdog IDs, `--sync-only` recovery, artifact audits, local mirror refresh, local notebook render, sacred-figure checks, and final archive/commit/push boundary exactly as written.

After setup validation, the user starts a normal OMP session from the
repository root and invokes:

```text
/skill:execute-plan .kilo/plans/1787522239980-bassez-annotation-benchmark-rerun.md
```

After a Mac restart, return to the repository and use OMP's session resume picker (or the saved session ID), then invoke the same skill/path. The resumed orchestrator must reconcile durable state before continuing.

## Out Of Scope

- Automatic `launchd`-driven OMP restart after Mac reboot.
- A webhook/mail listener that injects Bamboo completion into a closed conversation.
- Moving notebook rendering to Bamboo.
- Changing pipeline code solely to support OMP unless the smoke test exposes a concrete missing durable-status contract.
- Modifying `datasets.json`, package versions, publication scope, or sacred figure names.
