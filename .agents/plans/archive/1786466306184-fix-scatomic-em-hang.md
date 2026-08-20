# Fix: scATOMIC/cutoff.scATOMIC `em()` infinite EM loop hangs annotation tasks

## Context (root cause, confirmed by source + logs)

- `scATOMIC::automatic_threshold()` (pinned `d332cd5`) calls `cutoff.scATOMIC::em(score, "normal", "normal")` to fit a two-component normal mixture for per-cell-type confidence cutoffs.
- `cutoff.scATOMIC::em()` (fork pinned `a43fd5e1`, `R/em.R:111-143`) iterates `while (abs(lambda0 - mean(lambda)) > t)` with default `t = 1e-64` — a tolerance **below double-precision machine epsilon (~2.2e-16)** — and **no iteration cap**. The loop can only exit when two consecutive EM iterations produce bit-identical `mean(lambda)`.
- `automatic_threshold()` **jitters tail scores with `runif()`** before calling `em()`. On ill-conditioned score mixtures the M-step (`bbmle::mle2(mLL2, start, "Nelder-Mead")`, which repeatedly hits `maxit` and returns non-converged fits — the warning seen in the logs) yields different coefficients every iteration, so `mean(lambda)` never bit-stabilizes → the loop runs effectively forever, each iteration doing a full Nelder-Mead fit over up to ~10k cells.
- The worker's `withTimeout` interrupt fires (`Error in namedrop(args) : reached elapsed time limit` — `namedrop` is only the frame where the interrupt landed; `getAnywhere` confirms it is not in scATOMIC/cutoff.scATOMIC) but is swallowed inside the package stack, so the attempt continues until Slurm's 2h wall limit kills the task (no transient signature → correctly no requeue).
- Evidence: 2/512 tasks in array 4309053 hit "reached elapsed time limit"; the same sample (Pelka `C171_TA`, 10,709 cells, chunk_48 / task 378) hung identically in a previous array (log 4308096_378.err shows the same `mLL2`/`maxit` warnings + wall-limit kill). Deterministic for this sample; possible for any sample with an oscillating mixture fit.

## Decisions

- **Fix at worker level**: monkey-patch `cutoff.scATOMIC::em` in `2.1.1_process_chunk.R` right after `library(cutoff.scATOMIC)`. Do NOT fork/repin dependencies (keeps reproducibility story unchanged).
- The patch must override the **namespace binding** (`automatic_threshold` calls `cutoff.scATOMIC::em` explicitly via `::`): `unlockBinding` + `assignInNamespace`, and the patched function's `environment()` must be set to `asNamespace("cutoff.scATOMIC")` so internal helpers (`hash`, `startval`, `mLL`) resolve.
- Bounded loop: default `t = 1e-4`, cap `max_iter = 200`, `message()` on cap (goes to stderr → `.err` for observability). Numerically the threshold values are unchanged for well-behaved classes (1e-64 was never reachable in double precision).
- **Defer** the process-level watchdog for scATOMIC (subprocess + shell `timeout`) — out of scope; `withTimeout` + the 2h wall limit remain backstops.
- Docs: add a brief note to `docs/ARCHITECTURE.md` (annotation section) and one line to `AGENTS.md` (annotation bullet) referencing the upstream issue.
- GitHub issue: created manually by the user (gh CLI is not authenticated; plan mode cannot post). Subject + body ready in section "GitHub issue text".

## Tasks (ordered)

1. **Patch `src/4_cell_type_annotation/2.1.1_process_chunk.R`**: insert the patch block immediately after `library(cutoff.scATOMIC)` (line ~24). Use the code sketch below. Safe to edit while array 4309053 runs (R reads the script at task launch; a Slurm requeue would pick up the patch, which is fine).
2. **Docs**: `docs/ARCHITECTURE.md` annotation-pipeline note (worker monkey-patches `cutoff.scATOMIC::em` — bounded EM loop, upstream bug, issue link) + one line in `AGENTS.md` annotation bullet.
3. **Validation (HPC)**:
   a. Smoke test — run with the patch applied:
      ```bash
      ~/ECODA_paper/.pixi/envs/py-cuda13/bin/Rscript --vanilla -e '
      library(cutoff.scATOMIC)
      # <apply patch snippet from 2.1.1_process_chunk.R>
      set.seed(1)
      x <- c(runif(300, 0.9, 1), runif(300, 0.01, 0.1), rnorm(400, 0.5, 0.2))
      t0 <- Sys.time(); res <- cutoff.scATOMIC::em(x, "normal", "normal"); print(Sys.time() - t0); print(res$param)
      '
      ```
      Must return in seconds with the "EM loop capped" message (pre-patch it would hang).
   b. After array 4309053 leaves the scheduler: `./src/4_cell_type_annotation/2_submit_hpc_array.sh Pelka` → 511 tasks skip via the existing-feather check; chunk_48 must complete. Confirm `annotations_chunk_48.feather` exists and "chunk processing complete" appears in the new `.log`.
   c. Confirm the second affected task's chunk (the other file matching "reached elapsed time limit") also completes in the re-run.
   d. Then the standard post-run flow: sacct (once slurmdbd is back) → `./src/4_cell_type_annotation/2_submit_hpc_array.sh "" --sync-only 4309053` → `3_submit_merge.sh`.
4. **GitHub issue**: user creates on https://github.com/abelson-lab/scATOMIC using the text below. Optionally also file on https://github.com/inofechm/cutoff.scATOMIC (the fork owning `em()`).
5. **Task completion workflow** (AGENTS.md): archive this plan to `.kilo/plans/archive/`, `git add .`, commit, push.

## Code sketch (patch block for 2.1.1_process_chunk.R)

```r
# Bounded cutoff.scATOMIC::em(): upstream bug — EM loop with default t=1e-64
# (below machine epsilon) and no iteration cap never terminates on
# ill-conditioned score mixtures (automatic_threshold jitters scores with
# runif, so the mle2 M-step never bit-stabilizes). See abelson-lab/scATOMIC
# issue <NUMBER>.
em_bounded <- function(data, D1, D2, t = 1e-4, max_iter = 200) {
  data_name <- unlist(strsplit(deparse(match.call()), "="))[2]
  data_name <- sub(",.*$", "", gsub(" ", "", data_name))
  start <- as.list(startval(data, D1, D2))
  D1b <- hash[[D1]]
  D2b <- hash[[D2]]
  lambda0 <- 0
  iter <- 0
  with(start, {
    while (abs(lambda0 - mean(lambda)) > t) {
      iter <- iter + 1
      if (iter >= max_iter) {
        message(sprintf("cutoff.scATOMIC::em: EM loop capped at %d iterations (t=%g)", max_iter, t))
        break
      }
      lambda <- mean(lambda)
      lambda0 <- lambda
      distr1 <- lambda * D1b(data, mu1, sigma1)
      distr2 <- (1 - lambda) * D2b(data, mu2, sigma2)
      lambda <- distr1 / (distr1 + distr2)
      mLL2 <- function(mu1, sigma1, mu2, sigma2)
        return(mLL(mu1, sigma1, mu2, sigma2, lambda, data, D1b, D2b))
      start <- as.list(log(c(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2)))
      out <- bbmle::mle2(mLL2, start, "Nelder-Mead")
      coef <- out@coef
      coef_n <- names(coef)
      names(coef) <- NULL
      for (i in 1:4) assign(coef_n[i], exp(coef[i]))
    }
    out <- list(lambda = lambda, param = exp(out@coef), D1 = D1, D2 = D2,
                deviance = out@min, data = data, data_name = data_name, out = out, t = t)
    class(out) <- "em"
    return(out)
  })
}
em_ns <- asNamespace("cutoff.scATOMIC")
environment(em_bounded) <- em_ns
unlockBinding("em", em_ns)
assignInNamespace("em", em_bounded, em_ns)
lockBinding("em", em_ns)
message("Patched cutoff.scATOMIC::em with bounded EM loop (t=1e-4, max_iter=200)")
```

## GitHub issue text

**Subject**: `run_scATOMIC() hangs indefinitely on some samples: cutoff.scATOMIC::em() EM loop never terminates (t=1e-64, no iteration cap)`

**Body**:

```
## Summary
run_scATOMIC() occasionally never returns on certain samples: it consumes 100% CPU for hours without crashing, and only ends when the job's wall-clock limit kills the process. We traced the hang to the confidence-threshold step automatic_threshold() -> cutoff.scATOMIC::em().

## Environment
- scATOMIC @ d332cd5cf6a1ecef7c32d0adc4a862a4c47bcd95
- cutoff.scATOMIC fork @ a43fd5e1e8f0e3b71ec446970e4316f305939a17
- R 4.x, R.utils 2.13.0 (withTimeout wrapper), bbmle 2.x
- Reproduced deterministically on the same sample (colorectal cancer, ~10.7k cells) in two independent runs; intermittent on other samples.

## Symptom / log signature
  scATOMIC attempt 1 with 643 s timeout
  ... MAGIC warnings ...
  Error in namedrop(args) : reached elapsed time limit
  In addition: Warning messages:
  1: In asMethod(object) : sparse->dense coercion: allocating vector of size 1.7 GiB
  2: In bbmle::mle2(mLL2, start, "Nelder-Mead") : convergence failure: code=1 (iteration limit 'maxit' reached)

A wrapping withTimeout (setTimeLimit) fires, but the interrupt is swallowed inside the package's own error handling, so the call never returns and the process runs until the scheduler kills it at the wall limit. (The namedrop frame is not in scATOMIC or cutoff.scATOMIC; it is simply where the interrupt landed.)

## Root cause
1. cutoff.scATOMIC::em() (R/em.R) iterates
      while (abs(lambda0 - mean(lambda)) > t)
   with default t = 1e-64, a tolerance far below double-precision machine epsilon (~2.2e-16). The loop can only exit when two consecutive EM iterations produce bit-identical mean(lambda), which in practice never happens on real data.
2. automatic_threshold() (R/automatic_threshold.R) replaces tail scores with random jitter before calling em():
      score[index_1] <- runif(length(index_1), min = 0.9, max = 1)
      score[index_2] <- runif(length(index_2), min = 0.01, max = 0.1)
      mixed_model <- cutoff.scATOMIC::em(score, "normal", "normal")
   Because the input is randomized, the M-step (bbmle::mle2 Nelder-Mead, which here repeatedly hits the maxit limit and returns non-converged fits) yields slightly different coefficients every iteration, so mean(lambda) never bit-stabilizes and the while loop never terminates.
3. Each iteration runs a full 4-parameter Nelder-Mead fit over the per-class score vector (up to ~10k cells), so even a few hundred iterations take hours.

## Proposed fix
- cutoff.scATOMIC::em(): use a sane tolerance (e.g., t = 1e-4 or sqrt(.Machine$double.eps)) and add an iteration cap that warns on exit (e.g., max 200 iterations).
- automatic_threshold(): make the pre-em() score adjustment deterministic (fixed offsets instead of runif) so the fit is reproducible; optionally pass control = list(maxit = ...) to bbmle::mle2.

## Workaround (users)
Monkey-patch em() in the calling session: redefine it with a bounded loop (sane t + iteration cap) and install it into the cutoff.scATOMIC namespace via unlockBinding() + assignInNamespace(); the patched function must have environment set to asNamespace("cutoff.scATOMIC") so the internal helpers (hash, startval, mLL) resolve. Happy to share the snippet.
```

## Risks / notes

- The patched `em()` returns the last MLE estimates when the cap is hit — strictly better than never returning; thresholds for well-behaved classes are unchanged (1e-64 was never reachable).
- `lockBinding` after assignment keeps the namespace sealed for the rest of the session; the patch is applied once per worker R session.
- Editing `2.1.1_process_chunk.R` while the current array still runs is safe (no task re-reads it mid-run); do not touch `config_helper.R`/`datasets.json` (unchanged by this plan).
- Out of scope: process-level watchdog for scATOMIC; env/dependency changes; changes to `2.1_run_worker.sh`/`worker_retry.sh`.
