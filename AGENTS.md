# Agent Guardrails & Repository Conventions

> This document defines the core domain rules, architectural guardrails, and HPC execution invariants that AI agents must strictly follow in this repository.

---

## 1. Core Guardrails (Sacred Invariants)

### Paper Figures are Sacred
- Every plot saved with a filename starting with `Figure` or `Supp_fig` is a publication figure: **KEEP and FIX if broken — NEVER remove it**.
- **Benchmark Figure Hierarchy:**
  - **Main benchmark figure (`Figure 2A`):** Standard evaluation across all datasets using only the default/main parameter setting per method.
  - **Extended figure (`Supp fig 15`):** Extended non-standard methods (HiTME/scATOMIC annotations, cell-type pseudobulk, frequency-based composition).
  - **Parameter screening (`Supp fig 2`):** Main methods across parameter ranges (robustness check that default parameters are faithful).
- `ECODA_PB_combo_*` (ECODA + Pseudobulk distance combinations) are legacy/experimental and are not included in publication figures.

### No-Leakage Principle (Central Premise)
- Biological group labels (`Status`, `sample.origin`, `cond`, `Disease_Identity`, etc.) are **ground truth only**:
  - They must **NEVER** be passed as a design covariate, batch key, or input to preprocessing, DESeq2 normalization, batch correction, HVG selection, or embedding steps.
  - Batch correction must be strictly batch-only (e.g., no design protection argument in `removeBatchEffect`).
- **`DESeq2.normalize()` Semantics:**
  - **Benchmark mode (default):** `blind=TRUE`, `batch_col=NULL`, `correct_batch=FALSE` (unsupervised, design `~ 1`).
  - **Batch effect mode:** `batch_col=<col>`, `blind=FALSE`, `correct_batch=TRUE` (batch-only `limma::removeBatchEffect`).

### `datasets.json` Permissions
- `datasets.json` is the central ground truth for evaluated datasets, metadata columns, and view definitions.
- **Do not modify `datasets.json` without explicit user confirmation.**

### Reproducibility & Environments
- **Never drop defined package versions in `pixi.toml`** — preserving exact dependency versions is essential for reproducibility.
- Do not run full pipeline scripts (`.R`, `.py`, `.sh`) on large cohorts for minor validation checks unless explicitly requested by the user. Use the `_debug` dataset (Joanito 5-sample subset) for verification.

---

## 2. HPC Execution Invariants

### Cluster Access & Host Configuration
- **Cluster Endpoint:** `login1.bamboo.hpc.unige.ch` (SSH alias: `bamboo`, user: `halterc`).
- Remote commands and status checks can be executed directly via `ssh bamboo "<command>"`.

### Login Node Policy
- **Never execute heavy computation, preprocessing, or benchmarks on login nodes.**
- Login nodes are strictly for compiling, editing code, data staging (`1_stage_data.sh`), NAS synchronization, and submitting SLURM jobs/arrays.
- Long-running I/O operations (staging, NAS sync) should run inside persistent background sessions (`tmux` / `screen`).
- Use `debug-cpu` or `debug-gpu` partitions for short interactive checks.

### Data Protection & Safety
- **Never run recursive deletions (`rm -rf`)** on `$HOME/scratch` or `data/` directories without explicit user confirmation.
- **Job Monitoring:** Inspect runs with non-blocking Slurm queries (`squeue -u $USER`, `sacct -j <id> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS`).

### Repository vs. Scratch Directory Layout
- **The git repository clone lives at `$HOME/ECODA_paper`** — all git operations (`git pull`, `git status`, commits) and job submissions must run from `~/ECODA_paper`.
- **`$HOME/scratch/ECODA_paper` (`HPC_SCRATCH_DIR`) is data storage only** — it is not a git clone and does not carry tracked files.

### SLURM Submission Conventions
- **Working Directory:** All HPC bash scripts must source `src/slurm_config.sh` and set `cd "${PROJECT_ROOT}"`.
- **Spool Recovery (`BASH_SOURCE` safe):** Sbatch scripts copied to `/var/spool/slurmd/` must recover `SCRIPT_DIR` via:
  ```bash
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
      SCRIPT_DIR="$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/ {print $2}' | xargs dirname)"
  else
      SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  fi
  ```
- **Direct Environment Invocations:** Workers invoke `${PYTHON_BIN}` and `${PIXI_RSCRIPT}` (`pixi run --as-is -e py-cuda13 Rscript --vanilla`) directly. The `--as-is` flag prevents runtime lockfile modifications.
- **Environment Updates:** `setup_env_sbatch.sh` (worker node) and `refresh_env.sh` (login node) serialize via `logs/env_refresh.lock` and must never run while array jobs are active.

---

## 3. Task & Coding Conventions

- **Copy-Paste Ready Commands:** When providing commands or instructions to the user, always provide the full, copy-paste-ready commands including required environment sourcing (`source src/slurm_config.sh`).
- **Parallel Shells:** The user often operates multiple terminal sessions connected to different HPC login nodes. Highlight which commands can safely run in parallel.
- **Efficient Code Search:** Use targeted tools (`grep_search`, `git grep`). Never run un-scoped `grep -rn` across the entire workspace, as it scans the gitignored ~97 GB `data/` and `.pixi/` directories.
- **Plan Completion Workflow:** Whenever completing an implementation plan from `.kilo/plans/`:
  1. Move the completed plan to `.kilo/plans/archive/`.
  2. Stage the relevant modified files and the archived plan (`git add`).
  3. Commit and push the changes.

---

## 4. Domain Terminology

- **ECODA (Exploratory Compositional Data Analysis):** Uses CLR-transformed cell-type proportions for cohort-level patient stratification in an unsupervised setting.
- **CLR (Centered Log-Ratio):** Compositional data transformation: $\text{CLR}(x_i) = \log(x_i / g(x))$, where $g(x)$ is the geometric mean. Requires prior zero-imputation.
- **HVCs (Highly Variable Cell Types):** Cell types with the highest across-sample variance, selected for patient stratification.
- **Zero Imputation Strategies:** Count-based (`counts_zeros`, `counts_all`), percentage-based (`percentage_zeros`, `percentage_all`), and log-normal models (`multiLN`, `multiRepl` via `zCompositions`).
- **Pseudobulk:** Aggregating single-cell expression per sample followed by DESeq2 normalization.
- **Separation Metrics:** Evaluated in `src/utils/scoring_metrics.R` to quantify cluster recovery:
  - **ANOSIM:** Analysis of Similarities (`calc_sep_score()`).
  - **ARI:** Adjusted Rand Index (`clust_eval()`).
  - **Silhouette:** Silhouette width across biological classes (`calc_sil()`).
  - **Modularity:** Graph modularity across KNN graphs (`calc_modularity()`).
  - **LISI:** Local Inverse Simpson's Index (`calc_lisi()`).
- **Data Exchange Contracts:**
  - `.feather` files: Apache Arrow IPC format for cross-language distance matrices and sample embeddings.
  - `.rds` result bundles: Method output bundles saved atomically with MD5 sidecars (`checksums.md5`).
