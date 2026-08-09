# gloscope benchmark job failed: missing R packages HiTME, MOFA2

**Session ID:** ses_022e1cde5ffe8Ox5eOaOz1rFX4
**Created:** 07/08/2026, 18:46:11
**Updated:** 07/08/2026, 18:46:48

---

## User

what's wrong and how to fix? if it's just a bug due to the small sample size, it's fine, I'll test it later with a real (bigger) dataset.


(bamboo)-[halterc@login1 ECODA_paper]$ ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh   --ds_name _debug --methods gloscope,mofa,pseudobulk
Found 1 benchmark datasets.
NOTE: mofa/pseudobulk requested without prepare_pseudobulk; auto-prepending it.
=== Submitting R benchmark method arrays ===
Submitting prepare_pseudobulk array (1 datasets, partition=shared-cpu, 
  flags: --constraint=EPYC-7742 --cpus-per-task=16, mem=128G, throttle=500)
  prepare_pseudobulk array job ID: 4291184
=== Monitoring prepare_pseudobulk array 4291184 ===
prepare_pseudobulk array 4291184 finished. Checking task states...
Array Job 4291184 (prepare_pseudobulk): all tasks COMPLETED.
Submitting gloscope array (1 datasets, partition=shared-cpu, 
  flags: --constraint=EPYC-7742 --cpus-per-task=16, mem=128G, throttle=500)
  gloscope array job ID: 4291189
Submitting mofa array (1 datasets, partition=shared-cpu, 
  flags: --constraint=EPYC-7742 --cpus-per-task=16, mem=128G, throttle=500)
  mofa array job ID: 4291190
Submitting pseudobulk array (1 datasets, partition=shared-cpu, 
  flags: --constraint=EPYC-7742 --cpus-per-task=16, mem=128G, throttle=500)
  pseudobulk array job ID: 4291191
=== Monitoring benchmark array 4291189 ===
benchmark array 4291189 finished. Checking task states...
ERROR: Array Job 4291189 had non-COMPLETED tasks; NOT syncing to NAS.
JobID           JobName      State ExitCode 
------------ ---------- ---------- -------- 
4291189_1    benchmark+     FAILED      1:0 
4291189_1.b+      batch     FAILED      1:0 
4291189_1.e+     extern  COMPLETED      0:0 
(bamboo)-[halterc@login1 ECODA_paper]$ cat logs/5_benchmark_r_gloscope_4291189_1.err
 WARN cache for Repodata at [REDACTED_PATH]/.cache/rattler/cache/repodata is on a network/parallel filesystem (NFS/SMB/FUSE/BeeGFS/Lustre/GPFS/CephFS), redirected to /tmp/pixi-cache-halterc/repodata for this run. Set [cache.repodata] in config.toml or PIXI_CACHE_DIR to override, or [cache.netfs-redirect] = "never" to keep the original path.
Error in load_my_packages(my_packages) : 
  The following packages are missing from the pixi environment: 
HiTME, MOFA2

Please add them to your pixi.toml (e.g., as 'r-packagename') and run `pixi install`.
Calls: source ... source -> withVisible -> eval -> eval -> load_my_packages
In addition: There were 16 warnings (use warnings() to see them)
Execution halted

---

## Assistant

**Tool: grep**

**Tool: grep**

---

## Assistant

**Tool: read**

**Tool: read**

---

## Assistant

**Tool: bash**

**Tool: grep**

---

## Assistant

**Tool: read**

---

## Assistant

Not a sample-size bug — the job dies before touching any data, during package loading. Here's the diagnosis:

## What's wrong

The gloscope worker sources `load_all_functions.R` → `src/utils/imports.R`, whose `my_packages` list **hard-requires `HiTME` and `MOFA2` for every R script** (lines 21, 28). Those two packages are **not conda packages** — they're installed by the pixi `setup` task in `pixi.toml`:

- `HiTME` → `remotes::install_github('carmonalab/HiTME@87f8d414...')` (line 131)
- `MOFA2` → `remotes::install_version('MOFA2', version='1.20.2', repos=bioc_repos)` (line 147)

The last pixi.toml change is commit `3805890` ("Add R benchmark methods HPC pipeline..."), dated **today 18:13**. So your HPC clone's R library predates these additions — the `setup` task was never re-run there after pulling, and the worker fails at `require(HiTME)` / `require(MOFA2)` (exit code 1:0, nothing computed).

## How to fix

On the HPC **login node** (installs are allowed there), in the repo clone:

```bash
git pull                       # make sure pixi.toml is current
pixi install -e py-cuda13      # sync conda env (README HPC note, lines 73-77)
pixi run -e py-cuda13 setup    # installs the missing R packages
```

The setup task is idempotent (`requireNamespace` guards + remotes SHA check), so it only installs what's missing — but note it builds MOFA2 from Bioconductor source and pulls HiTME's deps (SignatuR/STACAS/ProjecTILs), so it may take a while. Then resubmit:

```bash
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name _debug --methods gloscope,mofa,pseudobulk
```

## Heads-up

The same fix is required **before any annotation array run** too — `src/4_cell_type_annotation/2.1.1_process_chunk.R` also does `library(HiTME)` (line 21). The error message's suggestion ("add to pixi.toml as 'r-packagename'") is misleading — those packages can't be `r-` conda deps; the setup task is the right mechanism.

---

