# Snakemake workflows unable to connect to network

- Source: https://hpc-community.unige.ch/t/3491

- Created: 2024-06-14T09:16:17.858Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Debajyoti.Sengupta (2024-06-14T09:16:17.905Z)

I am trying to run a snakemake workflow with the following slurm profile
```
verbose: false
printshellcmds: true
rerun-triggers: mtime
jobs: 50

use-apptainer: true
apptainer-args: >-
  --bind /srv,/home

executor: slurm
default-resources:
  slurm_account: golling
  slurm_partition: shared-cpu,private-dpnc-cpu,public-cpu
  runtime: 60
  threads: 4
  mem_mb: 20000
```
When I run a workflow with `snakemake -s workflow/windowwidth.smk --profile workflow/profiles/slurm/`, I get this error:
image
image2771×669 34.6 KB
[image2771×669 34.6 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/6164b0854a7b7905b7eee885db5983b21b552720.png)
```
Unable to test the validity of the given or guessed SLURM account 'golling' with sacctmgr: sacctmgr: error: slurm_persist_conn_open_without_init: failed to open persistent connection to host:lunihpcslurm1.admin.unige.ch:6819: Network is unreachable
sacctmgr: error: Sending PersistInit msg: Network is unreachable
```
I am not sure how to debug this. Looks like sacct is not working? Any clues?
Cheers,
Debajyoti.


## Post 2 by @Yann.Sagon (2024-06-14T09:43:53.336Z)

Hi, this is solved.
Best
