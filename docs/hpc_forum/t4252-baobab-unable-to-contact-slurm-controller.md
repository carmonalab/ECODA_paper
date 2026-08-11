# Baobab: unable to contact slurm controller

- Source: https://hpc-community.unige.ch/t/4252

- Created: 2026-03-17T13:29:20.233Z

- Tags: baobab, slurm

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2026-03-17T13:29:20.314Z)

## Primary informations
Username: coppinp
Cluster: Baobab
## Description
The slurm controller does not seem to be accessible:
> $ squeue --me
> slurm_load_jobs error: Unable to contact slurm controller (connect failure)
> $ sinfo
> slurm_load_partitions: Unable to contact slurm controller (connect failure)
Would it be possible to have a look?


## Post 2 by @Paul.Coppin (2026-03-17T13:48:18.702Z)

Seems to be back online now.
It may have just been a temporary issue.
