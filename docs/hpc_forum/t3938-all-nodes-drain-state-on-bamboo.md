# All nodes drain state on bamboo

- Source: https://hpc-community.unige.ch/t/3938

- Created: 2025-04-24T11:34:31.378Z

- Tags: bamboo

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Raphael.Rubino (2025-04-24T11:34:31.425Z)

## Primary informations
Username: rubinor
Cluster: bamboo
## Description
All nodes are in drain state on bamboo.
## Steps to Reproduce
`sinfo` on bamboo
image
image594×260 31 KB
[image594×260 31 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/880499ec82cdf5fcda84f02ce5f21e6a622d1aa7.png)


## Post 2 by @Yann.Sagon (2025-04-24T12:20:47.073Z)

Dear @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino)
this is fixed. [2025] Current issues on HPC Cluster - #10 by Yann.Sagon[[2025] Current issues on HPC Cluster - #10 by Yann.Sagon](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/10)


## Post 3 by @Raphael.Rubino (2025-04-25T09:03:23.474Z)

Thanks @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
Still on bamboo, gpu003 is currently unusable.
`scontrol` shows
`   Reason=gres/gpu GRES core specification 12-15 doesn't match socket boundaries. (Socket 0 is cores 0-16) [slurm@2025-04-25T10:52:59]`


## Post 4 by @Yann.Sagon (2025-04-25T09:08:35.827Z)

Yes, fixed. This is an annoying bug that will be fixed in the next Slurm release.
Best


## Post 5 by @Raphael.Rubino (2025-04-25T14:49:01.018Z)

Thanks @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
It seems that bamboo login node is unreachable outside unige network at the moment.
Best regards


## Post 6 by @Yann.Sagon (2025-04-25T15:12:02.773Z)

This is fixed. Time for the week-end!
