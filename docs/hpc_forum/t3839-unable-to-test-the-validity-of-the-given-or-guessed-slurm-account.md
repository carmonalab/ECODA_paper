# Unable to test the validity of the given or guessed SLURM account

- Source: https://hpc-community.unige.ch/t/3839

- Created: 2025-02-25T08:52:52.218Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Stephen.Mulligan (2025-02-25T08:52:52.259Z)

## Primary informations
Username: mulligas
Cluster: baobab, bamboo, yggdrasil
## Description
Using snakemake to submit jobs to slurm as part of a larger workflow but failing due to issues with slurm account
## Steps to Reproduce
submit any test job from a snakemake workflow from the slurm account “golling”
## Expected Result
the required resources should be requested and jobs submitted via slurm
## Actual Result
```
Unable to test the validity of the given or guessed SLURM account 'golling' with sacctmgr: sacctmgr: error: _open_persist_conn: failed to open persistent connection to host:lunihpcslurm1.admin.unige.ch:6819: Connection refused
sacctmgr: error: Sending PersistInit msg: Connection refused
```
not only an issue for me but also another member of my group, occurring as of this morning.


## Post 2 by @Adrien.Albert (2025-02-25T15:15:25.297Z)

Hello @Stephen.Mulligan[@Stephen.Mulligan](https://hpc-community.unige.ch/u/stephen.mulligan)
We had an issue this morning on slurm database. It have been fixed you should be able to run job again:
```
(baobab)-[root@login1~]$ su - mulligas
(baobab)-[mulligas@login1 ~]$ srun hostname
cpu001.baobab
```


## Post 3 by @Stephen.Mulligan (2025-02-25T15:17:43.147Z)

working again, thanks very much


## Post 4 by @Yann.Sagon (2025-03-14T10:16:51.164Z)

A post was merged into an existing topic: Issue with slurm[Issue with slurm](https://hpc-community.unige.ch/t/issue-with-slurm/3865/2)


## Post 5 by @Yann.Sagon (2025-03-14T10:14:52.733Z)

A post was split to a new topic: Issue with slurm[Issue with slurm](https://hpc-community.unige.ch/t/issue-with-slurm/3865)
