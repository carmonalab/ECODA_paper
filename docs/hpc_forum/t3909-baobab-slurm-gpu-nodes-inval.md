# [baobab][slurm] gpu nodes inval

- Source: https://hpc-community.unige.ch/t/3909

- Created: 2025-04-03T19:34:31.807Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2025-04-03T19:34:31.853Z)

## Primary informations
Username: falkiewi
Cluster: baobab
## Description
Dear @support[@support](https://hpc-community.unige.ch/groups/support),
There seems to be a problem with most of the GPU nodes on the cluster
```
(baobab)-[falkiewi@login1 ~]$ sinfo -R
REASON               USER      TIMESTAMP           NODELIST
issue-5894           root      2025-04-02T16:01:38 cpu242
issue-6543           root      2025-03-24T09:32:28 cpu319
issue-5958           root      2024-12-02T09:39:22 cpu246
gres/gpu GRES core s slurm     2025-04-03T15:48:07 gpu030
gres/gpu count repor slurm     2025-04-02T15:47:44 gpu002
gres/gpu GRES core s slurm     2025-04-03T15:48:07 gpu[004-009]
gres/gpu count repor slurm     2025-04-03T15:48:07 gpu010
gres/gpu GRES core s slurm     2025-04-02T15:47:44 gpu011
gres/gpu GRES core s slurm     2025-04-03T15:48:07 gpu[013-017,021-026,034-044,046-049]
gres/gpu GRES core s slurm     2025-04-03T15:48:07 gpu[018-020,027-029,031-033,045]
```
I kindly ask for a quick solution to the problem and thank you in advance!
Best regards,
Maciej Falkiewicz


## Post 2 by @Yann.Sagon (2025-04-04T06:35:29.526Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
maciej.falkiewicz:
> There seems to be a problem with most of the GPU nodes on the cluster
This is a bug in slurmd: 22498 – GRES cores doesn't match socket boundaries[22498 – GRES cores doesn't match socket boundaries](https://support.schedmd.com/show_bug.cgi?id=22498)
In the meantime, I’ll resume the nodes manually.
