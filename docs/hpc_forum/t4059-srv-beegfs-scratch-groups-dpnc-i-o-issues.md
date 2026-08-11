# /srv/beegfs/scratch/groups/dpnc/ I/O issues

- Source: https://hpc-community.unige.ch/t/4059

- Created: 2025-08-22T13:46:37.556Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Vilius.Cepaitis (2025-08-22T13:46:37.614Z)

Dear HPC team,
I’m having issues accessing shared scratch storage on baobab. Could you please investigate this?
Thank you very much in advance.
Best regards,
Vilius
```
➜ ~ cp /srv/beegfs/scratch/groups/dpnc/atlas/pileup/images/weakly-supervised-search_latest.sif ~/

cp: overwrite '/home/users/c/cepaitis/weakly-supervised-search_latest.sif'? y
cp: error reading '/srv/beegfs/scratch/groups/dpnc/atlas/pileup/images/weakly-supervised-search_latest.sif': Remote I/O error
```


## Post 2 by @Adrien.Albert (2025-08-22T16:45:04.299Z)

Hello,
The issue is related to [2025] Current issues on HPC Cluster - #20 by Yann.Sagon[[2025] Current issues on HPC Cluster - #20 by Yann.Sagon](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/20) and has been resolved.
