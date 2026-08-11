# [Storage][Bamboo]Symlink on scratch missing

- Source: https://hpc-community.unige.ch/t/3536

- Created: 2024-07-12T11:42:45.400Z

- Tags: bamboo

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-07-12T11:42:45.487Z)

Dear HPC users,
We have noticed the symlink `$HOME/scratch` was missing on bamboo.
It have been fixed.
```
(bamboo)-[root@login1 users]$ ll /home/users/a/alberta/
total 358
drwxr-xr-x 3 alberta hpc_users      1 Jun 26 11:18 ondemand
lrwxr-xr-x 1 root    root          35 Jul 12 13:14 scratch -> /srv/beegfs/scratch/users/a/alberta
```
