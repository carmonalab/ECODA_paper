# Private-gap-cpu nodes down/drained

- Source: https://hpc-community.unige.ch/t/3799

- Created: 2025-01-27T14:07:27.767Z

- Tags: bamboo

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Laure.Moinat (2025-01-27T14:07:27.812Z)

Username: moinatl
Custer: Bamboo
Dear team,
The nodes for the ‘private-gap-cpu’ partition seemed to be drained or down, is there a reason ? I have the impression that its the case since the end of December.
Best,
Laure


## Post 2 by @Yann.Sagon (2025-01-28T08:46:08.667Z)

Dear @Laure.Moinat[@Laure.Moinat](https://hpc-community.unige.ch/u/laure.moinat)
The partition is working and used:
```
(bamboo)-[root@admin1 ~]$ sinfo -p private-gap-cpu
PARTITION       AVAIL  TIMELIMIT  NODES  STATE NODELIST
private-gap-cpu    up 7-00:00:00      1    mix cpu007
```
What is the issue?


## Post 3 by @Laure.Moinat (2025-01-28T10:15:16.009Z)

Yesterday it wasn’t working, but today the nodes are available.
Thank you!
Best,
Laure
