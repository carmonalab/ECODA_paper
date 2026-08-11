# Nodes in drng state

- Source: https://hpc-community.unige.ch/t/3629

- Created: 2024-09-09T19:23:10.126Z

- Tags: bamboo

- Posts: 1

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Patrick.Remy (2024-09-09T19:23:10.181Z)

Hello,
there are quite a few nodes on bamboo in drng state:
```
(bamboo)-[remypa@login1 simulations]$ sinfo -R -l
Mon Sep 09 21:06:01 2024
REASON               USER         TIMESTAMP           STATE  NODELIST
Not responding       slurm(200)   2024-09-09T09:41:55 down*  cpu[014,038]
health_BEEGFS__tcp_c root(0)      2024-09-09T09:47:58 drng   cpu[018-021,023-024,026-027,029-030,033-034,036-037,039,041-043]
health_BEEGFS__tcp_c root(0)      2024-09-09T19:51:01 drng   cpu[017,028,035,040]
health_BEEGFS__tcp_c root(0)      2024-09-09T09:35:58 drain  cpu025
```
Thank you very much
Patrick
