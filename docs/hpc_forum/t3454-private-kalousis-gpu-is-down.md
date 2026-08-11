# `private-kalousis-gpu` is down

- Source: https://hpc-community.unige.ch/t/3454

- Created: 2024-05-21T15:21:21.512Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Imahn.Shekhzadeh1 (2024-05-21T15:21:21.558Z)

```
(base) (baobab)-[shekhza2@cpu322 benchmark]$ squeue -p private-kalousis-gpu
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          10159856 private-k ant-migr shekhza2 PD       0:00      1 (Nodes required for job are DOWN, DRAINED or reserved for jobs in higher priority partitions)
```


## Post 2 by @Adrien.Albert (2024-05-22T08:38:12.886Z)

Hi @Imahn.Shekhzadeh1[@Imahn.Shekhzadeh1](https://hpc-community.unige.ch/u/imahn.shekhzadeh1)
There’s no need to send or post a message about down or drained nodes. We check every day which nodes are ready to go back into production. We wait until all work has been completed before intervening on a node to avoid any impact on production :wink:
PS: In your sbatch, you can specify multiple partitions, allowing your jobs to start the first available node according to your priority.
