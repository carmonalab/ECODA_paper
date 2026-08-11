# New compute node installed cpu352.baobab

- Source: https://hpc-community.unige.ch/t/4188

- Created: 2026-01-09T15:08:25.717Z

- Tags: baobab

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Gael.Rossignol (2026-01-09T15:08:25.793Z)

Hello,
I pleased to announce we have installed a new cpu server on baobab: cpu352
```
(baobab)-[root@admin1 ~]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n cpu352
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
cpu352          1 shared-cpu            192   2:96:1 773000   873000 (null)                                             EPYC-9654,V12
cpu352          1 private-fpse-cpu      192   2:96:1 773000   873000 (null)                                             EPYC-9654,V12
```
Best regards,
