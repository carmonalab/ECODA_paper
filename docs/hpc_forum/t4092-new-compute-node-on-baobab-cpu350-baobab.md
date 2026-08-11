# New compute node on Baobab cpu350.baobab

- Source: https://hpc-community.unige.ch/t/4092

- Created: 2025-09-18T11:39:42.884Z

- Tags: baobab

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-09-18T11:39:43.002Z)

Dear users, please welcome a new family member: cpu350.baobab
```
(baobab)-[root@admin1 ~]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n cpu350
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
cpu350          1 shared-cpu            192   2:96:1 773000   900000 (null)                                             EPYC-9654,V12
cpu350          1 private-dpnc-cpu      192   2:96:1 773000   900000 (null)                                             EPYC-9654,V12
```
