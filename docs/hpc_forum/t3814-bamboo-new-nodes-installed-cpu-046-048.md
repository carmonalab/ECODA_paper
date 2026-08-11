# [bamboo] New nodes installed: cpu[046-048]

- Source: https://hpc-community.unige.ch/t/3814

- Created: 2025-02-07T12:37:50.133Z

- Tags: bamboo

- Posts: 1

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Gael.Rossignol (2025-02-07T12:37:50.233Z)

Dear users, we have installed 3 new compute nodes on Bamboo:
```
(bamboo)-[root@admin1 minion_id]$ sinfo -N --format="%10N %.6D %20P %.4c %.8z %.6m %.8d %G %.57f" -n cpu[046-048]
NODELIST    NODES PARTITION            CPUS    S:C:T MEMORY TMP_DISK GRES                                            AVAIL_FEATURES
cpu046          1 shared-cpu            128   2:64:1 512000   200000 (null)                                             EPYC-7763,V10
cpu046          1 private-cui-cpu       128   2:64:1 512000   200000 (null)                                             EPYC-7763,V10
cpu047          1 shared-cpu            128   2:64:1 512000   200000 (null)                                             EPYC-7763,V10
cpu047          1 private-cui-cpu       128   2:64:1 512000   200000 (null)                                             EPYC-7763,V10
cpu048          1 shared-cpu            128   2:64:1 512000   200000 (null)                                             EPYC-7763,V10
cpu048          1 private-cui-cpu       128   2:64:1 512000   200000 (null)                                             EPYC-7763,V10
```
Have a nice day.
