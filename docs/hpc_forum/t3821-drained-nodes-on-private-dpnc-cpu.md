# Drained nodes on `private-dpnc-cpu`

- Source: https://hpc-community.unige.ch/t/3821

- Created: 2025-02-12T16:31:25.382Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Andrea.Serpolla (2025-02-12T16:31:25.431Z)

## Primary informations
Username: serpolla
Cluster: Baobab
## Description
Dear HPC team,
I noticed many drained nodes for the `private-dpnc-cpu` queue (12 out of 17 currently).
I don’t know if there is any issue ongoing.
Below `sinfo` output:
```
(baobab)-[serpolla@login1 ~]$ sinfo -p private-dpnc-cpu
PARTITION        AVAIL  TIMELIMIT  NODES  STATE NODELIST
private-dpnc-cpu    up 7-00:00:00      2   drng cpu[226,277]
private-dpnc-cpu    up 7-00:00:00     10  drain cpu[084-088,210-211,213,227,229]
private-dpnc-cpu    up 7-00:00:00      1    mix cpu212
private-dpnc-cpu    up 7-00:00:00      4  alloc cpu[089-090,209,228]

(baobab)-[serpolla@login1 ~]$ sinfo -R -p private-dpnc-cpu
REASON               USER      TIMESTAMP           NODELIST
Kill task failed     root      2025-02-11T16:35:57 cpu[084,213]
health_BEEGFS__tcp_c root      2025-02-12T12:06:14 cpu[085,088,210-211,226]
Kill task failed     root      2025-02-11T16:30:56 cpu086
health_BEEGFS__tcp_c root      2025-02-12T12:06:15 cpu087
health_BEEGFS__tcp_c root      2025-02-12T11:51:13 cpu227
Kill task failed     root      2025-02-11T16:30:57 cpu[229,277]
```
Best,
Andrea


## Post 2 by @Gael.Rossignol (2025-03-06T17:18:50.622Z)

Dear Andrea,
No issues were at this special time, but we were facing some lack of storage space on cluster. Since your post storage has been released by users (thank you very much for your help!) and nodes running again in production.
Best regards,
