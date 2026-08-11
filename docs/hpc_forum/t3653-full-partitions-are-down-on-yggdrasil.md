# Full partitions are down on yggdrasil

- Source: https://hpc-community.unige.ch/t/3653

- Created: 2024-09-28T07:15:47.096Z

- Tags: yggdrasil

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthias.Kruckow (2024-09-28T07:15:47.134Z)

Currently, there are are partitions, where all nodes are down:
PARTITION              AVAIL  TIMELIMIT  NODES  STATE NODELIST
debug-gpu                 up      15:00      1  down* gpu001
private-euclid            up 7-00:00:00     10  down* cpu[125-134]
private-astro-cpu         up 7-00:00:00     18  down* cpu[123-124,135-150]
Is there another big problem arising like [2024] Current issues on HPC Cluster - #23 by Gael.Rossignol[[2024] Current issues on HPC Cluster - #23 by Gael.Rossignol](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/23)?
According to the monitoring of the effected nodes, they went suddenly down short before mid night (Sep.27-28), here an example of node cpu150[example of node cpu150](https://monitor.hpc.unige.ch/d/2SktsCHmf/host-overview-single?orgId=1&var-host=cpu150.yggdrasil&var-interval=$__auto_interval_interval&from=now-24h&to=now).


## Post 2 by @daniel.forerosanchez (2024-09-30T09:08:18.014Z)

Hi, Any news on this? It is quite hard to get some GPU time even for short (<15min) jobs.
Thanks!


## Post 3 by @Yann.Sagon (2024-09-30T09:57:51.967Z)

Please check the reason here: [2024] Current issues on HPC Cluster - #24 by Yann.Sagon[[2024] Current issues on HPC Cluster - #24 by Yann.Sagon](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/24)
