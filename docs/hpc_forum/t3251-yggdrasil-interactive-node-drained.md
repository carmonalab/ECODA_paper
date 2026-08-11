# Yggdrasil interactive node drained

- Source: https://hpc-community.unige.ch/t/3251

- Created: 2024-01-15T10:03:46.362Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yuzheng.Kang (2024-01-15T10:03:46.417Z)

## Primary informations
Username: $kang
Cluster: $Yggdrasil
## Description
Yggdrasil Public-interactive-cpu drained. Any expectation on when this node could be recovered?
## Steps to Reproduce
sinfo


## Post 2 by @Adrien.Albert (2024-01-17T17:04:21.801Z)

Hi @Yuzheng.Kang[@Yuzheng.Kang](https://hpc-community.unige.ch/u/yuzheng.kang)
For more information about draining nodes, run ‘sinfo -R’. In this case, the node has been drained for a scheduled reboot due to a scratch issue.
We check drained nodes everydays, and if the node is ready for production, we restore it as quickly as possible. However, some issues may require more time for resolution.
It is not necessary to contact us for updates, but rest assured that we will communicate if we believe an intervention/incident requires or takes more time than expected. In the meantime, we recommend use other clusters when some resources are unavailable.
Thank you for your understanding.
