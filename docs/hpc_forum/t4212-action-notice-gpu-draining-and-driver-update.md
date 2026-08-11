# [Action Notice] GPU Draining and Driver Update

- Source: https://hpc-community.unige.ch/t/4212

- Created: 2026-02-03T10:33:28.906Z

- Tags: all

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2026-02-03T10:33:29.003Z)

Dear users,
We would like to inform you that we are starting a rolling update of the GPU driver across the clusters.
The update has been successfully validated on gpu009.bamboo, and we expect it to resolve several user‑reported issues.
GPU nodes will be drained progressively, and each node will be updated and rebooted as soon as it becomes free.
We do not wait for all nodes to be drained before starting: the update is applied individually, as soon as each node becomes available, in order to minimize the impact on production.
Each node will be returned to service immediately after its update.
A follow‑up announcement will be made once the full maintenance is completed.
Thank you for yout understanding
