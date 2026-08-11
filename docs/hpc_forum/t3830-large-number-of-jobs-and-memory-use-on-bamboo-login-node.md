# Large number of jobs and memory use on bamboo login node

- Source: https://hpc-community.unige.ch/t/3830

- Created: 2025-02-20T10:16:40.939Z

- Tags: bamboo

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Raphael.Rubino (2025-02-20T10:16:40.988Z)

## Primary informations
Username: rubinor
Cluster: bamboo
## Description
Bamboo login node is slow due to high amount of RAM in use and a large number of jobs (1500+ jobs on login node from one user). Swap is also full, so I assume the jobs have peak RAM use larger than available RAM.
## Steps to Reproduce
Go to bamboo login node and check RAM use. Running htop on the login node will be quite slow already.
## Expected Result
Regular speed of login node.
## Actual Result
Slow login node.


## Post 2 by @Raphael.Rubino (2025-02-20T13:58:16.949Z)

The same is going on with yggdrasil login node now, lots of jobs running and all RAM being used.


## Post 3 by @Adrien.Albert (2025-02-20T14:26:34.796Z)

hello,
It should be resolved
