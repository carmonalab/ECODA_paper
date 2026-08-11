# [bamboo] gpu008 inter-GPU communication not working

- Source: https://hpc-community.unige.ch/t/4169

- Created: 2025-12-15T10:24:14.893Z

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ramon.CalvoGonzalez (2025-12-15T10:24:14.938Z)

## Primary informations
Username: calvogon
Cluster: bamboo
## Description
I am trying to run a job that requires inter-GPU communications. I know that my code works, since it runs on any other node of bamboo just fine. But when running on gpu008, the job is frozen and does not run. I think this has to do with the fact that gpu008 might have faulty inter-GPU communications, since running on a single GPU works fine.
## Steps to Reproduce
I’m trying to run jobs that require all 4 GPUs at the same time.
