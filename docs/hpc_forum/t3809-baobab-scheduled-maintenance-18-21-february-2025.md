# Baobab scheduled maintenance: 18-21 February 2025

- Source: https://hpc-community.unige.ch/t/3809

- Created: 2025-02-03T15:28:40.070Z

- Tags: baobab

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-02-03T15:28:40.192Z)

Dear users,
as just announced on the baobab-announce@ mailing list, we will be performing software and hardware maintenance on the Baobab HPC cluster from 18 to 21 February 2025, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after. What will be done during this maintenance
- All compute nodes will be reinstalled with a new major OS release, Rocky 9. Compute nodes are currently on Rocky 8.
- Upgrade Slurm to version 24.11.0.
- use NFSv4 for `/srv/fast`
- Migrate our physical admin servers to our new virtual machine cluster
- Update all servers to the latest security and bugfix releases.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Gael.Rossignol (2025-02-18T06:53:35.808Z)




## Post 3 by @Adrien.Albert (2025-02-24T16:29:51.414Z)
