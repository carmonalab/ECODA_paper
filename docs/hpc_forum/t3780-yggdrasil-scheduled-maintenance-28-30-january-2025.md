# Yggdrasil scheduled maintenance: 28-30 January 2025

- Source: https://hpc-community.unige.ch/t/3780

- Created: 2025-01-09T10:35:46.850Z

- Tags: yggdrasil

- Posts: 4

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-01-09T10:35:46.932Z)

Dear users,
We will be performing software and hardware maintenance on the Yggdrasil HPC cluster from 28 to 30 January 2025, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after. What will be done during this maintenance
- All compute nodes will be reinstalled with a new major OS release, Rocky 9. Compute nodes are currently on Rocky 8.
- Upgrade Slurm to version 24.11.0.
- Update all servers to the latest security and bugfix releases.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Matthias.Kruckow (2025-01-09T11:13:45.009Z)

Could you please resolve the conflicting information which cluster will get updated: the title states `Yggdrasil`, while the text states `Bamboo`. (The same holds for the email.)


## Post 3 by @Yann.Sagon (2025-01-09T12:33:29.907Z)

Indeed, thanks for the notification. The cluster that will get updated is Yggdrasil!


## Post 4 by @Gael.Rossignol (2025-01-28T07:44:50.097Z)




## Post 5 by @Yann.Sagon (2025-01-30T10:30:56.039Z)




## Post 6 by @Yann.Sagon (2025-01-30T10:36:38.104Z)

Dear users, the maintenance is now finished.
What has been done:
- reinstall all compute nodes and login node to Rocky 9.5
- upgrade slurm to version 24.11.1
- all storage servers reinstalled on Rocky 9.5
- reinstall the slurm server to Rocky 9.5
- update the admin server with the latest security and bug fixes
As always, please let us know if you find anything that does not work as expected.
Best regards,
HPC Team, Yann, Adrien, Gaël
