# Baobab scheduled maintenance: 02-05 December 2025

- Source: https://hpc-community.unige.ch/t/4134

- Created: 2025-10-28T17:28:13.104Z

- Tags: baobab

- Posts: 6

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-10-28T17:28:13.222Z)

Dear users,
as just announced on the baobab-announce@ mailing list, we will be performing software and hardware maintenance on the Baobab HPC cluster from 02 to 05 December 2025, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
We remember you that you our two other HPC clusters are available during the maintenance!
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after.
What will be done during this maintenance
- All compute nodes will be reinstalled with latest version of Rocky 9
- Upgrade Slurm to version 25.05.4[25.05.4](https://www.schedmd.com/slurm-version-25-05-4-is-now-available/).
- Update all servers to the latest security and bugfix releases.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Adrien.Albert (2025-12-02T07:47:47.266Z)




## Post 3 by @Gael.Rossignol (2025-12-05T11:36:04.770Z)




## Post 4 by @Gael.Rossignol (2025-12-05T11:36:36.589Z)

Dear users,
Maintenance on the Baobab HPC cluster has now been completed. However, an issue with OpenOnDemand remains, and we are actively working to resolve it as soon as possible.
Best regards,


## Post 5 by @Paul.Coppin (2025-12-05T12:20:25.286Z)

Dear Gael,
Thank you and the HPC team for completing this maintenance and update!
We are currently still unable to access the login node via ssh (ask for password and rejects the correct one). With the maintenance over, can we already access the login node again for access to scratch storage and Slurm job submission, or will this have to wait until the OpenOnDemand problem is solved?
All the best,
Paul


## Post 6 by @Ludovic.Dumoulin (2025-12-05T12:26:39.418Z)

Hello,
Same issue here, using Password or ssh key.
Best,
Ludovic


## Post 7 by @Gael.Rossignol (2025-12-05T13:04:36.396Z)

Hello,
Issue is now resolved, we forgot one step.
Sorry for inconvenience,
Best regards,


## Post 8 by @Adrien.Albert (2025-12-05T13:53:28.434Z)

Baobab OpenOndemand[Baobab OpenOndemand](https://openondemand.baobab.hpc.unige.ch/) is back again !
