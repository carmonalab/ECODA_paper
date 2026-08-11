# Baobab scheduled maintenance: 13-14 June 2024

- Source: https://hpc-community.unige.ch/t/3450

- Created: 2024-05-17T08:59:22.718Z

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-05-17T08:59:22.822Z)

Dear users,
As just announced on the baobab-announce@ mailing list, we will perform software and hardware maintenance on the Baobab HPC cluster on June 13th and 14th.
The maintenance will start at 08:00 +0100 and you will receive an email when the maintenance is finished.
The cluster will be completely unavailable during this time and you will not be able to access it at all (not even to retrieve files).
If you submit a job in the meantime, make sure that the expected wall time (duration) does not overlap with the start of the maintenance, or your job will be scheduled after the maintenance.
What will be done during this maintenance
- A brand new 128 core login node will be brought into production: login1.baobab.hpc.unige.ch.
- Upgrade OpenOnDemand to version 3.1.0.
- Replace the 13 old Infiniband switches with a new EDR model: this will increase the bandwidth between nodes and storage to 100G!
- Migrating our slurmdbd server (used for job accounting, job querying with sacct, etc) to a dedicated VM. This may cause some disruption to slurm on Yggdrasil as well, but running jobs shouldn’t be affected.
- Upgrade Slurm to version 23.11.6-1.
- Update all servers to the latest security and bugfix releases
- Re-install all nodes with latest security and bugfixes
- Various fixes (replace a battery, correct scripts etc)
As you may have noticed, we’ll be quite busy during this maintenance. We won’t be providing user support this week unless there is an emergency.
Thanks for your understanding.
Best regards,
The HPC Team


## Post 2 by @Yann.Sagon (2024-07-02T11:29:59.820Z)

Dear users, the maintenance is over.
Important information: the new url to reach the login node is login1.baobab.hpc.unige.ch
Take some time to update your preferences in your ssh/filezilla etc clients. The old login2.baobab.hpc.unige.ch url will still work for a month. The obsolete url baobab2.hpc.unige.ch is now deleted.
We were able to make the changes listed below with two modifications:
- Slurm version is now: 23.11.7
- We haven’t replaced the Infiniband switches as we haven’t received them yet.
Enjoy the new and powerful login node, test the latest version of OpenOnDemand, take advantage of the freshly installed compute node and, as always, do not hesitate to contact us if you find anything wrong.
We apologise for having announced the wrong Maintenance Day.
HPC Team
