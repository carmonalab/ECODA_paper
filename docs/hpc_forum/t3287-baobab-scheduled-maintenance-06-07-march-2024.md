# Baobab scheduled maintenance: 06-07 March 2024

- Source: https://hpc-community.unige.ch/t/3287

- Created: 2024-02-05T13:23:33.460Z

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-02-05T13:23:33.555Z)

Dear users,
as just announced on the baobab-announce@ mailing list, we perform software and hardware maintenance on the Baobab HPC cluster on March 06 and 07, 2024.
The maintenance will start at 08:00 +0100 and you will receive an email when the maintenance is finished.
The cluster will be completely unavailable during this time, with no access at all (not even to retrieve files).
If you submit a job in the meantime, make sure that the expected wall time (duration) does not overlap with the start of the maintenance, or your job will be scheduled after the maintenance.
What will be done during this maintenance:
- Re-install all nodes with latest Rocky8 (8.9)
- Upgrade Slurm to version 23.11.3
- Replace batteries on compute nodes
- Upgrade servers with latest security and bug fixes
Thank you for your understanding.
Best regards,
The HPC Team
