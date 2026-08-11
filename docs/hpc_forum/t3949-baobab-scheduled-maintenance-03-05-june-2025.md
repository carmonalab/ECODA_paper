# Baobab scheduled maintenance: 03-05 June 2025

- Source: https://hpc-community.unige.ch/t/3949

- Created: 2025-05-14T12:56:18.604Z

- Tags: baobab

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-05-14T12:56:18.704Z)

Dear users,
We will be performing software and hardware maintenance on the Baobab HPC cluster from 3 to 6 June 2025, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after. What will be done during this maintenance
- Reinstall all the compute nodes with latest bugfix and security patches (probably Rocky 9.6)
- Upgrade Slurm to version 24.11.5.
- Migrate our admin server to our VM cluster
- Update all servers to the latest security and bugfix releases.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Adrien.Albert (2025-06-03T07:06:48.561Z)




## Post 3 by @Adrien.Albert (2025-06-06T10:01:30.539Z)

Dear users,
The maintenance is now complete; thank you for your patience.
What has been done:
- Fixed BeeGFS tuning to improve performance
- Performed a minor Rocky Linux update
- Upgraded Slurm to version 24.11.5
- Limit resources on the login node (2 CPUs + 8 GB RAM):
- This helps keep the login node fast and responsive for everyone by preventing any single user from overloading it.
- Completed the final server migration to the virtual stack
- Reinstalled all compute nodes and updated all servers
- Applied various minor fixes
Unfortunately, we lost all Slurm queue information during the maintenance. All pending jobs have been lost.
We apologize for the inconvenience caused.


## Post 4 by @Adrien.Albert (2025-06-06T13:01:19.806Z)
