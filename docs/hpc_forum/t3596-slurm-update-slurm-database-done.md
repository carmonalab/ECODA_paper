# [Slurm] Update slurm database - done

- Source: https://hpc-community.unige.ch/t/3596

- Created: 2024-08-19T10:52:13.785Z

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-08-19T10:52:13.875Z)

Dear HPC Users,
We would like to inform you that we are currently performing an update on the SlurmDBD (Slurm Database Daemon). This update is necessary to ensure the ongoing performance,  security and carry out next maintenances of our clusters.
Important Notes:
- Start Time: The update began at 2024-08-19T11:30:00Z .
- Duration: Unknown
- Impact: Commands such as `sacct`, `sreport`, and `sacctmgr` querying the Slurm database will return errors. This is expected behavior, and normal operation will resume once the update is complete.
- No Impact on Running Jobs: Please note that this update will not affect any running jobs, and you can continue to schedule new jobs as usual. Only commands that interact with the Slurm database are affected.
We apologize for any inconvenience this may cause and appreciate your understanding.
Thank you for your cooperation.
Best regards,


## Post 2 by @Adrien.Albert (2024-08-19T13:10:06.579Z)

Dear HPC Users,
The slurmdbd has been succefully updated.
If you experience any abnormal behavior, please let us know.
