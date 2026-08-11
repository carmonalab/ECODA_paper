# [Slurm] Scheduled Reboot for VM Hosting SLURM Database

- Source: https://hpc-community.unige.ch/t/3697

- Created: 2024-10-23T08:56:35.550Z

- Tags: slurm

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-10-23T08:56:35.644Z)

Subject: Scheduled Reboot for VM Hosting Database and SLURM Upgrade
Dear Users,
We would like to inform you that the VM hosting the database for SLURM has been successfully upgraded and requires a reboot to apply the changes.
Details of the Intervention:
- Date/Time: 2024-10-23T09:00:00Z
- Expected Downtime: 5 minutes
- Impact on Commands:
  - During the reboot, commands such as `squeue` and `sacct` may not function properly.
  - Submitting jobs will not be affected.
  - Running jobs will NOT be impacted.
The service will be fully restored immediately after the reboot. We expect minimal disruption.
Thank you for your understanding.
Best regards,
HPC Team


## Post 2 by @Adrien.Albert (2024-10-23T09:07:29.940Z)

The upgrade have been successfully completed. All services are back online.
Best Regards,
–
HPC Team
