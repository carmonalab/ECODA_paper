# [2026-01-26] scheduled Action: Update server host slurmd database

- Source: https://hpc-community.unige.ch/t/4202

- Created: 2026-01-22T14:47:33.435Z

- Tags: all

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2026-01-22T14:47:33.561Z)

Dear  HPC users,
Please be advised that the server hosting the Slurm database (slurmdbd) is scheduled for maintenance on:
:date: Monday, 26 January 2026
:ten_o_clock: Starting at 10:00
During this maintenance window, the Slurm database will be temporarily unresponsive. As a result, any Slurm command that queries the database will fail. This includes, but is not limited to:
- `squeue`
- `sacct`
- `sreport`
- Any other Slurm command requiring access to the Slurm database
Job submission and scheduling may continue to function, but all accounting‑related queries will be unavailable until the maintenance is complete.
We will notify you once the service is fully restored.
Thank you for your understanding.
If you have any questions, feel free to reach out.


## Post 2 by @Adrien.Albert (2026-01-27T12:04:07.929Z)

Dear HPC Users,
2026-01-26T13:30:00Z:
The server has been successfully updated. All services are back in production, and no issues have been reported.
Thank you for your patience.
