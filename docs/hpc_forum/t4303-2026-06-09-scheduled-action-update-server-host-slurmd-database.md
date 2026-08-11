# [2026-06-09] scheduled Action: Update server host slurmd database

- Source: https://hpc-community.unige.ch/t/4303

- Created: 2026-06-08T15:05:17.068Z

- Tags: all

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2026-06-08T15:05:17.204Z)

Dear  HPC users,
Please be advised that the server hosting the Slurm database (slurmdbd) is scheduled for maintenance on:
:date: Tuesday, 09 June 2026
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


## Post 2 by @Adrien.Albert (2026-06-09T10:29:54.391Z)

Dear HPC Users,
2026-06-09T09:47:00Z:
The maintenance has been successfully completed. The Slurm database service is now fully operational, and all Slurm commands should respond as expected.
Despite our validation tests, if you encounter any issues, please do not hesitate to contact us.
Best regards,
HPC Team
