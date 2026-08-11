# [Slurm] Create public-interactive-gpu partition

- Source: https://hpc-community.unige.ch/t/3879

- Created: 2025-03-21T12:26:58.861Z

- Tags: bamboo, yggdrasil

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2025-03-21T12:26:58.946Z)

Dear Users,
Cluster: Bamboo, Yggdrasil
We would like to inform you about a recent change to the GPU partition on the system. The `gpu001` partition has been modified from `debug-gpu` to `public-interactive-gpu`. This change introduces the following constraints for all users:
### Key Changes:
- Partition Name: `gpu001` has been moved to the public-interactive-gpu partition.
- Constraints:
  - 1 Job per User: Each user is limited to running a single job at a time.
  - 1 GPU per Job: Each job will be allocated only one GPU.
  - XX CPU Max per Job: A maximum of 16 CPUs (Bamboo) and 2 CPUs (Yggdrasil) can be allocated to each job.
  - Max Walltime: The maximum allowed walltime for jobs in this partition is 04:00:00.
These changes are intended to help manage resources efficiently and ensure fair access to the GPU nodes.
Thank you for your understanding.
Best regards,
