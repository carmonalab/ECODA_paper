# Max number of jobs

- Source: https://hpc-community.unige.ch/t/3308

- Created: 2024-02-14T09:08:37.261Z

- Posts: 3

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jingze.Duan (2024-02-14T09:08:37.327Z)

Good morning,
My jobs on yggdrasil are pending for a few day. And I noticed that sometimes the same user are running about 30 jobs and using almost all the available shared gpu from the squeue list.
I also found this
https://hpc-community.unige.ch/t/max-number-of-running-and-pending-jobs-per-user/239[https://hpc-community.unige.ch/t/max-number-of-running-and-pending-jobs-per-user/239](https://hpc-community.unige.ch/t/max-number-of-running-and-pending-jobs-per-user/239)
It said the limitation of running and pending jobs per user is 10k. Could you mind my asking that is this limitation reasonable compared to the total number of the available shared gpu?
It’s just an idea. But how about setting a limit on the maximum percentage of the resource proportion? Like one user can only use no more than 60% of the total computational cpu/gpu. So that the computational resources can be cycled more healthily.
Thank you for your time. I’m looking for your insights.
Best,
Jingze


## Post 2 by @Matthias.Kruckow (2024-02-14T13:13:09.823Z)

The maximum number of jobs is the same for using cpus and gpus. For the cpus it makes more sense to have 10k as limit there, when thinking of jobs running on a single cpu.
Your request to better distribute the computational resources is already taken care of by the Job priorities[Job priorities](https://doc.eresearch.unige.ch/hpc/slurm#job_priorities) especially the `fairshare`. (If you see, that the fairshare doesn’t work correctly for the gpus, may give a more detailed report.)


## Post 3 by @Jingze.Duan (2024-02-14T14:26:40.821Z)

Thank you for your reply!
