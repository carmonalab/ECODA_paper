# Squeue unable to contact slurm controller

- Source: https://hpc-community.unige.ch/t/3769

- Created: 2024-12-17T13:17:18.611Z

- Posts: 2

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Cody.Cardenas (2024-12-17T13:17:18.648Z)

I just noticed this while creating a new job…
```
$ squeue
slurm_load_jobs error: Unable to contact slurm controller (connect failure)
```
I don’t have any jobs listed in my jobs folder, but there are also no `slurm-#######.out` files indicating that my most recent job failed in some way.
Is anyone else getting this behavior?
is it related to this issue[this issue](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2023-current-issues-on-hpc-cluster/2658) from last year?


## Post 2 by @Cody.Cardenas (2024-12-18T09:12:49.716Z)

If anyone lands here, I just checked and squeue seems like its working again.
