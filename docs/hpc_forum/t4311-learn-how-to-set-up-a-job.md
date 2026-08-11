# Learn how to set up a job

- Source: https://hpc-community.unige.ch/t/4311

- Created: 2026-06-12T07:59:42.734Z

- Tags: baobab

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yuv.Agarwal (2026-06-12T07:59:42.798Z)

Hi
I am a first time user of the Cluster.
I have to set up a job to run multiple variations of the same of the same CUDA code. (Around 300,000 variations of the same code - 1 of them takes about 15 mins to run.) I have arranged each variation in its own independent folder but as a first time user I don’t know what is the right way to assign this job so that it runs efficiently.
Looking for advice on how to do it effectively without hampering any queued task or jobs.
Best
Yuv


## Post 2 by @Yann.Sagon (2026-06-29T08:38:05.158Z)

Hi, you should have a look at job arrays maybe? hpc:slurm [eResearch Doc][hpc:slurm [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/slurm#job_array)
