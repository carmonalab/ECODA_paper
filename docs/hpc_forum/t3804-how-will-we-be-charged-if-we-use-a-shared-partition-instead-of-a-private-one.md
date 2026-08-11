# How will we be charged if we use a shared partition instead of a private one?

- Source: https://hpc-community.unige.ch/t/3804

- Created: 2025-01-30T02:21:21.106Z

- Tags: all

- Posts: 2

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2025-01-30T02:21:21.106Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) , I have a quick question about your suggestion on using shared-cpu. How would billing be handled in that case? Because I would like to avoid getting billed (which is the case with the private partition)


## Post 2 by @Yann.Sagon (2025-01-30T09:15:41.078Z)

Dear @Bharathkumar.Radhakrishnan[@Bharathkumar.Radhakrishnan](https://hpc-community.unige.ch/u/bharathkumar.radhakrishnan), I’ve created a new thread with your topic.
If your group has access to a private partition, they have a usage right for a number of hours per year. This number depends on the compute node (number of cores, GPUs etc).
You are then free to use the private partition or the shared partition for this number of hours. The only difference if you use the private partition is that you have a higher priority.
The main advantage is that you are not restricted to using your private nodes, but can access the three clusters and even the GPUs.
We are developing scripts to allow to check the usage and the amount of hours you have the right to use regarding the hardware your group owns.
More information usage limits[usage limits](https://doc.eresearch.unige.ch/hpc/hpc_clusters#usage_limits)
Best
