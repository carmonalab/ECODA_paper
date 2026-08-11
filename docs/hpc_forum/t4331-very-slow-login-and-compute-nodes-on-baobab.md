# Very slow login and compute nodes on Baobab

- Source: https://hpc-community.unige.ch/t/4331

- Created: 2026-06-29T12:24:47.381Z

- Tags: baobab

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Alexander.Froch (2026-06-29T12:24:47.446Z)

Hi,
I have some issues when accessing the login and compute node of Baobab. It seems that all of them are super slow, but the login node especially. Looking at htop on the login node, it seems that the swap is completely full again. Not sure if this could cause the issue ..
Cheers,
Alex


## Post 2 by @Gael.Rossignol (2026-06-29T12:42:41.660Z)

Dear Alexander,
This morning, an issue occurred on the scratch storage. It caused a lot of latency across the servers, but everything has now been resolved.
You can find more details here: [2026] Current issues on HPC Cluster - #23 by Yann.Sagon[[2026] Current issues on HPC Cluster - #23 by Yann.Sagon](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/23)
Best regards,
