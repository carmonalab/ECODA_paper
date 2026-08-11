# [baobab] No space left on device - /home

- Source: https://hpc-community.unige.ch/t/4289

- Created: 2026-05-01T06:42:19.606Z

- Tags: baobab

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2026-05-01T06:42:19.688Z)

Dear HPC Team,
Baobab’s /home is full; the cluster is not usable in this state.
```
beegfs_home                                         138T  138T   88G 100% /home
```
Kind regards,
Maciej Falkiewicz


## Post 2 by @Yann.Sagon (2026-05-04T08:28:40.939Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz) thanks for letting us know. Please see this announce for the action we’ll do as workaround. Home storage quota reduction on the Baobab cluster[Home storage quota reduction on the Baobab cluster](https://hpc-community.unige.ch/t/home-storage-quota-reduction-on-the-baobab-cluster/4290)


## Post 3 by @maciej.falkiewicz (2026-05-04T08:40:36.419Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) ,
Thank you for the reply.
While I expect the new policy to solve the problem for some time, in X amount of time, we will face the same problem again. At the same time, the current quota (and the new one even more) is far from satisfactory. There are training datasets that are bigger than that. Offloading to `scratch` is not an option either due to how slow it is.
I am also surprised that you didn’t take any action earlier to prevent the current situation.
Kind regards,
Maciej
