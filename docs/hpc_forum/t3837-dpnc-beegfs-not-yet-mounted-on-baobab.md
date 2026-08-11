# /dpnc/beegfs not yet mounted on Baobab

- Source: https://hpc-community.unige.ch/t/3837

- Created: 2025-02-21T16:44:30.549Z

- Tags: baobab

- Posts: 5

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2025-02-21T16:44:30.597Z)

Dear HPC,
It appears that /dpnc/beegfs is not yet mounted on the login-node following the OS migration of Baobab. Would it be possible to have a look at this?
Note that I have not checked the compute nodes, and that potentially it needs to be mounted for those as well.
All the best,
Paul


## Post 2 by @Adrien.Albert (2025-02-22T00:04:46.637Z)

Hello @Paul.Coppin[@Paul.Coppin](https://hpc-community.unige.ch/u/paul.coppin)
I checked and had no problems:
```
(baobab)-[alberta@login1 ~]$ df -t beegfs
Filesystem         1K-blocks          Used    Available Use% Mounted on
beegfs_dpnc     539049799680  535486640640   3563159040 100% /srv/beegfs/dpnc
beegfs_scratch 1617149399040 1320638278144 296511120896  82% /srv/beegfs/scratch
beegfs_home     147639165952  144832774144   2806391808  99% /home
```
Are you sure the issue is on login node ?


## Post 3 by @Paul.Coppin (2025-02-22T11:32:00.051Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
Thank you for the quick reply.
I mean the directory: /dpnc/beegfs
not /srv/beegfs
/dpnc/beegfs should point to a ~900TB BeeGFS file system managed by the DPNC department on Baobab, Yggdrasil, and Bamboo. It appears it has not yet been mounted on Baobab following the OS upgrade. If my memory serves correctly, I believe it is usually mounted via NSF through grid06.unige.ch (which acts as the BeeGFS client).
All the best,
Paul


## Post 4 by @Adrien.Albert (2025-02-22T20:18:41.144Z)

Hello,
Should be OK:
```
(baobab)-[root@login1 ~]$ ls /dpnc/beegfs/
ams  atlas  cta  dampe  etc  fast  herd  neutrino  pan  share  sys  usage  users
```


## Post 5 by @Paul.Coppin (2025-02-23T11:17:04.511Z)

Indeed, looks all good now. Thanks!
