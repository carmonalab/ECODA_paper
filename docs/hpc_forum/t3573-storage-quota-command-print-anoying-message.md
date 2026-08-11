# [storage] Quota command print anoying message

- Source: https://hpc-community.unige.ch/t/3573

- Created: 2024-07-26T15:16:43.896Z

- Tags: yggdrasil

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-07-26T15:16:43.990Z)

Dear HPC Users
You may have already notice this bug.
## Primary information:
Cluster: Yggdrasil (since the last maintenance in July)
## Describe the bug
Since the latest Beegfs version 7.4.4. The get/set quota print message announces BeeGFS Enterprise functionality even if sysNoEnterpriseFeatureMsg is set to true.
## To Reproduce
```
(cluster)-[root@login1 ~]$ beegfs-ctl --getquota --uid gandalf--mount=/srv/beegfs/scratch --connAuthFile=/etc/beegfs/connauthfile
--------------------------------------------------------------------------------
| BeeGFS Enterprise Feature                                                    |
|                                                                              |
| This beegfs-ctl mode configures a BeeGFS Enterprise Feature.                 |
|                                                                              |
| By downloading and/or installing BeeGFS, you have agreed to the EULA of      |
| BeeGFS: https://www.beegfs.io/docs/BeeGFS_EULA.txt                           |
|                                                                              |
| Please note that any use of Enterprise Features of BeeGFS for longer than    |
| the trial period of 60 (sixty) days requires a valid License & Support       |
| Agreement with the licensor of BeeGFS "ThinkParQ GmbH".                      |
|                                                                              |
| Contact: sales@thinkparq.com                                                 |
| Thank you for supporting BeeGFS development!                                 |
|                                                                              |
| If you are using BeeGFS in conformity with the EULA and do not wish to see   |
| this message in the future, you can set sysNoEnterpriseFeatureMsg to true in |
| beegfs-client.conf to disable it.                                            |
--------------------------------------------------------------------------------

Quota information for storage pool Default (ID: 1):

      user/group     ||           size          ||    chunk files    
     name     |  id  ||    used    |    hard    ||  used   |  hard   
--------------|------||------------|------------||---------|---------
       gandalf|76384||   80.00 KiB|   unlimited||       17| 10000000
```
## Expected behavior
No message printed
```
(cluster)-[root@login1 ~]$ beegfs-ctl --getquota --uid gandalf--mount=/srv/beegfs/scratch --connAuthFile=/etc/beegfs/connauthfile
Quota information for storage pool Default (ID: 1):

      user/group     ||           size          ||    chunk files    
     name     |  id  ||    used    |    hard    ||  used   |  hard   
--------------|------||------------|------------||---------|---------
       gandalf|76384||   80.00 KiB|   unlimited||       17| 10000000
```
## Resolution
I’ve posted on beegfs github: here[here](https://github.com/ThinkParQ/beegfs/issues/7) and we’re expecting a fix for the next release.
