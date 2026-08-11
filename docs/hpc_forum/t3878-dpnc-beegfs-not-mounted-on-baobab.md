# /dpnc/beegfs not mounted on Baobab

- Source: https://hpc-community.unige.ch/t/3878

- Created: 2025-03-20T14:53:52.807Z

- Posts: 11

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2025-03-20T14:53:52.844Z)

## Primary informations
Username: coppinp
Cluster: Baobab
## Description
/dpnc/beegfs is not mounted on the Baobab login node
(I didn’t check the compute nodes)
## Steps to Reproduce
(baobab)-[coppinp@login1 ~]$ ls /dpnc
(baobab)-[coppinp@login1 ~]$
## Expected Result
On Yggdrasil or Bamboo:
(yggdrasil)-[coppinp@login1 ~]$ ls /dpnc
beegfs
(yggdrasil)-[coppinp@login1 ~]$


## Post 2 by @Gael.Rossignol (2025-03-20T18:51:52.986Z)

Dear Paul,
Yes it’s possible that you had some problems to connect to this folder we had a scheduled a maintenance this afternoon :
baobab scheduled maintenance 20th of march 2025[baobab scheduled maintenance 20th of march 2025](https://hpc-community.unige.ch/t/baobab-scheduled-maintenance-20-march-2025/3868)
Sorry for inconvenance,
Best regards,


## Post 3 by @Paul.Coppin (2025-03-20T19:14:06.701Z)

Hi @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol)
Please note that this issue is still ongoing and already started yesterday evening.
```
(baobab)-[coppinp@login1 ~]$ ls /dpnc/beegfs
ls: cannot access '/dpnc/beegfs': No such file or directory
```
We knew about the maintenance and therefore did not contact you immediately, as we were thinking that the issue could perhaps be automatically resolved when Baobab came back online after the maintenance.
To give more info: /dpnc/beegfs is a BeeGFS filesystem with the BeeGFS client running on grid06. It is mounted on Baobab, Yggdrasil, and Bamboo as a NFS filesystem from grid06:/mnt/beegfs
On Yggdrasil and Bamboo we see the filesystem fine:
```
(yggdrasil)-[coppinp@login1 ~]$ ls /dpnc/beegfs
ams  atlas  cta  dampe  etc  fast  herd  neutrino  pan  share  sys  usage  users
```
but on the Baobab login node it seems that the mount is missing/broken. Again, I did not check the computing node. Can you have a look?
Cheers,
Paul


## Post 4 by @Gael.Rossignol (2025-03-20T20:04:26.502Z)

Dear Paul,
As explained in the previous post there are no issue now. Maybe you have wrong path :
```
(baobab)-[root@login1 ~]$ df -kh /srv/beegfs/dpnc/
Filesystem      Size  Used Avail Use% Mounted on
beegfs_dpnc     568T  500T   69T  88% /srv/beegfs/dpnc

(baobab)-[root@login1 ~]$ ls /srv/beegfs/dpnc/
groups
```
Best regards,


## Post 5 by @Paul.Coppin (2025-03-20T21:39:22.915Z)

Hi @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol)
Yes, /srv/beegfs/dpnc/ is working fine.
BUT I am talking about another filesytem: `/dpnc/beegfs`
`/dpnc/beegfs` is the correct path. It is a 900 TB server set up and maintained by Yann Meunier and myself from the DPNC. You can check on Yggdrasil and Bamboo that the mount is there and working fine. On Baobab the mount is not there. More details on this filesystem are provided in my previous message.
Thank you though for the late night reply!
All the best,
Paul


## Post 6 by @Gael.Rossignol (2025-03-21T06:36:06.615Z)

Dear Paul,
Sorry for confusion, I check now and server refuse connection :
```
(baobab)-[root@login1 dpnc]$ mount grid06.unige.ch:/mnt/beegfs /mnt
mount.nfs: access denied by server while mounting grid06.unige.ch:/mnt/beegfs
```
It’s working on nodes :
```
(baobab)-[root@gpu034 beegfs]$ df -kh .
Filesystem                   Size  Used Avail Use% Mounted on
grid06.unige.ch:/mnt/beegfs  916T  688T  229T  76% /dpnc/beegfs
```
Could you please check server side?
Best regards,


## Post 7 by @Paul.Coppin (2025-03-21T07:57:10.542Z)

Hi Gael,
I added
`129.194.9.190(rw,async,mp,no_subtree_check,fsid=200,crossmnt)`
to
`grid06:/etc/exports`
and ran
`exportfs -ra`
Can you try again? As to the cause, is it possible that the IP address of the baobab login node changed?
All the best,
Paul


## Post 8 by @Gael.Rossignol (2025-03-21T08:04:28.402Z)

Dear Paul,
It’s now working fine :slight_smile:
```
grid06.unige.ch:/mnt/beegfs on /mnt type nfs4 (rw,relatime,vers=4.2,rsize=1048576,wsize=1048576,namlen=255,hard,proto=tcp,timeo=600,retrans=2,sec=sys,clientaddr=129.194.9.190,local_lock=none,addr=192.33.218.206)
```
We did not change the ip of the login node since a long time.
Best regards,


## Post 9 by @Paul.Coppin (2025-03-21T08:24:29.986Z)

Dear Gael,
Ok, odd. I will check with Yann Meunier to try and understand if anything changed on our side that I was unaware of, and to prevent any issues in the future.
Thanks a lot for the help!
All the best,
Paul


## Post 10 by @Gael.Rossignol (2025-03-21T08:55:48.810Z)

Dear Paul,
I checked with Yann, and it seems we change route during previous maintenance because public interface has now 10g but we did not expect this result.
Sorry for this issue.
Best regards,


## Post 11 by @Paul.Coppin (2025-03-21T11:16:00.540Z)

Okay, now I understand better! Thanks :slight_smile:
