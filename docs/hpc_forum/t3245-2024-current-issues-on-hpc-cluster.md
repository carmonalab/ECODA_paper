# [2024] Current issues on HPC Cluster

- Source: https://hpc-community.unige.ch/t/3245

- Created: 2024-01-12T08:02:17.856Z

- Posts: 26

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-01-12T08:02:17.946Z)

# Yggdrasil: Scratch storage issue
Dear users,
there is an issue on scratch storage on Yggdrasil
- Our team is actively working on resolving this problem, and we apologize for any inconvenience caused. We will keep you updated on the progress.
Thank you for your understanding.
Best regards,
update 15.01.2024:
- We’ll stop the scratch storage completely until further notice.
- We’ll send today the disks to the vendor and they’ll try to rebuild the RAID in their lab.
- We’ll increase temporarily the disk quota $HOME on Yggdrasil from current 1TB to 3 TB in order to let you continue to work on the cluster.
update 16.01.2024:
- login1.yggdrasil restarted without scratch storage
- quota on home storage increased from 1TB to 3TB
- vendor received the disks, they are trying to rebuild the RAID. First attempt failed, they are first trying to clone the faulty disk.
update: 26.01.2024
- scratch storage restored.
- we removed corrupted files. If we removed files belonging to you, you’ll find in your home directory a file named `migration` with the list of erased files.
### Status : Solved :green_circle:
start: 2024-01-12T06:00:00Z
end:  2024-01-26T10:00:00Z


## Post 2 by @Yann.Sagon (2024-01-15T13:00:08.032Z)

# BAOBAB: fast storage down
Dear users,
the fast storage `/srv/fast` is currently down since the latest maintenance.
We need to go on site to put it back on production, this should be done today.
We apologize for the inconvenience.
Thank you for your understanding.
Best regards,
### Status : Solved :green_circle: :
start: 2024-01-10T02:26:00Z
end:2024-01-15T15:00:00Z


## Post 3 by @Yann.Sagon (2024-01-30T12:19:57.274Z)

# Baobab: recent software not working on legacy compute nodes and login node
Dear users,
Many software using the toolchain foss/2023b or GCCcore 13.2.0 aren’t working on legacy compute nodes or login node. The reason is they were compiled on a too recent compute node. We are rebuilding them to solve the issue.
Thanks for your understanding.
Best regards,
### Status : Solved :white_check_mark:
start: 2024-01-12T02:26:00Z
end: 2024-01-30T14:06:00Z


## Post 4 by @Yann.Sagon (2024-02-09T08:06:06.494Z)

# Baobab: scratch issue
Dear users,
one of our scratch server crashed due to too many open files. This happens when a buggy user application is opening too many files without closing them appropriately, there is not much we can do about it. If any user suspect his application may be the issue, please contact us.
Thank you for your understanding.
Best regards,
### Status : In Resolved :white_check_mark:
start: 2024-02-08T16:10:00Z
end:  2024-02-09T08:00:00Z


## Post 5 by @Adrien.Albert (2024-02-20T09:56:44.164Z)

# Baobab: scratch issue
Dear users,
One of our scratch server crashed due to too Unknown errors.  We are hardly working on it and keep you inform about the evolution.
Thank you for your understanding.
Best regards,
### Updates: 2024-02-21T12:15:00Z
We have several RaidSet in degraded state on scratch2 server. To avoid any more stress on storage, scratch have been disabled during the rebuilding of the RaidSet.
Before performing any action we have have been suspended all the jobs since 2024-02-20T14:00:00Z until 2024-02-21T11:00:00Z waiting for precise procedure from our Supplier. Then Jobs have been killed, compute and login node restarted.
A temporary scratch directory have been created in each home to hide as much as possible the impact of this issue.
### Updates: 2024-02-27T10:45:00Z
- scratch storage restored.
- No loss of data, all files created in scratch directory during the issue are present in the folder `scratch_during_scratch_failure`
### Status : In Resolved :white_check_mark:
start: 2024-02-20T09:45:00Z
end: 2024-02-27T10:45:00Z


## Post 6 by @Adrien.Albert (2024-02-21T13:19:13.394Z)




## Post 7 by @Adrien.Albert (2024-02-23T11:28:00.533Z)




## Post 8 by @Yann.Sagon (2024-03-21T09:52:18.096Z)

# Yggdrasil login node issue
Dear users,
We had a DNS issue on Yggdrasil. This prevented to use `salloc` and to connect to the allocated compute node.
### Status : Resolved :white_check_mark:
start: 2024-03-20T15:00:00Z
end: 2024-03-21T09:50:00Z


## Post 9 by @Adrien.Albert (2024-04-24T07:47:20.086Z)

# Yggdrasil login node crash
Dear users,
The connection node was blocked and remained unavailable for a few minutes.
The reason is currently unknown, but points strongly to a user process running on the connection node. We remind you that it is strictly and explicitly forbidden to run tasks on connection nodes.
If you think this is the cause of the problem, please contact us so that we can help you.
### Status : Resolved :white_check_mark:
start: 2024-04-24T07:25:00Z
end: 2024-04-24T08:37:00Z


## Post 10 by @Gael.Rossignol (2024-05-17T08:08:41.810Z)

# Baobab infiniband leaf12 crash
Dear users,
Yesterday we encounter some problems on the infiniband swich leaf12 impacting all servers connected to this switch ( gpu[030-045,047],cpu[312-321,332] ) acessing to storage home and scratch.
This leaf has been rebooted and all traffic is now working fine.
Sorry for inconvenience,
### Status : Resolved :white_check_mark:
start: 2024-05-16T13:00:00Z
end: 2024-05-16T14:00:00Z


## Post 11 by @Yann.Sagon (2024-06-03T07:08:39.047Z)

# Baobab: Slurm issue
Dear users,
Our slurm database (slurmdbd) crashed this week-end due to lack of disk space. We restored it quickly but we didn’t noticed slurm controller crashed as well.
Both services are now restored.
Thank you for your understanding.
Best regards,
---
### Status :  Resolved :white_check_mark:
Start: 2024-06-01T06:00:00Z
End:2024-06-03T07:00:00Z


## Post 12 by @Adrien.Albert (2024-06-04T14:05:11.475Z)

# Yggdrasil : Login node stuck
Dear users,
The login node (login1.yggdrasil) is blocked due to a computation process.
The server has been restarted and is available again.
Thank you for your understanding.
Best regards,
---
### Status :  Resolved :white_check_mark:
Start: 2024-06-04T13:40:00Z
End:2024-06-04T14:05:00Z


## Post 13 by @Adrien.Albert (2024-06-11T12:23:43.514Z)

# Yggdrasil : Login node stuck
Dear users,
The login node (login1.yggdrasil) is blocked due to a computation process.
The server has been restarted and is available again.
Thank you for your understanding.
Best regards,
---
### Status :  Resolved :white_check_mark:
Start: 2024-06-11T07:30:00Z
End:2024-06-11T12:15:00Z


## Post 14 by @Adrien.Albert (2024-06-17T08:18:27.421Z)

# Baobab : Storage latency
Dear users,
We are experiencing some latency on Baobab storage (home and scratch); the root cause is currently unknown. We are investigating
Thank you for your understanding.
Best regards,
---
### Status :  In progress :orange_circle:
Start: 2024-06-19T12:00:00Z
End:


## Post 15 by @Adrien.Albert (2024-06-19T12:51:33.613Z)

## All Clusters: GIO mount storage from NASAC
Dear User,
We have recently discovered issues with mounting storage spaces using the CIFS/SMB protocol from NASAC.
The problem comes from the update of Mellanox suite of packages, which manages our Infiniband (fast) network. Unfortunately, these packages disable the CIFS module, preventing CIFS/SMB storage mounts.
Since the Mellanox suite is crucial for the cluster’s operation, we cannot remove it. Currently, Mellanox doesn’t provides a workaround.
We are actively investigating a solution on our end to ensure you have the best possible experience. We appreciate your patience and will keep you updated on our progress.
Thank you for your understanding.
### Update:
2024-06-24T10:30:00Z→2024-06-24T11:45:00Z
- [HPC][Boabab] Monday 2024-06-24:12:30:00 Action on login node - #5 by Adrien.Albert[[HPC][Boabab] Monday 2024-06-24:12:30:00 Action on login node - #5 by Adrien.Albert](https://hpc-community.unige.ch/t/hpc-boabab-monday-2024-06-2430-00-action-on-login-node/3508/5)
2024-10-17T08:00:00Z
During the maintenance the patch cifs have been deployed.
```
(yggdrasil)-[alberta@login1 ~]$ gio  mount   smb://isis.unige.ch/nasac/hpc_exchange/backup < .credentials
Password required for share nasac on isis.unige.ch
User [alberta]: Domain [SAMBA]: Password: 

(yggdrasil)-[alberta@login1 ~]$ ps -edf |grep gvfs
alberta   594688       1  0 09:59 ?        00:00:00 /usr/libexec/gvfsd
alberta   594693       1  0 09:59 ?        00:00:00 /usr/libexec/gvfsd-fuse /run/user/401775/gvfs -f -o big_writes

(yggdrasil)-[alberta@login1 ~]$ ls /run/user/401775/gvfs
'smb-share:server=isis.unige.ch,share=nasac

(yggdrasil)-[alberta@login1 ~]$ ls /run/user/401775/gvfs/smb-share\:server\=isis.unige.ch\,share\=nasac/hpc_exchange/backup/
titi  toto
```
:warning:  Note that the module has been deactivated by mellanox (driver fiber provider) for a good reason, so we’re not immune to unexpected behaviour.
---
2024-12-04T23:00:00Z
Since Rocky9 was deployed on Bamboo during the last maintenance, the module now seems to be available again. You should be able to mount your NASAC share on both the login nodes and compute nodes.
Rocky9 will also be deployed on Baobab and Yggdrasil during the next maintenances.
## Status: Baobab: :orange_circle: | Yggdrasil: :orange_circle:  | Bamboo: :green_circle:
start: 2024-06-13T22:00:00Z
Patched: 2024-10-17T08:00:00Z


## Post 16 by @Gael.Rossignol (2024-06-24T11:55:01.551Z)

# Baobab : Storage scratch down
Dear users,
We are experiencing some problems on Baobab scratch storage; The root cause was the number of files open at the same time on the scratch. We had some issue to reboot server and get service up again but all is now resolved.
Thank you for your understanding.
Best regards,
---
### Status : Solved :green_circle:
Start: 2024-06-24T09:00:00Z
End: 2024-06-24T16:30:00Z


## Post 17 by @Yann.Sagon (2024-07-02T13:53:14.779Z)

# BAMBOO: Login node not reachable from outside of the university
Dear users,
it is not possible to connect to Bamboo login node from outside of the university. We’ll solve the issue.
In the meantime your options are:
- use the UNIGE VPN[UNIGE VPN](https://catalogue-si.unige.ch/en/vpn)
- connect through Baobab or Yggdrasil login node
Thank you for your understanding.
Best regards,
### Status : Solved :green_circle:
start: 2024-07-01T09:30:00Z
end: 2024-07-04T07:00:00Z


## Post 18 by @Adrien.Albert (2024-07-03T17:42:23.768Z)

# Baobab: DNS issue
Dear users,
We are currently experiencing DNS issue on Baobab, some of you have already encountered error messages as a result:
```
srun: error: _fwd_tree_get_addr: can't find address for host cpu002, check slurm.conf
srun: error: Task launch for StepId=11157586.0 failed on node cpu002: Can't find an address, check slurm.conf
srun: error: Application launch failed: Can't find an address, check slurm.conf
srun: Job step aborted
```
We are working on it to resolve it as soon as possible.
We apologize for any inconveniance caused
Thank you for your understanding.
### Status : Solved :green_circle:
start: 2024-07-01T16:30:00Z
end:2024-07-04T07:45:00Z


## Post 19 by @Adrien.Albert (2024-07-10T09:43:38.527Z)

# ALL Cluster: DPNC beegfs stuck
Dear DPNC users,
DPNC beegfs `/dpnc/beegfs` storage is currently blocked, resulting in all compute nodes, whose jobs uses this share, being drained from production. We have contacted the technical manager of this storage and are awaiting a solution.
In the meantime, we kindly ask you to wait before starting any job on this storage, and to cancel those currently using it.
Thank you for your understanding.
### Status :  Solved :green_circle:
start: 2024-07-10T07:00:00Z
end: 2024-07-10T10:03:00Z


## Post 20 by @Adrien.Albert (2024-07-10T09:44:08.242Z)
