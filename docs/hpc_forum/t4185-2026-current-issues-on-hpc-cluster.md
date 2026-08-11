# [2026] Current issues on HPC Cluster

- Source: https://hpc-community.unige.ch/t/4185

- Created: 2026-01-05T07:51:13.197Z

- Posts: 25

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-01-05T07:51:13.284Z)

## Description
Cluster: Baobab
Login node crashed with a kernel panic. We restarted it.
–
HPC Team
### Status : Resolved :green_circle:
start: 2026-01-02T18:53:00Z
end: 2026-01-05T07:50:00Z


## Post 2 by @Yann.Sagon (2026-01-05T08:31:03.828Z)

## Description
Cluster: Bamboo
Login node crashed with a kernel panic. We restarted it.
–
HPC Team
### Status : Resolved :green_circle:
start: 2025-12-31T07:17:00Z
end: 2026-01-05T07:55:00Z


## Post 3 by @Yann.Sagon (2026-01-05T08:33:45.621Z)

## Description
Cluster: Bamboo
Storage server was down for part of scratch storage.
–
HPC Team
### Status : solved :green_circle:
start: 2025-12-30T20:22:00Z
end: 2026-01-05T09:19:00Z


## Post 4 by @Adrien.Albert (2026-01-12T10:15:17.184Z)

## Description
Cluster: Bamboo
Login1 has crashed, server has been rebooted
### tatus : solved :green_circle:
start: 2026-01-12T10:06:00Z
end: 2026-01-11T10:15:00Z


## Post 5 by @Adrien.Albert (2026-01-15T09:44:06.856Z)

## Description
Cluster: Yggdrasil
We are experiencing issue with SLURM causing delays in job management. We are actively working to resolve this incident and limit its impact.
### Updates
The issue has been resolved and Slurm is back online.
### Status : Solved :green_circle:
start: 2026-01-15T10:05:00Z
end: 2026-01-15T10:30:00Z


## Post 7 by @Gael.Rossignol (2026-01-19T10:21:17.534Z)

## Description
Cluster: Yggdrasil
We are currently facing a power outage affecting the Yggdrasil cluster, which may result in multiple nodes being unreachable. We will provide an update once power has been fully restored.
### Status : Solved :green_circle:
start: 2026-01-19T10:00:00Z
end: 2026-01-19T10:30:00Z


## Post 8 by @Yann.Sagon (2026-01-19T13:22:14.275Z)

## Description
Cluster: all
Since the latest slurm update, some interactive jobs (using srun or salloc) are killed prematurely. For the user it appears as if the job reached its timelimit, but in the admin logs it is indicated that the job was killed due to inactivity timeout. We opened a case at schedMD.
### Status : Solved :green_circle:
start: 2025-12-19T10:00:00Z


## Post 9 by @Yann.Sagon (2026-01-29T11:01:37.350Z)

## Description
Cluster: baobab
We had to reboot login1.baobab because systemctl/dbus crashed.
### Status : Solved :green_circle:
start: 2026-01-29T10:45:00Z
end:2026-01-29T11:00:00Z


## Post 10 by @Gael.Rossignol (2026-02-12T08:26:04.632Z)

## Description
Cluster: yggdrasil
Due to a cut in the optical fiber, the IT and telephone networks are unavailable at the Ecogia and Sauverny sites. The OCSIN technicians who manage our optical fibers have been informed. No resolution time has been communicated for now.
## Edit:
Incident is now solved. Please not that cluster was not reachable but job were running.
### Status : Resolved :green_circle:
start: 2026-02-11T21:38:00Z
end: 2026-02-12T11:45:00Z


## Post 11 by @Adrien.Albert (2026-02-23T11:00:21.029Z)

Dear Users,
Cluster: Yggdrasil
### Description
An issue is currently affecting the Home storage on the Yggdrasil cluster.
Some users may encounter the following message:
“I/O remote error”
- Possible disruptions accessing the Home directory
- Risk of input/output errors during commands or job execution
- The cluster remains available, but some operations may fail
We are working to restore normal service as soon as possible.
An update will be posted as soon as we have more information.
We apologize for the delayed response.
Our team is operating with reduced staff this week and is focused on preparing the Bamboo maintenance starting tomorrow.
## Status: Resolved :green_circle:
Start: 2026-02-23T08:00:00Z
End: 2026-02-23T15:00:00Z


## Post 12 by @Adrien.Albert (2026-02-27T09:39:01.274Z)

Dear users,
Cluster: Bamboo
### Description
Following the recent maintenance on Bamboo, we upgraded the BeeGFS client from version 7.4.6 to 7.4.7.
After this update, we observed a change in behavior of the `beegfs-ctl --getquota` command.
When running:
```
(bamboo)-[alberta@login1 ~]$ beegfs-get-quota-home-scratch.sh -u alberta
home dir: /home/users/a/alberta
scratch dir: /srv/beegfs/scratch/users/a/alberta

          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
Unable to resolve given name: --connTCPRecvBufSize
home        | 
Unable to resolve given name: --connTCPRecvBufSize
scratch  
```
the command unexpectedly returns the following error:
```
Unable to resolve given name: --connTCPRecvBufSize
```
An issue has been opened with the BeeGFS team, and we are working to resolve this case as soon as possible.
## Status: WorkArround applied :blue_circle:
Start: 2026-02-26T15:00:00Z
End: 2026-02-27T10:00:00Z


## Post 13 by @Adrien.Albert (2026-03-02T12:59:52.889Z)

Dear users,
Cluster: Bamboo
### Description
Dear HPC users,
Following the last Bamboo maintenance[Bamboo maintenance](https://hpc-community.unige.ch/t/bamboo-scheduled-maintenance-24-27-february-2026/4216/6), the interactive RStudio application on Open OnDemand was not functioning correctly, causing sessions to terminate immediately after launch.
We have identified the cause of this behavior and applied a fix. The RStudio app should now operate as expected.
We also strongly suspect that this incident is related to a new behavior introduced by the recent upgrade of BeeGFS to version 7.4.7.
If you continue to experience any issues, please let us know.
## Status: Resolved :green_circle:
Start: 2026-02-26T15:00:00Z
End: 2026-03-02T10:00:00Z


## Post 14 by @Gael.Rossignol (2026-03-19T14:09:42.381Z)

Cluster: ALL
### Description
Dear users,
This message in only to inform that we need to update slurmdbd to latest version, this will no have any impact on the usage of clusters but some commands to monitor like sacct will be unavailable during few minutes.
Sorry for inconvenience,
## Status: Resolved :green_circle:
Start: 2026-03-19T14:00:00Z
End: 2026-03-19T16:00:00Z


## Post 15 by @Yann.Sagon (2026-03-27T14:41:08.639Z)

Cluster: Baobab
### Description
We we must perform an immediate emergency shutdown of the Baobab HPC cluster due to a critical cooling issue in the datacenter.
This action is required to protect the hardware and ensure the safety and integrity of your data.
All currently running jobs will unfortunately be stopped.
Our teams are actively working with the datacenter staff to resolve the situation as quickly as possible.
The Baobab HPC Team
## Status: Resolved :green_circle:
Start: 2026-03-27T14:00:00Z
End: 2026-03-27T17:00:00Z


## Post 16 by @Adrien.Albert (2026-04-05T10:18:59.989Z)

Cluster: Bamboo
## Description
An incident is currently ongoing on the Bamboo cluster. One of our storage server is experiencing an issue impacting the scratch and home storage.
## Status: Solved :green_circle:
Start: 2026-04-02T20:24:00Z
End: 2026-04-07T06:00:00Z


## Post 17 by @Gael.Rossignol (2026-04-17T09:18:31.814Z)

Cluster: Baobab
## Description
Hello,
We are currently experiencing an issue with login1.baobab. The server is not accessible at the moment and is currently rebooting.
Access should be restored within the next few minutes.
Best regards,
## Status: Solved :green_circle:
Start: 2026-04-17T09:00:00Z
End: 2026-04-17T09:30:00Z


## Post 18 by @Yann.Sagon (2026-04-30T15:14:45.350Z)

## Description
We have received notification of a critical security exploit. We are currently patching all our infrastructure immediately, including compute and login nodes.
We will reboot our three login nodes immediately and put the cluster into drain mode.
We will resume the cluster as soon as the nodes have been patched.
Best regards,
## Status: Resolved  :green_circle:
Start: 2026-04-30T15:00:00Z
End: 2026-05-12T22:00:00Z
Immediate maintenance – Emergency security patching[Immediate maintenance – Emergency security patching](https://hpc-community.unige.ch/t/immediate-maintenance-emergency-security-patching/4288)


## Post 19 by @Adrien.Albert (2026-05-11T09:38:01.502Z)

# Description
Cluster: All
Dear users,
We have been notified (again) of a critical security vulnerability. As a result, we must immediately patch our entire infrastructure, including compute and login nodes.
The three login nodes will be rebooted shortly, and the cluster will be put into drain mode.
The cluster will be resumed as soon as the nodes have been patched.
### Update:
For security reasons, we have disabled SSH access to compute nodes until they have all been patched.
We will re-enable this functionality as soon as possible.
## Status: Resolved :green_circle:
start: 2026-05-07T22:00:00Z→2026-05-11T07:00:00Z


## Post 20 by @Yann.Sagon (2026-05-15T15:47:25.283Z)

# Description
Cluster: All
Dear users,
We have been notified (again) of a critical security vulnerability in xdmod. As a result, we stop our xdmod instance until patched.
## Status: Patching done :green_circle:
start: 2026-05-14T22:00:00Z→2026-05-20T07:00:00Z


## Post 21 by @Adrien.Albert (2026-05-21T11:38:24.084Z)

# Description
Cluster: All
Dear users,
We have been notified (again and again) of a critical security vulnerability in on Nvidia driver. As a result, we have drained all gpu node to apply an update.
## Status: Resolved :green_circle:
start: 2026-05-21T10:00:00Z
end:
