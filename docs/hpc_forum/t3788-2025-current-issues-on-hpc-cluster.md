# [2025] Current issues on HPC Cluster

- Source: https://hpc-community.unige.ch/t/3788

- Created: 2025-01-20T09:06:04.069Z

- Tags: all

- Posts: 27

- Category: 9

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-01-20T09:06:04.178Z)

Dear users,
We have some problems with login1 on Bamboo, we are investigation to find a solution.
For the time being, running jobs are continuing, but it isn’t possible to login to the login node.
We’ll keep you informed, thank you for your understanding.
Kind regards
edit: problem solved, we’ll add the login node to our monitoring server.
### Status : Solved :green_circle:
start: 2025-01-18T10:00:00Z
end: 2025-01-20T08:00:00Z


## Post 2 by @Adrien.Albert (2025-01-24T16:14:38.886Z)

Dear users,
We are currently experiencing an issue with our storage server (home)  that may result in the following error message:
```
I/O remote error
```
We are investigating to determine the root cause.
### Status : Resolved :green_circle:
start: 2025-01-23T10:00:00Z
end:


## Post 3 by @Gael.Rossignol (2025-01-27T10:53:23.332Z)

Dear users,
Rhis morning, we had an issue with login node on Yggdrasil. Some process crash server but after reboot all is now up and running.
Root cause may be linked to I/O errors on storage.
Sorry for inconvenience.
### Status : Resolved :green_circle:
start: 2025-01-27T08:45:00Z
end: 2025-01-27T10:30:00Z


## Post 4 by @Adrien.Albert (2025-03-13T14:49:42.964Z)

## Description
We encountered an issue with OpenOnDemand authentication. While attempting to fix outsider access in collaboration with the Authentication team, the authentication rules were affected. As a result, some users with dual identities (Collaborator/Student) may have experienced account mismatches with the HPC system.
## Resolution
The fix has been rolled back, and the issue should no longer be present. We plan to test an alternative solution to allow outsider access to OpenOnDemand.
### Status : Resolved :green_circle:
start: 2025-03-12T08:45:00Z
end: 2025-03-13T07:30:00Z


## Post 5 by @Adrien.Albert (2025-03-13T23:30:39.973Z)

## Description
We encountered an issue with Slurm, the database is unreachable resulting errors executing slurm command (sinfo, sacct etc..)
## Resolution
Database service has been restarted.
### Status : Resolved :green_circle:
start: 2025-03-12T15:00:00Z
end: 2025-03-13T22:30:00Z


## Post 6 by @Adrien.Albert (2025-03-20T09:24:06.612Z)

## Description
We encountered an issue where Baobab login1 was stuck due to a Jobs user. The node became completely unreachable, preventing us from identifying the responsible user.
## Resolution
Node has been rebooted
### Status : Resolved :green_circle:
2025-03-20T09:20:00Z→2025-03-20T09:23:00Z


## Post 7 by @Yann.Sagon (2025-04-02T15:03:08.250Z)

## Description
We encountered an issue with Slurm on Baobab, the database is unreachable resulting errors executing slurm command (sinfo, sacct etc…) and jobs terminating too early with TIMEOUT reason.
## Resolution
There was a network issue between slurm and slurmdbd. It is now solved, thanks for your understanding.
### Status : Resolved :green_circle:
start: 2025-04-02T12:00:00Z
end: 2025-04-02T14:00:00Z


## Post 8 by @Yann.Sagon (2025-04-04T14:22:18.912Z)

## Description
We encountered an issue with admin servers and Slurm on Bamboo, the database is unreachable resulting errors executing slurm command (sinfo, sacct etc…) and jobs terminating too early with TIMEOUT reason.
## Resolution
We are working on it.
### Status : solved :green_circle:
start: 2025-04-04T14:00:00Z
start: 2025-04-04T14:54:00Z


## Post 9 by @Adrien.Albert (2025-04-09T11:52:54.387Z)

## Description
Cannot connect to yggdrasil[Cannot connect to yggdrasil](https://hpc-community.unige.ch/t/cannot-connect-to-yggdrasil/3912) HPC Support[HPC Support](https://hpc-community.unige.ch/c/hpc-support/13)
> I currently cannot connect to yggdrasil, using either X2Go or standard ssh.  I have confirmed with other users that this issue is not specific to me.  Is yggdrasil down right now?
### Status : solved :green_circle:
start: 2025-04-09T10:00:00Z
start: 2025-04-04T11:50:00Z


## Post 10 by @Yann.Sagon (2025-04-24T12:19:50.405Z)

## Description
We encountered an issue with infiniband network on Bamboo, the home and is unreachable and every node is in drain state.
## Resolution
There was an incoherent parameter in one of the configuration file since the latest maintenance and the effect was triggered today. It is now solved.
### Status : solved :green_circle:
start: 2025-04-24T09:30:00Z
start: 2025-04-24T12:15:00Z


## Post 11 by @Adrien.Albert (2025-04-24T16:05:15.523Z)

## Description
`Cluster: Baobab`
The login1 node was temporarily unavailable on  due to an issue caused by excessive memory usage from a user session.
As a result, the node became unresponsive and could not be accessed for a the period of time. The issue has since been resolved, and the system is back to normal.
### Status : solved :green_circle:
start: 2025-04-21T19:00:00Z
end : 2025-04-22T08:20:00Z


## Post 12 by @Adrien.Albert (2025-04-30T08:29:28.704Z)

## Description
`Cluster: baobab`
An unexpected behavior temporarily impacted the automatic sending of outgoing emails on one of our servers. The issue has been resolved, and email delivery is functioning normally. No email loss has been detected during the incident.
### Status : solved :green_circle:
start: 2025-04-29T09:00:00Z
end : 2025-04-30T08:21:00Z


## Post 13 by @Yann.Sagon (2025-05-26T13:22:21.265Z)

## Description
`Cluster: multi cluster`
We are reinstalling slurmdbd and for this reason commands such as sacct won’t work during the reinstallation. Thanks for your understanding.
### Status : solved :green_circle:
start: 2025-05-26T13:17:00Z
end : 2025-05-26T14:10:00Z


## Post 14 by @Adrien.Albert (2025-06-10T10:30:22.802Z)

## Description
`Cluster: Boabab`
We’ve identified a configuration issue with the `/tmp` directory on `login1`, which currently limits available space to only 15 GB for all users. This is causing latency and write errors for some process on `login node`.
To resolve this, we will perform a short maintenance and reboot the node at 13:30 today.
This intervention should be quick and will restore proper `/tmp` storage capacity.
Thank you for your understanding,
### Update:
2025-06-09T22:00:00Z Maintenance started at 14:00 PM and finished at 14:30 PM without any issue.
### Status : Resolved :green_circle:
start: 2025-06-10T06:00:00Z
end :2025-06-10T12:30:00Z


## Post 15 by @Yann.Sagon (2025-06-25T12:01:10.142Z)

## Description
`Cluster: multi cluster`
Under specific circumstances, it is possible that you will be able to use a resource such as a GPU that is already associated with another user. This may happen when you request resources using `salloc` and then connect to the compute node using `ssh` later. Alternatively, Slurm may allocate you a resource that is already in use by a previous job. This issue has already been resolved, but it will only be fully resolved once all resources are released by the running jobs.
### Status : solved :green_circle:
start: 2025-06-25T07:37:00Z
end : 2025-06-25T11:45:00Z


## Post 16 by @Adrien.Albert (2025-07-07T14:57:50.091Z)

## Description
Cluster: Bamboo
We are currently experiencing an issue with the Bamboo scratch and Home storage. This share is temporarily unavailable.
Our team is actively investigating the problem, and we will provide an update as soon as more information becomes available.
We apologize for the inconvenience and appreciate your patience.
### Status : alert :green_circle:
start: 2025-07-07T14:50:00Z
end : 2025-07-07T15:10:00Z
# Update
Filesystem impacted scratch:
start: 2025-07-07T20:00:00Z
end : 2025-07-09T14:20:00Z
The scratch has been re-activated. However, despite our investigation, the root cause of the issue has not yet been determined. We remain vigilant and will continue to monitor the situation closely.


## Post 17 by @Adrien.Albert (2025-07-31T14:37:21.655Z)

## Description
Cluster: Baobab
We are currently experiencing load on Boabab login node slowing the server. The server has been rebooted.
We kindly but strongly remind all users that only light, non-computational tasks must be run on the login node.
### Status : alert :green_circle:
start: 2025-07-31T14:20:00Z
end : 2025-07-31T14:35:00Z


## Post 18 by @Adrien.Albert (2025-08-04T10:08:25.983Z)

## Description
Cluster: Baobab
We experienced an issue on the scratch storage due to a timeout on one of the disks. This caused unresponsive data for a short period, which in turn led to a temporary crash of the BeeGFS service.
Some users may have noticed latency on the scratch or error messages during this time. The concerned disk has been removed from the pool, and the issue is now resolved.
Thank you for your understanding.
–
HPC Team
### Status : Resolved  :green_circle:
start: 2025-08-01T13:00:00Z
end: 2025-08-04T09:00:00Z


## Post 19 by @Gael.Rossignol (2025-08-12T07:35:22.164Z)

## Description
Cluster: Baobab
Last night, we encountered an issue affecting the home storage system, caused by an unexpected stoppage of the meta services. As a result, some jobs running at the time, as well as those triggered during this period, may have failed.
The services have now been fully restored, and all systems are operating normally. Our team is currently investigating the root cause to prevent similar occurrences in the future.
We appreciate your understanding.
–
HPC Team
### Status : Resolved  :green_circle:
start: 2025-08-11T15:15:00Z
end: 2025-08-12T08:23:00Z


## Post 20 by @Yann.Sagon (2025-08-22T14:56:04.921Z)

## Description
Cluster: Baobab
One of the scratch server and dpnc server crashed and we are unable to restart it. We need to go onsite to check what is going on.
In the meantime, scratch and dpnc storage are unavailable until further notice.
edit: this is fixed. I had to remove the RAID cards from the server and remove the backup battery to be able to restart the server again with the RAID card. It is now working as usual.
–
HPC Team
### Status : Resolved :green_circle:
start: 2025-08-22T12:15:00Z
end: 2025-08-22T16:43:00Z
