# SSH connection refused

- Source: https://hpc-community.unige.ch/t/3977

- Created: 2025-06-10T12:27:51.095Z

- Tags: baobab

- Posts: 3

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Sara.Ricossa (2025-06-10T12:27:51.140Z)

username: ricossa
cluster: baobab
Since today, I cannot seem to connect to the login node anymore. Either the ssh key is unrecognized, and I’m then asked for a password (the UNIGE account password doesn’t work), or it simply says connection refused.
```
ssh: connect to host login1.baobab.hpc.unige.ch port 22: Connection refused
```
Is this an ongoing issue, or something wrong on my side?
Thank you.


## Post 2 by @Adrien.Albert (2025-06-10T12:30:53.085Z)

Hi @Sara.Ricossa[@Sara.Ricossa](https://hpc-community.unige.ch/u/sara.ricossa)
Please refer to :
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/14) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: Boabab 
We’ve identified a configuration issue with the /tmp directory on login1, which currently limits available space to only 15 GB for all users. This is causing latency and write errors for some process on login node. 
To resolve this, we will perform a short maintenance and reboot the node at 13:30 today. 
This intervention should be quick and will restore proper /tmp storage capacity. 
Thank you for your understanding, 
Update:
2025-06-09T22:00:00Z (UTC) Maintenance s…


## Post 3 by @Sara.Ricossa (2025-06-10T12:31:38.019Z)

Hi! Thank you. I very much apologized. I checked the pinned posts, but didn’t see this one.
