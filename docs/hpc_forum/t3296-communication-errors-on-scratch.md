# Communication errors on scratch

- Source: https://hpc-community.unige.ch/t/3296

- Created: 2024-02-08T16:46:14.293Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Andrea.Serpolla (2024-02-08T16:46:14.349Z)

Hello,
I’m having troubles accessing files on `scratch`.
This is what I get if I try to `cat` a generic logfile of one of my running jobs:
```
(baobab)-[serpolla@cpu007 ~]$ cat /srv/beegfs/scratch/users/s/serpolla/__omissis__.out 
cat: /srv/beegfs/scratch/users/s/serpolla/__omissis__.out: Communication error on send
```
and this is what I get if I try to check my quota:
```
(baobab)-[serpolla@cpu007 ~]$ beegfs-get-quota-home-scratch.sh 
home dir: /home/users/s/serpolla
scratch dir: /srv/beegfs/scratch/users/s/serpolla

        user/group                 ||           size          ||    chunk files
storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
----------------------------|------||------------|------------||---------|---------
home        |       serpolla|427178||    4.84 GiB| 1024.00 GiB||    88421|unlimited
Could not gather quota information from all storage nodes.
Refusing to display incorrect information.
scratch     | (0) 17:44:02 Main [GetQuotaInfo.cpp:252] >> Unable to fetch quota information from storage target. Is the node offline? Storage target id: 1306
```
I’m logged in a shell of a compute node since on the login node my shell was just freezing trying to perform these operations.
Thank you for the help.
Andrea


## Post 2 by @Yann.Sagon (2024-02-09T08:09:53.604Z)

Hi,
this is solved: [2024] Current issues on HPC Cluster - #4 by Yann.Sagon[[2024] Current issues on HPC Cluster - #4 by Yann.Sagon](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/4)
