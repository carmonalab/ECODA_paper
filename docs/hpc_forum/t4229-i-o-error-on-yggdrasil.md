# I/O error on yggdrasil

- Source: https://hpc-community.unige.ch/t/4229

- Created: 2026-02-23T10:23:01.072Z

- Tags: yggdrasil

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Thibault.Garel (2026-02-23T10:23:01.137Z)

Dear HPC support team,
there seems to be a problem with the filesystem on yggdrasil, I am having this kind of error :
$ touch test_touch
touch: setting times of ‘test_touch’: Remote I/O error
Thanks.


## Post 2 by @Sebastien.Miche (2026-02-23T10:38:52.423Z)

Hello,
I am having the same issue whenever I try to save a file / script.
Best!


## Post 3 by @Adrien.Albert (2026-02-23T11:00:50.003Z)

Hello @Thibault.Garel[@Thibault.Garel](https://hpc-community.unige.ch/u/thibault.garel)
Thank you for your report, here the related post:
[2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/11) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear Users, 
Cluster: Yggdrasil 
Description
An issue is currently affecting the Home storage on the Yggdrasil cluster. 
Some users may encounter the following message: 
“I/O remote error” 

Possible disruptions accessing the Home directory
Risk of input/output errors during commands or job execution
The cluster remains available, but some operations may fail

We are working to restore normal service as soon as possible. 
An update will be posted as soon as we have more information. 
We apologiz…
