# Baobab: Remote I/O error with bash script

- Source: https://hpc-community.unige.ch/t/4058

- Created: 2025-08-22T13:20:55.933Z

- Tags: baobab

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Laure.Moinat (2025-08-22T13:20:55.996Z)

user: moinatl
Dear team,
It seems that Baobab has a problem making connections with certain files when running a bash script. I’m experiencing an error called ‘Remote I/O error,’ and I don’t understand why. Do you have an explanation?
Best,
Laure
PS: we are several to experience this error


## Post 2 by @Paul.Coppin (2025-08-22T13:39:10.994Z)

Hi,
we also see odd things:
`(baobab)-[coppinp@login1 ~]$ ls -l /srv/beegfs/scratch/groups/dpnc/`
`ls: cannot read symbolic link '/srv/beegfs/scratch/groups/dpnc/dampe': ``Remote I/O error`
from what I remember, ‘/srv/beegfs/scratch/groups/dpnc/dampe’ is a simlink to ‘srv/beegfs/dpnc/groups/dampe/’
but the directory is there and mounted:
`ls -l /srv/beegfs/dpnc/groups/dampe/`
`total 2`
`drwxrwx---  5 usr_s_dampe1 private_dpnc  3 Aug 13 11:53 prod`
`drwxrwxrwx 24 usr_s_dampe1 private_dpnc 27 Jul 15 10:03 public`
`drwxrwx---  8 usr_s_dampe1 private_dpnc  9 Jun  3 01:42 users`
Cheers,
Paul


## Post 3 by @Mathias.Elbaz (2025-08-22T13:52:19.863Z)

Hi,
Similar error here :
cp: error reading ‘/srv/beegfs/scratch/shares/sanchezf/…’: Remote I/O error
Best,
Mathias


## Post 4 by @Sophie.Slaats (2025-08-22T13:53:46.605Z)

Yes, the same error here.
hpc-remote-io-error
hpc-remote-io-error852×203 60.2 KB
[hpc-remote-io-error852×203 60.2 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/871d0525005e953a9917f13a1b1bd59862ae1ab3.png)


## Post 5 by @Yann.Sagon (2025-08-22T14:38:03.364Z)

Hello,
I’m checking. I’ll update this post.


## Post 6 by @Oscar.Dabrowski (2025-08-22T16:06:08.007Z)

Hello,
I am also experiencing I/O errors with my jobs on Baobab. It does not seem to be able to write data to scratch or when it does it is extremely slow. I also noticed that suddenly all my jobs (a dozen) crashed (almost at the same time) and the logs showed errors like “ Unable to write file … “
Best,
Oscar.


## Post 7 by @Adrien.Albert (2025-08-22T16:44:22.932Z)

Hello,
The issue is related to [2025] Current issues on HPC Cluster - #20 by Yann.Sagon[[2025] Current issues on HPC Cluster - #20 by Yann.Sagon](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/20) and has been resolved.
