# Issue :List Accessible shared directories not working

- Source: https://hpc-community.unige.ch/t/4075

- Created: 2025-08-31T10:58:37.271Z

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2025-08-30T07:49:55.769Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert) ,
when trying to rsync a file to Baobab, Yggdrasil, or Bamboo, I get the following error:
`protocol version mismatch -- is your shell clean?`
`(see the rsync manpage for an explanation)`
`rsync error: protocol incompatibility (code 2) at compat.c(626) [sender=3.4.1]`
The issues started this morning (yesterday everything still worked fine).
When doing some Googling, I saw that it could be related to: “your shell’s login outputting stuff on a non-interactive shell.” and indeed:
`ssh baobab "/bin/true" > testfile`
`pcoppin:HE$ cat testfile`
`======================================================================`
`                  Accessible Shared Directories                   `
`======================================================================`
`Directories                            Group`
`-------------------------------------  ------------`
`/srv/beegfs/scratch/groups/dpnc        private_dpnc`
`/srv/beegfs/scratch/groups/dpnc/dampe  share_dampe`
While this test should return an empty file.
Is there a way to output this message only on interactive shells?


## Post 2 by @michael.divia (2025-08-30T11:37:14.961Z)

We cannont `scp` eather anymore :
`(base) padi@TheBeast:~/Git/Projet de Bachelor/Project_Source$ scp -r CaloJetSSD divia@login1.baobab.hpc.unige.ch:/home/users/d/divia/`
`scp: Received message too long 171785533`
`scp: Ensure the remote shell produces no output for non-interactive sessions.`


## Post 3 by @Adrien.Albert (2025-08-31T10:58:37.271Z)

hi @Paul.Coppin[@Paul.Coppin](https://hpc-community.unige.ch/u/paul.coppin)
I reverted the changes while waiting for Monday. I had configured it to run only in interactive mode… well, I made a mistake.
My apologies for the inconvenience :pray:


## Post 4 by @Adrien.Albert (2025-09-01T15:43:28.668Z)

Hello @Paul.Coppin[@Paul.Coppin](https://hpc-community.unige.ch/u/paul.coppin)
I deployed the fix on Bamboo, and it works for me. Can you try it and let me know if it works for you?
```
[root@localhost ~]# rsync -avh slurm_conf.tgz alberta@login1.bamboo.hpc.unige.ch:
sending incremental file list

sent 52 bytes  received 12 bytes  42.67 bytes/sec
total size is 6.88K  speedup is 107.42
[root@localhost ~]# ssh alberta@login1.bamboo.hpc.unige.ch echo toto
toto
[root@localhost ~]# ssh alberta@login1.bamboo.hpc.unige.ch 
Last login: Mon Sep  1 17:41:15 2025 from 10.40.126.9
 ____                  _
|  _ \                | |
| |_) | __ _ _ __ ___ | |__   ___   ___
|  _ < / _` | '_ ` _ \| '_ \ / _ \ / _ \
| |_) | (_| | | | | | | |_) | (_) | (_) |
|____/ \__,_|_| |_| |_|_.__/ \___/ \___/
                 _             _      __ 
                | |           (_)    /_ |
                | | ___   __ _ _ _ __ | |
                | |/ _ \ / _` | | '_ \| |
                | | (_) | (_| | | | | | |
                |_|\___/ \__, |_|_| |_|_|
                          __/ |          
                         |___/  

 Documentation: https://doc.eresearch.unige.ch/hpc/start
 Forum: https://hpc-community.unige.ch/
 OpenOndemand: https://openondemand.baobab.hpc.unige.ch/
 support: https://doc.eresearch.unige.ch/hpc/start#support_-_get_help

======================================================================
                    Accessible Shared Directories                     
======================================================================

Directories                            Group
-------------------------------------  ------------
/srv/beegfs/scratch/groups/dpnc        private_dpnc
/srv/beegfs/scratch/groups/dpnc/atlas  share_atlas
/srv/beegfs/scratch/groups/rodem       share_rodem
/srv/fast/share/rodem/                 share_rodem
```


## Post 5 by @Paul.Coppin (2025-09-12T08:56:37.796Z)

Dear @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
Indeed, it works fine now.
Sorry for the very late reply! I just went into the profile-settings to turn on email notifications for all mentions to catch new message more quickly in the future.
All the best,
Paul


## Post 6 by @Adrien.Albert (2025-09-12T12:21:38.509Z)

Perfect ! Thank you for the feed back :pray:
