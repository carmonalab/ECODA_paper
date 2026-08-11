# Disk quota exceeded

- Source: https://hpc-community.unige.ch/t/3822

- Created: 2025-02-12T17:13:42.259Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Arthur.Freeman (2025-02-12T17:13:42.304Z)

## Primary informations
Username: freemana
Cluster: baobab
## Description
Can’t create files, I get an error saying I exceeded my quota despite making 800 GB of space.
## Steps to Reproduce
I was trying to connect using vscode remote SSH, but was getting an SCP error as it was unable to connect the required files to the remote. I’m able to connect manually by ssh, and if I try to touch a file, I get the error,
```
(baobab)-[freemana@login1 src]$ touch test
touch: cannot touch 'test': Disk quota exceeded
```
So I went ahead and deleted over 800 GB of checkpoints and logs I no longer needed.
The issue however persists.
```
(baobab)-[freemana@login1 src]$ beegfs-ctl --getquota --uid freemana --mount=/srv/beegfs/scratch --connAuthFile=/etc/beegfs/connauthfile

Quota information for storage pool Default (ID: 1):

      user/group     ||           size          ||    chunk files    
     name     |  id  ||    used    |    hard    ||  used   |  hard   
--------------|------||------------|------------||---------|---------
      freemana|388382||      0 Byte|   unlimited||        0| 1000000
```
## Expected Result
I need to be able to create files on the disk.
## Actual Result
I can’t create files even though I freed space, it says my disk quota was exceeded.
Cheers,
A. Freeman


## Post 2 by @Arthur.Freeman (2025-02-13T08:53:06.400Z)

Ok after trying again today, no idea what was done in the background, but I can now create files and vscode is able to connect.
Cheerio.


## Post 3 by @Yann.Sagon (2025-02-13T09:01:17.725Z)

Dear @Arthur.Freeman[@Arthur.Freeman](https://hpc-community.unige.ch/u/arthur.freeman)
Arthur.Freeman:
```
(baobab)-[freemana@login1 src]$ touch test
touch: cannot touch 'test': Disk quota exceeded
```
Where is located `src`? In scratch or home?
Arthur.Freeman:
```
(baobab)-[freemana@login1 src]$ beegfs-ctl --getquota --uid freemana --mount=/srv/beegfs/scratch --connAuthFile=/etc/beegfs/connauthfile

Quota information for storage pool Default (ID: 1):

      user/group     ||           size          ||    chunk files    
     name     |  id  ||    used    |    hard    ||  used   |  hard   
--------------|------||------------|------------||---------|---------
      freemana|388382||      0 Byte|   unlimited||        0| 1000000
```
Here you are looking at the quota on scratch, but we don’t have a quota on file size, only on number of files.


## Post 4 by @Adrien.Albert (2025-02-13T09:08:36.997Z)

hello
I recommend to use:
```
(baobab)-[root@admin1 ~]$ beegfs-get-quota-home-scratch.sh -u freemana
home dir: /home/users/f/freemana
scratch dir: /srv/beegfs/scratch/users/f/freemana

          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
home        |       freemana|388382||  260.92 GiB| 1024.00 GiB||    80487|unlimited
scratch     |       freemana|388382||      0 Byte|   unlimited||        0| 10000000
```
Has a running job generated and deleted data?
