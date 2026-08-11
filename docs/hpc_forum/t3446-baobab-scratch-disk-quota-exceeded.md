# Baobab: scratch: Disk quota exceeded

- Source: https://hpc-community.unige.ch/t/3446

- Created: 2024-05-13T15:50:05.695Z

- Posts: 12

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Malte.Algren (2024-05-13T15:50:05.740Z)

Hi,
When working on scratch, I am getting “Disk quota exceeded” just by simply moving files around. It doesn’t look like I have exceeded the limit of 10.000.000 chunks.
Hope you can help.
Cheers,
Malte
image
image1708×184 8.62 KB
[image1708×184 8.62 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/affd8493aeed712bca376218b84bdc91b2046d09.png)


## Post 2 by @Malte.Algren (2024-05-13T16:20:40.041Z)

Works again now
Cheers,
Malte


## Post 3 by @Paul.Coppin (2024-05-13T21:10:39.911Z)

I am getting the same error :confused:
scratch seems to be 99% filled. Could this be the reason?
```
(baobab)-[coppinp@login2 scratch]$ df -hT
beegfs_scratch                                       beegfs    1.5P  1.5P   20T  99% /srv/beegfs/scratch
```
Update, it seems that it was indeed because scratch filled up.
I created some space, and am currently moving some ~20TB of communal data used by our entire research group (DAMPE) to another storage. Hopefully this should help.


## Post 4 by @Debajyoti.Sengupta (2024-05-14T07:28:21.724Z)

Same here.
My job crashed with the same error.
Attaching my beegfs quota output below.
image
image943×212 7.74 KB
[image943×212 7.74 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/9a6a1f43da0fef6752204b684e42b59843d11a66.png)


## Post 5 by @Paul.Coppin (2024-05-14T19:31:05.141Z)

Update, since this afternoon it seems that there is 78T (5%) of free disk space.
Yet I still got jobs which failed at 21:05:08 this evening with the error:
`IOError: [Errno 122] Disk quota exceeded`
@HPC[@HPC](https://hpc-community.unige.ch/u/hpc) Any idea?


## Post 6 by @Malte.Algren (2024-05-15T05:56:14.392Z)

@Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
The problem still exists
Cheers,
Malte


## Post 7 by @Adrien.Albert (2024-05-15T07:09:41.285Z)

Hi @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren), @Paul.Coppin[@Paul.Coppin](https://hpc-community.unige.ch/u/paul.coppin)
Currently the scratch filesystem seems to have enough space both in use and inode
2024-05-15T07:04:00Z
```
(baobab)-[root@login2 ~]$ df -t beegfs -h
Filesystem      Size  Used Avail Use% Mounted on
beegfs_home     138T  116T   22T  85% /home
beegfs_scratch  1.5P  1.4P   83T  95% /srv/beegfs/scratch
```
however, we’re not immune to the possibility that the aggreagat of current jobs may temporarily fill the scratch.
Are you sure your error appeared after the space was freed up?


## Post 8 by @Malte.Algren (2024-05-15T07:23:58.813Z)

Well, I think that is the problem when running long jobs that save files from scratch. It might be temporarily filled and your job will crash.
But no at the moment I am not getting the error.
image
image1340×621 37.1 KB
[image1340×621 37.1 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/67b0085245e259d50e57c68176734335af4ad967.png)


## Post 9 by @Malte.Algren (2024-05-16T06:59:43.907Z)

Still seeing the issue where scripts can run for some time but eventually, I get “IOError: [Errno 122] Disk quota exceeded”.
Any way to add additional storage to the scratch?
image
image2001×859 55.2 KB
[image2001×859 55.2 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/6b597efd57bcd1c1b2122a86518337e73d99ff34.png)


## Post 10 by @Paul.Coppin (2024-05-16T07:06:05.306Z)

Same here, almost all my jobs suddenly failed again yesterday evening (2024-05-15) because of the disk quota error :confused:


## Post 11 by @Gael.Rossignol (2024-05-17T08:01:31.733Z)

Malte.Algren:
> Hi,
> When working on scratch, I am getting “Disk quota exceeded” just by simply moving files around. It doesn’t look like I have exceeded the limit of 10.000.000 chunks.
> Hope you can help.
> Cheers,
> Malte
> image
> image1708×184 8.62 KB
> [image1708×184 8.62 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/affd8493aeed712bca376218b84bdc91b2046d09.png)
Hello,
@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) found the problem! For share “private_dpnc” quota for number of files was reached :
```
(baobab)-[algren@admin1 ~]$ beegfs-get-quota-home-scratch.sh -g private_dpnc

          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
home        |   private_dpnc|  1014||   12.76 TiB|   unlimited|| 24719933|unlimited
scratch     |   private_dpnc|  1014||  489.86 TiB|  875.00 TiB|| 69942301| 70000000
```
As DPNC has dedicated storage I will update this value to 100000000.
Best regards,


## Post 12 by @Paul.Coppin (2024-05-17T12:42:20.508Z)

@Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) Thanks a lot for digging out the issue and updating the quota!
