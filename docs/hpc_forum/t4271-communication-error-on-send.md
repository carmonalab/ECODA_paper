# Communication error on send

- Source: https://hpc-community.unige.ch/t/4271

- Created: 2026-04-02T22:12:16.987Z

- Tags: bamboo

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @pablo.strasser1 (2026-04-02T22:12:17.041Z)

I’m unable to access my home on bamboo. When I login by ssh or when I cd to the directory, I get the following error:
Could not chdir to home directory /home/users/s/strassep: Communication error on send
-bash: /home/users/s/strassep/.bash_profile: Communication error on send


## Post 2 by @Raphael.Rubino (2026-04-03T05:22:56.745Z)

Same issue for me.
`ls: cannot access ‘/home/users/r/rubinor/’: Communication error on send`


## Post 3 by @Lucille.Delisle1 (2026-04-03T06:49:58.662Z)

The issue seems to be both for home and shared.
```
(bamboo)-[delislel@login1 /]$ ls /srv/beegfs/scratch/users/d/delislel 
ls: cannot access '/srv/beegfs/scratch/users/d/delislel': Communication error on send
(bamboo)-[delislel@login1 /]$ ls $HOME
ls: cannot access '/home/users/d/delislel': Communication error on send
(bamboo)-[delislel@login1 /]$ 
```
It seems similar to another issue that happened in March 2021: Current issues on Baobab and Yggdrasil - #37 by Luca.Capello[Current issues on Baobab and Yggdrasil - #37 by Luca.Capello](https://hpc-community.unige.ch/t/current-issues-on-baobab-and-yggdrasil/787/37)
Unfortunately this is Easter week-end… We probably need to use other clusters (no issue on baobab nor yggdrasil) up to Tuesday.


## Post 4 by @Lucille.Delisle1 (2026-04-03T07:06:23.971Z)

There is at least an issue with home02_storage01:
monitor.hpc.unige.ch[monitor.hpc.unige.ch](https://monitor.hpc.unige.ch/d/tf-bgennzz/beegfs-storage-server-bamboo?orgId=1&from=now-24h&to=now&var-storageID=home02_storage01)
### Grafana[Grafana](https://monitor.hpc.unige.ch/d/tf-bgennzz/beegfs-storage-server-bamboo?orgId=1&from=now-24h&to=now&var-storageID=home02_storage01)
image
image1252×801 95.6 KB
[image1252×801 95.6 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/8ea90dad0c47bf69dc66a9d21228fa7404a5876d.png)
And scratch02_storage01
monitor.hpc.unige.ch[monitor.hpc.unige.ch](https://monitor.hpc.unige.ch/d/tf-bgennzcv/beegfs-storage-server-bamboo-scratch?orgId=1&from=now-24h&to=now&var-storageID=scratch02_storage01)
### Grafana[Grafana](https://monitor.hpc.unige.ch/d/tf-bgennzcv/beegfs-storage-server-bamboo-scratch?orgId=1&from=now-24h&to=now&var-storageID=scratch02_storage01)
image
image1252×801 70.2 KB
[image1252×801 70.2 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/103985d51e41b36a49a5ccd02b38d2a9d1326283.png)


## Post 5 by @Ghasem.Hajianfar (2026-04-03T09:35:12.481Z)

Same error for me:
Could not chdir to home directory /home/users/: Communication error on send
/usr/bin/xauth:  error in locking authority file /home/users/.Xauthority
-bash: /home/users/.bash_profile: Communication error on send


## Post 6 by @Jaspreet.Saini (2026-04-03T14:49:30.685Z)

Same issue for me communication error.


## Post 7 by @Adrien.Albert (2026-04-05T10:10:00.148Z)

Super Cat GIFs  Tenor


## Post 8 by @Adrien.Albert (2026-04-05T10:46:29.065Z)

Dear HPC users,
There is currently a problem with the storage server. Since it is a hardware issue, I cannot intervene until Tuesday.
In the meantime, please use Baobab or Yggdrasil.
We apologize for the inconvenience caused.
Best Regards


## Post 9 by @Yann.Sagon (2026-04-07T06:23:17.922Z)

Hello this is solved. Of course this kind of s**t always happens during vacation :frowning:
[2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/16) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Cluster: Bamboo 
Description
An incident is currently ongoing on the Bamboo cluster. One of our storage server is experiencing an issue impacting the scratch and home storage. 
Status: Solved green_circle
Start: 2026-04-02T20:24:00Z (UTC) 
End: 2026-04-07T06:00:00Z (UTC)
