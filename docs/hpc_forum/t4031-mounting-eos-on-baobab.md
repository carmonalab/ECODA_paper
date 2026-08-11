# Mounting EOS on baobab

- Source: https://hpc-community.unige.ch/t/4031

- Created: 2025-08-04T13:34:36.943Z

- Posts: 8

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Tomke.Schroeer (2025-08-04T13:34:36.983Z)

what did you try:
I tried to mount eos on baobab by using the following commands:
> (baobab)-[schroeer@login1 ~]$ export EOS_MGM_URL=“root://eosatlas.cern.ch”
> (baobab)-[schroeer@login1 ~]$  export EOS_HOME=“/eos/user”
> (baobab)-[schroeer@login1 ~]$ eos fuse mount eosatlas
what didn’t work:
Mount fails
what was the error message:
```
250804 15:28:41 2037632 cryptossl_X509CreateProxy: unable to load EEC certificate from file: /home/users/s/schroeer/.globus/usercert.pem
warning: assuming you gave a relative path with respect to current working directory => mountpoint=eosatlas
===> Mountpoint   : /home/users/s/schroeer/eosatlas
===> Fuse-Options : fsname=eosatlas.cern.ch:
running eosxd /home/users/s/schroeer/eosatlas -ofsname=eosatlas.cern.ch:
error: failed mount, maybe still mounted? Check with df and eventually 'killall eosd'
```
when running `killall eosd`:
> eosd: no process found
What could be the issue?
Thanks in advance, Tomke


## Post 2 by @Adrien.Albert (2025-08-04T14:49:29.866Z)

Hi @Tomke.Schroeer[@Tomke.Schroeer](https://hpc-community.unige.ch/u/tomke.schroeer)
I don’t know your procedure but is it related to:
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#cvmfs)
### hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#cvmfs)


## Post 3 by @Tomke.Schroeer (2025-08-04T15:27:12.770Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
No, I have no issue accessing these. I was following the instructions given in this post:
eos-on-yggdrasil[eos-on-yggdrasil](https://hpc-community.unige.ch/t/eos-on-yggdrasil/1220)


## Post 4 by @Tomke.Schroeer (2025-08-07T07:32:37.880Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) the problem still persists, is there anything I could try to resolve it?


## Post 5 by @Adrien.Albert (2025-08-07T08:07:45.071Z)

I apologize, I understood that you no longer had any issuefollowing the procedure you shared.
I am not very familiar with EOS and the procedure I need to investigate.


## Post 6 by @Yann.Sagon (2025-08-12T09:43:47.493Z)

@Tomke.Schroeer[@Tomke.Schroeer](https://hpc-community.unige.ch/u/tomke.schroeer)
I’ve investigated more in depth the issue and I’m happy to announce you we found the issue.
First of all, you need to mount the eos filesystem to a location not in your home or scratch as they aren’t standard filesystems.
EOS on Yggdrasil[EOS on Yggdrasil](https://hpc-community.unige.ch/t/eos-on-yggdrasil/1220/3) HPC Announce[HPC Announce](https://hpc-community.unige.ch/c/hpc-announce/6)
> warning it isn’t working anymore if you try to mount the filesystem to “non standard filesystem”. 
For debugging purpose: 
(baobab)-[sagon@cpu001 ~]$ eosxd /home/sagon/eos -ofsname=eosatlas.cern.ch
[...]
fusermount: mounting over filesystem type 0x19830326 is forbidden  <<<<<< here >>>>>

The working solution is to use /tmp/<yourusername>/<yourmointpoint> for example as mounting point.
Second things, we had to create three directories to host the cache with correct rights as if those directories doesn’t exists, they are created automatically with the user/group of the first user running the command to mount the filesystem.
:warning: it isn’t working on Bamboo compute nodes yet due to firewall issue. We are awaiting a change from security team.


## Post 7 by @Tomke.Schroeer (2025-08-14T09:13:53.511Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon),
thanks for investigating! I am not sure I fully understand what I need to do, however, here is what I tried:
```
(baobab)-[schroeer@login1 schroeer]$cd /tmp/schroeer
(baobab)-[schroeer@login1 schroeer]$ export EOS_MGM_URL="root://eosatlas.cern.ch"
(baobab)-[schroeer@login1 schroeer]$ export EOS_HOME="/eos/atlas"
(baobab)-[schroeer@login1 schroeer]$ eos fuse mount eosatlas
```
giving me
```
250814 11:10:14 1727324 cryptossl_X509CreateProxy: unable to load EEC certificate from file: /home/users/s/schroeer/.globus/usercert.pem
warning: assuming you gave a relative path with respect to current working directory => mountpoint=eosatlas
===> Mountpoint   : /tmp/schroeer/eosatlas
===> Fuse-Options : fsname=eosatlas.cern.ch:
running eosxd /tmp/schroeer/eosatlas -ofsname=eosatlas.cern.ch:
error: failed mount, maybe still mounted? Check with df and eventually 'killall eosd'
```
when I then try he command explicitly I get
```
(baobab)-[schroeer@login1 schroeer]$ eosxd /tmp/schroeer/eosatlas -ofsname=eosatlas.cern.ch:
# fsname='eosatlas.cern.ch:'
# -o big_writes enabled
# no config file - running on default values
# no config file for local overwrites
# extracted remote mount dir from fsname is ''
# extracted connection host from fsname is 'eosatlas.cern.ch'
# enabling swapping inodes with md-cache in '/var/tmp/eos/fusex/md-cache/'
warning: sss keytabfile '/home/users/s/schroeer/.eos/fuse.sss.keytab' does not exist - disabling sss/oauth2
# File descriptor limit: 4096 soft, 4096 hard
# allowing max read-ahead buffers of 134217728 bytes
# allowing max write-back buffers of 134217728 bytes
# concurrent eosxd detect enabled, lock prefix /var/tmp/eos-403388/fusex/mount.-tmp-schroeer-eosatlas
mkdir: cannot create directory ‘/var/tmp/eos/fusex/md-cache’: Permission denied
mkdir: cannot create directory ‘/var/tmp/eos/fusex/cache’: Permission denied
mkdir: cannot create directory ‘/var/tmp/eos/fusex/cache’: Permission denied
mkdir: cannot create directory ‘/var/tmp/eos/fusex/credential-store’: Permission denied
error: failed to make path=/var/tmp/eos/fusex/md-cache/eosatlas.cern.ch-tmp-schroeer-eosatlas/899ea20e-78ee-11f0-a79f-74563cec13e8 RWX for root - errno=13
```
Do I need to do something in addition for the second point? Or what am I doing wrong?
Cheers, Tomke


## Post 8 by @Yann.Sagon (2025-08-14T12:22:24.742Z)

Dear @Tomke.Schroeer[@Tomke.Schroeer](https://hpc-community.unige.ch/u/tomke.schroeer), the issue is that you weren’t member of the hpc_users group as your gid is private_dpnc. I’ve added you to the hpc_users group. Can you please logout and login and give a new try?
