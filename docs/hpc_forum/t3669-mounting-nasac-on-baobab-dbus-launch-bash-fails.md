# Mounting nasac on baobab: dbus-launch bash fails

- Source: https://hpc-community.unige.ch/t/3669

- Created: 2024-10-02T13:06:56.306Z

- Tags: baobab

- Posts: 28

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthieu.Stigler (2024-10-02T13:06:56.341Z)

I tried to mount a new nasac share on baoab, and following instructions in hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#access_external_storage),
I run:
```
$dbus-launch bash
```
but this triggers the error message:
> dbus-daemon[154590]: Failed to start message bus: Failed to bind socket “/tmp/eb-240s9w15/dbus-BhXJEhGkTK”: No such file or directory
> EOF in dbus-launch reading address from bus daemon
This is happening on a cluster node on baobab.
Also, the documentation is not very clear on `smb://server_name/share_name`. Could you give a typical example in the documentation? I asked for a standard NASAC share and got `J'ai créé le partage `\isis.unige.ch\nasac\gsem\name`... what is `server_name`and`share_name` here?
Thanks!


## Post 2 by @Yann.Sagon (2024-10-04T13:01:37.021Z)

Dear @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
this is unfortunately a know issue [2024] Current issues on HPC Cluster - #15 by Adrien.Albert[[2024] Current issues on HPC Cluster - #15 by Adrien.Albert](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/15) still not solved.
If you want to copy the data from NASAC  to the cluster, you may want to try to  mount the cluster storage in your desktop/server and copy from there. Not very efficient I must admit:(
Best
Yann


## Post 3 by @Adrien.Albert (2024-10-04T13:12:46.306Z)

Hello @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
To briefly summarize, the issue is caused by a Mellanox driver update that disables the CIFS kernel module required for mounting smb share.
Here the mention in officiel Mellanox known issues 2657392[2657392](https://docs.nvidia.com/networking/display/mlnxofedv23103220lts/known+issues):
> Description: OFED installation caused CIFS to break in RHEL 8.4 and above. A dummy module was added so that CIFS would be disabled after the OFED installation in RHEL 8.4 and above.
We tried to apply a workaround based on the very limited information from Mellanox-NVIDIA, but it didn’t work. Mellanox-NVIDIA support has been mostly quiet on this issue.
I’m currently working on another clue, but I can’t guarantee a solution, as dealing with the kernel can be quite tricky.


## Post 4 by @Matthieu.Stigler (2024-10-07T07:31:52.016Z)

Thanks for following-up on this, I really appreciate! Hope a workaround can be found and especially that Mellanox responds to this!


## Post 5 by @Adrien.Albert (2024-10-08T22:28:07.586Z)

Hi @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
Few chance to get an answer from mellanox as they assume that is not their scope…
But the good news is :tada: mywork arround is working, we need to test the robustness of this patch.
If you are interested in testing this patch, please send an email to HPC.
PS: Note that the module has been deactivated by mellanox for a good reason, so we’re not immune to unexpected behaviour.


## Post 6 by @Alexis.Hervais-Adelman (2024-10-10T06:30:45.550Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
Thanks for the update in addressing this pernicious issue.
That’s good news for all users of the HPC who also depend upon the NASAC for sharing data across our teams.
I’d like to try your patch, but it’s not clear which HPC email address you would like us to write to (hpc or hpc-admin…?)
More generally - Given this appears to be an essential functionality, but you say it has been disabled by Mellanox for good reason, what alternative best practice can you propose for sharing data sustainably?
Thanks in advance!
Alexis


## Post 7 by @Adrien.Albert (2024-10-17T08:06:24.803Z)

Hi All
here the new about the work arround :wink:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/15) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> All Clusters: GIO mount storage from NASAC
Dear User, 
We have recently discovered issues with mounting storage spaces using the CIFS/SMB protocol from NASAC. 
The problem comes from the update of Mellanox suite of packages, which manages our Infiniband (fast) network. Unfortunately, these packages disable the CIFS module, preventing CIFS/SMB storage mounts. 
Since the Mellanox suite is crucial for the cluster’s operation, we cannot remove it. Currently, Mellanox doesn’t provides a workaround. 
W…


## Post 8 by @Matthieu.Stigler (2024-10-17T13:14:31.235Z)

great, thanks a lot Adrien!


## Post 9 by @Matthieu.Stigler (2024-10-17T13:15:55.444Z)

were you able to deploy the patch on the on demand docker images too, or is that still an issue? Thanks!


## Post 10 by @Adrien.Albert (2024-10-17T16:48:18.032Z)

As the patch could be unstable, we have limited its deployment to the login node to ensure there is no impact on production.
Matthieu.Stigler:
> Were you able to deploy the patch on the on-demand Docker images as well, or is that still an issue? Thanks!
I tried some basic adjustments, but unfortunately they were not successful.
As a workaround has been implemented, we won’t be exploring this fix any further. While I understand the convenience of directly mounting the share in Singularity images on compute nodes, the time investment required to avoid a simple data copy on the cluster is too significant. Furthermore, I was unable to find any relevant documentation on this approach.
As an alternative, we’ve discussed converting your CIFS share to NFS. This would allow us to easily mount your share across all compute nodes. Have you see with the storage team which is the best option?


## Post 11 by @Adrien.Albert (2024-12-09T10:55:24.190Z)

Hello @Matthieu.Stigler[@Matthieu.Stigler](https://hpc-community.unige.ch/u/matthieu.stigler)
Good News! :tada:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/15) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> All Clusters: GIO mount storage from NASAC
Dear User, 
We have recently discovered issues with mounting storage spaces using the CIFS/SMB protocol from NASAC. 
The problem comes from the update of Mellanox suite of packages, which manages our Infiniband (fast) network. Unfortunately, these packages disable the CIFS module, preventing CIFS/SMB storage mounts. 
Since the Mellanox suite is crucial for the cluster’s operation, we cannot remove it. Currently, Mellanox doesn’t provides a workaround. 
W…
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/15)
> 2024-12-04T23:00:00Z
> Since Rocky9 was deployed on Bamboo during the last maintenance, the module now seems to be available again. You should be able to mount your NASAC share on both the login nodes and compute nodes.
> Rocky9 will also be deployed on Baobab and Yggdrasil during the next maintenances.


## Post 12 by @Julien.Prados (2024-12-18T09:51:10.048Z)

Dear HPC Team,
We are trying to mount the NASAC on `bamboo` using the same procedure we are using on `baobab`. The procedure is working on `baobab` and we can access our data in `~/.gvfs/` ), but we are unable to localise the mount point when the procedure is repeated on `bamboo`.
More specifically, we are using these commands:
```
dbus-launch bash
gio mount 'smb://ISIS;prados@nasac-m2.unige.ch/m-BioinfoSupport'
```
The authentification process and mounting process seems to went well if I refer to running processes:
```
(bamboo)-[prados@login1 ~]$ ps -u prados
    PID TTY          TIME CMD
3749519 ?        00:00:00 sshd
3749520 pts/13   00:00:00 bash
3764679 pts/13   00:00:00 bash
3764683 ?        00:00:00 dbus-daemon
3764753 ?        00:00:00 gvfsd
3764760 ?        00:00:00 gvfs-udisks2-vo
3764768 ?        00:00:00 gvfsd-smb
3765448 pts/13   00:00:00 ps
```
However we are not able to find the mount point ! The files are note accessible in usual locations (`~/.gvfs`, and `/var/run/user/`), and the mount point is not listed by `findmnt`.
Could you please help us to localise the mount point ?
Thanks a lot,
Julien


## Post 14 by @Adrien.Albert (2024-12-18T23:11:45.306Z)

```
(bamboo)-[alberta@login1 ~]$ dbus-launch bash
(bamboo)-[alberta@login1 ~]$ gio  mount   smb://isis.unige.ch/nasac/hpc_exchange/backup < .credentials 
Authentication Required
Enter user and password for share “nasac” on “isis.unige.ch”:
User [alberta]: Domain [SAMBA]: Password: 

(bamboo)-[alberta@login1 ~]$ gio mount --list
Drive(0): SAMSUNG MZ7L3480HBLT-00A07
  Type: GProxyDrive (GProxyVolumeMonitorUDisks2)
Drive(1): SAMSUNG MZ7L3480HBLT-00A07
  Type: GProxyDrive (GProxyVolumeMonitorUDisks2)
Mount(0): nasac on isis.unige.ch -> smb://isis.unige.ch/nasac/
  Type: GDaemonMount
```
Following the documentation: hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc?s%5B%5D=gio#where_does_gio_mounts_my_data)
```
(bamboo)-[alberta@login1 ~]$  ps ux | grep -e '[g]vfsd-fuse'
alberta   526981  0.0  0.0 349216  2048 ?        Sl   Dec05   0:00 /usr/libexec/gvfsd-fuse /run/user/401775/gvfs -f

(bamboo)-[alberta@login1 ~]$ ls /run/user/401775/gvfs/smb-share\:server\=isis.unige.ch\,share\=nasac/hpc_exchange/backup/
titi  toto
```
With gio command:
```
(bamboo)-[alberta@login1 ~]$ gio mount --list smb://isis.unige.ch/nasac/ -i
[...]
Mount(0): nasac on isis.unige.ch -> smb://isis.unige.ch/nasac/
  Type: GDaemonMount
  default_location=smb://isis.unige.ch/nasac/hpc_exchange/backup <=== here the relevant information about mount point
  themed icons:  [folder-remote]  [folder]  [folder-remote-symbolic]  [folder-symbolic]
  symbolic themed icons:  [folder-remote-symbolic]  [folder-symbolic]  [folder-remote]  [folder]
  can_unmount=1
  can_eject=0
  is_shadowed=0
```


## Post 15 by @Yann.Sagon (2024-12-20T10:08:05.542Z)

One more things to check: be sure to use the correct gio and not the one provided for example by anaconda.
```
(bamboo)-[root@login1 ~]$ which gio
/usr/bin/gio
```


## Post 16 by @Lucille.Delisle1 (2025-03-03T09:59:35.437Z)

Dear all,
Now that bamboo and baobab have been updated, I observe the same behaviour on both.
```
(baobab)-[delislel@login1 ~]$ gio mount smb://nasac-m2.unige.ch/m-AndreyLab    
Authentication Required
Enter user and password for share “m-andreylab” on “nasac-m2.unige.ch”:
User [delislel]: 
Domain [SAMBA]: ISIS
Password: 
(baobab)-[delislel@login1 ~]$ ls .gvfs
ls: cannot access '.gvfs': No such file or directory
(baobab)-[delislel@login1 ~]$ ls /run/user/
240477
```
I cannot find where it is mounted but I see what is inside:
```
(baobab)-[delislel@login1 ~]$ gio list smb://nasac-m2.unige.ch/m-andreylab
#SHARE
```
With another nasac I can even go in subdirectories:
```
(baobab)-[delislel@login1 ~]$ gio list smb://nasac-m2.unige.ch/m-gherrera/GHerrera
directory1
directory2
LucilleDelisle
```
And I know that I can copy with:
```
(baobab)-[delislel@login1 ~]$ gio copy smb://nasac-m2.unige.ch/m-gherrera/LucilleDelisle/xxx/xxx/results_20241211/reports/X9_report-cutadapt.txt ./
```
But for the ‘#SHARE’ gio list is doing very strange thing:
```
(baobab)-[delislel@login1 ~]$  gio list smb://nasac-m2.unige.ch/m-andreylab/\#SHARE/\#SHARE/\#SHARE/\#SHARE/
#SHARE
```
Do you have a solution either to have a real mount point or at least to be able to list files where the directory is `#SHARE`?
Thanks


## Post 17 by @Adrien.Albert (2025-03-03T12:38:40.124Z)

Hi @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)
Indeed due to the maintenance on Baobab the workarround was no longer effective. Can you try again?


## Post 18 by @Adrien.Albert (2025-03-03T12:57:58.484Z)

also, did you check:
`ls /run//user/${UID}/gvfs`


## Post 19 by @Lucille.Delisle1 (2025-03-03T13:36:00.507Z)

I restarted from scratch:
```
(baobab)-[delislel@login1 ~]$ dbus-launch bash
(baobab)-[delislel@login1 ~]$ gio mount smb://nasac-m2.unige.ch/m-gherrera
Authentication Required
Enter user and password for share “m-gherrera” on “nasac-m2.unige.ch”:
User [delislel]: 
Domain [SAMBA]: ISIS
Password: 
(baobab)-[delislel@login1 ~]$ ls -alh .gvfs
ls: cannot access '.gvfs': No such file or directory
(baobab)-[delislel@login1 ~]$ ls -alh /run/user/$UID/gvfs
ls: cannot access '/run/user/313457/gvfs': No such file or directory
(baobab)-[delislel@login1 ~]$ 
```
All results are the same.


## Post 20 by @Adrien.Albert (2025-03-04T11:10:57.885Z)

Hi @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)
Could you give the outpout of the following command:
```
$ echo $XDG_RUNTIME_DIR
```
If empty set the variable to:
```
 $ export XDG_RUNTIME_DIR=/run/user/$(id -u)
```


## Post 21 by @Adrien.Albert (2025-03-04T11:38:19.632Z)

It seems the issue is only occuring on login node
```
Login node:
-------------
(baobab)-[alberta@login1 ~]$ dbus-launch bash
dbus[1533946]: Unable to set up transient service directory: XDG_RUNTIME_DIR "/run/user/401775" not available: No such file or directory
(baobab)-[alberta@login1 ~]$ 
exit

Compute node:
-------------
(baobab)-[alberta@login1 ~]$ salloc
salloc: Granted job allocation 15269387
salloc: Nodes cpu001 are ready for job
(baobab)-[alberta@cpu001 ~]$ dbus-launch bash
(baobab)-[alberta@cpu001 ~]$ gio mount smb://isis.unige.ch/nasac/faculty/alberta-smb
Authentication Required
Enter user and password for share “nasac” on “isis.unige.ch”:
User [alberta]: 
Domain [SAMBA]: ISIS
Password: 
(baobab)-[alberta@cpu001 ~]$ ls /run/user/401775/gvfs
'smb-share:server=isis.unige.ch,share=nasac'
```
XDG_RUNTIME_DIR does not initialise correctly on the login node.
```
XDG_RUNTIME_DIR "/run/user/401775" not available: No such file or directory
```
