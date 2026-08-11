# Issue with NASAC mount

- Source: https://hpc-community.unige.ch/t/4111

- Created: 2025-10-03T08:16:50.864Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicolas.Piron (2025-10-03T08:16:50.931Z)

Hi!
Username: pironn
Cluster: baobab
I used to transfer files between baobab and my team’s NAS and it’s not working anymore. The issue lies in the SMB share mount, and even though it was working two weeks ago (19.09.25), using the same commands today does not work anymore.
A possible source is the fact that I changed my UNIGE password, but I don’t get any error message when using my new password, while I get one when using the old one.
I normally use theses commands:
```
dbus-launch gio mount smb://nasac-m2.unige.ch/m-MyteamsNAS
```
```
# look in gvfs/ to see if it was correctly mounted
ls /run/user/$(id -u)/gvfs/
```
Usually this works and I’m able to access the NAS content in:
> /run/user/433544/gvfs/smb-share:server=nasac-m2.unige.ch,share=m-myteamsnas/
## Actual Result
I get the usual authentification request, and this seems to work as expected (no error message).
> Authentication Required
> Enter user and password for share “m-MyteamsNAS” on “nasac-m2.unige.ch”
> User [pironn]: pironn
> Domain [SAMBA]: isis.unige.ch
> Passeword: mynewpassword
Also I can connect from my personal laptop to the NAS without any problem using the same credentials.
But when I run the `ls` command to see if it is correctly mounted, the `gvfs/` folder is empty. And if I try to transfer some files, of course I get this error message:
```
rsync -rocvP $HOME/scratch/myproj /run/user/433544/gvfs/smb-share:server=nasac-m2.unige.ch,share=m-myteamsnas/myproject
```
> `sending incremental file list rsync: [Receiver] mkdir "/run/user/433544/gvfs/smb-share:server=nasac-m2.unige.ch,share=m-myteamsnas/myproject" failed: No such file or directory (2) rsync error: error in file IO (code 11) at main.c(783) [Receiver=3.2.3]`
NB: the issue is the same on Yggdrasil and Bamboo
Let me know if you can help me with this, thanks!


## Post 2 by @Yann.Sagon (2026-01-27T14:00:23.482Z)

Dear @Nicolas.Piron[@Nicolas.Piron](https://hpc-community.unige.ch/u/nicolas.piron) sorry, we missed your message. Apologies for Unanswered Posts — Notification Settings Updated[Apologies for Unanswered Posts — Notification Settings Updated](https://hpc-community.unige.ch/t/apologies-for-unanswered-posts-notification-settings-updated/4203)
We are aware that using gio isn’t reliable. Our advice is to be sure to have all the related process killed (dbus etc) and try again on a fresh session.
