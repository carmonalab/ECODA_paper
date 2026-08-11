# PuttyFatalError

- Source: https://hpc-community.unige.ch/t/3838

- Created: 2025-02-24T15:56:04.223Z

- Posts: 6

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Helena.Bach (2025-02-24T15:56:04.272Z)

am unable to access PuTTY anymore. When I open PuTTY and attempt to connect, I receive the following error message:
“Server’s host key did not match the signature supplied.”
I have verified that the IP address (login1.baobab.hpc.unige.ch), port number, and username are correct. However, I still cannot establish a connection.
Could you please assist me in resolving this issue? Let me know if you need any additional information.


## Post 2 by @Adrien.Albert (2025-02-24T16:32:33.904Z)

Hi Elena
On our side I do not see any error.
```
Feb 24 16:50:30 login1 sshd[1472643]: pam_sss(sshd:auth): authentication success; logname= uid=0 euid=0 tty=ssh ruser= rhost=XXX.XXX.XXX.XXX user=bachh
Feb 24 16:50:30 login1 sshd[1472641]: Accepted keyboard-interactive/pam for bachh from XXX.XXX.XXX.XXX port 58528 ssh2
Feb 24 16:50:30 login1 sshd[1472641]: pam_unix(sshd:session): session opened for user bachh(uid=404644) by bachh(uid=0)
```
Could you try to connect with SSH, and OpenOnDemand[OpenOnDemand](https://openondemand.baobab.hpc.unige.ch/) ?


## Post 3 by @Helena.Bach (2025-02-24T16:58:27.281Z)

Adrien.Albert:
> SSH
Thank you for your reply.
I can access Filezilla but not Putty (so it seems to be something specific to Putty). Last week I was able to connect. I have tried to connect OpenOnDemand but it keeps me waiting.


## Post 4 by @Adrien.Albert (2025-02-24T17:07:57.847Z)

Hi @Helena.Bach[@Helena.Bach](https://hpc-community.unige.ch/u/helena.bach)
I have connected through putty with my password and it’s working. Could you try too ?
Could try to change your ssh key  and connect with the new one:
https://my-account.unige.ch/main/home[https://my-account.unige.ch/main/home](https://my-account.unige.ch/main/home)


## Post 5 by @Camille.Serquet (2025-02-25T09:32:42.984Z)

Dear @Helena.Bach[@Helena.Bach](https://hpc-community.unige.ch/u/helena.bach) , @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert)
same issue here…
I am not even accessing the step to enter my password.
Did you manage to find a solution in the meantime ?
PuTTyError
PuTTyError305×152 2.03 KB
[PuTTyError305×152 2.03 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/94c4e0efe433dba4c47157921de244b7b02c71d5.png)


## Post 6 by @Camille.Serquet (2025-02-25T11:19:21.793Z)

Ended up deleting PuTTy and re-installing it. It works now !
My version was 0.74.0.0 (64bit) and now, it is 0.83.0.0 (64-bit).
Hope it helps !
Cheers
