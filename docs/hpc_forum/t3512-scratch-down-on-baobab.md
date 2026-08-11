# Scratch down on baobab

- Source: https://hpc-community.unige.ch/t/3512

- Created: 2024-06-24T15:02:31.514Z

- Posts: 6

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Samuel.Klein (2024-06-24T15:02:31.551Z)

Hi,
I cannot access scratch on baobab, neither can my colleagues.
Cheers,
Sam


## Post 2 by @Adrien.Albert (2024-06-24T15:46:32.093Z)

Dear @Samuel.Klein[@Samuel.Klein](https://hpc-community.unige.ch/u/samuel.klein)
We are aware of this problem. We have found the root cause and are checking the file system before making it accessible again, to avoid any further problems.
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/16) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Baobab : Storage scratch down
Dear users, 
We are experiencing some problems on Baobab scratch storage; The root cause was the number of files open at the same time on the scratch. We had some issue to reboot server and get service up again but all is now resolved. 
Thank you for your understanding. 
Best regards, 

Status : Solved green_circle
Start: 2024-06-24T09:00:00Z (UTC) 
End: 2024-06-24T16:30:00Z (UTC)
Thank you for your understanding
Best regards,


## Post 3 by @Debajyoti.Sengupta (2024-06-25T05:59:00.269Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert),
It looks like that the files in scratch are still inaccessible, or at least nothing shows up when I list screen in the directory.
EDIT: As of now (8:52 AM), it is back again.


## Post 4 by @Gael.Rossignol (2024-06-25T06:54:02.889Z)

Debajyoti.Sengupta:
> e, or at least nothing shows up when I list screen in the directory.
Hello,
Issue was closed yesterday for the cluster, but login1 was still not connected to scratch. I restart services this morning.
Sorry for inconvenience.


## Post 5 by @Anton.Hanke (2024-06-25T08:12:35.083Z)

Hello,
I am sorry to raise this again, but scratch is still slow/inaccessible for me and members of the group rightnow from the login node.
Best wishes,
Anton


## Post 6 by @Gael.Rossignol (2024-06-25T09:09:25.657Z)

Anton.Hanke:
> to raise this again, but scratch is still slow/inaccessible for me and members of the group rightnow from the login node.
> Best wishes,
> Anton
Hello,
Actually scratch has too much data and is very slow. But all is working fine.
We will check with users if we can remove some data.
Best regards,
