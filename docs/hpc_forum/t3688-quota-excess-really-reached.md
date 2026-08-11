# Quota excess really reached?

- Source: https://hpc-community.unige.ch/t/3688

- Created: 2024-10-16T10:39:33.429Z

- Tags: bamboo

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maxime.Juventin (2024-10-16T10:39:33.467Z)

## Primary informations
Username:juventin
Cluster: Bamboo
## Description
Hello,
I am having a 'Disk quota exceeded" error while running git clone on my Bamboo Home. I checked my quota with your `beegfs-get-quota-home-scratch.sh` script, it indicates I have 2To on my Home. However, `du -sh  /home/users/j/juventin` found only 13Go. How could be the difference be this important ? Is it related to the current issues on Bamboo storage or am I missing something ?
Thank you
Greetings


## Post 2 by @Adrien.Albert (2024-10-17T10:37:07.766Z)

Hi @Maxime.Juventin[@Maxime.Juventin](https://hpc-community.unige.ch/u/maxime.juventin)
We are working ont it:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/26) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear HPC Users, 
Bamboo cluster is currently experiencing issues with quota on home filesystem. The symptom are that the disk usage may be incorrect. We are investigating. 
Thank you for your understanding. 
Best Regards, 
Status : Solved green_circle
start: 2024-10-16T22:02:00Z (UTC) 
stop: 2024-11-07T23:12:00Z (UTC)


## Post 3 by @Maxime.Juventin (2024-11-08T11:57:41.852Z)

Hello,
Today I noticed the problem seems to be solved
Thank you :slight_smile:


## Post 4 by @Adrien.Albert (2024-11-08T12:44:25.860Z)

Hi @Maxime.Juventin[@Maxime.Juventin](https://hpc-community.unige.ch/u/maxime.juventin)
Indeed we have found and correct the issue this morning :wink:
