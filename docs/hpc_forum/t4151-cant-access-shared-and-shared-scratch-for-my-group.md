# Can't access shared and shared scratch for my group

- Source: https://hpc-community.unige.ch/t/4151

- Created: 2025-11-24T08:38:28.074Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matija.Trickovic (2025-11-24T08:38:28.121Z)

## Primary informations
Username: trickovi
Cluster: Yggdrasil
## Description
Recently, I cant access my group’s shared (/home/share/Trajkovski_Group/) and shared scratch (/srv/beegfs/scratch/shares/Trajkovski_Group/) with the same error message:
-bash: cd: /home/share/Trajkovski_Group/: Permission denied
-bash: cd: /srv/beegfs/scratch/shares/Trajkovski_Group/: Permission denied
## Steps to Reproduce
cd  /home/share/Trajkovski_Group/
cd  /srv/beegfs/scratch/shares/Trajkovski_Group/
## Expected Result
To be able to cd into them.
## Actual Result
Permission denied.
Thanks!


## Post 2 by @Matija.Trickovic (2025-11-27T10:20:10.801Z)

Hello,
The issue persisted even after cluster maintenance.
Best,
Matija


## Post 3 by @Gael.Rossignol (2025-11-27T12:54:22.059Z)

Dear Matija,
When did you access to this folder last time? It seems rights are not set to use this folder.
Thanks for informations,
