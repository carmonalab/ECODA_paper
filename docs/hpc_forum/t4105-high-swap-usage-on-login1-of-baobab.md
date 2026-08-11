# High SWAP Usage on Login1 of Baobab

- Source: https://hpc-community.unige.ch/t/4105

- Created: 2025-09-29T15:26:43.276Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Alexander.Froch (2025-09-29T15:26:43.336Z)

## Primary informations
Username: froch
Cluster: Baobab
## Description
The Login1 node is not usable for anything at the moment (not even proper login), because the SWAP is completely full again.
Due to that, I can’t even properly login on the login node and always get connection lost errors.
I know that changes were made to Baobab so that all users can only use 2 Cores and 2GBs of RAM. Is this maybe also possible for the SWAP?
Cheers,
Alex


## Post 2 by @pablo.strasser1 (2025-09-29T16:38:53.806Z)

Same cannot login or at least very difficult to login.
The load seem big to me.
Screenshot from 2025-09-29 18-41-26
Screenshot from 2025-09-29 18-41-262215×641 112 KB
[Screenshot from 2025-09-29 18-41-262215×641 112 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/d7787026df693344e55ef82798d1adaa2088b58c.png)
Just managed to login. Swap is full.


## Post 4 by @Yann.Sagon (2025-09-30T14:50:29.611Z)

Alexander.Froch:
> SWAP is completely full again.
Again? If you say so. I’m not aware of any SWAP issue we had in the past. We did have an issue on login1 but it was about the tmp space that was full. We fixed it today [2025] Current issues on HPC Cluster - #23 by Adrien.Albert[[2025] Current issues on HPC Cluster - #23 by Adrien.Albert](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/23) maybe this was the issue?


## Post 5 by @Alexander.Froch (2025-10-01T16:16:32.014Z)

Could be both. The only thing I was able to see was the completely full SWAP. This happened already before (I wrote emails the last time) and there it was the SWAP. Not sure if the two incidents are connected.
