# [bamboo] ssh not working

- Source: https://hpc-community.unige.ch/t/3519

- Created: 2024-07-01T09:35:15.672Z

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-07-01T09:35:15.708Z)

Dear @support[@support](https://hpc-community.unige.ch/groups/support) ,
I tried to connect to the new cluster with
```
ssh <user>@login1.bamboo.hpc.unige.ch
```
but this doesn’t seem to work.
Kind regards,
Maciej Falkiewicz


## Post 2 by @Gael.Rossignol (2024-07-01T12:06:40.390Z)

Dear Maciej,
We already have a lot of people connected sucessfully, I check you profile and nothing is wrong. Do you have any error messages?
Best regards,


## Post 3 by @Adrien.Deline (2024-07-01T12:42:12.426Z)

Dear Maciej and Gael,
I had a similar issue with a timeout when trying ssh connection.
It turned out that it is not working from “outside” the University, i.e. it seems we have to use a VPN or another server entry to be able to access bamboo.
This is not the case for Yggdrasil or Baobab by the way.
Best regards,
Adrien


## Post 4 by @maciej.falkiewicz (2024-07-01T13:11:11.893Z)

In light of what @Adrien.Deline[@Adrien.Deline](https://hpc-community.unige.ch/u/adrien.deline) wrote: why is Bamboo reachable only from UNIGE’s network? Is this a feature or a bug, i.e. should we get used to the VPN?
Best regards,
Maciej Falkiewicz


## Post 5 by @maciej.falkiewicz (2024-07-02T12:09:06.151Z)

@support[@support](https://hpc-community.unige.ch/groups/support) is it planned to enable access to Bamboo from the outside of UNIGE network?
Kind regards,
Maciej Falkiewicz


## Post 6 by @Yann.Sagon (2024-07-02T13:55:57.989Z)

Hi, indeed the login node isn’t reachable from outside the university. This is not done on purpose, we are working to fix it. Please check the status here:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/17) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> BAMBOO: Login node not reachable from outside of the university
Dear users, 
it is not possible to connect to Bamboo login node from outside of the university. We’ll solve the issue. 
In the meantime your options are: 

use the UNIGE VPN[UNIGE VPN](https://catalogue-si.unige.ch/en/vpn) 

connect through Baobab or Yggdrasil login node 

Thank you for your understanding. 
Best regards, 
Status : Solved green_circle
start: 2024-07-01T09:30:00Z (UTC) 
end: 2024-07-04T07:00:00Z (UTC)


## Post 7 by @Yann.Sagon (2024-07-04T09:07:05.073Z)

Hi, login1.bamboo is now fully available from outside unige.
Best
