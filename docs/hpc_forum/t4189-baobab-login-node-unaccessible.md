# Baobab login-node unaccessible

- Source: https://hpc-community.unige.ch/t/4189

- Created: 2026-01-12T14:03:01.982Z

- Tags: baobab

- Posts: 6

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2026-01-12T14:03:02.029Z)

## Primary informations
Username: coppinp
Cluster: Baobab
## Description
Hi all,
The Baobab login node seems to have gone down:
`pcoppin:~$ ping login1.baobab.hpc.unige.ch`
`PING login1.baobab.hpc.unige.ch (129.194.9.190): 56 data bytes`
`Request timeout for icmp_seq 0`
`Request timeout for icmp_seq 1`


## Post 2 by @Adrien.Albert (2026-01-12T14:04:18.944Z)

Hello @Paul.Coppin[@Paul.Coppin](https://hpc-community.unige.ch/u/paul.coppin) ;
Login1 has crashed and has been rebooted, it will be availbe within 5 min


## Post 3 by @Paul.Coppin (2026-01-12T14:04:54.016Z)

Thanks for the quick reply and response!


## Post 4 by @Davide.Piras (2026-01-12T14:49:39.759Z)

Hello, this seems to have happened again now, perhaps?
Thank you for your support!


## Post 5 by @Adrien.Albert (2026-01-12T15:00:44.450Z)

Yes, it’s happening again.
We suspect that a user is causing a kernel panic. We have identified what is triggering the issue, but we haven’t been able to identify the responsible user yet.


## Post 6 by @Davide.Piras (2026-01-12T15:02:09.779Z)

Ok, thank you for the reply!
