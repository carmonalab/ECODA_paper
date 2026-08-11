# [bamboo] gpu006 nvidia driver in a bad state

- Source: https://hpc-community.unige.ch/t/4259

- Created: 2026-03-20T07:17:43.308Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ramon.CalvoGonzalez (2026-03-20T07:17:43.367Z)

## Primary informations
Username: calvogon
Cluster: bamboo
## Description
The NVIDIA driver is in a bad state.
It indicates xid error 95 (which is a memory corruption error).
## Steps to Reproduce
`dmesg` showed Xid 95 uncontained errors.
Running `cuInit(0)`returns error code 999 (CUDA_ERROR_UNKNOWN)
I think the fix would be to either reboot the machine or to try to run nvidia-smi –gpu-reset.
On the other hand, what happened to gpu005 on bamboo? I can’t see it anymore.
Best,
Ramon.


## Post 2 by @Adrien.Albert (2026-03-20T09:40:10.252Z)

Hello @Ramon.CalvoGonzalez[@Ramon.CalvoGonzalez](https://hpc-community.unige.ch/u/ramon.calvogonzalez)
gpu006 has been put back into production and seems working again.
Regarding gpu005, this node is reserved for a central Unige project, for technical reason the node has been remove from slurm.
