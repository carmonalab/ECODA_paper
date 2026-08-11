# GPU nodes drained on Yggdrasil and Bamboo — issue-7530: security

- Source: https://hpc-community.unige.ch/t/4296

- Created: 2026-05-22T05:55:24.838Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Stephane.Bernardini (2026-05-22T05:55:24.934Z)

User: bernasb3
Cluster: : bamboo & yggdrasil
---
## Traces — Yggdrasil (all GPU nodes)
Node
State
Partitions
Reason
Drain set by
LastBusyTime
gpu001
IDLE+DRAIN
public-interactive-gpu
issue-7530: security
root@2026-05-21T13:21:57
2026-05-21T14:11:38
gpu003
IDLE+DRAIN
public-gpu, shared-gpu
issue-7530: security
root@2026-05-21T13:21:57
2026-05-21T23:54:09
gpu004
IDLE+DRAIN
public-gpu, shared-gpu
issue-7530: security
root@2026-05-21T13:21:57
2026-05-21T22:58:37
gpu005
MIXED+DRAIN
public-gpu, shared-gpu
issue-7530: security
root@2026-05-21T13:21:57
2026-05-19T09:03:26
gpu006
IDLE+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:21:57
2026-05-21T22:44:07
gpu007
IDLE+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:21:57
2026-05-21T23:54:09
gpu008
IDLE+DRAIN
public-gpu, shared-gpu
issue-7530: security
root@2026-05-21T13:21:57
2026-05-14T19:42:44
## Traces — Bamboo (all GPU nodes)
Node
State
Partitions
Reason
Drain set by
LastBusyTime
gpu001
IDLE+DRAIN
public-interactive-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-21T07:53:51
gpu002
MIXED+DRAIN
public-gpu, shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-17T07:21:32
gpu003
IDLE+DRAIN
public-gpu, shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-22T01:16:07
gpu004
IDLE+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-21T14:52:28
gpu006
MIXED+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-11T11:14:27
gpu007
IDLE+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-21T14:52:28
gpu008
IDLE+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-21T14:52:28
gpu009
IDLE+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-21T14:52:28
gpu010
DOWN+DRAIN+NOT_RESPONDING
shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-19T22:39:47
gpu011
IDLE+DRAIN
shared-gpu
issue-7530: security
root@2026-05-21T13:17:52
2026-05-22T00:08:34
## Pending jobs affected
Yggdrasil: 6 jobs pending (`ReqNodeNotAvail`)
Bamboo: 14 jobs pending (`ReqNodeNotAvail` + `Priority`)
## Note
Baobab GPU nodes are not affected and remain operational. The issue appears to be specific to Yggdrasil and Bamboo.


## Post 2 by @Adrien.Albert (2026-05-23T23:01:42.453Z)

hi @Stephane.Bernardini[@Stephane.Bernardini](https://hpc-community.unige.ch/u/stephane.bernardini)
All gpu node has been drained du to a critical security vulnerability by the HPC team.
7530 is the issue ID on our tasks tracker :wink:
Most of the Baobab nodes have been already patched.
Here the reference:
[2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/21) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: All 
Dear users, 
We have been notified (again and again) of a critical security vulnerability in on Nvidia driver. As a result, we have drained all gpu node to apply an update. 
Status: Patching in progress orange_circle
start: 2026-05-21T10:00:00Z (UTC)
