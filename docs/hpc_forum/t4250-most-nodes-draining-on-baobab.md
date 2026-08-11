# Most nodes draining on Baobab

- Source: https://hpc-community.unige.ch/t/4250

- Created: 2026-03-16T07:44:44.928Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2026-03-16T07:44:44.983Z)

## Primary informations
Username: coppinp
Cluster: baobab
## Description
Hi, very few nodes seem to be up on Baobab. Most are in the draining state. E.g.
```
sinfo
shared-cpu                up   12:00:00    120  drain cpu[084-090,193-195,197-202,205-213,216-217,220-225,237-242,244,247-276,279-280,282-284,289-299,304,319-340,342-352]

scontrol show node cpu280
Reason=health_BEEGFS__tcp_con_storage [root@2026-03-14T21:45:25]
```
At some point during the weekend /srv/beegfs/dpnc got filled. Could this be the reason?
(space has already been liberated in the meantime)
All the best,
Paul


## Post 2 by @Adrien.Albert (2026-03-16T10:05:22.359Z)

Hello   @Paul.Coppin[@Paul.Coppin](https://hpc-community.unige.ch/u/paul.coppin) Nodes are back into production :slight_smile:
