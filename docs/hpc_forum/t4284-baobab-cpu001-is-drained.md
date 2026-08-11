# [Baobab] cpu001 is drained

- Source: https://hpc-community.unige.ch/t/4284

- Created: 2026-04-22T14:31:33.891Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2026-04-22T14:31:33.969Z)

## Primary informations
Username: delislel
Cluster: baobab
## Description
cpu001 is currently drained and it is the only node of debug-cpu which is the default partition. Therefore jobs submitted without partition specification or with `-p debug-cpu` are not scheduled.
Thanks


## Post 2 by @Gael.Rossignol (2026-04-23T13:07:57.275Z)

Dear Lucille,
We had reserved the node to test the code. As you mentioned, it is currently the only one available in debug mode. I have now released the node.
Sorry for the inconvenience.
