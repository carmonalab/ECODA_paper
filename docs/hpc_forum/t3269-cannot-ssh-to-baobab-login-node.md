# Cannot ssh to Baobab login node

- Source: https://hpc-community.unige.ch/t/3269

- Created: 2024-01-24T16:04:55.866Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-01-24T16:04:55.928Z)

Hey,
I just noticed that the connection to Baobab is not working. My teammates confirmed it.
Kind regards,
Maciej Falkiewicz


## Post 2 by @Adrien.Albert (2024-01-24T16:19:50.092Z)

Hi,
I could not reproduce the issue. Could give the output of “ssh -vvv”


## Post 3 by @Yann.Sagon (2024-01-25T09:34:58.321Z)

Hi, it was probably related to high load on login2.baobab (see  Baobab Management node communication failed: admin1[Baobab Management node communication failed: admin1](https://hpc-community.unige.ch/t/baobab-management-node-communication-failed-admin1/3270) )
Fell free to comment.
