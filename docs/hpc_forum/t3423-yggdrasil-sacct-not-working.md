# [yggdrasil] sacct not working

- Source: https://hpc-community.unige.ch/t/3423

- Created: 2024-04-19T15:47:37.920Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-04-19T15:47:37.967Z)

If you are asking for help, try to provide information that can help us solve your issue, such as :
what did you try: `$ sacct`
what was the error message:
```
sacct: error: _conn_readable: persistent connection for fd 3 experienced error[104]: Connection reset by peer
sacct: error: slurm_persist_conn_open: No response to persist_init
sacct: error: Sending PersistInit msg: Resource temporarily unavailable
sacct: error: Problem talking to the database: Resource temporarily unavailable
```
Kind regards,
Maciej Falkiewicz


## Post 2 by @maciej.falkiewicz (2024-04-19T15:48:41.232Z)

Ok, the problem was solved seconds after I opened the topic.


## Post 3 by @Adrien.Albert (2024-04-21T18:22:18.845Z)




## Post 4 by @Adrien.Albert (2024-04-21T18:23:39.952Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
Noted, I close the topic!
