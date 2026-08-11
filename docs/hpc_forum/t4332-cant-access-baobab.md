# Can't access baobab

- Source: https://hpc-community.unige.ch/t/4332

- Created: 2026-06-30T09:39:15.310Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Vincent.Riechers (2026-06-30T09:39:15.369Z)

Primary Information
Username: riechers
Cluster: baobab
Hi HPC team,
I can’t get a working session on Baobab today.
SSH authentication to both login1 and login2 succeeds, but the session freezes after auth and times out:
```
Authenticated to login1.baobab.hpc.unige.ch using "publickey".
debug1: Entering interactive session.
Timeout, server login1.baobab.hpc.unige.ch not responding.
```
Same on login2, and a non-interactive command (`ssh baobab echo OK`) hangs too. Maybe you can check this?
Thanks for your help!
Vincent


## Post 2 by @Vincent.Riechers (2026-06-30T11:39:23.290Z)

UPDATE: Problem seems to be fixed now, thanks


## Post 3 by @Yann.Sagon (2026-07-01T12:15:19.645Z)

Vincent.Riechers:
> Same on login2
This host doesn’t exists anymore since many years:)


## Post 4 by @Vincent.Riechers (2026-07-01T13:03:46.381Z)

Thanks, I didn’t know. I tried login2 and now see that it lands on login1.
