# [bamboo] gpu006 is down

- Source: https://hpc-community.unige.ch/t/4198

- Created: 2026-01-19T13:54:46.940Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ramon.CalvoGonzalez (2026-01-19T13:54:46.983Z)

Dear HPC team,
My jobs to gpu006 are not being assigned, and with `scontrol` I can see that the node seems to be down:
Reason=health_ps___blocked
Best, Ramón.


## Post 2 by @Ramon.CalvoGonzalez (2026-01-19T14:35:35.257Z)

Also, gpu009 is down:
Reason=health_cuda___GPU_broken


## Post 3 by @Gael.Rossignol (2026-01-20T08:53:54.866Z)

Hello,
Nodes have been set again in production, we have a lot of issues with gpu and I’m checking if there are nvidia driver issues.
Nest regards,
