# [bamboo] gpu010 - cannot ssh to node

- Source: https://hpc-community.unige.ch/t/4300

- Created: 2026-06-04T20:26:38.724Z

- Tags: bamboo

- Posts: 3

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2026-06-04T20:26:38.798Z)

Dear HPC Team,
I have a job running on gpu010 and wanted to SSH to check resource utilization. Unfortunately, I got
```
$ ssh gpu010
Connection closed by 192.168.102.204 port 22
```
SSH to gpu030 works fine.
Thank you in advance,
Kind regards,
Maciej


## Post 2 by @Yann.Sagon (2026-07-03T09:57:14.194Z)

Dear  @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
sorry for the long delay. As you may have read, we’ve changed a lot of things about ssh connections those days. Is this still not working for you? I’ve tried to connect to gpu010 as user but the node is in use thus I can’t check.


## Post 3 by @maciej.falkiewicz (2026-07-13T09:23:57.682Z)

Works now. Thx @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) !
Best,
Maciej
