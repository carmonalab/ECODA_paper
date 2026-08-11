# Bamboo maintenance?

- Source: https://hpc-community.unige.ch/t/4136

- Created: 2025-10-29T15:52:25.403Z

- Tags: bamboo

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Laure.Moinat (2025-10-29T15:52:25.461Z)

user: moinatl
Dear HPC team,
Maybe I missed something about the maintenance regarding Bamboo (I only received the two announcements for Baobab and Yggdrasil). Still, when I try to launch a simulation on private-gapnl-cpu, I get the message ‘(ReqNodeNotAvail, Reserved for maintenance)’.
Is this normal?
Best,
Laure


## Post 2 by @Yann.Sagon (2025-10-31T08:45:50.936Z)

Dear Laure,
I’ve checked your jobs on Bamboo: they are all running right now and there isn’t a specify reservation on your partition. Please show us the full message and context if you have the issue again.
Best


## Post 3 by @Laure.Moinat (2025-10-31T11:19:47.208Z)

Dear Yann,
Thanks for your answer. I have launched a simulation with 7 days on the private partition, however I get this : JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
```
       2794101 private-g H480.con  moinatl PD       0:00      1 (ReqNodeNotAvail, Reserved for maintenance)

       2793484 private-g H460.con  moinatl  R 1-18:17:26      1 cpu046

       2793414 private-g H500.con  moinatl  R 1-20:24:39      1 cpu046

       2793483 public-cp H360.con  moinatl  R 1-18:19:12      3 cpu\[011,013,015\] 
```
Is this normal ?
Best,
Laure


## Post 4 by @Yann.Sagon (2025-10-31T12:47:47.690Z)

Oups indeed, I forgot to announce this extra maintenance, this is now done : Bamboo extra maintenance: 04th November 2025[Bamboo extra maintenance: 04th November 2025](https://hpc-community.unige.ch/t/bamboo-extra-maintenance-04th-november-2025/4137)
What you can do is to submit your job on Baobab or request less running time.
Best regards
Yann
